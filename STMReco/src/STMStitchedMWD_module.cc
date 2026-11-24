// Moving window deconvolution (MWD) run as a continuous stream over the STM digis
// produced by DAQ/src/STMDigisFromFragments_module.cc (raw or zero-suppressed).
//
// Differences from STMMovingWindowDeconvolution:
//  - Filter and peak-finder state persist across waveforms, EWTs and art events, so a pulse
//    straddling an EWT boundary is reconstructed once.
//  - Continuity is decided from the EWT and the sample index within the EWT. Contiguous ZS
//    regions within an EWT are already merged upstream, so a segment is only stitched to the
//    previous one if that reached the end of its window (STMEventHeader::rawLength) and this
//    one starts at index 0 of the next EWT. Otherwise the filter is reset and any peak in
//    progress is dropped:
//      raw: cold reset; no peaks are searched for in the first M samples (as online at the
//           start of a run), since the preceding signal is unknown (raw prescale gap).
//      zs : reset to the pedestal state (suppressed samples are below threshold), so peaks
//           are searched for immediately, as in the per-region algorithm.
//  - Filter and peak finding follow the online MWD (Mu2e-STMDAQ processing/mwd.cc): the
//    baseline (fcl baselineMean, else the STMEnergyCalib pedestal) is subtracted from the
//    ADCs, the threshold is -nsigma_cut * baselineSigma or fixedCut, peak finding starts
//    after M samples, and the peak time/EWT/height are those of the minimum.
//
// Output: STMPHDigiCollectionMap keyed by the STMEventHeader of the EWT containing the peak
// minimum; STMPHDigi::time() is the sample index within that EWT and energy() the rounded
// (negative) averaged MWD value, the same convention as the online PH digis. A peak that
// completes in a later art event is emitted in that later event, under the earlier header.
//
// Events must be processed in order on a single schedule.
//
// Algorithm: Claudia Alvarez-Garcia, Alex Keshavarzi, Mark Lancaster; Offline: Andy Edmonds,
// Pawel Plesniak. Streaming/stitching version adapted from STMMovingWindowDeconvolution.

// stdlib
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

// art
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline
#include "Offline/Mu2eUtilities/inc/STMUtils.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMEventHeader.hh"
#include "Offline/RecoDataProducts/inc/STMPHDigi.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/STMConditions/inc/STMEnergyCalib.hh"

// ROOT
#include "TTree.h"

namespace mu2e {

  namespace {

    // Streaming MWD filter + peak finder.
    class MWDStream {
    public:
      struct Params {
        unsigned M{0};
        unsigned L{0};
        double threshold{0};  // -nsigma * sigma, or a fixed cut
        double peakInit{0};   // initial peak minimum; a peak must go below this
      };

      struct Peak {
        STMEventHeader header;
        uint32_t index{0};  // sample index within the EWT of the peak minimum
        double height{0};   // averaged MWD value at the minimum (negative for a pulse)
      };

      explicit MWDStream(Params const& p)
        : _p(p), _dRing(p.M, 0.0), _diffRing(p.L, 0.0), _lowest(p.peakInit) {}

      void setConditions(double pedestal, double timeFactor) {
        _pedestal = pedestal;
        _timeFactor = timeFactor;
      }

      // Forget filter history; peaks are suppressed until M+L samples have been seen.
      // Returns true if a peak in progress was dropped.
      bool coldReset() {
        bool const dropped = _inPeak;
        _havePrev = false;
        _n = 0;
        _k = 0;
        _sumL = 0.0;
        _prevAvg = 0.0;
        clearPeak();
        return dropped;
      }

      // State equivalent to an infinitely long run of pedestal samples.
      void pedestalReset() {
        coldReset();
        _havePrev = true;
        _prevX = 0.0;
        _D = 0.0;
        std::fill(_dRing.begin(), _dRing.end(), 0.0);
        std::fill(_diffRing.begin(), _diffRing.end(), 0.0);
        _n = uint64_t(_p.M) + _p.L; // past the warm-up
      }

      bool inPeak() const { return _inPeak; }

      // Push one sample; appends to 'peaks' if a peak completes on this sample.
      void push(double adc, STMEventHeader const& header, uint32_t index,
                std::vector<Peak>& peaks) {
        // deconvolve (pedestal-subtracted, recursive)
        double const x = adc - _pedestal;
        _D = _havePrev ? (x - _timeFactor * _prevX + _D) : x;
        _prevX = x;
        _havePrev = true;

        // differentiate over M samples
        size_t const iM = _k % _p.M;
        double const diff = (_n >= _p.M) ? _D - _dRing[iM] : _D;
        _dRing[iM] = _D;

        // moving average over L samples
        size_t const iL = _k % _p.L;
        _sumL += diff - (_n >= _p.L ? _diffRing[iL] : 0.0);
        _diffRing[iL] = diff;
        if (iL == _p.L - 1) { // periodically recompute to stop rounding drift
          _sumL = 0.0;
          for (double d : _diffRing) _sumL += d;
        }
        double const avg = _sumL / _p.L;

        uint64_t const sampleIndex = _n; // samples since reset, before this one
        ++_k;
        ++_n;
        double const prevAvg = _prevAvg;
        _prevAvg = avg;
        if (sampleIndex < _p.M) return; // warm-up, as online (pf_index >= M)

        // peak finding, as online MWD::find_peaks: follow the minimum below threshold and
        // record its value, time and EWT; finalise when the signal rises back above threshold
        if (avg < _p.threshold && avg < prevAvg && avg < _lowest) {
          _lowest = avg;
          _inPeak = true;
          _peakHeader = header;
          _peakIndex = index;
        }
        if (_inPeak && avg > _p.threshold) {
          peaks.push_back(Peak{_peakHeader, _peakIndex, _lowest});
          clearPeak();
        }
      }

    private:
      void clearPeak() {
        _inPeak = false;
        _lowest = _p.peakInit;
      }

      Params _p;
      double _pedestal{0};
      double _timeFactor{0};

      // filter state
      bool _havePrev{false};
      double _prevX{0};
      double _D{0};
      std::vector<double> _dRing;     // last M deconvolved values
      std::vector<double> _diffRing;  // last L differentiated values
      double _sumL{0};
      double _prevAvg{0};
      uint64_t _n{0};                 // samples since last reset
      uint64_t _k{0};                 // ring index counter

      // peak state
      bool _inPeak{false};
      double _lowest{0};
      STMEventHeader _peakHeader;
      uint32_t _peakIndex{0};
    };

  } // namespace

  class STMStitchedMWD : public art::EDProducer {
  public:
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    struct Config {
      fhicl::Atom<art::InputTag> stmWaveformDigisMapTag{Name("stmWaveformDigisMapTag"), Comment("STMWaveformDigiCollectionMap from STMDigisFromFragments (rawHPGe, zsHPGe, rawLaBr or zsLaBr)")};
      fhicl::Atom<std::string> mode{Name("mode"), Comment("\"raw\" or \"zs\": how gaps between segments are treated")};
      fhicl::Atom<double> tau{Name("tau"), Comment("Decay constant of the waveform [ns]")};
      fhicl::Atom<double> M{Name("M"), Comment("Differentiation window [samples]")};
      fhicl::Atom<double> L{Name("L"), Comment("Averaging window [samples]")};
      fhicl::OptionalAtom<double> baselineMean{Name("baselineMean"), Comment("Baseline subtracted from the ADCs before deconvolution [ADC] (online manual_mean). If absent, the STMEnergyCalib pedestal is used")};
      fhicl::OptionalAtom<double> fixedCut{Name("fixedCut"), Comment("Fixed peak threshold on the averaged MWD signal [ADC] (online fixed_cut_value). If absent, -nsigma_cut * baselineSigma is used")};
      fhicl::OptionalAtom<double> nsigma_cut{Name("nsigma_cut"), Comment("Threshold in units of baselineSigma (online nsigma_cut)")};
      fhicl::OptionalAtom<double> baselineSigma{Name("baselineSigma"), Comment("Baseline sigma [ADC] (online manual_sigma)")};
      fhicl::Atom<double> peakMinInit{Name("peakMinInit"), Comment("Initial peak minimum; a peak must go below this (online peak_min_init)"), 0.0};
      fhicl::Atom<bool> makeTTreeEnergies{Name("makeTTreeEnergies"), Comment("Write a TTree of found peaks"), false};
      fhicl::Atom<double> TTreeEnergyCalib{Name("TTreeEnergyCalib"), Comment("Multiplier applied to energies in the TTree"), 1.0};
      fhicl::Atom<int> verbosityLevel{Name("verbosityLevel"), Comment("Verbosity level"), 0};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit STMStitchedMWD(Parameters const& conf);

  private:
    static double threshold(Config const& c);

    void produce(art::Event& event) override;
    void endJob() override;

    art::ProductToken<STMWaveformDigiCollectionMap> _digisToken;
    STMChannel _channel;
    bool _zsMode{false};
    double _tau{0};
    bool _haveManualBaseline{false};
    double _manualBaseline{0};
    unsigned _M{0};
    unsigned _L{0};
    int _verbosity{0};
    bool _makeTTree{false};
    double _ttreeCalib{1.0};

    ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h;
    MWDStream _mwd;

    // stream position: last sample fed
    bool _haveLast{false};
    STMEventHeader _lastHeader;
    uint64_t _lastEnd{0};     // index within the EWT one past the last sample fed

    // ttree
    TTree* _ttree{nullptr};
    uint64_t _tEWT{0};
    uint32_t _tTime{0};
    double _tE{0};
    uint32_t _tEventId{0};
    bool _tCrossEvent{false};

    // summary
    size_t _nEvents{0}, _nSegments{0}, _nSamples{0}, _nStitched{0}, _nGaps{0};
    size_t _nPeaks{0}, _nCrossEventPeaks{0}, _nDroppedPeaks{0};
  };

  STMStitchedMWD::STMStitchedMWD(Parameters const& conf)
    : art::EDProducer{conf},
      _digisToken(consumes<STMWaveformDigiCollectionMap>(conf().stmWaveformDigisMapTag())),
      _channel(STMUtils::getChannel(conf().stmWaveformDigisMapTag())),
      _tau(conf().tau()),
      _M(static_cast<unsigned>(conf().M())),
      _L(static_cast<unsigned>(conf().L())),
      _verbosity(conf().verbosityLevel()),
      _makeTTree(conf().makeTTreeEnergies()),
      _ttreeCalib(conf().TTreeEnergyCalib()),
      _mwd(MWDStream::Params{static_cast<unsigned>(conf().M()), static_cast<unsigned>(conf().L()),
                             threshold(conf()), conf().peakMinInit()}) {
    if (auto mean = conf().baselineMean()) {
      _haveManualBaseline = true;
      _manualBaseline = *mean;
    }
    std::string const mode = conf().mode();
    if (mode == "zs") _zsMode = true;
    else if (mode != "raw")
      throw cet::exception("Configuration") << "STMStitchedMWD: mode must be \"raw\" or \"zs\", got \"" << mode << "\"\n";
    if (_M == 0 || _L == 0 || _M < _L)
      throw cet::exception("Configuration") << "STMStitchedMWD: require 0 < L <= M (M=" << _M << ", L=" << _L << ")\n";
    if (_tau <= 0)
      throw cet::exception("Configuration") << "STMStitchedMWD: tau must be positive\n";

    produces<STMPHDigiCollectionMap>();

    if (_makeTTree) {
      art::ServiceHandle<art::TFileService> tfs;
      _ttree = tfs->make<TTree>("ttree", "STMStitchedMWD peaks");
      _ttree->Branch("ewt", &_tEWT, "ewt/l");
      _ttree->Branch("time", &_tTime, "time/i");
      _ttree->Branch("E", &_tE, "E/D");
      _ttree->Branch("eventId", &_tEventId, "eventId/i");
      _ttree->Branch("crossEvent", &_tCrossEvent, "crossEvent/O");
    }
  }

  double STMStitchedMWD::threshold(Config const& c) {
    if (auto cut = c.fixedCut()) return *cut;
    auto nsigma = c.nsigma_cut();
    auto sigma = c.baselineSigma();
    if (nsigma && sigma) return -(*nsigma) * (*sigma);
    throw cet::exception("Configuration") << "STMStitchedMWD: set either fixedCut, or both nsigma_cut and baselineSigma\n";
  }

  void STMStitchedMWD::produce(art::Event& event) {
    ++_nEvents;
    auto output = std::make_unique<STMPHDigiCollectionMap>();
    auto const& digisMap = *event.getValidHandle(_digisToken);

    STMEnergyCalib const& calib = _stmEnergyCalib_h.get(event.id());
    double const baseline = _haveManualBaseline ? _manualBaseline : calib.pedestal(_channel);
    _mwd.setConditions(baseline, 1.0 - calib.nsPerCt(_channel) / _tau); // tau_norm = 1 - tADC/tau

    std::vector<MWDStream::Peak> peaks;
    for (auto const& [header, digis] : digisMap) { // ordered by EWT
      bool firstInEWT = true;
      for (auto const& digi : digis) {
        auto const& adcs = digi.adcs();
        if (adcs.empty()) continue;
        ++_nSegments;
        uint64_t const start = _zsMode ? digi.trigTimeOffset() : 0;

        // Is this segment the direct continuation of the last sample fed? Contiguous segments
        // within an EWT are merged upstream, so this can only happen across an EWT boundary:
        // the previous segment reached the end of its window (if the length is known; older
        // files have rawLength 0) and this one starts the next window.
        bool const contiguous = _haveLast && firstInEWT
          && header.eventWindowTag() == _lastHeader.eventWindowTag() + 1
          && start == 0
          && (_lastHeader.rawLength() == 0 || _lastEnd == _lastHeader.rawLength());

        if (contiguous) {
          ++_nStitched;
        } else {
          if (_haveLast) ++_nGaps;
          if (_mwd.inPeak()) ++_nDroppedPeaks;
          // zs: suppressed samples are below threshold, so start from the pedestal state and
          //     search for peaks immediately (as the per-region algorithm did)
          // raw: real history is unknown, so wait for the filters to refill
          if (_zsMode) _mwd.pedestalReset();
          else _mwd.coldReset();
        }

        uint64_t index = start;
        for (int16_t adc : adcs) _mwd.push(adc, header, static_cast<uint32_t>(index++), peaks);
        _nSamples += adcs.size();
        _haveLast = true;
        _lastHeader = header;
        _lastEnd = index;
        firstInEWT = false;
      }
    }

    // Convert peaks to STMPHDigis
    for (auto const& peak : peaks) {
      // Same convention as the online PH: rounded, signed (negative) averaged MWD value.
      // Clip at the int16 range, since deconvolution of saturated pulses can overflow.
      double const h = std::clamp(peak.height, double(std::numeric_limits<int16_t>::min()),
                                  double(std::numeric_limits<int16_t>::max()));
      int16_t const energy = static_cast<int16_t>(std::lround(h));
      (*output)[peak.header].emplace_back(peak.index, energy);
      ++_nPeaks;
      bool const crossEvent = digisMap.find(peak.header) == digisMap.end();
      if (crossEvent) ++_nCrossEventPeaks;
      if (_makeTTree) {
        _tEWT = peak.header.eventWindowTag();
        _tTime = peak.index;
        _tE = energy * _ttreeCalib;
        _tEventId = event.id().event();
        _tCrossEvent = crossEvent;
        _ttree->Fill();
      }
      if (_verbosity > 3)
        mf::LogInfo("STMStitchedMWD") << _channel.name() << ": peak EWT " << peak.header.eventWindowTag()
                                      << " index " << peak.index << " energy " << energy
                                      << (crossEvent ? " (from previous art event)" : "");
    }
    if (_verbosity > 0)
      mf::LogInfo("STMStitchedMWD") << _channel.name() << ": event " << event.id() << ": "
                                    << peaks.size() << " peaks";

    event.put(std::move(output));
  }

  void STMStitchedMWD::endJob() {
    mf::LogInfo log("STMStitchedMWD");
    log << "============ STMStitchedMWD summary (" << _channel.name() << ", " << (_zsMode ? "zs" : "raw") << ") ============\n"
        << std::left
        << std::setw(32) << "\tEvents"                    << _nEvents << "\n"
        << std::setw(32) << "\tSegments"                  << _nSegments << "\n"
        << std::setw(32) << "\tSamples"                   << _nSamples << "\n"
        << std::setw(32) << "\tStitched across EWT"       << _nStitched << "\n"
        << std::setw(32) << "\tResets at gaps"            << _nGaps << "\n"
        << std::setw(32) << "\tPeaks"                     << _nPeaks << "\n"
        << std::setw(32) << "\tPeaks across art events"   << _nCrossEventPeaks << "\n"
        << std::setw(32) << "\tPeaks dropped at resets"   << _nDroppedPeaks + (_mwd.inPeak() ? 1 : 0) << "\n";
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::STMStitchedMWD)
