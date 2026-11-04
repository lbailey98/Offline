// Calibration source: each event fires Poisson(poissonMean) photons (exactly one if unset) from (x,y,z),
// uniformly inside a cone around the aim direction. Each photon's energy is drawn from the isotope's gamma
// lines weighted by their intensities; coincident gammas from one decay are not modelled.
//
// Built-in line tables by pdgId (intensities in photons per 100 decays):
//   Co-60  (1000270600): 1173.2 (99.85), 1332.5 (99.98) keV
//   Na-22  (1000110220):  511.0 (181.7), 1274.5 (99.94) keV
//   Cs-137 (1000551370):  661.7 (85.10) keV
//   Ba-133 (1000561330):   81.0 (34.06), 276.4 (7.16), 302.9 (18.34), 356.0 (61.94), 383.8 (8.94) keV
//   Eu-152 (1000631520):  121.8, 244.7, 344.3, 411.1, 443.9, 778.9, 867.4, 964.1, 1085.8, 1112.1, 1408.0 keV
// Other isotopes: set gammaLines_MeV and gammaIntensities.
//
// Author: Leo Liu

#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <vector>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/SeedService/inc/SeedService.hh"

namespace mu2e {

  class IsotopeSourceGun : public art::EDProducer {
  public:
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;

    struct Config {
      fhicl::Atom<double> x {
        Name("x"), Comment("Source x position [mm]") };
      fhicl::Atom<double> y {
        Name("y"), Comment("Source y position [mm]") };
      fhicl::Atom<double> z {
        Name("z"), Comment("Source z position [mm]") };
      fhicl::Atom<int> pdgId {
        Name("pdgId"),
        Comment("PDG code of the isotope.  Selects the built-in gamma spectrum table. "
                "Example: Co-60 = 1000270600.") };
      fhicl::OptionalAtom<double> aimX {
        Name("aimX"), Comment("Aim direction x (unnormalised). Default: 0") };
      fhicl::OptionalAtom<double> aimY {
        Name("aimY"), Comment("Aim direction y (unnormalised). Default: 0") };
      fhicl::OptionalAtom<double> aimZ {
        Name("aimZ"), Comment("Aim direction z (unnormalised). Default: 1") };
      fhicl::OptionalAtom<double> coneHalfAngle {
        Name("coneHalfAngle"),
        Comment("Half-angle of emission cone around aim [deg]. 0 = pencil beam. Default: 0") };
      fhicl::OptionalSequence<double> gammaLines_MeV {
        Name("gammaLines_MeV"),
        Comment("Override built-in table: gamma energies [MeV].") };
      fhicl::OptionalSequence<double> gammaIntensities {
        Name("gammaIntensities"),
        Comment("Override built-in table: absolute intensities (photons per 100 decays), "
                "one per entry in gammaLines_MeV.  Used as sampling weights.") };
      fhicl::OptionalAtom<double> poissonMean {
        Name("poissonMean"),
        Comment("Mean photons per event, Poisson distributed. Default exactly one photon per event") };
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit IsotopeSourceGun(const Parameters& conf);
    void produce(art::Event& event) override;

  private:
    struct GammaLine { double energy_MeV; double intensity; };

    CLHEP::Hep3Vector pos_;
    CLHEP::Hep3Vector aim_;
    double            cosHalfAngle_;
    int               pdgId_;
    std::vector<GammaLine> lines_;
    std::vector<double>    cdf_;    // normalised CDF for energy sampling
    double            poissonMean_;
    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat   randFlat_;
    std::unique_ptr<CLHEP::RandPoissonQ> randPoisson_;

    static std::vector<GammaLine> defaultLines(int pdgId);
    void buildCDF();
    CLHEP::Hep3Vector randomDirection();
  };

  // Built-in gamma line tables

  std::vector<IsotopeSourceGun::GammaLine>
  IsotopeSourceGun::defaultLines(int pdgId) {
    // intensities = photons per 100 decays (absolute)
    static const std::map<int, std::vector<GammaLine>> kTable = {
      { 1000270600, { {1.17320, 99.85}, {1.33250, 99.98} } },              // Co-60
      { 1000110220, { {0.51100, 181.7}, {1.27450, 99.94} } },              // Na-22
      { 1000551370, { {0.66170, 85.10} } },                                // Cs-137
      { 1000561330, { {0.08097, 34.06}, {0.27640,  7.16},                  // Ba-133
                      {0.30293, 18.34}, {0.35600, 61.94},
                      {0.38385,  8.94} } },
      { 1000631520, { {0.12182, 28.41}, {0.24470,  7.55},                  // Eu-152
                      {0.34430, 26.59}, {0.41110,  2.24},
                      {0.44390,  3.13}, {0.77890, 12.94},
                      {0.86740,  4.24}, {0.96410, 14.63},
                      {1.08580, 10.21}, {1.11210, 13.64},
                      {1.40800, 20.87} } },
    };

    auto it = kTable.find(pdgId);
    if (it != kTable.end()) return it->second;

    std::string known;
    for (auto& kv : kTable) known += "  " + std::to_string(kv.first) + "\n";
    throw cet::exception("BADINPUT")
      << "IsotopeSourceGun: pdgId=" << pdgId << " has no built-in table.\n"
      << "Known isotopes:\n" << known
      << "Supply gammaLines_MeV and gammaIntensities to use an unlisted isotope.\n";
  }

  void IsotopeSourceGun::buildCDF() {
    double total = 0;
    for (auto& l : lines_) total += l.intensity;
    cdf_.resize(lines_.size());
    double running = 0;
    for (size_t i = 0; i < lines_.size(); ++i) {
      running += lines_[i].intensity / total;
      cdf_[i] = running;
    }
    cdf_.back() = 1.0;   // guard against floating-point rounding
  }


  IsotopeSourceGun::IsotopeSourceGun(const Parameters& conf)
    : art::EDProducer(conf)
    , pos_(conf().x(), conf().y(), conf().z())
    , pdgId_(conf().pdgId())
    , poissonMean_(-1.0)
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randFlat_(eng_)
  {
    produces<GenParticleCollection>();

    // Aim direction
    const double ax = conf().aimX() ? *conf().aimX() : 0.0;
    const double ay = conf().aimY() ? *conf().aimY() : 0.0;
    const double az = conf().aimZ() ? *conf().aimZ() : 1.0;
    const double amag = std::sqrt(ax*ax + ay*ay + az*az);
    if (amag < std::numeric_limits<double>::epsilon())
      throw cet::exception("BADINPUT") << "IsotopeSourceGun: aim direction has zero magnitude\n";
    aim_.set(ax/amag, ay/amag, az/amag);

    const double halfDeg = conf().coneHalfAngle() ? *conf().coneHalfAngle() : 0.0;
    if (halfDeg < 0.0 || halfDeg > 180.0)
      throw cet::exception("BADINPUT") << "IsotopeSourceGun: coneHalfAngle must be in [0,180] deg\n";
    cosHalfAngle_ = std::cos(halfDeg * CLHEP::degree);

    // Gamma lines: explicit override takes priority
    std::vector<double> energies, intensities;
    const bool hasE = conf().gammaLines_MeV(energies);
    const bool hasI = conf().gammaIntensities(intensities);
    if (hasE != hasI)
      throw cet::exception("BADINPUT")
        << "IsotopeSourceGun: gammaLines_MeV and gammaIntensities must both be provided together.\n";
    if (hasE) {
      if (energies.empty())
        throw cet::exception("BADINPUT") << "IsotopeSourceGun: gammaLines_MeV must not be empty.\n";
      if (energies.size() != intensities.size())
        throw cet::exception("BADINPUT")
          << "IsotopeSourceGun: gammaLines_MeV and gammaIntensities must have the same length.\n";
      for (size_t i = 0; i < energies.size(); ++i)
        lines_.push_back({energies[i], intensities[i]});
    } else {
      lines_ = defaultLines(pdgId_);
    }
    buildCDF();

    if (conf().poissonMean()) {
      poissonMean_ = *conf().poissonMean();
      if (poissonMean_ <= 0.0)
        throw cet::exception("BADINPUT") << "IsotopeSourceGun: poissonMean must be > 0.\n";
      randPoisson_ = std::make_unique<CLHEP::RandPoissonQ>(eng_, poissonMean_);
    }

    std::string lineStr;
    for (auto& l : lines_)
      lineStr += "  " + std::to_string(l.energy_MeV * 1e3) + " keV  (" +
                 std::to_string(l.intensity) + "%)\n";
    mf::LogInfo("IsotopeSourceGun")
      << "pdgId=" << pdgId_
      << "  pos=(" << pos_.x() << ", " << pos_.y() << ", " << pos_.z() << ") mm"
      << "  aim=(" << aim_.x() << ", " << aim_.y() << ", " << aim_.z() << ")"
      << "  coneHalfAngle=" << halfDeg << " deg"
      << "  poissonMean=" << (poissonMean_ > 0 ? poissonMean_ : 1.0)
      << (poissonMean_ > 0 ? "" : " (fixed)")
      << "\n  Gamma lines (sampled by intensity weight):\n" << lineStr;
  }


  CLHEP::Hep3Vector IsotopeSourceGun::randomDirection() {
    const double cosTheta = cosHalfAngle_ + (1.0 - cosHalfAngle_) * randFlat_.fire();
    const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta*cosTheta));
    const double phi      = CLHEP::twopi * randFlat_.fire();

    CLHEP::Hep3Vector u, v;
    if (std::abs(aim_.x()) < 0.9)
      u = aim_.cross(CLHEP::Hep3Vector(1, 0, 0)).unit();
    else
      u = aim_.cross(CLHEP::Hep3Vector(0, 1, 0)).unit();
    v = aim_.cross(u);

    return sinTheta*std::cos(phi)*u + sinTheta*std::sin(phi)*v + cosTheta*aim_;
  }


  void IsotopeSourceGun::produce(art::Event& event) {
    const long nPhotons = randPoisson_ ? randPoisson_->fire() : 1;

    auto output = std::make_unique<GenParticleCollection>();
    output->reserve(nPhotons);

    for (long i = 0; i < nPhotons; ++i) {
      // Sample one gamma energy from the discrete spectrum via the CDF.
      const double u = randFlat_.fire();
      size_t idx = 0;
      while (idx < cdf_.size() - 1 && u > cdf_[idx]) ++idx;
      const double E_MeV = lines_[idx].energy_MeV;

      const CLHEP::Hep3Vector dir = randomDirection();
      CLHEP::HepLorentzVector mom(E_MeV * dir, E_MeV);
      output->emplace_back(PDGCode::gamma, GenId::particleGun, pos_, mom, 0.0);
    }

    event.put(std::move(output));
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::IsotopeSourceGun)
