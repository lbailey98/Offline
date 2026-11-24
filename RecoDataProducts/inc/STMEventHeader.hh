#ifndef RecoDataProducts_STMEventHeader_hh_
#define RecoDataProducts_STMEventHeader_hh_

#include "artdaq-core-mu2e/Data/Mu2eEventHeader.hh"

namespace mu2e
{
  class STMEventHeader
  {
  public:
    explicit STMEventHeader() {}
    STMEventHeader(uint64_t ewt, uint8_t mode, uint64_t adcClock, uint64_t dtcClock, uint16_t rawLength = 0) :
      _mu2eEventHeader(ewt, mode), _adcClock(adcClock), _dtcClock(dtcClock), _rawLength(rawLength) { }

    uint64_t eventWindowTag() const { return _mu2eEventHeader.EventWindowTag; }
    uint8_t eventMode() const { return _mu2eEventHeader.EventMode; }

    uint64_t adcClock() const { return _adcClock; }
    uint64_t dtcClock() const { return _dtcClock; }

    // Number of ADC samples in the full raw window for this EWT (from the raw header,
    // filled even when the raw waveform itself is prescaled). 0 if unknown.
    uint16_t rawLength() const { return _rawLength; }

    // Ordered by EWT, then use remaining header fields to further classify
    bool operator<(STMEventHeader const& other) const {
      if (eventWindowTag() != other.eventWindowTag()) {
        return eventWindowTag() < other.eventWindowTag();
      }
      if (eventMode() != other.eventMode()) {
        return eventMode() < other.eventMode();
      }
      if (adcClock() != other.adcClock()) {
        return adcClock() < other.adcClock();
      }
      return dtcClock() < other.dtcClock();
    }

  private:
    Mu2eEventHeader _mu2eEventHeader;
    uint64_t _adcClock{0};
    uint64_t _dtcClock{0};
    uint16_t _rawLength{0};
  };
}  // namespace mu2e


#endif
