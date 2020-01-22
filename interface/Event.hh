
#ifndef EVENT_hh
#define EVENT_hh 1

#include <vector>
#include <Rtypes.h>

namespace grinder{
  
  // ========================================================= Muon ========================================================= 
  class Muon{
    public:
    Float_t pt, eta, phi, charge, isLoose, isMedium, isTight, relIsoTrk, relIsoPF;
  };

  // ========================================================= Electron ========================================================= 
  class Electron{
    public:
    Float_t pt, eta, phi, charge;
  };

  // ========================================================= Jet ========================================================= 
  class Jet{
    public:
    Float_t pt, eta, phi, charge;
  };

  // ========================================================= Photon ========================================================= 
  class Photon{
    public:
    Float_t pt, eta, phi;
  };

  // ========================================================= Event ========================================================= 
  class Event {
    public: 
    UInt_t run;
    UInt_t lumi;
    ULong64_t event;
    UShort_t bunchCrossing;
  };

};

#endif
