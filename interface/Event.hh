
#ifndef EVENT_hh
#define EVENT_hh 1

#include <vector>
#include <Rtypes.h>

namespace grinder{
  // ========================================================= Photon ========================================================= 
  class Photon{
    public:
    Float_t pt, eta, phi;
    Bool_t  isLoose, isMedium, isTight;
    Float_t mva_value;
    Int_t   mva_category;
    Float_t sumChargedHadronPt, sumNeutralHadronEt, sumPhotonEt, sumPUPt;
  };
  // ========================================================= Electron ========================================================= 
  class Electron{
    public:
    Float_t pt, eta, phi, charge;
    Bool_t  isLoose, isMedium, isTight;
    Float_t sumChargedHadronPt, sumNeutralHadronEt, sumPhotonEt, sumPUPt;
  };
  // ========================================================= Muon ========================================================= 
  class Muon{
    public:
    Float_t pt, eta, phi, relIsoTrk, relIsoPF;
    Int_t charge;
    Bool_t isLoose, isMedium, isTight;
  };
  // ========================================================= Jet ========================================================= 
  class Jet{
    public:
    Float_t pt, eta, phi, charge;
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
