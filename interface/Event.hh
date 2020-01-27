
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
    Float_t ecalEnergyPreCorr, ecalEnergyPostCorr, energyScaleValue, energySigmaValue, energyScaleUp, energyScaleDown, energyScaleStatUp, energyScaleStatDown, energyScaleSystUp, energyScaleSystDown, energyScaleGainUp, energyScaleGainDown, energyScaleEtUp, energyScaleEtDown, energySigmaUp, energySigmaDown, energySigmaPhiUp, energySigmaPhiDown, energySigmaRhoUp, energySigmaRhoDown;
  };
  // ========================================================= Electron ========================================================= 
  class Electron{
    public:
    Float_t pt, eta, phi;
    Int_t charge;
    Bool_t  isLoose, isMedium, isTight;
    Float_t sumChargedHadronPt, sumNeutralHadronEt, sumPhotonEt, sumPUPt;
    Float_t ecalTrkEnergyPreCorr, ecalTrkEnergyPostCorr, energyScaleValue, energySigmaValue, energyScaleUp, energyScaleDown, energyScaleStatUp, energyScaleStatDown, energyScaleSystUp, energyScaleSystDown, energyScaleGainUp, energyScaleGainDown, energyScaleEtUp, energyScaleEtDown, energySigmaUp, energySigmaDown, energySigmaPhiUp, energySigmaPhiDown, energySigmaRhoUp, energySigmaRhoDown;
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
    Float_t pt, eta, phi, m, charge, area;
  };

  // ========================================================= Event ========================================================= 
  class Event {
    public: 
    UInt_t run;
    UInt_t lumi;
    ULong64_t event;
    UShort_t bunchCrossing;

    Float_t angular_pt_density;
  };

};

#endif
