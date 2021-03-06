
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
    Float_t ptRaw, etaRaw, phiRaw, mRaw, charge, area;
    Float_t pt, eta, phi, m;
    std::vector<Float_t> JEC_unc_v_u, JEC_unc_v_d;
    Float_t resolution, sf, sf_u, sf_d, getJet_pt;
    Bool_t isTight;

    Float_t pfDeepCSVJetTags_probb, pfDeepCSVJetTags_probbb, pfDeepCSVJetTags_probc, pfDeepCSVJetTags_probudsg;
    Float_t pfDeepFlavourJetTags_probb, pfDeepFlavourJetTags_probbb, pfDeepFlavourJetTags_problepb, pfDeepFlavourJetTags_probc, pfDeepFlavourJetTags_probuds, pfDeepFlavourJetTags_probg;
  };
  // ========================================================= MET ========================================================= 
  class MET{
    public:
    Float_t pt, eta, phi, gen_pt, gen_phi, significance;
    std::vector<Float_t> pt_unc_v, phi_unc_v;
    Bool_t Flag_goodVertices, Flag_globalSuperTightHalo2016Filter, Flag_HBHENoiseFilter, Flag_HBHENoiseIsoFilter;
    Bool_t Flag_EcalDeadCellTriggerPrimitiveFilter, Flag_BadPFMuonFilter, Flag_BadChargedCandidateFilter, Flag_eeBadScFilter, Flag_ecalBadCalibReducedMINIAODFilter;
  };

  // ========================================================= Event ========================================================= 
  class Event {
    public: 
    UInt_t run;
    UInt_t lumi;
    ULong64_t event;
    UShort_t bunchCrossing;

    Float_t angular_pt_density, angular_pt_density_central;
    Float_t weight, originalXWGTUP;
    std::vector<Float_t> weights, ps_weights;
    Int_t DicedMCNumInteractions, TrueMCNumInteractions, RecoNumInteractions;
    
    std::vector<Float_t> trigger_fires;
  };
  // ========================================================= Event Meta ========================================================= 
  class EventMetadata {
    public: 
    EventMetadata(){
      numEvents = 0;
      sumWeights = 0;
      originalXWGTUP = 0;
    }

    bool is_data;
    ULong64_t numEvents;
    long double sumWeights, originalXWGTUP;
    std::vector<std::string> selections_triggers_names;
  };
};












#endif














