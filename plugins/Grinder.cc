// -*- C++ -*-

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/PatCandidates/interface/Muon.h>     // pat::Muon
#include <DataFormats/PatCandidates/interface/Tau.h>      // pat::Tau
#include <DataFormats/PatCandidates/interface/Photon.h>   // pat::Tau
#include <DataFormats/PatCandidates/interface/Electron.h> 
#include <DataFormats/PatCandidates/interface/Jet.h> 
#include <DataFormats/PatCandidates/interface/MET.h> 
#include <DataFormats/METReco/interface/GenMET.h>

// electrons / photons
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

// muons
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

// jets
#include <DataFormats/JetReco/interface/GenJet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <JetMETCorrections/Objects/interface/JetCorrectionsRecord.h>
#include <JetMETCorrections/Modules/interface/JetResolution.h>

// weights
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// vertexes
#include <DataFormats/VertexReco/interface/VertexFwd.h>   // reco::VertexCollection
#include <DataFormats/VertexReco/interface/Vertex.h>      // reco::Vertex

// pile-up
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

// trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// other
#include <FWCore/Utilities/interface/EDMException.h>      // edm::Exception

#include <FWCore/ServiceRegistry/interface/Service.h>     // TFileService
#include <CommonTools/UtilAlgos/interface/TFileService.h> // TFileService

#include <TTree.h>

// Grinder
#include "Analysis/GRINDER/interface/Event.hh"
#include "Analysis/GRINDER/interface/Grinder_extensions.hh"
using namespace grinder;

class Grinder : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Grinder(const edm::ParameterSet&);
      ~Grinder();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      void beginRun(edm::Run const &, edm::EventSetup const &setup);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::Service<TFileService> fileService;

      edm::EDGetTokenT<edm::View<pat::Photon>>    photonToken;
      edm::EDGetTokenT<edm::View<pat::Muon>>      muonToken;
      edm::EDGetTokenT<edm::View<pat::Electron>>  electronToken;
      edm::EDGetTokenT<edm::View<pat::Jet>>       jetToken;
      edm::EDGetTokenT<edm::View<pat::MET>>       metToken;
      edm::EDGetTokenT<edm::View<pat::Tau>>       tauToken;
      edm::EDGetTokenT<edm::View<reco::GenJet>>   genJetToken;
      edm::EDGetTokenT<edm::TriggerResults>         triggerResultsToken;
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_L1T_token;
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_HLT_token;
      edm::EDGetTokenT<edm::TriggerResults>         metFilterResultsToken;

      edm::EDGetTokenT<reco::VertexCollection> primaryVerticesToken;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puSummaryToken;

      TTree *outTree, *outTreeMeta;
      grinder::Event  event;
      grinder::Event *event_ptr;
      grinder::EventMetadata  eventMeta;
      grinder::EventMetadata *eventMeta_ptr;

      grinder::Jet jet;
      grinder::Muon muon;
      grinder::Photon photon;
      grinder::Electron electron;
      grinder::MET met;

      std::vector<grinder::Jet>      jets;
      std::vector<grinder::Muon>     muons;
      std::vector<grinder::Photon>   photons;
      std::vector<grinder::Electron> electrons;

      std::vector<grinder::Jet>      *jets_ptr;
      std::vector<grinder::Muon>     *muons_ptr;
      std::vector<grinder::Photon>   *photons_ptr;
      std::vector<grinder::Electron> *electrons_ptr;
      grinder::MET                   *met_ptr;

      std::string era_label;

      // Photons
      std::string photon_loose_id_token, photon_medium_id_token, photon_tight_id_token, photon_mva_token, photon_mva_token_val, photon_mva_token_cat;
      EffectiveAreas effAreaChHadrons;
      EffectiveAreas effAreaNeuHadrons;
      EffectiveAreas effAreaPhotons;
      float superCluster_eta;
      edm::Handle< double > rho, rhoCentral;
      edm::EDGetTokenT<double> rhoToken, rhoCentralToken;

      // Electrons
      std::string electron_loose_id_token, electron_medium_id_token, electron_tight_id_token;
      const reco::GsfElectron::PflowIsolationVariables * pfIso;
      // Jets
      // ID
      double NHF, NEMF, CHF, MUF, CEMF, NumConst, NumNeutralParticles, CHM;
      // JEC
      std::string jet_type_label;
      std::unique_ptr<JetCorrectionUncertainty> jecUnc;
      std::vector< std::string >  jecUnc_names;
      std::vector<JetCorrectionUncertainty*> jecUnc_v;
      // JER
      std::unique_ptr<JME::JetResolution> jerResolution_ptr;
      std::unique_ptr<JME::JetResolution> jerScaleFactor_ptr;

      JME::JetResolution jerResolution;
      JME::JetResolutionScaleFactor jerScaleFactor;
      JME::JetParameters jerResolution_parameters, jerScaleFactor_parameters;
      // MET
      edm::EDGetTokenT< bool > ecalBadCalibFilterUpdate_token ;
      // Cuts
      double cut_photon_pt, cut_photon_eta, cut_electron_pt, cut_electron_eta, cut_muon_pt, cut_muon_eta;
      double cut_jet_pt, cut_jet_eta;
      // Triggers
      std::vector<int> trigger_indexes;
      std::vector<std::string> selections_triggers_names;
      edm::ParameterSetID prevTriggerParameterSetID;
      bool do_trigger_filtering;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Grinder::Grinder(const edm::ParameterSet& iConfig) : 
  effAreaChHadrons( (iConfig.getParameter<edm::FileInPath>("effAreaChHadFile")).fullPath() ),
  effAreaNeuHadrons( (iConfig.getParameter<edm::FileInPath>("effAreaNeuHadFile")).fullPath() ),
  effAreaPhotons( (iConfig.getParameter<edm::FileInPath>("effAreaPhoFile")).fullPath() )
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  era_label     = iConfig.getParameter<std::string>("era_label");

  electronToken         = consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons_token"));
  tauToken              = consumes<edm::View<pat::Tau>>(iConfig.getParameter<edm::InputTag>("taus_token"));
  rhoToken              = consumes<double>(iConfig.getParameter<edm::InputTag>("rho_token"));
  rhoCentralToken       = consumes<double>(iConfig.getParameter<edm::InputTag>("rho_central_token"));
  genJetToken           = consumes<edm::View<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genjets_token"));
  primaryVerticesToken  = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertex_token"));
  puSummaryToken        = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("puSummaryToken_token"));
  triggerResultsToken   = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults_token"));
  triggerPrescales_L1T_token = consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescales_L1T_token"));
  triggerPrescales_HLT_token = consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescales_HLT_token"));
  selections_triggers_names  = iConfig.getParameter< std::vector<std::string> >("selections_triggers_names");
  for(unsigned int i = 0; i < selections_triggers_names.size(); i++){
    eventMeta.selections_triggers_names.push_back( selections_triggers_names.at(i) );
    event.trigger_fires.push_back( 0 );
  }

  // read trigger options
  do_trigger_filtering = iConfig.getParameter<bool>("do_trigger_filtering");

  // read Photons options
  photonToken   = consumes<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons_token"));
  photon_loose_id_token  = iConfig.getParameter<std::string>("photon_loose_id_token");
  photon_medium_id_token = iConfig.getParameter<std::string>("photon_medium_id_token");
  photon_tight_id_token  = iConfig.getParameter<std::string>("photon_tight_id_token");
  photon_mva_token       = iConfig.getParameter<std::string>("photon_mva_token");
  photon_mva_token_val   = photon_mva_token + "Values";
  photon_mva_token_cat   = photon_mva_token + "Categories";

  // read Electron options
  electron_loose_id_token  = iConfig.getParameter<std::string>("electron_loose_id_token");
  electron_medium_id_token = iConfig.getParameter<std::string>("electron_medium_id_token");
  electron_tight_id_token  = iConfig.getParameter<std::string>("electron_tight_id_token");

  // read Muon options
  muonToken     = consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons_token"));

  // read Jet options
  jetToken      = consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets_token"));
  // jecUnc = new JetCorrectionUncertainty( iConfig.getParameter<std::string>("jet_JEC_Uncertainty_datafile_token") );
  // https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Recommendation_for_analysis
  jet_type_label = iConfig.getParameter<std::string>("jet_type_label");
  if(era_label == "2016"){
    jecUnc_names = {"Total", "SubTotalMC", "SubTotalAbsolute", "SubTotalScale", "SubTotalPt", "SubTotalRelative", "SubTotalPileUp", "FlavorQCD", "TimePtEta"};
  }
  if(era_label == "2017"){
    jecUnc_names = {"Total", "SubTotalMC", "SubTotalAbsolute", "SubTotalScale", "SubTotalPt", "SubTotalRelative", "SubTotalPileUp", "FlavorQCD", "TimePtEta"};
  }
  if(era_label == "2018"){
    jecUnc_names = {"Total", "SubTotalMC", "SubTotalAbsolute", "SubTotalScale", "SubTotalPt", "SubTotalRelative", "SubTotalPileUp", "FlavorQCD", "TimePtEta"};
  }

  for(auto item : jecUnc_names){
    jet.JEC_unc_v_u.push_back( 0.f );
    jet.JEC_unc_v_d.push_back( 0.f );
  }

  // read met options
  metToken              = consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets_token"));
  metFilterResultsToken = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterResults_token"));
  ecalBadCalibFilterUpdate_token = consumes< bool >(edm::InputTag("ecalBadCalibReducedMINIAODFilter"));

  // read cuts
  cut_photon_pt     = iConfig.getParameter<double>("cut_photon_pt");
  cut_photon_eta    = iConfig.getParameter<double>("cut_photon_eta");
  cut_electron_pt   = iConfig.getParameter<double>("cut_electron_pt");
  cut_electron_eta  = iConfig.getParameter<double>("cut_electron_eta");
  cut_muon_pt       = iConfig.getParameter<double>("cut_muon_pt");
  cut_muon_eta      = iConfig.getParameter<double>("cut_muon_eta");
  cut_jet_pt        = iConfig.getParameter<double>("cut_jet_pt");
  cut_jet_eta       = iConfig.getParameter<double>("cut_jet_eta");

  // setup output data
  outTree     = fileService->make<TTree>("Events", "Events");
  outTreeMeta = fileService->make<TTree>("EventsMeta", "EventsMeta");

  event_ptr     = &event;
  jets_ptr      = &jets;
  muons_ptr     = &muons;
  photons_ptr   = &photons;
  electrons_ptr = &electrons;
  met_ptr       = &met;

  outTree->Branch("Event",     &event_ptr);
  outTree->Branch("Photons",   &photons_ptr);
  outTree->Branch("Electrons", &electrons_ptr);
  outTree->Branch("Muons",     &muons_ptr);
  outTree->Branch("Jets",      &jets_ptr);
  outTree->Branch("MET",       &met_ptr);

  eventMeta_ptr = &eventMeta;

  outTreeMeta->Branch("EventMeta", &eventMeta_ptr);
}


Grinder::~Grinder(){}

//
// member functions
//

// ------------ method called for each event  ------------

void Grinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  // Remove Info from previous events ========================================================================================================
  jets.clear();
  muons.clear();
  photons.clear();
  electrons.clear();
  event.weights.clear();
  event.ps_weights.clear();
  met.pt_unc_v_u.clear();
  met.pt_unc_v_d.clear();
  met.phi_unc_v_u.clear();
  met.phi_unc_v_d.clear();

  // Set Event Info ========================================================================================================
  event.run   = iEvent.id().run();
  event.lumi  = iEvent.luminosityBlock();
  event.event = iEvent.id().event();
  if(iEvent.isRealData()) event.bunchCrossing = iEvent.bunchCrossing();
  eventMeta.numEvents += 1;

  iEvent.getByToken(rhoToken, rho);
  event.angular_pt_density = (*rho);

  iEvent.getByToken(rhoCentralToken, rhoCentral);
  event.angular_pt_density_central = (*rhoCentral);

  // weights ========================================================================================================
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW#Retrieving_the_weights
  if(not iEvent.isRealData()){
    edm::Handle<GenEventInfoProduct> genEvtInfo; 
    iEvent.getByLabel( "generator", genEvtInfo );

    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideDataFormatGeneratorInterface?redirectedfrom=CMS.SWGuideDataFormatGeneratorInterface
    // use only one weight 
    // const std::vector<double> & evtWeights = genEvtInfo->weights();
    event.weight          = genEvtInfo->weight();
    eventMeta.sumWeights += genEvtInfo->weight();

    // get LHEInfo if it is from external LHE generator
    edm::Handle<LHEEventProduct> LHEInfo ;
    bool product_exists = iEvent.getByLabel( "externalLHEProducer", LHEInfo ) ;
    
    if( product_exists ){
      event.originalXWGTUP      = LHEInfo->originalXWGTUP();
      eventMeta.originalXWGTUP += LHEInfo->originalXWGTUP();
      const std::vector<gen::WeightsInfo> & lhe_weights = LHEInfo->weights();
      for(unsigned int i=0, N_weights = lhe_weights.size(); i < N_weights; i++)
        event.weights.push_back( lhe_weights[i].wgt );
    } else event.originalXWGTUP = 0;

    // alternative PS event weights
    const std::vector<double> & ps_weights = genEvtInfo->weights();
    if(not ps_weights.empty() and ps_weights.size() > 1){
      for(unsigned int i=0, N_weights = ps_weights.size(); i < N_weights; i++)
        event.ps_weights.push_back( ps_weights[i] );
    }
  }

  // Trigger filter ========================================================================================================
  edm::Handle<edm::TriggerResults> triggerResults;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales_L1T;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales_HLT;

  iEvent.getByToken(triggerResultsToken,   triggerResults);
  iEvent.getByToken(triggerPrescales_L1T_token, triggerPrescales_L1T);
  iEvent.getByToken(triggerPrescales_HLT_token, triggerPrescales_HLT);

  // update indexes if trigger menu is changed
  if( triggerResults->parameterSetID() != prevTriggerParameterSetID ){
    prevTriggerParameterSetID = triggerResults->parameterSetID();
    trigger_indexes.clear();

    const edm::TriggerNames & names = iEvent.triggerNames(*triggerResults);
    for(int i = 0, N_triggers = selections_triggers_names.size(); i < N_triggers; i++){
      trigger_indexes.push_back( -1 );
      std::string trigger_name = selections_triggers_names[i];
      for(int index = 0, N_all_triggers = names.size(); index < N_all_triggers; index++){
        std::string name = names.triggerName(index);
        if( trigger_name != name ) continue;
        trigger_indexes[i] = index;
        break;
      }
      if( trigger_indexes[i] < 0 ) std::cout << "GRINDER::WARNING, Triger " << trigger_name << " not in trigger menu" << std::endl;
    }
  }

  // check if event is selected
  bool passed = false;
  for(int i = 0, N_triggers = selections_triggers_names.size(); i < N_triggers; i++){
    int index = trigger_indexes[ i ];
    std::cout << i << " " << index << std::endl;
    event.trigger_fires[i] = 0;
    if( index < 0 ) continue;
    if( triggerResults->accept( index ) ){
      event.trigger_fires[i] = triggerPrescales_L1T->getPrescaleForIndex(index) * triggerPrescales_HLT->getPrescaleForIndex(index);
      if( event.trigger_fires[i] > 0 ) passed = true;
    }
  }

  if(not passed and do_trigger_filtering) return;

  /*
  for (unsigned int i = 0, n = triggerResults->size(); i < n; ++i) {
    std::cout << "Trigger " << names.triggerName(i) << ", prescale " << triggerPrescales->getPrescaleForIndex(i) << ": " << (triggerResults->accept(i) ? "PASS" : "fail (or not run)") << std::endl;
  }
  */

  // Read primary vertices collection ========================================================================================================
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(primaryVerticesToken, vertices);
  if (vertices->size() == 0){
    edm::Exception excp(edm::errors::LogicError);
    excp << "Event contains zero good primary vertices.\n";
    excp.raise();
  }
  const reco::Vertex PrimaryVertex = vertices->front();
  event.RecoNumInteractions = vertices->size();

  // Pile-Up Info ========================================================================================================
  // https://twiki.cern.ch/twiki/bin/view/CMS/Pileup_MC_Information
  // uncertanties calc at the next step https://twiki.cern.ch/twiki/bin/view/CMS/PileupSystematicErrors
  if(not iEvent.isRealData()){
    edm::Handle<std::vector<PileupSummaryInfo> > puSummary;
    iEvent.getByToken(puSummaryToken, puSummary);

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    event.TrueMCNumInteractions = -404;
    for(PVI = puSummary->begin(); PVI != puSummary->end(); ++PVI) {
      if( PVI->getBunchCrossing() != 0 ) continue;
      event.DicedMCNumInteractions = PVI->getPU_NumInteractions();
      event.TrueMCNumInteractions  = PVI->getTrueNumInteractions();
      break;
    }
  }

  // Iterate over photons ========================================================================================================
  edm::Handle<edm::View<pat::Photon>> srcPhotons;
  iEvent.getByToken(photonToken, srcPhotons);
  for (unsigned i = 0; i < srcPhotons->size(); ++i){
    pat::Photon const &ph = srcPhotons->at(i);

    if( ph.pt() < cut_photon_pt or TMath::Abs( ph.eta() ) > cut_photon_eta ) continue;

    photon.pt  = ph.pt();
    photon.eta = ph.eta();
    photon.phi = ph.phi();

    // cut based ID 
    // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#Accessing_ID_result
    photon.isLoose  = ph.photonID( photon_loose_id_token  );
    photon.isMedium = ph.photonID( photon_medium_id_token );
    photon.isTight  = ph.photonID( photon_tight_id_token  );

    // MVA base ID
    // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#Accessing_MVA_variables
    photon.mva_value    = ph.userFloat( photon_mva_token_val );
    photon.mva_category = ph.userInt( photon_mva_token_cat );

    // ISO https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolationRun2
    // https://github.com/varuns23/phoJetAnalysis/blob/master/phoJetNtuplizer/plugins/phoJetNtuplizer_photons.cc
    // https://arxiv.org/pdf/1502.02702.pdf
    superCluster_eta = fabs( ph.superCluster()->eta());
    photon.sumChargedHadronPt = std::max( 0.0, ph.userFloat("phoChargedIsolation")       - (*rho) * effAreaChHadrons.getEffectiveArea( superCluster_eta )  );
    photon.sumNeutralHadronEt = std::max( 0.0, ph.userFloat("phoNeutralHadronIsolation") - (*rho) * effAreaNeuHadrons.getEffectiveArea( superCluster_eta ) );
    photon.sumPhotonEt        = std::max( 0.0, ph.userFloat("phoPhotonIsolation")        - (*rho) * effAreaPhotons.getEffectiveArea( superCluster_eta )    );
    // FIXME sumPUPt ???

    // Energy Scale and Smearing
    // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#Energy_Scale_and_Smearing
    photon.ecalEnergyPreCorr   = ph.userFloat("ecalEnergyPreCorr");
    photon.ecalEnergyPostCorr  = ph.userFloat("ecalEnergyPostCorr");
    photon.energyScaleValue    = ph.userFloat("energyScaleValue");
    photon.energySigmaValue    = ph.userFloat("energySigmaValue");
    photon.energyScaleUp       = ph.userFloat("energyScaleUp");
    photon.energyScaleDown     = ph.userFloat("energyScaleDown");
    photon.energyScaleStatUp   = ph.userFloat("energyScaleStatUp");
    photon.energyScaleStatDown = ph.userFloat("energyScaleStatDown");
    photon.energyScaleSystUp   = ph.userFloat("energyScaleSystUp");
    photon.energyScaleSystDown = ph.userFloat("energyScaleSystDown");
    photon.energyScaleGainUp   = ph.userFloat("energyScaleGainUp");
    photon.energyScaleGainDown = ph.userFloat("energyScaleGainDown");
    photon.energyScaleEtUp     = ph.userFloat("energyScaleEtUp");
    photon.energyScaleEtDown   = ph.userFloat("energyScaleEtDown");
    photon.energySigmaUp       = ph.userFloat("energySigmaUp");
    photon.energySigmaDown     = ph.userFloat("energySigmaDown");
    photon.energySigmaPhiUp    = ph.userFloat("energySigmaPhiUp");
    photon.energySigmaPhiDown  = ph.userFloat("energySigmaPhiDown");
    photon.energySigmaRhoUp    = ph.userFloat("energySigmaRhoUp");
    photon.energySigmaRhoDown  = ph.userFloat("energySigmaRhoDown");

    photons.emplace_back( photon );
  }

  // iterate over electrons ========================================================================================================
  edm::Handle<edm::View<pat::Electron>> srcElectrons;
  iEvent.getByToken(electronToken, srcElectrons);
  for (unsigned i = 0; i < srcElectrons->size(); ++i){
    pat::Electron const &el = srcElectrons->at(i);

    if( el.pt() < cut_electron_pt or TMath::Abs( el.eta() ) > cut_electron_eta ) continue;

    electron.pt     = el.pt();
    electron.eta    = el.eta();
    electron.phi    = el.phi();
    electron.charge = el.charge();

    // cut based ID
    // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#Accessing_ID_result
    electron.isLoose  = el.electronID( electron_loose_id_token  );
    electron.isMedium = el.electronID( electron_medium_id_token );
    electron.isTight  = el.electronID( electron_tight_id_token  );

    // ISO https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolationRun2
    pfIso = & ( el.pfIsolationVariables() );
    electron.sumChargedHadronPt = pfIso->sumChargedHadronPt;
    electron.sumNeutralHadronEt = pfIso->sumNeutralHadronEt;
    electron.sumPhotonEt        = pfIso->sumPhotonEt;
    electron.sumPUPt            = pfIso->sumPUPt;

    // Energy Scale and Smearing
    // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#Energy_Scale_and_Smearing
    electron.ecalTrkEnergyPreCorr  = el.userFloat("ecalTrkEnergyPreCorr");
    electron.ecalTrkEnergyPostCorr = el.userFloat("ecalTrkEnergyPostCorr");
    electron.energyScaleValue      = el.userFloat("energyScaleValue");
    electron.energySigmaValue      = el.userFloat("energySigmaValue");
    electron.energyScaleUp         = el.userFloat("energyScaleUp");
    electron.energyScaleDown       = el.userFloat("energyScaleDown");
    electron.energyScaleStatUp     = el.userFloat("energyScaleStatUp");
    electron.energyScaleStatDown   = el.userFloat("energyScaleStatDown");
    electron.energyScaleSystUp     = el.userFloat("energyScaleSystUp");
    electron.energyScaleSystDown   = el.userFloat("energyScaleSystDown");
    electron.energyScaleGainUp     = el.userFloat("energyScaleGainUp");
    electron.energyScaleGainDown   = el.userFloat("energyScaleGainDown");
    electron.energyScaleEtUp       = el.userFloat("energyScaleEtUp");
    electron.energyScaleEtDown     = el.userFloat("energyScaleEtDown");
    electron.energySigmaUp         = el.userFloat("energySigmaUp");
    electron.energySigmaDown       = el.userFloat("energySigmaDown");
    electron.energySigmaPhiUp      = el.userFloat("energySigmaPhiUp");
    electron.energySigmaPhiDown    = el.userFloat("energySigmaPhiDown");
    electron.energySigmaRhoUp      = el.userFloat("energySigmaRhoUp");
    electron.energySigmaRhoDown    = el.userFloat("energySigmaRhoDown");

    electrons.emplace_back( electron );
  }

  // iterate over muons ========================================================================================================
  edm::Handle<edm::View<pat::Muon>> srcMuons;
  iEvent.getByToken(muonToken, srcMuons);

  for (unsigned i = 0; i < srcMuons->size(); ++i){
    pat::Muon const &mu = srcMuons->at(i);

    if( mu.pt() < cut_muon_pt or TMath::Abs( mu.eta() ) > cut_muon_eta ) continue;

    muon.pt     = mu.pt();
    muon.eta    = mu.eta();
    muon.phi    = mu.phi();
    muon.charge = mu.charge();

    // ID
    muon.isLoose  = mu.passed( reco::Muon::CutBasedIdLoose  );
    muon.isMedium = mu.passed( reco::Muon::CutBasedIdMedium );
    muon.isTight  = mu.passed( reco::Muon::CutBasedIdTight  );

    // ISO
    // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Identification
    muon.relIsoTrk = (mu.pfIsolationR04().sumChargedHadronPt + std::max(0., mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5*mu.pfIsolationR04().sumPUPt))/mu.pt();
    muon.relIsoPF  = mu.isolationR03().sumPt/mu.pt();

    // SF
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffsRun2

    muons.emplace_back( muon );
  }

  // Gen jets ========================================================================================================
  edm::Handle<edm::View<reco::GenJet>> genJets;
  if(not iEvent.isRealData())
    iEvent.getByToken(genJetToken, genJets);

  // iterate over AK4 PF CHS jets
  edm::Handle<edm::View<pat::Jet>> srcJets;
  iEvent.getByToken(jetToken, srcJets);
  for (unsigned i = 0; i < srcJets->size(); ++i){
    pat::Jet const &j = srcJets->at(i);

    // RAW P4 vector 
    // https://twiki.cern.ch/twiki/bin/view/CMS/TopJME#Jets
    reco::Candidate::LorentzVector const &rawP4 = j.correctedP4("Uncorrected");

    // FIXME how to correct jets ???

    // pT scale corrections
    // here you must use the CORRECTED jet pt
    // https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Recommendation_for_analysis
    for(int i = jecUnc_v.size()-1; i >= 0; i--){
      JetCorrectionUncertainty * jecUnc_ptr = jecUnc_v[i];
      jecUnc_ptr->setJetEta( j.eta() );
      jecUnc_ptr->setJetPt(  j.pt()  );
      jet.JEC_unc_v_u[i] = jecUnc->getUncertainty(true);  // FIXME not filled !!!
      jecUnc_ptr->setJetEta( j.eta() );
      jecUnc_ptr->setJetPt(  j.pt()  ); 
      jet.JEC_unc_v_d[i] = jecUnc->getUncertainty(false); // FIXME not filled !!!

      std::cout << jecUnc->getUncertainty(true) << std::endl;
    }
    double JEC_uncertanty = std::max( jet.JEC_unc_v_u.at( 0 ), jet.JEC_unc_v_d.at( 0 ) );
    JEC_uncertanty = std::max(JEC_uncertanty, 0.);

    // pT resolution
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
    double JER_uncertanty = 0;
    if (not iEvent.isRealData()) {
      jerResolution_parameters.setJetPt(j.pt()).setJetEta(j.eta());
      jerScaleFactor_parameters.setJetEta(j.eta()).setRho( *rho );

      jet.resolution = jerResolution.getResolution( jerResolution_parameters );
      jet.sf         = jerScaleFactor.getScaleFactor(jerScaleFactor_parameters);
      jet.sf_u       = jerScaleFactor.getScaleFactor(jerScaleFactor_parameters, Variation::UP);
      jet.sf_d       = jerScaleFactor.getScaleFactor(jerScaleFactor_parameters, Variation::DOWN);

      // https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=54#Smearing_procedures
      // https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_18/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L236-L237
      reco::GenJet const * genJet = MatchGenJet( j, genJets, 3 * jet.resolution * j.pt() );
      if (genJet) {
        jet.getJet_pt = genJet->pt();
        double energy_factor = ( rawP4.pt() - genJet->pt() ) / rawP4.pt();
        JER_uncertanty = std::max( { (jet.sf - 1.) * energy_factor, (jet.sf_u - 1.) * energy_factor, (jet.sf_d - 1.) * energy_factor } );
      } else {
        jet.getJet_pt = -1;
        double max_unc = std::max( { TMath::Abs(jet.sf), TMath::Abs(jet.sf_u), TMath::Abs(jet.sf_d) } );
        JER_uncertanty = jet.resolution * std::sqrt( std::max(std::pow(max_unc, 2) - 1., 0.) );
      }
    }

    double jet_maxPt = j.pt() * ( 1. + JEC_uncertanty + JER_uncertanty );
    std::cout << JEC_uncertanty << " " << JER_uncertanty << std::endl;
    if( rawP4.pt()  < cut_jet_pt and jet_maxPt < cut_jet_pt ) continue;
    if( TMath::Abs(rawP4.eta()) > cut_jet_eta ) continue;

    jet.pt  = rawP4.pt();
    jet.eta = rawP4.eta();
    jet.phi = rawP4.phi();
    jet.m   = rawP4.mass();

    jet.charge = j.jetCharge();
    jet.area   = j.jetArea();

    // JetID
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetID
    NHF  = j.neutralHadronEnergyFraction();
    NEMF = j.neutralEmEnergyFraction();
    CHF  = j.chargedHadronEnergyFraction();
    MUF  = j.muonEnergyFraction();
    CEMF = j.chargedEmEnergyFraction();
    NumConst = j.chargedMultiplicity()+j.neutralMultiplicity();
    NumNeutralParticles = j.neutralMultiplicity();
    CHM  = j.chargedMultiplicity();

    double eta = rawP4.eta();
    if(era_label == "2016"){
      if      ( TMath::Abs( eta ) < 2.7 ) jet.isTight = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=2.7;
      else if ( TMath::Abs( eta ) < 3.0 ) jet.isTight = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 && abs(eta)>2.7 && abs(eta)<=3.0 );
      else                                jet.isTight = (NEMF<0.90 && NumNeutralParticles>10 && abs(eta)>3.0 );
    }
    else if(era_label == "2017"){
      if      ( TMath::Abs( eta ) < 2.7 ) jet.isTight = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0) || abs(eta)>2.4) && abs(eta)<=2.7;
      else if ( TMath::Abs( eta ) < 3.0 ) jet.isTight = (NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 && abs(eta)>2.7);
      else                                jet.isTight = (NEMF<0.90 && NumNeutralParticles>10);
    }
    else if(era_label == "2018"){
      if      ( TMath::Abs( eta ) < 2.6 ) jet.isTight = (abs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 ); 
      else if ( TMath::Abs( eta ) < 2.7 ) jet.isTight = (abs(eta)>2.6 && abs(eta)<=2.7 && CEMF<0.8 && CHM>0 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 ); 
      else if ( TMath::Abs( eta ) < 3.0 ) jet.isTight = (NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 && abs(eta)>2.7 && abs(eta)<=3.0 );
      else                                jet.isTight = (NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 && abs(eta)>3.0 );
    }

    // b-tagging
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#B_tagging
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    // DeepCSV
    jet.pfDeepCSVJetTags_probb    = j.bDiscriminator("pfDeepCSVJetTags:probb");
    jet.pfDeepCSVJetTags_probbb   = j.bDiscriminator("pfDeepCSVJetTags:probbb");
    jet.pfDeepCSVJetTags_probc    = j.bDiscriminator("pfDeepCSVJetTags_probc");
    jet.pfDeepCSVJetTags_probudsg = j.bDiscriminator("pfDeepCSVJetTags:probudsg");

    // DeepJet FIXME not working for b tag
    jet.pfDeepFlavourJetTags_probb    = j.bDiscriminator("pfDeepFlavourJetTags_probb");
    jet.pfDeepFlavourJetTags_probbb   = j.bDiscriminator("pfDeepFlavourJetTags_probbb");
    jet.pfDeepFlavourJetTags_problepb = j.bDiscriminator("pfDeepFlavourJetTags_problepb");
    jet.pfDeepFlavourJetTags_probc    = j.bDiscriminator("pfDeepFlavourJetTags:probc");
    jet.pfDeepFlavourJetTags_probuds  = j.bDiscriminator("pfDeepFlavourJetTags:probuds");
    jet.pfDeepFlavourJetTags_probg    = j.bDiscriminator("pfDeepFlavourJetTags:probg");

    // b-tagging SF
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    // https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration#Standalone
    // will use standalone

    // b-tagging uncertanties
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods
    // will use standalone

    jets.emplace_back( jet );
  }

  // iterate over MET ========================================================================================================
  Handle<View<pat::MET>> metHandle;
  iEvent.getByToken(metToken, metHandle);
  pat::MET const & srcMET = metHandle->front();

  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#ETmiss
  met.pt  = srcMET.shiftedPt(pat::MET::NoShift, pat::MET::Type1);
  met.phi = srcMET.shiftedPhi(pat::MET::NoShift, pat::MET::Type1);
  met.significance = srcMET.metSignificance();
  
  if (not iEvent.isRealData()) {
      met.gen_pt  = srcMET.genMET()->pt() ;
      met.gen_phi = srcMET.genMET()->phi();

      using Var = pat::MET::METUncertainty;
      for (Var const &var: {Var::JetEnUp, Var::JetEnDown, Var::JetResUp, Var::JetResDown, Var::UnclusteredEnUp, Var::UnclusteredEnDown}) {
        met.pt_unc_v_u.push_back(  srcMET.shiftedPt(var, pat::MET::Type1) ); 
        met.pt_unc_v_d.push_back(  srcMET.shiftedPt(var, pat::MET::Type1) );
        met.phi_unc_v_u.push_back( srcMET.shiftedPhi(var, pat::MET::Type1) ); 
        met.phi_unc_v_d.push_back( srcMET.shiftedPhi(var, pat::MET::Type1) ); 
      }
  }

  // store the MET filter results
  // "if your analysis is sensitive to the MET tails, we strongly recommend you to apply all recommended MET filters"
  // https://twiki.cern.ch/twiki/bin/view/CMS/MissingET
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#ETmiss_filters
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Analysis_Recommendations_for_ana
  edm::Handle<edm::TriggerResults> metFilterResults;
  iEvent.getByToken(metFilterResultsToken, metFilterResults);
  const edm::TriggerNames &filter_names = iEvent.triggerNames(*metFilterResults);

  met.Flag_goodVertices                        =  false ;
  met.Flag_globalSuperTightHalo2016Filter      =  false ;
  met.Flag_HBHENoiseFilter                     =  false ;
  met.Flag_HBHENoiseIsoFilter                  =  false ;
  met.Flag_EcalDeadCellTriggerPrimitiveFilter  =  false ;
  met.Flag_BadPFMuonFilter                     =  false ;
  met.Flag_BadChargedCandidateFilter           =  false ;
  met.Flag_eeBadScFilter                       =  false ;
  met.Flag_ecalBadCalibReducedMINIAODFilter    =  false ;
  for (unsigned int i = 0, n = metFilterResults->size(); i < n; ++i) {
    if( not metFilterResults->accept(i) ) continue;
         if( filter_names.triggerName(i) == "Flag_goodVertices" )                        met.Flag_goodVertices                       = true ; 
    else if( filter_names.triggerName(i) == "Flag_globalSuperTightHalo2016Filter" )      met.Flag_globalSuperTightHalo2016Filter     = true ;
    else if( filter_names.triggerName(i) == "Flag_HBHENoiseFilter" )                     met.Flag_HBHENoiseFilter                    = true ;
    else if( filter_names.triggerName(i) == "Flag_HBHENoiseIsoFilter" )                  met.Flag_HBHENoiseIsoFilter                 = true ;
    else if( filter_names.triggerName(i) == "Flag_EcalDeadCellTriggerPrimitiveFilter" )  met.Flag_EcalDeadCellTriggerPrimitiveFilter = true ; 
    else if( filter_names.triggerName(i) == "Flag_BadPFMuonFilter" )                     met.Flag_BadPFMuonFilter                    = true ; 
    else if( filter_names.triggerName(i) == "Flag_BadChargedCandidateFilter" )           met.Flag_BadChargedCandidateFilter          = true ; 
    else if( filter_names.triggerName(i) == "Flag_eeBadScFilter" )                       met.Flag_eeBadScFilter                      = true ;
    else if( filter_names.triggerName(i) == "Flag_ecalBadCalibReducedMINIAODFilter" )    met.Flag_ecalBadCalibReducedMINIAODFilter   = true ; // FIXME need to re-run filter
  }

  // re-run MET filter for "Flag_ecalBadCalibReducedMINIAODFilter"
  edm::Handle< bool > passecalBadCalibFilterUpdate ;
  iEvent.getByToken(ecalBadCalibFilterUpdate_token, passecalBadCalibFilterUpdate);
  met.Flag_ecalBadCalibReducedMINIAODFilter = (*passecalBadCalibFilterUpdate );

  // Iterate over GenParticle  ========================================================================================================
  // TODO skipped for this analyses

  // Fill the output tree
  outTree->Fill();
}

// ------------ method called when starting to processes a run  ------------
void Grinder::beginRun(edm::Run const &, edm::EventSetup const &setup){
  // Construct an object to obtain JEC uncertainty
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections?rev=137#JetCorUncertainties
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  setup.get<JetCorrectionsRecord>().get(jet_type_label, JetCorParColl); 
  for(std::string name : jecUnc_names){
    JetCorrectorParameters const & JetCorPar = (*JetCorParColl)[ name.c_str() ];
    jecUnc_v.push_back( new JetCorrectionUncertainty(JetCorPar) );
  }
  // Objects that provide jet energy resolution and its scale factors
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution#Jet_resolution
  // jerResolution_ptr.reset(  new JME::JetResolution(           std::move(JME::JetResolution::get(setup,            "AK4PFchs_pt"))));
  // jerScaleFactor_ptr.reset( new JME::JetResolutionScaleFactor(std::move(JME::JetResolutionScaleFactor::get(setup, "AK4PFchs"   ))));
  jerResolution  = JME::JetResolution::get(setup,            "AK4PFchs_pt" );
  jerScaleFactor = JME::JetResolutionScaleFactor::get(setup, "AK4PFchs"    );
}

// ------------ method called once each job just before starting event loop  ------------
void Grinder::beginJob(){
}

// ------------ method called once each job just after ending the event loop  ------------
void Grinder::endJob(){
  outTreeMeta->Fill();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Grinder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Grinder);










