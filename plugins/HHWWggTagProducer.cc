// -*- C++ -*-

// #define DEBUG_GRINDER 0

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/HHWWggTag.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// #include "flashgg/Taggers/interface/GlobalVariablesDumper.h"
#include "flashgg/DataFormats/interface/PDFWeightObject.h"

#include "flashgg/MicroAOD/interface/CutBasedDiPhotonObjectSelector.h"

#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
#include "flashgg/DataFormats/interface/WHLeptonicTag.h"

#include "flashgg/Taggers/interface/LeptonSelection.h"

#include <vector>
#include <algorithm>
#include "TGraph.h"
#include "TLorentzVector.h"

//============================== GRINDER
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

// genParticles
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>

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

// Grinder
#include "Analysis/GRINDER/interface/Event.hh"
using namespace grinder;

#include <TTree.h>

using namespace std;
using namespace edm;

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
  // https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
  const reco::GenJet * MatchGenJet(pat::Jet const &jet, edm::Handle<edm::View<reco::GenJet>> const &genJets, double maxDPt) {
    reco::GenJet const *matchedJet = nullptr;
    double minDR2 = std::numeric_limits<double>::infinity();
    double const maxDR2 = 0.04; // jetConeSize * jetConeSize / 4.;
      
    for (unsigned i = 0; i < genJets->size(); ++i){
      const reco::GenJet & genJet = genJets->at(i);
      double const dR2 = ROOT::Math::VectorUtil::DeltaR2(jet.p4(), genJet.p4());
      if (dR2 > maxDR2 or dR2 > minDR2) continue;
      if (std::abs(jet.pt() - genJet.pt()) > maxDPt) continue;
      minDR2 = dR2;
      matchedJet = & ( genJets->at(i) );
    }
    return matchedJet;
  }

namespace flashgg {
  class GrinderHHWWggTagProducer : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources>  {
  public:
    GrinderHHWWggTagProducer( const ParameterSet & );

    virtual void beginJob() override;
    virtual void beginRun(edm::Run const&,  edm::EventSetup const&) override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endRun(edm::Run const&, edm::EventSetup const&)override;
    virtual void endJob() override;

    //void produce( edm::Event &, const EventSetup & ) override;
    //void beginJob( edm::Run const & run, edm::EventSetup const & setup ) override;

    bool checkPassMVAs(const flashgg::Photon*& leading_photon, const flashgg::Photon*& subleading_photon, edm::Ptr<reco::Vertex>& diphoton_vertex);
    std::vector<edm::EDGetTokenT<edm::View<DiPhotonCandidate> > > diPhotonTokens_;
    std::string inputDiPhotonName_;

    std::string inputJetsName_;
    std::vector<std::string> inputJetsSuffixes_;
    unsigned int inputJetsCollSize_;

    // Adding Jets 
    std::vector<edm::EDGetTokenT<edm::View<flashgg::Jet> > > jetTokens_;

    EDGetTokenT<View<Photon> > photonToken_;
    Handle<View<flashgg::Photon> > photons;

    EDGetTokenT<View<DiPhotonCandidate> > diphotonToken_;
    Handle<View<flashgg::DiPhotonCandidate> > diphotons;

    EDGetTokenT<View<reco::Vertex> > vertexToken_;
    Handle<View<reco::Vertex> > vertex;

    EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
    Handle<View<reco::GenParticle> > genParticle;

    EDGetTokenT<View<Electron> > electronToken_;
    Handle<View<flashgg::Electron> > electrons;

    EDGetTokenT<View<Muon> > muonToken_;
    Handle<View<flashgg::Muon> > muons;

    EDGetTokenT<View<flashgg::Met> > METToken_;
    Handle<View<flashgg::Met> > METs;

    EDGetTokenT<View<DiPhotonMVAResult> > mvaResultToken_;
    Handle<View<flashgg::DiPhotonMVAResult> > mvaResults;

    Handle<View<reco::Vertex> > vertices;

    EDGetTokenT<double> rhoTag_;
    edm::EDGetTokenT<edm::TriggerResults> triggerRECO_;
    edm::EDGetTokenT<edm::TriggerResults> triggerPAT_;
    edm::EDGetTokenT<edm::TriggerResults> triggerFLASHggMicroAOD_;
    string systLabel_;
    edm::Handle<double>  rho;

    // std::vector<edm::EDGetTokenT<View<flashgg::Jet> > > JetToken_;


    std::vector< std::string > systematicsLabels;
    std::vector<std::string> inputDiPhotonSuffixes_;

    //---ID selector
    ConsumesCollector cc_;
    edm::InputTag pdfWeight_;
    edm::EDGetTokenT<std::vector<flashgg::PDFWeightObject> > pdfWeightToken_;
    GlobalVariablesComputer globalVariablesComputer_;
    // CutBasedDiPhotonObjectSelector idSelector_;

    //----output collection
    // auto_ptr<vector<HHWWggCandidate> > HHWWggColl_;

    // variables from WHLeptonicTagProducer
    double leptonPtThreshold_;
    double muonEtaThreshold_;
    double leadPhoOverMassThreshold_;
    double subleadPhoOverMassThreshold_;
    double MVAThreshold_;
    double deltaRMuonPhoThreshold_;
    double jetsNumberThreshold_;
    double jetPtThreshold_;
    double jetEtaThreshold_;
    double muPFIsoSumRelThreshold_;
    double PhoMVAThreshold_;
    double METThreshold_;
    bool useVertex0only_;
    double deltaRJetMuonThreshold_;
    double deltaRPhoLeadJet_;
    double deltaRPhoSubLeadJet_;

    double DeltaRTrkElec_;
    double TransverseImpactParam_;
    double LongitudinalImpactParam_;

    double deltaRPhoElectronThreshold_;
    double deltaMassElectronZThreshold_;

    vector<double> nonTrigMVAThresholds_;
    vector<double> nonTrigMVAEtaCuts_;

    double electronIsoThreshold_;
    double electronNumOfHitsThreshold_;
    vector<double> electronEtaThresholds_;
    bool useElectronMVARecipe_;
    bool useElectronLooseID_;
    // string bTag_;
    bool doHHWWggTagCutFlowAnalysis_;


    edm::InputTag genInfo_;
    edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
    
    
    // GRINDER PART ========================================>
      // GlobalVariablesDumper * globalVarsDumper_;
      std::string era_label;
      double lumiWeight;

      edm::Service<TFileService> fileService;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puSummaryToken;
      edm::EDGetTokenT<edm::View<reco::GenJet>>   genJetToken;
      edm::EDGetTokenT<GenEventInfoProduct>         genToken;

      TTree *outTree, *outTreeMeta;
      grinder::Event  event;
      grinder::Event *event_ptr;
      grinder::EventMetadata  eventMeta;
      grinder::EventMetadata *eventMeta_ptr;

      grinder::Jet jet;
      grinder::Muon muon;
      grinder::Photon photon;
      grinder::Electron electron;
      grinder::MET out_met;
      grinder::GenParticle particle;

      std::vector<grinder::Jet>      out_jets;
      std::vector<grinder::Muon>     out_muons;
      std::vector<grinder::Photon>   out_photons;
      std::vector<grinder::Electron> out_electrons;
      std::vector<grinder::GenParticle> out_particles;

      std::vector<grinder::Jet>      *jets_ptr;
      std::vector<grinder::Muon>     *muons_ptr;
      std::vector<grinder::Photon>   *photons_ptr;
      std::vector<grinder::Electron> *electrons_ptr;
      std::vector<grinder::GenParticle> *particles_ptr;
      grinder::MET                   *met_ptr;
      
      // Electrons
      std::string electron_loose_id_token, electron_medium_id_token, electron_tight_id_token;
      const reco::GsfElectron::PflowIsolationVariables * pfIso;
      // Jets
      // ID
      double NHF, NEMF, CHF, MUF, CEMF, NumConst, NumNeutralParticles, CHM;
      // JEC
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
      bool do_trigger_filtering;    // ============ ========================================>
  };

    //---standard
    GrinderHHWWggTagProducer::GrinderHHWWggTagProducer( const ParameterSet & iConfig):
      photonToken_( consumes<View<Photon> >( iConfig.getParameter<InputTag> ( "PhotonTag" ) ) ),
      diphotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
      vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
      genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
      muonToken_( consumes<View<Muon> >( iConfig.getParameter<InputTag> ( "MuonTag" ) ) ),
      METToken_( consumes<View<Met> >( iConfig.getParameter<InputTag> ( "METTag" ) ) ),
      mvaResultToken_( consumes<View<flashgg::DiPhotonMVAResult> >( iConfig.getParameter<InputTag> ( "MVAResultTag" ) ) ),
      rhoTag_( consumes<double>( iConfig.getParameter<InputTag>( "rhoTag" ) ) ),
      triggerRECO_( consumes<edm::TriggerResults>(iConfig.getParameter<InputTag>("RECOfilters") ) ),
      triggerPAT_( consumes<edm::TriggerResults>(iConfig.getParameter<InputTag>("PATfilters") ) ),
      triggerFLASHggMicroAOD_( consumes<edm::TriggerResults>( iConfig.getParameter<InputTag>("FLASHfilters") ) ),
      systLabel_( iConfig.getParameter<string> ( "SystLabel" ) ),
      cc_( consumesCollector() ), // need absence of comma on last entry 
      pdfWeight_( iConfig.getUntrackedParameter<edm::InputTag>("flashggPDFWeightObject", edm::InputTag("flashggPDFWeightObject") ) ),
      pdfWeightToken_( cc_.consumes<std::vector<flashgg::PDFWeightObject> >( pdfWeight_ ) ),
      globalVariablesComputer_(iConfig.getParameter<edm::ParameterSet>("globalVariables"), cc_)
      // idSelector_( iConfig.getParameter<ParameterSet> ( "idSelection" ), cc_ )
    {

      // FIXME
      // save true HH mass

      puSummaryToken        = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("puSummaryToken_token"));

      inputDiPhotonName_= iConfig.getParameter<std::string > ( "DiPhotonName" );
      inputDiPhotonSuffixes_= iConfig.getParameter<std::vector<std::string> > ( "DiPhotonSuffixes" );
      std::vector<edm::InputTag>  diPhotonTags;
      diPhotonTags.push_back( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) );
      for (auto & suffix : inputDiPhotonSuffixes_){ 
          diPhotonTags.push_back(edm::InputTag(inputDiPhotonName_, suffix));
      }

      for( auto & tag : diPhotonTags ) { diPhotonTokens_.push_back( consumes<edm::View<flashgg::DiPhotonCandidate> >( tag ) ); }

      bool breg = 0;

      inputJetsName_= iConfig.getParameter<std::string> ( "JetsName" );
      inputJetsCollSize_= iConfig.getParameter<unsigned int> ( "JetsCollSize" );
      inputJetsSuffixes_= iConfig.getParameter<std::vector<std::string> > ( "JetsSuffixes" );
      // cout << "inputJetsCollSize_ = " << inputJetsCollSize_ << endl;
      if (breg){
        std::vector<edm::InputTag>  jetTags; // With bregression on 
        for (auto & suffix : inputJetsSuffixes_) {
            if (!suffix.empty()) systematicsLabels.push_back(suffix);  //nominal is already put in the diphoton loop
            for (unsigned int i = 0; i < inputJetsCollSize_ ; i++) {
                  std::string bregtag = suffix;
                  bregtag.append(std::to_string(i));
                  if (breg) jetTags.push_back(edm::InputTag(inputJetsName_,bregtag)); // With bregression on 
            }         
        }

        for( auto & tag : jetTags ) { jetTokens_.push_back( consumes<edm::View<flashgg::Jet> >( tag ) ); } // With bregression on 
      }

      // Jets without bregression 
      if (!breg){
        auto jetTags = iConfig.getParameter<std::vector<edm::InputTag> > ( "JetTags" ); 
        for( auto & tag : jetTags ) { jetTokens_.push_back( consumes<edm::View<flashgg::Jet> >( tag ) ); }
      }


      genInfo_ = iConfig.getUntrackedParameter<edm::InputTag>( "genInfo", edm::InputTag("generator") );
      genInfoToken_ = consumes<GenEventInfoProduct>( genInfo_ );

      leptonPtThreshold_ = iConfig.getParameter<double>( "leptonPtThreshold");
      muonEtaThreshold_ = iConfig.getParameter<double>( "muonEtaThreshold");
      MVAThreshold_ = iConfig.getParameter<double>( "MVAThreshold");
      deltaRMuonPhoThreshold_ = iConfig.getParameter<double>( "deltaRMuonPhoThreshold");
      jetsNumberThreshold_ = iConfig.getParameter<double>( "jetsNumberThreshold");
      jetPtThreshold_ = iConfig.getParameter<double>( "jetPtThreshold");
      jetEtaThreshold_ = iConfig.getParameter<double>( "jetEtaThreshold");
      muPFIsoSumRelThreshold_ = iConfig.getParameter<double>( "muPFIsoSumRelThreshold");
      PhoMVAThreshold_ = iConfig.getParameter<double>( "PhoMVAThreshold");
      METThreshold_ = iConfig.getParameter<double>( "METThreshold");
      useVertex0only_              = iConfig.getParameter<bool>("useVertex0only");
      deltaRJetMuonThreshold_ = iConfig.getParameter<double>( "deltaRJetMuonThreshold");
      deltaRPhoLeadJet_ = iConfig.getParameter<double>( "deltaRPhoLeadJet");
      deltaRPhoSubLeadJet_ = iConfig.getParameter<double>( "deltaRPhoSubLeadJet");

      DeltaRTrkElec_ = iConfig.getParameter<double>( "DeltaRTrkElec");
      TransverseImpactParam_ = iConfig.getParameter<double>( "TransverseImpactParam");
      LongitudinalImpactParam_ = iConfig.getParameter<double>( "LongitudinalImpactParam");

      deltaRPhoElectronThreshold_ = iConfig.getParameter<double>( "deltaRPhoElectronThreshold");
      deltaMassElectronZThreshold_ = iConfig.getParameter<double>( "deltaMassElectronZThreshold");

      nonTrigMVAThresholds_ =  iConfig.getParameter<vector<double > >( "nonTrigMVAThresholds");
      nonTrigMVAEtaCuts_ =  iConfig.getParameter<vector<double > >( "nonTrigMVAEtaCuts");
      electronIsoThreshold_ = iConfig.getParameter<double>( "electronIsoThreshold");
      electronNumOfHitsThreshold_ = iConfig.getParameter<double>( "electronNumOfHitsThreshold");
      electronEtaThresholds_ = iConfig.getParameter<vector<double > >( "electronEtaThresholds");
      useElectronMVARecipe_=iConfig.getParameter<bool>("useElectronMVARecipe");
      useElectronLooseID_=iConfig.getParameter<bool>("useElectronLooseID");
      // bTag_ = iConfig.getParameter<string> ( "bTag");
      doHHWWggTagCutFlowAnalysis_ = iConfig.getParameter<bool>( "doHHWWggTagCutFlowAnalysis");

      // produces<vector<HHWWggTag>>();
      // for (auto & systname : systematicsLabels) { // to deal with systematics in producer 
      //     produces<vector<HHWWggTag>>(systname);
      // }
      // produces<vector<TagTruthBase>>();
      
      // ============ ========================================> GRIDNER PART
        // globalVarsDumper_ = new GlobalVariablesDumper( iConfig.getParameter<edm::ParameterSet>( "globalVariables" ), std::forward<edm::ConsumesCollector>(cc) );

        // wei = xsec["xs"]/float(totEvents)*self.targetLumi
        // wei *= xsec.get("br",1.)
        // wei *= xsec.get("kf",1.)
        lumiWeight           = iConfig.getParameter<double>( "lumiWeight" );
        eventMeta.sumWeights = lumiWeight;
        // #ifdef DEBUG_GRINDER
          std::cout << "Grinder ... lumiWeight = " << eventMeta.sumWeights << std::endl;
        // #endif

        genJetToken       = consumes<edm::View<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genjets_token"));
        era_label         = iConfig.getParameter<std::string>("era_label");
        eventMeta.is_data = iConfig.getParameter<bool>("is_data");
        genToken          = consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>( "generator_token") ) ;

        outTree     = fileService->make<TTree>("Events", "Events");
        outTreeMeta = fileService->make<TTree>("EventsMeta", "EventsMeta");

        event_ptr     = &event;
        jets_ptr      = &out_jets;
        muons_ptr     = &out_muons;
        photons_ptr   = &out_photons;
        electrons_ptr = &out_electrons;
        met_ptr       = &out_met;
        particles_ptr = &out_particles;

        outTree->Branch("Event",     &event_ptr);
        outTree->Branch("Photons",   &photons_ptr);
        outTree->Branch("Electrons", &electrons_ptr);
        outTree->Branch("Muons",     &muons_ptr);
        outTree->Branch("Jets",      &jets_ptr);
        outTree->Branch("MET",       &met_ptr);
        outTree->Branch("Particles",       &particles_ptr);

        eventMeta_ptr = & eventMeta;

        outTreeMeta->Branch("EventMeta", &eventMeta_ptr);
        
        // read Electron options
        // electronToken = consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons_token"));
        electronToken_ = consumes<View<Electron> >( iConfig.getParameter<InputTag> ( "ElectronTag" ) );
        electron_loose_id_token  = iConfig.getParameter<std::string>("electron_loose_id_token");
        electron_medium_id_token = iConfig.getParameter<std::string>("electron_medium_id_token");
        electron_tight_id_token  = iConfig.getParameter<std::string>("electron_tight_id_token");

        // read Jet options
        // jecUnc = new JetCorrectionUncertainty( iConfig.getParameter<std::string>("jet_JEC_Uncertainty_datafile_token") );
        // https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Recommendation_for_analysis
        if(era_label == "2016"){
          // jecUnc_names = {"Total", "SubTotalMC", "SubTotalAbsolute", "SubTotalScale", "SubTotalPt", "SubTotalRelative", "SubTotalPileUp", "FlavorQCD", "TimePtEta"};
          jecUnc_names = {"Uncertainty"};
        }
        if(era_label == "2017"){
          // jecUnc_names = {"Total", "SubTotalMC", "SubTotalAbsolute", "SubTotalScale", "SubTotalPt", "SubTotalRelative", "SubTotalPileUp", "FlavorQCD", "TimePtEta"};
          jecUnc_names = {"Uncertainty"};
        }
        if(era_label == "2018"){
          // jecUnc_names = {"Total", "SubTotalMC", "SubTotalAbsolute", "SubTotalScale", "SubTotalPt", "SubTotalRelative", "SubTotalPileUp", "FlavorQCD", "TimePtEta"};
          jecUnc_names = {"Uncertainty"};
        }

        for(auto item : jecUnc_names){
          jet.JEC_unc_v_u.push_back( 0.f );
          jet.JEC_unc_v_d.push_back( 0.f );
        }

    }

    bool GrinderHHWWggTagProducer::checkPassMVAs( const flashgg::Photon*& leading_photon, const flashgg::Photon*& subleading_photon, edm::Ptr<reco::Vertex>& diphoton_vertex){
      // MVA Check variables 
      double lp_mva_thresh = 0.07;
      double slp_mva_thresh = -0.03;

      bool lead_pass_TightPhoID = 0, sublead_pass_TightPhoID = 0;
      double lp_Hgg_MVA = -99, slp_Hgg_MVA = -99; 
      double leading_pho_eta = -99, sub_leading_pho_eta = -99;

      // Get MVA values wrt diphoton vertex
      lp_Hgg_MVA  = leading_photon->phoIdMvaDWrtVtx( diphoton_vertex ); 
      slp_Hgg_MVA = subleading_photon->phoIdMvaDWrtVtx( diphoton_vertex ); 

      // Get eta values
      leading_pho_eta = leading_photon->p4().eta();
      sub_leading_pho_eta = subleading_photon->p4().eta();

      // leading photon 
      // EB 
      if (( abs(leading_pho_eta) > 0) && ( abs(leading_pho_eta) < 1.4442)){
        // if (lead_pho_EG_MVA_ > 0.42) lead_pass_TightPhoID = 1; 
        if (lp_Hgg_MVA > lp_mva_thresh) lead_pass_TightPhoID = 1; 
      }

      // EE 
      else if (( abs(leading_pho_eta) > 1.566) && ( abs(leading_pho_eta) < 2.5)){
        // if (lead_pho_EG_MVA_ > 0.14) lead_pass_TightPhoID = 1;
        if (lp_Hgg_MVA > slp_mva_thresh) lead_pass_TightPhoID = 1;
      }

      // SubLeading Photon
      // EB 
      if (( abs(sub_leading_pho_eta) > 0) && ( abs(sub_leading_pho_eta) < 1.4442)){
        // if (sublead_pho_EG_MVA_ > 0.42) sublead_pass_TightPhoID = 1; 
        if (slp_Hgg_MVA > lp_mva_thresh) sublead_pass_TightPhoID = 1; 
      }

      // EE 
      else if (( abs(sub_leading_pho_eta) > 1.566) && ( abs(sub_leading_pho_eta) < 2.5)){
        // if (sublead_pho_EG_MVA_ > 0.14) sublead_pass_TightPhoID = 1;
        if (slp_Hgg_MVA > slp_mva_thresh) sublead_pass_TightPhoID = 1;
      }

      if (lead_pass_TightPhoID && sublead_pass_TightPhoID){
        return 1;
      }

      else{
        return 0; 
      }

    }

    void GrinderHHWWggTagProducer::analyze( const edm::Event &iEvent, const EventSetup & ) {
      /*
          TODO:
           V Remove Info from previous events
           V Set Event Info
           V weights
           V Read primary vertices collection
           V Pile-Up Info
           V electrons
           V muons
           V photons
           V genjets
           V jets
           V met
           V save objects weights
      */

      #ifdef DEBUG_GRINDER
        std::cout << "GrinderHHWWggTagProducer::analyze() ... " << std::endl;
      #endif
      // Remove Info from previous events ========================================================================================================
      out_jets.clear();
      out_muons.clear();
      out_photons.clear();
      out_electrons.clear();
      out_met.pt_unc_v.clear();
      out_met.phi_unc_v.clear();
      out_particles.clear();

      event.trigger_fires.clear();
      event.weights.clear();
      event.ps_weights.clear();
      event.flashgg_mc_weights.clear();
      event.flashgg_diphoton_weights.clear();

      // Get particle objects ======================================================================================================== 
      iEvent.getByToken( photonToken_, photons );
      iEvent.getByToken( electronToken_, electrons );
      iEvent.getByToken( muonToken_, muons );
      iEvent.getByToken( METToken_, METs );
      iEvent.getByToken( mvaResultToken_, mvaResults );
      iEvent.getByToken( vertexToken_, vertices );
      iEvent.getByToken( rhoTag_, rho);

      std::vector<edm::Ptr<flashgg::Muon> >     allGoodMuons     = selectAllMuons( muons->ptrs(), vertices->ptrs(), muonEtaThreshold_, leptonPtThreshold_, muPFIsoSumRelThreshold_ );
      std::vector<edm::Ptr<flashgg::Electron> > allGoodElectrons = selectStdAllElectrons( electrons->ptrs(), vertices->ptrs(), leptonPtThreshold_, electronEtaThresholds_, useElectronMVARecipe_, useElectronLooseID_, *rho, iEvent.isRealData() );


      #ifdef DEBUG_GRINDER
        std::cout << "  event content :" << std::endl;
        std::cout << "    N photons = "   << photons->size() << std::endl;
        std::cout << "    N electrons = " << electrons->size() << std::endl;
        std::cout << "    N muons = " << muons->size() << std::endl;
        std::cout << "    N good electrons = " << allGoodMuons.size() << std::endl;
        std::cout << "    N good muons = " << allGoodElectrons.size() << std::endl;
        std::cout << "    rho = " << *rho << std::endl;
      #endif
      if( (allGoodMuons.size() + allGoodElectrons.size()) < 1 ) return;

      // Set Event Info ============================================================================================================== TODO
      event.run   = iEvent.id().run();
      event.lumi  = iEvent.luminosityBlock();
      event.event = iEvent.id().event();
      event.RecoNumInteractions = vertices->size();
      event.angular_pt_density = (*rho);
      if(iEvent.isRealData()) event.bunchCrossing = iEvent.bunchCrossing();
      eventMeta.numEvents += 1;

      // iEvent.getByToken(rhoCentralToken, rhoCentral);
      // event.angular_pt_density_central = (*rhoCentral);

      #ifdef DEBUG_GRINDER
        std::cout << "Run/Lumi/Event " << event.run << "/" << event.lumi << "/" << event.event << std::endl;
      #endif

      // Read the generator-level particles ========================================================================================================
      if(not iEvent.isRealData()){
        Handle<View<reco::GenParticle>> genParticles;
        iEvent.getByToken(genParticleToken_, genParticles);

        #ifdef DEBUG_GRINDER
          cout << "Process GenParticles, size = " << genParticles->size() << endl;
        #endif

        for (reco::GenParticle const &p: *genParticles){
          #ifdef DEBUG_GRINDER
            cout << "   " << p.pdgId() << " " << p.mother(0) << " " << p.status() << endl;
          #endif
          
          particle.pt  = p.pt(); 
          particle.eta = p.eta();
          particle.phi = p.phi();
          particle.m   = p.mass();
          particle.pdg_id = p.pdgId();
          particle.status = p.status();
          out_particles.emplace_back( particle );
        }
      }

      // weights ========================================================================================================
      event.flashgg_weight = lumiWeight;

      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW#Retrieving_the_weights
      if(not iEvent.isRealData()){
        #ifdef DEBUG_GRINDER
          std::cout << "weights ... " << std::endl;
        #endif
        edm::Handle<GenEventInfoProduct> genEvtInfo; 
        iEvent.getByToken(genToken, genEvtInfo);

        // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideDataFormatGeneratorInterface?redirectedfrom=CMS.SWGuideDataFormatGeneratorInterface
        // use only one weight 
        // const std::vector<double> & evtWeights = genEvtInfo->weights();
        event.weight          = genEvtInfo->weight();
        // eventMeta.sumWeights += genEvtInfo->weight();

        // get LHEInfo if it is from external LHE generator
        edm::Handle<LHEEventProduct> LHEInfo ;
        bool product_exists = iEvent.getByLabel( "externalLHEProducer", LHEInfo ) ;
        
        if( product_exists ){
          event.originalXWGTUP      = LHEInfo->originalXWGTUP();
          eventMeta.originalXWGTUP += LHEInfo->originalXWGTUP();
          const std::vector<gen::WeightsInfo> & lhe_weights = LHEInfo->weights();
          for(unsigned int i=0, N_weights = lhe_weights.size(); i < N_weights; i++)
            event.weights.push_back( lhe_weights[i].wgt );
        } else event.originalXWGTUP = -11;

        // alternative PS event weights
        const std::vector<double> & ps_weights = genEvtInfo->weights();
        if(not ps_weights.empty() and ps_weights.size() > 1){
          for(unsigned int i=0, N_weights = ps_weights.size(); i < N_weights; i++)
            event.ps_weights.push_back( ps_weights[i] );
        }

        // flashgg events ...
        edm::Handle<vector<flashgg::PDFWeightObject> > WeightHandle;
        // product_exists = iEvent.getByLabel(pdfWeight_, WeightHandle);
        product_exists = iEvent.getByToken(pdfWeightToken_, WeightHandle);
        #ifdef DEBUG_GRINDER
          std::cout << "flashgg weights are exists? ... " << product_exists << std::endl;
        #endif
        if(product_exists){
          #ifdef DEBUG_GRINDER
            std::cout << "N weights handlers ... " << (*WeightHandle).size() << std::endl;
          #endif
          for( unsigned int weight_index = 0; weight_index < (*WeightHandle).size(); weight_index++ ){
            vector<uint16_t> compressed_weights = (*WeightHandle)[weight_index].pdf_weight_container; 
            vector<uint16_t> compressed_alpha_s_weights = (*WeightHandle)[weight_index].alpha_s_container; 
            vector<uint16_t> compressed_scale_weights = (*WeightHandle)[weight_index].qcd_scale_container;

            std::vector<float> uncompressed = (*WeightHandle)[weight_index].uncompress( compressed_weights );
            std::vector<float> uncompressed_alpha_s = (*WeightHandle)[weight_index].uncompress( compressed_alpha_s_weights );
            std::vector<float> uncompressed_scale = (*WeightHandle)[weight_index].uncompress( compressed_scale_weights );

            if ( (*WeightHandle)[weight_index].qcd_scale_container.size() == 0 ) {
              for ( unsigned int j = 0 ; j < 9 ; j++ ) 
                event.flashgg_mc_weights.push_back(0.);
            } else{
              for( unsigned int j=0; j<(*WeightHandle)[weight_index].qcd_scale_container.size();j++ )
                event.flashgg_mc_weights.push_back(uncompressed_scale[j]);
            }
            for( unsigned int j=0; j<(*WeightHandle)[weight_index].pdf_weight_container.size();j++ )
              event.flashgg_mc_weights.push_back(uncompressed[j]);
            for( unsigned int j=0; j<(*WeightHandle)[weight_index].alpha_s_container.size();j++ )
              event.flashgg_mc_weights.push_back(uncompressed_alpha_s[j]);

          #ifdef DEBUG_GRINDER
            std::cout << "i handler, N weights = ... " << weight_index << " " << uncompressed.size() << " " << uncompressed_alpha_s.size() << " " << uncompressed_scale.size() << std::endl;
          #endif
          }
        }
      } 

      // iterate over MET =============================================================================================================
      if( METs->size() != 1 ) { std::cout << "WARNING - #MET is not 1" << std::endl;}
      pat::MET const & srcMET  = METs->front();

      // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#ETmiss
      out_met.pt  = srcMET.shiftedPt (pat::MET::NoShift, pat::MET::Type1);
      out_met.phi = srcMET.shiftedPhi(pat::MET::NoShift, pat::MET::Type1);
      out_met.significance = srcMET.metSignificance();
      
      if (not iEvent.isRealData()) {
          out_met.gen_pt  = srcMET.genMET()->pt() ;
          out_met.gen_phi = srcMET.genMET()->phi();

          // https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/DataFormats/PatCandidates/interface/MET.h#L151-L168
          using Var = pat::MET::METUncertainty;
          for (Var const & var : {Var::JetEnUp, Var::JetEnDown, Var::JetResUp, Var::JetResDown, Var::UnclusteredEnUp, Var::UnclusteredEnDown}) {
            out_met.pt_unc_v.push_back(  srcMET.shiftedPt (var, pat::MET::Type1) ); 
            out_met.phi_unc_v.push_back( srcMET.shiftedPhi(var, pat::MET::Type1) ); 
            #ifdef DEBUG_GRINDER
              std::cout << "  met UNCERTANTIE : " << srcMET.shiftedPt (var, pat::MET::Type1) << " " << srcMET.shiftedPhi(var, pat::MET::Type1) << std::endl;
            #endif
          }
      }

      // ================================================================================================================================================ diphotons
      std::vector< edm::Ptr<flashgg::DiPhotonCandidate> > diphotons_sys;
      for( unsigned int i = 0; i < diPhotonTokens_.size(); ++i ){
      #ifdef DEBUG_GRINDER
        std::cout << "  load di-photons at" << i << std::endl;
      #endif

        Handle<View<flashgg::DiPhotonCandidate> > diphotons_vec_tmp;
        iEvent.getByToken( diPhotonTokens_[i], diphotons_vec_tmp );
        if(diphotons_vec_tmp->size() == 0) return;
        diphotons_sys.push_back( diphotons_vec_tmp->ptrAt( 0 ) );
      }

      edm::Ptr<flashgg::DiPhotonCandidate> diphotons_nom = diphotons_sys.at(0);

      // ================================================================================================================================================
      unsigned int jetCollectionIndex = diphotons_nom->jetCollectionIndex(); // always 0 for all diphotons/sys
      edm::Handle<edm::View<flashgg::Jet> > Jets_ ;
      iEvent.getByToken( jetTokens_[jetCollectionIndex], Jets_ );
      #ifdef DEBUG_GRINDER
        std::cout << "  dipho jetCollectionIndex = " << jetCollectionIndex << std::endl;
      #endif
      if( Jets_->size() < 1 ) return;

      // Pile-Up Info ========================================================================================================
      // https://twiki.cern.ch/twiki/bin/view/CMS/Pileup_MC_Information
      // uncertanties calc at the next step https://twiki.cern.ch/twiki/bin/view/CMS/PileupSystematicErrors
      #ifdef DEBUG_GRINDER
        std::cout << "Read Pile-Up Info ... " << std::endl;
      #endif
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

        // auto cache = globalVarsDumper->cache();
        globalVariablesComputer_.update(iEvent);
        auto cache = globalVariablesComputer_.cache();
        event.flashgg_puweight = cache.puweight;
        event.flashgg_nvtx     = cache.nvtx;
        event.flashgg_npu      = cache.npu;
        #ifdef DEBUG_GRINDER
          std::cout << "globalVarsDumper cashe info : " << std::endl;
          std::cout << "  rho = " << cache.rho << " " << *rho << std::endl;
          std::cout << "  nvtx = " << cache.nvtx << std::endl;
          std::cout << "  event = " << cache.event << std::endl;
          std::cout << "  lumi = " << cache.lumi << std::endl;
          std::cout << "  run = " << cache.run << std::endl;
          std::cout << "  npu = " << cache.npu << std::endl;
          std::cout << "  puweight = " << cache.puweight << std::endl;
          std::cout << "  processIndex = " << cache.processIndex << std::endl;
        #endif


        // Gen jets ========================================================================================================
        #ifdef DEBUG_GRINDER
          std::cout << "Gen jets ... " << std::endl;
        #endif
        edm::Handle<edm::View<reco::GenJet>> genJets;
        if(not iEvent.isRealData())
          iEvent.getByToken(genJetToken, genJets);
        
        // jets ===================================================
        #ifdef DEBUG_GRINDER
          std::cout << "Iterate over Jets ... " << std::endl;
        #endif
        for( unsigned int candIndex_outer = 0; candIndex_outer <  Jets_->size() ; candIndex_outer++ ) {
          const flashgg::Jet & j = Jets_->at( candIndex_outer );
          // RAW P4 vector 
          // https://twiki.cern.ch/twiki/bin/view/CMS/TopJME#Jets
          reco::Candidate::LorentzVector const &rawP4 = j.correctedP4("Uncorrected");
          
          jet.pt  = j.pt();
          jet.eta = j.eta();
          jet.phi = j.phi();
          jet.m   = j.mass();
          
          jet.ptRaw  = rawP4.pt();
          jet.etaRaw = rawP4.eta();
          jet.phiRaw = rawP4.phi();
          jet.mRaw   = rawP4.mass();

          jet.charge = j.jetCharge();
          jet.area   = j.jetArea();

          jet.PUJID      = j.puJetIdMVA();
          // pT scale corrections
          // here you must use the CORRECTED jet pt
          // https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Recommendation_for_analysis
          for(int i = jecUnc_v.size()-1; i >= 0; i--){
            JetCorrectionUncertainty * jecUnc_ptr = jecUnc_v[i];
            jecUnc_ptr->setJetEta( j.eta() );
            jecUnc_ptr->setJetPt(  j.pt()  );
            jet.JEC_unc_v_u[i] = jecUnc_ptr->getUncertainty(true);
            jecUnc_ptr->setJetEta( j.eta() );
            jecUnc_ptr->setJetPt(  j.pt()  ); 
            jet.JEC_unc_v_d[i] = jecUnc_ptr->getUncertainty(false);

            #ifdef DEBUG_GRINDER
              std::cout << "JEC uncertanties ... " << jet.JEC_unc_v_u[i] << " " << jet.JEC_unc_v_d[i] << std::endl;
            #endif
          }
          double JEC_uncertanty = std::max( jet.JEC_unc_v_u.at( 0 ), jet.JEC_unc_v_d.at( 0 ) );
          JEC_uncertanty = std::max(JEC_uncertanty, 0.);

          // pT resolution
          // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
          // double JER_uncertanty = 0;
          if (not iEvent.isRealData()) {
            jerResolution_parameters.setJetPt(j.pt()).setJetEta(j.eta());
            jerScaleFactor_parameters.setJetEta(j.eta()).setRho( *rho );

            // jet.resolution = jerResolution.getResolution( jerResolution_parameters );
            jet.resolution = jerResolution.getResolution( {{JME::Binning::JetPt, j.pt()}, {JME::Binning::JetEta, j.eta()}, {JME::Binning::Rho, *rho}} );

            jet.sf         = jerScaleFactor.getScaleFactor(jerScaleFactor_parameters);
            jet.sf_u       = jerScaleFactor.getScaleFactor(jerScaleFactor_parameters, Variation::UP);
            jet.sf_d       = jerScaleFactor.getScaleFactor(jerScaleFactor_parameters, Variation::DOWN);

            // https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Smearing_procedures
            // https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=54#Smearing_procedures
            // https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_18/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L236-L237
            reco::GenJet const * genJet = MatchGenJet( j, genJets, 3 * jet.resolution * j.pt() );
            if (genJet) {
              jet.getJet_pt = genJet->pt();
              double energy_factor = ( rawP4.pt() - genJet->pt() ) / rawP4.pt();
              // JER_uncertanty = std::max( { (jet.sf - 1.) * energy_factor, (jet.sf_u - 1.) * energy_factor, (jet.sf_d - 1.) * energy_factor } );
            } else {
              jet.getJet_pt = -1;
              double max_unc = std::max( { TMath::Abs(jet.sf), TMath::Abs(jet.sf_u), TMath::Abs(jet.sf_d) } );
              // JER_uncertanty = jet.resolution * std::sqrt( std::max(std::pow(max_unc, 2) - 1., 0.) );
            }
          }

          // double jet_maxPt = j.pt() * ( 1. + JEC_uncertanty + JER_uncertanty );
          // if( rawP4.pt()  < cut_jet_pt and jet_maxPt < cut_jet_pt ) continue;
          // if( TMath::Abs(rawP4.eta()) > cut_jet_eta ) continue;
          
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
          jet.isTight = false;
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
          jet.pfDeepCSVJetTags_probc    = j.bDiscriminator("pfDeepCSVJetTags:probc");
          jet.pfDeepCSVJetTags_probudsg = j.bDiscriminator("pfDeepCSVJetTags:probudsg");

          // DeepJet
          jet.pfDeepFlavourJetTags_probb    = j.bDiscriminator("mini_pfDeepFlavourJetTags:probb");
          jet.pfDeepFlavourJetTags_probbb   = j.bDiscriminator("mini_pfDeepFlavourJetTags:probbb");
          jet.pfDeepFlavourJetTags_problepb = j.bDiscriminator("mini_pfDeepFlavourJetTags:problepb");
          jet.pfDeepFlavourJetTags_probc    = j.bDiscriminator("mini_pfDeepFlavourJetTags:probc");
          jet.pfDeepFlavourJetTags_probuds  = j.bDiscriminator("mini_pfDeepFlavourJetTags:probuds");
          jet.pfDeepFlavourJetTags_probg    = j.bDiscriminator("mini_pfDeepFlavourJetTags:probg");

            #ifdef DEBUG_GRINDER
              std::cout << "jet.pfDeepFlavourJetTags_probb = "     << jet.pfDeepFlavourJetTags_probb << std::endl;
              std::cout << "jet.pfDeepFlavourJetTags_probbb = "    << jet.pfDeepFlavourJetTags_probbb << std::endl;
              std::cout << "jet.pfDeepFlavourJetTags_problepb = "  << jet.pfDeepFlavourJetTags_problepb << std::endl;
              std::cout << "jet.pfDeepFlavourJetTags_probc = "     << jet.pfDeepFlavourJetTags_probc << std::endl;
              std::cout << "jet.pfDeepFlavourJetTags_probuds = "   << jet.pfDeepFlavourJetTags_probuds << std::endl;
              std::cout << "jet.pfDeepFlavourJetTags_probg = "     << jet.pfDeepFlavourJetTags_probg << std::endl;
            #endif

          if (not iEvent.isRealData()) {
            jet.hadronFlavour = j.hadronFlavour();
            jet.partonFlavour = j.partonFlavour();
          }
          
          out_jets.emplace_back( jet );
        }
        
        // muons ===================================================                                                            
        // std::vector<edm::Ptr<flashgg::Muon> > goodMuons = selectMuons( muons->ptrs(), dipho, vertices->ptrs(), muonEtaThreshold_, leptonPtThreshold_, muPFIsoSumRelThreshold_, deltaRMuonPhoThreshold_, deltaRMuonPhoThreshold_ );

        // cout << "N good muons = " << allGoodMuons.size() << "/" << muons->ptrs().size() << endl;
        for(auto m : allGoodMuons){
          const Muon & mu = *m;
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
          muon.relIsoPF  = (mu.pfIsolationR04().sumChargedHadronPt + std::max(0., mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5*mu.pfIsolationR04().sumPUPt))/mu.pt();
          muon.relIsoTrk = mu.isolationR03().sumPt/mu.pt();

          muon.diphotons_veto.clear();
          for( auto diphoton_candidate : diphotons_sys ){
            float dRPhoLeadMuon    = deltaR( mu.eta(), mu.phi(), diphoton_candidate->leadingPhoton()->superCluster()->eta(),    diphoton_candidate->leadingPhoton()->superCluster()->phi()    );
            float dRPhoSubLeadMuon = deltaR( mu.eta(), mu.phi(), diphoton_candidate->subLeadingPhoton()->superCluster()->eta(), diphoton_candidate->subLeadingPhoton()->superCluster()->phi() ); 
            bool pass_diphoton_selection = ( dRPhoLeadMuon < deltaRMuonPhoThreshold_ || dRPhoSubLeadMuon < deltaRMuonPhoThreshold_ ); 
            muon.diphotons_veto.push_back( pass_diphoton_selection );
          }
          
          out_muons.emplace_back( muon );
        }
        
        // iterate over electrons ========================================================================================================
        // Electrons 
        //std::vector<edm::Ptr<Electron> > goodElectrons = selectStdElectrons( electrons->ptrs(), dipho, vertices->ptrs(), leptonPtThreshold_, electronEtaThresholds_, useElectronMVARecipe_, useElectronLooseID_, deltaRPhoElectronThreshold_, DeltaRTrkElec_, deltaMassElectronZThreshold_, *rho, iEvent.isRealData() );
        // std::vector<edm::Ptr<Electron> > goodElectrons = electrons->ptrs();

        // cout << "N good electrons = " << allGoodElectrons.size() << "/" << electrons->ptrs().size() << endl;
        for( auto e : allGoodElectrons ){
          const Electron & el = *e;
          // cout << "e pt, eta, phi, charge = " << el.pt() << " " << el.eta() << " " << el.phi() << " " << el.charge() << endl;
          electron.pt     = el.pt();
          electron.eta    = el.eta();
          electron.phi    = el.phi();
          electron.charge = el.charge();

          electron.diphotons_veto.clear();
          for( auto diphoton_candidate : diphotons_sys ){
            bool rho_veto = phoVeto(e, diphoton_candidate, deltaRPhoElectronThreshold_, DeltaRTrkElec_, deltaMassElectronZThreshold_);
            electron.diphotons_veto.push_back( rho_veto );
          }
          
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

          #ifdef DEBUG_GRINDER
            cout << "e loose,med,tight = " << electron.isLoose << " " << electron.isMedium << " " << electron.isTight << endl;
            cout << "e CH,NH,PH,PU = " << electron.sumChargedHadronPt << " " << electron.sumNeutralHadronEt << " " << electron.sumPhotonEt << " " << electron.sumPUPt << endl;
          #endif

          if(not iEvent.isRealData()){
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
            // electron.energyScaleEtUp       = el.userFloat("energyScaleEtUp"); //FIXME ???
            // electron.energyScaleEtDown     = el.userFloat("energyScaleEtDown"); //FIXME ???
            electron.energySigmaUp         = el.userFloat("energySigmaUp");
            electron.energySigmaDown       = el.userFloat("energySigmaDown");
            electron.energySigmaPhiUp      = el.userFloat("energySigmaPhiUp");
            electron.energySigmaPhiDown    = el.userFloat("energySigmaPhiDown");
            electron.energySigmaRhoUp      = el.userFloat("energySigmaRhoUp");
            electron.energySigmaRhoDown    = el.userFloat("energySigmaRhoDown");
            // cout << "some sys ... " << electron.ecalTrkEnergyPreCorr << " " << electron.ecalTrkEnergyPostCorr << " " << electron.energyScaleValue << " " << electron.energyScaleUp << " " << electron.energySigmaPhiUp << endl;

            // 2017 year : 
            // Requested UserFloat energyScaleEtUp is not available! Possible UserFloats are: 
            // ElectronMVAEstimatorRun2Fall17IsoV1Values ElectronMVAEstimatorRun2Fall17NoIsoV1Values ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values ElectronMVAEstimatorRun2Spring15Trig25nsV1Values ElectronMVAEstimatorRun2Spring15Trig50nsV1Values ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values ElectronMVAEstimatorRun2Spring16HZZV1Values ecalEnergyErrPostCorr ecalEnergyErrPreCorr ecalEnergyPostCorr ecalEnergyPreCorr ecalTrkEnergyErrPostCorr ecalTrkEnergyErrPreCorr ecalTrkEnergyPostCorr ecalTrkEnergyPreCorr energyScaleDown energyScaleGainDown energyScaleGainUp energyScaleStatDown energyScaleStatUp energyScaleSystDown energyScaleSystUp energyScaleUp energyScaleValue energySigmaDown energySigmaPhiDown energySigmaPhiUp energySigmaRhoDown energySigmaRhoUp energySigmaUp energySigmaValue energySmearNrSigma heepTrkPtIso 
          }
          
          out_electrons.emplace_back( electron );
        }
        
        // photons ===================================================
        edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( 0 );   
        for( auto diphoton_candidate : diphotons_sys ){
          edm::Ptr<reco::Vertex> diphoton_vertex = diphoton_candidate->vtx();
          
          const flashgg::Photon* leadPho    = diphoton_candidate->leadingPhoton();
          const flashgg::Photon* subleadPho = diphoton_candidate->subLeadingPhoton();
          // bool passMVAs = checkPassMVAs(leadPho, subleadPho, diphoton_vertex);

          #ifdef DEBUG_GRINDER
            auto dc_p4 = diphoton_candidate->p4();
            std::cout << "dipho mass vs ph+ph mass = " << dc_p4.mass() << " " << (leadPho->p4() + subleadPho->p4()).mass() << endl;
          #endif
          
          vector<const flashgg::Photon*> phs = { leadPho, subleadPho };
          for(const flashgg::Photon* ph : phs){
            photon.pt  = ph->pt();
            photon.eta = ph->eta();
            photon.phi = ph->phi();
            photon.mva_value = ph->phoIdMvaDWrtVtx( diphoton_vertex );
            out_photons.emplace_back( photon );

            #ifdef DEBUG_GRINDER
              std::cout << "Photon ... " << photon.pt << " " << photon.eta << " " << photon.phi << " " << photon.mva_value << std::endl;
              std::cout << "n weights = " << ph->weightListEnd() - ph->weightListBegin() << std::endl;
              for(auto it = ph->weightListBegin(); it != ph->weightListEnd(); ++it){
                std::cout << "with weights = " << ph->weight( *it ) << " " << (*it) << endl;
              }
            #endif
          }

            #ifdef DEBUG_GRINDER
              std::cout << "n weights = " << diphoton_candidate->weightListEnd() - diphoton_candidate->weightListBegin() << std::endl;
              for(auto it = diphoton_candidate->weightListBegin(); it != diphoton_candidate->weightListEnd(); ++it){
                std::cout << "with weights = " << diphoton_candidate->weight( *it ) << " " << (*it) << endl;
              }
            #endif
        }

        // weights ===================================================
        for(auto it = diphotons_nom->weightListBegin(); it != diphotons_nom->weightListEnd(); ++it){
          // std::cout << "with weights = " << diphotons_nom->weight( *it ) << " " << (*it) << endl;
          event.flashgg_diphoton_weights.push_back( diphotons_nom->weight( *it ) ); // Central
        }

        outTree->Fill();
        // cout << "outTree->Fill();" << endl;
        return;
      // ========================================================================================================
    }

    // ------------ method called when starting to processes a run  ------------
    void GrinderHHWWggTagProducer::beginRun( edm::Run const & run, edm::EventSetup const & setup ){
      #ifdef DEBUG_GRINDER
        std::cout << "Grinder::beginJob ... " << std::endl;
      #endif

      // Construct an object to obtain JEC uncertainty
      // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections?rev=137#JetCorUncertainties
      jecUnc_v.clear();
      edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
      setup.get<JetCorrectionsRecord>().get("AK4PF", JetCorParColl); 
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

      //#ifdef DEBUG_GRINDER
      //  JME::JetResolution::get(setup,            "AK4PFchs_pt" ).dump();
      //#endif
    }

    void GrinderHHWWggTagProducer::endRun(edm::Run const& run, edm::EventSetup const& setup){
      std::cout << "Grinder ... lumiWeight = " << eventMeta.sumWeights << std::endl;
      outTreeMeta->Fill();
    };

    // ------------ method called once each job just before starting event loop  ------------
    void GrinderHHWWggTagProducer::beginJob(){
    }

    // ------------ method called once each job just after ending the event loop  ------------
    void GrinderHHWWggTagProducer::endJob(){
    }

  } 

  typedef flashgg::GrinderHHWWggTagProducer GrinderFlashggHHWWggTagProducer;
  DEFINE_FWK_MODULE( GrinderFlashggHHWWggTagProducer );
  
  
  
  
  
  
