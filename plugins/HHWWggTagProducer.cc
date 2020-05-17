// Abe Tishelman-Charny
// November 2019
// Derived from HH->WWgg event dumper and HH->bbgg tagger 

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

#include "flashgg/MicroAOD/interface/CutBasedDiPhotonObjectSelector.h"

#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
#include "flashgg/DataFormats/interface/WHLeptonicTag.h"

#include "flashgg/Taggers/interface/LeptonSelection.h"

#include <vector>
#include <algorithm>
#include "TGraph.h"
#include "TLorentzVector.h"

// Grinder
#include "Analysis/GRINDER/interface/Event.hh"
using namespace grinder;
#include <TTree.h>

using namespace std;
using namespace edm;

namespace flashgg {
  class GrinderHHWWggTagProducer : public EDProducer
  {
  public:
    //---typedef
    typedef math::XYZTLorentzVector LorentzVector;

    //---ctors
    // HHWWggTagProducer();
    GrinderHHWWggTagProducer( const ParameterSet & );

    //---Outtree 
    edm::Service<TFileService> fs;
    
    // TH1F* indexes;
    // TH1F* btags;

  private:
    double genTotalWeight;
    bool checkPassMVAs(const flashgg::Photon*& leading_photon, const flashgg::Photon*& subleading_photon, edm::Ptr<reco::Vertex>& diphoton_vertex);
    void produce( edm::Event &, const EventSetup & ) override;
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

    bool hasGoodElec = false;
    bool hasGoodMuons = false;

    vector<double> nonTrigMVAThresholds_;
    vector<double> nonTrigMVAEtaCuts_;

    double electronIsoThreshold_;
    double electronNumOfHitsThreshold_;
    vector<double> electronEtaThresholds_;
    bool useElectronMVARecipe_;
    bool useElectronLooseID_;
    // string bTag_;
    double btagThresh_;
    bool doHHWWggTagCutFlowAnalysis_;


    edm::InputTag genInfo_;
    edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
    
    
    // GRINDER PART ========================================>
      edm::Service<TFileService> fileService;
      TTree *outTree, *outTreeMeta;
      grinder::Event  out_event;
      grinder::Event *event_ptr;
      grinder::EventMetadata  out_eventMeta;
      grinder::EventMetadata *eventMeta_ptr;

      grinder::Jet jet;
      grinder::Muon muon;
      grinder::Photon photon;
      grinder::Electron electron;
      grinder::MET out_met;

      std::vector<grinder::Jet>      out_jets;
      std::vector<grinder::Muon>     out_muons;
      std::vector<grinder::Photon>   out_photons;
      std::vector<grinder::Electron> out_electrons;

      std::vector<grinder::Jet>      *jets_ptr;
      std::vector<grinder::Muon>     *muons_ptr;
      std::vector<grinder::Photon>   *photons_ptr;
      std::vector<grinder::Electron> *electrons_ptr;
      grinder::MET                   *met_ptr;
    // ============ ========================================>
  };

  //---constructors
  // HHWWggTagProducer::HHWWggTagProducer( ):
  // photonToken_(),
  // diphotonToken_()
  // genParticleToken_(),
  // electronToken_(),
  // muonToken_(),
  // METToken_(),
  // cc_( consumesCollector() )
  // // idSelector_( ParameterSet(), cc_ )

  // {}

    //---standard
    GrinderHHWWggTagProducer::GrinderHHWWggTagProducer( const ParameterSet & pSet):
    photonToken_( consumes<View<Photon> >( pSet.getParameter<InputTag> ( "PhotonTag" ) ) ),
    diphotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( pSet.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
    vertexToken_( consumes<View<reco::Vertex> >( pSet.getParameter<InputTag> ( "VertexTag" ) ) ),
    genParticleToken_( consumes<View<reco::GenParticle> >( pSet.getParameter<InputTag> ( "GenParticleTag" ) ) ),
    electronToken_( consumes<View<Electron> >( pSet.getParameter<InputTag> ( "ElectronTag" ) ) ), 
    muonToken_( consumes<View<Muon> >( pSet.getParameter<InputTag> ( "MuonTag" ) ) ),
    METToken_( consumes<View<Met> >( pSet.getParameter<InputTag> ( "METTag" ) ) ),
    mvaResultToken_( consumes<View<flashgg::DiPhotonMVAResult> >( pSet.getParameter<InputTag> ( "MVAResultTag" ) ) ),
    rhoTag_( consumes<double>( pSet.getParameter<InputTag>( "rhoTag" ) ) ),
    triggerRECO_( consumes<edm::TriggerResults>(pSet.getParameter<InputTag>("RECOfilters") ) ),
    triggerPAT_( consumes<edm::TriggerResults>(pSet.getParameter<InputTag>("PATfilters") ) ),
    triggerFLASHggMicroAOD_( consumes<edm::TriggerResults>( pSet.getParameter<InputTag>("FLASHfilters") ) ),
    systLabel_( pSet.getParameter<string> ( "SystLabel" ) ),
    cc_( consumesCollector() ), // need absence of comma on last entry 
    globalVariablesComputer_(pSet.getParameter<edm::ParameterSet>("globalVariables"), cc_)
    // idSelector_( pSet.getParameter<ParameterSet> ( "idSelection" ), cc_ )

    {

      inputDiPhotonName_= pSet.getParameter<std::string > ( "DiPhotonName" );
      inputDiPhotonSuffixes_= pSet.getParameter<std::vector<std::string> > ( "DiPhotonSuffixes" );
      std::vector<edm::InputTag>  diPhotonTags;
      for (auto & suffix : inputDiPhotonSuffixes_){ 
          systematicsLabels.push_back(suffix);
          std::string inputName = inputDiPhotonName_;
          inputName.append(suffix);
          if (!suffix.empty()) diPhotonTags.push_back(edm::InputTag(inputName));
          else  diPhotonTags.push_back(edm::InputTag(inputDiPhotonName_));
      }
      for( auto & tag : diPhotonTags ) { diPhotonTokens_.push_back( consumes<edm::View<flashgg::DiPhotonCandidate> >( tag ) ); }

      bool breg = 0;

      inputJetsName_= pSet.getParameter<std::string> ( "JetsName" );
      inputJetsCollSize_= pSet.getParameter<unsigned int> ( "JetsCollSize" );
      inputJetsSuffixes_= pSet.getParameter<std::vector<std::string> > ( "JetsSuffixes" );
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
        auto jetTags = pSet.getParameter<std::vector<edm::InputTag> > ( "JetTags" ); 
        for( auto & tag : jetTags ) { jetTokens_.push_back( consumes<edm::View<flashgg::Jet> >( tag ) ); }
      }


      genInfo_ = pSet.getUntrackedParameter<edm::InputTag>( "genInfo", edm::InputTag("generator") );
      genInfoToken_ = consumes<GenEventInfoProduct>( genInfo_ );
      // numDiphoCand = fs->make<TH1F> ("numDiphoCand","numDiphoCand",10,0,10); 
      // diphoton_idx_h = fs->make<TH1F> ("diphoton_idx_h","diphoton_idx_h",20,0,20); 
      // diPhotons_size_h = fs->make<TH1F> ("diPhotons_size_h","diPhotons_size_h",20,0,20); 

      // indexes = fs->make<TH1F> ("indexes","indexes",5,0,5);
      // btags = fs->make<TH1F> ("btags","btags",100,0,1);

      // numEvents = fs->make<TH1F> ("numEvents","numEvents",1,0,10);

      // gen_weights = fs->make<TH1F> ("gen_weights","gen_weights",1000,-2,2);
      // vars = fs->make<TH1F> ("vars","vars",10,0,10);
      // cutFlow = fs->make<TH1F> ("cutFlow","Cut Flow",10,0,10);
      // WTags = fs->make<TH1F> ("WTags","W Tags",3,0,3);

      leptonPtThreshold_ = pSet.getParameter<double>( "leptonPtThreshold");
      muonEtaThreshold_ = pSet.getParameter<double>( "muonEtaThreshold");
      leadPhoOverMassThreshold_ = pSet.getParameter<double>( "leadPhoOverMassThreshold");
      subleadPhoOverMassThreshold_ = pSet.getParameter<double>( "subleadPhoOverMassThreshold");
      MVAThreshold_ = pSet.getParameter<double>( "MVAThreshold");
      deltaRMuonPhoThreshold_ = pSet.getParameter<double>( "deltaRMuonPhoThreshold");
      jetsNumberThreshold_ = pSet.getParameter<double>( "jetsNumberThreshold");
      jetPtThreshold_ = pSet.getParameter<double>( "jetPtThreshold");
      jetEtaThreshold_ = pSet.getParameter<double>( "jetEtaThreshold");
      muPFIsoSumRelThreshold_ = pSet.getParameter<double>( "muPFIsoSumRelThreshold");
      PhoMVAThreshold_ = pSet.getParameter<double>( "PhoMVAThreshold");
      METThreshold_ = pSet.getParameter<double>( "METThreshold");
      useVertex0only_              = pSet.getParameter<bool>("useVertex0only");
      deltaRJetMuonThreshold_ = pSet.getParameter<double>( "deltaRJetMuonThreshold");
      deltaRPhoLeadJet_ = pSet.getParameter<double>( "deltaRPhoLeadJet");
      deltaRPhoSubLeadJet_ = pSet.getParameter<double>( "deltaRPhoSubLeadJet");

      DeltaRTrkElec_ = pSet.getParameter<double>( "DeltaRTrkElec");
      TransverseImpactParam_ = pSet.getParameter<double>( "TransverseImpactParam");
      LongitudinalImpactParam_ = pSet.getParameter<double>( "LongitudinalImpactParam");

      deltaRPhoElectronThreshold_ = pSet.getParameter<double>( "deltaRPhoElectronThreshold");
      deltaMassElectronZThreshold_ = pSet.getParameter<double>( "deltaMassElectronZThreshold");

      nonTrigMVAThresholds_ =  pSet.getParameter<vector<double > >( "nonTrigMVAThresholds");
      nonTrigMVAEtaCuts_ =  pSet.getParameter<vector<double > >( "nonTrigMVAEtaCuts");
      electronIsoThreshold_ = pSet.getParameter<double>( "electronIsoThreshold");
      electronNumOfHitsThreshold_ = pSet.getParameter<double>( "electronNumOfHitsThreshold");
      electronEtaThresholds_ = pSet.getParameter<vector<double > >( "electronEtaThresholds");
      useElectronMVARecipe_=pSet.getParameter<bool>("useElectronMVARecipe");
      useElectronLooseID_=pSet.getParameter<bool>("useElectronLooseID");
      // bTag_ = pSet.getParameter<string> ( "bTag");
      btagThresh_ = pSet.getParameter<double>( "btagThresh");
      doHHWWggTagCutFlowAnalysis_ = pSet.getParameter<bool>( "doHHWWggTagCutFlowAnalysis");

      produces<vector<HHWWggTag>>();
      // for (auto & systname : systematicsLabels) { // to deal with systematics in producer 
      //     produces<vector<HHWWggTag>>(systname);
      // }
      produces<vector<TagTruthBase>>();
      
      // ============ ========================================> GRIDNER PART
        outTree     = fileService->make<TTree>("Events", "Events");
        outTreeMeta = fileService->make<TTree>("EventsMeta", "EventsMeta");

        event_ptr     = &out_event;
        jets_ptr      = &out_jets;
        muons_ptr     = &out_muons;
        photons_ptr   = &out_photons;
        electrons_ptr = &out_electrons;
        met_ptr       = &out_met;

        outTree->Branch("Event",     &event_ptr);
        outTree->Branch("Photons",   &photons_ptr);
        outTree->Branch("Electrons", &electrons_ptr);
        outTree->Branch("Muons",     &muons_ptr);
        outTree->Branch("Jets",      &jets_ptr);
        outTree->Branch("MET",       &met_ptr);

        eventMeta_ptr = &out_eventMeta;

        outTreeMeta->Branch("EventMeta", &eventMeta_ptr);
      // ============ ========================================>
    }

    bool GrinderHHWWggTagProducer::checkPassMVAs( const flashgg::Photon*& leading_photon, const flashgg::Photon*& subleading_photon, edm::Ptr<reco::Vertex>& diphoton_vertex){
      bool debug_mva = 0; 

      // MVA Check variables 
      double lp_mva_thresh = 0.07;
      double slp_mva_thresh = -0.03;

      bool lead_pass_TightPhoID = 0, sublead_pass_TightPhoID = 0;
      double lp_Hgg_MVA = -99, slp_Hgg_MVA = -99; 
      double leading_pho_eta = -99, sub_leading_pho_eta = -99;

      // Get MVA values wrt diphoton vertex
      lp_Hgg_MVA = leading_photon->phoIdMvaDWrtVtx( diphoton_vertex ); 
      slp_Hgg_MVA = subleading_photon->phoIdMvaDWrtVtx( diphoton_vertex ); 

      // Get eta values
      leading_pho_eta = leading_photon->p4().eta();
      sub_leading_pho_eta = subleading_photon->p4().eta();

      // Debug
      if(debug_mva){
        cout << "leading mva: " << lp_Hgg_MVA << endl;
        cout << "subleading mva: " << slp_Hgg_MVA << endl;
      }
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
        if(debug_mva) cout << "PASS MVA Selections" << endl;
        return 1;
    }

    else{
      if(debug_mva) cout << "FAIL MVA Selections" << endl;
      return 0; 
    }

    }

    void GrinderHHWWggTagProducer::produce( edm::Event &event, const EventSetup & )
    {
      cout << "HHWWggTagProducer::produce();" << endl;
      
      // cout << "[HHWWggTagProducer.cc] - Beginning of HHWWggTagProducer::produce" << endl;

      // update global variables
      // globalVariablesComputer_.update(event);

      // Get particle objects
      event.getByToken( photonToken_, photons );
      event.getByToken( diphotonToken_, diphotons );
      // event.getByToken( genParticleToken_, genParticle );
      event.getByToken( electronToken_, electrons );
      event.getByToken( muonToken_, muons );
      event.getByToken( METToken_, METs );
      event.getByToken( mvaResultToken_, mvaResults );
      event.getByToken( vertexToken_, vertices );
      event.getByToken( rhoTag_, rho);

      double rho_    = *rho;

      // Set cut booleans
      // std::vector<double> Cut_Results = {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; // Cut_Results[i] = 1: Event Passed Cut i 
      std::vector<double> Cut_Variables(20,0.0); // Cut_Results[i] = 1.0: Event Passed Cut i 
      // std::vector<double> Vertex_Variables(20,0.0); // Cut_Results[i] = 1.0: Event Passed Cut i 

      // Cut Variables 
      // double has_PS_Dipho = 0, pass_METfilters = 0, dipho_vertex_is_zero = 0, pass_leadPhoOverMassThreshold = 0, pass_subleadPhoOverMassThreshold = 0,
      //   pass_LeadPhoton_MVA = 0, pass_SubLeadPhoton_MVA = 0, pass_dipho_MVA = 0, number_passed_jetid = 0;
      // double dipho_vertex_is_zero = -999;
      // double SLW_Tag = 0.; // Semi-Leptonic W Tag  
      // double FLW_Tag = 0.; // Fully-Leptonic W Tag
      // double FHW_Tag = 0.; // Fully-Hadronic W Tag 
      // bool PS_dipho_tag = 0; // preselected diphoton 

      //---output collection
      // std::unique_ptr<vector<HHWWggCandidate> > HHWWggColl_( new vector<HHWWggCandidate> );
      // std::unique_ptr<vector<HHWWggTag> > tags( new vector<HHWWggTag> );
      // int n_METs = METs->size(); // Should be 1, but using as a way to obtain met four vector 
      int n_good_electrons = 0;
      int n_good_muons = 0;
      int n_good_leptons = 0;
      int n_good_jets = 0;
      bool hasHighbTag = 0;
      float btagVal = 0;
      // double dipho_MVA = -99;
      // double lead_pho_Hgg_MVA = -99, sublead_pho_Hgg_MVA = -99;
      // double CMS_hgg_mass = -99;
      // float bDiscriminatorValue = -2.;

      bool passMVAs = 0; // True if leading and subleading photons pass MVA selections 

      // Saved Objects after selections
      std::vector<flashgg::Jet> tagJets_;
      std::vector<flashgg::Muon> goodMuons_;
      std::vector<flashgg::Electron> goodElectrons_; 
      std::vector<flashgg::Met> theMET_;
      std::vector<flashgg::DiPhotonCandidate> diphoVector_;
      reco::GenParticle::Point genVertex;

      std::unique_ptr<vector<TagTruthBase> > truths( new vector<TagTruthBase> );
      edm::RefProd<vector<TagTruthBase> > rTagTruth = event.getRefBeforePut<vector<TagTruthBase> >();


      if (diphotons->size() == 0) return;
      edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( 0 );
          // Electrons 
          std::vector<edm::Ptr<Electron> > goodElectrons = selectStdElectrons( electrons->ptrs(), dipho, vertices->ptrs(), leptonPtThreshold_, electronEtaThresholds_,
                                                                            useElectronMVARecipe_,useElectronLooseID_,
                                                                            deltaRPhoElectronThreshold_,DeltaRTrkElec_,deltaMassElectronZThreshold_,
                                                                            rho_, event.isRealData() );
          // Muons                                                                   
          std::vector<edm::Ptr<flashgg::Muon> > goodMuons = selectMuons( muons->ptrs(), dipho, vertices->ptrs(), muonEtaThreshold_, leptonPtThreshold_,
          muPFIsoSumRelThreshold_, deltaRMuonPhoThreshold_, deltaRMuonPhoThreshold_ );
      
      unsigned int jetCollectionIndex = diphotons->at( 0 ).jetCollectionIndex();
      edm::Handle<edm::View<flashgg::Jet> > Jets_ ;
      event.getByToken( jetTokens_[jetCollectionIndex], Jets_ ); // testing 
      // Remove Info from previous events ========================================================================================================
        out_jets.clear();
        out_muons.clear();
        out_photons.clear();
        out_electrons.clear();
        out_event.weights.clear();
        out_event.ps_weights.clear();
        out_met.pt_unc_v.clear();
        out_met.phi_unc_v.clear();
        
        // jets ===================================================
        for( unsigned int candIndex_outer = 0; candIndex_outer <  Jets_->size() ; candIndex_outer++ ) {
          const edm::Ptr<flashgg::Jet> & j = Jets_->at( candIndex_outer );

          jet.pt  = j.pt();
          jet.eta = j.eta();
          jet.phi = j.phi();
          jet.m   = j.mass();
          
          jet.isTight = true; // all flashgg::Jet are tight?
          
          // b-tagging
          // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#B_tagging
          // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
          // DeepCSV
          jet.pfDeepCSVJetTags_probb    = j.bDiscriminator("pfDeepCSVJetTags:probb");
          jet.pfDeepCSVJetTags_probbb   = j.bDiscriminator("pfDeepCSVJetTags:probbb");
          jet.pfDeepCSVJetTags_probc    = j.bDiscriminator("pfDeepCSVJetTags:probc");
          jet.pfDeepCSVJetTags_probudsg = j.bDiscriminator("pfDeepCSVJetTags:probudsg");

          // DeepJet FIXME not working for b tag
          jet.pfDeepFlavourJetTags_probb    = j.bDiscriminator("pfDeepFlavourJetTags:probb");
          jet.pfDeepFlavourJetTags_probbb   = j.bDiscriminator("pfDeepFlavourJetTags:probbb");
          jet.pfDeepFlavourJetTags_problepb = j.bDiscriminator("pfDeepFlavourJetTags:problepb");
          jet.pfDeepFlavourJetTags_probc    = j.bDiscriminator("pfDeepFlavourJetTags:probc");
          jet.pfDeepFlavourJetTags_probuds  = j.bDiscriminator("pfDeepFlavourJetTags:probuds");
          jet.pfDeepFlavourJetTags_probg    = j.bDiscriminator("pfDeepFlavourJetTags:probg");
          
          out_jets.emplace_back( jet );
        }
        
        // muons ===================================================
        for(auto m : goodMuons){
          muon.pt     = m->pt();
          muon.eta    = m->eta();
          muon.phi    = m->phi();
          muon.charge = m->charge();
          
          out_muons.emplace_back( muon );
        }
        
        // electrons ===================================================
        for(auto e : goodElectrons){
          electron.pt     = e->pt();
          electron.eta    = e->eta();
          electron.phi    = e->phi();
          electron.charge = e->charge();
          
          out_electrons.emplace_back( electron );
        }
        
        // photons ===================================================
          edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( 0 );   
          edm::Ptr<reco::Vertex> diphoton_vertex = dipho->vtx();
          
          const flashgg::Photon* leadPho = dipho->leadingPhoton();
          const flashgg::Photon* subleadPho = dipho->subLeadingPhoton();
          passMVAs = checkPassMVAs(leadPho, subleadPho, diphoton_vertex);
          
          vector<const flashgg::Photon*> phs = { leadPho, subleadPho };
          for(const flashgg::Photon* ph : phs){
            photon.pt  = ph->pt();
            photon.eta = ph->eta();
            photon.phi = ph->phi();
            photon.mva_value = passMVAs;
            out_photons.emplace_back( photon );
          }
        
        // weights ===================================================
        
        // met ===================================================
        if( METs->size() != 1 ) { std::cout << "WARNING - #MET is not 1" << std::endl;}
        Ptr<flashgg::Met> theMET = METs->ptrAt( 0 );
        
        outTree->Fill();
        cout << "outTree->Fill();" << endl;
        return;
      // ========================================================================================================
    } 

  } 

  typedef flashgg::GrinderHHWWggTagProducer GrinderFlashggHHWWggTagProducer;
  DEFINE_FWK_MODULE( GrinderFlashggHHWWggTagProducer );
  
  
  
  
  
  
