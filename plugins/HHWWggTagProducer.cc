// -*- C++ -*-

#define DEBUG_GRINDER 1

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

    GrinderHHWWggTagProducer( const ParameterSet & );

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
    bool doHHWWggTagCutFlowAnalysis_;


    edm::InputTag genInfo_;
    edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
    
    
    // GRINDER PART ========================================>
      std::string era_label;

      edm::Service<TFileService> fileService;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puSummaryToken;

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
      
      // Jets
      // ID
      double NHF, NEMF, CHF, MUF, CEMF, NumConst, NumNeutralParticles, CHM;
      // Electrons
      std::string electron_loose_id_token, electron_medium_id_token, electron_tight_id_token;
      const reco::GsfElectron::PflowIsolationVariables * pfIso;
    // ============ ========================================>
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
      globalVariablesComputer_(iConfig.getParameter<edm::ParameterSet>("globalVariables"), cc_)
      // idSelector_( iConfig.getParameter<ParameterSet> ( "idSelection" ), cc_ )
    {
      puSummaryToken        = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("puSummaryToken_token"));

      inputDiPhotonName_= iConfig.getParameter<std::string > ( "DiPhotonName" );
      inputDiPhotonSuffixes_= iConfig.getParameter<std::vector<std::string> > ( "DiPhotonSuffixes" );
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

      produces<vector<HHWWggTag>>();
      // for (auto & systname : systematicsLabels) { // to deal with systematics in producer 
      //     produces<vector<HHWWggTag>>(systname);
      // }
      produces<vector<TagTruthBase>>();
      
      // ============ ========================================> GRIDNER PART
        era_label             = iConfig.getParameter<std::string>("era_label");
        out_eventMeta.is_data = iConfig.getParameter<bool>("is_data");

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
        
        // read Electron options
        // electronToken = consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons_token"));
        electronToken_ = consumes<View<Electron> >( iConfig.getParameter<InputTag> ( "ElectronTag" ) );
        electron_loose_id_token  = iConfig.getParameter<std::string>("electron_loose_id_token");
        electron_medium_id_token = iConfig.getParameter<std::string>("electron_medium_id_token");
        electron_tight_id_token  = iConfig.getParameter<std::string>("electron_tight_id_token");
      // ============ ========================================>
    }

    bool GrinderHHWWggTagProducer::checkPassMVAs( const flashgg::Photon*& leading_photon, const flashgg::Photon*& subleading_photon, edm::Ptr<reco::Vertex>& diphoton_vertex){
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

    void GrinderHHWWggTagProducer::produce( edm::Event &iEvent, const EventSetup & ) {
      /*
          TODO:
           V Remove Info from previous events
           - Set Event Info
           - weights
           - Read primary vertices collection
           V Pile-Up Info
           V- electrons
           V- muons
           - photons
           - genjets
           - jets
           V- met
      */

      #ifdef DEBUG_GRINDER
        std::cout << "GrinderHHWWggTagProducer::produce() ... " << std::endl;
      #endif
      // Remove Info from previous events ========================================================================================================
      out_jets.clear();
      out_muons.clear();
      out_photons.clear();
      out_electrons.clear();
      out_event.weights.clear();
      out_event.ps_weights.clear();
      out_met.pt_unc_v.clear();
      out_met.phi_unc_v.clear();

      // Set Event Info ============================================================================================================== TODO

      // Get particle objects ======================================================================================================== 
      iEvent.getByToken( photonToken_, photons );
      iEvent.getByToken( diphotonToken_, diphotons );
      iEvent.getByToken( electronToken_, electrons );
      iEvent.getByToken( muonToken_, muons );
      iEvent.getByToken( METToken_, METs );
      iEvent.getByToken( mvaResultToken_, mvaResults );
      iEvent.getByToken( vertexToken_, vertices );
      iEvent.getByToken( rhoTag_, rho);

      #ifdef DEBUG_GRINDER
        std::cout << "  event content :" << std::endl;
        std::cout << "    N photons = " << photons->size() << std::endl;
        std::cout << "    N diphotons = " << diphotons->size() << std::endl;
        std::cout << "    N electrons = " << electrons->size() << std::endl;
        std::cout << "    N muons = " << muons->size() << std::endl;
        std::cout << "    rho = " << *rho << std::endl;
      #endif

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
          }
      }
      out_met.Flag_goodVertices                        =  false ;
      out_met.Flag_globalSuperTightHalo2016Filter      =  false ;
      out_met.Flag_HBHENoiseFilter                     =  false ;
      out_met.Flag_HBHENoiseIsoFilter                  =  false ;
      out_met.Flag_EcalDeadCellTriggerPrimitiveFilter  =  false ;
      out_met.Flag_BadPFMuonFilter                     =  false ;
      out_met.Flag_BadChargedCandidateFilter           =  false ;
      out_met.Flag_eeBadScFilter                       =  false ;
      out_met.Flag_ecalBadCalibReducedMINIAODFilter    =  false ;

      // ================================================================================================================================================
      // FIXME count events 
      if (diphotons->size() == 0) return;
      edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( 0 );
      
      unsigned int jetCollectionIndex = dipho->jetCollectionIndex();
      edm::Handle<edm::View<flashgg::Jet> > Jets_ ;
      iEvent.getByToken( jetTokens_[jetCollectionIndex], Jets_ );
      #ifdef DEBUG_GRINDER
        std::cout << "  dipho jetCollectionIndex = " << jetCollectionIndex << std::endl;
      #endif

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
        out_event.TrueMCNumInteractions = -404;
        for(PVI = puSummary->begin(); PVI != puSummary->end(); ++PVI) {
          if( PVI->getBunchCrossing() != 0 ) continue;
          out_event.DicedMCNumInteractions = PVI->getPU_NumInteractions();
          out_event.TrueMCNumInteractions  = PVI->getTrueNumInteractions();
          break;
        }
      }
        
        // jets ===================================================
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
        std::vector<edm::Ptr<flashgg::Muon> > goodMuons = selectMuons( muons->ptrs(), dipho, vertices->ptrs(), muonEtaThreshold_, leptonPtThreshold_, muPFIsoSumRelThreshold_, deltaRMuonPhoThreshold_, deltaRMuonPhoThreshold_ );
        cout << "N good muons = " << goodMuons.size() << "/" << muons->ptrs().size() << endl;
        for(auto m : goodMuons){
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
          
          out_muons.emplace_back( muon );
        }
        
        // iterate over electrons ========================================================================================================
        // Electrons 
        // std::vector<edm::Ptr<Electron> > goodElectrons = selectStdElectrons( electrons->ptrs(), dipho, vertices->ptrs(), leptonPtThreshold_, electronEtaThresholds_, useElectronMVARecipe_, useElectronLooseID_, deltaRPhoElectronThreshold_, DeltaRTrkElec_, deltaMassElectronZThreshold_, *rho, iEvent.isRealData() );
        std::vector<edm::Ptr<Electron> > goodElectrons = electrons->ptrs();

        cout << "N good electrons = " << goodElectrons.size() << "/" << electrons->ptrs().size() << endl;
        for(auto e : goodElectrons){
          const Electron & el = *e;
          cout << "e pt, eta, phi, charge = " << el.pt() << " " << el.eta() << " " << el.phi() << " " << el.charge() << endl;
          electron.pt     = el.pt();
          electron.eta    = el.eta();
          electron.phi    = el.phi();
          electron.charge = el.charge();
          
          // cut based ID
          // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#Accessing_ID_result
          electron.isLoose  = el.electronID( electron_loose_id_token  );
          electron.isMedium = el.electronID( electron_medium_id_token );
          electron.isTight  = el.electronID( electron_tight_id_token  );
          cout << "e loose,med,tight = " << electron.isLoose << " " << electron.isMedium << " " << electron.isTight << endl;
          
          // ISO https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolationRun2
          pfIso = & ( el.pfIsolationVariables() );
          electron.sumChargedHadronPt = pfIso->sumChargedHadronPt;
          electron.sumNeutralHadronEt = pfIso->sumNeutralHadronEt;
          electron.sumPhotonEt        = pfIso->sumPhotonEt;
          electron.sumPUPt            = pfIso->sumPUPt;
          cout << "e CH,NH,PH,PU = " << electron.sumChargedHadronPt << " " << electron.sumNeutralHadronEt << " " << electron.sumPhotonEt << " " << electron.sumPUPt << endl;

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
            electron.energyScaleEtUp       = el.userFloat("energyScaleEtUp"); //FIXME ???
            electron.energyScaleEtDown     = el.userFloat("energyScaleEtDown"); //FIXME ???
            electron.energySigmaUp         = el.userFloat("energySigmaUp");
            electron.energySigmaDown       = el.userFloat("energySigmaDown");
            electron.energySigmaPhiUp      = el.userFloat("energySigmaPhiUp");
            electron.energySigmaPhiDown    = el.userFloat("energySigmaPhiDown");
            electron.energySigmaRhoUp      = el.userFloat("energySigmaRhoUp");
            electron.energySigmaRhoDown    = el.userFloat("energySigmaRhoDown");
            cout << "some sys ... " << electron.ecalTrkEnergyPreCorr << " " << electron.ecalTrkEnergyPostCorr << " " << electron.energyScaleValue << " " << electron.energyScaleUp << " " << electron.energySigmaPhiUp << endl;
          }
          
          out_electrons.emplace_back( electron );
        }
        
        // photons ===================================================
          edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( 0 );   
          edm::Ptr<reco::Vertex> diphoton_vertex = dipho->vtx();
          
          const flashgg::Photon* leadPho    = dipho->leadingPhoton();
          const flashgg::Photon* subleadPho = dipho->subLeadingPhoton();
          bool passMVAs = checkPassMVAs(leadPho, subleadPho, diphoton_vertex);
          
          vector<const flashgg::Photon*> phs = { leadPho, subleadPho };
          for(const flashgg::Photon* ph : phs){
            photon.pt  = ph->pt();
            photon.eta = ph->eta();
            photon.phi = ph->phi();
            photon.mva_value = passMVAs;
            out_photons.emplace_back( photon );
          }
        
        // weights ===================================================
        
        outTree->Fill();
        cout << "outTree->Fill();" << endl;
        return;
      // ========================================================================================================
    } 

  } 

  typedef flashgg::GrinderHHWWggTagProducer GrinderFlashggHHWWggTagProducer;
  DEFINE_FWK_MODULE( GrinderFlashggHHWWggTagProducer );
  
  
  
  
  
  
