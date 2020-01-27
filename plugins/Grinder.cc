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

#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include <DataFormats/VertexReco/interface/VertexFwd.h>   // reco::VertexCollection
#include <DataFormats/VertexReco/interface/Vertex.h>      // reco::Vertex

#include <FWCore/Utilities/interface/EDMException.h>      // edm::Exception

#include <FWCore/ServiceRegistry/interface/Service.h>     // TFileService
#include <CommonTools/UtilAlgos/interface/TFileService.h> // TFileService

#include <TTree.h>

// Grinder
#include "Analysis/GRINDER/interface/Event.hh"
using namespace grinder;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class Grinder : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Grinder(const edm::ParameterSet&);
      ~Grinder();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


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

      edm::EDGetTokenT<reco::VertexCollection> primaryVerticesToken;

      TTree *outTree;
      grinder::Event  event;
      grinder::Event *event_ptr;

      grinder::Jet jet;
      grinder::Muon muon;
      grinder::Photon photon;
      grinder::Electron electron;

      std::vector<grinder::Jet>      jets;
      std::vector<grinder::Muon>     muons;
      std::vector<grinder::Photon>   photons;
      std::vector<grinder::Electron> electrons;

      std::vector<grinder::Jet>      *jets_ptr;
      std::vector<grinder::Muon>     *muons_ptr;
      std::vector<grinder::Photon>   *photons_ptr;
      std::vector<grinder::Electron> *electrons_ptr;

      // Photons
      std::string photon_loose_id_token, photon_medium_id_token, photon_tight_id_token, photon_mva_token, photon_mva_token_val, photon_mva_token_cat;
      EffectiveAreas effAreaChHadrons;
      EffectiveAreas effAreaNeuHadrons;
      EffectiveAreas effAreaPhotons;
      float superCluster_eta;
      edm::Handle< double > rho;
      edm::EDGetTokenT<double> rhoToken;

      // Electrons
      std::string electron_loose_id_token, electron_medium_id_token, electron_tight_id_token;
      const reco::GsfElectron::PflowIsolationVariables * pfIso;
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

  muonToken     = consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons_token"));
  photonToken   = consumes<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons_token"));
  jetToken      = consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets_token"));
  metToken      = consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets_token"));
  electronToken = consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons_token"));
  tauToken      = consumes<edm::View<pat::Tau>>(iConfig.getParameter<edm::InputTag>("taus_token"));
  rhoToken      = consumes<double>(iConfig.getParameter<edm::InputTag>("rho_token"));

  primaryVerticesToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertex_token"));
  //FIXME hltPrescalesToken = consumes<pat::PackedTriggerPrescales>(cfg.getParameter<edm::InputTag>("hltPrescales"));
  //FIXME l1tPrescalesToken = consumes<pat::PackedTriggerPrescales>(cfg.getParameter<edm::InputTag>("l1tPrescales"));

  outTree = fileService->make<TTree>("Events", "Events");

  // read Photons options
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

  // setup output data
  event_ptr     = &event;
  jets_ptr      = &jets;
  muons_ptr     = &muons;
  photons_ptr   = &photons;
  electrons_ptr = &electrons;

  outTree->Branch("Event", &event_ptr);
  outTree->Branch("Photons", &photons_ptr);
  outTree->Branch("Electrons", &electrons_ptr);
  outTree->Branch("Muons", &muons_ptr);
  outTree->Branch("Jets", &jets_ptr);
}


Grinder::~Grinder(){
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void Grinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  // Set Event Info
  event.run   = iEvent.id().run();
  event.lumi  = iEvent.luminosityBlock();
  event.event = iEvent.id().event();
  if(iEvent.isRealData()) event.bunchCrossing = iEvent.bunchCrossing();

  iEvent.getByToken(rhoToken, rho);
  event.angular_pt_density = (*rho);

  // Read prescales
/*
  edm::Handle<pat::PackedTriggerPrescales> hltPrescales;
  edm::Handle<pat::PackedTriggerPrescales> l1tPrescales;
  iEvent.getByToken(hltPrescalesToken, hltPrescales);
  iEvent.getByToken(l1tPrescalesToken, l1tPrescales);
  
  event.prescale = hltPrescales->getPrescaleForIndex(t.second.index) * l1tPrescales->getPrescaleForIndex(t.second.index);
*/

  // Read primary vertices collection // FIXME
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(primaryVerticesToken, vertices);
  if (vertices->size() == 0){
    edm::Exception excp(edm::errors::LogicError);
    excp << "Event contains zero good primary vertices.\n";
    excp.raise();
  }
  const reco::Vertex PrimaryVertex = vertices->front();

  // iterate over photons
  edm::Handle<edm::View<pat::Photon>> srcPhotons;
  iEvent.getByToken(photonToken, srcPhotons);
  for (unsigned i = 0; i < srcPhotons->size(); ++i){
    pat::Photon const &ph = srcPhotons->at(i);

    photon.pt  = ph.pt();
    photon.eta = ph.eta();
    photon.phi = ph.phi();

    // if( photon.pt <  )

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
    superCluster_eta = fabs( ph.superCluster()->eta());
    photon.sumChargedHadronPt = std::max( 0.0, ph.userFloat("phoChargedIsolation")       - (*rho) * effAreaChHadrons.getEffectiveArea( superCluster_eta )  );
    photon.sumNeutralHadronEt = std::max( 0.0, ph.userFloat("phoNeutralHadronIsolation") - (*rho) * effAreaNeuHadrons.getEffectiveArea( superCluster_eta ) );
    photon.sumPhotonEt        = std::max( 0.0, ph.userFloat("phoPhotonIsolation")        - (*rho) * effAreaPhotons.getEffectiveArea( superCluster_eta )    );

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

  // iterate over electrons
  edm::Handle<edm::View<pat::Electron>> srcElectrons;
  iEvent.getByToken(electronToken, srcElectrons);
  for (unsigned i = 0; i < srcElectrons->size(); ++i){
    pat::Electron const &el = srcElectrons->at(i);

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

  // iterate over muons
  edm::Handle<edm::View<pat::Muon>> srcMuons;
  iEvent.getByToken(muonToken, srcMuons);

  for (unsigned i = 0; i < srcMuons->size(); ++i){
    pat::Muon const &mu = srcMuons->at(i);

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

  // iterate over jets TODO
  edm::Handle<edm::View<pat::Jet>> srcJets;
  iEvent.getByToken(jetToken, srcJets);
  for (unsigned i = 0; i < srcJets->size(); ++i){
    pat::Jet const &j = srcJets->at(i);

    // RAW P4 vector 
    // https://twiki.cern.ch/twiki/bin/view/CMS/TopJME#Jets
    reco::Candidate::LorentzVector const &rawP4 = j.correctedP4("Uncorrected");
    jet.pt  = rawP4.pt();
    jet.eta = rawP4.eta();
    jet.phi = rawP4.phi();
    jet.m   = rawP4.mass();

    jet.charge = j.jetCharge();
    jet.area   = j.jetArea();

    jets.emplace_back( jet );
  }

  // iterate over MET TODO
/*
*/


  // Fill the output tree
  outTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void Grinder::beginJob(){
}

// ------------ method called once each job just after ending the event loop  ------------
void Grinder::endJob(){
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










