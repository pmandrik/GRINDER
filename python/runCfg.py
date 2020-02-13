
import FWCore.ParameterSet.Config as cms

### Options
DEBUG_print_content = False
IS_DATA = True
YEAR_ERA = "2016" # DEFAULT
###

process = cms.Process("Demo")
process.path = cms.Path()
process.TFileService = cms.Service('TFileService', fileName = cms.string('grinder_out.root'))

process.load("FWCore.MessageService.MessageLogger_cfi")

### Input
process.source           = cms.Source("PoolSource")
process.source.fileNames = cms.untracked.vstring('file:/eos/cms/store/data/Run2016H/DoubleMuon/MINIAOD/03Feb2017_ver3-v1/50000/8643C759-9BEB-E611-A7CA-008CFA111270.root')
process.source.fileNames = cms.untracked.vstring('file:/eos/cms/store/data/Run2016B/DoubleMuon/MINIAOD/17Jul2018_ver2-v1/00000/3CB1477F-F98A-E811-B1BC-0CC47A4D760C.root')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

### Trigger options
if YEAR_ERA == "2016":
  pass

### Photons options
# https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Photon_ID_Working_Points_WP_defi
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#ID_information
if YEAR_ERA == "2016":
  # ID
  from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
  setupEgammaPostRecoSeq(process, runEnergyCorrections=False, era='2016-Legacy')
  photon_loose_id  = "cutBasedPhotonID-Fall17-94X-V1-loose"
  photon_medium_id = "cutBasedPhotonID-Fall17-94X-V1-medium"
  photon_tight_id  = "cutBasedPhotonID-Fall17-94X-V1-tight"
  photon_mva  = "PhotonMVAEstimatorRunIIFall17v1"
  # ISO
  effAreaChHad  = "RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt"
  effAreaNeuHad = "RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt"
  effAreaPho    = "RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt"
elif YEAR_ERA == "2017":
  # ID
  photon_loose_id  = "cutBasedPhotonID-Fall17-94X-V2-loose"
  photon_medium_id = "cutBasedPhotonID-Fall17-94X-V2-medium"
  photon_tight_id  = "cutBasedPhotonID-Fall17-94X-V2-tight"
  photon_mva  = "PhotonMVAEstimatorRunIIFall17v2"
  # ISO
  effAreaChHad  = "RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt"
  effAreaNeuHad = "RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt"
  effAreaPho    = "RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt"
# SF
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale
# Resolution / Correction

### Electrons options
if YEAR_ERA == "2016":
  electron_loose_id  = "cutBasedElectronID-Fall17-94X-V1-loose"
  electron_medium_id = "cutBasedElectronID-Fall17-94X-V1-medium"
  electron_tight_id  = "cutBasedElectronID-Fall17-94X-V1-tight"

### Jets options
jet_JEC_Uncertainty_datafile = "" # unused
if YEAR_ERA == "2016":
  pass
if YEAR_ERA == "2017":
  pass

### MET options
# set-up filter https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')
baddetEcallist = cms.vuint32( [872439604,872422825,872420274,872423218,872423215,872416066,872435036,872439336, 872420273,872436907,872420147,872439731,872436657,872420397,872439732,872439339, 872439603,872422436,872439861,872437051,872437052,872420649,872421950,872437185, 872422564,872421566,872421695,872421955,872421567,872437184,872421951,872421694, 872437056,872437057,872437313,872438182,872438951,872439990,872439864,872439609, 872437181,872437182,872437053,872436794,872436667,872436536,872421541,872421413, 872421414,872421031,872423083,872421439] )

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
process.load("Configuration.Geometry.GeometryECALHCAL_cff")

process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist, 
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
    )

### Extra paths
if DEBUG_print_content:
  process.content = cms.EDAnalyzer("EventContentAnalyzer")
  process.path += process.content

### Main path
process.grinderMain = cms.EDAnalyzer('Grinder',
  # events info
  era_label            = cms.string( YEAR_ERA ),
  primaryVertex_token  = cms.InputTag('offlineSlimmedPrimaryVertices'),
  puSummaryToken_token = cms.InputTag('addPileupInfo'),
  rho_token            = cms.InputTag("fixedGridRhoFastjetAll"),
  rho_central_token    = cms.InputTag("fixedGridRhoFastjetCentral"),
  # trigger
  triggerResults_token   = cms.InputTag("TriggerResults","","HLT"),
  triggerPrescales_L1T_token = cms.InputTag("patTrigger", "l1min"),
  triggerPrescales_HLT_token = cms.InputTag("patTrigger"),
  # photons
  photons_token          = cms.InputTag('slimmedPhotons'),
  photon_loose_id_token  = cms.string(photon_loose_id),
  photon_medium_id_token = cms.string(photon_medium_id),
  photon_tight_id_token  = cms.string(photon_tight_id),
  photon_mva_token       = cms.string(photon_mva),
  effAreaChHadFile       = cms.FileInPath( effAreaChHad ),
  effAreaNeuHadFile      = cms.FileInPath( effAreaNeuHad ),
  effAreaPhoFile         = cms.FileInPath( effAreaPho ),
  # electrons
  electrons_token          = cms.InputTag('slimmedElectrons'),
  electron_loose_id_token  = cms.string(electron_loose_id),
  electron_medium_id_token = cms.string(electron_medium_id),
  electron_tight_id_token  = cms.string(electron_tight_id),
  # muons
  muons_token    = cms.InputTag('slimmedMuons'),
  # jets
  jets_token     = cms.InputTag('slimmedJets'),
  jet_type_label = cms.string('AK4PF'),
  jet_JEC_Uncertainty_datafile_token = cms.InputTag( jet_JEC_Uncertainty_datafile ),
  # mets
  mets_token             = cms.InputTag('slimmedMETs'),
  metFilterResults_token = cms.InputTag("TriggerResults","","RECO") if IS_DATA else cms.InputTag("TriggerResults","","PAT"),
  # taus
  taus_token    = cms.InputTag('slimmedTaus'),
  # genjets
  genjets_token = cms.InputTag('slimmedGenJets'),
  # cuts are a little bit weaker than the cuts from the articles
  # https://arxiv.org/pdf/1804.02716.pdf
  cut_photon_pt  = cms.double(17.5),
  cut_photon_eta = cms.double(2.6),
  # https://arxiv.org/pdf/1812.10529.pdf
  cut_electron_pt   = cms.double(17.5),
  cut_electron_eta  = cms.double(2.6),
  # https://arxiv.org/pdf/1807.06325.pdf
  # https://arxiv.org/pdf/1706.09936.pdf
  cut_muon_pt   = cms.double(17.5),
  cut_muon_eta  = cms.double(2.5),
  #
  cut_jet_pt    = cms.double(20),
  cut_jet_eta   = cms.double(2.5),
)

process.content = cms.EDAnalyzer("EventContentAnalyzer")
# process.path += process.content
process.path += process.ecalBadCalibReducedMINIAODFilter * process.grinderMain














