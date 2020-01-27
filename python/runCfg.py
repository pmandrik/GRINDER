
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

### Extra paths
if DEBUG_print_content:
  process.content = cms.EDAnalyzer("EventContentAnalyzer")
  process.path += process.content

### Main path
process.grinderMain = cms.EDAnalyzer('Grinder',
  # events info
  primaryVertex_token = cms.InputTag('offlineSlimmedPrimaryVertices'),
  # photons
  photons_token = cms.InputTag('slimmedPhotons'),
  photon_loose_id_token = cms.string(photon_loose_id),
  photon_medium_id_token = cms.string(photon_medium_id),
  photon_tight_id_token = cms.string(photon_tight_id),
  photon_mva_token = cms.string(photon_mva),
  effAreaChHadFile = cms.FileInPath( effAreaChHad ),
  effAreaNeuHadFile= cms.FileInPath( effAreaNeuHad ),
  effAreaPhoFile   = cms.FileInPath( effAreaPho ),
  rho_token        = cms.InputTag("fixedGridRhoFastjetAll"),
  # electrons
  electrons_token = cms.InputTag('slimmedElectrons'),
  electron_loose_id_token = cms.string(electron_loose_id),
  electron_medium_id_token = cms.string(electron_medium_id),
  electron_tight_id_token = cms.string(electron_tight_id),
  # muons
  muons_token = cms.InputTag('slimmedMuons'),
  # jets
  jets_token = cms.InputTag('slimmedJets'),
  # mets
  mets_token = cms.InputTag('slimmedMETs'),
  # taus
  taus_token = cms.InputTag('slimmedTaus'),
)

process.content = cms.EDAnalyzer("EventContentAnalyzer")
# process.path += process.content
process.path += process.grinderMain









