
import FWCore.ParameterSet.Config as cms

### Options
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.register('isData',    True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,   '')
options.register('yearEra', '2016', VarParsing.multiplicity.singleton, VarParsing.varType.string, '')
options.parseArguments()

DEBUG_print_content = False
IS_DATA  = options.isData
YEAR_ERA = options.yearEra
###

print options.yearEra, options.isData

process = cms.Process("Demo")
process.path = cms.Path()
process.TFileService = cms.Service('TFileService', fileName = cms.string('grinder_out.root'))

process.load("FWCore.MessageService.MessageLogger_cfi")

### Input
# pathes from https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
# 2016 - DoubleEG
# 2017 - DoubleEG
# 2018 - EGamma
process.source           = cms.Source("PoolSource")
if IS_DATA : 
  if YEAR_ERA == "2016": process.source.fileNames = cms.untracked.vstring('file:/eos/cms/store/data/Run2016B/DoubleEG/MINIAOD/17Jul2018_ver2-v1/00000/06360CEF-9B8D-E811-8A07-008CFA1111B4.root')
  if YEAR_ERA == "2017": process.source.fileNames = cms.untracked.vstring('file:/eos/cms/store/data/Run2017C/DoubleEG/MINIAOD/31Mar2018-v1/00000/72D2823F-3B38-E811-9260-7CD30ABD2EE8.root')
  if YEAR_ERA == "2018": process.source.fileNames = cms.untracked.vstring('file:/eos/cms/store/data/Run2018C/EGamma/MINIAOD/17Sep2018-v1/00000/91003A57-CDF4-134E-8CAB-D2B86D95B1E5.root')
else:
  if YEAR_ERA == "2016": process.source.fileNames = cms.untracked.vstring('file:/eos/cms/store/mc/RunIISummer16MiniAODv3/...')
  if YEAR_ERA == "2017": process.source.fileNames = cms.untracked.vstring('file:/eos/cms/store/mc/RunIIFall17MiniAODv2/...')
  if YEAR_ERA == "2018": process.source.fileNames = cms.untracked.vstring('file:/eos/cms/store/mc/RunIIAutumn18MiniAOD/...')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

### Global Tag Options
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JecGlobalTag
GT_name = ""
if IS_DATA : 
  if YEAR_ERA == "2016": GT_name = "94X_dataRun2_v10"  # MiniAOD
  if YEAR_ERA == "2017": GT_name = "94X_dataRun2_v11"  # MiniAOD
  if YEAR_ERA == "2018": GT_name = "102X_dataRun2_v12" # 2018ABC
else:
  if YEAR_ERA == "2016": GT_name = "94X_mcRun2_asymptotic_v3"       # MiniAOD
  if YEAR_ERA == "2017": GT_name = "94X_mc2017_realistic_v17"       # MiniAOD
  if YEAR_ERA == "2018": GT_name = "102X_upgrade2018_realistic_v20" # MiniAOD

from Configuration.AlCa.GlobalTag import GlobalTag
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag = GlobalTag(process.GlobalTag, GT_name, '')
process.GlobalTag.globaltag = GT_name

### Trigger options
# https://twiki.cern.ch/twiki/bin/view/CMS/TriggerStudies
selections_triggers = []
if YEAR_ERA == "2016":
  selections_triggers = [ 'HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55',
                          'HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55',
                          'HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70',
                          'HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90',
                          'HLT_DoublePhoton60',
                          'HLT_DoublePhoton85']
if YEAR_ERA == "2017":
  selections_triggers = [ 'HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55',
                          'HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55',
                          'HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55',
                          'HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90',
                          'HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95',
                          'HLT_DoublePhoton70',
                          'HLT_DoublePhoton85']
if YEAR_ERA == "2018":
  selections_triggers = [ 'HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55',
                          'HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55',
                          'HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto',
                          'HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90',
                          'HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95',
                          'HLT_DoublePhoton70',
                          'HLT_DoublePhoton85',
                          'HLT_Photon100EE_TightID_TightIso']

### Setup Electron / Photon sequence
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#2018_MiniAOD
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
if YEAR_ERA == "2016":
  setupEgammaPostRecoSeq(process, runEnergyCorrections=False, era='2016-Legacy')
if YEAR_ERA == "2017":
  setupEgammaPostRecoSeq(process,runVID=False, era='2017-Nov17ReReco') #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
if YEAR_ERA == "2018":
  setupEgammaPostRecoSeq(process,era='2018-Prompt')  

### Photons options
# https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Photon_ID_Working_Points_WP_defi
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#ID_information
# https://twiki.cern.ch/twiki/bin/view/CMS/MultivariatePhotonIdentificationRun2#Recommended_MVA_Recipe_for_regul <- Fallv2
# https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
if YEAR_ERA == "2016":
  # ID
  photon_loose_id  = "cutBasedPhotonID-Spring16-V2p2-loose"
  photon_medium_id = "cutBasedPhotonID-Spring16-V2p2-medium"
  photon_tight_id  = "cutBasedPhotonID-Spring16-V2p2-tight"
  photon_mva  = "PhotonMVAEstimatorRunIIFall17v1"
  # ISO
  effAreaChHad  = "RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfChargedHadrons_90percentBased.txt"
  effAreaNeuHad = "RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased.txt"
  effAreaPho    = "RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased.txt"
if YEAR_ERA == "2017":
  # ID
  photon_loose_id  = "cutBasedPhotonID-Fall17-94X-V2-loose"
  photon_medium_id = "cutBasedPhotonID-Fall17-94X-V2-medium"
  photon_tight_id  = "cutBasedPhotonID-Fall17-94X-V2-tight"
  photon_mva  = "PhotonMVAEstimatorRunIIFall17v2"
  # ISO
  effAreaChHad  = "RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt"
  effAreaNeuHad = "RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt"
  effAreaPho    = "RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt"
if YEAR_ERA == "2018":
  # ID
  photon_loose_id  = "cutBasedPhotonID-Fall17-94X-V2-loose"
  photon_medium_id = "cutBasedPhotonID-Fall17-94X-V2-medium"
  photon_tight_id  = "cutBasedPhotonID-Fall17-94X-V2-tight"
  photon_mva  = "PhotonMVAEstimatorRunIIFall17v2"
  # ISO
  effAreaChHad  = "RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt"
  effAreaNeuHad = "RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt"
  effAreaPho    = "RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt"
  
  
### Electrons options
# https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_formats
if YEAR_ERA == "2016":
  electron_loose_id  = "cutBasedElectronID-Fall17-94X-V1-loose"
  electron_medium_id = "cutBasedElectronID-Fall17-94X-V1-medium"
  electron_tight_id  = "cutBasedElectronID-Fall17-94X-V1-tight"
if YEAR_ERA == "2017":
  electron_loose_id  = "cutBasedElectronID-Fall17-94X-V2-loose"
  electron_medium_id = "cutBasedElectronID-Fall17-94X-V2-medium"
  electron_tight_id  = "cutBasedElectronID-Fall17-94X-V2-tight"
if YEAR_ERA == "2018":
  electron_loose_id  = "cutBasedElectronID-Fall17-94X-V2-loose"
  electron_medium_id = "cutBasedElectronID-Fall17-94X-V2-medium"
  electron_tight_id  = "cutBasedElectronID-Fall17-94X-V2-tight"

# SF
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale
# Resolution / Correction

### Jets options
# re-apply JEC
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections?rev=139#CorrPatJets
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
  process,
  jetSource = cms.InputTag('slimmedJets'),
  labelName = 'UpdatedJEC', # -> updatedPatJetsUpdatedJEC
  jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')  # Do not forget 'L2L3Residual' on data!
)


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
  is_data              = cms.bool( IS_DATA ),
  era_label            = cms.string( YEAR_ERA ),
  primaryVertex_token  = cms.InputTag('offlineSlimmedPrimaryVertices'),
  puSummaryToken_token = cms.InputTag('addPileupInfo'),
  rho_token            = cms.InputTag("fixedGridRhoFastjetAll"),
  rho_central_token    = cms.InputTag("fixedGridRhoFastjetCentral"),
  # trigger
  triggerResults_token   = cms.InputTag("TriggerResults","","HLT"),
  triggerPrescales_L1T_token = cms.InputTag("patTrigger", "l1min"),
  triggerPrescales_HLT_token = cms.InputTag("patTrigger"),
  selections_triggers_names  = cms.vstring( selections_triggers ),
  do_trigger_filtering       = cms.bool(True),
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
  jets_token     = cms.InputTag('updatedPatJetsUpdatedJEC'), # cms.InputTag('slimmedJets'),
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
process.path += process.ecalBadCalibReducedMINIAODFilter * process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.grinderMain














