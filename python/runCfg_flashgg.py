#!/usr/bin/env cmsRun
print "======================================================> 0"

from importlib import import_module
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing
from flashgg.Systematics.SystematicDumperDefaultVariables import minimalVariables,minimalHistograms,minimalNonSignalVariables,systematicVariables
from flashgg.Systematics.SystematicDumperDefaultVariables import minimalVariablesHTXS,systematicVariablesHTXS
import os
from flashgg.MetaData.MetaConditionsReader import *
# from flashgg.Systematics.flashggDiPhotonSystematics_cfi import flashggDiPhotonSystematics

print "======================================================> 1"

# SYSTEMATICS SECTION
dropVBFInNonGold = False  # for 2015 only!

process = cms.Process("FLASHggSyst")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
# process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1 )

systlabels = [""]
phosystlabels = []
metsystlabels = []
jetsystlabels = []
elesystlabels = []
musystlabels  = []

from flashgg.MetaData.JobConfig import customize
customize.options.register('HHWWggTagsOnly',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'HHWWggTagsOnly'
                           )
customize.options.register('doHHWWggTag',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doHHWWggTag'
                           )
customize.options.register('doHHWWggTagCutFlow', # This saves all events for cutflow analysis
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doHHWWggTagCutFlow'
                           ),                       
customize.options.register('doBJetRegression',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doBJetRegression'
                           )
customize.options.register('doHTXS',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doHTXS'
                           )
customize.options.register('doMuFilter',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doMuFilter'
                           )
customize.options.register('doSystematics',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doSystematics'
                           )
customize.options.register('ignoreNegR9',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'ignoreNegR9'
                           )

print "======================================================> 2"

print "Printing defaults"
# import flashgg customization to check if we have signal or background
from flashgg.MetaData.JobConfig import customize
# set default options if needed
customize.setDefault("targetLumi", 1000. )

print "customize.parse() ... "
customize.parse()

print "MetaConditionsReader(customize.metaConditions) ... "
customize.metaConditions = MetaConditionsReader(customize.metaConditions)

print "======= Global Tag"
print customize.metaConditions
print customize.metaConditions['globalTags']
### Global Tag
from Configuration.AlCa.GlobalTag import GlobalTag
if customize.processId == "Data":
    process.GlobalTag.globaltag = str(customize.metaConditions['globalTags']['data'])
else:
    process.GlobalTag.globaltag = str(customize.metaConditions['globalTags']['MC'])


from flashgg.Systematics.SystematicsCustomize import *
if True : # jetSystematicsInputTags = createStandardSystematicsProducers(process , customize)
    print "Keep only di-photon systematics ... "
    process.load("flashgg/Taggers/flashggTagSequence_cfi")
    process.load("flashgg.Systematics.flashggDiPhotonSystematics_cfi")

    from flashgg.Taggers.flashggTagSequence_cfi import *
    process.flashggTagSequence = flashggPrepareTagSequence(process, customize.metaConditions)
    print "PATH =================", process.flashggTagSequence
    
    import flashgg.Systematics.flashggDiPhotonSystematics_cfi as diPhotons_syst
    diPhotons_syst.setupDiPhotonSystematics( process, customize )
    print "PATH =================", process.flashggTagSequence

    applyL1Prefiring = customizeForL1Prefiring(process, customize.metaConditions, customize.processId)

if True : # modifyTagSequenceForSystematics(process,jetSystematicsInputTags) # normally uncommented 
    process.flashggTagSequence.remove(process.flashggUnpackedJets) # to avoid unnecessary cloning
    process.flashggTagSequence.remove(process.flashggDifferentialPhoIdInputsCorrection) # Needs to be run before systematics
    from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet,massSearchReplaceAnyInputTag
    # massSearchReplaceAnyInputTag(process.flashggTagSequence,cms.InputTag("flashggDifferentialPhoIdInputsCorrection"),cms.InputTag("flashggDiPhotonSystematics"))
    process.flashggDiPhotonSystematics.src = cms.InputTag("flashggPreselectedDiPhotons")
    # massSearchReplaceAnyInputTag(process.flashggTagSequence,cms.InputTag("flashggPreselectedDiPhotons"),cms.InputTag("flashggDiPhotonSystematics"))
    from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag
    
    process.systematicsTagSequences = cms.Sequence()

print "Initial di-photon systematic:"
for pset in process.flashggDiPhotonSystematics.SystMethods:
  print  pset.Label.value()

print "======================================================> 3"


# process.load("flashgg/Taggers/flashggTagSequence_cfi")
# process.flashggTagSequence = flashggPrepareTagSequence(customize.metaConditions)

# needed for 0th vertex from microAOD
# HHWWgg: Want zeroeth vertex 
if customize.HHWWggTagsOnly:
    process.load("flashgg/MicroAOD/flashggDiPhotons_cfi")
    process.flashggDiPhotons.whichVertex = cms.uint32(0)
    process.flashggDiPhotons.useZerothVertexFromMicro = cms.bool(True)
    if customize.HHWWggTagsOnly: # not sure if this is needed for tthTagsOnly, but it is needed for HHWWgg 
        process.flashggDiPhotons.vertexProbMVAweightfile = "flashgg/MicroAOD/data/TMVAClassification_BDTVtxProb_SL_2016.xml" # Prob or Id ? 
        process.flashggDiPhotons.vertexIdMVAweightfile = "flashgg/MicroAOD/data/TMVAClassification_BDTVtxId_SL_2016.xml"

# if customize.HHWWggTagsOnly:
#     print'customizing for HHWWgg'
#     process.load("flashgg/MicroAOD/flashggDiPhotons_cfi")
#     process.flashggDiPhotons.whichVertex = cms.uint32(0)
#     process.flashggDiPhotons.useZerothVertexFromMicro = cms.bool(True)
#     process.flashggDiPhotons.vertexProbMVAweightfile = "flashgg/MicroAOD/data/TMVAClassification_BDTVtxId_SL_2016.xml"
#     process.flashggDiPhotons.vertexIdMVAweightfile = "flashgg/MicroAOD/data/TMVAClassification_BDTVtxId_SL_2016.xml"
                                                        
print 'here we print the tag sequence before'
print process.flashggTagSequence

#if customize.doHHWWggTag :
    #import flashgg.Systematics.HHWWggCustomize 
    #hhwwggc = flashgg.Systematics.HHWWggCustomize.HHWWggCustomize(process, customize, customize.metaConditions)
    #minimalVariables += hhwwggc.variablesToDump()
    #systematicVariables = hhwwggc.systematicVariables()
    
#if customize.doHHWWggTag:
    #import Analysis.GRINDER.HHWWggCustomize 
    #hhwwggc = Analysis.GRINDER.HHWWggCustomize.HHWWggCustomize(process, customize, customize.metaConditions)
    #minimalVariables += hhwwggc.variablesToDump()
    #systematicVariables = hhwwggc.systematicVariables()
    
if customize.doHHWWggTag:
    import Analysis.GRINDER.HHWWggCustomize 
    hhwwggc = Analysis.GRINDER.HHWWggCustomize.HHWWggCustomize(process, customize, customize.metaConditions)
    
    print "customizeTagSequencecustomizeTagSequencecustomizeTagSequencecustomizeTagSequencecustomizeTagSequence\n\n\n\n!!!!!!!!!"
    DEBUG_print_content = False
    IS_DATA  = False
    YEAR_ERA = "2017" # FIXME
    
    ### Electrons options ########################################################################################################################
    # https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_formats
    if YEAR_ERA == "2016":
      electron_loose_id  = "cutBasedElectronID-Fall17-94X-V1-loose"
      electron_medium_id = "cutBasedElectronID-Fall17-94X-V1-medium"
      electron_tight_id  = "cutBasedElectronID-Fall17-94X-V1-tight"
    if YEAR_ERA == "2017":
      electron_loose_id  = "cutBasedElectronID-Fall17-94X-V1-loose"  # FIXME cutBasedElectronID-Fall17-94X-V2-loose not available !!!
      electron_medium_id = "cutBasedElectronID-Fall17-94X-V1-medium"
      electron_tight_id  = "cutBasedElectronID-Fall17-94X-V1-tight"
    if YEAR_ERA == "2018":
      electron_loose_id  = "cutBasedElectronID-Fall17-94X-V2-loose"
      electron_medium_id = "cutBasedElectronID-Fall17-94X-V2-medium"
      electron_tight_id  = "cutBasedElectronID-Fall17-94X-V2-tight"

    # SF
    # https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale
    # Resolution / Correction

    
    import FWCore.ParameterSet.Config as cms
    from flashgg.Taggers.globalVariables_cff import globalVariables
    from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag # should include jet systematics 
    from flashgg.MicroAOD.flashggJets_cfi import  maxJetCollections, flashggDeepCSV

    print "UnpackedJetCollectionVInputTag = ", UnpackedJetCollectionVInputTag
    # cms.InputTag("flashggUnpackedJets","0"), cms.InputTag("flashggUnpackedJets","1"), cms.InputTag("flashggUnpackedJets","2"), cms.InputTag("flashggUnpackedJets","3"), cms.InputTag("flashggUnpackedJets","4"), cms.InputTag("flashggUnpackedJets","5"), cms.InputTag("flashggUnpackedJets","6"), cms.InputTag("flashggUnpackedJets","7"), cms.InputTag("flashggUnpackedJets","8"), cms.InputTag("flashggUnpackedJets","9"), cms.InputTag("flashggUnpackedJets","10"), cms.InputTag("flashggUnpackedJets","11")

    """
with weights = 0.928702 Central
with weights = 0.997022 FracRVWeightCentral
with weights = 0.928084 FracRVWeightDown01sigma
with weights = 0.929319 FracRVWeightUp01sigma
with weights = 1.0007   LooseMvaSFCentral
with weights = 0.928702 LooseMvaSFDown01sigma
with weights = 0.928702 LooseMvaSFUp01sigma
with weights = 0.990978 PreselSFCentral
with weights = 0.914777 PreselSFDown01sigma
with weights = 0.942728 PreselSFUp01sigma
with weights = 0.99133  TriggerWeightCentral
with weights = 0.926837 TriggerWeightDown01sigma
with weights = 0.930568 TriggerWeightUp01sigma
with weights = 0.947514 electronVetoSFCentral
with weights = 0.919505 electronVetoSFDown01sigma
with weights = 0.937929 electronVetoSFUp01sigma
0.997022 * 1.0007 * 0.990978 * 0.99133 * 0.947514 = 0.9287023385000381
    """

    diphoton_systematics_tags = cms.vstring("FNUFEBDown01sigma", "FNUFEBUp01sigma", "FNUFEEDown01sigma", "FNUFEEUp01sigma", "MCScaleGain1EBDown01sigma", "MCScaleGain1EBUp01sigma", "MCScaleGain6EBDown01sigma", "MCScaleGain6EBUp01sigma", "MCScaleHighR9EBDown01sigma", "MCScaleHighR9EBUp01sigma", "MCScaleHighR9EEDown01sigma", "MCScaleHighR9EEUp01sigma", "MCScaleLowR9EBDown01sigma", "MCScaleLowR9EBUp01sigma", "MCScaleLowR9EEDown01sigma", "MCScaleLowR9EEUp01sigma", "MCSmearHighR9EBPhiDown01sigma", "MCSmearHighR9EBPhiUp01sigma", "MCSmearHighR9EBRhoDown01sigma", "MCSmearHighR9EBRhoUp01sigma", "MCSmearHighR9EEPhiDown01sigma", "MCSmearHighR9EEPhiUp01sigma", "MCSmearHighR9EERhoDown01sigma", "MCSmearHighR9EERhoUp01sigma", "MCSmearLowR9EBPhiDown01sigma", "MCSmearLowR9EBPhiUp01sigma", "MCSmearLowR9EBRhoDown01sigma", "MCSmearLowR9EBRhoUp01sigma", "MCSmearLowR9EEPhiDown01sigma", "MCSmearLowR9EEPhiUp01sigma", "MCSmearLowR9EERhoDown01sigma", "MCSmearLowR9EERhoUp01sigma", "MaterialCentralBarrelDown01sigma", "MaterialCentralBarrelUp01sigma", "MaterialForwardDown01sigma", "MaterialForwardUp01sigma", "MaterialOuterBarrelDown01sigma", "MaterialOuterBarrelUp01sigma", "MvaShiftDown01sigma", "MvaShiftUp01sigma", "ShowerShapeHighR9EBDown01sigma", "ShowerShapeHighR9EBUp01sigma", "ShowerShapeHighR9EEDown01sigma", "ShowerShapeHighR9EEUp01sigma", "ShowerShapeLowR9EBDown01sigma", "ShowerShapeLowR9EBUp01sigma", "ShowerShapeLowR9EEDown01sigma", "ShowerShapeLowR9EEUp01sigma", "SigmaEOverEShiftDown01sigma", "SigmaEOverEShiftUp01sigma")
    
    IS_DATA  = True   # options.isData
    YEAR_ERA = "2017" # options.yearEra
    process.GrinderflashggHHWWggTag = cms.EDAnalyzer("GrinderFlashggHHWWggTagProducer",
                                    ### NEW VARIABLES =====================================================>
                                    is_data              = cms.bool( IS_DATA ),
                                    era_label            = cms.string( YEAR_ERA ),
                                    electron_loose_id_token  = cms.string(electron_loose_id),
                                    electron_medium_id_token = cms.string(electron_medium_id),
                                    electron_tight_id_token  = cms.string(electron_tight_id),
                                    puSummaryToken_token = cms.InputTag('slimmedAddPileupInfo'), # cms.InputTag('addPileupInfo'),
                                    genjets_token   = cms.InputTag('slimmedGenJets'),
                                    generator_token = cms.InputTag('generator'),
                                    lumiWeight      = cms.double(-999.0),
                                    ### OLD VARIABLES =====================================================>
                                    globalVariables=globalVariables,
                                    PhotonTag = cms.InputTag('flashggRandomizedPhotons'),
                                    DiPhotonTag = cms.InputTag('flashggDiPhotonSystematics'),
                                    # DiPhotonTag = cms.InputTag('flashggPreselectedDiPhotons'),
                                    # DiPhotonTag = cms.InputTag('flashggDiPhotons'),
                                    SystLabel = cms.string(""),
                                    JetsName = cms.string("bRegProducer"), # WHAT ???
                                    JetsCollSize = cms.uint32(maxJetCollections), #
                                    JetsSuffixes = cms.vstring(''), #nominal and systematic variations 
                                    DiPhotonName = cms.string('flashggDiPhotonSystematics'), # 
                                    VertexTag = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                    GenParticleTag         = cms.InputTag('flashggPrunedGenParticles'),
                                    ElectronTag            = cms.InputTag('flashggSelectedElectrons'),
                                    MuonTag                = cms.InputTag('flashggSelectedMuons'),
                                    METTag                 = cms.InputTag('flashggMets'),
                                    # METTag                 = cms.InputTag('flashggMetsCorr'), # RunIIFall17-3-2-0 contains these and NOT flashggMets
                                    JetTags                = UnpackedJetCollectionVInputTag, # 
                                    DiPhotonSuffixes = diphoton_systematics_tags, #nominal and systematic variations
                                    useVertex0only=cms.bool(False),
                                    MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                    leptonPtThreshold = cms.double(10.),
                                    muonEtaThreshold = cms.double(2.4),
                                    MVAThreshold = cms.double(0.0),                                                     
                                    deltaRMuonPhoThreshold = cms.double(0.4),
                                    jetsNumberThreshold = cms.double(99.), # originially 3. 
                                    jetPtThreshold = cms.double(25.),
                                    jetEtaThreshold= cms.double(2.4),
                                    deltaRPhoLeadJet = cms.double(0.4),
                                    deltaRPhoSubLeadJet = cms.double(0.4),
                                    muPFIsoSumRelThreshold = cms.double(0.15),
                                    deltaRJetMuonThreshold = cms.double(0.4),
                                    PuIDCutoffThreshold = cms.double(0.8),
                                    PhoMVAThreshold = cms.double(-0.9),
                                    METThreshold = cms.double(0.),
                                    DeltaRTrkElec = cms.double(.4),
                                    TransverseImpactParam = cms.double(0.02),
                                    LongitudinalImpactParam = cms.double(0.2),
                                    deltaRPhoElectronThreshold = cms.double(0.4), # was 1 
                                    deltaMassElectronZThreshold = cms.double(10.),                                  
                                    electronEtaThresholds=cms.vdouble(1.4442,1.566,2.5),
                                    nonTrigMVAThresholds = cms.vdouble(0.913286,0.805013,0.358969),
                                    nonTrigMVAEtaCuts = cms.vdouble(0.8,1.479,2.5),
                                    electronIsoThreshold = cms.double(0.15),
                                    electronNumOfHitsThreshold = cms.double(1),
                                    useElectronMVARecipe = cms.bool(False),
                                    useElectronLooseID = cms.bool(True),                                    
                                    rhoTag = cms.InputTag('fixedGridRhoFastjetAll'),
                                    RECOfilters = cms.InputTag('TriggerResults::RECO'),
                                    PATfilters = cms.InputTag('TriggerResults::PAT'),
                                    FLASHfilters = cms.InputTag('TriggerResults::FLASHggMicroAOD'),
                                    # bTag = cms.string(flashggDeepCSV),
                                    doHHWWggTagCutFlowAnalysis = cms.bool(False) # save events for cut flow analysis                                       
                                    )

    process.GrinderflashggHHWWggTag.doHHWWggTagCutFlowAnalysis = cms.bool(True)
    print'Removing single Higgs tags'
    process.flashggTagSequence.remove(process.flashggVBFTag)
    process.flashggTagSequence.remove(process.flashggTTHLeptonicTag)
    process.flashggTagSequence.remove(process.flashggTTHHadronicTag)
    process.flashggTagSequence.remove(process.flashggVHEtTag)
    process.flashggTagSequence.remove(process.flashggVHLooseTag)
    process.flashggTagSequence.remove(process.flashggVHTightTag)
    process.flashggTagSequence.remove(process.flashggVHMetTag)
    process.flashggTagSequence.remove(process.flashggWHLeptonicTag)
    process.flashggTagSequence.remove(process.flashggZHLeptonicTag)
    process.flashggTagSequence.remove(process.flashggVHLeptonicLooseTag)
    process.flashggTagSequence.remove(process.flashggVHHadronicTag)
    process.flashggTagSequence.remove(process.flashggVBFMVA)
    process.flashggTagSequence.remove(process.flashggVBFDiPhoDiJetMVA)
    process.flashggTagSequence.remove(process.flashggTTHDiLeptonTag)
    process.flashggTagSequence.remove(process.flashggUntagged)
    process.flashggTagSequence.remove(process.flashggUntagged)
    process.flashggTagSequence.remove(process.flashggTHQLeptonicTag)
    process.flashggTagSequence.remove(process.flashggPreselectedDiPhotons)
    process.flashggTagSequence.remove(process.flashggDiPhotonMVA)

    process.flashggTagSequence.remove(process.flashggPrefireDiPhotons)
    process.flashggTagSequence.remove(process.flashggVHhadMVA)

    process.flashggTagSequence.replace(process.flashggTagSorter, process.grinder_flashggHHWWggTagSequence)
    
    minimalVariables += []
    systematicVariables = []  
    
print process.flashggTagSequence
print "======================================================> 3"

#flashggPrefireDiPhotons+flashggPreselectedDiPhotons+flashggDiPhotonMVA+flashggVBFMVA+flashggVHhadMVA+flashggVBFDiPhoDiJetMVA+flashggUntagged+flashggVBFTag+flashggTTHDiLeptonTag+flashggTTHLeptonicTag+flashggTHQLeptonicTag+flashggTTHHadronicTag+flashggVHMetTag+flashggWHLeptonicTag+flashggZHLeptonicTag+flashggVHHadronicTag+flashggTagSorter


print 'here we print the tag sequence after'
print process.flashggTagSequence
print "customize.processId:",customize.processId

# Or use the official tool instead
useEGMTools(process)

# Only run systematics for signal events
# convention: ggh vbf wzh (wh zh) tth
signal_processes = ["ggh_","vbf_","wzh_","wh_","zh_","bbh_","thq_","thw_","tth_","HHTo2B2G","GluGluHToGG","VBFHToGG","VHToGG","ttHToGG","Acceptance","WWgg","HToAA"]
# ^^ WWgg present in HHWWgg signal samples 

# print'customize'
# print'checking customize options'
# print'customize.processId.count("ggF_X250_WWgg_qqlnugg") = ',customize.processId.count("ggF_X250_WWgg_qqlnugg")
# for thing in customize.processId.count(0):
    # print'thing = ',thing 
is_signal = reduce(lambda y,z: y or z, map(lambda x: customize.processId.count(x), signal_processes))
#if customize.processId.count("h_") or customize.processId.count("vbf_") or customize.processId.count("Acceptance") or customize.processId.count("hh_"): 

print "======================================================> 4"
if is_signal or True:
    print "Signal MC, so adding systematics and dZ"
    if customize.doHTXS:
        variablesToUse = minimalVariablesHTXS
    else:
        variablesToUse = minimalVariables

    if customize.doSystematics:
        for direction in ["Up","Down"]:
        # for direction in ["Up"]:
            phosystlabels.append("MvaShift%s01sigma" % direction)
#            phosystlabels.append("MvaLinearSyst%s01sigma" % direction)


################---- turning off to test MvaShift 
            phosystlabels.append("SigmaEOverEShift%s01sigma" % direction)
            phosystlabels.append("MaterialCentralBarrel%s01sigma" % direction)
            phosystlabels.append("MaterialOuterBarrel%s01sigma" % direction)
            phosystlabels.append("MaterialForward%s01sigma" % direction)
            phosystlabels.append("FNUFEB%s01sigma" % direction)
            phosystlabels.append("FNUFEE%s01sigma" % direction)
            phosystlabels.append("MCScaleGain6EB%s01sigma" % direction)
            phosystlabels.append("MCScaleGain1EB%s01sigma" % direction)
            jetsystlabels.append("JEC%s01sigma" % direction)
            jetsystlabels.append("JER%s01sigma" % direction)
            jetsystlabels.append("PUJIDShift%s01sigma" % direction)
            metsystlabels.append("metJecUncertainty%s01sigma" % direction)
            metsystlabels.append("metJerUncertainty%s01sigma" % direction)
            metsystlabels.append("metPhoUncertainty%s01sigma" % direction)
            metsystlabels.append("metUncUncertainty%s01sigma" % direction)
            variablesToUse.append("UnmatchedPUWeight%s01sigma[1,-999999.,999999.] := weight(\"UnmatchedPUWeight%s01sigma\")" % (direction,direction))
            variablesToUse.append("MvaLinearSyst%s01sigma[1,-999999.,999999.] := weight(\"MvaLinearSyst%s01sigma\")" % (direction,direction))
            variablesToUse.append("LooseMvaSF%s01sigma[1,-999999.,999999.] := weight(\"LooseMvaSF%s01sigma\")" % (direction,direction))
            variablesToUse.append("PreselSF%s01sigma[1,-999999.,999999.] := weight(\"PreselSF%s01sigma\")" % (direction,direction))
            variablesToUse.append("electronVetoSF%s01sigma[1,-999999.,999999.] := weight(\"electronVetoSF%s01sigma\")" % (direction,direction))
            variablesToUse.append("TriggerWeight%s01sigma[1,-999999.,999999.] := weight(\"TriggerWeight%s01sigma\")" % (direction,direction))
            variablesToUse.append("FracRVWeight%s01sigma[1,-999999.,999999.] := weight(\"FracRVWeight%s01sigma\")" % (direction,direction))
            # variablesToUse.append("FracRVNvtxWeight%s01sigma[1,-999999.,999999.] := weight(\"FracRVNvtxWeight%s01sigma\")" % (direction,direction)) # removed because not working for HHWWgg for some reason
            variablesToUse.append("ElectronWeight%s01sigma[1,-999999.,999999.] := weight(\"ElectronWeight%s01sigma\")" % (direction,direction))
            
            # 
            variablesToUse.append("MuonIDWeight%s01sigma[1,-999999.,999999.] := weight(\"Muon%sIDWeight%s01sigma\")" % (direction,str(customize.metaConditions["MUON_ID"]),direction))
            variablesToUse.append("ElectronIDWeight%s01sigma[1,-999999.,999999.] := weight(\"ElectronIDWeight%s01sigma\")" % (direction,direction))
            variablesToUse.append("ElectronRecoWeight%s01sigma[1,-999999.,999999.] := weight(\"ElectronRecoWeight%s01sigma\")" % (direction,direction))
            variablesToUse.append("MuonIsoWeight%s01sigma[1,-999999.,999999.] := weight(\"Muon%sISOWeight%s01sigma\")" % (direction,str(customize.metaConditions['MUON_ISO']),direction))
            
            #     variablesToUse.append("MuonWeight%s01sigma[1,-999999.,999999.] := weight(\"MuonWeight%s01sigma\")" % (direction,direction))
            #     variablesToUse.append("MuonMiniIsoWeight%s01sigma[1,-999999.,999999.] := weight(\"MuonMiniIsoWeight%s01sigma\")" % (direction,direction))

            # if os.environ["CMSSW_VERSION"].count("CMSSW_8_0"):
            #     variablesToUse.append("MuonWeight%s01sigma[1,-999999.,999999.] := weight(\"MuonWeight%s01sigma\")" % (direction,direction))
            #     variablesToUse.append("MuonMiniIsoWeight%s01sigma[1,-999999.,999999.] := weight(\"MuonMiniIsoWeight%s01sigma\")" % (direction,direction))
            # elif os.environ["CMSSW_VERSION"].count("CMSSW_9_4"):
            #     variablesToUse.append("MuonIDWeight%s01sigma[1,-999999.,999999.] := weight(\"Muon%sIDWeight%s01sigma\")" % (direction,MUON_ID,direction))
            #     variablesToUse.append("MuonIsoWeight%s01sigma[1,-999999.,999999.] := weight(\"Muon%sISOWeight%s01sigma\")" % (direction,MUON_ISO,direction))
            
            
            # 
            variablesToUse.append("JetBTagCutWeight%s01sigma[1,-999999.,999999.] := weight(\"JetBTagCutWeight%s01sigma\")" % (direction,direction))
            variablesToUse.append("JetBTagReshapeWeight%s01sigma[1,-999999.,999999.] := weight(\"JetBTagReshapeWeight%s01sigma\")" % (direction,direction))
            


            for r9 in ["HighR9","LowR9"]:
                for region in ["EB","EE"]:
                    phosystlabels.append("ShowerShape%s%s%s01sigma"%(r9,region,direction))
#                    phosystlabels.append("MCSmear%s%s%s01sigma" % (r9,region,direction))
                    phosystlabels.append("MCScale%s%s%s01sigma" % (r9,region,direction))
                    for var in ["Rho","Phi"]:
                        phosystlabels.append("MCSmear%s%s%s%s01sigma" % (r9,region,var,direction))
       
       
#############################################################
        
        systlabels += phosystlabels
        systlabels += jetsystlabels
        systlabels += metsystlabels
    customizeSystematicsForSignal(process)
elif customize.processId == "Data":
    print "Data, so turn off all shifts and systematics, with some exceptions"
    variablesToUse = minimalNonSignalVariables
    customizeSystematicsForData(process)
else:
    print "Background MC, so store mgg and central only"
    variablesToUse = minimalNonSignalVariables
    customizeSystematicsForBackground(process)

if customize.HHWWggTagsOnly:
    variablesToUse = minimalVariables
    if customize.processId == "Data":
        variablesToUse = minimalNonSignalVariables

print "--- Systematics  with independent collections ---"
print systlabels
print "-------------------------------------------------"
print "--- Variables to be dumped, including systematic weights ---"
print variablesToUse
print "------------------------------------------------------------"

print "======================================================> 5"

#from flashgg.Taggers.globalVariables_cff import globalVariables
#globalVariables.extraFloats.rho = cms.InputTag("rhoFixedGridAll")

# cloneTagSequenceForEachSystematic(process,systlabels,phosystlabels,jetsystlabels,jetSystematicsInputTags)
# cloneTagSequenceForEachSystematic(process,systlabels,phosystlabels,metsystlabels,jetsystlabels,jetSystematicsInputTags) # used in workspacestd 

# Dump an object called NoTag for untagged events in order to track QCD weights
# Will be broken if it's done for non-central values, so turn this on only for the non-syst tag sorter
# process.flashggTagSorter.CreateNoTag = True # MUST be after tag sequence cloning
# process.flashggTagSorter.CreateNoTag = False # MUST be after tag sequence cloning

###### Dumper section

from FWCore.ParameterSet.VarParsing import VarParsing

print "======================================================> 5a"
from flashgg.MetaData.samples_utils import SamplesManager
print "======================================================> 5b"

process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/user/spigazzi/flashgg/Era2017_RR-31Mar2018_v2/legacyRun2FullV1/DoubleEG/Era2017_RR-31Mar2018_v2-legacyRun2FullV1-v0-Run2017E-31Mar2018-v1/190606_095510/0000/myMicroAODOutputFile_639.root"))
process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_higgs/cmshgg/atishelm/flashgg/HHWWgg_v2-4/94X_mc2017-RunIIFall18/GluGluToHHTo_WWgg_qqlnu_node2/HHWWgg_v2-4-94X_mc2017-RunIIFall18-v0-atishelm-100000events_wPU_MINIAOD-5f646ecd4e1c7a39ab0ed099ff55ceb9/200429_093533/0000/myMicroAODOutputFile_1.root"))
process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_higgs/resonant_HH/RunII/MicroAOD/HHWWggSignal/HHWWgg_SM2017/94X_mc2017-RunIIFall18/GluGluToHHTo_WWgg_qqlnu_nodeSM/HHWWgg_SM2017-94X_mc2017-RunIIFall18-v0-atishelm-HHWWgg_SM2017_100000events_wPU_MINIAOD-5f646ecd4e1c7a39ab0ed099ff55ceb9/200724_084804/0000/myMicroAODOutputFile_32.root"))
#process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch:1094//store/user/spigazzi/flashgg/Era2017_RR-31Mar2018_v2/legacyRun2FullV1/VBFHToGG_M125_13TeV_amcatnlo_pythia8/Era2017_RR-31Mar2018_v2-legacyRun2FullV1-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/190703_112047/0000/myMicroAODOutputFile_68.root"))

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root"))

#if customize.options.WeightName :
    #lheProduct = customize.dataset[1]["LHESourceName"].split("_")
    ##print lheProduct
    #process.tagsDumper.LHEEventProduct = cms.untracked.InputTag( str(lheProduct[1]) , str(lheProduct[2]) , str(lheProduct[3]) )
    ##print process.tagsDumper.LHEEventProduct
    #process.tagsDumper.LHEWeightName = cms.untracked.string(customize.options.WeightName)

print "======================================================> 5a1"
#process.tagsDumper.NNLOPSWeightFile=cms.FileInPath("flashgg/Taggers/data/NNLOPS_reweight.root")
#process.tagsDumper.reweighGGHforNNLOPS = cms.untracked.bool(bool(customize.processId.count("ggh")))
#process.tagsDumper.classifierCfg.remap=cms.untracked.VPSet()

# debugging

print'customize = ',customize 
print'customize.datasetName() = ',customize.datasetName()

#
print "======================================================> 5a2"

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
hlt_paths = []
for dset in customize.metaConditions["TriggerPaths"]:
    print dset
    if not customize.datasetName() : continue
    if dset in customize.datasetName():
      hlt_paths.extend(customize.metaConditions["TriggerPaths"][dset])
print hlt_paths

print "!!!! 1"
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring(hlt_paths))
print "!!!! 2"
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
print "!!!! 3"
process.dataRequirements = cms.Sequence()
if customize.processId == "Data":
        process.dataRequirements += process.hltHighLevel
        
print "======================================================> 6"

# Split out prompt-fake or fake-fake
process.finalFilter = cms.Sequence()
if (customize.processId.count("qcd") or customize.processId.count("gjet")) and customize.processId.count("fake"):
    process.load("flashgg/Systematics/PromptFakeFilter_cfi")
    process.finalFilter += process.PromptFakeFilter
    if (customize.processId.count("promptfake")):
        process.PromptFakeFilter.doPromptFake = cms.bool(True)
        process.PromptFakeFilter.doFakeFake =cms.bool(False)
    elif (customize.processId.count("fakefake")):
        process.PromptFakeFilter.doPromptFake =cms.bool(False)
        process.PromptFakeFilter.doFakeFake =cms.bool(True)
    else:
        raise Exception,"Mis-configuration of python for prompt-fake filter"

# Met Filters
process.load('flashgg/Systematics/flashggMetFilters_cfi')

if customize.processId == "Data":
    metFilterSelector = "data"
    filtersInputTag = cms.InputTag("TriggerResults", "", "RECO")
else:
    metFilterSelector = "mc"
    filtersInputTag = cms.InputTag("TriggerResults", "", "PAT")

process.flashggMetFilters.requiredFilterNames = cms.untracked.vstring([filter.encode("ascii") for filter in customize.metaConditions["flashggMetFilters"][metFilterSelector]])
process.flashggMetFilters.filtersInputTag = filtersInputTag

# Split WH and ZH
process.genFilter = cms.Sequence()
if ((customize.processId.count("wh") or customize.processId.count("zh")) and not (customize.processId.count("powheg"))) and not customize.processId.count("wzh") :
    print "enabling vh filter!!!!!"
    process.load("flashgg/Systematics/VHFilter_cfi")
    process.genFilter += process.VHFilter
    process.VHFilter.chooseW = bool(customize.processId.count("wh"))
    process.VHFilter.chooseZ = bool(customize.processId.count("zh"))

# HHWWggTagsOnly requires zeroeth vertex, but not modifySystematicsWorkflowForttH
if customize.HHWWggTagsOnly or True:
    process.content = cms.EDAnalyzer("EventContentAnalyzer")
    process.p = cms.Path(process.dataRequirements*
                         process.flashggMetFilters*
                         process.genFilter*
                         process.flashggDiPhotons* # needed for 0th vertex from microAOD
                         process.flashggDifferentialPhoIdInputsCorrection*
                         process.flashggUnpackedJets*
                         process.flashggPrefireDiPhotons*
                         process.flashggPreselectedDiPhotons*
                         process.flashggDiPhotonSystematics*
                         process.flashggDiPhotonMVA*
                         #process.flashggMetSystematics*
                         #process.flashggMuonSystematics*process.flashggElectronSystematics*
                         #*process.jetSystematicsSequence)*
                         #process.content* 
                         process.flashggTagSequence)

if customize.doBJetRegression:
    print "DO bjets regressin ====================<<<<<<<<<"
    bregProducers = []
    doubleHTagProducers = []
    
    from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag
    from flashgg.Taggers.flashggbRegressionProducer_cfi import flashggbRegressionProducer
    recoJetCollections = UnpackedJetCollectionVInputTag

    jetsysts = cms.vstring()
    jetnames = cms.vstring()
    for jetsyst in [systlabels[0]]+jetsystlabels:
        jetsysts.append(jetsyst)
    for icoll,coll in enumerate(recoJetCollections):
        jetnames.append(coll.moduleLabel)
    producer = flashggbRegressionProducer.clone(JetSuffixes = jetsysts)
    producer.JetNames = jetnames
    producer.bRegressionWeightfile = cms.untracked.string(str(os.environ["CMSSW_BASE"]+customize.metaConditions['bRegression']['weightFile']))
    producer.y_mean = customize.metaConditions['bRegression']['y_mean']
    producer.y_std = customize.metaConditions['bRegression']['y_std']
    producer.year = cms.untracked.string(str(customize.metaConditions['bRegression']['year']))

    setattr(process,"bRegProducer",producer)
    bregProducers.append(producer)
    process.bregProducers = cms.Sequence(reduce(lambda x,y: x+y, bregProducers))
    process.p.replace(process.jetSystematicsSequence,process.jetSystematicsSequence*process.flashggUnpackedJets+process.bregProducers)
    
    
print "======================================================> 6"
if( not hasattr(process,"options") ): process.options = cms.untracked.PSet()
process.options.allowUnscheduled = cms.untracked.bool(True)

print "--- Dumping all modules: ---"
mns = process.p.moduleNames()
for mn in mns:
    module = getattr(process,mn)
    print str(module)

print "--- Dumping modules that take diphotons as input: ---"
mns = process.p.moduleNames()
for mn in mns:
    module = getattr(process,mn)
    if hasattr(module,"src") and type(module.src) == type(cms.InputTag("")) and module.src.value().count("DiPhoton"):
        print str(module),module.src
    elif hasattr(module,"DiPhotonTag"):
        print str(module),module.DiPhotonTag


def printSystematicInfo(process):
    vpsetlist = [process.flashggDiPhotonSystematics.SystMethods]
    print (14*"-"+" DUMPING SYSTEMATIC OVERVIEW "+14*"-")
    print "%20s %15s %20s" % ("Systematic","Central value?","Systematic shifts?")
    print 57*"-"
    printSystematicVPSet(vpsetlist)
    print (13*"-"+" DUMPING 2D SYSTEMATIC OVERVIEW "+12*"-")
    print "%20s %15s %20s" % ("Systematic","Central value?","Systematic shifts?")
    print 57*"-"
    vpsetlist2D  = [process.flashggDiPhotonSystematics.SystMethods2D]
    printSystematicVPSet(vpsetlist2D)
printSystematicInfo(process)


### Rerun microAOD sequence on top of microAODs using the parent dataset ???
if customize.useParentDataset:
    #print "\n\n\n Rerun microAOD sequence on top of microAODs using the parent dataset ??? <<<<<====="
    #runRivetSequence(process, customize.metaConditions)
    recalculatePDFWeights(process, customize.metaConditions)
    process.p.insert(0, process.content)
        
print "======================================================> 7"
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,era='2018-Prompt')  

customize(process)
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

print "Final path:"
print process.p

print "Final di-photon systematic:"
for pset in process.flashggDiPhotonSystematics.SystMethods:
  print  pset.Label.value()

print "======================================================> 7"







