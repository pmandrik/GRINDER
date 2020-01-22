
import FWCore.ParameterSet.Config as cms

### Options
DEBUG_print_content = False
###

process = cms.Process("Demo")
process.path = cms.Path()
process.TFileService = cms.Service('TFileService', fileName = cms.string('grinder_out.root'))

process.load("FWCore.MessageService.MessageLogger_cfi")

### Input
process.source           = cms.Source("PoolSource")
process.source.fileNames = cms.untracked.vstring('file:/eos/cms/store/data/Run2016H/DoubleMuon/MINIAOD/03Feb2017_ver3-v1/50000/8643C759-9BEB-E611-A7CA-008CFA111270.root')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

### Extra paths
if DEBUG_print_content:
  process.content = cms.EDAnalyzer("EventContentAnalyzer")
  process.path += process.content

### Main path
process.grinderMain = cms.EDAnalyzer('Grinder',
  muons_token = cms.InputTag('slimmedMuons'),
  photons_token = cms.InputTag('slimmedPhotons'),
  jets_token = cms.InputTag('slimmedJets'),
  mets_token = cms.InputTag('slimmedMETs'),
  electrons_token = cms.InputTag('slimmedElectrons'),
  taus_token = cms.InputTag('slimmedTaus'),

  primaryVertex_token = cms.InputTag('offlineSlimmedPrimaryVertices'),
)
process.path += process.grinderMain
