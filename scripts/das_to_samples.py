

def das_to_samples(txt):
  print "============================================="
  datasets = []
  for line in txt.split("\n"):
    if not "Dataset:" in line : continue
    datasets += [ line.split()[-1] ]

  for dataset in sorted(datasets):
    print dataset

das = """Dataset: /SinglePhoton/Run2016H-17Jul2018-v1/MINIAOD
Creation time: 2018-07-24 04:06:22 Physics group: NoGroup Status: VALID Type: data
Release, Blocks, Files, Runs, Configs, Parents, Children, Sites, Physics Groups XSDB Sources: dbs3 show Dataset: /SinglePhoton/Run2016G-17Jul2018-v1/MINIAOD
Creation time: 2018-07-18 18:11:01 Physics group: NoGroup Status: VALID Type: data
Release, Blocks, Files, Runs, Configs, Parents, Children, Sites, Physics Groups XSDB Sources: dbs3 show Dataset: /SinglePhoton/Run2016F-17Jul2018-v1/MINIAOD
Creation time: 2018-07-21 12:33:47 Physics group: NoGroup Status: VALID Type: data
Release, Blocks, Files, Runs, Configs, Parents, Children, Sites, Physics Groups XSDB Sources: dbs3 show Dataset: /SinglePhoton/Run2016E-17Jul2018-v1/MINIAOD
Creation time: 2018-07-18 05:38:33 Physics group: NoGroup Status: VALID Type: data
Release, Blocks, Files, Runs, Configs, Parents, Children, Sites, Physics Groups XSDB Sources: dbs3 show Dataset: /SinglePhoton/Run2016D-17Jul2018-v1/MINIAOD
Creation time: 2018-07-19 07:17:43 Physics group: NoGroup Status: VALID Type: data
Release, Blocks, Files, Runs, Configs, Parents, Children, Sites, Physics Groups XSDB Sources: dbs3 show Dataset: /SinglePhoton/Run2016C-17Jul2018-v1/MINIAOD
Creation time: 2018-07-18 15:17:43 Physics group: NoGroup Status: VALID Type: data
Release, Blocks, Files, Runs, Configs, Parents, Children, Sites, Physics Groups XSDB Sources: dbs3 show Dataset: /SinglePhoton/Run2016B-17Jul2018_ver2-v1/MINIAOD
Creation time: 2018-07-19 07:17:43 Physics group: NoGroup Status: VALID Type: data
Release, Blocks, Files, Runs, Configs, Parents, Children, Sites, Physics Groups XSDB Sources: dbs3 show Dataset: /SinglePhoton/Run2016B-17Jul2018_ver1-v1/MINIAOD
Creation time: 2018-08-02 01:04:57 Physics group: NoGroup Status: VALID Type: data
Release, Blocks, Files, Runs, Configs, Parents, Children, Sites, Physics Groups XSDB Sources: dbs3 show"""
das_to_samples(das)

das = """Dataset: /SinglePhoton/Run2017F-31Mar2018-v1/MINIAOD
Creation time: 2018-04-03 08:11:45 Physics group: NoGroup Status: VALID Type: data
Release, Blocks, Files, Runs, Configs, Parents, Children, Sites, Physics Groups XSDB Sources: dbs3 show Dataset: /SinglePhoton/Run2017E-31Mar2018-v1/MINIAOD
Creation time: 2018-04-03 05:08:12 Physics group: NoGroup Status: VALID Type: data
Release, Blocks, Files, Runs, Configs, Parents, Children, Sites, Physics Groups XSDB Sources: dbs3 show Dataset: /SinglePhoton/Run2017D-31Mar2018-v1/MINIAOD
Creation time: 2018-04-03 13:01:04 Physics group: NoGroup Status: VALID Type: data
Release, Blocks, Files, Runs, Configs, Parents, Children, Sites, Physics Groups XSDB Sources: dbs3 show Dataset: /SinglePhoton/Run2017C-31Mar2018-v1/MINIAOD
Creation time: 2018-04-04 08:02:36 Physics group: NoGroup Status: VALID Type: data
Release, Blocks, Files, Runs, Configs, Parents, Children, Sites, Physics Groups XSDB Sources: dbs3 show Dataset: /SinglePhoton/Run2017B-31Mar2018-v1/MINIAOD
Creation time: 2018-04-03 03:51:52 Physics group: NoGroup Status: VALID Type: data
Release, Blocks, Files, Runs, Configs, Parents, Children, Sites, Physics Groups XSDB Sources: dbs3 show"""
das_to_samples(das)


#####################################
from cmssw_das_client import *
def das_to_samples_v2( datasets_pattern ):
  print "============================================= v2"
  datasets = []
  answer = get_data( "dataset dataset=" + datasets_pattern )

  try: 
    data = answer['data']
  except : 
    print answer

  for item in answer['data']:
    datasets += [ item['dataset'][0]['name'] ]
  for dataset in sorted(datasets):
    print dataset

if False:
  # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable
  # https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
  print "\n\n\n DATA"
  das_to_samples_v2("/DoubleEG/Run2016*-17Jul2018*/MINIAOD")
  das_to_samples_v2("/DoubleEG/Run2017*-31Mar2018*/MINIAOD")
  das_to_samples_v2("/EGamma/Run2018*-17Sep2018*/MINIAOD")
  das_to_samples_v2("/EGamma/Run2018*-22Jan2019*/MINIAOD")

if False:
  # https://twiki.cern.ch/twiki/bin/view/CMS/Higgs2GGen
  print "\n\n\n MC Higgs -> yy", "ggF VBF VH ttH"
  print "2016 =========================================="
  das_to_samples_v2("/GluGluHToGG*125*13TeV*/RunIISummer16MiniAODv3-*/MINIAODSIM") # SM ggfusion
  das_to_samples_v2("/bbHToGG*125*13TeV*/RunIISummer16MiniAODv3-*/MINIAODSIM") # ttH
  das_to_samples_v2("/ttHJetToGG*125*13TeV*/RunIISummer16MiniAODv3-*/MINIAODSIM") # bbH https://arxiv.org/pdf/1409.5301.pdf
  das_to_samples_v2("/VBFHToGG*M125*13TeV*/RunIISummer16MiniAODv3-*/MINIAODSIM") # VBF
  das_to_samples_v2("/VHToGG*M125*13TeV*/RunIISummer16MiniAODv3-*/MINIAODSIM") # VH

  print "2017 =========================================="
  das_to_samples_v2("/GluGluHToGG*125*13TeV*/RunIIFall17MiniAODv2-*/MINIAODSIM") # SM ggfusion
  das_to_samples_v2("/bbHToGG*125*13TeV*/RunIIFall17MiniAODv2-*/MINIAODSIM") # ttH
  das_to_samples_v2("/ttHJetToGG*125*13TeV*/RunIIFall17MiniAODv2-*/MINIAODSIM") # bbH https://arxiv.org/pdf/1409.5301.pdf
  das_to_samples_v2("/VBFHToGG*M125*13TeV*/RunIIFall17MiniAODv2-*/MINIAODSIM") # VBF
  das_to_samples_v2("/VHToGG*M125*13TeV*/RunIIFall17MiniAODv2-*/MINIAODSIM") # VH

  print "2018 =========================================="
  das_to_samples_v2("/GluGluHToGG*125*13TeV*/RunIIAutumn18MiniAOD-*/MINIAODSIM") # SM ggfusion
  das_to_samples_v2("/bbHToGG*125*13TeV*/RunIIAutumn18MiniAOD-*/MINIAODSIM") # ttH
  das_to_samples_v2("/ttHJetToGG*125*13TeV*/RunIIAutumn18MiniAOD-*/MINIAODSIM") # bbH https://arxiv.org/pdf/1409.5301.pdf
  das_to_samples_v2("/VBFHToGG*M125*13TeV*/RunIIAutumn18MiniAOD-*/MINIAODSIM") # VBF
  das_to_samples_v2("/VHToGG*M125*13TeV*/RunIIAutumn18MiniAOD-*/MINIAODSIM") # VH


if True:
  das_to_samples_v2("/*HHTo*2G*/*/MINIAODSIM")





























