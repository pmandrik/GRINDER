
import os, sys
sys.path.insert(0, "/afs/cern.ch/user/p/pmandrik/public/PMANDRIK_LIBRARY")

def create_cfg(path, cfg_dic) :
  #os.system( "rm -rvf %s/* " % ( wdir ) )
  config_file_name = os.path.expandvars( path )
  config_file = open(config_file_name, 'w')
  
  # info about options https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile
  cfg_template = """
import os
from WMCore.Configuration import Configuration

config = Configuration()
config.section_("General")
config.General.requestName = {%% request_name  %%}
config.General.transferOutputs = True
config.section_("JobType")
config.JobType.pluginName = 'Analysis' # Specifies if this task is running an analysis ('Analysis') on an existing dataset or is running MC event generation ('PrivateMC')
config.JobType.psetName = os.path.expandvars( {%% input_run_cfg  %%} )
config.JobType.pyCfgParams = {% input_run_cfg_options  %}
config.JobType.allowUndistributedCMSSW = True
config.section_("Data")
config.Data.inputDataset = {%% dataset %%}
config.Data.inputDBS = 'global'
{% if is_data %}
config.Data.splitting = "LumiBased"
config.Data.unitsPerJob = {% units_per_job %}
config.Data.lumiMask = {%% lumi_mask %%}
{% else %}
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = {% units_per_job %}
{% endif %}
config.section_("Site")
config.Site.storageSite = "T2_RU_JINR"
"""

  import template_master as tm
  cfg_text = tm.parce_template( cfg_template, cfg_dic )

  print cfg_dic
  print cfg_text
  
  config_file.write( cfg_text )
  config_file.close()
  return cfg_text


def main() :
  # lumi_mask
  # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable
  # /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt
  # /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt
  # /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt
  def_lumi    = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt"
  def_dataset = "/DoubleEG/Run2016D-17Jul2018-v1/MINIAOD"

  import argparse
  parser = argparse.ArgumentParser(description='Create CRAB datacard for CMSSW job submition to GRID')
  parser.add_argument('--cfg_path',       default="${CMSSW_BASE}/src/Analysis/GRINDER/test/test_crab_cfg.py")
  parser.add_argument('--dataset',        default=def_dataset)
  parser.add_argument('--lumi_mask',      default=def_lumi)
  parser.add_argument('--input_run_cfg',  default="${CMSSW_BASE}/src/Analysis/GRINDER/python/runCfg.py") 
  parser.add_argument('--input_run_cfg_options',  default="['isData=True', 'yearEra=2016' ]") # [ dummy_key : \"dummy_value\"]"
  parser.add_argument('--prefix',   default="test_run_data2016D_v4")
  parser.add_argument('--is_data',  default=True)
  parser.add_argument('--units_per_job',  default=-1)

  cfg_dic = vars(parser.parse_args())
  cfg_dic["request_name"] = cfg_dic["prefix"]
  if cfg_dic["units_per_job"] < 0 : 
    if cfg_dic["is_data"] : cfg_dic["units_per_job"] = 400
    else :                  cfg_dic["units_per_job"] = 5
  
  create_cfg( cfg_dic["cfg_path"], cfg_dic )

if __name__ == "__main__" : sys.exit( main() )





