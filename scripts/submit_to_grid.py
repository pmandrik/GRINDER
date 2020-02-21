
import sys
sys.path.insert(0, "/afs/cern.ch/user/p/pmandrik/public/PMANDRIK_LIBRARY")

import template_master as tm

def create_cfg(wdir, prefix, dataset, options):
  #os.system( "rm -rvf %s/* " % ( wdir ) )
  #config_file_name = wdir+'/cfg_' + prefix + '.py'
  #config_file = open(config_file_name, 'w')
  
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
config.JobType.psetName = {%% input_run_cfg  %%}
config.JobType.pyCfgParams = {% input_run_cfg_options  %}
config.section_("Data")
config.Data.inputDataset = {%% dataset %%}
config.Data.inputDBS = 'global'
{% if is_DATA %}
{% else %}
{% fi %}
"""

  cfg_dic = {}
  cfg_dic["request_name"] = prefix
  cfg_dic["input_run_cfg"] = "${CMSSW_BASE}/src/Analysis/GRINDER/runCfg.py"
  cfg_dic["input_run_cfg_options"] = options # "[ dummy_key : \"dummy_value\" ]"
  cfg_dic["dataset"] = dataset

  cfg_text = tm.parce_template( cfg_template, cfg_dic )
  print cfg_text



create_cfg( "path", "prefix", "[ 'dummy_key=True', 'dummy_key2=-1' ]" )







