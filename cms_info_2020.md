#### RECOMENDED VERSIONS OF THE CMSSW
Take from tables:  
https://twiki.cern.ch/twiki/bin/view/CMS/PdmV  
https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable  

**-> CMSSW_10_2_18 <-** for all years!!!

Then do 
```shell
cmsrel CMSSW_10_2_18
mkdir -p CMSSW_10_2_18/src/Analysis; cd CMSSW_10_2_18/src/Analysis;
git clone https://github.com/pmandrik/GRINDER.git GRINDER; cd GRINDER
cmsnev
scram b
```

