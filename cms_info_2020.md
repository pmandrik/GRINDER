#### RECOMENDED VERSIONS OF THE CMSSW
Take from tables:  
https://twiki.cern.ch/twiki/bin/view/CMS/PdmV  
https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable  

**-> CMSSW_10_2_18 <-** for all years!!!

Then do:
**-> 2016 <-**
```shell
cmsrel CMSSW_10_2_18
cd CMSSW_10_2_18/src
cmsnev
git cms-init
# Extra packages for Photon ID from https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#102X
git cms-merge-topic cms-egamma:PhotonIDValueMapSpeedup1029
git cms-merge-topic cms-egamma:EgammaPostRecoTools
scram b -j 8
# extra packages for MET filters from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
git cms-addpkg RecoMET/METFilters
# GRIDNER packages
mkdir -p Analysis; cd Analysis;
git clone git@github.com:pmandrik/GRINDER.git GRINDER; cd GRINDER
scram b
# Extra data
cd data; ./download.sh; cd -
```

**-> 2017 <-**
```shell
cmsrel CMSSW_10_2_18
cd CMSSW_10_2_18/src
cmsnev
git cms-init
# Extra packages for Photon ID from https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#102XTODO
# extra packages for MET filters from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
git cms-addpkg RecoMET/METFilters
# GRIDNER packages
mkdir -p Analysis; cd Analysis;
git clone git@github.com:pmandrik/GRINDER.git GRINDER; cd GRINDER
scram b
# Extra data
cd data; ./download.sh; cd -
```

**-> 2018 <-**
```shell
cmsrel CMSSW_10_2_18
cd CMSSW_10_2_18/src
cmsnev
git cms-init
# Extra packages for Photon ID from https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#102XTODO
# extra packages for MET filters from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
git cms-addpkg RecoMET/METFilters
# GRIDNER packages
mkdir -p Analysis; cd Analysis;
git clone git@github.com:pmandrik/GRINDER.git GRINDER; cd GRINDER
scram b
# Extra data
cd data; ./download.sh; cd -
```

#### RUNNING IN THE GRID
Use *scripts/create_crab_cfg.py* in order to create crab cfg files.  













