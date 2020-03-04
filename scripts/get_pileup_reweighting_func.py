
# https://github.com/pfs/BJetEnergyPeak/blob/master/scripts/runPileupEstimation.py

import optparse
import os,sys
import pickle
import commands
import ROOT

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#Pile_up
MBXSEC=69000.
PUSCENARIOS={'nom':MBXSEC,'up':MBXSEC*(1.0 + 0.046),'down':MBXSEC*(1.0 - 0.046)}

"""
Be prepeared (!) using:
  cmsenv
  git cms-addpkg SimGeneral/MixingModule
maybe get some modules from https://github.com/cms-sw/cmssw/tree/master/SimGeneral/MixingModule/python
and run with
  python get_pileup_reweighting_func.py --year 2016
"""


def main():

    # configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--year',      
                      dest='year'  ,      
                      help='year',      
                      default='2016',
                      type='string')
    parser.add_option('--json',      
                      dest='inJson'  ,      
                      help='json file with processed runs',      
                      default='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt',
                      type='string')
    parser.add_option('--puJson',    
                      dest='puJson'  ,
                      help='pileup json file',      
                      default='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt',    
                      type='string')
    parser.add_option('--mixer',    
                      dest='mixer'  ,
                      help='mixer module',      
                      default='mix_2015_25ns_Startup_PoissonOOTPU_cfi',    
                      type='string')
    (opt, args) = parser.parse_args()

    # lumi_mask
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable
    # /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt
    # /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt
    # /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt

    if opt.year == "2016" :
      opt.inJson = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt"
      opt.puJson = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt"
      opt.mixer  = "mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi"
    if opt.year == "2017" :
      opt.inJson = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt"
      opt.puJson = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt"
      opt.mixer  = "mix_2017_25ns_UltraLegacy_PoissonOOTPU_cfi" # FIXME
    if opt.year == "2018" :
      opt.inJson = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt"
      opt.puJson = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt"
      opt.mixer  = "mix_2018_25ns_UltraLegacy_PoissonOOTPU_cfi" # FIXME

    print "loading pileup mixer module ... "
    import importlib
    mixer_module_name = "SimGeneral.MixingModule." + opt.mixer
    #mixer_module_name = "SimGeneral.MixingModule.MixingModule" #  + opt.mixer
    mixer_module      = importlib.import_module( mixer_module_name )
    print "ok ... "
    # simulated pileup
  
    NPUBINS=len(mixer_module.mix.input.nbPileupEvents.probValue)
    simPuH=ROOT.TH1F('simPuH','',NPUBINS,float(0),float(NPUBINS))
    for xbin in xrange(0,NPUBINS):
        probVal=mixer_module.mix.input.nbPileupEvents.probValue[xbin]
        simPuH.SetBinContent(xbin,probVal)
    simPuH.Scale(1./simPuH.Integral())

    # compute pileup in data assuming different xsec
    puWgts,puDists={},{}
    for scenario in PUSCENARIOS:
        print scenario, 'xsec=',PUSCENARIOS[scenario]
        cmd='pileupCalc.py -i %s --inputLumiJSON %s --calcMode true --minBiasXsec %f --maxPileupBin %d --numPileupBins %s Pileup.root'%(opt.inJson,opt.puJson,PUSCENARIOS[scenario],NPUBINS,NPUBINS)
        commands.getstatusoutput(cmd)

        fIn=ROOT.TFile.Open('Pileup.root')
        pileupH=fIn.Get('pileup')
        pileupH.Scale(1./pileupH.Integral())
        puDists[scenario]=ROOT.TGraph(pileupH)
        puDists[scenario].SetName('pu_'+scenario)

        pileupH.Divide(simPuH)
        puWgts[scenario]=ROOT.TGraph(pileupH)
        puWgts[scenario].SetName('puwgts_'+scenario)
        fIn.Close()
    commands.getstatusoutput('rm Pileup.root')

    # dump to pickle file                   
    cache='puweights_' + opt.year + '.root'
    cachefile=ROOT.TFile.Open(cache,'RECREATE')
    for scenario in puWgts:
        puWgts[scenario].Write()
        puDists[scenario].Write()
    cachefile.Close()
    print 'Produced normalization cache for pileup weights @ %s'%cache

"""
for execution from another script
"""

if __name__ == "__main__": main()




