
#ifndef Grinder_extensions_hh
#define Grinder_extensions_hh 1

namespace grinder{
  
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
  // https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
  const reco::GenJet * MatchGenJet(pat::Jet const &jet, edm::Handle<edm::View<reco::GenJet>> const &genJets, double maxDPt) {
    reco::GenJet const *matchedJet = nullptr;
    double minDR2 = std::numeric_limits<double>::infinity();
    double const maxDR2 = 0.04; // jetConeSize * jetConeSize / 4.;
      
    for (unsigned i = 0; i < genJets->size(); ++i){
      const reco::GenJet & genJet = genJets->at(i);
      double const dR2 = ROOT::Math::VectorUtil::DeltaR2(jet.p4(), genJet.p4());
      if (dR2 > maxDR2 or dR2 > minDR2) continue;
      if (std::abs(jet.pt() - genJet.pt()) > maxDPt) continue;
      minDR2 = dR2;
      matchedJet = & ( genJets->at(i) );
    }
    return matchedJet;
  }


};

#endif











