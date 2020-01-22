
// #include <PereFormat/interface/Muon.h>
// #include <PereFormat/interface/GenParticle.h>
// #include <PereFormat/interface/GenInfo.h>
#include "Analysis/GRINDER/interface/Event.hh"

#include <vector>

// Instantiate templates
template class std::vector<grinder::Muon>;
template class std::vector<grinder::Electron>;
template class std::vector<grinder::Photon>;
template class std::vector<grinder::Jet>;
class Event;
