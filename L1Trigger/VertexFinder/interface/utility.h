#ifndef __L1Trigger_VertexFinder_utility_h__
#define __L1Trigger_VertexFinder_utility_h__


#include <vector>


namespace l1tVertexFinder {

class AnalysisSettings;
class Stub;
class TP;

namespace utility {

// Count number of tracker layers a given list of stubs are in.
//
// By default uses the "reduced" layer ID if the configuration file requested it. However,
// you can insist on "normal" layer ID being used instead, ignoring the configuration file, by
// setting disableReducedLayerID = true.
//
// N.B. The "reduced" layer ID merges some layer IDs, so that no more than 8 ID are needed in any
// eta region, so as to simplify the firmware.
//
// N.B. You should set disableReducedLayerID = false when counting the number of layers on a tracking
// particle or how many layers it shares with an L1 track. Such counts by CMS convention use "normal" layer ID.
//
// By default, considers both PS+2S modules, but optionally considers only the PS ones if onlyPS = true.

unsigned int countLayers(const AnalysisSettings& settings, const std::vector<const Stub*>& stubs, bool disableReducedLayerID = false, bool onlyPS = false);


// Given a set of stubs (presumably on a reconstructed track candidate)
// return the best matching Tracking Particle (if any),
// the number of tracker layers in which one of the stubs matched one from this tracking particle,
// and the list of the subset of the stubs which match those on the tracking particle.

const TP* matchingTP(const AnalysisSettings& settings, const std::vector<Stub*>& vstubs, unsigned int& nMatchedLayersBest, std::vector<const Stub*>& matchedStubsBest);


} // end namespace utility
} // end namespace l1tVertexFinder

#endif
