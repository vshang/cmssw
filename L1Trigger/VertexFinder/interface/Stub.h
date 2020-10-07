#ifndef __L1Trigger_VertexFinder_Stub_h__
#define __L1Trigger_VertexFinder_Stub_h__


#include <array>
#include <map>
#include <set>

#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"

// TTStubAssociationMap.h forgets to two needed files, so must include them here ...
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"



class TrackerGeometry;

namespace l1tVertexFinder {

class AnalysisSettings;
class TP;

typedef edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_>> DetSetVec;
typedef edmNew::DetSet<TTStub<Ref_Phase2TrackerDigi_>> DetSet;
typedef edm::Ref<DetSetVec, TTStub<Ref_Phase2TrackerDigi_>> TTStubRef;
typedef edm::Ref<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_>>, TTCluster<Ref_Phase2TrackerDigi_>> TTClusterRef;
typedef TTStubAssociationMap<Ref_Phase2TrackerDigi_> TTStubAssMap;
typedef TTClusterAssociationMap<Ref_Phase2TrackerDigi_> TTClusterAssMap;


//=== Represents a Tracker stub (=pair of hits)

class Stub : public TTStubRef {

public:
  // Fill useful info about stub.
  Stub(const TTStubRef& ttStubRef,
       unsigned int index_in_vStubs,
       const AnalysisSettings& settings,
       const TrackerGeometry* trackerGeometry,
       const TrackerTopology* trackerTopology,
       const std::map<DetId, DetId>& geoDetIdMap);
  ~Stub() {}

  // Fill truth info with association from stub to tracking particles.
  // The 1st argument is a map relating TrackingParticles to TP.
  void fillTruth(const std::map<edm::Ptr<TrackingParticle>, const TP*>& translateTP, edm::Handle<TTStubAssMap> mcTruthTTStubHandle, edm::Handle<TTClusterAssMap> mcTruthTTClusterHandle);


  // === Functions for returning info about reconstructed stubs ===

  //--- Stub data and quantities derived from it ---

  // Stub coordinates (optionally after digitisation, if digitisation requested via cfg).
  // N.B. Digitisation is not run when the stubs are created, but later, after stubs are assigned to sectors.
  // Until then, these functions return the original coordinates.
  float phi() const { return phi_; }
  float r() const { return r_; }
  float z() const { return z_; }

  //--- Quantities common to all stubs in a given module ---

  // Module type: PS or 2S?
  bool psModule() const { return psModule_; }
  // Tracker layer ID number (1-6 = barrel layer; 11-15 = endcap A disk; 21-25 = endcap B disk)
  unsigned int layerId() const { return layerId_; }
  // Reduced layer ID (in range 1-7). This encodes the layer ID in only 3 bits (to simplify firmware) by merging some barrel layer and endcap disk layer IDs into a single ID.
  unsigned int layerIdReduced() const;

  //--- Truth info

  // Association of stub to tracking particles
  const std::set<const TP*>& assocTPs() const
  {
    return assocTPs_;
  }                                              // Return TPs associated to this stub. (Whether only TPs contributing to both clusters are returned is determined by "StubMatchStrict" config param.)
  const TP* assocTP() const { return assocTP_; } // If only one TP contributed to both clusters, this tells you which TP it is. Returns nullptr if none.

private:
  // Set info about the module that this stub is in.
  void setModuleInfo(const TrackerGeometry* trackerGeometry, const TrackerTopology* trackerTopology, const DetId& detId);

private:
  const AnalysisSettings* settings_; // configuration parameters.

  float phi_; // stub coords, optionally after digitisation.
  float r_;
  float z_;

  //--- Parameters common to all stubs in a given module.
  unsigned int idDet_;
  float moduleMinR_;
  float moduleMaxR_;
  float moduleMinPhi_;
  float moduleMaxPhi_;
  float moduleMinZ_;
  float moduleMaxZ_;
  bool psModule_;
  unsigned int layerId_;
  unsigned int endcapRing_;
  bool barrel_;
  float sigmaPerp_;
  float sigmaPar_;
  float stripPitch_;
  float stripLength_;
  unsigned int nStrips_;
  float sensorWidth_;
  bool outerModuleAtSmallerR_;
  //--- Truth info about stub.
  const TP* assocTP_;
  std::set<const TP*> assocTPs_;
  //--- Truth info about the two clusters that make up the stub
  std::array<const TP*, 2> assocTPofCluster_;
};

} // end namespace l1tVertexFinder

#endif
