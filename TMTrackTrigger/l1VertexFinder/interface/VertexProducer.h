#ifndef __TMTrackTrigger_VertexFinder_VertexProducer_h__
#define __TMTrackTrigger_VertexFinder_VertexProducer_h__

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"


#include <vector>
#include <map>
#include <string>


using namespace std;

class Settings;
class Histos;
class TrackFitGeneric;

class VertexProducer : public edm::EDProducer {

public:
  explicit VertexProducer(const edm::ParameterSet&);	
  ~VertexProducer(){}

private:

  typedef std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > TTTrackCollection;

  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

private:
  const edm::EDGetTokenT<TrackingParticleCollection> tpInputTag;
  typedef edmNew::DetSetVector< TTStub<Ref_Phase2TrackerDigi_> > DetSetVec;
  const edm::EDGetTokenT<DetSetVec> stubInputTag;
  typedef TTStubAssociationMap<Ref_Phase2TrackerDigi_>           TTStubAssMap;
  const edm::EDGetTokenT<TTStubAssMap> stubTruthInputTag;
  typedef TTClusterAssociationMap<Ref_Phase2TrackerDigi_>        TTClusterAssMap;
  const edm::EDGetTokenT<TTClusterAssMap> clusterTruthInputTag;
  const edm::EDGetTokenT<TTTrackCollection> l1TracksToken_;

  Settings *settings_;
  Histos   *hists_;
};

#endif

