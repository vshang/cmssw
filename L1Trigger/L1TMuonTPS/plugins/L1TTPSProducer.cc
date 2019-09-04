#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1TMuon/interface/L1MuCorrelatorHit.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
#include "L1Trigger/L1TMuonTPS/interface/L1TTPSCorrelator.h"

//
// class declaration
//

class L1TTPSProducer : public edm::stream::EDProducer<> {
   public:
      explicit L1TTPSProducer(const edm::ParameterSet&);
      ~L1TTPSProducer() override;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      void beginStream(edm::StreamID) override;
      void produce(edm::Event&, const edm::EventSetup&) override;
      void endStream() override;

  edm::EDGetTokenT<std::vector<L1MuCorrelatorHit> > srcStubs_;
  edm::EDGetTokenT<l1t::L1TkMuonParticle::L1TTTrackCollection> srcTracks_;
  std::unique_ptr<L1TTPSCorrelator> correlator_;
  double maxchi2_;
  unsigned int minStubs_;
  
  
  

};
L1TTPSProducer::L1TTPSProducer(const edm::ParameterSet& iConfig):
  srcStubs_(consumes<std::vector<L1MuCorrelatorHit> >(iConfig.getParameter<edm::InputTag>("srcStubs"))),
  srcTracks_(consumes<l1t::L1TkMuonParticle::L1TTTrackCollection>(iConfig.getParameter<edm::InputTag>("srcTracks"))),
  correlator_(new L1TTPSCorrelator(iConfig)),
  maxchi2_(iConfig.getParameter<double>("maxChi2")),
  minStubs_(iConfig.getParameter<unsigned int>("minStubs"))
{
  produces <std::vector<l1t::L1TkMuonParticle> >();
}


L1TTPSProducer::~L1TTPSProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}

 


//
// member functions
//

// ------------ method called to produce the data  ------------
void
L1TTPSProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   Handle<std::vector<L1MuCorrelatorHit> >stubHandle;
   Handle<l1t::L1TkMuonParticle::L1TTTrackCollection> trackHandle;
   iEvent.getByToken(srcStubs_,stubHandle);
   iEvent.getByToken(srcTracks_,trackHandle);
   L1MuCorrelatorHitRefVector stubs;
   for (uint i=0;i<stubHandle->size();++i) {
     L1MuCorrelatorHitRef r(stubHandle,i);
     if (r->bxNum()==0)
       stubs.push_back(r);
   }

   std::vector<edm::Ptr< l1t::L1TkMuonParticle::L1TTTrackType > > tracks;
   for (uint i=0;i<trackHandle->size();++i) {
     edm::Ptr< l1t::L1TkMuonParticle::L1TTTrackType > track(trackHandle, i);
     double chi2dof=track->getChi2()/(2*track->getStubRefs().size()-4);
     if (chi2dof>maxchi2_   || track->getStubRefs().size()<4) 
	continue;
     else
       tracks.push_back(track);
   }

   std::vector<l1t::L1TkMuonParticle> out = correlator_->process(tracks,stubs);
   std::unique_ptr<std::vector<l1t::L1TkMuonParticle> > out1 = std::make_unique<std::vector<l1t::L1TkMuonParticle> >(out); 
   iEvent.put(std::move(out1));

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
L1TTPSProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
L1TTPSProducer::endStream() {
}

void
L1TTPSProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TTPSProducer);
