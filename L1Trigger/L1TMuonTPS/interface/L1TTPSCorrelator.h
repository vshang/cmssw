#ifndef L1TTPSCORRELATOR_H
#define L1TTPSCORRELATOR_H

#include "DataFormats/L1TMuon/interface/L1MuCorrelatorHit.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "L1Trigger/L1TMuonTPS/interface/L1TTPSSectorProcessor.h"

class L1TTPSCorrelator {
 public:
  typedef std::vector<edm::Ptr< l1t::L1TkMuonParticle::L1TTTrackType > > TrackPtrVector;

  L1TTPSCorrelator(const edm::ParameterSet&);
  ~L1TTPSCorrelator();
  std::vector<l1t::L1TkMuonParticle> process(const TrackPtrVector& ,const L1MuCorrelatorHitRefVector&);
 private:
  
  std::vector<l1t::L1TkMuonParticle> clean(const std::vector<l1t::L1TkMuonParticle>&,const std::vector<l1t::L1TkMuonParticle>&,const std::vector<l1t::L1TkMuonParticle>&);

  int verbose_;
  std::vector<L1TTPSSectorProcessor> processor_;
};

#endif
