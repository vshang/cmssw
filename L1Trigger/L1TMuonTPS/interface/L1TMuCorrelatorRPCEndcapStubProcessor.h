//
// Combined correlator Hit producer
//Michalis Bachtis (UCLA)
//

#ifndef L1TMUCORRELATORRPCENDCAPSTUBPROCESSOR
#define L1TMUCORRELATORRPCENDCAPSTUBPROCESSOR

#include "DataFormats/L1TMuon/interface/L1MuCorrelatorHit.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "L1Trigger/L1TMuon/interface/GeometryTranslator.h"

class L1TMuCorrelatorRPCEndcapStubProcessor {


 public:
  L1TMuCorrelatorRPCEndcapStubProcessor();
  L1TMuCorrelatorRPCEndcapStubProcessor(const edm::ParameterSet&);
    ~L1TMuCorrelatorRPCEndcapStubProcessor();


    L1MuCorrelatorHitCollection makeStubs(const MuonDigiCollection<RPCDetId,RPCDigi>&,const L1TMuon::GeometryTranslator*,const edm::EventSetup& iSetup);

    
 private:
    L1MuCorrelatorHit buildStub(const RPCDetId& , const RPCDigi&,const L1TMuon::GeometryTranslator*);
  int minBX_;
  int maxBX_;
  double phiLSB_;
  double etaLSB_;

};


#endif
