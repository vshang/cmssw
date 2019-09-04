//
// Combined correlator Hit producer
//Michalis Bachtis (UCLA)
//

#ifndef L1TMUCORRELATORCSCSTUBPROCESSOR
#define L1TMUCORRELATORCSCSTUBPROCESSOR

#include "DataFormats/L1TMuon/interface/L1MuCorrelatorHit.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "L1Trigger/L1TMuon/interface/GeometryTranslator.h"

class L1TMuCorrelatorCSCStubProcessor {


 public:
  L1TMuCorrelatorCSCStubProcessor();
  L1TMuCorrelatorCSCStubProcessor(const edm::ParameterSet&);
    ~L1TMuCorrelatorCSCStubProcessor();


    L1MuCorrelatorHitCollection makeStubs(const MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi>&,const L1TMuon::GeometryTranslator*);

    
 private:
    L1MuCorrelatorHit buildStub(const CSCDetId& , const CSCCorrelatedLCTDigi&,const L1TMuon::GeometryTranslator*);
  int minBX_;
  int maxBX_;
  double phiLSB_;
  double etaLSB_;



};


#endif
