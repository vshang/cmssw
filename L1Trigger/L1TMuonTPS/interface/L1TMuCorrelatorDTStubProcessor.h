#ifndef L1TMUCORRELATORDTSTUBPROCESSOR
#define L1TMUCORRELATORDTSTUBPROCESSOR

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/L1TMuon/interface/L1MuCorrelatorHit.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondFormats/L1TObjects/interface/L1TMuonBarrelParams.h"
#include "CondFormats/DataRecord/interface/L1TMuonBarrelParamsRcd.h"

class L1MuDTTFMasks;


class L1TMuCorrelatorDTStubProcessor {
 public:
  L1TMuCorrelatorDTStubProcessor();
  L1TMuCorrelatorDTStubProcessor(const edm::ParameterSet&);
  
  ~L1TMuCorrelatorDTStubProcessor();


  L1MuCorrelatorHitCollection makeStubs(const L1MuDTChambPhContainer*,const L1MuDTChambThContainer*,const L1TMuonBarrelParams&);
  
 private:
  bool isGoodPhiStub(const L1MuDTChambPhDigi*); 
  L1MuCorrelatorHit buildStub(const L1MuDTChambPhDigi&,const L1MuDTChambThDigi*);
  L1MuCorrelatorHit buildStubNoEta(const L1MuDTChambPhDigi&);

  int calculateEta(uint, int,uint,uint);  
  int minPhiQuality_;
  int minBX_;
  int maxBX_;
  std::vector<int> eta1_;
  std::vector<int> eta2_;
  std::vector<int> eta3_;
  std::vector<int> coarseEta1_;
  std::vector<int> coarseEta2_;
  std::vector<int> coarseEta3_;
  std::vector<int> coarseEta4_;
  int verbose_;
  double phiLSB_;
  std::vector<double> bendingScale_;
  std::vector<int> phiOffset_;

  //    edm::ESHandle< L1TMuonBarrelParams > bmtfParamsHandle;
  //    L1MuDTTFMasks       masks_;





};


#endif
