
#include "L1Trigger/L1TMuonTPS/interface/L1TMuCorrelatorCSCStubProcessor.h"
#include "L1Trigger/L1TMuon/interface/MuonTriggerPrimitive.h"
#include <cmath>
#include <iostream> 
#include <string> 
#include <sstream> 

L1TMuCorrelatorCSCStubProcessor::L1TMuCorrelatorCSCStubProcessor():
  minBX_(-3),
  maxBX_(3)
{
 
} 



L1TMuCorrelatorCSCStubProcessor::L1TMuCorrelatorCSCStubProcessor(const edm::ParameterSet& iConfig):
  minBX_(iConfig.getParameter<int>("minBX")),
  maxBX_(iConfig.getParameter<int>("maxBX")),
  phiLSB_(iConfig.getParameter<double>("phiLSB")),
  etaLSB_(iConfig.getParameter<double>("etaLSB"))
{

} 



L1TMuCorrelatorCSCStubProcessor::~L1TMuCorrelatorCSCStubProcessor() {}




L1MuCorrelatorHit
L1TMuCorrelatorCSCStubProcessor::buildStub(const CSCDetId& detid,const CSCCorrelatedLCTDigi& digi, const L1TMuon::GeometryTranslator *translator) {

  int endcap = detid.endcap();
  int station = detid.station();
  int chamber =detid.chamber();
  int ring = detid.ring();

  L1TMuon::TriggerPrimitive primitive(detid,digi);
  const GlobalPoint& gp = translator->getGlobalPoint(primitive);


  int phi = int(gp.phi().value()/phiLSB_);
  int eta1 = int(gp.eta()/etaLSB_);

  int type=2;
  int wheel=0;;
  int sign = endcap==1 ? -1 : 1;

  if (ring==3)
    wheel = sign*3;
  if (ring==2)
    wheel = sign*4;
  if (ring==1)
    wheel = sign*5;

  int sector = fabs(chamber);
  int phiB =0;
  bool tag = 0;
  int bx=digi.getBX()-8;
  int quality=0;

  uint tfLayer=0;
  
  if (ring==1 && station==1)
    tfLayer=3;
  else if (station==1)
    tfLayer=8;
  else if (station==2)
    tfLayer=6;
  else if (station==3)
    tfLayer=7;
  else if (station==4)
    tfLayer=5;







  L1MuCorrelatorHit stub(wheel,sector,station,tfLayer,phi,phiB,tag,
			 bx,quality,eta1,1,type); 
  return stub;

  }






L1MuCorrelatorHitCollection 
L1TMuCorrelatorCSCStubProcessor::makeStubs(const MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi>& digis,const L1TMuon::GeometryTranslator *t) {
  L1MuCorrelatorHitCollection out;
   auto chamber = digis.begin();
   auto chend  = digis.end();
   for( ; chamber != chend; ++chamber ) {
     auto digi = (*chamber).second.first;
     auto dend = (*chamber).second.second;    
     for( ; digi != dend; ++digi ) {
       L1MuCorrelatorHit stub = buildStub((*chamber).first,*digi,t);
       if (stub.bxNum()>=minBX_&&stub.bxNum()<=maxBX_)
	 out.push_back(stub);
     }
   }
  return out;
}

