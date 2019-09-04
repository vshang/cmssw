#include "L1Trigger/L1TMuonTPS/interface/L1TMuCorrelatorRPCEndcapStubProcessor.h"
#include "L1Trigger/L1TMuon/interface/MuonTriggerPrimitive.h"
#include "L1Trigger/L1TTwinMux/interface/RPCHitCleaner.h"
#include <cmath>
#include <iostream> 
#include <string> 
#include <sstream> 

L1TMuCorrelatorRPCEndcapStubProcessor::L1TMuCorrelatorRPCEndcapStubProcessor():
  minBX_(-3),
  maxBX_(3)
{
 
} 



L1TMuCorrelatorRPCEndcapStubProcessor::L1TMuCorrelatorRPCEndcapStubProcessor(const edm::ParameterSet& iConfig):
  minBX_(iConfig.getParameter<int>("minBX")),
  maxBX_(iConfig.getParameter<int>("maxBX")),
  phiLSB_(iConfig.getParameter<double>("phiLSB")),
  etaLSB_(iConfig.getParameter<double>("etaLSB"))
{

} 



L1TMuCorrelatorRPCEndcapStubProcessor::~L1TMuCorrelatorRPCEndcapStubProcessor() {}




L1MuCorrelatorHit
L1TMuCorrelatorRPCEndcapStubProcessor::buildStub(const RPCDetId& detid,const RPCDigi& digi, const L1TMuon::GeometryTranslator *translator) {


  L1TMuon::TriggerPrimitive primitive(detid,digi);
  const GlobalPoint& gp = translator->getGlobalPoint(primitive);


  int phi = int(gp.phi().value()/phiLSB_);
  int eta1 = int(gp.eta()/etaLSB_);

  int type=3;
  int wheel=-(6-detid.ring())*detid.region();
  int sector =(detid.sector()-1)*6+detid.subsector();
  int station=detid.station();
  int phiB =0;
  bool tag = detid.trIndex();
  int bx=digi.bx();
  int quality=0;
  uint tfLayer=0;
  
  if (station==1)
    tfLayer=2;
  if (station==2)
    tfLayer=9;
  if (station==3)
    tfLayer=4;
  if (station==4)
    tfLayer=10;





  L1MuCorrelatorHit stub(wheel,sector,station,tfLayer,phi,phiB,tag,
			    bx,quality,eta1,0,type); 
  return stub;

  }






L1MuCorrelatorHitCollection 
L1TMuCorrelatorRPCEndcapStubProcessor::makeStubs(const MuonDigiCollection<RPCDetId,RPCDigi>& digis,const L1TMuon::GeometryTranslator *t,const edm::EventSetup& iSetup) {
  L1MuCorrelatorHitCollection out;

  RPCHitCleaner cleaner(digis);
  cleaner.run(iSetup,true);    
  RPCDigiCollection cleaned = cleaner.getRPCCollection();

   auto chamber = cleaned.begin();
   auto chend  = cleaned.end();
   for( ; chamber != chend; ++chamber ) {
       if((*chamber).first.region()==0)
	 continue;
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

