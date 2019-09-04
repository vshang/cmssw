#include "L1Trigger/L1TMuonTPS/interface/L1TMuCorrelatorDTStubProcessor.h"
#include <cmath>
#include "CondFormats/L1TObjects/interface/L1MuDTTFParameters.h"
#include "CondFormats/DataRecord/interface/L1MuDTTFParametersRcd.h"
#include "CondFormats/L1TObjects/interface/L1MuDTTFMasks.h"
#include "CondFormats/DataRecord/interface/L1MuDTTFMasksRcd.h"

#include <iostream> 
#include <string> 
#include <sstream> 

L1TMuCorrelatorDTStubProcessor::L1TMuCorrelatorDTStubProcessor():
  minPhiQuality_(0),
  minBX_(-3),
  maxBX_(3)
{
 
} 



L1TMuCorrelatorDTStubProcessor::L1TMuCorrelatorDTStubProcessor(const edm::ParameterSet& iConfig):
  minPhiQuality_(iConfig.getParameter<int>("minPhiQuality")),
  minBX_(iConfig.getParameter<int>("minBX")),
  maxBX_(iConfig.getParameter<int>("maxBX")),
  eta1_(iConfig.getParameter<std::vector<int> >("eta_1")),
  eta2_(iConfig.getParameter<std::vector<int> >("eta_2")),
  eta3_(iConfig.getParameter<std::vector<int> >("eta_3")),
  coarseEta1_(iConfig.getParameter<std::vector<int> >("coarseEta_1")),
  coarseEta2_(iConfig.getParameter<std::vector<int> >("coarseEta_2")),
  coarseEta3_(iConfig.getParameter<std::vector<int> >("coarseEta_3")),
  coarseEta4_(iConfig.getParameter<std::vector<int> >("coarseEta_4")),
  verbose_(iConfig.getParameter<int>("verbose")),
  phiLSB_(iConfig.getParameter<double>("phiLSB")),
  bendingScale_(iConfig.getParameter<std::vector<double> >("bendingScale"))
{

} 



L1TMuCorrelatorDTStubProcessor::~L1TMuCorrelatorDTStubProcessor() {}



bool L1TMuCorrelatorDTStubProcessor::isGoodPhiStub(const L1MuDTChambPhDigi * stub) {
  if (stub->code()<minPhiQuality_)
    return false;
  return true;
}




L1MuCorrelatorHit 
L1TMuCorrelatorDTStubProcessor::buildStub(const L1MuDTChambPhDigi& phiS,const L1MuDTChambThDigi* etaS) {
  
  L1MuCorrelatorHit stub = buildStubNoEta(phiS);


  //Now full eta
  int qeta1=0;
  int qeta2=0;
  int eta1=255;
  int eta2=255; 


  bool hasEta=false;
  for (uint i=0;i<7;++i) {
    if (etaS->position(i)==0)
      continue;
    if (!hasEta) {
      hasEta=true;
      eta1=calculateEta(i,etaS->whNum(),etaS->scNum(),etaS->stNum());
      if (etaS->quality(i)==1)
	qeta1=2;
      else
	qeta1=1;
    }
    else {
      eta2=calculateEta(i,etaS->whNum(),etaS->scNum(),etaS->stNum());
      if (etaS->quality(i)==1)
	qeta2=2;
      else
	qeta2=1;
    }
  }

  int eta=eta1;
  int qeta=qeta1;

  if (qeta2>0) {//both stubs->average
    eta=(eta1+eta2)>>1;
    qeta=0;
    stub.setEta(eta,qeta);
  }
  else if (qeta1>0) {//Good single stub
    eta=eta1;
    qeta=qeta1;
    stub.setEta(eta,qeta);
  }
 
  return stub;

}





L1MuCorrelatorHit
L1TMuCorrelatorDTStubProcessor::buildStubNoEta(const L1MuDTChambPhDigi& phiS) {
  int wheel = phiS.whNum();
  int abswheel = fabs(phiS.whNum());
  int sign  = wheel>0 ? 1: -1;
  int sector = phiS.scNum();
  int station = phiS.stNum();
  double globalPhi = (-180+sector*30)+phiS.phi()*30./2048.;
  if (globalPhi<-180)
    globalPhi+=360;
  if (globalPhi>180)
    globalPhi-=360;
  globalPhi = globalPhi*M_PI/180.;
  int phi = int(globalPhi/phiLSB_);
  double normPhiB = bendingScale_[phiS.stNum()-1]*phiS.phiB()*M_PI/(6*2048.);
  int phiB = int(normPhiB/phiLSB_);
  bool tag = (phiS.Ts2Tag()==1);
  int bx=phiS.bxNum();
  int quality=phiS.code();
  uint tfLayer=0;
  int eta=-255;
  if (station==1) {
    tfLayer=1; 
    eta=coarseEta1_[abswheel];
  }
  else if (station==2) {
    tfLayer=5;
    eta=coarseEta2_[abswheel];

  }
  else if (station==3) {
    tfLayer=8;
    eta=coarseEta3_[abswheel];

  }
  else if (station==4) {
    tfLayer=10;
    eta=coarseEta4_[abswheel];

  }
  else {
    tfLayer=0;
    eta=-255;
  }
  //Now full eta

  eta=eta*sign;
  L1MuCorrelatorHit stub(wheel,sector,station,tfLayer,phi,phiB,tag,
			 bx,quality,eta,0);
  return stub;

}





L1MuCorrelatorHitCollection 
L1TMuCorrelatorDTStubProcessor::makeStubs(const L1MuDTChambPhContainer* phiContainer,const L1MuDTChambThContainer* etaContainer,const L1TMuonBarrelParams& params) {

  L1MuCorrelatorHitCollection out;
  for (int bx=minBX_;bx<=maxBX_;bx++) {
    for (int wheel=-2;wheel<=2;wheel++) {
      for (uint sector=0;sector<12;sector++) {
	for (uint station=1;station<5;station++) {

	  bool hasEta=false;
	  L1MuDTChambThDigi const*  tseta   = etaContainer->chThetaSegm(wheel,station,sector,bx);
	  if (tseta) {
	    hasEta=true;
	  }


	  //	  printf("Wheel=%d LWheel=%d,%d sector=%d station=%d phiMask=%d etaMask=%d\n",wheel,lwheel1,lwheel2,sector,station,phiMask,etaMask);
	  //	  if (abs(wheel)==2 && station==1)
	  //	    continue;
	    

	  L1MuDTChambPhDigi const* high = phiContainer->chPhiSegm1(wheel,station,sector,bx);
	  if (high) {
	    if (high->code()>=minPhiQuality_) {
	      const L1MuDTChambPhDigi& stubPhi = *high;
	      if (hasEta) {
		out.push_back(buildStub(stubPhi,tseta));
	      }
	      else {
		out.push_back(buildStubNoEta(stubPhi));
	      }
	    }
	  }

	  L1MuDTChambPhDigi const* low = phiContainer->chPhiSegm2(wheel,station,sector,bx-1);
	  if (low) {
	    if (low->code()>=minPhiQuality_) {
	      const L1MuDTChambPhDigi& stubPhi = *low;
	      if (hasEta) {
		out.push_back(buildStub(stubPhi,tseta));
	      }
	      else {
		out.push_back(buildStubNoEta(stubPhi));
	      }
	    }
	  }
	}
      }
    }
  }
  return out;
}



int L1TMuCorrelatorDTStubProcessor::calculateEta(uint i, int wheel,uint sector,uint station) {
  int eta=0;
  if (wheel>0) {
	eta=7*wheel+3-i;
      }
  else if (wheel<0) {
	eta=7*wheel+i-3;
  }
  else {
    if (sector==0 || sector==3 ||sector==4 ||sector==7 ||sector==8 ||sector==11)
      eta=i-3;
    else
      eta=3-i;
  }


  if (station==1)
    eta=eta1_[eta+17];
  else if (station==2)
    eta=eta2_[eta+17];
  else 
    eta=eta3_[eta+17];

  return eta;



}



