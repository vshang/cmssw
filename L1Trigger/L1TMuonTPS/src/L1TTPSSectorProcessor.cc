#include "L1Trigger/L1TMuonTPS/interface/L1TTPSSectorProcessor.h"
#include "L1Trigger/L1TMuonTPS/interface/L1TPSLUTs.h"

L1TTPSSectorProcessor::L1TTPSSectorProcessor(const edm::ParameterSet& iConfig): 
  sectorNumber_(iConfig.getParameter<unsigned int>("sectorNumber")),
  barrelSectors_(iConfig.getParameter<std::vector<unsigned int> >("barrelSectors")),
  csc10DegreeChambers_(iConfig.getParameter<std::vector<unsigned int> >("csc10DegreeChambers")),
  csc20DegreeChambers_(iConfig.getParameter<std::vector<unsigned int> >("csc20DegreeChambers")),
  rpcEndcapChambers_(iConfig.getParameter<std::vector<unsigned int> >("rpcEndcapChambers")),
  iRpcChambers_(iConfig.getParameter<std::vector<unsigned int> >("iRpcChambers")),
  trackPhiLowerBound_(iConfig.getParameter<int>("phiLowerBound")),
  trackPhiUpperBound_(iConfig.getParameter<int>("phiUpperBound")),
  phiOffset_(iConfig.getParameter<int>("phiOffset"))
{
  //For a more readable python cfg file
  const edm::ParameterSet& algoSettings = iConfig.getParameter<edm::ParameterSet>("algoSettings");
  verbose_          = algoSettings.getParameter<int>("verbose");
  phiLSB_           = algoSettings.getParameter<double>("phiLSB");
  trackEtaLSB_      = algoSettings.getParameter<double>("etaLSB");
  trackEtaShift_    = algoSettings.getParameter<unsigned int>("etaShift");
  trackCurvLSB_     = algoSettings.getParameter<double>("curvLSB");
  vetoIndex_        = algoSettings.getParameter<std::vector<uint> >("vetoIndex");
  vetoPattern_      = algoSettings.getParameter<std::vector<uint> >("vetoPattern");
}

L1TTPSSectorProcessor::~L1TTPSSectorProcessor() {}

std::vector<l1t::L1TkMuonParticle> L1TTPSSectorProcessor::process(const TrackPtrVector& tracks,const L1MuCorrelatorHitRefVector& stubsAll) {
  //Output collection
  std::vector<l1t::L1TkMuonParticle> out;

  if (verbose_==1) 
    printf("----Filtering the stubs for the sector----\n");

  //First filter the stubs:
  L1MuCorrelatorHitRefVector stubsInSector;
  for (const auto& stub : stubsAll) {
    if (verbose_==1)
      printf("Candidate stub type=%d phiRegion=%d etaRegion=%d depthRegion=%d phi=%d\n",stub->type(),stub->phiRegion(),stub->etaRegion(),stub->depthRegion(),stub->phi());

    if (stub->type()<=1) {//barrel 
      for (const auto& sector : barrelSectors_) {
	if (stub->phiRegion()==int(sector)) {
	  if (verbose_==1)
	    printf("Passed stub type=%d phiRegion=%d etaRegion=%d depthRegion=%d phi=%d\n",stub->type(),stub->phiRegion(),stub->etaRegion(),stub->depthRegion(),stub->phi());
	  stubsInSector.push_back(stub);
	  break;
	}
      }
    }
    else if (stub->type()==2 && (stub->depthRegion()==1 || (stub->depthRegion()>1 && uint(fabs(stub->etaRegion()))==4))) {//CSC 10 degrees chambers in ME2/2,3/2,4/2
      for (const auto& sector : csc10DegreeChambers_) {
	if (stub->phiRegion()==int(sector)) {
	  if (verbose_==1)
	    printf("Passed stub type=%d phiRegion=%d etaRegion=%d depthRegion=%d phi=%d\n",stub->type(),stub->phiRegion(),stub->etaRegion(),stub->depthRegion(),stub->phi());

	  stubsInSector.push_back(stub);
	  break;
	}
      }
   }
    else if (stub->type()==2 && (stub->depthRegion()>1 && uint(fabs(stub->etaRegion()))==5)) {//CSC 
      for (const auto& sector : csc20DegreeChambers_) {
	if (stub->phiRegion()==int(sector)) {
	  if (verbose_==1)
	    printf("Passed stub type=%d phiRegion=%d etaRegion=%d depthRegion=%d phi=%d\n",stub->type(),stub->phiRegion(),stub->etaRegion(),stub->depthRegion(),stub->phi());
	  stubsInSector.push_back(stub);
	  break;
	}
      }
    }
    else if (stub->type()==3 &&uint(fabs(stub->etaRegion()))!=5) {//RPC
      for (const auto& sector : rpcEndcapChambers_) {
	if (stub->phiRegion()==int(sector)) {
	  if (verbose_==1)
	    printf("Passed stub type=%d phiRegion=%d etaRegion=%d depthRegion=%d phi=%d\n",stub->type(),stub->phiRegion(),stub->etaRegion(),stub->depthRegion(),stub->phi());
	  stubsInSector.push_back(stub);
	  break;
	}
      }
    }
    else if (stub->type()==3 &&uint(fabs(stub->etaRegion()))==5) {//iRPC
      for (const auto& sector : iRpcChambers_) {
	if (stub->phiRegion()==int(sector)) {
	  if (verbose_==1)
	    printf("Passed stub type=%d phiRegion=%d etaRegion=%d depthRegion=%d phi=%d\n",stub->type(),stub->phiRegion(),stub->etaRegion(),stub->depthRegion(),stub->phi());
	  stubsInSector.push_back(stub);
	  break;
	}
      }
    }


  }
  if (verbose_==1) {
    printf("---processing sector=%d - stubs to be processed=%d---\n",sectorNumber_,int(stubsInSector.size()));
  }

  //For faster emulator if there are no stubs at all return
  if(stubsInSector.size()==0)
    return out;

  //Now loop on the tracks
  for (const auto& track: tracks) {
    l1t::L1TkMuonParticle::LorentzVector vec(track->getMomentum().x(),
					track->getMomentum().y(),
					track->getMomentum().z(),
					track->getMomentum().mag());
    l1t::L1TkMuonParticle muon (vec,track);   

    if (muon.pt()<3.0)
      continue;

    //Set muon charge
    int charge=1;
    if (track->getRInv()<0)
      charge=-1;
    muon.setCharge(charge);
    //If not in nonant kill
    int globalPhi = muon.phi()/phiLSB_;
    int deltalow = deltaPhi(globalPhi,trackPhiLowerBound_);
    int deltahigh = deltaPhi(globalPhi,trackPhiUpperBound_);
    if(verbose_==1)
      printf("---> Testing track for sector with phi = %d in [%d %d]  and deltaPhi = %dm%d with eta=%f and pt=%f \n",int(muon.phi()/phiLSB_),trackPhiLowerBound_,trackPhiUpperBound_,deltalow,deltahigh,muon.eta(),muon.pt());
    
    if (globalPhi< trackPhiLowerBound_|| globalPhi>=trackPhiUpperBound_)
      continue;
    if(verbose_==1)
      printf("Passed\n");
    if(processTrack(muon,stubsInSector))
      out.push_back(muon);
  }

  std::vector<l1t::L1TkMuonParticle> cleaned = clean(out);
  return cleaned;
}

int L1TTPSSectorProcessor::deltaPhi(int phi1, int phi2) {
  //add offset code here 
  int delta = phi1-phi2;
  int pi = int(3.14159/phiLSB_);
  if (delta>pi)
    return delta-2*pi;
  if (delta<-pi)
    return delta+2*pi;
  return delta;
}


int L1TTPSSectorProcessor::stubPhi(const L1MuCorrelatorHitRef& stub ) {
  //add offset code here 
  return deltaPhi(stub->phi(),phiOffset_);
}

int L1TTPSSectorProcessor::trackPhi(const l1t::L1TkMuonParticle& track) {
  //add offset code here 
  return deltaPhi(int(track.phi()/phiLSB_),phiOffset_);
}

int L1TTPSSectorProcessor::trackEta(const l1t::L1TkMuonParticle& track) {
  //add offset code here 
  return int(track.eta()/trackEtaLSB_);
}

int L1TTPSSectorProcessor::trackCurv(const l1t::L1TkMuonParticle& track) {
  //add offset code here 

  return int(track.charge()/track.pt()/trackCurvLSB_);
}



L1TTPSSectorProcessor::PropagationInfo L1TTPSSectorProcessor::propagate(const l1t::L1TkMuonParticle& track,uint propIndex) {
  L1TTPSSectorProcessor::PropagationInfo out;
  //first calculate the index
  int eta = trackEta(track);
  uint absEta = uint(fabs(eta))>>trackEtaShift_;
  out.propagatedEta=eta>>trackEtaShift_;
 
  uint etaIndex = L1TTPS::etaIndex[absEta];

  if (verbose_==1) {
    printf("Calculating eta index , eta=%d |eta|=%d index=%d\n",eta,absEta,etaIndex);
  }

  double slope,resa,resb;
  double resa_phib = 0.0;
  double resb_phib = 0.0;
  uint hasPhiB;
  bool valid;

  switch(propIndex) {
  case 0:
    slope   = L1TTPS::slope_0[absEta];
    resa    = L1TTPS::resa_0[absEta];
    resb    = L1TTPS::resb_0[absEta];
    hasPhiB = 0;
    valid   = L1TTPS::valid_0[absEta];
    break;
  case 1:
    slope=L1TTPS::slope_1[absEta];
    resa =L1TTPS::resa_1[absEta];
    resb =L1TTPS::resb_1[absEta];
    resa_phib =L1TTPS::resa_phiB_1[absEta];
    resb_phib =L1TTPS::resb_phiB_1[absEta];
    hasPhiB = L1TTPS::hasPhiB_1[absEta];
    valid = L1TTPS::valid_1[absEta];
    break;
  case 2:
    slope   = L1TTPS::slope_2[absEta];
    resa    = L1TTPS::resa_2[absEta];
    resb    = L1TTPS::resb_2[absEta];
    hasPhiB = 0;
    valid   = L1TTPS::valid_2[absEta];
    break;
  case 3:
    slope   = L1TTPS::slope_3[absEta];
    resa    = L1TTPS::resa_3[absEta];
    resb    = L1TTPS::resb_3[absEta];
    hasPhiB = 0;
    valid   = L1TTPS::valid_3[absEta];
    break;
  case 4:
    slope   = L1TTPS::slope_4[absEta];
    resa    = L1TTPS::resa_4[absEta];
    resb    = L1TTPS::resb_4[absEta];
    hasPhiB = 0;
    valid   = L1TTPS::valid_4[absEta];
    break;
  case 5:
    slope=L1TTPS::slope_5[absEta];
    resa =L1TTPS::resa_5[absEta];
    resb =L1TTPS::resb_5[absEta];
    resa_phib =L1TTPS::resa_phiB_5[absEta];
    resb_phib =L1TTPS::resb_phiB_5[absEta];
    hasPhiB = L1TTPS::hasPhiB_5[absEta];
    valid = L1TTPS::valid_5[absEta];
    break;
  case 6:
    slope   = L1TTPS::slope_6[absEta];
    resa    = L1TTPS::resa_6[absEta];
    resb    = L1TTPS::resb_6[absEta];
    hasPhiB = 0;
    valid   = L1TTPS::valid_6[absEta];
    break;
  case 7:
    slope   = L1TTPS::slope_7[absEta];
    resa    = L1TTPS::resa_7[absEta];
    resb    = L1TTPS::resb_7[absEta];
    hasPhiB = 0;
    valid   = L1TTPS::valid_7[absEta];
    break;
  case 8:
    slope=L1TTPS::slope_8[absEta];
    resa =L1TTPS::resa_8[absEta];
    resb =L1TTPS::resb_8[absEta];
    resa_phib =L1TTPS::resa_phiB_8[absEta];
    resb_phib =L1TTPS::resb_phiB_8[absEta];
    hasPhiB = L1TTPS::hasPhiB_8[absEta];
    valid = L1TTPS::valid_8[absEta];
    break;
  case 9:
    slope   = L1TTPS::slope_9[absEta];
    resa    = L1TTPS::resa_9[absEta];
    resb    = L1TTPS::resb_9[absEta];
    hasPhiB = 0;
    valid   = L1TTPS::valid_9[absEta];
    break;
  case 10:
    slope=L1TTPS::slope_10[absEta];
    resa =L1TTPS::resa_10[absEta];
    resb =L1TTPS::resb_10[absEta];
    resa_phib =L1TTPS::resa_phiB_10[absEta];
    resb_phib =L1TTPS::resb_phiB_10[absEta];
    hasPhiB = L1TTPS::hasPhiB_10[absEta];
    valid = L1TTPS::valid_10[absEta];
    break;
  default:
    slope=0;
    resa =0;
    resb =0;
    hasPhiB = 0;
    valid = 0;
  }

  
  int phi = trackPhi(track);
  int curv = trackCurv(track);
  double prop = slope*curv;

  out.propagatedAngle = int(phi+prop);
  out.propagatedBend = hasPhiB ?  int(prop) : -1;
  //correct for boundary. Will be done automatically in firmware 
  int pi = M_PI/phiLSB_;
  if (out.propagatedAngle>pi)
    out.propagatedAngle = out.propagatedAngle-2*pi;
  if (out.propagatedAngle<-pi)
    out.propagatedAngle = out.propagatedAngle+2*pi;

  out.propagatorIndex = propIndex;
  out.propagatedSigmaAngle = uint(fabs(resa*curv)+resb);
  out.propagatedSigmaBend = hasPhiB ? uint(fabs(resa_phib*curv)+resb_phib) : 0.0;
  if (out.propagatedSigmaAngle>1023)
    out.propagatedSigmaAngle=1023;
  if (out.propagatedSigmaBend>1023)
    out.propagatedSigmaBend=1023;
  out.etaIndex=etaIndex;
  out.valid = valid;
  out.hasPhiB = hasPhiB;

  if (verbose_==1) {
    printf("Propagating index=%d eta=%d phi=%d curvature=%d angle=%d sigma=%d  bend=%d bend_sigma=%d eta index=%d valid=%d\n",propIndex,trackEta(track),phi,curv,out.propagatedAngle,out.propagatedSigmaAngle,out.propagatedBend,out.propagatedSigmaBend,out.etaIndex,out.valid);
  }

  return out;
}


uint L1TTPSSectorProcessor::match(l1t::L1TkMuonParticle& muon,const L1TTPSSectorProcessor::PropagationInfo& prop,const L1MuCorrelatorHitRefVector& stubs,uint& pattern) {
  if (!prop.valid)
    return 0; 
  //retrieve eta information for matching

  uint out=0;

  int bestStub=-1;
  uint maxDPhi=10000;
  
  L1MuCorrelatorHitRefVector stubsInLayer;
  for (uint i = 0;i<stubs.size();++i) {
    const L1MuCorrelatorHitRef& stub = stubs[i];
    if (stub->tfLayer() != prop.propagatorIndex)
      continue;
    if(verbose_==1)
      printf("Looking at stub at layer=%d etaRegion=%d depthRegion=%d type=%d\n",
	     stubs[i]->tfLayer(),
	     stubs[i]->etaRegion(),
	     stubs[i]->depthRegion(),
	     stubs[i]->type());
    
    //First Check eta
    uint etaCut = stub->etaQuality()==0 ? 100 :16;
    if (fabs(stub->eta()-prop.propagatedEta)>etaCut)
      continue;
    
    //Then check phi
    uint dPhi=fabs(deltaPhi(prop.propagatedAngle,stubPhi(stubs[i])));
    if (verbose_==1)
      printf("Found stub with phi =  %d phiB=%d sigma=%d  sigmaB=%d and deltaPhi=%d\n ",stubPhi(stubs[i]),stubs[i]->phiB(),prop.propagatedSigmaAngle,prop.propagatedSigmaBend,dPhi);
    if (dPhi<maxDPhi){
      maxDPhi=dPhi;
      bestStub=i;
    } 
  }

  if (bestStub!=-1) {
    uint deltaAngle = fabs(stubs[bestStub]->phi()-prop.propagatedAngle);
    uint deltaBend = fabs(stubs[bestStub]->phiB()-prop.propagatedBend);

    bool passAngle = deltaAngle<prop.propagatedSigmaAngle;
    bool passBend  = passAngle ? deltaBend<prop.propagatedSigmaBend : 0;

    

    if (verbose_==1)
      printf("Evaluating Best Stub dPhi=%d dBend=%d passed=%d %d\n",deltaAngle,deltaBend,passAngle,passBend);
    
    uint muonQuality = muon.quality();
    
    if (passAngle) {
      muonQuality+=8 - (deltaAngle<<4);
      muon.addStub(stubs[bestStub]);
      out+=1;
    }

    if (prop.hasPhiB) {
      pattern=(pattern<<2) | (passBend<<1) | passAngle;
      if (passBend) {
	muonQuality+=8 - (deltaBend<<4);
	out+=1;
      }
    }
    else {
      pattern=(pattern<<1) |passAngle;
    }      
  }

  return out;
}


  



bool L1TTPSSectorProcessor::processTrack(l1t::L1TkMuonParticle& muon,const L1MuCorrelatorHitRefVector& stubs) {
  //Propagate 12 times maximum and match
  muon.setQuality(288);
  uint pattern=0;
  uint nstubs=0;
  for (uint i=0;i<12;++i) {
    //propagate
    L1TTPSSectorProcessor::PropagationInfo prop = propagate(muon,i);
    nstubs+=match(muon,prop,stubs,pattern);
  }
  muon.setPattern(pattern);
  

  //Stubs requirements
  if (muon.getMatchedStubs().size()<2)
    return false;


  int eta = trackEta(muon);
  uint absEta = uint(fabs(eta))>>trackEtaShift_;
  uint etaIndex = L1TTPS::etaIndex[absEta];


  for (uint i =0 ;i<vetoIndex_.size();++i)
    if ((etaIndex==vetoIndex_[i]) && pattern==vetoPattern_[i])
      return false;

  
  return true;
}





std::vector<l1t::L1TkMuonParticle> L1TTPSSectorProcessor::clean(const std::vector<l1t::L1TkMuonParticle>& muons) {
  //    return muons;
    if (verbose_==1)
      printf("CROSS CLEANER IN SECTOR\n");
  std::vector<l1t::L1TkMuonParticle> out;
  if (muons.size()<=1)
    return muons;

  for (uint i=0;i<muons.size();++i) {
    if (verbose_==1)
      printf("->Muon with pt,eta,phi=%f %f %f\n",muons[i].pt(),muons[i].eta(),muons[i].phi());
    bool keep=true;
    const L1MuCorrelatorHitRefVector& stubs1 = muons[i].getMatchedStubs();
    if (verbose_==1)
      for (const auto& stub : stubs1)
	printf("stub %d %d %d %d %d\n",stub->etaRegion(),stub->phiRegion(),stub->depthRegion(),stub->id(),stub->phi());
    for (uint j=0;j<muons.size();++j) {
      if (i==j)
	continue;
      const L1MuCorrelatorHitRefVector& stubs2 = muons[j].getMatchedStubs();
      int overlap=0;
      for (const auto& stub1 : stubs1) {
	for (const auto& stub2 : stubs2) {
	  if ((*stub1)==(*stub2)) {
 	    overlap+=1;
	  }
	}
      }

      if ((muons[0].getMatchedStubs().size()-overlap)<2 && (muons[i].quality()<muons[j].quality())) {
	keep=false;
	break;
      }
    }
    if (keep)
      out.push_back(muons[i]);
  }


  if (verbose_==1) {
      printf("CLEANED MUONS\n");
      for (const auto& mu :out)
      printf("Muon with pt,eta,phi=%f %f %f\n",mu.pt(),mu.eta(),mu.phi());

  }


  return out;
}




