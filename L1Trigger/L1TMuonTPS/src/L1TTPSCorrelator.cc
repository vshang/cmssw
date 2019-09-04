#include "L1Trigger/L1TMuonTPS/interface/L1TTPSCorrelator.h"
#include "L1Trigger/L1TMuonTPS/interface/L1TPSLUTs.h"

L1TTPSCorrelator::L1TTPSCorrelator(const edm::ParameterSet& iConfig){

  verbose_          = iConfig.getParameter<int>("verbose");
  std::vector<edm::ParameterSet> info = iConfig.getParameter<std::vector<edm::ParameterSet> > ("sectors");
  for (const auto& set : info)
    processor_.push_back(L1TTPSSectorProcessor(set));
}

L1TTPSCorrelator::~L1TTPSCorrelator() {}

std::vector<l1t::L1TkMuonParticle> L1TTPSCorrelator::process(const TrackPtrVector& tracks,const L1MuCorrelatorHitRefVector& stubsAll) {
  std::map<uint,std::vector<l1t::L1TkMuonParticle> > sectorData;
  if (verbose_==1)
    printf("--------------- NEW EVENT ---------------\n");

  for (auto& proc : processor_) {
    std::vector<l1t::L1TkMuonParticle> data = proc.process(tracks,stubsAll);
    sectorData[proc.sector()] = data;
  }

  std::vector<l1t::L1TkMuonParticle> out;

  //Now clean the sectors
  for (const auto& proc : processor_) {
    uint sector = proc.sector();
    std::vector<l1t::L1TkMuonParticle> muons = sectorData[sector];
    std::vector<l1t::L1TkMuonParticle> muonsP;
    std::vector<l1t::L1TkMuonParticle> muonsN;
    if (sector==0) {
      if (sectorData.find(9)!=sectorData.end())
	muonsP = sectorData[8];
      if (sectorData.find(1)!=sectorData.end())
	muonsN = sectorData[1];
    }
    else if (sector==8) {
      if (sectorData.find(0)!=sectorData.end())
	muonsP = sectorData[0];
      if (sectorData.find(8)!=sectorData.end())
	muonsN = sectorData[7];
    }
    else {
      if (sectorData.find(sector+1)!=sectorData.end())
	muonsP = sectorData[sector+1];
      if (sectorData.find(sector-1)!=sectorData.end())
	muonsN = sectorData[sector-1];
    }

    std::vector<l1t::L1TkMuonParticle> tmp = clean(muons,muonsP,muonsN);
    std::copy(tmp.begin(), tmp.end(), std::back_inserter(out));
  }
  return out;
}


std::vector<l1t::L1TkMuonParticle> L1TTPSCorrelator::clean(const std::vector<l1t::L1TkMuonParticle>& central,const std::vector<l1t::L1TkMuonParticle>& before, const std::vector<l1t::L1TkMuonParticle>& after) {
  //  return central;

  std::vector<l1t::L1TkMuonParticle> tmp = before;
  //merge the other two
  std::copy (after.begin(), after.end(), std::back_inserter(tmp));

  std::vector<l1t::L1TkMuonParticle> out;
  for (uint i=0;i<central.size();++i) {
    bool keep=true;
    const L1MuCorrelatorHitRefVector& stubs1 = central[i].getMatchedStubs();
    for (uint j=0;j<tmp.size();++j) {
      const L1MuCorrelatorHitRefVector& stubs2 = tmp[j].getMatchedStubs();
      //check if stubs match
      bool overlap=false;
      for (const auto& stub1 : stubs1) {
	for (const auto& stub2 : stubs2) {
	  if (stub1==stub2) {
 	    overlap=true;
	    break;
	  }
	  
	}
      }
      if (overlap && (central[i].quality()<tmp[j].quality())) {
	keep=false;
	break;
      }
    }
    if (keep) {
      out.push_back(central[i]);
    }
  }
  return out;
}

