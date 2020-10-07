#include "L1Trigger/L1TMuonBayes/interface/MuonStubMakerBase.h"

#include "L1Trigger/L1TMuonBayes/interface/AngleConverterBase.h"
#include "L1Trigger/L1TMuonBayes/interface/Omtf/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonBayes/interface/Omtf/OMTFinput.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "L1Trigger/CSCCommonTrigger/interface/CSCConstants.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

///////////////////////////////////////
///////////////////////////////////////
MuonStubMakerBase::MuonStubMakerBase(): rpcClusterization(3, 2) {

}

///////////////////////////////////////
///////////////////////////////////////
void MuonStubMakerBase::initialize(const edm::ParameterSet& edmCfg, const edm::EventSetup& es, const ProcConfigurationBase* procConf, MuStubsInputTokens& muStubsInputTokens) {
  config = procConf;

  this->muStubsInputTokens = muStubsInputTokens;

  dropDTPrimitives = edmCfg.getParameter<bool>("dropDTPrimitives");
  dropRPCPrimitives = edmCfg.getParameter<bool>("dropRPCPrimitives");
  dropCSCPrimitives = edmCfg.getParameter<bool>("dropCSCPrimitives");

  if(edmCfg.exists("minDtPhQuality") ) {
    minDtPhQuality = edmCfg.getParameter<int>("minDtPhQuality");
  }

}
///////////////////////////////////////
///////////////////////////////////////
MuonStubMakerBase::~MuonStubMakerBase(){ }
///////////////////////////////////////
///////////////////////////////////////

void MuonStubMakerBase::loadAndFilterDigis(const edm::Event& event) {
  // Filter digis by dropping digis from selected (by cfg.py) subsystems
  if(!dropDTPrimitives){
    event.getByToken(muStubsInputTokens.inputTokenDTPh, dtPhDigis);
    event.getByToken(muStubsInputTokens.inputTokenDTTh, dtThDigis);
  }
  if(!dropRPCPrimitives) event.getByToken(muStubsInputTokens.inputTokenRPC, rpcDigis);
  if(!dropCSCPrimitives) event.getByToken(muStubsInputTokens.inputTokenCSC, cscDigis);
}


const void MuonStubMakerBase::buildInputForProcessor(MuonStubPtrs2D& muonStubsInLayers, unsigned int iProcessor,
    l1t::tftype type,
    int bxFrom, int bxTo) {
  LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" iProcessor "<<iProcessor<<std::endl;
  processDT( muonStubsInLayers, dtPhDigis.product(), dtThDigis.product(), iProcessor, type, false, bxFrom, bxTo);
  processCSC(muonStubsInLayers, cscDigis.product(), iProcessor, type, bxFrom, bxTo);
  processRPC(muonStubsInLayers, rpcDigis.product(), iProcessor, type, bxFrom, bxTo);
  //cout<<result<<endl;
}


//iProcessor counted from 0
/*int MuonStubMakerBase::getProcessorPhiZero(unsigned int iProcessor) {

  unsigned int nPhiBins = config->nPhiBins();

  int phiZero =  nPhiBins/6.*(iProcessor) + nPhiBins/24;  // "0" is 15degree moved cyclicaly to each processor, note [0,2pi]
  return config->foldPhi(phiZero);

}*/

void MuonStubMakerBase::processDT(MuonStubPtrs2D& muonStubsInLayers, const L1MuDTChambPhContainer *dtPhDigis,
    const L1MuDTChambThContainer *dtThDigis,
    unsigned int iProcessor, l1t::tftype procTyp, bool mergePhiAndTheta, int bxFrom, int bxTo)
{
  if(!dtPhDigis)
    return;

  for (const auto& digiIt: *dtPhDigis->getContainer()) {
    DTChamberId detid(digiIt.whNum(),digiIt.stNum(),digiIt.scNum()+1);

    ///Check it the data fits into given processor input range
    if(!acceptDigi(detid.rawId(), iProcessor, procTyp))
      continue;

    if (digiIt.bxNum() >= bxFrom && digiIt.bxNum() <= bxTo )
      addDTphiDigi(muonStubsInLayers, digiIt, dtThDigis, iProcessor, procTyp);
  }

  if(!mergePhiAndTheta) {
    for(auto& thetaDigi: (*(dtThDigis->getContainer()) ) ) {
      if(thetaDigi.bxNum() >= bxFrom && thetaDigi.bxNum() <= bxTo) {
        addDTetaStubs(muonStubsInLayers, thetaDigi, iProcessor, procTyp);
      }
    }
  }
  //std::cout<<__FUNCTION__<<":"<<__LINE__<<" iProcessor "<<iProcessor<<std::endl;
  //angleConverter->AngleConverterBase::getGlobalEta(dtThDigis, 0, 0);
}

////////////////////////////////////////////
////////////////////////////////////////////
void MuonStubMakerBase::processCSC(MuonStubPtrs2D& muonStubsInLayers, const CSCCorrelatedLCTDigiCollection *cscDigis,
    unsigned int iProcessor,
    l1t::tftype procTyp, int bxFrom, int bxTo) {

  int lctCentralBx = CSCConstants::LCT_CENTRAL_BX;
  lctCentralBx = 6; //TODO this was changed in CMSSW 10(?) to 8 if the data were generated with the previous CMSSW then you have to use 6

  if(!cscDigis)
    return;

  auto chamber = cscDigis->begin();
  auto chend  = cscDigis->end();
  for( ; chamber != chend; ++chamber ) {
    unsigned int rawid = (*chamber).first;
    ///Check it the data fits into given processor input range
    if(!acceptDigi(rawid, iProcessor, procTyp))
      continue;

    auto digi = (*chamber).second.first;
    auto dend = (*chamber).second.second;
    for( ; digi != dend; ++digi ) {
      ///Check if LCT trigger primitive has the right BX.
      int digiBx = digi->getBX()- lctCentralBx;

      if (digiBx >= bxFrom && digiBx <= bxTo )
        addCSCstubs(muonStubsInLayers, rawid, *digi, iProcessor, procTyp);
    }
  }
}

////////////////////////////////////////////
////////////////////////////////////////////
void MuonStubMakerBase::processRPC(MuonStubPtrs2D& muonStubsInLayers, const RPCDigiCollection *rpcDigis,
    unsigned int iProcessor,
    l1t::tftype procTyp, int bxFrom, int bxTo)
{
  if(!rpcDigis) return;
  //  std::cout <<" RPC HITS, processor : " << iProcessor << std::endl;

  const RPCDigiCollection & rpcDigiCollection = *rpcDigis;
  for (auto rollDigis : rpcDigiCollection) {
    RPCDetId roll = rollDigis.first;    
    unsigned int rawid = roll.rawId();

//    if(roll.region() != 0  &&  abs(roll.station()) >= 3 && roll.ring() == 1 ) {
//      //iRPC
//      for (auto pDigi=rollDigis.second.first; pDigi != rollDigis.second.second; pDigi++) {
//        LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" roll "<<roll
//            <<" strip "<<pDigi->strip()
//            <<" hasX "<<pDigi->hasX()<<" coordinateX "<<pDigi->coordinateX()<<" hasY "<<pDigi->hasY()<<" coordinateY "<<pDigi->coordinateY()
//            <<" bx "<<pDigi->bx()<<" time "<<pDigi->time()<<" irpc"<<std::endl;
//      }
//      //continue;
//    }

    if(!acceptDigi(rawid, iProcessor, procTyp))
      continue;

    ///To find the clusters we have to copy the digis in chamber to sort them (not optimal).
    //  for (auto tdigi = rollDigis.second.first; tdigi != rollDigis.second.second; tdigi++) { std::cout << "RPC DIGIS: " << roll.rawId()<< " "<<roll<<" digi: " << tdigi->strip() <<" bx: " << tdigi->bx() << std::endl; }
    std::vector<RPCDigi> digisCopy;
    //  std::copy_if(rollDigis.second.first, rollDigis.second.second, std::back_inserter(digisCopy), [](const RPCDigi & aDigi){return (aDigi.bx()==0);} );
    for (auto pDigi=rollDigis.second.first; pDigi != rollDigis.second.second; pDigi++) {
      if(pDigi->bx() >= bxFrom && pDigi->bx() <= bxTo ) {
        digisCopy.push_back( *pDigi);
      }
    }

    std::vector<RpcCluster> clusters = rpcClusterization.getClusters(roll, digisCopy);

    for (auto & cluster: clusters) {
      //LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<155<<" roll "<<roll<<" cluster: firstStrip "<<cluster.firstStrip<<" lastStrip "<<cluster.lastStrip<<" halfStrip "<<cluster.halfStrip()<<std::endl;
      addRPCstub(muonStubsInLayers, roll, cluster, iProcessor, procTyp);
    }
  }
}

//should it use GEMPadDigiClusterCollection???
void MuonStubMakerBase::processGEM(MuonStubPtrs2D& muonStubsInLayers, const GEMPadDigiCollection* gemDigis,
    unsigned int iProcessor,
    l1t::tftype procType, int bxFrom, int bxTo) {

  if(!gemDigis) return;

  for(auto chamber = gemDigis->begin() ; chamber != gemDigis->end(); ++chamber ) {
    for(auto digi = (*chamber).second.first ; digi != (*chamber).second.second; ++digi ) {
      //muon_primitives.emplace_back((*chamber).first, *digi);

      GEMDetId geDetId = (*chamber).first;
      unsigned int rawid = geDetId.rawId();

      if(!acceptDigi(rawid, iProcessor, procType))
        continue;

      //TODO implement
    }
  }

}
////////////////////////////////////////////
////////////////////////////////////////////
/*OMTFinput MuonStubMakerBase::buildInputForProcessor(const L1MuDTChambPhContainer *dtPhDigis,
    const L1MuDTChambThContainer *dtThDigis,
    const CSCCorrelatedLCTDigiCollection *cscDigis,
    const RPCDigiCollection *rpcDigis,
    unsigned int iProcessor,
    l1t::tftype procType,
    int bxFrom, int bxTo) {
  OMTFinput result(config);
  processDT(result, dtPhDigis, dtThDigis, iProcessor, type, bx);
  processCSC(result, cscDigis, iProcessor, type, bx);
  processRPC(result, rpcDigis, iProcessor, type, bx);
  //cout<<result<<endl;
  return result;
}*/
////////////////////////////////////////////
////////////////////////////////////////////
