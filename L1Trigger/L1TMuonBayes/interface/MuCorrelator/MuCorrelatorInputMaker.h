/*
 * MuonCorrelatorInputMaker.h
 *
 *  Created on: Jan 30, 2019
 *      Author: Karol Bunkowski kbunkow@cern.ch
 */

#ifndef MUCORRELATOR_MUONCORRELATORINPUTMAKER_H_
#define MUCORRELATOR_MUONCORRELATORINPUTMAKER_H_


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "L1Trigger/L1TMuonBayes/interface/MuonStubMakerBase.h"
#include "L1Trigger/L1TMuonBayes/interface/MuonStubsInput.h"
#include "L1Trigger/L1TMuonBayes/interface/AngleConverterBase.h"
#include "L1Trigger/L1TMuonBayes/interface/MuCorrelator/MuCorrelatorConfig.h"



namespace edm {
  class EventSetup;
}

class MuCorrelatorInputMaker : public MuonStubMakerBase {
public:
  MuCorrelatorInputMaker(const edm::ParameterSet& edmCfg, const edm::EventSetup& es, MuCorrelatorConfigPtr config, MuStubsInputTokens& muStubsInputTokens);
  virtual ~MuCorrelatorInputMaker();

private:

  virtual bool acceptDigi(uint32_t rawId, unsigned int iProcessor, l1t::tftype procType) {
    return true;
  }

  //the phi and eta digis are merged (even thought it is artificial)
  virtual void addDTphiDigi(MuonStubPtrs2D& muonStubsInLayers, const L1MuDTChambPhDigi& digi,
     const L1MuDTChambThContainer *dtThDigis,
     unsigned int iProcessor,
     l1t::tftype procTyp);

  virtual void addDTetaStubs(MuonStubPtrs2D& muonStubsInLayers, const L1MuDTChambThDigi& thetaDigi,
       unsigned int iProcessor, l1t::tftype procTyp);

  virtual void addCSCstubs(MuonStubPtrs2D& muonStubsInLayers,  unsigned int rawid, const CSCCorrelatedLCTDigi& digi,
     unsigned int iProcessor, l1t::tftype procTyp);

  virtual void addRPCstub(MuonStubPtrs2D& muonStubsInLayers, const RPCDetId& roll, const RpcCluster& cluster,
     unsigned int iProcessor, l1t::tftype procTyp);

  virtual void addStub(MuonStubPtrs2D& muonStubsInLayers, unsigned int iLayer, MuonStub& stub);

  //logic layer numbers
  uint32_t getLayerNumber(const DTChamberId& detid, bool eta = false) const;
  uint32_t getLayerNumber(const CSCDetId& detid, bool eta = false) const;
  uint32_t getLayerNumber(const RPCDetId& detid) const;

  MuCorrelatorConfigPtr config;

  AngleConverterBase angleConverter;
};

#endif /* INTERFACE_MUCORRELATOR_MUONCORRELATORINPUTMAKER_H_ */
