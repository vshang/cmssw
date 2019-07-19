#ifndef Phase2L1Trigger_DTTrigger_RPCIntegrator_cc
#define Phase2L1Trigger_DTTrigger_RPCIntegrator_cc

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/DTSuperLayerId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/DTDigi/interface/DTDigiCollection.h"
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTPhDigi.h"

#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"

#include "L1Trigger/DTPhase2Trigger/interface/analtypedefs.h"
//#include "L1Trigger/DTPhase2Trigger/interface/constants.h"


class RPCIntegrator {
    public:
        RPCIntegrator(const edm::ParameterSet& pset);
        ~RPCIntegrator();

        void initialise(const edm::EventSetup& iEventSetup);
        GlobalPoint getRPCGlobalPosition(RPCDetId rpcId, const RPCRecHit& rpcIt) const;
        void translateRPC(edm::Handle<RPCRecHitCollection> rpcRecHits);
        L1Phase2MuDTPhDigi* matchDTwithRPC(metaPrimitive* dt_metaprimitive);
        void confirmDT(std::vector<metaPrimitive> & dt_metaprimitives, double shift_back);

        std::vector<L1Phase2MuDTPhDigi> rpcRecHits_translated;

    private:
        //RPCRecHitCollection m_rpcRecHits;
        Bool_t m_debug;
        int m_max_quality_to_overwrite_t0;
        edm::ESHandle<RPCGeometry> m_rpcGeo;
        int m_dt_phi_granularity = 81920;
};
#endif
