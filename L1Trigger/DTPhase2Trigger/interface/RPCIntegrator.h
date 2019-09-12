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
        void finish();
        GlobalPoint getRPCGlobalPosition(RPCDetId rpcId, const RPCRecHit& rpcIt) const;
        void translateRPC(edm::Handle<RPCRecHitCollection> rpcRecHits);
        L1Phase2MuDTPhDigi* matchDTwithRPC(metaPrimitive* dt_metaprimitive, double shift_back);
        void confirmDT(std::vector<metaPrimitive> & dt_metaprimitives, double shift_back);
        void makeRPConlySegments();
        int getPhiBending(L1Phase2MuDTPhDigi* rpc_hit_1, L1Phase2MuDTPhDigi* rpc_hit_2);
        int getLocalX(int phi, int layer, int station);
        bool hasPosRF_rpc(int wh, int sec);
        void removeRPCHitsMatched();

        std::vector<L1Phase2MuDTPhDigi> rpcRecHits_translated;

    private:
        //RPCRecHitCollection m_rpcRecHits;
        Bool_t m_debug;
        int m_max_quality_to_overwrite_t0;
        int m_bx_window;
        bool m_storeAllRPCHits;
        edm::ESHandle<RPCGeometry> m_rpcGeo;
        int m_dt_phi_granularity = 81920;
        //R[stat][layer] - radius of rpc station/layer from center of CMS
        double R[2][2] = {{410.0, 444.8}, {492.7, 527.3}};
        //float m_radius_rb1_layer1 = 410.0, m_radius_rb1_layer2 = 444.8, m_radius_rb2_layer1 = 492.7, m_radius_rb2_layer2 = 527.3;
};
#endif
