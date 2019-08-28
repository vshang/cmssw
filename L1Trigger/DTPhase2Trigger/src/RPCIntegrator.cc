#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/ESInputTag.h"

#include "L1Trigger/DTPhase2Trigger/interface/RPCIntegrator.h"

#include <Geometry/Records/interface/MuonGeometryRecord.h>

RPCIntegrator::RPCIntegrator(const edm::ParameterSet& pset){
    m_debug = pset.getUntrackedParameter<Bool_t>("debug");
    if (m_debug) std::cout <<"RPCIntegrator constructor" << std::endl;
    m_max_quality_to_overwrite_t0 = pset.getUntrackedParameter<int>("max_quality_to_overwrite_t0");
    m_bx_window = pset.getUntrackedParameter<int>("bx_window");
    m_storeAllRPCHits = pset.getUntrackedParameter<bool>("storeAllRPCHits");
}

RPCIntegrator::~RPCIntegrator() {
    if (m_debug) std::cout <<"RPCIntegrator destructor" << std::endl;
}

void RPCIntegrator::initialise(const edm::EventSetup& iEventSetup) {
    if (m_debug) std::cout << "RPCIntegrator initialisation" << std::endl;
    
    if (m_debug) std::cout << "Getting RPC geometry" << std::endl;
    iEventSetup.get<MuonGeometryRecord>().get(m_rpcGeo);
}

GlobalPoint RPCIntegrator::getRPCGlobalPosition(RPCDetId rpcId, const RPCRecHit& rpcIt) const{

  RPCDetId rpcid = RPCDetId(rpcId);
  const LocalPoint& rpc_lp = rpcIt.localPosition();
  const GlobalPoint& rpc_gp = m_rpcGeo->idToDet(rpcid)->surface().toGlobal(rpc_lp);

  return rpc_gp;

}

void RPCIntegrator::translateRPC(edm::Handle<RPCRecHitCollection> rpcRecHits) {
    rpcRecHits_translated.clear();
    for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {
        // Retrieve RPC info and translate it to DT convention if needed
        int rpc_bx = rpcIt->BunchX();
        int rpc_time = int(rpcIt->time());
        RPCDetId rpcDetId = (RPCDetId)(*rpcIt).rpcId();
        if (m_debug) std::cout << "Getting RPC info from : " << rpcDetId << std::endl;
        int rpc_region = rpcDetId.region();
        if(rpc_region != 0 ) continue; // Region = 0 Barrel
        int rpc_wheel = rpcDetId.ring(); // In barrel, wheel is accessed via ring() method ([-2,+2])
        int rpc_dt_sector = rpcDetId.sector()-1; // DT sector:[0,11] while RPC sector:[1,12]
        int rpc_station = rpcDetId.station();
        int rpc_layer = rpcDetId.layer();

        if (m_debug) std::cout << "Getting RPC global point and translating to DT local coordinates" << std::endl;
        GlobalPoint rpc_gp = getRPCGlobalPosition(rpcDetId, *rpcIt);
        double rpc_global_phi = rpc_gp.phi();
        int rpc_localDT_phi = std::numeric_limits<int>::min();
        // Adaptation of https://github.com/cms-sw/cmssw/blob/master/L1Trigger/L1TTwinMux/src/RPCtoDTTranslator.cc#L349
        if (rpcDetId.sector() == 1) rpc_localDT_phi = int(rpc_global_phi * m_dt_phi_granularity);
        else {
            if (rpc_global_phi >= 0) rpc_localDT_phi = int((rpc_global_phi - (rpcDetId.sector() - 1) * Geom::pi() / 6.) * m_dt_phi_granularity);
            else rpc_localDT_phi = int((rpc_global_phi + (13 - rpcDetId.sector()) * Geom::pi() / 6.) * m_dt_phi_granularity);
        }
        int rpc_phiB = -100000; // single hit has no phiB, DT phiB ranges approx from -1500 to 1500
        int rpc_quality = -1; // to be decided
        int rpc_index = 0;
        int rpc_flag = 3; // only single hit for now
        if (m_debug) std::cout << "Creating DT TP out of RPC recHits" << std::endl;
        rpcRecHits_translated.push_back(L1Phase2MuDTPhDigi(rpc_bx,
                        rpc_wheel,
                        rpc_dt_sector,
                        rpc_station,
                        rpc_layer, //this would be the layer in the new dataformat
                        rpc_localDT_phi,
                        rpc_phiB,
                        rpc_quality,
                        rpc_index,
                        rpc_time,
                        -1, // signle hit --> no chi2
                        rpc_flag
                        ));
    }
}

L1Phase2MuDTPhDigi* RPCIntegrator::matchDTwithRPC(metaPrimitive* dt_metaprimitive, double shift_back) {
    int dt_bx = (int)round(dt_metaprimitive->t0/25.) - shift_back;
    DTChamberId dt_chId = DTChamberId(dt_metaprimitive->rawId);
    int sectorTP = dt_chId.sector();
    if (sectorTP == 13) sectorTP = 4;
    if (sectorTP == 14) sectorTP = 10;
    sectorTP = sectorTP - 1;
    L1Phase2MuDTPhDigi* bestMatch_rpcRecHit = NULL;
    auto bestMatch_rpcRecHit_idx = rpcRecHits_translated.begin();
    float min_dPhi = std::numeric_limits<float>::max();
    for (auto rpcRecHit_translated = rpcRecHits_translated.begin(); rpcRecHit_translated != rpcRecHits_translated.end(); rpcRecHit_translated++) {
        if (rpcRecHit_translated->whNum() == dt_chId.wheel() && rpcRecHit_translated->stNum() == dt_chId.station() && rpcRecHit_translated->scNum() == sectorTP && std::abs(rpcRecHit_translated->bxNum() - dt_bx) <= m_bx_window) {//FIXME improve DT/RPC matching
            //if (min_dPhi != std::numeric_limits<float>::max()){
            //    std::cout << "More than one match for DT/RPC" << dt_chId << std::endl;
            //}
            // Select the RPC hit closest in phi to the DT meta primitive
            if (std::abs(rpcRecHit_translated->phi() - dt_metaprimitive->phi) < min_dPhi){
                min_dPhi = std::abs(rpcRecHit_translated->phi() - dt_metaprimitive->phi);
                bestMatch_rpcRecHit = &*rpcRecHit_translated;
                bestMatch_rpcRecHit_idx = rpcRecHit_translated;
            }
        }
    }
    if (bestMatch_rpcRecHit and !m_storeAllRPCHits){
        rpcRecHits_translated.erase(bestMatch_rpcRecHit_idx);
    }
    return bestMatch_rpcRecHit;
}

void RPCIntegrator::confirmDT(std::vector<metaPrimitive> & dt_metaprimitives, double shift_back) {
    for (auto dt_metaprimitive = dt_metaprimitives.begin(); dt_metaprimitive != dt_metaprimitives.end(); dt_metaprimitive++) {
        L1Phase2MuDTPhDigi* bestMatch_rpcRecHit = matchDTwithRPC(&*dt_metaprimitive, shift_back);
        if (bestMatch_rpcRecHit) {
            (*dt_metaprimitive).rpcFlag = 4;
            if ((*dt_metaprimitive).quality < m_max_quality_to_overwrite_t0){
                (*dt_metaprimitive).t0 = bestMatch_rpcRecHit->t0() + 25 * shift_back; // Overwriting t0 will propagate to BX since it is defined by round((*metaPrimitiveIt).t0/25.)-shift_back
                                                                                      // but we need to add this shift back since all RPC chamber time is centered at 0 for prompt muon
                //(*dt_metaprimitive).phiB = bestMatch_rpcRecHit->phi() - dt_metaprimitive->phi; // Use to fine tune the phi matching window
                (*dt_metaprimitive).rpcFlag = 1;
            }
        }
    }
}
