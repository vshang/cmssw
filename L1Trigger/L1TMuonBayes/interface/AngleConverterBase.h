#ifndef ANGLECONVERTER_H
#define ANGLECONVERTER_H

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "L1Trigger/L1TMuonBayes/interface/ProcConfigurationBase.h"

#include <memory>

namespace edm {  
  class EventSetup;
}

class RPCGeometry;
class CSCGeometry;
class CSCLayer;
class DTGeometry;

class L1MuDTChambPhDigi;
class L1MuDTChambThDigi;
class L1MuDTChambThContainer;
class CSCCorrelatedLCTDigi;
class RPCDigi;

class DTChamberId;
class CSCDetId;
class RPCDetId;


struct EtaValue {
  int eta = 0;
  ///error of the eta measurement
  int etaSigma = 0;
  int quality = 0;

  int bx = 0;
  int timing = 0; //sub-bx timing, should be already in scale common for all muon subsystems
};

class AngleConverterBase {
  public:
    AngleConverterBase();
    virtual ~AngleConverterBase();

    ///Update the Geometry with current Event Setup
    virtual void checkAndUpdateGeometry(const edm::EventSetup&,  const ProcConfigurationBase* config);

    /// get phi of DT,CSC and RPC azimutal angle digi in processor scale, used by OMTF algorithm.
    /// in case of wrong phi returns OMTFConfiguration::instance()->nPhiBins
    /// phiZero - desired phi where the scale should start, should be in the desired scale, use getProcessorPhiZero to obtain it
    virtual int getProcessorPhi(int phiZero, l1t::tftype part, const L1MuDTChambPhDigi& digi) const;
    virtual int getProcessorPhi(int phiZero, l1t::tftype part, const CSCDetId& csc, const CSCCorrelatedLCTDigi& digi) const;

    virtual int getProcessorPhi(unsigned int iProcessor, l1t::tftype part, const RPCDetId& rollId, const unsigned int& digi) const;
    virtual int getProcessorPhi(int phiZero, l1t::tftype part, const RPCDetId& rollId, const unsigned int& digi1, const unsigned int& digi2) const;

    ///returns the eta position of the DT chamber
    ///(n.b. in the DT phi and eta segments are independent)
    virtual EtaValue getGlobalEtaDt(const DTChamberId& detId) const;
    
    //adds the eta segments from the thetaDigi to etaSegments
    virtual void getGlobalEta(const L1MuDTChambThDigi& thetaDigi, std::vector<EtaValue>& etaSegments);
    virtual std::vector<EtaValue> getGlobalEta(const L1MuDTChambThContainer* dtThDigis, int bxFrom, int bxTo);

    ///Convert local eta coordinate to global digital microGMT scale.
    virtual EtaValue getGlobalEta(const CSCDetId& detId, const CSCCorrelatedLCTDigi& aDigi);
    
    ///returns the eta position of the CSC chamber
    virtual EtaValue getGlobalEtaCsc(const CSCDetId& detId);

    ///Convert local eta coordinate to global digital microGMT scale.
    ///EtaValue::etaSigma is half of the strip
    virtual EtaValue getGlobalEta(unsigned int rawid, const unsigned int& aDigi);

    float cscChamberEtaSize(const CSCDetId& id);
  protected:

    ///Check orientation of strips in given CSC chamber
    virtual bool isCSCCounterClockwise(const std::unique_ptr<const CSCLayer>& layer) const;


    ///Find BTI group
    virtual const int findBTIgroup(const L1MuDTChambPhDigi &aDigi,
			   const L1MuDTChambThContainer *dtThDigis);
    
    // pointers to the current geometry records
    unsigned long long _geom_cache_id = 0;
    edm::ESHandle<RPCGeometry> _georpc;    
    edm::ESHandle<CSCGeometry> _geocsc;    
    edm::ESHandle<DTGeometry>  _geodt;    

    const ProcConfigurationBase* config = nullptr;
    ///Number of phi bins along 2Pi.
    unsigned int nPhiBins = 0;
   
  };

#endif
