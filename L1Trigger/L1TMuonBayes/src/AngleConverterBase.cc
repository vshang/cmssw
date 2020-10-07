#include "L1Trigger/L1TMuonBayes/interface/AngleConverterBase.h"
#include "L1Trigger/L1TMuonBayes/interface/Omtf/OMTFConfiguration.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCConstants.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCPatternLUT.h"

#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "L1Trigger/DTUtilities/interface/DTTrigGeom.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"

#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <cmath> 

namespace {
template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

int fixCscOffsetGeom(int offsetLoc) {
  // fix for CSC feo dependance from GlobalTag

  // dump of CSC offsets for MC global tag
  const std::vector<int> offCSC = { -154, -133, -17, -4, 4, 17, 133, 146, 154, 167, 
                            283, 296, 304, 317, 433, 446, 454, 467,  
                            583, 596, 604, 617, 733, 746, 754, 767,  883, 904};
  auto gep = std::lower_bound(offCSC.begin(), offCSC.end(), offsetLoc);
  int fixOff = (gep != offCSC.end()) ? *gep : *(gep-1); 
  if (gep != offCSC.begin() && std::abs(*(gep-1)-offsetLoc) < std::abs(fixOff-offsetLoc)  ) fixOff= *(gep-1);
  return fixOff;
}

}


AngleConverterBase::AngleConverterBase(): _geom_cache_id(0ULL) { }
///////////////////////////////////////
///////////////////////////////////////
AngleConverterBase::~AngleConverterBase() {  }
///////////////////////////////////////
///////////////////////////////////////
void AngleConverterBase::checkAndUpdateGeometry(const edm::EventSetup& es,  const ProcConfigurationBase* config) {
  const MuonGeometryRecord& geom = es.get<MuonGeometryRecord>();
  unsigned long long geomid = geom.cacheIdentifier();
  if( _geom_cache_id != geomid ) {
    geom.get(_georpc);  
    geom.get(_geocsc);    
    geom.get(_geodt);
    _geom_cache_id = geomid;
  }
  this->config = config;
  nPhiBins = config->nPhiBins();
}
///////////////////////////////////////
///////////////////////////////////////
int AngleConverterBase::getProcessorPhi(int phiZero, l1t::tftype part, const L1MuDTChambPhDigi &digi) const
{
  int dtPhiBins = 4096;

  double hsPhiPitch = 2*M_PI/nPhiBins; // width of phi Pitch, related to halfStrip at CSC station 2

  int sector  = digi.scNum()+1;   //NOTE: there is a inconsistency in DT sector numb. Thus +1 needed to get detector numb.
  //int wheel   = digi.whNum();
  //int station = digi.stNum();
  int phiDT   = digi.phi();

  //int offsetLoc = lround( ((ichamber-1)*M_PI/6+M_PI/12.)/hsPhiPitch );
  double scale = 1./dtPhiBins/hsPhiPitch;
  int scale_coeff = lround(scale* pow(2,11));
//  int phi = static_cast<int>(phiDT*scale) + offsetLoc;

  int ichamber = sector-1;
  if(ichamber > 6)
    ichamber = ichamber - 12;

  int offsetGlobal = (int)nPhiBins  * ichamber / 12;

  int phi = floor(phiDT*scale_coeff/pow(2,11)) + offsetGlobal - phiZero;

  //std::cout<<__FUNCTION__<<":"<<__LINE__<<" phiZero "<<phiZero<<" phiDT "<<phiDT<<" sector "<<sector<<" ichamber "<<ichamber<<" offsetGlobal "<<offsetGlobal<<" phi "<<phi<<" foldPhi(phi) "<<omtfConfig->foldPhi(phi)<<std::endl;
  return config->foldPhi(phi);
}
///////////////////////////////////////
///////////////////////////////////////
int AngleConverterBase::getProcessorPhi(int phiZero, l1t::tftype part, const CSCDetId & csc, const CSCCorrelatedLCTDigi &digi) const
{
  const double hsPhiPitch = 2*M_PI/nPhiBins;
  //
  // get offset for each chamber.
  // FIXME: These parameters depends on processor and chamber only so may be precomputed and put in map
  //
  const CSCChamber* chamber = _geocsc->chamber(csc);
  const CSCChamberSpecs* cspec = chamber->specs();
  const CSCLayer* layer = chamber->layer(3);
  int order = ( layer->centerOfStrip(2).phi() - layer->centerOfStrip(1).phi()> 0) ? 1 : -1;
  double stripPhiPitch = cspec->stripPhiPitch();
  double scale = fabs(stripPhiPitch/hsPhiPitch/2.);
  if ( fabs(scale-1.) < 0.0002) scale=1.;

  //double phiHalfStrip0 = layer->centerOfStrip(10).phi() - order*9*stripPhiPitch - order*stripPhiPitch/4.; //KB are the strips enumerated from 1??? - only then this has sense

  double phiHalfStrip0 = layer->centerOfStrip(1).phi() - order*stripPhiPitch/4.;

  int offsetLoc = lround( (phiHalfStrip0)/hsPhiPitch - phiZero);
  offsetLoc = config->foldPhi(offsetLoc);

  int halfStrip = digi.getStrip(); // returns halfStrip 0..159

  //FIXME: to be checked (only important for ME1/3) keep more bits for offset, truncate at the end

  // a quick fix for towards geometry changes due to global tag.
  // in case of MC tag fixOff shold be identical to offsetLoc 
  //int fixOff = fixCscOffsetGeom(offsetLoc); TODO does not work in correlator, i.e. when phiZero is always 0. .Fix this

  int fixOff = offsetLoc;

  int phi = fixOff + order*scale*halfStrip;

//  std::cout <<" phiZero "<<phiZero<<" hs: "<< halfStrip <<" phiHalfStrip0 "<<phiHalfStrip0<<" offset: " << offsetLoc <<" oder*scale: "<< order*scale
//            <<" phi: " <<phi<<" foldPhi(phi) "<<config->foldPhi(phi)<<" ("<<offsetLoc + order*scale*halfStrip<<")"<< std::endl;

  return config->foldPhi(phi);
}

///////////////////////////////////////
///////////////////////////////////////
int AngleConverterBase::getProcessorPhi(int phiZero, l1t::tftype part, const RPCDetId & rollId, const unsigned int &digi1, const unsigned int &digi2) const
{
  const double hsPhiPitch = 2*M_PI/nPhiBins;
  const int dummy = nPhiBins;
  const RPCRoll* roll = _georpc->roll(rollId);
  if (!roll) return dummy;

  double stripPhi1 = (roll->toGlobal(roll->centreOfStrip((int)digi1))).phi(); // note [-pi,pi]
  double stripPhi2 = (roll->toGlobal(roll->centreOfStrip((int)digi2))).phi(); // note [-pi,pi]
  // stripPhi from geometry is given in [-pi,pi] range.  

  // local angle in CSC halfStrip usnits
  int halfStrip = lround ( ( (stripPhi1+stripPhi2)/2.)/hsPhiPitch);
  halfStrip = config->foldPhi(halfStrip); //only for the case when the two strips are on different sides of phi = pi

//  LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" roll "<<rollId<<" cluster: firstStrip "<<digi1<<" stripPhi1 "<<stripPhi1
//      <<" lastStrip "<<digi2<<" stripPhi2 "<<stripPhi2<<" halfStrip "<<halfStrip<<std::endl;

  return config->foldPhi(halfStrip - phiZero);
}

int AngleConverterBase::getProcessorPhi(unsigned int iProcessor, l1t::tftype part, const RPCDetId & rollId, const unsigned int &digi) const
{

  const double hsPhiPitch = 2*M_PI/nPhiBins;
  const int dummy = nPhiBins;
  int processor = iProcessor+1;
  const RPCRoll* roll = _georpc->roll(rollId);
  if (!roll) return dummy;

  double phi15deg =  M_PI/3.*(processor-1)+M_PI/12.;                    // "0" is 15degree moved cyclicaly to each processor, note [0,2pi]
  double stripPhi = (roll->toGlobal(roll->centreOfStrip((int)digi))).phi(); // note [-pi,pi]

  // adjust [0,2pi] and [-pi,pi] to get deltaPhi difference properly
  switch (processor) {
    case 1: break;
    case 6: {phi15deg -= 2*M_PI; break; }
    default : {if (stripPhi < 0) stripPhi += 2*M_PI; break; }
  }

  // local angle in CSC halfStrip usnits
  int halfStrip = lround ( (stripPhi-phi15deg)/hsPhiPitch );
    
  return halfStrip;
}
///////////////////////////////////////
///////////////////////////////////////
EtaValue AngleConverterBase::getGlobalEtaDt(const DTChamberId& detId) const {
  // do not use this pointer for anything other than creating a trig geom
  std::unique_ptr<DTChamber> chamb(const_cast<DTChamber*>(_geodt->chamber(detId)));
  
  Local2DPoint chamberMiddleLP(0, 0);
  GlobalPoint chamberMiddleGP = chamb->toGlobal(chamberMiddleLP);
  chamb.release();

  const DTChamberId baseidNeigh(detId.wheel() + (detId.wheel() >= 0 ? -1 : +1), detId.station(), detId.sector());
  std::unique_ptr<DTChamber> chambNeigh(const_cast<DTChamber*>(_geodt->chamber(baseidNeigh)));
  GlobalPoint chambNeighMiddleGP = chambNeigh->toGlobal(chamberMiddleLP);
  chambNeigh.release();

  EtaValue etaValue = {
      config->etaToHwEta(chamberMiddleGP.eta() ),
      config->etaToHwEta( abs(chamberMiddleGP.eta() - chambNeighMiddleGP.eta() ) ) / 2,
      0, //quality
      0, //bx
      0 //timin
  };


  //std::cout<<__FUNCTION__<<":"<<__LINE__<<" rawid "<<detId.rawId()<<" baseid "<<detId<<" chamberMiddleGP.eta() "<<chamberMiddleGP.eta()<<" eta "<<etaValue.eta<<" etaSigma "<<etaValue.etaSigma<<std::endl;
  return etaValue;
}

///////////////////////////////////////
///////////////////////////////////////
void AngleConverterBase::getGlobalEta(const L1MuDTChambThDigi& thetaDigi, std::vector<EtaValue>& etaSegments) {
  const DTChamberId baseid(thetaDigi.whNum(), thetaDigi.stNum(), thetaDigi.scNum()+1);
  // do not use this pointer for anything other than creating a trig geom
  std::unique_ptr<DTChamber> chamb(const_cast<DTChamber*>(_geodt->chamber(baseid)));

  std::unique_ptr<DTTrigGeom> trig_geom( new DTTrigGeom(chamb.get(),false) );
  chamb.release(); // release it here so no one gets funny ideas
  // super layer 2 is the theta superlayer in a DT chamber
  // station 4 does not have a theta super layer
  // the BTI index from the theta trigger is an OR of some BTI outputs
  // so, we choose the BTI that's in the middle of the group
  // as the BTI that we get theta from
  // TODO:::::>>> need to make sure this ordering doesn't flip under wheel sign
  const int NBTI_theta = trig_geom->nCell(2);
  for(unsigned int btiGroup = 0; btiGroup < 7; ++btiGroup ) {
    if(thetaDigi.position(btiGroup) ) {
      unsigned btiActual = btiGroup * NBTI_theta/7 + NBTI_theta/14 + 1;
      DTBtiId thetaBTI = DTBtiId(baseid, 2, btiActual);
      GlobalPoint theta_gp = trig_geom->CMSPosition(thetaBTI);

      EtaValue etaValue = {
          config->etaToHwEta(theta_gp.eta() ),
          0,
          thetaDigi.quality(btiGroup),
          thetaDigi.bxNum(),
          0 //TODO what about sub-bx timing???
      };
      etaSegments.emplace_back(etaValue);


      //std::cout<<__FUNCTION__<<":"<<__LINE__<<" bx "<<thetaDigi.bxNum()<<" baseid "<<baseid<<" btiGroup "<<btiGroup<<" quality "<<thetaDigi.quality(btiGroup)<<" theta_gp.eta() "<<theta_gp.eta()<<" eta "<<etaValue.eta<<" etaSigma "<<etaValue.etaSigma<<std::endl;
    }
  }
}

std::vector<EtaValue> AngleConverterBase::getGlobalEta(const L1MuDTChambThContainer* dtThDigis, int bxFrom, int bxTo) {
  //std::cout<<__FUNCTION__<<":"<<__LINE__<<" dtThDigis size "<<dtThDigis->getContainer()->size()<<std::endl;

  std::vector<EtaValue> etaSegments;

  for(auto& thetaDigi: (*(dtThDigis->getContainer()) ) ) {
    if(thetaDigi.bxNum() >= bxFrom && thetaDigi.bxNum() <= bxTo) {
      getGlobalEta(thetaDigi, etaSegments);
    }
  }
  return etaSegments;
}

//just read from the drawing
float AngleConverterBase::cscChamberEtaSize(const CSCDetId& detId) {
  if(detId.station() == 1) {
    if(detId.ring() == 1) return (2.5 - 1.6)/2.; ///ME1/1 lower eta (b?, eta < ~2.1), muCorrelator eta bins 6-11 - but getGlobalEtaCsc(const CSCDetId& detId) gives the midle of the full chamber, so here we put the size of the full chamber
    if(detId.ring() == 2) return (1.7 - 1.2)/2.;
    if(detId.ring() == 3) return (1.12- 0.9)/2.;
    if(detId.ring() == 4) return (2.5 - 1.6)/2.; ///ME1/1 higher eta (a?, eta > ~2.1), muCorrelator eta bins 10-15
  }
  else if(detId.station() == 2) {
    if(detId.ring() == 1) return (2.5 - 1.6)/2.;
    if(detId.ring() == 2) return (1.6 - 1.0)/2.;
  }
  else if(detId.station() == 3) {
    if(detId.ring() == 1) return (2.5 - 1.7)/2.;
    if(detId.ring() == 2) return (1.7 - 1.1)/2.;
  }
  else if(detId.station() == 4) {
    if(detId.ring() == 1) return (2.45- 1.8)/2.;
    if(detId.ring() == 2) return (1.8 - 1.2)/2.;
  }
  return 0;
}

EtaValue AngleConverterBase::getGlobalEta(const CSCDetId& detId, const CSCCorrelatedLCTDigi &aDigi){
   ///Code taken from GeometryTranslator.
  ///Will be replaced by direct CSC phi local to global scale
  ///transformation as used in FPGA implementation

  
  // alot of this is transcription and consolidation of the CSC
  // global phi calculation code
  // this works directly with the geometry 
  // rather than using the old phi luts

  // we should change this to weak_ptrs at some point
  // requires introducing std::shared_ptrs to geometry
  std::unique_ptr<const CSCChamber> chamb(_geocsc->chamber(detId));
  std::unique_ptr<const CSCLayerGeometry> layer_geom(
						     chamb->layer(CSCConstants::KEY_ALCT_LAYER)->geometry()
						     );
  std::unique_ptr<const CSCLayer> layer(
					chamb->layer(CSCConstants::KEY_ALCT_LAYER)
					);

  const uint16_t keyWG = aDigi.getKeyWG();
/*  const uint16_t halfstrip = aDigi.getStrip();
  const uint16_t pattern = aDigi.getPattern();

  //const unsigned maxStrips = layer_geom->numberOfStrips();

  // so we can extend this later 
  // assume TMB2007 half-strips only as baseline
  double offset = 0.0;
  switch(1) {
  case 1:
    offset = CSCPatternLUT::get2007Position(pattern);
  }
  const unsigned halfstrip_offs = unsigned(0.5 + halfstrip + offset);
  const unsigned strip = halfstrip_offs/2 + 1; // geom starts from 1

  // the rough location of the hit at the ALCT key layer
  // we will refine this using the half strip information
  const LocalPoint coarse_lp = 
    layer_geom->stripWireGroupIntersection(strip,keyWG);  
  const GlobalPoint coarse_gp = layer->surface().toGlobal(coarse_lp);  
  
  // the strip width/4.0 gives the offset of the half-strip
  // center with respect to the strip center
  const double hs_offset = layer_geom->stripPhiPitch()/4.0;
  
  // determine handedness of the chamber
  const bool ccw = isCSCCounterClockwise(layer);
  // we need to subtract the offset of even half strips and add the odd ones
  const double phi_offset = ( ( halfstrip_offs%2 ? 1 : -1)*
			      ( ccw ? -hs_offset : hs_offset ) );
  
  // the global eta calculation uses the middle of the strip
  // so no need to increment it
  const GlobalPoint final_gp( GlobalPoint::Polar( coarse_gp.theta(),
						  (coarse_gp.phi().value() + 
						   phi_offset),
						  coarse_gp.mag() ) );*/

  const LocalPoint lpWg = layer_geom->localCenterOfWireGroup(keyWG);
  const GlobalPoint gpWg = layer->surface().toGlobal(lpWg);

  // release ownership of the pointers
  chamb.release();
  layer_geom.release();
  layer.release();

//  std::cout <<id<<" st: " << id.station()<< "ri: "<<id.ring()<<" eta: " <<  final_gp.eta() 
//           <<" etaCode_simple: " <<  etaVal2Code( final_gp.eta() )<< " KW: "<<keyWG <<" etaKeyWG2Code: "<<etaKeyWG2Code(id,keyWG)<< std::endl;
//  int station = (id.endcap()==1) ? id.station() : -id.station();
//  std::cout <<"ETA_CSC: " << station <<" "<<id.ring()<<" "<< final_gp.eta()<<" "<<keyWG <<" "<< etaKeyWG2Code(id,keyWG) << std::endl;

   EtaValue etaSegment = {
       config->etaToHwEta(gpWg.eta() ),
       0, //config->etaToHwEta(cscChamberEtaSize(id) ),
       0,
       aDigi.getBX(),
       0 //tming???
   };

   //std::cout<<__FUNCTION__<<":"<<__LINE__<<" csc "<<detId<<" eta "<<gpWg.eta()<<" etaHw "<<etaSegment.eta<<" etaSigma "<<etaSegment.etaSigma<<std::endl;
   return etaSegment;

}


//TODO the CSC ME1/1 has strips divided int tow part a nad b, so this function in principle ca include that, then it should also receive the roll number as parameter, off course implementation should be different then
EtaValue AngleConverterBase::getGlobalEtaCsc(const CSCDetId& detId) {
  std::unique_ptr<const CSCChamber> chamb(_geocsc->chamber(detId));

/*  std::unique_ptr<const CSCLayerGeometry> layer_geom(
                 chamb->layer(CSCConstants::KEY_ALCT_LAYER)->geometry()
                 );
  std::unique_ptr<const CSCLayer> layer(
          chamb->layer(CSCConstants::KEY_ALCT_LAYER)
          );*/

  Local2DPoint chamberMiddleLP(0, 0);
  GlobalPoint chamberMiddleGP = chamb->toGlobal(chamberMiddleLP);
  chamb.release();

  EtaValue etaValue = {
      config->etaToHwEta(chamberMiddleGP.eta() ),
      config->etaToHwEta(cscChamberEtaSize(detId) ),
      0,
      0, //bx
      0 //timnig
  };

  //std::cout<<__FUNCTION__<<":"<<__LINE__<<" rawid "<<detId.rawId()<<" detId "<<detId<<" chamberMiddleGP.eta() "<<chamberMiddleGP.eta()<<" eta "<<etaValue.eta<<" etaSigma "<<etaValue.etaSigma<<std::endl;
  return etaValue;

}

///////////////////////////////////////
///////////////////////////////////////
EtaValue AngleConverterBase::getGlobalEta(unsigned int rawid, const unsigned int &strip) {
  const RPCDetId id(rawid);

  std::unique_ptr<const RPCRoll>  roll(_georpc->roll(id));
  const LocalPoint lp = roll->centreOfStrip((int)strip);
  const GlobalPoint gp = roll->toGlobal(lp);

  int neighbRoll = 1; //neighbor roll in eta
  //roll->chamber()->nrolls() does not work
  if(id.region() == 0) {//barel
    if( id.station() == 2 && ( (abs(id.ring()) == 2 && id.layer() == 2 ) || (abs(id.ring()) != 2 && id.layer() == 1 ) ) ) { //three-roll chamber
      if(id.roll() == 2)
        neighbRoll = 1;
      else {
        neighbRoll = 2;
      }
    }
    else //two-roll chamber
      neighbRoll = (id.roll() == 1 ? 3 : 1 );
  }
  else {//endcap
    neighbRoll = id.roll() + (id.roll() == 1 ? +1 : -1);
  }
  roll.release();

  const RPCDetId idNeigh = RPCDetId(id.region(), id.ring(), id.station(), id.sector(), id.layer(), id.subsector(), neighbRoll );
  //std::cout<<__FUNCTION__<<":"<<__LINE__<<" rpc "<<id<<std::endl;
  //std::cout<<__FUNCTION__<<":"<<__LINE__<<" rpc "<<idNeigh<<std::endl;
  std::unique_ptr<const RPCRoll>  rollNeigh(_georpc->roll(idNeigh));
  const LocalPoint lpNeigh = rollNeigh->centreOfStrip((int)strip);
  const GlobalPoint gpNeigh = rollNeigh->toGlobal(lpNeigh);
  rollNeigh.release();


  EtaValue etaValue = {
      config->etaToHwEta(gp.eta()),
      config->etaToHwEta( abs(gp.eta() - gpNeigh.eta()) )/2, //half of the size of the strip in eta - not precise, but OK
      0
  };

  //std::cout<<__FUNCTION__<<":"<<__LINE__<<" rpc "<<id<<" eta "<<gp.eta()<<" etaHw "<<etaValue.eta<<" etaSigma "<<etaValue.etaSigma<<std::endl;
  return etaValue;
}
///////////////////////////////////////
///////////////////////////////////////
bool AngleConverterBase::
isCSCCounterClockwise(const std::unique_ptr<const CSCLayer>& layer) const {
  const int nStrips = layer->geometry()->numberOfStrips();
  const double phi1 = layer->centerOfStrip(1).phi();
  const double phiN = layer->centerOfStrip(nStrips).phi();
  return ( (std::abs(phi1 - phiN) < M_PI  && phi1 >= phiN) || 
	   (std::abs(phi1 - phiN) >= M_PI && phi1 < phiN)     );  
}
///////////////////////////////////////
///////////////////////////////////////
const int AngleConverterBase::findBTIgroup(const L1MuDTChambPhDigi &aDigi,
				       const L1MuDTChambThContainer *dtThDigis){

  int bti_group = -1;
  
  const L1MuDTChambThDigi *theta_segm = dtThDigis->chThetaSegm(aDigi.whNum(),
							       aDigi.stNum(),
							       aDigi.scNum(),
							       aDigi.bxNum());
  if(!theta_segm) return  bti_group;
  
  for(unsigned int i = 0; i < 7; ++i ){
    if(theta_segm->position(i) && bti_group<0) bti_group = i;
    ///If there are more than one theta digi we do not take is
    ///due to unresolvet ambiguity. In this case we take eta of the
    ///middle of the chamber.
    else if(theta_segm->position(i) && bti_group>-1) return -1;
  }
      
  return bti_group;
}
///////////////////////////////////////
///////////////////////////////////////
