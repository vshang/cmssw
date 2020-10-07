#ifndef OMTF_GOLDENPATTERNRESULTS_H
#define OMTF_GOLDENPATTERNRESULTS_H

#include <L1Trigger/L1TMuonBayes/interface/Omtf/OMTFConfiguration.h>
#include <L1Trigger/L1TMuonBayes/interface/StubResult.h>

#include <vector>
#include <ostream>

//result for one refHit of one GoldenPattern
class GoldenPatternResult {
public:

/*  struct LayerResult {
    PdfValueType pdfVal = 0;
    bool valid = false;

    int pdfBin = 0; //hit deltaPhi, bin=0 is reserved for no valid hit, see GoldenPatternBase::process1Layer1RefLayer, is obtained from it, see GoldenPatternBase::process1Layer1RefLayer
    int hit = 0; //original hit phi

    LayerResult(PdfValueType pdfVal, bool valid, int pdfBin, int hit) :
      pdfVal(pdfVal), valid(valid), pdfBin(pdfBin), hit(hit) {
    };
  };*/

private:
  bool valid = false;

  //number of the layer from which the reference hit originated
  int refLayer = 0;

  ///phi at the 2nd muon station (propagated refHitPhi)
  int phi = 0;

  ///eta at the 2nd muon station
  int eta = 0;

  ///Sum of pdfValues
  //omtfPdfValueType
  double pdfSum = 0;

  ///Number of fired layers - excluding banding layers
  unsigned int firedLayerCnt = 0;

  ///bits representing fired logicLayers (including banding layers),
  unsigned int firedLayerBits = 0;

  ///phi of the reference hits
  int refHitPhi = 0;

  static int finalizeFunction;

  double gpProbability1 = 0;

  double gpProbability2 = 0;

  StubResults stubResults;
public:
  void init(const OMTFConfiguration* omtfConfig);

  void reset();

  bool isValid() const {
    return valid;
  }

  void setValid(bool valid) {
    this->valid = valid;
  }

  void set(int refLayer, int phi, int eta, int refHitPhi);

  //void setLayerResult(unsigned int iLayer, LayerResult layerResult);

  void setStubResult(float pdfVal, bool valid, int pdfBin, int layer, MuonStubPtr stub);

  void setStubResult(int layer, StubResult& stubResult);

  int getRefLayer() const {
    return this->refLayer;
  }

  void setRefLayer(int refLayer) {
    this->refLayer = refLayer;
  }

  int getEta() const {
    return eta;
  }

  void setEta(int eta) {
    this->eta = eta;
  }

  unsigned int getFiredLayerBits() const {
    return firedLayerBits;
  }

  void setFiredLayerBits(unsigned int firedLayerBits) {
    this->firedLayerBits = firedLayerBits;
  }

  unsigned int getFiredLayerCnt() const {
    return firedLayerCnt;
  }

  void setFiredLayerCnt(unsigned int firedLayerCnt) {
    this->firedLayerCnt = firedLayerCnt;
  }

  /*
   * pdfValue from each layer
   */
/*  const std::vector<PdfValueType>& getPdfValues() const {
    return pdfValues;
  }*/

  /*
   * sum of the pdfValues in layers
   * if finalise2() it is product of the pdfValues
   */
  PdfValueType getPdfSum() const {
    return pdfSum;
  }

/*  const std::vector<int>& getHitPdfBins() const {
    return hitPdfBins;
  }

  const std::vector<int>& getHits() const {
    return hits;
  }*/

  const StubResults& getStubResults() const {
    return stubResults;
  }

  int getPhi() const {
    return phi;
  }

  void setPhi(int phi) {
    this->phi = phi;
  }

  int getRefHitPhi() const {
    return refHitPhi;
  }

  void setRefHitPhi(int refHitPhi) {
    this->refHitPhi = refHitPhi;
  }

  bool isLayerFired(unsigned int iLayer) const {
    return firedLayerBits & (1<<iLayer);
  }

  GoldenPatternResult() {
  };

  //dont use this in the pattern construction, since the myOmtfConfig is null then
  GoldenPatternResult(const OMTFConfiguration* omtfConfig);

  void set();

  void finalise() {
    if(finalizeFunction == 1)
      finalise1();
    else if(finalizeFunction == 2)
      finalise2();
    else if(finalizeFunction == 3)
      finalise3();
    else
      finalise0();
  }

  //version for the "normal" patterns, i.e. without pdfSum threshold
  void finalise0();

  //version for the patterns with pdfSum threshold
  void finalise1();

  //multiplication of PDF values instead of sum
  void finalise2();

  //for patterns generation
  void finalise3();

  //bool empty() const;

  friend std::ostream & operator << (std::ostream &out, const GoldenPatternResult & aResult);

  static void setFinalizeFunction(int finalizeFunction_) {
    finalizeFunction = finalizeFunction_;
  }

  double getGpProbability1() const {
    return gpProbability1;
  }

  void setGpProbability1(double probability1 = 0) {
    this->gpProbability1 = probability1;
  }

  double getGpProbability2() const {
    return gpProbability2;
  }

  void setGpProbability2(double probability2 = 0) {
    this->gpProbability2 = probability2;
  }

private:

  const OMTFConfiguration* omtfConfig =  nullptr;

};


#endif //OMTF_GOLDENPATTERNRESULTS_H
