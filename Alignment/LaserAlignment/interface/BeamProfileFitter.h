#ifndef LaserAlignment_BeamProfileFitter_h
#define LaserAlignment_BeamProfileFitter_h

/** \class BeamProfileFitter
 *  Fitting laser profiles from the beams in the Laser Alignment System
 *
 *  $Date: 2008/01/03 00:53:10 $
 *  $Revision: 1.8 $
 *  \author Maarten Thomas
 */

// Framework headers
#include "FWCore/Framework/interface/eventSetupGetImplementation.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"

// DetId
#include "DataFormats/DetId/interface/DetId.h"

// the result of the fit will be stored into the edm
#include "DataFormats/LaserAlignment/interface/LASBeamProfileFit.h"

// ROOT headers
#include "TVector3.h"
#include "TSpectrum.h"
class TH1D;

#include <string>

class BeamProfileFitter {
 public:
  /// default constructor
   BeamProfileFitter(edm::ParameterSet const& theConf, const edm::EventSetup* aSetup );
  
  /// default destructor
  ~BeamProfileFitter();
  
  /// fitting routine
  std::vector<LASBeamProfileFit> doFit( DetId theDetUnitId, TH1D* theHistogram, 
				        bool theSaveHistograms, 
				        int ScalingFactor, int theBeam, 
				        int theDisc, int theRing, 
				        int theSide, bool isTEC2TEC, 
				        bool & isGoodResult );

  /// the peakfinder
  std::vector<double> findPeakGaus(TH1D* theHist, int theDisc, int theRing);

 private:

  const edm::EventSetup* theSetup;

  bool theClearHistoAfterFit;
  bool theScaleHisto;
  double theMinSignalHeight;
  /// correct for the BS kink?
  bool theCorrectBSkink;
  /// systematic error on measured BS angles
  double theBSAnglesSystematic;
  
  /// function to calculate the error on phi
  Double_t phiError(TVector3 thePosition, TMatrix theCovarianceMatrix);
  
  /// function to calculate an angle between 0 and 2*Pi
  double angle(double theAngle);
};
#endif
