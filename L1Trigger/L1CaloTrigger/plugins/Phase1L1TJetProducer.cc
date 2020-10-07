// -*- C++ -*-
//
// Package: L1CaloTrigger
// Class: Phase1L1TJetProducer
//
/**\class Phase1L1TJetProducer Phase1L1TJetProducer.cc L1Trigger/L1CaloTrigger/plugin/Phase1L1TJetProducer.cc

Description: Produces jets with a phase-1 like sliding window algorithm using a collection of reco::Candidates in input.

*** INPUT PARAMETERS ***
  * etaBinning: vdouble with eta binning (allows non-homogeneous binning in eta)
  * nBinsPhi: uint32, number of bins in phi
  * phiLow: double, min phi (typically -pi)
  * phiUp: double, max phi (typically +pi)
  * jetIEtaSize: uint32, jet cluster size in ieta
  * jetIPhiSize: uint32, jet cluster size in iphi
  * seedPtThreshold: double, threshold of the seed tower
  * puSubtraction: bool, runs chunky doughnut pile-up subtraction, 9x9 jet only
  * outputCollectionName: string, tag for the output collection
  * vetoZeroPt: bool, controls whether jets with 0 pt should be save. 
    It matters if PU is ON, as you can get negative or zero pt jets after it.
  * inputCollectionTag: inputtag, collection of reco::candidates used as input to the algo

*/
//
// Original Simone Bologna
// Created: Mon Jul 02 2018
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCluster.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH2F.h"

#include <algorithm>

class Phase1L1TJetProducer : public edm::one::EDProducer<edm::one::SharedResources> {
   public:
      explicit Phase1L1TJetProducer(const edm::ParameterSet&);
      ~Phase1L1TJetProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);
      
      /// Finds the seeds in the caloGrid, seeds are saved in a vector that contain the index in the TH2F of each seed
      std::vector<std::tuple<int, int>> _findSeeds(const TH2F & caloGrid, float seedThreshold);
      
      std::vector<reco::CaloJet> _buildJetsFromSeedsWithPUSubtraction(const TH2F & caloGrid, const std::vector<std::tuple<int, int>> & seeds, bool killZeroPt);
      std::vector<reco::CaloJet> _buildJetsFromSeeds(const TH2F & caloGrid, const std::vector<std::tuple<int, int>> & seeds);
      
      void _subtract9x9Pileup(const TH2F & caloGrid, reco::CaloJet & jet);
      
      /// Get the energy of a certain tower while correctly handling phi periodicity in case of overflow
      float _getTowerEnergy(const TH2F & caloGrid, int iEta, int iPhi);
      
      reco::CaloJet _buildJetFromSeed(const TH2F & caloGrid, const std::tuple<int, int> & seed);

      // <3 handy method to fill the calogrid with whatever type
      template <class Container>
      void _fillCaloGrid(TH2F & caloGrid, const Container & triggerPrimitives);

      edm::EDGetTokenT<edm::View<reco::Candidate>> *_inputCollectionTag;
      // histogram containing our clustered inputs
      TH2F* _caloGrid;

      std::vector<double> _etaBinning;
      size_t _nBinsEta;
      unsigned int _nBinsPhi;
      double _phiLow;
      double _phiUp;
      unsigned int _jetIEtaSize;
      unsigned int _jetIPhiSize;
      double _seedPtThreshold;
      bool _puSubtraction;
      bool _vetoZeroPt;

      std::string _outputCollectionName;

};


Phase1L1TJetProducer::Phase1L1TJetProducer(const edm::ParameterSet& iConfig)
{

  // getting configuration settings
  this -> _etaBinning = iConfig.getParameter<std::vector<double> >("etaBinning");
  this -> _nBinsEta = this -> _etaBinning.size() - 1;
  this -> _nBinsPhi = iConfig.getParameter<unsigned int>("nBinsPhi");
  this -> _phiLow = iConfig.getParameter<double>("phiLow");
  this -> _phiUp = iConfig.getParameter<double>("phiUp");
  this -> _jetIEtaSize = iConfig.getParameter<unsigned int>("jetIEtaSize");
  this -> _jetIPhiSize = iConfig.getParameter<unsigned int>("jetIPhiSize");
  this -> _seedPtThreshold = iConfig.getParameter<double>("seedPtThreshold");
  this -> _puSubtraction = iConfig.getParameter<bool>("puSubtraction");
  this -> _outputCollectionName = iConfig.getParameter<std::string>("outputCollectionName");
  this -> _vetoZeroPt = iConfig.getParameter<bool>("vetoZeroPt");

  this -> _inputCollectionTag = new edm::EDGetTokenT< edm::View<reco::Candidate> >(consumes< edm::View<reco::Candidate> > (iConfig.getParameter< edm::InputTag >("inputCollectionTag")));
  this -> _caloGrid = new TH2F("caloGrid", "Calorimeter grid", this -> _nBinsEta, this-> _etaBinning.data(), this -> _nBinsPhi, this -> _phiLow, this -> _phiUp);
  this -> _caloGrid -> GetXaxis() -> SetTitle("#eta");
  this -> _caloGrid -> GetYaxis() -> SetTitle("#phi");
  produces<std::vector<reco::CaloJet> >( this->_outputCollectionName ).setBranchAlias(this->_outputCollectionName);
  
  
}

Phase1L1TJetProducer::~Phase1L1TJetProducer()
{
  delete this -> _inputCollectionTag;
  if (this -> _caloGrid) delete this -> _caloGrid;
}

float Phase1L1TJetProducer::_getTowerEnergy(const TH2F & caloGrid, int iEta, int iPhi)
{
  // We return the pt of a certain bin in the calo grid, taking account of the phi periodicity when overflowing (e.g. phi > phiSize), and returning 0 for the eta out of bounds

  int nBinsEta = caloGrid.GetNbinsX();
  int nBinsPhi = caloGrid.GetNbinsY();
  while (iPhi < 1) {
    iPhi += nBinsPhi;
  }
  while (iPhi > nBinsPhi) {
    iPhi -= nBinsPhi;
  }
  if (iEta < 1) {
    //std::cout << "Returning (pt, ieta, iphi): " << "(" << 0 << ", " << iEta << ", " << iPhi << ")" << std::endl;
    return 0;
  }
  if (iEta > nBinsEta) {
    //std::cout << "Returning (pt, ieta, iphi): " << "(" << 0 << ", " << iEta << ", " << iPhi << ")" << std::endl;
    return 0;
  }
    //std::cout << "Returning (pt, ieta, iphi): " << "(" << caloGrid.GetBinContent(iEta, iPhi) << ", " << iEta << ", " << iPhi << ")" << std::endl;
  return caloGrid.GetBinContent(iEta, iPhi);
}

void Phase1L1TJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  
  if (this -> _inputCollectionTag) {
    edm::Handle < edm::View< reco::Candidate > > inputCollectionHandle;
    iEvent.getByToken(*(this -> _inputCollectionTag), inputCollectionHandle);
    // dumping the data
    this -> _caloGrid -> Reset();
    this -> _fillCaloGrid<>(*(this -> _caloGrid), *inputCollectionHandle);

    //int nBinsX = this -> _caloGrid -> GetNbinsX();
    //int nBinsY = this -> _caloGrid -> GetNbinsY();
    //for (int iPhi = 1; iPhi <= nBinsY; iPhi++)
    //{
      // std::cout << "iPhi " << iPhi - 1 << " " << this -> _caloGrid -> GetYaxis() -> GetBinCenter(iPhi) << ": ";
      //for (int iEta = 1; iEta <= nBinsX; iEta++)
      //{
        // std::cout <<this -> _caloGrid -> GetBinContent(iEta, iPhi) << " ";
      //}
      // std::cout << std::endl;
    //}
    const auto seedsVector = this -> _findSeeds(*(this -> _caloGrid), this -> _seedPtThreshold); // seedPtThreshold = 6
    std::vector<reco::CaloJet> l1jetVector;
    if (this -> _puSubtraction)
    {
      l1jetVector = this -> _buildJetsFromSeedsWithPUSubtraction(*(this -> _caloGrid), seedsVector, this -> _vetoZeroPt);
    } else
    {
      l1jetVector = this -> _buildJetsFromSeeds(*(this -> _caloGrid), seedsVector);
    }

    // sort by pt
    std::sort(l1jetVector.begin(), l1jetVector.end(), 
      [](const reco::CaloJet& jet1, const reco::CaloJet& jet2) {
          return jet1.pt() > jet2.pt();   
      }
    );
    
    std::unique_ptr< std::vector<reco::CaloJet> > l1jetVectorPtr(new std::vector<reco::CaloJet>(l1jetVector));
    iEvent.put(std::move(l1jetVectorPtr), this -> _outputCollectionName);
  }

  return;

}

void Phase1L1TJetProducer::_subtract9x9Pileup(const TH2F & caloGrid, reco::CaloJet & jet) {
  //For each jet we compute the 4 side bands
  //for (reco::CaloJet & jet: jetCollection) {
  //for (auto jetIterator = jetCollection.begin(0); jetIterator != jetCollection.end(0); jetIterator++)
  // these variables host the total pt in each sideband and the total pileup contribution
  float topBandPt = 0;
  float leftBandPt = 0;
  float rightBandPt = 0;
  float bottomBandPt = 0;
  float pileUpEnergy;

  // hold the jet's x-y (and z, as I have to use it, even if 2D) location in the histo
  int xCenter, yCenter, zCenter;
  // Retrieving histo-coords for seed
  caloGrid.GetBinXYZ(caloGrid.FindFixBin(jet.eta(), jet.phi()), xCenter, yCenter, zCenter);
  
  // Computing pileup
  for (int x = -4; x <= 4; x++) {
    for (int y = 0; y < 3; y++) {
      // top band, I go up 5 squares to reach the bottom of the top band
      // +x scrolls from left to right, +y scrolls up
      topBandPt += this -> _getTowerEnergy(caloGrid, xCenter + x, yCenter + (5 + y));
      // left band, I go left 5 squares (-5) to reach the bottom of the top band
      // +x scrolls from bottom to top, +y scrolls left
      leftBandPt += this -> _getTowerEnergy(caloGrid, xCenter - (5 + y), yCenter + x);
      // right band, I go right 5 squares (+5) to reach the bottom of the top band
      // +x scrolls from bottom to top, +y scrolls right
      rightBandPt += this -> _getTowerEnergy(caloGrid, xCenter + (5 + y), yCenter + x);
      // right band, I go right 5 squares (+5) to reach the bottom of the top band
      // +x scrolls from bottom to top, +y scrolls right
      bottomBandPt += this -> _getTowerEnergy(caloGrid, xCenter + x, yCenter - (5 + y));
    }
  }
  // adding bands and removing the maximum band (equivalent to adding the three minimum bands)
  pileUpEnergy = topBandPt + leftBandPt + rightBandPt + bottomBandPt - std::max(topBandPt, std::max(leftBandPt, std::max(rightBandPt, bottomBandPt)));

  //preparing the new 4-momentum vector
  //math::PtEtaPhiMLorentzVector ptVector;
  reco::Candidate::PolarLorentzVector ptVector;
  // removing pu contribution
  float ptAfterPUSubtraction = jet.pt() - pileUpEnergy;
  ptVector.SetPt((ptAfterPUSubtraction > 0)? ptAfterPUSubtraction : 0); 
  ptVector.SetEta(jet.eta());
  ptVector.SetPhi(jet.phi());
  //updating the jet
  jet.setP4(ptVector);
  jet.setPileup(pileUpEnergy);
  return;
}

std::vector<std::tuple<int, int>> Phase1L1TJetProducer::_findSeeds(const TH2F & caloGrid, float seedThreshold) 
{
  int nBinsX = caloGrid.GetNbinsX();
  int nBinsY = caloGrid.GetNbinsY();

  std::vector<std::tuple<int, int>> seeds;

  int etaHalfSize = (int) this -> _jetIEtaSize/2;
  int phiHalfSize = (int) this -> _jetIPhiSize/2;

  // for each point of the grid check if it is a local maximum
  // to do so I take a point, and look if is greater than the points around it (in the 9x9 neighborhood)
  // to prevent mutual exclusion, I check greater or equal for points above and right to the one I am considering (including the top-left point)
  // to prevent mutual exclusion, I check greater for points below and left to the one I am considering (including the bottom-right point)

  for (int iPhi = 1; iPhi <= nBinsY; iPhi++)
  {
    for (int iEta = 1; iEta <= nBinsX; iEta++)
    {
      float centralPt = caloGrid.GetBinContent(iEta, iPhi);
      if (centralPt < seedThreshold) continue;
      bool isLocalMaximum = true;

      // Scanning through the grid centered on the seed
      for (int etaIndex = -etaHalfSize; etaIndex <= etaHalfSize; etaIndex++)
      {
        for (int phiIndex = -phiHalfSize; phiIndex <= phiHalfSize; phiIndex++)
        {
          if ((etaIndex == 0) && (phiIndex == 0)) continue;
          if (phiIndex > 0) {
            if (phiIndex > -etaIndex){
              isLocalMaximum = ((isLocalMaximum) && (centralPt > this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex)));
            } else {
              isLocalMaximum = ((isLocalMaximum) && (centralPt >= this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex)));
            }
          } else {
            if (phiIndex >= -etaIndex){
              isLocalMaximum = ((isLocalMaximum) && (centralPt > this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex)));
            } else {
              isLocalMaximum = ((isLocalMaximum) && (centralPt >= this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex)));
            }
          }
        }
      }
      if (isLocalMaximum)
      {
        seeds.emplace_back(std::make_tuple(iEta, iPhi));
      }

    }
  }

  return seeds;
}

reco::CaloJet Phase1L1TJetProducer::_buildJetFromSeed(const TH2F & caloGrid, const std::tuple<int, int> & seed) 
{
  int iEta = std::get<0>(seed);
  int iPhi = std::get<1>(seed);

  int etaHalfSize = (int) this -> _jetIEtaSize/2;
  int phiHalfSize = (int) this -> _jetIPhiSize/2;

  float ptSum = 0;
  // Scanning through the grid centered on the seed
  for (int etaIndex = -etaHalfSize; etaIndex <= etaHalfSize; etaIndex++)
  {
    for (int phiIndex = -phiHalfSize; phiIndex <= phiHalfSize; phiIndex++)
    {
      ptSum += this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex);
    }

  }

  // Creating a jet with eta phi centered on the seed and momentum equal to the sum of the pt of the components    
  //reco::Candidate::LorentzVector ptVector;
  reco::Candidate::PolarLorentzVector ptVector;
  ptVector.SetPt(ptSum);
  //ptVector.SetPtEtaPhiE(ptSum, caloGrid.GetXaxis() -> GetBinCenter(iEta), caloGrid.GetYaxis() -> GetBinCenter(iPhi), ptSum);
  ptVector.SetEta(caloGrid.GetXaxis() -> GetBinCenter(iEta));
  ptVector.SetPhi(caloGrid.GetYaxis() -> GetBinCenter(iPhi));
  // ptVector.SetEta(iEta);
  // ptVector.SetPhi(iPhi);
  reco::CaloJet jet;
  jet.setP4(ptVector);
  return jet;
}

std::vector<reco::CaloJet> Phase1L1TJetProducer::_buildJetsFromSeedsWithPUSubtraction(const TH2F & caloGrid, const std::vector<std::tuple<int, int>> & seeds, bool killZeroPt)
{

  // For each seed take a grid centered on the seed of the size specified by the user
  // Sum the pf in the grid, that will be the pt of the l1t jet. Eta and phi of the jet is taken from the seed.
  std::vector<reco::CaloJet> jets;
  for (const auto& seed: seeds)
  {
    reco::CaloJet jet = this -> _buildJetFromSeed(caloGrid, seed);
    this -> _subtract9x9Pileup(caloGrid, jet);
    //killing jets with 0 pt
    if ((this -> _vetoZeroPt) && (jet.pt() <= 0)) continue;
    jets.push_back(jet);
  }
  return jets;
}

std::vector<reco::CaloJet> Phase1L1TJetProducer::_buildJetsFromSeeds(const TH2F & caloGrid, const std::vector<std::tuple<int, int>> & seeds)
{
  // For each seed take a grid centered on the seed of the size specified by the user
  // Sum the pf in the grid, that will be the pt of the l1t jet. Eta and phi of the jet is taken from the seed.
  std::vector<reco::CaloJet> jets;
  for (const auto& seed: seeds)
  {
    reco::CaloJet jet = this -> _buildJetFromSeed(caloGrid, seed);
    jets.push_back(jet);
    //if(caloGrid.GetEntries() > 0)
    //  std::cout << std::fixed << std::setprecision(2) << jet.pt() << "\t" <<  
    //	jet.eta() << "\t" << jet.phi() << std::endl;
    
  }
  //if(caloGrid.GetEntries() > 0 && jets.size() == 0) std::cout << "0\t0\t0" << std::endl;
  //if(caloGrid.GetEntries() > 0) std::cout << " " << std::endl;
  return jets;
}


template <class Container>
void Phase1L1TJetProducer::_fillCaloGrid(TH2F & caloGrid, const Container & triggerPrimitives)
{
  //Filling the calo grid with the primitives
  for (auto primitiveIterator = triggerPrimitives.begin(); primitiveIterator != triggerPrimitives.end(); primitiveIterator++)
    {
      //if(primitiveIterator->eta() >= 0 && primitiveIterator->eta() < 1.5 && primitiveIterator->phi() >= 0 && primitiveIterator->phi() < 0.7)
      caloGrid.Fill((float) primitiveIterator -> eta(), (float) primitiveIterator -> phi(), (float) primitiveIterator -> pt());
    }
  return;
}

void
Phase1L1TJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(Phase1L1TJetProducer);
