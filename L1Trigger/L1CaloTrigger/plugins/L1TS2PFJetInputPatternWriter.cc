// -*- C++ -*-
//
// Package:    L1Trigger/L1CaloTrigger
// Class:      L1TS2PFJetInputPatternWriter
//
/**\class L1TS2PFJetInputPatternWriter L1TS2PFJetInputPatternWriter.cc L1Trigger/L1TCalorimeter/plugins/L1TS2PFJetInputPatternWriter.cc

   Description: 

   Implementation:

*/
//
// Original Author:  Aaron Bundock
//         Created:  Fri, 26 Jul 2018 14:20:25 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

//
// class declaration
//

class L1TS2PFJetInputPatternWriter : public edm::EDAnalyzer {
public:
  explicit L1TS2PFJetInputPatternWriter(const edm::ParameterSet&);
  ~L1TS2PFJetInputPatternWriter() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<l1t::PFCandidate>> pfToken_;
  std::string filename_;
  std::string outDir_;

  // constants
  unsigned nChan_;  // number of channels per quad
  unsigned nQuad_;
  unsigned nLink_;
  unsigned nHeaderFrames_;
  unsigned nPayloadFrames_;
  unsigned nClearFrames_;
  unsigned nFrame_;
  unsigned nFrameFile_;
  unsigned nEvents_;
  float ptLSB_;
  float phiLSB_;
  float etaLSB_;

  // data arranged by link and frame
  std::vector<std::vector<uint64_t>> data_;

  // data valid flags (just one per frame for now)
  std::vector<int> dataValid_;

  // map of towers onto links/frames
  std::map<int, int> map_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1TS2PFJetInputPatternWriter::L1TS2PFJetInputPatternWriter(const edm::ParameterSet& iConfig)
    : pfToken_(consumes<std::vector<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pfTag"))),
      filename_(iConfig.getUntrackedParameter<std::string>("filename")),
      outDir_(iConfig.getUntrackedParameter<std::string>("outDir")),
      nHeaderFrames_(iConfig.getUntrackedParameter<unsigned>("nHeaderFrames")),
      nPayloadFrames_(iConfig.getUntrackedParameter<unsigned>("nPayloadFrames")),
      nClearFrames_(iConfig.getUntrackedParameter<unsigned>("nClearFrames")) {
  //now do what ever initialization is needed

  // register what you consume and keep token for later access:

  nChan_ = 4;
  nQuad_ = 18;
  ptLSB_ = 0.25;
  etaLSB_ = 0.0043633231;
  phiLSB_ = 0.0043633231;

  nFrame_ = 0;
  nFrameFile_ = 0;
  nEvents_ = 0;

  nLink_ = nChan_ * nQuad_;
  data_.resize(nLink_);
  LogDebug("L1TDebug") << "Preparing for " << nLink_ << " links" << std::endl;
}

L1TS2PFJetInputPatternWriter::~L1TS2PFJetInputPatternWriter() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void L1TS2PFJetInputPatternWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  //count events
  nEvents_++;

  edm::Handle<std::vector<l1t::PFCandidate>> pfHandle;
  iEvent.getByToken(pfToken_, pfHandle);

  std::vector<l1t::PFCandidate> pfPartsA;
  std::vector<l1t::PFCandidate> pfPartsB;

  for (std::vector<l1t::PFCandidate>::const_iterator pfIt = pfHandle->begin(); pfIt != pfHandle->end(); pfIt++) {
    // select first two "small" regions for current fw
    if (pfIt->eta() >= 0 && pfIt->eta() < 0.75 && pfIt->phi() >= 0 && pfIt->phi() < 0.7)
      pfPartsA.push_back(*pfIt);
    if (pfIt->eta() >= 0.75 && pfIt->eta() < 1.5 && pfIt->phi() >= 0 && pfIt->phi() < 0.7)
      pfPartsB.push_back(*pfIt);
  }

  if (pfPartsA.empty() && pfPartsB.empty())
    return;

  if (nFrame_ == 0 || nFrameFile_ == 0) {
    //first empty frames
    while (nFrameFile_ < 14) {
      dataValid_.push_back(1);
      for (unsigned iQuad = 0; iQuad < nQuad_; ++iQuad) {
        for (unsigned iChan = 0; iChan < nChan_; ++iChan) {
          uint iLink = (iQuad * nChan_) + iChan;
          if (iLink == 0)
            data_.at(iLink).push_back(0);
          else
            data_.at(iLink).push_back(0);
          continue;
        }
      }
      nFrame_++;
      nFrameFile_++;
    }
  }

  // loop over frames
  for (unsigned iFrame = 0; iFrame < nPayloadFrames_; ++iFrame) {
    dataValid_.push_back(1);
    // loop over links
    for (unsigned iQuad = 0; iQuad < nQuad_; ++iQuad) {
      for (unsigned iChan = 0; iChan < nChan_; ++iChan) {
        // get tower ieta, iphi for link
        uint iLink = (iQuad * nChan_) + iChan;

        uint64_t data = 0;

        if ((nFrameFile_ % 13) == 1) {
          if (iLink < 24 && pfPartsA.size() > iLink) {
            data |= ((uint64_t)floor(pfPartsA.at(iLink).pt() / ptLSB_) & 0xffff);
            data |= ((uint64_t)floor(pfPartsA.at(iLink).phi() / phiLSB_) & 0x3ff) << 16;
            data |= ((uint64_t)floor(pfPartsA.at(iLink).eta() / etaLSB_) & 0x3ff) << 26;
          }
        }
        if ((nFrameFile_ % 13) == 2) {
          if (iLink < 24 && pfPartsB.size() > iLink) {
            data |= ((uint64_t)floor(pfPartsB.at(iLink).pt() / ptLSB_) & 0xffff);
            data |= ((uint64_t)floor(pfPartsB.at(iLink).phi() / phiLSB_) & 0x3ff) << 16;
            data |= ((uint64_t)floor((pfPartsB.at(iLink).eta() - 0.75) / etaLSB_) & 0x3ff) << 26;
          }
        }
        // add data to output
        data_.at(iLink).push_back(data);
      }
    }
    nFrame_++;
    nFrameFile_++;
    if (nFrame_ % 1015 == 0)
      nFrameFile_ = 0;
  }
}

// ------------ method called once each job just before starting event loop  ------------
void L1TS2PFJetInputPatternWriter::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void L1TS2PFJetInputPatternWriter::endJob() {
  //frames per event
  unsigned int framesPerEv = nHeaderFrames_ + nPayloadFrames_ + nClearFrames_;

  //frames per file
  unsigned int framesPerFile = 1015;

  //events per file
  unsigned int evPerFile = floor(framesPerFile / framesPerEv);

  //number of output files
  unsigned int nOutFiles = ceil((float)nEvents_ / (float)evPerFile);

  LogDebug("L1TDebug") << "Read " << nFrame_ << " frames" << std::endl;
  LogDebug("L1TDebug") << "Read " << nEvents_ << " events" << std::endl;
  LogDebug("L1TDebug") << "Writing " << nOutFiles << " files" << std::endl;
  LogDebug("L1TDebug") << "Output directory: ./" << outDir_ << "/" << std::endl;

  //files
  std::vector<std::ofstream> outFiles(nOutFiles);

  //make output files and write to them
  for (uint itFile = 0; itFile < nOutFiles; ++itFile) {
    std::stringstream outFilename;
    outFilename << outDir_ << "/" << filename_ << "_" << itFile << ".txt";
    outFiles[itFile] = std::ofstream(outFilename.str());
    LogDebug("L1TDebug") << "Writing to file: ./" << outFilename.str() << std::endl;
    std::cout << "Writing to file: ./" << outFilename.str() << std::endl;

    outFiles[itFile] << "Board SRNTY_TEST" << std::endl;

    // quad/chan numbers
    outFiles[itFile] << " Quad/Chan :      ";
    for (unsigned i = 0; i < nQuad_; ++i) {
      for (unsigned j = 0; j < nChan_; ++j) {
        outFiles[itFile] << "  q" << setfill('0') << setw(2) << i << "c" << j << "            ";
      }
    }
    outFiles[itFile] << std::endl;

    // link numbers
    outFiles[itFile] << "      Link :     ";
    for (unsigned i = 0; i < nQuad_; ++i) {
      for (unsigned j = 0; j < nChan_; ++j) {
        outFiles[itFile] << "    " << setfill('0') << setw(2) << (i * nChan_) + j << "             ";
      }
    }

    outFiles[itFile] << std::endl;

    // then the data
    unsigned iFileFrame = 0;
    for (unsigned iFrame = itFile * framesPerFile; iFrame < (itFile * framesPerFile + framesPerFile); ++iFrame) {
      if (iFrame <= nFrame_ && iFrame < (framesPerEv * nEvents_)) {
        outFiles[itFile] << "Frame " << std::dec << std::setw(4) << std::setfill('0') << iFileFrame << " : ";
        for (unsigned iQuad = 0; iQuad < nQuad_; ++iQuad) {
          for (unsigned iChan = 0; iChan < nChan_; ++iChan) {
            unsigned iLink = (iQuad * nChan_) + iChan;
            if (iLink < data_.size() && iFrame < data_.at(iLink).size()) {
              outFiles[itFile] << std::hex << ::std::setw(1) << dataValid_.at(iFrame) << "v" << std::hex
                               << std::setw(16) << std::setfill('0') << data_.at(iLink).at(iFrame) << " ";
            } else {
              outFiles[itFile] << std::hex << ::std::setw(1) << 0 << "v" << std::hex << std::setw(16)
                               << std::setfill('0') << 0 << " ";
            }
          }
        }
      }
      outFiles[itFile] << std::endl;
      iFileFrame++;
    }
    outFiles[itFile].close();
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1TS2PFJetInputPatternWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TS2PFJetInputPatternWriter);
