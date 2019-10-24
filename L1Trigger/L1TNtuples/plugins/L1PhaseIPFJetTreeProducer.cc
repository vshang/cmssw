// -*- C++ -*-
//
// Package:    UserCode/L1TriggerDPG
// Class:      L1PhaseIPFJetTreeProducer
// 
/**\class L1PhaseIPFJetTreeProducer L1PhaseIPFJetTreeProducer.cc UserCode/L1TriggerDPG/src/L1PhaseIPFJetTreeProducer.cc

   Description: Produce L1 Extra tree

   Implementation:

*/
//
// Original Author:  Alex Tapper
//         Created:  
// $Id: L1PhaseIPFJetTreeProducer.cc,v 1.5 2013/01/06 21:55:55 jbrooke Exp $
//
//


// system include files
#include <memory>

// framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// data formats
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFJet.h"

// ROOT output stuff
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"

#include "L1Trigger/L1TNtuples/interface/L1AnalysisPhaseIPFJet.h"

//
// class declaration
//

class L1PhaseIPFJetTreeProducer : public edm::EDAnalyzer {
public:
  explicit L1PhaseIPFJetTreeProducer(const edm::ParameterSet&);
  ~L1PhaseIPFJetTreeProducer();


private:
  virtual void beginJob(void) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

public:

  L1Analysis::L1AnalysisPhaseIPFJet* l1Extra;
  L1Analysis::L1AnalysisPhaseIPFJetDataFormat * l1ExtraData;

private:

  unsigned maxL1Extra_;

  // output file
  edm::Service<TFileService> fs_;

  // tree
  TTree * tree_;

  edm::EDGetTokenT<reco::GenJetCollection> genJets_;
  edm::EDGetTokenT<std::vector<l1t::PFJet>> ak4L1PF_;
  edm::EDGetTokenT<std::vector<reco::CaloJet>> phaseIL1PFJets_;

};

L1PhaseIPFJetTreeProducer::L1PhaseIPFJetTreeProducer(const edm::ParameterSet& iConfig){
  
  genJets_ = consumes<reco::GenJetCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genJetToken"));
  ak4L1PF_ = consumes<std::vector<l1t::PFJet> > (iConfig.getParameter<edm::InputTag>("ak4L1PF"));
  phaseIL1PFJets_ = consumes<std::vector<reco::CaloJet> > (iConfig.getParameter<edm::InputTag>("l1PhaseIPFJets"));
	
  maxL1Extra_ = iConfig.getParameter<unsigned int>("maxL1Extra");

  l1Extra     = new L1Analysis::L1AnalysisPhaseIPFJet();
  l1ExtraData = l1Extra->getData();

  // set up output
  tree_=fs_->make<TTree>("L1PhaseIPFJetTree", "L1PhaseIPFJetTree");
  tree_->Branch("L1PhaseIPFJet", "L1Analysis::L1AnalysisPhaseIPFJetDataFormat", &l1ExtraData, 32000, 3);

}


L1PhaseIPFJetTreeProducer::~L1PhaseIPFJetTreeProducer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
L1PhaseIPFJetTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

 l1Extra->Reset();


  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByToken(genJets_, genJets);

  if (genJets.isValid()){ 
    l1Extra->SetGenJet(genJets,maxL1Extra_);
  } else {
    edm::LogWarning("MissingProduct") << "Gen jets not found. Branch will not be filled" << std::endl;
  }


  edm::Handle<std::vector<l1t::PFJet>> ak4L1PFs;
  iEvent.getByToken(ak4L1PF_,ak4L1PFs);


  if (ak4L1PFs.isValid()){
    l1Extra->SetPFJet(ak4L1PFs, maxL1Extra_);
  } else {
    edm::LogWarning("MissingProduct") << "L1PhaseII PFJets not found. Branch will not be filled" << std::endl;
  }



  edm::Handle<std::vector<reco::CaloJet>> phaseIL1PFJets;
  iEvent.getByToken(phaseIL1PFJets_,phaseIL1PFJets);
  if (phaseIL1PFJets.isValid()){
    l1Extra->SetPhaseIPFJet(phaseIL1PFJets, maxL1Extra_);
  } else {
    edm::LogWarning("MissingProduct") << "L1PhaseIPFJet PFJets not found. Branch will not be filled" << std::endl;
  }




  tree_->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
void 
L1PhaseIPFJetTreeProducer::beginJob(void)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1PhaseIPFJetTreeProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1PhaseIPFJetTreeProducer);
