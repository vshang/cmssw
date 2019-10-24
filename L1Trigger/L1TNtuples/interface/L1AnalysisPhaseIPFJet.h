#ifndef __L1Analysis_L1AnalysisPhaseIPFJet_H__
#define __L1Analysis_L1AnalysisPhaseIPFJet_H__

//-------------------------------------------------------------------------------
// Created 02/03/2010 - A.C. Le Bihan
// 
// 
// Original code : UserCode/L1TriggerDPG/L1ExtraTreeProducer - Jim Brooke
//-------------------------------------------------------------------------------

#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisPhaseIPFJetDataFormat.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFJet.h"

namespace L1Analysis
{
  class L1AnalysisPhaseIPFJet 
  {
  public:
    L1AnalysisPhaseIPFJet();
    ~L1AnalysisPhaseIPFJet();
    void Reset() {l1extra_.Reset();}

    // Add new PFJet collections 
    void SetGenJet  (const edm::Handle<reco::GenJetCollection > genJets,    unsigned maxL1Extra);
    void SetPhaseIPFJet  (const edm::Handle< vector<reco::CaloJet> >  phaseIPFJets,    unsigned maxL1Extra);
    void SetPFJet  (const      edm::Handle<l1t::PFJetCollection>  PFJet,    unsigned maxL1Extra);

    L1AnalysisPhaseIPFJetDataFormat * getData() {return &l1extra_;}




  private :
    L1AnalysisPhaseIPFJetDataFormat l1extra_;
  }; 
}
#endif


