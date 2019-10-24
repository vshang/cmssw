#ifndef __L1Analysis_L1AnalysisPhaseIPFJetDataFormat_H__
#define __L1Analysis_L1AnalysisPhaseIPFJetDataFormat_H__

//-------------------------------------------------------------------------------
// Created 20/04/2010 - E. Conte, A.C. Le Bihan
// 
// 
// Original code : UserCode/L1TriggerDPG/L1ExtraTreeProducer - Jim Brooke
//-------------------------------------------------------------------------------


#include <vector>

namespace L1Analysis
{
  struct L1AnalysisPhaseIPFJetDataFormat
  {
    L1AnalysisPhaseIPFJetDataFormat(){Reset();};
    ~L1AnalysisPhaseIPFJetDataFormat(){};
    
    void Reset()
    {

      nPhaseIPFJets = 0;
      phaseIPFJetEt.clear();
      phaseIPFJetEta.clear();
      phaseIPFJetPhi.clear();
      phaseIPFJetHt = 0;

      nAK4PFJets = 0;
      ak4PFJetEt.clear();
      ak4PFJetEta.clear();
      ak4PFJetPhi.clear();
      ak4PFJetHt = 0;

      nGenJet = 0;     
      genJetPt.clear();
      genJetEta.clear();
      genJetPhi.clear();
      genJetM.clear();
      genJetHt = 0;


    }
 
    unsigned short int nPhaseIPFJets;
    std::vector<double> phaseIPFJetEt;
    std::vector<double> phaseIPFJetEta;
    std::vector<double> phaseIPFJetPhi;
    double phaseIPFJetHt;

    unsigned short int nAK4PFJets;
    std::vector<double> ak4PFJetEt;
    std::vector<double> ak4PFJetEta;
    std::vector<double> ak4PFJetPhi;
    double ak4PFJetHt;

    int nGenJet;
    std::vector<float> genJetPt;
    std::vector<float> genJetEta;
    std::vector<float> genJetPhi;
    std::vector<float> genJetM;
    double genJetHt;




  }; 
}
#endif


