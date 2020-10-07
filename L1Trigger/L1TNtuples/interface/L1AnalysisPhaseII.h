#ifndef __L1Analysis_L1AnalysisPhaseII_H__
#define __L1Analysis_L1AnalysisPhaseII_H__

//-------------------------------------------------------------------------------
// Created 02/03/2010 - A.C. Le Bihan
// 
// 
// Original code : UserCode/L1TriggerDPG/L1ExtraTreeProducer - Jim Brooke
//-------------------------------------------------------------------------------

#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"

#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1HFRingsFwd.h"
#include "DataFormats/L1Trigger/interface/L1HFRings.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/L1TVertex/interface/Vertex.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkGlbMuonParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkGlbMuonParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkPrimaryVertex.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkHTMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkHTMissParticleFwd.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkTauParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkTauParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TrkTauParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEGTauParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1CaloTkTauParticle.h"

//#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFJet.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h" 

#include "DataFormats/L1Trigger/interface/L1PFTau.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFTau.h"

#include "DataFormats/Phase2L1Taus/interface/L1HPSPFTau.h"
#include "DataFormats/Phase2L1Taus/interface/L1HPSPFTauFwd.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkBsCandidate.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkBsCandidateFwd.h"

#include "DataFormats/L1TMuon/interface/BayesMuCorrelatorTrack.h" 

#include "DataFormats/JetReco/interface/CaloJet.h"

#include "L1Trigger/L1TNtuples/interface/L1AnalysisPhaseIIDataFormat.h"



namespace L1Analysis
{
  class L1AnalysisPhaseII 
  {
  public:
    L1AnalysisPhaseII();
    ~L1AnalysisPhaseII();
    void Reset() {l1extra_.Reset();}

    // Fill DZ of Vertex, different algorithms
    void SetVertices(float z0Puppi, float z0VertexTDR, const edm::Handle<std::vector<l1t::Vertex> > l1vertices, const edm::Handle<std::vector<l1t::L1TkPrimaryVertex> > l1TkPrimaryVertex);


    // Old style objects (Phase I)
//    void SetJet  (const edm::Handle<l1t::JetBxCollection>    jet,  unsigned maxL1Extra);
//    void SetSum  (const edm::Handle<l1t::EtSumBxCollection>  sums, unsigned maxL1Extra);
    void SetMuon (const edm::Handle<l1t::MuonBxCollection>   muon, unsigned maxL1Extra);

    // Add new standalone objects 
    void SetEG   (const edm::Handle<l1t::EGammaBxCollection> EG,  const edm::Handle<l1t::EGammaBxCollection> EGHGC, unsigned maxL1Extra);
    void SetMuonKF (const edm::Handle<l1t::RegionalMuonCandBxCollection>   muonKF, unsigned maxL1Extra, unsigned int muonDetector);
    void SetCaloJet  (const edm::Handle<l1t::JetBxCollection>    calojet,  unsigned maxL1Extra, float caloJetHTT);
    void SetCaloTau  (const edm::Handle<l1t::TauBxCollection>    calotau,  unsigned maxL1Extra);


    // Add L1TrackTriggerObjects
    void SetTkEG   (const  edm::Handle<l1t::L1TkElectronParticleCollection>   tkEG,  const  edm::Handle<l1t::L1TkElectronParticleCollection>   tkEGHGC,    unsigned maxL1Extra);
    void SetTkEGV2   (const  edm::Handle<l1t::L1TkElectronParticleCollection>   tkEGV2, const  edm::Handle<l1t::L1TkElectronParticleCollection>   tkEGV2HGc,    unsigned maxL1Extra);
    void SetTkEM   (const  edm::Handle<l1t::L1TkEmParticleCollection>   tkEM,  const  edm::Handle<l1t::L1TkEmParticleCollection>   tkEMHGC,    unsigned maxL1Extra);
    void SetTkGlbMuon (const  edm::Handle<l1t::L1TkGlbMuonParticleCollection> TkGlbMuon,   unsigned maxL1Extra);
    void SetTkMuon (const  edm::Handle<l1t::L1TkMuonParticleCollection> TkMuon,   unsigned maxL1Extra);
    void SetTkMuonStubs (const  edm::Handle<l1t::L1TkMuonParticleCollection> TkMuon,   unsigned maxL1Extra,unsigned int muonDetector);
    void SetTkMuonStubsOMTF (const  edm::Handle<l1t::BayesMuCorrTrackBxCollection> TkMuonOMTF,   unsigned maxL1Extra,unsigned int muonDetector);

    void SetHSCP (const  edm::Handle<l1t::BayesMuCorrTrackBxCollection> TkMuonHSCP,   unsigned maxL1Extra);

    void SetTrkTau  (const  edm::Handle<l1t::L1TrkTauParticleCollection> tkTau, unsigned maxL1Extra);
    void SetCaloTkTau  (const  edm::Handle<l1t::L1CaloTkTauParticleCollection> CaloTkTau, unsigned maxL1Extra);
    void SetTkEGTau  (const  edm::Handle<l1t::L1TkEGTauParticleCollection> TkEGTau, unsigned maxL1Extra);


    void SetTkJet  (const  edm::Handle<l1t::L1TkJetParticleCollection>  tkTrackerJet,    unsigned maxL1Extra);
    void SetTkCaloJet  (const  edm::Handle<l1t::L1TkJetParticleCollection>  tkCaloJet,    unsigned maxL1Extra);
    void SetTkMET  (const  edm::Handle<l1t::L1TkEtMissParticleCollection> trackerMets);
    void SetTkMHT  (const  edm::Handle<l1t::L1TkHTMissParticleCollection> trackerMHTs);

    // Add new PFJet collections 
    void SetPFJet  (const      edm::Handle<l1t::PFJetCollection>  PFJet,    unsigned maxL1Extra);
    void SetL1METPF(const edm::Handle< std::vector<reco::PFMET> > l1MetPF);
    void SetPFObjects(const edm::Handle< vector<l1t::PFCandidate> > l1pfCandidates,unsigned maxL1Extra);
//    void SetPFJetForMET  (const      edm::Handle<l1t::PFJetCollection>  PFJet,    unsigned maxL1Extra); // this needs to be done better 

    // reco::caloJet collection for "Phase1L1Jets" ...
    void SetL1PfPhase1L1TJet(const      edm::Handle< std::vector<reco::CaloJet> >  l1L1PFPhase1L1Jet,    unsigned maxL1Extra);



    // Add new PFTau
    void SetPFTaus  (const      edm::Handle< std::vector<l1t::L1PFTau>>  l1pfTaus,    unsigned maxL1Extra);
    void SetNNTaus  (const      edm::Handle< std::vector<l1t::PFTau>>  l1nnTaus,    unsigned maxL1Extra);
    void SetHPSPFTaus  (const  edm::Handle<l1t::L1HPSPFTauCollection> HPSPFTau, unsigned maxL1Extra);
    void SetNNTauPFs  (const      edm::Handle< std::vector<l1t::PFTau>>  l1nnTauPFs,    unsigned maxL1Extra);

    void SetBsCands (const      edm::Handle< std::vector<l1t::L1TkBsCandidate>>  l1TkBs,    unsigned maxL1Extra, int kind);

    L1AnalysisPhaseIIDataFormat * getData() {return &l1extra_;}

    // DiObjects
    void SetDiMuonTk(const edm::Handle<l1t::L1TkMuonParticleCollection> muon, unsigned maxL1Extra); 

    static int transverseCoord(double cxa, double cya, double ra,
                       double cxb, double cyb, double rb,
                       double & xg1, double & yg1,
                       double & xg2, double & yg2)  dso_internal;

    // Computes z-coordinate on helix at given transverse coordinates
    static double zCoord(const GlobalVector& mom, const GlobalPoint& pos,
                        double r, double xc, double yc, double xg, double yg)  dso_internal;





  private :
    L1AnalysisPhaseIIDataFormat l1extra_;
    int tk_nFitParams_=4 ; // Harcoding this, choosing 4,
                              // to not have to store the chosen fitParams for all objects in this tree producer as a configuration.
                              // (it would be cleaner if all objects save the Z directly as well as the pointer to the track, or if
                              // it is clear that the default is 4 unless specifically stated)
  }; 
}
#endif


