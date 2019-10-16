import FWCore.ParameterSet.Config as cms

from L1Trigger.Phase2L1ParticleFlow.l1pfJetMet_cff import ak4PFL1PF, ak4PFL1PFCorrected
produceL1HPSPFTausPF = cms.Sequence(ak4PFL1PF + ak4PFL1PFCorrected)

l1pfNeutralCandidatesPF = cms.EDFilter("L1TPFCandSelector",
    src = cms.InputTag('l1pfCandidates:PF'),                                      
    cut = cms.string("pdgId = 22"),
    filter = cms.bool(False)                                           
)
produceL1HPSPFTausPF += l1pfNeutralCandidatesPF

kt6L1PFJetsNeutralsPF = ak4PFL1PF.clone(
    src          = cms.InputTag('l1pfNeutralCandidatesPF'), 
    jetAlgorithm = cms.string("Kt"),
    rParam       = cms.double(0.6),
    doRhoFastjet = cms.bool(True),
    Rho_EtaMax   = cms.double(3.0)
)
produceL1HPSPFTausPF += kt6L1PFJetsNeutralsPF

from L1Trigger.Phase2L1Taus.L1HPSPFTauProducerPF_cfi import L1HPSPFTauProducerPF
produceL1HPSPFTausPF += L1HPSPFTauProducerPF
