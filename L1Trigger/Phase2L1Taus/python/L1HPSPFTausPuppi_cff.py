import FWCore.ParameterSet.Config as cms

from L1Trigger.Phase2L1ParticleFlow.l1pfJetMet_cff import ak4PFL1Puppi, ak4PFL1PuppiCorrected
produceL1HPSPFTausPuppi = cms.Sequence(ak4PFL1Puppi + ak4PFL1PuppiCorrected)

l1pfNeutralCandidatesPuppi = cms.EDFilter("L1TPFCandSelector",
    src = cms.InputTag('l1pfCandidates:Puppi'),                                      
    cut = cms.string("pdgId = 22"),
    filter = cms.bool(False)                                           
)
produceL1HPSPFTausPuppi += l1pfNeutralCandidatesPuppi

from L1Trigger.Phase2L1ParticleFlow.l1pfJetMet_cff import ak4PFL1PF
kt6L1PFJetsNeutralsPuppi = ak4PFL1PF.clone(
    src          = cms.InputTag('l1pfNeutralCandidatesPuppi'), 
    jetAlgorithm = cms.string("Kt"),
    rParam       = cms.double(0.6),
    doRhoFastjet = cms.bool(True),
    Rho_EtaMax   = cms.double(3.0)
)
produceL1HPSPFTausPuppi += kt6L1PFJetsNeutralsPuppi

from L1Trigger.Phase2L1Taus.L1HPSPFTauProducerPuppi_cfi import L1HPSPFTauProducerPuppi
produceL1HPSPFTausPuppi += L1HPSPFTauProducerPuppi
