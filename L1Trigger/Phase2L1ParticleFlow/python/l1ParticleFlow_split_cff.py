import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.l1EGammaCrystalsProducer_cfi import l1EGammaCrystalsProducer
from L1Trigger.Phase2L1ParticleFlow.hgc3dClustersForPF_cff import *

l1ParticleFlow_prerequisites = cms.Sequence(
    l1EGammaCrystalsProducer + 
    hgc3dClustersForPF_STC
)

from L1Trigger.Phase2L1ParticleFlow.pfTracksFromL1Tracks_cfi import pfTracksFromL1Tracks
import L1Trigger.Phase2L1ParticleFlow.pfClustersFromHGC3DClusters_cfi
from L1Trigger.Phase2L1ParticleFlow.pfClustersFromL1EGClusters_cfi import pfClustersFromL1EGClusters
from L1Trigger.Phase2L1ParticleFlow.pfClustersFromCombinedCalo_cfi import pfClustersFromCombinedCalo
from L1Trigger.Phase2L1ParticleFlow.l1pfProducer_cfi import l1pfProducer

# Calorimeter part: ecal + hcal + hf only
pfClustersFromCombinedCaloHCal = pfClustersFromCombinedCalo.clone(
    hcalHGCTowers = [], 
    hcalDigisBarrel = True, hcalDigisHF = False,
    hadCorrector = cms.string("L1Trigger/Phase2L1ParticleFlow/data/hadcorr_barrel.root"),
    )
pfClustersFromCombinedCaloHF = pfClustersFromCombinedCalo.clone(
    ecalCandidates = [], hcalHGCTowers = [],
    hcalDigisBarrel = False, hcalDigisHF = True,
    hadCorrector = cms.string("L1Trigger/Phase2L1ParticleFlow/data/hfcorr.root"),
    resol = cms.PSet(
            etaBins = cms.vdouble( 3.500,  4.000,  4.500,  5.000),
            offset  = cms.vdouble( 0.580,  0.594,  0.355, -1.310),
            scale   = cms.vdouble( 0.174,  0.172,  0.220,  0.692),
            ptMin   = cms.vdouble( 5.000,  5.000,  5.000,  5.000),
            ptMax   = cms.vdouble(999999, 999999, 999999, 999999),
            kind    = cms.string('calo'),
    ))



# Calorimeter part: hgcal
hgc3DClustersNoNoise = cms.EDProducer("HGC3DClusterSimpleSelector",
    src = cms.InputTag("hgcalBackEndLayer2ProducerSTC","HGCalBackendLayer2Processor3DClustering"),
    cut = cms.string("coreShowerLength>3"),
    )
pfClustersFromHGC3DClusters = L1Trigger.Phase2L1ParticleFlow.pfClustersFromHGC3DClusters_cfi.pfClustersFromHGC3DClusters.clone(
    src = cms.InputTag("hgc3DClustersNoNoise"),
    corrector = cms.string("L1Trigger/Phase2L1ParticleFlow/data/hadcorr_HGCal3D_STC.root"),
    correctorEmfMax = cms.double(1.125),
    emId  = cms.string("hOverE < 0.25 && hOverE >= 0"),
    etMin = 1.0, 
    resol = cms.PSet(
        etaBins = cms.vdouble( 1.900,  2.200,  2.500,  2.800,  3.100),
        offset  = cms.vdouble( 0.731,  1.632,  0.282,  0.926,  1.122),
        scale   = cms.vdouble( 0.114,  0.065,  0.098,  0.075,  0.103),
        kind    = cms.string('calo')
    ),
)

l1ParticleFlow_calo = cms.Sequence(
    pfClustersFromL1EGClusters +
    pfClustersFromCombinedCaloHCal +
    pfClustersFromCombinedCaloHF +
    hgc3DClustersNoNoise +
    pfClustersFromHGC3DClusters
)


# PF in the barrel
l1pfProducerBarrel = l1pfProducer.clone(
    # inputs
    emClusters = [ cms.InputTag('pfClustersFromL1EGClusters') ],
    hadClusters = [ cms.InputTag('pfClustersFromCombinedCaloHCal:calibrated') ],
    # track-based PUPPI
    puppiUsingBareTracks = True,
    puppiDrMin = 0.07,
    puppiPtMax = 50.,
    vtxAlgo = "external",
    vtxCollection = cms.InputTag("VertexProducer","l1vertices"),
    # regionalize
    useRelativeRegionalCoordinates = cms.bool(False),
    trackRegionMode = cms.string("atCalo"),
    regions = cms.VPSet(
        cms.PSet(
            etaBoundaries = cms.vdouble(-1.5,1.5),
            phiSlices = cms.uint32(1),
            etaExtra = cms.double(0.3),
            phiExtra = cms.double(0.0)
        ),
    ),
)
l1ParticleFlow_pf_barrel = cms.Sequence(
    pfTracksFromL1Tracks +   
    l1pfProducerBarrel
)



# PF in HGCal
pfTracksFromL1TracksHGCal = pfTracksFromL1Tracks.clone(
    resolCalo = pfClustersFromHGC3DClusters.resol.clone(),
)
l1pfProducerHGCal = l1pfProducer.clone(
    # algo
    pfAlgo = "PFAlgo2HGC",
    # inputs
    tracks = cms.InputTag('pfTracksFromL1TracksHGCal'),
    emClusters  = [ ],  # EM clusters are not used (only added to NTuple for calibration/monitoring)
    hadClusters = [ cms.InputTag("pfClustersFromHGC3DClusters") ],
    # track-based PUPPI
    puppiDrMin = 0.04,
    puppiPtMax = 50.,
    puppiUsingBareTracks = True,
    vtxAlgo = "external",
    vtxCollection = cms.InputTag("VertexProducer","l1vertices"),
    # regionalize
    useRelativeRegionalCoordinates = cms.bool(False),
    trackRegionMode = cms.string("atCalo"),
    regions = cms.VPSet(
        cms.PSet(
            etaBoundaries = cms.vdouble(-3,-1.5),
            phiSlices = cms.uint32(1),
            etaExtra = cms.double(0.3),
            phiExtra = cms.double(0.0)
        ),
        cms.PSet(
            etaBoundaries = cms.vdouble(1.5,3.0),
            phiSlices = cms.uint32(1),
            etaExtra = cms.double(0.3),
            phiExtra = cms.double(0.0)
        ),
    ),
)

l1ParticleFlow_pf_hgcal = cms.Sequence(
    pfTracksFromL1TracksHGCal +   
    l1pfProducerHGCal
)



# PF in HF
l1pfProducerHF = l1pfProducer.clone(
    # inputs
    emClusters = [ ],
    hadClusters = [ cms.InputTag('pfClustersFromCombinedCaloHF:calibrated') ],
    # not really useful, but for consistency
    puppiDrMin = 0.1,
    puppiPtMax = 100.,
    vtxAlgo = "external",
    vtxCollection = cms.InputTag("VertexProducer","l1vertices"),
    # regionalize
    useRelativeRegionalCoordinates = cms.bool(False),
    trackRegionMode = cms.string("atCalo"),
    regions = cms.VPSet(
        cms.PSet(
            etaBoundaries = cms.vdouble(-5.5,-3),
            phiSlices = cms.uint32(1),
            etaExtra = cms.double(0.0),
            phiExtra = cms.double(0.0)
        ),
        cms.PSet(
            etaBoundaries = cms.vdouble(3,5.5),
            phiSlices = cms.uint32(1),
            etaExtra = cms.double(0.0),
            phiExtra = cms.double(0.0)
        ),
    )
)
l1ParticleFlow_pf_hf = cms.Sequence(
    l1pfProducerHF
)


# Merging all outputs
l1pfCandidates = cms.EDProducer("L1TPFCandMultiMerger",
    pfProducers = cms.VInputTag(
        cms.InputTag("l1pfProducerBarrel"), 
        cms.InputTag("l1pfProducerHGCal"),
        cms.InputTag("l1pfProducerHF")
    ),
    labelsToMerge = cms.vstring("Calo", "TK", "TKVtx", "PF", "Puppi"),
)

l1ParticleFlow_proper = cms.Sequence(
    l1ParticleFlow_calo +
    l1ParticleFlow_pf_barrel +
    l1ParticleFlow_pf_hgcal +
    l1ParticleFlow_pf_hf +
    l1pfCandidates
)

l1ParticleFlow = cms.Sequence(l1ParticleFlow_prerequisites + l1ParticleFlow_proper)
