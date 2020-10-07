import FWCore.ParameterSet.Config as cms

L1TkCaloJets = cms.EDProducer("L1TkCaloJetProducer",
    L1CentralJetInputTag = cms.InputTag("simCaloStage2Digis",""),    # new L1Jets
    L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
    TRK_ZMAX = cms.double(25.),         # max track z0 [cm]
    TRK_CHI2MAX = cms.double(100.),     # max track chi2
    TRK_PTMIN = cms.double(2.0),        # minimum track pt [GeV]
    TRK_ETAMAX = cms.double(2.5),       # maximum track eta
    TRK_NSTUBMIN = cms.int32(4),        # minimum number of stubs on track
    TRK_NSTUBPSMIN = cms.int32(2),      # minimum number of stubs in PS modules on track
    doPtComp = cms.bool( False ),	   # track-stubs PT compatibility cut
    doTightChi2 = cms.bool( False )    # chi2dof < 5 for tracks with PT > 10
)
