import FWCore.ParameterSet.Config as cms

TwoLayerJets = cms.EDProducer('TwoLayerJets',
	        L1TrackInputTag= cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),	
		L1Tk_nPar = cms.int32(5),         # use 4 or 5-parameter L1 track fit ??
		ZMAX = cms.double ( 15. ) ,
		PTMAX = cms.double( 200. ), # Maximumum track pT before saturation
		Etabins=cms.int32(24),
		Phibins=cms.int32(27),
		Zbins=cms.int32(60),
	 	TRK_PTMIN = cms.double(2.0),        # minimum track pt [GeV]
   	 	TRK_ETAMAX = cms.double(2.4),       # maximum track eta
		CHI2_MAX=cms.double(50.),	
		PromptBendConsistency=cms.double(1.75), #Bend Chi2 Cut for prompt tracks
		D0_CutNstubs4=cms.double(0.15),  #Flag for displaced tracks
		D0_CutNstubs5=cms.double(1.0),  #Flag for displaced tracks
		nPSStubsMin=cms.int32(-1),
		minTrkJetpT=cms.double(5.),
		LowpTJetMinTrackMultiplicity=cms.int32(2),
		HighpTJetMinTrackMultiplicity=cms.int32(3),	
		DisplacedAlgo=cms.bool(False),
		NStubs4DisplacedChi2_Loose=cms.double(5.0), ########Displaced track quality flags for loose/tight
		NStubs4Displacedbend_Loose=cms.double(3.0),	
		NStubs5DisplacedChi2_Loose=cms.double(3.0), ########Displaced track quality flags for loose/tight
		NStubs5Displacedbend_Loose=cms.double(3.0),	
		NStubs4DisplacedChi2_Tight=cms.double(3.0),
		NStubs4Displacedbend_Tight=cms.double(3.0),
		NStubs5DisplacedChi2_Tight=cms.double(2.5),
		NStubs5Displacedbend_Tight=cms.double(3.0)
)
