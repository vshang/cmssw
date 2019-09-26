import FWCore.ParameterSet.Config as cms
import math
from math import pi


l1TPSStubs = cms.EDProducer("L1TMuCorrelatorHitProducer",
    verbose = cms.int32(0),
    srcCSC = cms.InputTag("simCscTriggerPrimitiveDigis"),
    srcDT = cms.InputTag("simDtTriggerPrimitiveDigis"),
    srcRPC = cms.InputTag("simMuonRPCDigis"),
    CSC =cms.PSet(                            
        verbose = cms.int32(0),
        minBX              = cms.int32(0),                           
        maxBX              = cms.int32(0),         
        phiLSB             = cms.double(0.00019174760), 
        etaLSB             = cms.double(0.00585938), 
        ),
    DT = cms.PSet(                         
        verbose            = cms.int32(0),
        minPhiQuality      = cms.int32(0),
        minThetaQuality    = cms.int32(0),
        minBX              = cms.int32(0),                           
        maxBX              = cms.int32(0),                           
        disableMasks       = cms.bool(True),                               
        phiLSB             = cms.double(0.00019174760),
        eta_1         = cms.vint32(197,189,182,173,165,156,147,129,120,108,97,86,74,63,38,26,12,0,-12,-25,-38,-62,-74,-86,-98,-108,-120,-129,-148,-157,-165,-172,-182,-189,-197),
        eta_2         = cms.vint32(174,167,159,153,144,136,128,112,103,93,83,73,63,53,33,21,10,0,-10,-21,-32,-54,-63,-73,-84,-93,-102,-111,-128,-136,-144,-153,-160,-167,-175),
        eta_3         = cms.vint32(0,141,137,131,123,117,109,94,86,78,70,61,53,44,27,18,8,0,-8,-18,-26,-45,-53,-61,-70,-77,-86,-94,-109,-116,-123,-131,-137,-141,0),
        coarseEta_1 = cms.vint32(0,-102,-171),
        coarseEta_2 = cms.vint32(0,-85,-154),
        coarseEta_3 = cms.vint32(0,-68,-137),
        coarseEta_4 = cms.vint32(0,-68,-119),
        bendingScale       = cms.vdouble(16.0*1.7370158,16.0*2.6109661,16.0*4.4111160,16.0*11.660448),
   ),
   RPCBarrel = cms.PSet(                         
        verbose            = cms.int32(0),
        minPhiQuality      = cms.int32(0),
        minThetaQuality    = cms.int32(0),
        minBX              = cms.int32(0),                           
        maxBX              = cms.int32(0),                           
        disableMasks       = cms.bool(False),                               
        phiLSB             = cms.double(0.00019174760),
        etaLSB             = cms.double(0.00585938), 
        coarseEta_1 = cms.vint32(0,-102,-171),
        coarseEta_2 = cms.vint32(0,-85,-154),
        coarseEta_3 = cms.vint32(0,-68,-137),
        coarseEta_4 = cms.vint32(0,-68,-119)
   ),
   RPCEndcap = cms.PSet(                         
        verbose            = cms.int32(0),
        minBX              = cms.int32(0),                           
        maxBX              = cms.int32(0),                           
        phiLSB             = cms.double(0.00019174760),
        etaLSB             = cms.double(0.00585938), 
   )


)

tpsAlgoSettings = cms.PSet(
    verbose = cms.int32(0),
    phiLSB = cms.double(0.00019174760),
    etaLSB = cms.double(0.000366211),
    etaShift = cms.uint32(4),
    curvLSB  = cms.double(0.000122070),
    vetoIndex  =cms.vuint32(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,5,5,5,5,6),
    vetoPattern=cms.vuint32(3,9,10,18,120,320,48,192,3,5,6,9,10,11,12,18,44,48,88,864,1408,384,448,197,3,10,3,3,6,9,10,5)
)
        

        



l1StubMatchedMuons = cms.EDProducer("L1TTPSProducer",
    srcStubs = cms.InputTag("l1TPSStubs"),
    srcTracks = cms.InputTag("TTTracksFromTracklet:Level1TTTracks"),
    maxChi2 = cms.double(100000000),
    minStubs = cms.uint32(4),
    verbose = cms.int32(0),                                    
    sectors = cms.VPSet(
        cms.PSet(
            sectorNumber = cms.uint32(0),
            barrelSectors= cms.vuint32(0,1,2,3,4,5,6,7,8,9,10,11),
            csc10DegreeChambers=cms.vuint32(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36),
            csc20DegreeChambers=cms.vuint32(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),
            rpcEndcapChambers=cms.vuint32(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36),
            iRpcChambers=cms.vuint32(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36),
            phiLowerBound = cms.int32(-16384),
            phiUpperBound = cms.int32(16384),
            phiOffset=cms.int32(0),
            algoSettings = tpsAlgoSettings
        )
    )
)

l1TrackerPlusStubsSequence = cms.Sequence(l1TPSStubs*l1StubMatchedMuons)
