import FWCore.ParameterSet.Config as cms

from L1TriggerConfig.DTTPGConfigProducers.L1DTTPGConfigFromDB_cff import *
from L1Trigger.DTPhase2Trigger.HoughGrouping_cfi                  import HoughGrouping
from L1Trigger.DTPhase2Trigger.PseudoBayesGrouping_cfi            import PseudoBayesPattern

dtTriggerPhase2PrimitiveDigis = cms.EDProducer("DTTrigPhase2Prod",
                                               digiTag = cms.InputTag("CalibratedDigis"),
                                               trigger_with_sl = cms.untracked.int32(4),
                                               tanPhiTh = cms.untracked.double(1.),
                                               chi2Th = cms.untracked.double(0.01), #in cm^2
                                               chi2corTh = cms.untracked.double(0.1), #in cm^2
                                               do_correlation = cms.untracked.bool(True),
                                               dT0_correlate_TP = cms.untracked.double(25.),
                                               #minx_match_2digis = cms.untracked.double(2.1),
                                               minx_match_2digis = cms.untracked.double(1.),
                                               p2_df = cms.untracked.int32(2), #0 for phase-1, 1 for slice-test, 2 for phase-2 carlo-federica
                                               scenario = cms.untracked.int32(1), #0 for mc, 1 for data, 2 for slice test
                                               filter_cousins = cms.untracked.bool(True),
                                               apply_txt_ttrig_bc0 = cms.untracked.bool(False),

                                               ttrig_filename = cms.FileInPath('L1Trigger/DTPhase2Trigger/data/wire_rawId_ttrig.txt'),
                                               z_filename = cms.FileInPath('L1Trigger/DTPhase2Trigger/data/wire_rawId_z.txt'),
                                               shift_filename = cms.FileInPath('L1Trigger/DTPhase2Trigger/data/wire_rawId_x.txt'),
                                               grouping_code = cms.untracked.int32(0),       # 0 = initial grouping, 1 = Hough transform, 2 = PseudoBayes Approach
                                               min_phinhits_match_segment = cms.untracked.int32(8),
                                               min_dT0_match_segment = cms.untracked.double(12.5),
                                               minHits4Fit = cms.untracked.int32(4),
                                               #debugging
                                               debug = cms.untracked.bool(False),
                                               dump  = cms.untracked.bool(False),
                                               #RPC
                                               rpcRecHits = cms.untracked.InputTag("rpcRecHits"),
                                               useRPC = cms.untracked.bool(False)
                                               )

dtTriggerPhase2PrimitiveDigis.HoughGrouping      = HoughGrouping
dtTriggerPhase2PrimitiveDigis.PseudoBayesPattern = PseudoBayesPattern
