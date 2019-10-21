import FWCore.ParameterSet.Config as cms

from L1TriggerConfig.DTTPGConfigProducers.L1DTTPGConfigFromDB_cff import *
from L1Trigger.DTPhase2Trigger.HoughGrouping_cfi                  import HoughGrouping
from L1Trigger.DTPhase2Trigger.PseudoBayesGrouping_cfi            import PseudoBayesPattern

dtTriggerPhase2PrimitiveDigis = cms.EDProducer("DTTrigPhase2Prod",
                                               digiTag = cms.InputTag("CalibratedDigis"),
                                               trigger_with_sl = cms.untracked.int32(4),
                                               use_normal_chi2 = cms.untracked.bool(True),
                                               tanPhiTh = cms.untracked.double(999.), #Removing TanPsi filter in SliceTest, as FW doesnt use it
                                               #tanPhiTh = cms.untracked.double(1.),
                                               chi2Th = cms.untracked.double(0.03), #in cm^2 SliceTest JM 
                                               #chi2Th = cms.untracked.double(0.01), #in cm^2
                                               chi2corTh = cms.untracked.double(0.1), #in cm^2
                                               do_correlation = cms.untracked.bool(True),
                                               useBX_correlation = cms.untracked.bool(True),
                                               dT0_correlate_TP = cms.untracked.double(25.), 
                                               dBX_correlate_TP = cms.untracked.int32(0), 
                                               dTanPsi_correlate_TP = cms.untracked.double(620./4096.),
					       clean_chi2_correlation = cms.untracked.bool(True),
					       use_LSB = cms.untracked.bool(False),
					       tanPsi_precision = cms.untracked.double(1./4096.),
					       x_precision = cms.untracked.double(0.025),
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
                                               #Print hits and prims for tecno30
					       printPython  = cms.untracked.bool(False),
                                               printHits  = cms.untracked.bool(False),
                                               #RPC
                                               rpcRecHits = cms.untracked.InputTag("rpcRecHits"),
                                               useRPC = cms.untracked.bool(False),
                                               bx_window = cms.untracked.int32(1), # will look for RPC cluster within a bunch crossing window of 'dt.BX +- bx_window' 
                                               phi_window = cms.untracked.double(50.), # will look for RPC cluster within a phi window of +- phi_window in arbitrary coordinates (plot the value we cut on in RPCIntergator to fine tune it)
                                               max_quality_to_overwrite_t0 = cms.untracked.int32(9), # will use RPC  to set 't0' for TP with quality < 'max_quality_to_overwrite_t0'
                                               storeAllRPCHits = cms.untracked.bool(False)
                                               )

dtTriggerPhase2PrimitiveDigis.HoughGrouping      = HoughGrouping
dtTriggerPhase2PrimitiveDigis.PseudoBayesPattern = PseudoBayesPattern
