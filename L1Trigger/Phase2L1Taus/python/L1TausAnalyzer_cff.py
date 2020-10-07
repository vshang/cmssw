import FWCore.ParameterSet.Config as cms

from L1Trigger.Phase2L1Taus.L1TausAnalyzer_cfi import l1TausAnalysis
# TrkTaus
TrkTauEff                  = l1TausAnalysis.clone()
TrkTauEff.ObjectType       = cms.string("TrkTau")
TrkTauEff.GenEtaVisCutOff  = cms.double(2.3)
TrkTauEff.L1EtaCutOff      = cms.double(2.5)
TrkTauEff.AnalysisOption   = cms.string("Efficiency")
                           
TrkTauRate                 = l1TausAnalysis.clone()
TrkTauRate.ObjectType      = cms.string("TrkTau")
TrkTauRate.GenEtaVisCutOff = cms.double(2.3)
TrkTauRate.L1EtaCutOff     = cms.double(2.5)
TrkTauRate.AnalysisOption  = cms.string("Rate")
                           
# TkEGTaus                 
TkEGEff                    = l1TausAnalysis.clone()
TkEGEff.ObjectType         = cms.string("TkEG")
TkEGEff.AnalysisOption     = cms.string("Efficiency")
                           
TkEGRate                   = l1TausAnalysis.clone()
TkEGRate.ObjectType        = cms.string("TkEG")
TkEGRate.AnalysisOption    = cms.string("Rate")
                           
# CaloTkTaus               
CaloTkEff                  = l1TausAnalysis.clone()
CaloTkEff.ObjectType       = cms.string("CaloTk")
CaloTkEff.AnalysisOption   = cms.string("Efficiency")
CaloTkEff.DRMatching       = cms.double(0.2)

CaloTkRate                 = l1TausAnalysis.clone()
CaloTkRate.ObjectType      = cms.string("CaloTk")
CaloTkRate.AnalysisOption  = cms.string("Rate")
