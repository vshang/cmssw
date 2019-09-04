import FWCore.ParameterSet.Config as cms
process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')   
process.load("Configuration.StandardSequences.MagneticField_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/scratch3/MTF/data/190401/singleMu0.root'
    ),
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '90X_upgrade2023_realistic_v9', '')

#Some stuff that we do not care###################
###############################################
process.load('L1Trigger.L1TMuonBarrel.fakeBmtfParams_cff')
process.esProd = cms.EDAnalyzer("EventSetupRecordDataGetter",
   toGet = cms.VPSet(
      cms.PSet(record = cms.string('L1TMuonBarrelParamsRcd'),
               data = cms.vstring('L1TMuonBarrelParams'))
                   ),
   verbose = cms.untracked.bool(True)
)
process.fakeBmtfParams.fwVersion = cms.uint32(2)
process.fakeBmtfParams.BX_max = cms.int32(2)
process.fakeBmtfParams.BX_min = cms.int32(-2)
maskenable      = '000000000000'
maskdisable     = '111111111111'
process.fakeBmtfParams.mask_phtf_st1        = cms.vstring(maskenable, maskenable, maskenable, maskenable, maskenable, maskenable, maskenable)
process.fakeBmtfParams.mask_phtf_st2        = cms.vstring(maskenable,  maskenable, maskenable, maskenable, maskenable, maskenable, maskenable)
process.fakeBmtfParams.mask_phtf_st3        = cms.vstring(maskenable,  maskenable, maskenable, maskenable, maskenable, maskenable, maskenable)
process.fakeBmtfParams.mask_phtf_st4        = cms.vstring(maskenable,  maskenable, maskenable, maskenable, maskenable, maskenable, maskenable)
process.fakeBmtfParams.mask_ettf_st1        = cms.vstring(maskenable, maskenable, maskenable, maskenable, maskenable, maskenable, maskenable)
process.fakeBmtfParams.mask_ettf_st2        = cms.vstring(maskenable,  maskenable, maskenable, maskenable, maskenable, maskenable, maskenable)
process.fakeBmtfParams.mask_ettf_st3        = cms.vstring(maskenable,  maskenable, maskenable, maskenable, maskenable, maskenable, maskenable)
##################################################
################################################


# Load Tracker Plus Stubs
process.load('L1Trigger.L1TTrackMatch.L1TTrackerPlusStubs_cfi')
#print('Loaded tracker plus stubs')


#Configure for Phase 2



#print('Configuring for phase 2')
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),

)

  
#Path that makes stubs and then runs Tracker+Stubs
#print('Running tracker+stubs')
process.p = cms.Path(process.l1TPSStubs)
process.e = cms.EndPath(process.out)
process.schedule = cms.Schedule(process.p,process.e)

#print('End of python script')
