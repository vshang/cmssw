import FWCore.ParameterSet.Config as cms

# Defines the L1 Emulator sequence for simulation use-case subsystem emulators
# run on the results of previous (in the hardware chain) subsystem emulator:
#
#     SimL1Emulator = cms.Sequence(...)
#
# properly configured for the current Era (e.g. Run1, 2015, or 2016).  Also
# configures event setup producers appropriate to the current Era, to handle
# conditions which are not yet available in the GT.
#
# Author List
# Jim Brooke, 24 April 2008
# Vasile Mihai Ghete, 2009
# Jim Brooke, Michael Mulhearn, 2015
# Vladimir Rekovic 2016,2017

# Notes on Inputs:

# ECAL TPG emulator and HCAL TPG run in the simulation sequence in order to be able
# to use unsuppressed digis produced by ECAL and HCAL simulation, respectively
# in Configuration/StandardSequences/python/Digi_cff.py
# SimCalorimetry.Configuration.SimCalorimetry_cff
# which calls
# SimCalorimetry.Configuration.ecalDigiSequence_cff
# SimCalorimetry.Configuration.hcalDigiSequence_cff

#
# At the moment, there is no emulator available for upgrade HF Trigger Primitives,
# so these missing (required!) inputs are presently ignored by downstream modules.
#

from L1Trigger.Configuration.SimL1TechnicalTriggers_cff import *

from L1Trigger.L1TCalorimeter.simDigis_cff import *
from L1Trigger.L1TMuon.simDigis_cff import *
from L1Trigger.L1TGlobal.simDigis_cff import *

# define a core which can be extented in customizations:
SimL1EmulatorCore = cms.Sequence(
    SimL1TCalorimeter +
    SimL1TMuon +
    SimL1TechnicalTriggers +
    SimL1TGlobal
    )

SimL1Emulator = cms.Sequence( SimL1EmulatorCore )

#
# Emulators are configured from DB (GlobalTags)
#

from L1Trigger.L1TGlobal.GlobalParameters_cff import *

# 2017 EMTF and TwinMux emulator use payloads from DB, not yet in GT,
# soon to be removed when availble in GTs
from L1Trigger.L1TTwinMux.fakeTwinMuxParams_cff import *

_phase2_siml1emulator = SimL1EmulatorTask.copy()

# ########################################################################
# ########################################################################
#
# Phase-2 
#
# ########################################################################
# ########################################################################

# ########################################################################
# Phase-2 Trigger Primitives
# ########################################################################

# HGCAL TP 
# ########################################################################
from  L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff import *
_phase2_siml1emulator.add(hgcalTriggerPrimitivesTask)
 
# ########################################################################
# Phase 2 L1T
# ########################################################################

# Barrel and EndCap EGamma
# ########################################################################

from L1Trigger.L1CaloTrigger.L1EGammaCrystalsEmulatorProducer_cfi import *
_phase2_siml1emulator.add(L1EGammaClusterEmuProducer)

from L1Trigger.L1CaloTrigger.l1EGammaEEProducer_cfi import *
_phase2_siml1emulator.add(l1EGammaEEProducer)

# ########################################################################
# Phase-2 L1T - TrackTrigger dependent modules
# ########################################################################

# Tk + StandaloneObj, including L1TkPrimaryVertex
# ########################################################################
from L1Trigger.L1TTrackMatch.L1TkObjectProducers_cff import *

_phase2_siml1emulator.add(L1TkPrimaryVertex)

_phase2_siml1emulator.add(L1TkElectronsCrystal)
_phase2_siml1emulator.add(L1TkElectronsLooseCrystal)
_phase2_siml1emulator.add(L1TkElectronsEllipticMatchCrystal)
_phase2_siml1emulator.add(L1TkIsoElectronsCrystal)
_phase2_siml1emulator.add(L1TkPhotonsCrystal)

_phase2_siml1emulator.add(L1TkElectronsHGC)
_phase2_siml1emulator.add(L1TkElectronsEllipticMatchHGC)
_phase2_siml1emulator.add(L1TkIsoElectronsHGC)
_phase2_siml1emulator.add(L1TkPhotonsHGC)

_phase2_siml1emulator.add( L1TkMuons )

# PF Candidates
# ########################################################################
from L1Trigger.Phase2L1ParticleFlow.l1ParticleFlow_cff import *
_phase2_siml1emulator.add(l1ParticleFlowTask)

# PF JetMET
# ########################################################################
from L1Trigger.Phase2L1ParticleFlow.l1pfJetMet_cff import *
# Describe here l1PFJets Task
# ###############################
l1PFJetsTask = cms.Task(
  ak4PFL1Calo , ak4PFL1PF , ak4PFL1Puppi ,
  ak4PFL1CaloCorrected , ak4PFL1PFCorrected , ak4PFL1PuppiCorrected)
_phase2_siml1emulator.add(l1PFJetsTask)
# Describe here l1PFMets Task
# ###############################
l1PFMetsTask = cms.Task(l1PFMetCalo , l1PFMetPF , l1PFMetPuppi)
_phase2_siml1emulator.add(l1PFMetsTask)

# --> add modules
#%% # Barrel EGamma
#%% # ########################################################################
from L1Trigger.L1CaloTrigger.L1EGammaCrystalsEmulatorProducer_cfi import *
_phase2_siml1emulator.add(L1EGammaClusterEmuProducer)

from L1Trigger.L1CaloTrigger.l1EGammaEEProducer_cfi import *
_phase2_siml1emulator.add(l1EGammaEEProducer)

# Tk + StandaloneObj, including L1TkPrimaryVertex
# ########################################################################
from L1Trigger.L1TTrackMatch.L1TkObjectProducers_cff import *

_phase2_siml1emulator.add(L1TkPrimaryVertex)

_phase2_siml1emulator.add(L1TkElectronsCrystal)
_phase2_siml1emulator.add(L1TkElectronsLooseCrystal)
_phase2_siml1emulator.add(L1TkElectronsEllipticMatchCrystal)
_phase2_siml1emulator.add(L1TkIsoElectronsCrystal)
_phase2_siml1emulator.add(L1TkPhotonsCrystal)

_phase2_siml1emulator.add(L1TkElectronsHGC)
_phase2_siml1emulator.add(L1TkElectronsEllipticMatchHGC)
_phase2_siml1emulator.add(L1TkIsoElectronsHGC)
_phase2_siml1emulator.add(L1TkPhotonsHGC)

_phase2_siml1emulator.add( L1TkMuons )


# PFTaus(HPS)
# ########################################################################
from L1Trigger.L1CaloTrigger.Phase1L1TJetProducer_cfi import Phase1L1TJetProducer
l1pfPhase1L1TJetProducer = Phase1L1TJetProducer.clone()
phase2_SimL1Emulator += l1pfPhase1L1TJetProducer

# PFTaus(HPS)
# ########################################################################
from L1Trigger.Phase2L1Taus.L1PFTauProducer_cff import L1PFTauProducer
l1pfTauProducer = L1PFTauProducer.clone()
l1pfTauProducer.L1PFObjects = cms.InputTag("l1pfCandidates","PF")
l1pfTauProducer.L1Neutrals = cms.InputTag("l1pfCandidates")
phase2_SimL1Emulator += l1pfTauProducer

from L1Trigger.Phase2L1Taus.L1HPSPFTausPF_cff import *
phase2_SimL1Emulator += produceL1HPSPFTausPF

from L1Trigger.Phase2L1Taus.L1HPSPFTausPuppi_cff import *
phase2_SimL1Emulator += produceL1HPSPFTausPuppi

# NNTaus
# ########################################################################
from L1Trigger.Phase2L1Taus.L1NNTauProducer_cff import *
l1NNTauProducer = L1NNTauProducer.clone()
l1NNTauProducer.L1PFObjects = cms.InputTag("l1pfCandidates","PF")
l1NNTauProducerPuppi = L1NNTauProducerPuppi.clone()
l1NNTauProducerPuppi.L1PFObjects = cms.InputTag("l1pfCandidates","Puppi")
phase2_SimL1Emulator += l1NNTauProducer
phase2_SimL1Emulator += l1NNTauProducerPuppi

# NNTaus
# ########################################################################
from L1Trigger.L1TTrackMatch.L1TkBsCandidateProducer_cfi import *
l1TkBsCandidates = L1TkBsCandidates.clone()
l1TkBsCandidatesLooseWP = L1TkBsCandidatesLooseWP.clone()
l1TkBsCandidatesTightWP = L1TkBsCandidatesTightWP.clone()
phase2_SimL1Emulator += l1TkBsCandidates
phase2_SimL1Emulator += l1TkBsCandidatesLooseWP
phase2_SimL1Emulator += l1TkBsCandidatesTightWP

from Configuration.Eras.Modifier_phase2_trigger_cff import phase2_trigger
from Configuration.Eras.Modifier_phase2_trackerV14_cff import phase2_trackerV14
(phase2_trigger & phase2_trackerV14).toReplaceWith( SimL1EmulatorTask , _phase2_siml1emulator)
