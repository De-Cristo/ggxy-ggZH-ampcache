# Pythia8 hadronizer for the LO companion sample (Born-only events: bornonly 1, LOevents 1):
# no POWHEG emission, so no POWHEG veto; the shower starts at SCALUP = sqrt(s-hat).
# CP5 tune (Run 3), H -> bb.  Configuration of the LO validation sample.
import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunesRun3ECM13p6TeV.PythiaCP5Settings_cfi import *

generator = cms.EDFilter("Pythia8HadronizerFilter",
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13600.),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        processParameters = cms.vstring(
            '25:m0 = 125.0',
            '25:onMode = off',
            '25:onIfMatch = 5 -5',
        ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CP5Settings',
                                    'processParameters'),
    ),
)

ProductionFilterSequence = cms.Sequence(generator)
