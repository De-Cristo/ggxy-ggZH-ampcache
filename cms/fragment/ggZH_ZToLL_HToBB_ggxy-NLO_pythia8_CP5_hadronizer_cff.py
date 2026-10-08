# Pythia8 hadronizer for the POWHEG-BOX-V2 ggxy_ggZH NLO gridpack (gg -> Z(l+l-) H, H -> bb):
# CP5 tune (Run 3), POWHEG emission veto with POWHEG:nFinal = 3 (l+, l-, H at the Born level;
# the value used for the validation samples, without Pythia errors on the final-state count).
# This is the configuration of the validation samples (LHE input, step GEN,NANOGEN);
# ggZH_ZToLL_HToBB_ggxy-NLO_pythia8_CP5_cff.py is the same with the ExternalLHEProducer.
import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunesRun3ECM13p6TeV.PythiaCP5Settings_cfi import *
from Configuration.Generator.Pythia8PowhegEmissionVetoSettings_cfi import *
# from Configuration.Generator.PSweightsPythia.PythiaPSweightsSettings_cfi import *   # optional PS weights

generator = cms.EDFilter("Pythia8HadronizerFilter",
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13600.),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        pythia8PowhegEmissionVetoSettingsBlock,
        # pythia8PSweightsSettingsBlock,
        processParameters = cms.vstring(
            'POWHEG:nFinal = 3',      # l+ l- H
            '25:m0 = 125.0',
            '25:onMode = off',
            '25:onIfMatch = 5 -5',    # H -> bb; change for other Higgs decays
        ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CP5Settings',
                                    'pythia8PowhegEmissionVetoSettings',
                                    # 'pythia8PSweightsSettings',
                                    'processParameters'),
    ),
)

ProductionFilterSequence = cms.Sequence(generator)
