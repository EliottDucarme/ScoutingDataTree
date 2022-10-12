import FWCore.ParameterSet.Config as cms

"The HLT objects saved in the tree are the ones matched with the filter names from 2018 HLT menu"
"DYTreeProducer::SavedFilterCondition should be updated if you want to run it on 2016 or 17!"

# -- usage: usage: cmsRun ProduceTree.py sampleType=<sample type>

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.register('sampleType',
                  "none", # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.string,         # string, int, or float
                  "Sample type defined in DYScouting/TreeProducer/NtuplerArgument.py")

options.parseArguments()

print "input sample type = ", options.sampleType

from DYScouting.TreeProducer.NtuplerArgument import GetArgument
theExampleFile, theGlobalTag, isMC, isMiniAOD = GetArgument( options.sampleType )

print "   [example file] ", theExampleFile
print "   [global tag]   ", theGlobalTag
print "   [isMC]         ", isMC
print "   [isMiniAOD]    ", isMiniAOD
print ""

process = cms.Process("TreeProducer")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(theExampleFile), # -- @ KNU
    secondaryFileNames = cms.untracked.vstring(),
    # lumisToProcess = cms.untracked.VLuminosityBlockRange('258158:1-258158:1786'),
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = theGlobalTag

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')

process.TFileService = cms.Service("TFileService",
  fileName = cms.string("ntuple.root"),
  closeFileFast = cms.untracked.bool(False),
)

# -- produce PAT trigger object (to use the trigger information with a pre-defined format)
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi")
process.patTrigger.onlyStandAlone = True # -- produce triggerObjectStandAlone objects only

# -- extract IterL3MuonCandidateNoVtx object (TriggerObjectStandAlone) from patTrigger
process.IterL3MuonCandidatesNoVtx = cms.EDProducer("TriggerObjectFilterByCollection",
    src = cms.InputTag("patTrigger"),
    collections = cms.vstring("hltIterL3MuonCandidatesNoVtx"),
)

# -- for the extrapolation of offlie muon to 2nd muon station
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")

##--- l1 stage2 digis ---
process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )
process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')

from DYScouting.TreeProducer.L1SeedList import GetL1SeedList

process.DYTree = cms.EDAnalyzer('DYTreeProducer',
  L1Muon                = cms.untracked.InputTag("gmtStage2Digis", "Muon", "RECO"),
  globalAlgBlk          = cms.untracked.InputTag("gtStage2Digis"),
  l1tAlgBlkInputTag     = cms.InputTag("gtStage2Digis"), # -- for L1TGlobalUtil
  l1tExtBlkInputTag     = cms.InputTag("gtStage2Digis"), # -- for L1TGlobalUtil
  ReadPrescalesFromFile = cms.bool( False ),             # -- for L1TGlobalUtil
  triggerResults        = cms.untracked.InputTag("TriggerResults", "", "HLT"),
  triggerEvent          = cms.untracked.InputTag("hltTriggerSummaryAOD"), # -- for the trigger objects in AOD: will be skipped if the collection is not available (e.g. RAW)

  SCDimuonVtx      = cms.untracked.InputTag("hltScoutingMuonPackerCalo",              "displacedVtx", "HLT"),
  SCPixelVtx       = cms.untracked.InputTag("hltScoutingPrimaryVertexPacker",         "primaryVtx",   "HLT"),
  SCPixelVtxNearMu = cms.untracked.InputTag("hltScoutingPrimaryVertexPackerCaloMuon", "primaryVtx",   "HLT"),

  SCMuon           = cms.untracked.InputTag("hltScoutingMuonPackerCalo"),
  SCCaloJet        = cms.untracked.InputTag("hltScoutingCaloPacker"),
  SCCaloMETPhi     = cms.untracked.InputTag("hltScoutingCaloPacker", "caloMetPhi", "HLT"),
  SCCaloMETPt      = cms.untracked.InputTag("hltScoutingCaloPacker", "caloMetPt", "HLT"),
  SCRho            = cms.untracked.InputTag("hltScoutingCaloPacker", "rho", "HLT"),

  # triggerObject_L3MuonNoVtx = cms.untracked.InputTag("IterL3MuonCandidatesNoVtx"),
)

# RAW data tier: RAWtoDigi step is needed to retrieve L1 information
process.p = cms.Path(process.patTrigger + process.IterL3MuonCandidatesNoVtx + process.gtStage2Digis + process.DYTree)
