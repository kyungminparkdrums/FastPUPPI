import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar 

def LazyVar(expr, valtype, doc=None, precision=-1):
    return Var(expr, valtype, doc, precision, lazyEval=True)

process = cms.Process("L1Dump", eras.Phase2C17I13M9)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/cms/store/cmst3/group/l1tr/vcamagni/Zee/GEN-SIM-DIGI-RAW/m20/2STEPS/140X_Phase2Spring24_FEVTDEBUGHLT_v0/prod_140X_6516663_37.root')
    #fileNames = cms.untracked.vstring('file:/eos/cms/store/cmst3/group/l1tr/kypark/REPO_Scouting_Hackathon_December/NuGunAllEta_PU200/FP/v140X/perfNano_7807692_6607.root')
)

#process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2026D88_cff')
process.load('Configuration.Geometry.GeometryExtended2026D95Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D95_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi') # needed for HGCAL_noise_fC
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '125X_mcRun4_realistic_v2', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_mcRun4_realistic_v3', '')

process.l1tTrackSelectionProducer.processSimulatedTracks = False # these would need stubs, and are not used anyway

process.deps = cms.Task(
    #process.l1tTkMuonsGmt,
    #process.l1tSAMuonsGmt,
    #process.l1tGTTInputProducer,
    #process.l1tTrackSelectionProducer,
    #process.l1tVertexFinderEmulator,
    #process.L1TLayer1TaskInputsTask,
    #process.L1TLayer1Task,
    process.l1tLayer2EG,
    #process.L1TPFJetsEmulationTask,
    #process.L1TPFJetsExtendedTask,
    #process.L1TBJetsTask
)

process.egDump = cms.EDAnalyzer("L1CTL2EgammaBinaryDumper",
                srcEle = cms.InputTag("l1tLayer2EG", "L1CtTkElectron"),
                srcEm = cms.InputTag("l1tLayer2EG", "L1CtTkEm"),
                #srcPuppi = cms.InputTag("l1tLayer2Deregionizer", "Puppi"),
                interleaveOutputs = cms.bool(False), # False = first 12 photons, then electrons; True = pho1, ele1, pho2, ele2, ...
                outName = cms.string("egamma.dump"))

process.p_dumps = cms.EndPath(
    process.egDump
)
process.p_dumps.associate(process.deps)

process.phoTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("l1tLayer2EG","L1CtTkEm"),
    cut = cms.string(""),
    name = cms.string("Pho"),
    doc = cms.string("Photons (TkEm) from CTL2"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table
    variables = cms.PSet(
        pt   = Var("pt",  float),
        eta  = Var("eta", float),
        phi  = Var("phi", float),
        mass = Var("0", float),
        quality = LazyVar("hwQual", int, doc="quality (TBD)"),
        trkIsol = LazyVar("trkIsol", float),
        trkIsolPV = LazyVar("trkIsolPV", float),
        puppiIsol = LazyVar("puppiIsol", float),
        puppiIsolPV = LazyVar("puppiIsolPV", float),
        hwPt   = LazyVar("hwPt", int),
        hwEta  = LazyVar("hwEta", int),
        hwPhi  = LazyVar("hwPhi", int),
        hwQual  = LazyVar("hwQual", int),
    )
)

process.eleTable = process.phoTable.clone(
    src = "l1tLayer2EG:L1CtTkElectron",
    cut = "",
    name = "Ele",
    doc = "TkElectrons from CTL2",
    variables = dict(
        z0   = LazyVar("trkzVtx",  float),
        idScore = LazyVar("idScore",  float),
    )
)
process.genPho = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("genParticles"),
    cut = cms.string("pdgId == 22 && status == 1 && (pt > 10 || isPromptFinalState())"),
)
process.genEle = process.genPho.clone(
    cut = "pdgId == 11 && status == 1",
)
process.genPhoTable = cms.EDProducer("SimpleGenParticleFlatTableProducer",
    src = cms.InputTag("genPho"),
    cut = cms.string(""),
    name = cms.string("GenPho"),
    doc = cms.string("gen photons (pt > 10 || prompt)"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table
    variables = cms.PSet(
        pt  = Var("pt",  float),
        phi = Var("phi", float),
        eta  = Var("eta", float),
        mass  = Var("mass", float),
        z0   = Var("vz",  float, doc="Production point along the beam axis"),
        dxy   = Var("vertex.Rho",  float, doc="transverse distance of production point from the beam axis"),
        isPrompt  = Var("statusFlags().isPrompt()", int, doc="Prompt"),
        isFromTau  = Var("statusFlags().isDirectPromptTauDecayProduct()", int, doc="Electron from prompt tau decay"),
        motherId  = Var("? motherRef.isNonnull() ? motherRef.pdgId : 0", int, doc="Mother particle ID"),
    )
)
process.genEleTable = process.genPhoTable.clone(
    src = "genEle",
    cut = "",
    name = "GenEle",
    doc = "gen photons (pt > 10 || prompt)",
)
process.phoMCMatch = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = process.phoTable.src,                         # final reco collection
    matched     = cms.InputTag("genPho"),     # final mc-truth particle collection
    mcPdgId     = cms.vint32(22),               # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),              # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.2),              # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),    # False = just match input in order; True = pick lowest deltaR pair first
)
process.eleMCMatch = process.phoMCMatch.clone(
    src         = process.eleTable.src,  # final reco collection
    matched     = "genEle",              # final mc-truth particle collection
    mcPdgId     = [11],                  # one or more PDG ID (13 = mu); absolute values (see below)
)

process.phoMCTable = cms.EDProducer("CandMCMatchTableProducer",
    src     = process.phoTable.src,
    mcMap   = cms.InputTag("phoMCMatch"),
    objName = process.phoTable.name,
    objType = cms.string("Other"),
    branchName = cms.string("GenPho"),
    docString = cms.string("MC matching"),
)
process.eleMCTable = process.phoMCTable.clone(
    src     = process.eleTable.src,
    mcMap   = "eleMCMatch",
    objName = process.eleTable.name,
    branchName = "GenEle",
)
process.p_pho = cms.Path(process.phoTable)
process.p_phoMC = cms.Path(process.genPho + process.genPhoTable + process.phoMCMatch + process.phoMCTable)
process.p_ele = cms.Path(process.eleTable)
process.p_eleMC = cms.Path(process.genEle + process.genEleTable + process.eleMCMatch + process.eleMCTable)

process.outnano = cms.OutputModule("NanoAODOutputModule",
    fileName = cms.untracked.string("l1Nano.root"),
    outputCommands = cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*Table_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("ZLIB"),
)
process.end = cms.EndPath(process.outnano)

def noPU():
    for X in "", "EM", "Raw", "EMRaw":
        pfc = getattr(process, 'pfClustersFromHGC3DClusters'+X, None)
        if not pfc: continue
        pfc.emVsPUID.wp = "-1.0"

def noNano():
    process.schedule = cms.Schedule(process.p_dumps)

