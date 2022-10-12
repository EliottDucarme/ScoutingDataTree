

def GetArgument(sampleType):
    exampleEDMFile = ''
    globalTag = ''
    isMC = ''
    isMiniAOD = ''

    if sampleType == "ScoutingMuon2018":
        # exampleEDMFile = '/store/data/Run2018A/ScoutingCaloMuon/RAW/v1/000/316/377/00000/6AAB14C9-B058-E811-B4CC-FA163EBBDB1E.root' # -- 9MB with 4114 events
        exampleEDMFile = '/store/data/Run2018D/ScoutingCaloMuon/RAW/v1/000/325/170/00000/29C6126D-C33D-2744-B4E8-92637A8EDE74.root'
        globalTag = '101X_dataRun2_HLT_v7'
        isMC = False
        isMiniAOD = False

    elif sampleType == "Data2018_MINIAOD":
        exampleEDMFile = '/store/data/Run2018C/SingleMuon/MINIAOD/17Sep2018-v1/270000/CBAACD5A-9EAF-9E46-BC1C-49C2EB361154.root'
        globalTag = '102X_dataRun2_Sep2018Rereco_v1'
        isMC = False
        isMiniAOD = True

    elif sampleType == "Data2018Prompt_MINIAOD":
        exampleEDMFile = '/store/data/Run2018C/SingleMuon/MINIAOD/17Sep2018-v1/270000/CBAACD5A-9EAF-9E46-BC1C-49C2EB361154.root'
        globalTag = '102X_dataRun2_Prompt_v11'
        isMC = False
        isMiniAOD = True

    elif sampleType == "MC2018": # -- Autumn18
         # exampleEDMFile = '/store/mc/RunIIAutumn18DRPremix/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v2/120000/362F2012-0EFD-6C45-8BFF-DED4B10C7F4F.root'
         exampleEDMFile = '/store/mc/RunIIAutumn18DRPremix/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/AODSIM/102X_upgrade2018_realistic_v15-v4/110000/8497240D-E938-D040-9688-2695C710FBFE.root' # -- QCD pT20-30, test
         globalTag = '102X_upgrade2018_realistic_v15'
         isMC = True
         isMiniAOD = False

    elif sampleType == "MC2018UL": # -- AOD
         exampleEDMFile = '/store/mc/RunIISummer19UL18RECO/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/100004/F83D534F-17F0-884A-91D4-E42F43132438.root'
         globalTag = '106X_upgrade2018_realistic_v11_L1v1'
         isMC = True
         isMiniAOD = False

    elif sampleType == "MC2018_MINIAOD":
         exampleEDMFile = '/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/A128C694-51C2-8744-BC1E-B5BD88CB8B8C.root'
         globalTag = '102X_upgrade2018_realistic_v15'
         isMC = True
         isMiniAOD = True

    elif sampleType == "MC2018UL_MINIAOD":
         exampleEDMFile = 'file:/afs/cern.ch/work/k/kplee/private/ROOTFile_Test/106X/UL18/MINIAODSIM_DYLL_M50_aMCNLO_Summer19UL18.root' # -- lxplus
         globalTag = '106X_upgrade2018_realistic_v11_L1v1'
         isMC = True
         isMiniAOD = True

    elif sampleType == "Data2018UL_MINIAOD":
         exampleEDMFile = 'file:/afs/cern.ch/work/k/kplee/private/ROOTFile_Test/106X/UL18/MINIAOD_EGamma_12Nov2019_UL2018.root' # -- lxplus
         # exampleEDMFile = 'file:/afs/cern.ch/work/k/kplee/private/ROOTFile_Test/106X/UL18/MINIAOD_DoubleMuon_12Nov2019_UL2018.root' # -- lxplus, for test
         globalTag = '106X_dataRun2_v28'
         isMC = False
         isMiniAOD = True

    else:
        print "[GetArgument] type = %s is not available sample type" % (sampleType)


    return exampleEDMFile, globalTag, isMC, isMiniAOD