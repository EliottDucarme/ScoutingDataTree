from CRABClient.UserUtilities import config
config = config()

config.General.requestName = ''
config.General.workArea = 'CRABDir'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../test/ProduceTree.py'
config.JobType.pyCfgParams = ['sampleType=ScoutingMuon2018']
# config.JobType.numCores = 4
# config.JobType.maxMemoryMB = 2500
# config.JobType.maxJobRuntimeMin = 2000

config.Data.inputDataset = ''
# config.Data.useParent = True

config.Data.inputDBS = 'global'
# config.Data.splitting = 'Automatic'
config.Data.splitting = 'LumiBased'
# config.Data.unitsPerJob = 20 # -- too many jobs ( > 1000)
config.Data.unitsPerJob = 50

config.Data.publication = False
# config.Data.ignoreLocality = True
# config.Site.whitelist = ['T2_KR_*', 'T2_US_*', 'T2_IT_*', 'T3_IT_*'] # -- mandatory for ignoreLocality option
config.Site.storageSite = 'T2_BE_IIHE'

config.Data.lumiMask = './JSON/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'

version = 'v02_2nd'
config.Data.outLFNDirBase = '/store/user/kplee/DYScoutingTree_%s' % version

config.JobType.allowUndistributedCMSSW = True

# 'MultiCRAB' part
if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand

    config.General.requestName = 'ScoutingCaloMuon_Run2018Av1_GoldenJSON'
    config.Data.inputDataset = '/ScoutingCaloMuon/Run2018A-v1/RAW'
    crabCommand('submit', config = config)

    config.General.requestName = 'ScoutingCaloMuon_Run2018Bv1_GoldenJSON'
    config.Data.inputDataset = '/ScoutingCaloMuon/Run2018B-v1/RAW'
    crabCommand('submit', config = config)

    config.General.requestName = 'ScoutingCaloMuon_Run2018Cv1_GoldenJSON'
    config.Data.inputDataset = '/ScoutingCaloMuon/Run2018C-v1/RAW'
    crabCommand('submit', config = config)

    config.General.requestName = 'ScoutingCaloMuon_Run2018Dv1_GoldenJSON'
    config.Data.inputDataset = '/ScoutingCaloMuon/Run2018D-v1/RAW'
    crabCommand('submit', config = config)
