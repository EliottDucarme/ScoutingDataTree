from CRABClient.UserUtilities import config
config = config()

config.General.requestName = ''
config.General.workArea = 'CRABDir'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../test/ProduceTree.py'
config.JobType.pyCfgParams = ['sampleType=ScoutingMuon2022']
# config.JobType.numCores = 4
# config.JobType.maxMemoryMB = 2500
# config.JobType.maxJobRuntimeMin = 600

config.Data.inputDataset = ''
# config.Data.useParent = True

config.Data.inputDBS = 'global'
# config.Data.splitting = 'Automatic'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 1 
# config.Data.unitsPerJob = 50

config.Data.publication = False
# config.Data.ignoreLocality = True
# config.Site.whitelist = ['T2_KR_*', 'T2_US_*', 'T2_IT_*', 'T3_IT_*'] # -- mandatory for ignoreLocality option
config.Site.storageSite = 'T2_BE_IIHE'

config.Data.lumiMask = './JSON/Cert_Collisions2022_355100_362760_Golden.json'

version = 'v01'
config.Data.outLFNDirBase = '/store/user/educarme/DYScoutingRun3Tree_%s' % version

config.JobType.allowUndistributedCMSSW = True

# config.General.requestName = 'ScoutingPFRun3_Run2022Fv1_GoldenJSON'
config.Data.inputDataset = '/ScoutingPFRun3/Run2022F-v1/RAW'
config.Data.runRange = '361438-361447'


# 'MultiCRAB' part
# if __name__ == '__main__':

# from CRABAPI.RawCommand import crabCommand

#   config.General.requestName = 'ScoutingPFRun3_Run2022Av1_GoldenJSON'
#   config.Data.inputDataset = '/ScoutingPFRun3/Run2022A-v1/RAW'
#   crabCommand('submit', config = config)

#   config.General.requestName = 'ScoutingPFRun3_Run2022Bv1_GoldenJSON'
#   config.Data.inputDataset = '/ScoutingPFRun3/Run2022B-v1/RAW'
#   crabCommand('submit', config = config)

#   config.General.requestName = 'ScoutingPFRun3_Run2022Cv1_GoldenJSON'
#   config.Data.inputDataset = '/ScoutingPFRun3/Run2022C-v1/RAW'
#   crabCommand('submit', config = config)

#   config.General.requestName = 'ScoutingPFRun3_Run2022Dv1_GoldenJSON'
#   config.Data.inputDataset = '/ScoutingPFRun3/Run2022C-v1/RAW'
#   crabCommand('submit', config = config)

#   config.General.requestName = 'ScoutingPFRun3_Run2022Ev1_GoldenJSON'
#   config.Data.inputDataset = '/ScoutingPFRun3/Run2022C-v1/RAW'
#   crabCommand('submit', config = config)

#   config.General.requestName = 'ScoutingPFRun3_Run2022Fv1_GoldenJSON'
#   config.Data.inputDataset = '/ScoutingPFRun3/Run2022F-v1/RAW'
#   crabCommand('submit', config = config)

#   config.General.requestName = 'ScoutingPFRun3_Run2022Gv1_GoldenJSON'
#   config.Data.inputDataset = '/ScoutingPFRun3/Run2022C-v1/RAW'
#   crabCommand('submit', config = config)
