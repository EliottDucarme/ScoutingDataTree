# NanoAOD-like ntuple producer for the scouting data

## Setup (for the first time)

```bash
# -- if it is not CentOS7 (e.g. lxplus)
# cmssw-el7 # -- use CentOS7 env. (for this CMSSW & ARCH)

cmsrel CMSSW_10_6_27
cd CMSSW_10_6_27/src
cmsenv

git clone git@github.com:KyeongPil-Lee/ScoutingDataTree.git
# https://github.com/KyeongPil-Lee/ScoutingDataTree.git # -- https

scram b

cd ScoutingDataTree/TreeProducer/test
voms-proxy-init --voms cms # -- to load exmaple RAW file
cmsRun ProduceTree.py sampleType=ScoutingMuon2018 >&ProduceTree.log&
```

### submit via CRAB

```bash
cd ../CRAB
python crabcfg_2018Data.py
```

