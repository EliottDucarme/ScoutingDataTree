// -- Flat tree producer with the nanoAOD-like format for the scouting data
// -- based on https://github.com/KyeongPil-Lee/DYScouting/tree/master/TreeProducer
// -- author: Kyeongpil Lee (ULB)

// -- system include files
#include <memory>


// -- frameworks
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// -- triggers
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

// -- scouting
#include "DataFormats/Scouting/interface/ScoutingMuon.h"
#include "DataFormats/Scouting/interface/ScoutingCaloJet.h"
#include "DataFormats/Scouting/interface/ScoutingVertex.h"

// -- tracking
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

// -- others
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

// -- ROOT includes
#include <TTree.h>
#include <TLorentzVector.h>

using namespace std;
using namespace reco;
using namespace edm;

class ScoutingDataTreeProducer : public edm::EDAnalyzer
{
public:
  explicit ScoutingDataTreeProducer(const edm::ParameterSet&);
  ~ScoutingDataTreeProducer();

private:
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun( const edm::Run&, const edm::EventSetup& );
  virtual void endRun( const edm::Run&, const edm::EventSetup& );

  void Init();
  void Make_Branch();

  void Fill_L1( const edm::Event&, const edm::EventSetup& );
  void Fill_HLT( const edm::Event & );
  // void Fill_L3MuonNoVtx( const edm::Event& );
  void Fill_SCDimuonVtx( const edm::Event& );
  void Fill_SCPixelVtx( const edm::Event& );
  void Fill_SCPixelVtxNearMu( const edm::Event& );
  void Fill_SCMuon( const edm::Event& );
  void Fill_SCCaloJet( const edm::Event& );

  void SetTrue_HLTBitInfo( const std::string& );
  void Set_L1BitAndPrescaleInfo();

  void GetMuonIndex_AssociatedToVertex(edm::Event& iEvent, const ScoutingVertex& vtx, int theVtxIndex, int& index1_mu, int& index2_mu);


  // -- tokens
  edm::EDGetTokenT< l1t::MuonBxCollection >       t_L1Muon_;

  edm::EDGetTokenT< BXVector<GlobalAlgBlk> >      t_globalAlgBlk_;

  edm::EDGetTokenT< edm::TriggerResults >         t_triggerResults_;
  // edm::EDGetTokenT< trigger::TriggerEvent >       t_triggerEvent_;

  edm::EDGetTokenT<std::vector<ScoutingVertex> >  t_SCDimuonVtx_;
  edm::EDGetTokenT<std::vector<ScoutingVertex> >  t_SCPixelVtx_;
  edm::EDGetTokenT<std::vector<ScoutingVertex> >  t_SCPixelVtxNearMu_;

  edm::EDGetTokenT<std::vector<ScoutingMuon> >    t_SCMuon_;
  edm::EDGetTokenT<std::vector<ScoutingCaloJet> > t_SCCaloJet_;

  edm::EDGetTokenT< double >  t_SCCaloMETPhi_;
  edm::EDGetTokenT< double >  t_SCCaloMETPt_;
  edm::EDGetTokenT< double >  t_SCRho_;


  // -- variable for L1 information
  l1t::L1TGlobalUtil   *L1GtUtils_;

  // -- IterL3MuonCandidateNoVtx object
  // edm::EDGetTokenT< std::vector<pat::TriggerObjectStandAlone> > t_trigObj_L3MuonNoVtx_;

  // -- debugging switch
  bool debug_ = false;

  // -- tree
  TTree* ntuple_;

  // -- array size
  static const int arrSize_ = 2000;

  // --  event informations
  int run_;
  int luminosityBlock_;
  unsigned long long event_;

  // -- L1 flags
  bool   L1_DoubleMu_15_7_;
  bool   L1_DoubleMu4p5_SQ_OS_dR_Max1p2_;

  // -- HLT flags
  bool   DST_DoubleMu1_noVtx_CaloScouting_;
  bool   DST_DoubleMu3_noVtx_CaloScouting_Monitoring_;
  bool   DST_DoubleMu3_noVtx_CaloScouting_;
  bool   HLT_IsoMu24_;

  // -- HLT objects
  // -- use float (instead of double) to sync. with nanoAOD MC
  int   nTrigObj_;
  float TrigObj_pt[arrSize_];
  float TrigObj_eta[arrSize_];
  float TrigObj_phi[arrSize_];

  // -- dimuon vertex information (@ HLT)
  int   nSCDimuonVtx_;
  float SCDimuonVtx_x_[arrSize_];
  float SCDimuonVtx_y_[arrSize_];
  float SCDimuonVtx_z_[arrSize_];
  float SCDimuonVtx_xErr_[arrSize_];
  float SCDimuonVtx_yErr_[arrSize_];
  float SCDimuonVtx_zErr_[arrSize_];
  float SCDimuonVtx_chi2_[arrSize_];
  int   SCDimuonVtx_nDOF_[arrSize_];
  int   SCDimuonVtx_muonIndex1_[arrSize_];
  int   SCDimuonVtx_muonIndex2_[arrSize_];
  int   SCDimuonVtx_isValid_[arrSize_];

  // -- pixel vertex information from full tracking @ HLT
  int   nSCPixelVtx_;
  float SCPixelVtx_x_[arrSize_];
  float SCPixelVtx_y_[arrSize_];
  float SCPixelVtx_z_[arrSize_];
  float SCPixelVtx_xErr_[arrSize_];
  float SCPixelVtx_yErr_[arrSize_];
  float SCPixelVtx_zErr_[arrSize_];
  float SCPixelVtx_chi2_[arrSize_];
  int   SCPixelVtx_nDOF_[arrSize_];
  int   SCPixelVtx_isValid_[arrSize_];

  // -- pixel vertex made by tracks near L3 muons @ HLT
  int nSCPixelVtxNearMu_;
  float SCPixelVtxNearMu_x_[arrSize_];
  float SCPixelVtxNearMu_y_[arrSize_];
  float SCPixelVtxNearMu_z_[arrSize_];
  float SCPixelVtxNearMu_xErr_[arrSize_];
  float SCPixelVtxNearMu_yErr_[arrSize_];
  float SCPixelVtxNearMu_zErr_[arrSize_];
  float SCPixelVtxNearMu_chi2_[arrSize_];
  int   SCPixelVtxNearMu_nDOF_[arrSize_];
  int   SCPixelVtxNearMu_isValid_[arrSize_];

  // -- muon information (scouting)
  int nSCMuon_;
  float SCMuon_pt_[arrSize_];
  float SCMuon_eta_[arrSize_];
  float SCMuon_phi_[arrSize_];
  float SCMuon_charge_[arrSize_];

  int SCMuon_nPixelHit_[arrSize_];
  int SCMuon_nStripHit_[arrSize_];
  int SCMuon_nTrackerLayer_[arrSize_];
  int SCMuon_nMuonHit_[arrSize_];
  int SCMuon_nMatchedStation_[arrSize_];
  float SCMuon_normChi2_[arrSize_];
  float SCMuon_dxy_[arrSize_];
  float SCMuon_dz_[arrSize_];
  float SCMuon_trkIso_[arrSize_];

  // -- jet information
  int nSCCaloJet_;
  double SCCaloJet_pt_[arrSize_];
  double SCCaloJet_eta_[arrSize_];
  double SCCaloJet_phi_[arrSize_];
  double SCCaloJet_m_[arrSize_];
  double SCCaloJet_jetArea_[arrSize_];
  double SCCaloJet_maxEInEmTowers_[arrSize_];
  double SCCaloJet_maxEInHadTowers_[arrSize_];
  double SCCaloJet_hadEnergyInHB_[arrSize_];
  double SCCaloJet_hadEnergyInHE_[arrSize_];
  double SCCaloJet_hadEnergyInHF_[arrSize_];
  double SCCaloJet_emEnergyInEB_[arrSize_];
  double SCCaloJet_emEnergyInEE_[arrSize_];
  double SCCaloJet_emEnergyInHF_[arrSize_];
  double SCCaloJet_towersArea_[arrSize_];
  double SCCaloJet_mvaDiscriminator_[arrSize_];
  double SCCaloJet_btagDiscriminator_[arrSize_];

  // -- MET, Calo info.
  double SCCaloMET_phi_;
  double SCCaloMET_pt_;
  double SCRho_;

  // int nL1Muon_;
  // double L1Muon_pt_[arrSize_];
  // double L1Muon_eta_[arrSize_];
  // double L1Muon_phi_[arrSize_];
  // double L1Muon_charge_[arrSize_];
  // double L1Muon_quality_[arrSize_];

  // int nL3MuonNoVtx_;
  // double L3MuonNoVtx_pt_[arrSize_];
  // double L3MuonNoVtx_eta_[arrSize_];
  // double L3MuonNoVtx_phi_[arrSize_];
};

ScoutingDataTreeProducer::ScoutingDataTreeProducer(const edm::ParameterSet& iConfig):
t_L1Muon_              ( consumes< l1t::MuonBxCollection  >         (iConfig.getUntrackedParameter<edm::InputTag>("L1Muon")) ),
t_globalAlgBlk_        ( consumes< BXVector< GlobalAlgBlk > >       (iConfig.getUntrackedParameter<edm::InputTag>("globalAlgBlk")) ),
t_triggerResults_      ( consumes< edm::TriggerResults >            (iConfig.getUntrackedParameter<edm::InputTag>("triggerResults")) ),
t_SCDimuonVtx_         ( consumes< std::vector<ScoutingVertex> >    (iConfig.getUntrackedParameter<edm::InputTag>("SCDimuonVtx")) ),
t_SCPixelVtx_          ( consumes< std::vector<ScoutingVertex> >    (iConfig.getUntrackedParameter<edm::InputTag>("SCPixelVtx")) ),
t_SCPixelVtxNearMu_    ( consumes< std::vector<ScoutingVertex> >    (iConfig.getUntrackedParameter<edm::InputTag>("SCPixelVtxNearMu")) ),
t_SCMuon_              ( consumes< std::vector<ScoutingMuon> >      (iConfig.getUntrackedParameter<edm::InputTag>("SCMuon")) ),
t_SCCaloJet_           ( consumes< std::vector<ScoutingCaloJet> >   (iConfig.getUntrackedParameter<edm::InputTag>("SCCaloJet")) ),
t_SCCaloMETPhi_        ( consumes< double >                         (iConfig.getUntrackedParameter<edm::InputTag>("SCCaloMETPhi")) ),
t_SCCaloMETPt_         ( consumes< double >                         (iConfig.getUntrackedParameter<edm::InputTag>("SCCaloMETPt")) ),
t_SCRho_               ( consumes< double >                         (iConfig.getUntrackedParameter<edm::InputTag>("SCRho")) )
// t_trigObj_L3MuonNoVtx_ ( consumes< std::vector<pat::TriggerObjectStandAlone> >   (iConfig.getUntrackedParameter<edm::InputTag>("triggerObject_L3MuonNoVtx")) ),
{
  L1GtUtils_  = new l1t::L1TGlobalUtil(iConfig, consumesCollector());
}

ScoutingDataTreeProducer::~ScoutingDataTreeProducer()
{

}

void ScoutingDataTreeProducer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  Init();

  run_             = iEvent.id().run();
  luminosityBlock_ = iEvent.id().luminosityBlock();
  event_           = iEvent.id().event();

  // -- fill each object
  Fill_L1(iEvent, iSetup);
  Fill_HLT(iEvent);
  // Fill_L3MuonNoVtx(iEvent);
  Fill_SCDimuonVtx(iEvent);
  Fill_SCPixelVtx(iEvent);
  Fill_SCPixelVtxNearMu(iEvent);
  Fill_SCMuon(iEvent);
  Fill_SCCaloJet(iEvent);

  // -- MET variables
  edm::Handle<double> h_SCCaloMETPhi;
  iEvent.getByToken(t_SCCaloMETPhi_, h_SCCaloMETPhi);
  if( h_SCCaloMETPhi.isValid() ) SCCaloMET_phi_ = *(h_SCCaloMETPhi.product());

  edm::Handle<double> h_SCCaloMETPt;
  iEvent.getByToken(t_SCCaloMETPt_, h_SCCaloMETPt);
  if( h_SCCaloMETPt.isValid() ) SCCaloMET_pt_ = *(h_SCCaloMETPt.product());

  edm::Handle<double> h_SCRho;
  iEvent.getByToken(t_SCRho_, h_SCRho);
  if( h_SCRho.isValid() ) SCRho_ = *(h_SCRho.product());

  ntuple_->Fill();
}

void ScoutingDataTreeProducer::beginJob()
{
  edm::Service<TFileService> fs;
  ntuple_ = fs->make<TTree>("Events","Events");

  Make_Branch();
}

void ScoutingDataTreeProducer::Init()
{
  // --  event information
  run_             = -999;
  luminosityBlock_ = -999;
  event_           = 0;

  // -- L1 flags
  L1_DoubleMu_15_7_ = false;
  L1_DoubleMu4p5_SQ_OS_dR_Max1p2_ = false;

  // -- HLT flags
  DST_DoubleMu1_noVtx_CaloScouting_ = false;
  DST_DoubleMu3_noVtx_CaloScouting_Monitoring_ = false;
  DST_DoubleMu3_noVtx_CaloScouting_ = false;
  HLT_IsoMu24_ = false;

  // -- HLT objects
  nTrigObj_ = -999;
  for(Int_t i=0; i<arrSize_; i++) {
    TrigObj_pt[i] = -999;
    TrigObj_eta[i] = -999;
    TrigObj_phi[i] = -999;
  }


  // -- vertex information (@ HLT)
  nSCDimuonVtx_      = -999;
  nSCPixelVtx_       = -999;
  nSCPixelVtxNearMu_ = -999;
  for(Int_t i=0; i<arrSize_; i++) {
    SCDimuonVtx_x_[i] = -999;
    SCDimuonVtx_y_[i] = -999;
    SCDimuonVtx_z_[i] = -999;
    SCDimuonVtx_xErr_[i] = -999;
    SCDimuonVtx_yErr_[i] = -999;
    SCDimuonVtx_zErr_[i] = -999;
    SCDimuonVtx_chi2_[i] = -999;
    SCDimuonVtx_nDOF_[i] = -999;
    SCDimuonVtx_muonIndex1_[i] = -999;
    SCDimuonVtx_muonIndex2_[i] = -999;
    SCDimuonVtx_isValid_[i] = 0;

    SCPixelVtx_x_[i] = -999;
    SCPixelVtx_y_[i] = -999;
    SCPixelVtx_z_[i] = -999;
    SCPixelVtx_xErr_[i] = -999;
    SCPixelVtx_yErr_[i] = -999;
    SCPixelVtx_zErr_[i] = -999;
    SCPixelVtx_chi2_[i] = -999;
    SCPixelVtx_nDOF_[i] = -999;
    SCPixelVtx_isValid_[i] = 0;

    SCPixelVtxNearMu_x_[i] = -999;
    SCPixelVtxNearMu_y_[i] = -999;
    SCPixelVtxNearMu_z_[i] = -999;
    SCPixelVtxNearMu_xErr_[i] = -999;
    SCPixelVtxNearMu_yErr_[i] = -999;
    SCPixelVtxNearMu_zErr_[i] = -999;
    SCPixelVtxNearMu_chi2_[i] = -999;
    SCPixelVtxNearMu_nDOF_[i] = -999;
    SCPixelVtxNearMu_isValid_[i] = 0;
  }

  // -- muon information
  nSCMuon_ = -999;
  for(Int_t i=0; i<arrSize_; i++)
  {
    SCMuon_pt_[i] = -999;
    SCMuon_eta_[i] = -999;
    SCMuon_phi_[i] = -999;
    SCMuon_charge_[i] = -999;

    SCMuon_nPixelHit_[i] = -999;
    SCMuon_nStripHit_[i] = -999;
    SCMuon_nTrackerLayer_[i] = -999;
    SCMuon_nMuonHit_[i] = -999;
    SCMuon_nMatchedStation_[i] = -999;
    SCMuon_normChi2_[i] = -999;
    SCMuon_dxy_[i] = -999;
    SCMuon_dz_[i] = -999;
    SCMuon_trkIso_[i] = -999;
  }

  // -- jet information
  nSCCaloJet_ = -999;
  for(Int_t i=0; i<arrSize_; i++)
  {
    SCCaloJet_pt_[i] = -999;
    SCCaloJet_eta_[i] = -999;
    SCCaloJet_phi_[i] = -999;
    SCCaloJet_m_[i] = -999;
    SCCaloJet_jetArea_[i] = -999;
    SCCaloJet_maxEInEmTowers_[i] = -999;
    SCCaloJet_maxEInHadTowers_[i] = -999;
    SCCaloJet_hadEnergyInHB_[i] = -999;
    SCCaloJet_hadEnergyInHE_[i] = -999;
    SCCaloJet_hadEnergyInHF_[i] = -999;
    SCCaloJet_emEnergyInEB_[i] = -999;
    SCCaloJet_emEnergyInEE_[i] = -999;
    SCCaloJet_emEnergyInHF_[i] = -999;
    SCCaloJet_towersArea_[i] = -999;
    SCCaloJet_mvaDiscriminator_[i] = -999;
    SCCaloJet_btagDiscriminator_[i] = -999;
  }

  SCCaloMET_phi_ = -999;
  SCCaloMET_pt_ = -999;
  SCRho_ = -999;

  // nL1Muon_ = -999;
  // for( int i=0; i<arrSize_; i++)
  // {
  //   L1Muon_pt_[i] = -999;
  //   L1Muon_eta_[i] = -999;
  //   L1Muon_phi_[i] = -999;
  //   L1Muon_charge_[i] = -999;
  //   L1Muon_quality_[i] = -999;
  // }

  // nL3MuonNoVtx_ = -999;
  // for( int i=0; i<arrSize_; i++)
  // {
  //   L3MuonNoVtx_pt_[i] = -999;
  //   L3MuonNoVtx_eta_[i] = -999;
  //   L3MuonNoVtx_phi_[i] = -999;
  // }

}

void ScoutingDataTreeProducer::Make_Branch()
{
  ntuple_->Branch("run",&run_,"run/I");
  ntuple_->Branch("luminosityBlock",&luminosityBlock_,"luminosityBlock/I");
  ntuple_->Branch("event",&event_,"event/l"); // -- unsigned long long -- //

  ntuple_->Branch("L1_DoubleMu_15_7",               &L1_DoubleMu_15_7_,               "L1_DoubleMu_15_7/O");
  ntuple_->Branch("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2_, "L1_DoubleMu4p5_SQ_OS_dR_Max1p2/O");

  ntuple_->Branch("DST_DoubleMu1_noVtx_CaloScouting",            &DST_DoubleMu1_noVtx_CaloScouting_,            "DST_DoubleMu1_noVtx_CaloScouting/O");
  ntuple_->Branch("DST_DoubleMu3_noVtx_CaloScouting_Monitoring", &DST_DoubleMu3_noVtx_CaloScouting_Monitoring_, "DST_DoubleMu3_noVtx_CaloScouting_Monitoring/O");
  ntuple_->Branch("DST_DoubleMu3_noVtx_CaloScouting",            &DST_DoubleMu3_noVtx_CaloScouting_,            "DST_DoubleMu3_noVtx_CaloScouting/O");
  ntuple_->Branch("HLT_IsoMu24",                                 &HLT_IsoMu24_,                                 "HLT_IsoMu24/O");

  // ntuple_->Branch("nTrigObj", &nTrigObj_, "nTrigObj/I");
  // ntuple_->Branch("TrigObj_pt",  &TrigObj_pt_, "TrigObj_pt[nTrigObj]/F");
  // ntuple_->Branch("TrigObj_eta", &TrigObj_eta_, "TrigObj_eta[nTrigObj]/F");
  // ntuple_->Branch("TrigObj_phi", &TrigObj_phi_, "TrigObj_phi[nTrigObj]/F");

  ntuple_->Branch("nSCDimuonVtx", &nSCDimuonVtx_, "nSCDimuonVtx/I");
  ntuple_->Branch("SCDimuonVtx_x", &SCDimuonVtx_x_, "SCDimuonVtx_x[nSCDimuonVtx]/F");
  ntuple_->Branch("SCDimuonVtx_y", &SCDimuonVtx_y_, "SCDimuonVtx_y[nSCDimuonVtx]/F");
  ntuple_->Branch("SCDimuonVtx_z", &SCDimuonVtx_z_, "SCDimuonVtx_z[nSCDimuonVtx]/F");
  ntuple_->Branch("SCDimuonVtx_xErr", &SCDimuonVtx_xErr_, "SCDimuonVtx_xErr[nSCDimuonVtx]/F");
  ntuple_->Branch("SCDimuonVtx_yErr", &SCDimuonVtx_yErr_, "SCDimuonVtx_yErr[nSCDimuonVtx]/F");
  ntuple_->Branch("SCDimuonVtx_zErr", &SCDimuonVtx_zErr_, "SCDimuonVtx_zErr[nSCDimuonVtx]/F");
  ntuple_->Branch("SCDimuonVtx_chi2", &SCDimuonVtx_chi2_, "SCDimuonVtx_chi2[nSCDimuonVtx]/F");
  ntuple_->Branch("SCDimuonVtx_nDOF", &SCDimuonVtx_nDOF_, "SCDimuonVtx_nDOF[nSCDimuonVtx]/I");
  ntuple_->Branch("SCDimuonVtx_muonIndex1", &SCDimuonVtx_muonIndex1_, "SCDimuonVtx_muonIndex1[nSCDimuonVtx]/I");
  ntuple_->Branch("SCDimuonVtx_muonIndex2", &SCDimuonVtx_muonIndex2_, "SCDimuonVtx_muonIndex2[nSCDimuonVtx]/I");
  ntuple_->Branch("SCDimuonVtx_isValid", &SCDimuonVtx_isValid_, "SCDimuonVtx_isValid[nSCDimuonVtx]/I");


  ntuple_->Branch("nSCPixelVtx", &nSCPixelVtx_, "nSCPixelVtx/I");
  ntuple_->Branch("SCPixelVtx_x", &SCPixelVtx_x_, "SCPixelVtx_x[nSCPixelVtx]/F");
  ntuple_->Branch("SCPixelVtx_y", &SCPixelVtx_y_, "SCPixelVtx_y[nSCPixelVtx]/F");
  ntuple_->Branch("SCPixelVtx_z", &SCPixelVtx_z_, "SCPixelVtx_z[nSCPixelVtx]/F");
  ntuple_->Branch("SCPixelVtx_xErr", &SCPixelVtx_xErr_, "SCPixelVtx_xErr[nSCPixelVtx]/F");
  ntuple_->Branch("SCPixelVtx_yErr", &SCPixelVtx_yErr_, "SCPixelVtx_yErr[nSCPixelVtx]/F");
  ntuple_->Branch("SCPixelVtx_zErr", &SCPixelVtx_zErr_, "SCPixelVtx_zErr[nSCPixelVtx]/F");
  ntuple_->Branch("SCPixelVtx_chi2", &SCPixelVtx_chi2_, "SCPixelVtx_chi2[nSCPixelVtx]/F");
  ntuple_->Branch("SCPixelVtx_nDOF", &SCPixelVtx_nDOF_, "SCPixelVtx_nDOF[nSCPixelVtx]/I");
  ntuple_->Branch("SCPixelVtx_isValid", &SCPixelVtx_isValid_, "SCPixelVtx_isValid[nSCPixelVtx]/I");


  ntuple_->Branch("nSCPixelVtxNearMu", &nSCPixelVtxNearMu_, "nSCPixelVtxNearMu/I");
  ntuple_->Branch("SCPixelVtxNearMu_x", &SCPixelVtxNearMu_x_, "SCPixelVtxNearMu_x[nSCPixelVtxNearMu]/F");
  ntuple_->Branch("SCPixelVtxNearMu_y", &SCPixelVtxNearMu_y_, "SCPixelVtxNearMu_y[nSCPixelVtxNearMu]/F");
  ntuple_->Branch("SCPixelVtxNearMu_z", &SCPixelVtxNearMu_z_, "SCPixelVtxNearMu_z[nSCPixelVtxNearMu]/F");
  ntuple_->Branch("SCPixelVtxNearMu_xErr", &SCPixelVtxNearMu_xErr_, "SCPixelVtxNearMu_xErr[nSCPixelVtxNearMu]/F");
  ntuple_->Branch("SCPixelVtxNearMu_yErr", &SCPixelVtxNearMu_yErr_, "SCPixelVtxNearMu_yErr[nSCPixelVtxNearMu]/F");
  ntuple_->Branch("SCPixelVtxNearMu_zErr", &SCPixelVtxNearMu_zErr_, "SCPixelVtxNearMu_zErr[nSCPixelVtxNearMu]/F");
  ntuple_->Branch("SCPixelVtxNearMu_chi2", &SCPixelVtxNearMu_chi2_, "SCPixelVtxNearMu_chi2[nSCPixelVtxNearMu]/F");
  ntuple_->Branch("SCPixelVtxNearMu_nDOF", &SCPixelVtxNearMu_nDOF_, "SCPixelVtxNearMu_nDOF[nSCPixelVtxNearMu]/I");
  ntuple_->Branch("SCPixelVtxNearMu_isValid", &SCPixelVtxNearMu_isValid_, "SCPixelVtxNearMu_isValid[nSCPixelVtxNearMu]/I");


  ntuple_->Branch("nSCMuon", &nSCMuon_, "nSCMuon/I");
  ntuple_->Branch("SCMuon_pt", &SCMuon_pt_, "SCMuon_pt[nSCMuon]/F");
  ntuple_->Branch("SCMuon_eta", &SCMuon_eta_, "SCMuon_eta[nSCMuon]/F");
  ntuple_->Branch("SCMuon_phi", &SCMuon_phi_, "SCMuon_phi[nSCMuon]/F");
  ntuple_->Branch("SCMuon_charge", &SCMuon_charge_, "SCMuon_charge[nSCMuon]/F");

  ntuple_->Branch("SCMuon_normChi2", &SCMuon_normChi2_, "SCMuon_normChi2[nSCMuon]/F");
  ntuple_->Branch("SCMuon_nTrackerLayer", &SCMuon_nTrackerLayer_, "SCMuon_nTrackerLayer[nSCMuon]/I");
  ntuple_->Branch("SCMuon_nPixelHit", &SCMuon_nPixelHit_, "SCMuon_nPixelHit[nSCMuon]/I");
  ntuple_->Branch("SCMuon_nStripHit", &SCMuon_nStripHit_, "SCMuon_nStripHit[nSCMuon]/I");
  ntuple_->Branch("SCMuon_nMuonHit", &SCMuon_nMuonHit_, "SCMuon_nMuonHit[nSCMuon]/I");
  ntuple_->Branch("SCMuon_nMatchedStation", &SCMuon_nMatchedStation_, "SCMuon_nMatchedStation[nSCMuon]/I");
  ntuple_->Branch("SCMuon_dxy", &SCMuon_dxy_, "SCMuon_dxy[nSCMuon]/F");
  ntuple_->Branch("SCMuon_dz",  &SCMuon_dz_,  "SCMuon_dz[nSCMuon]/F");
  ntuple_->Branch("SCMuon_trkIso", &SCMuon_trkIso_, "SCMuon_trkIso[nSCMuon]/F");


  ntuple_->Branch("nSCCaloJet", &nSCCaloJet_, "nSCCaloJet/I");
  ntuple_->Branch("SCCaloJet_pt",  &SCCaloJet_pt_,  "SCCaloJet_pt[nSCCaloJet]/F");
  ntuple_->Branch("SCCaloJet_eta", &SCCaloJet_eta_, "SCCaloJet_eta[nSCCaloJet]/F");
  ntuple_->Branch("SCCaloJet_phi", &SCCaloJet_phi_, "SCCaloJet_phi[nSCCaloJet]/F");
  ntuple_->Branch("SCCaloJet_m",   &SCCaloJet_m_,   "SCCaloJet_m[nSCCaloJet]/F");
  ntuple_->Branch("SCCaloJet_jetArea",        &SCCaloJet_jetArea_,           "SCCaloJet_jetArea[nSCCaloJet]/F");
  ntuple_->Branch("SCCaloJet_maxEInEmTowers", &SCCaloJet_maxEInEmTowers_,    "SCCaloJet_maxEInEmTowers[nSCCaloJet]/F");
  ntuple_->Branch("SCCaloJet_maxEInHadTowers", &SCCaloJet_maxEInHadTowers_, "SCCaloJet_maxEInHadTowers_[nSCCaloJet]/F");
  ntuple_->Branch("SCCaloJet_hadEnergyInHB",   &SCCaloJet_hadEnergyInHB_,   "SCCaloJet_hadEnergyInHB_[nSCCaloJet]/F");
  ntuple_->Branch("SCCaloJet_hadEnergyInHE",   &SCCaloJet_hadEnergyInHE_,   "SCCaloJet_hadEnergyInHE_[nSCCaloJet]/F");
  ntuple_->Branch("SCCaloJet_hadEnergyInHF",   &SCCaloJet_hadEnergyInHF_,   "SCCaloJet_hadEnergyInHF_[nSCCaloJet]/F");
  ntuple_->Branch("SCCaloJet_emEnergyInEB",    &SCCaloJet_emEnergyInEB_,    "SCCaloJet_emEnergyInEB_[nSCCaloJet]/F");
  ntuple_->Branch("SCCaloJet_emEnergyInEE",    &SCCaloJet_emEnergyInEE_,    "SCCaloJet_emEnergyInEE_[nSCCaloJet]/F");
  ntuple_->Branch("SCCaloJet_emEnergyInHF",    &SCCaloJet_emEnergyInHF_,    "SCCaloJet_emEnergyInHF_[nSCCaloJet]/F");
  ntuple_->Branch("SCCaloJet_towersArea",      &SCCaloJet_towersArea_,      "SCCaloJet_towersArea_[nSCCaloJet]/F");
  ntuple_->Branch("SCCaloJet_mvaDiscriminator",  &SCCaloJet_mvaDiscriminator_,  "SCCaloJet_mvaDiscriminator_[nSCCaloJet]/F");
  ntuple_->Branch("SCCaloJet_btagDiscriminator", &SCCaloJet_btagDiscriminator_, "SCCaloJet_btagDiscriminator_[nSCCaloJet]/F");

  ntuple_->Branch("SCCaloMET_phi", &SCCaloMET_phi_, "SCCaloMET_phi/F");
  ntuple_->Branch("SCCaloMET_pt",  &SCCaloMET_pt_,  "SCCaloMET_pt/F");
  ntuple_->Branch("SCRho",         &SCRho_,         "SCRho/F");


  // ntuple_->Branch("nL1Muon",        &nL1Muon_,        "nL1Muon/I");
  // ntuple_->Branch("L1Muon_pt",      &L1Muon_pt_,      "L1Muon_pt[nL1Muon]/F");
  // ntuple_->Branch("L1Muon_eta",     &L1Muon_eta_,     "L1Muon_eta[nL1Muon]/F");
  // ntuple_->Branch("L1Muon_phi",     &L1Muon_phi_,     "L1Muon_phi[nL1Muon]/F");
  // ntuple_->Branch("L1Muon_charge",  &L1Muon_charge_,  "L1Muon_charge[nL1Muon]/F");
  // ntuple_->Branch("L1Muon_quality", &L1Muon_quality_, "L1Muon_quality[nL1Muon]/F");

  // ntuple_->Branch("nL3MuonNoVtx",        &nL3MuonNoVtx_,        "nL3MuonNoVtx/I");
  // ntuple_->Branch("L3MuonNoVtx_pt",      &L3MuonNoVtx_pt_,      "L3MuonNoVtx_pt[nL3MuonNoVtx]/F");
  // ntuple_->Branch("L3MuonNoVtx_eta",     &L3MuonNoVtx_eta_,     "L3MuonNoVtx_eta[nL3MuonNoVtx]/F");
  // ntuple_->Branch("L3MuonNoVtx_phi",     &L3MuonNoVtx_phi_,     "L3MuonNoVtx_phi[nL3MuonNoVtx]/F");
}

void ScoutingDataTreeProducer::Fill_HLT(const edm::Event &iEvent)
{
  edm::Handle<edm::TriggerResults>  h_triggerResults;
  iEvent.getByToken(t_triggerResults_, h_triggerResults);

  // -- save trigger bit
  edm::TriggerNames triggerNames = iEvent.triggerNames(*h_triggerResults);
  for(unsigned int i_trig=0; i_trig<triggerNames.size(); ++i_trig)
  {
    if( h_triggerResults->accept(i_trig) ) // -- if accepted,
      SetTrue_HLTBitInfo(triggerNames.triggerName(i_trig));

  } // -- end of iteration over all trigger names -- //

  // -- save trigger object (if available)
  // edm::Handle<trigger::TriggerEvent> h_triggerEvent;
  // if( iEvent.getByToken(t_triggerEvent_, h_triggerEvent) ) // -- run only when TriggerEvent content is available (e.g. RAW tier (Scouting data) doesn't have TriggerEvent so below lines will not run. But it will run for usual data or MC in AOD format)
  // {
  //   const trigger::size_type nFilter(h_triggerEvent->sizeFilters());
  //   for( trigger::size_type i_filter=0; i_filter<nFilter; i_filter++)
  //   {
  //     std::string filterName = h_triggerEvent->filterTag(i_filter).encode();

  //     if( SavedFilterCondition(filterName) )
  //     {
  //       trigger::Keys objectKeys = h_triggerEvent->filterKeys(i_filter);
  //       const trigger::TriggerObjectCollection& triggerObjects(h_triggerEvent->getObjects());

  //       for( trigger::size_type i_key=0; i_key<objectKeys.size(); i_key++)
  //       {
  //         trigger::size_type objKey = objectKeys.at(i_key);
  //         const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);

  //         vec_filterName_.push_back( filterName );
  //         vec_HLTObj_pt_.push_back( triggerObj.pt() );
  //         vec_HLTObj_eta_.push_back( triggerObj.eta() );
  //         vec_HLTObj_phi_.push_back( triggerObj.phi() );
  //       }
  //     } // -- end of if( muon filters )-- //
  //   } // -- end of filter iteration -- //
  // } // -- end of if(token is available)
 
}

void ScoutingDataTreeProducer::SetTrue_HLTBitInfo(const std::string &pathName) {
  // cout << "[SetTrue_HLTBitInfo] pathName = " << pathName << endl;

  // -- check up to "_v": without this, the other triggers can be mixed (e.g. IsoMu24_XX_v triggers -> can be added in IsoMu24)
  if( pathName.find("DST_DoubleMu1_noVtx_CaloScouting_v") != std::string::npos )            DST_DoubleMu1_noVtx_CaloScouting_ = true;
  if( pathName.find("DST_DoubleMu3_noVtx_CaloScouting_Monitoring_v") != std::string::npos ) DST_DoubleMu3_noVtx_CaloScouting_Monitoring_ = true;
  if( pathName.find("DST_DoubleMu3_noVtx_CaloScouting_v") != std::string::npos )            DST_DoubleMu3_noVtx_CaloScouting_ = true;
  if( pathName.find("HLT_IsoMu24_v") != std::string::npos )                                 HLT_IsoMu24_ = true;
}

void ScoutingDataTreeProducer::Fill_SCDimuonVtx( const edm::Event& iEvent )
{
  Handle<std::vector<ScoutingVertex> > h_SCDimuonVtx;
  iEvent.getByToken(t_SCDimuonVtx_, h_SCDimuonVtx);

  int _nSCDimuonVtx = 0;
  if( h_SCDimuonVtx.isValid() )
  {
    // cout << "h_SCDimuonVtx->size() = " << h_SCDimuonVtx->size() << endl;

    for(unsigned int i_vtx=0; i_vtx<h_SCDimuonVtx->size(); ++i_vtx)
    {
      const ScoutingVertex &SCDimuonVtx = (*h_SCDimuonVtx)[i_vtx];

      SCDimuonVtx_x_[i_vtx]    = SCDimuonVtx.x();
      SCDimuonVtx_y_[i_vtx]    = SCDimuonVtx.y();
      SCDimuonVtx_z_[i_vtx]    = SCDimuonVtx.z();
      SCDimuonVtx_xErr_[i_vtx]    = SCDimuonVtx.zError();
      SCDimuonVtx_yErr_[i_vtx]    = SCDimuonVtx.yError();
      SCDimuonVtx_zErr_[i_vtx]    = SCDimuonVtx.zError();
      SCDimuonVtx_chi2_[i_vtx]    = SCDimuonVtx.chi2();
      SCDimuonVtx_nDOF_[i_vtx]    = SCDimuonVtx.ndof();
      SCDimuonVtx_isValid_[i_vtx] = SCDimuonVtx.isValidVtx();

      int index1_mu, index2_mu;
      GetMuonIndex_AssociatedToVertex(iEvent, vtx, i_vtx, index1_mu, index2_mu);
      SCDimuonVtx_muonIndex1_[i_vtx] = index1_mu;
      SCDimuonVtx_muonIndex2_[i_vtx] = index2_mu;

      _nSCDimuonVtx++;
    }

    nSCDimuonVtx_ = _nSCDimuonVtx;
  }

}

void ScoutingDataTreeProducer::Fill_SCPixelVtx( const edm::Event& iEvent )
{
  Handle<std::vector<ScoutingVertex> > h_SCPixelVtx;
  iEvent.getByToken(t_SCPixelVtx_, h_SCPixelVtx);

  int _nSCPixelVtx = 0;
  if( h_SCPixelVtx.isValid() )
  {
    for(unsigned int i_vtx=0; i_vtx<h_SCPixelVtx->size(); ++i_vtx)
    {
      const ScoutingVertex &scoutingVertex = (*h_SCPixelVtx)[i_vtx];

      SCPixelVtx_x_[i_vtx]    = scoutingVertex.x();
      SCPixelVtx_y_[i_vtx]    = scoutingVertex.y();
      SCPixelVtx_z_[i_vtx]    = scoutingVertex.z();
      SCPixelVtx_xErr_[i_vtx]    = scoutingVertex.zError();
      SCPixelVtx_yErr_[i_vtx]    = scoutingVertex.yError();
      SCPixelVtx_zErr_[i_vtx]    = scoutingVertex.zError();
      SCPixelVtx_chi2_[i_vtx]    = scoutingVertex.chi2();
      SCPixelVtx_nDOF_[i_vtx]    = scoutingVertex.ndof();
      SCPixelVtx_isValid_[i_vtx] = scoutingVertex.isValidVtx();

      _nSCPixelVtx++;
    }

    nSCPixelVtx_ = _nSCPixelVtx;
  }

}

void ScoutingDataTreeProducer::Fill_SCPixelVtxNearMu( const edm::Event& iEvent )
{
  Handle<std::vector<ScoutingVertex> > h_SCPixelVtxNearMu;
  iEvent.getByToken(t_SCPixelVtxNearMu_, h_SCPixelVtxNearMu);

  int _nSCPixelVtxNearMu = 0;
  if( h_SCPixelVtxNearMu.isValid() )
  {
    for(unsigned int i_vtx=0; i_vtx<h_SCPixelVtxNearMu->size(); ++i_vtx)
    {
      const ScoutingVertex &scoutingVertex = (*h_SCPixelVtxNearMu)[i_vtx];

      SCPixelVtxNearMu_x_[i_vtx]    = scoutingVertex.x();
      SCPixelVtxNearMu_y_[i_vtx]    = scoutingVertex.y();
      SCPixelVtxNearMu_z_[i_vtx]    = scoutingVertex.z();
      SCPixelVtxNearMu_xErr_[i_vtx]    = scoutingVertex.zError();
      SCPixelVtxNearMu_yErr_[i_vtx]    = scoutingVertex.yError();
      SCPixelVtxNearMu_zErr_[i_vtx]    = scoutingVertex.zError();
      SCPixelVtxNearMu_chi2_[i_vtx]    = scoutingVertex.chi2();
      SCPixelVtxNearMu_nDOF_[i_vtx]    = scoutingVertex.ndof();
      SCPixelVtxNearMu_isValid_[i_vtx] = scoutingVertex.isValidVtx();

      _nSCPixelVtxNearMu++;
    }

    nSCPixelVtxNearMu_ = _nSCPixelVtxNearMu;
  }

}

void ScoutingDataTreeProducer::Fill_SCMuon( const edm::Event& iEvent )
{
  Handle<std::vector<ScoutingMuon> > h_SCMuon;
  iEvent.getByToken(t_SCMuon_, h_SCMuon);

  int _nSCMuon = 0;
  if( h_SCMuon.isValid() )
  {
    // cout << "h_SCMuon->size() = " << h_SCMuon->size() << endl;

    for(unsigned int i_mu=0; i_mu<h_SCMuon->size(); ++i_mu)
    {
      const ScoutingMuon& muon = (*h_SCMuon)[i_mu];

      SCMuon_pt_[i_mu]  = muon.pt();
      SCMuon_eta_[i_mu] = muon.eta();
      SCMuon_phi_[i_mu] = muon.phi();
      SCMuon_charge_[i_mu] = muon.charge();

      SCMuon_nPixelHit_[i_mu]       = muon.nValidPixelHits();
      SCMuon_nStripHit_[i_mu]       = muon.nValidStripHits();
      SCMuon_nTrackerLayer_[i_mu]   = muon.nTrackerLayersWithMeasurement();
      SCMuon_nMuonHit_[i_mu]        = muon.nValidMuonHits();
      SCMuon_nMatchedStation_[i_mu] = muon.nMatchedStations();
      SCMuon_normChi2_[i_mu]        = muon.ndof() > 0. ? muon.chi2() / muon.ndof() : 1e4;

      SCMuon_dxy_[i_mu] = muon.dxy();
      SCMuon_dz_[i_mu]  = muon.dz();

      SCMuon_trkIso_[i_mu] = muon.trackIso();      

      // cout << "[Scouting muon: isolation] (ECAL, HCAL) = (" << muon.ecalIso() << ", " << muon.hcalIso() << ")" << endl; 

      _nSCMuon++;
    }

    nSCMuon_ = _nSCMuon;
  }

}

void ScoutingDataTreeProducer::Fill_SCCaloJet( const edm::Event& iEvent )
{
  Handle<std::vector<ScoutingCaloJet> > h_SCCaloJet;
  iEvent.getByToken(t_SCCaloJet_, h_SCCaloJet);

  int _nSCCaloJet = 0;
  if( h_SCCaloJet.isValid() )
  {
    for(unsigned int i_jet=0; i_jet<h_SCCaloJet->size(); ++i_jet)
    {
      const ScoutingCaloJet& caloJet = (*h_SCCaloJet)[i_jet];

      SCCaloJet_pt_[i_jet]   = caloJet.pt();
      SCCaloJet_eta_[i_jet]  = caloJet.eta();
      SCCaloJet_phi_[i_jet]  = caloJet.phi();
      SCCaloJet_m_[i_jet]    = caloJet.m();

      SCCaloJet_jetArea_[i_jet]         = caloJet.jetArea();
      SCCaloJet_maxEInEmTowers_[i_jet]  = caloJet.maxEInEmTowers();
      SCCaloJet_maxEInHadTowers_[i_jet] = caloJet.maxEInHadTowers();
      SCCaloJet_hadEnergyInHB_[i_jet]   = caloJet.hadEnergyInHB();
      SCCaloJet_hadEnergyInHE_[i_jet]   = caloJet.hadEnergyInHE();
      SCCaloJet_hadEnergyInHF_[i_jet]   = caloJet.hadEnergyInHF();
      SCCaloJet_emEnergyInEB_[i_jet]    = caloJet.emEnergyInEB();
      SCCaloJet_emEnergyInEE_[i_jet]    = caloJet.emEnergyInEE();
      SCCaloJet_emEnergyInHF_[i_jet]    = caloJet.emEnergyInHF();
      SCCaloJet_towersArea_[i_jet]      = caloJet.towersArea();

      SCCaloJet_mvaDiscriminator_[i_jet]  = caloJet.mvaDiscriminator();
      SCCaloJet_btagDiscriminator_[i_jet] = caloJet.btagDiscriminator();

      _nSCCaloJet++;
    }

    nSCCaloJet_ = _nSCCaloJet;
  }

}

void ScoutingDataTreeProducer::Fill_L1( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  // -- seed information
  L1GtUtils_->retrieveL1(iEvent, iSetup, t_globalAlgBlk_);

  Set_L1BitAndPrescaleInfo();

  // -- L1 muon information
  // edm::Handle<l1t::MuonBxCollection> h_L1Muon;
  // if( iEvent.getByToken(t_L1Muon_, h_L1Muon) )
  // {
  //   int _nL1Muon = 0;
  //   for(int ibx = h_L1Muon->getFirstBX(); ibx<=h_L1Muon->getLastBX(); ++ibx)
  //   {
  //     if(ibx != 0) continue; // -- only take when ibx == 0 -- //
  //     for(auto it=h_L1Muon->begin(ibx); it!=h_L1Muon->end(ibx); it++)
  //     {
  //       l1t::MuonRef ref_L1Mu(h_L1Muon, distance(h_L1Muon->begin(h_L1Muon->getFirstBX()), it) );

  //       L1Muon_pt_[_nL1Muon]      = ref_L1Mu->pt();
  //       L1Muon_eta_[_nL1Muon]     = ref_L1Mu->eta();
  //       L1Muon_phi_[_nL1Muon]     = ref_L1Mu->phi();
  //       L1Muon_charge_[_nL1Muon]  = ref_L1Mu->charge();
  //       L1Muon_quality_[_nL1Muon] = ref_L1Mu->hwQual();

  //       _nL1Muon++;
  //     }
  //   }
  //   nL1Muon_ = _nL1Muon;
  // }

}

void ScoutingDataTreeProducer::Set_L1BitAndPrescaleInfo() {
  L1GtUtils_->getFinalDecisionByName("L1_DoubleMu_15_7", L1_DoubleMu_15_7_);
  // L1GtUtils_->getPrescaleByName("L1_DoubleMu_15_7",      L1_DoubleMu_15_7_PS_ );

  L1GtUtils_->getFinalDecisionByName("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", L1_DoubleMu4p5_SQ_OS_dR_Max1p2_);
  // L1GtUtils_->getPrescaleByName("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", L1_DoubleMu4p5_SQ_OS_dR_Max1p2_PS_);
}

void ScoutingDataTreeProducer::GetMuonIndex_AssociatedToVertex(edm::Event& iEvent, const ScoutingVertex& vtx, int theVtxIndex, int& index1_mu, int& index2_mu) {
  index1_mu = -1;
  index2_mu = -1;

  edm::Handle< std::vector<ScoutingMuon> > h_SCMuon;
  iEvent.getByToken( t_SCMuon_, h_SCMuon );

  for(unsigned int i_mu=0; i_mu<h_SCMuon->size(); ++i_mu) {
    const auto& mu = (*h_SCMuon)[i_mu];

    std::vector<int> vec_vtxIndex = mu.vtxIndx();

    for(int i_vtx : vec_vtxIndex) {
      if( theVtxIndex == i_vtx ) {
        if( index1_mu == -1 )      index1_mu = i_mu;
        else if( index2_mu == -1 ) index2_mu = i_mu;
        else
          std::cout << "(Index1_mu, index2_mu) = (" << index1_mu << ", " << index2_mu << "): already filled --> the " << i_mu << "th muon will be ignored" << std::endl;

        if( index1_mu != -1 && index2_mu != -1 ) break; // -- to speed up
      }
    } // -- iteration over vertex index associated with the given muon
  } // -- iteration over muons  
}

void ScoutingDataTreeProducer::endJob() {}
void ScoutingDataTreeProducer::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {}
void ScoutingDataTreeProducer::endRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {}

DEFINE_FWK_MODULE(ScoutingDataTreeProducer);
