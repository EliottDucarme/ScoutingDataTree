// -- Flat tree producer with the nanoAOD-like format for the scouting data
// -- based on https://github.com/KyeongPil-Lee/DYScouting/tree/master/TreeProducer
// -- author: Kyeongpil Lee (ULB)
// -- update for RunIII by Eliott Ducarme (ULB)

// -- system include files
#include <memory>

// -- frameworks
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// -- triggers
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"


// -- scouting
#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"

// -- others
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

// -- ROOT includes
#include <TTree.h>
#include <TLorentzVector.h>

using namespace std;
using namespace reco;
using namespace edm;

class ScoutingDataTreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>{
public:
  explicit ScoutingDataTreeProducer(const edm::ParameterSet&);
  ~ScoutingDataTreeProducer();

  // static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
  void analyze( const edm::Event&, const edm::EventSetup& );
  void beginJob();
  void endJob();
  void beginRun( const edm::Run&, const edm::EventSetup& );
  void endRun( const edm::Run&, const edm::EventSetup& );

  void Init();
  void Make_Branch();

  void Fill_L1( const edm::Event&, const edm::EventSetup& );
  void Fill_HLT( const edm::Event & );
  // void Fill_L3MuonNoVtx( const edm::Event& );
  void Fill_SCPrimaryVtx( const edm::Event& );
  void Fill_SCDisplacedVtx( const edm::Event& );
  void Fill_SCMuon( const edm::Event& );
  void Fill_SCPFJet( const edm::Event& );

  void SetTrue_HLTBitInfo( const std::string& );
  void Set_L1BitAndPrescaleInfo();

  void GetMuonIndex_AssociatedToVertex(const edm::Event& iEvent, const Run3ScoutingVertex& vtx, int theVtxIndex, int& index1_mu, int& index2_mu);

  const edm::InputTag triggerResultsTag;

  // -- tokens
  const edm::EDGetTokenT< l1t::MuonBxCollection >       t_L1Muon_;
  const edm::EDGetTokenT< BXVector<GlobalAlgBlk> >      t_globalAlgBlk_;
  const edm::EDGetTokenT< edm::TriggerResults >         t_triggerResults_;
  //const  edm::EDGetTokenT< trigger::TriggerEvent >       t_triggerEvent_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> >  t_SCPrimaryVtx_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> >  t_SCDisplacedVtx_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon> >    t_SCMuon_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPFJet> > t_SCPFJet_;
  const edm::EDGetTokenT< double >  t_SCMETPhi_;
  const edm::EDGetTokenT< double >  t_SCMETPt_;
  const edm::EDGetTokenT< double >  t_SCRho_;



  // -- variable for L1 information
  edm::InputTag                algInputTag_;
  edm::InputTag                extInputTag_;
  edm::EDGetToken              algToken_;
  std::vector<bool>            l1Result_;
  std::unique_ptr<l1t::L1TGlobalUtil>   L1GtUtils_;

  // -- IterL3MuonCandidateNoVtx object
  // edm::EDGetTokenT< std::vector<pat::TriggerObjectStandAlone> > t_trigObj_L3MuonNoVtx_;

  // -- debugging switch
  bool debug_ = false;

  // -- tree
  TTree* ntuple_;

  // -- array size
  static const int arrSize_ = 2000;

  // --  event informations
  unsigned int run_;
  unsigned int luminosityBlock_;
  unsigned long long event_;

  // -- L1 flags
  bool   L1_DoubleMu_15_7_;
  bool   L1_DoubleMu4p5_SQ_OS_dR_Max1p2_;
  bool   L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7_;
  bool   L1_DoubleMu8_SQ_;

  // -- HLT flags
  bool   DST_Physics_;
  bool   DST_Run3_PFScoutingPixelTracking_;
  bool   DST_HLTMuon_Run3_PFScoutingPixelTracking_;
  bool   HLT_IsoMu24_;

  // -- HLT objects
  // -- use float (instead of double) to sync. with nanoAOD MC
  unsigned int   nTrigObj_;
  float TrigObj_pt[arrSize_];
  float TrigObj_eta[arrSize_];
  float TrigObj_phi[arrSize_];

  // -- primary vertex information
  unsigned int   nSCPrimaryVtx_;
  float SCPrimaryVtx_x_[arrSize_];
  float SCPrimaryVtx_y_[arrSize_];
  float SCPrimaryVtx_z_[arrSize_];
  float SCPrimaryVtx_xErr_[arrSize_];
  float SCPrimaryVtx_yErr_[arrSize_];
  float SCPrimaryVtx_zErr_[arrSize_];
  float SCPrimaryVtx_chi2_[arrSize_];
  int   SCPrimaryVtx_nDOF_[arrSize_];
  int   SCPrimaryVtx_muonIndex1_[arrSize_];
  int   SCPrimaryVtx_muonIndex2_[arrSize_];
  bool  SCPrimaryVtx_isValid_[arrSize_];

  // -- displaced vertex information
  unsigned int nSCDisplacedVtx_;
  float SCDisplacedVtx_x_[arrSize_];
  float SCDisplacedVtx_y_[arrSize_];
  float SCDisplacedVtx_z_[arrSize_];
  float SCDisplacedVtx_xErr_[arrSize_];
  float SCDisplacedVtx_yErr_[arrSize_];
  float SCDisplacedVtx_zErr_[arrSize_];
  float SCDisplacedVtx_chi2_[arrSize_];
  int   SCDisplacedVtx_nDOF_[arrSize_];
  int   SCDisplacedVtx_isValid_[arrSize_];

  // -- muon information (scouting)
  unsigned int nSCMuon_;
  float SCMuon_pt_[arrSize_];
  float SCMuon_eta_[arrSize_];
  float SCMuon_phi_[arrSize_];
  float SCMuon_mass_[arrSize_];
  float SCMuon_charge_[arrSize_];

  int SCMuon_nPixelHit_[arrSize_];
  int SCMuon_nStripHit_[arrSize_];
  int SCMuon_nTrackerLayer_[arrSize_];
  int SCMuon_nMuonHit_[arrSize_];
  int SCMuon_nMatchedStation_[arrSize_];
  int SCMuon_nDOF_[arrSize_];
  float SCMuon_chi2_[arrSize_];
  float SCMuon_dxy_[arrSize_];
  float SCMuon_dz_[arrSize_];
  float SCMuon_trackIso_[arrSize_];
  float SCMuon_ecalIso_[arrSize_];
  float SCMuon_hcalIso_[arrSize_];

  float SCMuon_isGlobal_[arrSize_];
  float SCMuon_isTracker_[arrSize_];

  // -- jet information
  unsigned int nSCPFJet_;
  double SCPFJet_pt_[arrSize_];
  double SCPFJet_eta_[arrSize_];
  double SCPFJet_phi_[arrSize_];
  double SCPFJet_m_[arrSize_];
  double SCPFJet_jetArea_[arrSize_];
  double SCPFJet_maxEInEmTowers_[arrSize_];
  double SCPFJet_maxEInHadTowers_[arrSize_];
  double SCPFJet_hadEnergyInHB_[arrSize_];
  double SCPFJet_hadEnergyInHE_[arrSize_];
  double SCPFJet_hadEnergyInHF_[arrSize_];
  double SCPFJet_emEnergyInEB_[arrSize_];
  double SCPFJet_emEnergyInEE_[arrSize_];
  double SCPFJet_emEnergyInHF_[arrSize_];
  double SCPFJet_towersArea_[arrSize_];
  double SCPFJet_mvaDiscriminator_[arrSize_];
  double SCPFJet_btagDiscriminator_[arrSize_];

  // -- MET, Calo info.
  double SCMET_phi_;
  double SCMET_pt_;
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
triggerResultsTag       (iConfig.getUntrackedParameter<edm::InputTag>("triggerResults")),
t_L1Muon_               ( consumes< l1t::MuonBxCollection  >         (iConfig.getUntrackedParameter<edm::InputTag>("L1Muon")) ),
t_globalAlgBlk_         ( consumes< BXVector< GlobalAlgBlk > >       (iConfig.getUntrackedParameter<edm::InputTag>("globalAlgBlk")) ),
t_triggerResults_       ( consumes< edm::TriggerResults >            (iConfig.getUntrackedParameter<edm::InputTag>("triggerResults")) ),
t_SCPrimaryVtx_         ( consumes< std::vector<Run3ScoutingVertex> >    (iConfig.getUntrackedParameter<edm::InputTag>("SCPrimaryVtx")) ),
t_SCDisplacedVtx_       ( consumes< std::vector<Run3ScoutingVertex> >    (iConfig.getUntrackedParameter<edm::InputTag>("SCDisplacedVtx")) ),
t_SCMuon_               ( consumes< std::vector<Run3ScoutingMuon> >      (iConfig.getUntrackedParameter<edm::InputTag>("SCMuon")) ),
t_SCPFJet_            ( consumes< std::vector<Run3ScoutingPFJet> >   (iConfig.getUntrackedParameter<edm::InputTag>("SCPFJet")) ),
t_SCMETPhi_         ( consumes< double >                         (iConfig.getUntrackedParameter<edm::InputTag>("SCMETPhi")) ),
t_SCMETPt_          ( consumes< double >                         (iConfig.getUntrackedParameter<edm::InputTag>("SCMETPt")) ),
t_SCRho_                ( consumes< double >                         (iConfig.getUntrackedParameter<edm::InputTag>("SCRho")) )
// t_trigObj_L3MuonNoVtx_ ( consumes< std::vector<pat::TriggerObjectStandAlone> >   (iConfig.getParameter<edm::InputTag>("triggerObject_L3MuonNoVtx")) ),
{
  usesResource("TFileService");
  algInputTag_ = iConfig.getUntrackedParameter<edm::InputTag>("AlgInputTag");
  extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
  algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
  L1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
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
  Fill_SCPrimaryVtx(iEvent);
  Fill_SCDisplacedVtx(iEvent);
  Fill_SCMuon(iEvent);
  Fill_SCPFJet(iEvent);

  // -- MET variables
  edm::Handle<double> h_SCMETPhi;
  iEvent.getByToken(t_SCMETPhi_, h_SCMETPhi);
  if( h_SCMETPhi.isValid() ) SCMET_phi_ = *(h_SCMETPhi.product());

  edm::Handle<double> h_SCMETPt;
  iEvent.getByToken(t_SCMETPt_, h_SCMETPt);
  if( h_SCMETPt.isValid() ) SCMET_pt_ = *(h_SCMETPt.product());

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
  L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7_ = false;
  L1_DoubleMu8_SQ_ = false;

  // -- HLT flags
  DST_Physics_ = false;
  DST_Run3_PFScoutingPixelTracking_ = false;
  DST_HLTMuon_Run3_PFScoutingPixelTracking_ = false;
  HLT_IsoMu24_ = false;

  // -- HLT objects
  nTrigObj_ = -999;
  for(Int_t i=0; i<arrSize_; i++) {
    TrigObj_pt[i] = -999;
    TrigObj_eta[i] = -999;
    TrigObj_phi[i] = -999;
  }


  // -- vertex information (@ HLT)
  nSCPrimaryVtx_      = -999;
  nSCDisplacedVtx_ = -999;
  for(Int_t i=0; i<arrSize_; i++) {
    SCPrimaryVtx_x_[i] = -999;
    SCPrimaryVtx_y_[i] = -999;
    SCPrimaryVtx_z_[i] = -999;
    SCPrimaryVtx_xErr_[i] = -999;
    SCPrimaryVtx_yErr_[i] = -999;
    SCPrimaryVtx_zErr_[i] = -999;
    SCPrimaryVtx_chi2_[i] = -999;
    SCPrimaryVtx_nDOF_[i] = -999;
    SCPrimaryVtx_muonIndex1_[i] = -999;
    SCPrimaryVtx_muonIndex2_[i] = -999;
    SCPrimaryVtx_isValid_[i] = 0;

    SCDisplacedVtx_x_[i] = -999;
    SCDisplacedVtx_y_[i] = -999;
    SCDisplacedVtx_z_[i] = -999;
    SCDisplacedVtx_xErr_[i] = -999;
    SCDisplacedVtx_yErr_[i] = -999;
    SCDisplacedVtx_zErr_[i] = -999;
    SCDisplacedVtx_chi2_[i] = -999;
    SCDisplacedVtx_nDOF_[i] = -999;
    SCDisplacedVtx_isValid_[i] = 0;
  }

  // -- muon information
  nSCMuon_ = -999;
  for(Int_t i=0; i<arrSize_; i++)
  {
    SCMuon_pt_[i] = -999;
    SCMuon_eta_[i] = -999;
    SCMuon_phi_[i] = -999;
    SCMuon_mass_[i] = -999;
    SCMuon_charge_[i] = -999;

    SCMuon_nPixelHit_[i] = -999;
    SCMuon_nStripHit_[i] = -999;
    SCMuon_nTrackerLayer_[i] = -999;
    SCMuon_nMuonHit_[i] = -999;
    SCMuon_nMatchedStation_[i] = -999;
    SCMuon_nDOF_[i] = -999;
    SCMuon_chi2_[i] = -999;
    SCMuon_dxy_[i] = -999;
    SCMuon_dz_[i] = -999;
    SCMuon_trackIso_[i] = -999;
    SCMuon_hcalIso_[i] = -999;
    SCMuon_ecalIso_[i] = -999;
    
    SCMuon_isGlobal_[i] = 0;
    SCMuon_isTracker_[i] = 0;
  }

  // -- jet information
  nSCPFJet_ = -999;
  for(Int_t i=0; i<arrSize_; i++)
  {
    SCPFJet_pt_[i] = -999;
    SCPFJet_eta_[i] = -999;
    SCPFJet_phi_[i] = -999;
    SCPFJet_m_[i] = -999;
    SCPFJet_jetArea_[i] = -999;
    // SCPFJet_maxEInEmTowers_[i] = -999;
    // SCPFJet_maxEInHadTowers_[i] = -999;
    // SCPFJet_hadEnergyInHB_[i] = -999;
    // SCPFJet_hadEnergyInHE_[i] = -999;
    // SCPFJet_hadEnergyInHF_[i] = -999;
    // SCPFJet_emEnergyInEB_[i] = -999;
    // SCPFJet_emEnergyInEE_[i] = -999;
    // SCPFJet_emEnergyInHF_[i] = -999;
    // SCPFJet_towersArea_[i] = -999;
    SCPFJet_mvaDiscriminator_[i] = -999;
    // SCPFJet_btagDiscriminator_[i] = -999;
  }

  SCMET_phi_ = -999;
  SCMET_pt_ = -999;
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
  ntuple_->Branch("run",&run_,"run/i");
  ntuple_->Branch("luminosityBlock",&luminosityBlock_,"luminosityBlock/i");
  ntuple_->Branch("event",&event_,"event/l"); // -- unsigned long long -- //

  ntuple_->Branch("L1_DoubleMu_15_7",               &L1_DoubleMu_15_7_,               "L1_DoubleMu_15_7/O");
  ntuple_->Branch("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2_, "L1_DoubleMu4p5_SQ_OS_dR_Max1p2/O");
  ntuple_->Branch("L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7", &L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7_, "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7/O");
  ntuple_->Branch("L1_DoubleMu8_SQ", &L1_DoubleMu8_SQ_, "L1_DoubleMu8_SQ/O");

  ntuple_->Branch("DST_HLTMuon_Run3_PFScoutingPixelTracking",            &DST_HLTMuon_Run3_PFScoutingPixelTracking_,            "DST_HLTMuon_Run3_PFScoutingPixelTracking/O");
  ntuple_->Branch("DST_Run3_PFScoutingPixelTracking",            &DST_Run3_PFScoutingPixelTracking_,            "DST_Run3_PFScoutingPixelTracking/O");
  ntuple_->Branch("DST_Physics",            &DST_Physics_,            "DST_Physics/O");

  // ntuple_->Branch("nTrigObj", &nTrigObj_, "nTrigObj/I");
  // ntuple_->Branch("TrigObj_pt",  &TrigObj_pt_, "TrigObj_pt[nTrigObj]/F");
  // ntuple_->Branch("TrigObj_eta", &TrigObj_eta_, "TrigObj_eta[nTrigObj]/F");
  // ntuple_->Branch("TrigObj_phi", &TrigObj_phi_, "TrigObj_phi[nTrigObj]/F");

  ntuple_->Branch("nSCPrimaryVtx", &nSCPrimaryVtx_, "nSCPrimaryVtx/i");
  ntuple_->Branch("SCPrimaryVtx_x", &SCPrimaryVtx_x_, "SCPrimaryVtx_x[nSCPrimaryVtx]/F");
  ntuple_->Branch("SCPrimaryVtx_y", &SCPrimaryVtx_y_, "SCPrimaryVtx_y[nSCPrimaryVtx]/F");
  ntuple_->Branch("SCPrimaryVtx_z", &SCPrimaryVtx_z_, "SCPrimaryVtx_z[nSCPrimaryVtx]/F");
  ntuple_->Branch("SCPrimaryVtx_xErr", &SCPrimaryVtx_xErr_, "SCPrimaryVtx_xErr[nSCPrimaryVtx]/F");
  ntuple_->Branch("SCPrimaryVtx_yErr", &SCPrimaryVtx_yErr_, "SCPrimaryVtx_yErr[nSCPrimaryVtx]/F");
  ntuple_->Branch("SCPrimaryVtx_zErr", &SCPrimaryVtx_zErr_, "SCPrimaryVtx_zErr[nSCPrimaryVtx]/F");
  ntuple_->Branch("SCPrimaryVtx_chi2", &SCPrimaryVtx_chi2_, "SCPrimaryVtx_chi2[nSCPrimaryVtx]/F");
  ntuple_->Branch("SCPrimaryVtx_nDOF", &SCPrimaryVtx_nDOF_, "SCPrimaryVtx_nDOF[nSCPrimaryVtx]/I");
  ntuple_->Branch("SCPrimaryVtx_muonIndex1", &SCPrimaryVtx_muonIndex1_, "SCPrimaryVtx_muonIndex1[nSCPrimaryVtx]/I");
  ntuple_->Branch("SCPrimaryVtx_muonIndex2", &SCPrimaryVtx_muonIndex2_, "SCPrimaryVtx_muonIndex2[nSCPrimaryVtx]/I");
  ntuple_->Branch("SCPrimaryVtx_isValid", &SCPrimaryVtx_isValid_, "SCPrimaryVtx_isValid[nSCPrimaryVtx]/O");

  ntuple_->Branch("nSCDisplacedVtx", &nSCDisplacedVtx_, "nSCDisplacedVtx/i");
  ntuple_->Branch("SCDisplacedVtx_x", &SCDisplacedVtx_x_, "SCDisplacedVtx_x[nSCDisplacedVtx]/F");
  ntuple_->Branch("SCDisplacedVtx_y", &SCDisplacedVtx_y_, "SCDisplacedVtx_y[nSCDisplacedVtx]/F");
  ntuple_->Branch("SCDisplacedVtx_z", &SCDisplacedVtx_z_, "SCDisplacedVtx_z[nSCDisplacedVtx]/F");
  ntuple_->Branch("SCDisplacedVtx_xErr", &SCDisplacedVtx_xErr_, "SCDisplacedVtx_xErr[nSCDisplacedVtx]/F");
  ntuple_->Branch("SCDisplacedVtx_yErr", &SCDisplacedVtx_yErr_, "SCDisplacedVtx_yErr[nSCDisplacedVtx]/F");
  ntuple_->Branch("SCDisplacedVtx_zErr", &SCDisplacedVtx_zErr_, "SCDisplacedVtx_zErr[nSCDisplacedVtx]/F");
  ntuple_->Branch("SCDisplacedVtx_chi2", &SCDisplacedVtx_chi2_, "SCDisplacedVtx_chi2[nSCDisplacedVtx]/F");
  ntuple_->Branch("SCDisplacedVtx_nDOF", &SCDisplacedVtx_nDOF_, "SCDisplacedVtx_nDOF[nSCDisplacedVtx]/I");
  ntuple_->Branch("SCDisplacedVtx_isValid", &SCDisplacedVtx_isValid_, "SCDisplacedVtx_isValid[nSCDisplacedVtx]/I");


  ntuple_->Branch("nSCMuon", &nSCMuon_, "nSCMuon/i");
  ntuple_->Branch("SCMuon_pt", &SCMuon_pt_, "SCMuon_pt[nSCMuon]/F");
  ntuple_->Branch("SCMuon_eta", &SCMuon_eta_, "SCMuon_eta[nSCMuon]/F");
  ntuple_->Branch("SCMuon_phi", &SCMuon_phi_, "SCMuon_phi[nSCMuon]/F");
  ntuple_->Branch("SCMuon_mass", &SCMuon_mass_, "SCMuon_mass[nSCMuon]/F");
  ntuple_->Branch("SCMuon_charge", &SCMuon_charge_, "SCMuon_charge[nSCMuon]/F");

  ntuple_->Branch("SCMuon_chi2", &SCMuon_chi2_, "SCMuon_chi2[nSCMuon]/F");
  ntuple_->Branch("SCMuon_nTrackerLayer", &SCMuon_nTrackerLayer_, "SCMuon_nTrackerLayer[nSCMuon]/I");
  ntuple_->Branch("SCMuon_nPixelHit", &SCMuon_nPixelHit_, "SCMuon_nPixelHit[nSCMuon]/I");
  ntuple_->Branch("SCMuon_nStripHit", &SCMuon_nStripHit_, "SCMuon_nStripHit[nSCMuon]/I");
  ntuple_->Branch("SCMuon_nMuonHit", &SCMuon_nMuonHit_, "SCMuon_nMuonHit[nSCMuon]/I");
  ntuple_->Branch("SCMuon_nMatchedStation", &SCMuon_nMatchedStation_, "SCMuon_nMatchedStation[nSCMuon]/I");
  ntuple_->Branch("SCMuon_nDOF", &SCMuon_nDOF_, "SCMuon_nDOF[nSCMuon]/I");
  ntuple_->Branch("SCMuon_dxy", &SCMuon_dxy_, "SCMuon_dxy[nSCMuon]/F");
  ntuple_->Branch("SCMuon_dz",  &SCMuon_dz_,  "SCMuon_dz[nSCMuon]/F");
  ntuple_->Branch("SCMuon_trackIso", &SCMuon_trackIso_, "SCMuon_trackIso[nSCMuon]/F");
  ntuple_->Branch("SCMuon_ecalIso", &SCMuon_ecalIso_, "SCMuon_ecalIso[nSCMuon]/F");
  ntuple_->Branch("SCMuon_hcalIso", &SCMuon_hcalIso_, "SCMuon_hcalIso[nSCMuon]/F");
  ntuple_->Branch("SCMuon_isGlobal", &SCMuon_isGlobal_, "SCMuon_isGlobal[nSCMuon]/O");
  ntuple_->Branch("SCMuon_isTracker", &SCMuon_isTracker_, "SCMuon_isTracker[nSCMuon]/O");


  ntuple_->Branch("nSCPFJet", &nSCPFJet_, "nSCPFJet/i");
  ntuple_->Branch("SCPFJet_pt",  &SCPFJet_pt_,  "SCPFJet_pt[nSCPFJet]/F");
  ntuple_->Branch("SCPFJet_eta", &SCPFJet_eta_, "SCPFJet_eta[nSCPFJet]/F");
  ntuple_->Branch("SCPFJet_phi", &SCPFJet_phi_, "SCPFJet_phi[nSCPFJet]/F");
  ntuple_->Branch("SCPFJet_m",   &SCPFJet_m_,   "SCPFJet_m[nSCPFJet]/F");
  ntuple_->Branch("SCPFJet_jetArea",        &SCPFJet_jetArea_,           "SCPFJet_jetArea[nSCPFJet]/F");
  // ntuple_->Branch("SCPFJet_maxEInEmTowers", &SCPFJet_maxEInEmTowers_,    "SCPFJet_maxEInEmTowers[nSCPFJet]/F");
  // ntuple_->Branch("SCPFJet_maxEInHadTowers", &SCPFJet_maxEInHadTowers_, "SCPFJet_maxEInHadTowers_[nSCPFJet]/F");
  // ntuple_->Branch("SCPFJet_hadEnergyInHB",   &SCPFJet_hadEnergyInHB_,   "SCPFJet_hadEnergyInHB_[nSCPFJet]/F");
  // ntuple_->Branch("SCPFJet_hadEnergyInHE",   &SCPFJet_hadEnergyInHE_,   "SCPFJet_hadEnergyInHE_[nSCPFJet]/F");
  // ntuple_->Branch("SCPFJet_hadEnergyInHF",   &SCPFJet_hadEnergyInHF_,   "SCPFJet_hadEnergyInHF_[nSCPFJet]/F");
  // ntuple_->Branch("SCPFJet_emEnergyInEB",    &SCPFJet_emEnergyInEB_,    "SCPFJet_emEnergyInEB_[nSCPFJet]/F");
  // ntuple_->Branch("SCPFJet_emEnergyInEE",    &SCPFJet_emEnergyInEE_,    "SCPFJet_emEnergyInEE_[nSCPFJet]/F");
  // ntuple_->Branch("SCPFJet_emEnergyInHF",    &SCPFJet_emEnergyInHF_,    "SCPFJet_emEnergyInHF_[nSCPFJet]/F");
  // ntuple_->Branch("SCPFJet_towersArea",      &SCPFJet_towersArea_,      "SCPFJet_towersArea_[nSCPFJet]/F");
  ntuple_->Branch("SCPFJet_mvaDiscriminator",  &SCPFJet_mvaDiscriminator_,  "SCPFJet_mvaDiscriminator_[nSCPFJet]/F");
  // ntuple_->Branch("SCPFJet_btagDiscriminator", &SCPFJet_btagDiscriminator_, "SCPFJet_btagDiscriminator_[nSCPFJet]/F");

  ntuple_->Branch("SCMET_phi", &SCMET_phi_, "SCMET_phi/F");
  ntuple_->Branch("SCMET_pt",  &SCMET_pt_,  "SCMET_pt/F");
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
  if( pathName.find("DST_Physics_v") != std::string::npos )                                 DST_Physics_ = true;
  if( pathName.find("DST_Run3_PFScoutingPixelTracking_v") != std::string::npos ) DST_Run3_PFScoutingPixelTracking_ = true;
  if( pathName.find("DST_HLTMuon_Run3_PFScoutingPixelTracking_v") != std::string::npos )    DST_HLTMuon_Run3_PFScoutingPixelTracking_ = true;
  if( pathName.find("HLT_IsoMu24_v") != std::string::npos )                                 HLT_IsoMu24_ = true;
}

void ScoutingDataTreeProducer::Fill_SCPrimaryVtx( const edm::Event& iEvent )
{
  Handle<std::vector<Run3ScoutingVertex> > h_SCPrimaryVtx;
  iEvent.getByToken(t_SCPrimaryVtx_, h_SCPrimaryVtx);

  int _nSCPrimaryVtx = 0;
  if( h_SCPrimaryVtx.isValid() )
  {
    // cout << "h_SCPrimaryVtx->size() = " << h_SCPrimaryVtx->size() << endl;

    for(unsigned int i_vtx=0; i_vtx<h_SCPrimaryVtx->size(); ++i_vtx)
    {
      const Run3ScoutingVertex &SCPrimaryVtx = (*h_SCPrimaryVtx)[i_vtx];

      SCPrimaryVtx_x_[i_vtx]    = SCPrimaryVtx.x();
      SCPrimaryVtx_y_[i_vtx]    = SCPrimaryVtx.y();
      SCPrimaryVtx_z_[i_vtx]    = SCPrimaryVtx.z();
      SCPrimaryVtx_xErr_[i_vtx]    = SCPrimaryVtx.zError();
      SCPrimaryVtx_yErr_[i_vtx]    = SCPrimaryVtx.yError();
      SCPrimaryVtx_zErr_[i_vtx]    = SCPrimaryVtx.zError();
      SCPrimaryVtx_chi2_[i_vtx]    = SCPrimaryVtx.chi2();
      SCPrimaryVtx_nDOF_[i_vtx]    = SCPrimaryVtx.ndof();
      SCPrimaryVtx_isValid_[i_vtx] = SCPrimaryVtx.isValidVtx();

      int index1_mu, index2_mu;
      GetMuonIndex_AssociatedToVertex(iEvent, SCPrimaryVtx, i_vtx, index1_mu, index2_mu);
      SCPrimaryVtx_muonIndex1_[i_vtx] = index1_mu;
      SCPrimaryVtx_muonIndex2_[i_vtx] = index2_mu;

      _nSCPrimaryVtx++;
    }

    nSCPrimaryVtx_ = _nSCPrimaryVtx;
  }

}

void ScoutingDataTreeProducer::Fill_SCDisplacedVtx( const edm::Event& iEvent )
{
  Handle<std::vector<Run3ScoutingVertex> > h_SCDisplacedVtx;
  iEvent.getByToken(t_SCDisplacedVtx_, h_SCDisplacedVtx);

  int _nSCDisplacedVtx = 0;
  if( h_SCDisplacedVtx.isValid() )
  {
    for(unsigned int i_vtx=0; i_vtx<h_SCDisplacedVtx->size(); ++i_vtx)
    {
      const Run3ScoutingVertex &Run3ScoutingVertex = (*h_SCDisplacedVtx)[i_vtx];

      SCDisplacedVtx_x_[i_vtx]    = Run3ScoutingVertex.x();
      SCDisplacedVtx_y_[i_vtx]    = Run3ScoutingVertex.y();
      SCDisplacedVtx_z_[i_vtx]    = Run3ScoutingVertex.z();
      SCDisplacedVtx_xErr_[i_vtx]    = Run3ScoutingVertex.zError();
      SCDisplacedVtx_yErr_[i_vtx]    = Run3ScoutingVertex.yError();
      SCDisplacedVtx_zErr_[i_vtx]    = Run3ScoutingVertex.zError();
      SCDisplacedVtx_chi2_[i_vtx]    = Run3ScoutingVertex.chi2();
      SCDisplacedVtx_nDOF_[i_vtx]    = Run3ScoutingVertex.ndof();
      SCDisplacedVtx_isValid_[i_vtx] = Run3ScoutingVertex.isValidVtx();

      _nSCDisplacedVtx++;
    }

    nSCDisplacedVtx_ = _nSCDisplacedVtx;
  }

}

void ScoutingDataTreeProducer::Fill_SCMuon( const edm::Event& iEvent )
{
  Handle<std::vector<Run3ScoutingMuon> > h_SCMuon;
  iEvent.getByToken(t_SCMuon_, h_SCMuon);

  int _nSCMuon = 0;
  // if( h_SCMuon.isGlobalMuon() || h_SCMuon.isTrackerMuon() )
  // {
    // cout << "h_SCMuon->size() = " << h_SCMuon->size() << endl;

    for(unsigned int i_mu=0; i_mu<h_SCMuon->size(); ++i_mu)
    {
      const Run3ScoutingMuon& muon = (*h_SCMuon)[i_mu];

      SCMuon_mass_[i_mu] = muon.m();
      SCMuon_pt_[i_mu]  = muon.pt();
      SCMuon_eta_[i_mu] = muon.eta();
      SCMuon_phi_[i_mu] = muon.phi();
      SCMuon_charge_[i_mu] = muon.charge();
      SCMuon_isGlobal_[i_mu] = muon.isGlobalMuon();
      SCMuon_isTracker_[i_mu] = muon.isTrackerMuon();

      SCMuon_trackIso_[i_mu] = muon.trackIso();
      SCMuon_ecalIso_[i_mu] = muon.ecalIso();
      SCMuon_hcalIso_[i_mu] = muon.hcalIso();

      SCMuon_nPixelHit_[i_mu]       = muon.nValidPixelHits();
      SCMuon_nStripHit_[i_mu]       = muon.nValidStripHits();
      SCMuon_nTrackerLayer_[i_mu]   = muon.nTrackerLayersWithMeasurement();
      SCMuon_nMuonHit_[i_mu]        = muon.nValidRecoMuonHits();
      SCMuon_nMatchedStation_[i_mu] = muon.nRecoMuonMatchedStations();
      SCMuon_chi2_[i_mu]            = muon.trk_chi2();
      SCMuon_nDOF_[i_mu]            = muon.trk_ndof();

      SCMuon_dxy_[i_mu] = muon.trk_dxy();
      SCMuon_dz_[i_mu]  = muon.trk_dz();

      // cout << "[Scouting muon: isolation] (ECAL, HCAL) = (" << muon.ecalIso() << ", " << muon.hcalIso() << ")" << endl; 

      _nSCMuon++;
    }

    nSCMuon_ = _nSCMuon;
  // }

}

void ScoutingDataTreeProducer::Fill_SCPFJet( const edm::Event& iEvent )
{
  Handle<vector<Run3ScoutingPFJet> > h_SCPFJet;
  iEvent.getByToken(t_SCPFJet_, h_SCPFJet);

  int _nSCPFJet = 0;
  // if( h_SCPFJet.isValid() )
  // {
    for(unsigned int i_jet=0; i_jet<h_SCPFJet->size(); ++i_jet)
    {
      const Run3ScoutingPFJet& PFJet = (*h_SCPFJet)[i_jet];

      SCPFJet_pt_[i_jet]   = PFJet.pt();
      SCPFJet_eta_[i_jet]  = PFJet.eta();
      SCPFJet_phi_[i_jet]  = PFJet.phi();
      SCPFJet_m_[i_jet]    = PFJet.m();

      SCPFJet_jetArea_[i_jet]         = PFJet.jetArea();
      // SCPFJet_maxEInEmTowers_[i_jet]  = PFJet.maxEInEmTowers();
      // SCPFJet_maxEInHadTowers_[i_jet] = PFJet.maxEInHadTowers();
      // SCPFJet_hadEnergyInHB_[i_jet]   = PFJet.hadEnergyInHB();
      // SCPFJet_hadEnergyInHE_[i_jet]   = PFJet.hadEnergyInHE();
      // SCPFJet_hadEnergyInHF_[i_jet]   = PFJet.hadEnergyInHF();
      // SCPFJet_emEnergyInEB_[i_jet]    = PFJet.emEnergyInEB();
      // SCPFJet_emEnergyInEE_[i_jet]    = PFJet.emEnergyInEE();
      // SCPFJet_emEnergyInHF_[i_jet]    = PFJet.emEnergyInHF();
      // SCPFJet_towersArea_[i_jet]      = PFJet.towersArea();

      SCPFJet_mvaDiscriminator_[i_jet]  = PFJet.mvaDiscriminator();
      // SCPFJet_btagDiscriminator_[i_jet] = PFJet.btagDiscriminator();

      _nSCPFJet++;
    }

    nSCPFJet_ = _nSCPFJet;
  // }

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

  L1GtUtils_->getFinalDecisionByName("L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7", L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7_);
  L1GtUtils_->getFinalDecisionByName("L1_DoubleMu8_SQ", L1_DoubleMu8_SQ_);
  
}

void ScoutingDataTreeProducer::GetMuonIndex_AssociatedToVertex(const edm::Event& iEvent, const Run3ScoutingVertex& vtx, int theVtxIndex, int& index1_mu, int& index2_mu) {
  index1_mu = -1;
  index2_mu = -1;

  edm::Handle< std::vector<Run3ScoutingMuon> > h_SCMuon;
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
