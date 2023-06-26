#ifndef TFLATMCOUT
#define TFLATMCOUT

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "ReturnRes.hh"
#include "TDataProcess.h"
#include "TFile.h"
#include "THyphiAttributes.h"
#include "TRandom3.h"
#include "TTree.h"

class DataML
{
protected:
  const THyphiAttributes& att;

public:
  TTree* out_tree;
  TTree* out_tree_signal;
  TTree* out_tree_background;

  int SplitTree;

  DataML(const THyphiAttributes& att_) : att(att_), out_tree(nullptr), out_tree_signal(nullptr), out_tree_background(nullptr) {}
  DataML(const THyphiAttributes& att_, TTree* outT) :
                                    att(att_), out_tree(outT), out_tree_signal(nullptr), out_tree_background(nullptr) {}
  DataML(const THyphiAttributes& att_, TTree* outT, TTree* outTs, TTree* outTb) :
                                    att(att_), out_tree(outT), out_tree_signal(outTs), out_tree_background(outTb) {}
  virtual ~DataML() { out_tree = nullptr; out_tree_signal = nullptr; out_tree_background = nullptr; }
  virtual void FillEvent(FullRecoEvent& RecoEvent) = 0;
};

class DataML_momfit : public DataML
{
public:
  Float_t b_tx, b_ty, b_ptnx, b_ptny, b_vx, b_vy, b_vz, b_x, b_y, b_z, b_pt, b_phi, b_theta;
  Float_t a_tx, a_ty, a_ptnx, a_ptny, a_vx, a_vy, a_vz, a_x, a_y, a_z, a_pt, a_phi, a_theta;
  Float_t Dphi;
  Float_t poq, qop;
  Float_t ptoq, qopt;
  Float_t q;
  Float_t tof;
  Float_t psb_z;
  Int_t psce_psbe;
  Int_t pdg;
  DataML_momfit(const THyphiAttributes& att, TTree* outT);
  void FillEvent(FullRecoEvent& RecoEvent) override;
};

class DataML_Hyp : public DataML
{
public:
  Int_t Pattern; /// 1 = Simulation / 2 = KFParticle_real / 3 = KFParticle_cut / 4 = KFParticle / 5 = KFParticle_Mass / 6 = LorentzVector
  
  //Mother:
  Int_t PDG;
  Int_t N_Mother;
  Double32_t Chi2ndf;
  Int_t NDF;
  Double32_t MomX;
  Double32_t MomY;
  Double32_t MomZ;
  Double32_t E;
  Double32_t PrimVtx_PosX;
  Double32_t PrimVtx_PosY;
  Double32_t PrimVtx_PosZ;
  Double32_t DecayVtx_PosX;
  Double32_t DecayVtx_PosY;
  Double32_t DecayVtx_PosZ;
  Double32_t Dist_RealReconsVtx_X;
  Double32_t Dist_RealReconsVtx_Y;
  Double32_t Dist_RealReconsVtx_Z;
  Double32_t Dist_MotherPrimVtx;
  Double32_t Angle_MotherPrimVtx;
  Double32_t InvMass;
  Double32_t ErrInvMass;
  Int_t ErrGetMass;
  Double32_t LifeTime;
  Double32_t ErrLifeTime;
  Int_t ErrGetLifeTime;
  Int_t Mother_IsFromHyp; // 0-> No; 1-> Yes

  //Daughters:
  Int_t Id_Fragment;
  Double32_t MomX_Fragment;
  Double32_t MomY_Fragment;
  Double32_t MomZ_Fragment;
  Double32_t E_Fragment;
  Double32_t Chi2ndf_Fragment;
  Int_t NDF_Fragment;
  Double32_t Pvalue_Fragment;
  Double32_t Angle_MotherFragment;
  Int_t Fragment_IsFromHyp; // 0-> No; 1-> Yes
  Int_t RealPDG_Fragment;
  
  Int_t Id_Pion;
  Double32_t MomX_Pion;
  Double32_t MomY_Pion;
  Double32_t MomZ_Pion;
  Double32_t E_Pion;
  Double32_t Chi2ndf_Pion;
  Int_t NDF_Pion;
  Double32_t Pvalue_Pion;
  Double32_t Angle_MotherPion;
  Int_t NHitsMDC_Pion;
  Int_t NHitsMinifiber_Pion;
  Int_t N_Pion;
  Int_t Pion_IsFromHyp; // 0-> No; 1-> Yes
  Int_t RealPDG_Pion;

  Double32_t Dist_Daughters;
  Double32_t ArmPod_Qt;
  Double32_t ArmPod_Alfa;

  DataML_Hyp(const THyphiAttributes& att, TTree* outT, TTree* outTs, TTree* outTb);
  void FillEvent(FullRecoEvent& RecoEvent) override;
};

//typedef TDataProcess<FullRecoEvent, MCAnaEventG4Sol> TDataProcessInterface;
template<class Out>
using TDataProcessInterface = TDataProcess<FullRecoEvent, Out>;


template <class Out>
class TFlatMCOutputML final : public TDataProcessInterface<Out>
{
public:
  const THyphiAttributes& att;

  TFlatMCOutputML(const THyphiAttributes& attr);
  ~TFlatMCOutputML();

  // int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator()(FullRecoEvent& RecoEvent, Out* OutTree) override;

private:
  int Exec(FullRecoEvent& RecoEvent, Out* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;

  int FlattenOut(FullRecoEvent& RecoEvent);

  TString namefileFlat;
  TFile* f_flat = nullptr;
  TTree* t_flat = nullptr;
  TTree* t_flatS = nullptr;
  TTree* t_flatB = nullptr;

  bool SplitTree = false;

  DataML* data_out = nullptr;

  struct LocalHists
  {
  };
  LocalHists LocalHisto;
};

#endif
