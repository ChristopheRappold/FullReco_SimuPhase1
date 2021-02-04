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

typedef TDataProcess<FullRecoEvent, MCAnaEventG4Sol> TDataProcessInterface;

class TFlatMCOutputML final : public TDataProcessInterface
{
public:
  const THyphiAttributes& att;

  TFlatMCOutputML(const THyphiAttributes& attr);
  ~TFlatMCOutputML();

  // int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator()(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree) override;

private:
  int Exec(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;

  int FlattenOut(FullRecoEvent& RecoEvent);

  TString namefileFlat;
  TFile* f_flat = nullptr;
  TTree* t_flat = nullptr;

  Float_t b_tx, b_ty, b_vx, b_vy, b_vz, b_x, b_y, b_z, b_pt, b_phi, b_theta;
  Float_t a_tx, a_ty, a_vx, a_vy, a_vz, a_x, a_y, a_z, a_pt, a_phi, a_theta;
  Float_t poq, qop;
  Float_t q;
  Float_t tof;
  Float_t psb_z;
  Int_t psce_psbe;
  Int_t pdg;

  struct LocalHists
  {
  };
  LocalHists LocalHisto;
};

#endif
