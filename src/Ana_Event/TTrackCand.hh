#ifndef TRACKCAND_h
#define TRACKCAND_h

#include "TObject.h"
#include "TString.h"
#include <vector>
#include <array>
#include "TLorentzVector.h"
#include "TVector3.h"

class TTrackCand : public TObject
{

  public:

  Int_t MCTrackID = -1;
  Int_t TrackID = -1;

  Int_t NSet = 0.;
  std::vector<UInt_t> SetLayerID = {};
  std::vector<UInt_t> SetHitID = {};

  Int_t FitStatus = -1.;
  Int_t Charge  = -999;
  Double32_t Chi2_C = -1.;
  Double32_t Chi2_L = -1.;

  std::array<Double32_t,5>  Seed_Par = {0., 0., 0., 0., 0.};
  std::array<Double32_t,25> Seed_Cov = {-9999., -9999., -9999., -9999., -9999.,
					-9999., -9999., -9999., -9999., -9999.,
					-9999., -9999., -9999., -9999., -9999.,
					-9999., -9999., -9999., -9999., -9999.,
					-9999., -9999., -9999., -9999., -9999.};

  Double_t Mass;
  Int_t pdgcode;
  Int_t MC_status;
  TLorentzVector MomMass;
  TVector3 Pos;
  TVector3 Mom;

  Int_t BarId;
  TVector3 Pos_PS;
  TVector3 Mom_PS;
  Float_t dE_PS;
  Float_t dx_PS;

  Float_t Beta;
  Double_t PathLength;
  Double_t TOF;

  Int_t NCent;
  Int_t Nmfib;

  TTrackCand() = default;
  TTrackCand(const TTrackCand& H);
  ~TTrackCand();

  virtual void Clear(Option_t* option = "");
  ClassDef(TTrackCand, 2)
};

#endif
