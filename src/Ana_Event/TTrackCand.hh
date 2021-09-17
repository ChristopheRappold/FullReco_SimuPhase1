#ifndef TRACKCAND_h
#define TRACKCAND_h

#include "TObject.h"
#include "TString.h"
#include <vector>
#include <array>

class TTrackCand : public TObject
{

  public:

  Int_t MCTrackID = -1;
  Int_t TrackID = -1;

  Int_t NSet = 0.;
  std::vector<UInt_t> SetLayerID = {};
  std::vector<UInt_t> SetHitID = {};

  Int_t FitStatus = -1.;
  Int_t Charge  = 0.;
  Double32_t Chi2_C = -1.;
  Double32_t Chi2_L = -1.;

  std::array<Double32_t,5>  Seed_Par = {0., 0., 0., 0., 0.};
  std::array<Double32_t,25> Seed_Cov = {-9999., -9999., -9999., -9999., -9999.,
					-9999., -9999., -9999., -9999., -9999.,
					-9999., -9999., -9999., -9999., -9999.,
					-9999., -9999., -9999., -9999., -9999.,
					-9999., -9999., -9999., -9999., -9999.};

  TTrackCand() = default;
  TTrackCand(const TTrackCand& H) = default;
  ~TTrackCand();

  virtual void Clear(Option_t* option = "");
  ClassDef(TTrackCand, 1)
};

#endif
