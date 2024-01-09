#include "TTrackCand.hh"

ClassImp(TTrackCand)

TTrackCand::TTrackCand(const TTrackCand& H)
{
  MCTrackID = H.MCTrackID;
  TrackID = H.TrackID;
  NSet = H.NSet;
  SetLayerID = H.SetLayerID;
  SetHitID = H.SetHitID;
  FitStatus = H.FitStatus;
  Charge = H.Charge;
  Chi2_C = H.Chi2_C;
  Chi2_L = H.Chi2_L;
  for(int i=0; i<5; ++i)
    Seed_Par[i] = H.Seed_Par[i];
  for(int i=0; i<25; ++i)
    Seed_Cov[i] = H.Seed_Cov[i];
  Mass = H.Mass;
  pdgcode = H.pdgcode;
  MC_status = H.MC_status;
  MomMass = H.MomMass;
  Pos = H.Pos;
  Mom = H.Mom;
  BarId = H.BarId;
  Pos_PS = H.Pos_PS;
  Mom_PS = H.Mom_PS;
  dE_PS = H.dE_PS;
  dx_PS = H.dx_PS;
  Beta = H.Beta;
  PathLength = H.PathLength;
  TOF = H.TOF;
  NCent = H.NCent;
  Nmfib = H.Nmfib;
}

TTrackCand::~TTrackCand()
{ }

void TTrackCand::Clear(Option_t *option)
{
  MCTrackID = -1;
  TrackID = -1;
  NSet = 0;
  SetLayerID.clear();
  SetHitID.clear();
  FitStatus = -1;
  Charge = -999;
  Chi2_C = -1.;
  Chi2_L = -1.;
  Seed_Par = {0.,0.,0.,0.,0.};
  Seed_Cov = {-9999., -9999., -9999., -9999., -9999.,
	      -9999., -9999., -9999., -9999., -9999.,
	      -9999., -9999., -9999., -9999., -9999.,
	      -9999., -9999., -9999., -9999., -9999.,
	      -9999., -9999., -9999., -9999., -9999.};
  Mass = -1.;
  pdgcode = -1;
  MC_status = -1;
  MomMass.SetXYZM(0.,0.,0.,0.);
  Pos.SetXYZ(0.,0.,0.);
  Mom.SetXYZ(0.,0.,0.);
  BarId = -1;
  Pos_PS.SetXYZ(0.,0.,0.);
  Mom_PS.SetXYZ(0.,0.,0.);
  dE_PS = -1.;
  dx_PS = -1.;
  Beta = -1.;
  PathLength = -1.;
  TOF = -1.;
  NCent = -1;
  Nmfib = -1;
}
