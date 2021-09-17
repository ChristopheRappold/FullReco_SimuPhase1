#include "TTrackCand.hh"

ClassImp(TTrackCand)

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
  Charge = 0.;
  Chi2_C = -1.;
  Chi2_L = -1.;
  Seed_Par = {0.,0.,0.,0.,0.};
  Seed_Cov = {-9999., -9999., -9999., -9999., -9999.,
	      -9999., -9999., -9999., -9999., -9999.,
	      -9999., -9999., -9999., -9999., -9999.,
	      -9999., -9999., -9999., -9999., -9999.,
	      -9999., -9999., -9999., -9999., -9999.};
}
