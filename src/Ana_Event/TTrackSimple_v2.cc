#include "TTrackSimple_v2.hh"

using namespace std;

ClassImp(TTrackSimple)

TTrackSimple::TTrackSimple():type(""),MC_status(0),Chi2(-1.),Chi2_X(-1.),Chi2_Y(-1.),Mass(-1.),pdgcode(0),MomMass(0.,0.,0.,0.),Mom(0.,0.,0.)
{
}

TTrackSimple::TTrackSimple(const TTrackSimple& H)
{
  type=H.type;
  MC_status=H.MC_status;
  Chi2=H.Chi2;
  Chi2_X=H.Chi2_X;
  Chi2_Y=H.Chi2_Y;
  Mass=H.Mass;
  pdgcode=H.pdgcode;
  MomMass=H.MomMass;
  Mom=H.Mom;
}

//TTrackSimple& TTrackSimple::operator=(const TTrackSimple& H)
//{
//   TTrackSimple temp(H);
//   return temp;
//}
void TTrackSimple::SetPara(const TString& name,Int_t MC,Double_t chi2,Double_t chi2_x, Double_t chi2_y,Double_t mass, Int_t pdg,const TVector3& vec)
{
  type=name;
  MC_status=MC;
  Chi2 = chi2;
  Chi2_X = chi2_x;
  Chi2_Y = chi2_y;
  Mass = mass;
  pdgcode = pdg;
  MomMass.SetXYZM(vec.X(),vec.Y(),vec.Z(),mass);
  Mom = vec;

}

TTrackSimple::~TTrackSimple()
{
}

void TTrackSimple::Clear(Option_t *option)
{

  MomMass.SetXYZT(0.,0.,0.,0.);
  Mom.SetXYZ(0.,0.,0.);
}

