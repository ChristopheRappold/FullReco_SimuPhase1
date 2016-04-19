#include "TTrackSimple.hh"

using namespace std;

ClassImp(TTrackSimple)

TTrackSimple::TTrackSimple():Chi2(-1.),Mass(-1.),MomMass(0.,0.,0.,0.)
{
}

TTrackSimple::TTrackSimple(const TTrackSimple& H)
{
  Chi2=H.Chi2;
  Mass=H.Mass;
  MomMass=H.MomMass;
  
}

//TTrackSimple& TTrackSimple::operator=(const TTrackSimple& H)
//{
//   TTrackSimple temp(H);
//   return temp;
//}


TTrackSimple::~TTrackSimple()
{
}

void TTrackSimple::Clear(Option_t *option)
{

  MomMass.SetXYZT(0.,0.,0.,0.);

}

