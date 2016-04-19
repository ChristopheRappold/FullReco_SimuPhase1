#include "THypernucleus.hh"

using namespace std;

ClassImp(THypernucleus)

THypernucleus::THypernucleus():MomMass(0.,0.,0.,0.),Vtx(0.,0.,0.,0.),MomMassD1(0.,0.,0.,0.),MomMassD2(0.,0.,0.,0.),MomMassD3(0.,0.,0.,0.)
{
}

THypernucleus::THypernucleus(const THypernucleus& H)
{
  type=H.type;
  pattern=H.pattern;
  Ndecay=H.Ndecay;
  //   Chi2=H.Chi2;
  Pvalue=H.Pvalue;
  InvMass=H.InvMass;
  Dist=H.Dist;
  MomMass=H.MomMass;
  Vtx=H.Vtx;
  MomMassD1=H.MomMassD1;
  MomMassD2=H.MomMassD2;
  MomMassD3=H.MomMassD3;

}

//THypernucleus& THypernucleus::operator=(const THypernucleus& H)
//{
//   THypernucleus temp(H);
//   return temp;
//}


THypernucleus::~THypernucleus()
{
}

void THypernucleus::Clear(Option_t *option)
{

  MomMass.SetXYZT(0.,0.,0.,0.);
  Vtx.SetXYZT(0.,0.,0.,0.);
  MomMassD1.SetXYZT(0.,0.,0.,0.);
  MomMassD2.SetXYZT(0.,0.,0.,0.);
  MomMassD3.SetXYZT(0.,0.,0.,0.);

}

