#include "TMcParticle.hh"
#include <iostream>
ClassImp(TMcParticle)

using namespace std;

TMcParticle::TMcParticle():type(""),Mc_id(-1),Mother_id(-1),Pdg(0),Charge(0),MomMass(0.,0.,0.,0.),Vtx(0.,0.,0.,0.),Weigth(1.),GeoAcc(true)
{
}

TMcParticle::TMcParticle(const TMcParticle& M)
{
  type = M.type;
  Mc_id = M.Mc_id;
  Mother_id = M.Mother_id;
  Pdg = M.Pdg;
  Charge = M.Charge;
  MomMass = M.MomMass;
  Vtx = M.Vtx;
  Weigth = M.Weigth;
  GeoAcc=M.GeoAcc;
}


TMcParticle::~TMcParticle()
{
}

void TMcParticle::Clear(Option_t *option)
{

  type = "";
  Charge = 0;
  Pdg = 0;
  Mc_id = -1;
  Mother_id = -1;
  MomMass.SetXYZT(0.,0.,0.,0.);
  Vtx.SetXYZT(0.,0.,0.,0.);
  Weigth = 1.;
  GeoAcc=true;
}

void TMcParticle::Print(Option_t* option) const
{

  std::cout<<"MC particle :"<<type<<" ["<<Pdg<<"] "<<" Id ("<<Mc_id<<";"<<Mother_id<<") "<<std::endl;
  std::cout<<" P = "<<MomMass.X()<<" "<<MomMass.Y()<<" "<<MomMass.Z()<<" "<<MomMass.M()<< " Mag: "<<MomMass.P()<<std::endl;
  std::cout<<" Vtx = "<<Vtx.X()<<" "<<Vtx.Y()<<" "<<Vtx.Z()<<" "<<Vtx.T()<<std::endl;
  

}
