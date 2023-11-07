#include "THyphiTrack_v4.hh"
#include "Riostream.h"

using namespace std;

ClassImp(THyphiTrack)
  
THyphiTrack::THyphiTrack():type(""),MC_status(0),Chi2(-1.),Chi2_X(-1.),Chi2_Y(-1.),Mass(-1.),pdgcode(0),MomMass(0.,0.,0.,0.),Mom(0.,0.,0.),Pos_PS(0.,0.,0.),Mom_PS(0.,0.,0.),Charge(0),BarId(0),dE(0.),dx(0.),Beta(0),RefPoint(-999.,-999.,-999.),Pval2(-1),TofsBar(-1),PathLength(-999),TOF(-999),MomIni(-999,-999,-999),RChiIni(-1),PathLengthIni(-999.),TOFIni(-999.),BetaIni(-999),MassIni(-999.),Sim2Vtx(-999.,-999.,-999.,-999.),NCent(0),iterNum(-1),Nmfib(0)//,state(5,1),cov(5,5)
{
  
  //   std::cout<<"THyphiTrack Constructor"<<std::endl;
  //   std::cout<<"come 11 "<<state.GetNrows()<<" "<<state.GetNcols()<<std::endl;
  //   std::cout<<"come 21 "<<cov.GetNrows()<<" "<<cov.GetNcols()<<std::endl;
  for(int i=0;i<6;i++)
    {
      State[i] = -9999;
      for(int j=0;j<6;j++)
	{
	  Cov[i][j]=-9999;
	}
    }
  for(int i=0; i<17; ++i){
    for(int j=0; j<3; ++j){
      ResMDC[i][j] = -999.;
      WeightMDC[i][j] = -999.;
    }
  }
  for(int i=0; i<9; ++i){
      ResFiber[i] = -999.;
      WeightFiber[i] = -999.;
  }
  for(int i=0; i<6; ++i){
      ResMiniFiber[i] = -999.;
      WeightMiniFiber[i] = -999.;
  }
  for(int i=0; i<2; ++i){
    ResPSCE[i] = -999.;
    WeightPSCE[i] = -999.;
  }

}

THyphiTrack::THyphiTrack(const THyphiTrack& H)
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

  Pos_PS=H.Pos_PS;
  Mom_PS=H.Mom_PS;

  Charge=H.Charge;
  BarId =H.BarId;
  dE = H.dE;
  dx = H.dx;
  Beta = H.Beta;
  RefPoint = H.RefPoint;
  Pval2  = H.Pval2;
  TofsBar= H.TofsBar;

  PathLength = H.PathLength;
  TOF = H.TOF;

  MomIni = H.MomIni;
  RChiIni = H.RChiIni;
  PathLengthIni = H.PathLengthIni;
  TOFIni = H.TOFIni;
  BetaIni = H.BetaIni;
  MassIni = H.MassIni;

  for(int i=0;i<6;i++)
    {
      State[i] = H.State[i];
      for(int j=0;j<6;j++)
        Cov[i][j]=H.Cov[i][j];
    }

  NCent  = H.NCent;
  iterNum= H.iterNum;
  Nmfib = H.Nmfib;

  for(int i=0; i<17; ++i){
    for(int j=0; j<3; ++j){
      ResMDC[i][j] = H.ResMDC[i][j];
      WeightMDC[i][j] = H.WeightMDC[i][j];
    }
  }
  for(int i=0; i<9; ++i){
    ResFiber[i] = H.ResFiber[i];
    WeightFiber[i] = H.WeightFiber[i];
  }
  for(int i=0; i<6; ++i){
    ResMiniFiber[i] = H.ResMiniFiber[i];
    WeightMiniFiber[i] = H.WeightMiniFiber[i];
  }
  for(int i=0; i<2; ++i){
    ResPSCE[i] = H.ResPSCE[i];
    WeightPSCE[i] = H.WeightPSCE[i];
  }
}

//THyphiTrack& THyphiTrack::operator=(const THyphiTrack& H)
//{
//   THyphiTrack temp(H);
//   return temp;
//}

THyphiTrack::~THyphiTrack()
{ }

void THyphiTrack::Clear(Option_t *option)
{

  MomMass.SetXYZT(0.,0.,0.,0.);
  Mom.SetXYZ(0.,0.,0.);

  Pos_PS.SetXYZ(0.,0.,0.);
  Mom_PS.SetXYZ(0.,0.,0.);

  type="";
  MC_status=-999;
  Chi2 = -999.;
  Chi2_X = -999.;
  Chi2_Y = -999.;
  Mass = -999.;
  pdgcode = 0;

  Charge=-999;
  BarId=-999;
  dE = 0.;
  dx = 0.;
  Beta=-999.;
  RefPoint.SetXYZ(-999.,-999.,-999.);
  Pval2 = -999;
  TofsBar=-1;

  PathLength = -999;
  TOF = -999;

  MomIni.SetXYZ(-999,-999,-999);
  RChiIni = -1;
  PathLengthIni = -999.;
  TOFIni = -999.;
  BetaIni = -999.;
  MassIni = -999.;

  for(int i=0;i<6;i++)
    {
      State[i] = -9999;
      for(int j=0;j<6;j++)
	{
	  Cov[i][j]=-9999;
	}
    }

  NCent = 0;
  iterNum = -1;
  Nmfib = 0;

  for(int i=0; i<17; ++i){
    for(int j=0; j<3; ++j){
      ResMDC[i][j] = -999.;
      WeightMDC[i][j] = -999.;
    }
  }
  for(int i=0; i<9; ++i){
    ResFiber[i] = -999.;
    WeightFiber[i] = -999.;
  }
  for(int i=0; i<6; ++i){
    ResMiniFiber[i] = -999.;
    WeightMiniFiber[i] = -999.;
  }
  for(int i=0; i<2; ++i){
    ResPSCE[i] = -999.;
    WeightPSCE[i] = -999.;
  }

}

