#include "THyphiTrack_v2.hh"
#include "Riostream.h"

using namespace std;

ClassImp(THyphiTrack)
  
THyphiTrack::THyphiTrack():type(""),MC_status(0),Chi2(-1.),Chi2_X(-1.),Chi2_Y(-1.),Mass(-1.),pdgcode(0),MomMass(0.,0.,0.,0.),Mom(0.,0.,0.),Charge(0),BarId(0),Beta(0),HitTR1(-999,-999,-999),HitTR2(-999,-999,-999),HitDC1(-999,-999,-999),HitDC1p(-999,-999,-999),NumDC1(-1,-1,-1),HitDC2(-999,-999,-999),NumDC2(-1,-1,-1),Pval2(-1),TofsBar(-1),PathLength(-999),TOF(-999),MomIni(-999,-999,-999),RChiIni(-1),BetaIni(-999)//,state(5,1),cov(5,5)
{
  
  //   std::cout<<"THyphiTrack Constructor"<<std::endl;
  //   std::cout<<"come 11 "<<state.GetNrows()<<" "<<state.GetNcols()<<std::endl;
  //   std::cout<<"come 21 "<<cov.GetNrows()<<" "<<cov.GetNcols()<<std::endl;
  for(int i=0;i<5;i++)
    {
      State[i] = -9999;
      for(int j=0;j<5;j++)
	{
	  Cov[i][j]=-9999;
	}
    }
  
  MomTof.SetXYZ(-999,-999,-999);
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
  
  Charge=H.Charge;
  BarId =H.BarId;
  dE = H.dE;
  Beta = H.Beta;
  RefPoint = H.RefPoint;
  HitTR1 = H.HitTR1;
  HitTR2 = H.HitTR2;
  HitDC1 = H.HitDC1;
  NumDC1 = H.NumDC1;
  HitDC1p = H.HitDC1p;
  HitDC2 = H.HitDC2;
  NumDC2 = H.NumDC2;
  HitTOF = H.HitTOF;
  Pval2  = H.Pval2;
  TofsBar= H.TofsBar;

  PathLength = H.PathLength;
  TOF = H.TOF;

  MomIni = H.MomIni;
  RChiIni = H.RChiIni;
  //  PathLengthIni = H.PathLengthIni;
  BetaIni = H.BetaIni;
  //  state.ResizeTo(5,1);
  //  cov.ResizeTo(5,5);
  //   std::cout<<"come 1 "<<H.state.GetNrows()<<" * "<<H.state.GetNcols()<<" == "<<state.GetNrows()<<" "<<state.GetNcols()<<std::endl;
  //   std::cout<<"come 2 "<<H.cov.GetNrows()<<" * "<<H.cov.GetNcols()<<" ==" <<cov.GetNrows()<<" "<<cov.GetNcols()<<std::endl;
  
  //  state = H.state;
  //  cov = H.cov;
  //   for(int i=0;i<5;i++)
  //     {
  //       std::cout<<"come 2 "<<std::endl;
  //       state[i][0] = H.state[i][0];
  //       for(int j=0;j<5;j++)
  // 	cov[i][j] = H.cov[i][j];
  //     }
  
  /////////////////////////
  for(int i=0;i<5;i++)
    {
      State[i] = H.State[i];
      for(int j=0;j<5;j++)
	{
	  Cov[i][j]=H.Cov[i][j];
	}
    }
  ///////////////////////
  MomTof = H.MomTof;
  ////////////////////////

}

//THyphiTrack& THyphiTrack::operator=(const THyphiTrack& H)
//{
//   THyphiTrack temp(H);
//   return temp;
//}
void THyphiTrack::SetPara(const TString& name,Int_t MC,Double_t chi2,Double_t chi2_x, Double_t chi2_y,Double_t mass, Int_t pdg,const TVector3& vec)
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
void THyphiTrack::SetPara(const TString& name,Int_t MC,Double_t chi2,Double_t chi2_x, Double_t chi2_y,Double_t mass, Int_t pdg,const TVector3& vec,Int_t charge,Int_t barid,Double_t beta,const TVector3& hit_tr1,const TVector3& hit_tr2, const TVector3 hit_dc2,const TVector3 num_dc2,Double_t pval2,Int_t tofsbar )
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

  Charge=charge;
  BarId=barid;
  Beta=beta;
  HitTR1=hit_tr1;
  HitTR2=hit_tr2;
  HitDC2=hit_dc2;
  NumDC2 = num_dc2;
  Pval2 =pval2;
  TofsBar=tofsbar;

}

THyphiTrack::~THyphiTrack()
{
}

void THyphiTrack::Clear(Option_t *option)
{

  MomMass.SetXYZT(0.,0.,0.,0.);
  Mom.SetXYZ(0.,0.,0.);

  RefPoint.SetXYZ(-999.,-999.,-999.);

  HitTR1.SetXYZ(-999.,-999.,-999.);
  HitTR2.SetXYZ(-999.,-999.,-999.);
  HitDC1.SetXYZ(-999.,-999.,-999.);
  HitDC1p.SetXYZ(-999.,-999.,-999.);
  HitDC2.SetXYZ(-999.,-999.,-999.);

  NumDC1.SetXYZ(-1,-1,-1);
  NumDC2.SetXYZ(-1,-1,-1);

  HitTOF.SetXYZ(-999.,-999.,-999.);

  type="";
  MC_status=-999;
  Chi2 = -999.;
  Chi2_X = -999.;
  Chi2_Y = -999.;
  Mass = -999.;
  pdgcode = -999;

  Charge=-999;
  BarId=-999;
  dE = -999;
  Beta=-999.;
  Pval2 = -999;
  TofsBar=-1;

  PathLength = -999;
  TOF = -999;

  MomIni.SetXYZ(-999,-999,-999);
  RChiIni = -1;
  //  PathLengthIni =-999;
  BetaIni =-999;
  //  std::cout<<"come 3 "<<std::endl;
  //  state.Zero();
  //  cov.Zero();

//   for(int i=0;i<5;i++)
//     {
//       std::cout<<"come 4 "<<std::endl;
//       state[i][0] = 0;
//       for(int j=0;j<5;j++)
// 	cov[i][j] = 0;
//     }

  //////////////////////////////
  for(int i=0;i<5;i++)
    {
      State[i] = -9999;
      for(int j=0;j<5;j++)
	{
	  Cov[i][j]=-9999;
	}
    }
  //////////////////////////////
  MomTof.SetXYZ(-999,-999,-999);

}

