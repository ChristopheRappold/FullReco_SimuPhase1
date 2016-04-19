#include "TMath.h"
#include "GFHypFieldMap_new.h"
#include "Riostream.h"

namespace genfit {

double F_scale(double x,double p0,double p1,double p0b,double p1b)
{
  if(x<60. && x>200.)
    return 1.;

  double v_fit = TMath::Exp( p0+p1*x);
  double v_new = TMath::Exp( p0b+p1b*x);
  if(TMath::Abs(v_fit)<1e-9)
    return 1.;
  else
    return v_new/v_fit;
}

GFHypFieldMap_new::GFHypFieldMap_new(bool normalized_,double fac,bool Telsa,bool DoFieldScaling,double p1_scaling,bool moreScale,TVirtualMagField* field):normalized(normalized_),scale(DoFieldScaling),more_scaling(moreScale),factor(fac),FieldMap(field)
{
  if(Telsa)
    telsa = 0.1;
  else
    telsa = 1.;

  p0 = 1.67525;
  p1 = -3.82014e-2;

  p1b = p1_scaling;
  p0b=p0+(p1-p1b)*60.;

}

TVector3 GFHypFieldMap_new::get(const TVector3& Pos) const
{
  TVector3 temp_B;
  Double_t B[3];
  Double_t pos[3];
  
  pos[0]=Pos.X();
  pos[1]=Pos.Y();
  pos[2]=Pos.Z();
  
  FieldMap->Field(pos,B);
  
  temp_B.SetXYZ(B[0]*telsa,B[1]*telsa,B[2]*telsa);
  
  // if(scale)
  //   {
  //     double scaling_factor=F_scale(Pos.Z(),p0,p1,p0b,p1b);
  //     temp_B*=scaling_factor;
  //   }

  if(scale)
    {
      Double_t xl = pos[0]; 
      Double_t zl = pos[2];
      const Double_t Aladin_Angle = -5.6;
      const Double_t fAngleRot = Aladin_Angle*TMath::Pi()/180.;
      Double_t xl2 = xl*TMath::Cos(fAngleRot)-zl*TMath::Sin(fAngleRot);
      Double_t zl2 = zl*TMath::Cos(fAngleRot)+xl*TMath::Sin(fAngleRot);
	  
      double scaling_factor=F_scale(zl2,p0,p1,p0b,p1b);
      temp_B*=scaling_factor;
      
      if(more_scaling)
	{	  
	  if(zl2>184.)
	    {
	      double new_by = TMath::Exp( p0b+p1b*zl2);
	      temp_B.SetY(new_by);
	    }
	}
    }

  if(normalized)
    temp_B*=factor/0.750533;

  //std::cout<<" GFHypFieldMap_new "<<" normalized"<<Pos.X()<<" "<<Pos.Z()<<" B: "<<temp_B.X()<<" "<<temp_B.Y()<<" "<<temp_B.Z()<<std::endl;

  
  return temp_B;
}


void GFHypFieldMap_new::get(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz) const 
{
  Double_t B_field[3];
  Double_t pos[3]={posX,posY,posZ};
  
  FieldMap->Field(pos,B_field);
  
  Double_t temp_B[3]={B_field[0]*telsa,B_field[1]*telsa,B_field[2]*telsa};
  
  // if(scale)
  //   {
  //     double scaling_factor=F_scale(Pos.Z(),p0,p1,p0b,p1b);
  //     temp_B*=scaling_factor;
  //   }

  if(scale)
    {
      Double_t xl = pos[0]; 
      Double_t zl = pos[2];
      const Double_t Aladin_Angle = -5.6;
      const Double_t fAngleRot = Aladin_Angle*TMath::Pi()/180.;
      Double_t xl2 = xl*TMath::Cos(fAngleRot)-zl*TMath::Sin(fAngleRot);
      Double_t zl2 = zl*TMath::Cos(fAngleRot)+xl*TMath::Sin(fAngleRot);
	  
      double scaling_factor=F_scale(zl2,p0,p1,p0b,p1b);
      temp_B[0] *= scaling_factor;
      temp_B[1] *= scaling_factor;
      temp_B[2] *= scaling_factor;
      
      if(more_scaling)
	{	  
	  if(zl2>184.)
	    {
	      double new_by = TMath::Exp( p0b+p1b*zl2);
	      temp_B[1]=(new_by);
	    }
	}
    }

  if(normalized)
    {
      temp_B[0] *= factor/0.750533;
      temp_B[1] *= factor/0.750533;
      temp_B[2] *= factor/0.750533;
    }
  //std::cout<<" GFHypFieldMap_new "<<" normalized"<<Pos.X()<<" "<<Pos.Z()<<" B: "<<temp_B.X()<<" "<<temp_B.Y()<<" "<<temp_B.Z()<<std::endl;

  Bx = temp_B[0];
  By = temp_B[1];
  Bz = temp_B[2];
  //return temp_B;
}


// TVector3 GFHypFieldMap_new::get(const ROOT::Math::SVector<double,3>& Pos) const
// {
//   TVector3 temp_B;
//   Double_t B[3];
//   Double_t pos[3];
  
//   pos[0]=Pos[0];
//   pos[1]=Pos[1];
//   pos[2]=Pos[2];
  
//   FieldMap->Field(pos,B);
  
//   temp_B.SetXYZ(B[0]*telsa,B[1]*telsa,B[2]*telsa);
  
//   if(scale)
//     {
//       double scaling_factor=F_scale(Pos[2],p0,p1,p0b,p1b);
//       temp_B*=scaling_factor;
//     }

//   if(normalized)
//     temp_B*=factor/0.750533;
  
//   return temp_B;
// }

}
