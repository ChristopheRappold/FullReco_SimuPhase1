#include "TMath.h"
#include "GFHypFieldMap.h"
#include "Riostream.h"

typedef std::map<int, std::map<int,std::map<int,double> > > map_magneticfield;


GFHypFieldMap::GFHypFieldMap(bool _fieldmap, bool _simu, bool _rot,const map_magneticfield& field):FieldMap(field)
{

  if(_fieldmap==true && _simu == false)
    {
      FieldMap_X = 80;
      FieldMap_Y = 27;
      FieldMap_Z[0] = -185; // cm
      FieldMap_Z[1] = 185; // cm

      double m=100.;

      Aladin_field = true;
      aladin_rotation = _rot;
      Aladin_fieldmap = true;
    
      Aladin_width = 1.56*m;
      Aladin_length = 1.76*m; // Field length is 1.4 m 1.7
      Field_length = 1.48*m;
      Aladin_gap = 0.5*m;
      Aladin_angle = -7.3;
      if(_simu==false)
        Aladin_angle= -5.6;
      
      if(aladin_rotation)
        Aladin_angle = 0.0;
      
      DistanceToTarget = 1.5*m;
  
      DistanceFromtargetToAladinCenter = DistanceToTarget + Aladin_length/2.0;
      fBx=0.;
      fBy=0.75;//  -0.7 Telsa
      //fBy=-0.7;//  -0.7 Telsa
      fBz=0.;
      Field_length = 1.5*2*m; // Field length is 1.4 m 1.7

      /*if(aladin_rot)
	{
	translation_T2M trans(DistanceFromtargetToAladinCenter);
	std::vector<double> steping_for_rot;
	steping_for_rot.resize(steping.size());
	std::transform(steping.begin(),steping.end(),steping_for_rot.begin(),trans);
	std::swap(stepZinField,steping_for_rot);
	fZmin = stepZinField.front();
	fZmax = stepZinField.back();
	}
	else
	{*/
      fZmin = 0.89*(DistanceFromtargetToAladinCenter-Field_length/2.);
      fZmax = 1.11*(DistanceFromtargetToAladinCenter+Field_length/2.);
      //}

    }
  else
    {
      Aladin_field = _simu;
      aladin_rotation = _rot;
      Aladin_fieldmap=false;
      double Aladin_angle_old = 0.0;
      if(Aladin_field)
	{
	  double m=100.;
	  //std::cout<<"Setup Aladin field !"<<std::endl;
	  Aladin_width = 1.56*m;
	  Aladin_length = 1.76*m; // Field length is 1.4 m 1.7
	  Field_length = 1.48*m;
	  Aladin_gap = 0.5*m;
	  //Aladin_angle = -7.3;
	  Aladin_angle = -7.3;
	  if(aladin_rotation)
	    {
	      Aladin_angle = 0.000;
	      Aladin_angle_old = -7.3;
	    }
	  DistanceToTarget = 1.5*m;
	  //G4double Yoke_thickness = 0.5*m;
	  
	  DistanceFromtargetToAladinCenter = DistanceToTarget + Aladin_length/2.0;
	  fBx=0.;
	  fBy=0.7;//  .7 Telsa
	  fBz=0.;
	  
	  //fZmin = DistanceToTarget-m/10.;
	  double temp_Z = 1./TMath::Cos(Aladin_angle*TMath::Pi()/180.)*Aladin_length/2.+(Aladin_width/2.-Aladin_length/2.*TMath::Tan(Aladin_angle*TMath::Pi()/180.))*TMath::Cos(Aladin_angle*TMath::Pi()/180.);
	  fZmin = DistanceFromtargetToAladinCenter-temp_Z;
	  fZmax = DistanceFromtargetToAladinCenter+temp_Z;
	}
      else
        {
	  fXmin=0.;
	  fXmax=0.;
	  fYmin=0.;
	  fYmax=0.;
	  fZmin=0.;
	  fZmax=0.;
	  
	  fBx=0.;
	  fBy=0.;
	  fBz=0.;
	  
	  if(Aladin_field)
	    {
	      std::cout<<"!> Incompatibility Field is OFF but asking for ALADIN FIELD ON !"<<std::endl;
	    }
	}

    }

}


TVector3 GFHypFieldMap::get(const TVector3& Pos) const
{
  TVector3 temp_B(0.,0.,0.);
  /*	std::cout<<"get field";
	Pos.Print();
	std::cout<<"field map:("<<fXmin<<"/"<<fXmax<<");("<<fYmin<<"/"<<fYmax<<");("<<fZmin<<"/"<<fZmax<<");"<<std::endl;*/
  if(Aladin_field==false)
    {
      if(Pos.X()>fXmin && Pos.X()<fXmax)
	if(Pos.Y()>fYmin && Pos.Y()<fYmax)
	  if(Pos.Z()>fZmin && Pos.Z()<fZmax)
	    {
	      temp_B.SetXYZ(fBx,fBy,fBz);
	      //std::cout<<"field setup";
	    }
      //temp_B.Print();
    }
  else
    {
      //G4double DistanceFromtargetToAladinCenter = 2.0*m;
      //
      double z_shifted = Pos.Z() - DistanceFromtargetToAladinCenter;
      double z = z_shifted*cos(-Aladin_angle*TMath::DegToRad()) +	Pos.X()*sin(-Aladin_angle*TMath::DegToRad());
      double x = -1.0*z_shifted*sin(-Aladin_angle*TMath::DegToRad()) + Pos.X()*cos(-Aladin_angle*TMath::DegToRad());
      double y = Pos.Y();

      if(aladin_rotation)
	{
	  z = Pos.Z();
	  x = Pos.X();
	  y = Pos.Y();
	}

      //if(TMath::Abs(z)<=Field_length/2.0 && TMath::Abs(x)<Aladin_width/2.0 && TMath::Abs(y)<Aladin_gap/2.)
      if(FieldMap_Z[0]<z && z<FieldMap_Z[1] && TMath::Abs(x)<FieldMap_X && TMath::Abs(y)<FieldMap_Y)
	{
	  if(Aladin_fieldmap==false)
	    {
	      //G4cout << Point[2] << " " << z_shifted << G4endl;
	      //G4cout << x << " " << y << " " << z << G4endl;
	      temp_B.SetXYZ(fBx,fBy,fBz);
	    }   
	  else
	    {
	      int zz = static_cast<int>(z);
	      int yy = static_cast<int>(y);
	      int xx = static_cast<int>(x);
	      std::map<int, std::map<int,std::map<int,double> > >::const_iterator it_B_x = FieldMap.find(xx);
	      double B_map = 0.;
	      if(it_B_x == FieldMap.end())
		std::cout<<"E> Genfit FieldMap requested a x value out of the field ! ["<<xx<<" "<<yy<<" "<<zz<<"]"<<std::endl;
	      else
		{
		  std::map<int,std::map<int,double> >::const_iterator it_B_y = it_B_x->second.find(yy);
		  if(it_B_y == it_B_x->second.end())
		    std::cout<<"E> Genfit FieldMap requested a y value out of the field ! ["<<xx<<" "<<yy<<" "<<zz<<"]"<<std::endl;
		  else
		    {
		      std::map<int,double>::const_iterator it_B_z = it_B_y->second.find(zz);
		      if(it_B_z == it_B_y->second.end())
			std::cout<<"E> Genfit FieldMap requested a z value out of the field ! ["<<xx<<" "<<yy<<" "<<zz<<"]"<<std::endl;
		      else
			{
			  B_map = it_B_z->second;
			}
		    } 
		} 
	      temp_B.SetXYZ(fBx,B_map,fBz);
	      //double B_map = FieldMap[xx][yy][zz];
	    }
	}
      else 
	{
	  temp_B.SetXYZ(0.,0.,0.);
	}
    }
  
  return temp_B;
}
