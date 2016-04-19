#ifndef GFHypFIELD_H
#define GFHypFIELD_H

#include"GFAbsBField.h"
#include <map>

typedef std::map<int, std::map<int,std::map<int,double> > > map_magneticfield;

class GFHypFieldMap : public GFAbsBField{
 public:
  GFHypFieldMap(bool _fielmap, bool type_simu, bool _rot,const map_magneticfield& field);

  TVector3 get(const TVector3& pos) const;

 private:

  //int LoadFieldMap(std::string type);

  /** Limits of the field region **/
  Double_t fXmin;   
  Double_t fXmax;
  Double_t fYmin;
  Double_t fYmax;
  Double_t fZmin;
  Double_t fZmax;
  int FieldMap_X;
  int FieldMap_Y;
  int FieldMap_Z[2];

  /** Field components inside the field region **/
  Double_t fBx;
  Double_t fBy;
  Double_t fBz;

  bool Aladin_field;
  bool aladin_rotation;
  bool Aladin_fieldmap;


  double Aladin_width;// = 1.56*m;
  double Aladin_length;// = 1.76*m; // Field length is 1.4 m 1.7
  double Field_length;// = 1.48*m;
  double Aladin_gap;// = 0.5*m;
  double Aladin_angle;// = -7.3*degree;
  //  Aladin_angle = -0.0*degree;
  double DistanceToTarget;// = 1.5*m;
  //G4double Yoke_thickness = 0.5*m;
  //
  double DistanceFromtargetToAladinCenter;// =	DistanceToTarget + Aladin_length/2.0;

  //std::map<int, std::map<int,std::map<int,double> > > FieldMap;

  const map_magneticfield& FieldMap;

};

#endif
