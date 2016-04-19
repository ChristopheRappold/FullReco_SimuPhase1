#ifndef GFHypFIELDNEW_H
#define GFHypFIELDNEW_H

#include "AbsBField.h"
#include "TVirtualMagField.h"

namespace genfit {


class GFHypFieldMap_new : public AbsBField 
{
 public:
  GFHypFieldMap_new(bool s,double f, bool Telsa,bool DoFS,double p1_new,bool ms,TVirtualMagField* field);

  virtual TVector3 get(const TVector3& pos) const;
  //TVector3 get(const ROOT::Math::SVector<double,3>& pos) const;
  virtual void get(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz) const ;

private:

  bool normalized;
  bool scale;
  bool more_scaling;

  double factor;
  double telsa;
  // exp(p0b + p1b*x) with p0b = p0 + (p1-p1b)*x0 and x0 = 60.
  double p0b;
  double p1b;
  // exp(p0+p1*x) : tail scaling 
  double p0;
  double p1;

  //TF1* f_fit;
  //TF1* f_new;

  TVirtualMagField* const FieldMap;

};

}
#endif
