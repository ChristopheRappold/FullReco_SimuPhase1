#ifndef HYPDCHRECOHIT_H
#define HYPDCHRECOHIT_H 1

// Hyp includes
#include "GFRecoHitIfc.h"
#include "GFPlanarHitPolicy.h"

// c++ headers
#include <ostream> 

class HypDchCylinderHit;


typedef GFRecoHitIfc<GFPlanarHitPolicy> PlanarRecoHit;

class HypDchRecoHit : public PlanarRecoHit {

public:

  /** Constructors  **/
  HypDchRecoHit();
  //  HypDchRecoHit(double r, double wireposx, double angle, double z,double sigr);
  //HypDchRecoHit(const HypDchCylinderHit* cylHit);
  HypDchRecoHit(const TVector3& end1,const TVector3& end2,double sin_a,double cos_a,double distance,double dist_sigr,bool m);
  /** Destructor  **/
  virtual ~HypDchRecoHit();
  
  /** Public method clone() **/
  virtual GFAbsRecoHit* clone();
  
  /** Public method residualVector(...) **/
  TMatrixT<double> residualVector(const GFAbsTrackRep* stateVector,const TMatrixT<double>& state,const GFDetPlane& d);

  /** Operations **/
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);

  virtual void  Print();
  //const std::string& getPolicyName() {return PolicyName;}
  int getNparHit() {return NparHitRep;}

private:

  static const int NparHitRep = 1;
  //static const std::string PolicyName;

  double _sangle;  ///<  sine of the angle of the wire
  double _cangle;  ///<  cosine of the angle of the wire
  double x_wirepos; ///<  x coordinate of wire in detector plane
  double y_wirepos;
  bool mirror;

  ClassDef(HypDchRecoHit,1)

};

#endif
