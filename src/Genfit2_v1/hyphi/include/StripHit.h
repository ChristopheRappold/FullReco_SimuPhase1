#ifndef STRIPHIT_HH
#define STRIPHIT_HH

// Base Class Headers ----------------
#include <assert.h>

#include "GFRecoHitIfc.h"
#include "GFPlanarHitPolicy.h"

// Collaborating Class Headers -------
#include <ostream> // remove if you do not need streaming op

// Collaborating Class Declarations --

typedef GFRecoHitIfc<GFPlanarHitPolicy> PlanarRecoHit;

class StripHit : public PlanarRecoHit {
public:

  // Constructors/Destructors ---------
  StripHit();
  StripHit(TVector3 point,double res,int proj);
  StripHit(TVector3 point,double res,double Plane_angle);
  StripHit(TVector3 point,double res,double Plane_angle,double Angle_det);

  virtual ~StripHit();

  virtual GFAbsRecoHit* clone();
  
  // Operations ----------------------
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);
  void Print();


private:

  // Private Data Members ------------
  static const int NparHitRep = 1;

  // Private Methods -----------------

public:
  ClassDef(StripHit,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
