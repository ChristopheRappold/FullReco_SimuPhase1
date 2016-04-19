
#ifndef HYPPIXHIT_HH
#define HYPPIXHIT_HH

// Base Class Headers ----------------
#include "GFRecoHitIfc.h"
#include "GFPlanarHitPolicy.h"
#include "GFDetPlane.h"

// Collaborating Class Headers -------
#include <ostream> // remove if you do not need streaming op


// Collaborating Class Declarations --

typedef GFRecoHitIfc<GFPlanarHitPolicy> PlanarRecoHit;

class HypPixHit : public PlanarRecoHit {
public:

  // Constructors/Destructors ---------
  HypPixHit();
  HypPixHit(TVector3 point,double xres,double yres);
  HypPixHit(const GFDetPlane& pl,double xres,double yres);

  virtual ~HypPixHit();

  virtual GFAbsRecoHit* clone();
  
  // Operations ----------------------
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);


private:
  static const int NparHitRep = 2;

  // Private Methods -----------------

public:
  ClassDef(HypPixHit,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
