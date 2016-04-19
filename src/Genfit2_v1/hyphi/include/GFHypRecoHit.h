#ifndef GFHYPRECOHIT_HH
#define GFHYPRECOHIT_HH

// Base Class Headers ----------------
#include "GFRecoHitIfc.h"
#include "GFPlanarHitPolicy.h"
#include "TMath.h"
// Collaborating Class Headers -------
//#include <ostream> // remove if you do not need streaming op

// Collaborating Class Declarations --
//class CbmMCPoint;

typedef GFRecoHitIfc<GFPlanarHitPolicy> PlanarRecoHit;

class GFHypRecoHit : public PlanarRecoHit {
public:

  // Constructors/Destructors ---------
  GFHypRecoHit();
  GFHypRecoHit(double x, double y, double z, double sigx, double sigy,double angle);
  GFHypRecoHit(const TVector3& pos, double sigx, double sigy,double angle);
  //GFHypRecoHit(CbmMCPoint* point);
  GFHypRecoHit& operator= ( GFHypRecoHit& hit);

  virtual ~GFHypRecoHit();

  virtual GFAbsRecoHit* clone() ;
  
  // Operations ----------------------
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);

  //friend bool operator== (const GFHypRecoHit& h1,const GFHypRecoHit& h2);

/*  virtual double residualScalar(AbsTrackRep* stateVector,
                                const TMatrixT<double>& state);
*/
  int getNparHit() {return NparHitRep;}
  virtual void Print() ;
  
private:

  // Private Data Members ------------
  static const int NparHitRep = 2;
  //static const std::string PolicyName;
  // Private Methods -----------------

public:
  ClassDef(GFHypRecoHit,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
