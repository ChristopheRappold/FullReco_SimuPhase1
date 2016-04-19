#ifndef HYPDCHRECOHIT2_H
#define HYPDCHRECOHIT2_H 1

// Hyp includes
#include "GFRecoHitIfc.h"
#include "GFWireHitPolicy.h"

typedef GFRecoHitIfc<GFWireHitPolicy> WireRecoHit;

class HypDchRecoHit2 : public WireRecoHit {

public:

  HypDchRecoHit2();
  HypDchRecoHit2(const TVector3& end1,const TVector3& end2,double distance,double dist_sigr,double max);
  /** Default Constructor for dch **/
  //HypDchRecoHit2(const HypDchCylinderHit* cylHit);

  /** virtual destructor **/
  virtual ~HypDchRecoHit2();

  /** Public method clone(...) creates a clone of argument**/
   virtual GFAbsRecoHit* clone(){return new HypDchRecoHit2(*this);};

  /** Public method Print() **/
  virtual void  Print();

 // Operations ----------------------
  TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);
  int getNparHit() {return NparHitRep;}

 protected:
  static const int NparHitRep = 7;
  //static const std::string PolicyName;
 public:
  ClassDef(HypDchRecoHit2,2)

};

#endif
