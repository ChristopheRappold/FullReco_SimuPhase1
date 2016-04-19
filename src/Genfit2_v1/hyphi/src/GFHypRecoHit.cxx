// This Class' Header ------------------
#include "GFHypRecoHit.h"

// C/C++ Headers ----------------------
#include <assert.h>

// Collaborating Class Headers --------
//#include "CbmMCPoint.h"
#include "LSLTrackRep.h"
#include "RKTrackRep.h"
#include "RKTrackRepXY.h"
//#include "GeaneTrackRep.h"
#include "GFDetPlane.h"

// Class Member definitions -----------

ClassImp(GFHypRecoHit)

//const std::string PolicyName = "PlanarHitPolicy";

GFHypRecoHit::~GFHypRecoHit()
{}

GFHypRecoHit::GFHypRecoHit()  : PlanarRecoHit(NparHitRep)
{}

GFHypRecoHit::GFHypRecoHit(double x, double y, double z, double sigx, double sigy,double angle) : PlanarRecoHit(NparHitRep)
{
  fHitCoord[0][0] = 0.;//x;
  fHitCoord[1][0] = 0.;//y;
  TVector3 temp_u(1.,0.,0.);
  temp_u.RotateY(angle);
  fPolicy.setDetPlane(GFDetPlane(TVector3(x,y,z),
		       temp_u,
		       TVector3(0,1,0)));
  fHitCov[0][0] = sigx*sigx;
  fHitCov[1][1] = sigy*sigy;
}

GFHypRecoHit::GFHypRecoHit(const TVector3& pos, double sigx, double sigy, double angle) : PlanarRecoHit(NparHitRep)
{
  fHitCoord[0][0] = 0.;//pos.X();
  fHitCoord[1][0] = 0.;//pos.Y();
  TVector3 temp_u(1.,0.,0.);
  temp_u.RotateY(angle);
  fPolicy.setDetPlane(GFDetPlane(TVector3(pos.X(),pos.Y(),pos.Z()),
		       temp_u,
		       TVector3(0,1,0)));
  fHitCov[0][0] = sigx*sigx;
  fHitCov[1][1] = sigy*sigy;
}
/*
GFHypRecoHit::GFHypRecoHit(CbmMCPoint* point)
  : PlanarRecoHit(NparHitRep)
{
  fHitCoord[0][0] = point->GetX();
  fHitCoord[1][0] = point->GetY();
  setDetPlane(DetPlane(TVector3(0,0,point->GetZ()),
		       TVector3(1,0,0),
		       TVector3(0,1,0)));
  double sigx=0.1;
  double sigy=0.1;
  fHitCov[0][0] = sigx*sigx;
  fHitCov[1][1] = sigy*sigy;
}
*/
GFHypRecoHit& GFHypRecoHit::operator= ( GFHypRecoHit& hit) 
{
  //assert(hit.getPolicyName() != fPolicyName);

  //_HMatrix = hit._HMatrix; 
  fHitCoord=hit.getRawHitCoord();
  fHitCov=hit.getRawHitCov();
  GFDetPlane pl(  hit.getDetPlane(NULL));
  fPolicy.setDetPlane(pl);

  return *this;
}

GFAbsRecoHit* GFHypRecoHit::clone()
{
  return new GFHypRecoHit(*this);
}


TMatrixT<double> GFHypRecoHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
  assert(stateVector!=NULL);
  if (dynamic_cast<const LSLTrackRep*>(stateVector) != NULL || dynamic_cast<const RKTrackRepXY*>(stateVector) != NULL) {
    // LSLTrackRep (x,y,x',y',q/p)
    TMatrixT<double> _HMatrix(NparHitRep,5);
    _HMatrix[0][0] = 1.;
    _HMatrix[0][1] = 0.;
    _HMatrix[0][2] = 0.;
    _HMatrix[0][3] = 0.;
    _HMatrix[0][4] = 0.;

    _HMatrix[1][0] = 0.;
    _HMatrix[1][1] = 1.;
    _HMatrix[1][2] = 0.;
    _HMatrix[1][3] = 0.;
    _HMatrix[1][4] = 0.;
    return _HMatrix;
    }
  else if (dynamic_cast<const RKTrackRep*>(stateVector) != NULL ) 
  {
    TMatrixT<double> HMatrix(2,5);
    
    HMatrix[0][0] = 0.;
    HMatrix[0][1] = 0.;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 1.;
    HMatrix[0][4] = 0.;
    
    HMatrix[1][0] = 0.;
    HMatrix[1][1] = 0.;
    HMatrix[1][2] = 0.;
    HMatrix[1][3] = 0.;
    HMatrix[1][4] = 1.;

    return HMatrix;
  }  
  /*else if (dynamic_cast<const GeaneTrackRep*>(stateVector) != NULL) {
    // Uses TrackParP (q/p,v',w',v,w)
    // coordinates are defined by detplane!
    _HMatrix.ResizeTo(NparHitRep,5);

    _HMatrix[0][0] = 0.;
    _HMatrix[0][1] = 0.;
    _HMatrix[0][2] = 0.;
    _HMatrix[0][3] = 1.;
    _HMatrix[0][4] = 0.;

    _HMatrix[1][0] = 0.;
    _HMatrix[1][1] = 0.;
    _HMatrix[1][2] = 0.;
    _HMatrix[1][3] = 0.;
    _HMatrix[1][4] = 1.;
    }*/
  else 
    {
      std::cerr << "GFHypRecoHit can only handle state"
		<< " vectors of type LSLTrackRep -> abort" // or GeaneTrackRep -> abort" 
		<< std::endl;
      throw;
    }
 
}

void GFHypRecoHit::Print()
{
  std::cout<<"GFHypRecoHit :"<<getPolicyName()<<std::endl;
  fHitCoord.Print();
  fHitCov.Print();
  fPolicy.Print();

}

/*bool operator== (const GFHypRecoHit& h1,const GFHypRecoHit& h2)
{
  if(h1.getDetPlane(NULL) == h2.getDetPlane(NULL) )
    {
      bool check_cov = true;
      for(int row=0;row<2;row++)
	for(int col=0;col<2;col++)
	  if(TMath::Abs(h1.getRawHitCov()[row][col]-h2.getRawHitCov()[row][col])>1e-9)
	    check_cov=false;

      return check_cov;
    }
  else
    return false;
}
*/

/*double 
GFHypRecoHit::residualScalar(AbsTrackRep* stateVector,
			    const TMatrixT<double>& state)
{
  throw;
}

*/
