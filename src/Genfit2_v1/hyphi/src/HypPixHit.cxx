// This Class' Header ------------------
#include "HypPixHit.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "RKTrackRep.h"
#include "GFDetPlane.h"
//#include "GFRectFinitePlane.h"
#include "TRandom.h"
#include "TMath.h"

// Class Member definitions -----------

ClassImp(HypPixHit)

HypPixHit::~HypPixHit()
{}

HypPixHit::HypPixHit()
  : PlanarRecoHit(NparHitRep)
{}

HypPixHit::HypPixHit(const GFDetPlane& pl,double x_res,double y_res) : PlanarRecoHit(NparHitRep){

  fPolicy.setDetPlane(pl);
  fHitCoord[0][0] = 0.;//gRandom->Gaus(0.,res);
  fHitCoord[1][0] = 0.;//gRandom->Gaus(0.,y_res);//gRandom->Gaus(0.,res);
  fHitCov[0][0] = x_res*x_res;
  fHitCov[1][1] = y_res*y_res;

}

HypPixHit::HypPixHit(TVector3 point,double x_res,double y_res)
  : PlanarRecoHit(NparHitRep){

  fHitCov[0][0] = x_res*x_res;
  fHitCov[1][1] = y_res*y_res;
  GFDetPlane d;

  d.setO(point);
  d.setNormal(TVector3(point.X(),point.Y(),0.));
  //d.setNormal(TVector3(0.,0.,1.));
  
  
    fHitCoord[0][0] = 0.;//gRandom->Gaus(0,res);
    fHitCoord[1][0] = 0.;//gRandom->Gaus(0,y_res);

  //  GFAbsFinitePlane * fin = new GFRectFinitePlane(-20,20,-3,3);
  //  d.setFinitePlane(fin);

  fPolicy.setDetPlane(d);
}

GFAbsRecoHit* 
HypPixHit::clone(){
  return new HypPixHit(*this);
}


TMatrixT<double> HypPixHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
  if ( dynamic_cast<const RKTrackRep*>(stateVector) != NULL ) {
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
  else {
    std::cerr << "HypPixHit can only handle state"
	      << " vectors of type RKTrackRep -> abort" 
	      << std::endl;
    throw;
  }
}
