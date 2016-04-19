
// The DCH point is not defined by the (x, y, z) coordinates, but we ONLY know the
// (x, y) coordinates of the wire center (not the z!); moreover we know the drift
// radius. In order to find the MEASURED POINT to be used in the Kalman fit a 
// preliminary propagation must be performed to the point of closest approach to
// the firing wire.
// The plane is a virtual detector plane defined via:
// dj = wvi - vpf     (U coordinate in DetPlane)
// dk = wiredirection (V coordinate in DetPlane)
//
// input: 7 entries
// 0-1-2 ==> x,y,z of the 1st extremity of the firing wire;
// 3-4-5 ==> x,y,z of the 2nd extremity of the firing wire;
// 6     ==> drift radius;

// Hyp Headers 
#include "HypDchRecoHit2.h"
//#include "HypDchCylinderHit.h"
#include "GFDetPlane.h"

#include "RKTrackRep.h"
#include "LSLTrackRep.h"
#include "SlTrackRep.h"
//ROOT includes
#include "TMath.h"

// C/C++ Headers 
#include <assert.h>
#include <exception>
#include <iostream> 
#include <cmath>

using std::cout;
using std::endl;


ClassImp(HypDchRecoHit2)

//const std::string PolicyName ="WireHitPolicy";

HypDchRecoHit2::~HypDchRecoHit2()
{;}

HypDchRecoHit2::HypDchRecoHit2(): WireRecoHit(NparHitRep)
{;}

//------------------------------------------------------------------------------
HypDchRecoHit2::HypDchRecoHit2(const TVector3& end1,const TVector3& end2,double distance,double dist_sigr,double max_dist): WireRecoHit(NparHitRep)
{
  //_sangle = TMath::Sin(angle);
  //_cangle = TMath::Cos(angle);
 
  if(end1.Z()!=end2.Z())
    {
      std::cout<<"E> HypDchRecoHit2 not same z";
      end1.Print();
      end2.Print();
    }

 fHitCoord[0][0] = end1.X();
 fHitCoord[1][0] = end1.Y();
 fHitCoord[2][0] = end1.Z();
 fHitCoord[3][0] = end2.X();
 fHitCoord[4][0] = end2.Y();
 fHitCoord[5][0] = end2.Z();
 fHitCoord[6][0] = distance;

  for(int i = 0; i < NparHitRep; i++) 
    for(int j = 0; j < NparHitRep; j++) 
      fHitCov[i][j] = 0.;

  fHitCov[6][6] = dist_sigr;//*dist_sigr;

  fPolicy.setMaxDistance(max_dist);

  //std::cout<<" HypDchRecoHit2 : Detplane "<<std::endl<<"wire end :";
  
  //end1.Print();
  //end2.Print();


  //TVector3 V(end2-end1),W(0.,0.,1.);
  //V.Unit();
  //std::cout<<"V:";
  //V.Print();

  //V*=1./V.Mag();
  //std::cout<<"V unit:";
  //V.Print();

  //TVector3 U = V.Cross(W);

  //std::cout<<"U:";
  //U.Print();
  //setDetPlane(DetPlane(TVector3(0,0,z),// position of wi in glo coord. at y=0
  //TVector3(_cangle,_sangle,0),
  //TVector3(-_sangle,_cangle,0)));
  //fPolicy.fDetPlane=GFDetPlane(end1,U,V);


  //cylHit->Print();
  //Print();
}

//------------------------------------------------------------------------------
void  HypDchRecoHit2::Print()
{
  
  cout<<"HypDchRecoHit2 :"<<getPolicyName()<<endl;
  fHitCoord.Print();
  fPolicy.Print();

  // std::cout<<"hit HMatrix:";getHMatrix().Print();
  //std::cout<<"hit DetPlane:";fDetPlane.Print();
  //   std::cout<<"\n hitCov:";getHitCov(getDetPlane(0)).Print();
  
}

TMatrixT<double> HypDchRecoHit2::getHMatrix(const GFAbsTrackRep* stateVector)
{
  // if (dynamic_cast<const GeaneTrackRep*>(stateVector) != NULL) 
  //   {
  //     TMatrixT<double> HMatrix(1,5);
      
  //     HMatrix[0][0] = 0.;
  //     HMatrix[0][1] = 0.;
  //     HMatrix[0][2] = 0.;
  //     HMatrix[0][3] = 1.;
  //     HMatrix[0][4] = 0.;
  //     return HMatrix;
  //   }
  if (dynamic_cast<const LSLTrackRep*>(stateVector) != NULL) 
    {
      // LSLTrackRep (x,y,x',y',q/p)
       TMatrixT<double> HMatrix(1,5);

      HMatrix[0][0] = 1.;//_cangle;
      HMatrix[0][1] = 0.;//_sangle;
      HMatrix[0][2] = 0.;
      HMatrix[0][3] = 0.;
      HMatrix[0][4] = 0.;
      return HMatrix;
    }
  else if(dynamic_cast<const SlTrackRep*>(stateVector) != NULL) 
    {
      // SlTrackRep (x,y,x',y')
      TMatrixT<double> HMatrix(1,4);
      
      HMatrix[0][0] = 1.;//_cangle;
      HMatrix[0][1] = 0.;//_sangle;
      HMatrix[0][2] = 0.;
      HMatrix[0][3] = 0.;
      return HMatrix;
      //     fHMatrix[1][0] = _sangle;
      //     fHMatrix[1][1] = _cangle;
      //     fHMatrix[1][2] = 0.;
      //     fHMatrix[1][3] = 0.;
      
    }
  else if ((dynamic_cast<const RKTrackRep*>(stateVector) != NULL)) {
    //I know, since this is the same everytime, it could be done in the
    //the constructor, but I do it here anyway, to make clear that in the
    //case of several track-reps per hit, it would have to be done here
    //    fHMatrix.ResizeTo(NparHitRep,5);
    TMatrixT<double> HMatrix(1,5);
    
    HMatrix[0][0] = 0.;
    HMatrix[0][1] = 0.;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 1.;
    HMatrix[0][4] = 0.;
    
    return HMatrix;
  }
  else 
    {
      std::cerr << "HypDchRecoHit2 can only handle state"
		<< " vectors of type Hyp{LSL,Sl}TrackRep -> abort" << std::endl;
      throw;
    }
 
}

