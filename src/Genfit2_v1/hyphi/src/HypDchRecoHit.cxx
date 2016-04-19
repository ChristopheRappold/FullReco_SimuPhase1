// Hyp Headers 
#include "HypDchRecoHit.h"
//#include "HypDchCylinderHit.h"
#include "LSLTrackRep.h"
#include "SlTrackRep.h"
//#include "GeaneTrackRep.h"
#include "GFDetPlane.h"
#include "RKTrackRep.h"
#include "RKTrackRepXY.h"

// C/C++ Headers 
#include <assert.h>
#include "math.h"

//ROOT includes
#include "TMath.h"

ClassImp(HypDchRecoHit)

//#define DEBUG_DCH

//const std::string PolicyName ="PlanarHitPolicy";

HypDchRecoHit::~HypDchRecoHit()
{}

HypDchRecoHit::HypDchRecoHit(): PlanarRecoHit(NparHitRep)//,"PlanarHitPolicy")
{}

//HypDchRecoHit::HypDchRecoHit(double r, double wireposx, double angle, double z, double sigr): PlanarRecoHit(NparHitRep), _wirepos(wireposx)
HypDchRecoHit::HypDchRecoHit(const TVector3& end1,const TVector3& end2,double cos_a,double sin_a,double distance,double dist_sigr,bool mirroring): PlanarRecoHit(NparHitRep),mirror(mirroring)//,_sangle(sin_a),_cangle(cos_a)
{
  //_sangle=TMath::Sin(angle);
  //_cangle=TMath::Cos(angle);

  fHitCoord[0][0] = distance;
  TVector3 V(end2-end1),W(0.,0.,1.),Y(0.,1.,0.),V2(V);
  V*=1./V.Mag();
  
  TVector3 U = V.Cross(W);

  V2.SetZ(0.);
  
  double angle = U.Phi();//V2.Angle(Y);
#ifdef DEBUG_DCH
  U.Print();
#endif
  //std::cout<<"angle : "<<angle*180./TMath::Pi()<<std::endl;

  //_sangle=-1.*TMath::Sin(angle);
  _sangle=TMath::Sin(angle);
  _cangle=TMath::Cos(angle);

  //assert(TMath::Abs(_sangle-sin_a)<1e-5 && TMath::Abs(_cangle-cos_a)<1e-5);

  //setDetPlane(DetPlane(TVector3(0,0,z),// position of wi in glo coord. at y=0
  //TVector3(_cangle,_sangle,0),
  //TVector3(-_sangle,_cangle,0)));
  fPolicy.setDetPlane(GFDetPlane(end1,U,V));
  
  x_wirepos = end1.X();
  y_wirepos = end1.Y();

  fHitCov[0][0] = dist_sigr;//*dist_sigr;
}



/*HypDchRecoHit::HypDchRecoHit(const HypDchCylinderHit* cylHit): PlanarRecoHit(NparHitRep)
{
  //  std::cout<<"I create HypDchRecoHit out of a CylinderHit"<<std::endl;
  _sangle = TMath::Sin(cylHit->GetWireAngle());
  _cangle = TMath::Cos(cylHit->GetWireAngle());
  fHitCoord[0][0] = cylHit->GetDistance();
  Double_t z = cylHit->GetWireZcoordGlobal();
  _wirepos = cylHit->GetWireXcoordLocal();
  fPolicy.setDetPlane(DetPlane(TVector3(0,0,z),// position of wi in glo coord. at y=0
		       TVector3(_cangle,_sangle,0),
		       TVector3(-_sangle,_cangle,0)));
  Double_t sigr = cylHit->GetDistanceError();
  fHitCov[0][0] = sigr*sigr;

//   cylHit->Print();
//   this->Print();
}*/


GFAbsRecoHit* HypDchRecoHit::clone()
{
  return new HypDchRecoHit(*this);
}


TMatrixT<double> HypDchRecoHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
#ifdef DEBUG_DCH1
  std::cout<<"I set up the H matrix now with NparHitRep = "<<NparHitRep<<std::endl;
#endif
  assert(stateVector!=NULL);
  if (dynamic_cast<const LSLTrackRep*>(stateVector) != NULL || dynamic_cast<const RKTrackRepXY*>(stateVector) != NULL)  
  {
    // LSLTrackRep (x,y,x',y',q/p)
    // SlTrackRep (x,y,x',y')
    TMatrixT<double> HMatrix(1,5);

    HMatrix[0][0] = _cangle;
    HMatrix[0][1] = _sangle;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 0.;
    HMatrix[0][4] = 0.;

    return HMatrix;

    //     fHMatrix[1][0] = _sangle;
    //     fHMatrix[1][1] = _cangle;
    //     fHMatrix[1][2] = 0.;
    //     fHMatrix[1][3] = 0.;
    //     fHMatrix[1][4] = 0.;
  }
  else if(dynamic_cast<const SlTrackRep*>(stateVector) != NULL) 
  {
    // SlTrackRep (x,y,x',y')
    TMatrixT<double> HMatrix(1,4);

    HMatrix[0][0] = _cangle;
    HMatrix[0][1] = _sangle;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 0.;
    
    //     fHMatrix[1][0] = _sangle;
    //     fHMatrix[1][1] = _cangle;
    //     fHMatrix[1][2] = 0.;
    //     fHMatrix[1][3] = 0.;
    
    return HMatrix;
  }
  else if (dynamic_cast<const RKTrackRep*>(stateVector) != NULL) 
  {
    TMatrixT<double> HMatrix(1,5);

    HMatrix[0][0] = 0.;
    HMatrix[0][1] = 0.;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 1.;
    HMatrix[0][4] = 0.;
    //HMatrix[0][3] = _cangle;
    //HMatrix[0][4] = _sangle;

    return HMatrix;
  }
  else
  {
    std::cerr << "HypDchRecoHit can only handle state" << " vectors of type LSLTrackRep  -> abort" << std::endl;
    throw;
  }
}

TMatrixT<double> HypDchRecoHit::residualVector(const GFAbsTrackRep* stateVector, const TMatrixT<double>& state, const GFDetPlane& d)
{
#ifdef DEBUG_DCH
  std::cout<<" residu Vector HypDchRecoHit :"<<std::endl;
#endif
  
    //setHMatrix(stateVector,state);
  TMatrixT<double> H=getHMatrix(stateVector);

#ifdef DEBUG_DCH
  std::cout<<" H done ";H.Print();
#endif
  //return ( getHitCoord(d) - (_HMatrix*state ));


  int repDim=stateVector->getDim();
  TMatrixT<double> s(repDim,1);
  s[0][0]= x_wirepos;
  s[1][0]= y_wirepos;
  
  if(dynamic_cast<const RKTrackRep*>(stateVector) != NULL)
    {
      s[0][0]= 0.;
      s[1][0]= 0.;
      s[3][0]= 0.;//x_wirepos;
      s[4][0]= 0.;//y_wirepos;
    }


#ifdef DEBUG_DCH
  std::cout<<" s done "<<std::endl;
  state.Print();
  s.Print();
#endif

  TMatrixT<double> trackpos=(H*(state-s));

#ifdef DEBUG_DCH
  std::cout<<" trackpos done "<<std::endl;
#endif  

  TMatrixT<double> res = getHitCoord(d)-trackpos[0][0];
  if(mirror)
    {
      Double_t res1 = (getHitCoord(d))[0][0]-trackpos[0][0];
      Double_t res2 = -1.*(getHitCoord(d))[0][0]-trackpos[0][0];
      //TMatrixT<double> res(1,1);
      if(TMath::Abs(res1)>TMath::Abs(res2))
	res[0][0]=res2;
      else 
	res[0][0]=res1;
    }

#ifdef DEBUG_DCH
  std::cout<<" res ";res.Print();
#endif  

  return res;
}

void  HypDchRecoHit::Print()
{
  std::cout<<"hitCoord:"; fHitCoord.Print();
  //  std::cout<<"hit HMatrix:";getHMatrix().Print();
  std::cout<<"hit DetPlane:";getDetPlane(0).Print();
  std::cout<<"\n hitCov:";getHitCov(getDetPlane(0)).Print();
  
}
