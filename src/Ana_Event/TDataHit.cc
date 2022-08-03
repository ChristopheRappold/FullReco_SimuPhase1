#include "TDataHit.hh"
#include "Riostream.h"
//#define DEBUG_MEM

using namespace std;

ClassImp(TDataHit)


TDataHit::TDataHit():name("default_hit"),LayerID(-1),HitID(-1),TrackID(-1),Hit(-999.,-999.,-999.),Charge(0),Pdg(0),Brho(0.),MagnetInteraction(0.),Time(0.),Energy(0.),TrackLength(0.)
{ 

}


TDataHit::~TDataHit() 
{ 
#ifdef DEBUG_MEM
  cout<<"HitDet Destructor :"<<name<<endl;
#endif   
}


TDataHit::TDataHit(const TDataHit& M):LayerID(M.LayerID),HitID(M.HitID),TrackID(M.TrackID),Hit(M.Hit),Charge(M.Charge),Pdg(M.Pdg),Brho(M.Brho),MagnetInteraction(M.MagnetInteraction),Time(M.Time),Energy(M.Energy),TrackLength(M.TrackLength)
{
  name = "const_copy";
  name+=M.name;
#ifdef DEBUG_MEM
  cout<<"HitDet const Copy Constructor :"<<name<<endl;
#endif
  
  //HitPos = M.HitPos;
  //Id = M.Id;
  //Qdc = M.Qdc;
  //TimeCh=M.TimeCh;
  //TimeNs = M.TimeNs;
}

TDataHit& TDataHit::operator=(const TDataHit& M) 
{
  name = M.name;
  name+="op=";
  LayerID = M.LayerID;
  HitID = M.HitID;
  TrackID = M.TrackID;
  Hit.SetXYZ(M.Hit.X(),M.Hit.Y(),M.Hit.Z());
  Charge = M.Charge;
  Pdg = M.Pdg;
  Brho = M.Brho;
  MagnetInteraction = M.MagnetInteraction;
  Time = M.Time;
  Energy = M.Energy;
  TrackLength = M.TrackLength;

  return *this;
}



void TDataHit::Clear(Option_t * /*option*/)
{
#ifdef DEBUG_MEM
  cout<<"HitDet Clearing:"<<name<<endl;
#endif
  name="";
  LayerID = -1;
  HitID = -1;
  TrackID = -1;
  Hit.SetXYZ(-999.,-999.,-999.);
  Charge = 0;
  Pdg = 0;
  Brho = 0;
  MagnetInteraction = 0.;
  Time = 0.;
  Energy = 0.;
  TrackLength = 0.;
}
