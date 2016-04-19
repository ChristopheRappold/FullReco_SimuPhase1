#include "TMcHit.hh"
#include "Riostream.h"
//#define DEBUG_MEM

using namespace std;

ClassImp(TMcHit)


TMcHit::TMcHit():name("default_hit"),LayerID(-1),HitID(-1),MCHit(-999.,-999.,-999.),Hit(-999.,-999.,-999.),MC_id(0),Charge(0),Pdg(0),Brho(0.),MagnetInteraction(0.),MCparticle()
{ 


}
// TMcHit::TMcHit(TString n,double xx,double yy,double zz,double tt,double EE,int ppid):name(n),x(xx),y(yy),z(zz),t(tt),E(EE),PID(ppid) 
// { 

// #ifdef DEBUG_MEM
//   cout<<"Hit Constructor :"<<name<<endl;
// #endif
// }

TMcHit::~TMcHit() 
{ 
#ifdef DEBUG_MEM
  cout<<"HitDet Destructor :"<<name<<endl;
#endif   
}


TMcHit::TMcHit(const TMcHit& M):LayerID(M.LayerID),HitID(M.HitID),MCHit(M.MCHit),Hit(M.Hit),MC_id(M.MC_id),Charge(M.Charge),Pdg(M.Pdg),Brho(M.Brho),MagnetInteraction(M.MagnetInteraction),MCparticle(M.MCparticle)
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

TMcHit& TMcHit::operator=(const TMcHit& M) 
{
  name = M.name;
  name+="op=";
  LayerID = M.LayerID;
  HitID = M.HitID;
  MCHit.SetXYZ(M.MCHit.X(),M.MCHit.Y(),M.MCHit.Z());
  Hit.SetXYZ(M.Hit.X(),M.Hit.Y(),M.Hit.Z());
  MC_id = M.MC_id;
  Charge = M.Charge;
  Pdg = M.Pdg;
  Brho = M.Brho;
  MagnetInteraction = M.MagnetInteraction;
  MCparticle.SetXYZT(M.MCparticle.X(),M.MCparticle.Y(),M.MCparticle.Z(),M.MCparticle.T());
    
  return *this;
}

// TBuffer &operator<<(TBuffer &buf, const TMcHit *obj)
// {
//    ((TMcHit*)obj)->Streamer(buf);
//    return buf;
// }


// TMcHit* TMcHit::CloneD(const char * name)
// {
//   const TMcHit& hit_Ref = *this;
//   TMcHit* temp = new TMcHit(hit_Ref);
//   return temp;

// }

void TMcHit::Clear(Option_t * /*option*/)
{
#ifdef DEBUG_MEM
  cout<<"HitDet Clearing:"<<name<<endl;
#endif
  name="";
  LayerID = -1;
  HitID = -1;
  MCHit.SetXYZ(-999.,-999.,-999.);
  Hit.SetXYZ(-999.,-999.,-999.);
  MC_id = 0;
  Charge = 0;
  Pdg = 0;
  Brho = 0;
  MagnetInteraction = 0.;
  MCparticle.SetXYZT(0.,0.,0.,0.);
  
}
