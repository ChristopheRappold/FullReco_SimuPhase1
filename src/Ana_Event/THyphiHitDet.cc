#include "THyphiHitDet.hh"
#include "Riostream.h"
//#define DEBUG_MEM

using namespace std;

ClassImp(THyphiHitDet)


THyphiHitDet::THyphiHitDet() 
{ 
  name ="default_hit";

  HitPos.SetXYZ(-999,-999,-999);
  Id.clear();
  Qdc.clear();
  TimeCh.clear();
  TimeNs=-999;

}
// THyphiHitDet::THyphiHitDet(TString n,double xx,double yy,double zz,double tt,double EE,int ppid):name(n),x(xx),y(yy),z(zz),t(tt),E(EE),PID(ppid) 
// { 

// #ifdef DEBUG_MEM
//   cout<<"Hit Constructor :"<<name<<endl;
// #endif
// }

THyphiHitDet::~THyphiHitDet() 
{ 
#ifdef DEBUG_MEM
  cout<<"HitDet Destructor :"<<name<<endl;
#endif   
}

THyphiHitDet::THyphiHitDet(THyphiHitDet& M) 
{
  name = "copy";
  name +=M.name;
#ifdef DEBUG_MEM
  cout<<"HitDet Copy Constructor :"<<name<<endl;
#endif
  
}

THyphiHitDet::THyphiHitDet(const THyphiHitDet& M):HitPos(M.HitPos),Id(M.Id),Qdc(M.Qdc),TimeCh(M.TimeCh),TimeNs(M.TimeNs)
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

THyphiHitDet& THyphiHitDet::operator=(const THyphiHitDet& M) 
{
  name = M.name;
  name+="op=";
  HitPos = M.HitPos;
  Id = M.Id;
  Qdc = M.Qdc;
  TimeCh=M.TimeCh;
  TimeNs = M.TimeNs;
    
  return *this;
}

TBuffer &operator<<(TBuffer &buf, const THyphiHitDet *obj)
{
   ((THyphiHitDet*)obj)->Streamer(buf);
   return buf;
}


// THyphiHitDet* THyphiHitDet::CloneD(const char * name)
// {
//   const THyphiHitDet& hit_Ref = *this;
//   THyphiHitDet* temp = new THyphiHitDet(hit_Ref);
//   return temp;

// }

void THyphiHitDet::Clear(Option_t * /*option*/)
{
#ifdef DEBUG_MEM
  cout<<"HitDet Clearing:"<<name<<endl;
#endif
  name="";
  HitPos.SetXYZ(-999,-999,-999);
  Id.clear();
  Qdc.clear();
  TimeCh.clear();
  TimeNs=-999;
  
}
