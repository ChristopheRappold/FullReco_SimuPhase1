#include "THyphiSimpleHit.hh"
#include "Riostream.h"
//#define DEBUG_MEM

using namespace std;

ClassImp(THyphiSimpleHit)


THyphiSimpleHit::THyphiSimpleHit() 
{ 
  name ="default_hit";
  x = -99999.;
  y = -99999.;
  z = -99999.;
  t = -99999.;
  E = -99999.;
  PID=-1;
}
THyphiSimpleHit::THyphiSimpleHit(TString n,double xx,double yy,double zz,double tt,double EE,int ppid):name(n),x(xx),y(yy),z(zz),t(tt),E(EE),PID(ppid) 
{ 

#ifdef DEBUG_MEM
  cout<<"SimpleHit Constructor :"<<name<<endl;
#endif
}

THyphiSimpleHit::~THyphiSimpleHit() 
{ 
#ifdef DEBUG_MEM
  cout<<"SimpleHit Destructor :"<<name<<endl;
#endif   
}

THyphiSimpleHit::THyphiSimpleHit(THyphiSimpleHit& M) 
{
  name = "copy";
  name +=M.name;
#ifdef DEBUG_MEM
  cout<<"SimpleHit Copy Constructor :"<<name<<endl;
#endif
  x=M.x;
  y=M.y;
  z=M.z;
  t=M.t;
  E=M.E;
  PID=M.PID;
  
}

THyphiSimpleHit::THyphiSimpleHit(const THyphiSimpleHit& M) 
{
  name = "const_copy";
  name+=M.name;
#ifdef DEBUG_MEM
    cout<<"SimpleHit const Copy Constructor :"<<name<<endl;
#endif
  x=M.x;
  y=M.y;
  z=M.z;
  t=M.t;
  E=M.E;
  PID=M.PID;
  
}

THyphiSimpleHit& THyphiSimpleHit::operator=(const THyphiSimpleHit& M) 
{
  name = M.name;
  name+="op=";

  x=M.x;
  y=M.y;
  z=M.z;
  t=M.t;
  E=M.E;
  PID=M.PID;
  
    
   return *this;
}

TBuffer &operator<<(TBuffer &buf, const THyphiSimpleHit *obj)
{
   ((THyphiSimpleHit*)obj)->Streamer(buf);
   return buf;
}


// THyphiSimpleHit* THyphiSimpleHit::CloneD(const char * name)
// {
//   const THyphiSimpleHit& hit_Ref = *this;
//   THyphiSimpleHit* temp = new THyphiSimpleHit(hit_Ref);
//   return temp;

// }

void THyphiSimpleHit::Clear(Option_t * /*option*/)
{
#ifdef DEBUG_MEM
  cout<<"SimpleHit Clearing:"<<name<<endl;
#endif
  name="";
  x = -99999.;
  y = -99999.;
  z = -99999.;
  t = -99999.;
  E = -99999.;
  PID=-1;
  
}
