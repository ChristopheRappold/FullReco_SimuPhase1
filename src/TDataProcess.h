#ifndef TDATAPROCESS
#define TDATAPROCESS
#include <string>
#include "Ana_Hist.hh"
#include "Debug.hh"

template <typename T0,typename T1>
class TDataProcess 
{
  public :
  const std::string signature;
  Ana_Hist* AnaHisto;

  TDataProcess():signature("proc_default"),AnaHisto(NULL) {;}
  TDataProcess(const std::string name):signature(name),AnaHisto(NULL) {;}
  virtual ~TDataProcess() { AnaHisto = 0;}

  virtual int Init(Ana_Hist* h) { AnaHisto = h; return 0;} 
  virtual int operator() (T0& t1,T1* t2) = 0;
  private :
  virtual int Exec(T0& t1, T1* t2) = 0;

  virtual int SoftExit(int) = 0;

};
/*
template <typename T0,typename T1>
class TDataProcessTR0Finding : TDataProcess<T0,T1>
{
public :

  TDataProcessTR0Finding():TDataProcess<T0,T1>("TR0_find_default") {;}
  ~TDataProcessTR0Finding() {;}

}

template <typename T0,typename T1>
class TDataProcessPreTracking : TDataProcess<T0,T1>
{
public :

  TDataProcessPreTracking():TDataProcess<T0,T1>("PreTracking_default") {;}
  ~TDataProcessPreTracking() {;}

}

template <typename T0,typename T1>
class TDataProcessTrackFinding : TDataProcess<T0,T1>
{
public :

  TDataProcessTrackFinding():TDataProcess<T0,T1>("TrackFinding_default") {;}
  ~TDataProcessTrackFinding() {;}

}

template <typename T0,typename T1>
class TDataProcessTrackFitting : TDataProcess<T0,T1>
{
public :

  TDataProcessTrackFitting():TDataProcess<T0,T1>("TrackFitting_default") {;}
  ~TDataProcessTrackFitting() {;}

}

template <typename T0,typename T1>
class TDataProcessPhys : TDataProcess<T0,T1>
{
public :

  TDataProcessPhys():TDataProcess<T0,T1>("Phys_default") {;}
  ~TDataProcessPhys() {;}

}*/

#endif
