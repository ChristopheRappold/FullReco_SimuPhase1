#ifndef TDATAPROCESS
#define TDATAPROCESS
#include <string>
#include "Ana_Hist.hh"
#include "Debug.hh"
#include "ReturnRes.hh"

template <typename T0,typename T1>
class TDataProcess 
{
  public :
  const std::string signature;
  Ana_Hist* AnaHisto;

  TDataProcess():signature("proc_default"),AnaHisto(nullptr) {;}
  explicit TDataProcess(const std::string& name):signature(name),AnaHisto(nullptr) {;}
  virtual ~TDataProcess() { AnaHisto = nullptr;}

  virtual int Init(Ana_Hist* h) { AnaHisto = h; SelectHists(); return 0;} 
  virtual void InitMT() = 0;
  virtual ReturnRes::InfoM operator() (T0& t1,T1* t2) = 0;
  private :
  virtual int Exec(T0& t1, T1* t2) = 0;
  
  virtual ReturnRes::InfoM SoftExit(int) = 0;
  virtual void SelectHists() = 0;
  
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
