#ifndef TDATABUILDER
#define TDATABUILDER
#include <string>
#include "Ana_Hist.hh"
#include "Debug.hh"

//template <typename T0,typename T1,typename T2,typename T3>

class TDataBuilder
{
  public :
  const std::string signature;
  Ana_Hist* AnaHisto;

  TDataBuilder():signature("build_default"),AnaHisto(nullptr) {;}
  explicit TDataBuilder(const std::string name):signature(name),AnaHisto(nullptr) {;}
  virtual ~TDataBuilder() { AnaHisto = nullptr;}

  virtual int Init(Ana_Hist* h) { AnaHisto = h; return 0;}

  //virtual int Exec(THypHi_Event *event,std::vector<TUTracker_Event*> *UTrackerEvents,FullRecoEvent& RecoEvent) = 0;
  //virtual int Exec(THyphiEvent_Phys* event,FullRecoEvent& RecoEvent) = 0;
  private :
    //virtual int Exec() = 0;
  virtual int SoftExit(int) = 0;

};

#endif 
