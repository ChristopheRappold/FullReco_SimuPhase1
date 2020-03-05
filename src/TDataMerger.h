#ifndef TDATAMERGER
#define TDATAMERGER
#include <string>
#include "Ana_Hist.hh"
#include "Debug.hh"
#include "ReturnRes.hh"

//template <typename T0,typename T1,typename T2,typename T3>

class TDataMerger
{
  public :
  const std::string signature;
  Ana_Hist* AnaHisto;

  TDataMerger():signature("merger_default"),AnaHisto(nullptr) {;}
  explicit TDataMerger(const std::string& name):signature(name),AnaHisto(nullptr) {;}
  virtual ~TDataMerger() { AnaHisto = nullptr;}

  virtual int Init(Ana_Hist* h) { AnaHisto = h; SelectHists() ; return 0;}

  //virtual int Exec(THypHi_Event *event,std::vector<TUTracker_Event*> *UTrackerEvents,FullRecoEvent& RecoEvent) = 0;
  //virtual int Exec(THyphiEvent_Phys* event,FullRecoEvent& RecoEvent) = 0;
  private :
    //virtual int Exec() = 0;
  virtual ReturnRes::InfoM SoftExit(int) = 0;
  virtual void SelectHists() = 0;
};

#endif 
