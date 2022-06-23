#ifndef TCHECKFIELD 
#define TCHECKFIELD 

#include "TDataProcess.h"
#include "Debug.hh"
#include "FullRecoEvent.hh"

#include "Ana_Event/MCAnaEventG4Sol.hh"

#include "THyphiAttributes.h"

template<class Out>
using TDataProcessInterface = TDataProcess<FullRecoEvent,Out>;

template<class Out>
class CheckField final :  public TDataProcessInterface<Out>
{
  public :
  const THyphiAttributes& att;


  CheckField(const THyphiAttributes& attr);
  ~CheckField() final = default;

  void InitMT() final;
  ReturnRes::InfoM operator() (FullRecoEvent& RecoEvent,Out* OutTree) final;
 private:
  int Exec(FullRecoEvent& RecoEvent,Out* OutTree) final;
  ReturnRes::InfoM SoftExit(int) final;
  void SelectHists() final;
  int Check();

  bool done;

  struct LocalHists
  {
    TH2F* FieldXY[3];
    TH2F* FieldXZ[3];
    TH2F* FieldYZ[3];
    TH2F* FieldXYmax[3];
    TH2F* FieldXZmax[3];
    TH2F* FieldYZmax[3];
    TH2F* FieldXY_n[3];
    TH2F* FieldXZ_n[3];
    TH2F* FieldYZ_n[3];
  };
  LocalHists LocalHisto;

};


#endif
