#ifndef TCHECKFIELD 
#define TCHECKFIELD 

#include "TDataProcess.h"
#include "Debug.hh"
#include "FullRecoEvent.hh"

#include "Ana_Event/MCAnaEventG4Sol.hh"

#include "THyphiAttributes.h"

typedef TDataProcess<FullRecoEvent,MCAnaEventG4Sol> TDataProcessInterface;

class CheckField :  public TDataProcessInterface
{
  public :
  const THyphiAttributes& att;


  CheckField(const THyphiAttributes& attr);
  ~CheckField() final = default;

  //int Init(Ana_Hist* h);
  int operator() (FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree) final;
 private:
  int Exec(FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree) final;
  int SoftExit(int) final;
  int Check();

  bool done;
  
};


#endif
