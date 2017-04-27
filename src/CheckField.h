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
  ~CheckField();

  //int Init(Ana_Hist* h);
  int operator() (FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree) override;
 private:
  int Exec(FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree) override;
  int SoftExit(int) override;
  int Check();

  bool done;
  
};


#endif
