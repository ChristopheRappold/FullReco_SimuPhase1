#ifndef TCHECKRZ
#define TCHECKRZ

#include "ReturnRes.hh"
#include "TDataProcess.h"
#include "Debug.hh"
#include "FullRecoEvent.hh"

#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TRandom3.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"

#include "THyphiAttributes.h"

typedef TDataProcess<FullRecoEvent,MCAnaEventG4Sol> TDataProcessInterface;

class TCheckRZ final :  public TDataProcessInterface
{
  public :
  const THyphiAttributes& att;

  TCheckRZ(const THyphiAttributes& attr);
  ~TCheckRZ();

  //int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator() (FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree) override;
 private:
  int Exec(FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;

  int FinderTrack(FullRecoEvent& RecoEvent);

  struct LocalHists
  {
    TH2F* h_RZ;
    TH2F* h_MDC_Z_residu;
  };
  LocalHists LocalHisto;
};


#endif
