#ifndef TBAYESFINDER 
#define TBAYESFINDER 

#include "TDataProcess.h"
#include "Debug.hh"
#include "FullRecoEvent.hh"

#include "Ana_Event/MCAnaEventG4Sol.hh"

#include "THyphiAttributes.h"

typedef TDataProcess<FullRecoEvent,MCAnaEventG4Sol> TDataProcessInterface;

class TBayesFinder :  public TDataProcessInterface
{
  public :
  const THyphiAttributes& att;


  TBayesFinder(const THyphiAttributes& attr);
  ~TBayesFinder();

  //int Init(Ana_Hist* h);
  int operator() (FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree) override;
 private:
  int Exec(FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree) override;
  int SoftExit(int) override;
  int FinderTrack(FullRecoEvent& RecoEvent);

  std::vector<double> radiusCDC;
  //TVector3 Plane_time;

};


#endif
