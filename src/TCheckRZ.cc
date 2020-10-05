#include "TCheckRZ.h"

#include <tuple>

using namespace std;
using namespace G4Sol;

TCheckRZ::TCheckRZ(const THyphiAttributes& attribut):TDataProcessInterface("Check_R_Z"),att(attribut)
{

}

TCheckRZ::~TCheckRZ()
{

}

int TCheckRZ::operator() (FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree)
{

  int result_finder = Exec(RecoEvent,OutTree);

  return SoftExit(result_finder);
}

int TCheckRZ::Exec(FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree)
{
  return FinderTrack(RecoEvent);
}

int TCheckRZ::SoftExit(int result_full)
{
  return ReturnRes::Fine;
}


void TCheckRZ::SelectHists()
{

  LocalHisto.h_RZ           = AnaHisto->CloneAndRegister(AnaHisto->h_RK);
  LocalHisto.h_MDC_Z_residu = AnaHisto->CloneAndRegister(AnaHisto->h_MDC_Z_residu);

}



int TCheckRZ::FinderTrack(FullRecoEvent& RecoEvent)
{



  
  return 0;
}
