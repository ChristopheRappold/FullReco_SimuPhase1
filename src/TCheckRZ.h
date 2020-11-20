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

constexpr double dl_max(int idDet) {
  switch(idDet){
  case 0:
  case 1:
  case 2:
  case 3:
  case 4:
    return 0.2;
    break;
  case 5:
  case 6:
  case 7:
  case 8:
  case 9:
  case 10:
    return 0.3;
    break;
  case 11:
  case 12:
  case 13:
  case 14:
  case 15:
  case 16:
    return 0.4;
    break;
  default:
    return 0.;
    break;
  }
};




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
  //                                 1  2  3  4  5   6     7     8      9      10      11      12      13     14     15     16     17
  //std::vector<double> correctBias = {0.,0.,0.,0.,0.,-1.4,-0.017,-0.8852,0.069,-0.6512,0.3412,-0.1674,0.4637,0.1999,0.5361,0.3536,0.5534};
  std::vector<double> correctBias = {0.,0.,0.,0.,0.,-1.56,-0.1328,-1.013,-0.0871,-0.8192,0.1454,-0.3784,0.2278,-0.05195,0.2511,0.04834,0.2226};
  //std::vector<double> correctBias = {0.,0.,0., 0.,0.,0., 0.,0.,0., 0.,0.,0., 0.,0.,0., 0.,0.};

  bool ChangeMiniFiber = false;
  int MDCWireType = 0; // 0 -> MeasurementWire2 / 1 -> ProlateSpace
  bool MDCBiasCorr = true;

  struct LocalHists
  {
    TH1F* h_RZStats;
    TH2F* h_RZ;
    TH1F* h_RZfit_mom;
    TH2F* h_RZfit_Chi2;
    TH2F* h_XYfit_miniF;
    TH2F* h_MDC_Z_residu;
    TH2F* h_MDC_R_residu;

  };
  LocalHists LocalHisto;
};


#endif
