#ifndef TCHECKFIBERXUV
#define TCHECKFIBERXUV

#include "ReturnRes.hh"
#include "TDataProcess.h"
#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "THyphiAttributes.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"

#include "TFile.h"
//#include "TRandom3.h"
#include "TTree.h"
//#include "TH1I.h"
//#include "TH2I.h"
#include "TVector2.h"
#include "TMath.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <cmath>

template<class Out>
using TDataProcessInterface = TDataProcess<FullRecoEvent,Out>;

template<class Out>
class TCheckFiberXUV final :  public TDataProcessInterface<Out>
{
  public :
  const THyphiAttributes& att;

  TCheckFiberXUV(const THyphiAttributes& attr);
  ~TCheckFiberXUV();

  //int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator() (FullRecoEvent& RecoEvent,Out* OutTree) override;
 private:
  int Exec(FullRecoEvent& RecoEvent,Out* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;
  int CheckHitXUV(const FullRecoEvent& RecoEvent);

  void FindHitXUV(const std::vector<genfit::AbsMeasurement*>& hitx, const std::vector<genfit::AbsMeasurement*>& hitu,
                    const std::vector<genfit::AbsMeasurement*>& hitv, int id_det);

  void FindSingleHitXUVId(const std::vector<genfit::AbsMeasurement*>& hitx, const std::vector<genfit::AbsMeasurement*>& hitu,
                            const std::vector<genfit::AbsMeasurement*>& hitv, int id_det, int id_track);

  std::vector<genfit::AbsMeasurement*> Clusterize(const std::vector<std::unique_ptr<genfit::AbsMeasurement>>& hit);
  std::vector<genfit::AbsMeasurement*> Clusterize(const std::vector<genfit::AbsMeasurement*>& hit);

  std::vector<std::vector<std::tuple<double, double, double> > >   vect_realCombHit = {{}, {}, {}, {}, {}, {}, {}};
  std::vector<std::vector<TVector2> >                              vect_HitXY = {{}, {}, {}, {}, {}, {}, {}};
  std::vector<std::vector<std::tuple<TVector2, int> > >            vect_SingleHitXYId = {{}, {}, {}, {}, {}, {}, {}};
  std::vector<std::vector<std::tuple<double, double, double> > >   vect_CombHit = {{}, {}, {}, {}, {}, {}, {}};
  std::vector<std::vector<std::tuple<TVector2, double, int> > >    vect_realHitXYAngId = {{}, {}, {}, {}, {}, {}, {}};
  std::vector<std::unordered_map<int,std::tuple<std::vector<size_t>,std::vector<size_t>,std::vector<size_t> > > > vect_realCombIdHit = {};

  std::vector<double> cut_d = {};

  double fiber_resolution = 0.015; // in cm
  double fiber_width = 0.11; // in cm

  std::vector<int> id_detector = {G4Sol::FiberD1_x, G4Sol::FiberD2_x, G4Sol::FiberD3_x,
                                          G4Sol::MiniFiberD1_x1, G4Sol::MiniFiberD1_x2,
                                                    G4Sol::FiberD4_x, G4Sol::FiberD5_x};

  double ang[7][3] = {
    {  0.,  30., -30.},  // UFT1
    {  0.,  30., -30.},  // UFT2 
    {  0.,  30., -30.},  // UFT3
    {  0., -60.,  60.},  // MFT1
    {  0.,  60., -60.},  // MFT2
    { 30., -30.,   0.},  // DFT1
    {  0.,  30., -30.}}; // DFT2

  int id_mid[7] = {1, 1, 1, 1, 2, 1, 1};

  struct LocalHists
  {
    TH1F* h_ResidualFiberHitX[7];
    TH1F* h_ResidualFiberHitY[7];
    TH1F* h_ResidualFiberHitR[7];
    TH2F* h_ResidualFiberHitXY[7];
    TH2F* h_ResidualFiberHitX_Angle[7];
    TH2F* h_ResidualFiberHitY_Angle[7];
    TH2F* h_ResidualFiberHitX_HitX[7];
    TH2F* h_ResidualFiberHitY_HitY[7];
    TH2F* h_ResidualFiberHitR_Angle[7];
    TH1F* h_EfficiencyFiberHit[7];
    TH1F* h_ResidualSingleFiberHitX[7];
    TH1F* h_ResidualSingleFiberHitY[7];
    TH1F* h_ResidualSingleFiberHitR[7];
    TH2F* h_ResidualSingleFiberHitX_Angle[7];
    TH2F* h_ResidualSingleFiberHitY_Angle[7];
    TH2F* h_ResidualSingleFiberHitR_Angle[7];
  };
  LocalHists LocalHisto;
};


#endif
