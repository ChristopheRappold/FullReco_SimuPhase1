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
#include "TRandom3.h"
#include "TLinearFitter.h"

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

  void FindHitXUV_v1(const std::vector<genfit::AbsMeasurement*>& hitx, const std::vector<genfit::AbsMeasurement*>& hitu,
                      const std::vector<genfit::AbsMeasurement*>& hitv, int id_det);

  void FindHitXUV_v2(const std::vector<genfit::AbsMeasurement*>& hitx, const std::vector<genfit::AbsMeasurement*>& hitu,
                      const std::vector<genfit::AbsMeasurement*>& hitv, int id_det);

  void FindHitXUV_v3(const std::vector<genfit::AbsMeasurement*>& hitx, const std::vector<genfit::AbsMeasurement*>& hitu,
                      const std::vector<genfit::AbsMeasurement*>& hitv, int id_det);

  void FindHitXUV_v4(const std::vector<genfit::AbsMeasurement*>& hitx, const std::vector<genfit::AbsMeasurement*>& hitu,
                      const std::vector<genfit::AbsMeasurement*>& hitv, int id_det);

  void FindSingleHitXUVId_v1(const std::vector<genfit::AbsMeasurement*>& hitx, const std::vector<genfit::AbsMeasurement*>& hitu,
                              const std::vector<genfit::AbsMeasurement*>& hitv, int id_det, int id_track);

  void FindSingleHitXUVId_v2(const std::vector<genfit::AbsMeasurement*>& hitx, const std::vector<genfit::AbsMeasurement*>& hitu,
                              const std::vector<genfit::AbsMeasurement*>& hitv, int id_det, int id_track);

  void FindSingleHitXUVId_v3(const std::vector<genfit::AbsMeasurement*>& hitx, const std::vector<genfit::AbsMeasurement*>& hitu,
                              const std::vector<genfit::AbsMeasurement*>& hitv, int id_det, int id_track);

  void FindSingleHitXUVId_v4(const std::vector<genfit::AbsMeasurement*>& hitx, const std::vector<genfit::AbsMeasurement*>& hitu,
                              const std::vector<genfit::AbsMeasurement*>& hitv, int id_det, int id_track);

  double FindHitReal_dvalue(double hitx, double hitu, double hitv, int id_det);

  double d_function1(int id_det, double hitx, double hity);
  double d_function2(int id_det, double posx);

  std::vector<genfit::AbsMeasurement*> Clusterize(const std::vector<std::unique_ptr<genfit::AbsMeasurement>>& hit);
  std::vector<genfit::AbsMeasurement*> Clusterize(const std::vector<genfit::AbsMeasurement*>& hit);

  std::vector<std::vector<std::tuple<double, double, double> > >   vect_realCombHit = {{}, {}, {}, {}, {}, {}, {}};
  std::vector<std::vector<TVector2> >                              vect_HitXY = {{}, {}, {}, {}, {}, {}, {}};
  std::vector<std::vector<std::tuple<TVector2, int> > >            vect_SingleHitXYId = {{}, {}, {}, {}, {}, {}, {}};
  std::vector<std::vector<std::tuple<double, double, double> > >   vect_CombHit = {{}, {}, {}, {}, {}, {}, {}};
  std::vector<std::vector<std::tuple<double, double, double> > >   vect_SingleCombHit = {{}, {}, {}, {}, {}, {}, {}};
  std::vector<std::vector<std::tuple<TVector2, double, int> > >    vect_realHitXYAngId = {{}, {}, {}, {}, {}, {}, {}};
  std::vector<std::unordered_map<int,std::tuple<std::vector<size_t>,std::vector<size_t>,std::vector<size_t> > > > vect_realCombIdHit = {};


  //std::vector<double> cut_d = {0.4, 0.4, 2., 0.4, 0.4, 0.6, 0.6}; //in cm
  std::vector<double> cut_d = {1., 1., 2., 0.4, 0.4, 0.4, 0.4}; //in cm
  std::vector<double> cut_diff_d = {1., 1., 1., 0.3, 0.3, 0.35, 0.35}; //in cm

  double fiber_resolution = 0.015; // in cm
  double fiber_width = 0.11; // in cm

  double Zpos_target = 196.12; // in cm
  std::vector<double> Zpos_Fiber = {154.77, 171.17, 199.67, 226.93, 230.93, 396.00, 405.63}; //in cm

  std::vector< std::tuple<double, double, double, double, double> > param_d_funct1 = // p0, p1, p2, p3, phi_offset
               {std::make_tuple(     0.,     0.,     0.,     0.,   0.),  // UFT1
                std::make_tuple(     0.,     0.,     0.,     0.,   0.),  // UFT2
                std::make_tuple(-0.0006, 2.1169, 0.0014, -0.001,  11.),  // UFT3
                std::make_tuple(-0.0010, -0.785,  0.054,  -0.18, -30.),  // MFT1
                std::make_tuple(-0.0009, -0.786,  0.029,  -0.25,  30.),  // MFT2
                std::make_tuple( 0.0002,     2.,     0.,     0., -13.),  // DFT1
                std::make_tuple(-0.0010,   2.21,     0.,     0.,  13.)}; // DFT2

  std::vector< std::tuple<double, double, double, double> > param_d_funct2 = // p0, p1, p2, p3
               {std::make_tuple(     0.,       0.,       0.,       0.),  // UFT1
                std::make_tuple(     0.,       0.,       0.,       0.),  // UFT2
                std::make_tuple( 0.0035,  0.31200,  -0.0017,   0.0552),  // UFT3
                std::make_tuple( 0.0010, -0.02370,       0.,       0.),  // MFT1
                std::make_tuple( 0.0008, -0.02147,       0.,       0.),  // MFT2
                std::make_tuple( 0.0012,  0.01001,       0.,       0.),  // DFT1
                std::make_tuple( 0.0001,  0.00985,       0.,       0.)}; // DFT2

  std::vector<int> id_detector = {G4Sol::FiberD1_x, G4Sol::FiberD2_x, G4Sol::FiberD3_x,
                                            G4Sol::MiniFiberD1_x, G4Sol::MiniFiberD2_x,
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

  bool ifcut_MFTtheta_UFT3 = false;
  double maxMFTtheta = 18.;
  double minMFTtheta = 8.;

  TVector2 realIP;
  double randIPXY = 0.015; // in cm

  TRandom3* rand;

  //TLinearFitter *lf = new TLinearFitter(2);
  //lf->SetFormula("1 ++ x ++ y ++ xy");
  //lf->StoreData(kFALSE);

  struct LocalHists
  {
    TH1F* h_ResidualFiberHitX[7];
    TH1F* h_ResidualFiberHitY[7];
    TH1F* h_ResidualFiberHitR[7];
    TH2F* h_ResidualFiberHitXY[7];
    TH2F* h_ResidualFiberHitX_Theta[7];
    TH2F* h_ResidualFiberHitY_Theta[7];
    TH2F* h_ResidualFiberHitX_HitX[7];
    TH2F* h_ResidualFiberHitY_HitY[7];
    TH2F* h_ResidualFiberHitR_Theta[7];
    TH1F* h_ResidualSingleFiberHitX[7];
    TH1F* h_ResidualSingleFiberHitY[7];
    TH1F* h_ResidualSingleFiberHitR[7];
    TH2F* h_ResidualSingleFiberHitX_Theta[7];
    TH2F* h_ResidualSingleFiberHitY_Theta[7];
    TH2F* h_ResidualSingleFiberHitR_Theta[7];
    TH2F* h_EfficiencyFiberHit;
    TH2F* h_EfficiencySingleFiberHit;
    TH2F* h_EfficiencyFiberHit_Theta[7];
    TH2F* h_EfficiencyFiberHit_dvalue[7];
    TH2F* h_EfficiencyFiberHit_mult[7];
    TH2F* h_EfficiencySingleFiberHit_Theta[7];
    TH2F* h_EfficiencySingleFiberHit_dvalue[7];
    TH2F* h_NumFiberHit_GoodReco[7];
    TH2F* h_NumFiberHit_Ghost[7];
    TH1F* h_FiberHit_dvalue[7];
    TH1F* h_FiberHitSingle_dvalue[7];
    TH1F* h_FiberHitReal_dvalue[7];
    TH2F* h_FiberHitReal_dvalue_Theta[7];
    TH2F* h_FiberHitReal_dvalue_Phi[7];
    TH2F* h_FiberHitReal_dvalue_Theta03_Phi[7];
    TH2F* h_FiberHitReal_dvalue_Theta310_Phi[7];
    TH2F* h_FiberHitReal_dvalue_Theta1020_Phi[7];
    TH2F* h_FiberHitReal_dvalue_HitX[7];
    TH2F* h_FiberHitReal_dvalue_HitY[7];
    TH2F* h_FiberHitReal_dvalue_PosX[7];
    TH2F* h_FiberHitReal_dvalue_tanThetacosPhi[7];
    TH2F* h_FiberHitReal_dvalue_dfunction[7];
    TH1F* h_FiberHit_Residualdvalue[7];
    TH2F* h_FiberHit_Residualdvalue_Realdvalue[7];
  };
  LocalHists LocalHisto;
};


#endif
