#ifndef TPRIMARYVERTEX
#define TPRIMARYVERTEX

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"

#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "ReturnRes.hh"
#include "TDataProcess.h"
#include "THyphiAttributes.h"

#include "TVector3.h"
#include "TRandom3.h"

//typedef TDataProcess<FullRecoEvent, MCAnaEventG4Sol> TDataProcessInterface;
// template<class Out>
// struct TDataProcessInterfaceImp { using type = TDataProcess<FullRecoEvent, Out>; } ;

// template<class Out>
// using TDataProcessInterface = typename TDataProcessInterfaceImp<Out>::type;

template<class Out>
using TDataProcessInterface = TDataProcess<FullRecoEvent, Out>;


class PrimVtxTrack 
{
  public:
    PrimVtxTrack() { };
    ~PrimVtxTrack() = default;

    void SetX(double _x) { x = _x; };
    void SetY(double _y) { y = _y; };
    void SetA(double _a) { a = _a; };
    void SetB(double _b) { b = _b; };
    void SetChi2NDF(double _chi2ndf) { chi2ndf = _chi2ndf; };
    void SetFTHit(size_t i_ft, double hit_pos) { ft_hits.emplace_back(std::make_tuple(i_ft, hit_pos)); };

    double GetX() { return x; };
    double GetY() { return y; };
    double GetA() { return a; };
    double GetB() { return b; };
    double GetChi2NDF() { return chi2ndf; };
    std::vector<std::tuple<size_t,double>> GetFTHits() { return ft_hits; };

    //bool IsFTHit(size_t i_ft);
    double GetTheta();
    double GetPhi();
    size_t GetNHits() { return ft_hits.size(); };

  private:
    double x = -999.;
    double y = -999.; //z = mid of target;
    double a = -999.;
    double b = -999.;
    double chi2ndf = -999.;
    std::vector<std::tuple<size_t,double>> ft_hits = {};
};


template<class Out>
class TPrimaryVertex final : public TDataProcessInterface<Out> //TDataProcess<FullRecoEvent, Out> //TDataProcessInterface
{
public:
  const THyphiAttributes& att;

  TPrimaryVertex(const THyphiAttributes& attr);
  ~TPrimaryVertex();

  // int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator()(FullRecoEvent& RecoEvent, Out* OutTree) override;

private:
  int Exec(FullRecoEvent& RecoEvent, Out* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;

  int FinderPrimaryVertex(FullRecoEvent& RecoEvent);

  void CloseDist_TrackTrack(PrimVtxTrack& Track1, PrimVtxTrack& Track2, double& distance, TVector3& midpoint);

  double CloseDist_TrackPoint(PrimVtxTrack& Track, TVector3& Point);

  double f_function(PrimVtxTrack& Track, std::vector<double>& PosXYZ);

  double V_function(std::vector<double>& f_vector);

  void SpaceDiscretization(double& Xi, double& Xf, double& Yi, double& Yf, size_t& NstepsXY, double& Zi, double& Zf, size_t& NstepsZ,
                              size_t& border, std::vector<std::vector<double> >& PosXYZ);
  
  void BeamTracksFinder(std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > >& ListHits, std::vector<PrimVtxTrack>& BeamTracks_All);

  void BeamTracksSelector(std::vector<PrimVtxTrack>& BeamTracks_All, std::vector<PrimVtxTrack>& BeamTracks);

  void PrimaryTracksFinder(std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > >& ListHits, std::vector<PrimVtxTrack>& PrimaryTracks_All);

  void PrimaryTracksFinder_v2(std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > >& ListHits, PrimVtxTrack& BeamTrack,
                                  std::vector<PrimVtxTrack>& PrimaryTracks_All);

  void RealPrimaryTracksFinder(std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > >& ListHits,
                                std::vector<std::vector<int> >& ListHitsToTracks,
                                std::unordered_map<int, std::vector<std::vector<SimHit> > >& TrackDAFSim,
                                  std::vector<PrimVtxTrack>& PrimaryTracks_All);

  void PrimaryTracksSelector(std::vector<PrimVtxTrack>& PrimaryTracks_All, std::vector<PrimVtxTrack>& BeamTracks,
                                  std::vector<PrimVtxTrack>& PrimaryTracks);

  void InteractionPointFinder(std::vector<PrimVtxTrack>& BeamTracks, std::vector<PrimVtxTrack>& PrimaryTracks,
                                TVector3& IP_recons, std::vector<double>& f_values_IP, size_t& i_BeamTracks_Vmax);

  void CovarianceMatrix(PrimVtxTrack& BeamTrack, std::vector<PrimVtxTrack>& PrimaryTracks,
                          TVector3& IP_average, std::vector<double>& f_values_IP,
                              std::vector<std::vector<double> >& CovMatrix);


  void FindHitXUV_v3(const std::vector<genfit::AbsMeasurement*>& hitx, const std::vector<genfit::AbsMeasurement*>& hitu,
                      const std::vector<genfit::AbsMeasurement*>& hitv, int id_det);

  void FindHitXUV_v4(const std::vector<genfit::AbsMeasurement*>& hitx, const std::vector<genfit::AbsMeasurement*>& hitu,
                      const std::vector<genfit::AbsMeasurement*>& hitv, int id_det);

  void FindCombXUV_d(const std::vector<genfit::AbsMeasurement*>& hitx, const std::vector<genfit::AbsMeasurement*>& hitu,
                      const std::vector<genfit::AbsMeasurement*>& hitv, int id_det,
                        std::map<double,std::tuple<double,double,double> >& CombXUV_diffd);

  void SelectCombXUV_d(std::map<double,std::tuple<double,double,double> >& CombXUV_diffd);

  double d_function1(int id_det, double hitx, double hity);
  double d_function2(int id_det, double posx);

  void Get_dvalueThetaPhi_CombXUV(int id_det, std::tuple<double,double,double> CombXUV);

  std::vector<genfit::AbsMeasurement*> Clusterize(const std::vector<std::unique_ptr<genfit::AbsMeasurement>>& hit);
  std::vector<genfit::AbsMeasurement*> Clusterize(const std::vector<genfit::AbsMeasurement*>& hit);

  void BeamTracking(const std::vector<std::tuple<double,double,double>>& CombXUV_UFT1,
                    const std::vector<std::tuple<double,double,double>>& CombXUV_UFT2,
                      std::vector<PrimVtxTrack>& BeamTracks);

  void PrimaryTracking(const std::vector<std::tuple<double,double,double>>& CombXUV_UFT3,
                        const std::vector<std::tuple<double,double,double>>& CombXUV_MFT1,
                        const std::vector<std::tuple<double,double,double>>& CombXUV_MFT2,
                          std::vector<PrimVtxTrack>& PrimaryTracks);

  bool PrimVtxTracking(PrimVtxTrack& PrimaryVtxTrack);

  bool CloseToBeam(std::tuple<double,double,double,double> BeamEllipsePar, PrimVtxTrack& PrimaryVtxTrack);

  PrimVtxTrack TrackFitting(int nlayer, double* w, double* z, double* angle, double* s);

  bool GaussJordan( double **a, int n, double *b, int *indxr, int *indxc, int *ipiv );

  std::tuple<double,double,double,double> BeamEllipseParam(PrimVtxTrack& BeamTrack);

  std::vector<std::vector<TVector2> >                            vect_HitXY = {{}, {}, {}, {}, {}, {}, {}};
  std::vector<std::vector<std::tuple<double, double, double> > > vect_CombHit = {{}, {}, {}, {}, {}, {}, {}};

  //std::vector<double> cut_d = {0.4, 0.4, 2., 0.4, 0.4, 0.6, 0.6}; //in cm
  std::vector<double> cut_d = {0.4, 0.4, 2., 0.4, 0.4, 0.4, 0.4}; //in cm
  std::vector<double> cut_diff_d = {0.35, 0.35, 1., 0.3, 0.3, 0.35, 0.35}; //in cm

  double fiber_resolution = 0.015; // in cm
  double fiber_width = 0.11; // in cm

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

  bool SelectCombXUV_flag[7] = {true, true, false, false, false, true, true};
  bool ifcut_MFTtheta_UFT3 = false;
  double maxMFTtheta = 18.;
  double minMFTtheta = 8.;

  double randIPXY = fiber_resolution;

  TRandom3* rand;

  bool RealPrimTrack = false;

  bool onlyMFT = true;

  //Experimental setup parameters about target
  TVector3 target_pos;
  TVector3 target_size;

  //Cut conditions on reconstructed beam tracks
  bool ifCut_MaxChi2ndf_BeamTracks = false;
  double MaxChi2ndf_BeamTracks = 10.; //Change !

  bool ifCut_MaxTheta_BeamTracks = false;
  double MaxTheta_BeamTracks = 2.; //Change !

  bool ifCut_CloseDistToTgtCenter_BeamTracks = false;
  double MaxCloseDistToTgtCenter_BeamTracks = std::sqrt(2)/2. + 0.1; // in target_size.Mag() units, which is 3 cm


  //Cut conditions on reconstructed primary tracks
  bool ifCut_MaxChi2ndf_PrimaryTracks = true;
  double MaxChi2ndf_PrimaryTracks = 10.; //Change !

  bool ifCut_MaxCloseDistToBeam_PrimaryTracks = true;
  double MaxCloseDistToBeam_PrimaryTracks = 0.5;

  bool ifCut_MaxDistZBeamToTarget_PrimaryTracks = true;
  double MaxDistZBeamToTarget_PrimaryTracks = 0.6; //in target_size.Mag() units, which is 3 cm, from target_pos.Z()

  bool ifCut_MinTheta_PrimaryTracks = true;
  double MinTheta_PrimaryTracks = 5.;

  bool ifCut_MaxTheta_PrimaryTracks = true;
  double MaxTheta_PrimaryTracks = 20.;

  bool ifCut_CloseDistToTgtCenter_PrimaryTracks = true;
  double MaxCloseDistToTgtCenter_PrimaryTracks = 1.1; // in target_size.Mag() units, which is 3 cm


/*
  bool ifCut_FlagPSBhit_PrimaryTracks = false;
  bool FlagPSBhit_PrimaryTracks = true; //Change !

  bool ifCut_FlagPSBEhit_PrimaryTracks = false;
  bool FlagPSBEhit_PrimaryTracks = true; //Change !

  bool ifCut_FlagPSFEhit_PrimaryTracks = false;
  bool FlagPSFEhit_PrimaryTracks = true; //Change !

  bool ifCut_FlagPShit_PrimaryTracks = false;
  bool FlagPShit_PrimaryTracks = true; //Change !
*/

  //Parameters for space discretization
  size_t NstepsdiscretXY      = 5;
  size_t NstepsdiscretZ       = 16;
  size_t Nstepsdiscretbox     = 5;
  size_t nTimesDiscretization = 5;
  double boxDistBeamXY        = 0.4; //in cm

  //Parameter for v function method
  double k_factor       = 30.; // multiplies the beam track f function value
  bool ifCut_min_V_value = false;
  double min_V_value = 0.9; //change

  size_t i_BeamTracks_Vmax = 99.;

  //Parameters for covariance matrix calculation
  double min_f_value = 0.5;
  double sigma_FT = fiber_resolution;

  //Results for primary vertex Finder
  TVector3 IP_real;
  TVector3 IP_average;
  TVector3 IP_recons;


  struct LocalHists
  {
    TH2F* h_NFiberMult;
    TH2F* h_NCombsXUV_UFT12;
    
    TH2F* h_NCombsXUV[5];
    TH2F* h_CombsXUV_dvalue_theta[5];
    TH2F* h_CombsXUV_dvalue_phi[5];
    TH2F* h_CombsXUV_dvalue_phi_theta5[5];
    TH2F* h_CombsXUV_dvalue_phi_theta10[5];

    TH2F* h_NHits_PrimaryTracks;

    TH2F* h_nTrackCandidates;
    TH2F* h_DistanceBeamTracks;
    TH2F* h_PosZBeamTracks;
    TH2F* h_thetaTracks;
    TH2F* h_chi2ndfTracks;

    TH1F* h_fvalues;

    TH1F* h_InteractionPointPosX;
    TH1F* h_InteractionPointPosY;
    TH1F* h_InteractionPointPosZ;

    TH2F* h_InteractionPointDistance_V_value;

    TH1F* h_InteractionPointDistanceX;
    TH1F* h_InteractionPointDistanceY;
    TH1F* h_InteractionPointDistanceZ;

    TH1F* h_InteractionPointDistanceX_pull;
    TH1F* h_InteractionPointDistanceY_pull;
    TH1F* h_InteractionPointDistanceZ_pull;

    TH1F* h_CovarianceSigmaX;
    TH1F* h_CovarianceSigmaY;
    TH1F* h_CovarianceSigmaZ;

    TH1F* h_PrimVtxstats;
    TH2F* h_PrimStatus;
  };
  LocalHists LocalHisto;
};

#endif
