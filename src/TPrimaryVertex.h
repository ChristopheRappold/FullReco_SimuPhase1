#ifndef TPRIMARYVERTEX
#define TPRIMARYVERTEX

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"

#include "HitAna/FiberHitXUV.hh"
#include "HitAna/FiberTrackAna.hh"
#include "HitAna/ConstantParameter.hh"
#include "HitAna/ParaManager.hh"
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

  void CloseDist_TrackTrack(FiberTrackAna* Track1, FiberTrackAna* Track2, double& distance, TVector3& midpoint);

  void CloseDist_TrackPoint(FiberTrackAna* Track, TVector3& Point, double& distance);

  double f_function(FiberTrackAna* Track, std::vector<double>& PosXYZ);

  double V_function(std::vector<double>& f_vector);

  void SpaceDiscretization(double& Xi, double& Xf, double& Yi, double& Yf, size_t& NstepsXY, double& Zi, double& Zf, size_t& NstepsZ,
                              size_t& border, std::vector<std::vector<double> >& PosXYZ);
  
  void BeamTracksFinder(std::vector<std::vector<FiberHitXUV*> > FiberXUVCont, std::vector<FiberTrackAna*>& BeamTracks_All);

  void BeamTracksSelector(std::vector<FiberTrackAna*>& BeamTracks_All, std::vector<FiberTrackAna*>& BeamTracks);

  void PrimaryTracksFinder(std::vector<std::vector<FiberHitXUV*> > FiberXUVCont, std::vector<FiberTrackAna*>& PrimaryTracks_All);

  void PrimaryTracksSelector(std::vector<FiberTrackAna*>& PrimaryTracks_All, std::vector<FiberTrackAna*>& BeamTracks,
                                  std::vector<FiberTrackAna*>& PrimaryTracks);

  void InteractionPointFinder(std::vector<FiberTrackAna*>& BeamTracks, std::vector<FiberTrackAna*>& PrimaryTracks,
                                TVector3& IP_recons, std::vector<double>& f_values_IP, size_t& i_BeamTracks_Vmax);

  void CovarianceMatrix(FiberTrackAna* BeamTrack, std::vector<FiberTrackAna*>& PrimaryTracks,
                          TVector3& IP_average, std::vector<double>& f_values_IP,
                              std::vector<std::vector<double> >& CovMatrix);


  TRandom3* rand;

  //Experimental setup parameters about target
  std::unique_ptr<ParaManager> par;
  TVector3 target_pos;
  TVector3 target_size;

  //Cut conditions on reconstructed beam tracks
  bool ifCut_MaxChi2ndf_BeamTracks = false;
  double MaxChi2ndf_BeamTracks = 3.; //Change !

  bool ifCut_MaxTheta_BeamTracks = false;
  double MaxTheta_BeamTracks = 1; //Change !

  bool ifCut_CloseDistToTgtCenter_BeamTracks = true;
  double MaxCloseDistToTgtCenter_BeamTracks = 0.7; // in target_size.Mag() units, which is 30 mm


  //Cut conditions on reconstructed primary tracks
  bool ifCut_MaxChi2ndf_PrimaryTracks = false;
  double MaxChi2ndf_PrimaryTracks = 3.; //Change !

  bool ifCut_MaxCloseDistToBeam_PrimaryTracks = false;
  double MaxCloseDistToBeam_PrimaryTracks = 10.;

  bool ifCut_CloseDistToTgtCenter_PrimaryTracks = true;
  double MaxCloseDistToTgtCenter_PrimaryTracks = 1.1; // in target_size.Mag() units, which is 30 mm

  bool ifCut_FlagPSBhit_PrimaryTracks = false;
  bool FlagPSBhit_PrimaryTracks = true; //Change !

  bool ifCut_FlagPSBEhit_PrimaryTracks = false;
  bool FlagPSBEhit_PrimaryTracks = true; //Change !

  bool ifCut_FlagPSFEhit_PrimaryTracks = false;
  bool FlagPSFEhit_PrimaryTracks = true; //Change !

  bool ifCut_FlagPShit_PrimaryTracks = false;
  bool FlagPShit_PrimaryTracks = true; //Change !

  //Parameters for space discretization
  size_t NstepsdiscretXY      = 5;
  size_t NstepsdiscretZ       = 11;
  size_t Nstepsdiscretbox     = 5;
  size_t nTimesDiscretization = 6;
  double boxDistBeamXY        = 4.;

  //Parameter for v function method
  double k_factor       = 3.; // multiplies the beam track f function value

  size_t i_BeamTracks_Vmax = 99.;


  double randInteractionPointXY = 0.1; // in mm

  //Parameters for covariance matrix calculation
  double min_f_value = 0.5;
  double sigma_FT; // Change !

  //Results for primary vertex Finder
  TVector3 IP_real;
  TVector3 IP_average;
  TVector3 IP_recons;



  struct LocalHists
  {
    TH2F* h_nTrackCandidates;
    TH2F* h_DistanceBeamTracks;
    TH2F* h_PosZBeamTracks;
    TH2F* h_thetaTracks;

    TH1F* h_fvalues;

    TH1F* h_InteractionPointDistance;
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
