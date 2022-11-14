#include "TPrimaryVertex.h"
/*
#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"
#include "FullRecoEvent.hh"
#include "ReturnRes.hh"

#include "TVector3.h"
*/


//#define DEBUG_PRIMVTX
#define VERTEX_RECONS_CHECK
#define COVARIANCE_MATRIX

using namespace std;
using namespace G4Sol;

template <class Out>
TPrimaryVertex<Out>::TPrimaryVertex(const THyphiAttributes& attribut)
    : TDataProcessInterface<Out>("PrimaryVertexReco"), att(attribut)
{
  att._logger->info("TPrimaryVertex::TPrimaryVertex");

  rand = new TRandom3();
 
  par = std::make_unique<ParaManager>(att.map_ParamFiles);
  target_pos.SetXYZ(par->fiber_tgt_pos_x, par->fiber_tgt_pos_y, par->fiber_tgt_pos_z);
  target_size.SetXYZ(par->fiber_tgt_size_x, par->fiber_tgt_size_y, par->fiber_tgt_size_z);
  sigma_FT = par->fiber_res;
}


template <class Out>
TPrimaryVertex<Out>::~TPrimaryVertex() {}


template <class Out>
void TPrimaryVertex<Out>::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }


template <class Out>
ReturnRes::InfoM TPrimaryVertex<Out>::operator()(FullRecoEvent& RecoEvent, Out* OutTree)
{

  int result_finder = Exec(RecoEvent, OutTree);

  return SoftExit(result_finder);
}


template <class Out>
int TPrimaryVertex<Out>::Exec(FullRecoEvent& RecoEvent, Out* OutTree) { return FinderPrimaryVertex(RecoEvent); }


template <class Out>
ReturnRes::InfoM TPrimaryVertex<Out>::SoftExit(int result_full) {
   
  if(result_full == -1)
    {
      att._logger->debug("No enough Beam tracks for primary vertex recons");
      LocalHisto.h_PrimVtxstats->Fill("BeamTracks=0", 1.);
      return ReturnRes::Fine;
    }

  else if(result_full == -2)
    {
      att._logger->debug("No enough Primary tracks for primary vertex recons");
      LocalHisto.h_PrimVtxstats->Fill("PrimaryTracks=0", 1.);
      return ReturnRes::Fine;
    }

  else if(result_full == -3)
    {
      att._logger->debug("No simulated hypernucleus");
      LocalHisto.h_PrimVtxstats->Fill("Mother=0", 1.);
      //return ReturnRes::PrimVtxError;
      return ReturnRes::Fine;
    }

  LocalHisto.h_PrimVtxstats->Fill("Fine", 1.);

  return ReturnRes::Fine; 
}


template <class Out>
void TPrimaryVertex<Out>::SelectHists()
{
  LocalHisto.h_nTrackCandidates   = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_nTrackCandidates);
  LocalHisto.h_DistanceBeamTracks = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DistanceBeamTracks);
  LocalHisto.h_PosZBeamTracks     = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_PosZBeamTracks);
  LocalHisto.h_thetaTracks        = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_thetaTracks);

  LocalHisto.h_fvalues           = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_fvalues);

  LocalHisto.h_InteractionPointDistance  = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_InteractionPointDistance);
  LocalHisto.h_InteractionPointDistanceX = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_InteractionPointDistanceX);
  LocalHisto.h_InteractionPointDistanceY = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_InteractionPointDistanceY);
  LocalHisto.h_InteractionPointDistanceZ = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_InteractionPointDistanceZ);

  LocalHisto.h_InteractionPointDistanceX_pull = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_InteractionPointDistanceX_pull);
  LocalHisto.h_InteractionPointDistanceY_pull = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_InteractionPointDistanceY_pull);
  LocalHisto.h_InteractionPointDistanceZ_pull = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_InteractionPointDistanceZ_pull);

  LocalHisto.h_CovarianceSigmaX = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_CovarianceSigmaX);
  LocalHisto.h_CovarianceSigmaY = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_CovarianceSigmaY);
  LocalHisto.h_CovarianceSigmaZ = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_CovarianceSigmaZ);

  LocalHisto.h_PrimVtxstats = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_PrimVtxstats);
  LocalHisto.h_PrimStatus = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_PrimStatus);
}

template <class Out>
int TPrimaryVertex<Out>::FinderPrimaryVertex(FullRecoEvent& RecoEvent)
{

  LocalHisto.h_PrimStatus->Fill("PrimVtx_acc", 1., 1);

  IP_real.SetXYZ(RecoEvent.InteractionPoint[0], RecoEvent.InteractionPoint[1], RecoEvent.InteractionPoint[2]);

  // Exp Beam Tracks finding

  std::vector<FiberTrackAna*> BeamTracks_All;
  BeamTracksFinder(RecoEvent.FiberXUVCont, BeamTracks_All);

  LocalHisto.h_nTrackCandidates->Fill(BeamTracks_All.size(), "BeamTracks_All", 1.);
  for(size_t i = 0; i < BeamTracks_All.size(); ++i)
    {
      double Theta_BeamTracks = std::sqrt( std::pow(BeamTracks_All[i]->GetA(),2) + std::pow(BeamTracks_All[i]->GetB(),2) );
      LocalHisto.h_thetaTracks->Fill(Theta_BeamTracks, "BeamTracks_All", 1.);
    }

  std::vector<FiberTrackAna*> BeamTracks;
  BeamTracksSelector(BeamTracks_All, BeamTracks);

  LocalHisto.h_nTrackCandidates->Fill(BeamTracks.size(), "BeamTracks", 1.);

  if(BeamTracks.size() == 0)
    return -1;


  // Primary Tracks finding

  std::vector<FiberTrackAna*> PrimaryTracks_All;
  PrimaryTracksFinder(RecoEvent.FiberXUVCont, PrimaryTracks_All);

  LocalHisto.h_nTrackCandidates->Fill(PrimaryTracks_All.size(), "PrimaryTracks_All", 1.);
  for(size_t i = 0; i < PrimaryTracks_All.size(); ++i)
    {
      double Theta_PrimaryTracks = std::sqrt( std::pow(PrimaryTracks_All[i]->GetA(),2) + std::pow(PrimaryTracks_All[i]->GetB(),2) );
      LocalHisto.h_thetaTracks->Fill(Theta_PrimaryTracks, "PrimaryTracks_All", 1.);

      for(size_t j = 0; j < BeamTracks.size(); ++j)
        {
          double temp_CloseDistToBeam_PrimaryTracks;
          TVector3 temp_midpoint;
          CloseDist_TrackTrack(PrimaryTracks_All[i], BeamTracks[j], temp_CloseDistToBeam_PrimaryTracks, temp_midpoint);

          LocalHisto.h_DistanceBeamTracks->Fill(temp_CloseDistToBeam_PrimaryTracks, "PrimaryTracks_All", 1.);
          LocalHisto.h_PosZBeamTracks->Fill(temp_midpoint.Z(), "PrimaryTracks_All", 1.);
        }
    }
  
  std::vector<FiberTrackAna*> PrimaryTracks;
  PrimaryTracksSelector(PrimaryTracks_All, BeamTracks, PrimaryTracks);

  LocalHisto.h_nTrackCandidates->Fill(PrimaryTracks.size(), "PrimaryTracks", 1.);
  if(PrimaryTracks.size() == 0)
    return -2;


  std::vector<double> f_values_IP(PrimaryTracks.size() + 1, 0.);

  InteractionPointFinder(BeamTracks, PrimaryTracks, IP_recons, f_values_IP, i_BeamTracks_Vmax);
  RecoEvent.PrimVtxRecons = IP_recons;

#ifdef VERTEX_RECONS_CHECK

  for(size_t i = 0; i < f_values_IP.size(); ++i)
    LocalHisto.h_fvalues->Fill(f_values_IP[i], 1.);

  TVector3 IP_res(IP_real - IP_recons);

  LocalHisto.h_InteractionPointDistance->Fill(IP_res.Mag(), 1.);
  LocalHisto.h_InteractionPointDistanceX->Fill(IP_res.X(), 1.);
  LocalHisto.h_InteractionPointDistanceY->Fill(IP_res.Y(), 1.);
  LocalHisto.h_InteractionPointDistanceZ->Fill(IP_res.Z(), 1.);

#endif

#ifdef COVARIANCE_MATRIX

  std::vector<std::vector<double> > CovMatrix;
  CovarianceMatrix(BeamTracks[i_BeamTracks_Vmax], PrimaryTracks, IP_average, f_values_IP, CovMatrix);

  RecoEvent.CovMatrix_IP = {CovMatrix[0][0],
                            CovMatrix[1][0], CovMatrix[1][1],
                            CovMatrix[2][0], CovMatrix[2][1], CovMatrix[2][2]};

  LocalHisto.h_InteractionPointDistanceX_pull->Fill(IP_res.X() / std::sqrt(TMath::Max(CovMatrix[0][0], 1e-10)), 1.);
  LocalHisto.h_InteractionPointDistanceY_pull->Fill(IP_res.Y() / std::sqrt(TMath::Max(CovMatrix[1][1], 1e-10)), 1.);
  LocalHisto.h_InteractionPointDistanceZ_pull->Fill(IP_res.Z() / std::sqrt(TMath::Max(CovMatrix[2][2], 1e-10)), 1.);


  LocalHisto.h_CovarianceSigmaX->Fill(std::sqrt(CovMatrix[0][0]), 1.);
  LocalHisto.h_CovarianceSigmaY->Fill(std::sqrt(CovMatrix[1][1]), 1.);
  LocalHisto.h_CovarianceSigmaZ->Fill(std::sqrt(CovMatrix[2][2]), 1.);

#endif


  return 0;
}



template <class Out>
void TPrimaryVertex<Out>::CloseDist_TrackTrack(FiberTrackAna* Track1, FiberTrackAna* Track2,
                                      double& distance, TVector3& midpoint)
{
  std::vector<double> Track1Hit{Track1->GetXtgt(), Track1->GetYtgt(), target_pos.Z()};
  std::vector<double> Track1Par{Track1->GetA(), Track1->GetB(), 1.};

  std::vector<double> Track2Hit{Track2->GetXtgt(), Track2->GetYtgt(), target_pos.Z()};
  std::vector<double> Track2Par{Track2->GetA(), Track2->GetB(), 1.};

  std::vector<double> n{Track1Par[1] * Track2Par[2] - Track1Par[2] * Track2Par[1],
                        Track1Par[2] * Track2Par[0] - Track1Par[0] * Track2Par[2],
                        Track1Par[0] * Track2Par[1] - Track1Par[1] * Track2Par[0]};
  std::vector<double> n1{Track1Par[1] * n[2] - Track1Par[2] * n[1],
                         Track1Par[2] * n[0] - Track1Par[0] * n[2],
                         Track1Par[0] * n[1] - Track1Par[1] * n[0]};
  std::vector<double> n2{Track2Par[1] * n[2] - Track2Par[2] * n[1],
                         Track2Par[2] * n[0] - Track2Par[0] * n[2],
                         Track2Par[0] * n[1] - Track2Par[1] * n[0]};

  std::vector<double> c1 = {-999., -999., -999.};
  std::vector<double> c2 = { -99.,  -99.,  -99.};

  for(size_t i = 0; i < 3; ++i)
    {
      c1[i] = Track1Hit[i] + ((Track2Hit[1] - Track1Hit[1]) * n2[0] + (Track2Hit[2] - Track1Hit[2]) * n2[1] +
                              (Track2Hit[3] - Track1Hit[3]) * n2[2]) /
                                  (Track1Par[0] * n2[0] + Track1Par[1] * n2[1] + Track1Par[2] * n2[2]) * Track1Par[i];
      c2[i] = Track2Hit[i] + -((Track2Hit[1] - Track1Hit[1]) * n1[0] + (Track2Hit[2] - Track1Hit[2]) * n1[1] +
                               (Track2Hit[3] - Track1Hit[3]) * n1[2]) /
                                  (Track2Par[0] * n1[0] + Track2Par[1] * n1[1] + Track2Par[2] * n1[2]) * Track2Par[i];
    }

  distance = std::sqrt((c2[0] - c1[0]) * (c2[0] - c1[0]) + (c2[1] - c1[1]) * (c2[1] - c1[1]) + (c2[2] - c1[2]) * (c2[2] - c1[2]));
  midpoint.SetXYZ( (c1[0] + c2[0])/2., (c1[1] + c2[1])/2. , (c1[2] + c2[2])/2. );

  return;
}


template <class Out>
void TPrimaryVertex<Out>::CloseDist_TrackPoint(FiberTrackAna* Track, TVector3& Point, double& distance)
{
  TVector3 TrackHit(Track->GetXtgt(), Track->GetYtgt(), target_pos.Z());
  TVector3 TrackPar(Track->GetA(), Track->GetB(), 1.);

  TVector3 TrackHit_Point(Point - TrackHit);
  TVector3 TrackHit_Point_x_TrackPar(TrackHit_Point.Cross(TrackPar));

  double mag_TrackHit_Point_x_TrackPar = TrackHit_Point_x_TrackPar.Mag();
  double mag_TrackPar = TrackPar.Mag();

  distance = mag_TrackHit_Point_x_TrackPar / mag_TrackPar;

  return;
}


template <class Out>
double TPrimaryVertex<Out>::f_function(FiberTrackAna* Track, std::vector<double>& PosXYZ)
{
  double slope_x     = Track->GetA();
  double intercept_x = Track->GetXtgt() - slope_x * target_pos.Z();
  double slope_y     = Track->GetB();
  double intercept_y = Track->GetYtgt() - slope_y * target_pos.Z();

  double distanceStepXY = 2. * boxDistBeamXY / static_cast<double>(NstepsdiscretXY - 1);
  double sigma2         = std::pow(distanceStepXY, 2.) / 12.;

  double f = std::exp(-0.5 *
                 (std::pow((PosXYZ[0] - slope_x * PosXYZ[2] - intercept_x), 2.) +
                  std::pow((PosXYZ[1] - slope_y * PosXYZ[2] - intercept_y), 2.)) /
                 sigma2);
  return f;
}


template <class Out>
double TPrimaryVertex<Out>::V_function(std::vector<double>& f_vector)
{
  double sum_f  = 0;
  double sum_f2 = 0;
  double v      = 0.;

  for(size_t i = 0; i < f_vector.size() - 1; ++i)
    {
      sum_f += f_vector[i];
      sum_f2 += std::pow(f_vector[i], 2.);
    }

  if((sum_f > 1.E-9) && (sum_f2 > 1.E-9))
    {
      v = k_factor * f_vector[-1] + sum_f -
          (std::pow(k_factor * f_vector[-1], 2.) + sum_f2) /
              (k_factor * f_vector[-1] + sum_f);
    }

  return v;
}


template <class Out>
void TPrimaryVertex<Out>::SpaceDiscretization(double& Xi, double& Xf, double& Yi, double& Yf, size_t& NstepsXY,
                                                double& Zi, double& Zf, size_t& NstepsZ, size_t& border,
                                                    std::vector<std::vector<double> >& PosXYZ)
{
  // if border == 0 -> Borders included
  // if border == 1 -> Borders not included

  if((border != 0) && (border != 1))
    {
      border = 0;
      cout << "Problem with parameter border (must be 0 or 1) -> Solved by including borders\n";
    }

  PosXYZ = {};

  for(size_t i = 0 + border; i < NstepsXY - border; ++i)
    {
      double PosX = Xi + i * (Xf - Xi) / static_cast<double>(NstepsXY - 1);

      for(size_t j = 0 + border; j < NstepsXY - border; ++j)
        {
          double PosY = Yi + j * (Yf - Yi) / static_cast<double>(NstepsXY - 1);

          for(size_t k = 0 + border; k < NstepsZ - border; ++k)
            {
              double PosZ = Zi + k * (Zf - Zi) / static_cast<double>(NstepsZ - 1);

              std::vector<double> temp_PosXYZ{PosX, PosY, PosZ};
              PosXYZ.emplace_back(temp_PosXYZ);
            }
        }
    }

  return;
}

template <class Out>
void TPrimaryVertex<Out>::BeamTracksFinder(std::vector<std::vector<FiberHitXUV*> > FiberXUVCont,
                                              std::vector<FiberTrackAna*>& BeamTracks_All)
{

  for(size_t i = 0; i < FiberXUVCont[0].size(); ++i)
    {
      for(size_t j = 0; j < FiberXUVCont[1].size(); ++j)
        {
            std::vector<FiberHitXUV*> buf_xuv;

            buf_xuv.emplace_back(FiberXUVCont[0][i]);
            buf_xuv.emplace_back(FiberXUVCont[1][j]);

            FiberTrackAna *track = new FiberTrackAna(buf_xuv, par.get());
            BeamTracks_All.emplace_back(track);
        }
    }

  return;
}


template <class Out>
void TPrimaryVertex<Out>::BeamTracksSelector(std::vector<FiberTrackAna*>& BeamTracks_All,
                                              std::vector<FiberTrackAna*>& BeamTracks)
{

/*
  Possible cuts:
  - upper limit to chi2/ndf
  - upper limit to theta angle
  - upper limit to distance to target center == flag inside target
*/

  for(size_t i = 0; i < BeamTracks_All.size(); ++i)
    {
      double Chi2ndf_BeamTracks = BeamTracks_All[i]->GetChi2();
      if( ifCut_MaxChi2ndf_BeamTracks && (Chi2ndf_BeamTracks > MaxChi2ndf_BeamTracks) )
        continue;

      double Theta_BeamTracks = std::sqrt( std::pow(BeamTracks_All[i]->GetA(),2) + std::pow(BeamTracks_All[i]->GetB(),2) );

      LocalHisto.h_thetaTracks->Fill(Theta_BeamTracks, "BeamTracks_All", 1.);

      if( ifCut_MaxTheta_BeamTracks && (Theta_BeamTracks > MaxTheta_BeamTracks) )
        continue;

      double CloseDistToTgtCenter_BeamTracks;
      CloseDist_TrackPoint(BeamTracks_All[i], target_pos, CloseDistToTgtCenter_BeamTracks);
      if( ifCut_CloseDistToTgtCenter_BeamTracks && (CloseDistToTgtCenter_BeamTracks > target_size.Mag()*MaxCloseDistToTgtCenter_BeamTracks) )
        continue;

      BeamTracks.emplace_back(BeamTracks_All[i]);
    }

  return;
}


template <class Out>
void TPrimaryVertex<Out>::PrimaryTracksFinder(std::vector<std::vector<FiberHitXUV*> > FiberXUVCont,
                                              std::vector<FiberTrackAna*>& PrimaryTracks_All)
{

  for(size_t i = 0; i < FiberXUVCont[3].size(); ++i)
    {
      for(size_t j = 0; j < FiberXUVCont[4].size(); ++j)
        {
          for(size_t k = 0; k < FiberXUVCont[2].size(); ++k)
            {
              std::vector<FiberHitXUV*> buf_xuv;

              buf_xuv.emplace_back(FiberXUVCont[2][k]);
              buf_xuv.emplace_back(FiberXUVCont[3][i]);
              buf_xuv.emplace_back(FiberXUVCont[4][j]);

              FiberTrackAna *track = new FiberTrackAna(buf_xuv, par.get());
              PrimaryTracks_All.emplace_back(track);
            }
        }
    }

  return;
}


template <class Out>
void TPrimaryVertex<Out>::PrimaryTracksSelector(std::vector<FiberTrackAna*>& PrimaryTracks_All, std::vector<FiberTrackAna*>& BeamTracks,
                                                    std::vector<FiberTrackAna*>& PrimaryTracks)
{

/*
  Possible cuts:
  - upper limit to chi2/ndf
  - upper limit to distance to beam track
  - upper limit to distance to target center == flag inside target
  - flag PSB hit
  - flag PSFE hit
  - flag PSBE hit
  - flag PS hit
*/

  for(size_t i = 0; i < PrimaryTracks_All.size(); ++i)
    {
      double Chi2ndf_PrimaryTracks = PrimaryTracks_All[i]->GetChi2();
      if( ifCut_MaxChi2ndf_PrimaryTracks && (Chi2ndf_PrimaryTracks > MaxChi2ndf_PrimaryTracks) )
        continue;

      double CloseDistToBeam_PrimaryTracks = 9999.;
      for(size_t j = 0; j < BeamTracks.size(); ++j)
        {
          double temp_CloseDistToBeam_PrimaryTracks;
          TVector3 temp_midpoint;
          CloseDist_TrackTrack(PrimaryTracks_All[i], BeamTracks[j], temp_CloseDistToBeam_PrimaryTracks, temp_midpoint);

          if(temp_CloseDistToBeam_PrimaryTracks < CloseDistToBeam_PrimaryTracks)
            CloseDistToBeam_PrimaryTracks = temp_CloseDistToBeam_PrimaryTracks;
        }
      if( ifCut_MaxCloseDistToBeam_PrimaryTracks && (CloseDistToBeam_PrimaryTracks > MaxCloseDistToBeam_PrimaryTracks) )
        continue;

      double CloseDistToTgtCenter_PrimaryTracks;
      CloseDist_TrackPoint(PrimaryTracks_All[i], target_pos, CloseDistToTgtCenter_PrimaryTracks);
      if( ifCut_CloseDistToTgtCenter_PrimaryTracks && (CloseDistToTgtCenter_PrimaryTracks > target_size.Mag()*MaxCloseDistToTgtCenter_PrimaryTracks) )
        continue;

      if (ifCut_FlagPSBhit_PrimaryTracks && (FlagPSBhit_PrimaryTracks == PrimaryTracks_All[i]->IsFlagPSB()))
        continue;

      if (ifCut_FlagPSFEhit_PrimaryTracks && (FlagPSFEhit_PrimaryTracks == PrimaryTracks_All[i]->IsFlagPSFE()))
        continue;

      if (ifCut_FlagPSBEhit_PrimaryTracks && (FlagPSBEhit_PrimaryTracks == PrimaryTracks_All[i]->IsFlagPSBE()))
        continue;

      if (ifCut_FlagPShit_PrimaryTracks && (FlagPShit_PrimaryTracks ==
                                             (PrimaryTracks_All[i]->IsFlagPSB() || PrimaryTracks_All[i]->IsFlagPSFE() || PrimaryTracks_All[i]->IsFlagPSBE() ) ) )
        continue;

      PrimaryTracks.emplace_back(PrimaryTracks_All[i]);
    }

  return;
}


template <class Out>
void TPrimaryVertex<Out>::InteractionPointFinder(std::vector<FiberTrackAna*>& BeamTracks, std::vector<FiberTrackAna*>& PrimaryTracks,
                                                  TVector3& InteractionPointRecons, std::vector<double>& f_values_IP, size_t& i_BeamTracks_Vmax)
{
  double V_total = 0.;

  for(size_t l = 0; l < BeamTracks.size(); ++l)
    {
      TVector3 temp_InteractionPointRecons(-99., -99., -99.);

      std::vector<double> temp_f_values_IP(PrimaryTracks.size() + 1, 0.);
      std::vector<double> temp_f(PrimaryTracks.size() + 1, 0.);

      double V_temp = 0.;
      double V_new  = 0.;

      double Xi            = BeamTracks[l]->GetXtgt() - boxDistBeamXY;
      double Xf            = BeamTracks[l]->GetXtgt() + boxDistBeamXY;
      double distanceStepX = (Xf - Xi) / static_cast<double>(NstepsdiscretXY - 1);

      double Yi            = BeamTracks[l]->GetYtgt() - boxDistBeamXY;
      double Yf            = BeamTracks[l]->GetYtgt() + boxDistBeamXY;
      double distanceStepY = (Yf - Yi) / static_cast<double>(NstepsdiscretXY - 1);

      double Zi            = target_pos.Z() - target_size.Z()/2.;
      double Zf            = target_pos.Z() + target_size.Z()/2.;
      double distanceStepZ = (Zf - Zi) / static_cast<double>(NstepsdiscretZ - 1);

      size_t border = 0;
      std::vector<std::vector<double> > PosXYZ{};
      SpaceDiscretization(Xi, Xf, Yi, Yf, NstepsdiscretXY, Zi, Zf, NstepsdiscretZ, border, PosXYZ);

      border = 1;

      for(size_t k = 0; k < nTimesDiscretization; ++k)
        {
          if(k != 0)
            {
              Xi            = temp_InteractionPointRecons.X() - distanceStepX;
              Xf            = temp_InteractionPointRecons.X() + distanceStepX;
              distanceStepX = (Xf - Xi) / static_cast<double>(Nstepsdiscretbox - 1);

              Yi            = temp_InteractionPointRecons.Y() - distanceStepY;
              Yf            = temp_InteractionPointRecons.Y() + distanceStepY;
              distanceStepY = (Yf - Yi) / static_cast<double>(Nstepsdiscretbox - 1);

              Zi            = temp_InteractionPointRecons.Z() - distanceStepZ;
              Zf            = temp_InteractionPointRecons.Z() + distanceStepZ;
              distanceStepZ = (Zf - Zi) / static_cast<double>(Nstepsdiscretbox - 1);

              SpaceDiscretization(Xi, Xf, Yi, Yf, Nstepsdiscretbox, Zi, Zf, Nstepsdiscretbox, border, PosXYZ);
            }

          for(size_t i = 0; i < PosXYZ.size(); ++i)
            {
              for(size_t j = 0; j < PrimaryTracks.size(); ++j)
                temp_f[j] = f_function(PrimaryTracks[j], PosXYZ[i]);

              temp_f[PrimaryTracks.size()] = f_function(BeamTracks[l], PosXYZ[i]);

              V_new = V_function(temp_f);

              if(V_new > V_temp)
                {
                  V_temp = V_new;
                  temp_f_values_IP = temp_f;
                  temp_InteractionPointRecons.SetXYZ(PosXYZ[i][0], PosXYZ[i][1], PosXYZ[i][2]);
                }
            }
        }

      if(V_temp > V_total)
        {
          V_total = V_temp;
          f_values_IP = temp_f_values_IP;
          InteractionPointRecons.SetXYZ(temp_InteractionPointRecons.X(), temp_InteractionPointRecons.Y(), temp_InteractionPointRecons.Z());
          i_BeamTracks_Vmax = l;
        }
    }

  return;
}


template <class Out>
void TPrimaryVertex<Out>::CovarianceMatrix(FiberTrackAna* BeamTrack, std::vector<FiberTrackAna*>& PrimaryTracks,
                                            TVector3& IP_average, std::vector<double>& f_values_IP,
                                                std::vector<std::vector<double> >& CovMatrix)
{
  std::vector<FiberTrackAna*> PrimaryTracks_Decisive;

  for(size_t i = 0; i < PrimaryTracks.size(); ++i)
    if(f_values_IP[i] > min_f_value)
      PrimaryTracks_Decisive.emplace_back(PrimaryTracks[i]);

  LocalHisto.h_nTrackCandidates->Fill(PrimaryTracks_Decisive.size(), "PrimaryTracks_Decisive", 1.);

  FiberTrackAna* const_BeamTrack = BeamTrack;
  std::vector<FiberTrackAna*> const_PrimaryTracks_Decisive = PrimaryTracks_Decisive;

  TVector3 temp_IP(-99., -99., -99.);
  std::vector<double> temp_f_values_IP(PrimaryTracks_Decisive.size() + 1, 0.);
  size_t temp_i_BeamTracks_Vmax = 99;

  std::vector<TVector3> variations_IP;
  double nDimensions = 6. * PrimaryTracks_Decisive.size() + 4.;

  TVector3 temp_IP_average;


  for(size_t variationsign = 0; variationsign <= 1; ++variationsign)
    {
      for(size_t idDir = 0; idDir <= 1; ++idDir)
        {
          for(size_t idFT = 0; idFT <= 2; ++idFT) // 0-> UFT3, 1-> MFT1, 2-> MFT2
            {
              for(size_t idPrimaryTrack_Decisive = 0; idPrimaryTrack_Decisive < PrimaryTracks_Decisive.size(); ++idPrimaryTrack_Decisive)
                {
                  std::vector<FiberHitXUV*> var_cont_xuv = PrimaryTracks_Decisive[idPrimaryTrack_Decisive]->GetContXUV();

                  FiberHitXUV* var_hit_xuv = new FiberHitXUV(var_cont_xuv[idFT]->GetPosX() +
                                                                                (1-idDir) * std::pow(-1., static_cast<double>(variationsign)) * std::sqrt(nDimensions) * sigma_FT,
                                                            var_cont_xuv[idFT]->GetPosY() + 
                                                                                (1-idDir) * std::pow(-1., static_cast<double>(variationsign)) * std::sqrt(nDimensions) * sigma_FT,
                                                            var_cont_xuv[idFT]->GetD(), var_cont_xuv[idFT]->GetHitX(), var_cont_xuv[idFT]->GetHitU(), var_cont_xuv[idFT]->GetHitV());

                  var_cont_xuv.erase(var_cont_xuv.begin() + idFT);
                  var_cont_xuv.emplace_back(var_hit_xuv);

                  FiberTrackAna* var_FiberTrackAna = new FiberTrackAna(var_cont_xuv, par.get());

                  PrimaryTracks_Decisive.erase(PrimaryTracks_Decisive.begin() + idPrimaryTrack_Decisive);
                  PrimaryTracks_Decisive.emplace_back(var_FiberTrackAna);

                  std::vector<FiberTrackAna*> BeamTracks = {BeamTrack};
                  InteractionPointFinder(BeamTracks, PrimaryTracks_Decisive, temp_IP, temp_f_values_IP, temp_i_BeamTracks_Vmax);
                  
                  variations_IP.emplace_back(temp_IP);
                  temp_IP_average += temp_IP;
                  PrimaryTracks_Decisive = const_PrimaryTracks_Decisive;
                }
            }

          for(size_t idFT = 0; idFT <=1; ++idFT) // 0-> UFT1, 1-> UFT2
            {
              std::vector<FiberHitXUV*> var_cont_xuv = BeamTrack->GetContXUV();

              FiberHitXUV* var_hit_xuv = new FiberHitXUV(var_cont_xuv[idFT]->GetPosX() +
                                                                            (1-idDir) * std::pow(-1., static_cast<double>(variationsign)) * std::sqrt(nDimensions) * sigma_FT,
                                                        var_cont_xuv[idFT]->GetPosY() + 
                                                                            (1-idDir) * std::pow(-1., static_cast<double>(variationsign)) * std::sqrt(nDimensions) * sigma_FT,
                                                        var_cont_xuv[idFT]->GetD(), var_cont_xuv[idFT]->GetHitX(), var_cont_xuv[idFT]->GetHitU(), var_cont_xuv[idFT]->GetHitV());

              var_cont_xuv.erase(var_cont_xuv.begin() + idFT);
              var_cont_xuv.emplace_back(var_hit_xuv);

              FiberTrackAna* var_FiberTrackAna = new FiberTrackAna(var_cont_xuv, par.get());

              std::vector<FiberTrackAna*> BeamTracks = {var_FiberTrackAna};
              InteractionPointFinder(BeamTracks, PrimaryTracks_Decisive, temp_IP, temp_f_values_IP, temp_i_BeamTracks_Vmax);
              
              variations_IP.emplace_back(temp_IP);
              temp_IP_average += temp_IP;
              BeamTrack = const_BeamTrack;
            }
        }
    }

  if(std::abs(variations_IP.size() - 2. * nDimensions) > 0.01)
      std::cout << "Error with dimensions covariance matrix\t" << variations_IP.size() << "\t" << 2. * nDimensions << "\n";

  IP_average = temp_IP_average * (1. / static_cast<double>(variations_IP.size()));

  std::vector<double> zero_vect(3, 0.);
  std::vector<std::vector<double> > zero_matrix(3, zero_vect);

  CovMatrix = zero_matrix;

  for(size_t i = 0; i < variations_IP.size(); ++i)
    {

      TVector3 temp_var_IP(variations_IP[i] - temp_IP_average);

      for(size_t j = 0; j < 3; ++j)
        for(size_t k = 0; k < 3; ++k)
          CovMatrix[j][k] += temp_var_IP(j) * temp_var_IP(k);
    }

  for(size_t j = 0; j < CovMatrix.size(); ++j)
    for(size_t k = 0; k < CovMatrix.size(); ++k)
      CovMatrix[j][k] /= static_cast<double>(variations_IP.size());

  return;
}



template class TPrimaryVertex<MCAnaEventG4Sol>;
template class TPrimaryVertex<Ana_WasaEvent>;
