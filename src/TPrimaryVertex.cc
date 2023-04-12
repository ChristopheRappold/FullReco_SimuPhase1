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
//#define COVARIANCE_MATRIX

using namespace std;
using namespace G4Sol;

template <class Out>
TPrimaryVertex<Out>::TPrimaryVertex(const THyphiAttributes& attribut)
    : TDataProcessInterface<Out>("PrimaryVertexReco"), att(attribut)
{
  att._logger->info("TPrimaryVertex::TPrimaryVertex");

  rand = new TRandom3();

  if(att.G4_simu == true)
    {
      //target_pos.SetXYZ(att.Target_PositionX, att.Target_PositionY, att.Target_PositionZ);
      //target_size.SetXYZ(att.Target_Size, att.Target_Size, att.Target_Size);
      //std::cout << "target_pos: " << att.Target_PositionX << "\t" << att.Target_PositionY << "\t" << att.Target_PositionZ << "\n";
      //std::cout << "target_size: " << att.Target_Size << "\n";
      target_pos.SetXYZ( 0., 0., 196.12); //in cm //Change to att.
      target_size.SetXYZ(3., 3.,   3.  ); //in cm //Change to att.

    }
  else
    {
      target_pos.SetXYZ( 0., 0., 196.12); //in cm //Change to att.
      target_size.SetXYZ(3., 3.,   3.  ); //in cm //Change to att.
    }

  RealPrimTrack = att.PV_RealPrimTrack;
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
      LocalHisto.h_PrimVtxstats->Fill("N_BeamTracks=0", 1.);
      return ReturnRes::Fine;
    }

  else if(result_full == -2)
    {
      att._logger->debug("No enough Primary tracks for primary vertex recons");
      LocalHisto.h_PrimVtxstats->Fill("N_PrimaryTracks=0", 1.);
      return ReturnRes::Fine;
    }

  else if(result_full == -3)
    {
      att._logger->debug("No proper PrimVtxTracks");
      LocalHisto.h_PrimVtxstats->Fill("V<Vmin", 1.);
      //return ReturnRes::PrimVtxError;
      return ReturnRes::Fine;
    }
  
  else if(result_full == -4)
    {
      att._logger->debug("More than 1 Beam tracks for primary vertex recons");
      LocalHisto.h_PrimVtxstats->Fill("N_BeamTracks>1", 1.);
      return ReturnRes::Fine;
    }
  LocalHisto.h_PrimVtxstats->Fill("Fine", 1.);

  return ReturnRes::Fine; 
}


template <class Out>
void TPrimaryVertex<Out>::SelectHists()
{

  LocalHisto.h_NFiberMult      = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_NFiberMult);
  LocalHisto.h_NCombsXUV_UFT12 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_NCombsXUV_UFT12);

  for(size_t i = 0; i < 5; ++i)
    LocalHisto.h_NCombsXUV[i]  = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_NCombsXUV[i]);

  for(size_t i = 2; i < 5; ++i)
    { 
      LocalHisto.h_CombsXUV_dvalue_theta[i]       = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_CombsXUV_dvalue_theta[i]);
      LocalHisto.h_CombsXUV_dvalue_phi[i]         = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_CombsXUV_dvalue_phi[i]);
      LocalHisto.h_CombsXUV_dvalue_phi_theta5[i]  = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_CombsXUV_dvalue_phi_theta5[i]);
      LocalHisto.h_CombsXUV_dvalue_phi_theta10[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_CombsXUV_dvalue_phi_theta10[i]);
    }

  LocalHisto.h_NHits_PrimaryTracks = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_NHits_PrimaryTracks);

  LocalHisto.h_nTrackCandidates   = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_nTrackCandidates);
  LocalHisto.h_DistanceBeamTracks = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DistanceBeamTracks);
  LocalHisto.h_PosZBeamTracks     = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_PosZBeamTracks);
  LocalHisto.h_thetaTracks        = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_thetaTracks);
  LocalHisto.h_chi2ndfTracks      = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_chi2ndfTracks);

  LocalHisto.h_fvalues           = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_fvalues);

  LocalHisto.h_InteractionPointPosX = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_InteractionPointPosX);
  LocalHisto.h_InteractionPointPosY = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_InteractionPointPosY);
  LocalHisto.h_InteractionPointPosZ = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_InteractionPointPosZ);

  LocalHisto.h_InteractionPointDistance_V_value  = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_InteractionPointDistance_V_value);
  
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

  IP_real.SetXYZ(RecoEvent.InteractionPoint[0], RecoEvent.InteractionPoint[1], RecoEvent.InteractionPoint[2]); //Simulation

  //Get beam tracks
  std::vector<PrimVtxTrack> BeamTracks_All = {};
  BeamTracksFinder(RecoEvent.ListHits, BeamTracks_All);

  LocalHisto.h_nTrackCandidates->Fill(BeamTracks_All.size(), "BeamTracks_All", 1.);
  for(size_t i = 0; i < BeamTracks_All.size(); ++i)
    {
      LocalHisto.h_thetaTracks->Fill(BeamTracks_All[i].GetTheta(), "BeamTracks_All", 1.);
      LocalHisto.h_chi2ndfTracks->Fill(BeamTracks_All[i].GetChi2NDF(), "BeamTracks_All", 1.);

#ifdef DEBUG_PRIMVTX
      std::cout << "BeamTracks_All[i].x: "       << BeamTracks_All[i].GetX()       << "\n"; 
      std::cout << "BeamTracks_All[i].y: "       << BeamTracks_All[i].GetY()       << "\n"; 
      std::cout << "BeamTracks_All[i].a: "       << BeamTracks_All[i].GetA()       << "\n"; 
      std::cout << "BeamTracks_All[i].b: "       << BeamTracks_All[i].GetB()       << "\n"; 
      std::cout << "BeamTracks_All[i].chi2ndf: " << BeamTracks_All[i].GetChi2NDF() << "\n\n"; 
#endif

    }

  std::vector<PrimVtxTrack> BeamTracks = {};
  BeamTracksSelector(BeamTracks_All, BeamTracks);
  //RecoEvent.BeamTracks = BeamTracks;

  LocalHisto.h_nTrackCandidates->Fill(BeamTracks.size(), "BeamTracks", 1.);

  if(BeamTracks.size() == 0)
    return -1;

  if(BeamTracks.size() > 1)
    return -4;

  double Beam_IPrecons_x = 0.;
  double Beam_IPrecons_y = 0.;

  for(size_t i = 0; i < BeamTracks.size(); ++i)
    {
      Beam_IPrecons_x += BeamTracks[i].GetX();
      Beam_IPrecons_y += BeamTracks[i].GetY();

      LocalHisto.h_thetaTracks->Fill(BeamTracks[i].GetTheta(), "BeamTracks", 1.);
      LocalHisto.h_chi2ndfTracks->Fill(BeamTracks[i].GetChi2NDF(), "BeamTracks", 1.);
    }

    Beam_IPrecons_x /= static_cast<double>(BeamTracks.size());
    Beam_IPrecons_y /= static_cast<double>(BeamTracks.size());

    IP_recons.SetX(Beam_IPrecons_x);
    IP_recons.SetY(Beam_IPrecons_y);

  //Get primary tracks

  std::vector<PrimVtxTrack> PrimaryTracks_All = {};
  if(att.PV_RealXUVComb == false)
    //PrimaryTracksFinder(RecoEvent.ListHits, PrimaryTracks_All);
    PrimaryTracksFinder_v2(RecoEvent.ListHits, BeamTracks[0], PrimaryTracks_All);
  else if(att.PV_RealXUVComb == true)
    RealPrimaryTracksFinder(RecoEvent.ListHits, RecoEvent.ListHitsToTracks, RecoEvent.TrackDAFSim, PrimaryTracks_All);

  LocalHisto.h_nTrackCandidates->Fill(PrimaryTracks_All.size(), "PrimaryTracks_All", 1.);
  for(size_t i = 0; i < PrimaryTracks_All.size(); ++i)
    {
      LocalHisto.h_thetaTracks->Fill(PrimaryTracks_All[i].GetTheta(), "PrimaryTracks_All", 1.);
      LocalHisto.h_chi2ndfTracks->Fill(PrimaryTracks_All[i].GetChi2NDF(), "PrimaryTracks_All", 1.);

      for(size_t j = 0; j < BeamTracks.size(); ++j)
        {
          double temp_CloseDistToBeam_PrimaryTracks;
          TVector3 temp_midpoint;
          CloseDist_TrackTrack(PrimaryTracks_All[i], BeamTracks[j], temp_CloseDistToBeam_PrimaryTracks, temp_midpoint);

          LocalHisto.h_DistanceBeamTracks->Fill(temp_CloseDistToBeam_PrimaryTracks, "PrimaryTracks_All", 1.);
          LocalHisto.h_PosZBeamTracks->Fill(temp_midpoint.Z(), "PrimaryTracks_All", 1.);
        }

#ifdef DEBUG_PRIMVTX
      std::cout << "PrimaryTracks_All[i].x: "       << PrimaryTracks_All[i].GetX()       << "\n"; 
      std::cout << "PrimaryTracks_All[i].y: "       << PrimaryTracks_All[i].GetY()       << "\n"; 
      std::cout << "PrimaryTracks_All[i].a: "       << PrimaryTracks_All[i].GetA()       << "\n"; 
      std::cout << "PrimaryTracks_All[i].b: "       << PrimaryTracks_All[i].GetB()       << "\n"; 
      std::cout << "PrimaryTracks_All[i].chi2ndf: " << PrimaryTracks_All[i].GetChi2NDF() << "\n\n"; 
#endif
    }

  std::vector<PrimVtxTrack> PrimaryTracks;
  PrimaryTracksSelector(PrimaryTracks_All, BeamTracks, PrimaryTracks);

  LocalHisto.h_nTrackCandidates->Fill(PrimaryTracks.size(), "PrimaryTracks", 1.);
  for(size_t i = 0; i < PrimaryTracks.size(); ++i)
    {
      LocalHisto.h_thetaTracks->Fill(PrimaryTracks[i].GetTheta(), "PrimaryTracks", 1.);
      LocalHisto.h_chi2ndfTracks->Fill(PrimaryTracks[i].GetChi2NDF(), "PrimaryTracks", 1.);

      for(size_t j = 0; j < BeamTracks.size(); ++j)
        {
          double temp_CloseDistToBeam_PrimaryTracks;
          TVector3 temp_midpoint;
          CloseDist_TrackTrack(PrimaryTracks[i], BeamTracks[j], temp_CloseDistToBeam_PrimaryTracks, temp_midpoint);

          LocalHisto.h_DistanceBeamTracks->Fill(temp_CloseDistToBeam_PrimaryTracks, "PrimaryTracks", 1.);
          LocalHisto.h_PosZBeamTracks->Fill(temp_midpoint.Z(), "PrimaryTracks", 1.);
        }

#ifdef DEBUG_PRIMVTX
      std::cout << "PrimaryTracks[i].x: "       << PrimaryTracks[i].GetX()       << "\n"; 
      std::cout << "PrimaryTracks[i].y: "       << PrimaryTracks[i].GetY()       << "\n"; 
      std::cout << "PrimaryTracks[i].a: "       << PrimaryTracks[i].GetA()       << "\n"; 
      std::cout << "PrimaryTracks[i].b: "       << PrimaryTracks[i].GetB()       << "\n"; 
      std::cout << "PrimaryTracks[i].chi2ndf: " << PrimaryTracks[i].GetChi2NDF() << "\n\n"; 
#endif
      
    }

  if(PrimaryTracks.size() == 0)
    return -2;


  std::vector<double> f_values_IP(PrimaryTracks.size() + 1, 0.);

  InteractionPointFinder(BeamTracks, PrimaryTracks, IP_recons, f_values_IP, i_BeamTracks_Vmax);

  if(i_BeamTracks_Vmax == 999)
    return -3;

  RecoEvent.PrimVtxRecons = IP_recons;

#ifdef VERTEX_RECONS_CHECK

  for(size_t i = 0; i < f_values_IP.size(); ++i)
    LocalHisto.h_fvalues->Fill(f_values_IP[i], 1.);

  LocalHisto.h_InteractionPointPosX->Fill(IP_recons.X(), 1.);
  LocalHisto.h_InteractionPointPosY->Fill(IP_recons.Y(), 1.);
  LocalHisto.h_InteractionPointPosZ->Fill(IP_recons.Z(), 1.);

  TVector3 IP_res(IP_real - IP_recons);

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
void TPrimaryVertex<Out>::CloseDist_TrackTrack(PrimVtxTrack& Track1, PrimVtxTrack& Track2,
                                      double& distance, TVector3& midpoint)
{
  std::vector<double> Track1Hit{Track1.GetX(), Track1.GetY(), target_pos.Z()};
  std::vector<double> Track1Par{Track1.GetA(), Track1.GetB(), 1.};

  std::vector<double> Track2Hit{Track2.GetX(), Track2.GetY(), target_pos.Z()};
  std::vector<double> Track2Par{Track2.GetA(), Track2.GetB(), 1.};

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
double TPrimaryVertex<Out>::CloseDist_TrackPoint(PrimVtxTrack& Track, TVector3& Point)
{
  TVector3 TrackHit(Track.GetX(), Track.GetY(), target_pos.Z());
  TVector3 TrackPar(Track.GetA(), Track.GetB(), 1.);

  TVector3 TrackHit_Point(Point - TrackHit);
  TVector3 TrackHit_Point_x_TrackPar(TrackHit_Point.Cross(TrackPar));

  double mag_TrackHit_Point_x_TrackPar = TrackHit_Point_x_TrackPar.Mag();
  double mag_TrackPar = TrackPar.Mag();

  return mag_TrackHit_Point_x_TrackPar / mag_TrackPar;
}


template <class Out>
double TPrimaryVertex<Out>::f_function(PrimVtxTrack& Track, std::vector<double>& PosXYZ)
{
  double slope_x     = Track.GetA();
  double intercept_x = Track.GetX() - slope_x * target_pos.Z();
  double slope_y     = Track.GetB();
  double intercept_y = Track.GetY() - slope_y * target_pos.Z();

  double sigma2 = std::pow(boxDistBeamXY/2., 2.) / 12.;

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
void TPrimaryVertex<Out>::BeamTracksFinder(std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > >& ListHits,
                                              std::vector<PrimVtxTrack>& BeamTracks_All)
{

  std::map<double,std::tuple<double,double,double> > empty_map;
  std::vector<std::map<double,std::tuple<double,double,double> > > vect_CombXUV_diffd = {empty_map, empty_map};

  if(att.G4_simu == true) //Simulation
    {
      PrimVtxTrack tmp_BeamTrack;
      tmp_BeamTrack.SetX(rand->Uniform(IP_real.X() - randIPXY, IP_real.X() + randIPXY));
      tmp_BeamTrack.SetY(rand->Uniform(IP_real.Y() - randIPXY, IP_real.Y() + randIPXY));
      tmp_BeamTrack.SetA(rand->Uniform(0. - randIPXY/(Zpos_Fiber[1] - Zpos_Fiber[0]), 0. + randIPXY/(Zpos_Fiber[1] - Zpos_Fiber[0])));
      tmp_BeamTrack.SetB(rand->Uniform(0. - randIPXY/(Zpos_Fiber[1] - Zpos_Fiber[0]), 0. + randIPXY/(Zpos_Fiber[1] - Zpos_Fiber[0])));
      tmp_BeamTrack.SetChi2NDF(0.);
      BeamTracks_All.emplace_back(tmp_BeamTrack);

      LocalHisto.h_NCombsXUV[0]->Fill(0., "All", 1.);
      LocalHisto.h_NCombsXUV[0]->Fill(0., "Selected", 1.);

      LocalHisto.h_NCombsXUV[1]->Fill(0., "All", 1.);
      LocalHisto.h_NCombsXUV[1]->Fill(0., "Selected", 1.);
    }
  else //Exp data
    {
      std::vector<size_t> FiberMult(6, 0);
      for(int i_det = 0; i_det < 2; ++i_det)
        {
          int i = id_detector[i_det];

          std::vector<genfit::AbsMeasurement*> hitx = Clusterize(ListHits[i]);
          std::vector<genfit::AbsMeasurement*> hitu = Clusterize(ListHits[i+1]);
          std::vector<genfit::AbsMeasurement*> hitv = Clusterize(ListHits[i+2]);

          FiberMult[i_det*3+0] = hitx.size();
          FiberMult[i_det*3+1] = hitu.size();
          FiberMult[i_det*3+2] = hitv.size();

          FindCombXUV_d(hitx, hitu, hitv, i_det, vect_CombXUV_diffd[i_det]);
          LocalHisto.h_NCombsXUV[i_det]->Fill(vect_CombXUV_diffd[i_det].size(), "All", 1.);

          if(SelectCombXUV_flag[i_det])
            {
              SelectCombXUV_d(vect_CombXUV_diffd[i_det]);
              LocalHisto.h_NCombsXUV[i_det]->Fill(vect_CombXUV_diffd[i_det].size(), "Selected", 1.);
            }
        }

      LocalHisto.h_NFiberMult->Fill(FiberMult[0], "UFT1_x", 1.);
      LocalHisto.h_NFiberMult->Fill(FiberMult[1], "UFT1_u", 1.);
      LocalHisto.h_NFiberMult->Fill(FiberMult[2], "UFT1_v", 1.);
      LocalHisto.h_NFiberMult->Fill(FiberMult[3], "UFT2_x", 1.);
      LocalHisto.h_NFiberMult->Fill(FiberMult[4], "UFT2_u", 1.);
      LocalHisto.h_NFiberMult->Fill(FiberMult[5], "UFT2_v", 1.);
      
      LocalHisto.h_NCombsXUV_UFT12->Fill(vect_CombXUV_diffd[0].size(), vect_CombXUV_diffd[1].size(), 1.);

      std::vector<std::tuple<double,double,double>> CombXUV_UFT1 = {};
      std::vector<std::tuple<double,double,double>> CombXUV_UFT2 = {};

      for(std::map<double,std::tuple<double,double,double>>::const_iterator it=vect_CombXUV_diffd[0].cbegin(); it!=vect_CombXUV_diffd[0].cend(); ++it)
        CombXUV_UFT1.emplace_back(it->second);


      for(std::map<double,std::tuple<double,double,double>>::const_iterator it=vect_CombXUV_diffd[1].cbegin(); it!=vect_CombXUV_diffd[1].cend(); ++it)
        CombXUV_UFT2.emplace_back(it->second);


      if(CombXUV_UFT1.size() == 0 || CombXUV_UFT2.size() == 0)
        return;

      BeamTracking(CombXUV_UFT1, CombXUV_UFT2, BeamTracks_All);
    }

  return;
}


template <class Out>
void TPrimaryVertex<Out>::BeamTracksSelector(std::vector<PrimVtxTrack>& BeamTracks_All,
                                              std::vector<PrimVtxTrack>& BeamTracks)
{

/*
  Possible cuts:
  - upper limit to chi2/ndf
  - upper limit to theta angle
  - upper limit to distance to target center == flag inside target
*/

  for(size_t i = 0; i < BeamTracks_All.size(); ++i)
    {
      if( ifCut_MaxChi2ndf_BeamTracks && (BeamTracks_All[i].GetChi2NDF() > MaxChi2ndf_BeamTracks) )
        continue;

      if( ifCut_MaxTheta_BeamTracks && (BeamTracks_All[i].GetTheta() > MaxTheta_BeamTracks) )
        continue;

      double CloseDistToTgtCenter_BeamTracks = CloseDist_TrackPoint(BeamTracks_All[i], target_pos);
      if( ifCut_CloseDistToTgtCenter_BeamTracks && (CloseDistToTgtCenter_BeamTracks > target_size.Mag()*MaxCloseDistToTgtCenter_BeamTracks) )
        continue;

      BeamTracks.emplace_back(BeamTracks_All[i]);
    }

  return;
}


template <class Out>
void TPrimaryVertex<Out>::PrimaryTracksFinder(std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > >& ListHits,
                                              std::vector<PrimVtxTrack>& PrimaryTracks_All)
{

  std::map<double,std::tuple<double,double,double> > empty_map;
  std::vector<std::map<double,std::tuple<double,double,double> > > vect_CombXUV_diffd = {empty_map, empty_map, empty_map, empty_map, empty_map};

  std::vector<size_t> FiberMult(9, 0);

  for(int i_det = 2; i_det < 5; ++i_det)
    {
      int i = id_detector[i_det];

      std::vector<genfit::AbsMeasurement*> hitx = Clusterize(ListHits[i]);
      std::vector<genfit::AbsMeasurement*> hitu = Clusterize(ListHits[i+1]);
      std::vector<genfit::AbsMeasurement*> hitv = Clusterize(ListHits[i+2]);

      FiberMult[(i_det-2)*3+0] = hitx.size();
      FiberMult[(i_det-2)*3+1] = hitu.size();
      FiberMult[(i_det-2)*3+2] = hitv.size();

      FindCombXUV_d(hitx, hitu, hitv, i_det, vect_CombXUV_diffd[i_det]);
      LocalHisto.h_NCombsXUV[i_det]->Fill(vect_CombXUV_diffd[i_det].size(), "All", 1.);

      if(SelectCombXUV_flag[i_det])
        {
          SelectCombXUV_d(vect_CombXUV_diffd[i_det]);
          LocalHisto.h_NCombsXUV[i_det]->Fill(vect_CombXUV_diffd[i_det].size(), "Selected", 1.);
        }

      for(std::map<double,std::tuple<double,double,double>>::const_iterator it=vect_CombXUV_diffd[i_det].cbegin(); it!=vect_CombXUV_diffd[i_det].cend(); ++it)
        Get_dvalueThetaPhi_CombXUV(i_det, it->second);
    }

  LocalHisto.h_NFiberMult->Fill(FiberMult[0], "UFT3_x", 1.);
  LocalHisto.h_NFiberMult->Fill(FiberMult[1], "UFT3_u", 1.);
  LocalHisto.h_NFiberMult->Fill(FiberMult[2], "UFT3_v", 1.);
  LocalHisto.h_NFiberMult->Fill(FiberMult[3], "MFT1_x", 1.);
  LocalHisto.h_NFiberMult->Fill(FiberMult[4], "MFT1_u", 1.);
  LocalHisto.h_NFiberMult->Fill(FiberMult[5], "MFT1_v", 1.);
  LocalHisto.h_NFiberMult->Fill(FiberMult[6], "MFT2_x", 1.);
  LocalHisto.h_NFiberMult->Fill(FiberMult[7], "MFT2_u", 1.);
  LocalHisto.h_NFiberMult->Fill(FiberMult[8], "MFT2_v", 1.);

  std::vector<std::tuple<double,double,double>> CombXUV_UFT3 = {};
  std::vector<std::tuple<double,double,double>> CombXUV_MFT1 = {};
  std::vector<std::tuple<double,double,double>> CombXUV_MFT2 = {};

  for(std::map<double,std::tuple<double,double,double>>::const_iterator it=vect_CombXUV_diffd[2].cbegin(); it!=vect_CombXUV_diffd[2].cend(); ++it)
      CombXUV_UFT3.emplace_back(it->second);

  for(std::map<double,std::tuple<double,double,double>>::const_iterator it=vect_CombXUV_diffd[3].cbegin(); it!=vect_CombXUV_diffd[3].cend(); ++it)
      CombXUV_MFT1.emplace_back(it->second);

  for(std::map<double,std::tuple<double,double,double>>::const_iterator it=vect_CombXUV_diffd[4].cbegin(); it!=vect_CombXUV_diffd[4].cend(); ++it)
      CombXUV_MFT2.emplace_back(it->second);

  if(CombXUV_UFT3.size() == 0 || CombXUV_MFT1.size() == 0 || CombXUV_MFT2.size() == 0)
    return;

  PrimaryTracking(CombXUV_UFT3, CombXUV_MFT1, CombXUV_MFT2, PrimaryTracks_All);

  return;
}


template <class Out>
void TPrimaryVertex<Out>::PrimaryTracksFinder_v2(std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > >& ListHits,
                                                  PrimVtxTrack& BeamTrack,
                                                  std::vector<PrimVtxTrack>& PrimaryTracks_All)
{

  std::vector<std::vector<std::vector<genfit::AbsMeasurement*>>> hits_FT = {{}, {}, {}};

  for(int i_det = 2; i_det < 5; ++i_det)
    {
      int i = id_detector[i_det];

      std::vector<genfit::AbsMeasurement*> hitx = Clusterize(ListHits[i]);
      std::vector<genfit::AbsMeasurement*> hitu = Clusterize(ListHits[i+1]);
      std::vector<genfit::AbsMeasurement*> hitv = Clusterize(ListHits[i+2]);

      hits_FT[(i_det-2)].emplace_back(hitx);
      hits_FT[(i_det-2)].emplace_back(hitu);
      hits_FT[(i_det-2)].emplace_back(hitv);
    }

  LocalHisto.h_NFiberMult->Fill(hits_FT[0][0].size(), "UFT3_x", 1.);
  LocalHisto.h_NFiberMult->Fill(hits_FT[0][1].size(), "UFT3_u", 1.);
  LocalHisto.h_NFiberMult->Fill(hits_FT[0][2].size(), "UFT3_v", 1.);
  LocalHisto.h_NFiberMult->Fill(hits_FT[1][0].size(), "MFT1_x", 1.);
  LocalHisto.h_NFiberMult->Fill(hits_FT[1][1].size(), "MFT1_u", 1.);
  LocalHisto.h_NFiberMult->Fill(hits_FT[1][2].size(), "MFT1_v", 1.);
  LocalHisto.h_NFiberMult->Fill(hits_FT[2][0].size(), "MFT2_x", 1.);
  LocalHisto.h_NFiberMult->Fill(hits_FT[2][1].size(), "MFT2_u", 1.);
  LocalHisto.h_NFiberMult->Fill(hits_FT[2][2].size(), "MFT2_v", 1.);


  std::tuple<double,double,double,double> BeamEllPar = BeamEllipseParam(BeamTrack); //a_ell, b_ell, h_ell, k_ell


  size_t n_cand = 0;
  size_t n_fitted = 0;
  size_t n_closebeam = 0;

  for(size_t d1=0; d1<2; ++d1)
  {
    if(onlyMFT)
      ++d1;
    for(size_t d2=d1+1; d2<3; ++d2)
      for(size_t d1lay1=0; d1lay1<2; ++d1lay1)
        for(size_t d1lay2=d1lay1+1; d1lay2<3; ++d1lay2)
          for(size_t d2lay1=0; d2lay1<2; ++d2lay1)
            for(size_t d2lay2=d2lay1+1; d2lay2<3; ++d2lay2)
              for(size_t hit1=0; hit1<hits_FT[d1][d1lay1].size(); ++hit1)
                for(size_t hit2=0; hit2<hits_FT[d1][d1lay2].size(); ++hit2)
                  for(size_t hit3=0; hit3<hits_FT[d2][d2lay1].size(); ++hit3)
                    for(size_t hit4=0; hit4<hits_FT[d2][d2lay2].size(); ++hit4)
                      {
                        PrimVtxTrack tmp_PrimaryTrack;
                        tmp_PrimaryTrack.SetFTHit((d1+2)*3+d1lay1, (hits_FT[d1][d1lay1][hit1]->getRawHitCoords())[0]);
                        tmp_PrimaryTrack.SetFTHit((d1+2)*3+d1lay2, (hits_FT[d1][d1lay2][hit2]->getRawHitCoords())[0]);
                        tmp_PrimaryTrack.SetFTHit((d2+2)*3+d2lay1, (hits_FT[d2][d2lay1][hit3]->getRawHitCoords())[0]);
                        tmp_PrimaryTrack.SetFTHit((d2+2)*3+d2lay2, (hits_FT[d2][d2lay2][hit4]->getRawHitCoords())[0]);

                        ++n_cand;

                        if(PrimVtxTracking(tmp_PrimaryTrack) == false)
                          continue;

                        ++n_fitted;

                        if(CloseToBeam(BeamEllPar, tmp_PrimaryTrack) == false)
                          continue;

                        ++n_closebeam;

                        PrimaryTracks_All.emplace_back(tmp_PrimaryTrack);

                      }
  }

  LocalHisto.h_nTrackCandidates->Fill(n_cand, "PrimaryTracks_n_cand", 1.);
  LocalHisto.h_nTrackCandidates->Fill(n_fitted, "PrimaryTracks_n_fitted", 1.);
  LocalHisto.h_nTrackCandidates->Fill(n_closebeam, "PrimaryTracks_n_closebeam", 1.);

  return;
}


template <class Out>
void TPrimaryVertex<Out>::RealPrimaryTracksFinder(std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > >& ListHits,
                                                    std::vector<std::vector<int> >& ListHitsToTracks,
                                                    std::unordered_map<int, std::vector<std::vector<SimHit> > >& TrackDAFSim,
                                                      std::vector<PrimVtxTrack>& PrimaryTracks_All)
{
  if(RealPrimTrack)
    {
      for(std::pair<int, std::vector<std::vector<SimHit> > > it_track : TrackDAFSim)
        {
          int i_closest_det = 0;
          size_t nHits = 0;

          for(int i_det = 4; i_det >= 2; --i_det)
            {
              int i = id_detector[i_det];

              if(!it_track.second[i+2].empty())
                {
                  ++nHits;
                  i_closest_det = i+2;
                }

              if(!it_track.second[i+1].empty())
                {
                  ++nHits;
                  i_closest_det = i+1;
                }

              if(!it_track.second[i].empty())
                {
                  ++nHits;
                  i_closest_det = i;
                }
            }

          LocalHisto.h_NHits_PrimaryTracks->Fill(nHits, "RealPrimTrack", 1.);

          if(nHits < 8)
            continue;

          PrimVtxTrack tmp_RealPrimTrack;

          tmp_RealPrimTrack.SetA(it_track.second[i_closest_det][0].momX / it_track.second[i_closest_det][0].momZ);
          tmp_RealPrimTrack.SetB(it_track.second[i_closest_det][0].momY / it_track.second[i_closest_det][0].momZ);

          tmp_RealPrimTrack.SetX(it_track.second[i_closest_det][0].hitX - tmp_RealPrimTrack.GetA() * (it_track.second[i_closest_det][0].hitZ - target_pos.Z()));
          tmp_RealPrimTrack.SetY(it_track.second[i_closest_det][0].hitY - tmp_RealPrimTrack.GetB() * (it_track.second[i_closest_det][0].hitZ - target_pos.Z()));

          tmp_RealPrimTrack.SetChi2NDF(0.);

          PrimaryTracks_All.emplace_back(tmp_RealPrimTrack);
        }

      return;
    }

  std::vector<std::tuple<double,double,double>> CombXUV_UFT3 = {};
  std::vector<std::tuple<double,double,double>> CombXUV_MFT1 = {};
  std::vector<std::tuple<double,double,double>> CombXUV_MFT2 = {};

  for(int i_det = 2; i_det < 5; ++i_det)
    {
      int i = id_detector[i_det];

      for(std::pair<int, std::vector<std::vector<SimHit> > > it_track : TrackDAFSim)
        {
          const int id_track = it_track.first;

          double tmp_hit_x = 0.;
          double tmp_hit_u = 0.;
          double tmp_hit_v = 0.;

          int counter_x = 0;
          int counter_u = 0;
          int counter_v = 0;

          std::vector<size_t> tmp_idhit_x = {};
          std::vector<size_t> tmp_idhit_u = {};
          std::vector<size_t> tmp_idhit_v = {};

          for(size_t j = 0; j < ListHits[i].size(); ++j)
            {
              if(ListHitsToTracks[i][j] == id_track)
                {
                  TVectorD& vpos_x = ListHits[i][j]->getRawHitCoords();
                  tmp_hit_x += vpos_x[0];
                  ++counter_x;
                  tmp_idhit_x.emplace_back(j);
                }
            }

          for(size_t j = 0; j < ListHits[i+1].size(); ++j)
            {
              if(ListHitsToTracks[i+1][j] == id_track)
                {
                  TVectorD& vpos_u = ListHits[i+1][j]->getRawHitCoords();
                  tmp_hit_u += vpos_u[0];
                  ++counter_u;
                  tmp_idhit_u.emplace_back(j);
                }
            }

          for(size_t j = 0; j < ListHits[i+2].size(); ++j)
            {
              if(ListHitsToTracks[i+2][j] == id_track)
                {
                  TVectorD& vpos_v = ListHits[i+2][j]->getRawHitCoords();
                  tmp_hit_v += vpos_v[0];
                  ++counter_v;
                  tmp_idhit_v.emplace_back(j);
                }
            }

          if((counter_x == 0) || (counter_u == 0) || (counter_v == 0))
            continue;

          tmp_hit_x /= static_cast<double>(counter_x);
          tmp_hit_u /= static_cast<double>(counter_u);
          tmp_hit_v /= static_cast<double>(counter_v);

          if(i_det == 2)
            CombXUV_UFT3.emplace_back(std::make_tuple(tmp_hit_x, tmp_hit_u, tmp_hit_v));
          else if(i_det == 3)
            CombXUV_MFT1.emplace_back(std::make_tuple(tmp_hit_x, tmp_hit_u, tmp_hit_v));
          else if(i_det == 4)
            CombXUV_MFT2.emplace_back(std::make_tuple(tmp_hit_x, tmp_hit_u, tmp_hit_v));

          Get_dvalueThetaPhi_CombXUV(i_det, std::make_tuple(tmp_hit_x, tmp_hit_u, tmp_hit_v));
        }

      if(i_det == 2)
        LocalHisto.h_NCombsXUV[i_det]->Fill(CombXUV_UFT3.size(), "All", 1.);
      else if(i_det == 3)
        LocalHisto.h_NCombsXUV[i_det]->Fill(CombXUV_MFT1.size(), "All", 1.);
      else if(i_det == 4)
        LocalHisto.h_NCombsXUV[i_det]->Fill(CombXUV_MFT2.size(), "All", 1.);
    }

  if(CombXUV_UFT3.size() == 0 || CombXUV_MFT1.size() == 0 || CombXUV_MFT2.size() == 0)
    return;

  PrimaryTracking(CombXUV_UFT3, CombXUV_MFT1, CombXUV_MFT2, PrimaryTracks_All);

  return;
}


template <class Out>
void TPrimaryVertex<Out>::PrimaryTracksSelector(std::vector<PrimVtxTrack>& PrimaryTracks_All, std::vector<PrimVtxTrack>& BeamTracks,
                                                    std::vector<PrimVtxTrack>& PrimaryTracks)
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
      if( ifCut_MaxChi2ndf_PrimaryTracks && (PrimaryTracks_All[i].GetChi2NDF() > MaxChi2ndf_PrimaryTracks) )
        continue;

      if( ifCut_MinTheta_PrimaryTracks && (PrimaryTracks_All[i].GetTheta() < MinTheta_PrimaryTracks) )
        continue;

      if( ifCut_MaxTheta_PrimaryTracks && (PrimaryTracks_All[i].GetTheta() > MaxTheta_PrimaryTracks) )
        continue;

      double CloseDistToBeam_PrimaryTracks = 9999.;
      double DistZBeamToTarget_PrimaryTracks = 9999.;
      for(size_t j = 0; j < BeamTracks.size(); ++j)
        {
          double temp_CloseDistToBeam_PrimaryTracks;
          TVector3 temp_midpoint;
          CloseDist_TrackTrack(PrimaryTracks_All[i], BeamTracks[j], temp_CloseDistToBeam_PrimaryTracks, temp_midpoint);

          if(temp_CloseDistToBeam_PrimaryTracks < CloseDistToBeam_PrimaryTracks)
            CloseDistToBeam_PrimaryTracks = temp_CloseDistToBeam_PrimaryTracks;
          if(std::fabs(target_pos.Z()-temp_midpoint.Z()) < DistZBeamToTarget_PrimaryTracks)
            DistZBeamToTarget_PrimaryTracks = temp_CloseDistToBeam_PrimaryTracks;        
        }
      if( ifCut_MaxCloseDistToBeam_PrimaryTracks && (CloseDistToBeam_PrimaryTracks > MaxCloseDistToBeam_PrimaryTracks) )
        continue;

      if( ifCut_MaxDistZBeamToTarget_PrimaryTracks && (DistZBeamToTarget_PrimaryTracks > MaxDistZBeamToTarget_PrimaryTracks * target_size.Z()))
        continue;

      double CloseDistToTgtCenter_PrimaryTracks = CloseDist_TrackPoint(PrimaryTracks_All[i], target_pos);
      if( ifCut_CloseDistToTgtCenter_PrimaryTracks && (CloseDistToTgtCenter_PrimaryTracks > target_size.Mag()*MaxCloseDistToTgtCenter_PrimaryTracks) )
        continue;
/*
      if (ifCut_FlagPSBhit_PrimaryTracks && (FlagPSBhit_PrimaryTracks == PrimaryTracks_All[i]->IsFlagPSB()))
        continue;

      if (ifCut_FlagPSFEhit_PrimaryTracks && (FlagPSFEhit_PrimaryTracks == PrimaryTracks_All[i]->IsFlagPSFE()))
        continue;

      if (ifCut_FlagPSBEhit_PrimaryTracks && (FlagPSBEhit_PrimaryTracks == PrimaryTracks_All[i]->IsFlagPSBE()))
        continue;

      if (ifCut_FlagPShit_PrimaryTracks && (FlagPShit_PrimaryTracks ==
                                             (PrimaryTracks_All[i]->IsFlagPSB() || PrimaryTracks_All[i]->IsFlagPSFE() || PrimaryTracks_All[i]->IsFlagPSBE() ) ) )
        continue;
*/
      PrimaryTracks.emplace_back(PrimaryTracks_All[i]);
    }

  return;
}


template <class Out>
void TPrimaryVertex<Out>::InteractionPointFinder(std::vector<PrimVtxTrack>& BeamTracks, std::vector<PrimVtxTrack>& PrimaryTracks,
                                                  TVector3& InteractionPointRecons, std::vector<double>& f_values_IP, size_t& i_BeamTracks_Vmax)
{
  double V_total = 0.;
  std::vector<double> temp_f_values_IP(PrimaryTracks.size() + 1, 0.);

  for(size_t l = 0; l < BeamTracks.size(); ++l)
    {
      TVector3 temp_InteractionPointRecons(-99., -99., -99.);

      std::vector<double> temp_f(PrimaryTracks.size() + 1, 0.);

      double V_temp = 0.;

      double Xi            = BeamTracks[l].GetX() - boxDistBeamXY;
      double Xf            = BeamTracks[l].GetX() + boxDistBeamXY;
      double distanceStepX = (Xf - Xi) / static_cast<double>(NstepsdiscretXY - 1);

      double Yi            = BeamTracks[l].GetY() - boxDistBeamXY;
      double Yf            = BeamTracks[l].GetY() + boxDistBeamXY;
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

              PosXYZ.clear();
              SpaceDiscretization(Xi, Xf, Yi, Yf, Nstepsdiscretbox, Zi, Zf, Nstepsdiscretbox, border, PosXYZ);
            }

          for(size_t i = 0; i < PosXYZ.size(); ++i)
            {
              for(size_t j = 0; j < PrimaryTracks.size(); ++j)
                temp_f[j] = f_function(PrimaryTracks[j], PosXYZ[i]);

              temp_f[PrimaryTracks.size()] = f_function(BeamTracks[l], PosXYZ[i]);

              double V_new = V_function(temp_f);

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

  if(ifCut_min_V_value && (V_total < min_V_value))
    {
      f_values_IP.clear();
      InteractionPointRecons.SetXYZ(-999., -999., -999.);
      i_BeamTracks_Vmax = 999;

      return;
    }

  if(InteractionPointRecons.X() < target_pos.X() - target_size.X()/2.)
    InteractionPointRecons.SetX(target_pos.X() - target_size.X()/2.);
  else if(InteractionPointRecons.X() > target_pos.X() + target_size.X()/2.)
    InteractionPointRecons.SetX(target_pos.X() + target_size.X()/2.);

  if(InteractionPointRecons.Y() < target_pos.Y() - target_size.Y()/2.)
    InteractionPointRecons.SetY(target_pos.Y() - target_size.Y()/2.);
  else if(InteractionPointRecons.Y() > target_pos.Y() + target_size.Y()/2.)
    InteractionPointRecons.SetY(target_pos.Y() + target_size.Y()/2.);

  if(InteractionPointRecons.Z() < target_pos.Z() - target_size.Z()/2.)
    InteractionPointRecons.SetZ(target_pos.Z() - target_size.Z()/2.);
  else if(InteractionPointRecons.Z() > target_pos.Z() + target_size.Z()/2.)
    InteractionPointRecons.SetZ(target_pos.Z() + target_size.Z()/2.);

  TVector3 IP_res(IP_real - InteractionPointRecons);
  LocalHisto.h_InteractionPointDistance_V_value->Fill(IP_res.Mag(), V_total, 1.);

  return;
}


/*
template <class Out>
void TPrimaryVertex<Out>::CovarianceMatrix(PrimVtxTrack& BeamTrack, std::vector<PrimVtxTrack>& PrimaryTracks,
                                            TVector3& IP_average, std::vector<double>& f_values_IP,
                                                std::vector<std::vector<double> >& CovMatrix)
{
  std::vector<PrimVtxTrack> PrimaryTracks_Decisive;

  for(size_t i = 0; i < PrimaryTracks.size(); ++i)
    if(f_values_IP[i] > min_f_value)
      PrimaryTracks_Decisive.emplace_back(PrimaryTracks[i]);

  LocalHisto.h_nTrackCandidates->Fill(PrimaryTracks_Decisive.size(), "PrimaryTracks_Decisive", 1.);

  PrimVtxTrack const_BeamTrack = BeamTrack;
  std::vector<PrimVtxTrack> const_PrimaryTracks_Decisive = PrimaryTracks_Decisive;

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
*/

/*
template<class Out>
void TPrimaryVertex<Out>::FindHitXUV_v3(const std::vector<genfit::AbsMeasurement*>& hitx,
                                          const std::vector<genfit::AbsMeasurement*>& hitu,
                                          const std::vector<genfit::AbsMeasurement*>& hitv,
                                            int id_det)
{
  std::set<int> flag_x;
  std::set<int> flag_u;
  std::set<int> flag_v;

  int nx = hitx.size();
  int nu = hitu.size();
  int nv = hitv.size();
  if((nx==0)||(nu==0)||(nv==0))
    return;

  double ang_u = ang[id_det][id_mid[id_det]] * M_PI / 180.;

  while(true)
    {
      double hit_xu_x = -999.;
      double hit_xu_y =  999.;
      double hit_xv_y = -999.;

      double diff_d = 999.;

      double hit_u = -999.;
      double hit_v = -999.;

      int used_i = -999;
      int used_j = -999;
      int used_k = -999;

      for(int i = 0; i < nx ; ++i)
        {
          if(flag_x.size()>0 && flag_x.find(i)!=flag_x.end())
            continue;

          TVectorD& vpos_x = hitx[i]->getRawHitCoords();
          double pos_x = vpos_x[0];

          for(int j = 0; j < nu ; ++j)
            {
              if(flag_u.size()>0 && flag_u.find(j)!=flag_u.end())
                continue;
              
              TVectorD& vpos_u = hitu[j]->getRawHitCoords();
              double pos_u = vpos_u[0];
              double tmp_hit_u = pos_u;
              pos_u = pos_u / std::cos(ang_u);

              double tmp_hit_xu_y = 999.;
              if(ang_u>0.)
                tmp_hit_xu_y = -std::tan((M_PI/2. - ang_u))*pos_x + std::tan((M_PI/2. -ang_u))*pos_u;
              else
                tmp_hit_xu_y =  std::tan((M_PI/2. + ang_u))*pos_x - std::tan((M_PI/2. +ang_u))*pos_u;

              for(int k = 0; k < nv; ++k)
                {
                  if(flag_v.size()>0 && flag_v.find(k)!=flag_v.end())
                    continue;

                  TVectorD& vpos_v = hitv[k]->getRawHitCoords();
                  double pos_v = vpos_v[0];
                  double tmp_hit_v = pos_v;
                  pos_v = pos_v / std::cos(ang_u);

                  double tmp_hit_xv_y = -999.;
                  if(ang_u>0)
                    tmp_hit_xv_y =  std::tan((M_PI/2. - ang_u))*pos_x - std::tan((M_PI/2. - ang_u))*pos_v;
                  else
                    tmp_hit_xv_y = -std::tan((M_PI/2. + ang_u))*pos_x + std::tan((M_PI/2. + ang_u))*pos_v;

                  double tmp_x = pos_x + (tmp_hit_xu_y-tmp_hit_xv_y)/2.*std::tan(ang_u/2.);
                  double tmp_y = (tmp_hit_xu_y+tmp_hit_xv_y)/2.;

                  double tmp_d = tmp_hit_xu_y-tmp_hit_xv_y;
                  double expected_d = d_function1(id_det, tmp_x, tmp_y);
                  //double expected_d = d_function2(id_det, pos_x);
                  double tmp_diff_d = std::fabs(tmp_d - expected_d);

                  if(tmp_diff_d < diff_d)
                    {
                      diff_d = tmp_diff_d;
                      hit_xu_x = pos_x; 
                      hit_xu_y = tmp_hit_xu_y;
                      hit_xv_y = tmp_hit_xv_y;
                      hit_u = tmp_hit_u;
                      hit_v = tmp_hit_v;
                      used_i = i;
                      used_j = j;
                      used_k = k;
                    }
                }
            }
        }

        if( (std::fabs(hit_xu_x+999.) < 1.) || (std::fabs(hit_xu_y-999.) < 1.) || (std::fabs(hit_xv_y+999.) < 1.) )
          break;

        if(diff_d > cut_diff_d[id_det])
          break;

        LocalHisto.h_FiberHit_dvalue[id_det]->Fill(hit_xu_y-hit_xv_y, 1.);

        double buf_x = hit_xu_x + (hit_xu_y-hit_xv_y)/2.*std::tan(ang_u/2.);
        double buf_y = (hit_xu_y+hit_xv_y)/2.;

        TVector2 tmp_hitxy(buf_x, buf_y);
        vect_HitXY[id_det].emplace_back(tmp_hitxy);

        vect_CombHit[id_det].emplace_back(std::make_tuple(hit_xu_x, hit_u, hit_v));

        flag_x.insert(used_i);
        flag_u.insert(used_j);
        flag_v.insert(used_k);  
    }

  return;
}


template<class Out>
void TPrimaryVertex<Out>::FindHitXUV_v4(const std::vector<genfit::AbsMeasurement*>& hitx,
                                          const std::vector<genfit::AbsMeasurement*>& hitu,
                                          const std::vector<genfit::AbsMeasurement*>& hitv,
                                            int id_det)
{
  int nx = hitx.size();
  int nu = hitu.size();
  int nv = hitv.size();
  if((nx==0)||(nu==0)||(nv==0))
    return;

  double ang_u = ang[id_det][id_mid[id_det]] * M_PI / 180.;

  for(int i = 0; i < nx ; ++i)
    {
      TVectorD& vpos_x = hitx[i]->getRawHitCoords();
      double pos_x = vpos_x[0];

      for(int j = 0; j < nu ; ++j)
        {
          TVectorD& vpos_u = hitu[j]->getRawHitCoords();
          double pos_u = vpos_u[0];
          pos_u = pos_u / std::cos(ang_u);

          double tmp_hit_xu_y = 999.;
          if(ang_u>0.)
            tmp_hit_xu_y = -std::tan((M_PI/2. - ang_u))*pos_x + std::tan((M_PI/2. -ang_u))*pos_u;
          else
            tmp_hit_xu_y =  std::tan((M_PI/2. + ang_u))*pos_x - std::tan((M_PI/2. +ang_u))*pos_u;

          for(int k = 0; k < nv; ++k)
            {
              TVectorD& vpos_v = hitv[k]->getRawHitCoords();
              double pos_v = vpos_v[0];
              pos_v = pos_v / std::cos(ang_u);

              double tmp_hit_xv_y = -999.;
              if(ang_u>0)
                tmp_hit_xv_y =  std::tan((M_PI/2. - ang_u))*pos_x - std::tan((M_PI/2. - ang_u))*pos_v;
              else
                tmp_hit_xv_y = -std::tan((M_PI/2. + ang_u))*pos_x + std::tan((M_PI/2. + ang_u))*pos_v;

              double tmp_x = pos_x + (tmp_hit_xu_y-tmp_hit_xv_y)/2.*std::tan(ang_u/2.);
              double tmp_y = (tmp_hit_xu_y+tmp_hit_xv_y)/2.;

              double tmp_d = tmp_hit_xu_y-tmp_hit_xv_y;
              double expected_d = d_function1(id_det, tmp_x, tmp_y);
              //double expected_d = d_function2(id_det, pos_x);
              double tmp_diff_d = std::fabs(tmp_d - expected_d);

              if(tmp_diff_d < cut_diff_d[id_det])
                {
                  LocalHisto.h_FiberHit_dvalue[id_det]->Fill(tmp_hit_xu_y-tmp_hit_xv_y, 1.);

                  TVector2 tmp_hitxy(tmp_x, tmp_y);
                  vect_HitXY[id_det].emplace_back(tmp_hitxy);

                  vect_CombHit[id_det].emplace_back(std::make_tuple(vpos_x[0], vpos_u[0], vpos_v[0]));
                }
            }
        }
    }

  return;
}
*/

template<class Out>
void TPrimaryVertex<Out>::FindCombXUV_d(const std::vector<genfit::AbsMeasurement*>& hitx, const std::vector<genfit::AbsMeasurement*>& hitu,
                                          const std::vector<genfit::AbsMeasurement*>& hitv, int id_det,
                                            std::map<double,std::tuple<double,double,double> >& CombXUV_diffd)
{
  int nx = hitx.size();
  int nu = hitu.size();
  int nv = hitv.size();
  if((nx==0)||(nu==0)||(nv==0))
    return;

  double ang_u = ang[id_det][id_mid[id_det]] * M_PI / 180.;

  for(int i = 0; i < nx ; ++i)
    {
      TVectorD& vpos_x = hitx[i]->getRawHitCoords();
      double pos_x = vpos_x[0];

      for(int j = 0; j < nu ; ++j)
        {
          TVectorD& vpos_u = hitu[j]->getRawHitCoords();
          double pos_u = vpos_u[0] / std::cos(ang_u);

          double tmp_hit_xu_y = 999.;
          if(ang_u>0.)
            tmp_hit_xu_y = -std::tan((M_PI/2. - ang_u))*pos_x + std::tan((M_PI/2. -ang_u))*pos_u;
          else
            tmp_hit_xu_y =  std::tan((M_PI/2. + ang_u))*pos_x - std::tan((M_PI/2. +ang_u))*pos_u;

          for(int k = 0; k < nv; ++k)
            {
              TVectorD& vpos_v = hitv[k]->getRawHitCoords();
              double pos_v = vpos_v[0] / std::cos(ang_u);

              double tmp_hit_xv_y = -999.;
              if(ang_u>0)
                tmp_hit_xv_y =  std::tan((M_PI/2. - ang_u))*pos_x - std::tan((M_PI/2. - ang_u))*pos_v;
              else
                tmp_hit_xv_y = -std::tan((M_PI/2. + ang_u))*pos_x + std::tan((M_PI/2. + ang_u))*pos_v;

              double tmp_x = pos_x + (tmp_hit_xu_y-tmp_hit_xv_y)/2.*std::tan(ang_u/2.);
              double tmp_y = (tmp_hit_xu_y+tmp_hit_xv_y)/2.;

              double tmp_d = tmp_hit_xu_y-tmp_hit_xv_y;
              double expected_d = d_function1(id_det, tmp_x, tmp_y);
              //double expected_d = d_function2(id_det, pos_x);
              double tmp_diff_d = std::fabs(tmp_d - expected_d);

              if(tmp_diff_d < cut_diff_d[id_det])
                {
                  //std::cout << "tmp_x: " << tmp_x << "\t tmp_y: " << tmp_y << "\n";
                  CombXUV_diffd.insert(std::pair<double,std::tuple<double,double,double>>(tmp_diff_d, std::make_tuple(vpos_x[0], vpos_u[0], vpos_v[0])));
/*
                  if(id_det != 0 && id_det != 1)
                    {
                      double x = tmp_x - IP_recons.X();
                      double y = tmp_y - IP_recons.Y();
                      double z = Zpos_Fiber[id_det] - target_pos.Z();

                      double theta = TMath::ATan2( std::sqrt(x * x + y * y), z) * 180. / M_PI;
                      LocalHisto.h_CombsXUV_dvalue_theta[id_det]->Fill(theta, tmp_d, 1.);

                      double phi = TMath::ATan2( y , x ) * 180. / M_PI;
                      LocalHisto.h_CombsXUV_dvalue_phi[id_det]->Fill(phi, tmp_d, 1.);
                    }
*/  
                }
            }
        }
    }

  return;
}


template<class Out>
void TPrimaryVertex<Out>::SelectCombXUV_d(std::map<double,std::tuple<double,double,double> >& CombXUV_diffd)
{
  if(CombXUV_diffd.size() == 0 || CombXUV_diffd.size() == 1)
    return;

  std::set<double> used_x;
  std::set<double> used_u;
  std::set<double> used_v;

  for(std::map<double,std::tuple<double,double,double>>::const_iterator it=CombXUV_diffd.cbegin(); it!=CombXUV_diffd.cend();)
    {
      if(used_x.size()>0 && used_x.find(get<0>(it->second))!=used_x.end())
          it = CombXUV_diffd.erase(it);
      else if(used_u.size()>0 && used_u.find(get<1>(it->second))!=used_u.end())
          it = CombXUV_diffd.erase(it);
      else if(used_v.size()>0 && used_v.find(get<2>(it->second))!=used_v.end())
          it = CombXUV_diffd.erase(it);
      else
        {
          used_x.insert(get<0>(it->second));
          used_u.insert(get<1>(it->second));
          used_v.insert(get<2>(it->second));
          ++it;
        }
    }

  return;
}



template<class Out>
std::vector<genfit::AbsMeasurement*> TPrimaryVertex<Out>::Clusterize(const std::vector<std::unique_ptr<genfit::AbsMeasurement>>& hit)
{
  std::vector<genfit::AbsMeasurement*> hit_cluster = {};
  for(size_t i = 0; i < hit.size(); ++i)
    hit_cluster.emplace_back(hit[i]->clone());

  if(hit_cluster.size() < 2)
    return hit_cluster;

  std::set<int> flag_used;

  for(size_t i = 0; i < hit_cluster.size(); ++i)
    {
      if(flag_used.size()>0 && flag_used.find(i)!=flag_used.end())
        continue;

      TVectorD& vpos = hit_cluster[i]->getRawHitCoords();
      double pos_cluster = vpos[0];
      double edge_left  = pos_cluster;
      double edge_right = pos_cluster;
      size_t size_cluster = 1;

      for(size_t j = i+1; j < hit_cluster.size(); ++j)
        {
          if(flag_used.size()>0 && flag_used.find(j)!=flag_used.end())
            continue;

          TVectorD& tmp_vpos = hit_cluster[j]->getRawHitCoords();

          if(((pos_cluster-tmp_vpos[0])>0.) && std::fabs(edge_left-tmp_vpos[0])<0.0826)
            {
              pos_cluster = (pos_cluster * size_cluster + tmp_vpos[0]) / static_cast<double>(size_cluster + 1);
              flag_used.insert(j);
              edge_left = tmp_vpos[0];
              ++size_cluster;
              j=i;
            }
          else if(((tmp_vpos[0]-pos_cluster)>0.) && std::fabs(tmp_vpos[0]-edge_right)<0.0826)
            {
              pos_cluster = (pos_cluster * size_cluster + tmp_vpos[0]) / static_cast<double>(size_cluster + 1);
              flag_used.insert(j);
              edge_right = tmp_vpos[0];
              ++size_cluster;
              j=i;
            }
        }

      vpos[0] = pos_cluster;
      hit_cluster[i]->setRawHitCoords(vpos);
    }

  for(int i = (int)hit_cluster.size()-1; i >= 0; --i)
    {
      if(flag_used.size()>0 && flag_used.find(i)!=flag_used.end())
        {
          delete hit_cluster[i];
          hit_cluster.erase(hit_cluster.begin() + i);
        }
    }

  return hit_cluster;
}


template<class Out>
std::vector<genfit::AbsMeasurement*> TPrimaryVertex<Out>::Clusterize(const std::vector<genfit::AbsMeasurement*>& hit)
{
  std::vector<genfit::AbsMeasurement*> hit_cluster = {};
  for(size_t i = 0; i < hit.size(); ++i)
    hit_cluster.emplace_back(hit[i]->clone());

  if(hit_cluster.size() < 2)
    return hit_cluster;

  std::set<int> flag_used;

  for(size_t i = 0; i < hit_cluster.size(); ++i)
    {
      if(flag_used.size()>0 && flag_used.find(i)!=flag_used.end())
        continue;

      TVectorD& vpos = hit_cluster[i]->getRawHitCoords();
      double pos_cluster = vpos[0];
      double edge_left  = pos_cluster;
      double edge_right = pos_cluster;
      size_t size_cluster = 1;

      for(size_t j = i+1; j < hit_cluster.size(); ++j)
        {
          if(flag_used.size()>0 && flag_used.find(j)!=flag_used.end())
            continue;

          TVectorD& tmp_vpos = hit_cluster[j]->getRawHitCoords();

          if(((pos_cluster-tmp_vpos[0])>0.) && std::fabs(edge_left - tmp_vpos[0])<0.0826)
            {
              double tmp_pos_cluster = (pos_cluster * size_cluster + tmp_vpos[0]) / (size_cluster + 1);
              pos_cluster = tmp_pos_cluster;

              flag_used.insert(j);
              edge_left = tmp_vpos[0];
              ++size_cluster;
              j=i;
            }
          else if(((tmp_vpos[0]-pos_cluster)>0.) && std::fabs(tmp_vpos[0]-edge_right)<0.0826)
            {
              double tmp_pos_cluster = (pos_cluster * size_cluster + tmp_vpos[0]) / (size_cluster + 1);
              pos_cluster = tmp_pos_cluster;

              flag_used.insert(j);
              edge_right = tmp_vpos[0];
              ++size_cluster;
              j=i;
            }
        }

      vpos[0] = pos_cluster;
      hit_cluster[i]->setRawHitCoords(vpos);
    }

  for(int i = (int)hit_cluster.size()-1; i >= 0; --i)
    {
      if(flag_used.size()>0 && flag_used.find(i)!=flag_used.end())
        {
          delete hit_cluster[i];
          hit_cluster.erase(hit_cluster.begin() + i);
        }
    }

  return hit_cluster;
}


template<class Out>
double TPrimaryVertex<Out>::d_function1(int id_det, double hitx, double hity)
{
  if(id_det == 0 || id_det == 1)
    return 0.;

  double x = hitx - IP_recons.X();
  double y = hity - IP_recons.Y();
  double z = Zpos_Fiber[id_det] - target_pos.Z();

  double theta = TMath::ATan2( std::sqrt(x * x + y * y), z);
  double phi = TMath::ATan2( y , x );

  if(theta*180./M_PI < 5.)
    {
      double d = d_function2(id_det, hitx);
      return d;
    }

  double f = std::tan(theta) * std::cos(phi + std::get<4>(param_d_funct1[id_det]) * M_PI / 180.);

  double d = std::get<0>(param_d_funct1[id_det]) + std::get<1>(param_d_funct1[id_det]) * f + std::get<2>(param_d_funct1[id_det]) * f * f
                  + std::get<3>(param_d_funct1[id_det]) * f * f * f; //d_function_1
  return d;
}


template<class Out>
double TPrimaryVertex<Out>::d_function2(int id_det, double posx)
{
  double x = posx - IP_recons.X();

  double d = std::get<0>(param_d_funct2[id_det]) + std::get<1>(param_d_funct2[id_det]) * x + std::get<2>(param_d_funct2[id_det]) * x * x
                + std::get<3>(param_d_funct2[id_det]) * x * x * x; //d_function_2
  return d;
}


template<class Out>
void TPrimaryVertex<Out>::Get_dvalueThetaPhi_CombXUV(int id_det, std::tuple<double,double,double> CombXUV)
{
  if(id_det == 0 || id_det == 1)
    return;

  double ang_u = ang[id_det][id_mid[id_det]] * M_PI / 180.;

  double pos_x = get<0>(CombXUV);
  double pos_u = get<1>(CombXUV) / std::cos(ang_u);
  double pos_v = get<2>(CombXUV) / std::cos(ang_u);

  double tmp_hit_xu_y = 999.;
  if(ang_u>0.)
    tmp_hit_xu_y = -std::tan((M_PI/2. - ang_u))*pos_x + std::tan((M_PI/2. -ang_u))*pos_u;
  else
    tmp_hit_xu_y =  std::tan((M_PI/2. + ang_u))*pos_x - std::tan((M_PI/2. +ang_u))*pos_u;


  double tmp_hit_xv_y = -999.;
  if(ang_u>0)
    tmp_hit_xv_y =  std::tan((M_PI/2. - ang_u))*pos_x - std::tan((M_PI/2. - ang_u))*pos_v;
  else
    tmp_hit_xv_y = -std::tan((M_PI/2. + ang_u))*pos_x + std::tan((M_PI/2. + ang_u))*pos_v;

  double tmp_x = pos_x + (tmp_hit_xu_y-tmp_hit_xv_y)/2.*std::tan(ang_u/2.);
  double tmp_y = (tmp_hit_xu_y+tmp_hit_xv_y)/2.;

  double tmp_d = tmp_hit_xu_y-tmp_hit_xv_y;

  double x = tmp_x - IP_recons.X();
  double y = tmp_y - IP_recons.Y();
  double z = Zpos_Fiber[id_det] - target_pos.Z();

  double theta = TMath::ATan2( std::sqrt(x * x + y * y), z) * 180. / M_PI;
  double phi = TMath::ATan2( y , x ) * 180. / M_PI;
  
  LocalHisto.h_CombsXUV_dvalue_theta[id_det]->Fill(theta, tmp_d, 1.);
  LocalHisto.h_CombsXUV_dvalue_phi[id_det]->Fill(phi, tmp_d, 1.);

  if(theta >= 5.)
    LocalHisto.h_CombsXUV_dvalue_phi_theta5[id_det]->Fill(phi, tmp_d, 1.);
  if(theta >= 10.)
    LocalHisto.h_CombsXUV_dvalue_phi_theta10[id_det]->Fill(phi, tmp_d, 1.);

  return;
}


template<class Out>
void TPrimaryVertex<Out>::BeamTracking(const std::vector<std::tuple<double,double,double>>& CombXUV_UFT1,
                                        const std::vector<std::tuple<double,double,double>>& CombXUV_UFT2,
                                          std::vector<PrimVtxTrack>& BeamTracks)
{
  int num_UFT1 = CombXUV_UFT1.size();
  int num_UFT2 = CombXUV_UFT2.size();

  int nlayer = 6;

  std::vector<double> w;
  std::vector<double> z;
  std::vector<double> angle;
  std::vector<double> s;

  w.emplace_back(    fiber_resolution );
  w.emplace_back(    fiber_resolution );
  w.emplace_back(    fiber_resolution );
  z.emplace_back(    Zpos_Fiber[0] - 0.4);
  z.emplace_back(    Zpos_Fiber[0]      );
  z.emplace_back(    Zpos_Fiber[0] + 0.4);
  angle.emplace_back(ang[0][0] * M_PI / 180.);
  angle.emplace_back(ang[0][1] * M_PI / 180.);
  angle.emplace_back(ang[0][2] * M_PI / 180.);

  w.emplace_back(    fiber_resolution );
  w.emplace_back(    fiber_resolution );
  w.emplace_back(    fiber_resolution );
  z.emplace_back(    Zpos_Fiber[1] - 0.4);
  z.emplace_back(    Zpos_Fiber[1]      );
  z.emplace_back(    Zpos_Fiber[1] + 0.4);
  angle.emplace_back(ang[1][0] * M_PI / 180.);
  angle.emplace_back(ang[1][1] * M_PI / 180.);
  angle.emplace_back(ang[1][2] * M_PI / 180.);

  for(int i=0; i<num_UFT1; ++i)
    {
      for(int j=0; j<num_UFT2; ++j)
        {
          s.emplace_back(get<0>(CombXUV_UFT1[i]));
          s.emplace_back(get<1>(CombXUV_UFT1[i]));
          s.emplace_back(get<2>(CombXUV_UFT1[i]));

          s.emplace_back(get<0>(CombXUV_UFT2[j]));
          s.emplace_back(get<1>(CombXUV_UFT2[j]));
          s.emplace_back(get<2>(CombXUV_UFT2[j]));

          PrimVtxTrack tmp_BeamTrack = TrackFitting(nlayer, &w[0], &z[0], &angle[0], &s[0]);

          tmp_BeamTrack.SetFTHit(0, get<0>(CombXUV_UFT1[i]));
          tmp_BeamTrack.SetFTHit(1, get<1>(CombXUV_UFT1[i]));
          tmp_BeamTrack.SetFTHit(2, get<2>(CombXUV_UFT1[i]));

          tmp_BeamTrack.SetFTHit(3, get<0>(CombXUV_UFT2[j]));
          tmp_BeamTrack.SetFTHit(4, get<1>(CombXUV_UFT2[j]));
          tmp_BeamTrack.SetFTHit(5, get<2>(CombXUV_UFT2[j]));

          BeamTracks.emplace_back(tmp_BeamTrack);

          s.clear();
        }
    }

  return;
}


template<class Out>
void TPrimaryVertex<Out>::PrimaryTracking(const std::vector<std::tuple<double,double,double>>& CombXUV_UFT3,
                                            const std::vector<std::tuple<double,double,double>>& CombXUV_MFT1,
                                            const std::vector<std::tuple<double,double,double>>& CombXUV_MFT2,
                                              std::vector<PrimVtxTrack>& PrimaryTracks)
{
  int num_UFT3 = CombXUV_UFT3.size();
  int num_MFT1 = CombXUV_MFT1.size();
  int num_MFT2 = CombXUV_MFT2.size();

  int nlayer = 9;

  std::vector<double> w;
  std::vector<double> z;
  std::vector<double> angle;
  std::vector<double> s;

  w.emplace_back(    fiber_resolution );
  w.emplace_back(    fiber_resolution );
  w.emplace_back(    fiber_resolution );
  z.emplace_back(    Zpos_Fiber[2] - 0.4);
  z.emplace_back(    Zpos_Fiber[2]      );
  z.emplace_back(    Zpos_Fiber[2] + 0.4);
  angle.emplace_back(ang[2][0] * M_PI / 180.);
  angle.emplace_back(ang[2][1] * M_PI / 180.);
  angle.emplace_back(ang[2][2] * M_PI / 180.);

  w.emplace_back(    fiber_resolution );
  w.emplace_back(    fiber_resolution );
  w.emplace_back(    fiber_resolution );
  z.emplace_back(    Zpos_Fiber[3] - 0.4);
  z.emplace_back(    Zpos_Fiber[3]      );
  z.emplace_back(    Zpos_Fiber[3] + 0.4);
  angle.emplace_back(ang[3][0] * M_PI / 180.);
  angle.emplace_back(ang[3][1] * M_PI / 180.);
  angle.emplace_back(ang[3][2] * M_PI / 180.);

  w.emplace_back(    fiber_resolution );
  w.emplace_back(    fiber_resolution );
  w.emplace_back(    fiber_resolution );
  z.emplace_back(    Zpos_Fiber[4] - 0.4);
  z.emplace_back(    Zpos_Fiber[4]      );
  z.emplace_back(    Zpos_Fiber[4] + 0.4);
  angle.emplace_back(ang[4][0] * M_PI / 180.);
  angle.emplace_back(ang[4][1] * M_PI / 180.);
  angle.emplace_back(ang[4][2] * M_PI / 180.);

  for(int i=0; i<num_UFT3; ++i)
    {
      for(int j=0; j<num_MFT1; ++j)
        {
          for(int k=0; k<num_MFT2; ++k)
            {
              s.emplace_back(get<0>(CombXUV_UFT3[i]));
              s.emplace_back(get<1>(CombXUV_UFT3[i]));
              s.emplace_back(get<2>(CombXUV_UFT3[i]));

              s.emplace_back(get<0>(CombXUV_MFT1[j]));
              s.emplace_back(get<1>(CombXUV_MFT1[j]));
              s.emplace_back(get<2>(CombXUV_MFT1[j]));

              s.emplace_back(get<0>(CombXUV_MFT2[k]));
              s.emplace_back(get<2>(CombXUV_MFT2[k]));
              s.emplace_back(get<1>(CombXUV_MFT2[k]));

              PrimVtxTrack tmp_PrimaryTrack = TrackFitting(nlayer, &w[0], &z[0], &angle[0], &s[0]);

              tmp_PrimaryTrack.SetFTHit( 6, get<0>(CombXUV_UFT3[i]));
              tmp_PrimaryTrack.SetFTHit( 7, get<1>(CombXUV_UFT3[i]));
              tmp_PrimaryTrack.SetFTHit( 8, get<2>(CombXUV_UFT3[i]));

              tmp_PrimaryTrack.SetFTHit( 9, get<0>(CombXUV_MFT1[j]));
              tmp_PrimaryTrack.SetFTHit(10, get<1>(CombXUV_MFT1[j]));
              tmp_PrimaryTrack.SetFTHit(11, get<2>(CombXUV_MFT1[j]));

              tmp_PrimaryTrack.SetFTHit(12, get<0>(CombXUV_MFT1[k]));
              tmp_PrimaryTrack.SetFTHit(13, get<1>(CombXUV_MFT1[k]));
              tmp_PrimaryTrack.SetFTHit(14, get<2>(CombXUV_MFT1[k]));

              PrimaryTracks.emplace_back(tmp_PrimaryTrack);

              s.clear();
            }
        }
    }

  return;
}


template<class Out>
bool TPrimaryVertex<Out>::PrimVtxTracking(PrimVtxTrack& PrimaryVtxTrack)
{

  std::vector<std::tuple<size_t,double>> hits_FT =  PrimaryVtxTrack.GetFTHits();
  size_t nlayer = hits_FT.size();

  std::vector<double> w;
  std::vector<double> z;
  std::vector<double> angle;
  std::vector<double> s;

  for(size_t i = 0; i < hits_FT.size(); ++i)
    {
      size_t i_det = get<0>(hits_FT[i])/3;
      size_t i_lay = get<0>(hits_FT[i])%3;

      if(get<0>(hits_FT[i]) == 13)
        ++i_lay;
      else if(get<0>(hits_FT[i]) == 14)
        --i_lay;

      w.emplace_back(    fiber_resolution );
      z.emplace_back(    Zpos_Fiber[i_det] + 0.4 * (i_lay - 1));
      angle.emplace_back(ang[i_det][i_lay] * M_PI / 180.);
      s.emplace_back(    get<1>(hits_FT[i]));
    }

  PrimVtxTrack tmp_PrimaryTrack = TrackFitting(nlayer, &w[0], &z[0], &angle[0], &s[0]);

  if(tmp_PrimaryTrack.GetX() < -990.) //Fitting failed
    return false;

  PrimaryVtxTrack.SetX(tmp_PrimaryTrack.GetX());
  PrimaryVtxTrack.SetY(tmp_PrimaryTrack.GetY());
  PrimaryVtxTrack.SetA(tmp_PrimaryTrack.GetA());
  PrimaryVtxTrack.SetB(tmp_PrimaryTrack.GetB());

  if(PrimaryVtxTrack.GetNHits() == 4)
    PrimaryVtxTrack.SetChi2NDF(0.);
  else
    PrimaryVtxTrack.SetChi2NDF(tmp_PrimaryTrack.GetChi2NDF());

  return true;
}


template<class Out>
bool TPrimaryVertex<Out>::CloseToBeam(std::tuple<double,double,double,double> BeamEllipsePar, PrimVtxTrack& PrimaryVtxTrack)
{

  double a_ell = get<0>(BeamEllipsePar);
  double b_ell = get<1>(BeamEllipsePar);
  double h_ell = get<2>(BeamEllipsePar);
  double k_ell = get<3>(BeamEllipsePar);

  double m = PrimaryVtxTrack.GetB() / PrimaryVtxTrack.GetA();
  double s = PrimaryVtxTrack.GetY() - PrimaryVtxTrack.GetX() * m - k_ell;

  double radicand = std::pow(a_ell*a_ell*m*s - b_ell*b_ell*h_ell, 2)
                      - (b_ell*b_ell + a_ell*a_ell*m*m)*(s*s*a_ell + h_ell*h_ell*b_ell-a_ell*a_ell*b_ell*b_ell);

  if (radicand < 0.)
    return false;
  else
    return true;
}


template<class Out>
PrimVtxTrack TPrimaryVertex<Out>::TrackFitting(int nlayer, double* w, double* z, double* angle, double* s)
{
  PrimVtxTrack tmp_PrimVtxTrack;

  size_t n = nlayer;
  double ct[n],st[n],wt[n];

  for( std::size_t i=0; i<n; ++i ){
    double ww = w[i];
    double aa = angle[i];
    wt[i] = 1./(ww*ww);
    ct[i] = std::cos(aa);
    st[i] = std::sin(aa);
  }

  double matrx[16], *mtp[4], fitp[4];
  mtp[0]=&matrx[0]; mtp[1]=&matrx[4]; mtp[2]=&matrx[8]; mtp[3]=&matrx[12];

  for( int i=0; i<4; ++i ){
    fitp[i]=0.0;
    for( int j=0; j<4; ++j ){
      mtp[i][j]=0.0;
    }
  }

  for( std::size_t i=0; i<n; ++i ){
    double ww=wt[i], zz=z[i], ss=s[i], ctt=ct[i], stt=st[i];
    mtp[0][0] += ww*ctt*ctt;
    mtp[0][1] += ww*zz*ctt*ctt;
    mtp[0][2] += ww*ctt*stt;
    mtp[0][3] += ww*zz*ctt*stt;
    mtp[1][1] += ww*zz*zz*ctt*ctt;
    mtp[1][2] += ww*zz*ctt*stt;
    mtp[1][3] += ww*zz*zz*ctt*stt;
    mtp[2][2] += ww*stt*stt;
    mtp[2][3] += ww*zz*stt*stt;
    mtp[3][3] += ww*zz*zz*stt*stt;

    fitp[0] += ww*ss*ctt;
    fitp[1] += ww*zz*ss*ctt;
    fitp[2] += ww*ss*stt;
    fitp[3] += ww*zz*ss*stt;
  }
  mtp[1][0]=mtp[0][1]; mtp[2][0]=mtp[0][2]; mtp[3][0]=mtp[0][3];
  mtp[2][1]=mtp[1][2]; mtp[3][1]=mtp[1][3]; mtp[3][2]=mtp[2][3];

  std::vector<int> indxc(n), indxd(n), ipiv(n);

  if( GaussJordan( mtp, 4, fitp, &indxc[0],
        &indxd[0], &ipiv[0] )==false ){
    std::cerr <<  "Fitting fails" << std::endl;
    return tmp_PrimVtxTrack;
  }

  tmp_PrimVtxTrack.SetX(fitp[0] + fitp[1] * target_pos.Z());
  tmp_PrimVtxTrack.SetY(fitp[2] + fitp[3] * target_pos.Z());
  tmp_PrimVtxTrack.SetA(fitp[1]);
  tmp_PrimVtxTrack.SetB(fitp[3]);

  double chisqr = 0.;
  for( std::size_t i=0; i<n; ++i ){
    double ww=wt[i], zz=z[i];
    double scal=(fitp[0]+fitp[1]*zz)*ct[i] + (fitp[2]+fitp[3]*zz)*st[i];
    chisqr += ww*(s[i]-scal)*(s[i]-scal);
  }

  chisqr /= n-4.;
  tmp_PrimVtxTrack.SetChi2NDF(chisqr);

  return tmp_PrimVtxTrack;
}


template<class Out>
bool TPrimaryVertex<Out>::GaussJordan( double **a, int n, double *b, int *indxr, int *indxc, int *ipiv )
{

  for( int j=0; j<n; ++j ) ipiv[j]=0;
  for( int i=0; i<n; ++i ){
    double big=0.0;
    int irow=-1, icol=-1;
    for( int j=0; j<n; ++j )
      if( ipiv[j]!=1 )
        for( int k=0; k<n; ++k ){
          if( ipiv[k]==0 ){
            if( std::fabs(a[j][k])>=big ){
              big=std::fabs(a[j][k]);
              irow=j; icol=k;
            }
          }
          else if( ipiv[k]>1 ){
            return false;
          }
        }
    ++(ipiv[icol]);

    if( irow!=icol ){
      for( int k=0; k<n; ++k ){
        double ta=a[irow][k];
        a[irow][k]=a[icol][k];
        a[icol][k]=ta;
      }
      double tb=b[irow];
      b[irow]=b[icol];
      b[icol]=tb;
    }

    indxr[i]=irow; indxc[i]=icol;

    if(a[icol][icol]==0.0 || std::fabs(a[icol][icol])<1e-30 || std::isinf(a[icol][icol]) ){
      return false;
    }
    double pivinv=1./a[icol][icol];
    a[icol][icol]=1.;
    for(int k=0; k<n; ++k) a[icol][k]*=pivinv;
    b[icol]*=pivinv;
    for( int k=0; k<n; ++k ){
      if(k!=icol){
        double d=a[k][icol];
        a[k][icol]=0.;
        for( int l=0; l<n; ++l ) a[k][l] -= a[icol][l]*d;
        b[k] -= b[icol]*d;
      }
    }
  }

  for(int l=n-1; l>=0; --l){
    if( indxr[l]!=indxc[l] ){
      for(int k=0; k<n; ++k ){
        double t=a[k][indxr[l]];
        a[k][indxr[l]]=a[k][indxc[l]];
        a[k][indxc[l]]=t;
      }
    }
  }
  return true;
}


template<class Out>
std::tuple<double,double,double,double> TPrimaryVertex<Out>::BeamEllipseParam(PrimVtxTrack& BeamTrack)
{
  return std::make_tuple(50.*fiber_resolution, 50.*fiber_resolution, BeamTrack.GetX(), BeamTrack.GetY());//Change!
}


double PrimVtxTrack::GetTheta()
{
  if(a > -990. && b > -990.)
    return TMath::ATan2(std::sqrt( a*a + b*b ), 1.) * 180. / M_PI;
  else
    return -999.;
}

double PrimVtxTrack::GetPhi()
{
  if(a > -990. && b > -990.)
    return TMath::ATan2(b, a) * 180. / M_PI;
  else
    return -999.;
}




template class TPrimaryVertex<MCAnaEventG4Sol>;
template class TPrimaryVertex<Ana_WasaEvent>;
