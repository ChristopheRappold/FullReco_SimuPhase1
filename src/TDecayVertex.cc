#include "TDecayVertex.h"

#include "FullRecoEvent.hh"
#include "ReturnRes.hh"

#include <tuple>
#include <string>

#include "TLorentzVector.h"
#include "TVector3.h"

//#define DEBUG_DECAYVTX

#define REAL_PIONS_CHECK
#define CUT_PIONS_CHECK
#define CENTROID_METHOD
#define VFUNCTION_METHOD

using namespace std;
using namespace G4Sol;

TDecayVertex::TDecayVertex(const THyphiAttributes& attribut)
    : TDataProcessInterface("DecayVertexReco"), att(attribut)
{
  //rand = new TRandom3(); CHANGE
  //PDG_fromName pid_fromName ; //Is it necessary ?
}

TDecayVertex::~TDecayVertex() {}

void TDecayVertex::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }

ReturnRes::InfoM TDecayVertex::operator()(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
{

  int result_finder = Exec(RecoEvent, OutTree);

  return SoftExit(result_finder);
}

int TDecayVertex::Exec(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree) { return FinderDecayVertex(RecoEvent); }

ReturnRes::InfoM TDecayVertex::SoftExit(int result_full) {
  
  if(result_full == -1)
    {
      att._logger->debug("No real/reconstructed fragment tracks for decay vertex");
      LocalHisto.h_DecayVtxstats->Fill("FragmentTracks=0", 1.);
      //return ReturnRes::DecayVtxError; //Check in the future
      return ReturnRes::Fine;
    }

  else if(result_full == -2)
    {
      att._logger->debug("No real pion tracks for decay vertex");
      LocalHisto.h_DecayVtxstats->Fill("RealPionTracks=0", 1.);
      //return ReturnRes::DecayVtxError; //Check in the future
      return ReturnRes::Fine;
    }

  else if(result_full == -3)
    {
      att._logger->debug("No cut pion tracks for decay vertex");
      LocalHisto.h_DecayVtxstats->Fill("CutPionTracks=0", 1.);
      //return ReturnRes::DecayVtxError; //Check in the future
      return ReturnRes::Fine;
    }

  else if(result_full == -4)
    {
      att._logger->debug("No pion tracks reconstructed for decay vertex");
      LocalHisto.h_DecayVtxstats->Fill("PionTracks=0", 1.);
      //return ReturnRes::DecayVtxError; //Check in the future
      return ReturnRes::Fine;
    }

  LocalHisto.h_DecayVtxstats->Fill("Fine", 1.);

  return ReturnRes::Fine; 
}



void TDecayVertex::SelectHists()
{
  LocalHisto.h_Pt_fragments = AnaHisto->CloneAndRegister(AnaHisto->h_Pt_fragments);
  LocalHisto.h_Pz_fragments = AnaHisto->CloneAndRegister(AnaHisto->h_Pz_fragments);

  LocalHisto.h_Pt_pions = AnaHisto->CloneAndRegister(AnaHisto->h_Pt_pions);
  LocalHisto.h_Pz_pions = AnaHisto->CloneAndRegister(AnaHisto->h_Pz_pions);
  LocalHisto.h_Chi2ndf_pions = AnaHisto->CloneAndRegister(AnaHisto->h_Chi2ndf_pions);

  LocalHisto.h_Pt_realpions = AnaHisto->CloneAndRegister(AnaHisto->h_Pt_realpions);
  LocalHisto.h_Pz_realpions = AnaHisto->CloneAndRegister(AnaHisto->h_Pz_realpions);

  LocalHisto.h_Pt_cutpions = AnaHisto->CloneAndRegister(AnaHisto->h_Pt_cutpions);
  LocalHisto.h_Pz_cutpions = AnaHisto->CloneAndRegister(AnaHisto->h_Pz_cutpions);

  LocalHisto.h_Nrealpions = AnaHisto->CloneAndRegister(AnaHisto->h_Nrealpions);
  LocalHisto.h_Ncutpions = AnaHisto->CloneAndRegister(AnaHisto->h_Ncutpions);
  LocalHisto.h_Npions = AnaHisto->CloneAndRegister(AnaHisto->h_Npions);

  LocalHisto.h_Closedist_Distance = AnaHisto->CloneAndRegister(AnaHisto->h_Closedist_Distance);
  LocalHisto.h_Closedist_PosZ = AnaHisto->CloneAndRegister(AnaHisto->h_Closedist_PosZ);
  LocalHisto.h_Dist_DecayTrackPrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_Dist_DecayTrackPrimVtx);


  LocalHisto.h_Closedist_realDistance = AnaHisto->CloneAndRegister(AnaHisto->h_Closedist_realDistance);
  LocalHisto.h_Closedist_realPosZ = AnaHisto->CloneAndRegister(AnaHisto->h_Closedist_realPosZ);
  LocalHisto.h_Dist_realDecayTrackPrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_Dist_realDecayTrackPrimVtx);

  LocalHisto.h_Closedist_cutDistance = AnaHisto->CloneAndRegister(AnaHisto->h_Closedist_cutDistance);
  LocalHisto.h_Closedist_cutPosZ = AnaHisto->CloneAndRegister(AnaHisto->h_Closedist_cutPosZ);
  LocalHisto.h_Dist_cutDecayTrackPrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_Dist_cutDecayTrackPrimVtx);


  LocalHisto.h_DecayVertexDistance = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistance);
  LocalHisto.h_DecayVertexDistanceX = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceX);
  LocalHisto.h_DecayVertexDistanceY = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceY);
  LocalHisto.h_DecayVertexDistanceZ = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceZ);

  LocalHisto.h_DecayVertexDistance_centroid = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistance_centroid);
  LocalHisto.h_DecayVertexDistanceX_centroid = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceX_centroid);
  LocalHisto.h_DecayVertexDistanceY_centroid = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceY_centroid);
  LocalHisto.h_DecayVertexDistanceZ_centroid = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceZ_centroid);

  LocalHisto.h_DecayVertexDistance_2centroid_average = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistance_2centroid_average);
  LocalHisto.h_DecayVertexDistanceX_2centroid_average = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceX_2centroid_average);
  LocalHisto.h_DecayVertexDistanceY_2centroid_average = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceY_2centroid_average);
  LocalHisto.h_DecayVertexDistanceZ_2centroid_average = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceZ_2centroid_average);

  LocalHisto.h_DecayVertexDistance_2centroid_closest = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistance_2centroid_closest);
  LocalHisto.h_DecayVertexDistanceX_2centroid_closest = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceX_2centroid_closest);
  LocalHisto.h_DecayVertexDistanceY_2centroid_closest = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceY_2centroid_closest);
  LocalHisto.h_DecayVertexDistanceZ_2centroid_closest = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceZ_2centroid_closest);

  LocalHisto.h_DecayVertexDistance_2centroid_IPCheck = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistance_2centroid_IPCheck);
  LocalHisto.h_DecayVertexDistanceX_2centroid_IPCheck = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceX_2centroid_IPCheck);
  LocalHisto.h_DecayVertexDistanceY_2centroid_IPCheck = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceY_2centroid_IPCheck);
  LocalHisto.h_DecayVertexDistanceZ_2centroid_IPCheck = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceZ_2centroid_IPCheck);

  LocalHisto.h_DecayVertexrealDistance = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexrealDistance);
  LocalHisto.h_DecayVertexrealDistanceX = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexrealDistanceX);
  LocalHisto.h_DecayVertexrealDistanceY = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexrealDistanceY);
  LocalHisto.h_DecayVertexrealDistanceZ = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexrealDistanceZ);

  LocalHisto.h_DecayVertexcutDistance = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexcutDistance);
  LocalHisto.h_DecayVertexcutDistanceX = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexcutDistanceX);
  LocalHisto.h_DecayVertexcutDistanceY = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexcutDistanceY);
  LocalHisto.h_DecayVertexcutDistanceZ = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexcutDistanceZ);

  LocalHisto.h_DecayVertexPosZ_real = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexPosZ_real);
  LocalHisto.h_DecayVertexPosZ_vfunction = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexPosZ_vfunction);
  LocalHisto.h_DecayVertexPosZ_centroid = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexPosZ_centroid);

  LocalHisto.h_HypInvariantMass = AnaHisto->CloneAndRegister(AnaHisto->h_HypInvariantMass);
  
  LocalHisto.h_DecayVtxstats = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVtxstats);
}

int TDecayVertex::FinderDecayVertex(FullRecoEvent& RecoEvent)
{

  double DecayVertex_real_X = RecoEvent.DecayVertex[0];
  double DecayVertex_real_Y = RecoEvent.DecayVertex[1];
  double DecayVertex_real_Z = RecoEvent.DecayVertex[2];

  LocalHisto.h_DecayVertexPosZ_real->Fill(DecayVertex_real_Z, 1.);

  //Fragment tracks
  RealTracksFinder(RecoEvent.TrackDAFSim, He3_pdg, 0, RecoEvent.FragmentTracks);

  if(RecoEvent.FragmentTracks.size() == 0)
    return -1;

  for(size_t i = 0; i < RecoEvent.FragmentTracks.size(); ++i)
    {
      LocalHisto.h_Pt_fragments->Fill(sqrt(pow(RecoEvent.FragmentTracks[i].Hit_MomEnergy.Px(),2.)
                                            + pow(RecoEvent.FragmentTracks[i].Hit_MomEnergy.Py(),2.)), 1.);
      LocalHisto.h_Pz_fragments->Fill(RecoEvent.FragmentTracks[i].Hit_MomEnergy.Pz(), 1.);
    }
  

#ifdef REAL_PIONS_CHECK

  //Real pion tracks
  std::vector<DecayTrackInfo> RealPionTracks {};
  RealTracksFinder(RecoEvent.TrackDAFSim, pi_pdg, 0, RealPionTracks);

  LocalHisto.h_Nrealpions->Fill(RealPionTracks.size(), 1.);
  if(RealPionTracks.size() == 0)
    return -2;

  double closedist_realdistance = 0.;
  TVector3 closedist_realpos;
  double dist_realDecayTrackPrimVtx = 0.;

  std::unordered_map<int, InfoInit>::iterator itr_real;
  bool ifDaughter_real = false;

  for(size_t i = 0; i < RealPionTracks.size(); ++i)
    {
      LocalHisto.h_Pt_realpions->Fill(sqrt(pow(RealPionTracks[i].Hit_MomEnergy.Px(),2.)
                                        + pow(RealPionTracks[i].Hit_MomEnergy.Py(),2.)), 1.);
      LocalHisto.h_Pz_realpions->Fill(RealPionTracks[i].Hit_MomEnergy.Pz(), 1.);

      CloseDist(RecoEvent.FragmentTracks[0], RealPionTracks[i], closedist_realdistance, closedist_realpos);
      LocalHisto.h_Closedist_realDistance->Fill(closedist_realdistance, 1.);
      LocalHisto.h_Closedist_realPosZ->Fill(closedist_realpos.Z(), 1.);

      Dist_DecayTrackPrimVtx(RealPionTracks[i], RecoEvent.PrimVtxRecons, dist_realDecayTrackPrimVtx);
      LocalHisto.h_Dist_realDecayTrackPrimVtx->Fill(dist_realDecayTrackPrimVtx, "All", 1.);

      for(itr_real = RecoEvent.DaughtersTrackDAFInit.begin(); itr_real != RecoEvent.DaughtersTrackDAFInit.end(); ++itr_real)
        {
          if(itr_real->first == RealPionTracks[i].Id)
            ifDaughter_real = true;
        }

      if(ifDaughter_real)
        LocalHisto.h_Dist_realDecayTrackPrimVtx->Fill(dist_realDecayTrackPrimVtx, "Daughters", 1.);
      else
        LocalHisto.h_Dist_realDecayTrackPrimVtx->Fill(dist_realDecayTrackPrimVtx, "Primaries", 1.);
    }

  //Decay vertex reconstruction
  TVector3 DecayVertexReconsreal;

  TrackstoDecayVertex(RecoEvent.FragmentTracks, RealPionTracks, RecoEvent.PrimVtxRecons, DecayVertexReconsreal);

  /*
  std::cout << "X:\t" << DecayVertex_real_X << "\t" << RecoEvent.DecayVtxRecons.X() << "\n";
  std::cout << "Y:\t" << DecayVertex_real_Y << "\t" << RecoEvent.DecayVtxRecons.Y() << "\n";
  std::cout << "Z:\t" << DecayVertex_real_Z << "\t" << RecoEvent.DecayVtxRecons.Z() << "\n";
  std::cout << "\n";
  */


  double realdistance  = sqrt(pow((DecayVertex_real_X - DecayVertexReconsreal.X()), 2.) +
                              pow((DecayVertex_real_Y - DecayVertexReconsreal.Y()), 2.) +
                              pow((DecayVertex_real_Z - DecayVertexReconsreal.Z()), 2.));
  double realdistanceX = DecayVertex_real_X - DecayVertexReconsreal.X();
  double realdistanceY = DecayVertex_real_Y - DecayVertexReconsreal.Y();
  double realdistanceZ = DecayVertex_real_Z - DecayVertexReconsreal.Z();

  LocalHisto.h_DecayVertexrealDistance->Fill(realdistance, 1.);
  LocalHisto.h_DecayVertexrealDistanceX->Fill(realdistanceX, 1.);
  LocalHisto.h_DecayVertexrealDistanceY->Fill(realdistanceY, 1.);
  LocalHisto.h_DecayVertexrealDistanceZ->Fill(realdistanceZ, 1.);

#endif

#ifdef CUT_PIONS_CHECK

 //Cut pion tracks
  std::vector<DecayTrackInfo> CutPionTracks {};
  RealTracksFinder(RecoEvent.TrackDAFSim, pi_pdg, 1, CutPionTracks);

  LocalHisto.h_Ncutpions->Fill(CutPionTracks.size(), 1.);
  if(CutPionTracks.size() == 0)
    return -3;

  double closedist_cutdistance = 0.;
  TVector3 closedist_cutpos;
  double dist_cutDecayTrackPrimVtx = 0.;

  std::unordered_map<int, InfoInit>::iterator itr_cut;
  bool ifDaughter_cut = false;

  for(size_t i = 0; i < CutPionTracks.size(); ++i)
    {
      LocalHisto.h_Pt_cutpions->Fill(sqrt(pow(CutPionTracks[i].Hit_MomEnergy.Px(),2.)
                                        + pow(CutPionTracks[i].Hit_MomEnergy.Py(),2.)), 1.);
      LocalHisto.h_Pz_cutpions->Fill(CutPionTracks[i].Hit_MomEnergy.Pz(), 1.);

      CloseDist(RecoEvent.FragmentTracks[0], CutPionTracks[i], closedist_cutdistance, closedist_cutpos);
      LocalHisto.h_Closedist_cutDistance->Fill(closedist_cutdistance, 1.);
      LocalHisto.h_Closedist_cutPosZ->Fill(closedist_cutpos.Z(), 1.);

      Dist_DecayTrackPrimVtx(CutPionTracks[i], RecoEvent.PrimVtxRecons, dist_cutDecayTrackPrimVtx);
      LocalHisto.h_Dist_cutDecayTrackPrimVtx->Fill(dist_cutDecayTrackPrimVtx, "All", 1.);

      for(itr_cut = RecoEvent.DaughtersTrackDAFInit.begin(); itr_cut != RecoEvent.DaughtersTrackDAFInit.end(); ++itr_cut)
        {
          if(itr_cut->first == CutPionTracks[i].Id)
            ifDaughter_cut = true;
        }

      if(ifDaughter_cut)
        LocalHisto.h_Dist_cutDecayTrackPrimVtx->Fill(dist_cutDecayTrackPrimVtx, "Daughters", 1.);
      else
        LocalHisto.h_Dist_cutDecayTrackPrimVtx->Fill(dist_cutDecayTrackPrimVtx, "Primaries", 1.);
    }

  //Decay vertex reconstruction
  TVector3 DecayVertexReconscut;

  TrackstoDecayVertex(RecoEvent.FragmentTracks, CutPionTracks, RecoEvent.PrimVtxRecons, DecayVertexReconscut);

  /*
  std::cout << "X:\t" << DecayVertex_cut_X << "\t" << RecoEvent.DecayVtxRecons.X() << "\n";
  std::cout << "Y:\t" << DecayVertex_cut_Y << "\t" << RecoEvent.DecayVtxRecons.Y() << "\n";
  std::cout << "Z:\t" << DecayVertex_cut_Z << "\t" << RecoEvent.DecayVtxRecons.Z() << "\n";
  std::cout << "\n";
  */


  double cutdistance  = sqrt(pow((DecayVertex_real_X - DecayVertexReconscut.X()), 2.) +
                              pow((DecayVertex_real_Y - DecayVertexReconscut.Y()), 2.) +
                              pow((DecayVertex_real_Z - DecayVertexReconscut.Z()), 2.));
  double cutdistanceX = DecayVertex_real_X - DecayVertexReconscut.X();
  double cutdistanceY = DecayVertex_real_Y - DecayVertexReconscut.Y();
  double cutdistanceZ = DecayVertex_real_Z - DecayVertexReconscut.Z();

  LocalHisto.h_DecayVertexcutDistance->Fill(cutdistance, 1.);
  LocalHisto.h_DecayVertexcutDistanceX->Fill(cutdistanceX, 1.);
  LocalHisto.h_DecayVertexcutDistanceY->Fill(cutdistanceY, 1.);
  LocalHisto.h_DecayVertexcutDistanceZ->Fill(cutdistanceZ, 1.);

#endif

  //Pion tracks
  PionTracksFinder(RecoEvent.DAF_results, RecoEvent.PionTracks);

  LocalHisto.h_Npions->Fill(RecoEvent.PionTracks.size(), 1.);
  if(RecoEvent.PionTracks.size() == 0)
    return -4;

  double closedist_distance = 0.;
  TVector3 temp_closedist_pos;

  double newclosest_distance = 1000.;
  TVector3 new_closedist_pos;

  std::vector<TVector3> vect_closedist_pos;
  double dist_DecayTrackPrimVtx = 0.;

  std::unordered_map<int, InfoInit>::iterator itr_recons;
  bool ifDaughter_recons = false;

  for(size_t i = 0; i < RecoEvent.PionTracks.size(); ++i)
    {
      LocalHisto.h_Pt_pions->Fill(sqrt(pow(RecoEvent.PionTracks[i].Hit_MomEnergy.Px(),2.)
                                        + pow(RecoEvent.PionTracks[i].Hit_MomEnergy.Py(),2.)), 1.);
      LocalHisto.h_Pz_pions->Fill(RecoEvent.PionTracks[i].Hit_MomEnergy.Pz(), 1.);

      CloseDist(RecoEvent.FragmentTracks[0], RecoEvent.PionTracks[i], closedist_distance, temp_closedist_pos);
      vect_closedist_pos.emplace_back(temp_closedist_pos);

      LocalHisto.h_Closedist_Distance->Fill(closedist_distance, 1.);
      LocalHisto.h_Closedist_PosZ->Fill(temp_closedist_pos.Z(), 1.);

      Dist_DecayTrackPrimVtx(RecoEvent.PionTracks[i], RecoEvent.PrimVtxRecons, dist_DecayTrackPrimVtx);
      LocalHisto.h_Dist_DecayTrackPrimVtx->Fill(dist_DecayTrackPrimVtx, "All", 1.);

      for(itr_recons = RecoEvent.DaughtersTrackDAFInit.begin(); itr_recons != RecoEvent.DaughtersTrackDAFInit.end(); ++itr_recons)
        {
          if(itr_recons->first == RecoEvent.PionTracks[i].Id)
            ifDaughter_recons = true;
        }

      if(ifDaughter_recons)
        LocalHisto.h_Dist_DecayTrackPrimVtx->Fill(dist_DecayTrackPrimVtx, "Daughters", 1.);
      else
        LocalHisto.h_Dist_DecayTrackPrimVtx->Fill(dist_DecayTrackPrimVtx, "Primaries", 1.);

      if(closedist_distance < newclosest_distance)
        {
          newclosest_distance = closedist_distance;
          new_closedist_pos = temp_closedist_pos;
        }
    }


#ifdef CENTROID_METHOD

  TVector3 closedist_pos;

  for(size_t i = 0; i < vect_closedist_pos.size(); ++i)
    {
      closedist_pos += vect_closedist_pos[i];
    }
  closedist_pos *= (1./vect_closedist_pos.size());

  double distance_centroid  = sqrt(pow((DecayVertex_real_X - closedist_pos.X()), 2.) +
                                   pow((DecayVertex_real_Y - closedist_pos.Y()), 2.) +
                                   pow((DecayVertex_real_Z - closedist_pos.Z()), 2.));
  double distanceX_centroid = DecayVertex_real_X - closedist_pos.X();
  double distanceY_centroid = DecayVertex_real_Y - closedist_pos.Y();
  double distanceZ_centroid = DecayVertex_real_Z - closedist_pos.Z();

  LocalHisto.h_DecayVertexDistance_centroid->Fill(distance_centroid, 1.);
  LocalHisto.h_DecayVertexDistanceX_centroid->Fill(distanceX_centroid, 1.);
  LocalHisto.h_DecayVertexDistanceY_centroid->Fill(distanceY_centroid, 1.);
  LocalHisto.h_DecayVertexDistanceZ_centroid->Fill(distanceZ_centroid, 1.);

  LocalHisto.h_DecayVertexPosZ_centroid->Fill(closedist_pos.Z(), 1.);

  if(RecoEvent.PionTracks.size() > 1)
    {
      LocalHisto.h_DecayVertexDistance_2centroid_average->Fill(distance_centroid, 1.);
      LocalHisto.h_DecayVertexDistanceX_2centroid_average->Fill(distanceX_centroid, 1.);
      LocalHisto.h_DecayVertexDistanceY_2centroid_average->Fill(distanceY_centroid, 1.);
      LocalHisto.h_DecayVertexDistanceZ_2centroid_average->Fill(distanceZ_centroid, 1.);

      double distance_2centroid_closest  = sqrt(pow((DecayVertex_real_X - new_closedist_pos.X()), 2.) +
                                                pow((DecayVertex_real_Y - new_closedist_pos.Y()), 2.) +
                                                pow((DecayVertex_real_Z - new_closedist_pos.Z()), 2.));
      double distanceX_2centroid_closest = DecayVertex_real_X - new_closedist_pos.X();
      double distanceY_2centroid_closest = DecayVertex_real_Y - new_closedist_pos.Y();
      double distanceZ_2centroid_closest = DecayVertex_real_Z - new_closedist_pos.Z();

      LocalHisto.h_DecayVertexDistance_2centroid_closest->Fill(distance_2centroid_closest, 1.);
      LocalHisto.h_DecayVertexDistanceX_2centroid_closest->Fill(distanceX_2centroid_closest, 1.);
      LocalHisto.h_DecayVertexDistanceY_2centroid_closest->Fill(distanceY_2centroid_closest, 1.);
      LocalHisto.h_DecayVertexDistanceZ_2centroid_closest->Fill(distanceZ_2centroid_closest, 1.);


      double dist_2DecayTrackPrimVtx = 0.;

      double closedist_2distance = 0.;
      TVector3 temp_closedist_2pos;

      double newclosest_2distance = 1000.;
      TVector3 new_closedist_2pos;

      for(size_t i = 0; i < RecoEvent.PionTracks.size(); ++i)
        {
          Dist_DecayTrackPrimVtx(RecoEvent.PionTracks[i], RecoEvent.PrimVtxRecons, dist_2DecayTrackPrimVtx);

          if(dist_2DecayTrackPrimVtx <= MinDist_DecayTracksPrimVtx)
            continue;

          CloseDist(RecoEvent.FragmentTracks[0], RecoEvent.PionTracks[i], closedist_2distance, temp_closedist_2pos);

          if(closedist_2distance < newclosest_2distance)
            {
              newclosest_2distance = closedist_2distance;
              new_closedist_2pos = temp_closedist_2pos;
            }

        }

      double distance_2centroid_IPCheck  = sqrt(pow((DecayVertex_real_X - new_closedist_2pos.X()), 2.) +
                                                pow((DecayVertex_real_Y - new_closedist_2pos.Y()), 2.) +
                                                pow((DecayVertex_real_Z - new_closedist_2pos.Z()), 2.));
      double distanceX_2centroid_IPCheck = DecayVertex_real_X - new_closedist_2pos.X();
      double distanceY_2centroid_IPCheck = DecayVertex_real_Y - new_closedist_2pos.Y();
      double distanceZ_2centroid_IPCheck = DecayVertex_real_Z - new_closedist_2pos.Z();

      LocalHisto.h_DecayVertexDistance_2centroid_IPCheck->Fill(distance_2centroid_IPCheck, 1.);
      LocalHisto.h_DecayVertexDistanceX_2centroid_IPCheck->Fill(distanceX_2centroid_IPCheck, 1.);
      LocalHisto.h_DecayVertexDistanceY_2centroid_IPCheck->Fill(distanceY_2centroid_IPCheck, 1.);
      LocalHisto.h_DecayVertexDistanceZ_2centroid_IPCheck->Fill(distanceZ_2centroid_IPCheck, 1.);
    }

#endif

#ifdef VFUNCTION_METHOD

  //Decay vertex reconstruction
  TVector3 DecayVertexRecons;

  TrackstoDecayVertex(RecoEvent.FragmentTracks, RecoEvent.PionTracks, RecoEvent.PrimVtxRecons, DecayVertexRecons);
  RecoEvent.DecayVtxRecons.SetXYZ(DecayVertexRecons.X(), DecayVertexRecons.Y(), DecayVertexRecons.Z());

  /*
  std::cout << "X:\t" << DecayVertex_real_X << "\t" << RecoEvent.DecayVtxRecons.X() << "\n";
  std::cout << "Y:\t" << DecayVertex_real_Y << "\t" << RecoEvent.DecayVtxRecons.Y() << "\n";
  std::cout << "Z:\t" << DecayVertex_real_Z << "\t" << RecoEvent.DecayVtxRecons.Z() << "\n";
  std::cout << "\n";
  */


  double distance  = sqrt(pow((DecayVertex_real_X - DecayVertexRecons.X()), 2.) +
                          pow((DecayVertex_real_Y - DecayVertexRecons.Y()), 2.) +
                          pow((DecayVertex_real_Z - DecayVertexRecons.Z()), 2.));
  double distanceX = DecayVertex_real_X - DecayVertexRecons.X();
  double distanceY = DecayVertex_real_Y - DecayVertexRecons.Y();
  double distanceZ = DecayVertex_real_Z - DecayVertexRecons.Z();

  LocalHisto.h_DecayVertexDistance->Fill(distance, 1.);
  LocalHisto.h_DecayVertexDistanceX->Fill(distanceX, 1.);
  LocalHisto.h_DecayVertexDistanceY->Fill(distanceY, 1.);
  LocalHisto.h_DecayVertexDistanceZ->Fill(distanceZ, 1.);

  LocalHisto.h_DecayVertexPosZ_vfunction->Fill(DecayVertexRecons.Z(), 1.);

#endif

  //Hypernucleus reconstruction
  TLorentzVector Hypernucleus_MomEnergy = RecoEvent.FragmentTracks[0].Hit_MomEnergy + RecoEvent.PionTracks[0].Hit_MomEnergy;  
  LocalHisto.h_HypInvariantMass->Fill(Hypernucleus_MomEnergy.M(), 1.);

  return 0;
}


void TDecayVertex::RealTracksFinder(std::unordered_map<int, std::vector<std::vector<SimHit> > >& TrackDAFSim,
                                        int& pdgParticle, int cutConditions,
                                        std::vector<DecayTrackInfo>& RealTracks)
{

/*
  char temp_namePart[namePart.size() + 1];
  namePart.copy(temp_namePart, namePart.size() + 1);
  temp_namePart[namePart.size()] = '\0';


  auto tempPart = TDatabasePDG::Instance()->GetParticle(temp_namePart);
  int pdgParticle = tempPart->PdgCode();
*/

  std::unordered_map<int, std::vector<std::vector<SimHit> > >::iterator itr;
  for(itr = TrackDAFSim.begin(); itr != TrackDAFSim.end(); ++itr)
    {
      size_t iDetFirst = -1;
      
      size_t nHits_MDC = 0;
      size_t nHits_MiniFiber = 0;
      size_t nHits_PSCE = 0;
      size_t nHits_PSBE = 0;

      for(size_t iDet = 0; iDet < itr->second.size(); ++iDet)
        {
          if(itr->second[iDet].size() == 0)
            continue;
          if(iDetFirst == -1)
            iDetFirst = iDet;
          if(iDet >= G4Sol::MG01 && iDet <= G4Sol::MG17)
            nHits_MDC += 1;
          if(iDet >= G4Sol::MiniFiberD1_x1 && iDet <= G4Sol::MiniFiberD1_v2)
            nHits_MiniFiber += 1;
          if(iDet == G4Sol::PSCE)
            nHits_PSCE += 1;
          if(iDet == G4Sol::PSBE)
            nHits_PSBE += 1;
        }

      if(iDetFirst == -1)
        continue;

      if(itr->second[iDetFirst][0].pdg == pdgParticle)
        {
          DecayTrackInfo temp_track;

          temp_track.Pdg = itr->second[iDetFirst][0].pdg;
          temp_track.Id = itr->first;

          TVector3 temp_Hit_Pos(TVector3(itr->second[iDetFirst][0].hitX, itr->second[iDetFirst][0].hitY, itr->second[iDetFirst][0].hitZ));
          temp_track.Hit_Pos = temp_Hit_Pos;

          TLorentzVector temp_Hit_MomEnergy;
          temp_Hit_MomEnergy.SetXYZM(itr->second[iDetFirst][0].momX, itr->second[iDetFirst][0].momY,
                                     itr->second[iDetFirst][0].momZ, itr->second[iDetFirst][0].mass);
          temp_track.Hit_MomEnergy = temp_Hit_MomEnergy;

          if(cutConditions == 0)
            {
              RealTracks.emplace_back(temp_track);
            }
          else if((cutConditions == 1) && (nHits_MDC >= 6) && (nHits_MiniFiber >= 4) && ((nHits_PSCE != 0) || (nHits_PSBE != 0)))
            {
              RealTracks.emplace_back(temp_track);
            }
        }
    }
}


void TDecayVertex::PionTracksFinder(std::unordered_map<int, ResSolDAF>& DAF_results,
                                    std::vector<DecayTrackInfo>& PionTracks)
{
  std::unordered_map<int, ResSolDAF>::iterator itr;
  for(itr = DAF_results.begin(); itr != DAF_results.end(); ++itr)
    {
      if((itr->second.charge == -1) && (itr->second.chi2 / itr->second.ndf < 3.) && (itr->second.Ncentral >= 6))
        {
          LocalHisto.h_Chi2ndf_pions->Fill(itr->second.chi2 / itr->second.ndf, 1.);

          DecayTrackInfo temp_track;

          temp_track.Pdg = itr->second.pdg_guess;
          temp_track.Id = itr->first;
          temp_track.Chi2 = itr->second.chi2;

          TVector3 temp_Hit_Pos = TVector3(itr->second.posX, itr->second.posY, itr->second.posZ);
          temp_track.Hit_Pos = temp_Hit_Pos;

          TLorentzVector temp_Hit_MomEnergy;
          temp_Hit_MomEnergy.SetXYZM(itr->second.momX, itr->second.momY, itr->second.momZ, pi_mass);
          temp_track.Hit_MomEnergy = temp_Hit_MomEnergy;

          PionTracks.emplace_back(temp_track);
        }
    }
}

void TDecayVertex::CloseDist(DecayTrackInfo& FragmentTrack, DecayTrackInfo& PionTrack, double& distance, TVector3& centroid)
{
  TVector3 n(FragmentTrack.Hit_MomEnergy.Py() * PionTrack.Hit_MomEnergy.Pz() - FragmentTrack.Hit_MomEnergy.Pz() * PionTrack.Hit_MomEnergy.Py(),
             FragmentTrack.Hit_MomEnergy.Pz() * PionTrack.Hit_MomEnergy.Px() - FragmentTrack.Hit_MomEnergy.Px() * PionTrack.Hit_MomEnergy.Pz(),
             FragmentTrack.Hit_MomEnergy.Px() * PionTrack.Hit_MomEnergy.Py() - FragmentTrack.Hit_MomEnergy.Py() * PionTrack.Hit_MomEnergy.Px());

  TVector3 n1(FragmentTrack.Hit_MomEnergy.Py() * n.Z() - FragmentTrack.Hit_MomEnergy.Pz() * n.Y(),
              FragmentTrack.Hit_MomEnergy.Pz() * n.X() - FragmentTrack.Hit_MomEnergy.Px() * n.Z(),
              FragmentTrack.Hit_MomEnergy.Px() * n.Y() - FragmentTrack.Hit_MomEnergy.Py() * n.X());

  TVector3 n2(PionTrack.Hit_MomEnergy.Py() * n.Z() - PionTrack.Hit_MomEnergy.Pz() * n.Y(),
              PionTrack.Hit_MomEnergy.Pz() * n.X() - PionTrack.Hit_MomEnergy.Px() * n.Z(),
              PionTrack.Hit_MomEnergy.Px() * n.Y() - PionTrack.Hit_MomEnergy.Py() * n.X());

  double FragmentFactor = ((PionTrack.Hit_Pos.X() - FragmentTrack.Hit_Pos.X()) * n2.X()
                            + (PionTrack.Hit_Pos.Y() - FragmentTrack.Hit_Pos.Y()) * n2.Y()
                            + (PionTrack.Hit_Pos.Z() - FragmentTrack.Hit_Pos.Z()) * n2.Z()) /
                                  (FragmentTrack.Hit_MomEnergy.Px() * n2.X() + FragmentTrack.Hit_MomEnergy.Py() * n2.Y()
                                   + FragmentTrack.Hit_MomEnergy.Pz() * n2.Z());

  double PionFactor = -((PionTrack.Hit_Pos.X() - FragmentTrack.Hit_Pos.X()) * n1.X()
                         + (PionTrack.Hit_Pos.Y() - FragmentTrack.Hit_Pos.Y()) * n1.Y()
                         + (PionTrack.Hit_Pos.Z() - FragmentTrack.Hit_Pos.Z()) * n1.Z()) /
                                (PionTrack.Hit_MomEnergy.Px() * n1.X() + PionTrack.Hit_MomEnergy.Py() * n1.Y()
                                 + PionTrack.Hit_MomEnergy.Pz() * n1.Z());
  
  TVector3 c1(FragmentTrack.Hit_Pos.X() + FragmentFactor * FragmentTrack.Hit_MomEnergy.Px(),
              FragmentTrack.Hit_Pos.Y() + FragmentFactor * FragmentTrack.Hit_MomEnergy.Py(),
              FragmentTrack.Hit_Pos.Z() + FragmentFactor * FragmentTrack.Hit_MomEnergy.Pz());
  
  TVector3 c2(PionTrack.Hit_Pos.X() + PionFactor * PionTrack.Hit_MomEnergy.Px(),
              PionTrack.Hit_Pos.Y() + PionFactor * PionTrack.Hit_MomEnergy.Py(),
              PionTrack.Hit_Pos.Z() + PionFactor * PionTrack.Hit_MomEnergy.Pz());;

  distance = (c2 - c1).Mag();
  centroid = (c1 + c2);
  centroid *= 0.5;
}

double TDecayVertex::f_function(DecayTrackInfo& DecayTrack, TVector3& PosXYZ)
{
  double slope_x     = DecayTrack.Hit_MomEnergy.Px() / DecayTrack.Hit_MomEnergy.Pz();
  double intercept_x = DecayTrack.Hit_Pos.X() - slope_x * DecayTrack.Hit_Pos.Z();
  double slope_y     = DecayTrack.Hit_MomEnergy.Py() / DecayTrack.Hit_MomEnergy.Pz();
  double intercept_y = DecayTrack.Hit_Pos.Y() - slope_y * DecayTrack.Hit_Pos.Z();

  double distanceStepX = 2. * boxDistXY / static_cast<double>(NstepsdiscretXY - 1); //Change(?)
  double sigma2        = pow(distanceStepX, 2.) / 12.;

  double f = exp(-0.5 *
                 (pow((PosXYZ.X() - slope_x * PosXYZ.Z() - intercept_x), 2.) +
                  pow((PosXYZ.Y() - slope_y * PosXYZ.Z() - intercept_y), 2.)) /
                 sigma2);
  return f;
}

double TDecayVertex::V_function(std::vector<double>& f_vector)
{
  double sum_f  = 0;
  double sum_f2 = 0;
  double v      = 0.;

  for(size_t i = 1; i < f_vector.size(); ++i)
    {
      sum_f += f_vector[i];
      sum_f2 += pow(f_vector[i], 2.);
    }

  if((sum_f > 1.E-9) && (sum_f2 > 1.E-9))
    {
      v = k_factor * f_vector[0] + sum_f -
          (k_factor * pow(f_vector[0], 2.) + sum_f2) / (k_factor * f_vector[0] + sum_f);
    }

  return v;
}

void TDecayVertex::SpaceDiscretization(double& Xi, double& Xf, size_t& NstepsX, double& Yi, double& Yf,
                                         size_t& NstepsY, double& Zi, double& Zf, size_t& NstepsZ, size_t& border,
                                         std::vector<TVector3>& PosXYZ)
{
  // if border == 0 -> Borders included
  // if border == 1 -> Borders not included

  if((border != 0) && (border != 1))
    {
      border = 0;
      cout << "Problem with parameter border (must be 0 or 1) -> Solved by including borders\n";
    }

  PosXYZ = {};

  for(size_t i = 0 + border; i < NstepsX - border; ++i)
    {
      double PosX = Xi + i * (Xf - Xi) / static_cast<double>(NstepsX - 1);

      for(size_t j = 0 + border; j < NstepsY - border; ++j)
        {
          double PosY = Yi + j * (Yf - Yi) / static_cast<double>(NstepsY - 1);

          for(size_t k = 0 + border; k < NstepsZ - border; ++k)
            {
              double PosZ = Zi + k * (Zf - Zi) / static_cast<double>(NstepsZ - 1);

              TVector3 temp_PosXYZ(PosX, PosY, PosZ);
              PosXYZ.emplace_back(temp_PosXYZ);
            }
        }
    }
}

void TDecayVertex::TrackstoDecayVertex(std::vector<DecayTrackInfo>& FragmentTracks, std::vector<DecayTrackInfo>& PionTracks,
                                        TVector3& PrimVtxRecons, TVector3& DecayVertexRecons)
{
  std::vector<double> temp_f(FragmentTracks.size() + PionTracks.size(), 0.);

  double V    = 0.;
  double Vnew = 0.;

  double Xi            = PrimVtxRecons.X() - boxDistXY;
  double Xf            = PrimVtxRecons.X() + boxDistXY;
  double distanceStepX = (Xf - Xi) / static_cast<double>(NstepsdiscretXY - 1);

  double Yi            = PrimVtxRecons.Y() - boxDistXY;
  double Yf            = PrimVtxRecons.Y() + boxDistXY;
  double distanceStepY = (Yf - Yi) / static_cast<double>(NstepsdiscretXY - 1);

  double Zi            = Zo_target;
  double Zf            = Zo_minifibers;
  double distanceStepZ = (Zf - Zi) / static_cast<double>(NstepsdiscretZ - 1);

  size_t border = 0;
  std::vector<TVector3> PosXYZ{};
  SpaceDiscretization(Xi, Xf, NstepsdiscretXY, Yi, Yf, NstepsdiscretXY, Zi, Zf, NstepsdiscretZ, border, PosXYZ);

  border = 1;

  for(size_t k = 0; k < nTimesDiscretization; ++k)
    {
      if(k != 0)
        {
          Xi            = DecayVertexRecons.X() - distanceStepX;
          Xf            = DecayVertexRecons.X() + distanceStepX;
          distanceStepX = (Xf - Xi) / static_cast<double>(Nstepsdiscretbox - 1);

          Yi            = DecayVertexRecons.Y() - distanceStepY;
          Yf            = DecayVertexRecons.Y() + distanceStepY;
          distanceStepY = (Yf - Yi) / static_cast<double>(Nstepsdiscretbox - 1);

          Zi            = DecayVertexRecons.Z() - distanceStepZ;
          Zf            = DecayVertexRecons.Z() + distanceStepZ;
          distanceStepZ = (Zf - Zi) / static_cast<double>(Nstepsdiscretbox - 1);

          SpaceDiscretization(Xi, Xf, Nstepsdiscretbox, Yi, Yf, Nstepsdiscretbox, Zi, Zf, Nstepsdiscretbox, border,
                              PosXYZ);
        }

      for(size_t i = 0; i < PosXYZ.size(); ++i)
        {
          TVector3 temp_PosXYZ = PosXYZ[i];
          
          for(size_t j = 0; j < FragmentTracks.size(); ++j)
            {
              temp_f[j] = f_function(FragmentTracks[j], temp_PosXYZ);
            }

          for(size_t j = 0; j < PionTracks.size(); ++j)
            {
              temp_f[FragmentTracks.size()+j] = f_function(PionTracks[j], temp_PosXYZ);
            }

          Vnew = V_function(temp_f);

          if(Vnew > V)
            {
              V = Vnew;

              DecayVertexRecons.SetXYZ(temp_PosXYZ.X(), temp_PosXYZ.Y(), temp_PosXYZ.Z());
            }
        }
    }
}

void TDecayVertex::Dist_DecayTrackPrimVtx(DecayTrackInfo& PionTrack, TVector3& PrimVtxRecons, double& distance)
{
  TVector3 u(PionTrack.Hit_MomEnergy.Px(), PionTrack.Hit_MomEnergy.Py(), PionTrack.Hit_MomEnergy.Pz());
  TVector3 PA = PrimVtxRecons - PionTrack.Hit_Pos;

  distance = (PA.Cross(u)).Mag()/u.Mag();
}

