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
      LocalHisto.h_DecayVtxstats->Fill("N_FragmentTracks=0", 1.);
      //return ReturnRes::DecayVtxError; //Check in the future
      return ReturnRes::Fine;
    }

  else if(result_full == -2)
    {
      att._logger->debug("No real pion tracks for decay vertex");
      LocalHisto.h_DecayVtxstats->Fill("N_RealPionTracks=0", 1.);
      //return ReturnRes::DecayVtxError; //Check in the future
      return ReturnRes::Fine;
    }

  else if(result_full == -3)
    {
      att._logger->debug("No cut pion tracks for decay vertex");
      LocalHisto.h_DecayVtxstats->Fill("N_CutPionTracks=0", 1.);
      //return ReturnRes::DecayVtxError; //Check in the future
      return ReturnRes::Fine;
    }

  else if(result_full == -4)
    {
      att._logger->debug("No pion tracks reconstructed for decay vertex");
      LocalHisto.h_DecayVtxstats->Fill("N_PionTracks=0", 1.);
      //return ReturnRes::DecayVtxError; //Check in the future
      return ReturnRes::Fine;
    }

  else if(result_full == -5)
    {
      att._logger->debug("No Si_MotherTrack reconstructed for decay vertex");
      LocalHisto.h_DecayVtxstats->Fill("N_Si_MotherTracks=0", 1.);
      //return ReturnRes::DecayVtxError; //Check in the future
      return ReturnRes::Fine;
    }

  else if(result_full == -6)
    {
      att._logger->debug("No MotherTracks reconstructed for decay vertex");
      LocalHisto.h_DecayVtxstats->Fill("N_MotherTracks=0", 1.);
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
  LocalHisto.h_Dist_FragmentTrackPrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_Dist_FragmentTrackPrimVtx);

  LocalHisto.h_Pt_pions = AnaHisto->CloneAndRegister(AnaHisto->h_Pt_pions);
  LocalHisto.h_Pz_pions = AnaHisto->CloneAndRegister(AnaHisto->h_Pz_pions);
  LocalHisto.h_Chi2ndf_pions = AnaHisto->CloneAndRegister(AnaHisto->h_Chi2ndf_pions);

  LocalHisto.h_Pt_cutpions = AnaHisto->CloneAndRegister(AnaHisto->h_Pt_cutpions);
  LocalHisto.h_Pz_cutpions = AnaHisto->CloneAndRegister(AnaHisto->h_Pz_cutpions);

  LocalHisto.h_Nrealpions = AnaHisto->CloneAndRegister(AnaHisto->h_Nrealpions);
  LocalHisto.h_Ncutpions = AnaHisto->CloneAndRegister(AnaHisto->h_Ncutpions);
  LocalHisto.h_Npions = AnaHisto->CloneAndRegister(AnaHisto->h_Npions);


  LocalHisto.h_Closedist_Distance = AnaHisto->CloneAndRegister(AnaHisto->h_Closedist_Distance);
  LocalHisto.h_Closedist_PosZ = AnaHisto->CloneAndRegister(AnaHisto->h_Closedist_PosZ);
  LocalHisto.h_Dist_DecayTrackPrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_Dist_DecayTrackPrimVtx);

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

  LocalHisto.h_DecayVertexDistance_KFPart = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistance_KFPart);
  LocalHisto.h_DecayVertexDistanceX_KFPart = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceX_KFPart);
  LocalHisto.h_DecayVertexDistanceY_KFPart = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceY_KFPart);
  LocalHisto.h_DecayVertexDistanceZ_KFPart = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceZ_KFPart);

  LocalHisto.h_DecayVertexDistance_KFPart_PrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistance_KFPart_PrimVtx);
  LocalHisto.h_DecayVertexDistanceX_KFPart_PrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceX_KFPart_PrimVtx);
  LocalHisto.h_DecayVertexDistanceY_KFPart_PrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceY_KFPart_PrimVtx);
  LocalHisto.h_DecayVertexDistanceZ_KFPart_PrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceZ_KFPart_PrimVtx);

  LocalHisto.h_DecayVertexDistance_KFPart_PrimVtx_Mass = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistance_KFPart_PrimVtx_Mass);
  LocalHisto.h_DecayVertexDistanceX_KFPart_PrimVtx_Mass = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceX_KFPart_PrimVtx_Mass);
  LocalHisto.h_DecayVertexDistanceY_KFPart_PrimVtx_Mass = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceY_KFPart_PrimVtx_Mass);
  LocalHisto.h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass);

  LocalHisto.h_DecayVertexcutDistance = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexcutDistance);
  LocalHisto.h_DecayVertexcutDistanceX = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexcutDistanceX);
  LocalHisto.h_DecayVertexcutDistanceY = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexcutDistanceY);
  LocalHisto.h_DecayVertexcutDistanceZ = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexcutDistanceZ);

/*
  LocalHisto.h_DecayVertexcutDistance_KFPart = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexcutDistance_KFPart);
  LocalHisto.h_DecayVertexcutDistanceX_KFPart = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexcutDistanceX_KFPart);
  LocalHisto.h_DecayVertexcutDistanceY_KFPart = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexcutDistanceY_KFPart);
  LocalHisto.h_DecayVertexcutDistanceZ_KFPart = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexcutDistanceZ_KFPart);
*/

  LocalHisto.h_DecayVertexcutDistance_KFPart_PrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexcutDistance_KFPart_PrimVtx);
  LocalHisto.h_DecayVertexcutDistanceX_KFPart_PrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexcutDistanceX_KFPart_PrimVtx);
  LocalHisto.h_DecayVertexcutDistanceY_KFPart_PrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexcutDistanceY_KFPart_PrimVtx);
  LocalHisto.h_DecayVertexcutDistanceZ_KFPart_PrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexcutDistanceZ_KFPart_PrimVtx);


  LocalHisto.h_DecayVertexPosZ_real = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexPosZ_real);
  LocalHisto.h_DecayVertexPosZ_vfunction = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexPosZ_vfunction);
  LocalHisto.h_DecayVertexPosZ_centroid = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexPosZ_centroid);
  LocalHisto.h_DecayVertexPosZ_KFPart = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexPosZ_KFPart);
  LocalHisto.h_DecayVertexPosZ_AllVfunc = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexPosZ_AllVfunc);
  LocalHisto.h_DecayVertexPosZ_AllCentroid = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexPosZ_AllCentroid);
  LocalHisto.h_DecayVertexPosZ_AllKFPart = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexPosZ_AllKFPart);

//HISTO FOR POSSIBLE CUTS ON MOTHER TRACK
  LocalHisto.h_N_MotherTracks = AnaHisto->CloneAndRegister(AnaHisto->h_N_MotherTracks);
  LocalHisto.h_Dist_DaughterTracks = AnaHisto->CloneAndRegister(AnaHisto->h_Dist_DaughterTracks);
  LocalHisto.h_Angle_MotherFragment = AnaHisto->CloneAndRegister(AnaHisto->h_Angle_MotherFragment);
  LocalHisto.h_Angle_MotherPion = AnaHisto->CloneAndRegister(AnaHisto->h_Angle_MotherPion);
  LocalHisto.h_Chi2ndf_MotherTracks = AnaHisto->CloneAndRegister(AnaHisto->h_Chi2ndf_MotherTracks);
  LocalHisto.h_Dist_MotherTrackPrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_Dist_MotherTrackPrimVtx);
  LocalHisto.h_Theta_MotherTrackPrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_Theta_MotherTrackPrimVtx);
  LocalHisto.h_DecayVertexPosZ_KFPart_PrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexPosZ_KFPart_PrimVtx);
  LocalHisto.h_Hyp_ArmenterosPodolanski = AnaHisto->CloneAndRegister(AnaHisto->h_Hyp_ArmenterosPodolanski);
  LocalHisto.h_Hyp_CutArmenterosPodolanski = AnaHisto->CloneAndRegister(AnaHisto->h_Hyp_CutArmenterosPodolanski);
  
  LocalHisto.h_HypInvariantMass = AnaHisto->CloneAndRegister(AnaHisto->h_HypInvariantMass);
  LocalHisto.h_HypErrorInvariantMass = AnaHisto->CloneAndRegister(AnaHisto->h_HypErrorInvariantMass);

  LocalHisto.h_Hyp_RealLifeTime = AnaHisto->CloneAndRegister(AnaHisto->h_Hyp_RealLifeTime);
  LocalHisto.h_HypLifeTime_PrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_HypLifeTime_PrimVtx);
  LocalHisto.h_HypErrorLifeTime_PrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_HypErrorLifeTime_PrimVtx);
  LocalHisto.h_HypcutLifeTime_PrimVtx = AnaHisto->CloneAndRegister(AnaHisto->h_HypcutLifeTime_PrimVtx);

  LocalHisto.h_HypInvariantMassCheck = AnaHisto->CloneAndRegister(AnaHisto->h_HypInvariantMassCheck);
  LocalHisto.h_HypInvariantErrorMassCheck = AnaHisto->CloneAndRegister(AnaHisto->h_HypInvariantErrorMassCheck);

  LocalHisto.h_HypInvariantMass_LorentzVect = AnaHisto->CloneAndRegister(AnaHisto->h_HypInvariantMass_LorentzVect);
  LocalHisto.h_HypInvariantMass_CutLorentzVect = AnaHisto->CloneAndRegister(AnaHisto->h_HypInvariantMass_CutLorentzVect);

  LocalHisto.h_EffPosZ_real = AnaHisto->CloneAndRegister(AnaHisto->h_EffPosZ_real);
  LocalHisto.h_EffPosZ_preKF = AnaHisto->CloneAndRegister(AnaHisto->h_EffPosZ_preKF);
  LocalHisto.h_EffPosZ_postKF = AnaHisto->CloneAndRegister(AnaHisto->h_EffPosZ_postKF);
  LocalHisto.h_EffPosZ_preKFPart = AnaHisto->CloneAndRegister(AnaHisto->h_EffPosZ_preKFPart);
  LocalHisto.h_EffPosZ_postKFPart = AnaHisto->CloneAndRegister(AnaHisto->h_EffPosZ_postKFPart);

/*
  LocalHisto.h_N_Si_MotherTracks = AnaHisto->CloneAndRegister(AnaHisto->h_N_Si_MotherTracks);

  LocalHisto.h_DecayVertexDistance_AllVfunc = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistance_AllVfunc);
  LocalHisto.h_DecayVertexDistanceX_AllVfunc = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceX_AllVfunc);
  LocalHisto.h_DecayVertexDistanceY_AllVfunc = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceY_AllVfunc);
  LocalHisto.h_DecayVertexDistanceZ_AllVfunc = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceZ_AllVfunc);

  LocalHisto.h_DecayVertexDistance_AllCentroid = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistance_AllCentroid);
  LocalHisto.h_DecayVertexDistanceX_AllCentroid = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceX_AllCentroid);
  LocalHisto.h_DecayVertexDistanceY_AllCentroid = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceY_AllCentroid);
  LocalHisto.h_DecayVertexDistanceZ_AllCentroid = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceZ_AllCentroid);

  LocalHisto.h_DecayVertexDistance_AllKFPart = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistance_AllKFPart);
  LocalHisto.h_DecayVertexDistanceX_AllKFPart = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceX_AllKFPart);
  LocalHisto.h_DecayVertexDistanceY_AllKFPart = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceY_AllKFPart);
  LocalHisto.h_DecayVertexDistanceZ_AllKFPart = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVertexDistanceZ_AllKFPart);
*/
  LocalHisto.h_DecayVtxstats = AnaHisto->CloneAndRegister(AnaHisto->h_DecayVtxstats);
}

int TDecayVertex::FinderDecayVertex(FullRecoEvent& RecoEvent)
{

  double DecayVertex_real_X = RecoEvent.DecayVertex[0];
  double DecayVertex_real_Y = RecoEvent.DecayVertex[1];
  double DecayVertex_real_Z = RecoEvent.DecayVertex[2];

  LocalHisto.h_DecayVertexPosZ_real->Fill(DecayVertex_real_Z, 1.);
  LocalHisto.h_EffPosZ_real->Fill(DecayVertex_real_Z, 1.);
  LocalHisto.h_Hyp_RealLifeTime->Fill(RecoEvent.Hyp_LifeTime, 1.);

  StudyCaseSelector(att.StudyCase, Fragment_pdg);

  //Primary vertex KFParticle initialization
  KFParticleSIMD KFPart_PrimVtx;
  KFPart_PrimaryVertex(RecoEvent.PrimVtxRecons, RecoEvent.CovMatrix_IP, KFPart_PrimVtx);
  const KFParticleSIMD* pointer_PrimVtx = &KFPart_PrimVtx;

  //Fragment tracks
  std::vector<KFParticle> FragmentTracks_All {};
  RealTracksFinder(RecoEvent.TrackDAFSim, Fragment_pdg, No_cutconditions, FragmentTracks_All);
  if(FragmentTracks_All.size() == 0)
    return -1;

  std::vector<KFParticle> FragmentTracks {};
  FragmentSelector(FragmentTracks_All, RecoEvent.PrimVtxRecons, FragmentTracks);
  if(FragmentTracks.size() == 0)
    return -1;

  double dist_FragmentTrackPrimVtx;
  double theta_FragmentTrackPrimVtx;

  TLorentzVector Fragment_LorentzVector;
  Fragment_LorentzVector.SetPxPyPzE(FragmentTracks[0].GetPx(),FragmentTracks[0].GetPy(),FragmentTracks[0].GetPz(),FragmentTracks[0].GetE());

  for(size_t i = 0; i < FragmentTracks.size(); ++i)
    {
      LocalHisto.h_Pt_fragments->Fill(FragmentTracks[i].GetPt(), 1.);
      LocalHisto.h_Pz_fragments->Fill(FragmentTracks[i].GetPz(), 1.);

      ThetaDist_TrackPrimVtx(FragmentTracks_All[i], RecoEvent.PrimVtxRecons, theta_FragmentTrackPrimVtx, dist_FragmentTrackPrimVtx);
      LocalHisto.h_Dist_FragmentTrackPrimVtx->Fill(dist_FragmentTrackPrimVtx, 1.);
    }

#ifdef REAL_PIONS_CHECK

  //Real pion tracks
  std::vector<KFParticle> RealPionTracks_All {};
  RealTracksFinder(RecoEvent.TrackDAFSim, pi_pdg, No_cutconditions, RealPionTracks_All);
  LocalHisto.h_Nrealpions->Fill(RealPionTracks_All.size(), 1.);
  if(RealPionTracks_All.size() == 0)
    return -2;

  std::vector<KFParticle> RealPionTracks {};
  PionSelector(RealPionTracks_All, RecoEvent.PrimVtxRecons, RealPionTracks);
  if(RealPionTracks.size() == 0)
    return -2;

/*
  double closedist_realdistance = 0.;
  TVector3 closedist_realpos;
  double dist_realDecayTrackPrimVtx;
  double theta_realDecayTrackPrimVtx;

  std::unordered_map<int, InfoInit>::iterator itr_real;
  bool ifDaughter_real = false;

  for(size_t i = 0; i < RealPionTracks.size(); ++i)
    {
      LocalHisto.h_Pt_realpions->Fill(RealPionTracks[i].GetPt(), 1.);
      LocalHisto.h_Pz_realpions->Fill(RealPionTracks[i].GetPz(), 1.);

      CloseDist(FragmentTracks[0], RealPionTracks[i], closedist_realdistance, closedist_realpos);
      LocalHisto.h_Closedist_realDistance->Fill(closedist_realdistance, 1.);
      LocalHisto.h_Closedist_realPosZ->Fill(closedist_realpos.Z(), 1.);

      ThetaDist_TrackPrimVtx(RealPionTracks[i], RecoEvent.PrimVtxRecons, theta_realDecayTrackPrimVtx, dist_realDecayTrackPrimVtx);
      LocalHisto.h_Dist_realDecayTrackPrimVtx->Fill(dist_realDecayTrackPrimVtx, "All", 1.);

      for(itr_real = RecoEvent.DaughtersTrackDAFInit.begin(); itr_real != RecoEvent.DaughtersTrackDAFInit.end(); ++itr_real)
        {
          if(itr_real->first == RealPionTracks[i].Id())
            ifDaughter_real = true;
        }

      if(ifDaughter_real)
        LocalHisto.h_Dist_realDecayTrackPrimVtx->Fill(dist_realDecayTrackPrimVtx, "Daughters", 1.);
      else
        LocalHisto.h_Dist_realDecayTrackPrimVtx->Fill(dist_realDecayTrackPrimVtx, "Primaries", 1.);
    }

  //Decay vertex reconstruction
  TVector3 DecayVertexReconsreal;

  TrackstoDecayVertex(FragmentTracks, RealPionTracks, RecoEvent.PrimVtxRecons, DecayVertexReconsreal);

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
*/
#endif

#ifdef CUT_PIONS_CHECK

 //Cut pion tracks
  std::vector<KFParticle> CutPionTracks_All {};
  RealTracksFinder(RecoEvent.TrackDAFSim, pi_pdg, Yes_cutconditions, CutPionTracks_All);

  LocalHisto.h_Ncutpions->Fill(CutPionTracks_All.size(), 1.);
  if(CutPionTracks_All.size() == 0)
    return -3;

  std::vector<KFParticle> CutPionTracks {};
  PionSelector(CutPionTracks_All, RecoEvent.PrimVtxRecons, CutPionTracks);
  if(CutPionTracks.size() == 0)
    return -3;

  double closedist_cutdistance = 0.;
  TVector3 closedist_cutpos;
  double theta_cutDecayTrackPrimVtx;
  double dist_cutDecayTrackPrimVtx;

  std::unordered_map<int, InfoInit>::iterator itr_cut;
  bool ifDaughter_cut = false;

  TLorentzVector CutPion_LorentzVector;

  for(size_t i = 0; i < CutPionTracks.size(); ++i)
    {
      LocalHisto.h_Pt_cutpions->Fill(CutPionTracks[i].GetPt(), 1.);
      LocalHisto.h_Pz_cutpions->Fill(CutPionTracks[i].GetPz(), 1.);

      CloseDist(FragmentTracks[0], CutPionTracks[i], closedist_cutdistance, closedist_cutpos);
      LocalHisto.h_Closedist_cutDistance->Fill(closedist_cutdistance, 1.);
      LocalHisto.h_Closedist_cutPosZ->Fill(closedist_cutpos.Z(), 1.);

      ThetaDist_TrackPrimVtx(CutPionTracks[i], RecoEvent.PrimVtxRecons, theta_cutDecayTrackPrimVtx, dist_cutDecayTrackPrimVtx);
      LocalHisto.h_Dist_cutDecayTrackPrimVtx->Fill(dist_cutDecayTrackPrimVtx, "All", 1.);

      for(itr_cut = RecoEvent.DaughtersTrackDAFInit.begin(); itr_cut != RecoEvent.DaughtersTrackDAFInit.end(); ++itr_cut)
        {
          if(itr_cut->first == CutPionTracks[i].Id())
            ifDaughter_cut = true;
        }

      if(ifDaughter_cut)
        {
          LocalHisto.h_Dist_cutDecayTrackPrimVtx->Fill(dist_cutDecayTrackPrimVtx, "Daughters", 1.);

          CutPion_LorentzVector.SetPxPyPzE(CutPionTracks[i].GetPx(),CutPionTracks[i].GetPy(),CutPionTracks[i].GetPz(),CutPionTracks[i].GetE());
        }
      else
        LocalHisto.h_Dist_cutDecayTrackPrimVtx->Fill(dist_cutDecayTrackPrimVtx, "Primaries", 1.);
    }

  //Decay vertex reconstruction
  TVector3 DecayVertexReconscut;

  TrackstoDecayVertex(FragmentTracks, CutPionTracks, RecoEvent.PrimVtxRecons, DecayVertexReconscut);

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


  std::vector<KFParticle> CutMotherTracks_PrimVtx_All;
  std::vector<std::tuple<size_t, size_t>> RefCutDaughtersTracks_PrimVtx_All;
  MotherTracksRecons(FragmentTracks, CutPionTracks, pointer_PrimVtx, CutMotherTracks_PrimVtx_All, RefCutDaughtersTracks_PrimVtx_All);

  std::vector<KFParticle> CutMotherTracks_PrimVtx;
  std::vector<std::tuple<size_t, size_t>> RefCutDaughtersTracks_PrimVtx;
  MotherSelector(CutMotherTracks_PrimVtx_All, RefCutDaughtersTracks_PrimVtx_All, FragmentTracks, CutPionTracks, RecoEvent.PrimVtxRecons,
                 CutMotherTracks_PrimVtx, RefCutDaughtersTracks_PrimVtx);

  for(size_t i = 0; i < CutMotherTracks_PrimVtx.size(); ++i)
    {
      double cutdistance_KFPart_PrimVtx  = sqrt(pow((DecayVertex_real_X - CutMotherTracks_PrimVtx[i].GetX()), 2.) +
                              pow((DecayVertex_real_Y - CutMotherTracks_PrimVtx[i].GetY()), 2.) +
                              pow((DecayVertex_real_Z - CutMotherTracks_PrimVtx[i].GetZ()), 2.));
      double cutdistanceX_KFPart_PrimVtx = DecayVertex_real_X - CutMotherTracks_PrimVtx[i].GetX();
      double cutdistanceY_KFPart_PrimVtx = DecayVertex_real_Y - CutMotherTracks_PrimVtx[i].GetY();
      double cutdistanceZ_KFPart_PrimVtx = DecayVertex_real_Z - CutMotherTracks_PrimVtx[i].GetZ();

      LocalHisto.h_DecayVertexcutDistance_KFPart_PrimVtx->Fill(cutdistance_KFPart_PrimVtx, 1.);
      LocalHisto.h_DecayVertexcutDistanceX_KFPart_PrimVtx->Fill(cutdistanceX_KFPart_PrimVtx, 1.);
      LocalHisto.h_DecayVertexcutDistanceY_KFPart_PrimVtx->Fill(cutdistanceY_KFPart_PrimVtx, 1.);
      LocalHisto.h_DecayVertexcutDistanceZ_KFPart_PrimVtx->Fill(cutdistanceZ_KFPart_PrimVtx, 1.);

      LocalHisto.h_EffPosZ_preKF->Fill(CutMotherTracks_PrimVtx[i].GetZ(), 1.);

      double cutTimeLife_PrimVtx = CutMotherTracks_PrimVtx[i].GetLifeTime() / c_light_speed_cmps; //in ps
      LocalHisto.h_HypcutLifeTime_PrimVtx->Fill(cutTimeLife_PrimVtx, 1.);

      size_t temp_id_fragment = get<0>(RefCutDaughtersTracks_PrimVtx[i]);
      size_t temp_id_pion = get<1>(RefCutDaughtersTracks_PrimVtx[i]);

      float Vertex[3] = {CutMotherTracks_PrimVtx[i].GetX(), CutMotherTracks_PrimVtx[i].GetY(), CutMotherTracks_PrimVtx[i].GetZ()};
      FragmentTracks[temp_id_fragment].TransportToPoint(Vertex);
      CutPionTracks[temp_id_pion].TransportToPoint(Vertex);

      float armenterosQtAlfa[2] = {0.};
      CutMotherTracks_PrimVtx[i].GetArmenterosPodolanski_FromMother(FragmentTracks[temp_id_fragment], CutPionTracks[temp_id_pion], armenterosQtAlfa);
      LocalHisto.h_Hyp_CutArmenterosPodolanski->Fill(armenterosQtAlfa[1], armenterosQtAlfa[0], 1.);
    }

#endif

  //Pion tracks
  std::vector<KFParticle> PionTracks_All {};
  PionTracksFinder(RecoEvent.DAF_results, PionTracks_All);
  LocalHisto.h_Npions->Fill(PionTracks_All.size(), 1.);
  if(PionTracks_All.size() == 0)
    return -4;

  std::vector<KFParticle> PionTracks {};
  PionSelector(PionTracks_All, RecoEvent.PrimVtxRecons, PionTracks);
  if(PionTracks.size() == 0)
    return -4;

  double closedist_distance = 0.;
  TVector3 temp_closedist_pos;

  double newclosest_distance = 1000.;
  TVector3 new_closedist_pos;

  std::vector<TVector3> vect_closedist_pos;
  double theta_DecayTrackPrimVtx;
  double dist_DecayTrackPrimVtx;

  std::unordered_map<int, InfoInit>::iterator itr_recons;
  bool ifDaughter_recons = false;

  TLorentzVector Pion_LorentzVector;

  for(size_t i = 0; i < PionTracks.size(); ++i)
    {
      LocalHisto.h_Pt_pions->Fill(sqrt(pow(PionTracks[i].GetPx(),2.)
                                        + pow(PionTracks[i].GetPy(),2.)), 1.);
      LocalHisto.h_Pz_pions->Fill(PionTracks[i].GetPz(), 1.);

      CloseDist(FragmentTracks[0], PionTracks[i], closedist_distance, temp_closedist_pos);
      vect_closedist_pos.emplace_back(temp_closedist_pos);

      LocalHisto.h_Closedist_Distance->Fill(closedist_distance, 1.);
      LocalHisto.h_Closedist_PosZ->Fill(temp_closedist_pos.Z(), 1.);

      ThetaDist_TrackPrimVtx(PionTracks[i], RecoEvent.PrimVtxRecons, theta_DecayTrackPrimVtx, dist_DecayTrackPrimVtx);
      LocalHisto.h_Dist_DecayTrackPrimVtx->Fill(dist_DecayTrackPrimVtx, "All", 1.);

      for(itr_recons = RecoEvent.DaughtersTrackDAFInit.begin(); itr_recons != RecoEvent.DaughtersTrackDAFInit.end(); ++itr_recons)
        {
          if(itr_recons->first == PionTracks[i].Id())
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

          Pion_LorentzVector.SetPxPyPzE(PionTracks[i].GetPx(),PionTracks[i].GetPy(),PionTracks[i].GetPz(),PionTracks[i].GetE());
        }
    }

    TLorentzVector Mother_LorentzVector;
    Mother_LorentzVector = Fragment_LorentzVector + Pion_LorentzVector;
    LocalHisto.h_HypInvariantMass_LorentzVect->Fill(Mother_LorentzVector.M(), 1.);

    TLorentzVector CutMother_LorentzVector;
    CutMother_LorentzVector = Fragment_LorentzVector + CutPion_LorentzVector;   
    LocalHisto.h_HypInvariantMass_CutLorentzVect->Fill(CutMother_LorentzVector.M(), 1.);


#ifdef CENTROID_METHOD

  TVector3 closedist_pos;

  for(size_t i = 0; i < vect_closedist_pos.size(); ++i)
    {
      closedist_pos += vect_closedist_pos[i];
    }
  closedist_pos *= (1./static_cast<double>(vect_closedist_pos.size()));

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

#endif

#ifdef VFUNCTION_METHOD

  //Decay vertex reconstruction
  TVector3 DecayVertexRecons;

  TrackstoDecayVertex(FragmentTracks, PionTracks, RecoEvent.PrimVtxRecons, DecayVertexRecons);
  RecoEvent.DecayVtxRecons.SetXYZ(DecayVertexRecons.X(), DecayVertexRecons.Y(), DecayVertexRecons.Z());

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
  std::vector<KFParticle> MotherTracks_All;
  std::vector<std::tuple<size_t, size_t>> RefDaughtersTracks_All;
  ifSet_ProductionVertex = 0;
  ifSet_MassConstraint = 0;
  MotherTracksRecons(FragmentTracks, PionTracks, nullptr, MotherTracks_All, RefDaughtersTracks_All);

  std::vector<KFParticle> MotherTracks;
  std::vector<std::tuple<size_t, size_t>> RefDaughtersTracks;
  MotherSelector(MotherTracks_All, RefDaughtersTracks_All, FragmentTracks, PionTracks, RecoEvent.PrimVtxRecons,
                 MotherTracks, RefDaughtersTracks);

  for(size_t i = 0; i < MotherTracks.size(); ++i)
    {
      double distance_KFPart  = sqrt(pow((DecayVertex_real_X - MotherTracks[i].GetX()), 2.) +
                              pow((DecayVertex_real_Y - MotherTracks[i].GetY()), 2.) +
                              pow((DecayVertex_real_Z - MotherTracks[i].GetZ()), 2.));
      double distanceX_KFPart = DecayVertex_real_X - MotherTracks[i].GetX();
      double distanceY_KFPart = DecayVertex_real_Y - MotherTracks[i].GetY();
      double distanceZ_KFPart = DecayVertex_real_Z - MotherTracks[i].GetZ();

      LocalHisto.h_DecayVertexDistance_KFPart->Fill(distance_KFPart, 1.);
      LocalHisto.h_DecayVertexDistanceX_KFPart->Fill(distanceX_KFPart, 1.);
      LocalHisto.h_DecayVertexDistanceY_KFPart->Fill(distanceY_KFPart, 1.);
      LocalHisto.h_DecayVertexDistanceZ_KFPart->Fill(distanceZ_KFPart, 1.);

      LocalHisto.h_DecayVertexPosZ_KFPart->Fill(MotherTracks[i].GetZ(), 1.);
    }

/*
  std::vector<KFParticle> CutMotherTracks_All;
  std::vector<std::tuple<size_t, size_t>> RefCutDaughtersTracks_All;
  MotherTracksRecons(FragmentTracks, CutPionTracks, nullptr, CutMotherTracks_All, RefCutDaughtersTracks_All);

  std::vector<KFParticle> CutMotherTracks;
  std::vector<std::tuple<size_t, size_t>> RefCutDaughtersTracks;
  MotherSelector(CutMotherTracks_All, RefCutDaughtersTracks_All, FragmentTracks, CutPionTracks, RecoEvent.PrimVtxRecons,
                 CutMotherTracks, RefCutDaughtersTracks);

  for(size_t i = 0; i < CutMotherTracks.size(); ++i)
    {
      double cutdistance_KFPart  = sqrt(pow((DecayVertex_real_X - CutMotherTracks[i].GetX()), 2.) +
                              pow((DecayVertex_real_Y - CutMotherTracks[i].GetY()), 2.) +
                              pow((DecayVertex_real_Z - CutMotherTracks[i].GetZ()), 2.));
      double cutdistanceX_KFPart = DecayVertex_real_X - CutMotherTracks[i].GetX();
      double cutdistanceY_KFPart = DecayVertex_real_Y - CutMotherTracks[i].GetY();
      double cutdistanceZ_KFPart = DecayVertex_real_Z - CutMotherTracks[i].GetZ();

      LocalHisto.h_DecayVertexcutDistance_KFPart->Fill(cutdistance_KFPart, 1.);
      LocalHisto.h_DecayVertexcutDistanceX_KFPart->Fill(cutdistanceX_KFPart, 1.);
      LocalHisto.h_DecayVertexcutDistanceY_KFPart->Fill(cutdistanceY_KFPart, 1.);
      LocalHisto.h_DecayVertexcutDistanceZ_KFPart->Fill(cutdistanceZ_KFPart, 1.);
    }
*/

  std::vector<KFParticle> MotherTracks_PrimVtx_All;
  std::vector<std::tuple<size_t, size_t>> RefDaughtersTracks_PrimVtx_All;
  ifSet_ProductionVertex = 1;
  ifSet_MassConstraint = 0;
  MotherTracksRecons(FragmentTracks, PionTracks, pointer_PrimVtx, MotherTracks_PrimVtx_All, RefDaughtersTracks_PrimVtx_All);

  std::vector<KFParticle> MotherTracks_PrimVtx;
  std::vector<std::tuple<size_t, size_t>> RefDaughtersTracks_PrimVtx;
  MotherSelector(MotherTracks_PrimVtx_All, RefDaughtersTracks_PrimVtx_All, FragmentTracks, PionTracks, RecoEvent.PrimVtxRecons,
                 MotherTracks_PrimVtx, RefDaughtersTracks_PrimVtx);


  for(size_t i = 0; i < MotherTracks_PrimVtx.size(); ++i)
    {
      double distance_KFPart_PrimVtx  = sqrt(pow((DecayVertex_real_X - MotherTracks_PrimVtx[i].GetX()), 2.) +
                              pow((DecayVertex_real_Y - MotherTracks_PrimVtx[i].GetY()), 2.) +
                              pow((DecayVertex_real_Z - MotherTracks_PrimVtx[i].GetZ()), 2.));
      double distanceX_KFPart_PrimVtx = DecayVertex_real_X - MotherTracks_PrimVtx[i].GetX();
      double distanceY_KFPart_PrimVtx = DecayVertex_real_Y - MotherTracks_PrimVtx[i].GetY();
      double distanceZ_KFPart_PrimVtx = DecayVertex_real_Z - MotherTracks_PrimVtx[i].GetZ();

      LocalHisto.h_DecayVertexDistance_KFPart_PrimVtx->Fill(distance_KFPart_PrimVtx, 1.);
      LocalHisto.h_DecayVertexDistanceX_KFPart_PrimVtx->Fill(distanceX_KFPart_PrimVtx, 1.);
      LocalHisto.h_DecayVertexDistanceY_KFPart_PrimVtx->Fill(distanceY_KFPart_PrimVtx, 1.);
      LocalHisto.h_DecayVertexDistanceZ_KFPart_PrimVtx->Fill(distanceZ_KFPart_PrimVtx, 1.);

      //Check the possible cuts
      size_t temp_id_fragment = get<0>(RefDaughtersTracks_PrimVtx[i]);
      size_t temp_id_pion = get<1>(RefDaughtersTracks_PrimVtx[i]);

      double Closedist_DaughterTracks;
      TVector3 Centroid_DaughtersTracks;
      CloseDist(FragmentTracks[temp_id_fragment], PionTracks[temp_id_pion], Closedist_DaughterTracks, Centroid_DaughtersTracks);
      LocalHisto.h_Dist_DaughterTracks->Fill(Closedist_DaughterTracks, MotherTracks_PrimVtx[i].GetMass(), 1.);

      double angle_MotherFragment = MotherTracks_PrimVtx[i].GetAngle(FragmentTracks[temp_id_fragment]) * 180. / M_PI;
      LocalHisto.h_Angle_MotherFragment->Fill(angle_MotherFragment, MotherTracks_PrimVtx[i].GetMass(), 1.);

      double angle_MotherPion = MotherTracks_PrimVtx[i].GetAngle(PionTracks[temp_id_pion]) * 180. / M_PI;
      LocalHisto.h_Angle_MotherPion->Fill(angle_MotherPion, MotherTracks_PrimVtx[i].GetMass(), 1.);

      double temp_chi2ndf = MotherTracks_PrimVtx[i].GetChi2() / static_cast<double>(MotherTracks_PrimVtx[i].GetNDF());
      LocalHisto.h_Chi2ndf_MotherTracks->Fill(temp_chi2ndf, MotherTracks_PrimVtx[i].GetMass(), 1.);

      double theta_MotherTrackPrimVtx;
      double dist_MotherTrackPrimVtx;
      ThetaDist_TrackPrimVtx(MotherTracks_PrimVtx[i], RecoEvent.PrimVtxRecons, theta_MotherTrackPrimVtx, dist_MotherTrackPrimVtx);
      LocalHisto.h_Dist_MotherTrackPrimVtx->Fill(dist_MotherTrackPrimVtx, MotherTracks_PrimVtx[i].GetMass(), 1.);
      LocalHisto.h_Theta_MotherTrackPrimVtx->Fill(theta_MotherTrackPrimVtx, MotherTracks_PrimVtx[i].GetMass(), 1.);

      LocalHisto.h_EffPosZ_postKFPart->Fill(MotherTracks_PrimVtx[i].GetZ(), 1.);
      LocalHisto.h_DecayVertexPosZ_KFPart_PrimVtx->Fill(MotherTracks_PrimVtx[i].GetZ(), MotherTracks_PrimVtx[i].GetMass(), 1.);

      LocalHisto.h_N_MotherTracks->Fill(MotherTracks_PrimVtx.size(), MotherTracks_PrimVtx[i].GetMass(), 1.);

      LocalHisto.h_HypInvariantMass->Fill(MotherTracks_PrimVtx[i].GetMass(), 1.);
      LocalHisto.h_HypErrorInvariantMass->Fill(MotherTracks_PrimVtx[i].GetErrMass(), 1.);

      float m_getmass;
      float error_getmass;
      int status_getmass = MotherTracks_PrimVtx[i].GetMass(m_getmass, error_getmass);

      LocalHisto.h_HypInvariantMassCheck->Fill(m_getmass, status_getmass, 1.);
      LocalHisto.h_HypInvariantErrorMassCheck->Fill(error_getmass, status_getmass, 1.);


      double TimeLife_PrimVtx = MotherTracks_PrimVtx[i].GetLifeTime() / c_light_speed_cmps; //in ps
      LocalHisto.h_HypLifeTime_PrimVtx->Fill(TimeLife_PrimVtx, 1.);

      double ErrorTimeLife_PrimVtx = MotherTracks_PrimVtx[i].GetErrLifeTime() / c_light_speed_cmps; //in ps
      LocalHisto.h_HypErrorLifeTime_PrimVtx->Fill(ErrorTimeLife_PrimVtx, 1.);

      float Vertex[3] = {MotherTracks_PrimVtx[i].GetX(), MotherTracks_PrimVtx[i].GetY(), MotherTracks_PrimVtx[i].GetZ()};
      FragmentTracks[temp_id_fragment].TransportToPoint(Vertex);
      PionTracks[temp_id_pion].TransportToPoint(Vertex);

      float armenterosQtAlfa[2] = {0.};
      MotherTracks_PrimVtx[i].GetArmenterosPodolanski_FromMother(FragmentTracks[temp_id_fragment], PionTracks[temp_id_pion], armenterosQtAlfa);
      LocalHisto.h_Hyp_ArmenterosPodolanski->Fill(armenterosQtAlfa[1], armenterosQtAlfa[0], 1.);
    }
  

  std::vector<KFParticle> MotherTracks_PrimVtx_Mass_All;
  std::vector<std::tuple<size_t, size_t>> RefDaughtersTracks_PrimVtx_Mass_All;
  ifSet_ProductionVertex = 1;
  ifSet_MassConstraint = 1;
  MotherTracksRecons(FragmentTracks, PionTracks, pointer_PrimVtx, MotherTracks_PrimVtx_Mass_All, RefDaughtersTracks_PrimVtx_Mass_All);

  std::vector<KFParticle> MotherTracks_PrimVtx_Mass;
  std::vector<std::tuple<size_t, size_t>> RefDaughtersTracks_PrimVtx_Mass;
  MotherSelector(MotherTracks_PrimVtx_Mass_All, RefDaughtersTracks_PrimVtx_Mass_All, FragmentTracks, PionTracks, RecoEvent.PrimVtxRecons,
                 MotherTracks_PrimVtx_Mass, RefDaughtersTracks_PrimVtx_Mass);


  for(size_t i = 0; i < MotherTracks_PrimVtx_Mass.size(); ++i)
    {
      double distance_KFPart_PrimVtx_Mass  = sqrt(pow((DecayVertex_real_X - MotherTracks_PrimVtx_Mass[i].GetX()), 2.) +
                              pow((DecayVertex_real_Y - MotherTracks_PrimVtx_Mass[i].GetY()), 2.) +
                              pow((DecayVertex_real_Z - MotherTracks_PrimVtx_Mass[i].GetZ()), 2.));
      double distanceX_KFPart_PrimVtx_Mass = DecayVertex_real_X - MotherTracks_PrimVtx_Mass[i].GetX();
      double distanceY_KFPart_PrimVtx_Mass = DecayVertex_real_Y - MotherTracks_PrimVtx_Mass[i].GetY();
      double distanceZ_KFPart_PrimVtx_Mass = DecayVertex_real_Z - MotherTracks_PrimVtx_Mass[i].GetZ();

      LocalHisto.h_DecayVertexDistance_KFPart_PrimVtx_Mass->Fill(distance_KFPart_PrimVtx_Mass, 1.);
      LocalHisto.h_DecayVertexDistanceX_KFPart_PrimVtx_Mass->Fill(distanceX_KFPart_PrimVtx_Mass, 1.);
      LocalHisto.h_DecayVertexDistanceY_KFPart_PrimVtx_Mass->Fill(distanceY_KFPart_PrimVtx_Mass, 1.);
      LocalHisto.h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass->Fill(distanceZ_KFPart_PrimVtx_Mass, 1.);

/*
      LocalHisto.h_DecayVertexPosZ_KFPart_PrimVtx->Fill(MotherTracks_PrimVtx[i].GetZ(), 1.);

      LocalHisto.h_N_MotherTracks->Fill(MotherTracks_PrimVtx.size(), MotherTracks_PrimVtx[i].GetMass(), 1.);

      ThetaDist_TrackPrimVtx(MotherTracks_PrimVtx[i], RecoEvent.PrimVtxRecons, theta_MotherTrackPrimVtx, dist_MotherTrackPrimVtx);
      LocalHisto.h_Dist_MotherTrackPrimVtx->Fill(dist_MotherTrackPrimVtx, MotherTracks_PrimVtx[i].GetMass(), 1.);
      LocalHisto.h_Theta_MotherTrackPrimVtx->Fill(theta_MotherTrackPrimVtx, MotherTracks_PrimVtx[i].GetMass(), 1.);

      LocalHisto.h_HypInvariantMass->Fill(MotherTracks_PrimVtx[i].GetMass(), 1.);
      LocalHisto.h_HypInvariantMassCheck->Fill(MotherTracks_PrimVtx[i].GetMass(), 1.);
      LocalHisto.h_HypErrorInvariantMass->Fill(MotherTracks_PrimVtx[i].GetErrMass(), 1.);

      double TimeLife_PrimVtx = MotherTracks_PrimVtx[i].GetLifeTime() / c_light_speed_cmps; //in ps
      LocalHisto.h_HypLifeTime_PrimVtx->Fill(TimeLife_PrimVtx, 1.);

      double ErrorTimeLife_PrimVtx = MotherTracks_PrimVtx[i].GetErrLifeTime() / c_light_speed_cmps; //in ps
      LocalHisto.h_HypErrorLifeTime_PrimVtx->Fill(ErrorTimeLife_PrimVtx, 1.);

      size_t temp_id_fragment = get<0>(RefDaughtersTracks_PrimVtx[i]);
      size_t temp_id_pion = get<1>(RefDaughtersTracks_PrimVtx[i]);

      float Vertex[3] = {MotherTracks_PrimVtx[i].GetX(), MotherTracks_PrimVtx[i].GetY(), MotherTracks_PrimVtx[i].GetZ()};
      FragmentTracks[temp_id_fragment].TransportToPoint(Vertex);
      PionTracks[temp_id_pion].TransportToPoint(Vertex);

      float armenterosQtAlfa[2] = {0.};
      MotherTracks_PrimVtx[i].GetArmenterosPodolanski_FromMother(FragmentTracks[temp_id_fragment], PionTracks[temp_id_pion], armenterosQtAlfa);
      LocalHisto.h_Hyp_ArmenterosPodolanski->Fill(armenterosQtAlfa[1], armenterosQtAlfa[0], 1.);
*/
    }



/*
  if(MotherTracks_PrimVtx.size() == 0)
    return -5;


  if(DecayVertexRecons.Z() > Z_plane_Si2 + Dist_to_Silicons)
    {
      LocalHisto.h_N_Si_MotherTracks->Fill(0.5, 1.);

      KFParticle Si_MotherTrack;
      MotherTrackSiliconHits(RecoEvent.PrimVtxRecons, DecayVertexRecons, RecoEvent.Hits_Si1, RecoEvent.Hits_Si2, MotherTracks_PrimVtx[0], Si_MotherTrack);
      
      if(Si_MotherTrack.GetP() < 1.e-4)
        return -6;
                  
      LocalHisto.h_N_Si_MotherTracks->Fill(2.5, 1.);

      std::vector<KFParticle> AllTracks;
      AllTracks.emplace_back(PionTracks[0]);
      AllTracks.emplace_back(FragmentTracks[0]);
      AllTracks.emplace_back(Si_MotherTrack);

      TVector3 AllVfunc_DecayVertexRecons;
      AllTrackstoDecayVertex_Vfunction(AllTracks, closedist_pos, AllVfunc_DecayVertexRecons);

      double AllVfunc_distance  = sqrt(pow((DecayVertex_real_X - AllVfunc_DecayVertexRecons.X()), 2.) +
                                         pow((DecayVertex_real_Y - AllVfunc_DecayVertexRecons.Y()), 2.) +
                                          pow((DecayVertex_real_Z - AllVfunc_DecayVertexRecons.Z()), 2.));
      double AllVfunc_distanceX = DecayVertex_real_X - AllVfunc_DecayVertexRecons.X();
      double AllVfunc_distanceY = DecayVertex_real_Y - AllVfunc_DecayVertexRecons.Y();
      double AllVfunc_distanceZ = DecayVertex_real_Z - AllVfunc_DecayVertexRecons.Z();

      LocalHisto.h_DecayVertexDistance_AllVfunc->Fill(AllVfunc_distance, 1.);
      LocalHisto.h_DecayVertexDistanceX_AllVfunc->Fill(AllVfunc_distanceX, 1.);
      LocalHisto.h_DecayVertexDistanceY_AllVfunc->Fill(AllVfunc_distanceY, 1.);
      LocalHisto.h_DecayVertexDistanceZ_AllVfunc->Fill(AllVfunc_distanceZ, 1.);

      LocalHisto.h_DecayVertexPosZ_AllVfunc->Fill(AllVfunc_DecayVertexRecons.Z(), 1.);



      TVector3 AllCentroid_DecayVertexRecons;
      AllTrackstoDecayVertex_Centroids(AllTracks, AllCentroid_DecayVertexRecons);

      double AllCentroid_distance  = sqrt(pow((DecayVertex_real_X - AllCentroid_DecayVertexRecons.X()), 2.) +
                                            pow((DecayVertex_real_Y - AllCentroid_DecayVertexRecons.Y()), 2.) +
                                             pow((DecayVertex_real_Z - AllCentroid_DecayVertexRecons.Z()), 2.));
      double AllCentroid_distanceX = DecayVertex_real_X - AllCentroid_DecayVertexRecons.X();
      double AllCentroid_distanceY = DecayVertex_real_Y - AllCentroid_DecayVertexRecons.Y();
      double AllCentroid_distanceZ = DecayVertex_real_Z - AllCentroid_DecayVertexRecons.Z();

      LocalHisto.h_DecayVertexDistance_AllCentroid->Fill(AllCentroid_distance, 1.);
      LocalHisto.h_DecayVertexDistanceX_AllCentroid->Fill(AllCentroid_distanceX, 1.);
      LocalHisto.h_DecayVertexDistanceY_AllCentroid->Fill(AllCentroid_distanceY, 1.);
      LocalHisto.h_DecayVertexDistanceZ_AllCentroid->Fill(AllCentroid_distanceZ, 1.);

      LocalHisto.h_DecayVertexPosZ_AllCentroid->Fill(AllCentroid_DecayVertexRecons.Z(), 1.);



      KFParticle KFPart_SecondaryVertex;

      KFParticle KFPart_Daughters [] = {PionTracks[0], FragmentTracks[0]};
      KFParticle *pointer_Daughters;
      pointer_Daughters = &KFPart_Daughters;
      int KFPart_n_daughters = 2;
      KFParticle KFPart_Parent [] = {Si_MotherTrack};
      KFParticle *pointer_Parent;
      pointer_Parent = &KFPart_Parent;
      double KFPart_mass = 0.; //Hyp_mass; for mass hypothesis

      KFPart_SecondaryVertex.Construct(&pointer_Daughters, KFPart_n_daughters, &pointer_Parent, KFPart_mass);

      TVector3 AllKFPart_DecayVertexRecons(KFPart_SecondaryVertex.GetX(), KFPart_SecondaryVertex.GetY(), KFPart_SecondaryVertex.GetY());

      double AllKFPart_distance  = sqrt(pow((DecayVertex_real_X - AllKFPart_DecayVertexRecons.X()), 2.) +
                                            pow((DecayVertex_real_Y - AllKFPart_DecayVertexRecons.Y()), 2.) +
                                             pow((DecayVertex_real_Z - AllKFPart_DecayVertexRecons.Z()), 2.));
      double AllKFPart_distanceX = DecayVertex_real_X - AllKFPart_DecayVertexRecons.X();
      double AllKFPart_distanceY = DecayVertex_real_Y - AllKFPart_DecayVertexRecons.Y();
      double AllKFPart_distanceZ = DecayVertex_real_Z - AllKFPart_DecayVertexRecons.Z();

      LocalHisto.h_DecayVertexDistance_AllKFPart->Fill(AllKFPart_distance, 1.);
      LocalHisto.h_DecayVertexDistanceX_AllKFPart->Fill(AllKFPart_distanceX, 1.);
      LocalHisto.h_DecayVertexDistanceY_AllKFPart->Fill(AllKFPart_distanceY, 1.);
      LocalHisto.h_DecayVertexDistanceZ_AllKFPart->Fill(AllKFPart_distanceZ, 1.);

      LocalHisto.h_DecayVertexPosZ_AllKFPart->Fill(AllKFPart_DecayVertexRecons.Z(), 1.);

    }
*/

  return 0;
}



void TDecayVertex::StudyCaseSelector(std::string StudyCase, int& Fragment_pdg)
{
  if(StudyCase.compare("H3L") == 0)
    {
      Fragment_pdg = He3_pdg;
    }
  else if(StudyCase.compare("H4L") == 0)
    {
      Fragment_pdg = He4_pdg;
    }
  else if(StudyCase.compare("nnL") == 0)
    {
      Fragment_pdg = deuteron_pdg;
    }
  else if(StudyCase.compare("lambda") == 0)
    {
      Fragment_pdg = proton_pdg;
    }
  else if(StudyCase.compare("background") == 0)
    {

    }

  Hyp_pdg = pid_fromName(StudyCase);
  Hyp_charge = TDatabasePDG::Instance()->GetParticle(Hyp_pdg)->Charge()/3.;
  Hyp_mass = TDatabasePDG::Instance()->GetParticle(Hyp_pdg)->Mass();
}


void TDecayVertex::RealTracksFinder(std::unordered_map<int, std::vector<std::vector<SimHit> > >& TrackDAFSim,
                                        int& pdgParticle, int& cutConditions,
                                        std::vector<KFParticle>& RealTracks)
{

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
          double temp_fP[] = {itr->second[iDetFirst][0].hitX, itr->second[iDetFirst][0].hitY, itr->second[iDetFirst][0].hitZ, // X, Y, Z
                              itr->second[iDetFirst][0].momX, itr->second[iDetFirst][0].momY, itr->second[iDetFirst][0].momZ}; // Px, Py, Pz

          double temp_fC[] = {1.e-8,
                                 0., 1.e-8,
                                 0.,    0., 1.e-8,
                                 0.,    0.,    0., 1.e-8,
                                 0.,    0.,    0.,    0., 1.e-8,
                                 0.,    0.,    0.,    0.,    0., 1.e-8}; //Change!

          int temp_charge = TDatabasePDG::Instance()->GetParticle(pdgParticle)->Charge()/3.;
          double temp_mass = TDatabasePDG::Instance()->GetParticle(pdgParticle)->Mass();

          KFParticle temp_particle;
          temp_particle.Create(temp_fP, temp_fC, temp_charge, temp_mass);
          temp_particle.SetId(itr->first);

          temp_particle.SetNDF(5); //Change !
          temp_particle.SetChi2(5.); // Change !
          
          //for(int iF = 0; iF < 10; ++iF)
            //temp_particle.SetFieldCoeff(0., iF);
          temp_particle.SetField(0.);

          if(cutConditions == 0)
              RealTracks.emplace_back(temp_particle);
          else if((cutConditions == 1) && (nHits_MDC >= 6) && (nHits_MiniFiber >= 4) && ((nHits_PSCE != 0) || (nHits_PSBE != 0)))
              RealTracks.emplace_back(temp_particle);
        }
    }
}

void TDecayVertex::FragmentSelector(std::vector<KFParticle>& FragmentTracks_All, TVector3& PrimVtxRecons, std::vector<KFParticle>& FragmentTracks)
{
/*
  Possible cuts:
  - upper limit to chi2/ndf
  - lower limit to d with IP
*/
  if(FragmentTracks_All.size() == 0)
    std::cout << "Error: N_FragmentTracks_All=0 in FragmentSelector\n";

  double temp_chi2ndf;
  double temp_theta;
  double temp_dist;

  for(size_t i = 0; i < FragmentTracks_All.size(); ++i)
    {
      temp_chi2ndf = FragmentTracks_All[i].GetChi2() / static_cast<double>(FragmentTracks_All[i].GetNDF());
      if( (ifCut_MaxChi2ndf_FragmentTracks == 1) && (temp_chi2ndf > MaxChi2ndf_FragmentTracks) )
        continue;

      ThetaDist_TrackPrimVtx(FragmentTracks_All[i], PrimVtxRecons, temp_theta, temp_dist);
      if( (ifCut_MinDist_FragmentTracksPrimVtx == 1) && (temp_dist < MinDist_FragmentTracksPrimVtx) )
        continue;

      FragmentTracks.emplace_back(FragmentTracks_All[i]);
    }
}



void TDecayVertex::PionTracksFinder(std::unordered_map<int, ResSolDAF>& DAF_results,
                                    std::vector<KFParticle>& PionTracks)
{
  std::unordered_map<int, ResSolDAF>::iterator itr;
  for(itr = DAF_results.begin(); itr != DAF_results.end(); ++itr)
    {
      if((itr->second.charge == -1) && (itr->second.Ncentral > att.KF_NbCentralCut))
        {
          LocalHisto.h_Chi2ndf_pions->Fill(itr->second.chi2 / itr->second.ndf, 1.);

          double temp_fP[] = {itr->second.posX, itr->second.posY, itr->second.posZ, // X, Y, Z
                              itr->second.momX, itr->second.momY, itr->second.momZ}; // Px, Py, Pz

          double temp_fC[] = {itr->second.cov_matrix[0][0],
                              itr->second.cov_matrix[1][0], itr->second.cov_matrix[1][1],
                              itr->second.cov_matrix[2][0], itr->second.cov_matrix[2][1], itr->second.cov_matrix[2][2],
                              itr->second.cov_matrix[3][0], itr->second.cov_matrix[3][1], itr->second.cov_matrix[3][2], itr->second.cov_matrix[3][3],
                              itr->second.cov_matrix[4][0], itr->second.cov_matrix[4][1], itr->second.cov_matrix[4][2], itr->second.cov_matrix[4][3], itr->second.cov_matrix[4][4],
                              itr->second.cov_matrix[5][0], itr->second.cov_matrix[5][1], itr->second.cov_matrix[5][2], itr->second.cov_matrix[5][3], itr->second.cov_matrix[5][4], itr->second.cov_matrix[5][5]};

          int temp_charge = TDatabasePDG::Instance()->GetParticle(pi_pdg)->Charge()/3.;
          double temp_mass = TDatabasePDG::Instance()->GetParticle(pi_pdg)->Mass();

          KFParticle temp_particle;

          temp_particle.Create(temp_fP, temp_fC, temp_charge, temp_mass);
          temp_particle.SetId(itr->first);
          temp_particle.SetField(0.);

          temp_particle.SetNDF(itr->second.ndf);
          temp_particle.SetChi2(itr->second.chi2);

          PionTracks.emplace_back(temp_particle);
        }
    }
}

void TDecayVertex::PionSelector(std::vector<KFParticle>& PionTracks_All, TVector3& PrimVtxRecons, std::vector<KFParticle>& PionTracks)
{
/*
  Possible cuts:
  - upper limit to chi2/ndf
  - lower limit to d with IP
*/
  if(PionTracks_All.size() == 0)
    std::cout << "Error: N_PionTracks_All=0 in PionSelector\n";

  double temp_chi2ndf;
  double temp_theta;
  double temp_dist;

  for(size_t i = 0; i < PionTracks_All.size(); ++i)
    {
      temp_chi2ndf = PionTracks_All[i].GetChi2() / static_cast<double>(PionTracks_All[i].GetNDF());
      if( (ifCut_MaxChi2ndf_PionTracks == 1) && (temp_chi2ndf > MaxChi2ndf_PionTracks) )
        continue;

      ThetaDist_TrackPrimVtx(PionTracks_All[i], PrimVtxRecons, temp_theta, temp_dist);
      if( (ifCut_MinDist_PionTracksPrimVtx == 1) && (temp_dist < MinDist_PionTracksPrimVtx) )
        continue;

      PionTracks.emplace_back(PionTracks_All[i]);
    }
}


void TDecayVertex::CloseDist(KFParticle& FragmentTrack, KFParticle& PionTrack, double& distance, TVector3& centroid)
{
  TVector3 n(FragmentTrack.GetPy() * PionTrack.GetPz() - FragmentTrack.GetPz() * PionTrack.GetPy(),
             FragmentTrack.GetPz() * PionTrack.GetPx() - FragmentTrack.GetPx() * PionTrack.GetPz(),
             FragmentTrack.GetPx() * PionTrack.GetPy() - FragmentTrack.GetPy() * PionTrack.GetPx());

  TVector3 n1(FragmentTrack.GetPy() * n.Z() - FragmentTrack.GetPz() * n.Y(),
              FragmentTrack.GetPz() * n.X() - FragmentTrack.GetPx() * n.Z(),
              FragmentTrack.GetPx() * n.Y() - FragmentTrack.GetPy() * n.X());

  TVector3 n2(PionTrack.GetPy() * n.Z() - PionTrack.GetPz() * n.Y(),
              PionTrack.GetPz() * n.X() - PionTrack.GetPx() * n.Z(),
              PionTrack.GetPx() * n.Y() - PionTrack.GetPy() * n.X());

  double FragmentFactor = ((PionTrack.GetX() - FragmentTrack.GetX()) * n2.X()
                            + (PionTrack.GetY() - FragmentTrack.GetY()) * n2.Y()
                            + (PionTrack.GetZ() - FragmentTrack.GetZ()) * n2.Z()) /
                                  (FragmentTrack.GetPx() * n2.X() + FragmentTrack.GetPy() * n2.Y()
                                   + FragmentTrack.GetPz() * n2.Z());

  double PionFactor = -((PionTrack.GetX() - FragmentTrack.GetX()) * n1.X()
                         + (PionTrack.GetY() - FragmentTrack.GetY()) * n1.Y()
                         + (PionTrack.GetZ() - FragmentTrack.GetZ()) * n1.Z()) /
                                (PionTrack.GetPx() * n1.X() + PionTrack.GetPy() * n1.Y()
                                 + PionTrack.GetPz() * n1.Z());
  
  TVector3 c1(FragmentTrack.GetX() + FragmentFactor * FragmentTrack.GetPx(),
              FragmentTrack.GetY() + FragmentFactor * FragmentTrack.GetPy(),
              FragmentTrack.GetZ() + FragmentFactor * FragmentTrack.GetPz());
  
  TVector3 c2(PionTrack.GetX() + PionFactor * PionTrack.GetPx(),
              PionTrack.GetY() + PionFactor * PionTrack.GetPy(),
              PionTrack.GetZ() + PionFactor * PionTrack.GetPz());;

  distance = (c2 - c1).Mag();
  centroid = (c1 + c2);
  centroid *= 0.5;
}

double TDecayVertex::f_function(KFParticle& DecayTrack, TVector3& PosXYZ)
{
  double slope_x     = DecayTrack.GetPx() / DecayTrack.GetPz();
  double intercept_x = DecayTrack.GetX() - slope_x * DecayTrack.GetZ();
  double slope_y     = DecayTrack.GetPy() / DecayTrack.GetPz();
  double intercept_y = DecayTrack.GetY() - slope_y * DecayTrack.GetZ();

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

void TDecayVertex::TrackstoDecayVertex(std::vector<KFParticle>& FragmentTracks, std::vector<KFParticle>& PionTracks,
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

void TDecayVertex::ThetaDist_TrackPrimVtx(KFParticle& Track, TVector3& PrimVtxRecons, double& theta, double& distance)
{
  TVector3 u(Track.GetPx(), Track.GetPy(), Track.GetPz());
  TVector3 P(Track.GetX(), Track.GetY(), Track.GetZ());
  TVector3 PA = P - PrimVtxRecons;

  theta = PA.Angle(u) * 180. / M_PI;
  distance = (PA.Cross(u)).Mag()/u.Mag();
}

void TDecayVertex::KFPart_PrimaryVertex(TVector3& PrimVtxRecons, std::array<double,6> Cov_PrimVtx, KFParticleSIMD& temp_PrimVtx)
{
  const float_v fP_PrimVtx [] = {PrimVtxRecons.X(), PrimVtxRecons.Y(), PrimVtxRecons.Z(), 0., 0., 0.};
  
  const float_v fC_PrimVtx [] = {Cov_PrimVtx[0],
                                 Cov_PrimVtx[1], Cov_PrimVtx[2],
                                 Cov_PrimVtx[3], Cov_PrimVtx[4], Cov_PrimVtx[5],
                                             0.,             0.,             0., 0.,
                                             0.,             0.,             0., 0., 0.,
                                             0.,             0.,             0., 0., 0., 0.};

  int_v fQ_PrimVtx = 0;
  float_v fMass_PrimVtx = 0.;

  temp_PrimVtx.Create(fP_PrimVtx, fC_PrimVtx, fQ_PrimVtx, fMass_PrimVtx);
  temp_PrimVtx.SetField(0.);
}

void TDecayVertex::MotherTracksRecons(std::vector<KFParticle>& FragmentTracks, std::vector<KFParticle>& PionTracks,
                                      const KFParticleSIMD* pointer_PrimVtx, std::vector<KFParticle>& MotherTracks,
                                      std::vector<std::tuple<size_t, size_t>>& RefDaughtersTracks)
{
  float_v ds[2] = {0.f,0.f};
  float_v dsdr[4][6];

  for(size_t i = 0; i < FragmentTracks.size(); ++i)
    {
      for(size_t j = 0; j < PionTracks.size(); ++j)
        {
          KFParticleSIMD particleSIMD1(FragmentTracks[i]);    // the same particle is copied to each SIMD element
          KFParticleSIMD particleSIMD2(PionTracks[j]);

          particleSIMD1.GetDStoParticle( particleSIMD2, ds, dsdr );
          particleSIMD1.TransportToDS(ds[0], dsdr[0]);
          particleSIMD2.TransportToDS(ds[1], dsdr[3]);

          const KFParticleSIMD* vDaughtersPointer[2] = {&particleSIMD1, &particleSIMD2};

          KFParticleSIMD mother;
          mother.SetConstructMethod(KFPart_fConstructMethod);
          mother.Construct(vDaughtersPointer, 2);

          if( (ifSet_ProductionVertex == 1) && (pointer_PrimVtx != nullptr) )
            mother.SetProductionVertex(*pointer_PrimVtx);

          if( ifSet_MassConstraint == 1)
            mother.SetMassConstraint(Hyp_mass);

          mother.TransportToDecayVertex();

          KFParticle temp_MotherTrack;
          mother.GetKFParticle(temp_MotherTrack, 0);
          temp_MotherTrack.SetField(0.);

          MotherTracks.emplace_back(temp_MotherTrack);
          RefDaughtersTracks.emplace_back(std::make_tuple(i,j));
        }
    }
}

void TDecayVertex::MotherSelector(std::vector<KFParticle>& MotherTracks_All, std::vector<std::tuple<size_t, size_t>>& RefDaughtersTracks_All,
                                  std::vector<KFParticle>& FragmentTracks, std::vector<KFParticle>& PionTracks, TVector3& PrimVtxRecons,
                                  std::vector<KFParticle>& MotherTracks, std::vector<std::tuple<size_t, size_t>>& RefDaughtersTracks)
{
/*
  Possible cuts:
  - upper limit d between daughters
  - upper limit to theta between mother and fragment
  - upper limit to theta between mother and pion
  - upper limit to chi2/ndf
  - upper limit to d with IP
  - upper limit to theta with IP-SV
  - upper and lower limit to Z_SV
  - cut on Armenteros-Podolanski plot
*/
  size_t temp_id_fragment;
  size_t temp_id_pion;

  double Closedist_DaughterTracks;
  TVector3 Centroid_DaughtersTracks;

  double angle_MotherFragment;
  double angle_MotherPion;
  double temp_chi2ndf;

  double Dist_MotherTrackPrimVtx;
  double Angle_MotherTrackPrimVtx;

  double PosZ_DecayVertex;

  float armenterosQtAlfa[2] = {0.};

  for(size_t i = 0; i < MotherTracks_All.size(); ++i)
    {
      std::tie(temp_id_fragment, temp_id_pion) = RefDaughtersTracks_All[i];

      CloseDist(FragmentTracks[temp_id_fragment], PionTracks[temp_id_pion], Closedist_DaughterTracks, Centroid_DaughtersTracks);
      if( (ifCut_MaxClosedist_DaughterTracks == 1) && (Closedist_DaughterTracks > MaxClosedist_DaughterTracks) )
        continue;

      angle_MotherFragment = MotherTracks_All[i].GetAngle(FragmentTracks[temp_id_fragment]) * 180. / M_PI;
      if( (ifCut_MaxAngle_MotherFragment == 1) && (angle_MotherFragment > MaxAngle_MotherFragment) )
        continue;

      angle_MotherPion = MotherTracks_All[i].GetAngle(PionTracks[temp_id_pion]) * 180. / M_PI;
      if( (ifCut_MaxAngle_MotherPion == 1) && (angle_MotherPion > MaxAngle_MotherPion) )
        continue;

      temp_chi2ndf = MotherTracks_All[i].GetChi2() / static_cast<double>(MotherTracks_All[i].GetNDF());
      if( (ifCut_MaxChi2ndf == 1) && (temp_chi2ndf < MaxChi2ndf) )
        continue;

      ThetaDist_TrackPrimVtx(MotherTracks_All[i], PrimVtxRecons, Dist_MotherTrackPrimVtx, Angle_MotherTrackPrimVtx);
      if( (ifCut_MaxDist_MotherTrackPrimVtx == 1) && (Dist_MotherTrackPrimVtx > MaxDist_MotherTrackPrimVtx) )
        continue;

      if( (ifCut_MaxAngle_MotherTrackPrimVtx == 1) && (Angle_MotherTrackPrimVtx > MaxAngle_MotherTrackPrimVtx))
        continue;

      PosZ_DecayVertex = MotherTracks_All[i].GetZ();
      if( (ifCut_MaxPosZ_DecayVertex == 1) && (PosZ_DecayVertex > MaxPosZ_DecayVertex) )
        continue;

      if( (ifCut_MinPosZ_DecayVertex == 1) && (PosZ_DecayVertex < MinPosZ_DecayVertex) )
        continue;

      MotherTracks_All[i].GetArmenterosPodolanski_FromMother(FragmentTracks[temp_id_fragment], PionTracks[temp_id_pion], armenterosQtAlfa);
      if ( ifCut_ArmenterosPodolanski ) //Change !
        continue;


      MotherTracks.emplace_back(MotherTracks_All[i]);
      RefDaughtersTracks.emplace_back(RefDaughtersTracks_All[i]);
    }
}


void TDecayVertex::MotherTrack_SiHit(TVector3& PrimVtxRecons, TVector3& DecayVtxRecons, double& Z_plane,
                                      std::vector<std::vector<double> >& Hits_Si, TVector3& Mother_Sihit)
{
  TVector3 IP_SV = DecayVtxRecons - PrimVtxRecons;

  double track_parameter = (Z_plane - PrimVtxRecons.Z()) / IP_SV.Z();
  double X_track = PrimVtxRecons.X() + track_parameter * IP_SV.X();
  double Y_track = PrimVtxRecons.Y() + track_parameter * IP_SV.Y();

  TVector3 Hit_track (X_track, Y_track, Z_plane);

  double new_dist = 50.;

  for(size_t i = 0; i < Hits_Si.size(); ++i)
    {
      if(Hits_Si[i][0] < Min_EnergyDeposition_Si)
        continue;

      TVector3 temp_hit (Hits_Si[i][1], Hits_Si[i][2], Hits_Si[i][3]);
      TVector3 temp_Vect_dist = Hit_track - temp_hit;
      double temp_dist = temp_Vect_dist.Mag();

      if(temp_dist < new_dist)
      {
        new_dist = temp_dist;
        Mother_Sihit = temp_hit;
      }
    }

  if(new_dist > Max_dist_Mother_SiHit)
    Mother_Sihit.SetXYZ(0.,0.,0.);
}

void TDecayVertex::MotherTrackSiliconHits(TVector3& PrimVtxRecons, TVector3& DecayVtxRecons, std::vector<std::vector<double> >& Hits_Si1,
                                            std::vector<std::vector<double> >& Hits_Si2, KFParticle& Mother, KFParticle& Si_MotherTrack)
{
  TVector3 Vect_Zero;

  TVector3 Mother_Hit_Si1;
  MotherTrack_SiHit(PrimVtxRecons, DecayVtxRecons, Z_plane_Si1, Hits_Si1, Mother_Hit_Si1);
  if(Mother_Hit_Si1 == Vect_Zero)
    return;

  TVector3 Mother_Hit_Si2;
  MotherTrack_SiHit(PrimVtxRecons, DecayVtxRecons, Z_plane_Si2, Hits_Si2, Mother_Hit_Si2);
  if(Mother_Hit_Si2 == Vect_Zero)
    return;

  TVector3 Si_MotherMom = Mother_Hit_Si2 - Mother_Hit_Si1;
  double sigma2         = pow(widthStrip_Si1, 2.) / 12.;

  double temp_fP[] = {Mother_Hit_Si2.X(), Mother_Hit_Si2.Y(), Mother_Hit_Si2.Z(), // X, Y, Z
                      Si_MotherMom.X(), Si_MotherMom.Y(), Si_MotherMom.Z()}; // Px, Py, Pz

  double temp_fC[] = {sigma2,
                          0., sigma2,
                          0.,     0., sigma2,
                          0.,     0.,     0., 1.e-6,
                          0.,     0.,     0.,    0., 1.e-6,
                          0.,     0.,     0.,    0.,    0., 1.e-6}; //Change momentum values!

  int temp_charge = Hyp_charge; //Change!
  const double temp_mass = Mother.GetMass(); //Change!

  Si_MotherTrack.Create(temp_fP, temp_fC, temp_charge, temp_mass);
  Si_MotherTrack.SetField(0.);
}

void TDecayVertex::AllTrackstoDecayVertex_Vfunction(std::vector<KFParticle>& AllTracks, TVector3& Old_DecayVertexRecons, TVector3& DecayVertexRecons)
{
  std::vector<double> temp_f(AllTracks.size(), 0.);

  double V    = 0.;
  double Vnew = 0.;

  double Xi            = Old_DecayVertexRecons.X() - boxXYZ;
  double Xf            = Old_DecayVertexRecons.X() + boxXYZ;
  double distanceStepX = (Xf - Xi) / static_cast<double>(NstepsdiscretboxXYZ - 1);

  double Yi            = Old_DecayVertexRecons.Y() - boxXYZ;
  double Yf            = Old_DecayVertexRecons.Y() + boxXYZ;
  double distanceStepY = (Yf - Yi) / static_cast<double>(NstepsdiscretboxXYZ - 1);

  double Zi            = Old_DecayVertexRecons.Z() - boxXYZ;
  double Zf            = Old_DecayVertexRecons.Z() + boxXYZ;
  double distanceStepZ = (Zf - Zi) / static_cast<double>(NstepsdiscretboxXYZ - 1);

  size_t border = 0;
  std::vector<TVector3> PosXYZ{};
  SpaceDiscretization(Xi, Xf, NstepsdiscretboxXYZ, Yi, Yf, NstepsdiscretboxXYZ, Zi, Zf, NstepsdiscretboxXYZ, border, PosXYZ);

  border = 1;

  for(size_t k = 0; k < nTimesBoxXYZ; ++k)
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

          SpaceDiscretization(Xi, Xf, NstepsdiscretboxXYZ, Yi, Yf, NstepsdiscretboxXYZ, Zi, Zf, NstepsdiscretboxXYZ, border,
                              PosXYZ);
        }

      for(size_t i = 0; i < PosXYZ.size(); ++i)
        {
          TVector3 temp_PosXYZ = PosXYZ[i];
          
          for(size_t j = 0; j < AllTracks.size(); ++j)
            {
              temp_f[j] = f_function(AllTracks[j], temp_PosXYZ);
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

void TDecayVertex::AllTrackstoDecayVertex_Centroids(std::vector<KFParticle>& AllTracks, TVector3& DecayVertexRecons)
{
  std::vector<TVector3> Vect_Centroids;

  double temp_distance;
  TVector3 temp_centroid;

  for(size_t i = 0; i < AllTracks.size()-1; ++i)
    {
      for(size_t j = i+1; j < AllTracks.size(); ++j)
        {
            CloseDist(AllTracks[i], AllTracks[j], temp_distance, temp_centroid);
            Vect_Centroids.emplace_back(temp_centroid);
        }
    }

  for(size_t i = 0; i < Vect_Centroids.size(); ++i)
    DecayVertexRecons += Vect_Centroids[i];

  DecayVertexRecons *= 1. / static_cast<double>(Vect_Centroids.size());
}
