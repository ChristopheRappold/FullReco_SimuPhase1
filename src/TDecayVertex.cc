#include "TDecayVertex.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"
#include "FullRecoEvent.hh"
#include "ReturnRes.hh"

#include <tuple>
#include <string>
#include <cmath>

#include "TLorentzVector.h"
#include "TVector3.h"

//#define DEBUG_DECAYVTX

#define REAL_PIONS_CHECK
//#define CUT_PIONS_CHECK
#define CENTROID_METHOD
#define VFUNCTION_METHOD 

using namespace std;
using namespace G4Sol;

template <class Out>
TDecayVertex<Out>::TDecayVertex(const THyphiAttributes& attribut, int pi_type) //pi_type=0 -> pi-  ; pi_type=1 -> pi+
    : TDataProcessInterface<Out>("DecayVertex_Reco_" + std::to_string(pi_type)), att(attribut), fieldWASA(att.Field)
  {
    Zo_target = -att.Target_Size; // in cm
    Zf_target =  att.Target_Size; // in cm
    Zo_minifibers = att.fiber_mft1_pos_z; // in cm

    pion_type = pi_type;
    if(pi_type==0) pi_pdg = piminus_pdg;
    if(pi_type==1) pi_pdg = piplus_pdg;

    StudyCaseSelector_Hyp(att.StudyCase);
    Hyp_charge = TDatabasePDG::Instance()->GetParticle(Hyp_pdg)->Charge()/3.;
    Hyp_mass = TDatabasePDG::Instance()->GetParticle(Hyp_pdg)->Mass();


  }

template <class Out>
TDecayVertex<Out>::~TDecayVertex() {}

template <class Out>
void TDecayVertex<Out>::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }

template <class Out>
ReturnRes::InfoM TDecayVertex<Out>::operator()(FullRecoEvent& RecoEvent, Out* OutTree)
{
  int result_finder = Exec(RecoEvent, OutTree);
  return SoftExit(result_finder);
}

template <class Out>
int TDecayVertex<Out>::Exec(FullRecoEvent& RecoEvent, Out* OutTree)
{
  int res_finder = FinderDecayVertex(RecoEvent); 

  for(auto i_Hyp : RecoEvent.Hyp_Vect)
  {
    if(std::isfinite(i_Hyp.Chi2ndf) == false)
      continue;
    if(std::isfinite(i_Hyp.Dist_Daughters) == false)
      continue;
    if(std::isfinite(i_Hyp.MomE_Fragment.Px()) == false)
      continue;

    THypernucleus* OutHyp = dynamic_cast<THypernucleus*>(OutTree->fHyp->ConstructedAt(OutTree->fHyp->GetEntries()));

    OutHyp->Pattern              = i_Hyp.Pattern;

    OutHyp->PDG                  = i_Hyp.PDG;
    OutHyp->N_Mother             = i_Hyp.N_Mother;
    OutHyp->Chi2ndf              = i_Hyp.Chi2ndf;
    OutHyp->NDF                  = i_Hyp.NDF;
    OutHyp->MomE                 = i_Hyp.MomE;
    OutHyp->PrimVtx              = i_Hyp.PrimVtx;
    OutHyp->DecayVtx             = i_Hyp.DecayVtx;
    OutHyp->Dist_RealReconsVtx   = i_Hyp.Dist_RealReconsVtx;
    OutHyp->Dist_MotherPrimVtx   = i_Hyp.Dist_MotherPrimVtx;
    OutHyp->Angle_MotherPrimVtx  = i_Hyp.Angle_MotherPrimVtx;
    OutHyp->InvMass              = i_Hyp.InvMass;
    OutHyp->ErrInvMass           = i_Hyp.ErrInvMass;
    OutHyp->ErrGetMass           = i_Hyp.ErrGetMass;
    OutHyp->LifeTime             = i_Hyp.LifeTime;
    OutHyp->ErrLifeTime          = i_Hyp.ErrLifeTime;
    OutHyp->ErrGetLifeTime       = i_Hyp.ErrGetLifeTime;
    OutHyp->Mother_IsFromHyp     = i_Hyp.Mother_IsFromHyp;

    OutHyp->PDG_Fragment         = i_Hyp.PDG_Fragment;
    OutHyp->Id_Fragment          = i_Hyp.Id_Fragment;
    OutHyp->MomE_Fragment        = i_Hyp.MomE_Fragment;
    OutHyp->Chi2ndf_Fragment     = i_Hyp.Chi2ndf_Fragment;
    OutHyp->NDF_Fragment         = i_Hyp.NDF_Fragment;
    OutHyp->Pvalue_Fragment      = i_Hyp.Pvalue_Fragment;
    OutHyp->Angle_MotherFragment = i_Hyp.Angle_MotherFragment;
    OutHyp->Fragment_IsFromHyp   = i_Hyp.Fragment_IsFromHyp;

    OutHyp->PDG_Pion             = i_Hyp.PDG_Pion;
    OutHyp->Id_Pion              = i_Hyp.Id_Pion;
    OutHyp->MomE_Pion            = i_Hyp.MomE_Pion;
    OutHyp->Chi2ndf_Pion         = i_Hyp.Chi2ndf_Pion;
    OutHyp->NDF_Pion             = i_Hyp.NDF_Pion;
    OutHyp->Pvalue_Pion          = i_Hyp.Pvalue_Pion;
    OutHyp->Angle_MotherPion     = i_Hyp.Angle_MotherPion;
    OutHyp->NHitsMDC_Pion        = i_Hyp.NHitsMDC_Pion;
    OutHyp->NHitsMinifiber_Pion  = i_Hyp.NHitsMinifiber_Pion;
    OutHyp->N_Pion               = i_Hyp.N_Pion;
    OutHyp->Pion_IsFromHyp       = i_Hyp.Pion_IsFromHyp;

    OutHyp->Dist_Daughters       = i_Hyp.Dist_Daughters;
    OutHyp->ArmPod_Qt            = i_Hyp.ArmPod_Qt;
    OutHyp->ArmPod_Alfa          = i_Hyp.ArmPod_Alfa;
  }

  OutTree->Nhyp = OutTree->fHyp->GetEntries();
  
  return res_finder;
}

template <class Out>
ReturnRes::InfoM TDecayVertex<Out>::SoftExit(int result_full) {
  
  if(result_full == -1)
    {
      att._logger->debug("No real/reconstructed fragment tracks for decay vertex");
      LocalHisto.h_DecayVtxstats[pion_type]->Fill("N_FragmentTracks=0", 1.);
      //return ReturnRes::DecayVtxError; //Check in the future
      return ReturnRes::Fine;
    }

  else if(result_full == -2)
    {
      att._logger->debug("No real pion tracks for decay vertex");
      LocalHisto.h_DecayVtxstats[pion_type]->Fill("N_RealPionTracks=0", 1.);
      //return ReturnRes::DecayVtxError; //Check in the future
      return ReturnRes::Fine;
    }

  else if(result_full == -3)
    {
      att._logger->debug("No cut pion tracks for decay vertex");
      LocalHisto.h_DecayVtxstats[pion_type]->Fill("N_CutPionTracks=0", 1.);
      //return ReturnRes::DecayVtxError; //Check in the future

      return ReturnRes::Fine;
    }

  else if(result_full == -4)
    {
      att._logger->debug("No pion tracks reconstructed for decay vertex");
      LocalHisto.h_DecayVtxstats[pion_type]->Fill("N_PionTracks=0", 1.);
      //return ReturnRes::DecayVtxError; //Check in the future
      return ReturnRes::Fine;
    }

  else if(result_full == -5)
    {
      att._logger->debug("No Si_MotherTrack reconstructed for decay vertex");
      LocalHisto.h_DecayVtxstats[pion_type]->Fill("N_Si_MotherTracks=0", 1.);
      //return ReturnRes::DecayVtxError; //Check in the future
      return ReturnRes::Fine;
    }

  else if(result_full == -6)
    {
      att._logger->debug("No MotherTracks reconstructed for decay vertex");
      LocalHisto.h_DecayVtxstats[pion_type]->Fill("N_MotherTracks=0", 1.);
      //return ReturnRes::DecayVtxError; //Check in the future
      return ReturnRes::Fine;
    }

  LocalHisto.h_DecayVtxstats[pion_type]->Fill("Fine", 1.);

  return ReturnRes::Fine; 
}



template <class Out>
void TDecayVertex<Out>::SelectHists()
{
  for(size_t i = 0; i < 2; ++i)
    {
      LocalHisto.h_P_fragments[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_P_fragments[i]);
      LocalHisto.h_Pt_fragments[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Pt_fragments[i]);
      LocalHisto.h_Pz_fragments[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Pz_fragments[i]);
      LocalHisto.h_Dist_FragmentTrackPrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Dist_FragmentTrackPrimVtx[i]);

      LocalHisto.h_P_pions[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_P_pions[i]);
      LocalHisto.h_Pt_pions[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Pt_pions[i]);
      LocalHisto.h_Pz_pions[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Pz_pions[i]);
      LocalHisto.h_Chi2ndf_pions[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Chi2ndf_pions[i]);

      LocalHisto.h_Pt_cutpions[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Pt_cutpions[i]);
      LocalHisto.h_Pz_cutpions[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Pz_cutpions[i]);

      LocalHisto.h_Nrealpions[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Nrealpions[i]);
      LocalHisto.h_Ncutpions[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Ncutpions[i]);
      LocalHisto.h_Npions[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Npions[i]);


      LocalHisto.h_Closedist_Distance[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Closedist_Distance[i]);
      LocalHisto.h_Closedist_PosZ[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Closedist_PosZ[i]);
      LocalHisto.h_Dist_DecayTrackPrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Dist_DecayTrackPrimVtx[i]);

      LocalHisto.h_Closedist_cutDistance[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Closedist_cutDistance[i]);
      LocalHisto.h_Closedist_cutPosZ[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Closedist_cutPosZ[i]);
      LocalHisto.h_Dist_cutDecayTrackPrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Dist_cutDecayTrackPrimVtx[i]);


      LocalHisto.h_DecayVertexDistance[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistance[i]);
      LocalHisto.h_DecayVertexDistanceX[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceX[i]);
      LocalHisto.h_DecayVertexDistanceY[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceY[i]);
      LocalHisto.h_DecayVertexDistanceZ[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceZ[i]);

      LocalHisto.h_DecayVertexDistance_centroid[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistance_centroid[i]);
      LocalHisto.h_DecayVertexDistanceX_centroid[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceX_centroid[i]);
      LocalHisto.h_DecayVertexDistanceY_centroid[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceY_centroid[i]);
      LocalHisto.h_DecayVertexDistanceZ_centroid[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceZ_centroid[i]);

      LocalHisto.h_DecayVertexDistance_KFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistance_KFPart[i]);
      LocalHisto.h_DecayVertexDistanceX_KFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceX_KFPart[i]);
      LocalHisto.h_DecayVertexDistanceY_KFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceY_KFPart[i]);
      LocalHisto.h_DecayVertexDistanceZ_KFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceZ_KFPart[i]);

      LocalHisto.h_DecayVertexDistance_KFPart_PrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistance_KFPart_PrimVtx[i]);
      LocalHisto.h_DecayVertexDistanceX_KFPart_PrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceX_KFPart_PrimVtx[i]);
      LocalHisto.h_DecayVertexDistanceY_KFPart_PrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceY_KFPart_PrimVtx[i]);
      LocalHisto.h_DecayVertexDistanceZ_KFPart_PrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceZ_KFPart_PrimVtx[i]);

      LocalHisto.h_DecayVertexDistance_KFPart_PrimVtx_Mass[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistance_KFPart_PrimVtx_Mass[i]);
      LocalHisto.h_DecayVertexDistanceX_KFPart_PrimVtx_Mass[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceX_KFPart_PrimVtx_Mass[i]);
      LocalHisto.h_DecayVertexDistanceY_KFPart_PrimVtx_Mass[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceY_KFPart_PrimVtx_Mass[i]);
      LocalHisto.h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass[i]);

      LocalHisto.h_DecayVertexcutDistance[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexcutDistance[i]);
      LocalHisto.h_DecayVertexcutDistanceX[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexcutDistanceX[i]);
      LocalHisto.h_DecayVertexcutDistanceY[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexcutDistanceY[i]);
      LocalHisto.h_DecayVertexcutDistanceZ[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexcutDistanceZ[i]);

    /*
      LocalHisto.h_DecayVertexcutDistance_KFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexcutDistance_KFPart[i]);
      LocalHisto.h_DecayVertexcutDistanceX_KFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexcutDistanceX_KFPart[i]);
      LocalHisto.h_DecayVertexcutDistanceY_KFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexcutDistanceY_KFPart[i]);
      LocalHisto.h_DecayVertexcutDistanceZ_KFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexcutDistanceZ_KFPart[i]);
    */

      LocalHisto.h_DecayVertexcutDistance_KFPart_PrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexcutDistance_KFPart_PrimVtx[i]);
      LocalHisto.h_DecayVertexcutDistanceX_KFPart_PrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexcutDistanceX_KFPart_PrimVtx[i]);
      LocalHisto.h_DecayVertexcutDistanceY_KFPart_PrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexcutDistanceY_KFPart_PrimVtx[i]);
      LocalHisto.h_DecayVertexcutDistanceZ_KFPart_PrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexcutDistanceZ_KFPart_PrimVtx[i]);


      LocalHisto.h_DecayVertexPosZ_real[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexPosZ_real[i]);
      LocalHisto.h_DecayVertexPosZ_vfunction[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexPosZ_vfunction[i]);
      LocalHisto.h_DecayVertexPosZ_centroid[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexPosZ_centroid[i]);
      LocalHisto.h_DecayVertexPosZ_KFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexPosZ_KFPart[i]);
      LocalHisto.h_DecayVertexPosZ_AllVfunc[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexPosZ_AllVfunc[i]);
      LocalHisto.h_DecayVertexPosZ_AllCentroid[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexPosZ_AllCentroid[i]);
      LocalHisto.h_DecayVertexPosZ_AllKFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexPosZ_AllKFPart[i]);

    //HISTO FOR POSSIBLE CUTS ON MOTHER TRACK
      LocalHisto.h_N_MotherTracks[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_N_MotherTracks[i]);
      LocalHisto.h_Dist_DaughterTracks[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Dist_DaughterTracks[i]);
      LocalHisto.h_Angle_MotherFragment[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Angle_MotherFragment[i]);
      LocalHisto.h_Angle_MotherPion[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Angle_MotherPion[i]);
      LocalHisto.h_Chi2ndf_MotherTracks[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Chi2ndf_MotherTracks[i]);
      LocalHisto.h_Dist_MotherTrackPrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Dist_MotherTrackPrimVtx[i]);
      LocalHisto.h_Theta_MotherTrackPrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Theta_MotherTrackPrimVtx[i]);
      LocalHisto.h_DecayVertexPosZ_KFPart_PrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexPosZ_KFPart_PrimVtx[i]);
      LocalHisto.h_DecayFragmentMomZ_KFPart_PrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayFragmentMomZ_KFPart_PrimVtx[i]);
      LocalHisto.h_DecayPionMomZ_KFPart_PrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayPionMomZ_KFPart_PrimVtx[i]);
      LocalHisto.h_Hyp_ArmenterosPodolanski[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Hyp_ArmenterosPodolanski[i]);
      LocalHisto.h_Hyp_CutArmenterosPodolanski[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Hyp_CutArmenterosPodolanski[i]);
      
      LocalHisto.h_HypInvariantMass[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HypInvariantMass[i]);
      LocalHisto.h_HypInvariantMass_Z05[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HypInvariantMass_Z05[i]);
      LocalHisto.h_HypInvariantMass_Z10[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HypInvariantMass_Z10[i]);
      LocalHisto.h_HypInvariantMass_Z15[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HypInvariantMass_Z15[i]);
      LocalHisto.h_HypInvariantMass_Z20[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HypInvariantMass_Z20[i]);
      LocalHisto.h_HypErrorInvariantMass[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HypErrorInvariantMass[i]);

      LocalHisto.h_Hyp_RealLifeTime[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Hyp_RealLifeTime[i]);
      LocalHisto.h_HypLifeTime_PrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HypLifeTime_PrimVtx[i]);
      LocalHisto.h_HypErrorLifeTime_PrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HypErrorLifeTime_PrimVtx[i]);
      LocalHisto.h_HypcutLifeTime_PrimVtx[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HypcutLifeTime_PrimVtx[i]);

      LocalHisto.h_HypInvariantMassCheck[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HypInvariantMassCheck[i]);
      LocalHisto.h_HypInvariantErrorMassCheck[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HypInvariantErrorMassCheck[i]);

      LocalHisto.h_HypInvariantMass_LorentzVect[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HypInvariantMass_LorentzVect[i]);
      LocalHisto.h_HypInvariantMass_CutLorentzVect[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HypInvariantMass_CutLorentzVect[i]);

      LocalHisto.h_EffPosZ_real[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EffPosZ_real[i]);
      LocalHisto.h_EffPosZ_preKF[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EffPosZ_preKF[i]);
      LocalHisto.h_EffPosZ_postKF[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EffPosZ_postKF[i]);
      LocalHisto.h_EffPosZ_preKFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EffPosZ_preKFPart[i]);
      LocalHisto.h_EffPosZ_postKFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EffPosZ_postKFPart[i]);

      LocalHisto.h_EffPosZPosR_real      [i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EffPosZPosR_real[i]);
      LocalHisto.h_EffPosZPosR_postKFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EffPosZPosR_postKFPart[i]);

    /*
      LocalHisto.h_N_SiHits_ReconsTracks[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_N_SiHits_ReconsTracks[i]);

      LocalHisto.h_N_Si_MotherTracks[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_N_Si_MotherTracks[i]);

      LocalHisto.h_DecayVertexDistance_AllVfunc[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistance_AllVfunc[i]);
      LocalHisto.h_DecayVertexDistanceX_AllVfunc[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceX_AllVfunc[i]);
      LocalHisto.h_DecayVertexDistanceY_AllVfunc[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceY_AllVfunc[i]);
      LocalHisto.h_DecayVertexDistanceZ_AllVfunc[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceZ_AllVfunc[i]);

      LocalHisto.h_DecayVertexDistance_AllCentroid[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistance_AllCentroid[i]);
      LocalHisto.h_DecayVertexDistanceX_AllCentroid[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceX_AllCentroid[i]);
      LocalHisto.h_DecayVertexDistanceY_AllCentroid[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceY_AllCentroid[i]);
      LocalHisto.h_DecayVertexDistanceZ_AllCentroid[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceZ_AllCentroid[i]);

      LocalHisto.h_DecayVertexDistance_AllKFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistance_AllKFPart[i]);
      LocalHisto.h_DecayVertexDistanceX_AllKFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceX_AllKFPart[i]);
      LocalHisto.h_DecayVertexDistanceY_AllKFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceY_AllKFPart[i]);
      LocalHisto.h_DecayVertexDistanceZ_AllKFPart[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVertexDistanceZ_AllKFPart[i]);
    */
      LocalHisto.h_DecayVtxstats[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayVtxstats[i]);
    }

}

template <class Out>
int TDecayVertex<Out>::FinderDecayVertex(FullRecoEvent& RecoEvent)
{
  TVector3 InteractionPoint_real(RecoEvent.InteractionPoint[0], RecoEvent.InteractionPoint[1], RecoEvent.InteractionPoint[2]);
  TVector3 DecayVertex_real(RecoEvent.DecayVertex[0], RecoEvent.DecayVertex[1], RecoEvent.DecayVertex[2]);

  LocalHisto.h_DecayVertexPosZ_real[pion_type]->Fill(DecayVertex_real.Z() - att.Target_PositionZ, 1.);
  LocalHisto.h_EffPosZ_real[pion_type]->Fill(DecayVertex_real.Z() - att.Target_PositionZ, 1.);
  double PosR_real = std::sqrt(std::pow(DecayVertex_real.X() - att.Target_PositionX, 2.)
                                + std::pow(DecayVertex_real.Y() - att.Target_PositionY, 2.));
  LocalHisto.h_EffPosZPosR_real[pion_type]->Fill(DecayVertex_real.Z() - att.Target_PositionZ, PosR_real, 1.);
  LocalHisto.h_Hyp_RealLifeTime[pion_type]->Fill(RecoEvent.Hyp_LifeTime, 1.);

  //Primary vertex KFParticle initialization
  KFParticleSIMD KFPart_PrimVtx_real;
  std::array<double,6> CovMatrix_IP_real = {0.01,
                                              0., 0.01,
                                              0.,   0., 0.1}; //temp value

  KFPart_PrimaryVertex(InteractionPoint_real, CovMatrix_IP_real, KFPart_PrimVtx_real);
  const KFParticleSIMD* pointer_PrimVtx_real = &KFPart_PrimVtx_real;

  KFParticleSIMD KFPart_PrimVtx;
  KFPart_PrimaryVertex(RecoEvent.PrimVtxRecons, RecoEvent.CovMatrix_IP, KFPart_PrimVtx);
  const KFParticleSIMD* pointer_PrimVtx = &KFPart_PrimVtx;


  //Fragment tracks

  if(recons_from_FRS_MDC == 2)
    FragmentMDCTracksFinder(RecoEvent.DAF_results, Fragment_pdg, RecoEvent.FragmentTracks);

  std::vector<KFParticle> FragmentTracks_All {};
  std::vector<KFFitInfo> Fragment_FitInfo {};
  
  for(size_t i = 0; i < RecoEvent.FragmentTracks.size(); ++i)
    {
      double temp_fP[] = {RecoEvent.FragmentTracks[i].GetPos().X(),  // x
                          RecoEvent.FragmentTracks[i].GetPos().Y(),  // y
                          RecoEvent.FragmentTracks[i].GetPos().Z(),  // z
                          RecoEvent.FragmentTracks[i].GetMom().X(),  // px
                          RecoEvent.FragmentTracks[i].GetMom().Y(),  // py
                          RecoEvent.FragmentTracks[i].GetMom().Z()}; // pz

      double temp_fC[] = {RecoEvent.FragmentTracks[i].GetCovMatrix()[0],   //x-x
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[1],   //x-y
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[2],   //y-y
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[3],   //x-z
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[4],   //y-z
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[5],   //z-z
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[6],   //x-px
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[7],   //y-px
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[8],   //z-px
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[9],   //px-px
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[10],  //x-py
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[11],  //y-py
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[12],  //z-py
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[13],  //px-py
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[14],  //py-py
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[15],  //x-pz
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[16],  //y-pz
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[17],  //z-pz
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[18],  //px-pz
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[19],  //py-pz
                          RecoEvent.FragmentTracks[i].GetCovMatrix()[20]}; //pz-pz

      int temp_charge = TDatabasePDG::Instance()->GetParticle(RecoEvent.FragmentTracks[i].GetPID())->Charge()/3.;
      double temp_mass = TDatabasePDG::Instance()->GetParticle(RecoEvent.FragmentTracks[i].GetPID())->Mass();

      KFParticle temp_particle;
      temp_particle.Create(temp_fP, temp_fC, temp_charge, temp_mass);
      temp_particle.SetId(RecoEvent.FragmentTracks[i].GetTrackID());
      temp_particle.SetPDG(RecoEvent.FragmentTracks[i].GetPID());

      temp_particle.SetNDF(1); //CHECK Change !
      temp_particle.SetChi2(RecoEvent.FragmentTracks[i].GetChi2NDF()); //CHECK Change !

      temp_particle.SetField(att.Field);

      KFFitInfo temp_FitInfo;
      temp_FitInfo.Pvalue = 0.; //CHECK Change!
      temp_FitInfo.NHitsMDC = 0;
      temp_FitInfo.NHitsMinifiber = 0;


      FragmentTracks_All.emplace_back(temp_particle);
      Fragment_FitInfo.emplace_back(temp_FitInfo);
    };

  if(FragmentTracks_All.size() == 0)
    return -1;
  
  std::vector<KFParticle> FragmentTracks {};
  FragmentSelector(FragmentTracks_All, RecoEvent.PrimVtxRecons, FragmentTracks);
  if(FragmentTracks.size() == 0)
    return -1;

  if(att.G4_simu)
    {
      std::unordered_map<int, InfoInit>::iterator itr_fragment;
      for(size_t i = 0; i < FragmentTracks.size(); ++i)
        {
          for(itr_fragment = RecoEvent.DaughtersTrackDAFInit.begin(); itr_fragment != RecoEvent.DaughtersTrackDAFInit.end(); ++itr_fragment)
            {
              if(itr_fragment->first == FragmentTracks[i].Id())
                ref_RealFragment = i;
            }
        }
    }
  else
    ref_RealFragment = 0; //CHECK Change!

  if(ifOnlyRealFragment && (ref_RealFragment == -1))
    return -1;

  double dist_FragmentTrackPrimVtx;
  double theta_FragmentTrackPrimVtx;

  TLorentzVector Fragment_LorentzVector;
  Fragment_LorentzVector.SetPxPyPzE(FragmentTracks[ref_RealFragment].GetPx(),FragmentTracks[ref_RealFragment].GetPy(),FragmentTracks[ref_RealFragment].GetPz(),FragmentTracks[ref_RealFragment].GetE());

  for(size_t i = 0; i < FragmentTracks.size(); ++i)
    {
      LocalHisto.h_P_fragments[pion_type]->Fill(FragmentTracks[i].GetP(), 1.);
      LocalHisto.h_Pt_fragments[pion_type]->Fill(FragmentTracks[i].GetPt(), 1.);
      LocalHisto.h_Pz_fragments[pion_type]->Fill(FragmentTracks[i].GetPz(), 1.);

      ThetaDist_TrackPrimVtx(FragmentTracks[i], RecoEvent.PrimVtxRecons, theta_FragmentTrackPrimVtx, dist_FragmentTrackPrimVtx);
      LocalHisto.h_Dist_FragmentTrackPrimVtx[pion_type]->Fill(dist_FragmentTrackPrimVtx, 1.);
    }

#ifdef REAL_PIONS_CHECK

  //Real pion tracks
  std::vector<KFParticle> RealPionTracks_All {};
  std::vector<KFFitInfo> RealPion_FitInfo {};
  RealTracksFinder(RecoEvent.TrackDAFSim, pi_pdg, No_cutconditions, RealPionTracks_All, RealPion_FitInfo);
  LocalHisto.h_Nrealpions[pion_type]->Fill(RealPionTracks_All.size(), 1.);
  if(RealPionTracks_All.size() == 0)
    return -2;

  std::vector<KFParticle> RealPionTracks {};
  PionSelector(RealPionTracks_All, RecoEvent.PrimVtxRecons, RealPionTracks, RealPion_FitInfo);
  if(RealPionTracks.size() == 0)
    return -2;

  double dist_RealDaughters;
  TVector3 centroid_RealDaughters;

  std::unordered_map<int, InfoInit>::iterator itr_real;

  for(size_t i = 0; i < RealPionTracks.size(); ++i)
    {
      Hyp temp_Hyp_Real;

      temp_Hyp_Real.Pattern = 1;

      temp_Hyp_Real.PDG = Hyp_pdg;
      temp_Hyp_Real.N_Mother = RealPionTracks.size();
      temp_Hyp_Real.MomE = RecoEvent.Mother_MomE;
      temp_Hyp_Real.PrimVtx = InteractionPoint_real;
      temp_Hyp_Real.DecayVtx.SetXYZ(RecoEvent.DecayVertex[0], RecoEvent.DecayVertex[1], RecoEvent.DecayVertex[2]);
      temp_Hyp_Real.Dist_RealReconsVtx.SetXYZ(0.,0.,0.);
      temp_Hyp_Real.Dist_MotherPrimVtx = 0.;
      temp_Hyp_Real.Angle_MotherPrimVtx = 0.;
      temp_Hyp_Real.InvMass = Hyp_mass;
      temp_Hyp_Real.ErrInvMass = 0.;
      temp_Hyp_Real.ErrGetMass = 0;
      temp_Hyp_Real.LifeTime = RecoEvent.Hyp_LifeTime;
      temp_Hyp_Real.ErrLifeTime = 0.;
      temp_Hyp_Real.ErrGetLifeTime = 0;

      TVector3 Mother_Mom(RecoEvent.Mother_MomE.Px(), RecoEvent.Mother_MomE.Py(), RecoEvent.Mother_MomE.Pz());

      temp_Hyp_Real.PDG_Fragment = FragmentTracks[ref_RealFragment].GetPDG();
      temp_Hyp_Real.Id_Fragment = FragmentTracks[ref_RealFragment].Id();
      temp_Hyp_Real.MomE_Fragment.SetPxPyPzE(FragmentTracks[ref_RealFragment].GetPx(), FragmentTracks[ref_RealFragment].GetPy(), FragmentTracks[ref_RealFragment].GetPz(), FragmentTracks[ref_RealFragment].GetE());
      temp_Hyp_Real.Chi2ndf_Fragment = FragmentTracks[ref_RealFragment].GetChi2() / static_cast<double>(FragmentTracks[ref_RealFragment].GetNDF());
      temp_Hyp_Real.NDF_Fragment = FragmentTracks[ref_RealFragment].GetNDF();
      temp_Hyp_Real.Pvalue_Fragment = Fragment_FitInfo[ref_RealFragment].Pvalue;
      TVector3 Fragment_Mom(FragmentTracks[ref_RealFragment].GetPx(), FragmentTracks[ref_RealFragment].GetPy(), FragmentTracks[ref_RealFragment].GetPz());
      temp_Hyp_Real.Angle_MotherFragment = Mother_Mom.Angle(Fragment_Mom) * 180. / M_PI;
      temp_Hyp_Real.Fragment_IsFromHyp = 1;

      temp_Hyp_Real.PDG_Pion = RealPionTracks[i].GetPDG();
      temp_Hyp_Real.Id_Pion = RealPionTracks[i].Id();
      temp_Hyp_Real.MomE_Pion.SetPxPyPzE(RealPionTracks[i].GetPx(), RealPionTracks[i].GetPy(), RealPionTracks[i].GetPz(), RealPionTracks[i].GetE());
      temp_Hyp_Real.Chi2ndf_Pion = RealPionTracks[i].GetChi2() / static_cast<double>(RealPionTracks[i].GetNDF());
      temp_Hyp_Real.NDF_Pion = RealPionTracks[i].GetNDF();
      temp_Hyp_Real.Pvalue_Pion = RealPion_FitInfo[i].Pvalue;
      TVector3 Pion_Mom(RealPionTracks[i].GetPx(), RealPionTracks[i].GetPy(), RealPionTracks[i].GetPz());
      temp_Hyp_Real.Angle_MotherPion = Mother_Mom.Angle(Pion_Mom) * 180. / M_PI;
      temp_Hyp_Real.NHitsMDC_Pion = RealPion_FitInfo[i].NHitsMDC;
      temp_Hyp_Real.NHitsMinifiber_Pion = RealPion_FitInfo[i].NHitsMinifiber;
      temp_Hyp_Real.N_Pion = RealPionTracks.size();
      temp_Hyp_Real.Pion_IsFromHyp = 0;

      for(itr_real = RecoEvent.DaughtersTrackDAFInit.begin(); itr_real != RecoEvent.DaughtersTrackDAFInit.end(); ++itr_real)
        {
          if(itr_real->first == RealPionTracks[i].Id())
            temp_Hyp_Real.Pion_IsFromHyp = 1;
        }

      temp_Hyp_Real.Mother_IsFromHyp = temp_Hyp_Real.Fragment_IsFromHyp * temp_Hyp_Real.Pion_IsFromHyp;

      CloseDist(FragmentTracks[ref_RealFragment], RealPionTracks[i], dist_RealDaughters, centroid_RealDaughters);
      temp_Hyp_Real.Dist_Daughters = dist_RealDaughters;

      float armenterosQtAlfa[2] = {0.};
      FragmentTracks[ref_RealFragment].GetArmenterosPodolanski(FragmentTracks[ref_RealFragment], RealPionTracks[i], armenterosQtAlfa);
      temp_Hyp_Real.ArmPod_Qt = armenterosQtAlfa[0];
      temp_Hyp_Real.ArmPod_Alfa = armenterosQtAlfa[1];

      RecoEvent.Hyp_Vect.emplace_back(temp_Hyp_Real);
    }

  std::vector<KFParticle> RealMotherTracks_All;
  std::vector<std::tuple<size_t, size_t>> RefRealDaughtersTracks_All;
  ifSet_ProductionVertex = true;
  ifSet_MassConstraint = false;
  MotherTracksRecons(FragmentTracks, RealPionTracks, pointer_PrimVtx_real, RealMotherTracks_All, RefRealDaughtersTracks_All);

  std::vector<KFParticle> RealMotherTracks;
  std::vector<std::tuple<size_t, size_t>> RefRealDaughtersTracks;
  MotherSelector(RealMotherTracks_All, RefRealDaughtersTracks_All, FragmentTracks, RealPionTracks, InteractionPoint_real,
                 RealMotherTracks, RefRealDaughtersTracks);

  double theta_realDecayTrackPrimVtx;
  double dist_realDecayTrackPrimVtx;
  float RealMother_InvMass;
  float RealMother_ErrInvMass;
  float RealMother_LifeTime;
  float RealMother_ErrLifeTime;

  for(size_t i = 0; i < RealMotherTracks.size(); ++i)
    {
      Hyp temp_Hyp_Real;

      temp_Hyp_Real.Pattern = 2;

/*
      if(std::isnan(RealMotherTracks[i].GetMass()))
      {
        std::cout << "0DcyPosX: " << RealMotherTracks[i].GetX()    << "\n";
        std::cout << "0DcyPosY: " << RealMotherTracks[i].GetY()    << "\n";
        std::cout << "0DcyPosZ: " << RealMotherTracks[i].GetZ()    << "\n";
        std::cout << "0FraPosX: " << FragmentTracks[get<0>(RefRealDaughtersTracks[i])].GetX() << "\n";
        std::cout << "0FraPosY: " << FragmentTracks[get<0>(RefRealDaughtersTracks[i])].GetY() << "\n";
        std::cout << "0FraPosZ: " << FragmentTracks[get<0>(RefRealDaughtersTracks[i])].GetZ() << "\n";
        std::cout << "0PioPosX: " << RealPionTracks[get<1>(RefRealDaughtersTracks[i])].GetX() << "\n";
        std::cout << "0PioPosY: " << RealPionTracks[get<1>(RefRealDaughtersTracks[i])].GetY() << "\n";
        std::cout << "0PioPosZ: " << RealPionTracks[get<1>(RefRealDaughtersTracks[i])].GetZ() << "\n";
      }
*/
      temp_Hyp_Real.PDG = RealMotherTracks[i].GetPDG();
      temp_Hyp_Real.N_Mother = RealMotherTracks.size();
      temp_Hyp_Real.Chi2ndf = RealMotherTracks[i].GetChi2() / static_cast<double>(RealMotherTracks[i].GetNDF());
      temp_Hyp_Real.NDF = RealMotherTracks[i].GetNDF();
      temp_Hyp_Real.MomE.SetPxPyPzE(RealMotherTracks[i].GetPx(),RealMotherTracks[i].GetPy(),RealMotherTracks[i].GetPz(),RealMotherTracks[i].GetE());
      RealMotherTracks[i].TransportToProductionVertex();
      temp_Hyp_Real.PrimVtx.SetXYZ(RealMotherTracks[i].GetX(),RealMotherTracks[i].GetY(),RealMotherTracks[i].GetZ());
      RealMotherTracks[i].TransportToDecayVertex();
      temp_Hyp_Real.DecayVtx.SetXYZ(RealMotherTracks[i].GetX(),RealMotherTracks[i].GetY(),RealMotherTracks[i].GetZ());
      temp_Hyp_Real.Dist_RealReconsVtx.SetXYZ(DecayVertex_real.X()-RealMotherTracks[i].GetX(),DecayVertex_real.Y()-RealMotherTracks[i].GetY(),DecayVertex_real.Z()-RealMotherTracks[i].GetZ());
      ThetaDist_TrackPrimVtx(RealMotherTracks[i], InteractionPoint_real, theta_realDecayTrackPrimVtx, dist_realDecayTrackPrimVtx);
      temp_Hyp_Real.Dist_MotherPrimVtx = dist_realDecayTrackPrimVtx;
      temp_Hyp_Real.Angle_MotherPrimVtx = theta_realDecayTrackPrimVtx;
      temp_Hyp_Real.ErrGetMass = RealMotherTracks[i].GetMass(RealMother_InvMass, RealMother_ErrInvMass);
      temp_Hyp_Real.InvMass = RealMother_InvMass;
      temp_Hyp_Real.ErrInvMass = RealMother_ErrInvMass;
      temp_Hyp_Real.ErrGetLifeTime = RealMotherTracks[i].GetLifeTime(RealMother_LifeTime, RealMother_ErrLifeTime);
      temp_Hyp_Real.LifeTime = RealMother_LifeTime / c_light_speed_cmps;
      temp_Hyp_Real.ErrLifeTime = RealMother_ErrLifeTime / c_light_speed_cmps;

      TVector3 Mother_Mom(RealMotherTracks[i].GetPx(), RealMotherTracks[i].GetPy(), RealMotherTracks[i].GetPz());

      size_t temp_id_fragment = get<0>(RefRealDaughtersTracks[i]);
      size_t temp_id_pion = get<1>(RefRealDaughtersTracks[i]);

      temp_Hyp_Real.PDG_Fragment = FragmentTracks[temp_id_fragment].GetPDG();
      temp_Hyp_Real.Id_Fragment = FragmentTracks[temp_id_fragment].Id();
      temp_Hyp_Real.MomE_Fragment.SetPxPyPzE(FragmentTracks[temp_id_fragment].GetPx(), FragmentTracks[temp_id_fragment].GetPy(), FragmentTracks[temp_id_fragment].GetPz(), FragmentTracks[temp_id_fragment].GetE());
      temp_Hyp_Real.Chi2ndf_Fragment = FragmentTracks[temp_id_fragment].GetChi2() / static_cast<double>(FragmentTracks[temp_id_fragment].GetNDF());      
      temp_Hyp_Real.NDF_Fragment = FragmentTracks[temp_id_fragment].GetNDF();
      temp_Hyp_Real.Pvalue_Fragment = Fragment_FitInfo[temp_id_fragment].Pvalue;
      TVector3 Fragment_Mom(FragmentTracks[temp_id_fragment].GetPx(), FragmentTracks[temp_id_fragment].GetPy(), FragmentTracks[temp_id_fragment].GetPz());
      temp_Hyp_Real.Angle_MotherFragment = Mother_Mom.Angle(Fragment_Mom) * 180. / M_PI;
      temp_Hyp_Real.Fragment_IsFromHyp = 0;
      if(FragmentTracks[temp_id_fragment].Id() == FragmentTracks[ref_RealFragment].Id())
        temp_Hyp_Real.Fragment_IsFromHyp = 1;


      temp_Hyp_Real.PDG_Pion = RealPionTracks[temp_id_pion].GetPDG();
      temp_Hyp_Real.Id_Pion = RealPionTracks[temp_id_pion].Id();
      temp_Hyp_Real.MomE_Pion.SetPxPyPzE(RealPionTracks[temp_id_pion].GetPx(), RealPionTracks[temp_id_pion].GetPy(), RealPionTracks[temp_id_pion].GetPz(), RealPionTracks[temp_id_pion].GetE());
      temp_Hyp_Real.Chi2ndf_Pion = RealPionTracks[temp_id_pion].GetChi2() / static_cast<double>(RealPionTracks[temp_id_pion].GetNDF());
      temp_Hyp_Real.NDF_Pion = RealPionTracks[temp_id_pion].GetNDF();
      temp_Hyp_Real.Pvalue_Pion = RealPion_FitInfo[temp_id_pion].Pvalue;
      TVector3 Pion_Mom(RealPionTracks[temp_id_pion].GetPx(), RealPionTracks[temp_id_pion].GetPy(), RealPionTracks[temp_id_pion].GetPz());
      temp_Hyp_Real.Angle_MotherPion = Mother_Mom.Angle(Pion_Mom) * 180. / M_PI;
      temp_Hyp_Real.NHitsMDC_Pion = RealPion_FitInfo[temp_id_pion].NHitsMDC;
      temp_Hyp_Real.NHitsMinifiber_Pion = RealPion_FitInfo[temp_id_pion].NHitsMinifiber;
      temp_Hyp_Real.N_Pion = RealPionTracks.size();
      temp_Hyp_Real.Pion_IsFromHyp = 0;

      for(itr_real = RecoEvent.DaughtersTrackDAFInit.begin(); itr_real != RecoEvent.DaughtersTrackDAFInit.end(); ++itr_real)
        {
          if(itr_real->first == RealPionTracks[temp_id_pion].Id())
            temp_Hyp_Real.Pion_IsFromHyp = 1;
        }

      temp_Hyp_Real.Mother_IsFromHyp = temp_Hyp_Real.Fragment_IsFromHyp * temp_Hyp_Real.Pion_IsFromHyp;

      float Vertex[3] = {RealMotherTracks[i].GetX(), RealMotherTracks[i].GetY(), RealMotherTracks[i].GetZ()};
      FragmentTracks[temp_id_fragment].TransportToPoint(Vertex);
      RealPionTracks[temp_id_pion].TransportToPoint(Vertex);

      CloseDist(FragmentTracks[temp_id_fragment], RealPionTracks[temp_id_pion], dist_RealDaughters, centroid_RealDaughters);
      temp_Hyp_Real.Dist_Daughters = dist_RealDaughters;

      float armenterosQtAlfa[2] = {0.};
      RealMotherTracks[i].GetArmenterosPodolanski_FromMother(FragmentTracks[temp_id_fragment], RealPionTracks[temp_id_pion], armenterosQtAlfa);
      temp_Hyp_Real.ArmPod_Qt = armenterosQtAlfa[0];
      temp_Hyp_Real.ArmPod_Alfa = armenterosQtAlfa[1];


/*
      if(std::isnan(RealMotherTracks[i].GetMass()))
      {
        std::cout << "1DcyPosX: " << RealMotherTracks[i].GetX()    << "\n";
        std::cout << "1DcyPosY: " << RealMotherTracks[i].GetY()    << "\n";
        std::cout << "1DcyPosZ: " << RealMotherTracks[i].GetZ()    << "\n";
        std::cout << "1FraPosX: " << FragmentTracks[temp_id_fragment].GetX() << "\n";
        std::cout << "1FraPosY: " << FragmentTracks[temp_id_fragment].GetY() << "\n";
        std::cout << "1FraPosZ: " << FragmentTracks[temp_id_fragment].GetZ() << "\n";
        std::cout << "1PioPosX: " << RealPionTracks[temp_id_pion].GetX() << "\n";
        std::cout << "1PioPosY: " << RealPionTracks[temp_id_pion].GetY() << "\n";
        std::cout << "1PioPosZ: " << RealPionTracks[temp_id_pion].GetZ() << "\n";
      }
*/




      RecoEvent.Hyp_Vect.emplace_back(temp_Hyp_Real);
    }

#endif

#ifdef CUT_PIONS_CHECK

 //Cut pion tracks
  std::vector<KFParticle> CutPionTracks_All {};
  std::vector<KFFitInfo> CutPion_FitInfo {};
  RealTracksFinder(RecoEvent.TrackDAFSim, pi_pdg, Yes_cutconditions, CutPionTracks_All, CutPion_FitInfo);

  LocalHisto.h_Ncutpions[pion_type]->Fill(CutPionTracks_All.size(), 1.);

  if(CutPionTracks_All.size() == 0)
    return -3;

  std::vector<KFParticle> CutPionTracks {};
  PionSelector(CutPionTracks_All, RecoEvent.PrimVtxRecons, CutPionTracks, CutPion_FitInfo);
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
      LocalHisto.h_Pt_cutpions[pion_type]->Fill(CutPionTracks[i].GetPt(), 1.);
      LocalHisto.h_Pz_cutpions[pion_type]->Fill(CutPionTracks[i].GetPz(), 1.);

      CloseDist(FragmentTracks[ref_RealFragment], CutPionTracks[i], closedist_cutdistance, closedist_cutpos);
      LocalHisto.h_Closedist_cutDistance[pion_type]->Fill(closedist_cutdistance, 1.);
      LocalHisto.h_Closedist_cutPosZ[pion_type]->Fill(closedist_cutpos.Z() - att.Target_PositionZ, 1.);

      ThetaDist_TrackPrimVtx(CutPionTracks[i], RecoEvent.PrimVtxRecons, theta_cutDecayTrackPrimVtx, dist_cutDecayTrackPrimVtx);
      LocalHisto.h_Dist_cutDecayTrackPrimVtx[pion_type]->Fill(dist_cutDecayTrackPrimVtx, "All", 1.);

      for(itr_cut = RecoEvent.DaughtersTrackDAFInit.begin(); itr_cut != RecoEvent.DaughtersTrackDAFInit.end(); ++itr_cut)
        {
          if(itr_cut->first == CutPionTracks[i].Id())
            ifDaughter_cut = true;
        }

      if(ifDaughter_cut)
        {
          LocalHisto.h_Dist_cutDecayTrackPrimVtx[pion_type]->Fill(dist_cutDecayTrackPrimVtx, "Daughters", 1.);

          CutPion_LorentzVector.SetPxPyPzE(CutPionTracks[i].GetPx(),CutPionTracks[i].GetPy(),CutPionTracks[i].GetPz(),CutPionTracks[i].GetE());
        }
      else
        LocalHisto.h_Dist_cutDecayTrackPrimVtx[pion_type]->Fill(dist_cutDecayTrackPrimVtx, "Primaries", 1.);
    }

  //Decay vertex reconstruction
  TVector3 DecayVertexReconscut;

  TrackstoDecayVertex(FragmentTracks, CutPionTracks, RecoEvent.PrimVtxRecons, DecayVertexReconscut);

  double cutdistance  = (DecayVertex_real - DecayVertexReconscut).Mag();
  double cutdistanceX = DecayVertex_real.X() - DecayVertexReconscut.X();
  double cutdistanceY = DecayVertex_real.Y() - DecayVertexReconscut.Y();
  double cutdistanceZ = DecayVertex_real.Z() - DecayVertexReconscut.Z();

  LocalHisto.h_DecayVertexcutDistance[pion_type]->Fill(cutdistance, 1.);
  LocalHisto.h_DecayVertexcutDistanceX[pion_type]->Fill(cutdistanceX, 1.);
  LocalHisto.h_DecayVertexcutDistanceY[pion_type]->Fill(cutdistanceY, 1.);
  LocalHisto.h_DecayVertexcutDistanceZ[pion_type]->Fill(cutdistanceZ, 1.);


  std::vector<KFParticle> CutMotherTracks_PrimVtx_All;
  std::vector<std::tuple<size_t, size_t>> RefCutDaughtersTracks_PrimVtx_All;
  ifSet_ProductionVertex = true;
  ifSet_MassConstraint = false;
  MotherTracksRecons(FragmentTracks, CutPionTracks, pointer_PrimVtx, CutMotherTracks_PrimVtx_All, RefCutDaughtersTracks_PrimVtx_All);

  std::vector<KFParticle> CutMotherTracks_PrimVtx;
  std::vector<std::tuple<size_t, size_t>> RefCutDaughtersTracks_PrimVtx;
  MotherSelector(CutMotherTracks_PrimVtx_All, RefCutDaughtersTracks_PrimVtx_All, FragmentTracks, CutPionTracks, RecoEvent.PrimVtxRecons,
                 CutMotherTracks_PrimVtx, RefCutDaughtersTracks_PrimVtx);

  float CutMother_InvMass;
  float CutMother_ErrInvMass;
  float CutMother_LifeTime;
  float CutMother_ErrLifeTime;

  double dist_CutDaughters;
  TVector3 centroid_CutDaughters;

  for(size_t i = 0; i < CutMotherTracks_PrimVtx.size(); ++i)
    {
      Hyp temp_Hyp_Cut;

      temp_Hyp_Cut.Pattern = 3;

      temp_Hyp_Cut.PDG = CutMotherTracks_PrimVtx[i].GetPDG();
      temp_Hyp_Cut.N_Mother = CutMotherTracks_PrimVtx.size();
      temp_Hyp_Cut.Chi2ndf = CutMotherTracks_PrimVtx[i].GetChi2() / static_cast<double>(CutMotherTracks_PrimVtx[i].GetNDF());
      temp_Hyp_Cut.NDF = CutMotherTracks_PrimVtx[i].GetNDF();
      temp_Hyp_Cut.MomE.SetPxPyPzE(CutMotherTracks_PrimVtx[i].GetPx(),CutMotherTracks_PrimVtx[i].GetPy(),CutMotherTracks_PrimVtx[i].GetPz(),CutMotherTracks_PrimVtx[i].GetE());
      CutMotherTracks_PrimVtx[i].TransportToProductionVertex();
      temp_Hyp_Cut.PrimVtx.SetXYZ(CutMotherTracks_PrimVtx[i].GetX(),CutMotherTracks_PrimVtx[i].GetY(),CutMotherTracks_PrimVtx[i].GetZ());
      CutMotherTracks_PrimVtx[i].TransportToDecayVertex();
      temp_Hyp_Cut.DecayVtx.SetXYZ(CutMotherTracks_PrimVtx[i].GetX(),CutMotherTracks_PrimVtx[i].GetY(),CutMotherTracks_PrimVtx[i].GetZ());
      temp_Hyp_Cut.Dist_RealReconsVtx.SetXYZ(DecayVertex_real.X()-CutMotherTracks_PrimVtx[i].GetX(),DecayVertex_real.Y()-CutMotherTracks_PrimVtx[i].GetY(),DecayVertex_real.Z()-CutMotherTracks_PrimVtx[i].GetZ());
      ThetaDist_TrackPrimVtx(CutMotherTracks_PrimVtx[i], InteractionPoint_real, theta_cutDecayTrackPrimVtx, dist_cutDecayTrackPrimVtx);
      temp_Hyp_Cut.Dist_MotherPrimVtx = dist_cutDecayTrackPrimVtx;
      temp_Hyp_Cut.Angle_MotherPrimVtx = theta_cutDecayTrackPrimVtx;
      temp_Hyp_Cut.ErrGetMass = CutMotherTracks_PrimVtx[i].GetMass(CutMother_InvMass, CutMother_ErrInvMass);
      temp_Hyp_Cut.InvMass = CutMother_InvMass;
      temp_Hyp_Cut.ErrInvMass = CutMother_ErrInvMass;
      temp_Hyp_Cut.ErrGetLifeTime = CutMotherTracks_PrimVtx[i].GetLifeTime(CutMother_LifeTime, CutMother_ErrLifeTime);
      temp_Hyp_Cut.LifeTime = CutMother_LifeTime / c_light_speed_cmps;
      temp_Hyp_Cut.ErrLifeTime = CutMother_ErrLifeTime / c_light_speed_cmps;

      TVector3 Mother_Mom(CutMotherTracks_PrimVtx[i].GetPx(), CutMotherTracks_PrimVtx[i].GetPy(), CutMotherTracks_PrimVtx[i].GetPz());

      size_t temp_id_fragment = get<0>(RefCutDaughtersTracks_PrimVtx[i]);
      size_t temp_id_pion = get<1>(RefCutDaughtersTracks_PrimVtx[i]);

      temp_Hyp_Cut.PDG_Fragment = FragmentTracks[temp_id_fragment].GetPDG();
      temp_Hyp_Cut.Id_Fragment = FragmentTracks[temp_id_fragment].Id();
      temp_Hyp_Cut.MomE_Fragment.SetPxPyPzE(FragmentTracks[temp_id_fragment].GetPx(), FragmentTracks[temp_id_fragment].GetPy(), FragmentTracks[temp_id_fragment].GetPz(), FragmentTracks[temp_id_fragment].GetE());
      temp_Hyp_Cut.Chi2ndf_Fragment = FragmentTracks[temp_id_fragment].GetChi2() / static_cast<double>(FragmentTracks[temp_id_fragment].GetNDF());      
      temp_Hyp_Cut.NDF_Fragment = FragmentTracks[temp_id_fragment].GetNDF();
      temp_Hyp_Cut.Pvalue_Fragment = Fragment_FitInfo[temp_id_fragment].Pvalue;
      TVector3 Fragment_Mom(FragmentTracks[temp_id_fragment].GetPx(), FragmentTracks[temp_id_fragment].GetPy(), FragmentTracks[temp_id_fragment].GetPz());
      temp_Hyp_Cut.Angle_MotherFragment = Mother_Mom.Angle(Fragment_Mom) * 180. / M_PI;
      temp_Hyp_Cut.Fragment_IsFromHyp = 0;
      if(FragmentTracks[temp_id_fragment].Id() == FragmentTracks[ref_RealFragment].Id())
        temp_Hyp_Cut.Fragment_IsFromHyp = 1;        

      temp_Hyp_Cut.PDG_Pion = CutPionTracks[temp_id_pion].GetPDG();
      temp_Hyp_Cut.Id_Pion = CutPionTracks[temp_id_pion].Id();
      temp_Hyp_Cut.MomE_Pion.SetPxPyPzE(CutPionTracks[temp_id_pion].GetPx(), CutPionTracks[temp_id_pion].GetPy(), CutPionTracks[temp_id_pion].GetPz(), CutPionTracks[temp_id_pion].GetE());
      temp_Hyp_Cut.Chi2ndf_Pion = CutPionTracks[temp_id_pion].GetChi2() / static_cast<double>(CutPionTracks[temp_id_pion].GetNDF());
      temp_Hyp_Cut.NDF_Pion = CutPionTracks[temp_id_pion].GetNDF();
      temp_Hyp_Cut.Pvalue_Pion = CutPion_FitInfo[temp_id_pion].Pvalue;
      TVector3 Pion_Mom(CutPionTracks[temp_id_pion].GetPx(), CutPionTracks[temp_id_pion].GetPy(), CutPionTracks[temp_id_pion].GetPz());
      temp_Hyp_Cut.Angle_MotherPion = Mother_Mom.Angle(Pion_Mom) * 180. / M_PI;
      temp_Hyp_Cut.NHitsMDC_Pion = CutPion_FitInfo[temp_id_pion].NHitsMDC;
      temp_Hyp_Cut.NHitsMinifiber_Pion = CutPion_FitInfo[temp_id_pion].NHitsMinifiber;
      temp_Hyp_Cut.N_Pion = CutPionTracks.size();
      temp_Hyp_Cut.Pion_IsFromHyp = 0;

      for(itr_cut = RecoEvent.DaughtersTrackDAFInit.begin(); itr_cut != RecoEvent.DaughtersTrackDAFInit.end(); ++itr_cut)
        {
          if(itr_cut->first == CutPionTracks[temp_id_pion].Id())
            temp_Hyp_Cut.Pion_IsFromHyp = 1;
        }

      temp_Hyp_Cut.Mother_IsFromHyp = temp_Hyp_Cut.Fragment_IsFromHyp * temp_Hyp_Cut.Pion_IsFromHyp;

      float Vertex[3] = {CutMotherTracks_PrimVtx[i].GetX(), CutMotherTracks_PrimVtx[i].GetY(), CutMotherTracks_PrimVtx[i].GetZ()};
      FragmentTracks[temp_id_fragment].TransportToPoint(Vertex);
      CutPionTracks[temp_id_pion].TransportToPoint(Vertex);

      CloseDist(FragmentTracks[temp_id_fragment], CutPionTracks[temp_id_pion], dist_CutDaughters, centroid_CutDaughters);
      temp_Hyp_Cut.Dist_Daughters = dist_CutDaughters;

      float armenterosQtAlfa[2] = {0.};
      CutMotherTracks_PrimVtx[i].GetArmenterosPodolanski_FromMother(FragmentTracks[temp_id_fragment], CutPionTracks[temp_id_pion], armenterosQtAlfa);
      temp_Hyp_Cut.ArmPod_Qt = armenterosQtAlfa[0];
      temp_Hyp_Cut.ArmPod_Alfa = armenterosQtAlfa[1];

      RecoEvent.Hyp_Vect.emplace_back(temp_Hyp_Cut);


      double cutdistance_KFPart_PrimVtx  = std::sqrt(std::pow((DecayVertex_real.X() - CutMotherTracks_PrimVtx[i].GetX()), 2.) +
                              std::pow((DecayVertex_real.Y() - CutMotherTracks_PrimVtx[i].GetY()), 2.) +
                              std::pow((DecayVertex_real.Z() - CutMotherTracks_PrimVtx[i].GetZ()), 2.));
      double cutdistanceX_KFPart_PrimVtx = DecayVertex_real.X() - CutMotherTracks_PrimVtx[i].GetX();
      double cutdistanceY_KFPart_PrimVtx = DecayVertex_real.Y() - CutMotherTracks_PrimVtx[i].GetY();
      double cutdistanceZ_KFPart_PrimVtx = DecayVertex_real.Z() - CutMotherTracks_PrimVtx[i].GetZ();

      LocalHisto.h_DecayVertexcutDistance_KFPart_PrimVtx[pion_type]->Fill(cutdistance_KFPart_PrimVtx, 1.);
      LocalHisto.h_DecayVertexcutDistanceX_KFPart_PrimVtx[pion_type]->Fill(cutdistanceX_KFPart_PrimVtx, 1.);
      LocalHisto.h_DecayVertexcutDistanceY_KFPart_PrimVtx[pion_type]->Fill(cutdistanceY_KFPart_PrimVtx, 1.);
      LocalHisto.h_DecayVertexcutDistanceZ_KFPart_PrimVtx[pion_type]->Fill(cutdistanceZ_KFPart_PrimVtx, 1.);

      LocalHisto.h_EffPosZ_preKF[pion_type]->Fill(CutMotherTracks_PrimVtx[i].GetZ() - att.Target_PositionZ, 1.);

      double cutTimeLife_PrimVtx = CutMotherTracks_PrimVtx[i].GetLifeTime() / c_light_speed_cmps; //in ps
      LocalHisto.h_HypcutLifeTime_PrimVtx[pion_type]->Fill(cutTimeLife_PrimVtx, 1.);

      LocalHisto.h_Hyp_CutArmenterosPodolanski[pion_type]->Fill(armenterosQtAlfa[1], armenterosQtAlfa[0], 1.);
    }

    TLorentzVector CutMother_LorentzVector;
    CutMother_LorentzVector = Fragment_LorentzVector + CutPion_LorentzVector;   
    LocalHisto.h_HypInvariantMass_CutLorentzVect[pion_type]->Fill(CutMother_LorentzVector.M(), 1.);

#endif

  //Pion tracks
  std::vector<KFParticle> PionTracks_All {};
  std::vector<KFFitInfo> Pion_FitInfo {};
  PionTracksFinder(RecoEvent.DAF_results, PionTracks_All, Pion_FitInfo);
  LocalHisto.h_Npions[pion_type]->Fill(PionTracks_All.size(), 1.);
  if(PionTracks_All.size() == 0)
    return -4;

  std::vector<KFParticle> PionTracks {};
  PionSelector(PionTracks_All, RecoEvent.PrimVtxRecons, PionTracks, Pion_FitInfo);
  if(PionTracks.size() == 0)
    return -4;

  double closedist_distance = 0.;
  TVector3 temp_closedist_pos;

  double newclosest_distance = 1000.;
  TVector3 new_closedist_pos;
  size_t new_pionid;

  std::vector<TVector3> vect_closedist_pos;
  double theta_DecayTrackPrimVtx;
  double dist_DecayTrackPrimVtx;

  std::unordered_map<int, InfoInit>::iterator itr_recons;
  bool ifDaughter_recons = false;

  TLorentzVector Pion_LorentzVector;

  for(size_t i = 0; i < PionTracks.size(); ++i)
    {
      LocalHisto.h_P_pions[pion_type]->Fill(PionTracks[i].GetP(), 1.);
      LocalHisto.h_Pt_pions[pion_type]->Fill(PionTracks[i].GetPt(), 1.);
      LocalHisto.h_Pz_pions[pion_type]->Fill(PionTracks[i].GetPz(), 1.);

      CloseDist(FragmentTracks[ref_RealFragment], PionTracks[i], closedist_distance, temp_closedist_pos);
      vect_closedist_pos.emplace_back(temp_closedist_pos);

      LocalHisto.h_Closedist_Distance[pion_type]->Fill(closedist_distance, 1.);
      LocalHisto.h_Closedist_PosZ[pion_type]->Fill(temp_closedist_pos.Z() - att.Target_PositionZ, 1.);

      ThetaDist_TrackPrimVtx(PionTracks[i], RecoEvent.PrimVtxRecons, theta_DecayTrackPrimVtx, dist_DecayTrackPrimVtx);
      LocalHisto.h_Dist_DecayTrackPrimVtx[pion_type]->Fill(dist_DecayTrackPrimVtx, "All", 1.);

      for(itr_recons = RecoEvent.DaughtersTrackDAFInit.begin(); itr_recons != RecoEvent.DaughtersTrackDAFInit.end(); ++itr_recons)
        {
          if(itr_recons->first == PionTracks[i].Id())
            ifDaughter_recons = true;
        }

      if(ifDaughter_recons)
        LocalHisto.h_Dist_DecayTrackPrimVtx[pion_type]->Fill(dist_DecayTrackPrimVtx, "Daughters", 1.);
      else
        LocalHisto.h_Dist_DecayTrackPrimVtx[pion_type]->Fill(dist_DecayTrackPrimVtx, "Primaries", 1.);

      if(closedist_distance < newclosest_distance)
        {
          newclosest_distance = closedist_distance;
          new_closedist_pos = temp_closedist_pos;
          new_pionid = i;

          Pion_LorentzVector.SetPxPyPzE(PionTracks[i].GetPx(),PionTracks[i].GetPy(),PionTracks[i].GetPz(),PionTracks[i].GetE());
        }
    }

    double dist_DaughtersLV;
    TVector3 centroid_DaughtersLV;

    if(Pion_LorentzVector.P() > 0.0001)
      {
        TLorentzVector Mother_LorentzVector;
        Mother_LorentzVector = Fragment_LorentzVector + Pion_LorentzVector;
        LocalHisto.h_HypInvariantMass_LorentzVect[pion_type]->Fill(Mother_LorentzVector.M(), 1.);

        Hyp temp_Hyp_LV;

        temp_Hyp_LV.Pattern = 6;

        temp_Hyp_LV.PDG = Hyp_pdg;
        temp_Hyp_LV.N_Mother = 1; //Change!
        temp_Hyp_LV.Chi2ndf = -1.; //Change!
        temp_Hyp_LV.NDF = -999; //Change!
        temp_Hyp_LV.MomE.SetPxPyPzE(Mother_LorentzVector.Px(),Mother_LorentzVector.Py(),Mother_LorentzVector.Pz(),Mother_LorentzVector.E());
        temp_Hyp_LV.PrimVtx.SetXYZ(0.,0.,0.); //Change!
        temp_Hyp_LV.DecayVtx.SetXYZ(temp_closedist_pos.X(),temp_closedist_pos.Y(),temp_closedist_pos.Z());
        temp_Hyp_LV.Dist_RealReconsVtx.SetXYZ(DecayVertex_real.X()-temp_closedist_pos.X(),DecayVertex_real.Y()-temp_closedist_pos.Y(),DecayVertex_real.Z()-temp_closedist_pos.Z());
        temp_Hyp_LV.Dist_MotherPrimVtx = -1.; //Change!
        temp_Hyp_LV.Angle_MotherPrimVtx = -1.; //Change!
        temp_Hyp_LV.ErrGetMass = -1; //Change!
        temp_Hyp_LV.InvMass = Mother_LorentzVector.M();
        temp_Hyp_LV.ErrInvMass = -1.; //Change!
        temp_Hyp_LV.ErrGetLifeTime = -1; //Change!
        temp_Hyp_LV.LifeTime = -1.; //Change!
        temp_Hyp_LV.ErrLifeTime = -1.; //Change!

        TVector3 Mother_Mom(Mother_LorentzVector.Px(), Mother_LorentzVector.Py(), Mother_LorentzVector.Pz());

        temp_Hyp_LV.PDG_Fragment = FragmentTracks[ref_RealFragment].GetPDG();
        temp_Hyp_LV.Id_Fragment = FragmentTracks[ref_RealFragment].Id();
        temp_Hyp_LV.MomE_Fragment.SetPxPyPzE(FragmentTracks[ref_RealFragment].GetPx(), FragmentTracks[ref_RealFragment].GetPy(), FragmentTracks[ref_RealFragment].GetPz(), FragmentTracks[ref_RealFragment].GetE());
        temp_Hyp_LV.Chi2ndf_Fragment = FragmentTracks[ref_RealFragment].GetChi2() / static_cast<double>(FragmentTracks[ref_RealFragment].GetNDF());      
        temp_Hyp_LV.NDF_Fragment = FragmentTracks[ref_RealFragment].GetNDF();
        temp_Hyp_LV.Pvalue_Fragment = Fragment_FitInfo[ref_RealFragment].Pvalue;
        TVector3 Fragment_Mom(FragmentTracks[ref_RealFragment].GetPx(), FragmentTracks[ref_RealFragment].GetPy(), FragmentTracks[ref_RealFragment].GetPz());
        temp_Hyp_LV.Angle_MotherFragment = Mother_Mom.Angle(Fragment_Mom) * 180. / M_PI;
        temp_Hyp_LV.Fragment_IsFromHyp = 1;

        temp_Hyp_LV.PDG_Pion = PionTracks[new_pionid].GetPDG();
        temp_Hyp_LV.Id_Pion = PionTracks[new_pionid].Id();
        temp_Hyp_LV.MomE_Pion.SetPxPyPzE(PionTracks[new_pionid].GetPx(), PionTracks[new_pionid].GetPy(), PionTracks[new_pionid].GetPz(), PionTracks[new_pionid].GetE());
        temp_Hyp_LV.Chi2ndf_Pion = PionTracks[new_pionid].GetChi2() / static_cast<double>(PionTracks[new_pionid].GetNDF());
        temp_Hyp_LV.NDF_Pion = PionTracks[new_pionid].GetNDF();
        temp_Hyp_LV.Pvalue_Pion = Pion_FitInfo[new_pionid].Pvalue;
        TVector3 Pion_Mom(PionTracks[new_pionid].GetPx(), PionTracks[new_pionid].GetPy(), PionTracks[new_pionid].GetPz());
        temp_Hyp_LV.Angle_MotherPion = Mother_Mom.Angle(Pion_Mom) * 180. / M_PI;
        temp_Hyp_LV.NHitsMDC_Pion = Pion_FitInfo[new_pionid].NHitsMDC;
        temp_Hyp_LV.NHitsMinifiber_Pion = Pion_FitInfo[new_pionid].NHitsMinifiber;
        temp_Hyp_LV.N_Pion = 1;
        temp_Hyp_LV.Pion_IsFromHyp = 0;

        for(itr_recons = RecoEvent.DaughtersTrackDAFInit.begin(); itr_recons != RecoEvent.DaughtersTrackDAFInit.end(); ++itr_recons)
          {
            if(itr_recons->first == PionTracks[new_pionid].Id())
              temp_Hyp_LV.Pion_IsFromHyp = 1;
          }

        temp_Hyp_LV.Mother_IsFromHyp = temp_Hyp_LV.Fragment_IsFromHyp * temp_Hyp_LV.Pion_IsFromHyp;

        CloseDist(FragmentTracks[ref_RealFragment], PionTracks[new_pionid], dist_DaughtersLV, centroid_DaughtersLV);
        temp_Hyp_LV.Dist_Daughters = dist_DaughtersLV;

        temp_Hyp_LV.ArmPod_Qt = -1.; //Change!
        temp_Hyp_LV.ArmPod_Alfa = -1.; //Change!

        RecoEvent.Hyp_Vect.emplace_back(temp_Hyp_LV);
      }

#ifdef CENTROID_METHOD

  TVector3 closedist_pos;

  for(size_t i = 0; i < vect_closedist_pos.size(); ++i)
    {
      closedist_pos += vect_closedist_pos[i];
    }
  closedist_pos *= (1./static_cast<double>(vect_closedist_pos.size()));

  double distance_centroid  = std::sqrt(std::pow((DecayVertex_real.X() - closedist_pos.X()), 2.) +
                                   std::pow((DecayVertex_real.Y() - closedist_pos.Y()), 2.) +
                                   std::pow((DecayVertex_real.Z() - closedist_pos.Z()), 2.));
  double distanceX_centroid = DecayVertex_real.X() - closedist_pos.X();
  double distanceY_centroid = DecayVertex_real.Y() - closedist_pos.Y();
  double distanceZ_centroid = DecayVertex_real.Z() - closedist_pos.Z();

  LocalHisto.h_DecayVertexDistance_centroid[pion_type]->Fill(distance_centroid, 1.);
  LocalHisto.h_DecayVertexDistanceX_centroid[pion_type]->Fill(distanceX_centroid, 1.);
  LocalHisto.h_DecayVertexDistanceY_centroid[pion_type]->Fill(distanceY_centroid, 1.);
  LocalHisto.h_DecayVertexDistanceZ_centroid[pion_type]->Fill(distanceZ_centroid, 1.);

  LocalHisto.h_DecayVertexPosZ_centroid[pion_type]->Fill(closedist_pos.Z() - att.Target_PositionZ, 1.);

#endif

#ifdef VFUNCTION_METHOD

  //Decay vertex reconstruction
  TVector3 DecayVertexRecons;

  TrackstoDecayVertex(FragmentTracks, PionTracks, RecoEvent.PrimVtxRecons, DecayVertexRecons);
  RecoEvent.DecayVtxRecons.SetXYZ(DecayVertexRecons.X(), DecayVertexRecons.Y(), DecayVertexRecons.Z());

  double distance  = (DecayVertex_real - DecayVertexRecons).Mag();
  double distanceX = DecayVertex_real.X() - DecayVertexRecons.X();
  double distanceY = DecayVertex_real.Y() - DecayVertexRecons.Y();
  double distanceZ = DecayVertex_real.Z() - DecayVertexRecons.Z();

  LocalHisto.h_DecayVertexDistance[pion_type]->Fill(distance, 1.);
  LocalHisto.h_DecayVertexDistanceX[pion_type]->Fill(distanceX, 1.);
  LocalHisto.h_DecayVertexDistanceY[pion_type]->Fill(distanceY, 1.);
  LocalHisto.h_DecayVertexDistanceZ[pion_type]->Fill(distanceZ, 1.);

  LocalHisto.h_DecayVertexPosZ_vfunction[pion_type]->Fill(DecayVertexRecons.Z() - att.Target_PositionZ, 1.);

#endif


  //Hypernucleus reconstruction
  std::vector<KFParticle> MotherTracks_All;
  std::vector<std::tuple<size_t, size_t>> RefDaughtersTracks_All;
  ifSet_ProductionVertex = false;
  ifSet_MassConstraint = false;
  MotherTracksRecons(FragmentTracks, PionTracks, nullptr, MotherTracks_All, RefDaughtersTracks_All);

  std::vector<KFParticle> MotherTracks;
  std::vector<std::tuple<size_t, size_t>> RefDaughtersTracks;
  MotherSelector(MotherTracks_All, RefDaughtersTracks_All, FragmentTracks, PionTracks, RecoEvent.PrimVtxRecons,
                 MotherTracks, RefDaughtersTracks);

  for(size_t i = 0; i < MotherTracks.size(); ++i)
    {
      double distance_KFPart  = std::sqrt(std::pow((DecayVertex_real.X() - MotherTracks[i].GetX()), 2.) +
                              std::pow((DecayVertex_real.Y() - MotherTracks[i].GetY()), 2.) +
                              std::pow((DecayVertex_real.Z() - MotherTracks[i].GetZ()), 2.));
      double distanceX_KFPart = DecayVertex_real.X() - MotherTracks[i].GetX();
      double distanceY_KFPart = DecayVertex_real.Y() - MotherTracks[i].GetY();
      double distanceZ_KFPart = DecayVertex_real.Z() - MotherTracks[i].GetZ();

      LocalHisto.h_DecayVertexDistance_KFPart[pion_type]->Fill(distance_KFPart, 1.);
      LocalHisto.h_DecayVertexDistanceX_KFPart[pion_type]->Fill(distanceX_KFPart, 1.);
      LocalHisto.h_DecayVertexDistanceY_KFPart[pion_type]->Fill(distanceY_KFPart, 1.);
      LocalHisto.h_DecayVertexDistanceZ_KFPart[pion_type]->Fill(distanceZ_KFPart, 1.);

      LocalHisto.h_DecayVertexPosZ_KFPart[pion_type]->Fill(MotherTracks[i].GetZ() - att.Target_PositionZ, 1.);
    }


  std::vector<KFParticle> MotherTracks_PrimVtx_All;
  std::vector<std::tuple<size_t, size_t>> RefDaughtersTracks_PrimVtx_All;
  ifSet_ProductionVertex = true;
  ifSet_MassConstraint = false;
  MotherTracksRecons(FragmentTracks, PionTracks, pointer_PrimVtx, MotherTracks_PrimVtx_All, RefDaughtersTracks_PrimVtx_All);

  std::vector<KFParticle> MotherTracks_PrimVtx;
  std::vector<std::tuple<size_t, size_t>> RefDaughtersTracks_PrimVtx;
  MotherSelector(MotherTracks_PrimVtx_All, RefDaughtersTracks_PrimVtx_All, FragmentTracks, PionTracks, RecoEvent.PrimVtxRecons,
                 MotherTracks_PrimVtx, RefDaughtersTracks_PrimVtx);

  float Mother_InvMass;
  float Mother_ErrInvMass;
  float Mother_LifeTime;
  float Mother_ErrLifeTime;

  double dist_Daughters;
  TVector3 centroid_Daughters;

  for(size_t i = 0; i < MotherTracks_PrimVtx.size(); ++i)
    {
      Hyp temp_Hyp;

      temp_Hyp.Pattern = 4;

      temp_Hyp.PDG = MotherTracks_PrimVtx[i].GetPDG();
      temp_Hyp.N_Mother = MotherTracks_PrimVtx.size();
      temp_Hyp.Chi2ndf = MotherTracks_PrimVtx[i].GetChi2() / static_cast<double>(MotherTracks_PrimVtx[i].GetNDF());
      temp_Hyp.NDF = MotherTracks_PrimVtx[i].GetNDF();
      temp_Hyp.MomE.SetPxPyPzE(MotherTracks_PrimVtx[i].GetPx(),MotherTracks_PrimVtx[i].GetPy(),MotherTracks_PrimVtx[i].GetPz(),MotherTracks_PrimVtx[i].GetE());
      MotherTracks_PrimVtx[i].TransportToProductionVertex();
      temp_Hyp.PrimVtx.SetXYZ(MotherTracks_PrimVtx[i].GetX(),MotherTracks_PrimVtx[i].GetY(),MotherTracks_PrimVtx[i].GetZ());
      MotherTracks_PrimVtx[i].TransportToDecayVertex();
      temp_Hyp.DecayVtx.SetXYZ(MotherTracks_PrimVtx[i].GetX(),MotherTracks_PrimVtx[i].GetY(),MotherTracks_PrimVtx[i].GetZ());
      temp_Hyp.Dist_RealReconsVtx.SetXYZ(DecayVertex_real.X()-MotherTracks_PrimVtx[i].GetX(),DecayVertex_real.Y()-MotherTracks_PrimVtx[i].GetY(),DecayVertex_real.Z()-MotherTracks_PrimVtx[i].GetZ());
      ThetaDist_TrackPrimVtx(MotherTracks_PrimVtx[i], InteractionPoint_real, theta_DecayTrackPrimVtx, dist_DecayTrackPrimVtx);
      temp_Hyp.Dist_MotherPrimVtx = dist_DecayTrackPrimVtx;
      temp_Hyp.Angle_MotherPrimVtx = theta_DecayTrackPrimVtx;
      temp_Hyp.ErrGetMass = MotherTracks_PrimVtx[i].GetMass(Mother_InvMass, Mother_ErrInvMass);
      temp_Hyp.InvMass = Mother_InvMass;
      temp_Hyp.ErrInvMass = Mother_ErrInvMass;
      temp_Hyp.ErrGetLifeTime = MotherTracks_PrimVtx[i].GetLifeTime(Mother_LifeTime, Mother_ErrLifeTime);
      temp_Hyp.LifeTime = Mother_LifeTime / c_light_speed_cmps;
      temp_Hyp.ErrLifeTime = Mother_ErrLifeTime / c_light_speed_cmps;

      TVector3 Mother_Mom(MotherTracks_PrimVtx[i].GetPx(), MotherTracks_PrimVtx[i].GetPy(), MotherTracks_PrimVtx[i].GetPz());

      size_t temp_id_fragment = get<0>(RefDaughtersTracks_PrimVtx[i]);
      size_t temp_id_pion = get<1>(RefDaughtersTracks_PrimVtx[i]);

      temp_Hyp.PDG_Fragment = FragmentTracks[temp_id_fragment].GetPDG();
      temp_Hyp.Id_Fragment = FragmentTracks[temp_id_fragment].Id();
      temp_Hyp.MomE_Fragment.SetPxPyPzE(FragmentTracks[temp_id_fragment].GetPx(), FragmentTracks[temp_id_fragment].GetPy(), FragmentTracks[temp_id_fragment].GetPz(), FragmentTracks[temp_id_fragment].GetE());
      temp_Hyp.Chi2ndf_Fragment = FragmentTracks[temp_id_fragment].GetChi2() / static_cast<double>(FragmentTracks[temp_id_fragment].GetNDF());      
      temp_Hyp.NDF_Fragment = FragmentTracks[temp_id_fragment].GetNDF();
      temp_Hyp.Pvalue_Fragment = Fragment_FitInfo[temp_id_fragment].Pvalue;
      TVector3 Fragment_Mom(FragmentTracks[temp_id_fragment].GetPx(), FragmentTracks[temp_id_fragment].GetPy(), FragmentTracks[temp_id_fragment].GetPz());
      temp_Hyp.Angle_MotherFragment = Mother_Mom.Angle(Fragment_Mom) * 180. / M_PI;
      temp_Hyp.Fragment_IsFromHyp = 0;
      if(FragmentTracks[temp_id_fragment].Id() == FragmentTracks[ref_RealFragment].Id())
        temp_Hyp.Fragment_IsFromHyp = 1;

      temp_Hyp.PDG_Pion = PionTracks[temp_id_pion].GetPDG();
      temp_Hyp.Id_Pion = PionTracks[temp_id_pion].Id();
      temp_Hyp.MomE_Pion.SetPxPyPzE(PionTracks[temp_id_pion].GetPx(), PionTracks[temp_id_pion].GetPy(), PionTracks[temp_id_pion].GetPz(), PionTracks[temp_id_pion].GetE());
      temp_Hyp.Chi2ndf_Pion = PionTracks[temp_id_pion].GetChi2() / static_cast<double>(PionTracks[temp_id_pion].GetNDF());
      temp_Hyp.NDF_Pion = PionTracks[temp_id_pion].GetNDF();
      temp_Hyp.Pvalue_Pion = Pion_FitInfo[temp_id_pion].Pvalue;
      TVector3 Pion_Mom(PionTracks[temp_id_pion].GetPx(), PionTracks[temp_id_pion].GetPy(), PionTracks[temp_id_pion].GetPz());
      temp_Hyp.Angle_MotherPion = Mother_Mom.Angle(Pion_Mom) * 180. / M_PI;
      temp_Hyp.NHitsMDC_Pion = Pion_FitInfo[temp_id_pion].NHitsMDC;
      temp_Hyp.NHitsMinifiber_Pion = Pion_FitInfo[temp_id_pion].NHitsMinifiber;
      temp_Hyp.N_Pion = PionTracks.size();
      temp_Hyp.Pion_IsFromHyp = 0;

      for(itr_recons = RecoEvent.DaughtersTrackDAFInit.begin(); itr_recons != RecoEvent.DaughtersTrackDAFInit.end(); ++itr_recons)
        {
          if(itr_recons->first == PionTracks[temp_id_pion].Id())
            temp_Hyp.Pion_IsFromHyp = 1;
        }

      temp_Hyp.Mother_IsFromHyp = temp_Hyp.Fragment_IsFromHyp * temp_Hyp.Pion_IsFromHyp;

      float Vertex[3] = {MotherTracks_PrimVtx[i].GetX(), MotherTracks_PrimVtx[i].GetY(), MotherTracks_PrimVtx[i].GetZ()};
      FragmentTracks[temp_id_fragment].TransportToPoint(Vertex);
      PionTracks[temp_id_pion].TransportToPoint(Vertex);

      CloseDist(FragmentTracks[temp_id_fragment], PionTracks[temp_id_pion], dist_Daughters, centroid_Daughters);
      temp_Hyp.Dist_Daughters = dist_Daughters;

      float armenterosQtAlfa[2] = {0.};
      MotherTracks_PrimVtx[i].GetArmenterosPodolanski_FromMother(FragmentTracks[temp_id_fragment], PionTracks[temp_id_pion], armenterosQtAlfa);
      temp_Hyp.ArmPod_Qt = armenterosQtAlfa[0];
      temp_Hyp.ArmPod_Alfa = armenterosQtAlfa[1];

      RecoEvent.Hyp_Vect.emplace_back(temp_Hyp);


      double distance_KFPart_PrimVtx  = std::sqrt(std::pow((DecayVertex_real.X() - MotherTracks_PrimVtx[i].GetX()), 2.) +
                              std::pow((DecayVertex_real.Y() - MotherTracks_PrimVtx[i].GetY()), 2.) +
                              std::pow((DecayVertex_real.Z() - MotherTracks_PrimVtx[i].GetZ()), 2.));
      double distanceX_KFPart_PrimVtx = DecayVertex_real.X() - MotherTracks_PrimVtx[i].GetX();
      double distanceY_KFPart_PrimVtx = DecayVertex_real.Y() - MotherTracks_PrimVtx[i].GetY();
      double distanceZ_KFPart_PrimVtx = DecayVertex_real.Z() - MotherTracks_PrimVtx[i].GetZ();

      LocalHisto.h_DecayVertexDistance_KFPart_PrimVtx[pion_type]->Fill(distance_KFPart_PrimVtx, 1.);
      LocalHisto.h_DecayVertexDistanceX_KFPart_PrimVtx[pion_type]->Fill(distanceX_KFPart_PrimVtx, 1.);
      LocalHisto.h_DecayVertexDistanceY_KFPart_PrimVtx[pion_type]->Fill(distanceY_KFPart_PrimVtx, 1.);
      LocalHisto.h_DecayVertexDistanceZ_KFPart_PrimVtx[pion_type]->Fill(distanceZ_KFPart_PrimVtx, 1.);

      //Check the possible cuts
      double Closedist_DaughterTracks;
      TVector3 Centroid_DaughtersTracks;
      CloseDist(FragmentTracks[temp_id_fragment], PionTracks[temp_id_pion], Closedist_DaughterTracks, Centroid_DaughtersTracks);
      LocalHisto.h_Dist_DaughterTracks[pion_type]->Fill(Closedist_DaughterTracks, MotherTracks_PrimVtx[i].GetMass(), 1.);

      double angle_MotherFragment = MotherTracks_PrimVtx[i].GetAngle(FragmentTracks[temp_id_fragment]) * 180. / M_PI;
      LocalHisto.h_Angle_MotherFragment[pion_type]->Fill(angle_MotherFragment, MotherTracks_PrimVtx[i].GetMass(), 1.);

      double angle_MotherPion = MotherTracks_PrimVtx[i].GetAngle(PionTracks[temp_id_pion]) * 180. / M_PI;
      LocalHisto.h_Angle_MotherPion[pion_type]->Fill(angle_MotherPion, MotherTracks_PrimVtx[i].GetMass(), 1.);

      double temp_chi2ndf = MotherTracks_PrimVtx[i].GetChi2() / static_cast<double>(MotherTracks_PrimVtx[i].GetNDF());
      LocalHisto.h_Chi2ndf_MotherTracks[pion_type]->Fill(temp_chi2ndf, MotherTracks_PrimVtx[i].GetMass(), 1.);

      double theta_MotherTrackPrimVtx;
      double dist_MotherTrackPrimVtx;
      ThetaDist_TrackPrimVtx(MotherTracks_PrimVtx[i], RecoEvent.PrimVtxRecons, theta_MotherTrackPrimVtx, dist_MotherTrackPrimVtx);
      LocalHisto.h_Dist_MotherTrackPrimVtx[pion_type]->Fill(dist_MotherTrackPrimVtx, MotherTracks_PrimVtx[i].GetMass(), 1.);
      LocalHisto.h_Theta_MotherTrackPrimVtx[pion_type]->Fill(theta_MotherTrackPrimVtx, MotherTracks_PrimVtx[i].GetMass(), 1.);

      LocalHisto.h_EffPosZ_postKFPart[pion_type]->Fill(MotherTracks_PrimVtx[i].GetZ() - att.Target_PositionZ, 1.);
      double PosR_DecayVertex = std::sqrt(std::pow(MotherTracks_PrimVtx[i].GetX() - att.Target_PositionX, 2.)
                                            + std::pow(MotherTracks_PrimVtx[i].GetY() - att.Target_PositionY, 2.));
      LocalHisto.h_EffPosZPosR_postKFPart[pion_type]->Fill(MotherTracks_PrimVtx[i].GetZ() - att.Target_PositionZ, PosR_DecayVertex, 1.);
      LocalHisto.h_DecayVertexPosZ_KFPart_PrimVtx[pion_type]->Fill(MotherTracks_PrimVtx[i].GetZ() - att.Target_PositionZ, MotherTracks_PrimVtx[i].GetMass(), 1.);

      LocalHisto.h_DecayFragmentMomZ_KFPart_PrimVtx[pion_type]->Fill(FragmentTracks[temp_id_fragment].GetPz(), MotherTracks_PrimVtx[i].GetMass(), 1.);
      LocalHisto.h_DecayPionMomZ_KFPart_PrimVtx[pion_type]->Fill(PionTracks[temp_id_pion].GetPz(), MotherTracks_PrimVtx[i].GetMass(), 1.);

      LocalHisto.h_N_MotherTracks[pion_type]->Fill(MotherTracks_PrimVtx.size(), MotherTracks_PrimVtx[i].GetMass(), 1.);

/*
      std::cout << "Mother mass: " << MotherTracks_PrimVtx[i].GetMass() << "\n";
      std::cout << "Mother ErrM: " << MotherTracks_PrimVtx[i].GetErrMass() << "\n";
      std::cout << "Mother Ener: " << MotherTracks_PrimVtx[i].GetE() << "\n";
      std::cout << "Mother MomP: " << MotherTracks_PrimVtx[i].GetP() << "\n\n";
*/
/*
      if(std::isnan(MotherTracks_PrimVtx[i].GetMass()))
      {
        std::cout << "---------InvMass: " << MotherTracks_PrimVtx[i].GetMass() << "\n";
        std::cout << "---------RelPosX: " << DecayVertex_real.X()    << "\n";
        std::cout << "---------RelPosY: " << DecayVertex_real.Y()    << "\n";
        std::cout << "---------RelPosZ: " << DecayVertex_real.Z()    << "\n";
        std::cout << "---------DcyPosX: " << MotherTracks_PrimVtx[i].GetX()    << "\n";
        std::cout << "---------DcyPosY: " << MotherTracks_PrimVtx[i].GetY()    << "\n";
        std::cout << "---------DcyPosZ: " << MotherTracks_PrimVtx[i].GetZ()    << "\n";
        std::cout << "---------FraPosX: " << FragmentTracks[temp_id_fragment].GetX() << "\n";
        std::cout << "---------FraPosY: " << FragmentTracks[temp_id_fragment].GetY() << "\n";
        std::cout << "---------FraPosZ: " << FragmentTracks[temp_id_fragment].GetZ() << "\n";
        std::cout << "---------PioPosX: " << PionTracks[temp_id_pion].GetX() << "\n";
        std::cout << "---------PioPosY: " << PionTracks[temp_id_pion].GetY() << "\n";
        std::cout << "---------PioPosZ: " << PionTracks[temp_id_pion].GetZ() << "\n";
        std::cout << "---------FraMomX: " << FragmentTracks[temp_id_fragment].GetPx() << "\n";
        std::cout << "---------FraMomY: " << FragmentTracks[temp_id_fragment].GetPy() << "\n";
        std::cout << "---------FraMomZ: " << FragmentTracks[temp_id_fragment].GetPz() << "\n";
        std::cout << "---------PioMomX: " << PionTracks[temp_id_pion].GetPx() << "\n";
        std::cout << "---------PioMomY: " << PionTracks[temp_id_pion].GetPy() << "\n";
        std::cout << "---------PioMomZ: " << PionTracks[temp_id_pion].GetPz() << "\n";
        std::cout << "---------NFragms: " << FragmentTracks.size() << "\n";
        std::cout << "---------NPions : " << PionTracks.size() << "\n";
        std::cout << "---------FragInd: " << temp_id_fragment << "\n";
        std::cout << "---------PionInd: " << temp_id_pion << "\n";
      }
*/
      LocalHisto.h_HypInvariantMass[pion_type]->Fill(MotherTracks_PrimVtx[i].GetMass(), 1.);
      LocalHisto.h_HypErrorInvariantMass[pion_type]->Fill(MotherTracks_PrimVtx[i].GetErrMass(), 1.);

      if((MotherTracks_PrimVtx[i].GetZ() - att.Target_PositionZ >=  5.)
          && (MotherTracks_PrimVtx[i].GetZ() - att.Target_PositionZ <= 27.))
            LocalHisto.h_HypInvariantMass_Z05[pion_type]->Fill(MotherTracks_PrimVtx[i].GetMass(), 1.);
      if((MotherTracks_PrimVtx[i].GetZ() - att.Target_PositionZ >= 10.)
          && (MotherTracks_PrimVtx[i].GetZ() - att.Target_PositionZ <= 27.))
            LocalHisto.h_HypInvariantMass_Z10[pion_type]->Fill(MotherTracks_PrimVtx[i].GetMass(), 1.);
      if((MotherTracks_PrimVtx[i].GetZ() - att.Target_PositionZ >= 15.)
          && (MotherTracks_PrimVtx[i].GetZ() - att.Target_PositionZ <= 27.))
            LocalHisto.h_HypInvariantMass_Z15[pion_type]->Fill(MotherTracks_PrimVtx[i].GetMass(), 1.);
      if((MotherTracks_PrimVtx[i].GetZ() - att.Target_PositionZ >= 20.)
          && (MotherTracks_PrimVtx[i].GetZ() - att.Target_PositionZ <= 27.))
            LocalHisto.h_HypInvariantMass_Z20[pion_type]->Fill(MotherTracks_PrimVtx[i].GetMass(), 1.);

      float m_getmass;
      float error_getmass;
      int status_getmass = MotherTracks_PrimVtx[i].GetMass(m_getmass, error_getmass);

      LocalHisto.h_HypInvariantMassCheck[pion_type]->Fill(m_getmass, status_getmass, 1.);
      LocalHisto.h_HypInvariantErrorMassCheck[pion_type]->Fill(error_getmass, status_getmass, 1.);


      double TimeLife_PrimVtx = MotherTracks_PrimVtx[i].GetLifeTime() / c_light_speed_cmps; //in ps
      LocalHisto.h_HypLifeTime_PrimVtx[pion_type]->Fill(TimeLife_PrimVtx, 1.);

      double ErrorTimeLife_PrimVtx = MotherTracks_PrimVtx[i].GetErrLifeTime() / c_light_speed_cmps; //in ps
      LocalHisto.h_HypErrorLifeTime_PrimVtx[pion_type]->Fill(ErrorTimeLife_PrimVtx, 1.);

      LocalHisto.h_Hyp_ArmenterosPodolanski[pion_type]->Fill(armenterosQtAlfa[1], armenterosQtAlfa[0], 1.);
    }
  
  std::vector<KFParticle> MotherTracks_PrimVtx_Mass_All;
  std::vector<std::tuple<size_t, size_t>> RefDaughtersTracks_PrimVtx_Mass_All;
  ifSet_ProductionVertex = true;
  ifSet_MassConstraint = true;
  MotherTracksRecons(FragmentTracks, PionTracks, pointer_PrimVtx, MotherTracks_PrimVtx_Mass_All, RefDaughtersTracks_PrimVtx_Mass_All);

  std::vector<KFParticle> MotherTracks_PrimVtx_Mass;
  std::vector<std::tuple<size_t, size_t>> RefDaughtersTracks_PrimVtx_Mass;
  MotherSelector(MotherTracks_PrimVtx_Mass_All, RefDaughtersTracks_PrimVtx_Mass_All, FragmentTracks, PionTracks, RecoEvent.PrimVtxRecons,
                 MotherTracks_PrimVtx_Mass, RefDaughtersTracks_PrimVtx_Mass);

  float MotherMass_InvMass;
  float MotherMass_ErrInvMass;
  float MotherMass_LifeTime;
  float MotherMass_ErrLifeTime;

  double dist_DaughtersMass;
  TVector3 centroid_DaughtersMass;

  for(size_t i = 0; i < MotherTracks_PrimVtx_Mass.size(); ++i)
    {
      Hyp temp_Hyp_Mass;

      temp_Hyp_Mass.Pattern = 5;

      temp_Hyp_Mass.PDG = MotherTracks_PrimVtx_Mass[i].GetPDG();
      temp_Hyp_Mass.N_Mother = MotherTracks_PrimVtx_Mass.size();
      temp_Hyp_Mass.Chi2ndf = MotherTracks_PrimVtx_Mass[i].GetChi2() / static_cast<double>(MotherTracks_PrimVtx_Mass[i].GetNDF());
      temp_Hyp_Mass.NDF = MotherTracks_PrimVtx_Mass[i].GetNDF();
      temp_Hyp_Mass.MomE.SetPxPyPzE(MotherTracks_PrimVtx_Mass[i].GetPx(),MotherTracks_PrimVtx_Mass[i].GetPy(),MotherTracks_PrimVtx_Mass[i].GetPz(),MotherTracks_PrimVtx_Mass[i].GetE());
      MotherTracks_PrimVtx_Mass[i].TransportToProductionVertex();
      temp_Hyp_Mass.PrimVtx.SetXYZ(MotherTracks_PrimVtx_Mass[i].GetX(),MotherTracks_PrimVtx_Mass[i].GetY(),MotherTracks_PrimVtx_Mass[i].GetZ());
      MotherTracks_PrimVtx_Mass[i].TransportToDecayVertex();
      temp_Hyp_Mass.DecayVtx.SetXYZ(MotherTracks_PrimVtx_Mass[i].GetX(),MotherTracks_PrimVtx_Mass[i].GetY(),MotherTracks_PrimVtx_Mass[i].GetZ());
      temp_Hyp_Mass.Dist_RealReconsVtx.SetXYZ(DecayVertex_real.X()-MotherTracks_PrimVtx_Mass[i].GetX(),DecayVertex_real.Y()-MotherTracks_PrimVtx_Mass[i].GetY(),DecayVertex_real.Z()-MotherTracks_PrimVtx_Mass[i].GetZ());
      ThetaDist_TrackPrimVtx(MotherTracks_PrimVtx_Mass[i], InteractionPoint_real, theta_DecayTrackPrimVtx, dist_DecayTrackPrimVtx);
      temp_Hyp_Mass.Dist_MotherPrimVtx = dist_DecayTrackPrimVtx;
      temp_Hyp_Mass.Angle_MotherPrimVtx = theta_DecayTrackPrimVtx;
      temp_Hyp_Mass.ErrGetMass = MotherTracks_PrimVtx_Mass[i].GetMass(MotherMass_InvMass, MotherMass_ErrInvMass);
      temp_Hyp_Mass.InvMass = MotherMass_InvMass;
      temp_Hyp_Mass.ErrInvMass = MotherMass_ErrInvMass;
      temp_Hyp_Mass.ErrGetLifeTime = MotherTracks_PrimVtx_Mass[i].GetLifeTime(MotherMass_LifeTime, MotherMass_ErrLifeTime);
      temp_Hyp_Mass.LifeTime = MotherMass_LifeTime / c_light_speed_cmps;
      temp_Hyp_Mass.ErrLifeTime = MotherMass_ErrLifeTime / c_light_speed_cmps;

      TVector3 Mother_Mom(MotherTracks_PrimVtx_Mass[i].GetPx(), MotherTracks_PrimVtx_Mass[i].GetPy(), MotherTracks_PrimVtx_Mass[i].GetPz());

      size_t temp_id_fragment = get<0>(RefDaughtersTracks_PrimVtx_Mass[i]);
      size_t temp_id_pion = get<1>(RefDaughtersTracks_PrimVtx_Mass[i]);

      temp_Hyp_Mass.PDG_Fragment = FragmentTracks[temp_id_fragment].GetPDG();
      temp_Hyp_Mass.Id_Fragment = FragmentTracks[temp_id_fragment].Id();
      temp_Hyp_Mass.MomE_Fragment.SetPxPyPzE(FragmentTracks[temp_id_fragment].GetPx(), FragmentTracks[temp_id_fragment].GetPy(), FragmentTracks[temp_id_fragment].GetPz(), FragmentTracks[temp_id_fragment].GetE());
      temp_Hyp_Mass.Chi2ndf_Fragment = FragmentTracks[temp_id_fragment].GetChi2() / static_cast<double>(FragmentTracks[temp_id_fragment].GetNDF());      
      temp_Hyp_Mass.NDF_Fragment = FragmentTracks[temp_id_fragment].GetNDF();
      temp_Hyp_Mass.Pvalue_Fragment = Fragment_FitInfo[temp_id_fragment].Pvalue;
      TVector3 Fragment_Mom(FragmentTracks[temp_id_fragment].GetPx(), FragmentTracks[temp_id_fragment].GetPy(), FragmentTracks[temp_id_fragment].GetPz());
      temp_Hyp_Mass.Angle_MotherFragment = Mother_Mom.Angle(Fragment_Mom) * 180. / M_PI;
      temp_Hyp_Mass.Fragment_IsFromHyp = 0;
      if(FragmentTracks[temp_id_fragment].Id() == FragmentTracks[ref_RealFragment].Id())
        temp_Hyp_Mass.Fragment_IsFromHyp = 1;

      temp_Hyp_Mass.PDG_Pion = PionTracks[temp_id_pion].GetPDG();
      temp_Hyp_Mass.Id_Pion = PionTracks[temp_id_pion].Id();
      temp_Hyp_Mass.MomE_Pion.SetPxPyPzE(PionTracks[temp_id_pion].GetPx(), PionTracks[temp_id_pion].GetPy(), PionTracks[temp_id_pion].GetPz(), PionTracks[temp_id_pion].GetE());
      temp_Hyp_Mass.Chi2ndf_Pion = PionTracks[temp_id_pion].GetChi2() / static_cast<double>(PionTracks[temp_id_pion].GetNDF());
      temp_Hyp_Mass.NDF_Pion = PionTracks[temp_id_pion].GetNDF();
      temp_Hyp_Mass.Pvalue_Pion = Pion_FitInfo[temp_id_pion].Pvalue;
      TVector3 Pion_Mom(PionTracks[temp_id_pion].GetPx(), PionTracks[temp_id_pion].GetPy(), PionTracks[temp_id_pion].GetPz());
      temp_Hyp_Mass.Angle_MotherPion = Mother_Mom.Angle(Pion_Mom) * 180. / M_PI;
      temp_Hyp_Mass.NHitsMDC_Pion = Pion_FitInfo[temp_id_pion].NHitsMDC;
      temp_Hyp_Mass.NHitsMinifiber_Pion = Pion_FitInfo[temp_id_pion].NHitsMinifiber;
      temp_Hyp_Mass.N_Pion = PionTracks.size();
      temp_Hyp_Mass.Pion_IsFromHyp = 0;

      for(itr_recons = RecoEvent.DaughtersTrackDAFInit.begin(); itr_recons != RecoEvent.DaughtersTrackDAFInit.end(); ++itr_recons)
        {
          if(itr_recons->first == PionTracks[temp_id_pion].Id())
            temp_Hyp_Mass.Pion_IsFromHyp = 1;
        }

      temp_Hyp_Mass.Mother_IsFromHyp = temp_Hyp_Mass.Fragment_IsFromHyp * temp_Hyp_Mass.Pion_IsFromHyp;

      float Vertex[3] = {MotherTracks_PrimVtx_Mass[i].GetX(), MotherTracks_PrimVtx_Mass[i].GetY(), MotherTracks_PrimVtx_Mass[i].GetZ()};
      FragmentTracks[temp_id_fragment].TransportToPoint(Vertex);
      PionTracks[temp_id_pion].TransportToPoint(Vertex);

      CloseDist(FragmentTracks[temp_id_fragment], PionTracks[temp_id_pion], dist_DaughtersMass, centroid_DaughtersMass);
      temp_Hyp_Mass.Dist_Daughters = dist_DaughtersMass;

      float armenterosQtAlfa[2] = {0.};
      MotherTracks_PrimVtx_Mass[i].GetArmenterosPodolanski_FromMother(FragmentTracks[temp_id_fragment], PionTracks[temp_id_pion], armenterosQtAlfa);
      temp_Hyp_Mass.ArmPod_Qt = armenterosQtAlfa[0];
      temp_Hyp_Mass.ArmPod_Alfa = armenterosQtAlfa[1];

      RecoEvent.Hyp_Vect.emplace_back(temp_Hyp_Mass);


      double distance_KFPart_PrimVtx_Mass  = std::sqrt(std::pow((DecayVertex_real.X() - MotherTracks_PrimVtx_Mass[i].GetX()), 2.) +
                              std::pow((DecayVertex_real.Y() - MotherTracks_PrimVtx_Mass[i].GetY()), 2.) +
                              std::pow((DecayVertex_real.Z() - MotherTracks_PrimVtx_Mass[i].GetZ()), 2.));
      double distanceX_KFPart_PrimVtx_Mass = DecayVertex_real.X() - MotherTracks_PrimVtx_Mass[i].GetX();
      double distanceY_KFPart_PrimVtx_Mass = DecayVertex_real.Y() - MotherTracks_PrimVtx_Mass[i].GetY();
      double distanceZ_KFPart_PrimVtx_Mass = DecayVertex_real.Z() - MotherTracks_PrimVtx_Mass[i].GetZ();

      LocalHisto.h_DecayVertexDistance_KFPart_PrimVtx_Mass[pion_type]->Fill(distance_KFPart_PrimVtx_Mass, 1.);
      LocalHisto.h_DecayVertexDistanceX_KFPart_PrimVtx_Mass[pion_type]->Fill(distanceX_KFPart_PrimVtx_Mass, 1.);
      LocalHisto.h_DecayVertexDistanceY_KFPart_PrimVtx_Mass[pion_type]->Fill(distanceY_KFPart_PrimVtx_Mass, 1.);
      LocalHisto.h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass[pion_type]->Fill(distanceZ_KFPart_PrimVtx_Mass, 1.);
    }



  if(MotherTracks_PrimVtx.size() == 0)
    return -5;

/*
  //Silicon HITS study of reconstructed mother and daughters
  std::vector<std::tuple<int, std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>>>> Fragment_SiHits;
  std::vector<std::tuple<int, std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>>>> Pion_SiHits;
  std::vector<std::tuple<int, std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>>>> Mother_SiHits;

  MotherDaughtersTrack_SiHits(FragmentTracks, PionTracks, MotherTracks_PrimVtx, RefDaughtersTracks_PrimVtx,
                              RecoEvent.Hits_Si1, RecoEvent.Hits_Si2, Fragment_SiHits, Pion_SiHits, Mother_SiHits);

  for(size_t i = 0; i < Fragment_SiHits.size(); ++i)
    {
      int temp_Nhits_Si1 = get<0>(get<1>(Fragment_SiHits[i])).size();
      int temp_Nhits_Si2 = get<1>(get<1>(Fragment_SiHits[i])).size();

      for(size_t k = 0; k < 2; ++k)
        {
          if(get<k>(get<1>(Fragment_SiHits[i]))[0][0] == -1.)
            {
              if(k == 0)
                temp_Nhits_Si1 -= 1;
              else if(k == 1)
                temp_Nhits_Si2 -= 1;
            }
        }

      LocalHisto.h_N_SiHits_ReconsTracks[pion_type]->Fill(temp_Nhits_Si1, "Fragment_Si1", 1.);
      LocalHisto.h_N_SiHits_ReconsTracks[pion_type]->Fill(temp_Nhits_Si2, "Fragment_Si2", 1.);
    }

  for(size_t i = 0; i < Pion_SiHits.size(); ++i)
    {
      int temp_Nhits_Si1 = get<0>(get<1>(Pion_SiHits[i])).size();
      int temp_Nhits_Si2 = get<1>(get<1>(Pion_SiHits[i])).size();

      for(size_t k = 0; k < 2; ++k)
        {
          if(get<k>(get<1>(Pion_SiHits[i]))[0][0] == -1.)
            {
              if(k == 0)
                temp_Nhits_Si1 -= 1;
              else if(k == 1)
                temp_Nhits_Si2 -= 1;
            }
        }

      LocalHisto.h_N_SiHits_ReconsTracks[pion_type]->Fill(temp_Nhits_Si1, "Pion_Si1", 1.);
      LocalHisto.h_N_SiHits_ReconsTracks[pion_type]->Fill(temp_Nhits_Si2, "Pion_Si2", 1.);
    }

  for(size_t i = 0; i < Mother_SiHits.size(); ++i)
    {
      int temp_Nhits_Si1 = get<0>(get<1>(Mother_SiHits[i])).size();
      int temp_Nhits_Si2 = get<1>(get<1>(Mother_SiHits[i])).size();

      for(size_t k = 0; k < 2; ++k)
        {
          if(get<k>(get<1>(Mother_SiHits[i]))[0][0] == -1.)
            {
              if(k == 0)
                temp_Nhits_Si1 -= 1;
              else if(k == 1)
                temp_Nhits_Si2 -= 1;
            }
        }

      LocalHisto.h_N_SiHits_ReconsTracks[pion_type]->Fill(temp_Nhits_Si1, "Mother_Si1", 1.);
      LocalHisto.h_N_SiHits_ReconsTracks[pion_type]->Fill(temp_Nhits_Si2, "Mother_Si2", 1.);
    }
*/

  //Silicon SIGNALS study of reconstructed mother and daughters
/*
  std::vector<std::tuple<int, std::vector<std::vector<std::vector<double>>>>> Fragment_SiHits {};
  std::vector<std::tuple<int, std::vector<std::vector<std::vector<double>>>>> Pion_SiHits {};
  std::vector<std::tuple<int, std::vector<std::vector<std::vector<double>>>>> Mother_SiHits {};

  MotherDaughtersTrack_SiHits2(FragmentTracks, PionTracks, MotherTracks_PrimVtx, RefDaughtersTracks_PrimVtx,
                                RecoEvent.HitsX_Si1, RecoEvent.HitsY_Si1, RecoEvent.HitsX_Si2, RecoEvent.HitsY_Si2,
                                Fragment_SiHits, Pion_SiHits, Mother_SiHits);

  for(size_t i = 0; i < Fragment_SiHits.size(); ++i)
    {
      double temp_Nhits_Si1 = (get<1>(Fragment_SiHits[i])[0].size() + get<1>(Fragment_SiHits[i])[1].size()) / 2.;
      double temp_Nhits_Si2 = (get<1>(Fragment_SiHits[i])[2].size() + get<1>(Fragment_SiHits[i])[3].size()) / 2.;

      for(size_t k = 0; k < 4; ++k)
        {
          if(get<1>(Fragment_SiHits[i])[k][0][0] == -1.)
            {
              if((k == 0) || (k == 1))
                temp_Nhits_Si1 -= 0.5;
              else if((k == 2) || (k == 3))
                temp_Nhits_Si2 -= 0.5;
            }
        }

      LocalHisto.h_N_SiHits_ReconsTracks[pion_type]->Fill(temp_Nhits_Si1, "Fragment_Si1", 1.);
      LocalHisto.h_N_SiHits_ReconsTracks[pion_type]->Fill(temp_Nhits_Si2, "Fragment_Si2", 1.);
    }

  for(size_t i = 0; i < Pion_SiHits.size(); ++i)
    {
      double temp_Nhits_Si1 = (get<1>(Pion_SiHits[i])[0].size() + get<1>(Pion_SiHits[i])[1].size()) / 2.;
      double temp_Nhits_Si2 = (get<1>(Pion_SiHits[i])[2].size() + get<1>(Pion_SiHits[i])[3].size()) / 2.;

      for(size_t k = 0; k < 4; ++k)
        {
          if(get<1>(Pion_SiHits[i])[k][0][0] == -1.)
            {
              if((k == 0) || (k == 1))
                temp_Nhits_Si1 -= 0.5;
              else if((k == 2) || (k == 3))
                temp_Nhits_Si2 -= 0.5;
            }
        }

      LocalHisto.h_N_SiHits_ReconsTracks[pion_type]->Fill(temp_Nhits_Si1, "Pion_Si1", 1.);
      LocalHisto.h_N_SiHits_ReconsTracks[pion_type]->Fill(temp_Nhits_Si2, "Pion_Si2", 1.);
    }

  for(size_t i = 0; i < Mother_SiHits.size(); ++i)
    {
      double temp_Nhits_Si1 = (get<1>(Mother_SiHits[i])[0].size() + get<1>(Mother_SiHits[i])[1].size()) / 2.;
      double temp_Nhits_Si2 = (get<1>(Mother_SiHits[i])[2].size() + get<1>(Mother_SiHits[i])[3].size()) / 2.;

      for(size_t k = 0; k < 4; ++k)
        {
          if(get<1>(Mother_SiHits[i])[k][0][0] == -1.)
            {
              if((k == 0) || (k == 1))
                temp_Nhits_Si1 -= 0.5;
              else if((k == 2) || (k == 3))
                temp_Nhits_Si2 -= 0.5;
            }
        }

      LocalHisto.h_N_SiHits_ReconsTracks[pion_type]->Fill(temp_Nhits_Si1, "Mother_Si1", 1.);
      LocalHisto.h_N_SiHits_ReconsTracks[pion_type]->Fill(temp_Nhits_Si2, "Mother_Si2", 1.);
    }
*/

/*    NOT USEFUL, NOT WORKING
      LocalHisto.h_N_Si_MotherTracks[pion_type]->Fill(0.5, 1.);

      KFParticle Si_MotherTrack;
      MotherTrackSiliconHits(RecoEvent.PrimVtxRecons, DecayVertexRecons, RecoEvent.Hits_Si1, RecoEvent.Hits_Si2, MotherTracks_PrimVtx[0], Si_MotherTrack);
      
      if(Si_MotherTrack.GetP() < 1.e-4)
        return -6;
                  
      LocalHisto.h_N_Si_MotherTracks[pion_type]->Fill(2.5, 1.);

      std::vector<KFParticle> AllTracks;
      AllTracks.emplace_back(PionTracks[0]);
      AllTracks.emplace_back(FragmentTracks[0]);
      AllTracks.emplace_back(Si_MotherTrack);

      TVector3 AllVfunc_DecayVertexRecons;
      AllTrackstoDecayVertex_Vfunction(AllTracks, closedist_pos, AllVfunc_DecayVertexRecons);

      double AllVfunc_distance  = std::sqrt(std::pow((DecayVertex_real.X() - AllVfunc_DecayVertexRecons.X()), 2.) +
                                         std::pow((DecayVertex_real.Y() - AllVfunc_DecayVertexRecons.Y()), 2.) +
                                          std::pow((DecayVertex_real.Z() - AllVfunc_DecayVertexRecons.Z()), 2.));
      double AllVfunc_distanceX = DecayVertex_real.X() - AllVfunc_DecayVertexRecons.X();
      double AllVfunc_distanceY = DecayVertex_real.Y() - AllVfunc_DecayVertexRecons.Y();
      double AllVfunc_distanceZ = DecayVertex_real.Z() - AllVfunc_DecayVertexRecons.Z();

      LocalHisto.h_DecayVertexDistance_AllVfunc[pion_type]->Fill(AllVfunc_distance, 1.);
      LocalHisto.h_DecayVertexDistanceX_AllVfunc[pion_type]->Fill(AllVfunc_distanceX, 1.);
      LocalHisto.h_DecayVertexDistanceY_AllVfunc[pion_type]->Fill(AllVfunc_distanceY, 1.);
      LocalHisto.h_DecayVertexDistanceZ_AllVfunc[pion_type]->Fill(AllVfunc_distanceZ, 1.);

      LocalHisto.h_DecayVertexPosZ_AllVfunc[pion_type]->Fill(AllVfunc_DecayVertexRecons.Z() - att.Target_PositionZ, 1.);



      TVector3 AllCentroid_DecayVertexRecons;
      AllTrackstoDecayVertex_Centroids(AllTracks, AllCentroid_DecayVertexRecons);

      double AllCentroid_distance  = std::sqrt(std::pow((DecayVertex_real.X() - AllCentroid_DecayVertexRecons.X()), 2.) +
                                            std::pow((DecayVertex_real.Y() - AllCentroid_DecayVertexRecons.Y()), 2.) +
                                             std::pow((DecayVertex_real.Z() - AllCentroid_DecayVertexRecons.Z()), 2.));
      double AllCentroid_distanceX = DecayVertex_real.X() - AllCentroid_DecayVertexRecons.X();
      double AllCentroid_distanceY = DecayVertex_real.Y() - AllCentroid_DecayVertexRecons.Y();
      double AllCentroid_distanceZ = DecayVertex_real.Z() - AllCentroid_DecayVertexRecons.Z();

      LocalHisto.h_DecayVertexDistance_AllCentroid[pion_type]->Fill(AllCentroid_distance, 1.);
      LocalHisto.h_DecayVertexDistanceX_AllCentroid[pion_type]->Fill(AllCentroid_distanceX, 1.);
      LocalHisto.h_DecayVertexDistanceY_AllCentroid[pion_type]->Fill(AllCentroid_distanceY, 1.);
      LocalHisto.h_DecayVertexDistanceZ_AllCentroid[pion_type]->Fill(AllCentroid_distanceZ, 1.);

      LocalHisto.h_DecayVertexPosZ_AllCentroid[pion_type]->Fill(AllCentroid_DecayVertexRecons.Z() - att.Target_PositionZ, 1.);



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

      double AllKFPart_distance  = std::sqrt(std::pow((DecayVertex_real.X() - AllKFPart_DecayVertexRecons.X()), 2.) +
                                            std::pow((DecayVertex_real.Y() - AllKFPart_DecayVertexRecons.Y()), 2.) +
                                             std::pow((DecayVertex_real.Z() - AllKFPart_DecayVertexRecons.Z()), 2.));
      double AllKFPart_distanceX = DecayVertex_real.X() - AllKFPart_DecayVertexRecons.X();
      double AllKFPart_distanceY = DecayVertex_real.Y() - AllKFPart_DecayVertexRecons.Y();
      double AllKFPart_distanceZ = DecayVertex_real.Z() - AllKFPart_DecayVertexRecons.Z();

      LocalHisto.h_DecayVertexDistance_AllKFPart[pion_type]->Fill(AllKFPart_distance, 1.);
      LocalHisto.h_DecayVertexDistanceX_AllKFPart[pion_type]->Fill(AllKFPart_distanceX, 1.);
      LocalHisto.h_DecayVertexDistanceY_AllKFPart[pion_type]->Fill(AllKFPart_distanceY, 1.);
      LocalHisto.h_DecayVertexDistanceZ_AllKFPart[pion_type]->Fill(AllKFPart_distanceZ, 1.);

      LocalHisto.h_DecayVertexPosZ_AllKFPart[pion_type]->Fill(AllKFPart_DecayVertexRecons.Z() - att.Target_PositionZ, 1.);
*/


  return 0;
}



template <class Out>
void TDecayVertex<Out>::StudyCaseSelector_Hyp(std::string StudyCase)
{
  if(StudyCase.compare("H3L") == 0)
    {
      Hyp_pdg = H3L_pdg;
      recons_from_FRS_MDC = 1;
    }
  else if(StudyCase.compare("H4L") == 0)
    {
      Hyp_pdg = H4L_pdg;
      recons_from_FRS_MDC = 1;
    }
  else if(StudyCase.compare("nnL") == 0)
    {
      Hyp_pdg = nnL_pdg;
      recons_from_FRS_MDC = 1;
    }
  else if(StudyCase.compare("lambda") == 0)
    {
      Hyp_pdg = lambda_pdg;
      //FragmentPDG = proton_pdg;
      recons_from_FRS_MDC = 2;
    }
  else if(StudyCase.compare("background_H3L") == 0)
    {
      Hyp_pdg = H3L_pdg;
      recons_from_FRS_MDC = 1;
    }
  else if(StudyCase.compare("background_H4L") == 0)
    {
      Hyp_pdg = H4L_pdg;
      recons_from_FRS_MDC = 1;
    }

  return;
}


template <class Out>
void TDecayVertex<Out>::RealTracksFinder(std::unordered_map<int, std::vector<std::vector<SimHit> > >& TrackDAFSim,
                                        int& pdgParticle, int& cutConditions,
                                        std::vector<KFParticle>& RealTracks,
                                        std::vector<KFFitInfo>& Vect_FitInfo)
{
  std::unordered_map<int, std::vector<std::vector<SimHit> > >::iterator itr;
  for(itr = TrackDAFSim.begin(); itr != TrackDAFSim.end(); ++itr)
    {
      int iDetFirst = -1;
      
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
          if(iDet >= G4Sol::MiniFiberD1_x && iDet <= G4Sol::MiniFiberD2_u)
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

          double temp_fC[] = {1.e-1,
                                 0., 1.e-1,
                                 0.,    0., 1.e-1,
                                 0.,    0.,    0., 1.e-2,
                                 0.,    0.,    0.,    0., 1.e-2,
                                 0.,    0.,    0.,    0.,    0., 1.e-2}; //Change!

          int temp_charge = TDatabasePDG::Instance()->GetParticle(pdgParticle)->Charge()/3.;
          double temp_mass = TDatabasePDG::Instance()->GetParticle(pdgParticle)->Mass();

          KFParticle temp_particle;
          temp_particle.Create(temp_fP, temp_fC, temp_charge, temp_mass);
          temp_particle.SetId(itr->first);

          temp_particle.SetNDF(5); //Change !
          temp_particle.SetChi2(5.); // Change !

          temp_particle.SetField(att.Field);

          KFFitInfo temp_FitInfo;
          temp_FitInfo.Pvalue = 0.; //Change!
          temp_FitInfo.NHitsMDC = nHits_MDC;
          temp_FitInfo.NHitsMinifiber = nHits_MiniFiber;

          if((cutConditions == 1) && ((nHits_MDC < 6) || (nHits_MiniFiber < 3) || ((nHits_PSCE == 0) && (nHits_PSBE == 0))))
            continue;

          RealTracks.emplace_back(temp_particle);
          Vect_FitInfo.emplace_back(temp_FitInfo);
        }
    }

  return;
}


template <class Out>
void TDecayVertex<Out>::FragmentMDCTracksFinder(std::unordered_map<int, ResSolDAF>& DAF_results, int& fragment_pdg,
                                                    std::vector<FragmentTrack>& FragmentMDCTracks)
{
  //int temp_charge = TDatabasePDG::Instance()->GetParticle(fragment_pdg)->Charge()/3.;

  std::unordered_map<int, ResSolDAF>::iterator itr;
  for(itr = DAF_results.begin(); itr != DAF_results.end(); ++itr)
    {
      if((itr->second.pdg_guess == fragment_pdg) && (itr->second.Ncentral > att.KF_NbCentralCut))
        {
          TVector3 tmp_pos = TVector3(itr->second.posX, itr->second.posY, itr->second.posZ);
          TVector3 tmp_mom = TVector3(itr->second.momX, itr->second.momY, itr->second.momZ);
          std::vector<double> tmp_covmatrix = {1.e-8,
                                                  0., 1.e-8,
                                                  0.,    0., 1.e-8,
                                                  0.,    0.,    0., 1.e-8,
                                                  0.,    0.,    0.,    0., 1.e-8,
                                                  0.,    0.,    0.,    0.,    0., 1.e-8}; //Change!

          FragmentTrack tmp_FragmentTrack;
          tmp_FragmentTrack.SetPos(tmp_pos);
          tmp_FragmentTrack.SetMom(tmp_mom);
          tmp_FragmentTrack.SetCovMatrix(tmp_covmatrix);
          tmp_FragmentTrack.SetTOT(1.); //CHECK Change!
          tmp_FragmentTrack.SetTime(1.); //CHECK Change!
          tmp_FragmentTrack.SetChi2NDF(1.); //CHECK Change!
          tmp_FragmentTrack.SetIsBest(true);
          tmp_FragmentTrack.SetPID(fragment_pdg);

          FragmentMDCTracks.emplace_back(tmp_FragmentTrack);
        }
    }

  return;
}

/*
template <class Out>
void TDecayVertex<Out>::FragmentMDCTracksFinder(std::unordered_map<int, ResSolDAF>& DAF_results, int& fragment_pdg,
                                    std::vector<KFParticle>& FragmentMDCTracks)
{
  int temp_charge = TDatabasePDG::Instance()->GetParticle(fragment_pdg)->Charge()/3.;
  double temp_mass = TDatabasePDG::Instance()->GetParticle(fragment_pdg)->Mass();

  std::unordered_map<int, ResSolDAF>::iterator itr;
  for(itr = DAF_results.begin(); itr != DAF_results.end(); ++itr)
    {
      if((itr->second.pdg_guess == fragment_pdg) && (itr->second.Ncentral > att.KF_NbCentralCut))
        {
          double temp_fP[] = {itr->second.posX, itr->second.posY, itr->second.posZ, // X, Y, Z
                              itr->second.momX, itr->second.momY, itr->second.momZ}; // Px, Py, Pz

          double temp_fC[] = {itr->second.cov_matrix[0][0],
                              itr->second.cov_matrix[1][0], itr->second.cov_matrix[1][1],
                              itr->second.cov_matrix[2][0], itr->second.cov_matrix[2][1], itr->second.cov_matrix[2][2],
                              itr->second.cov_matrix[3][0], itr->second.cov_matrix[3][1], itr->second.cov_matrix[3][2], itr->second.cov_matrix[3][3],
                              itr->second.cov_matrix[4][0], itr->second.cov_matrix[4][1], itr->second.cov_matrix[4][2], itr->second.cov_matrix[4][3], itr->second.cov_matrix[4][4],
                              itr->second.cov_matrix[5][0], itr->second.cov_matrix[5][1], itr->second.cov_matrix[5][2], itr->second.cov_matrix[5][3], itr->second.cov_matrix[5][4], itr->second.cov_matrix[5][5]};

          KFParticle temp_particle;

          temp_particle.Create(temp_fP, temp_fC, temp_charge, temp_mass);
          temp_particle.SetId(itr->first);

          temp_particle.SetField(att.Field);
          temp_particle.SetNDF(itr->second.ndf);
          temp_particle.SetChi2(itr->second.chi2);

          FragmentMDCTracks.emplace_back(temp_particle);
        }
    }

  return;
}
*/

template <class Out>
void TDecayVertex<Out>::FragmentSelector(std::vector<KFParticle>& FragmentTracks_All, TVector3& PrimVtxRecons, std::vector<KFParticle>& FragmentTracks)
{
/*
  Possible cuts:
  - upper limit to chi2/ndf
  - lower limit to d with IP
*/
  if(FragmentTracks_All.size() == 0)
    {
      std::cout << "Error: N_FragmentTracks_All=0 in FragmentSelector\n";
      return;
    }

  double temp_chi2ndf;
  double temp_theta;
  double temp_dist;


  for(size_t i = 0; i < FragmentTracks_All.size(); ++i)
    {
      temp_chi2ndf = FragmentTracks_All[i].GetChi2() / static_cast<double>(FragmentTracks_All[i].GetNDF());
      if( ifCut_MaxChi2ndf_FragmentTracks && (temp_chi2ndf > MaxChi2ndf_FragmentTracks) )
        continue;

      ThetaDist_TrackPrimVtx(FragmentTracks_All[i], PrimVtxRecons, temp_theta, temp_dist);
      if( ifCut_MinDist_FragmentTracksPrimVtx && (temp_dist < MinDist_FragmentTracksPrimVtx) )
        continue;

      if( ifCut_MinMomZ_FragmentTracks && (FragmentTracks_All[i].GetPz() < MinMomZ_FragmentTracks) )
        continue;

      double theta_FragmentTracks = acos(FragmentTracks_All[i].GetPz() / FragmentTracks_All[i].GetP()) * 180. / M_PI;
      if(ifCut_MaxTheta_FragmentTracks)
        if( ((recons_from_FRS_MDC == 1) && (theta_FragmentTracks > MaxTheta_FragmentTracks)) || ((recons_from_FRS_MDC == 2) && (theta_FragmentTracks > MaxTheta_FragmentMDCTracks)) )
          continue;

      FragmentTracks.emplace_back(FragmentTracks_All[i]);
    }

  return;
}


template <class Out>
void TDecayVertex<Out>::PionTracksFinder(std::unordered_map<int, ResSolDAF>& DAF_results,
                                          std::vector<KFParticle>& PionTracks,
                                          std::vector<KFFitInfo>& Vect_FitInfo)
{
  std::unordered_map<int, ResSolDAF>::iterator itr;
  for(itr = DAF_results.begin(); itr != DAF_results.end(); ++itr)
    {
      if((itr->second.pdg_guess == pi_pdg) && (itr->second.Ncentral > att.KF_NbCentralCut))
        {
          LocalHisto.h_Chi2ndf_pions[pion_type]->Fill(itr->second.chi2 / itr->second.ndf, 1.);

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
          temp_particle.SetPDG(pi_pdg);

          temp_particle.SetField(att.Field);

          temp_particle.SetNDF(itr->second.ndf);
          temp_particle.SetChi2(itr->second.chi2);

          PionTracks.emplace_back(temp_particle);

          KFFitInfo temp_FitInfo;
          temp_FitInfo.Pvalue = itr->second.pvalue;
          temp_FitInfo.NHitsMDC = itr->second.Ncentral;
          temp_FitInfo.NHitsMinifiber = itr->second.Nmfiber;

          Vect_FitInfo.emplace_back(temp_FitInfo);
        }
    }

  return;
}

template <class Out>
void TDecayVertex<Out>::PionSelector(std::vector<KFParticle>& PionTracks_All, TVector3& PrimVtxRecons, std::vector<KFParticle>& PionTracks, std::vector<KFFitInfo>& Vect_FitInfo)
{
/*
  Possible cuts:
  - upper limit to chi2/ndf
  - lower limit to d with IP
  - lower limit to Pz
*/
  if(PionTracks_All.size() == 0)
    {
      std::cout << "Error: N_PionTracks_All=0 in PionSelector\n";
      return;
    }

  double temp_chi2ndf;
  double temp_theta;
  double temp_dist;

  std::vector<KFFitInfo> copy_Vect_FitInfo;
  copy_Vect_FitInfo = Vect_FitInfo;
  Vect_FitInfo.clear();

  for(size_t i = 0; i < PionTracks_All.size(); ++i)
    {
      temp_chi2ndf = PionTracks_All[i].GetChi2() / static_cast<double>(PionTracks_All[i].GetNDF());
      if( ifCut_MaxChi2ndf_PionTracks && (temp_chi2ndf > MaxChi2ndf_PionTracks) )
        continue;

      ThetaDist_TrackPrimVtx(PionTracks_All[i], PrimVtxRecons, temp_theta, temp_dist);
      if( ifCut_MinDist_PionTracksPrimVtx && (temp_dist < MinDist_PionTracksPrimVtx) )
        continue;

      if( ifCut_MinMomZ_PionTracks && (PionTracks_All[i].GetPz() < MinMomZ_PionTracks) )
        continue;

      PionTracks.emplace_back(PionTracks_All[i]);
      Vect_FitInfo.emplace_back(copy_Vect_FitInfo[i]);
    }

  return;
}


template <class Out>
void TDecayVertex<Out>::CloseDist(KFParticle& FragmentTrack, KFParticle& PionTrack, double& distance, TVector3& centroid)
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

  return;
}

template <class Out>
double TDecayVertex<Out>::f_function(KFParticle& DecayTrack, TVector3& PosXYZ)
{
  double slope_x     = DecayTrack.GetPx() / DecayTrack.GetPz();
  double intercept_x = DecayTrack.GetX() - slope_x * DecayTrack.GetZ();
  double slope_y     = DecayTrack.GetPy() / DecayTrack.GetPz();
  double intercept_y = DecayTrack.GetY() - slope_y * DecayTrack.GetZ();

  double distanceStepX = 2. * boxDistXY / static_cast<double>(NstepsdiscretXY - 1); //Change(?)
  double sigma2        = std::pow(distanceStepX, 2.) / 12.;

  double f = exp(-0.5 *
                 (std::pow((PosXYZ.X() - slope_x * PosXYZ.Z() - intercept_x), 2.) +
                  std::pow((PosXYZ.Y() - slope_y * PosXYZ.Z() - intercept_y), 2.)) /
                 sigma2);
  return f;
}

template <class Out>
double TDecayVertex<Out>::V_function(std::vector<double>& f_vector)
{
  double sum_f  = 0;
  double sum_f2 = 0;
  double v      = 0.;

  for(size_t i = 1; i < f_vector.size(); ++i)
    {
      sum_f += f_vector[i];
      sum_f2 += std::pow(f_vector[i], 2.);
    }

  if((sum_f > 1.E-9) && (sum_f2 > 1.E-9))
      v = k_factor * f_vector[0] + sum_f -
          (k_factor * std::pow(f_vector[0], 2.) + sum_f2) / (k_factor * f_vector[0] + sum_f);

  return v;
}

template <class Out>
void TDecayVertex<Out>::SpaceDiscretization(double& Xi, double& Xf, size_t& NstepsX, double& Yi, double& Yf,
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

  return;
}

template <class Out>
void TDecayVertex<Out>::TrackstoDecayVertex(std::vector<KFParticle>& FragmentTracks, std::vector<KFParticle>& PionTracks,
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
              temp_f[j] = f_function(FragmentTracks[j], temp_PosXYZ);

          for(size_t j = 0; j < PionTracks.size(); ++j)
              temp_f[FragmentTracks.size()+j] = f_function(PionTracks[j], temp_PosXYZ);

          Vnew = V_function(temp_f);

          if(Vnew > V)
            {
              V = Vnew;
              DecayVertexRecons.SetXYZ(temp_PosXYZ.X(), temp_PosXYZ.Y(), temp_PosXYZ.Z());
            }
        }
    }

  return;
}

template <class Out>
void TDecayVertex<Out>::ThetaDist_TrackPrimVtx(KFParticle& Track, TVector3& PrimVtxRecons, double& theta, double& distance)
{
  TVector3 u(Track.GetPx(), Track.GetPy(), Track.GetPz());
  TVector3 P(Track.GetX(), Track.GetY(), Track.GetZ());
  TVector3 PA = P - PrimVtxRecons;

  theta = PA.Angle(u) * 180. / M_PI;
  distance = (PA.Cross(u)).Mag()/u.Mag();

  return;
}

template <class Out>
void TDecayVertex<Out>::KFPart_PrimaryVertex(TVector3& PrimVtxRecons, std::array<double,6> Cov_PrimVtx, KFParticleSIMD& temp_PrimVtx)
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
  temp_PrimVtx.SetField(fieldWASA);

  return;
}

/* //WORKING FINE FOR HOMOGENEUS FIELD
template <class Out>
void TDecayVertex<Out>::MotherTracksRecons(std::vector<KFParticle>& FragmentTracks, std::vector<KFParticle>& PionTracks,
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

          std::cout << "Fragment Mom: " << FragmentTracks[i].GetP() << "\n";
          std::cout << "Pion Mom: " << PionTracks[j].GetP() << "\n\n";


          particleSIMD1.GetDStoParticle( particleSIMD2, ds, dsdr );
          particleSIMD1.TransportToDS(ds[0], dsdr[0]);
          particleSIMD2.TransportToDS(ds[1], dsdr[3]);

          const KFParticleSIMD* vDaughtersPointer[2] = {&particleSIMD1, &particleSIMD2};

          KFParticleSIMD mother;
          mother.SetConstructMethod(KFPart_fConstructMethod);
          mother.Construct(vDaughtersPointer, 2);

          if( ifSet_ProductionVertex && (pointer_PrimVtx != nullptr) )
            mother.SetProductionVertex(*pointer_PrimVtx);

          if( ifSet_MassConstraint )
            mother.SetMassConstraint(Hyp_mass);

          mother.TransportToDecayVertex();

          KFParticle temp_MotherTrack;
          mother.GetKFParticle(temp_MotherTrack, 0);

          MotherTracks.emplace_back(temp_MotherTrack);
          RefDaughtersTracks.emplace_back(std::make_tuple(i,j));
        }
    }

  return;
}
*/

template <class Out>
void TDecayVertex<Out>::MotherTracksRecons(std::vector<KFParticle>& FragmentTracks, std::vector<KFParticle>& PionTracks,
                                      const KFParticleSIMD* pointer_PrimVtx, std::vector<KFParticle>& MotherTracks,
                                      std::vector<std::tuple<size_t, size_t>>& RefDaughtersTracks)
{
  float_v ds[2] = {0.f, 0.f};
  float_v dsdr[4][6];

  for(size_t i = 0; i < FragmentTracks.size(); ++i)
    {
      for(size_t j = 0; j < PionTracks.size(); ++j)
        {
          KFParticleSIMD particleSIMD1(FragmentTracks[i]);    // the same particle is copied to each SIMD element
          KFParticleSIMD particleSIMD2(PionTracks[j]);

          particleSIMD1.SetField(fieldWASA);
          particleSIMD2.SetField(fieldWASA);

          particleSIMD1.GetDStoParticle( particleSIMD2, ds, dsdr );
          particleSIMD1.TransportToDS(ds[0], dsdr[0]);
          particleSIMD2.TransportToDS(ds[1], dsdr[3]);

          const KFParticleSIMD* vDaughtersPointer[2] = {&particleSIMD1, &particleSIMD2};
/*
            std::cout << "M0FraPosX: " << particleSIMD1.GetX()    << "\n";
            std::cout << "M0FraPosY: " << particleSIMD1.GetY()    << "\n";
            std::cout << "M0FraPosZ: " << particleSIMD1.GetZ()    << "\n";
            std::cout << "M0FraMomX: " << particleSIMD1.GetPx()    << "\n";
            std::cout << "M0FraMomY: " << particleSIMD1.GetPy()    << "\n";
            std::cout << "M0FraMomZ: " << particleSIMD1.GetPz()    << "\n";

            std::cout << "M0PioPosX: " << particleSIMD2.GetX()    << "\n";
            std::cout << "M0PioPosY: " << particleSIMD2.GetY()    << "\n";
            std::cout << "M0PioPosZ: " << particleSIMD2.GetZ()    << "\n";
            std::cout << "M0PioMomX: " << particleSIMD2.GetPx()    << "\n";
            std::cout << "M0PioMomY: " << particleSIMD2.GetPy()    << "\n";
            std::cout << "M0PioMomZ: " << particleSIMD2.GetPz()    << "\n";
*/
          KFParticleSIMD mother;
          mother.SetConstructMethod(KFPart_fConstructMethod);
          mother.Construct(vDaughtersPointer, 2);
          mother.SetField(fieldWASA);
/*
          if(std::isnan(mother.GetMass()[0]))
          {
            std::cout << "M0DcyPosX: " << mother.GetX()    << "\n";
            std::cout << "M0DcyPosY: " << mother.GetY()    << "\n";
            std::cout << "M0DcyPosZ: " << mother.GetZ()    << "\n";
            std::cout << "M0SumDMas: " << mother.GetSumDaughterMass()    << "\n";
          }
*/
          if( ifSet_ProductionVertex && (pointer_PrimVtx != nullptr) )
            mother.SetProductionVertex(*pointer_PrimVtx);

          if( ifSet_MassConstraint )
            mother.SetMassConstraint(Hyp_mass);
/*
          if(std::isnan(mother.GetMass()[0]))
          {
            std::cout << "M1DcyPosX: " << mother.GetX()    << "\n";
            std::cout << "M1DcyPosY: " << mother.GetY()    << "\n";
            std::cout << "M1DcyPosZ: " << mother.GetZ()    << "\n";
            std::cout << "M1SumDMas: " << mother.GetSumDaughterMass()    << "\n";
          }
*/
          mother.TransportToDecayVertex();

          if(ifRemoveNaNMother && (std::isnan(mother.GetMass()[0])))
            continue;

          KFParticle temp_MotherTrack;
          mother.GetKFParticle(temp_MotherTrack, 0);

          temp_MotherTrack.SetField(att.Field);
          temp_MotherTrack.SetPDG(Hyp_pdg);

          MotherTracks.emplace_back(temp_MotherTrack);
          RefDaughtersTracks.emplace_back(std::make_tuple(i,j));
        }
    }

  return;
}

template <class Out>
void TDecayVertex<Out>::MotherSelector(std::vector<KFParticle>& MotherTracks_All, std::vector<std::tuple<size_t, size_t>>& RefDaughtersTracks_All,
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
      if( ifCut_MaxClosedist_DaughterTracks && (Closedist_DaughterTracks > MaxClosedist_DaughterTracks) )
        continue;

      angle_MotherFragment = MotherTracks_All[i].GetAngle(FragmentTracks[temp_id_fragment]) * 180. / M_PI;
      if( ifCut_MaxAngle_MotherFragment && (angle_MotherFragment > MaxAngle_MotherFragment) )
        continue;

      angle_MotherPion = MotherTracks_All[i].GetAngle(PionTracks[temp_id_pion]) * 180. / M_PI;
      if( ifCut_MaxAngle_MotherPion && (angle_MotherPion > MaxAngle_MotherPion) )
        continue;

      temp_chi2ndf = MotherTracks_All[i].GetChi2() / static_cast<double>(MotherTracks_All[i].GetNDF());
      if( ifCut_MaxChi2ndf && (temp_chi2ndf < MaxChi2ndf) )
        continue;

      ThetaDist_TrackPrimVtx(MotherTracks_All[i], PrimVtxRecons, Dist_MotherTrackPrimVtx, Angle_MotherTrackPrimVtx);
      if( ifCut_MaxDist_MotherTrackPrimVtx && (Dist_MotherTrackPrimVtx > MaxDist_MotherTrackPrimVtx) )
        continue;

      if( ifCut_MaxAngle_MotherTrackPrimVtx && (Angle_MotherTrackPrimVtx > MaxAngle_MotherTrackPrimVtx))
        continue;

      PosZ_DecayVertex = MotherTracks_All[i].GetZ() - att.Target_PositionZ;
      if( ifCut_MaxPosZ_DecayVertex && (PosZ_DecayVertex > MaxPosZ_DecayVertex) )
        continue;

      if( ifCut_MinPosZ_DecayVertex && (PosZ_DecayVertex < MinPosZ_DecayVertex) )
        continue;

      MotherTracks_All[i].GetArmenterosPodolanski_FromMother(FragmentTracks[temp_id_fragment], PionTracks[temp_id_pion], armenterosQtAlfa);
      if ( ifCut_ArmenterosPodolanski ) //Change !
        continue;


      MotherTracks.emplace_back(MotherTracks_All[i]);
      RefDaughtersTracks.emplace_back(RefDaughtersTracks_All[i]);
    }

  return;
}

/*
template <class Out>
void TDecayVertex<Out>::SiHitsFinder(KFParticle& Track, std::vector<std::vector<double> >& Hits_Si,
                                  std::vector<std::vector<double> >& Track_Sihit)
{
  if(Hits_Si.size() == 0)
    return;

  Track_Sihit.clear();
  Z_plane = Hits_Si[0][3];
  std::vector<double> Empty_Sihit = {-1., -1., -1., -1.}; // (E, X, Y, Z)

  double track_parameter = (Z_plane - Track.GetZ()) / Track.GetPz();
  float X_track = Track.GetX() + track_parameter * Track.GetPx();
  float Y_track = Track.GetY() + track_parameter * Track.GetPy();

  float TrackHit_Aprox[3] = {X_track, Y_track, Z_plane};
  Track.TransportToPoint(TrackHit_Aprox);

  for(size_t i = 0; i < Hits_Si.size(); ++i)
    {
      if(ifCut_MinEnergyDeposition_SiHit && (Hits_Si[i][0] < MinEnergyDeposition_SiHit))
        continue;

      double temp_dist = std::sqrt(std::pow(Hits_Si[i][1] - Track.GetX(), 2.) + std::pow(Hits_Si[i][2] - Track.GetY(), 2.));
      if(ifCut_MaxDist_SiHit && (temp_dist > MaxDist_SiHit))
        continue;

      Track_Sihit.emplace_back(Hits_Si[i]);
    }

  if(Track_Sihit.size() == 0)
    Track_Sihit.emplace_back(Empty_Sihit);

  return;
}

template <class Out>
void TDecayVertex<Out>::MotherDaughtersTrack_SiHits(std::vector<KFParticle>& FragmentTracks, std::vector<KFParticle>& PionTracks,
                    std::vector<KFParticle>& MotherTracks, std::vector<std::tuple<size_t, size_t>>& RefDaughtersTracks,
                    std::vector<std::vector<double> >& Hits_Si1, std::vector<std::vector<double> >& Hits_Si2,
                    std::vector<std::tuple<int, std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>>>>& Fragment_SiHits,
                    std::vector<std::tuple<int, std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>>>>& Pion_SiHits,
                    std::vector<std::tuple<int, std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>>>>& Mother_SiHits)
{
  std::vector<std::vector<double> > Empty_Sihit = {{-1., -1., -1., -1.}};
  std::vector<std::vector<double> > temp_Track_Sihit1;
  std::vector<std::vector<double> > temp_Track_Sihit2;

  for(size_t i = 0; i < MotherTracks.size(); ++i)
    {
      MotherTracks[i].TransportToDecayVertex();

      size_t temp_id_fragment = get<0>(RefDaughtersTracks[i]);
      size_t temp_id_pion = get<1>(RefDaughtersTracks[i]);

      if((MotherTracks[i].GetZ() > Zf_target + Dist_to_Target) && (MotherTracks[i].GetZ() < Z_plane_Si1 - Dist_to_Silicons))
        {
          SiHitsFinder(FragmentTracks[temp_id_fragment], Hits_Si1, temp_Track_Sihit1);
          SiHitsFinder(FragmentTracks[temp_id_fragment], Hits_Si2, temp_Track_Sihit2);
          if((temp_Track_Sihit1[0][0] != Empty_Sihit[0][0]) || (temp_Track_Sihit2[0][0] != Empty_Sihit[0][0]))
            {
              auto temp_fragmenthit = std::make_tuple(temp_Track_Sihit1, temp_Track_Sihit2);
              Fragment_SiHits.emplace_back(std::make_tuple(i, temp_fragmenthit));
            }

          SiHitsFinder(PionTracks[temp_id_pion], Hits_Si1, temp_Track_Sihit1);
          SiHitsFinder(PionTracks[temp_id_pion], Hits_Si2, temp_Track_Sihit2);
          if((temp_Track_Sihit1[0][0] != Empty_Sihit[0][0]) || (temp_Track_Sihit2[0][0] != Empty_Sihit[0][0]))
            {
              auto temp_pionhit = std::make_tuple(temp_Track_Sihit1, temp_Track_Sihit2);
              Pion_SiHits.emplace_back(std::make_tuple(i, temp_pionhit));
            }
        }

      else if((MotherTracks[i].GetZ() > Z_plane_Si1 + Dist_to_Silicons) && (MotherTracks[i].GetZ() < Z_plane_Si2 - Dist_to_Silicons))
        {
          SiHitsFinder(MotherTracks[i], Hits_Si1, temp_Track_Sihit1);
          if(temp_Track_Sihit1[0][0] != Empty_Sihit[0][0])
            {
              auto temp_motherhit = std::make_tuple(temp_Track_Sihit1, Empty_Sihit);
              Mother_SiHits.emplace_back(std::make_tuple(i, temp_motherhit));
            }


          SiHitsFinder(FragmentTracks[temp_id_fragment], Hits_Si2, temp_Track_Sihit2);
          if(temp_Track_Sihit2[0][0] != Empty_Sihit[0][0])
            {
              auto temp_fragmenthit = std::make_tuple(Empty_Sihit, temp_Track_Sihit2);
              Fragment_SiHits.emplace_back(std::make_tuple(i, temp_fragmenthit));
            }

          SiHitsFinder(PionTracks[temp_id_pion], Hits_Si2, temp_Track_Sihit2);
          if(temp_Track_Sihit2[0][0] != Empty_Sihit[0][0])
            {
              auto temp_pionhit = std::make_tuple(Empty_Sihit, temp_Track_Sihit2);
              Pion_SiHits.emplace_back(std::make_tuple(i, temp_pionhit));
            }
        }

      else if((MotherTracks[i].GetZ() > Z_plane_Si2 + Dist_to_Silicons) && (MotherTracks[i].GetZ() < Zo_minifibers - Dist_to_Minifibers))
        {
          SiHitsFinder(MotherTracks[i], Hits_Si1, temp_Track_Sihit1);
          SiHitsFinder(MotherTracks[i], Hits_Si2, temp_Track_Sihit2);
          if((temp_Track_Sihit1[0][0] != Empty_Sihit[0][0]) || (temp_Track_Sihit2[0][0] != Empty_Sihit[0][0]))
            {
              auto temp_motherhit = std::make_tuple(temp_Track_Sihit1, temp_Track_Sihit2);
              Mother_SiHits.emplace_back(std::make_tuple(i, temp_motherhit));
            }
        }
    }

  return;
}


template <class Out>
void TDecayVertex<Out>::SiHitsFinder2(KFParticle& Track, int idSilicon, int stripDirection, //Strip direction: 1 -> X, 2 -> Y
                                  std::vector<std::tuple<double, size_t> >& Hits_Si, std::vector<std::vector<double> >& Track_Sihit)
{
  if(Hits_Si.size() == 0)
    return;

  std::vector<double> Empty_Sihit = {-1., -1., -1., -1.}; // (E, XorY, sigmaXorY, Z)
  Track_Sihit.clear();

  if(idSilicon == 1)
    {
      if(stripDirection == 1)
        Z_plane = Z_plane_Si1x;
      else if(stripDirection == 2)
        Z_plane = Z_plane_Si1y;

      widthStrip = widthStrip_Si1;
      lenghtSi   = lenghtSi_Si1;

      sigma = sigma_Si1;
    }

  else if(idSilicon == 2)
    {
      if(stripDirection == 1)
        Z_plane = Z_plane_Si2x;
      else if(stripDirection == 2)
        Z_plane = Z_plane_Si2y;

      widthStrip = widthStrip_Si2;
      lenghtSi   = lenghtSi_Si2;
      sigma = sigma_Si2;
    }

  nStrips = static_cast<int>(lenghtSi / widthStrip);

  double track_parameter = (Z_plane - Track.GetZ()) / Track.GetPz();
  float X_track = Track.GetX() + track_parameter * Track.GetPx();
  float Y_track = Track.GetY() + track_parameter * Track.GetPy();

  float TrackHit_Aprox[3] = {X_track, Y_track, Z_plane};
  Track.TransportToPoint(TrackHit_Aprox);

  for(size_t i = 0; i < Hits_Si.size(); ++i)
    {
      if(ifCut_MinEnergyDeposition_SiHit && (get<0>(Hits_Si[i]) < MinEnergyDeposition_SiHit))
        continue;

      double posHit = -50.;
      double posTrack = -100.;

      if(stripDirection == 1)
        {
          posHit = -lenghtSi / 2. + (get<1>(Hits_Si[i])%nStrips + 0.5) * widthStrip;
          posTrack = Track.GetX();
        }
      else if(stripDirection == 2)
        {
          posHit = +lenghtSi / 2. - (get<1>(Hits_Si[i])%nStrips + 0.5) * widthStrip;
          posTrack = Track.GetY();
        }

      double temp_dist = std::abs(posHit - posTrack);
      if(ifCut_MaxDist_SiHit && (temp_dist > MaxDist_SiHit))
        continue;

      std::vector<double> temp_Sihit = {get<0>(Hits_Si[i]), posHit, sigma, Z_plane}; // (E, XorY, sigmaXorY, Z)
      Track_Sihit.emplace_back(temp_Sihit);
    }

  if(Track_Sihit.size() == 0)
    Track_Sihit.emplace_back(Empty_Sihit);

  return;
}

template <class Out>
void TDecayVertex<Out>::MotherDaughtersTrack_SiHits2(std::vector<KFParticle>& FragmentTracks, std::vector<KFParticle>& PionTracks,
                    std::vector<KFParticle>& MotherTracks, std::vector<std::tuple<size_t, size_t>>& RefDaughtersTracks,
                    std::vector<std::tuple<double, size_t> >& HitsX_Si1, std::vector<std::tuple<double, size_t> >& HitsY_Si1,
                    std::vector<std::tuple<double, size_t> >& HitsX_Si2, std::vector<std::tuple<double, size_t> >& HitsY_Si2,
                    std::vector<std::tuple<int, std::vector<std::vector<std::vector<double>>>>>& Fragment_SiHits,
                    std::vector<std::tuple<int, std::vector<std::vector<std::vector<double>>>>>& Pion_SiHits,
                    std::vector<std::tuple<int, std::vector<std::vector<std::vector<double>>>>>& Mother_SiHits)
{
  std::vector<std::vector<double> > Empty_Sihit = {{-1., -1., -1., -1.}};

  std::vector<std::vector<double> > temp_TrackHits_Si1X;
  std::vector<std::vector<double> > temp_TrackHits_Si1Y;
  std::vector<std::vector<double> > temp_TrackHits_Si2X;
  std::vector<std::vector<double> > temp_TrackHits_Si2Y;

  for(size_t i = 0; i < MotherTracks.size(); ++i)
    {
      MotherTracks[i].TransportToDecayVertex();

      std::cout << "HELLO\n";
      size_t temp_id_fragment = get<0>(RefDaughtersTracks[i]);
      size_t temp_id_pion = get<1>(RefDaughtersTracks[i]);

      if((MotherTracks[i].GetZ() > Zf_target + Dist_to_Target) && (MotherTracks[i].GetZ() < Z_plane_Si1 - Dist_to_Silicons))
        {
                std::cout << "HELLO1\n";

          SiHitsFinder2(FragmentTracks[temp_id_fragment], 1, 1, HitsX_Si1, temp_TrackHits_Si1X);
          SiHitsFinder2(FragmentTracks[temp_id_fragment], 1, 2, HitsY_Si1, temp_TrackHits_Si1Y);
          SiHitsFinder2(FragmentTracks[temp_id_fragment], 2, 1, HitsX_Si2, temp_TrackHits_Si2X);
          SiHitsFinder2(FragmentTracks[temp_id_fragment], 2, 2, HitsY_Si2, temp_TrackHits_Si2Y);
          if((temp_TrackHits_Si1X[0][0] != Empty_Sihit[0][0]) || (temp_TrackHits_Si1Y[0][0] != Empty_Sihit[0][0])
                || (temp_TrackHits_Si2X[0][0] != Empty_Sihit[0][0]) || (temp_TrackHits_Si2Y[0][0] != Empty_Sihit[0][0]))
            {
              std::vector<std::vector<std::vector<double>>> temp_fragmenthit =
                                                    {temp_TrackHits_Si1X, temp_TrackHits_Si1Y, temp_TrackHits_Si2X, temp_TrackHits_Si2Y};
              Fragment_SiHits.emplace_back(std::make_tuple(i, temp_fragmenthit));
            }

          SiHitsFinder2(PionTracks[temp_id_pion], 1, 1, HitsX_Si1, temp_TrackHits_Si1X);
          SiHitsFinder2(PionTracks[temp_id_pion], 1, 2, HitsY_Si1, temp_TrackHits_Si1Y);
          SiHitsFinder2(PionTracks[temp_id_pion], 2, 1, HitsX_Si2, temp_TrackHits_Si2X);
          SiHitsFinder2(PionTracks[temp_id_pion], 2, 2, HitsY_Si2, temp_TrackHits_Si2Y);
          if((temp_TrackHits_Si1X[0][0] != Empty_Sihit[0][0]) || (temp_TrackHits_Si1Y[0][0] != Empty_Sihit[0][0])
                || (temp_TrackHits_Si2X[0][0] != Empty_Sihit[0][0]) || (temp_TrackHits_Si2Y[0][0] != Empty_Sihit[0][0]))
            {
              std::vector<std::vector<std::vector<double>>> temp_pionhit =
                                                    {temp_TrackHits_Si1X, temp_TrackHits_Si1Y, temp_TrackHits_Si2X, temp_TrackHits_Si2Y};
              Pion_SiHits.emplace_back(std::make_tuple(i, temp_pionhit));
            }
        }

      else if((MotherTracks[i].GetZ() > Z_plane_Si1 + Dist_to_Silicons) && (MotherTracks[i].GetZ() < Z_plane_Si2 - Dist_to_Silicons))
        {
                std::cout << "HELLO2\n";

          SiHitsFinder2(MotherTracks[i], 1, 1, HitsX_Si1, temp_TrackHits_Si1X);
          SiHitsFinder2(MotherTracks[i], 1, 2, HitsY_Si1, temp_TrackHits_Si1Y);
          if((temp_TrackHits_Si1X[0][0] != Empty_Sihit[0][0]) || (temp_TrackHits_Si1Y[0][0] != Empty_Sihit[0][0]))
            {
              std::vector<std::vector<std::vector<double>>> temp_motherhit =
                                                    {temp_TrackHits_Si1X, temp_TrackHits_Si1Y, Empty_Sihit, Empty_Sihit};
              Mother_SiHits.emplace_back(std::make_tuple(i, temp_motherhit));
            }

          SiHitsFinder2(FragmentTracks[temp_id_fragment], 2, 1, HitsX_Si2, temp_TrackHits_Si2X);
          SiHitsFinder2(FragmentTracks[temp_id_fragment], 2, 2, HitsY_Si2, temp_TrackHits_Si2Y);
          if((temp_TrackHits_Si2X[0][0] != Empty_Sihit[0][0]) || (temp_TrackHits_Si2Y[0][0] != Empty_Sihit[0][0]))
            {
              std::vector<std::vector<std::vector<double>>> temp_fragmenthit =
                                                    {Empty_Sihit, Empty_Sihit, temp_TrackHits_Si2X, temp_TrackHits_Si2Y};
              Fragment_SiHits.emplace_back(std::make_tuple(i, temp_fragmenthit));
            }

          SiHitsFinder2(PionTracks[temp_id_pion], 2, 1, HitsX_Si2, temp_TrackHits_Si2X);
          SiHitsFinder2(PionTracks[temp_id_pion], 2, 2, HitsY_Si2, temp_TrackHits_Si2Y);
          if((temp_TrackHits_Si2X[0][0] != Empty_Sihit[0][0]) || (temp_TrackHits_Si2Y[0][0] != Empty_Sihit[0][0]))
            {
              std::vector<std::vector<std::vector<double>>> temp_pionhit =
                                                    {Empty_Sihit, Empty_Sihit, temp_TrackHits_Si2X, temp_TrackHits_Si2Y};
              Pion_SiHits.emplace_back(std::make_tuple(i, temp_pionhit));
            }
        }

      else if((MotherTracks[i].GetZ() > Z_plane_Si2 + Dist_to_Silicons) && (MotherTracks[i].GetZ() < Zo_minifibers - Dist_to_Minifibers))
        {
                std::cout << "HELLO3\n";

          SiHitsFinder2(MotherTracks[i], 1, 1, HitsX_Si1, temp_TrackHits_Si1X);
          SiHitsFinder2(MotherTracks[i], 1, 2, HitsY_Si1, temp_TrackHits_Si1Y);
          SiHitsFinder2(MotherTracks[i], 2, 1, HitsX_Si2, temp_TrackHits_Si2X);
          SiHitsFinder2(MotherTracks[i], 2, 2, HitsY_Si2, temp_TrackHits_Si2Y);

                std::cout << "HELLO4\n";

          if((temp_TrackHits_Si1X[0][0] != Empty_Sihit[0][0]) || (temp_TrackHits_Si1Y[0][0] != Empty_Sihit[0][0])
                || (temp_TrackHits_Si2X[0][0] != Empty_Sihit[0][0]) || (temp_TrackHits_Si2Y[0][0] != Empty_Sihit[0][0]))
            {
                    std::cout << "HELLO5\n";

              std::vector<std::vector<std::vector<double>>> temp_motherhit =
                                                    {temp_TrackHits_Si1X, temp_TrackHits_Si1Y, temp_TrackHits_Si2X, temp_TrackHits_Si2Y};
              Mother_SiHits.emplace_back(std::make_tuple(i, temp_motherhit));
            }
        }

              std::cout << "HELLO6\n";

    }

          std::cout << "HELLO9\n";


  return;
}
*/
/*
  TVector3 Mother_Hit_Si1;
  MotherTrack_SiHit(PrimVtxRecons, DecayVtxRecons, Z_plane_Si1, Hits_Si1, Mother_Hit_Si1);
  if(Mother_Hit_Si1 == Vect_Zero)
    return;

  TVector3 Mother_Hit_Si2;
  MotherTrack_SiHit(PrimVtxRecons, DecayVtxRecons, Z_plane_Si2, Hits_Si2, Mother_Hit_Si2);
  if(Mother_Hit_Si2 == Vect_Zero)
    return;

  TVector3 Si_MotherMom = Mother_Hit_Si2 - Mother_Hit_Si1;
  double sigma2         = std::pow(widthStrip_Si1, 2.) / 12.;

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
  
  for(int iF = 0; iF < 10; ++iF)
    Si_MotherTrack.SetFieldCoeff(fieldMDCParameters[iF], iF);
*/


template <class Out>
void TDecayVertex<Out>::AllTrackstoDecayVertex_Vfunction(std::vector<KFParticle>& AllTracks, TVector3& Old_DecayVertexRecons, TVector3& DecayVertexRecons)
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
              temp_f[j] = f_function(AllTracks[j], temp_PosXYZ);

          Vnew = V_function(temp_f);

          if(Vnew > V)
            {
              V = Vnew;
              DecayVertexRecons.SetXYZ(temp_PosXYZ.X(), temp_PosXYZ.Y(), temp_PosXYZ.Z());
            }
        }
    }
}

template <class Out>
void TDecayVertex<Out>::AllTrackstoDecayVertex_Centroids(std::vector<KFParticle>& AllTracks, TVector3& DecayVertexRecons)
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


template class TDecayVertex<MCAnaEventG4Sol>;
template class TDecayVertex<Ana_WasaEvent>;
