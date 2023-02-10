#include "TPrimaryVertex.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"
#include "FullRecoEvent.hh"
#include "ReturnRes.hh"

#include <tuple>

#include "TVector3.h"

//#define DEBUG_PRIMVTX

//#define RECONS_HITS_MULTIPLICITY
//#define HIT_RECONS_CHECK
//#define TRACK_RECONS_CHECK
//#define MOTHER_DAUGHTERS_CHECK
#define VERTEX_RECONS_CHECK
#define COVARIANCE_MATRIX
//#define DECAY_VERTEX

using namespace std;
using namespace G4Sol;

template <class Out>
TPrimaryVertex<Out>::TPrimaryVertex(const THyphiAttributes& attribut)
    : TDataProcessInterface<Out>("PrimaryVertexReco"), att(attribut), SiliconHitsSD_Si1(1), SiliconHitsSD_Si2(2),
        SiliconHitsSD_Si3(3), SiliconHitsSD_Si4(4)
{
  TranslationZ_Target_System(att.Target_PositionZ);
  rand = new TRandom3();
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
int TPrimaryVertex<Out>::Exec(FullRecoEvent& RecoEvent, Out* ) { return FinderPrimaryVertex(RecoEvent); }

template <class Out>
ReturnRes::InfoM TPrimaryVertex<Out>::SoftExit(int result_full) {
   
  if(result_full == -1)
    {
      att._logger->debug("No enough candidate tracks for primary vertex recons");
      LocalHisto.h_PrimVtxstats->Fill("CandidateTracks=0", 1.);
      //return ReturnRes::PrimVtxError; //Check in the future
      return ReturnRes::Fine;
    }

  else if(result_full == -2)
    {
      att._logger->debug("No enough candidate tracks for decay vertex recons");
      LocalHisto.h_PrimVtxstats->Fill("CandidateDecayTracks=0", 1.);
      return ReturnRes::Fine;
    }

  else if(result_full == -3)
    {
      att._logger->debug("No simulated hypernucleus");
      LocalHisto.h_PrimVtxstats->Fill("Mother=0", 1.);
      //return ReturnRes::PrimVtxError;
      return ReturnRes::Fine;
    }

  else if(result_full == -4)
    {
      att._logger->debug("No enough hits in Silicons");
      LocalHisto.h_PrimVtxstats->Fill("SiliconHits<2", 1.);
      //return ReturnRes::PrimVtxError;
      return ReturnRes::Fine;
    }

  LocalHisto.h_PrimVtxstats->Fill("Fine", 1.);

  return ReturnRes::Fine; 
}





template <class Out>
void TPrimaryVertex<Out>::SelectHists()
{
  LocalHisto.h_HitMultiplicity_Si1          = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HitMultiplicity_Si1);
  LocalHisto.h_HitMultiplicityRecons_Si1    = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HitMultiplicityRecons_Si1);
  LocalHisto.h_HitMultiplicityDiff_Si1      = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HitMultiplicityDiff_Si1);
  LocalHisto.h_HitMultiplicityDiffNHits_Si1 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HitMultiplicityDiffNHits_Si1);

  LocalHisto.h_EnergyDiffStrips_Si1          = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EnergyDiffStrips_Si1);
  LocalHisto.h_nEventsGoodrecons_Si1         = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_nEventsGoodrecons_Si1);
  LocalHisto.h_nEventsGhost_Si1              = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_nEventsGhost_Si1);
  LocalHisto.h_nEventsGoodreconsGhost_Si1    = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_nEventsGoodreconsGhost_Si1);
  LocalHisto.h_nEventsRealGoodrecons_Si1     = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_nEventsRealGoodrecons_Si1);
  LocalHisto.h_nEventsRealRejectPad_Si1      = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_nEventsRealRejectPad_Si1);

  LocalHisto.h_HitMultiplicity_Si2          = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HitMultiplicity_Si2);
  LocalHisto.h_HitMultiplicityRecons_Si2    = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HitMultiplicityRecons_Si2);
  LocalHisto.h_HitMultiplicityDiff_Si2      = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HitMultiplicityDiff_Si2);
  LocalHisto.h_HitMultiplicityDiffNHits_Si2 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_HitMultiplicityDiffNHits_Si2);

  LocalHisto.h_EnergyDiffStrips_Si2          = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EnergyDiffStrips_Si2);
  LocalHisto.h_nEventsGoodrecons_Si2         = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_nEventsGoodrecons_Si2);
  LocalHisto.h_nEventsGhost_Si2              = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_nEventsGhost_Si2);
  LocalHisto.h_nEventsGoodreconsGhost_Si2    = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_nEventsGoodreconsGhost_Si2);
  LocalHisto.h_nEventsRealGoodrecons_Si2     = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_nEventsRealGoodrecons_Si2);
  LocalHisto.h_nEventsRealRejectPad_Si2      = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_nEventsRealRejectPad_Si2);

  LocalHisto.h_MFCheck_Theta_MomSi1MomSi2 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_MFCheck_Theta_MomSi1MomSi2);
  LocalHisto.h_MFCheck_Dist_MomSi1HitSi2  = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_MFCheck_Dist_MomSi1HitSi2);
  LocalHisto.h_MFCheck_Dist_MomSi2HitSi1  = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_MFCheck_Dist_MomSi2HitSi1);
  LocalHisto.h_MFCheck_Dist_MomSi1HitIP   = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_MFCheck_Dist_MomSi1HitIP);
  LocalHisto.h_MFCheck_Dist_MomSi2HitIP   = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_MFCheck_Dist_MomSi2HitIP);

  LocalHisto.h_EnergyStripEnergyTotalReal = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EnergyStripEnergyTotalReal);
  LocalHisto.h_EnergyStripEnergyTotal     = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EnergyStripEnergyTotal);
  LocalHisto.h_EnergyDiffSilicons         = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EnergyDiffSilicons);

  LocalHisto.h_EnergyDepositionMother    = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EnergyDepositionMother);
  LocalHisto.h_EnergyDepositionDaughters = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EnergyDepositionDaughters);

  LocalHisto.h_nTrackCandidates   = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_nTrackCandidates);
  LocalHisto.h_DistanceBeamTracks = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DistanceBeamTracks);
  LocalHisto.h_PosZBeamTracks     = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_PosZBeamTracks);
  LocalHisto.h_thetaTracks        = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_thetaTracks);
  LocalHisto.h_thetaResol         = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_thetaResol);

  LocalHisto.h_Acc_ThetaCandidates = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Acc_ThetaCandidates);
  LocalHisto.h_Acc_ThetaAllReal    = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Acc_ThetaAllReal);

  LocalHisto.h_nCandidatesRealTracks          = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_nCandidatesRealTracks);
  LocalHisto.h_nCandidatesRealTracks_IfRecons = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_nCandidatesRealTracks_IfRecons);

  LocalHisto.h_nHypernucleiTrack = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_nHypernucleiTrack);
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

  LocalHisto.h_IP_DecayDistance  = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_IP_DecayDistance);
  LocalHisto.h_IP_DecayDistanceX = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_IP_DecayDistanceX);
  LocalHisto.h_IP_DecayDistanceY = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_IP_DecayDistanceY);
  LocalHisto.h_IP_DecayDistanceZ = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_IP_DecayDistanceZ);

  LocalHisto.h_DecayPositionDistance  = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayPositionDistance);
  LocalHisto.h_DecayPositionDistanceX = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayPositionDistanceX);
  LocalHisto.h_DecayPositionDistanceY = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayPositionDistanceY);
  LocalHisto.h_DecayPositionDistanceZ = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_DecayPositionDistanceZ);

  LocalHisto.h_PrimVtxstats = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_PrimVtxstats);
  LocalHisto.h_PrimStatus = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_PrimStatus);
}

template <class Out>
int TPrimaryVertex<Out>::FinderPrimaryVertex(FullRecoEvent& RecoEvent)
{
  LocalHisto.h_PrimStatus->Fill("Si1x", RecoEvent.Si_HitsEnergyLayer[5].size(), 1);
  LocalHisto.h_PrimStatus->Fill("Si1y", RecoEvent.Si_HitsEnergyLayer[4].size(), 1);
  LocalHisto.h_PrimStatus->Fill("Si2x", RecoEvent.Si_HitsEnergyLayer[7].size(), 1);
  LocalHisto.h_PrimStatus->Fill("Si2y", RecoEvent.Si_HitsEnergyLayer[6].size(), 1);
  LocalHisto.h_PrimStatus->Fill(
      "Si1xy", (RecoEvent.Si_HitsEnergyLayer[5].size() < 2 && RecoEvent.Si_HitsEnergyLayer[4].size() < 2), 1);
  LocalHisto.h_PrimStatus->Fill(
      "Si2xy", (RecoEvent.Si_HitsEnergyLayer[7].size() < 2 && RecoEvent.Si_HitsEnergyLayer[6].size() < 2), 1);
  LocalHisto.h_PrimStatus->Fill("Si_all", 1., 1);

  if(Si_4layers)
    if(((RecoEvent.Si_HitsEnergyLayer[0].size() < 2 && RecoEvent.Si_HitsEnergyLayer[1].size() < 2) ||
      (RecoEvent.Si_HitsEnergyLayer[2].size() < 2 && RecoEvent.Si_HitsEnergyLayer[3].size() < 2)) &&
      ((RecoEvent.Si_HitsEnergyLayer[5].size() < 2 && RecoEvent.Si_HitsEnergyLayer[4].size() < 2) ||
      (RecoEvent.Si_HitsEnergyLayer[7].size() < 2 && RecoEvent.Si_HitsEnergyLayer[6].size() < 2)))
        return -4;

  if((RecoEvent.Si_HitsEnergyLayer[5].size() < 2 && RecoEvent.Si_HitsEnergyLayer[4].size() < 2) ||
     (RecoEvent.Si_HitsEnergyLayer[7].size() < 2 && RecoEvent.Si_HitsEnergyLayer[6].size() < 2))
        return -4;

  LocalHisto.h_PrimStatus->Fill("Si_acc", 1., 1);

  // Check all real particles
  size_t nRealTracksCounter = 0;

  for(auto& it: RecoEvent.TrackDAFInit)
    {
      if(it.second.charge == 0)
        continue;

      if(it.second.posZ >= Z_plane_Si1)
        continue;

      nRealTracksCounter += 1;

      double temp_thetaAllReal = std::atan(std::sqrt(std::pow(it.second.momX, 2.) + std::pow(it.second.momY, 2.))
                                            / it.second.momZ) * 180. / M_PI;
      LocalHisto.h_Acc_ThetaAllReal->Fill(temp_thetaAllReal, 1.);
    }

  // Reconstruct the real hits from the simulation
  std::vector<std::tuple<double, double, double, size_t, double, double, std::string> > HitEnergyPosXYreal_Si1{};
  std::vector<std::tuple<size_t,TVector3,TVector3>> HitIdMomPos_Si1{};
  simulHitstoRealHits(RecoEvent, HitEnergyPosXYreal_Si1, G4Sol::Si1y_SD_pad, G4Sol::Si1x_SD_pad, LocalHisto.h_EnergyDiffStrips_Si1, HitIdMomPos_Si1);

  std::vector<std::tuple<double, double, double, size_t, double, double, std::string> > HitEnergyPosXYreal_Si2{};
  std::vector<std::tuple<size_t,TVector3,TVector3>> HitIdMomPos_Si2{};
  simulHitstoRealHits(RecoEvent, HitEnergyPosXYreal_Si2, G4Sol::Si2y_SD_pad, G4Sol::Si2x_SD_pad, LocalHisto.h_EnergyDiffStrips_Si1, HitIdMomPos_Si2);


  std::vector<std::tuple<double, double, double, size_t, double, double, std::string> > HitEnergyPosXYreal_Si3{};
  std::vector<std::tuple<size_t,TVector3,TVector3>> HitIdMomPos_Si3{};

  std::vector<std::tuple<double, double, double, size_t, double, double, std::string> > HitEnergyPosXYreal_Si4{};
  std::vector<std::tuple<size_t,TVector3,TVector3>> HitIdMomPos_Si4{};
  
  if(Si_4layers)
  {
    simulHitstoRealHits(RecoEvent, HitEnergyPosXYreal_Si3, G4Sol::Si1x, G4Sol::Si1y, LocalHisto.h_EnergyDiffStrips_Si2, HitIdMomPos_Si3);
    simulHitstoRealHits(RecoEvent, HitEnergyPosXYreal_Si4, G4Sol::Si2x, G4Sol::Si2y, LocalHisto.h_EnergyDiffStrips_Si2, HitIdMomPos_Si4);
  }

  //Check magnetic field effects
  MFcheck(HitIdMomPos_Si1, HitIdMomPos_Si2, RecoEvent.InteractionPoint);
  
  if(Si_4layers)
    MFcheck(HitIdMomPos_Si3, HitIdMomPos_Si4, RecoEvent.InteractionPoint);


  std::vector<std::tuple<double, size_t> > HitEnergyLayerX_Si1{};
  for(auto it_HitSi : RecoEvent.Si_HitsEnergyLayer[5])
    HitEnergyLayerX_Si1.emplace_back(std::make_tuple(it_HitSi.second, it_HitSi.first));

  std::vector<std::tuple<double, size_t> > HitEnergyLayerY_Si1{};
  for(auto it_HitSi : RecoEvent.Si_HitsEnergyLayer[4])
    HitEnergyLayerY_Si1.emplace_back(std::make_tuple(it_HitSi.second, it_HitSi.first));

  std::vector<std::tuple<double, size_t> > HitEnergyLayerX_Si2{};
  for(auto it_HitSi : RecoEvent.Si_HitsEnergyLayer[7])
    HitEnergyLayerX_Si2.emplace_back(std::make_tuple(it_HitSi.second, it_HitSi.first));

  std::vector<std::tuple<double, size_t> > HitEnergyLayerY_Si2{};
  for(auto it_HitSi : RecoEvent.Si_HitsEnergyLayer[6])
    HitEnergyLayerY_Si2.emplace_back(std::make_tuple(it_HitSi.second, it_HitSi.first));


  std::vector<std::tuple<double, size_t> > HitEnergyLayerX_Si3{};
  std::vector<std::tuple<double, size_t> > HitEnergyLayerY_Si3{};
  std::vector<std::tuple<double, size_t> > HitEnergyLayerX_Si4{};
  std::vector<std::tuple<double, size_t> > HitEnergyLayerY_Si4{};

  if(Si_4layers)
    {
      for(auto it_HitSi : RecoEvent.Si_HitsEnergyLayer[0])
        HitEnergyLayerX_Si3.emplace_back(std::make_tuple(it_HitSi.second, it_HitSi.first));

      for(auto it_HitSi : RecoEvent.Si_HitsEnergyLayer[1])
        HitEnergyLayerY_Si3.emplace_back(std::make_tuple(it_HitSi.second, it_HitSi.first));

      for(auto it_HitSi : RecoEvent.Si_HitsEnergyLayer[2])
        HitEnergyLayerX_Si4.emplace_back(std::make_tuple(it_HitSi.second, it_HitSi.first));

      for(auto it_HitSi : RecoEvent.Si_HitsEnergyLayer[3])
        HitEnergyLayerY_Si4.emplace_back(std::make_tuple(it_HitSi.second, it_HitSi.first));
    }


#ifdef DEBUG_PRIMVTX
  std::cout << "Size HitEnergyLayerX_Si1: " << HitEnergyLayerX_Si1.size() << "\n";
  std::cout << "Size HitEnergyLayerY_Si1: " << HitEnergyLayerY_Si1.size() << "\n";
  std::cout << "Size HitEnergyLayerX_Si2: " << HitEnergyLayerX_Si2.size() << "\n";
  std::cout << "Size HitEnergyLayerY_Si2: " << HitEnergyLayerY_Si2.size() << "\n\n";
  
  if(Si_4layers)
    {
      std::cout << "Size HitEnergyLayerX_Si3: " << HitEnergyLayerX_Si3.size() << "\n";
      std::cout << "Size HitEnergyLayerY_Si3: " << HitEnergyLayerY_Si3.size() << "\n";
      std::cout << "Size HitEnergyLayerX_Si4: " << HitEnergyLayerX_Si4.size() << "\n";
      std::cout << "Size HitEnergyLayerY_Si4: " << HitEnergyLayerY_Si4.size() << "\n\n";
    }
#endif

  //Combine strips if set in the header
  SiliconHitsSD_Si1.CombineStrips(HitEnergyLayerX_Si1);
  SiliconHitsSD_Si1.CombineStrips(HitEnergyLayerY_Si1);
  SiliconHitsSD_Si2.CombineStrips(HitEnergyLayerX_Si2);
  SiliconHitsSD_Si2.CombineStrips(HitEnergyLayerY_Si2);

  if(Si_4layers)
    {
      SiliconHitsSD_Si3.CombineStrips(HitEnergyLayerX_Si3);
      SiliconHitsSD_Si3.CombineStrips(HitEnergyLayerY_Si3);
      SiliconHitsSD_Si4.CombineStrips(HitEnergyLayerX_Si4);
      SiliconHitsSD_Si4.CombineStrips(HitEnergyLayerY_Si4);
    }

  //Copy silicon signals for later checking the produced by daughters/mother
  RecoEvent.HitsX_Si1 = HitEnergyLayerX_Si1;
  RecoEvent.HitsY_Si1 = HitEnergyLayerY_Si1;
  RecoEvent.HitsX_Si2 = HitEnergyLayerX_Si2;
  RecoEvent.HitsY_Si2 = HitEnergyLayerY_Si2;

/* FOR PREVIOUS SILICON GEOMETRY
  // Obtain the hits from the energy signals
  std::vector<std::vector<double> > HitEnergyPosXY_Si1{};
  SiliconHits_Si1.SignalstoHits(HitEnergyLayerX_Si1, HitEnergyLayerY_Si1, HitEnergyPosXY_Si1);
  RecoEvent.Hits_Si1 = HitEnergyPosXY_Si1;

  std::vector<std::vector<double> > HitEnergyPosXY_Si2{};
  SiliconHits_Si2.SignalstoHits(HitEnergyLayerX_Si2, HitEnergyLayerY_Si2, HitEnergyPosXY_Si2);
  RecoEvent.Hits_Si2 = HitEnergyPosXY_Si2;

  FOR HAMATATSU SILICON GEOMETRY
  // Obtain the hits from the energy signals
  std::vector<std::vector<double> > HitEnergyPosXY_Si1{};
  SiliconHitsSD_Si1.SignalstoHits_SD(HitEnergyLayerX_Si1, HitEnergyLayerY_Si1, HitEnergyPosXY_Si1);
  RecoEvent.Hits_Si1 = HitEnergyPosXY_Si1;

  std::vector<std::vector<double> > HitEnergyPosXY_Si2{};
  SiliconHitsSD_Si2.SignalstoHits_SD(HitEnergyLayerX_Si2, HitEnergyLayerY_Si2, HitEnergyPosXY_Si2);
  RecoEvent.Hits_Si2 = HitEnergyPosXY_Si2;
*/

  //FOR NEW PAD SILICON GEOMETRY
  // Obtain the hits from the energy signals
  std::vector<std::vector<double> > HitEnergyPosXY_Si1{};
  SiliconHitsSD_Si1.SignalstoHits_SDpad(HitEnergyLayerX_Si1, HitEnergyLayerY_Si1, HitEnergyPosXY_Si1);
  RecoEvent.Hits_Si1 = HitEnergyPosXY_Si1;

  std::vector<std::vector<double> > HitEnergyPosXY_Si2{};
  SiliconHitsSD_Si2.SignalstoHits_SDpad(HitEnergyLayerX_Si2, HitEnergyLayerY_Si2, HitEnergyPosXY_Si2);
  RecoEvent.Hits_Si2 = HitEnergyPosXY_Si2;

  std::vector<std::vector<double> > HitEnergyPosXY_Si3{};
  std::vector<std::vector<double> > HitEnergyPosXY_Si4{};

  if(Si_4layers)
    {
      SiliconHitsSD_Si3.SignalstoHits_SDpad(HitEnergyLayerX_Si3, HitEnergyLayerY_Si3, HitEnergyPosXY_Si3);
      RecoEvent.Hits_Si3 = HitEnergyPosXY_Si3;

      SiliconHitsSD_Si4.SignalstoHits_SDpad(HitEnergyLayerX_Si4, HitEnergyLayerY_Si4, HitEnergyPosXY_Si4);
      RecoEvent.Hits_Si4 = HitEnergyPosXY_Si4;
    }


#ifdef HIT_RECONS_CHECK

  size_t nGoodrecons_Si1 = 0;
  nGoodEventsCounter(HitEnergyPosXY_Si1, HitEnergyPosXYreal_Si1, widthStrip_Si1, nGoodrecons_Si1);
  int nGhost_Si1    = HitEnergyPosXY_Si1.size() - nGoodrecons_Si1;

  LocalHisto.h_HitMultiplicity_Si1->Fill(HitEnergyPosXYreal_Si1.size(), 1.);
  LocalHisto.h_HitMultiplicityRecons_Si1->Fill(HitEnergyPosXY_Si1.size(), 1.);
  LocalHisto.h_HitMultiplicityDiff_Si1->Fill(
      static_cast<int>(HitEnergyPosXY_Si1.size()) - static_cast<int>(HitEnergyPosXYreal_Si1.size()), 1.);
  LocalHisto.h_HitMultiplicityDiffNHits_Si1->Fill(static_cast<int>(HitEnergyPosXY_Si1.size()) -
                                                      static_cast<int>(HitEnergyPosXYreal_Si1.size()),
                                                  HitEnergyPosXYreal_Si1.size(), 1.);

  LocalHisto.h_nEventsGoodrecons_Si1->Fill(nGoodrecons_Si1, 1.);
  LocalHisto.h_nEventsGhost_Si1->Fill(nGhost_Si1, 1.);
  LocalHisto.h_nEventsGoodreconsGhost_Si1->Fill(nGoodrecons_Si1, nGhost_Si1, 1.);
  LocalHisto.h_nEventsRealGoodrecons_Si1->Fill(HitEnergyPosXYreal_Si1.size(), nGoodrecons_Si1, 1.);
  LocalHisto.h_nEventsRealRejectPad_Si1->Fill(HitEnergyPosXYreal_Si1.size(), SiliconHitsSD_Si1.CountRejectPad(), 1.);
  SiliconHitsSD_Si1.Clear();

  size_t nGoodrecons_Si2 = 0;
  nGoodEventsCounter(HitEnergyPosXY_Si2, HitEnergyPosXYreal_Si2, widthStrip_Si2, nGoodrecons_Si2);
  int nGhost_Si2    = HitEnergyPosXY_Si2.size() - nGoodrecons_Si2;

  LocalHisto.h_HitMultiplicity_Si2->Fill(HitEnergyPosXYreal_Si2.size(), 1.);
  LocalHisto.h_HitMultiplicityRecons_Si2->Fill(HitEnergyPosXY_Si2.size(), 1.);
  LocalHisto.h_HitMultiplicityDiff_Si2->Fill(
      static_cast<int>(HitEnergyPosXY_Si2.size()) - static_cast<int>(HitEnergyPosXYreal_Si2.size()), 1.);
  LocalHisto.h_HitMultiplicityDiffNHits_Si2->Fill(static_cast<int>(HitEnergyPosXY_Si2.size()) -
                                                      static_cast<int>(HitEnergyPosXYreal_Si2.size()),
                                                  HitEnergyPosXYreal_Si2.size(), 1.);

  LocalHisto.h_nEventsGoodrecons_Si2->Fill(nGoodrecons_Si2, 1.);
  LocalHisto.h_nEventsGhost_Si2->Fill(nGhost_Si2, 1.);
  LocalHisto.h_nEventsGoodreconsGhost_Si2->Fill(nGoodrecons_Si2, nGhost_Si2, 1.);
  LocalHisto.h_nEventsRealGoodrecons_Si2->Fill(HitEnergyPosXYreal_Si2.size(), nGoodrecons_Si2, 1.);
  LocalHisto.h_nEventsRealRejectPad_Si2->Fill(HitEnergyPosXYreal_Si2.size(), SiliconHitsSD_Si2.CountRejectPad(), 1.);
  SiliconHitsSD_Si2.Clear();

#endif

  double InteractionPoint_real_X = RecoEvent.InteractionPoint[0];
  double InteractionPoint_real_Y = RecoEvent.InteractionPoint[1];
  double InteractionPoint_real_Z = RecoEvent.InteractionPoint[2];

  double InteractionPoint_rand_X =
      rand->Uniform(InteractionPoint_real_X - randInteractionPointXY, InteractionPoint_real_X + randInteractionPointXY);
  double InteractionPoint_rand_Y =
      rand->Uniform(InteractionPoint_real_Y - randInteractionPointXY, InteractionPoint_real_Y + randInteractionPointXY);

  double EnergyBeamHit1 = 0.1;  // random
  double EnergyBeamHit2 = 0.01; // random

  std::vector<double> BeamHit1{EnergyBeamHit1, InteractionPoint_rand_X, InteractionPoint_rand_Y,
                               InteractionPoint_real_Z};
  std::vector<double> BeamHit2{EnergyBeamHit2, InteractionPoint_rand_X, InteractionPoint_rand_Y,
                               InteractionPoint_real_Z + 1.};

  // Track finding

  std::vector<std::vector<std::vector<double> > > CandidateTracks{};
  HitstoTracks(HitEnergyPosXY_Si1, HitEnergyPosXY_Si2, BeamHit1, BeamHit2, CandidateTracks);
  
  if(Si_4layers)
    HitstoTracks(HitEnergyPosXY_Si3, HitEnergyPosXY_Si4, BeamHit1, BeamHit2, CandidateTracks);


#ifdef TRACK_RECONS_CHECK

  // Real tracks
  std::vector<std::vector<std::vector<double> > > RealTracks{};
  RealHitstoRealTracks(HitEnergyPosXYreal_Si1, HitEnergyPosXYreal_Si2, RealTracks);
  
  if(Si_4layers)
    RealHitstoRealTracks(HitEnergyPosXYreal_Si3, HitEnergyPosXYreal_Si4, RealTracks);

  LocalHisto.h_nCandidatesRealTracks->Fill(CandidateTracks.size(), nRealTracksCounter, 1.);
#endif

#ifdef MOTHER_DAUGHTERS_CHECK

  auto it_mother = RecoEvent.TrackMother.begin();
  if(it_mother == RecoEvent.TrackMother.end())
    return -3;

  int mother_trackID = std::get<0>(it_mother->second);
  
  std::vector<int> daughters_TrackID {};
  for(auto it_daughters = RecoEvent.TrackMother.begin(); it_daughters != RecoEvent.TrackMother.end(); ++it_daughters)
    {
      daughters_TrackID.emplace_back(it_daughters->first);
    }

  std::vector<std::vector<std::vector<double>>> motherTracks{};
  std::vector<std::vector<std::vector<double>>> daughtersTracks{};

  for(size_t i = 0; i < RealTracks.size(); ++i)
    {
      if(std::abs(RealTracks[i][0][7] - mother_trackID) < 0.01)
        motherTracks.emplace_back(RealTracks[i]);

      for(size_t j = 0; j < daughters_TrackID.size(); ++j)
        {
          if(std::abs(RealTracks[i][0][7] - daughters_TrackID[j]) < 0.01)
            daughtersTracks.emplace_back(RealTracks[i]);
        }
    }


  if(motherTracks.size() == 1)
    {
      LocalHisto.h_EnergyDepositionMother->Fill(motherTracks[0][0][0] + motherTracks[0][1][0], "Total", 1.);
      LocalHisto.h_EnergyDepositionMother->Fill(motherTracks[0][0][4], "Si1 X", 1.);
      LocalHisto.h_EnergyDepositionMother->Fill(motherTracks[0][0][5], "Si1 Y", 1.);
      LocalHisto.h_EnergyDepositionMother->Fill(motherTracks[0][1][4], "Si2 X", 1.);
      LocalHisto.h_EnergyDepositionMother->Fill(motherTracks[0][1][5], "Si2 Y", 1.);

      for(size_t i = 0; i < CandidateTracks.size(); ++i)
        {
          if((std::abs(motherTracks[0][0][1] - CandidateTracks[i][0][1]) < widthStrip_Si1 / 2.) &&
             (std::abs(motherTracks[0][0][2] - CandidateTracks[i][0][2]) < widthStrip_Si1 / 2.) &&
             (std::abs(motherTracks[0][1][1] - CandidateTracks[i][1][1]) < widthStrip_Si2 / 2.) &&
             (std::abs(motherTracks[0][1][2] - CandidateTracks[i][1][2]) < widthStrip_Si2 / 2.))
            {
              LocalHisto.h_nHypernucleiTrack->Fill("Recons H3L", 1.);
            }
        }

      double distance = 0.;
      double z        = 0.;

      CloseDist(BeamHit1, BeamHit2, motherTracks[0][0], motherTracks[0][1], distance, z);

      LocalHisto.h_DistanceBeamTracks->Fill(distance, "Mother", 1.);
      LocalHisto.h_PosZBeamTracks->Fill(z, "Mother", 1.);
    }

  for(size_t j = 0; j < daughtersTracks.size(); ++j)
    {
      LocalHisto.h_EnergyDepositionDaughters->Fill(daughtersTracks[j][0][0] + daughtersTracks[j][1][0], "Total", 1.);
      LocalHisto.h_EnergyDepositionDaughters->Fill(daughtersTracks[j][0][4], "Si1 X", 1.);
      LocalHisto.h_EnergyDepositionDaughters->Fill(daughtersTracks[j][0][5], "Si1 Y", 1.);
      LocalHisto.h_EnergyDepositionDaughters->Fill(daughtersTracks[j][1][4], "Si2 X", 1.);
      LocalHisto.h_EnergyDepositionDaughters->Fill(daughtersTracks[j][1][5], "Si2 Y", 1.);


      for(size_t i = 0; i < CandidateTracks.size(); ++i)
        {
          if((std::abs(daughtersTracks[j][0][1] - CandidateTracks[i][0][1]) < widthStrip_Si1 / 2.) &&
             (std::abs(daughtersTracks[j][0][2] - CandidateTracks[i][0][2]) < widthStrip_Si1 / 2.) &&
             (std::abs(daughtersTracks[j][1][1] - CandidateTracks[i][1][1]) < widthStrip_Si2 / 2.) &&
             (std::abs(daughtersTracks[j][1][2] - CandidateTracks[i][1][2]) < widthStrip_Si2 / 2.))
            {
              if(daughtersTracks[j][0][6] == 2)
                {
                  LocalHisto.h_nHypernucleiTrack->Fill("Recons He3", 1.);
                }

              else if(daughtersTracks[j][0][6] == 3)
                {
                  LocalHisto.h_nHypernucleiTrack->Fill("Recons pi-", 1.);
                }
            }
        }

      double distance = 0.;
      double z        = 0.;

      CloseDist(BeamHit1, BeamHit2, daughtersTracks[j][0], daughtersTracks[j][1], distance, z);

      LocalHisto.h_DistanceBeamTracks->Fill(distance, "Daughters", 1.);
      LocalHisto.h_PosZBeamTracks->Fill(z, "Daughters", 1.);

      double thetadaughters = std::atan(std::sqrt(std::pow((daughtersTracks[j][1][2] - daughtersTracks[j][0][2]), 2.) +
                                        std::pow((daughtersTracks[j][1][1] - daughtersTracks[j][0][1]), 2.)) /
                                   (Z_plane_Si2 - Z_plane_Si1)) * 180. / M_PI;
      LocalHisto.h_thetaTracks->Fill(thetadaughters, "Daughters", 1.);
    }

#endif

#ifdef TRACK_RECONS_CHECK

  for(size_t i = 0; i < CandidateTracks.size(); ++i)
    {
      double trackPosX_Si1 = CandidateTracks[i][0][1];
      double trackPosY_Si1 = CandidateTracks[i][0][2];

      double trackPosX_Si2 = CandidateTracks[i][1][1];
      double trackPosY_Si2 = CandidateTracks[i][1][2];

      double TotalEnergy = CandidateTracks[i][0][0] + CandidateTracks[i][1][0];

      for(size_t j = 0; j < HitEnergyLayerX_Si1.size(); ++j)
        {
          double HitX_Si1 = -lenghtSi_Si1 / 2. + (get<1>(HitEnergyLayerX_Si1[j])%nStrips_Si1 + 0.5) * widthStrip_Si1;

          if(std::abs(trackPosX_Si1 - HitX_Si1) <= widthStrip_Si1 / 2.)
            {
              for(size_t k = 0; k < HitEnergyLayerY_Si1.size(); ++k)
                {
                  double HitY_Si1 = +lenghtSi_Si1 / 2. - (get<1>(HitEnergyLayerY_Si1[k])%nStrips_Si1 + 0.5) * widthStrip_Si1;

                  if(std::abs(trackPosY_Si1 - HitY_Si1) <= widthStrip_Si1 / 2.)
                    {
                      for(size_t l = 0; l < HitEnergyLayerX_Si2.size(); ++l)
                        {
                          double HitX_Si2 =
                              -lenghtSi_Si2 / 2. + (get<1>(HitEnergyLayerX_Si2[l])%nStrips_Si2 + 0.5) * widthStrip_Si2;

                          if(std::abs(trackPosX_Si2 - HitX_Si2) <= widthStrip_Si2 / 2.)
                            {
                              for(size_t m = 0; m < HitEnergyLayerY_Si2.size(); ++m)
                                {
                                  double HitY_Si2 =
                                      +lenghtSi_Si2 / 2. - (get<1>(HitEnergyLayerY_Si2[m])%nStrips_Si2 + 0.5) * widthStrip_Si2;

                                  if(std::abs(trackPosY_Si2 - HitY_Si2) <= widthStrip_Si2 / 2.)
                                    {
                                      LocalHisto.h_EnergyStripEnergyTotal->Fill(get<0>(HitEnergyLayerX_Si1[j]),
                                                                                TotalEnergy, 1.);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

  for(size_t i = 0; i < RealTracks.size(); ++i)
    {
      double distance = 0.;
      double z        = 0.;

      CloseDist(BeamHit1, BeamHit2, RealTracks[i][0], RealTracks[i][1], distance, z);

      LocalHisto.h_DistanceBeamTracks->Fill(distance, "Real track", 1.);
      LocalHisto.h_PosZBeamTracks->Fill(z, "Real track", 1.);

      LocalHisto.h_EnergyStripEnergyTotalReal->Fill(RealTracks[i][0][4], RealTracks[i][0][0] + RealTracks[i][1][0], 1.);
    }

  // Good tracks
  size_t nGoodTracks = 0;
  std::vector<size_t> goodCandidateTracks(CandidateTracks.size(), 0);
  nGoodTracksCounter(CandidateTracks, RealTracks, nGoodTracks, goodCandidateTracks);

  size_t nForwardTracks = 0;
  std::vector<size_t> forwardCandidateTracks(CandidateTracks.size(), 0);
  nForwardTracksCounter(CandidateTracks, nForwardTracks, forwardCandidateTracks);

  LocalHisto.h_nTrackCandidates->Fill(CandidateTracks.size(), "Total", 1.);
  LocalHisto.h_nTrackCandidates->Fill(nGoodTracks, "Good", 1.);
  LocalHisto.h_nTrackCandidates->Fill(nForwardTracks, "Forward", 1.);

#endif

  std::vector<double> InteractionPointRecons(3, 0.);
  std::vector<double> InteractionPointAverage(3, 0.);

  if(CandidateTracks.size() == 0)
    return -1;

  std::vector<double> f_values_IP(CandidateTracks.size() + 1, 0.);

  TrackstoVertexPosition(CandidateTracks, BeamHit1, BeamHit2, InteractionPointRecons, f_values_IP);
  RecoEvent.PrimVtxRecons.SetXYZ(InteractionPointRecons[0],InteractionPointRecons[1],InteractionPointRecons[2]);

#ifdef VERTEX_RECONS_CHECK

  LocalHisto.h_nCandidatesRealTracks_IfRecons->Fill(CandidateTracks.size(), nRealTracksCounter, 1.);

  for(size_t i = 0; i < f_values_IP.size(); ++i)
    LocalHisto.h_fvalues->Fill(f_values_IP[i], 1.);

  double distance  = std::sqrt(std::pow((InteractionPoint_real_X - InteractionPointRecons[0]), 2.) +
                          std::pow((InteractionPoint_real_Y - InteractionPointRecons[1]), 2.) +
                          std::pow((InteractionPoint_real_Z - InteractionPointRecons[2]), 2.));
  double distanceX = InteractionPoint_real_X - InteractionPointRecons[0];
  double distanceY = InteractionPoint_real_Y - InteractionPointRecons[1];
  double distanceZ = InteractionPoint_real_Z - InteractionPointRecons[2];

  LocalHisto.h_InteractionPointDistance->Fill(distance, 1.);
  LocalHisto.h_InteractionPointDistanceX->Fill(distanceX, 1.);
  LocalHisto.h_InteractionPointDistanceY->Fill(distanceY, 1.);
  LocalHisto.h_InteractionPointDistanceZ->Fill(distanceZ, 1.);

#endif

#ifdef COVARIANCE_MATRIX

  std::vector<std::vector<double> > CovMatrix;
  CovarianceMatrix(CandidateTracks, BeamHit1, BeamHit2, InteractionPointAverage, f_values_IP, CovMatrix);

  RecoEvent.CovMatrix_IP = {CovMatrix[0][0],
                            CovMatrix[1][0], CovMatrix[1][1],
                            CovMatrix[2][0], CovMatrix[2][1], CovMatrix[2][2]};

  LocalHisto.h_InteractionPointDistanceX_pull->Fill(distanceX / std::sqrt(TMath::Max(CovMatrix[0][0], 1e-10)), 1.);
  LocalHisto.h_InteractionPointDistanceY_pull->Fill(distanceY / std::sqrt(TMath::Max(CovMatrix[1][1], 1e-10)), 1.);
  LocalHisto.h_InteractionPointDistanceZ_pull->Fill(distanceZ / std::sqrt(TMath::Max(CovMatrix[2][2], 1e-10)), 1.);


  LocalHisto.h_CovarianceSigmaX->Fill(std::sqrt(CovMatrix[0][0]), 1.);
  LocalHisto.h_CovarianceSigmaY->Fill(std::sqrt(CovMatrix[1][1]), 1.);
  LocalHisto.h_CovarianceSigmaZ->Fill(std::sqrt(CovMatrix[2][2]), 1.);

#endif

#ifdef DECAY_VERTEX
  
  double DecayPosition_real_X = std::get<1>(it_mother->second);
  double DecayPosition_real_Y = std::get<2>(it_mother->second);
  double DecayPosition_real_Z = std::get<3>(it_mother->second);

  double IP_DecayDistanceX = DecayPosition_real_X - InteractionPoint_real_X;
  double IP_DecayDistanceY = DecayPosition_real_Y - InteractionPoint_real_Y;
  double IP_DecayDistanceZ = DecayPosition_real_Z - InteractionPoint_real_Z;

  LocalHisto.h_IP_DecayDistance->Fill(
      std::sqrt(std::pow(IP_DecayDistanceX, 2.) + std::pow(IP_DecayDistanceY, 2.) + std::pow(IP_DecayDistanceZ, 2.)), 1.);
  LocalHisto.h_IP_DecayDistanceX->Fill(IP_DecayDistanceX, 1.);
  LocalHisto.h_IP_DecayDistanceY->Fill(IP_DecayDistanceY, 1.);
  LocalHisto.h_IP_DecayDistanceZ->Fill(IP_DecayDistanceZ, 1.);

  std::vector<std::vector<std::vector<double> > > CandidateDecayTracks{};

  HitstoDecayTracks(HitEnergyPosXY_Si1, HitEnergyPosXY_Si2, BeamHit1, BeamHit2, InteractionPointRecons,
                    CandidateDecayTracks);
  LocalHisto.h_nTrackCandidates->Fill(CandidateDecayTracks.size(), "Decay", 1.);


  if(CandidateDecayTracks.size() == 0)
    {
      return -2;
    }

  std::vector<double> DecayPositionRecons(3, 0.);
  DecayTrackstoDecayPosition(CandidateDecayTracks, InteractionPointRecons, DecayPositionRecons);

  double distanceDecayPosition  = std::sqrt(std::pow((DecayPosition_real_X - DecayPositionRecons[0]), 2.) +
                                      std::pow((DecayPosition_real_Y - DecayPositionRecons[1]), 2.) +
                                      std::pow((DecayPosition_real_Z - DecayPositionRecons[2]), 2.));
  double distanceDecayPositionX = DecayPosition_real_X - DecayPositionRecons[0];
  double distanceDecayPositionY = DecayPosition_real_Y - DecayPositionRecons[1];
  double distanceDecayPositionZ = DecayPosition_real_Z - DecayPositionRecons[2];

  LocalHisto.h_DecayPositionDistance->Fill(distanceDecayPosition, 1.);
  LocalHisto.h_DecayPositionDistanceX->Fill(distanceDecayPositionX, 1.);
  LocalHisto.h_DecayPositionDistanceY->Fill(distanceDecayPositionY, 1.);
  LocalHisto.h_DecayPositionDistanceZ->Fill(distanceDecayPositionZ, 1.);

#ifdef MOTHER_DAUGHTERS_CHECK

  if(motherTracks.size() == 1)
    {
      for(size_t i = 0; i < CandidateDecayTracks.size(); ++i)
        {
          if((std::abs(motherTracks[0][0][1] - CandidateDecayTracks[i][0][1]) < widthStrip_Si1 / 2.) &&
             (std::abs(motherTracks[0][0][2] - CandidateDecayTracks[i][0][2]) < widthStrip_Si1 / 2.) &&
             (std::abs(motherTracks[0][1][1] - CandidateDecayTracks[i][1][1]) < widthStrip_Si2 / 2.) &&
             (std::abs(motherTracks[0][1][2] - CandidateDecayTracks[i][1][2]) < widthStrip_Si2 / 2.))
            {
              LocalHisto.h_nHypernucleiTrack->Fill("Decay H3L", 1.);
            }
        }
    }

  for(size_t j = 0; j < daughtersTracks.size(); ++j)
    {
      for(size_t i = 0; i < CandidateDecayTracks.size(); ++i)
        {
          if((std::abs(daughtersTracks[j][0][1] - CandidateDecayTracks[i][0][1]) < widthStrip_Si1 / 2.) &&
             (std::abs(daughtersTracks[j][0][2] - CandidateDecayTracks[i][0][2]) < widthStrip_Si1 / 2.) &&
             (std::abs(daughtersTracks[j][1][1] - CandidateDecayTracks[i][1][1]) < widthStrip_Si2 / 2.) &&
             (std::abs(daughtersTracks[j][1][2] - CandidateDecayTracks[i][1][2]) < widthStrip_Si2 / 2.))
            {
              if(daughtersTracks[j][0][6] == 2)
                {
                  LocalHisto.h_nHypernucleiTrack->Fill("Decay He3", 1.);
                }

              else if(daughtersTracks[j][0][6] == 3)
                {
                  LocalHisto.h_nHypernucleiTrack->Fill("Decay pi-", 1.);
                }
            }
        }
    }

#endif

#endif

  return 0;
}

template <class Out>
void TPrimaryVertex<Out>::simulHitstoRealHits(
    FullRecoEvent& REvent,
    std::vector<std::tuple<double, double, double, size_t, double, double, std::string> >& HitEnergyPosXYreal,
    int id_det_x, int id_det_y, TH1F* h_Diff, std::vector<std::tuple<size_t,TVector3,TVector3>>& HitIdMomPos)
{
  for(auto it_track : REvent.TrackDAFSim)
    {
      auto it_SiX = it_track.second[id_det_x];
      auto it_SiY = it_track.second[id_det_y];

      if(it_SiX.size() == 0 || it_SiY.size() == 0)
        continue;

      double meanPosX = 0.;
      double EnergyX  = 0.;
      double meanPosY = 0.;
      double EnergyY  = 0.;

      double meanPosYinXstrips = 0.;
      double meanPosZinXstrips = 0.;

      double meanPosXinYstrips = 0.;
      double meanPosZinYstrips = 0.;

      for(size_t idX = 0; idX < it_SiX.size(); ++idX)
        {
          if(it_SiX[idX].Eloss > EnergyThreshold)
            {
              meanPosX += it_SiX[idX].hitX * it_SiX[idX].Eloss;
              EnergyX += it_SiX[idX].Eloss;

              meanPosYinXstrips += it_SiX[idX].hitY * it_SiX[idX].Eloss;
              meanPosZinXstrips += it_SiX[idX].hitZ * it_SiX[idX].Eloss;
            }
        }

      for(size_t idY = 0; idY < it_SiY.size(); ++idY)
        {
          if(it_SiY[idY].Eloss > EnergyThreshold)
            {
              meanPosY += it_SiY[idY].hitY * it_SiY[idY].Eloss;
              EnergyY += it_SiY[idY].Eloss;

              meanPosXinYstrips += it_SiY[idY].hitX * it_SiY[idY].Eloss;
              meanPosZinYstrips += it_SiY[idY].hitZ * it_SiY[idY].Eloss;
            }
        }

      meanPosX /= EnergyX;
      meanPosY /= EnergyY;

      meanPosYinXstrips /= EnergyX;
      meanPosZinXstrips /= EnergyX;

      meanPosXinYstrips /= EnergyY;
      meanPosZinYstrips /= EnergyY;

      std::string namePart(TDatabasePDG::Instance()->GetParticle(it_SiX[0].pdg)->GetName());

      std::tuple<double, double, double, size_t, double, double, std::string> tempEnergyPosXYreal(
          EnergyX + EnergyY, meanPosX, meanPosY, it_track.first, EnergyX, EnergyY, namePart);
      HitEnergyPosXYreal.emplace_back(tempEnergyPosXYreal);

      h_Diff->Fill(EnergyX - EnergyY, 1.);

      //Check magnetic field effects on silicon tracking

      TVector3 HitSiX(meanPosX, meanPosYinXstrips, meanPosZinXstrips);
      TVector3 HitSiY(meanPosXinYstrips, meanPosY, meanPosZinYstrips);

      TVector3 MomSi = HitSiY - HitSiX;
      std::tuple<size_t,TVector3,TVector3> tempHitIdMomPos(it_track.first, MomSi, HitSiX);
      HitIdMomPos.emplace_back(tempHitIdMomPos);
    }
}

template <class Out>
void TPrimaryVertex<Out>::MFcheck(std::vector<std::tuple<size_t,TVector3,TVector3>>& HitIdMomPos_Si1,
                             std::vector<std::tuple<size_t,TVector3,TVector3>>& HitIdMomPos_Si2,
                             std::array<double,3> InteractionPoint)
{
  TVector3 Real_IP(InteractionPoint[0], InteractionPoint[1], InteractionPoint[2]);

  for(size_t i = 0; i < HitIdMomPos_Si1.size(); ++i)
    {
      size_t temp_idSi1 = get<0>(HitIdMomPos_Si1[i]);

      for(size_t j = 0; j < HitIdMomPos_Si2.size(); ++j)
        {
          size_t temp_idSi2 = get<0>(HitIdMomPos_Si2[j]);
          if (temp_idSi1 != temp_idSi2)
            continue;

          TVector3 temp_MomSi1 = get<1>(HitIdMomPos_Si1[i]);
          TVector3 temp_PosSi1 = get<2>(HitIdMomPos_Si1[i]);

          TVector3 temp_MomSi2 = get<1>(HitIdMomPos_Si2[j]);
          TVector3 temp_PosSi2 = get<2>(HitIdMomPos_Si2[j]);

          //std::cout << "MomSi1: " << temp_MomSi1.X()/temp_MomSi1.Z() << "\t" << temp_MomSi1.Y()/temp_MomSi1.Z() << "\n";
          //std::cout << "MomSi2: " << temp_MomSi2.X()/temp_MomSi2.Z() << "\t" << temp_MomSi2.Y()/temp_MomSi2.Z() << "\n\n";

          double Theta_MomSi1MomSi2 = temp_MomSi1.Angle(temp_MomSi2) * 180. / M_PI;
          LocalHisto.h_MFCheck_Theta_MomSi1MomSi2->Fill(Theta_MomSi1MomSi2, 1.);

          TVector3 PosSi1_PosSi2 = temp_PosSi1 - temp_PosSi2;
          double Dist_MomSi1HitSi2 = (PosSi1_PosSi2.Cross(temp_MomSi1)).Mag() / temp_MomSi1.Mag();
          LocalHisto.h_MFCheck_Dist_MomSi1HitSi2->Fill(Dist_MomSi1HitSi2, 1.);

          TVector3 PosSi2_PosSi1 = temp_PosSi2 - temp_PosSi1;
          double Dist_MomSi2HitSi1 = (PosSi2_PosSi1.Cross(temp_MomSi2)).Mag() / temp_MomSi2.Mag();
          LocalHisto.h_MFCheck_Dist_MomSi2HitSi1->Fill(Dist_MomSi2HitSi1, 1.);

          TVector3 PosSi1_IP = temp_PosSi1 - Real_IP;
          double Dist_MomSi1HitIP = (PosSi1_IP.Cross(temp_MomSi1)).Mag() / temp_MomSi1.Mag();
          LocalHisto.h_MFCheck_Dist_MomSi1HitIP->Fill(Dist_MomSi1HitIP, 1.);

          TVector3 PosSi2_IP = temp_PosSi2 - Real_IP;
          double Dist_MomSi2HitIP = (PosSi2_IP.Cross(temp_MomSi2)).Mag() / temp_MomSi2.Mag();
          LocalHisto.h_MFCheck_Dist_MomSi2HitIP->Fill(Dist_MomSi2HitIP, 1.);
        }
    }
}

/*
void SiliconHits::SignalstoHits(std::vector<std::tuple<double, size_t> >& HitEnergyLayerX,
                                std::vector<std::tuple<double, size_t> >& HitEnergyLayerY,
                                std::vector<std::vector<double> >& HitEnergyPosXY)
{
  std::vector<double> countHitLayerX(HitEnergyLayerX.size(), 0.);
  std::vector<double> countHitLayerY(HitEnergyLayerY.size(), 0.);

  // Events with multiplicityX = multiplicityY = 1
  if((HitEnergyLayerX.size() > 0.1) && (HitEnergyLayerY.size() > 0.1))
    {
      for(size_t i = 0; i < HitEnergyLayerX.size(); ++i)
        {
          for(size_t j = 0; j < HitEnergyLayerY.size(); ++j)
            {
              double EnergyDiff = std::abs(get<0>(HitEnergyLayerX[i]) - get<0>(HitEnergyLayerY[j]));
              double HitX       = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[i]) + 0.5) * widthStrip;
              double HitY       = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[j]) + 0.5) * widthStrip;

              if(EnergyDiff < MaxEnergyDiffStrips)
                {
                  double HitEnergy = get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerY[j]);
                  std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                  HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                  countHitLayerX[i] += 1.;
                  countHitLayerY[j] += 1.;
                }
            }
        }
    }

#ifdef RECONS_HITS_MULTIPLICITY

  for(int i = HitEnergyLayerX.size() - 1; i >= 0; --i)
    {
      if(countHitLayerX[i] > 0.1)
        {
          HitEnergyLayerX.erase(HitEnergyLayerX.begin() + i);
          countHitLayerX.erase(countHitLayerX.begin() + i);
        }
    }

  for(int j = HitEnergyLayerY.size() - 1; j >= 0; --j)
    {
      if(countHitLayerY[j] > 0.1)
        {
          HitEnergyLayerY.erase(HitEnergyLayerY.begin() + j);
          countHitLayerY.erase(countHitLayerY.begin() + j);
        }
    }

  // Events with multiplicityX = 1 and multiplicityY = 2
  if((HitEnergyLayerX.size() > 1.1) && (HitEnergyLayerY.size() > 0.1))
    {
      for(size_t i = 0; i < HitEnergyLayerX.size() - 1; ++i)
        {
          for(size_t j = i + 1; j < HitEnergyLayerX.size(); ++j)
            {
              for(size_t k = 0; k < HitEnergyLayerY.size(); ++k)
                {
                  double EnergyDiff =
                      std::abs(get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j]) - get<0>(HitEnergyLayerY[k]));
                  if(EnergyDiff < MaxEnergyDiffStrips)
                    {
                      if((size_type_abs(get<1>(HitEnergyLayerX[i]), get<1>(HitEnergyLayerX[j])) == 1) &&
                         ((get<0>(HitEnergyLayerX[i]) < MaxEnergyMultiplicity) ||
                          (get<0>(HitEnergyLayerX[j]) < MaxEnergyMultiplicity)))
                        {
                          double HitX =
                              -lenghtSi / 2. + ((get<0>(HitEnergyLayerX[i]) * get<1>(HitEnergyLayerX[i]) +
                                                 get<0>(HitEnergyLayerX[j]) * get<1>(HitEnergyLayerX[j])) /
                                                    (get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j])) +
                                                0.5) *
                                                   widthStrip;
                          double HitY = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[k]) + 0.5) * widthStrip;
                          double HitEnergy =
                              get<0>(HitEnergyLayerY[k]) + get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j]);
                          std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                          countHitLayerX[i] += 0.5;
                          countHitLayerX[j] += 0.5;
                          countHitLayerY[k] += 1.;
                        }
                      else
                        {
                          double HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[i]) + 0.5) * widthStrip;
                          double HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[k]) + 0.5) * widthStrip;
                          double HitEnergy = get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerY[k]) / 2.;
                          std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                          HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[j]) + 0.5) * widthStrip;
                          HitEnergy = get<0>(HitEnergyLayerX[j]) + get<0>(HitEnergyLayerY[k]) / 2.;
                          std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                          countHitLayerX[i] += 1.;
                          countHitLayerX[j] += 1.;
                          countHitLayerY[k] += 2.;
                        }
                    }
                }
            }
        }

      for(int i = HitEnergyLayerX.size() - 1; i >= 0; --i)
        {
          if(countHitLayerX[i] > 0.1)
            {
              HitEnergyLayerX.erase(HitEnergyLayerX.begin() + i);
              countHitLayerX.erase(countHitLayerX.begin() + i);
            }
        }

      for(int j = HitEnergyLayerY.size() - 1; j >= 0; --j)
        {
          if(countHitLayerY[j] > 0.1)
            {
              HitEnergyLayerY.erase(HitEnergyLayerY.begin() + j);
              countHitLayerY.erase(countHitLayerY.begin() + j);
            }
        }
    }

  // Events with multiplicityX = 2 and multiplicityY = 1
  if((HitEnergyLayerX.size() > 0.1) && (HitEnergyLayerY.size() > 1.1))
    {
      for(size_t i = 0; i < HitEnergyLayerY.size() - 1; ++i)
        {
          for(size_t j = i + 1; j < HitEnergyLayerY.size(); ++j)
            {
              for(size_t k = 0; k < HitEnergyLayerX.size(); ++k)
                {
                  double EnergyDiff =
                      std::abs(get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j]) - get<0>(HitEnergyLayerX[k]));
                  if(EnergyDiff < MaxEnergyDiffStrips)
                    {
                      if((size_type_abs(get<1>(HitEnergyLayerY[i]), get<1>(HitEnergyLayerY[j])) == 1) &&
                         ((get<0>(HitEnergyLayerY[i]) < MaxEnergyMultiplicity) ||
                          (get<0>(HitEnergyLayerY[j]) < MaxEnergyMultiplicity)))
                        {
                          double HitX = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[k]) + 0.5) * widthStrip;
                          double HitY = -lenghtSi / 2. + ((get<0>(HitEnergyLayerY[i]) * get<1>(HitEnergyLayerY[i]) +
                                                          get<0>(HitEnergyLayerY[j]) * get<1>(HitEnergyLayerY[j])) /
                                                             (get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j])) +
                                                         0.5) *
                                                            widthStrip;
                          double HitEnergy =
                              get<0>(HitEnergyLayerX[k]) + get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j]);
                          std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                          countHitLayerY[i] += 0.5;
                          countHitLayerY[j] += 0.5;
                          countHitLayerX[k] += 1.;
                        }
                      else
                        {
                          double HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[k]) + 0.5) * widthStrip;
                          double HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[i]) + 0.5) * widthStrip;
                          double HitEnergy = get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerX[k]) / 2.;
                          std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                          HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[j]) + 0.5) * widthStrip;
                          HitEnergy = get<0>(HitEnergyLayerY[j]) + get<0>(HitEnergyLayerX[k]) / 2.;
                          std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                          countHitLayerY[i] += 1.;
                          countHitLayerY[j] += 1.;
                          countHitLayerX[k] += 2.;
                        }
                    }
                }
            }
        }

      for(int i = HitEnergyLayerX.size() - 1; i >= 0; --i)
        {
          if(countHitLayerX[i] > 0.1)
            {
              HitEnergyLayerX.erase(HitEnergyLayerX.begin() + i);
              countHitLayerX.erase(countHitLayerX.begin() + i);
            }
        }

      for(int j = HitEnergyLayerY.size() - 1; j >= 0; --j)
        {
          if(countHitLayerY[j] > 0.1)
            {
              HitEnergyLayerY.erase(HitEnergyLayerY.begin() + j);
              countHitLayerY.erase(countHitLayerY.begin() + j);
            }
        }
    }

  // Events with multiplicityX = 1 and multiplicityY = 3
  if((HitEnergyLayerX.size() > 2.1) && (HitEnergyLayerY.size() > 0.1))
    {
      for(size_t i = 0; i < HitEnergyLayerX.size() - 2; ++i)
        {
          for(size_t j = i + 1; j < HitEnergyLayerX.size() - 1; ++j)
            {
              for(size_t k = j + 1; k < HitEnergyLayerX.size(); ++k)
                {
                  for(size_t l = 0; l < HitEnergyLayerY.size(); ++l)
                    {
                      double EnergyDiff = std::abs(get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j]) +
                                              get<0>(HitEnergyLayerX[k]) - get<0>(HitEnergyLayerY[l]));
                      if(EnergyDiff < MaxEnergyDiffStrips)
                        {
                          if((size_type_abs(get<1>(HitEnergyLayerX[i]), get<1>(HitEnergyLayerX[j])) == 1) &&
                             ((get<0>(HitEnergyLayerX[i]) < MaxEnergyMultiplicity) ||
                              (get<0>(HitEnergyLayerX[j]) < MaxEnergyMultiplicity)))
                            {
                              double HitX =
                                  -lenghtSi / 2. + ((get<0>(HitEnergyLayerX[i]) * get<1>(HitEnergyLayerX[i]) +
                                                     get<0>(HitEnergyLayerX[j]) * get<1>(HitEnergyLayerX[j])) /
                                                        (get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j])) +
                                                    0.5) *
                                                       widthStrip;
                              double HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[l]) + 0.5) * widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j]) +
                                                 get<0>(HitEnergyLayerY[l]) / 2.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[k]) + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerX[k]) + get<0>(HitEnergyLayerY[l]) / 2.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              countHitLayerX[i] += 0.5;
                              countHitLayerX[j] += 0.5;
                              countHitLayerX[k] += 1.;
                              countHitLayerY[l] += 2.;
                            }
                          else if((size_type_abs(get<1>(HitEnergyLayerX[j]), get<1>(HitEnergyLayerX[k])) == 1) &&
                                  ((get<0>(HitEnergyLayerX[j]) < MaxEnergyMultiplicity) ||
                                   (get<0>(HitEnergyLayerX[k]) < MaxEnergyMultiplicity)))
                            {
                              double HitX =
                                  -lenghtSi / 2. + ((get<0>(HitEnergyLayerX[j]) * get<1>(HitEnergyLayerX[j]) +
                                                     get<0>(HitEnergyLayerX[k]) * get<1>(HitEnergyLayerX[k])) /
                                                        (get<0>(HitEnergyLayerX[j]) + get<0>(HitEnergyLayerX[k])) +
                                                    0.5) *
                                                       widthStrip;
                              double HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[l]) + 0.5) * widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerX[j]) + get<0>(HitEnergyLayerX[k]) +
                                                 get<0>(HitEnergyLayerY[l]) / 2.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[i]) + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerY[l]) / 2.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              countHitLayerX[i] += 1.;
                              countHitLayerX[j] += 0.5;
                              countHitLayerX[k] += 0.5;
                              countHitLayerY[l] += 2.;
                            }
                          else
                            {
                              double HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[i]) + 0.5) * widthStrip;
                              double HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[l]) + 0.5) * widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerY[l]) / 3.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[j]) + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerX[j]) + get<0>(HitEnergyLayerY[l]) / 3.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[k]) + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerX[k]) + get<0>(HitEnergyLayerY[l]) / 3.;
                              std::vector<double> tempEnergyPosXY3{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY3);

                              countHitLayerX[i] += 1.;
                              countHitLayerX[j] += 1.;
                              countHitLayerX[k] += 1.;
                              countHitLayerY[l] += 3.;
                            }
                        }
                    }
                }
            }
        }

      for(int i = HitEnergyLayerX.size() - 1; i >= 0; --i)
        {
          if(countHitLayerX[i] > 0.1)
            {
              HitEnergyLayerX.erase(HitEnergyLayerX.begin() + i);
              countHitLayerX.erase(countHitLayerX.begin() + i);
            }
        }

      for(int j = HitEnergyLayerY.size() - 1; j >= 0; --j)
        {
          if(countHitLayerY[j] > 0.1)
            {
              HitEnergyLayerY.erase(HitEnergyLayerY.begin() + j);
              countHitLayerY.erase(countHitLayerY.begin() + j);
            }
        }
    }

  // Events with multiplicityX = 3 and multiplicityY = 1
  if((HitEnergyLayerX.size() > 0.1) && (HitEnergyLayerY.size() > 2.1))
    {
      for(size_t i = 0; i < HitEnergyLayerY.size() - 2; ++i)
        {
          for(size_t j = i + 1; j < HitEnergyLayerY.size() - 1; ++j)
            {
              for(size_t k = j + 1; k < HitEnergyLayerY.size(); ++k)
                {
                  for(size_t l = 0; l < HitEnergyLayerX.size(); ++l)
                    {
                      double EnergyDiff = std::abs(get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j]) +
                                              get<0>(HitEnergyLayerY[k]) - get<0>(HitEnergyLayerX[l]));
                      if(EnergyDiff < MaxEnergyDiffStrips)
                        {
                          if((size_type_abs(get<1>(HitEnergyLayerY[i]), get<1>(HitEnergyLayerY[j])) == 1) &&
                             ((get<0>(HitEnergyLayerY[i]) < MaxEnergyMultiplicity) ||
                              (get<0>(HitEnergyLayerY[j]) < MaxEnergyMultiplicity)))
                            {
                              double HitX = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[l]) + 0.5) * widthStrip;
                              double HitY =
                                  -lenghtSi / 2. + ((get<0>(HitEnergyLayerY[i]) * get<1>(HitEnergyLayerY[i]) +
                                                    get<0>(HitEnergyLayerY[j]) * get<1>(HitEnergyLayerY[j])) /
                                                       (get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j])) +
                                                   0.5) *
                                                      widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j]) +
                                                 get<0>(HitEnergyLayerX[l]) / 2.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[k]) + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerY[k]) + get<0>(HitEnergyLayerX[l]) / 2.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              countHitLayerY[i] += 0.5;
                              countHitLayerY[j] += 0.5;
                              countHitLayerY[k] += 1.;
                              countHitLayerX[l] += 2.;
                            }
                          else if((size_type_abs(get<1>(HitEnergyLayerY[j]), get<1>(HitEnergyLayerY[k])) == 1) &&
                                  ((get<0>(HitEnergyLayerY[j]) < MaxEnergyMultiplicity) ||
                                   (get<0>(HitEnergyLayerY[k]) < MaxEnergyMultiplicity)))
                            {
                              double HitX = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[l]) + 0.5) * widthStrip;
                              double HitY =
                                  -lenghtSi / 2. + ((get<0>(HitEnergyLayerY[j]) * get<1>(HitEnergyLayerY[j]) +
                                                    get<0>(HitEnergyLayerY[k]) * get<1>(HitEnergyLayerY[k])) /
                                                       (get<0>(HitEnergyLayerY[j]) + get<0>(HitEnergyLayerY[k])) +
                                                   0.5) *
                                                      widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerY[j]) + get<0>(HitEnergyLayerY[k]) +
                                                 get<0>(HitEnergyLayerX[l]) / 2.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[i]) + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerX[l]) / 2.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              countHitLayerY[i] += 1.;
                              countHitLayerY[j] += 0.5;
                              countHitLayerY[k] += 0.5;
                              countHitLayerX[l] += 2.;
                            }
                          else
                            {
                              double HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[l]) + 0.5) * widthStrip;
                              double HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[i]) + 0.5) * widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerX[l]) / 3.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[j]) + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerY[j]) + get<0>(HitEnergyLayerX[l]) / 3.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[k]) + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerY[k]) + get<0>(HitEnergyLayerX[l]) / 3.;
                              std::vector<double> tempEnergyPosXY3{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY3);

                              countHitLayerY[i] += 1.;
                              countHitLayerY[j] += 1.;
                              countHitLayerY[k] += 1.;
                              countHitLayerX[l] += 3.;
                            }
                        }
                    }
                }
            }
        }

      for(int i = HitEnergyLayerX.size() - 1; i >= 0; --i)
        {
          if(countHitLayerX[i] > 0.1)
            {
              HitEnergyLayerX.erase(HitEnergyLayerX.begin() + i);
              countHitLayerX.erase(countHitLayerX.begin() + i);
            }
        }

      for(int j = HitEnergyLayerY.size() - 1; j >= 0; --j)
        {
          if(countHitLayerY[j] > 0.1)
            {
              HitEnergyLayerY.erase(HitEnergyLayerY.begin() + j);
              countHitLayerY.erase(countHitLayerY.begin() + j);
            }
        }
    }

  if(HitEnergyPosXY.size() >= 2)
    {
      for(size_t i = 0; i < HitEnergyPosXY.size() - 1; ++i)
        {
          for(size_t j = i + 1; j < HitEnergyPosXY.size(); ++j)
            {
              if((HitEnergyPosXY[i][0] == HitEnergyPosXY[j][0]) && (HitEnergyPosXY[i][1] == HitEnergyPosXY[j][1]) &&
                 (HitEnergyPosXY[i][2] == HitEnergyPosXY[j][2]))
                {
                  HitEnergyPosXY.erase(HitEnergyPosXY.begin() + j);
                  --j;
                }
            }
        }
    }

#endif
}
*/
/*
bool SiliconHits_SD::IsValidCuadrant(size_t& layerX_id, size_t& layerY_id)
{
  //Cuadrant I
  if((layerX_id >= nStrips/2) && (layerX_id < nStrips) && (layerY_id >= 2*nStrips) && (layerY_id < 5*nStrips/2))
    return true;

  //Cuadrant II
  else if((layerX_id >= 0) && (layerX_id < nStrips/2) && (layerY_id >= 0) && (layerY_id < nStrips/2))
    return true;

  //Cuadrant III
  else if((layerX_id >= 2*nStrips) && (layerX_id < 5*nStrips/2) && (layerY_id >= nStrips/2) && (layerY_id < nStrips))
    return true;
    //Cuadrant I
  else if((layerX_id >= 5*nStrips/2) && (layerX_id < 3*nStrips) && (layerY_id >= 5*nStrips/2) && (layerY_id < 3*nStrips))
    return true;

  else
    {
      nRejectCuadrant += 1;
      return false;
    }
}

void SiliconHits_SD::CombineStrips(std::vector<std::tuple<double, size_t> >& HitEnergyLayer)
{
  if(combineStrips == 1)
    return;

  if((HitEnergyLayer.size() == 0) || (HitEnergyLayer.size() == 1))
    return;

  std::vector<std::tuple<double, size_t> > temp_HitEnergyLayer;

  std::vector<size_t> usedLayer {};
  std::vector<size_t>::iterator it_usedLayer;

  for(size_t i = 0; i < HitEnergyLayer.size()-1; ++i)
    {
      it_usedLayer = find(usedLayer.begin(), usedLayer.end(), get<1>(HitEnergyLayer[i]));
      if(it_usedLayer != usedLayer.end())
        continue;

      double temp_HitEnergy = get<0>(HitEnergyLayer[i]);
      size_t temp_HitLayer = (get<1>(HitEnergyLayer[i])%(nStrips*combineStrips))/combineStrips;
      
      if(get<1>(HitEnergyLayer[i]) >= 2*nStrips*combineStrips)
        temp_HitLayer += 2*nStrips;

      int temp_remainder = get<1>(HitEnergyLayer[i])%combineStrips;

      for(size_t j = i+1; j < HitEnergyLayer.size(); ++j)
        {
          it_usedLayer = find(usedLayer.begin(), usedLayer.end(), get<1>(HitEnergyLayer[j]));
          if(it_usedLayer != usedLayer.end())
            continue;

          int temp_diflayers = get<1>(HitEnergyLayer[i]) - get<1>(HitEnergyLayer[j]);
          if((temp_diflayers <= temp_remainder) && (temp_diflayers > temp_remainder - combineStrips))
            {
              temp_HitEnergy += get<0>(HitEnergyLayer[j]);
              usedLayer.emplace_back(j);
            }
        }

        temp_HitEnergyLayer.emplace_back(std::make_tuple(temp_HitEnergy, temp_HitLayer));
    }

    HitEnergyLayer.clear();
    HitEnergyLayer = temp_HitEnergyLayer;

    return;
}


void SiliconHits_SD::SignalstoHits_SD(std::vector<std::tuple<double, size_t> > HitEnergyLayerX,
                                      std::vector<std::tuple<double, size_t> > HitEnergyLayerY,
                                      std::vector<std::vector<double> >& HitEnergyPosXY)
{
  std::vector<double> countHitLayerX(HitEnergyLayerX.size(), 0.);
  std::vector<double> countHitLayerY(HitEnergyLayerY.size(), 0.);

  // Events with multiplicityX = multiplicityY = 1
  if((HitEnergyLayerX.size() > 0.1) && (HitEnergyLayerY.size() > 0.1))
    {
      for(size_t i = 0; i < HitEnergyLayerX.size(); ++i)
        {
          for(size_t j = 0; j < HitEnergyLayerY.size(); ++j)
            {
              double EnergyDiff = std::abs(get<0>(HitEnergyLayerX[i]) - get<0>(HitEnergyLayerY[j]));

              if(EnergyDiff > MaxEnergyDiffStrips)
                continue;

              if(ifveto && (std::find(inactwiresX.begin(), inactwiresX.end(), get<1>(HitEnergyLayerX[i])) != inactwiresX.end()))
                continue;

              if(ifveto && (std::find(inactwiresY.begin(), inactwiresY.end(), get<1>(HitEnergyLayerY[j])) != inactwiresY.end()))
                continue;

              if(IsValidCuadrant(get<1>(HitEnergyLayerX[i]), get<1>(HitEnergyLayerY[j])))
                {
                  double HitX       = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[i])%nStrips + 0.5) * widthStrip;
                  double HitY       = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[j])%nStrips + 0.5) * widthStrip;
                  double HitEnergy = get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerY[j]);

                  if(restrict_actlenght && ((HitX < -actlenghtX/2.) || (HitX > actlenghtX/2.) || (HitY < -actlenghtY/2.) || (HitY > actlenghtY/2.)))
                    continue;

                  if(restrict_gapcenter && (HitX > -gapcenter/2.) && (HitX < gapcenter/2.))
                    continue;

                  if(restrict_gapcenter && (HitY > -gapcenter/2.) && (HitY < gapcenter/2.))
                    continue;

                  std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                  HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                  countHitLayerX[i] += 1.;
                  countHitLayerY[j] += 1.;
                }
            }
        }
    }

#ifdef RECONS_HITS_MULTIPLICITY
    //If used, it has to be added:
    // - restrict_actlenght conditions
    // - inactwires conditions

  for(int i = HitEnergyLayerX.size() - 1; i >= 0; --i)
    {
      if(countHitLayerX[i] > 0.1)
        {
          HitEnergyLayerX.erase(HitEnergyLayerX.begin() + i);
          countHitLayerX.erase(countHitLayerX.begin() + i);
        }
    }

  for(int j = HitEnergyLayerY.size() - 1; j >= 0; --j)
    {
      if(countHitLayerY[j] > 0.1)
        {
          HitEnergyLayerY.erase(HitEnergyLayerY.begin() + j);
          countHitLayerY.erase(countHitLayerY.begin() + j);
        }
    }

  // Events with multiplicityX = 1 and multiplicityY = 2
  if((HitEnergyLayerX.size() > 1.1) && (HitEnergyLayerY.size() > 0.1))
    {
      for(size_t i = 0; i < HitEnergyLayerX.size() - 1; ++i)
        {
          for(size_t j = i + 1; j < HitEnergyLayerX.size(); ++j)
            {
              for(size_t k = 0; k < HitEnergyLayerY.size(); ++k)
                {
                  double EnergyDiff =
                      std::abs(get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j]) - get<0>(HitEnergyLayerY[k]));
                  if(EnergyDiff > MaxEnergyDiffStrips)
                    continue;

                  if((IsValidCuadrant(get<1>(HitEnergyLayerX[i]), get<1>(HitEnergyLayerY[k])))
                    && (IsValidCuadrant(get<1>(HitEnergyLayerX[j]), get<1>(HitEnergyLayerY[k]))))
                    {
                         if((size_type_abs(get<1>(HitEnergyLayerX[i]), get<1>(HitEnergyLayerX[j])) == 1) &&
                         ((get<0>(HitEnergyLayerX[i]) < MaxEnergyMultiplicity) ||
                          (get<0>(HitEnergyLayerX[j]) < MaxEnergyMultiplicity)))
                        {
                          double HitX =
                              -lenghtSi / 2. + ((get<0>(HitEnergyLayerX[i]) * (get<1>(HitEnergyLayerX[i])%nStrips) +
                                                 get<0>(HitEnergyLayerX[j]) * (get<1>(HitEnergyLayerX[j])%nStrips)) /
                                                    (get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j])) +
                                                0.5) *
                                                   widthStrip;
                          double HitY = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[k])%nStrips + 0.5) * widthStrip;
                          double HitEnergy =
                              get<0>(HitEnergyLayerY[k]) + get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j]);

                          std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                          countHitLayerX[i] += 0.5;
                          countHitLayerX[j] += 0.5;
                          countHitLayerY[k] += 1.;
                        }
                      else
                        {
                          double HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[i])%nStrips + 0.5) * widthStrip;
                          double HitY      = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[k])%nStrips + 0.5) * widthStrip;
                          double HitEnergy = get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerY[k]) / 2.;

                          std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                          HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[j])%nStrips + 0.5) * widthStrip;
                          HitEnergy = get<0>(HitEnergyLayerX[j]) + get<0>(HitEnergyLayerY[k]) / 2.;
                          std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                          countHitLayerX[i] += 1.;
                          countHitLayerX[j] += 1.;
                          countHitLayerY[k] += 2.;
                        }
                    }
                }
            }
        }

      for(int i = HitEnergyLayerX.size() - 1; i >= 0; --i)
        {
          if(countHitLayerX[i] > 0.1)
            {
              HitEnergyLayerX.erase(HitEnergyLayerX.begin() + i);
              countHitLayerX.erase(countHitLayerX.begin() + i);
            }
        }

      for(int j = HitEnergyLayerY.size() - 1; j >= 0; --j)
        {
          if(countHitLayerY[j] > 0.1)
            {
              HitEnergyLayerY.erase(HitEnergyLayerY.begin() + j);
              countHitLayerY.erase(countHitLayerY.begin() + j);
            }
        }
    }

  // Events with multiplicityX = 2 and multiplicityY = 1
  if((HitEnergyLayerX.size() > 0.1) && (HitEnergyLayerY.size() > 1.1))
    {
      for(size_t i = 0; i < HitEnergyLayerY.size() - 1; ++i)
        {
          for(size_t j = i + 1; j < HitEnergyLayerY.size(); ++j)
            {
              for(size_t k = 0; k < HitEnergyLayerX.size(); ++k)
                {
                  double EnergyDiff =
                      std::abs(get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j]) - get<0>(HitEnergyLayerX[k]));
                  if(EnergyDiff > MaxEnergyDiffStrips)
                    continue;

                  if((IsValidCuadrant(get<1>(HitEnergyLayerX[k]), get<1>(HitEnergyLayerY[i]))) && 
                      (IsValidCuadrant(get<1>(HitEnergyLayerX[k]), get<1>(HitEnergyLayerY[j]))))
                    {
                      if((size_type_abs(get<1>(HitEnergyLayerY[i]), get<1>(HitEnergyLayerY[j])) == 1) &&
                         ((get<0>(HitEnergyLayerY[i]) < MaxEnergyMultiplicity) ||
                          (get<0>(HitEnergyLayerY[j]) < MaxEnergyMultiplicity)))
                        {
                          double HitX = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[k])%nStrips + 0.5) * widthStrip;
                          double HitY = +lenghtSi / 2. - ((get<0>(HitEnergyLayerY[i]) * (get<1>(HitEnergyLayerY[i])%nStrips) +
                                                          get<0>(HitEnergyLayerY[j]) * (get<1>(HitEnergyLayerY[j])%nStrips)) /
                                                             (get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j])) +
                                                         0.5) *
                                                            widthStrip;
                          double HitEnergy =
                              get<0>(HitEnergyLayerX[k]) + get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j]);
                          std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                          countHitLayerY[i] += 0.5;
                          countHitLayerY[j] += 0.5;
                          countHitLayerX[k] += 1.;
                        }
                      else
                        {
                          double HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[k])%nStrips + 0.5) * widthStrip;
                          double HitY      = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[i])%nStrips + 0.5) * widthStrip;
                          double HitEnergy = get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerX[k]) / 2.;
                          std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                          HitY      = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[j])%nStrips + 0.5) * widthStrip;
                          HitEnergy = get<0>(HitEnergyLayerY[j]) + get<0>(HitEnergyLayerX[k]) / 2.;
                          std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                          countHitLayerY[i] += 1.;
                          countHitLayerY[j] += 1.;
                          countHitLayerX[k] += 2.;
                        }
                    }
                }
            }
        }

      for(int i = HitEnergyLayerX.size() - 1; i >= 0; --i)
        {
          if(countHitLayerX[i] > 0.1)
            {
              HitEnergyLayerX.erase(HitEnergyLayerX.begin() + i);
              countHitLayerX.erase(countHitLayerX.begin() + i);
            }
        }

      for(int j = HitEnergyLayerY.size() - 1; j >= 0; --j)
        {
          if(countHitLayerY[j] > 0.1)
            {
              HitEnergyLayerY.erase(HitEnergyLayerY.begin() + j);
              countHitLayerY.erase(countHitLayerY.begin() + j);
            }
        }
    }

  // Events with multiplicityX = 1 and multiplicityY = 3
  if((HitEnergyLayerX.size() > 2.1) && (HitEnergyLayerY.size() > 0.1))
    {
      for(size_t i = 0; i < HitEnergyLayerX.size() - 2; ++i)
        {
          for(size_t j = i + 1; j < HitEnergyLayerX.size() - 1; ++j)
            {
              for(size_t k = j + 1; k < HitEnergyLayerX.size(); ++k)
                {
                  for(size_t l = 0; l < HitEnergyLayerY.size(); ++l)
                    {
                      double EnergyDiff = std::abs(get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j]) +
                                              get<0>(HitEnergyLayerX[k]) - get<0>(HitEnergyLayerY[l]));
                      if(EnergyDiff > MaxEnergyDiffStrips)
                        continue;

                      if((IsValidCuadrant(get<1>(HitEnergyLayerX[i]), get<1>(HitEnergyLayerY[l])))
                        && (IsValidCuadrant(get<1>(HitEnergyLayerX[j]), get<1>(HitEnergyLayerY[l])))
                        && (IsValidCuadrant(get<1>(HitEnergyLayerX[k]), get<1>(HitEnergyLayerY[l]))))
                        {
                          if((size_type_abs(get<1>(HitEnergyLayerX[i]), get<1>(HitEnergyLayerX[j])) == 1) &&
                             ((get<0>(HitEnergyLayerX[i]) < MaxEnergyMultiplicity) ||
                              (get<0>(HitEnergyLayerX[j]) < MaxEnergyMultiplicity)))
                            {
                              double HitX =
                                  -lenghtSi / 2. + ((get<0>(HitEnergyLayerX[i]) * (get<1>(HitEnergyLayerX[i])%nStrips) +
                                                     get<0>(HitEnergyLayerX[j]) * (get<1>(HitEnergyLayerX[j]))%nStrips) /
                                                        (get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j])) +
                                                    0.5) *
                                                       widthStrip;
                              double HitY      = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[l])%nStrips + 0.5) * widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j]) +
                                                 get<0>(HitEnergyLayerY[l]) / 2.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[k])%nStrips + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerX[k]) + get<0>(HitEnergyLayerY[l]) / 2.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              countHitLayerX[i] += 0.5;
                              countHitLayerX[j] += 0.5;
                              countHitLayerX[k] += 1.;
                              countHitLayerY[l] += 2.;
                            }
                          else if((size_type_abs(get<1>(HitEnergyLayerX[j]), get<1>(HitEnergyLayerX[k])) == 1) &&
                                  ((get<0>(HitEnergyLayerX[j]) < MaxEnergyMultiplicity) ||
                                   (get<0>(HitEnergyLayerX[k]) < MaxEnergyMultiplicity)))
                            {
                              double HitX =
                                  -lenghtSi / 2. + ((get<0>(HitEnergyLayerX[j]) * (get<1>(HitEnergyLayerX[j])%nStrips) +
                                                     get<0>(HitEnergyLayerX[k]) * (get<1>(HitEnergyLayerX[k])%nStrips)) /
                                                        (get<0>(HitEnergyLayerX[j]) + get<0>(HitEnergyLayerX[k])) +
                                                    0.5) *
                                                       widthStrip;
                              double HitY      = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[l])%nStrips + 0.5) * widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerX[j]) + get<0>(HitEnergyLayerX[k]) +
                                                 get<0>(HitEnergyLayerY[l]) / 2.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[i])%nStrips + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerY[l]) / 2.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              countHitLayerX[i] += 1.;
                              countHitLayerX[j] += 0.5;
                              countHitLayerX[k] += 0.5;
                              countHitLayerY[l] += 2.;
                            }
                          else
                            {
                              double HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[i])%nStrips + 0.5) * widthStrip;
                              double HitY      = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[l])%nStrips + 0.5) * widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerY[l]) / 3.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[j])%nStrips + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerX[j]) + get<0>(HitEnergyLayerY[l]) / 3.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[k])%nStrips + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerX[k]) + get<0>(HitEnergyLayerY[l]) / 3.;
                              std::vector<double> tempEnergyPosXY3{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY3);

                              countHitLayerX[i] += 1.;
                              countHitLayerX[j] += 1.;
                              countHitLayerX[k] += 1.;
                              countHitLayerY[l] += 3.;
                            }
                        }
                    }
                }
            }
        }

      for(int i = HitEnergyLayerX.size() - 1; i >= 0; --i)
        {
          if(countHitLayerX[i] > 0.1)
            {
              HitEnergyLayerX.erase(HitEnergyLayerX.begin() + i);
              countHitLayerX.erase(countHitLayerX.begin() + i);
            }
        }

      for(int j = HitEnergyLayerY.size() - 1; j >= 0; --j)
        {
          if(countHitLayerY[j] > 0.1)
            {
              HitEnergyLayerY.erase(HitEnergyLayerY.begin() + j);
              countHitLayerY.erase(countHitLayerY.begin() + j);
            }
        }
    }

  // Events with multiplicityX = 3 and multiplicityY = 1
  if((HitEnergyLayerX.size() > 0.1) && (HitEnergyLayerY.size() > 2.1))
    {
      for(size_t i = 0; i < HitEnergyLayerY.size() - 2; ++i)
        {
          for(size_t j = i + 1; j < HitEnergyLayerY.size() - 1; ++j)
            {
              for(size_t k = j + 1; k < HitEnergyLayerY.size(); ++k)
                {
                  for(size_t l = 0; l < HitEnergyLayerX.size(); ++l)
                    {
                      double EnergyDiff = std::abs(get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j]) +
                                              get<0>(HitEnergyLayerY[k]) - get<0>(HitEnergyLayerX[l]));
                      if(EnergyDiff > MaxEnergyDiffStrips)
                        continue;

                      if((IsValidCuadrant(get<1>(HitEnergyLayerX[l]), get<1>(HitEnergyLayerY[i]))) && 
                      (IsValidCuadrant(get<1>(HitEnergyLayerX[l]), get<1>(HitEnergyLayerY[j]))) &&
                      (IsValidCuadrant(get<1>(HitEnergyLayerX[l]), get<1>(HitEnergyLayerY[k]))))
                        {
                          if((size_type_abs(get<1>(HitEnergyLayerY[i]), get<1>(HitEnergyLayerY[j])) == 1) &&
                             ((get<0>(HitEnergyLayerY[i]) < MaxEnergyMultiplicity) ||
                              (get<0>(HitEnergyLayerY[j]) < MaxEnergyMultiplicity)))
                            {
                              double HitX = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[l])%nStrips + 0.5) * widthStrip;
                              double HitY =
                                  +lenghtSi / 2. - ((get<0>(HitEnergyLayerY[i]) * (get<1>(HitEnergyLayerY[i])%nStrips) +
                                                    get<0>(HitEnergyLayerY[j]) * (get<1>(HitEnergyLayerY[j])%nStrips)) /
                                                       (get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j])) +
                                                   0.5) *
                                                      widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j]) +
                                                 get<0>(HitEnergyLayerX[l]) / 2.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[k])%nStrips + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerY[k]) + get<0>(HitEnergyLayerX[l]) / 2.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              countHitLayerY[i] += 0.5;
                              countHitLayerY[j] += 0.5;
                              countHitLayerY[k] += 1.;
                              countHitLayerX[l] += 2.;
                            }
                          else if((size_type_abs(get<1>(HitEnergyLayerY[j]), get<1>(HitEnergyLayerY[k])) == 1) &&
                                  ((get<0>(HitEnergyLayerY[j]) < MaxEnergyMultiplicity) ||
                                   (get<0>(HitEnergyLayerY[k]) < MaxEnergyMultiplicity)))
                            {
                              double HitX = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[l])%nStrips + 0.5) * widthStrip;
                              double HitY =
                                  +lenghtSi / 2. - ((get<0>(HitEnergyLayerY[j]) * (get<1>(HitEnergyLayerY[j])%nStrips +
                                                    get<0>(HitEnergyLayerY[k]) * (get<1>(HitEnergyLayerY[k])%nStrips)) /
                                                       (get<0>(HitEnergyLayerY[j]) + get<0>(HitEnergyLayerY[k])) +
                                                   0.5) *
                                                      widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerY[j]) + get<0>(HitEnergyLayerY[k]) +
                                                 get<0>(HitEnergyLayerX[l]) / 2.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[i])%nStrips + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerX[l]) / 2.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              countHitLayerY[i] += 1.;
                              countHitLayerY[j] += 0.5;
                              countHitLayerY[k] += 0.5;
                              countHitLayerX[l] += 2.;
                            }
                          else
                            {
                              double HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[l])%nStrips + 0.5) * widthStrip;
                              double HitY      = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[i])%nStrips + 0.5) * widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerX[l]) / 3.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[j])%nStrips + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerY[j]) + get<0>(HitEnergyLayerX[l]) / 3.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[k])%nStrips + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerY[k]) + get<0>(HitEnergyLayerX[l]) / 3.;
                              std::vector<double> tempEnergyPosXY3{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY3);

                              countHitLayerY[i] += 1.;
                              countHitLayerY[j] += 1.;
                              countHitLayerY[k] += 1.;
                              countHitLayerX[l] += 3.;
                            }
                        }
                    }
                }
            }
        }

      for(int i = HitEnergyLayerX.size() - 1; i >= 0; --i)
        {
          if(countHitLayerX[i] > 0.1)
            {
              HitEnergyLayerX.erase(HitEnergyLayerX.begin() + i);
              countHitLayerX.erase(countHitLayerX.begin() + i);
            }
        }

      for(int j = HitEnergyLayerY.size() - 1; j >= 0; --j)
        {
          if(countHitLayerY[j] > 0.1)
            {
              HitEnergyLayerY.erase(HitEnergyLayerY.begin() + j);
              countHitLayerY.erase(countHitLayerY.begin() + j);
            }
        }
    }

  if(HitEnergyPosXY.size() >= 2)
    {
      for(size_t i = 0; i < HitEnergyPosXY.size() - 1; ++i)
        {
          for(size_t j = i + 1; j < HitEnergyPosXY.size(); ++j)
            {
              if((HitEnergyPosXY[i][0] == HitEnergyPosXY[j][0]) && (HitEnergyPosXY[i][1] == HitEnergyPosXY[j][1]) &&
                 (HitEnergyPosXY[i][2] == HitEnergyPosXY[j][2]))
                {
                  HitEnergyPosXY.erase(HitEnergyPosXY.begin() + j);
                  --j;
                }
            }
        }
    }

#endif
}
*/

void SiliconHits_SDpad::TranslationZ_Si_System(double Target_PosZ)
{
  double Dist_to_TransZ = Target_PosZ - 25.; // Reference target PosZ = 25.
  Z_plane += Dist_to_TransZ;
}

bool SiliconHits_SDpad::IsSamePad(size_t& layerX_id, size_t& layerY_id)
{
  size_t iPad_X = layerX_id / nStrips;
  size_t iPad_Y = layerY_id / nStrips;

  if(iPad_X == iPad_Y)
    return true;

  else
    {
      nRejectPad += 1;
      return false;
    }
}

void SiliconHits_SDpad::originPad(size_t& layerX_id, double& originPadX, double& originPadY)
{
  size_t temp_Pad = layerX_id / nStrips / 2;
  size_t temp_Row = temp_Pad / PadDistribution[0];
  size_t temp_Col = temp_Pad % PadDistribution[0];

  originPadX =  lenghtSi * (temp_Col - PadDistribution[0]/2.) + gapcenter * (temp_Col - PadDistribution[0]/2. + 0.5);
  originPadY = -lenghtSi * (temp_Row - PadDistribution[1]/2.) - gapcenter * (temp_Row - PadDistribution[1]/2. + 0.5);

 #ifdef DEBUG_PRIMVTX
    std::cout << "Si_station: " << Si_station << "\n";
    std::cout << "temp_Pad: " << temp_Pad << "\n";
    std::cout << "temp_Col: " << temp_Col << "\n";
    std::cout << "temp_Row: " << temp_Row << "\n";
    std::cout << "PadDistribution: " << PadDistribution[0] << ", " << PadDistribution[1] << "\n";
    std::cout << "origin_X: " << originPadX << "\n";
    std::cout << "origin_Y: " << originPadY << "\n\n";
  #endif
}

void SiliconHits_SDpad::CombineStrips(std::vector<std::tuple<double, size_t> >& HitEnergyLayer)
{
  if(combineStrips == 1)
    return;

  if((HitEnergyLayer.size() == 0) || (HitEnergyLayer.size() == 1))
    return;

  std::vector<std::tuple<double, size_t> > temp_HitEnergyLayer;

  std::vector<size_t> usedLayer {};
  std::vector<size_t>::iterator it_usedLayer;

  for(size_t i = 0; i < HitEnergyLayer.size()-1; ++i)
    {
      it_usedLayer = find(usedLayer.begin(), usedLayer.end(), get<1>(HitEnergyLayer[i]));
      if(it_usedLayer != usedLayer.end())
        continue;

      double temp_HitEnergy = get<0>(HitEnergyLayer[i]);
      size_t temp_HitLayer = (get<1>(HitEnergyLayer[i])%(nStrips*combineStrips))/combineStrips;
      
      for(size_t k = 1; k < nPads; ++k)
        if((get<1>(HitEnergyLayer[i]) >= k*2*nStrips*combineStrips) && (get<1>(HitEnergyLayer[i]) < (k+0.5)*2*nStrips*combineStrips))
          temp_HitLayer += k*2*nStrips;

      for(size_t j = i+1; j < HitEnergyLayer.size(); ++j)
        {
          it_usedLayer = find(usedLayer.begin(), usedLayer.end(), get<1>(HitEnergyLayer[j]));
          if(it_usedLayer != usedLayer.end())
            continue;

          size_t temp_diflayers = size_type_abs(get<1>(HitEnergyLayer[i]), get<1>(HitEnergyLayer[j]));
          if(temp_diflayers >= combineStrips)
            continue;

          if(((get<1>(HitEnergyLayer[i])%(nStrips*combineStrips))/combineStrips)
              != ((get<1>(HitEnergyLayer[j])%(nStrips*combineStrips))/combineStrips))
                continue;

          temp_HitEnergy += get<0>(HitEnergyLayer[j]);
          usedLayer.emplace_back(get<1>(HitEnergyLayer[j]));
        }

        temp_HitEnergyLayer.emplace_back(std::make_tuple(temp_HitEnergy, temp_HitLayer));
    }

    //Check last entry: vector.size()-1
    it_usedLayer = find(usedLayer.begin(), usedLayer.end(), get<1>(HitEnergyLayer[HitEnergyLayer.size()-1]));
    if(it_usedLayer == usedLayer.end())
    {
      double temp_HitEnergy = get<0>(HitEnergyLayer[HitEnergyLayer.size()-1]);
      size_t temp_HitLayer = (get<1>(HitEnergyLayer[HitEnergyLayer.size()-1])%(nStrips*combineStrips))/combineStrips;

      for(size_t k = 1; k < nPads; ++k)
        if((get<1>(HitEnergyLayer[HitEnergyLayer.size()-1]) >= k*2*nStrips*combineStrips)
          && (get<1>(HitEnergyLayer[HitEnergyLayer.size()-1]) < (k+0.5)*2*nStrips*combineStrips))
            temp_HitLayer += k*2*nStrips;

      temp_HitEnergyLayer.emplace_back(std::make_tuple(temp_HitEnergy, temp_HitLayer));
    }

    HitEnergyLayer.clear();
    HitEnergyLayer = temp_HitEnergyLayer;

    return;
}


void SiliconHits_SDpad::SignalstoHits_SDpad(std::vector<std::tuple<double, size_t> > HitEnergyLayerX,
                                      std::vector<std::tuple<double, size_t> > HitEnergyLayerY,
                                      std::vector<std::vector<double> >& HitEnergyPosXY)
{
  std::vector<double> countHitLayerX(HitEnergyLayerX.size(), 0.);
  std::vector<double> countHitLayerY(HitEnergyLayerY.size(), 0.);

  double originPadX;
  double originPadY;

  // Events with multiplicityX = multiplicityY = 1
  if((HitEnergyLayerX.size() > 0.1) && (HitEnergyLayerY.size() > 0.1))
    {
      for(size_t i = 0; i < HitEnergyLayerX.size(); ++i)
        {
          if(ifveto && (std::find(inactwiresX.begin(), inactwiresX.end(), get<1>(HitEnergyLayerX[i])) != inactwiresX.end()))
            continue;

          for(size_t j = 0; j < HitEnergyLayerY.size(); ++j)
            {
              if(ifveto && (std::find(inactwiresY.begin(), inactwiresY.end(), get<1>(HitEnergyLayerY[j])) != inactwiresY.end()))
                continue;

              double EnergyDiff = std::abs(get<0>(HitEnergyLayerX[i]) - get<0>(HitEnergyLayerY[j]));
              if(EnergyDiff > MaxEnergyDiffStrips)
                continue;

              if(IsSamePad(get<1>(HitEnergyLayerX[i]), get<1>(HitEnergyLayerY[j])))
                {
                  originPad(get<1>(HitEnergyLayerX[i]), originPadX, originPadY);

                  double HitX_prime = originPadX + (get<1>(HitEnergyLayerX[i])%nStrips + 0.5) * widthStrip;
                  double HitY_prime = originPadY - (get<1>(HitEnergyLayerY[j])%nStrips + 0.5) * widthStrip;

                  if(ifreflectionX)
                    HitX_prime *= -1.;

                  if(ifreflectionY)
                    HitY_prime *= -1.;

                  if(ifphirot)
                  {
                    double HitX_rot = HitX_prime * std::cos(phirot * M_PI / 180.) - HitY_prime * std::sin(phirot * M_PI / 180.);
                    double HitY_rot = HitX_prime * std::sin(phirot * M_PI / 180.) + HitY_prime * std::cos(phirot * M_PI / 180.);

                    HitX_prime = HitX_rot;
                    HitY_prime = HitY_rot;
                  }

                  //if(Si_station == 1)
                      //std::cout << "ReconsPosX : " << HitX_prime << "\tReconsPosY : " << HitY_prime << std::endl;

                  double HitX = HitX_prime;
                  double HitY = HitY_prime;
                  double HitEnergy = get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerY[j]);

#ifdef DEBUG_PRIMVTX
                  std::cout << "HitX: " << HitX << "\n";
                  std::cout << "LayerX: " << get<1>(HitEnergyLayerX[i])%nStrips << "\n";
                  std::cout << "HitY: " << HitY << "\n";
                  std::cout << "LayerY: " << get<1>(HitEnergyLayerY[j])%nStrips << "\n\n";
#endif

                  if(restrict_actlenght && ((HitX < -actlenghtX/2.) || (HitX > actlenghtX/2.) || (HitY < -actlenghtY/2.) || (HitY > actlenghtY/2.)))
                    continue;

                  if(restrict_gapcenter && (HitX > -gapcenter/2.) && (HitX < gapcenter/2.))
                    continue;

                  if(restrict_gapcenter && (HitY > -gapcenter/2.) && (HitY < gapcenter/2.))
                    continue;

                  std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                  HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                  countHitLayerX[i] += 1.;
                  countHitLayerY[j] += 1.;
                }
            }
        }
    }

#ifdef RECONS_HITS_MULTIPLICITY
    //If used, it has to be added:
    // - restrict_actlenght conditions
    // - inactwires conditions
    // - changes from hamamatsu to pad distribution

  for(int i = HitEnergyLayerX.size() - 1; i >= 0; --i)
    {
      if(countHitLayerX[i] > 0.1)
        {
          HitEnergyLayerX.erase(HitEnergyLayerX.begin() + i);
          countHitLayerX.erase(countHitLayerX.begin() + i);
        }
    }

  for(int j = HitEnergyLayerY.size() - 1; j >= 0; --j)
    {
      if(countHitLayerY[j] > 0.1)
        {
          HitEnergyLayerY.erase(HitEnergyLayerY.begin() + j);
          countHitLayerY.erase(countHitLayerY.begin() + j);
        }
    }

  // Events with multiplicityX = 1 and multiplicityY = 2
  if((HitEnergyLayerX.size() > 1.1) && (HitEnergyLayerY.size() > 0.1))
    {
      for(size_t i = 0; i < HitEnergyLayerX.size() - 1; ++i)
        {
          for(size_t j = i + 1; j < HitEnergyLayerX.size(); ++j)
            {
              for(size_t k = 0; k < HitEnergyLayerY.size(); ++k)
                {
                  double EnergyDiff =
                      std::abs(get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j]) - get<0>(HitEnergyLayerY[k]));
                  if(EnergyDiff > MaxEnergyDiffStrips)
                    continue;

                  if((IsSamePad(get<1>(HitEnergyLayerX[i]), get<1>(HitEnergyLayerY[k])))
                    && (IsSamePad(get<1>(HitEnergyLayerX[j]), get<1>(HitEnergyLayerY[k]))))
                    {
                         if((size_type_abs(get<1>(HitEnergyLayerX[i]), get<1>(HitEnergyLayerX[j])) == 1) &&
                         ((get<0>(HitEnergyLayerX[i]) < MaxEnergyMultiplicity) ||
                          (get<0>(HitEnergyLayerX[j]) < MaxEnergyMultiplicity)))
                        {
                          double HitX =
                              -lenghtSi / 2. + ((get<0>(HitEnergyLayerX[i]) * (get<1>(HitEnergyLayerX[i])%nStrips) +
                                                 get<0>(HitEnergyLayerX[j]) * (get<1>(HitEnergyLayerX[j])%nStrips)) /
                                                    (get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j])) +
                                                0.5) *
                                                   widthStrip;
                          double HitY = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[k])%nStrips + 0.5) * widthStrip;
                          double HitEnergy =
                              get<0>(HitEnergyLayerY[k]) + get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j]);

                          std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                          countHitLayerX[i] += 0.5;
                          countHitLayerX[j] += 0.5;
                          countHitLayerY[k] += 1.;
                        }
                      else
                        {
                          double HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[i])%nStrips + 0.5) * widthStrip;
                          double HitY      = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[k])%nStrips + 0.5) * widthStrip;
                          double HitEnergy = get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerY[k]) / 2.;

                          std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                          HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[j])%nStrips + 0.5) * widthStrip;
                          HitEnergy = get<0>(HitEnergyLayerX[j]) + get<0>(HitEnergyLayerY[k]) / 2.;
                          std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                          countHitLayerX[i] += 1.;
                          countHitLayerX[j] += 1.;
                          countHitLayerY[k] += 2.;
                        }
                    }
                }
            }
        }

      for(int i = HitEnergyLayerX.size() - 1; i >= 0; --i)
        {
          if(countHitLayerX[i] > 0.1)
            {
              HitEnergyLayerX.erase(HitEnergyLayerX.begin() + i);
              countHitLayerX.erase(countHitLayerX.begin() + i);
            }
        }

      for(int j = HitEnergyLayerY.size() - 1; j >= 0; --j)
        {
          if(countHitLayerY[j] > 0.1)
            {
              HitEnergyLayerY.erase(HitEnergyLayerY.begin() + j);
              countHitLayerY.erase(countHitLayerY.begin() + j);
            }
        }
    }

  // Events with multiplicityX = 2 and multiplicityY = 1
  if((HitEnergyLayerX.size() > 0.1) && (HitEnergyLayerY.size() > 1.1))
    {
      for(size_t i = 0; i < HitEnergyLayerY.size() - 1; ++i)
        {
          for(size_t j = i + 1; j < HitEnergyLayerY.size(); ++j)
            {
              for(size_t k = 0; k < HitEnergyLayerX.size(); ++k)
                {
                  double EnergyDiff =
                      std::abs(get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j]) - get<0>(HitEnergyLayerX[k]));
                  if(EnergyDiff > MaxEnergyDiffStrips)
                    continue;

                  if((IsSamePad(get<1>(HitEnergyLayerX[k]), get<1>(HitEnergyLayerY[i]))) && 
                      (IsSamePad(get<1>(HitEnergyLayerX[k]), get<1>(HitEnergyLayerY[j]))))
                    {
                      if((size_type_abs(get<1>(HitEnergyLayerY[i]), get<1>(HitEnergyLayerY[j])) == 1) &&
                         ((get<0>(HitEnergyLayerY[i]) < MaxEnergyMultiplicity) ||
                          (get<0>(HitEnergyLayerY[j]) < MaxEnergyMultiplicity)))
                        {
                          double HitX = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[k])%nStrips + 0.5) * widthStrip;
                          double HitY = +lenghtSi / 2. - ((get<0>(HitEnergyLayerY[i]) * (get<1>(HitEnergyLayerY[i])%nStrips) +
                                                          get<0>(HitEnergyLayerY[j]) * (get<1>(HitEnergyLayerY[j])%nStrips)) /
                                                             (get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j])) +
                                                         0.5) *
                                                            widthStrip;
                          double HitEnergy =
                              get<0>(HitEnergyLayerX[k]) + get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j]);
                          std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                          countHitLayerY[i] += 0.5;
                          countHitLayerY[j] += 0.5;
                          countHitLayerX[k] += 1.;
                        }
                      else
                        {
                          double HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[k])%nStrips + 0.5) * widthStrip;
                          double HitY      = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[i])%nStrips + 0.5) * widthStrip;
                          double HitEnergy = get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerX[k]) / 2.;
                          std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                          HitY      = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[j])%nStrips + 0.5) * widthStrip;
                          HitEnergy = get<0>(HitEnergyLayerY[j]) + get<0>(HitEnergyLayerX[k]) / 2.;
                          std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                          HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                          countHitLayerY[i] += 1.;
                          countHitLayerY[j] += 1.;
                          countHitLayerX[k] += 2.;
                        }
                    }
                }
            }
        }

      for(int i = HitEnergyLayerX.size() - 1; i >= 0; --i)
        {
          if(countHitLayerX[i] > 0.1)
            {
              HitEnergyLayerX.erase(HitEnergyLayerX.begin() + i);
              countHitLayerX.erase(countHitLayerX.begin() + i);
            }
        }

      for(int j = HitEnergyLayerY.size() - 1; j >= 0; --j)
        {
          if(countHitLayerY[j] > 0.1)
            {
              HitEnergyLayerY.erase(HitEnergyLayerY.begin() + j);
              countHitLayerY.erase(countHitLayerY.begin() + j);
            }
        }
    }

  // Events with multiplicityX = 1 and multiplicityY = 3
  if((HitEnergyLayerX.size() > 2.1) && (HitEnergyLayerY.size() > 0.1))
    {
      for(size_t i = 0; i < HitEnergyLayerX.size() - 2; ++i)
        {
          for(size_t j = i + 1; j < HitEnergyLayerX.size() - 1; ++j)
            {
              for(size_t k = j + 1; k < HitEnergyLayerX.size(); ++k)
                {
                  for(size_t l = 0; l < HitEnergyLayerY.size(); ++l)
                    {
                      double EnergyDiff = std::abs(get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j]) +
                                              get<0>(HitEnergyLayerX[k]) - get<0>(HitEnergyLayerY[l]));
                      if(EnergyDiff > MaxEnergyDiffStrips)
                        continue;

                      if((IsSamePad(get<1>(HitEnergyLayerX[i]), get<1>(HitEnergyLayerY[l])))
                        && (IsSamePad(get<1>(HitEnergyLayerX[j]), get<1>(HitEnergyLayerY[l])))
                        && (IsSamePad(get<1>(HitEnergyLayerX[k]), get<1>(HitEnergyLayerY[l]))))
                        {
                          if((size_type_abs(get<1>(HitEnergyLayerX[i]), get<1>(HitEnergyLayerX[j])) == 1) &&
                             ((get<0>(HitEnergyLayerX[i]) < MaxEnergyMultiplicity) ||
                              (get<0>(HitEnergyLayerX[j]) < MaxEnergyMultiplicity)))
                            {
                              double HitX =
                                  -lenghtSi / 2. + ((get<0>(HitEnergyLayerX[i]) * (get<1>(HitEnergyLayerX[i])%nStrips) +
                                                     get<0>(HitEnergyLayerX[j]) * (get<1>(HitEnergyLayerX[j]))%nStrips) /
                                                        (get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j])) +
                                                    0.5) *
                                                       widthStrip;
                              double HitY      = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[l])%nStrips + 0.5) * widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerX[j]) +
                                                 get<0>(HitEnergyLayerY[l]) / 2.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[k])%nStrips + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerX[k]) + get<0>(HitEnergyLayerY[l]) / 2.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              countHitLayerX[i] += 0.5;
                              countHitLayerX[j] += 0.5;
                              countHitLayerX[k] += 1.;
                              countHitLayerY[l] += 2.;
                            }
                          else if((size_type_abs(get<1>(HitEnergyLayerX[j]), get<1>(HitEnergyLayerX[k])) == 1) &&
                                  ((get<0>(HitEnergyLayerX[j]) < MaxEnergyMultiplicity) ||
                                   (get<0>(HitEnergyLayerX[k]) < MaxEnergyMultiplicity)))
                            {
                              double HitX =
                                  -lenghtSi / 2. + ((get<0>(HitEnergyLayerX[j]) * (get<1>(HitEnergyLayerX[j])%nStrips) +
                                                     get<0>(HitEnergyLayerX[k]) * (get<1>(HitEnergyLayerX[k])%nStrips)) /
                                                        (get<0>(HitEnergyLayerX[j]) + get<0>(HitEnergyLayerX[k])) +
                                                    0.5) *
                                                       widthStrip;
                              double HitY      = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[l])%nStrips + 0.5) * widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerX[j]) + get<0>(HitEnergyLayerX[k]) +
                                                 get<0>(HitEnergyLayerY[l]) / 2.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[i])%nStrips + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerY[l]) / 2.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              countHitLayerX[i] += 1.;
                              countHitLayerX[j] += 0.5;
                              countHitLayerX[k] += 0.5;
                              countHitLayerY[l] += 2.;
                            }
                          else
                            {
                              double HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[i])%nStrips + 0.5) * widthStrip;
                              double HitY      = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[l])%nStrips + 0.5) * widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerX[i]) + get<0>(HitEnergyLayerY[l]) / 3.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[j])%nStrips + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerX[j]) + get<0>(HitEnergyLayerY[l]) / 3.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[k])%nStrips + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerX[k]) + get<0>(HitEnergyLayerY[l]) / 3.;
                              std::vector<double> tempEnergyPosXY3{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY3);

                              countHitLayerX[i] += 1.;
                              countHitLayerX[j] += 1.;
                              countHitLayerX[k] += 1.;
                              countHitLayerY[l] += 3.;
                            }
                        }
                    }
                }
            }
        }

      for(int i = HitEnergyLayerX.size() - 1; i >= 0; --i)
        {
          if(countHitLayerX[i] > 0.1)
            {
              HitEnergyLayerX.erase(HitEnergyLayerX.begin() + i);
              countHitLayerX.erase(countHitLayerX.begin() + i);
            }
        }

      for(int j = HitEnergyLayerY.size() - 1; j >= 0; --j)
        {
          if(countHitLayerY[j] > 0.1)
            {
              HitEnergyLayerY.erase(HitEnergyLayerY.begin() + j);
              countHitLayerY.erase(countHitLayerY.begin() + j);
            }
        }
    }

  // Events with multiplicityX = 3 and multiplicityY = 1
  if((HitEnergyLayerX.size() > 0.1) && (HitEnergyLayerY.size() > 2.1))
    {
      for(size_t i = 0; i < HitEnergyLayerY.size() - 2; ++i)
        {
          for(size_t j = i + 1; j < HitEnergyLayerY.size() - 1; ++j)
            {
              for(size_t k = j + 1; k < HitEnergyLayerY.size(); ++k)
                {
                  for(size_t l = 0; l < HitEnergyLayerX.size(); ++l)
                    {
                      double EnergyDiff = std::abs(get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j]) +
                                              get<0>(HitEnergyLayerY[k]) - get<0>(HitEnergyLayerX[l]));
                      if(EnergyDiff > MaxEnergyDiffStrips)
                        continue;

                      if((IsSamePad(get<1>(HitEnergyLayerX[l]), get<1>(HitEnergyLayerY[i]))) && 
                      (IsSamePad(get<1>(HitEnergyLayerX[l]), get<1>(HitEnergyLayerY[j]))) &&
                      (IsSamePad(get<1>(HitEnergyLayerX[l]), get<1>(HitEnergyLayerY[k]))))
                        {
                          if((size_type_abs(get<1>(HitEnergyLayerY[i]), get<1>(HitEnergyLayerY[j])) == 1) &&
                             ((get<0>(HitEnergyLayerY[i]) < MaxEnergyMultiplicity) ||
                              (get<0>(HitEnergyLayerY[j]) < MaxEnergyMultiplicity)))
                            {
                              double HitX = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[l])%nStrips + 0.5) * widthStrip;
                              double HitY =
                                  +lenghtSi / 2. - ((get<0>(HitEnergyLayerY[i]) * (get<1>(HitEnergyLayerY[i])%nStrips) +
                                                    get<0>(HitEnergyLayerY[j]) * (get<1>(HitEnergyLayerY[j])%nStrips)) /
                                                       (get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j])) +
                                                   0.5) *
                                                      widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerY[j]) +
                                                 get<0>(HitEnergyLayerX[l]) / 2.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[k])%nStrips + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerY[k]) + get<0>(HitEnergyLayerX[l]) / 2.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              countHitLayerY[i] += 0.5;
                              countHitLayerY[j] += 0.5;
                              countHitLayerY[k] += 1.;
                              countHitLayerX[l] += 2.;
                            }
                          else if((size_type_abs(get<1>(HitEnergyLayerY[j]), get<1>(HitEnergyLayerY[k])) == 1) &&
                                  ((get<0>(HitEnergyLayerY[j]) < MaxEnergyMultiplicity) ||
                                   (get<0>(HitEnergyLayerY[k]) < MaxEnergyMultiplicity)))
                            {
                              double HitX = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[l])%nStrips + 0.5) * widthStrip;
                              double HitY =
                                  +lenghtSi / 2. - ((get<0>(HitEnergyLayerY[j]) * (get<1>(HitEnergyLayerY[j])%nStrips +
                                                    get<0>(HitEnergyLayerY[k]) * (get<1>(HitEnergyLayerY[k])%nStrips)) /
                                                       (get<0>(HitEnergyLayerY[j]) + get<0>(HitEnergyLayerY[k])) +
                                                   0.5) *
                                                      widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerY[j]) + get<0>(HitEnergyLayerY[k]) +
                                                 get<0>(HitEnergyLayerX[l]) / 2.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[i])%nStrips + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerX[l]) / 2.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              countHitLayerY[i] += 1.;
                              countHitLayerY[j] += 0.5;
                              countHitLayerY[k] += 0.5;
                              countHitLayerX[l] += 2.;
                            }
                          else
                            {
                              double HitX      = -lenghtSi / 2. + (get<1>(HitEnergyLayerX[l])%nStrips + 0.5) * widthStrip;
                              double HitY      = +lenghtSi / 2. - (get<1>(HitEnergyLayerY[i])%nStrips + 0.5) * widthStrip;
                              double HitEnergy = get<0>(HitEnergyLayerY[i]) + get<0>(HitEnergyLayerX[l]) / 3.;
                              std::vector<double> tempEnergyPosXY{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY);

                              HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[j])%nStrips + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerY[j]) + get<0>(HitEnergyLayerX[l]) / 3.;
                              std::vector<double> tempEnergyPosXY2{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY2);

                              HitY      = -lenghtSi / 2. + (get<1>(HitEnergyLayerY[k])%nStrips + 0.5) * widthStrip;
                              HitEnergy = get<0>(HitEnergyLayerY[k]) + get<0>(HitEnergyLayerX[l]) / 3.;
                              std::vector<double> tempEnergyPosXY3{HitEnergy, HitX, HitY, Z_plane};
                              HitEnergyPosXY.emplace_back(tempEnergyPosXY3);

                              countHitLayerY[i] += 1.;
                              countHitLayerY[j] += 1.;
                              countHitLayerY[k] += 1.;
                              countHitLayerX[l] += 3.;
                            }
                        }
                    }
                }
            }
        }

      for(int i = HitEnergyLayerX.size() - 1; i >= 0; --i)
        {
          if(countHitLayerX[i] > 0.1)
            {
              HitEnergyLayerX.erase(HitEnergyLayerX.begin() + i);
              countHitLayerX.erase(countHitLayerX.begin() + i);
            }
        }

      for(int j = HitEnergyLayerY.size() - 1; j >= 0; --j)
        {
          if(countHitLayerY[j] > 0.1)
            {
              HitEnergyLayerY.erase(HitEnergyLayerY.begin() + j);
              countHitLayerY.erase(countHitLayerY.begin() + j);
            }
        }
    }

  if(HitEnergyPosXY.size() >= 2)
    {
      for(size_t i = 0; i < HitEnergyPosXY.size() - 1; ++i)
        {
          for(size_t j = i + 1; j < HitEnergyPosXY.size(); ++j)
            {
              if((HitEnergyPosXY[i][0] == HitEnergyPosXY[j][0]) && (HitEnergyPosXY[i][1] == HitEnergyPosXY[j][1]) &&
                 (HitEnergyPosXY[i][2] == HitEnergyPosXY[j][2]))
                {
                  HitEnergyPosXY.erase(HitEnergyPosXY.begin() + j);
                  --j;
                }
            }
        }
    }

#endif
}


template <class Out>
void TPrimaryVertex<Out>::TranslationZ_Target_System(double Target_PosZ)
{
  double Dist_to_TransZ = Target_PosZ - 25.; // Reference target PosZ = 25.
  Zo_target += Dist_to_TransZ;
  Zf_target += Dist_to_TransZ;
/* 
  Z_plane_Si1 += Dist_to_TransZ;
  Z_plane_Si2 += Dist_to_TransZ;

  SiliconHitsSD_Si1.TranslationZ_Si_System(Target_PosZ);
  SiliconHitsSD_Si2.TranslationZ_Si_System(Target_PosZ);
  
  if(Si_4layers)
    {
      SiliconHitsSD_Si3.TranslationZ_Si_System(Target_PosZ);
      SiliconHitsSD_Si4.TranslationZ_Si_System(Target_PosZ);
    }
*/
}

template <class Out>
void TPrimaryVertex<Out>::CloseDist(std::vector<double>& BeamHit1, std::vector<double>& BeamHit2,
                               std::vector<double>& TrackHit1, std::vector<double>& TrackHit2, double& distance,
                               double& z)
{
  std::vector<double> BeamPar{BeamHit2[1] - BeamHit1[1], BeamHit2[2] - BeamHit1[2], BeamHit2[3] - BeamHit1[3]};
  std::vector<double> TrackPar{TrackHit2[1] - TrackHit1[1], TrackHit2[2] - TrackHit1[2], TrackHit2[3] - TrackHit1[3]};

  std::vector<double> n{BeamPar[1] * TrackPar[2] - BeamPar[2] * TrackPar[1],
                        BeamPar[2] * TrackPar[0] - BeamPar[0] * TrackPar[2],
                        BeamPar[0] * TrackPar[1] - BeamPar[1] * TrackPar[0]};
  std::vector<double> n1{BeamPar[1] * n[2] - BeamPar[2] * n[1], BeamPar[2] * n[0] - BeamPar[0] * n[2],
                         BeamPar[0] * n[1] - BeamPar[1] * n[0]};
  std::vector<double> n2{TrackPar[1] * n[2] - TrackPar[2] * n[1], TrackPar[2] * n[0] - TrackPar[0] * n[2],
                         TrackPar[0] * n[1] - TrackPar[1] * n[0]};

  std::vector<double> c1(3, 0.);
  std::vector<double> c2(3, 0.);

  for(size_t i = 0; i < 3; ++i)
    {
      c1[i] = BeamHit1[i + 1] + ((TrackHit1[1] - BeamHit1[1]) * n2[0] + (TrackHit1[2] - BeamHit1[2]) * n2[1] +
                                 (TrackHit1[3] - BeamHit1[3]) * n2[2]) /
                                    (BeamPar[0] * n2[0] + BeamPar[1] * n2[1] + BeamPar[2] * n2[2]) * BeamPar[i];
      c2[i] = TrackHit1[i + 1] + -((TrackHit1[1] - BeamHit1[1]) * n1[0] + (TrackHit1[2] - BeamHit1[2]) * n1[1] +
                                   (TrackHit1[3] - BeamHit1[3]) * n1[2]) /
                                     (TrackPar[0] * n1[0] + TrackPar[1] * n1[1] + TrackPar[2] * n1[2]) * TrackPar[i];
    }

  distance =
      std::sqrt((c2[0] - c1[0]) * (c2[0] - c1[0]) + (c2[1] - c1[1]) * (c2[1] - c1[1]) + (c2[2] - c1[2]) * (c2[2] - c1[2]));
  z = (c1[2] + c2[2]) / 2.;
}

template <class Out>
double TPrimaryVertex<Out>::f_function(std::vector<double>& Hit1, std::vector<double>& Hit2, std::vector<double>& PosXYZ)
{
  double slope_x     = (Hit2[1] - Hit1[1]) / (Hit2[3] - Hit1[3]);
  double intercept_x = Hit2[1] - slope_x * Hit2[3];
  double slope_y     = (Hit2[2] - Hit1[2]) / (Hit2[3] - Hit1[3]);
  double intercept_y = Hit2[2] - slope_y * Hit2[3];

  double distanceStepX = 2. * boxDistBeamXY / static_cast<double>(NstepsdiscretXY - 1);
  double sigma2        = std::pow(distanceStepX, 2.) / 12.;

  double f = exp(-0.5 *
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
      v = k_factor * f_vector[f_vector.size() - 1] + sum_f -
          (k_factor * std::pow(f_vector[f_vector.size() - 1], 2.) + sum_f2) /
              (k_factor * f_vector[f_vector.size() - 1] + sum_f);
    }

  return v;
}

template <class Out>
void TPrimaryVertex<Out>::SpaceDiscretization(double& Xi, double& Xf, size_t& NstepsX, double& Yi, double& Yf,
                                         size_t& NstepsY, double& Zi, double& Zf, size_t& NstepsZ, size_t& border,
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

  for(size_t i = 0 + border; i < NstepsX - border; ++i)
    {
      double PosX = Xi + i * (Xf - Xi) / static_cast<double>(NstepsX - 1);

      for(size_t j = 0 + border; j < NstepsY - border; ++j)
        {
          double PosY = Yi + j * (Yf - Yi) / static_cast<double>(NstepsY - 1);

          for(size_t k = 0 + border; k < NstepsZ - border; ++k)
            {
              double PosZ = Zi + k * (Zf - Zi) / static_cast<double>(NstepsZ - 1);

              std::vector<double> temp_PosXYZ{PosX, PosY, PosZ};
              PosXYZ.emplace_back(temp_PosXYZ);
            }
        }
    }
}

template <class Out>
void TPrimaryVertex<Out>::HitstoTracks(std::vector<std::vector<double> >& HitEnergyPosXY_Si1,
                                  std::vector<std::vector<double> >& HitEnergyPosXY_Si2,
                                  std::vector<double>& BeamHit1, std::vector<double>& BeamHit2,
                                  std::vector<std::vector<std::vector<double> > >& CandidateTracks)
{
  double distance = 0.;
  double z        = 0.;

  for(size_t i = 0; i < HitEnergyPosXY_Si1.size(); ++i)
    {
      for(size_t j = 0; j < HitEnergyPosXY_Si2.size(); ++j)
        {
          CloseDist(BeamHit1, BeamHit2, HitEnergyPosXY_Si1[i], HitEnergyPosXY_Si2[j], distance, z);

          if((z > Zo_target - ErrorDistTarget) && (z < Zf_target + ErrorDistTarget) &&
             (distance < MaxClosestDistance) &&
             (std::abs(HitEnergyPosXY_Si1[i][0] - HitEnergyPosXY_Si2[j][0]) < MaxEnergyDiffSilicons))
            {
              std::vector<std::vector<double> > temp_CandidateTracks{HitEnergyPosXY_Si1[i], HitEnergyPosXY_Si2[j]};
              CandidateTracks.emplace_back(temp_CandidateTracks);

              double temp_thetaCand = std::atan(std::sqrt(std::pow((HitEnergyPosXY_Si1[i][1] - HitEnergyPosXY_Si2[j][1]), 2.)
                                                          + std::pow((HitEnergyPosXY_Si1[i][2] - HitEnergyPosXY_Si2[j][2]), 2.))
                                                  / (HitEnergyPosXY_Si1[i][3] - HitEnergyPosXY_Si2[j][3])) * 180. / M_PI;
              LocalHisto.h_Acc_ThetaCandidates->Fill(temp_thetaCand, 1.);
            }
        }
    }
}

template <class Out>
void TPrimaryVertex<Out>::TrackstoVertexPosition(std::vector<std::vector<std::vector<double> > >& CandidateTracks,
                                            std::vector<double>& BeamHit1, std::vector<double>& BeamHit2,
                                            std::vector<double>& InteractionPointRecons,
                                            std::vector<double>& f_values_IP)
{
  double InteractionPoint_rand_X = (BeamHit1[1] + BeamHit2[1]) / 2.;
  double InteractionPoint_rand_Y = (BeamHit1[2] + BeamHit2[2]) / 2.;

  std::vector<double> temp_f(CandidateTracks.size() + 1, 0.);

  double V    = 0.;
  double Vnew = 0.;

  double Xi            = InteractionPoint_rand_X - boxDistBeamXY;
  double Xf            = InteractionPoint_rand_X + boxDistBeamXY;
  double distanceStepX = (Xf - Xi) / static_cast<double>(NstepsdiscretXY - 1);

  double Yi            = InteractionPoint_rand_Y - boxDistBeamXY;
  double Yf            = InteractionPoint_rand_Y + boxDistBeamXY;
  double distanceStepY = (Yf - Yi) / static_cast<double>(NstepsdiscretXY - 1);

  double Zi            = Zo_target;
  double Zf            = Zf_target;
  double distanceStepZ = (Zf - Zi) / static_cast<double>(NstepsdiscretZ - 1);

  size_t border = 0;
  std::vector<std::vector<double> > PosXYZ{};
  SpaceDiscretization(Xi, Xf, NstepsdiscretXY, Yi, Yf, NstepsdiscretXY, Zi, Zf, NstepsdiscretZ, border, PosXYZ);

  border = 1;

  for(size_t k = 0; k < nTimesDiscretization; ++k)
    {
      if(k != 0)
        {
          Xi            = InteractionPointRecons[0] - distanceStepX;
          Xf            = InteractionPointRecons[0] + distanceStepX;
          distanceStepX = (Xf - Xi) / static_cast<double>(Nstepsdiscretbox - 1);

          Yi            = InteractionPointRecons[1] - distanceStepY;
          Yf            = InteractionPointRecons[1] + distanceStepY;
          distanceStepY = (Yf - Yi) / static_cast<double>(Nstepsdiscretbox - 1);

          Zi            = InteractionPointRecons[2] - distanceStepZ;
          Zf            = InteractionPointRecons[2] + distanceStepZ;
          distanceStepZ = (Zf - Zi) / static_cast<double>(Nstepsdiscretbox - 1);

          SpaceDiscretization(Xi, Xf, Nstepsdiscretbox, Yi, Yf, Nstepsdiscretbox, Zi, Zf, Nstepsdiscretbox, border,
                              PosXYZ);
        }

      for(size_t i = 0; i < PosXYZ.size(); ++i)
        {
          std::vector<double> temp_PosXYZ = PosXYZ[i];
          for(size_t j = 0; j < CandidateTracks.size(); ++j)
            {
              temp_f[j] = f_function(CandidateTracks[j][0], CandidateTracks[j][1], temp_PosXYZ);
            }
          temp_f[CandidateTracks.size()] = f_function(BeamHit1, BeamHit2, temp_PosXYZ);

          Vnew = V_function(temp_f);

          if(Vnew > V)
            {
              V = Vnew;

              f_values_IP = temp_f;

              InteractionPointRecons[0] = temp_PosXYZ[0];
              InteractionPointRecons[1] = temp_PosXYZ[1];
              InteractionPointRecons[2] = temp_PosXYZ[2];
            }
        }
    }
}

template <class Out>
void TPrimaryVertex<Out>::CovarianceMatrix(std::vector<std::vector<std::vector<double> > >& CandidateTracks,
                                      std::vector<double>& BeamHit1, std::vector<double>& BeamHit2,
                                      std::vector<double>& InteractionPointAverage, std::vector<double>& f_values_IP,
                                      std::vector<std::vector<double> >& CovMatrix)
{
  std::vector<std::vector<std::vector<double> > > DecisiveTracks;

  for(size_t i = 0; i < CandidateTracks.size(); ++i)
    {
      if(f_values_IP[i] > min_f_value)
        {
          DecisiveTracks.emplace_back(CandidateTracks[i]);
        }
    }

  LocalHisto.h_nTrackCandidates->Fill(DecisiveTracks.size(), "Decisive", 1.);

  std::vector<double> temp_IP(3, 0.);
  std::vector<double> temp_f_values_IP(DecisiveTracks.size() + 1, 0.);

  std::vector<std::vector<double> > variations_IP{};
  double nDimensions = 4. * DecisiveTracks.size() + 2.;

  double average_X_IP = 0.;
  double average_Y_IP = 0.;
  double average_Z_IP = 0.;

  for(size_t idSilicon = 0; idSilicon <= 1; ++idSilicon)
    {
      for(size_t variationsign = 0; variationsign <= 1; ++variationsign)
        {
          for(size_t idStrip = 1; idStrip <= 2; ++idStrip)
            {
              for(size_t idDecisiveTrack = 0; idDecisiveTrack < DecisiveTracks.size(); ++idDecisiveTrack)
                {
                  double value_not_variation = DecisiveTracks[idDecisiveTrack][idSilicon][idStrip];
                  DecisiveTracks[idDecisiveTrack][idSilicon][idStrip] += std::pow(-1., static_cast<double>(variationsign)) *
                                                                         std::sqrt(nDimensions) *
                                                                         sigma_Si[idSilicon];

                  TrackstoVertexPosition(DecisiveTracks, BeamHit1, BeamHit2, temp_IP, temp_f_values_IP);
                  variations_IP.emplace_back(temp_IP);

                  average_X_IP += temp_IP[0];
                  average_Y_IP += temp_IP[1];
                  average_Z_IP += temp_IP[2];

                  DecisiveTracks[idDecisiveTrack][idSilicon][idStrip] = value_not_variation;
                }
            }

          double value_not_variation2 = BeamHit1[idSilicon + 1];
          double value_not_variation3 = BeamHit2[idSilicon + 1];

          BeamHit1[idSilicon + 1] +=
              std::pow(-1., static_cast<double>(variationsign)) * std::sqrt(nDimensions) * sigma_beam;
          BeamHit2[idSilicon + 1] +=
              std::pow(-1., static_cast<double>(variationsign)) * std::sqrt(nDimensions) * sigma_beam;

          TrackstoVertexPosition(DecisiveTracks, BeamHit1, BeamHit2, temp_IP, temp_f_values_IP);
          variations_IP.emplace_back(temp_IP);

          average_X_IP += temp_IP[0];
          average_Y_IP += temp_IP[1];
          average_Z_IP += temp_IP[2];

          BeamHit1[idSilicon + 1] = value_not_variation2;
          BeamHit2[idSilicon + 1] = value_not_variation3;
        }
    }

  if(std::abs(variations_IP.size() - 2. * nDimensions) > 0.01)
    {
      std::cout << "Error with dimensions covariance matrix\t" << variations_IP.size() << "\t" << 2. * nDimensions
                << "\n";
    }

  average_X_IP = average_X_IP / static_cast<double>(variations_IP.size());
  average_Y_IP = average_Y_IP / static_cast<double>(variations_IP.size());
  average_Z_IP = average_Z_IP / static_cast<double>(variations_IP.size());

  InteractionPointAverage[0] = average_X_IP;
  InteractionPointAverage[1] = average_Y_IP;
  InteractionPointAverage[2] = average_Z_IP;

  std::vector<double> zero_vect(3, 0.);
  std::vector<std::vector<double> > zero_matrix(3, zero_vect);

  CovMatrix = zero_matrix;

  for(size_t i = 0; i < variations_IP.size(); ++i)
    {
      std::vector<double> temp_var_IP{variations_IP[i][0] - average_X_IP, variations_IP[i][1] - average_Y_IP,
                                      variations_IP[i][2] - average_Z_IP};

      for(size_t j = 0; j < temp_var_IP.size(); ++j)
        {
          for(size_t k = 0; k < temp_var_IP.size(); ++k)
            {
              CovMatrix[j][k] += temp_var_IP[j] * temp_var_IP[k];
            }
        }
    }

  for(size_t j = 0; j < CovMatrix.size(); ++j)
    {
      for(size_t k = 0; k < CovMatrix.size(); ++k)
        {
          CovMatrix[j][k] /= static_cast<double>(variations_IP.size());
        }
    }
}

template <class Out>
void TPrimaryVertex<Out>::HitstoDecayTracks(std::vector<std::vector<double> >& HitEnergyPosXY_Si1,
                                       std::vector<std::vector<double> >& HitEnergyPosXY_Si2,
                                       std::vector<double>& BeamHit1, std::vector<double>& BeamHit2,
                                       std::vector<double>& InteractionPointRecons,
                                       std::vector<std::vector<std::vector<double> > >& CandidateDecayTracks)
{
  double distance = 0.;
  double z        = 0.;

  for(size_t i = 0; i < HitEnergyPosXY_Si1.size(); ++i)
    {
      for(size_t j = 0; j < HitEnergyPosXY_Si2.size(); ++j)
        {
          CloseDist(BeamHit1, BeamHit2, HitEnergyPosXY_Si1[i], HitEnergyPosXY_Si2[j], distance, z);

          if((z > InteractionPointRecons[2] + MinDistIPDecay) && (z < Z_plane_Si1 - MinDistIPDecay) &&
             (distance < 1.5 * MaxClosestDistance) &&
             (std::abs(HitEnergyPosXY_Si1[i][0] - HitEnergyPosXY_Si2[j][0]) < MaxEnergyDiffSilicons))
            {
              std::vector<std::vector<double> > temp_CandidateDecayTracks{HitEnergyPosXY_Si1[i], HitEnergyPosXY_Si2[j]};
              CandidateDecayTracks.emplace_back(temp_CandidateDecayTracks);
            }
        }
    }
}

template <class Out>
void TPrimaryVertex<Out>::DecayTrackstoDecayPosition(std::vector<std::vector<std::vector<double> > >& CandidateTracks,
                                                std::vector<double>& InteractionPointRecons,
                                                std::vector<double>& DecayPositionRecons)
{
  std::vector<double> temp_f(CandidateTracks.size(), 0.);

  double V    = 0.;
  double Vnew = 0.;

  double Xi            = InteractionPointRecons[0] - boxDistDecayXY;
  double Xf            = InteractionPointRecons[0] + boxDistDecayXY;
  double distanceStepX = (Xf - Xi) / static_cast<double>(NstepsdiscretDecay - 1);

  double Yi            = InteractionPointRecons[1] - boxDistDecayXY;
  double Yf            = InteractionPointRecons[1] + boxDistDecayXY;
  double distanceStepY = (Yf - Yi) / static_cast<double>(NstepsdiscretDecay - 1);

  double Zi            = InteractionPointRecons[2] + MinDistIPDecay;
  double Zf            = Z_plane_Si1 - MinDistIPDecay;
  double distanceStepZ = (Zf - Zi) / static_cast<double>(NstepsdiscretDecay - 1);

  size_t border = 0;
  std::vector<std::vector<double> > PosXYZ{};
  SpaceDiscretization(Xi, Xf, NstepsdiscretDecay, Yi, Yf, NstepsdiscretDecay, Zi, Zf, NstepsdiscretDecay, border,
                      PosXYZ);

  border = 1;

  for(size_t k = 0; k < nTimesDiscretization; ++k)
    {
      if(k != 0)
        {
          Xi            = DecayPositionRecons[0] - distanceStepX;
          Xf            = DecayPositionRecons[0] + distanceStepX;
          distanceStepX = (Xf - Xi) / static_cast<double>(Nstepsdiscretbox - 1);

          Yi            = DecayPositionRecons[1] - distanceStepY;
          Yf            = DecayPositionRecons[1] + distanceStepY;
          distanceStepY = (Yf - Yi) / static_cast<double>(Nstepsdiscretbox - 1);

          Zi            = DecayPositionRecons[2] - distanceStepZ;
          Zf            = DecayPositionRecons[2] + distanceStepZ;
          distanceStepZ = (Zf - Zi) / static_cast<double>(Nstepsdiscretbox - 1);

          SpaceDiscretization(Xi, Xf, Nstepsdiscretbox, Yi, Yf, Nstepsdiscretbox, Zi, Zf, Nstepsdiscretbox, border,
                              PosXYZ);
        }

      for(size_t i = 0; i < PosXYZ.size(); ++i)
        {
          std::vector<double> temp_PosXYZ = PosXYZ[i];
          double alpha                    = std::atan(std::sqrt(std::pow((InteractionPointRecons[0] - temp_PosXYZ[0]), 2.) +
                                   std::pow((InteractionPointRecons[1] - temp_PosXYZ[1]), 2.)) /
                              (InteractionPointRecons[2] - temp_PosXYZ[2])) * 180. / M_PI;

          for(size_t j = 0; j < CandidateTracks.size(); ++j)
            {
              temp_f[j] = f_function(CandidateTracks[j][0], CandidateTracks[j][1], temp_PosXYZ);
            }

          Vnew = V_function(temp_f) * exp(-k_alpha_factor * alpha);

          if(Vnew > V)
            {
              V = Vnew;

              DecayPositionRecons[0] = temp_PosXYZ[0];
              DecayPositionRecons[1] = temp_PosXYZ[1];
              DecayPositionRecons[2] = temp_PosXYZ[2];
            }
        }
    }
}

template <class Out>
void TPrimaryVertex<Out>::nGoodEventsCounter(
    std::vector<std::vector<double> >& HitEnergyPosXY,
    std::vector<std::tuple<double, double, double, size_t, double, double, std::string> >& HitEnergyPosXYreal,
    double& widthStrip, size_t& nGoodrecons)
{
  for(size_t i = 0; i < HitEnergyPosXY.size(); ++i)
    {
      for(size_t j = 0; j < HitEnergyPosXYreal.size(); ++j)
        {
          if((std::abs(HitEnergyPosXY[i][1] - get<1>(HitEnergyPosXYreal[j])) < widthStrip / 2.) &&
             (std::abs(HitEnergyPosXY[i][2] - get<2>(HitEnergyPosXYreal[j])) < widthStrip / 2.))
            {
              nGoodrecons += 1;
            }
        }
    }
}

template <class Out>
void TPrimaryVertex<Out>::RealHitstoRealTracks(
    std::vector<std::tuple<double, double, double, size_t, double, double, std::string> >& HitEnergyPosXYreal_Si1,
    std::vector<std::tuple<double, double, double, size_t, double, double, std::string> >& HitEnergyPosXYreal_Si2,
    std::vector<std::vector<std::vector<double> > >& RealTracks)
{
  Z_plane_Si1 = 27.; // in cm
  Z_plane_Si2 = 30.; // in cm

  double particletype;

  for(size_t i = 0; i < HitEnergyPosXYreal_Si1.size(); ++i)
    {
      for(size_t j = 0; j < HitEnergyPosXYreal_Si2.size(); ++j)
        {
          if(get<3>(HitEnergyPosXYreal_Si1[i]) == get<3>(HitEnergyPosXYreal_Si2[j]))
            {
              if((get<6>(HitEnergyPosXYreal_Si1[i]).compare("H3L") == 0) ||
                 (get<6>(HitEnergyPosXYreal_Si1[i]).compare("H4L") == 0) ||
                 (get<6>(HitEnergyPosXYreal_Si1[i]).compare("He5L") == 0))
                {
                  particletype = 1.;
                }
              else if(get<6>(HitEnergyPosXYreal_Si1[i]).compare("He3") == 0)
                {
                  particletype = 2.;
                }
              else if(get<6>(HitEnergyPosXYreal_Si1[i]).compare("pi-") == 0)
                {
                  particletype = 3.;
                }
              else
                {
                  particletype = 0.;
                }

              std::vector<double> temp_PosHit1{get<0>(HitEnergyPosXYreal_Si1[i]),
                                               get<1>(HitEnergyPosXYreal_Si1[i]),
                                               get<2>(HitEnergyPosXYreal_Si1[i]),
                                               Z_plane_Si1,
                                               get<4>(HitEnergyPosXYreal_Si1[i]),
                                               get<5>(HitEnergyPosXYreal_Si1[i]),
                                               particletype,
                                               static_cast<double>(get<3>(HitEnergyPosXYreal_Si1[i]))};
              std::vector<double> temp_PosHit2{get<0>(HitEnergyPosXYreal_Si2[j]),
                                               get<1>(HitEnergyPosXYreal_Si2[j]),
                                               get<2>(HitEnergyPosXYreal_Si2[j]),
                                               Z_plane_Si2,
                                               get<4>(HitEnergyPosXYreal_Si2[j]),
                                               get<5>(HitEnergyPosXYreal_Si2[j]),
                                               particletype,
                                               static_cast<double>(get<3>(HitEnergyPosXYreal_Si2[j]))};

              std::vector<std::vector<double> > temp_RealTracks{temp_PosHit1, temp_PosHit2};
              RealTracks.emplace_back(temp_RealTracks);

              LocalHisto.h_EnergyDiffSilicons->Fill(
                  get<0>(HitEnergyPosXYreal_Si1[i]) - get<0>(HitEnergyPosXYreal_Si2[j]), "Real track", 1.);
            }
        }
    }
}

template <class Out>
void TPrimaryVertex<Out>::nGoodTracksCounter(std::vector<std::vector<std::vector<double> > >& CandidateTracks,
                                        std::vector<std::vector<std::vector<double> > >& RealTracks,
                                        size_t& nGoodTracks, std::vector<size_t>& goodCandidateTracks)
{
  for(size_t i = 0; i < RealTracks.size(); ++i)
    {
      double realPosX_Si1 = RealTracks[i][0][1];
      double realPosY_Si1 = RealTracks[i][0][2];

      double realPosX_Si2 = RealTracks[i][1][1];
      double realPosY_Si2 = RealTracks[i][1][2];

      double thetareal    = std::atan(std::sqrt(std::pow((realPosY_Si2 - realPosY_Si1), 2.) + std::pow((realPosX_Si2 - realPosX_Si1), 2.)) /
                              (Z_plane_Si2 - Z_plane_Si1)) * 180. / M_PI;

      LocalHisto.h_thetaTracks->Fill(thetareal, "Real", 1.);

      for(size_t j = 0; j < CandidateTracks.size(); ++j)
        {
          double candPosX_Si1 = CandidateTracks[j][0][1];
          double candPosY_Si1 = CandidateTracks[j][0][2];

          double candPosX_Si2 = CandidateTracks[j][1][1];
          double candPosY_Si2 = CandidateTracks[j][1][2];

          double DiffEnergySilicons = CandidateTracks[j][0][0] - CandidateTracks[j][1][0];
          double theta = std::atan(std::sqrt(std::pow((candPosY_Si2 - candPosY_Si1), 2.) + std::pow((candPosX_Si2 - candPosX_Si1), 2.)) /
                              (Z_plane_Si2 - Z_plane_Si1)) * 180. / M_PI;

          if((std::abs(realPosX_Si1 - candPosX_Si1) < widthStrip_Si1 / 2.) &&
             (std::abs(realPosY_Si1 - candPosY_Si1) < widthStrip_Si1 / 2.) &&
             (std::abs(realPosX_Si2 - candPosX_Si2) < widthStrip_Si2 / 2.) &&
             (std::abs(realPosY_Si2 - candPosY_Si2) < widthStrip_Si2 / 2.))
            {
              nGoodTracks += 1;
              goodCandidateTracks[j] = 1;

              double thetaResol = thetareal - theta;

              LocalHisto.h_EnergyDiffSilicons->Fill(DiffEnergySilicons, "Good With-energy", 1.);
              LocalHisto.h_thetaTracks->Fill(theta, "Recons", 1.);
              LocalHisto.h_thetaResol->Fill(thetaResol, 1.);
            }
          else
            LocalHisto.h_EnergyDiffSilicons->Fill(DiffEnergySilicons, "False With-energy", 1.);
        }
    }
}

template <class Out>
void TPrimaryVertex<Out>::nForwardTracksCounter(std::vector<std::vector<std::vector<double> > >& CandidateTracks,
                                           size_t& nForwardTracks, std::vector<size_t>& forwardCandidateTracks)
{

  nForwardTracks = 0;

  for(size_t i = 0; i < CandidateTracks.size(); ++i)
    {
      double candPosX_Si1 = CandidateTracks[i][0][1];
      double candPosY_Si1 = CandidateTracks[i][0][2];

      double candPosX_Si2 = CandidateTracks[i][1][1];
      double candPosY_Si2 = CandidateTracks[i][1][2];

      double theta = std::atan(std::sqrt(std::pow((candPosY_Si2 - candPosY_Si1), 2.) + std::pow((candPosX_Si2 - candPosX_Si1), 2.)) /
                          (Z_plane_Si2 - Z_plane_Si1)) * 180. / M_PI;

      if(theta <= 5.)
        {
          nForwardTracks += 1;
          forwardCandidateTracks[i] = 1;
        }
    }
}

template class TPrimaryVertex<MCAnaEventG4Sol>;
template class TPrimaryVertex<Ana_WasaEvent>;
