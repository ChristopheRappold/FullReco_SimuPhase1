#ifndef TPRIMARYVERTEX
#define TPRIMARYVERTEX

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "ReturnRes.hh"
#include "TDataProcess.h"
#include "THyphiAttributes.h"
#include "TRandom3.h"

//typedef TDataProcess<FullRecoEvent, MCAnaEventG4Sol> TDataProcessInterface;
// template<class Out>
// struct TDataProcessInterfaceImp { using type = TDataProcess<FullRecoEvent, Out>; } ;

// template<class Out>
// using TDataProcessInterface = typename TDataProcessInterfaceImp<Out>::type;

template<class Out>
using TDataProcessInterface = TDataProcess<FullRecoEvent, Out>;

/*
class SiliconHits
{
  double Z_plane_Si1    = 27.;  // in cm
  double widthStrip_Si1 = 0.03; // in cm
  double lenghtSi_Si1   = 4.;   // in cm

  double Z_plane_Si2    = 30.;  // in cm
  double widthStrip_Si2 = 0.03; // in cm
  double lenghtSi_Si2   = 6.;   // in cm

  double MaxEnergyDiffStrips = 0.2;   // in MeV

  double MaxEnergyMultiplicity = 0.03; // in MeV

  double Z_plane;
  double widthStrip;
  double lenghtSi;
  int nStrips;

public:
  SiliconHits(size_t idSilicon)
  {
    if(idSilicon == 1)
      {
        Z_plane    = Z_plane_Si1;
        widthStrip = widthStrip_Si1;
        lenghtSi   = lenghtSi_Si1;
      }
    else if(idSilicon == 2)
      {
        Z_plane    = Z_plane_Si2;
        widthStrip = widthStrip_Si2;
        lenghtSi   = lenghtSi_Si2;
      }
    else
      std::cout << "Wrong idSilicon: it must be 1 for Silicon1 or 2 for Silicon2\n";

    nStrips = static_cast<int>(lenghtSi / widthStrip);
  }

  size_t size_type_abs(size_t a, size_t b)
  {
    if(a >= b)
      return a - b;
    else
      return b - a;
  }

  void SignalstoHits(std::vector<std::tuple<double, size_t> >& HitEnergyLayerX,
                      std::vector<std::tuple<double, size_t> >& HitEnergyLayerY,
                      std::vector<std::vector<double> >& HitEnergyPosXY);
};
*/

class SiliconHits_SD
{
  double Z_plane_Si1x    = 27.; // in cm
  double Z_plane_Si1y    = 27.05; // in cm
  size_t combineStrips_Si1 = 1; // power of 2
  double widthStrip_Si1 = 0.019 * combineStrips_Si1; // in cm
  double lenghtSi_Si1   = 9.728; // in cm
  double thicknessSi_Si1 = 0.032; // in cm
  bool restrict_actlenght_Si1 = true;
  double actlenghtX_Si1 = 3.648; // in cm
  double actlenghtY_Si1 = 3.648; // in cm
  bool restrict_gapcenter_Si1 = true;
  double gapcenter_Si1 = 0.05; //in cm
  bool ifveto_Si1 = false;
  std::vector<int> inactwiresX_Si1 {254, 255, 256, 257, 1278, 1279, 1280, 1281}; 
  std::vector<int> inactwiresY_Si1 {254, 255, 256, 257, 1278, 1279, 1280, 1281};

  double Z_plane_Si2x    = 30.; // in cm
  double Z_plane_Si2y    = 30.05; // in cm
  size_t combineStrips_Si2 = 2; // power of 2
  double widthStrip_Si2 = 0.019 * combineStrips_Si2; // in cm
  double lenghtSi_Si2   = 9.728; // in cm
  double thicknessSi_Si2 = 0.032; // in cm
  bool restrict_actlenght_Si2 = true;
  double actlenghtX_Si2 = 7.296; // in cm
  double actlenghtY_Si2 = 7.296; // in cm
  bool restrict_gapcenter_Si2 = true;
  double gapcenter_Si2 = 0.05; //in cm
  bool ifveto_Si2 = false;
  std::vector<int> inactwiresX_Si2 {254, 255, 256, 257, 1278, 1279, 1280, 1281};
  std::vector<int> inactwiresY_Si2 {254, 255, 256, 257, 1278, 1279, 1280, 1281};

  double MaxEnergyDiffStrips = 0.4;   // in MeV
  double MaxEnergyMultiplicity = 0.06; // in MeV


  double Z_plane;
  double widthStrip;
  double lenghtSi;
  bool restrict_actlenght;
  double actlenghtX;
  double actlenghtY;
  bool restrict_gapcenter;
  double gapcenter;
  size_t combineStrips;
  bool ifveto;
  std::vector<int> inactwiresX; 
  std::vector<int> inactwiresY;
  int nStrips;

  int nRejectCuadrant = 0;

public:
  SiliconHits_SD(size_t idSilicon)
  {
    if(idSilicon == 1)
      {
        Z_plane    = (Z_plane_Si1x + Z_plane_Si1y)/2.;
        widthStrip = widthStrip_Si1;
        lenghtSi   = lenghtSi_Si1;
        restrict_actlenght = restrict_actlenght_Si1;
        actlenghtX = actlenghtX_Si1;
        actlenghtY = actlenghtY_Si1;
        restrict_gapcenter = restrict_gapcenter_Si1;
        gapcenter = gapcenter_Si1;
        combineStrips = combineStrips_Si1;
        ifveto = ifveto_Si1;
        inactwiresX = inactwiresX_Si1;
        inactwiresY = inactwiresY_Si1;
      }
    else if(idSilicon == 2)
      {
        Z_plane    = (Z_plane_Si2x + Z_plane_Si2y)/2.;
        widthStrip = widthStrip_Si2;
        lenghtSi   = lenghtSi_Si2;
        restrict_actlenght = restrict_actlenght_Si2;
        actlenghtX = actlenghtX_Si2;
        actlenghtY = actlenghtY_Si2;
        restrict_gapcenter = restrict_gapcenter_Si2;
        gapcenter = gapcenter_Si2;
        combineStrips = combineStrips_Si2;
        ifveto = ifveto_Si2;
        inactwiresX = inactwiresX_Si2;
        inactwiresY = inactwiresY_Si2;
      }
    else
      std::cout << "Wrong idSilicon: it must be 1 for Silicon1 or 2 for Silicon2\n";

    nStrips = static_cast<int>(lenghtSi / widthStrip);
  }

  size_t size_type_abs(size_t a, size_t b)
  {
    if(a >= b)
      return a - b;
    else
      return b - a;
  }

  void SignalstoHits_SD(std::vector<std::tuple<double, size_t> > HitEnergyLayerX,
                        std::vector<std::tuple<double, size_t> > HitEnergyLayerY,
                        std::vector<std::vector<double> >& HitEnergyPosXY);

  bool IsValidCuadrant(size_t& layerX_id, size_t& layerY_id);

  void CombineStrips(std::vector<std::tuple<double, size_t> >& HitEnergyLayer);

  int CountRejectCuadrant()
  {
    return nRejectCuadrant;
  }

  void Clear()
  {
    nRejectCuadrant = 0;
  }
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
  void simulHitstoRealHits(
      FullRecoEvent& REvent,
      std::vector<std::tuple<double, double, double, size_t, double, double, std::string> >& HitEnergyPosXYreal,
      int id_det_x, int id_det_y, TH1F* h_Diff, std::vector<std::tuple<size_t,TVector3,TVector3>>& HitIdMomPos);

  void MFcheck(std::vector<std::tuple<size_t,TVector3,TVector3>>& HitIdMomPos_Si1,
               std::vector<std::tuple<size_t,TVector3,TVector3>>& HitIdMomPos_Si2,
               std::array<double,3> InteractionPoint);

  void nGoodEventsCounter(
      std::vector<std::vector<double> >& HitEnergyPosXY,
      std::vector<std::tuple<double, double, double, size_t, double, double, std::string> >& HitEnergyPosXYreal,
      double& widthStrip, size_t& nGoodrecons);

  void RealHitstoRealTracks(
      std::vector<std::tuple<double, double, double, size_t, double, double, std::string> >& HitEnergyPosXYreal_Si1,
      std::vector<std::tuple<double, double, double, size_t, double, double, std::string> >& HitEnergyPosXYreal_Si2,
      std::vector<std::vector<std::vector<double> > >& RealTracks);

  void nGoodTracksCounter(std::vector<std::vector<std::vector<double> > >& CandidateTracks,
                          std::vector<std::vector<std::vector<double> > >& RealTracks, size_t& nGoodTracks,
                          std::vector<size_t>& goodCandidateTracks);

  void nForwardTracksCounter(std::vector<std::vector<std::vector<double> > >& CandidateTracks, size_t& nForwardTracks,
                             std::vector<size_t>& forwardCandidateTracks);

  void CloseDist(std::vector<double>& BeamHit1, std::vector<double>& BeamHit2, std::vector<double>& TrackHit1,
                 std::vector<double>& TrackHit2, double& distance, double& z);

  double f_function(std::vector<double>& Hit1, std::vector<double>& Hit2, std::vector<double>& PosXYZ);

  double V_function(std::vector<double>& f_vector);

  void SpaceDiscretization(double& Xi, double& Xf, size_t& NstepsX, double& Yi, double& Yf, size_t& NstepsY, double& Zi,
                           double& Zf, size_t& NstepsZ, size_t& border, std::vector<std::vector<double> >& PosXYZ);
  
  void HitstoTracks(std::vector<std::vector<double> >& HitEnergyPosXY_Si1,
                    std::vector<std::vector<double> >& HitEnergyPosXY_Si2, std::vector<double>& BeamHit1,
                    std::vector<double>& BeamHit2, std::vector<std::vector<std::vector<double> > >& CandidateTracks);

  void TrackstoVertexPosition(std::vector<std::vector<std::vector<double> > >& CandidateTracks,
                              std::vector<double>& BeamHit1, std::vector<double>& BeamHit2,
                              std::vector<double>& InteractionPointRecons, std::vector<double>& f_values_IP);

  void CovarianceMatrix(std::vector<std::vector<std::vector<double> > >& CandidateTracks, std::vector<double>& BeamHit1,
                        std::vector<double>& BeamHit2, std::vector<double>& InteractionPointAverage,
                        std::vector<double>& f_values_IP, std::vector<std::vector<double> >& CovMatrix);

  void HitstoDecayTracks(std::vector<std::vector<double> >& HitEnergyPosXY_Si1,
                         std::vector<std::vector<double> >& HitEnergyPosXY_Si2, std::vector<double>& BeamHit1,
                         std::vector<double>& BeamHit2, std::vector<double>& InteractionPointRecons,
                         std::vector<std::vector<std::vector<double> > >& CandidateDecayTracks);

  void DecayTrackstoDecayPosition(std::vector<std::vector<std::vector<double> > >& CandidateTracks,
                                  std::vector<double>& InteractionPointRecons,
                                  std::vector<double>& DecayPositionRecons);

  //SiliconHits SiliconHits_Si1;
  //SiliconHits SiliconHits_Si2;

  SiliconHits_SD SiliconHitsSD_Si1;
  SiliconHits_SD SiliconHitsSD_Si2;

  double Zo_target = 24.5; // in cm
  double Zf_target = 25.5; // in cm

/*
  double Z_plane_Si1    = 27.;  // in cm
  double widthStrip_Si1 = 0.03; // in cm
  double sigma_Si1      = widthStrip_Si1 / std::sqrt(12.);
  double lenghtSi_Si1    = 4.;   // in cm

  double Z_plane_Si2     = 30.;  // in cm
  double widthStrip_Si2  = 0.03; // in cm
  double sigma_Si2       = widthStrip_Si2 / std::sqrt(12.);
  double lenghtSi_Si2    = 6.;   // in cm
*/

  double Z_plane_Si1x   = 27.; // in cm
  double Z_plane_Si1y   = 27.05; // in cm
  double Z_plane_Si1    = (Z_plane_Si1x + Z_plane_Si1y)/2.;
  double widthStrip_Si1 = 0.019; // in cm
  double sigma_Si1      = widthStrip_Si1 / std::sqrt(12.);
  double lenghtSi_Si1    = 9.728;   // in cm
  int nStrips_Si1 = static_cast<int>(lenghtSi_Si1 / widthStrip_Si1);

  double Z_plane_Si2x    = 30.; // in cm
  double Z_plane_Si2y    = 30.05; // in cm
  double Z_plane_Si2     = (Z_plane_Si2x + Z_plane_Si2y)/2.;  // in cm
  double widthStrip_Si2  = 0.038; // in cm
  double sigma_Si2       = widthStrip_Si2 / std::sqrt(12.);
  double lenghtSi_Si2    = 9.728;   // in cm
  int nStrips_Si2 = static_cast<int>(lenghtSi_Si2 / widthStrip_Si2);


  double ErrorDistTarget    = 0.2;  // in cm
  double MaxClosestDistance = 0.04; // in cm

  size_t NstepsdiscretXY      = 5;
  size_t NstepsdiscretZ       = 11;
  size_t Nstepsdiscretbox     = 5;
  size_t nTimesDiscretization = 6;
  double boxDistBeamXY        = 0.4;

  double boxDistDecayXY     = 0.75; // in cm
  size_t NstepsdiscretDecay = 11;

  double k_factor       = 3.; // multiplies the beam track f function value
  double k_alpha_factor = 5.;
  double MinDistIPDecay = 0.15;

  double randInteractionPointXY = 0.01; // in cm

  double min_f_value = 0.5;
  std::vector<double> sigma_Si{sigma_Si1, sigma_Si2};
  double sigma_beam = 2. * randInteractionPointXY / std::sqrt(12.);

  double EnergyThreshold     = 0.001; // in MeV
  double MaxEnergyDiffSilicons = 0.3;   // in MeV

  TRandom3* rand;

  struct LocalHists
  {
    TH1F* h_HitMultiplicity_Si1;
    TH1F* h_HitMultiplicityRecons_Si1;
    TH1F* h_HitMultiplicityDiff_Si1;
    TH2F* h_HitMultiplicityDiffNHits_Si1;

    TH1F* h_EnergyDiffStrips_Si1;

    TH1F* h_nEventsGoodrecons_Si1;
    TH1F* h_nEventsGhost_Si1;
    TH2F* h_nEventsGoodreconsGhost_Si1;
    TH2F* h_nEventsRealGoodrecons_Si1;
    TH2F* h_nEventsRealRejectCuadrant_Si1;


    TH1F* h_HitMultiplicity_Si2;
    TH1F* h_HitMultiplicityRecons_Si2;
    TH1F* h_HitMultiplicityDiff_Si2;
    TH2F* h_HitMultiplicityDiffNHits_Si2;

    TH1F* h_EnergyDiffStrips_Si2;

    TH1F* h_nEventsGoodrecons_Si2;
    TH1F* h_nEventsGhost_Si2;
    TH2F* h_nEventsGoodreconsGhost_Si2;
    TH2F* h_nEventsRealGoodrecons_Si2;
    TH2F* h_nEventsRealRejectCuadrant_Si2;

    TH1F* h_MFCheck_Theta_MomSi1MomSi2;
    TH1F* h_MFCheck_Dist_MomSi1HitSi2;
    TH1F* h_MFCheck_Dist_MomSi2HitSi1;
    TH1F* h_MFCheck_Dist_MomSi1HitIP;
    TH1F* h_MFCheck_Dist_MomSi2HitIP;

    TH2F* h_EnergyStripEnergyTotalReal;
    TH2F* h_EnergyStripEnergyTotal;
    TH2F* h_EnergyDiffSilicons;

    TH2F* h_EnergyDepositionMother;
    TH2F* h_EnergyDepositionDaughters;

    TH2F* h_nTrackCandidates;
    TH2F* h_DistanceBeamTracks;
    TH2F* h_PosZBeamTracks;
    TH2F* h_thetaTracks;
    TH1F* h_thetaResol;

    TH1F* h_nHypernucleiTrack;
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

    TH1F* h_IP_DecayDistance;
    TH1F* h_IP_DecayDistanceX;
    TH1F* h_IP_DecayDistanceY;
    TH1F* h_IP_DecayDistanceZ;

    TH1F* h_DecayPositionDistance;
    TH1F* h_DecayPositionDistanceX;
    TH1F* h_DecayPositionDistanceY;
    TH1F* h_DecayPositionDistanceZ;

    TH1F* h_PrimVtxstats;
    TH2F* h_PrimStatus;
  };
  LocalHists LocalHisto;
};

#endif
