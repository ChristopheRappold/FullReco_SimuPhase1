#ifndef TPRIMARYVERTEX
#define TPRIMARYVERTEX

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "ReturnRes.hh"
#include "TDataProcess.h"
#include "THyphiAttributes.h"
#include "TRandom3.h"

typedef TDataProcess<FullRecoEvent, MCAnaEventG4Sol> TDataProcessInterface;

class SiliconHits
{
  double Z_plane_Si1    = 27.;  // in cm
  double widthStrip_Si1 = 0.03; // in cm
  double lenghtSi_Si1   = 4.;   // in cm

  double Z_plane_Si2    = 30.;  // in cm
  double widthStrip_Si2 = 0.03; // in cm
  double lenghtSi_Si2   = 6.;   // in cm

  double EnergyThreshold     = 0.001; // in MeV
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

class TPrimaryVertex final : public TDataProcessInterface
{
public:
  const THyphiAttributes& att;

  TPrimaryVertex(const THyphiAttributes& attr);
  ~TPrimaryVertex();

  // int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator()(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree) override;

private:
  int Exec(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;

  int FinderPrimaryVertex(FullRecoEvent& RecoEvent);
  void simulHitstoRealHits(
      FullRecoEvent& REvent,
      std::vector<std::tuple<double, double, double, size_t, double, double, std::string> >& HitEnergyPosXYreal,
      int id_det_x, int id_det_y, TH1F* h_Diff);

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

  SiliconHits SiliconHits_Si1;
  SiliconHits SiliconHits_Si2;

  double Zo_target = 24.5; // in cm
  double Zf_target = 25.5; // in cm
  double Xo_target = -1.;
  double Xf_target = 1.;
  double Yo_target = -1.;
  double Yf_target = +1.;

  double Z_plane_Si1    = 27.;  // in cm
  double widthStrip_Si1 = 0.03; // in cm
  double sigma_Si1      = widthStrip_Si1 / std::sqrt(12.);
  double thicknessSi_Si1 = 0.03; // cm
  double lenghtSi_Si1    = 4.;   // in cm
  int nStrips_Si1        = (int)(lenghtSi_Si1 / widthStrip_Si1);

  double Z_plane_Si2     = 30.;  // in cm
  double widthStrip_Si2  = 0.03; // in cm
  double sigma_Si2       = widthStrip_Si2 / std::sqrt(12.);
  double thicknessSi_Si2 = 0.03; // cm
  double lenghtSi_Si2    = 6.;   // in cm
  int nStrips_Si2        = (int)(lenghtSi_Si2 / widthStrip_Si2);

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

  double EnergyThreshold       = 0.001; // in MeV
  double MaxEnergyDiffStrips   = 0.2;   // in MeV
  double MaxEnergyDiffSilicons = 0.1;   // in MeV

  double MaxDistTracks = 0.2;

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

    TH1F* h_HitMultiplicity_Si2;
    TH1F* h_HitMultiplicityRecons_Si2;
    TH1F* h_HitMultiplicityDiff_Si2;
    TH2F* h_HitMultiplicityDiffNHits_Si2;

    TH1F* h_EnergyDiffStrips_Si2;

    TH1F* h_nEventsGoodrecons_Si2;
    TH1F* h_nEventsGhost_Si2;
    TH2F* h_nEventsGoodreconsGhost_Si2;
    TH2F* h_nEventsRealGoodrecons_Si2;

    TH2F* h_EnergyStripEnergyTotalReal;
    TH2F* h_EnergyStripEnergyTotal;
    TH2F* h_EnergyDiffSilicons;

    TH2F* h_EnergyDepositionMother;
    TH2F* h_EnergyDepositionDaughters;

    TH2F* h_nTrackCandidates;
    TH2F* h_DistanceBeamTracks;
    TH2F* h_PosZBeamTracks;
    TH2F* h_thetaTracks;

    TH1F* h_nHypernucleiTrack;
    TH1F* h_fvalues;

    TH1F* h_InteractionPointDistance;
    TH1F* h_InteractionPointDistanceX;
    TH1F* h_InteractionPointDistanceY;
    TH1F* h_InteractionPointDistanceZ;

    TH1F* h_InteractionPointDistanceX_pull;
    TH1F* h_InteractionPointDistanceY_pull;
    TH1F* h_InteractionPointDistanceZ_pull;

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
