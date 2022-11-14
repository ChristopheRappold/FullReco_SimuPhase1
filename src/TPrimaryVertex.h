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


class SiliconHits_SDpad
{

/* For 2 layers: 1p2pXiYj
  double Z_plane_Si1x    = 27.; // in cm
  double Z_plane_Si1y    = 27.05; // in cm
  size_t combineStrips_Si1 = 2; // power of 2
  double widthStrip_Si1 = 0.008 * combineStrips_Si1; // in cm
  double lenghtSi_Si1   = 2.048; // in cm
  double thicknessSi_Si1 = 0.0285; // in cm
  bool restrict_actlenght_Si1 = false;
  double actlenghtX_Si1 = 50.; // in cm
  double actlenghtY_Si1 = 50.; // in cm
  bool restrict_gapcenter_Si1 = false;
  double gapcenter_Si1 = 0.05; //in cm
  bool ifveto_Si1 = false;
  std::vector<int> inactwiresX_Si1 {};
  std::vector<int> inactwiresY_Si1 {};
  std::vector<size_t> PadDistribution_Si1 {2, 2}; // {nPads_X, n_Pads_Y}
  bool ifreflectionX_Si1 = false;
  bool ifreflectionY_Si1 = false;
  bool ifphirot_Si1 = false;
  double phirot_Si1 = 0.; //in degrees

  double Z_plane_Si2x    = 28.; // in cm
  double Z_plane_Si2y    = 28.05; // in cm
  size_t combineStrips_Si2 = 2; // power of 2
  double widthStrip_Si2 = 0.008 * combineStrips_Si2; // in cm
  double lenghtSi_Si2   = 2.048; // in cm
  double thicknessSi_Si2 = 0.0285; // in cm
  bool restrict_actlenght_Si2 = false;
  double actlenghtX_Si2 = 50.; // in cm
  double actlenghtY_Si2 = 50.; // in cm
  bool restrict_gapcenter_Si2 = false;
  double gapcenter_Si2 = 0.05; //in cm
  bool ifveto_Si2 = false;
  std::vector<int> inactwiresX_Si2 {};
  std::vector<int> inactwiresY_Si2 {};
  std::vector<size_t> PadDistribution_Si2 {2, 3}; // {nPads_X, n_Pads_Y}
  bool ifreflectionX_Si2 = false;
  bool ifreflectionY_Si2 = false;
  bool ifphirot_Si2 = false;
  double phirot_Si2 = 0.; //in degrees
*/

// For 4 layers: pad1234a
  double Z_plane_Si1x    = 27.4; // in cm
  double Z_plane_Si1y    = 27.45; // in cm
  size_t combineStrips_Si1 = 4; // power of 2
  double widthStrip_Si1 = 0.008 * combineStrips_Si1; // in cm
  double lenghtSi_Si1   = 2.048; // in cm
  double thicknessSi_Si1 = 0.0285; // in cm
  bool restrict_actlenght_Si1 = false;
  double actlenghtX_Si1 = 50.; // in cm
  double actlenghtY_Si1 = 50.; // in cm
  bool restrict_gapcenter_Si1 = false;
  double gapcenter_Si1 = 0.188*2.; //in cm
  bool ifveto_Si1 = false;
  std::vector<int> inactwiresX_Si1 {};
  std::vector<int> inactwiresY_Si1 {};
  std::vector<size_t> PadDistribution_Si1 {2, 2}; // {nPads_X, n_Pads_Y}
  bool ifreflectionX_Si1 = false;
  bool ifreflectionY_Si1 = false;
  bool ifphirot_Si1 = true;
  double phirot_Si1 = -45.; //in degrees

  double Z_plane_Si2x    = 28.2; // in cm
  double Z_plane_Si2y    = 28.25; // in cm
  size_t combineStrips_Si2 = 4; // power of 2
  double widthStrip_Si2 = 0.008 * combineStrips_Si2; // in cm
  double lenghtSi_Si2   = 2.048; // in cm
  double thicknessSi_Si2 = 0.0285; // in cm
  bool restrict_actlenght_Si2 = false;
  double actlenghtX_Si2 = 50.; // in cm
  double actlenghtY_Si2 = 50.; // in cm
  bool restrict_gapcenter_Si2 = false;
  double gapcenter_Si2 = 0.188*2.; //in cm
  bool ifveto_Si2 = false;
  std::vector<int> inactwiresX_Si2 {};
  std::vector<int> inactwiresY_Si2 {};
  std::vector<size_t> PadDistribution_Si2 {2, 2}; // {nPads_X, n_Pads_Y}
  bool ifreflectionX_Si2 = false;
  bool ifreflectionY_Si2 = false;
  bool ifphirot_Si2 = true;
  double phirot_Si2 = -45.; //in degrees


  double Z_plane_Si3x    = 26.475; // in cm
  double Z_plane_Si3y    = 26.525; // in cm
  size_t combineStrips_Si3 = 1; // power of 2
  double widthStrip_Si3 = 0.008 * combineStrips_Si3; // in cm
  double lenghtSi_Si3   = 2.048; // in cm
  double thicknessSi_Si3 = 0.0285; // in cm
  bool restrict_actlenght_Si3 = false;
  double actlenghtX_Si3 = 50.; // in cm
  double actlenghtY_Si3 = 50.; // in cm
  bool restrict_gapcenter_Si3 = false;
  double gapcenter_Si3 = 0.188*2.; //in cm
  bool ifveto_Si3 = false;
  std::vector<int> inactwiresX_Si3 {};
  std::vector<int> inactwiresY_Si3 {};
  std::vector<size_t> PadDistribution_Si3 {1, 1}; // {nPads_X, n_Pads_Y}
  bool ifreflectionX_Si3 = false;
  bool ifreflectionY_Si3 = true;
  bool ifphirot_Si3 = false;
  double phirot_Si3 = 0.; //in degrees

  double Z_plane_Si4x    = 27.075; // in cm
  double Z_plane_Si4y    = 27.125; // in cm
  size_t combineStrips_Si4 = 1; // power of 2
  double widthStrip_Si4 = 0.008 * combineStrips_Si4; // in cm
  double lenghtSi_Si4   = 2.048; // in cm
  double thicknessSi_Si4 = 0.0285; // in cm
  bool restrict_actlenght_Si4 = false;
  double actlenghtX_Si4 = 50.; // in cm
  double actlenghtY_Si4 = 50.; // in cm
  bool restrict_gapcenter_Si4 = false;
  double gapcenter_Si4 = 0.188*2.; //in cm
  bool ifveto_Si4 = false;
  std::vector<int> inactwiresX_Si4 {};
  std::vector<int> inactwiresY_Si4 {};
  std::vector<size_t> PadDistribution_Si4 {1, 1}; // {nPads_X, n_Pads_Y}
  bool ifreflectionX_Si4 = false;
  bool ifreflectionY_Si4 = true;
  bool ifphirot_Si4 = false;
  double phirot_Si4 = 0.; //in degrees

  double MaxEnergyDiffStrips = 0.2;   // in MeV
  double MaxEnergyMultiplicity = 0.06; // in MeV

  size_t Si_station;
  double Z_plane;
  double widthStrip;
  double lenghtSi;
  double thicknessSi;
  bool restrict_actlenght;
  double actlenghtX;
  double actlenghtY;
  bool restrict_gapcenter;
  double gapcenter;
  size_t combineStrips;
  bool ifveto;
  std::vector<int> inactwiresX; 
  std::vector<int> inactwiresY;
  std::vector<size_t> PadDistribution;
  bool ifreflectionX;
  bool ifreflectionY;
  bool ifphirot;
  double phirot;
  int nStrips;
  size_t nPads;

  int nRejectPad = 0;

public:
  SiliconHits_SDpad(size_t idSilicon)
  {
    if(idSilicon == 1)
      {
        Z_plane    = (Z_plane_Si1x + Z_plane_Si1y)/2.;
        widthStrip = widthStrip_Si1;
        lenghtSi   = lenghtSi_Si1;
        thicknessSi = thicknessSi_Si1;
        restrict_actlenght = restrict_actlenght_Si1;
        actlenghtX = actlenghtX_Si1;
        actlenghtY = actlenghtY_Si1;
        restrict_gapcenter = restrict_gapcenter_Si1;
        gapcenter = gapcenter_Si1;
        combineStrips = combineStrips_Si1;
        ifveto = ifveto_Si1;
        inactwiresX = inactwiresX_Si1;
        inactwiresY = inactwiresY_Si1;
        PadDistribution = PadDistribution_Si1;
        ifreflectionX = ifreflectionX_Si1;
        ifreflectionY = ifreflectionY_Si1;
        ifphirot = ifphirot_Si1;
        phirot = phirot_Si1;
      }
    else if(idSilicon == 2)
      {
        Z_plane    = (Z_plane_Si2x + Z_plane_Si2y)/2.;
        widthStrip = widthStrip_Si2;
        lenghtSi   = lenghtSi_Si2;
        thicknessSi = thicknessSi_Si2;
        restrict_actlenght = restrict_actlenght_Si2;
        actlenghtX = actlenghtX_Si2;
        actlenghtY = actlenghtY_Si2;
        restrict_gapcenter = restrict_gapcenter_Si2;
        gapcenter = gapcenter_Si2;
        combineStrips = combineStrips_Si2;
        ifveto = ifveto_Si2;
        inactwiresX = inactwiresX_Si2;
        inactwiresY = inactwiresY_Si2;
        PadDistribution = PadDistribution_Si2;
        ifreflectionX = ifreflectionX_Si2;
        ifreflectionY = ifreflectionY_Si2;
        ifphirot = ifphirot_Si2;
        phirot = phirot_Si2;
      }
    else if(idSilicon == 3)
      {
        Z_plane    = (Z_plane_Si3x + Z_plane_Si3y)/2.;
        widthStrip = widthStrip_Si3;
        lenghtSi   = lenghtSi_Si3;
        thicknessSi = thicknessSi_Si3;
        restrict_actlenght = restrict_actlenght_Si3;
        actlenghtX = actlenghtX_Si3;
        actlenghtY = actlenghtY_Si3;
        restrict_gapcenter = restrict_gapcenter_Si3;
        gapcenter = gapcenter_Si3;
        combineStrips = combineStrips_Si3;
        ifveto = ifveto_Si3;
        inactwiresX = inactwiresX_Si3;
        inactwiresY = inactwiresY_Si3;
        PadDistribution = PadDistribution_Si3;
        ifreflectionX = ifreflectionX_Si3;
        ifreflectionY = ifreflectionY_Si3;
        ifphirot = ifphirot_Si3;
        phirot = phirot_Si3;
      }
    else if(idSilicon == 4)
      {
        Z_plane    = (Z_plane_Si4x + Z_plane_Si4y)/2.;
        widthStrip = widthStrip_Si4;
        lenghtSi   = lenghtSi_Si4;
        thicknessSi = thicknessSi_Si4;
        restrict_actlenght = restrict_actlenght_Si4;
        actlenghtX = actlenghtX_Si4;
        actlenghtY = actlenghtY_Si4;
        restrict_gapcenter = restrict_gapcenter_Si4;
        gapcenter = gapcenter_Si4;
        combineStrips = combineStrips_Si4;
        ifveto = ifveto_Si4;
        inactwiresX = inactwiresX_Si4;
        inactwiresY = inactwiresY_Si4;
        PadDistribution = PadDistribution_Si4;
        ifreflectionX = ifreflectionX_Si4;
        ifreflectionY = ifreflectionY_Si4;
        ifphirot = ifphirot_Si4;
        phirot = phirot_Si4;
      }
    else
      std::cout << "Wrong idSilicon: it must be 1,2,3,4 for Silicon1,2,3,4\n";

    Si_station = idSilicon;
    nStrips = static_cast<int>(lenghtSi / widthStrip);
    nPads = PadDistribution[0] * PadDistribution[1];
  }

  size_t size_type_abs(size_t a, size_t b)
  {
    if(a >= b)
      return a - b;
    else
      return b - a;
  }

  void SignalstoHits_SDpad(std::vector<std::tuple<double, size_t> > HitEnergyLayerX,
                        std::vector<std::tuple<double, size_t> > HitEnergyLayerY,
                        std::vector<std::vector<double> >& HitEnergyPosXY);

  bool IsSamePad(size_t& layerX_id, size_t& layerY_id);

  void originPad(size_t& layerX_id, double& originPadX, double& originPadY);

  void CombineStrips(std::vector<std::tuple<double, size_t> >& HitEnergyLayer);

  int CountRejectPad()
  {
    return nRejectPad;
  }

  void Clear()
  {
    nRejectPad = 0;
  }

  void TranslationZ_Si_System(double Target_PosZ);

};


/* For Hamamatsu silicons
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
*/

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

  void TranslationZ_Target_System(double Target_PosZ);

  SiliconHits_SDpad SiliconHitsSD_Si1;
  SiliconHits_SDpad SiliconHitsSD_Si2;

  bool Si_4layers = true;
  SiliconHits_SDpad SiliconHitsSD_Si3;
  SiliconHits_SDpad SiliconHitsSD_Si4;

  double Zo_target = 24.; // in cm
  double Zf_target = 26.; // in cm

  double Z_plane_Si1x   = 27.4; // in cm
  double Z_plane_Si1y   = 27.45; // in cm
  double Z_plane_Si1    = (Z_plane_Si1x + Z_plane_Si1y)/2.;
  size_t combineStrips_Si1 = 4; // power of 2
  double widthStrip_Si1 = 0.008 * combineStrips_Si1; // in cm
  double sigma_Si1      = widthStrip_Si1 / std::sqrt(12.);
  double lenghtSi_Si1   = 2.048; // in cm
  int nStrips_Si1 = static_cast<int>(lenghtSi_Si1 / widthStrip_Si1);

  double Z_plane_Si2x    = 28.2; // in cm
  double Z_plane_Si2y    = 28.25; // in cm
  double Z_plane_Si2     = (Z_plane_Si2x + Z_plane_Si2y)/2.;  // in cm
  size_t combineStrips_Si2 = 4; // power of 2
  double widthStrip_Si2 = 0.008 * combineStrips_Si2; // in cm
  double sigma_Si2       = widthStrip_Si2 / std::sqrt(12.);
  double lenghtSi_Si2   = 2.048; // in cm
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
    TH2F* h_nEventsRealRejectPad_Si1;


    TH1F* h_HitMultiplicity_Si2;
    TH1F* h_HitMultiplicityRecons_Si2;
    TH1F* h_HitMultiplicityDiff_Si2;
    TH2F* h_HitMultiplicityDiffNHits_Si2;

    TH1F* h_EnergyDiffStrips_Si2;

    TH1F* h_nEventsGoodrecons_Si2;
    TH1F* h_nEventsGhost_Si2;
    TH2F* h_nEventsGoodreconsGhost_Si2;
    TH2F* h_nEventsRealGoodrecons_Si2;
    TH2F* h_nEventsRealRejectPad_Si2;

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

    TH1F* h_Acc_ThetaCandidates;
    TH1F* h_Acc_ThetaAllReal;
    
    TH2F* h_nCandidatesRealTracks;
    TH2F* h_nCandidatesRealTracks_IfRecons;

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
