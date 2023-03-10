#include "TBuildWASACalibrationLayerPlane.h"

#include <list>
#include <map>
#include <set>
#include <vector>

//#define DEBUG_BUILD
//#define DEBUG_BUILD3

using namespace std;

TBuildWASACalibrationLayerPlane::TBuildWASACalibrationLayerPlane(const THyphiAttributes& attribut)
    : TDataBuilder("build_det"), att(attribut)
{
  att._logger->info("TBuildWASACalibrationLayerPlane::TBuildWASACalibrationLayerPlane");

  par = std::make_unique<ParaManager>(att.map_ParamFiles);

  std::string volMDCfirst = gGeoManager->GetVolume("INNER")->GetNode(0)->GetVolume()->GetName();
  if(volMDCfirst == "MD01")
    offsetGeoNameID_MDC = 0;
  if(volMDCfirst == "SOL")
    offsetGeoNameID_MDC = 1;

  auto listNodes = gGeoManager->GetVolume("INNER")->GetNodes();
  int index_lastMDC = -1, index_firstPSCE = -1;
  for(int i=0;i<listNodes->GetEntries();++i)
    {
      std::string tempName(listNodes->At(i)->GetName());
      if(tempName == "MD17_1")
        index_lastMDC = i;
      if(tempName == "PSCE_1")
        index_firstPSCE = i;
    }
  offsetGeoNameID_PSCE = index_firstPSCE - index_lastMDC + offsetGeoNameID_MDC + 17 -1;
}

TBuildWASACalibrationLayerPlane::~TBuildWASACalibrationLayerPlane() {}

#ifdef ROOT6
ReturnRes::InfoM TBuildWASACalibrationLayerPlane::operator()(const EventWASAUnpack& event, FullRecoEvent& RecoEvent,
                                                             Ana_WasaEvent* OutTree)
{
  int result = Exec(event, RecoEvent, OutTree);

  return SoftExit(result);
}
#else
ReturnRes::InfoM TBuildWASACalibrationLayerPlane::operator()(const EventWASAUnpack& event, FullRecoEvent& RecoEvent,
                                                             Ana_WasaEvent* OutTree)
{
  int result = Exec(event, RecoEvent, OutTree);

  return SoftExit(result);
}

ReturnRes::InfoM TBuildWASACalibrationLayerPlane::operator()(const EventWASAUnpack& event,
                                                          FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
{
  return ReturnRes::BuildError;
}

#endif
void TBuildWASACalibrationLayerPlane::SelectHists()
{
  LocalHisto.h_Builderstats = AnaHisto->CloneAndRegister(AnaHisto->h_Builderstats);

  for(int i = 0; i < 7; ++i)
    {
      for(int j = 0; j < 3; ++j)
        {
          LocalHisto.h10[i][j]         = AnaHisto->CloneAndRegister(AnaHisto->h10[i][j]);
          LocalHisto.h11[i][j]         = AnaHisto->CloneAndRegister(AnaHisto->h11[i][j]);
          LocalHisto.h12[i][j]         = AnaHisto->CloneAndRegister(AnaHisto->h12[i][j]);
          LocalHisto.h13[i][j]         = AnaHisto->CloneAndRegister(AnaHisto->h13[i][j]);
          LocalHisto.h14[i][j]         = AnaHisto->CloneAndRegister(AnaHisto->h14[i][j]);
          LocalHisto.h15[i][j]         = AnaHisto->CloneAndRegister(AnaHisto->h15[i][j]);
          LocalHisto.hfiber_13_0[i][j] = AnaHisto->CloneAndRegister(AnaHisto->hfiber_13_0[i][j]);
          LocalHisto.hfiber_13_1[i][j] = AnaHisto->CloneAndRegister(AnaHisto->hfiber_13_1[i][j]);
          LocalHisto.hfiber_13_2[i][j] = AnaHisto->CloneAndRegister(AnaHisto->hfiber_13_2[i][j]);
          LocalHisto.hfiber_13_3[i][j] = AnaHisto->CloneAndRegister(AnaHisto->hfiber_13_3[i][j]);
          LocalHisto.hfiber_13_4[i][j] = AnaHisto->CloneAndRegister(AnaHisto->hfiber_13_4[i][j]);
          LocalHisto.h75[i][j]         = AnaHisto->CloneAndRegister(AnaHisto->h75[i][j]);
        }
      LocalHisto.h16[i]   = AnaHisto->CloneAndRegister(AnaHisto->h16[i]);
      LocalHisto.h17[i]   = AnaHisto->CloneAndRegister(AnaHisto->h17[i]);
      LocalHisto.h17_2[i] = AnaHisto->CloneAndRegister(AnaHisto->h17_2[i]);
    }
  LocalHisto.h18_3_1 = AnaHisto->CloneAndRegister(AnaHisto->h18_3_1);
  LocalHisto.h18_3_2 = AnaHisto->CloneAndRegister(AnaHisto->h18_3_2);
  LocalHisto.h18_3_3 = AnaHisto->CloneAndRegister(AnaHisto->h18_3_3);
  LocalHisto.h18_3_4 = AnaHisto->CloneAndRegister(AnaHisto->h18_3_4);
  LocalHisto.h18_3_5 = AnaHisto->CloneAndRegister(AnaHisto->h18_3_5);
  LocalHisto.h18_3_6 = AnaHisto->CloneAndRegister(AnaHisto->h18_3_6);
  LocalHisto.h18_3_7 = AnaHisto->CloneAndRegister(AnaHisto->h18_3_7);
  LocalHisto.h18_3_8 = AnaHisto->CloneAndRegister(AnaHisto->h18_3_8);
  for(int i=0; i<3 ; ++i)
    for(int j=0; j<3; ++j)
      for(int k=0; k<2; ++k)
        LocalHisto.h51[i][j][k] = AnaHisto->CloneAndRegister(AnaHisto->h51[i][j][k]);

  //MFT12
  LocalHisto.hfiber_1_1     = AnaHisto->CloneAndRegister(AnaHisto->hfiber_1_1);
  LocalHisto.hfiber_1_2     = AnaHisto->CloneAndRegister(AnaHisto->hfiber_1_2);
  LocalHisto.hfiber_1_3     = AnaHisto->CloneAndRegister(AnaHisto->hfiber_1_3);
  LocalHisto.hfiber_1_4     = AnaHisto->CloneAndRegister(AnaHisto->hfiber_1_4);
  LocalHisto.hfiber_1_5     = AnaHisto->CloneAndRegister(AnaHisto->hfiber_1_5);
  LocalHisto.hfiber_1_6     = AnaHisto->CloneAndRegister(AnaHisto->hfiber_1_6);
  LocalHisto.hfiber_1_7     = AnaHisto->CloneAndRegister(AnaHisto->hfiber_1_7);
  LocalHisto.hfiber_1_9     = AnaHisto->CloneAndRegister(AnaHisto->hfiber_1_9);
  LocalHisto.hfiber_2_1_1   = AnaHisto->CloneAndRegister(AnaHisto->hfiber_2_1_1);
  LocalHisto.hfiber_2_1_2   = AnaHisto->CloneAndRegister(AnaHisto->hfiber_2_1_2);
  LocalHisto.hfiber_2_2_1   = AnaHisto->CloneAndRegister(AnaHisto->hfiber_2_2_1);
  LocalHisto.hfiber_2_2_2   = AnaHisto->CloneAndRegister(AnaHisto->hfiber_2_2_2);
  LocalHisto.hfiber_2_3     = AnaHisto->CloneAndRegister(AnaHisto->hfiber_2_3);
  LocalHisto.hfiber_3_0     = AnaHisto->CloneAndRegister(AnaHisto->hfiber_3_0);
  LocalHisto.hfiber_3_0_2   = AnaHisto->CloneAndRegister(AnaHisto->hfiber_3_0_2);
  LocalHisto.hfiber_6_1     = AnaHisto->CloneAndRegister(AnaHisto->hfiber_6_1);
  LocalHisto.hfiber_6_2     = AnaHisto->CloneAndRegister(AnaHisto->hfiber_6_2);
  LocalHisto.hfiber_6_3     = AnaHisto->CloneAndRegister(AnaHisto->hfiber_6_3);
  LocalHisto.hfiber_6_4     = AnaHisto->CloneAndRegister(AnaHisto->hfiber_6_4);
  LocalHisto.hfiber_12_1_1  = AnaHisto->CloneAndRegister(AnaHisto->hfiber_12_1_1);
  LocalHisto.hfiber_12_2_1  = AnaHisto->CloneAndRegister(AnaHisto->hfiber_12_2_1);
  LocalHisto.hfiber_12_3_1  = AnaHisto->CloneAndRegister(AnaHisto->hfiber_12_3_1);
  LocalHisto.hfiber_12_1_2  = AnaHisto->CloneAndRegister(AnaHisto->hfiber_12_1_2);
  LocalHisto.hfiber_12_2_2  = AnaHisto->CloneAndRegister(AnaHisto->hfiber_12_2_2);
  LocalHisto.hfiber_12_3_2  = AnaHisto->CloneAndRegister(AnaHisto->hfiber_12_3_2);

  //DFT12
  LocalHisto.hfiber_4_1   = AnaHisto->CloneAndRegister(AnaHisto->hfiber_4_1);
  LocalHisto.hfiber_4_2_1 = AnaHisto->CloneAndRegister(AnaHisto->hfiber_4_2_1);
  LocalHisto.hfiber_4_3_1 = AnaHisto->CloneAndRegister(AnaHisto->hfiber_4_3_1);
  LocalHisto.hfiber_4_4_1 = AnaHisto->CloneAndRegister(AnaHisto->hfiber_4_4_1);
  LocalHisto.hfiber_4_5_1 = AnaHisto->CloneAndRegister(AnaHisto->hfiber_4_5_1);
  LocalHisto.hfiber_4_2_2 = AnaHisto->CloneAndRegister(AnaHisto->hfiber_4_2_2);
  LocalHisto.hfiber_4_3_2 = AnaHisto->CloneAndRegister(AnaHisto->hfiber_4_3_2);
  LocalHisto.hfiber_4_4_2 = AnaHisto->CloneAndRegister(AnaHisto->hfiber_4_4_2);
  LocalHisto.hfiber_4_5_2 = AnaHisto->CloneAndRegister(AnaHisto->hfiber_4_5_2);
  LocalHisto.hfiber_4_1_3 = AnaHisto->CloneAndRegister(AnaHisto->hfiber_4_1_3);
  LocalHisto.hfiber_4_2_3 = AnaHisto->CloneAndRegister(AnaHisto->hfiber_4_2_3);
  LocalHisto.hfiber_4_3_3 = AnaHisto->CloneAndRegister(AnaHisto->hfiber_4_3_3);
  LocalHisto.hfiber_4_4_3 = AnaHisto->CloneAndRegister(AnaHisto->hfiber_4_4_3);
  LocalHisto.hfiber_4_5_3 = AnaHisto->CloneAndRegister(AnaHisto->hfiber_4_5_3);
  LocalHisto.hfiber_5_1   = AnaHisto->CloneAndRegister(AnaHisto->hfiber_5_1);
  LocalHisto.hfiber_5_2   = AnaHisto->CloneAndRegister(AnaHisto->hfiber_5_2);
  LocalHisto.hfiber_5_3   = AnaHisto->CloneAndRegister(AnaHisto->hfiber_5_3);
  LocalHisto.hfiber_5_4   = AnaHisto->CloneAndRegister(AnaHisto->hfiber_5_4);
  LocalHisto.hfiber_5_5   = AnaHisto->CloneAndRegister(AnaHisto->hfiber_5_5);
  LocalHisto.hfiber_5_6   = AnaHisto->CloneAndRegister(AnaHisto->hfiber_5_6);
  LocalHisto.hfiber_5_7   = AnaHisto->CloneAndRegister(AnaHisto->hfiber_5_7);

  LocalHisto.hpsb_0_1 = AnaHisto->CloneAndRegister(AnaHisto->hpsb_0_1);
  LocalHisto.hpsb_0_2 = AnaHisto->CloneAndRegister(AnaHisto->hpsb_0_2);
  LocalHisto.hpsb_0_3 = AnaHisto->CloneAndRegister(AnaHisto->hpsb_0_3);
  LocalHisto.hpsb_0_4 = AnaHisto->CloneAndRegister(AnaHisto->hpsb_0_4);
  LocalHisto.hpsb_1_1 = AnaHisto->CloneAndRegister(AnaHisto->hpsb_1_1);
  for(int i = 0; i < 46; ++i)
    {
      LocalHisto.hpsb_2[i] = AnaHisto->CloneAndRegister(AnaHisto->hpsb_2[i]);
      LocalHisto.hpsb_3[i] = AnaHisto->CloneAndRegister(AnaHisto->hpsb_3[i]);
      LocalHisto.hpsb_4[i] = AnaHisto->CloneAndRegister(AnaHisto->hpsb_4[i]);
    }
  LocalHisto.h76 = AnaHisto->CloneAndRegister(AnaHisto->h76);

  LocalHisto.hpsfe_0_1 = AnaHisto->CloneAndRegister(AnaHisto->hpsfe_0_1);
  LocalHisto.hpsfe_0_2 = AnaHisto->CloneAndRegister(AnaHisto->hpsfe_0_2);
  LocalHisto.hpsfe_0_3 = AnaHisto->CloneAndRegister(AnaHisto->hpsfe_0_3);
  LocalHisto.hpsfe_0_4 = AnaHisto->CloneAndRegister(AnaHisto->hpsfe_0_4);

  LocalHisto.hpsbe_0_1 = AnaHisto->CloneAndRegister(AnaHisto->hpsbe_0_1);
  LocalHisto.hpsbe_0_2 = AnaHisto->CloneAndRegister(AnaHisto->hpsbe_0_2);
  LocalHisto.hpsbe_0_3 = AnaHisto->CloneAndRegister(AnaHisto->hpsbe_0_3);
  LocalHisto.hpsbe_0_4 = AnaHisto->CloneAndRegister(AnaHisto->hpsbe_0_4);
  LocalHisto.hpsbe_1_0 = AnaHisto->CloneAndRegister(AnaHisto->hpsbe_1_0);

  LocalHisto.ht0_0_1 = AnaHisto->CloneAndRegister(AnaHisto->ht0_0_1);
  LocalHisto.ht0_0_2 = AnaHisto->CloneAndRegister(AnaHisto->ht0_0_2);
  LocalHisto.ht0_0_3 = AnaHisto->CloneAndRegister(AnaHisto->ht0_0_3);
  LocalHisto.ht0_0_4 = AnaHisto->CloneAndRegister(AnaHisto->ht0_0_4);
  for(int i = 0; i < 28; ++i)
    LocalHisto.ht0_1[i] = AnaHisto->CloneAndRegister(AnaHisto->ht0_1[i]);

  LocalHisto.hmdc_0_1 = AnaHisto->CloneAndRegister(AnaHisto->hmdc_0_1);
  LocalHisto.hmdc_0_2 = AnaHisto->CloneAndRegister(AnaHisto->hmdc_0_2);
  LocalHisto.hmdc_0_3 = AnaHisto->CloneAndRegister(AnaHisto->hmdc_0_3);
  LocalHisto.hmdc_0_4 = AnaHisto->CloneAndRegister(AnaHisto->hmdc_0_4);
  LocalHisto.hmdc_0_5 = AnaHisto->CloneAndRegister(AnaHisto->hmdc_0_5);
  LocalHisto.hmdc_0_6 = AnaHisto->CloneAndRegister(AnaHisto->hmdc_0_6);
  LocalHisto.hmdc_0_9 = AnaHisto->CloneAndRegister(AnaHisto->hmdc_0_9);
  for(int i = 0; i < 17; ++i)
    {
      LocalHisto.hmdc_1[i]   = AnaHisto->CloneAndRegister(AnaHisto->hmdc_1[i]);
      LocalHisto.hmdc_2[i]   = AnaHisto->CloneAndRegister(AnaHisto->hmdc_2[i]);
      LocalHisto.hmdc_2_2[i] = AnaHisto->CloneAndRegister(AnaHisto->hmdc_2_2[i]);
      LocalHisto.hmdc_3_2[i] = AnaHisto->CloneAndRegister(AnaHisto->hmdc_3_2[i]);
      LocalHisto.hmdc_3[i]   = AnaHisto->CloneAndRegister(AnaHisto->hmdc_3[i]);
      LocalHisto.hmdc_3_2[i] = AnaHisto->CloneAndRegister(AnaHisto->hmdc_3_2[i]);
      LocalHisto.hmdc_3_3[i] = AnaHisto->CloneAndRegister(AnaHisto->hmdc_3_3[i]);
    }

  LocalHisto.hmwdc_1_1 = AnaHisto->CloneAndRegister(AnaHisto->hmwdc_1_1);
  LocalHisto.hmwdc_1_2 = AnaHisto->CloneAndRegister(AnaHisto->hmwdc_1_2);
  LocalHisto.hmwdc_1_3 = AnaHisto->CloneAndRegister(AnaHisto->hmwdc_1_3);
  LocalHisto.hmwdc_1_4 = AnaHisto->CloneAndRegister(AnaHisto->hmwdc_1_4);
  LocalHisto.hmwdc_1_5 = AnaHisto->CloneAndRegister(AnaHisto->hmwdc_1_5);
  LocalHisto.hmwdc_1_6 = AnaHisto->CloneAndRegister(AnaHisto->hmwdc_1_6);

  LocalHisto.hs4sci_1_1 = AnaHisto->CloneAndRegister(AnaHisto->hs4sci_1_1);
  LocalHisto.hs4sci_1_2 = AnaHisto->CloneAndRegister(AnaHisto->hs4sci_1_2);
  LocalHisto.hs4sci_1_3 = AnaHisto->CloneAndRegister(AnaHisto->hs4sci_1_3);
  LocalHisto.hs4sci_1_4 = AnaHisto->CloneAndRegister(AnaHisto->hs4sci_1_4);
  LocalHisto.hs4sci_2_1 = AnaHisto->CloneAndRegister(AnaHisto->hs4sci_2_1);
  LocalHisto.hs4sci_2_2 = AnaHisto->CloneAndRegister(AnaHisto->hs4sci_2_2);
  LocalHisto.hs4sci_2_3 = AnaHisto->CloneAndRegister(AnaHisto->hs4sci_2_3);
  LocalHisto.hs4sci_2_4 = AnaHisto->CloneAndRegister(AnaHisto->hs4sci_2_4);

  LocalHisto.hopt_1_1 = AnaHisto->CloneAndRegister(AnaHisto->hopt_1_1);
  LocalHisto.hopt_1_2 = AnaHisto->CloneAndRegister(AnaHisto->hopt_1_2);
  LocalHisto.hopt_1_3 = AnaHisto->CloneAndRegister(AnaHisto->hopt_1_3);
  LocalHisto.hopt_1_4 = AnaHisto->CloneAndRegister(AnaHisto->hopt_1_4);
  LocalHisto.hopt_2_1 = AnaHisto->CloneAndRegister(AnaHisto->hopt_2_1);
  LocalHisto.hopt_2_2 = AnaHisto->CloneAndRegister(AnaHisto->hopt_2_2);
  LocalHisto.hopt_2_3 = AnaHisto->CloneAndRegister(AnaHisto->hopt_2_3);
  LocalHisto.hopt_2_4 = AnaHisto->CloneAndRegister(AnaHisto->hopt_2_4);

  LocalHisto.htrig_0 = AnaHisto->CloneAndRegister(AnaHisto->htrig_0);
  LocalHisto.htrig_1 = AnaHisto->CloneAndRegister(AnaHisto->htrig_1);
  LocalHisto.htrig_2 = AnaHisto->CloneAndRegister(AnaHisto->htrig_2);
  LocalHisto.htrig_3 = AnaHisto->CloneAndRegister(AnaHisto->htrig_3);
  LocalHisto.htrig_4 = AnaHisto->CloneAndRegister(AnaHisto->htrig_4);
}

ReturnRes::InfoM TBuildWASACalibrationLayerPlane::SoftExit(int return_build)
{
  if(return_build == -1)
    {
      att._logger->warn("!> Multiplicity > 2 on Start : event rejected");
      LocalHisto.h_Builderstats->Fill("start M>2", 1);
      return ReturnRes::MultiS2_Start;
    }
  else if(return_build == -2)
    {
      att._logger->warn("!> TDC Timing Start cut : event rejected");
      LocalHisto.h_Builderstats->Fill("start Timing cut", 1);
      return ReturnRes::StartTimingCut;
    }
  else if(return_build == -3)
    {
      att._logger->warn("!> Chamber Hit > 1000 : event rejected");
      LocalHisto.h_Builderstats->Fill("chamber hit>1000", 1);
      return ReturnRes::ChamberHitLimit;
    }
  else if(return_build == -4)
    {
      att._logger->warn("!> Trigger Cut: event rejected");
      LocalHisto.h_Builderstats->Fill("Trigger cut", 1);
      return ReturnRes::TriggerCut;
    }
  else if(return_build == -9)
    {
      att._logger->warn("!> No Beam : event rejected");
      LocalHisto.h_Builderstats->Fill("No Beam", 1);
      return ReturnRes::NoBeam;
    }
  else if(return_build != 0)
    {
      att._logger->warn("Error in Build Detector !");
      LocalHisto.h_Builderstats->Fill("Error", 1);
      return ReturnRes::BuildError;
    }
  LocalHisto.h_Builderstats->Fill("start Ok", 1);

  return ReturnRes::Fine;
}

#ifdef ROOT6
int TBuildWASACalibrationLayerPlane::Exec(const EventWASAUnpack& event, FullRecoEvent& RecoEvent,
                                          Ana_WasaEvent* OutTree)
#else
int TBuildWASACalibrationLayerPlane::Exec(const EventWASAUnpack& event, FullRecoEvent& RecoEvent,
                                          Ana_WasaEvent* OutTree)
#endif
{
  //att._logger->info("Start BuildWASACalibration!");

  // Sub classes of calibrations
  // Fiber : Input event.s2fiber

  // MDC : Input event.s2mdc

  // PSB - PSBE - PSFE - T0 : Input event.s2tq1 event.s2tq2

  // CsI

  // S4 - FRSTPC

  // After calibration :
  // creation of the genfit measurements
  // collection of additional information that does not go into genfit measurements (i.e ToF, dE, ion optics ...)

  // Save calibration into Ana_WasaEvent for possible restart after this point.

  //  ana tmp  //////////////////////////////////////////////////////////////////////////////////

  RecoEvent.ListHits.resize(G4Sol::SIZEOF_G4SOLDETTYPE);
  RecoEvent.ListHitsInfo.resize(G4Sol::SIZEOF_G4SOLDETTYPE);

  RecoEvent.InteractionPoint[0] = 0.;
  RecoEvent.InteractionPoint[1] = 0.;
  RecoEvent.InteractionPoint[2] = att.Target_PositionZ;

  RecoEvent.PrimVtxRecons.SetXYZ(0., 0., att.Target_PositionZ);

  // Trigger //////////////////
  bool trig[16];
  for(int i=0; i<32; ++i)
    {
      double val = event.s4tq->v775[i];
      if(val>0)
        LocalHisto.htrig_1->Fill(i, val);
    }
  for(int i=0; i<16; ++i)
  {
    trig[i] = event.s4tq->trig[i];
    int offset = 0;
    if(i>7)
      offset = 8;
    if(trig[i])
      {
        LocalHisto.htrig_0->Fill(i);
        double val = event.s4tq->v775[i+offset];
        LocalHisto.htrig_2->Fill(i+offset, val);
      }
  }

  if( par->trig_main  && !trig[3] ) return -4;
  if( par->trig_clock && !trig[5] ) return -4;
  if( par->trig_t0    && !(trig[7] && trig[9]) ) return -4;
  if( par->trig_sc41  && !trig[1] ) return -4;

  for(int i=0; i<16; ++i)
    {
      int offset = 0;
      if(i>7)
        offset = 8;
      if(trig[i])
        {
          LocalHisto.htrig_3->Fill(i);
          double val = event.s4tq->v775[i+offset];
          LocalHisto.htrig_4->Fill(i+offset, val);
        }
    }


  // T0 ana //////////////////////////////////
  std::vector<T0HitAna*> T0HitCont;
  for(int i = 0; i < 28; ++i)
    {
      bool flag_u = false;
      bool flag_d = false;
      int i_u     = -1;
      int i_d     = -1;
      for(int k = 0; k < event.s2tq2->nhit_t0[i][0]; ++k)
        {
          if(par->t0_tcut_min < event.s2tq2->tdc_t0[i][0][k] && event.s2tq2->tdc_t0[i][0][k] < par->t0_tcut_max)
            {
              flag_u = true;
              i_u    = k;
              break;
            }
        }
      for(int k = 0; k < event.s2tq2->nhit_t0[i][1]; ++k)
        {
          if(par->t0_tcut_min < event.s2tq2->tdc_t0[i][1][k] && event.s2tq2->tdc_t0[i][1][k] < par->t0_tcut_max)
            {
              flag_d = true;
              i_d    = k;
              break;
            }
        }
      if(flag_u && flag_d)
        {
          int t_u           = event.s2tq2->tdc_t0[i][0][i_u];
          int t_d           = event.s2tq2->tdc_t0[i][1][i_d];
          int q_u           = event.s2tq2->qdc_t0[i][0];
          int q_d           = event.s2tq2->qdc_t0[i][1];
          T0HitAna* hit_ana = new T0HitAna(i, t_u, t_d, q_u, q_d, par.get());
          T0HitCont.emplace_back(hit_ana);
        }
    }

  T0HitAna* hit_t0_main = nullptr;
  double time_ref = 9999.;

  for(int i = 0; i < (int)T0HitCont.size(); ++i)
    {
      LocalHisto.ht0_0_1->Fill(T0HitCont[i]->GetTU());
      LocalHisto.ht0_0_2->Fill(T0HitCont[i]->GetTD());
      LocalHisto.ht0_0_4->Fill(T0HitCont[i]->GetSeg());
      LocalHisto.ht0_1[T0HitCont[i]->GetSeg()]->Fill(T0HitCont[i]->GetTime());

      double time_buf = std::fabs(T0HitCont[i]->GetTime());
      if(time_buf < time_ref && time_buf < par->cut_t0_time)
      {
        hit_t0_main = T0HitCont[i];
        time_ref = time_buf;
      }
    }

    //std::cout << time_ref << "\n";
  LocalHisto.ht0_0_3->Fill(T0HitCont.size());




  // PSB ana //////////////////////////////////
  std::vector<PSBHitAna*> PSBHitCont;
  for(int i = 0; i < 46; ++i)
    {
      bool flag_u = false;
      bool flag_d = false;
      int i_u     = -1;
      int i_d     = -1;
      for(int k = 0; k < event.s2tq1->nhit_psb[i][0]; ++k)
        {
          if(par->psb_tcut_min < event.s2tq1->tdc_psb[i][0][k] && event.s2tq1->tdc_psb[i][0][k] < par->psb_tcut_max)
            {
              flag_u = true;
              i_u    = k;
              break;
            }
        }
      for(int k = 0; k < event.s2tq1->nhit_psb[i][1]; ++k)
        {
          if(par->psb_tcut_min < event.s2tq1->tdc_psb[i][1][k] && event.s2tq1->tdc_psb[i][1][k] < par->psb_tcut_max)
            {
              flag_d = true;
              i_d    = k;
              break;
            }
        }
      if(flag_u && flag_d)
        {
          double t_t0 = 0.;
          if(hit_t0_main) t_t0 = hit_t0_main->GetTime();

          int t_u            = event.s2tq1->tdc_psb[i][0][i_u];
          int t_d            = event.s2tq1->tdc_psb[i][1][i_d];
          int q_u            = event.s2tq1->qdc_psb[i][0];
          int q_d            = event.s2tq1->qdc_psb[i][1];
          PSBHitAna* hit_ana = new PSBHitAna(i, t_u, t_d, q_u, q_d, par.get(), t_t0);
          PSBHitCont.emplace_back(hit_ana);
        }
    }

  for(const auto& hit_PSB : PSBHitCont)
    {
      int hitID      = hit_PSB->GetSeg();
      TGeoMatrix* g1 = gGeoManager->GetVolume("INNER")->GetNode(offsetGeoNameID_PSCE + hitID)->GetMatrix(); // PSCE
      TGeoMatrix* g2 = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix();           // INNER
      TGeoMatrix* g3 = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();           // MFLD
      TGeoHMatrix H1(*g1), H2(*g2), H3(*g3);
      TGeoHMatrix H = H2 * H1;
      H             = H3 * H;
#ifdef DEBUG_BUILD2
      H.Print();
#endif
      TGeoHMatrix Hsf("Hsf"); // PSCE inner surface
      Hsf.SetDz(-0.4);
      H             = H * Hsf;
      double* shift = H.GetTranslation();

      TVector3 o(shift[0], shift[1], shift[2]), phidir(shift[0], shift[1], 0), zdir(0., 0., 1.);
      phidir     = phidir.Unit();
      TVector3 u = zdir.Cross(phidir);
      TVector3 v = zdir;
      genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));

      TVectorD hitCoords(2);
      hitCoords(0) = 0.;
      hitCoords(1) = hit_PSB->GetZ() * 0.1; // mm to cm
      TMatrixDSym hitCov(2);
      //hitCov(0, 0)     = TMath::Sq(par->psb_res_phi * 0.1);//0.1; // to be adjusted resolution_psce * resolution_psce;
      //hitCov(1, 1)     = TMath::Sq(par->psb_res_z   * 0.1);//0.1; // to be adjusted resolution_psce_z * resolution_psce_z;
      hitCov(0, 0)     = std::pow(par->psb_res_phi * 0.1, 2.);//0.1; // to be adjusted resolution_psce * resolution_psce;
      hitCov(1, 1)     = std::pow(par->psb_res_z   * 0.1, 2.);//0.1; // to be adjusted resolution_psce_z * resolution_psce_z;
      auto measurement = std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, G4Sol::PSCE, hitID, nullptr);
      dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

      RecoEvent.ListHits[G4Sol::PSCE].emplace_back(measurement.release());

      MeasurementInfo measinfo(-9999., hit_PSB->GetTime(), -9999.); //CHECK include dE
      RecoEvent.ListHitsInfo[G4Sol::PSCE].emplace_back(measinfo);

      LocalHisto.h76->Fill(hitID, atan2(shift[1], shift[0]));
    }

  for(int i = 0; i < (int)PSBHitCont.size(); ++i)
    {
      LocalHisto.hpsb_0_1->Fill(PSBHitCont[i]->GetTU());
      LocalHisto.hpsb_0_3->Fill(PSBHitCont[i]->GetPhi());
      LocalHisto.hpsb_0_4->Fill(PSBHitCont[i]->GetSeg(), PSBHitCont[i]->GetPhi());
      LocalHisto.hpsb_1_1->Fill(PSBHitCont[i]->GetZ());
      LocalHisto.hpsb_2[PSBHitCont[i]->GetSeg()]->Fill(PSBHitCont[i]->GetTime());
      LocalHisto.hpsb_4[PSBHitCont[i]->GetSeg()]->Fill(PSBHitCont[i]->GetZ());
    }
  LocalHisto.hpsb_0_2->Fill(PSBHitCont.size());

  // PSFE ana //////////////////////////////////
  std::vector<PSFEHitAna*> PSFEHitCont;
  for(int i = 0; i < 44; ++i)
    for(int j = 0; j < event.s2tq2->nhit_psfe[i]; ++j)
      {
        if(par->psfe_tcut_min < event.s2tq2->tdc_psfe[i][j] && event.s2tq2->tdc_psfe[i][j] < par->psfe_tcut_max)
          {
            double t_t0 = 0.;
            if(hit_t0_main) t_t0 = hit_t0_main->GetTime();

            int t_buf           = event.s2tq2->tdc_psfe[i][j];
            int q_buf           = event.s2tq2->qdc_psfe[i];
            PSFEHitAna* hit_ana = new PSFEHitAna(i, t_buf, q_buf, par.get(), t_t0);
            PSFEHitCont.emplace_back(hit_ana);
          }
      }

  for(const auto& hit_PSFE : PSFEHitCont)
    {
      int hitID = hit_PSFE->GetSeg();
      TGeoMatrix* g1   = gGeoManager->GetVolume("PSB")->GetNode(hitID)->GetMatrix();    // PSFE   -> In Geometry file PSBE/PSFE flipped
      TGeoMatrix* g1_1 = gGeoManager->GetVolume("MFLD")->GetNode("PSB_1")->GetMatrix(); // PSFE box
      TGeoMatrix* g2   = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix();       // INNER
      TGeoMatrix* g3   = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();       // MFLD
      TGeoHMatrix H1(*g1), H1_1(*g1_1), H2(*g2), H3(*g3);
      TGeoHMatrix H = H1_1 * H1;
      H             = H2 * H;

      H = H3 * H;
#ifdef DEBUG_BUILD2
      H.Print();
#endif

      double* shift     = H.GetTranslation();
      double* local_rot = H.GetRotationMatrix();

      TVector3 v(local_rot[0], local_rot[3], local_rot[6]);
      // v is at the left border of the bar -> rotate 3.75 degree to be at the center of the bar
      v.RotateZ(-3.75 * TMath::DegToRad());
      v = v.Unit();

      TVector3 u(v.Y(), -v.X(), 0.);
      TVector3 o(shift[0], shift[1], shift[2]);
      genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));

      // TVectorD hitCoords(2);
      // hitCoords(0) = gRandom->Uniform(-3.75, 3.75); // to be adjusted // phi ! be aware ! not u-dim
      // hitCoords(1) = gRandom->Uniform(6., 22.);     // to be adjusted   // r -> v dir

      // TMatrixDSym hitCov(2);
      // hitCov(0, 0)     = TMath::Sq(2 * hitCoords(1) * TMath::Sin(3.75 * TMath::DegToRad())) / 12.;
      // hitCov(1, 1)     = TMath::Sq(22. - 6.) / 12.;

      TVectorD hitCoords(1);
      hitCoords(0) = 0.;
      TMatrixDSym hitCov(1);
      //hitCov(0, 0) = TMath::Sq(par->psfe_res_phi * 0.1); // mm to cm
      hitCov(0, 0) = std::pow(par->psfe_res_phi * 0.1, 2.); // mm to cm

      auto measurement = std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, G4Sol::PSFE, hitID, nullptr);
      dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

      RecoEvent.ListHits[G4Sol::PSFE].emplace_back(measurement.release());

      MeasurementInfo measinfo(-9999., hit_PSFE->GetTime(), -9999.);
      RecoEvent.ListHitsInfo[G4Sol::PSFE].emplace_back(measinfo);
    }

  for(int i = 0; i < (int)PSFEHitCont.size(); ++i)
    {
      LocalHisto.hpsfe_0_1->Fill(PSFEHitCont[i]->GetT());
      LocalHisto.hpsfe_0_3->Fill(PSFEHitCont[i]->GetPhi());
      LocalHisto.hpsfe_0_4->Fill(PSFEHitCont[i]->GetSeg(), PSFEHitCont[i]->GetPhi());
    }
  LocalHisto.hpsfe_0_2->Fill(PSFEHitCont.size());

  // PSBE ana //////////////////////////////////
  std::vector<PSBEHitAna*> PSBEHitCont;
  for(int i = 0; i < 38; ++i)
    {
      for(int j = 0; j < event.s2tq2->nhit_psbe[i]; ++j)
        {
          LocalHisto.hpsbe_1_0->Fill(i, event.s2tq2->tdc_psbe[i][j]);
          if(par->psbe_tcut_min < event.s2tq2->tdc_psbe[i][j] && event.s2tq2->tdc_psbe[i][j] < par->psbe_tcut_max)
            {
              int t_buf           = event.s2tq2->tdc_psbe[i][j];
              int q_buf           = event.s2tq2->qdc_psbe[i];
              PSBEHitAna* hit_ana = new PSBEHitAna(i, t_buf, q_buf, par.get());
              PSBEHitCont.emplace_back(hit_ana);
            }
        }
    }

  for(const auto& hit_PSBE : PSBEHitCont)
    {
      int hitID = hit_PSBE->GetSeg();

      TGeoMatrix* g1   = gGeoManager->GetVolume("PSF")->GetNode(hitID)->GetMatrix();    // PSBE    -> In geometry file, PSFE/PSBE are flipped 
      TGeoMatrix* g1_1 = gGeoManager->GetVolume("MFLD")->GetNode("PSF_1")->GetMatrix(); // PSBE box
      TGeoMatrix* g2   = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix();       // INNNER
      TGeoMatrix* g3   = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();       // MFLD
      TGeoHMatrix H1(*g1), H1_1(*g1_1), H2(*g2), H3(*g3);
      TGeoHMatrix H = H1_1 * H1;
      H             = H2 * H;

      H = H3 * H;
#ifdef DEBUG_BUILD2
      H.Print();
#endif

      double* shift     = H.GetTranslation();
      double* local_rot = H.GetRotationMatrix();

      TVector3 v(local_rot[0], local_rot[3], local_rot[6]);
      // v is at the left border of the bar -> rotate 3.75 degree to be at the center of the bar
      v.RotateZ(-3.75 * TMath::DegToRad());
      v = v.Unit();

      TVector3 u(v.Y(), -v.X(), 0.);
      TVector3 o(shift[0], shift[1], shift[2]);
      genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));

      // TVectorD hitCoords(2);
      // hitCoords(0) = gRandom->Uniform(-3.75, 3.75); // to be adjusted // phi ! be aware ! not u-dim
      // hitCoords(1) = gRandom->Uniform(6., 22.);     // to be adjusted   // r -> v dir

      // TMatrixDSym hitCov(2);
      // hitCov(0, 0)     = TMath::Sq(2 * hitCoords(1) * TMath::Sin(3.75 * TMath::DegToRad())) / 12.;
      // hitCov(1, 1)     = TMath::Sq(22. - 6.) / 12.;
      TVectorD hitCoords(1);
      hitCoords(0) = 0.;
      TMatrixDSym hitCov(1);
      hitCov(0, 0) = pow(par->psfe_res_phi * mm2cm, 2); // must be corrected

      auto measurement = std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, G4Sol::PSBE, hitID, nullptr);
      dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

      RecoEvent.ListHits[G4Sol::PSBE].emplace_back(measurement.release());
    }

  for(int i = 0; i < (int)PSBEHitCont.size(); ++i)
    {
      LocalHisto.hpsbe_0_1->Fill(PSBEHitCont[i]->GetT());
      LocalHisto.hpsbe_0_3->Fill(PSBEHitCont[i]->GetPhi());
      LocalHisto.hpsbe_0_4->Fill(PSBEHitCont[i]->GetSeg(), PSBEHitCont[i]->GetPhi());
    }
  LocalHisto.hpsbe_0_2->Fill(PSBEHitCont.size());


  // MDC ana //////////////////////////////////
  std::vector<std::vector<MDCHitAna*> > MDCHitCont;
  for(int i = 0; i < 17; ++i)
    {
      std::vector<MDCHitAna*> buf_v;
      MDCHitCont.emplace_back(buf_v);
    }
  for(int i = 0; i < (int)event.s2mdc->mdchit.size(); ++i)
    {
      double t_t0 = 0.;
      if(hit_t0_main) t_t0 = hit_t0_main->GetTime();

      MDCHitAna* hit_ana = new MDCHitAna(event.s2mdc->mdchit[i], par.get(), event.s2mdc->tref[0], t_t0);
      if(!hit_ana->IsValid())
        {
          delete hit_ana;
          continue;
        }
      MDCHitCont[hit_ana->GetLay()].emplace_back(hit_ana);
      // hit_ana->Print();
    }

  for(int MDCLayerID = 0; MDCLayerID < 17; ++MDCLayerID)
    for(const auto& hit_MDC : MDCHitCont[MDCLayerID])
      {
        int hitID = hit_MDC->GetWir();

        TGeoMatrix* g1 =
            gGeoManager->GetVolume("INNER")->GetNode(MDCLayerID + offsetGeoNameID_MDC)->GetVolume()->GetNode(hitID)->GetMatrix(); // ME,
                                                                                                                // MG
        TGeoShape* tempShape = gGeoManager->GetVolume("INNER")->GetNode(MDCLayerID + offsetGeoNameID_MDC)->GetVolume()->GetShape();
        TGeoMatrix* g2       = gGeoManager->GetVolume("INNER")->GetNode(MDCLayerID + offsetGeoNameID_MDC)->GetMatrix(); // MD
        TGeoMatrix* g3       = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix();               // INNER
        TGeoMatrix* g4       = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();               // MFLD
        TGeoHMatrix H1(*g1), H2(*g2), H3(*g3), H4(*g4);
        TGeoHMatrix H = H2 * H1;
        H             = H3 * H;
        H             = H4 * H;
        double* shift = H.GetTranslation();
        TGeoHMatrix w1("w1");
        TGeoHMatrix w2("w2");
        Double_t minZ, maxZ;
        tempShape->GetAxisRange(3, minZ, maxZ); // length of wire
        w1.SetDz(minZ);
        w2.SetDz(maxZ);
        TGeoHMatrix Hw1 = H * w1;
        TGeoHMatrix Hw2 = H * w2;
#ifdef DEBUG_BUILD2
        H.Print();
        Hw1.Print();
        Hw2.Print();
#endif
        double* edge1 = Hw1.GetTranslation();
        double* edge2 = Hw2.GetTranslation();

        double dl = hit_MDC->GetDl() * 0.1 ; // mm to cm

        double dlmax = 0;
        switch(MDCLayerID + 1)
          {
          case 1:
          case 2:
          case 3:
          case 4:
          case 5:
            dlmax = 0.2;
            break;
          case 6:
          case 7:
          case 8:
          case 9:
          case 10:
          case 11:
            dlmax = 0.3;
            break;
          case 12:
          case 13:
          case 14:
          case 15:
          case 16:
          case 17:
            dlmax = 0.4;
            break;
          default:
            att._logger->warn("Error in WireMeasurement !");
            break;
          }

        TVectorD hitCoords(7);
        hitCoords(0) = edge1[0];
        hitCoords(1) = edge1[1];
        hitCoords(2) = edge1[2];
        hitCoords(3) = edge2[0];
        hitCoords(4) = edge2[1];
        hitCoords(5) = edge2[2];
        hitCoords(6) = dl;
        if(edge1[2] > edge2[2])
          {
            hitCoords(0) = edge2[0];
            hitCoords(1) = edge2[1];
            hitCoords(2) = edge2[2];
            hitCoords(3) = edge1[0];
            hitCoords(4) = edge1[1];
            hitCoords(5) = edge1[2];
          }

        //printf("MDC Edge1z: %f , Edge2: %f \n", edge1[2], edge2[2]);
        TMatrixDSym hitCov(7);
        //hitCov(6, 6) = TMath::Sq(hit_MDC->GetRes()*0.1); // to be adjucted resolution_dl * resolution_dl;
        hitCov(6, 6) = std::pow(hit_MDC->GetRes()*0.1, 2.); // to be adjucted resolution_dl * resolution_dl;
        auto measurement =
            std::make_unique<genfit::WireMeasurement>(hitCoords, hitCov, G4Sol::MG01 + MDCLayerID, hitID, nullptr);
        dynamic_cast<genfit::WireMeasurement*>(measurement.get())->setLeftRightResolution(1);
        dynamic_cast<genfit::WireMeasurement*>(measurement.get())->setMaxDistance(dlmax);

        RecoEvent.ListHits[G4Sol::MG01 + MDCLayerID].emplace_back(measurement.release());

        MeasurementInfo measinfo(hit_MDC->GetTot(), -9999., -9999.);
        RecoEvent.ListHitsInfo[G4Sol::MG01 + MDCLayerID].emplace_back(measinfo);
      }

  for(int i = 0; i < 17; ++i)
    {
      for(int j = 0; j < (int)MDCHitCont[i].size(); ++j)
        {
          LocalHisto.hmdc_0_1->Fill(MDCHitCont[i][j]->GetTL());
          LocalHisto.hmdc_0_2->Fill(MDCHitCont[i][j]->GetTL(), MDCHitCont[i][j]->GetTot());
          LocalHisto.hmdc_0_3->Fill(MDCHitCont[i][j]->GetR());
          LocalHisto.hmdc_0_4->Fill(MDCHitCont[i][j]->GetPhi());
          LocalHisto.hmdc_1[i]->Fill(MDCHitCont[i][j]->GetTL());
          LocalHisto.hmdc_2[i]->Fill(MDCHitCont[i][j]->GetDt());
          LocalHisto.hmdc_3[i]->Fill(MDCHitCont[i][j]->GetDl());
          if(i == 16)
            {
              LocalHisto.hmdc_0_5->Fill(MDCHitCont[i][j]->GetTL());
              LocalHisto.hmdc_0_6->Fill(MDCHitCont[i][j]->GetTL(), MDCHitCont[i][j]->GetTot());
              LocalHisto.hmdc_0_9->Fill(MDCHitCont[i][j]->GetWir(), MDCHitCont[i][j]->GetPhi());
            }
        }
    }

  // Fiber ana //////////////////////////////////
  std::vector<std::vector<std::vector<FiberHitAna*> > > FiberHitCont;
  std::vector<std::vector<std::vector<FiberHitAna*> > > FiberHitClCont;

  for(int i = 0; i < 7; ++i)
    {
      std::vector<FiberHitAna*> buf_v;
      std::vector<std::vector<FiberHitAna*> > buf_vv;
      for(int j = 0; j < 3; ++j)
        {
          buf_vv.emplace_back(buf_v);
        }
      FiberHitCont.emplace_back(buf_vv);
      FiberHitClCont.emplace_back(buf_vv);
    }

  int t_r = event.s2fiber->tref[3];

  for(int i = 0; i < (int)event.s2fiber->fiberhit.size(); ++i)
    {
      double t_t0 = 0.;
      if(hit_t0_main) t_t0 = hit_t0_main->GetTime();
      FiberHitAna* hit_ana = new FiberHitAna(event.s2fiber->fiberhit[i], par.get(), t_r, t_t0);

      if(!hit_ana->IsValid())
        {
          delete hit_ana;
          continue;
        }
      FiberHitCont[hit_ana->GetDet()][hit_ana->GetLay()].emplace_back(hit_ana);
    }

  FiberAnalyzer* fiberana = new FiberAnalyzer();
  FiberHitClCont          = fiberana->Clusterize(FiberHitCont);

  for(int i=0; i<7; ++i){
    for(int j=0; j<3; ++j){
      //std::cout << "\nsize before : " << FiberHitCont[i][j].size() << std::endl;
      //std::cout << "size after  : " << FiberHitClCont[i][j].size() << std::endl;
      for(int k=0; k<(int)FiberHitClCont[i][j].size(); ++k){
        //FiberHitClCont[i][j][k]->Print();
        FiberHitAna *hit = FiberHitClCont[i][j][k];
        LocalHisto.hfiber_13_0[i][j]->Fill(hit->GetTL());
        LocalHisto.hfiber_13_1[i][j]->Fill(hit->GetTime());
        double t_buf = 0;
        if(hit_t0_main) t_buf = hit_t0_main->GetTime();
        LocalHisto.hfiber_13_2[i][j]->Fill(hit->GetTime() + t_buf);
        LocalHisto.hfiber_13_3[i][j]->Fill(hit->GetFib(), hit->GetTL());
        LocalHisto.hfiber_13_4[i][j]->Fill(hit->GetFib(), hit->GetTime());
      }
    }
  }

  for(int i=2; i<5; ++i){
    for(int j=0; j<3; ++j){
      for(int k=0; k<(int)FiberHitClCont[i][j].size(); ++k){
        LocalHisto.h51[i-2][j][0]->Fill(FiberHitClCont[i][j][k]->GetFib(), FiberHitClCont[i][j][k]->GetTOT());
      }
    }
  }

    //temp
    std::vector< std::vector<FiberHitXUV*> > FiberXUVCont = fiberana->FindHit(FiberHitClCont, par.get());
    std::map< std::string, std::vector<FiberTrackAna*> > FiberTrackCont;

    // MFT12
    int nt_mft12 = 0;
    int nt_mft12_xuv = 0;
    if(FiberXUVCont[3].size()>0 && FiberXUVCont[4].size()>0 && !par->flag_mft12_allcombi){
      std::vector<FiberTrackAna*> buf_track;
      for(int i=0; i<(int)FiberXUVCont[3].size(); ++i){
        for(int j=0; j<(int)FiberXUVCont[4].size(); ++j){
          std::vector<FiberHitXUV*>   buf_xuv;
          buf_xuv.emplace_back(FiberXUVCont[3][i]);
          buf_xuv.emplace_back(FiberXUVCont[4][j]);
          FiberTrackAna *track = new FiberTrackAna(buf_xuv, par.get());
          if(par->flag_mft12_posang){
            track->CorrectMFT(par.get());
            double buf_x = track->GetXmft();
            double buf_y = track->GetYmft();
            double buf_a = track->GetA();
            double buf_b = track->GetB();
            if( fabs(buf_x * 0.003 - buf_a)>0.3 || fabs(buf_y * 0.003 - buf_b)>0.3 ){ delete track; continue;}
          }
          if(!par->flag_mft12_combi)                      buf_track.emplace_back(track);
          else if(track->GetChi2() < par->cut_chi2_mft12) buf_track.emplace_back(track);
        }
      }
/*
      if(par->flag_mft12_xuv_psb){
        for(int i=0; i<(int)buf_track.size(); ++i){
          for(int j=0; j<(int)PSBHitCont.size(); ++j){
            double a_fiber = buf_track[i]->GetA();
            double b_fiber = buf_track[i]->GetB();
            double x_fiber = buf_track[i]->GetX();
            double y_fiber = buf_track[i]->GetY();
            double phi_psb   = PSBHitCont[j]->GetPhi() + par->psb_rot_z*Deg2Rad;
            double r_psb     = PSBHitCont[j]->GetR();
            double z_psb     = PSBHitCont[j]->GetZ();

            double par_a = pow(a_fiber, 2) + pow(b_fiber, 2);
            double par_b = a_fiber * (x_fiber - par->psb_pos_x)+ b_fiber * (y_fiber - par->psb_pos_y);
            double par_c = pow(x_fiber - par->psb_pos_x, 2) + pow(y_fiber - par->psb_pos_y, 2) - pow(r_psb, 2);
            double z_fiber  = (-par_b + sqrt( pow(par_b,2) - par_a * par_c)) / par_a;

            double fiber_x_buf = x_fiber + a_fiber * z_fiber - par->psb_pos_x;
            double fiber_y_buf = y_fiber + b_fiber * z_fiber - par->psb_pos_y;
            double phi_fiber = atan2(fiber_y_buf, fiber_x_buf);
            if( fabs(fiberana->CalcPhiDif(phi_psb, phi_fiber)) < par->cut_psb_phi ){
              if(fabs( (z_fiber - par->psb_pos_z) - z_psb)<par->cut_psb_z){
                if(buf_track[i]->IsFlagPSB()){
                  double phi_psb_buf = buf_track[i]->GetPSBHit()->GetPhi();
                  if( fabs( fiberana->CalcPhiDif(phi_psb, phi_fiber) ) < fabs( fiberana->CalcPhiDif(phi_psb_buf, phi_fiber) ) ){
                    buf_track[i]->SetSegPSB(PSBHitCont[j]->GetSeg());
                    buf_track[i]->SetPSBHit(PSBHitCont[j]);
                  }
                }
                else{
                  buf_track[i]->SetFlagPSB();
                  buf_track[i]->SetSegPSB(PSBHitCont[j]->GetSeg());
                  buf_track[i]->SetPSBHit(PSBHitCont[j]);
                }
              }
            }
          }
        }
        int num_buf2 = buf_track.size();
        for(int i=num_buf2-1; i>=0; --i){
          if(!buf_track[i]->IsFlagPSB()){
            delete buf_track[i];
            buf_track.erase(buf_track.begin() + i);
          }
        }
      }
*/
      if(par->flag_dup_mft12_xuv && (int)buf_track.size()>0) buf_track = fiberana->DeleteDup(buf_track);
      FiberTrackCont["mft12"] = buf_track;

      nt_mft12     = FiberTrackCont["mft12"].size();
      nt_mft12_xuv = FiberTrackCont["mft12"].size();
      for(auto v: FiberTrackCont["mft12"]){
        for(int i=0; i<6; ++i){
          if(par->flag_dup_mft12_xuv) v->GetContHit().at(i)->SetUsed();
        }
        //std::cout << "-- 1st cor"  << std::endl;
        v->CorrectMFT(par.get());
        //std::cout << "-- 2nd cor"  << std::endl;
        //v->CorrectMFT(par);
        //std::cout << "-- 3rd cor"  << std::endl;
        //v->CorrectMFT(par);
        v->SetPosL();
        //if(par->flag_mftcor_xy) v->CorrectMFTXY(par.get());
      }
    }

    int num_combi_mft12 = 1;
    for(int i=3; i<5; ++i){
      for(int j=0; j<3; ++j){
        num_combi_mft12 *= ( (int)FiberHitClCont[i][j].size() + 1 );
      }
    }
    LocalHisto.hfiber_1_5->Fill((double)num_combi_mft12*1e-3);
    LocalHisto.hfiber_1_6->Fill((double)num_combi_mft12*1e-3);
    LocalHisto.hfiber_1_7->Fill((double)num_combi_mft12*1e-6);
    //std::cout << "num_combi_mft12 : " << num_combi_mft12 << std::endl;
    if(par->flag_debug) std::cout << "- num_combi_mft12 : " << num_combi_mft12 << std::endl;

    if(par->flag_mft12_combi && num_combi_mft12<par->cut_mft12_combi){
      LocalHisto.hfiber_1_1->Fill(num_combi_mft12);
      LocalHisto.hfiber_1_2->Fill((double)num_combi_mft12*1e-3);
      LocalHisto.hfiber_1_3->Fill((double)num_combi_mft12*1e-3);
      LocalHisto.hfiber_1_4->Fill((double)num_combi_mft12*1e-6);
      //std::vector< std::vector<int> > hit_combi;
      std::vector<FiberTrackAna*> buf_track;
      for(int a=-1; a<(int)FiberHitClCont[3][0].size(); ++a){
        for(int b=-1; b<(int)FiberHitClCont[3][1].size(); ++b){
          for(int c=-1; c<(int)FiberHitClCont[3][2].size(); ++c){
            for(int d=-1; d<(int)FiberHitClCont[4][0].size(); ++d){
              for(int e=-1; e<(int)FiberHitClCont[4][1].size(); ++e){
                for(int f=-1; f<(int)FiberHitClCont[4][2].size(); ++f){
                  std::vector<FiberHitAna*> buf_hit;
                  int count = 0;
                  if(a>-1 && !FiberHitClCont[3][0][a]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[3][0][a]); count++;}
                  if(b>-1 && !FiberHitClCont[3][1][b]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[3][1][b]); count++;}
                  if(c>-1 && !FiberHitClCont[3][2][c]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[3][2][c]); count++;}
                  if(d>-1 && !FiberHitClCont[4][0][d]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[4][0][d]); count++;}
                  if(e>-1 && !FiberHitClCont[4][1][e]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[4][1][e]); count++;}
                  if(f>-1 && !FiberHitClCont[4][2][f]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[4][2][f]); count++;}
                  if(count<4) continue;
                  //if(count!=4) continue;
                  FiberTrackAna *track = new FiberTrackAna(buf_hit, par.get());
                  track->SetFlagCombi();
                  if(par->flag_mft12_posang){
                    track->CorrectMFTCombi(par.get());
                    //track->CorrectMFTCombi(par);
                    double buf_x = track->GetXmft();
                    double buf_y = track->GetYmft();
                    double buf_a = track->GetA();
                    double buf_b = track->GetB();
                    if( fabs(buf_x * 0.003 - buf_a)>0.3 || fabs(buf_y * 0.003 - buf_b)>0.3 ){ delete track; continue;}
                    //if( pow(buf_x, 2.) + pow(buf_y, 2.) > pow(100., 2.)){ delete track; continue;}
                  }
                  switch(track->GetNlayer()){
                    case 4: buf_track.emplace_back(track); break;
                    case 5: if(track->GetChi2() < par->cut_chi2_mft12) buf_track.emplace_back(track); else delete track; break;
                    case 6: if(track->GetChi2() < par->cut_chi2_mft12) buf_track.emplace_back(track); else delete track; break;
                    default: break;
                  }
                }
              }
            }
          }
        }
      }

      LocalHisto.hfiber_2_1_1->Fill(buf_track.size());
      LocalHisto.hfiber_2_1_2->Fill((double)buf_track.size()/1000.);

      for(int i=0; i<(int)buf_track.size(); ++i){
        for(int j=0; j<(int)PSBHitCont.size(); ++j){
          double a_fiber = buf_track[i]->GetA();
          double b_fiber = buf_track[i]->GetB();
          double x_fiber = buf_track[i]->GetX();
          double y_fiber = buf_track[i]->GetY();
          //double phi_fiber = atan2(b_fiber, a_fiber);
          double phi_psb   = PSBHitCont[j]->GetPhi() + par->psb_rot_z*Deg2Rad;
          double r_psb     = PSBHitCont[j]->GetR();
          double z_psb     = PSBHitCont[j]->GetZ();

          double par_a = pow(a_fiber, 2) + pow(b_fiber, 2);
          double par_b = a_fiber * (x_fiber - par->psb_pos_x)+ b_fiber * (y_fiber - par->psb_pos_y);
          double par_c = pow(x_fiber - par->psb_pos_x, 2) + pow(y_fiber - par->psb_pos_y, 2) - pow(r_psb, 2);
          double z_fiber  = (-par_b + sqrt( pow(par_b,2) - par_a * par_c)) / par_a;

          double fiber_x_buf = x_fiber + a_fiber * z_fiber - par->psb_pos_x;
          double fiber_y_buf = y_fiber + b_fiber * z_fiber - par->psb_pos_y;
          double phi_fiber = atan2(fiber_y_buf, fiber_x_buf);
          LocalHisto.hfiber_6_1->Fill( fiberana->CalcPhiDif(phi_psb, phi_fiber) );
          LocalHisto.hfiber_6_2->Fill( (z_fiber - par->psb_pos_z) - z_psb );
          LocalHisto.hfiber_6_3->Fill( phi_fiber * Rad2Deg, phi_psb * Rad2Deg );
          LocalHisto.hfiber_6_4->Fill( (z_fiber - par->psb_pos_z),  z_psb );
          if( fabs(fiberana->CalcPhiDif(phi_psb, phi_fiber)) < par->cut_psb_phi ){
            if(fabs( (z_fiber - par->psb_pos_z) - z_psb)<par->cut_psb_z){
              if(buf_track[i]->IsFlagPSB()){
                double phi_psb_buf = buf_track[i]->GetPSBHit()->GetPhi();
                if( fabs( fiberana->CalcPhiDif(phi_psb, phi_fiber) ) < fabs( fiberana->CalcPhiDif(phi_psb_buf, phi_fiber) ) ){
                  buf_track[i]->SetSegPSB(PSBHitCont[j]->GetSeg());
                  buf_track[i]->SetPSBHit(PSBHitCont[j]);
                }
              }
              else{
                buf_track[i]->SetFlagPSB();
                buf_track[i]->SetSegPSB(PSBHitCont[j]->GetSeg());
                buf_track[i]->SetPSBHit(PSBHitCont[j]);
              }
            }
          }

        }
      }

      int num_buf = buf_track.size();
      for(int i=num_buf-1; i>=0; --i){
        if(!buf_track[i]->IsFlagPSB()){
          delete buf_track[i];
          buf_track.erase(buf_track.begin() + i);
        }
      }

      LocalHisto.hfiber_2_3->Fill(buf_track.size());

      if(par->flag_dup_mft12_combi){
        std::vector<FiberTrackAna*> buf_track_tmp = fiberana->DeleteDupCombi(buf_track);
        buf_track = buf_track_tmp;
      }

      for(auto v: buf_track){
        for(auto v2: v->GetContHit()){ v2->SetUsed(); }
        v->CorrectMFTCombi(par.get());
        //v->CorrectMFTCombi(par);
        //v->CorrectMFTCombi(par);
        v->SetPosL();
        v->DelFlagPSB();
        //if(par->flag_mftcor_xy) v->CorrectMFTXY(par);
      }

      LocalHisto.hfiber_2_2_1->Fill(buf_track.size());
      LocalHisto.hfiber_2_2_2->Fill((double)buf_track.size()/1000.);

      for(auto v : buf_track){
        FiberTrackCont["mft12"].emplace_back(v);
      }

      FiberTrackCont["mft12"] = fiberana->DeleteSame(FiberTrackCont["mft12"]);
      if(par->flag_debug) std::cout << "- before Inclusive : " << FiberTrackCont["mft12"].size() << std::endl;
      if(par->flag_mft12_inclusive) FiberTrackCont["mft12"] = fiberana->DeleteInclusive(FiberTrackCont["mft12"]);
      if(par->flag_debug) std::cout << "- after  Inclusive : " << FiberTrackCont["mft12"].size() << std::endl;
      nt_mft12 = FiberTrackCont["mft12"].size();

    }

    if(par->flag_mft12_pair){
      std::vector<std::pair<FiberHitAna*, FiberHitAna*> > pair_x;
      std::vector<std::pair<FiberHitAna*, FiberHitAna*> > pair_u;
      std::vector<std::pair<FiberHitAna*, FiberHitAna*> > pair_v;
      for(auto v1: FiberHitClCont[3][0]){
        for(auto v2: FiberHitClCont[4][0]){
          double pos1 = v1->GetPos();
          double pos2 = v2->GetPos();
          if( fabs(pos1 - pos2) < 20 ) pair_x.emplace_back(std::make_pair(v1, v2));
        }
      }
      for(auto v1: FiberHitClCont[3][1]){
        for(auto v2: FiberHitClCont[4][2]){
          double pos1 = v1->GetPos();
          double pos2 = v2->GetPos();
          if( fabs(pos1 - pos2) < 20 ) pair_u.emplace_back(std::make_pair(v1, v2));
        }
      }
      for(auto v1: FiberHitClCont[3][2]){
        for(auto v2: FiberHitClCont[4][1]){
          double pos1 = v1->GetPos();
          double pos2 = v2->GetPos();
          if( fabs(pos1 - pos2) < 20 ) pair_v.emplace_back(std::make_pair(v1, v2));
        }
      }

      int num_xp = pair_x.size();
      int num_up = pair_u.size();
      int num_vp = pair_v.size();
      int num_combi_pair = (num_xp+1) * (num_up+1) * (num_vp+1);
      if(par->flag_debug) std::cout << "- num_combi_pair_mft12 : " << num_combi_pair << std::endl;
      LocalHisto.hfiber_1_9->Fill((double)num_combi_pair*1e-3);

      if(num_combi_pair < par->cut_mft12_combi){

        std::vector<FiberTrackAna*> buf_track;
        for(int a=-1; a<num_xp; ++a){
          for(int b=-1; b<num_up; ++b){
            for(int c=-1; c<num_vp; ++c){
              std::vector<FiberHitAna*> buf_hit;
              int count = 0;
              if( a>-1 && !pair_x[a].first->IsUsed() && !pair_x[a].second->IsUsed() ){
                buf_hit.emplace_back(pair_x[a].first); buf_hit.emplace_back(pair_x[a].second); count+=2; }
              if( b>-1 && !pair_u[b].first->IsUsed() && !pair_u[b].second->IsUsed() ){
                buf_hit.emplace_back(pair_u[b].first); buf_hit.emplace_back(pair_u[b].second); count+=2; }
              if( c>-1 && !pair_v[c].first->IsUsed() && !pair_v[c].second->IsUsed() ){
                buf_hit.emplace_back(pair_v[c].first); buf_hit.emplace_back(pair_v[c].second); count+=2; }
              if(count<4) continue;
              FiberTrackAna *track = new FiberTrackAna(buf_hit, par.get());
              track->SetFlagPair();
              if(par->flag_mft12_posang){
                track->CorrectMFTCombi(par.get());
                //track->CorrectMFTCombi(par);
                double buf_x = track->GetXmft();
                double buf_y = track->GetYmft();
                double buf_a = track->GetA();
                double buf_b = track->GetB();
                if( fabs(buf_x * 0.003 - buf_a)>0.3 || fabs(buf_y * 0.003 - buf_b)>0.3 ){ delete track; continue;}
                //if( pow(buf_x, 2.) + pow(buf_y, 2.) > pow(100., 2.)){ delete track; continue;}
              }
              switch(track->GetNlayer()){
                case 4: buf_track.emplace_back(track); break;
                case 5: if(track->GetChi2() < par->cut_chi2_mft12) buf_track.emplace_back(track); else delete track; break;
                case 6: if(track->GetChi2() < par->cut_chi2_mft12) buf_track.emplace_back(track); else delete track; break;
                default: break;
              }

            }
          }
        }

        for(int i=0; i<(int)buf_track.size(); ++i){
          for(int j=0; j<(int)PSBHitCont.size(); ++j){
            double a_fiber = buf_track[i]->GetA();
            double b_fiber = buf_track[i]->GetB();
            double x_fiber = buf_track[i]->GetX();
            double y_fiber = buf_track[i]->GetY();
            //double phi_fiber = atan2(b_fiber, a_fiber);
            double phi_psb   = PSBHitCont[j]->GetPhi() + par->psb_rot_z*Deg2Rad;
            double r_psb     = PSBHitCont[j]->GetR();
            double z_psb     = PSBHitCont[j]->GetZ();

            double par_a = pow(a_fiber, 2) + pow(b_fiber, 2);
            double par_b = a_fiber * (x_fiber - par->psb_pos_x)+ b_fiber * (y_fiber - par->psb_pos_y);
            double par_c = pow(x_fiber - par->psb_pos_x, 2) + pow(y_fiber - par->psb_pos_y, 2) - pow(r_psb, 2);
            double z_fiber  = (-par_b + sqrt( pow(par_b,2) - par_a * par_c)) / par_a;

            double fiber_x_buf = x_fiber + a_fiber * z_fiber - par->psb_pos_x;
            double fiber_y_buf = y_fiber + b_fiber * z_fiber - par->psb_pos_y;
            double phi_fiber = atan2(fiber_y_buf, fiber_x_buf);
            if( fabs(fiberana->CalcPhiDif(phi_psb, phi_fiber)) < par->cut_psb_phi ){
              if(fabs( (z_fiber - par->psb_pos_z) - z_psb)<par->cut_psb_z){
                if(buf_track[i]->IsFlagPSB()){
                  double phi_psb_buf = buf_track[i]->GetPSBHit()->GetPhi();
                  if( fabs( fiberana->CalcPhiDif(phi_psb, phi_fiber) ) < fabs( fiberana->CalcPhiDif(phi_psb_buf, phi_fiber) ) ){
                    buf_track[i]->SetSegPSB(PSBHitCont[j]->GetSeg());
                    buf_track[i]->SetPSBHit(PSBHitCont[j]);
                  }
                }
                else{
                  buf_track[i]->SetFlagPSB();
                  buf_track[i]->SetSegPSB(PSBHitCont[j]->GetSeg());
                  buf_track[i]->SetPSBHit(PSBHitCont[j]);
                }
              }
            }

          }
        }

        int num_buf = buf_track.size();
        for(int i=num_buf-1; i>=0; --i){
          if(!buf_track[i]->IsFlagPSB()){
            delete buf_track[i];
            buf_track.erase(buf_track.begin() + i);
          }
        }

        ////////////////////////////
        //  delete dup ???? ////////
        ////////////////////////////

        for(auto v: buf_track){
          for(auto v2: v->GetContHit()){ v2->SetUsed(); }
          v->CorrectMFTCombi(par.get());
          //v->CorrectMFTCombi(par);
          //v->CorrectMFTCombi(par);
          v->SetPosL();
          v->DelFlagPSB();
          //if(par->flag_mftcor_xy) v->CorrectMFTXY(par);
        }

        for(auto v : buf_track){
          FiberTrackCont["mft12"].emplace_back(v);
        }

        FiberTrackCont["mft12"] = fiberana->DeleteSame(FiberTrackCont["mft12"]);
        if(par->flag_debug) std::cout << "- before Inclusive : " << FiberTrackCont["mft12"].size() << std::endl;
        if(par->flag_mft12_inclusive) FiberTrackCont["mft12"] = fiberana->DeleteInclusive(FiberTrackCont["mft12"]);
        if(par->flag_debug) std::cout << "- after  Inclusive : " << FiberTrackCont["mft12"].size() << std::endl;
        nt_mft12 = FiberTrackCont["mft12"].size();

      }

    }

    for(auto v: FiberTrackCont["mft12"]){
      v->SortContHit();
    }

    for(auto track: FiberTrackCont["mft12"]){
      double x1 = -9999.;
      double x2 = -9999.;
      double u1 = -9999.;
      double u2 = -9999.;
      double v1 = -9999.;
      double v2 = -9999.;
      for(auto hit: track->GetContHit()){
        double pos = hit->GetPos();
        if(hit->GetDid()==30) x1 = pos;
        if(hit->GetDid()==31) u1 = pos;
        if(hit->GetDid()==32) v1 = pos;
        if(hit->GetDid()==40) x2 = pos;
        if(hit->GetDid()==41) v2 = pos;
        if(hit->GetDid()==42) u2 = pos;
      }
      if(x1!=-9999. && x2!=-9999.) LocalHisto.hfiber_12_1_1->Fill(x1, x2);
      if(u1!=-9999. && u2!=-9999.) LocalHisto.hfiber_12_2_1->Fill(u1, u2);
      if(v1!=-9999. && v2!=-9999.) LocalHisto.hfiber_12_3_1->Fill(v1, v2);
      if(track->IsFlagPair()){
        if(x1!=-9999. && x2!=-9999.) LocalHisto.hfiber_12_1_2->Fill(x1, x2);
        if(u1!=-9999. && u2!=-9999.) LocalHisto.hfiber_12_2_2->Fill(u1, u2);
        if(v1!=-9999. && v2!=-9999.) LocalHisto.hfiber_12_3_2->Fill(v1, v2);
      }
    }

    LocalHisto.hfiber_3_0->Fill(nt_mft12);
    LocalHisto.hfiber_3_0_2->Fill(nt_mft12_xuv);

    if(par->flag_debug) std::cout << "- mft12 end" << std::endl;


  RecoEvent.FiberHitClCont = FiberHitClCont;

  int TypeDet[7][3] = {{G4Sol::FiberD1_x, G4Sol::FiberD1_u, G4Sol::FiberD1_v},
                       {G4Sol::FiberD2_x, G4Sol::FiberD2_u, G4Sol::FiberD2_v},
                       {G4Sol::FiberD3_x, G4Sol::FiberD3_u, G4Sol::FiberD3_v},
                       {G4Sol::MiniFiberD1_x, G4Sol::MiniFiberD1_u, G4Sol::MiniFiberD1_v},
                       {G4Sol::MiniFiberD2_x, G4Sol::MiniFiberD2_v, G4Sol::MiniFiberD2_u},
                       {G4Sol::FiberD4_v, G4Sol::FiberD4_u, G4Sol::FiberD4_x},
                       {G4Sol::FiberD5_x, G4Sol::FiberD5_u, G4Sol::FiberD5_v}};

  std::string nameVolDet[7][3] = {{"FiberD1_log_x", "FiberD1_log_u", "FiberD1_log_v"},
                                  {"FiberD2_log_x", "FiberD2_log_u", "FiberD2_log_v"},
                                  {"FiberD3_log_x", "FiberD3_log_u", "FiberD3_log_v"},
                                  {"MiniFiberD1_log_x", "MiniFiberD1_log_u", "MiniFiberD1_log_v"},
                                  {"MiniFiberD2_log_x", "MiniFiberD2_log_v", "MiniFiberD2_log_u"},
                                  {"FiberD4_log_v", "FiberD4_log_u", "FiberD4_log_x"},
                                  {"FiberD5_log_x", "FiberD5_log_u", "FiberD5_log_v"}};

  std::string nameFiber[7][3] = {{"FiberD1_Core_log_x", "FiberD1_Core_log_u", "FiberD1_Core_log_v"},
                                  {"FiberD2_Core_log_x", "FiberD2_Core_log_u", "FiberD2_Core_log_v"},
                                  {"FiberD3_Core_log_x", "FiberD3_Core_log_u", "FiberD3_Core_log_v"},
                                  {"MiniFiberD1_Core_log_x", "MiniFiberD1_Core_log_u", "MiniFiberD1_Core_log_v"},
                                  {"MiniFiberD2_Core_log_x", "MiniFiberD2_Core_log_v", "MiniFiberD2_Core_log_u"},
                                  {"FiberD4_Core_log_v", "FiberD4_Core_log_u", "FiberD4_Core_log_x"},
                                  {"FiberD5_Core_log_x", "FiberD5_Core_log_u", "FiberD5_Core_log_v"}};

  std::string nameMotherDet[7] = {"FiberD1_log_0",     "FiberD2_log_0", "FiberD3_log_0",
                                  "MiniFiberD1_log_0", "MiniFiberD2_log_0",
                                  "FiberD4_log_0",     "FiberD5_log_0"};

  for(int i = 0; i < 7; ++i)
    {
      for(int j = 0; j < 3; ++j)
        {
          // std::cout << "\nsize before : " << FiberHitCont[i][j].size() << std::endl;
          // std::cout << "size after  : " << FiberHitClCont[i][j].size() << std::endl;
          for(int k = 0; k < (int)FiberHitClCont[i][j].size(); ++k)
            {
              // FiberHitClCont[i][j][k]->Print();
              std::string volumeName = nameVolDet[i][j];
              std::string motherName = nameMotherDet[i];
              std::string fiberName = nameFiber[i][j];

              int hitID = FiberHitClCont[i][j][k]->GetClFib(); // CHECK

              TGeoMatrix* g1 =
                  gGeoManager->GetVolume(volumeName.c_str())->GetNode( (fiberName + Form("_%d", hitID*2)).c_str() )->GetMatrix(); //->GetNode(hitID * 2 + 1)->GetMatrix(); // fiber core
              TGeoMatrix* g1_pair =
                  gGeoManager->GetVolume(volumeName.c_str())->GetNode( (fiberName + Form("_%d", hitID*2 +1)).c_str() )->GetMatrix(); // minifiber pair core
              TGeoHMatrix H1(*g1), H1_pair(*g1_pair);
              Double_t* center_hit = H1.GetTranslation();
              Double_t* center_pair = H1_pair.GetTranslation();
              Double_t center_both[3];
              center_both[0] = (center_hit[0] + center_pair[0]) / 2.;
              center_both[1] = (center_hit[1] + center_pair[1]) / 2.;
              center_both[2] = (center_hit[2] + center_pair[2]) / 2.;

              H1.SetTranslation(center_both);

              TGeoMatrix* g2 = gGeoManager->GetVolume("MFLD")
                                   ->GetNode(motherName.c_str())
                                   ->GetVolume()
                                   ->GetNode((volumeName + "_0").c_str())
                                   ->GetMatrix(); // fiber layer
              TGeoMatrix* g3 =
                  gGeoManager->GetVolume("MFLD")->GetNode(motherName.c_str())->GetMatrix(); // fiber station
              TGeoMatrix* g4 = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();     // MFLD
              TGeoHMatrix H2(*g2), H3(*g3), H4(*g4);
              TGeoHMatrix H = H2 * H1;
              H             = H3 * H;
              H             = H4 * H;
              TGeoHMatrix w1("w1");
              TGeoHMatrix w2("w2");
              w1.SetDz(-10);
              w2.SetDz(10);
              TGeoHMatrix Hw1 = H * w1;
              TGeoHMatrix Hw2 = H * w2;
#ifdef DEBUG_BUILD2
              H.Print();
              Hw1.Print();
              Hw2.Print();
#endif
              double* edge1 = Hw1.GetTranslation();
              double* edge2 = Hw2.GetTranslation();
              // std::cout << "edge1[0] : " << edge1[0] << std::endl;
              // std::cout << "edge1[1] : " << edge1[1] << std::endl;
              // std::cout << "edge1[2] : " << edge1[2] << std::endl;
              // std::cout << "edge2[0] : " << edge2[0] << std::endl;
              // std::cout << "edge2[1] : " << edge2[1] << std::endl;
              // std::cout << "edge2[2] : " << edge2[2] << std::endl;
              double* shift = H.GetTranslation();
              // std::cout << "shift[0] : " << shift[0] << std::endl;
              // std::cout << "shift[1] : " << shift[1] << std::endl;
              // std::cout << "shift[2] : " << shift[2] << std::endl;
              //TVector3 o(0., 0., shift[2]), zdir(0., 0., 1.);
              TVector3 o(0., 0., FiberHitClCont[i][j][k]->GetZ()*0.1), zdir(0., 0., 1.);
              TVector3 fiber_dir(edge2[0] - edge1[0], edge2[1] - edge1[1], edge2[2] - edge1[2]);
              fiber_dir  = fiber_dir.Unit();
              TVector3 u = fiber_dir.Cross(zdir);
              TVector3 v = fiber_dir;
              genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));

              TVectorD hitCoords(1);
              hitCoords(0) = FiberHitClCont[i][j][k]->GetPos()*0.1; //u.Dot(TVector3(shift[0], shift[1], 0));
              TMatrixDSym hitCov(1);
              //hitCov(0, 0) = TMath::Sq(FiberHitClCont[i][j][k]->GetRes()*0.1); // mm to cm // to be adjusted resolution_fiber * resolution_fiber; // Check with or without Square root?
              hitCov(0, 0) = std::pow(FiberHitClCont[i][j][k]->GetRes()*0.1, 2.); // mm to cm
              auto measurement =
                  std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, TypeDet[i][j], hitID, nullptr);
              dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

              RecoEvent.ListHits[TypeDet[i][j]].emplace_back(measurement.release());

              MeasurementInfo measinfo(FiberHitClCont[i][j][k]->GetTOT(), FiberHitClCont[i][j][k]->GetTime(), -9999.);
              RecoEvent.ListHitsInfo[TypeDet[i][j]].emplace_back(measinfo);
            }
        }
    }

  for(int i = 0; i < 7; ++i)
    {
      for(int j = 0; j < 3; ++j)
        {
          LocalHisto.h13[i][j]->Fill(FiberHitCont[i][j].size());
          LocalHisto.h14[i][j]->Fill(FiberHitClCont[i][j].size());

          for(int k = 0; k < (int)FiberHitCont[i][j].size(); ++k)
            {
                            //printf("%d %d %.8f\n", i, j, FiberHitCont[i][j][k]->GetPos());

              LocalHisto.h10[i][j]->Fill(FiberHitCont[i][j][k]->GetFib(), FiberHitCont[i][j][k]->GetTL());
              LocalHisto.h11[i][j]->Fill(FiberHitCont[i][j][k]->GetPos());
            }

          for(int k = 0; k < (int)FiberHitClCont[i][j].size(); ++k)
            {
              LocalHisto.h12[i][j]->Fill(FiberHitClCont[i][j][k]->GetPos());
              LocalHisto.h15[i][j]->Fill(FiberHitClCont[i][j][k]->GetClsize());
              LocalHisto.h75[i][j]->Fill(FiberHitClCont[i][j][k]->GetClFib(), FiberHitClCont[i][j][k]->GetPos());
            }
        }
    }



    for(int i=0; i<7; ++i)
      {
        LocalHisto.h17[i]->Fill( FiberXUVCont[i].size() );
        for(int j=0; j<(int)FiberXUVCont[i].size(); ++j)
          {
            LocalHisto.h16[i]->Fill(FiberXUVCont[i][j]->GetPosX(), FiberXUVCont[i][j]->GetPosY());
            LocalHisto.h17_2[i]->Fill(FiberXUVCont[i][j]->GetD());
          }
      }



    // dft12
    int nt_dft12 = 0;

    if(FiberXUVCont[5].size()>0 && FiberXUVCont[6].size()>0){
      std::vector<FiberTrackAna*> buf_track;
      for(int i=0; i<(int)FiberXUVCont[5].size(); ++i){
        for(int j=0; j<(int)FiberXUVCont[6].size(); ++j){
          std::vector<FiberHitXUV*>   buf_xuv;
          buf_xuv.emplace_back(FiberXUVCont[5][i]);
          buf_xuv.emplace_back(FiberXUVCont[6][j]);
          FiberTrackAna *track = new FiberTrackAna(buf_xuv, par.get());
          if(par->flag_dft12_cut){
            double tot_mean = track->GetTOT();
            double x_buf = track->GetXdet();
            double y_buf = track->GetYdet();
            double a_buf = track->GetA() * 1000;
            double b_buf = track->GetB() * 1000;
            bool flag_cut = false;
            if( pow(x_buf/60., 2.) + pow(y_buf/40., 2) > 1. ) flag_cut = true;
            if( fabs( x_buf * 23./50. - a_buf) > 20. ) flag_cut = true;
            if( fabs( y_buf * 35./60. - b_buf) > 30. ) flag_cut = true;
            if( tot_mean < par->cut_dft12_tot_max ) flag_cut = true;
            if(flag_cut){ delete track; continue; }
          }
          if(!par->flag_dft12_combi)                      buf_track.emplace_back(track);
          else if(track->GetChi2() < par->cut_chi2_dft12) buf_track.emplace_back(track);
        }
      }
      if((int)buf_track.size()>0) FiberTrackCont["dft12"] = fiberana->DeleteDup(buf_track);
      nt_dft12     = FiberTrackCont["dft12"].size();
    }

    LocalHisto.hfiber_4_1->Fill(nt_dft12);
    for(int i=0; i<nt_dft12; ++i){
      FiberTrackAna *track = FiberTrackCont["dft12"][i];
      double tot_mean = track->GetTOT();
      if(nt_dft12==1){
        LocalHisto.hfiber_4_2_1->Fill(track->GetXdet(), track->GetYdet());
        LocalHisto.hfiber_4_3_1->Fill(track->GetXdet(), track->GetA()*1000);
        LocalHisto.hfiber_4_4_1->Fill(track->GetYdet(), track->GetB()*1000);
        LocalHisto.hfiber_4_5_1->Fill(tot_mean);
      }
      if(nt_dft12>1){
        LocalHisto.hfiber_4_2_2->Fill(track->GetXdet(), track->GetYdet());
        LocalHisto.hfiber_4_3_2->Fill(track->GetXdet(), track->GetA()*1000);
        LocalHisto.hfiber_4_4_2->Fill(track->GetYdet(), track->GetB()*1000);
        LocalHisto.hfiber_4_5_2->Fill(tot_mean);
      }
    }


    for(auto v: FiberTrackCont["dft12"]){
      for(int i=0; i<6; ++i){
        v->GetContHit().at(i)->SetUsed();
      }
    }

    int num_combi_dft12 = 1;
    for(int i=5; i<7; ++i){
      for(int j=0; j<3; ++j){
        num_combi_dft12 *= ( (int)FiberHitClCont[i][j].size() + 1 );
      }
    }
    LocalHisto.hfiber_5_5->Fill((double)num_combi_dft12*1e-3);
    LocalHisto.hfiber_5_6->Fill((double)num_combi_dft12*1e-3);
    LocalHisto.hfiber_5_7->Fill((double)num_combi_dft12*1e-6);
    if(par->flag_debug) std::cout << "- num_combi_dft12 : " << num_combi_dft12 << std::endl;

    if(par->flag_dft12_combi && num_combi_dft12<par->cut_dft12_combi){
      LocalHisto.hfiber_5_1->Fill(num_combi_dft12);
      LocalHisto.hfiber_5_2->Fill((double)num_combi_dft12*1e-3);
      LocalHisto.hfiber_5_3->Fill((double)num_combi_dft12*1e-3);
      LocalHisto.hfiber_5_4->Fill((double)num_combi_dft12*1e-6);
      std::vector<FiberTrackAna*> buf_track;
      for(int a=-1; a<(int)FiberHitClCont[5][0].size(); ++a){
        for(int b=-1; b<(int)FiberHitClCont[5][1].size(); ++b){
          for(int c=-1; c<(int)FiberHitClCont[5][2].size(); ++c){
            for(int d=-1; d<(int)FiberHitClCont[6][0].size(); ++d){
              for(int e=-1; e<(int)FiberHitClCont[6][1].size(); ++e){
                for(int f=-1; f<(int)FiberHitClCont[6][2].size(); ++f){
                  std::vector<FiberHitAna*> buf_hit;
                  int count = 0;
                  if(a>-1 && !FiberHitClCont[5][0][a]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[5][0][a]); count++;}
                  if(b>-1 && !FiberHitClCont[5][1][b]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[5][1][b]); count++;}
                  if(c>-1 && !FiberHitClCont[5][2][c]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[5][2][c]); count++;}
                  if(d>-1 && !FiberHitClCont[6][0][d]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[6][0][d]); count++;}
                  if(e>-1 && !FiberHitClCont[6][1][e]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[6][1][e]); count++;}
                  if(f>-1 && !FiberHitClCont[6][2][f]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[6][2][f]); count++;}
                  if(count<4) continue;
                  //if(count!=4) continue;
                  FiberTrackAna *track = new FiberTrackAna(buf_hit, par.get());

                  if(par->flag_dft12_cut){
                    double tot_mean = track->GetTOT();
                    double x_buf = track->GetXdet();
                    double y_buf = track->GetYdet();
                    double a_buf = track->GetA() * 1000;
                    double b_buf = track->GetB() * 1000;
                    bool flag_cut = false;
                    if( pow(x_buf/60., 2.) + pow(y_buf/40., 2) > 1. ) flag_cut = true;
                    if( fabs( x_buf * 23./50. - a_buf) > 20. ) flag_cut = true;
                    if( fabs( y_buf * 35./60. - b_buf) > 30. ) flag_cut = true;
                    if( tot_mean < 75. ) flag_cut = true;
                    if(flag_cut){ delete track; continue; }
                  }

                  switch(track->GetNlayer()){
                    case 4: buf_track.emplace_back(track); break;
                    case 5: if(track->GetChi2() < par->cut_chi2_dft12) buf_track.emplace_back(track); else delete track; break;
                    case 6: if(track->GetChi2() < par->cut_chi2_dft12) buf_track.emplace_back(track); else delete track; break;
                    default: break;
                  }
                }
              }
            }
          }
        }
      }
      if(par->flag_debug) std::cout << "- dft12 combi" << std::endl;

      std::vector<FiberTrackAna*> buf_track_tmp = fiberana->DeleteDupCombi(buf_track);
      if(par->flag_debug) std::cout << "- dft12 dup" << std::endl;
      buf_track = buf_track_tmp;

      for(auto v: buf_track){
        for(auto v2: v->GetContHit()){ v2->SetUsed(); }
      }

      for(auto v : buf_track){
        FiberTrackCont["dft12"].emplace_back(v);
      }
      nt_dft12 = FiberTrackCont["dft12"].size();

    }

    LocalHisto.hfiber_4_1_3->Fill(nt_dft12);
    int buf_i_dft12 = -1;
    double buf_diff_dft12 = -9999.;
    for(int i=0; i<nt_dft12; ++i){
      FiberTrackAna *track = FiberTrackCont["dft12"][i];
      double tot_mean  = track->GetTOT();
      double time_mean = track->GetTime();
      LocalHisto.hfiber_4_2_3->Fill(track->GetXdet(), track->GetYdet());
      LocalHisto.hfiber_4_3_3->Fill(track->GetXdet(), track->GetA()*1000);
      LocalHisto.hfiber_4_4_3->Fill(track->GetYdet(), track->GetB()*1000);
      LocalHisto.hfiber_4_5_3->Fill(tot_mean);
      double dist_dft12 =
        pow( (tot_mean  - par->cut_dft12_tot_mean)  / par->cut_dft12_tot_sig , 2.) +
        pow( (time_mean - par->cut_dft12_time_mean) / par->cut_dft12_time_sig, 2.);
      if(buf_diff_dft12<0 || dist_dft12 < buf_diff_dft12){
        buf_i_dft12 = i;
        buf_diff_dft12 = dist_dft12;
      }
    }

    if(buf_i_dft12>-1) FiberTrackCont["dft12"][buf_i_dft12]->SetBest();


    LocalHisto.h18_3_1->Fill(nt_dft12);
    if(nt_dft12>0){
      LocalHisto.h18_3_2->Fill(FiberTrackCont["dft12"][0]->GetChi2());
      LocalHisto.h18_3_3->Fill(FiberTrackCont["dft12"][0]->GetXtgt(), FiberTrackCont["dft12"][0]->GetA()*1000);
      LocalHisto.h18_3_4->Fill(FiberTrackCont["dft12"][0]->GetYtgt(), FiberTrackCont["dft12"][0]->GetB()*1000);
      LocalHisto.h18_3_5->Fill(FiberTrackCont["dft12"][0]->GetXtgt());
      LocalHisto.h18_3_6->Fill(FiberTrackCont["dft12"][0]->GetYtgt());
      LocalHisto.h18_3_7->Fill(FiberTrackCont["dft12"][0]->GetA()*1000);
      LocalHisto.h18_3_8->Fill(FiberTrackCont["dft12"][0]->GetB()*1000);
    }

    RecoEvent.FiberTrackCont = FiberTrackCont;

    if(par->flag_debug) std::cout << "- dft12 end" << std::endl;


    //  S4 Scintillators  //
    S4SciHitAna *s4hit = new S4SciHitAna(event.s4tq, par.get(), att.StudyCase);

    RecoEvent.FragmentPID = s4hit->GetPID();

    LocalHisto.hs4sci_1_1->Fill(s4hit->GetdE_sc31());
    LocalHisto.hs4sci_1_2->Fill(s4hit->GetdE_sc41());
    LocalHisto.hs4sci_1_3->Fill(s4hit->GetdE_sc42_high());
    LocalHisto.hs4sci_1_4->Fill(s4hit->GetdE_sc42_low());
    LocalHisto.hs4sci_2_1->Fill(s4hit->GetTOF_sc3141() , s4hit->GetdE_sc31());
    LocalHisto.hs4sci_2_2->Fill(s4hit->GetTOF_sc3141() , s4hit->GetdE_sc41());
    LocalHisto.hs4sci_2_3->Fill(s4hit->GetTOF_sc3141() , s4hit->GetdE_sc42_high());
    LocalHisto.hs4sci_2_4->Fill(s4hit->GetTOF_sc3141() , s4hit->GetdE_sc42_low());
    //s4hit->Print();

    /// --- MWDC HIT ---- ///////////////////////////////////////////////
    int nt_mwdc = 0;
    //std::vector<MWDCHitAna *> MWDCHitCont;
    //MWDCHitCont.clear();
    //for(auto mwdc_hit : s4mwdc->mwdchit){
    //  MWDCHitAna *hit_ana = new MWDCHitAna(mwdc_hit, par, s4mwdc->tref[0]);
    //  MWDCHitCont.emplace_back(hit_ana);
    //}

    float mwdc_x = -999.;
    float mwdc_a = -999.;
    float mwdc_y = -999.;
    float mwdc_b = -999.;
    float mwdc_chi2 = -999.;
    MWDCTracking mwdc_track(par.get());// = new MWDCTracking(par);
    // for(auto hit_mwdc : MWDCHitCont){
    //   mwdc_track.StoreHit(hit_mwdc, par); /// only leading time cut
    // }

    for(auto mwdc_hit : event.s4mwdc->mwdchit){
      float lt_tmp = mwdc_hit.t_leading - event.s4mwdc->tref[0];
      float tot_tmp = mwdc_hit.t_trailing - mwdc_hit.t_leading;
      // if(tot_tmp<200) continue;
      // if(lt_tmp>-600 || lt_tmp<-850) continue;
      if(lt_tmp < par->mwdc_lt_valid_min || lt_tmp > par->mwdc_lt_valid_max) continue;
      float mwdc_length = par->mwdc_dtdxconversion(mwdc_hit.i_plane, lt_tmp);
      float wirepos_tmp = (par->mwdc_plane_sign[mwdc_hit.i_plane])*(5.0)*((float)(mwdc_hit.i_wire)-(par->mwdc_center_id[mwdc_hit.i_plane]));
      mwdc_track.StoreHit(mwdc_hit.i_plane,
          -wirepos_tmp,
          fabs(mwdc_length),
          tot_tmp,
          lt_tmp,
          -par->mwdc_plane_angle[mwdc_hit.i_plane],
          par->mwdc_zpos[mwdc_hit.i_plane]);
    }

    mwdc_track.SetAssumedResolution(par->mwdc_assumed_plane_resolution);
    mwdc_track.SetMinPlaneEnebaled(par->mwdc_min_plane_with_hit_for_tracking);   // 12,13,14,15,16 or something like that
    mwdc_track.SetMaxHitCombination(par->mwdc_max_hit_combination_for_tracking); // around 10. small number for online analysis!

    // bool tracking_status = mwdc_track.Tracking();
    int tracking_status = mwdc_track.Tracking_PairLR();  /// which is better? seems this is faster
    if(tracking_status==1){
      mwdc_chi2 = mwdc_track.GetChi2();
      if(mwdc_chi2 >= 0.0 && mwdc_chi2 < par->mwdc_tracking_chi2cut_max){
        mwdc_x = mwdc_track.GetX() + mwdc_track.GetA() * (par->dist_focS4 - par->dist_MWDC_zref)
          + par->mwdc_shift_x_alignment;
        mwdc_a = mwdc_track.GetA() * 1000. + par->mwdc_shift_a_alignment;
        mwdc_y = mwdc_track.GetY() + mwdc_track.GetB() * (par->dist_focS4 - par->dist_MWDC_zref)
          + par->mwdc_shift_y_alignment;
        mwdc_b = mwdc_track.GetB() * 1000. + par->mwdc_shift_b_alignment;
        nt_mwdc++;
        //for(int i=0; i<16; ++i){ //CHECK LATER
        //  res_mwdc[i] = mwdc_track.GetResidual(i);
        //  hmwdc_2[i]->Fill(res_mwdc[i]);
        //}
        LocalHisto.hmwdc_1_1->Fill(mwdc_x);
        LocalHisto.hmwdc_1_2->Fill(mwdc_y);
        LocalHisto.hmwdc_1_3->Fill(mwdc_a);
        LocalHisto.hmwdc_1_4->Fill(mwdc_b);
        LocalHisto.hmwdc_1_5->Fill(mwdc_chi2);

        mwdc_track.SetX(mwdc_x);
        mwdc_track.SetY(mwdc_y);
        mwdc_track.SetA(mwdc_a);
        mwdc_track.SetB(mwdc_b);


        //t_mwdc_x.emplace_back(mwdc_x); //Check include later to output root file
        //t_mwdc_y.emplace_back(mwdc_y);
        //t_mwdc_a.emplace_back(mwdc_a);
        //t_mwdc_b.emplace_back(mwdc_b);
        //t_mwdc_chi2.emplace_back(mwdc_chi2);
      }
    }

    RecoEvent.MWDCTracks.emplace_back(mwdc_track);

    //t_nt_mwdc = nt_mwdc;
    LocalHisto.hmwdc_1_6->Fill(nt_mwdc);
    //  MWDC end  ///////////////////////////////////////////



  // clear /////////////////////////////////////

  for(int i = 0; i < 7; ++i)
    {
      for(int j = 0; j < 3; ++j)
        {
          int num = (int)FiberHitCont[i][j].size();
          for(int k = 0; k < num; ++k)
            {
              delete FiberHitCont[i][j].back();
              FiberHitCont[i][j].pop_back();
            }
          int num2 = (int)FiberHitClCont[i][j].size();
          for(int k = 0; k < num2; ++k)
            {
              //delete FiberHitClCont[i][j].back(); //CHECK when removing FiberHitClCont from RecoEvent
              //FiberHitClCont[i][j].pop_back();
            }
        }
    }

  {
    int num = (int)PSBHitCont.size();
    for(int i = 0; i < num; ++i)
      {
        delete PSBHitCont.back();
        PSBHitCont.pop_back();
      }
  }

  for(int i = 0; i < 17; ++i)
    {
      int num = (int)MDCHitCont[i].size();
      for(int j = 0; j < num; ++j)
        {
          delete MDCHitCont[i].back();
          MDCHitCont[i].pop_back();
        }
    }

  // #ifdef DEBUG_BUILD2
  //                   std::cout << "fiber" << std::endl;
  //                   std::string tmpName = orderDetName.find(TypeDet)->second;
  //                   std::cout << "name : " << tmpName << std::endl;
  //                   std::cout << "LayerID : " << LayerID << std::endl;
  //                   std::cout << "HitPosX : " << hit.HitPosX << std::endl;
  //                   std::cout << "HitPosY : " << hit.HitPosY << std::endl;
  //                   std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
  //                   gGeoManager->GetVolume(volumeName.c_str())->GetNode(LayerID * 2 + 1)->Print();
  //                   gGeoManager->GetVolume(volumeName.c_str())->GetNode(LayerID * 2 + 1)->GetMatrix()->Print();
  //                   gGeoManager->GetVolume("MFLD")
  //                       ->GetNode(motherName.c_str())
  //                       ->GetVolume()
  //                       ->GetNode((volumeName + "_0").c_str())
  //                       ->Print();
  //                   gGeoManager->GetVolume("MFLD")
  //                       ->GetNode(motherName.c_str())
  //                       ->GetVolume()
  //                       ->GetNode((volumeName + "_0").c_str())
  //                       ->GetMatrix()
  //                       ->Print();
  //                   gGeoManager->GetVolume("MFLD")->GetNode(motherName.c_str())->Print();
  //                   gGeoManager->GetVolume("MFLD")->GetNode(motherName.c_str())->GetMatrix()->Print();
  // #endif
  //                   TGeoMatrix* g1 =
  //                       gGeoManager->GetVolume(volumeName.c_str())->GetNode(LayerID * 2 + 1)->GetMatrix(); // fiber
  //                       core
  //                   TGeoMatrix* g2 = gGeoManager->GetVolume("MFLD")
  //                                        ->GetNode(motherName.c_str())
  //                                        ->GetVolume()
  //                                        ->GetNode((volumeName + "_0").c_str())
  //                                        ->GetMatrix(); // fiber layer
  //                   TGeoMatrix* g3 =
  //                       gGeoManager->GetVolume("MFLD")->GetNode(motherName.c_str())->GetMatrix(); // fiber station
  //                   TGeoMatrix* g4 = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();     // MFLD
  //                   TGeoHMatrix H1(*g1), H2(*g2), H3(*g3), H4(*g4);
  //                   TGeoHMatrix H = H2 * H1;
  //                   H             = H3 * H;
  //                   H             = H4 * H;
  //                   TGeoHMatrix w1("w1");
  //                   TGeoHMatrix w2("w2");
  //                   w1.SetDz(-10);
  //                   w2.SetDz(10);
  //                   TGeoHMatrix Hw1 = H * w1;
  //                   TGeoHMatrix Hw2 = H * w2;
  // #ifdef DEBUG_BUILD2
  //                   H.Print();
  //                   Hw1.Print();
  //                   Hw2.Print();
  // #endif
  //                   double* edge1 = Hw1.GetTranslation();
  //                   double* edge2 = Hw2.GetTranslation();
  //                   double* shift = H.GetTranslation();
  //                   TVector3 o(0., 0., shift[2]), zdir(0., 0., 1.);
  //                   TVector3 fiber_dir(edge2[0] - edge1[0], edge2[1] - edge1[1], edge2[2] - edge1[2]);
  //                   fiber_dir  = fiber_dir.Unit();
  //                   TVector3 u = fiber_dir.Cross(zdir);
  //                   TVector3 v = fiber_dir;
  //                   genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));

  //                   TVectorD hitCoords(1);
  //                   hitCoords(0) = u.Dot(TVector3(shift[0], shift[1], 0));
  //                   TMatrixDSym hitCov(1);
  //                   hitCov(0, 0) = resolution_fiber * resolution_fiber;
  //                   measurement =
  //                       std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID,
  //                       nullptr);
  //                   dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

  //                   hitCoordsTree(0) = hit.HitPosX;
  //                   hitCoordsTree(1) = hit.HitPosY;
  //                   hitCoordsTree(2) = hit.HitPosZ;
  //                 }

  //               else if(IsWire(TypeDet))
  //                 {
  // #ifdef DEBUG_BUILD2
  //                   std::cout << "wire" << std::endl;
  //                   std::string tmpName = orderDetName.find(TypeDet)->second;
  //                   std::cout << "name : " << tmpName << std::endl;
  //                   std::cout << "LayerID : " << LayerID << std::endl;
  //                   std::cout << "HitPosX : " << hit.HitPosX << std::endl;
  //                   std::cout << "HitPosY : " << hit.HitPosY << std::endl;
  //                   std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
  //                   gGeoManager->GetVolume("INNER")
  //                       ->GetNode(TypeDet - G4Sol::MG01 + 1)
  //                       ->GetVolume()
  //                       ->GetNode(LayerID - 1)
  //                       ->Print();
  //                   gGeoManager->GetVolume("INNER")
  //                       ->GetNode(TypeDet - G4Sol::MG01 + 1)
  //                       ->GetVolume()
  //                       ->GetNode(LayerID - 1)
  //                       ->GetMatrix()
  //                       ->Print();
  //                   gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 1)->Print();
  //                   gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 1)->GetMatrix()->Print();
  //                   gGeoManager->GetVolume("MFLD")->GetNode(0)->Print();
  //                   gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix()->Print();
  //                   gGeoManager->GetVolume("WASA")->GetNode(0)->Print();
  //                   gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix()->Print();
  // #endif
  //                   TGeoMatrix* g1 = gGeoManager->GetVolume("INNER")
  //                                        ->GetNode(TypeDet - G4Sol::MG01 + 1)
  //                                        ->GetVolume()
  //                                        ->GetNode(LayerID - 1)
  //                                        ->GetMatrix(); // ME, MG
  //                   TGeoShape* tempShape =
  //                       gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 1)->GetVolume()->GetShape();
  //                   TGeoMatrix* g2 =
  //                       gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 1)->GetMatrix(); // MD
  //                   TGeoMatrix* g3 = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix();             // INNER
  //                   TGeoMatrix* g4 = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();             // MFLD
  //                   TGeoHMatrix H1(*g1), H2(*g2), H3(*g3), H4(*g4);
  //                   TGeoHMatrix H = H2 * H1;
  //                   H             = H3 * H;
  //                   H             = H4 * H;
  //                   double* shift = H.GetTranslation();
  //                   TGeoHMatrix w1("w1");
  //                   TGeoHMatrix w2("w2");
  //                   Double_t minZ, maxZ;
  //                   tempShape->GetAxisRange(3, minZ, maxZ);
  //                   w1.SetDz(minZ);
  //                   w2.SetDz(maxZ);
  //                   TGeoHMatrix Hw1 = H * w1;
  //                   TGeoHMatrix Hw2 = H * w2;
  // #ifdef DEBUG_BUILD2
  //                   H.Print();
  //                   Hw1.Print();
  //                   Hw2.Print();
  // #endif
  //                   double* edge1 = Hw1.GetTranslation();
  //                   double* edge2 = Hw2.GetTranslation();

  //                   TVector3 x1(shift[0], shift[1], shift[2]);
  //                   TVector3 p1(edge2[0] - edge1[0], edge2[1] - edge1[1], edge2[2] - edge1[2]);
  //                   TVector3 x2(hit.HitPosX, hit.HitPosY, hit.HitPosZ);
  //                   TVector3 p2(hit.MomX, hit.MomY, hit.MomZ);
  //                   double dl = CloseDist(x1, x2, p1, p2);

  //                   double dlmax = 0;
  //                   switch(TypeDet - G4Sol::MG01 + 1)
  //                     {
  //                     case 1:
  //                     case 2:
  //                     case 3:
  //                     case 4:
  //                     case 5:
  //                       dlmax = 0.2;
  //                       break;
  //                     case 6:
  //                     case 7:
  //                     case 8:
  //                     case 9:
  //                     case 10:
  //                     case 11:
  //                       dlmax = 0.3;
  //                       break;
  //                     case 12:
  //                     case 13:
  //                     case 14:
  //                     case 15:
  //                     case 16:
  //                     case 17:
  //                       dlmax = 0.4;
  //                       break;
  //                     default:
  //                       att._logger->warn("Error in WireMeasurement !");
  //                       break;
  //                     }
  //                   double temp_dl = gRandom->Gaus(dl, resolution_dl);
  //                   bool doneRand  = false;
  //                   while(doneRand)
  //                     {
  //                       if(temp_dl < 0 || temp_dl > dlmax)
  //                         temp_dl = gRandom->Gaus(dl, resolution_dl);
  //                       else
  //                         doneRand = true;
  //                     }
  //                   // if(temp_dl<0)     dl = 0;
  //                   // if(temp_dl>dlmax) dl = dlmax;

  //                   TVectorD hitCoords(7);
  //                   hitCoords(0) = edge1[0];
  //                   hitCoords(1) = edge1[1];
  //                   hitCoords(2) = edge1[2];
  //                   hitCoords(3) = edge2[0];
  //                   hitCoords(4) = edge2[1];
  //                   hitCoords(5) = edge2[2];
  //                   hitCoords(6) = temp_dl;
  //                   if(edge1[2] > edge2[2])
  //                     {
  //                       hitCoords(0) = edge2[0];
  //                       hitCoords(1) = edge2[1];
  //                       hitCoords(2) = edge2[2];
  //                       hitCoords(3) = edge1[0];
  //                       hitCoords(4) = edge1[1];
  //                       hitCoords(5) = edge1[2];
  //                     }
  //                   TMatrixDSym hitCov(7);
  //                   hitCov(6, 6) = resolution_dl * resolution_dl;
  //                   measurement =
  //                       std::make_unique<genfit::WireMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID, nullptr);
  //                   dynamic_cast<genfit::WireMeasurement*>(measurement.get())->setLeftRightResolution(0);
  //                   dynamic_cast<genfit::WireMeasurement*>(measurement.get())->setMaxDistance(dlmax);

  //                   hitCoordsTree(0) = hit.HitPosX;
  //                   hitCoordsTree(1) = hit.HitPosY;
  //                   hitCoordsTree(2) = hit.HitPosZ;
  //                 }

  return 0;
}
