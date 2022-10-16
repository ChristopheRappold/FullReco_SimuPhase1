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
      if(tempName == "MDC17_1")
	index_lastMDC = i;
      if(tempName == "PSCE_1")
	index_firstPSCE = i;
    }
  offsetGeoNameID_PSCE = index_firstPSCE - index_lastMDC + offsetGeoNameID_MDC + 17;

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
          LocalHisto.h10[i][j] = AnaHisto->CloneAndRegister(AnaHisto->h10[i][j]);
          LocalHisto.h11[i][j] = AnaHisto->CloneAndRegister(AnaHisto->h11[i][j]);
          LocalHisto.h12[i][j] = AnaHisto->CloneAndRegister(AnaHisto->h12[i][j]);
          LocalHisto.h13[i][j] = AnaHisto->CloneAndRegister(AnaHisto->h13[i][j]);
          LocalHisto.h14[i][j] = AnaHisto->CloneAndRegister(AnaHisto->h14[i][j]);
          LocalHisto.h15[i][j] = AnaHisto->CloneAndRegister(AnaHisto->h15[i][j]);
        }
      LocalHisto.h16[i]   = AnaHisto->CloneAndRegister(AnaHisto->h16[i]);
      LocalHisto.h17[i]   = AnaHisto->CloneAndRegister(AnaHisto->h17[i]);
      LocalHisto.h17_2[i] = AnaHisto->CloneAndRegister(AnaHisto->h17_2[i]);
    }

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
  for(int i = 0; i < 17; ++i)
    {
      LocalHisto.hmdc_1[i]   = AnaHisto->CloneAndRegister(AnaHisto->hmdc_1[i]);
      LocalHisto.hmdc_2[i]   = AnaHisto->CloneAndRegister(AnaHisto->hmdc_2[i]);
      LocalHisto.hmdc_2_2[i] = AnaHisto->CloneAndRegister(AnaHisto->hmdc_2_2[i]);
      LocalHisto.hmdc_2_3[i] = AnaHisto->CloneAndRegister(AnaHisto->hmdc_2_3[i]);
      LocalHisto.hmdc_3[i]   = AnaHisto->CloneAndRegister(AnaHisto->hmdc_3[i]);
      LocalHisto.hmdc_3_2[i] = AnaHisto->CloneAndRegister(AnaHisto->hmdc_3_2[i]);
      LocalHisto.hmdc_3_3[i] = AnaHisto->CloneAndRegister(AnaHisto->hmdc_3_3[i]);
    }
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
  att._logger->info("Start Build Calibration of Wasa Detector !");
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
          int t_u            = event.s2tq1->tdc_psb[i][0][i_u];
          int t_d            = event.s2tq1->tdc_psb[i][1][i_d];
          int q_u            = event.s2tq1->qdc_psb[i][0];
          int q_d            = event.s2tq1->qdc_psb[i][1];
          PSBHitAna* hit_ana = new PSBHitAna(i, t_u, t_d, q_u, q_d, par.get());
          PSBHitCont.emplace_back(hit_ana);
        }
    }

  for(const auto& hit_PSB : PSBHitCont)
    {
      int hitID      = hit_PSB->GetSeg();
      TGeoMatrix* g1 = gGeoManager->GetVolume("INNER")->GetNode(offsetGeoNameID_PSCE + hitID)->GetMatrix(); // PSCE
      TGeoMatrix* g2 = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix();           // INNNER
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
      hitCov(0, 0)     = TMath::Sq(par->psb_res_phi * 0.1);//0.1; // to be adjusted resolution_psce * resolution_psce;
      hitCov(1, 1)     = TMath::Sq(par->psb_res_z   * 0.1);//0.1; // to be adjusted resolution_psce_z * resolution_psce_z;
      auto measurement = std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, G4Sol::PSCE, hitID, nullptr);
      dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

      RecoEvent.ListHits[G4Sol::PSCE].emplace_back(measurement.release());
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
            int t_buf           = event.s2tq2->tdc_psfe[i][j];
            int q_buf           = event.s2tq2->qdc_psfe[i];
            PSFEHitAna* hit_ana = new PSFEHitAna(i, t_buf, q_buf, par.get());
            PSFEHitCont.emplace_back(hit_ana);
          }
      }

  for(const auto& hit_PSFE : PSFEHitCont)
    {
      int hitID = hit_PSFE->GetSeg();

      TGeoMatrix* g1   = gGeoManager->GetVolume("PSF")->GetNode(hitID)->GetMatrix();    // PSCE
      TGeoMatrix* g1_1 = gGeoManager->GetVolume("MFLD")->GetNode("PSF_1")->GetMatrix(); // PSB box
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
      hitCov(0, 0) = TMath::Sq(par->psfe_res_phi * 0.1); // mm to cm

      auto measurement = std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, G4Sol::PSFE, hitID, nullptr);
      dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

      RecoEvent.ListHits[G4Sol::PSFE].emplace_back(measurement.release());
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

      TGeoMatrix* g1   = gGeoManager->GetVolume("PSB")->GetNode(hitID)->GetMatrix();    // PSCE
      TGeoMatrix* g1_1 = gGeoManager->GetVolume("MFLD")->GetNode("PSB_1")->GetMatrix(); // PSB box
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
  for(int i = 0; i < (int)T0HitCont.size(); ++i)
    {
      LocalHisto.ht0_0_1->Fill(T0HitCont[i]->GetTU());
      LocalHisto.ht0_0_2->Fill(T0HitCont[i]->GetTD());
      LocalHisto.ht0_0_4->Fill(T0HitCont[i]->GetSeg());
      LocalHisto.ht0_1[T0HitCont[i]->GetSeg()]->Fill(T0HitCont[i]->GetTime());
    }
  LocalHisto.ht0_0_3->Fill(T0HitCont.size());

  // MDC ana //////////////////////////////////
  std::vector<std::vector<MDCHitAna*> > MDCHitCont;
  for(int i = 0; i < 17; ++i)
    {
      std::vector<MDCHitAna*> buf_v;
      MDCHitCont.emplace_back(buf_v);
    }
  for(int i = 0; i < (int)event.s2mdc->mdchit.size(); ++i)
    {
      MDCHitAna* hit_ana = new MDCHitAna(event.s2mdc->mdchit[i], par.get(), event.s2mdc->tref[0]);
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
        TMatrixDSym hitCov(7);
        hitCov(6, 6) = TMath::Sq(hit_MDC->GetRes()*0.1); // to be adjucted resolution_dl * resolution_dl;
        auto measurement =
            std::make_unique<genfit::WireMeasurement>(hitCoords, hitCov, G4Sol::MG01 + MDCLayerID, hitID, nullptr);
        dynamic_cast<genfit::WireMeasurement*>(measurement.get())->setLeftRightResolution(1);
        dynamic_cast<genfit::WireMeasurement*>(measurement.get())->setMaxDistance(dlmax);

        RecoEvent.ListHits[G4Sol::MG01 + MDCLayerID].emplace_back(measurement.release());
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

  // make container

  int t_r = event.s2fiber->tref[3];

  for(int i = 0; i < (int)event.s2fiber->fiberhit.size(); ++i)
    {
      FiberHitAna* hit_ana = new FiberHitAna(event.s2fiber->fiberhit[i], par.get(), t_r);
      if(!hit_ana->IsValid())
        {
          delete hit_ana;
          continue;
        }
      FiberHitCont[hit_ana->GetDet()][hit_ana->GetLay()].emplace_back(hit_ana);
    }

  for(int i = 0; i < 7; ++i)
    {
      for(int j = 0; j < 3; ++j)
        {
          for(int k = 0; k < (int)FiberHitCont[i][j].size(); ++k)
            {
              LocalHisto.h10[i][j]->Fill(FiberHitCont[i][j][k]->GetFib(), FiberHitCont[i][j][k]->GetTL());
              LocalHisto.h11[i][j]->Fill(FiberHitCont[i][j][k]->GetPos());
            }
        }
    }

  FiberAnalyzer* fiberana = new FiberAnalyzer();
  FiberHitClCont          = fiberana->Clusterize(FiberHitCont);

  int TypeDet[7][3] = {{G4Sol::FiberD1_x, G4Sol::FiberD1_u, G4Sol::FiberD1_v},
                       {G4Sol::FiberD2_x, G4Sol::FiberD2_u, G4Sol::FiberD2_v},
                       {G4Sol::FiberD3_x, G4Sol::FiberD3_u, G4Sol::FiberD3_v},
                       {G4Sol::MiniFiberD1_x1, G4Sol::MiniFiberD1_u1, G4Sol::MiniFiberD1_v1},
                       {G4Sol::MiniFiberD1_x2, G4Sol::MiniFiberD1_u2, G4Sol::MiniFiberD1_v2},
                       {G4Sol::FiberD4_x, G4Sol::FiberD4_u, G4Sol::FiberD4_v},
                       {G4Sol::FiberD5_x, G4Sol::FiberD5_u, G4Sol::FiberD5_v}};

  std::string nameVolDet[7][3] = {{"FiberD1_log_x", "FiberD1_log_u", "FiberD1_log_v"},
                                  {"FiberD2_log_x", "FiberD2_log_u", "FiberD2_log_v"},
                                  {"FiberD3_log_x", "FiberD3_log_u", "FiberD3_log_v"},
                                  {"MiniFiberD1_log_x1", "MiniFiberD1_log_u1", "MiniFiberD1_log_v1"},
                                  {"MiniFiberD1_log_x2", "MiniFiberD1_log_u2", "MiniFiberD1_log_v2"},
                                  {"FiberD4_log_x", "FiberD4_log_u", "FiberD4_log_v"},
                                  {"FiberD5_log_x", "FiberD5_log_u", "FiberD5_log_v"}};

  std::string nameMotherDet[7] = {"FiberD1_log_0",     "FiberD2_log_0", "FiberD3_log_0", "MiniFiberD1_log_0",
                                  "MiniFiberD1_log_0", "FiberD4_log_0", "FiberD5_log_0"};

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

              int hitID = FiberHitClCont[i][j][k]->GetClFib(); // todo check

              TGeoMatrix* g1 =
                  gGeoManager->GetVolume(volumeName.c_str())->GetNode(hitID * 2 + 1)->GetMatrix(); // fiber core
              TGeoMatrix* g2 = gGeoManager->GetVolume("MFLD")
                                   ->GetNode(motherName.c_str())
                                   ->GetVolume()
                                   ->GetNode((volumeName + "_0").c_str())
                                   ->GetMatrix(); // fiber layer
              TGeoMatrix* g3 =
                  gGeoManager->GetVolume("MFLD")->GetNode(motherName.c_str())->GetMatrix(); // fiber station
              TGeoMatrix* g4 = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();     // MFLD
              TGeoHMatrix H1(*g1), H2(*g2), H3(*g3), H4(*g4);
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
              TVector3 o(0., 0., shift[2]), zdir(0., 0., 1.);
              TVector3 fiber_dir(edge2[0] - edge1[0], edge2[1] - edge1[1], edge2[2] - edge1[2]);
              fiber_dir  = fiber_dir.Unit();
              TVector3 u = fiber_dir.Cross(zdir);
              TVector3 v = fiber_dir;
              genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));

              TVectorD hitCoords(1);
              hitCoords(0) = FiberHitClCont[i][j][k]->GetPos()*0.1; //u.Dot(TVector3(shift[0], shift[1], 0));
              TMatrixDSym hitCov(1);
              hitCov(0, 0) = TMath::Sq(FiberHitClCont[i][j][k]->GetRes()*0.1); // mm to cm // to be adjusted resolution_fiber * resolution_fiber;
              auto measurement =
                  std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, TypeDet[i][j], hitID, nullptr);
              dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

              RecoEvent.ListHits[TypeDet[i][j]].emplace_back(measurement.release());
            }
        }
    }

  for(int i = 0; i < 7; ++i)
    {
      for(int j = 0; j < 3; ++j)
        {
          for(int k = 0; k < (int)FiberHitClCont[i][j].size(); ++k)
            {
              LocalHisto.h12[i][j]->Fill(FiberHitClCont[i][j][k]->GetPos());
            }
        }
    }

  for(int i = 0; i < 7; ++i)
    {
      for(int j = 0; j < 3; ++j)
        {
          LocalHisto.h13[i][j]->Fill(FiberHitCont[i][j].size());
          LocalHisto.h14[i][j]->Fill(FiberHitClCont[i][j].size());
          for(int k = 0; k < (int)FiberHitClCont[i][j].size(); ++k)
            {
              LocalHisto.h15[i][j]->Fill(FiberHitClCont[i][j][k]->GetClsize());
            }
        }
    }

  std::vector<std::vector<FiberHitXUV*> > FiberXUVCont = fiberana->FindHit(FiberHitClCont, par.get());
  for(int i = 0; i < 7; ++i)
    {
      for(int j = 0; j < (int)FiberXUVCont[i].size(); ++j)
        {
          LocalHisto.h16[i]->Fill(FiberXUVCont[i][j]->GetPosX(), FiberXUVCont[i][j]->GetPosY());
          LocalHisto.h17_2[i]->Fill(FiberXUVCont[i][j]->GetD());
        }
      LocalHisto.h17[i]->Fill(FiberXUVCont[i].size());
    }

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
              delete FiberHitClCont[i][j].back();
              FiberHitClCont[i][j].pop_back();
            }
        }
    }
  for(int i = 0; i < 7; ++i)
    {
      int num = (int)FiberXUVCont[i].size();
      for(int j = 0; j < num; ++j)
        {
          delete FiberXUVCont[i].back();
          FiberXUVCont[i].pop_back();
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
