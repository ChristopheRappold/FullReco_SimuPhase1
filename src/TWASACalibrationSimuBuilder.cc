#include "TWASACalibrationSimuBuilder.h"

#include "Debug.hh"

#include <list>
#include <map>
#include <set>
#include <vector>

//#define DEBUG_BUILD
//#define DEBUG_BUILD2

using namespace std;

TWASACalibrationSimuBuilder::TWASACalibrationSimuBuilder(const THyphiAttributes& attribut) : TDataBuilder("build_det"), att(attribut)
{
  att._logger->info("TWASACalibrationSimuBuilder::TWASACalibrationSimuBuilder");

  std::vector<std::string> tempName = {
    "HypHI_InSi_log0", "HypHI_InSi_log1", "HypHI_InSi_log2", "HypHI_InSi_log3",
    "TR1_log","TR2_log",
    "Si1_Strip_log_x", "Si1_Strip_log_y", "Si2_Strip_log_x", "Si2_Strip_log_y",
    "SD1_Strip_log_u", "SD1_Strip_log_v", "SD2_Strip_log_u", "SD2_Strip_log_v",
    "SD1pad_Strip_log_u", "SD1pad_Strip_log_v", "SD2pad_Strip_log_u", "SD2pad_Strip_log_v",
    "TO_Counter",
    "FiberD1_Core_log_x", "FiberD1_Core_log_u", "FiberD1_Core_log_v",
    "FiberD2_Core_log_x", "FiberD2_Core_log_u", "FiberD2_Core_log_v",
    "FiberD3_Core_log_x", "FiberD3_Core_log_u", "FiberD3_Core_log_v",
    "MiniFiberD1_Core_log_x1", "MiniFiberD1_Core_log_u1", "MiniFiberD1_Core_log_v1",
    "MiniFiberD1_Core_log_x2", "MiniFiberD1_Core_log_u2", "MiniFiberD1_Core_log_v2",
    "MiniFiberD1_Core_log_x", "MiniFiberD1_Core_log_u", "MiniFiberD1_Core_log_v",
    "MiniFiberD2_Core_log_x", "MiniFiberD2_Core_log_v", "MiniFiberD2_Core_log_u",
    "FiberD4_Core_log_v", "FiberD4_Core_log_u", "FiberD4_Core_log_x",
    "FiberD5_Core_log_x", "FiberD5_Core_log_u", "FiberD5_Core_log_v",
    "PSFE",
    "MG01", "MG02", "MG03", "MG04", "MG05", "MG06", "MG07", "MG08", "MG09", "MG10", "MG11", "MG12",
    "MG13", "MG14", "MG15", "MG16", "MG17",
    "PSCE", "PSBE",
    "CDC_log0", "CDC_log1", "CDC_log2", "CDC_log3", "CDC_log4", "CDC_log5", "CDC_log6", "CDC_log7",
    "CDC_log8", "CDC_log9", "CDC_log10", "CDC_log11", "CDC_log12", "CDC_log13", "CDC_log14",
    "CDH_log",
    "HypHI_TrackFwd_log", "HypHI_TrackFwd_logDummy1", "HypHI_TrackFwd_logDummy2",
    "HypHI_RPC_l_log", "HypHI_RPC_h_log",
    "FMF2_log"};

  for(size_t iName = 0; iName < att.InputPar.nameDet->size(); ++iName)
    {
      std::string nameDetTemp(att.InputPar.nameDet->at(iName));
      for(size_t iTypeDet = 0; iTypeDet < tempName.size(); ++iTypeDet)
        {
          if(nameDetTemp == tempName[iTypeDet])
            {
              orderDetectors.insert(std::make_pair(iName, iTypeDet));
              orderDetName.insert(std::make_pair(iTypeDet, nameDetTemp));
            }
        }
    }
  
  std::string volMDCfirst = gGeoManager->GetVolume("INNER")->GetNode(0)->GetVolume()->GetName();
  if(volMDCfirst == "MD01")      offsetGeoNameID_MDC = 0;
  else if(volMDCfirst == "SOL")  offsetGeoNameID_MDC = 1;
  else                           offsetGeoNameID_MDC = 0;

  auto listNodes = gGeoManager->GetVolume("INNER")->GetNodes();
  int index_lastMDC = -1, index_firstPSCE = -1;
  for(int i=0;i<listNodes->GetEntries();++i)
    {
      std::string tempNameList(listNodes->At(i)->GetName());
      if(tempNameList == "MD17_1") index_lastMDC = i;
      else if(tempNameList == "PSCE_1") index_firstPSCE = i;
      else if(tempNameList == "MDC_1") ++newGeoExp;
      else if(tempNameList == "PSCEall_1") ++newGeoExp;
    }

  if(index_lastMDC == -1 && index_firstPSCE == -1) offsetGeoNameID_PSCE = 0;
  else offsetGeoNameID_PSCE = index_firstPSCE - index_lastMDC + offsetGeoNameID_MDC -1;


  att._logger->info("Builder: geometry type: {} | offsets MDC {} PSCE {}",newGeoExp,offsetGeoNameID_MDC,offsetGeoNameID_PSCE);

}

TWASACalibrationSimuBuilder::~TWASACalibrationSimuBuilder() {}
#ifdef ROOT6
ReturnRes::InfoM TWASACalibrationSimuBuilder::operator()(const TG4Sol_Event& event,
                                                         const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits,
                                                         FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
{
  int result = Exec(event, hits, RecoEvent, OutTree);

  return SoftExit(result);
}
#else
ReturnRes::InfoM TWASACalibrationSimuBuilder::operator()(const TG4Sol_Event& event,
                                                         const std::vector<TClonesArray*>& hits,
                                                         FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
{
  int result = Exec(event, hits, RecoEvent, OutTree);

  return SoftExit(result);
}

ReturnRes::InfoM TWASACalibrationSimuBuilder::operator()(const TG4Sol_Event& event,
                                                         const std::vector<TClonesArray*>& hits,
                                                         FullRecoEvent& RecoEvent, Ana_WasaEvent* OutTree)
{
  return ReturnRes::BuildError;
}

#endif
void TWASACalibrationSimuBuilder::SelectHists()
{
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

  //PSB
  LocalHisto.h76 = AnaHisto->CloneAndRegister(AnaHisto->h76);

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


  LocalHisto.h_stats = AnaHisto->CloneAndRegister(AnaHisto->h_stats);
}

ReturnRes::InfoM TWASACalibrationSimuBuilder::SoftExit(int return_build)
{
  if(return_build == -1)
    {
      att._logger->warn("!> Multiplicity > 2 on Start : event rejected");
      LocalHisto.h_stats->Fill("start M>2", 1);
      return ReturnRes::MultiS2_Start;
    }
  else if(return_build == -2)
    {
      att._logger->warn("!> TDC Timing Start cut : event rejected");
      LocalHisto.h_stats->Fill("start Timing cut", 1);
      return ReturnRes::StartTimingCut;
    }
  else if(return_build == -3)
    {
      att._logger->warn("!> Chamber Hit > 1000 : event rejected");
      LocalHisto.h_stats->Fill("chamber hit>1000", 1);
      return ReturnRes::ChamberHitLimit;
    }
  else if(return_build == -9)
    {
      att._logger->warn("!> No Beam : event rejected");
      LocalHisto.h_stats->Fill("No Beam", 1);
      return ReturnRes::NoBeam;
    }
  else if(return_build != 0)
    {
      att._logger->warn("Error in Build Detector !");
      LocalHisto.h_stats->Fill("Error", 1);
      return ReturnRes::BuildError;
    }
  LocalHisto.h_stats->Fill("start Ok", 1);

  return ReturnRes::Fine;
}
#ifdef ROOT6
int TWASACalibrationSimuBuilder::Exec(const TG4Sol_Event& event, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits,
                                      FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
#else
int TWASACalibrationSimuBuilder::Exec(const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits,
                                      FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
#endif
{
  // int NumFilled = 0;
  //  std::cout<<"TBuidlDetectorLayerPlaneDAF : Proton = "<<event.n_Proton<<endl;

#ifdef DEBUG_BUILD
  auto printW = [](const auto a, const int width, bool lr = true) -> std::string {
    std::stringstream ss;
    ss << std::fixed;
    if(lr)
      ss << std::right;
    else
      ss << std::left;
    ss.fill(' ');    // fill space around displayed #
    ss.width(width); // set  width around displayed #
    ss << a;
    return ss.str();
  };
  auto printFixed = [](const double a, const int decDigits, const int width) -> std::string {
    std::stringstream ss;
    ss << std::fixed << std::right;
    ss.fill(' ');            // fill space around displayed #
    ss.width(width);         // set  width around displayed #
    ss.precision(decDigits); // set # places after decimal
    ss << a;
    return ss.str();
  };

  for(size_t id = 0; id < event.BeamNames.size(); ++id)
    {
      int trackID = event.BeamTrackID[id];
      att._logger->debug("beam : {} # {}", event.BeamNames[id], trackID);

      for(auto& det : hits)
        {
#ifdef ROOT6
          for(auto hit : *det)
            {
              if(hit.TrackID == trackID)
                att._logger->debug("Branch: {} {} {} {} {} {}", printW(det->GetBranchName(), 18, false),
                                   printFixed(hit.HitPosX, 4, 6), printFixed(hit.HitPosY, 4, 6),
                                   printFixed(hit.HitPosZ, 4, 6), hit.LayerID, hit.Pname);
            }
#else
          for(int j = 0; j < det->GetEntries(); ++j)
            {
              TG4Sol_Hit* hit = dynamic_cast<TG4Sol_Hit*>(det->At(j));
              if(hit->TrackID == trackID)
                att._logger->debug("Branch: {} hit:{} {} {} mom:{} {} {} | {} {}", printW(det->GetName(), 18, false),
                                   printFixed(hit->HitPosX, 4, 9), printFixed(hit->HitPosY, 4, 9),
                                   printFixed(hit->HitPosZ, 4, 9), printFixed(hit->MomX, 4, 8),
                                   printFixed(hit->MomY, 4, 8), printFixed(hit->MomZ, 4, 8), printW(hit->LayerID, 2),
                                   printW(hit->Pname, 6));
            }
#endif
        }
    }
  for(size_t id = 0; id < event.DaughterNames.size(); ++id)
    {
      int trackID = event.DaughterTrackID[id];
      att._logger->debug("decayed : {} # {}", event.DaughterNames[id], trackID);

      for(auto& det : hits)
        {
#ifdef ROOT6
          for(auto hit : *det)
            {
              // hit.Print();
              if(hit.TrackID == trackID)
                att._logger->debug("Branch: {} {} {} {} {} {}", det->GetBranchName(), hit.HitPosX, hit.HitPosY,
                                   hit.HitPosZ, hit.LayerID, hit.Pname);
            }
#else
          for(int j = 0; j < det->GetEntries(); ++j)
            {
              // hit.Print();
              TG4Sol_Hit* hit = dynamic_cast<TG4Sol_Hit*>(det->At(j));
              if(hit->TrackID == trackID)
                att._logger->debug("Branch: {} hit:{} {} {} mom:{} {} {} | {} {}", printW(det->GetName(), 18, false),
                                   printFixed(hit->HitPosX, 4, 9), printFixed(hit->HitPosY, 4, 9),
                                   printFixed(hit->HitPosZ, 4, 9), printFixed(hit->MomX, 4, 8),
                                   printFixed(hit->MomY, 4, 8), printFixed(hit->MomZ, 4, 8), printW(hit->LayerID, 2),
                                   printW(hit->Pname, 6));
            }

#endif
        }
    }
#endif

  std::string nameMother(event.MotherName);
  int id_mother = event.MotherTrackID;

  double Mother_totalMomentum = sqrt(pow(event.MotherMomentumAtDecay_X, 2.) + pow(event.MotherMomentumAtDecay_Y, 2.) +
                                     pow(event.MotherMomentumAtDecay_Z, 2.));
  double Mother_Energy        = pow(event.MotherMass, 2.) + pow(Mother_totalMomentum, 2.);

  RecoEvent.Mother_MomE.SetPxPyPzE(event.MotherMomentumAtDecay_X, event.MotherMomentumAtDecay_Y,
                                   event.MotherMomentumAtDecay_Z, Mother_Energy);

  double gamma_factor    = sqrt(1. + pow(Mother_totalMomentum / event.MotherMass, 2.));
  RecoEvent.Hyp_LifeTime = event.DecayTime * 1000. / gamma_factor; // in ps in the CM frame

  RecoEvent.InteractionPoint[0] = event.InteractionPoint_X;
  RecoEvent.InteractionPoint[1] = event.InteractionPoint_Y;
  RecoEvent.InteractionPoint[2] = event.InteractionPoint_Z;

  RecoEvent.PrimVtxRecons.SetXYZ(event.InteractionPoint_X, event.InteractionPoint_Y, event.InteractionPoint_Z);
  RecoEvent.CovMatrix_IP = {0.01,
                              0., 0.01,
                              0.,   0., 0.1}; //temp value

  RecoEvent.DecayVertex[0] = event.DecayVertex_X;
  RecoEvent.DecayVertex[1] = event.DecayVertex_Y;
  RecoEvent.DecayVertex[2] = event.DecayVertex_Z;

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

  double time_res_t0counter = att.t0_timeres; // ns

  for(size_t index = 0; index < event.BeamTrackID.size(); ++index)
    {
      TMcParticle* OutParticle =
          dynamic_cast<TMcParticle*>(OutTree->fMC_Particle->ConstructedAt(OutTree->fMC_Particle->GetEntries()));
      OutParticle->type      = event.BeamNames[index];
      OutParticle->Mc_id     = event.BeamTrackID[index];
      OutParticle->Mother_id = -1;
      OutParticle->Pdg       = pid_fromName(event.BeamNames[index]);
      OutParticle->Charge    = event.BeamCharges[index];
      OutParticle->MomMass.SetXYZM(event.BeamMomentums_X[index], event.BeamMomentums_Y[index],
                                   event.BeamMomentums_Z[index], event.BeamMasses[index]);
      OutParticle->Vtx.SetXYZT(event.InteractionPoint_X, event.InteractionPoint_Y, event.InteractionPoint_Z, 0.);
      OutParticle->Weigth = event.BeamCharges[index] == 0 ? 0. : 1.;
      OutParticle->GeoAcc = 1.;

      if(event.BeamNames[index] == nameMother)
        {
          id_mother = event.BeamTrackID[index];
          continue;
        }
      if(event.BeamCharges[index] == 0)
        continue;
      int TrackID = event.BeamTrackID[index];

      InfoInit tempInit;
      tempInit.pdg = pid_fromName(event.BeamNames[index]);
      tempInit.charge = event.BeamCharges[index];
      tempInit.time = gRandom->Gaus(0., time_res_t0counter);
      tempInit.posX = event.InteractionPoint_X;
      tempInit.posY = event.InteractionPoint_Y;
      tempInit.posZ = event.InteractionPoint_Z;
      tempInit.momX = event.BeamMomentums_X[index];
      tempInit.momY = event.BeamMomentums_Y[index];
      tempInit.momZ = event.BeamMomentums_Z[index];
      tempInit.mass = event.BeamMasses[index];
      RecoEvent.TrackDAFInitSim.insert(std::make_pair(TrackID, tempInit));

      std::vector<std::vector<SimHit> > tempSetSimHit(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackDAFSim.insert(std::make_pair(TrackID, tempSetSimHit));

      //std::vector<InfoPar> tempSetInfo(G4Sol::SIZEOF_G4SOLDETTYPE);
      //RecoEvent.TrackInfo.insert(std::make_pair(TrackID, tempSetInfo));
    }
  for(size_t index = 0; index < event.DaughterTrackID.size(); ++index)
    {
      TMcParticle* OutParticle =
          dynamic_cast<TMcParticle*>(OutTree->fMC_Particle->ConstructedAt(OutTree->fMC_Particle->GetEntries()));
      OutParticle->type      = event.DaughterNames[index];
      OutParticle->Mc_id     = event.DaughterTrackID[index];
      OutParticle->Mother_id = id_mother;
      OutParticle->Pdg       = pid_fromName(event.DaughterNames[index]);
      OutParticle->Charge    = event.DaughterCharges[index];
      OutParticle->MomMass.SetXYZM(event.DaughterMomentums_X[index], event.DaughterMomentums_Y[index],
                                   event.DaughterMomentums_Z[index], event.DaughterMasses[index]);
      OutParticle->Vtx.SetXYZT(event.DecayVertex_X, event.DecayVertex_Y, event.DecayVertex_Z, event.DecayTime);
      OutParticle->Weigth = event.DaughterCharges[index] == 0 ? 0. : 1.;
      OutParticle->GeoAcc = 1.;

      int TrackID = event.DaughterTrackID[index];

      InfoInit tempInit;
      tempInit.pdg = pid_fromName(event.DaughterNames[index]);
      tempInit.charge = event.DaughterCharges[index];
      tempInit.time = gRandom->Gaus(event.DecayTime, time_res_t0counter);
      tempInit.posX = event.DecayVertex_X;
      tempInit.posY = event.DecayVertex_Y;
      tempInit.posZ = event.DecayVertex_Z;
      tempInit.momX = event.DaughterMomentums_X[index];
      tempInit.momY = event.DaughterMomentums_Y[index];
      tempInit.momZ = event.DaughterMomentums_Z[index];
      tempInit.mass = event.DaughterMasses[index];
      RecoEvent.TrackDAFInitSim.insert(std::make_pair(TrackID, tempInit));

      RecoEvent.DaughtersTrackDAFInit.insert(std::make_pair(TrackID, tempInit));

      std::vector<std::vector<SimHit> > tempSetSimHit(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackDAFSim.insert(std::make_pair(TrackID, tempSetSimHit));

      RecoEvent.TrackMother.insert(
          std::make_pair(TrackID, std::make_tuple(event.MotherTrackID, event.DecayVertex_X, event.DecayVertex_Y,
                                                  event.DecayVertex_Z, event.DecayTime)));
    }

  OutTree->Nmc = OutTree->fMC_Particle->GetEntries();

  RecoEvent.ListHits.resize(G4Sol::SIZEOF_G4SOLDETTYPE);
  RecoEvent.ListHitsInfo.resize(G4Sol::SIZEOF_G4SOLDETTYPE);
  RecoEvent.ListHitsToTracks.resize(G4Sol::SIZEOF_G4SOLDETTYPE);
  RecoEvent.OldListHits.resize(G4Sol::SIZEOF_G4SOLDETTYPE);
  //RecoEvent.Si_HitsEnergyLayer.resize(8);

  auto fillOutHit = [](TClonesArray* out, const TG4Sol_Hit& hit, int PDG, double charge, const TVectorD& hitR,
                       int LayerID, int HitID) {
    TMcHit* OutHit  = dynamic_cast<TMcHit*>(out->ConstructedAt(out->GetEntries()));
    OutHit->name    = hit.Pname;
    OutHit->LayerID = LayerID;
    OutHit->HitID   = HitID;
    OutHit->TrackID = hit.TrackID;
    OutHit->MCHit.SetXYZ(hit.HitPosX, hit.HitPosY, hit.HitPosZ);
    OutHit->Hit.SetXYZ(hitR(0), hitR(1), hitR(2));
    OutHit->MC_id  = hit.TrackID;
    OutHit->Pdg    = PDG;
    OutHit->Charge = charge;
    OutHit->MCparticle.SetXYZM(hit.MomX, hit.MomY, hit.MomZ, hit.Mass);
    OutHit->Brho              = 3.10715497 * OutHit->MCparticle.P() / charge;
    OutHit->MagnetInteraction = 0.;
    OutHit->Time              = hit.Time;
    OutHit->Energy            = hit.Energy;
    OutHit->TrackLength       = hit.TrackLength;
    OutHit->HitLength         = TMath::Sqrt(TMath::Sq(hit.ExitPosX-hit.HitPosX)+TMath::Sq(hit.ExitPosY-hit.HitPosY)+TMath::Sq(hit.ExitPosZ-hit.HitPosZ));

    // std::cout<<" Out> LayerID:"<<LayerID<<" "<<HitID<<std::endl;
  };

  for(size_t iDet = 0; iDet < hits.size(); ++iDet)
    {
#ifdef ROOT6
      TTreeReaderArray<TG4Sol_Hit>* tempHits = hits[iDet];
      std::string nameTempBr(tempHits->GetBranchName());
#else
      TClonesArray* tempHits = hits[iDet];
      std::string nameTempBr(tempHits->GetName());
#endif

      //printf("iDet: %d , nameTempBr:  , TypeDet: %d \n", iDet, G4Sol::FiberD5_v, TypeDet);

      auto tempPair                                       = orderDetectors.find(iDet);
      G4Sol::SolDet TypeDet                               = G4Sol::SolDet(tempPair->second);
      std::unique_ptr<genfit::AbsMeasurement> measurement = nullptr;
      MeasurementInfo measinfo;

      //printf("iDet: %d , nameTempBr:  , TypeDet: %d \n", iDet, G4Sol::FiberD5_v, TypeDet);

#ifdef DEBUG_BUILD
      att._logger->debug("iDet # {} {} {}", iDet, nameTempBr, TypeDet);
#endif

      std::set<int> used_fiber_pair;
      // double resolution_wire = 0.01;
      // double resolution_wire_z = 0.1;
      // double resolution_planar = 0.05; // cm
      double resolution_wire   = 1;
      double resolution_wire_z = 10;
      double resolution_planar = 1;    // cm
      double resolution_dl     = 0.02; // cm
      double resolution_fiber  = 0.015;
      double resolution_psce   = 1.1; // 3.8/sqrt(12.)
      double resolution_psce_z = 1.0;
      double time_res_psb      = att.psb_timeres; // ns
      double dE_res_psb        = 0.1;   // in % of dE
      double time_res_psbe     = 0.150; // ns
      double dE_res_psbe       = 0.1;   // in % of dE
      double time_res_psfe     = 0.150; // ns
      double dE_res_psfe       = 0.1;   // in % of dE
      double time_res_fiber    = 0.150; // ns
      double time_res_mdc      = 0.150; // ns
      double time_res          = 0.150; // ns
      if(nameTempBr == "FMF2_log" || nameTempBr == "HypHI_TrackFwd_log")
        {
          // continue;
#ifdef ROOT6
          for(auto it_hit = tempHits->begin(), it_hit_end = tempHits->end(); it_hit != it_hit_end; ++it_hit)
#else
          for(int it_hit = 0; it_hit < tempHits->GetEntries(); ++it_hit)
#endif
            {
#ifdef ROOT6
              auto hit          = *it_hit;
              int indexInBranch = std::distance(tempHits->begin(), it_hit);
#else
              const TG4Sol_Hit& hit = *(dynamic_cast<TG4Sol_Hit*>(tempHits->At(it_hit)));
              int indexInBranch     = it_hit;
#endif
              int TrackID = hit.TrackID;
              int LayerID = hit.LayerID;
#ifdef DEBUG_BUILD
              att._logger->debug(" hit#{} {} {} {}", indexInBranch, hit.Pname, hit.TrackID, hit.LayerID);
#endif
              auto tempTrack = RecoEvent.TrackDAFSim.find(TrackID);
              if(tempTrack == RecoEvent.TrackDAFSim.end())
                continue;

              TVectorD hitCoordsTree(3);
              hitCoordsTree(0) = gRandom->Gaus(hit.HitPosX, resolution_planar);
              hitCoordsTree(1) = gRandom->Gaus(hit.HitPosY, resolution_planar);
              hitCoordsTree(2) = hit.HitPosZ;

              int pdg_code = pid_fromName(hit.Pname);

              auto tempTrackSimLayers = RecoEvent.TrackDAFSim.find(TrackID);
              SimHit tempHitSim;
              tempHitSim.layerID     = LayerID;
              tempHitSim.hitX        = hit.HitPosX;
              tempHitSim.hitY        = hit.HitPosY;
              tempHitSim.hitZ        = hit.HitPosZ;
              tempHitSim.momX        = hit.MomX;
              tempHitSim.momY        = hit.MomY;
              tempHitSim.momZ        = hit.MomZ;
              tempHitSim.pdg         = pdg_code;
              tempHitSim.mass        = hit.Mass;
              tempHitSim.Eloss       = hit.Energy;
              tempHitSim.time        = gRandom->Gaus(hit.Time, time_res);
              tempHitSim.tracklength = hit.TrackLength;
              tempHitSim.hitlength   = TMath::Sqrt(TMath::Sq(hit.ExitPosX-hit.HitPosX)+TMath::Sq(hit.ExitPosY-hit.HitPosY)+TMath::Sq(hit.ExitPosZ-hit.HitPosZ));
              tempTrackSimLayers->second[TypeDet + LayerID].emplace_back(tempHitSim);

              auto PDG_particle = TDatabasePDG::Instance()->GetParticle(pdg_code);
              if(PDG_particle == nullptr)
                {
                  att._logger->error("E> PDG not found !");
                  continue;
                }
              const double charge = PDG_particle->Charge() / 3.;

              if(TypeDet + LayerID >= G4Sol::TrFwd0 && TypeDet + LayerID <= G4Sol::TrFwd2)
                fillOutHit(OutTree->FwdTracker, hit, pdg_code, charge, hitCoordsTree, TypeDet + LayerID, 0.);
              if(TypeDet + LayerID >= G4Sol::FMF2Stop0 && TypeDet + LayerID <= G4Sol::FMF2Stop2)
                fillOutHit(OutTree->FMF2, hit, pdg_code, charge, hitCoordsTree, TypeDet + LayerID, 0.);
            }
        }
      else
        {
#ifdef ROOT6
          for(auto it_hit = tempHits->begin(), it_hit_end = tempHits->end(); it_hit != it_hit_end; ++it_hit)
#else
          for(int it_hit = 0; it_hit < tempHits->GetEntries(); ++it_hit)
#endif
            {
#ifdef ROOT6
              auto hit          = *it_hit;
              int indexInBranch = std::distance(tempHits->begin(), it_hit);
#else
              const TG4Sol_Hit& hit = *(dynamic_cast<TG4Sol_Hit*>(tempHits->At(it_hit)));
              // int indexInBranch     = it_hit;
#endif
              int TrackID = hit.TrackID;
              int LayerID = hit.LayerID;
#ifdef DEBUG_BUILD
              att._logger->debug(" hit#{} {} {} {}", indexInBranch, hit.Pname, hit.TrackID, hit.LayerID);
#endif
              auto tempTrack = RecoEvent.TrackDAFSim.find(TrackID);
              if(tempTrack == RecoEvent.TrackDAFSim.end())
                continue;

              TVectorD hitCoordsTree(3);

              if(Fiber_removefragment_flag && hit.Pname == "He3")
                continue;

/*
              if(IsSilicon(TypeDet))
                {
                  // void simulHitstoSignals(TTreeReaderArray<TG4Sol_Hit>* DetHits,
                  // std::vector<std::tuple<double,size_t>>& HitEnergyLayer)
                  const double EnergyThreshold = 0.001; // MeV
#ifdef DEBUG_BUILD3
                  std::cout << "Silicon \n";
                  std::string tempName = orderDetName.find(TypeDet)->second;
                  std::cout << " name : " << tempName << "\n";
                  std::cout << " LayerID : " << LayerID << "\n";
                  std::cout << " Energy : " << hit.Energy << "\n";
                  std::cout << " PosZ : " << hit.HitPosZ << "\n";
#endif
                  if(hit.Energy < EnergyThreshold)
                    continue;

                  int idSi      = TypeDet - G4Sol::Si1x;
                  auto it_SiHit = RecoEvent.Si_HitsEnergyLayer[idSi].find(LayerID);
                  if(it_SiHit != RecoEvent.Si_HitsEnergyLayer[idSi].end())
                    it_SiHit->second += hit.Energy;
                  else
                    RecoEvent.Si_HitsEnergyLayer[idSi].insert(std::make_pair(LayerID, hit.Energy));
                }

              else if(IsSilicon_SD(TypeDet))
                {
                  // void simulHitstoSignals(TTreeReaderArray<TG4Sol_Hit>* DetHits,
                  // std::vector<std::tuple<double,size_t>>& HitEnergyLayer)
                  const double EnergyThreshold = 0.001; // MeV
#ifdef DEBUG_BUILD
                  std::cout << "Silicon \n";
                  std::string tempName = orderDetName.find(TypeDet)->second;
                  std::cout << " name : " << tempName << "\n";
                  std::cout << " LayerID :" << LayerID << "\n";
                  std::cout << " Energy :" << hit.Energy << "\n";
#endif
                  if(hit.Energy < EnergyThreshold)
                    continue;

                  int idSi      = TypeDet - G4Sol::Si1x_SD;
                  auto it_SiHit = RecoEvent.Si_HitsEnergyLayer[idSi].find(LayerID);
                  if(it_SiHit != RecoEvent.Si_HitsEnergyLayer[idSi].end())
                    it_SiHit->second += hit.Energy;
                  else
                    RecoEvent.Si_HitsEnergyLayer[idSi].insert(std::make_pair(LayerID, hit.Energy));
                }

	      else if(IsSilicon_SD_pad(TypeDet))
                {
                  // void simulHitstoSignals(TTreeReaderArray<TG4Sol_Hit>* DetHits,
                  // std::vector<std::tuple<double,size_t>>& HitEnergyLayer)
                  const double EnergyThreshold = 0.001; // MeV
#ifdef DEBUG_BUILD
                  std::cout << "Silicon pad \n";
                  std::cout << "SiliconID: " << TypeDet << "\n";
                  std::string tempName = orderDetName.find(TypeDet)->second;
                  std::cout << " name : " << tempName << "\n";
                  std::cout << " LayerID :" << LayerID << "\n";
                  std::cout << " Energy :" << hit.Energy << "\n";
#endif
                  if(hit.Energy < EnergyThreshold)
                    continue;

                  int idSi      = TypeDet - G4Sol::Si1x_SD; // id 0 1 2 3 : first and second station [4 5 6 7] was for Hamamatsu geo -> not use anymore / 8 9 10 11 : pad config
                  auto it_SiHit = RecoEvent.Si_HitsEnergyLayer[idSi].find(LayerID);
                  if(it_SiHit != RecoEvent.Si_HitsEnergyLayer[idSi].end())
                    it_SiHit->second += hit.Energy;
                  else
                    RecoEvent.Si_HitsEnergyLayer[idSi].insert(std::make_pair(LayerID, hit.Energy));

                  //if(idSi == 4)
                    //std::cout << "HitPosX : " << hit.HitPosX << "\tHitPosY : " << hit.HitPosY << std::endl;
                }
*/
              if(IsPSCE(TypeDet))
                {
                  //LayerID--;

                  std::array<double,3> shift;
                  TVector3 o,u,v;

                  if(newGeoExp==2)
                    {
#ifdef DEBUG_BUILD2
                      std::cout << "PSC" << std::endl;
                      std::string tmpName = orderDetName.find(TypeDet)->second;
                      std::cout << "name : " << tmpName << std::endl;
                      std::cout << "LayerID: " << LayerID <<" : " <<TypeDet - G4Sol::MG01 + offsetGeoNameID_PSCE + LayerID - 1<<std::endl;
                      std::cout << "HitPosX: " << hit.HitPosX << std::endl;
                      std::cout << "HitPosY: " << hit.HitPosY << std::endl;
                      std::cout << "HitPosZ: " << hit.HitPosZ << std::endl;
                      gGeoManager->GetVolume("PSCEall")->GetNode(LayerID - 1)->Print();
                      gGeoManager->GetVolume("PSCEall")->GetNode(LayerID - 1)->GetMatrix()->Print();
                      gGeoManager->GetVolume("INNER")->GetNode(1)->Print();
                      gGeoManager->GetVolume("INNER")->GetNode(1)->GetMatrix()->Print();
                      gGeoManager->GetVolume("MFLD")->GetNode(0)->Print();
                      gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix()->Print();
                      gGeoManager->GetVolume("WASA")->GetNode(0)->Print();
                      gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix()->Print();

#endif
                      TGeoMatrix* g0 = gGeoManager->GetVolume("PSCEall")->GetNode(LayerID - 1)->GetMatrix(); // PSCEbar
                      TGeoMatrix* g1 = gGeoManager->GetVolume("INNER")->GetNode(1)->GetMatrix(); // PSCEall
                      TGeoMatrix* g2 = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix(); // INNER
                      TGeoMatrix* g3 = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix(); // MFLD
                      TGeoHMatrix H0(*g0), H1(*g1), H2(*g2), H3(*g3);

                      TGeoHMatrix H = H1 * H0;
                      H             = H2 * H;
                      H             = H3 * H;
                      TGeoHMatrix Hinsf1("Hsf1"); // PSCE inner surface
                      Hinsf1.SetDz(-0.4);
                      Hinsf1.SetDx(-27.5);
                      TGeoHMatrix Hinsf2("Hsf2"); // PSCE inner surface
                      Hinsf2.SetDz(-0.4);
                      Hinsf2.SetDx(27.5);
                      TGeoHMatrix       Hsf1      = H * Hinsf1;
                      TGeoHMatrix       Hsf2      = H * Hinsf2;
#ifdef DEBUG_BUILD2
                      H.Print();
                      Hsf1.Print();
                      Hsf2.Print();
#endif
                      double* shift1 = Hsf1.GetTranslation();
                      double* shift2 = Hsf2.GetTranslation();

                      shift[0] = 0.5*(shift1[0]+shift2[0]);
                      shift[1] = 0.5*(shift1[1]+shift2[1]);
                      shift[2] = 0.5*(shift1[2]+shift2[2]);


#ifdef DEBUG_BUILD2
                      TVector3 Bar_dir(shift2[0]-shift1[0],shift2[1]-shift1[1],shift2[2]-shift1[2]);
                      Bar_dir = Bar_dir.Unit();
                      att._logger->debug("center bar : {} {} {} | Bar dir : {} {} {} unit {} {} {}",shift[0],shift[1],shift[2],shift2[0]-shift1[0],shift2[1]-shift1[1],shift2[2]-shift1[2], Bar_dir.X(), Bar_dir.Y(),Bar_dir.Z());
#endif
                      o = TVector3(shift[0], shift[1], shift[2]);
                      TVector3 phidir(shift[0], shift[1], 0), zdir(0., 0., 1.);
                      phidir     = phidir.Unit();
                      u = zdir.Cross(phidir);
                      v = zdir;
                    }
                  else
                    {
#ifdef DEBUG_BUILD2
                      std::cout << "PSC" << std::endl;
                      std::string tmpName = orderDetName.find(TypeDet)->second;
                      std::cout << "name: " << tmpName << std::endl;
                      std::cout << "LayerID: " << LayerID <<" : " <<TypeDet - G4Sol::MG01 + offsetGeoNameID_PSCE + LayerID - 1<<std::endl;
                      std::cout << "HitPosX: " << hit.HitPosX << std::endl;
                      std::cout << "HitPosY: " << hit.HitPosY << std::endl;
                      std::cout << "HitPosZ: " << hit.HitPosZ << std::endl;
                      gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_PSCE + LayerID - 1)->Print();
                      gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_PSCE + LayerID - 1)->GetMatrix()->Print();
                      gGeoManager->GetVolume("MFLD")->GetNode(0)->Print();
                      gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix()->Print();
                      gGeoManager->GetVolume("WASA")->GetNode(0)->Print();
                      gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix()->Print();
#endif
                      TGeoMatrix* g1 = gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_PSCE + LayerID - 1)->GetMatrix(); // PSCE
                      TGeoMatrix* g2 = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix(); // INNER
                      TGeoMatrix* g3 = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix(); // MFLD
                      TGeoHMatrix H1(*g1), H2(*g2), H3(*g3);
                      TGeoHMatrix H = H2 * H1;
                      H             = H3 * H;

#ifdef DEBUG_BUILD2
                      H.Print();
#endif
                      TGeoHMatrix Hsf("Hsf"); // PSCE inner surface
                      Hsf.SetDz(-0.4);
                      H             = H * Hsf;
                      double* shift_temp = H.GetTranslation();
                      shift = {shift_temp[0], shift_temp[1], shift_temp[2]};
                      o = TVector3(shift[0], shift[1], shift[2]);
                      TVector3 phidir(shift[0], shift[1], 0), zdir(0., 0., 1.);
                      phidir     = phidir.Unit();
                      u = zdir.Cross(phidir);
                      v = zdir;
                    }

                  genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));

                  TVectorD hitCoords(2);
                  hitCoords(0) = 0.;
                  hitCoords(1) = gRandom->Gaus(hit.HitPosZ - shift[2], resolution_psce_z);
                  TMatrixDSym hitCov(2);
                  hitCov(0, 0) = resolution_psce * resolution_psce;
                  hitCov(1, 1) = resolution_psce_z * resolution_psce_z;
                  measurement =
                      std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID, nullptr);
                  dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

                  measinfo.SetInfo(-9999., gRandom->Gaus(hit.Time, time_res_psb), gRandom->Gaus(hit.Energy, hit.Energy * dE_res_psb));
                  measinfo.SetPDG(pid_fromName(hit.Pname));

                  hitCoordsTree(0) = hit.HitPosX;
                  hitCoordsTree(1) = hit.HitPosY;
                  hitCoordsTree(2) = hit.HitPosZ;

                  LocalHisto.h76->Fill(LayerID, atan2(shift[1], shift[0]));
                }
              else if(IsPSBE(TypeDet))
                {
#ifdef DEBUG_BUILD2
                  std::cout << "PSBE" << std::endl;
                  std::string tmpName = orderDetName.find(TypeDet)->second;
                  std::cout << "name : " << tmpName << std::endl;
                  std::cout << "LayerID : " << LayerID << std::endl;
                  std::cout << "HitPosX : " << hit.HitPosX << std::endl;
                  std::cout << "HitPosY : " << hit.HitPosY << std::endl;
                  std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
                  gGeoManager->GetVolume("PSB")->GetNode(LayerID - 1)->Print();
                  gGeoManager->GetVolume("PSB")->GetNode(LayerID - 1)->GetMatrix()->Print();
                  gGeoManager->GetVolume("MFLD")->GetNode(0)->Print();
                  gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix()->Print();
                  gGeoManager->GetVolume("WASA")->GetNode(0)->Print();
                  gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix()->Print();
#endif
                  TGeoMatrix* g1   = gGeoManager->GetVolume("PSB")->GetNode(LayerID - 1)->GetMatrix(); // PSCE
                  TGeoMatrix* g1_1 = gGeoManager->GetVolume("MFLD")->GetNode("PSB_1")->GetMatrix();    // PSB box
                  TGeoMatrix* g2   = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix();          // INNER
                  TGeoMatrix* g3   = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();          // MFLD
                  TGeoHMatrix H1(*g1), H1_1(*g1_1), H2(*g2), H3(*g3);
                  TGeoHMatrix H = H1_1 * H1;
//                  H             = H2 * H;

                  H = H3 * H;
#ifdef DEBUG_BUILD2
                  H.Print();
#endif

                  double* shift     = H.GetTranslation();
                  double* local_rot = H.GetRotationMatrix();

                  TVector3 v(local_rot[0], local_rot[3], local_rot[6]);
                  // v is at the left border of the bar -> rotate 4.09 degree to be at the center of the bar
//                  v.RotateZ(-4.09 * TMath::DegToRad());
                  v = v.Unit();
                  TVector3 u(v.Y(), -v.X(), 0.);

                  TGeoTubeSeg* shapePSBE = dynamic_cast<TGeoTubeSeg*>(gGeoManager->GetVolume("PSBE")->GetShape());
                  double mid_r = 0.5*(shapePSBE->GetRmin()+shapePSBE->GetRmax());
                  shift[0] = mid_r*std::cos(v.Phi());
                  shift[1] = mid_r*std::sin(v.Phi());

                  TVector3 o(shift[0], shift[1], shift[2]);

                  //u.SetXYZ(1,0,0);
                  //v.SetXYZ(0,1,0);
                  //o.SetXYZ(hit.HitPosX,hit.HitPosY,hit.HitPosZ);

                  genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));

                  TVectorD hitCoords(2);
                  //hitCoords(0) = gRandom->Uniform(-3.75, 3.75); // phi ! be aware ! not u-dim
                  //hitCoords(1) = gRandom->Uniform(6., 22.);     // r -> v dir
                  hitCoords(0) = 0.; // phi ! be aware ! not u-dim
                  hitCoords(1) = 0.;     // r -> v dir

                  TMatrixDSym hitCov(2);
                  //hitCov(0, 0) = TMath::Sq(2 * hitCoords(1) * TMath::Sin(3.75 * TMath::DegToRad())) / 12.;
                  //hitCov(1, 1) = TMath::Sq(22. - 6.) / 12.;
                  //hitCov(0, 0) = 0.1;
                  //hitCov(1, 1) = 0.1;
                  hitCov(0, 0) = 0.5744*0.5744;
                  hitCov(1, 1) = 2.956*2.956;
                  hitCov(0, 1) = 0.00075487717;
                  hitCov(1, 0) = 0.00075487717;

                  measurement =
                      std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID, nullptr);
                  dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

                  measinfo.SetInfo(-9999., gRandom->Gaus(hit.Time, time_res_psbe), gRandom->Gaus(hit.Energy, hit.Energy * dE_res_psbe));
                  measinfo.SetPDG(pid_fromName(hit.Pname));

                  hitCoordsTree(0) = hit.HitPosX;
                  hitCoordsTree(1) = hit.HitPosY;
                  hitCoordsTree(2) = hit.HitPosZ;
                }
              else if(IsPSFE(TypeDet))
                {
#ifdef DEBUG_BUILD2
                  std::cout << "PSFE" << std::endl;
                  std::string tmpName = orderDetName.find(TypeDet)->second;
                  std::cout << "name : " << tmpName << std::endl;
                  std::cout << "LayerID : " << LayerID << std::endl;
                  std::cout << "HitPosX : " << hit.HitPosX << std::endl;
                  std::cout << "HitPosY : " << hit.HitPosY << std::endl;
                  std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
                  gGeoManager->GetVolume("PSF")->GetNode(LayerID - 1)->Print();
                  gGeoManager->GetVolume("PSF")->GetNode(LayerID - 1)->GetMatrix()->Print();
                  gGeoManager->GetVolume("MFLD")->GetNode(0)->Print();
                  gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix()->Print();
                  gGeoManager->GetVolume("WASA")->GetNode(0)->Print();
                  gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix()->Print();
#endif
                  TGeoMatrix* g1   = gGeoManager->GetVolume("PSF")->GetNode(LayerID - 1)->GetMatrix(); // PSCE
                  TGeoMatrix* g1_1 = gGeoManager->GetVolume("MFLD")->GetNode("PSF_1")->GetMatrix();    // PSB box
                  TGeoMatrix* g2   = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix();          // INNER
                  TGeoMatrix* g3   = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();          // MFLD
                  TGeoHMatrix H1(*g1), H1_1(*g1_1), H2(*g2), H3(*g3);
                  TGeoHMatrix H = H1_1 * H1;
//                  H             = H2 * H;

                  H = H3 * H;
#ifdef DEBUG_BUILD2
                  H.Print();
#endif

                  double* shift     = H.GetTranslation();
                  double* local_rot = H.GetRotationMatrix();

                  TVector3 v(local_rot[0], local_rot[3], local_rot[6]);
                  // v is at the left border of the bar -> rotate 3.75 degree to be at the center of the bar
//                  v.RotateZ(-3.75 * TMath::DegToRad());
                  v = v.Unit();

                  TVector3 u(v.Y(), -v.X(), 0.);

                  TGeoTubeSeg* shapePSFE = dynamic_cast<TGeoTubeSeg*>(gGeoManager->GetVolume("PSFE")->GetShape());
                  double mid_r = 0.5*(shapePSFE->GetRmin()+shapePSFE->GetRmax());
                  shift[0] = mid_r*std::cos(v.Phi());
                  shift[1] = mid_r*std::sin(v.Phi());

                  TVector3 o(shift[0], shift[1], shift[2]);
                  genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));

                  TVectorD hitCoords(2);
                  //hitCoords(0) = gRandom->Uniform(-3.75, 3.75); // phi ! be aware ! not u-dim
                  //hitCoords(1) = gRandom->Uniform(6., 22.);     // r -> v dir
                  hitCoords(0) = 0.; // phi ! be aware ! not u-dim
                  hitCoords(1) = 0.;     // r -> v dir

                  TMatrixDSym hitCov(2);
                  //hitCov(0, 0) = TMath::Sq(2 * hitCoords(1) * TMath::Sin(3.75 * TMath::DegToRad())) / 12.;
                  //hitCov(1, 1) = TMath::Sq(22. - 6.) / 12.;
                  hitCov(0, 0) = 0.6328*0.6328;
                  hitCov(1, 1) = 2.033*2.033;
                  hitCov(0, 1) = 0.00051441990;
                  hitCov(1, 0) = 0.00051441990;

                  measurement =
                      std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID, nullptr);
                  dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

                  measinfo.SetInfo(-9999., gRandom->Gaus(hit.Time, time_res_psfe), gRandom->Gaus(hit.Energy, hit.Energy * dE_res_psfe));
                  measinfo.SetPDG(pid_fromName(hit.Pname));

                  hitCoordsTree(0) = hit.HitPosX;
                  hitCoordsTree(1) = hit.HitPosY;
                  hitCoordsTree(2) = hit.HitPosZ;
                }
              else if(IsFiber(TypeDet))
                {
                  //if(IsFiberU_Vetoed(TypeDet))
                  //  continue;

#ifdef DEBUG_BUILD2
                  std::cout << "fiber" << std::endl;
                  std::string tmpName = orderDetName.find(TypeDet)->second;
                  std::cout << "name : " << tmpName << std::endl;
                  std::cout << "LayerID (hitID in Simu) : " << LayerID << std::endl;
                  std::cout << "HitID   (hitID in Data) : " << LayerID/2 << std::endl;
                  std::cout << "HitPosX : " << hit.HitPosX << std::endl;
                  std::cout << "HitPosY : " << hit.HitPosY << std::endl;
                  std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
#endif

                  int i_fiber = -1;
                  int i_layer = -1;
                  switch(TypeDet)
                    {
                      case G4Sol::FiberD1_x:     { i_fiber = 0; i_layer = 0; } ; break;
                      case G4Sol::FiberD1_u:     { i_fiber = 0; i_layer = 1; } ; break;
                      case G4Sol::FiberD1_v:     { i_fiber = 0; i_layer = 2; } ; break;
                      case G4Sol::FiberD2_x:     { i_fiber = 1; i_layer = 0; } ; break;
                      case G4Sol::FiberD2_u:     { i_fiber = 1; i_layer = 1; } ; break;
                      case G4Sol::FiberD2_v:     { i_fiber = 1; i_layer = 2; } ; break;
                      case G4Sol::FiberD3_x:     { i_fiber = 2; i_layer = 0; } ; break;
                      case G4Sol::FiberD3_u:     { i_fiber = 2; i_layer = 1; } ; break;
                      case G4Sol::FiberD3_v:     { i_fiber = 2; i_layer = 2; } ; break;
                      case G4Sol::MiniFiberD1_x: { i_fiber = 3; i_layer = 0; } ; break;
                      case G4Sol::MiniFiberD1_u: { i_fiber = 3; i_layer = 1; } ; break;
                      case G4Sol::MiniFiberD1_v: { i_fiber = 3; i_layer = 2; } ; break;
                      case G4Sol::MiniFiberD2_x: { i_fiber = 4; i_layer = 0; } ; break;
                      case G4Sol::MiniFiberD2_v: { i_fiber = 4; i_layer = 1; } ; break;
                      case G4Sol::MiniFiberD2_u: { i_fiber = 4; i_layer = 2; } ; break;
                      case G4Sol::FiberD4_v:     { i_fiber = 5; i_layer = 0; } ; break;
                      case G4Sol::FiberD4_u:     { i_fiber = 5; i_layer = 1; } ; break;
                      case G4Sol::FiberD4_x:     { i_fiber = 5; i_layer = 2; } ; break;
                      case G4Sol::FiberD5_x:     { i_fiber = 6; i_layer = 0; } ; break;
                      case G4Sol::FiberD5_u:     { i_fiber = 6; i_layer = 1; } ; break;
                      case G4Sol::FiberD5_v:     { i_fiber = 6; i_layer = 2; } ; break;

                      default:  std::cerr << "something wrong" << std::endl; break;
                    }


                  int pdg_code = pid_fromName(hit.Pname);
                  if(pdg_code == 0)
                    att._logger->debug("!> Builder : pdg_code = 0 ! {}", hit.Pname);

                  FiberHitAna* hit_ana = new FiberHitAna(i_fiber, i_layer, LayerID/2,
                                  gRandom->Gaus(hit.Time, time_res_fiber), hit.Energy, pdg_code, TrackID, att.par.get());
                  FiberHitCont[hit_ana->GetDet()][hit_ana->GetLay()].emplace_back(hit_ana);

                  auto tempTrackSimFibers = RecoEvent.TrackDAFSim.find(TrackID);
                  SimHit tempHitSim;
                  tempHitSim.layerID      = LayerID;
                  tempHitSim.hitX         = hit.HitPosX;
                  tempHitSim.hitY         = hit.HitPosY;
                  tempHitSim.hitZ         = hit.HitPosZ;
                  tempHitSim.momX         = hit.MomX;
                  tempHitSim.momY         = hit.MomY;
                  tempHitSim.momZ         = hit.MomZ;
                  tempHitSim.pdg          = pdg_code;
                  tempHitSim.mass         = hit.Mass;
                  tempHitSim.Eloss        = hit.Energy;
                  tempHitSim.time         = hit.Time;
                  tempHitSim.tracklength  = hit.TrackLength;
                  tempHitSim.hitlength    = TMath::Sqrt(TMath::Sq(hit.ExitPosX-hit.HitPosX)+TMath::Sq(hit.ExitPosY-hit.HitPosY)+TMath::Sq(hit.ExitPosZ-hit.HitPosZ));
                  tempTrackSimFibers->second[TypeDet].emplace_back(tempHitSim);

                  auto PDG_particle = TDatabasePDG::Instance()->GetParticle(pdg_code);
                  if(PDG_particle == nullptr)
                    {
                      att._logger->error("E> PDG not found !");
                      continue;
                    }
                  const double charge = PDG_particle->Charge() / 3.;

                  hitCoordsTree(0) = hit.HitPosX;
                  hitCoordsTree(1) = hit.HitPosY;
                  hitCoordsTree(2) = hit.HitPosZ;
                  
                  fillOutHit(OutTree->Fiber, hit, pdg_code, charge, hitCoordsTree, TypeDet, LayerID);

                  continue;
                }
              else if(IsWire(TypeDet))
                {

                  double* shift = nullptr;
                  double* edge1 = nullptr;
                  double* edge2 = nullptr;

                  if(newGeoExp==2)
                    {
#ifdef DEBUG_BUILD2
                      std::cout << "wire" << std::endl;
                      std::string tmpName = orderDetName.find(TypeDet)->second;
                      std::cout << "name : " << tmpName << std::endl;
                      std::cout << "LayerID : " << LayerID <<" : "<< TypeDet - G4Sol::MG01 + 1 <<std::endl;
                      std::cout << "HitPosX : " << hit.HitPosX << std::endl;
                      std::cout << "HitPosY : " << hit.HitPosY << std::endl;
                      std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
                      gGeoManager->GetVolume("MDC")->GetNode(TypeDet - G4Sol::MG01)->GetVolume()->GetNode(LayerID - 1)->GetVolume()->GetNode(0)->Print();
                      gGeoManager->GetVolume("MDC")->GetNode(TypeDet - G4Sol::MG01)->GetVolume()->GetNode(LayerID - 1)->GetVolume()->GetNode(0)->GetMatrix()->Print();

                      gGeoManager->GetVolume("MDC")->GetNode(TypeDet - G4Sol::MG01)->GetVolume()->GetNode(LayerID - 1)->Print();
                      gGeoManager->GetVolume("MDC")->GetNode(TypeDet - G4Sol::MG01)->GetVolume()->GetNode(LayerID - 1)->GetMatrix()->Print();
                      gGeoManager->GetVolume("MDC")->GetNode(TypeDet - G4Sol::MG01)->Print();
                      gGeoManager->GetVolume("MDC")->GetNode(TypeDet - G4Sol::MG01)->GetMatrix()->Print();
                      gGeoManager->GetVolume("INNER")->GetNode(0)->Print();
                      gGeoManager->GetVolume("INNER")->GetNode(0)->GetMatrix()->Print();
                      gGeoManager->GetVolume("MFLD")->GetNode(0)->Print();
                      gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix()->Print();
                      std::cout<<"WASA geo\n";
                      gGeoManager->GetVolume("WASA")->GetNode(0)->Print();
                      gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix()->Print();
#endif
                      TGeoMatrix* g0 = gGeoManager->GetVolume("MDC")->GetNode(TypeDet - G4Sol::MG01)->GetVolume()->GetNode(LayerID - 1)->GetMatrix(); // ME, MG
                      TGeoShape* tempShape =
                        gGeoManager->GetVolume("MDC")->GetNode(TypeDet - G4Sol::MG01)->GetVolume()->GetShape();
                      TGeoMatrix* g1 =
                        gGeoManager->GetVolume("MDC")->GetNode(TypeDet - G4Sol::MG01)->GetMatrix(); // MD
                      TGeoMatrix* g2 = gGeoManager->GetVolume("INNER")->GetNode(0)->GetMatrix(); // MDC
                      TGeoMatrix* g3 = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix(); // INNER
                      TGeoMatrix* g4 = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix(); // MFLD
                      TGeoHMatrix H0(*g0), H1(*g1), H2(*g2), H3(*g3), H4(*g4);
                      TGeoHMatrix H = H1 * H0;
                      H             = H2 * H;
                      H             = H3 * H;
                      H             = H4 * H;
                      shift = H.GetTranslation();
                      TGeoHMatrix w1("w1");
                      TGeoHMatrix w2("w2");
                      Double_t minZ, maxZ;
                      tempShape->GetAxisRange(3, minZ, maxZ);
                      w1.SetDz(minZ);
                      w2.SetDz(maxZ);
                      TGeoHMatrix Hw1 = H * w1;
                      TGeoHMatrix Hw2 = H * w2;
#ifdef DEBUG_BUILD2
                      H.Print();
                      Hw1.Print();
                      Hw2.Print();
#endif
                      edge1 = Hw1.GetTranslation();
                      edge2 = Hw2.GetTranslation();
                    }
                  else
                    {
#ifdef DEBUG_BUILD2
                      std::cout << "wire" << std::endl;
                      std::string tmpName = orderDetName.find(TypeDet)->second;
                      std::cout << "name : " << tmpName << std::endl;
                      std::cout << "LayerID : " << LayerID <<" : "<< TypeDet - G4Sol::MG01 + 1 <<std::endl;
                      std::cout << "HitPosX : " << hit.HitPosX << std::endl;
                      std::cout << "HitPosY : " << hit.HitPosY << std::endl;
                      std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
                      gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_MDC)->GetVolume()->GetNode(LayerID - 1)->Print();
                      gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_MDC)->GetVolume()->GetNode(LayerID - 1)->GetMatrix()->Print();
                      gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_MDC)->Print();
                      gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_MDC)->GetMatrix()->Print();
                      gGeoManager->GetVolume("MFLD")->GetNode(0)->Print();
                      gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix()->Print();
                      gGeoManager->GetVolume("WASA")->GetNode(0)->Print();
                      gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix()->Print();
#endif
                      TGeoMatrix* g1 = gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_MDC)->GetVolume()->GetNode(LayerID - 1)->GetMatrix(); // ME, MG
                      TGeoShape* tempShape =
                        gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_MDC)->GetVolume()->GetShape();
                      TGeoMatrix* g2 =
                        gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_MDC)->GetMatrix(); // MD
                      TGeoMatrix* g3 = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix();             // INNER
                      TGeoMatrix* g4 = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();             // MFLD
                      TGeoHMatrix H1(*g1), H2(*g2), H3(*g3), H4(*g4);
                      TGeoHMatrix H = H2 * H1;
                      H             = H3 * H;
                      H             = H4 * H;
                      shift = H.GetTranslation();
                      TGeoHMatrix w1("w1");
                      TGeoHMatrix w2("w2");
                      Double_t minZ, maxZ;
                      tempShape->GetAxisRange(3, minZ, maxZ);
                      w1.SetDz(minZ);
                      w2.SetDz(maxZ);
                      TGeoHMatrix Hw1 = H * w1;
                      TGeoHMatrix Hw2 = H * w2;
#ifdef DEBUG_BUILD2
                      H.Print();
                      Hw1.Print();
                      Hw2.Print();
#endif
                      edge1 = Hw1.GetTranslation();
                      edge2 = Hw2.GetTranslation();
                    }

                  TVector3 x1(shift[0], shift[1], shift[2]);
                  TVector3 p1(edge2[0] - edge1[0], edge2[1] - edge1[1], edge2[2] - edge1[2]);
                  TVector3 x2(hit.HitPosX, hit.HitPosY, hit.HitPosZ);
                  TVector3 p2(hit.MomX, hit.MomY, hit.MomZ);
                  //double dl2 = CloseDist(x1, x2, p1, p2);
                  TVector3 ClosestPointWire, ClosestPointTrack;
                  double dl  = closestDistanceApproach(x1, x2, p1, p2,ClosestPointWire,ClosestPointTrack);
#ifdef DEBUG_BUILD2
                  //att._logger->debug("Wire Closest distance : dl {} | dl_2 {}",dl, dl2);
                  att._logger->debug("Wire Closest distance : dl {}",dl);

                  att._logger->debug("geometry wire : o {} {} {}, wire_dir {} {} {}",x1.X(),x1.Y(),x1.Z(), p1.X(),p1.Y(),p1.Z());
                  att._logger->debug("sim hit : {} {} {}, mom {} {} {}",x2.X(),x2.Y(),x2.Z(), p2.X(),p2.Y(),p2.Z());
                  att._logger->debug("mid hit : {} {} {}",0.5*(hit.HitPosX+hit.ExitPosX),0.5*(hit.HitPosY+hit.ExitPosY),0.5*(hit.HitPosZ+hit.ExitPosZ));
                  att._logger->debug("closest distance :sim vs wire {} / on fiber {} {} {} / on track {} {} {}",dl, ClosestPointWire.X(),ClosestPointWire.Y(),ClosestPointWire.Z(),ClosestPointTrack.X(),ClosestPointTrack.Y(),ClosestPointTrack.Z());
                  att._logger->debug("Diff {} {} {} | {}",0.5*(hit.HitPosX+hit.ExitPosX)-ClosestPointTrack.X(),0.5*(hit.HitPosY+hit.ExitPosY)-ClosestPointTrack.Y(),0.5*(hit.HitPosZ+hit.ExitPosZ)-ClosestPointTrack.Z(), TMath::Sqrt(TMath::Sq(0.5*(hit.HitPosX+hit.ExitPosX)-ClosestPointTrack.X())+TMath::Sq(0.5*(hit.HitPosY+hit.ExitPosY)-ClosestPointTrack.Y())+TMath::Sq(0.5*(hit.HitPosZ+hit.ExitPosZ)-ClosestPointTrack.Z())));
#endif
                  double dlmax = 0;
                  switch(TypeDet - G4Sol::MG01 + 1)
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
                  double temp_dl = gRandom->Gaus(dl, resolution_dl);
                  bool doneRand  = false;
                  while(doneRand)
                    {
                      if(temp_dl < 0 || temp_dl > dlmax)
                        temp_dl = gRandom->Gaus(dl, resolution_dl);
                      else
                        doneRand = true;
                    }
                  // if(temp_dl<0)     dl = 0;
                  // if(temp_dl>dlmax) dl = dlmax;

                  TVectorD hitCoords(7);
                  hitCoords(0) = edge1[0];
                  hitCoords(1) = edge1[1];
                  hitCoords(2) = edge1[2];
                  hitCoords(3) = edge2[0];
                  hitCoords(4) = edge2[1];
                  hitCoords(5) = edge2[2];
                  hitCoords(6) = temp_dl;
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
                  hitCov(6, 6) = resolution_dl * resolution_dl;
                  measurement =
                      std::make_unique<genfit::WireMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID, nullptr);
                  dynamic_cast<genfit::WireMeasurement*>(measurement.get())->setLeftRightResolution(0);
                  dynamic_cast<genfit::WireMeasurement*>(measurement.get())->setMaxDistance(dlmax);

                  measinfo.SetInfo(-9999., gRandom->Gaus(hit.Time, time_res_mdc), hit.Energy); //Change check TOT?
                  measinfo.SetPDG(pid_fromName(hit.Pname));

                  hitCoordsTree(0) = ClosestPointTrack.X();
                  hitCoordsTree(1) = ClosestPointTrack.Y();
                  hitCoordsTree(2) = ClosestPointTrack.Z();

                  //std::cout << "MDCHit of track: " << TrackID << "\n";
                  //std::cout << "MDC Layer: " << TypeDet-G4Sol::MG01+1 << " hit: " << LayerID << " TrackID: " << TrackID << "\n";
                }
              else
                {
                  continue;
                  // std::cout << "else"  << std::endl;
                  TVectorD hitCoords(3);
                  // std::cout << "HitPosX : " << hit.HitPosX << std::endl;
                  // std::cout << "HitPosY : " << hit.HitPosX << std::endl;

                  hitCoords(0) = gRandom->Gaus(hit.HitPosX, resolution_wire);
                  hitCoords(1) = gRandom->Gaus(hit.HitPosY, resolution_wire);
                  hitCoords(2) = gRandom->Gaus(hit.HitPosZ, resolution_wire_z);

                  TMatrixDSym hitCov(3);
                  hitCov(0, 0) = resolution_wire   * resolution_wire;
                  hitCov(1, 1) = resolution_wire   * resolution_wire;
                  hitCov(2, 2) = resolution_wire_z * resolution_wire_z;

                  measurement = std::make_unique<genfit::ProlateSpacepointMeasurement>(hitCoords, hitCov, int(TypeDet),
                                                                                       LayerID, nullptr);
                  const TVector3 WireDir(0., 0., 1.);
                  dynamic_cast<genfit::ProlateSpacepointMeasurement*>(measurement.get())
                      ->setLargestErrorDirection(WireDir);

                  hitCoordsTree(0) = hitCoords(0);
                  hitCoordsTree(1) = hitCoords(1);
                  hitCoordsTree(2) = hitCoords(2);
                }

              RecoEvent.ListHits[TypeDet].emplace_back(measurement.release());
              RecoEvent.ListHitsInfo[TypeDet].emplace_back(measinfo);
              RecoEvent.ListHitsToTracks[TypeDet].emplace_back(TrackID);

              int pdg_code = pid_fromName(hit.Pname);
              if(pdg_code == 0)
                att._logger->debug("!> Builder : pdg_code = 0 ! {}", hit.Pname);

              auto tempTrackSim = RecoEvent.TrackDAFSim.find(TrackID);
              SimHit tempHitSim;
              tempHitSim.layerID     = LayerID;
              tempHitSim.hitX        = hit.HitPosX;
              tempHitSim.hitY        = hit.HitPosY;
              tempHitSim.hitZ        = hit.HitPosZ;
              tempHitSim.momX        = hit.MomX;
              tempHitSim.momY        = hit.MomY;
              tempHitSim.momZ        = hit.MomZ;
              tempHitSim.pdg         = pdg_code;
              tempHitSim.mass        = hit.Mass;
              tempHitSim.Eloss       = hit.Energy;
              tempHitSim.time        = hit.Time;
              tempHitSim.tracklength = hit.TrackLength;
              tempHitSim.hitlength    = TMath::Sqrt(TMath::Sq(hit.ExitPosX-hit.HitPosX)+TMath::Sq(hit.ExitPosY-hit.HitPosY)+TMath::Sq(hit.ExitPosZ-hit.HitPosZ));
              tempTrackSim->second[TypeDet].emplace_back(tempHitSim);

              auto PDG_particle = TDatabasePDG::Instance()->GetParticle(pdg_code);
              if(PDG_particle == nullptr)
                {
                  att._logger->error("E> PDG not found !");
                  continue;
                }
              const double charge = PDG_particle->Charge() / 3.;

              if(TypeDet >= G4Sol::InSi0 && TypeDet <= G4Sol::InSi3)
                fillOutHit(OutTree->InSi, hit, pdg_code, charge, hitCoordsTree, TypeDet, LayerID);

              if(TypeDet >= G4Sol::TR1 && TypeDet <= G4Sol::TR2)
                fillOutHit(OutTree->TR, hit, pdg_code, charge, hitCoordsTree, TypeDet, LayerID);

              if(TypeDet >= G4Sol::CDC_layer0 && TypeDet <= G4Sol::CDC_layer14)
                fillOutHit(OutTree->CDC, hit, pdg_code, charge, hitCoordsTree, TypeDet, LayerID);

              if(TypeDet >= G4Sol::MG01 && TypeDet <= G4Sol::MG17)
                fillOutHit(OutTree->CDC, hit, pdg_code, charge, hitCoordsTree, TypeDet, LayerID);

              if(TypeDet == G4Sol::CDHBar)
                fillOutHit(OutTree->CDH, hit, pdg_code, charge, hitCoordsTree, TypeDet, LayerID);

              if(TypeDet == G4Sol::RPC_l || TypeDet == G4Sol::RPC_h)
                fillOutHit(OutTree->RPC, hit, pdg_code, charge, hitCoordsTree, TypeDet, LayerID);

              if(TypeDet == G4Sol::PSFE)
                fillOutHit(OutTree->PSFE, hit, pdg_code, charge, hitCoordsTree, TypeDet, LayerID);

              if(TypeDet == G4Sol::PSCE)
                fillOutHit(OutTree->PSCE, hit, pdg_code, charge, hitCoordsTree, TypeDet, LayerID);

              if(TypeDet == G4Sol::PSBE)
                fillOutHit(OutTree->PSBE, hit, pdg_code, charge, hitCoordsTree, TypeDet, LayerID);
            }
        }
    }

/*
  TString tmp_pionminus = TDatabasePDG::Instance()->GetParticle(-211)->GetName();
  TString tmp_pionplus  = TDatabasePDG::Instance()->GetParticle( 211)->GetName();
  TString tmp_proton    = TDatabasePDG::Instance()->GetParticle(2212)->GetName();
  TString tmp_kaonplus  = TDatabasePDG::Instance()->GetParticle(321)->GetName();

  std::cout << "Pi- name: " << tmp_pionminus << "\n";
  std::cout << "Pi+ name: " << tmp_pionplus << "\n";
  std::cout << "Proton name: " << tmp_proton << "\n";
  std::cout << "Kaon name: " << tmp_kaonplus << "\n";
*/

  // Fiber ana //////////////////////////////////
  FiberAnalyzer* fiberana = new FiberAnalyzer();
  FiberHitClCont          = fiberana->Clusterize(FiberHitCont, att.WF_perfect);

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

              int hitID = FiberHitClCont[i][j][k]->GetClFib();

#ifdef DEBUG_BUILD2
              std::cout << "Fiber:" << i << "   Layer: " << j << "   HitID: " << hitID << "\n"; 
#endif

              TGeoMatrix* g1 =
                  gGeoManager->GetVolume(volumeName.c_str())->GetNode( (fiberName + Form("_%d", hitID*2)).c_str() )->GetMatrix(); // fiber core

              TGeoMatrix* g1_pair =
                  gGeoManager->GetVolume(volumeName.c_str())->GetNode( (fiberName + Form("_%d", hitID*2 +1)).c_str() )->GetMatrix(); // fiber pair core

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
              // double* shift = H.GetTranslation();
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

              MeasurementInfo measinfo(FiberHitClCont[i][j][k]->GetTOT(), FiberHitClCont[i][j][k]->GetTime(), FiberHitClCont[i][j][k]->GetdE());
              measinfo.SetPDG(FiberHitClCont[i][j][k]->GetPDG());
              RecoEvent.ListHitsInfo[TypeDet[i][j]].emplace_back(measinfo);

              RecoEvent.ListHitsToTracks[TypeDet[i][j]].emplace_back(FiberHitClCont[i][j][k]->GetSimTrackID());

              //if(i == 3 || i == 4)
                //std::cout << "FiberHit of track: " << FiberHitClCont[i][j][k]->GetSimTrackID() << "\n";

              //std::cout << "FiberCl Det: " << i << " layer: " << j << " hit: " << hitID << " TrackID: " << FiberHitClCont[i][j][k]->GetSimTrackID() << "\n";
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

  for(int i=0; i<7; ++i){
    for(int j=0; j<3; ++j){
      for(int k=0; k<(int)FiberHitClCont[i][j].size(); ++k){
        FiberHitAna *hit = FiberHitClCont[i][j][k];
        LocalHisto.hfiber_13_0[i][j]->Fill(hit->GetTL());
        LocalHisto.hfiber_13_1[i][j]->Fill(hit->GetTime());
        LocalHisto.hfiber_13_2[i][j]->Fill(hit->GetTime());
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


  // Fiber Track Analysis
  std::map< std::string, std::vector<FiberTrackAna*> > FiberTrackCont
                  = fiberana->FiberTracking(FiberHitClCont, att.par.get(), RecoEvent.ListHits[G4Sol::PSCE]);

  RecoEvent.FiberTrackCont = FiberTrackCont;

  int num_combi_mft12 = 1;
  for(int i=3; i<5; ++i)
    for(int j=0; j<3; ++j)
      num_combi_mft12 *= ( (int)FiberHitClCont[i][j].size() + 1 );

  if(att.par->flag_mft12_combi && num_combi_mft12<att.par->cut_mft12_combi){
    LocalHisto.hfiber_1_1->Fill(num_combi_mft12);
    LocalHisto.hfiber_1_2->Fill((double)num_combi_mft12*1e-3);
    LocalHisto.hfiber_1_3->Fill((double)num_combi_mft12*1e-3);
    LocalHisto.hfiber_1_4->Fill((double)num_combi_mft12*1e-6);
  }

  LocalHisto.hfiber_1_5->Fill((double)num_combi_mft12*1e-3);
  LocalHisto.hfiber_1_6->Fill((double)num_combi_mft12*1e-3);
  LocalHisto.hfiber_1_7->Fill((double)num_combi_mft12*1e-6);
  if(att.par->flag_debug) std::cout << "- num_combi_mft12 : " << num_combi_mft12 << std::endl;

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

  LocalHisto.hfiber_3_0->Fill(FiberTrackCont["mft12"].size());

  int num_combi_dft12 = 1;
  for(int i=5; i<7; ++i)
    for(int j=0; j<3; ++j)
      num_combi_dft12 *= ( (int)FiberHitClCont[i][j].size() + 1 );

  if(att.par->flag_dft12_combi && num_combi_dft12<att.par->cut_dft12_combi){
    LocalHisto.hfiber_5_1->Fill(num_combi_dft12);
    LocalHisto.hfiber_5_2->Fill((double)num_combi_dft12*1e-3);
    LocalHisto.hfiber_5_3->Fill((double)num_combi_dft12*1e-3);
    LocalHisto.hfiber_5_4->Fill((double)num_combi_dft12*1e-6);
  }

  LocalHisto.hfiber_5_5->Fill((double)num_combi_dft12*1e-3);
  LocalHisto.hfiber_5_6->Fill((double)num_combi_dft12*1e-3);
  LocalHisto.hfiber_5_7->Fill((double)num_combi_dft12*1e-6);
  if(att.par->flag_debug) std::cout << "- num_combi_dft12 : " << num_combi_dft12 << std::endl;

  for(size_t i=0; i<FiberTrackCont["dft12"].size(); ++i){
    FiberTrackAna *track = FiberTrackCont["dft12"][i];
    LocalHisto.hfiber_4_2_3->Fill(track->GetXdet(), track->GetYdet());
    LocalHisto.hfiber_4_3_3->Fill(track->GetXdet(), track->GetA()*1000);
    LocalHisto.hfiber_4_4_3->Fill(track->GetYdet(), track->GetB()*1000);
    LocalHisto.hfiber_4_5_3->Fill(track->GetTOT());
  }

  LocalHisto.h18_3_1->Fill(FiberTrackCont["dft12"].size());
  if(FiberTrackCont["dft12"].size()>0){
    LocalHisto.h18_3_2->Fill(FiberTrackCont["dft12"][0]->GetChi2());
    LocalHisto.h18_3_3->Fill(FiberTrackCont["dft12"][0]->GetXtgt(), FiberTrackCont["dft12"][0]->GetA()*1000);
    LocalHisto.h18_3_4->Fill(FiberTrackCont["dft12"][0]->GetYtgt(), FiberTrackCont["dft12"][0]->GetB()*1000);
    LocalHisto.h18_3_5->Fill(FiberTrackCont["dft12"][0]->GetXtgt());
    LocalHisto.h18_3_6->Fill(FiberTrackCont["dft12"][0]->GetYtgt());
    LocalHisto.h18_3_7->Fill(FiberTrackCont["dft12"][0]->GetA()*1000);
    LocalHisto.h18_3_8->Fill(FiberTrackCont["dft12"][0]->GetB()*1000);
  }



  OutTree->Field       = att.Field_Strength;

  OutTree->NInSi       = OutTree->InSi->GetEntries();
  OutTree->NTr         = OutTree->TR->GetEntries();
  OutTree->NFiber      = OutTree->Fiber->GetEntries();
  OutTree->NCdc        = OutTree->CDC->GetEntries();
  OutTree->NCdh        = OutTree->CDH->GetEntries();
  OutTree->NRpc        = OutTree->RPC->GetEntries();
  OutTree->NFwdtracker = OutTree->FwdTracker->GetEntries();
  OutTree->NFmf2       = OutTree->FMF2->GetEntries();
  OutTree->NPsbe       = OutTree->PSBE->GetEntries();
  OutTree->NPsfe       = OutTree->PSFE->GetEntries();
  OutTree->NPsce       = OutTree->PSCE->GetEntries();

#ifdef DEBUG_BUILD
  att._logger->debug("done !");
#endif

#ifdef DEBUG_BUILD
  /*
  att._logger->debug(" DAF Hit :");
  for(auto track : RecoEvent.TrackDAF)
    {
      att._logger->debug("TrackID # {} hit_id [", track.first);
      std::vector<std::stringstream> s1(track.second.size() / 8 + 1);
      std::vector<std::stringstream> s2(track.second.size() / 8 + 1);
      std::vector<std::stringstream> s3(track.second.size() / 8 + 1);
      for(size_t i = 0; i < track.second.size(); ++i)
        {
          s1[i / 8] << printW(G4Sol::nameLiteralDet.begin()[i], 14) << ", ";
          s2[i / 8] << printW(i, 14) << ", ";
          s3[i / 8] << printW(track.second[i], 14) << ", ";
        }
      for(size_t i = 0; i < s1.size(); ++i)
        {
          att._logger->debug("idDet:{}", s1[i].str());
          att._logger->debug("stat :{}", s3[i].str());
        }

      att._logger->debug("] ");
    }

  for(auto track : RecoEvent.TrackInfo)
    {
      att._logger->debug("TrackID #{} PID [", track.first);
      std::vector<std::stringstream> s1(track.second.size() / 8 + 1);
      std::vector<std::stringstream> s2(track.second.size() / 8 + 1);
      std::vector<std::stringstream> s3(track.second.size() / 8 + 1);
      for(size_t i = 0; i < track.second.size(); ++i)
        {
          s1[i / 8] << printW(G4Sol::nameLiteralDet.begin()[i], 14) << ", ";
          s2[i / 8] << printW(i, 14) << ", ";
          s3[i / 8] << printW(track.second[i].pdg, 14) << ", ";
        }
      for(size_t i = 0; i < s1.size(); ++i)
        {
          att._logger->debug("idDet:{}", s1[i].str());
          att._logger->debug("stat :{}", s3[i].str());
        }

      att._logger->debug("] ");
    }
*/
  for(size_t i_det = 0; i_det < RecoEvent.ListHits.size(); ++i_det)
    att._logger->debug("RecoListHits det {} size {}", i_det, RecoEvent.ListHits[i_det].size());

#endif

  return 0;
}

double TWASACalibrationSimuBuilder::CloseDist(const TVector3& Xin, const TVector3& Xout, const TVector3& Pin, const TVector3& Pout)
{
  double ui = Pin.x() / Pin.z(), vi = Pin.y() / Pin.z();
  double uo = Pout.x() / Pout.z(), vo = Pout.y() / Pout.z();
  double xi = Xin.x(), yi = Xin.y();
  double xo = Xout.x() + Pout.x() / Pout.z() * (Xin.z() - Xout.z());
  double yo = Xout.y() + Pout.y() / Pout.z() * (Xin.z() - Xout.z());

  double z  = ((xi - xo) * (uo - ui) + (yi - yo) * (vo - vi)) / ((uo - ui) * (uo - ui) + (vo - vi) * (vo - vi));
  double x1 = xi + ui * z, y1 = yi + vi * z;
  double x2 = xo + uo * z, y2 = yo + vo * z;

  return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

double TWASACalibrationSimuBuilder::closestDistanceApproach(const TVector3& point1, const TVector3& point2, const TVector3& dir1, const TVector3& dir2, TVector3 &closestPoint1, TVector3 &closestPoint2)
{
  TVector3 v0 = point1 - point2;
  double d1 = dir1.Dot(dir1);
  double d12 = dir1.Dot(dir2);
  double d2 = dir2.Dot(dir2);
  double d1_v0= dir1.Dot(v0);
  double d2_v0= dir2.Dot(v0);
  double denom = d1 * d2 - d12 * d12;
  double distance = 0;
  if (TMath::Abs(denom)> 1e-5)
    {
      double sc = (d12 * d2_v0 - d2 * d1_v0) / denom;
      double tc = (d1  * d2_v0 - d12 * d1_v0) / denom;
      closestPoint1 = point1;
      closestPoint1 += dir1 * sc;
      closestPoint2 = point2;
      closestPoint2 += dir2 * tc;
      distance = (closestPoint1 - closestPoint2).Mag();
    }
  else
    {
      closestPoint1 = point1;
      closestPoint2 = point2;
      closestPoint2 += (d1_v0 / d2) * dir2 ;
      distance = (point1 - closestPoint2).Mag();
    }

  return distance;
}
