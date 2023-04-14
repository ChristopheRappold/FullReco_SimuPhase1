#include "TBuildDetectorLayerPlaneDAF.h"

#include "Debug.hh"
#include "TGeoManager.h"

#include <list>
#include <map>
#include <set>
#include <vector>

//#define DEBUG_BUILD
//#define DEBUG_BUILD2

using namespace std;

TBuildDetectorLayerPlaneDAF::TBuildDetectorLayerPlaneDAF(const THyphiAttributes& attribut) : TDataBuilder("build_det"), att(attribut)
{
  att._logger->info("TBuildDetectorLayerPlaneDAF::TBuildDetectorLayerPlaneDAF");

  par = std::make_unique<ParaManager>(att.map_ParamFiles);

  std::vector<std::string> tempName = {"HypHI_InSi_log0", "HypHI_InSi_log1", "HypHI_InSi_log2", "HypHI_InSi_log3",
    "TR1_log","TR2_log","Si1_Strip_log_x", "Si1_Strip_log_y", "Si2_Strip_log_x", "Si2_Strip_log_y",
    "SD1_Strip_log_u", "SD1_Strip_log_v", "SD2_Strip_log_u", "SD2_Strip_log_v",
    "SD1pad_Strip_log_u", "SD1pad_Strip_log_v", "SD2pad_Strip_log_u", "SD2pad_Strip_log_v",
    "TO_Counter",
    "FiberD1_Core_log_x", "FiberD1_Core_log_u", "FiberD1_Core_log_v",
    "FiberD2_Core_log_x", "FiberD2_Core_log_u", "FiberD2_Core_log_v",
    "FiberD3_Core_log_x", "FiberD3_Core_log_u", "FiberD3_Core_log_v",
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
  if(volMDCfirst == "MD01") offsetGeoNameID_MDC = 0;
  if(volMDCfirst == "SOL")  offsetGeoNameID_MDC = 1;

  auto listNodes = gGeoManager->GetVolume("INNER")->GetNodes();
  int index_lastMDC = -1, index_firstPSCE = -1;
  for(int i=0;i<listNodes->GetEntries();++i)
    {
      std::string tempName(listNodes->At(i)->GetName());
      if(tempName == "MD17_1") index_lastMDC = i;
      if(tempName == "PSCE_1") index_firstPSCE = i;
    }

  offsetGeoNameID_PSCE = index_firstPSCE - index_lastMDC + offsetGeoNameID_MDC -1;
}

TBuildDetectorLayerPlaneDAF::~TBuildDetectorLayerPlaneDAF() {}
#ifdef ROOT6
ReturnRes::InfoM TBuildDetectorLayerPlaneDAF::operator()(const TG4Sol_Event& event,
                                                         const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits,
                                                         FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
{
  int result = Exec(event, hits, RecoEvent, OutTree);

  return SoftExit(result);
}
#else
ReturnRes::InfoM TBuildDetectorLayerPlaneDAF::operator()(const TG4Sol_Event& event,
                                                         const std::vector<TClonesArray*>& hits,
                                                         FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
{
  int result = Exec(event, hits, RecoEvent, OutTree);

  return SoftExit(result);
}

ReturnRes::InfoM TBuildDetectorLayerPlaneDAF::operator()(const TG4Sol_Event& event,
                                                         const std::vector<TClonesArray*>& hits,
                                                         FullRecoEvent& RecoEvent, Ana_WasaEvent* OutTree)
{
  return ReturnRes::BuildError;
}

#endif
void TBuildDetectorLayerPlaneDAF::SelectHists()
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

ReturnRes::InfoM TBuildDetectorLayerPlaneDAF::SoftExit(int return_build)
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
int TBuildDetectorLayerPlaneDAF::Exec(const TG4Sol_Event& event, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits,
                                      FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
#else
int TBuildDetectorLayerPlaneDAF::Exec(const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits,
                                      FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
#endif
{
  int NumFilled = 0;
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
      tempInit.charge = event.BeamCharges[index];
      tempInit.posX = event.InteractionPoint_X;
      tempInit.posY = event.InteractionPoint_Y;
      tempInit.posZ = event.InteractionPoint_Z;
      tempInit.momX = event.BeamMomentums_X[index];
      tempInit.momY = event.BeamMomentums_Y[index];
      tempInit.momZ = event.BeamMomentums_Z[index];
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
      tempInit.charge = event.DaughterCharges[index];
      tempInit.posX = event.DecayVertex_X;
      tempInit.posY = event.DecayVertex_Y;
      tempInit.posZ = event.DecayVertex_Z;
      tempInit.momX = event.DaughterMomentums_X[index];
      tempInit.momY = event.DaughterMomentums_Y[index];
      tempInit.momZ = event.DaughterMomentums_Z[index];
      RecoEvent.TrackDAFInitSim.insert(std::make_pair(TrackID, tempInit));

      RecoEvent.DaughtersTrackDAFInit.insert(std::make_pair(TrackID, tempInit));

      std::vector<std::vector<SimHit> > tempSetSimHit(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackDAFSim.insert(std::make_pair(TrackID, tempSetSimHit));

      //std::vector<InfoPar> tempSetInfo(G4Sol::SIZEOF_G4SOLDETTYPE);
      //RecoEvent.TrackInfo.insert(std::make_pair(TrackID, tempSetInfo));

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
      double time_res          = 0.150; // ns
      if(nameTempBr == "FMF2_log" || nameTempBr == "HypHI_TrackFwd_log")
        {
          // continue;
#ifdef ROOT6
          for(auto it_hit = tempHits->begin(), it_hit_end = tempHits->end(); it_hit != it_hit_end; ++it_hit)
#else
          for(size_t it_hit = 0; it_hit < tempHits->GetEntries(); ++it_hit)
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

              TVectorD hitCoords(2);
              hitCoords(0) = gRandom->Gaus(hit.HitPosX, resolution_planar);
              hitCoords(1) = gRandom->Gaus(hit.HitPosY, resolution_planar);
              // hitCoords(2) = gRandom->Gaus(hit.HitPosZ, resolution);

              TVectorD hitCoordsTree(3);
              hitCoordsTree(0) = hitCoords(0);
              hitCoordsTree(1) = hitCoords(1);
              hitCoordsTree(2) = hit.HitPosZ;

              TMatrixDSym hitCov(2);
              hitCov(0, 0) = resolution_planar * resolution_planar;
              hitCov(1, 1) = resolution_planar * resolution_planar;
              // hitCov(2, 2) = resolution * resolution;
              // measurement = std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet) + LayerID, 0,
              // nullptr);

              // TVector3 o(0.,0.,hit.HitPosZ), u(1.,0.,0.), v(0.,1.,0.);
              // genfit::SharedPlanePtr plane(new genfit::DetPlane(o,u,v));
              // dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

              // RecoEvent.ListHits[TypeDet + LayerID].emplace_back(measurement.release());
              // int indexHit = RecoEvent.ListHits[TypeDet + LayerID].size() - 1;

              // tempTrack->second[TypeDet + LayerID] = indexHit;

              int pdg_code = pid_fromName(hit.Pname);

              auto tempTrackSimLayers = RecoEvent.TrackDAFSim.find(TrackID);
              SimHit tempHitSim;
              tempHitSim.layerID = LayerID;
              tempHitSim.hitX    = hit.HitPosX;
              tempHitSim.hitY    = hit.HitPosY;
              tempHitSim.hitZ    = hit.HitPosZ;
              tempHitSim.momX    = hit.MomX;
              tempHitSim.momY    = hit.MomY;
              tempHitSim.momZ    = hit.MomZ;
              tempHitSim.pdg     = pdg_code;
              tempHitSim.mass    = hit.Mass;
              tempHitSim.Eloss   = hit.Energy;
              tempHitSim.time    = gRandom->Gaus(hit.Time, time_res);
              tempHitSim.length  = hit.TrackLength;
              tempTrackSimLayers->second[TypeDet + LayerID].emplace_back(tempHitSim);
/*
              auto tempTrackInfo                              = RecoEvent.TrackInfo.find(TrackID);
              tempTrackInfo->second[TypeDet + LayerID].pdg    = pdg_code;
              tempTrackInfo->second[TypeDet + LayerID].momX   = hit.MomX;
              tempTrackInfo->second[TypeDet + LayerID].momY   = hit.MomY;
              tempTrackInfo->second[TypeDet + LayerID].momZ   = hit.MomZ;
              tempTrackInfo->second[TypeDet + LayerID].mass   = hit.Mass;
              tempTrackInfo->second[TypeDet + LayerID].Eloss  = hit.Energy;
              tempTrackInfo->second[TypeDet + LayerID].time   = hit.Time;
              tempTrackInfo->second[TypeDet + LayerID].length = hit.TrackLength;
*/
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
          for(size_t it_hit = 0; it_hit < tempHits->GetEntries(); ++it_hit)
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

              if(Fiber_removefragment_flag && hit.Pname == "He3")
                continue;

              // std::cout << "\nTypeDet : " << TypeDet << std::endl;
              // if(IsPlanar(TypeDet))
              //{
              //  std::cout << "planar"  << std::endl;
              //  TVectorD hitCoords(2);
              //  hitCoords(0) = gRandom->Gaus(hit.HitPosX, resolution_planar);
              //  hitCoords(1) = gRandom->Gaus(hit.HitPosY, resolution_planar);
              //  //hitCoords(2) = gRandom->Gaus(hit.HitPosZ, resolution);
              //  TMatrixDSym hitCov(2);
              //  hitCov(0, 0) = resolution_planar * resolution_planar;
              //  hitCov(1, 1) = resolution_planar * resolution_planar;
              //  //hitCov(2, 2) = resolution * resolution;
              //  measurement = std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID,
              //  nullptr);

              //  TVector3 o(0.,0.,hit.HitPosZ), u(1.,0.,0.), v(0.,1.,0.);
              //  genfit::SharedPlanePtr plane(new genfit::DetPlane(o,u,v));
              //  dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

              //  hitCoordsTree(0) = hitCoords(0);
              //  hitCoordsTree(1) = hitCoords(1);
              //  hitCoordsTree(2) = hit.HitPosZ;
              //}
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
#ifdef DEBUG_BUILD2
                  std::cout << "PSC" << std::endl;
                  std::string tmpName = orderDetName.find(TypeDet)->second;
                  std::cout << "name : " << tmpName << std::endl;
                  std::cout << "LayerID : " << LayerID-1 <<" : " <<TypeDet - G4Sol::MG01 + offsetGeoNameID_PSCE + LayerID - 1<<std::endl;
                  std::cout << "HitPosX : " << hit.HitPosX << std::endl;
                  std::cout << "HitPosY : " << hit.HitPosY << std::endl;
                  std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
                  gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_PSCE + LayerID - 1)->Print();
                  gGeoManager->GetVolume("INNER")
                      ->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_PSCE + LayerID - 1)
                      ->GetMatrix()
                      ->Print();
                  gGeoManager->GetVolume("MFLD")->GetNode(0)->Print();
                  gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix()->Print();
                  gGeoManager->GetVolume("WASA")->GetNode(0)->Print();
                  gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix()->Print();
#endif
                  TGeoMatrix* g1 = gGeoManager->GetVolume("INNER")
                                       ->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_PSCE + LayerID - 1)
                                       ->GetMatrix();                                       // PSCE
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
                  double* shift = H.GetTranslation();
                  TVector3 o(shift[0], shift[1], shift[2]), phidir(shift[0], shift[1], 0), zdir(0., 0., 1.);
                  phidir     = phidir.Unit();
                  TVector3 u = zdir.Cross(phidir);
                  TVector3 v = zdir;
                  genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));

                  TVectorD hitCoords(2);
                  hitCoords(0) = 0.;
                  hitCoords(1) = gRandom->Gaus(hit.HitPosZ - shift[2], resolution_psce_z);
                  TMatrixDSym hitCov(2);
                  hitCov(0, 0) = resolution_psce * resolution_psce;
                  hitCov(1, 1) = resolution_psce_z * resolution_psce_z;
                  measurement =
                      std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID-1, nullptr);
                  dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

                  measinfo.SetInfo(-9999., gRandom->Gaus(hit.Time, time_res), hit.Energy);

                  hitCoordsTree(0) = hit.HitPosX;
                  hitCoordsTree(1) = hit.HitPosY;
                  hitCoordsTree(2) = hit.HitPosZ;

                  LocalHisto.h76->Fill(LayerID-1, atan2(shift[1], shift[0]));
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
                  TGeoMatrix* g2   = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix();          // INNNER
                  TGeoMatrix* g3   = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();          // MFLD
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

                  TVectorD hitCoords(2);
                  hitCoords(0) = gRandom->Uniform(-3.75, 3.75); // phi ! be aware ! not u-dim
                  hitCoords(1) = gRandom->Uniform(6., 22.);     // r -> v dir

                  TMatrixDSym hitCov(2);
                  hitCov(0, 0) = TMath::Sq(2 * hitCoords(1) * TMath::Sin(3.75 * TMath::DegToRad())) / 12.;
                  hitCov(1, 1) = TMath::Sq(22. - 6.) / 12.;
                  measurement =
                      std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID, nullptr);
                  dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

                  measinfo.SetInfo(-9999., gRandom->Gaus(hit.Time, time_res), hit.Energy);

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
                  TGeoMatrix* g2   = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix();          // INNNER
                  TGeoMatrix* g3   = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();          // MFLD
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

                  TVectorD hitCoords(2);
                  hitCoords(0) = gRandom->Uniform(-3.75, 3.75); // phi ! be aware ! not u-dim
                  hitCoords(1) = gRandom->Uniform(6., 22.);     // r -> v dir

                  TMatrixDSym hitCov(2);
                  hitCov(0, 0) = TMath::Sq(2 * hitCoords(1) * TMath::Sin(3.75 * TMath::DegToRad())) / 12.;
                  hitCov(1, 1) = TMath::Sq(22. - 6.) / 12.;
                  measurement =
                      std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID, nullptr);
                  dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

                  measinfo.SetInfo(-9999., gRandom->Gaus(hit.Time, time_res), hit.Energy);

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

                  FiberHitAna* hit_ana = new FiberHitAna(i_fiber, i_layer, LayerID/2,
                                  gRandom->Gaus(hit.Time, time_res), hit.Energy, TrackID, par.get());
                  FiberHitCont[hit_ana->GetDet()][hit_ana->GetLay()].emplace_back(hit_ana);

                  int pdg_code = pid_fromName(hit.Pname);
                  if(pdg_code == 0)
                    att._logger->debug("!> Builder : pdg_code = 0 ! {}", hit.Pname);

                  auto tempTrackSimFibers = RecoEvent.TrackDAFSim.find(TrackID);
                  SimHit tempHitSim;
                  tempHitSim.layerID = LayerID;
                  tempHitSim.hitX    = hit.HitPosX;
                  tempHitSim.hitY    = hit.HitPosY;
                  tempHitSim.hitZ    = hit.HitPosZ;
                  tempHitSim.momX    = hit.MomX;
                  tempHitSim.momY    = hit.MomY;
                  tempHitSim.momZ    = hit.MomZ;
                  tempHitSim.pdg     = pdg_code;
                  tempHitSim.mass    = hit.Mass;
                  tempHitSim.Eloss   = hit.Energy;
                  tempHitSim.time    = hit.Time;
                  tempHitSim.length  = hit.TrackLength;
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

/*
                  //if(IsFiberU_Vetoed(TypeDet))
                  //  continue;
                  int leftright = LayerID%2;
                  if(used_fiber_pair.size()>0 && used_fiber_pair.find(LayerID+std::pow(-1,leftright))!=used_fiber_pair.end())
                    continue;

                  string volumeName;
                  string motherName;
                  int i_fiber;
                  int i_lay;
                  switch(TypeDet)
                    {
                    case G4Sol::FiberD1_x: {
                      volumeName = "FiberD1_log_x";
                      motherName = "FiberD1_log_0";
                      i_fiber = 0;
                      i_lay   = 0; } ;
                      break;
                    case G4Sol::FiberD1_u: {
                      volumeName = "FiberD1_log_u";
                      motherName = "FiberD1_log_0";
                      i_fiber = 0;
                      i_lay   = 1; } ;
                      break;
                    case G4Sol::FiberD1_v: {
                      volumeName = "FiberD1_log_v";
                      motherName = "FiberD1_log_0";
                      i_fiber = 0;
                      i_lay   = 2; } ;
                      break;
                    case G4Sol::FiberD2_x: {
                      volumeName = "FiberD2_log_x";
                      motherName = "FiberD2_log_0";
                      i_fiber = 1;
                      i_lay   = 0; } ;
                      break;
                    case G4Sol::FiberD2_u: {
                      volumeName = "FiberD2_log_u";
                      motherName = "FiberD2_log_0";
                      i_fiber = 1;
                      i_lay   = 1; } ;
                      break;
                    case G4Sol::FiberD2_v: {
                      volumeName = "FiberD2_log_v";
                      motherName = "FiberD2_log_0";
                      i_fiber = 1;
                      i_lay   = 2; } ;
                      break;
                    case G4Sol::FiberD3_x: {
                      volumeName = "FiberD3_log_x";
                      motherName = "FiberD3_log_0";
                      i_fiber = 2;
                      i_lay   = 0; } ;
                      break;
                    case G4Sol::FiberD3_u: {
                      volumeName = "FiberD3_log_u";
                      motherName = "FiberD3_log_0";
                      i_fiber = 2;
                      i_lay   = 1; } ;
                      break;
                    case G4Sol::FiberD3_v: {
                      volumeName = "FiberD3_log_v";
                      motherName = "FiberD3_log_0";
                      i_fiber = 2;
                      i_lay   = 2; } ;
                      break;
                    case G4Sol::FiberD4_x: {
                      volumeName = "FiberD4_log_x";
                      motherName = "FiberD4_log_0";
                      i_fiber = 5;
                      i_lay   = 0; } ;
                      break;
                    case G4Sol::FiberD4_u: {
                      volumeName = "FiberD4_log_u";
                      motherName = "FiberD4_log_0";
                      i_fiber = 5;
                      i_lay   = 1; } ;
                      break;
                    case G4Sol::FiberD4_v: {
                      volumeName = "FiberD4_log_v";
                      motherName = "FiberD4_log_0";
                      i_fiber = 5;
                      i_lay   = 2; } ;
                      break;
                    case G4Sol::FiberD5_x: {
                      volumeName = "FiberD5_log_x";
                      motherName = "FiberD5_log_0";
                      i_fiber = 6;
                      i_lay   = 0; } ;
                      break;
                    case G4Sol::FiberD5_u: {
                      volumeName = "FiberD5_log_u";
                      motherName = "FiberD5_log_0";
                      i_fiber = 6;
                      i_lay   = 1; } ;
                      break;
                    case G4Sol::FiberD5_v: {
                      volumeName = "FiberD5_log_v";
                      motherName = "FiberD5_log_0";
                      i_fiber = 6;
                      i_lay   = 2; } ;
                      break;
                    default:
                      std::cerr << "something wrong" << std::endl;
                      break;
                    }

#ifdef DEBUG_BUILD2
                  std::cout << "fiber" << std::endl;
                  std::string tmpName = orderDetName.find(TypeDet)->second;
                  std::cout << "name : " << tmpName << std::endl;
                  std::cout << "LayerID : " << LayerID << std::endl;
                  std::cout << "HitPosX : " << hit.HitPosX << std::endl;
                  std::cout << "HitPosY : " << hit.HitPosY << std::endl;
                  std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
                  gGeoManager->GetVolume(volumeName.c_str())->GetNode(LayerID * 2 + 1)->Print();
                  gGeoManager->GetVolume(volumeName.c_str())->GetNode(LayerID * 2 + 1)->GetMatrix()->Print();
                  gGeoManager->GetVolume("MFLD")
                      ->GetNode(motherName.c_str())
                      ->GetVolume()
                      ->GetNode((volumeName + "_0").c_str())
                      ->Print();
                  gGeoManager->GetVolume("MFLD")
                      ->GetNode(motherName.c_str())
                      ->GetVolume()
                      ->GetNode((volumeName + "_0").c_str())
                      ->GetMatrix()
                      ->Print();
                  gGeoManager->GetVolume("MFLD")->GetNode(motherName.c_str())->Print();
                  gGeoManager->GetVolume("MFLD")->GetNode(motherName.c_str())->GetMatrix()->Print();
#endif
                  TGeoMatrix* g1 =
                      gGeoManager->GetVolume(volumeName.c_str())->GetNode(LayerID * 2 + 1)->GetMatrix(); // fiber hit core
                  TGeoMatrix* g1_pair =
                      gGeoManager->GetVolume(volumeName.c_str())->GetNode((LayerID + std::pow(-1, leftright)) * 2 + 1)->GetMatrix(); // fiber pair core
                  TGeoHMatrix H1(*g1), H1_pair(*g1_pair);
                  Double_t* center_hit = H1.GetTranslation();
                  Double_t* center_pair = H1_pair.GetTranslation();
                  Double_t center_both[3];
                  center_both[0] = (center_hit[0] + center_pair[0]) / 2.;
                  center_both[1] = (center_hit[1] + center_pair[1]) / 2.;
                  center_both[2] = (center_hit[2] + center_pair[2]) / 2.;

                  //if(Fiber_moveXlayer_flag && i_lay == 0)
                  //  center_both[0] += Fiber_moveXlayer_stepsize * Fiber_moveXlayer_ntimes;

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
                  double* shift = H.GetTranslation();
                  TVector3 o(0., 0., shift[2]), zdir(0., 0., 1.);
                  TVector3 fiber_dir(edge2[0] - edge1[0], edge2[1] - edge1[1], edge2[2] - edge1[2]);
                  fiber_dir  = fiber_dir.Unit();
                  TVector3 u = fiber_dir.Cross(zdir);
                  TVector3 v = fiber_dir;
                  genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));

                  TVectorD hitCoords(1);
                  hitCoords(0) = u.Dot(TVector3(shift[0], shift[1], 0));
                  TMatrixDSym hitCov(1);
                  hitCov(0, 0) = resolution_fiber * resolution_fiber;
                  measurement =
                      std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID, nullptr);
                  dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

                  measinfo.SetInfo(-9999., gRandom->Gaus(hit.Time, time_res), hit.Energy);

                  hitCoordsTree(0) = hit.HitPosX;
                  hitCoordsTree(1) = hit.HitPosY;
                  hitCoordsTree(2) = hit.HitPosZ;

                  //if(Fiber_moveXlayer_flag && i_lay == 0)
                  //  hitCoordsTree(0) += Fiber_moveXlayer_stepsize * Fiber_moveXlayer_ntimes;

                  used_fiber_pair.insert(LayerID);
                }
              else if(IsFiberM(TypeDet))
                {
                  int leftright = LayerID%2;
                  if(used_fiber_pair.size()>0 && used_fiber_pair.find(LayerID+std::pow(-1,leftright))!=used_fiber_pair.end())
                    continue;

                  string volumeName;
                  string motherName;
                  int i_fiber;
                  int i_lay;
                  switch(TypeDet)
                    {
                    case G4Sol::MiniFiberD1_x: {
                      volumeName = "MiniFiberD1_log_x";
                      motherName = "MiniFiberD1_log_0";
                      i_fiber = 3;
                      i_lay   = 0; } ;
                      break;
                    case G4Sol::MiniFiberD1_u: {
                      volumeName = "MiniFiberD1_log_u";
                      motherName = "MiniFiberD1_log_0";
                      i_fiber = 3;
                      i_lay   = 1; } ;
                      break;
                    case G4Sol::MiniFiberD1_v: {
                      volumeName = "MiniFiberD1_log_v";
                      motherName = "MiniFiberD1_log_0";
                      i_fiber = 3;
                      i_lay   = 2; } ;
                      break;
                    case G4Sol::MiniFiberD2_x: {
                      volumeName = "MiniFiberD2_log_x";
                      motherName = "MiniFiberD2_log_0";
                      i_fiber = 4;
                      i_lay   = 0; } ;
                      break;
                    case G4Sol::MiniFiberD2_v: {
                      volumeName = "MiniFiberD2_log_v";
                      motherName = "MiniFiberD2_log_0";
                      i_fiber = 4;
                      i_lay   = 1; } ;
                      break;
                    case G4Sol::MiniFiberD2_u: {
                      volumeName = "MiniFiberD2_log_u";
                      motherName = "MiniFiberD2_log_0";
                      i_fiber = 4;
                      i_lay   = 2; } ;
                      break;
                    default:
                      std::cerr << "something wrong" << std::endl;
                      break;
                    }

#ifdef DEBUG_BUILD2
                  std::cout << "fiberM" << std::endl;
                  std::string tmpName = orderDetName.find(TypeDet)->second;
                  std::cout << "name : " << tmpName << std::endl;
                  std::cout << "LayerID : " << LayerID << std::endl;
                  std::cout << "HitPosX : " << hit.HitPosX << std::endl;
                  std::cout << "HitPosY : " << hit.HitPosY << std::endl;
                  std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
                  // gGeoManager->GetVolume(volumeName.c_str())->Print();
                  gGeoManager->GetVolume(volumeName.c_str())->GetNode(LayerID * 2 + 1)->Print();
                  gGeoManager->GetVolume(volumeName.c_str())->GetNode(LayerID * 2 + 1)->GetMatrix()->Print();
                  gGeoManager->GetVolume("MFLD")
                      ->GetNode(motherName.c_str())
                      ->GetVolume()
                      ->GetNode((volumeName + "_0").c_str())
                      ->Print();
                  gGeoManager->GetVolume("MFLD")
                      ->GetNode(motherName.c_str())
                      ->GetVolume()
                      ->GetNode((volumeName + "_0").c_str())
                      ->GetMatrix()
                      ->Print();
                  gGeoManager->GetVolume("MFLD")->GetNode(motherName.c_str())->Print();
                  gGeoManager->GetVolume("MFLD")->GetNode(motherName.c_str())->GetMatrix()->Print();
#endif
                  TGeoMatrix* g1 =
                      gGeoManager->GetVolume(volumeName.c_str())->GetNode(LayerID * 2 + 1)->GetMatrix(); // minifiber hit core
                  TGeoMatrix* g1_pair =
                      gGeoManager->GetVolume(volumeName.c_str())->GetNode((LayerID + std::pow(-1, leftright)) * 2 + 1)->GetMatrix(); // minifiber pair core
                  TGeoHMatrix H1(*g1), H1_pair(*g1_pair);
                  Double_t* center_hit = H1.GetTranslation();
                  Double_t* center_pair = H1_pair.GetTranslation();
                  Double_t center_both[3];
                  center_both[0] = (center_hit[0] + center_pair[0]) / 2.;
                  center_both[1] = (center_hit[1] + center_pair[1]) / 2.;
                  center_both[2] = (center_hit[2] + center_pair[2]) / 2.;

                  //if(Fiber_moveXlayer_flag && i_lay == 0)
                  //  center_both[0] += Fiber_moveXlayer_stepsize * Fiber_moveXlayer_ntimes;

                  H1.SetTranslation(center_both);

                  TGeoMatrix* g2 = gGeoManager->GetVolume("MFLD")
                                       ->GetNode(motherName.c_str())
                                       ->GetVolume()
                                       ->GetNode((volumeName + "_0").c_str())
                                       ->GetMatrix(); // minifiber layer
                  TGeoMatrix* g3 =
                      gGeoManager->GetVolume("MFLD")->GetNode(motherName.c_str())->GetMatrix(); // minifiber station
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
                  TVector3 o(0., 0., shift[2]), zdir(0., 0., 1.);
                  TVector3 fiber_dir(edge2[0] - edge1[0], edge2[1] - edge1[1], edge2[2] - edge1[2]);
                  fiber_dir  = fiber_dir.Unit();
                  TVector3 u = fiber_dir.Cross(zdir);
                  TVector3 v = fiber_dir;
                  genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));

                  TVectorD hitCoords(1);
                  hitCoords(0) = u.Dot(TVector3(shift[0], shift[1], 0));
                  TMatrixDSym hitCov(1);
                  hitCov(0, 0) = resolution_fiber * resolution_fiber;
                  measurement =
                      std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID, nullptr);
                  dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

                  measinfo.SetInfo(-9999., gRandom->Gaus(hit.Time, time_res), hit.Energy);

                  hitCoordsTree(0) = hit.HitPosX;
                  hitCoordsTree(1) = hit.HitPosY;
                  hitCoordsTree(2) = hit.HitPosZ;

                  //if(Fiber_moveXlayer_flag && i_lay == 0)
                  //  hitCoordsTree(0) += Fiber_moveXlayer_stepsize * Fiber_moveXlayer_ntimes;
                  
                  used_fiber_pair.insert(LayerID);
                }
*/
              else if(IsWire(TypeDet))
                {
#ifdef DEBUG_BUILD2
                  std::cout << "wire" << std::endl;
                  std::string tmpName = orderDetName.find(TypeDet)->second;
                  std::cout << "name : " << tmpName << std::endl;
                  std::cout << "LayerID : " << LayerID <<" : "<< TypeDet - G4Sol::MG01 + 1 <<std::endl;
                  std::cout << "HitPosX : " << hit.HitPosX << std::endl;
                  std::cout << "HitPosY : " << hit.HitPosY << std::endl;
                  std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
                  gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_MDC)->Print();
                  gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_MDC)->GetMatrix()->Print();
                  gGeoManager->GetVolume("INNER")
                                ->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_MDC)
                                ->GetVolume()
                                ->GetNode(LayerID - 1)
                                ->Print();
                  gGeoManager->GetVolume("INNER")
                                ->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_MDC)
                                ->GetVolume()
                                ->GetNode(LayerID - 1)
                                ->GetMatrix()
                                ->Print();
                  gGeoManager->GetVolume("MFLD")->GetNode(0)->Print();
                  gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix()->Print();
                  gGeoManager->GetVolume("WASA")->GetNode(0)->Print();
                  gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix()->Print();
#endif
                  TGeoMatrix* g1 = gGeoManager->GetVolume("INNER")
		                              ->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_MDC)
		                              ->GetVolume()
                                       ->GetNode(LayerID - 1)
                                       ->GetMatrix(); // ME, MG
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
                  double* shift = H.GetTranslation();
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
                  double* edge1 = Hw1.GetTranslation();
                  double* edge2 = Hw2.GetTranslation();

                  TVector3 x1(shift[0], shift[1], shift[2]);
                  TVector3 p1(edge2[0] - edge1[0], edge2[1] - edge1[1], edge2[2] - edge1[2]);
                  TVector3 x2(hit.HitPosX, hit.HitPosY, hit.HitPosZ);
                  TVector3 p2(hit.MomX, hit.MomY, hit.MomZ);
                  double dl = CloseDist(x1, x2, p1, p2);

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

                  measinfo.SetInfo(-9999., gRandom->Gaus(hit.Time, time_res), hit.Energy); //Change check TOT?

                  hitCoordsTree(0) = hit.HitPosX;
                  hitCoordsTree(1) = hit.HitPosY;
                  hitCoordsTree(2) = hit.HitPosZ;
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
                  hitCov(0, 0) = resolution_wire * resolution_wire;
                  hitCov(1, 1) = resolution_wire * resolution_wire;
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

              //int indexHit = RecoEvent.ListHits[TypeDet].size() - 1;
              //tempTrack->second[TypeDet] = indexHit;

              int pdg_code = pid_fromName(hit.Pname);
              if(pdg_code == 0)
                att._logger->debug("!> Builder : pdg_code = 0 ! {}", hit.Pname);

              //printf("Fibers id[%d, %d] ; TypeDet: %d \n", G4Sol::FiberD1_x, G4Sol::FiberD5_v, TypeDet);
              //printf("TrackID: %d ; hit.trackID: %d \n", TrackID, hit.TrackID);

              //for(auto iter = RecoEvent.TrackDAFSim.begin(); iter != RecoEvent.TrackDAFSim.end(); ++iter)
              //  printf("Map TrackIDs: %d \n", iter->first);


              auto tempTrackSim = RecoEvent.TrackDAFSim.find(TrackID);
              SimHit tempHitSim;
              tempHitSim.layerID = LayerID;
              tempHitSim.hitX    = hit.HitPosX;
              tempHitSim.hitY    = hit.HitPosY;
              tempHitSim.hitZ    = hit.HitPosZ;
              tempHitSim.momX    = hit.MomX;
              tempHitSim.momY    = hit.MomY;
              tempHitSim.momZ    = hit.MomZ;
              tempHitSim.pdg     = pdg_code;
              tempHitSim.mass    = hit.Mass;
              tempHitSim.Eloss   = hit.Energy;
              tempHitSim.time    = hit.Time;
              tempHitSim.length  = hit.TrackLength;
              tempTrackSim->second[TypeDet].emplace_back(tempHitSim);
/*
              auto tempTrackInfo                    = RecoEvent.TrackInfo.find(TrackID);
              tempTrackInfo->second[TypeDet].pdg    = pdg_code;
              tempTrackInfo->second[TypeDet].momX   = hit.MomX;
              tempTrackInfo->second[TypeDet].momY   = hit.MomY;
              tempTrackInfo->second[TypeDet].momZ   = hit.MomZ;
              tempTrackInfo->second[TypeDet].mass   = hit.Mass;
              tempTrackInfo->second[TypeDet].Eloss  = hit.Energy;
              tempTrackInfo->second[TypeDet].time   = hit.Time;
              tempTrackInfo->second[TypeDet].length = hit.TrackLength;
*/
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
/*
              if(TypeDet >= G4Sol::FiberD1_x && TypeDet <= G4Sol::FiberD3_v)
                fillOutHit(OutTree->Fiber, hit, pdg_code, charge, hitCoordsTree, TypeDet, LayerID);
              if(TypeDet >= G4Sol::MiniFiberD1_x && TypeDet <= G4Sol::MiniFiberD2_u)
                fillOutHit(OutTree->Fiber, hit, pdg_code, charge, hitCoordsTree, TypeDet, LayerID);
              if(TypeDet >= G4Sol::FiberD4_v && TypeDet <= G4Sol::FiberD5_v)
                fillOutHit(OutTree->Fiber, hit, pdg_code, charge, hitCoordsTree, TypeDet, LayerID);
*/
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


  // Fiber ana //////////////////////////////////
  FiberAnalyzer* fiberana = new FiberAnalyzer();
  FiberHitClCont          = fiberana->Clusterize(FiberHitCont);

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

              MeasurementInfo measinfo(FiberHitClCont[i][j][k]->GetTOT(), FiberHitClCont[i][j][k]->GetTime(), FiberHitClCont[i][j][k]->GetdE());
              RecoEvent.ListHitsInfo[TypeDet[i][j]].emplace_back(measinfo);

              RecoEvent.ListHitsToTracks[TypeDet[i][j]].emplace_back(FiberHitClCont[i][j][k]->GetSimTrackID());
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
      //std::cout << "\nsize before : " << FiberHitCont[i][j].size() << std::endl;
      //std::cout << "size after  : " << FiberHitClCont[i][j].size() << std::endl;
      for(int k=0; k<(int)FiberHitClCont[i][j].size(); ++k){
        //FiberHitClCont[i][j][k]->Print();
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
                  = fiberana->FiberTracking(FiberHitClCont, par.get(), RecoEvent.ListHits[G4Sol::PSCE]);

  RecoEvent.FiberTrackCont = FiberTrackCont;

  int num_combi_mft12 = 1;
  for(int i=3; i<5; ++i)
    for(int j=0; j<3; ++j)
      num_combi_mft12 *= ( (int)FiberHitClCont[i][j].size() + 1 );

  if(par->flag_mft12_combi && num_combi_mft12<par->cut_mft12_combi){
    LocalHisto.hfiber_1_1->Fill(num_combi_mft12);
    LocalHisto.hfiber_1_2->Fill((double)num_combi_mft12*1e-3);
    LocalHisto.hfiber_1_3->Fill((double)num_combi_mft12*1e-3);
    LocalHisto.hfiber_1_4->Fill((double)num_combi_mft12*1e-6);
  }

  LocalHisto.hfiber_1_5->Fill((double)num_combi_mft12*1e-3);
  LocalHisto.hfiber_1_6->Fill((double)num_combi_mft12*1e-3);
  LocalHisto.hfiber_1_7->Fill((double)num_combi_mft12*1e-6);
  if(par->flag_debug) std::cout << "- num_combi_mft12 : " << num_combi_mft12 << std::endl;

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

  if(par->flag_dft12_combi && num_combi_dft12<par->cut_dft12_combi){
    LocalHisto.hfiber_5_1->Fill(num_combi_dft12);
    LocalHisto.hfiber_5_2->Fill((double)num_combi_dft12*1e-3);
    LocalHisto.hfiber_5_3->Fill((double)num_combi_dft12*1e-3);
    LocalHisto.hfiber_5_4->Fill((double)num_combi_dft12*1e-6);
  }

  LocalHisto.hfiber_5_5->Fill((double)num_combi_dft12*1e-3);
  LocalHisto.hfiber_5_6->Fill((double)num_combi_dft12*1e-3);
  LocalHisto.hfiber_5_7->Fill((double)num_combi_dft12*1e-6);
  if(par->flag_debug) std::cout << "- num_combi_dft12 : " << num_combi_dft12 << std::endl;

  for(int i=0; i<FiberTrackCont["dft12"].size(); ++i){
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

/*
  std::vector< std::vector<FiberHitXUV*> > FiberXUVCont = fiberana->FindHit(FiberHitClCont, par.get());

  for(int i=0; i<7; ++i)
    {
      LocalHisto.h17[i]->Fill( FiberXUVCont[i].size() );
      for(int j=0; j<(int)FiberXUVCont[i].size(); ++j)
        {
          LocalHisto.h16[i]->Fill(FiberXUVCont[i][j]->GetPosX(), FiberXUVCont[i][j]->GetPosY());
          LocalHisto.h17_2[i]->Fill(FiberXUVCont[i][j]->GetD());
        }
    }


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

    if(par->flag_dup_mft12_xuv && (int)buf_track.size()>0) buf_track = fiberana->DeleteDup(buf_track);
    FiberTrackCont["mft12"] = buf_track;

    nt_mft12     = FiberTrackCont["mft12"].size();
    nt_mft12_xuv = FiberTrackCont["mft12"].size();
    for(auto v: FiberTrackCont["mft12"]){
      for(int i=0; i<6; ++i){
        if(par->flag_dup_mft12_xuv) v->GetContHit().at(i)->SetUsed();
      }
      v->CorrectMFT(par.get());
      v->SetPosL();
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
  if(par->flag_debug) std::cout << "- num_combi_mft12 : " << num_combi_mft12 << std::endl;

  if(par->flag_mft12_combi && num_combi_mft12<par->cut_mft12_combi){
    LocalHisto.hfiber_1_1->Fill(num_combi_mft12);
    LocalHisto.hfiber_1_2->Fill((double)num_combi_mft12*1e-3);
    LocalHisto.hfiber_1_3->Fill((double)num_combi_mft12*1e-3);
    LocalHisto.hfiber_1_4->Fill((double)num_combi_mft12*1e-6);
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
                FiberTrackAna *track = new FiberTrackAna(buf_hit, par.get());
                track->SetFlagCombi();
                if(par->flag_mft12_posang){
                  track->CorrectMFTCombi(par.get());
                  double buf_x = track->GetXmft();
                  double buf_y = track->GetYmft();
                  double buf_a = track->GetA();
                  double buf_b = track->GetB();
                  if( fabs(buf_x * 0.003 - buf_a)>0.3 || fabs(buf_y * 0.003 - buf_b)>0.3 ){ delete track; continue;}
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
      for(int j=0; j<(int)RecoEvent.ListHits[G4Sol::PSCE].size(); ++j){
        double a_fiber = buf_track[i]->GetA();
        double b_fiber = buf_track[i]->GetB();
        double x_fiber = buf_track[i]->GetX();
        double y_fiber = buf_track[i]->GetY();
        double phi_psb   = GetPSB_Phi(RecoEvent.ListHits[G4Sol::PSCE][j]->getHitId()) + att.psb_rot_z*Deg2Rad;
        double r_psb     = GetPSB_R(RecoEvent.ListHits[G4Sol::PSCE][j]->getHitId());
        double z_psb     = RecoEvent.ListHits[G4Sol::PSCE][j]->getRawHitCoords()[1]*10.; //in mm

        double par_a = pow(a_fiber, 2) + pow(b_fiber, 2);
        double par_b = a_fiber * (x_fiber - att.psb_pos_x)+ b_fiber * (y_fiber - att.psb_pos_y);
        double par_c = pow(x_fiber - att.psb_pos_x, 2) + pow(y_fiber - att.psb_pos_y, 2) - pow(r_psb, 2);
        double z_fiber  = (-par_b + sqrt( pow(par_b,2) - par_a * par_c)) / par_a;

        double fiber_x_buf = x_fiber + a_fiber * z_fiber - att.psb_pos_x;
        double fiber_y_buf = y_fiber + b_fiber * z_fiber - att.psb_pos_y;
        double phi_fiber = atan2(fiber_y_buf, fiber_x_buf);
        LocalHisto.hfiber_6_1->Fill( fiberana->CalcPhiDif(phi_psb, phi_fiber) );
        LocalHisto.hfiber_6_2->Fill( (z_fiber - att.psb_pos_z) - z_psb );
        LocalHisto.hfiber_6_3->Fill( phi_fiber * Rad2Deg, phi_psb * Rad2Deg );
        LocalHisto.hfiber_6_4->Fill( (z_fiber - att.psb_pos_z),  z_psb );
        if( fabs(fiberana->CalcPhiDif(phi_psb, phi_fiber)) < att.cut_psb_phi ){
          if(fabs( (z_fiber - att.psb_pos_z) - z_psb)<att.cut_psb_z){
            if(buf_track[i]->IsFlagPSB()){
                if( fabs( fiberana->CalcPhiDif(phi_psb, phi_fiber) ) < fabs( buf_track[i]->GetPSBDifPhi() ) ){
                  buf_track[i]->SetSegPSB(RecoEvent.ListHits[G4Sol::PSCE][j]->getHitId());
                  buf_track[i]->SetPSBDifZ(z_psb - (z_fiber - att.psb_pos_z));
                  buf_track[i]->SetPSBDifPhi(fiberana->CalcPhiDif(phi_psb, phi_fiber));
                }
            }
            else{
              buf_track[i]->SetFlagPSB();
              buf_track[i]->SetSegPSB(RecoEvent.ListHits[G4Sol::PSCE][j]->getHitId());
              buf_track[i]->SetPSBDifZ(z_psb - (z_fiber - att.psb_pos_z));
              buf_track[i]->SetPSBDifPhi(fiberana->CalcPhiDif(phi_psb, phi_fiber));
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
      v->SetPosL();
      v->DelFlagPSB();
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
              double buf_x = track->GetXmft();
              double buf_y = track->GetYmft();
              double buf_a = track->GetA();
              double buf_b = track->GetB();
              if( fabs(buf_x * 0.003 - buf_a)>0.3 || fabs(buf_y * 0.003 - buf_b)>0.3 ){ delete track; continue;}
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
        for(int j=0; j<(int)RecoEvent.ListHits[G4Sol::PSCE].size(); ++j){
          double a_fiber = buf_track[i]->GetA();
          double b_fiber = buf_track[i]->GetB();
          double x_fiber = buf_track[i]->GetX();
          double y_fiber = buf_track[i]->GetY();
          double phi_psb   = GetPSB_Phi(RecoEvent.ListHits[G4Sol::PSCE][j]->getHitId()) + att.psb_rot_z*Deg2Rad;
          double r_psb     = GetPSB_R(RecoEvent.ListHits[G4Sol::PSCE][j]->getHitId());
          double z_psb     = RecoEvent.ListHits[G4Sol::PSCE][j]->getRawHitCoords()[1]*10.; //in mm

          double par_a = pow(a_fiber, 2) + pow(b_fiber, 2);
          double par_b = a_fiber * (x_fiber - att.psb_pos_x)+ b_fiber * (y_fiber - att.psb_pos_y);
          double par_c = pow(x_fiber - att.psb_pos_x, 2) + pow(y_fiber - att.psb_pos_y, 2) - pow(r_psb, 2);
          double z_fiber  = (-par_b + sqrt( pow(par_b,2) - par_a * par_c)) / par_a;

          double fiber_x_buf = x_fiber + a_fiber * z_fiber - att.psb_pos_x;
          double fiber_y_buf = y_fiber + b_fiber * z_fiber - att.psb_pos_y;
          double phi_fiber = atan2(fiber_y_buf, fiber_x_buf);
          if( fabs(fiberana->CalcPhiDif(phi_psb, phi_fiber)) < att.cut_psb_phi ){
            if(fabs( (z_fiber - att.psb_pos_z) - z_psb)<att.cut_psb_z){
              if(buf_track[i]->IsFlagPSB()){
                if( fabs( fiberana->CalcPhiDif(phi_psb, phi_fiber) ) < fabs( buf_track[i]->GetPSBDifPhi() ) ){
                  buf_track[i]->SetSegPSB(RecoEvent.ListHits[G4Sol::PSCE][j]->getHitId());
                  buf_track[i]->SetPSBDifZ(z_psb - (z_fiber - att.psb_pos_z));
                  buf_track[i]->SetPSBDifPhi(fiberana->CalcPhiDif(phi_psb, phi_fiber));
                }
              }

              else{
                buf_track[i]->SetFlagPSB();
                buf_track[i]->SetSegPSB(RecoEvent.ListHits[G4Sol::PSCE][j]->getHitId());
                buf_track[i]->SetPSBDifZ(z_psb - (z_fiber - att.psb_pos_z));
                buf_track[i]->SetPSBDifPhi(fiberana->CalcPhiDif(phi_psb, phi_fiber));
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
        v->SetPosL();
        v->DelFlagPSB();
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


  // DFT12
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

  if(par->flag_debug) std::cout << "- dft12 end" << std::endl;
*/


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

double CloseDist(const TVector3& Xin, const TVector3& Xout, const TVector3& Pin, const TVector3& Pout)
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
/*
double TBuildDetectorLayerPlaneDAF::GetPSB_R(int _seg)
{
  double _r = -999.;

  if(_seg<23){
    if (0 == (_seg % 2)) { // inner PSB
      _r   = 217.;
    } else { // outer PSB
      _r   = 227.75;
    }
  }
  else{
    if (0 == ((_seg-23) % 2)) { // inner PSB
      _r   = 217.;
    } else { // outer PSB
      _r   = 227.75;
    }
  }
  return _r;
}


double TBuildDetectorLayerPlaneDAF::GetPSB_Phi(int _seg)
{
  double _phi = -999.;

  if(_seg<23){
    if (0 == (_seg % 2)) { // inner PSB
      _phi = TMath::Pi() * (9.15 + 14.7 * ((double)(_seg / 2))) / 180.0;
    } else { // outer PSB
      _phi = TMath::Pi() * (16.5 + 14.7 * ((double)((_seg - 1) / 2))) / 180.0;
    }
  }
  else{
    if (0 == ((_seg-23) % 2)) { // inner PSB
      _phi = TMath::Pi() * (189.15 + 14.7 * ((float)((_seg-23) / 2))) / 180.0;
    } else { // outer PSB
      _phi = TMath::Pi() * (196.5 + 14.7 * ((float)(((_seg-23) - 1) / 2))) / 180.0;
    }
  }

  _phi += TMath::Pi()/2.;
  if( _phi > TMath::Pi() ) _phi -= 2*TMath::Pi();

  return _phi;
}
*/