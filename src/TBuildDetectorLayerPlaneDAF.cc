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
/*
  att.FiberXUV_cutd[0] = par->fiber_uft1_cut_d;
  att.FiberXUV_cutd[1] = par->fiber_uft2_cut_d;
  att.FiberXUV_cutd[2] = par->fiber_uft3_cut_d;
  att.FiberXUV_cutd[3] = par->fiber_mft_cut_d;
  att.FiberXUV_cutd[4] = par->fiber_mft_cut_d;
  att.FiberXUV_cutd[5] = par->fiber_dft1_cut_d;
  att.FiberXUV_cutd[6] = par->fiber_dft2_cut_d;
*/
  std::vector<std::string> tempName = {"HypHI_InSi_log0", "HypHI_InSi_log1", "HypHI_InSi_log2", "HypHI_InSi_log3",
    "TR1_log","TR2_log","Si1_Strip_log_x", "Si1_Strip_log_y", "Si2_Strip_log_x", "Si2_Strip_log_y",
    "SD1_Strip_log_u", "SD1_Strip_log_v", "SD2_Strip_log_u", "SD2_Strip_log_v",
    "SD1pad_Strip_log_u", "SD1pad_Strip_log_v", "SD2pad_Strip_log_u", "SD2pad_Strip_log_v",
    "MiniFiberD1_Core_log_x1", "MiniFiberD1_Core_log_u1", "MiniFiberD1_Core_log_v1",
    "MiniFiberD1_Core_log_x2", "MiniFiberD1_Core_log_u2", "MiniFiberD1_Core_log_v2",
    "FiberD1_Core_log_x", "FiberD1_Core_log_u", "FiberD1_Core_log_v",
    "FiberD2_Core_log_x", "FiberD2_Core_log_u", "FiberD2_Core_log_v",
    "FiberD3_Core_log_x", "FiberD3_Core_log_u", "FiberD3_Core_log_v",
    "FiberD4_Core_log_x", "FiberD4_Core_log_u", "FiberD4_Core_log_v",
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
void TBuildDetectorLayerPlaneDAF::SelectHists() { LocalHisto.h_stats = AnaHisto->CloneAndRegister(AnaHisto->h_stats); }

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

      std::vector<int> tempSetHit(G4Sol::SIZEOF_G4SOLDETTYPE, -1);
      RecoEvent.TrackDAF.insert(std::make_pair(TrackID, tempSetHit));

      InfoInit tempInit;
      tempInit.charge = event.BeamCharges[index];
      tempInit.posX = event.InteractionPoint_X;
      tempInit.posY = event.InteractionPoint_Y;
      tempInit.posZ = event.InteractionPoint_Z;
      tempInit.momX = event.BeamMomentums_X[index];
      tempInit.momY = event.BeamMomentums_Y[index];
      tempInit.momZ = event.BeamMomentums_Z[index];
      RecoEvent.TrackDAFInit.insert(std::make_pair(TrackID, tempInit));

      std::vector<std::vector<SimHit> > tempSetSimHit(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackDAFSim.insert(std::make_pair(TrackID, tempSetSimHit));

      std::vector<InfoPar> tempSetInfo(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackInfo.insert(std::make_pair(TrackID, tempSetInfo));
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

      std::vector<int> tempSetHit(G4Sol::SIZEOF_G4SOLDETTYPE, -1);
      RecoEvent.TrackDAF.insert(std::make_pair(TrackID, tempSetHit));

      InfoInit tempInit;
      tempInit.charge = event.DaughterCharges[index];
      tempInit.posX = event.DecayVertex_X;
      tempInit.posY = event.DecayVertex_Y;
      tempInit.posZ = event.DecayVertex_Z;
      tempInit.momX = event.DaughterMomentums_X[index];
      tempInit.momY = event.DaughterMomentums_Y[index];
      tempInit.momZ = event.DaughterMomentums_Z[index];
      RecoEvent.TrackDAFInit.insert(std::make_pair(TrackID, tempInit));
      RecoEvent.DaughtersTrackDAFInit.insert(std::make_pair(TrackID, tempInit));

      std::vector<std::vector<SimHit> > tempSetSimHit(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackDAFSim.insert(std::make_pair(TrackID, tempSetSimHit));

      std::vector<InfoPar> tempSetInfo(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackInfo.insert(std::make_pair(TrackID, tempSetInfo));

      RecoEvent.TrackMother.insert(
          std::make_pair(TrackID, std::make_tuple(event.MotherTrackID, event.DecayVertex_X, event.DecayVertex_Y,
                                                  event.DecayVertex_Z, event.DecayTime)));
    }
  OutTree->Nmc = OutTree->fMC_Particle->GetEntries();

  RecoEvent.ListHits.resize(G4Sol::SIZEOF_G4SOLDETTYPE);
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
      auto tempPair                                       = orderDetectors.find(iDet);
      G4Sol::SolDet TypeDet                               = G4Sol::SolDet(tempPair->second);
      std::unique_ptr<genfit::AbsMeasurement> measurement = nullptr;

#ifdef DEBUG_BUILD
      att._logger->debug("iDet # {} {} {}", iDet, nameTempBr, TypeDet);
#endif

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
              auto tempTrack = RecoEvent.TrackDAF.find(TrackID);
              if(tempTrack == RecoEvent.TrackDAF.end())
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

              auto tempTrackInfo                              = RecoEvent.TrackInfo.find(TrackID);
              tempTrackInfo->second[TypeDet + LayerID].pdg    = pdg_code;
              tempTrackInfo->second[TypeDet + LayerID].momX   = hit.MomX;
              tempTrackInfo->second[TypeDet + LayerID].momY   = hit.MomY;
              tempTrackInfo->second[TypeDet + LayerID].momZ   = hit.MomZ;
              tempTrackInfo->second[TypeDet + LayerID].mass   = hit.Mass;
              tempTrackInfo->second[TypeDet + LayerID].Eloss  = hit.Energy;
              tempTrackInfo->second[TypeDet + LayerID].time   = hit.Time;
              tempTrackInfo->second[TypeDet + LayerID].length = hit.TrackLength;

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
              auto tempTrack = RecoEvent.TrackDAF.find(TrackID);
              if(tempTrack == RecoEvent.TrackDAF.end())
                continue;

              TVectorD hitCoordsTree(3);
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
                  std::cout << "LayerID : " << LayerID <<" : " <<TypeDet - G4Sol::MG01 + offsetGeoNameID_PSCE + LayerID - 1<<std::endl;
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
                  TGeoMatrix* g2 = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix(); // INNNER
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
                      std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID, nullptr);
                  dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

                  hitCoordsTree(0) = hit.HitPosX;
                  hitCoordsTree(1) = hit.HitPosY;
                  hitCoordsTree(2) = hit.HitPosZ;
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

                  hitCoordsTree(0) = hit.HitPosX;
                  hitCoordsTree(1) = hit.HitPosY;
                  hitCoordsTree(2) = hit.HitPosZ;
                }
              else if(IsFiberU(TypeDet))
                {
                  //if(IsFiberU_Vetoed(TypeDet))
                  //  continue;
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
                      gGeoManager->GetVolume(volumeName.c_str())->GetNode(LayerID * 2 + 1)->GetMatrix(); // fiber core
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

                  hitCoordsTree(0) = hit.HitPosX;
                  hitCoordsTree(1) = hit.HitPosY;
                  hitCoordsTree(2) = hit.HitPosZ;
                }
              else if(IsFiberM(TypeDet))
                {
                  string volumeName;
                  string motherName;
                  int i_fiber;
                  int i_lay;
                  switch(TypeDet)
                    {
                    case G4Sol::MiniFiberD1_x1: {
                      volumeName = "MiniFiberD1_log_x1";
                      motherName = "MiniFiberD1_log_0";
                      i_fiber = 3;
                      i_lay   = 0; } ;
                      break;
                    case G4Sol::MiniFiberD1_u1: {
                      volumeName = "MiniFiberD1_log_u1";
                      motherName = "MiniFiberD1_log_0";
                      i_fiber = 3;
                      i_lay   = 1; } ;
                      break;
                    case G4Sol::MiniFiberD1_v1: {
                      volumeName = "MiniFiberD1_log_v1";
                      motherName = "MiniFiberD1_log_0";
                      i_fiber = 3;
                      i_lay   = 2; } ;
                      break;
                    case G4Sol::MiniFiberD1_x2: {
                      volumeName = "MiniFiberD1_log_x2";
                      motherName = "MiniFiberD2_log_0";
                      i_fiber = 4;
                      i_lay   = 0; } ;
                      break;
                    case G4Sol::MiniFiberD1_u2: {
                      volumeName = "MiniFiberD1_log_u2";
                      motherName = "MiniFiberD2_log_0";
                      i_fiber = 4;
                      i_lay   = 1; } ;
                      break;
                    case G4Sol::MiniFiberD1_v2: {
                      volumeName = "MiniFiberD1_log_v2";
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
                      gGeoManager->GetVolume(volumeName.c_str())->GetNode(LayerID * 2 + 1)->GetMatrix(); // minifiber core
                  TGeoMatrix* g2 = gGeoManager->GetVolume("MFLD")
                                       ->GetNode(motherName.c_str())
                                       ->GetVolume()
                                       ->GetNode((volumeName + "_0").c_str())
                                       ->GetMatrix(); // minifiber layer
                  TGeoMatrix* g3 =
                      gGeoManager->GetVolume("MFLD")->GetNode(motherName.c_str())->GetMatrix(); // minifiber station
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
                  hitCoords(0) = u.Dot(TVector3(shift[0], shift[1], 0));
                  TMatrixDSym hitCov(1);
                  hitCov(0, 0) = resolution_fiber * resolution_fiber;
                  measurement =
                      std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID, nullptr);
                  dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

                  hitCoordsTree(0) = hit.HitPosX;
                  hitCoordsTree(1) = hit.HitPosY;
                  hitCoordsTree(2) = hit.HitPosZ;
                }
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
                  gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_MDC)->Print();
                  gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + offsetGeoNameID_MDC)->GetMatrix()->Print();
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
              RecoEvent.ListHitsToTracks[TypeDet].emplace_back(TrackID);
              int indexHit = RecoEvent.ListHits[TypeDet].size() - 1;

              tempTrack->second[TypeDet] = indexHit;

              int pdg_code = pid_fromName(hit.Pname);
              if(pdg_code == 0)
                att._logger->debug("!> Builder : pdg_code = 0 ! {}", hit.Pname);

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

              auto tempTrackInfo                    = RecoEvent.TrackInfo.find(TrackID);
              tempTrackInfo->second[TypeDet].pdg    = pdg_code;
              tempTrackInfo->second[TypeDet].momX   = hit.MomX;
              tempTrackInfo->second[TypeDet].momY   = hit.MomY;
              tempTrackInfo->second[TypeDet].momZ   = hit.MomZ;
              tempTrackInfo->second[TypeDet].mass   = hit.Mass;
              tempTrackInfo->second[TypeDet].Eloss  = hit.Energy;
              tempTrackInfo->second[TypeDet].time   = hit.Time;
              tempTrackInfo->second[TypeDet].length = hit.TrackLength;

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

              if(TypeDet >= G4Sol::MiniFiberD1_x1 && TypeDet <= G4Sol::FiberD5_v)
                fillOutHit(OutTree->Fiber, hit, pdg_code, charge, hitCoordsTree, TypeDet, LayerID);

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

  for(size_t i_det = 0; i_det < RecoEvent.ListHits.size(); ++i_det) // const auto& det : RecoEvent.ListHits)
    {
      att._logger->debug("RecoListHits det {} size {}", i_det, RecoEvent.ListHits[i_det].size());
      //   for(const auto& hit : det)
      //     if(hit != nullptr)
      //       hit->Print();
    }
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
