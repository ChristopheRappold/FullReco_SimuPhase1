#include "TBuildDetectorLayerPlaneDAF_MT.h"

#include "Debug.hh"
#include "TGeoManager.h"

#include <list>
#include <map>
#include <set>
#include <vector>

//#define DEBUG_BUILD

using namespace std;

// const std::string PDG_fromName::ElName2[] = {
//     "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
//     "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
//     "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba",
//     "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re",
//     "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu",
//     "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds"};

TBuildDetectorLayerPlaneDAF_MT::TBuildDetectorLayerPlaneDAF_MT(const THyphiAttributes& attribut)
    : TDataBuilder("build_det"), att(attribut)
{
  att._logger->info("TBuildDetectorLayerPlaneDAF_MT::TBuildDetectorLayerPlaneDAF_MT");

  std::vector<std::string> tempName = {"HypHI_InSi_log0",
                                       "HypHI_InSi_log1",
                                       "HypHI_InSi_log2",
                                       "HypHI_InSi_log3",
                                       "TR1_log",
                                       "TR2_log",
                                       "FiberD1_Core_log_x",
                                       "FiberD1_Core_log_u",
                                       "FiberD1_Core_log_v",
                                       "FiberD2_Core_log_x",
                                       "FiberD2_Core_log_u",
                                       "FiberD2_Core_log_v",
                                       "FiberD3_Core_log_x",
                                       "FiberD3_Core_log_u",
                                       "FiberD3_Core_log_v",
                                       "FiberD4_Core_log_x",
                                       "FiberD4_Core_log_u",
                                       "FiberD4_Core_log_v",
                                       "FiberD5_Core_log_x",
                                       "FiberD5_Core_log_u",
                                       "FiberD5_Core_log_v",
                                       "MiniFiberD1_Core_log_x",
                                       "MiniFiberD1_Core_log_u",
                                       "MiniFiberD1_Core_log_v",
                                       "MiniFiberD2_Core_log_x",
                                       "MiniFiberD2_Core_log_u",
                                       "MiniFiberD2_Core_log_v",
                                       "PSFE",
                                       "MG01",
                                       "MG02",
                                       "MG03",
                                       "MG04",
                                       "MG05",
                                       "MG06",
                                       "MG07",
                                       "MG08",
                                       "MG09",
                                       "MG10",
                                       "MG11",
                                       "MG12",
                                       "MG13",
                                       "MG14",
                                       "MG15",
                                       "MG16",
                                       "MG17",
                                       "PSCE",
                                       "PSBE",
                                       "CDC_log0",
                                       "CDC_log1",
                                       "CDC_log2",
                                       "CDC_log3",
                                       "CDC_log4",
                                       "CDC_log5",
                                       "CDC_log6",
                                       "CDC_log7",
                                       "CDC_log8",
                                       "CDC_log9",
                                       "CDC_log10",
                                       "CDC_log11",
                                       "CDC_log12",
                                       "CDC_log13",
                                       "CDC_log14",
                                       "CDH_log",
                                       "HypHI_TrackFwd_log",
                                       "HypHI_TrackFwd_logDummy1",
                                       "HypHI_TrackFwd_logDummy2",
                                       "HypHI_RPC_l_log",
                                       "HypHI_RPC_h_log",
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
}

TBuildDetectorLayerPlaneDAF_MT::~TBuildDetectorLayerPlaneDAF_MT() {}
#ifdef ROOT6
ReturnRes::InfoM TBuildDetectorLayerPlaneDAF_MT::operator()(const TG4Sol_Event& event,
                                                            const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits,
                                                            FullRecoEvent& RecoEvent)
{
  int result = Exec(event, hits, RecoEvent);

  return SoftExit(result);
}
#else
ReturnRes::InfoM TBuildDetectorLayerPlaneDAF_MT::operator()(const TG4Sol_Event& event,
                                                            const std::vector<TClonesArray*>& hits,
                                                            FullRecoEvent& RecoEvent)
{
  int result = Exec(event, hits, RecoEvent);

  return SoftExit(result);
}

#endif

void TBuildDetectorLayerPlaneDAF_MT::SelectHists()
{
  LocalHisto.h_stats = AnaHisto->CloneAndRegister(AnaHisto->h_stats);
}

ReturnRes::InfoM TBuildDetectorLayerPlaneDAF_MT::SoftExit(int return_build)
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
int TBuildDetectorLayerPlaneDAF_MT::Exec(const TG4Sol_Event& event,
                                         const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits,
                                         FullRecoEvent& RecoEvent)
#else
int TBuildDetectorLayerPlaneDAF_MT::Exec(const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits,
                                         FullRecoEvent& RecoEvent)
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

  for(size_t index = 0; index < event.BeamTrackID.size(); ++index)
    {
      OutParticle OParticle;
      OParticle.type      = event.BeamNames[index];
      OParticle.Mc_id     = event.BeamTrackID[index];
      OParticle.Mother_id = -1;
      OParticle.Pdg       = pid_fromName(event.BeamNames[index]);
      OParticle.Charge    = event.BeamCharges[index];
      OParticle.MomMass   = {event.BeamMomentums_X[index], event.BeamMomentums_Y[index], event.BeamMomentums_Z[index],
                           event.BeamMasses[index]};
      OParticle.Vtx       = {event.InteractionPoint_X, event.InteractionPoint_Y, event.InteractionPoint_Z, 0.};
      OParticle.Weigth    = event.BeamCharges[index] == 0 ? 0. : 1.;
      OParticle.GeoAcc    = 1.;
      RecoEvent.ToDumpParticles.emplace_back(OParticle);

      if(event.BeamNames[index] == nameMother)
        continue;
      int TrackID = event.BeamTrackID[index];

      std::vector<int> tempSetHit(G4Sol::SIZEOF_G4SOLDETTYPE, -1);
      RecoEvent.TrackDAF.insert(std::make_pair(event.BeamTrackID[index], tempSetHit));

      std::vector<std::vector<SimHit> > tempSetSimHit(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackDAFSim.insert(std::make_pair(TrackID, tempSetSimHit));

      std::vector<InfoPar> tempSetInfo(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackInfo.insert(std::make_pair(TrackID, tempSetInfo));
    }

  for(size_t index = 0; index < event.DaughterTrackID.size(); ++index)
    {

      OutParticle OParticle;
      OParticle.type      = event.DaughterNames[index];
      OParticle.Mc_id     = event.DaughterTrackID[index];
      OParticle.Mother_id = event.MotherTrackID;
      OParticle.Pdg       = pid_fromName(event.DaughterNames[index]);
      OParticle.Charge    = event.DaughterCharges[index];
      OParticle.MomMass   = {event.DaughterMomentums_X[index], event.DaughterMomentums_Y[index],
                           event.DaughterMomentums_Z[index], event.DaughterMasses[index]};
      OParticle.Vtx       = {event.DecayVertex_X, event.DecayVertex_Y, event.DecayVertex_Z, event.DecayTime};
      OParticle.Weigth    = event.DaughterCharges[index] == 0 ? 0. : 1.;
      OParticle.GeoAcc    = 1.;

      RecoEvent.ToDumpParticles.emplace_back(OParticle);

      int TrackID = event.DaughterTrackID[index];

      std::vector<int> tempSetHit(G4Sol::SIZEOF_G4SOLDETTYPE, -1);
      RecoEvent.TrackDAF.insert(std::make_pair(TrackID, tempSetHit));

      std::vector<std::vector<SimHit> > tempSetSimHit(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackDAFSim.insert(std::make_pair(TrackID, tempSetSimHit));

      std::vector<InfoPar> tempSetInfo(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackInfo.insert(std::make_pair(TrackID, tempSetInfo));

      RecoEvent.TrackMother.insert(
          std::make_pair(TrackID, std::make_tuple(event.MotherTrackID, event.DecayVertex_X, event.DecayVertex_Y,
                                                  event.DecayVertex_Z, event.DecayTime)));
    }

  RecoEvent.ListHits.resize(G4Sol::SIZEOF_G4SOLDETTYPE);

  auto fillOutHit = [](std::vector<OutHit>& out, const TG4Sol_Hit& hit, int PDG, double charge, const TVectorD& hitR,
                       int LayerID, int HitID) {
    OutHit OHit;
    OHit.name       = hit.Pname;
    OHit.LayerID    = LayerID;
    OHit.HitID      = HitID;
    OHit.MCHit      = {hit.HitPosX, hit.HitPosY, hit.HitPosZ};
    OHit.Hit        = {hitR(0), hitR(1), hitR(2)};
    OHit.MC_id      = hit.TrackID;
    OHit.Pdg        = PDG;
    OHit.Charge     = charge;
    OHit.MCparticle = {hit.MomX, hit.MomY, hit.MomZ, hit.Mass};
    OHit.Brho       = 3.10715497 *
                std::sqrt(OHit.MCparticle[0] * OHit.MCparticle[0] + OHit.MCparticle[1] * OHit.MCparticle[1] +
                          OHit.MCparticle[2] * OHit.MCparticle[2]) /
                charge;
    OHit.MagnetInteraction = 0.;

    out.emplace_back(OHit);
    // std::cout<<" Out> LayerID:"<<LayerID<<" "<<HitID<<std::endl;
  };

  RecoEvent.ToDumpHits.resize(G4Sol::SIZEOF_G4SOLDETTYPE);

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
      double resolution_fiber  = 0.0144;
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

              // TVector3 o(0., 0., hit.HitPosZ), u(1., 0., 0.), v(0., 1., 0.);
              // genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));
              // dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

              // RecoEvent.ListHits[TypeDet + LayerID].emplace_back(measurement.release());
              // int indexHit = RecoEvent.ListHits[TypeDet + LayerID].size() - 1;

              // tempTrack->second[TypeDet + LayerID] = indexHit;

              int pdg_code = pid_fromName(hit.Pname);

              auto tempTrackSim                               = RecoEvent.TrackDAFSim.find(TrackID);
	      SimHit tempSimHit;
              tempSimHit.layerID = LayerID;
              tempSimHit.hitX    = hit.HitPosX;
              tempSimHit.hitY    = hit.HitPosY;
              tempSimHit.hitZ    = hit.HitPosZ;
              tempSimHit.momX    = hit.MomX;
              tempSimHit.momY    = hit.MomY;
              tempSimHit.momZ    = hit.MomZ;
              tempSimHit.pdg     = pdg_code;
              tempSimHit.mass    = hit.Mass;
              tempSimHit.Eloss   = hit.Energy;
              tempSimHit.time    = gRandom->Gaus(hit.Time, time_res);
              tempSimHit.length  = hit.TrackLength;
	      tempTrackSim->second[TypeDet + LayerID].emplace_back(tempSimHit);
				   
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

              fillOutHit(RecoEvent.ToDumpHits[TypeDet + LayerID], hit, pdg_code, charge, hitCoordsTree,
                         TypeDet + LayerID, 0.);
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
              // if(IsPlanar(TypeDet))
              //   {
              //     TVectorD hitCoords(2);
              //     hitCoords(0) = gRandom->Gaus(hit.HitPosX, resolution_planar);
              //     hitCoords(1) = gRandom->Gaus(hit.HitPosY, resolution_planar);
              //     // hitCoords(2) = gRandom->Gaus(hit.HitPosZ, resolution);
              //     TMatrixDSym hitCov(2);
              //     hitCov(0, 0) = resolution_planar * resolution_planar;
              //     hitCov(1, 1) = resolution_planar * resolution_planar;
              //     // hitCov(2, 2) = resolution * resolution;
              //     measurement =
              //         std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID, nullptr);

              //     TVector3 o(0., 0., hit.HitPosZ), u(1., 0., 0.), v(0., 1., 0.);
              //     genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));
              //     dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

              //     hitCoordsTree(0) = hitCoords(0);
              //     hitCoordsTree(1) = hitCoords(1);
              //     hitCoordsTree(2) = hit.HitPosZ;
              //   }

              if(IsPSCE(TypeDet))
                {
#ifdef DEBUG_BUILD
                  std::cout << "PSB" << std::endl;
                  std::string tmpName = orderDetName.find(TypeDet)->second;
                  std::cout << "name : " << tmpName << std::endl;
                  std::cout << "LayerID : " << LayerID << std::endl;
                  std::cout << "HitPosX : " << hit.HitPosX << std::endl;
                  std::cout << "HitPosY : " << hit.HitPosY << std::endl;
                  std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
                  gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 4 + LayerID - 1)->Print();
                  gGeoManager->GetVolume("INNER")
                      ->GetNode(TypeDet - G4Sol::MG01 + 4 + LayerID - 1)
                      ->GetMatrix()
                      ->Print();
                  gGeoManager->GetVolume("MFLD")->GetNode(0)->Print();
                  gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix()->Print();
                  gGeoManager->GetVolume("WASA")->GetNode(0)->Print();
                  gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix()->Print();
#endif
                  TGeoMatrix* g1 = gGeoManager->GetVolume("INNER")
                                       ->GetNode(TypeDet - G4Sol::MG01 + 4 + LayerID - 1)
                                       ->GetMatrix();                                       // PSCE
                  TGeoMatrix* g2 = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix(); // INNNER
                  TGeoMatrix* g3 = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix(); // MFLD
                  TGeoHMatrix H1(*g1), H2(*g2), H3(*g3);
                  TGeoHMatrix H = H2 * H1;
                  H             = H3 * H;
#ifdef DEBUG_BUILD
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
              else if(IsFiberU(TypeDet))
                {
                  if(IsFiberU_Vetoed(TypeDet))
                    continue;
                  string volumeName;
                  switch(TypeDet)
                    {
                    case G4Sol::FiberD1_x:
                      volumeName = "FiberD1_log_x";
                      break;
                    case G4Sol::FiberD1_u:
                      volumeName = "FiberD1_log_u";
                      break;
                    case G4Sol::FiberD1_v:
                      volumeName = "FiberD1_log_v";
                      break;
                    case G4Sol::FiberD2_x:
                      volumeName = "FiberD2_log_x";
                      break;
                    case G4Sol::FiberD2_u:
                      volumeName = "FiberD2_log_u";
                      break;
                    case G4Sol::FiberD2_v:
                      volumeName = "FiberD2_log_v";
                      break;
                    case G4Sol::FiberD3_x:
                      volumeName = "FiberD3_log_x";
                      break;
                    case G4Sol::FiberD3_u:
                      volumeName = "FiberD3_log_u";
                      break;
                    case G4Sol::FiberD3_v:
                      volumeName = "FiberD3_log_v";
                      break;
                    case G4Sol::FiberD4_x:
                      volumeName = "FiberD4_log_x";
                      break;
                    case G4Sol::FiberD4_u:
                      volumeName = "FiberD4_log_u";
                      break;
                    case G4Sol::FiberD4_v:
                      volumeName = "FiberD4_log_v";
                      break;
                    case G4Sol::FiberD5_x:
                      volumeName = "FiberD5_log_x";
                      break;
                    case G4Sol::FiberD5_u:
                      volumeName = "FiberD5_log_u";
                      break;
                    case G4Sol::FiberD5_v:
                      volumeName = "FiberD5_log_v";
                      break;
                    case G4Sol::MiniFiberD1_x:
                      volumeName = "MiniFiberD1_log_x";
                      break;
                    case G4Sol::MiniFiberD1_u:
                      volumeName = "MiniFiberD1_log_u";
                      break;
                    case G4Sol::MiniFiberD1_v:
                      volumeName = "MiniFiberD1_log_v";
                      break;
                    case G4Sol::MiniFiberD2_x:
                      volumeName = "MiniFiberD2_log_x";
                      break;
                    case G4Sol::MiniFiberD2_u:
                      volumeName = "MiniFiberD2_log_u";
                      break;
                    case G4Sol::MiniFiberD2_v:
                      volumeName = "MiniFiberD2_log_v";
                      break;
                    default:
                      std::cerr << "something wrong" << std::endl;
                      break;
                    }
                  string motherName;
                  switch(TypeDet)
                    {
                    case G4Sol::FiberD1_x:
                      motherName = "FiberD1_log_0";
                      break;
                    case G4Sol::FiberD1_u:
                      motherName = "FiberD1_log_0";
                      break;
                    case G4Sol::FiberD1_v:
                      motherName = "FiberD1_log_0";
                      break;
                    case G4Sol::FiberD2_x:
                      motherName = "FiberD2_log_0";
                      break;
                    case G4Sol::FiberD2_u:
                      motherName = "FiberD2_log_0";
                      break;
                    case G4Sol::FiberD2_v:
                      motherName = "FiberD2_log_0";
                      break;
                    case G4Sol::FiberD3_x:
                      motherName = "FiberD3_log_0";
                      break;
                    case G4Sol::FiberD3_u:
                      motherName = "FiberD3_log_0";
                      break;
                    case G4Sol::FiberD3_v:
                      motherName = "FiberD3_log_0";
                      break;
                    case G4Sol::FiberD4_x:
                      motherName = "FiberD4_log_0";
                      break;
                    case G4Sol::FiberD4_u:
                      motherName = "FiberD4_log_0";
                      break;
                    case G4Sol::FiberD4_v:
                      motherName = "FiberD4_log_0";
                      break;
                    case G4Sol::FiberD5_x:
                      motherName = "FiberD5_log_0";
                      break;
                    case G4Sol::FiberD5_u:
                      motherName = "FiberD5_log_0";
                      break;
                    case G4Sol::FiberD5_v:
                      motherName = "FiberD5_log_0";
                      break;
                    case G4Sol::MiniFiberD1_x:
                      motherName = "MiniFiberD1_log_0";
                      break;
                    case G4Sol::MiniFiberD1_u:
                      motherName = "MiniFiberD1_log_0";
                      break;
                    case G4Sol::MiniFiberD1_v:
                      motherName = "MiniFiberD1_log_0";
                      break;
                    case G4Sol::MiniFiberD2_x:
                      motherName = "MiniFiberD2_log_0";
                      break;
                    case G4Sol::MiniFiberD2_u:
                      motherName = "MiniFiberD2_log_0";
                      break;
                    case G4Sol::MiniFiberD2_v:
                      motherName = "MiniFiberD2_log_0";
                      break;
                    default:
                      std::cerr << "something wrong" << std::endl;
                      break;
                    }
#ifdef DEBUG_BUILD
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
#ifdef DEBUG_BUILD
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
                  hitCoords(0) = u.Dot(TVector3(shift[0], 0, 0));
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
#ifdef DEBUG_BUILD
                  std::cout << "wire" << std::endl;
                  std::string tmpName = orderDetName.find(TypeDet)->second;
                  std::cout << "name : " << tmpName << std::endl;
                  std::cout << "LayerID : " << LayerID << std::endl;
                  std::cout << "HitPosX : " << hit.HitPosX << std::endl;
                  std::cout << "HitPosY : " << hit.HitPosY << std::endl;
                  std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
                  gGeoManager->GetVolume("INNER")
                      ->GetNode(TypeDet - G4Sol::MG01 + 1)
                      ->GetVolume()
                      ->GetNode(LayerID - 1)
                      ->Print();
                  gGeoManager->GetVolume("INNER")
                      ->GetNode(TypeDet - G4Sol::MG01 + 1)
                      ->GetVolume()
                      ->GetNode(LayerID - 1)
                      ->GetMatrix()
                      ->Print();
                  gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 1)->Print();
                  gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 1)->GetMatrix()->Print();
                  gGeoManager->GetVolume("MFLD")->GetNode(0)->Print();
                  gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix()->Print();
                  gGeoManager->GetVolume("WASA")->GetNode(0)->Print();
                  gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix()->Print();
#endif
                  TGeoMatrix* g1 = gGeoManager->GetVolume("INNER")
                                       ->GetNode(TypeDet - G4Sol::MG01 + 1)
                                       ->GetVolume()
                                       ->GetNode(LayerID - 1)
                                       ->GetMatrix(); // ME, MG
                  TGeoMatrix* g2 =
                      gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 1)->GetMatrix(); // MD
                  TGeoMatrix* g3 = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix();             // INNER
                  TGeoMatrix* g4 = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();             // MFLD
                  TGeoHMatrix H1(*g1), H2(*g2), H3(*g3), H4(*g4);
                  TGeoHMatrix H = H2 * H1;
                  H             = H3 * H;
                  H             = H4 * H;
                  double* shift = H.GetTranslation();
                  TGeoHMatrix w1("w1");
                  TGeoHMatrix w2("w2");
                  w1.SetDz(-10);
                  w2.SetDz(10);
                  TGeoHMatrix Hw1 = H * w1;
                  TGeoHMatrix Hw2 = H * w2;
#ifdef DEBUG_BUILD
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
                  double dl    = CloseDist(x1, x2, p1, p2);
                  dl           = gRandom->Gaus(dl, resolution_dl);
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
                  if(dl < 0)
                    dl = 0;
                  if(dl > dlmax)
                    dl = dlmax;

                  TVectorD hitCoords(7);
                  hitCoords(0) = edge1[0];
                  hitCoords(1) = edge1[1];
                  hitCoords(2) = edge1[2];
                  hitCoords(3) = edge2[0];
                  hitCoords(4) = edge2[1];
                  hitCoords(5) = edge2[2];
                  hitCoords(6) = dl;
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
                  TVectorD hitCoords(3);

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
              int indexHit = RecoEvent.ListHits[TypeDet].size() - 1;

              tempTrack->second[TypeDet] = indexHit;

              int pdg_code = pid_fromName(hit.Pname);

              auto tempTrackSim                     = RecoEvent.TrackDAFSim.find(TrackID);
	      SimHit tempSimHit;
              tempSimHit.layerID = LayerID;
              tempSimHit.hitX    = hit.HitPosX;
              tempSimHit.hitY    = hit.HitPosY;
              tempSimHit.hitZ    = hit.HitPosZ;
              tempSimHit.momX    = hit.MomX;
              tempSimHit.momY    = hit.MomY;
              tempSimHit.momZ    = hit.MomZ;
              tempSimHit.pdg     = pdg_code;
              tempSimHit.mass    = hit.Mass;
              tempSimHit.Eloss   = hit.Energy;
              tempSimHit.time    = hit.Time;
              tempSimHit.length  = hit.TrackLength;
              tempTrackSim->second[TypeDet].emplace_back(tempSimHit);

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

              fillOutHit(RecoEvent.ToDumpHits[TypeDet], hit, pdg_code, charge, hitCoordsTree, TypeDet, LayerID);
            }
        }
    }

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
  for(const auto& det : RecoEvent.ListHits)
    {
      for(const auto& hit : det)
        hit->Print();
    }
#endif
  return 0;
}

double TBuildDetectorLayerPlaneDAF_MT::CloseDist(const TVector3& Xin, const TVector3& Xout, const TVector3& Pin, const TVector3& Pout)
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