#include "TBuildDetectorLayerPlaneDAF_ZMQ.h"

#include "Debug.hh"

#include <list>
#include <map>
#include <set>
#include <vector>

//#define DEBUG_BUILD

using namespace std;
using namespace ZMQ;

TBuildDetectorLayerPlaneDAF_ZMQ::TBuildDetectorLayerPlaneDAF_ZMQ(const THyphiAttributes& attribut)
    : TDataBuilder("build_det"), att(attribut)
{
  att._logger->info("TBuildDetectorLayerPlaneDAF_ZMQ::TBuildDetectorLayerPlaneDAF_ZMQ");

  std::vector<std::string> tempName = {"HypHI_InSi_log0",
                                       "HypHI_InSi_log1",
                                       "HypHI_InSi_log2",
                                       "HypHI_InSi_log3",
                                       "TR1_log",
                                       "TR2_log",
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
            orderDetectors.insert(std::make_pair(iName, iTypeDet));
        }
    }
}

TBuildDetectorLayerPlaneDAF_ZMQ::~TBuildDetectorLayerPlaneDAF_ZMQ() {}
#ifdef ROOT6
ReturnRes::InfoM TBuildDetectorLayerPlaneDAF_ZMQ::
operator()(const TG4Sol_Event& event, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits, DataBuilderOut& RecoEvent)
{
  int result = Exec(event, hits, RecoEvent);

  return SoftExit(result);
}
#else
ReturnRes::InfoM TBuildDetectorLayerPlaneDAF_ZMQ::
operator()(const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits, DataBuilderOut& RecoEvent)
{
  int result = Exec(event, hits, RecoEvent);

  return SoftExit(result);
}

#endif

void TBuildDetectorLayerPlaneDAF_ZMQ::SelectHists()
{
  LocalHisto.h_stats = AnaHisto->CloneAndRegister(AnaHisto->h_stats);
}

ReturnRes::InfoM TBuildDetectorLayerPlaneDAF_ZMQ::SoftExit(int return_build)
{
  switch(return_build)
    {
    case 0:
      {
	LocalHisto.h_stats->Fill("start Ok", 1);
	return ReturnRes::Fine;
      }
    case -1:
      {
	att._logger->warn("!> Multiplicity > 2 on Start : event rejected");
	LocalHisto.h_stats->Fill("start M>2", 1);
	return ReturnRes::MultiS2_Start;
      }
    case -2:
      {
	att._logger->warn("!> TDC Timing Start cut : event rejected");
	LocalHisto.h_stats->Fill("start Timing cut", 1);
	return ReturnRes::StartTimingCut;
      }
    case -3:
      {
	att._logger->warn("!> Chamber Hit > 1000 : event rejected");
	LocalHisto.h_stats->Fill("chamber hit>1000", 1);
	return ReturnRes::ChamberHitLimit;
      }
    case -9:
      {
	att._logger->warn("!> No Beam : event rejected");
	LocalHisto.h_stats->Fill("No Beam", 1);
	return ReturnRes::NoBeam;
      }
    default:
      {
	att._logger->warn("Error in Build Detector !");
	LocalHisto.h_stats->Fill("Error", 1);
	return ReturnRes::BuildError;
      }
    }
}
#ifdef ROOT6
int TBuildDetectorLayerPlaneDAF_ZMQ::Exec(const TG4Sol_Event& event,
                                          const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits,
                                          DataBuilderOut& RecoEvent)
#else
int TBuildDetectorLayerPlaneDAF_ZMQ::Exec(const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits,
                                          DataBuilderOut& RecoEvent)
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

  // std::string nameMother(event.MotherName);
  // int id_mother = event.MotherTrackID;
  att._logger->info("Builder> beam size :{}", event.BeamTrackID.size());
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

      RecoEvent.DumpParticles.emplace_back(OParticle);

      // if(event.BeamNames[index] == nameMother)
      //   continue;
      // int TrackID = event.BeamTrackID[index];

      // std::vector<int> tempSetHit(G4Sol::SIZEOF_G4SOLDETTYPE, -1);
      // RecoEvent.TrackDAF.insert(std::make_pair(event.BeamTrackID[index], tempSetHit));

      // std::vector<SimHit> tempSetSimHit(G4Sol::SIZEOF_G4SOLDETTYPE);
      // RecoEvent.TrackDAFSim.insert(std::make_pair(TrackID, tempSetSimHit));

      // std::vector<InfoPar> tempSetInfo(G4Sol::SIZEOF_G4SOLDETTYPE);
      // RecoEvent.TrackInfo.insert(std::make_pair(TrackID, tempSetInfo));
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

      RecoEvent.DumpParticles.emplace_back(OParticle);

      // int TrackID = event.DaughterTrackID[index];

      // std::vector<int> tempSetHit(G4Sol::SIZEOF_G4SOLDETTYPE, -1);
      // RecoEvent.TrackDAF.insert(std::make_pair(TrackID, tempSetHit));

      // std::vector<SimHit> tempSetSimHit(G4Sol::SIZEOF_G4SOLDETTYPE);
      // RecoEvent.TrackDAFSim.insert(std::make_pair(TrackID, tempSetSimHit));

      // std::vector<InfoPar> tempSetInfo(G4Sol::SIZEOF_G4SOLDETTYPE);
      // RecoEvent.TrackInfo.insert(std::make_pair(TrackID, tempSetInfo));

      // RecoEvent.TrackMother.insert(
      //     std::make_pair(TrackID, std::make_tuple(event.MotherTrackID, event.DecayVertex_X, event.DecayVertex_Y,
      //                                             event.DecayVertex_Z, event.DecayTime)));
    }

  // RecoEvent.ListHits.resize(G4Sol::SIZEOF_G4SOLDETTYPE);

  auto fillOutHit = [](std::vector<OutHit>& out, const TG4Sol_Hit& hit, int PDG, double charge, double timeR,
                       const TVectorD& hitR, int LayerID, int HitID) {
    OutHit OHit;
    OHit.name       = hit.Pname;
    OHit.LayerID    = LayerID;
    OHit.HitID      = HitID;
    OHit.MCHit      = {hit.HitPosX, hit.HitPosY, hit.HitPosZ};
    OHit.Hit        = {hitR(0), hitR(1), hitR(2)};
    OHit.Energy     = hit.Energy;
    OHit.MCTime     = hit.Time;
    OHit.Time       = timeR;
    OHit.Length     = hit.TrackLength;
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

  RecoEvent.DumpHits.resize(G4Sol::SIZEOF_G4SOLDETTYPE);

  for(size_t iDet = 0; iDet < hits.size(); ++iDet)
    {
#ifdef ROOT6
      TTreeReaderArray<TG4Sol_Hit>* tempHits = hits[iDet];
      std::string nameTempBr(tempHits->GetBranchName());
#else
      TClonesArray* tempHits = hits[iDet];
      std::string nameTempBr(tempHits->GetName());
#endif
      auto tempPair         = orderDetectors.find(iDet);
      G4Sol::SolDet TypeDet = G4Sol::SolDet(tempPair->second);
      // std::unique_ptr<genfit::AbsMeasurement> measurement = nullptr;

#ifdef DEBUG_BUILD
      att._logger->debug("iDet # {} {} {}", iDet, nameTempBr, TypeDet);
#endif

      double resolution_wire   = 0.01;
      double resolution_wire_z = 0.1;
      double resolution_planar = 0.05;  // cm
      double time_res          = 0.150; // ns
      if(nameTempBr == "FMF2_log" || nameTempBr == "HypHI_TrackFwd_log")
        {
#ifdef ROOT6
          for(auto it_hit = tempHits->begin(), it_hit_end = tempHits->end(); it_hit != it_hit_end; ++it_hit)
#else
          for(int it_hit = 0; it_hit < tempHits->GetEntries(); ++it_hit)
#endif
            {
#ifdef ROOT6
              auto hit = *it_hit;
#else
              const TG4Sol_Hit& hit = *(dynamic_cast<TG4Sol_Hit*>(tempHits->At(it_hit)));
#endif
              int LayerID = hit.LayerID;

#ifdef DEBUG_BUILD
#ifdef ROOT6
              int indexInBranch = std::distance(tempHits->begin(), it_hit);
#else
              int indexInBranch = it_hit;
#endif
              int TrackID = hit.TrackID;
              att._logger->debug(" hit#{} {} {} {}", indexInBranch, hit.Pname, hit.TrackID, hit.LayerID);
#endif
              // auto tempTrack = RecoEvent.TrackDAF.find(TrackID);
              // if(tempTrack == RecoEvent.TrackDAF.end())
              //   continue;

              TVectorD hitCoords(2);
              hitCoords(0) = gRandom->Gaus(hit.HitPosX, resolution_planar);
              hitCoords(1) = gRandom->Gaus(hit.HitPosY, resolution_planar);
              // hitCoords(2) = gRandom->Gaus(hit.HitPosZ, resolution);

              TVectorD hitCoordsTree(3);
              hitCoordsTree(0) = hitCoords(0);
              hitCoordsTree(1) = hitCoords(1);
              hitCoordsTree(2) = hit.HitPosZ;

              // TMatrixDSym hitCov(2);
              // hitCov(0, 0) = resolution_planar * resolution_planar;
              // hitCov(1, 1) = resolution_planar * resolution_planar;
              // // hitCov(2, 2) = resolution * resolution;
              // measurement =
              //     std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet) + LayerID, 0, nullptr);

              // TVector3 o(0., 0., hit.HitPosZ), u(1., 0., 0.), v(0., 1., 0.);
              // genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));
              // dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

              // RecoEvent.ListHits[TypeDet + LayerID].emplace_back(measurement.release());
              // int indexHit = RecoEvent.ListHits[TypeDet + LayerID].size() - 1;

              // tempTrack->second[TypeDet + LayerID] = indexHit;

              int pdg_code = pid_fromName(hit.Pname);
              double time  = gRandom->Gaus(hit.Time, time_res);

              // auto tempTrackSim                               = RecoEvent.TrackDAFSim.find(TrackID);
              // tempTrackSim->second[TypeDet + LayerID].layerID = LayerID;
              // tempTrackSim->second[TypeDet + LayerID].hitX    = hit.HitPosX;
              // tempTrackSim->second[TypeDet + LayerID].hitY    = hit.HitPosY;
              // tempTrackSim->second[TypeDet + LayerID].hitZ    = hit.HitPosZ;
              // tempTrackSim->second[TypeDet + LayerID].momX    = hit.MomX;
              // tempTrackSim->second[TypeDet + LayerID].momY    = hit.MomY;
              // tempTrackSim->second[TypeDet + LayerID].momZ    = hit.MomZ;
              // tempTrackSim->second[TypeDet + LayerID].pdg     = pdg_code;
              // tempTrackSim->second[TypeDet + LayerID].mass    = hit.Mass;
              // tempTrackSim->second[TypeDet + LayerID].Eloss   = hit.Energy;
              // tempTrackSim->second[TypeDet + LayerID].time    = gRandom->Gaus(hit.Time, time_res);
              // tempTrackSim->second[TypeDet + LayerID].length  = hit.TrackLength;

              // auto tempTrackInfo                              = RecoEvent.TrackInfo.find(TrackID);
              // tempTrackInfo->second[TypeDet + LayerID].pdg    = pdg_code;
              // tempTrackInfo->second[TypeDet + LayerID].momX   = hit.MomX;
              // tempTrackInfo->second[TypeDet + LayerID].momY   = hit.MomY;
              // tempTrackInfo->second[TypeDet + LayerID].momZ   = hit.MomZ;
              // tempTrackInfo->second[TypeDet + LayerID].mass   = hit.Mass;
              // tempTrackInfo->second[TypeDet + LayerID].Eloss  = hit.Energy;
              // tempTrackInfo->second[TypeDet + LayerID].time   = hit.Time;
              // tempTrackInfo->second[TypeDet + LayerID].length = hit.TrackLength;

              auto PDG_particle = TDatabasePDG::Instance()->GetParticle(pdg_code);
              if(PDG_particle == nullptr)
                {
                  att._logger->error("E> PDG not found !");
                  continue;
                }
              const double charge = PDG_particle->Charge() / 3.;

              fillOutHit(RecoEvent.DumpHits[TypeDet + LayerID], hit, pdg_code, charge, time, hitCoordsTree,
                         TypeDet + LayerID, 0.);
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
              auto hit = *it_hit;
#else
              const TG4Sol_Hit& hit = *(dynamic_cast<TG4Sol_Hit*>(tempHits->At(it_hit)));
#endif

              int LayerID = hit.LayerID;
#ifdef DEBUG_BUILD
#ifdef ROOT6
              int indexInBranch = std::distance(tempHits->begin(), it_hit);
#else
              int indexInBranch = it_hit;
#endif
              int TrackID = hit.TrackID;
              att._logger->debug(" hit#{} {} {} {}", indexInBranch, hit.Pname, hit.TrackID, hit.LayerID);
#endif
              // auto tempTrack = RecoEvent.TrackDAF.find(TrackID);
              // if(tempTrack == RecoEvent.TrackDAF.end())
              //   continue;

              TVectorD hitCoordsTree(3);
              if(IsPlanar(TypeDet))
                {
                  TVectorD hitCoords(2);
                  hitCoords(0) = gRandom->Gaus(hit.HitPosX, resolution_planar);
                  hitCoords(1) = gRandom->Gaus(hit.HitPosY, resolution_planar);
                  // hitCoords(2) = gRandom->Gaus(hit.HitPosZ, resolution);
                  // TMatrixDSym hitCov(2);
                  // hitCov(0, 0) = resolution_planar * resolution_planar;
                  // hitCov(1, 1) = resolution_planar * resolution_planar;
                  // // hitCov(2, 2) = resolution * resolution;
                  // measurement =
                  //     std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID, nullptr);

                  // TVector3 o(0., 0., hit.HitPosZ), u(1., 0., 0.), v(0., 1., 0.);
                  // genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));
                  // dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

                  hitCoordsTree(0) = hitCoords(0);
                  hitCoordsTree(1) = hitCoords(1);
                  hitCoordsTree(2) = hit.HitPosZ;
                }
              else
                {
                  TVectorD hitCoords(3);

                  hitCoords(0) = gRandom->Gaus(hit.HitPosX, resolution_wire);
                  hitCoords(1) = gRandom->Gaus(hit.HitPosY, resolution_wire);
                  hitCoords(2) = gRandom->Gaus(hit.HitPosZ, resolution_wire_z);

                  // TMatrixDSym hitCov(3);
                  // hitCov(0, 0) = resolution_wire * resolution_wire;
                  // hitCov(1, 1) = resolution_wire * resolution_wire;
                  // hitCov(2, 2) = resolution_wire_z * resolution_wire_z;

                  // measurement = std::make_unique<genfit::ProlateSpacepointMeasurement>(hitCoords, hitCov,
                  // int(TypeDet),
                  //                                                                      LayerID, nullptr);
                  // const TVector3 WireDir(0., 0., 1.);
                  // dynamic_cast<genfit::ProlateSpacepointMeasurement*>(measurement.get())
                  //     ->setLargestErrorDirection(WireDir);

                  hitCoordsTree(0) = hitCoords(0);
                  hitCoordsTree(1) = hitCoords(1);
                  hitCoordsTree(2) = hitCoords(2);
                }

              // RecoEvent.ListHits[TypeDet].emplace_back(measurement.release());
              // int indexHit = RecoEvent.ListHits[TypeDet].size() - 1;

              // tempTrack->second[TypeDet] = indexHit;

              int pdg_code = pid_fromName(hit.Pname);
              double time  = hit.Time;
              // auto tempTrackSim                     = RecoEvent.TrackDAFSim.find(TrackID);
              // tempTrackSim->second[TypeDet].layerID = LayerID;
              // tempTrackSim->second[TypeDet].hitX    = hit.HitPosX;
              // tempTrackSim->second[TypeDet].hitY    = hit.HitPosY;
              // tempTrackSim->second[TypeDet].hitZ    = hit.HitPosZ;
              // tempTrackSim->second[TypeDet].momX    = hit.MomX;
              // tempTrackSim->second[TypeDet].momY    = hit.MomY;
              // tempTrackSim->second[TypeDet].momZ    = hit.MomZ;
              // tempTrackSim->second[TypeDet].pdg     = pdg_code;
              // tempTrackSim->second[TypeDet].mass    = hit.Mass;
              // tempTrackSim->second[TypeDet].Eloss   = hit.Energy;
              // tempTrackSim->second[TypeDet].time    = hit.Time;
              // tempTrackSim->second[TypeDet].length  = hit.TrackLength;

              // auto tempTrackInfo                    = RecoEvent.TrackInfo.find(TrackID);
              // tempTrackInfo->second[TypeDet].pdg    = pdg_code;
              // tempTrackInfo->second[TypeDet].momX   = hit.MomX;
              // tempTrackInfo->second[TypeDet].momY   = hit.MomY;
              // tempTrackInfo->second[TypeDet].momZ   = hit.MomZ;
              // tempTrackInfo->second[TypeDet].mass   = hit.Mass;
              // tempTrackInfo->second[TypeDet].Eloss  = hit.Energy;
              // tempTrackInfo->second[TypeDet].time   = hit.Time;
              // tempTrackInfo->second[TypeDet].length = hit.TrackLength;

              auto PDG_particle = TDatabasePDG::Instance()->GetParticle(pdg_code);
              if(PDG_particle == nullptr)
                {
                  att._logger->error("E> PDG not found !");
                  continue;
                }
              const double charge = PDG_particle->Charge() / 3.;

              fillOutHit(RecoEvent.DumpHits[TypeDet], hit, pdg_code, charge, time, hitCoordsTree, TypeDet, LayerID);
            }
        }
    }

#ifdef DEBUG_BUILD
  att._logger->debug("done !");
#endif

  // #ifdef DEBUG_BUILD
  //   att._logger->debug(" DAF Hit :");
  //   for(auto track : RecoEvent.TrackDAF)
  //     {
  //       att._logger->debug("TrackID # {} hit_id [", track.first);
  //       std::vector<std::stringstream> s1(track.second.size() / 20 + 1);
  //       std::vector<std::stringstream> s2(track.second.size() / 20 + 1);
  //       std::vector<std::stringstream> s3(track.second.size() / 20 + 1);
  //       for(size_t i = 0; i < track.second.size(); ++i)
  //         {
  //           s1[i / 20] << printW(G4Sol::nameLiteralDet.begin()[i], 6) << ", ";
  //           s2[i / 20] << printW(i, 6) << ", ";
  //           s3[i / 20] << printW(track.second[i], 6) << ", ";
  //         }
  //       for(size_t i = 0; i < s1.size(); ++i)
  //         {
  //           att._logger->debug("idDet:{}", s1[i].str());
  //           att._logger->debug("stat :{}", s3[i].str());
  //         }

  //       att._logger->debug("] ");
  //     }

  //   for(auto track : RecoEvent.TrackInfo)
  //     {
  //       att._logger->debug("TrackID #{} PID [", track.first);
  //       std::vector<std::stringstream> s1(track.second.size() / 20 + 1);
  //       std::vector<std::stringstream> s2(track.second.size() / 20 + 1);
  //       std::vector<std::stringstream> s3(track.second.size() / 20 + 1);
  //       for(size_t i = 0; i < track.second.size(); ++i)
  //         {
  //           s1[i / 20] << printW(G4Sol::nameLiteralDet.begin()[i], 6) << ", ";
  //           s2[i / 20] << printW(i, 6) << ", ";
  //           s3[i / 20] << printW(track.second[i].pdg, 6) << ", ";
  //         }
  //       for(size_t i = 0; i < s1.size(); ++i)
  //         {
  //           att._logger->debug("idDet:{}", s1[i].str());
  //           att._logger->debug("stat :{}", s3[i].str());
  //         }

  //       att._logger->debug("] ");
  //     }
  //   for(const auto& det : RecoEvent.ListHits)
  //     {
  //       for(const auto& hit : det)
  //         hit->Print();
  //     }
  // #endif
  return 0;
}
