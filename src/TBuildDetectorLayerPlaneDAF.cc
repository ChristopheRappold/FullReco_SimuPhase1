#include "TBuildDetectorLayerPlaneDAF.h"
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <vector>

#include "Debug.hh"

//#define DEBUG_BUILD

using namespace std;

const std::string PDG_fromName::ElName2[] = {"n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",
                                             "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
                                             "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
                                             "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
                                             "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",
                                             "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am",
                                             "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds"};

TBuildDetectorLayerPlaneDAF::TBuildDetectorLayerPlaneDAF(const THyphiAttributes& attribut) : TDataBuilder("build_det"), att(attribut)
{
  cout << "TBuildDetectorLayerPlaneDAF::TBuildDetectorLayerPlaneDAF" << endl;

  std::vector<std::string> tempName = {"HypHI_InSi_log0", "HypHI_InSi_log1", "HypHI_InSi_log2", "HypHI_InSi_log3",
				       "TR1_log","TR2_log",
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
            orderDetectors.insert(std::make_pair(iName, iTypeDet));
        }
    }
}

TBuildDetectorLayerPlaneDAF::~TBuildDetectorLayerPlaneDAF() {}
#ifdef ROOT6
int TBuildDetectorLayerPlaneDAF::operator()(const TG4Sol_Event& event, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits,
                                            FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
{
  int result = Exec(event, hits, RecoEvent, OutTree);

  return SoftExit(result);
}
#else
int TBuildDetectorLayerPlaneDAF::operator()(const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits, FullRecoEvent& RecoEvent,
                                            MCAnaEventG4Sol* OutTree)
{
  int result = Exec(event, hits, RecoEvent, OutTree);

  return SoftExit(result);
}

#endif
int TBuildDetectorLayerPlaneDAF::SoftExit(int return_build)
{
  if(return_build == -1)
    {
#ifndef QUIET
      std::cout << "!> Multiplicity > 2 on Start : event rejected" << std::endl;
#endif
      AnaHisto->h_stats->Fill("start M>2", 1);
      return -1;
    }
  else if(return_build == -2)
    {
#ifndef QUIET
      std::cout << "!> TDC Timing Start cut : event rejected" << std::endl;
#endif
      AnaHisto->h_stats->Fill("start Timing cut", 1);
      return -1;
    }
  else if(return_build == -3)
    {
#ifndef QUIET
      std::cout << "!> Chamber Hit > 1000 : event rejected" << std::endl;
#endif
      AnaHisto->h_stats->Fill("chamber hit>1000", 1);
      return -1;
    }
  else if(return_build == -9)
    {
#ifndef QUIET
      std::cout << "!> No Beam : event rejected" << std::endl;
#endif
      AnaHisto->h_stats->Fill("No Beam", 1);
      return -1;
    }
  else if(return_build != 0)
    {
#ifndef QUIET
      std::cout << "Error in Build Detector !" << std::endl;
#endif
      AnaHisto->h_stats->Fill("Error", 1);
      return -1;
    }
  AnaHisto->h_stats->Fill("start Ok", 1);

  return return_build;
}
#ifdef ROOT6
int TBuildDetectorLayerPlaneDAF::Exec(const TG4Sol_Event& event, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits,
                                      FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
#else
int TBuildDetectorLayerPlaneDAF::Exec(const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits, FullRecoEvent& RecoEvent,
                                      MCAnaEventG4Sol* OutTree)
#endif
{
  int NumFilled = 0;
//  std::cout<<"TBuidlDetectorLayerPlaneDAF : Proton = "<<event.n_Proton<<endl;
#ifdef DEBUG_BUILD
  for(size_t id = 0; id < event.BeamNames.size(); ++id)
    {
      int trackID = event.BeamTrackID[id];
      std::cout << "beam : " << event.BeamNames[id] << " #" << trackID << "\n";

      for(auto& det : hits)
        {
#ifdef ROOT6
          for(auto hit : *det)
            {
              if(hit.TrackID == trackID)
                std::cout << "Branch: " << det->GetBranchName() << " " << hit.HitPosX << " " << hit.HitPosY << " " << hit.HitPosZ << " "
                          << " " << hit.LayerID << " " << hit.Pname << "\n";
            }
#else
          for(int j = 0; j < det->GetEntries(); ++j)
            {
              TG4Sol_Hit* hit = dynamic_cast<TG4Sol_Hit*>(det->At(j));
              if(hit->TrackID == trackID)
                std::cout << "Branch: " << det->GetName() << " hit:" << hit->HitPosX << " " << hit->HitPosY << " " << hit->HitPosZ
                          << " mom:" << hit->MomX << " " << hit->MomY << " " << hit->MomZ << " | "
                          << " " << hit->LayerID << " " << hit->Pname << "\n";
            }
#endif
        }
    }
  for(size_t id = 0; id < event.DaughterNames.size(); ++id)
    {
      int trackID = event.DaughterTrackID[id];
      std::cout << "decayed : " << event.DaughterNames[id] << " #" << trackID << "\n";

      for(auto& det : hits)
        {
#ifdef ROOT6
          for(auto hit : *det)
            {
              // hit.Print();
              if(hit.TrackID == trackID)
                std::cout << "Branch: " << det->GetBranchName() << " " << hit.HitPosX << " " << hit.HitPosY << " " << hit.HitPosZ << " "
                          << hit.LayerID << " " << hit.Pname << "\n";
            }
#else
          for(int j = 0; j < det->GetEntries(); ++j)
            {
              // hit.Print();
              TG4Sol_Hit* hit = dynamic_cast<TG4Sol_Hit*>(det->At(j));
              if(hit->TrackID == trackID)
                std::cout << "Branch: " << det->GetName() << " hit:" << hit->HitPosX << " " << hit->HitPosY << " " << hit->HitPosZ
                          << " mom:" << hit->MomX << " " << hit->MomY << " " << hit->MomZ << " | "
                          << " " << hit->LayerID << " " << hit->Pname << "\n";
            }

#endif
        }
    }
#endif
  
  auto pid_fromName = [] (const std::string& name) {
    if(name=="proton")
      return 2212;
    if(name=="He3")
      return 10003;
    if(name=="pi-")
      return -211;
    if(name=="pi+")
      return 211;
    if(name=="neutron")
      return 2112;
    if(name=="kaon0")
      return 311;
    if(name=="kaon+")
      return 321;
    if(name=="kaon-")
      return -321;
    if(name=="H3L")
      return 20001;
    if(name=="pi0")
      return 111;
    if(name=="triton")
      return 10001;
    if(name=="deuteron")
      return 10000;
    if(name=="alpha")
      return 10002;

    return 0;
  };
  

  std::string nameMother(event.MotherName);
  int id_mother = event.MotherTrackID;

  for(size_t index = 0; index < event.BeamTrackID.size(); ++index)
    {
      TMcParticle* OutParticle = dynamic_cast<TMcParticle*>(OutTree->fMC_Particle->ConstructedAt(OutTree->fMC_Particle->GetEntries()));
      OutParticle->type = event.BeamNames[index];
      OutParticle->Mc_id = event.BeamTrackID[index];
      OutParticle->Mother_id = -1;
      OutParticle->Pdg = pid_fromName(event.BeamNames[index]);
      OutParticle->Charge = event.BeamCharges[index];
      OutParticle->MomMass.SetXYZM(event.BeamMomentums_X[index], event.BeamMomentums_Y[index], event.BeamMomentums_Z[index],
                                   event.BeamMasses[index]);
      OutParticle->Vtx.SetXYZT(event.InteractionPoint_X, event.InteractionPoint_Y, event.InteractionPoint_Z, 0.);
      OutParticle->Weigth = event.BeamCharges[index] == 0 ? 0. : 1.;
      OutParticle->GeoAcc = 1.;

      if(event.BeamNames[index] == nameMother)
        continue;
      int TrackID = event.BeamTrackID[index];

      std::vector<int> tempSetHit(G4Sol::SIZEOF_G4SOLDETTYPE, -1);
      RecoEvent.TrackDAF.insert(std::make_pair(event.BeamTrackID[index], tempSetHit));

      std::vector<SimHit> tempSetSimHit(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackDAFSim.insert(std::make_pair(TrackID, tempSetSimHit));

      std::vector<InfoPar> tempSetInfo(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackInfo.insert(std::make_pair(TrackID, tempSetInfo));
    }
  for(size_t index = 0; index < event.DaughterTrackID.size(); ++index)
    {

      TMcParticle* OutParticle = dynamic_cast<TMcParticle*>(OutTree->fMC_Particle->ConstructedAt(OutTree->fMC_Particle->GetEntries()));
      OutParticle->type = event.DaughterNames[index];
      OutParticle->Mc_id = event.DaughterTrackID[index];
      OutParticle->Mother_id = event.MotherTrackID;
      OutParticle->Pdg = pid_fromName(event.DaughterNames[index]);
      OutParticle->Charge = event.DaughterCharges[index];
      OutParticle->MomMass.SetXYZM(event.DaughterMomentums_X[index], event.DaughterMomentums_Y[index], event.DaughterMomentums_Z[index],
                                   event.DaughterMasses[index]);
      OutParticle->Vtx.SetXYZT(event.DecayVertex_X, event.DecayVertex_Y, event.DecayVertex_Z, event.DecayTime);
      OutParticle->Weigth = event.DaughterCharges[index] == 0 ? 0. : 1.;
      OutParticle->GeoAcc = 1.;

      int TrackID = event.DaughterTrackID[index];

      std::vector<int> tempSetHit(G4Sol::SIZEOF_G4SOLDETTYPE, -1);
      RecoEvent.TrackDAF.insert(std::make_pair(TrackID, tempSetHit));

      std::vector<SimHit> tempSetSimHit(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackDAFSim.insert(std::make_pair(TrackID, tempSetSimHit));

      std::vector<InfoPar> tempSetInfo(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackInfo.insert(std::make_pair(TrackID, tempSetInfo));

      RecoEvent.TrackMother.insert(std::make_pair(TrackID,std::make_tuple(event.MotherTrackID, event.DecayVertex_X, event.DecayVertex_Y, event.DecayVertex_Z, event.DecayTime)));
    }
  OutTree->Nmc = OutTree->fMC_Particle->GetEntries();

  RecoEvent.ListHits.resize(G4Sol::SIZEOF_G4SOLDETTYPE);

  auto fillOutHit = [](TClonesArray* out, const TG4Sol_Hit& hit, int PDG, double charge, const TVectorD& hitR, int LayerID, int HitID) {
    TMcHit* OutHit = dynamic_cast<TMcHit*>(out->ConstructedAt(out->GetEntries()));
    OutHit->name = hit.Pname;
    OutHit->LayerID = LayerID;
    OutHit->HitID = HitID;
    OutHit->MCHit.SetXYZ(hit.HitPosX, hit.HitPosY, hit.HitPosZ);
    OutHit->Hit.SetXYZ(hitR(0), hitR(1), hitR(2));
    OutHit->MC_id = hit.TrackID;
    OutHit->Pdg = PDG;
    OutHit->Charge = charge;
    OutHit->MCparticle.SetXYZM(hit.MomX, hit.MomY, hit.MomZ, hit.Mass);
    OutHit->Brho = 3.10715497 * OutHit->MCparticle.P() / charge;
    OutHit->MagnetInteraction = 0.;
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
      auto tempPair = orderDetectors.find(iDet);
      G4Sol::SolDet TypeDet = G4Sol::SolDet(tempPair->second);
      std::unique_ptr<genfit::AbsMeasurement> measurement = nullptr;

#ifdef DEBUG_BUILD
      std::cout << "iDet #" << iDet << " " << nameTempBr << " " << TypeDet << std::endl;
#endif
      double resolution = 0.01;

      if(nameTempBr == "FMF2_log" || nameTempBr == "HypHI_TrackFwd_log")
        {
#ifdef ROOT6
          for(auto it_hit = tempHits->begin(), it_hit_end = tempHits->end(); it_hit != it_hit_end; ++it_hit)
#else
          for(size_t it_hit = 0; it_hit < tempHits->GetEntries(); ++it_hit)
#endif
            {
#ifdef ROOT6
              auto hit = *it_hit;
              int indexInBranch = std::distance(tempHits->begin(), it_hit);
#else
              const TG4Sol_Hit& hit = *(dynamic_cast<TG4Sol_Hit*>(tempHits->At(it_hit)));
              int indexInBranch = it_hit;
#endif
              int TrackID = hit.TrackID;
              int LayerID = hit.LayerID;
#ifdef DEBUG_BUILD
              std::cout << " hit#" << indexInBranch << " " << hit.Pname << " " << hit.TrackID << " " << hit.LayerID << std::endl;
#endif
              auto tempTrack = RecoEvent.TrackDAF.find(TrackID);
              if(tempTrack == RecoEvent.TrackDAF.end())
                continue;

              TVectorD hitCoords(3);
              hitCoords(0) = gRandom->Gaus(hit.HitPosX, resolution);
              hitCoords(1) = gRandom->Gaus(hit.HitPosY, resolution);
              hitCoords(2) = gRandom->Gaus(hit.HitPosZ, resolution);
              TMatrixDSym hitCov(3);
              hitCov(0, 0) = resolution * resolution;
              hitCov(1, 1) = resolution * resolution;
              hitCov(2, 2) = resolution * resolution;
              measurement = std::make_unique<genfit::SpacepointMeasurement>(hitCoords, hitCov, int(TypeDet) + LayerID, 0, nullptr);

              RecoEvent.ListHits[TypeDet + LayerID].emplace_back(measurement.release());
              int indexHit = RecoEvent.ListHits[TypeDet + LayerID].size() - 1;

              tempTrack->second[TypeDet + LayerID] = indexHit;

              int pdg_code = pid_fromName(hit.Pname);

              auto tempTrackSim = RecoEvent.TrackDAFSim.find(TrackID);
              tempTrackSim->second[TypeDet + LayerID].layerID = LayerID;
              tempTrackSim->second[TypeDet + LayerID].hitX = hit.HitPosX;
              tempTrackSim->second[TypeDet + LayerID].hitY = hit.HitPosY;
              tempTrackSim->second[TypeDet + LayerID].hitZ = hit.HitPosZ;
              tempTrackSim->second[TypeDet + LayerID].momX = hit.MomX;
              tempTrackSim->second[TypeDet + LayerID].momY = hit.MomY;
              tempTrackSim->second[TypeDet + LayerID].momZ = hit.MomZ;
              tempTrackSim->second[TypeDet + LayerID].pdg = pdg_code;
              tempTrackSim->second[TypeDet + LayerID].mass = hit.Mass;
              tempTrackSim->second[TypeDet + LayerID].Eloss = hit.Energy;
              tempTrackSim->second[TypeDet + LayerID].time = gRandom->Gaus(hit.Time, time_res);
              tempTrackSim->second[TypeDet + LayerID].length = hit.TrackLength;

              auto tempTrackInfo = RecoEvent.TrackInfo.find(TrackID);
              tempTrackInfo->second[TypeDet + LayerID].pdg = pdg_code;
              tempTrackInfo->second[TypeDet + LayerID].momX = hit.MomX;
              tempTrackInfo->second[TypeDet + LayerID].momY = hit.MomY;
              tempTrackInfo->second[TypeDet + LayerID].momZ = hit.MomZ;
              tempTrackInfo->second[TypeDet + LayerID].mass = hit.Mass;
              tempTrackInfo->second[TypeDet + LayerID].Eloss = hit.Energy;
              tempTrackInfo->second[TypeDet + LayerID].time = hit.Time;
              tempTrackInfo->second[TypeDet + LayerID].length = hit.TrackLength;

              auto PDG_particle = TDatabasePDG::Instance()->GetParticle(pdg_code);
              if(PDG_particle == nullptr)
                {
                  std::cout << "E> PDG not found !" << std::endl;
                  continue;
                }
              const double charge = PDG_particle->Charge() / 3.;

              if(TypeDet + LayerID >= G4Sol::TrFwd0 && TypeDet + LayerID <= G4Sol::TrFwd2)
                fillOutHit(OutTree->FwdTracker, hit, pdg_code, charge, hitCoords, TypeDet + LayerID, 0.);
              if(TypeDet + LayerID >= G4Sol::FMF2Stop0 && TypeDet + LayerID <= G4Sol::FMF2Stop2)
                fillOutHit(OutTree->FMF2, hit, pdg_code, charge, hitCoords, TypeDet + LayerID, 0.);
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
              auto hit = *it_hit;
              int indexInBranch = std::distance(tempHits->begin(), it_hit);
#else
              const TG4Sol_Hit& hit = *(dynamic_cast<TG4Sol_Hit*>(tempHits->At(it_hit)));
              int indexInBranch = it_hit;
#endif
              int TrackID = hit.TrackID;
              int LayerID = hit.LayerID;
#ifdef DEBUG_BUILD
              std::cout << " hit#" << indexInBranch << " " << hit.Pname << " " << hit.TrackID << " " << hit.LayerID << std::endl;
#endif
              auto tempTrack = RecoEvent.TrackDAF.find(TrackID);
              if(tempTrack == RecoEvent.TrackDAF.end())
                continue;
              TVectorD hitCoords(3);
              hitCoords(0) = gRandom->Gaus(hit.HitPosX, resolution);
              hitCoords(1) = gRandom->Gaus(hit.HitPosY, resolution);
              hitCoords(2) = gRandom->Gaus(hit.HitPosZ, resolution);
              TMatrixDSym hitCov(3);
              hitCov(0, 0) = resolution * resolution;
              hitCov(1, 1) = resolution * resolution;
              hitCov(2, 2) = resolution * resolution;
              measurement = std::make_unique<genfit::SpacepointMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID, nullptr);

              RecoEvent.ListHits[TypeDet].emplace_back(measurement.release());
              int indexHit = RecoEvent.ListHits[TypeDet].size() - 1;

              tempTrack->second[TypeDet] = indexHit;

              int pdg_code = pid_fromName(hit.Pname);

              auto tempTrackSim = RecoEvent.TrackDAFSim.find(TrackID);
              tempTrackSim->second[TypeDet].layerID = LayerID;
              tempTrackSim->second[TypeDet].hitX = hit.HitPosX;
              tempTrackSim->second[TypeDet].hitY = hit.HitPosY;
              tempTrackSim->second[TypeDet].hitZ = hit.HitPosZ;
              tempTrackSim->second[TypeDet].momX = hit.MomX;
              tempTrackSim->second[TypeDet].momY = hit.MomY;
              tempTrackSim->second[TypeDet].momZ = hit.MomZ;
              tempTrackSim->second[TypeDet].pdg = pdg_code;
              tempTrackSim->second[TypeDet].mass = hit.Mass;
              tempTrackSim->second[TypeDet].Eloss = hit.Energy;
              tempTrackSim->second[TypeDet].time = hit.Time;
              tempTrackSim->second[TypeDet].length = hit.TrackLength;

              auto tempTrackInfo = RecoEvent.TrackInfo.find(TrackID);
              tempTrackInfo->second[TypeDet].pdg = pdg_code;
              tempTrackInfo->second[TypeDet].momX = hit.MomX;
              tempTrackInfo->second[TypeDet].momY = hit.MomY;
              tempTrackInfo->second[TypeDet].momZ = hit.MomZ;
              tempTrackInfo->second[TypeDet].mass = hit.Mass;
              tempTrackInfo->second[TypeDet].Eloss = hit.Energy;
              tempTrackInfo->second[TypeDet].time = hit.Time;
              tempTrackInfo->second[TypeDet].length = hit.TrackLength;

              auto PDG_particle = TDatabasePDG::Instance()->GetParticle(pdg_code);
              if(PDG_particle == nullptr)
                {
                  std::cout << "E> PDG not found !" << std::endl;
                  continue;
                }
              const double charge = PDG_particle->Charge() / 3.;

              if(TypeDet >= G4Sol::InSi0 && TypeDet <= G4Sol::InSi3)
                fillOutHit(OutTree->InSi, hit, pdg_code, charge, hitCoords, TypeDet, LayerID);

              if(TypeDet >= G4Sol::TR1 && TypeDet <= G4Sol::TR2)
                fillOutHit(OutTree->TR, hit, pdg_code, charge, hitCoords, TypeDet, LayerID);

              if(TypeDet >= G4Sol::CDC_layer0 && TypeDet <= G4Sol::CDC_layer14)
                fillOutHit(OutTree->CDC, hit, pdg_code, charge, hitCoords, TypeDet, LayerID);

              if(TypeDet >= G4Sol::MG01 && TypeDet <= G4Sol::MG17)
                fillOutHit(OutTree->CDC, hit, pdg_code, charge, hitCoords, TypeDet, LayerID);

              if(TypeDet == G4Sol::CDHBar)
                fillOutHit(OutTree->CDH, hit, pdg_code, charge, hitCoords, TypeDet, LayerID);

              if(TypeDet == G4Sol::RPC_l || TypeDet == G4Sol::RPC_h)
                fillOutHit(OutTree->RPC, hit, pdg_code, charge, hitCoords, TypeDet, LayerID);

              if(TypeDet == G4Sol::PSFE)
                fillOutHit(OutTree->PSFE, hit, pdg_code, charge, hitCoords, TypeDet, LayerID);

              if(TypeDet == G4Sol::PSCE)
                fillOutHit(OutTree->PSCE, hit, pdg_code, charge, hitCoords, TypeDet, LayerID);

              if(TypeDet == G4Sol::PSBE)
                fillOutHit(OutTree->PSBE, hit, pdg_code, charge, hitCoords, TypeDet, LayerID);
            }
        }
    }

  OutTree->Field = att.Field_Strength;

  OutTree->NInSi = OutTree->InSi->GetEntries();
  OutTree->NTr = OutTree->TR->GetEntries();
  OutTree->NCdc = OutTree->CDC->GetEntries();
  OutTree->NCdh = OutTree->CDH->GetEntries();
  OutTree->NRpc = OutTree->RPC->GetEntries();
  OutTree->NFwdtracker = OutTree->FwdTracker->GetEntries();
  OutTree->NFmf2 = OutTree->FMF2->GetEntries();
  OutTree->NPsbe = OutTree->PSBE->GetEntries();
  OutTree->NPsfe = OutTree->PSFE->GetEntries();
  OutTree->NPsce = OutTree->PSCE->GetEntries();

#ifdef DEBUG_BUILD
  std::cout << "done !" << std::endl;
#endif

#ifdef DEBUG_BUILD
  cout << " DAF Hit :" << endl;
  for(auto track : RecoEvent.TrackDAF)
    {
      std::cout << "TrackID #" << track.first << " hit_id [";
      for(auto id_hit : track.second)
        std::cout << " " << id_hit << ", ";
      std::cout << "] " << std::endl;
    }
  for(auto track : RecoEvent.TrackInfo)
    {
      std::cout << "TrackID #" << track.first << " PID [";
      for(auto id_hit : track.second)
        std::cout << " :" << id_hit.pdg << ", ";
      std::cout << "] " << std::endl;
    }
// for(const auto& det : RecoEvent.ListHits)
//   {
//     for(const auto& hit : det)
// 	hit->Print();
//   }
#endif
  return 0;
}
