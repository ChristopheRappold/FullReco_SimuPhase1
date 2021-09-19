#include "TBuildRestarter.h"

#include "Ana_Event/TMcHit.hh"
#include "Ana_Event/TMcParticle.hh"
#include "Ana_Event/TTrackCand.hh"
#include "Debug.hh"
#include "TGeoManager.h"

#include <list>
#include <map>
#include <set>
#include <vector>

//#define DEBUG_BUILD

using namespace std;

TBuildRestarter::TBuildRestarter(const THyphiAttributes& attribut) : TDataBuilder("build_restart"), att(attribut)
{
  att._logger->info("TBuildRestarter::TBuildRestarter");
}

TBuildRestarter::~TBuildRestarter() {}
#ifdef ROOT6
ReturnRes::InfoM TBuildRestarter::operator()(const MCAnaEventG4Sol& event, FullRecoEvent& RecoEvent,
                                             MCAnaEventG4Sol* OutTree)
{
  int result = Exec(event, hits, RecoEvent, OutTree);

  return SoftExit(result);
}
#else
ReturnRes::InfoM TBuildRestarter::operator()(MCAnaEventG4Sol* event, FullRecoEvent& RecoEvent,
                                             MCAnaEventG4Sol* OutTree)
{
  int result = Exec(event, RecoEvent, OutTree);

  return SoftExit(result);
}

#endif
void TBuildRestarter::SelectHists() { LocalHisto.h_stats = AnaHisto->CloneAndRegister(AnaHisto->h_stats); }

ReturnRes::InfoM TBuildRestarter::SoftExit(int return_build)
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
int TBuildRestarter::Exec(const MCAnaEventG4Sol& event, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
#else
int TBuildRestarter::Exec(MCAnaEventG4Sol* event, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
#endif
{
  int NumFilled = 0;
  //  std::cout<<"TBuidlDetectorLayerPlaneDAF : Proton = "<<event->n_Proton<<endl;
  int MotherTrackID = -1;
  for(size_t it_MC = 0; it_MC < event->fMC_Particle->GetEntries(); ++it_MC)
    {
      const TMcParticle& MCpar = *(dynamic_cast<TMcParticle*>(event->fMC_Particle->At(it_MC)));

      if(MCpar.Mother_id == -1)
        {
          RecoEvent.InteractionPoint[0] = MCpar.Vtx.X();
          RecoEvent.InteractionPoint[1] = MCpar.Vtx.Y();
          RecoEvent.InteractionPoint[2] = MCpar.Vtx.Z();
        }
      else
        {
          RecoEvent.DecayVertex[0] = MCpar.Vtx.X();
          RecoEvent.DecayVertex[1] = MCpar.Vtx.Y();
          RecoEvent.DecayVertex[2] = MCpar.Vtx.Z();

          MotherTrackID = MCpar.Mother_id;
        }
    }

  for(size_t it_MC = 0; it_MC < event->fMC_Particle->GetEntries(); ++it_MC)
    {
      const TMcParticle& MCpar = *(dynamic_cast<TMcParticle*>(event->fMC_Particle->At(it_MC)));

      int TrackID = MCpar.Mc_id;

      if(MCpar.Mc_id == MotherTrackID)
        continue;

      if(MCpar.Charge == 0)
	continue;

      std::vector<int> tempSetHit(G4Sol::SIZEOF_G4SOLDETTYPE, -1);
      RecoEvent.TrackDAF.insert(std::make_pair(TrackID, tempSetHit));

      InfoInit tempInit;
      tempInit.posX = MCpar.Vtx.X();
      tempInit.posY = MCpar.Vtx.Y();
      tempInit.posZ = MCpar.Vtx.Z();
      tempInit.momX = MCpar.MomMass.Px();
      tempInit.momY = MCpar.MomMass.Py();
      tempInit.momZ = MCpar.MomMass.Pz();
      RecoEvent.TrackDAFInit.insert(std::make_pair(TrackID, tempInit));
      if(MCpar.Mother_id != -1)
        {
          RecoEvent.TrackMother.insert(std::make_pair(
              TrackID, std::make_tuple(MotherTrackID, MCpar.Vtx.X(), MCpar.Vtx.Y(), MCpar.Vtx.Z(), MCpar.Vtx.T())));

          RecoEvent.DaughtersTrackDAFInit.insert(std::make_pair(TrackID, tempInit));
        }
      std::vector<std::vector<SimHit> > tempSetSimHit(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackDAFSim.insert(std::make_pair(TrackID, tempSetSimHit));

      std::vector<InfoPar> tempSetInfo(G4Sol::SIZEOF_G4SOLDETTYPE);
      RecoEvent.TrackInfo.insert(std::make_pair(TrackID, tempSetInfo));
    }

  RecoEvent.ListHits.resize(G4Sol::SIZEOF_G4SOLDETTYPE);
  RecoEvent.ListHitsToTracks.resize(G4Sol::SIZEOF_G4SOLDETTYPE);
  RecoEvent.OldListHits.resize(G4Sol::SIZEOF_G4SOLDETTYPE);
  RecoEvent.Si_HitsEnergyLayer.resize(4);

  std::unique_ptr<genfit::AbsMeasurement> measurement = nullptr;

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

  for(size_t it_Si = 0; it_Si < event->NInSi; ++it_Si)
    {

      const TMcHit& SiHit = *(dynamic_cast<TMcHit*>(event->InSi->At(it_Si)));
      auto TypeDet        = G4Sol::SolDet(SiHit.LayerID);

      // if(IsSilicon(TypeDet))
      // 	{
      // void simulHitstoSignals(TTreeReaderArray<TG4Sol_Hit>* DetHits, std::vector<std::tuple<double,size_t>>&
      // HitEnergyLayer)
      const double EnergyThreshold = 0.001; // MeV
#ifdef DEBUG_BUILD2
      std::cout << "Silicon \n";
      std::string tempName = orderDetName.find(TypeDet)->second;
      std::cout << " name : " << tempName << "\n";
      std::cout << " LayerID :" << LayerID << "\n";
      std::cout << " Energy :" << hit.Energy << "\n";
#endif
      if(SiHit.Energy < EnergyThreshold)
        continue;

      int idSi      = TypeDet - G4Sol::Si1x;
      auto it_SiHit = RecoEvent.Si_HitsEnergyLayer[idSi].find(SiHit.HitID);
      if(it_SiHit != RecoEvent.Si_HitsEnergyLayer[idSi].end())
        it_SiHit->second += SiHit.Energy;
      else
        RecoEvent.Si_HitsEnergyLayer[idSi].insert(std::make_pair(SiHit.HitID, SiHit.Energy));
    }

  att._logger->debug("NMC   : {} {}",event->Nmc,event->fMC_Particle->GetEntries());
  att._logger->debug("NPSCE : {} {}",event->NPsce,event->PSCE->GetEntries());
  att._logger->debug("NPSBE : {} {}",event->NPsbe,event->PSBE->GetEntries());
  att._logger->debug("NFiber: {} {}",event->NFiber,event->Fiber->GetEntries());
  att._logger->debug("NMDC  : {} {}",event->NCdc,event->CDC->GetEntries());

  for(size_t it_PSCE = 0; it_PSCE < event->PSCE->GetEntries(); ++it_PSCE)
    {

      TMcHit* PSCEHit = dynamic_cast<TMcHit*>(event->PSCE->At(it_PSCE));

      auto TypeDet          = G4Sol::SolDet(PSCEHit->LayerID);
      int TrackID           = PSCEHit->TrackID;

      auto tempTrack = RecoEvent.TrackDAF.find(TrackID);
      if(tempTrack == RecoEvent.TrackDAF.end())
        continue;
        //   else if(IsPSCE(TypeDet))
        // {
#ifdef DEBUG_BUILD2
      std::cout << "PSC" << std::endl;
      std::string tmpName = orderDetName.find(TypeDet)->second;
      std::cout << "name : " << tmpName << std::endl;
      std::cout << "PSCEHit->HitID : " << PSCEHit->HitID << std::endl;
      std::cout << "HitPosX : " << hit.HitPosX << std::endl;
      std::cout << "HitPosY : " << hit.HitPosY << std::endl;
      std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
      gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 4 + PSCEHit->HitID - 1)->Print();
      gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 4 + PSCEHit->HitID - 1)->GetMatrix()->Print();
      gGeoManager->GetVolume("MFLD")->GetNode(0)->Print();
      gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix()->Print();
      gGeoManager->GetVolume("WASA")->GetNode(0)->Print();
      gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix()->Print();
#endif
      TGeoMatrix* g1 =
          gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 4 + PSCEHit->HitID - 1)->GetMatrix(); // PSCE
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
      hitCoords(1) = gRandom->Gaus(PSCEHit->MCHit.Z() - shift[2], resolution_psce_z);
      TMatrixDSym hitCov(2);
      hitCov(0, 0) = resolution_psce * resolution_psce;
      hitCov(1, 1) = resolution_psce_z * resolution_psce_z;
      measurement =
          std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), PSCEHit->HitID, nullptr);
      dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

      RecoEvent.ListHits[TypeDet].emplace_back(measurement.release());
      RecoEvent.ListHitsToTracks[TypeDet].emplace_back(TrackID);
      int indexHit = RecoEvent.ListHits[TypeDet].size() - 1;

      tempTrack->second[TypeDet] = indexHit;

      auto tempTrackSim = RecoEvent.TrackDAFSim.find(TrackID);
      SimHit tempHitSim;
      tempHitSim.layerID = PSCEHit->HitID;
      tempHitSim.hitX    = PSCEHit->MCHit.X();
      tempHitSim.hitY    = PSCEHit->MCHit.Y();
      tempHitSim.hitZ    = PSCEHit->MCHit.Z();
      tempHitSim.momX    = PSCEHit->MCparticle.Px();
      tempHitSim.momY    = PSCEHit->MCparticle.Py();
      tempHitSim.momZ    = PSCEHit->MCparticle.Pz();
      tempHitSim.pdg     = PSCEHit->Pdg;
      tempHitSim.mass    = PSCEHit->MCparticle.M();
      tempHitSim.Eloss   = PSCEHit->Energy;
      tempHitSim.time    = PSCEHit->Time;
      tempHitSim.length  = PSCEHit->TrackLength;
      tempTrackSim->second[TypeDet].emplace_back(tempHitSim);

      auto tempTrackInfo                    = RecoEvent.TrackInfo.find(TrackID);
      tempTrackInfo->second[TypeDet].pdg    = PSCEHit->Pdg;
      tempTrackInfo->second[TypeDet].momX   = PSCEHit->MCparticle.Px();
      tempTrackInfo->second[TypeDet].momY   = PSCEHit->MCparticle.Py();
      tempTrackInfo->second[TypeDet].momZ   = PSCEHit->MCparticle.Pz();
      tempTrackInfo->second[TypeDet].mass   = PSCEHit->MCparticle.M();
      tempTrackInfo->second[TypeDet].Eloss  = PSCEHit->Energy;
      tempTrackInfo->second[TypeDet].time   = PSCEHit->Time;
      tempTrackInfo->second[TypeDet].length = PSCEHit->TrackLength;
    }

  for(size_t it_PSBE = 0; it_PSBE < event->NPsbe; ++it_PSBE)
    {

      const TMcHit& PSBEHit = *(dynamic_cast<TMcHit*>(event->PSBE->At(it_PSBE)));
      auto TypeDet          = G4Sol::SolDet(PSBEHit.LayerID);
      int TrackID           = PSBEHit.TrackID;

      auto tempTrack = RecoEvent.TrackDAF.find(TrackID);
      if(tempTrack == RecoEvent.TrackDAF.end())
        continue;

        // else if(IsPSBE(TypeDet))
        //   {
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
      TGeoMatrix* g1   = gGeoManager->GetVolume("PSB")->GetNode(PSBEHit.HitID - 1)->GetMatrix(); // PSCE
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
          std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), PSBEHit.HitID, nullptr);
      dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

      RecoEvent.ListHits[TypeDet].emplace_back(measurement.release());
      RecoEvent.ListHitsToTracks[TypeDet].emplace_back(TrackID);
      int indexHit = RecoEvent.ListHits[TypeDet].size() - 1;

      tempTrack->second[TypeDet] = indexHit;

      auto tempTrackSim = RecoEvent.TrackDAFSim.find(TrackID);
      SimHit tempHitSim;
      tempHitSim.layerID = PSBEHit.HitID;
      tempHitSim.hitX    = PSBEHit.MCHit.X();
      tempHitSim.hitY    = PSBEHit.MCHit.Y();
      tempHitSim.hitZ    = PSBEHit.MCHit.Z();
      tempHitSim.momX    = PSBEHit.MCparticle.Px();
      tempHitSim.momY    = PSBEHit.MCparticle.Py();
      tempHitSim.momZ    = PSBEHit.MCparticle.Pz();
      tempHitSim.pdg     = PSBEHit.Pdg;
      tempHitSim.mass    = PSBEHit.MCparticle.M();
      tempHitSim.Eloss   = PSBEHit.Energy;
      tempHitSim.time    = PSBEHit.Time;
      tempHitSim.length  = PSBEHit.TrackLength;
      tempTrackSim->second[TypeDet].emplace_back(tempHitSim);

      auto tempTrackInfo                    = RecoEvent.TrackInfo.find(TrackID);
      tempTrackInfo->second[TypeDet].pdg    = PSBEHit.Pdg;
      tempTrackInfo->second[TypeDet].momX   = PSBEHit.MCparticle.Px();
      tempTrackInfo->second[TypeDet].momY   = PSBEHit.MCparticle.Py();
      tempTrackInfo->second[TypeDet].momZ   = PSBEHit.MCparticle.Pz();
      tempTrackInfo->second[TypeDet].mass   = PSBEHit.MCparticle.M();
      tempTrackInfo->second[TypeDet].Eloss  = PSBEHit.Energy;
      tempTrackInfo->second[TypeDet].time   = PSBEHit.Time;
      tempTrackInfo->second[TypeDet].length = PSBEHit.TrackLength;
    }

  for(size_t it_Fiber = 0; it_Fiber < event->NFiber; ++it_Fiber)
    {

      const TMcHit& FiberHit = *(dynamic_cast<TMcHit*>(event->Fiber->At(it_Fiber)));
      auto TypeDet           = G4Sol::SolDet(FiberHit.LayerID);
      int TrackID            = FiberHit.TrackID;

      auto tempTrack = RecoEvent.TrackDAF.find(TrackID);
      if(tempTrack == RecoEvent.TrackDAF.end())
        continue;

      if(restart::IsFiberU(TypeDet))
        {
          if(restart::IsFiberU_Vetoed(TypeDet))
            continue;
          std::string volumeName;
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
            // case G4Sol::MiniFiberD1_x1 :  volumeName = "MiniFiberD1_log_x1"; break;
            // case G4Sol::MiniFiberD1_u1 :  volumeName = "MiniFiberD1_log_u1"; break;
            // case G4Sol::MiniFiberD1_v1 :  volumeName = "MiniFiberD1_log_v1"; break;
            // case G4Sol::MiniFiberD1_x2 :  volumeName = "MiniFiberD1_log_x2"; break;
            // case G4Sol::MiniFiberD1_u2 :  volumeName = "MiniFiberD1_log_u2"; break;
            // case G4Sol::MiniFiberD1_v2 :  volumeName = "MiniFiberD1_log_v2"; break;
            default:
              std::cerr << "something wrong" << std::endl;
              break;
            }
          std::string motherName;
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
            // case G4Sol::MiniFiberD1_x1 :  motherName = "MiniFiberD1_log_0"; break;
            // case G4Sol::MiniFiberD1_u1 :  motherName = "MiniFiberD1_log_0"; break;
            // case G4Sol::MiniFiberD1_v1 :  motherName = "MiniFiberD1_log_0"; break;
            // case G4Sol::MiniFiberD1_x2 :  motherName = "MiniFiberD1_log_0"; break;
            // case G4Sol::MiniFiberD1_u2 :  motherName = "MiniFiberD1_log_0"; break;
            // case G4Sol::MiniFiberD1_v2 :  motherName = "MiniFiberD1_log_0"; break;
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
          gGeoManager->GetVolume(volumeName.c_str())->GetNode(FiberHit.HitID * 2 + 1)->Print();
          gGeoManager->GetVolume(volumeName.c_str())->GetNode(FiberHit.HitID * 2 + 1)->GetMatrix()->Print();
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
              gGeoManager->GetVolume(volumeName.c_str())->GetNode(FiberHit.HitID * 2 + 1)->GetMatrix(); // fiber core
          TGeoMatrix* g2 = gGeoManager->GetVolume("MFLD")
                               ->GetNode(motherName.c_str())
                               ->GetVolume()
                               ->GetNode((volumeName + "_0").c_str())
                               ->GetMatrix();                                                        // fiber layer
          TGeoMatrix* g3 = gGeoManager->GetVolume("MFLD")->GetNode(motherName.c_str())->GetMatrix(); // fiber station
          TGeoMatrix* g4 = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();                  // MFLD
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
              std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), FiberHit.HitID, nullptr);
          dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

          RecoEvent.ListHits[TypeDet].emplace_back(measurement.release());
          RecoEvent.ListHitsToTracks[TypeDet].emplace_back(TrackID);
          int indexHit = RecoEvent.ListHits[TypeDet].size() - 1;

          tempTrack->second[TypeDet] = indexHit;

          auto tempTrackSim = RecoEvent.TrackDAFSim.find(TrackID);
          SimHit tempHitSim;
          tempHitSim.layerID = FiberHit.HitID;
          tempHitSim.hitX    = FiberHit.MCHit.X();
          tempHitSim.hitY    = FiberHit.MCHit.Y();
          tempHitSim.hitZ    = FiberHit.MCHit.Z();
          tempHitSim.momX    = FiberHit.MCparticle.Px();
          tempHitSim.momY    = FiberHit.MCparticle.Py();
          tempHitSim.momZ    = FiberHit.MCparticle.Pz();
          tempHitSim.pdg     = FiberHit.Pdg;
          tempHitSim.mass    = FiberHit.MCparticle.M();
          tempHitSim.Eloss   = FiberHit.Energy;
          tempHitSim.time    = FiberHit.Time;
          tempHitSim.length  = FiberHit.TrackLength;
          tempTrackSim->second[TypeDet].emplace_back(tempHitSim);

          auto tempTrackInfo                    = RecoEvent.TrackInfo.find(TrackID);
          tempTrackInfo->second[TypeDet].pdg    = FiberHit.Pdg;
          tempTrackInfo->second[TypeDet].momX   = FiberHit.MCparticle.Px();
          tempTrackInfo->second[TypeDet].momY   = FiberHit.MCparticle.Py();
          tempTrackInfo->second[TypeDet].momZ   = FiberHit.MCparticle.Pz();
          tempTrackInfo->second[TypeDet].mass   = FiberHit.MCparticle.M();
          tempTrackInfo->second[TypeDet].Eloss  = FiberHit.Energy;
          tempTrackInfo->second[TypeDet].time   = FiberHit.Time;
          tempTrackInfo->second[TypeDet].length = FiberHit.TrackLength;
        }

      else if(restart::IsFiberM(TypeDet))
        {
          string volumeName;
          switch(TypeDet)
            {
            case G4Sol::MiniFiberD1_x1:
              volumeName = "MiniFiberD1_log_x1";
              break;
            case G4Sol::MiniFiberD1_u1:
              volumeName = "MiniFiberD1_log_u1";
              break;
            case G4Sol::MiniFiberD1_v1:
              volumeName = "MiniFiberD1_log_v1";
              break;
            case G4Sol::MiniFiberD1_x2:
              volumeName = "MiniFiberD1_log_x2";
              break;
            case G4Sol::MiniFiberD1_u2:
              volumeName = "MiniFiberD1_log_u2";
              break;
            case G4Sol::MiniFiberD1_v2:
              volumeName = "MiniFiberD1_log_v2";
              break;
            default:
              std::cerr << "something wrong" << std::endl;
              break;
            }
          string motherName;
          switch(TypeDet)
            {
            case G4Sol::MiniFiberD1_x1:
              motherName = "MiniFiberD1_log_0";
              break;
            case G4Sol::MiniFiberD1_u1:
              motherName = "MiniFiberD1_log_0";
              break;
            case G4Sol::MiniFiberD1_v1:
              motherName = "MiniFiberD1_log_0";
              break;
            case G4Sol::MiniFiberD1_x2:
              motherName = "MiniFiberD1_log_0";
              break;
            case G4Sol::MiniFiberD1_u2:
              motherName = "MiniFiberD1_log_0";
              break;
            case G4Sol::MiniFiberD1_v2:
              motherName = "MiniFiberD1_log_0";
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
              gGeoManager->GetVolume(volumeName.c_str())->GetNode(FiberHit.HitID * 2 + 1)->GetMatrix(); // fiber core
          TGeoMatrix* g2 = gGeoManager->GetVolume("MFLD")
                               ->GetNode(motherName.c_str())
                               ->GetVolume()
                               ->GetNode((volumeName + "_0").c_str())
                               ->GetMatrix();                                                        // fiber layer
          TGeoMatrix* g3 = gGeoManager->GetVolume("MFLD")->GetNode(motherName.c_str())->GetMatrix(); // fiber station
          TGeoMatrix* g4 = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();                  // MFLD
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
              std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), FiberHit.HitID, nullptr);
          dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

          RecoEvent.ListHits[TypeDet].emplace_back(measurement.release());
          RecoEvent.ListHitsToTracks[TypeDet].emplace_back(TrackID);
          int indexHit = RecoEvent.ListHits[TypeDet].size() - 1;

          tempTrack->second[TypeDet] = indexHit;

          auto tempTrackSim = RecoEvent.TrackDAFSim.find(TrackID);
          SimHit tempHitSim;
          tempHitSim.layerID = FiberHit.HitID;
          tempHitSim.hitX    = FiberHit.MCHit.X();
          tempHitSim.hitY    = FiberHit.MCHit.Y();
          tempHitSim.hitZ    = FiberHit.MCHit.Z();
          tempHitSim.momX    = FiberHit.MCparticle.Px();
          tempHitSim.momY    = FiberHit.MCparticle.Py();
          tempHitSim.momZ    = FiberHit.MCparticle.Pz();
          tempHitSim.pdg     = FiberHit.Pdg;
          tempHitSim.mass    = FiberHit.MCparticle.M();
          tempHitSim.Eloss   = FiberHit.Energy;
          tempHitSim.time    = FiberHit.Time;
          tempHitSim.length  = FiberHit.TrackLength;
          tempTrackSim->second[TypeDet].emplace_back(tempHitSim);

          auto tempTrackInfo                    = RecoEvent.TrackInfo.find(TrackID);
          tempTrackInfo->second[TypeDet].pdg    = FiberHit.Pdg;
          tempTrackInfo->second[TypeDet].momX   = FiberHit.MCparticle.Px();
          tempTrackInfo->second[TypeDet].momY   = FiberHit.MCparticle.Py();
          tempTrackInfo->second[TypeDet].momZ   = FiberHit.MCparticle.Pz();
          tempTrackInfo->second[TypeDet].mass   = FiberHit.MCparticle.M();
          tempTrackInfo->second[TypeDet].Eloss  = FiberHit.Energy;
          tempTrackInfo->second[TypeDet].time   = FiberHit.Time;
          tempTrackInfo->second[TypeDet].length = FiberHit.TrackLength;
        }
    }

  for(size_t it_MDC = 0; it_MDC < event->NCdc; ++it_MDC)
    {

      const TMcHit& MDCHit = *(dynamic_cast<TMcHit*>(event->CDC->At(it_MDC)));
      auto TypeDet         = G4Sol::SolDet(MDCHit.LayerID);
      int TrackID          = MDCHit.TrackID;

      auto tempTrack = RecoEvent.TrackDAF.find(TrackID);
      if(tempTrack == RecoEvent.TrackDAF.end())
        continue;

        // else if(IsWire(TypeDet)){
#ifdef DEBUG_BUILD2
      std::cout << "wire" << std::endl;
      std::string tmpName = orderDetName.find(TypeDet)->second;
      std::cout << "name : " << tmpName << std::endl;
      std::cout << "LayerID : " << LayerID << std::endl;
      std::cout << "HitPosX : " << hit.HitPosX << std::endl;
      std::cout << "HitPosY : " << hit.HitPosY << std::endl;
      std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
      gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 1)->GetVolume()->GetNode(LayerID - 1)->Print();
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
                           ->GetNode(MDCHit.HitID - 1)
                           ->GetMatrix(); // ME, MG
      TGeoShape* tempShape =
          gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 1)->GetVolume()->GetShape();
      TGeoMatrix* g2 = gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 1)->GetMatrix(); // MD
      TGeoMatrix* g3 = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix();                          // INNER
      TGeoMatrix* g4 = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();                          // MFLD
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
      TVector3 x2(MDCHit.MCHit.X(), MDCHit.MCHit.Y(), MDCHit.MCHit.Z());
      TVector3 p2(MDCHit.MCparticle.Px(), MDCHit.MCparticle.Py(), MDCHit.MCparticle.Pz());
      double dl = restart::CloseDist(x1, x2, p1, p2);

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
      measurement  = std::make_unique<genfit::WireMeasurement>(hitCoords, hitCov, int(TypeDet), MDCHit.HitID, nullptr);
      dynamic_cast<genfit::WireMeasurement*>(measurement.get())->setLeftRightResolution(0);
      dynamic_cast<genfit::WireMeasurement*>(measurement.get())->setMaxDistance(dlmax);

      RecoEvent.ListHits[TypeDet].emplace_back(measurement.release());
      RecoEvent.ListHitsToTracks[TypeDet].emplace_back(TrackID);
      int indexHit = RecoEvent.ListHits[TypeDet].size() - 1;

      tempTrack->second[TypeDet] = indexHit;

      auto tempTrackSim = RecoEvent.TrackDAFSim.find(TrackID);
      SimHit tempHitSim;
      tempHitSim.layerID = MDCHit.HitID;
      tempHitSim.hitX    = MDCHit.MCHit.X();
      tempHitSim.hitY    = MDCHit.MCHit.Y();
      tempHitSim.hitZ    = MDCHit.MCHit.Z();
      tempHitSim.momX    = MDCHit.MCparticle.Px();
      tempHitSim.momY    = MDCHit.MCparticle.Py();
      tempHitSim.momZ    = MDCHit.MCparticle.Pz();
      tempHitSim.pdg     = MDCHit.Pdg;
      tempHitSim.mass    = MDCHit.MCparticle.M();
      tempHitSim.Eloss   = MDCHit.Energy;
      tempHitSim.time    = MDCHit.Time;
      tempHitSim.length  = MDCHit.TrackLength;
      tempTrackSim->second[TypeDet].emplace_back(tempHitSim);

      auto tempTrackInfo                   = RecoEvent.TrackInfo.find(TrackID);
      tempTrackInfo->second[TypeDet].pdg   = MDCHit.Pdg;
      tempTrackInfo->second[TypeDet].momX  = MDCHit.MCparticle.Px();
      tempTrackInfo->second[TypeDet].momY  = MDCHit.MCparticle.Py();
      tempTrackInfo->second[TypeDet].momZ  = MDCHit.MCparticle.Pz();
      tempTrackInfo->second[TypeDet].mass  = MDCHit.MCparticle.M();
      tempTrackInfo->second[TypeDet].Eloss = MDCHit.Energy;
      tempTrackInfo->second[TypeDet].time  = MDCHit.Time;
    }

  for(size_t it_TrFwd = 0; it_TrFwd < event->NFwdtracker; ++it_TrFwd)
    {

      const TMcHit& TrFwdHit = *(dynamic_cast<TMcHit*>(event->FwdTracker->At(it_TrFwd)));
      auto TypeDet           = G4Sol::SolDet(TrFwdHit.LayerID);
      int TrackID            = TrFwdHit.TrackID;

      auto tempTrack = RecoEvent.TrackDAF.find(TrackID);
      if(tempTrack == RecoEvent.TrackDAF.end())
        continue;

      TVectorD hitCoords(2);
      hitCoords(0) = gRandom->Gaus(TrFwdHit.MCHit.X(), resolution_planar);
      hitCoords(1) = gRandom->Gaus(TrFwdHit.MCHit.Y(), resolution_planar);
      // hitCoords(2) = gRandom->Gaus(TrFwdHit.MCHit.Z, resolution);

      TMatrixDSym hitCov(2);
      hitCov(0, 0) = resolution_planar * resolution_planar;
      hitCov(1, 1) = resolution_planar * resolution_planar;
      // hitCov(2, 2) = resolution * resolution;
      // measurement = std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet) + LayerID, 0,
      // nullptr);

      // TVector3 o(0.,0.,TrFwdHit.MCHit.Z), u(1.,0.,0.), v(0.,1.,0.);
      // genfit::SharedPlanePtr plane(new genfit::DetPlane(o,u,v));
      // dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

      // RecoEvent.ListHits[TypeDet + LayerID].emplace_back(measurement.release());
      // int indexHit = RecoEvent.ListHits[TypeDet + LayerID].size() - 1;

      // tempTrack->second[TypeDet + LayerID] = indexHit;

      // RecoEvent.ListHits[TypeDet].emplace_back(measurement.release());
      // RecoEvent.ListHitsToTracks[TypeDet].emplace_back(TrackID);
      // int indexHit = RecoEvent.ListHits[TypeDet].size() - 1;

      // tempTrack->second[TypeDet] = indexHit;

      auto tempTrackSim = RecoEvent.TrackDAFSim.find(TrackID);
      SimHit tempHitSim;
      tempHitSim.layerID = TrFwdHit.HitID;
      tempHitSim.hitX    = TrFwdHit.MCHit.X();
      tempHitSim.hitY    = TrFwdHit.MCHit.Y();
      tempHitSim.hitZ    = TrFwdHit.MCHit.Z();
      tempHitSim.momX    = TrFwdHit.MCparticle.Px();
      tempHitSim.momY    = TrFwdHit.MCparticle.Py();
      tempHitSim.momZ    = TrFwdHit.MCparticle.Pz();
      tempHitSim.pdg     = TrFwdHit.Pdg;
      tempHitSim.mass    = TrFwdHit.MCparticle.M();
      tempHitSim.Eloss   = TrFwdHit.Energy;
      tempHitSim.time    = TrFwdHit.Time;
      tempHitSim.length  = TrFwdHit.TrackLength;
      tempTrackSim->second[TypeDet].emplace_back(tempHitSim);

      auto tempTrackInfo                   = RecoEvent.TrackInfo.find(TrackID);
      tempTrackInfo->second[TypeDet].pdg   = TrFwdHit.Pdg;
      tempTrackInfo->second[TypeDet].momX  = TrFwdHit.MCparticle.Px();
      tempTrackInfo->second[TypeDet].momY  = TrFwdHit.MCparticle.Py();
      tempTrackInfo->second[TypeDet].momZ  = TrFwdHit.MCparticle.Pz();
      tempTrackInfo->second[TypeDet].mass  = TrFwdHit.MCparticle.M();
      tempTrackInfo->second[TypeDet].Eloss = TrFwdHit.Energy;
      tempTrackInfo->second[TypeDet].time  = TrFwdHit.Time;
    }

  for(size_t it_FMF2 = 0; it_FMF2 < event->NFmf2; ++it_FMF2)
    {

      const TMcHit& FMF2Hit = *(dynamic_cast<TMcHit*>(event->FwdTracker->At(it_FMF2)));
      auto TypeDet          = G4Sol::SolDet(FMF2Hit.LayerID);
      int TrackID           = FMF2Hit.TrackID;

      auto tempTrack = RecoEvent.TrackDAF.find(TrackID);
      if(tempTrack == RecoEvent.TrackDAF.end())
        continue;

      TVectorD hitCoords(2);
      hitCoords(0) = gRandom->Gaus(FMF2Hit.MCHit.X(), resolution_planar);
      hitCoords(1) = gRandom->Gaus(FMF2Hit.MCHit.Y(), resolution_planar);
      // hitCoords(2) = gRandom->Gaus(FMF2Hit.MCHit.Z, resolution);

      TMatrixDSym hitCov(2);
      hitCov(0, 0) = resolution_planar * resolution_planar;
      hitCov(1, 1) = resolution_planar * resolution_planar;
      // hitCov(2, 2) = resolution * resolution;
      // measurement = std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet) + LayerID, 0,
      // nullptr);

      // TVector3 o(0.,0.,FMF2Hit.MCHit.Z), u(1.,0.,0.), v(0.,1.,0.);
      // genfit::SharedPlanePtr plane(new genfit::DetPlane(o,u,v));
      // dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

      // RecoEvent.ListHits[TypeDet + LayerID].emplace_back(measurement.release());
      // int indexHit = RecoEvent.ListHits[TypeDet + LayerID].size() - 1;

      // tempTrack->second[TypeDet + LayerID] = indexHit;

      // RecoEvent.ListHits[TypeDet].emplace_back(measurement.release());
      // RecoEvent.ListHitsToTracks[TypeDet].emplace_back(TrackID);
      // int indexHit = RecoEvent.ListHits[TypeDet].size() - 1;

      // tempTrack->second[TypeDet] = indexHit;

      auto tempTrackSim = RecoEvent.TrackDAFSim.find(TrackID);
      SimHit tempHitSim;
      tempHitSim.layerID = FMF2Hit.HitID;
      tempHitSim.hitX    = FMF2Hit.MCHit.X();
      tempHitSim.hitY    = FMF2Hit.MCHit.Y();
      tempHitSim.hitZ    = FMF2Hit.MCHit.Z();
      tempHitSim.momX    = FMF2Hit.MCparticle.Px();
      tempHitSim.momY    = FMF2Hit.MCparticle.Py();
      tempHitSim.momZ    = FMF2Hit.MCparticle.Pz();
      tempHitSim.pdg     = FMF2Hit.Pdg;
      tempHitSim.mass    = FMF2Hit.MCparticle.M();
      tempHitSim.Eloss   = FMF2Hit.Energy;
      tempHitSim.time    = FMF2Hit.Time;
      tempHitSim.length  = FMF2Hit.TrackLength;
      tempTrackSim->second[TypeDet].emplace_back(tempHitSim);

      auto tempTrackInfo                   = RecoEvent.TrackInfo.find(TrackID);
      tempTrackInfo->second[TypeDet].pdg   = FMF2Hit.Pdg;
      tempTrackInfo->second[TypeDet].momX  = FMF2Hit.MCparticle.Px();
      tempTrackInfo->second[TypeDet].momY  = FMF2Hit.MCparticle.Py();
      tempTrackInfo->second[TypeDet].momZ  = FMF2Hit.MCparticle.Pz();
      tempTrackInfo->second[TypeDet].mass  = FMF2Hit.MCparticle.M();
      tempTrackInfo->second[TypeDet].Eloss = FMF2Hit.Energy;
      tempTrackInfo->second[TypeDet].time  = FMF2Hit.Time;
    }


  for(size_t it_trCand = 0; it_trCand < event->NtrackCand;++it_trCand)
    {
      const TTrackCand& TrackC = *(dynamic_cast<TTrackCand*>(OutTree->TrackCand->At(it_trCand)));

      std::vector<int> sortedHits;

      for(size_t it_set = 0; it_set<TrackC.NSet;++it_set)
	{
	  RecoEvent.IdHitsToMeasurement.emplace_back(TrackC.SetLayerID[it_set],TrackC.SetHitID[it_set],TrackC.TrackID);
	  int tempId = RecoEvent.IdHitsToMeasurement.size() - 1;
	  sortedHits.emplace_back(tempId);

	}

      TMatrixD tempCov(5,5,TrackC.Seed_Cov.data());
      std::vector<double> tempPar(TrackC.Seed_Par.begin(), TrackC.Seed_Par.end());
      RecoEvent.TracksFound.emplace_back(sortedHits, TrackC.FitStatus, TrackC.Charge, tempPar, tempCov, TrackC.Chi2_C, TrackC.Chi2_L);

    }




  for(size_t it_track = 0; it_track < event->Ntrack; ++it_track)
    {
      const THyphiTrack& OutTrack = *(dynamic_cast<THyphiTrack*>(OutTree->fTrack->At(it_track)));
      int TrackID                 = -1;
      int Decay                   = 0;
      if(OutTrack.MC_status >= 10000)
        {
          Decay   = 1;
          TrackID = OutTrack.MC_status - 10000;
        }
      else
        TrackID = OutTrack.MC_status;

      ResSolDAF FitRes;

      FitRes.chi2      = OutTrack.Chi2;
      FitRes.ndf       = OutTrack.Chi2_X;
      FitRes.firstHit  = OutTrack.Chi2_Y;
      FitRes.mass      = OutTrack.Mass;
      FitRes.pdg_guess = OutTrack.pdgcode;
      FitRes.pdg_ini   = OutTrack.pdgcode;

      FitRes.momX = OutTrack.MomMass.Px();
      FitRes.momY = OutTrack.MomMass.Py();
      FitRes.momZ = OutTrack.MomMass.Pz();

      FitRes.lastHit = OutTrack.BarId;
      FitRes.charge  = OutTrack.Charge;
      // OutTrack.dE     = TInfo.second[FitRes.lastHit].Eloss;
      FitRes.beta = OutTrack.Beta;
      FitRes.posX = OutTrack.RefPoint.X();
      FitRes.posY = OutTrack.RefPoint.Y();
      FitRes.posZ = OutTrack.RefPoint.Z();

      FitRes.pvalue      = OutTrack.Pval2;
      FitRes.path_length = OutTrack.PathLength;
      FitRes.tof         = OutTrack.TOF;

      FitRes.momX_init = OutTrack.MomIni.X();
      FitRes.momY_init = OutTrack.MomIni.Y();
      FitRes.momZ_init = OutTrack.MomIni.Z();

      FitRes.beta2 = OutTrack.BetaIni;
      FitRes.mass2 = OutTrack.MassIni;
      FitRes.tof2  = OutTrack.TOFIni;

      for(int row = 0; row < 6; row++)
        for(int col = 0; col < 6; col++)
          FitRes.cov_matrix[row][col] = OutTrack.Cov[row][col];

      FitRes.Ncentral = OutTrack.NCent;
      FitRes.Nmfiber  = OutTrack.Nmfib;
      FitRes.iterNum  = OutTrack.iterNum;
      for(int i = 0; i < 17; ++i)
        {
          for(int j = 0; j < 3; ++j)
            {
              FitRes.ResMDC[i][j]    = OutTrack.ResMDC[i][j];
              FitRes.WeightMDC[i][j] = OutTrack.WeightMDC[i][j];
            }
        }
      for(int i = 0; i < 9; ++i)
        {
          FitRes.ResFiber[i]    = OutTrack.ResFiber[i];
          FitRes.WeightFiber[i] = OutTrack.WeightFiber[i];
        }
      for(int i = 0; i < 6; ++i)
        {
          FitRes.ResMiniFiber[i]    = OutTrack.ResMiniFiber[i];
          FitRes.WeightMiniFiber[i] = OutTrack.WeightMiniFiber[i];
        }
      for(int i = 0; i < 2; ++i)
        {
          FitRes.ResPSCE[i] = OutTrack.ResPSCE[i];
        }

      RecoEvent.DAF_results.insert(std::make_pair(TrackID, FitRes));
    }

  return 0;
}

double restart::CloseDist(const TVector3& Xin, const TVector3& Xout, const TVector3& Pin, const TVector3& Pout)
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
