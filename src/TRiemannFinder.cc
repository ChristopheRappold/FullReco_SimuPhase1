#include "TRiemannFinder.h"

#include "FullRecoEvent.hh"
#include "KalmanFittedStateOnPlane.h"
#include "KalmanFitterInfo.h"
#include "PlanarMeasurement.h"
#include "ReturnRes.hh"
#include "StateOnPlane.h"

#include "tricktrack/RiemannFit.h"

#include <set>
#include <sstream>
#include <tuple>

#define DEBUG_RIEMANNFINDER

using namespace std;
using namespace G4Sol;

TRiemannFinder::TRiemannFinder(const THyphiAttributes& attribut) : TDataProcessInterface("RiemannFinder"), att(attribut)
{

  std::string temp_name_out       = att.Config.Get<std::string>("Output_Namefile");
  std::string temp_file_base_name = temp_name_out.substr(0, temp_name_out.find_last_of('.'));

  temp_file_base_name += "FinderTest.root";
  namefileFinder = temp_file_base_name;

  f_finder = new TFile(namefileFinder, "RECREATE");
  t_finder = new TTree("RiemannFinder", "Finder output tracks");

  t_finder->SetDirectory(f_finder);
  // t_finder->AutoSave();

  h_XY          = new TH2I("XY", "XY", 200, -30, 30, 200, -30, 30);
  h_XYSim       = new TH2I("XYSim", "XYSim", 200, -30, 30, 200, -30, 30);
  h_XYReco      = new TH2I("XYReco", "XYReco", 200, -30, 30, 200, -30, 30);
  h_RPhi        = new TH2I("RPhi", "RPhi", 200, 0, 30, 360, -180, 180);
  h_RPhiSim     = new TH2I("RPhiSim", "RPhiSim", 200, 0, 30, 360, -180, 180);
  h_XYConformal = new TH2I("XYConformal", "XYConformal", 300, -0.3, 0.3, 300, -0.3, 0.3);
  // h_XYRiemann = new TH2I("XYRiemann","XYRiemann",300,-0.3,0.3,300,-0.3,0.3);
  h_XYRiemann       = new TH2I("XYRiemann", "XYRiemann", 300, 0, 0.5, 360, 0, 360);
  h_XYTracks        = new TH2I("XYTracks", "XYTracks", 200, -30, 30, 200, -30, 30);
  h_XYRiemannTracks = new TH2I("XYRiemannTracks", "XYRiemannTracks", 300, -.3, .3, 300, -0.3, 0.3);

  t_finder->Branch("XY", "TH2I", &h_XY, 12800, 0);
  t_finder->Branch("XYSim", "TH2I", &h_XYSim, 12800, 0);
  t_finder->Branch("XYReco", "TH2I", &h_XYReco, 12800, 0);
  t_finder->Branch("RPhi", "TH2I", &h_RPhi, 12800, 0);
  t_finder->Branch("RPhiSim", "TH2I", &h_RPhiSim, 12800, 0);
  t_finder->Branch("XYConformal", "TH2I", &h_XYConformal, 12800, 0);
  t_finder->Branch("XYRiemann", "TH2I", &h_XYRiemann, 12800, 0);
  t_finder->Branch("XYTracks", "TH2I", &h_XYTracks, 12800, 0);
  t_finder->Branch("XYRiemannTracks", "TH2I", &h_XYRiemannTracks, 12800, 0);

  LastFrontWall = att.Wasa_Side == 0 ? G4Sol::PSFE : G4Sol::PSBE;

  att._logger->info("RiemannFinder : tree out set ");

}

TRiemannFinder::~TRiemannFinder()
{
  f_finder->cd();
  t_finder->Write();
  // f_finder->Write();
  f_finder->Close();

  if(f_finder != nullptr)
    {
      f_finder->Delete();
      f_finder = nullptr;
    }
}

void TRiemannFinder::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }

ReturnRes::InfoM TRiemannFinder::operator()(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
{

  int result_finder = Exec(RecoEvent, OutTree);

  return SoftExit(result_finder);
}

int TRiemannFinder::Exec(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree) { return FinderTrack(RecoEvent); }

ReturnRes::InfoM TRiemannFinder::SoftExit(int result_full) { return ReturnRes::Fine; }

void TRiemannFinder::SelectHists() {}


void TRiemannFinder::BuildKDTree(const FullRecoEvent& RecoEvent, tricktrack::FKDTree<RPhiHit, double, 3>& KDtree, std::vector<RPhiHit>& TempHits)
{
  int ntrack = -1;
  unsigned int nCount = 0;

  for(auto it_trackDAF : RecoEvent.TrackDAF)
    {
      ntrack++;

      const int id_track  = it_trackDAF.first;
      auto it_ListHits    = it_trackDAF.second;
      auto it_ListHitsSim = RecoEvent.TrackDAFSim.find(id_track);
      auto it_Info        = RecoEvent.TrackInfo.find(id_track);
#ifdef DEBUG_RIEMANNFINDER
      att._logger->debug("track {} id {}", ntrack, id_track);
#endif
      std::vector<std::array<double, 3> > SimHits;
      int FirstId = -1;

      for(int id_det = 0; id_det < it_ListHits.size(); ++id_det)
        {
          int id_hit = it_ListHits[id_det];

          if(id_hit < 0)
            continue;
          if(FirstId == -1)
            {
	      FirstId = id_det;
#ifdef DEBUG_RIEMANNFINDER
	      att._logger->debug("pid {}",it_Info->second[FirstId].pdg);
#endif
	    }
          if(id_det >= G4Sol::MG01 && id_det <= G4Sol::MG17)
            {

              auto SimHit = it_ListHitsSim->second[id_det][0];
              h_XYSim->Fill(SimHit.hitX, SimHit.hitY, id_track);
              h_RPhiSim->Fill(TMath::Sqrt(TMath::Sq(SimHit.hitX) + TMath::Sq(SimHit.hitY)),
                              TMath::ATan2(SimHit.hitY, SimHit.hitX) * TMath::RadToDeg(), id_track );

              SimHits.emplace_back(std::array<double, 3>{SimHit.hitX, SimHit.hitY, SimHit.hitZ});

              genfit::WireMeasurement* currentHit =
                  dynamic_cast<genfit::WireMeasurement*>(RecoEvent.ListHits[id_det][id_hit].get());

              auto HitCoord = currentHit->getRawHitCoords();
              auto HitCov   = currentHit->getRawHitCov();

              TVector3 wire_side1(HitCoord[0], HitCoord[1], HitCoord[2]);
              TVector3 wire_side2(HitCoord[3], HitCoord[4], HitCoord[5]);
              //double DriftDist    = HitCoord[6];
              //double ErrDriftDist = HitCov[6][6];
              const int idMDC     = id_det - G4Sol::MG01;
              double R_size       = dl_max(idMDC) * 0.5;

#ifdef DEBUG_RIEMANNFINDER
              double Pt    = TMath::Sqrt(TMath::Sq(SimHit.momX) + TMath::Sq(SimHit.momY));
	      double Rd = TMath::Hypot(SimHit.hitX, SimHit.hitY);
	      double Phid = TMath::ATan2(SimHit.hitY, SimHit.hitX) * TMath::RadToDeg();
              double Phi   = TMath::ATan2(SimHit.momY, SimHit.momX) * TMath::RadToDeg();
              double Theta = TMath::ATan2(Pt, SimHit.momZ) * TMath::RadToDeg();
              att._logger->debug(" layer ID : {} / Hit Id {} / Pt : {} / Phi : {} / Theta : {} | R {} Phi {} | ID {}", idMDC, id_hit, Pt, Phi, Theta, Rd, Phid, SimHit.layerID);
#endif
              auto f_propagWire = [](const TVector3& w1, const TVector3& w2, double relativePos, TVector3& Pos) {
                TVector3 Dir = w2 - w1;
                Dir *= 1. / (w2.Z() - w1.Z());
                double posZ = Pos.Z();
                if(TMath::Abs(relativePos) <= 1.)
                  posZ = w2.Z() * (relativePos) + w1.Z() * (1 - relativePos);

                Pos = Dir * (posZ - w1.Z()) + w1;
              };

              double meanX = 0.5 * (wire_side1.X() + wire_side2.X());
              double meanY = 0.5 * (wire_side1.Y() + wire_side2.Y());
              if(TMath::Abs(meanX - wire_side1.X()) > 0.01)
                {
                  TVector3 midPos(0., 0., 0.);
                  // double RelPos = idMDC / static_cast<double>(G4Sol::MG17);
                  f_propagWire(wire_side1, wire_side2, 0., midPos);
                  meanX = midPos.X();
                  meanY = midPos.Y();
                  redo_layer.insert(id_det);
		  redo_layerToTrackID[id_det][id_hit] = id_track;
		  continue;
                }

              double PhiXY1 = wire_side1.Phi() * TMath::RadToDeg();
              double PhiXY2 = wire_side2.Phi() * TMath::RadToDeg();

              double RXY1 = wire_side1.Perp();
              double RXY2 = wire_side2.Perp();

              double dPhiXY1 = TMath::ASin(R_size / RXY1);
              //double dPhiXY2 = TMath::ASin(R_size / RXY2);
              double dRXY1   = R_size;
              //double dRXY2   = R_size;

              // if(TMath::Abs(PhiXY1-PhiXY2)>0.1)
              // 	att._logger->warn("!> Parallel wire but different Phi ! {} {}",PhiXY1, PhiXY2);

              double R = TMath::Sq(meanX) + TMath::Sq(meanY);
              h_XY->Fill(meanX, meanY, id_track);
              h_XYConformal->Fill(meanX / R, meanY / R, id_track );

              double Xp = TMath::Sqrt(R) * TMath::Cos(PhiXY1 * TMath::DegToRad()) / (1. + R);
              double Yp = TMath::Sqrt(R) * TMath::Sin(PhiXY1 * TMath::DegToRad()) / (1. + R);
              double Zp = R / (1. + R);

              h_XYRiemann->Fill(TMath::Sqrt(Xp * Xp + Yp * Yp), (TMath::ATan2(Yp, Xp)+TMath::Pi()) * TMath::RadToDeg(),
                                id_track);

              TempHits.emplace_back(
				    RPhiHit(TMath::Sqrt(Xp * Xp + Yp * Yp), (TMath::ATan2(Yp, Xp)+TMath::Pi()) * TMath::RadToDeg(), Zp, nCount));
              ++nCount;
              TempHitToAllHits.emplace_back(std::make_tuple(id_det, id_hit, id_track));
	      TempHitXYZ.emplace_back(TVector3(meanX,meanY,0.));

              for(int it1 = 0; it1 < 1000; ++it1)
                {
                  // double tempPhi1 = gRandom->Uniform(PhiXY1-dPhiXY1,PhiXY1+dPhiXY1);
                  // double tempPhi2 = gRandom->Uniform(PhiXY2-dPhiXY2,PhiXY2+dPhiXY2);

                  double tempPhi1 = gRandom->Uniform(PhiXY1, PhiXY2);
                  double tempR1   = gRandom->Uniform(RXY1, RXY2);

                  double tempPhiMin = TMath::Min(tempPhi1 - dPhiXY1, tempPhi1 + dPhiXY1);
                  double tempPhiMax = TMath::Min(tempPhi1 - dPhiXY1, tempPhi1 + dPhiXY1);
                  for(int it = 0; it < 1000; ++it)
                    h_RPhi->Fill(gRandom->Uniform(tempR1 - dRXY1, tempR1 + dRXY1),
                                 gRandom->Uniform(tempPhiMin, tempPhiMax), id_track);
                }
            }

          if(id_det == G4Sol::PSCE)
            {

              auto SimHit = it_ListHitsSim->second[id_det][0];
              h_XYSim->Fill(SimHit.hitX, SimHit.hitY, id_track);
              h_RPhiSim->Fill(TMath::Sqrt(TMath::Sq(SimHit.hitX) + TMath::Sq(SimHit.hitY)),
                              TMath::ATan2(SimHit.hitY, SimHit.hitX) * TMath::RadToDeg(), id_track);
#ifdef DEBUG_RIEMANNFINDER
              att._logger->debug(" layer ID : {} / R : {} / Phi : {} / Theta : {} | XYZ {} {} {} | ID {} / {}", id_det,
                                 TMath::Sqrt(TMath::Sq(SimHit.hitX) + TMath::Sq(SimHit.hitY)),
                                 TMath::ATan2(SimHit.hitY, SimHit.hitX) * TMath::RadToDeg(), -1., SimHit.hitX, SimHit.hitY, SimHit.hitZ, SimHit.layerID, it_ListHitsSim->second[id_det].size());
#endif
              SimHits.emplace_back(std::array<double, 3>{SimHit.hitX, SimHit.hitY, SimHit.hitZ});

              genfit::PlanarMeasurement* currentHit =
                  dynamic_cast<genfit::PlanarMeasurement*>(RecoEvent.ListHits[id_det][id_hit].get());
              auto dummyState   = genfit::StateOnPlane();
              auto currentPlane = currentHit->constructPlane(genfit::StateOnPlane());
              auto HitCoord     = currentHit->getRawHitCoords();
              auto HitCov       = currentHit->getRawHitCov();

              TVector2 tempLocal;
              if(HitCoord.GetNrows() == 1)
                tempLocal.Set(HitCoord[0], 0.);
              else
                tempLocal.Set(HitCoord[0], HitCoord[1]);

              auto currentLabPos = currentPlane->toLab(tempLocal);
              double meanX       = currentLabPos.X();
              double meanY       = currentLabPos.Y();
              double R           = TMath::Sq(meanX) + TMath::Sq(meanY);
              h_XY->Fill(meanX, meanY, id_track );
              h_XYConformal->Fill(meanX / R, meanY / R, id_track );

              double Phi1       = TMath::ATan2(meanY, meanX) * TMath::RadToDeg();
              double dPhiXY1    = TMath::ASin(TMath::Sqrt(HitCov[0][0] / R));
              double dRXY1      = 0.5;
              double tempPhiMin = TMath::Min(Phi1 - dPhiXY1, Phi1 + dPhiXY1);
              double tempPhiMax = TMath::Min(Phi1 - dPhiXY1, Phi1 + dPhiXY1);

              for(int it1 = 0; it1 < 1000; ++it1)
                h_RPhi->Fill(gRandom->Uniform(R - dRXY1, R + dRXY1), gRandom->Uniform(tempPhiMin, tempPhiMax),
                             id_track);

              double Xp = TMath::Sqrt(R) * TMath::Cos(Phi1 * TMath::DegToRad()) / (1. + R);
              double Yp = TMath::Sqrt(R) * TMath::Sin(Phi1 * TMath::DegToRad()) / (1. + R);
              double Zp = R / (1. + R);

              h_XYRiemann->Fill(TMath::Sqrt(Xp * Xp + Yp * Yp), (TMath::ATan2(Yp, Xp)+TMath::Pi()) * TMath::RadToDeg(),
                                id_track);

              TempHits.emplace_back(
				    RPhiHit(TMath::Sqrt(Xp * Xp + Yp * Yp), (TMath::ATan2(Yp, Xp)+TMath::Pi()) * TMath::RadToDeg(), Zp, nCount));
              ++nCount;
              TempHitToAllHits.emplace_back(std::make_tuple(id_det, id_hit, id_track));
	      TempHitXYZ.emplace_back(TVector3(meanX,meanY,currentLabPos.Z()));
	    }

	  if(id_det == LastFrontWall)
            {
              auto SimHit = it_ListHitsSim->second[id_det][0];
#ifdef DEBUG_RIEMANNFINDER
              att._logger->debug(" layer ID : {} / R : {} / Phi : {} / Theta : {} / Pt {} / Phi {} | XYZ {} {} {} | ID {} / {}", id_det,
                                 TMath::Sqrt(TMath::Sq(SimHit.hitX) + TMath::Sq(SimHit.hitY)),
                                 TMath::ATan2(SimHit.hitY, SimHit.hitX) * TMath::RadToDeg(), -1., TMath::Hypot(SimHit.momX,SimHit.momY), TMath::ATan2(SimHit.momY, SimHit.momX) * TMath::RadToDeg(), SimHit.hitX, SimHit.hitY, SimHit.hitZ,  SimHit.layerID, it_ListHitsSim->second[id_det].size());
#endif
              SimHits.emplace_back(std::array<double, 3>{SimHit.hitX, SimHit.hitY, SimHit.hitZ});
              TempHitPSBE.emplace_back(std::make_tuple(id_det, id_hit, id_track));
	    }
	  if(id_det == G4Sol::PSFE)
            {
              auto SimHit = it_ListHitsSim->second[id_det][0];
#ifdef DEBUG_RIEMANNFINDER
              att._logger->debug(" PSCE layer ID : {} / Pt : {} / Phi : {} / Theta : {}", id_det,
                                 TMath::Sqrt(TMath::Sq(SimHit.hitX) + TMath::Sq(SimHit.hitY)),
                                 TMath::ATan2(SimHit.hitY, SimHit.hitX) * TMath::RadToDeg(), -1.);
#endif
              //SimHits.emplace_back(std::array<double, 3>{SimHit.hitX, SimHit.hitY, SimHit.hitZ});
              //TempHitPSBE.emplace_back(std::make_tuple(id_det, id_hit, id_track));
	    }

        }

      // if(SimHits.size()>2)
      // 	{
      // 	  tricktrack::Matrix3xNd riemannHits = tricktrack::Matrix3xNd::Random(3,SimHits.size());
      // 	  for(size_t idH = 0; idH< SimHits.size();++idH)
      // 	    riemannHits.col(idH) << SimHits[idH][0], SimHits[idH][1], SimHits[idH][2];

      // 	  tricktrack::Matrix3Nd hits_cov = tricktrack::Matrix3Nd::Identity(3*SimHits.size(),3*SimHits.size());
      // 	  auto h = tricktrack::Helix_fit(riemannHits, hits_cov, 1./(0.30286*1000.), true, true);

      // 	  //att._logger->debug("Fit: para [phi {}, Tip {}, pt {}, cotanTheta {}, Zip {}] , chi2_line {},
      // chi2_circle {}, q {} ", h.par(0), h.par(1), h.par(2), h.par(3), h.par(4), h.chi2_line, h.chi2_circle, h.q);
      // 	  att._logger->debug("Sim Hit : {} {}", it_Info->second[FirstId].pdg,
      // TMath::Sqrt(TMath::Sq(it_Info->second[FirstId].momX) + TMath::Sq(it_Info->second[FirstId].momY)));
      // 	}
      // else
      // 	att._logger->debug("Fit : not enough hit");
    }
#ifdef DEBUG_RIEMANNFINDER
  for(size_t i_h = 0; i_h < TempHits.size(); ++i_h)
    att._logger->debug("RPhi hit# {} : {}, {}, {}, {} | LayerID {} hit {} track {}", i_h, TempHits[i_h][0], TempHits[i_h][1], TempHits[i_h][2],
                       TempHits[i_h].getId(), std::get<0>(TempHitToAllHits[TempHits[i_h].getId()])-G4Sol::MG01,
                       std::get<1>(TempHitToAllHits[TempHits[i_h].getId()]),std::get<2>(TempHitToAllHits[TempHits[i_h].getId()]) );

  att._logger->debug("KDPoint done ");
#endif

  KDtree.build(TempHits);
  // att._logger->debug("KDTree build done: {}", KDtree.size());
  // auto TreeIds = KDtree.getTheIds();
  // std::stringstream ss_tree;
  // ss_tree<<" ids of KDtree :";
  // for(size_t id_ids = 0; id_ids< TreeIds.size(); ++id_ids)
  //   ss_tree<<" #"<<id_ids<<": "<<TreeIds[id_ids];
  // att._logger->debug(ss_tree.str());

}

void TRiemannFinder::BuildTrackCand(const FullRecoEvent& RecoEvent, tricktrack::FKDTree<RPhiHit, double, 3>& KDtree, const std::vector<RPhiHit>& TempHits, std::vector<RTrack>& newTracksCand)
{

  std::vector<unsigned int> reverseTempHitIndex(TempHits.size());

#ifdef DEBUG_RIEMANNFINDER
  att._logger->debug("ReverseId");
#endif
  for(size_t idtempHit = 0;idtempHit<reverseTempHitIndex.size();++idtempHit)
    {
      reverseTempHitIndex[TempHits[idtempHit].getId()]=idtempHit;
#ifdef DEBUG_RIEMANNFINDER
      att._logger->debug("{} -> {}",TempHits[idtempHit].getId(),idtempHit);
#endif
    }

#ifdef DEBUG_RIEMANNFINDER
  int NtrackC = 0;
#endif

  std::vector<std::set<int> > TrackCandidates;

  for(int iR = 0; iR < 10; ++iR)
    {
      std::vector<unsigned int> initialHits;

      double minR = iR * 0.02;
      double maxR = (iR + 1) * 0.02;

      RPhiHit minHit(minR, 0., 0.8), maxHit(maxR, 360.,1.2);

      //initialHits.resize(0);
      //minHit.setDimension(0, minR);
      //maxHit.setDimension(0, maxR);
      KDtree.search(minHit, maxHit, initialHits);

#ifdef DEBUG_RIEMANNFINDER
      for(size_t i_h = 0; i_h < TempHits.size(); ++i_h)
	att._logger->debug("RPhi hit# {} : {}, {}, {}, {} | LayerID {} hit {} track {}", i_h, TempHits[i_h][0], TempHits[i_h][1], TempHits[i_h][2],
			   TempHits[i_h].getId(), std::get<0>(TempHitToAllHits[TempHits[i_h].getId()])-G4Sol::MG01,
			   std::get<1>(TempHitToAllHits[TempHits[i_h].getId()]),std::get<2>(TempHitToAllHits[TempHits[i_h].getId()]) );
#endif

      if(initialHits.size() > 0)
        {
#ifdef DEBUG_RIEMANNFINDER
	  att._logger->debug("scan R: {} {} | found {}",minR,maxR, initialHits.size());
#endif
          for(auto TempInitHit : initialHits)
            {
	      auto initHit = reverseTempHitIndex[TempInitHit];
              double iniR   = TempHits[initHit][0];
              double iniPhi = TempHits[initHit][1];
	      double iniZ   = TempHits[initHit][2];
#ifdef DEBUG_RIEMANNFINDER
              att._logger->debug("Found IniHit : {} [{} {} {}] {} | {} {} {}", initHit, iniR, iniPhi, iniZ, TempHits[initHit].getId(),
				 std::get<0>(TempHitToAllHits[TempHits[initHit].getId()]),std::get<1>(TempHitToAllHits[TempHits[initHit].getId()]),std::get<2>(TempHitToAllHits[TempHits[initHit].getId()])
				 );
#endif
              RPhiHit minHit2(iniR*0.95, iniPhi - 3.9, iniZ*0.9), maxHit2(iniR+(0.02+0.01*iR), iniPhi+ 3.9, iniZ*1.1);

              // minHit2.setDimension(0, iniR * 0.95);
              // minHit2.setDimension(1, iniPhi - 3.5);
              // // minHit2.setDimension(2,iniZ*0.90);

              // maxHit2.setDimension(0, iniR+0.02);
              // maxHit2.setDimension(1, iniPhi + 3.5);
              // maxHit2.setDimension(2,iniZ*1.10);
              std::vector<unsigned int> trackC;
              KDtree.search(minHit2, maxHit2, trackC);

#ifdef DEBUG_RIEMANNFINDER
              att._logger->debug("track Candidate: {} | found adjacent {} | cuts : Min [{} {} {}] Max [{} {} {}]", NtrackC++, trackC.size(),minHit2[0],minHit2[1],minHit[2], maxHit2[0], maxHit2[1],maxHit2[2]);
	      bool noFound=true;
	      for(auto Temp_idT : trackC)
                {
		  auto idT = reverseTempHitIndex[Temp_idT];

		  att._logger->debug(" {} / {} - {} {} {}", idT, TempHits[idT].getId(), std::get<0>(TempHitToAllHits[TempHits[idT].getId()]),std::get<1>(TempHitToAllHits[TempHits[idT].getId()]),std::get<2>(TempHitToAllHits[TempHits[idT].getId()]));
		  if(idT==initHit)
		    noFound=false;
		}
	      if(noFound == true)
		att._logger->debug("initHit no found in the proximity search !");
#endif
              if(trackC.size() < 2)
                continue;

              std::set<int> trackSet;
              for(auto iH : trackC)
                trackSet.insert(iH);

              bool IsInsert = false;
              for(auto& prevT : TrackCandidates)
                {
#ifdef DEBUG_RIEMANNFINDER
		  std::stringstream sst_pr;
		  sst_pr<<" --- ";
		  for(auto iS : prevT)
		    sst_pr<<iS<<" ";
		  att._logger->debug(sst_pr.str());
#endif
                  std::set<int> res;
                  std::set_intersection(prevT.begin(), prevT.end(), trackSet.begin(), trackSet.end(),
                                        std::inserter(res, res.begin()));
                  if(res.size() > 0)
                    {
#ifdef DEBUG_RIEMANNFINDER
		      std::stringstream sst_res;
		      sst_res<<" >> ";
		      for(auto iRset : res)
			sst_res<<iRset<<" ";
		      att._logger->debug(sst_res.str());
#endif
                      prevT.insert(trackSet.begin(), trackSet.end());
                      IsInsert = true;
                      break;
                    }
                }
              if(IsInsert == false)
                {
                  TrackCandidates.emplace_back(trackSet);
                }
            }
        }
    }

#ifdef DEBUG_RIEMANNFINDER
  att._logger->debug("-- Before merge --");
  int nTc = 0;
  for(auto trackC : TrackCandidates)
    {
      std::stringstream ssT;
      ssT << " track  " << nTc++ << " : ";
      for(auto t_hit : trackC)
        {
	  ssT << t_hit << " [";
	  auto [id_det1, id_hit1, id_track1] = TempHitToAllHits[t_hit];
	  ssT<<id_det1 - G4Sol::MG01<<", "<<id_hit1<<", "<<id_track1<<" ] ";
	}
      att._logger->debug(ssT.str());
    }
  att._logger->debug("Plots tracks :");
#endif

  auto f_extract_WireXY = [](const auto& ListHits, int id_det, int id_hit, bool DoCov=false, TMatrixD* cov=nullptr, double* dl = nullptr) {
    genfit::WireMeasurement* currentHit1 = dynamic_cast<genfit::WireMeasurement*>(ListHits[id_det][id_hit].get());

    auto HitCoord = currentHit1->getRawHitCoords();
    auto HitCov   = currentHit1->getRawHitCov();

    TVector3 wire_side1(HitCoord[0], HitCoord[1], HitCoord[2]);
    TVector3 wire_side2(HitCoord[3], HitCoord[4], HitCoord[5]);

    double meanX = 0.5 * (wire_side1.X() + wire_side2.X());
    double meanY = 0.5 * (wire_side1.Y() + wire_side2.Y());
    if(DoCov)
      {
	*dl = HitCoord[6];
	cov->ResizeTo(2,2);
	cov->Zero();
	(*cov)[0][0] = HitCov[6][6];
	(*cov)[1][1] = HitCov[6][6];
      }
    return std::make_tuple(meanX, meanY);
  };
  auto f_extract_PSBXY = [](const auto& ListHits, int id_det, int id_hit, bool DoCov=false, TMatrixD* cov = nullptr) {
    genfit::PlanarMeasurement* currentHit = dynamic_cast<genfit::PlanarMeasurement*>(ListHits[id_det][id_hit].get());
    auto dummyState                       = genfit::StateOnPlane();
    auto currentPlane                     = currentHit->constructPlane(genfit::StateOnPlane());
    auto HitCoord                         = currentHit->getRawHitCoords();
    auto HitCov                           = currentHit->getRawHitCov();

    TVector2 tempLocal;
    if(HitCoord.GetNrows() == 1)
      tempLocal.Set(HitCoord[0], 0.);
    else
      tempLocal.Set(HitCoord[0], HitCoord[1]);

    auto currentLabPos = currentPlane->toLab(tempLocal);
    double meanX       = currentLabPos.X();
    double meanY       = currentLabPos.Y();
    if(DoCov==true)
      {
	TMatrixDSym hitCovXY(2);
	hitCovXY.Zero();
	hitCovXY(0,0) = HitCov(0,0);
	hitCovXY(1,1) = 1./12.;
	TMatrixD rotMat(2,2);
	rotMat.Zero();
	tempLocal.Set(1.,0.);
	currentLabPos = currentPlane->toLab(tempLocal);
	rotMat[0][0] = currentLabPos.X();
	rotMat[0][1] = currentLabPos.Y();
	tempLocal.Set(0.,1.);
	currentLabPos = currentPlane->toLab(tempLocal);
	rotMat[1][0] = currentLabPos.X();
	rotMat[1][1] = currentLabPos.Y();

	TMatrixD newCov (rotMat,TMatrixD::kMult,hitCovXY);
	rotMat.T();
	newCov *= rotMat;
	cov->ResizeTo(3,3);
	cov->Zero();
	(*cov)[0][0] = newCov[0][0];
	(*cov)[1][0] = newCov[1][0];
	(*cov)[0][1] = newCov[0][1];
	(*cov)[1][1] = newCov[1][1];
	(*cov)[2][2] = HitCov(1,1);
      }
    return std::make_tuple(meanX, meanY);
  };

  int idTrack = 0;
  for(auto trackC : TrackCandidates)
    {
#ifdef DEBUG_RIEMANNFINDER
      att._logger->debug("track {}", idTrack);
#endif
      auto itR_hitT = trackC.begin();
      bool done    = false;
      while(done != true)
        // r(auto it_hitT = trackC.begin(), it_end = trackC.end(); it_hitT != it_end; ++it_hitT)
        {
          auto itR_hitT2 = itR_hitT;
          itR_hitT2++;
	  auto it_hitT = *itR_hitT;
	  auto it_hitT2 = *itR_hitT2;

          auto [id_det1, id_hit1, id_track1] = TempHitToAllHits[it_hitT];
          auto [id_det2, id_hit2, id_track2] = TempHitToAllHits[it_hitT2];

#ifdef DEBUG_RIEMANNFINDER
          att._logger->debug("  -- {} {} ; [{} {}] [{} {}]", it_hitT, it_hitT2, id_det1, id_hit1, id_det2, id_hit2);
#endif
          auto [meanX1, meanY1] = id_det1 == G4Sol::PSCE ? f_extract_PSBXY(RecoEvent.ListHits, id_det1, id_hit1)
                                                         : f_extract_WireXY(RecoEvent.ListHits, id_det1, id_hit1);
          auto [meanX2, meanY2] = id_det2 == G4Sol::PSCE ? f_extract_PSBXY(RecoEvent.ListHits, id_det2, id_hit2)
                                                         : f_extract_WireXY(RecoEvent.ListHits, id_det2, id_hit2);

          int minBin = h_XYTracks->GetXaxis()->FindBin(TMath::Min(meanX1, meanX2));
          int maxBin = h_XYTracks->GetXaxis()->FindBin(TMath::Max(meanX1, meanX2));
          for(int iBin = minBin; iBin < maxBin; ++iBin)
            {
              double tempX = h_XYTracks->GetXaxis()->GetBinCenter(iBin);
              double tempY = (meanY2 - meanY1) / (meanX2 - meanX1) * (tempX - meanX1) + meanY1;
              h_XYTracks->Fill(tempX, tempY, idTrack + 1);
            }

          itR_hitT2++;
          if(itR_hitT2 == trackC.end())
            done = true;
          itR_hitT++;
        }
      idTrack++;
    }

  TempCovXYZ.resize(TempHitToAllHits.size(), {0.,0.,0.,0.,0.,0.,0.,0.,0.});
  for(auto trackC : TrackCandidates)
    {
      if(trackC.size()<3)
	{
	  newTracksCand.emplace_back(trackC, 0.,0.,0., -1., 0.,0.,0);
	  continue;
	}
      int nWire = 0;
      int nPSB = 0;
      for(auto it_hit : trackC)
	if(auto [id_det1, id_hit1, id_track1] = TempHitToAllHits[it_hit] ; id_det1 != G4Sol::PSCE)
	  ++nWire;
	else
	  ++nPSB;
      //int totalHit = 2*nWire + nPSB;
      int totalHit = nWire + nPSB;
#ifdef DEBUG_RIEMANNFINDER
      att._logger->debug(" --> Circle fit : nHits {}",totalHit);
#endif
      tricktrack::Matrix3xNd riemannHits = tricktrack::Matrix3xNd::Random(3,totalHit);
      tricktrack::Matrix3Nd hits_cov = tricktrack::Matrix3Nd::Identity(3*totalHit,3*totalHit);
      tricktrack::Matrix3d tempCov = tricktrack::Matrix3d::Identity(3,3);

      int nhit = 0;
      double tempZ = 0.1;
      double minPhi = 380.;
      double maxPhi = -380.;

      for(auto it_hit : trackC)
	{
	  auto [id_det1, id_hit1, id_track1] = TempHitToAllHits[it_hit];
	  TMatrixD CovXY(3,3);
	  CovXY.Zero();
	  double dl = 0;
	  const int idMDC     = id_det1 - G4Sol::MG01;
	  double R_size       = dl_max(idMDC) * 0.5;
	  auto [meanX1, meanY1] = id_det1 == G4Sol::PSCE ? f_extract_PSBXY(RecoEvent.ListHits, id_det1, id_hit1,true,&CovXY)
	    : f_extract_WireXY(RecoEvent.ListHits, id_det1, id_hit1,true,&CovXY,&dl);

	  double tempR  = TMath::Sqrt(TMath::Sq(meanX1)+TMath::Sq(meanY1));
	  double tempPhi = TMath::ATan2(meanY1,meanX1);
	  double dPhi = TMath::ASin(R_size/tempR);
	  if(minPhi> tempPhi-dPhi)
	    minPhi = tempPhi-dPhi;
	  if(maxPhi<tempPhi+dPhi)
	    maxPhi = tempPhi+dPhi;

	  if(id_det1 != G4Sol::PSCE)
	    {
	      // double newX1 = tempR*TMath::Cos(tempPhi+dPhi);
	      // double newY1 = tempR*TMath::Sin(tempPhi+dPhi);

	      // double newX2 = tempR*TMath::Cos(tempPhi-dPhi);
	      // double newY2 = tempR*TMath::Sin(tempPhi-dPhi);

	      // riemannHits.col(nhit) << newX1, newY1, 0.;
	      // tempCov(0,0) = CovXY[0][0]; tempCov(0,1) = CovXY[0][1];
	      // tempCov(1,0) = CovXY[1][0]; tempCov(1,1) = CovXY[1][1];
	      // hits_cov.block<3,3>(3*nhit,3*nhit) = tempCov;
	      // ++nhit;
	      // riemannHits.col(nhit) << newX2, newY2, 0.;
	      // hits_cov.block<3,3>(3*nhit,3*nhit) = tempCov;
	      // ++nhit;

	      riemannHits.col(nhit) << meanX1, meanY1, tempZ;
	      tempCov(0,0) = R_size *0.5; tempCov(0,1) = 0.;
	      tempCov(1,0) = 0.; tempCov(1,1) = R_size*0.5;
	      hits_cov(nhit,nhit) = tempCov(0,0);
	      hits_cov(nhit+totalHit,nhit+totalHit) =  tempCov(1,1);
	      hits_cov(nhit,nhit+totalHit) = tempCov(0,1);
	      hits_cov(totalHit+nhit,nhit) = tempCov(0,1);

	      TempCovXYZ[it_hit] = {tempCov(0,0),tempCov(1,0),tempCov(2,0),tempCov(0,1),tempCov(1,1),tempCov(2,1),tempCov(0,2),tempCov(1,2),tempCov(2,2)};

	      ++nhit;
	      tempZ += 1.;
	    }
	  else
	    {
	      riemannHits.col(nhit) << meanX1, meanY1, tempZ;
	      tempCov(0,0) = CovXY[0][0]; tempCov(0,1) = CovXY[0][1];
	      tempCov(1,0) = CovXY[1][0]; tempCov(1,1) = CovXY[1][1];
	      hits_cov(nhit,nhit) = tempCov(0,0);
	      hits_cov(nhit+totalHit,nhit+totalHit) =  tempCov(1,1);
	      hits_cov(nhit,nhit+totalHit) = tempCov(0,1);
	      hits_cov(totalHit+nhit,nhit) = tempCov(0,1);

	      TempCovXYZ[it_hit] = {tempCov(0,0),tempCov(1,0),tempCov(2,0),tempCov(0,1),tempCov(1,1),tempCov(2,1),tempCov(0,2),tempCov(1,2),tempCov(2,2)};

	      ++nhit;
	      tempZ += 1.;
	    }
	}

      u_int nC = riemannHits.cols();
      tricktrack::VectorNd rad = (riemannHits.block(0, 0, 2, nC).colwise().norm());
#ifdef DEBUG_RIEMANNFINDER
      std::stringstream radOut;
      radOut<<rad;
      att._logger->debug("radial matrix : ");
      att._logger->debug(radOut.str());
#endif
      // Fast_fit gives back (X0, Y0, R, theta) w/o errors, using only 3 points.
      const tricktrack::Vector4d fast_fit = tricktrack::Fast_fit(riemannHits);

#ifdef DEBUG_RIEMANNFINDER
      std::stringstream radOut1;
      radOut1<<fast_fit;
      att._logger->debug("fast fit :");
      att._logger->debug(radOut1.str());
#endif
      tricktrack::circle_fit circle = tricktrack::Circle_fit(riemannHits.block(0, 0, 2, nC), hits_cov.block(0, 0, 2 * nC, 2 * nC),
							     fast_fit, rad, true, true);

      newTracksCand.emplace_back(trackC, circle.par(0), circle.par(1), circle.par(2), circle.chi2, minPhi, maxPhi, circle.q);
      //auto h = tricktrack::Helix_fit(riemannHits, hits_cov, 1./(0.30286*1000.), true, true);
#ifdef DEBUG_RIEMANNFINDER
      att._logger->debug("Fit: para [X0 {}, Y0 {}, R0 {}] , chi2 {}, q {}, minPhi {}, maxPhi {} ", circle.par(0), circle.par(1), circle.par(2), circle.chi2, circle.q, minPhi, maxPhi);
#endif
    }

}

void TRiemannFinder::AddStereoWire(const FullRecoEvent& RecoEvent, std::vector<RTrack>& newTracksCand)
{

  auto f_propagWire = [](const TVector3& w1, const TVector3& w2, TVector3& Pos) {
    TVector3 Dir = w2 - w1;
    Dir *= 1. / (w2.Z() - w1.Z());
    double posZ = Pos.Z();
    Pos = Dir * (posZ - w1.Z()) + w1;
  };

#ifdef DEBUG_RIEMANNFINDER
  att._logger->debug("redo layer step : {}",redo_layer.size());
#endif

  for(auto layerRedo : redo_layer)
    {
#ifdef DEBUG_RIEMANNFINDER
      att._logger->debug("#############");
      att._logger->debug("Redo Layer {}", layerRedo - G4Sol::MG01);
#endif
      for(size_t id_hit = 0; id_hit < RecoEvent.ListHits[layerRedo].size(); ++id_hit)
        {
	  const int idMDCL     = layerRedo - G4Sol::MG01;
	  double R_sizeL       = dl_max(idMDCL) * 0.5;

	  genfit::WireMeasurement* currentHit1 = dynamic_cast<genfit::WireMeasurement*>(RecoEvent.ListHits[layerRedo][id_hit].get());

	  auto HitCoord = currentHit1->getRawHitCoords();

	  TVector3 wire_side1(HitCoord[0], HitCoord[1], HitCoord[2]);
	  TVector3 wire_side2(HitCoord[3], HitCoord[4], HitCoord[5]);

	  double Phi_S1 = wire_side1.Phi()+TMath::Pi();
	  double dPhi_S1 = TMath::ASin(R_sizeL/wire_side1.Perp());

	  double Phi_S2 = wire_side2.Phi()+TMath::Pi();
	  double dPhi_S2 = TMath::ASin(R_sizeL/wire_side2.Perp());

	  double minPhi = TMath::Min(Phi_S1-dPhi_S1,Phi_S2-dPhi_S2);
	  double maxPhi = TMath::Max(Phi_S1+dPhi_S1,Phi_S2+dPhi_S2);

	  if(minPhi > maxPhi)
	    std::swap(minPhi,maxPhi);

	  double alphaX = (wire_side2.X()-wire_side1.X()) / (wire_side2.Z() - wire_side1.Z());
	  double alphaY = (wire_side2.Y()-wire_side1.Y()) / (wire_side2.Z() - wire_side1.Z());

	  bool isAugmented = false;
	  if(minPhi < 0.2*TMath::Pi() || maxPhi > 1.8*TMath::Pi())
	    {
	      minPhi += 0.5*TMath::Pi();
	      maxPhi += 0.5*TMath::Pi();

	      if(minPhi > TMath::TwoPi())
		minPhi -= TMath::TwoPi();
	      if(maxPhi > TMath::TwoPi())
		maxPhi -= TMath::TwoPi();
	      isAugmented = true;
	      if(minPhi > maxPhi)
		std::swap(minPhi,maxPhi);

	    }
#ifdef DEBUG_RIEMANNFINDER
	  att._logger->debug(" -- hit : {}, Phi1 {} Phi2 {}, dPhi1 {}, dPhi2 {}, minPhi {}, maxPhi {} | Augmented ? {}", id_hit, Phi_S1*TMath::RadToDeg(), Phi_S2*TMath::RadToDeg(), dPhi_S1*TMath::RadToDeg(), dPhi_S2*TMath::RadToDeg(), minPhi*TMath::RadToDeg(), maxPhi*TMath::RadToDeg(), isAugmented);
	  int nTc = 0;
#endif
	  for(auto& trackC : newTracksCand)
	    {
	      std::set<int> toAdd ;
#ifdef DEBUG_RIEMANNFINDER
	      std::stringstream ssT1;
#endif
	      double TrackMinPhi = trackC.minPhi+TMath::Pi();
	      double TrackMaxPhi = trackC.maxPhi+TMath::Pi();
	      if(isAugmented)
		{
		  TrackMinPhi += 0.5*TMath::Pi();
		  TrackMaxPhi += 0.5*TMath::Pi();
		  if(TrackMinPhi > TMath::TwoPi())
		    TrackMinPhi -= TMath::TwoPi();
		  if(TrackMaxPhi > TMath::TwoPi())
		    TrackMaxPhi -= TMath::TwoPi();
		  if(TrackMinPhi > TrackMaxPhi)
		    std::swap(TrackMinPhi,TrackMaxPhi);
		}
#ifdef DEBUG_RIEMANNFINDER
	      ssT1 << " track  " << nTc++ << " [ "<<TrackMinPhi*TMath::RadToDeg()<<", "<<TrackMaxPhi*TMath::RadToDeg()<<" ] : ";
	      for(auto t_hit : trackC.hits)
		ssT1 << t_hit << " [ "<<std::get<0>(TempHitToAllHits[t_hit]) - G4Sol::MG01<<"] ";
	      att._logger->debug(ssT1.str());
#endif
	      // for(auto it_hitT = trackC.hits.begin(), it_hitT2 = ++(trackC.hits.begin()), it_hitend = trackC.hits.end(); it_hitT2 != it_hitend;++it_hitT, ++it_hitT2)
	      // 	{

	      // 	  auto [id_det1, id_hit1] = TempHitToAllHits[*it_hitT];
	      // 	  auto [id_det2, id_hit2] = TempHitToAllHits[*it_hitT2];

	      // 	  if(id_det1 > layerRedo || layerRedo > id_det2)
	      // 	    continue;

	      // 	  att._logger->debug(" In Track : {} {}", id_det1- G4Sol::MG01, id_det2- G4Sol::MG01);
	      // 	  //att._logger->debug("  -- {} {} ; [{} {}] [{} {}]", *it_hitT, *it_hitT2, id_det1, id_hit1, id_det2, id_hit2);

	      // 	  auto [meanX1, meanY1] = id_det1 == G4Sol::PSCE ? f_extract_PSBXY(RecoEvent.ListHits, id_det1, id_hit1)
	      // 	    : f_extract_WireXY(RecoEvent.ListHits, id_det1, id_hit1);
	      // 	  auto [meanX2, meanY2] = id_det2 == G4Sol::PSCE ? f_extract_PSBXY(RecoEvent.ListHits, id_det2, id_hit2)
	      // 	    : f_extract_WireXY(RecoEvent.ListHits, id_det2, id_hit2);

	      // 	  const int idMDC     = id_det2 - G4Sol::MG01;
	      // 	  double R_size       = dl_max(idMDC) * 0.5;

	      // 	  double Phi1 = TMath::ATan2(meanY1,meanX1);
	      // 	  double R1 = TMath::Sqrt(TMath::Sq(meanX1)+TMath::Sq(meanY1));

	      // 	  double Phi2 = TMath::ATan2(meanY2,meanX2);
	      // 	  double R2 = TMath::Sqrt(TMath::Sq(meanX2)+TMath::Sq(meanY2));
	      // 	  double dPhi2 = TMath::ASin(R_size / R2);

	      // 	  att._logger->debug("   -- Phi Before {}, Phi After {}, check {} {}", Phi1, Phi2, Phi2+dPhi2, Phi2-dPhi2);
	      if(trackC.chi2<0)
		continue;
	      // 	    {
	      // 	      att._logger->debug("In Phi checker chi2<0 {} | min {}, max {}",trackC.chi2,minPhi, maxPhi);
	      // 	      if(Phi2 + dPhi2 < minPhi || Phi2 - dPhi2 > maxPhi)
	      // 		continue;

	      // 	      double dir_t = (meanY2 - meanY1) / (meanX2 - meanX1);
	      // 	      double posZwire = wire_side1.Z() - (dir_t*wire_side1.X() - wire_side1.Y()) / (dir_t*alphaX - alphaY) + (dir_t*meanX1 - meanY1) / (dir_t*alphaX - alphaY);

	      // 	      double comPhi = Phi1 + TMath::ATan2(R2*TMath::Sin(Phi2-Phi1), R1 + R2*TMath::Cos(Phi2-Phi1));
	      // 	      double comZ = (TMath::Tan(comPhi)*alphaX - alphaY)*wire_side1.Z() + wire_side1.Y() - TMath::Tan(comPhi)*wire_side1.X();

	      // 	      att._logger->debug(" Hit layer belong to track ! {} {}", dir_t, posZwire);

	      // 	      TVector3 comPos (0.,0.,posZwire);
	      // 	      f_propagWire(wire_side1,wire_side2,comPos);

	      // 	      TempHitToAllHits.emplace_back(std::make_tuple(layerRedo,id_hit));
	      // 	      TempHitXYZ.emplace_back(comPos);

	      // 	      size_t newId = TempHitToAllHits.size()-1;
	      // 	      toAdd.insert(newId);
	      // 	    }
	      // 	  else
	      //     {

	      if(TrackMaxPhi < minPhi || TrackMinPhi > maxPhi)
		continue;

	      double withCir_a = TMath::Sq(alphaX) + TMath::Sq(alphaY);
	      double withCir_b = alphaX*(wire_side1.X()-trackC.x0)+alphaY*(wire_side1.Y()-trackC.y0);
	      double withCir_c = TMath::Sq(wire_side1.X()-trackC.x0) + TMath::Sq(wire_side1.Y()-trackC.y0) - TMath::Sq(trackC.r0);

	      double CirZ1 = -TMath::Sqrt(TMath::Sq(withCir_b) - withCir_a*withCir_c) - withCir_b;
	      CirZ1 /= withCir_a;
	      double CirZ2 = TMath::Sqrt(TMath::Sq(withCir_b) - withCir_a*withCir_c) - withCir_b;
	      CirZ2 /= withCir_a;

	      double CirZ = TMath::Abs(CirZ1) < TMath::Abs(CirZ2) ? CirZ1 : CirZ2;

	      double diffZ = wire_side2.Z()-wire_side1.Z() ;

#ifdef DEBUG_RIEMANNFINDER
	      att._logger->debug("In Cirle checker chi2>0 {} | CirZ1 {} CirZ2 {} | DiffZ {}",trackC.chi2, CirZ1, CirZ2, diffZ);
#endif
	      // if(CirZ2 < 0.)
	      // 	continue;
	      // if(CirZ1 > diffZ)
	      // 	continue;
	      // if(CirZ1 < 0. && CirZ2 > diffZ)
	      // 	continue;

	      //double CirZ = 0.;
	      // if(CirZ1 < 0. && CirZ2 <= diffZ)
	      // 	CirZ = CirZ2;
	      // else if(CirZ1 <= diffZ && CirZ2 > diffZ)
	      // 	CirZ = CirZ1;
	      // else
	      // 	att._logger->debug("E> Circle / StereoWire : two solutions valid ! {} {}",CirZ1,CirZ2);

	      if(CirZ < 0. || CirZ > diffZ)
		continue;

	      CirZ += wire_side1.Z();
	      TVector3 comPos (0.,0.,CirZ);
	      f_propagWire(wire_side1,wire_side2,comPos);
#ifdef DEBUG_RIEMANNFINDER
	      att._logger->debug("Circle Extrap : {} {} {} | {} {}",comPos.X(),comPos.Y(),comPos.Z(), TMath::Sq(trackC.x0-comPos.X())+TMath::Sq(trackC.y0-comPos.Y()), TMath::Sq(trackC.r0));

	      TVector3 comPosW1 = wire_side1 - comPos;
	      TVector3 DirWire = wire_side2 - wire_side1;
	      double WireLength = DirWire.Mag();
	      DirWire = DirWire.Unit();
	      TVector3 PerpComWire = comPosW1 - (comPosW1.Dot(DirWire))*DirWire;
	      att._logger->debug("Check : dist on wire {} | WireL {} | Prep dist {}",comPosW1.Dot(DirWire),WireLength,PerpComWire.Mag());
#endif
	      TempHitToAllHits.emplace_back(std::make_tuple(layerRedo,id_hit,redo_layerToTrackID[layerRedo][id_hit]));
	      TempHitXYZ.emplace_back(comPos);

	      size_t newId = TempHitToAllHits.size()-1;
	      toAdd.insert(newId);
		      //}
	      //  }

#ifdef DEBUG_RIEMANNFINDER
	      att._logger->debug("update redo : {}",toAdd.size());
#endif
	      trackC.hits.insert(toAdd.begin(), toAdd.end());
	    }
        }
    }
}

void TRiemannFinder::AddEndCap(const FullRecoEvent& RecoEvent, std::vector<RTrack>& newTracksCand)
{
  for(auto trackC : newTracksCand)
    {
      bool hasPSCE = false;
      for(auto hitT : trackC.hits)
	if(auto [id_det1, id_hit1, id_track1] = TempHitToAllHits[hitT]; id_det1 == G4Sol::PSCE)
	  hasPSCE = true;

      if(hasPSCE)
	continue;

      for(size_t id_LastWall = 0; id_LastWall<TempHitPSBE.size();++id_LastWall)
	{
	  auto [id_det, id_hit, id_track] = TempHitPSBE[id_LastWall];

	  genfit::PlanarMeasurement* currentHit = dynamic_cast<genfit::PlanarMeasurement*>(RecoEvent.ListHits[id_det][id_hit].get());
	  auto dummyState                       = genfit::StateOnPlane();
	  auto currentPlane                     = currentHit->constructPlane(genfit::StateOnPlane());
	  auto HitCoord                         = currentHit->getRawHitCoords();
	  //auto HitCov                           = currentHit->getRawHitCov();

	  auto v = currentPlane->getV();
	  double PhiBar = v.Phi();

	  double minPhiBar = PhiBar - 3.5*TMath::DegToRad();
	  double maxPhiBar = PhiBar + 3.5*TMath::DegToRad();

	  if(trackC.minPhi>maxPhiBar || trackC.maxPhi < minPhiBar)
	    continue;

	  auto O = currentPlane->getO();

	  // cross point in PSBE
	  double a = v.Y() / v.X();
	  double b = O.Y() - a*O.X();
	  double Xhit1 = ((b-trackC.y0)*a - trackC.x0 - TMath::Sqrt(TMath::Sq((b-trackC.y0)*a - trackC.x0) - (1.+a*a)*(TMath::Sq(trackC.x0)+TMath::Sq(b-trackC.y0)-TMath::Sq(trackC.r0))))/(1.+a*a);
	  double Xhit2 = ((b-trackC.y0)*a - trackC.x0 + TMath::Sqrt(TMath::Sq((b-trackC.y0)*a - trackC.x0) - (1.+a*a)*(TMath::Sq(trackC.x0)+TMath::Sq(b-trackC.y0)-TMath::Sq(trackC.r0))))/(1.+a*a);

	  double Yhit1 = a*Xhit1+b;
	  double Yhit2 = a*Xhit1+b;

	  double R1 = TMath::Hypot(Xhit1,Yhit1);
	  double Phi1 = TMath::ATan2(Yhit1,Xhit1);
	  double R2 = TMath::Hypot(Xhit2,Yhit2);
	  double Phi2 = TMath::ATan2(Yhit2,Xhit2);

	  auto InPSWall = [&minPhiBar,&maxPhiBar](double R, double Phi)->bool {
	    if(R < 6. || R> 22.)
	      return false;

	    if(Phi>maxPhiBar || Phi<minPhiBar)
	      return false;

	    return true;
	  };

	  bool isHit1 = InPSWall(R1,Phi1);
	  bool isHit2 = InPSWall(R2,Phi2);

	  if((isHit1 == false && isHit2 == false) || (isHit1 == true && isHit2 == true))
	    {
	      att._logger->debug("E> CheckPSBE failed ! {} {} / {} {} =! r 6. 22. phi {} {}",R1,Phi1,R2,Phi2,minPhiBar,maxPhiBar);
	      continue;
	    }
	  if(isHit1)
	    TempHitXYZ.emplace_back(Xhit1,Yhit1,O.Z());
	  if(isHit2)
	    TempHitXYZ.emplace_back(Xhit2,Yhit2,O.Z());

	  TempHitToAllHits.emplace_back(TempHitPSBE[id_LastWall]);
	  trackC.hits.insert(TempHitToAllHits.size()-1);
	}
    }
}

int TRiemannFinder::FinderTrack(FullRecoEvent& RecoEvent)
{

  tricktrack::FKDTree<RPhiHit, double, 3> KDtree;
  std::vector<RPhiHit> TempHits;

  redo_layer.clear();
  redo_layerToTrackID.clear();
  TempHitToAllHits.clear();
  TempHitXYZ.clear();
  TempCovXYZ.clear();
  TempHitPSBE.clear();

  BuildKDTree(RecoEvent, KDtree, TempHits);

#ifdef DEBUG_RIEMANNFINDER
  for(size_t i_h = 0; i_h < TempHits.size(); ++i_h)
    att._logger->debug("RPhi hit# {} : {}, {}, {}, {} | LayerID {} hit {} track {}", i_h, TempHits[i_h][0], TempHits[i_h][1], TempHits[i_h][2],
                       TempHits[i_h].getId(), std::get<0>(TempHitToAllHits[TempHits[i_h].getId()])-G4Sol::MG01,
                       std::get<1>(TempHitToAllHits[TempHits[i_h].getId()]),std::get<2>(TempHitToAllHits[TempHits[i_h].getId()]) );
#endif

  std::vector<RTrack> newTracksCand;
  BuildTrackCand(RecoEvent, KDtree, TempHits, newTracksCand);

  AddStereoWire(RecoEvent, newTracksCand);

  AddEndCap(RecoEvent, newTracksCand);

  int idTrack = 0;
  for(auto trackC : newTracksCand)
    {
#ifdef DEBUG_RIEMANNFINDER
      att._logger->debug("track {}", idTrack);
      std::stringstream ssT1;
#endif

      for(auto itR_hitT = trackC.hits.begin(), it_end = trackC.hits.end(); itR_hitT != it_end; ++itR_hitT)
        {
	  auto it_hitT = *itR_hitT;
	  auto vecXYZ = TempHitXYZ[it_hitT];

#ifdef DEBUG_RIEMANNFINDER
          auto [id_det1, id_hit1, id_track1] = TempHitToAllHits[it_hitT];
	  ssT1<<" "<<id_det1 - G4Sol::MG01<<" ["<<id_hit1<<", "<<id_track1<<"]";
#endif

	  h_XYReco->Fill(vecXYZ.X(), vecXYZ.Y(), idTrack + 1);
        }

#ifdef DEBUG_RIEMANNFINDER
      att._logger->debug(ssT1.str());
      ++idTrack;
#endif
    }




  
  // std::set<std::set<int> > TrackCandidates2;
  // for(auto t1 = TrackCandidates.begin(), t1_end = TrackCandidates.end(); t1 != t1_end; ++t1)
  //   {
  //     auto t1c = t1;
  //     for(auto t2 = t1c, t2_end = TrackCandidates.end();t2 != t2_end; ++t2)
  // 	{
  // 	  for(auto idH : *t1)
  // 	    if(auto it_ret = t2->find(idH); it_ret != t2->end())
  // 	      {
  // 		std::set<int> t1n = *t1;
  // 		t1n.insert(t2->begin(), t2->end());
  // 		TrackCandidates2.insert(t1n);
  // 		TrackCandidates.remove(t2);
  // 		break;
  // 	      }

  // 	  TrackCandidates2.insert(*t1);
  // 	}
  //   }

  // att._logger->debug("----");
  // nTc = 0;
  // for(auto trackC : TrackCandidates2)
  //   {
  //     std::stringstream ssT;
  //     ssT<<" track  "<<nTc++<<" : ";
  //     for(auto t_hit : trackC)
  // 	ssT<<t_hit<<" ";
  //     att._logger->debug(ssT.str());
  //   }

  t_finder->Fill();

  h_XYSim->Reset();
  h_XY->Reset();
  h_XYReco->Reset();
  h_RPhiSim->Reset();
  h_RPhi->Reset();
  h_XYConformal->Reset();
  h_XYRiemann->Reset();
  h_XYTracks->Reset();
  return 0;
}
