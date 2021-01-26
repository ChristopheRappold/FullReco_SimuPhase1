#include "TFlatMCOutputML.h"

#include "FullRecoEvent.hh"
#include "ReturnRes.hh"

#include <tuple>

//#define DEBUG_FLATTEN

TFlatMCOutputML::TFlatMCOutputML(const THyphiAttributes& attribut)
    : TDataProcessInterface("FlatMCOutputML"), att(attribut), namefileFlat(att.FlatML_namefile)
{
  f_flat = new TFile(namefileFlat, "RECREATE");
  t_flat = new TTree("MCDataML_Tree", "Flat data event tracks");

  t_flat->Branch("b_tx", &b_tx, "b_tx/F");
  t_flat->Branch("b_ty", &b_ty, "b_ty/F");
  t_flat->Branch("b_vx", &b_vx, "b_vx/F");
  t_flat->Branch("b_vy", &b_vy, "b_vy/F");
  t_flat->Branch("b_vz", &b_vz, "b_vz/F");
  t_flat->Branch("b_x", &b_x, "b_x/F");
  t_flat->Branch("b_y", &b_y, "b_y/F");
  t_flat->Branch("b_z", &b_z, "b_z/F");

  t_flat->Branch("a_tx", &a_tx, "a_tx/F");
  t_flat->Branch("a_ty", &a_ty, "a_ty/F");
  t_flat->Branch("a_vx", &a_vx, "a_vx/F");
  t_flat->Branch("a_vy", &a_vy, "a_vy/F");
  t_flat->Branch("a_vz", &a_vz, "a_vz/F");
  t_flat->Branch("a_x", &a_x, "a_x/F");
  t_flat->Branch("a_y", &a_y, "a_y/F");
  t_flat->Branch("a_z", &a_z, "a_z/F");

  t_flat->Branch("poq", &poq, "poq/F");
  t_flat->Branch("qop", &qop, "qop/F");
  t_flat->Branch("tof", &tof, "tof/F");
  t_flat->Branch("psce_psbe", &psce_psbe, "psce_psbe/I");
  t_flat->Branch("psb_z", &psb_z, "psb_z/F");
  t_flat->Branch("q", &q, "q/F");
  t_flat->Branch("pdg", &pdg, "pdg/I");

  t_flat->SetDirectory(f_flat);
  // t_flat->AutoSave();

  att._logger->info("FlatMC : tree out set ");
}

TFlatMCOutputML::~TFlatMCOutputML()
{
  f_flat->cd();
  t_flat->Write();
  // f_flat->Write();
  f_flat->Close();
  if(f_flat != nullptr)
    {
      f_flat->Delete();
      f_flat = nullptr;
    }
}

void TFlatMCOutputML::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }

ReturnRes::InfoM TFlatMCOutputML::operator()(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
{
  int result_finder = Exec(RecoEvent, OutTree);

  return SoftExit(result_finder);
}

int TFlatMCOutputML::Exec(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree) { return FlattenOut(RecoEvent); }

ReturnRes::InfoM TFlatMCOutputML::SoftExit(int result_full) { return ReturnRes::Fine; }

void TFlatMCOutputML::SelectHists()
{
  // LocalHisto.h_RZStats      = AnaHisto->CloneAndRegister(AnaHisto->h_RZStats);
}

int TFlatMCOutputML::FlattenOut(FullRecoEvent& RecoEvent)
{
  int ntrack = -1;
  for(auto it_trackInfo : RecoEvent.TrackInfo)
    {
      ntrack++;
#ifdef DEBUG_FLATTEN
      att._logger->debug("track {}", ntrack);
#endif

      const int id_track  = it_trackInfo.first;
      auto it_ListHits    = RecoEvent.TrackDAF.find(id_track);
      auto it_ListHitsSim = RecoEvent.TrackDAFSim.find(id_track);

      std::tuple<int, int> id_before_mag = std::make_tuple(-1, -1);
      std::tuple<int, int> id_after_mag  = std::make_tuple(-1, -1);
      std::tuple<int, int> id_psb        = std::make_tuple(-1, -1);

      for(int id_det = 0; id_det < it_ListHits->second.size(); ++id_det)
        {
          int id_hit = it_ListHits->second[id_det];

          if(id_hit < 0)
            continue;

#ifdef DEBUG_FLATTEN
          att._logger->debug("id_det {} id_hit {}", id_det, id_hit);
#endif
          if(id_det >= G4Sol::MG01 && id_det <= G4Sol::MG17)
            {
#ifdef DEBUG_FLATTEN
              att._logger->debug("MDC: {}", id_det);
#endif
              if(id_det > std::get<0>(id_after_mag))
                {
                  // auto temp_hit = RecoEvent.ListHits[id_det][id_hit].get();
                  // id_after_mag = std::make_tuple(id_det,temp_hit->getHitId());

                  id_after_mag = std::make_tuple(id_det, id_hit);
#ifdef DEBUG_FLATTEN
                  att._logger->debug("update! {}", id_det);
#endif
                }
            }

          if(id_det >= G4Sol::MiniFiberD1_x1 && id_det <= G4Sol::MiniFiberD1_v2)
            {
#ifdef DEBUG_FLATTEN
              att._logger->debug("MiniFiber: ");
#endif
              if(id_det > std::get<0>(id_before_mag))
                {
                  // auto temp_hit = RecoEvent.ListHits[id_det][id_hit].get();
                  // id_before_mag = std::make_tuple(id_det,temp_hit->getHitId());

                  id_before_mag = std::make_tuple(id_det, id_hit);
#ifdef DEBUG_FLATTEN
                  att._logger->debug("update! {}", id_det);
#endif
                }
            }
          if(id_det == G4Sol::PSCE)
            {
#ifdef DEBUG_FLATTEN
              att._logger->debug("PSCE: ");
#endif
              if(id_det > std::get<0>(id_psb))
                {
                  // auto temp_hit = RecoEvent.ListHits[id_det][id_hit].get();
                  // id_psb = std::make_tuple(id_det,temp_hit->getHitId());

                  id_psb = std::make_tuple(id_det, id_hit);
#ifdef DEBUG_FLATTEN
                  att._logger->debug("update! {}", id_det);
#endif
                }

              psce_psbe = 0;
            }
          G4Sol::SolDet LastFrontWall = att.Wasa_Side == 0 ? G4Sol::PSFE : G4Sol::PSBE;
          if(id_det == LastFrontWall)
            {
#ifdef DEBUG_FLATTEN
              att._logger->debug("LastWall: ");
#endif
              if(id_det > std::get<0>(id_psb))
                {
                  // auto temp_hit = RecoEvent.ListHits[id_det][id_hit].get();
                  // id_psb = std::make_tuple(id_det,temp_hit->getHitId());

                  id_psb = std::make_tuple(id_det, id_hit);
#ifdef DEBUG_FLATTEN
                  att._logger->debug("update! {}", id_det);
#endif
                }

              psce_psbe = 1;
            }
        }

#ifdef DEBUG_FLATTEN
      att._logger->debug("id last fiber : {} {} | id last MDC {} {} | id PSB {} {}", std::get<0>(id_before_mag),
                         std::get<1>(id_before_mag), std::get<0>(id_after_mag), std::get<1>(id_after_mag),
                         std::get<0>(id_psb), std::get<1>(id_psb));
#endif

      if(std::get<0>(id_psb) == -1 || std::get<0>(id_before_mag) == -1 || std::get<0>(id_after_mag) == -1)
        continue;

#ifdef DEBUG_FLATTEN
      att._logger->debug("to fill ");
#endif
      std::vector<double> mom;
      int temp_pdg = 0;

      auto it_hitBeforeSim = it_ListHitsSim->second[std::get<0>(id_before_mag)];
      // if(it_hitBeforeSim.layerID==std::get<1>(id_before_mag))
      //{
      b_tx = it_hitBeforeSim.momX / it_hitBeforeSim.momZ;
      b_ty = it_hitBeforeSim.momY / it_hitBeforeSim.momZ;

      double temp_mom = TMath::Sqrt(TMath::Sq(it_hitBeforeSim.momX) + TMath::Sq(it_hitBeforeSim.momY) +
                                    TMath::Sq(it_hitBeforeSim.momZ));

      b_vx = it_hitBeforeSim.momX / temp_mom;
      b_vy = it_hitBeforeSim.momY / temp_mom;
      b_vz = it_hitBeforeSim.momZ / temp_mom;

      b_x = it_hitBeforeSim.hitX;
      b_y = it_hitBeforeSim.hitY;
      b_z = it_hitBeforeSim.hitZ;

      mom.emplace_back(temp_mom);

      temp_pdg = it_hitBeforeSim.pdg;
      //}
      // else
      // 	att._logger->error("FlatML : Sim Before LayerID different from id_before_mag.id_hit {}
      // {}",it_hitBeforeSim.layerID, std::get<1>(id_before_mag));

      auto it_hitAfterSim = it_ListHitsSim->second[std::get<0>(id_after_mag)];
      // if(it_hitAfterSim.layerID==std::get<1>(id_after_mag))
      //	{
      a_tx = it_hitAfterSim.momX / it_hitAfterSim.momZ;
      a_ty = it_hitAfterSim.momY / it_hitAfterSim.momZ;

      temp_mom =
          TMath::Sqrt(TMath::Sq(it_hitAfterSim.momX) + TMath::Sq(it_hitAfterSim.momY) + TMath::Sq(it_hitAfterSim.momZ));

      a_vx = it_hitAfterSim.momX / temp_mom;
      a_vy = it_hitAfterSim.momY / temp_mom;
      a_vz = it_hitAfterSim.momZ / temp_mom;

      a_x = it_hitAfterSim.hitX;
      a_y = it_hitAfterSim.hitY;
      a_z = it_hitAfterSim.hitZ;

      mom.emplace_back(temp_mom);
      //}
      // else
      // 	att._logger->error("FlatML : Sim After LayerID different from id_after_mag.id_hit {}
      // {}",it_hitAfterSim.layerID, std::get<1>(id_after_mag));

      auto it_hitPSB = it_ListHitsSim->second[std::get<0>(id_psb)];
      // if(it_hitPSB.layerID==std::get<1>(id_psb))
      //{
      psb_z = it_hitPSB.hitZ;
      tof   = it_hitPSB.time;
      //}
      // else
      // 	att._logger->error("FlatML : Sim PSB LayerID different from id_psb.id_hit {} {}",it_hitPSB.layerID,
      // std::get<1>(id_psb));

      auto PDG_particle = TDatabasePDG::Instance()->GetParticle(temp_pdg);
      double charge     = PDG_particle->Charge() / 3.;
      q                 = charge;
      pdg               = temp_pdg;
      poq               = mom[0] / charge;
      qop               = 1. / poq;

      t_flat->Fill();

#ifdef DEBUG_FLATTEN
      att._logger->debug("filled");
#endif

      b_tx = 0., b_ty = 0., b_vx = 0., b_vy = 0., b_vz = 0., b_x = 0., b_y = 0., b_z = 0.;
      a_tx = 0., a_ty = 0., a_vx = 0., a_vy = 0., a_vz = 0., a_x = 0., a_y = 0., a_z = 0.;
      poq = 0., qop = 0.;
      q         = 0.;
      tof       = 0.;
      psb_z     = 0.;
      psce_psbe = -1;
      pdg       = 0;
    }

  return 0;
}
