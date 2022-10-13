#include "TFlatMCOutputML.h"

#include "FullRecoEvent.hh"
#include "ReturnRes.hh"

#include <tuple>

//#define DEBUG_FLATTEN
DataML_momfit::DataML_momfit(const THyphiAttributes& att_, TTree* outT) : DataML(att_, outT)
{
  out_tree->Branch("b_tx", &b_tx, "b_tx/F");
  out_tree->Branch("b_ty", &b_ty, "b_ty/F");
  out_tree->Branch("b_ptnx", &b_ptnx, "b_ptnx/F");
  out_tree->Branch("b_ptny", &b_ptny, "b_ptny/F");
  out_tree->Branch("b_vx", &b_vx, "b_vx/F");
  out_tree->Branch("b_vy", &b_vy, "b_vy/F");
  out_tree->Branch("b_vz", &b_vz, "b_vz/F");
  out_tree->Branch("b_x", &b_x, "b_x/F");
  out_tree->Branch("b_y", &b_y, "b_y/F");
  out_tree->Branch("b_z", &b_z, "b_z/F");

  out_tree->Branch("b_pt", &b_pt, "b_pt/F");
  out_tree->Branch("b_phi", &b_phi, "b_phi/F");
  out_tree->Branch("b_theta", &b_theta, "b_theta/F");

  out_tree->Branch("a_tx", &a_tx, "a_tx/F");
  out_tree->Branch("a_ty", &a_ty, "a_ty/F");
  out_tree->Branch("a_ptnx", &a_ptnx, "a_ptnx/F");
  out_tree->Branch("a_ptny", &a_ptny, "a_ptny/F");
  out_tree->Branch("a_vx", &a_vx, "a_vx/F");
  out_tree->Branch("a_vy", &a_vy, "a_vy/F");
  out_tree->Branch("a_vz", &a_vz, "a_vz/F");
  out_tree->Branch("a_x", &a_x, "a_x/F");
  out_tree->Branch("a_y", &a_y, "a_y/F");
  out_tree->Branch("a_z", &a_z, "a_z/F");

  out_tree->Branch("a_pt", &a_pt, "a_pt/F");
  out_tree->Branch("a_phi", &a_phi, "a_phi/F");
  out_tree->Branch("a_theta", &a_theta, "a_theta/F");

  out_tree->Branch("Dphi",&Dphi, "Dphi/F");
  out_tree->Branch("poq", &poq, "poq/F");
  out_tree->Branch("qop", &qop, "qop/F");
  out_tree->Branch("ptoq", &ptoq, "ptoq/F");
  out_tree->Branch("qopt", &qopt, "qopt/F");
  out_tree->Branch("tof", &tof, "tof/F");
  out_tree->Branch("psce_psbe", &psce_psbe, "psce_psbe/I");
  out_tree->Branch("psb_z", &psb_z, "psb_z/F");
  out_tree->Branch("q", &q, "q/F");
  out_tree->Branch("pdg", &pdg, "pdg/I");
}

void DataML_momfit::FillEvent(FullRecoEvent& REvent)
{
  int ntrack = -1;
  for(auto it_trackInfo : REvent.TrackInfo)
    {
      ntrack++;
#ifdef DEBUG_FLATTEN
      att._logger->debug("track {}", ntrack);
#endif

      const int id_track  = it_trackInfo.first;
      auto it_ListHits    = REvent.TrackDAF.find(id_track);
      auto it_ListHitsSim = REvent.TrackDAFSim.find(id_track);

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
                  // auto temp_hit = REvent.ListHits[id_det][id_hit].get();
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
                  // auto temp_hit = REvent.ListHits[id_det][id_hit].get();
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
                  // auto temp_hit = REvent.ListHits[id_det][id_hit].get();
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
                  // auto temp_hit = REvent.ListHits[id_det][id_hit].get();
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

      auto it_hitBeforeSim = it_ListHitsSim->second[std::get<0>(id_before_mag)][0];
      // if(it_hitBeforeSim.layerID==std::get<1>(id_before_mag))
      //{
      b_tx = it_hitBeforeSim.momX / it_hitBeforeSim.momZ;
      b_ty = it_hitBeforeSim.momY / it_hitBeforeSim.momZ;

      double temp_mom = TMath::Sqrt(TMath::Sq(it_hitBeforeSim.momX) + TMath::Sq(it_hitBeforeSim.momY) +
                                    TMath::Sq(it_hitBeforeSim.momZ));

      b_vx = it_hitBeforeSim.momX / temp_mom;
      b_vy = it_hitBeforeSim.momY / temp_mom;
      b_vz = it_hitBeforeSim.momZ / temp_mom;

      b_pt    = TMath::Sqrt(TMath::Sq(it_hitBeforeSim.momX) + TMath::Sq(it_hitBeforeSim.momY));

      a_ptnx  = it_hitBeforeSim.momX / b_pt;
      a_ptny  = it_hitBeforeSim.momY / b_pt;

      b_phi   = TMath::ATan2(it_hitBeforeSim.momY, it_hitBeforeSim.momX);
      b_theta = TMath::ASin(a_pt / temp_mom);

      b_x = it_hitBeforeSim.hitX;
      b_y = it_hitBeforeSim.hitY;
      b_z = it_hitBeforeSim.hitZ;

      mom.emplace_back(temp_mom);

      temp_pdg = it_hitBeforeSim.pdg;
      //}
      // else
      // 	logger->error("FlatML : Sim Before LayerID different from id_before_mag.id_hit {}
      // {}",it_hitBeforeSim.layerID, std::get<1>(id_before_mag));

      auto it_hitAfterSim = it_ListHitsSim->second[std::get<0>(id_after_mag)][0];
      // if(it_hitAfterSim.layerID==std::get<1>(id_after_mag))
      //	{
      a_tx = it_hitAfterSim.momX / it_hitAfterSim.momZ;
      a_ty = it_hitAfterSim.momY / it_hitAfterSim.momZ;

      temp_mom =
          TMath::Sqrt(TMath::Sq(it_hitAfterSim.momX) + TMath::Sq(it_hitAfterSim.momY) + TMath::Sq(it_hitAfterSim.momZ));

      a_vx = it_hitAfterSim.momX / temp_mom;
      a_vy = it_hitAfterSim.momY / temp_mom;
      a_vz = it_hitAfterSim.momZ / temp_mom;

      a_pt    = TMath::Sqrt(TMath::Sq(it_hitAfterSim.momX) + TMath::Sq(it_hitAfterSim.momY));
      a_ptnx  = it_hitAfterSim.momX / a_pt;
      a_ptny  = it_hitAfterSim.momY / a_pt;
      
      a_phi   = TMath::ATan2(it_hitAfterSim.momY, it_hitAfterSim.momX);
      a_theta = TMath::ASin(a_pt / temp_mom);

      a_x = it_hitAfterSim.hitX;
      a_y = it_hitAfterSim.hitY;
      a_z = it_hitAfterSim.hitZ;

      mom.emplace_back(temp_mom);
      //}
      // else
      // 	logger->error("FlatML : Sim After LayerID different from id_after_mag.id_hit {}
      // {}",it_hitAfterSim.layerID, std::get<1>(id_after_mag));

      auto it_hitPSB = it_ListHitsSim->second[std::get<0>(id_psb)][0];
      // if(it_hitPSB.layerID==std::get<1>(id_psb))
      //{
      psb_z = it_hitPSB.hitZ;
      tof   = it_hitPSB.time;
      //}
      // else
      // 	logger->error("FlatML : Sim PSB LayerID different from id_psb.id_hit {} {}",it_hitPSB.layerID,
      // std::get<1>(id_psb));

      auto PDG_particle = TDatabasePDG::Instance()->GetParticle(temp_pdg);
      double charge     = PDG_particle->Charge() / 3.;       
      q                 = charge;
      pdg               = temp_pdg;
      Dphi              = (b_phi-a_phi)*TMath::RadToDeg() + 360.*((b_phi-a_phi)<-TMath::Pi()) - 360.*((b_phi-a_phi)>TMath::Pi());
      poq               = mom[0] / charge;
      ptoq              = b_pt / charge;
      qop               = 1. / poq;
      qopt              = 1. /ptoq;

      out_tree->Fill();

#ifdef DEBUG_FLATTEN
      att._logger->debug("filled");
#endif

      b_tx = 0.; b_ty = 0.; b_ptnx = 0.; b_ptny = 0.; b_vx = 0.; b_vy = 0.; b_vz = 0.; b_x = 0.; b_y = 0.; b_z = 0.; b_pt = 0.; b_phi = 0.;
      b_theta = 0.;
      a_tx = 0.; a_ty = 0.; a_ptnx = 0.; a_ptny = 0.; a_vx = 0.; a_vy = 0.; a_vz = 0.; a_x = 0.; a_y = 0.; a_z = 0.; a_pt = 0.; a_phi = 0.;
      a_theta = 0.;
      Dphi = 0.;
      poq = 0.; qop = 0.;
      ptoq = 0.; qopt = 0.;
      q         = 0.;
      tof       = 0.;
      psb_z     = 0.;
      psce_psbe = -1;
      pdg       = 0;
    }
}

TFlatMCOutputML::TFlatMCOutputML(const THyphiAttributes& attribut)
    : TDataProcessInterface("FlatMCOutputML"), att(attribut), namefileFlat(att.FlatML_namefile)
{
  f_flat = new TFile(namefileFlat, "RECREATE");
  t_flat = new TTree("MCDataML_Tree", "Flat data event tracks");

  data_out = nullptr;
  if(att.DataML_Out == "DataML_momfit")
    data_out = new DataML_momfit(att, t_flat);
  else
    att._logger->error("FlatMC: selection of the data out for the ML study failed {}, check your config file !!",
                       att.DataML_Out);

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

  delete data_out;

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
  data_out->FillEvent(RecoEvent);

  return 0;
}
