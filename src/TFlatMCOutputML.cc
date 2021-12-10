#include "TFlatMCOutputML.h"

#include "FullRecoEvent.hh"
#include "ReturnRes.hh"

#include <tuple>

//#define DEBUG_FLATTEN
DataML_momfit::DataML_momfit(const THyphiAttributes& att_, TTree* outT) : DataML(att_, outT)
{
  out_tree->Branch("b_tx", &b_tx, "b_tx/F");
  out_tree->Branch("b_ty", &b_ty, "b_ty/F");
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
  out_tree->Branch("a_vx", &a_vx, "a_vx/F");
  out_tree->Branch("a_vy", &a_vy, "a_vy/F");
  out_tree->Branch("a_vz", &a_vz, "a_vz/F");
  out_tree->Branch("a_x", &a_x, "a_x/F");
  out_tree->Branch("a_y", &a_y, "a_y/F");
  out_tree->Branch("a_z", &a_z, "a_z/F");

  out_tree->Branch("a_pt", &a_pt, "a_pt/F");
  out_tree->Branch("a_phi", &a_phi, "a_phi/F");
  out_tree->Branch("a_theta", &a_theta, "a_theta/F");

  out_tree->Branch("poq", &poq, "poq/F");
  out_tree->Branch("qop", &qop, "qop/F");
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
      poq               = mom[0] / charge;
      qop               = 1. / poq;

      out_tree->Fill();

#ifdef DEBUG_FLATTEN
      att._logger->debug("filled");
#endif

      b_tx = 0., b_ty = 0., b_vx = 0., b_vy = 0., b_vz = 0., b_x = 0., b_y = 0., b_z = 0., b_pt = 0., b_phi = 0.,
      b_theta = 0.;
      a_tx = 0., a_ty = 0., a_vx = 0., a_vy = 0., a_vz = 0., a_x = 0., a_y = 0., a_z = 0., a_pt = 0., a_phi = 0.,
      a_theta = 0.;
      poq = 0., qop = 0.;
      q         = 0.;
      tof       = 0.;
      psb_z     = 0.;
      psce_psbe = -1;
      pdg       = 0;
    }
}

DataML_Hyp::DataML_Hyp(const THyphiAttributes& att_, TTree* outT, TTree* outTs, TTree* outTb) : DataML(att_, outT, outTs, outTb)
{
  out_tree->Branch("Pattern", &Pattern, "Pattern/I");

  out_tree->Branch("PDG", &PDG, "PDG/I");
  out_tree->Branch("N_Mother", &N_Mother, "N_Mother/I");
  out_tree->Branch("Chi2ndf", &Chi2ndf, "Chi2ndf/d");
  out_tree->Branch("NDF", &NDF, "NDF/I");
  out_tree->Branch("MomX", &MomX, "MomX/d");
  out_tree->Branch("MomY", &MomY, "MomY/d");
  out_tree->Branch("MomZ", &MomZ, "MomZ/d");
  out_tree->Branch("E", &E, "E/d");
  out_tree->Branch("PrimVtx_PosX", &PrimVtx_PosX, "PrimVtx_PosX/d");
  out_tree->Branch("PrimVtx_PosY", &PrimVtx_PosY, "PrimVtx_PosY/d");
  out_tree->Branch("PrimVtx_PosZ", &PrimVtx_PosZ, "PrimVtx_PosZ/d");
  out_tree->Branch("DecayVtx_PosX", &DecayVtx_PosX, "DecayVtx_PosX/d");
  out_tree->Branch("DecayVtx_PosY", &DecayVtx_PosY, "DecayVtx_PosY/d");
  out_tree->Branch("DecayVtx_PosZ", &DecayVtx_PosZ, "DecayVtx_PosZ/d");
  out_tree->Branch("Dist_RealReconsVtx_X", &Dist_RealReconsVtx_X, "Dist_RealReconsVtx_X/d");
  out_tree->Branch("Dist_RealReconsVtx_Y", &Dist_RealReconsVtx_Y, "Dist_RealReconsVtx_Y/d");
  out_tree->Branch("Dist_RealReconsVtx_Z", &Dist_RealReconsVtx_Z, "Dist_RealReconsVtx_Z/d");
  out_tree->Branch("Dist_MotherPrimVtx", &Dist_MotherPrimVtx, "Dist_MotherPrimVtx/d");
  out_tree->Branch("Angle_MotherPrimVtx", &Angle_MotherPrimVtx, "Angle_MotherPrimVtx/d");
  out_tree->Branch("InvMass", &InvMass, "InvMass/d");
  out_tree->Branch("ErrInvMass", &ErrInvMass, "ErrInvMass/d");
  out_tree->Branch("ErrGetMass", &ErrGetMass, "ErrGetMass/I");
  out_tree->Branch("LifeTime", &LifeTime, "LifeTime/d");
  out_tree->Branch("ErrLifeTime", &ErrLifeTime, "ErrLifeTime/d");
  out_tree->Branch("ErrGetLifeTime", &ErrGetLifeTime, "ErrGetLifeTime/I");
  out_tree->Branch("Mother_IsFromHyp", &Mother_IsFromHyp, "Mother_IsFromHyp/I");

  out_tree->Branch("Id_Fragment", &Id_Fragment, "Id_Fragment/I");
  out_tree->Branch("MomX_Fragment", &MomX_Fragment, "MomX_Fragment/d");
  out_tree->Branch("MomY_Fragment", &MomY_Fragment, "MomY_Fragment/d");
  out_tree->Branch("MomZ_Fragment", &MomZ_Fragment, "MomZ_Fragment/d");
  out_tree->Branch("E_Fragment", &E_Fragment, "E_Fragment/d");
  out_tree->Branch("Chi2ndf_Fragment", &Chi2ndf_Fragment, "Chi2ndf_Fragment/d");
  out_tree->Branch("NDF_Fragment", &NDF_Fragment, "NDF_Fragment/I");
  out_tree->Branch("Pvalue_Fragment", &Pvalue_Fragment, "Pvalue_Fragment/d");
  out_tree->Branch("Angle_MotherFragment", &Angle_MotherFragment, "Angle_MotherFragment/d");
  out_tree->Branch("Fragment_IsFromHyp", &Fragment_IsFromHyp, "Fragment_IsFromHyp/I");

  out_tree->Branch("Id_Pion", &Id_Pion, "Id_Pion/I");
  out_tree->Branch("MomX_Pion", &MomX_Pion, "MomX_Pion/d");
  out_tree->Branch("MomY_Pion", &MomY_Pion, "MomY_Pion/d");
  out_tree->Branch("MomZ_Pion", &MomZ_Pion, "MomZ_Pion/d");
  out_tree->Branch("E_Pion", &E_Pion, "E_Pion/d");
  out_tree->Branch("Chi2ndf_Pion", &Chi2ndf_Pion, "Chi2ndf_Pion/d");
  out_tree->Branch("NDF_Pion", &NDF_Pion, "NDF_Pion/I");
  out_tree->Branch("Pvalue_Pion", &Pvalue_Pion, "Pvalue_Pion/d");
  out_tree->Branch("Angle_MotherPion", &Angle_MotherPion, "Angle_MotherPion/d");
  out_tree->Branch("NHitsMDC_Pion", &NHitsMDC_Pion, "NHitsMDC_Pion/I");
  out_tree->Branch("NHitsMinifiber_Pion", &NHitsMinifiber_Pion, "NHitsMinifiber_Pion/I");
  out_tree->Branch("N_Pion", &N_Pion, "N_Pion/I");
  out_tree->Branch("Pion_IsFromHyp", &Pion_IsFromHyp, "Pion_IsFromHyp/I");

  out_tree->Branch("Dist_Daughters", &Dist_Daughters, "Dist_Daughters/d");
  out_tree->Branch("ArmPod_Qt", &ArmPod_Qt, "ArmPod_Qt/d");
  out_tree->Branch("ArmPod_Alfa", &ArmPod_Alfa, "ArmPod_Alfa/d");



  out_tree_signal->Branch("Pattern", &Pattern, "Pattern/I");

  out_tree_signal->Branch("PDG", &PDG, "PDG/I");
  out_tree_signal->Branch("N_Mother", &N_Mother, "N_Mother/I");
  out_tree_signal->Branch("Chi2ndf", &Chi2ndf, "Chi2ndf/d");
  out_tree_signal->Branch("NDF", &NDF, "NDF/I");
  out_tree_signal->Branch("MomX", &MomX, "MomX/d");
  out_tree_signal->Branch("MomY", &MomY, "MomY/d");
  out_tree_signal->Branch("MomZ", &MomZ, "MomZ/d");
  out_tree_signal->Branch("E", &E, "E/d");
  out_tree_signal->Branch("PrimVtx_PosX", &PrimVtx_PosX, "PrimVtx_PosX/d");
  out_tree_signal->Branch("PrimVtx_PosY", &PrimVtx_PosY, "PrimVtx_PosY/d");
  out_tree_signal->Branch("PrimVtx_PosZ", &PrimVtx_PosZ, "PrimVtx_PosZ/d");
  out_tree_signal->Branch("DecayVtx_PosX", &DecayVtx_PosX, "DecayVtx_PosX/d");
  out_tree_signal->Branch("DecayVtx_PosY", &DecayVtx_PosY, "DecayVtx_PosY/d");
  out_tree_signal->Branch("DecayVtx_PosZ", &DecayVtx_PosZ, "DecayVtx_PosZ/d");
  out_tree_signal->Branch("Dist_RealReconsVtx_X", &Dist_RealReconsVtx_X, "Dist_RealReconsVtx_X/d");
  out_tree_signal->Branch("Dist_RealReconsVtx_Y", &Dist_RealReconsVtx_Y, "Dist_RealReconsVtx_Y/d");
  out_tree_signal->Branch("Dist_RealReconsVtx_Z", &Dist_RealReconsVtx_Z, "Dist_RealReconsVtx_Z/d");
  out_tree_signal->Branch("Dist_MotherPrimVtx", &Dist_MotherPrimVtx, "Dist_MotherPrimVtx/d");
  out_tree_signal->Branch("Angle_MotherPrimVtx", &Angle_MotherPrimVtx, "Angle_MotherPrimVtx/d");
  out_tree_signal->Branch("InvMass", &InvMass, "InvMass/d");
  out_tree_signal->Branch("ErrInvMass", &ErrInvMass, "ErrInvMass/d");
  out_tree_signal->Branch("ErrGetMass", &ErrGetMass, "ErrGetMass/I");
  out_tree_signal->Branch("LifeTime", &LifeTime, "LifeTime/d");
  out_tree_signal->Branch("ErrLifeTime", &ErrLifeTime, "ErrLifeTime/d");
  out_tree_signal->Branch("ErrGetLifeTime", &ErrGetLifeTime, "ErrGetLifeTime/I");
  out_tree_signal->Branch("Mother_IsFromHyp", &Mother_IsFromHyp, "Mother_IsFromHyp/I");

  out_tree_signal->Branch("Id_Fragment", &Id_Fragment, "Id_Fragment/I");
  out_tree_signal->Branch("MomX_Fragment", &MomX_Fragment, "MomX_Fragment/d");
  out_tree_signal->Branch("MomY_Fragment", &MomY_Fragment, "MomY_Fragment/d");
  out_tree_signal->Branch("MomZ_Fragment", &MomZ_Fragment, "MomZ_Fragment/d");
  out_tree_signal->Branch("E_Fragment", &E_Fragment, "E_Fragment/d");
  out_tree_signal->Branch("Chi2ndf_Fragment", &Chi2ndf_Fragment, "Chi2ndf_Fragment/d");
  out_tree_signal->Branch("NDF_Fragment", &NDF_Fragment, "NDF_Fragment/I");
  out_tree_signal->Branch("Pvalue_Fragment", &Pvalue_Fragment, "Pvalue_Fragment/d");
  out_tree_signal->Branch("Angle_MotherFragment", &Angle_MotherFragment, "Angle_MotherFragment/d");
  out_tree_signal->Branch("Fragment_IsFromHyp", &Fragment_IsFromHyp, "Fragment_IsFromHyp/I");

  out_tree_signal->Branch("Id_Pion", &Id_Pion, "Id_Pion/I");
  out_tree_signal->Branch("MomX_Pion", &MomX_Pion, "MomX_Pion/d");
  out_tree_signal->Branch("MomY_Pion", &MomY_Pion, "MomY_Pion/d");
  out_tree_signal->Branch("MomZ_Pion", &MomZ_Pion, "MomZ_Pion/d");
  out_tree_signal->Branch("E_Pion", &E_Pion, "E_Pion/d");
  out_tree_signal->Branch("Chi2ndf_Pion", &Chi2ndf_Pion, "Chi2ndf_Pion/d");
  out_tree_signal->Branch("NDF_Pion", &NDF_Pion, "NDF_Pion/I");
  out_tree_signal->Branch("Pvalue_Pion", &Pvalue_Pion, "Pvalue_Pion/d");
  out_tree_signal->Branch("Angle_MotherPion", &Angle_MotherPion, "Angle_MotherPion/d");
  out_tree_signal->Branch("NHitsMDC_Pion", &NHitsMDC_Pion, "NHitsMDC_Pion/I");
  out_tree_signal->Branch("NHitsMinifiber_Pion", &NHitsMinifiber_Pion, "NHitsMinifiber_Pion/I");
  out_tree_signal->Branch("N_Pion", &N_Pion, "N_Pion/I");
  out_tree_signal->Branch("Pion_IsFromHyp", &Pion_IsFromHyp, "Pion_IsFromHyp/I");

  out_tree_signal->Branch("Dist_Daughters", &Dist_Daughters, "Dist_Daughters/d");
  out_tree_signal->Branch("ArmPod_Qt", &ArmPod_Qt, "ArmPod_Qt/d");
  out_tree_signal->Branch("ArmPod_Alfa", &ArmPod_Alfa, "ArmPod_Alfa/d");



  out_tree_background->Branch("Pattern", &Pattern, "Pattern/I");

  out_tree_background->Branch("PDG", &PDG, "PDG/I");
  out_tree_background->Branch("N_Mother", &N_Mother, "N_Mother/I");
  out_tree_background->Branch("Chi2ndf", &Chi2ndf, "Chi2ndf/d");
  out_tree_background->Branch("NDF", &NDF, "NDF/I");
  out_tree_background->Branch("MomX", &MomX, "MomX/d");
  out_tree_background->Branch("MomY", &MomY, "MomY/d");
  out_tree_background->Branch("MomZ", &MomZ, "MomZ/d");
  out_tree_background->Branch("E", &E, "E/d");
  out_tree_background->Branch("PrimVtx_PosX", &PrimVtx_PosX, "PrimVtx_PosX/d");
  out_tree_background->Branch("PrimVtx_PosY", &PrimVtx_PosY, "PrimVtx_PosY/d");
  out_tree_background->Branch("PrimVtx_PosZ", &PrimVtx_PosZ, "PrimVtx_PosZ/d");
  out_tree_background->Branch("DecayVtx_PosX", &DecayVtx_PosX, "DecayVtx_PosX/d");
  out_tree_background->Branch("DecayVtx_PosY", &DecayVtx_PosY, "DecayVtx_PosY/d");
  out_tree_background->Branch("DecayVtx_PosZ", &DecayVtx_PosZ, "DecayVtx_PosZ/d");
  out_tree_background->Branch("Dist_RealReconsVtx_X", &Dist_RealReconsVtx_X, "Dist_RealReconsVtx_X/d");
  out_tree_background->Branch("Dist_RealReconsVtx_Y", &Dist_RealReconsVtx_Y, "Dist_RealReconsVtx_Y/d");
  out_tree_background->Branch("Dist_RealReconsVtx_Z", &Dist_RealReconsVtx_Z, "Dist_RealReconsVtx_Z/d");
  out_tree_background->Branch("Dist_MotherPrimVtx", &Dist_MotherPrimVtx, "Dist_MotherPrimVtx/d");
  out_tree_background->Branch("Angle_MotherPrimVtx", &Angle_MotherPrimVtx, "Angle_MotherPrimVtx/d");
  out_tree_background->Branch("InvMass", &InvMass, "InvMass/d");
  out_tree_background->Branch("ErrInvMass", &ErrInvMass, "ErrInvMass/d");
  out_tree_background->Branch("ErrGetMass", &ErrGetMass, "ErrGetMass/I");
  out_tree_background->Branch("LifeTime", &LifeTime, "LifeTime/d");
  out_tree_background->Branch("ErrLifeTime", &ErrLifeTime, "ErrLifeTime/d");
  out_tree_background->Branch("ErrGetLifeTime", &ErrGetLifeTime, "ErrGetLifeTime/I");
  out_tree_background->Branch("Mother_IsFromHyp", &Mother_IsFromHyp, "Mother_IsFromHyp/I");

  out_tree_background->Branch("Id_Fragment", &Id_Fragment, "Id_Fragment/I");
  out_tree_background->Branch("MomX_Fragment", &MomX_Fragment, "MomX_Fragment/d");
  out_tree_background->Branch("MomY_Fragment", &MomY_Fragment, "MomY_Fragment/d");
  out_tree_background->Branch("MomZ_Fragment", &MomZ_Fragment, "MomZ_Fragment/d");
  out_tree_background->Branch("E_Fragment", &E_Fragment, "E_Fragment/d");
  out_tree_background->Branch("Chi2ndf_Fragment", &Chi2ndf_Fragment, "Chi2ndf_Fragment/d");
  out_tree_background->Branch("NDF_Fragment", &NDF_Fragment, "NDF_Fragment/I");
  out_tree_background->Branch("Pvalue_Fragment", &Pvalue_Fragment, "Pvalue_Fragment/d");
  out_tree_background->Branch("Angle_MotherFragment", &Angle_MotherFragment, "Angle_MotherFragment/d");
  out_tree_background->Branch("Fragment_IsFromHyp", &Fragment_IsFromHyp, "Fragment_IsFromHyp/I");

  out_tree_background->Branch("Id_Pion", &Id_Pion, "Id_Pion/I");
  out_tree_background->Branch("MomX_Pion", &MomX_Pion, "MomX_Pion/d");
  out_tree_background->Branch("MomY_Pion", &MomY_Pion, "MomY_Pion/d");
  out_tree_background->Branch("MomZ_Pion", &MomZ_Pion, "MomZ_Pion/d");
  out_tree_background->Branch("E_Pion", &E_Pion, "E_Pion/d");
  out_tree_background->Branch("Chi2ndf_Pion", &Chi2ndf_Pion, "Chi2ndf_Pion/d");
  out_tree_background->Branch("NDF_Pion", &NDF_Pion, "NDF_Pion/I");
  out_tree_background->Branch("Pvalue_Pion", &Pvalue_Pion, "Pvalue_Pion/d");
  out_tree_background->Branch("Angle_MotherPion", &Angle_MotherPion, "Angle_MotherPion/d");
  out_tree_background->Branch("NHitsMDC_Pion", &NHitsMDC_Pion, "NHitsMDC_Pion/I");
  out_tree_background->Branch("NHitsMinifiber_Pion", &NHitsMinifiber_Pion, "NHitsMinifiber_Pion/I");
  out_tree_background->Branch("N_Pion", &N_Pion, "N_Pion/I");
  out_tree_background->Branch("Pion_IsFromHyp", &Pion_IsFromHyp, "Pion_IsFromHyp/I");

  out_tree_background->Branch("Dist_Daughters", &Dist_Daughters, "Dist_Daughters/d");
  out_tree_background->Branch("ArmPod_Qt", &ArmPod_Qt, "ArmPod_Qt/d");
  out_tree_background->Branch("ArmPod_Alfa", &ArmPod_Alfa, "ArmPod_Alfa/d");
}

void DataML_Hyp::FillEvent(FullRecoEvent& REvent)
{
  for(size_t i = 0; i < REvent.Hyp_Vect.size(); ++i)
    {
      //if(REvent.Hyp_Vect[i].Pattern != 4) //Save only KFParticle reconstructed mothers
      //  continue;

      if(std::isfinite(REvent.Hyp_Vect[i].Chi2ndf) == false)
        continue;
      if(std::isfinite(REvent.Hyp_Vect[i].Dist_Daughters) == false)
        continue;
      if(std::isfinite(REvent.Hyp_Vect[i].MomE_Fragment.Px()) == false)
        continue;
      

      Pattern = REvent.Hyp_Vect[i].Pattern; 
      
      PDG  = REvent.Hyp_Vect[i].PDG;
      N_Mother = REvent.Hyp_Vect[i].N_Mother;
      Chi2ndf = REvent.Hyp_Vect[i].Chi2ndf;
      NDF = REvent.Hyp_Vect[i].NDF;
      MomX = REvent.Hyp_Vect[i].MomE.Px();
      MomY = REvent.Hyp_Vect[i].MomE.Py();
      MomZ = REvent.Hyp_Vect[i].MomE.Pz();
      E = REvent.Hyp_Vect[i].MomE.E();
      PrimVtx_PosX = REvent.Hyp_Vect[i].PrimVtx.X();
      PrimVtx_PosY = REvent.Hyp_Vect[i].PrimVtx.Y();
      PrimVtx_PosZ = REvent.Hyp_Vect[i].PrimVtx.Z();
      DecayVtx_PosX = REvent.Hyp_Vect[i].DecayVtx.X();
      DecayVtx_PosY = REvent.Hyp_Vect[i].DecayVtx.Y();
      DecayVtx_PosZ = REvent.Hyp_Vect[i].DecayVtx.Z();
      Dist_RealReconsVtx_X = REvent.Hyp_Vect[i].Dist_RealReconsVtx.X();
      Dist_RealReconsVtx_Y = REvent.Hyp_Vect[i].Dist_RealReconsVtx.Y();
      Dist_RealReconsVtx_Z = REvent.Hyp_Vect[i].Dist_RealReconsVtx.Z();
      Dist_MotherPrimVtx = REvent.Hyp_Vect[i].Dist_MotherPrimVtx;
      Angle_MotherPrimVtx = REvent.Hyp_Vect[i].Angle_MotherPrimVtx;
      InvMass = REvent.Hyp_Vect[i].InvMass;
      ErrInvMass = REvent.Hyp_Vect[i].ErrInvMass;
      ErrGetMass = REvent.Hyp_Vect[i].ErrGetMass;
      LifeTime = REvent.Hyp_Vect[i].LifeTime;
      ErrLifeTime = REvent.Hyp_Vect[i].ErrLifeTime;
      ErrGetLifeTime = REvent.Hyp_Vect[i].ErrGetLifeTime;
      Mother_IsFromHyp = REvent.Hyp_Vect[i].Mother_IsFromHyp;

      Id_Fragment = REvent.Hyp_Vect[i].Id_Fragment;
      MomX_Fragment = REvent.Hyp_Vect[i].MomE_Fragment.Px();
      MomY_Fragment = REvent.Hyp_Vect[i].MomE_Fragment.Py();
      MomZ_Fragment = REvent.Hyp_Vect[i].MomE_Fragment.Pz();
      E_Fragment = REvent.Hyp_Vect[i].MomE_Fragment.E();
      Chi2ndf_Fragment = REvent.Hyp_Vect[i].Chi2ndf_Fragment;
      NDF_Fragment = REvent.Hyp_Vect[i].NDF_Fragment;
      Pvalue_Fragment = REvent.Hyp_Vect[i].Pvalue_Fragment;
      Angle_MotherFragment = REvent.Hyp_Vect[i].Angle_MotherFragment;
      Fragment_IsFromHyp = REvent.Hyp_Vect[i].Fragment_IsFromHyp;
      
      Id_Pion = REvent.Hyp_Vect[i].Id_Pion;
      MomX_Pion = REvent.Hyp_Vect[i].MomE_Pion.Px();
      MomY_Pion = REvent.Hyp_Vect[i].MomE_Pion.Py();
      MomZ_Pion = REvent.Hyp_Vect[i].MomE_Pion.Pz();
      E_Pion = REvent.Hyp_Vect[i].MomE_Pion.E();
      Chi2ndf_Pion = REvent.Hyp_Vect[i].Chi2ndf_Pion;
      NDF_Pion = REvent.Hyp_Vect[i].NDF_Pion;
      Pvalue_Pion = REvent.Hyp_Vect[i].Pvalue_Pion;
      Angle_MotherPion = REvent.Hyp_Vect[i].Angle_MotherPion;
      NHitsMDC_Pion = REvent.Hyp_Vect[i].NHitsMDC_Pion;
      NHitsMinifiber_Pion = REvent.Hyp_Vect[i].NHitsMinifiber_Pion;
      N_Pion = REvent.Hyp_Vect[i].N_Pion;
      Pion_IsFromHyp = REvent.Hyp_Vect[i].Pion_IsFromHyp;

      Dist_Daughters = REvent.Hyp_Vect[i].Dist_Daughters;
      ArmPod_Qt = REvent.Hyp_Vect[i].ArmPod_Qt;
      ArmPod_Alfa = REvent.Hyp_Vect[i].ArmPod_Alfa;


      out_tree->Fill();

      if(SplitTree)
        {
          if(Mother_IsFromHyp == 1)
            out_tree_signal->Fill();
          else if(Mother_IsFromHyp == 0)
            out_tree_background->Fill();
        }
    }
}


TFlatMCOutputML::TFlatMCOutputML(const THyphiAttributes& attribut)
    : TDataProcessInterface("FlatMCOutputML"), att(attribut), namefileFlat(att.FlatML_namefile)
{
  f_flat = new TFile(namefileFlat, "RECREATE");
  t_flat = new TTree("MCDataML_Tree", "Flat data");
  t_flatS = new TTree("MCDataML_TreeS", "Flat signal data");
  t_flatB = new TTree("MCDataML_TreeB", "Flat background data");

  data_out = nullptr;
  if(att.DataML_Out == "DataML_momfit")
    data_out = new DataML_momfit(att, t_flat);
  else if(att.DataML_Out == "DataML_Hyp")
    data_out = new DataML_Hyp(att, t_flat, t_flatS, t_flatB);
  else if(att.DataML_Out == "DataML_Hyp_Split")
    {
      data_out = new DataML_Hyp(att, t_flat, t_flatS, t_flatB);
      SplitTree = true;
    }
  else
    att._logger->error("FlatMC: selection of the data out for the ML study failed {}, check your config file !!",
                       att.DataML_Out);

  t_flat ->SetDirectory(f_flat);
  t_flatS->SetDirectory(f_flat);
  t_flatB->SetDirectory(f_flat);
  // t_flat->AutoSave();

  att._logger->info("FlatMC : tree out set ");
}

TFlatMCOutputML::~TFlatMCOutputML()
{
  f_flat->cd();
  t_flat->Write();

  if(SplitTree)
    {
      t_flatS->Write();
      t_flatB->Write();
    }
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
