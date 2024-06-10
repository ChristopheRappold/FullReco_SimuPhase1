#include "TImproveHitsFiber.h"

#include "Ana_Event/Ana_WasaEvent.hh"
#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "FullRecoEvent.hh"
#include "KalmanFittedStateOnPlane.h"
#include "KalmanFitterInfo.h"
#include "PlanarMeasurement.h"
#include "ProlateSpacepointMeasurement.h"
#include "ReturnRes.hh"
#include "StateOnPlane.h"
#include "TLinearFitter.h"
#include "TRandom3.h"
#include "UnscentedTrans.hh"

#include <iterator>
#include <set>
#include <sstream>
#include <tuple>

//#define DEBUG_RIEMANNFINDER

using namespace std;
using namespace G4Sol;

template <class Out>
TImproveHitsFiber<Out>::TImproveHitsFiber(const THyphiAttributes& attribut): TDataProcessInterface<Out>("ImproveHitsFiber"), att(attribut)
{

  // LastFrontWall = att.Wasa_Side == 0 ? G4Sol::PSFE : G4Sol::PSBE;

  // res_h.time.resize(1,1);

  att._logger->info("ImproveHits : set ");
}

template <class Out>
TImproveHitsFiber<Out>::~TImproveHitsFiber()
{
}

template <class Out>
void TImproveHitsFiber<Out>::InitMT()
{
  att._logger->error("E> Not supposed to be multithreaded !");
}

template <class Out>
ReturnRes::InfoM TImproveHitsFiber<Out>::operator()(FullRecoEvent& RecoEvent, Out* OutTree)
{
  int result_finder = Exec(RecoEvent, OutTree);

  return SoftExit(result_finder);
}

template <class Out>
int TImproveHitsFiber<Out>::Exec(FullRecoEvent& RecoEvent, Out* OutTree)
{

  int res = ImproveFiber(RecoEvent);
  return res;
}

template <class Out>
ReturnRes::InfoM TImproveHitsFiber<Out>::SoftExit(int)
{
  return ReturnRes::Fine;
}

template <class Out>
void TImproveHitsFiber<Out>::SelectHists()
{
  LocalHisto.h_ImproveHitsResidus  = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ImproveHitsResidus);
  LocalHisto.h_ImproveHitsResidus2 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ImproveHitsResidus2);
}

template <class Out>
int TImproveHitsFiber<Out>::ImproveFiber(FullRecoEvent& RecoEvent)
{

  for(const auto& it_trackDAF : RecoEvent.TrackDAF)
    {
      const int id_track         = it_trackDAF.first;
      const auto& it_ListHits    = it_trackDAF.second;
      const auto& it_ListHitsSim = RecoEvent.TrackDAFSim.find(id_track)->second;

      // auto it_paramFitRZ = RecoEvent.paramFitRZ.find(id_track);
      // if(it_paramFitRZ == RecoEvent.paramFitRZ.end())
      // 	{
      // 	  att._logger->debug("E> ImproveHit : paramFitRZ for idTrack {} does not exist !",id_track);
      // 	  continue;
      // 	}
      // att._logger->debug(" ImproveHit : track {} paramRZ : {} {} ",id_track,
      // it_paramFitRZ->second.p0,it_paramFitRZ->second.p1);

      std::vector<Improve::miniF> new_miniFibers[2];

      std::string TrackLog("");
      for(int id_det = 0; id_det < it_ListHits.size(); ++id_det)
        {
          int id_hit = it_ListHits[id_det];

          if(id_hit < 0)
            continue;

          TrackLog += G4Sol::nameLiteralDet.begin()[id_det];
          TrackLog += "_hit";
          TrackLog += std::to_string(id_hit);
          TrackLog += "->";
        }
      att._logger->debug("track #{} : {}", id_track, TrackLog);

      for(int id_det = 0; id_det < it_ListHits.size(); ++id_det)
        {
          int id_hit = it_ListHits[id_det];

          if(id_hit < 0)
            continue;
          if(!(id_det >= G4Sol::MiniFiberD1_x &&
               id_det <= G4Sol::MiniFiberD2_u)) // || id_det == G4Sol::PSFE || id_det == G4Sol::PSBE))
            continue;

          att._logger->debug(" Improve hit > det : {} {}, hit {} ", id_det, G4Sol::nameLiteralDet.begin()[id_det],
                             id_hit);

          RecoEvent.OldListHits[id_det].resize(RecoEvent.ListHits[id_det].size());
          if(id_det >= G4Sol::MiniFiberD1_x && id_det <= G4Sol::MiniFiberD1_v)
            {
              att._logger->debug("miniFiber 1");
              new_miniFibers[0].emplace_back(
					     ImproveFiber2(RecoEvent.ListHits, id_det, id_hit)); //,it_paramFitRZ->second));
            }
          if(id_det >= G4Sol::MiniFiberD2_x && id_det <= G4Sol::MiniFiberD2_u)
            {
              att._logger->debug("miniFiber 2");
              new_miniFibers[1].emplace_back(
					     ImproveFiber2(RecoEvent.ListHits, id_det, id_hit)); //,it_paramFitRZ->second));
            }

          for(int id_mini = 0; id_mini < 2; ++id_mini)
            if(new_miniFibers[id_mini].size() == 3)
              {
                att._logger->debug(" case 3: minifiber {}", id_mini + 1);
                for(int id_M = 0; id_M < 3; ++id_M)
                  {
                    int id_det = new_miniFibers[id_mini][id_M].det;
                    int id_hit = new_miniFibers[id_mini][id_M].hit;
                    att._logger->debug(
				       "  det {} hit {} posR1 [ {}, {}, {} ] posR2 [ {}, {}, {} ]",
				       G4Sol::nameLiteralDet.begin()[id_det], new_miniFibers[id_mini][id_M].hit,
				       new_miniFibers[id_mini][id_M].xyz[0].X(), new_miniFibers[id_mini][id_M].xyz[0].Y(),
				       new_miniFibers[id_mini][id_M].xyz[0].Z(), new_miniFibers[id_mini][id_M].xyz[1].X(),
				       new_miniFibers[id_mini][id_M].xyz[1].Y(), new_miniFibers[id_mini][id_M].xyz[1].Z());

                    if(it_ListHitsSim[id_det].size() != 0)
                      {
                        std::array<double, 3> sim_mean_xyz{0., 0., 0.};
                        double size_mean = 0.;
                        for(const auto& SimH : it_ListHitsSim[id_det])
                          {
                            sim_mean_xyz[0] += 0.5 * (SimH.hitX + SimH.hitXexit);
                            sim_mean_xyz[1] += 0.5 * (SimH.hitY + SimH.hitYexit);
                            sim_mean_xyz[2] += 0.5 * (SimH.hitZ + SimH.hitZexit);
                            size_mean += 1.;
                          }

                        sim_mean_xyz[0] /= size_mean;
                        sim_mean_xyz[1] /= size_mean;
                        sim_mean_xyz[2] /= size_mean;

                        att._logger->debug("SImHit : det {} hit {}, [{}, {}, {}]",
                                           G4Sol::nameLiteralDet.begin()[id_det], id_hit, sim_mean_xyz[0],
                                           sim_mean_xyz[1], sim_mean_xyz[2]);
                      }
                  }

                // auto temp01_11 = new_miniFibers[id_mini][0].xyz[0]-new_miniFibers[id_mini][1].xyz[0] ;
                // auto temp01_12 = new_miniFibers[id_mini][0].xyz[0]-new_miniFibers[id_mini][1].xyz[1] ;

                // auto temp01_21 = new_miniFibers[id_mini][0].xyz[1]-new_miniFibers[id_mini][1].xyz[0] ;
                // auto temp01_22 = new_miniFibers[id_mini][0].xyz[1]-new_miniFibers[id_mini][1].xyz[1] ;

                auto [x12, y12] = Intersection(new_miniFibers[id_mini][0], new_miniFibers[id_mini][1]);
                auto [x23, y23] = Intersection(new_miniFibers[id_mini][1], new_miniFibers[id_mini][2]);
                auto [x13, y13] = Intersection(new_miniFibers[id_mini][0], new_miniFibers[id_mini][2]);

                TLinearFitter fitXYZ(2, "pol1");
                fitXYZ.StoreData(false);

                std::array<double, 3> posZ{
		  0.5 * (new_miniFibers[id_mini][0].xyz[0].Z() + new_miniFibers[id_mini][1].xyz[0].Z()),
		  0.5 * (new_miniFibers[id_mini][1].xyz[0].Z() + new_miniFibers[id_mini][2].xyz[0].Z()),
		  0.5 * (new_miniFibers[id_mini][0].xyz[0].Z() + new_miniFibers[id_mini][2].xyz[0].Z())};

                std::array<double, 3> posX{x12, x23, x13};
                std::array<double, 3> posY{y12, y23, y13};

                att._logger->debug(" Fiber Intersection12 : {} {} {}", x12, y12, posZ[0]);
                att._logger->debug(" Fiber Intersection23 : {} {} {}", x23, y23, posZ[1]);
                att._logger->debug(" Fiber Intersection13 : {} {} {}", x13, y13, posZ[2]);

                fitXYZ.AssignData(3, 1, posZ.data(), posX.data());
                int statusFitXZ = fitXYZ.Eval();
                TVectorD parametersXZ;
                fitXYZ.GetParameters(parametersXZ); // X = [0]+Z*[1];
                att._logger->debug("{} {} {} | {} {} {}", posZ[0], posZ[1], posZ[2], posY[0], posY[1], posY[2]);
                TLinearFitter fitXYZ2(2, "pol1");
                fitXYZ2.StoreData(false);
                fitXYZ2.AssignData(3, 1, posZ.data(), posY.data());
                int statusFitYZ = fitXYZ2.Eval();
                TVectorD parametersYZ;
                fitXYZ2.GetParameters(parametersYZ); // Y = [0]+Z*[1];

                att._logger->debug(" XZ : {} {} ", parametersXZ[0], parametersXZ[1]);
                att._logger->debug(" YZ : {} {} ", parametersYZ[0], parametersYZ[1]);

                // std::array<double,4> difIall = { temp01_11.Mag(), temp01_12.Mag(), temp01_21.Mag(), temp01_22.Mag()};
                // auto it_min = TMath::LocMin(difIall.begin(),difIall.end());
                // int id_min = std::distance(difIall.begin(),it_min);
                // if(id_min == 1 || id_min ==3)
                // att._logger->error("E> ImproveHit : MiniFiber hit coincidence with different R on different layers
                // (case size=3)| track {}", id_track); else
		// int id_min_red = id_min/2;
		// auto temp02_11 =
		// new_miniFibers[id_mini][0].xyz[id_min_red]-new_miniFibers[id_mini][2].xyz[id_min_red] ; auto
		// temp12_11 = new_miniFibers[id_mini][1].xyz[id_min_red]-new_miniFibers[id_mini][2].xyz[id_min_red] ;

		// double difI1 = temp02_11.Mag();
		// double difI2 = temp02_11.Mag();

		// if(TMath::Abs(difI1-*it_min)> 0.5 || TMath::Abs(difI2-*it_min) > 0.5 ) //cm
		//   {
		//     att._logger->error("E> ImproveHit : MiniFiber third hit at the same positions than others |
		//     track {}", id_track);
		//   }
		// else
		//   {
		for(int i_xuv = 0; i_xuv < 3; ++i_xuv)
		  {
		    int id_det = new_miniFibers[id_mini][i_xuv].det;
		    int id_hit = new_miniFibers[id_mini][i_xuv].hit;

		    genfit::PlanarMeasurement* currentHit =
		      dynamic_cast<genfit::PlanarMeasurement*>(RecoEvent.ListHits[id_det][id_hit].get());
		    auto dummyState   = genfit::StateOnPlane();
		    auto currentPlane = currentHit->constructPlane(genfit::StateOnPlane());
		    auto HitCoord     = currentHit->getRawHitCoords();
		    auto HitCov       = currentHit->getRawHitCov();

		    double posZ = new_miniFibers[id_mini][i_xuv].xyz[0].Z();

		    TVectorD hitCoordsNew(2);
		    // hitCoordsNew(0) = new_miniFibers[id_mini][i_xuv].xyz[id_min_red].X();
		    // hitCoordsNew(1) = new_miniFibers[id_mini][i_xuv].xyz[id_min_red].Y();
		    hitCoordsNew(0) = parametersXZ[1] * posZ + parametersXZ[0];
		    hitCoordsNew(1) = parametersYZ[1] * posZ + parametersYZ[0];
		    TMatrixDSym hitCovNew(2);
		    hitCovNew(0, 0) = HitCov(0, 0);
		    hitCovNew(0, 1) = 0.;
		    hitCovNew(1, 0) = 0.;
		    hitCovNew(1, 1) = 1.2 * HitCov(0, 0);

		    TVector3 o(0, 0, posZ);
		    TVector3 u(1., 0., 0.), v(0., 1., 0.);
		    genfit::SharedPlanePtr newPlane(new genfit::DetPlane(o, u, v));

		    RecoEvent.OldListHits[id_det][id_hit] =
		      std::make_unique<genfit::PlanarMeasurement>(hitCoordsNew, hitCovNew, id_det, id_hit, nullptr);
		    dynamic_cast<genfit::PlanarMeasurement*>(RecoEvent.OldListHits[id_det][id_hit].get())
		      ->setPlane(newPlane);

		    RecoEvent.OldListHits[id_det][id_hit].swap(RecoEvent.ListHits[id_det][id_hit]);

		    if(it_ListHitsSim[id_det].size() != 0)
		      {
			std::array<double, 3> sim_mean_xyz{0., 0., 0.};
			double size_mean = 0.;
			for(const auto& SimH : it_ListHitsSim[id_det])
			  {
			    sim_mean_xyz[0] += 0.5 * (SimH.hitX + SimH.hitXexit);
			    sim_mean_xyz[1] += 0.5 * (SimH.hitY + SimH.hitYexit);
			    sim_mean_xyz[2] += 0.5 * (SimH.hitZ + SimH.hitZexit);
			    size_mean += 1.;
			  }

			sim_mean_xyz[0] /= size_mean;
			sim_mean_xyz[1] /= size_mean;
			sim_mean_xyz[2] /= size_mean;

			att._logger->debug("SImHit : det {} hit {}, [{}, {}, {}] | Est : [{}, {}, {}]",
					   G4Sol::nameLiteralDet.begin()[id_det], id_hit, sim_mean_xyz[0],
					   sim_mean_xyz[1], sim_mean_xyz[2], hitCoordsNew(0), hitCoordsNew(1), posZ);

			TString nameTemp1(G4Sol::nameLiteralDet.begin()[id_det]);
			LocalHisto.h_ImproveHitsResidus->Fill(nameTemp1 + "_x", sim_mean_xyz[0] - hitCoordsNew(0),
							      1.);
			LocalHisto.h_ImproveHitsResidus->Fill(nameTemp1 + "_y", sim_mean_xyz[1] - hitCoordsNew(1),
							      1.);
			LocalHisto.h_ImproveHitsResidus->Fill(nameTemp1 + "_z", sim_mean_xyz[2] - posZ, 1.);
			LocalHisto.h_ImproveHitsResidus->Fill(
							      nameTemp1 + "_r",
							      TMath::Hypot(sim_mean_xyz[0] - hitCoordsNew(0), sim_mean_xyz[1] - hitCoordsNew(1)), 1.);
			LocalHisto.h_ImproveHitsResidus->Fill(
							      nameTemp1 + "_r",
							      TMath::Sqrt(TMath::Sq(sim_mean_xyz[0] - hitCoordsNew(0)) +
									  TMath::Sq(sim_mean_xyz[1] - hitCoordsNew(1)) +
									  TMath::Sq(sim_mean_xyz[2] - posZ)),
							      1.);

			// double Rest = it_paramFitRZ->second.extrapR(posZ);

			// LocalHisto.h_ImproveHitsResidus->Fill(nameTemp1+"_fitRZ",Rest-TMath::Hypot(hitCoordsNew(0),hitCoordsNew(1)),1.);

			LocalHisto.h_ImproveHitsResidus2->Fill(nameTemp1 + "_x", sim_mean_xyz[0] - hitCoordsNew(0),
							       1.);
			LocalHisto.h_ImproveHitsResidus2->Fill(nameTemp1 + "_y", sim_mean_xyz[1] - hitCoordsNew(1),
							       1.);
			LocalHisto.h_ImproveHitsResidus2->Fill(nameTemp1 + "_z", sim_mean_xyz[2] - posZ, 1.);
			LocalHisto.h_ImproveHitsResidus2->Fill(
							       nameTemp1 + "_r",
							       TMath::Hypot(sim_mean_xyz[0] - hitCoordsNew(0), sim_mean_xyz[1] - hitCoordsNew(1)), 1.);
			LocalHisto.h_ImproveHitsResidus2->Fill(
							       nameTemp1 + "_r",
							       TMath::Sqrt(TMath::Sq(sim_mean_xyz[0] - hitCoordsNew(0)) +
									   TMath::Sq(sim_mean_xyz[1] - hitCoordsNew(1)) +
									   TMath::Sq(sim_mean_xyz[2] - posZ)),
							       1.);
			// LocalHisto.h_ImproveHitsResidus2->Fill(nameTemp1+"_fitRZ",Rest-TMath::Hypot(hitCoordsNew(0),hitCoordsNew(1)),1.);
		      }

		    //}
		  }
              }
            else if(new_miniFibers[id_mini].size() == 2)
              {
                att._logger->debug(" case 2: minifiber {}", new_miniFibers[id_mini].size());
                for(int id_M = 0; id_M < 2; ++id_M)
                  {
                    int id_det = new_miniFibers[id_mini][id_M].det;
                    int id_hit = new_miniFibers[id_mini][id_M].hit;
                    att._logger->debug(
				       "  det {} hit {} posR1 [ {}, {}, {} ] posR2 [ {}, {}, {} ]",
				       G4Sol::nameLiteralDet.begin()[id_det], new_miniFibers[id_mini][id_M].hit,
				       new_miniFibers[id_mini][id_M].xyz[0].X(), new_miniFibers[id_mini][id_M].xyz[0].Y(),
				       new_miniFibers[id_mini][id_M].xyz[0].Z(), new_miniFibers[id_mini][id_M].xyz[1].X(),
				       new_miniFibers[id_mini][id_M].xyz[1].Y(), new_miniFibers[id_mini][id_M].xyz[1].Z());

                    if(it_ListHitsSim[id_det].size() != 0)
                      {
                        std::array<double, 3> sim_mean_xyz{0., 0., 0.};
                        double size_mean = 0.;
                        for(const auto& SimH : it_ListHitsSim[id_det])
                          {
                            sim_mean_xyz[0] += 0.5 * (SimH.hitX + SimH.hitXexit);
                            sim_mean_xyz[1] += 0.5 * (SimH.hitY + SimH.hitYexit);
                            sim_mean_xyz[2] += 0.5 * (SimH.hitZ + SimH.hitZexit);
                            size_mean += 1.;
                          }

                        sim_mean_xyz[0] /= size_mean;
                        sim_mean_xyz[1] /= size_mean;
                        sim_mean_xyz[2] /= size_mean;

                        att._logger->debug("SImHit : det {} hit {}, [{}, {}, {}]",
                                           G4Sol::nameLiteralDet.begin()[id_det], id_hit, sim_mean_xyz[0],
                                           sim_mean_xyz[1], sim_mean_xyz[2]);
                      }
                  }

                // auto temp01_11 = new_miniFibers[id_mini][0].xyz[0]-new_miniFibers[id_mini][1].xyz[0] ;
                // auto temp01_12 = new_miniFibers[id_mini][0].xyz[0]-new_miniFibers[id_mini][1].xyz[1] ;

                // auto temp01_21 = new_miniFibers[id_mini][0].xyz[1]-new_miniFibers[id_mini][1].xyz[0] ;
                // auto temp01_22 = new_miniFibers[id_mini][0].xyz[1]-new_miniFibers[id_mini][1].xyz[1] ;

                auto [x12, y12] = Intersection(new_miniFibers[id_mini][0], new_miniFibers[id_mini][1]);

                // std::array<double,4> difIall = { temp01_11.Mag(), temp01_12.Mag(), temp01_21.Mag(), temp01_22.Mag()};
                // auto it_min = TMath::LocMin(difIall.begin(),difIall.end());
                // int id_min = std::distance(difIall.begin(),it_min);
                // if(id_min == 1 || id_min ==3)
                //   att._logger->error("E> ImproveHit : MiniFiber hit coincidence with different R on different layers
                //   (case size=2)| track {}", id_track);
                // else
		// int id_min_red = id_min/2;
		for(int i_xuv = 0; i_xuv < 2; ++i_xuv)
		  {
		    int id_det = new_miniFibers[id_mini][i_xuv].det;
		    int id_hit = new_miniFibers[id_mini][i_xuv].hit;

		    genfit::PlanarMeasurement* currentHit =
		      dynamic_cast<genfit::PlanarMeasurement*>(RecoEvent.ListHits[id_det][id_hit].get());
		    auto dummyState   = genfit::StateOnPlane();
		    auto currentPlane = currentHit->constructPlane(genfit::StateOnPlane());
		    auto HitCoord     = currentHit->getRawHitCoords();
		    auto HitCov       = currentHit->getRawHitCov();

		    TVectorD hitCoordsNew(2);
		    // hitCoordsNew(0) = new_miniFibers[id_mini][i_xuv].xyz[id_min_red].X();
		    // hitCoordsNew(1) = new_miniFibers[id_mini][i_xuv].xyz[id_min_red].Y();
		    hitCoordsNew(0) = x12;
		    hitCoordsNew(1) = y12;
		    TMatrixDSym hitCovNew(2);
		    hitCovNew(0, 0) = HitCov(0, 0);
		    hitCovNew(0, 1) = 0.;
		    hitCovNew(1, 0) = 0.;
		    hitCovNew(1, 1) = 1.2 * HitCov(0, 0);

		    TVector3 o(0, 0, new_miniFibers[id_mini][i_xuv].xyz[0].Z());
		    TVector3 u(1., 0., 0.), v(0., 1., 0.);
		    genfit::SharedPlanePtr newPlane(new genfit::DetPlane(o, u, v));

		    RecoEvent.OldListHits[id_det][id_hit] =
		      std::make_unique<genfit::PlanarMeasurement>(hitCoordsNew, hitCovNew, id_det, id_hit, nullptr);
		    dynamic_cast<genfit::PlanarMeasurement*>(RecoEvent.OldListHits[id_det][id_hit].get())
		      ->setPlane(newPlane);

		    RecoEvent.OldListHits[id_det][id_hit].swap(RecoEvent.ListHits[id_det][id_hit]);

		    if(it_ListHitsSim[id_det].size() != 0)
		      {
			std::array<double, 3> sim_mean_xyz{0., 0., 0.};
			double size_mean = 0.;
			for(const auto& SimH : it_ListHitsSim[id_det])
			  {
			    sim_mean_xyz[0] += 0.5 * (SimH.hitX + SimH.hitXexit);
			    sim_mean_xyz[1] += 0.5 * (SimH.hitY + SimH.hitYexit);
			    sim_mean_xyz[2] += 0.5 * (SimH.hitZ + SimH.hitZexit);
			    size_mean += 1.;
			  }

			sim_mean_xyz[0] /= size_mean;
			sim_mean_xyz[1] /= size_mean;
			sim_mean_xyz[2] /= size_mean;

			att._logger->debug("SImHit : det {} hit {}, [{}, {}, {}] | Est : [{}, {}, {}]",
					   G4Sol::nameLiteralDet.begin()[id_det], id_hit, sim_mean_xyz[0],
					   sim_mean_xyz[1], sim_mean_xyz[2], hitCoordsNew(0), hitCoordsNew(1),
					   new_miniFibers[id_mini][i_xuv].xyz[0].Z());

			TString nameTemp1(G4Sol::nameLiteralDet.begin()[id_det]);
			LocalHisto.h_ImproveHitsResidus->Fill(nameTemp1 + "_x", sim_mean_xyz[0] - hitCoordsNew(0),
							      1.);
			LocalHisto.h_ImproveHitsResidus->Fill(nameTemp1 + "_y", sim_mean_xyz[1] - hitCoordsNew(1),
							      1.);
			LocalHisto.h_ImproveHitsResidus->Fill(
							      nameTemp1 + "_z", sim_mean_xyz[2] - new_miniFibers[id_mini][i_xuv].xyz[0].Z(), 1.);
			LocalHisto.h_ImproveHitsResidus->Fill(
							      nameTemp1 + "_r",
							      TMath::Hypot(sim_mean_xyz[0] - hitCoordsNew(0), sim_mean_xyz[1] - hitCoordsNew(1)), 1.);
			LocalHisto.h_ImproveHitsResidus->Fill(
							      nameTemp1 + "_r",
							      TMath::Sqrt(TMath::Sq(sim_mean_xyz[0] - hitCoordsNew(0)) +
									  TMath::Sq(sim_mean_xyz[1] - hitCoordsNew(1)) +
									  TMath::Sq(sim_mean_xyz[2] - new_miniFibers[id_mini][i_xuv].xyz[0].Z())),
							      1.);

			// double Rest = it_paramFitRZ->second.extrapR(o.Z());

			// LocalHisto.h_ImproveHitsResidus->Fill(nameTemp1+"_fitRZ",Rest-TMath::Hypot(hitCoordsNew(0),hitCoordsNew(1)),1.);

			LocalHisto.h_ImproveHitsResidus2->Fill(nameTemp1 + "_x", sim_mean_xyz[0] - hitCoordsNew(0),
							       1.);
			LocalHisto.h_ImproveHitsResidus2->Fill(nameTemp1 + "_y", sim_mean_xyz[1] - hitCoordsNew(1),
							       1.);
			LocalHisto.h_ImproveHitsResidus2->Fill(
							       nameTemp1 + "_z", sim_mean_xyz[2] - new_miniFibers[id_mini][i_xuv].xyz[0].Z(), 1.);
			LocalHisto.h_ImproveHitsResidus2->Fill(
							       nameTemp1 + "_r",
							       TMath::Hypot(sim_mean_xyz[0] - hitCoordsNew(0), sim_mean_xyz[1] - hitCoordsNew(1)), 1.);
			LocalHisto.h_ImproveHitsResidus2->Fill(
							       nameTemp1 + "_r",
							       TMath::Sqrt(TMath::Sq(sim_mean_xyz[0] - hitCoordsNew(0)) +
									   TMath::Sq(sim_mean_xyz[1] - hitCoordsNew(1)) +
									   TMath::Sq(sim_mean_xyz[2] - new_miniFibers[id_mini][i_xuv].xyz[0].Z())),
							       1.);
			// LocalHisto.h_ImproveHitsResidus->Fill(nameTemp1+"_fitRZ",Rest-TMath::Hypot(hitCoordsNew(0),hitCoordsNew(1)),1.);
		      }
		  }
              }
        }
    }

  return 0;
}

template <class Out>
Improve::miniF TImproveHitsFiber<Out>::ImproveFiber2(const std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > >& ListHits, int id_det, int id_hit)
{

  genfit::PlanarMeasurement* currentHit = dynamic_cast<genfit::PlanarMeasurement*>(ListHits[id_det][id_hit].get());
  auto dummyState                       = genfit::StateOnPlane();
  auto currentPlane                     = currentHit->constructPlane(genfit::StateOnPlane());
  auto HitCoord                         = currentHit->getRawHitCoords();
  auto HitCov                           = currentHit->getRawHitCov();

  auto posPlane = currentPlane->getO();
  auto planeU   = currentPlane->getU();
  auto planeV   = currentPlane->getV();

  att._logger->debug("ImproveFiber : HitCoord : {} | O [{}, {}, {}] U [{}, {}, {}] V [{}, {}, {}]", HitCoord(0),
                     posPlane.X(), posPlane.Y(), posPlane.Z(), planeU.X(), planeU.Y(), planeU.Z(), planeV.X(),
                     planeV.Y(), planeV.Z());

  double res = std::sqrt(HitCov(0, 0));

  TVector3 posXYZ_1 = currentPlane->toLab({HitCoord(0) + res, 0});
  TVector3 posXYZ_2 = currentPlane->toLab({HitCoord(0) - res, 0});

  return {id_det, id_hit, {posXYZ_1, planeV}, {posXYZ_2, planeV}, res};
};

template <class Out>
std::tuple<double, double> TImproveHitsFiber<Out>::Intersection(const Improve::miniF& fiber1,
                                                                const Improve::miniF& fiber2)
{
  std::vector<std::tuple<double, double> > intersections;
  intersections.push_back(Inter(fiber1.xyz[0], fiber1.xyz[1], fiber2.xyz[0], fiber2.xyz[1]));
  intersections.push_back(Inter(fiber1.xyz_low[0], fiber1.xyz_low[1], fiber2.xyz[0], fiber2.xyz[1]));
  intersections.push_back(Inter(fiber1.xyz[0], fiber1.xyz[1], fiber2.xyz[0], fiber2.xyz_low[1]));
  intersections.push_back(Inter(fiber1.xyz_low[0], fiber1.xyz_low[1], fiber2.xyz[0], fiber2.xyz_low[1]));

  double mean_x = 0., mean_y = 0.;
  for(auto [x1, y1] : intersections)
    {
      mean_x += x1;
      mean_y += y1;
    }
  mean_x /= 4.;
  mean_y /= 4.;

  return std::make_tuple(mean_x, mean_y);
};

template <class Out>
std::tuple<double, double> TImproveHitsFiber<Out>::Inter(const TVector3& p1, const TVector3& d1, const TVector3& p2,
                                                         const TVector3& d2)
{
  double cross  = d2.Y() * d1.X() - d2.X() * d1.Y();
  double cross2 = TMath::Min(1. / cross, 10000000.);
  double t1     = (p2.X() - p1.X()) * d2.Y() - (p2.Y() - p1.Y()) * d2.X();
  t1 *= cross2;
  double s1 = (p2.X() - p1.X()) * d1.Y() - (p2.Y() - p1.Y()) * d1.X();
  s1 *= cross2;

  att._logger->debug("inter : t1 {} s1 {} | [{} {}] [{} {}]", t1, s1, p1.X() + t1 * d1.X(), p1.Y() + t1 * d1.Y(),
                     p2.X() + s1 * d2.X(), p2.Y() + s1 * d2.Y());

  return std::make_tuple(p1.X() + t1 * d1.X(), p1.Y() + t1 * d1.Y());
}

template class TImproveHitsFiber<MCAnaEventG4Sol>;
template class TImproveHitsFiber<Ana_WasaEvent>;
