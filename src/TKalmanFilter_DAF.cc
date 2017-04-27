#include "TKalmanFilter_DAF.h"

#include "TRandom3.h"

#include "KalmanFittedStateOnPlane.h"
#include "KalmanFitterInfo.h"
#include "StateOnPlane.h"

#include "Exception.h"

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"

#include "TGeoBBox.h"
#include "TGeoNode.h"
#include "TObjArray.h"

#include <algorithm>
#include <memory>
#include <numeric>

//#define DEBUG_KALMAN2
//#define DEBUG_KALMAN
//#define ROTATION_KALMAN

//#define DISPLAY

using namespace std;
using namespace G4Sol;

TKalmanFilter_DAF::TKalmanFilter_DAF(const THyphiAttributes& attribut) : TDataProcessInterface("mom_fit_kalman"), att(attribut)
{

  const int nIter = 7;        // max number of iterations
  const double dPVal = 1.e-3; // convergence criterion

  // const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedAverage;
  // For Ref
  // const genfit::eMultipleMeasurementHandling mmHandling = genfit::unweightedClosestToReference;
  // const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedClosestToReference;
  // const genfit::eMultipleMeasurementHandling mmHandling = genfit::unweightedClosestToReferenceWire;
  // const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedClosestToReferenceWire;

  // const genfit::eMultipleMeasurementHandling mmHandling = genfit::unweightedClosestToPrediction;
  // const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedClosestToPrediction;
  const genfit::eMultipleMeasurementHandling mmHandling = genfit::unweightedClosestToPredictionWire;
  // const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedClosestToPredictionWire;

  Fitter_rescue = new genfit::KalmanFitterRefTrack(nIter, dPVal); //,1e4);
  Fitter_rescue->setMultipleMeasurementHandling(mmHandling);
  Fitter_rescue->setMinIterations(2);
  Fitter_rescue->setMaxIterations(nIter);

  // Fitter = new genfit::KalmanFitter(6,dPVal,1e-3,true);//
  // Fitter = new genfit::DAF(true);
  Fitter = new genfit::KalmanFitterRefTrack(nIter, dPVal); //,1e4);
  Fitter->setMinIterations(2);
  Fitter->setMultipleMeasurementHandling(mmHandling);
  Fitter->setMaxIterations(nIter);

  // Fitter->setDebugLvl(10);

  rep = new genfit::RKTrackRep();

  Vtracks = new genfit::Track(rep, TVector3(0., 0., 0.), TVector3(0., 0., 1.0));

#ifdef DISPLAY
  display = genfit::EventDisplay::getInstance();
  display->reset();
  std::cout << " Display On :" << display << std::endl;
#endif
}

TKalmanFilter_DAF::~TKalmanFilter_DAF()
{
  if(Vtracks != NULL)
    {
      delete Vtracks;
      Vtracks = 0;
    }

  if(Fitter != NULL)
    {
      delete Fitter;
      Fitter = 0;
    }
  if(Fitter_rescue != NULL)
    {
      delete Fitter_rescue;
      Fitter_rescue = 0;
    }

  for(auto& elem : list_Plane)
    if(elem != 0)
      {
        delete elem;
        elem = 0;
      }
}

int TKalmanFilter_DAF::operator()(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
{

  int result_mom = Exec(RecoEvent, OutTree);

  return SoftExit(result_mom);
}

int TKalmanFilter_DAF::Exec(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
{

  //  int result_kalman = Kalman_Filter(RecoEvent);
  int result_kalman = Kalman_Filter_FromTrack(RecoEvent);

  for(auto TrackRes : RecoEvent.DAF_results)
    {
      int TrackID = TrackRes.first;
      const ResSolDAF& FitRes = TrackRes.second;
      auto TInfo = RecoEvent.TrackInfo.find(TrackID);
      auto IsDecay = RecoEvent.TrackMother.find(TrackID);
      int Decay = IsDecay == RecoEvent.TrackMother.end() ? 0 : 1;

      THyphiTrack* OutTrack = dynamic_cast<THyphiTrack*>(OutTree->fTrack->ConstructedAt(OutTree->fTrack->GetEntries()));
      auto PDG_particle = TDatabasePDG::Instance()->GetParticle(FitRes.pdg_guess);
      OutTrack->type = PDG_particle->GetName();
      OutTrack->MC_status = TrackID + 10000 * Decay;
      OutTrack->Chi2 = FitRes.chi2;
      OutTrack->Chi2_X = FitRes.ndf;
      OutTrack->Chi2_Y = FitRes.firstHit;
      OutTrack->Mass = FitRes.mass;
      OutTrack->pdgcode = FitRes.pdg_guess;
      OutTrack->MomMass.SetXYZM(FitRes.momX, FitRes.momY, FitRes.momZ, FitRes.mass);
      OutTrack->Mom.SetXYZ(FitRes.momX, FitRes.momY, FitRes.momZ);

      OutTrack->BarId = FitRes.lastHit;
      OutTrack->Charge = FitRes.charge;
      OutTrack->dE = TInfo->second[FitRes.lastHit].Eloss;
      OutTrack->Beta = FitRes.beta;
      OutTrack->RefPoint.SetXYZ(FitRes.posX, FitRes.posY, FitRes.posZ);
      OutTrack->Pval2 = FitRes.pvalue;
      OutTrack->PathLength = FitRes.path_length;
      OutTrack->TOF = FitRes.tof;

      OutTrack->MomIni.SetXYZ(FitRes.momX_init, FitRes.momY_init, FitRes.momZ_init);
      OutTrack->BetaIni = FitRes.beta2;
      OutTrack->MassIni = FitRes.mass2;
      OutTrack->TOFIni = FitRes.tof2;
      OutTrack->PathLengthIni = TInfo->second[FitRes.lastHit].length; // FitRes.path_length2;
      OutTrack->RChiIni = FitRes.fitter;
      if(Decay == 1)
        OutTrack->Sim2Vtx.SetXYZT(std::get<1>(IsDecay->second), std::get<2>(IsDecay->second), std::get<3>(IsDecay->second),
                                  std::get<4>(IsDecay->second));

      OutTrack->State[0] = FitRes.posX;
      OutTrack->State[1] = FitRes.posY;
      OutTrack->State[2] = FitRes.posZ;
      OutTrack->State[3] = FitRes.momX;
      OutTrack->State[4] = FitRes.momY;
      OutTrack->State[5] = FitRes.momZ;
      for(int row = 0; row < 6; row++)
        for(int col = 0; col < 6; col++)
          OutTrack->Cov[row][col] = FitRes.cov_matrix[row][col];
    }
  OutTree->Ntrack = OutTree->fTrack->GetEntries();

  // #ifdef DEBUG_DAISUKE
  //       std::cout<<"HitTR1 ("<<it_mom_pv->second[13]<<","<<it_mom_pv->second[14]<<","<<it_mom_pv->second[15]<<")"<<std::endl;
  //       std::cout<<"HitTR2 ("<<it_mom_pv->second[16]<<","<<it_mom_pv->second[17]<<","<<it_mom_pv->second[18]<<")"<<std::endl;
  //       std::cout<<"HitDC2 ("<<it_mom_pv->second[19]<<","<<it_mom_pv->second[20]<<","<<it_mom_pv->second[21]<<")"<<std::endl;
  //       std::cout<<"NumDC2 ("<<it_mom_pv->second[22]<<","<<it_mom_pv->second[23]<<","<<it_mom_pv->second[24]<<")"<<std::endl;

  //       std::cout<<"Pval2 = "<<it_mom_pv->second[27]<<std::endl;
  //       std::cout<<"TofsBar = "<<it_mom_pv->second[28]<<std::endl;
  //       std::cout<<"BarId = "<<it_mom_pv->second[25]<<std::endl;
  //       std::cout<<"Charge = "<<it_mom_pv->second[6]<<std::endl;
  //       std::cout<<"Beta = "<<it_mom_pv->second[26]<<std::endl;
  // #endif
  //     }

  return result_kalman;
}

int TKalmanFilter_DAF::SoftExit(int result_full) { return result_full; }

int TKalmanFilter_DAF::Kalman_Filter_FromTrack(FullRecoEvent& RecoEvent)
{
  auto printW = [](const auto a, const int width) -> std::string {
    std::stringstream ss;
    ss << std::fixed << std::right;
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

  unsigned int Nb_trackSeed = RecoEvent.TrackDAF.size();

  int ntrack = -1;
  for(auto it_trackInfo : RecoEvent.TrackInfo)
    {
      ntrack++;

#ifdef DEBUG_KALMAN
      std::cout << "start #" << ntrack; //<<std::endl;
// std::cout<<" | Pointsize "<<Vtracks->getNumPoints()<<" Measurements "<<Vtracks->getNumPointsWithMeasurement()<<" Rep
// "<<Vtracks->getNumReps()<<std::endl;
// std::vector<std::string> name_det;
#endif

      const int id_track = it_trackInfo.first;
      auto it_ListHits = RecoEvent.TrackDAF.find(id_track);

      std::set<std::tuple<double, int, int> > id_dets;

      for(size_t id_det = 0; id_det < it_ListHits->second.size(); ++id_det)
        {
          int id_hit = it_ListHits->second[id_det];
          if(id_hit < 0)
            continue;

          genfit::AbsMeasurement* currentHit = RecoEvent.ListHits[id_det][id_hit].get();
          id_dets.insert(std::make_tuple(currentHit->getRawHitCoords()(2), id_det, id_hit));
        }

      // if(id_firstDet != std::get<1>(*firstHit) )
      // 	{
      // 	  std::cout<<"!> Kalman: firstHit not matched !"<<id_firstDet<<" "<<std::get<1>(*firstHit)<<"\n";
      // 	  auto Decay = RecoEvent.TrackMother.find(id_track);
      // 	  bool IsDecay = Decay != RecoEvent.TrackMother.end();
      // 	  std::cout<<"track #"<<id_track<<" "<<PDG<<" "<<IsDecay<<"\n";
      // 	  for(auto ii : id_dets)
      // 	    std::cout<<"| :"<<std::get<0>(ii)<<" "<<std::get<1>(ii)<<" "<<std::get<2>(ii)<<"\n";
      // 	  std::cout<<" +---\n";
      // 	}

      // for(size_t idTempDet = 0; idTempDet < it_ListHits->second.size(); ++idTempDet)
      //   if(it_ListHits->second[idTempDet] != -1)
      //     {
      //       if(id_firstDet == -1 && idTempDet != G4Sol::PSFE)
      //         id_firstDet = idTempDet;

      //       ++nb_ValidHits;
      //       // break;
      //     }

      // if(id_firstDet == -1)
      //   continue;

      if(id_dets.size() < 3) //(nb_ValidHits < 3)
        continue;

      auto f_LastHitIsValid = [](const auto& it_ListHits, std::set<G4Sol::SolDet> listToTest) {

        for(auto it_det : listToTest)
          if(it_ListHits->second[it_det] >= 0)
            return it_det;

        return G4Sol::SIZEOF_G4SOLDETTYPE;

      };
      G4Sol::SolDet lastValidHit = f_LastHitIsValid(it_ListHits, {G4Sol::PSCE, G4Sol::PSFE, G4Sol::RPC_l, G4Sol::RPC_h});

      if(lastValidHit == G4Sol::SIZEOF_G4SOLDETTYPE)
        continue;

      const auto firstHit = id_dets.cbegin();
      const auto lastHit = id_dets.crbegin();

      int id_firstDet = std::get<1>(*firstHit);
      int id_lastHit = std::get<1>(*lastHit); // std::get<0>(lastHit);

      if(TMath::Abs(it_trackInfo.second[id_firstDet].pdg) < 1e-2)
        continue;

#ifdef DEBUG_KALMAN
      std::cout << " Id_track: " << id_track << std::endl;
#endif

      const InfoPar track_state(it_trackInfo.second[id_firstDet]);

      const int PDG = static_cast<int>(track_state.pdg);
      auto PDG_particle = TDatabasePDG::Instance()->GetParticle(PDG);
      if(PDG_particle == nullptr)
        {
          std::cout << "E> PDG not found !" << std::endl;
          continue;
        }
      const double charge = PDG_particle->Charge() / 3.;

      const int PID = charge * charge;

#ifdef DEBUG_KALMAN
      std::cout << " PID:" << PID << endl;
#endif

      const TVector3 init_p(track_state.momX, track_state.momY, track_state.momZ);

      const double seed_Mom_Mag = init_p.Mag();
      if(TMath::Abs(seed_Mom_Mag) < 1e-9)
        {
          std::cout << "E> Seed Momemtum with TVector3 is zero ! " << PDG_particle->GetName() << " Mom:" << seed_Mom_Mag << "\n";
          auto printW = [](const auto a, const int width) -> std::string {
            std::stringstream ss;
            ss << std::fixed << std::right;
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

          std::cout << "TrackID #" << it_ListHits->first << " hit_id :\n";
          std::vector<std::stringstream> s1(it_ListHits->second.size() / 20 + 1);
          std::vector<std::stringstream> s2(it_ListHits->second.size() / 20 + 1);
          for(size_t i = 0; i < it_ListHits->second.size(); ++i)
            {
              s1[i / 20] << printW(i, 3) << ", ";
              s2[i / 20] << printW(it_ListHits->second[i], 3) << ", ";
            }
          for(size_t i = 0; i < s1.size(); ++i)
            {
              std::cout << "idDet:" << s1[i].str() << "\n";
              std::cout << "stat :" << s2[i].str() << "\n";
            }

          std::cout << "Track Info :\n";
          std::vector<std::stringstream> s11(it_ListHits->second.size() / 10 + 1);
          std::vector<std::stringstream> s22(it_ListHits->second.size() / 10 + 1);
          std::vector<std::stringstream> s33(it_ListHits->second.size() / 10 + 1);
          std::vector<std::stringstream> s44(it_ListHits->second.size() / 10 + 1);
          std::vector<std::stringstream> s55(it_ListHits->second.size() / 10 + 1);
          std::vector<std::stringstream> s66(it_ListHits->second.size() / 10 + 1);
          for(size_t i = 0; i < it_trackInfo.second.size(); ++i)
            {
              s11[i / 10] << printW(i, 9) << ", ";
              s22[i / 10] << printFixed(it_trackInfo.second[i].pdg, 3, 9) << ", ";
              s33[i / 10] << printFixed(it_trackInfo.second[i].momX, 3, 9) << ", ";
              s44[i / 10] << printFixed(it_trackInfo.second[i].momY, 3, 9) << ", ";
              s55[i / 10] << printFixed(it_trackInfo.second[i].momZ, 3, 9) << ", ";
              s66[i / 10] << printFixed(it_trackInfo.second[i].time, 3, 9) << ", ";
            }
          for(size_t i = 0; i < s1.size(); ++i)
            {
              std::cout << "idDet:" << s11[i].str() << "\n";
              std::cout << "pdg  :" << s22[i].str() << "\n";
              std::cout << "momX :" << s33[i].str() << "\n";
              std::cout << "momY :" << s44[i].str() << "\n";
              std::cout << "momZ :" << s55[i].str() << "\n";
              std::cout << "time :" << s66[i].str() << "\n";
            }

          // for(auto id_hit : track.second)
          // 	std::cout<<" "<<id_hit<<", ";
          // std::cout<<"] "<<std::endl;

          continue;
        }

#ifdef DEBUG_KALMAN
      std::cout << "init mom : ";
      init_p.Print();
#endif

      // const double mom_res = .200; // GeV
      // double new_P = gRandom->Gaus(seed_Mom_Mag,mom_res);
      TVector3 seed_p(init_p);
      // seed_p.SetMag(new_P);

      genfit::AbsTrackRep* REP = new genfit::RKTrackRep(PDG);

      genfit::AbsMeasurement* FirstHit = RecoEvent.ListHits[id_firstDet][it_ListHits->second[id_firstDet]].get();
      TVectorD& HitrawRef = FirstHit->getRawHitCoords();
      const TVector3 init_point(HitrawRef[0], HitrawRef[1], HitrawRef[2]);

      // TVector3 init_point(track_state[9],track_state[10],track_state[11]);

      TMatrixDSym CovM(6);
      // double resolution = 0.1;
      for(int i = 0; i < 6; ++i)
        CovM(i, i) = 1000.;

      genfit::MeasuredStateOnPlane stateRef(REP);
      // TEP->setPDG(infoTrack.PDG);
      REP->setPosMomCov(stateRef, init_point, seed_p, CovM);
      // remember original initial state
      const genfit::StateOnPlane stateRefOrig(stateRef.getState(), stateRef.getPlane(), REP, stateRef.getAuxInfo());

      // create track
      TVectorD seedState(6);
      TMatrixDSym seedCov(6);
      REP->get6DStateCov(stateRef, seedState, seedCov);

#ifdef DEBUG_KALMAN
      std::cout << "init_p:" << std::endl;
      init_p.Print();
      init_p.Unit().Print();
      std::cout << "seed_p:" << std::endl;
      seed_p.Print();
      std::cout << "init_point :" << std::endl;
      init_point.Print();
      std::cout << " SEED :" << std::endl;
      seedState.Print();
      std::cout << "--" << std::endl;
      seedCov.Print();
// std::cout<<"init_plane :"<<std::endl;
// init_plane.Print();
#endif

      if(REP == NULL)
        std::cout << "E> no Rep build" << std::endl;
#ifdef DEBUG_KALMAN
      else
        std::cout << "rep done" << std::endl;
#endif

// std::auto_ptr<genfit::Track> fitTrack(new genfit::Track(REP, seedState, seedCov));
#ifndef DISPLAY
      std::unique_ptr<genfit::Track> fitTrack = std::make_unique<genfit::Track>(REP, seedState, seedCov);
#else
      std::shared_ptr<genfit::Track> fitTrack = std::make_shared<genfit::Track>(REP, seedState, seedCov);
#endif

      // std::tuple<int, double> lastHit = std::make_tuple(-1, -9999.);

      for(auto ids : id_dets)
        {
          int id_det = std::get<1>(ids);
          int id_hit = std::get<2>(ids);

          genfit::AbsMeasurement* currentHit = RecoEvent.ListHits[id_det][id_hit].get();
#ifdef DEBUG_KALMAN
          std::cout << "Loop insertPoint: #det:" << id_det << " #hit:" << id_hit << " " << currentHit << "\n";
          currentHit->Print();
#endif
          // if(std::get<1>(lastHit) < currentHit->getRawHitCoords()(2))
          //   {
          //     lastHit = std::make_tuple(id_det, currentHit->getRawHitCoords()(2));
          //   }

          // Vtracks->insertPoint(new genfit::TrackPoint(RecoEvent.ListHitsDAF[id_track][ii],Vtracks));
          fitTrack->insertPoint(new genfit::TrackPoint(currentHit, fitTrack.get(), false));
          // id_lastHit = id_det;
        }

      InfoPar track_stateLast(it_trackInfo.second[id_lastHit]);

      // assert(fitTrack->checkConsistency());
      fitTrack->checkConsistency();
#ifdef DEBUG_KALMAN2
      std::cout << "track n'" << ntrack << std::endl;
#endif
      int whoDoFit = 1;
      try
        {
          Fitter->processTrack(fitTrack.get(), true);
        }
      catch(genfit::Exception& e1)
        {

          std::cout << "*** FITTER EXCEPTION *** Rescue fitter take place !" << std::endl;
          std::cout << e1.what() << std::endl;
          AnaHisto->h_stats->Fill("Exc:RescueFit", 1.);

          try
            {
              whoDoFit = 2;
              Fitter_rescue->processTrack(fitTrack.get(), true);
            }
          catch(genfit::Exception& e2)
            {
              std::cout << "*** FITTER Rescue EXCEPTION *** cancel fit !" << std::endl;
              std::cout << e2.what() << std::endl;
              AnaHisto->h_stats->Fill("Exc:FitFail", 1.);

              continue;
            }
        }

#ifdef DISPLAY
      display->addEvent(fitTrack.get()); // Vtracks);
#endif

#ifdef DEBUG_KALMAN2
      std::cout << "SUCESSFULL FIT!" << std::endl;
      std::cout << "track n'" << ntrack << " / " << fitTrack->getFitStatus(REP)->isFitConverged() << std::endl;
      fitTrack->Print();
#endif

      if(fitTrack->getFitStatus(REP)->isFitConverged())
        {
          assert(fitTrack->hasKalmanFitStatus() == true);
          // assert(fitTrack->checkConsistency());
          fitTrack->checkConsistency();
          genfit::TrackPoint* tp = fitTrack->getPointWithMeasurementAndFitterInfo(0, REP);
          if(tp == NULL)
            {
              std::cout << "Track has no TrackPoint with fitterInfo! \n";
              AnaHisto->h_stats->Fill("Exc:NoTrackPointInfo", 1.);
              continue;
            }

          genfit::KalmanFittedStateOnPlane kfsop(*(dynamic_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(REP))->getBackwardUpdate()));

          try
            {
              REP->extrapolateToPlane(kfsop, stateRefOrig.getPlane());
            }
          catch(genfit::Exception& e)
            {
              std::cerr << "Exception, next track" << std::endl;
              AnaHisto->h_stats->Fill("Exc:FitBacktoOrig", 1.);
              std::cerr << e.what();
              continue;
            }

          ResSolDAF tempResults;
          tempResults.charge = charge;
          tempResults.pdg_ini = PDG;
          tempResults.pdg_guess = PDG;

          const TVectorD& referenceState = stateRefOrig.getState();
          const TVectorD& stateFit = kfsop.getState();
          TVector3 p3, posRef;
          TMatrixDSym covFit(6);

          kfsop.getPosMomCov(posRef, p3, covFit);

          tempResults.momX_init = init_p.X();
          tempResults.momY_init = init_p.Y();
          tempResults.momZ_init = init_p.Z();

          tempResults.momX = p3.X();
          tempResults.momY = p3.Y();
          tempResults.momZ = p3.Z();

          tempResults.posX = posRef.X();
          tempResults.posY = posRef.Y();
          tempResults.posZ = posRef.Z();

          for(int row = 0; row < 6; ++row)
            for(int col = 0; col < 6; ++col)
              tempResults.cov_matrix[row][col] = covFit(row, col);

          const double chi2 = whoDoFit == 1 ? Fitter->getChiSqu(fitTrack.get(), REP)
                                            : Fitter_rescue->getChiSqu(fitTrack.get(), REP); // Vtracks->getChiSqu();

          const double ndf = whoDoFit == 1 ? Fitter->getNdf(fitTrack.get(), REP)
                                           : Fitter_rescue->getNdf(fitTrack.get(), REP); // static_cast<int>(Vtracks->getNDF());

          tempResults.chi2 = chi2;
          tempResults.ndf = ndf;
          tempResults.fitter = whoDoFit;
#ifdef DEBUG_KALMAN2
          std::cout << " / chi2 =" << chi2 << " / ndf =" << ndf << std::endl;
#endif

          // double p_value = TMath::Prob(chi2,ndf);
          // double p_value = MathKalman::Prob()(chi2,ndf);
          // double p_value2= 1.-p_value;
          const double p_value =
              whoDoFit == 1 ? Fitter->getPVal(fitTrack.get(), REP) : Fitter_rescue->getPVal(fitTrack.get(), REP); // TMath::Prob(chi2,ndf);
          const double p_value2 = 1. - p_value;

          tempResults.pvalue = p_value;

          AnaHisto->h_pv->Fill(p_value2);
          AnaHisto->h_chi2->Fill(chi2);

// h_chi2_smooth->Fill(Vtracks[itr]->getChiSquSmooth());

#ifdef DEBUG_KALMAN2
          std::cout << " / 1-p_value =" << 1 - p_value << std::endl;
#endif

          // n_status[static_cast<int>(charge)+1]++;

          if(charge > 0)
            {
              AnaHisto->hd_chi[0]->Fill(chi2 / ndf);
              AnaHisto->hd_pv[0]->Fill(p_value2);
            }
          else
            {
              AnaHisto->hd_chi[1]->Fill(chi2 / ndf);
              AnaHisto->hd_pv[1]->Fill(p_value2);
            }

          // GFDetPlane InitPlane(RecoEvent.ListHitsDAF[id_track][TR1X][0]->getDetPlane(rep));
          // TVector3 TR1x_Pos(Vtracks->getPos(InitPlane));

          // double DistToTR1= 0;
          // try
          //   {
          //     const genfit::MeasuredStateOnPlane  TR1PlaneState(fitTrack->getFittedState(find_tr1->second,REP,true));
          //     TVector3 TR1x_Pos(TR1PlaneState.getPos());
          //     HitTR1 = TR1x_Pos;
          //     DistToTR1=(TR1x_Pos-Plane_time).Mag();
          //   }
          // catch(genfit::Exception& e)
          //   {
          //     std::cerr<<"Exception in TR1PlaneState, next track"<<std::endl;
          //     std::cerr << e.what();
          //     continue;
          //   }

          const double p = p3.Mag();

          // if(att.G4_simu==false)
          // DistToTR1+=5.; // distance from Start to TR1

          double Path_length = -999;
          double Path_time = -999;
          try
            {
              Path_length = fitTrack->getTrackLen(REP); // getTotalLength();
              // double Path_length = Vtracks->getTotalLength();
              // double Path_lengthB = Vtracks->getTotalLengthBack();
            }
          catch(genfit::Exception& e)
            {
              std::cerr << e.what();
              std::cout << "could not get TrackLen ! \n";
              continue;
            }
          try
            {
              Path_time = fitTrack->getTOF(REP);
              // Path_length = fitTrack->getTrackLen(REP);//getTotalLength();
              // double Path_length = Vtracks->getTotalLength();
              // double Path_lengthB = Vtracks->getTotalLengthBack();
            }
          catch(genfit::Exception& e)
            {
              std::cerr << e.what();
              std::cout << "could not get TOF! \n";
              AnaHisto->h_stats->Fill("Exc:TOF", 1.);
              continue;
            }

          // double Path_timeB = Vtracks->getTotalTimeBack();
          // double Path_lengthMean = Path_length + Path_lengthB;
          double Path_lengthMean = Path_length;
          // Path_lengthMean/=2.;
          Path_lengthMean += init_point.Mag(); // DistToTR1;
                                               // 	  if(ndf==1)
                                               // 	    Path_lengthMean=Path_length;

          AnaHisto->h_Path->Fill(Path_length);
          // AnaHisto->h_Path_Back->Fill(Path_lengthB);
          AnaHisto->h_MeanPath->Fill(Path_lengthMean);
// AnaHisto->h_dpath->Fill(Path_lengthB-Path_length);

// double rapidity = 0.5 * TMath::Log(( TMath::Sqrt(rep_length->getMass()*rep_length->getMass() + p3.Mag2()) +p3.Z())/(
// TMath::Sqrt(rep_length->getMass()*rep_length->getMass() + p3.Mag2()) - p3.Z()));
// for(unsigned int vol = 0 ; vol < att.name_GeoVolumes.size() ; ++vol)
//   {
//     AnaHisto->Material_XX0_y[vol]->Fill(rapidity,rep_length->getXX0(att.name_GeoVolumes[vol]));
//     AnaHisto->Material_dE_y[vol]->Fill(rapidity,rep_length->getDE(att.name_GeoVolumes[vol]));
//   }

#ifdef VERBOSE_EVE

          std::cout << "$$$$$$$$$$$$-----------------------$$$$$$$$$$$$$$$" << std::endl;
          std::cout << "DAF Length :ListTrackPoints :" << std::endl;
          std::cout << "TEveLine* EveTrack_DAF1 = new TEveLine(TrackRepList_DAF," << rep_length->EveListTrackPoint.size() << ");"
                    << std::endl;
          for(unsigned int i = 0; i < rep_length->EveListTrackPoint.size(); ++i)
            {
              std::cout << std::setprecision(10) << "EveTrack_DAF1->SetPoint(" << i << "," << rep_length->EveListTrackPoint[i].X() << ","
                        << rep_length->EveListTrackPoint[i].Y() << "," << rep_length->EveListTrackPoint[i].Z() << ");" << std::endl;
            }

          std::cout << "Fitted Hit :" << std::endl;
          std::cout << "TEvePointSet* EveHitFitter1 = new TEvePointSet(HitFitter1," << 5 << ");" << std::endl;
          std::cout << std::setprecision(10) << "EveHitFitter1->SetPoint(" << 0 << "," << HitTR1.X() << "," << HitTR1.Y() << ","
                    << HitTR1.Z() << ");" << std::endl;
          std::cout << std::setprecision(10) << "EveHitFitter1->SetPoint(" << 1 << "," << HitDC1.X() << "," << HitDC1.Y() << ","
                    << HitDC1.Z() << ");" << std::endl;
          std::cout << std::setprecision(10) << "EveHitFitter1->SetPoint(" << 2 << "," << HitTR2.X() << "," << HitTR2.Y() << ","
                    << HitTR2.Z() << ");" << std::endl;
          std::cout << std::setprecision(10) << "EveHitFitter1->SetPoint(" << 3 << "," << HitDC2.X() << "," << HitDC2.Y() << ","
                    << HitDC2.Z() << ");" << std::endl;
          std::cout << std::setprecision(10) << "EveHitFitter1->SetPoint(" << 4 << "," << HitTOF.X() << "," << HitTOF.Y() << ","
                    << HitTOF.Z() << ");" << std::endl;

          std::cout << "DAF Track :" << std::endl;
          std::cout << "TEveRecTrack *DAF_track_rc = new TEveRecTrack();" << std::endl;
          std::cout << "DAF_track_rc->fV.Set(" << HitTR1.X() << "," << HitTR1.Y() << "," << HitTR1.Z() << ");" << std::endl;
          std::cout << "DAF_track_rc->fP.Set(" << p3.X() << "," << p3.Y() << "," << p3.Z() << ");" << std::endl;
          // std::cout<<"int sign = -"<<list_MC[current_idMC]->Charge<<std::endl;
          std::cout << "DAF_track_rc->fSign = -" << charge << std::endl;
          std::cout << "$$$$$$$$$$$$-----------------------$$$$$$$$$$$$$$$" << std::endl;

#endif

          const double dE = track_stateLast.Eloss;
          const double time_of_flight = track_stateLast.time; // Vtracks->getTimeOfFlight();//+6.;
          const double time_fromFirstToLastHit = track_stateLast.time - track_state.time;

          const double beta =
              1. / 29.9792458 * Path_lengthMean / time_of_flight; /// beta = dL[m]/(c[m/s]* dT[s]) = dL[mm]/(300[mm/ns] * dT[ns])
          const double beta_FirstToLast = 1. / 29.9792458 * Path_length / time_fromFirstToLastHit;
          // double beta = 1./30.*new_Total_Length/time_of_flight; ///beta = dL[m]/(c[m/s]* dT[s]) = dL[mm]/(300[mm/ns] * dT[ns])
          // double gamma = 1./sqrt(1.-beta*beta);
          const double mass = p * sqrt(1. / (beta * beta) - 1.);
          const double mass2 = p * sqrt(1. / (beta_FirstToLast * beta_FirstToLast) - 1.);

          if(beta < 0.)
            {
              cout << "E> Beta Negative ! " << Path_lengthMean << " " << time_of_flight << " | " << Path_length << " "
                   << endl; //<<Path_lengthB<<" "<<endl;
              cout << charge << " " << p << " " << chi2 << " " << ndf << " " << p_value2 << endl;
            }

#ifdef DEBUG_KALMAN2
          std::cout << "charge: " << charge
                    << " / time = " << time_of_flight /*Detectors_UTr["TOF"][id_track][3]*/ << " / Path = " << Path_length
                    << " PathMean = " << Path_lengthMean << " / Beta = " << beta << " / Mass = " << mass << std::endl;
#endif

          AnaHisto->h_beta->Fill(beta);
          AnaHisto->h_Mass_All->Fill(mass);
          AnaHisto->h_Mass_charge_All->Fill(mass, charge);
          AnaHisto->h_beta_mom->Fill(p, beta);
          AnaHisto->h_path_tof->Fill(Path_lengthMean / 30., time_of_flight);

          AnaHisto->h_pv_mom->Fill(p, p_value2);
          AnaHisto->h_pv_beta->Fill(beta, p_value2);
          AnaHisto->h_pv_mass->Fill(mass, p_value2);

          const std::vector<int> hist_to_pdg = {2212, -211, 211, -321, 321};
          for(size_t it_pdg = 0; it_pdg < hist_to_pdg.size(); ++it_pdg)
            {
              if(hist_to_pdg[it_pdg] == PDG)
                AnaHisto->h_mom_res[it_pdg]->Fill(init_p.Mag(), std::fabs(init_p.Mag() - p3.Mag()) / init_p.Mag());
            }

          if(p_value2 < .75)
            {
              AnaHisto->h_beta2->Fill(beta);
              AnaHisto->h_Mass_All2->Fill(mass);
              AnaHisto->h_beta_mom2->Fill(p, beta);
              AnaHisto->h_Mass_charge_All2->Fill(mass, charge);
            }
          if(p_value2 < .1)
            {
              AnaHisto->h_beta3->Fill(beta);
              AnaHisto->h_Mass_All3->Fill(mass);
              AnaHisto->h_beta_mom3->Fill(p, beta);
              AnaHisto->h_Mass_charge_All3->Fill(mass, charge);
              AnaHisto->h_mom_tof_cut->Fill(p, time_of_flight);
              AnaHisto->h_path_mom_cut->Fill(Path_lengthMean / 30., p);
              AnaHisto->h_path_tof_cut->Fill(Path_lengthMean / 30., time_of_flight);
            }
          tempResults.firstHit = id_firstDet;
          tempResults.lastHit = id_lastHit;
          tempResults.dE = dE;
          tempResults.path_time = Path_time;
          tempResults.path_length = Path_lengthMean;
          tempResults.path_length2 = Path_length;
          tempResults.tof = time_of_flight;
          tempResults.tof2 = time_fromFirstToLastHit;
          tempResults.beta = beta;
          tempResults.beta2 = beta_FirstToLast;
          tempResults.mass = mass;
          tempResults.mass2 = mass2;

          double m_range[4] = {0.9383, 3.72738, 0.1396, 2.809};
          double m_charge[4] = {1., 2., -1., 2.};

          for(int i = 0; i < 4; i++)
            // if (TMath::Abs(mass-m_range[i])<0.25*m_range[i] && TMath::Abs(charge-m_charge[i])<0.1)
            if(TMath::Abs(charge - m_charge[i]) < 0.1)
              {
                // double dmass = TMath::Abs(mass-m_range[i])/m_range[i];
                AnaHisto->h_Mass[i]->Fill(mass);
                AnaHisto->h_chi2_particle[i]->Fill(chi2);
                AnaHisto->h_pv_particle[i]->Fill(p_value2);
              }

#ifdef DEBUG_KALMAN
          std::cout << "Temp :" << ntrack << " " << p_value2 << " " << mass << std::endl;
          std::cout << "MomRef :" << id_track << " " << p_value << " " << mass << std::endl;
#endif

          RecoEvent.DAF_results.insert(std::make_pair(id_track, tempResults));
// RecoEvent.Sigma[ntrack]=CovSigma;
#ifdef DEBUG_KALMAN
          std::cout << "Kalman pval " << p_value2 << " charge " << charge << std::endl;
#endif
        }
    }

#ifdef DISPLAY
  display->setOptions("ABDEFGHMPT"); // G show geometry
  // if (matFX) display->setOptions("ABDEFGHMPT"); // G show geometry
  display->open();
#endif

  return 0;
}
