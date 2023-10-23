#include "TKalmanFilter_DAF_PID.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"
#include "KalmanFittedStateOnPlane.h"
#include "KalmanFitterInfo.h"
#include "StateOnPlane.h"
#include "TRandom3.h"
//#include "IO.h"
#include "Exception.h"
#include "KalmanFitStatus.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TMath.h"
#include "TMinuit.h"
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
template<class Out>
TKalmanFilter_DAF_PID<Out>::TKalmanFilter_DAF_PID(const THyphiAttributes& attribut)
    : TDataProcessInterface<Out>("mom_fit_kalman_pid"), att(attribut)
{

  const int nIter    = 10;    // max number of iterations
  const double dPVal = 1.e-3; // convergence criterion

  // genfit::errorOut.rdbuf(errorOut.rdbuf());
  // genfit::printOut.rdbuf(printOut.rdbuf());

  // const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedAverage;
  // For Ref
  // const genfit::eMultipleMeasurementHandling mmHandling = genfit::unweightedClosestToReference;
  // const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedClosestToReference;
  // const genfit::eMultipleMeasurementHandling mmHandling = genfit::unweightedClosestToReferenceWire;
  // const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedClosestToReferenceWire;

  // const genfit::eMultipleMeasurementHandling mmHandling = genfit::unweightedClosestToPrediction;
  // const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedClosestToPrediction;
  //const genfit::eMultipleMeasurementHandling mmHandling = genfit::unweightedClosestToPredictionWire;
  // const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedClosestToPredictionWire;

  //Fitter_rescue = new genfit::DAF();
  Fitter_rescue = new genfit::KalmanFitterRefTrack(nIter, dPVal); //,1e4);
  //Fitter_rescue->setMultipleMeasurementHandling(mmHandling);
  Fitter_rescue->setMinIterations(2);
  Fitter_rescue->setMaxIterations(nIter);

  Fitter     = nullptr;
  Fitter_pid = nullptr;

  if(att.KF_Kalman || att.KF_KalmanSqrt)
    {
      Fitter     = new genfit::KalmanFitter(10, dPVal, 1e3, att.KF_KalmanSqrt); //
      Fitter_pid = new genfit::KalmanFitter(10, dPVal, 1e3, att.KF_KalmanSqrt); //
    }
  else if(att.KF_KalmanRef)
    {
      Fitter     = new genfit::KalmanFitterRefTrack(nIter, dPVal); //,1e4);
      Fitter_pid = new genfit::KalmanFitterRefTrack(nIter, dPVal); //,1e4);
    //Fitter->setMultipleMeasurementHandling(mmHandling);
    }
  else if(att.KF_DAF || att.KF_DAFRef)
    {
      Fitter     = new genfit::DAF(att.KF_DAFRef);
      Fitter_pid = new genfit::DAF(att.KF_DAFRef);
    }
  else
    att._logger->error("E> Kalman Filter fitter not set correctly !");

  Fitter->setMinIterations(3);
  Fitter->setMaxIterations(nIter);

  Fitter_pid->setMinIterations(3);
  Fitter_pid->setMaxIterations(nIter);

  Nb_CentralCut = att.KF_NbCentralCut;
  Nb_MiniFiberCut = att.KF_NbMiniFiberCut;

  if(att.PID_CutorProb == true)
    {
      CutPID.insert(std::pair<int,CutMomBeta*>(-211,new CutMomBeta(0.13957,0.03,true,true)));
      CutPID.insert(std::pair<int,CutMomBeta*>(211,new CutMomBeta(0.13957,0.5,true,false)));
      //CutPID.insert(std::pair<int,CutMomBeta*>(211,new CutMomBeta(0.13957,0.009,true,false)));
      CutPID.insert(std::pair<int,CutMomBeta*>(321,new CutMomBeta(0.4936,0.12,true,false)));
      CutPID.insert(std::pair<int,CutMomBeta*>(2212,new CutMomBeta(0.938272013,0.2,true,false)));
      CutPID.insert(std::pair<int,CutMomBeta*>(-991,new CutMomBeta(0.4936,0.02,true,true)));
      CutPID[-991]->SetMargingLimit(.4);
      CutPID.insert(std::pair<int,CutMomBeta*>(991,new CutMomBeta(1.90925,0.02,true,false)));
      CutPID[991]->SetMargingLimit(.4);
    }

  // Fitter->setDebugLvl(10);

  rep = new genfit::RKTrackRep();

  Vtracks = new genfit::Track(rep, TVector3(0., 0., 0.), TVector3(0., 0., 1.0));

#ifdef DISPLAY
  display = genfit::EventDisplay::getInstance();
  display->reset();
  att._logger->info(" Display On :{}", fmt::ptr(display));
#endif
}

template<class Out>
TKalmanFilter_DAF_PID<Out>::~TKalmanFilter_DAF_PID()
{
  if(Vtracks != nullptr)
    {
      delete Vtracks;
      Vtracks = nullptr;
    }

  if(Fitter != nullptr)
    {
      delete Fitter;
      Fitter = nullptr;
    }
  if(Fitter_pid != nullptr)
    {
      delete Fitter_pid;
      Fitter_pid = nullptr;
    }
  if(Fitter_rescue != nullptr)
    {
      delete Fitter_rescue;
      Fitter_rescue = nullptr;
    }

  for(auto& elem : list_Plane)
    if(elem != nullptr)
      {
        delete elem;
        elem = nullptr;
      }
}

template<class Out>
void TKalmanFilter_DAF_PID<Out>::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }

template<class Out>
ReturnRes::InfoM TKalmanFilter_DAF_PID<Out>::operator()(FullRecoEvent& RecoEvent, Out* OutTree)
{

  int result_mom = Exec(RecoEvent, OutTree);

  ++Nb_event;
#ifdef DISPLAY
  if(Nb_event == att.NEvent)
    {
      display->setOptions("ABDEFGHMPT"); // G show geometry
      // if (matFX) display->setOptions("ABDEFGHMPT"); // G show geometry
      display->open();
    }
#endif

  return SoftExit(result_mom);
}

template<class Out>
void TKalmanFilter_DAF_PID<Out>::SelectHists()
{

  LocalHisto.h_stats         = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_stats);
  LocalHisto.h_statsLess3Mes = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_statsLess3Mes);
  LocalHisto.h_statsInvalid  = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_statsInvalid);
  LocalHisto.h_pv            = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_pv);
  LocalHisto.h_chi2          = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_chi2);
  for(size_t i = 0; i < 2; ++i)
    {
      LocalHisto.hd_chi[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->hd_chi[i]);
      LocalHisto.hd_pv[i]  = this->AnaHisto->CloneAndRegister(this->AnaHisto->hd_pv[i]);
    }
  LocalHisto.h_Path             = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Path);
  LocalHisto.h_MeanPath         = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_MeanPath);
  LocalHisto.h_beta             = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_beta);
  LocalHisto.h_beta2            = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_beta2);
  LocalHisto.h_beta3            = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_beta3);
  LocalHisto.h_Mass_All         = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Mass_All);
  LocalHisto.h_Mass_All2        = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Mass_All2);
  LocalHisto.h_Mass_All3        = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Mass_All3);
  LocalHisto.h_Mass_charge_All  = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Mass_charge_All);
  LocalHisto.h_Mass_charge_All2 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Mass_charge_All2);
  LocalHisto.h_Mass_charge_All3 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Mass_charge_All3);
  LocalHisto.h_beta_mom         = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_beta_mom);
  LocalHisto.h_beta_mom2        = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_beta_mom2);
  LocalHisto.h_beta_mom3        = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_beta_mom3);
  LocalHisto.h_beta_momcharge   = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_beta_momcharge);
  LocalHisto.h_beta_momcharge2  = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_beta_momcharge2);
  LocalHisto.h_beta_momcharge3  = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_beta_momcharge3);
  LocalHisto.h_pv_mom           = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_pv_mom);
  LocalHisto.h_pv_beta          = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_pv_beta);
  LocalHisto.h_pv_mass          = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_pv_mass);
  LocalHisto.h_path_tof         = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_path_tof);
  LocalHisto.h_mom_tof_cut      = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_mom_tof_cut);
  LocalHisto.h_path_mom_cut     = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_path_mom_cut);
  LocalHisto.h_path_tof_cut     = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_path_tof_cut);
  for(size_t i = 0; i < 4; ++i)
    {
      LocalHisto.h_Mass[i]          = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_Mass[i]);
      LocalHisto.h_chi2_particle[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_chi2_particle[i]);
      LocalHisto.h_pv_particle[i]   = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_pv_particle[i]);
    }
  for(size_t i = 0; i < 5; ++i)
    {
      LocalHisto.h_mom_res[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_mom_res[i]);
      for(size_t j = 0; j < 10; ++j)
        {
          LocalHisto.h_ResPull[i][j]        = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResPull[i][j]);
          LocalHisto.h_ResPull_normal[i][j] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResPull_normal[i][j]);
        }
    }
  LocalHisto.h_total_dE = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_total_dE);

  for(size_t i = 0;i < 9 ;++i)
    LocalHisto.h_ResFiber[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResFiber[i]);

  for(size_t i = 0;i < 6 ;++i)
    LocalHisto.h_ResMiniFiber[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResMiniFiber[i]);

  for(size_t i = 0;i < 17 ;++i)
    for(size_t j =0 ; j<3;++j)
      LocalHisto.h_ResMDC[i][j] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResMDC[i][j]);

  for(size_t i = 0;i < 2 ;++i)
    LocalHisto.h_ResPSCE[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResPSCE[i]);
  
}
  

template<class Out>
int TKalmanFilter_DAF_PID<Out>::Exec(FullRecoEvent& RecoEvent, Out* OutTree)
{

  //  int result_kalman = Kalman_Filter(RecoEvent);
  int result_kalman = Kalman_Filter_FromTrack(RecoEvent);

  for(auto TrackRes : RecoEvent.DAF_results)
    {
      int TrackID             = TrackRes.first;
      const ResSolDAF& FitRes = TrackRes.second;
      auto TInfo              = RecoEvent.TrackInfo.find(TrackID);
      auto IsDecay            = RecoEvent.TrackMother.find(TrackID);
      int Decay               = IsDecay == RecoEvent.TrackMother.end() ? 0 : 1;

      THyphiTrack* OutTrack = dynamic_cast<THyphiTrack*>(OutTree->fTrack->ConstructedAt(OutTree->fTrack->GetEntries()));
//    auto PDG_particle     = TDatabasePDG::Instance()->GetParticle(FitRes.pdg_guess);
      auto PDG_particle     = TDatabasePDG::Instance()->GetParticle(TInfo->second[FitRes.lastHit].pdg);
      OutTrack->type        = PDG_particle->GetName();
      OutTrack->MC_status   = TrackID + 10000 * Decay;
      OutTrack->Chi2        = FitRes.chi2;
      OutTrack->Chi2_X      = FitRes.ndf;
      OutTrack->Chi2_Y      = FitRes.firstHit;
      OutTrack->Mass        = FitRes.mass;
      OutTrack->pdgcode     = FitRes.pdg_guess;
      OutTrack->MomMass.SetXYZM(FitRes.momX, FitRes.momY, FitRes.momZ, FitRes.mass);
      OutTrack->Mom.SetXYZ(FitRes.momX, FitRes.momY, FitRes.momZ);

      OutTrack->BarId  = FitRes.lastHit;
      OutTrack->Charge = FitRes.charge;
      OutTrack->dE     = TInfo->second[FitRes.lastHit].Eloss;
      OutTrack->Beta   = FitRes.beta;
      OutTrack->RefPoint.SetXYZ(FitRes.posX, FitRes.posY, FitRes.posZ);
      OutTrack->Pval2      = FitRes.pvalue;
      OutTrack->PathLength = FitRes.path_length;
      OutTrack->TOF        = FitRes.tof;

      OutTrack->MomIni.SetXYZ(FitRes.momX_init, FitRes.momY_init, FitRes.momZ_init);
      OutTrack->BetaIni       = FitRes.beta2;
      OutTrack->MassIni       = FitRes.mass2;
      OutTrack->TOFIni        = FitRes.tof2;
      OutTrack->PathLengthIni = TInfo->second[FitRes.lastHit].length; // FitRes.path_length2;
      OutTrack->RChiIni       = FitRes.fitter;
      if(Decay == 1)
        OutTrack->Sim2Vtx.SetXYZT(std::get<1>(IsDecay->second), std::get<2>(IsDecay->second),
                                  std::get<3>(IsDecay->second), std::get<4>(IsDecay->second));

      OutTrack->State[0] = FitRes.posX;
      OutTrack->State[1] = FitRes.posY;
      OutTrack->State[2] = FitRes.posZ;
      OutTrack->State[3] = FitRes.momX;
      OutTrack->State[4] = FitRes.momY;
      OutTrack->State[5] = FitRes.momZ;
      for(int row = 0; row < 6; row++)
        for(int col = 0; col < 6; col++)
          OutTrack->Cov[row][col] = FitRes.cov_matrix[row][col];

      OutTrack->NCent   = FitRes.Ncentral;
      OutTrack->Nmfib  = FitRes.Nmfiber;
      OutTrack->iterNum = FitRes.iterNum;
      for(int i = 0; i < 17; ++i)
        {
          for(int j = 0; j < 3; ++j)
            {
              OutTrack->ResMDC[i][j]    = FitRes.ResMDC[i][j];
              OutTrack->WeightMDC[i][j] = FitRes.WeightMDC[i][j];
            }
        }
      for(int i = 0; i < 9; ++i)
        {
          OutTrack->ResFiber[i]    = FitRes.ResFiber[i];
          OutTrack->WeightFiber[i] = FitRes.WeightFiber[i];
        }
      for(int i=0; i<6; ++i)
        {
          OutTrack->ResMiniFiber[i] = FitRes.ResMiniFiber[i];
          OutTrack->WeightMiniFiber[i] = FitRes.WeightMiniFiber[i];
        }
      for(int i = 0; i < 2; ++i)
        {
          OutTrack->ResPSCE[i] = FitRes.ResPSCE[i];
        }
    }

  OutTree->Ntrack = OutTree->fTrack->GetEntries();

  // #ifdef DEBUG_DAISUKE
  //       std::cout<<"HitTR1
  //       ("<<it_mom_pv->second[13]<<","<<it_mom_pv->second[14]<<","<<it_mom_pv->second[15]<<")"<<std::endl;
  //       std::cout<<"HitTR2
  //       ("<<it_mom_pv->second[16]<<","<<it_mom_pv->second[17]<<","<<it_mom_pv->second[18]<<")"<<std::endl;
  //       std::cout<<"HitDC2
  //       ("<<it_mom_pv->second[19]<<","<<it_mom_pv->second[20]<<","<<it_mom_pv->second[21]<<")"<<std::endl;
  //       std::cout<<"NumDC2
  //       ("<<it_mom_pv->second[22]<<","<<it_mom_pv->second[23]<<","<<it_mom_pv->second[24]<<")"<<std::endl;

  //       std::cout<<"Pval2 = "<<it_mom_pv->second[27]<<std::endl;
  //       std::cout<<"TofsBar = "<<it_mom_pv->second[28]<<std::endl;
  //       std::cout<<"BarId = "<<it_mom_pv->second[25]<<std::endl;
  //       std::cout<<"Charge = "<<it_mom_pv->second[6]<<std::endl;
  //       std::cout<<"Beta = "<<it_mom_pv->second[26]<<std::endl;
  // #endif

  return result_kalman;
}

template<class Out>
ReturnRes::InfoM TKalmanFilter_DAF_PID<Out>::SoftExit(int result_full) { return ReturnRes::Fine; }

template<class Out>
int TKalmanFilter_DAF_PID<Out>::Kalman_Filter_FromTrack(FullRecoEvent& RecoEvent)
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

  // unsigned int Nb_trackSeed = RecoEvent.TrackDAF.size();

  int ntrack = -1;
  for(auto it_trackInfo : RecoEvent.TrackInfo)
    {
      ntrack++;

#ifdef DEBUG_KALMAN
      att._logger->debug("start #{}", ntrack); //<<std::endl;
                                               // std::cout<<" | Pointsize "<<Vtracks->getNumPoints()<<" Measurements
                                               // "<<Vtracks->getNumPointsWithMeasurement()<<" Rep
                                               // "<<Vtracks->getNumReps()<<std::endl;
                                               // std::vector<std::string> name_det;
#endif

      const int id_track = it_trackInfo.first;
      auto it_ListHits   = RecoEvent.TrackDAF.find(id_track);

      // auto getZpos = [](genfit::AbsMeasurement* m) {
      //   TVectorD& HitrawRef = m->getRawHitCoords();
      //   if(HitrawRef.GetNrows() == 2)
      //     {
      //       genfit::StateOnPlane dummy;
      //       genfit::SharedPlanePtr plane = dynamic_cast<genfit::PlanarMeasurement*>(m)->constructPlane(dummy);
      //       TVector3 tempO(plane->getO());
      //       return tempO.Z();
      //     }
      //   else if(HitrawRef.GetNrows() == 3)
      //     return HitrawRef[2];
      //   else
      //     {
      //       fmt::print("E> rawref not proper ! {}", HitrawRef.GetNrows());
      //       HitrawRef.Print();
      //       return -999.;
      //     }
      // };
      std::set<std::tuple<double, int, int> > id_dets;
      double total_dE = 0.;
      int n_Central   = 0;
      int n_MiniFiber = 0;

      //std::cout << "\nTrackID: " << it_ListHits->first << "\n";

      for(size_t id_det = 0; id_det < it_ListHits->second.size(); ++id_det)
        {
          int id_hit = it_ListHits->second[id_det];
          if(id_hit < 0)
            continue;
          if(id_det >= G4Sol::Si1x && id_det <= G4Sol::Si2y)
            continue;
          if(id_det >= G4Sol::Si1x_SD && id_det <= G4Sol::Si2y_SD_pad)
            continue;
          if(id_det >= G4Sol::MG01 && id_det <= G4Sol::MG17)
            ++n_Central;
          if(id_det>=G4Sol::MiniFiberD1_x && id_det<=G4Sol::MiniFiberD2_u)
            ++n_MiniFiber;
          // genfit::AbsMeasurement* currentHit = RecoEvent.ListHits[id_det][id_hit].get();

          total_dE += it_trackInfo.second[id_det].Eloss;   //CHECK Change!!

          id_dets.insert(std::make_tuple(id_det, id_det, id_hit));

          //std::cout << "DetID: " << id_det << "  HitID: " << currentHit->getHitId() << "\n";
        }

      if(id_dets.size() < 3) //(nb_ValidHits < 3)
        {
#ifdef DEBUG_KALMAN
          att._logger->debug("!> less than 3 measurements: {}", id_dets.size());
#endif
          LocalHisto.h_stats->Fill("Less3Mes", 1);
          int idPDG = 0;
          for(size_t i = 0; i < it_trackInfo.second.size(); ++i)
            if(it_trackInfo.second[i].pdg != idPDG)
              {
                idPDG = it_trackInfo.second[i].pdg;
                break;
              }
          if(idPDG==0)
            {
              att._logger->debug("!> less than 3 measurements: PDG = 0 : id# {}" , id_track);
              for(size_t i = 0; i < it_trackInfo.second.size(); ++i)
                att._logger->debug("!> #{} : {} {} {}",i,it_trackInfo.second[i].pdg, it_trackInfo.second[i].mass);
            }
          std::string namePDG = std::to_string(idPDG);
          LocalHisto.h_statsLess3Mes->Fill(namePDG.c_str(), "Less3Mes", 1.);

          continue;
        }
      if(n_Central < Nb_CentralCut)
        {
          LocalHisto.h_stats->Fill("Less3MesCentral", 1);
          int idPDG = 0;
          for(size_t i = 0; i < it_trackInfo.second.size(); ++i)
            if(it_trackInfo.second[i].pdg != idPDG)
              {
                idPDG = it_trackInfo.second[i].pdg;
                break;
              }
          std::string namePDG = std::to_string(idPDG);
          LocalHisto.h_statsLess3Mes->Fill(namePDG.c_str(), "Less3MesCentral", 1.);
          continue;
        }
      if(n_MiniFiber< Nb_MiniFiberCut)
        {
          LocalHisto.h_stats->Fill("LessMiniFiber<4",1);
          int idPDG = 0;
          for(size_t i = 0; i < it_trackInfo.second.size(); ++i)
            if(it_trackInfo.second[i].pdg!=idPDG)
            {
              idPDG = it_trackInfo.second[i].pdg;
              break;
            }
          std::string namePDG = std::to_string(idPDG);
          LocalHisto.h_statsLess3Mes->Fill(namePDG.c_str(),"LessMiniFiber<4",1.);
          continue;
        }

      auto f_LastHitIsValid = [](const auto& it_ListHits, const std::set<G4Sol::SolDet>& listToTest) {
        for(auto it_det : listToTest)
          if(it_ListHits->second[it_det] >= 0)
            return it_det;

        return G4Sol::SIZEOF_G4SOLDETTYPE;
      };

      G4Sol::SolDet LastFrontWall = att.Wasa_Side == 0 ? G4Sol::PSFE : G4Sol::PSBE;
      G4Sol::SolDet lastValidHit =
          f_LastHitIsValid(it_ListHits, {G4Sol::PSCE, LastFrontWall, G4Sol::RPC_l, G4Sol::RPC_h});

      if(lastValidHit == G4Sol::SIZEOF_G4SOLDETTYPE)
        {
          int idPDG = 0;
          for(size_t i = 0; i < it_trackInfo.second.size(); ++i)
            if(it_trackInfo.second[i].pdg != idPDG)
              {
                idPDG = it_trackInfo.second[i].pdg;
                break;
              }

#ifdef DEBUG_KALMAN
          att._logger->debug("!> LastValidHit not found !");
          std::vector<std::stringstream> s1(it_trackInfo.second.size() / 20 + 1);
          std::vector<std::stringstream> s2(it_trackInfo.second.size() / 20 + 1);
          std::vector<std::stringstream> s3(it_trackInfo.second.size() / 20 + 1);
          for(size_t i = 0; i < it_trackInfo.second.size(); ++i)
            {
              s1[i / 20] << printW(G4Sol::nameLiteralDet.begin()[i], 6) << ", ";
              s2[i / 20] << printW(i, 6) << ", ";
              s3[i / 20] << printW(it_trackInfo.second[i].pdg, 6) << ", ";
            }
          for(size_t i = 0; i < s1.size(); ++i)
            {
              att._logger->debug("Det  :{}", s1[i].str());
              att._logger->debug("idDet:{}", s2[i].str());
              att._logger->debug("stat :{}", s3[i].str());
            }
#endif
          LocalHisto.h_stats->Fill("NoValidLastHit", 1.);

          std::string namePDG = std::to_string(idPDG);
          for(auto it_det = it_ListHits->second.crbegin(), it_det_end = it_ListHits->second.crend();
              it_det != it_det_end; ++it_det)
            if(*it_det >= 0)
              {
                int id_det_inv = it_ListHits->second.size() - 1 - std::distance(it_ListHits->second.crbegin(), it_det);
                LocalHisto.h_statsInvalid->Fill(namePDG.c_str(), G4Sol::nameLiteralDet.begin()[id_det_inv], 1.);
                break;
              }
          continue;
        }

      auto firstHit   = id_dets.cbegin();
      int id_firstDet = std::get<1>(*firstHit);

      if(TMath::Abs(it_trackInfo.second[id_firstDet].pdg) < 1e-2)
        {
#ifdef DEBUG_KALMAN
          att._logger->debug("!> pdgcode = 0 {}", it_trackInfo.second[id_firstDet].pdg);
#endif
          continue;
        }
#ifdef DEBUG_KALMAN
      att._logger->debug(" Id_track: {}", id_track);
#endif

      const InfoPar track_state(it_trackInfo.second[id_firstDet]);

      genfit::AbsMeasurement* tempHit = RecoEvent.ListHits[id_firstDet][std::get<2>(*firstHit)].get();
      TVectorD& tempHitrawRef         = tempHit->getRawHitCoords();

      int PDG     = -211; //pi minus
      if(att.G4_simu && att.WF_perfect)
        PDG     = static_cast<int>(track_state.pdg);

      auto PDG_particle = TDatabasePDG::Instance()->GetParticle(PDG);
      
      LocalHisto.h_stats->Fill("Beginning Kalman", 1.);

      std::string temp_namePDG = std::to_string(PDG);
      LocalHisto.h_statsLess3Mes->Fill(temp_namePDG.c_str(),"Beginning Kalman", 1.);

      LocalHisto.h_total_dE->Fill(PDG_particle->GetName(), total_dE, 1.);

      if(PDG_particle == nullptr)
        {
          att._logger->error("E> PDG not found !");
          continue;
        }
      const double charge = PDG_particle->Charge() / 3.;

      if(charge <= 1)
        {
          std::set<std::tuple<double, int, int> > temp_id_dets;
          for(auto id_det : id_dets)
            {
              if(std::get<1>(id_det) <= G4Sol::RPC_h)
                temp_id_dets.insert(id_det);
            }

          std::swap(id_dets, temp_id_dets);
        }

      //std::cout << "AfterCuts YES \n\n";

      firstHit             = id_dets.cbegin();
      id_firstDet          = std::get<1>(*firstHit);
      const auto lastHit   = id_dets.crbegin();
      const int id_lastDet = std::get<1>(*lastHit); // std::get<0>(lastHit);

      //std::cout << "FirstHit: " << id_firstDet << "  #det: " << G4Sol::nameLiteralDet.begin()[id_firstDet] << "\n";
      //std::cout << " LastHit: " << id_lastDet  << "  #det: " << G4Sol::nameLiteralDet.begin()[id_lastDet]  << "\n";

#ifdef DEBUG_KALMAN
      for(auto idet : id_dets)
        att._logger->debug("det:{} [{}] at Z:{}", std::get<1>(idet), std::get<2>(idet), std::get<0>(idet));

      const int PID = charge * charge;
      att._logger->debug(" PID:{} {} {}", PID, charge, PDG);
#endif

      const InfoPar track_stateLast(it_trackInfo.second[id_lastDet]);

      // Forward

      auto it_init   = RecoEvent.TrackDAFInit.find(id_track);
      double init_px = it_init->second.momX;
      double init_py = it_init->second.momY;
      double init_pz = it_init->second.momZ;

      TVector3 init_p(init_px, init_py, init_pz);

#ifdef DEBUG_KALMAN
      att._logger->debug("init mom : ");
      init_p.Print();
#endif

      // Forward
      genfit::AbsTrackRep* REP = new genfit::RKTrackRep(PDG, 1);
      // REP->setDebugLvl(2);

      double init_x = it_init->second.posX;
      double init_y = it_init->second.posY;
      double init_z = it_init->second.posZ;

      const TVector3 init_point(init_x, init_y, init_z);

      TMatrixDSym CovM(6);
      for(int i = 0; i < 6; ++i)
        CovM(i, i) = 1000.;

      genfit::MeasuredStateOnPlane stateRef(REP);
      REP->setPosMomCov(stateRef, init_point, init_p, CovM);
      // remember original initial state
      const genfit::StateOnPlane stateRefOrig(stateRef.getState(), stateRef.getPlane(), REP, stateRef.getAuxInfo());

      // create track
      TVectorD seedState(6);
      TMatrixDSym seedCov(6);
      REP->get6DStateCov(stateRef, seedState, seedCov);

#ifdef DEBUG_KALMAN
      att._logger->debug("init_p:");
      init_p.Print();
      init_p.Unit().Print();
      att._logger->debug("seed_p:");
      seed_p.Print();
      att._logger->debug("init_point :");
      init_point.Print();
      att._logger->debug(" SEED :");
      seedState.Print();
      att._logger->debug("--");
      seedCov.Print();
      // std::cout<<"init_plane :"<<std::endl;
      // init_plane.Print();
#endif

      if(REP == NULL)
        att._logger->error("E> no Rep build");
#ifdef DEBUG_KALMAN
      else
        att._logger->debug("rep done");
#endif

#ifndef DISPLAY
      std::unique_ptr<genfit::Track> fitTrack = std::make_unique<genfit::Track>(REP, seedState, seedCov);
#else
      std::shared_ptr<genfit::Track> fitTrack = std::make_shared<genfit::Track>(REP, seedState, seedCov);
#endif

      for(auto ids : id_dets)
        {
          int id_det = std::get<1>(ids);
          int id_hit = std::get<2>(ids);

          genfit::AbsMeasurement* currentHit = RecoEvent.ListHits[id_det][id_hit].get();

          //std::cout << "Loop insertPoint: #det: " << G4Sol::nameLiteralDet.begin()[id_det] << " #hit: " << id_hit << "\n";

#ifdef DEBUG_KALMAN
          att._logger->debug("Loop insertPoint: #det:{} #hit {} {}", G4Sol::nameLiteralDet.begin()[id_det], id_hit,
                             fmt::ptr(currentHit));
          currentHit->Print();
#endif

          fitTrack->insertPoint(new genfit::TrackPoint(currentHit, fitTrack.get(), false));
        }


      fitTrack->checkConsistency();
#ifdef DEBUG_KALMAN2
      att._logger->debug("track n'{}", ntrack);
#endif
      int whoDoFit = 1;
      try
        {
          Fitter->processTrack(fitTrack.get(), true);
        }
      catch(genfit::Exception& e1)
        {

          att._logger->info("*** FITTER EXCEPTION *** Rescue fitter take place !");
          att._logger->info(e1.what());
          LocalHisto.h_stats->Fill("Exc:RescueFit", 1.);

          try
            {
              whoDoFit = 2;
              Fitter_rescue->processTrack(fitTrack.get(), true);
            }
          catch(genfit::Exception& e2)
            {
              att._logger->info("*** FITTER Rescue EXCEPTION *** cancel fit !");
              att._logger->info(e2.what());
              LocalHisto.h_stats->Fill("Exc:FitFail", 1.);
              
              continue;
            }
        }
      
#ifdef DISPLAY
      display->addEvent(fitTrack.get()); // Vtracks);
#endif

#ifdef DEBUG_KALMAN2
      att._logger->debug("SUCCESSFUL FIT!");
      att._logger->debug("track n'{} / {}", ntrack, fitTrack->getFitStatus(REP)->isFitConverged());
      fitTrack->Print();
#endif

      if(fitTrack->getFitStatus(REP)->isFitConverged())
        {
          assert(fitTrack->hasKalmanFitStatus() == true);
          // assert(fitTrack->checkConsistency());
          
          LocalHisto.h_stats->Fill("Converged Kalman", 1.);
          LocalHisto.h_statsLess3Mes->Fill(temp_namePDG.c_str(),"Converged Kalman", 1.);

          fitTrack->checkConsistency();
          genfit::TrackPoint* tp = fitTrack->getPointWithMeasurementAndFitterInfo(0, REP);
          if(tp == NULL)
            {
              att._logger->info("Track has no TrackPoint with fitterInfo!");
              LocalHisto.h_stats->Fill("Exc:NoTrackPointInfo", 1.);
              continue;
            }

          genfit::KalmanFittedStateOnPlane* kfsop = nullptr;
#ifdef DISPLAY
          kfsop = (dynamic_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(REP))->getBackwardUpdate());
#else
          try
            {
              kfsop = (dynamic_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(REP))->getBackwardUpdate());
              // Forward
              REP->extrapolateToPlane(*kfsop, stateRefOrig.getPlane());
            }
          catch(genfit::Exception& e)
            {
              // std::cerr << "Exception, next track" << std::endl;
              LocalHisto.h_stats->Fill("Exc:FitBacktoOrig", 1.);
              // std::cerr << e.what();
              continue;
            }
#endif

          assert(kfsop != nullptr);

          TVector3 p3, posRef;
          TMatrixDSym covFit(6);
          kfsop->getPosMomCov(posRef, p3, covFit);

          TVector3 posBegin;
          try{
            genfit::TrackPoint* tp_tmp = fitTrack->getPointWithMeasurementAndFitterInfo(0, REP);
            genfit::KalmanFitterInfo* kfi = static_cast<genfit::KalmanFitterInfo*>(tp_tmp->getFitterInfo(REP));
            const genfit::MeasuredStateOnPlane& buf_state = kfi->getFittedState();
            TVector3 v_pos;
            TVector3 v_mom;
            TMatrixDSym v_cov;
            buf_state.getPosMomCov(v_pos, v_mom, v_cov);
            posBegin = v_pos;
          }
          catch(genfit::Exception& e){
            std::cerr<<"Exception posBegin, next track"<<std::endl;
            std::cerr << e.what();
            continue;
          }
          
          //std::cout << Form("- posBegin   : (%.2f, %.2f, %.2f)",   posBegin.x(),   posBegin.y(),   posBegin.z()) << std::endl;
          //std::cout << Form("- init_point : (%.2f, %.2f, %.2f)", init_point.x(), init_point.y(), init_point.z()) << std::endl;

          double p = p3.Mag();
          double Path_length = -999.;
          try
            {
              Path_length = fitTrack->getTrackLen(REP);
            }
          catch(genfit::Exception& e)
            {
              att._logger->info("could not get TrackLen ! ");
              att._logger->info(e.what());
              LocalHisto.h_stats->Fill("Exc:TrackLen", 1.);
              continue;
            }

          //std::cout << Form("- TOF : %.4f",   track_stateLast.time) << std::endl;

          double init_time = it_init->second.time;

          double Path_lengthMean = Path_length;
          Path_lengthMean += (posBegin-init_point).Mag();
          double time_of_flight          = track_stateLast.time - init_time;
          double beta = 1. / 29.9792458 * Path_lengthMean / time_of_flight;
          double mass  = p * sqrt(1. / (beta * beta) - 1.);

          int PDG_guess = 2212;
          if(att.PID_CutorProb == false)
            {
              if(att.cut_pi->IsInside(mass*mass, p) && kfsop->getCharge()>0) PDG_guess =  211;
              if(att.cut_pi->IsInside(mass*mass, p) && kfsop->getCharge()<0) PDG_guess = -211;
            }
          else
            {
              if(kfsop->getCharge()>0) PDG_guess =  ProbPIDAssign_Pos(p, beta);
              if(kfsop->getCharge()<0) PDG_guess =  ProbPIDAssign_Neg(p, beta);
            }

          if(std::abs(PDG_guess) == 991) PDG_guess = 2212;


          // GENFIT PID
          genfit::AbsTrackRep* rep_pid = new genfit::RKTrackRep(PDG_guess, 1);

          // smeared start state
          genfit::MeasuredStateOnPlane state_pid(rep_pid);

          TVector3 posM_pid = posRef;
          TVector3 momM_pid = p3;

          state_pid.setPosMomCov(posM_pid, momM_pid, covFit);

          // create track
          TVectorD seedState_pid(6);
          TMatrixDSym seedCov_pid(6);
          state_pid.get6DStateCov(seedState_pid, seedCov_pid);

          if(rep_pid == NULL)
            att._logger->error("E> no Rep build");
#ifdef DEBUG_KALMAN
          else
            att._logger->debug("rep done");
#endif

#ifndef DISPLAY
          std::unique_ptr<genfit::Track> fitTrack_pid = std::make_unique<genfit::Track>(rep_pid, seedState_pid, seedCov_pid);
#else
          std::shared_ptr<genfit::Track> fitTrack_pid = std::make_shared<genfit::Track>(rep_pid, seedState_pid, seedCov_pid);
#endif

          for(auto ids : id_dets)
            {
              int id_det = std::get<1>(ids);
              int id_hit = std::get<2>(ids);

              genfit::AbsMeasurement* currentHit = RecoEvent.ListHits[id_det][id_hit].get();

              //std::cout << "Loop insertPoint: #det: " << G4Sol::nameLiteralDet.begin()[id_det] << " #hit: " << id_hit << "\n";

#ifdef DEBUG_KALMAN
              att._logger->debug("Loop insertPoint 2nd fit: #det:{} #hit {} {}", G4Sol::nameLiteralDet.begin()[id_det], id_hit,
                                 fmt::ptr(currentHit));
              currentHit->Print();
#endif

              fitTrack_pid->insertPoint(new genfit::TrackPoint(currentHit, fitTrack_pid.get(), false));
            }


          fitTrack_pid->checkConsistency();

          int whoDoFit_pid = 1;
          try
            {
              Fitter_pid->processTrack(fitTrack_pid.get(), true);
            }
          catch(genfit::Exception& e1)
            {

              att._logger->info("*** FITTER EXCEPTION *** Rescue fitter take place - 2nd fit!");
              att._logger->info(e1.what());
              LocalHisto.h_stats->Fill("Exc:RescueFit-2", 1.);

              try
                {
                  whoDoFit_pid = 2;
                  Fitter_rescue->processTrack(fitTrack_pid.get(), true);
                }
              catch(genfit::Exception& e2)
                {
                  att._logger->info("*** FITTER Rescue EXCEPTION *** cancel fit - 2nd fit!");
                  att._logger->info(e2.what());
                  LocalHisto.h_stats->Fill("Exc:FitFail-2", 1.);
                  
                  continue;
                }
            }
          
#ifdef DISPLAY
          display->addEvent(fitTrack_pid.get());
#endif

#ifdef DEBUG_KALMAN2
          att._logger->debug("SUCCESSFUL FIT - 2nd!");
          att._logger->debug("track n'{} / {}", ntrack, fitTrack_pid->getFitStatus(rep_pid)->isFitConverged());
          fitTrack_pid->Print();
#endif


          if(fitTrack_pid->getFitStatus(rep_pid)->isFitConverged())
            {
              assert(fitTrack_pid->hasKalmanFitStatus() == true);
              // assert(fitTrack->checkConsistency());
              
              LocalHisto.h_stats->Fill("Converged Kalman - 2nd", 1.);
              LocalHisto.h_statsLess3Mes->Fill(temp_namePDG.c_str(),"Converged Kalman - 2nd", 1.);

              fitTrack_pid->checkConsistency();
              genfit::TrackPoint* tp_pid = fitTrack_pid->getPointWithMeasurementAndFitterInfo(0, rep_pid);
              if(tp_pid == NULL)
                {
                  att._logger->info("Track has no TrackPoint with fitterInfo - 2nd!");
                  LocalHisto.h_stats->Fill("Exc:NoTrackPointInfo-2", 1.);
                  continue;
                }

              genfit::KalmanFittedStateOnPlane* kfsop_pid = nullptr;

#ifdef DISPLAY
              kfsop_pid = (dynamic_cast<genfit::KalmanFitterInfo*>(tp_pid->getFitterInfo(rep_pid))->getBackwardUpdate());
#else
              try
                {
                  kfsop_pid = (dynamic_cast<genfit::KalmanFitterInfo*>(tp_pid->getFitterInfo(rep_pid))->getBackwardUpdate());
                  // Forward
                  rep_pid->extrapolateToPlane(*kfsop_pid, state_pid.getPlane());
                }
              catch(genfit::Exception& e)
                {
                  // std::cerr << "Exception, next track" << std::endl;
                  LocalHisto.h_stats->Fill("Exc:FitBacktoOrig-2", 1.);
                  // std::cerr << e.what();
                  continue;
                }
#endif

              assert(kfsop_pid != nullptr);

              TVector3 p3_pid, posRef_pid;
              TMatrixDSym covFit_pid(6);
              kfsop_pid->getPosMomCov(posRef_pid, p3_pid, covFit_pid);

              TVector3 posBegin_pid;
              try{
                genfit::TrackPoint* tp_tmp = fitTrack_pid->getPointWithMeasurementAndFitterInfo(0, rep_pid);
                genfit::KalmanFitterInfo* kfi = static_cast<genfit::KalmanFitterInfo*>(tp_tmp->getFitterInfo(rep_pid));
                const genfit::MeasuredStateOnPlane& buf_state = kfi->getFittedState();
                TVector3 v_pos;
                TVector3 v_mom;
                TMatrixDSym v_cov;
                buf_state.getPosMomCov(v_pos, v_mom, v_cov);
                posBegin_pid = v_pos;
              }
              catch(genfit::Exception& e){
                std::cerr<<"Exception posBegin, next track"<<std::endl;
                std::cerr << e.what();
                continue;
              }

              //std::cout << Form("- posBeginPID: (%.3f, %.3f, %.3f)",   posBegin_pid.x(),   posBegin_pid.y(),   posBegin_pid.z()) << std::endl;

              double p_pid = p3_pid.Mag();
              double Path_length_pid = -999.;
              try
                {
                  Path_length_pid = fitTrack_pid->getTrackLen(rep_pid);
                }
              catch(genfit::Exception& e)
                {
                  att._logger->info("could not get TrackLen - 2nd! ");
                  att._logger->info(e.what());
                  LocalHisto.h_stats->Fill("Exc:TrackLen-2", 1.);
                  continue;
                }

              double Path_lengthMean_pid = Path_length_pid;
              Path_lengthMean_pid += (posBegin_pid-init_point).Mag();

              double time_of_flight_pid          = track_stateLast.time - init_time;
              double beta_pid = 1. / 29.9792458 * Path_lengthMean_pid / time_of_flight_pid;
              double mass_pid  = p_pid * sqrt(1. / (beta_pid * beta_pid) - 1.);

              int PDG_guess_pid = 2212;
              if(att.PID_CutorProb == false)
                {
                  if(att.cut_pi->IsInside(mass_pid*mass_pid, p_pid) && kfsop_pid->getCharge()>0) PDG_guess_pid =  211;
                  if(att.cut_pi->IsInside(mass_pid*mass_pid, p_pid) && kfsop_pid->getCharge()<0) PDG_guess_pid = -211;
                }
              else
                {
                  if(kfsop_pid->getCharge()>0) PDG_guess_pid =  ProbPIDAssign_Pos(p, beta);
                  if(kfsop_pid->getCharge()<0) PDG_guess_pid =  ProbPIDAssign_Neg(p, beta);
                }


              ResSolDAF tempResults;
              tempResults.charge    = kfsop_pid->getCharge();
              tempResults.pdg_ini   = PDG_guess;


              // Forward
              tempResults.momX_init = it_init->second.momX;
              tempResults.momY_init = it_init->second.momY;
              tempResults.momZ_init = it_init->second.momZ;

              tempResults.momX = p3_pid.X();
              tempResults.momY = p3_pid.Y();
              tempResults.momZ = p3_pid.Z();

              tempResults.posX = posRef_pid.X();
              tempResults.posY = posRef_pid.Y();
              tempResults.posZ = posRef_pid.Z();

              for(int row = 0; row < 6; ++row)
                for(int col = 0; col < 6; ++col)
                  tempResults.cov_matrix[row][col] = covFit_pid(row, col);

              const double chi2_pid = whoDoFit_pid == 1
                                      ? Fitter_pid->getChiSqu(fitTrack_pid.get(), rep_pid)
                                      : Fitter_rescue->getChiSqu(fitTrack_pid.get(), rep_pid);

              const double ndf_pid =
                  whoDoFit_pid == 1 ? Fitter_pid->getNdf(fitTrack_pid.get(), rep_pid)
                                    : Fitter_rescue->getNdf(fitTrack_pid.get(), rep_pid);

              tempResults.chi2    = chi2_pid;
              tempResults.ndf     = ndf_pid;
              tempResults.fitter  = whoDoFit_pid;
              tempResults.iterNum = (int)(fitTrack_pid->getKalmanFitStatus()->getNumIterations());
#ifdef DEBUG_KALMAN2
              att._logger->debug(" / chi2 ={} / ndf ={}", chi2_pid, ndf_pid);
#endif

              {
                unsigned int np_pid = fitTrack_pid->getNumPointsWithMeasurement();
                for(size_t i = 0; i < np_pid; ++i)
                  {
                    genfit::TrackPoint* tp_tmp    = fitTrack_pid->getPointWithMeasurementAndFitterInfo(i);
                      if(tp_tmp==nullptr)
                        {
                          att._logger->debug("E> TrackPoint null - 2nd! #i {} / {}\n",i,np_pid);
                          continue;
                        }
                    genfit::KalmanFitterInfo* kfi = static_cast<genfit::KalmanFitterInfo*>(tp_tmp->getFitterInfo(rep_pid));
                    std::vector<double> weights   = kfi->getWeights();
                    int id_det                    = tp_tmp->getRawMeasurement(0)->getDetId();
                    // std::cout << "id_det : " << id_det << std::endl;
                    if(id_det >= G4Sol::MG01 && id_det <= G4Sol::MG17)
                      {
                        double res_min = 9999999;
                        for(size_t j = 0; j < weights.size() ; ++j)
                          {
                            const genfit::MeasurementOnPlane& residual = kfi->getResidual(j, false, true);
                            const TVectorD& resid(residual.getState());
                            double res = resid(0);
                            if(fabs(res) < fabs(res_min))
                              res_min = res;
                            // std::cout << "res : " << res << std::endl;
                            LocalHisto.h_ResMDC[id_det - G4Sol::MG01][j]->Fill(res);
                            tempResults.ResMDC[id_det - G4Sol::MG01][j]    = res;
                            tempResults.WeightMDC[id_det - G4Sol::MG01][j] = weights[j];
                          }
                        LocalHisto.h_ResMDC[id_det - G4Sol::MG01][2]->Fill(res_min);
                        tempResults.ResMDC[id_det - G4Sol::MG01][2] = res_min;
                      }
                    else if(id_det >= G4Sol::FiberD1_x && id_det <= G4Sol::FiberD3_v)
                      {
                        const genfit::MeasurementOnPlane& residual = kfi->getResidual(0, false, true);
                        const TVectorD& resid(residual.getState());
                        double res = resid(0);
                        // std::cout << "res : " << res << std::endl;
                        LocalHisto.h_ResFiber[id_det - G4Sol::FiberD1_x]->Fill(res);
                        tempResults.ResFiber[id_det - G4Sol::FiberD1_x]    = res;
                        tempResults.WeightFiber[id_det - G4Sol::FiberD1_x] = weights[0];
                      }
                    else if(id_det>=G4Sol::MiniFiberD1_x && id_det<=G4Sol::MiniFiberD2_u)
                      {
                        const genfit::MeasurementOnPlane& residual = kfi->getResidual(0, false, true);
                        const TVectorD& resid(residual.getState());
                        double res = resid(0);
                        //std::cout << "res : " << res << std::endl;
                        LocalHisto.h_ResMiniFiber[id_det-G4Sol::MiniFiberD1_x]->Fill(res);
                        tempResults.ResMiniFiber[id_det-G4Sol::MiniFiberD1_x] = res;
                        tempResults.WeightMiniFiber[id_det-G4Sol::MiniFiberD1_x] = weights[0];
                      }
                    else if(id_det == G4Sol::PSCE)
                      {
                        const genfit::MeasurementOnPlane& residual = kfi->getResidual(0, false, true);
                        const TVectorD& resid(residual.getState());
                        double res  = resid(0);
                        double resz = resid(1);
                        // std::cout << "res : " << res << std::endl;
                        LocalHisto.h_ResPSCE[0]->Fill(res);
                        LocalHisto.h_ResPSCE[1]->Fill(resz);
                        tempResults.ResPSCE[0]    = res;
                        tempResults.ResPSCE[1]    = resz;
                        tempResults.WeightPSCE[0] = weights[0];
                      }
                    else
                      continue;
                  }
              }

              const double p_value_pid = whoDoFit_pid == 1
                                         ? Fitter_pid->getPVal(fitTrack_pid.get(), rep_pid)
                                         : Fitter_rescue->getPVal(fitTrack_pid.get(), rep_pid); // TMath::Prob(chi2,ndf);
              const double p_value2_pid = 1. - p_value_pid;

              tempResults.pvalue = p_value_pid;

              LocalHisto.h_pv->Fill(p_value2_pid);
              LocalHisto.h_chi2->Fill(chi2_pid);


#ifdef DEBUG_KALMAN2
              att._logger->debug(" / 1-p_value ={}", 1 - p_value_pid);
#endif

              if(kfsop_pid->getCharge() > 0)
                {
                  LocalHisto.hd_chi[0]->Fill(chi2_pid / ndf_pid);
                  LocalHisto.hd_pv[0]->Fill(p_value2_pid);
                }
              else
                {
                  LocalHisto.hd_chi[1]->Fill(chi2_pid / ndf_pid);
                  LocalHisto.hd_pv[1]->Fill(p_value2_pid);
                }

              double Path_time_pid   = -999;
              try
                {
                  Path_time_pid = fitTrack_pid->getTOF(rep_pid);
                }
              catch(genfit::Exception& e)
                {
                  att._logger->info("could not get TOF - 2nd! ");
                  att._logger->info(e.what());
                  LocalHisto.h_stats->Fill("Exc:TOF-2", 1.);
                  continue;
                }

              LocalHisto.h_Path->Fill(Path_length_pid);
              LocalHisto.h_MeanPath->Fill(Path_lengthMean_pid);

              double dE_pid                      = track_stateLast.Eloss;
              double time_fromFirstToLastHit_pid = track_stateLast.time - track_state.time;

              double beta_FirstToLast_pid = 1. / 29.9792458 * Path_length_pid / time_fromFirstToLastHit_pid;
              double mass2_pid = p * sqrt(1. / (beta_FirstToLast_pid * beta_FirstToLast_pid) - 1.);

              if(beta_pid < 0.)
                {
                  att._logger->error("E> Beta Negative - 2nd! {} {} | {}", Path_lengthMean_pid, time_of_flight_pid, Path_length_pid);
                  att._logger->info("{} {} {} {} {}", kfsop_pid->getCharge(), p_pid, chi2_pid, ndf_pid, p_value2_pid);
                }

    #ifdef DEBUG_KALMAN2
              att._logger->debug("charge: {} / time = {} / Path = {} PathMean = {} / Beta = {} / Mass = {}",
                                  kfsop_pid->getCharge(), time_of_flight_pid, Path_length_pid, Path_lengthMean_pid,
                                  beta_pid, mass_pid);
    #endif

              LocalHisto.h_beta->Fill(beta_pid);
              LocalHisto.h_Mass_All->Fill(mass_pid);
              LocalHisto.h_Mass_charge_All->Fill(mass_pid, kfsop_pid->getCharge());
              LocalHisto.h_beta_mom->Fill(p_pid, beta_pid);
              LocalHisto.h_beta_momcharge->Fill(p_pid*kfsop_pid->getCharge(), beta_pid);
              LocalHisto.h_path_tof->Fill(Path_lengthMean_pid , time_of_flight_pid);

              LocalHisto.h_pv_mom->Fill(p_pid, p_value2_pid);
              LocalHisto.h_pv_beta->Fill(beta_pid, p_value2_pid);
              LocalHisto.h_pv_mass->Fill(mass_pid, p_value2_pid);

              // Forward
              const TVectorD& referenceState = stateRefOrig.getState();
              const TVectorD& stateFit_pid = kfsop_pid->getState();

              const TMatrixDSym& covPull_pid = kfsop_pid->getCov();

              const std::vector<int> hist_to_pdg = {2212, -211, 211, -321, 321};
              for(size_t it_pdg = 0; it_pdg < hist_to_pdg.size(); ++it_pdg)
                {
                  if(hist_to_pdg[it_pdg] == PDG_guess_pid)
                    {
                      LocalHisto.h_mom_res[it_pdg]->Fill(momM_pid.Mag(), (momM_pid.Mag() - p3_pid.Mag()) / momM_pid.Mag());
                      LocalHisto.h_ResPull[it_pdg][0]->Fill((kfsop_pid->getCharge() / stateFit_pid[0] - momM_pid.Mag()));
                      LocalHisto.h_ResPull_normal[it_pdg][0]->Fill(momM_pid.Mag(), (kfsop_pid->getCharge() / stateFit_pid[0] - momM_pid.Mag()) /
                                                                                     momM_pid.Mag());
                      for(int i_Res = 1; i_Res < 5; ++i_Res)
                        {
                          LocalHisto.h_ResPull[it_pdg][i_Res]->Fill((stateFit_pid[i_Res] - referenceState[i_Res]));
                          LocalHisto.h_ResPull_normal[it_pdg][i_Res]->Fill(
                              referenceState[i_Res], (stateFit_pid[i_Res] - referenceState[i_Res]) / referenceState[i_Res]);
                        }
                      for(int i_Pull = 0; i_Pull < 5; ++i_Pull)
                        LocalHisto.h_ResPull[it_pdg][i_Pull + 5]->Fill((stateFit_pid[i_Pull] - referenceState[i_Pull]) /
                                                                       sqrt(covPull_pid[i_Pull][i_Pull]));
                    }
                }

              if(p_value2_pid < .75)
                {
                  LocalHisto.h_beta2->Fill(beta_pid);
                  LocalHisto.h_Mass_All2->Fill(mass_pid);
                  LocalHisto.h_beta_mom2->Fill(p_pid, beta_pid);
                  LocalHisto.h_beta_momcharge2->Fill(p_pid*kfsop_pid->getCharge(), beta_pid);
                  LocalHisto.h_Mass_charge_All2->Fill(mass_pid, kfsop_pid->getCharge());
                }
              if(p_value2_pid < .1)
                {
                  LocalHisto.h_beta3->Fill(beta_pid);
                  LocalHisto.h_Mass_All3->Fill(mass_pid);
                  LocalHisto.h_beta_mom3->Fill(p_pid, beta_pid);
                  LocalHisto.h_beta_momcharge3->Fill(p_pid*kfsop_pid->getCharge(), beta_pid);
                  LocalHisto.h_Mass_charge_All3->Fill(mass_pid, kfsop_pid->getCharge());
                  LocalHisto.h_mom_tof_cut->Fill(p_pid, time_of_flight_pid);
                  LocalHisto.h_path_mom_cut->Fill(Path_lengthMean_pid / 30., p_pid);
                  LocalHisto.h_path_tof_cut->Fill(Path_lengthMean_pid / 30., time_of_flight_pid);
                }

              tempResults.firstHit     = id_firstDet;
              tempResults.lastHit      = id_lastDet;
              tempResults.Ncentral     = n_Central;
              tempResults.Nmfiber      = n_MiniFiber;
              tempResults.dE           = dE_pid;
              tempResults.path_time    = Path_time_pid;
              tempResults.path_length  = Path_lengthMean_pid;
              tempResults.path_length2 = Path_length_pid;
              tempResults.tof          = time_of_flight_pid;
              tempResults.tof2         = time_fromFirstToLastHit_pid;
              tempResults.beta         = beta_pid;
              tempResults.beta2        = beta_FirstToLast_pid;
              tempResults.mass         = mass_pid;
              tempResults.mass2        = mass2_pid;
              tempResults.pdg_guess    = PDG_guess_pid;

              // double m_range[4]  = {0.9383, 3.72738, 0.1396, 2.809};
              double m_charge[4] = {1., 2., -1., 2.};

              for(int i = 0; i < 4; i++)
                if(TMath::Abs(kfsop_pid->getCharge() - m_charge[i]) < 0.1)
                  {
                    LocalHisto.h_Mass[i]->Fill(mass_pid);
                    LocalHisto.h_chi2_particle[i]->Fill(chi2_pid);
                    LocalHisto.h_pv_particle[i]->Fill(p_value2_pid);
                  }

    #ifdef DEBUG_KALMAN
              att._logger->debug("Temp :{} {} {}", ntrack, p_value2_pid, mass_pid);
              att._logger->debug("MomRef :{} {} {}", id_track, p_value_pid, mass_pid);
              att._logger->debug("Kalman pval {} charge {}", p_value2_pid, kfsop_pid->getCharge());

    #endif

              RecoEvent.DAF_results.insert(std::make_pair(id_track, tempResults));

            } // if second fit converged
        } // if first fit converged
    } //loop over tracks
  
  return 0;
}


template<class Out>
int TKalmanFilter_DAF_PID<Out>::ProbPIDAssign_Pos(double momenta, double beta)
{
  double prob_beta_Pi = (*(CutPID[211]))(momenta,beta);
  double diff_beta_Pi = (*(CutPID[211])).GetDiff(momenta,beta);
  double chi2_beta_Pi = (*(CutPID[211])).GetChi2(momenta,beta);

  double prob_beta_K = (*(CutPID[321]))(momenta,beta);
  double diff_beta_K = (*(CutPID[321])).GetDiff(momenta,beta);
  double chi2_beta_K = (*(CutPID[321])).GetChi2(momenta,beta);

  double prob_beta_Pr = (*(CutPID[2212]))(momenta,beta);
  double diff_beta_Pr = (*(CutPID[2212])).GetDiff(momenta,beta);
  double chi2_beta_Pr = (*(CutPID[2212])).GetChi2(momenta,beta);

  double prob_beta_Limit = (*(CutPID[991]))(momenta,beta);
  double diff_beta_Limit = (*(CutPID[991])).GetDiff(momenta,beta);
  double chi2_beta_Limit = (*(CutPID[991])).GetChi2(momenta,beta);

  std::vector<int> indexToPdg = {211,321,2212,991};

  double chi2_beta_array[4] = {chi2_beta_Pi, chi2_beta_K, chi2_beta_Pr, chi2_beta_Limit};
  double diff_beta_array[4] = {diff_beta_Pi, diff_beta_K, diff_beta_Pr, diff_beta_Limit};
  double prob_beta_array[4] = {prob_beta_Pi, prob_beta_K, prob_beta_Pr, prob_beta_Limit};
  //std::cout<<" Pi:"<<prob_beta_Pi<<" K+:"<<prob_beta_K<<" p:"<<prob_beta_Pr;

  double Tot_prob = prob_beta_Pi + prob_beta_K + prob_beta_Pr + prob_beta_Limit;
  //std::cout<<" / tot :"<<Tot_prob<<std::endl;
  double prob_pid_array[4];
  prob_pid_array[0] = prob_beta_Pi / Tot_prob;
  prob_pid_array[1]= prob_beta_K / Tot_prob;
  prob_pid_array[2] = prob_beta_Pr / Tot_prob;
  prob_pid_array[3] = prob_beta_Limit / Tot_prob;

  int id_max_pid = TMath::LocMax(4,prob_pid_array);

  int pid_new = indexToPdg[id_max_pid];
  //double prob_pid = prob_pid_array[id_max_pid];
  double prob_beta = prob_beta_array[id_max_pid];
  //double chi2_beta = chi2_beta_array[id_max_pid];
  //double diff_beta = diff_beta_array[id_max_pid];

  if(prob_beta < att.PID_minProb) pid_new = 0;

  return pid_new;
}


template<class Out>
int TKalmanFilter_DAF_PID<Out>::ProbPIDAssign_Neg(double momenta, double beta)
{
  double prob_beta_Pi = (*(CutPID[-211]))(momenta,beta);
  double diff_beta_Pi = (*(CutPID[-211])).GetDiff(momenta,beta);
  double chi2_beta_Pi = (*(CutPID[-211])).GetChi2(momenta,beta);

  double prob_beta_Limit = (*(CutPID[-991]))(momenta,beta);
  double diff_beta_Limit = (*(CutPID[-991])).GetDiff(momenta,beta);
  double chi2_beta_Limit = (*(CutPID[-991])).GetChi2(momenta,beta);

  std::vector<int> indexToPdg = {-211,-991};

  double chi2_beta_array[2] = {chi2_beta_Pi, chi2_beta_Limit};
  double diff_beta_array[2] = {diff_beta_Pi, diff_beta_Limit};
  double prob_beta_array[2] = {prob_beta_Pi, prob_beta_Limit};
  //std::cout<<" Pi:"<<prob_beta_Pi<;

  double Tot_prob = prob_beta_Pi + prob_beta_Limit ;
  //std::cout<<" / tot :"<<Tot_prob<<std::endl;
  double prob_pid_array[2];
  prob_pid_array[0] = prob_beta_Pi / Tot_prob;
  prob_pid_array[1]= prob_beta_Limit / Tot_prob;

  int id_max_pid = TMath::LocMax(2,prob_pid_array);

  int pid_new = indexToPdg[id_max_pid];
  //double prob_pid = prob_pid_array[id_max_pid];
  double prob_beta = prob_beta_array[id_max_pid];
  //double chi2_beta = chi2_beta_array[id_max_pid];
  //double diff_beta = diff_beta_array[id_max_pid];

  if(prob_beta < att.PID_minProb) pid_new = 0;

  return pid_new;
}


template class TKalmanFilter_DAF_PID<MCAnaEventG4Sol>;
template class TKalmanFilter_DAF_PID<Ana_WasaEvent>;
