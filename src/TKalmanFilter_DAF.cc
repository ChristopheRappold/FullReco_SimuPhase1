#include "TKalmanFilter_DAF.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"
#include "FullRecoEvent.hh"
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
TKalmanFilter_DAF<Out>::TKalmanFilter_DAF(const THyphiAttributes& attribut)
    : TDataProcessInterface<Out>("mom_fit_kalman"), att(attribut)
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

  Fitter = nullptr;

  if(att.KF_Kalman || att.KF_KalmanSqrt)
    Fitter = new genfit::KalmanFitter(10, dPVal, 1e3, att.KF_KalmanSqrt); //
  else if(att.KF_KalmanRef)
    Fitter = new genfit::KalmanFitterRefTrack(nIter, dPVal); //,1e4);
  //Fitter->setMultipleMeasurementHandling(mmHandling);
  else if(att.KF_DAF || att.KF_DAFRef)
    Fitter = new genfit::DAF(att.KF_DAFRef);
  else
    att._logger->error("E> Kalman Filter fitter not set correctly !");

  if(att.KF_G4e)
    {
      att._logger->error("!> Start KF_G4e !");
      std::string nameField = att.Wasa_FieldMapName;
      double signDir = att.Wasa_Side ? 1.0 : -1.0;
      std::array<double,3> fieldValues = {att.Field_Strength*signDir*0.1, 0., 0.}; 
      // std::unordered_map<std::string, std::string> G4eConfig = {{"FullProp","0"},{"G4eBasf2PhysicsList","1"},{"ExactJac","1"},
      // 								//,{"TrackingVerbose","10"},{"ControlVerbose","2"},{"G4eVerbose","4"}};
      // 								//{"StepLength","25 mm"},
      // 								//{"MagFieldDiff","0.01",},
      // 								{"MaxEnergyLoss","0.0005"}

      std::unordered_map<std::string, std::string> G4eConfig = {{"FullProp",att.G4e_FullProp},
								{"G4eBasf2PhysicsList",att.G4e_Basf2List},
								{"ExactJac",att.G4e_ExactJac},
								//,{"TrackingVerbose","10"},{"ControlVerbose","2"},
								{"G4eVerbose",att.G4e_Verbose},
								//{"StepLength","25 mm"},
								//{"MagFieldDiff","0.01",},
								{"MaxEnergyLoss", att.G4e_MaxEnergyLoss}};

      G4eMag = std::make_unique<G4eManager>(G4eConfig,fieldValues, true, nameField);
      G4eMag->InitG4e();
    }


  
  Fitter->setMinIterations(3);
  Fitter->setMaxIterations(nIter);

  Nb_CentralCut = att.KF_NbCentralCut;
  Nb_MiniFiberCut = att.KF_NbMiniFiberCut;

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
TKalmanFilter_DAF<Out>::~TKalmanFilter_DAF()
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
void TKalmanFilter_DAF<Out>::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }

template<class Out>
ReturnRes::InfoM TKalmanFilter_DAF<Out>::operator()(FullRecoEvent& RecoEvent, Out* OutTree)
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
void TKalmanFilter_DAF<Out>::SelectHists()
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
int TKalmanFilter_DAF<Out>::Exec(FullRecoEvent& RecoEvent, Out* OutTree)
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
      auto PDG_particle     = TDatabasePDG::Instance()->GetParticle(FitRes.pdg_guess);
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
ReturnRes::InfoM TKalmanFilter_DAF<Out>::SoftExit(int ) { return ReturnRes::Fine; }

template<class Out>
int TKalmanFilter_DAF<Out>::Kalman_Filter_FromTrack(FullRecoEvent& RecoEvent)
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
      auto it_ListHitsSim = RecoEvent.TrackDAFSim.find(id_track);

      auto getZpos = [](genfit::AbsMeasurement* m) {
        TVectorD& HitrawRef = m->getRawHitCoords();
        if(HitrawRef.GetNrows() == 2 || HitrawRef.GetNrows() == 1)
          {
            genfit::StateOnPlane dummy;
            genfit::SharedPlanePtr plane = dynamic_cast<genfit::PlanarMeasurement*>(m)->constructPlane(dummy);
            TVector3 tempO(plane->getO());
            return tempO.Z();
          }
        else if(HitrawRef.GetNrows() == 3)
          return HitrawRef[2];
	else if(HitrawRef.GetNrows() == 7)
	  {
	    return 0.5*(HitrawRef[2]+HitrawRef[5]);
	  }
        else
          {
            fmt::print("E> rawref not proper ! {}", HitrawRef.GetNrows());
            HitrawRef.Print();
            return -999.;
          }
      };
      std::set<std::tuple<double, int, int, double> > id_dets;
      double total_dE = 0.;
      int n_Central   = 0;
      int n_MiniFiber = 0;
      // std::cout << "it_ListHits->second.size : " << it_ListHits->second.size() << std::endl;
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
          if( (id_det>=G4Sol::MiniFiberD1_x1 && id_det<=G4Sol::MiniFiberD1_v2) || (id_det>=G4Sol::MiniFiberD1_x && id_det<= G4Sol::MiniFiberD2_u))
            ++n_MiniFiber;
          genfit::AbsMeasurement* currentHit = RecoEvent.ListHits[id_det][id_hit].get();

          total_dE += it_trackInfo.second[id_det].Eloss;
	  //std::cout << "id_det : " << id_det << std::endl;
	  //std::cout << "id_hit : " << id_hit << std::endl;
	  //id_dets.insert(std::make_tuple(getZpos(currentHit), id_det, id_hit));
          id_dets.insert(std::make_tuple(id_det, id_det, id_hit, getZpos(currentHit)));
          // std::cout << " id_det : " << id_det << std::endl;
        }
      // std::cout << "n_Central : " << n_Central << std::endl;
      // std::cout << "id_dets : " << id_dets.size() << std::endl;

      // if(id_firstDet != std::get<1>(*firstHit) )
      //  {
      //    std::cout<<"!> Kalman: firstHit not matched !"<<id_firstDet<<" "<<std::get<1>(*firstHit)<<"\n";
      //    auto Decay = RecoEvent.TrackMother.find(id_track);
      //    bool IsDecay = Decay != RecoEvent.TrackMother.end();
      //    std::cout<<"track #"<<id_track<<" "<<PDG<<" "<<IsDecay<<"\n";
      //    for(auto ii : id_dets)
      //      std::cout<<"| :"<<std::get<0>(ii)<<" "<<std::get<1>(ii)<<" "<<std::get<2>(ii)<<"\n";
      //    std::cout<<" +---\n";
      //  }

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
		{
		  att._logger->debug("!> #{} : {} {} {}",i,it_trackInfo.second[i].pdg, it_trackInfo.second[i].mass);
		}
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
      // TVectorD& tempHitrawRef         = tempHit->getRawHitCoords();
      // std::cout << "tempHitrawRef 0 : " << tempHitrawRef[0] << std::endl;
      // std::cout << "tempHitrawRef 1 : " << tempHitrawRef[1] << std::endl;
      // const TVector3 firstPos(tempHitrawRef[0], tempHitrawRef[1], getZpos(tempHit));
      //const TVector3 firstPos(0, 0, 0);
      auto it_hitFirstSim = it_ListHitsSim->second[id_firstDet][0];
      const TVector3 firstPos(it_hitFirstSim.hitX, it_hitFirstSim.hitY, it_hitFirstSim.hitZ);

      const int PDG     = static_cast<int>(track_state.pdg);
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

      // std::cout<<"id_dets:"<<id_dets.size()<<" ";
      if(charge <= 1)
        {
          std::set<std::tuple<double, int, int, double> > temp_id_dets;
          for(auto id_det : id_dets)
            {
              if(std::get<1>(id_det) <= G4Sol::RPC_h)
                temp_id_dets.insert(id_det);
            }

          std::swap(id_dets, temp_id_dets);
        }

      // cout<<"after :"<<id_dets.size()<<"\n";

      firstHit             = id_dets.cbegin();
      id_firstDet          = std::get<1>(*firstHit);
      const auto lastHit   = id_dets.crbegin();
      const int id_lastDet = std::get<1>(*lastHit); // std::get<0>(lastHit);

#ifdef DEBUG_KALMAN
      for(auto idet : id_dets)
        att._logger->debug("det:{} / {} / [{}] at Z:{}", std::get<1>(idet),G4Sol::nameLiteralDet.begin()[std::get<1>(idet)],  std::get<2>(idet), std::get<3>(idet));

      const int PID = charge * charge;
      att._logger->debug(" PID:{} {} {}", PID, charge, PDG);
#endif

      const InfoPar track_stateLast(it_trackInfo.second[id_lastDet]);

      // Forward
      // TVector3 init_p(track_state.momX, track_state.momY, track_state.momZ);
      // init_p.SetMag(0.5);
      auto it_init   = RecoEvent.TrackDAFInit.find(id_track);
      double init_px = gRandom->Gaus(it_init->second.momX, it_init->second.momX * 0.2);
      double init_py = gRandom->Gaus(it_init->second.momY, it_init->second.momY * 0.2);
      double init_pz = gRandom->Gaus(it_init->second.momZ, it_init->second.momZ * 0.2);
      TVector3 init_p(init_px, init_py, init_pz);
      // TVector3 init_p(track_stateLast.momX, track_stateLast.momY, track_stateLast.momZ);

      double seed_Mom_Mag = init_p.Mag();
      if(TMath::Abs(seed_Mom_Mag) < 1e-9)
        {
	  //#ifdef DEBUG_KALMAN
          att._logger->debug("!> Seed Momemtum with TVector3 is zero ! correcting ");
	  //#endif
          auto tempLastHit = id_dets.crbegin();
          ++tempLastHit;
          auto lastHit2 = tempLastHit;
          const InfoPar track_stateBeforeLast(it_trackInfo.second[std::get<1>(*lastHit2)]);
          const TVector3 init_p2(track_stateBeforeLast.momX, track_stateBeforeLast.momY, track_stateBeforeLast.momZ);

          if(TMath::Abs(init_p2.Mag()) < 1e-7)
            {
              att._logger->error("E> Seed Momemtum with TVector3 is zero ! {} Mom:{}", PDG_particle->GetName(),
                                 seed_Mom_Mag);
              att._logger->error("TrackID #{} hit_id :", it_ListHits->first);
              std::vector<std::stringstream> s1(it_ListHits->second.size() / 20 + 1);
              std::vector<std::stringstream> s2(it_ListHits->second.size() / 20 + 1);
              for(size_t i = 0; i < it_ListHits->second.size(); ++i)
                {
                  s1[i / 20] << printW(i, 3) << ", ";
                  s2[i / 20] << printW(it_ListHits->second[i], 3) << ", ";
                }
              for(size_t i = 0; i < s1.size(); ++i)
                {
                  att._logger->error("idDet:{}", s1[i].str());
                  att._logger->error("stat :{}", s2[i].str());
                }

              att._logger->error("Track Info :");
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
                  att._logger->error("idDet:{}", s11[i].str());
                  att._logger->error("pdg  :{}", s22[i].str());
                  att._logger->error("momX :{}", s33[i].str());
                  att._logger->error("momY :{}", s44[i].str());
                  att._logger->error("momZ :{}", s55[i].str());
                  att._logger->error("time :{}", s66[i].str());
                }

              // for(auto id_hit : track.second)
              //  std::cout<<" "<<id_hit<<", ";
              // std::cout<<"] "<<std::endl;
              continue;
            }
          else
            {
              init_p.SetXYZ(init_p2.X(), init_p2.Y(), init_p2.Z());
              seed_Mom_Mag = init_p.Mag();
#ifdef DEBUG_KALMAN
              att._logger->debug("!> Reset init_p N#{}", Nb_event);

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
              for(auto idet : id_dets)
                {
                  att._logger->debug("det:{} [{}] at Z:{}", std::get<1>(idet), std::get<2>(idet), std::get<0>(idet));
                }
#endif
            }
        }
#ifdef DEBUG_KALMAN
      att._logger->debug("init mom : ");
      init_p.Print();
#endif

      const double mom_res = .0500;
      double new_P         = gRandom->Gaus(seed_Mom_Mag, mom_res * seed_Mom_Mag);
      TVector3 seed_p(init_p);
      seed_p.SetMag(new_P);

      // Forward
      genfit::AbsTrackRep* REP = att.KF_G4e ? new genfit::G4eTrackRep(*G4eMag,PDG) : new genfit::RKTrackRep(PDG, 1);
      // REP->setDebugLvl(2);
      // genfit::AbsTrackRep* REP = new genfit::RKTrackRep(PDG,-1);

      // std::cout<<" -- :"<<id_firstDet<<" "<<std::get<2>(*firstHit)<<" \n";
      // std::cout<<" -- :"<<id_lastDet<<" "<<std::get<2>(*lastHit)<<" \n";

      // Forward
      //genfit::AbsMeasurement* FirstHit = RecoEvent.ListHits[id_firstDet][std::get<2>(*firstHit)].get();
      //TVectorD& HitrawRef              = FirstHit->getRawHitCoords();

      //const TVector3 init_point(HitrawRef[0], HitrawRef[1], getZpos(FirstHit));

      double init_x = gRandom->Gaus(it_init->second.posX, 1.0);
      double init_y = gRandom->Gaus(it_init->second.posY, 1.0);
      double init_z = gRandom->Gaus(it_init->second.posZ, 3.0);
      const TVector3 init_point(init_x, init_y, init_z);
      //const TVector3 init_point(it_hitFirstSim.hitX, it_hitFirstSim.hitY, it_hitFirstSim.hitZ);

      // TVector3 init_point(track_state[9],track_state[10],track_state[11]);

      TMatrixDSym CovM(6);
      // double resolution = 0.1;
      for(int i = 0; i < 6; ++i)
        CovM(i, i) = 1000.;

      ////Backward
      // genfit::MeasuredStateOnPlane stateToFinalRef(REP);
      // TVector3 momFinalRef(track_state.momX,track_state.momY,track_state.momZ);
      // if(momFinalRef.Mag()<1e-8)
      //  att._logger->error("E> momFinalRef null !{} {} {}", track_state.momX, track_state.momY, track_state.momZ);
      // REP->setPosMomCov(stateToFinalRef, firstPos, momFinalRef, CovM);

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

        // std::auto_ptr<genfit::Track> fitTrack(new genfit::Track(REP, seedState, seedCov));
#ifndef DISPLAY
      std::unique_ptr<genfit::Track> fitTrack = std::make_unique<genfit::Track>(REP, seedState, seedCov);
#else
      std::shared_ptr<genfit::Track> fitTrack = std::make_shared<genfit::Track>(REP, seedState, seedCov);
#endif

      // std::tuple<int, double> lastHit = std::make_tuple(-1, -9999.);

      //// Backward
      // std::set<std::tuple<double,int,int> > temp_id_det2;
      // for(auto ids = id_dets.crbegin(), ids_end = id_dets.crend(); ids != ids_end ; ++ids)
      //  temp_id_det2.insert(*ids);
      // std::swap(id_dets,temp_id_det2);

      for(auto ids : id_dets)
        {
          int id_det = std::get<1>(ids);
          int id_hit = std::get<2>(ids);

          genfit::AbsMeasurement* currentHit = RecoEvent.ListHits[id_det][id_hit].get();
#ifdef DEBUG_KALMAN
          att._logger->debug("Loop insertPoint: #det:{} #hit {} {}", G4Sol::nameLiteralDet.begin()[id_det], id_hit,
                             fmt::ptr(currentHit));
          currentHit->Print();
#endif
          // if(std::get<1>(lastHit) < currentHit->getRawHitCoords()(2))
          //   {
          //     lastHit = std::make_tuple(id_det, currentHit->getRawHitCoords()(2));
          //   }

          // Vtracks->insertPoint(new genfit::TrackPoint(RecoEvent.ListHitsDAF[id_track][ii],Vtracks));
          fitTrack->insertPoint(new genfit::TrackPoint(currentHit, fitTrack.get(), false));
          // id_lastDet = id_det;
        }

      // InfoPar track_stateLast(it_trackInfo.second[id_lastDet]);

      // assert(fitTrack->checkConsistency());
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
                  // REP->extrapolateToPlane(*kfsop, stateToFinalRef.getPlane());
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

              ResSolDAF tempResults;
              tempResults.charge    = charge;
              tempResults.pdg_ini   = PDG;
              tempResults.pdg_guess = PDG;
              // Forward
              const TVectorD& referenceState = stateRefOrig.getState();
              // const TVectorD& referenceState = stateToFinalRef.getState();
              const TVectorD& stateFit = kfsop->getState();

              const TMatrixDSym& covPull = kfsop->getCov();

              TVector3 p3, posRef;
              TMatrixDSym covFit(6);

              kfsop->getPosMomCov(posRef, p3, covFit);
              // std::cout << Form("pos pi  : %.2f ,  %.2f , %.2f",posRef.X(), posRef.Y(), posRef.Z()) << std::endl;

              // Forward
              tempResults.momX_init = init_p.X();
              tempResults.momY_init = init_p.Y();
              tempResults.momZ_init = init_p.Z();
              // tempResults.momX_init = track_state.momX;
              // tempResults.momY_init = track_state.momY;
              // tempResults.momZ_init = track_state.momZ;

              tempResults.momX = p3.X();
              tempResults.momY = p3.Y();
              tempResults.momZ = p3.Z();

              tempResults.posX = posRef.X();
              tempResults.posY = posRef.Y();
              tempResults.posZ = posRef.Z();

              for(int row = 0; row < 6; ++row)
                for(int col = 0; col < 6; ++col)
                  tempResults.cov_matrix[row][col] = covFit(row, col);

              const double chi2 = whoDoFit == 1
                                      ? Fitter->getChiSqu(fitTrack.get(), REP)
                                      : Fitter_rescue->getChiSqu(fitTrack.get(), REP); // Vtracks->getChiSqu();

              const double ndf =
                  whoDoFit == 1 ? Fitter->getNdf(fitTrack.get(), REP)
                                : Fitter_rescue->getNdf(fitTrack.get(), REP); // static_cast<int>(Vtracks->getNDF());

              tempResults.chi2    = chi2;
              tempResults.ndf     = ndf;
              tempResults.fitter  = whoDoFit;
              tempResults.iterNum = (int)(fitTrack->getKalmanFitStatus()->getNumIterations());
#ifdef DEBUG_KALMAN2
              att._logger->debug(" / chi2 ={} / ndf ={}", chi2, ndf);
#endif

              {
                unsigned int np = fitTrack->getNumPointsWithMeasurement();
                //std::cout << "np : " << np << std::endl;
                for(unsigned int i = 0; i < np; ++i)
		  {
                    genfit::TrackPoint* tp_tmp    = fitTrack->getPointWithMeasurementAndFitterInfo(i);//, REP);
		    if(tp_tmp==nullptr)
		      {
			att._logger->debug("E> TrackPoint null ! #i {} / {}\n",i,np);
			continue;
		      }
		    genfit::KalmanFitterInfo* kfi = static_cast<genfit::KalmanFitterInfo*>(tp_tmp->getFitterInfo(REP));
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
                            // double weight = residual.getWeight();
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
                    else if(id_det>=G4Sol::MiniFiberD1_x && id_det<=G4Sol::MiniFiberD2_u){
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
                        // double resz  = -999;
                        // std::cout << "res : " << res << std::endl;
                        LocalHisto.h_ResPSCE[0]->Fill(res);
                        LocalHisto.h_ResPSCE[1]->Fill(resz);
                        tempResults.ResPSCE[0]    = res;
                        tempResults.ResPSCE[1]    = resz;
                        tempResults.WeightPSCE[0] = weights[0];
                      }
                    else
                      {
                        continue;
                      }
                  }
              }

              // double p_value = TMath::Prob(chi2,ndf);
              // double p_value = MathKalman::Prob()(chi2,ndf);
              // double p_value2= 1.-p_value;
              const double p_value = whoDoFit == 1
                                         ? Fitter->getPVal(fitTrack.get(), REP)
                                         : Fitter_rescue->getPVal(fitTrack.get(), REP); // TMath::Prob(chi2,ndf);
              const double p_value2 = 1. - p_value;

              tempResults.pvalue = p_value;

              LocalHisto.h_pv->Fill(p_value2);
              LocalHisto.h_chi2->Fill(chi2);

              // h_chi2_smooth->Fill(Vtracks[itr]->getChiSquSmooth());

#ifdef DEBUG_KALMAN2
              att._logger->debug(" / 1-p_value ={}", 1 - p_value);
#endif

              // n_status[static_cast<int>(charge)+1]++;

              if(charge > 0)
                {
                  LocalHisto.hd_chi[0]->Fill(chi2 / ndf);
                  LocalHisto.hd_pv[0]->Fill(p_value2);
                }
              else
                {
                  LocalHisto.hd_chi[1]->Fill(chi2 / ndf);
                  LocalHisto.hd_pv[1]->Fill(p_value2);
                }

              // GFDetPlane InitPlane(RecoEvent.ListHitsDAF[id_track][TR1X][0]->getDetPlane(rep));
              // TVector3 TR1x_Pos(Vtracks->getPos(InitPlane));

              // double DistToTR1= 0;
              // try
              //   {
              //     const genfit::MeasuredStateOnPlane
              //     TR1PlaneState(fitTrack->getFittedState(find_tr1->second,REP,true)); TVector3
              //     TR1x_Pos(TR1PlaneState.getPos()); HitTR1 = TR1x_Pos; DistToTR1=(TR1x_Pos-Plane_time).Mag();
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
              double Path_time   = -999;
              try
                {
                  Path_length = fitTrack->getTrackLen(REP); // getTotalLength();
                  // Path_length = -fitTrack->getTrackLen(REP); // getTotalLength();
                  // double Path_length = Vtracks->getTotalLength();
                  // double Path_lengthB = Vtracks->getTotalLengthBack();
                }
              catch(genfit::Exception& e)
                {
                  att._logger->info("could not get TrackLen ! ");
                  att._logger->info(e.what());
                  LocalHisto.h_stats->Fill("Exc:TrackLen", 1.);
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
                  att._logger->info("could not get TOF! ");
                  att._logger->info(e.what());
                  LocalHisto.h_stats->Fill("Exc:TOF", 1.);
                  continue;
                }

              // double Path_timeB = Vtracks->getTotalTimeBack();
              // double Path_lengthMean = Path_length + Path_lengthB;
              double Path_lengthMean = Path_length;
              // Path_lengthMean/=2.;
              Path_lengthMean += TMath::Sqrt(TMath::Sq(firstPos.X()-att.Target_PositionX)+TMath::Sq(firstPos.Y()-att.Target_PositionY)+TMath::Sq(firstPos.Z()-att.Target_PositionZ)); // DistToTR1;
                                                 // 	  if(ndf==1)
                                                 // 	    Path_lengthMean=Path_length;

              LocalHisto.h_Path->Fill(Path_length);
              // LocalHisto.h_Path_Back->Fill(Path_lengthB);
              LocalHisto.h_MeanPath->Fill(Path_lengthMean);
              // LocalHisto.h_dpath->Fill(Path_lengthB-Path_length);

              // double rapidity = 0.5 * TMath::Log(( TMath::Sqrt(rep_length->getMass()*rep_length->getMass() +
              // p3.Mag2()) +p3.Z())/( TMath::Sqrt(rep_length->getMass()*rep_length->getMass() + p3.Mag2()) - p3.Z()));
              // for(unsigned int vol = 0 ; vol < att.name_GeoVolumes.size() ; ++vol)
              //   {
              //     LocalHisto.Material_XX0_y[vol]->Fill(rapidity,rep_length->getXX0(att.name_GeoVolumes[vol]));
              //     LocalHisto.Material_dE_y[vol]->Fill(rapidity,rep_length->getDE(att.name_GeoVolumes[vol]));
              //   }

#ifdef VERBOSE_EVE

              att._logger->debug("$$$$$$$$$$$$-----------------------$$$$$$$$$$$$$$$");
              att._logger->debug("DAF Length :ListTrackPoints :");
              att._logger->debug("TEveLine* EveTrack_DAF1 = new TEveLine(TrackRepList_DAF,{});",
                                 rep_length->EveListTrackPoint.size());

              for(unsigned int i = 0; i < rep_length->EveListTrackPoint.size(); ++i)
                {
                  att._logger->debug("EveTrack_DAF1->SetPoint({}, {:0.10f},{:0.10f},{:0.10f});", i,
                                     rep_length->EveListTrackPoint[i].X(), rep_length->EveListTrackPoint[i].Y(),
                                     rep_length->EveListTrackPoint[i].Z());
                }

              att._logger->debug("Fitted Hit :");
              att._logger->debug("TEvePointSet* EveHitFitter1 = new TEvePointSet(HitFitter1,5)");
              att._logger->debug("EveHitFitter1->SetPoint({}, {:0.10f},{:0.10f},{:0.10f});", 0, HitTR1.X(), HitTR1.Y(),
                                 HitTR1.Z());
              att._logger->debug("EveHitFitter1->SetPoint({}, {:0.10f},{:0.10f},{:0.10f});", 1, HitDC1.X(), HitDC1.Y(),
                                 HitDC1.Z());
              att._logger->debug("EveHitFitter1->SetPoint({}, {:0.10f},{:0.10f},{:0.10f});", 2, HitTR2.X(), HitTR2.Y(),
                                 HitTR2.Z());
              att._logger->debug("EveHitFitter1->SetPoint({}, {:0.10f},{:0.10f},{:0.10f});", 3, HitDC2.X(), HitDC2.Y(),
                                 HitDC2.Z());
              att._logger->debug("EveHitFitter1->SetPoint({}, {:0.10f},{:0.10f},{:0.10f});", 4, HitTOF.X(), HitTOF.Y(),
                                 HitTOF.Z());

              att._logger->debug("DAF Track :");
              att._logger->debug("TEveRecTrack *DAF_track_rc = new TEveRecTrack();");
              att._logger->debug("DAF_track_rc->fV.Set({},{},{});", HitTR1.X(), HitTR1.Y(), HitTR1.Z());
              att._logger->debug("DAF_track_rc->fP.Set({},{},{});", p3.X(), p3.Y(), p3.Z());
              // std::cout<<"int sign = -"<<list_MC[current_idMC]->Charge<<std::endl;
              att._logger->debug("DAF_track_rc->fSign = -{};", charge);
              att._logger->debug("$$$$$$$$$$$$-----------------------$$$$$$$$$$$$$$$");

#endif

              const double dE                      = track_stateLast.Eloss;
              const double time_of_flight          = track_stateLast.time; // Vtracks->getTimeOfFlight();//+6.;
              const double time_fromFirstToLastHit = track_stateLast.time - track_state.time;

              const double beta = 1. / 29.9792458 * Path_lengthMean /
                                  time_of_flight; /// beta = dL[m]/(c[m/s]* dT[s]) = dL[mm]/(300[mm/ns] * dT[ns])
              const double beta_FirstToLast = 1. / 29.9792458 * Path_length / time_fromFirstToLastHit;
              // double beta = 1./30.*new_Total_Length/time_of_flight; ///beta = dL[m]/(c[m/s]* dT[s]) =
              // dL[mm]/(300[mm/ns] * dT[ns]) double gamma = 1./sqrt(1.-beta*beta);
              const double mass  = p * sqrt(1. / (beta * beta) - 1.);
              const double mass2 = p * sqrt(1. / (beta_FirstToLast * beta_FirstToLast) - 1.);

              if(beta < 0.)
                {
                  att._logger->error("E> Beta Negative ! {} {} | {}", Path_lengthMean, time_of_flight, Path_length);
                  //<<Path_lengthB<<" "<<endl;
                  att._logger->info("{} {} {} {} {}", charge, p, chi2, ndf, p_value2);
                }

#ifdef DEBUG_KALMAN2
              att._logger->debug("charge: {} / time = {} / Path = {} PathMean = {} / Beta = {} / Mass = {}", charge,
                                 time_of_flight /*Detectors_UTr["TOF"][id_track][3]*/, Path_length, Path_lengthMean,
                                 beta, mass);
#endif

              LocalHisto.h_beta->Fill(beta);
              LocalHisto.h_Mass_All->Fill(mass);
              LocalHisto.h_Mass_charge_All->Fill(mass, charge);
              LocalHisto.h_beta_mom->Fill(p, beta);
              LocalHisto.h_path_tof->Fill(Path_lengthMean / 30., time_of_flight);

              LocalHisto.h_pv_mom->Fill(p, p_value2);
              LocalHisto.h_pv_beta->Fill(beta, p_value2);
              LocalHisto.h_pv_mass->Fill(mass, p_value2);

              const std::vector<int> hist_to_pdg = {2212, -211, 211, -321, 321};
              for(size_t it_pdg = 0; it_pdg < hist_to_pdg.size(); ++it_pdg)
                {
                  if(hist_to_pdg[it_pdg] == PDG)
                    {
                      LocalHisto.h_mom_res[it_pdg]->Fill(init_p.Mag(), (init_p.Mag() - p3.Mag()) / init_p.Mag());

                      LocalHisto.h_ResPull[it_pdg][0]->Fill((charge / stateFit[0] - init_p.Mag()));
                      LocalHisto.h_ResPull_normal[it_pdg][0]->Fill(init_p.Mag(), (charge / stateFit[0] - init_p.Mag()) /
                                                                                     init_p.Mag());
                      for(int i_Res = 1; i_Res < 5; ++i_Res)
                        {
                          LocalHisto.h_ResPull[it_pdg][i_Res]->Fill((stateFit[i_Res] - referenceState[i_Res]));
                          LocalHisto.h_ResPull_normal[it_pdg][i_Res]->Fill(
                              referenceState[i_Res], (stateFit[i_Res] - referenceState[i_Res]) / referenceState[i_Res]);
                        }
                      for(int i_Pull = 0; i_Pull < 5; ++i_Pull)
                        LocalHisto.h_ResPull[it_pdg][i_Pull + 5]->Fill((stateFit[i_Pull] - referenceState[i_Pull]) /
                                                                       sqrt(covPull[i_Pull][i_Pull]));
                    }
                }

              if(p_value2 < .75)
                {
                  LocalHisto.h_beta2->Fill(beta);
                  LocalHisto.h_Mass_All2->Fill(mass);
                  LocalHisto.h_beta_mom2->Fill(p, beta);
                  LocalHisto.h_Mass_charge_All2->Fill(mass, charge);
                }
              if(p_value2 < .1)
                {
                  LocalHisto.h_beta3->Fill(beta);
                  LocalHisto.h_Mass_All3->Fill(mass);
                  LocalHisto.h_beta_mom3->Fill(p, beta);
                  LocalHisto.h_Mass_charge_All3->Fill(mass, charge);
                  LocalHisto.h_mom_tof_cut->Fill(p, time_of_flight);
                  LocalHisto.h_path_mom_cut->Fill(Path_lengthMean / 30., p);
                  LocalHisto.h_path_tof_cut->Fill(Path_lengthMean / 30., time_of_flight);
                }
              tempResults.firstHit     = id_firstDet;
              tempResults.lastHit      = id_lastDet;
              tempResults.Ncentral     = n_Central;
              tempResults.Nmfiber      = n_MiniFiber;
              tempResults.dE           = dE;
              tempResults.path_time    = Path_time;
              tempResults.path_length  = Path_lengthMean;
              tempResults.path_length2 = Path_length;
              tempResults.tof          = time_of_flight;
              tempResults.tof2         = time_fromFirstToLastHit;
              tempResults.beta         = beta;
              tempResults.beta2        = beta_FirstToLast;
              tempResults.mass         = mass;
              tempResults.mass2        = mass2;

              // double m_range[4]  = {0.9383, 3.72738, 0.1396, 2.809};
              double m_charge[4] = {1., 2., -1., 2.};

              for(int i = 0; i < 4; i++)
                // if (TMath::Abs(mass-m_range[i])<0.25*m_range[i] && TMath::Abs(charge-m_charge[i])<0.1)
                if(TMath::Abs(charge - m_charge[i]) < 0.1)
                  {
                    // double dmass = TMath::Abs(mass-m_range[i])/m_range[i];
                    LocalHisto.h_Mass[i]->Fill(mass);
                    LocalHisto.h_chi2_particle[i]->Fill(chi2);
                    LocalHisto.h_pv_particle[i]->Fill(p_value2);
                  }

#ifdef DEBUG_KALMAN
              att._logger->debug("Temp :{} {} {}", ntrack, p_value2, mass);
              att._logger->debug("MomRef :{} {} {}", id_track, p_value, mass);
#endif

              RecoEvent.DAF_results.insert(std::make_pair(id_track, tempResults));
              // RecoEvent.Sigma[ntrack]=CovSigma;
#ifdef DEBUG_KALMAN
              att._logger->debug("Kalman pval {} charge {}", p_value2, charge);
#endif
            }
        }
  
  return 0;
}


template class TKalmanFilter_DAF<MCAnaEventG4Sol>;
template class TKalmanFilter_DAF<Ana_WasaEvent>;
