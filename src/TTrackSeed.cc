#include "TTrackSeed.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"
#include "FullRecoEvent.hh"
#include "KalmanFittedStateOnPlane.h"
#include "KalmanFitterInfo.h"
#include "PlanarMeasurement.h"
#include "ProlateSpacepointMeasurement.h"
#include "ReturnRes.hh"
#include "StateOnPlane.h"

#include <set>
#include <sstream>
#include <tuple>

//#define DEBUG_RIEMANNFINDER

using namespace std;
using namespace G4Sol;

template<class Out>
TTrackSeed<Out>::TTrackSeed(const THyphiAttributes& attribut) : TDataProcessInterface<Out>("TrackSeed"), att(attribut)
{

  LastFrontWall = att.Wasa_Side == 0 ? G4Sol::PSFE : G4Sol::PSBE;

  //res_h.time.resize(1,1);

  att._logger->info("TrackSeed : set ");
}

template<class Out>
TTrackSeed<Out>::~TTrackSeed()
{ }

template<class Out>
void TTrackSeed<Out>::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }

template<class Out>
ReturnRes::InfoM TTrackSeed<Out>::operator()(FullRecoEvent& RecoEvent, Out* OutTree)
{

  int result_finder = Exec(RecoEvent, OutTree);

  return SoftExit(result_finder);
}

template<class Out>
int TTrackSeed<Out>::Exec(FullRecoEvent& RecoEvent, Out* OutTree)
{

  std::vector<Seed::RTrack> newTracksCand;
  int res = TrackFit(RecoEvent, newTracksCand);

  for(int idT = 0; auto& trackC : newTracksCand)
    {
      int charge = trackC.helix.q;
      RecoEvent.TracksFound.emplace_back(trackC.sortedHits,true,charge, trackC.helix.Par,
					 trackC.helix.Cov, trackC.helix.chi2_circle,trackC.helix.chi2_line);

      const auto it_firsthit = trackC.sortedHits.begin();
      auto [_1,_2,id_track] = TempHitToAllHits[*it_firsthit];

      auto it_InitInfo = RecoEvent.TrackDAFInit.find(id_track);
      it_InitInfo->second.charge = charge;
      TVector3 tempMom(it_InitInfo->second.momX,it_InitInfo->second.momY,it_InitInfo->second.momZ);
      tempMom.Unit();
      tempMom.SetMag(trackC.helix.Par[2]/tempMom.Perp());

      it_InitInfo->second.momX = tempMom.X();
      it_InitInfo->second.momY = tempMom.Y();
      it_InitInfo->second.momZ = tempMom.Z();

      TTrackCand* OutTrack = dynamic_cast<TTrackCand*>(OutTree->TrackCand->ConstructedAt(OutTree->TrackCand->GetEntries()));

      OutTrack->FitStatus = trackC.toRefit;
      OutTrack->Charge = charge;

      for(size_t id_par = 0; id_par<trackC.helix.Par.size();++id_par)
	OutTrack->Seed_Par[id_par] = trackC.helix.Par[id_par];

      auto* DataCov = trackC.helix.Cov.GetMatrixArray();
      for(int id_cov = 0; id_cov<trackC.helix.Cov.GetNoElements();++id_cov)
	OutTrack->Seed_Cov[id_cov] = DataCov[id_cov];

      OutTrack->Chi2_C = trackC.helix.chi2_circle;
      OutTrack->Chi2_L = trackC.helix.chi2_line;

      for(size_t id_dethit = 0; id_dethit<trackC.sortedHits.size();++id_dethit)
	{
	  auto [id_det,id_hit,id_track] = TempHitToAllHits[id_dethit];
	  OutTrack->SetLayerID.emplace_back(id_det);
	  OutTrack->SetHitID.emplace_back(id_hit);
	  OutTrack->MCTrackID = id_track;
	}
      OutTrack->NSet = OutTrack->SetLayerID.size();
      OutTrack->TrackID = idT;
      ++idT;
    }

  for(auto [id_det,id_hit,id_track] : TempHitToAllHits)
    RecoEvent.IdHitsToMeasurement.emplace_back(id_det,id_hit,id_track);

  return res;
}

template<class Out>
ReturnRes::InfoM TTrackSeed<Out>::SoftExit(int ) { return ReturnRes::Fine; }

template<class Out>
void TTrackSeed<Out>::SelectHists()
{
  LocalHisto.h_SeedRiemannChi2    = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_SeedRiemannChi2);
  LocalHisto.h_SeedRiemannResidus = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_SeedRiemannResidus);
}

template<class Out>
void TTrackSeed<Out>::SortAndPosZ(std::vector<Seed::RTrack>& newTracksCand)
{

  for(auto& trackC : newTracksCand)
    {
      std::vector<std::vector<int> > tempSort(G4Sol::PSBE - G4Sol::MG01+1);
      att._logger->debug("Sort {}",tempSort.size());
      std::vector<int> hitWithZ;
      for(auto hitT : trackC.hits)
        {
          auto [id_det1, id_hit1, id_track1] = TempHitToAllHits[hitT];
          tempSort[id_det1 - G4Sol::MG01].emplace_back(hitT);
          if(TempHitXYZ[hitT].Z() > 0.1)
            hitWithZ.emplace_back(hitT);
        }

      //#ifdef DEBUG_RIEMANNFINDER
      //std::stringstream ssT;
      //int iddet = 0;
      for( auto hitIds : tempSort)
        {
	  //ssT<<" id "<<iddet++<<" [";
	  for(auto hitId : hitIds)
	    {
	      trackC.sortedHits.emplace_back(hitId);
	      //auto [id_det1, id_hit1, id_track1] = TempHitToAllHits[hitId];
	      //ssT<<hitId<<" ("<<id_det1<<" "<<id_hit1<<" "<<TempHitXYZ[hitId].Z()<<") ";
	    }
	  //ssT<<" ]";
	}
#ifdef DEBUG_RIEMANNFINDER
      att._logger->debug(ssT.str());
      att._logger->debug("HitsWithZ {}",hitWithZ.size());
#endif
      if(hitWithZ.size() >= 2)
        {
          TMatrixD Z(hitWithZ.size(), 2);
          TMatrixD R(hitWithZ.size(), 1);

          for(size_t i = 0; i < hitWithZ.size(); ++i)
            {
              Z(i, 0) = TempHitXYZ[hitWithZ[i]].Z();
              Z(i, 1) = 1.;
              R(i, 0) = TempHitXYZ[hitWithZ[i]].Perp();
            }

          TMatrixD ZtZ(Z, TMatrixD::kTransposeMult, Z);
          TMatrixD ZtR(Z, TMatrixD::kTransposeMult, R);
	  ZtZ.Invert();
	  //TMatrixD A = ZtZ * ZtR;
	  TMatrixD A(ZtZ, TMatrixD::kMult, ZtR); // R = A(0,0) * z + A(1,0);

          auto f_line = [&A](double R) { return (R - A(1, 0)) / TMath::Max(A(0, 0), 1.e-5); };

          for(auto hitT : trackC.hits)
            //if(TempHitXYZ[hitT].Z() <= 0.1)
              TempHitXYZ[hitT].SetZ(f_line(TempHitXYZ[hitT].Perp()));

          trackC.toRefit = true;
        }
    }
}

template<class Out>
void TTrackSeed<Out>::FitterRiemann(std::vector<Seed::RTrack>& newTracksCand)
{

  for(auto& trackC : newTracksCand)
    {
      size_t totalHit                    = trackC.sortedHits.size();
      tricktrack::Matrix3xNd riemannHits = tricktrack::Matrix3xNd::Random(3, totalHit);
      tricktrack::Matrix3Nd hits_cov     = tricktrack::Matrix3Nd::Identity(3 * totalHit, 3 * totalHit);

      int nhit = 0;
      for(auto it_hit : trackC.sortedHits)
        {
          double meanX1 = TempHitXYZ[it_hit].X();
          double meanY1 = TempHitXYZ[it_hit].Y();
          double meanZ1 = TempHitXYZ[it_hit].Z();

          riemannHits.col(nhit) << meanX1, meanY1, meanZ1;

          hits_cov(nhit, nhit)                               = TempCovXYZ[it_hit][0];
          hits_cov(nhit + totalHit, nhit + totalHit)         = TempCovXYZ[it_hit][4];
          hits_cov(nhit + 2 * totalHit, nhit + 2 * totalHit) = TempCovXYZ[it_hit][8];

	  hits_cov(nhit, nhit + totalHit)                    = TempCovXYZ[it_hit][1];
          hits_cov(totalHit + nhit, nhit)                    = TempCovXYZ[it_hit][3];

	  hits_cov(nhit, nhit + 2 * totalHit)                = TempCovXYZ[it_hit][2];
          hits_cov(2 * totalHit + nhit, nhit)                = TempCovXYZ[it_hit][6];

          hits_cov(nhit + totalHit, nhit + 2 * totalHit)     = TempCovXYZ[it_hit][5];
          hits_cov(2 * totalHit + nhit, nhit + totalHit)     = TempCovXYZ[it_hit][7];

	  //att._logger->debug("TempCov {} [ {} {} {}, {} {} {}, {} {} {}] ",it_hit, TempCovXYZ[it_hit][0],TempCovXYZ[it_hit][1],TempCovXYZ[it_hit][2],TempCovXYZ[it_hit][3],TempCovXYZ[it_hit][4],TempCovXYZ[it_hit][5],TempCovXYZ[it_hit][6],TempCovXYZ[it_hit][7],TempCovXYZ[it_hit][8]);

          ++nhit;
        }
#ifdef DEBUG_RIEMANNFINDER
      std::stringstream ssT;
      ssT <<"\n"<<riemannHits<<"\n";
      ssT<< hits_cov ;
      att._logger->debug("Fit :");
      att._logger->debug(ssT.str());
#endif

      //const tricktrack::Vector4d fast_fit = tricktrack::Fast_fit(riemannHits);

//       u_int nC                 = riemannHits.cols();
//       tricktrack::VectorNd rad = (riemannHits.block(0, 0, 2, nC).colwise().norm());
// #ifdef DEBUG_RIEMANNFINDER
//       std::stringstream radOut1;
//       radOut1 << fast_fit;
//       att._logger->debug("fast fit :");
//       att._logger->debug(radOut1.str());
// #endif
//       tricktrack::circle_fit circle = tricktrack::Circle_fit(
//           riemannHits.block(0, 0, 2, nC), hits_cov.block(0, 0, 2 * nC, 2 * nC), fast_fit, rad, true, true);

      res_h = tricktrack::Helix_fit(riemannHits, hits_cov, 1./(2.99792458e-3), true, true);

      for(int i = 0; i < 5; ++i)
        {
          trackC.helix.Par[i]    = res_h.par(i);
          trackC.helix.Cov[i][0] = res_h.cov(i, 0);
          trackC.helix.Cov[i][1] = res_h.cov(i, 1);
          trackC.helix.Cov[i][2] = res_h.cov(i, 2);
          trackC.helix.Cov[i][3] = res_h.cov(i, 3);
          trackC.helix.Cov[i][4] = res_h.cov(i, 4);
        }
      trackC.helix.q           = res_h.q;
      trackC.helix.chi2_circle = res_h.chi2_circle;
      trackC.helix.chi2_line   = res_h.chi2_line;
    }
}

template<class Out>
int TTrackSeed<Out>::TrackFit(FullRecoEvent& RecoEvent, std::vector<Seed::RTrack>& newTracksCand)
{

  for(auto it_trackDAF : RecoEvent.TrackDAF)
    {
      const int id_track  = it_trackDAF.first;
      auto it_ListHits    = it_trackDAF.second;

      const auto& it_fitRZ = RecoEvent.paramFitRZ.find(id_track);
      if(it_fitRZ == RecoEvent.paramFitRZ.end())
	{
	  att._logger->debug(" TrackSeed : track id {} did not pass the fitRZ -> no Riemmann Fit possible ",id_track);
	  continue;
	}

      std::set<int> hits;

      int n_MDC=0,n_PSCE=0, n_PSBE=0;
      for(int id_det = 0; id_det < it_ListHits.size(); ++id_det)
        {
          int id_hit = it_ListHits[id_det];

          if(id_hit < 0)
            continue;
	  if(id_det >= G4Sol::MG01 && id_det<= G4Sol::MG17)
	    ++n_MDC;
	  if( id_det == G4Sol::PSCE)
	    ++n_PSCE;
	  if(id_det == PSBE)
	    ++n_PSBE;
	}

      att._logger->debug("Track {}: MDC {} PSCE {} PSBE {}",id_track,n_MDC,n_PSCE, n_PSBE);
      for(int id_det = 0; id_det < it_ListHits.size(); ++id_det)
        {
          int id_hit = it_ListHits[id_det];

          if(id_hit < 0)
            continue;
	  if( !((id_det >= G4Sol::MG01 && id_det<= G4Sol::MG17) || id_det == G4Sol::PSCE || id_det == PSBE))
	    continue;

	  TempHitToAllHits.push_back(std::make_tuple(id_det, id_hit, id_track));
	  hits.insert(TempHitToAllHits.size()-1);

	  if(id_det >= G4Sol::MG01 && id_det<= G4Sol::MG17)
	    {
	      TMatrixD tempCov;
	      double temp_dl;
	      auto [status, meanX, meanY, meanZ] = extract_WireXY(RecoEvent.ListHits, id_det, id_hit, true, &tempCov, &temp_dl);

	      if(status!=true)
		att._logger->error("E> TrackSeed : extract Wire measurement data failed ! det:{} {}  hit:{} track:{}",id_det-G4Sol::MG01,G4Sol::nameLiteralDet.begin()[id_det],id_hit,id_track);

	      TempHitXYZ.emplace_back(meanX,meanY,meanZ);
	      std::array<double,9> tempC = {tempCov(0, 0), tempCov(1, 0), tempCov(2, 0), tempCov(0, 1), tempCov(1, 1), tempCov(2, 1), tempCov(0, 2), tempCov(1, 2), tempCov(2, 2)};
	      TempCovXYZ.emplace_back(tempC);
	    }
	  else
	    {
	      TMatrixD tempCov;
	      auto [meanX, meanY, meanZ] = extract_PSBXY(RecoEvent.ListHits, id_det, id_hit, true, &tempCov);

	      TempHitXYZ.emplace_back(meanX,meanY,meanZ);
	      std::array<double,9> tempC = {tempCov(0, 0), tempCov(1, 0), tempCov(2, 0), tempCov(0, 1), tempCov(1, 1), tempCov(2, 1), tempCov(0, 2), tempCov(1, 2), tempCov(2, 2)};
	      TempCovXYZ.emplace_back(tempC);
	    }
	}
      newTracksCand.emplace_back(hits, 0., 0., 0., -1., 0., 0., 0);
    }

  SortAndPosZ(newTracksCand);
  // for(int i=0; const auto& hit_tuple : TempHitToAllHits)
  //   {
  //     auto [id_det, id_hit, id_track] = hit_tuple;
  //     att._logger->debug("Sorter : TempHitToAllHits [{}] : det {} (from MDC1 {}), hit {}, track {}",i++, G4Sol::nameLiteralDet.begin()[id_det], id_det-G4Sol::MG01,id_hit,id_track);
  //   }

  // att._logger->debug("Check TempCov After Sorting: {}",TempCovXYZ.size());
  // for(int iC=0; auto trackC : newTracksCand)
  //   {
  //     att._logger->debug("track {}",iC++);
  //     for(auto it_hitT : trackC.sortedHits)
  //       {
  //         const auto& CovXYZ  = TempCovXYZ[it_hitT];

  // 	  att._logger->debug(": TempCov {} [ {} {} {}, {} {} {}, {} {} {}] ",it_hitT, CovXYZ[0],CovXYZ[1],CovXYZ[2],CovXYZ[3],CovXYZ[4],CovXYZ[5],CovXYZ[6],CovXYZ[7],CovXYZ[8]);

  // 	}

  //   }
  // att._logger->debug("---");

  FitterRiemann(newTracksCand);

  att._logger->debug("Fill histos ");

#ifdef DEBUG_RIEMANNFINDER
  int idTrack = 0;
#endif

  for(auto trackC : newTracksCand)
    {
#ifdef DEBUG_RIEMANNFINDER
      att._logger->debug("track {}", idTrack);
      std::stringstream ssT1;
      std::stringstream ssT2;
#endif

      LocalHisto.h_SeedRiemannChi2->Fill("chi2_pre",trackC.chi2,1.);
      if(trackC.toRefit)
	{
	  LocalHisto.h_SeedRiemannChi2->Fill("chi2_helix",trackC.helix.chi2_circle,1.);
	  LocalHisto.h_SeedRiemannChi2->Fill("chi2_line",trackC.helix.chi2_line,1.);
	}

#ifdef DEBUG_RIEMANNFINDER
      att._logger->debug(ssT1.str());
      att._logger->debug(ssT2.str());
      att._logger->debug("refit ? {} chi2 C:{} L:{}",trackC.toRefit, trackC.helix.chi2_circle, trackC.helix.chi2_line);
      ++idTrack;
#endif
    }
  return 0;
}


template<class Out>
std::tuple<int, double, double, double> TTrackSeed<Out>::extract_WireXY(const std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > >& ListHits, int id_det, int id_hit, bool DoCov, TMatrixD* cov, double* dl)
{
  att._logger->debug("extract WireXY: {} {} {}",fmt::ptr(dynamic_cast<genfit::WirePointMeasurement*>(ListHits[id_det][id_hit].get())), fmt::ptr(dynamic_cast<genfit::ProlateSpacepointMeasurement*>(ListHits[id_det][id_hit].get())), fmt::ptr(dynamic_cast<genfit::WireMeasurement*>(ListHits[id_det][id_hit].get())) );


  bool extractDone  = false;
  double meanX = 0., meanY = 0., meanZ = 0.;
  genfit::WirePointMeasurement* currentHit1 = dynamic_cast<genfit::WirePointMeasurement*>(ListHits[id_det][id_hit].get());

  if(currentHit1 != nullptr)
    {
      att._logger->debug("extract WirePointMeasurement");

      extractDone = true;
      auto HitCoord = currentHit1->getRawHitCoords();
      auto HitCov   = currentHit1->getRawHitCov();

      TVector3 wire_side1(HitCoord[0], HitCoord[1], HitCoord[2]);
      TVector3 wire_side2(HitCoord[3], HitCoord[4], HitCoord[5]);

      meanX = 0.5 * (wire_side1.X() + wire_side2.X());
      meanY = 0.5 * (wire_side1.Y() + wire_side2.Y());
      meanZ = HitCoord[7];
      if(DoCov)
	{
	  *dl = HitCoord[6];
	  cov->ResizeTo(3, 3);
	  cov->Zero();
	  (*cov)[0][0] = HitCov[6][6];
	  (*cov)[1][1] = HitCov[6][6];
	  (*cov)[2][2] = HitCov[7][7];
	}
    }

  genfit::ProlateSpacepointMeasurement* currentHit2 = dynamic_cast<genfit::ProlateSpacepointMeasurement*>(ListHits[id_det][id_hit].get());

  if(currentHit2 != nullptr)
    {
      att._logger->debug("extract ProlateSpacepointMeasurement");

      extractDone = true;
      auto HitCoord = currentHit2->getRawHitCoords();
      auto HitCov   = currentHit2->getRawHitCov();

      meanX = HitCoord[0];
      meanY = HitCoord[1];
      meanZ = HitCoord[2];
      if(DoCov)
	{
	  cov->ResizeTo(3, 3);
	  cov->Zero();
	  for(int i=0;i<3;++i)
	    for(int j=0;j<3;++j)
	      (*cov)[i][j] = HitCov[i][j];
	}
    }

  return std::make_tuple(extractDone, meanX, meanY,meanZ);
}

template<class Out>
std::tuple<double,double,double> TTrackSeed<Out>::extract_PSBXY(const  std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > >& ListHits, int id_det, int id_hit, bool DoCov, TMatrixD* cov)
{
  att._logger->debug("extract PSBXY");

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
  double meanZ       = currentLabPos.Z();
  if(DoCov == true)
    {
      TMatrixDSym hitCovXY(2);
      hitCovXY.Zero();
      hitCovXY(0, 0) = HitCov(0, 0);
      hitCovXY(1, 1) = 1. / 12.;
      TMatrixD rotMat(2, 2);
      rotMat.Zero();
      tempLocal.Set(1., 0.);
      currentLabPos = currentPlane->toLab(tempLocal);
      rotMat[0][0]  = currentLabPos.X();
      rotMat[0][1]  = currentLabPos.Y();
      tempLocal.Set(0., 1.);
      currentLabPos = currentPlane->toLab(tempLocal);
      rotMat[1][0]  = currentLabPos.X();
      rotMat[1][1]  = currentLabPos.Y();

      TMatrixD newCov(rotMat, TMatrixD::kMult, hitCovXY);
      rotMat.T();
      newCov *= rotMat;
      cov->ResizeTo(3, 3);
      cov->Zero();
      (*cov)[0][0] = newCov[0][0];
      (*cov)[1][0] = newCov[1][0];
      (*cov)[0][1] = newCov[0][1];
      (*cov)[1][1] = newCov[1][1];
      (*cov)[2][2] = HitCov(1, 1);
    }
  return std::make_tuple(meanX, meanY, meanZ);
}




template class TTrackSeed<MCAnaEventG4Sol>;
template class TTrackSeed<Ana_WasaEvent>;
