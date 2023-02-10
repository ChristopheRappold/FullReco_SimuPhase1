#include "TCheckRZ.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"
#include "FullRecoEvent.hh"
#include "KalmanFittedStateOnPlane.h"
#include "KalmanFitterInfo.h"
#include "PlanarMeasurement.h"
#include "ReturnRes.hh"
#include "StateOnPlane.h"
#include "TLinearFitter.h"

#include <tuple>

//#define DEBUG_CHECKRZ

using namespace std;
using namespace G4Sol;

template<class Out>
TCheckRZ<Out>::TCheckRZ(const THyphiAttributes& attribut) : TDataProcessInterface<Out>("Check_R_Z"), att(attribut)
{
  ChangeMiniFiber = att.RZ_ChangeMiniFiber;

  if(att.RZ_MDCProlate)
    MDCWireType = 0;
  if(att.RZ_MDCWire2)
    MDCWireType = 1;

  MDCBiasCorr = att.RZ_MDCBiasCorr;
  if(MDCBiasCorr == false)
    {
      for(auto id : correctBiasPSCE)
        id = 0.;
      for(auto id : correctBiasPSEndCap)
        id = 0.;
    }
  // else
  //   {
  //     for(size_t i=0;i<correctBiasPSCE.size();++i)
  //       correctBiasPSEndCap[i] += correctBiasPSCE[i];
  //   }
}

template<class Out>
TCheckRZ<Out>::~TCheckRZ() {}

template<class Out>
void TCheckRZ<Out>::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }

template<class Out>
ReturnRes::InfoM TCheckRZ<Out>::operator()(FullRecoEvent& RecoEvent, Out* OutTree)
{

  int result_finder = Exec(RecoEvent, OutTree);

  return SoftExit(result_finder);
}

template<class Out>
int TCheckRZ<Out>::Exec(FullRecoEvent& RecoEvent, Out* ) { return FinderTrack(RecoEvent); }

template<class Out>
ReturnRes::InfoM TCheckRZ<Out>::SoftExit(int ) { return ReturnRes::Fine; }

template<class Out>
void TCheckRZ<Out>::SelectHists()
{

  LocalHisto.h_RZStats      = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_RZStats);
  LocalHisto.h_RZ           = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_RZ);
  LocalHisto.h_XYfit_miniF  = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_XYfit_miniF);
  LocalHisto.h_RZfit_mom    = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_RZfit_mom);
  LocalHisto.h_RZfit_Chi2   = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_RZfit_Chi2);
  LocalHisto.h_MDC_Z_residu = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_MDC_Z_residu);
  LocalHisto.h_MDC_R_residu = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_MDC_R_residu);
  LocalHisto.h_MDC_Z_pull   = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_MDC_Z_pull);
  LocalHisto.h_MDC_R_pull   = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_MDC_R_pull);
}

template<class Out>
int TCheckRZ<Out>::FinderTrack(FullRecoEvent& RecoEvent)
{

  // auto printW = [](const auto a, const int width) -> std::string {
  //   std::stringstream ss;
  //   ss << std::fixed << std::right;
  //   ss.fill(' ');    // fill space around displayed #
  //   ss.width(width); // set  width around displayed #
  //   ss << a;
  //   return ss.str();
  // };
  // auto printFixed = [](const double a, const int decDigits, const int width) -> std::string {
  //   std::stringstream ss;
  //   ss << std::fixed << std::right;
  //   ss.fill(' ');            // fill space around displayed #
  //   ss.width(width);         // set  width around displayed #
  //   ss.precision(decDigits); // set # places after decimal
  //   ss << a;
  //   return ss.str();
  // };

  // unsigned int Nb_trackSeed = RecoEvent.TrackDAF.size();

  for(size_t id_det = G4Sol::MG01; id_det <= G4Sol::MG17; ++id_det)
    RecoEvent.OldListHits[id_det].resize(RecoEvent.ListHits[id_det].size());

  for(size_t id_det = G4Sol::MiniFiberD1_x1; id_det <= G4Sol::MiniFiberD1_v2; ++id_det)
    RecoEvent.OldListHits[id_det].resize(RecoEvent.ListHits[id_det].size());

  RecoEvent.OldListHits[G4Sol::PSBE].resize(RecoEvent.ListHits[G4Sol::PSBE].size());

  int ntrack = -1;
  for(auto it_trackInfo : RecoEvent.TrackInfo)
    {
      ntrack++;

#ifdef DEBUG_CHECKRZ
      att._logger->debug("start check RZ #{}", ntrack); //<<std::endl;
      // std::cout<<" | Pointsize "<<Vtracks->getNumPoints()<<" Measurements
      // "<<Vtracks->getNumPointsWithMeasurement()<<" Rep
      // "<<Vtracks->getNumReps()<<std::endl;
      // std::vector<std::string> name_det;
#endif
      LocalHisto.h_RZStats->Fill("allRZ", 1.);
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

      std::set<std::tuple<int, int> > id_dets_MDC;
      std::set<std::tuple<int, int> > id_dets_MiniFiber;
      std::set<std::tuple<int, int> > id_dets_PSCE;

      G4Sol::SolDet LastFrontWall = att.Wasa_Side == 0 ? G4Sol::PSFE : G4Sol::PSBE;

      std::set<std::tuple<int, int> > id_dets_PSEndcap;

      int n_Central   = 0;
      int n_MiniFiber = 0;
      // std::cout << "it_ListHits->second.size : " << it_ListHits->second.size() << std::endl;
      for(size_t id_det = 0; id_det < it_ListHits->second.size(); ++id_det)
        {
          int id_hit = it_ListHits->second[id_det];
          if(id_hit < 0)
            continue;
          if(id_det >= G4Sol::MG01 && id_det <= G4Sol::MG17)
            {
              ++n_Central;
              id_dets_MDC.insert(std::make_tuple(id_det, id_hit));
            }
          if(id_det >= G4Sol::MiniFiberD1_x1 && id_det <= G4Sol::MiniFiberD1_v2)
            {
              ++n_MiniFiber;
              id_dets_MiniFiber.insert(std::make_tuple(id_det, id_hit));
            }
          if(id_det == G4Sol::PSCE)
            id_dets_PSCE.insert(std::make_tuple(id_det, id_hit));
          if(id_det == LastFrontWall)
            id_dets_PSEndcap.insert(std::make_tuple(id_det, id_hit));

#ifdef DEBUG_CHECKRZ
          att._logger->debug("id_det / id_hit: {} / {}", id_det, id_hit);
#endif
          // id_dets.insert(std::make_tuple(getZpos(currentHit), id_det, id_hit));
          id_dets.insert(std::make_tuple(id_det, id_det, id_hit));
          // std::cout << " id_det : " << id_det << std::endl;
        }

      if(id_dets.size() < 3) //(nb_ValidHits < 3)
        {
#ifdef DEBUG_CHECKRZ
          att._logger->debug("!> less than 3 measurements: {}", id_dets.size());
#endif
          LocalHisto.h_RZStats->Fill("Rej<3Mes", 1.);
          continue;
        }

      auto f_LastHitIsValid = [](const auto& it_ListHits, const std::set<G4Sol::SolDet>& listToTest) {
        for(auto it_det : listToTest)
          if(it_ListHits->second[it_det] >= 0)
            return it_det;

        return G4Sol::SIZEOF_G4SOLDETTYPE;
      };

      // G4Sol::SolDet LastFrontWall = att.Wasa_Side == 0 ? G4Sol::PSFE : G4Sol::PSBE;
      G4Sol::SolDet lastValidHit = f_LastHitIsValid(it_ListHits, {G4Sol::PSCE, LastFrontWall});

      if(lastValidHit == G4Sol::SIZEOF_G4SOLDETTYPE)
        {
          LocalHisto.h_RZStats->Fill("NotPSCE", 1.);
          continue;
        }
      if(n_MiniFiber < 4)
        {
          LocalHisto.h_RZStats->Fill("CutFiber<4", 1.);
          continue;
        }

      LocalHisto.h_RZStats->Fill("Good", 1.);

      // auto firstHit   = id_dets.cbegin();
      // int id_firstDet = std::get<1>(*firstHit);

      // const InfoPar track_state(it_trackInfo.second[id_firstDet]);

      auto it_ListHitsSim = RecoEvent.TrackDAFSim.find(id_track);

      // firstHit             = id_dets.cbegin();
      // id_firstDet          = std::get<1>(*firstHit);
      // const auto lastHit   = id_dets.crbegin();
      // const int id_lastDet = std::get<1>(*lastHit); // std::get<0>(lastHit);

      std::vector<std::array<double, 3> > Hits_MiniFiber_RZ;
      std::vector<std::array<double, 3> > Hits_PSLast_RZ;

      auto f_extractHitLabs = [](const auto& id_dets, auto& Hits_RZ, const auto& ListHits) {
        for(auto ids : id_dets)
          {
            int id_det = std::get<0>(ids);
            int id_hit = std::get<1>(ids);

            genfit::PlanarMeasurement* currentHit =
                dynamic_cast<genfit::PlanarMeasurement*>(ListHits[id_det][id_hit].get());
            auto dummyState   = genfit::StateOnPlane();
            auto currentPlane = currentHit->constructPlane(genfit::StateOnPlane());
            auto HitCoord     = currentHit->getRawHitCoords();
            auto HitCov       = currentHit->getRawHitCov();
            // std::cout<<"f_extractHitLabs :\n";
            // HitCoord.Print();
            // HitCov.Print();

            TVector2 tempLocal;
            if(HitCoord.GetNrows() == 1)
              tempLocal.Set(HitCoord[0], 0.);
            else
              tempLocal.Set(HitCoord[0], HitCoord[1]);

            auto currentLabPos = currentPlane->toLab(tempLocal);

            TMatrixD currentLabCov(3, 3);
            currentLabCov.Zero();

            double scale    = TMath::Sqrt(HitCov.GetNrows());
            auto tempLabPos = decltype(currentLabPos)(0., 0., 0.);

            std::vector<decltype(tempLabPos)> listLabPos;

            size_t nSize = HitCov.GetNrows();
            std::vector<std::array<double, 2> > Errs(nSize * 2, {0., 0.});
            Errs[0] = {scale * TMath::Sqrt(HitCov[0][0]), 0.};
            Errs[1] = {-scale * TMath::Sqrt(HitCov[0][0]), 0.};
            if(nSize == 2)
              {
                Errs[2] = {0., scale * TMath::Sqrt(HitCov[1][1])};
                Errs[3] = {0., -scale * TMath::Sqrt(HitCov[1][1])};
              }
            for(auto Err : Errs)
              {
                TVector2 tempLocal_err(HitCoord[0] + Err[0], 0.);
                if(nSize == 2)
                  tempLocal_err.SetY(HitCoord[1] + Err[1]);

                auto currentLabPos_err = currentPlane->toLab(tempLocal_err);

                tempLabPos += currentLabPos_err;
                listLabPos.emplace_back(currentLabPos_err);
              }

            tempLabPos *= 1. / (double)listLabPos.size();
            // tempLabPos.Print();
            for(auto tempVecPos : listLabPos)
              {
                // std::cout<<"f_extract : "<<tempVecPos[0]<<" "<<tempVecPos[1]<<" "<<tempVecPos[2]<<"\n";
                TMatrixD temp1(3, 1);
                for(int iM = 0; iM < temp1.GetNrows(); ++iM)
                  temp1(iM, 0) = (TMath::Floor(1000 * tempVecPos[iM]) * 0.001) - tempLabPos[iM];
                auto temp2 = temp1;
                temp2.T();
                currentLabCov += temp1 * temp2;
              }
            currentLabCov *= 1. / (double)listLabPos.size();

            auto SubCov = currentLabCov.GetSub(0, 1, 0, 1);

            auto f_sqrtMat2x2 = [](TMatrixD& m) -> TMatrixD {
              double det_m   = m.Determinant();
              double trace_m = m(0, 0) + m(1, 1);

              double Norm = TMath::Sqrt(trace_m + 2 * TMath::Sqrt(det_m));

              TMatrixD U(2, 2);
              U.UnitMatrix();

              auto mSqrt = (m + TMath::Sqrt(Norm) * U);
              mSqrt *= 1. / Norm;
              return mSqrt;
            };

            TMatrixD SubCovSqrt = f_sqrtMat2x2(SubCov);
            SubCovSqrt *= TMath::Sqrt(2 * SubCov.GetNcols());
            std::vector<double> Rs;
            double mean_R = 0.;
            for(int iC = 0; iC < SubCovSqrt.GetNcols(); ++iC)
              {
                double x1 = tempLabPos.X() + SubCovSqrt(0, iC);
                double y1 = tempLabPos.Y() + SubCovSqrt(0, iC);
                double r1 = TMath::Sqrt(x1 * x1 + y1 * y1);
                Rs.emplace_back(r1);
                mean_R += r1;
                double x2 = tempLabPos.X() - SubCovSqrt(0, iC);
                double y2 = tempLabPos.Y() - SubCovSqrt(0, iC);
                double r2 = TMath::Sqrt(x2 * x2 + y2 * y2);
                Rs.emplace_back(r2);
                mean_R += r2;
              }
            mean_R /= (double)Rs.size();
            double cov_R = 0.;
            for(auto R : Rs)
              cov_R += TMath::Sq(R - mean_R);
            cov_R /= (double)Rs.size();

            // std::cout<<"f_extract : Cov\n";
            // currentLabCov.Print();
            std::array<double, 3> temp_hitRZ = {mean_R, tempLabPos.Z(), TMath::Sqrt(cov_R)};
            // TMath::Sqrt(currentLabCov(0, 0) + currentLabCov(1, 1) + 2 * currentLabCov(0, 1))};
            // tempLabPos.Perp()*TMath::Sqrt(currentLabCov(0,0)*tempLabPos.X()*tempLabPos.X()+currentLabCov(1,1)*tempLabPos.Y()*tempLabPos.Y()+2.*currentLabCov(0,1)*tempLabPos.X()*tempLabPos.Y())};

            Hits_RZ.emplace_back(temp_hitRZ);
          }
      };

      auto f_extractPSEndcapHit = [](const auto& id_dets, auto& Hits_RZ, const auto& ListHits, const auto& TrackFiber,
                                     const auto& TrackFiberCov) {
        for(auto ids : id_dets)
          {
            int id_det = std::get<0>(ids);
            int id_hit = std::get<1>(ids);

            genfit::PlanarMeasurement* currentHit =
                dynamic_cast<genfit::PlanarMeasurement*>(ListHits[id_det][id_hit].get());
            auto dummyState   = genfit::StateOnPlane();
            auto currentPlane = currentHit->constructPlane(genfit::StateOnPlane());
            auto HitCoord     = currentHit->getRawHitCoords();
            auto HitCov       = currentHit->getRawHitCov();
            // std::cout<<"f_extractHitLabs :\n";

            TVector3 Orig = currentPlane->getO();

            auto f_sqrtMat = [](const TMatrixD& m) -> TMatrixD {
              TVectorD EigVal;
              TMatrixD EigVec = m.EigenVectors(EigVal);
              TMatrixD m_new(m.GetNrows(), m.GetNcols());
              m_new.UnitMatrix();
              for(int i = 0; i < m.GetNrows(); ++i)
                m_new(i, i) = EigVal(i);
              m_new.Sqrt();
              auto EigVecInv = EigVec;
              EigVecInv.Invert();
              auto mSqrt = EigVec * m_new * EigVecInv;
              return mSqrt;
            };

            TMatrixD SubCovSqrt = f_sqrtMat(TrackFiberCov);
            SubCovSqrt *= TMath::Sqrt(2 * TrackFiberCov.GetNcols());
            std::vector<double> Rs;
            double mean_R = 0.;
            for(int iC = 0; iC < SubCovSqrt.GetNcols(); ++iC)
              {
                double XPosAtZ = (TrackFiber(2) + SubCovSqrt(2, iC)) * Orig.Z() + (TrackFiber(0) + SubCovSqrt(0, iC));
                double YPosAtZ = (TrackFiber(3) + SubCovSqrt(3, iC)) * Orig.Z() + (TrackFiber(1) + SubCovSqrt(1, iC));

                double RPosAtZ = TMath::Sqrt(XPosAtZ * XPosAtZ + YPosAtZ * YPosAtZ);

                Rs.emplace_back(RPosAtZ);
                mean_R += RPosAtZ;

                double XPosAtZ2 = (TrackFiber(2) - SubCovSqrt(2, iC)) * Orig.Z() + (TrackFiber(0) - SubCovSqrt(0, iC));
                double YPosAtZ2 = (TrackFiber(3) - SubCovSqrt(3, iC)) * Orig.Z() + (TrackFiber(1) - SubCovSqrt(1, iC));

                double RPosAtZ2 = TMath::Sqrt(XPosAtZ2 * XPosAtZ2 + YPosAtZ2 * YPosAtZ2);
                Rs.emplace_back(RPosAtZ2);
                mean_R += RPosAtZ2;
              }

            mean_R /= (double)Rs.size();
            double cov_R = 0.;
            for(auto R : Rs)
              cov_R += TMath::Sq(R - mean_R);
            cov_R /= (double)Rs.size();

            std::array<double, 3> temp_hitRZ = {mean_R, Orig.Z(), TMath::Sqrt(cov_R)};

            Hits_RZ.emplace_back(temp_hitRZ);

            // double CovRAtZ = TMath::Sqrt(miniTrackCov(0,0)*XPosAtZ*XPosAtZ/RPosAtZ/RPosAtZ +
            // miniTrackCov(1,1)*YPosAtZ*YPosAtZ/RPosAtZ/RPosAtZ+2.*miniTrackCov(1,0)*XPosAtZ*YPosAtZ/RPosAtZ/RPosAtZ);
            // double CovRAtZ = TMath::Sqrt(miniTrackCov(0, 0) + miniTrackCov(1, 1) + 2. * miniTrackCov(1, 0));
          }
      };

      auto f_extractTrackMF = [](const auto& id_dets, const auto& ListHits, auto& Res, auto& Cov,
                                 auto& ResOnly1) -> std::tuple<int, int> {
        auto f_extractTrackMiniF = [](const auto& id_dets, const auto& ListHits, const std::vector<double>& div,
                                      auto& Res) -> int {
          TMatrixD A(4, 4); //  least square def : A.(x0 y0 tx ty)^T = U // { x = alpha_x*z+x0 , y = alpha_y*z+y0 }
          TVectorD U(4);
          A.Zero();
          U.Zero();

          int id_div = 0;
          for(auto ids : id_dets)
            {

              int id_det = std::get<0>(ids);
              int id_hit = std::get<1>(ids);

              genfit::PlanarMeasurement* currentHit =
                  dynamic_cast<genfit::PlanarMeasurement*>(ListHits[id_det][id_hit].get());
              auto dummyState   = genfit::StateOnPlane();
              auto currentPlane = currentHit->constructPlane(genfit::StateOnPlane());
              auto HitCoord     = currentHit->getRawHitCoords();
              auto HitCov       = currentHit->getRawHitCov();
              auto currentO     = currentPlane->getO();
              auto currentU     = currentPlane->getU();
              auto currentV     = currentPlane->getV();

              double sinTheta  = -currentV.X();
              double cosTheta  = currentU.X();
              double sinTheta2 = sinTheta * sinTheta;
              double cosTheta2 = cosTheta * cosTheta;
              double z         = currentO.Z();
              double z2        = z * z;
              double Cov2      = HitCov(0, 0);
              double u         = HitCoord(0) + div[id_div];

              A(0, 0) += cosTheta2 * z / Cov2;
              A(0, 1) += cosTheta * sinTheta * z / Cov2;
              A(0, 2) += cosTheta2 * z2 / Cov2;
              A(0, 3) += cosTheta * sinTheta * z2 / Cov2;

              A(1, 0) += cosTheta2 / Cov2;
              A(1, 1) += cosTheta * sinTheta / Cov2;
              A(1, 2) += cosTheta2 * z / Cov2;
              A(1, 3) += cosTheta * sinTheta * z / Cov2;

              A(2, 0) += cosTheta * sinTheta * z / Cov2;
              A(2, 1) += sinTheta2 * z / Cov2;
              A(2, 2) += cosTheta * sinTheta * z2 / Cov2;
              A(2, 3) += sinTheta2 * z2 / Cov2;

              A(3, 0) += cosTheta * sinTheta / Cov2;
              A(3, 1) += sinTheta2 / Cov2;
              A(3, 2) += cosTheta * sinTheta * z / Cov2;
              A(3, 3) += sinTheta2 * z / Cov2;

              U(0) += cosTheta * z * u / Cov2;
              U(1) += cosTheta * u / Cov2;
              U(2) += sinTheta * z * u / Cov2;
              U(3) += sinTheta * u / Cov2;

              ++id_div;
            }
          double det_A  = 0;
          double maxInA = A.Max();
          A *= 1. / maxInA;
          A.Invert(&det_A);

          if(TMath::Abs(det_A) < 1e-30)
            return 0;

          A *= 1. / maxInA;

          Res = A * U;
          return 1;
        };

        std::vector<double> div_null(id_dets.size(), 0.);
        int status = f_extractTrackMiniF(id_dets, ListHits, div_null, ResOnly1);

        std::vector<TVectorD> miniTrackAll(id_dets.size() * 2, TVectorD(4));
        int id_div       = 0;
        size_t statusAll = 0;
        for(auto ids : id_dets)
          {

            int id_det = std::get<0>(ids);
            int id_hit = std::get<1>(ids);

            genfit::PlanarMeasurement* currentHit =
                dynamic_cast<genfit::PlanarMeasurement*>(ListHits[id_det][id_hit].get());
            auto dummyState   = genfit::StateOnPlane();
            auto currentPlane = currentHit->constructPlane(genfit::StateOnPlane());
            auto HitCoord     = currentHit->getRawHitCoords();
            auto HitCov       = currentHit->getRawHitCov();

            double Cov2  = HitCov(0, 0);
            double scale = TMath::Sqrt(id_dets.size());
            std::vector<double> div_pos(id_dets.size(), 0.), div_neg(id_dets.size(), 0.);
            div_pos[id_div] = scale * TMath::Sqrt(Cov2);
            div_neg[id_div] = -scale * TMath::Sqrt(Cov2);

            miniTrackAll[id_div].Zero();
            int statusP = f_extractTrackMiniF(id_dets, ListHits, div_pos, miniTrackAll[id_div]);

            miniTrackAll[id_div + id_dets.size()].Zero();
            int statusN = f_extractTrackMiniF(id_dets, ListHits, div_neg, miniTrackAll[id_div + id_dets.size()]);
            ++id_div;
            if(statusN != 1 || statusP != 1)
              break;
            statusAll += statusN + statusP;
          }
        if(statusAll != miniTrackAll.size())
          return std::make_tuple(status, 0);

        Res.Zero();
        for(auto tempV : miniTrackAll)
          Res += tempV;

        Res *= 1. / (double)(miniTrackAll.size());

        Cov.Zero();
        for(auto tempV : miniTrackAll)
          {
            TMatrixD temp1(4, 1);
            for(int iM = 0; iM < temp1.GetNrows(); ++iM)
              temp1(iM, 0) = tempV[iM] - Res[iM];
            auto temp2 = temp1;
            temp2.T();
            Cov += temp1 * temp2;
          }

        Cov *= 1. / (double)miniTrackAll.size();
        return std::make_tuple(status, 1);
      };

      TVectorD miniTrack(4);
      TVectorD miniTrackOnly1(4);
      TMatrixD miniTrackCov(4, 4);
      miniTrack.Zero();
      auto [statusOnly1, statusWithCov] =
          f_extractTrackMF(id_dets_MiniFiber, RecoEvent.ListHits, miniTrack, miniTrackCov, miniTrackOnly1);

      if(statusWithCov == 0)
        {
          if(statusOnly1 == 0)
            {
              LocalHisto.h_RZStats->Fill("MiniF_FitFailedBoth", 1.);
              continue;
            }
          LocalHisto.h_RZStats->Fill("MiniF_FitnoCov", 1.);

          miniTrack = miniTrackOnly1;
          miniTrackCov.Zero();
          miniTrackCov(0, 0) = 0.1 * 0.1;
        }

      TVector3 momRefAtminiF, PosRefAtminiF;
      {
        auto ids = id_dets_MiniFiber.begin();

        int id_det  = std::get<0>(*ids);
        auto it_hit = it_ListHitsSim->second[id_det][0];
        PosRefAtminiF.SetXYZ(it_hit.hitX, it_hit.hitY, it_hit.hitZ);
        momRefAtminiF.SetXYZ(it_hit.momX, it_hit.momY, it_hit.momZ);
      }
      TVectorD miniTrackRef(4);
      miniTrackRef(2) = momRefAtminiF.X() / momRefAtminiF.Z();
      miniTrackRef(3) = momRefAtminiF.Y() / momRefAtminiF.Z();
      miniTrackRef(0) = PosRefAtminiF.X() - miniTrackRef(2) * PosRefAtminiF.Z();
      miniTrackRef(1) = PosRefAtminiF.Y() - miniTrackRef(3) * PosRefAtminiF.Z();

      static const std::vector<std::string> nameXYfit = {"alpha_x",      "alpha_y",      "x0",      "y0",
                                                         "alpha_x_pull", "alpha_y_pull", "x0_pull", "y0_pull"};

      LocalHisto.h_XYfit_miniF->Fill(nameXYfit[0].c_str(), miniTrack(2) - miniTrackRef(2), 1.);
      LocalHisto.h_XYfit_miniF->Fill(nameXYfit[1].c_str(), miniTrack(3) - miniTrackRef(3), 1.);
      // LocalHisto.h_XYfit_miniF->Fill(nameXYfit[2].c_str(), miniTrack(0) - miniTrackRef(0), 1.);
      // LocalHisto.h_XYfit_miniF->Fill(nameXYfit[3].c_str(), miniTrack(1) - miniTrackRef(1), 1.);
      LocalHisto.h_XYfit_miniF->Fill(nameXYfit[2].c_str(),
                                     (miniTrack(2) * PosRefAtminiF.Z() + miniTrack(0)) -
                                         (miniTrackRef(2) * PosRefAtminiF.Z() + miniTrackRef(0)),
                                     1.);
      LocalHisto.h_XYfit_miniF->Fill(nameXYfit[3].c_str(),
                                     (miniTrack(3) * PosRefAtminiF.Z() + miniTrack(1)) -
                                         (miniTrackRef(3) * PosRefAtminiF.Z() + miniTrackRef(1)),
                                     1.);

      LocalHisto.h_XYfit_miniF->Fill(nameXYfit[4].c_str(),
                                     (miniTrack(2) - miniTrackRef(2)) / TMath::Sqrt(miniTrackCov(2, 2)), 1.);
      LocalHisto.h_XYfit_miniF->Fill(nameXYfit[5].c_str(),
                                     (miniTrack(3) - miniTrackRef(3)) / TMath::Sqrt(miniTrackCov(3, 3)), 1.);
      LocalHisto.h_XYfit_miniF->Fill(nameXYfit[6].c_str(),
                                     (miniTrack(0) - miniTrackRef(0)) / TMath::Sqrt(miniTrackCov(0, 0)), 1.);
      LocalHisto.h_XYfit_miniF->Fill(nameXYfit[7].c_str(),
                                     (miniTrack(1) - miniTrackRef(1)) / TMath::Sqrt(miniTrackCov(1, 1)), 1.);

      LocalHisto.h_XYfit_miniF->Fill("Dalpha_x", miniTrack(2) - miniTrackOnly1(2), 1.);
      LocalHisto.h_XYfit_miniF->Fill("Dalpha_y", miniTrack(3) - miniTrackOnly1(3), 1.);
      LocalHisto.h_XYfit_miniF->Fill("Dx0", miniTrack(0) - miniTrackOnly1(0), 1.);
      LocalHisto.h_XYfit_miniF->Fill("Dy0", miniTrack(1) - miniTrackOnly1(1), 1.);

#ifdef DEBUG_CHECKRZ
      att._logger->debug("MiniTrack :");
      att._logger->debug("zx-plane : x = alpha_x*z+x0 : alpha_x={} x0={}", miniTrack(2), miniTrack(0));
      att._logger->debug("zy-plane : y = alpha_y*z+y0 : alpha_y={} y0={}", miniTrack(3), miniTrack(1));

      att._logger->debug("Ref zx-plane : x = alpha_x*z+x0 : alpha_x={} x0={}", miniTrackRef(2), miniTrackRef(0));
      att._logger->debug("Ref zy-plane : y = alpha_y*z+y0 : alpha_y={} y0={}", miniTrackRef(3), miniTrackRef(1));

#endif
      int EndCap_Case = 0;
      if(id_dets_PSCE.size() != 0)
        f_extractHitLabs(id_dets_PSCE, Hits_PSLast_RZ, RecoEvent.ListHits);
      else if(id_dets_PSEndcap.size() != 0)
        {
          f_extractPSEndcapHit(id_dets_PSEndcap, Hits_PSLast_RZ, RecoEvent.ListHits, miniTrack, miniTrackCov);
          EndCap_Case = 1;
        }
      else
        att._logger->error("E> CheckRZ does not have PSCE and LastPSE even after selection !");

#ifdef DEBUG_CHECKRZ
      att._logger->debug("Extracted Hit in Lab ref");

      att._logger->debug("MiniFiber {}", Hits_MiniFiber_RZ.size());
      for(auto tempHit : Hits_MiniFiber_RZ)
        att._logger->debug(" [Z, R, errR] : {}, {}, {}", tempHit[1], tempHit[0], tempHit[2]);

      for(auto ids : id_dets_MiniFiber)
        {
          int id_det  = std::get<0>(ids);
          auto it_hit = it_ListHitsSim->second[id_det];
          att._logger->debug("Hit in {} : [x, y, z] : {} {} {} | [Z, R] {}, {}", G4Sol::nameLiteralDet.begin()[id_det],
                             it_hit.hitX, it_hit.hitY, it_hit.hitZ,
                             TMath::Sqrt(it_hit.hitX * it_hit.hitX + it_hit.hitY * it_hit.hitY), it_hit.hitZ);
        }

      att._logger->debug("PSCE {}", Hits_PSLast_RZ.size());
      for(auto tempHit : Hits_PSLast_RZ)
        att._logger->debug(" [Z, R, errR] : {}, {}, {}", tempHit[1], tempHit[0], tempHit[2]);

      for(auto ids : id_dets_PSCE)
        {
          int id_det  = std::get<0>(ids);
          auto it_hit = it_ListHitsSim->second[id_det];
          att._logger->debug("Hit in {} : [x, y, z] : {} {} {} | [Z, R] {}, {}", G4Sol::nameLiteralDet.begin()[id_det],
                             it_hit.hitX, it_hit.hitY, it_hit.hitZ,
                             TMath::Sqrt(it_hit.hitX * it_hit.hitX + it_hit.hitY * it_hit.hitY), it_hit.hitZ);
        }
#endif

      TVector3 momAtFiber;
      for(auto ids : id_dets_MiniFiber)
        {
          int id_det  = std::get<0>(ids);
          auto it_hit = it_ListHitsSim->second[id_det][0];
          // Hits_MiniFiber_RZ.emplace_back(std::array<double,3>{TMath::Sqrt(it_hit.hitX*it_hit.hitX+it_hit.hitY*it_hit.hitY),
          // it_hit.hitZ, 0.1});

          double XPosAtZ = miniTrack(2) * it_hit.hitZ + miniTrack(0);
          double YPosAtZ = miniTrack(3) * it_hit.hitZ + miniTrack(1);

          double RPosAtZ = TMath::Sqrt(XPosAtZ * XPosAtZ + YPosAtZ * YPosAtZ);
          // double CovRAtZ = TMath::Sqrt(miniTrackCov(0,0)*XPosAtZ*XPosAtZ/RPosAtZ/RPosAtZ +
          // miniTrackCov(1,1)*YPosAtZ*YPosAtZ/RPosAtZ/RPosAtZ+2.*miniTrackCov(1,0)*XPosAtZ*YPosAtZ/RPosAtZ/RPosAtZ);
          double CovRAtZ = TMath::Sqrt(miniTrackCov(0, 0) + miniTrackCov(1, 1) + 2. * miniTrackCov(1, 0));

          Hits_MiniFiber_RZ.emplace_back(std::array<double, 3>{RPosAtZ, it_hit.hitZ, CovRAtZ});

          momAtFiber.SetXYZ(it_hit.momX, it_hit.momY, it_hit.momZ);
        }

      // for(auto ids : id_dets_PSCE)
      //   {
      //     int id_det = std::get<0>(ids);
      //     auto it_hit = it_ListHitsSim->second[id_det];
      //     Hits_PSLast_RZ.emplace_back(std::array<double,3>{TMath::Sqrt(it_hit.hitX*it_hit.hitX+it_hit.hitY*it_hit.hitY),
      //     it_hit.hitZ,0.1});
      //   }

      TLinearFitter fitRZ(2, "pol1");
      fitRZ.StoreData(false);
      std::unique_ptr<Double_t[]> Zall(new Double_t[Hits_MiniFiber_RZ.size() + Hits_PSLast_RZ.size()]);
      std::unique_ptr<Double_t[]> Rall(new Double_t[Hits_MiniFiber_RZ.size() + Hits_PSLast_RZ.size()]);
      std::unique_ptr<Double_t[]> Eall(new Double_t[Hits_MiniFiber_RZ.size() + Hits_PSLast_RZ.size()]);

      for(size_t iD = 0; iD < Hits_MiniFiber_RZ.size(); ++iD)
        {
          Rall.get()[iD] = Hits_MiniFiber_RZ[iD][0];
          Zall.get()[iD] = Hits_MiniFiber_RZ[iD][1];
          LocalHisto.h_RZ->Fill(Hits_MiniFiber_RZ[iD][1] / 100., Hits_MiniFiber_RZ[iD][0] / 100.);
          Eall.get()[iD] = Hits_MiniFiber_RZ[iD][2];
        }
      for(size_t iD = 0; iD < Hits_PSLast_RZ.size(); ++iD)
        {
          Rall.get()[iD + Hits_MiniFiber_RZ.size()] = Hits_PSLast_RZ[iD][0];
          Zall.get()[iD + Hits_MiniFiber_RZ.size()] = Hits_PSLast_RZ[iD][1];
          LocalHisto.h_RZ->Fill(Hits_PSLast_RZ[iD][1] / 100., Hits_PSLast_RZ[iD][0] / 100.);
          Eall.get()[iD + Hits_MiniFiber_RZ.size()] = Hits_PSLast_RZ[iD][2];
          if(EndCap_Case == 1)
            {
              if(5.8 <= Hits_PSLast_RZ[iD][0] && Hits_PSLast_RZ[iD][0] <= 22.2)
                LocalHisto.h_RZStats->Fill("PSEndCap_InsideBar", 1);
              else
                LocalHisto.h_RZStats->Fill("PSEndCap_OutsideBar", 1);
            }
        }

      fitRZ.AssignData(Hits_MiniFiber_RZ.size() + Hits_PSLast_RZ.size(), 1, Zall.get(), Rall.get(), Eall.get());
      int statusFitRZ = fitRZ.Eval();
      if(statusFitRZ == 1)
        {
          LocalHisto.h_RZStats->Fill("RZFit_Failed", 1.);
          continue;
        }

      TVectorD parametersRZ;
      fitRZ.GetParameters(parametersRZ); // R = [0]+Z*[1];
      LocalHisto.h_RZfit_mom->Fill(momAtFiber.Pt() / momAtFiber.Pz() - parametersRZ[1]);

      double Chi2_RZ = fitRZ.GetChisquare();
      double Ndf_RZ  = fitRZ.GetNumberFreeParameters();
      LocalHisto.h_RZfit_Chi2->Fill("chi2", Chi2_RZ, 1.);
      LocalHisto.h_RZfit_Chi2->Fill("chi2/ndf", Chi2_RZ / Ndf_RZ, 1.);
      if(Chi2_RZ / Ndf_RZ < 3)
        LocalHisto.h_RZStats->Fill("chi2Cut", 1.);

#ifdef DEBUG_CHECKRZ
      att._logger->debug("Fit in RZ from Fiber & PSCE");
      att._logger->debug("Fit chi2 : {}", Chi2_RZ);
      att._logger->debug("Fit R = [0]+Z*[1] : {} {}", parametersRZ[0], parametersRZ[1]);
#endif

      for(auto ids : id_dets_MDC)
        {
          int id_det = std::get<0>(ids);
          int id_hit = std::get<1>(ids);

          genfit::WireMeasurement* currentHit =
              dynamic_cast<genfit::WireMeasurement*>(RecoEvent.ListHits[id_det][id_hit].get());

          auto HitCoord = currentHit->getRawHitCoords();
          auto HitCov   = currentHit->getRawHitCov();

          TVector3 wire_side1(HitCoord[0], HitCoord[1], HitCoord[2]);
          TVector3 wire_side2(HitCoord[3], HitCoord[4], HitCoord[5]);
          double DriftDist    = HitCoord[6];
          double ErrDriftDist = HitCov[6][6];
          const int idMDC     = id_det - G4Sol::MG01;
          double R_size       = dl_max(idMDC) * 0.5;

          std::vector<std::array<double, 2> > list_parametersRZ;
          // std::vector<std::tuple<double,double> > WireDiffs = {{0.,0.}};
          std::vector<std::tuple<double, double> > WireDiffs = {
              {DriftDist, 0.}, {-DriftDist, 0.}, {0., DriftDist}, {0., -DriftDist}};

          for(auto WireDiff : WireDiffs)
            {
              double Rwire_side1 = wire_side1.Perp() - R_size + TMath::Sqrt(WireDiffs.size()) * std::get<1>(WireDiff);
              double Rwire_side2 = wire_side2.Perp() - R_size + TMath::Sqrt(WireDiffs.size()) * std::get<0>(WireDiff);
              double temp_slope  = (Rwire_side2 - Rwire_side1) / (wire_side2.Z() - wire_side1.Z());
              if(TMath::Abs(temp_slope) < 1e-7)
                temp_slope = 0.;
              double temp_org = Rwire_side1 - temp_slope * wire_side1.Z();
              TVectorD parameterWireRZ(2);
              parameterWireRZ(0) = temp_org;
              parameterWireRZ(1) = temp_slope;

#ifdef DEBUG_CHECKRZ
              att._logger->debug("Check in MDC {}", G4Sol::nameLiteralDet.begin()[id_det]);
              att._logger->debug("Wire : [{}, {}, {}] -> [{}, {}, {}] | R: {} - {}", wire_side1[0], wire_side1[1],
                                 wire_side1[2], wire_side2[0], wire_side2[1], wire_side2[2], wire_side1.Perp(),
                                 wire_side2.Perp());
              att._logger->debug("Wire init R= [0]+Z*[1] {} {}", temp_org, temp_slope);
              att._logger->debug("Wire R= [0]+Z*[1] {} {}", parameterWireRZ[0], parameterWireRZ[1]);
#endif

              TMatrixD systemCross(2, 2);
              systemCross(0, 0) = systemCross(1, 0) = 1.;
              systemCross(0, 1)                     = -parameterWireRZ(1);
              systemCross(1, 1)                     = -parametersRZ(1);
#ifdef DEBUG_CHECKRZ
              systemCross.Print();
#endif
              double detC        = 0;
              double maxSysCross = systemCross.Max();
              systemCross *= 1. / maxSysCross;
              systemCross.Invert(&detC);

              if(TMath::Abs(detC) < 1e-25)
                {
                  att._logger->error("E> CheckRZ : SystemCross invert failed ! {} # size {}", detC, id_dets_MDC.size());
                  att._logger->error("E> (0, 0) = (1, 0) = 1.");
                  att._logger->error("E> (0, 1) = {} . (1,1) = {}", -parameterWireRZ(1), -parametersRZ(1));
                }
#ifdef DEBUG_CHECKRZ
              att._logger->debug("Invert :");
              systemCross.Print();
#endif
              systemCross *= 1. / maxSysCross;

              TVectorD Org(2);
              Org(0)              = parameterWireRZ(0);
              Org(1)              = parametersRZ(0);
              TVectorD CrossPoint = systemCross * Org;

              std::array<double, 2> temp_RZ = {CrossPoint[1], CrossPoint[0]};
              list_parametersRZ.emplace_back(temp_RZ);
            }

          std::array<double, 2> Mean_CrossRZ = {0., 0.};
          for(auto [tempZ, tempR] : list_parametersRZ)
            {
              Mean_CrossRZ[0] += tempZ;
              Mean_CrossRZ[1] += tempR;
            }

          Mean_CrossRZ[0] *= 1. / (double)list_parametersRZ.size();
          Mean_CrossRZ[1] *= 1. / (double)list_parametersRZ.size();

          TMatrixD Cov_CrossRZ(2, 2);
          Cov_CrossRZ.Zero();
          for(auto [tempZ, tempR] : list_parametersRZ)
            {
              Cov_CrossRZ(0, 0) += (tempZ - Mean_CrossRZ[0]) * (tempZ - Mean_CrossRZ[0]);
              Cov_CrossRZ(1, 0) += (tempZ - Mean_CrossRZ[0]) * (tempR - Mean_CrossRZ[1]);
              Cov_CrossRZ(0, 1) += (tempZ - Mean_CrossRZ[0]) * (tempR - Mean_CrossRZ[1]);
              Cov_CrossRZ(1, 1) += (tempR - Mean_CrossRZ[1]) * (tempR - Mean_CrossRZ[1]);
            }

          Cov_CrossRZ *= 1. / (double)list_parametersRZ.size();

          auto HitSimulated = it_ListHitsSim->second[id_det][0];

#ifdef DEBUG_CHECKRZ
          att._logger->debug("Hit in MDC : {} {} {}, R {} | CalculatedRZ : [R, Z]: {} {}", HitSimulated.hitX,
                             HitSimulated.hitY, HitSimulated.hitZ,
                             TMath::Sqrt(HitSimulated.hitX * HitSimulated.hitX + HitSimulated.hitY * HitSimulated.hitY),
                             Mean_CrossRZ[1], Mean_CrossRZ[0]);
#endif

          LocalHisto.h_MDC_Z_residu->Fill(
              ((Mean_CrossRZ[0] -
                (correctBiasPSCE[id_det - G4Sol::MG01] + EndCap_Case * correctBiasPSEndCap[id_det - G4Sol::MG01])) -
               HitSimulated.hitZ),
              id_det - G4Sol::MG01 + EndCap_Case * 20, 1.);
          LocalHisto.h_MDC_R_residu->Fill((Mean_CrossRZ[1] - TMath::Sqrt(HitSimulated.hitX * HitSimulated.hitX +
                                                                         HitSimulated.hitY * HitSimulated.hitY)),
                                          id_det - G4Sol::MG01 + EndCap_Case * 20, 1.);
          LocalHisto.h_MDC_Z_pull->Fill(((Mean_CrossRZ[0] - (correctBiasPSCE[id_det - G4Sol::MG01] +
                                                             EndCap_Case * correctBiasPSEndCap[id_det - G4Sol::MG01])) -
                                         HitSimulated.hitZ) /
                                            TMath::Sqrt(Cov_CrossRZ(0, 0)),
                                        id_det - G4Sol::MG01 + EndCap_Case * 20, 1.);
          LocalHisto.h_MDC_R_pull->Fill((Mean_CrossRZ[1] - TMath::Sqrt(HitSimulated.hitX * HitSimulated.hitX +
                                                                       HitSimulated.hitY * HitSimulated.hitY)) /
                                            TMath::Sqrt(Cov_CrossRZ(1, 1)),
                                        id_det - G4Sol::MG01 + EndCap_Case * 20, 1.);

          if(MDCWireType == 0)
            {
              TVector3 WireDir = wire_side2 - wire_side1;
              TVector3 VecWire = WireDir;
              VecWire *= 1. / (wire_side2.Z() - wire_side1.Z());
              TVector3 PosAtWire = wire_side1 + VecWire * (Mean_CrossRZ[0] - wire_side1.Z());

              TVectorD hitCoordsNew(3);
              hitCoordsNew(0) = PosAtWire.X();
              hitCoordsNew(1) = PosAtWire.Y();
              hitCoordsNew(2) = Mean_CrossRZ[0] - (correctBiasPSCE[id_det - G4Sol::MG01] +
                                                   EndCap_Case * correctBiasPSEndCap[id_det - G4Sol::MG01]);

              TMatrixDSym hitCovNew(3);
              hitCovNew.Zero();
              hitCovNew(0, 0) = Cov_CrossRZ(1, 1);
              hitCovNew(1, 1) = Cov_CrossRZ(1, 1);
              hitCovNew(2, 2) = Cov_CrossRZ(0, 0);
              // hitCovNew(1, 2) = Cov_CrossRZ(1,0);
              // hitCovNew(2, 1) = Cov_CrossRZ(1,0);
              // hitCovNew(0, 2) = Cov_CrossRZ(1,0);
              // hitCovNew(2, 0) = Cov_CrossRZ(1,0);

              double det                = 0.;
              TMatrixDSym hitCovNewCopy = hitCovNew;
              double maxHiCovNew        = hitCovNewCopy.Max();
              hitCovNewCopy *= 1. / maxHiCovNew;
              hitCovNewCopy.Invert(&det);
              if(TMath::Abs(det) < 1.e-40)
                {
                  hitCovNewCopy *= 1. / maxHiCovNew;
                  att._logger->error("E> CheckRZ : HitCovNew for prolateSpacepoint not inversable ! det:{} max val:{}",
                                     det, maxHiCovNew);
                  att._logger->error("E> Cov : 0:{} 1:{} 2:{}", Cov_CrossRZ(1, 1), Cov_CrossRZ(1, 1),
                                     Cov_CrossRZ(0, 0));
                  hitCovNewCopy.Print();
                  att._logger->error("E> Inv:");
                  hitCovNew.Print();
                  att._logger->error(
                      "E> Hit in MDC : {} {} {}, R {} | CalculatedRZ : [R, Z]: {} {}", HitSimulated.hitX,
                      HitSimulated.hitY, HitSimulated.hitZ,
                      TMath::Sqrt(HitSimulated.hitX * HitSimulated.hitX + HitSimulated.hitY * HitSimulated.hitY),
                      Mean_CrossRZ[1], Mean_CrossRZ[0]);

                  TMatrixD Cov_CrossRZ_Err(2, 2);
                  Cov_CrossRZ_Err.Zero();
                  for(auto [tempZ, tempR] : list_parametersRZ)
                    {
                      att._logger->error("E> {} {} {} {},", tempZ, Mean_CrossRZ[0], tempR, Mean_CrossRZ[1]);
                      Cov_CrossRZ_Err(0, 0) += (tempZ - Mean_CrossRZ[0]) * (tempZ - Mean_CrossRZ[0]);
                      Cov_CrossRZ_Err(1, 0) += (tempZ - Mean_CrossRZ[0]) * (tempR - Mean_CrossRZ[1]);
                      Cov_CrossRZ_Err(0, 1) += (tempZ - Mean_CrossRZ[0]) * (tempR - Mean_CrossRZ[1]);
                      Cov_CrossRZ_Err(1, 1) += (tempR - Mean_CrossRZ[1]) * (tempR - Mean_CrossRZ[1]);
                    }
                  Cov_CrossRZ_Err.Print();
                  Cov_CrossRZ *= 1. / (double)list_parametersRZ.size();
                  Cov_CrossRZ_Err.Print();
                  att._logger->error("E> DriftDist {} ", DriftDist);
                  int indexE = 0;
                  att._logger->error("E> MiniFiber size {}, PSCE size {}", Hits_MiniFiber_RZ.size(),
                                     Hits_PSLast_RZ.size());
                  for(auto vecF : Hits_MiniFiber_RZ)
                    {
                      att._logger->error("MiniFiber : {}, {}, {}", vecF[0], vecF[1], vecF[2]);
                    }
                  for(size_t iD = 0; iD < Hits_PSLast_RZ.size(); ++iD)
                    {
                      att._logger->error("PSCE id#{} : {}, {}, {}", iD, Hits_PSLast_RZ[iD][0], Hits_PSLast_RZ[iD][1],
                                         Hits_PSLast_RZ[iD][2]);
                    }

                  att._logger->error("Fit in RZ from Fiber & PSCE");
                  att._logger->error("Fit chi2 : {}", Chi2_RZ);
                  att._logger->error("Fit R = [0]+Z*[1] : {} {}", parametersRZ[0], parametersRZ[1]);

                  for(auto WireDiff : WireDiffs)
                    {
                      double Rwire_side1 =
                          wire_side1.Perp() - R_size + TMath::Sqrt(WireDiffs.size()) * std::get<1>(WireDiff);
                      double Rwire_side2 =
                          wire_side2.Perp() - R_size + TMath::Sqrt(WireDiffs.size()) * std::get<0>(WireDiff);
                      double temp_slope = (Rwire_side2 - Rwire_side1) / (wire_side2.Z() - wire_side1.Z());
                      if(TMath::Abs(temp_slope) < 1e-7)
                        temp_slope = 0.;
                      double temp_org = Rwire_side1 - temp_slope * wire_side1.Z();
                      TVectorD parameterWireRZ(2);
                      parameterWireRZ(0) = temp_org;
                      parameterWireRZ(1) = temp_slope;
                      att._logger->error("WireDiff : #{}", indexE++);
                      att._logger->error("Check in MDC {}", G4Sol::nameLiteralDet.begin()[id_det]);
                      att._logger->error("Rwire :1:{} 2:{} R_size:{} Diff1{} Diff2{}", Rwire_side1, Rwire_side2, R_size,
                                         std::get<1>(WireDiff), std::get<0>(WireDiff));
                      att._logger->error("Wire : [{}, {}, {}] -> [{}, {}, {}] | R: {} - {}", wire_side1[0],
                                         wire_side1[1], wire_side1[2], wire_side2[0], wire_side2[1], wire_side2[2],
                                         wire_side1.Perp(), wire_side2.Perp());
                      att._logger->error("Wire init R= [0]+Z*[1] {} {}", temp_org, temp_slope);
                      att._logger->error("Wire R= [0]+Z*[1] {} {}", parameterWireRZ[0], parameterWireRZ[1]);

                      TMatrixD systemCross(2, 2);
                      systemCross(0, 0) = systemCross(1, 0) = 1.;
                      systemCross(0, 1)                     = -parameterWireRZ(1);
                      systemCross(1, 1)                     = -parametersRZ(1);

                      systemCross.Print();

                      double detC        = 0;
                      double maxSysCross = systemCross.Max();
                      systemCross *= 1. / maxSysCross;
                      systemCross.Invert(&detC);

                      if(TMath::Abs(detC) < 1e-25)
                        {
                          att._logger->error("E> CheckRZ : SystemCross invert failed ! {} # size {}", detC,
                                             id_dets_MDC.size());
                          att._logger->error("E> (0, 0) = (1, 0) = 1.");
                          att._logger->error("E> (0, 1) = {} . (1,1) = {}", -parameterWireRZ(1), -parametersRZ(1));
                        }
                      att._logger->debug("Invert :");
                      systemCross.Print();

                      systemCross *= 1. / maxSysCross;

                      TVectorD Org(2);
                      Org(0)              = parameterWireRZ(0);
                      Org(1)              = parametersRZ(0);
                      TVectorD CrossPoint = systemCross * Org;

                      std::array<double, 2> temp_RZ = {CrossPoint[1], CrossPoint[0]};
                      att._logger->error("E> temp_RZ {} {}", temp_RZ[0], temp_RZ[1]);
                      // list_parametersRZ.emplace_back(temp_RZ);
                    }
                  exit(-1);
                }

              RecoEvent.OldListHits[id_det][id_hit] = std::make_unique<genfit::ProlateSpacepointMeasurement>(
                  hitCoordsNew, hitCovNew, currentHit->getDetId(), currentHit->getHitId(), nullptr);
              dynamic_cast<genfit::ProlateSpacepointMeasurement*>(RecoEvent.OldListHits[id_det][id_hit].get())
                  ->setLargestErrorDirection(WireDir);
            }
          else if(MDCWireType == 1)
            {
              TVectorD hitCoordsNew(8);
              hitCoordsNew(0) = wire_side1[0];
              hitCoordsNew(1) = wire_side1[1];
              hitCoordsNew(2) = wire_side1[2];
              hitCoordsNew(3) = wire_side2[0];
              hitCoordsNew(4) = wire_side2[1];
              hitCoordsNew(5) = wire_side2[2];
              hitCoordsNew(6) = DriftDist;
              hitCoordsNew(7) =
                  Mean_CrossRZ[0] -
                  (correctBiasPSCE[id_det - G4Sol::MG01] + EndCap_Case * correctBiasPSEndCap[id_det - G4Sol::MG01]) -
                  wire_side1[2];
              TMatrixDSym hitCovNew(8);
              hitCovNew(6, 6) = ErrDriftDist;
              hitCovNew(7, 7) = Cov_CrossRZ(0, 0);

              RecoEvent.OldListHits[id_det][id_hit] = std::make_unique<genfit::WirePointMeasurement>(
                  hitCoordsNew, hitCovNew, currentHit->getDetId(), currentHit->getHitId(), nullptr);
            }

          RecoEvent.OldListHits[id_det][id_hit].swap(RecoEvent.ListHits[id_det][id_hit]);
        }

      LocalHisto.h_RZStats->Fill("RZ_Succeed", 1.);

      if(id_dets_PSEndcap.size() != 0)
        {
          int index_RZ = 0;
          for(auto ids : id_dets_PSEndcap)
            {
              int id_det = std::get<0>(ids);
              int id_hit = std::get<1>(ids);

              genfit::PlanarMeasurement* currentHit =
                  dynamic_cast<genfit::PlanarMeasurement*>(RecoEvent.ListHits[id_det][id_hit].get());
              auto dummyState   = genfit::StateOnPlane();
              auto currentPlane = currentHit->constructPlane(genfit::StateOnPlane());
              auto HitCoord     = currentHit->getRawHitCoords();
              auto HitCov       = currentHit->getRawHitCov();
              // std::cout<<"f_extractHitLabs :\n";

              const double phi_max = 7.5 * TMath::DegToRad();
              TVector3 Orig        = currentPlane->getO();
              TVector3 V           = currentPlane->getV();
              TVector3 U           = currentPlane->getU();
              genfit::SharedPlanePtr planeNew(new genfit::DetPlane(Orig, U, V));

              double R_est    = Hits_PSLast_RZ[index_RZ][0];
              double RCov_est = Hits_PSLast_RZ[index_RZ][2];
              double Phi_Uni  = HitCoord(0) * TMath::DegToRad();

              TVectorD hitCoordsNew(2);
              hitCoordsNew(0) = R_est * TMath::Sin(Phi_Uni); // u
              hitCoordsNew(1) = R_est;                       // v

              TMatrixDSym hitCovNew(2);
              hitCovNew.Zero();
              hitCovNew(0, 0) = TMath::Sq(2 * R_est * TMath::Sin(0.5 * phi_max)) / 12.;
              hitCovNew(1, 1) = RCov_est;

              RecoEvent.OldListHits[id_det][id_hit] = std::make_unique<genfit::PlanarMeasurement>(
                  hitCoordsNew, hitCovNew, currentHit->getDetId(), currentHit->getHitId(), nullptr);
              dynamic_cast<genfit::PlanarMeasurement*>(RecoEvent.OldListHits[id_det][id_hit].get())->setPlane(planeNew);

              RecoEvent.OldListHits[id_det][id_hit].swap(RecoEvent.ListHits[id_det][id_hit]);
            }
        }

      if(ChangeMiniFiber)
        {
          TVector3 Xdir(1., 0., 0.), Ydir(0., 1., 0.);
          for(auto ids : id_dets_MiniFiber)
            {
              int id_det = std::get<0>(ids);
              int id_hit = std::get<1>(ids);

              genfit::PlanarMeasurement* currentHit =
                  dynamic_cast<genfit::PlanarMeasurement*>(RecoEvent.ListHits[id_det][id_hit].get());

              auto dummyState   = genfit::StateOnPlane();
              auto currentPlane = currentHit->constructPlane(genfit::StateOnPlane());
              auto currentO     = currentPlane->getO();

              double z = currentO.Z();
              TVector3 o(0., 0., z);

              double XPosAtZ = miniTrack(2) * z + miniTrack(0);
              double YPosAtZ = miniTrack(3) * z + miniTrack(1);

              genfit::SharedPlanePtr planeNew(new genfit::DetPlane(o, Xdir, Ydir));

              TVectorD hitCoordsNew(2);
              hitCoordsNew(0) = XPosAtZ;
              hitCoordsNew(1) = YPosAtZ;

              TMatrixDSym hitCovNew(2);
              hitCovNew(0, 0) = miniTrackCov(0, 0);
              hitCovNew(0, 1) = miniTrackCov(0, 1);
              hitCovNew(1, 0) = miniTrackCov(1, 0);
              hitCovNew(1, 1) = miniTrackCov(1, 1);

              RecoEvent.OldListHits[id_det][id_hit] = std::make_unique<genfit::PlanarMeasurement>(
                  hitCoordsNew, hitCovNew, currentHit->getDetId(), currentHit->getHitId(), nullptr);
              dynamic_cast<genfit::PlanarMeasurement*>(RecoEvent.OldListHits[id_det][id_hit].get())->setPlane(planeNew);

              RecoEvent.OldListHits[id_det][id_hit].swap(RecoEvent.ListHits[id_det][id_hit]);
            }
        }
    }

  return 0;
}

template class TCheckRZ<MCAnaEventG4Sol>;
template class TCheckRZ<Ana_WasaEvent>;
