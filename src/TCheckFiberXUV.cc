#include "TCheckFiberXUV.h"

//#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"
//#include "FullRecoEvent.hh"
//#include "KalmanFittedStateOnPlane.h"
//#include "KalmanFitterInfo.h"
//#include "PlanarMeasurement.h"
#include "AbsMeasurement.h"
//#include "ReturnRes.hh"
#include "StateOnPlane.h"

//#include <numeric>
//#include <set>
#include <sstream>
//#include <tuple>
//#include <list>

using namespace std;
using namespace G4Sol;

template<class Out>
TCheckFiberXUV<Out>::TCheckFiberXUV(const THyphiAttributes& attribut) : TDataProcessInterface<Out>("CheckFiberXUV"), att(attribut)
{
	cut_d = att.FiberXUV_cutd;
}

template<class Out>
TCheckFiberXUV<Out>::~TCheckFiberXUV()
{

}

template<class Out>
void TCheckFiberXUV<Out>::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }

template<class Out>
ReturnRes::InfoM TCheckFiberXUV<Out>::operator()(FullRecoEvent& RecoEvent, Out* OutTree)
{

  int result_finder = Exec(RecoEvent, OutTree);

  return SoftExit(result_finder);
}

template<class Out>
int TCheckFiberXUV<Out>::Exec(FullRecoEvent& RecoEvent, Out* OutTree)
{
  return CheckHitXUV(RecoEvent);
}

template<class Out>
ReturnRes::InfoM TCheckFiberXUV<Out>::SoftExit(int result_full) { return ReturnRes::Fine; }

template<class Out>
void TCheckFiberXUV<Out>::SelectHists()
{
  for(size_t i = 0; i < 7; ++i)
	  {
	  	LocalHisto.h_ResidualFiberHitX[i]       = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitX[i]);
	  	LocalHisto.h_ResidualFiberHitY[i]       = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitY[i]);
	  	LocalHisto.h_ResidualFiberHitXY[i]      = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitXY[i]);
	  	LocalHisto.h_ResidualFiberHitX_Angle[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitX_Angle[i]);
	  	LocalHisto.h_ResidualFiberHitY_Angle[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitY_Angle[i]);
	  }
}

template<class Out>
int TCheckFiberXUV<Out>::CheckHitXUV(const FullRecoEvent& RecoEvent)
{

/*
	std::vector<genfit::PlanarMeasurement*> 


  886            genfit::PlanarMeasurement* currentHit1 =
  887:               dynamic_cast<genfit::WireMeasurement*>(RecoEvent.ListHits[layerRedo][id_hit].get());
  888  
  889            auto HitCoord = currentHit1->getRawHitCoords();
*/

/*
  std::vector< std::vector< std::vector<FiberHitAna*> > > FiberHitCont;
  std::vector< std::vector< std::vector<FiberHitAna*> > > FiberHitClCont;
  for(int i=0; i<7; ++i){
    std::vector<FiberHitAna*> buf_v;
    std::vector< std::vector<FiberHitAna*> > buf_vv;
    for(int j=0; j<3; ++j){
      buf_vv.emplace_back(buf_v);
    }
    FiberHitCont.emplace_back(buf_vv);
    FiberHitClCont.emplace_back(buf_vv);
  }					  



    // make container
    for(int i=0; i<(int)s2fiber->fiberhit.size(); ++i){
      FiberHitAna *hit_ana = new FiberHitAna(s2fiber->fiberhit[i], par, t_r);
      if(!hit_ana->IsValid()){ delete hit_ana; continue; }
      FiberHitCont[hit_ana->GetDet()][hit_ana->GetLay()].emplace_back(hit_ana);
    }


    FiberAnalyzer *fiberana = new FiberAnalyzer();
    FiberHitClCont = fiberana->Clusterize(FiberHitCont);


    std::vector< std::vector<FiberHitXUV*> > FiberXUVCont = fiberana->FindHit(FiberHitClCont, par);
    for(int i=0; i<7; ++i){
      for(int j=0; j<(int)FiberXUVCont[i].size(); ++j){
        h16[i]->Fill(FiberXUVCont[i][j]->GetPosX(), FiberXUVCont[i][j]->GetPosY());
        h17_2[i]->Fill(FiberXUVCont[i][j]->GetD());
      }
      h17[i]->Fill( FiberXUVCont[i].size() );
    }

*/

	std::vector<int> id_detector = {G4Sol::FiberD1_x, G4Sol::FiberD2_x, G4Sol::FiberD3_x, G4Sol::MiniFiberD1_x1, G4Sol::MiniFiberD1_x2, G4Sol::FiberD4_x, G4Sol::FiberD5_x};
  
  //Get recons hits
  for(int tmp_det = 0; tmp_det < id_detector.size(); ++tmp_det)
	  {
	  	int i = id_detector[tmp_det];
	  	if((i == G4Sol::FiberD1_x) || (i == G4Sol::FiberD2_x)) //Skip UFT12
	  		continue;
	  	if((RecoEvent.ListHits[i].size() != 1) || (RecoEvent.ListHits[i+1].size() != 1) || (RecoEvent.ListHits[i+2].size() != 1)) //Only 1-hit events
	  		continue;

	  	TVector2 tmp_HitXY(-9999., -9999.);
	  	FindHitXUV(RecoEvent.ListHits[i][0].get(), RecoEvent.ListHits[i+1][0].get(), RecoEvent.ListHits[i+2][0].get(), tmp_det, &tmp_HitXY);
	  	vect_HitXY[tmp_det].emplace_back(tmp_HitXY);
	  }

	//Get real hits
  for(int tmp_det = 0; tmp_det < id_detector.size(); ++tmp_det)
	  {
	  	int i = id_detector[tmp_det];
	  	for(auto it_track : RecoEvent.TrackDAFSim) //Get real hits
		  	{
		  		auto it_tmpFiber = it_track.second[i+1];
				  if(it_tmpFiber.size() == 0)
			      continue;

			    if(it_tmpFiber[0].momZ < 0.)
			      continue;

			    double meanPosX = 0.;
			    double meanPosY = 0.;
			    double energy = 0.;

			    for(size_t idX = 0; idX < it_tmpFiber.size(); ++idX)
				    {
			        if(it_tmpFiber[idX].Eloss < EnergyThreshold)
			        	continue;

		          meanPosX += it_tmpFiber[idX].hitX * it_tmpFiber[idX].Eloss;
		          meanPosY += it_tmpFiber[idX].hitY * it_tmpFiber[idX].Eloss;
		          energy   += it_tmpFiber[idX].Eloss;
			      }

		      meanPosX /= energy;
		      meanPosY /= energy;

	  			TVector2 tmp_realHitXY(meanPosX, meanPosY);
	  			vect_realHitXY[tmp_det].emplace_back(tmp_realHitXY);

			    double theta = TMath::ATan2( std::sqrt(it_tmpFiber[0].momX * it_tmpFiber[0].momX + it_tmpFiber[0].momY * it_tmpFiber[0].momY), 
			    															it_tmpFiber[0].momZ) * 180. / M_PI;
			    vect_realHitAng[tmp_det].emplace_back(theta);
		  	} 	
		}

	//Compare real-recons hits
	for(size_t i = 0; i < id_detector.size(); ++i)
		{
			if(vect_HitXY[i].size() == 0)
				continue;

			double Hitdistance = 9999.;
			size_t closest_j;
			for(size_t j = 0; j < vect_realHitXY[i].size(); ++j)
				{
					double tmp_Hitdistance = (vect_realHitXY[i][j] - vect_HitXY[i][0]).Mod();
					if(tmp_Hitdistance < Hitdistance)
						{
							Hitdistance = tmp_Hitdistance;
							closest_j = j;
						}
				}

			TVector2 resHit = vect_realHitXY[i][closest_j] - vect_HitXY[i][0];

      LocalHisto.h_ResidualFiberHitX[i]      ->Fill(resHit.X(), 1.);
      LocalHisto.h_ResidualFiberHitY[i]      ->Fill(resHit.Y(), 1.);
      LocalHisto.h_ResidualFiberHitXY[i]     ->Fill(resHit.X(), resHit.Y(), 1.);
      LocalHisto.h_ResidualFiberHitX_Angle[i]->Fill(resHit.X(), vect_realHitAng[i][0], 1.);
      LocalHisto.h_ResidualFiberHitY_Angle[i]->Fill(resHit.Y(), vect_realHitAng[i][0], 1.);
		}

  return 0;
}


template<class Out>
void TCheckFiberXUV<Out>::FindHitXUV(genfit::AbsMeasurement* hitx, genfit::AbsMeasurement* hitu, genfit::AbsMeasurement* hitv, int id_det, TVector2* HitXY)
{
	double ang_uv = ang[id_det][1] * M_PI / 180.;

	auto vpos_x = hitx->getRawHitCoords();
	auto vpos_u = hitu->getRawHitCoords();
	auto vpos_v = hitv->getRawHitCoords();

	double pos_x = vpos_x[0];
	double pos_u = vpos_u[0];
	double pos_v = vpos_v[0];

  pos_u = pos_u / std::cos(ang_uv);
  pos_v = pos_v / std::cos(ang_uv);

  double hit_xu_x = pos_x;

  double hit_xu_y = 999.;
  if(ang_uv>0) hit_xu_y = -std::tan((M_PI/2. - ang_uv))*pos_x + std::tan((M_PI/2. -ang_uv))*pos_u;
  else         hit_xu_y =  std::tan((M_PI/2. + ang_uv))*pos_x - std::tan((M_PI/2. +ang_uv))*pos_u;

  double hit_xv_y = -999.;
  if(ang_uv>0) hit_xv_y =  std::tan((M_PI/2. - ang_uv))*pos_x - std::tan((M_PI/2. - ang_uv))*pos_v;
  else         hit_xv_y = -std::tan((M_PI/2. + ang_uv))*pos_x + std::tan((M_PI/2. + ang_uv))*pos_v;

	if((std::fabs(hit_xu_y-999.)<1.)||(std::fabs(hit_xv_y+999.)<1.)) return;
	if(std::fabs(hit_xu_y-hit_xv_y) > cut_d[id_det]) return;

  double buf_x = hit_xu_x + (hit_xu_y-hit_xv_y)/2.*std::tan(ang_uv/2.);
  double buf_y = (hit_xu_y+hit_xv_y)/2.;

  HitXY->Set(buf_x, buf_y);
}


template class TCheckFiberXUV<MCAnaEventG4Sol>;
template class TCheckFiberXUV<Ana_WasaEvent>;
