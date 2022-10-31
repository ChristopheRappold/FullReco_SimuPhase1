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
#include <set>
#include <sstream>
//#include <tuple>
//#include <list>

//#define DEBUG_FIBERXUV
//#define SINGLE_FIBERXUV

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
	  	LocalHisto.h_ResidualFiberHitX[i]       			= this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitX[i]);
	  	LocalHisto.h_ResidualFiberHitY[i]       			= this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitY[i]);
	  	LocalHisto.h_ResidualFiberHitR[i]       			= this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitR[i]);
	  	LocalHisto.h_ResidualFiberHitXY[i]      			= this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitXY[i]);
	  	LocalHisto.h_ResidualFiberHitX_Angle[i] 			= this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitX_Angle[i]);
	  	LocalHisto.h_ResidualFiberHitY_Angle[i] 			= this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitY_Angle[i]);
	  	LocalHisto.h_ResidualFiberHitX_HitX[i] 		  	= this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitX_HitX[i]);
	  	LocalHisto.h_ResidualFiberHitY_HitY[i] 		  	= this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitY_HitY[i]);
	  	LocalHisto.h_ResidualFiberHitR_Angle[i]				= this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitR_Angle[i]);
	  	LocalHisto.h_EfficiencyFiberHit[i]      			= this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EfficiencyFiberHit[i]);
	  	LocalHisto.h_ResidualSingleFiberHitX[i]       = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualSingleFiberHitX[i]);
	  	LocalHisto.h_ResidualSingleFiberHitY[i]       = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualSingleFiberHitY[i]);
	  	LocalHisto.h_ResidualSingleFiberHitR[i]       = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualSingleFiberHitR[i]);
	  	LocalHisto.h_ResidualSingleFiberHitX_Angle[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualSingleFiberHitX_Angle[i]);
	  	LocalHisto.h_ResidualSingleFiberHitY_Angle[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualSingleFiberHitY_Angle[i]);
	  	LocalHisto.h_ResidualSingleFiberHitR_Angle[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualSingleFiberHitR_Angle[i]);
	  }
}

template<class Out>
int TCheckFiberXUV<Out>::CheckHitXUV(const FullRecoEvent& RecoEvent)
{
	//Get 3-fiber correct combination & real hit position
  for(int i_det = 0; i_det < id_detector.size(); ++i_det)
	  {
	  	int i = id_detector[i_det];
		  std::unordered_map<int,std::tuple<std::vector<size_t>,std::vector<size_t>,std::vector<size_t> > > tmp_realCombIdHit;
		  vect_realCombIdHit.emplace_back(tmp_realCombIdHit);

	  	if((i == G4Sol::FiberD1_x) || (i == G4Sol::FiberD2_x)) //Skip UFT12
		  		continue;

		  int imid = id_mid[i_det];

	  	for(std::pair<int, std::vector<std::vector<SimHit> > > it_track : RecoEvent.TrackDAFSim)
		    {
		      const int id_track = it_track.first;

			  	double tmp_hit_x = 0.;
			  	double tmp_hit_u = 0.;
			  	double tmp_hit_v = 0.;

			  	int counter_x = 0;
			  	int counter_u = 0;
			  	int counter_v = 0;

			  	std::vector<size_t> tmp_idhit_x = {};
			  	std::vector<size_t> tmp_idhit_u = {};
			  	std::vector<size_t> tmp_idhit_v = {};

		      for(size_t j = 0; j < RecoEvent.ListHits[i].size(); ++j)
			      {
							if(RecoEvent.ListHitsToTracks[i][j] == id_track)
								{
									TVectorD& vpos_x = RecoEvent.ListHits[i][j]->getRawHitCoords();
									tmp_hit_x += vpos_x[0];
									++counter_x;
									tmp_idhit_x.emplace_back(j);
									//std::cout << "ListHits x index j: " << j << "\n";
								}
						}

		      for(size_t j = 0; j < RecoEvent.ListHits[i+1].size(); ++j)
			      {
							if(RecoEvent.ListHitsToTracks[i+1][j] == id_track)
								{
									TVectorD& vpos_u = RecoEvent.ListHits[i+1][j]->getRawHitCoords();
									tmp_hit_u += vpos_u[0];
									++counter_u;
									tmp_idhit_u.emplace_back(j);
									//std::cout << "ListHits u index j: " << j << "\n";
								}
						}

		      for(size_t j = 0; j < RecoEvent.ListHits[i+2].size(); ++j)
			      {
							if(RecoEvent.ListHitsToTracks[i+2][j] == id_track)
								{
									TVectorD& vpos_v = RecoEvent.ListHits[i+2][j]->getRawHitCoords();
									tmp_hit_v += vpos_v[0];
									++counter_v;
									tmp_idhit_v.emplace_back(j);
									//std::cout << "ListHits v index j: " << j << "\n";
								}
						}

					if((counter_x == 0) || (counter_u == 0) || (counter_v == 0))
						continue;

					tmp_hit_x /= counter_x;
					tmp_hit_u /= counter_u;
					tmp_hit_v /= counter_v;

					vect_realCombHit[i_det].emplace_back(std::make_tuple(tmp_hit_x, tmp_hit_u, tmp_hit_v));
					vect_realCombIdHit[i_det].insert(std::make_pair(id_track, std::make_tuple(tmp_idhit_x, tmp_idhit_u, tmp_idhit_v)));

		  		std::vector<SimHit> it_tmpFiber = it_track.second[i+imid];
			    if(it_tmpFiber[0].momZ < 0.)
			      continue;

			    double meanPosX = 0.;
			    double meanPosY = 0.;
			    double energy = 0.;

			    for(size_t j = 0; j < it_tmpFiber.size(); ++j)
				    {
		          meanPosX += it_tmpFiber[j].hitX * it_tmpFiber[j].Eloss;
		          meanPosY += it_tmpFiber[j].hitY * it_tmpFiber[j].Eloss;
		          energy   += it_tmpFiber[j].Eloss;
			      }

		      meanPosX /= energy;
		      meanPosY /= energy;

	  			TVector2 tmp_realHitXY(meanPosX, meanPosY);

			    double theta = TMath::ATan2( std::sqrt(it_tmpFiber[0].momX * it_tmpFiber[0].momX + it_tmpFiber[0].momY * it_tmpFiber[0].momY), 
			    															it_tmpFiber[0].momZ) * 180. / M_PI;

	  			vect_realHitXYAngId[i_det].emplace_back(std::make_tuple(tmp_realHitXY, theta, id_track));
		  	}
		}


  //Get recons hits
  for(int i_det = 0; i_det < id_detector.size(); ++i_det)
	  {
	  	int i = id_detector[i_det];
	  	if((i == G4Sol::FiberD1_x) || (i == G4Sol::FiberD2_x)) //Skip UFT12
	  		continue;

#ifdef DEBUG_FIBERXUV
		  std::cout << "Before Clusterize()\n";
#endif

		  std::vector<genfit::AbsMeasurement*> hitx = Clusterize(RecoEvent.ListHits[i]);
		  std::vector<genfit::AbsMeasurement*> hitu = Clusterize(RecoEvent.ListHits[i+1]);
		  std::vector<genfit::AbsMeasurement*> hitv = Clusterize(RecoEvent.ListHits[i+2]);

#ifdef DEBUG_FIBERXUV
		  std::cout << "Before FindHitXUV()\n";
#endif

	  	FindHitXUV(hitx, hitu, hitv, i_det);

#ifdef SINGLE_FIBERXUV
	#ifdef DEBUG_FIBERXUV
		  std::cout << "Before FindSingleHitXUVId()\n";
	#endif

	  	for(std::pair<int,std::tuple<std::vector<size_t>,std::vector<size_t>,std::vector<size_t>>> it_CombHit : vect_realCombIdHit[i_det])
		  	{
		  		const int id_track = it_CombHit.first;

		  		std::vector<genfit::AbsMeasurement*> single_hitx = {};
		  		std::vector<genfit::AbsMeasurement*> single_hitu = {};
		  		std::vector<genfit::AbsMeasurement*> single_hitv = {};

#ifdef DEBUG_FIBERXUV
		  		std::cout << "(std::get<0>(it_CombHit.second)).size(): " << (std::get<0>(it_CombHit.second)).size() << "\n";
		  		std::cout << "(std::get<1>(it_CombHit.second)).size(): " << (std::get<1>(it_CombHit.second)).size() << "\n";
		  		//std::cout << "(std::get<2>(it_CombHit.second)).size(): " << (std::get<2>(it_CombHit.second)).size() << "\n";

		  		std::cout << "RecoEvent.ListHits[i].size(): "   << RecoEvent.ListHits[i].size() << "\n";
		  		std::cout << "RecoEvent.ListHits[i+1].size(): " << RecoEvent.ListHits[i+1].size() << "\n";
		  		//std::cout << "RecoEvent.ListHits[i+2].size(): " << RecoEvent.ListHits[i+2].size() << "\n";
#endif

		  		for(size_t j = 0; j < (std::get<0>(it_CombHit.second)).size(); ++j)
		  		{
		  			std::cout << "Inside x loop: j = " << j << "\n";
			  		std::cout << "std::get<0>(it_CombHit.second)[j]: " << std::get<0>(it_CombHit.second)[j] << "\n";
						single_hitx.emplace_back(RecoEvent.ListHits[i][(std::get<0>(it_CombHit.second))[j]]->clone());
		  		}

		  		for(size_t j = 0; j < (std::get<1>(it_CombHit.second)).size(); ++j)
		  		{
		  			std::cout << "Inside u loop: j = " << j << "\n";
			  		std::cout << "std::get<1>(it_CombHit.second)[j]: " << std::get<1>(it_CombHit.second)[j] << "\n";
						single_hitu.emplace_back(RecoEvent.ListHits[i+1][(std::get<1>(it_CombHit.second))[j]]->clone());
		  		}

		  		for(size_t j = 0; j < (std::get<2>(it_CombHit.second)).size(); ++j)
		  		{
		  			//std::cout << "Inside v loop: j = " << j << "\n";
			  		//std::cout << "std::get<2>(it_CombHit.second)[j]: " << std::get<2>(it_CombHit.second)[j] << "\n";
						single_hitv.emplace_back(RecoEvent.ListHits[i+2][(std::get<2>(it_CombHit.second))[j]]->clone());
		  		}

		  		std::vector<genfit::AbsMeasurement*> cluster_single_hitx = Clusterize(single_hitx);
		  		std::vector<genfit::AbsMeasurement*> cluster_single_hitu = Clusterize(single_hitu);
		  		std::vector<genfit::AbsMeasurement*> cluster_single_hitv = Clusterize(single_hitv);

		  		FindSingleHitXUVId(cluster_single_hitx, cluster_single_hitu, cluster_single_hitv, i_det, id_track);
		  	}

#endif

	  }


#ifdef DEBUG_FIBERXUV
	std::cout << "Before get real hits\n";
#endif


#ifdef DEBUG_FIBERXUV
	std::cout << "Before compare real-recons hits\n";
#endif

	//Compare real-recons hits
	for(size_t i_det = 0; i_det < id_detector.size(); ++i_det)
		{

#ifdef DEBUG_FIBERXUV
			std::cout << "Fiber: " << i_det << "  NumRealHit: " << vect_realHitXYAngId[i_det].size() << "   NumRecoHit: " << vect_HitXY[i_det].size() << "\n";
#endif

			for(size_t j = 0; j < vect_realHitXYAngId[i_det].size(); ++j)
				{
					if(vect_HitXY[i_det].size() == 0)
						break;

					double Hitdistance = 99999.;
					size_t closest_i;
					for(size_t i = 0; i < vect_HitXY[i_det].size(); ++i)
						{
							double tmp_Hitdistance = (std::get<0>(vect_realHitXYAngId[i_det][j]) - vect_HitXY[i_det][i]).Mod();
							if(tmp_Hitdistance < Hitdistance)
								{
									Hitdistance = tmp_Hitdistance;
									closest_i = i;
								}
						}

					TVector2 resHit = std::get<0>(vect_realHitXYAngId[i_det][j]) - vect_HitXY[i_det][closest_i];

		      LocalHisto.h_ResidualFiberHitX[i_det]      ->Fill(resHit.X(), 1.);
		      LocalHisto.h_ResidualFiberHitY[i_det]      ->Fill(resHit.Y(), 1.);
		      LocalHisto.h_ResidualFiberHitR[i_det]      ->Fill(resHit.Mod(), 1.);
		      LocalHisto.h_ResidualFiberHitXY[i_det]     ->Fill(resHit.X(), resHit.Y(), 1.);
		      LocalHisto.h_ResidualFiberHitX_Angle[i_det]->Fill(std::get<1>(vect_realHitXYAngId[i_det][j]), resHit.X(), 1.);
		      LocalHisto.h_ResidualFiberHitY_Angle[i_det]->Fill(std::get<1>(vect_realHitXYAngId[i_det][j]), resHit.Y(), 1.);
		      LocalHisto.h_ResidualFiberHitX_HitX[i_det] ->Fill(std::get<0>(vect_realHitXYAngId[i_det][j]).X(), resHit.X(), 1.);
		      LocalHisto.h_ResidualFiberHitY_HitY[i_det] ->Fill(std::get<0>(vect_realHitXYAngId[i_det][j]).Y(), resHit.Y(), 1.);
		      LocalHisto.h_ResidualFiberHitR_Angle[i_det]->Fill(std::get<1>(vect_realHitXYAngId[i_det][j]), resHit.Mod(), 1.);
				}

			vect_HitXY[i_det].clear();

#ifdef SINGLE_FIBERXUV

			for(size_t i = 0; i < vect_SingleHitXYId[i_det].size(); ++i)
				{
					const int id_track = std::get<1>(vect_SingleHitXYId[i_det][i]);

					for(size_t j = 0; j < vect_realHitXYAngId[i_det].size(); ++j)
						{
							if(std::get<2>(vect_realHitXYAngId[i_det][j]) != id_track)
								continue;

							TVector2 resSingleHit = std::get<0>(vect_realHitXYAngId[i_det][j]) - std::get<0>(vect_SingleHitXYId[i_det][i]);

							LocalHisto.h_ResidualSingleFiberHitX[i_det]      ->Fill(resSingleHit.X(), 1.);
				      LocalHisto.h_ResidualSingleFiberHitY[i_det]      ->Fill(resSingleHit.Y(), 1.);
				      LocalHisto.h_ResidualSingleFiberHitR[i_det]      ->Fill(resSingleHit.Mod(), 1.);
				      LocalHisto.h_ResidualSingleFiberHitX_Angle[i_det]->Fill(std::get<1>(vect_realHitXYAngId[i_det][j]), resSingleHit.X(), 1.);
				      LocalHisto.h_ResidualSingleFiberHitY_Angle[i_det]->Fill(std::get<1>(vect_realHitXYAngId[i_det][j]), resSingleHit.Y(), 1.);
				      LocalHisto.h_ResidualSingleFiberHitR_Angle[i_det]->Fill(std::get<1>(vect_realHitXYAngId[i_det][j]), resSingleHit.Mod(), 1.);
						}
				}

			vect_SingleHitXYId[i_det].clear();

#endif

			vect_realHitXYAngId[i_det].clear();

			for(size_t i = 0; i < vect_realCombHit[i_det].size(); ++i)
				{
					LocalHisto.h_EfficiencyFiberHit[i_det]->Fill(0., 1.);

					for(size_t j = 0; j < vect_CombHit[i_det].size(); ++j)
						{
							double res_x = std::fabs(std::get<0>(vect_realCombHit[i_det][i]) - std::get<0>(vect_CombHit[i_det][j]));
							double res_u = std::fabs(std::get<1>(vect_realCombHit[i_det][i]) - std::get<1>(vect_CombHit[i_det][j]));
							double res_v = std::fabs(std::get<2>(vect_realCombHit[i_det][i]) - std::get<2>(vect_CombHit[i_det][j]));

							if((res_x < fiber_width) && (res_u < fiber_width) && (res_v < fiber_width))
							{
								LocalHisto.h_EfficiencyFiberHit[i_det]->Fill(1., 1.);
								break;
							}
						}
				}

			vect_realCombHit[i_det].clear();
			vect_CombHit[i_det].clear();
		}

#ifdef DEBUG_FIBERXUV
	std::cout << "TCheckFiberXUV End\n";
#endif

  return 0;
}


template<class Out>
void TCheckFiberXUV<Out>::FindHitXUV(const std::vector<genfit::AbsMeasurement*>& hitx,
																			const std::vector<genfit::AbsMeasurement*>& hitu,
																			const std::vector<genfit::AbsMeasurement*>& hitv,
																				int id_det)
{
  std::set<int> flag_x;
  std::set<int> flag_u;
  std::set<int> flag_v;

  int nx = hitx.size();
  int nu = hitu.size();
  int nv = hitv.size();
	if((nx==0)||(nu==0)||(nv==0))
		return;

	double ang_u = ang[id_det][id_mid[id_det]] * M_PI / 180.;

	for(int i = 0; i < nx ; ++i)
	{
		if(flag_x.size()>0 && flag_x.find(i)!=flag_x.end())
			continue;

		TVectorD& vpos_x = hitx[i]->getRawHitCoords();
		double pos_x = vpos_x[0];

	  double hit_xu_x = pos_x;
	  double hit_xu_y = 999.;
	  double hit_xv_y = -999.;

	  double tmp_hit_u;
	  double tmp_hit_v;

	  int used_j = -999;
	  int used_k = -999;

	  for(int j = 0; j < nu ; ++j)
		  {
        if(flag_u.size()>0 && flag_u.find(j)!=flag_u.end())
        	continue;
				
				TVectorD& vpos_u = hitu[j]->getRawHitCoords();
				double pos_u = vpos_u[0];
				tmp_hit_u = pos_u;
				pos_u = pos_u / std::cos(ang_u);

				double tmp_hit_xu_y = 999.;
				if(ang_u>0.)
					tmp_hit_xu_y = -std::tan((M_PI/2. - ang_u))*pos_x + std::tan((M_PI/2. -ang_u))*pos_u;
				else
					tmp_hit_xu_y =  std::tan((M_PI/2. + ang_u))*pos_x - std::tan((M_PI/2. +ang_u))*pos_u;

        for(int k = 0; k < nv; ++k)
        	{
          	if(flag_v.size()>0 && flag_v.find(k)!=flag_v.end())
          		continue;

          	TVectorD& vpos_v = hitv[k]->getRawHitCoords();
						double pos_v = vpos_v[0];
						tmp_hit_v = pos_v;
					  pos_v = pos_v / std::cos(ang_u);

					  double tmp_hit_xv_y = -999.;
					  if(ang_u>0)
					  	tmp_hit_xv_y =  std::tan((M_PI/2. - ang_u))*pos_x - std::tan((M_PI/2. - ang_u))*pos_v;
					  else
					  	tmp_hit_xv_y = -std::tan((M_PI/2. + ang_u))*pos_x + std::tan((M_PI/2. + ang_u))*pos_v;

					  if(std::fabs(hit_xu_y-hit_xv_y) > std::fabs(tmp_hit_xu_y-tmp_hit_xv_y))
						  {
		            hit_xu_y = tmp_hit_xu_y;
		            hit_xv_y = tmp_hit_xv_y;
		            used_j = j;
		            used_k = k;
						  }
	  			}
  		}

		if((std::fabs(hit_xu_y-999.)<1.)||(std::fabs(hit_xv_y+999.)<1.))
			continue;
		if(std::fabs(hit_xu_y-hit_xv_y) > cut_d[id_det])
			continue;

	  double buf_x = hit_xu_x + (hit_xu_y-hit_xv_y)/2.*std::tan(ang_u/2.);
	  double buf_y = (hit_xu_y+hit_xv_y)/2.;

	  TVector2 tmp_hitxy(buf_x, buf_y);
	  vect_HitXY[id_det].emplace_back(tmp_hitxy);

	  vect_CombHit[id_det].emplace_back(std::make_tuple(pos_x, tmp_hit_u, tmp_hit_v));

	  flag_x.insert(i);
	  flag_u.insert(used_j);
	  flag_v.insert(used_k);
	}

	return;
}


template<class Out>
void TCheckFiberXUV<Out>::FindSingleHitXUVId(const std::vector<genfit::AbsMeasurement*>& hitx,
																							const std::vector<genfit::AbsMeasurement*>& hitu,
																							const std::vector<genfit::AbsMeasurement*>& hitv,
																								int id_det, int id_track)
{
  std::set<int> flag_x;
  std::set<int> flag_u;
  std::set<int> flag_v;

  int nx = hitx.size();
  int nu = hitu.size();
  int nv = hitv.size();
	if((nx==0)||(nu==0)||(nv==0))
		return;

	double ang_u = ang[id_det][1] * M_PI / 180.;

	for(int i = 0; i < nx ; ++i)
	{
		if(flag_x.size()>0 && flag_x.find(i)!=flag_x.end())
			continue;

		TVectorD& vpos_x = hitx[i]->getRawHitCoords();
		double pos_x = vpos_x[0];

	  double hit_xu_x = pos_x;
	  double hit_xu_y = 999.;
	  double hit_xv_y = -999.;

	  int used_j = -999;
	  int used_k = -999;

	  for(int j = 0; j < nu ; ++j)
		  {
        if(flag_u.size()>0 && flag_u.find(j)!=flag_u.end())
        	continue;
				
				TVectorD& vpos_u = hitu[j]->getRawHitCoords();
				double pos_u = vpos_u[0];
				pos_u = pos_u / std::cos(ang_u);

				double tmp_hit_xu_y = 999.;
				if(ang_u>0.)
					tmp_hit_xu_y = -std::tan((M_PI/2. - ang_u))*pos_x + std::tan((M_PI/2. -ang_u))*pos_u;
				else
					tmp_hit_xu_y =  std::tan((M_PI/2. + ang_u))*pos_x - std::tan((M_PI/2. +ang_u))*pos_u;

        for(int k = 0; k < nv; ++k)
        	{
          	if(flag_v.size()>0 && flag_v.find(k)!=flag_v.end())
          		continue;

          	TVectorD& vpos_v = hitv[k]->getRawHitCoords();
						double pos_v = vpos_v[0];
					  pos_v = pos_v / std::cos(ang_u);

					  double tmp_hit_xv_y = -999.;
					  if(ang_u>0)
					  	tmp_hit_xv_y =  std::tan((M_PI/2. - ang_u))*pos_x - std::tan((M_PI/2. - ang_u))*pos_v;
					  else
					  	tmp_hit_xv_y = -std::tan((M_PI/2. + ang_u))*pos_x + std::tan((M_PI/2. + ang_u))*pos_v;

					  if(std::fabs(hit_xu_y-hit_xv_y) > std::fabs(tmp_hit_xu_y-tmp_hit_xv_y))
						  {
		            hit_xu_y = tmp_hit_xu_y;
		            hit_xv_y = tmp_hit_xv_y;
		            used_j = j;
		            used_k = k;
						  }
	  			}
  		}

		if((std::fabs(hit_xu_y-999.)<1.)||(std::fabs(hit_xv_y+999.)<1.))
			continue;
		if(std::fabs(hit_xu_y-hit_xv_y) > cut_d[id_det])
			continue;

	  double buf_x = hit_xu_x + (hit_xu_y-hit_xv_y)/2.*std::tan(ang_u/2.);
	  double buf_y = (hit_xu_y+hit_xv_y)/2.;

	  TVector2 tmp_hitxy(buf_x, buf_y);
	  vect_SingleHitXYId[id_det].emplace_back(std::make_tuple(tmp_hitxy, id_track));

	  flag_x.insert(i);
	  flag_u.insert(used_j);
	  flag_v.insert(used_k);
	}

	return;
}


template<class Out>
std::vector<genfit::AbsMeasurement*> TCheckFiberXUV<Out>::Clusterize(const std::vector<std::unique_ptr<genfit::AbsMeasurement>>& hit)
{
	std::vector<genfit::AbsMeasurement*> hit_cluster = {};
	for(size_t i = 0; i < hit.size(); ++i)
		hit_cluster.emplace_back(hit[i]->clone());

	if(hit_cluster.size() < 2)
		return hit_cluster;

  std::set<int> flag_used;

  for(size_t i = 0; i < hit_cluster.size(); ++i)
	  {
			if(flag_used.size()>0 && flag_used.find(i)!=flag_used.end())
				continue;

			TVectorD& vpos = hit_cluster[i]->getRawHitCoords();
			double pos_cluster = vpos[0];
			double edge_left  = pos_cluster;
			double edge_right = pos_cluster;
			size_t size_cluster = 1;

      for(size_t j = i+1; j < hit_cluster.size(); ++j)
	      {
          if(flag_used.size()>0 && flag_used.find(j)!=flag_used.end())
          	continue;

          TVectorD& tmp_vpos = hit_cluster[j]->getRawHitCoords();

	        if(((pos_cluster-tmp_vpos[0])>0.) && std::fabs(edge_left - tmp_vpos[0])<0.0826)
	          {
	          	double tmp_pos_cluster = (pos_cluster * size_cluster + tmp_vpos[0]) / (size_cluster + 1);
	          	pos_cluster = tmp_pos_cluster;

	            flag_used.insert(j);
	            edge_left = tmp_vpos[0];
	            ++size_cluster;
	            j=i;
	          }
          else if(((tmp_vpos[0]-pos_cluster)>0.) && std::fabs(tmp_vpos[0]-edge_right)<0.0826)
	          {
	          	double tmp_pos_cluster = (pos_cluster * size_cluster + tmp_vpos[0]) / (size_cluster + 1);
	          	pos_cluster = tmp_pos_cluster;

	            flag_used.insert(j);
	            edge_right = tmp_vpos[0];
	            ++size_cluster;
	            j=i;
	          }
	      }

			vpos[0] = pos_cluster;
      hit_cluster[i]->setRawHitCoords(vpos);
	  }

  for(int i = (int)hit_cluster.size()-1; i >= 0; --i)
    {
	    if(flag_used.size()>0 && flag_used.find(i)!=flag_used.end())
		    {
		      delete hit_cluster[i];
		      hit_cluster.erase(hit_cluster.begin() + i);
		    }
  	}

	return hit_cluster;
}


template<class Out>
std::vector<genfit::AbsMeasurement*> TCheckFiberXUV<Out>::Clusterize(const std::vector<genfit::AbsMeasurement*>& hit)
{
	std::vector<genfit::AbsMeasurement*> hit_cluster = {};
	for(size_t i = 0; i < hit.size(); ++i)
		hit_cluster.emplace_back(hit[i]->clone());

	if(hit_cluster.size() < 2)
		return hit_cluster;

  std::set<int> flag_used;

  for(size_t i = 0; i < hit_cluster.size(); ++i)
	  {
			if(flag_used.size()>0 && flag_used.find(i)!=flag_used.end())
				continue;

			TVectorD& vpos = hit_cluster[i]->getRawHitCoords();
			double pos_cluster = vpos[0];
			double edge_left  = pos_cluster;
			double edge_right = pos_cluster;
			size_t size_cluster = 1;

      for(size_t j = i+1; j < hit_cluster.size(); ++j)
	      {
          if(flag_used.size()>0 && flag_used.find(j)!=flag_used.end())
          	continue;

          TVectorD& tmp_vpos = hit_cluster[j]->getRawHitCoords();

	        if(((pos_cluster-tmp_vpos[0])>0.) && std::fabs(edge_left - tmp_vpos[0])<0.0826)
	          {
	          	double tmp_pos_cluster = (pos_cluster * size_cluster + tmp_vpos[0]) / (size_cluster + 1);
	          	pos_cluster = tmp_pos_cluster;

	            flag_used.insert(j);
	            edge_left = tmp_vpos[0];
	            ++size_cluster;
	            j=i;
	          }
          else if(((tmp_vpos[0]-pos_cluster)>0.) && std::fabs(tmp_vpos[0]-edge_right)<0.0826)
	          {
	          	double tmp_pos_cluster = (pos_cluster * size_cluster + tmp_vpos[0]) / (size_cluster + 1);
	          	pos_cluster = tmp_pos_cluster;

	            flag_used.insert(j);
	            edge_right = tmp_vpos[0];
	            ++size_cluster;
	            j=i;
	          }
	      }

			vpos[0] = pos_cluster;
      hit_cluster[i]->setRawHitCoords(vpos);
	  }

  for(int i = (int)hit_cluster.size()-1; i >= 0; --i)
    {
	    if(flag_used.size()>0 && flag_used.find(i)!=flag_used.end())
		    {
		      delete hit_cluster[i];
		      hit_cluster.erase(hit_cluster.begin() + i);
		    }
  	}

	return hit_cluster;
}

template class TCheckFiberXUV<MCAnaEventG4Sol>;
template class TCheckFiberXUV<Ana_WasaEvent>;
