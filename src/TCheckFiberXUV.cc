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
#define SINGLE_FIBERXUV
#define DFUNCTION_FIBERXUV

using namespace std;
using namespace G4Sol;

template<class Out>
TCheckFiberXUV<Out>::TCheckFiberXUV(const THyphiAttributes& attribut) : TDataProcessInterface<Out>("CheckFiberXUV"), att(attribut)
{
  rand = new TRandom3();
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
	LocalHisto.h_EfficiencyFiberHit 	    = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EfficiencyFiberHit);
	LocalHisto.h_EfficiencySingleFiberHit = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EfficiencySingleFiberHit);

  for(size_t i = 2; i < 7; ++i)
	  {
	  	LocalHisto.h_ResidualFiberHitX[i]       					 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitX[i]);
	  	LocalHisto.h_ResidualFiberHitY[i]       					 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitY[i]);
	  	LocalHisto.h_ResidualFiberHitR[i]     		  			 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitR[i]);
	  	LocalHisto.h_ResidualFiberHitXY[i]		      			 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitXY[i]);
	  	LocalHisto.h_ResidualFiberHitX_Theta[i]			 			 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitX_Theta[i]);
	  	LocalHisto.h_ResidualFiberHitY_Theta[i] 					 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitY_Theta[i]);
	  	LocalHisto.h_ResidualFiberHitX_HitX[i] 				  	 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitX_HitX[i]);
	  	LocalHisto.h_ResidualFiberHitY_HitY[i]		 		  	 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitY_HitY[i]);
	  	LocalHisto.h_ResidualFiberHitR_Theta[i]						 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberHitR_Theta[i]);
	  	LocalHisto.h_ResidualSingleFiberHitX[i]			       = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualSingleFiberHitX[i]);
	  	LocalHisto.h_ResidualSingleFiberHitY[i]     		   = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualSingleFiberHitY[i]);
	  	LocalHisto.h_ResidualSingleFiberHitR[i]			       = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualSingleFiberHitR[i]);
	  	LocalHisto.h_ResidualSingleFiberHitX_Theta[i]			 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualSingleFiberHitX_Theta[i]);
	  	LocalHisto.h_ResidualSingleFiberHitY_Theta[i] 		 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualSingleFiberHitY_Theta[i]);
	  	LocalHisto.h_ResidualSingleFiberHitR_Theta[i] 		 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualSingleFiberHitR_Theta[i]);
	  	LocalHisto.h_EfficiencyFiberHit_Theta[i]    			 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EfficiencyFiberHit_Theta[i]);
	  	LocalHisto.h_EfficiencyFiberHit_dvalue[i]    			 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EfficiencyFiberHit_dvalue[i]);
	  	LocalHisto.h_EfficiencyFiberHit_mult[i]    	  		 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EfficiencyFiberHit_mult[i]);
	  	LocalHisto.h_EfficiencySingleFiberHit_Theta[i]     = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EfficiencySingleFiberHit_Theta[i]);
	  	LocalHisto.h_EfficiencySingleFiberHit_dvalue[i]    = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_EfficiencySingleFiberHit_dvalue[i]);
	  	LocalHisto.h_NumFiberHit_GoodReco[i]		      		 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_NumFiberHit_GoodReco[i]);
	  	LocalHisto.h_NumFiberHit_Ghost[i]      		  			 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_NumFiberHit_Ghost[i]);
	  	LocalHisto.h_FiberHit_dvalue[i]      				  	   = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_FiberHit_dvalue[i]);
	  	LocalHisto.h_FiberHitSingle_dvalue[i]		      	   = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_FiberHitSingle_dvalue[i]);
	  	LocalHisto.h_FiberHitReal_dvalue[i]      		  		 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_FiberHitReal_dvalue[i]);
	  	LocalHisto.h_FiberHitReal_dvalue_Theta[i]     		 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_FiberHitReal_dvalue_Theta[i]);
	  	LocalHisto.h_FiberHitReal_dvalue_Phi[i]  					 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_FiberHitReal_dvalue_Phi[i]);
	  	LocalHisto.h_FiberHitReal_dvalue_Theta03_Phi[i]    = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_FiberHitReal_dvalue_Theta03_Phi[i]);
	  	LocalHisto.h_FiberHitReal_dvalue_Theta310_Phi[i]   = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_FiberHitReal_dvalue_Theta310_Phi[i]);
	  	LocalHisto.h_FiberHitReal_dvalue_Theta1020_Phi[i]  = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_FiberHitReal_dvalue_Theta1020_Phi[i]);
	  	LocalHisto.h_FiberHitReal_dvalue_HitX[i]      		 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_FiberHitReal_dvalue_HitX[i]);
	  	LocalHisto.h_FiberHitReal_dvalue_HitY[i]      		 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_FiberHitReal_dvalue_HitY[i]);
	  	LocalHisto.h_FiberHitReal_dvalue_PosX[i]           = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_FiberHitReal_dvalue_PosX[i]);
	  	LocalHisto.h_FiberHitReal_dvalue_tanThetacosPhi[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_FiberHitReal_dvalue_tanThetacosPhi[i]);
	  	LocalHisto.h_FiberHitReal_dvalue_dfunction[i]      = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_FiberHitReal_dvalue_dfunction[i]);
	  	LocalHisto.h_FiberHit_Residualdvalue[i]            = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_FiberHit_Residualdvalue[i]);
	  	LocalHisto.h_FiberHit_Residualdvalue_Realdvalue[i] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_FiberHit_Residualdvalue_Realdvalue[i]);
	  }
}

template<class Out>
int TCheckFiberXUV<Out>::CheckHitXUV(const FullRecoEvent& RecoEvent)
{
	realIP.Set(rand->Uniform(RecoEvent.InteractionPoint[0] - randIPXY, RecoEvent.InteractionPoint[0] + randIPXY),
							rand->Uniform(RecoEvent.InteractionPoint[1] - randIPXY, RecoEvent.InteractionPoint[1] + randIPXY));

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
									//std::cout << "id_track: " << id_track << "\t ListHits x index j: " << j << "\n";
									//std::cout << "id_track: " << id_track << "\t vpos_x[0]: " << vpos_x[0] <<  "\n";
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
									//std::cout << "id_track: " << id_track << "\t ListHits u index j: " << j << "\n";
									//std::cout << "id_track: " << id_track << "\t vpos_u[0]: " << vpos_u[0] <<  "\n";
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
									//std::cout << "id_track: " << id_track << "\t ListHits v index j: " << j << "\n";
									//std::cout << "id_track: " << id_track << "\t vpos_v[0]: " << vpos_v[0] <<  "\n";
								}
						}

					if((counter_x == 0) || (counter_u == 0) || (counter_v == 0))
						continue;

					tmp_hit_x /= static_cast<double>(counter_x);
					tmp_hit_u /= static_cast<double>(counter_u);
					tmp_hit_v /= static_cast<double>(counter_v);

					//std::cout << "hitX: " << tmp_hit_x << "\t hitU: " << tmp_hit_u << "\t hitV: " << tmp_hit_v << "\n";
					//std::cout << "id_track: " << id_track <<"\t tmp_idhit_u.size(): " << tmp_idhit_u.size() << "\n";

		  		std::vector<SimHit> it_tmpFiber = it_track.second[i+imid];
			    //if(it_tmpFiber[0].momZ < 0.)
			   //  continue;

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
/* To select different MFT areas
		      if(meanPosX > 0.)
		      	continue;
		      if(meanPosY > std::fabs(std::tan(30. * M_PI / 180.) * meanPosX))
		      	continue;
*/
	  			TVector2 tmp_realHitXY(meanPosX, meanPosY);

			    double theta = TMath::ATan2( std::sqrt(it_tmpFiber[0].momX * it_tmpFiber[0].momX + it_tmpFiber[0].momY * it_tmpFiber[0].momY), 
			    															it_tmpFiber[0].momZ) * 180. / M_PI;

			    double phi = TMath::ATan2( it_tmpFiber[0].momY , it_tmpFiber[0].momX )  * 180. / M_PI;;

			    if(i_det == 4)
				    {
				    	std::swap(tmp_hit_u  ,   tmp_hit_v);
				    	std::swap(tmp_idhit_u, tmp_idhit_u);
				    }
			    else if(i_det == 5)
				    {
				    	std::swap(tmp_hit_x  ,   tmp_hit_v);
				    	std::swap(tmp_idhit_x, tmp_idhit_u);
				    }

					vect_realCombHit[i_det].emplace_back(std::make_tuple(tmp_hit_x, tmp_hit_u, tmp_hit_v));
					vect_realCombIdHit[i_det].insert(std::make_pair(id_track, std::make_tuple(tmp_idhit_x, tmp_idhit_u, tmp_idhit_v)));
	  			vect_realHitXYAngId[i_det].emplace_back(std::make_tuple(tmp_realHitXY, theta, id_track));

  				double tmp_dvalue = FindHitReal_dvalue(tmp_hit_x, tmp_hit_u, tmp_hit_v, i_det);

				  //std::cout << "Real_dvalue: " << std::fabs(tmp_dvalue) << "\n";
					LocalHisto.h_FiberHitReal_dvalue[i_det]                 ->Fill(tmp_dvalue, 1.);
				  LocalHisto.h_FiberHitReal_dvalue_Theta[i_det]           ->Fill(theta, tmp_dvalue, 1.);
				  LocalHisto.h_FiberHitReal_dvalue_Phi[i_det]             ->Fill(phi, tmp_dvalue, 1.);
				  LocalHisto.h_FiberHitReal_dvalue_HitX[i_det]            ->Fill(meanPosX, tmp_dvalue, 1.);
				  LocalHisto.h_FiberHitReal_dvalue_HitY[i_det]            ->Fill(meanPosY, tmp_dvalue, 1.);
				  LocalHisto.h_FiberHitReal_dvalue_PosX[i_det]            ->Fill(tmp_hit_x - realIP.X(), tmp_dvalue, 1.);
				  LocalHisto.h_FiberHitReal_dvalue_tanThetacosPhi[i_det]  ->Fill(std::tan(theta*M_PI/180.)*std::cos((phi+std::get<4>(param_d_funct1[i_det]))*M_PI/180.),
																																						tmp_dvalue, 1.);

				  double d_funct = d_function1(i_det, meanPosX, meanPosY);
				  //double d_funct = d_function2(i_det, tmp_hit_x);

				  //if(i_det == 3)
				  //	std::cout << "Real_d: " << tmp_dvalue << "\t Reco_d: " << d_funct << "\n";

				  LocalHisto.h_FiberHitReal_dvalue_dfunction[i_det]->Fill(d_funct, tmp_dvalue, 1.);
				  LocalHisto.h_FiberHit_Residualdvalue[i_det]->Fill(tmp_dvalue - d_funct, 1.);
					LocalHisto.h_FiberHit_Residualdvalue_Realdvalue[i_det]->Fill(tmp_dvalue, tmp_dvalue - d_funct, 1.);

				  if((0. <= theta ) && (theta < 3.))
				  	LocalHisto.h_FiberHitReal_dvalue_Theta03_Phi[i_det]->Fill(phi, tmp_dvalue, 1.);
				  else if((3. <= theta ) && (theta < 10.))
				  	LocalHisto.h_FiberHitReal_dvalue_Theta310_Phi[i_det]->Fill(phi, tmp_dvalue, 1.);
				  else if((10. <= theta ) && (theta < 20.))
				  	LocalHisto.h_FiberHitReal_dvalue_Theta1020_Phi[i_det]->Fill(phi, tmp_dvalue, 1.);
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

	  	FindHitXUV_v4(hitx, hitu, hitv, i_det);

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

		  		//std::cout << "(std::get<0>(it_CombHit.second)).size(): " << (std::get<0>(it_CombHit.second)).size() << "\n";
		  		//std::cout << "(std::get<1>(it_CombHit.second)).size(): " << (std::get<1>(it_CombHit.second)).size() << "\n";
		  		//std::cout << "(std::get<2>(it_CombHit.second)).size(): " << (std::get<2>(it_CombHit.second)).size() << "\n";

		  		//std::cout << "RecoEvent.ListHits[i].size(): "   << RecoEvent.ListHits[i].size() << "\n";
		  		//std::cout << "RecoEvent.ListHits[i+1].size(): " << RecoEvent.ListHits[i+1].size() << "\n";
		  		//std::cout << "RecoEvent.ListHits[i+2].size(): " << RecoEvent.ListHits[i+2].size() << "\n";

		  		for(size_t j = 0; j < (std::get<0>(it_CombHit.second)).size(); ++j)
		  		{
		  			//std::cout << "Inside x loop: j = " << j << "\n";
			  		//std::cout << "std::get<0>(it_CombHit.second)[j]: " << std::get<0>(it_CombHit.second)[j] << "\n";
						single_hitx.emplace_back(RecoEvent.ListHits[i][(std::get<0>(it_CombHit.second))[j]]->clone());
		  		}

		  		for(size_t j = 0; j < (std::get<1>(it_CombHit.second)).size(); ++j)
		  		{
		  			//std::cout << "Inside u loop: i_det: " <<  i_det << "\t id_track: " << id_track << "\t j = " << j << "\n";
			  		//std::cout << "std::get<1>(it_CombHit.second)[j]: " << std::get<1>(it_CombHit.second)[j] << "\n";
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

		  		FindSingleHitXUVId_v4(cluster_single_hitx, cluster_single_hitu, cluster_single_hitv, i_det, id_track);
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
	for(size_t i_det = 2; i_det < id_detector.size(); ++i_det)
		{
			//std::cout << "Fiber: " << i_det << "  NumRealHit: " << vect_realHitXYAngId[i_det].size() << "   NumRecoHit: " << vect_HitXY[i_det].size() << "\n";

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

					if((i_det == 2) && ifcut_MFTtheta_UFT3 &&
						((std::get<1>(vect_realHitXYAngId[i_det][j]) < minMFTtheta) || (std::get<1>(vect_realHitXYAngId[i_det][j]) > maxMFTtheta)))
							continue;

					TVector2 resHit = std::get<0>(vect_realHitXYAngId[i_det][j]) - vect_HitXY[i_det][closest_i];

		      LocalHisto.h_ResidualFiberHitX[i_det]      ->Fill(resHit.X(), 1.);
		      LocalHisto.h_ResidualFiberHitY[i_det]      ->Fill(resHit.Y(), 1.);
		      LocalHisto.h_ResidualFiberHitR[i_det]      ->Fill(resHit.Mod(), 1.);
		      LocalHisto.h_ResidualFiberHitXY[i_det]     ->Fill(resHit.X(), resHit.Y(), 1.);
		      LocalHisto.h_ResidualFiberHitX_Theta[i_det]->Fill(std::get<1>(vect_realHitXYAngId[i_det][j]), resHit.X(), 1.);
		      LocalHisto.h_ResidualFiberHitY_Theta[i_det]->Fill(std::get<1>(vect_realHitXYAngId[i_det][j]), resHit.Y(), 1.);
		      LocalHisto.h_ResidualFiberHitX_HitX[i_det] ->Fill(std::get<0>(vect_realHitXYAngId[i_det][j]).X(), resHit.X(), 1.);
		      LocalHisto.h_ResidualFiberHitY_HitY[i_det] ->Fill(std::get<0>(vect_realHitXYAngId[i_det][j]).Y(), resHit.Y(), 1.);
		      LocalHisto.h_ResidualFiberHitR_Theta[i_det]->Fill(std::get<1>(vect_realHitXYAngId[i_det][j]), resHit.Mod(), 1.);
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
/*
							if(resSingleHit.X() > fiber_width || resSingleHit.Y() > fiber_width)
								{
									std::cout << "RealSingX: " << (std::get<0>(vect_realHitXYAngId[i_det][j])).X() <<
																"\t RecoSingX: " << (std::get<0>(vect_SingleHitXYId[i_det][i])).X() << "\n";
									std::cout << "RealSingY: " << (std::get<0>(vect_realHitXYAngId[i_det][j])).Y() <<
																"\t RecoSingY: " << (std::get<0>(vect_SingleHitXYId[i_det][i])).Y() << "\n";
								}
*/
							LocalHisto.h_ResidualSingleFiberHitX[i_det]      ->Fill(resSingleHit.X(), 1.);
				      LocalHisto.h_ResidualSingleFiberHitY[i_det]      ->Fill(resSingleHit.Y(), 1.);
				      LocalHisto.h_ResidualSingleFiberHitR[i_det]      ->Fill(resSingleHit.Mod(), 1.);
				      LocalHisto.h_ResidualSingleFiberHitX_Theta[i_det]->Fill(std::get<1>(vect_realHitXYAngId[i_det][j]), resSingleHit.X(), 1.);
				      LocalHisto.h_ResidualSingleFiberHitY_Theta[i_det]->Fill(std::get<1>(vect_realHitXYAngId[i_det][j]), resSingleHit.Y(), 1.);
				      LocalHisto.h_ResidualSingleFiberHitR_Theta[i_det]->Fill(std::get<1>(vect_realHitXYAngId[i_det][j]), resSingleHit.Mod(), 1.);
						}
				}

			vect_SingleHitXYId[i_det].clear();
#endif

			vect_realCombIdHit[i_det].clear();

/*
			if(vect_realCombHit[i_det].size() == 2)
				{
					int i = id_detector[i_det];
				  std::vector<genfit::AbsMeasurement*> hitx = Clusterize(RecoEvent.ListHits[i]);
				  std::vector<genfit::AbsMeasurement*> hitu = Clusterize(RecoEvent.ListHits[i+1]);
				  std::vector<genfit::AbsMeasurement*> hitv = Clusterize(RecoEvent.ListHits[i+2]);

					std::cout << "i_det: " << i_det << "\n";

		      for(size_t j = 0; j < RecoEvent.ListHits[i].size(); ++j)
			      std::cout << "ListHitsX: " << (RecoEvent.ListHits[i][j]->getRawHitCoords())[0] << "\n";

		      for(size_t j = 0; j < hitx.size(); ++j)
			      std::cout << "ClusHitsX: " << (hitx[j]->getRawHitCoords())[0] << "\n";

		      for(size_t j = 0; j < RecoEvent.ListHits[i+1].size(); ++j)
			      std::cout << "ListHitsU: " << (RecoEvent.ListHits[i+1][j]->getRawHitCoords())[0] << "\n";

		      for(size_t j = 0; j < hitu.size(); ++j)
			      std::cout << "ClusHitsU: " << (hitu[j]->getRawHitCoords())[0] << "\n";

		      for(size_t j = 0; j < RecoEvent.ListHits[i+2].size(); ++j)
			      std::cout << "ListHitsV: " << (RecoEvent.ListHits[i+2][j]->getRawHitCoords())[0] << "\n";

		      for(size_t j = 0; j < hitv.size(); ++j)
			      std::cout << "ClusHitsV: " << (hitv[j]->getRawHitCoords())[0] << "\n";


					for(size_t i = 0; i < vect_realCombHit[i_det].size(); ++i)
						std::cout << "RealX: " << std::get<0>(vect_realCombHit[i_det][i]) << "\t RealU: " << std::get<1>(vect_realCombHit[i_det][i])
										<< "\t RealV: " << std::get<2>(vect_realCombHit[i_det][i]) << "\n";

					for(size_t i = 0; i < vect_CombHit[i_det].size(); ++i)
						std::cout << "RecoX: " << std::get<0>(vect_CombHit[i_det][i]) << "\t RecoU: " << std::get<1>(vect_CombHit[i_det][i])
										<< "\t RecoV: " << std::get<2>(vect_CombHit[i_det][i]) << "\n";
				}
*/

			int n_correct = 0;
		  std::set<int> flag_used;
		  std::set<int> flag_used_Single;

			for(size_t i = 0; i < vect_realCombHit[i_det].size(); ++i)
				{
					LocalHisto.h_EfficiencyFiberHit->Fill(i_det, "All", 1.);
					LocalHisto.h_EfficiencyFiberHit_Theta[i_det]->Fill(std::get<1>(vect_realHitXYAngId[i_det][i]), "All", 1.);

  				double tmp_realdvalue = FindHitReal_dvalue(std::get<0>(vect_realCombHit[i_det][i]), std::get<1>(vect_realCombHit[i_det][i]),
  																											std::get<2>(vect_realCombHit[i_det][i]), i_det);

					LocalHisto.h_EfficiencyFiberHit_dvalue[i_det]->Fill(tmp_realdvalue, "All", 1.);
					LocalHisto.h_EfficiencyFiberHit_mult[i_det]->Fill(vect_realCombHit[i_det].size(), "All", 1.);

					for(size_t j = 0; j < vect_CombHit[i_det].size(); ++j)
						{
							if(flag_used.size()>0 && flag_used.find(j)!=flag_used.end())
								continue;

							double res_x = std::fabs(std::get<0>(vect_realCombHit[i_det][i]) - std::get<0>(vect_CombHit[i_det][j]));
							double res_u = std::fabs(std::get<1>(vect_realCombHit[i_det][i]) - std::get<1>(vect_CombHit[i_det][j]));
							double res_v = std::fabs(std::get<2>(vect_realCombHit[i_det][i]) - std::get<2>(vect_CombHit[i_det][j]));

							if((res_x < fiber_width) && (res_u < fiber_width) && (res_v < fiber_width))
								{
									LocalHisto.h_EfficiencyFiberHit->Fill(i_det, "Acc", 1.);
									LocalHisto.h_EfficiencyFiberHit_Theta[i_det]->Fill(std::get<1>(vect_realHitXYAngId[i_det][i]), "Acc", 1.);

	  							double tmp_recodvalue = FindHitReal_dvalue(std::get<0>(vect_CombHit[i_det][j]), std::get<1>(vect_CombHit[i_det][j]),
  																														std::get<2>(vect_CombHit[i_det][j]), i_det);

									LocalHisto.h_EfficiencyFiberHit_dvalue[i_det]->Fill(tmp_recodvalue, "Acc", 1.);
									LocalHisto.h_EfficiencyFiberHit_mult[i_det]->Fill(vect_realCombHit[i_det].size(), "Acc", 1.);

									++n_correct;
									flag_used.insert(j);
									break;
								}
						}

					LocalHisto.h_EfficiencySingleFiberHit->Fill(i_det, "All", 1.);
					LocalHisto.h_EfficiencySingleFiberHit_Theta[i_det]->Fill(std::get<1>(vect_realHitXYAngId[i_det][i]), "All", 1.);
					LocalHisto.h_EfficiencySingleFiberHit_dvalue[i_det]->Fill(tmp_realdvalue, "All", 1.);

					for(size_t j = 0; j < vect_SingleCombHit[i_det].size(); ++j)
						{
							if(flag_used_Single.size()>0 && flag_used_Single.find(j)!=flag_used_Single.end())
								continue;

							double res_x = std::fabs(std::get<0>(vect_realCombHit[i_det][i]) - std::get<0>(vect_SingleCombHit[i_det][j]));
							double res_u = std::fabs(std::get<1>(vect_realCombHit[i_det][i]) - std::get<1>(vect_SingleCombHit[i_det][j]));
							double res_v = std::fabs(std::get<2>(vect_realCombHit[i_det][i]) - std::get<2>(vect_SingleCombHit[i_det][j]));

							if((res_x < fiber_width) && (res_u < fiber_width) && (res_v < fiber_width))
								{
									LocalHisto.h_EfficiencySingleFiberHit->Fill(i_det, "Acc", 1.);
									LocalHisto.h_EfficiencySingleFiberHit_Theta[i_det]->Fill(std::get<1>(vect_realHitXYAngId[i_det][i]), "Acc", 1.);

	  							double tmp_recodvalue = FindHitReal_dvalue(std::get<0>(vect_SingleCombHit[i_det][j]), std::get<1>(vect_SingleCombHit[i_det][j]),
  																															std::get<2>(vect_SingleCombHit[i_det][j]), i_det);

									LocalHisto.h_EfficiencySingleFiberHit_dvalue[i_det]->Fill(tmp_recodvalue, "Acc", 1.);

									flag_used_Single.insert(j);
									break;
								}
						}

				}

			//std::cout << "i_det: " << i_det << "\t vect_realCombHit[i_det].size(): " << vect_realCombHit[i_det].size() << "\t n_correct: " << n_correct << "\n";
			//std::cout << "i_det: " << i_det << "\t Not_correct: " << vect_realCombHit[i_det].size() - n_correct << "\n";
			int nghosts = vect_CombHit[i_det].size() - n_correct;

			LocalHisto.h_NumFiberHit_GoodReco[i_det]->Fill(vect_realCombHit[i_det].size(), n_correct, 1.);
			LocalHisto.h_NumFiberHit_Ghost[i_det]   ->Fill(vect_realCombHit[i_det].size(), nghosts, 1.);

			vect_realHitXYAngId[i_det].clear();
			vect_realCombHit[i_det].clear();
			vect_CombHit[i_det].clear();
			vect_SingleCombHit[i_det].clear();
		}

#ifdef DEBUG_FIBERXUV
	std::cout << "TCheckFiberXUV End\n";
#endif

  return 0;
}


template<class Out>
void TCheckFiberXUV<Out>::FindHitXUV_v1(const std::vector<genfit::AbsMeasurement*>& hitx,
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

		  double hit_u = -999.;
		  double hit_v = -999.;

		  int used_j = -999;
		  int used_k = -999;

		  for(int j = 0; j < nu ; ++j)
			  {
	        if(flag_u.size()>0 && flag_u.find(j)!=flag_u.end())
	        	continue;
					
					TVectorD& vpos_u = hitu[j]->getRawHitCoords();
					double pos_u = vpos_u[0];
					double tmp_hit_u = pos_u;
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
							double tmp_hit_v = pos_v;
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
			            hit_u = tmp_hit_u;
			            hit_v = tmp_hit_v;
			            used_j = j;
			            used_k = k;
							  }
		  			}
	  		}

			if((std::fabs(hit_xu_y-999.)<1.)||(std::fabs(hit_xv_y+999.)<1.))
				continue;

			if(std::fabs(hit_xu_y-hit_xv_y) > cut_d[id_det])
				continue;

			LocalHisto.h_FiberHit_dvalue[id_det]->Fill(hit_xu_y-hit_xv_y, 1.);

		  double buf_x = hit_xu_x + (hit_xu_y-hit_xv_y)/2.*std::tan(ang_u/2.);
		  double buf_y = (hit_xu_y+hit_xv_y)/2.;

		  TVector2 tmp_hitxy(buf_x, buf_y);
		  vect_HitXY[id_det].emplace_back(tmp_hitxy);

		  vect_CombHit[id_det].emplace_back(std::make_tuple(pos_x, hit_u, hit_v));

		  flag_x.insert(i);
		  flag_u.insert(used_j);
		  flag_v.insert(used_k);
		}

	return;
}


template<class Out>
double TCheckFiberXUV<Out>::FindHitReal_dvalue(double hitx, double hitu, double hitv, int id_det)
{
	double ang_u = ang[id_det][id_mid[id_det]] * M_PI / 180.;

	double pos_x = hitx;
  double hit_xu_y = 999.;
  double hit_xv_y = -999.;

	double pos_u = hitu;
	pos_u = pos_u / std::cos(ang_u);

	if(ang_u>0.)
		hit_xu_y = -std::tan((M_PI/2. - ang_u))*pos_x + std::tan((M_PI/2. -ang_u))*pos_u;
	else
		hit_xu_y =  std::tan((M_PI/2. + ang_u))*pos_x - std::tan((M_PI/2. +ang_u))*pos_u;

	double pos_v = hitv;
  pos_v = pos_v / std::cos(ang_u);

  if(ang_u>0)
  	hit_xv_y =  std::tan((M_PI/2. - ang_u))*pos_x - std::tan((M_PI/2. - ang_u))*pos_v;
  else
  	hit_xv_y = -std::tan((M_PI/2. + ang_u))*pos_x + std::tan((M_PI/2. + ang_u))*pos_v;

	return hit_xu_y-hit_xv_y;
}



template<class Out>
void TCheckFiberXUV<Out>::FindHitXUV_v2(const std::vector<genfit::AbsMeasurement*>& hitx,
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

	while(true)
		{
			double hit_xu_x = -999.;
		  double hit_xu_y =  999.;
		  double hit_xv_y = -999.;

		  double hit_u = -999.;
		  double hit_v = -999.;

		  int used_i = -999;
		  int used_j = -999;
		  int used_k = -999;

			for(int i = 0; i < nx ; ++i)
				{
					if(flag_x.size()>0 && flag_x.find(i)!=flag_x.end())
						continue;

					TVectorD& vpos_x = hitx[i]->getRawHitCoords();
					double pos_x = vpos_x[0];

				  for(int j = 0; j < nu ; ++j)
					  {
			        if(flag_u.size()>0 && flag_u.find(j)!=flag_u.end())
			        	continue;
							
							TVectorD& vpos_u = hitu[j]->getRawHitCoords();
							double pos_u = vpos_u[0];
							double tmp_hit_u = pos_u;
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
									double tmp_hit_v = pos_v;
								  pos_v = pos_v / std::cos(ang_u);

								  double tmp_hit_xv_y = -999.;
								  if(ang_u>0)
								  	tmp_hit_xv_y =  std::tan((M_PI/2. - ang_u))*pos_x - std::tan((M_PI/2. - ang_u))*pos_v;
								  else
								  	tmp_hit_xv_y = -std::tan((M_PI/2. + ang_u))*pos_x + std::tan((M_PI/2. + ang_u))*pos_v;

								  if(std::fabs(hit_xu_y-hit_xv_y) > std::fabs(tmp_hit_xu_y-tmp_hit_xv_y))
									  {
									  	hit_xu_x = pos_x;	
					            hit_xu_y = tmp_hit_xu_y;
					            hit_xv_y = tmp_hit_xv_y;
					            hit_u = tmp_hit_u;
					            hit_v = tmp_hit_v;
					            used_i = i;
					            used_j = j;
					            used_k = k;
									  }
				  			}
			  		}
			  }

				if( (std::fabs(hit_xu_x+999.) < 1.) || (std::fabs(hit_xu_y-999.) < 1.) || (std::fabs(hit_xv_y+999.) < 1.) )
					break;

				if(std::fabs(hit_xu_y-hit_xv_y) > cut_d[id_det])
					break;

				LocalHisto.h_FiberHit_dvalue[id_det]->Fill(hit_xu_y-hit_xv_y, 1.);

			  double buf_x = hit_xu_x + (hit_xu_y-hit_xv_y)/2.*std::tan(ang_u/2.);
			  double buf_y = (hit_xu_y+hit_xv_y)/2.;

			  TVector2 tmp_hitxy(buf_x, buf_y);
			  vect_HitXY[id_det].emplace_back(tmp_hitxy);

			  vect_CombHit[id_det].emplace_back(std::make_tuple(hit_xu_x, hit_u, hit_v));

			  flag_x.insert(used_i);
			  flag_u.insert(used_j);
			  flag_v.insert(used_k);	
		}

	return;
}


template<class Out>
void TCheckFiberXUV<Out>::FindHitXUV_v3(const std::vector<genfit::AbsMeasurement*>& hitx,
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

	while(true)
		{
			double hit_xu_x = -999.;
		  double hit_xu_y =  999.;
		  double hit_xv_y = -999.;

		  double diff_d = 999.;

		  double hit_u = -999.;
		  double hit_v = -999.;

		  int used_i = -999;
		  int used_j = -999;
		  int used_k = -999;

			for(int i = 0; i < nx ; ++i)
				{
					if(flag_x.size()>0 && flag_x.find(i)!=flag_x.end())
						continue;

					TVectorD& vpos_x = hitx[i]->getRawHitCoords();
					double pos_x = vpos_x[0];

				  for(int j = 0; j < nu ; ++j)
					  {
			        if(flag_u.size()>0 && flag_u.find(j)!=flag_u.end())
			        	continue;
							
							TVectorD& vpos_u = hitu[j]->getRawHitCoords();
							double pos_u = vpos_u[0];
							double tmp_hit_u = pos_u;
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
									double tmp_hit_v = pos_v;
								  pos_v = pos_v / std::cos(ang_u);

								  double tmp_hit_xv_y = -999.;
								  if(ang_u>0)
								  	tmp_hit_xv_y =  std::tan((M_PI/2. - ang_u))*pos_x - std::tan((M_PI/2. - ang_u))*pos_v;
								  else
								  	tmp_hit_xv_y = -std::tan((M_PI/2. + ang_u))*pos_x + std::tan((M_PI/2. + ang_u))*pos_v;

								  double tmp_x = pos_x + (tmp_hit_xu_y-tmp_hit_xv_y)/2.*std::tan(ang_u/2.);
								  double tmp_y = (tmp_hit_xu_y+tmp_hit_xv_y)/2.;

								  double tmp_d = tmp_hit_xu_y-tmp_hit_xv_y;
								  double expected_d = d_function1(id_det, tmp_x, tmp_y);
								  //double expected_d = d_function2(id_det, pos_x);
								  double tmp_diff_d = std::fabs(tmp_d - expected_d);

								  if(tmp_diff_d < diff_d)
									  {
									  	diff_d = tmp_diff_d;
									  	hit_xu_x = pos_x;	
					            hit_xu_y = tmp_hit_xu_y;
					            hit_xv_y = tmp_hit_xv_y;
					            hit_u = tmp_hit_u;
					            hit_v = tmp_hit_v;
					            used_i = i;
					            used_j = j;
					            used_k = k;
									  }
				  			}
			  		}
			  }

				if( (std::fabs(hit_xu_x+999.) < 1.) || (std::fabs(hit_xu_y-999.) < 1.) || (std::fabs(hit_xv_y+999.) < 1.) )
					break;

				if(diff_d > cut_diff_d[id_det])
					break;

				LocalHisto.h_FiberHit_dvalue[id_det]->Fill(hit_xu_y-hit_xv_y, 1.);

			  double buf_x = hit_xu_x + (hit_xu_y-hit_xv_y)/2.*std::tan(ang_u/2.);
			  double buf_y = (hit_xu_y+hit_xv_y)/2.;

			  TVector2 tmp_hitxy(buf_x, buf_y);
			  vect_HitXY[id_det].emplace_back(tmp_hitxy);

			  vect_CombHit[id_det].emplace_back(std::make_tuple(hit_xu_x, hit_u, hit_v));

			  flag_x.insert(used_i);
			  flag_u.insert(used_j);
			  flag_v.insert(used_k);	
		}

	return;
}


template<class Out>
void TCheckFiberXUV<Out>::FindHitXUV_v4(const std::vector<genfit::AbsMeasurement*>& hitx,
																					const std::vector<genfit::AbsMeasurement*>& hitu,
																					const std::vector<genfit::AbsMeasurement*>& hitv,
																						int id_det)
{
  int nx = hitx.size();
  int nu = hitu.size();
  int nv = hitv.size();
	if((nx==0)||(nu==0)||(nv==0))
		return;

	double ang_u = ang[id_det][id_mid[id_det]] * M_PI / 180.;

	for(int i = 0; i < nx ; ++i)
		{
			TVectorD& vpos_x = hitx[i]->getRawHitCoords();
			double pos_x = vpos_x[0];

		  for(int j = 0; j < nu ; ++j)
			  {
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
	          	TVectorD& vpos_v = hitv[k]->getRawHitCoords();
							double pos_v = vpos_v[0];
						  pos_v = pos_v / std::cos(ang_u);

						  double tmp_hit_xv_y = -999.;
						  if(ang_u>0)
						  	tmp_hit_xv_y =  std::tan((M_PI/2. - ang_u))*pos_x - std::tan((M_PI/2. - ang_u))*pos_v;
						  else
						  	tmp_hit_xv_y = -std::tan((M_PI/2. + ang_u))*pos_x + std::tan((M_PI/2. + ang_u))*pos_v;

						  double tmp_x = pos_x + (tmp_hit_xu_y-tmp_hit_xv_y)/2.*std::tan(ang_u/2.);
						  double tmp_y = (tmp_hit_xu_y+tmp_hit_xv_y)/2.;

						  double tmp_d = tmp_hit_xu_y-tmp_hit_xv_y;
						  double expected_d = d_function1(id_det, tmp_x, tmp_y);
						  //double expected_d = d_function2(id_det, pos_x);
						  double tmp_diff_d = std::fabs(tmp_d - expected_d);

						  if(tmp_diff_d < cut_diff_d[id_det])
							  {
									LocalHisto.h_FiberHit_dvalue[id_det]->Fill(tmp_hit_xu_y-tmp_hit_xv_y, 1.);

								  TVector2 tmp_hitxy(tmp_x, tmp_y);
								  vect_HitXY[id_det].emplace_back(tmp_hitxy);

								  vect_CombHit[id_det].emplace_back(std::make_tuple(vpos_x[0], vpos_u[0], vpos_v[0]));
							  }
		  			}
	  		}
	  }

	return;
}



template<class Out>
void TCheckFiberXUV<Out>::FindSingleHitXUVId_v1(const std::vector<genfit::AbsMeasurement*>& hitx,
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

		double hit_u = -999.;
		double hit_v = -999.;

	  int used_j = -999;
	  int used_k = -999;

	  for(int j = 0; j < nu ; ++j)
		  {
        if(flag_u.size()>0 && flag_u.find(j)!=flag_u.end())
        	continue;
				
				TVectorD& vpos_u = hitu[j]->getRawHitCoords();
				double pos_u = vpos_u[0];
				double tmp_hit_u = pos_u;
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
						double tmp_hit_v = pos_v;
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
		            hit_u = tmp_hit_u;
		            hit_v = tmp_hit_v;
		            used_j = j;
		            used_k = k;
						  }
	  			}
  		}

		if((std::fabs(hit_xu_y-999.)<1.)||(std::fabs(hit_xv_y+999.)<1.))
			continue;

		if(std::fabs(hit_xu_y-hit_xv_y) > cut_d[id_det])
			continue;

		LocalHisto.h_FiberHitSingle_dvalue[id_det]->Fill(hit_xu_y-hit_xv_y, 1.);

	  double buf_x = hit_xu_x + (hit_xu_y-hit_xv_y)/2.*std::tan(ang_u/2.);
	  double buf_y = (hit_xu_y+hit_xv_y)/2.;

	  TVector2 tmp_hitxy(buf_x, buf_y);
	  vect_SingleHitXYId[id_det].emplace_back(std::make_tuple(tmp_hitxy, id_track));

		vect_SingleCombHit[id_det].emplace_back(std::make_tuple(hit_xu_x, hit_u, hit_v));

	  flag_x.insert(i);
	  flag_u.insert(used_j);
	  flag_v.insert(used_k);
	}

	return;
}


template<class Out>
void TCheckFiberXUV<Out>::FindSingleHitXUVId_v2(const std::vector<genfit::AbsMeasurement*>& hitx,
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

	double ang_u = ang[id_det][id_mid[id_det]] * M_PI / 180.;

	while(true)
		{
			double hit_xu_x = -999.;
		  double hit_xu_y =  999.;
		  double hit_xv_y = -999.;

			double hit_u = -999.;
			double hit_v = -999.;

		  int used_i = -999;
		  int used_j = -999;
		  int used_k = -999;

			for(int i = 0; i < nx ; ++i)
				{
					if(flag_x.size()>0 && flag_x.find(i)!=flag_x.end())
						continue;

					TVectorD& vpos_x = hitx[i]->getRawHitCoords();
					double pos_x = vpos_x[0];

				  for(int j = 0; j < nu ; ++j)
					  {
			        if(flag_u.size()>0 && flag_u.find(j)!=flag_u.end())
			        	continue;
							
							TVectorD& vpos_u = hitu[j]->getRawHitCoords();
							double pos_u = vpos_u[0];
							double tmp_hit_u = pos_u;
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
									double tmp_hit_v = pos_v;
								  pos_v = pos_v / std::cos(ang_u);

								  double tmp_hit_xv_y = -999.;
								  if(ang_u>0)
								  	tmp_hit_xv_y =  std::tan((M_PI/2. - ang_u))*pos_x - std::tan((M_PI/2. - ang_u))*pos_v;
								  else
								  	tmp_hit_xv_y = -std::tan((M_PI/2. + ang_u))*pos_x + std::tan((M_PI/2. + ang_u))*pos_v;

								  if(std::fabs(hit_xu_y-hit_xv_y) > std::fabs(tmp_hit_xu_y-tmp_hit_xv_y))
									  {
									  	hit_xu_x = pos_x;	
					            hit_xu_y = tmp_hit_xu_y;
					            hit_xv_y = tmp_hit_xv_y;
					            hit_u = tmp_hit_u;
					            hit_v = tmp_hit_v;
					            used_i = i;
					            used_j = j;
					            used_k = k;
									  }
				  			}
			  		}
			  }

				if( (std::fabs(hit_xu_x+999.) < 1.) || (std::fabs(hit_xu_y-999.) < 1.) || (std::fabs(hit_xv_y+999.) < 1.) )
					break;

				if(std::fabs(hit_xu_y-hit_xv_y) > cut_d[id_det])
					break;

				LocalHisto.h_FiberHitSingle_dvalue[id_det]->Fill(hit_xu_y-hit_xv_y, 1.);

			  double buf_x = hit_xu_x + (hit_xu_y-hit_xv_y)/2.*std::tan(ang_u/2.);
			  double buf_y = (hit_xu_y+hit_xv_y)/2.;

			  TVector2 tmp_hitxy(buf_x, buf_y);
	  		vect_SingleHitXYId[id_det].emplace_back(std::make_tuple(tmp_hitxy, id_track));

	  		vect_SingleCombHit[id_det].emplace_back(std::make_tuple(hit_xu_x, hit_u, hit_v));

			  flag_x.insert(used_i);
			  flag_u.insert(used_j);
			  flag_v.insert(used_k);	
		}

	return;
}


template<class Out>
void TCheckFiberXUV<Out>::FindSingleHitXUVId_v3(const std::vector<genfit::AbsMeasurement*>& hitx,
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

	double ang_u = ang[id_det][id_mid[id_det]] * M_PI / 180.;

	while(true)
		{
			double hit_xu_x = -999.;
		  double hit_xu_y =  999.;
		  double hit_xv_y = -999.;

		  double diff_d = 999.;

			double hit_u = -999.;
			double hit_v = -999.;

		  int used_i = -999;
		  int used_j = -999;
		  int used_k = -999;

			for(int i = 0; i < nx ; ++i)
				{
					if(flag_x.size()>0 && flag_x.find(i)!=flag_x.end())
						continue;

					TVectorD& vpos_x = hitx[i]->getRawHitCoords();
					double pos_x = vpos_x[0];

				  for(int j = 0; j < nu ; ++j)
					  {
			        if(flag_u.size()>0 && flag_u.find(j)!=flag_u.end())
			        	continue;
							
							TVectorD& vpos_u = hitu[j]->getRawHitCoords();
							double pos_u = vpos_u[0];
							double tmp_hit_u = pos_u;
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
									double tmp_hit_v = pos_v;
								  pos_v = pos_v / std::cos(ang_u);

								  double tmp_hit_xv_y = -999.;
								  if(ang_u>0)
								  	tmp_hit_xv_y =  std::tan((M_PI/2. - ang_u))*pos_x - std::tan((M_PI/2. - ang_u))*pos_v;
								  else
								  	tmp_hit_xv_y = -std::tan((M_PI/2. + ang_u))*pos_x + std::tan((M_PI/2. + ang_u))*pos_v;

								  double tmp_x = pos_x + (tmp_hit_xu_y-tmp_hit_xv_y)/2.*std::tan(ang_u/2.);
								  double tmp_y = (tmp_hit_xu_y+tmp_hit_xv_y)/2.;

								  double tmp_d = tmp_hit_xu_y-tmp_hit_xv_y;
								  double expected_d = d_function1(id_det, tmp_x, tmp_y);
								  //double expected_d = d_function2(id_det, pos_x);
								  double tmp_diff_d = std::fabs(tmp_d - expected_d);


								  if(tmp_diff_d < diff_d)
									  {
									  	diff_d = tmp_diff_d;
									  	hit_xu_x = pos_x;	
					            hit_xu_y = tmp_hit_xu_y;
					            hit_xv_y = tmp_hit_xv_y;
					            hit_u = tmp_hit_u;
					            hit_v = tmp_hit_v;
					            used_i = i;
					            used_j = j;
					            used_k = k;
									  }
				  			}
			  		}
			  }

				if( (std::fabs(hit_xu_x+999.) < 1.) || (std::fabs(hit_xu_y-999.) < 1.) || (std::fabs(hit_xv_y+999.) < 1.) )
					break;

				if(diff_d > cut_diff_d[id_det])
					break;

				LocalHisto.h_FiberHitSingle_dvalue[id_det]->Fill(hit_xu_y-hit_xv_y, 1.);

			  double buf_x = hit_xu_x + (hit_xu_y-hit_xv_y)/2.*std::tan(ang_u/2.);
			  double buf_y = (hit_xu_y+hit_xv_y)/2.;

			  TVector2 tmp_hitxy(buf_x, buf_y);
	  		vect_SingleHitXYId[id_det].emplace_back(std::make_tuple(tmp_hitxy, id_track));

	  		vect_SingleCombHit[id_det].emplace_back(std::make_tuple(hit_xu_x, hit_u, hit_v));

			  flag_x.insert(used_i);
			  flag_u.insert(used_j);
			  flag_v.insert(used_k);	
		}

	return;
}


template<class Out>
void TCheckFiberXUV<Out>::FindSingleHitXUVId_v4(const std::vector<genfit::AbsMeasurement*>& hitx,
																									const std::vector<genfit::AbsMeasurement*>& hitu,
																									const std::vector<genfit::AbsMeasurement*>& hitv,
																										int id_det, int id_track)
{
  int nx = hitx.size();
  int nu = hitu.size();
  int nv = hitv.size();
	if((nx==0)||(nu==0)||(nv==0))
		return;

	double ang_u = ang[id_det][id_mid[id_det]] * M_PI / 180.;

	for(int i = 0; i < nx ; ++i)
		{
			TVectorD& vpos_x = hitx[i]->getRawHitCoords();
			double pos_x = vpos_x[0];

		  for(int j = 0; j < nu ; ++j)
			  {
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
	          	TVectorD& vpos_v = hitv[k]->getRawHitCoords();
							double pos_v = vpos_v[0];
						  pos_v = pos_v / std::cos(ang_u);

						  double tmp_hit_xv_y = -999.;
						  if(ang_u>0)
						  	tmp_hit_xv_y =  std::tan((M_PI/2. - ang_u))*pos_x - std::tan((M_PI/2. - ang_u))*pos_v;
						  else
						  	tmp_hit_xv_y = -std::tan((M_PI/2. + ang_u))*pos_x + std::tan((M_PI/2. + ang_u))*pos_v;

						  double tmp_x = pos_x + (tmp_hit_xu_y-tmp_hit_xv_y)/2.*std::tan(ang_u/2.);
						  double tmp_y = (tmp_hit_xu_y+tmp_hit_xv_y)/2.;

						  double tmp_d = tmp_hit_xu_y-tmp_hit_xv_y;
						  double expected_d = d_function1(id_det, tmp_x, tmp_y);
						  //double expected_d = d_function2(id_det, pos_x);
						  double tmp_diff_d = std::fabs(tmp_d - expected_d);

						  if(tmp_diff_d < cut_diff_d[id_det])
							  {
									LocalHisto.h_FiberHitSingle_dvalue[id_det]->Fill(tmp_hit_xu_y-tmp_hit_xv_y, 1.);

								  TVector2 tmp_hitxy(tmp_x, tmp_y);
								  vect_SingleHitXYId[id_det].emplace_back(tmp_hitxy, id_track);

								  vect_SingleCombHit[id_det].emplace_back(std::make_tuple(vpos_x[0], vpos_u[0], vpos_v[0]));
							  }
		  			}
	  		}
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
	    //std::cout << "PrePos: " << vpos[0] << "\n";

      for(size_t j = i+1; j < hit_cluster.size(); ++j)
	      {
          if(flag_used.size()>0 && flag_used.find(j)!=flag_used.end())
          	continue;

          TVectorD& tmp_vpos = hit_cluster[j]->getRawHitCoords();

	        if(((pos_cluster-tmp_vpos[0])>0.) && std::fabs(edge_left-tmp_vpos[0])<0.0826)
	          {
	          	pos_cluster = (pos_cluster * size_cluster + tmp_vpos[0]) / static_cast<double>(size_cluster + 1);
	            flag_used.insert(j);
	            edge_left = tmp_vpos[0];
	            ++size_cluster;
	            j=i;
	            //std::cout << "PrePos: " << tmp_vpos[0] << "\n";
	          }
          else if(((tmp_vpos[0]-pos_cluster)>0.) && std::fabs(tmp_vpos[0]-edge_right)<0.0826)
	          {
	          	pos_cluster = (pos_cluster * size_cluster + tmp_vpos[0]) / static_cast<double>(size_cluster + 1);
	            flag_used.insert(j);
	            edge_right = tmp_vpos[0];
	            ++size_cluster;
	            j=i;
	            //std::cout << "PrePos: " << tmp_vpos[0] << "\n";
	          }
	      }

			vpos[0] = pos_cluster;
	    //std::cout << "PostPos: " << pos_cluster << "\t vpos[0]: " << vpos[0] << "\n";
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

  //for(size_t i = 0; i < hit_cluster.size(); ++i)
  	//std::cout << "Vect_PostPos: " << (hit_cluster[i]->getRawHitCoords())[0] << "\n";
  //std::cout << "Presize: " << hit.size() << "\t Postsize: " << hit_cluster.size() << "\n";
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


template<class Out>
double TCheckFiberXUV<Out>::d_function1(int id_det, double hitx, double hity)
{
	double x = hitx - realIP.X();
	double y = hity - realIP.Y();
	double z = Zpos_Fiber[id_det] - Zpos_target;

	double theta = TMath::ATan2( std::sqrt(x * x + y * y), z);
	double phi = TMath::ATan2( y , x );

	if(theta*180./M_PI < 5.)
		{
			double d = d_function2(id_det, hitx);
			return d;
		}

	double f = std::tan(theta) * std::cos(phi + std::get<4>(param_d_funct1[id_det]) * M_PI / 180.);

	double d = std::get<0>(param_d_funct1[id_det]) + std::get<1>(param_d_funct1[id_det]) * f + std::get<2>(param_d_funct1[id_det]) * f * f
									+ std::get<3>(param_d_funct1[id_det]) * f * f * f; //d_function_1
	return d;
}

template<class Out>
double TCheckFiberXUV<Out>::d_function2(int id_det, double posx)
{
	double x = posx - realIP.X();

	double d = std::get<0>(param_d_funct2[id_det]) + std::get<1>(param_d_funct2[id_det]) * x + std::get<2>(param_d_funct2[id_det]) * x * x
								+ std::get<3>(param_d_funct2[id_det]) * x * x * x; //d_function_2
	return d;
}

template class TCheckFiberXUV<MCAnaEventG4Sol>;
template class TCheckFiberXUV<Ana_WasaEvent>;
