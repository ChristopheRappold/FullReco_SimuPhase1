#include "TKalmanFilter_FRS.h"
#include "MathematicalTools.hh"

#include "GFTrack.h"
#include "RKTrackRep.h"
#include "GFKalman.h"
#include "GFException.h"
#include "GFHypRecoHit.h"


#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TGeoMatrix.h"

#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "TObjArray.h"

#include <algorithm>
#include <numeric>
//#include "Ana_Event/TTrackSimple.hh"

// #define DEBUG_KALMAN2
// #define DEBUG_KALMAN
//#define ROTATION_KALMAN

//#define OLD 1
//#define DEBUG_DAISUKE

//#define OLDHITS

using namespace std;


TKalmanFilter_FRS::TKalmanFilter_FRS(const THyphiAttributes& attribut):TDataProcessInterfaceMC("mom_fit_kalmanMC"),att(attribut)
{
  //int nb_iteration = 2;
  //Fitter.setLazy(0); // tell the fitter to skip hits if error occurs  
  //Fitter.setNumIterations(nb_iteration);

  //bool hermit_interpolation = true;
  //bool adding_noise_multiscattering = true;
  //bool B_interpolation = true;

  rep = new RKTrackRep();
  rep_length = new RKTrackRep();

  Vtracks = new GFTrack (rep);
  Vtracks->setSmoothing();
  //Vtracks->SetException(false);

  // if(att.EventGen==true)
  //   {
  //     TGeoNode* world = gGeoManager->GetNode(0);
  //     TGeoNode* setup = world->GetDaughter(0);
  //     TObjArray* list_node = setup->GetNodes();
  //     TString name_node_tr0x("TR0BoxLog_0");
  //     TGeoNodeMatrix* temp_node = (TGeoNodeMatrix*)list_node->FindObject(name_node_tr0x);
  //     TGeoMatrix* mat_vol = temp_node->GetMatrix();
  //     assert(mat_vol!=0);
  //     const Double_t* orig = mat_vol->GetTranslation();
  //     Plane_time.SetXYZ(orig[0],orig[1],orig[2]);
  //   }
  // else
  //   {
  //     TVector3 Start_Pos(0.,0.,att.TOFStart_Z);
  //     Rotating(Start_Pos,0.1);
  //     //Start_Pos *= 0.1;
  //     Plane_time=Start_Pos;
  //   }


  //std::vector<GFRectFinitePlane*> finiteplane;

  LayerDetectorType name_plane[10] = {TR0X,TR1X,DC1X,TR2X,DC2X,DC2XP,TOFP,DC2Y,DC2YP,TFWN};
  
  TGeoNode* world = gGeoManager->GetNode(0);
  TGeoNode* setup = world->GetDaughter(0);
  TObjArray* list_node = setup->GetNodes();
  
  TString name_node[10] = { "TR0_1","TR1_1","DC1_1","TR2_1","DC2_1","DC3_1","TOFP_1","DC2STOP_1","STOP_1","STOP2_1"};
  std::vector<TGeoMatrix*> mat_vol;
  for(int i=0;i<10;++i)  
    {
      TGeoNodeMatrix* temp_node = (TGeoNodeMatrix*)list_node->FindObject(name_node[i]);
      mat_vol.push_back(temp_node->GetMatrix());
      assert(mat_vol.back()!=0);
    }
  
  list_Plane.resize(SIZEOF_LAYERDETTYPE,0);
  for(int i=0;i<10;++i)
    {
      const Double_t* orig = mat_vol[i]->GetTranslation();
      const Double_t* rot = mat_vol[i]-> GetRotationMatrix();
      TVector3 O(orig[0],orig[1],orig[2]);
      TVector3 U(rot[0],rot[3],rot[6]);
      TVector3 V(rot[1],rot[4],rot[7]);
      list_Plane[name_plane[i]]= new GFDetPlane(O,U,V);
    }



  Kalman.setLazy(0); // tell the fitter to skip hits if error occurs  
  Kalman.setNumIterations(3);
  //Kalman.setNumIterations(2);

}

TKalmanFilter_FRS::~TKalmanFilter_FRS()
{
  if(Vtracks!=NULL)
    {
      delete Vtracks;
      Vtracks= 0;
    }
  if(rep_length!=NULL)
    {
      delete rep_length;
      rep_length=0;
    }
  for(unsigned int i=0;i<list_Plane.size();++i)
    if(list_Plane[i]!=0)
      {
	delete list_Plane[i];
	list_Plane[i] = 0;
      }
}

int TKalmanFilter_FRS::operator() (FullRecoEvent& RecoEvent,MCAnaEvent* OutTree)
{

  int result_mom = Exec(RecoEvent,OutTree);

  return SoftExit(result_mom);
}

int TKalmanFilter_FRS::Exec(FullRecoEvent& RecoEvent,MCAnaEvent* OutTree)
{


  //  int result_kalman = Kalman_Filter(RecoEvent);
  int result_kalman = Kalman_Filter_FromTrack(RecoEvent);


  //OutTree->Ntracks=0;
  
#ifdef DEBUG_DAISUKE
  std::cout<<"-------------------------------------------------------- "<<std::endl;
#endif
  for(map_mom3::const_iterator it_mom_pv = RecoEvent.Mom_Particle.begin(); it_mom_pv!=RecoEvent.Mom_Particle.end();++it_mom_pv)
    {
#ifdef DEBUG_DAISUKE
      std::cout<<"********************* "<<std::endl;
#endif

      //TTrackSimple track_temp;
      THyphiTrack track_temp;
      track_temp.Chi2 = it_mom_pv->second[7];
      track_temp.Chi2_Y = it_mom_pv->second[8];
      track_temp.Mass = it_mom_pv->first.mass;
      track_temp.MomMass.SetXYZM(it_mom_pv->second[0],it_mom_pv->second[1],it_mom_pv->second[2],it_mom_pv->first.mass);
      track_temp.Mom.SetXYZ(it_mom_pv->second[0],it_mom_pv->second[1],it_mom_pv->second[2]);
    
      track_temp.Charge = it_mom_pv->second[6];
      track_temp.BarId = it_mom_pv->second[25];
      track_temp.Beta =it_mom_pv->second[26];
      track_temp.HitTR1.SetXYZ(it_mom_pv->second[13],it_mom_pv->second[14],it_mom_pv->second[15]);
      track_temp.HitTR2.SetXYZ(it_mom_pv->second[16],it_mom_pv->second[17],it_mom_pv->second[18]);
      track_temp.HitDC2.SetXYZ(it_mom_pv->second[19],it_mom_pv->second[20],it_mom_pv->second[21]);
      track_temp.NumDC2.SetXYZ(it_mom_pv->second[22],it_mom_pv->second[23],it_mom_pv->second[24]);
      track_temp.Pval2 = it_mom_pv->second[27];
      track_temp.TofsBar = it_mom_pv->second[28];
      track_temp.RefPoint.SetXYZ(it_mom_pv->second[3],it_mom_pv->second[4],it_mom_pv->second[5]);
      track_temp.HitTOF.SetXYZ(it_mom_pv->second[29],it_mom_pv->second[30],it_mom_pv->second[31]);
      track_temp.dE=it_mom_pv->second[12];

      track_temp.HitDC1.SetXYZ(it_mom_pv->second[32],it_mom_pv->second[33],it_mom_pv->second[34]);
      track_temp.HitDC1p.SetXYZ(it_mom_pv->second[35],it_mom_pv->second[36],it_mom_pv->second[37]);
       
      track_temp.MomIni.SetXYZ(it_mom_pv->second[38],it_mom_pv->second[39],it_mom_pv->second[40]);
      track_temp.RChiIni = it_mom_pv->second[41];
      track_temp.BetaIni= it_mom_pv->second[42];
      track_temp.Chi2_X = it_mom_pv->second[9];//FOR TEST it_mom_pv->second[99];//it_mom_pv->second[43];
      track_temp.PathLength = it_mom_pv->second[10];
      track_temp.TOF = it_mom_pv->second[11];
      //        if(it_mom_pv->second[6]>0)
      //        std::cout<<it_mom_pv->second[6]<<" BarId="<<it_mom_pv->second[25]<<" "<<it_mom_pv->second[43]<<std::endl;


      //////////////////////////

      // cout<<" before Tree"<<endl<<"state: ";
      for(int row=0;row<6;row++)
	{
	  track_temp.State[row] = it_mom_pv->second[44+row];
	  // 	   cout<<it_mom_pv->second[44+row]<<" ";
	}
      // cout<<endl<<"cov:";
	
      for(int row=0;row<6;row++)
	for(int col=0;col<6;col++)
	  {
	    track_temp.Cov[row][col] = it_mom_pv->second[50+row*6+col];
	    //cout<<it_mom_pv->second[49+row*5+col]<<endl;
	  }
      //       cout<<"*****"<<endl;

      track_temp.MomTof.SetXYZ(it_mom_pv->second[87],it_mom_pv->second[88],it_mom_pv->second[89]);

      track_temp.Residu_X_TR1 = it_mom_pv->second[90];
      track_temp.Residu_X_TR2 = it_mom_pv->second[91];
      track_temp.Residu_X_DC2 = it_mom_pv->second[92];
      track_temp.Residu_X_TOF = it_mom_pv->second[93];
      track_temp.Residu_Y_TR1 = it_mom_pv->second[94];
      track_temp.Residu_Y_TR2 = it_mom_pv->second[95];
      track_temp.Residu_Y_DC2 = it_mom_pv->second[96];
      track_temp.Residu_Y_TOF = it_mom_pv->second[97];

      track_temp.MC_status = it_mom_pv->second[98];

      OutTree->Add_Track(track_temp);



#ifdef DEBUG_DAISUKE
      std::cout<<"HitTR1 ("<<it_mom_pv->second[13]<<","<<it_mom_pv->second[14]<<","<<it_mom_pv->second[15]<<")"<<std::endl;
      std::cout<<"HitTR2 ("<<it_mom_pv->second[16]<<","<<it_mom_pv->second[17]<<","<<it_mom_pv->second[18]<<")"<<std::endl;
      std::cout<<"HitDC2 ("<<it_mom_pv->second[19]<<","<<it_mom_pv->second[20]<<","<<it_mom_pv->second[21]<<")"<<std::endl;
      std::cout<<"NumDC2 ("<<it_mom_pv->second[22]<<","<<it_mom_pv->second[23]<<","<<it_mom_pv->second[24]<<")"<<std::endl;
      
      std::cout<<"Pval2 = "<<it_mom_pv->second[27]<<std::endl;
      std::cout<<"TofsBar = "<<it_mom_pv->second[28]<<std::endl;
      std::cout<<"BarId = "<<it_mom_pv->second[25]<<std::endl;
      std::cout<<"Charge = "<<it_mom_pv->second[6]<<std::endl;
      std::cout<<"Beta = "<<it_mom_pv->second[26]<<std::endl;
#endif
    }
#ifdef DEBUG_DAISUKE
  std::cout<<"-------------------------------------------------------- "<<std::endl;
#endif

  return result_kalman;
}

int TKalmanFilter_FRS::SoftExit(int result_full)
{
  return result_full;    
}


int TKalmanFilter_FRS::Kalman_Filter_FromTrack(FullRecoEvent& RecoEvent)
{

  double x_TOF;  
  double y_TOF;
  double z_TOF;
  double E_TOF;

  int TofsBar=-1;
  TVector3 mom_Track;  
  TVector3 mom_Track_ini;  
  double mass_Track_ini;
  double beta_Track_ini;
  double rchi_Track_ini;
  int charge_Track_ini;
  int bar_Track;
  TVector3 HitTR1;
  TVector3 HitTR2;
  TVector3 HitDC2;
  TVector3 HitDC1;
  TVector3 HitDC1p;
  TVector3 NumDC2;
  TVector3 HitTOF;
  int ProtonNewOrOld;
  int NeigHit;
  double CoinWindow;


  double best_pv[32];
  double best_mass[32];
  double best_beta[32];
  double best_plen[32];
  double best_tof[32];
  double best_mom[32];
  double best_dE[32];
  int best_charge[32];
  bool best_tofs_coin[32];
  int best_layer[32];
      
  double best_rchi[32];
  double best_mass_ini[32];
  double best_beta_ini[32];
  double best_mom_ini[32];
  double best_dE_ini[32];
  bool best_tofs_coin_ini[32];
  int best_charge_ini[32];
  int best_layer_ini[32];
  
  double best_pv_pi[32];
  double best_mass_pi[32];
  double best_beta_pi[32];
  double best_mom_pi[32];
  bool best_tofs_coin_pi[32];
  int best_layer_pi[32];
  
  double best_rchi_pi[32];
  double best_mass_ini_pi[32];
  double best_beta_pi_ini[32];
  double best_mom_pi_ini[32];
  bool best_tofs_coin_pi_ini[32];
  int best_layer_pi_ini[32];
  
  for(int bar=0;bar<32;bar++)
    {
      best_pv[bar]=-9999;
      best_mass[bar]=-1;
      best_beta[bar]=-1;
      best_plen[bar]=-1;
      best_tof[bar]=-1;
      best_mom[bar]=-1;
      best_dE[bar]=-1;
      best_tofs_coin[bar]=false;
      best_charge[bar]=0;
      best_rchi[bar]=9999;
      
      best_mass_ini[bar]=-1;
      best_beta_ini[bar]=-1;
      best_mom_ini[bar]=-1;
      best_dE_ini[bar]=-1;
      best_tofs_coin_ini[bar]=false;
      best_charge_ini[bar]=0;
      
      best_pv_pi[bar]=-9999;
      best_mass_pi[bar]=-1;
      best_beta_pi[bar]=-1;
      best_mom_pi[bar]=-1;
      best_tofs_coin_pi[bar]=false;
      
      best_rchi_pi[bar]=9999;
      best_mass_ini_pi[bar]=-1;
      best_beta_pi_ini[bar]=-1;
      best_mom_pi_ini[bar]=-1;
      best_tofs_coin_pi_ini[bar]=false;
      
    }


  GFDetPlane* PlaneBeforeTR1 = list_Plane[TR1X];

  
  

  std::vector<int> n_status(5,0); // id 0 => q=-1 
  std::vector<int> n_status_rej(5,0); // id 0 => q=-1 
  //std::vector<TMatrixT<double> > MCtruth;

  //std::cout<<"vtracks max size "<<Vtracks.max_size()<<std::endl;

  //double reso_y[5] = {0.038,0.01,0.038,0.01,5.};
  //double reso_x[5] = {0.038,0.01,0.038,0.01,1.5};
  //int nb_det = 5;
  //int nb_track=-1;
  //std::map<int,double> DistToTR1;
  if(RecoEvent.seed_kalman_tracks_fromD2.size()>500)
    std::cout<<"vtracks size ="<<RecoEvent.seed_kalman_tracks_fromD2.size()<<std::endl;
  
  
  if(RecoEvent.seed_kalman_tracks_fromD2.size()>2000)
    return -1;

	
#ifdef DEBUG_KALMAN
  std::cout<<"1 Kalman : RecoEvent.TrackDAF.size() = "<<RecoEvent.TrackDAF.size()<<endl;//" RecoEvent.pre_tracks.size()="<<RecoEvent.pre_tracks.size()<<std::endl;
#endif
  int ntrack = -1;
  for(std::map<int,std::vector<double> >::const_iterator it_track = RecoEvent.TrackDAF.begin(),it_track_end = RecoEvent.TrackDAF.end();it_track!=it_track_end;++it_track)
    {
      ntrack++;
      Vtracks->clear();
      
#ifdef DEBUG_KALMAN
      //#ifdef DEBUG_KALMAN2
      std::cout<<"start #"<<ntrack;//<<std::endl;
      //Vtracks->trackReps->Print();
      std::cout<<" | hitsize "<<Vtracks->fHits.size()<<" "<<Vtracks->fBookkeeping.size()<<" "<<Vtracks->fRepAtHit.size()<<" "<<Vtracks->fHits.size()<<std::endl;
      //Vtracks->fCand.Print();
      std::vector<std::string> name_det;
#endif
      
      int id_track = it_track->first;
#ifdef DEBUG_KALMAN2
      std::cout<<" Id_track: "<<id_track<<std::endl;
#endif
      //std::vector<std::vector<GFHypRecoHit*> > tmp_hits(RecoEvent.ListHitsDAF[id_track]);
      std::map<int,int> ListNamePlane;
      int temp_det_i=0;
      for(unsigned int ii = 0; ii<RecoEvent.ListHitsDAF[id_track].size();++ii)
	{
	  //std::cout<<" LayerDet :"<<ii<<" "<<RecoEvent.ListHitsDAF[id_track][ii].size()<<std::endl;
	  if(RecoEvent.ListHitsDAF[id_track][ii].size()!=0)
	    {
	      ListNamePlane[ii]=temp_det_i;
	      for(unsigned int temp_single_hit = 0;temp_single_hit<RecoEvent.ListHitsDAF[id_track][ii].size();++temp_single_hit)
		Vtracks->addHit(RecoEvent.ListHitsDAF[id_track][ii][temp_single_hit],ii,temp_single_hit,0.,ii);
	      //std::cout<<" HitDet :"<<temp_det_i<<std::endl; 
	      
	      ++temp_det_i;
	    }
	}
      

      /*int temp_layer_id = 0;
      for(unsigned int ii = 0; ii<RecoEvent.ListIdHitsDAFCluster[id_track].size();++ii)
	if(RecoEvent.ListIdHitsDAFCluster[id_track][ii].size()!=0)
	  {
	    ListNamePlane[ii]=temp_det_i;
	    std::cout<<" Layer :"<< temp_det_i<<std::endl;
	    for(unsigned int temp_cluster_hit = 0;temp_cluster_hit<RecoEvent.ListIdHitsDAFCluster[id_track][ii].size();++temp_cluster_hit)
	      {
		std::cout<<" Cluster id :"<<temp_layer_id<<" Hits Id:";
		for(unsigned int temp_single_hit = 0;temp_single_hit< RecoEvent.ListIdHitsDAFCluster[id_track][ii][temp_cluster_hit].size();++temp_single_hit)
		  {
		    int id_hit = RecoEvent.ListIdHitsDAFCluster[id_track][ii][temp_cluster_hit][temp_single_hit];
		    Vtracks->addHit(RecoEvent.ListHitsDAF[id_track][ii][id_hit],ii,id_hit,0.,temp_layer_id);
		    std::cout<<" "<<id_hit;
		  }
		std::cout<<std::endl;
		++temp_layer_id;
	      }
	    ++temp_det_i;
	    }*/

#ifdef DEBUG_KALMAN
      std::cout<<" NPlane in KalmanDAF from ListHitsDAF :"<<temp_det_i<<endl;//" "<< temp_layer_id<<endl;
      Vtracks->fCand.Print();
      
#endif
      double time_of_flight = 0;
      int PID =-10;
      bool BigTofSide=false;
      std::vector<double> track_seed(it_track->second);
      std::vector<double> track_state((RecoEvent.Detectors_UTr.find("Track")->second).find(id_track)->second);
      
      PID = static_cast<int>(track_state[6]);
      bar_Track = static_cast<int>(track_state[7]);
      TofsBar   = static_cast<int>(track_state[8]);
      CoinWindow = track_state[9];
      if(PID>0)
	{
	  time_of_flight = (RecoEvent.Detectors_UTr.find("TOF")->second).find(id_track)->second[3];
	  x_TOF = (RecoEvent.Detectors_UTr.find("TOF")->second).find(id_track)->second[0];
	  y_TOF = (RecoEvent.Detectors_UTr.find("TOF")->second).find(id_track)->second[1];
	  z_TOF = (RecoEvent.Detectors_UTr.find("TOF")->second).find(id_track)->second[2];
	  E_TOF = (RecoEvent.Detectors_UTr.find("TOF")->second).find(id_track)->second[4];
	}
      else
	{
	  time_of_flight = (RecoEvent.Detectors_UTr.find("TOF")->second).find(id_track)->second[3];
	  x_TOF = (RecoEvent.Detectors_UTr.find("TOF")->second).find(id_track)->second[0];
	  y_TOF = (RecoEvent.Detectors_UTr.find("TOF")->second).find(id_track)->second[1];
	  z_TOF = (RecoEvent.Detectors_UTr.find("TOF")->second).find(id_track)->second[2];
	  E_TOF = (RecoEvent.Detectors_UTr.find("TOF")->second).find(id_track)->second[4];
	  BigTofSide=true;
	  //std::cout<<"Kalman pion filled "<<" "<<PID<<std::endl;
	}
	
#ifdef DEBUG_KALMAN
	std::cout<<" "<<PID<<endl;
#endif
	ProtonNewOrOld = 0;
	
	if(!BigTofSide)
	  {
	    NeigHit=track_state[10];
	  }
	else
	  NeigHit=0;
	
	TVector3 init_p(track_state[0],track_state[1],track_state[2]);

#ifdef DEBUG_KALMAN
	std::cout<<"init_p : "; init_p.Print();
#endif
   
	mom_Track=init_p;
	mom_Track_ini = init_p;
	beta_Track_ini=track_state[3];
	mass_Track_ini=track_state[4];
	rchi_Track_ini=track_state[5];
	charge_Track_ini=PID;
	
	
	//if(!BigTofSide&&PID==1)
	//       std::cout<<"Filled "<<((RecoEvent.Detectors_UTr.find("TR1")->second).find(id_track)->second[3])<<" "<<
	// 	PID<<" BigTof "<<BigTofSide<<" NewOrOld="<<(RecoEvent.Detectors_UTr.find("TOFInfo")->second).find(id_track)->second[1]<<std::endl;
	
	//       if(PID>0)
	// 	std::cout<<"--------------- "<<((RecoEvent.Detectors_UTr.find("TR1")->second).find(id_track)->second[3])<<" "<< time_of_flight<<std::endl;



#ifdef DEBUG_KALMAN
	std::cout<<"Momentum seed ("<<track_state[0]<<","<<track_state[1]<<", "<<track_state[2]<<") TOF = "<<time_of_flight<<std::endl;
#endif
	// #ifdef OLD
	// 	    for(int layer = 0;layer<z.size();layer++)
	// 	      {
	// 		std::cout<<" z "<<z[layer]<<" ";
	// 	      }
	// 	    std::cout<<" "<<z.size()<<std::endl;
	// #endif
      
	    
	double charge = 0.;
	double mass = 0.;
	int PDG = 0;
	{
	  if(PID == -10 )
	    cout<<"E> No pid found on TOF+ data !"<<endl;
	  else if(PID == 1)
	    {
	      charge = 1.;
	      
	      if(init_p.Mag()<=4.5)
		{
		  mass = 0.9383;
		  PDG = 2212;
		  AnaHisto->h_detection_eff_proton->Fill("proton",1.);
		}
	      else if(init_p.Mag()>4.5 && init_p.Mag()<=7.5)
		{
		  mass = 1.875613;
		  PDG = 10000;
		  AnaHisto->h_detection_eff_proton->Fill("deutron",1.);
		}
	      else if(init_p.Mag()>7.5)
		{
		  mass =  2.80925;
		  PDG = 10001;
		  AnaHisto->h_detection_eff_proton->Fill("triton",1.);
		}
	      else
		{
		  mass = 0.9383;
		  PDG = 2212;
		  std::cout<<"strange "<<__LINE__<<std::endl;
		}
	    }
	  else if(PID == 2)
	    {
	      cout<<"E> No matched pid !"<<PID<<endl;
	      continue;
	    }
	  else if(PID == 4)
	    {
	      charge = 2.;
	      if(init_p.Mag()<9.4)
		{
		  AnaHisto->h_detection_eff_alpha->Fill("He3",1.);
		  mass =  3.0160293*0.931478;
		  PDG = 10003;
		}
	      else
		{
		  mass = 3.72738;
		  PDG = 10002;
		  AnaHisto->h_detection_eff_alpha->Fill("He4",1.);
		}
	      
	    }
	  else if(PID == 9)
	    {
	      AnaHisto->h_detection_eff_alpha->Fill("Li6",1.);
	      charge = 3.;
	      mass = 5.6015194;
	      PDG = 10004;
	    }
	  else if(PID == -1)
	    {
	      AnaHisto->h_detection_eff_pion->Fill("pi-",1.);
	      charge = -1.;
	      mass = 0.1396;
	      PDG = -211;
	    }
	  else
	    {
	      cout<<"E> No matched pid !"<<PID<<endl;
	      continue;
	    }
	}
	
	//double DistToTR1;
	
	//DistToTR1=TMath::Sqrt(track_seed[3]/10.*track_seed[3]/10.+track_seed[4]/10.*track_seed[4]/10.+(track_seed[5]/10.-att.TOFStart_Z/10.0)*(track_seed[5]/10.-att.TOFStart_Z/10.0));
	
	//if(att.G4_simu==false)
	//  DistToTR1+=5; // distance from Start to TR1 
	
	//GFDetPlane init_plane((temp_listhits.front())->getDetPlane(Vtracks->getCardinalRep()));
	TVector3 init_point(track_seed[3],track_seed[4],track_seed[5]);

#ifdef DEBUG_KALMAN
	std::cout<<"init_p:"<<std::endl;
	init_p.Print();
	init_p.Unit().Print();
	std::cout<<"init_point :"<<std::endl;
	init_point.Print();
	//std::cout<<"init_plane :"<<std::endl;
	//init_plane.Print();
	
#endif
	double seed_Mom_Mag = init_p.Mag();
	if(TMath::Abs(seed_Mom_Mag)<1e-9)
	  {
	    std::cout<<"E> Seed Momemtum with TVector3 is zero !"<<std::endl;
	    //delete Vtracks;
	    //Vtracks = 0 ;
	    continue;
	  }
	
	//bool hermit_interpolation = true;
	//bool adding_noise_multiscattering = true;
	//bool B_interpolation = true;
	//bool adding_noise_multiscattering = false;//material correction
	//bool hermit_interpolation = true;
	//bool adding_noise_multiscattering = true;
	//bool B_interpolation = true;
	//AbsTrackRep* rep = 0;
	if(att.G4_simu==false)
	  {
	    
	    //rep = new HypLSLTrackFieldMapRep_v2 (init_point, init_p,posErr,momErr,charge,mass,att.StepingInField,hermit_interpolation,rotation,att.Aladin_Angle*TMath::Pi()/180.,adding_noise_multiscattering,B_interpolation);
	    //cout<<charge<<std::endl;
	    //rep = new HypLSLTrackFieldMapRep_v2 (att.StepingInField,hermit_interpolation,adding_noise_multiscattering,B_interpolation);
	    dynamic_cast<RKTrackRep*>(rep)->Init(init_point,init_p,PDG);
	  }
	else // SIMULATION 
	  {
	    //rep = new HypLSLTrackFieldMapRep_v2 (init_point, init_p,posErr,momErr,charge,mass,att.StepingInField,hermit_interpolation,rotation,adding_noise_multiscattering,B_interpolation,true);
	    //dynamic_cast<HypLSLTrackFieldMapRep_v2*>(rep)->init(init_point,init_p, posErr,momErr,sState_MC);
	    //rep = new HypLSLTrackFieldMapRep_v2 (att.StepingInField,hermit_interpolation,adding_noise_multiscattering,B_interpolation);
	    dynamic_cast<RKTrackRep*>(rep)->Init(init_point,init_p,PDG);
	  }

	if(rep==NULL)
	  std::cout<<"E> no Rep build"<<std::endl;
#ifdef DEBUG_KALMAN	
	else	
	  std::cout<<"rep done"<<std::endl;
#endif

	//Vtracks->setShift(DistToTR1);
	Vtracks->setShift(0.);
	Vtracks->setTOF(time_of_flight);
	
	//Track* Vtracks = new Track (rep,DistToTR1,time_of_flight);
	//std::cout<<"DistToTR1 "<<DistToTR1<<std::endl;
	
#ifdef DEBUG_KALMAN	  
	std::cout<<"Tkalman ref plane:";
	Vtracks->getCardinalRep()->getReferencePlane().Print();
	
     	std::cout<<"test Rep start point"<<std::endl;
	std::cout<<(rep->getReferencePlane()).getO().Z();
	rep->getState().Print();
	//rep->getStartState().Print();
	rep->getCov().Print();
	std::cout<<" -- --"<<std::endl;
#endif

      
#ifdef DEBUG_KALMAN2
      std::cout<<"track n'"<<ntrack<<std::endl;
#endif

      try
	{
	  Kalman.processTrack(Vtracks);
	  //Daf.processTrack(Vtracks);
	  //fitter.smoothing(Vtracks[itr]);
	}
      catch (GFException& e)
	{
	  std::cout<<"*** FITTER EXCEPTION ***"<<std::endl;
	  std::cout<<e.what()<<std::endl;
	  Vtracks->getTrackRep(0)->setStatusFlag(1);
	}

      TVector3 p3,posRef;
      TMatrixT<double> covRef;
      try
	{
	  Vtracks->getPosMomCov(*PlaneBeforeTR1,posRef,p3,covRef);
	}
      catch (GFException& e)
	{
	  std::cout<<"Vtracks->getPosMomCov ==> Exception "<<std::endl;
	  std::cout<<e.what()<<std::endl;
	  Vtracks->getTrackRep(0)->setStatusFlag(1);
	}
      
#ifdef DEBUG_KALMAN2
      std::cout<<"SUCESSFULL FIT!"<<std::endl;
      std::cout<<"track n'"<<ntrack<< " / "<<Vtracks->getTrackRep(0)->getStatusFlag()<<std::endl;
#endif
       
      if(Vtracks->getTrackRep(0)->getStatusFlag()==0)
	{
	  
#ifdef DEBUG_KALMAN2
	  Vtracks->getTrackRep(0)->Print();
#endif
	  
	  
#ifdef DEBUG_KALMAN2
	  TVector3 p3_old;
	  if(att.back_tracking)
	    {
	      GFDetPlane plane(TVector3(0,0,40.-att.Aladin_CenterToTarget/10.),TVector3(1,0,0),TVector3(0,1,0));
	      p3_old=Vtracks->getTrackRep(0)->getMom(plane);
	    }
	  else
	    p3_old=Vtracks->getMom();
	  double p_old=p3_old.Mag();
	  std::cout<<" mom before ROT :"; p3_old.Print();

#endif

#ifdef DEBUG_KALMAN2
	  std::cout<< "mom ="<<p_old ;
	  p3_old.Print();
#endif
 

	  double chi2=Vtracks->getChiSqu();
	  //int ndf=static_cast<int>(Vtracks->getNDF());
	  double ndf=Vtracks->getNDF();
#ifdef DEBUG_KALMAN2
	  std::cout<<" / chi2 ="<<chi2<<" / ndf ="<<ndf<<std::endl;
#endif
 
	  //double p_value = TMath::Prob(chi2,ndf);
	  double p_value = MathKalmanFRS::Prob()(chi2,ndf);
	  double p_value2= 1.-p_value;
 
#ifdef DEBUG_DAISUKE
	  double p_value_new = TMath::Prob(chi2,2*4-5);
	  double p_value2_new= 1.-p_value_new;
#endif
	  AnaHisto->h_pv->Fill(p_value2);
	  AnaHisto->h_chi2->Fill(chi2);

	  //h_chi2_smooth->Fill(Vtracks[itr]->getChiSquSmooth());
 
#ifdef DEBUG_KALMAN2
	  std::cout<< " / p_value ="<<1-p_value<<std::endl;
#endif

	  double charge = Vtracks->getCharge();
	  n_status[static_cast<int>(charge)+1]++;

#ifdef DEBUG_DAISUKE
	  std::cout<< " / p_value ="<<1-p_value<< " / p_value_new ="<<1-p_value_new<<std::endl;
#endif
	  if(charge>0)
	    {
	      AnaHisto->hd_chi[0]->Fill(chi2/ndf);
	      AnaHisto->hd_pv[0]->Fill(p_value2);
#ifdef DEBUG_DAISUKE
	      std::cout<<"TOFp / chi2 ="<<chi2<<" / ndf ="<<ndf<<std::endl;
#endif
	    }

	  else
	    {
	      AnaHisto->hd_chi[1]->Fill(chi2/ndf);
	      AnaHisto->hd_pv[1]->Fill(p_value2);
#ifdef DEBUG_DAISUKE
	      std::cout<<"BigTof / chi2 ="<<chi2<<" / ndf ="<<ndf<<std::endl;
#endif
	    }


#ifdef DEBUG_DAISUKE
	  std::cout<<"TKALMAN FILTER Come  "<<__LINE__<<std::endl;
#endif

	  GFDetPlane InitPlane(RecoEvent.ListHitsDAF[id_track][TR1X][0]->getDetPlane(rep));

	  TVector3 TR1x_Pos(Vtracks->getPos(InitPlane));
	  HitTR1 = TR1x_Pos;

	  double DistToTR1=(TR1x_Pos-Plane_time).Mag();

	  TVector3 temp_pos,temp_mom;
	  if( GFTools::getSmoothedPosMom(Vtracks,0,ListNamePlane[TR1X],temp_pos,temp_mom) )
	    {
	      HitTR1=temp_pos;
	      DistToTR1=(HitTR1-Plane_time).Mag();
	    }
	  if( GFTools::getSmoothedPosMom(Vtracks,0,ListNamePlane[DC1X],temp_pos,temp_mom) )
	    HitDC1=temp_pos;
	  if( GFTools::getSmoothedPosMom(Vtracks,0,ListNamePlane[TR2X],temp_pos,temp_mom) )
	    HitTR2=temp_pos;
	  if(charge < 0 )
	    {
	      if( GFTools::getSmoothedPosMom(Vtracks,0,ListNamePlane[DC2X],temp_pos,temp_mom) )
		HitDC2=temp_pos;
	      if( GFTools::getSmoothedPosMom(Vtracks,0,ListNamePlane[TOFP],temp_pos,temp_mom) )
		HitTOF=temp_pos;
	    }
	  else
	    {
	      if( GFTools::getSmoothedPosMom(Vtracks,0,ListNamePlane[DC2Y],temp_pos,temp_mom) )
		HitDC2=temp_pos;
	      if( GFTools::getSmoothedPosMom(Vtracks,0,ListNamePlane[DC2YP],temp_pos,temp_mom) )
		HitTOF=temp_pos;
	    }
	  HitDC1p=temp_pos;
	  

	  
	  double p = p3.Mag();
	  
	  //if(att.G4_simu==false)
	  //DistToTR1+=5.; // distance from Start to TR1 

	  double Path_length = Vtracks->getTotalLength();
	  double Path_lengthB = Vtracks->getTotalLengthBack();

	  double Path_time = Vtracks->getTotalTime();
	  double Path_timeB = Vtracks->getTotalTimeBack();
	  //double Path_lengthMean = Path_length + Path_lengthB;
	  double Path_lengthMean = Path_lengthB + Path_lengthB;
	  Path_lengthMean/=2.;
	  Path_lengthMean+= DistToTR1;
	  // 	  if(ndf==1)
	  // 	    Path_lengthMean=Path_length;
	    
	  AnaHisto->h_Path->Fill(Path_length);
	  AnaHisto->h_Path_Back->Fill(Path_lengthB);
	  AnaHisto->h_MeanPath->Fill(Path_lengthMean);
	  AnaHisto->h_dpath->Fill(Path_lengthB-Path_length);

	  //TMatrixT<double> statePredDummy(5,1);
	  double new_Length = 0.;
	  double new_Time = 0.;

	  LayerDetectorType s_type[10] = {TR0X,TR1X,DC1X,TR2X,DC2X,DC2XP,TOFP,DC2Y,DC2YP,TFWN};
	  int NplaneDetector = 7;
	  if(charge>0)
	    {
	      NplaneDetector = 10;
	    } //=======
// 	  LayerDetectorType s_type[6] = {TR1X,DC1X,TR2X,DC2X,DC2XP,TOFP};
// 	  if(charge<0)
// 	    {
// 	      s_type[3]=DC2Y;
// 	      s_type[4]=DC2YP;
// 	      s_type[5]=TFWN;
//	    }
// >>>>>>> start git repo from code in home/gsi
  
	  dynamic_cast<RKTrackRep*>(rep_length)->Init(HitTR1,p3,PDG);

	  for(int plane_id = 0;plane_id<NplaneDetector;++plane_id)
// =======
// 	  for(int plane_id = 0;plane_id<6;++plane_id)
// >>>>>>> start git repo from code in home/gsi
	    {
	      //std::cout<<" DetName :"<<name_type[i]<<std::endl<<"Pos:"<<std::endl;
	      //GFDetPlane goal_pl(TVector3(1.,1.,z[i]),TVector3(1.,0.,0.),TVector3(0.,1.,0.));
	      
	      TMatrixT<double> state_Length(5,1);
	      TMatrixT<double> cov_Length(5,5);
	      GFDetPlane* goal_pl_Length = list_Plane[s_type[plane_id]];
	      
	      assert(goal_pl_Length!=0);
	      
	      try
		{
		  rep_length->extrapolate(*goal_pl_Length,state_Length,cov_Length);
		} 
	      catch(GFException& e)
		{
		  e.what();
		  std::cerr<<"Error extrapolating in the Length calculation"<<std::endl;
		  std::cerr<<"Going to plane #"<<plane_id<<" "<<s_type[plane_id]<<std::endl;
		  break;
		}
	      rep_length->setData(state_Length, *goal_pl_Length, &cov_Length);
	    }

	  double rapidity = 0.5 * TMath::Log(( TMath::Sqrt(rep_length->getMass()*rep_length->getMass() + p3.Mag2()) +p3.Z())/( TMath::Sqrt(rep_length->getMass()*rep_length->getMass() + p3.Mag2()) - p3.Z())); 
	  for(unsigned int vol = 0 ; vol < att.name_GeoVolumes.size() ; ++vol)
	    {
	      AnaHisto->Material_XX0_y[vol]->Fill(rapidity,rep_length->getXX0(att.name_GeoVolumes[vol]));
	      AnaHisto->Material_dE_y[vol]->Fill(rapidity,rep_length->getDE(att.name_GeoVolumes[vol]));
	    }



#ifdef VERBOSE_EVE

	  std::cout<<"$$$$$$$$$$$$-----------------------$$$$$$$$$$$$$$$"<<std::endl;
	  std::cout<<"DAF Length :ListTrackPoints :"<<std::endl;
	  std::cout<<"TEveLine* EveTrack_FRS1 = new TEveLine(TrackRepList_FRS,"<<rep_length->EveListTrackPoint.size()<<");"<<std::endl;
	  for(unsigned int i = 0;i<rep_length->EveListTrackPoint.size();++i)
	    {
	      std::cout<<std::setprecision(10)<<"EveTrack_FRS1->SetPoint("<<i<<","<<rep_length->EveListTrackPoint[i].X()<<","<<rep_length->EveListTrackPoint[i].Y()<<","<<rep_length->EveListTrackPoint[i].Z()<<");"<<std::endl;
	    }	
      
	      
	  std::cout<<"Fitted Hit :"<<std::endl;
	  std::cout<<"TEvePointSet* EveHitFitter1 = new TEvePointSet(HitFitter1,"<<5<<");"<<std::endl;
	  std::cout<<std::setprecision(10)<<"EveHitFitter1->SetPoint("<<0<<","<<HitTR1.X()<<","<<HitTR1.Y()<<","<<HitTR1.Z()<<");"<<std::endl;
	  std::cout<<std::setprecision(10)<<"EveHitFitter1->SetPoint("<<1<<","<<HitDC1.X()<<","<<HitDC1.Y()<<","<<HitDC1.Z()<<");"<<std::endl;
	  std::cout<<std::setprecision(10)<<"EveHitFitter1->SetPoint("<<2<<","<<HitTR2.X()<<","<<HitTR2.Y()<<","<<HitTR2.Z()<<");"<<std::endl;
	  std::cout<<std::setprecision(10)<<"EveHitFitter1->SetPoint("<<3<<","<<HitDC2.X()<<","<<HitDC2.Y()<<","<<HitDC2.Z()<<");"<<std::endl;
	  std::cout<<std::setprecision(10)<<"EveHitFitter1->SetPoint("<<4<<","<<HitTOF.X()<<","<<HitTOF.Y()<<","<<HitTOF.Z()<<");"<<std::endl;

	  std::cout<<"DAF Track :"<<std::endl;
	  std::cout<<"TEveRecTrack *DAF_track_rc = new TEveRecTrack();"<<std::endl;
	  std::cout<<"DAF_track_rc->fV.Set("<<HitTR1.X()<<","<<HitTR1.Y()<<","<<HitTR1.Z()<<");"<<std::endl;
	  std::cout<<"DAF_track_rc->fP.Set("<<p3.X()<<","<<p3.Y()<<","<<p3.Z()<<");"<<std::endl;
	  //std::cout<<"int sign = -"<<list_MC[current_idMC]->Charge<<std::endl;
	  std::cout<<"DAF_track_rc->fSign = -"<<charge<<std::endl;
	  std::cout<<"$$$$$$$$$$$$-----------------------$$$$$$$$$$$$$$$"<<std::endl;

#endif


	  new_Length = rep_length->getTotalLength();
	  new_Time = rep_length->getTotalTime();

	  double new_Length_check = rep_length->getTotalLengthBack();
	  if(TMath::Abs(new_Length_check)>1e-4)
	    std::cout<<" E> Strange New_Length_check != 0 "<<new_Length_check<<std::endl;

	  {
	    dynamic_cast<RKTrackRep*>(rep_length)->Init(HitTR1,p3,PDG);
	    GFDetPlane* goal_pl_Length = list_Plane[TR0X];
	    TMatrixT<double> state_Length(5,1);
	    TMatrixT<double> cov_Length(5,5);

	    rep_length->switchDirection();
	    try
	      {
		rep_length->extrapolate(*goal_pl_Length,state_Length,cov_Length);
	      } 
	    catch(GFException& e)
	      {
		e.what();
		std::cerr<<"Error extrapolating in the Length calculation TR0 TR1"<<std::endl;
	      }
	  }
	   
	  double new_Length_TR1toTR0 = -rep_length->getTotalLengthBack();
	  double new_Time_TR1toTR0 = -rep_length->getTotalTimeBack();
	  double new_Length_TR1toTR0_check = rep_length->getTotalLength();
	  if(TMath::Abs(new_Length_TR1toTR0_check)>1e-4)
	    std::cout<<" E> Strange New_Length_TR1toTR0_check != 0 "<<new_Length_TR1toTR0_check<<std::endl;

	  double new_Total_Length = new_Length_TR1toTR0+new_Length;
	  double new_Total_Time = new_Time_TR1toTR0+new_Time;
	  
	  AnaHisto->h_Path_newL->Fill(new_Total_Length);
	  AnaHisto->h_Path_LTr1toTOF->Fill(new_Length);
	  AnaHisto->h_Path_LTr0toTr1->Fill(new_Length_TR1toTR0);
	  AnaHisto->h_Path_check1->Fill(new_Length_check);
	  AnaHisto->h_Path_check2->Fill(new_Length_TR1toTR0_check);
	  
	  AnaHisto->h_dpath_new->Fill(Path_lengthMean-new_Total_Length,0);
	  AnaHisto->h_dpath_new->Fill(Path_length-new_Length,1);
	  AnaHisto->h_dpath_new->Fill(Path_lengthB-new_Length,2);
	  AnaHisto->h_dpath_new->Fill(30.*beta_Track_ini*Vtracks->getTimeOfFlight()-new_Total_Length,3);
	  
	  AnaHisto->h_dpath_new->Fill(DistToTR1-new_Length_TR1toTR0,5);
	  
	  //try
	      //bin/{
	      //new_Length = DistToTR1 + Vtracks->getCardinalRep()->extrapolate(LastDetPlane,statePredDummy);
	      //}
	      //catch (FitterException& e)
	      //{
	      //std::cout<<"*** EXTRAPOLATION EXCEPTION Length Calculation ***"<<std::endl;
	      //std::cout<<e.what()<<std::endl;
	      //}
	  
	  //AnaHisto->h_dpath2_chi2->Fill(Path_lengthMean-new_Length,chi2);
	  //AnaHisto->h_dpath2_pv->Fill(Path_lengthMean-new_Length,p_value2);

	  double time_of_flight= Vtracks->getTimeOfFlight();//+6.;
	  
	  AnaHisto->h_Time->Fill(new_Total_Time);
	  AnaHisto->h_Time_Estimation->Fill(time_of_flight - new_Total_Time,0);
	  AnaHisto->h_Time_Estimation->Fill(time_of_flight - Path_time,1);
	  AnaHisto->h_Time_Estimation->Fill(time_of_flight - Path_timeB,2);
	  AnaHisto->h_Time_Estimation->Fill(time_of_flight - 0.5*(Path_time+Path_timeB),3);

	  //double beta = 1./30.*Path_lengthMean/time_of_flight; ///beta = dL[m]/(c[m/s]* dT[s]) = dL[mm]/(300[mm/ns] * dT[ns])
	  double beta = 1./30.*new_Total_Length/time_of_flight; ///beta = dL[m]/(c[m/s]* dT[s]) = dL[mm]/(300[mm/ns] * dT[ns])
	  //double gamma = 1./sqrt(1.-beta*beta);
	  double mass = p*sqrt(1./(beta*beta)-1.);

	  if(beta < 0.)
	    {
	      cout<<"E> Beta Negative ! "<<Path_lengthMean<<" "<<time_of_flight<<" | "<<Path_length<<" "<<Path_lengthB<<" "<<endl;
	      cout<<charge<<" "<<p<<" "<<chi2<<" "<<ndf<<" "<<p_value2<<endl;
	    }
 
#ifdef DEBUG_KALMAN2
	  std::cout<<"charge: "<<charge<<" / time = "<<time_of_flight /*Detectors_UTr["TOF"][id_track][3]*/ <<" / Path = "<<Path_length<<" / Path back "<<Path_lengthB<<" PathMean = "<<Path_lengthMean<<" New Length = "<<new_Length <<" DTR1STart"<< DistToTR1 <<" / Beta = "<<beta<<" / Mass = "<<mass<<std::endl;
#endif
 
 
 
	  AnaHisto->h_beta->Fill(beta);
	  AnaHisto->h_Mass_All->Fill(mass);
	  AnaHisto->h_Mass_charge_All->Fill(mass,charge);
	  AnaHisto->h_beta_mom->Fill(p,beta);
	  AnaHisto->h_path_tof->Fill(Path_lengthMean/30.,time_of_flight);
	  
	  AnaHisto->h_pv_mom->Fill(p,p_value2);
	  AnaHisto->h_pv_beta->Fill(beta,p_value2);
	  AnaHisto->h_pv_mass->Fill(mass,p_value2);

#ifdef DEBUG_KALMAN
	  std::cout<<"TKalmanFilter_FRS "<<__LINE__<<std::endl;
#endif	  
	  if(p_value2<.75)
	    {
#ifdef DEBUG_KALMAN
	      std::cout<<"TKalmanFilter_FRS "<<__LINE__<<std::endl;
#endif	  

	      AnaHisto->h_beta2->Fill(beta);
	      AnaHisto->h_Mass_All2->Fill(mass);
	      AnaHisto->h_beta_mom2->Fill(p,beta);
	      AnaHisto->h_Mass_charge_All2->Fill(mass,charge);
	    }
	  if(p_value2<.1)
	    {
#ifdef DEBUG_KALMAN
	      std::cout<<"TKalmanFilter_FRS "<<__LINE__<<std::endl;
#endif	  

	      AnaHisto->h_beta3->Fill(beta);
	      AnaHisto->h_Mass_All3->Fill(mass);
	      AnaHisto->h_beta_mom3->Fill(p,beta);
	      AnaHisto->h_Mass_charge_All3->Fill(mass,charge);
	      AnaHisto->h_mom_tof_cut->Fill(p,time_of_flight);
	      AnaHisto->h_path_mom_cut->Fill(Path_lengthMean/30.,p);
	      AnaHisto->h_path_tof_cut->Fill(Path_lengthMean/30.,time_of_flight);
	    }
 
	  //if(charge<0.)
	  //  mass=0.139;
 
	  //std::vector<double> mom_vec(13,0.);
	  //std::vector<double> mom_vec(32,0.);
	  //std::vector<double> mom_vec(43,0.);
	  //std::vector<double> mom_vec(77,0.); // 30 12 1010 
	  //std::vector<double> mom_vec(85,0.); // 30 12 1010 
	  std::vector<double> mom_vec(100,0.); // 30 12 1010 
	  mom_vec[0]=p3.X();
	  mom_vec[1]=p3.Y();
	  mom_vec[2]=p3.Z();
#ifdef DEBUG_KALMAN
	  std::cout<<"TKalmanFilter_FRS "<<__LINE__<<std::endl;
#endif	  

#ifdef OLD
#ifndef ROTATION_KALMAN
	  TVector3 temp_Pos = Vtracks->getPos();
	  Rotating(temp_Pos,1.,true);
	  mom_vec[3]=temp_Pos.X();
	  mom_vec[4]=temp_Pos.Y();
	  mom_vec[5]=temp_Pos.Z();//+att.Aladin_CenterToTarget/10.;

	  

#else
	  TVector3 point_ref(Vtracks->getPos());
	  mom_vec[3]=point_ref.X();
	  mom_vec[4]=point_ref.Y();
	  mom_vec[5]=point_ref.Z(); //cm !!
	  //	  std::cout<<"def "<<std::endl;
	  //std::cout<<"point("<<mom_vec[5]*10<<", "<<mom_vec[3]*10<<")"<<std::endl;
	  AnaHisto->hd_theWorld->Fill(mom_vec[5]*10,mom_vec[3]*10);
#endif
#endif
	  mom_vec[3]=posRef.X(); // cm
	  mom_vec[4]=posRef.Y(); // cm
	  mom_vec[5]=posRef.Z(); // cm

	  mom_vec[6]=charge;
	  mom_vec[7]=chi2;
	  mom_vec[8]=ndf;

	  mom_vec[9]=new_Total_Time;//new_Total_Length;

	  mom_vec[10]=new_Total_Length;//Path_lengthMean;//Path_lengthB;

	  mom_vec[11]=time_of_flight;
#ifdef DEBUG_DAISUKE
	  std::cout<<"p 1 = "<<Path_length<<" p2 ="<<Path_lengthB<<std::endl;
#endif

#ifdef DEBUG_DAISUKE
	  std::cout<<"TKalmanFilter_FRS "<<__LINE__<<" etofsize"<<E_TOF.size()<<" itr ="<<itr<<" x tof "<<x_TOF.size()<<" vtracks size="<<Vtracks.size()<<std::endl;
#endif	  

	  mom_vec[12]=E_TOF;

#ifdef DEBUG_DAISUKE
	  std::cout<<"TKalmanFilter_FRS "<<__LINE__<<std::endl;
#endif	  
	  
	  //cout<<" ListNamePlane :"<<ListNamePlane.size()<<endl;

	  /*GFDetPlane planeExtrap(TVector3(0,0,41.11),TVector3(1,0,0),TVector3(0,1,0));
	  HitTR1=Vtracks->getTrackRep(0)->getPos(planeExtrap);
	  Rotating(HitTR1,true);
	  
	  planeExtrap.setO(TVector3(0,0,721.14502));
	  HitTR2=Vtracks->getTrackRep(0)->getPos(planeExtrap);
	  Rotating(HitTR2,true);
	  
	  planeExtrap.setO(TVector3(0,0,355.601));
	  HitDC2=Vtracks->getTrackRep(0)->getPos(planeExtrap);
	  Rotating(HitDC2,true);*/
	  
	  mom_vec[13] = HitTR1.X();
	  mom_vec[14] = HitTR1.Y();
	  mom_vec[15] = HitTR1.Z();

	  mom_vec[16] = HitTR2.X();
	  mom_vec[17] = HitTR2.Y();
	  mom_vec[18] = HitTR2.Z();
	  
	  mom_vec[19] = HitDC2.X();
	  mom_vec[20] = HitDC2.Y();
	  mom_vec[21] = HitDC2.Z();

	  mom_vec[22] = NumDC2.X();
	  mom_vec[23] = NumDC2.Y();
	  mom_vec[24] = NumDC2.Z();

	  mom_vec[25] = bar_Track;
	  mom_vec[26]= beta;
	  mom_vec[27]= p_value2;
	  mom_vec[28]= TofsBar;

	  mom_vec[29] = x_TOF;
	  mom_vec[30] = y_TOF;
	  mom_vec[31] = z_TOF;


	  mom_vec[32] = -1;//HitDC1.X();
	  mom_vec[33] = -1;//HitDC1.Y();
	  mom_vec[34] = -1;//HitDC1.Z();
	  mom_vec[35] = -1;//HitDC1p.X();
	  mom_vec[36] = -1;//HitDC1p.Y();
	  mom_vec[37] = -1;//HitDC1p.Z();


	  mom_vec[38] = mom_Track_ini.X();
	  mom_vec[39] = mom_Track_ini.Y();
	  mom_vec[40] = mom_Track_ini.Z();

	  mom_vec[41] = rchi_Track_ini;

	  mom_vec[42] = beta_Track_ini;

	  mom_vec[43] = CoinWindow;


#ifdef DEBUG_DAISUKE
	  std::cout<<" -----------------------------------"<<std::endl;
	  std::cout<<"mom vec ("<<mom_vec[3]*10<<","<<mom_vec[4]*10<<","<<mom_vec[5]*10<<")"<<std::endl;
	  std::cout<<"mom vec ("<<mom_vec[13]<<","<<mom_vec[14]<<","<<mom_vec[15]<<")"<<std::endl;
#endif
	  //double m_range[4] ={ 0.9383,3.72738,0.1396,2.809};
	  double m_charge[4] ={ 1.,2.,-1.,2.};

#ifdef DEBUG_DASUKE
	  std::cout<<"TKalmanFilter_FRS "<<__LINE__<<std::endl;
#endif	  
	  
	  for (int i=0;i<4;i++)
	    //if (TMath::Abs(mass-m_range[i])<0.25*m_range[i] && TMath::Abs(charge-m_charge[i])<0.1)
	    if (TMath::Abs(charge-m_charge[i])<0.1)
	      {
		//double dmass = TMath::Abs(mass-m_range[i])/m_range[i];
		AnaHisto->h_Mass[i]->Fill(mass);
		AnaHisto->h_chi2_particle[i]->Fill(chi2);
		AnaHisto->h_pv_particle[i]->Fill(p_value2);
	      }
 
	  MomRef temp_mom_ref(ntrack,p_value2,mass,x_TOF,y_TOF);
#ifdef DEBUG_KALMAN
	  std::cout<<"Temp :"<<ntrack<<" "<<p_value2<<" "<<mass<<std::endl;
	  std::cout<<"MomRef :"<<temp_mom_ref.id_track<<" "<<temp_mom_ref.p_value<<" "<<temp_mom_ref.mass<<std::endl;
#endif

	  /*std::vector<double> CovSigma(5,0.);
	  for(int sigma_index =0;sigma_index<5;sigma_index++)
	  CovSigma[sigma_index]=TMath::Sqrt(Vtracks->getCardinalRep()->getCovElem(sigma_index,sigma_index));*/

#ifdef LATER
	  if(att.Debug_FRS==true && false)
	    {
	      std::vector<std::vector<std::vector<double> > > weights(Daf.getWeights());
	      std::vector<std::vector<double> > chi2_planes(Daf.getChi2Planes());
	      assert(weights[0].size() == ListNamePlane.size());
	      /*if(weights[0].size() != ListNamePlane.size())
		{
		std::cout<<"E> Weight != PlaneName "<<weights[0].size()<<" "<<ListNamePlane.size()<<" Id_track: "<<id_track<<std::endl;
		for(std::map<int,int>::const_iterator it_nameplane = ListNamePlane.begin(), it_nameplane_end = ListNamePlane.end();it_nameplane!=it_nameplane_end;++it_nameplane)
		{
		std::cout<<it_nameplane->first<<" "<<it_nameplane->second<<std::endl;
		}
		assert(weights[0].size() == ListNamePlane.size());
		}
	      */
	      int index_charge_1 = charge < 0 ? 1+static_cast<int>(charge) : static_cast<int>(charge) ;
	      for(std::map<int,int>::const_iterator it_nameplane = ListNamePlane.begin(), it_nameplane_end = ListNamePlane.end();it_nameplane!=it_nameplane_end;++it_nameplane)
		{
		  
		  std::vector<double> temp_W(weights[0][it_nameplane->second]);
		  if(temp_W.size()>1)
		    {
		      double temp_max = *std::max_element(temp_W.begin(),temp_W.end());
		      double temp_min = *std::min_element(temp_W.begin(),temp_W.end());
		      double init = 0.;
		      double temp_sum = std::accumulate(temp_W.begin(),temp_W.end(),init);
		      
		      AnaHisto->h_weigth_max[index_charge_1][it_nameplane->first]->Fill(temp_max);
		      AnaHisto->h_weigth_min[index_charge_1][it_nameplane->first]->Fill(temp_min);
		      AnaHisto->h_weigth_sum[index_charge_1][it_nameplane->first]->Fill(temp_sum);
		      AnaHisto->h_chi2_det[index_charge_1][it_nameplane->first]->Fill(chi2_planes[0][it_nameplane->second]);
		    }
		  else if(temp_W.size()==1)
		    {
		      //std::cout<<"hist:"<<AnaHisto->h_weigth_one[index_charge_1][it_nameplane->first]<<" "<<temp_W.front()<<std::endl;
		      AnaHisto->h_weigth_one[index_charge_1][it_nameplane->first]->Print();
		      AnaHisto->h_weigth_one[index_charge_1][it_nameplane->first]->Fill(temp_W.front());
		      AnaHisto->h_chi2_det_one[index_charge_1][it_nameplane->first]->Fill(chi2_planes[0][it_nameplane->second]);
		    }
		}
	      //int index_charge = charge;
	      //if(charge<0)
	      //index_charge=0;

	      /*for(int det = 0;det<4;++det)
		{
		  AnaHisto->h_ResidualX[det][index_charge]->Fill(mom_vec[77+det]);
		  AnaHisto->h_ResidualY[det][index_charge]->Fill(mom_vec[84+det]);
		  if(det==3)
		    {
		      int index_tof = charge > 0 ? det : det+1;
		      AnaHisto->h_ResidualX[index_tof][index_charge]->Fill(mom_vec[77+det]);
		      AnaHisto->h_ResidualY[index_tof][index_charge]->Fill(mom_vec[84+det]);
		    }
		    }*/	  
	    }

	  std::vector< std::vector<int>* > planes;
	  Vtracks->getHitsByPlane(planes);
	  int nPlanes = planes.size();
	  	  
	  std::vector<GFDafHit*> eff_hits;
	  
	  for(int i_plane=0; i_plane<nPlanes; i_plane++) 
	    {
	      
	      std::vector<GFAbsRecoHit*> hits;
	    
	      for(unsigned int j_hit=0; j_hit<planes.at(i_plane)->size(); j_hit++) 
		{
		  hits.push_back(Vtracks->getHit(planes.at(i_plane)->at(j_hit)) );
		}

	      GFDafHit* eff_hit = new GFDafHit(hits);
	      eff_hit->setWeights(weights[0][i_plane]);
	      eff_hits.push_back(eff_hit);
	  }
	  
	  std::vector<std::vector<double> > result_residual;
	  std::vector<std::vector<double> > sigma_residual;
	  
	  GFTools::getResidualsDAF(Vtracks,0,result_residual,sigma_residual,eff_hits);
	  for(unsigned int ii= 0 ;ii<eff_hits.size();++ii)
	    {
	      delete eff_hits[ii];
	      eff_hits[ii]= 0;
	    }
	  eff_hits.clear();

	  mom_vec[77] = result_residual[ListNamePlane[TR1X]][0]/TMath::Sqrt(sigma_residual[ListNamePlane[TR1X]][0]); //tr1x
	  mom_vec[78] = result_residual[ListNamePlane[TR2X]][0]/TMath::Sqrt(sigma_residual[ListNamePlane[TR2X]][0]); //tr2x
	  mom_vec[79] = result_residual[ListNamePlane[DC2X]][0]/TMath::Sqrt(sigma_residual[ListNamePlane[DC2X]][0]); //dc2x
	  if(charge >0)
	    mom_vec[80] = result_residual[ListNamePlane[TOFP]][0]/TMath::Sqrt(sigma_residual[ListNamePlane[TOFP]][0]); //tofx
	  else
	    mom_vec[80] = result_residual[ListNamePlane[TFWN]][0]/TMath::Sqrt(sigma_residual[ListNamePlane[TFWN]][0]); //tofx

	  mom_vec[81] = result_residual[ListNamePlane[TR1Y]][0]/TMath::Sqrt(sigma_residual[ListNamePlane[TR1Y]][0]); //tr1y
	  mom_vec[82] = result_residual[ListNamePlane[TR2Y]][0]/TMath::Sqrt(sigma_residual[ListNamePlane[TR2Y]][0]); //tr2y
	  mom_vec[83] = result_residual[ListNamePlane[DC2Y]][0]/TMath::Sqrt(sigma_residual[ListNamePlane[DC2Y]][0]); //dc2y
	  if(charge >0)
	    mom_vec[84] = result_residual[ListNamePlane[TOFP]][1]/TMath::Sqrt(sigma_residual[ListNamePlane[TOFP]][1]); //tofx
	  else
	    mom_vec[84] = result_residual[ListNamePlane[TFWN]][1]/TMath::Sqrt(sigma_residual[ListNamePlane[TFWN]][1]); //tofx
#endif 
	  	  

	  TVector3 vec_track(mom_vec[0],mom_vec[1],mom_vec[2]);
	  //TVector3 point_ref_cm(mom_vec[3]*10,mom_vec[4]*10,mom_vec[5]*10);
	  TVector3 point_ref_cm(mom_vec[3],mom_vec[4],mom_vec[5]);
	  MT_Line track_scand;
	  track_scand.From_SlopeAndPoint(vec_track,point_ref_cm);

	  double tpos[3]={0,0,0};
	  bool temp = track_scand.Get_XY_atZ(tpos);
	  // double cutfuncl = TofsBar*3.0-18.0;
	  // double cutfunch = TofsBar*3.0-12.0;
	  // bool tofs_coin = false;
	  // if(tpos[0]>cutfuncl &&tpos[0]<cutfunch)
	  bool tofs_coin=true;
	  mom_vec[98] = tofs_coin; 
	  mom_vec[99] = NeigHit;
	  if(p_value2<0.9999)
	    {

	      //TMatrixT<double> State(5,1);
	      //TMatrixT<double> Cov(5,5);
	      //State = Vtracks->getCardinalRep()->getState();
	      //Cov=Vtracks->getCardinalRep()->getCov();
	      // for(int row=0;row<5;row++)
	      // 	{
	      // 	  mom_vec[44+row] = State[row][0];
	      // 	}
	      // for(int row=0;row<5;row++)
	      // 	for(int col=0;col<5;col++)
	      // 	  {
	      // 	    mom_vec[49+row*5+col] = Cov[row][col];
	      // 	  }
	      // TVector3 p_back(0.,0.,0.);
	      /*
		try
		{
		p_back=Vtracks->getTrackRep(0)->getMom(LastDetPlane);
		}
		catch (FitterException& e)
		{
		std::cout<<"*** EXTRAPOLATION EXCEPTION ***"<<std::endl;
		std::cout<<e.what()<<std::endl;
		}
	      */
	      for(int row=0;row<3;++row)
		mom_vec[44+row] = posRef[row];
	      for(int row=0;row<3;++row)
		mom_vec[47+row] = p3[row];

	      for(int row=0;row<6;++row)			  
		for(int col=0;col<6;++col)
		  mom_vec[50+row*6+col] = covRef[row][col];
	      bool test_smooth = false;
	      if(charge<0)
		test_smooth = GFTools::getSmoothedPosMom(Vtracks,0,ListNamePlane[TOFP],temp_pos,temp_mom);
	      else
		test_smooth = GFTools::getSmoothedPosMom(Vtracks,0,ListNamePlane[DC2YP],temp_pos,temp_mom);

	      TVector3 p_back;
	      if(test_smooth)
		p_back = temp_mom;
	      
	      //mom_vec[74] = p_back.X();
	      //mom_vec[75] = p_back.Y();
	      //mom_vec[76] = p_back.Z();
	      mom_vec[87] = p_back.X();
	      mom_vec[88] = p_back.Y();
	      mom_vec[89] = p_back.Z();
	      
	      RecoEvent.Mom_Particle[temp_mom_ref]=mom_vec;
	      //RecoEvent.Sigma[ntrack]=CovSigma;
#ifdef DEBUG_KALMAN	      
	      std::cout<<"Kalman pval "<<p_value2<<" charge "<<charge<<std::endl;
#endif
	    }
	  //std::cout<<" Mass = "<<mass<<" pval "<<p_value<<std::endl;
	  //	  std::cout<<mom_Track[itr].Z()<<" " <<mom_vec[2]<<std::endl;
	  AnaHisto->h_MomCor[0]->Fill(mom_Track.X(),mom_vec[0]);
	  AnaHisto->h_MomCor[1]->Fill(mom_Track.Y(),mom_vec[1]);
	  AnaHisto->h_MomCor[2]->Fill(mom_Track.Z(),mom_vec[2]);


	  if(charge>0)
	    {
	      // 	      if(charge==1 && ProtonNewOrOld[itr]==0)
	      //std::cout<<"Bar Id "<<bar_Track<<" NewOrOld"<<ProtonNewOrOld<<" Neig "<<NeigHit<<" charge = "<<charge<<" Pvalue"<<p_value<<std::endl;
	      //int RecalBarId = TMath::Nint((164.96 - x_TOF[itr])/29.5);
	      if(bar_Track>=32 || (ProtonNewOrOld==1&&charge==1) || NeigHit!=0) 
		{
		  //std::cout<<"end #"<<ntrack<<std::endl;
		  
		  continue;//<-Old Momentum 1 New Momentum 0
		}
	      //cout<<"passed and just be filled !"<<endl;
	      int RecalBarId=1;//bar_Track;
	      
	      if(charge>0)
		{
		  //std::cout<<"x_TOF "<<x_TOF[itr]<<" e_TOF"<<E_TOF[itr]<<std::endl;
		  // 		  if(RecalBarId!=bar_Track[itr]+1)
		  // 		    std::cout<<"Bar = "<<RecalBarId<<" "<<bar_Track[itr] -1<<std::endl;
		  if(TofsBar<0)
		    std::cout<<" Tof Start Bar not defined id="<<TofsBar<<std::endl;

		  AnaHisto->hd_TrackFrom[0]->Fill(TofsBar,tpos[0]);

		  //std::cout<<"Mass="<<mass<<" Charge="<<charge<<std::endl;
		  if(RecalBarId<0 || RecalBarId>31||p_value>1)
		    std::cout<<"Daisuke "<<charge<<" "<<RecalBarId<<" "<<p_value<<std::endl;

		  if(charge>0 &&best_pv[RecalBarId]<p_value)
		    {
		      best_pv[RecalBarId]=p_value;
		      best_mass[RecalBarId]=mass;
		      best_beta[RecalBarId]=beta;
		      best_plen[RecalBarId]=Path_lengthMean;
		      best_tof[RecalBarId]=time_of_flight;
		      best_mom[RecalBarId]=sqrt(mom_vec[0]*mom_vec[0]+mom_vec[1]*mom_vec[1]+mom_vec[2]*mom_vec[2]);
		      best_dE[RecalBarId] = mom_vec[12];
		      best_tofs_coin[RecalBarId] = tofs_coin;
		      best_charge[RecalBarId] = charge*charge;
		      best_layer[RecalBarId] = 4;
		      //cout<<mass<<endl;
		    }
		  //cout<<" best :"<<best_charge[RecalBarId]<<endl;
		}
		
	      if(charge_Track_ini>0 && rchi_Track_ini < best_rchi[RecalBarId ])
		{
		  best_rchi[RecalBarId] = rchi_Track_ini;
		  best_mass_ini[RecalBarId]=mass_Track_ini;
		  best_mom_ini[RecalBarId]=mom_Track_ini.Mag();
		  best_beta_ini[RecalBarId]=beta_Track_ini;
		  best_tofs_coin_ini[RecalBarId] = tofs_coin;
		  best_dE_ini[RecalBarId] = mom_vec[12];
		  best_charge_ini[RecalBarId] = charge_Track_ini;
		  best_layer_ini[RecalBarId] = 4.;
		}
	      
	    }
	  else
	    {
	      
	      AnaHisto->hd_TrackFrom[1]->Fill(TofsBar,tpos[0]);
	      //int RecalBarId =  (819.28 -(x_TOF -1679.26)/TMath::Cos(22.8*TMath::DegToRad()))/105;
	      int RecalBarId = 1;
// =======
// 	      int RecalBarId =  (819.28 -(x_TOF -1679.26)/TMath::Cos(22.8*TMath::DegToRad()))/105;
// >>>>>>> start git repo from code in home/gsi
	      if(TMath::Abs(charge+1)<0.1)
		{
		  //std::cout<<"x_TOF "<<x_TOF[itr]<<" e_TOF"<<E_TOF[itr]<<std::endl;
		  
		  //std::cout<<RecalBarId<<std::endl;//" "<<bar_Track[itr] -1<<std::endl;
		  //std::cout<<"p_value = "<<p_value<<std::endl;
		  if(RecalBarId<0 || RecalBarId>18||p_value>1)
		    std::cout<<"Daisuke "<<charge<<" "<<RecalBarId<<" "<<p_value<<std::endl;
		  else
		    {
		      if(TMath::Abs(charge+1)<0.1 && best_pv_pi[RecalBarId]<p_value)
			{
			  best_pv_pi[RecalBarId]=p_value;
			  best_mass_pi[RecalBarId]=mass;
			  //cout<<mass<<endl;
			  best_beta_pi[RecalBarId]=beta;
			  best_mom_pi[RecalBarId]=sqrt(mom_vec[0]*mom_vec[0]+mom_vec[1]*mom_vec[1]+mom_vec[2]*mom_vec[2]);
			  best_tofs_coin_pi[RecalBarId] = tofs_coin;
			  best_layer_pi[RecalBarId] = 4.;
			}
		    }
		}
	      if(rchi_Track_ini < best_rchi_pi[RecalBarId ])
		{
		  best_rchi_pi[RecalBarId] = rchi_Track_ini;
		  best_mass_ini_pi[RecalBarId]=mass_Track_ini;
		  best_mom_pi_ini[RecalBarId]=mom_Track_ini.Mag();
		  best_beta_pi_ini[RecalBarId]=beta_Track_ini;
		  best_tofs_coin_pi_ini[RecalBarId] = tofs_coin;
		  best_layer_pi_ini[RecalBarId] = 4.;
		}
	    }

	} //StatusFlag == 0
      else
	{
	  n_status_rej[static_cast<int>(Vtracks->getCharge())+1]++;
	}
      //delete Vtracks;
      //Vtracks = 0;
    }// for(unsigned int itr=0;itr<Vtracks.size();++itr)

  Vtracks->clear();
#ifdef OLDHITS
  for(unsigned int ii = 0;ii<handlerhit.size();++ii)
    if(handlerhit[ii]!=0)
      delete handlerhit[ii];
  handlerhit.clear();
#endif

#ifdef OUTBAR
  cout<<"Fit OK q -1 ="<<n_status[0] << " | q 1 ="<<n_status[2]<<" | q 2 ="<<n_status[3]<<" | q 3 ="<<n_status[4]<<endl;
    cout<<"Fit RJ q -1 ="<<n_status_rej[0] << " | q 1 ="<<n_status_rej[2]<<" | q 2 ="<<n_status_rej[3]<<" | q 3 ="<<n_status_rej[4]<<endl;
  for(int bar=0;bar<32;bar++)
    {
      if(best_charge[bar]>0.)
	{
	  cout<<" bar #"<<bar<<endl;
	  cout<<" pv :"<<best_pv[bar]<<" mass:"<<best_mass[bar]<<" beta:"<<best_beta[bar]<<" plen:"<<best_plen[bar]<<" tofs:"<<best_tof[bar]<<endl;
	  cout<<" mom:"<<best_mom[bar]<<" E:"<<best_dE[bar]<<" Tofscoin:"<<best_tofs_coin[bar]<<" Q:"<<best_charge[bar]<<" rchi:"<<best_rchi[bar]<<endl;
	  cout<<" iniQ:"<<best_charge_ini[bar]<<endl;
	  cout<<"================="<<endl;
	}
      if(best_pv_pi[bar]>-1.)
	{
	  cout<<"bar #"<<bar<<endl;
	  cout<<" Pi-"<<endl;
	  cout<<"pv "<<best_pv_pi[bar]<<" mom:"<<best_mom_pi[bar]<<endl;
	  cout<<"================="<<endl;
	}
    }

#endif
    



    for(int bar=0;bar<32;bar++)
      {
	if(best_mom[bar]>0)
	  {
	    if(best_charge[bar]==1)
	      {
		//	      cout<<"TKalman "<<__LINE__<<endl;
		AnaHisto->h_MassZ1BarBest->Fill(best_mass[bar],bar);
		// 	      if(best_layer[bar]==3)
		// 		AnaHisto->h_MassZ1BarBest3->Fill(best_mass[bar],bar);
		AnaHisto->hd_MassPv[0]->Fill(best_mass[bar],1-best_pv[bar]);
		//	  std::cout<<"best mass "<<bar<<" "<<best_mass[bar]<<endl;
		//	      AnaHisto->hd_BetaMom_Tp[bar]->Fill(best_mom[bar],best_beta[bar]);
		AnaHisto->hd_BetaMom_Tp_all->Fill(best_mom[bar],best_beta[bar]);
		//	      cout<<"TKalman "<<__LINE__<<endl;
	      
		if(best_tofs_coin[bar])
		  {
		    AnaHisto->h_MassZ1BarBestTofs->Fill(best_mass[bar],bar);
		    // 		  if(best_layer[bar]==3)
		    // 		    AnaHisto->h_MassZ1BarBestTofs3->Fill(best_mass[bar],bar);
		    AnaHisto->hd_BetaMom_Tp_allTofs->Fill(best_mom[bar],best_beta[bar]);

		    if(best_mass[bar]>0.5 && best_mass[bar]<1.5)
		      {
			double tote = sqrt(best_mom[bar]*best_mom[bar]+0.9383*0.9383);
			double calc_beta = best_mom[bar]/tote;
		      
			double calc_time = best_plen[bar]/calc_beta/30;
		      
			AnaHisto->hd_TimeCalTp[bar]->Fill(calc_time-best_tof[bar]);
		      }
		    //	      cout<<"TKalman "<<__LINE__<<endl;
		  }
	      
	      
		//std::cout<<"Tofp Mom ="<<best_mom[bar]<<" beta ="<<best_beta[bar]<<std::endl;
		// #ifdef DEBUG_DAISUKE
		// 	      std::cout<<"Cal = "<<calc_time<<" Mes="<<best_tof[bar]<<" Calc-Mes="<<calc_time-best_tof[bar]<<std::endl;
		// #endif
	      }
	    else if(best_charge[bar]==4)
	      {
		AnaHisto->h_MassZ2BarBest->Fill(best_mass[bar],bar);
		// 	      if(best_layer[bar]==3)
		// 		AnaHisto->h_MassZ2BarBest3->Fill(best_mass[bar],bar);
		if(best_tofs_coin[bar])
		  {
		    AnaHisto->h_MassZ2BarBestTofs->Fill(best_mass[bar],bar);
		    // 		  if(best_layer[bar]==3)
		    // 		    AnaHisto->h_MassZ2BarBestTofs3->Fill(best_mass[bar],bar);
		    //	      cout<<"TKalman "<<__LINE__<<endl;
		  }

	      }

	    AnaHisto->hd_PoQ_E->Fill(best_mom[bar]/sqrt(best_charge[bar]),best_dE[bar]);
	    // 	  if(best_dE[bar]>40)
	    // 	    std::cout<<"alpha come"<<std::endl;
	    // 	  if(best_layer[bar]==3)
	    // 	    AnaHisto->hd_PoQ_E3->Fill(best_mom[bar]/sqrt(best_charge[bar]),best_dE[bar]);
	    if(best_tofs_coin[bar])
	      {
		AnaHisto->hd_PoQ_ETofs->Fill(best_mom[bar]/sqrt(best_charge[bar]),best_dE[bar]);
		// 	      if(best_layer[bar]==3)
		// 		AnaHisto->hd_PoQ_ETofs3->Fill(best_mom[bar]/sqrt(best_charge[bar]),best_dE[bar]);
	      }

	  }
	//	      cout<<"TKalman "<<__LINE__<<endl;
	if(best_mom_ini[bar]>0)
	  {
	    if(best_charge_ini[bar]==1)
	      {
		AnaHisto->h_MassZ1BarBest_ini->Fill(best_mass_ini[bar],bar);
		// 	      if(best_layer_ini[bar]==3)
		// 		AnaHisto->h_MassZ1BarBest_ini3->Fill(best_mass_ini[bar],bar);
		AnaHisto->hd_BetaMom_Tp_ini_all->Fill(best_mom_ini[bar],best_beta_ini[bar]);
		if(best_tofs_coin_ini[bar])
		  {
		    AnaHisto->h_MassZ1BarBest_iniTofs->Fill(best_mass_ini[bar],bar);
		    // 		  if(best_layer_ini[bar]==3)
		    // 		    AnaHisto->h_MassZ1BarBest_iniTofs3->Fill(best_mass_ini[bar],bar);

		    AnaHisto->hd_BetaMom_Tp_ini_allTofs->Fill(best_mom_ini[bar],best_beta_ini[bar]);
		  }
	      }
	    else if(best_charge_ini[bar]==4)
	      {
		AnaHisto->h_MassZ2BarBest_ini->Fill(best_mass_ini[bar],bar);
		// 	      if(best_layer_ini[bar]==3)
		// 		AnaHisto->h_MassZ2BarBest_ini3->Fill(best_mass_ini[bar],bar);
		if(best_tofs_coin_ini[bar])
		  {
		    AnaHisto->h_MassZ2BarBest_iniTofs->Fill(best_mass_ini[bar],bar);
		    //		  cout<<"TKalman "<<__LINE__<<endl;
		    // 		  if(best_layer_ini[bar]==3)
		    // 		    AnaHisto->h_MassZ2BarBest_iniTofs3->Fill(best_mass_ini[bar],bar);
		  }
	      }

	    AnaHisto->hd_PoQ_E_ini->Fill(best_mom_ini[bar]/sqrt(best_charge_ini[bar]),best_dE_ini[bar]);
	    // 	  if(best_layer_ini[bar]==3)
	    // 	    AnaHisto->hd_PoQ_E_ini3->Fill(best_mom_ini[bar]/sqrt(best_charge_ini[bar]),best_dE_ini[bar]);
	    //	  std::cout<<"best mass "<<bar<<" "<<best_mass[bar]<<endl;
	    if(best_tofs_coin_ini[bar])
	      {
		AnaHisto->hd_PoQ_E_iniTofs->Fill(best_mom_ini[bar]/sqrt(best_charge_ini[bar]),best_dE_ini[bar]);
		// 	      if(best_layer_ini[bar]==3)
		// 		AnaHisto->hd_PoQ_E_iniTofs3->Fill(best_mom_ini[bar]/sqrt(best_charge_ini[bar]),best_dE_ini[bar]);
	      }
	  }


      

	//	      cout<<"TKalman "<<__LINE__<<endl;
	if(best_mom_pi[bar]>0)
	  {
	    AnaHisto->h_MassZ1BarBest_pi->Fill(best_mass_pi[bar],bar);
	    // 	  if(best_layer_pi[bar]==3)
	    // 	    AnaHisto->h_MassZ1BarBest_pi3->Fill(best_mass_pi[bar],bar);
	    AnaHisto->hd_MassPv[1]->Fill(best_mass_pi[bar],1-best_pv_pi[bar]);
	    //	  std::cout<<"best mass "<<bar<<" "<<best_mass[bar]<<endl;
	    //	  AnaHisto->hd_BetaMom_Bg[bar]->Fill(best_mom_pi[bar],best_beta_pi[bar]);
	    AnaHisto->hd_BetaMom_Bg_all->Fill(best_mom_pi[bar],best_beta_pi[bar]);
	    //	      cout<<"TKalman "<<__LINE__<<endl;
	    if(best_tofs_coin_pi[bar])
	      {
		AnaHisto->h_MassZ1BarBest_piTofs->Fill(best_mass_pi[bar],bar);
		// 	      if(best_layer_pi[bar]==3)
		// 		AnaHisto->h_MassZ1BarBest_piTofs3->Fill(best_mass_pi[bar],bar);
		AnaHisto->hd_BetaMom_Bg_allTofs->Fill(best_mom_pi[bar],best_beta_pi[bar]);
	      }

#ifdef DEBUG_DAISUKE
	    std::cout<<"Bg Mom ="<<best_mom_pi[bar]<<" beta ="<<best_beta_pi[bar]<<std::endl;
#endif
	  }

	if(best_mom_pi_ini[bar]>0)
	  {
	    AnaHisto->h_MassZ1BarBest_ini_pi->Fill(best_mass_ini_pi[bar],bar);
	    AnaHisto->hd_BetaMom_Bg_ini_all->Fill(best_mom_pi_ini[bar],best_beta_pi_ini[bar]);
	    // 	  if(best_layer_pi_ini[bar]==3)
	    // 	    AnaHisto->h_MassZ1BarBest_ini_pi3->Fill(best_mass_ini_pi[bar],bar);
	    //	      cout<<"TKalman "<<__LINE__<<endl;
	    if(best_tofs_coin_pi_ini[bar])


	      {
		AnaHisto->h_MassZ1BarBest_ini_piTofs->Fill(best_mass_ini_pi[bar],bar);
		// 	      if(best_layer_pi_ini[bar]==3)
		// 		AnaHisto->h_MassZ1BarBest_ini_piTofs3->Fill(best_mass_ini_pi[bar],bar);
		AnaHisto->hd_BetaMom_Bg_ini_allTofs->Fill(best_mom_pi_ini[bar],best_beta_pi_ini[bar]);
	      }
	    //	  std::cout<<"best mass "<<bar<<" "<<best_mass[bar]<<endl;
	  }
      }

    return 0;

}
