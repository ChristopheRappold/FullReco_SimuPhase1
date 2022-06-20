#include "TBayesFinder.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "ReturnRes.hh"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TGeoMatrix.h"
#include "TGeoTube.h"

#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "TObjArray.h"

#include <iomanip>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <memory>
#include <tuple>

using namespace std;
using namespace G4Sol;
using namespace BayesFind;

#define DEBUG_BAYES

template<class Out>
TBayesFinder<Out>::TBayesFinder(const THyphiAttributes& attribut):TDataProcessInterface<Out>("bayes_finder"),att(attribut),ME( 17, std::vector<TGeoNodeMatrix*>()),LayerGeo(17, std::vector<DataLayer>())
{
  for(auto it : *gGeoManager->GetListOfVolumes() )
    {
      TString name(it->GetName());
      
      if(name.Contains("MD"))
	{
	  att._logger->debug("{} {}",name.Data(), fmt::ptr(it));
	  TGeoTubeSeg* tube = dynamic_cast<TGeoTubeSeg*>( dynamic_cast<TGeoVolume*>(it)->GetShape());
	  radiusCDC.emplace_back(0.5*(tube->GetRmax()+tube->GetRmin()));
	  att._logger->debug("{} {}",name.Data(), radiusCDC.back());
	}
    }


  

  auto N00 = gGeoManager->GetListOfNodes();
  TGeoNodeMatrix* WASA_1 = (TGeoNodeMatrix*)N00->At(0);
  TString nameWASA = WASA_1->GetName();
  if(nameWASA.Contains("WASA"))
    {
      std::cout<<"WASA:\n";
      WASA_1->GetMatrix()->Print();
      auto N01 = WASA_1->GetNodes();

      TGeoNodeMatrix* MFLD_1 = (TGeoNodeMatrix*)N01->At(0);
      std::cout<<"MFLD:\n";
      MFLD_1->GetMatrix()->Print();
  
      auto N02 = MFLD_1->GetNodes();

      //N02->Print();
      TGeoNodeMatrix* INNER_1 = (TGeoNodeMatrix*)N02->At(0);
      std::cout<<"INNER:\n";
      INNER_1->GetMatrix()->Print();
      TGeoHMatrix tempMasterInner(*INNER_1->GetMatrix()), tempMasterMag(*MFLD_1->GetMatrix()); 
      TGeoHMatrix ToMaster = tempMasterInner*tempMasterMag;
      std::cout<<"ToINNER:\n:";
      ToMaster.Print();
  
      auto N03 = INNER_1->GetNodes();
      //N03->Print();

      TGeoNodeMatrix* MD[17]; 
      unsigned int Nsize[17];
      std::vector<TGeoHMatrix> MatMD;
      for(int i=0;i<17;++i)
	{
	  MD[i] = nullptr;
	  MD[i] = (TGeoNodeMatrix*)N03->At(1+i);
	  Nsize[i] = MD[i]->GetNodes()->GetEntries()-2;

	  std::cout<<i<<"MD \n";
	  //MD[i]->GetMatrix()->Print();
	  TGeoHMatrix temp1(*MD[i]->GetMatrix());
	  
	  MatMD.emplace_back(temp1*ToMaster);

	  MatMD.back().Print();
	}

      //ME.resize( 17, std::vector<TGeoNodeMatrix*>());
      //LayerGeo.resize(17, std::vector<DataLayer>());
  
      for( auto i = 0 ; i < ME.size(); ++i)
	{
	  att._logger->debug("layer:{} nb channels: {}", i, Nsize[i]);
	  ME[i] = std::vector<TGeoNodeMatrix*>(Nsize[i],nullptr);
	  LayerGeo[i].resize(Nsize[i]);
	  auto N04 = MD[i]->GetNodes();
	  for( auto j = 0; j<Nsize[i]; ++j)
	    {
	      //att._logger->debug("ch#{} / {}",j, Nsize[i]);
	      ME[i][j] = (TGeoNodeMatrix*)N04->At(j);

	      auto TempTube = dynamic_cast<TGeoTube*>(ME[i][j]->GetVolume()->GetShape());

	      Double_t TempCenter[3] = {0.,0.,0.}; 
	      Double_t TempMinC[3] = {0.,0.,-TempTube->GetDz()};
	      Double_t TempMaxC[3] = {0.,0.,TempTube->GetDz()};

	      Double_t CenterInMaster[3] = {0.,0.,0.};
	      Double_t MinInMaster[3] = {0.,0.,0.};
	      Double_t MaxInMaster[3] = {0.,0.,0.};
	      ME[i][j]->LocalToMaster(TempCenter, CenterInMaster);
	      ME[i][j]->LocalToMaster(TempMinC, MinInMaster);
	      ME[i][j]->LocalToMaster(TempMaxC, MaxInMaster);

	      Double_t CenterInLab[3] = {0.,0.,0.};
	      Double_t MinInLab[3] = {0.,0.,0.};
	      Double_t MaxInLab[3] = {0.,0.,0.};

	      MatMD[i].LocalToMaster(CenterInMaster, CenterInLab);
	      MatMD[i].LocalToMaster(MinInMaster, MinInLab);
	      MatMD[i].LocalToMaster(MaxInMaster, MaxInLab);
	  
	      LayerGeo[i][j].cenX = CenterInLab[0];
	      LayerGeo[i][j].cenY = CenterInLab[1];
	      LayerGeo[i][j].cenZ = CenterInLab[2];
	      
	      LayerGeo[i][j].minX = MinInLab[0];
	      LayerGeo[i][j].minY = MinInLab[1];
	      LayerGeo[i][j].minZ = MinInLab[2];
	      
	      LayerGeo[i][j].maxX = MaxInLab[0];
	      LayerGeo[i][j].maxY = MaxInLab[1];
	      LayerGeo[i][j].maxZ = MaxInLab[2];

	      LayerGeo[i][j].radius = TempTube->GetRmax();
	    }
	}
    }
    
  
}
template<class Out>
TBayesFinder<Out>::~TBayesFinder()
{

}

template<class Out>
void TBayesFinder<Out>::InitMT() {att._logger->error("E> Not supposed to be multithreaded !"); }

template<class Out>
ReturnRes::InfoM TBayesFinder<Out>::operator() (FullRecoEvent& RecoEvent,Out* OutTree)
{

  int result_finder = Exec(RecoEvent,OutTree);

  return SoftExit(result_finder);
}

template<class Out>
int TBayesFinder<Out>::Exec(FullRecoEvent& RecoEvent,Out* OutTree)
{
  return FinderTrack(RecoEvent);
  
}

template<class Out>
ReturnRes::InfoM TBayesFinder<Out>::SoftExit(int result_full)
{
  return ReturnRes::Fine;
}

template<class Out>
void TBayesFinder<Out>::SelectHists()
{
    LocalHisto.h_xy          = AnaHisto->CloneAndRegister(AnaHisto->h_xy);
    LocalHisto.h_PxPy        = AnaHisto->CloneAndRegister(AnaHisto->h_PxPy);
    LocalHisto.h_xy_extrap   = AnaHisto->CloneAndRegister(AnaHisto->h_xy_extrap);
    LocalHisto.h_PxPy_extrap = AnaHisto->CloneAndRegister(AnaHisto->h_PxPy_extrap);
    LocalHisto.h_TrackFindingStat = AnaHisto->CloneAndRegister(AnaHisto->h_TrackFindingStat);
  for(int i=0;i<3;++i)
    LocalHisto.h_SolenoidGeo[i]         = AnaHisto->CloneAndRegister(AnaHisto->h_SolenoidGeo[i]);
}


template<class Out>
int TBayesFinder<Out>::FinderTrack(FullRecoEvent& RecoEvent)
{
  int index_ell = 0;
  for(double radius : radiusCDC)
    if(index_ell<17)
      {
	att._logger->error("GEO: {} {}",fmt::ptr(AnaHisto->geoSolenoid[index_ell]), radius);
	AnaHisto->geoSolenoid[index_ell] = new TEllipse(0,0,radius);
	AnaHisto->geoSolenoid[index_ell]->SetFillStyle(2);
	++index_ell;
      }

  
#ifdef DEBUG_BAYES
  cout<<"Finder DAF Hit :"<<endl;
  auto printW = [](const auto a, const int width) -> std::string {
		  std::stringstream ss;
		  ss << std::fixed << std::right;
		  ss.fill(' ');        // fill space around displayed #
		  ss.width(width);     // set  width around displayed #
		  ss << a;
		  return ss.str();
		};
  for(auto track : RecoEvent.TrackDAF)
    {

      std::cout<<"TrackID #"<<track.first<<" hit_id :\n";
      std::vector<std::stringstream> s1(track.second.size()/20+1);
      std::vector<std::stringstream> s2(track.second.size()/20+1);
      for(size_t i = 0; i<track.second.size(); ++i)
	{
	  s1[i/20] << printW(G4Sol::nameLiteralDet.begin()[i],6)<<", ";
	  s2[i/20] << printW(track.second[i],6)<<", ";
	}
      for(size_t i = 0; i<s1.size();++i)
	{
	  std::cout<<"idDet:"<<s1[i].str()<<"\n";
	  std::cout<<"stat :"<<s2[i].str()<<"\n";
	}
      // for(auto id_hit : track.second)
      // 	std::cout<<" "<<id_hit<<", ";
      // std::cout<<"] "<<std::endl;
    }
  for(auto track : RecoEvent.TrackInfo)
    {
      std::cout<<"TrackID #"<<track.first<<" PID :\n";
      std::vector<std::stringstream> s1(track.second.size()/20+1);
      std::vector<std::stringstream> s2(track.second.size()/20+1);
      for(size_t i = 0; i<track.second.size(); ++i)
	{
	  s1[i/20] << printW(G4Sol::nameLiteralDet.begin()[i],6)<<", ";
	  s2[i/20] << printW(track.second[i].pdg,6)<<", ";
	}
      for(size_t i = 0; i<s1.size();++i)
	{
	  std::cout<<"idDet:"<<s1[i].str()<<"\n";
	  std::cout<<"PDG  :"<<s2[i].str()<<"\n";
	}
      // for(auto id_hit : track.second)
      // 	std::cout<<" :"<<id_hit.pdg<<", ";
      // std::cout<<"] "<<std::endl;
    }


  std::cout<<"Start Finder:\n";
#endif

  std::vector< std::vector<std::tuple<SimHit,int> > > AllHits(G4Sol::SIZEOF_G4SOLDETTYPE);
  for(auto it_track : RecoEvent.TrackDAFSim)
    {
      for(size_t id_det = 0 ; id_det < it_track.second.size() ; ++id_det )
	{
	  if(it_track.second[id_det][0].layerID>=0)
	    AllHits[id_det].emplace_back(std::make_tuple(it_track.second[id_det][0], it_track.first) );
	}
    }

  auto DrawMG = [](const auto& ListHits, const auto& LayerG, TH2F* h, TRandom3& rand) {
    
    for(SolDetIter i = G4Sol::MG01 ; i != G4Sol::MG17;++i)
      {
	for(auto it_hit : ListHits[i.pos])
	  {
	    SimHit hit( std::get<0>(it_hit));
	    double xH = hit.hitX;
	    double yH = hit.hitY;
	    int WireID = hit.layerID;
	    int LayID = i.pos - G4Sol::MG01;
	    //LayID -= initLayID;
	    for(int ran = 0;ran<1000;++ran)
	      {
		double temp_x,temp_y;
		rand.Circle(temp_x,temp_y,LayerG[LayID][WireID].radius);
		h->Fill(xH+temp_x,yH+temp_y);
	      }
	  }
      }
  };

  DrawMG(AllHits,LayerGeo,LocalHisto.h_SolenoidGeo[0],rand);

#ifdef DEBUG_BAYES
  auto printSimHit = [](const SimHit& hit) {
		       std::cout<<"     layerID: "<<hit.layerID<<" pdg:"<<hit.pdg<<"\n";
		       std::cout<<"     hit: ["<<hit.hitX<<", "<<hit.hitY<<", "<<hit.hitZ<<"] \n";
		     };
  for(size_t id_det = 0; id_det < AllHits.size() ; ++id_det)
    if(AllHits[id_det].size()>0)
      {
	std::cout<<"detector: #"<<id_det<<" "<<G4Sol::nameLiteralDet.begin()[id_det]<<" "<<AllHits[id_det].size()<<"\n";
	for(auto it_hit : AllHits[id_det])
	  {
	    printSimHit(std::get<0>(it_hit));
	    std::cout<<"     track id:"<<std::get<1>(it_hit)<<"\n";
	    std::cout<<"     ----\n";
	  }
      }
#endif

  std::vector< std::vector<int> > FoundTrack;
  
  for(auto it_det : RecoEvent.TrackInfo)
    {
      int id_track = it_det.first;

      auto ListValidHits = RecoEvent.TrackDAF.find(id_track)->second;

      if(ListValidHits[G4Sol::PSCE]<0 && ListValidHits[G4Sol::PSFE]<0)
	continue;

      std::vector<int> SetHits (G4Sol::SIZEOF_G4SOLDETTYPE,-1);
      
      auto ListHits = RecoEvent.TrackDAFSim.find(id_track)->second[0];

#ifdef DEBUG_BAYES
      auto printValid = [](auto ListValidHits) {
			  for(auto id_hit : ListValidHits)
			    std::cout<<" :"<< id_hit <<", ";
			};
      printValid(ListValidHits);
      cout<<"\n Just CDC:\n";
      auto printMG = [](auto List) {
		       for(SolDetIter i = G4Sol::MG01 ; i != G4Sol::PSCE;++i)
			 std::cout<<List[i.pos]<<", ";
		     };
      printMG(ListValidHits);
#endif
           
      const double charge = -1.;
      double ch =charge*att.Field_Strength/3.10715497;
      
      InfoPar info;
      
      double init_x = 0, init_y = 0, init_momX = 0, init_momY = 0;

      G4Sol::SolDet LastFrontWall = att.Wasa_Side == 0 ? G4Sol::PSFE : G4Sol::PSBE ;
      const std::string nameFinalHit[2] = {"PSCE","PSBFE"};

      int idFinal = 0;
      if(ListValidHits[G4Sol::PSCE]>=0)
	{
	  init_x = ListHits[G4Sol::PSCE].hitX;
	  init_y = ListHits[G4Sol::PSCE].hitY;
	  init_momX = ListHits[G4Sol::PSCE].momX;
	  init_momY = ListHits[G4Sol::PSCE].momY;
	  info = it_det.second[G4Sol::PSCE];
	  idFinal = 0;
	}
      else if(ListValidHits[LastFrontWall]>=0)
	{
	  init_x = ListHits[LastFrontWall].hitX;
	  init_y = ListHits[LastFrontWall].hitY;
	  init_momX = ListHits[LastFrontWall].momX;
	  init_momY = ListHits[LastFrontWall].momY;
	  info = it_det.second[LastFrontWall];
	  idFinal = 1;
	}

#ifdef DEBUG_BAYES
      auto printInfo = [](const InfoPar& a,double ch)
		       {
			 std::cout<<"pdg: "<<a.pdg<<" ["<<a.momX<<","<<a.momY<<","<<a.momZ<<"] t:"<<a.time<<"\n";
			 std::cout<<"Pt:"<<TMath::Sqrt(a.momX*a.momX+a.momY*a.momY)<<" R:"<<TMath::Sqrt(a.momX*a.momX+a.momY*a.momY)/ch<<"\n";
		       };
	  std::cout<<"track: "<<id_track<<" ";
	  printInfo(info,ch);
#endif

      LocalHisto.h_xy->Fill(init_x,init_y);
      LocalHisto.h_PxPy->Fill(init_momX,init_momY);
      
      double new_x[2]={0.,0.}, new_y[2]={0.,0.}, new_momX[2]={0.,0.}, new_momY[2]={0.,0.};
      
      
#ifdef DEBUG_BAYES
      auto printFixed = [](const double a, const int decDigits, const int width) -> std::string {
			  std::stringstream ss;
			  ss << std::fixed << std::right;
			  ss.fill(' ');        // fill space around displayed #
			  ss.width(width);     // set  width around displayed #
			  ss.precision(decDigits); // set # places after decimal
			  ss << a;
			  return ss.str();
			};
      std::cout<<"Init: {"<<printFixed(init_x,5,10)<<"} \n      {"<<printFixed(init_y,5,10)<<"} \n      {"<<printFixed(init_momX,5,10)<<"} \n      {"<<printFixed(init_momY,5,10)<<"} \n";
#endif
      // int n_ValidCDC = 0;
      // for(int idCDC = G4Sol::MG17; idCDC >= G4Sol::MG01; --idCDC)
      // 	if(ListValidHits[idCDC]>=0)
      // 	  ++n_ValidCDC;
      
      for(int idCDC = G4Sol::MG17; idCDC >= G4Sol::MG01; --idCDC)
	//if(ListValidHits[idCDC]>=0)
	if(AllHits[idCDC].size()>0)
	  {

	    double best_x = 0., best_y = 0., best_momX = 0., best_momY = 0.;
	    bool best_found =  false;
	    int id_hit = 0;
	    for(auto it_hit_track : AllHits[idCDC])
	      {
		SimHit it_hit = std::get<0>(it_hit_track);

#ifdef DEBUG_BAYES
		std::cout<<"MG "<<idCDC-G4Sol::MG01+1<<" "<<radiusCDC[idCDC-G4Sol::MG01]<<" hit#"<<id_hit<<" \n";
#endif
		LocalHisto.h_xy->Fill(it_hit.hitX, it_hit.hitY);
		LocalHisto.h_PxPy->Fill(it_hit.momX, it_hit.momY);
	    
		double Px = init_momX; //= ListHits[idCDC].momX;
		double Py = init_momY; //= ListHits[idCDC].momY;
		double x  = init_x*1.e-2; //= ListHits[idCDC].hitX*1e-2;
		double y  = init_y*1.e-2; //= ListHits[idCDC].hitY*1e-2;
	    
		double Rn = radiusCDC[idCDC-G4Sol::MG01]*1e-2;

		double a = Px*Px/ch/ch + Py*Py/ch/ch + x*Py/ch - y*Px/ch; //Px/ch+Py/ch; 
		double b = 2.*(x*Px+y*Py)/ch;
		double c = x*x + y*y - Rn*Rn;

		double delta = b*b - 4.*a*c;
		double r1 = 0.5*(-b + TMath::Sqrt(delta)) / a;
		double r2 = 0.5*(-b - TMath::Sqrt(delta)) / a;


		//double dphi2_1 =-(-ch*Px*x-ch*Py*y+ch*TMath::Sqrt(ch*Px*y*y*y+(-Px*Px-ch*Py*x)*y*y+(-ch*Px*Rn*Rn+2*Px*Py*x+ch*Px*x*x)*y-ch*Py*x*x*x-Py*Py*x*x+ch*Py*Rn*Rn*x+(Px*Px+Py*Py)*Rn*Rn))/(ch*Px*y-ch*Py*x-Py*Py-Px*Px);
		//double dphi2_2 =(ch*Px*x+ch*Py*y+ch*TMath::Sqrt(ch*Px*y*y*y+(-Px*Px-ch*Py*x)*y*y+(-ch*Px*Rn*Rn+2*Px*Py*x+ch*Px*x*x)*y-ch*Py*x*x*x-Py*Py*x*x+ch*Py*Rn*Rn*x+(Px*Px+Py*Py)*Rn*Rn))/(ch*Px*y-ch*Py*x-Py*Py-Px*Px);

	    
#ifdef DEBUG_BAYES
		std::cout<<" a dPhi^2 + b dPhi + c = 0 \n"<<a<<" "<<b<<" "<<c<<"\n";
		cout<<"delta:"<<delta<<" "<<r1<<" "<<r2<<"\n";
		//cout<<"dphi2: "<<dphi2_1<<" "<<dphi2_2<<"\n";
#endif
		double dPhi = TMath::Abs(r1) < TMath::Abs(r2) ? r1 : r2;

		auto f_applyTracking = [](double dPhi, std::tuple<double, double, double, double> vec, double ch) {
				     
					 double x,y,Px,Py;
					 std::tie(x,y,Px,Py) = vec;
					 return std::make_tuple(x+(1.-TMath::Cos(dPhi))*Py/ch + TMath::Sin(dPhi)*Px/ch,
								y-(1.-TMath::Cos(dPhi))*Px/ch + TMath::Sin(dPhi)*Py/ch,
								TMath::Cos(dPhi)*Px+TMath::Sin(dPhi)*Py,
								TMath::Cos(dPhi)*Py-TMath::Sin(dPhi)*Px);
				       };


		std::tie(new_x[0], new_y[0], new_momX[0], new_momY[0]) = f_applyTracking(r1,std::make_tuple(x,y,Px,Py),ch);
		std::tie(new_x[1], new_y[1], new_momX[1], new_momY[1]) = f_applyTracking(r2,std::make_tuple(x,y,Px,Py),ch);

		double temp_x = 0., temp_y = 0., temp_momX = 0., temp_momY = 0.;
		std::tie(temp_x, temp_y, temp_momX, temp_momY) = f_applyTracking(dPhi,std::make_tuple(x,y,Px,Py),ch);

		// init_x *= 1e2;
		// init_y *= 1e2;
	    
#ifdef DEBUG_BAYES
		std::cout<<"new X :"<<printFixed(it_hit.hitX,5,10)<<" {"<<printFixed(new_x[0]*1e2,5,10)<<", "<<printFixed(new_x[1]*1e2,5,10)<<"} \n";
		std::cout<<"new Y :"<<printFixed(it_hit.hitY,5,10)<<" {"<<printFixed(new_y[0]*1e2,5,10)<<", "<<printFixed(new_y[1]*1e2,5,10)<<"} \n";
		std::cout<<"new Px:"<<printFixed(it_hit.momX,5,10)<<" {"<<printFixed(new_momX[0] ,5,10)<<", "<<printFixed(new_momX[1] ,5,10)<<"} \n";
		std::cout<<"new Py:"<<printFixed(it_hit.momY,5,10)<<" {"<<printFixed(new_momY[0] ,5,10)<<", "<<printFixed(new_momY[1] ,5,10)<<"} \n";

		std::cout<<"diff :"<<printFixed(it_hit.hitX-temp_x*1e2,5,10)<<" "<<printFixed(it_hit.hitY-temp_y*1e2,5,10)<<"\n";
	    
		std::cout<<"Pt:"<<TMath::Sqrt(it_hit.momX*it_hit.momX+it_hit.momY*it_hit.momY)
			 <<"R: "<<TMath::Sqrt(it_hit.momX*it_hit.momX+it_hit.momY*it_hit.momY)/ch<<"\n";
#endif
		if(TMath::Sqrt( (it_hit.hitX-temp_x*1e2)*(it_hit.hitX-temp_x*1e2) + (it_hit.hitY-temp_y*1e2)*(it_hit.hitY-temp_y*1e2)) < 1. /*cm*/ )
		  {
#ifdef DEBUG_BAYES
		    std::cout<<" close hit !\n";
#endif
		    best_x = it_hit.hitX;
		    best_y = it_hit.hitY;
		    best_momX = temp_momX;
		    best_momY = temp_momY;
		    best_found = true;
		    SetHits[idCDC] = id_hit;
		  }
		++id_hit;
	      }

	    if(best_found==true)
	      {
		init_x = best_x;
		init_y = best_y;
		init_momX = best_momX;
		init_momY = best_momY;

		LocalHisto.h_xy_extrap->Fill(init_x, init_y);
		LocalHisto.h_PxPy_extrap->Fill(init_momX, init_momY);
	      }
	  }

      auto checkSameTrack = [](const auto ListValidHits, const auto SetHits) {
			      bool same = true;
			      int nb_diff = 0;
			      for(SolDetIter i = G4Sol::MG01; i != G4Sol::PSCE;++i)
				if(ListValidHits[i.pos]!=SetHits[i.pos])
				  {
				    same = false;
				    ++nb_diff;
				  }
			      if(same)
				nb_diff = -1;
			      
			      return nb_diff;
			    };
      
      int sameTrack = checkSameTrack(ListValidHits,SetHits);
      std::string nameX ( TDatabasePDG::Instance()->GetParticle(info.pdg)->GetName());
      nameX += nameFinalHit[idFinal];
      LocalHisto.h_TrackFindingStat->Fill(nameX.c_str(),sameTrack,1);
      LocalHisto.h_TrackFindingStat->Fill(nameX.c_str(),0.,1);
      
#ifdef DEBUG_BAYES
      std::cout<<"results:"<<sameTrack <<" \n";
      std::cout<<"Sim:   "; printMG(ListValidHits); std::cout<<"\n";
      std::cout<<"Found: "; printMG(SetHits); std::cout<<"\n";
      std::cout<<"--------\n";
#endif
    }

  



  
  return 0;
}

template class TBayesFinder<MCAnaEventG4Sol>;
