#include "TCheckFiberTrack.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"
#include "FullRecoEvent.hh"
#include "FullRecoEventZMQ.hh"
#include "KalmanFittedStateOnPlane.h"
#include "KalmanFitterInfo.h"
#include "PlanarMeasurement.h"
#include "ReturnRes.hh"
#include "StateOnPlane.h"

#include "TVector3.h"

#include <numeric>
#include <set>
#include <sstream>
#include <tuple>
#include <list>

using namespace std;
using namespace G4Sol;

template<class Out>
TCheckFiberTrack<Out>::TCheckFiberTrack(const THyphiAttributes& attribut) : TDataProcessInterface<Out>("CheckFiberTrack"), att(attribut)
{

}

template<class Out>
TCheckFiberTrack<Out>::~TCheckFiberTrack()
{ }

template<class Out>
void TCheckFiberTrack<Out>::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }

template<class Out>
ReturnRes::InfoM TCheckFiberTrack<Out>::operator()(FullRecoEvent& RecoEvent, Out* OutTree)
{

  int result_finder = Exec(RecoEvent, OutTree);

  return SoftExit(result_finder);
}

template<class Out>
int TCheckFiberTrack<Out>::Exec(FullRecoEvent& RecoEvent, Out* )
{
  return CheckTrackFinding(RecoEvent);
}

template<class Out>
ReturnRes::InfoM TCheckFiberTrack<Out>::SoftExit(int ) { return ReturnRes::Fine; }

template<class Out>
void TCheckFiberTrack<Out>::SelectHists()
{
  LocalHisto.h_ResidualFiberX = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberX);
  LocalHisto.h_ResidualFiberY = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberY);
  for(int i=0; i<11; ++i)
    for(int j=0; j<2; ++j)
      {
	LocalHisto.h_ResidualFiberX_Angle[i][j] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberX_Angle[i][j]);
	LocalHisto.h_ResidualFiberY_Angle[i][j] = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberY_Angle[i][j]);
      }
  LocalHisto.h_ResidualFiberDzDphi = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberDzDphi);
  LocalHisto.h_ResidualFiberDzDtheta = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberDzDtheta);

}

template<class Out>
int TCheckFiberTrack<Out>::CheckTrackFinding(const FullRecoEvent& RecoEvent)
{

  int ntrack = -1;
  for(auto it_trackInfo : RecoEvent.TrackInfo)
    {
      ntrack++;

#ifdef DEBUG_FIBERTRACK
      att._logger->debug("start #{}", ntrack); //<<std::endl;
                                               // std::cout<<" | Pointsize "<<Vtracks->getNumPoints()<<" Measurements
                                               // "<<Vtracks->getNumPointsWithMeasurement()<<" Rep
                                               // "<<Vtracks->getNumReps()<<std::endl;
                                               // std::vector<std::string> name_det;
#endif

      const int id_track = it_trackInfo.first;
      auto it_ListHits   = RecoEvent.TrackDAF.find(id_track);
      auto it_ListHitsSim = RecoEvent.TrackDAFSim.find(id_track);

#ifdef DEBUG_FIBERTRACK

      auto printW = [](const auto a, const int width) -> std::string {
	std::stringstream ss;
	ss << std::fixed << std::right;
	ss.fill(' ');    // fill space around displayed #
	ss.width(width); // set  width around displayed #
	ss << a;
	return ss.str();
      };

      std::vector<std::stringstream> s1(it_trackInfo.second.size() / 10 + 1);
      std::vector<std::stringstream> s2(it_trackInfo.second.size() / 10 + 1);
      std::vector<std::stringstream> s3(it_trackInfo.second.size() / 10 + 1);
      for(size_t i = 0; i < it_trackInfo.second.size(); ++i)
	{
	  s1[i / 10] << printW(G4Sol::nameLiteralDet.begin()[i], 10) << ", ";
	  s2[i / 10] << printW(i, 10) << ", ";
	  s3[i / 10] << printW(it_trackInfo.second[i].pdg, 10) << ", ";
	}
      for(size_t i = 0; i < s1.size(); ++i)
	{
	  att._logger->debug("Det  :{}", s1[i].str());
	  att._logger->debug("idDet:{}", s2[i].str());
	  att._logger->debug("stat :{}", s3[i].str());
	}
#endif
	  
      bool id_FirstFiberX = it_ListHitsSim->second[G4Sol::FiberD3_x].size() > 0 ;
      bool id_FirstFiberU = it_ListHitsSim->second[G4Sol::FiberD3_u].size() > 0 ;
      bool id_FirstFiberV = it_ListHitsSim->second[G4Sol::FiberD3_v].size() > 0 ;

      std::list<int> id_dets = {G4Sol::MiniFiberD1_x1,  G4Sol::MiniFiberD1_u1,  G4Sol::MiniFiberD1_v1,  G4Sol::MiniFiberD1_x2,  G4Sol::MiniFiberD1_u2,  G4Sol::MiniFiberD1_v2, G4Sol::PSFE, G4Sol::PSCE, G4Sol::PSBE, G4Sol::FiberD4_x, G4Sol::FiberD4_u, G4Sol::FiberD4_v, G4Sol::FiberD5_x, G4Sol::FiberD5_u, G4Sol::FiberD5_v};

      int id_FirstFiber = -1;
      if(id_FirstFiberX)
	{
	  id_FirstFiber = G4Sol::FiberD3_x;
	  id_dets.push_front(G4Sol::FiberD3_v);
	  id_dets.push_front(G4Sol::FiberD3_u);
	}
      else if(id_FirstFiberU)
	{
	  id_FirstFiber = G4Sol::FiberD3_u;
	  id_dets.push_front(G4Sol::FiberD3_v);
	}
      else if(id_FirstFiberV)
	id_FirstFiber = G4Sol::FiberD3_v;
      
      //att._logger->debug("first fiber : {} {} {} | {}",id_FirstFiberX,id_FirstFiberU,id_FirstFiberV, id_FirstFiber); 
      
      if(id_FirstFiber == -1)
	{
#ifdef DEBUG_FIBERTRACK
	  att._logger->debug("E> should have fiber D3 hit ! ");

#endif

	  continue;
	}

      std::array<double,6> LineParams = {it_ListHitsSim->second[id_FirstFiber][0].hitX,
					 it_ListHitsSim->second[id_FirstFiber][0].hitY,
					 it_ListHitsSim->second[id_FirstFiber][0].hitZ,
					 it_ListHitsSim->second[id_FirstFiber][0].momX,
					 it_ListHitsSim->second[id_FirstFiber][0].momY,
					 it_ListHitsSim->second[id_FirstFiber][0].momZ};
      
      auto PDG_particle = TDatabasePDG::Instance()->GetParticle(it_ListHitsSim->second[id_FirstFiber][0].pdg);
      const double charge = PDG_particle->Charge() / 3.;

      auto f_line = [](const std::array<double, 6>& params, double posZ, std::array<double, 3>& hitExtrap)
      {
	if(std::abs(params[2]-posZ)<1e-3)
	  {
	    hitExtrap[0] = params[0];
	    hitExtrap[1] = params[1];
	    hitExtrap[2] = params[2];
	  }
	else
	  {
	    hitExtrap[0] = params[0] + ( params[3] / params[5] ) * ( posZ - params[2]);
	    hitExtrap[1] = params[1] + ( params[4] / params[5] ) * ( posZ - params[2]);
	    hitExtrap[2] = posZ;
	  }
      };

      double theta = TMath::ATan2(TMath::Sqrt(TMath::Sq(LineParams[3])+TMath::Sq(LineParams[4])),LineParams[5])*TMath::RadToDeg();

      std::array<double, 3> hitExp = {0.,0.,0.};
      TVector3 beforeV(LineParams[3],LineParams[4],LineParams[5]), afterV;
      for(int id_det : id_dets)
	{
	  //att._logger->debug(" det: {} : size {} / condition : {}",G4Sol::nameLiteralDet.begin()[id_det], it_ListHitsSim->second[id_det].size(), id_det < id_FirstFiber);
	  
	  // if(id_det < id_FirstFiber)
	  //   continue;

	  if(it_ListHitsSim->second[id_det].size()>0)
	    {
	      
	      f_line(LineParams, it_ListHitsSim->second[id_det][0].hitZ, hitExp);
	      TString nameT(G4Sol::nameLiteralDet.begin()[id_det]);
	      //att._logger->debug("    -> hitZ: {}, hit xy {},{} : Exp xy {},{}", it_ListHitsSim->second[id_det][0].hitZ, it_ListHitsSim->second[id_det][0].hitX, it_ListHitsSim->second[id_det][0].hitY, hitExp[0], hitExp[1]);

	      const TString nameTemp[] = {"FiberD3_u","FiberD3_v","MiniFiberD1_x","MiniFiberD1_u","MiniFiberD1_v",
					  "MiniFiberD2_x","MiniFiberD2_u","MiniFiberD2_v","PSCE","PSFE","PSCE","PSBE"};

	      int temp_i = 0;
	      for(int i=0;i<11;++i)
		if(nameT == nameTemp[i])
		  {
		    temp_i = i;
		    break;
		  }

	      int temp_j = charge < 0 ? 0 : 1;
	      
	      LocalHisto.h_ResidualFiberX->Fill(nameT.Data(), charge*(hitExp[0] - it_ListHitsSim->second[id_det][0].hitX),1.);
	      LocalHisto.h_ResidualFiberY->Fill(nameT.Data(), charge*(hitExp[1] - it_ListHitsSim->second[id_det][0].hitY),1.);

	      LocalHisto.h_ResidualFiberX_Angle[temp_i][temp_j]->Fill(theta, (hitExp[0] - it_ListHitsSim->second[id_det][0].hitX),1.);
	      LocalHisto.h_ResidualFiberY_Angle[temp_i][temp_j]->Fill(theta, (hitExp[1] - it_ListHitsSim->second[id_det][0].hitY),1.);

	      afterV.SetXYZ(it_ListHitsSim->second[id_det][0].momX,it_ListHitsSim->second[id_det][0].momY,it_ListHitsSim->second[id_det][0].momZ);
	      double Dphi = (afterV.Phi() - beforeV.Phi())*TMath::RadToDeg();
	      double Dtheta = (afterV.Theta() - beforeV.Theta())*TMath::RadToDeg();
	      LocalHisto.h_ResidualFiberDzDphi->Fill(it_ListHitsSim->second[id_det][0].hitZ-LineParams[2], Dphi);
	      LocalHisto.h_ResidualFiberDzDtheta->Fill(it_ListHitsSim->second[id_det][0].hitZ-LineParams[2], Dtheta);

	    }
	}
    }

  return 0;
}

template class TCheckFiberTrack<MCAnaEventG4Sol>;
template class TCheckFiberTrack<Ana_WasaEvent>;

