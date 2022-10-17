#include "TCheckFiberTrack.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"
#include "FullRecoEvent.hh"
#include "KalmanFittedStateOnPlane.h"
#include "KalmanFitterInfo.h"
#include "PlanarMeasurement.h"
#include "ReturnRes.hh"
#include "StateOnPlane.h"

#include <numeric>
#include <set>
#include <sstream>
#include <tuple>

using namespace std;
using namespace G4Sol;

template<class Out>
TCheckFiberTrack<Out>::TCheckFiberTrack(const THyphiAttributes& attribut) : TDataProcessInterface<Out>("CheckFiberTrack"), att(attribut)
{

}

template<class Out>
TCheckFiberTrack<Out>::~TCheckFiberTrack()
{

}

template<class Out>
void TCheckFiberTrack<Out>::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }

template<class Out>
ReturnRes::InfoM TCheckFiberTrack<Out>::operator()(FullRecoEvent& RecoEvent, Out* OutTree)
{

  int result_finder = Exec(RecoEvent, OutTree);

  return SoftExit(result_finder);
}

template<class Out>
int TCheckFiberTrack<Out>::Exec(FullRecoEvent& RecoEvent, Out* OutTree)
{
  return CheckTrackFinding(RecoEvent);
}

template<class Out>
ReturnRes::InfoM TCheckFiberTrack<Out>::SoftExit(int result_full) { return ReturnRes::Fine; }

template<class Out>
void TCheckFiberTrack<Out>::SelectHists()
{
  LocalHisto.h_ResidualFiberX = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberX);
  LocalHisto.h_ResidualFiberY = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualFiberY);

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


      bool id_FirstFiberX = it_ListHitsSim->second[G4Sol::FiberD3_x].size() > 0 ;
      bool id_FirstFiberU = it_ListHitsSim->second[G4Sol::FiberD3_u].size() > 0 ;
      bool id_FirstFiberV = it_ListHitsSim->second[G4Sol::FiberD3_v].size() > 0 ;

      int id_FirstFiber = -1;
      if(id_FirstFiberX)
	id_FirstFiber = G4Sol::FiberD3_x;
      else if(id_FirstFiberU)
	id_FirstFiber = G4Sol::FiberD3_u;
      else if(id_FirstFiberV)
	id_FirstFiber = G4Sol::FiberD3_v;

      if(id_FirstFiber == -1)
	{
#ifdef DEBUG_FIBERTRACK
	  att._logger->debug("E> should have fiber D3 hit ! ");

	  auto printW = [](const auto a, const int width) -> std::string {
	    std::stringstream ss;
	    ss << std::fixed << std::right;
	    ss.fill(' ');    // fill space around displayed #
	    ss.width(width); // set  width around displayed #
	    ss << a;
	    return ss.str();
	  };

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

	  continue;
	}

      std::array<double,6> LineParams = {it_ListHitsSim->second[id_FirstFiber][0].hitX,
					 it_ListHitsSim->second[id_FirstFiber][0].hitY,
					 it_ListHitsSim->second[id_FirstFiber][0].hitZ,
					 it_ListHitsSim->second[id_FirstFiber][0].momX,
					 it_ListHitsSim->second[id_FirstFiber][0].momY,
					 it_ListHitsSim->second[id_FirstFiber][0].momZ};


      auto f_line = [](const std::array<double, 6>& params, double posZ, std::array<double, 3>& hitExtrap)
      {
	if(std::abs(params[3]-posZ)<1e-3)
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

      std::vector<int> id_dets = {G4Sol::FiberD3_x, G4Sol::FiberD3_u, G4Sol::FiberD3_v, G4Sol::MiniFiberD1_x1,  G4Sol::MiniFiberD1_u1,  G4Sol::MiniFiberD1_v1,  G4Sol::MiniFiberD1_x2,  G4Sol::MiniFiberD1_u2,  G4Sol::MiniFiberD1_v2, G4Sol::PSFE};

      std::array<double, 3> hitExp = {0.,0.,0.};
      for(int id_det : id_dets)
	{
	  if(id_det < id_FirstFiber)
	    continue;

	  if(it_ListHitsSim->second[id_det].size()>0)
	    {
	      f_line(LineParams, it_ListHitsSim->second[id_det][0].hitZ, hitExp);
	      LocalHisto.h_ResidualFiberX->Fill(G4Sol::nameLiteralDet.begin()[id_det], hitExp[0] - it_ListHitsSim->second[id_det][0].hitX,1.);
	      LocalHisto.h_ResidualFiberY->Fill(G4Sol::nameLiteralDet.begin()[id_det], hitExp[1] - it_ListHitsSim->second[id_det][0].hitY,1.);
	    }

	}
    }
  return 0;
}

template class TCheckFiberTrack<MCAnaEventG4Sol>;
template class TCheckFiberTrack<Ana_WasaEvent>;
