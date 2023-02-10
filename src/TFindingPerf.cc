#include "TFindingPerf.h"

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

#define DEBUG_RIEMANNFINDER

const char TPerf::GenNameHit::Chars[] = { 'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};
const char TPerf::GenNameHit::Nums[] = { '1','2','3','4','5','6','7','8','9','0'};

std::size_t TPerf::Levenshtein_distance2(std::string_view string_a, std::string_view string_b)
{
    const auto size_a = string_a.size()/2;
    const auto size_b = string_b.size()/2;

    std::vector<std::size_t> distances(size_b + 1);
    std::iota(distances.begin(), distances.end(), std::size_t{0});

    for (std::size_t i = 0; i < size_a; ++i)
      {
        std::size_t previous_distance = 0;
        for (std::size_t j = 0; j < size_b; ++j)
	  {
        //cout<<string_a.substr(2*i,2)<<" "<<string_b.substr(2*j,2)<<"\n";
	    distances[j + 1] = std::min({
                std::exchange(previous_distance, distances[j + 1]) + (string_a.substr(2*i,2) == string_b.substr(2*j,2) ? 0 : 1),
                distances[j] + 1,
                distances[j + 1] + 1
	      });
        }
    }
    return distances[size_b];
}

std::size_t TPerf::Levenshtein_distance3(const std::vector<std::string>& string_a, const std::vector<std::string>& string_b)
{
    const auto size_a = string_a.size();
    const auto size_b = string_b.size();

    std::vector<std::size_t> distances(size_b + 1);
    std::iota(distances.begin(), distances.end(), std::size_t{0});

    for (std::size_t i = 0; i < size_a; ++i)
      {
        std::size_t previous_distance = 0;
        for (std::size_t j = 0; j < size_b; ++j)
	  {
	    distances[j + 1] = std::min({
                std::exchange(previous_distance, distances[j + 1]) + (string_a[i] == string_b[j] ? 0 : 1),
                distances[j] + 1,
                distances[j + 1] + 1
	      });
        }
      }
    return distances[size_b];
}

std::size_t TPerf::Levenshtein_distance4(const std::set<std::string>& string_a, const std::set<std::string>& string_b)
{
    const auto size_a = string_a.size();
    const auto size_b = string_b.size();

    std::vector<std::size_t> distances(size_b + 1);
    std::iota(distances.begin(), distances.end(), std::size_t{0});

    for (auto it_a = string_a.begin(), it_a_end = string_a.end();it_a != it_a_end;++it_a)
      {
        std::size_t previous_distance = 0;
        for (auto it_b = string_b.begin(), it_b_end = string_b.end();it_b != it_b_end;++it_b)
	  {
	    int j = std::distance(string_b.begin(), it_b);
	    distances[j + 1] = std::min({
                std::exchange(previous_distance, distances[j + 1]) + (*it_a == *it_b ? 0 : 1),
                distances[j] + 1,
                distances[j + 1] + 1
	      });
        }
      }
    return distances[size_b];
}


std::size_t TPerf::Levenshtein_distance5(const std::set<std::string>& string_a, const std::set<std::string>& string_b)
{
    const auto size_a = string_a.size();
    const auto size_b = string_b.size();

    std::vector<std::size_t> distances_bb(size_b + 1);
    std::vector<std::size_t> distances_b(size_b + 1);
    std::vector<std::size_t> distances_c(size_b + 1);
    std::iota(distances_b.begin(), distances_b.end(), std::size_t{0});

    for (auto it_a = string_a.begin(), it_a_end = string_a.end();it_a != it_a_end;++it_a)
      {
        //std::size_t previous_distance = 0;
	int i = std::distance(string_a.begin(), it_a);
	distances_c[0] = i+1;
	for (auto it_b = string_b.begin(), it_b_end = string_b.end();it_b != it_b_end;++it_b)
	  {
	    int j = std::distance(string_b.begin(), it_b);
	    distances_c[j+1]=distances_b[j]+ 1*(*it_a != *it_b);

	    if (i > 0 && j > 0 && *(std::prev(it_a)) == *it_b && *it_a == *(std::prev(it_b)) && distances_c[j + 1] > distances_bb[j - 1] + 1)
	      distances_c[j + 1] = distances_bb[j - 1] + 1;
	    /* deletion */
	    if (distances_c[j + 1] > distances_b[j + 1] + 1)
	      distances_c[j + 1] = distances_b[j + 1] + 1;
	    /* insertion */
	    if (distances_c[j + 1] > distances_c[j] + 1)
	      distances_c[j + 1] = distances_c[j] + 1;

	  }

	std::swap(distances_bb,distances_b);
	std::swap(distances_b,distances_c);

      }
    return distances_b[size_b];
}



using namespace std;
using namespace G4Sol;

template<class Out>
TFindingPerf<Out>::TFindingPerf(const THyphiAttributes& attribut) : TDataProcessInterface<Out>("FindingPerf"), att(attribut)
{

}

template<class Out>
TFindingPerf<Out>::~TFindingPerf()
{

}

template<class Out>
void TFindingPerf<Out>::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }

template<class Out>
ReturnRes::InfoM TFindingPerf<Out>::operator()(FullRecoEvent& RecoEvent, Out* OutTree)
{

  int result_finder = Exec(RecoEvent, OutTree);

  return SoftExit(result_finder);
}

template<class Out>
int TFindingPerf<Out>::Exec(FullRecoEvent& RecoEvent, Out* )
{
  return CheckTrackFinding(RecoEvent);
}

template<class Out>
ReturnRes::InfoM TFindingPerf<Out>::SoftExit(int ) { return ReturnRes::Fine; }

template<class Out>
void TFindingPerf<Out>::SelectHists()
{
  LocalHisto.h_PerfFinder =  this->AnaHisto->CloneAndRegister(this->AnaHisto->h_PerfFinder);
  LocalHisto.h_PerfFinderLevenshtein = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_PerfFinderLevenshtein);
  // auto Yaxis = LocalHisto.h_PerfFinder->GetYaxis();
  // Yaxis->SetBinLabel(1,"NotFoundReal");
  // Yaxis->SetBinLabel(2,"FoundReal");
  // Yaxis->SetBinLabel(3,"FoundRealClones");
  // Yaxis->SetBinLabel(4,"Found_>90%");
  // Yaxis->SetBinLabel(5,"Found_>70%");
  // Yaxis->SetBinLabel(6,"Found_>50%");
  // Yaxis->SetBinLabel(7,"Found_70%90%");
  // Yaxis->SetBinLabel(8,"Found_>70%");
  // Yaxis->SetBinLabel(9,"Found_>50%");
  // Yaxis->SetBinLabel(10,"Found_50%70%");
  // Yaxis->SetBinLabel(11,"Found_>50%");
  // Yaxis->SetBinLabel(12,"SharedReal");
  // Yaxis->SetBinLabel(13,"Shared_>90%");
  // Yaxis->SetBinLabel(14,"Shared_70%90%");
  // Yaxis->SetBinLabel(15,"SharedGhost");
  // Yaxis->SetBinLabel(16,"Shared_50%70%");
  // Yaxis->SetBinLabel(17,"SharedGhost");
  // Yaxis->SetBinLabel(18,"Shared_30%50%");
  // Yaxis->SetBinLabel(19,"SharedGhost");
  // Yaxis->SetBinLabel(20,"Shared_10%30%");
  // Yaxis->SetBinLabel(21,"SharedGhost");
  // Yaxis->SetBinLabel(22,"Shared_<10%");

}

template<class Out>
int TFindingPerf<Out>::CheckTrackFinding(const FullRecoEvent& RecoEvent)
{

  auto minDet = G4Sol::MG01;
  auto maxDet = G4Sol::PSBE;

  std::unordered_map<int, std::map<int,std::vector<int>>> realTracks;

  using key_t1 = std::tuple<int, int>;
  struct key_hash : public std::unary_function<key_t1, std::size_t>
  {
    std::size_t operator()(const key_t1& k) const
    {
      return std::get<0>(k) ^ std::get<1>(k);
    }
  };


  using DetHitmap_t = std::unordered_map<key_t1, std::string, key_hash,std::equal_to<std::tuple<int,int>>>;
  DetHitmap_t nameHits;

  TPerf::GenNameHit genNames;

  att._logger->debug("ListNameHits : {} ",nameHits.size());
  for(int idDet = minDet; idDet <= maxDet;++idDet)
    {
      for(int idHit = 0; idHit < RecoEvent.ListHitsToTracks[idDet].size();++idHit)
	{
	  nameHits.insert({std::make_tuple(idDet,idHit),genNames.nextName()});
	  att._logger->debug(" -> {} {}",idDet, idHit);
	  auto mes = RecoEvent.ListHits[idDet][idHit].get();

	  att._logger->debug("nameHit : {} {} # {} / idT {} / [{} {}] = {}",idDet,G4Sol::nameLiteralDet.begin()[idDet], idHit, RecoEvent.ListHitsToTracks[idDet][idHit], mes->getDetId(), mes->getHitId(), nameHits[std::make_tuple(idDet,idHit)]);

	  auto it_idT = realTracks.find(RecoEvent.ListHitsToTracks[idDet][idHit]);
	  if(it_idT == realTracks.end())
	    {
	      std::map<int, std::vector<int> > tempT;
	      std::vector<int> tempH = {idHit};
	      tempT.insert({idDet,tempH});
	      realTracks.insert({RecoEvent.ListHitsToTracks[idDet][idHit],tempT});
	    }
	  else
	    {
	      auto it_detT = it_idT->second.find(idDet);
	      if(it_detT == it_idT->second.end())
		{
		  std::vector<int> tempH = {idHit};
		  it_idT->second.insert({idDet,tempH});
		}
	      else
		it_detT->second.emplace_back(idHit);
	    }
	}
    }


  std::unordered_map<int,std::set<std::string>> nameRealTracks;
  for(const auto& [idTrack, Rtrack] : realTracks)
    {
      std::set<std::string> nameTrack;
      for(const auto& [idDet, idHits] : Rtrack)
	for(const auto& idHit : idHits)
	  nameTrack.insert(nameHits[std::make_tuple(idDet,idHit)]);

      nameRealTracks.insert({idTrack,nameTrack});
    }

  std::vector<std::unordered_map<int,std::vector<std::tuple<int,int>>>> FoundToReal;
  std::vector<std::set<std::string>> nameTrackFound;

  for(const auto& trackC : RecoEvent.TracksFound)
    {
      std::unordered_map<int,std::vector<std::tuple<int,int>>> trackIDs;
      std::set<std::string> nameTrack;
      for(int id_hit : trackC.orderedHitIds)
	{
	  const auto& Ids = RecoEvent.IdHitsToMeasurement[id_hit];
	  int realTrackId = RecoEvent.ListHitsToTracks[Ids.id_det][Ids.id_hit];

	  nameTrack.insert(nameHits[std::make_tuple(Ids.id_det,Ids.id_hit)]);

	  auto it_idT = trackIDs.find(realTrackId);
	  if(it_idT == trackIDs.end())
	    {
	      std::vector<std::tuple<int,int>> tempVec;
	      tempVec.emplace_back(Ids.id_det,Ids.id_hit);
	      trackIDs.insert({realTrackId,tempVec});
	    }
	  else
	    it_idT->second.emplace_back(Ids.id_det,Ids.id_hit);
	}
      nameTrackFound.push_back(nameTrack);
      FoundToReal.push_back(trackIDs);
    }


  att._logger->debug("Real Tracks: {}",nameRealTracks.size());
  for(const auto& [idT, nT] : nameRealTracks)
    att._logger->debug("{}: {}",idT,std::accumulate(nT.begin(),nT.end(),std::string(),std::plus<std::string>()));

  att._logger->debug("Found Tracks: {}",nameTrackFound.size());
  for(const auto& nameT : nameTrackFound)
    att._logger->debug("{}",std::accumulate(nameT.begin(),nameT.end(),std::string(),std::plus<std::string>()));

  for(const auto& [idT, nT] : nameRealTracks)
    {
      int best_dist = 999;
      for(const auto& nameT : nameTrackFound)
	{
	  int dist = TPerf::Levenshtein_distance5(nameT,nT);
	  att._logger->debug("Dist Levenshtein : {} = {} / {}: {}",dist,std::accumulate(nameT.begin(),nameT.end(),std::string(),std::plus<std::string>()),idT,std::accumulate(nT.begin(),nT.end(),std::string(),std::plus<std::string>()));

	  int max_dist = TMath::Max(nameT.size(),nT.size());
	  if(dist >= max_dist)
	    {
	      att._logger->debug(" -> completely diff ! max {} dist {}",max_dist,dist);
	      continue;
	    }
	  double perC = 1. - dist / static_cast<double>(nT.size());
	  if(perC >= 0.9)
	    att._logger->debug(" -> Found ! dist {} real {} % = {}",dist,nT.size(),perC);
	  else
	    att._logger->debug(" -> Shared : dist {} real {} % = {}",dist,nT.size(),perC);

	  if(dist<=best_dist)
	    best_dist = dist;

	}
      LocalHisto.h_PerfFinderLevenshtein->Fill(best_dist,nT.size());
    }


  struct Results {
    int nFound = 0;
    int nShared = 0;
    int nFound90per = 0; // > 90%
    int nFound70per = 0; // [70%, 90%]
    int nFound50per = 0; // [50%, 70%]
    int nShared90per = 0; // > 90%
    int nShared70per = 0; // [70%, 90%]
    int nShared50per = 0; // [50%, 70%]
    int nShared30per = 0; // [30%, 50%]
    int nShared10per = 0; // [10%, 30%]
    int nShared0per = 0;
  };

  auto f_setPer = [](double Per, bool FS, Results& res) {
    if(FS == true)
      {
	if(Per>=0.9)
	  res.nFound90per += 1;
	else if(Per >= 0.7)
	  res.nFound70per += 1;
	else if(Per >= 0.5)
	  res.nFound50per += 1;
      }
    else
      {
	if(Per>=0.9)
	  res.nShared90per += 1;
	else if(Per >= 0.7)
	  res.nShared70per += 1;
	else if(Per >= 0.5)
	  res.nShared50per += 1;
	else if(Per >= 0.3)
	  res.nShared30per += 1;
	else if(Per >= 0.1)
	  res.nShared10per += 1;
	else
	  res.nShared0per += 1;
      }
  };

  std::unordered_map<int, Results > TrackSummary;

  for(const auto& trackC_toReal : FoundToReal)
    {
      att._logger->debug("-- FoundTrack :");
      for(const auto& [idTrack, vecDetHits] : trackC_toReal)
	{
	  att._logger->debug("part of Track {},",idTrack);
	  for(const auto& [idDet,idHit] : vecDetHits)
	    att._logger->debug("dethit {} {}",idDet,idHit);
	}


      if(trackC_toReal.size()==1)
	{
	  att._logger->debug("trackC_toReal == 1");

	  auto itT = trackC_toReal.begin();
	  int idFound = itT->first;
	  if(auto itS = TrackSummary.find(idFound); itS == TrackSummary.end())
	    {
	      TrackSummary.insert({idFound, {1,0,0,0,0,0,0,0,0,0} });
	      att._logger->debug("new Summary, {}",idFound);
	    }
	  else
	    {
	      itS->second.nFound += 1;
	      att._logger->debug("add+1 Summary, {}",idFound);
	    }


	  const auto& it_realT = realTracks.find(idFound);
	  int nHit = 0;
	  int nHitFound = 0;
	  for(const auto& [Found_Det,Found_Hit] : itT->second)
	    {
	      att._logger->debug("find hits :{} {}",Found_Det,Found_Hit);
	      auto itDet = it_realT->second.find(Found_Det);
	      for(auto realHit : itDet->second)
		{
		  att._logger->debug(" -- real hit : {}",realHit);
		  if(realHit == Found_Hit)
		    ++nHitFound;
		}
	    }

	  for(const auto& itDet : it_realT->second)
	    nHit += itDet.second.size();

	  double perFound = nHitFound / static_cast<double>(nHit);
	  att._logger->debug("Stats found : {} {} / {} ",nHitFound,nHit,perFound);
	  const auto& TempRes = TrackSummary.find(idFound);
	  f_setPer(perFound,true,TempRes->second);
	}
      else
	{
	  att._logger->debug("trackC_toReal > 1 -> sharing hits");

	  std::unordered_map<int, std::tuple<int,int> > trackShares; // trackID -> {shared, total}
	  for(const auto& [idTrack, DetHitId] : trackC_toReal)
	    {
	      if(const auto& itS = TrackSummary.find(idTrack); itS == TrackSummary.end())
		{
		  TrackSummary.insert({idTrack,{0,1,0,0,0,0,0,0,0,0}});
		  att._logger->debug("new Summary shared, {}",idTrack);
		}
	      else
		{
		  itS->second.nShared += 1;
		  att._logger->debug("add+1 Summary shared, {}",idTrack);
		}

	      const auto& it_realT = realTracks.find(idTrack);
	      int nHit = 0;
	      int nHitShared = 0;
	      for(const auto& [Found_Det,Found_Hit] : DetHitId)
		{
		  auto itDet = it_realT->second.find(Found_Det);
		  for(auto realHit : itDet->second)
		    if(realHit == Found_Hit)
		      ++nHitShared;
		}
	      for(const auto& itDet : it_realT->second)
		nHit += itDet.second.size();

	      trackShares.insert({idTrack, {nHitShared, nHit}});
	      double perFound = nHitShared / static_cast<double>(nHit);
	      att._logger->debug("Stats found : {} {} / {} ",nHitShared,nHit,perFound);

	    }

	  for(const auto& [idTrack, Nhits] : trackShares)
	    {
	      const auto& itS = TrackSummary.find(idTrack);
	      if(itS == TrackSummary.end())
		{
		  att._logger->warn("Perf : Shared track summary should have been found ! {}", idTrack);
		  continue;
		}
	      double perC = std::get<0>(Nhits) / static_cast<double>(std::get<1>(Nhits));
	      f_setPer(perC, false, itS->second);
	    }
	}
    }
  int n_realTracks = realTracks.size();
  for(int n_i = 1;n_i<n_realTracks;++n_i)
    LocalHisto.h_PerfFinder->Fill("TotalReal",n_i,n_realTracks);

  LocalHisto.h_PerfFinder->Fill("TotalReal_dN",n_realTracks,n_realTracks);

  for(const auto& [idTrack,Res] : TrackSummary)
    {
      att._logger->debug("Stats : trackID {} / {}",idTrack, n_realTracks);
      att._logger->debug(" - Found  > {} | 90%={} 70%={} 50%={}",Res.nFound,Res.nFound90per,Res.nFound70per,Res.nFound50per);
      att._logger->debug(" - Shared > {} | 90%={} 70%={} 50%={} 30%={} 10%={} 0%={}",Res.nShared,Res.nShared90per,Res.nShared70per,Res.nShared50per, Res.nShared30per, Res.nShared10per, Res.nShared0per);

      if(Res.nFound == 0)
	LocalHisto.h_PerfFinder->Fill("NotFoundReal",n_realTracks,1.);
      else
	{
	  if(Res.nFound == 1)
	    LocalHisto.h_PerfFinder->Fill("FoundReal",n_realTracks,1.);
	  else
	    LocalHisto.h_PerfFinder->Fill("FoundRealClones",n_realTracks,1.);

	  if(Res.nFound90per>0)
	    {
	      LocalHisto.h_PerfFinder->Fill("Found_>90%",n_realTracks,Res.nFound90per);
	      LocalHisto.h_PerfFinder->Fill("Found_>70%",n_realTracks,Res.nFound90per);
	      LocalHisto.h_PerfFinder->Fill("Found_>50%",n_realTracks,Res.nFound90per);
	    }
	  if(Res.nFound70per>0)
	    {
	      LocalHisto.h_PerfFinder->Fill("Found_70%90%",n_realTracks,Res.nFound70per);
	      LocalHisto.h_PerfFinder->Fill("Found_>70%",n_realTracks,Res.nFound70per);
	      LocalHisto.h_PerfFinder->Fill("Found_>50%",n_realTracks,Res.nFound70per);
	    }
	  if(Res.nFound50per>0)
	    {
	      LocalHisto.h_PerfFinder->Fill("Found_50%70%",n_realTracks,Res.nFound50per);
	      LocalHisto.h_PerfFinder->Fill("Found_>50%",n_realTracks,Res.nFound50per);
	    }
	}

      if(Res.nShared!=0)
	{
	  LocalHisto.h_PerfFinder->Fill("SharedReal",n_realTracks,Res.nShared);
	  if(Res.nShared90per>0)
	    {
	      LocalHisto.h_PerfFinder->Fill("Shared_>90%",n_realTracks,Res.nShared90per);
	    }
	  if(Res.nShared70per>0)
	    {
	      LocalHisto.h_PerfFinder->Fill("Shared_70%90%",n_realTracks,Res.nShared70per);
	      if(Res.nFound)
		LocalHisto.h_PerfFinder->Fill("SharedGhost",n_realTracks,Res.nShared50per);
	    }
	  if(Res.nShared50per>0)
	    {
	      LocalHisto.h_PerfFinder->Fill("Shared_50%70%",n_realTracks,Res.nShared50per);
	      if(Res.nFound)
		LocalHisto.h_PerfFinder->Fill("SharedGhost",n_realTracks,Res.nShared50per);
	    }
	  if(Res.nShared30per>0)
	    {
	      LocalHisto.h_PerfFinder->Fill("Shared_30%50%",n_realTracks,Res.nShared30per);
	      if(Res.nFound)
		LocalHisto.h_PerfFinder->Fill("SharedGhost",n_realTracks,Res.nShared50per);
	    }
	  if(Res.nShared10per>0)
	    {
	      LocalHisto.h_PerfFinder->Fill("Shared_10%30%",n_realTracks,Res.nShared10per);
	      if(Res.nFound)
		LocalHisto.h_PerfFinder->Fill("SharedGhost",n_realTracks,Res.nShared50per);
	    }
	  if(Res.nShared0per>0)
	    {
	      LocalHisto.h_PerfFinder->Fill("Shared_<10%",n_realTracks,Res.nShared0per);
	    }
	}
    }

  return 0;
}

template class TFindingPerf<MCAnaEventG4Sol>;
template class TFindingPerf<Ana_WasaEvent>;
