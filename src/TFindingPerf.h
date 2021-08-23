#ifndef TFINDINGPERF
#define TFINDINGPERF

#include "ReturnRes.hh"
#include "TDataProcess.h"
#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "THyphiAttributes.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"

#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TH2I.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>


namespace TPerf {

  struct GenNameHit{

    static const char Chars[];
    static const char Nums[];

    int nextChar = 0;
    int nextNum  = 0;

    std::string&& nextName()
    {
      int currentC = nextChar;
      int currentN = nextNum;
      if(currentC == 25)
	{
	  nextChar = 0;
	  ++nextNum;
	}
      else
	++nextChar;

      std::string tempN (1,Chars[currentC]);
      tempN.append(1,Nums[currentN]);
      return std::move(tempN);
    }
  };

  //  std::size_t Levenshtein_distance(std::string_view string_a, std::string_view string_b);

// std::size_t Levenshtein_distance(std::string_view string_a, std::string_view string_b)
// {
//     const auto size_a = string_a.size();
//     const auto size_b = string_b.size();

//     std::vector<std::size_t> distances(size_b + 1);
//     std::iota(distances.begin(), distances.end(), std::size_t{0});

//     for (std::size_t i = 0; i < size_a; ++i)
//       {
//         std::size_t previous_distance = 0;
//         for (std::size_t j = 0; j < size_b; ++j)
// 	  {
// 	    distances[j + 1] = std::min({
//                 std::exchange(previous_distance, distances[j + 1]) + (string_a[i] == string_b[j] ? 0 : 1),
//                 distances[j] + 1,
//                 distances[j + 1] + 1
// 	      });
//         }
//     }
//     return distances[size_b];
// };

  std::size_t Levenshtein_distance2(std::string_view string_a, std::string_view string_b);
  std::size_t Levenshtein_distance3(const std::vector<std::string>& string_a, const std::vector<std::string>& string_b);
  std::size_t Levenshtein_distance4(const std::set<std::string>& string_a, const std::set<std::string>& string_b);
  std::size_t Levenshtein_distance5(const std::set<std::string>& string_a, const std::set<std::string>& string_b);
// {
//     const auto size_a = string_a.size();
//     const auto size_b = string_b.size();

//     std::vector<std::size_t> distances(size_b + 1);
//     std::iota(distances.begin(), distances.end(), std::size_t{0});

//     for (std::size_t i = 0; i < size_a; ++i)
//       {
//         std::size_t previous_distance = 0;
//         for (std::size_t j = 0; j < size_b; ++j)
// 	  {
// 	    distances[j + 1] = std::min({
//                 std::exchange(previous_distance, distances[j + 1]) + (string_a.substr(2*i,2) == string_b.substr(2*i,2) ? 0 : 1),
//                 distances[j] + 1,
//                 distances[j + 1] + 1
// 	      });
//         }
//     }
//     return distances[size_b];
// };

};

typedef TDataProcess<FullRecoEvent,MCAnaEventG4Sol> TDataProcessInterface;


class TFindingPerf final :  public TDataProcessInterface
{
  public :
  const THyphiAttributes& att;

  TFindingPerf(const THyphiAttributes& attr);
  ~TFindingPerf();

  //int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator() (FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree) override;
 private:
  int Exec(FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;
  int CheckTrackFinding(const FullRecoEvent& RecoEvent);

  struct LocalHists
  {
    TH2F* h_PerfFinder;
    TH2F* h_PerfFinderLevenshtein;
  };
  LocalHists LocalHisto;
};


#endif
