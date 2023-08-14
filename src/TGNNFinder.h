#ifndef TGNNFINDER
#define TGNNFINDER

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"

#include "HitAna/ConstantParameter.hh"
#include "HitAna/FiberHitAna.hh"
#include "HitAna/PSBHitAna.hh"
#include "HitAna/PSFEHitAna.hh"
#include "HitAna/FiberTrackAna.hh"
#include "HitAna/FiberAnalyzer.hh"
#include "HitAna/TrackHit.hh"

#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "ReturnRes.hh"
#include "TDataProcess.h"
#include "THyphiAttributes.h"
#include "HitGnnStruct.hh"

#include "TVector3.h"
#include "TRandom3.h"

#include "torch/torch.h"
//#include "torch/script.h"


template<class Out>
using TDataProcessInterface = TDataProcess<FullRecoEvent, Out>;


template<class Out>
class TGNNFinder final : public TDataProcessInterface<Out> //TDataProcess<FullRecoEvent, Out> //TDataProcessInterface
{
public:
  const THyphiAttributes& att;

  TGNNFinder(const THyphiAttributes& attr);
  ~TGNNFinder();

  // int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator()(FullRecoEvent& RecoEvent, Out* OutTree) override;

private:
  int Exec(FullRecoEvent& RecoEvent, Out* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;

  int FinderGNN(FullRecoEvent& RecoEvent);

  static bool did_order(HitGnn a, HitGnn b){
    if(a.did != b.did) return a.did < b.did;
    else if(a.did>6){
      if (a.lid != b.lid ) return a.lid < b.lid;
      else                 return a.de  > b.de;
    }
    else{
      if (a.s != b.s ) return a.s  < b.s;
      else             return a.de > b.de;
    }
  }


  HitGnn make_cluster(HitGnn hit_a, HitGnn hit_b, double lid, double s, double dl, double phi, double wang1, double wang2);

  bool not_dup(torch::Tensor input_node, std::set<int> label_g1, std::set<int> label_g2);
  std::set<int> get_edge_id_g(torch::Tensor edge_index, std::set<int> label_g);
  std::tuple< std::map<int, std::set<int> >, torch::Tensor, torch::Tensor, bool > get_label_g(
      torch::Tensor input_node, torch::Tensor edge_index, torch::Tensor gnn_label_n, torch::Tensor gnn_label_e,
      std::map<int, std::set<int> > gnn_label_g, torch::Tensor gnn_label_g_bce, torch::Tensor gnn_pred_e,
      int id=-1, int id_s=-1, int id_d=-1, int g_s=0);

  //torch::jit::script::Module model;


};

#endif
