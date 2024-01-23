//#include <torch/torch.h>
#include <torch/script.h>
#include "TGNNFinder.h"
#include "HitGnnStruct.hh"

using namespace std;
using namespace G4Sol;

template<class Out>
TGNNFinder<Out>::TGNNFinder(const THyphiAttributes& attribut) : TDataProcessInterface<Out>("GNN_Find"), att(attribut)
{
  att._logger->info("TGNNFinder::TGNNFinder");
  //torch::jit::script::Module model;
  try{
    model = torch::jit::load("scripted_model/scripted_model_m7_l4_5_h180_480_b25_e50_d040_lr000004_v04_div0005.pt");
                             //"scripted_model/scripted_model_m7_l4_5_h180_480_b25_e50_d040_lr000004_v04_div0005.pt"
    std::cout << "model is loaded" << std::endl;
  }
  catch(const c10::Error& e){
    std::cerr << "error loading the model\n";
    return;
  }
  model.eval();

  if(att.GNN_Text)
    ifs_gnn.open(att.GNN_Node);
}

template<class Out>
TGNNFinder<Out>::~TGNNFinder() {}

template<class Out>
void TGNNFinder<Out>::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }

template<class Out>
ReturnRes::InfoM TGNNFinder<Out>::operator()(FullRecoEvent& RecoEvent, Out* OutTree)
{

  int result_finder = Exec(RecoEvent, OutTree);

  return SoftExit(result_finder);
}

template<class Out>
int TGNNFinder<Out>::Exec(FullRecoEvent& RecoEvent, Out* ) {
 
  if(att.GNN_Text)
    return FinderGNNText(RecoEvent);
  else
    return FinderGNN(RecoEvent);
}

template<class Out>
ReturnRes::InfoM TGNNFinder<Out>::SoftExit(int result_full) { return ReturnRes::Fine; }


template<class Out>
void TGNNFinder<Out>::SelectHists() {}


template<class Out>
int TGNNFinder<Out>::FinderGNN(FullRecoEvent& RecoEvent)
{

  double mm2cm = 0.1;

  //torch::jit::script::Module model;
  //try{
  //  model = torch::jit::load("scripted_model/scripted_model_m7_l4_5_h180_480_b25_e50_d040_lr000004_v04_div0005.pt");
  //                           //"scripted_model/scripted_model_m7_l4_5_h180_480_b25_e50_d040_lr000004_v04_div0005.pt"
  //  std::cout << "model is loaded" << std::endl;
  //}
  //catch(const c10::Error& e){
  //  std::cerr << "error loading the model\n";
  //  return -1;
  //}
  //model.eval();
  auto bce_loss = torch::nn::BCELoss();

  std::vector<HitGnn> hitgnn_cont;

  std::map<int, std::set<int> > gnn_label_g;
  for(auto fiber_hit_det: RecoEvent.FiberHitClCont){
    for(auto fiber_hit_lay: fiber_hit_det){
      for(auto fiber_hit: fiber_hit_lay){
        if(fiber_hit->GetDet()!=3 && fiber_hit->GetDet()!=4) continue;
        HitGnn buf;
        buf.did  = (fiber_hit->GetDet() - 3)*3 + fiber_hit->GetLay() + 1;
        buf.lid  = fiber_hit->GetFib();
        buf.z = fiber_hit->GetZ() * mm2cm;
        buf.s    = fiber_hit->GetPos() * mm2cm;
        buf.phi  = 0;
        buf.r    = 0;
        buf.dl   = 0;
        buf.dle  = 0.2 * mm2cm; //tmp
        buf.wang1 = 90. * TMath::DegToRad();
        buf.wang2 = 90. * TMath::DegToRad() + fiber_hit->GetAng();
        buf.time = 0;
        buf.de   = fiber_hit->GetTOT()/70.;
        buf.clsize = fiber_hit->GetClsize();
        buf.cls = fiber_hit->GetFib();
        buf.cle = fiber_hit->GetFib();
        hitgnn_cont.push_back(buf);
      }
    }
  }
  att._logger->debug("TGNNFinder: Fiber  hit end {}",hitgnn_cont.size());

  for(auto mdc_hit_lay: RecoEvent.MDCHitCont){
    for(auto mdc_hit: mdc_hit_lay){
      HitGnn buf;
      buf.did  = mdc_hit->GetLay() + 7;
      buf.lid  = mdc_hit->GetWir();;
      buf.time = 0;
      buf.clsize = 1;
      buf.cls = mdc_hit->GetWir();;
      buf.cle = mdc_hit->GetWir();;
      std::array<double, 6> edge;
      edge = mdc_hit->GetEdge();
      TVector3 pw(edge[3]-edge[0], edge[4]-edge[1], edge[5]-edge[2]);
      double wang1 = pw.Theta();
      double wang2 = pw.Phi();
      if(wang1<1.0e-5){
        wang1 = 0.;
        wang2 = 0.;
      }
      buf.wang1 = wang1;
      buf.wang2 = wang2;
      buf.z = (edge[2] + edge[5]) / 2. * mm2cm;
      buf.s     = 0.;
      buf.r     = sqrt( pow((edge[0]+edge[3])/2.,2.) + pow((edge[1]+edge[4])/2.,2.)  ) * mm2cm;
      buf.phi   = atan2( (edge[1]+edge[4])/2, (edge[0]+edge[3])/2. );
      buf.dl    = mdc_hit->GetDl() * mm2cm;
      buf.dle   = 0.3 * mm2cm; //tmp
      buf.de    = (double) (mdc_hit->GetTot()) / 250.;
      hitgnn_cont.push_back(buf);
    }
  }
  att._logger->debug("TGNNFinder: MDC    hit end {}",hitgnn_cont.size());


  for(auto psb_hit: RecoEvent.PSBHitCont){
    HitGnn buf;
    buf.did  = 24;
    buf.lid  = psb_hit->GetSeg();;
    buf.time = 0;
    //buf.de   = 0;
    buf.clsize = 1;
    buf.cls = psb_hit->GetSeg();;
    buf.cle = psb_hit->GetSeg();;
    std::array<double, 6> edge;
    edge = psb_hit->GetEdge();
    TVector3 pw(edge[3]-edge[0], edge[4]-edge[1], edge[5]-edge[2]);
    double wang1 = pw.Theta();
    double wang2 = pw.Phi();
    if(wang1<1.0e-5){
      wang1 = 0.;
      wang2 = 0.;
    }
    buf.wang1 = wang1;
    buf.wang2 = wang2;
    buf.z = (edge[2] + edge[5]) / 2. * mm2cm;
    buf.s     = 0.;
    buf.r     = sqrt( pow((edge[0]+edge[3])/2.,2.) + pow((edge[1]+edge[4])/2.,2.)  ) * mm2cm;
    buf.phi   = atan2( (edge[1]+edge[4])/2, (edge[0]+edge[3])/2. );
    buf.dl    = 0.;
    buf.dle   = 1.1;
    buf.de    = 1.0;
    hitgnn_cont.push_back(buf);
  }
  att._logger->debug("TGNNFinder: PSB    hit end {}",hitgnn_cont.size());


  for(auto psfe_hit: RecoEvent.PSFEHitCont){
    HitGnn buf;
    buf.did  = 25;
    buf.lid  = psfe_hit->GetSeg();;
    buf.z = psfe_hit->GetZ() * mm2cm;
    buf.phi  = psfe_hit->GetPhi();
    buf.time = 0;
    //buf.de   = 0;
    buf.clsize = 1;
    buf.cls = psfe_hit->GetSeg();;
    buf.cle = psfe_hit->GetSeg();;
    buf.wang1 = 90. * TMath::DegToRad();
    buf.wang2 = psfe_hit->GetPhi();
    buf.r = 0.;
    buf.s = 0.;
    buf.dl = 0.;
    buf.dle = 0.7;
    buf.de = 1.0;
    hitgnn_cont.push_back(buf);
  }
  att._logger->debug("TGNNFinder: PSFE   hit end {}",hitgnn_cont.size());


  // delete overlap
  std::sort(hitgnn_cont.begin(), hitgnn_cont.end(), did_order);
  std::vector<HitGnn> hitgnn_cont_ov;
  std::map< double, std::set<double> > hit_map;
  for(int i=0; i<(int)hitgnn_cont.size(); ++i){
    if( hit_map.find(hitgnn_cont[i].did)==hit_map.end() ){
      std::set<double> buf_set;
      buf_set.insert(hitgnn_cont[i].lid);
      hit_map[hitgnn_cont[i].did] = buf_set;
      hitgnn_cont_ov.emplace_back(hitgnn_cont[i]);
    }
    else{
      if( hit_map[hitgnn_cont[i].did].find(hitgnn_cont[i].lid)==hit_map[hitgnn_cont[i].did].end() ){
        hit_map[hitgnn_cont[i].did].insert(hitgnn_cont[i].lid);
        hitgnn_cont_ov.emplace_back(hitgnn_cont[i]);
      }
      else continue;
    }
  }
  hitgnn_cont = hitgnn_cont_ov;

  att._logger->debug("TGNNFinder: delete overlap {}",hitgnn_cont.size());

  //  clustering
  std::sort(hitgnn_cont.begin(), hitgnn_cont.end(), did_order);
  std::vector<HitGnn> hitgnn_cont_cl;
  std::vector<HitGnn> hitgnn_cont_cl_layer;
  for(int i=1; i<26; ++i){
    hitgnn_cont_cl_layer.clear();

    for(int j=0; j<(int)hitgnn_cont.size(); ++j){
      if(hitgnn_cont[j].did < i) continue;
      if(hitgnn_cont[j].did > i) break;
      if(hitgnn_cont[j].did==i) hitgnn_cont_cl_layer.emplace_back(hitgnn_cont[j]);
    }

    if(i<=6){
      for(auto hit: hitgnn_cont_cl_layer){
        hitgnn_cont_cl.emplace_back(hit);
      }
    }

    else{
      int size_layer = (int)hitgnn_cont_cl_layer.size();

      for(int j=0; j<size_layer; ++j){
        //std::cout << Form("i : %d,  j : %d,  size_layer : %d",i, j, size_layer) << std::endl;
        HitGnn hit_a = hitgnn_cont_cl_layer[j];
        //std::cout << Form("- hit_a  did : %d,  lid : %.1f", (int)hit_a.did, hit_a.lid) << std::endl;
        if(j == size_layer-1){
          continue;
        }
        else{
          HitGnn hit_b = hitgnn_cont_cl_layer[j+1];
          //std::cout << Form("- hit_b  did : %d,  lid : %.1f", (int)hit_b.did, hit_b.lid) << std::endl;

          if(abs(hit_a.cle - hit_b.cls)==1){
            if(hit_a.clsize==1 && hit_b.clsize==1){
              double lid = (hit_a.lid + hit_b.lid) / 2;
              double dl_a = hit_a.dl;
              double dl_b = hit_b.dl;
              if(dl_a + dl_b == 0){
                dl_a = 1;
                dl_b = 1;
              }
              double phi_a = hit_a.phi;
              double phi_b = hit_b.phi;
              if(      (phi_a - phi_b)>TMath::Pi() ) phi_a -= TMath::Pi()*2;
              else if( (phi_b - phi_a)>TMath::Pi() ) phi_b -= TMath::Pi()*2;
              double phi = (phi_a*dl_b + phi_b*dl_a) / (dl_a + dl_b);
              if(phi<-TMath::Pi()) phi += TMath::Pi()*2;
              if(phi> TMath::Pi()) phi -= TMath::Pi()*2;
              double wang2_a = hit_a.wang2;
              double wang2_b = hit_b.wang2;
              if(      (wang2_a - wang2_b)>TMath::Pi() ) wang2_a -= TMath::Pi()*2;
              else if( (wang2_b - wang2_a)>TMath::Pi() ) wang2_b -= TMath::Pi()*2;
              double wang2 = (wang2_a*dl_b + wang2_b*dl_a) / (dl_a + dl_b);
              if(wang2<-TMath::Pi()) wang2 += TMath::Pi()*2;
              if(wang2> TMath::Pi()) wang2 -= TMath::Pi()*2;
              double wang1_a = hit_a.wang1;
              double wang1_b = hit_b.wang1;
              double wang1 = (wang1_a*dl_b + wang1_b*dl_a) / (dl_a + dl_b);
              double dl = 0;
              double s  = 0;
              HitGnn hit_buf = make_cluster(hit_a, hit_b, lid, s, dl, phi, wang1, wang2);
              hitgnn_cont_cl_layer.erase(hitgnn_cont_cl_layer.begin() + j );
              hitgnn_cont_cl_layer.erase(hitgnn_cont_cl_layer.begin() + j );
              hitgnn_cont_cl_layer.insert(hitgnn_cont_cl_layer.begin()+ j  , hit_buf);
              j--;
              size_layer--;
            }
            else if(hit_a.clsize + hit_b.clsize < 4){
              double lid = (hit_a.lid*hit_a.clsize + hit_b.lid*hit_b.clsize) / (hit_a.clsize + hit_b.clsize);
              double phi_a = hit_a.phi;
              double phi_b = hit_b.phi;
              if(      (phi_a - phi_b)>TMath::Pi() ) phi_a -= TMath::Pi()*2;
              else if( (phi_b - phi_a)>TMath::Pi() ) phi_b -= TMath::Pi()*2;
              double phi = (phi_a*hit_a.clsize + phi_b*hit_b.clsize) / (hit_a.clsize + hit_b.clsize);
              if(phi<-TMath::Pi()) phi += TMath::Pi()*2;
              if(phi> TMath::Pi()) phi -= TMath::Pi()*2;
              double wang2_a = hit_a.wang2;
              double wang2_b = hit_b.wang2;
              if(      (wang2_a - wang2_b)>TMath::Pi() ) wang2_a -= TMath::Pi()*2;
              else if( (wang2_b - wang2_a)>TMath::Pi() ) wang2_b -= TMath::Pi()*2;
              double wang2 = (wang2_a*hit_a.clsize + wang2_b*hit_b.clsize) / (hit_a.clsize + hit_b.clsize);
              if(wang2<-TMath::Pi()) wang2 += TMath::Pi()*2;
              if(wang2> TMath::Pi()) wang2 -= TMath::Pi()*2;
              double wang1_a = hit_a.wang1;
              double wang1_b = hit_b.wang1;
              double wang1 = (wang1_a*hit_a.clsize + wang1_b*hit_b.clsize) / (hit_a.clsize + hit_b.clsize);
              double dl = 0;
              double s  = 0;
              HitGnn hit_buf = make_cluster(hit_a, hit_b, lid, s, dl, phi, wang1, wang2);
              hitgnn_cont_cl_layer.erase(hitgnn_cont_cl_layer.begin() + j );
              hitgnn_cont_cl_layer.erase(hitgnn_cont_cl_layer.begin() + j );
              hitgnn_cont_cl_layer.insert(hitgnn_cont_cl_layer.begin()+ j  , hit_buf);
              j--;
              size_layer--;
            }
          }
          else{
            continue;
          }
        }
      }

      for(auto hit: hitgnn_cont_cl_layer){
        hitgnn_cont_cl.emplace_back(hit);
      }

    }

  }
  hitgnn_cont = hitgnn_cont_cl;

  att._logger->debug("TGNNFinder: clustering     {}",hitgnn_cont.size());


  torch::Tensor input_node = torch::full({(int)hitgnn_cont.size(), 13}, -999, torch::TensorOptions().dtype(torch::kFloat));
  for(int i=0; i<(int)hitgnn_cont.size(); ++i){
    int dtype = -1;
    int layer_buf = hitgnn_cont[i].did;
    if(     layer_buf<=6)                dtype = 0;
    else if(layer_buf>6 && layer_buf<=23)dtype = 1;
    else if(layer_buf==24)               dtype = 2;
    else if(layer_buf==25)               dtype = 3;
    input_node[i][ 0] = hitgnn_cont[i].did;
    input_node[i][ 1] = cos(hitgnn_cont[i].phi);
    input_node[i][ 2] = sin(hitgnn_cont[i].phi);
    input_node[i][ 3] = hitgnn_cont[i].wang1;
    input_node[i][ 4] = cos(hitgnn_cont[i].wang2);
    input_node[i][ 5] = sin(hitgnn_cont[i].wang2);
    input_node[i][ 6] = hitgnn_cont[i].r;
    input_node[i][ 7] = hitgnn_cont[i].z;
    input_node[i][ 8] = hitgnn_cont[i].s;
    input_node[i][ 9] = hitgnn_cont[i].dl;
    input_node[i][10] = hitgnn_cont[i].dle;
    //input_node[i][11] = hitgnn_cont[i].de;
    input_node[i][11] = hitgnn_cont[i].clsize;
    input_node[i][12] = dtype;
  }

  std::vector<int> gnn_src;
  std::vector<int> gnn_dst;

  for(int i=0; i<(int)hitgnn_cont.size()-1; ++i){
    for(int j=i+1; j<(int)hitgnn_cont.size(); ++j){
      int layer_i  = hitgnn_cont[i].did;
      int layer_j  = hitgnn_cont[j].did;
      double phi_i = hitgnn_cont[i].phi;
      double phi_j = hitgnn_cont[j].phi;

      if( (layer_i<=6 && layer_j<=6) && (layer_i!=layer_j) ){
        gnn_src.emplace_back(i);
        gnn_dst.emplace_back(j);
        gnn_src.emplace_back(j);
        gnn_dst.emplace_back(i);
      }

      if( (layer_i>0 && layer_i<=6 && layer_j>6 && layer_j<19) ||
          (layer_j>0 && layer_j<=6 && layer_i>6 && layer_i<19) ){
        gnn_src.emplace_back(i);
        gnn_dst.emplace_back(j);
        gnn_src.emplace_back(j);
        gnn_dst.emplace_back(i);
      }

      if( (layer_i>6 && layer_i<=24 && layer_j>6 && layer_j<=24) &&
          (layer_i!=layer_j && (abs(phi_i-phi_j)<TMath::Pi()/4 || abs(phi_i-phi_j)>TMath::Pi()*7/4)) ){
        gnn_src.emplace_back(i);
        gnn_dst.emplace_back(j);
        gnn_src.emplace_back(j);
        gnn_dst.emplace_back(i);
      }

      if( ( (layer_i==25 && layer_j>6 && layer_j<=23) || (layer_j==25 && layer_i>6 && layer_i<=23) ) &&
          (layer_i!=layer_j && (abs(phi_i-phi_j)<TMath::Pi()/2 or abs(phi_i-phi_j)>TMath::Pi()*3/2) ) ){
        gnn_src.emplace_back(i);
        gnn_dst.emplace_back(j);
        gnn_src.emplace_back(j);
        gnn_dst.emplace_back(i);
      }

    }
  }

  torch::Tensor edge_index_buf = torch::full({2, (int)gnn_src.size()}, -999, torch::TensorOptions().dtype(torch::kInt));
  for(int i=0; i<(int)gnn_src.size(); ++i){
    edge_index_buf[0][i] = gnn_src[i];
    edge_index_buf[1][i] = gnn_dst[i];
  }

  if(edge_index_buf.size(1)==0) return -1;
  //std::cout << "input_node : " << input_node << std::endl;
  //std::cout << "edge_index : " << edge_index_buf << std::endl;

  std::vector<torch::jit::IValue> inputs;
  inputs.push_back(input_node);
  inputs.push_back(edge_index_buf);

  auto gnn_output = model.forward(inputs).toTuple();
  auto gnn_output_n = gnn_output->elements()[0].toTensor();
  auto gnn_output_e = gnn_output->elements()[1].toTensor();

  att._logger->debug("TGNNFinder: GNN output");

  //std::cout << "gnn_output_n : " << gnn_output_n << std::endl;
  //std::cout << "gnn_output_e : " << gnn_output_e << std::endl;

  auto gnn_pred_e_buf = at::softmax(gnn_output_e, 1).index({torch::indexing::Slice(), 1});

  //std::cout << "gnn_pred_e_buf : " << gnn_pred_e_buf << std::endl;

  // bidir
  torch::Tensor edge_index = torch::full({2, (int)gnn_src.size()/2}, -999, torch::TensorOptions().dtype(torch::kInt));
  for(int i=0; i<(int)gnn_src.size(); ++i){
    if(i%2!=0) continue;
    edge_index[0][i/2] = gnn_src[i];
    edge_index[1][i/2] = gnn_dst[i];
  }

  auto gnn_pred_e = torch::full({(int)gnn_src.size()/2, 1}, -999, torch::TensorOptions().dtype(torch::kFloat));
  for(int i=0; i<(int)gnn_src.size()/2; ++i){
    gnn_pred_e[i] = (gnn_pred_e_buf[i*2] + gnn_pred_e_buf[i*2 +1])/2.;
  }

  auto gnn_pred_n = at::softmax(gnn_output_n, 1).index({torch::indexing::Slice(), 1});

  //std::cout << "gnn_pred_n : " << gnn_pred_n << std::endl;


  // clustering
  double threshold_n = 0.5;
  double threshold_e = 0.5;

  torch::Tensor gnn_label_e = gnn_pred_e > 100;
  torch::Tensor gnn_label_n = gnn_pred_n > threshold_n;

  torch::Tensor gnn_label_g_bce = torch::full({(int)gnn_src.size()/2, 1}, 0, torch::TensorOptions().dtype(torch::kFloat));

  std::tuple< std::map<int, std::set<int> >, torch::Tensor, torch::Tensor, bool > ret = get_label_g(
      input_node, edge_index, gnn_label_n, gnn_label_e, gnn_label_g, gnn_label_g_bce, gnn_pred_e);

  gnn_label_g     = std::get<0>(ret);
  gnn_label_e     = std::get<1>(ret);
  gnn_label_g_bce = std::get<2>(ret);

  double gloss_last = bce_loss(gnn_pred_e, gnn_label_g_bce).item<double>();
  double gloss_first = gloss_last;

  auto gnn_pred_order = gnn_pred_e.argsort(0, true);

  std::vector<int> gnn_list_ov2;


  for(int i=0; i<(int)gnn_pred_order.size(0); ++i ){
    int id   = gnn_pred_order[i].item<int>();
    int id_s = edge_index[0][id].item<int>();
    int id_d = edge_index[1][id].item<int>();
    if(gnn_pred_e[id].item<double>() < threshold_e) break;
    //if(par->flag_debug_gnn) std::cout << Form("- %3d <-> %3d : pred %12.5e",
    //    edge_index[0][id].item<int>(), edge_index[1][id].item<int>(), gnn_pred_e[id].item<double>()) << std::endl;

    if(gnn_label_n[id_s].item<int>() && gnn_label_n[id_d].item<int>()){
      gnn_list_ov2.emplace_back(id);
      //if(par->flag_debug_gnn) std::cout << "- ov2" << std::endl;
      continue;
    }

    ret = get_label_g(input_node, edge_index, gnn_label_n, gnn_label_e, gnn_label_g, gnn_label_g_bce, gnn_pred_e, id, id_s, id_d);
    std::map<int, std::set<int> > gnn_label_g_buf     = std::get<0>(ret);
    torch::Tensor                 gnn_label_e_buf     = std::get<1>(ret);
    torch::Tensor                 gnn_label_g_bce_buf = std::get<2>(ret);
    bool                          gnn_flag_in         = std::get<3>(ret);

    if(gnn_flag_in){
      //if(par->flag_debug_gnn) std::cout << "- included" << std::endl;
      gnn_label_g     = gnn_label_g_buf;
      gnn_label_e     = gnn_label_e_buf;
      gnn_label_g_bce = gnn_label_g_bce_buf;
      continue;
    }

    double gloss = bce_loss(gnn_pred_e, gnn_label_g_bce_buf).item<double>();

    if(gloss < gloss_last){
      //if(par->flag_debug_gnn) std::cout << Form("- update : %12.5e -> %12.5e", gloss_last, gloss) << std::endl;
      gloss_last = gloss;
      gnn_label_g     = gnn_label_g_buf;
      gnn_label_e     = gnn_label_e_buf;
      gnn_label_g_bce = gnn_label_g_bce_buf;
    }

  }

  //if(par->flag_debug_gnn) std::cout << "- gnn_label_g 1st\n" << gnn_label_g << std::endl;

  for(auto id: gnn_list_ov2){
    int id_s = edge_index[0][id].item<int>();
    int id_d = edge_index[1][id].item<int>();
    if(gnn_pred_e[id].item<double>() < threshold_e) break;
    //if(par->flag_debug_gnn) std::cout << Form("- [ov2] %3d <-> %3d : pred %12.5e",
    //    edge_index[0][id].item<int>(), edge_index[1][id].item<int>(), gnn_pred_e[id].item<double>()) << std::endl;

    std::set<int> list_s;
    std::set<int> list_d;
    for(auto l: gnn_label_g){
      if(l.second.size()>0 && l.second.find(id_s)!=l.second.end()) list_s.insert(l.first);
      if(l.second.size()>0 && l.second.find(id_d)!=l.second.end()) list_d.insert(l.first);
    }

    for(auto i: list_s){
      ret = get_label_g(input_node, edge_index, gnn_label_n, gnn_label_e, gnn_label_g, gnn_label_g_bce, gnn_pred_e, id, id_s, id_d, i);
      std::map<int, std::set<int> > gnn_label_g_buf     = std::get<0>(ret);
      torch::Tensor                 gnn_label_e_buf     = std::get<1>(ret);
      torch::Tensor                 gnn_label_g_bce_buf = std::get<2>(ret);
      bool                          gnn_flag_in         = std::get<3>(ret);

      if(gnn_flag_in){
        //if(par->flag_debug_gnn) std::cout << "- included" << std::endl;
        gnn_label_g     = gnn_label_g_buf;
        gnn_label_e     = gnn_label_e_buf;
        gnn_label_g_bce = gnn_label_g_bce_buf;
        continue;
      }

      double gloss = bce_loss(gnn_pred_e, gnn_label_g_bce_buf).item<double>();

      if(gloss < gloss_last){
        //if(par->flag_debug_gnn) std::cout << Form("- update : %12.5e -> %12.5e", gloss_last, gloss) << std::endl;
        gloss_last = gloss;
        gnn_label_g     = gnn_label_g_buf;
        gnn_label_e     = gnn_label_e_buf;
        gnn_label_g_bce = gnn_label_g_bce_buf;
      }
    }

    for(auto i: list_d){
      ret = get_label_g(input_node, edge_index, gnn_label_n, gnn_label_e, gnn_label_g, gnn_label_g_bce, gnn_pred_e, id, id_d, id_s, i);
      std::map<int, std::set<int> > gnn_label_g_buf     = std::get<0>(ret);
      torch::Tensor                 gnn_label_e_buf     = std::get<1>(ret);
      torch::Tensor                 gnn_label_g_bce_buf = std::get<2>(ret);
      bool                          gnn_flag_in         = std::get<3>(ret);

      if(gnn_flag_in){
        //if(par->flag_debug_gnn) std::cout << "- included" << std::endl;
        gnn_label_g     = gnn_label_g_buf;
        gnn_label_e     = gnn_label_e_buf;
        gnn_label_g_bce = gnn_label_g_bce_buf;
        continue;
      }

      double gloss = bce_loss(gnn_pred_e, gnn_label_g_bce_buf).item<double>();

      if(gloss < gloss_last){
        //if(par->flag_debug_gnn) std::cout << Form("- update : %12.5e -> %12.5e", gloss_last, gloss) << std::endl;
        gloss_last = gloss;
        gnn_label_g     = gnn_label_g_buf;
        gnn_label_e     = gnn_label_e_buf;
        gnn_label_g_bce = gnn_label_g_bce_buf;
      }
    }

  }

  //if(par->flag_debug_gnn) std::cout << "- gnn_label_g 2nd\n" << gnn_label_g << std::endl;

  // delete inclusion
  std::map<int, std::set<int> > gnn_label_g_buf = gnn_label_g;
  for(auto g1: gnn_label_g){
    for(auto g2: gnn_label_g){
      if(g1.first!=g2.first && gnn_label_g_buf.find(g1.first)!=gnn_label_g_buf.end() && gnn_label_g_buf.find(g2.first)!=gnn_label_g_buf.end() ){
        bool flag = true;
        for(auto i: g1.second){
          if(g2.second.find(i)==g2.second.end()) flag = false;
        }
        if(flag){
          //if(par->flag_debug_gnn) std::cout << Form("- delete inclusion : %d", g1.first) << std::endl;
          gnn_label_g_buf.erase(g1.first);
        }
      }
    }
  }
  gnn_label_g = gnn_label_g_buf;

  att._logger->debug("TGNNFinder: gloss {:.5e} -> {:.5e}",gloss_first, gloss_last);
  for(auto g: gnn_label_g){
    std::string group_buf;
    std::for_each(g.second.begin(), g.second.end(), [&group_buf](int num) {
        if (!group_buf.empty()){ group_buf += ", ";}
        group_buf += std::to_string(num);
        });
    att._logger->debug("TGNNFinder: {:3d} : {}",g.first, group_buf);
  }

  att._logger->debug("TGNNFinder: TrackHit start");
  std::vector<TrackHit*> TrackHitCont;

  for(auto g: gnn_label_g){

    int num_fiber = 0;
    int num_mdc   = 0;
    bool flag_psb  = false;
    bool flag_psfe = false;
    for(auto id: g.second){
      int layer = hitgnn_cont[id].did;
      if(1<=layer && layer<=6 ) ++num_fiber;
      if(7<=layer && layer<=23) ++num_mdc;
      if(layer==24) flag_psb  = true;
      if(layer==25) flag_psfe = true;
    }
    //att._logger->debug("TGNNFinder: Fiber:{}, MDC:{}, PSB:{}, PSFE:{}",num_fiber, num_mdc, flag_psb, flag_psfe);
    if(!flag_psb && !flag_psfe) continue;
    if(num_fiber <= 1 ) continue;
    if(num_mdc   <= 3 ) continue;
    att._logger->debug("TGNNFinder: Fiber:{}, MDC:{}, PSB:{}, PSFE:{}",num_fiber, num_mdc, flag_psb, flag_psfe);

    TrackHit *track_hit = new TrackHit();
    PSBHitAna *psb_hit = nullptr;
    PSFEHitAna *psfe_hit = nullptr;
    for(auto id: g.second){
      int layer = hitgnn_cont[id].did;
      int hitid = hitgnn_cont[id].lid;

      if(1<=layer && layer<=6){
        int det = (layer-1)/3 + 3;
        int lay = (layer-1)%3;
        for(auto hit: RecoEvent.FiberHitClCont[det][lay]){
          if(fabs(hit->GetFib()-hitid)<1.5){
            track_hit->AddFiber(layer-1, hit->GetClFib());
          }
        }
      }

      if(7<=layer && layer<=23){
        int lay = layer-7;
        for(auto hit: RecoEvent.MDCHitCont[lay]){
          double dif = hit->GetWir() - hitid;
          if(fabs(dif) < 2){
            hit->SetDif(dif);
            track_hit->SetMDCdif(lay, hit->GetWir(), dif);
          }
        }
      }

      if(layer==24){
        double r = 1000;
        for(auto hit: RecoEvent.PSBHitCont){
          if( (fabs(hit->GetSeg() - hitid) < 0.6) && (hit->GetR() < r) ){
            psb_hit = hit;
            r = hit->GetR();
          }
        }
        track_hit->AddPSB(psb_hit->GetSeg());
        track_hit->SetFlagPSB();
      }

      if(layer==25){
        double dif_min = 1000;
        for(auto hit: RecoEvent.PSFEHitCont){
          double dif = fabs(hit->GetSeg() - hitid);
          if( (fabs(dif) < 0.6) && (dif < dif_min) ){
            psfe_hit = hit;
            dif_min = dif;
          }
        }
        track_hit->AddPSFE(psfe_hit->GetSeg());
        track_hit->SetFlagPSFE();
      }

    }

    double x1 = att.Target_PositionX;
    double y1 = att.Target_PositionY;
    double z1 = att.Target_PositionZ;
    double x2 = -9999.;
    double y2 = -9999.;
    double z2 = -9999.;
    if(track_hit->IsFlagPSB()){
      x2 = ( psb_hit->GetR() * cos(psb_hit->GetPhi()) + att.psb_pos_x ) * mm2cm;
      y2 = ( psb_hit->GetR() * sin(psb_hit->GetPhi()) + att.psb_pos_y ) * mm2cm;
      z2 = ( psb_hit->GetZ()                          + att.psb_pos_z ) * mm2cm;
    }
    else if(track_hit->IsFlagPSFE()){
      x2 = (psfe_hit->GetRmax() + psfe_hit->GetRmin())/2. * cos(psfe_hit->GetPhi()) * mm2cm;
      y2 = (psfe_hit->GetRmax() + psfe_hit->GetRmin())/2. * sin(psfe_hit->GetPhi()) * mm2cm;
      z2 = psfe_hit->GetZ()                                                         * mm2cm;
    }

    double a_buf = (x2 - x1) / (z2 - z1);
    double b_buf = (y2 - y1) / (z2 - z1);
    track_hit->SetTrackA(a_buf);
    track_hit->SetTrackB(b_buf);
    //std::cout << "x1 : " << x1 << std::endl;
    //std::cout << "y1 : " << y1 << std::endl;
    //std::cout << "z1 : " << z1 << std::endl;
    //std::cout << "x2 : " << x2 << std::endl;
    //std::cout << "y2 : " << y2 << std::endl;
    //std::cout << "z2 : " << z2 << std::endl;
    //std::cout << "a : " << a_buf << std::endl;
    //std::cout << "b : " << b_buf << std::endl;

    //track_hit->SetMDCLayHitCont();
    track_hit->DeleteDupMDC();
    //track_hit->SetDidCh();


    TrackHitCont.emplace_back(track_hit);
  }

  att._logger->debug("TGNNFinder: TrackHit {}", (int)TrackHitCont.size());


  for(size_t i = 0; i < TrackHitCont.size(); ++i)
  {
    std::vector<int> tempSetHit(G4Sol::SIZEOF_G4SOLDETTYPE, -1);
    std::vector<InfoPar> tempSetInfo(G4Sol::SIZEOF_G4SOLDETTYPE);

    //FiberHits
    for(int j = 0; j < 6; ++j)
    {
      if(TrackHitCont[i]->GetFiberHit()[j] == -1) continue;

      bool flag_found = false;
      for(int k = 0; k < (int)RecoEvent.ListHits[G4Sol::MiniFiberD1_x + j].size(); ++k)
      {
        //printf("TrackHit: %d ; ListHists: %d \n", TrackHitCont[i]->GetFiberHit()[j], RecoEvent.ListHits[G4Sol::MiniFiberD1_x + j][k]->getHitId());
        if(TrackHitCont[i]->GetFiberHit()[j] == RecoEvent.ListHits[G4Sol::MiniFiberD1_x + j][k]->getHitId())
        {
          tempSetHit[G4Sol::MiniFiberD1_x + j] = k;
          InfoPar tmp_infopar;
          tmp_infopar.pdg    = -211;
          //tmp_infopar.momX   = hit.MomX;
          //tmp_infopar.momY   = hit.MomY;
          //tmp_infopar.momZ   = hit.MomZ;
          //tmp_infopar.mass   = hit.Mass;
          tmp_infopar.Eloss  = RecoEvent.ListHitsInfo[G4Sol::MiniFiberD1_x + j][k].dE;
          tmp_infopar.time   = RecoEvent.ListHitsInfo[G4Sol::MiniFiberD1_x + j][k].time;
          tmp_infopar.TOT    = RecoEvent.ListHitsInfo[G4Sol::MiniFiberD1_x + j][k].TOT;
          //tmp_infopar.length = hit.TrackLength;
          tempSetInfo[G4Sol::MiniFiberD1_x + j] = tmp_infopar;
          flag_found = true;
        }
      }

      if(!flag_found) printf("Error: %d-layer Fiber Hit not found in GnnFinder\n", j);
    }

    //MDCHits
    for(int j = 0; j < 17; ++j)
    {
      if(TrackHitCont[i]->GetMDCHit()[j] == -1) continue;

      bool flag_found = false;
      for(int k = 0; k < (int)RecoEvent.ListHits[G4Sol::MG01 + j].size(); ++k)
      {
        if(TrackHitCont[i]->GetMDCHit()[j] == RecoEvent.ListHits[G4Sol::MG01 + j][k]->getHitId())
        {
          tempSetHit[G4Sol::MG01 + j] = k;
          InfoPar tmp_infopar;
          tmp_infopar.pdg    = -211;
          //tmp_infopar.momX;
          //tmp_infopar.momY;
          //tmp_infopar.momZ;
          //tmp_infopar.mass;
          tmp_infopar.Eloss  = RecoEvent.ListHitsInfo[G4Sol::MG01 + j][k].dE;
          tmp_infopar.time   = RecoEvent.ListHitsInfo[G4Sol::MG01 + j][k].time;
          tmp_infopar.TOT    = RecoEvent.ListHitsInfo[G4Sol::MG01 + j][k].TOT;
          //tmp_infopar.length = hit.TrackLength;
          tempSetInfo[G4Sol::MG01 + j] = tmp_infopar;
          flag_found = true;
        }
      }

      if(!flag_found)
        printf("Error: %d-layer MDC Hit not found in GnnFinder\n", j);
    }

    //PSBHit
    if(att.WF_PSBHits)
      if(TrackHitCont[i]->IsFlagPSB())
      {
        bool flag_found = false;
        for(int j = 0; j < (int)RecoEvent.ListHits[G4Sol::PSCE].size(); ++j)
        {
          if(TrackHitCont[i]->GetPSBHit() == RecoEvent.ListHits[G4Sol::PSCE][j]->getHitId())
          {
            tempSetHit[G4Sol::PSCE] = j;
            InfoPar tmp_infopar;
            tmp_infopar.pdg    = -211;
            //tmp_infopar.momX   = hit.MomX;
            //tmp_infopar.momY   = hit.MomY;
            //tmp_infopar.momZ   = hit.MomZ;
            //tmp_infopar.mass   = hit.Mass;
            tmp_infopar.Eloss  = RecoEvent.ListHitsInfo[G4Sol::PSCE][j].dE;
            tmp_infopar.time   = RecoEvent.ListHitsInfo[G4Sol::PSCE][j].time;
            //tmp_infopar.TOT    = RecoEvent.ListHitsInfo[G4Sol::PSCE][j].TOT;
            //tmp_infopar.length = hit.TrackLength;
            tempSetInfo[G4Sol::PSCE] = tmp_infopar;
            flag_found = true;
          }
        }

        if(!flag_found)
          std::cout << "Error: PSB Hit not found in GnnFinder\n";
      }

    //PSFEHit
    if(att.WF_PSFEHits)
      if(TrackHitCont[i]->IsFlagPSFE())
      {
        bool flag_found = false;
        for(int j = 0; j < (int)RecoEvent.ListHits[G4Sol::PSFE].size(); ++j)
        {
          if(TrackHitCont[i]->GetPSFEHit() == RecoEvent.ListHits[G4Sol::PSFE][j]->getHitId())
          {
            tempSetHit[G4Sol::PSFE] = j;
            InfoPar tmp_infopar;
            tmp_infopar.pdg    = -211;
            //tmp_infopar.momX   = hit.MomX;
            //tmp_infopar.momY   = hit.MomY;
            //tmp_infopar.momZ   = hit.MomZ;
            //tmp_infopar.mass   = hit.Mass;
            tmp_infopar.Eloss  = RecoEvent.ListHitsInfo[G4Sol::PSFE][j].dE;
            tmp_infopar.time   = RecoEvent.ListHitsInfo[G4Sol::PSFE][j].time;
            //tmp_infopar.TOT    = RecoEvent.ListHitsInfo[G4Sol::PSFE][j].TOT;
            //tmp_infopar.length = hit.TrackLength;
            tempSetInfo[G4Sol::PSFE] = tmp_infopar;
            flag_found = true;
          }
        }

        if(!flag_found)
          std::cout << "Error: PSFE Hit not found in GnnFinder\n";
      }

    RecoEvent.TrackDAF.insert(std::make_pair(i, tempSetHit));
    RecoEvent.TrackInfo.insert(std::make_pair(i, tempSetInfo));

    //TVector3 tmp_closepoint_Pos;
    //double tmp_closepoint_Dist;
    //if(!RecoEvent.FragmentTracks.empty())
    //  CloseDist(RecoEvent.FragmentTracks[0], TrackHitCont[i], tmp_closepoint_Dist, tmp_closepoint_Pos); //CHECK!!

    InfoInit tempInit;
    tempInit.charge = -1;
    //tempInit.posX = tmp_closepoint_Pos.X();
    //tempInit.posY = tmp_closepoint_Pos.Y();
    //tempInit.posZ = tmp_closepoint_Pos.Z();
    tempInit.posX = att.Target_PositionX;
    tempInit.posY = att.Target_PositionY;
    tempInit.posZ = att.Target_PositionZ;

    //std::cout << "InitPoint: " << tempInit.posX << "  " << tempInit.posY << "  " << tempInit.posZ << "\n";
    double tmp_momZ = 1.;
    tempInit.momX = tmp_momZ*TrackHitCont[i]->GetTrackA();
    tempInit.momY = tmp_momZ*TrackHitCont[i]->GetTrackB();
    tempInit.momZ = tmp_momZ;

    RecoEvent.TrackDAFInit.insert(std::make_pair(i, tempInit));
  }

  att._logger->debug("TGNNFinder: FinderGNN end");
  return 0;
}


template<class Out>
int TGNNFinder<Out>::FinderGNNText(FullRecoEvent& RecoEvent)
{
  att._logger->debug("TGNNFinder: FinderGNNText start");
  std::vector<TrackHit*> TrackHitCont;

  //std::cout << "Text : " << att.GNN_Node << std::endl;
  if(!ifs_gnn.is_open()){
    std::cout << "GNN clust file not open" << std::endl;
    return -1;
  }

  std::vector< std::map<int, double> > node_clust;

  std::string temp_line;
  std::getline(ifs_gnn,temp_line);
  if(ifs_gnn.eof()){
    std::cout << "GNN file end" << std::endl;
    //std::ofstream ofs_err(gnn_err_filename.c_str());
    //ofs_err << Form("GNN file end :  %s", gnn_err_filename.c_str()) << std::endl;
    //ofs_err << Form("last ev      :  %d", ev) << std::endl;
    return -1;
  }
  std::stringstream stream(temp_line);
  int num_cl;
  int num_hit;
  int    layer;
  double hitid;
  stream >> num_cl;
  att._logger->debug("TGNNFinder: num_cl {}",num_cl);
  for(int i=0; i<num_cl; ++i){
    std::map<int, double> buf_map;
    node_clust.emplace_back(buf_map);
    stream >> num_hit;
    for(int j=0; j<num_hit; ++j){
      stream >> layer >> hitid;
      if(layer==-1){
        att._logger->debug("TGNNFinder: GNN ev {}",(int)hitid);
        //if(par->flag_debug) std::cout << Form("- ev : %d %.0f",ev, hitid) <<std::endl;
        //if(hitid!=ev){
        //  std::cout << Form("Event ID mismatch : %d %.0f",ev, hitid) << std::endl;
        //  std::ofstream ofs_err(gnn_err_filename.c_str());
        //  ofs_err << Form("Event ID mismatch : ev %d  gnn %.0f",ev, hitid) << std::endl;
        //  return -1;
        //}
      }
      else (node_clust.back())[layer] = hitid;
    }
  }

  const double mm2cm = 0.1;
  for(auto cl: node_clust){
    int num_fiber = 0;
    int num_mdc   = 0;
    bool flag_psb  = false;
    bool flag_psfe = false;
    for(auto v: cl){
      int layer = v.first;
      if(1<=layer && layer<=6 ) ++num_fiber;
      if(7<=layer && layer<=23) ++num_mdc;
      if(layer==24) flag_psb  = true;
      if(layer==25) flag_psfe = true;
    }
    if(!flag_psb && !flag_psfe) continue;
    if(num_fiber <= 1) continue;
    if(num_mdc   <= 3) continue;
    att._logger->debug("TGNNFinder: Fiber {}, MDC {}, PSB {} PSFE {}",num_fiber, num_mdc, flag_psb, flag_psfe);

    TrackHit *track_hit = new TrackHit();
    PSBHitAna *psb_hit = nullptr;
    PSFEHitAna *psfe_hit = nullptr;
    for(auto v: cl){
      int layer = v.first;
      int hitid = v.second;

      if(1<=layer && layer<=6){
        int det = (layer-1)/3 + 3;
        int lay = (layer-1)%3;
        for(auto hit: RecoEvent.FiberHitClCont[det][lay]){
          if(fabs(hit->GetFib()-hitid)<1.5){
            track_hit->AddFiber(layer-1, hit->GetClFib());
          }
        }
      }

      if(7<=layer && layer<=23){
        int lay = layer-7;
        for(auto hit: RecoEvent.MDCHitCont[lay]){
          double dif = hit->GetWir() - hitid;
          if(fabs(dif) < 2){
            hit->SetDif(dif);
            track_hit->SetMDCdif(lay, hit->GetWir(), dif);
          }
        }
      }

      if(layer==24){
        double r = 1000;
        for(auto hit: RecoEvent.PSBHitCont){
          if( (fabs(hit->GetSeg() - hitid) < 0.6) && (hit->GetR() < r) ){
            psb_hit = hit;
            r = hit->GetR();
          }
        }
        track_hit->AddPSB(psb_hit->GetSeg());
        track_hit->SetFlagPSB();
      }

      if(layer==25){
        double dif_min = 1000;
        for(auto hit: RecoEvent.PSFEHitCont){
          double dif = fabs(hit->GetSeg() - hitid);
          if( (fabs(dif) < 0.6) && (dif < dif_min) ){
            psfe_hit = hit;
            dif_min = dif;
          }
        }
        track_hit->AddPSFE(psfe_hit->GetSeg());
        track_hit->SetFlagPSFE();
      }

    }

    double x1 = att.Target_PositionX;
    double y1 = att.Target_PositionY;
    double z1 = att.Target_PositionZ;
    double x2 = -9999.;
    double y2 = -9999.;
    double z2 = -9999.;
    if(track_hit->IsFlagPSB()){
      x2 = ( psb_hit->GetR() * cos(psb_hit->GetPhi()) + att.psb_pos_x ) * mm2cm;
      y2 = ( psb_hit->GetR() * sin(psb_hit->GetPhi()) + att.psb_pos_y ) * mm2cm;
      z2 = ( psb_hit->GetZ()                          + att.psb_pos_z ) * mm2cm;
    }
    else if(track_hit->IsFlagPSFE()){
      x2 = (psfe_hit->GetRmax() + psfe_hit->GetRmin())/2. * cos(psfe_hit->GetPhi()) * mm2cm;
      y2 = (psfe_hit->GetRmax() + psfe_hit->GetRmin())/2. * sin(psfe_hit->GetPhi()) * mm2cm;
      z2 = psfe_hit->GetZ()                                                         * mm2cm;
    }

    double a_buf = (x2 - x1) / (z2 - z1);
    double b_buf = (y2 - y1) / (z2 - z1);
    track_hit->SetTrackA(a_buf);
    track_hit->SetTrackB(b_buf);

    //track_hit->SetMDCLayHitCont();
    track_hit->DeleteDupMDC();
    //track_hit->SetDidCh();

    TrackHitCont.emplace_back(track_hit);
  }

  att._logger->debug("TGNNFinder: TrackHit {}", (int)TrackHitCont.size());

  for(size_t i = 0; i < TrackHitCont.size(); ++i)
  {
    std::vector<int> tempSetHit(G4Sol::SIZEOF_G4SOLDETTYPE, -1);
    std::vector<InfoPar> tempSetInfo(G4Sol::SIZEOF_G4SOLDETTYPE);

    //FiberHits
    for(int j = 0; j < 6; ++j)
    {
      if(TrackHitCont[i]->GetFiberHit()[j] == -1) continue;

      bool flag_found = false;
      for(int k = 0; k < (int)RecoEvent.ListHits[G4Sol::MiniFiberD1_x + j].size(); ++k)
      {
        //printf("TrackHit: %d ; ListHists: %d \n", TrackHitCont[i]->GetFiberHit()[j], RecoEvent.ListHits[G4Sol::MiniFiberD1_x + j][k]->getHitId());
        if(TrackHitCont[i]->GetFiberHit()[j] == RecoEvent.ListHits[G4Sol::MiniFiberD1_x + j][k]->getHitId())
        {
          tempSetHit[G4Sol::MiniFiberD1_x + j] = k;
          InfoPar tmp_infopar;
          tmp_infopar.pdg    = -211;
          //tmp_infopar.momX   = hit.MomX;
          //tmp_infopar.momY   = hit.MomY;
          //tmp_infopar.momZ   = hit.MomZ;
          //tmp_infopar.mass   = hit.Mass;
          tmp_infopar.Eloss  = RecoEvent.ListHitsInfo[G4Sol::MiniFiberD1_x + j][k].dE;
          tmp_infopar.time   = RecoEvent.ListHitsInfo[G4Sol::MiniFiberD1_x + j][k].time;
          tmp_infopar.TOT    = RecoEvent.ListHitsInfo[G4Sol::MiniFiberD1_x + j][k].TOT;
          //tmp_infopar.length = hit.TrackLength;
          tempSetInfo[G4Sol::MiniFiberD1_x + j] = tmp_infopar;
          flag_found = true;
        }
      }

      if(!flag_found) printf("Error: %d-layer Fiber Hit not found in GnnFinder\n", j);
    }

    //MDCHits
    for(int j = 0; j < 17; ++j)
    {
      if(TrackHitCont[i]->GetMDCHit()[j] == -1) continue;

      bool flag_found = false;
      for(int k = 0; k < (int)RecoEvent.ListHits[G4Sol::MG01 + j].size(); ++k)
      {
        if(TrackHitCont[i]->GetMDCHit()[j] == RecoEvent.ListHits[G4Sol::MG01 + j][k]->getHitId())
        {
          tempSetHit[G4Sol::MG01 + j] = k;
          InfoPar tmp_infopar;
          tmp_infopar.pdg    = -211;
          //tmp_infopar.momX;
          //tmp_infopar.momY;
          //tmp_infopar.momZ;
          //tmp_infopar.mass;
          tmp_infopar.Eloss  = RecoEvent.ListHitsInfo[G4Sol::MG01 + j][k].dE;
          tmp_infopar.time   = RecoEvent.ListHitsInfo[G4Sol::MG01 + j][k].time;
          tmp_infopar.TOT    = RecoEvent.ListHitsInfo[G4Sol::MG01 + j][k].TOT;
          //tmp_infopar.length = hit.TrackLength;
          tempSetInfo[G4Sol::MG01 + j] = tmp_infopar;
          flag_found = true;
        }
      }

      if(!flag_found)
        printf("Error: %d-layer MDC Hit not found in GnnFinder\n", j);
    }

    //PSBHit
    if(att.WF_PSBHits)
      if(TrackHitCont[i]->IsFlagPSB())
      {
        bool flag_found = false;
        for(int j = 0; j < (int)RecoEvent.ListHits[G4Sol::PSCE].size(); ++j)
        {
          if(TrackHitCont[i]->GetPSBHit() == RecoEvent.ListHits[G4Sol::PSCE][j]->getHitId())
          {
            tempSetHit[G4Sol::PSCE] = j;
            InfoPar tmp_infopar;
            tmp_infopar.pdg    = -211;
            //tmp_infopar.momX   = hit.MomX;
            //tmp_infopar.momY   = hit.MomY;
            //tmp_infopar.momZ   = hit.MomZ;
            //tmp_infopar.mass   = hit.Mass;
            tmp_infopar.Eloss  = RecoEvent.ListHitsInfo[G4Sol::PSCE][j].dE;
            tmp_infopar.time   = RecoEvent.ListHitsInfo[G4Sol::PSCE][j].time;
            //tmp_infopar.TOT    = RecoEvent.ListHitsInfo[G4Sol::PSCE][j].TOT;
            //tmp_infopar.length = hit.TrackLength;
            tempSetInfo[G4Sol::PSCE] = tmp_infopar;
            flag_found = true;
          }
        }

        if(!flag_found)
          std::cout << "Error: PSB Hit not found in GnnFinder\n";
      }

    //PSFEHit
    if(att.WF_PSFEHits)
      if(TrackHitCont[i]->IsFlagPSFE())
      {
        bool flag_found = false;
        for(int j = 0; j < (int)RecoEvent.ListHits[G4Sol::PSFE].size(); ++j)
        {
          if(TrackHitCont[i]->GetPSFEHit() == RecoEvent.ListHits[G4Sol::PSFE][j]->getHitId())
          {
            tempSetHit[G4Sol::PSFE] = j;
            InfoPar tmp_infopar;
            tmp_infopar.pdg    = -211;
            //tmp_infopar.momX   = hit.MomX;
            //tmp_infopar.momY   = hit.MomY;
            //tmp_infopar.momZ   = hit.MomZ;
            //tmp_infopar.mass   = hit.Mass;
            tmp_infopar.Eloss  = RecoEvent.ListHitsInfo[G4Sol::PSFE][j].dE;
            tmp_infopar.time   = RecoEvent.ListHitsInfo[G4Sol::PSFE][j].time;
            //tmp_infopar.TOT    = RecoEvent.ListHitsInfo[G4Sol::PSFE][j].TOT;
            //tmp_infopar.length = hit.TrackLength;
            tempSetInfo[G4Sol::PSFE] = tmp_infopar;
            flag_found = true;
          }
        }

        if(!flag_found)
          std::cout << "Error: PSFE Hit not found in GnnFinder\n";
      }

    RecoEvent.TrackDAF.insert(std::make_pair(i, tempSetHit));
    RecoEvent.TrackInfo.insert(std::make_pair(i, tempSetInfo));

    //TVector3 tmp_closepoint_Pos;
    //double tmp_closepoint_Dist;
    //if(!RecoEvent.FragmentTracks.empty())
    //  CloseDist(RecoEvent.FragmentTracks[0], TrackHitCont[i], tmp_closepoint_Dist, tmp_closepoint_Pos); //CHECK!!

    InfoInit tempInit;
    tempInit.charge = -1;
    //tempInit.posX = tmp_closepoint_Pos.X();
    //tempInit.posY = tmp_closepoint_Pos.Y();
    //tempInit.posZ = tmp_closepoint_Pos.Z();
    tempInit.posX = att.Target_PositionX;
    tempInit.posY = att.Target_PositionY;
    tempInit.posZ = att.Target_PositionZ;

    //std::cout << "InitPoint: " << tempInit.posX << "  " << tempInit.posY << "  " << tempInit.posZ << "\n";
    double tmp_momZ = 1.;
    tempInit.momX = tmp_momZ*TrackHitCont[i]->GetTrackA();
    tempInit.momY = tmp_momZ*TrackHitCont[i]->GetTrackB();
    tempInit.momZ = tmp_momZ;

    RecoEvent.TrackDAFInit.insert(std::make_pair(i, tempInit));
  }



  att._logger->debug("TGNNFinder: FinderGNNText end");
  return 0;
}


template<class Out>
HitGnn TGNNFinder<Out>::make_cluster(HitGnn hit_a, HitGnn hit_b, double lid, double s, double dl, double phi, double wang1, double wang2){

  int clsize_sum = hit_a.clsize + hit_b.clsize;

  HitGnn hit_buf;
  hit_buf.did     = hit_a.did;
  hit_buf.lid     = lid;
  hit_buf.phi     = phi;
  hit_buf.wang1   = wang1;
  hit_buf.wang2   = wang2;
  hit_buf.r       = (hit_a.r * hit_a.clsize) + (hit_b.r * hit_b.clsize) / clsize_sum;
  hit_buf.z       = (hit_a.z * hit_a.clsize) + (hit_b.z * hit_b.clsize) / clsize_sum;
  hit_buf.s       = s;
  hit_buf.dl      = dl;
  hit_buf.dle     = hit_a.dle;
  hit_buf.de      = hit_a.de + hit_b.de;
  hit_buf.clsize  = clsize_sum;
  hit_buf.cls     = hit_a.cls;
  hit_buf.cle     = hit_b.cle;

  return hit_buf;

}

template<class Out>
bool  TGNNFinder<Out>::not_dup(torch::Tensor input_node, std::set<int> label_g1, std::set<int> label_g2){
  bool flag_dup = true;
  bool flag_psb = false;
  bool flag_psfe = false;
  for(auto l1: label_g1){
    for(auto l2: label_g2){
      if( (l1 != l2) && (input_node[l1][0].item<int>() == input_node[l2][0]).item<int>()  ) flag_dup = false;
      if( (input_node[l1][0].item<int>()==24) || (input_node[l2][0].item<int>()==24) ) flag_psb  = true;
      if( (input_node[l1][0].item<int>()==25) || (input_node[l2][0].item<int>()==25) ) flag_psfe = true;
    }
  }
  if(flag_psb && flag_psfe) flag_dup = false;
  //if(par->flag_debug_gnn && !flag_dup){
  //  std::cout << "- dup" << std::endl;
  //  std::cout << "label_g1 : " << label_g1 << std::endl;
  //  std::cout << "label_g2 : " << label_g2 << std::endl;
  //}
  return flag_dup;
}

template<class Out>
std::set<int> TGNNFinder<Out>::get_edge_id_g(torch::Tensor edge_index, std::set<int> label_g){
  std::set<int> id_g;
  for(int i=0; i<(int)edge_index.size(1); ++i){
    int id_s = edge_index[0][i].item<int>();
    int id_d = edge_index[1][i].item<int>();
    if(label_g.size()>0 && label_g.find(id_s)!=label_g.end() && label_g.find(id_d)!=label_g.end()) id_g.insert(i);
  }
  return id_g;
}


template<class Out>
std::tuple< std::map<int, std::set<int> >, torch::Tensor, torch::Tensor, bool > TGNNFinder<Out>::get_label_g(
    torch::Tensor input_node, torch::Tensor edge_index, torch::Tensor gnn_label_n, torch::Tensor gnn_label_e,
    std::map<int, std::set<int> > gnn_label_g, torch::Tensor gnn_label_g_bce, torch::Tensor gnn_pred_e, int id, int id_s, int id_d, int g_s){

  std::tuple< std::map<int, std::set<int> >, torch::Tensor, torch::Tensor, bool > ret;

  torch::Tensor src = edge_index.index({0, torch::indexing::Slice()});
  torch::Tensor dst = edge_index.index({1, torch::indexing::Slice()});

  // initialize
  if(id == -1){
    std::map<int, std::set<int> > gnn_label_g_buf;
    for(int i=0; i<input_node.size(0); ++i){
      std::set<int> s{i};
      gnn_label_g_buf[i] = s;
    }
    torch::Tensor gnn_label_g_bce_buf = torch::zeros( {gnn_label_g_bce.size(0), 1}, torch::TensorOptions().dtype(torch::kFloat) );
    ret = std::make_tuple( gnn_label_g_buf, gnn_label_e, gnn_label_g_bce_buf, false);
    return ret;
  }

  torch::Tensor gnn_label_e_buf;
  torch::Tensor gnn_label_g_bce_buf;
  std::map<int, std::set<int> > gnn_label_g_buf = gnn_label_g;
  bool flag_in = false;

  gnn_label_e_buf = gnn_label_e.clone();
  gnn_label_g_bce_buf = gnn_label_g_bce.clone();

  gnn_label_e_buf[id] = 1;

  if(gnn_label_n[id_s].item<int>()!=1 && gnn_label_n[id_d].item<int>()!=1){
    //if(par->flag_debug_gnn) std::cout<< Form("- 0-0 edge : %3d and %3d", id_s, id_d) << std::endl;
    int id_s_l = -1;
    int id_d_l = -2;
    for(auto l: gnn_label_g){
      if( (l.second).size()>0 && (l.second).find(id_s)!=(l.second).end() ) id_s_l = l.first;
      if( (l.second).size()>0 && (l.second).find(id_d)!=(l.second).end() ) id_d_l = l.first;
    }
    if(id_s_l == id_d_l) flag_in = true;
    else if(id_s_l != id_d_l && not_dup(input_node, gnn_label_g[id_s_l], gnn_label_g[id_d_l]) ){
      int id_small = id_s_l < id_d_l ? id_s_l : id_d_l;
      int id_large = id_s_l < id_d_l ? id_d_l : id_s_l;
      std::set_union(gnn_label_g[id_small].begin(), gnn_label_g[id_small].end(),
          gnn_label_g[id_large].begin(), gnn_label_g[id_large].end(),
          std::inserter(gnn_label_g_buf[id_small], gnn_label_g_buf[id_small].begin()));
      gnn_label_g_buf.erase(id_large);
      std::set<int> id_g = get_edge_id_g(edge_index, gnn_label_g_buf[id_small]);
      for(auto id: id_g){
        gnn_label_g_bce_buf[id] = 1;
      }
    }
  }

  if(gnn_label_n[id_s].item<int>()==1 && gnn_label_n[id_d].item<int>()!=1){
    //if(par->flag_debug_gnn) std::cout<< Form("- 1-0 edge : %3d into %3d", id_s, id_d) << std::endl;
    int id_l = -1;
    for(auto l: gnn_label_g){
      if( (l.second).size()>0 && (l.second).find(id_d)!=(l.second).end() ){
        id_l = l.first;
        break;
      }
    }
    if( gnn_label_g[id_l].size()>0 && gnn_label_g[id_l].find(id_s)!=gnn_label_g[id_l].end() ) flag_in = true;
    if( not_dup(input_node, gnn_label_g[id_s], gnn_label_g[id_l]) ){ // ##### check is_connect
      std::set_union(gnn_label_g[id_l].begin(), gnn_label_g[id_l].end(),
          gnn_label_g[id_s].begin(), gnn_label_g[id_s].end(),
          std::inserter(gnn_label_g_buf[id_l], gnn_label_g_buf[id_l].begin()));
      std::set<int> id_g = get_edge_id_g(edge_index, gnn_label_g_buf[id_l]);
      for(auto id: id_g){
        gnn_label_g_bce_buf[id] = 1;
      }
    }
  }


  if(gnn_label_n[id_d].item<int>()==1 && gnn_label_n[id_s].item<int>()!=1){
    //if(par->flag_debug_gnn) std::cout<< Form("- 1-0 edge : %3d into %3d", id_d, id_s) << std::endl;
    int id_l = -1;
    for(auto l: gnn_label_g){
      if( (l.second).size()>0 && (l.second).find(id_s)!=(l.second).end() ){
        id_l = l.first;
        break;
      }
    }
    if( gnn_label_g[id_l].size()>0 && gnn_label_g[id_l].find(id_d)!=gnn_label_g[id_l].end() ) flag_in = true;
    if( not_dup(input_node, gnn_label_g[id_d], gnn_label_g[id_l]) ){ // ##### check is_connect
      std::set_union(gnn_label_g[id_l].begin(), gnn_label_g[id_l].end(),
          gnn_label_g[id_d].begin(), gnn_label_g[id_d].end(),
          std::inserter(gnn_label_g_buf[id_l], gnn_label_g_buf[id_l].begin()));
      std::set<int> id_g = get_edge_id_g(edge_index, gnn_label_g_buf[id_l]);
      for(auto id: id_g){
        gnn_label_g_bce_buf[id] = 1;
      }
    }
  }

  if(gnn_label_n[id_s].item<int>()==1 && gnn_label_n[id_d].item<int>()==1 && gnn_label_n[g_s].item<int>()!=1){
    //if(par->flag_debug_gnn) std::cout<< Form("- 1-1 edge : %3d into %3d", id_d, g_s) << std::endl;
    if( gnn_label_g[g_s].size()>0 && gnn_label_g[g_s].find(id_d)!=gnn_label_g[g_s].end() ) flag_in = true;
    if( not_dup(input_node, gnn_label_g[g_s], gnn_label_g[id_d]) ){ // ##### check is_connect
      std::set_union(gnn_label_g[g_s].begin(), gnn_label_g[g_s].end(),
          gnn_label_g[id_d].begin(), gnn_label_g[id_d].end(),
          std::inserter(gnn_label_g_buf[g_s], gnn_label_g_buf[g_s].begin()));
      std::set<int> id_g = get_edge_id_g(edge_index, gnn_label_g_buf[g_s]);
      for(auto id: id_g){
        gnn_label_g_bce_buf[id] = 1;
      }
    }
  }

  ret = std::make_tuple( gnn_label_g_buf, gnn_label_e, gnn_label_g_bce_buf, flag_in);

  return ret;
}





template class TGNNFinder<MCAnaEventG4Sol>;
template class TGNNFinder<Ana_WasaEvent>;
