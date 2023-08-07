#include <torch/torch.h>
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
  //try{
  //  model = torch::jit::load("scripted_model/scripted_model_m7_l4_5_h180_480_b25_e50_d040_lr000004_v04_div0005.pt");
  //                           //"scripted_model/scripted_model_m7_l4_5_h180_480_b25_e50_d040_lr000004_v04_div0005.pt"
  //  std::cout << "model is loaded" << std::endl;
  //}
  //catch(const c10::Error& e){
  //  std::cerr << "error loading the model\n";
  //  return;
  //}
  //model.eval();
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
int TGNNFinder<Out>::Exec(FullRecoEvent& RecoEvent, Out* ) { return FinderGNN(RecoEvent); }

template<class Out>
ReturnRes::InfoM TGNNFinder<Out>::SoftExit(int result_full) { return ReturnRes::Fine; }


  template<class Out>
void TGNNFinder<Out>::SelectHists() {}


template<class Out>
int TGNNFinder<Out>::FinderGNN(FullRecoEvent& RecoEvent)
{

  double mm2cm = 0.1;

  torch::Tensor x = torch::full({3, 3}, 1.5, torch::TensorOptions().dtype(torch::kFloat));
  std::cout << x << std::endl;

  torch::jit::script::Module model;
  try{
    model = torch::jit::load("scripted_model/scripted_model_m7_l4_5_h180_480_b25_e50_d040_lr000004_v04_div0005.pt");
                             //"scripted_model/scripted_model_m7_l4_5_h180_480_b25_e50_d040_lr000004_v04_div0005.pt"
    std::cout << "model is loaded" << std::endl;
  }
  catch(const c10::Error& e){
    std::cerr << "error loading the model\n";
    return -1;
  }
  model.eval();

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
  att._logger->debug(Form("TGNNFinder: Fiber hit end %d",hitgnn_cont.size()));

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
  att._logger->debug(Form("TGNNFinder: MDC   hit end %d",hitgnn_cont.size()));


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
  att._logger->debug(Form("TGNNFinder: PSB   hit end %d",hitgnn_cont.size()));


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
  att._logger->debug(Form("TGNNFinder: PSFE  hit end %d",hitgnn_cont.size()));


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

  att._logger->debug(Form("TGNNFinder: delete overlap %d",hitgnn_cont.size()));

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

  att._logger->debug(Form("TGNNFinder: clustering     %d",hitgnn_cont.size()));


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

  std::vector<torch::jit::IValue> inputs;
  inputs.push_back(input_node);
  inputs.push_back(edge_index_buf);

  auto gnn_output = model.forward(inputs).toTuple();
  auto gnn_output_n = gnn_output->elements()[0].toTensor();
  auto gnn_output_e = gnn_output->elements()[1].toTensor();

  auto gnn_pred_e_buf = at::softmax(gnn_output_e, 1).index({torch::indexing::Slice(), 1});

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



  att._logger->debug("TGNNFinder: FinderGNN end");
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



template class TGNNFinder<MCAnaEventG4Sol>;
template class TGNNFinder<Ana_WasaEvent>;
