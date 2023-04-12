#include "ConstantParameter.hh"
#include "FiberAnalyzer.hh"
#include "FiberHitAna.hh"
#include "FiberHitXUV.hh"
#include <math.h>
#include <iostream>
#include <set>
#include <algorithm>
#include "TMath.h"

FiberAnalyzer::FiberAnalyzer(){}
FiberAnalyzer::~FiberAnalyzer() {}


std::vector< std::vector< std::vector< FiberHitAna* > > > FiberAnalyzer::Clusterize(std::vector< std::vector< std::vector< FiberHitAna* > > > &cont){

  std::vector< std::vector< std::vector< FiberHitAna* > > > buf_cont;
  for(int i=0; i<7; ++i){
    std::vector<FiberHitAna*> buf_v;
    std::vector< std::vector<FiberHitAna*> > buf_vv;
    for(int j=0; j<3; ++j){
      buf_vv.emplace_back(buf_v);
    }
    buf_cont.emplace_back(buf_vv);
  }

  for(int i=0; i<7; ++i){
    for(int j=0; j<3; ++j){
      for(int k=0; k<(int)cont[i][j].size(); ++k){
        FiberHitAna *buf_hit = new FiberHitAna(cont[i][j][k]);
        buf_cont[i][j].emplace_back(buf_hit);
        //cont[i][j][k]->Print();
      }
    }
  }

  for(int i=0; i<7; ++i){
    for(int j=0; j<3; ++j){
      if(buf_cont[i][j].size()<2) continue;
      std::set<int> flag_used;
      //std::cout << "\nsize before : " << buf_cont[i][j].size() << std::endl;
      for(size_t k=0; k<buf_cont[i][j].size(); ++k)
        {
          if(flag_used.size()>0 && flag_used.find(k)!=flag_used.end()) continue;
          double pos_cluster = buf_cont[i][j][k]->GetPos();
          double edge_left  = pos_cluster;
          double edge_right = pos_cluster;
          for(size_t l = k+1; l < buf_cont[i][j].size(); ++l)
            {

              //printf("%d %d %.8f\n" , i, j, fabs(buf_cont[i][j][k]->GetPos()-buf_cont[i][j][l]->GetPos()));

              if(flag_used.size()>0 && flag_used.find(l)!=flag_used.end()) continue;

              if( (buf_cont[i][j][l]->GetPos()>=edge_left &&  buf_cont[i][j][l]->GetPos()<=edge_right)
                || fabs(edge_left - buf_cont[i][j][l]->GetPos())<1.1 || fabs(edge_left - buf_cont[i][j][l]->GetPos())-1.1 <1e-4
                || fabs(edge_right - buf_cont[i][j][l]->GetPos())<1.1 || fabs(edge_right - buf_cont[i][j][l]->GetPos())-1.1 <1e-4)
                {
                                //printf("%d %d %.8f %.8f\n" , i, j, fabs(edge_right-buf_cont[i][j][k]->GetPos()), fabs(edge_right-buf_cont[i][j][l]->GetPos()));
                                //printf("%d %d %.8f %.8f\n" , i, j, fabs(edge_left-buf_cont[i][j][l]->GetPos()), fabs(edge_right-buf_cont[i][j][l]->GetPos()));

                  buf_cont[i][j][k]->Add(buf_cont[i][j][l]);
                  flag_used.insert(l);
                  if(edge_left >= buf_cont[i][j][l]->GetPos()) edge_left = buf_cont[i][j][l]->GetPos();
                  if(edge_right <= buf_cont[i][j][l]->GetPos()) edge_right = buf_cont[i][j][l]->GetPos();
                  l=k;
                }

              /*
              if(((((buf_cont[i][j][k]->GetPos()-buf_cont[i][j][l]->GetPos())>0.) || fabs(buf_cont[i][j][k]->GetPos()-buf_cont[i][j][l]->GetPos())<1e-4))
                && fabs(edge_left - buf_cont[i][j][l]->GetPos())<=1.1)
                {
                  
                printf("%d %d %.8f\n" , i, j, fabs(edge_left - buf_cont[i][j][l]->GetPos()));

                  buf_cont[i][j][k]->Add(buf_cont[i][j][l]);
                  flag_used.insert(l);
                  if(edge_left >= buf_cont[i][j][l]->GetPos()) edge_left = buf_cont[i][j][l]->GetPos();
                  l=k;
                }
              else if(((((buf_cont[i][j][l]->GetPos()-buf_cont[i][j][k]->GetPos())>0.) || fabs(buf_cont[i][j][k]->GetPos()-buf_cont[i][j][l]->GetPos())<1e-4))
                && fabs(buf_cont[i][j][l]->GetPos()-edge_right)<=1.1)
                {
                  buf_cont[i][j][k]->Add(buf_cont[i][j][l]);
                  flag_used.insert(l);
                  if(edge_right <= buf_cont[i][j][l]->GetPos()) edge_right = buf_cont[i][j][l]->GetPos();
                                        printf("%d %d %.8f\n" , i, j, fabs(edge_right - buf_cont[i][j][l]->GetPos()));

                  l=k;

                }
                */
            } // for l
        } // for k

      for(int k=(int)buf_cont[i][j].size()-1; k>=0; --k){
        if(flag_used.size()>0 && flag_used.find(k)!=flag_used.end()){
          //std::cout << "k : " << k << std::endl;
          delete buf_cont[i][j][k];
          buf_cont[i][j].erase(buf_cont[i][j].begin() + k);
        }
      }
      //std::cout << "size after : " << buf_cont[i][j].size() << std::endl;

    } // for j
  } // for i


  return buf_cont;

}


std::vector< std::vector< FiberHitXUV* > > FiberAnalyzer::FindHit(std::vector< std::vector< std::vector< FiberHitAna* > > > &cont, ParaManager *par){

  // ang > 0 is U (only here)
  int id_x[7] = {0, 0, 0, 0, 0, 2, 0};
  int id_u[7] = {1, 1, 1, 2, 1, 0, 1};
  int id_v[7] = {2, 2, 2, 1, 2, 1, 2};

  std::vector< std::vector< FiberHitXUV* > > buf_cont;
  for(int i=0; i<7; ++i){
    std::vector<FiberHitXUV*> buf_v;
    buf_cont.emplace_back(buf_v);
  }

  for(int det=0; det<7; ++det){
    int id = 0;
    std::set<int> flag_x;
    std::set<int> flag_u;
    std::set<int> flag_v;
    int nh0 = cont[det][0].size();
    int nh1 = cont[det][1].size();
    int nh2 = cont[det][2].size();
    if((nh0==0)||(nh1==0)||(nh2==0)) continue;
    int ix=-1, iu=-1, iv=-1, nx=-1, nu=-1, nv=-1;
    ix = id_x[det];
    iu = id_u[det];
    iv = id_v[det];
    nx = cont[det][ix].size();
    nu = cont[det][iu].size();
    nv = cont[det][iv].size();
    double ang_uv = cont[det][iu][0]->GetAng();
    double cut_d = -999.;
    switch(det){
      case 0: cut_d = par->fiber_uft1_cut_d; break;
      case 1: cut_d = par->fiber_uft2_cut_d; break;
      case 2: cut_d = par->fiber_uft3_cut_d; break;
      case 3: cut_d = par->fiber_mft_cut_d ; break;
      case 4: cut_d = par->fiber_mft_cut_d ; break;
      case 5: cut_d = par->fiber_dft1_cut_d; break;
      case 6: cut_d = par->fiber_dft2_cut_d; break;
      default: break;
    }

    for(int i=0; i<nx; ++i){
      //if(flag_x.size()>0 && flag_x.find(i)!=flag_x.end()) continue;
      double pos_x = cont[det][ix][i]->GetPos();
      double hit_xu_x = pos_x;
      double hit_xu_y = 999;
      double hit_xv_y = -999;
      int used_j = -999;
      int used_k = -999;

      for(int j=0; j<nu; ++j){
        //if(flag_u.size()>0 && flag_u.find(j)!=flag_u.end()) continue;
        double pos_u = cont[det][iu][j]->GetPos();
        pos_u = pos_u/cos(ang_uv);
        double temp_hit_xu_y = 999;
        if(ang_uv>0) temp_hit_xu_y = -tan((TMath::Pi()/2. - ang_uv))*pos_x + tan((TMath::Pi()/2. -ang_uv))*pos_u;
        else         temp_hit_xu_y =  tan((TMath::Pi()/2. + ang_uv))*pos_x - tan((TMath::Pi()/2. +ang_uv))*pos_u;

        for(int k=0; k<nv; ++k){
          //if(flag_v.size()>0 && flag_v.find(k)!=flag_v.end()) continue;
          double pos_v = cont[det][iv][k]->GetPos();
          pos_v = pos_v/cos(ang_uv);
          double temp_hit_xv_y = -999;
          if(ang_uv>0) temp_hit_xv_y =  tan((TMath::Pi()/2. - ang_uv))*pos_x - tan((TMath::Pi()/2. - ang_uv))*pos_v;
          else         temp_hit_xv_y = -tan((TMath::Pi()/2. + ang_uv))*pos_x + tan((TMath::Pi()/2. + ang_uv))*pos_v;
          //if(fabs(hit_xu_y-hit_xv_y) > fabs(temp_hit_xu_y-temp_hit_xv_y)){
          //  hit_xu_y = temp_hit_xu_y;
          //  hit_xv_y = temp_hit_xv_y;
          //  used_j = j;
          //  used_k = k;
          //}
          hit_xu_y = temp_hit_xu_y;
          hit_xv_y = temp_hit_xv_y;
          used_j = j;
          used_k = k;
          double buf_x = hit_xu_x + (hit_xu_y-hit_xv_y)/2.*tan(ang_uv/2.);
          double buf_y = (hit_xu_y+hit_xv_y)/2.;
          double buf_d = hit_xu_y-hit_xv_y;
          if((fabs(hit_xu_y-999)<1.)||(fabs(hit_xv_y+999)<1.)) continue;
          if(fabs(hit_xu_y-hit_xv_y) > cut_d) continue;
          FiberHitXUV *buf_hit = new FiberHitXUV(buf_x, buf_y, buf_d, cont[det][ix][i], cont[det][iu][used_j], cont[det][iv][used_k], id++);
          buf_cont[det].emplace_back(buf_hit);

        } // for k
      } // for j
      //if((fabs(hit_xu_y-999)<1.)||(fabs(hit_xv_y+999)<1.)) continue;
      //if(fabs(hit_xu_y-hit_xv_y) > cut_d) continue;
      //double buf_x = hit_xu_x + (hit_xu_y-hit_xv_y)/2.*tan(ang_uv/2.);
      //double buf_y = (hit_xu_y+hit_xv_y)/2.;
      //double buf_d = hit_xu_y-hit_xv_y;
      //FiberHitXUV *buf_hit = new FiberHitXUV(buf_x, buf_y, buf_d, cont[det][ix][i], cont[det][iu][used_j], cont[det][iv][used_k], id++);
      //buf_cont[det].emplace_back(buf_hit);

      flag_x.insert(i);
      flag_u.insert(used_j);
      flag_v.insert(used_k);

    } // for i

    switch(det){
      case 0:
      case 1:
        buf_cont[det] = FiberAnalyzer::DeleteDupXUV(buf_cont[det]);
        break;

      case 2:
        if(par->flag_dup_xuv_uft3) buf_cont[det] = FiberAnalyzer::DeleteDupXUV(buf_cont[det]);
        break;

      case 3:
      case 4:
        if(par->flag_dup_xuv_mft12) buf_cont[det] = FiberAnalyzer::DeleteDupXUV(buf_cont[det]);
        break;

      case 5:
      case 6:
        buf_cont[det] = FiberAnalyzer::DeleteDupXUV(buf_cont[det]);
        break;

      default: break;
    }

  } // for det

  return buf_cont;

}


std::vector< FiberTrackAna* > FiberAnalyzer::DeleteDup( std::vector< FiberTrackAna* > &cont){

  std::vector< FiberTrackAna* > buf_cont;
  std::sort(cont.begin(), cont.end(), [](auto const& a, auto const& b) { return a->GetChi2() < b->GetChi2(); });
  int num_det = cont[0]->GetContXUV().size();
  //std::cout << "num_det : " << num_det << std::endl;
  std::vector< std::set<int> > used_hit;
  for(int i=0; i<num_det; ++i){
    std::set<int> buf_set;
    used_hit.emplace_back(buf_set);
  }
  int num_cont = cont.size();
  for(int i=0; i<num_cont; ++i){
    bool flag = true;
    std::vector<FiberHitXUV*> buf_xuv = cont[i]->GetContXUV();
    for(int j=0; j<(int)buf_xuv.size(); ++j){
      int id = buf_xuv[j]->GetID();
      if(used_hit[j].size()>0 && used_hit[j].find(id)!=used_hit[j].end()) flag = false;
    }
    if(flag){
      buf_cont.emplace_back(cont[i]);
      for(int j=0; j<(int)buf_xuv.size(); ++j){
        int id = buf_xuv[j]->GetID();
        used_hit[j].insert(id);
      }
    }
    else delete cont[i];
  }

  return buf_cont;

}

std::vector< FiberTrackAna* > FiberAnalyzer::DeleteDupCombi( std::vector< FiberTrackAna* > &cont){

  std::vector< FiberTrackAna* > buf_cont;
  if(cont.size()==0) return buf_cont;

  int det_first = cont.front()->GetContHit().at(0)->GetDet();
  int id_offset = -999;
  if(     det_first==0 || det_first==1) id_offset = 0; // UFT12
  else if(det_first==3 || det_first==4) id_offset = 3; // MFT12
  else if(det_first==5 || det_first==6) id_offset = 5; // DFT12

  std::sort(cont.begin(), cont.end(), [](auto const& a, auto const& b) {
      return a->GetNlayer() != b->GetNlayer() ? a->GetNlayer() > b->GetNlayer() : a->GetChi2() < b->GetChi2(); });
  int num_lay = 6;
  //std::cout << "num_det : " << num_det << std::endl;
  std::vector< std::set<int> > used_hit;
  for(int i=0; i<num_lay; ++i){
    std::set<int> buf_set;
    used_hit.emplace_back(buf_set);
  }
  int num_cont = cont.size();
  for(int i=0; i<num_cont; ++i){
    bool flag = true;
    std::vector<FiberHitAna*> buf_hit = cont[i]->GetContHit();
    for(int j=0; j<(int)buf_hit.size(); ++j){
      int fib = buf_hit[j]->GetFib();
      int det = buf_hit[j]->GetDet();
      int lay = buf_hit[j]->GetLay();
      int buf_id = (det - id_offset)*3 + lay;
      if(used_hit[buf_id].size()>0 && used_hit[buf_id].find(fib)!=used_hit[buf_id].end()) flag = false;
    }
    if(flag){
      buf_cont.emplace_back(cont[i]);
      for(int j=0; j<(int)buf_hit.size(); ++j){
        int fib = buf_hit[j]->GetFib();
        int det = buf_hit[j]->GetDet();
        int lay = buf_hit[j]->GetLay();
        int buf_id = (det-id_offset)*3 + lay;
        used_hit[buf_id].insert(fib);
      }
    }
    else delete cont[i];
  }

  return buf_cont;

}


std::vector< FiberTrackAna* > FiberAnalyzer::DeleteSame( std::vector< FiberTrackAna* > &cont){

  std::vector< FiberTrackAna* > buf_cont;
  std::sort(cont.begin(), cont.end(), [](auto const& a, auto const& b) { return a->IsFlagXUV() > b->IsFlagXUV(); });
  for(int i=0; i<(int)cont.size(); ++i){
    bool flag = true;
    bool flag_same = false;
    for(int j=(int)cont.size()-1; j>i; --j){
      if( cont[i]->GetNlayer()==6 && cont[j]->GetNlayer()==6 && cont[i]->IsFlagXUV()!=cont[j]->IsFlagXUV() ){
        flag_same = true;
        for(int k=0; k<6; ++k){
          if(cont[i]->GetHit(k)->GetFib() != cont[j]->GetHit(k)->GetFib()) flag_same = false;
        }
        if(flag_same){
          flag = cont[i]->IsFlagXUV();
          if(flag){
            delete cont[j];
            cont.erase(cont.begin() + j);
          }
          else break;
        }
      }
    }
    if(!flag_same) buf_cont.emplace_back(cont[i]);
    else if(flag){
      buf_cont.emplace_back(cont[i]);
    }
    else delete cont[i];
  }

  for(int i=0; i<(int)buf_cont.size(); ++i){
    bool flag_same = false;
    for(int j=(int)buf_cont.size()-1; j>i; --j){
      if( (buf_cont[i]->GetNlayer()==4 && buf_cont[j]->GetNlayer()==4) || (buf_cont[i]->GetNlayer()==5 && buf_cont[j]->GetNlayer()==5) ){
        flag_same = true;
        for(int k=0; k<buf_cont[i]->GetNlayer(); ++k){
          if(buf_cont[i]->GetHit(k)->GetDid() != buf_cont[j]->GetHit(k)->GetDid()) flag_same = false;
          if(buf_cont[i]->GetHit(k)->GetFib() != buf_cont[j]->GetHit(k)->GetFib()) flag_same = false;
        }
        if(flag_same){
          delete buf_cont[j];
          buf_cont.erase(buf_cont.begin() + j);
        }
      }
    }
  }

  return buf_cont;

}


std::vector< FiberTrackAna* > FiberAnalyzer::DeleteInclusive( std::vector< FiberTrackAna* > &cont){

  std::vector< FiberTrackAna* > buf_cont;
  std::sort(cont.begin(), cont.end(), [](auto const& a, auto const& b) { return a->GetNlayer() > b->GetNlayer(); });

  for(int i=0; i<(int)cont.size(); ++i){
    bool flag = false;
    for(int j=(int)cont.size()-1; j>i; --j){
      if( cont[i]->GetNlayer() >  cont[j]->GetNlayer() ){
        if(cont[i]->IsInclusive(cont[j])){
          delete cont[j];
          cont.erase(cont.begin() + j);
        }

      }
      else if( cont[i]->GetNlayer() <  cont[j]->GetNlayer() ){
        if(cont[j]->IsInclusive(cont[i])){
          flag = true;
          break;
        }
      }
    }
    if(flag){
      delete cont[i];
      cont.erase(cont.begin() + i);
      i--;
    }
    else buf_cont.emplace_back(cont[i]);

  }

  return buf_cont;

}


std::vector< TrackHit* > FiberAnalyzer::DeleteInclusiveTrackHit( std::vector< TrackHit* > &cont){

  std::vector< TrackHit* > buf_cont;
  std::sort(cont.begin(), cont.end(), [](auto const& a, auto const& b) { return a->GetNum() > b->GetNum(); });

  for(int i=0; i<(int)cont.size(); ++i){
    bool flag = false;
    for(int j=(int)cont.size()-1; j>i; --j){
      if( cont[i]->GetNum() >  cont[j]->GetNum() ){
        if(cont[i]->IsInclusive(cont[j])){
          delete cont[j];
          cont.erase(cont.begin() + j);
        }

      }
      else if( cont[i]->GetNum() <  cont[j]->GetNum() ){
        if(cont[j]->IsInclusive(cont[i])){
          flag = true;
          break;
        }
      }
    }
    if(flag){
      delete cont[i];
      cont.erase(cont.begin() + i);
      i--;
    }
    else buf_cont.emplace_back(cont[i]);

  }

  return buf_cont;

}

std::vector< TrackHit* > FiberAnalyzer::DeleteDupTrackHit( std::vector< TrackHit* > &cont){

  std::vector< TrackHit* >  buf_cont;
  if(cont.size()==0) return buf_cont;
  int num_cont = cont.size();

  //std::sort(cont.begin(), cont.end(), [](auto const& a, auto const& b) {
  //    return a->GetNumMDC() != b->GetNumMDC() ? a->GetNumMDC() > b->GetNumMDC() :
  //    (a->GetNumFiber() != b->GetNumFiber() ? a->GetNumFiber() > b->GetNumFiber() :
  //     (a->GetTrack()->GetChi2() == b->GetTrack()->GetChi2() ? a->GetTrack()->IsFlagXUV() > b->GetTrack()->IsFlagXUV() :
  //     a->GetTrack()->GetChi2() < b->GetTrack()->GetChi2() ) )
  //    ;});
  std::sort(cont.begin(), cont.end(), [](auto const& a, auto const& b) {
      return a->GetNumFiber() != b->GetNumFiber() ? a->GetNumFiber() > b->GetNumFiber() :
      (a->GetNumMDC() != b->GetNumMDC() ? a->GetNumMDC() > b->GetNumMDC() :
       a->GetChi2NDF() < b->GetChi2NDF()) ;});

  std::map< int, std::set<int> > used_hit_fib;
  std::map< int, std::set<int> > used_hit_mdc;

  for(int i=0; i<num_cont; ++i){
    TrackHit *track_hit =cont[i];
    std::vector<int> _fiberhit = track_hit->GetFiberHit();
    std::vector<int> _mdchit_bestdif = track_hit->GetMDCHit();
    bool flag = true;

    for(size_t i=0;i<6;++i)
    {
      if(_fiberhit[i]==-1) continue;
      int fib = _fiberhit[i];
      int did = i;
      if(used_hit_fib[did].size()>0 && used_hit_fib[did].find(fib)!=used_hit_fib[did].end()) flag = false;
    }
    for(size_t i=0;i<17;++i)
    {
      if(_mdchit_bestdif[i]==-1) continue;
      int wir = _mdchit_bestdif[i];
      int did = i;
      if(used_hit_mdc[did].size()>0 && used_hit_mdc[did].find(wir)!=used_hit_mdc[did].end()) flag = false;
    }

    if(flag){
      buf_cont.emplace_back(cont[i]);
      for(size_t i=0;i<6;++i)
      {
        if(_fiberhit[i]==-1) continue;
        int fib = _fiberhit[i];
        int did = i;
        used_hit_fib[did].insert(fib);
      }
      for(size_t i=0;i<17;++i)
      {
        if(_mdchit_bestdif[i]==-1) continue;
        int wir = _mdchit_bestdif[i];
        int did = i;
        used_hit_mdc[did].insert(wir);
      }

    }
    else delete cont[i];
  }

  return buf_cont;

}


std::vector< TrackHit* > FiberAnalyzer::DeleteDupTrackHitMDC( std::vector< TrackHit* > &cont){

  std::vector< TrackHit* >  buf_cont;
  if(cont.size()==0) return buf_cont;

  int num_cont = cont.size();

  std::sort(cont.begin(), cont.end(), [](auto const& a, auto const& b) {
      return a->GetNumMDC() != b->GetNumMDC() ? a->GetNumMDC() > b->GetNumMDC() :
      (a->GetNumFiber() != b->GetNumFiber() ? a->GetNumFiber() > b->GetNumFiber() :
       a->GetChi2NDF() < b->GetChi2NDF() )
      ;});

  std::map< int, std::set<int> > used_hit_mdc;

  for(int i=0; i<num_cont; ++i){
    TrackHit *track_hit =cont[i];
    std::vector<int> _mdchit_bestdif = track_hit->GetMDCHit();
    bool flag = true;


    for(size_t i=0;i<17;++i)
    {
      if(_mdchit_bestdif[i]==-1) continue;
      int wir = _mdchit_bestdif[i];
      int did = i;
      if(used_hit_mdc[did].size()>0 && used_hit_mdc[did].find(wir)!=used_hit_mdc[did].end()) flag = false;
    }

    if(flag){
      buf_cont.emplace_back(cont[i]);
      for(size_t i=0;i<17;++i)
      {
        if(_mdchit_bestdif[i]==-1) continue;
        int wir = _mdchit_bestdif[i];
        int did = i;
        used_hit_mdc[did].insert(wir);
      }
    }
    else delete cont[i];
  }

  return buf_cont;

}


std::vector< FiberHitXUV* > FiberAnalyzer::DeleteDupXUV( std::vector< FiberHitXUV* > &cont){

  std::vector< FiberHitXUV* > buf_cont;
  if(cont.size()==0) return buf_cont;

  std::sort(cont.begin(), cont.end(), [](auto const& a, auto const& b) {
      return a->GetD() <  b->GetD();});

  std::map< int, std::set<int> > used_hit;

  int num_cont = cont.size();

  for(int i=0; i<num_cont; ++i){
    FiberHitXUV *hit = cont[i];
    bool flag = true;

    for(int j=0; j<3; ++j){
      int fib = hit->GetHit(j)->GetFib();
      int did = hit->GetHit(j)->GetDid();
      if(used_hit[did].size()>0 && used_hit[did].find(fib)!=used_hit[did].end()) flag = false;
    }

    if(flag){
      buf_cont.emplace_back(cont[i]);
      for(int j=0; j<3; ++j){
        int fib = hit->GetHit(j)->GetFib();
        int did = hit->GetHit(j)->GetDid();
        used_hit[did].insert(fib);
      }

    }
    else delete cont[i];
  }

  return buf_cont;

}


TVector3 FiberAnalyzer::GetVertexPoint( const TVector3 & Xin, const TVector3 & Xout,
    const TVector3 & Pin, const TVector3 & Pout,
    double & dist )
{
  double xi=Xin.x(), yi=Xin.y(), xo=Xout.x(), yo=Xout.y();
  double ui=Pin.x()/Pin.z(), vi=Pin.y()/Pin.z();
  double uo=Pout.x()/Pout.z(), vo=Pout.y()/Pout.z();

  double z=((xi-xo)*(uo-ui)+(yi-yo)*(vo-vi))/
    ((uo-ui)*(uo-ui)+(vo-vi)*(vo-vi));
  double x1=xi+ui*z, y1=yi+vi*z;
  double x2=xo+uo*z, y2=yo+vo*z;
  dist=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

  return TVector3( 0.5*(x1+x2), 0.5*(y1+y2), z);
}


double FiberAnalyzer::CalcPhiDif(const double phi1, const double phi2){

  double diff = phi1 - phi2;
  while( diff >  TMath::Pi() ){
    diff -= 2*TMath::Pi();
  }
  while( diff < -TMath::Pi() ){
    diff += 2*TMath::Pi();
  }

  return diff;
}

double FiberAnalyzer::GetPSB_R(int _seg)
{
  double _r = -999.;

  if(_seg<23){
    if (0 == (_seg % 2)) { // inner PSB
      _r   = 217.;
    } else { // outer PSB
      _r   = 227.75;
    }
  }
  else{
    if (0 == ((_seg-23) % 2)) { // inner PSB
      _r   = 217.;
    } else { // outer PSB
      _r   = 227.75;
    }
  }
  return _r;
}


double FiberAnalyzer::GetPSB_Phi(int _seg)
{
  double _phi = -999.;

  if(_seg<23){
    if (0 == (_seg % 2)) { // inner PSB
      _phi = TMath::Pi() * (9.15 + 14.7 * ((double)(_seg / 2))) / 180.0;
    } else { // outer PSB
      _phi = TMath::Pi() * (16.5 + 14.7 * ((double)((_seg - 1) / 2))) / 180.0;
    }
  }
  else{
    if (0 == ((_seg-23) % 2)) { // inner PSB
      _phi = TMath::Pi() * (189.15 + 14.7 * ((float)((_seg-23) / 2))) / 180.0;
    } else { // outer PSB
      _phi = TMath::Pi() * (196.5 + 14.7 * ((float)(((_seg-23) - 1) / 2))) / 180.0;
    }
  }

  _phi += TMath::Pi()/2.;
  if( _phi > TMath::Pi() ) _phi -= 2*TMath::Pi();

  return _phi;
}


std::map< std::string, std::vector<FiberTrackAna*> > FiberAnalyzer::FiberTracking( std::vector< std::vector< std::vector< FiberHitAna* > > >& FiberHitClCont,
                                                                      ParaManager *par, std::vector<std::unique_ptr<genfit::AbsMeasurement> >& ListHits_PSB)
{
  std::map< std::string, std::vector<FiberTrackAna*> > FiberTrackCont;
  std::vector< std::vector<FiberHitXUV*> > FiberXUVCont = FindHit(FiberHitClCont, par);


  // UFT12
  int nt_uft12 = 0;

  if(FiberXUVCont[0].size()>0 && FiberXUVCont[1].size()>0)
    {
      std::vector<FiberTrackAna*> buf_track;
      for(int i=0; i<(int)FiberXUVCont[0].size(); ++i)
        {
          for(int j=0; j<(int)FiberXUVCont[1].size(); ++j)
            {
              std::vector<FiberHitXUV*>   buf_xuv;
              buf_xuv.emplace_back(FiberXUVCont[0][i]);
              buf_xuv.emplace_back(FiberXUVCont[1][j]);
              FiberTrackAna *track = new FiberTrackAna(buf_xuv, par);

              if(!par->flag_uft12_combi)                      buf_track.emplace_back(track);
              else if(track->GetChi2() < par->cut_chi2_uft12) buf_track.emplace_back(track);
            }
        }
      if((int)buf_track.size()>0) FiberTrackCont["uft12"] = DeleteDup(buf_track);
      nt_uft12     = FiberTrackCont["uft12"].size();
    }

  int buf_i_uft12 = -1;
  double buf_diff_uft12 = -9999.;
  for(int i=0; i<nt_uft12; ++i)
    {
      FiberTrackAna *track = FiberTrackCont["uft12"][i];
      double time_mean = track->GetTime();
      double dist_uft12 = fabs(time_mean);
      if(buf_diff_uft12<0 || dist_uft12 < buf_diff_uft12)
        {
          buf_i_uft12 = i;
          buf_diff_uft12 = dist_uft12;
        }
    }

  if(buf_i_uft12>-1) FiberTrackCont["uft12"][buf_i_uft12]->SetBest();

  if(par->flag_debug) std::cout << "- uft12 end" << std::endl;


  // UFT3MFT12
  if(FiberXUVCont[2].size()>0 && FiberXUVCont[3].size()>0 && FiberXUVCont[4].size()>0
                                  && par->flag_dup_xuv_mft12 && par->flag_dup_xuv_uft3)
    {
      std::vector<FiberTrackAna*> buf_track;
      for(int h=0; h<(int)FiberXUVCont[2].size(); ++h)
        for(int i=0; i<(int)FiberXUVCont[3].size(); ++i)
          for(int j=0; j<(int)FiberXUVCont[4].size(); ++j)
            {
              std::vector<FiberHitXUV*>   buf_xuv;
              buf_xuv.emplace_back(FiberXUVCont[2][h]);
              buf_xuv.emplace_back(FiberXUVCont[3][i]);
              buf_xuv.emplace_back(FiberXUVCont[4][j]);
              FiberTrackAna *track = new FiberTrackAna(buf_xuv, par);
              //track->CorrectMFT(par);
              buf_track.emplace_back(track);
            }

      FiberTrackCont["uft3mft12"] = DeleteDup(buf_track);
    }


  // MFT12
  //int nt_mft12 = 0;
  //int nt_mft12_xuv = 0;
  if(FiberXUVCont[3].size()>0 && FiberXUVCont[4].size()>0 && !par->flag_mft12_allcombi){
    std::vector<FiberTrackAna*> buf_track;
    for(int i=0; i<(int)FiberXUVCont[3].size(); ++i){
      for(int j=0; j<(int)FiberXUVCont[4].size(); ++j){
        std::vector<FiberHitXUV*>   buf_xuv;
        buf_xuv.emplace_back(FiberXUVCont[3][i]);
        buf_xuv.emplace_back(FiberXUVCont[4][j]);
        FiberTrackAna *track = new FiberTrackAna(buf_xuv, par);
        if(par->flag_mft12_posang){
          track->CorrectMFT(par);
          double buf_x = track->GetXmft();
          double buf_y = track->GetYmft();
          double buf_a = track->GetA();
          double buf_b = track->GetB();
          if( fabs(buf_x * 0.003 - buf_a)>0.3 || fabs(buf_y * 0.003 - buf_b)>0.3 ){ delete track; continue;}
        }
        if(!par->flag_mft12_combi)                      buf_track.emplace_back(track);
        else if(track->GetChi2() < par->cut_chi2_mft12) buf_track.emplace_back(track);
      }
    }

    if(par->flag_dup_mft12_xuv && (int)buf_track.size()>0) buf_track = DeleteDup(buf_track);
    FiberTrackCont["mft12"] = buf_track;

    //nt_mft12     = FiberTrackCont["mft12"].size();
    //nt_mft12_xuv = FiberTrackCont["mft12"].size();
    for(auto v: FiberTrackCont["mft12"]){
      for(int i=0; i<6; ++i){
        if(par->flag_dup_mft12_xuv) v->GetContHit().at(i)->SetUsed();
      }
      v->CorrectMFT(par);
      v->SetPosL();
    }
  }

  int num_combi_mft12 = 1;
  for(int i=3; i<5; ++i){
    for(int j=0; j<3; ++j){
      num_combi_mft12 *= ( (int)FiberHitClCont[i][j].size() + 1 );
    }
  }

  if(par->flag_debug) std::cout << "- num_combi_mft12 : " << num_combi_mft12 << std::endl;

  if(par->flag_mft12_combi && num_combi_mft12<par->cut_mft12_combi){

    std::vector<FiberTrackAna*> buf_track;
    for(int a=-1; a<(int)FiberHitClCont[3][0].size(); ++a){
      for(int b=-1; b<(int)FiberHitClCont[3][1].size(); ++b){
        for(int c=-1; c<(int)FiberHitClCont[3][2].size(); ++c){
          for(int d=-1; d<(int)FiberHitClCont[4][0].size(); ++d){
            for(int e=-1; e<(int)FiberHitClCont[4][1].size(); ++e){
              for(int f=-1; f<(int)FiberHitClCont[4][2].size(); ++f){
                std::vector<FiberHitAna*> buf_hit;
                int count = 0;
                if(a>-1 && !FiberHitClCont[3][0][a]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[3][0][a]); count++;}
                if(b>-1 && !FiberHitClCont[3][1][b]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[3][1][b]); count++;}
                if(c>-1 && !FiberHitClCont[3][2][c]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[3][2][c]); count++;}
                if(d>-1 && !FiberHitClCont[4][0][d]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[4][0][d]); count++;}
                if(e>-1 && !FiberHitClCont[4][1][e]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[4][1][e]); count++;}
                if(f>-1 && !FiberHitClCont[4][2][f]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[4][2][f]); count++;}
                if(count<4) continue;
                FiberTrackAna *track = new FiberTrackAna(buf_hit, par);
                track->SetFlagCombi();
                if(par->flag_mft12_posang){
                  track->CorrectMFTCombi(par);
                  double buf_x = track->GetXmft();
                  double buf_y = track->GetYmft();
                  double buf_a = track->GetA();
                  double buf_b = track->GetB();
                  if( fabs(buf_x * 0.003 - buf_a)>0.3 || fabs(buf_y * 0.003 - buf_b)>0.3 ){ delete track; continue;}
                }
                switch(track->GetNlayer()){
                  case 4: buf_track.emplace_back(track); break;
                  case 5: if(track->GetChi2() < par->cut_chi2_mft12) buf_track.emplace_back(track); else delete track; break;
                  case 6: if(track->GetChi2() < par->cut_chi2_mft12) buf_track.emplace_back(track); else delete track; break;
                  default: break;
                }
              }
            }
          }
        }
      }
    }

    for(int i=0; i<(int)buf_track.size(); ++i){
      for(int j=0; j<(int)ListHits_PSB.size(); ++j){
        double a_fiber = buf_track[i]->GetA();
        double b_fiber = buf_track[i]->GetB();
        double x_fiber = buf_track[i]->GetX();
        double y_fiber = buf_track[i]->GetY();
        double phi_psb   = GetPSB_Phi(ListHits_PSB[j]->getHitId()) + par->psb_rot_z*Deg2Rad;
        double r_psb     = GetPSB_R(ListHits_PSB[j]->getHitId());
        double z_psb     = ListHits_PSB[j]->getRawHitCoords()[1]*10.; //in mm

        double par_a = pow(a_fiber, 2) + pow(b_fiber, 2);
        double par_b = a_fiber * (x_fiber - par->psb_pos_x)+ b_fiber * (y_fiber - par->psb_pos_y);
        double par_c = pow(x_fiber - par->psb_pos_x, 2) + pow(y_fiber - par->psb_pos_y, 2) - pow(r_psb, 2);
        double z_fiber  = (-par_b + sqrt( pow(par_b,2) - par_a * par_c)) / par_a;

        double fiber_x_buf = x_fiber + a_fiber * z_fiber - par->psb_pos_x;
        double fiber_y_buf = y_fiber + b_fiber * z_fiber - par->psb_pos_y;
        double phi_fiber = atan2(fiber_y_buf, fiber_x_buf);

        if( fabs(CalcPhiDif(phi_psb, phi_fiber)) < par->cut_psb_phi ){
          if(fabs( (z_fiber - par->psb_pos_z) - z_psb)<par->cut_psb_z){
            if(buf_track[i]->IsFlagPSB()){
                if( fabs( CalcPhiDif(phi_psb, phi_fiber) ) < fabs( buf_track[i]->GetPSBDifPhi() ) ){
                  buf_track[i]->SetSegPSB(ListHits_PSB[j]->getHitId());
                  buf_track[i]->SetPSBDifZ(z_psb - (z_fiber - par->psb_pos_z));
                  buf_track[i]->SetPSBDifPhi(CalcPhiDif(phi_psb, phi_fiber));
                }
            }
            else{
              buf_track[i]->SetFlagPSB();
              buf_track[i]->SetSegPSB(ListHits_PSB[j]->getHitId());
              buf_track[i]->SetPSBDifZ(z_psb - (z_fiber - par->psb_pos_z));
              buf_track[i]->SetPSBDifPhi(CalcPhiDif(phi_psb, phi_fiber));
            }
          }
        }

      }
    }

    int num_buf = buf_track.size();
    for(int i=num_buf-1; i>=0; --i){
      if(!buf_track[i]->IsFlagPSB()){
        delete buf_track[i];
        buf_track.erase(buf_track.begin() + i);
      }
    }

    if(par->flag_dup_mft12_combi){
      std::vector<FiberTrackAna*> buf_track_tmp = DeleteDupCombi(buf_track);
      buf_track = buf_track_tmp;
    }

    for(auto v: buf_track){
      for(auto v2: v->GetContHit()){ v2->SetUsed(); }
      v->CorrectMFTCombi(par);
      v->SetPosL();
      v->DelFlagPSB();
    }

    for(auto v : buf_track){
      FiberTrackCont["mft12"].emplace_back(v);
    }

    FiberTrackCont["mft12"] = DeleteSame(FiberTrackCont["mft12"]);
    if(par->flag_debug) std::cout << "- before Inclusive : " << FiberTrackCont["mft12"].size() << std::endl;
    if(par->flag_mft12_inclusive) FiberTrackCont["mft12"] = DeleteInclusive(FiberTrackCont["mft12"]);
    if(par->flag_debug) std::cout << "- after  Inclusive : " << FiberTrackCont["mft12"].size() << std::endl;
    //nt_mft12 = FiberTrackCont["mft12"].size();

  }

  if(par->flag_mft12_pair){
    std::vector<std::pair<FiberHitAna*, FiberHitAna*> > pair_x;
    std::vector<std::pair<FiberHitAna*, FiberHitAna*> > pair_u;
    std::vector<std::pair<FiberHitAna*, FiberHitAna*> > pair_v;
    for(auto v1: FiberHitClCont[3][0]){
      for(auto v2: FiberHitClCont[4][0]){
        double pos1 = v1->GetPos();
        double pos2 = v2->GetPos();
        if( fabs(pos1 - pos2) < 20 ) pair_x.emplace_back(std::make_pair(v1, v2));
      }
    }
    for(auto v1: FiberHitClCont[3][1]){
      for(auto v2: FiberHitClCont[4][2]){
        double pos1 = v1->GetPos();
        double pos2 = v2->GetPos();
        if( fabs(pos1 - pos2) < 20 ) pair_u.emplace_back(std::make_pair(v1, v2));
      }
    }
    for(auto v1: FiberHitClCont[3][2]){
      for(auto v2: FiberHitClCont[4][1]){
        double pos1 = v1->GetPos();
        double pos2 = v2->GetPos();
        if( fabs(pos1 - pos2) < 20 ) pair_v.emplace_back(std::make_pair(v1, v2));
      }
    }

    int num_xp = pair_x.size();
    int num_up = pair_u.size();
    int num_vp = pair_v.size();
    int num_combi_pair = (num_xp+1) * (num_up+1) * (num_vp+1);
    if(par->flag_debug) std::cout << "- num_combi_pair_mft12 : " << num_combi_pair << std::endl;

    if(num_combi_pair < par->cut_mft12_combi){

      std::vector<FiberTrackAna*> buf_track;
      for(int a=-1; a<num_xp; ++a){
        for(int b=-1; b<num_up; ++b){
          for(int c=-1; c<num_vp; ++c){
            std::vector<FiberHitAna*> buf_hit;
            int count = 0;
            if( a>-1 && !pair_x[a].first->IsUsed() && !pair_x[a].second->IsUsed() ){
              buf_hit.emplace_back(pair_x[a].first); buf_hit.emplace_back(pair_x[a].second); count+=2; }
            if( b>-1 && !pair_u[b].first->IsUsed() && !pair_u[b].second->IsUsed() ){
              buf_hit.emplace_back(pair_u[b].first); buf_hit.emplace_back(pair_u[b].second); count+=2; }
            if( c>-1 && !pair_v[c].first->IsUsed() && !pair_v[c].second->IsUsed() ){
              buf_hit.emplace_back(pair_v[c].first); buf_hit.emplace_back(pair_v[c].second); count+=2; }
            if(count<4) continue;
            FiberTrackAna *track = new FiberTrackAna(buf_hit, par);
            track->SetFlagPair();
            if(par->flag_mft12_posang){
              track->CorrectMFTCombi(par);
              double buf_x = track->GetXmft();
              double buf_y = track->GetYmft();
              double buf_a = track->GetA();
              double buf_b = track->GetB();
              if( fabs(buf_x * 0.003 - buf_a)>0.3 || fabs(buf_y * 0.003 - buf_b)>0.3 ){ delete track; continue;}
            }
            switch(track->GetNlayer()){
              case 4: buf_track.emplace_back(track); break;
              case 5: if(track->GetChi2() < par->cut_chi2_mft12) buf_track.emplace_back(track); else delete track; break;
              case 6: if(track->GetChi2() < par->cut_chi2_mft12) buf_track.emplace_back(track); else delete track; break;
              default: break;
            }

          }
        }
      }

      for(int i=0; i<(int)buf_track.size(); ++i){
        for(int j=0; j<(int)ListHits_PSB.size(); ++j){
          double a_fiber = buf_track[i]->GetA();
          double b_fiber = buf_track[i]->GetB();
          double x_fiber = buf_track[i]->GetX();
          double y_fiber = buf_track[i]->GetY();
          double phi_psb   = GetPSB_Phi(ListHits_PSB[j]->getHitId()) + par->psb_rot_z*Deg2Rad;
          double r_psb     = GetPSB_R(ListHits_PSB[j]->getHitId());
          double z_psb     = ListHits_PSB[j]->getRawHitCoords()[1]*10.; //in mm

          double par_a = pow(a_fiber, 2) + pow(b_fiber, 2);
          double par_b = a_fiber * (x_fiber - par->psb_pos_x)+ b_fiber * (y_fiber - par->psb_pos_y);
          double par_c = pow(x_fiber - par->psb_pos_x, 2) + pow(y_fiber - par->psb_pos_y, 2) - pow(r_psb, 2);
          double z_fiber  = (-par_b + sqrt( pow(par_b,2) - par_a * par_c)) / par_a;

          double fiber_x_buf = x_fiber + a_fiber * z_fiber - par->psb_pos_x;
          double fiber_y_buf = y_fiber + b_fiber * z_fiber - par->psb_pos_y;
          double phi_fiber = atan2(fiber_y_buf, fiber_x_buf);
          if( fabs(CalcPhiDif(phi_psb, phi_fiber)) < par->cut_psb_phi ){
            if(fabs( (z_fiber - par->psb_pos_z) - z_psb)<par->cut_psb_z){
              if(buf_track[i]->IsFlagPSB()){
                if( fabs( CalcPhiDif(phi_psb, phi_fiber) ) < fabs( buf_track[i]->GetPSBDifPhi() ) ){
                  buf_track[i]->SetSegPSB(ListHits_PSB[j]->getHitId());
                  buf_track[i]->SetPSBDifZ(z_psb - (z_fiber - par->psb_pos_z));
                  buf_track[i]->SetPSBDifPhi(CalcPhiDif(phi_psb, phi_fiber));
                }
              }

              else{
                buf_track[i]->SetFlagPSB();
                buf_track[i]->SetSegPSB(ListHits_PSB[j]->getHitId());
                buf_track[i]->SetPSBDifZ(z_psb - (z_fiber - par->psb_pos_z));
                buf_track[i]->SetPSBDifPhi(CalcPhiDif(phi_psb, phi_fiber));
              }
            }
          }

        }
      }

      int num_buf = buf_track.size();
      for(int i=num_buf-1; i>=0; --i){
        if(!buf_track[i]->IsFlagPSB()){
          delete buf_track[i];
          buf_track.erase(buf_track.begin() + i);
        }
      }

      ////////////////////////////
      //  delete dup ???? ////////
      ////////////////////////////

      for(auto v: buf_track){
        for(auto v2: v->GetContHit()){ v2->SetUsed(); }
        v->CorrectMFTCombi(par);
        v->SetPosL();
        v->DelFlagPSB();
      }

      for(auto v : buf_track){
        FiberTrackCont["mft12"].emplace_back(v);
      }

      FiberTrackCont["mft12"] = DeleteSame(FiberTrackCont["mft12"]);
      if(par->flag_debug) std::cout << "- before Inclusive : " << FiberTrackCont["mft12"].size() << std::endl;
      if(par->flag_mft12_inclusive) FiberTrackCont["mft12"] = DeleteInclusive(FiberTrackCont["mft12"]);
      if(par->flag_debug) std::cout << "- after  Inclusive : " << FiberTrackCont["mft12"].size() << std::endl;
      //nt_mft12 = FiberTrackCont["mft12"].size();

    }

  }

  for(auto v: FiberTrackCont["mft12"]){
    v->SortContHit();
  }

  if(par->flag_debug) std::cout << "- mft12 end" << std::endl;


  // DFT12
  int nt_dft12 = 0;

  if(FiberXUVCont[5].size()>0 && FiberXUVCont[6].size()>0){
    std::vector<FiberTrackAna*> buf_track;
    for(int i=0; i<(int)FiberXUVCont[5].size(); ++i){
      for(int j=0; j<(int)FiberXUVCont[6].size(); ++j){
        std::vector<FiberHitXUV*>   buf_xuv;
        buf_xuv.emplace_back(FiberXUVCont[5][i]);
        buf_xuv.emplace_back(FiberXUVCont[6][j]);
        FiberTrackAna *track = new FiberTrackAna(buf_xuv, par);
        if(par->flag_dft12_cut){
          double tot_mean = track->GetTOT();
          double x_buf = track->GetXdet();
          double y_buf = track->GetYdet();
          double a_buf = track->GetA() * 1000;
          double b_buf = track->GetB() * 1000;
          bool flag_cut = false;
          if( pow(x_buf/60., 2.) + pow(y_buf/40., 2) > 1. ) flag_cut = true;
          if( fabs( x_buf * 23./50. - a_buf) > 20. ) flag_cut = true;
          if( fabs( y_buf * 35./60. - b_buf) > 30. ) flag_cut = true;
          if( tot_mean < par->cut_dft12_tot_max ) flag_cut = true;
          if(flag_cut){ delete track; continue; }
        }
        if(!par->flag_dft12_combi)                      buf_track.emplace_back(track);
        else if(track->GetChi2() < par->cut_chi2_dft12) buf_track.emplace_back(track);
      }
    }
    if((int)buf_track.size()>0) FiberTrackCont["dft12"] = DeleteDup(buf_track);
    nt_dft12     = FiberTrackCont["dft12"].size();
  }

  for(auto v: FiberTrackCont["dft12"]){
    for(int i=0; i<6; ++i){
      v->GetContHit().at(i)->SetUsed();
    }
  }

  int num_combi_dft12 = 1;
  for(int i=5; i<7; ++i){
    for(int j=0; j<3; ++j){
      num_combi_dft12 *= ( (int)FiberHitClCont[i][j].size() + 1 );
    }
  }

  if(par->flag_debug) std::cout << "- num_combi_dft12 : " << num_combi_dft12 << std::endl;

  if(par->flag_dft12_combi && num_combi_dft12<par->cut_dft12_combi){

    std::vector<FiberTrackAna*> buf_track;
    for(int a=-1; a<(int)FiberHitClCont[5][0].size(); ++a){
      for(int b=-1; b<(int)FiberHitClCont[5][1].size(); ++b){
        for(int c=-1; c<(int)FiberHitClCont[5][2].size(); ++c){
          for(int d=-1; d<(int)FiberHitClCont[6][0].size(); ++d){
            for(int e=-1; e<(int)FiberHitClCont[6][1].size(); ++e){
              for(int f=-1; f<(int)FiberHitClCont[6][2].size(); ++f){
                std::vector<FiberHitAna*> buf_hit;
                int count = 0;
                if(a>-1 && !FiberHitClCont[5][0][a]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[5][0][a]); count++;}
                if(b>-1 && !FiberHitClCont[5][1][b]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[5][1][b]); count++;}
                if(c>-1 && !FiberHitClCont[5][2][c]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[5][2][c]); count++;}
                if(d>-1 && !FiberHitClCont[6][0][d]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[6][0][d]); count++;}
                if(e>-1 && !FiberHitClCont[6][1][e]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[6][1][e]); count++;}
                if(f>-1 && !FiberHitClCont[6][2][f]->IsUsed() ) {buf_hit.emplace_back(FiberHitClCont[6][2][f]); count++;}
                if(count<4) continue;
                FiberTrackAna *track = new FiberTrackAna(buf_hit, par);

                if(par->flag_dft12_cut){
                  double tot_mean = track->GetTOT();
                  double x_buf = track->GetXdet();
                  double y_buf = track->GetYdet();
                  double a_buf = track->GetA() * 1000;
                  double b_buf = track->GetB() * 1000;
                  bool flag_cut = false;
                  if( pow(x_buf/60., 2.) + pow(y_buf/40., 2) > 1. ) flag_cut = true;
                  if( fabs( x_buf * 23./50. - a_buf) > 20. ) flag_cut = true;
                  if( fabs( y_buf * 35./60. - b_buf) > 30. ) flag_cut = true;
                  if( tot_mean < 75. ) flag_cut = true;
                  if(flag_cut){ delete track; continue; }
                }

                switch(track->GetNlayer()){
                  case 4: buf_track.emplace_back(track); break;
                  case 5: if(track->GetChi2() < par->cut_chi2_dft12) buf_track.emplace_back(track); else delete track; break;
                  case 6: if(track->GetChi2() < par->cut_chi2_dft12) buf_track.emplace_back(track); else delete track; break;
                  default: break;
                }
              }
            }
          }
        }
      }
    }
    if(par->flag_debug) std::cout << "- dft12 combi" << std::endl;

    std::vector<FiberTrackAna*> buf_track_tmp = DeleteDupCombi(buf_track);
    if(par->flag_debug) std::cout << "- dft12 dup" << std::endl;
    buf_track = buf_track_tmp;

    for(auto v: buf_track){
      for(auto v2: v->GetContHit()){ v2->SetUsed(); }
    }

    for(auto v : buf_track){
      FiberTrackCont["dft12"].emplace_back(v);
    }
    nt_dft12 = FiberTrackCont["dft12"].size();

  }

  int buf_i_dft12 = -1;
  double buf_diff_dft12 = -9999.;
  for(int i=0; i<nt_dft12; ++i){
    FiberTrackAna *track = FiberTrackCont["dft12"][i];
    double tot_mean  = track->GetTOT();
    double time_mean = track->GetTime();
    double dist_dft12 =
      pow( (tot_mean  - par->cut_dft12_tot_mean)  / par->cut_dft12_tot_sig , 2.) +
      pow( (time_mean - par->cut_dft12_time_mean) / par->cut_dft12_time_sig, 2.);
    if(buf_diff_dft12<0 || dist_dft12 < buf_diff_dft12){
      buf_i_dft12 = i;
      buf_diff_dft12 = dist_dft12;
    }
  }

  if(buf_i_dft12>-1) FiberTrackCont["dft12"][buf_i_dft12]->SetBest();

  if(par->flag_debug) std::cout << "- dft12 end" << std::endl;

  return FiberTrackCont;
}

