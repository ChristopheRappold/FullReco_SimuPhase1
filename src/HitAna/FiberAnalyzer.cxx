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
      for(size_t k=0; k<buf_cont[i][j].size(); ++k){
        if(flag_used.size()>0 && flag_used.find(k)!=flag_used.end()) continue;
        double pos_cluster = buf_cont[i][j][k]->GetPos();
        double edge_left  = pos_cluster;
        double edge_right = pos_cluster;
        for(size_t l = k+1; l < buf_cont[i][j].size(); ++l){
          if(flag_used.size()>0 && flag_used.find(l)!=flag_used.end()) continue;
          if(((buf_cont[i][j][k]->GetPos()-buf_cont[i][j][l]->GetPos())>0.) && fabs(edge_left - buf_cont[i][j][l]->GetPos())<1.1){
            buf_cont[i][j][k]->Add(buf_cont[i][j][l]);
            flag_used.insert(l);
            edge_left = buf_cont[i][j][l]->GetPos();
            l=k;
          }
          else if(((buf_cont[i][j][l]->GetPos()-buf_cont[i][j][k]->GetPos())>0.) && fabs(buf_cont[i][j][l]->GetPos()-edge_right)<1.1){
            buf_cont[i][j][k]->Add(buf_cont[i][j][l]);
            flag_used.insert(l);
            edge_right = buf_cont[i][j][l]->GetPos();
            l=k;
          }
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
