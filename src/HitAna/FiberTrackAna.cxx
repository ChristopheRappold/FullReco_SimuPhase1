#include "FiberTrackAna.hh"
#include "ConstantParameter.hh"
#include "FiberHitXUV.hh"
#include "TString.h"
#include "TProfile2D.h"
#include <iostream>
#include <numeric>
#include <math.h>

FiberTrackAna::FiberTrackAna(std::vector<FiberHitXUV*> &cont, ParaManager *par){
  Tracking(cont, par);
  for(int i=0; i<(int)cont.size(); ++i){
    _cont_xuv.emplace_back(cont[i]);
    _cont_hit.emplace_back(cont[i]->GetHit(0));
    _cont_hit.emplace_back(cont[i]->GetHit(1));
    _cont_hit.emplace_back(cont[i]->GetHit(2));
  }
  std::sort(_cont_hit.begin(), _cont_hit.end(), [](auto const& a, auto const& b) { return a->GetZOrg() < b->GetZOrg(); });
  SetResidual();
  SetTOT();
  SetTime();
  for(auto hit: _cont_hit){ _cont_did_fib[hit->GetDid()] = hit->GetFib(); }
  _flag_xuv = true;
  _valid = true;
}

FiberTrackAna::FiberTrackAna(std::vector<FiberHitAna*> &cont, ParaManager *par){
  Tracking(cont,par);
  for(int i=0; i<(int)cont.size(); ++i){
    _cont_hit.emplace_back(cont[i]);
  }
  std::sort(_cont_hit.begin(), _cont_hit.end(), [](auto const& a, auto const& b) { return a->GetZOrg() < b->GetZOrg(); });
  //SetResidual();
  SetTOT();
  SetTime();
  for(auto hit: _cont_hit){ _cont_did_fib[hit->GetDid()] = hit->GetFib(); }
  _valid = true;
}

FiberTrackAna::FiberTrackAna(FiberHitXUV *hit_xuv, PSBHitAna *hit_psb, ParaManager *par){
  Tracking(hit_xuv, hit_psb, par);
  _cont_xuv.emplace_back(hit_xuv);
  _cont_hit.emplace_back(hit_xuv->GetHit(0));
  _cont_hit.emplace_back(hit_xuv->GetHit(1));
  _cont_hit.emplace_back(hit_xuv->GetHit(2));
  SetTOT();
  SetTime();
  for(auto hit: _cont_hit){ _cont_did_fib[hit->GetDid()] = hit->GetFib(); }
  _valid = true;
}

FiberTrackAna::FiberTrackAna(double x_tgt, double y_tgt, PSBHitAna *hit_psb, ParaManager *par){
  Tracking(x_tgt, y_tgt, hit_psb, par);
  _valid = true;
}

FiberTrackAna::~FiberTrackAna() {}


void FiberTrackAna::Print(void){
  std::cout << "-- track -----"   << std::endl;
  std::cout << "x    : " << _x    << std::endl;
  std::cout << "y    : " << _y    << std::endl;
  std::cout << "a    : " << _a    << std::endl;
  std::cout << "b    : " << _b    << std::endl;
  std::cout << "chi2 : " << _chi2 << std::endl;
  std::cout << "--------------"   << std::endl;
}

void FiberTrackAna::Tracking(std::vector<FiberHitXUV*> &cont, ParaManager *par){

  int num = cont.size();
  std::vector<double> w;
  std::vector<double> z;
  std::vector<double> angle;
  std::vector<double> s;
  for(int i=0; i<num; ++i){
    w.emplace_back(    cont[i]->GetHit0()->GetRes() );
    w.emplace_back(    cont[i]->GetHit1()->GetRes() );
    w.emplace_back(    cont[i]->GetHit2()->GetRes() );
    z.emplace_back(    cont[i]->GetHit0()->GetZ()   );
    z.emplace_back(    cont[i]->GetHit1()->GetZ()   );
    z.emplace_back(    cont[i]->GetHit2()->GetZ()   );
    angle.emplace_back(cont[i]->GetHit0()->GetAng() );
    angle.emplace_back(cont[i]->GetHit1()->GetAng() );
    angle.emplace_back(cont[i]->GetHit2()->GetAng() );
    s.emplace_back(    cont[i]->GetHit0()->GetPos() );
    s.emplace_back(    cont[i]->GetHit1()->GetPos() );
    s.emplace_back(    cont[i]->GetHit2()->GetPos() );
  }

  int nlayer = 3 * num;
  _nlayer =  nlayer;
  _z_mean = std::accumulate(z.begin(), z.end(), 0.0) / z.size();

  TrackFitting(nlayer, &w[0], &z[0], &angle[0], &s[0]);

  _xtgt = _x + _a * par->target_pos_z;
  _ytgt = _y + _b * par->target_pos_z;
  _xdet = _x + _a * _z_mean;
  _ydet = _y + _b * _z_mean;
  _xmft = _x + _a * (par->fiber_mft1_pos_z + par->fiber_mft2_pos_z)/2.;
  _ymft = _y + _b * (par->fiber_mft1_pos_z + par->fiber_mft2_pos_z)/2.;

}

void FiberTrackAna::Tracking(std::vector<FiberHitAna*> &cont, ParaManager *par){

  int num = cont.size();
  if(num<4){
    std::cout << "Num of Layer small : " << num << std::endl;
    return;
  }
  std::vector<double> w;
  std::vector<double> z;
  std::vector<double> angle;
  std::vector<double> s;
  for(int i=0; i<num; ++i){
    w.emplace_back(    cont[i]->GetRes() );
    z.emplace_back(    cont[i]->GetZ()   );
    angle.emplace_back(cont[i]->GetAng() );
    s.emplace_back(    cont[i]->GetPos() );
  }

  int nlayer = num;
  _nlayer =  nlayer;
  _z_mean = std::accumulate(z.begin(), z.end(), 0.0) / z.size();

  TrackFitting(nlayer, &w[0], &z[0], &angle[0], &s[0]);

  _xtgt = _x + _a * par->target_pos_z;
  _ytgt = _y + _b * par->target_pos_z;
  _xdet = _x + _a * _z_mean;
  _ydet = _y + _b * _z_mean;
  _xmft = _x + _a * (par->fiber_mft1_pos_z + par->fiber_mft2_pos_z)/2.;
  _ymft = _y + _b * (par->fiber_mft1_pos_z + par->fiber_mft2_pos_z)/2.;

}

void FiberTrackAna::Tracking(FiberHitXUV *hit_xuv, PSBHitAna *hit_psb, ParaManager *par){

  double x_uft3 = hit_xuv->GetPosX();
  double y_uft3 = hit_xuv->GetPosY();
  double z_uft3 = par->fiber_uft3_pos_z;

  double buf_r   = hit_psb->GetR();
  double buf_phi = hit_psb->GetPhi();
  double x_psb = buf_r * cos(buf_phi) + par->psb_pos_x;
  double y_psb = buf_r * sin(buf_phi) + par->psb_pos_y;
  double z_psb = hit_psb->GetZ()      + par->psb_pos_z;

  _a = (x_psb - x_uft3) / (z_psb - z_uft3);
  _x = x_uft3 - _a * z_uft3;
  _b = (y_psb - y_uft3) / (z_psb - z_uft3);
  _y = y_uft3 - _b * z_uft3;
  _chi2 = -9999.;

  _nlayer = 3;
  _z_mean = (z_uft3 + z_psb) / 2.;

  _xtgt = _x + _a * par->target_pos_z;
  _ytgt = _y + _b * par->target_pos_z;
  _xdet = _x + _a * _z_mean;
  _ydet = _y + _b * _z_mean;
  _xmft = _x + _a * (par->fiber_mft1_pos_z + par->fiber_mft2_pos_z)/2.;
  _ymft = _y + _b * (par->fiber_mft1_pos_z + par->fiber_mft2_pos_z)/2.;

}

void FiberTrackAna::Tracking(double x_tgt, double y_tgt, PSBHitAna *hit_psb, ParaManager *par){

  double z_tgt = par->target_pos_z;

  double buf_r   = hit_psb->GetR();
  double buf_phi = hit_psb->GetPhi();
  double x_psb = buf_r * cos(buf_phi) + par->psb_pos_x;
  double y_psb = buf_r * sin(buf_phi) + par->psb_pos_y;
  double z_psb = hit_psb->GetZ()      + par->psb_pos_z;

  _a = (x_psb - x_tgt) / (z_psb - z_tgt);
  _x = x_tgt - _a * z_tgt;
  _b = (y_psb - y_tgt) / (z_psb - z_tgt);
  _y = y_tgt - _b * z_tgt;
  _chi2 = -9999.;

  _nlayer = 0;
  _z_mean = (z_tgt + z_psb) / 2.;

  _xtgt = _x + _a * par->target_pos_z;
  _ytgt = _y + _b * par->target_pos_z;
  _xdet = _x + _a * _z_mean;
  _ydet = _y + _b * _z_mean;
  _xmft = _x + _a * (par->fiber_mft1_pos_z + par->fiber_mft2_pos_z)/2.;
  _ymft = _y + _b * (par->fiber_mft1_pos_z + par->fiber_mft2_pos_z)/2.;

}

bool FiberTrackAna::TrackFitting(int nlayer, double* w, double* z, double* angle, double* s){

  std::size_t n = nlayer;

  //if(n < TrLocalMinNHits ) return status_ = false;

  double ct[n],st[n];

  for( std::size_t i=0; i<n; ++i ){
    double ww = w[i];
    double aa = angle[i];
    w[i] = 1./(ww*ww);
    ct[i] = cos(aa);
    st[i] = sin(aa);
  }

  double matrx[16], *mtp[4], fitp[4];
  mtp[0]=&matrx[0]; mtp[1]=&matrx[4]; mtp[2]=&matrx[8]; mtp[3]=&matrx[12];

  for( int i=0; i<4; ++i ){
    fitp[i]=0.0;
    for( int j=0; j<4; ++j ){
      mtp[i][j]=0.0;
    }
  }

  for( std::size_t i=0; i<n; ++i ){
    double ww=w[i], zz=z[i], ss=s[i], ctt=ct[i], stt=st[i];
    mtp[0][0] += ww*ctt*ctt;
    mtp[0][1] += ww*zz*ctt*ctt;
    mtp[0][2] += ww*ctt*stt;
    mtp[0][3] += ww*zz*ctt*stt;
    mtp[1][1] += ww*zz*zz*ctt*ctt;
    mtp[1][2] += ww*zz*ctt*stt;
    mtp[1][3] += ww*zz*zz*ctt*stt;
    mtp[2][2] += ww*stt*stt;
    mtp[2][3] += ww*zz*stt*stt;
    mtp[3][3] += ww*zz*zz*stt*stt;

    fitp[0] += ww*ss*ctt;
    fitp[1] += ww*zz*ss*ctt;
    fitp[2] += ww*ss*stt;
    fitp[3] += ww*zz*ss*stt;
  }
  mtp[1][0]=mtp[0][1]; mtp[2][0]=mtp[0][2]; mtp[3][0]=mtp[0][3];
  mtp[2][1]=mtp[1][2]; mtp[3][1]=mtp[1][3]; mtp[3][2]=mtp[2][3];

  std::vector<int> indxc(n), indxd(n), ipiv(n);

  if( GaussJordan( mtp, 4, fitp, &indxc[0],
        &indxd[0], &ipiv[0] )==false ){
    std::cerr <<  "Fitting fails" << std::endl;
    return false;
  }

  _x = fitp[0];
  _y = fitp[2];
  _a = fitp[1];
  _b = fitp[3];

  double chisqr = 0.;
  for( std::size_t i=0; i<n; ++i ){
    double ww=w[i], zz=z[i];
    double scal=(_x+_a*zz)*ct[i] + (_y+_b*zz)*st[i];
    chisqr += ww*(s[i]-scal)*(s[i]-scal);
  }
  if(n>4) chisqr /= n-4.;
  else    chisqr = 9999.;
  _chi2 = chisqr;

  return true;
}


bool FiberTrackAna::GaussJordan( double **a, int n, double *b, int *indxr, int *indxc, int *ipiv ){

  for( int j=0; j<n; ++j ) ipiv[j]=0;
  for( int i=0; i<n; ++i ){
    double big=0.0;
    int irow=-1, icol=-1;
    for( int j=0; j<n; ++j )
      if( ipiv[j]!=1 )
        for( int k=0; k<n; ++k ){
          if( ipiv[k]==0 ){
            if( fabs(a[j][k])>=big ){
              big=fabs(a[j][k]);
              irow=j; icol=k;
            }
          }
          else if( ipiv[k]>1 ){
            std::cerr << "Singular Matrix" << std::endl;
            return false;
          }
        }
    ++(ipiv[icol]);

    if( irow!=icol ){
      for( int k=0; k<n; ++k ){
        double ta=a[irow][k];
        a[irow][k]=a[icol][k];
        a[icol][k]=ta;
      }
      double tb=b[irow];
      b[irow]=b[icol];
      b[icol]=tb;
    }

    indxr[i]=irow; indxc[i]=icol;

    if(a[icol][icol]==0.0 || fabs(a[icol][icol])<1e-30 || std::isinf(a[icol][icol]) ){
      std::cerr << "Singular Matrix"  << std::endl;
      return false;
    }
    double pivinv=1./a[icol][icol];
    a[icol][icol]=1.;
    for(int k=0; k<n; ++k) a[icol][k]*=pivinv;
    b[icol]*=pivinv;
    for( int k=0; k<n; ++k ){
      if(k!=icol){
        double d=a[k][icol];
        a[k][icol]=0.;
        for( int l=0; l<n; ++l ) a[k][l] -= a[icol][l]*d;
        b[k] -= b[icol]*d;
      }
    }
  }

  for(int l=n-1; l>=0; --l){
    if( indxr[l]!=indxc[l] ){
      for(int k=0; k<n; ++k ){
        double t=a[k][indxr[l]];
        a[k][indxr[l]]=a[k][indxc[l]];
        a[k][indxc[l]]=t;
      }
    }
  }
  return true;
}

void FiberTrackAna::SetResidual(){

  int layer = 0;
  for(auto v : _cont_xuv){
    double ang0  = v->GetHit0()->GetAng();
    double posz0 = v->GetHit0()->GetZ();
    double pos_track0 = (_x + _a*posz0) * cos(ang0) + (_y + _b*posz0) * sin(ang0);
    double pos_hit0   = v->GetHit0()->GetPos();
    double ch_hit0    = v->GetHit0()->GetFib();
    _cont_res[layer] = pos_hit0 - pos_track0;
    _cont_pos_res[layer][ch_hit0] = pos_hit0 - pos_track0;
    _cont_res_ch_pos[layer] = {(pos_hit0 - pos_track0), ch_hit0, pos_hit0};
    layer++;

    double ang1  = v->GetHit1()->GetAng();
    double posz1 = v->GetHit1()->GetZ();
    double pos_track1 = (_x + _a*posz1) * cos(ang1) + (_y + _b*posz1) * sin(ang1);
    double pos_hit1   = v->GetHit1()->GetPos();
    double ch_hit1    = v->GetHit1()->GetFib();
    _cont_res[layer] = pos_hit1 - pos_track1;
    _cont_pos_res[layer][ch_hit1] = pos_hit1 - pos_track1;
    _cont_res_ch_pos[layer] = {(pos_hit1 - pos_track1), ch_hit1, pos_hit1};
    layer++;

    double ang2  = v->GetHit2()->GetAng();
    double posz2 = v->GetHit2()->GetZ();
    double pos_track2 = (_x + _a*posz2) * cos(ang2) + (_y + _b*posz2) * sin(ang2);
    double pos_hit2   = v->GetHit2()->GetPos();
    double ch_hit2    = v->GetHit2()->GetFib();
    _cont_res[layer] = pos_hit2 - pos_track2;
    _cont_pos_res[layer][ch_hit2] = pos_hit2 - pos_track2;
    _cont_res_ch_pos[layer] = {(pos_hit2 - pos_track2), ch_hit2, pos_hit2};
    layer++;
  }

}

void FiberTrackAna::SetTOT(){

  double tot = 0;
  int    num = 0;
  for(auto hit: _cont_hit){
    tot += hit->GetTOT();
    num++;
  }
  tot /= num;

  _tot = tot;

}

void FiberTrackAna::SetTime(){

  double time = 0;
  int    num = 0;
  for(auto hit: _cont_hit){
    time += hit->GetTime();
    num++;
  }
  time /= num;

  _time = time;

}

void FiberTrackAna::CorrectMFT(ParaManager *par){

  for(auto v: _cont_xuv){
    int det = v->GetHit0()->GetDet();
    if(det==3 || det==4){
      for(int i=0; i<3; ++i){
        int lay = v->GetHit(i)->GetLay();
        int fib = v->GetHit(i)->GetFib();
        int seg = -1;
        if(fib<128) seg = 0;
        else        seg = 1;
        double r_ring = 210/2.; //mm
        //double cor_min = par->fiber_mft_cor[det-3][lay][seg][0];
        //double cor_max = par->fiber_mft_cor[det-3][lay][seg][1];
        double buf_z      = v->GetHit(i)->GetZ();
        double buf_z_org  = v->GetHit(i)->GetZOrg();
        double buf_ang = v->GetHit(i)->GetAng();
        double buf_x = _x + _a * buf_z;
        double buf_y = _y + _b * buf_z;
        double hit_s = v->GetHit(i)->GetPosOrg();
        //double hit_s =      buf_x * cos(buf_ang) + buf_y * sin(buf_ang);
        double hit_l = -1 * buf_x * sin(buf_ang) + buf_y * cos(buf_ang);
        double buf_cent_s =      par->fiber_mft1_pos_x * cos(buf_ang) + par->fiber_mft1_pos_y * sin(buf_ang);
        double buf_cent_l = -1 * par->fiber_mft1_pos_x * sin(buf_ang) + par->fiber_mft1_pos_y * cos(buf_ang);
        if(det==4){
          buf_cent_s =      par->fiber_mft2_pos_x * cos(buf_ang) + par->fiber_mft2_pos_y * sin(buf_ang);
          buf_cent_l = -1 * par->fiber_mft2_pos_x * sin(buf_ang) + par->fiber_mft2_pos_y * cos(buf_ang);
        }
        hit_s -= buf_cent_s;
        hit_l -= buf_cent_l;
        // double l_buf =  sqrt( pow(r_ring, 2.) - pow(hit_s, 2.) );
        //double cor_buf =  cor_max - (fabs(hit_s) - 12.6)/(83.0 - 12.6) * (cor_max - cor_min);
        //double r = (pow(l_buf,2) + pow(cor_buf,2)) / (cor_buf * 2);
        //double theta = asin(l_buf/r);
        //double theta_cor = asin(hit_l/r);
        //double cor = r*cos(theta_cor) - r*cos(theta);
        //double theta_cor_pos = asin( (cor_max - cor_min)/70.4 );
        //double cor_pos = cor * tan(theta_cor_pos);
        //if(seg==0) cor_pos *= -1;
        double p[3];
        for(int j=0; j<3; ++j){
          p[j] = par->fiber_mft_cor_par[det-3][lay][seg][j];
        }
        double cor = pow(hit_s, 2) / p[0] + pow(hit_l, 2)/p[1] + p[2];

        //std::cout << "\ni : " << i << std::endl;
        //std::cout << "buf_z : " << buf_z << std::endl;
        //std::cout << "hit_s : " << hit_s  << std::endl;
        //std::cout << "hit_l : " << hit_l  << std::endl;
        //std::cout << "hit : " << v->GetHit(i)->GetPos() - buf_cent_s << std::endl;
        //std::cout << "cor : " << cor << std::endl;
        //std::cout << "cor_pos : " << cor_pos << std::endl;
        v->GetHit(i)->SetZ(buf_z_org + cor);
        //v->GetHit(i)->SetPos( v->GetHit(i)->GetPosOrg() + cor_pos );
      }
    }
    else continue;
  }

  Tracking(_cont_xuv, par);
  SetResidual();

}

void FiberTrackAna::CorrectMFTCombi(ParaManager *par){

  for(auto v: _cont_hit){
    int det = v->GetDet();
    if(det==3 || det==4){
      int lay = v->GetLay();
      int fib = v->GetFib();
      int seg = -1;
      if(fib<128) seg = 0;
      else        seg = 1;
      double r_ring = 210/2.; //mm
      //double cor_min = par->fiber_mft_cor[det-3][lay][seg][0];
      //double cor_max = par->fiber_mft_cor[det-3][lay][seg][1];
      double buf_z      = v->GetZ();
      double buf_z_org  = v->GetZOrg();
      double buf_ang = v->GetAng();
      double buf_x = _x + _a * buf_z;
      double buf_y = _y + _b * buf_z;
      double hit_s = v->GetPosOrg();
      //double hit_s =      buf_x * cos(buf_ang) + buf_y * sin(buf_ang);
      double hit_l = -1 * buf_x * sin(buf_ang) + buf_y * cos(buf_ang);
      double buf_cent_s =      par->fiber_mft1_pos_x * cos(buf_ang) + par->fiber_mft1_pos_y * sin(buf_ang);
      double buf_cent_l = -1 * par->fiber_mft1_pos_x * sin(buf_ang) + par->fiber_mft1_pos_y * cos(buf_ang);
      if(det==4){
        buf_cent_s =      par->fiber_mft2_pos_x * cos(buf_ang) + par->fiber_mft2_pos_y * sin(buf_ang);
        buf_cent_l = -1 * par->fiber_mft2_pos_x * sin(buf_ang) + par->fiber_mft2_pos_y * cos(buf_ang);
      }
      hit_s -= buf_cent_s;
      hit_l -= buf_cent_l;
      double l_buf =  sqrt( pow(r_ring, 2.) - pow(hit_s, 2.) );

      double p[3];
      for(int j=0; j<3; ++j){
        p[j] = par->fiber_mft_cor_par[det-3][lay][seg][j];
      }
      double cor = pow(hit_s, 2) / p[0] + pow(hit_l, 2)/p[1] + p[2];
      if(fabs(hit_l) > l_buf) cor = p[2];

      v->SetZ(buf_z_org + cor);
    }
    else continue;
  }

  Tracking(_cont_hit, par);

}

void FiberTrackAna::SetPosL(){

  for(auto hit: _cont_hit){
    // int lay = hit->GetLay();
    // int fib = hit->GetFib();
    double buf_z   = hit->GetZ();
    double buf_ang = hit->GetAng();
    double buf_x = _x + _a * buf_z;
    double buf_y = _y + _b * buf_z;
    double pos_l = -1 * buf_x * sin(buf_ang) + buf_y * cos(buf_ang);
    hit->SetPosL(pos_l);
  }

}
/*
void FiberTrackAna::CorrectMFTXY(ParaManager *par){

  for(auto v: _cont_xuv){
    int det = v->GetHit0()->GetDet();
    if(det==3 || det==4){
      for(int i=0; i<3; ++i){
        FiberHitAna *hit = v->GetHit(i);
        int lay = hit->GetLay();
        int fib = hit->GetFib();
        //double pos_s = hit->GetPos();
        double pos_l = hit->GetPosL();

        TProfile2D* prof = (TProfile2D*)par->fiber_mft_cor_xy->Get(Form("h47_3[%d]_pxy",(det-3)*3 + lay));
        int xbin = prof->GetXaxis()->GetNbins();
        int ybin = prof->GetYaxis()->GetNbins();
        double xmax = prof->GetXaxis()->GetXmax();
        double ymax = prof->GetYaxis()->GetXmax();
        double xmin = prof->GetXaxis()->GetXmin();
        double ymin = prof->GetYaxis()->GetXmin();
        double xstep = (xmax - xmin)/xbin;
        double ystep = (ymax - ymin)/ybin;

        //int ix = (pos_s - xmin) /xstep + 1;
        int ix = (fib   - xmin) /xstep + 1;
        int iy = (pos_l - ymin) /ystep + 1;
        double val = prof->GetBinContent(ix, iy);

        double cor_pos = -1 * val;
        //std::cout << "cor_pos : " << cor_pos << std::endl;

        hit->SetPos( hit->GetPosOrg() + cor_pos );
      }
    }
    else continue;
  }

  Tracking(_cont_xuv, par);
  SetResidual();

}
*/
void FiberTrackAna::SortContHit(){

  std::sort(_cont_hit.begin(), _cont_hit.end(), [](auto const& a, auto const& b) { return a->GetZ() < b->GetZ(); });

}

bool FiberTrackAna::IsInclusive(FiberTrackAna* track){

  bool flag = false;

  if(_nlayer <= track->GetNlayer()) flag = false;
  else {
    bool flag_tmp = true;
    for(auto cont : track->GetDidFib()){
      int did = cont.first;
      int fib = cont.second;;
      if( _cont_did_fib.find(did) == _cont_did_fib.end() || _cont_did_fib[did]!=fib) flag_tmp = false;
    }
    if(flag_tmp) flag = true;
  }

  return flag;

}



