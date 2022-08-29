#include "FiberTrackAna.hh"
#include "ConstantParameter.hh"
#include "FiberHitXUV.hh"
#include "TString.h"
#include <iostream>
#include <math.h>

FiberTrackAna::FiberTrackAna(std::vector<FiberHitXUV*> &cont, ParaManager *par){
  Tracking(cont);
  _xtgt = _x + _a*par->fiber_tgt_pos_z;
  _ytgt = _y + _b*par->fiber_tgt_pos_z;
  for(int i=0; i<(int)cont.size(); ++i){
    _cont_xuv.emplace_back(cont[i]);
  }
  SetResidual();
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

void FiberTrackAna::Tracking(std::vector<FiberHitXUV*> &cont){

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

  TrackFitting(nlayer, &w[0], &z[0], &angle[0], &s[0]);

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
  chisqr /= n-4.;
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
    _cont_res_ch_pos[layer] = {(pos_hit0 - pos_track0), ch_hit0, pos_hit0};
    layer++;

    double ang1  = v->GetHit1()->GetAng();
    double posz1 = v->GetHit1()->GetZ();
    double pos_track1 = (_x + _a*posz1) * cos(ang1) + (_y + _b*posz1) * sin(ang1);
    double pos_hit1   = v->GetHit1()->GetPos();
    double ch_hit1    = v->GetHit1()->GetFib();
    _cont_res[layer] = pos_hit1 - pos_track1;
    _cont_res_ch_pos[layer] = {(pos_hit1 - pos_track1), ch_hit1, pos_hit1};
    layer++;

    double ang2  = v->GetHit2()->GetAng();
    double posz2 = v->GetHit2()->GetZ();
    double pos_track2 = (_x + _a*posz2) * cos(ang2) + (_y + _b*posz2) * sin(ang2);
    double pos_hit2   = v->GetHit2()->GetPos();
    double ch_hit2    = v->GetHit2()->GetFib();
    _cont_res[layer] = pos_hit2 - pos_track2;
    _cont_res_ch_pos[layer] = {(pos_hit2 - pos_track2), ch_hit2, pos_hit2};
    layer++;
  }

}
