#include "MWDCTracking.hh"
// #include "TFRSParameter.h"
#include "TMinuit.h"
#include "TMatrixFSym.h"

using namespace std;

MWDCTracking::MWDCTracking(ParaManager *par){
  //
  mwdc_tracking_assumption_resoltion = 0.5; //initial value. can be changed
  max_hit_combination=10;// initial value. can be changed
  min_plane_enabled=6;   // initial value. can be changed
  min_uv_plane_enabled=3;   // initial value. can be changed
  max_pair_combination=10000; // initial value. can be changed
  tracking_status=-1;
  chi2_track=99999.9; x_track=-99999.9; y_track=-99999.9; a_track=-99999.9; b_track=-99999.9;
  chi2_temp=99999.9; x_temp=-99999.9; y_temp=-99999.9; a_temp=-99999.9; b_temp=-99999.9;
  for(int i_plane=0;i_plane<16;i_plane++){
    LR_fixed[i_plane]=0;
    numhit[i_plane]=0; // <<== this must be here
    plane_enabled_temp[i_plane]=0;
    plane_enabled_track[i_plane]=0;
    theta[i_plane] =-99999.9;
    zplane[i_plane]=-99999.9;
    residual_track[i_plane]=-99999.9;
    residual_temp[i_plane]=-99999.9;
    time_wMaxToT[i_plane] = -9999;
    tot_max[i_plane] = -99999.9;
    num_pair[i_plane/2]=0; //<<== this must be here
    pair_enabled_temp[i_plane/2]=0;
    for(int i_hit=0; i_hit<64; i_hit++){
      length[i_plane][i_hit]  = -99999.9;
      wirepos[i_plane][i_hit] = -99999.9;
    }
  }
}

MWDCTracking::~MWDCTracking(){}


//--------------------------------------
//     AMSC tracking fit functions
//
// par[i]  >>  i: 0:x0, 1:a0, 2:y0, 3:b0 /////
// temp[i_plane][m] >>  m: 0:theta, 1:wpos, m:zplane, 3:drift length,
//                         4:factor in chi2 (=1./(dof*sigma*sigma), or 0 if plane disabled)
//

static double temp[16][5]={};

double distance_wire_track(int i_plane, Float_t par0, Float_t par1, Float_t par2, Float_t par3){
  return (abs(par0*cos(temp[i_plane][0])+par2*sin(temp[i_plane][0]) -temp[i_plane][1] + (par1*cos(temp[i_plane][0])+par3*sin(temp[i_plane][0]))*temp[i_plane][2]))/sqrt(1.0+pow((par1*cos(temp[i_plane][0])+par3*sin(temp[i_plane][0])), 2)) ;
}

double drift_length(int i_plane){
  return temp[i_plane][3];
}

void CylindricalChi2Funct(int& nDim, double* gout, double& result, double par[], int flg){
  double ans;  ///// par >> 0:x0, 1:y0, 2:a0, 3:b0 /////
  ans = 0.0;
  for(int i_plane=0; i_plane<16; i_plane++){
    if(0.0!=temp[i_plane][4]){
      ans += temp[i_plane][4] * pow((distance_wire_track(i_plane,par[0],par[1],par[2],par[3]) - drift_length(i_plane)), 2) ;
    }
  }
  result = ans;
}

void MWDCTracking::SetAssumedResolution(float setval){
  mwdc_tracking_assumption_resoltion = setval;
}


//-------------------------------------------------
//   Functions for MWDCTracking::Tracking

void MWDCTracking::StoreHit(int i_plane, Float_t wirepos_hit, Float_t length_hit, Float_t tot_hit, Float_t time_hit,  Float_t theta_hit, Float_t zplane_hit){
  if(0<= i_plane && i_plane<16){
    length[i_plane][ (numhit[i_plane]) ]  = length_hit  ;
    wirepos[i_plane][ (numhit[i_plane]) ] = wirepos_hit ;
    tot[i_plane][ (numhit[i_plane]) ]  = tot_hit  ;
    time[i_plane][ (numhit[i_plane]) ]  = time_hit  ;
    theta[i_plane]  = theta_hit  ;
    zplane[i_plane] = zplane_hit ;
    numhit[i_plane]++;
  }else{
    cout << "Wrong i_plane value (" << i_plane <<") is set..." << endl;
  }
  return ;
}

void MWDCTracking::StoreHit(MWDCHitAna *a, ParaManager *par){
  if(0<=a->GetPlane() && a->GetPlane()<16){
    if(a->GetLT()>par->mwdc_lt_valid_min && a->GetLT()<par->mwdc_lt_valid_max){
      length[a->GetPlane()][(numhit[a->GetPlane()])] = a->GetDl();;
      wirepos[a->GetPlane()][(numhit[a->GetPlane()])] = (par->mwdc_plane_sign[a->GetPlane()])*(5.0)*((float)(a->GetWire())-(par->mwdc_center_id[a->GetPlane()]));
      tot[a->GetPlane()][(numhit[a->GetPlane()])] = a->GetTOT();
      time[a->GetPlane()][(numhit[a->GetPlane()])] = a->GetLT();
      theta[a->GetPlane()] = -par->mwdc_plane_angle[a->GetPlane()];
      zplane[a->GetPlane()] = par->mwdc_zpos[a->GetPlane()];
      numhit[a->GetPlane()]++;
    }
  }
  else{
    std::cout<< "Wrong i_plane value ("<< a->GetPlane() <<") is set..."<<std::endl;
  }
  return ;
}

void MWDCTracking::SetMaxHitCombination(int setvalue){
  max_hit_combination = setvalue;
  return ;
}

void MWDCTracking::SetMinPlaneEnebaled(int setvalue){
  min_plane_enabled = setvalue;
  return ;
}



//-------------------------------------------------
//   Functions for MWDCTracking::Tracking

int MWDCTracking::Tracking(void){
  // 1: done
  // 0: give up because of too many combination
  //-1: give up because # of planes with hits is small.

  //------ set plane_enabled[] ----
  for(int i_plane=0; i_plane<16; i_plane++){
    if(0<numhit[i_plane]){
      plane_enabled_temp[i_plane]=1; // presently, just require all planes with hits to participate
    }
  }
  // temp test
  //  plane_enabled_temp[2]=0;    plane_enabled_temp[3]=0;  plane_enabled_temp[4]=0; plane_enabled_temp[5]=0;
  //      plane_enabled_temp[10]=0;   plane_enabled_temp[11]=0;  plane_enabled_temp[12]=0; plane_enabled_temp[13]=0;

  //------ condition for giving up -----
  if(NumHitCombinationTemp()>max_hit_combination){
    //  cout << "too many combination... "<< endl;
    return 0;
  }

  if(NumEnabledPlaneTemp()<min_plane_enabled){
    //cout << "num of hit planes is not enough... "<< endl;
    //printf("%d, %d \n",NumEnabledPlaneTemp(),min_plane_enabled);
    return -1;
  }


  //------try all combination by inverse matrix method------
  chi2_track = 999999.9;
  SetHitCombinationTemp0();
  for(int i=0; i<NumHitCombinationTemp(); i++){
    SetLRCombinationTemp0();
    for(int j=0; j<NumLRCombinationTemp(); j++){
      int status_tracking_matrix;
      status_tracking_matrix = TrackingInverseMatrix();
      if(1!=status_tracking_matrix){
        cout << "enabled plane choice seems bad..."<< endl;
        exit(-1);
      }
      if(chi2_temp < chi2_track){ // overwrite result (*_track)
        x_track = x_temp;
        y_track = y_temp;
        a_track = a_temp;
        b_track = b_temp;
        chi2_track = chi2_temp;
        CopyHitCombination(i_hit_used_temp, i_hit_used_track);
        CopyLRCombination(LR_used_temp, LR_used_track);
        CopyPlaneEnabled(plane_enabled_temp, plane_enabled_track);
        for(int i_plane=0;i_plane<16; i_plane++){
          residual_track[i_plane] = residual_temp[i_plane];
        }
      }
      LRCombinationTempIncrement(); // _temp is incremented (in enabled planes)
    }//j-end
    HitCombinationTempIncrement();// _temp is incremented (in enabled planes)
  }//i-end


  /* // comment out for cosmic test

  // ------- fit again (chi2 calculation = cylindrically) -------
  // ------- use previous result for initial values.
  //
  x_temp = x_track;
  y_temp = y_track;
  a_temp = a_track;
  b_temp = b_track;
  CopyHitCombination(i_hit_used_track, i_hit_used_temp);
  CopyLRCombination(LR_used_track, LR_used_temp);
  CopyPlaneEnabled(plane_enabled_track, plane_enabled_temp);
  TrackingFitCylindrical();
  x_track = x_temp;
  y_track = y_temp;
  a_track = a_temp;
  b_track = b_temp;
  cout << "chi2_before=" << chi2_track <<"   / chi2_new="<< chi2_temp <<endl;
  chi2_track = chi2_temp;
  CopyHitCombination(i_hit_used_temp, i_hit_used_track);
  CopyLRCombination(LR_used_temp, LR_used_track);
  CopyPlaneEnabled(plane_enabled_temp, plane_enabled_track);
  for(int i_plane=0;i_plane<16; i_plane++){
  residual_track[i_plane] = residual_temp[i_plane];
  }
  */


  return 1;
}



//-------------------------------------------------
//   Functions for MWDCTracking::Tracking

int MWDCTracking::Tracking_PairLR(void){
  // 1: done
  // 0: give up because of too many combination
  //-1: give up because # of planes with hits is small.

  //------ set plane_enabled[] ----
  for(int i_plane=0; i_plane<16; i_plane++){
    if(0<numhit[i_plane]){
      plane_enabled_temp[i_plane]=1; // presently, just require all planes with hits to participate
    }
  }

  //--------------------------------------------
  for(int i_pair=0; i_pair<8; i_pair++){
    pair_enabled_temp[i_pair]=1; // firstly, enable all pairs
  }
  // pair_enabled_temp[1]=0;  plane_enabled_temp[2]=0;  plane_enabled_temp[3]=0;
  // pair_enabled_temp[5]=0;  plane_enabled_temp[10]=0; plane_enabled_temp[11]=0;
  pair_enabled_temp[0]=0; // plane_enabled_temp[0]=0;  plane_enabled_temp[1]=0;
  pair_enabled_temp[4]=0; // plane_enabled_temp[8]=0;  plane_enabled_temp[9]=0;

  SetPairTemp();
  //        fprintf(stderr,"NumPairCombination = %d \n",NumPairCombinationTemp());

  //------ condition for giving up -----
  if(NumEnabledPlaneTemp()<min_plane_enabled){
    //cout << "num of hit planes is not enough... "<< endl;
    //  printf("%d, %d \n",NumEnabledPlaneTemp(),min_plane_enabled);
    return -1;
  }
  if(NumPairCombinationTemp()>max_pair_combination){
    //  cout << "too many combination... "<< endl;
    return 0;
  }
  if(0==NumPairCombinationTemp()){
    return -1;
  }

  //------try all combination by inverse matrix method------
  chi2_track = 999999.9;
  SetPairCombinationTemp0();
  for(int i=0; i<NumPairCombinationTemp(); i++){
    SetHitLRfromPair();
    if(NumEnabledPlaneTemp() >= min_plane_enabled){
      int status_tracking_matrix;
      status_tracking_matrix = TrackingInverseMatrix();
      if(1!=status_tracking_matrix){
        cout << "enabled plane choice seems bad..."<< endl;
        exit(-1);
      }
      if(chi2_temp < chi2_track){ // overwrite result (*_track)
        x_track = x_temp;
        y_track = y_temp;
        a_track = a_temp;
        b_track = b_temp;
        chi2_track = chi2_temp;
        CopyHitCombination(i_hit_used_temp, i_hit_used_track);
        CopyLRCombination(LR_used_temp, LR_used_track);
        CopyPlaneEnabled(plane_enabled_temp, plane_enabled_track);
        for(int i_plane=0;i_plane<16; i_plane++){
          residual_track[i_plane] = residual_temp[i_plane];
        }
      }
    }
    PairCombinationTempIncrement();// _temp is incremented (in enabled planes)
  }//i-end
  return 1;
}


// tracking method with FitCylindrical method
int MWDCTracking::Tracking_FitCylindrical(void){
  // 1: done
  // 0: give up because of too many combination
  //-1: give up because # of planes with hits is small.

  //------ set plane_enabled[] ----

  for(int i_plane=0; i_plane<16; i_plane++){
    if(0<numhit[i_plane]){
      plane_enabled_temp[i_plane]=1; // presently, just require all planes with hits to participate
    }
  }
  //------ condition for giving up -----
  if(NumHitCombinationTemp()>max_hit_combination){
    //  cout << "too many combination... "<< endl;
    return 0;
  }
  if(NumEnabledPlaneTemp()<min_plane_enabled){
    //cout << "num of hit planes is not enough... "<< endl;
    return -1;
  }

  //------try all combination by FitCylindrical method------
  chi2_track = 999999.9;
  SetHitCombinationTemp0();
  for(int i=0; i<NumHitCombinationTemp(); i++){
    SetLRCombinationTemp0();
    //for(int j=0; j<NumLRCombinationTemp(); j++){ //}
    int status_tracking_matrix;
    int status_tracking_fit = TrackingFitCylindrical();
    if(1!=status_tracking_fit){
      cout << "enabled plane choice seems bad..."<< endl;
      exit(-1);
    }
    if(chi2_temp < chi2_track){ // overwrite result (*_track)
      x_track = x_temp;
      y_track = y_temp;
      a_track = a_temp;
      b_track = b_temp;
      chi2_track = chi2_temp;
      CopyHitCombination(i_hit_used_temp, i_hit_used_track);
      CopyPlaneEnabled(plane_enabled_temp, plane_enabled_track);
      for(int i_plane=0;i_plane<16; i_plane++){
        residual_track[i_plane] = residual_temp[i_plane];
      }
    }
    HitCombinationTempIncrement();// _temp is incremented (in enabled planes)
  }//i-end
  return 1;
}

int MWDCTracking::Tracking_LRfixed(void){
  //TMWDCParameter* mwdc = (TMWDCParameter*) GetParameter("MWDCPar");
  // 1: done
  // 0: give up because of too many combination
  //-1: give up because # of planes with hits is small.
  //-2: give up because # of uv planes with hits is small.

  //------ set plane_enabled[] ----
  int nuvhit=0;
  for(int i_plane=0; i_plane<16; i_plane++){
    if(0<numhit[i_plane]){
      plane_enabled_temp[i_plane]=1; // presently, just require all planes with hits to participate
      if(i_plane>3&&i_plane<8) nuvhit++;
      if(i_plane>11&&i_plane<16) nuvhit++;
    }
  }
  //  plane_enabled_temp[2]=0; plane_enabled_temp[3]=0; plane_enabled_temp[8]=0; plane_enabled_temp[9]=0;

  //------ condition for giving up -----
  if(NumHitCombinationTemp()>max_hit_combination){
    //  cout << "too many combination... "<< endl;
    return 0;
  }

  if(NumEnabledPlaneTemp()<min_plane_enabled){
    //  cout << "num of hit planes is not enough... "<< endl;
    return -1;
  }

  if(nuvhit<min_uv_plane_enabled){
    //  cout << "num of hit planes is not enough... "<< endl;
    return -2;
  }

  //------try all combination by inverse matrix method------
  chi2_track = 999999.9;
  SetHitCombinationTemp0();
  for(int i=0; i<NumHitCombinationTemp(); i++){
    SetLRCombinationTemp0();
#if 1
    bool WRONGCOMB=false;
    int plane_enabled_save[16];
    CopyPlaneEnabled(plane_enabled_temp, plane_enabled_save);
    for(int i_pair=0;i_pair<8;i_pair++){
      if(numhit[i_pair*2]>0&&numhit[i_pair*2+1]>0&&false==LR_fixed[i_pair*2]&&false==LR_fixed[i_pair*2+1]){
        double tmp1=wirepos[i_pair*2][i_hit_used_temp[i_pair*2]];
        double tmp2=wirepos[i_pair*2+1][i_hit_used_temp[i_pair*2+1]];
        //if(length[i_pair*2][i_hit_used_temp[i_pair*2]]<0.2||
        //length[i_pair*2+1][i_hit_used_temp[i_pair*2+1]]<0.2)
        //continue;
        // if(0||i_pair!=0&&i_pair!=1&&i_pair!=4&&i_pair!=5)
        //  continue;
        double addlength=length[i_pair*2][i_hit_used_temp[i_pair*2]];
        addlength+=length[i_pair*2+1][i_hit_used_temp[i_pair*2+1]];
        if( addlength>4.){
          plane_enabled_temp[i_pair*2]=0;
          plane_enabled_temp[i_pair*2+1]=0;
          // WRONGCOMB=true;
          // break;
        }
        if(fabs(tmp1-tmp2)<=4.){
          LR_fixed[i_pair*2]=true;
          LR_fixed[i_pair*2+1]=true;
          plane_enabled_temp[i_pair*2]=1;
          plane_enabled_temp[i_pair*2+1]=1;
          if(tmp1-tmp2<0){
            LR_used_temp[i_pair*2]=0;
            LR_used_temp[i_pair*2+1]=1;
          }else{
            LR_used_temp[i_pair*2]=1;
            LR_used_temp[i_pair*2+1]=0;
          }
        }else{
          plane_enabled_temp[i_pair*2]=0;
          plane_enabled_temp[i_pair*2+1]=0;
          //    WRONGCOMB=true;
          //    break;
        }
      }
    }
    if(NumHitCombinationTemp()>max_hit_combination) continue;
    nuvhit=0;
    for(int i_plane=0; i_plane<16; i_plane++){
      if(0<plane_enabled_temp[i_plane]){
        if(i_plane>3&&i_plane<8) nuvhit++;
        if(i_plane>11&&i_plane<16) nuvhit++;
      }
    }
    if(nuvhit<min_uv_plane_enabled) continue;
    //    if(WRONGCOMB) continue;
#endif

    for(int j=0; j<NumLRCombinationTemp(); j++){
      int status_tracking_matrix;
      status_tracking_matrix = TrackingInverseMatrix();
      if(1!=status_tracking_matrix){
        cout << "enabled plane choice seems bad..."<< endl;
        exit(-1);
      }
      if(chi2_temp < chi2_track){ // overwrite result (*_track)
        x_track = x_temp;
        y_track = y_temp;
        a_track = a_temp;
        b_track = b_temp;
        chi2_track = chi2_temp;
        CopyHitCombination(i_hit_used_temp, i_hit_used_track);
        CopyLRCombination(LR_used_temp, LR_used_track);
        CopyPlaneEnabled(plane_enabled_temp, plane_enabled_track);
        for(int i_plane=0;i_plane<16; i_plane++){
          residual_track[i_plane] = residual_temp[i_plane];
        }
      }
      LRCombinationTempIncrement(); // _temp is incremented (in enabled planes)
    }//j-end
    HitCombinationTempIncrement();// _temp is incremented (in enabled planes)
  }//i-end

  //if(chi2_track>10000) return -3;
  return 1;
}

int MWDCTracking::ResidualExCalc(void){
  double w_hit_temp;
  CopyHitCombination(i_hit_used_track, i_hit_used_temp);
  CopyLRCombination(LR_used_track, LR_used_temp);
  CopyPlaneEnabled(plane_enabled_track, plane_enabled_temp);
  for(int i_plane=0;i_plane<16;i_plane++){
    if(plane_enabled_temp[i_plane]==1){
      plane_enabled_temp[i_plane]=0;
      TrackingInverseMatrix();
      w_hit_temp = wirepos[i_plane][(i_hit_used_temp[i_plane])] + (1.0-2.0*((float)(LR_used_temp[i_plane])))*length[i_plane][(i_hit_used_temp[i_plane])];
      residual_track_ex[i_plane] = ((x_temp+a_temp*zplane[i_plane])*cos(theta[i_plane])+(y_temp+b_temp*zplane[i_plane])*sin(theta[i_plane])-w_hit_temp);
      //for(int j_plane=0;j_plane<16;j_plane++)cout<<" "<<plane_enabled_temp[j_plane];
      //cout<<endl;
      //cout<<i_plane<<" "<<residual_woSelfData[i_plane]<<endl;
      plane_enabled_temp[i_plane]=1;
    }
    else residual_track_ex[i_plane] = -99999.9;
  }

  return 0;
}

int MWDCTracking::ResidualExCalc2(void){
  double w_hit_temp;
  CopyHitCombination(i_hit_used_track, i_hit_used_temp);
  CopyLRCombination(LR_used_track, LR_used_temp);
  CopyPlaneEnabled(plane_enabled_track, plane_enabled_temp);
  for(int i_plane=0;i_plane<16;i_plane++){
    if(plane_enabled_temp[i_plane]==1){
      for(int j_plane=0; j_plane<16; j_plane++){
        if(1==CheckSameType(i_plane,j_plane)){
          plane_enabled_temp[j_plane]=0;
        }
      }
      TrackingInverseMatrix();
      w_hit_temp = wirepos[i_plane][(i_hit_used_temp[i_plane])] + (1.0-2.0*((float)(LR_used_temp[i_plane])))*length[i_plane][(i_hit_used_temp[i_plane])];
      residual_track_ex[i_plane] = ((x_temp+a_temp*zplane[i_plane])*cos(theta[i_plane])+(y_temp+b_temp*zplane[i_plane])*sin(theta[i_plane])-w_hit_temp);
      //for(int j_plane=0;j_plane<16;j_plane++)cout<<" "<<plane_enabled_temp[j_plane];
      //cout<<endl;
      //cout<<i_plane<<" "<<residual_woSelfData[i_plane]<<endl;
      CopyPlaneEnabled(plane_enabled_track, plane_enabled_temp);
    }
    else residual_track_ex[i_plane] = -99999.9;
  }

  return 0;
}

int MWDCTracking::ResidualExCalc3(void){
  double w_hit_track, w_deduced_by_other;
  CopyHitCombination(i_hit_used_track, i_hit_used_temp);
  CopyLRCombination(LR_used_track, LR_used_temp);
  CopyPlaneEnabled(plane_enabled_track, plane_enabled_temp);
  for(int i_plane=0;i_plane<16;i_plane++){
    //printf("a%d\n",NumHitCombinationTemp());
    if(plane_enabled_temp[i_plane]==1 && (2000 >= (unsigned int)(NumHitCombinationTemp()))){
      for(int j_plane=0; j_plane<16; j_plane++){
        if(1==CheckSameType(i_plane,j_plane)){
          plane_enabled_temp[j_plane]=0;
        }
      }
      //printf("b%d\n",NumHitCombinationTemp());
      //------try all combination by inverse matrix method------
      double chi2_excalc3 = 999999.9;
      double a_excalc3,b_excalc3,x_excalc3,y_excalc3;
      int i_hit_used_excalc3[16];
      int LR_used_excalc3[16];
      SetHitCombinationTemp0();
      for(int i=0; i<NumHitCombinationTemp(); i++){
        SetLRCombinationTemp0();
        for(int j=0; j<NumLRCombinationTemp(); j++){
          //  printf("j=%d\n",j);
          int status_tracking_matrix;
          status_tracking_matrix = TrackingInverseMatrix();
          if(1!=status_tracking_matrix){
            cout << "enabled plane choice seems bad..."<< endl;
            exit(-1);
          }
          if(chi2_temp < chi2_excalc3){ // overwrite result (*_track)
            x_excalc3 = x_temp;
            y_excalc3 = y_temp;
            a_excalc3 = a_temp;
            b_excalc3 = b_temp;
            chi2_excalc3 = chi2_temp;
            CopyHitCombination(i_hit_used_temp, i_hit_used_excalc3);
            CopyLRCombination(LR_used_temp, LR_used_excalc3);
          }
          LRCombinationTempIncrement(); // _temp is incremented (in enabled planes)
        }//j-end
        HitCombinationTempIncrement();// _temp is incremented (in enabled planes)
      }//i-end

      // wire coordinate (w) deduced by others
      w_deduced_by_other = ((x_excalc3+a_excalc3*zplane[i_plane])*cos(theta[i_plane])+(y_excalc3+b_excalc3*zplane[i_plane])*sin(theta[i_plane]));

      // wire coordinate (w) for ith-plane (i_hit, LR "_track" should be used !!)
      w_hit_track = wirepos[i_plane][(i_hit_used_track[i_plane])] + (1.0-2.0*((float)(LR_used_track[i_plane])))*length[i_plane][(i_hit_used_track[i_plane])];

      // residual exclusive (method3)
      residual_track_ex[i_plane] = (w_deduced_by_other - w_hit_track);

      // set back "plane_enabled" for next plane calculation
      CopyPlaneEnabled(plane_enabled_track, plane_enabled_temp);

    }else{
      residual_track_ex[i_plane] = -99999.9;
    }
  }
  return 0;
}




int MWDCTracking::CheckSameType(int i_plane, int j_plane){
  int type_id[16]={0,0,0,0,1,1,2,2,0,0,0,0,1,1,2,2};
  if(0<= i_plane && i_plane<16 && 0<= j_plane && j_plane<16 ){
    if(type_id[i_plane]==type_id[j_plane]){
      return 1; // same type
    }else{
      return 0; // different type
    }
  }else{
    return -1; // error
  }
  return -1; // error
}


int MWDCTracking::FindTime_wMaxToT(void){
  for(int i_plane=0; i_plane<16; i_plane++){
    for(int ihit =0; ihit < numhit[i_plane]; ihit++){
      if(tot[i_plane][ihit]>tot_max[i_plane]){
        tot_max[i_plane] = tot[i_plane][ihit];
        time_wMaxToT[i_plane] = time[i_plane][ihit];
      }
    }
  }
  return 1;
}
int MWDCTracking::TrackingInverseMatrix(){
  //  1 = done
  //  0 = no inverse matrix exist

  // calculation of matrix elements (note p.63)
  Float_t m00=0.0;
  Float_t m01=0.0;
  Float_t m02=0.0;
  Float_t m03=0.0;
  Float_t m11=0.0;
  Float_t m12=0.0;
  Float_t m13=0.0;
  Float_t m22=0.0;
  Float_t m23=0.0;
  Float_t m33=0.0;
  Float_t v[4] = {0.0, 0.0 ,0.0 ,0.0}; // vector on right-hand side
  Float_t w_hit[16];
  for(int i_plane=0; i_plane<16; i_plane++){
    if(1==plane_enabled_temp[i_plane]){
      w_hit[i_plane] = wirepos[i_plane][(i_hit_used_temp[i_plane])] + (1.0-2.0*((float)(LR_used_temp[i_plane])))*length[i_plane][(i_hit_used_temp[i_plane])];
      m00 += cos(theta[i_plane])*cos(theta[i_plane]);
      m01 += zplane[i_plane]*cos(theta[i_plane])*cos(theta[i_plane]);
      m02 += sin(theta[i_plane])*cos(theta[i_plane]);
      m03 += zplane[i_plane]*sin(theta[i_plane])*cos(theta[i_plane]);
      m11 += zplane[i_plane]*zplane[i_plane]*cos(theta[i_plane])*cos(theta[i_plane]);
      m12 += zplane[i_plane]*sin(theta[i_plane])*cos(theta[i_plane]);
      m13 += zplane[i_plane]*zplane[i_plane]*sin(theta[i_plane])*cos(theta[i_plane]);
      m22 += sin(theta[i_plane])*sin(theta[i_plane]);
      m23 += zplane[i_plane]*sin(theta[i_plane])*sin(theta[i_plane]);
      m33 += zplane[i_plane]*zplane[i_plane]*sin(theta[i_plane])*sin(theta[i_plane]);
      v[0] += w_hit[i_plane]*cos(theta[i_plane]);
      v[1] += w_hit[i_plane]*zplane[i_plane]*cos(theta[i_plane]);
      v[2] += w_hit[i_plane]*sin(theta[i_plane]);
      v[3] += w_hit[i_plane]*zplane[i_plane]*sin(theta[i_plane]);
    }
  }


  // define matrix, set values
  TMatrixFSym M1(4);
  Float_t *M1_element = M1.GetMatrixArray();
  *(M1_element+(4*0+0)) = m00;
  *(M1_element+(4*0+1)) = m01;
  *(M1_element+(4*0+2)) = m02;
  *(M1_element+(4*0+3)) = m03;
  *(M1_element+(4*1+0)) = m01;
  *(M1_element+(4*1+1)) = m11;
  *(M1_element+(4*1+2)) = m12;
  *(M1_element+(4*1+3)) = m13;
  *(M1_element+(4*2+0)) = m02;
  *(M1_element+(4*2+1)) = m12;
  *(M1_element+(4*2+2)) = m22;
  *(M1_element+(4*2+3)) = m23;
  *(M1_element+(4*3+0)) = m03;
  *(M1_element+(4*3+1)) = m13;
  *(M1_element+(4*3+2)) = m23;
  *(M1_element+(4*3+3)) = m33;

  // Invert and check determinant
  Double_t det;
  M1.InvertFast(&det);
  // no inverted matrix
  if(0==det){
    cout << "no inverse matrix ..." <<endl ;
    chi2_temp = 999999.9;
    return 0;
  }

  // calculate (x,a,y,b) and chi2
  x_temp = 0.0; a_temp = 0.0; y_temp = 0.0; b_temp = 0.0; chi2_temp = 0.0;
  for(int i=0; i<4; i++){
    x_temp += (*(M1_element+(4*0+i))) * v[i] ;
    a_temp += (*(M1_element+(4*1+i))) * v[i] ;
    y_temp += (*(M1_element+(4*2+i))) * v[i] ;
    b_temp += (*(M1_element+(4*3+i))) * v[i] ;
  }

  Float_t dof;
  dof = (Float_t)(NumEnabledPlaneTemp()-4); // 4 is num. of parameter
  for(int i_plane=0;i_plane<16;i_plane++){
    if(1==plane_enabled_temp[i_plane]){
      residual_temp[i_plane] = ((x_temp+a_temp*zplane[i_plane])*cos(theta[i_plane])+(y_temp+b_temp*zplane[i_plane])*sin(theta[i_plane])-w_hit[i_plane]);
      chi2_temp += pow((residual_temp[i_plane])/(mwdc_tracking_assumption_resoltion),2)/dof;
    }else{
      residual_temp[i_plane]=99999.9;
    }
  }
  return 1;
}


int MWDCTracking::TrackingFitCylindrical(void){
  //  1 = done
  //  0 = fit failed
  // -1 = hit information missing

  Float_t dof;
  dof = (Float_t)(NumEnabledPlaneTemp()-4); // 4 is num. of parameter

  //-----Set information in temp[][]-----
  for(int i_plane=0; i_plane<16; i_plane++){
    if(1==plane_enabled_temp[i_plane]){
      temp[i_plane][0] = theta[i_plane];
      temp[i_plane][1] = wirepos[i_plane][i_hit_used_temp[i_plane]];
      temp[i_plane][2] = zplane[i_plane];
      temp[i_plane][3] = length[i_plane][i_hit_used_temp[i_plane]];
      temp[i_plane][4] = 1.0/(dof*(mwdc_tracking_assumption_resoltion)*(mwdc_tracking_assumption_resoltion));
    }else{
      temp[i_plane][4] = 0.0;
    }
  }

  //----Minimization----
  TMinuit minimizer(4);
  minimizer.SetFCN(CylindricalChi2Funct);
  minimizer.SetPrintLevel(-1);
  minimizer.DefineParameter(0,"x0",x_temp,0.01,-300.,300.);
  minimizer.DefineParameter(1,"a0",a_temp,0.01,-2.,2.);
  minimizer.DefineParameter(2,"y0",y_temp,0.01,-300.,300.);
  minimizer.DefineParameter(3,"b0",b_temp,0.01,-2.,2.);
  minimizer.Migrad();
  //minimizer.Release(0);

  //------ result -------
  minimizer.GetParameter(0, x_temp, xerr_temp);
  minimizer.GetParameter(1, a_temp, aerr_temp);
  minimizer.GetParameter(2, y_temp, yerr_temp);
  minimizer.GetParameter(3, b_temp, berr_temp);

  //----- residual_temp[], chi2_temp-----
  chi2_temp=0.0;
  for(int i_plane=0; i_plane<16; i_plane++){
    if(1==plane_enabled_temp[i_plane]){
      residual_temp[i_plane] = 1.0 * (distance_wire_track(i_plane,x_temp,a_temp,y_temp,b_temp) - drift_length(i_plane));
      chi2_temp += temp[i_plane][4] * pow(residual_temp[i_plane], 2) ;
    }else{
      residual_temp[i_plane] = 99999.9;
    }
  }

  return 1;
}



//-------------------------------------------------
//   Functions to treat
//     plane_enable[], i_hit_used[], LR_used[]

int MWDCTracking::NumEnabledPlaneTemp(){
  int count=0;
  for(int i_plane=0; i_plane<16; i_plane++){
    if(1==plane_enabled_temp[i_plane]){
      count++;
    }
  }
  return count;
}

int MWDCTracking::NumLRCombinationTemp(){
  int count=0;
  for(int i_plane=0; i_plane<16; i_plane++){
    if(1==plane_enabled_temp[i_plane]&&0==LR_fixed[i_plane]){
      count++;
    }
  }
  return (0x0001 << count); // 2^n
}


int MWDCTracking::NumHitCombinationTemp(){
  int count=1;
  for(int i_plane=0; i_plane<16; i_plane++){
    if(1==plane_enabled_temp[i_plane]){
      count*=numhit[i_plane];
    }
  }
  return count;
}


int MWDCTracking::NumPairCombinationTemp(){
  int count=1;
  for(int i_pair=0; i_pair<8; i_pair++){
    if(1==pair_enabled_temp[i_pair]){
      count*=num_pair[i_pair];
    }
  }
  return count;
}

int MWDCTracking::HitCombinationTempIncrement(){ // go to next combination
  // 1 : done
  // 0 : can not increment any more...
  for(int i_plane=0; i_plane<16; i_plane++){
    if(1==plane_enabled_temp[i_plane]){
      if(i_hit_used_temp[i_plane]==(numhit[i_plane]-1)){
        i_hit_used_temp[i_plane]=0; // (kuriagari) continue
      }else{
        i_hit_used_temp[i_plane]+=1;
        return 1;
      }
    }
  }
  return 0;
}

int MWDCTracking::SetHitLRfromPair(){
  for(int i_pair=0; i_pair<8; i_pair++){
    int i_plane1, i_plane2;
    i_plane1 = 2*i_pair;
    i_plane2 = 2*i_pair+1;
    if(1==pair_enabled_temp[i_pair]){
      if(0==num_pair[i_pair]){
        fprintf(stderr,"SetHitLRfromPair() was executed, but there num_pair[i_pair=%d]=0\n",i_pair);
        fprintf(stderr,"This should not happen...\n");
        exit(-1);
      }
      int ihit1,ihit2,ilr1,ilr2;
      ihit1 = i_hit1_pair[i_pair][(i_pair_used_temp[i_pair])];
      ihit2 = i_hit2_pair[i_pair][(i_pair_used_temp[i_pair])];
      ilr1  = i_lr1_pair[i_pair][(i_pair_used_temp[i_pair])];
      ilr2  = i_lr2_pair[i_pair][(i_pair_used_temp[i_pair])];
      //plane1
      if(-1!=ihit1){
        i_hit_used_temp[i_plane1] = ihit1;
        LR_used_temp[i_plane1]  = ilr1;
        plane_enabled_temp[i_plane1] = 1;
      }else{
        i_hit_used_temp[i_plane1] = ihit1;
        LR_used_temp[i_plane1]  = ilr1;
        plane_enabled_temp[i_plane1] = 0;
      }
      //plane2
      if(-1!=ihit2){
        i_hit_used_temp[i_plane2] = ihit2;
        LR_used_temp[i_plane2]  = ilr2;
        plane_enabled_temp[i_plane2] = 1;
      }else{
        i_hit_used_temp[i_plane2] = ihit2;
        LR_used_temp[i_plane2]  = ilr2;
        plane_enabled_temp[i_plane2] = 0;
      }
    }else{ // case for i_pair disabled
      plane_enabled_temp[i_plane1] = 0;
      plane_enabled_temp[i_plane2] = 0;
    }
  }
  return 1;
}

int MWDCTracking::SetPairTemp(){
  for(int i_pair=0; i_pair<8; i_pair++){
    num_pair[i_pair]=0;
    if(1==pair_enabled_temp[i_pair]){
      int i_plane1, i_plane2;
      i_plane1 = 2*i_pair;
      i_plane2 = 2*i_pair+1;
      if(1==plane_enabled_temp[i_plane1] && 1==plane_enabled_temp[i_plane2]){
        for(int i_hit1=0; i_hit1<numhit[i_plane1]; i_hit1++){
          for(int i_lr1=0; i_lr1<2; i_lr1++){
            for(int i_hit2=0; i_hit2<numhit[i_plane2]; i_hit2++){
              for(int i_lr2=0; i_lr2<2; i_lr2++){
                // check condition
                bool condition_for_pair;
                double pos_plane1, pos_plane2;
                pos_plane1 = wirepos[i_plane1][i_hit1] + (1.0-2.0*i_lr1)*length[i_plane1][i_hit1];
                pos_plane2 = wirepos[i_plane2][i_hit2] + (1.0-2.0*i_lr2)*length[i_plane2][i_hit2];
                condition_for_pair = (2.0 > fabs(pos_plane1 - pos_plane2));
                if(condition_for_pair){
                  int i = num_pair[i_pair];
                  if(i>=32) continue;
                  i_hit1_pair[i_pair][i]= i_hit1;
                  i_hit2_pair[i_pair][i]= i_hit2;
                  i_lr1_pair[i_pair][i] = i_lr1;
                  i_lr2_pair[i_pair][i] = i_lr2;
                  num_pair[i_pair]  += 1;
                }
              }
            }
          }
        }
      }else if( 1==plane_enabled_temp[i_plane1] && 0==plane_enabled_temp[i_plane2] ){
        for(int i_hit1=0; i_hit1<numhit[i_plane1]; i_hit1++){
          for(int i_lr1=0; i_lr1<2; i_lr1++){
            int i = num_pair[i_pair];
            if(i>=32) continue;
            i_hit1_pair[i_pair][i]= i_hit1;
            i_hit2_pair[i_pair][i]= -1;
            i_lr1_pair[i_pair][i] = i_lr1;
            i_lr2_pair[i_pair][i] = -1;
            num_pair[i_pair]  += 1;
          }
        }
      }else if ( 0==plane_enabled_temp[i_plane1] && 1==plane_enabled_temp[i_plane2] ){
        for(int i_hit2=0; i_hit2<numhit[i_plane2]; i_hit2++){
          for(int i_lr2=0; i_lr2<2; i_lr2++){
            int i = num_pair[i_pair];
            if(i>=32) continue;
            i_hit2_pair[i_pair][i]= i_hit2;
            i_hit1_pair[i_pair][i]= -1;
            i_lr2_pair[i_pair][i] = i_lr2;
            i_lr1_pair[i_pair][i] = -1;
            num_pair[i_pair]  += 1;
          }
        }
      }else{
        num_pair[i_pair] = 0;
      }
    }
    if(num_pair[i_pair]>=32) pair_enabled_temp[i_pair] = 0;
  }
  return 1;
}


int MWDCTracking::PairCombinationTempIncrement(){ // go to next combination
  // 1 : done
  // 0 : can not increment any more...
  for(int i_pair=0; i_pair<8; i_pair++){
    if(0<num_pair[i_pair]){
      if(i_pair_used_temp[i_pair]==(num_pair[i_pair]-1)){
        i_pair_used_temp[i_pair]=0; // (kuriagari) continue
      }else{
        i_pair_used_temp[i_pair]+=1;
        return 1;
      }
    }
  }
  return 0;
}


int MWDCTracking::SetPairCombinationTemp0(){ // go to next combination
  for(int i_pair=0; i_pair<8; i_pair++){
    if(1==pair_enabled_temp[i_pair]){
      i_pair_used_temp[i_pair]=0;
    }
  }
  return 0;
}



int MWDCTracking::LRCombinationTempIncrement(){
  // 1 : done
  // 0 : can not increment any more...
  for(int i_plane=0; i_plane<16; i_plane++){
    if(0==LR_fixed[i_plane]&&1==plane_enabled_temp[i_plane]){
      if(1==LR_used_temp[i_plane]){
        LR_used_temp[i_plane]=0; // kuriagari, continue
      }else{
        LR_used_temp[i_plane]=1;
        return 1;
      }
    }
  }
  return 0;
}

void MWDCTracking::SetHitCombinationTemp0(){
  for(int i_plane=0; i_plane<16; i_plane++){
    i_hit_used_temp[i_plane]=0;
  }
}

void MWDCTracking::SetLRCombinationTemp0(){
  for(int i_plane=0; i_plane<16; i_plane++){
    LR_used_temp[i_plane]=0;
  }
  return ;
}

void MWDCTracking::CopyHitCombination(int* array_src, int* array_tgt){
  for(int i_plane=0; i_plane<16; i_plane++){
    array_tgt[i_plane] = array_src[i_plane];
  }
  return ;
}

void MWDCTracking::CopyLRCombination(int* array_src, int* array_tgt){
  for(int i_plane=0; i_plane<16; i_plane++){
    array_tgt[i_plane] = array_src[i_plane];
  }
  return ;
}

void MWDCTracking::CopyPlaneEnabled(int* array_src, int* array_tgt){
  for(int i_plane=0; i_plane<16; i_plane++){
    array_tgt[i_plane] = array_src[i_plane];
  }
  return ;
}


//-----------------
//    Get****()

Float_t MWDCTracking::GetX(void){
  return x_track;
}

Float_t MWDCTracking::GetY(void){
  return y_track;
}

Float_t MWDCTracking::GetA(void){
  return a_track;
}

Float_t MWDCTracking::GetB(void){
  return b_track;
}


Float_t MWDCTracking::GetXError(void){
  return xerr_track;
}

Float_t MWDCTracking::GetYError(void){
  return yerr_track;
}

Float_t MWDCTracking::GetAError(void){
  return aerr_track;
}

Float_t MWDCTracking::GetBError(void){
  return berr_track;
}

Float_t MWDCTracking::GetTime_wMaxToT(int iplane){
  return time_wMaxToT[iplane];
}

int MWDCTracking::GetTrackingStatus(void){
  return tracking_status;
}

Float_t MWDCTracking::GetChi2(void){
  return chi2_track;
}

Float_t MWDCTracking::GetResidual(int i_plane){
  return residual_track[i_plane];
}

Float_t MWDCTracking::GetExclusiveResidual(int i_plane){
  return residual_track_ex[i_plane];
}

void MWDCTracking::GetLRCombination(int* lr_combi){
  for(int i_plane=0; i_plane<16; i_plane++){
    lr_combi[i_plane] = LR_used_track[i_plane];
  }
  return ;
}

void MWDCTracking::GetHitCombination(int* hit_combi){
  for(int i_plane=0; i_plane<16; i_plane++){
    hit_combi[i_plane] = i_hit_used_track[i_plane];
  }
  return ;
}

int MWDCTracking::GetNumEnabledPlane(void){
  int count=0;
  for(int i_plane=0; i_plane<16; i_plane++){
    if(1==plane_enabled_track[i_plane]){
      count++;
    }
  }
  return count;
}

int MWDCTracking::GetNumLRCombination(void){
  return (0x0001 << GetNumEnabledPlane()); // 2^n
}

int MWDCTracking::GetNumHitCombination(void){
  int count=1;
  for(int i_plane=0; i_plane<16; i_plane++){
    if(1==plane_enabled_track[i_plane]){
      count*=numhit[i_plane];
    }
  }
  return count;
}



int MWDCTracking::GetiHit(int i_plane){
  if(1==plane_enabled_track[i_plane]){
    return i_hit_used_track[i_plane];
  }
  return -1;
}

int MWDCTracking::GetLR(int i_plane){
  if(1==plane_enabled_track[i_plane]){
    return LR_used_track[i_plane];
  }
  return -1;
}

Float_t MWDCTracking::GetUsedWirepos(int i_plane){
  if(1==plane_enabled_track[i_plane]){
    int tmp_i_hit = i_hit_used_track[i_plane];
    return wirepos[i_plane][tmp_i_hit];
  }
  return -999.9;
}

Float_t MWDCTracking::GetUsedLength(int i_plane){
  if(1==plane_enabled_track[i_plane]){
    int tmp_i_hit = i_hit_used_track[i_plane];
    return length[i_plane][tmp_i_hit];
  }
  return -999.9;
}
