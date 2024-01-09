#ifndef TKALMANFILTERDAFPID
#define TKALMANFILTERDAFPID

#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "TDataProcess.h"
//#include "Ana_EventNew.hh"
//#include "Ana_Event/Ana_EventNew_v10.hh"
//#include "Ana_Event/Ana_EventNew_v16.hh"

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "THyphiAttributes.h"
//#include "Genfit/HypLSLTrackFieldMapRep_v2.h"
//#include "Genfit_v3/core/GFTools.h"
//#include "Genfit_v3/core/GFDaf.h"
//#include "Genfit_v3/core/GFTrack.h"
//#include "Genfit_v3/RKTrackRep/RKTrackRep.h"

//#include "GFTools.h"
#include "AbsKalmanFitter.h"
#include "DAF.h"
#include "EventDisplay.h"
#include "KalmanFitter.h"
#include "KalmanFitterRefTrack.h"
#include "RKTrackRep.h"
#include "Track.h"

#include "TMath.h"
#include "TString.h"
#include "TCutG.h"

//#include "GFRaveVertexFactory.h"

//#include "Math/ProbFunc.h"
template<class Out>
using TDataProcessInterface = TDataProcess<FullRecoEvent, Out>;
// namespace MathKalman
// {
//   struct Prob {
//     double operator() (double chi2,double ndf)
//     {
//       if (ndf <= 0)
//  return 0; // Set CL to zero in case ndf<=0

//       if (chi2 <= 0)
//  {
//    if (chi2 < 0)
//      return 0;
//    else
//      return 1;
//  }

//       return ::ROOT::Math::chisquared_cdf_c(chi2,ndf);
//     }
//   };

//}


class CutMomBeta 
{

private :
  // mom-beta curve : beta = 1/sqrt(mass*mass/mom/mom + 1)
  const double mass; 
  double marging;
  double marging_limit;
  double mass_up2;
  double mass_down2;
  const bool norm;
  const bool transl_beta;
  inline double CalBetaDown(double mom_test) const
  {return transl_beta == false ? 1./TMath::Sqrt(mass_up2/mom_test/mom_test+1.) : 1./TMath::Sqrt(mass*mass/mom_test/mom_test+1.) - marging_limit;}
  inline double CalBetaUp(double mom_test) const
  {return transl_beta == false ? 1./TMath::Sqrt(mass_down2/mom_test/mom_test+1.) : 1./TMath::Sqrt(mass*mass/mom_test/mom_test+1.) + marging;}

  
public :
  double operator () (double mom_test,double beta_test) const 
  {
    
    double beta_down_marge = CalBetaDown(mom_test);//1./TMath::Sqrt(mass_up2/mom_test/mom_test+1.);
    double beta_up_marge = CalBetaUp(mom_test);//1./TMath::Sqrt(mass_down2/mom_test/mom_test+1.);

    double beta_th = 1./TMath::Sqrt(mass*mass/mom_test/mom_test+1.);

    double diff = beta_test-beta_th;
    double sigma = 0.;

    //std::cout<<"CutMomBeta("<<mass<<","<<marging<<"): mom="<<mom_test<<" beta="<<beta_test<<" Bup="<<beta_up_marge<<" Bdown:"<<beta_down_marge<<" Bth:"<<beta_th<<" diff:"<<diff<<" ";

    if(diff>0.)
      {
        sigma = beta_up_marge-beta_th;
        //double arg = -0.5*diff*diff ;///beta_up_marge/beta_up_marge;
        //std::cout<<" sigma:"<<sigma<<" prob:"<<TMath::Gaus(beta_test,beta_th,sigma,norm)<<std::endl;
        return TMath::Gaus(beta_test,beta_th,sigma,norm);//diff*diff/beta_up_marge;
      }
    else
      {
        sigma = beta_th-beta_down_marge;
        //double arg = -0.5*diff*diff; //beta_down_marge/beta_down_marge;
        //std::cout<<" sigma:"<<sigma<<" prob:"<<TMath::Gaus(beta_test,beta_th,sigma,norm)<<std::endl;
        return TMath::Gaus(beta_test,beta_th,sigma,norm);
      }
    // if(beta_down_marge < beta_test && beta_test < beta_up_marge)
    //   return true;
    // else
    //   return false;
  }
  double GetDiff( double mom_test, double beta_test) const
  {
    double beta_th = 1./TMath::Sqrt(mass*mass/mom_test/mom_test+1.);
    return beta_test-beta_th;    
  }
  double GetChi2( double mom_test, double beta_test) const
  {
    double beta_down_marge = CalBetaDown(mom_test);//1./TMath::Sqrt(mass*mass*(1.+marging)*(1.+marging)/mom_test/mom_test+1.);
    double beta_up_marge = CalBetaUp(mom_test);//1./TMath::Sqrt(mass*mass*(1.-marging)*(1.-marging)/mom_test/mom_test+1.);

    double beta_th = 1./TMath::Sqrt(mass*mass/mom_test/mom_test+1.);

    double diff = beta_test-beta_th;
    if(diff>0.)
      return diff*diff/beta_up_marge/beta_up_marge;
    else
      return diff*diff/beta_down_marge/beta_down_marge;
  }
  double GetProb(double mom_test, double beta_test) const
  {
    double beta_th = 1./TMath::Sqrt(mass*mass/mom_test/mom_test+1.);

    double beta_down_marge = CalBetaDown(mom_test);//1./TMath::Sqrt(mass*mass*(1.+marging)*(1.+marging)/mom_test/mom_test+1.);
    double beta_up_marge = CalBetaUp(mom_test);//1./TMath::Sqrt(mass*mass*(1.-marging)*(1.-marging)/mom_test/mom_test+1.);
    
    double diff = beta_test-beta_th;
    double sigma = 0.;
    if(diff>0)
      sigma = beta_up_marge-beta_th;
    else
      sigma = beta_th-beta_down_marge;

    return TMath::Gaus(beta_test,beta_th,sigma,true);//TMath::Exp(-0.5*diff*diff/sigma/sigma);
  }


  CutMomBeta():mass(1.),marging(0.),norm(true),transl_beta(false) {}
  CutMomBeta(double M,double marge,bool nn=true,bool translbeta =false):mass(M),marging(marge),marging_limit(marge),norm(nn),transl_beta(translbeta)
  {
    mass_up2=TMath::Sq(mass*(1.+marging));
    mass_down2=TMath::Sq(mass*(1.-marging));
  }
  ~CutMomBeta() {}
  
  void SetMarging(double marge) 
  {
    marging=marge;     
    mass_up2=TMath::Sq(mass*(1.+marging));
    mass_down2=TMath::Sq(mass*(1.-marging));
  }
  void SetMargingLimit(double margeL)
  {
    marging_limit=margeL;
    mass_up2=TMath::Sq(mass*(1.+marging_limit));
  }
  double GetMass() const { return mass;  }
};



template<class Out>
class TKalmanFilter_DAF_PID final : public TDataProcessInterface<Out>
{
public:
  const THyphiAttributes& att;

  explicit TKalmanFilter_DAF_PID(const THyphiAttributes& attr);
  ~TKalmanFilter_DAF_PID() final;

  // int Init(Ana_Hist* h);
  ReturnRes::InfoM operator()(FullRecoEvent& RecoEvent, Out* OutTree) final;
  void InitMT() final;

private:
  int Exec(FullRecoEvent& RecoEvent, Out* OutTree) final;
  ReturnRes::InfoM SoftExit(int) final;
  void SelectHists() final;

  int Kalman_Filter_FromTrack(FullRecoEvent& RecoEvent);

  int Nb_CentralCut;
  int Nb_MiniFiberCut;

  std::map<int,CutMomBeta*> CutPID;
  double PID_maxMom_piPlus = 10000.;//1.5;
  double PID_maxMom_kaonPlus = 10000.;// 2.5;
  double PID_maxMom_proton = 10000.;//3.5;

  int ProbPIDAssign_Pos(double momenta, double beta);
  int ProbPIDAssign_Neg(double momenta, double beta);

  // genfit::DAF* Fitter;
  genfit::AbsKalmanFitter* Fitter;
  genfit::AbsKalmanFitter* Fitter_pid;
  genfit::AbsKalmanFitter* Fitter_rescue;
  genfit::Track* Vtracks;
  genfit::RKTrackRep* rep;
  genfit::EventDisplay* display;
  // genfit::GFRaveVertexFactory* vertexFactory;

  std::vector<genfit::DetPlane*> list_Plane;

  TVector3 Plane_time;

  int Nb_event = 0;
  struct LocalHists
  {
    TH1I* h_stats;
    TH2I* h_statsLess3Mes;
    TH2I* h_statsInvalid;

    TH1F* h_pv;
    TH1F* h_chi2;
    TH1F* hd_chi[2];
    TH1F* hd_pv[2];

    TH1F* h_Path;
    TH1F* h_MeanPath;

    TH1F* h_beta;
    TH1F* h_beta2;
    TH1F* h_beta3;

    TH1F* h_Mass_All;
    TH1F* h_Mass_All2;
    TH1F* h_Mass_All3;

    TH2F* h_Mass_charge_All;
    TH2F* h_Mass_charge_All2;
    TH2F* h_Mass_charge_All3;

    TH2F* h_beta_mom;
    TH2F* h_beta_mom2;
    TH2F* h_beta_mom3;
    TH2F* h_beta_momcharge;
    TH2F* h_beta_momcharge2;
    TH2F* h_beta_momcharge3;

    TH2F* h_pv_mom;
    TH2F* h_pv_beta;
    TH2F* h_pv_mass;

    TH2F* h_path_tof;

    TH2F* h_mom_tof_cut;
    TH2F* h_path_mom_cut;
    TH2F* h_path_tof_cut;

    TH1F* h_Mass[4];
    TH1F* h_chi2_particle[4];
    TH1F* h_pv_particle[4];

    TH2F* h_mom_res[5];
    TH1D* h_ResPull[5][10];
    TH2F* h_ResPull_normal[5][10];

    TH2F* h_total_dE;

    TH1F* h_ResFiber[9];
    TH1F* h_ResMiniFiber[6];
    TH1F* h_ResMDC[17][3];
    TH1F* h_ResPSCE[2];
    
  };
  LocalHists LocalHisto;
};

#endif
