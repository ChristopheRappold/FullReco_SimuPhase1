#include "TWASAFinder.h"

//#define DEBUG_WASAFINDER

using namespace std;
using namespace G4Sol;

template<class Out>
TWASAFinder<Out>::TWASAFinder(const THyphiAttributes& attribut) : TDataProcessInterface<Out>("WASA_Reco"), att(attribut)
{
  att._logger->info("TWASAFinder::TWASAFinder");

}

template<class Out>
TWASAFinder<Out>::~TWASAFinder() {}

template<class Out>
void TWASAFinder<Out>::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }

template<class Out>
ReturnRes::InfoM TWASAFinder<Out>::operator()(FullRecoEvent& RecoEvent, Out* OutTree)
{

  int result_finder = Exec(RecoEvent, OutTree);

  return SoftExit(result_finder);
}

template<class Out>
int TWASAFinder<Out>::Exec(FullRecoEvent& RecoEvent, Out* OutTree) { return FinderWASA(RecoEvent); }

template<class Out>
ReturnRes::InfoM TWASAFinder<Out>::SoftExit(int result_full) { return ReturnRes::Fine; }

template<class Out>
void TWASAFinder<Out>::SelectHists()
{
  LocalHisto.h23_1      = this->AnaHisto->CloneAndRegister(this->AnaHisto->h23_1);
  LocalHisto.h23_2      = this->AnaHisto->CloneAndRegister(this->AnaHisto->h23_2);
  LocalHisto.h24_1      = this->AnaHisto->CloneAndRegister(this->AnaHisto->h24_1);
  LocalHisto.h24_2      = this->AnaHisto->CloneAndRegister(this->AnaHisto->h24_2);
  LocalHisto.h24_2_1    = this->AnaHisto->CloneAndRegister(this->AnaHisto->h24_2_1);
  LocalHisto.h24_2_2    = this->AnaHisto->CloneAndRegister(this->AnaHisto->h24_2_2);
  for(size_t i = 0; i < 17; ++i)
    {
      LocalHisto.h24_3[i]      = this->AnaHisto->CloneAndRegister(this->AnaHisto->h24_3[i]);
      LocalHisto.h24_4[i]      = this->AnaHisto->CloneAndRegister(this->AnaHisto->h24_4[i]);
      LocalHisto.h24_2_3[i]    = this->AnaHisto->CloneAndRegister(this->AnaHisto->h24_2_3[i]);
      LocalHisto.h24_2_4[i]    = this->AnaHisto->CloneAndRegister(this->AnaHisto->h24_2_4[i]);
      LocalHisto.hmdc_2_2[i]   = this->AnaHisto->CloneAndRegister(this->AnaHisto->hmdc_2_2[i]);
      LocalHisto.hmdc_3_2[i]   = this->AnaHisto->CloneAndRegister(this->AnaHisto->hmdc_3_2[i]);
    }
}

template<class Out>
int TWASAFinder<Out>::FinderWASA(FullRecoEvent& RecoEvent)
{
  FiberAnalyzer* fiberana = new FiberAnalyzer();
  std::string track_type = "mft12";

  // Fiber PSB ////
  for(int i=0; i<RecoEvent.FiberTrackCont[track_type].size(); ++i)
    {
      for(int j=0; j<(int)RecoEvent.ListHits[G4Sol::PSCE].size(); ++j)
        {
          if(RecoEvent.FiberTrackCont[track_type][i]->GetNlayer()>4 && RecoEvent.FiberTrackCont[track_type][i]->GetChi2()>30) continue;
          double a_fiber = RecoEvent.FiberTrackCont[track_type][i]->GetA();
          double b_fiber = RecoEvent.FiberTrackCont[track_type][i]->GetB();
          double x_fiber = RecoEvent.FiberTrackCont[track_type][i]->GetX();
          double y_fiber = RecoEvent.FiberTrackCont[track_type][i]->GetY();
          //double phi_fiber = atan2(b_fiber, a_fiber);
          double phi_psb   = GetPSB_Phi(RecoEvent.ListHits[G4Sol::PSCE][j]->getHitId())    + att.psb_rot_z*Deg2Rad;
          double r_psb     = GetPSB_R(RecoEvent.ListHits[G4Sol::PSCE][j]->getHitId());
          double z_psb     = RecoEvent.ListHits[G4Sol::PSCE][j]->getRawHitCoords()[1]*10.; //in mm

          double par_a = pow(a_fiber, 2) + pow(b_fiber, 2);
          double par_b = a_fiber * (x_fiber - att.psb_pos_x)+ b_fiber * (y_fiber - att.psb_pos_y);
          double par_c = pow(x_fiber - att.psb_pos_x, 2) + pow(y_fiber - att.psb_pos_y, 2) - pow(r_psb, 2);
          double z_fiber  = (-par_b + sqrt( pow(par_b,2) - par_a * par_c)) / par_a;

          double fiber_x_buf = x_fiber + a_fiber * z_fiber - att.psb_pos_x;
          double fiber_y_buf = y_fiber + b_fiber * z_fiber - att.psb_pos_y;
          double phi_fiber = atan2(fiber_y_buf, fiber_x_buf);

          LocalHisto.h23_1->Fill(phi_fiber, phi_psb);
          LocalHisto.h23_2->Fill(z_fiber - att.psb_pos_z, z_psb);

          if( fabs(fiberana->CalcPhiDif(phi_psb, phi_fiber)) < att.cut_psb_phi )
            {
              if(fabs( (z_fiber - att.psb_pos_z) - z_psb)<att.cut_psb_z)
                {
                  if(RecoEvent.FiberTrackCont[track_type][i]->IsFlagPSB())
                    {
                      //double phi_psb_buf = RecoEvent.FiberTrackCont[track_type][i]->GetPSBHit()->GetPhi();
                      if( fabs( fiberana->CalcPhiDif(phi_psb, phi_fiber) ) < fabs( RecoEvent.FiberTrackCont[track_type][i]->GetPSBDifPhi() ) )
                        {
                          RecoEvent.FiberTrackCont[track_type][i]->SetSegPSB(RecoEvent.ListHits[G4Sol::PSCE][j]->getHitId());
                          //RecoEvent.FiberTrackCont[track_type][i]->SetPSBHit(PSBHitCont[j]);
                          RecoEvent.FiberTrackCont[track_type][i]->SetPSBDifZ(z_psb - (z_fiber - att.psb_pos_z));
                          RecoEvent.FiberTrackCont[track_type][i]->SetPSBDifPhi(fiberana->CalcPhiDif(phi_psb, phi_fiber));
                        }
                    }
                  else
                    {
                      RecoEvent.FiberTrackCont[track_type][i]->SetFlagPSB();
                      RecoEvent.FiberTrackCont[track_type][i]->SetSegPSB(RecoEvent.ListHits[G4Sol::PSCE][j]->getHitId());
                      //RecoEvent.FiberTrackCont[track_type][i]->SetPSBHit(PSBHitCont[j]);
                      RecoEvent.FiberTrackCont[track_type][i]->SetPSBDifZ(z_psb - (z_fiber - att.psb_pos_z));
                      RecoEvent.FiberTrackCont[track_type][i]->SetPSBDifPhi(fiberana->CalcPhiDif(phi_psb, phi_fiber));
                    }
                }
            }
        }
    }

  // Fiber MDC ////
  std::vector<TrackHit*> TrackHitCont;


  for(int i=0; i<RecoEvent.FiberTrackCont[track_type].size(); ++i)
    {
      TrackHit *track_hit = new TrackHit();
      for(int j=0; j<17; ++j)
        {
          for(int k=0; k<(int)RecoEvent.ListHits[G4Sol::MG01+j].size(); ++k)
            {
              if(RecoEvent.FiberTrackCont[track_type][i]->GetNlayer()>4 && RecoEvent.FiberTrackCont[track_type][i]->GetChi2()>30) continue;
              //if(par->flag_genfit_psbonly && !FiberTrackCont[track_type][i]->IsFlagPSB()) continue;
              if( (!RecoEvent.FiberTrackCont[track_type][i]->IsFlagPSB()) && (!RecoEvent.FiberTrackCont[track_type][i]->IsFlagPSFE()) ) continue;
              if(RecoEvent.FiberTrackCont[track_type][i]->IsFlagPSB() ) track_hit->SetFlagPSB();
              if(RecoEvent.FiberTrackCont[track_type][i]->IsFlagPSFE()) track_hit->SetFlagPSFE();
              //if(MDCHitCont[j][k]->GetDpf()==0) continue;
              //if(MDCHitCont[j][k]->GetDpf()!=0) continue;
              //if(MDCHitCont[j][k]->GetDpf()!=0 && j!=13) continue;
              double a_fiber = RecoEvent.FiberTrackCont[track_type][i]->GetA();
              double b_fiber = RecoEvent.FiberTrackCont[track_type][i]->GetB();
              double x_fiber = RecoEvent.FiberTrackCont[track_type][i]->GetX();
              double y_fiber = RecoEvent.FiberTrackCont[track_type][i]->GetY();
              //double phi_fiber = atan2(b_fiber, a_fiber);
              //double phi_mdc   = MDCHitCont[j][k]->GetPhi();

              double mdc_edge_1_x = RecoEvent.ListHits[G4Sol::MG01+j][k]->getRawHitCoords()[0]*10.;
              double mdc_edge_1_y = RecoEvent.ListHits[G4Sol::MG01+j][k]->getRawHitCoords()[1]*10.;
              double mdc_edge_1_z = RecoEvent.ListHits[G4Sol::MG01+j][k]->getRawHitCoords()[2]*10.;
              double mdc_edge_2_x = RecoEvent.ListHits[G4Sol::MG01+j][k]->getRawHitCoords()[3]*10.;
              double mdc_edge_2_y = RecoEvent.ListHits[G4Sol::MG01+j][k]->getRawHitCoords()[4]*10.;
              double mdc_edge_2_z = RecoEvent.ListHits[G4Sol::MG01+j][k]->getRawHitCoords()[5]*10.;
              TVector3 pos_mdc(mdc_edge_1_x, mdc_edge_1_y, mdc_edge_1_z);
              TVector3 ang_mdc(mdc_edge_2_x-mdc_edge_1_x, mdc_edge_2_y-mdc_edge_1_y, mdc_edge_2_z-mdc_edge_1_z);
              TVector3 pos_fiber(x_fiber + a_fiber*mdc_edge_1_z, y_fiber + b_fiber*mdc_edge_1_z, 0.);
              TVector3 ang_fiber(a_fiber, b_fiber, 1.);
              //pos_mdc.Print();
              double dist = 9999.;
              TVector3 buf_vertex = fiberana->GetVertexPoint( pos_fiber, pos_mdc, ang_fiber, ang_mdc, dist);
              double buf_vtz = buf_vertex.z() + mdc_edge_1_z;
              if(buf_vtz < mdc_edge_1_z || buf_vtz > mdc_edge_2_z) continue;
              //std::cout << "vtz : " << buf_vtz << std::endl;
              double fiber_x_vtx = x_fiber + a_fiber * buf_vtz;
              double fiber_y_vtx = y_fiber + b_fiber * buf_vtz;
              double phi_fiber = atan2(fiber_y_vtx, fiber_x_vtx);
              double mdc_x_vtx = mdc_edge_1_x + (mdc_edge_2_x-mdc_edge_1_x)/(mdc_edge_2_z-mdc_edge_1_z)*(buf_vtz-mdc_edge_1_z);
              double mdc_y_vtx = mdc_edge_1_y + (mdc_edge_2_y-mdc_edge_1_y)/(mdc_edge_2_z-mdc_edge_1_z)*(buf_vtz-mdc_edge_1_z);
              double phi_mdc = atan2(mdc_y_vtx, mdc_x_vtx);
              //if(fabs(phi_mdc)>0.1) continue;
              //std::cout << "phi_mdc_c : " << phi_mdc_c << std::endl;
              //std::cout << "phi_mdc : " << phi_mdc << std::endl;
              LocalHisto.h24_1->Fill(   phi_fiber, phi_mdc);
              LocalHisto.h24_3[j]->Fill(phi_fiber, phi_mdc);
              LocalHisto.h24_2_1->Fill(   fiberana->CalcPhiDif(phi_mdc, phi_fiber) );
              LocalHisto.h24_2_3[j]->Fill(fiberana->CalcPhiDif(phi_mdc, phi_fiber) );
              if( fabs( fiberana->CalcPhiDif(phi_fiber, phi_mdc) ) < att.cut_phi_fm )
                {
                  LocalHisto.h24_2->Fill(   phi_fiber, phi_mdc);
                  LocalHisto.h24_4[j]->Fill(phi_fiber, phi_mdc);
                  LocalHisto.h24_2_2->Fill(   fiberana->CalcPhiDif(phi_mdc, phi_fiber) );
                  LocalHisto.h24_2_4[j]->Fill(fiberana->CalcPhiDif(phi_mdc, phi_fiber) );
                  //LocalHisto.hmdc_2_2[j]->Fill(MDCHitCont[j][k]->GetDt());
                  LocalHisto.hmdc_3_2[j]->Fill(RecoEvent.ListHits[G4Sol::MG01+j][k]->getRawHitCoords()[6]*10.); //in mm
                  //MDCHitCont[j][k]->SetDif(fiberana->CalcPhiDif(phi_mdc, phi_fiber) );
                  track_hit->SetMDCdif(j, RecoEvent.ListHits[G4Sol::MG01+j][k]->getHitId(), fabs( fiberana->CalcPhiDif(phi_fiber, phi_mdc) ));
                }
            }
        }

      //track_hit->SetMDCLayHitCont(); //CHECK Maybe include later
      track_hit->DeleteDupMDC();
      //track_hit->SetDidCh(); //CHECK Maybe include later

      for( auto v : RecoEvent.FiberTrackCont[track_type][i]->GetContHit() )
        track_hit->AddFiber((v->GetDet()-3)*3 + v->GetLay(), v->GetClFib());

      if(track_hit->IsFlagPSB())       track_hit->AddPSB(RecoEvent.FiberTrackCont[track_type][i]->GetSegPSB());
      else if(track_hit->IsFlagPSFE()) track_hit->AddPSFE(RecoEvent.FiberTrackCont[track_type][i]->GetSegPSFE());

      track_hit->SetTrack(RecoEvent.FiberTrackCont[track_type][i], att.Target_PositionZ);
      TrackHitCont.emplace_back(track_hit);
    }

if(att.flag_dup_trackhit)     TrackHitCont = fiberana->DeleteDupTrackHit(TrackHitCont);
if(att.flag_dup_trackhit_mdc) TrackHitCont = fiberana->DeleteDupTrackHitMDC(TrackHitCont);
if(att.flag_trackhit_inclusive) TrackHitCont = fiberana->DeleteInclusiveTrackHit(TrackHitCont);


for(int i = 0; i < TrackHitCont.size(); ++i)
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
                //continue; //CHECK Remove when FiberID is same for data and simu
                tempSetHit[G4Sol::MiniFiberD1_x + j] = k;
                InfoPar tmp_infopar;
                tmp_infopar.pdg    = -211;
                //tmp_infopar.momX   = hit.MomX;
                //tmp_infopar.momY   = hit.MomY;
                //tmp_infopar.momZ   = hit.MomZ;
                //tmp_infopar.mass   = hit.Mass;
                //tmp_infopar.Eloss  = hit.Energy;
                //tmp_infopar.time   = hit.Time;
                //tmp_infopar.length = hit.TrackLength;
                tempSetInfo[G4Sol::MiniFiberD1_x + j] = tmp_infopar;
                flag_found = true;
              }
          }

        if(!flag_found) printf("Error: %d-layer Fiber Hit not found in WasaFinder\n", j);
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
                //tmp_infopar.momX   = hit.MomX;
                //tmp_infopar.momY   = hit.MomY;
                //tmp_infopar.momZ   = hit.MomZ;
                //tmp_infopar.mass   = hit.Mass;
                //tmp_infopar.Eloss  = hit.Energy;
                //tmp_infopar.time   = hit.Time;
                //tmp_infopar.length = hit.TrackLength;
                tempSetInfo[G4Sol::MG01 + j] = tmp_infopar;
                flag_found = true;
              }
          }

        if(!flag_found)
          printf("Error: %d-layer MDC Hit not found in WasaFinder\n", j);
      }

    //PSBHit
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
                //tmp_infopar.Eloss  = hit.Energy;
                //tmp_infopar.time   = hit.Time;
                //tmp_infopar.length = hit.TrackLength;
                tempSetInfo[G4Sol::PSCE] = tmp_infopar;

                flag_found = true;
              }
          }

        if(!flag_found)
          std::cout << "Error: PSB Hit not found in WasaFinder\n";
      }

    //PSFEHit
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
                //tmp_infopar.Eloss  = hit.Energy;
                //tmp_infopar.time   = hit.Time;
                //tmp_infopar.length = hit.TrackLength;
                tempSetInfo[G4Sol::PSFE] = tmp_infopar;

                flag_found = true;
              }
          }

        if(!flag_found)
          std::cout << "Error: PSFE Hit not found in WasaFinder\n";
      }

    RecoEvent.TrackDAF.insert(std::make_pair(i, tempSetHit));
    RecoEvent.TrackInfo.insert(std::make_pair(i, tempSetInfo));


    InfoInit tempInit;
    tempInit.charge = -1;
    tempInit.posX = TrackHitCont[i]->GetTrackX() *0.1;
    tempInit.posY = TrackHitCont[i]->GetTrackY() *0.1;
    tempInit.posZ = TrackHitCont[i]->GetTrackZ() *0.1;

    double tmp_momZ = 1.;
    tempInit.momX = tmp_momZ*TrackHitCont[i]->GetTrackA();
    tempInit.momY = tmp_momZ*TrackHitCont[i]->GetTrackB();
    tempInit.momZ = tmp_momZ;

    RecoEvent.TrackDAFInit.insert(std::make_pair(i, tempInit));
  }

  return 0;
}

template<class Out>
double TWASAFinder<Out>::GetPSB_R(int _seg)
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


template<class Out>
double TWASAFinder<Out>::GetPSB_Phi(int _seg)
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


template class TWASAFinder<MCAnaEventG4Sol>;
template class TWASAFinder<Ana_WasaEvent>;
