#include "TFragmentFinder.h"

//#define DEBUG_FragmentFinder

using namespace std;
using namespace G4Sol;

template<class Out>
TFragmentFinder<Out>::TFragmentFinder(const THyphiAttributes& attribut) : TDataProcessInterface<Out>("Fragment_Reco"), att(attribut)
{
  att._logger->info("TFragmentFinder::TFragmentFinder");

}

template<class Out>
TFragmentFinder<Out>::~TFragmentFinder() {}

template<class Out>
void TFragmentFinder<Out>::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }

template<class Out>
ReturnRes::InfoM TFragmentFinder<Out>::operator()(FullRecoEvent& RecoEvent, Out* OutTree)
{

  int result_finder = Exec(RecoEvent, OutTree);

  return SoftExit(result_finder);
}

template<class Out>
int TFragmentFinder<Out>::Exec(FullRecoEvent& RecoEvent, Out* OutTree) { return FinderFragment(RecoEvent); }

template<class Out>
ReturnRes::InfoM TFragmentFinder<Out>::SoftExit(int result_full) { return ReturnRes::Fine; }

template<class Out>
void TFragmentFinder<Out>::SelectHists()
{
  LocalHisto.hopt_1_1      = this->AnaHisto->CloneAndRegister(this->AnaHisto->hopt_1_1);
  LocalHisto.hopt_1_2      = this->AnaHisto->CloneAndRegister(this->AnaHisto->hopt_1_2);
  LocalHisto.hopt_1_3      = this->AnaHisto->CloneAndRegister(this->AnaHisto->hopt_1_3);
  LocalHisto.hopt_1_4      = this->AnaHisto->CloneAndRegister(this->AnaHisto->hopt_1_4);

  LocalHisto.hopt_2_1      = this->AnaHisto->CloneAndRegister(this->AnaHisto->hopt_2_1);
  LocalHisto.hopt_2_2      = this->AnaHisto->CloneAndRegister(this->AnaHisto->hopt_2_2);
  LocalHisto.hopt_2_3      = this->AnaHisto->CloneAndRegister(this->AnaHisto->hopt_2_3);
  LocalHisto.hopt_2_4      = this->AnaHisto->CloneAndRegister(this->AnaHisto->hopt_2_4);
}

template<class Out>
int TFragmentFinder<Out>::FinderFragment(FullRecoEvent& RecoEvent)
{
  StudyCaseSelector_Fr(att.StudyCase, Fragment_pdg);

  if(recons_from_FRS_MDC == 2)
    {
      //FragmentMDCTracksFinder(RecoEvent.DAF_results, Fragment_pdg, RecoEvent.FragmentTracks);
      LocalHisto.hopt_1_4->Fill(RecoEvent.FragmentTracks.size());
      return 0;
    }

  if(att.G4_simu)
    {
      RealFragmentFinder(RecoEvent.TrackDAFSim, Fragment_pdg, RecoEvent.FragmentTracks);
      LocalHisto.hopt_1_4->Fill(RecoEvent.FragmentTracks.size());

      for(int i=0; i<RecoEvent.FragmentTracks.size(); ++i)
        {
          LocalHisto.hopt_1_1->Fill(RecoEvent.FragmentTracks[i].GetMom().X() / RecoEvent.FragmentTracks[i].GetMom().Z());
          LocalHisto.hopt_1_2->Fill(RecoEvent.FragmentTracks[i].GetMom().Y() / RecoEvent.FragmentTracks[i].GetMom().Z());
          LocalHisto.hopt_1_3->Fill(RecoEvent.FragmentTracks[i].GetMom().Mag());
        }

      return 0;
    }

  if(RecoEvent.MWDCTracks.size() == 1 && RecoEvent.FragmentPID != -999)
    {
      for(int i=0; i<RecoEvent.FiberTrackCont["dft12"].size(); ++i)
        {
          //std::cout << "s4PID : " << s4hit->GetPID() << std::endl;
          OpticsMom *optics = new OpticsMom( RecoEvent.FiberTrackCont["dft12"][i], RecoEvent.MWDCTracks[0], RecoEvent.FragmentPID, att.optics_par, att.optics_s2z);
          if( !optics->IsValid() )
            {
              delete optics;
              continue;
            }
          if( RecoEvent.FiberTrackCont["dft12"][i]->IsBest() ) optics->SetBest();
          //optics->Print();
          LocalHisto.hopt_1_1->Fill(optics->GetA2Rec());
          LocalHisto.hopt_1_2->Fill(optics->GetB2Rec());
          LocalHisto.hopt_1_3->Fill(optics->GetMom());

          LocalHisto.hopt_2_1->Fill(optics->GetA2(), optics->GetA2Rec());
          LocalHisto.hopt_2_2->Fill(optics->GetB2(), optics->GetB2Rec());
          LocalHisto.hopt_2_3->Fill(optics->GetA2Rec() - optics->GetA2());
          LocalHisto.hopt_2_4->Fill(optics->GetB2Rec() - optics->GetB2());


          FragmentTrack tmp_FragmentTrack;

          tmp_FragmentTrack.SetPos(TVector3(optics->GetX2()*0.1, optics->GetY2()*0.1, optics->GetZ2()*0.1));
          tmp_FragmentTrack.SetMom(optics->GetMomV());

          std::vector<double> tmp_covmatrix = {1.e-8,
                                                  0., 1.e-8,
                                                  0.,    0., 1.e-8,
                                                  0.,    0.,    0., 1.e-8,
                                                  0.,    0.,    0.,    0., 1.e-8,
                                                  0.,    0.,    0.,    0.,    0., 1.e-8}; //Change!
          tmp_FragmentTrack.SetCovMatrix(tmp_covmatrix); //CHECK Change!
          tmp_FragmentTrack.SetTOT(RecoEvent.FiberTrackCont["dft12"][i]->GetTOT());
          tmp_FragmentTrack.SetTime(RecoEvent.FiberTrackCont["dft12"][i]->GetTime());
          tmp_FragmentTrack.SetChi2NDF(RecoEvent.FiberTrackCont["dft12"][i]->GetChi2());
          tmp_FragmentTrack.SetIsBest(optics->IsBest());
          tmp_FragmentTrack.SetPID(optics->GetPID());
          tmp_FragmentTrack.SetTrackID(i);

          RecoEvent.FragmentTracks.emplace_back(tmp_FragmentTrack);
        }
    }

  LocalHisto.hopt_1_4->Fill(RecoEvent.FragmentTracks.size());

  return 0;
}


template <class Out>
void TFragmentFinder<Out>::StudyCaseSelector_Fr(std::string StudyCase, int& Fragment_pdg)
{
  if(StudyCase.compare("H3L") == 0)
    {
      Fragment_pdg = He3_pdg;
      recons_from_FRS_MDC = 1;
    }
  else if(StudyCase.compare("H4L") == 0)
    {
      Fragment_pdg = He4_pdg;
      recons_from_FRS_MDC = 1;
    }
  else if(StudyCase.compare("nnL") == 0)
    {
      Fragment_pdg = deuteron_pdg;
      recons_from_FRS_MDC = 1;
    }
  else if(StudyCase.compare("lambda") == 0)
    {
      Fragment_pdg = proton_pdg;
      recons_from_FRS_MDC = 2;
    }
  else if(StudyCase.compare("background_H3L") == 0)
    {
      Fragment_pdg = He3_pdg;
      recons_from_FRS_MDC = 1;
    }
  else if(StudyCase.compare("background_H4L") == 0)
    {
      Fragment_pdg = He4_pdg;
      recons_from_FRS_MDC = 1;
    }
  
  return;
}


template <class Out>
void TFragmentFinder<Out>::RealFragmentFinder(std::unordered_map<int, std::vector<std::vector<SimHit> > >& TrackDAFSim,
                                               int& fragment_pdg, std::vector<FragmentTrack>& RealFragmentTracks)
{
  std::unordered_map<int, std::vector<std::vector<SimHit> > >::iterator itr;
  for(itr = TrackDAFSim.begin(); itr != TrackDAFSim.end(); ++itr)
    {
      size_t iDetFirst = -1;

      for(size_t iDet = G4Sol::FiberD4_v; iDet < G4Sol::FiberD5_v; ++iDet)
        {
          if(itr->second[iDet].size() == 0)
            continue;

          iDetFirst = iDet;
          break;
        }

      if(iDetFirst == -1)
        continue;

      if(itr->second[iDetFirst][0].pdg == fragment_pdg)
        {

          TVector3 tmp_pos = TVector3(itr->second[iDetFirst][0].hitX, itr->second[iDetFirst][0].hitY, itr->second[iDetFirst][0].hitZ);
          TVector3 tmp_mom = TVector3(itr->second[iDetFirst][0].momX, itr->second[iDetFirst][0].momY, itr->second[iDetFirst][0].momZ);
          std::vector<double> tmp_covmatrix = {1.e-2,
                                                  0., 1.e-2,
                                                  0.,    0., 1.e-2,
                                                  0.,    0.,    0., 1.e-3,
                                                  0.,    0.,    0.,    0., 1.e-3,
                                                  0.,    0.,    0.,    0.,    0., 1.e-3}; //Change!
          

          FragmentTrack tmp_FragmentTrack;
          tmp_FragmentTrack.SetPos(tmp_pos);
          tmp_FragmentTrack.SetMom(tmp_mom);
          tmp_FragmentTrack.SetCovMatrix(tmp_covmatrix);
          tmp_FragmentTrack.SetTOT(1.); //CHECK Change!
          tmp_FragmentTrack.SetTime(1.); //CHECK Change!
          tmp_FragmentTrack.SetChi2NDF(1.); //CHECK Change!
          tmp_FragmentTrack.SetIsBest(true);
          tmp_FragmentTrack.SetPID(fragment_pdg);
          tmp_FragmentTrack.SetTrackID(itr->first);

          RealFragmentTracks.emplace_back(tmp_FragmentTrack);
        }
    }

  return;
}

/*
template <class Out>
void TFragmentFinder<Out>::FragmentMDCTracksFinder(std::unordered_map<int, ResSolDAF>& DAF_results, int& fragment_pdg,
                                                    std::vector<FragmentTrack>& FragmentMDCTracks)
{
  int temp_charge = TDatabasePDG::Instance()->GetParticle(fragment_pdg)->Charge()/3.;

  std::unordered_map<int, ResSolDAF>::iterator itr;
  for(itr = DAF_results.begin(); itr != DAF_results.end(); ++itr)
    {
      if((itr->second.charge == temp_charge) && (itr->second.Ncentral > att.KF_NbCentralCut))
        {
          TVector3 tmp_pos = TVector3(itr->second.posX, itr->second.posY, itr->second.posZ);
          TVector3 tmp_mom = TVector3(itr->second.momX, itr->second.momY, itr->second.momZ);
          std::vector<double> tmp_covmatrix = {1.e-8,
                                                  0., 1.e-8,
                                                  0.,    0., 1.e-8,
                                                  0.,    0.,    0., 1.e-8,
                                                  0.,    0.,    0.,    0., 1.e-8,
                                                  0.,    0.,    0.,    0.,    0., 1.e-8}; //Change!

          FragmentTrack tmp_FragmentTrack;
          tmp_FragmentTrack.SetPos(tmp_pos);
          tmp_FragmentTrack.SetMom(tmp_mom);
          tmp_FragmentTrack.SetCovMatrix(tmp_covmatrix);
          tmp_FragmentTrack.SetTOT(1.); //CHECK Change!
          tmp_FragmentTrack.SetTime(1.); //CHECK Change!
          tmp_FragmentTrack.SetChi2NDF(1.); //CHECK Change!
          tmp_FragmentTrack.SetIsBest(true);
          tmp_FragmentTrack.SetPID(fragment_pdg);

          FragmentMDCTracks.emplace_back(tmp_FragmentTrack);
        }
    }

  return;
}
*/

template class TFragmentFinder<MCAnaEventG4Sol>;
template class TFragmentFinder<Ana_WasaEvent>;
