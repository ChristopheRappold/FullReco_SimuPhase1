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

  if(RecoEvent.MWDCTracks.size() == 1 && RecoEvent.FragmentPID != -999){
    for(int i=0; i<RecoEvent.FiberTrackCont["dft12"].size(); ++i){
      //std::cout << "s4PID : " << s4hit->GetPID() << std::endl;
      OpticsMom *optics = new OpticsMom( RecoEvent.FiberTrackCont["dft12"][i], RecoEvent.MWDCTracks[0], RecoEvent.FragmentPID, att.optics_par, att.optics_s2z);
      if( !optics->IsValid() ){
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
      //tmp_FragmentTrack.SetCovMatrix(); CHECK
      tmp_FragmentTrack.SetTOT(RecoEvent.FiberTrackCont["dft12"][i]->GetTOT());
      tmp_FragmentTrack.SetTime(RecoEvent.FiberTrackCont["dft12"][i]->GetTime());
      tmp_FragmentTrack.SetChi2NDF(RecoEvent.FiberTrackCont["dft12"][i]->GetChi2());
      tmp_FragmentTrack.SetIsBest(optics->IsBest());
      tmp_FragmentTrack.SetPID(optics->GetPID());

      RecoEvent.FragmentTracks.emplace_back(tmp_FragmentTrack);
    }
  }

  LocalHisto.hopt_1_4->Fill(RecoEvent.FragmentTracks.size());

  return 0;
}

template class TFragmentFinder<MCAnaEventG4Sol>;
template class TFragmentFinder<Ana_WasaEvent>;
