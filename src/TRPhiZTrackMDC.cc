#include "TRPhiZTrackMDC.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"
#include "FullRecoEvent.hh"
#include "FullRecoEventZMQ.hh"
#include "KalmanFittedStateOnPlane.h"
#include "KalmanFitterInfo.h"
#include "PlanarMeasurement.h"
#include "ReturnRes.hh"
#include "StateOnPlane.h"

#include "TVector3.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "TComplex.h"

#include <numeric>
#include <set>
#include <sstream>
#include <tuple>
#include <list>
#include <map>
#include <unordered_map>

#include "TGraphErrors.h"
#include "TFitResult.h"

using namespace std;
using namespace G4Sol;

//#define DEBUGRPZ

template<class Out>
TRPhiZTrackMDC<Out>::TRPhiZTrackMDC(const THyphiAttributes& attribut) : TDataProcessInterface<Out>("RPhiZTrackMDC"), att(attribut)
{

  OutputEvents = att.RPZ_OutputEvents;
  std::string temp_name_out       = att.Config.Get<std::string>("Output_Namefile");
  std::string temp_file_base_name = temp_name_out.substr(0, temp_name_out.find_last_of('.'));

  temp_file_base_name += "PhiZTest.root";
  namefilePhiZ = temp_file_base_name;

  RZfit = att.RPZ_RZfit;
  MDCWireType = att.RPZ_MDCWireType;

  if(OutputEvents)
    {
      att._logger->info("RF: OutputEvent set on !");
      f_phiZ = new TFile(namefilePhiZ, "RECREATE");
      t_phiZ = new TTree("PhiZTree", "Finder output tracks");

      t_phiZ->SetDirectory(f_phiZ);
      // t_phiZ->AutoSave();

      mg_trackPhiZ = new TMultiGraph;
      mg_trackPhiZ->SetNameTitle("PhiZTrack","PhiZTrack");
      mg_trackPhiZNoCorr = new TMultiGraph;
      mg_trackPhiZNoCorr->SetNameTitle("PhiZTrackNoCorr","PhiZTrackNoCorr");
      mg_trackRZ = new TMultiGraph;
      mg_trackRZ->SetNameTitle("RZTrack","RZTrack");
      tempEvent = 0;
      t_phiZ->Branch("PhiZTrack", "TMultiGraph", &mg_trackPhiZ, 12800, 0);
      t_phiZ->Branch("PhiZTrackNoCorr", "TMultiGraph", &mg_trackPhiZNoCorr, 12800, 0);
      t_phiZ->Branch("RZTrack", "TMultiGraph", &mg_trackRZ, 12800, 0);
      t_phiZ->Branch("nEvent",&tempEvent);
    }
  else
    {
      mg_trackPhiZ = new TMultiGraph;
      mg_trackPhiZ->SetNameTitle("PhiZTrack","PhiZTrack");
      mg_trackPhiZNoCorr = new TMultiGraph;
      mg_trackPhiZNoCorr->SetNameTitle("PhiZTrackNoCorr","PhiZTrackNoCorr");
      mg_trackRZ = new TMultiGraph;
      mg_trackRZ->SetNameTitle("RZTrack","RZTrack");
    }
}

template<class Out>
TRPhiZTrackMDC<Out>::~TRPhiZTrackMDC()
{
  if(OutputEvents)
    {
      f_phiZ->cd();
      t_phiZ->Write();
      // f_phiZ->Write();
      f_phiZ->Close();

      if(f_phiZ != nullptr)
	{
	  f_phiZ->Delete();
	  f_phiZ = nullptr;
	}
    }
}

template<class Out>
void TRPhiZTrackMDC<Out>::InitMT() { att._logger->error("E> Not supposed to be multithreaded !"); }

template<class Out>
ReturnRes::InfoM TRPhiZTrackMDC<Out>::operator()(FullRecoEvent& RecoEvent, Out* OutTree)
{

  int result_phiZ = Exec(RecoEvent, OutTree);

  return SoftExit(result_phiZ);
}

template<class Out>
int TRPhiZTrackMDC<Out>::Exec(FullRecoEvent& RecoEvent, Out* )
{
  return CheckTrackFinding(RecoEvent);
}

template<class Out>
ReturnRes::InfoM TRPhiZTrackMDC<Out>::SoftExit(int ) { return ReturnRes::Fine; }

template<class Out>
void TRPhiZTrackMDC<Out>::SelectHists()
{

  LocalHisto.h_RPhiZMDC_Chi2 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_RPhiZMDC_Chi2);
  LocalHisto.h_RPhiZMDC_Status = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_RPhiZMDC_Status);
  LocalHisto.h_ResidualMDC_dZ1 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualMDC_dZ1);
  LocalHisto.h_ResidualMDC_dZ2 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualMDC_dZ2);
  LocalHisto.h_RPhiZMDC_Sigma = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_RPhiZMDC_Sigma);
  LocalHisto.h_RPhiZMDC_Sigma2 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_RPhiZMDC_Sigma2);

  LocalHisto.h_ResidualMDC_dZ_PSB = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualMDC_dZ_PSB);
  LocalHisto.h_ResidualMDC_dZ_PSBE = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualMDC_dZ_PSBE);
  LocalHisto.h_ResidualMDC_dZ_PSFE = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualMDC_dZ_PSFE);
  LocalHisto.h_ResidualMDC_dZ_More6 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualMDC_dZ_More6);

  LocalHisto.h_PullMDC_dZ1 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_PullMDC_dZ1);
  LocalHisto.h_PullMDC_dZ2 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_PullMDC_dZ2);

  LocalHisto.h_PullMDC_dZ_PSB = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_PullMDC_dZ_PSB);
  LocalHisto.h_PullMDC_dZ_PSBE = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_PullMDC_dZ_PSBE);
  LocalHisto.h_PullMDC_dZ_PSFE = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_PullMDC_dZ_PSFE);
  LocalHisto.h_PullMDC_dZ_More6 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_PullMDC_dZ_More6);

}

template<class Out>
int TRPhiZTrackMDC<Out>::CheckTrackFinding(FullRecoEvent& RecoEvent)
{

  int ntrack = -1;
  for(auto it_trackInfo : RecoEvent.TrackInfo)
    {
      ntrack++;
      const int id_track = it_trackInfo.first;
      auto it_ListHits   = RecoEvent.TrackDAF.find(id_track);
      auto it_ListHitsSim = RecoEvent.TrackDAFSim.find(id_track);

      std::vector<std::tuple<int,int> > id_detMFT;
      std::vector<std::tuple<int,int> > id_detMG;
      std::vector<std::tuple<int,int> > id_detPSB;
      std::vector<std::tuple<int,int> > id_detPSFE;
      std::vector<std::tuple<int,int> > id_detPSBE;
      int index_mg = 0;

      att._logger->debug(" track {} TreeEvent {}",ntrack,tempEvent);
      LocalHisto.h_RPhiZMDC_Status->Fill("Status","AllTrack",1.);
      double meanPhi = 0;
      TComplex sumAngles(0,0);

      for(size_t id_det = 0; id_det < it_ListHits->second.size(); ++id_det)
        {
          int id_hit = it_ListHits->second[id_det];
          if(id_hit < 0)
            continue;

	  if(!((id_det>=G4Sol::MG01 && id_det<=G4Sol::MG17) || id_det==G4Sol::PSCE ||id_det==G4Sol::PSBE||id_det==G4Sol::PSFE))
	    continue;

	  auto TempSeg = RecoEvent.SegmentHit1Ds[id_det][id_hit];

	  TVector3 afterV;
	  for(size_t iS=0;iS<TempSeg.size();++iS)
	    {
	      double tempAbs = TMath::Hypot(TempSeg[iS][0],TempSeg[iS][1]);
	      TComplex tempAngles(TempSeg[iS][0]/tempAbs,TempSeg[iS][1]/tempAbs);
	      sumAngles += tempAngles;
	    }
	}

      meanPhi = TMath::ATan2(sumAngles.Im(), sumAngles.Re());

      for(size_t id_det = 0; id_det < it_ListHits->second.size(); ++id_det)
        {
          int id_hit = it_ListHits->second[id_det];
          if(id_hit < 0)
            continue;

	  if((id_det >= G4Sol::MiniFiberD1_x &&
	      id_det <= G4Sol::MiniFiberD2_u)) // || id_det == G4Sol::PSFE || id_det == G4Sol::PSBE))
	    {
	      id_detMFT.push_back({id_det,id_hit});
	    }

	  if(!((id_det>=G4Sol::MG01 && id_det<=G4Sol::MG17) || id_det==G4Sol::PSCE ||id_det==G4Sol::PSBE||id_det==G4Sol::PSFE))
	    continue;

	  if(id_det>=G4Sol::MG01 && id_det<=G4Sol::MG17)
	    id_detMG.push_back({index_mg,id_det});
	  if(id_det==G4Sol::PSCE)
	    id_detPSB.push_back({index_mg,id_det});
	  if(id_det==G4Sol::PSBE)
	    id_detPSBE.push_back({index_mg,id_det});
	  if(id_det==G4Sol::PSFE)
	    id_detPSFE.push_back({index_mg,id_det});

	  att._logger->debug("layer found : {}, {}, size: {}",id_det,G4Sol::nameLiteralDet.begin()[id_det], RecoEvent.SegmentHit1Ds[id_det][id_hit].size());
	  auto TempSeg = RecoEvent.SegmentHit1Ds[id_det][id_hit];

	  TVector3 afterV;

	  auto graphDetPhiZ = std::make_unique<TGraphErrors>();
	  auto graphDetPhiZ_NoCorr = std::make_unique<TGraphErrors>();
	  auto graphDetRZ = std::make_unique<TGraphErrors>();

	  double previousPhi = TMath::ATan2(TempSeg[0][1],TempSeg[0][0]);

	  for(size_t iS=0;iS<TempSeg.size();++iS)
	    {
	      afterV.SetXYZ(TempSeg[iS][0],TempSeg[iS][1],TempSeg[iS][2]);
	      double tempPhi = afterV.Phi();
	      //att._logger->debug("TempPhi : {}, previousPhi {} : Diff {}",tempPhi,previousPhi,tempPhi-previousPhi);
	      if(TMath::Abs(tempPhi - previousPhi) > 1.5 ) // ~ pi/2
		tempPhi += tempPhi<previousPhi ? 2*TMath::Pi() : -2*TMath::Pi() ;

	      //att._logger->debug("After Phi corr : tempPhi {} meanPhi {}" , tempPhi,meanPhi);
	      if(TMath::Abs(tempPhi - meanPhi)> 2.5) // ~ 3pi/2
		tempPhi += tempPhi<meanPhi ? 2*TMath::Pi() : -2*TMath::Pi() ;

	      //att._logger->debug("After Phi corr2 : tempPhi {} meanPhi {}" , tempPhi,meanPhi);
	      double tempPhi2 = tempPhi-meanPhi;

	      previousPhi = tempPhi;

	      graphDetPhiZ->SetPoint(iS,TempSeg[iS][2],tempPhi2);
	      graphDetPhiZ->SetPointError(iS,TempSeg[iS][5],TempSeg[iS][4]);

	      graphDetPhiZ_NoCorr->SetPoint(iS,TempSeg[iS][2],tempPhi2+meanPhi);
	      graphDetPhiZ_NoCorr->SetPointError(iS,TempSeg[iS][5],TempSeg[iS][4]);

	      graphDetRZ->SetPoint(iS,TempSeg[iS][2],afterV.Perp());
	      graphDetRZ->SetPointError(iS,TempSeg[iS][5],TempSeg[iS][3]);
	    }

	  graphDetPhiZ->SetMarkerStyle(27+id_det-G4Sol::MG01);
	  graphDetPhiZ->SetMarkerColor(30+id_det-G4Sol::MG01);
	  graphDetPhiZ->SetLineColor(30+id_det-G4Sol::MG01);
	  TString nameTemp1(G4Sol::nameLiteralDet.begin()[id_det]);
	  nameTemp1+="_PhiZ_";
	  nameTemp1+=id_hit;
	  graphDetPhiZ->SetNameTitle(nameTemp1,nameTemp1);

	  graphDetPhiZ_NoCorr->SetMarkerStyle(27+id_det-G4Sol::MG01);
	  graphDetPhiZ_NoCorr->SetMarkerColor(30+id_det-G4Sol::MG01);
	  graphDetPhiZ_NoCorr->SetLineColor(30+id_det-G4Sol::MG01);
	  TString nameTemp1_1(G4Sol::nameLiteralDet.begin()[id_det]);
	  nameTemp1_1+="_PhiZNoCorr_";
	  nameTemp1_1+=id_hit;
	  graphDetPhiZ_NoCorr->SetNameTitle(nameTemp1_1,nameTemp1_1);

	  graphDetRZ->SetMarkerStyle(27+id_det-G4Sol::MG01);
	  graphDetRZ->SetMarkerColor(30+id_det-G4Sol::MG01);
	  graphDetRZ->SetLineColor(30+id_det-G4Sol::MG01);
	  TString nameTemp2(G4Sol::nameLiteralDet.begin()[id_det]);
	  nameTemp2+="_RZ_";
	  nameTemp2+=id_hit;
	  graphDetRZ->SetNameTitle(nameTemp2,nameTemp2);

	  //att._logger->debug("nameDet : {}",nameTemp.Data());
	  //graphDetPhiZ->Print();

	  mg_trackPhiZ->Add(graphDetPhiZ.release(),"PL");
	  mg_trackPhiZNoCorr->Add(graphDetPhiZ_NoCorr.release(),"PL");
	  mg_trackRZ->Add(graphDetRZ.release(),"PL");

	  ++index_mg;
	}

      unsigned int Nb_hit = id_detMG.size() + id_detPSB.size();
      std::string tempNbhit = std::to_string(Nb_hit);
      if(id_detPSB.size()>0)
	LocalHisto.h_RPhiZMDC_Status->Fill("AllMDC+PSB",tempNbhit.c_str(),1.);
      else
	LocalHisto.h_RPhiZMDC_Status->Fill("AllMDConly",tempNbhit.c_str(),1.);

      if( Nb_hit < 3)
	{
	  auto* listG1 = mg_trackPhiZ->GetListOfGraphs();
	  listG1->Delete();
	  auto* listG1_1 = mg_trackPhiZNoCorr->GetListOfGraphs();
	  listG1_1->Delete();
	  auto* listG2 = mg_trackRZ->GetListOfGraphs();
	  listG2->Delete();

	  continue;
	}

      LocalHisto.h_RPhiZMDC_Status->Fill("Status","GoodTrack",1.);
      // for(int id_det : id_detMG)
      // 	att._logger->debug("list MDC : {}",id_det);

      auto* listPZ = mg_trackPhiZ->GetListOfGraphs();
      auto* listRZ = mg_trackRZ->GetListOfGraphs();

      auto graphInterPhiZ = std::make_unique<TGraphErrors>();
      auto graphInterRZ = std::make_unique<TGraphErrors>();

      graphInterPhiZ->SetNameTitle("IntersectionPhiZ","IntersectionPhiZ");
      graphInterRZ->SetNameTitle("IntersectionRZ","IntersectionRZ");

      std::map<int,std::vector<std::tuple<double,double,double,double> > > MG_posZestimated;
      // if(id_detPSB.size()>0)
      // 	std::cout<<"!> PSCE : "<<id_detPSB.size()<<"x: "<<dynamic_cast<TGraph*>(listRZ->At(std::get<0>(id_detPSB[0])))->GetMean(1)<<"y: "<<dynamic_cast<TGraph*>(listRZ->At(std::get<0>(id_detPSB[0])))->GetMean(2)<<"xe: "<<dynamic_cast<TGraph*>(listRZ->At(std::get<0>(id_detPSB[0])))->GetRMS(1)<<"ye: "<<dynamic_cast<TGraphErrors*>(listRZ->At(std::get<0>(id_detPSB[0])))->GetErrorY(0)<<"\n";

      for(size_t iM = 1; iM<id_detMG.size();++iM)
	{
	  int id_mgN1 = std::get<0>(id_detMG[iM-1]);
	  int id_mgN2 = std::get<0>(id_detMG[iM]);
	  int id_detN1 = std::get<1>(id_detMG[iM-1]);
	  int id_detN2 = std::get<1>(id_detMG[iM]);

	  att._logger->debug("check intersection : {}, {} | {} {} | {}",id_mgN1, id_mgN2, G4Sol::nameLiteralDet.begin()[id_detN1], G4Sol::nameLiteralDet.begin()[id_detN2], id_detN2 - id_detN1 );

	  if(id_detN2 - id_detN1 != 1 && id_detN2 - id_detN1 !=2)
	    continue;

	  auto [ret, z_intersect, phi_intersect, r_intersect1, r_intersect2, z_e, phi_e, r1_e, r2_e] = MDC_PRZ::f_interset(listPZ, listRZ, id_mgN1, id_mgN2);

	  if(ret!=1)
	    continue;

	  graphInterPhiZ->SetPoint(graphInterPhiZ->GetN(),z_intersect,phi_intersect);
	  graphInterPhiZ->SetPointError(graphInterPhiZ->GetN()-1,z_e,phi_e);

	  graphInterRZ->SetPoint(graphInterRZ->GetN(),z_intersect, r_intersect1);
	  graphInterRZ->SetPointError(graphInterRZ->GetN()-1,z_e, r1_e);
	  graphInterRZ->SetPoint(graphInterRZ->GetN(),z_intersect, r_intersect2);
	  graphInterRZ->SetPointError(graphInterRZ->GetN()-1,z_e, r2_e);

	  att._logger->debug("   -> Interset done : z_inter {} phi_intersect {} | r1_intersect {} r2_interset {}",z_intersect, phi_intersect,r_intersect1,r_intersect2);

	  int lowId, highId;
	  if(r_intersect1<r_intersect2)
	    {
	      lowId = id_detN1;
	      highId = id_detN2;
	    }
	  else
	    {
	      lowId = id_detN2;
	      highId = id_detN1;
	    }
	  auto it_N1 = MG_posZestimated.find(lowId);
	  if(it_N1!=MG_posZestimated.end())
	    it_N1->second.push_back({z_intersect,r_intersect1,z_e,r1_e});
	  else
	    MG_posZestimated.insert({lowId,std::vector<std::tuple<double,double,double,double> >{{z_intersect,r_intersect1,z_e,r1_e}}});

	  auto it_N2 = MG_posZestimated.find(highId);
	  if(it_N2!=MG_posZestimated.end())
	    it_N2->second.push_back({z_intersect,r_intersect2,z_e,r2_e});
	  else
	    MG_posZestimated.insert({highId,std::vector<std::tuple<double,double,double,double> >{{z_intersect,r_intersect2,z_e,r2_e}}});

	}


      for(size_t iM = 2; iM<id_detMG.size();++iM)
	{
	  int id_mgN1 = std::get<0>(id_detMG[iM-2]);
	  int id_mgN2 = std::get<0>(id_detMG[iM]);
	  int id_detN1 = std::get<1>(id_detMG[iM-2]);
	  int id_detN2 = std::get<1>(id_detMG[iM]);

	  att._logger->debug("check intersection : {}, {} | {} {} | {}",id_mgN1, id_mgN2, G4Sol::nameLiteralDet.begin()[id_detN1], G4Sol::nameLiteralDet.begin()[id_detN2], id_detN2 - id_detN1 );

	  if(id_detN2 - id_detN1 != 2)
	    continue;

	  auto [ret, z_intersect, phi_intersect, r_intersect1, r_intersect2,z_e,phi_e,r1_e,r2_e] = MDC_PRZ::f_interset(listPZ, listRZ, id_mgN1, id_mgN2);

	  if(ret!=1)
	    continue;

	  graphInterPhiZ->SetPoint(graphInterPhiZ->GetN(),z_intersect,phi_intersect);
	  graphInterPhiZ->SetPointError(graphInterPhiZ->GetN()-1,z_e,phi_e);
	  graphInterRZ->SetPoint(graphInterRZ->GetN(),z_intersect, r_intersect1);
	  graphInterRZ->SetPointError(graphInterRZ->GetN()-1,z_e, r1_e);
	  graphInterRZ->SetPoint(graphInterRZ->GetN(),z_intersect, r_intersect2);
	  graphInterRZ->SetPointError(graphInterRZ->GetN()-1,z_e, r2_e);

	  att._logger->debug("   -> Interset done : z_inter {} phi_intersect {} | r1_intersect {} r2_interset {}",z_intersect, phi_intersect,r_intersect1,r_intersect2);
	  int lowId, highId;
	  if(r_intersect1<r_intersect2)
	    {
	      lowId = id_detN1;
	      highId = id_detN2;
	    }
	  else
	    {
	      lowId = id_detN2;
	      highId = id_detN1;
	    }
	  auto it_N1 = MG_posZestimated.find(lowId);
	  if(it_N1!=MG_posZestimated.end())
	    it_N1->second.push_back({z_intersect,r_intersect1,z_e,r1_e});
	  else
	    MG_posZestimated.insert({lowId,std::vector<std::tuple<double,double,double,double> >{{z_intersect,r_intersect1,z_e,r1_e}}});

	  auto it_N2 = MG_posZestimated.find(highId);
	  if(it_N2!=MG_posZestimated.end())
	    it_N2->second.push_back({z_intersect,r_intersect2,z_e,r2_e});
	  else
	    MG_posZestimated.insert({highId,std::vector<std::tuple<double,double,double,double> >{{z_intersect,r_intersect2,z_e,r2_e}}});

	}


      if(id_detPSB.size()>0)
	{
	  std::vector<std::tuple<double,double,double,double> > posAllPSB;
	  int id_det = -1;
	  for(auto [id1,id2] : id_detPSB)
	    {
	      id_det =id2;
	      att._logger->debug("!> PSCE : {} x: {} y: {} xe: {} ye: {}",id1,dynamic_cast<TGraph*>(listRZ->At(std::get<0>(id_detPSB[0])))->GetMean(1),dynamic_cast<TGraph*>(listRZ->At(std::get<0>(id_detPSB[0])))->GetMean(2),dynamic_cast<TGraph*>(listRZ->At(std::get<0>(id_detPSB[0])))->GetRMS(1),dynamic_cast<TGraphErrors*>(listRZ->At(std::get<0>(id_detPSB[0])))->GetErrorY(0));
	      posAllPSB.push_back({dynamic_cast<TGraph*>(listRZ->At(id1))->GetMean(1), dynamic_cast<TGraph*>(listRZ->At(id1))->GetMean(2), dynamic_cast<TGraph*>(listRZ->At(id1))->GetRMS(1), dynamic_cast<TGraphErrors*>(listRZ->At(id1))->GetErrorY(0)});
	    }

	  std::tuple<double,double,double,double> posPSB = {0.,0.,0.,0.};
	  for(auto [temp_z,temp_r,temp_zerr,temp_rerr] : posAllPSB)
	    {
	      std::get<0>(posPSB) += temp_z;
	      std::get<1>(posPSB) += temp_r;
	      std::get<2>(posPSB) += temp_zerr;
	      std::get<3>(posPSB) += temp_rerr;
	    }
	  double posPSBsize = static_cast<double>(posAllPSB.size());
	  std::get<0>(posPSB) /= posPSBsize;
	  std::get<1>(posPSB) /= posPSBsize;
	  std::get<2>(posPSB) /= posPSBsize;
	  std::get<3>(posPSB) /= posPSBsize;

	  graphInterRZ->SetPoint(graphInterRZ->GetN(),std::get<0>(posPSB),std::get<1>(posPSB));
	  graphInterRZ->SetPointError(graphInterRZ->GetN()-1,std::get<2>(posPSB),std::get<3>(posPSB));

	  MG_posZestimated.insert({id_det,std::vector<std::tuple<double,double,double,double> >{{std::get<0>(posPSB),std::get<1>(posPSB),std::get<2>(posPSB),std::get<3>(posPSB)}}});

	}

      if(id_detMFT.size()>0)
	{
	  auto f_extractXYZ_Fiber = [](const auto& ListHits, int id_det, int id_hit) {
	    genfit::PlanarMeasurement* currentHit = dynamic_cast<genfit::PlanarMeasurement*>(ListHits[id_det][id_hit].get());
	    auto dummyState                       = genfit::StateOnPlane();
	    auto currentPlane                     = currentHit->constructPlane(genfit::StateOnPlane());
	    auto HitCoord                         = currentHit->getRawHitCoords();
	    auto HitCov                           = currentHit->getRawHitCov();

	    if(HitCoord.GetNrows() != 2)
	      return std::make_tuple(-1,0.,0.,0.,0.,0.);

	    //auto posPlane = currentPlane->getO();
	    //auto planeU = currentPlane->getU();
	    //auto planeV = currentPlane->getV();

	    //att._logger->debug("ImproveFiber : HitCoord : {} | O [{}, {}, {}] U [{}, {}, {}] V [{}, {}, {}]",HitCoord(0),posPlane.X(),posPlane.Y(),posPlane.Z(),planeU.X(),planeU.Y(),planeU.Z(),planeV.X(),planeV.Y(),planeV.Z());

	    //double res = std::sqrt(HitCov(0,0));

	    TVector3 posXYZ = currentPlane->toLab({HitCoord(0),HitCoord(1)});

	    return std::make_tuple(1,posXYZ.X(),posXYZ.Y(),posXYZ.Z(), std::sqrt(HitCov(0,0)),std::sqrt(HitCov(1,1)));

	  };
	  int id_detF1 = -1, id_detF2 = -1;
	  std::array<double,5> xyz_f1, xyz_f2;
	  for(auto [id_detF,id_hitF] : id_detMFT)
	    {
	      if((id_detF >= G4Sol::MiniFiberD1_x && id_detF <= G4Sol::MiniFiberD1_v))
		{
		  auto [res_f,x_f,y_f,z_f,xer_f,yer_f] = f_extractXYZ_Fiber(RecoEvent.ListHits,id_detF,id_hitF);
		  if(res_f == 1)
		    {
		      xyz_f1 = {x_f,y_f,z_f,xer_f*10.,yer_f*10.};
		      id_detF1 = id_detF;
		    }
		}
	      if((id_detF >= G4Sol::MiniFiberD2_x && id_detF <= G4Sol::MiniFiberD2_u))
		{
		  auto [res_f,x_f,y_f,z_f,xer_f,yer_f] = f_extractXYZ_Fiber(RecoEvent.ListHits,id_detF,id_hitF);
		  if(res_f == 1)
		    {
		      xyz_f2 = {x_f,y_f,z_f,xer_f*10.,yer_f*10.};
		      id_detF2 = id_detF;
		    }
		}
	    }

	  if(id_detF1!=-1)
	    {
	      MG_posZestimated.insert({id_detF1,std::vector<std::tuple<double,double,double,double> >{{xyz_f1[2], TMath::Hypot(xyz_f1[0],xyz_f1[1]), 0.1, std::sqrt(TMath::Sq(xyz_f1[0]*xyz_f1[3])/(TMath::Sq(xyz_f1[0])+TMath::Sq(xyz_f1[1]))+TMath::Sq(xyz_f1[1]*xyz_f1[4])/(TMath::Sq(xyz_f1[0])+TMath::Sq(xyz_f1[1])))}}});

	      graphInterRZ->SetPoint(graphInterRZ->GetN(),xyz_f1[2], TMath::Hypot(xyz_f1[0],xyz_f1[1]));
	      graphInterRZ->SetPointError(graphInterRZ->GetN()-1,.5, std::sqrt(TMath::Sq(xyz_f1[0]*xyz_f1[3])/(TMath::Sq(xyz_f1[0])+TMath::Sq(xyz_f1[1]))+TMath::Sq(xyz_f1[1]*xyz_f1[4])/(TMath::Sq(xyz_f1[0])+TMath::Sq(xyz_f1[1]))));

	    }
	  if(id_detF2!=-1)
	    {
	      MG_posZestimated.insert({id_detF2,std::vector<std::tuple<double,double,double,double> >{{xyz_f2[2], TMath::Hypot(xyz_f2[0],xyz_f2[1]), 0.1, std::sqrt(TMath::Sq(xyz_f2[0]*xyz_f2[3])/(TMath::Sq(xyz_f2[0])+TMath::Sq(xyz_f2[1]))+TMath::Sq(xyz_f2[1]*xyz_f2[4])/(TMath::Sq(xyz_f2[0])+TMath::Sq(xyz_f2[1])))}}});
	      graphInterRZ->SetPoint(graphInterRZ->GetN(),xyz_f2[2], TMath::Hypot(xyz_f2[0],xyz_f2[1]));
	      graphInterRZ->SetPointError(graphInterRZ->GetN()-1,.5, std::sqrt(TMath::Sq(xyz_f2[0]*xyz_f2[3])/(TMath::Sq(xyz_f2[0])+TMath::Sq(xyz_f2[1]))+TMath::Sq(xyz_f2[1]*xyz_f2[4])/(TMath::Sq(xyz_f2[0])+TMath::Sq(xyz_f2[1]))));
	    }
	}



      auto graphResultZ = std::make_unique<TGraphErrors>();
      graphResultZ->SetNameTitle("IntersectionResultZ","IntersectionResultZ");

      for(const auto& it_MG : MG_posZestimated)
	{
	  att._logger->debug("E> CheckIntersectionResultZ : wrong r_intersect ! {} | size {}", G4Sol::nameLiteralDet.begin()[it_MG.first], MG_posZestimated.size());
	  //if(MG_posZestimated.size() == 1)
	  //  continue;

	  double z_mean = 0.;
	  double r_check = 0;
	  double dem_Z = 0;
	  double dem_R = 0;
	  double sig_Z = 0;
	  double sig_R = 0;
	  for(auto [tempZ, tempR, tempZ_e, tempR_e] : it_MG.second)
	    {
	      z_mean += tempZ/TMath::Sq(tempZ_e);
	      r_check += tempR/TMath::Sq(tempR_e);
	      dem_Z += 1./TMath::Sq(tempZ_e);
	      dem_R += 1./TMath::Sq(tempR_e);
	      sig_Z += TMath::Sq(tempZ_e);
	      sig_R += TMath::Sq(tempR_e);
	    }
	  z_mean /= dem_Z;//static_cast<double>(it_MG.second.size());
	  r_check /= dem_R;//static_cast<double>(it_MG.second.size());

	  sig_Z /= static_cast<double>(it_MG.second.size());
	  sig_R /= static_cast<double>(it_MG.second.size());

	  // if(TMath::Abs(r_check-std::get<1>(it_MG.second[0]))>0.5)
	  //   {
	  for(auto [tempZ, tempR, tempZ_e, tempR_e] : it_MG.second)
	    att._logger->debug("E> z {} r {} | +- z {} +- r {}",tempZ, tempR, tempZ_e, tempR_e);
	  // }

	  att._logger->debug("mean> z {} r {} | +- z {} +- r {}",z_mean, r_check, TMath::Sqrt(sig_Z), TMath::Sqrt(sig_R));
	  graphResultZ->SetPoint(graphResultZ->GetN(), z_mean,r_check);
	  graphResultZ->SetPointError(graphResultZ->GetN()-1, TMath::Sqrt(sig_Z), TMath::Sqrt(sig_R));
	}

      auto graphResultZ_Clean = std::make_unique<TGraphErrors>();
      graphResultZ_Clean->SetNameTitle("IntersectionResultZ_Clean","IntersectionResultZ_Clean");

      std::set<int> indexToExclude_graphR;

      att._logger->debug(" IntersectionResultZ_Clean {}",graphResultZ->GetN());
      // int id_middle = graphResultZ->GetN()/2;

      // att._logger->debug(" from middle to end");
      // for(int i=id_middle; i < graphResultZ->GetN()-1;++i)
      // 	{
      // 	  double tempZ1 = 0,tempR1=0;
      // 	  graphResultZ->GetPoint(i,tempZ1,tempR1);
      // 	  double errZ1=graphResultZ->GetErrorX(i);
      // 	  double errR1=graphResultZ->GetErrorY(i);
      // 	  double tempZ0 = 0,tempR0=0, errZ0;
      // 	  graphResultZ->GetPoint(i-1,tempZ0,tempR0);
      // 	  errZ0=graphResultZ->GetErrorX(i+1);
      // 	  double tempZ2 = 0,tempR2=0, errZ2;
      // 	  graphResultZ->GetPoint(i+1,tempZ2,tempR2);
      // 	  errZ2=graphResultZ->GetErrorX(i+1);

      // 	  std::array<double,2> ABvec = {tempZ1-tempZ0,tempR1-tempR0};
      // 	  std::array<double,2> BCvec = {tempZ2-tempZ1,tempZ2-tempR1};

      // 	  double AngleABC = (ABvec[0]*BCvec[0]+ABvec[1]*BCvec[1])/std::hypot(ABvec[0],ABvec[1])/std::hypot(BCvec[0],BCvec[1]);

      // 	  att._logger->debug(" #{} : Z0 {} R0 {} errZ0 {} | Z1 {} R1 {} errZ1 {} errR1 {} | Z2 {} R2 {} errZ2 {} | diff : {} <> {} / {} <> {} | cos(angleABC) {}",i,tempZ0, tempR0,errZ0,tempZ1,tempR1,errZ1,errR1,tempZ2, tempR2,errZ2, tempZ1+errZ1/2,tempZ2-errZ2/2, tempZ0+errZ0/2,tempZ1-errZ1/2, AngleABC);

      // 	  if(tempR2 < tempR1 || tempR1 < tempR0)
      // 	    att._logger->error("E> RPhiZ : graphResultZ R[{}] < R[{}] ! R[{}] = {} / R[{}] = {}",i+1,i,i, tempR1, i+1, tempR2);

      // 	  if((tempZ1 > tempZ2 || tempZ0 > tempZ1) && TMath::Abs(AngleABC)<0.1 )
      // 	    {
      // 	      //graphResultZ_Clean->SetPoint(i,tempZ1,tempR1);
      // 	      //graphResultZ_Clean->SetPointError(graphResultZ_Clean->GetN()-1,errZ1,errR1);
      // 	      att._logger->debug(" -> Exclude insert {}",i);
      // 	      indexToExclude_graphR.insert(i);
      // 	    }
      // 	}

      // att._logger->debug(" from middle to first");
      // for(int i=1; i < id_middle; ++i)
      // 	{
      // 	  double tempZ1 = 0,tempR1=0;
      // 	  graphResultZ->GetPoint(i,tempZ1,tempR1);
      // 	  double errZ1=graphResultZ->GetErrorX(i);
      // 	  double errR1=graphResultZ->GetErrorY(i);
      // 	  double tempZ0 = 0,tempR0=0, errZ0;
      // 	  graphResultZ->GetPoint(i-1,tempZ0,tempR0);
      // 	  errZ0=graphResultZ->GetErrorX(i+1);
      // 	  double tempZ2 = 0,tempR2=0, errZ2;
      // 	  graphResultZ->GetPoint(i+1,tempZ2,tempR2);
      // 	  errZ2=graphResultZ->GetErrorX(i+1);

      // 	  std::array<double,2> ABvec = {tempZ1-tempZ0,tempR1-tempR0};
      // 	  std::array<double,2> BCvec = {tempZ2-tempZ1,tempZ2-tempR1};

      // 	  double AngleABC = (ABvec[0]*BCvec[0]+ABvec[1]*BCvec[1])/std::hypot(ABvec[0],ABvec[1])/std::hypot(BCvec[0],BCvec[1]);

      // 	  att._logger->debug(" #{} : Z0 {} R0 {} errZ0 {} | Z1 {} R1 {} errZ1 {} errR1 {} | Z2 {} R2 {} errZ2 {} | diff : {} <> {} / {} <> {} | cos(angleABC) {}",i,tempZ0, tempR0,errZ0,tempZ1,tempR1,errZ1,errR1,tempZ2, tempR2,errZ2, tempZ1+errZ1/2,tempZ2-errZ2/2, tempZ0+errZ0/2,tempZ1-errZ1/2, AngleABC);

      // 	  if(tempR2 < tempR1 || tempR1 < tempR0)
      // 	    att._logger->error("E> RPhiZ : graphResultZ R[{}] < R[{}] ! R[{}] = {} / R[{}] = {}",i+1,i,i, tempR1, i+1, tempR2);

      // 	  if((tempZ1 > tempZ2 || tempZ0 > tempZ1) && TMath::Abs(AngleABC)<0.1 )
      // 	    {
      // 	      //graphResultZ_Clean->SetPoint(i,tempZ1,tempR1);
      // 	      //graphResultZ_Clean->SetPointError(graphResultZ_Clean->GetN()-1,errZ1,errR1);
      // 	      att._logger->debug(" -> Exclude insert {}",i);
      // 	      indexToExclude_graphR.insert(i);
      // 	    }
      // 	}

      if(graphResultZ->GetN()>0)
	{
	  double CorrelationFactorGraphRes = graphResultZ->GetCorrelationFactor();
	  double tempM = CorrelationFactorGraphRes * std::sqrt(1. - CorrelationFactorGraphRes*CorrelationFactorGraphRes);
	  double tempB = graphResultZ->GetMean(2) - tempM*graphResultZ->GetMean(1); // R = m*z+b -> z = R/m - b/m

	  att._logger->debug(" Correlation line : M={} B={}",tempM, tempB);

	  if(TMath::Abs(tempM)>1e-6)
	    {
	      std::vector<double> diff_XRes;
	      for(int i_resZ=0;i_resZ<graphResultZ->GetN();++i_resZ)
		{
		  double tempZ1 = 0,tempR1=0;
		  graphResultZ->GetPoint(i_resZ,tempZ1,tempR1);
		  double errZ1=graphResultZ->GetErrorX(i_resZ);
		  //double errR1=graphResultZ->GetErrorY(i_resZ);

		  double tempZ_fromLine1 = tempR1/tempM - tempB/tempM;

		  double diffZ1 = TMath::Abs((tempZ1+errZ1) - tempZ_fromLine1);
		  double diffZ2 = TMath::Abs((tempZ1-errZ1) - tempZ_fromLine1);

		  double diffZ = diffZ1 < diffZ2 ? diffZ1 : diffZ2;

		  att._logger->debug("  Diff Z #{}: {} | {} ",i_resZ,diffZ1, diffZ2);
		  diff_XRes.push_back(diffZ);
		}

	      std::sort(diff_XRes.begin(), diff_XRes.end());
	      std::array<double,3> TempQuantile;
	      std::array<double,3> TempProbQ = {0.25,0.5,0.75};
	      TMath::Quantiles(diff_XRes.size(),3,diff_XRes.data(), TempQuantile.data(), TempProbQ.data());

	      att._logger->debug("Quantiles : Median {}, Q1 {}, Q3 {}",TempQuantile[1],TempQuantile[0],TempQuantile[2]);

	      double tempIQR = TempQuantile[2]-TempQuantile[0];
	      double lowFence = TempQuantile[0] - 1.5*tempIQR;
	      double highFence = TempQuantile[2] + 1.5*tempIQR;

	      att._logger->debug("Fences : {} low {} high {}",tempIQR,lowFence,highFence);
	      for(int i_resZ=0;i_resZ<graphResultZ->GetN();++i_resZ)
		{
		  double tempZ1 = 0,tempR1=0;
		  graphResultZ->GetPoint(i_resZ,tempZ1,tempR1);
		  double tempZ_fromLine1 = tempR1/tempM - tempB/tempM;
		  double errZ1=graphResultZ->GetErrorX(i_resZ);
		  double diffZ1 = TMath::Abs((tempZ1+errZ1) - tempZ_fromLine1);
		  double diffZ2 = TMath::Abs((tempZ1-errZ1) - tempZ_fromLine1);

		  double diffZ0 = diffZ1 < diffZ2 ? diffZ1 : diffZ2;
		  double diffZ = tempZ1-tempZ_fromLine1;

		  //double errZ1=graphResultZ->GetErrorX(i_resZ);
		  att._logger->debug(" Check Outliers {}: diffZ0 {} diffZ {} <> {} | {}",i_resZ,diffZ0,diffZ,lowFence, highFence);

		  if(diffZ0 < lowFence)
		    indexToExclude_graphR.insert(i_resZ);
		  if(diffZ0 > highFence)
		    indexToExclude_graphR.insert(i_resZ);
		}

	      att._logger->debug("index size {}",indexToExclude_graphR.size());
	      for(auto i : indexToExclude_graphR)
		att._logger->debug("index : {}",i);

	      if(indexToExclude_graphR.size()!=0)
		{
		  for(int i=0;i<graphResultZ->GetN();++i)
		    {
		      auto it_i = indexToExclude_graphR.find(i);
		      double tempZ1 = 0,tempR1=0;
		      graphResultZ->GetPoint(i,tempZ1,tempR1);

		      if(it_i == indexToExclude_graphR.end() || tempZ1<245. || (i == graphResultZ->GetN()-1 && id_detPSB.size()>0))
			{
			  double errZ1=graphResultZ->GetErrorX(i);
			  double errR1=graphResultZ->GetErrorY(i);
			  graphResultZ_Clean->SetPoint(graphResultZ_Clean->GetN(),tempZ1,tempR1);
			  graphResultZ_Clean->SetPointError(graphResultZ_Clean->GetN()-1,errZ1,errR1);
			}
		    }
		}
	    }
	}
      att._logger->debug("Result cleaning :");
      //graphResultZ_Clean->Print();

      std::vector<std::tuple<double,double,double,double,double,double> > graphFits(17);
      std::tuple<int,double,double> paraFitRZ;

      if(graphInterRZ->GetN() >=2)
	LocalHisto.h_RPhiZMDC_Status->Fill("Status","GoodTrack_InterRZ>=2",1.);
      else
	LocalHisto.h_RPhiZMDC_Status->Fill("Status","GoodTrack_InterRZ<2",1.);

      std::string Nb_InterRZstr = std::to_string(graphInterRZ->GetN());
      LocalHisto.h_RPhiZMDC_Status->Fill("InterRZ_GetN",Nb_InterRZstr.c_str(),1.);

      if(graphResultZ_Clean->GetN()>=2)
	LocalHisto.h_RPhiZMDC_Status->Fill("Status","GoodTrack_InterRZClean>=2",1.);
      else
	LocalHisto.h_RPhiZMDC_Status->Fill("Status","GoodTrack_InterRZClean<2",1.);

      std::string Nb_ResultRZstr = std::to_string(graphResultZ->GetN());
      LocalHisto.h_RPhiZMDC_Status->Fill("ResultRZ_GetN",Nb_ResultRZstr.c_str(),1.);
      std::string Nb_ResultRZCleanstr = std::to_string(graphResultZ_Clean->GetN());
      LocalHisto.h_RPhiZMDC_Status->Fill("ResultRZClean_GetN",Nb_ResultRZCleanstr.c_str(),1.);

      if(RZfit && graphInterRZ->GetN() >=2)
	{
	  att._logger->debug("RZfit : IntersectRZ:{} | mg_trackRZ: {}", graphInterRZ->GetN(), mg_trackRZ->GetListOfGraphs()->GetEntries());

	  double CorrelationFactorGraphRes = graphResultZ->GetCorrelationFactor();
	  double tempM = CorrelationFactorGraphRes * std::sqrt(1. - CorrelationFactorGraphRes*CorrelationFactorGraphRes);
	  double tempB = graphResultZ->GetMean(2) - tempM*graphResultZ->GetMean(1); // R = m*z+b -> z = R/m - b/m
	  pol.SetParameter(0,tempB);
	  pol.SetParameter(1,tempM);

	  double chi2_sC = 0.;
	  for(int iRZ=0;iRZ<graphResultZ->GetN();++iRZ)
	    {
	      double tempZ1 = 0,tempR1=0;
	      graphResultZ->GetPoint(iRZ,tempZ1,tempR1);
	      double errZ1=graphResultZ->GetErrorX(iRZ);
	      double errR1=graphResultZ->GetErrorY(iRZ);

	      chi2_sC += TMath::Sq(tempR1 - (tempM*tempZ1+tempB))/(TMath::Sq(errR1) + TMath::Sq(errZ1)*tempM);
	    }
	  att._logger->debug("chi2 from correlation : {}", chi2_sC);

	  TFitResultPtr fitR = graphInterRZ->Fit(&pol,"SNBQ");

	  TFitResult* Rf = fitR.Get();
	  double p0 = Rf->Parameter(0);
	  double p1 = Rf->Parameter(1);
	  TMatrixD cov = Rf->GetCovarianceMatrix();
	  double chi2  = Rf->Chi2();
	  double ndf = Rf->Ndf();
	  int res_fit = 1;
	  if(!cov.IsSymmetric())
	    {
	      att._logger->error("E> fit did not work ! RPhiZ > fitRZ : chi2 {} ndf {} prob {}| p0 {} p1 {} | cov {} {} {} {}",chi2,ndf,Rf->Prob(), p0, p1, cov(0,0),cov(0,1),cov(1,0),cov(1,1));
	      graphInterRZ->Print();
	    }
	  LocalHisto.h_RPhiZMDC_Chi2->Fill("Chi2_1",chi2,1.);
	  LocalHisto.h_RPhiZMDC_Chi2->Fill("Chi2r_1",chi2/ndf,1.);
	  att._logger->debug("RPhiZ > fitRZ : chi2 {} ndf {} prob {}| p0 {} p1 {}",chi2,ndf,Rf->Prob(), p0, p1);
	  if(graphResultZ_Clean->GetN()>=2)
	    {
	      TFitResultPtr fitR_C = graphResultZ_Clean->Fit(&pol,"SNBQ");
	      TFitResult* Rf_C = fitR_C.Get();
	      double p0_C = Rf_C->Parameter(0);
	      double p1_C = Rf_C->Parameter(1);
	      //paraFitRZ = {1,p0,p1};
	      TMatrixD cov_C = Rf_C->GetCovarianceMatrix();
	      double chi2_C  = Rf_C->Chi2();
	      double ndf_C = Rf_C->Ndf();
	      LocalHisto.h_RPhiZMDC_Chi2->Fill("Chi2_2",chi2_C,1.);
	      LocalHisto.h_RPhiZMDC_Chi2->Fill("Chi2r_2",chi2_C/ndf_C,1.);
	      att._logger->debug("RPhiZ > fitRZ 2 : {} | {} = p0 {} p1 {}",chi2_C, ndf_C,Rf_C->Prob(), p0_C, p1_C);
	      if(chi2_C/ndf_C < chi2/ndf)
		{
		  att._logger->debug("RPhiZ > fitRZ replaced with ResultZ_Clean !");
		  res_fit = 2;
		  p0 = p0_C;
		  p1 = p1_C;
		  cov = cov_C;
		  chi2 = chi2_C;
		  ndf = ndf_C;

		  LocalHisto.h_RPhiZMDC_Status->Fill("BetterFit","Clean",1.);
		}
	      else
		LocalHisto.h_RPhiZMDC_Status->Fill("BetterFit","Old",1.);

	      LocalHisto.h_RPhiZMDC_Status->Fill("BetterFit","All",1.);

	    }
	  // else if(chi2 > chi2_sC)
	  //   {
	  //     LocalHisto.h_RPhiZMDC_Status->Fill("BetterFit","Corr",1.);
	  //     LocalHisto.h_RPhiZMDC_Chi2->Fill("Chi2CR_1",chi2_sC,1.);
	  //     LocalHisto.h_RPhiZMDC_Chi2->Fill("Chi2CRr_1",chi2_sC/ndf,1.);
	  //     p0 = tempB;
	  //     p1 = tempM;
	  //     cov(0,0) = TMath::Sq(graphResultZ->GetRMS(1));
	  //     cov(1,0) = graphResultZ->GetCovariance();
	  //     cov(0,1) = graphResultZ->GetCovariance();
	  //     cov(1,1) = TMath::Sq(graphResultZ->GetRMS(2));
	  //     chi2 = chi2_sC;
	  //   }

	  paraFitRZ = {res_fit,p0,p1};


	  RecoEvent.paramFitRZ.insert({id_track,{p0,p1,cov}});

	  auto graphRZfit1 = std::make_unique<TGraphErrors>();
	  graphRZfit1->SetNameTitle("FitIntRZ","FitIntRZ");
	  auto graphRZfit2 = std::make_unique<TGraphErrors>();
	  graphRZfit2->SetNameTitle("FitErrRZ","FitErrRZ");
	  auto graphPZfit2 = std::make_unique<TGraphErrors>();
	  graphPZfit2->SetNameTitle("FitErrPZ","FitErrPZ");

	  for( auto [id_mg, id_det] : id_detMG)
	    {
	      auto [res, Zin, Rin, meanZ, meanR, errZ, errR] = MDC_PRZ::f_intersetTrack(listRZ,id_mg, p0,p1,cov);
	      graphRZfit1->SetPoint(graphRZfit1->GetN(),Zin,Rin);

	      graphRZfit2->SetPoint(graphRZfit2->GetN(),meanZ,meanR);
	      graphRZfit2->SetPointError(graphRZfit2->GetN()-1,errZ,errR);

	      double meanPhi = dynamic_cast<TGraphErrors*>(listPZ->At(id_mg))->Eval(meanZ);
	      double errPhi1 = dynamic_cast<TGraphErrors*>(listPZ->At(id_mg))->Eval(meanZ+errZ);
	      double errPhi2 = dynamic_cast<TGraphErrors*>(listPZ->At(id_mg))->Eval(meanZ-errZ);
	      double errPhi = 0.5*(TMath::Sq(meanPhi-errPhi1)+TMath::Sq(meanPhi-errPhi2));
	      errPhi = TMath::Sqrt(errPhi);

	      graphPZfit2->SetPoint(graphPZfit2->GetN(),meanZ,meanPhi);
	      graphPZfit2->SetPointError(graphPZfit2->GetN()-1,errZ,errPhi);

	      int id_det_reduced = id_det - G4Sol::MG01;
	      graphFits[id_det_reduced] = {meanZ,meanR,errZ,errR,meanPhi,errPhi};
	    }

	  graphRZfit1->SetMarkerStyle(23);
	  graphRZfit1->SetMarkerSize(1.2);
	  graphRZfit1->SetMarkerColor(6);

	  graphRZfit2->SetMarkerStyle(23);
	  graphRZfit2->SetMarkerSize(1.2);
	  graphRZfit2->SetMarkerColor(7);

	  graphPZfit2->SetMarkerStyle(23);
	  graphPZfit2->SetMarkerSize(1.2);
	  graphPZfit2->SetMarkerColor(7);

	  if(OutputEvents && graphRZfit1->GetN()!=0)
	    {
	      mg_trackRZ->Add(graphRZfit1.release(),"P");
	      mg_trackRZ->Add(graphRZfit2.release(),"P");
	      mg_trackPhiZ->Add(graphPZfit2.release(),"P");

	    }
	}
      else
	paraFitRZ = {-1,0,0};

      if(OutputEvents && graphResultZ->GetN()!=0)
	{
	  graphResultZ->SetMarkerStyle(22);
	  graphResultZ->SetMarkerSize(1.2);
	  graphResultZ->SetMarkerColor(4);

	  graphResultZ_Clean->SetMarkerStyle(25);
	  graphResultZ_Clean->SetMarkerSize(1.5);
	  graphResultZ_Clean->SetMarkerColor(8);

	  mg_trackRZ->Add(graphResultZ.release(),"P");
	  mg_trackRZ->Add(graphResultZ_Clean.release(),"P");
	}

      if(OutputEvents && graphInterPhiZ->GetN()!=0)
	{
	  graphInterPhiZ->SetMarkerStyle(20);
	  graphInterPhiZ->SetMarkerSize(1.2);
	  graphInterPhiZ->SetMarkerColor(3);
	  graphInterRZ->SetMarkerStyle(20);
	  graphInterRZ->SetMarkerSize(1.2);
	  graphInterRZ->SetMarkerColor(3);

	  mg_trackPhiZ->Add(graphInterPhiZ.release(),"P");
	  mg_trackRZ->Add(graphInterRZ.release(),"P");
	}

      for(size_t id_det = G4Sol::MG01; id_det <= G4Sol::MG17; ++id_det)
        RecoEvent.OldListHits[id_det].resize(RecoEvent.ListHits[id_det].size());

      if(std::get<0>(paraFitRZ) >= 1)
	{
	  if(MDCWireType == 0)
	    {
	      for(size_t iM = 0; iM<id_detMG.size();++iM)
		{

		  int id_detN = std::get<1>(id_detMG[iM]);
		  int id_hit = it_ListHits->second[id_detN];

		  att._logger->debug("Create new WireType 0 : {}, {}",id_detN,G4Sol::nameLiteralDet.begin()[id_detN]);

		  int id_det_reduced = id_detN - G4Sol::MG01;
		  auto [meanZ,meanR,errZ,errR,meanPhi,errPhi] = graphFits[id_det_reduced];

		  att._logger->debug(" meanZ {}  errZ {} | meanR {} errR {} | meanPhi {} errPhi {}",meanZ,errZ,meanR,errR,meanPhi,errPhi);

		  auto TempSeg = RecoEvent.SegmentHit1Ds[id_detN][id_hit];

		  TVector3 wire_side1(TempSeg[0][0],TempSeg[0][1],TempSeg[0][2]);
		  TVector3 wire_side2(TempSeg[TempSeg.size()-1][0],TempSeg[TempSeg.size()-1][1],TempSeg[TempSeg.size()-1][2]);

		  TVector3 WireDir = wire_side2 - wire_side1;
		  TVector3 VecWire = WireDir;
		  VecWire *= 1. / (wire_side2.Z() - wire_side1.Z());

		  //TVector3 PosAtWire = wire_side1 + VecWire * (Mean_CrossRZ[0] - wire_side1.Z());

		  ROOT::Math::Cylindrical3D<double> vecCyl(meanR,meanZ,meanPhi);
		  ROOT::Math::Cartesian3D<double> vecCar(vecCyl);


		  TVectorD hitCoordsNew(3);
		  hitCoordsNew(0) = vecCar.x();
		  hitCoordsNew(1) = vecCar.y();
		  hitCoordsNew(2) = vecCar.z();

		  TMatrixDSym hitCovNew(3);
		  hitCovNew.Zero();
		  hitCovNew(0, 0) = errR;
		  hitCovNew(1, 1) = errZ;
		  hitCovNew(2, 2) = errPhi;

		  hitCovNew *= TMath::Sqrt(2 * hitCovNew.GetNcols());

		  std::vector<ROOT::Math::XYZVector> vecCarAll;
		  ROOT::Math::XYZVector meanVecCar(0.,0.,0.);
		  for(int iC = 0; iC < hitCovNew.GetNcols(); ++iC)
		    {
		      ROOT::Math::RhoZPhiVector vecCyltemp(meanR-hitCovNew(0,iC),meanZ-hitCovNew(1,iC),meanPhi-hitCovNew(2,iC));

		      ROOT::Math::XYZVector vecCartemp(vecCyltemp);
		      vecCarAll.push_back(vecCartemp);
		      meanVecCar += vecCartemp;
		      ROOT::Math::RhoZPhiPoint vecCyltemp2(meanR+hitCovNew(0,iC),meanZ+hitCovNew(1,iC),meanPhi+hitCovNew(2,iC));

		      ROOT::Math::XYZVector vecCartemp2(vecCyltemp);
		      vecCarAll.push_back(vecCartemp2);
		      meanVecCar += vecCartemp2;
		    }

		  meanVecCar *= 1./static_cast<double>(vecCarAll.size());
		  att._logger->debug("E> meanVecCar : {} {} {}", meanVecCar.x(), meanVecCar.y(), meanVecCar.z());

		  TMatrixD hitCov(3,3);
		  hitCov.Zero();
		  for(auto tempV : vecCarAll)
		    {
		      TMatrixD temp1(3, 1);
		      temp1(0, 0) = tempV.x() - meanVecCar.x();
		      temp1(1, 0) = tempV.y() - meanVecCar.y();
		      temp1(2, 0) = tempV.z() - meanVecCar.z();
		      auto temp2 = temp1;
		      temp2.T();
		      hitCov += temp1 * temp2;
		    }

		  hitCov *= 1./static_cast<double>(vecCarAll.size());

		  TMatrixDSym hitCovSymNew(3);
		  hitCovSymNew.Zero();
		  hitCovSymNew(0, 0) = hitCov(0, 0);
		  hitCovSymNew(0, 1) = 0;//hitCov(0, 1);
		  hitCovSymNew(0, 2) = 0;//hitCov(0, 2);
		  hitCovSymNew(1, 0) = 0;//hitCov(1, 0);
		  hitCovSymNew(1, 1) = hitCov(1, 1);
		  hitCovSymNew(1, 2) = hitCov(1, 2);
		  hitCovSymNew(2, 0) = 0;//hitCov(2, 0);
		  hitCovSymNew(2, 1) = hitCov(2, 1);
		  hitCovSymNew(2, 2) = hitCov(2, 2);

		  hitCovSymNew.Print();
		  RecoEvent.OldListHits[id_detN][id_hit] = std::make_unique<genfit::ProlateSpacepointMeasurement>(
														  hitCoordsNew, hitCovSymNew, id_detN, id_hit, nullptr);
		  dynamic_cast<genfit::ProlateSpacepointMeasurement*>(RecoEvent.OldListHits[id_detN][id_hit].get())
		    ->setLargestErrorDirection(WireDir);

		  RecoEvent.OldListHits[id_detN][id_hit].swap(RecoEvent.ListHits[id_detN][id_hit]);

		}
	    }
	  else if(MDCWireType == 1)
	    {
	      for(size_t iM = 0; iM<id_detMG.size();++iM)
		{

		  int id_detN = std::get<1>(id_detMG[iM]);
		  int id_hit = it_ListHits->second[id_detN];

		  att._logger->debug("Create new WireType 1 : {}, {}",id_detN,G4Sol::nameLiteralDet.begin()[id_detN]);

		  genfit::WireMeasurement* currentHit1 = dynamic_cast<genfit::WireMeasurement*>(RecoEvent.ListHits[id_detN][id_hit].get());
		  auto HitCoord = currentHit1->getRawHitCoords();
		  auto HitCov   = currentHit1->getRawHitCov();

		  TVector3 wire_side1_org(HitCoord[0], HitCoord[1], HitCoord[2]);
		  TVector3 wire_side2_org(HitCoord[3], HitCoord[4], HitCoord[5]);

		  att._logger->debug("Wires 1 Org : [{}, {}, {}] Wires 2 Org: [{}, {}, {}]",wire_side1_org.X(),wire_side1_org.Y(),wire_side1_org.Z(),wire_side2_org.X(),wire_side2_org.Y(),wire_side2_org.Z());

		  int id_det_reduced = id_detN - G4Sol::MG01;
		  auto [meanZ,meanR,errZ,errR,meanPhi,errPhi] = graphFits[id_det_reduced];

		  auto TempSeg = RecoEvent.SegmentHit1Ds[id_detN][id_hit];

		  TVector3 wire_side1(TempSeg[0][0],TempSeg[0][1],TempSeg[0][2]);
		  TVector3 wire_side2(TempSeg[TempSeg.size()-1][0],TempSeg[TempSeg.size()-1][1],TempSeg[TempSeg.size()-1][2]);

		  att._logger->debug("Wires 1 New : [{}, {}, {}] Wires 2 New: [{}, {}, {}]",wire_side1.X(),wire_side1.Y(),wire_side1.Z(),wire_side2.X(),wire_side2.Y(),wire_side2.Z());

		  TVector3 WireDir = wire_side2 - wire_side1;
		  TVector3 VecWire = WireDir;
		  VecWire *= 1. / (wire_side2.Z() - wire_side1.Z());

		  double meanX =  meanR * std::cos(meanPhi);
		  double meanY = meanR * std::sin(meanPhi);

		  TVector3 PosAtWire = wire_side1 + VecWire * (meanZ - wire_side1.Z());
		  double newDist = TMath::Hypot(PosAtWire.X()-meanX,PosAtWire.Y()-meanY);

		  TVectorD hitCoordsNew(8);
		  hitCoordsNew(0) = wire_side1[0];
		  hitCoordsNew(1) = wire_side1[1];
		  hitCoordsNew(2) = wire_side1[2];
		  hitCoordsNew(3) = wire_side2[0];
		  hitCoordsNew(4) = wire_side2[1];
		  hitCoordsNew(5) = wire_side2[2];
		  hitCoordsNew(6) = newDist;
		  hitCoordsNew(7) = TMath::Abs(meanZ - wire_side1[2]);

		  att._logger->debug("Wire local r, pos : {} | {} {} {}",hitCoordsNew(6),hitCoordsNew(7), meanZ, errZ);

		  TMatrixDSym hitCovNew(8);
		  const double resolution_dl  = 0.02 ;
		  hitCovNew(6, 6) =  resolution_dl* resolution_dl;
		  hitCovNew(7, 7) = errZ*errZ;

		  RecoEvent.OldListHits[id_detN][id_hit] =
		    std::make_unique<genfit::WirePointMeasurement>(hitCoordsNew, hitCovNew, id_detN, id_hit, nullptr);

		  RecoEvent.OldListHits[id_detN][id_hit].swap(RecoEvent.ListHits[id_detN][id_hit]);

		}
	    }
	}


      auto simTrackPhiZ = std::make_unique<TGraphErrors>();
      simTrackPhiZ->SetNameTitle("simTrackPhiZ","simTrackPhiZ");
      auto simTrackRZ = std::make_unique<TGraphErrors>();
      simTrackRZ->SetNameTitle("simTrackRZ","simTrackRZ");
      auto simTrackZR = std::make_unique<TGraphErrors>();
      simTrackRZ->SetNameTitle("simTrackZR","simTrackZR");

      auto realTrackPhiZ = std::make_unique<TGraphErrors>();
      realTrackPhiZ->SetNameTitle("realTrackPhiZ","realTrackPhiZ");
      auto realTrackRZ = std::make_unique<TGraphErrors>();
      realTrackRZ->SetNameTitle("realTrackRZ","realTrackRZ");
      auto realTrackZR = std::make_unique<TGraphErrors>();
      realTrackZR->SetNameTitle("realTrackZR","realTrackZR");

      std::map<double, tuple< double,double,double,double,double,int> > orderZRPhiTrackR;
      std::map<double, tuple< double,double,double,double,double,int> > orderZRPhiTrackS;

      TVector3 afterV;

      for(size_t idet = 0; idet<it_ListHitsSim->second.size();++idet)
	{
	  if(it_ListHitsSim->second[idet].size()==0)
	    continue;
	  if(!((idet>=G4Sol::MG01 && idet<=G4Sol::MG17) || idet==G4Sol::PSCE ||idet==G4Sol::PSBE||idet==G4Sol::PSFE))
	    continue;

	  afterV.SetXYZ(0.5*(it_ListHitsSim->second[idet][0].hitX+it_ListHitsSim->second[idet][0].hitXexit),
			0.5*(it_ListHitsSim->second[idet][0].hitY+it_ListHitsSim->second[idet][0].hitYexit),
			0.5*(it_ListHitsSim->second[idet][0].hitZ+it_ListHitsSim->second[idet][0].hitZexit));

	  double dx = 0.5*(it_ListHitsSim->second[idet][0].hitX-it_ListHitsSim->second[idet][0].hitXexit);
	  double dy = 0.5*(it_ListHitsSim->second[idet][0].hitY-it_ListHitsSim->second[idet][0].hitYexit);
	  double dz = TMath::Abs(0.5*(it_ListHitsSim->second[idet][0].hitZ-it_ListHitsSim->second[idet][0].hitZexit));

	  att._logger->debug(" sim track : {} {} {} | hit {} exit {} mid {}",afterV.Z(), afterV.Phi(), G4Sol::nameLiteralDet.begin()[idet],it_ListHitsSim->second[idet][0].hitZ, it_ListHitsSim->second[idet][0].hitZexit, it_ListHitsSim->second[idet][0].hitZmid);

	  orderZRPhiTrackS.insert({afterV.Z(), {afterV.Phi()-meanPhi,
		afterV.Perp(),
		dz,
		TMath::Sqrt(TMath::Sq(afterV.Y()*dx)+TMath::Sq(afterV.X()*dy))/afterV.Perp2(),// dPhi
		TMath::Sqrt(TMath::Sq(afterV.X()*dx)+TMath::Sq(afterV.Y()*dy))/afterV.Perp(), // dR
		idet,
	      }
	    });

	  afterV.SetXYZ(it_ListHitsSim->second[idet][0].hitXmid,
			it_ListHitsSim->second[idet][0].hitYmid,
			it_ListHitsSim->second[idet][0].hitZmid);

	  orderZRPhiTrackR.insert({afterV.Z(), {afterV.Phi()-meanPhi,
		afterV.Perp(),
		dz,
		TMath::Sqrt(TMath::Sq(afterV.Y()*dx)+TMath::Sq(afterV.X()*dy))/afterV.Perp2(),// dPhi
		TMath::Sqrt(TMath::Sq(afterV.X()*dx)+TMath::Sq(afterV.Y()*dy))/afterV.Perp(), // dR
		idet,
	      }
	    });
	}

      for(auto [tempZ, tempPhiRDzdRDPhi] : orderZRPhiTrackS)
	{
	  simTrackPhiZ->SetPoint(simTrackPhiZ->GetN(), tempZ, std::get<0>(tempPhiRDzdRDPhi));
	  simTrackPhiZ->SetPointError(simTrackPhiZ->GetN()-1, std::get<2>(tempPhiRDzdRDPhi), std::get<3>(tempPhiRDzdRDPhi));

	  simTrackRZ->SetPoint(simTrackRZ->GetN(), tempZ, std::get<1>(tempPhiRDzdRDPhi));
	  simTrackRZ->SetPointError(simTrackRZ->GetN()-1, std::get<2>(tempPhiRDzdRDPhi), std::get<4>(tempPhiRDzdRDPhi));

	  simTrackZR->SetPoint(simTrackZR->GetN(), std::get<1>(tempPhiRDzdRDPhi), tempZ);
	  simTrackZR->SetPointError(simTrackZR->GetN()-1, std::get<4>(tempPhiRDzdRDPhi), std::get<2>(tempPhiRDzdRDPhi));
	}

      for(auto [tempZ, tempPhiRDzdRDPhi] : orderZRPhiTrackS)
	{
	  if(std::get<5>(tempPhiRDzdRDPhi)>=G4Sol::MG01 && std::get<5>(tempPhiRDzdRDPhi)<=G4Sol::MG17)
	    {
	      int id_red = std::get<5>(tempPhiRDzdRDPhi)-G4Sol::MG01;

	      double tempR = std::get<1>(graphFits[id_red]);
	      double errZ = std::get<2>(graphFits[id_red]);
	      double errR = std::get<3>(graphFits[id_red]);
	      double SimZ = simTrackZR->Eval(tempR);
	      //double SimZErr1 = simTrackZR->Eval(tempR+simTrackZR->GetRMS(1));
	      //double SimZErr2 = simTrackZR->Eval(tempR-simTrackZR->GetRMS(1));
	      double SimZErr1 = simTrackZR->Eval(tempR+errR);
	      double SimZErr2 = simTrackZR->Eval(tempR-errR);
	      double SimZ2 = gRandom->Uniform(TMath::Min(SimZErr1,SimZErr2),TMath::Max(SimZErr1,SimZErr2));

	      double Std_SimZerr = TMath::Sqrt(TMath::Sq(SimZ-SimZErr1)+TMath::Sq(SimZ-SimZErr2))/TMath::Sqrt(0.5);
	      LocalHisto.h_ResidualMDC_dZ1->Fill(id_red,(SimZ2-std::get<0>(graphFits[id_red])));// /TMath::Sqrt(TMath::Sq(std::get<1>(tempPhiRDzdRDPhi))+TMath::Sq(std::get<2>(graphFits[id_red]))));
	      LocalHisto.h_PullMDC_dZ1->Fill(id_red,(SimZ2-std::get<0>(graphFits[id_red]))/Std_SimZerr);// /TMath::Sqrt(TMath::Sq(std::get<1>(tempPhiRDzdRDPhi))+TMath::Sq(std::get<2>(graphFits[id_red]))));
	      LocalHisto.h_RPhiZMDC_Sigma->Fill(id_red,Std_SimZerr);// /TMath::Sqrt(TMath::Sq(std::get<1>(tempPhiRDzdRDPhi))+TMath::Sq(std::get<2>(graphFits[id_red]))));
	      LocalHisto.h_RPhiZMDC_Sigma2->Fill(id_red,errZ);// /TMath::Sqrt(TMath::Sq(std::get<1>(tempPhiRDzdRDPhi))+TMath::Sq(std::get<2>(graphFits[id_red]))));

	      if(id_detMG.size()>=6)
		{
		  LocalHisto.h_ResidualMDC_dZ_More6->Fill(id_red,(SimZ2-std::get<0>(graphFits[id_red])));
		  LocalHisto.h_PullMDC_dZ_More6->Fill(id_red,(SimZ2-std::get<0>(graphFits[id_red]))/Std_SimZerr);
		}
	      if(id_detPSB.size()>0)
		{
		  LocalHisto.h_ResidualMDC_dZ_PSB->Fill(id_red,(SimZ2-std::get<0>(graphFits[id_red])));
		  LocalHisto.h_PullMDC_dZ_PSB->Fill(id_red,(SimZ2-std::get<0>(graphFits[id_red]))/Std_SimZerr);
		}
	      if(id_detPSBE.size()>0)
		{
		  LocalHisto.h_ResidualMDC_dZ_PSBE->Fill(id_red,(SimZ2-std::get<0>(graphFits[id_red])));
		  LocalHisto.h_PullMDC_dZ_PSBE->Fill(id_red,(SimZ2-std::get<0>(graphFits[id_red]))/Std_SimZerr);
		}
	      if(id_detPSFE.size()>0)
		{
		  LocalHisto.h_ResidualMDC_dZ_PSFE->Fill(id_red,(SimZ2-std::get<0>(graphFits[id_red])));
		  LocalHisto.h_PullMDC_dZ_PSFE->Fill(id_red,(SimZ2-std::get<0>(graphFits[id_red]))/Std_SimZerr);
		}
	    }
	}
      for(auto [tempZ, tempPhiRDzdRDPhi] : orderZRPhiTrackR)
	{
	  realTrackPhiZ->SetPoint(realTrackPhiZ->GetN(), tempZ, std::get<0>(tempPhiRDzdRDPhi));
	  realTrackPhiZ->SetPointError(realTrackPhiZ->GetN()-1, std::get<2>(tempPhiRDzdRDPhi), std::get<3>(tempPhiRDzdRDPhi));

	  realTrackRZ->SetPoint(realTrackRZ->GetN(), tempZ, std::get<1>(tempPhiRDzdRDPhi));
	  realTrackRZ->SetPointError(realTrackRZ->GetN()-1, std::get<2>(tempPhiRDzdRDPhi), std::get<4>(tempPhiRDzdRDPhi));

	  realTrackZR->SetPoint(realTrackZR->GetN(), std::get<1>(tempPhiRDzdRDPhi), tempZ);
	  realTrackZR->SetPointError(realTrackZR->GetN()-1, std::get<4>(tempPhiRDzdRDPhi), std::get<2>(tempPhiRDzdRDPhi));
	}

      for(auto [tempZ, tempPhiRDzdRDPhi] : orderZRPhiTrackR)
	{
	  if(std::get<5>(tempPhiRDzdRDPhi)>=G4Sol::MG01 && std::get<5>(tempPhiRDzdRDPhi)<=G4Sol::MG17)
	    {
	      int id_red = std::get<5>(tempPhiRDzdRDPhi)-G4Sol::MG01;

	      double tempR = std::get<1>(graphFits[id_red]);
	      double SimZ = realTrackZR->Eval(tempR);
	      //double SimZErr1 = realTrackZR->Eval(tempR+realTrackZR->GetRMS(1));
	      //double SimZErr2 = realTrackZR->Eval(tempR-realTrackZR->GetRMS(1));

	      //double Std_SimZerr = TMath::Sqrt(TMath::Sq(SimZ-SimZErr1)+TMath::Sq(SimZ-SimZErr2));
	      LocalHisto.h_ResidualMDC_dZ2->Fill(id_red,(SimZ-std::get<0>(graphFits[id_red])));
	      LocalHisto.h_PullMDC_dZ2->Fill(id_red,(SimZ-std::get<0>(graphFits[id_red]))/TMath::Sqrt(TMath::Sq(std::get<1>(tempPhiRDzdRDPhi))+TMath::Sq(std::get<2>(graphFits[id_red]))));
	    }
	}



      if(OutputEvents && id_detMG.size()>0)
	{
	  simTrackPhiZ->SetMarkerStyle(20);
	  simTrackPhiZ->SetMarkerColor(2);
	  simTrackPhiZ->SetLineColor(2);

	  simTrackRZ->SetMarkerStyle(20);
	  simTrackRZ->SetMarkerColor(2);
	  simTrackRZ->SetLineColor(2);

	  realTrackPhiZ->SetMarkerStyle(21);
	  realTrackPhiZ->SetMarkerColor(6);
	  realTrackPhiZ->SetLineColor(6);

	  realTrackRZ->SetMarkerStyle(21);
	  realTrackRZ->SetMarkerColor(6);
	  realTrackRZ->SetLineColor(6);

	  mg_trackPhiZ->Add(simTrackPhiZ.release(),"PL");
	  mg_trackRZ->Add(simTrackRZ.release(),"PL");

	  mg_trackPhiZ->Add(realTrackPhiZ.release(),"PL");
	  mg_trackRZ->Add(realTrackRZ.release(),"PL");

	  t_phiZ->Fill();

	  auto* listG1 = mg_trackPhiZ->GetListOfGraphs();
	  listG1->Delete();
	  auto* listG1_1 = mg_trackPhiZNoCorr->GetListOfGraphs();
	  listG1_1->Delete();
	  auto* listG2 = mg_trackRZ->GetListOfGraphs();
	  listG2->Delete();
	}
      else
	{
	  auto* listG1 = mg_trackPhiZ->GetListOfGraphs();
	  listG1->Delete();
	  auto* listG1_1 = mg_trackPhiZNoCorr->GetListOfGraphs();
	  listG1_1->Delete();
	  auto* listG2 = mg_trackRZ->GetListOfGraphs();
	  listG2->Delete();
	}

      att._logger->debug("RPhiZ MDC > @@@");
    }

  ++tempEvent;

  return 0;
}

template class TRPhiZTrackMDC<MCAnaEventG4Sol>;
template class TRPhiZTrackMDC<Ana_WasaEvent>;



int MDC_PRZ::f_checkDiscontinuity(TGraphErrors* g, double limit)
{
  for (int i = 1; i < g->GetN(); ++i)
    {
      double x1, y1, x2, y2;
      g->GetPoint(i-1, x1, y1);
      g->GetPoint(i, x2, y2);
      if (TMath::Abs(y2 - y1) > limit)
	return i;
    }
  return -1;
};

std::tuple<int,double,double> MDC_PRZ::f_calInterset(std::array<double,4>& x, std::array<double,4> y)
{
  double denom = (x[0] - x[1]) * (y[2] - y[3]) - (y[0] - y[1]) * (x[2] - x[3]);
  // att._logger->debug("denom : {}",denom);
  if (TMath::Abs(denom) < 1e-3)
    {
      return std::make_tuple(0,0.,0.);
      //continue;
    }
  double x_num   = (x[0] * y[1] - y[0] * x[1]) * (x[2] - x[3]) - (x[0] - x[1]) * (x[2] * y[3] - y[2] * x[3]);
  double y_num = (x[0] * y[1] - y[0] * x[1]) * (y[2] - y[3]) - (y[0] - y[1]) * (x[2] * y[3] - y[2] * x[3]);
  double x_intersect   = x_num / denom;
  double y_intersect = y_num / denom;

  return std::make_tuple(1,x_intersect,y_intersect);//,r_intersect1,r_intersect2);
};

void MDC_PRZ::f_MeanCov(const std::vector<std::tuple<double,double>>& ResZPR, TVectorD& mean, TMatrixD& cov)
{
  mean.ResizeTo(2);
  mean(0) = 0.;
  mean(1) = 0.;

  for(auto [z,pr] : ResZPR)
    {
      mean(0) += z;
      mean(1) += pr;
    }

  mean *= 1./static_cast<double>(ResZPR.size());
  cov.ResizeTo(2,2);
  cov.Zero();

  for(auto [z,pr] : ResZPR)
    {
      TMatrixD temp1(2, 1);
      temp1(0, 0) = mean(0) - z;
      temp1(1, 0) = mean(1) - pr;

      auto temp2 = temp1;
      temp2.T();
      cov += temp1 * temp2;
    }
  cov *= 1. / static_cast<double>(ResZPR.size());
};


std::tuple<int,double,double,double,double,double,double,double,double> MDC_PRZ::f_interset(const TList* listPZ, const TList* listRZ, int id_1, int id_2)
{

  TGraphErrors* graph1 = dynamic_cast<TGraphErrors*>(listPZ->At(id_1));
  TGraphErrors* graph2 = dynamic_cast<TGraphErrors*>(listPZ->At(id_2));
  TGraphErrors* graph1RZ = dynamic_cast<TGraphErrors*>(listRZ->At(id_1));
  TGraphErrors* graph2RZ = dynamic_cast<TGraphErrors*>(listRZ->At(id_2));

  // int idDisc_1 = MDC_PRZ::f_checkDiscontinuity(graph1, TMath::PiOver2());
  // int idDisc_2 = MDC_PRZ::f_checkDiscontinuity(graph2, TMath::PiOver2());

  std::pair<int,int> minmaxPoint1, minmaxPoint2;
  // if(idDisc_1==-1)
  minmaxPoint1 = {0, graph1->GetN()-1};
  // else
  //   minmaxPoint1 = {0, idDisc_1 -1};

  // if(idDisc_2==-1)
  minmaxPoint2 = {0, graph2->GetN()-1};
  // else
  //   minmaxPoint2 = {0, idDisc_2 -1};


  double z1 = graph1->GetPointX(minmaxPoint1.first);
  double phi1 = graph1->GetPointY(minmaxPoint1.first);
  double z2 = graph1->GetPointX(minmaxPoint1.second);
  double phi2 = graph1->GetPointY(minmaxPoint1.second);
  double z3 = graph2->GetPointX(minmaxPoint2.first);
  double phi3 = graph2->GetPointY(minmaxPoint2.first);
  double z4 = graph2->GetPointX(minmaxPoint2.second);
  double phi4 = graph2->GetPointY(minmaxPoint2.second);

  // att._logger->debug("graph1 : {}",graph1->GetName());
  // graph1->Print();
  // att._logger->debug("graph2 : {}", graph2->GetName());
  // graph1->Print();
  // att._logger->debug("coordinates g1 : ({}, {}) ; ({}, {}) | g2: ({}, {}) ; ({}, {}) ",z1,phi1,z2,phi2,z3,phi3,z4,phi4);

  std::array<double,4> allZ = {z1,z2,z3,z4}, allPhi = {phi1,phi2,phi3,phi4};


  auto [res,z_intersect,phi_intersect] = MDC_PRZ::f_calInterset(allZ,allPhi);
  if(res == 0)
    return std::make_tuple(0, 0.,0.,0.,0.,0.,0.,0.,0.);

  if(z_intersect<0.99*z1 || z_intersect > 1.01 *z2)
    return std::make_tuple(-1, 0.,0.,0.,0.,0.,0.,0.,0.);


  double phi1Er = graph1->GetErrorY(minmaxPoint1.first);
  double phi2Er = graph1->GetErrorY(minmaxPoint1.second);
  double phi3Er = graph2->GetErrorY(minmaxPoint2.first);
  double phi4Er = graph2->GetErrorY(minmaxPoint2.second);

  std::vector<std::tuple<double,double> > ResZPhi;
  std::vector<std::array<double,4> > ErrAllPhi = { {-phi1Er,0.,0.,0.},
						   { phi1Er,0.,0.,0.},
						   {0.,-phi2Er,0.,0.},
						   {0., phi2Er,0.,0.},
						   {0.,0.,-phi3Er,0.},
						   {0.,0., phi3Er,0.},
						   {0.,0.,0.,-phi4Er},
						   {0.,0.,0., phi4Er}};

  for(size_t iR = 0; iR < ErrAllPhi.size();++iR)
    {
      std::array<double,4> TempAllPhi = allPhi;
      for(size_t iE = 0 ; iE<ErrAllPhi[iR].size();++iE)
	TempAllPhi[iE] += ErrAllPhi[iR][iE];

      auto [res, tempZ, tempPhi] = MDC_PRZ::f_calInterset(allZ,TempAllPhi);
      ResZPhi.push_back({tempZ,tempPhi});
    }


  TVectorD meanAll(2);
  TMatrixD covAll(2,2);

  MDC_PRZ::f_MeanCov(ResZPhi,meanAll,covAll);


  //double r_intersect1 = graph1RZ->Eval(z_intersect);
  //double r_intersect2 = graph2RZ->Eval(z_intersect);

  double r_intersect1_low = graph1RZ->Eval(z_intersect-TMath::Sqrt(covAll(0,0)));
  double r_intersect1_high = graph1RZ->Eval(z_intersect+TMath::Sqrt(covAll(0,0)));
  double r_intersect2_low = graph2RZ->Eval(z_intersect-TMath::Sqrt(covAll(0,0)));
  double r_intersect2_high = graph2RZ->Eval(z_intersect+TMath::Sqrt(covAll(0,0)));

  double r1 = 0.5*(r_intersect1_low+r_intersect1_high);
  double r2 = 0.5*(r_intersect2_low+r_intersect2_high);
  double r1_e = 0.5*(TMath::Sq(r1-r_intersect1_low)+TMath::Sq(r1-r_intersect1_high));
  double r2_e = 0.5*(TMath::Sq(r2-r_intersect2_low)+TMath::Sq(r2-r_intersect2_high));

  double RErrGraph1 = graph1RZ->GetErrorY(0);
  double RErrGraph2 = graph2RZ->GetErrorY(0);

  r1_e = TMath::Sqrt(TMath::Sq(r1_e)+TMath::Sq(RErrGraph1));
  r2_e = TMath::Sqrt(TMath::Sq(r2_e)+TMath::Sq(RErrGraph2));

  // std::cout<<"Check Error : z "<< z_intersect<<" phi "<< phi_intersect <<" | z "<< meanAll(0)<< "phi "<< meanAll(1)<<"| +- z "<< TMath::Sqrt(covAll(0,0))<<" +- phi"<<TMath::Sqrt(covAll(1,1)) <<"\n";
  // std::cout<<" r1 e "<<r_intersect1_low<< " "<<r_intersect1_high<<"\n";
  // std::cout<<" r2 e "<<r_intersect2_low<< " "<<r_intersect2_high<<"\n";
  // std::cout<<"r mean "<<r1<<" "<<r2<<"\n";
  // std::cout<<"r std  "<<r1_e<<" "<<r2_e<<"\n";
  //return std::make_tuple(1,z_intersect,phi_intersect,r_intersect1,r_intersect2, );
  return std::make_tuple(1, meanAll(0), meanAll(1) ,r1 ,r2, TMath::Sqrt(covAll(0,0)), TMath::Sqrt(covAll(1,1)), TMath::Sqrt(r1_e),TMath::Sqrt(r2_e));

  // att._logger->debug(" intersect: {}, {}",z_intersect,phi_intersect);
};

std::tuple<int,double,double,double,double,double,double> MDC_PRZ::f_intersetTrack(const TList* listRZ, int id_1, double p0, double p1, const TMatrixD& covP01)
{

  TGraphErrors* graph1RZ = dynamic_cast<TGraphErrors*>(listRZ->At(id_1));
  double min_Z = graph1RZ->GetPointX(0);
  double max_Z = graph1RZ->GetPointX(graph1RZ->GetN()-1);

  double TrackR_minZ = p0 + min_Z*p1;
  double TrackR_maxZ = p0 + max_Z*p1;

#ifdef DEBUGRPZ
  std::cout<<"In tube :"<<" "<<id_1<<" "<<p0<<" "<<p1<<"\n";
  covP01.Print();
#endif
  auto f_TrackTube = [](double p0,double p1, const TMatrixD& cov, double z) -> std::tuple<double,double> {

    auto f_sqrtMat2x2 = [](const TMatrixD& m) -> TMatrixD {
      //double det_m   = m.Determinant();
      //double trace_m = m(0, 0) + m(1, 1);

      //double Norm = TMath::Sqrt(det_m);

      //TMatrixD U(2, 2);
      //U.UnitMatrix();

      //auto mSqrt = (m + Norm * U);
      //mSqrt *= 1. / det_m;
      //return mSqrt;

      TVectorD EigVal;
      TMatrixD EigVec = m.EigenVectors(EigVal);
      TMatrixD m_new(m.GetNrows(), m.GetNcols());
      m_new.UnitMatrix();
      for(int i = 0; i < m.GetNrows(); ++i)
	m_new(i, i) = EigVal(i);
      m_new.Sqrt();
      auto EigVecInv = EigVec;
      EigVecInv.Invert();
      auto mSqrt = EigVec * m_new * EigVecInv;
      return mSqrt;
    };

    // if(!cov.IsSymmetric())
    //   {
    // 	std::cout<<"E> f_TrackTube : cov not symmetric ! ";
    // 	cov.Print();
    //   }
    TMatrixD SubCovSqrt = f_sqrtMat2x2(cov);
    SubCovSqrt *= TMath::Sqrt(2 * cov.GetNcols());

#ifdef DEBUGRPZ
    SubCovSqrt.Print();
#endif

    std::vector<double> Rs;
    double mean_R = 0.;
    for(int iC = 0; iC < SubCovSqrt.GetNcols(); ++iC)
      {
	double p0_eP = p0 + SubCovSqrt(0, iC);
	double p1_eP = p1 + SubCovSqrt(1, iC);

	double p0_eN = p0 - SubCovSqrt(0, iC);
	double p1_eN = p1 - SubCovSqrt(1, iC);

#ifdef DEBUGRPZ
	std::cout<<" iC "<<iC<<" +:"<<p0_eP<<" "<<p1_eP<<" -:"<<p0_eN<<" "<<p1_eN<<"\n";
#endif
	mean_R += p0_eP + z*p1_eP;
	Rs.push_back(p0_eP + z*p1_eP);
	mean_R += p0_eN + z*p1_eN;
	Rs.push_back(p0_eN + z*p1_eN);
      }
    mean_R *= 1./static_cast<double>(Rs.size());

    double Rstd = 0;
    for(auto Rsig : Rs)
      {
	Rstd += TMath::Sq(mean_R - Rsig);
#ifdef DEBUGRPZ
	std::cout<<" Rs :"<<Rsig<<"\n";
#endif
      }
    Rstd *= 1./static_cast<double>(Rs.size());

    return std::make_tuple(mean_R,TMath::Sqrt(Rstd));
  };

#ifdef DEBUGRPZ
  std::cout<<"intersetTrack : graph1RZ :"<< graph1RZ->GetName() << "\n";
  //graph1RZ->Print();
#endif
  int res, iP1 = 0;
  double z_intersect, r_intersect;
  for(int iP = 1; iP < graph1RZ->GetN(); ++iP)
    {
      double z1 = graph1RZ->GetPointX(iP-1);
      double R1 = graph1RZ->GetPointY(iP-1);
      double z2 = graph1RZ->GetPointX(iP);
      double R2 = graph1RZ->GetPointY(iP);

#ifdef DEBUGRPZ
      //graph1RZ->Print();
      std::cout<<" coordinates g1RZ :("<<z1<<", "<<R1<<") ; ("<<z2<<", "<<R2<<") track:"<<min_Z<<", "<<TrackR_minZ<<") ; ("<<max_Z<<", "<<TrackR_maxZ<<") \n";
#endif
      std::array<double,4> allZ = {min_Z,max_Z,z1,z2}, allR = {TrackR_minZ,TrackR_maxZ,R1,R2};

      std::tie(res,z_intersect,r_intersect) = MDC_PRZ::f_calInterset(allZ,allR);
      if(res == 0)
	continue;
      else
	{
	  if(z_intersect< z1 || z_intersect>z2)
	    if(z_intersect > min_Z && z_intersect < max_Z)
	      continue;

	  iP1 = iP-1;
	  break;
	}
    }

  if(res==0)
    std::cout<<"E> not interset between Pol1 and RZ wire !\n";
#ifdef DEBUGRPZ
  else
    std::cout<<"Find Intersect : graphMG: "<<graph1RZ->GetName()<<" "<<z_intersect<<" "<<r_intersect<<"\n";
#endif

  auto [Rm1, R1Er] = f_TrackTube(p0,p1,covP01,min_Z);
  auto [Rm2, R2Er] = f_TrackTube(p0,p1,covP01,max_Z);

#ifdef DEBUGRPZ
  std::cout<<"!> tubeTrack "<<Rm1<<" "<<R1Er<<" "<<Rm2<<" "<<R2Er<<"\n";
#endif

  double R3Er = graph1RZ->GetErrorY(iP1);
  double R4Er = graph1RZ->GetErrorY(iP1+1);

  std::vector<std::tuple<double,double> > ResZR;
  std::vector<std::array<double,4> > ErrAllR = { {-R1Er,0.,0.,0.},
						 { R1Er,0.,0.,0.},
						 {0.,-R2Er,0.,0.},
						 {0., R2Er,0.,0.},
						 {0.,0.,-R3Er,0.},
						 {0.,0., R3Er,0.},
						 {0.,0.,0.,-R4Er},
						 {0.,0.,0., R4Er}};

  std::array<double,4> allZ = {min_Z,max_Z,graph1RZ->GetPointX(iP1),graph1RZ->GetPointX(iP1+1)},
                       allR = {  Rm1,  Rm2,graph1RZ->GetPointY(iP1),graph1RZ->GetPointY(iP1+1)};

  std::tuple<double, double, double> best_resRZ = {10000.,0.,0.};

  for(size_t iR = 0; iR < ErrAllR.size();++iR)
    {
      int resFinal = 0;
      //for(int iP = TMath::Max(0, iP1-1); iP < TMath::Min(graph1RZ->GetN(),iP1+1) ; ++iP)
      //{
      // int iP = iP1;
      // allZ[2] = graph1RZ->GetPointX(iP);
      // allZ[3] = graph1RZ->GetPointX(iP+1);

      // allR[2] = graph1RZ->GetPointY(iP);
      // allR[3] = graph1RZ->GetPointY(iP+1);

      std::array<double,4> TempAllR = allR;

      for(size_t iE = 0 ; iE<ErrAllR[iR].size();++iE)
	TempAllR[iE] += ErrAllR[iR][iE];

      auto [res, tempZ, tempR] = MDC_PRZ::f_calInterset(allZ,TempAllR);
#ifdef DEBUGRPZ
      std::cout<<"loop intersect iR:"<<iR<<" "<<res<<" "<<tempZ<<" "<<tempR<<" | [z1 z2] < "<<allZ[2]<<" "<<allZ[3]<<"\n";

      if(TMath::Abs(z_intersect-tempZ)> 2 || TMath::Abs(r_intersect-tempR)> 2)
       	std::cout<<"E> interset and mean/error calculation different :"<< iR<<" | "<<z_intersect<<" "<<tempZ<<" | "<<r_intersect<<" "<<tempR<<" | "<<
       	  allZ[0]<<" "<<allZ[1]<<" "<<allZ[2]<<" "<<allZ[3]<<" / "<<TempAllR[0]<<" "<<TempAllR[1]<<" "<<TempAllR[2]<<" "<<TempAllR[3]<<"\n";
#endif
      if(res == 0)
	{
	  resFinal+=10;
	  //continue;
	}
      else
	{
	  double temp_diffZ = TMath::Abs(z_intersect-tempZ);
	  if(temp_diffZ<std::get<0>(best_resRZ))
	    best_resRZ = {temp_diffZ,tempZ,tempR};
	  if(temp_diffZ > 10.)
	    {
	      ++resFinal;
	      continue;
	    }

#ifdef DEBUGRPZ
	  std::cout<<"  -- Save values "<<res<<" "<<tempZ<<" "<<tempR<<"\n";
#endif
	  ResZR.push_back({tempZ,tempR});
	  resFinal = -1;
	  //break;
	}


#ifdef DEBUGRPZ
      if(resFinal != -1)
	std::cout<<"E> segment RZ intersection with errors failed ! ["<<iR<<"] "<<graph1RZ->GetName()<<" "<<resFinal<<" "<<TMath::Max(0, iP1-10)<<" "<< TMath::Min(graph1RZ->GetN(),iP1+10)<<"\n";
#endif
    }
  if(ResZR.size()==0)
    ResZR.push_back({std::get<1>(best_resRZ),std::get<2>(best_resRZ)});

#ifdef DEBUGRPZ
  std::cout<<" !> out :"<<ResZR.size()<<"\n";
#endif
  TVectorD meanAll(2);
  TMatrixD covAll(2,2);

  MDC_PRZ::f_MeanCov(ResZR,meanAll,covAll);

#ifdef DEBUGRPZ
  if(TMath::Abs(z_intersect-meanAll(0))>0.5 || TMath::Abs(r_intersect-meanAll(1))>0.5)
    std::cout<<"E> interset and mean/error calculation different :"<<z_intersect<<" "<<meanAll(0)<<" | "<<r_intersect<<" "<<meanAll(1)<<"\n";
#endif
  return std::make_tuple(1, z_intersect,r_intersect, meanAll(0), meanAll(1), TMath::Sqrt(covAll(0,0)), TMath::Sqrt(covAll(1,1)));

};
