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

#include <numeric>
#include <set>
#include <sstream>
#include <tuple>
#include <list>

#include "TGraphErrors.h"
#include "TFitResult.h"

using namespace std;
using namespace G4Sol;

template<class Out>
TRPhiZTrackMDC<Out>::TRPhiZTrackMDC(const THyphiAttributes& attribut) : TDataProcessInterface<Out>("CheckFiberTrack"), att(attribut)
{

  OutputEvents = att.RF_OutputEvents;
  std::string temp_name_out       = att.Config.Get<std::string>("Output_Namefile");
  std::string temp_file_base_name = temp_name_out.substr(0, temp_name_out.find_last_of('.'));

  temp_file_base_name += "PhiZTest.root";
  namefilePhiZ = temp_file_base_name;

  RZfit = att.RPZ_RZfit;

  if(OutputEvents)
    {
      att._logger->info("RF: OutputEvent set on !");
      f_phiZ = new TFile(namefilePhiZ, "RECREATE");
      t_phiZ = new TTree("PhiZTree", "Finder output tracks");

      t_phiZ->SetDirectory(f_phiZ);
      // t_phiZ->AutoSave();

      mg_trackPhiZ = new TMultiGraph;
      mg_trackPhiZ->SetNameTitle("PhiZTrack","PhiZTrack");
      mg_trackRZ = new TMultiGraph;
      mg_trackRZ->SetNameTitle("RZTrack","RZTrack");
      tempEvent = 0;
      t_phiZ->Branch("PhiZTrack", "TMultiGraph", &mg_trackPhiZ, 12800, 0);
      t_phiZ->Branch("RZTrack", "TMultiGraph", &mg_trackRZ, 12800, 0);
      t_phiZ->Branch("nEvent",&tempEvent);
    }
  else
    {
      mg_trackPhiZ = new TMultiGraph;
      mg_trackPhiZ->SetNameTitle("PhiZTrack","PhiZTrack");
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

  LocalHisto.h_ResidualMDC_dZ1 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualMDC_dZ1);
  LocalHisto.h_ResidualMDC_dZ2 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualMDC_dZ2);
  LocalHisto.h_RPhiZMDC_Sigma = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_RPhiZMDC_Sigma);

  LocalHisto.h_ResidualMDC_dZ_PSB = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualMDC_dZ_PSB);
  LocalHisto.h_ResidualMDC_dZ_PSBE = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualMDC_dZ_PSBE);
  LocalHisto.h_ResidualMDC_dZ_PSFE = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualMDC_dZ_PSFE);
  LocalHisto.h_ResidualMDC_dZ_More6 = this->AnaHisto->CloneAndRegister(this->AnaHisto->h_ResidualMDC_dZ_More6);

}

template<class Out>
int TRPhiZTrackMDC<Out>::CheckTrackFinding(const FullRecoEvent& RecoEvent)
{

  int ntrack = -1;
  for(auto it_trackInfo : RecoEvent.TrackInfo)
    {
      ntrack++;
      const int id_track = it_trackInfo.first;
      auto it_ListHits   = RecoEvent.TrackDAF.find(id_track);
      auto it_ListHitsSim = RecoEvent.TrackDAFSim.find(id_track);

      std::vector<std::tuple<int,int> > id_detMG;
      std::vector<std::tuple<int,int> > id_detPSB;
      std::vector<std::tuple<int,int> > id_detPSFE;
      std::vector<std::tuple<int,int> > id_detPSBE;
      int index_mg = 0;
      for(size_t id_det = 0; id_det < it_ListHits->second.size(); ++id_det)
        {
          int id_hit = it_ListHits->second[id_det];
          if(id_hit < 0)
            continue;

	  if(RecoEvent.SegmentHit1Ds[id_det][id_hit].size()==0)
	    continue;

	  if(id_det>=G4Sol::MG01 && id_det<=G4Sol::MG17)
	    {
	      id_detMG.push_back({index_mg,id_det});
	      // att._logger->debug("MDC found : {}, {}",id_det,G4Sol::nameLiteralDet.begin()[id_det]);
	    }
	  if(id_det==G4Sol::PSCE)
	    id_detPSB.push_back({index_mg,id_det});
	  if(id_det==G4Sol::PSBE)
	    id_detPSBE.push_back({index_mg,id_det});
	  if(id_det==G4Sol::PSFE)
	    id_detPSFE.push_back({index_mg,id_det});


	  auto TempSeg = RecoEvent.SegmentHit1Ds[id_det][id_hit];

	  TVector3 afterV;

	  auto graphDetPhiZ = std::make_unique<TGraphErrors>();
	  auto graphDetRZ = std::make_unique<TGraphErrors>();
	  for(size_t iS=0;iS<TempSeg.size();++iS)
	    {
	      afterV.SetXYZ(TempSeg[iS][0],TempSeg[iS][1],TempSeg[iS][2]);
	      double tempPhi = afterV.Phi();

	      graphDetPhiZ->SetPoint(iS,TempSeg[iS][2],tempPhi);
	      graphDetPhiZ->SetPointError(iS,TempSeg[iS][5],TempSeg[iS][4]);

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
	  mg_trackRZ->Add(graphDetRZ.release(),"PL");

	  ++index_mg;
	}

      // for(int id_det : id_detMG)
      // 	att._logger->debug("list MDC : {}",id_det);

      auto* listPZ = mg_trackPhiZ->GetListOfGraphs();
      auto* listRZ = mg_trackRZ->GetListOfGraphs();

      auto graphInterPhiZ = std::make_unique<TGraphErrors>();
      auto graphInterRZ = std::make_unique<TGraphErrors>();

      graphInterPhiZ->SetNameTitle("IntersectionPhiZ","IntersectionPhiZ");
      graphInterRZ->SetNameTitle("IntersectionRZ","IntersectionRZ");

      std::unordered_map<int,std::vector<std::tuple<double,double,double,double> > > MG_posZestimated;
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

	  if(ret==0)
	    continue;

	  graphInterPhiZ->SetPoint(graphInterPhiZ->GetN(),z_intersect,phi_intersect);
	  graphInterPhiZ->SetPointError(graphInterPhiZ->GetN()-1,z_e,phi_e);

	  graphInterRZ->SetPoint(graphInterRZ->GetN(),z_intersect, r_intersect1);
	  graphInterRZ->SetPointError(graphInterRZ->GetN()-1,z_e, r1_e);
	  graphInterRZ->SetPoint(graphInterRZ->GetN(),z_intersect, r_intersect2);
	  graphInterRZ->SetPointError(graphInterRZ->GetN()-1,z_e, r2_e);

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

	  if(ret==0)
	    continue;

	  graphInterPhiZ->SetPoint(graphInterPhiZ->GetN(),z_intersect,phi_intersect);
	  graphInterPhiZ->SetPointError(graphInterPhiZ->GetN()-1,z_e,phi_e);
	  graphInterRZ->SetPoint(graphInterRZ->GetN(),z_intersect, r_intersect1);
	  graphInterRZ->SetPointError(graphInterRZ->GetN()-1,z_e, r1_e);
	  graphInterRZ->SetPoint(graphInterRZ->GetN(),z_intersect, r_intersect2);
	  graphInterRZ->SetPointError(graphInterRZ->GetN()-1,z_e, r2_e);

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

      std::vector<std::tuple<double,double,double,double> > graphFits(17);
      if(RZfit && graphInterRZ->GetN() !=0)
	{
	  att._logger->debug("RZfit : IntersectRZ:{} | mg_trackRZ: {}", graphInterRZ->GetN(), mg_trackRZ->GetListOfGraphs()->GetEntries());
	  if(id_detPSB.size()>0)
	    {
	      std::vector<std::tuple<double,double,double,double> > posAllPSB;
	      for(auto [id1,id2] : id_detPSB)
		{
		  //std::cout<<"!> PSCE : "<<id_detPSB.size()<<"x: "<<dynamic_cast<TGraph*>(listRZ->At(std::get<0>(id_detPSB[0])))->GetMean(1)<<"y: "<<dynamic_cast<TGraph*>(listRZ->At(std::get<0>(id_detPSB[0])))->GetMean(2)<<"xe: "<<dynamic_cast<TGraph*>(listRZ->At(std::get<0>(id_detPSB[0])))->GetRMS(1)<<"ye: "<<dynamic_cast<TGraphErrors*>(listRZ->At(std::get<0>(id_detPSB[0])))->GetErrorY(0)<<"\n";
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
	    }

	  TFitResultPtr fitR = graphInterRZ->Fit("pol1","SN0Q");
	  TFitResult* Rf = fitR.Get();
	  double p0 = Rf->Parameter(0);
	  double p1 = Rf->Parameter(1);
	  TMatrixD cov = Rf->GetCovarianceMatrix();
	  //double chi2  = fitR->Chi2();

	  auto graphRZfit1 = std::make_unique<TGraphErrors>();
	  graphRZfit1->SetNameTitle("FitIntRZ","FitIntRZ");
	  auto graphRZfit2 = std::make_unique<TGraphErrors>();
	  graphRZfit2->SetNameTitle("FitErrRZ","FitErrRZ");

	  for( auto [id_mg, id_det] : id_detMG)
	    {
	      auto [res, Zin, Rin, meanZ, meanR, errZ, errR] = MDC_PRZ::f_intersetTrack(listRZ,id_mg, p0,p1,cov);
	      graphRZfit1->SetPoint(graphRZfit1->GetN(),Zin,Rin);

	      graphRZfit2->SetPoint(graphRZfit2->GetN(),meanZ,meanR);
	      graphRZfit2->SetPointError(graphRZfit2->GetN()-1,errZ,errR);

	      int id_det_reduced = id_det - G4Sol::MG01;
	      graphFits[id_det_reduced] = {meanZ,meanR,errZ,errR};
	    }

	  graphRZfit1->SetMarkerStyle(23);
	  graphRZfit1->SetMarkerSize(1.2);
	  graphRZfit1->SetMarkerColor(6);

	  graphRZfit2->SetMarkerStyle(23);
	  graphRZfit2->SetMarkerSize(1.2);
	  graphRZfit2->SetMarkerColor(7);

	  if(OutputEvents && graphRZfit1->GetN()!=0)
	    {
	      mg_trackRZ->Add(graphRZfit1.release(),"P");
	      mg_trackRZ->Add(graphRZfit2.release(),"P");
	    }
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

      auto graphResultZ = std::make_unique<TGraphErrors>();
      graphResultZ->SetNameTitle("IntersectionResultZ","IntersectionResultZ");

      for(const auto& it_MG : MG_posZestimated)
	{
	  if(MG_posZestimated.size() == 1)
	    continue;

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

	  if(TMath::Abs(r_check-std::get<1>(it_MG.second[0]))>0.5)
	    {
	      att._logger->debug("E> CheckIntersectionResultZ : wrong r_intersect ! {}", G4Sol::nameLiteralDet.begin()[it_MG.first]);
	      for(auto [tempZ, tempR, tempZ_e, tempR_e] : it_MG.second)
		att._logger->debug("E> z {} r {} | +- z {} +- r {}",tempZ, tempR, tempZ_e, tempR_e);
	    }

	  graphResultZ->SetPoint(graphResultZ->GetN(), z_mean,r_check);
	  graphResultZ->SetPointError(graphResultZ->GetN()-1, TMath::Sqrt(sig_Z), TMath::Sqrt(sig_R));
	}

      graphResultZ->SetMarkerStyle(22);
      graphResultZ->SetMarkerSize(1.2);
      graphResultZ->SetMarkerColor(4);

      if(OutputEvents && graphResultZ->GetN()!=0)
	mg_trackRZ->Add(graphResultZ.release(),"P");

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

	  afterV.SetXYZ(0.5*(it_ListHitsSim->second[idet][0].hitX+it_ListHitsSim->second[idet][0].hitXexit),
			0.5*(it_ListHitsSim->second[idet][0].hitY+it_ListHitsSim->second[idet][0].hitYexit),
			0.5*(it_ListHitsSim->second[idet][0].hitZ+it_ListHitsSim->second[idet][0].hitZexit));

	  double dx = 0.5*(it_ListHitsSim->second[idet][0].hitX-it_ListHitsSim->second[idet][0].hitXexit);
	  double dy = 0.5*(it_ListHitsSim->second[idet][0].hitY-it_ListHitsSim->second[idet][0].hitYexit);
	  double dz = TMath::Abs(0.5*(it_ListHitsSim->second[idet][0].hitZ-it_ListHitsSim->second[idet][0].hitZexit));

	  orderZRPhiTrackS.insert({afterV.Z(), {afterV.Phi(),
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

	  orderZRPhiTrackR.insert({afterV.Z(), {afterV.Phi(),
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
	      double SimZ = simTrackZR->Eval(tempR);
	      double SimZErr1 = simTrackZR->Eval(tempR+simTrackZR->GetRMS(1));
	      double SimZErr2 = simTrackZR->Eval(tempR-simTrackZR->GetRMS(1));
	      double SimZ2 = gRandom->Uniform(TMath::Min(SimZErr1,SimZErr2),TMath::Max(SimZErr1,SimZErr2));

	      double Std_SimZerr = TMath::Sqrt(TMath::Sq(SimZ-SimZErr1)+TMath::Sq(SimZ-SimZErr2));
	      LocalHisto.h_ResidualMDC_dZ1->Fill(id_red,(SimZ2-std::get<0>(graphFits[id_red]))/Std_SimZerr);// /TMath::Sqrt(TMath::Sq(std::get<1>(tempPhiRDzdRDPhi))+TMath::Sq(std::get<2>(graphFits[id_red]))));

	      LocalHisto.h_RPhiZMDC_Sigma->Fill(id_red,Std_SimZerr);// /TMath::Sqrt(TMath::Sq(std::get<1>(tempPhiRDzdRDPhi))+TMath::Sq(std::get<2>(graphFits[id_red]))));

	      if(id_detMG.size()>=6)
		LocalHisto.h_ResidualMDC_dZ_More6->Fill(id_red,(SimZ2-std::get<0>(graphFits[id_red]))/Std_SimZerr);
	      if(id_detPSB.size()>0)
		LocalHisto.h_ResidualMDC_dZ_PSB->Fill(id_red,(SimZ2-std::get<0>(graphFits[id_red]))/Std_SimZerr);
	      if(id_detPSBE.size()>0)
		LocalHisto.h_ResidualMDC_dZ_PSBE->Fill(id_red,(SimZ2-std::get<0>(graphFits[id_red]))/Std_SimZerr);
	      if(id_detPSFE.size()>0)
		LocalHisto.h_ResidualMDC_dZ_PSFE->Fill(id_red,(SimZ2-std::get<0>(graphFits[id_red]))/Std_SimZerr);
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
	      LocalHisto.h_ResidualMDC_dZ2->Fill(id_red,(SimZ-std::get<0>(graphFits[id_red]))/TMath::Sqrt(TMath::Sq(std::get<1>(tempPhiRDzdRDPhi))+TMath::Sq(std::get<2>(graphFits[id_red]))));
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
	  auto* listG2 = mg_trackRZ->GetListOfGraphs();
	  listG2->Delete();
	}
      else
	{
	  auto* listG1 = mg_trackPhiZ->GetListOfGraphs();
	  listG1->Delete();
	  auto* listG2 = mg_trackRZ->GetListOfGraphs();
	  listG2->Delete();
	}

      att._logger->debug("@@@");
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

  int idDisc_1 = MDC_PRZ::f_checkDiscontinuity(graph1, TMath::PiOver2());
  int idDisc_2 = MDC_PRZ::f_checkDiscontinuity(graph2, TMath::PiOver2());

  std::pair<int,int> minmaxPoint1, minmaxPoint2;
  if(idDisc_1==-1)
    minmaxPoint1 = {0, graph1->GetN()-1};
  else
    minmaxPoint1 = {0, idDisc_1 -1};

  if(idDisc_2==-1)
    minmaxPoint2 = {0, graph2->GetN()-1};
  else
    minmaxPoint2 = {0, idDisc_2 -1};


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




    //std::cout<<"In tube :"<<p0<<" "<<p1<<"\n";
    //cov.Print();

    TMatrixD SubCovSqrt = f_sqrtMat2x2(cov);
    SubCovSqrt *= TMath::Sqrt(2 * cov.GetNcols());

    //SubCovSqrt.Print();

    std::vector<double> Rs;
    double mean_R = 0.;
    for(int iC = 0; iC < SubCovSqrt.GetNcols(); ++iC)
      {
	double p0_eP = p0 + SubCovSqrt(0, iC);
	double p1_eP = p1 + SubCovSqrt(1, iC);

	double p0_eN = p0 - SubCovSqrt(0, iC);
	double p1_eN = p1 - SubCovSqrt(1, iC);

	//std::cout<<" iC "<<iC<<" +:"<<p0_eP<<" "<<p1_eP<<" -:"<<p0_eN<<" "<<p1_eN<<"\n";

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
	//std::cout<<" Rs :"<<Rsig<<"\n";
      }
    Rstd *= 1./static_cast<double>(Rs.size());

    return std::make_tuple(mean_R,TMath::Sqrt(Rstd));
  };

  //std::cout<<"intersetTrack : graph1RZ :"<< graph1RZ->GetName() << "\n";
  //graph1RZ->Print();

  int res, iP1 = 0;
  double z_intersect, r_intersect;
  for(int iP = 1; iP < graph1RZ->GetN(); ++iP)
    {
      double z1 = graph1RZ->GetPointX(iP-1);
      double R1 = graph1RZ->GetPointY(iP-1);
      double z2 = graph1RZ->GetPointX(iP);
      double R2 = graph1RZ->GetPointY(iP);

      // att._logger->debug("graph1 : {}",graph1->GetName());
      // att._logger->debug("graph2 : {}", graph2->GetName());
      // graph1->Print();
      // att._logger->debug("coordinates g1 : ({}, {}) ; ({}, {}) | g2: ({}, {}) ; ({}, {}) ",z1,phi1,z2,phi2,z3,phi3,z4,phi4);
      //std::cout<<" coordinates g1RZ :("<<z1<<", "<<R1<<") ; ("<<z2<<", "<<R2<<") track:"<<min_Z<<", "<<TrackR_minZ<<") ; ("<<max_Z<<", "<<TrackR_maxZ<<") \n";
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

  //if(res==0)
  //   std::cout<<"E> not interset between Pol1 and RZ wire !\n";
  //else
  //std::cout<<"Find Intersect : graphMG: "<<graph1RZ->GetName()<<" "<<z_intersect<<" "<<r_intersect<<"\n";
  
  auto [Rm1, R1Er] = f_TrackTube(p0,p1,covP01,min_Z);
  auto [Rm2, R2Er] = f_TrackTube(p0,p1,covP01,max_Z);

  //std::cout<<"!> tubeTrack "<<Rm1<<" "<<R1Er<<" "<<Rm2<<" "<<R2Er<<"\n";

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
      //std::cout<<"loop intersect iR:"<<iR<<" iP:"<<iP<<" "<<res<<" "<<tempZ<<" "<<tempR<<" < z < "<<allZ[2]<<" "<<allZ[3]<<"\n";

      // if(TMath::Abs(z_intersect-tempZ)> 2 || TMath::Abs(r_intersect-tempR)> 2)
      // 	std::cout<<"E> interset and mean/error calculation different :"<< iR<<" | "<<z_intersect<<" "<<tempZ<<" | "<<r_intersect<<" "<<tempR<<" | "<<
      // 	  allZ[0]<<" "<<allZ[1]<<" "<<allZ[2]<<" "<<allZ[3]<<" / "<<TempAllR[0]<<" "<<TempAllR[1]<<" "<<TempAllR[2]<<" "<<TempAllR[3]<<"\n";

      if(res == 0)
	{
	  resFinal+=10;
	  //continue;
	}
      else
	{
	  if(TMath::Abs(z_intersect-tempZ) > 10.)
	    {
	      ++resFinal;
	      continue;
	    }

	  ResZR.push_back({tempZ,tempR});
	  resFinal = -1;
	  //break;
	}


      // if(resFinal != -1)
      // 	std::cout<<"E> segment RZ intersection with errors failed ! ["<<iR<<"] "<<graph1RZ->GetName()<<" "<<resFinal<<" "<<TMath::Max(0, iP1-10)<<" "<< TMath::Min(graph1RZ->GetN(),iP1+10)<<"\n";

    }
  //std::cout<<" !> out :"<<ResZR.size()<<"\n";

  TVectorD meanAll(2);
  TMatrixD covAll(2,2);

  MDC_PRZ::f_MeanCov(ResZR,meanAll,covAll);

  // if(TMath::Abs(z_intersect-meanAll(0))>0.5 || TMath::Abs(r_intersect-meanAll(1))>0.5)
  //   std::cout<<"E> interset and mean/error calculation different :"<<z_intersect<<" "<<meanAll(0)<<" | "<<r_intersect<<" "<<meanAll(1)<<"\n";

  return std::make_tuple(1, z_intersect,r_intersect, meanAll(0), meanAll(1), TMath::Sqrt(covAll(0,0)), TMath::Sqrt(covAll(1,1)));

};
