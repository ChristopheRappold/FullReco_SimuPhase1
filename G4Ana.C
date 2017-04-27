#include "./src/Ana_Event/MCAnaEventG4Sol.hh"

#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TDatabasePDG.h"

#define ROOT6
#ifdef ROOT6
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#endif

#include "TVector3.h"

struct DataTreeMC
{
  int MC_id;
  int Mother_id;
  int pdg;
  TLorentzVector MomMass;
};




class Timer
{
  public:
  Timer() : beg_(clock_::now()) {}
  void reset() { beg_ = clock_::now(); }
  double elapsed() const { return std::chrono::duration_cast<second_>(clock_::now() - beg_).count(); }

  private:
  typedef std::chrono::high_resolution_clock clock_;
  typedef std::chrono::duration<double, std::ratio<1> > second_;
  std::chrono::time_point<clock_> beg_;
};

void G4Ana(const std::set<std::string>& ParticleList = std::set<std::string>(), const std::string& nameList = "./rMinBias_r1.0.root",
           const std::string& outfile = "")
{

  TDatabasePDG::Instance()->AddParticle("deuteron", "deuteron", 1.875613, kTRUE, 0., 1. * 3., "Ions", 10000);
  TDatabasePDG::Instance()->AddParticle("triton", "triton", 2.80925, kTRUE, 0., 1. * 3., "Ions", 10001);
  TDatabasePDG::Instance()->AddParticle("alpha", "alpha", 3.727417, kTRUE, 0., 2. * 3., "Ions", 10002);
  TDatabasePDG::Instance()->AddParticle("He3", "He3", 2.80923, kTRUE, 0., 2. * 3., "Ions", 10003);
  TDatabasePDG::Instance()->AddParticle("Li6", "Li6", 5.6015194, kTRUE, 0., 3. * 3., "Ions", 10004);
  TDatabasePDG::Instance()->AddParticle("H3L", "H3L", 2.99114, kFALSE, 0., 1. * 3., "Ions", 20001);
  TDatabasePDG::Instance()->AddParticle("H4L", "H4L", 3.9225, kFALSE, 0., 1. * 3., "Ions", 20002);
  TDatabasePDG::Instance()->AddParticle("He5L", "He5L", 4.8399, kFALSE, 0., 1. * 3., "Ions", 20003);
  
  auto ploting = [](TH1F& h, const std::string& nameC, std::vector<double> range = {10, -10}, std::string option = "") {
    TCanvas* c = new TCanvas(nameC.c_str(), nameC.c_str(), 500, 500);
    c->cd();
    if(range[0] < range[1])
      h.GetXaxis()->SetRangeUser(range[0], range[1]);
    if(option.empty() == false)
      h.Draw(option.c_str());
    else
      h.Draw();
    c->Draw();
  };

  auto plotingArray = [](std::unordered_map<int, TH2F*>& h2, const std::string& nameC, const std::string& nameOpt = "colz") {
    // int index = 0;
    for(auto& hist2d : h2)
      {
        if(hist2d.second->GetEntries() <= 0)
          continue;

        std::string nameCtemp(nameC);
        auto PDG_particle = TDatabasePDG::Instance()->GetParticle(hist2d.first);
        nameCtemp += PDG_particle->GetName(); // std::tostring(index);

        TCanvas* c = new TCanvas(nameCtemp.c_str(), nameCtemp.c_str(), 500, 500);
        c->cd();
        hist2d.second->Draw(nameOpt.c_str());
        c->Draw();
        //++index;
      }
  };

  auto plotingAcceptance = [&ParticleList](std::unordered_map<int, std::vector<TH1F*> >& h1, const std::string& nameC, bool AllIn = false) {
    // int index = 0;
    for(auto& hist1d : h1)
      {
	auto PDG_particle = TDatabasePDG::Instance()->GetParticle(hist1d.first);
	if(std::fabs(PDG_particle->Charge())<1e-1)
	  continue;
	std::string namePar(PDG_particle->GetName());
	
	if(ParticleList.size()>0)
	  {
	    auto it_findPar = ParticleList.find(namePar);
	    if(it_findPar == ParticleList.end())
	      continue;
	  }
	if(hist1d.second[0]->GetEntries()<=0)
	  continue;
	
	std::string nameCtemp(nameC);

	nameCtemp+= namePar;//std::tostring(index);
	for(auto hist : hist1d.second)
	  hist->Sumw2();

	//std::string nameDet [] = {"_CDH","_RPC","_FMF2","_AllDet"};
	
	TCanvas* c = new TCanvas(nameCtemp.c_str(),nameCtemp.c_str(),500,500);
	if(AllIn==false)
	  {
	    c->Divide(2,2);
	    // c->cd(1);
	    // hist1d.second[0]->Draw("e1");
	    // c->cd(2);
	    // hist1d.second[1]->Draw("e1");
	    // c->cd(3);
	    for(size_t i=1; i<hist1d.second.size();++i)
	      {
		c->cd(i);
		std::string nameAcc(hist1d.second[i]->GetName());
		nameAcc+="GeoAcceptance";
		//TH1F* htemp1 = new TH1F(nameAcc.c_str(),nameAcc.c_str(),hist1d.second[0]->GetNbinsX(),hist1d.second[0]->GetXaxis()->GetXmin(),hist1d.second[0]->GetXaxis()->GetXmax());
		for(int n = 1; n<=hist1d.second[0]->GetNbinsX();++n)
		  {
		    Int_t Nall = hist1d.second[0]->GetBinContent(n);
		    Int_t Nacc = hist1d.second[i]->GetBinContent(n);
		    if(Nall<=0)
		      {
			Nacc=0;
			hist1d.second[i]->SetBinContent(n,Nacc);
		      }
		    if(Nall<Nacc)
		      {
			Nacc = Nall;
			hist1d.second[i]->SetBinContent(n,Nacc);
		      }
		  }
		if(hist1d.second[i]->GetEntries()>1)
		  {
		    TH1F* h_acc = (TH1F*) hist1d.second[i]->Clone();
		    h_acc->GetXaxis()->SetRangeUser(hist1d.second[0]->GetXaxis()->GetXmin(),hist1d.second[0]->GetXaxis()->GetXmax());
		    //TGraphAsymmErrors* g_acc = new TGraphAsymmErrors(hist1d.second[0],hist1d.second[1]);
		    h_acc->SetNameTitle(nameAcc.c_str(),nameAcc.c_str());
		    h_acc->Divide(hist1d.second[0]);
		    h_acc->Draw("e1");
		  }
	      }
	  }
	else
	  {
	    c->cd();
	    for(size_t i=hist1d.second.size()-1; i>=1;--i)
	      {
		std::string nameAcc(hist1d.second[i]->GetName());
		nameAcc+="GeoAcceptance";
		//TH1F* htemp1 = new TH1F(nameAcc.c_str(),nameAcc.c_str(),hist1d.second[0]->GetNbinsX(),hist1d.second[0]->GetXaxis()->GetXmin(),hist1d.second[0]->GetXaxis()->GetXmax());
		for(int n = 1; n<=hist1d.second[0]->GetNbinsX();++n)
		  {
		    Int_t Nall = hist1d.second[0]->GetBinContent(n);
		    Int_t Nacc = hist1d.second[i]->GetBinContent(n);
		    if(Nall<=0)
		      {
			Nacc=0;
			hist1d.second[i]->SetBinContent(n,Nacc);
		      }
		    if(Nall<Nacc)
		      {
			Nacc = Nall;
			hist1d.second[i]->SetBinContent(n,Nacc);
		      }
		  }
		if(hist1d.second[i]->GetEntries()>1)
		  {
		    TH1F* h_acc = (TH1F*) hist1d.second[i]->Clone();
		    h_acc->GetXaxis()->SetRangeUser(hist1d.second[0]->GetXaxis()->GetXmin(),hist1d.second[0]->GetXaxis()->GetXmax());
		    //TGraphAsymmErrors* g_acc = new TGraphAsymmErrors(hist1d.second[0],hist1d.second[1]);
		    h_acc->SetNameTitle(nameAcc.c_str(),nameAcc.c_str());
		    h_acc->Divide(hist1d.second[0]);
		    if(i==hist1d.second.size()-1)
		      h_acc->Draw("e1");
		    else
		      {
			h_acc->SetLineColor(i);
			h_acc->Draw("same e1");
		      }
		  }
	      }

	  }
	c->Draw();
	//++index;
      }
  };

  
  auto ploting2D = [](TH2F& h, const std::string& nameC)
    {
      TCanvas* c = new TCanvas(nameC.c_str(),nameC.c_str(),500,500);
      c->cd();
      h.Draw("colz");
      c->Draw();
    };
  
  //TFile* f = new TFile(nameList.c_str());
  TChain* Chain = new TChain("T");
  std::cout << "Files :" << nameList << std::endl;
  if(nameList.find(".root") != std::string::npos)
    {
      std::cout<<"Load from single file "<<nameList<<std::endl;
      int temp_nb = Chain->AddFile(nameList.c_str());
      std::cout<<" Loaded "<<temp_nb<<" files "<<std::endl;
    }
  else
    {
      std::cout << "Adding Chain from List files" << std::endl;
      std::ifstream List(nameList.c_str());
      std::string infiles;
      int nb_file = 0;
      while(std::getline(List,infiles))
	{
	  std::cout<<infiles<<std::endl;
	  int temp_nb = Chain->AddFile(infiles.c_str());
          nb_file += temp_nb;
        }
      std::cout << " Loaded " << nb_file << " files " << std::endl;
    }

  MCAnaEventG4Sol* event = new MCAnaEventG4Sol();

  Chain->SetBranchAddress("MCAnaEventG4Sol", &event);

  // TH1F* tofP = new TH1F("tofplus","tofplus",32,0,32);
  // TH2F* tofP_bar = new TH2F("tofplus_bar","tofplus_bar",32,0,32,10,0,10);
  // std::unordered_map<std::string,TH2F*> tofP_barMomDiff;
  // std::unordered_map<std::string,TH2F*> tofP_barMomParticle;
  // std::unordered_map<std::string,TH2F*> tofP_barTR1x;
  // std::unordered_map<std::string,TH2F*> tofP_barTR2x;
  // std::unordered_map<std::string,TH2F*> tofP_barTR2xMom;
  // std::unordered_map<std::string,TH2F*> tofP_barTR2xTR1x;

  // std::vector< std::vector<std::string> > nameParticleInit = {
  // {"proton","proton"},{"deuteron","deuteron"},{"triton","triton"},{"kaon+","kaonP"},{"pi+","piP"},
  // 							       {"alpha","alpha"},{"He3","He3"},{"e+","electron"},{"pi-","piN"},{"He6[0.0]","He6"},{"Li6[0.0]","Li6"}
  // };

  // std::vector< std::vector<std::string> > nameParticle;
  // if(nameSelected.empty())
  //   nameParticle = nameParticleInit;
  // else
  //   {
  //     for(auto& nameTemp : nameParticleInit)
  // 	if(nameTemp[0]==nameSelected)
  // 	  nameParticle.emplace_back(nameTemp);
  //   }
  // assert(nameParticle.size()!=0);

  // for(auto& Particle : nameParticle)
  //   {
  //     std::string nameH("h_BarMom_");
  //     nameH+=Particle[1];
  //     std::vector<double> range = {0.,10.};
  //     if(Particle[0] == "alpha")
  // 	{range[0] = 8.; range[1] = 16.;}
  //     else if(Particle[0] == "He3")
  // 	{range[0] = 5.; range[1] = 15.;}
  //     else if(Particle[1] == "He6" || Particle[1] == "Li6")
  // 	{range[0] = 14.; range[1] = 24.;}

  //     tofP_barMomParticle.insert({Particle[0],new TH2F(nameH.c_str(),nameH.c_str(),32,0,32,200,range[0],range[1])});

  //     //tofP_barMomDiff.insert({key,new TH2F(nameHD.c_str(),nameHD.c_str(),32,0,32,100,-5,5)});
  //     for(int i=0;i<32;++i)
  // 	{
  // 	  std::string nameHD("h_BarMomDiff_");
  // 	  std::string key(Particle[0]);
  // 	  key+=std::to_string(i);
  // 	  nameHD+=Particle[1];
  // 	  nameHD+=std::to_string(i);
  // 	  tofP_barMomDiff.insert({key,new TH2F(nameHD.c_str(),nameHD.c_str(),100,range[0],range[1],1000,-1,0)});
  // 	}

  //     std::string nameH2("h_BarTR1x_");
  //     nameH2+=Particle[1];
  //     tofP_barTR1x.insert({Particle[0],new TH2F(nameH2.c_str(),nameH2.c_str(),32,0,32,250,0,250)});

  //     std::string nameH3("h_BarTR2x_");
  //     nameH3+=Particle[1];
  //     tofP_barTR2x.insert({Particle[0],new TH2F(nameH3.c_str(),nameH3.c_str(),32,0,32,420,0,420)});

  //     for(int i=0;i<32;++i)
  // 	{
  // 	  std::string nameH4("h_TR2xMom_");
  // 	  std::string key(Particle[0]);
  // 	  key+=std::to_string(i);
  // 	  nameH4+=Particle[1];
  // 	  nameH4+=std::to_string(i);
  // 	  tofP_barTR2xMom.insert({key,new TH2F(nameH4.c_str(),nameH4.c_str(),420,0,420,200,-0.2,0.2)});//range[0],range[1])});
  // 	}

  //     for(int i=0;i<32;++i)
  // 	{
  // 	  std::string nameH5("h_TR2xTR1x_");
  // 	  std::string key(Particle[0]);
  // 	  key+=std::to_string(i);
  // 	  nameH5+=Particle[1];
  // 	  nameH5+=std::to_string(i);
  // 	  tofP_barTR2xTR1x.insert({key,new TH2F(nameH5.c_str(),nameH5.c_str(),420,0,420,250,0,250)});//range[0],range[1])});
  // 	}

  //     TableFinding.insert({Particle[0],std::vector< std::vector< std::set<int> > >(32,std::vector< std::set<int>
  //     >(420,std::set<int>()))});
  //   }

  TH2F* h_Acceptance = new TH2F("Acceptance","Acceptance",20,0,20,10,0,10);

  std::unordered_map<int, std::vector<TH1F*> > h_MomAcc; 
  std::unordered_map<int, std::vector<TH1F*> > h_MomAccDecay; 
  std::unordered_map<int, std::vector<TH1F*> > h_AngleAcc; 
  std::unordered_map<int, std::vector<TH1F*> > h_AngleAccDecay; 
  std::unordered_map<int, std::vector<TH1F*> > h_RapidityAcc; 

  std::unordered_map<int, TH2F*> h_MomAccReco; 

  std::unordered_map<std::string, std::vector<double> > ParticleBinMM = {{"pi-",{0.,2.}},{"pi+",{0.,2.}},{"K+",{0.,3.}},{"K-",{0.,3.}},{"proton",{0.,5.}},
									{"neutron",{0.,5.}},{"He3",{5.,10.}},{"deuteron",{4.,8}},{"triton",{7.,12.}},
									{"alpha",{9.,12.}}};

  std::unordered_map<std::string, std::vector<double> > ParticleBinAA = {{"pi-",{0.,90.}},{"pi+",{0.,90.}},{"K+",{0.,90.}},{"K-",{0.,90.}},{"proton",{0.,90.}},
									{"neutron",{0.,90.}},{"He3",{0.,10.}},{"deuteron",{0.,10.}},{"triton",{0.,10.}},
									{"alpha",{0.,10.}}};

  std::unordered_map<std::string, std::vector<double> > ParticleBinYY = {{"pi-",{0.,5.}},{"pi+",{0.,5.}},{"K+",{0.,5.}},{"K-",{0.,5.}},{"proton",{0.,5.}},
									{"neutron",{0.,5.}},{"He3",{0.,5.}},{"deuteron",{0.,5.}},{"triton",{0.,5.}},
									{"alpha",{0.,5.}}};

  
  
  auto f_createAccHist = [] (std::unordered_map<int, std::vector<TH1F*> >& h_Acc,const std::unordered_map<std::string, std::vector<double> >& ParBin, int  pdg, std::string nameH) -> TH1F*
    {
      auto it_Acc = h_Acc.find(pdg);
      if(it_Acc==h_Acc.end())
	{
	  auto PDG_particle = TDatabasePDG::Instance()->GetParticle(pdg);
	  std::string namePar(PDG_particle->GetName());
		  
	  std::string nameParHist(PDG_particle->GetName());
	  std::replace(nameParHist.begin(),nameParHist.end(),'+','P');
	  std::replace(nameParHist.begin(),nameParHist.end(),'-','N');
	  std::string nameAll(nameParHist);
	  nameAll+=nameH;
	  nameAll+="_All";
	  std::string nameSel(nameParHist);
	  nameSel+=nameH;
	  nameSel+="_Acc";
	  std::string nameSel2 [7] = {nameSel+"CDH",nameSel+"RPC",nameSel+"FMF2",nameSel+"PSCE",nameSel+"PSBE",nameSel+"PSFE",nameSel+"Det"};
	  std::vector<TH1F*> tempHist(8,nullptr);
	  double binMin = 0.;
	  double binMax = 10.;
	  auto it_binMinMax = ParBin.find(namePar);
	  if(it_binMinMax != ParBin.end())
	    {
	      binMin = it_binMinMax->second[0];
	      binMax = it_binMinMax->second[1];
	    }
	  tempHist[0] = new TH1F(nameAll.c_str(),nameAll.c_str(),20,binMin,binMax);
	  tempHist[1] = new TH1F(nameSel2[0].c_str(),nameSel2[0].c_str(),20,binMin,binMax);
	  tempHist[2] = new TH1F(nameSel2[1].c_str(),nameSel2[1].c_str(),20,binMin,binMax);
	  tempHist[3] = new TH1F(nameSel2[2].c_str(),nameSel2[2].c_str(),20,binMin,binMax);
	  tempHist[4] = new TH1F(nameSel2[3].c_str(),nameSel2[3].c_str(),20,binMin,binMax);
	  tempHist[5] = new TH1F(nameSel2[4].c_str(),nameSel2[4].c_str(),20,binMin,binMax);
	  tempHist[6] = new TH1F(nameSel2[5].c_str(),nameSel2[5].c_str(),20,binMin,binMax);
	  tempHist[7] = new TH1F(nameSel2[6].c_str(),nameSel2[6].c_str(),20,binMin,binMax);
	  h_Acc.insert(std::make_pair(pdg,tempHist));
	  return tempHist[0];
	}
      else
	return it_Acc->second[0];
    };
    
  
  
  
  const auto Entries = Chain->GetEntries();
  std::cout << " Entries :" << Entries << std::endl;
  int timing = 0;
  // reader.SetEntriesRange(0,10);
  for(Long64_t iEntry = 0; iEntry < Entries; ++iEntry)
    {
      // cout<<" #"<<iEntry<<" "<<Entries<<" "<<iEntry/Entries<<endl;
      if(static_cast<int>(static_cast<double>(iEntry) / static_cast<double>(Entries) * 10) == timing)
        {
          std::cout << "Processing :" << timing * 10 << "% \n";
          ++timing;
        }
      Chain->GetEntry(iEntry);

      std::vector<DataTreeMC> TempData(event->Nmc);
      
      for(Int_t id = 0; id< event->fMC_Particle->GetEntries();++id)
	{
	  const TMcParticle& MCpar = *(dynamic_cast<TMcParticle*>(event->fMC_Particle->At(id)));
	  TempData[id].MC_id = MCpar.Mc_id;
	  TempData[id].Mother_id = MCpar.Mother_id;
	  TempData[id].pdg = MCpar.Pdg;
	  TempData[id].MomMass = MCpar.MomMass;
	  auto PDG_particle = TDatabasePDG::Instance()->GetParticle(MCpar.Pdg);
	  h_Acceptance->Fill(PDG_particle->GetName(),"All",1);

	  TH1F* hMomAll = f_createAccHist(h_MomAcc,ParticleBinMM,MCpar.Pdg,"Mom");
	  hMomAll->Fill(MCpar.MomMass.P());
	  TH1F* hAngleAll = f_createAccHist(h_AngleAcc,ParticleBinAA,MCpar.Pdg,"Theta");
	  hAngleAll->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
	  TH1F* hYAll = f_createAccHist(h_RapidityAcc,ParticleBinYY,MCpar.Pdg,"Y");
	  hYAll->Fill(MCpar.MomMass.Rapidity());
	  if(MCpar.Mother_id>=0)
	    {
	      h_Acceptance->Fill(PDG_particle->GetName(),"Decay",1);
	      TH1F* hMomAll2 = f_createAccHist(h_MomAccDecay,ParticleBinMM,MCpar.Pdg,"Decay_Mom");
	      hMomAll2->Fill(MCpar.MomMass.P());
	      TH1F* hAngleAll2 = f_createAccHist(h_AngleAccDecay,ParticleBinAA,MCpar.Pdg,"Decay_Theta");
	      hAngleAll2->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
	    }
	}



      std::unordered_map<int,int> idDet;
      for(Int_t id = 0; id< event->FMF2->GetEntries();++id)
	{
	  const TMcHit& MChit = *(dynamic_cast<TMcHit*>(event->FMF2->At(id)));
	  int TrackID = MChit.MC_id;
	  for(auto MCpar : TempData)
	    if(TrackID == MCpar.MC_id)
	      {
		auto it_det = idDet.find(TrackID);
		if(it_det == idDet.end())
		  {
		    auto PDG_particle = TDatabasePDG::Instance()->GetParticle(MChit.Pdg);
		    h_Acceptance->Fill(PDG_particle->GetName(),"FMF2",1);
		    //h_Acceptance->Fill(PDG_particle->GetName(),"Det",1);
		    idDet.insert(std::make_pair(TrackID, MChit.LayerID));
		    auto it_momAcc = h_MomAcc.find(MChit.Pdg);
		    it_momAcc->second[3]->Fill(MCpar.MomMass.P());
		    //it_momAcc->second[7]->Fill(MCpar.MomMass.P());
		    auto it_angleAcc = h_AngleAcc.find(MChit.Pdg);
		    it_angleAcc->second[3]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
		    //it_angleAcc->second[7]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
		    auto it_YAcc = h_RapidityAcc.find(MChit.Pdg);
		    it_YAcc->second[3]->Fill(MCpar.MomMass.Rapidity());
		    //it_YAcc->second[7]->Fill(MCpar.MomMass.Rapidity());
		    if(MCpar.Mother_id>=0)
		      {
			//h_Acceptance->Fill(PDG_particle->GetName(),"Decay_Det",1);
			auto it_momAcc2 = h_MomAccDecay.find(MChit.Pdg);
			it_momAcc2->second[3]->Fill(MCpar.MomMass.P());
			//it_momAcc2->second[7]->Fill(MCpar.MomMass.P());
			auto it_angleAcc2 = h_AngleAccDecay.find(MChit.Pdg);
			it_angleAcc2->second[3]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
			//it_angleAcc2->second[7]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
		      }
		  }
	      }
	}
	
      
      for(Int_t id = 0; id < event->RPC->GetEntries(); ++id)
        {
          const TMcHit& MChit = *(dynamic_cast<TMcHit*>(event->RPC->At(id)));
          int TrackID = MChit.MC_id;
          for(auto& MCpar : TempData)
            if(TrackID == MCpar.MC_id)
              {
                MCpar.FinalHit.insert(std::make_pair("RPC", MChit.MCHit));
                auto PDG_particle = TDatabasePDG::Instance()->GetParticle(MChit.Pdg);
                h_Acceptance->Fill(PDG_particle->GetName(), "RPC", 1);
                h_Acceptance->Fill(PDG_particle->GetName(), "Det", 1);
                auto it_momAcc = h_MomAcc.find(MChit.Pdg);
                it_momAcc->second[2]->Fill(MCpar.MomMass.P());
                it_momAcc->second[7]->Fill(MCpar.MomMass.P());
                auto it_angleAcc = h_AngleAcc.find(MChit.Pdg);
                it_angleAcc->second[2]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
                it_angleAcc->second[7]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
                auto it_YAcc = h_RapidityAcc.find(MChit.Pdg);
                it_YAcc->second[2]->Fill(MCpar.MomMass.Rapidity());
                it_YAcc->second[7]->Fill(MCpar.MomMass.Rapidity());
                if(MCpar.Mother_id >= 0)
                  {
                    h_Acceptance->Fill(PDG_particle->GetName(), "Decay_Det", 1);
                    auto it_momAcc2 = h_MomAccDecay.find(MChit.Pdg);
                    it_momAcc2->second[2]->Fill(MCpar.MomMass.P());
                    it_momAcc2->second[7]->Fill(MCpar.MomMass.P());
                    auto it_angleAcc2 = h_AngleAccDecay.find(MChit.Pdg);
                    it_angleAcc2->second[2]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
                    it_angleAcc2->second[7]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
		    if(DaughterFragInFRS)
		      {
			h_Acceptance->Fill(PDG_particle->GetName(), "Decay_DetCoinFRS", 1);
			auto it_momAcc3 = h_MomAccDecayCoinFRS.find(MChit.Pdg);
			it_momAcc3->second[2]->Fill(MCpar.MomMass.P());
			it_momAcc3->second[7]->Fill(MCpar.MomMass.P());
			auto it_angleAcc3 = h_AngleAccDecayCoinFRS.find(MChit.Pdg);
			it_angleAcc3->second[2]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
			it_angleAcc3->second[7]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
		      }
		    
                  }
              }

      for(Int_t id = 0; id< event->PSCE->GetEntries();++id)
	{
	  const TMcHit& MChit = *(dynamic_cast<TMcHit*>(event->PSCE->At(id)));
	  int TrackID = MChit.MC_id;
	  for(auto MCpar : TempData)
	    if(TrackID == MCpar.MC_id)
	      {
		auto PDG_particle = TDatabasePDG::Instance()->GetParticle(MChit.Pdg);
		h_Acceptance->Fill(PDG_particle->GetName(),"PSCE",1);
		h_Acceptance->Fill(PDG_particle->GetName(),"Det",1);
		auto it_momAcc = h_MomAcc.find(MChit.Pdg);
		it_momAcc->second[4]->Fill(MCpar.MomMass.P());
		it_momAcc->second[7]->Fill(MCpar.MomMass.P());
		auto it_angleAcc = h_AngleAcc.find(MChit.Pdg);
		it_angleAcc->second[4]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
		it_angleAcc->second[7]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
		auto it_YAcc = h_RapidityAcc.find(MChit.Pdg);
		it_YAcc->second[4]->Fill(MCpar.MomMass.Rapidity());
		it_YAcc->second[7]->Fill(MCpar.MomMass.Rapidity());
		if(MCpar.Mother_id>=0)
		  {
		    h_Acceptance->Fill(PDG_particle->GetName(),"Decay_Det",1);
		    auto it_momAcc2 = h_MomAccDecay.find(MChit.Pdg);
		    it_momAcc2->second[4]->Fill(MCpar.MomMass.P());
		    it_momAcc2->second[7]->Fill(MCpar.MomMass.P());
		    auto it_angleAcc2 = h_AngleAccDecay.find(MChit.Pdg);
		    it_angleAcc2->second[4]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
		    it_angleAcc2->second[7]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
		  }


	      }
	  
	}
      
      for(Int_t id = 0; id< event->PSBE->GetEntries();++id)
	{
	  const TMcHit& MChit = *(dynamic_cast<TMcHit*>(event->PSBE->At(id)));
	  int TrackID = MChit.MC_id;
	  for(auto MCpar : TempData)
	    if(TrackID == MCpar.MC_id)
	      {
		auto PDG_particle = TDatabasePDG::Instance()->GetParticle(MChit.Pdg);
		h_Acceptance->Fill(PDG_particle->GetName(),"PSBE",1);
		//h_Acceptance->Fill(PDG_particle->GetName(),"Det",1);
		auto it_momAcc = h_MomAcc.find(MChit.Pdg);
		it_momAcc->second[5]->Fill(MCpar.MomMass.P());
		//it_momAcc->second[7]->Fill(MCpar.MomMass.P());
		auto it_angleAcc = h_AngleAcc.find(MChit.Pdg);
		it_angleAcc->second[5]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
		//it_angleAcc->second[7]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
		auto it_YAcc = h_RapidityAcc.find(MChit.Pdg);
		it_YAcc->second[5]->Fill(MCpar.MomMass.Rapidity());
		//it_YAcc->second[7]->Fill(MCpar.MomMass.Rapidity());
		if(MCpar.Mother_id>=0)
		  {
		    //h_Acceptance->Fill(PDG_particle->GetName(),"Decay_Det",1);
		    auto it_momAcc2 = h_MomAccDecay.find(MChit.Pdg);
		    it_momAcc2->second[5]->Fill(MCpar.MomMass.P());
		    //it_momAcc2->second[7]->Fill(MCpar.MomMass.P());
		    auto it_angleAcc2 = h_AngleAccDecay.find(MChit.Pdg);
		    it_angleAcc2->second[5]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
		    //it_angleAcc2->second[7]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
      for(Int_t id = 0; id < event->CDH->GetEntries(); ++id)
        {
          const TMcHit& MChit = *(dynamic_cast<TMcHit*>(event->CDH->At(id)));
          int TrackID = MChit.MC_id;
          for(auto& MCpar : TempData)
            if(TrackID == MCpar.MC_id)
              {
                MCpar.FinalHit.insert(std::make_pair("CDH", MChit.MCHit));

                auto PDG_particle = TDatabasePDG::Instance()->GetParticle(MChit.Pdg);
                h_Acceptance->Fill(PDG_particle->GetName(), "CDH", 1);
                h_Acceptance->Fill(PDG_particle->GetName(), "Det", 1);
                h_Acceptance->Fill(PDG_particle->GetName(), "Det2", 1);
                auto it_momAcc = h_MomAcc.find(MChit.Pdg);
                it_momAcc->second[1]->Fill(MCpar.MomMass.P());
                it_momAcc->second[7]->Fill(MCpar.MomMass.P());
                auto it_angleAcc = h_AngleAcc.find(MChit.Pdg);
                it_angleAcc->second[1]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
                it_angleAcc->second[7]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
                auto it_YAcc = h_RapidityAcc.find(MChit.Pdg);
                it_YAcc->second[1]->Fill(MCpar.MomMass.Rapidity());
                it_YAcc->second[7]->Fill(MCpar.MomMass.Rapidity());
                if(MCpar.Mother_id >= 0)
                  {
                    h_Acceptance->Fill(PDG_particle->GetName(), "Decay_Det", 1);
                    h_Acceptance->Fill(PDG_particle->GetName(), "Decay_Det2", 1);
                    auto it_momAcc2 = h_MomAccDecay.find(MChit.Pdg);
                    it_momAcc2->second[1]->Fill(MCpar.MomMass.P());
                    it_momAcc2->second[7]->Fill(MCpar.MomMass.P());
                    auto it_angleAcc2 = h_AngleAccDecay.find(MChit.Pdg);
                    it_angleAcc2->second[1]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
                    it_angleAcc2->second[7]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
		    if(DaughterFragInFRS)
		      {
			h_Acceptance->Fill(PDG_particle->GetName(), "Decay_DetCoinFRS", 1);
			auto it_momAcc3 = h_MomAccDecayCoinFRS.find(MChit.Pdg);
			it_momAcc3->second[1]->Fill(MCpar.MomMass.P());
			it_momAcc3->second[7]->Fill(MCpar.MomMass.P());
			auto it_angleAcc3 = h_AngleAccDecayCoinFRS.find(MChit.Pdg);
			it_angleAcc3->second[1]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
			it_angleAcc3->second[7]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
		      }
		  }
	      }
	  
	}

      
      for(Int_t id = 0; id< event->PSFE->GetEntries();++id)
	{
	  const TMcHit& MChit = *(dynamic_cast<TMcHit*>(event->PSFE->At(id)));
	  int TrackID = MChit.MC_id;
	  for(auto MCpar : TempData)
	    if(TrackID == MCpar.MC_id)
	      {
		auto PDG_particle = TDatabasePDG::Instance()->GetParticle(MChit.Pdg);
		h_Acceptance->Fill(PDG_particle->GetName(),"PSFE",1);
		//h_Acceptance->Fill(PDG_particle->GetName(),"Det",1);
		auto it_momAcc = h_MomAcc.find(MChit.Pdg);
		it_momAcc->second[6]->Fill(MCpar.MomMass.P());
		//it_momAcc->second[7]->Fill(MCpar.MomMass.P());
		auto it_angleAcc = h_AngleAcc.find(MChit.Pdg);
		it_angleAcc->second[6]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
		//it_angleAcc->second[7]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
		auto it_YAcc = h_RapidityAcc.find(MChit.Pdg);
		it_YAcc->second[6]->Fill(MCpar.MomMass.Rapidity());
		//it_YAcc->second[7]->Fill(MCpar.MomMass.Rapidity());
		if(MCpar.Mother_id>=0)
		  {
		    //h_Acceptance->Fill(PDG_particle->GetName(),"Decay_Det",1);
		    auto it_momAcc2 = h_MomAccDecay.find(MChit.Pdg);
		    it_momAcc2->second[6]->Fill(MCpar.MomMass.P());
		    //it_momAcc2->second[7]->Fill(MCpar.MomMass.P());
		    auto it_angleAcc2 = h_AngleAccDecay.find(MChit.Pdg);
		    it_angleAcc2->second[6]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
		    //it_angleAcc2->second[7]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
		  }
	      }
	  
	}





      
      for(Int_t id = 0; id< event->fTrack->GetEntries();++id)
	{
	  const THyphiTrack& RecoTrack = *(dynamic_cast<THyphiTrack*>(event->fTrack->At(id)));
	  int TrackID = RecoTrack.MC_status;
	  for(auto MCpar :TempData)
	    if(TrackID == MCpar.MC_id)
	      {
		if(RecoTrack.Pval2<0.001)
		  continue;
		auto it_2DHist = h_MomAccReco.find(MCpar.pdg);
		if(it_2DHist == h_MomAccReco.end())
		  {
		    auto PDG_particle = TDatabasePDG::Instance()->GetParticle(MCpar.pdg);
		    std::string namePar(PDG_particle->GetName());
		    std::string nameParHist(namePar);
		    std::replace(nameParHist.begin(),nameParHist.end(),'+','P');
		    std::replace(nameParHist.begin(),nameParHist.end(),'-','N');
		    nameParHist+="_RecoMom";
		    double binMin = 0.;
		    double binMax = 10.;
		    auto it_binMinMax = ParticleBinMM.find(namePar);
		    if(it_binMinMax != ParticleBinMM.end())
		      {
			binMin = it_binMinMax->second[0];
			binMax = it_binMinMax->second[1];
		      }
		    TH2F* h_Reco = new TH2F(nameParHist.c_str(),nameParHist.c_str(),25,binMin,binMax,1000,0,0.1);
		    h_MomAccReco.insert(std::make_pair(MCpar.pdg,h_Reco));
		    h_Reco->Fill(MCpar.MomMass.P(),std::fabs(RecoTrack.MomMass.P()-MCpar.MomMass.P())/MCpar.MomMass.P());
		  }
		else
		  it_2DHist->second->Fill(MCpar.MomMass.P(),std::fabs(RecoTrack.MomMass.P()-MCpar.MomMass.P())/MCpar.MomMass.P());
	      }
	}
      
      //cout<<"\n";
    }

  if(outfile.empty())
    {
      ploting2D(*h_Acceptance,"c1");

      plotingAcceptance(h_MomAcc,"c1_MomAcc",true);

      plotingAcceptance(h_AngleAcc,"c1_AngleAcc",true);

      plotingAcceptance(h_RapidityAcc,"c1_RapidityAcc",true);

      plotingAcceptance(h_MomAccDecay,"c1_MomAccDecay",true);
      
      plotingAcceptance(h_AngleAccDecay,"c1_AngleAccDecay",true);
      
      plotingArray(h_MomAccReco,"c1_Reco","candle");
    }
  else
    {
      TFile* fout = new TFile(outfile.c_str(),"RECREATE");
      fout->cd();
      h_Acceptance->Write();

      auto f_DoEff = [] (std::vector<TH1F*>& hist1d) -> std::vector<TH1F*> {
	
	std::vector<TH1F*> h_out(hist1d.size()-1,nullptr);
	for(size_t i=hist1d.size()-1; i>=1;--i)
	  {
	    std::string nameAcc(hist1d[i]->GetName());
	    nameAcc+="GeoAcceptance";
	    //TH1F* htemp1 = new TH1F(nameAcc.c_str(),nameAcc.c_str(),hist1d.second[0]->GetNbinsX(),hist1d.second[0]->GetXaxis()->GetXmin(),hist1d.second[0]->GetXaxis()->GetXmax());
	    for(int n = 1; n<=hist1d[0]->GetNbinsX();++n)
	      {
		Int_t Nall = hist1d[0]->GetBinContent(n);
		Int_t Nacc = hist1d[i]->GetBinContent(n);
		if(Nall<=0)
		  {
		    Nacc=0;
		    hist1d[i]->SetBinContent(n,Nacc);
		  }
		if(Nall<Nacc)
		  {
			Nacc = Nall;
			hist1d[i]->SetBinContent(n,Nacc);
		  }
	      }
	    if(hist1d[i]->GetEntries()>1)
	      {
		h_out[i-1] = (TH1F*) hist1d[i]->Clone();
		h_out[i-1]->GetXaxis()->SetRangeUser(hist1d[0]->GetXaxis()->GetXmin(),hist1d[0]->GetXaxis()->GetXmax());
		//TGraphAsymmErrors* g_acc = new TGraphAsymmErrors(hist1d.second[0],hist1d.second[1]);
		h_out[i-1]->SetNameTitle(nameAcc.c_str(),nameAcc.c_str());
		h_out[i-1]->Divide(hist1d[0]);
		if(i!=hist1d.size()-1)
		  h_out[i-1]->SetLineColor(i);
	      }
	  }
	return h_out;
     };
      

      for(auto it_momAcc : h_MomAcc)
	{
	  int pdg = it_momAcc.first;
	  auto PDG_particle = TDatabasePDG::Instance()->GetParticle(pdg);
	  std::string nameDir(PDG_particle->GetName());

	  auto it_angleAcc = h_AngleAcc.find(pdg);
	  auto it_rapAcc = h_RapidityAcc.find(pdg);
	  auto it_momRes = h_MomAccReco.find(pdg);

	  std::vector<TH1F*> momEff = f_DoEff(it_momAcc.second);
	  std::vector<TH1F*> angleEff = f_DoEff(it_angleAcc->second);
	  std::vector<TH1F*> rapEff = f_DoEff(it_rapAcc->second);
       	  TDirectory* dirTemp = fout->mkdir(nameDir.c_str(),nameDir.c_str());
       	  dirTemp->cd();
	  for(auto htemp : momEff)
	    if(htemp!=nullptr)
	      htemp->Write();

	  for(auto htemp : angleEff)
	    if(htemp!=nullptr)
	      htemp->Write();

	  for(auto htemp : rapEff)
	    if(htemp!=nullptr)
	      htemp->Write();

	  if(it_momRes!=h_MomAccReco.end())
	    it_momRes->second->Write();
	}
      
      fout->Close();
      // auto writeToFile = [&fout,&nameParticle](std::unordered_map<int,TH2F*>& h2, const std::string& nameDir)
      //  	{
      // 	  TDirectory* dirTemp = fout->mkdir(nameDir.c_str(),nameDir.c_str());
      // 	  dirTemp->cd();
      // 	  std::vector<TDirectory*> subDirTemp (nameParticle.size(),nullptr);

      // 	  for(auto it_nameP = nameParticle.begin(), it_namePend = nameParticle.end();it_nameP!=it_namePend;++it_nameP)
      // 	    {
      // 	      size_t index = std::distance(nameParticle.begin(),it_nameP);
      // 	      subDirTemp[index] = dirTemp->mkdir((*it_nameP)[0].c_str());
      // 	    }

      // 	  for(auto& hist2d : h2)
      // 	    {
      // 	      if(hist2d.second->GetEntries()<=0)
      // 		continue;

      // 	      for(auto it_Dir : subDirTemp)
      // 		{
      // 		  std::string nameSubDir(it_Dir->GetName());
      // 		  std::string::size_type found = hist2d.first.find(nameSubDir);
      // 		  if(found == std::string::npos)
      // 		    continue;
      // 		  it_Dir->cd();
      // 		  hist2d.second->SetDrawOption("colz");
      // 		  hist2d.second->Write();
      // 		}
      // 	    }
      // 	};

      // tofP->Write();
      // tofP_bar->Write();
      // writeToFile(tofP_barTR2xMom,"Tr2Bar_PxPz");
      // writeToFile(tofP_barTR1x,"Tr1");
      // writeToFile(tofP_barTR2x,"Tr2");
      // writeToFile(tofP_barMomParticle,"TOFP_Mom");
      // writeToFile(tofP_barMomDiff,"TOFP_MomDiff");
      // writeToFile(tofP_barTR2xTR1x,"TOFP_Tr2Tr1");

      // fout->Close();
    }

}

