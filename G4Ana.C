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
  std::unordered_map<std::string, TVector3> FinalHit;
};

struct ParticleD
{
  TLorentzVector MomMass;
  TLorentzVector Vtx;
  int pdg;
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
  auto AddToDatabasePDG = [](TObjArray* array) {
			  for(auto temp: *array)
			    {
			      TParticlePDG* temp1 = dynamic_cast<TParticlePDG*>(temp);
			      TParticlePDG* old = TDatabasePDG::Instance()->GetParticle(temp1->PdgCode());
			      if(old == nullptr)
				TDatabasePDG::Instance()->AddParticle(temp1->GetName(),temp1->GetTitle(), temp1->Mass(), temp1->Stable(), temp1->Width(), temp1->Charge(), temp1->ParticleClass(), temp1->PdgCode());
			    }
			};

  
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
        if(std::fabs(PDG_particle->Charge()) < 1e-1)
          continue;
        std::string namePar(PDG_particle->GetName());

        if(ParticleList.size() > 0)
          {
            auto it_findPar = ParticleList.find(namePar);
            if(it_findPar == ParticleList.end())
              continue;
          }
        if(hist1d.second[0]->GetEntries() <= 0)
          continue;

        std::string nameCtemp(nameC);

        nameCtemp += namePar; // std::tostring(index);
        for(auto hist : hist1d.second)
          hist->Sumw2();

        // std::string nameDet [] = {"_CDH","_RPC","_FMF2","_AllDet"};

        TCanvas* c = new TCanvas(nameCtemp.c_str(), nameCtemp.c_str(), 500, 500);
        if(AllIn == false)
          {
            c->Divide(2, 2);
            // c->cd(1);
            // hist1d.second[0]->Draw("e1");
            // c->cd(2);
            // hist1d.second[1]->Draw("e1");
            // c->cd(3);
            for(size_t i = 1; i < hist1d.second.size(); ++i)
              {
                c->cd(i);
                std::string nameAcc(hist1d.second[i]->GetName());
                nameAcc += "GeoAcceptance";
                // TH1F* htemp1 = new
                // TH1F(nameAcc.c_str(),nameAcc.c_str(),hist1d.second[0]->GetNbinsX(),hist1d.second[0]->GetXaxis()->GetXmin(),hist1d.second[0]->GetXaxis()->GetXmax());
                for(int n = 1; n <= hist1d.second[0]->GetNbinsX(); ++n)
                  {
                    Int_t Nall = hist1d.second[0]->GetBinContent(n);
                    Int_t Nacc = hist1d.second[i]->GetBinContent(n);
                    if(Nall <= 0)
                      {
                        Nacc = 0;
                        hist1d.second[i]->SetBinContent(n, Nacc);
                      }
                    if(Nall < Nacc)
                      {
                        Nacc = Nall;
                        hist1d.second[i]->SetBinContent(n, Nacc);
                      }
                  }
                if(hist1d.second[i]->GetEntries() > 1)
                  {
                    TH1F* h_acc = (TH1F*)hist1d.second[i]->Clone();
                    h_acc->GetXaxis()->SetRangeUser(hist1d.second[0]->GetXaxis()->GetXmin(), hist1d.second[0]->GetXaxis()->GetXmax());
                    // TGraphAsymmErrors* g_acc = new TGraphAsymmErrors(hist1d.second[0],hist1d.second[1]);
                    h_acc->SetNameTitle(nameAcc.c_str(), nameAcc.c_str());
                    h_acc->Divide(hist1d.second[0]);
                    h_acc->Draw("e1");
                  }
              }
          }
        else
          {
            c->cd();
            for(size_t i = hist1d.second.size() - 1; i >= 1; --i)
              {
                std::string nameAcc(hist1d.second[i]->GetName());
                nameAcc += "GeoAcceptance";
                // TH1F* htemp1 = new
                // TH1F(nameAcc.c_str(),nameAcc.c_str(),hist1d.second[0]->GetNbinsX(),hist1d.second[0]->GetXaxis()->GetXmin(),hist1d.second[0]->GetXaxis()->GetXmax());
                for(int n = 1; n <= hist1d.second[0]->GetNbinsX(); ++n)
                  {
                    Int_t Nall = hist1d.second[0]->GetBinContent(n);
                    Int_t Nacc = hist1d.second[i]->GetBinContent(n);
                    if(Nall <= 0)
                      {
                        Nacc = 0;
                        hist1d.second[i]->SetBinContent(n, Nacc);
                      }
                    if(Nall < Nacc)
                      {
                        Nacc = Nall;
                        hist1d.second[i]->SetBinContent(n, Nacc);
                      }
                  }
                if(hist1d.second[i]->GetEntries() > 1)
                  {
                    TH1F* h_acc = (TH1F*)hist1d.second[i]->Clone();
                    h_acc->GetXaxis()->SetRangeUser(hist1d.second[0]->GetXaxis()->GetXmin(), hist1d.second[0]->GetXaxis()->GetXmax());
                    // TGraphAsymmErrors* g_acc = new TGraphAsymmErrors(hist1d.second[0],hist1d.second[1]);
                    h_acc->SetNameTitle(nameAcc.c_str(), nameAcc.c_str());
                    h_acc->Divide(hist1d.second[0]);
                    if(i == hist1d.second.size() - 1)
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

  auto ploting2D = [](TH2F& h, const std::string& nameC) {
    TCanvas* c = new TCanvas(nameC.c_str(), nameC.c_str(), 500, 500);
    c->cd();
    h.Draw("colz");
    c->Draw();
  };

  auto plotingCDCFinal = [&ParticleList](std::unordered_map<std::string, std::vector<TH2F*> >& h_CDCFinalHit, const std::string& name) {
    for(auto it_hist : h_CDCFinalHit)
      {
        std::string nameC(name);
	std::cout<<"name FinalHit key:"<<it_hist.first<<"\n";
        // int pdg_code = std::stoi(it_hist.first, nullptr);
        auto it_posPDG = it_hist.first.find_first_of("-1234567890");
        std::string nameDet = it_hist.first.substr(0, it_posPDG);

        std::string namePDG = it_hist.first.substr(it_posPDG, std::string::npos);
        int pdg_code = std::stoi(namePDG);

	auto PDG_particle = TDatabasePDG::Instance()->GetParticle(pdg_code);

	std::string namePar(PDG_particle->GetName());
        if(ParticleList.size() > 0)
          {
            auto it_findPar = ParticleList.find(namePar);
            if(it_findPar == ParticleList.end())
              continue;
          }
	
	
	std::replace(namePar.begin(), namePar.end(), '+', 'P');
        std::replace(namePar.begin(), namePar.end(), '-', 'N');

	nameC += nameDet;
        nameC += namePar;
	auto it_multi = it_hist.first.find("Multiplicity");
	if(it_multi!=std::string::npos)
	  nameC += "Multiplicity";
	auto it_decay = it_hist.first.find("Decay");
	if(it_decay!=std::string::npos)
	  nameC += "Decay";
	auto it_frs = it_hist.first.find("CoinFRS");
	if(it_frs != std::string::npos)
	  nameC += "CoinFRS";
	
	TCanvas* c = new TCanvas(nameC.c_str(), nameC.c_str(), 600, 600);
        c->Divide(2, 2);
        for(size_t id = 0; id < it_hist.second.size(); ++id)
          {
            c->cd(1 + id);
            it_hist.second[id]->DrawNormalized("colz");
          }
        c->Draw();
      }
  };

  // TFile* f = new TFile(nameList.c_str());
  TChain* Chain = new TChain("T");
  std::cout << "Files :" << nameList << std::endl;
  if(nameList.find(".root") != std::string::npos)
    {
      TFile* f_in = new TFile(nameList.c_str());
      TObjArray* addPDG = (TObjArray*) f_in->Get("additionalPDG");
      if(addPDG!=nullptr)
	AddToDatabasePDG(addPDG);
      f_in->Close();
      
      std::cout << "Load from single file " << nameList << std::endl;
      int temp_nb = Chain->AddFile(nameList.c_str());
      std::cout << " Loaded " << temp_nb << " files " << std::endl;

    }
  else
    {
      std::cout << "Adding Chain from List files" << std::endl;
      std::ifstream List(nameList.c_str());
      std::string infiles;
      int nb_file = 0;
      while(std::getline(List, infiles))
        {	  
	  std::cout << infiles << std::endl;
	  TFile* f_in = new TFile(infiles.c_str());
	  TObjArray* addPDG = (TObjArray*) f_in->Get("additionalPDG");
	  if(addPDG!=nullptr)
	    AddToDatabasePDG(addPDG);
	  f_in->Close();

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

  std::unordered_map< std::string,std::vector<double> > binInvMass = { {"H3L",{2.5,3.5}}, {"H4L", {3.5,4.5}}, {"nnL",{1.6,2.6}}, {"Lambda",{0.5,1.5}}};

  
  TH2F* h_Acceptance = new TH2F("Acceptance", "Acceptance", 20, 0, 20, 10, 0, 10);
  TH1F* h_Vtx = new TH1F("Vtx","Vtx",1000,0,10);
  TH1F* h_VtxMix = new TH1F("VtxMix","VtxMix",1000,0,10);
  TH2F* h_RZVtx = new TH2F("RZVtx","RZVtx",200,0,10,200,0,200);

  TH1F* h_Vtx_X = new TH1F("Res_Xvtx","Res_Xvtx",1000,-2,2);
  TH1F* h_Vtx_Y = new TH1F("Res_Yvtx","Res_Yvtx",1000,-1,1);
  TH1F* h_Vtx_Z = new TH1F("Res_Zvtx","Res_Zvtx",1000,-1,1);

  TH1F* h_InvMass = nullptr;
  std::unordered_map<int,TH1F*> h_InvMassMix;
  std::unordered_map<int,TH2F*> h_InvMassBrhoMix;
  TH2F* h_InvMassVtx = nullptr;
  TH2F* h_InvMassVtxMix = nullptr;
  TH2F* h_InvMassBrho = nullptr;
  
  std::vector<TH2F*> h_InvMassBrho_multi; 

  for(auto it_binInv : binInvMass)
    {
      if( nameList.find(it_binInv.first) != std::string::npos)
	{
	  h_InvMass = new TH1F("InvMass","InvMass",1000,it_binInv.second[0],it_binInv.second[1]);
	  h_InvMassVtx = new TH2F("InvMassVtx","InvMassVtx",1000,it_binInv.second[0],it_binInv.second[1],200,0,1);
	  h_InvMassBrho = new TH2F("InvMassBrho","InvMassBrho",1000,it_binInv.second[0],it_binInv.second[1],200,10,20);
	  for(size_t i=0; i<20; ++i)
	    {
	      TString nameTemp("InvMassBrho_MoreMulti_");
	      nameTemp+=i;
	      h_InvMassBrho_multi.emplace_back( new TH2F(nameTemp,nameTemp,1000,it_binInv.second[0],it_binInv.second[1],200,10,20) );      
	    }
	}
    }
  h_InvMassMix.insert(std::make_pair(10003,new TH1F("h_InvMassMix_H3L","h_InvMassMix_H3L",1000, binInvMass["H3L"][0], binInvMass["H3L"][1])));
  h_InvMassMix.insert(std::make_pair(10002,new TH1F("h_InvMassMix_H4L","h_InvMassMix_H4L",1000, binInvMass["H4L"][0], binInvMass["H4L"][1])));
  h_InvMassMix.insert(std::make_pair(10000,new TH1F("h_InvMassMix_nnL","h_InvMassMix_nnL",1000, binInvMass["nnL"][0], binInvMass["nnL"][1])));
  h_InvMassMix.insert(std::make_pair(2212,new TH1F("h_InvMassMix_Lambda","h_InvMassMix_Lambda",1000, binInvMass["Lambda"][0], binInvMass["Lambda"][1])));

  h_InvMassBrhoMix.insert(std::make_pair(10003,new TH2F("h_InvMassBrhoMix_H3L","h_InvMassBrhoMix_H3L",1000, binInvMass["H3L"][0], binInvMass["H3L"][1],200,10,20)));
  h_InvMassBrhoMix.insert(std::make_pair(10002,new TH2F("h_InvMassBrhoMix_H4L","h_InvMassBrhoMix_H4L",1000, binInvMass["H4L"][0], binInvMass["H4L"][1],200,10,20)));
  h_InvMassBrhoMix.insert(std::make_pair(10000,new TH2F("h_InvMassBrhoMix_nnL","h_InvMassBrhoMix_nnL",1000, binInvMass["nnL"][0], binInvMass["nnL"][1],200,10,20)));
  h_InvMassBrhoMix.insert(std::make_pair(2212,new TH2F("h_InvMassBrhoMix_Lambda","h_InvMassBrhoMix_Lambda",1000, binInvMass["Lambda"][0], binInvMass["Lambda"][1],200,10,20)));

  h_InvMassVtxMix = new TH2F("InvMassVtxMix","InvMassVtxMix",3000,1.6,4.6,200,0,10);
  
  std::unordered_map<int, std::vector<TH1F*> > h_MomAcc;
  std::unordered_map<int, std::vector<TH1F*> > h_MomAccDecay;
  std::unordered_map<int, std::vector<TH1F*> > h_MomAccDecayCoinFRS;
  std::unordered_map<int, std::vector<TH1F*> > h_AngleAcc;
  std::unordered_map<int, std::vector<TH1F*> > h_AngleAccDecay;
  std::unordered_map<int, std::vector<TH1F*> > h_AngleAccDecayCoinFRS;
  std::unordered_map<int, std::vector<TH1F*> > h_RapidityAcc;

  std::unordered_map<std::string, std::vector<TH2F*> > h_FinalHit_CDC;

  std::unordered_map<int, TH2F*> h_MomAccReco;
  std::unordered_map<int, std::vector<TH2F*> > h_MomAccReco_multi;

  std::unordered_map<std::string, std::vector<double> > ParticleBinMM = {
      {"pi-", {0., 2.}},     {"pi+", {0., 2.}},  {"K+", {0., 3.}},      {"K-", {0., 3.}},      {"proton", {0., 5.}},
      {"neutron", {0., 5.}}, {"He3", {5., 10.}}, {"deuteron", {4., 8}}, {"triton", {7., 12.}}, {"alpha", {9., 12.}}};

  std::unordered_map<std::string, std::vector<double> > ParticleBinAA = {
      {"pi-", {0., 90.}},     {"pi+", {0., 90.}}, {"K+", {0., 90.}},       {"K-", {0., 90.}},     {"proton", {0., 90.}},
      {"neutron", {0., 90.}}, {"He3", {0., 10.}}, {"deuteron", {0., 10.}}, {"triton", {0., 10.}}, {"alpha", {0., 10.}}};

  std::unordered_map<std::string, std::vector<double> > ParticleBinYY = {
      {"pi-", {0., 5.}},     {"pi+", {0., 5.}}, {"K+", {0., 5.}},       {"K-", {0., 5.}},     {"proton", {0., 5.}},
      {"neutron", {0., 5.}}, {"He3", {0., 5.}}, {"deuteron", {0., 5.}}, {"triton", {0., 5.}}, {"alpha", {0., 5.}}};

  // std::vector<std::string> ParticleBinFF = {"pi-","pi+","K+","K-","proton","neutron","He3","deuteron","triton","alpha"};

  auto f_createAccHist = [](std::unordered_map<int, std::vector<TH1F*> >& h_Acc,
                            const std::unordered_map<std::string, std::vector<double> >& ParBin, int pdg, std::string nameH) -> TH1F* {
    auto it_Acc = h_Acc.find(pdg);
    if(it_Acc == h_Acc.end())
      {
        auto PDG_particle = TDatabasePDG::Instance()->GetParticle(pdg);
        std::string namePar(PDG_particle->GetName());

        std::string nameParHist(PDG_particle->GetName());
        std::replace(nameParHist.begin(), nameParHist.end(), '+', 'P');
        std::replace(nameParHist.begin(), nameParHist.end(), '-', 'N');
        std::string nameAll(nameParHist);
        nameAll += nameH;
        nameAll += "_All";
        std::string nameSel(nameParHist);
        nameSel += nameH;
        nameSel += "_Acc";
        std::string nameSel2[] = {nameSel + "CDH",  nameSel + "RPC",  nameSel + "FMF2", nameSel + "PSCE",
                                   nameSel + "PSBE", nameSel + "PSFE", nameSel + "Det", nameSel + "Det2"};
        std::vector<TH1F*> tempHist(9, nullptr);
        double binMin = 0.;
        double binMax = 10.;
        auto it_binMinMax = ParBin.find(namePar);
        if(it_binMinMax != ParBin.end())
          {
            binMin = it_binMinMax->second[0];
            binMax = it_binMinMax->second[1];
          }
        tempHist[0] = new TH1F(nameAll.c_str(), nameAll.c_str(), 20, binMin, binMax);
        tempHist[1] = new TH1F(nameSel2[0].c_str(), nameSel2[0].c_str(), 20, binMin, binMax);
        tempHist[2] = new TH1F(nameSel2[1].c_str(), nameSel2[1].c_str(), 20, binMin, binMax);
        tempHist[3] = new TH1F(nameSel2[2].c_str(), nameSel2[2].c_str(), 20, binMin, binMax);
        tempHist[4] = new TH1F(nameSel2[3].c_str(), nameSel2[3].c_str(), 20, binMin, binMax);
        tempHist[5] = new TH1F(nameSel2[4].c_str(), nameSel2[4].c_str(), 20, binMin, binMax);
        tempHist[6] = new TH1F(nameSel2[5].c_str(), nameSel2[5].c_str(), 20, binMin, binMax);
        tempHist[7] = new TH1F(nameSel2[6].c_str(), nameSel2[6].c_str(), 20, binMin, binMax);
        tempHist[8] = new TH1F(nameSel2[7].c_str(), nameSel2[7].c_str(), 20, binMin, binMax);
        h_Acc.insert(std::make_pair(pdg, tempHist));
        return tempHist[0];
      }
    else
      return it_Acc->second[0];
  };

  auto f_createCDCFinal = [](std::unordered_map<std::string, std::vector<TH2F*> >& h_CDCFinal,
                             const std::string& nameH) -> std::vector<TH2F*>& {
    auto it_CDCFinal = h_CDCFinal.find(nameH);
    if(it_CDCFinal == h_CDCFinal.end())
      {
        std::string nameHist = nameH;
        std::string tempNameHist = nameH + "_X";
        std::vector<TH2F*> h_vecH(3, nullptr);
        h_vecH[0] = new TH2F(tempNameHist.c_str(), tempNameHist.c_str(), 20, 0, 20, 100, -50, 50);
        tempNameHist = nameH + "_Y";
        h_vecH[1] = new TH2F(tempNameHist.c_str(), tempNameHist.c_str(), 20, 0, 20, 100, -50, 50);
        tempNameHist = nameH + "_Z";
        h_vecH[2] = new TH2F(tempNameHist.c_str(), tempNameHist.c_str(), 20, 0, 20, 100, 60, 160);

        auto it_ret = h_CDCFinal.insert(std::make_pair(nameHist, h_vecH));

        return it_ret.first->second;
      }
    else
      {
        return it_CDCFinal->second;
      }
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
      std::vector<std::tuple<ParticleD,double> > daugthers;
      std::vector<std::tuple<ParticleD,double> > fake_daugthers;
      int pdg_frag = -1;
      for(Int_t id = 0; id < event->fMC_Particle->GetEntries(); ++id)
        {
          const TMcParticle& MCpar = *(dynamic_cast<TMcParticle*>(event->fMC_Particle->At(id)));
          TempData[id].MC_id = MCpar.Mc_id;
          TempData[id].Mother_id = MCpar.Mother_id;
          TempData[id].pdg = MCpar.Pdg;
          TempData[id].MomMass = MCpar.MomMass;
          auto PDG_particle = TDatabasePDG::Instance()->GetParticle(MCpar.Pdg);

	  h_Acceptance->Fill(PDG_particle->GetName(), "All", 1);

          TH1F* hMomAll = f_createAccHist(h_MomAcc, ParticleBinMM, MCpar.Pdg, "Mom");
          hMomAll->Fill(MCpar.MomMass.P());
          TH1F* hAngleAll = f_createAccHist(h_AngleAcc, ParticleBinAA, MCpar.Pdg, "Theta");
          hAngleAll->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
          TH1F* hYAll = f_createAccHist(h_RapidityAcc, ParticleBinYY, MCpar.Pdg, "Y");
          hYAll->Fill(MCpar.MomMass.Rapidity());
          if(MCpar.Mother_id >= 0)
            {
              h_Acceptance->Fill(PDG_particle->GetName(), "Decay", 1);
              TH1F* hMomAll2 = f_createAccHist(h_MomAccDecay, ParticleBinMM, MCpar.Pdg, "Decay_Mom");
              hMomAll2->Fill(MCpar.MomMass.P());
              TH1F* hAngleAll2 = f_createAccHist(h_AngleAccDecay, ParticleBinAA, MCpar.Pdg, "Decay_Theta");
              hAngleAll2->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());

              h_Acceptance->Fill(PDG_particle->GetName(), "DecayCoinFRS", 1);
	      TH1F* hMomAll3 = f_createAccHist(h_MomAccDecayCoinFRS, ParticleBinMM, MCpar.Pdg, "DecayCoinFRS_Mom");
              hMomAll3->Fill(MCpar.MomMass.P());
              TH1F* hAngleAll3 = f_createAccHist(h_AngleAccDecayCoinFRS, ParticleBinAA, MCpar.Pdg, "DecayCoinFRS_Theta");
              hAngleAll3->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());

	      if(MCpar.Pdg >= 2000)
		{
		  pdg_frag =  MCpar.Pdg;
		  ParticleD d_frag;
		  d_frag.MomMass = MCpar.MomMass;
		  d_frag.Vtx = MCpar.Vtx;
		  d_frag.pdg = pdg_frag;
		  daugthers.emplace_back(std::make_tuple(d_frag,MCpar.Charge));
		}
	    }
	  else
	    {
	      if(MCpar.Pdg >= 2000)
		{
		  ParticleD d_frag;
		  d_frag.MomMass = MCpar.MomMass;
		  d_frag.Vtx = MCpar.Vtx;
		  d_frag.pdg = MCpar.Pdg;
		  fake_daugthers.emplace_back(std::make_tuple(d_frag,MCpar.Charge));
		}
	    }
	}

      std::unordered_map<int, int> idDet;
      bool DaughterFragInFRS = false; 
      for(Int_t id = 0; id < event->FMF2->GetEntries(); ++id)
        {
          const TMcHit& MChit = *(dynamic_cast<TMcHit*>(event->FMF2->At(id)));
          int TrackID = MChit.MC_id;
          for(auto MCpar : TempData)
            if(TrackID == MCpar.MC_id)
              {

                auto it_det = idDet.find(TrackID);
                if(it_det == idDet.end())
                  {
                    MCpar.FinalHit.insert(std::make_pair("FMF2", MChit.MCHit));

                    auto PDG_particle = TDatabasePDG::Instance()->GetParticle(MChit.Pdg);
                    h_Acceptance->Fill(PDG_particle->GetName(), "FMF2", 1);
                    // h_Acceptance->Fill(PDG_particle->GetName(),"Det",1);
                    idDet.insert(std::make_pair(TrackID, MChit.LayerID));

                    auto it_momAcc = h_MomAcc.find(MChit.Pdg);
                    it_momAcc->second[3]->Fill(MCpar.MomMass.P());
                    // it_momAcc->second[7]->Fill(MCpar.MomMass.P());
                    auto it_angleAcc = h_AngleAcc.find(MChit.Pdg);
                    it_angleAcc->second[3]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
                    // it_angleAcc->second[7]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
                    auto it_YAcc = h_RapidityAcc.find(MChit.Pdg);
                    it_YAcc->second[3]->Fill(MCpar.MomMass.Rapidity());
                    // it_YAcc->second[7]->Fill(MCpar.MomMass.Rapidity());
                    if(MCpar.Mother_id >= 0)
                      {
                        // h_Acceptance->Fill(PDG_particle->GetName(),"Decay_Det",1);
                        auto it_momAcc2 = h_MomAccDecay.find(MChit.Pdg);
                        it_momAcc2->second[3]->Fill(MCpar.MomMass.P());
                        // it_momAcc2->second[7]->Fill(MCpar.MomMass.P());
                        auto it_angleAcc2 = h_AngleAccDecay.find(MChit.Pdg);
                        it_angleAcc2->second[3]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());

			if(MCpar.pdg >= 10000)
			  {
			    if(TMath::Abs(MChit.MCHit.X())<7 && TMath::Abs(MChit.MCHit.Y())<7)
			      {
				DaughterFragInFRS = true;
				h_Acceptance->Fill(PDG_particle->GetName(), "FMF2DecayCoinFRS", 1);

				auto it_momAcc3 = h_MomAccDecayCoinFRS.find(MChit.Pdg);
				it_momAcc3->second[3]->Fill(MCpar.MomMass.P());
				// it_momAcc2->second[7]->Fill(MCpar.MomMass.P());
				auto it_angleAcc3 = h_AngleAccDecayCoinFRS.find(MChit.Pdg);
				it_angleAcc3->second[3]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
			      }
			  }		
			// it_angleAcc2->second[7]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
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
        }

      for(Int_t id = 0; id < event->PSCE->GetEntries(); ++id)
        {
          const TMcHit& MChit = *(dynamic_cast<TMcHit*>(event->PSCE->At(id)));
          int TrackID = MChit.MC_id;
          for(auto& MCpar : TempData)
            if(TrackID == MCpar.MC_id)
              {
                MCpar.FinalHit.insert(std::make_pair("PSCE", MChit.MCHit));

                auto PDG_particle = TDatabasePDG::Instance()->GetParticle(MChit.Pdg);
                h_Acceptance->Fill(PDG_particle->GetName(), "PSCE", 1);
                h_Acceptance->Fill(PDG_particle->GetName(), "Det", 1);
                h_Acceptance->Fill(PDG_particle->GetName(), "Det2", 1);
                auto it_momAcc = h_MomAcc.find(MChit.Pdg);
                it_momAcc->second[4]->Fill(MCpar.MomMass.P());
                it_momAcc->second[7]->Fill(MCpar.MomMass.P());
                it_momAcc->second[8]->Fill(MCpar.MomMass.P());
                auto it_angleAcc = h_AngleAcc.find(MChit.Pdg);
                it_angleAcc->second[4]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
                it_angleAcc->second[7]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
                it_angleAcc->second[8]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
                auto it_YAcc = h_RapidityAcc.find(MChit.Pdg);
                it_YAcc->second[4]->Fill(MCpar.MomMass.Rapidity());
                it_YAcc->second[7]->Fill(MCpar.MomMass.Rapidity());
                it_YAcc->second[8]->Fill(MCpar.MomMass.Rapidity());
                if(MCpar.Mother_id >= 0)
                  {
                    h_Acceptance->Fill(PDG_particle->GetName(), "Decay_Det", 1);
                    h_Acceptance->Fill(PDG_particle->GetName(), "Decay_Det2", 1);
                    auto it_momAcc2 = h_MomAccDecay.find(MChit.Pdg);
                    it_momAcc2->second[4]->Fill(MCpar.MomMass.P());
                    it_momAcc2->second[7]->Fill(MCpar.MomMass.P());
                    it_momAcc2->second[8]->Fill(MCpar.MomMass.P());
                    auto it_angleAcc2 = h_AngleAccDecay.find(MChit.Pdg);
                    it_angleAcc2->second[4]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
                    it_angleAcc2->second[7]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
                    it_angleAcc2->second[8]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());

		    if(DaughterFragInFRS)
		      {
			h_Acceptance->Fill(PDG_particle->GetName(), "Decay_DetCoinFRS", 1);
			auto it_momAcc3 = h_MomAccDecayCoinFRS.find(MChit.Pdg);
			it_momAcc3->second[4]->Fill(MCpar.MomMass.P());
			it_momAcc3->second[7]->Fill(MCpar.MomMass.P());
			it_momAcc3->second[8]->Fill(MCpar.MomMass.P());
			auto it_angleAcc3 = h_AngleAccDecayCoinFRS.find(MChit.Pdg);
			it_angleAcc3->second[4]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
			it_angleAcc3->second[7]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
			it_angleAcc3->second[8]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
		      }
		  }
              }
        }

      for(Int_t id = 0; id < event->PSBE->GetEntries(); ++id)
        {
          const TMcHit& MChit = *(dynamic_cast<TMcHit*>(event->PSBE->At(id)));
          int TrackID = MChit.MC_id;
          for(auto MCpar : TempData)
            if(TrackID == MCpar.MC_id)
              {
                auto PDG_particle = TDatabasePDG::Instance()->GetParticle(MChit.Pdg);
                h_Acceptance->Fill(PDG_particle->GetName(), "PSBE", 1);
                // h_Acceptance->Fill(PDG_particle->GetName(),"Det",1);
                auto it_momAcc = h_MomAcc.find(MChit.Pdg);
                it_momAcc->second[5]->Fill(MCpar.MomMass.P());
                // it_momAcc->second[7]->Fill(MCpar.MomMass.P());
                auto it_angleAcc = h_AngleAcc.find(MChit.Pdg);
                it_angleAcc->second[5]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
                // it_angleAcc->second[7]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
                auto it_YAcc = h_RapidityAcc.find(MChit.Pdg);
                it_YAcc->second[5]->Fill(MCpar.MomMass.Rapidity());
                // it_YAcc->second[7]->Fill(MCpar.MomMass.Rapidity());
                if(MCpar.Mother_id >= 0)
                  {
                    // h_Acceptance->Fill(PDG_particle->GetName(),"Decay_Det",1);
                    auto it_momAcc2 = h_MomAccDecay.find(MChit.Pdg);
                    it_momAcc2->second[5]->Fill(MCpar.MomMass.P());
                    // it_momAcc2->second[7]->Fill(MCpar.MomMass.P());
                    auto it_angleAcc2 = h_AngleAccDecay.find(MChit.Pdg);
                    it_angleAcc2->second[5]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
                    // it_angleAcc2->second[7]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
                  }
              }
        }

      for(Int_t id = 0; id < event->PSFE->GetEntries(); ++id)
        {
          const TMcHit& MChit = *(dynamic_cast<TMcHit*>(event->PSFE->At(id)));
          int TrackID = MChit.MC_id;
          for(auto MCpar : TempData)
            if(TrackID == MCpar.MC_id)
              {
                auto PDG_particle = TDatabasePDG::Instance()->GetParticle(MChit.Pdg);
                h_Acceptance->Fill(PDG_particle->GetName(), "PSFE", 1);
                h_Acceptance->Fill(PDG_particle->GetName(),"Det2",1);
                auto it_momAcc = h_MomAcc.find(MChit.Pdg);
                it_momAcc->second[6]->Fill(MCpar.MomMass.P());
                it_momAcc->second[8]->Fill(MCpar.MomMass.P());
                auto it_angleAcc = h_AngleAcc.find(MChit.Pdg);
                it_angleAcc->second[6]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
                it_angleAcc->second[8]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
                auto it_YAcc = h_RapidityAcc.find(MChit.Pdg);
                it_YAcc->second[6]->Fill(MCpar.MomMass.Rapidity());
                it_YAcc->second[8]->Fill(MCpar.MomMass.Rapidity());
                if(MCpar.Mother_id >= 0)
                  {
                    h_Acceptance->Fill(PDG_particle->GetName(),"Decay_Det2",1);
                    auto it_momAcc2 = h_MomAccDecay.find(MChit.Pdg);
                    it_momAcc2->second[6]->Fill(MCpar.MomMass.P());
                    it_momAcc2->second[8]->Fill(MCpar.MomMass.P());
                    auto it_angleAcc2 = h_AngleAccDecay.find(MChit.Pdg);
                    it_angleAcc2->second[6]->Fill(MCpar.MomMass.Theta() * TMath::RadToDeg());
                    it_angleAcc2->second[8]->Fill(MCpar.MomMass.Theta()*TMath::RadToDeg());
                  }
              }
        }

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

      std::unordered_map<int, int> TrackMulti;
      for(Int_t id = 0; id < event->CDC->GetEntries(); ++id)
        {
          const TMcHit& MChit = *(dynamic_cast<TMcHit*>(event->CDC->At(id)));
          int TrackID = MChit.MC_id;
          auto it_trackM = TrackMulti.find(TrackID);
          if(it_trackM != TrackMulti.end())
            it_trackM->second += 1;
          else
            TrackMulti.insert(std::make_pair(TrackID, 1));
        }

      for(Int_t id = 0; id < event->CDC->GetEntries(); ++id)
        {
          const TMcHit& MChit = *(dynamic_cast<TMcHit*>(event->CDC->At(id)));
          int TrackID = MChit.MC_id;
          for(auto MCpar : TempData)
            if(TrackID == MCpar.MC_id)
              {
                auto PDG_particle = TDatabasePDG::Instance()->GetParticle(MChit.Pdg);
                auto it_trackMulti = TrackMulti.find(TrackID);
                if(it_trackMulti == TrackMulti.end())
                  std::cout << "!> trackMulti : unfound trackID :" << TrackID << "\n";

                for(auto finalHit_Det : MCpar.FinalHit)
                  {
                    std::string nameHist = finalHit_Det.first;
                    nameHist += std::to_string(MChit.Pdg);

                    auto it_FinalHit_CDC = f_createCDCFinal(h_FinalHit_CDC, nameHist);
                    it_FinalHit_CDC[0]->Fill(MChit.LayerID-6, finalHit_Det.second.X());
                    it_FinalHit_CDC[1]->Fill(MChit.LayerID-6, finalHit_Det.second.Y());
                    it_FinalHit_CDC[2]->Fill(MChit.LayerID-6, finalHit_Det.second.Z());

                    std::string nameHistM(nameHist);
                    nameHistM += "MultiplicityCDC";
                    auto it_FinalHit_CDCMulti = f_createCDCFinal(h_FinalHit_CDC, nameHistM);
                    it_FinalHit_CDCMulti[0]->Fill(it_trackMulti->second, finalHit_Det.second.X());
                    it_FinalHit_CDCMulti[1]->Fill(it_trackMulti->second, finalHit_Det.second.Y());
                    it_FinalHit_CDCMulti[2]->Fill(it_trackMulti->second, finalHit_Det.second.Z());

                    if(MCpar.Mother_id >= 0)
                      {
                        nameHist += "_Decay";
                        auto it_FinalHit_CDC_Decay = f_createCDCFinal(h_FinalHit_CDC, nameHist);
                        it_FinalHit_CDC_Decay[0]->Fill(MChit.LayerID-6, finalHit_Det.second.X());
			it_FinalHit_CDC_Decay[1]->Fill(MChit.LayerID-6, finalHit_Det.second.Y());
			it_FinalHit_CDC_Decay[2]->Fill(MChit.LayerID-6, finalHit_Det.second.Z());

                        nameHistM += "_Decay";
                        auto it_FinalHit_CDCMulti_Decay = f_createCDCFinal(h_FinalHit_CDC, nameHistM);
                        it_FinalHit_CDCMulti_Decay[0]->Fill(it_trackMulti->second, finalHit_Det.second.X());
                        it_FinalHit_CDCMulti_Decay[1]->Fill(it_trackMulti->second, finalHit_Det.second.Y());
                        it_FinalHit_CDCMulti_Decay[2]->Fill(it_trackMulti->second, finalHit_Det.second.Z());
			if(DaughterFragInFRS)
			  {
			    nameHist += "_DecayCoinFRS";
			    auto it_FinalHit_CDC_DecayFRS = f_createCDCFinal(h_FinalHit_CDC, nameHist);
			    it_FinalHit_CDC_DecayFRS[0]->Fill(MChit.LayerID-6, finalHit_Det.second.X());
			    it_FinalHit_CDC_DecayFRS[1]->Fill(MChit.LayerID-6, finalHit_Det.second.Y());
			    it_FinalHit_CDC_DecayFRS[2]->Fill(MChit.LayerID-6, finalHit_Det.second.Z());
			    
			    nameHistM += "_DecayCoinFRS";
			    auto it_FinalHit_CDCMulti_DecayFRS = f_createCDCFinal(h_FinalHit_CDC, nameHistM);
			    it_FinalHit_CDCMulti_DecayFRS[0]->Fill(it_trackMulti->second, finalHit_Det.second.X());
			    it_FinalHit_CDCMulti_DecayFRS[1]->Fill(it_trackMulti->second, finalHit_Det.second.Y());
			    it_FinalHit_CDCMulti_DecayFRS[2]->Fill(it_trackMulti->second, finalHit_Det.second.Z());

			  }
			
                      }
                  }
              }
        }
      
      for(Int_t id = 0; id < event->fTrack->GetEntries(); ++id)
        {
          const THyphiTrack& RecoTrack = *(dynamic_cast<THyphiTrack*>(event->fTrack->At(id)));
          int TrackID = RecoTrack.MC_status;
	  if(RecoTrack.Sim2Vtx.Z()>0 && TrackID>2000)
	    TrackID -= 10000;

          for(auto MCpar : TempData)
            if(TrackID == MCpar.MC_id)
              {
		auto PDG_particle = TDatabasePDG::Instance()->GetParticle(MCpar.pdg);
		if(MCpar.pdg<1000)
		  {
		    ParticleD d_pi;
		    //TLorentzVector p4;
		    d_pi.MomMass.SetXYZM(RecoTrack.MomMass.Px(),RecoTrack.MomMass.Py(),RecoTrack.MomMass.Pz(),PDG_particle->Mass());
		    d_pi.Vtx.SetXYZT(RecoTrack.RefPoint.X(),RecoTrack.RefPoint.Y(),RecoTrack.RefPoint.Z(),0.);
		    d_pi.pdg = MCpar.pdg;
		    daugthers.emplace_back(std::make_tuple(d_pi,PDG_particle->Charge()/3.));
		  }
                // if(RecoTrack.Pval2 < 0.000000001)
                //   continue;

		// auto f_extrap = [](const TVector3& dir,const TVector3& point, const double z, TVector3& point_in_z)
		// 		{
		// 		  TVector3 dir_temp(dir);
		// 		  if(TMath::Abs(dir.Z()-1.)>1e-9)
		// 		    dir_temp*=1./dir_temp.Z();
				  
		// 		  //cout<<" debug Line 3D :"<<"("<<dir_temp.X()<<" "<<dir_temp.Y()<<" "<<dir_temp.Z()<<")"<<endl;
				  
		// 		  point_in_z.SetX(dir_temp.X()*(z-point.Z())+point.X());
		// 		  point_in_z.SetY(dir_temp.Y()*(z-point.Z())+point.Y());
		// 		  point_in_z.SetZ(z);
		// 		};


		auto it_2DHist = h_MomAccReco.find(MCpar.pdg);
                if(it_2DHist == h_MomAccReco.end())
                  {
                    std::string namePar(PDG_particle->GetName());
                    std::string nameParHist(namePar);
                    std::replace(nameParHist.begin(), nameParHist.end(), '+', 'P');
                    std::replace(nameParHist.begin(), nameParHist.end(), '-', 'N');
                    nameParHist += "_RecoMom";
                    double binMin = 0.;
                    double binMax = 10.;
                    auto it_binMinMax = ParticleBinMM.find(namePar);
                    if(it_binMinMax != ParticleBinMM.end())
                      {
                        binMin = it_binMinMax->second[0];
                        binMax = it_binMinMax->second[1];
                      }
		    TH2F* h_Reco = new TH2F(nameParHist.c_str(), nameParHist.c_str(), 25, binMin, binMax, 1000, 0, 0.1);
                    h_MomAccReco.insert(std::make_pair(MCpar.pdg, h_Reco));
                    h_Reco->Fill(MCpar.MomMass.P(), std::fabs(RecoTrack.MomMass.P() - MCpar.MomMass.P()) / MCpar.MomMass.P());
                  }
                else
                  it_2DHist->second->Fill(MCpar.MomMass.P(), std::fabs(RecoTrack.MomMass.P() - MCpar.MomMass.P()) / MCpar.MomMass.P());
              }
        }

      auto f_Svtx = [](const ParticleD& d1, const ParticleD& d2)
		    {
		      double v1[3],v2[3],p1[3],p2[3];
		      double vec1[3],vec2[3];
		      
		      v1[0] = d1.MomMass.Px();
		      v1[1] = d1.MomMass.Py();
		      v1[2] = d1.MomMass.Pz();
		      p1[0] = d1.Vtx.X();
		      p1[1] = d1.Vtx.Y();
		      p1[2] = d1.Vtx.Z();
					      
		      v2[0] = d2.MomMass.Px();
		      v2[1] = d2.MomMass.Py();
		      v2[2] = d2.MomMass.Pz();
		      p2[0] = d2.Vtx.X();
		      p2[1] = d2.Vtx.Y();
		      p2[2] = d2.Vtx.Z();
					      
		      vec1[0]=v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
		      vec1[1]=-(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
		      vec1[2]=v1[0]*(p1[0]-p2[0])+v1[1]*(p1[1]-p2[1])+v1[2]*(p1[2]-p2[2]);
					      
		      vec2[0]=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
		      vec2[1]=-(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
		      vec2[2]=v2[0]*(p1[0]-p2[0])+v2[1]*(p1[1]-p2[1])+v2[2]*(p1[2]-p2[2]);
					      
		      auto f_solve = [](double v1[3], double v2[3]) {

				       double det = v1[0]*v2[1]-v1[1]*v2[0];
				       if (std::abs(det)<1e-5)
					 {
					   return std::make_tuple(-1,0.,0.);
					 }
				       else
					 {
					   double x=(v1[1]*v2[2]-v1[2]*v2[1])/det;
					   double y=(v1[2]*v2[0]-v1[0]*v2[2])/det;
					   return std::make_tuple(1,x,y);
					 }
				     };

		      int res=0;
		      double x=0,y=0, dist=0;
		      TVector3 point, poca1, poca2;
		      std::tie(res,x,y) = f_solve(vec1,vec2);
		      if(res==1)
			{
			  for (int i=0; i< 3; i++)
			    {
			      poca1[i]=x*v1[i]+p1[i];
			      poca2[i]=y*v2[i]+p2[i];
			      point[i]=(poca1[i]+poca2[i])/2.;
			      dist += (poca1[i]-poca2[i])*(poca1[i]-poca2[i]);
			    }
			  dist = std::sqrt(dist);
			}
		      return std::make_tuple(res,dist,point,poca1,poca2); 
		    };
      

      
      if(daugthers.size()>1 && h_InvMass != nullptr)
	{
	  TLorentzVector v4_mother;
	  ParticleD v4_piN;
	  bool have_pi = false;
	  for(auto v4_d : daugthers)
	    {
	      v4_mother += std::get<0>(v4_d).MomMass;
	      if(std::get<1>(v4_d)<0.)
		{
		  have_pi = true;
		  v4_piN = std::get<0>(v4_d);
		}
	    }
	  //v4_mother+=d_frag;
	  
	  if(daugthers.size()==2)
	    {
	      TVector3 Svtx,Poca1,Poca2;
	      int res;
	      double dist;
	      std::tie(res,dist,Svtx,Poca1,Poca2) = f_Svtx(std::get<0>(daugthers[0]),std::get<0>(daugthers[1]));
	      if(res==1)
		{
		  h_Vtx->Fill(dist);
		  h_RZVtx->Fill(Svtx.Perp(),Svtx.Z());
		  
		  h_InvMassVtx->Fill(v4_mother.M(),dist);
		  for(auto v4_d : daugthers)
		    if(std::get<0>(v4_d).MomMass.M() > .5)
		      {
			h_Vtx_X->Fill(Svtx.X()-std::get<0>(v4_d).Vtx.X());
			h_Vtx_Y->Fill(Svtx.Y()-std::get<0>(v4_d).Vtx.Y());
			h_Vtx_Z->Fill(Svtx.Z()-std::get<0>(v4_d).Vtx.Z());
		      }
		}
	    }
	  
	  h_InvMass->Fill(v4_mother.M());

	  for(auto v4_d : daugthers)
	    if(std::get<0>(v4_d).MomMass.M() > .5)
	      h_InvMassBrho->Fill(v4_mother.M(),3.10715497 * std::get<0>(v4_d).MomMass.P() / std::get<1>(v4_d));
	  //h_InvMassBrho->Fill(v4_mother.M(),3.10715497 * d_frag.P() / 2.);

	  if(fake_daugthers.size()>0 && have_pi==true)
	    {
	      for(auto v4_fake : fake_daugthers)
		{
		  TLorentzVector v4_m = v4_piN.MomMass;
		  v4_m += std::get<0>(v4_fake).MomMass;
		  auto it_hist = h_InvMassMix.find(std::get<0>(v4_fake).pdg);
		  if(it_hist != h_InvMassMix.end())
		    it_hist->second->Fill(v4_m.M());

		  auto it_hist_brho = h_InvMassBrhoMix.find(std::get<0>(v4_fake).pdg);
		  if(it_hist_brho != h_InvMassBrhoMix.end())
		    it_hist_brho->second->Fill(v4_m.M(),3.10715497 * std::get<0>(v4_fake).MomMass.P() / std::get<1>(v4_fake));
		  //else
		  //  std::cout<<"!> fake inv mass pdg not found !"<<std::get<0>(v4_fake).pdg<<std::endl;

		  TVector3 Svtx,Poca1,Poca2;
		  int res;
		  double dist;
		  std::tie(res,dist,Svtx,Poca1,Poca2) = f_Svtx(v4_piN,std::get<0>(v4_fake));
		  if(res==1)
		    {
		      h_VtxMix->Fill(dist);
		      h_InvMassVtxMix->Fill(v4_m.M(),dist);
		    }
		}
	    }	     
	}
      
      // cout<<"\n";
    }

  if(outfile.empty())
    {
      ploting2D(*h_Acceptance, "c1");

      //plotingAcceptance(h_MomAcc, "c1_MomAcc", true);

      //plotingAcceptance(h_AngleAcc, "c1_AngleAcc", true);

      //plotingAcceptance(h_RapidityAcc, "c1_RapidityAcc", true);

      //plotingAcceptance(h_MomAccDecay, "c1_MomAccDecay", true);

      //plotingAcceptance(h_AngleAccDecay, "c1_AngleAccDecay", true);

      plotingArray(h_MomAccReco, "c1_Reco", "candle");

      //plotingCDCFinal(h_FinalHit_CDC, "c2_");

      TCanvas* c_Inv = new TCanvas("c_Inv","c_Inv",1000,500); 
      c_Inv->Divide(2,1);
      c_Inv->cd(1);
      h_InvMass->Draw();
      c_Inv->cd(2);
      h_InvMassBrho->Draw();
    }
  else
    {
      TFile* fout = new TFile(outfile.c_str(), "RECREATE");
      fout->cd();
      h_Acceptance->Write();
      if(h_InvMass != nullptr)
	{
	  h_InvMass->Write();
	  h_Vtx->Write();
	  h_RZVtx->Write();
	  h_InvMassVtx->Write();
	  h_InvMassBrho->Write();
	  h_Vtx_X->Write();
	  h_Vtx_Y->Write();
	  h_Vtx_Z->Write();
	  for(auto it_hist : h_InvMassMix)
	    it_hist.second->Write();
	  for(auto it_hist : h_InvMassBrhoMix)
	    it_hist.second->Write();

	  h_VtxMix->Write();
	  h_InvMassVtxMix->Write();

	}
      auto f_DoEff = [](std::vector<TH1F*>& hist1d) -> std::vector<TH1F*> {

        std::vector<TH1F*> h_out(hist1d.size() - 1, nullptr);
	hist1d[0]->Sumw2();
        for(size_t i = hist1d.size() - 1; i >= 1; --i)
          {
            std::string nameAcc(hist1d[i]->GetName());
            nameAcc += "GeoAcceptance";
            // TH1F* htemp1 = new
            // TH1F(nameAcc.c_str(),nameAcc.c_str(),hist1d.second[0]->GetNbinsX(),hist1d.second[0]->GetXaxis()->GetXmin(),hist1d.second[0]->GetXaxis()->GetXmax());
            for(int n = 1; n <= hist1d[0]->GetNbinsX(); ++n)
              {
                Int_t Nall = hist1d[0]->GetBinContent(n);
                Int_t Nacc = hist1d[i]->GetBinContent(n);
                if(Nall <= 0)
                  {
                    Nacc = 0;
                    hist1d[i]->SetBinContent(n, Nacc);
                  }
                if(Nall < Nacc)
                  {
                    Nacc = Nall;
                    hist1d[i]->SetBinContent(n, Nacc);
                  }
              }
            if(hist1d[i]->GetEntries() > 1)
              {
                h_out[i - 1] = (TH1F*)hist1d[i]->Clone();
		h_out[i - 1]->Sumw2();
		h_out[i - 1]->GetXaxis()->SetRangeUser(hist1d[0]->GetXaxis()->GetXmin(), hist1d[0]->GetXaxis()->GetXmax());
                // TGraphAsymmErrors* g_acc = new TGraphAsymmErrors(hist1d.second[0],hist1d.second[1]);
                h_out[i - 1]->SetNameTitle(nameAcc.c_str(), nameAcc.c_str());
                h_out[i - 1]->Divide(hist1d[0]);
                if(i != hist1d.size() - 1)
                  h_out[i - 1]->SetLineColor(i);
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

	  int nb_hist = 0;

	  for(auto htemp : momEff)
            if(htemp != nullptr)
	      if(htemp->GetEntries()>0)
	      ++nb_hist;

	  for(auto htemp : angleEff)
            if(htemp != nullptr)
	      if( htemp->GetEntries()>0)
	      ++nb_hist;

          for(auto htemp : rapEff)
            if(htemp != nullptr)
	      if(htemp->GetEntries()>0)
	      ++nb_hist;
          if(it_momRes != h_MomAccReco.end())
	    if(it_momRes->second->GetEntries()>0)
	      ++nb_hist;

	  if( nb_hist == 0)
	    continue;
	  
	  TDirectory* dirTemp = fout->mkdir(nameDir.c_str(), nameDir.c_str());
          dirTemp->cd();
          for(auto htemp : momEff)
            if(htemp != nullptr)
	      if( htemp->GetEntries()>0)
              htemp->Write();

          for(auto htemp : angleEff)
            if(htemp != nullptr)
	      if(htemp->GetEntries()>0)
              htemp->Write();

          for(auto htemp : rapEff)
            if(htemp != nullptr)
	      if(htemp->GetEntries()>0)
              htemp->Write();

          if(it_momRes != h_MomAccReco.end())
	    if(it_momRes->second->GetEntries()>0)
	      it_momRes->second->Write();
        }

      for(auto it_momAcc : h_MomAccDecay)
        {
          int pdg = it_momAcc.first;
          auto PDG_particle = TDatabasePDG::Instance()->GetParticle(pdg);
          std::string nameDir(PDG_particle->GetName());
	  nameDir+="Decay";
          auto it_angleAcc = h_AngleAccDecay.find(pdg);

          std::vector<TH1F*> momEff = f_DoEff(it_momAcc.second);
          std::vector<TH1F*> angleEff = f_DoEff(it_angleAcc->second);

	  int nb_hist = 0;

	  for(auto htemp : momEff)
            if(htemp != nullptr)
	      if(htemp->GetEntries()>0)
	      ++nb_hist;

	  for(auto htemp : angleEff)
            if(htemp != nullptr)
	      if( htemp->GetEntries()>0)
	      ++nb_hist;

	  if( nb_hist == 0)
	    continue;

          TDirectory* dirTemp = fout->mkdir(nameDir.c_str(), nameDir.c_str());
          dirTemp->cd();

	  for(auto htemp : momEff)
            if(htemp != nullptr)
              if( htemp->GetEntries()>0)
		htemp->Write();

          for(auto htemp : angleEff)
            if(htemp != nullptr)
              if( htemp->GetEntries()>0)
		htemp->Write();

        }

      for(auto it_momAcc : h_MomAccDecayCoinFRS)
        {
          int pdg = it_momAcc.first;
          auto PDG_particle = TDatabasePDG::Instance()->GetParticle(pdg);
          std::string nameDir(PDG_particle->GetName());
	  nameDir+="DecayCoinFRS";
          auto it_angleAcc = h_AngleAccDecay.find(pdg);

          std::vector<TH1F*> momEff = f_DoEff(it_momAcc.second);
          std::vector<TH1F*> angleEff = f_DoEff(it_angleAcc->second);

	  int nb_hist = 0;

	  for(auto htemp : momEff)
            if(htemp != nullptr)
	      if(htemp->GetEntries()>0)
	      ++nb_hist;

	  for(auto htemp : angleEff)
            if(htemp != nullptr)
	      if( htemp->GetEntries()>0)
	      ++nb_hist;

	  if( nb_hist == 0)
	    continue;
	  
          TDirectory* dirTemp = fout->mkdir(nameDir.c_str(), nameDir.c_str());
          dirTemp->cd();

	  for(auto htemp : momEff)
            if(htemp != nullptr)
	      if( htemp->GetEntries()>0)
		htemp->Write();

          for(auto htemp : angleEff)
            if(htemp != nullptr)
	      if( htemp->GetEntries()>0)
		htemp->Write();

        }

      for(auto it_hist : h_FinalHit_CDC)
	{

	  auto it_posPDG = it_hist.first.find_first_of("-1234567890");
	  std::string nameDet = it_hist.first.substr(0, it_posPDG);
	      
	  std::string namePDG = it_hist.first.substr(it_posPDG, std::string::npos);
	  int pdg_code = std::stoi(namePDG);
	  auto PDG_particle = TDatabasePDG::Instance()->GetParticle(pdg_code);
	      
	  std::string namePar(PDG_particle->GetName());

	  // std::replace(namePar.begin(), namePar.end(), '+', 'P');
	  // std::replace(namePar.begin(), namePar.end(), '-', 'N');

	  auto getDir = [](TFile* out_file, const std::string& name_dir) {
			  if(!out_file->GetDirectory(name_dir.c_str()))
			    return out_file->mkdir(name_dir.c_str());
			  else
			    return out_file->GetDirectory(name_dir.c_str());
			};
	  std::string nameDir (namePar);
	      
	  auto it_decay = it_hist.first.find("Decay");
	  if(it_decay!=std::string::npos)
	    nameDir += "Decay";
	  auto it_frs = it_hist.first.find("CoinFRS");
	  if(it_frs != std::string::npos)
	    nameDir += "CoinFRS";

	  TDirectory* dirTemp = getDir(fout,nameDir);
	  dirTemp->cd();
	  for(size_t id = 0; id < it_hist.second.size(); ++id)
	    {
	      it_hist.second[id]->Write();
	    }
	      
	  
	

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

  auto f_delete_map = [](auto& maps) {
  			for(auto it : maps)
			  for(auto it1 : it.second)
			    {
			      it1->Reset();
			      it1->Delete();
			    }
		      };
  
  f_delete_map(h_MomAcc);
  f_delete_map(h_MomAccDecay);
  f_delete_map(h_MomAccDecayCoinFRS);
  f_delete_map(h_AngleAcc);
  f_delete_map(h_AngleAccDecay);
  f_delete_map(h_AngleAccDecayCoinFRS);
  f_delete_map(h_RapidityAcc);

  f_delete_map(h_FinalHit_CDC);
  for(auto it : h_MomAccReco)
    {
      it.second->Reset();
      it.second->Delete();
    }
  f_delete_map(h_MomAccReco_multi);
  for(auto it : h_InvMassBrho_multi)
    {
      it->Reset();
      it->Delete();
    }

  h_Acceptance->Reset();
  if(h_InvMass != nullptr)
    {
      h_Vtx->Reset();
      h_RZVtx->Reset();
      h_InvMass->Reset();
      h_InvMassVtx->Reset();
      h_InvMassBrho->Reset();
      h_Vtx_X->Reset();
      h_Vtx_Y->Reset();
      h_Vtx_Z->Reset();
      for(auto& it : h_InvMassMix)
	it.second->Reset();
      for(auto& it : h_InvMassBrhoMix)
	it.second->Reset();
    }
  h_Acceptance->Delete();
  if(h_InvMass != nullptr)
    {
      h_Vtx->Delete();
      h_RZVtx->Delete();
      h_InvMass->Delete();
      h_InvMassVtx->Delete();
      h_InvMassBrho->Delete();
      h_Vtx_X->Delete();
      h_Vtx_Y->Delete();
      h_Vtx_Z->Delete();
    }
}
