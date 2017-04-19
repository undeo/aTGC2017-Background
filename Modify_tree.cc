#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TFile.h>
#include <TString.h>
#include <TSystem.h>
#include <TMath.h>
#include <iostream>

void merge_trees(std::vector<TString> files, TChain &chain)
{
	for(unsigned int i=0; i<files.size(); i++)
		chain.Add(files[i]);
}

void fill_tree_with_cuts(TTree &oldtree, TTree &tree, TString ch)
{
	Double_t jet_pt,jet_tau2tau1,Mjpruned,W_pt,deltaR_LeptonWJet,deltaPhi_WJetMet,deltaPhi_WJetWlep,pfMET,METCUT;
	Int_t nbtag;
	METCUT	= ch=="ele" ? 80. : 40;

	oldtree.SetBranchAddress("jet_pt",&jet_pt);
	oldtree.SetBranchAddress("jet_tau2tau1",&jet_tau2tau1);
	oldtree.SetBranchAddress("Mjpruned",&Mjpruned);
	oldtree.SetBranchAddress("W_pt",&W_pt);
	oldtree.SetBranchAddress("deltaR_LeptonWJet",&deltaR_LeptonWJet);
	oldtree.SetBranchAddress("deltaPhi_WJetMet",&deltaPhi_WJetMet);
	oldtree.SetBranchAddress("deltaPhi_WJetWlep",&deltaPhi_WJetWlep);
	oldtree.SetBranchAddress("nbtag",&nbtag);
	oldtree.SetBranchAddress("pfMET",&pfMET);

	Long64_t nEntries	= oldtree.GetEntries();

	int used_events		= 0;

	for(unsigned int i=0; i<nEntries; i++)
	{
		oldtree.GetEntry(i);
		if(jet_pt>200. && jet_tau2tau1<0.6 && Mjpruned<150. && Mjpruned>40. && W_pt>200. && fabs(deltaR_LeptonWJet)>TMath::Pi()/2 && fabs(deltaPhi_WJetMet)>2. && fabs(deltaPhi_WJetWlep)>2. && nbtag==0 && pfMET>METCUT)
		{
			tree.Fill();
			if(used_events%1000==0 && i!=0)
				std::cout<<used_events<<"/"<<i<<"/"<<nEntries<<", "<<(float)i/(float)nEntries*100<<"%"<<", used: "<<(float)used_events/(float)i*100<<"%"<<std::endl;
			used_events++;
		}
	}
	std::cout<<"using "<<used_events<<" events of total "<<nEntries<<std::endl;
}

void save_tree(TString name, TTree &tree, TString ch)
{
	TString path	= "../InputTrees/"+ch(0,2)+"/";
	TFile * fileOut	= TFile::Open(path+"tree_"+name+"_"+ch(0,2)+".root","recreate");
	tree.Write();
	fileOut->Close();
	std::cout<<"saved "<<name<<" as "<<path+"tree_"+name+"_"+ch(0,2)+".root"<<std::endl;
}

void make_trees(TString ch)
{

	std::vector<TString> WJets_names;	
	std::vector<TString> TTbar_names;
	std::vector<TString> STop_names;
	std::vector<TString> WW_names;
	std::vector<TString> WZ_names;
	std::vector<TString> data_names;
	
	WJets_names.push_back("WJets_HT-100To200-tot_"+ch+".root");
	WJets_names.push_back("WJets_HT-200To400-tot_"+ch+".root");
	WJets_names.push_back("WJets_HT-400To600_"+ch+".root");
	WJets_names.push_back("WJets_HT-600To800_"+ch+".root");
	WJets_names.push_back("WJets_HT-800To1200-tot_"+ch+".root");
	WJets_names.push_back("WJets_HT-1200To2500-tot_"+ch+".root");
	WJets_names.push_back("WJets_HT-2500ToInf-tot_"+ch+".root");

	TTbar_names.push_back("ttbar-powheg-tot_"+ch+".root");

	STop_names.push_back("s-ch_"+ch+".root");
	STop_names.push_back("tW-ch-antitop_"+ch+".root");
	STop_names.push_back("tW-ch-top_"+ch+".root");
	STop_names.push_back("t-ch-tot_"+ch+".root");

	WW_names.push_back("WW-tot_"+ch+".root");
	
	WZ_names.push_back("WZ_"+ch+".root");

	data_names.push_back("data-RunD_"+ch+".root");


	TChain WJets_old("BasicTree");
	merge_trees(WJets_names,WJets_old);
	TTree * WJets_new	= WJets_old.CloneTree(0);
	std::cout<<"reding W+jets..."<<std::endl;
	fill_tree_with_cuts(WJets_old,*WJets_new,ch);

	std::cout<<"reading TTbar..."<<std::endl;
	TChain TTbar_old("BasicTree");
	merge_trees(TTbar_names,TTbar_old);
	TTree * TTbar_new	= TTbar_old.CloneTree(0);
	fill_tree_with_cuts(TTbar_old,*TTbar_new,ch);	

	std::cout<<"reading STop..."<<std::endl;
	TChain STop_old("BasicTree");
	merge_trees(STop_names,STop_old);
	TTree * STop_new	= STop_old.CloneTree(0);
	fill_tree_with_cuts(STop_old,*STop_new,ch);

	std::cout<<"reading WW..."<<std::endl;
	TChain WW_old("BasicTree");
	merge_trees(WW_names,WW_old);
	TTree * WW_new		= WW_old.CloneTree(0);
	fill_tree_with_cuts(WW_old,*WW_new,ch);
	
	std::cout<<"reading WZ..."<<std::endl;
	TChain WZ_old("BasicTree");
	merge_trees(WZ_names,WZ_old);
	TTree * WZ_new		= WZ_old.CloneTree(0);
	fill_tree_with_cuts(WZ_old,*WZ_new,ch);
	
	std::cout<<"reading data..."<<std::endl;
	TChain data_old("treeDumper/BasicTree");
	merge_trees(data_names,data_old);
	TTree * data_new	= data_old.CloneTree(0);
	fill_tree_with_cuts(data_old,*data_new,ch);


	save_tree("WJets",*WJets_new,ch);
	save_tree("TTbar",*TTbar_new,ch);
	save_tree("STop",*STop_new,ch);
	save_tree("WW",*WW_new,ch);
	save_tree("WZ",*WZ_new,ch);
	save_tree("data",*data_new,ch);
			
}

void Modify_tree()
{
	std::cout<<"reading electron trees"<<std::endl;
	make_trees("ele");
	std::cout<<"reading muon trees"<<std::endl;
	make_trees("mu");
}

