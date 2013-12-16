//***** Simple macro to compute systematics on di-jet tagged categories from PU-jetID efficiency.
// Inputs needed: 
// (1) file containing mini-tree with jet related variables for ggh, vbf, wzh, tth
// (2) data/MC scale factors derived from Z->mumu, e.g.: AnalysisScripts/aux/JetIdScaleFactor_ZmumuJets_12fb_hcp.root

#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TH1D.h"
#include "TCanvas.h"

int etaregion(float jeteta){
  int region=-1;
  if (fabs(jeteta)<2.50 )                     region = 0; //TK
  if (fabs(jeteta)>2.50 && fabs(jeteta)<2.75) region = 1; //HEin
  if (fabs(jeteta)>2.75 && fabs(jeteta)<3.0)  region = 2; //HEout
  if (fabs(jeteta)>3.00 )                     region = 3; //HF
  return region;
}


void makeJetIdEffUncertainties(string infilename, int mass=125, bool isMVA=true)
{
  
  TFile *inFile = TFile::Open(infilename.c_str());

  //--- process names
  vector<string> procnames;
  procnames.push_back("ggH");
  procnames.push_back("qqH");
  procnames.push_back("wzH");
  procnames.push_back("ttH");

  vector<int> procid;
  procid.push_back(0);
  procid.push_back(1);
  procid.push_back(2);
  procid.push_back(3);
  
  //--- number of categories
  int ncats = 14; // MVA - 8TeV
  if (!isMVA) ncats = 16; // CiC - TeV

  //--- book histograms
  TH1F *h[4][20];
  TH1F *hs[4][20];
  char histo[100];
  for (int iproc=0; iproc < 4; iproc++){
    for (int icat=0; icat < ncats; icat++){
      // -- histograms
      sprintf(histo,"h_%s_cat%d", (procnames.at(iproc)).c_str(), icat);
      h[iproc][icat]= new TH1F(histo, histo, 100, -1,1);
      // -- histograms after re-weighting for PU-jetID
      sprintf(histo,"hs_%s_cat%d",(procnames.at(iproc)).c_str(), icat);
      hs[iproc][icat]= new TH1F(histo, histo, 100, -1,1);
    }
  }

  //--- load data/MC scale factors
  TFile *f = TFile::Open("../../../AnalysisScripts/aux/JetIdScaleFactor_ZmumuJets_legacypaper.root"); 
  //TFile *f = TFile::Open("../../../AnalysisScripts/aux//JetIdScaleFactor_ZmumuJets_12fb_hcp.root");
  TH1F *hJetIdScaleFactor[4];
  hJetIdScaleFactor[0]=(TH1F*)f->Get("hJetIdScaleFactor_TK");
  hJetIdScaleFactor[1]=(TH1F*)f->Get("hJetIdScaleFactor_HEin");
  hJetIdScaleFactor[2]=(TH1F*)f->Get("hJetIdScaleFactor_HEout");
  hJetIdScaleFactor[3]=(TH1F*)f->Get("hJetIdScaleFactor_HF");

  float w1 = 1.; // weight for leading jet
  float w2 = 1.; // weight for sub-leading jet
  float w  = 1.;
  int proc = 0;
  int cat  = 0;
  int r1, r2;  

  //--- load trees
  vector<TTree*> trees;
  trees.push_back((TTree*)inFile->Get(Form("vbf_trees/ggh_m%d_8TeV",mass)));
  trees.push_back((TTree*)inFile->Get(Form("vbf_trees/vbf_m%d_8TeV",mass)));
  trees.push_back((TTree*)inFile->Get(Form("vbf_trees/wzh_m%d_8TeV",mass)));
  trees.push_back((TTree*)inFile->Get(Form("vbf_trees/tth_m%d_8TeV",mass)));

  //  float leadPtOverM, subleadPtOverM;
  float leadJPt, subleadJPt;
  float leadJEta, subleadJEta;
  //float vbfMVA, diphotonMVA;
  int category ;
  float weight;
  //float zepp;
  //float MJJ;
  //float deltaPhiJJGamGam;
  //float deltaEtaJJ;

  TTree *tree;      

  for (int itree =0; itree < trees->size(); itree++){

    tree = trees.at(itree);

    //tree->SetBranchAddress("sampleType",&sampleType);
    //tree->SetBranchAddress("leadPtOverM",&leadPtOverM);
    //tree->SetBranchAddress("subleadPtOverM",&subleadPtOverM);
    tree->SetBranchAddress("leadJPt",&leadJPt);
    tree->SetBranchAddress("subleadJPt",&subleadJPt);
    tree->SetBranchAddress("leadJEta",&leadJEta);
    tree->SetBranchAddress("subleadJEta",&subleadJEta);
    //tree->SetBranchAddress("MJJ",&MJJ);
    //tree->SetBranchAddress("deltaEtaJJ",&deltaEtaJJ);
    //tree->SetBranchAddress("deltaPhiJJGamGam",&deltaPhiJJGamGam);
    //tree->SetBranchAddress("zepp",&zepp);
    //tree->SetBranchAddress("vbfMVA",&vbfMVA);
    //tree->SetBranchAddress("diphotonMVA",&diphotonMVA);
    tree->SetBranchAddress("category",&category);
    tree->SetBranchAddress("weight",&weight);
    

    //--- loop over entries
    for (int ientry = 0 ; ientry < tree->GetEntries(); ientry++ ){
      
      tree->GetEntry(ientry);
      
      if (ientry%10000==0) cout << "Analyzing entry:" << ientry << endl;
      
      //-- 
      if (!(leadJPt>0 && subleadJPt>0)) continue; 
      
      w1=1;
      w2=1;
      w =1;
      r1=-1;
      r2=-1;
        
      //-- select event in category 4 or 5 (di-jet tagged categories)
      bool isVBF = false;
      if (isMVA) 
	isVBF = (category==5 || category==6 || category==7);
      else
	isVBF = (category==8 || category==9 );

      if ( isVBF )  { 
	//-- find the jet eta pseudorapidity region
	r1=etaregion(leadJEta);
	r2=etaregion(subleadJEta);

	if (r1!=-1) {
	  w1 = hJetIdScaleFactor[r1]->GetBinContent(hJetIdScaleFactor[r1]->FindBin(leadJPt));
	  if (leadJPt>100) w1 = hJetIdScaleFactor[r1]->GetBinContent(hJetIdScaleFactor[r1]->GetNbinsX());
	}
	if (r2!=-1) {
	  w2 = hJetIdScaleFactor[r2]->GetBinContent(hJetIdScaleFactor[r2]->FindBin(subleadJPt));
	  if (subleadJPt>100) w2 = hJetIdScaleFactor[r2]->GetBinContent(hJetIdScaleFactor[r2]->GetNbinsX());
	}
	
	// -- to be conservative if the data/MC scale factor is > 1, do not correct. 
	if (w1>1) w1=1.;
	if (w2>1) w2=1.;
	w  = w1*w2;
	
	//-- fill histograms with weights
	proc = procid.at(itree);
	h[proc][category]->Fill(0, weight);
	hs[proc][category]->Fill(0, (w*weight));
      }
    }
  }  

  cout << "proc     cat    N(nominal)     N(effW)     diff(%) "  << endl; 
  float r;
  int catmin = 5;
  int catmax = 7;
  if (!isMVA){
    catmin = 8;
    catmax=9;
  }


  for (int iproc=0; iproc < 4; iproc++){
    for (int icat=catmin; icat < (catmax+1); icat++){
      r=0;
      if (hs[iproc][icat]-> GetSumOfWeights()>0 && h[iproc][icat]-> GetSumOfWeights()>0)
	r = hs[iproc][icat]-> GetSumOfWeights()/h[iproc][icat]-> GetSumOfWeights()-1 ;
      cout << (procnames.at(iproc)).c_str() << "   " 
	   << icat << "   " 
	   << setprecision(3) << fixed << h[iproc][icat] -> GetSumOfWeights() << "   " 
	   << setprecision(3) << fixed << hs[iproc][icat]-> GetSumOfWeights() << "   "
	   << setprecision(2) << fixed << r*100 << "%" 
	   << endl;
    }
    cout << endl;
  }
  
  
}
