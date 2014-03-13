///==== include ====

#include "TFile.h"
#include "TChain.h"
#include "TMinuit.h"
#include <sstream>
#include <iostream>
#include "TLorentzVector.h"
#include "TMVA/Factory.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif

using namespace std;
 
// --------- MAIN -------------------

int main(int argc, char** argv)
{ 

  // options
  bool  useDiphotonPt = true;
  bool  usePhotonsPt  = true;
  
  // Chain
  TChain* tree = new TChain("vtxTree/ggh_m125_14TeV");
  tree->Add("/afs/cern.ch/work/m/malberti/HGG/Upgrade/CMSSW_6_1_2/src/h2gglobe/AnalysisScripts/vertexOpt.root");
  
  
  // Declaration of leaf types
  TLorentzVector *pho1=0;
  TLorentzVector *pho2=0;

  float logsumpt2 ;
  float ptbal;
  float ptasym;
  float limpulltoconv;
  float nconv;
  float tReco1;
  float tReco2;
  float dtof1;
  float dtof2;
  bool isClosestToGen;

    
  tree->SetBranchAddress("logsumpt2",    &logsumpt2);
  tree->SetBranchAddress("ptbal",        &ptbal);
  tree->SetBranchAddress("ptasym",       &ptasym);
  tree->SetBranchAddress("limpulltoconv",&limpulltoconv);
  tree->SetBranchAddress("nconv",        &nconv);
  tree->SetBranchAddress("tReco1",       &tReco1);
  tree->SetBranchAddress("tReco2",       &tReco2);
  tree->SetBranchAddress("dtof1",        &dtof1);
  tree->SetBranchAddress("dtof2",        &dtof2);
  tree->SetBranchAddress("pho1",         &pho1);
  tree->SetBranchAddress("pho2",         &pho2);


  // Create a new root output file.
  string outputFileName = argv[1];

  TFile* outputFile = TFile::Open((outputFileName+".root").c_str(), "RECREATE" );
  TMVA::Factory* factory = new TMVA::Factory(outputFileName.c_str(), outputFile, 
					     "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );


  // -- variables

  factory->AddVariable( "logsumpt2");
  factory->AddVariable( "ptbal" );
  factory->AddVariable( "ptasym" );
  //factory->AddVariable( "dtof:=(tReco1-dtof1-(tReco2-dtof2))" );
  //factory->AddVariable( "deta:=abs(pho1.Eta()-pho2.Eta())");
  //  factory->AddVariable( "limpulltoconv" );
  //  factory->AddVariable( "nconv" );
  
  
  // -- spectators
  //factory->AddSpectator("diphoM");
    

  //event weights per tree (see below for setting event-wise weights)
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;
  

  // ====== To give different trees for training and testing, do as follows:
  //TFile *tmp = new TFile("tmvatrainingtemp.root","RECREATE");
  //tmp->ls();
   
  //TTree *trainingSignal     = tree->CopyTree("isClosestToGen==1");
  //TTree *trainingBackground = tree->CopyTree("isClosestToGen==0");
  //TTree *testSignal     = tree->CopyTree("isClosestToGen==1");
  //TTree *testBackground = tree->CopyTree("isClosestToGen==0");
  //factory->AddSignalTree( trainingSignal, signalWeight, "Training" );
  //factory->AddSignalTree( testSignal,     signalWeight,  "Test" );
  //factory->AddBackgroundTree( trainingBackground, backgroundWeight, "Training" );
  //factory->AddBackgroundTree( testBackground,     backgroundWeight,  "Test" );
    
  // ====== register trees ====================================================
  //
  // the following method is the prefered one:
  // you can add an arbitrary number of signal or background trees
  factory->AddSignalTree    ( tree,     signalWeight     );
  factory->AddBackgroundTree( tree,     backgroundWeight );
   

  // Apply additional cuts on the signal and background samples (can be different)
  TCut mycuts = "pho1.Pt()>25 && pho2.Pt()>25 && isClosestToGen==1"; 
  TCut mycutb = "pho1.Pt()>25 && pho2.Pt()>25 && isClosestToGen==0"; 

   
  // tell the factory to use all remaining events in the trees after training for testing:
  factory->PrepareTrainingAndTestTree( mycuts, mycutb,
				       "SplitMode=Random:NormMode=NumEvents:!V" );
   
   
  // Boosted Decision Trees: use BDTG ( Gradient Boost )
  factory->BookMethod( TMVA::Types::kBDT, "BDTG",
		       "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5:MaxDepth=3" );
  //"!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=15:MaxDepth=5" );
   
   
  // book Cuts
  //factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
  //"H:!V:FitMethod=GA:CutRangeMin[0]=20:CutRangeMax[0]=500:CutRangeMin[1]=20:CutRangeMax[1]=500:VarProp=FSmart:VarProp[4]=FMin:EffSel:Steps=30:Cycles=3:PopSize=500:SC_steps=10:SC_rate=5:SC_factor=0.95" );

  // ---- Now you can tell the factory to train, test, and evaluate the MVAs

  // Train MVAs using the set of training events
  factory->TrainAllMethods();

  // ---- Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

  // ----- Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();

  // --------------------------------------------------------------

  // Save the output
  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

  delete factory;
}
