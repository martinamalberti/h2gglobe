///==== include ====

#include "TFile.h"
#include "TChain.h"
#include "TMinuit.h"
#include <sstream>
#include <iostream>
#include "TLorentzVector.h"
#include "TMVA/Factory.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/MethodCategory.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif

using namespace std;
 
// --------- MAIN -------------------

int main(int argc, char** argv)
{ 
  // TTree
  string inputfile = argv[1];
  TChain* tree = new TChain("vtxTree/ggh_m125_14TeV");
  //tree->Add("~/eos/cms//store/cmst3/user/malberti/HIGGS/Upgrade/vtxPU140_testTiming_50ps/histograms_vertexOpt.root");
  tree->Add(inputfile.c_str());

  // options
  bool  useTiming    = atoi(argv[2]);
  bool  useDeltaEta  = atoi(argv[3]);
  
  
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

    
  tree->SetBranchAddress("logsumpt2",     &logsumpt2);
  tree->SetBranchAddress("ptbal",         &ptbal);
  tree->SetBranchAddress("ptasym",        &ptasym);
  tree->SetBranchAddress("limpulltoconv", &limpulltoconv);
  tree->SetBranchAddress("nconv",         &nconv);
  tree->SetBranchAddress("tReco1",        &tReco1);
  tree->SetBranchAddress("tReco2",        &tReco2);
  tree->SetBranchAddress("dtof1",         &dtof1);
  tree->SetBranchAddress("dtof2",         &dtof2);
  tree->SetBranchAddress("pho1",          &pho1);
  tree->SetBranchAddress("pho2",          &pho2);
  tree->SetBranchAddress("isClosestToGen",&isClosestToGen);


  // Create a new root output file.
  string outputFileName = argv[4];

  TFile* outputFile = TFile::Open((outputFileName+".root").c_str(), "RECREATE" );
  TMVA::Factory* factory = new TMVA::Factory(outputFileName.c_str(), outputFile, 
					     "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );


  // -- variables
  factory->AddVariable( "logsumpt2");
  factory->AddVariable( "ptbal" );
  factory->AddVariable( "ptasym" );
  factory->AddVariable( "limpulltoconv" );
  if (useTiming){
    factory->AddVariable( "dtof:=(tReco1-dtof1-(tReco2-dtof2))" );
    if (useDeltaEta){
      factory->AddVariable( "deta:=abs(pho1.Eta()-pho2.Eta())");
    }
  }

  // -- spectators
  factory->AddVariable( "nconv" );
  //factory->AddSpectator("diphoM:=((pho1+pho2).M())");
  factory->AddSpectator("diphoPt:=(sqrt(pow(pho1.Px()+pho2.Px(),2)+pow(pho1.Py()+pho2.Py(),2)))");
    

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
  //				       "nTrain_Signal=1000:nTrain_Background=1000:nTest_Signal=1000:nTest_Background=1000:SplitMode=Random:NormMode=NumEvents:!V" );

   
   
  // Boosted Decision Trees: use BDTG ( Gradient Boost )
  factory->BookMethod( TMVA::Types::kBDT, "BDTG",
		       "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5:MaxDepth=3" );



  // --- Categorised classifier
  TMVA::MethodCategory* mcat = 0;

  // The variable sets
  TString theCat1Vars = "logsumpt2:ptbal:ptasym";
  TString theCat2Vars = "logsumpt2:ptbal:ptasym:limpulltoconv";
  if (useTiming){
    if (useDeltaEta){
      theCat1Vars = "logsumpt2:ptbal:ptasym:dtof:deta";
      theCat2Vars = "logsumpt2:ptbal:ptasym:limpulltoconv:dtof:deta";
    }
    else{
      theCat1Vars = "logsumpt2:ptbal:ptasym:dtof";
      theCat2Vars = "logsumpt2:ptbal:ptasym:limpulltoconv:dtof";
    }
  }

  // BDTG with categories
  TMVA::MethodBase* bdtgCat = factory->BookMethod( TMVA::Types::kCategory, "BDTGCat","" );
  mcat = dynamic_cast<TMVA::MethodCategory*>(bdtgCat);
  mcat->AddMethod( "nconv==0",theCat1Vars, TMVA::Types::kBDT,
		   "Category_BDTG_1","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5:MaxDepth=3" );
  mcat->AddMethod( "nconv>0",theCat1Vars, TMVA::Types::kBDT,
		   "Category_BDTG_2","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5:MaxDepth=3" );
		   


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
