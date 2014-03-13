#ifndef __SimpleVertexAnalysis__
#define __SimpleVertexAnalysis__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis/interface/PhotonAnalysis.h"
#include "PhotonAnalysis/interface/StatAnalysis.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"

#include "EnergySmearer.h"
#include "EfficiencySmearer.h"
#include "DiPhoEfficiencySmearer.h"
#include "KFactorSmearer.h"
#include <iostream>
#include <fstream>
#include "math.h"
#include "TRandom.h"


// ------------------------------------------------------------------------------------
class SimpleVertexAnalysis : public StatAnalysis 
{
 public:
    
    SimpleVertexAnalysis();
    virtual ~SimpleVertexAnalysis();
    
    virtual const std::string & name() const { return name_; };
    
    // LoopAll analysis interface implementation
    void Init(LoopAll&);
    void Term(LoopAll&);
    
    void ReducedOutputTree(LoopAll &l, TTree * outputTree);
    void GetBranches(TTree *t, std::set<TBranch *>& s ) ;

    void FillReductionVariables(LoopAll& l, int jentry){};
    bool SelectEventsReduction(LoopAll&, int){};
    bool SkimEvents(LoopAll& l, int jentry){};


    virtual bool AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, 
			      float & mass, float & evweight, int & category, int & diphoton_id,
			      bool & isCorrectVertex, float &kinematic_bdtout,
			      bool isSyst=false, 
			      float syst_shift=0., bool skipSelection=false,
			      BaseGenLevelSmearer *genSys=0, BaseSmearer *phoSys=0, BaseDiPhotonSmearer * diPhoSys=0); 
    
    TString minBiasRefName;
    int storeNVert;
    bool runOnReducedNtuples;
    double sigmaT;

private:

    TRandom *rnd;

    std::vector<std::string> vtxVarNames_;
    std::vector<float> vtxVars_;
    
    TFile * uFile_;
    TTree * uTree_;
    TH1 * hMinBias_, *hHiggs_, *hMinBiasRef_[4];
    bool isClosestToGen_, passCiC_;
    int   nPU_, nVert_, itype_;
    float evWeight_, ksprob_[4];
    vector<float> alpha02_, alpha03_, alpha04_,alpha05_;
    TLorentzVector *pho1_;
    TLorentzVector *pho2_;
    TLorentzVector *dipho_;
    float dtof1_, dtof2_;
    float tof1_, tof2_;
    float tReco1_;
    float tReco2_;
    float vtxZ_;
    
    TTree *evTree_;
    Float_t dZTrue_, zRMS_, zTrue_;
    Int_t category_, nConv_;
    Float_t mTrue_, mTrueVtx_;
    vector<float> MVA_;
    vector<float> dZ_;
    vector<float> diphoM_;
    vector<float> diphoPt_;
    
};

#endif


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
