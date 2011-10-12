#ifndef __PHOTONANALYSIS__
#define __PHOTONANALYSIS__

#include "BaseAnalysis.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"
#include "VertexAnalysis/interface/PhotonInfo.h"

#include "TriggerSelection.h"
#include "EnergySmearer.h"
#include "MassResolution.h"

// ------------------------------------------------------------------------------------
class PhotonAnalysis : public BaseAnalysis 
{
public:
	
	PhotonAnalysis();
	virtual ~PhotonAnalysis();
	
	virtual const std::string & name() const { return name_; };
	
	// LoopAll analysis interface implementation
	virtual void Init(LoopAll&);
	virtual void Term(LoopAll&);
	
	virtual void ReducedOutputTree(LoopAll &l, TTree *);
	virtual void GetBranches(TTree *, std::set<TBranch *>& );
	
	virtual void FillReductionVariables(LoopAll& l, int jentry);   
	virtual bool SelectEventsReduction(LoopAll&, int);

	virtual bool SkimEvents(LoopAll&, int);
	virtual bool SelectEvents(LoopAll&, int);
	virtual void ResetAnalysis();
	virtual void Analysis(LoopAll&, Int_t);
	
	// Public parameters to be read from config file
	VertexAlgoParameters vtxAlgoParams;	 
	std::vector<std::string> vtxVarNames;
	bool useDefaultVertex;
	float forcedRho;
	
	bool doTriggerSelection; 
	std::vector<TriggerSelection> triggerSelections;
	
	// Preselection indexes
	float presel_scet1, presel_scet2, presel_maxeta;
	float presel_ecaliso_eb, presel_ecaliso_ee, presel_sieie_eb, presel_sieie_ee, presel_hoe;

	bool doEcorrectionSmear, doEcorrectionSyst;
	
	EnergySmearer::energySmearingParameters eSmearDataPars;
	std::string scale_offset_file;
	float scale_offset_EBHighR9         ;
	float scale_offset_EBLowR9          ;
	float scale_offset_EEHighR9         ;
	float scale_offset_EELowR9          ;
	float scale_offset_error_EBHighR9   ;
	float scale_offset_error_EBLowR9    ;
	float scale_offset_error_EEHighR9   ;
	float scale_offset_error_EELowR9    ;

	EnergySmearer::energySmearingParameters eSmearPars;
	float smearing_sigma_EBHighR9       ;
	float smearing_sigma_EBLowR9        ;
	float smearing_sigma_EEHighR9       ;
	float smearing_sigma_EELowR9        ;
	float smearing_sigma_error_EBHighR9 ;
	float smearing_sigma_error_EBLowR9  ;
	float smearing_sigma_error_EEHighR9 ;
	float smearing_sigma_error_EELowR9  ;

	std::vector<int> pho_acc;
	std::vector<int> pho_presel;
	std::vector<int> pho_presel_lead;
	std::vector<float> pho_et;
	// Other options
	bool runStatAnalysis;
	bool reRunCiC;
        TString puHist, puMap;//name of pileup reweighting histogram
	bool applyPtoverM;
	float leadEtCut;
	float subleadEtCut;
	std::string massResolutionFileName;

	enum BkgCategory{promptprompt,promptfake,fakefake};
	bool keepPP, keepPF, keepFF;

protected:
	void PreselectPhotons(LoopAll& l, int jentry);
	void loadPuMap(const char * fname, TDirectory * dir);
	void loadPuWeights(int typid, TDirectory * dir);

	std::string name_;
	
	// Vertex analysis
	HggVertexAnalyzer vtxAna_;
	HggVertexFromConversions vtxConv_;

	// Mass Resolution
	MassResolution *massResolutionCalculator;	
	std::map<int, vector<double> > weights;
	int trigCounter_;

	EnergySmearer *eScaleDataSmearer ; // corrections for energy scale data
	EnergySmearer *eScaleSmearer;      // corrections for energy scale  MC
	EnergySmearer *eCorrSmearer;      // corrections for energy scale  MC
	std::vector<float> corrected_pho_energy;
	
	
};

#endif
