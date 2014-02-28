#include "../interface/SimpleVertexAnalysis.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
SimpleVertexAnalysis::SimpleVertexAnalysis()  
{
    name_ = "SimpleVertexAnalysis";

    minBiasRefName = "aux/minBiasRef.root";
    storeNVert = 10;
}

// ----------------------------------------------------------------------------------------------------
SimpleVertexAnalysis::~SimpleVertexAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void SimpleVertexAnalysis::Term(LoopAll& l) 
{
    //// hMinBiasSpecturm_->Write();
    //// hHiggsSpecturm_->Write();
    //// hMinBiasSpecturm_->SetDirectory(0);
    //// hHiggsSpecturm_->SetDirectory(0);
}

// ----------------------------------------------------------------------------------------------------
void SimpleVertexAnalysis::Init(LoopAll& l) 
{
    doSystematics = false;

    if( minBiasRefName != "" ) {
	TFile * mbRef = TFile::Open(minBiasRefName);
	hMinBiasRef_ = (TH1*)mbRef->Get("minBiasSpecturm")->Clone("hMinBiasRef");
	hMinBiasRef_->SetDirectory(0);
	mbRef->Close();
    } else {
	hMinBiasRef_ = 0;
    }
        
    // per vertex tree
    pho1_=0, pho2_=0, dipho_=0;
    l.InitTrees("vtxTree");
    l.BookExternalTreeBranch("nVert",   &nVert_   , "vtxTree");
    l.BookExternalTreeBranch("nPU",     &nPU_     , "vtxTree");
    l.BookExternalTreeBranch("evweight",&evWeight_, "vtxTree");
    l.BookExternalTreeBranch("ksprob",&ksprob_, "vtxTree");
    l.BookExternalTreeBranch("isClosestToGen",&isClosestToGen_, "vtxTree");
    l.BookExternalTreeBranch("passCiC",&passCiC_, "vtxTree");
    l.BookExternalTreeBranch("pho1",&pho1_,32000,0, "vtxTree");
    l.BookExternalTreeBranch("pho2",&pho2_,32000,0, "vtxTree");
    l.BookExternalTreeBranch("dipho",&dipho_,32000,0, "vtxTree");
    l.BookExternalTreeBranch("itype",&itype_, "vtxTree");
    
    vtxVarNames_.push_back("ptvtx"), vtxVarNames_.push_back("ptasym"), vtxVarNames_.push_back("ptratio"), 
	vtxVarNames_.push_back("ptbal"), vtxVarNames_.push_back("logsumpt2"), vtxVarNames_.push_back("ptmax3"), 
	vtxVarNames_.push_back("ptmax"), vtxVarNames_.push_back("nchthr"), vtxVarNames_.push_back("sumtwd"),
	vtxVarNames_.push_back("pulltoconv"), vtxVarNames_.push_back("limpulltoconv"), 
	vtxVarNames_.push_back("nch"), vtxVarNames_.push_back("nconv"), vtxVarNames_.push_back("nlegs"), 
	vtxVarNames_.push_back("mva"); 
    vtxVars_.resize( vtxVarNames_.size() );
    for( size_t iv=0; iv<vtxVarNames_.size(); ++iv ) {
	l.BookExternalTreeBranch(vtxVarNames_[iv].c_str(),&vtxVars_[iv],"vtxTree");
    }

    hMinBiasSpecturm_ = new TH1F("minBiasSpecturm","minBiasSpecturm;;p_{T} (GeV/c)",200,0,20);
    hMinBiasSpecturm_->SetDirectory(0);
    hHiggsSpecturm_   = new TH1F("higgsSpecturm","higgsSpecturm;;p_{T} (GeV/c)",200,0,20);
    hHiggsSpecturm_->SetDirectory(0); 
    
    // per event tree
    MVA_.resize(storeNVert);
    dZ_.resize(storeNVert);
    diphoM_.resize(storeNVert);
    diphoPt_.resize(storeNVert);
    
    l.InitTrees("evtTree");
    l.BookExternalTreeBranch("itype",&itype_, "evtTree");
    l.BookExternalTreeBranch("dZTrue",&dZTrue_, "evtTree");
    l.BookExternalTreeBranch("zTrue",&zTrue_, "evtTree");
    l.BookExternalTreeBranch("zRMS",&zRMS_, "evtTree");
    l.BookExternalTreeBranch("category",&category_,"category/I", "evtTree");
    l.BookExternalTreeBranch("mTrue",&mTrue_, "evtTree");
    l.BookExternalTreeBranch("mTrueVtx",&mTrueVtx_, "evtTree");
    l.BookExternalTreeBranch("nVert",&nVert_, "evtTree");
    l.BookExternalTreeBranch("nPU",&nPU_, "evtTree");
    l.BookExternalTreeBranch("evWeight",&evWeight_, "evtTree");
    l.BookExternalTreeBranch("nConv",&nConv_, "evtTree");
    for(int i=0;i<storeNVert; i++){
    	l.BookExternalTreeBranch(Form("MVA%d",i),&MVA_[i], "evtTree");
    	l.BookExternalTreeBranch(Form("dZ%d",i),&dZ_[i], "evtTree");
	l.BookExternalTreeBranch(Form("diphoM%d",i),&diphoM_[i], "evtTree");
	l.BookExternalTreeBranch(Form("diphoPt%d",i),&diphoPt_[i], "evtTree");
    }
    
    StatAnalysis::Init(l);
}

void SimpleVertexAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
    PhotonAnalysis::ReducedOutputTree(l,0);
}

// ----------------------------------------------------------------------------------------------------
bool SimpleVertexAnalysis::AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, 
					      float & mass, float & evweight, int & category, int & diphoton_id,
					      bool & isCorrectVertex, float &kinematic_bdtout,
					      bool isSyst, 
					      float syst_shift, bool skipSelection,
					      BaseGenLevelSmearer *genSys, 
					      BaseSmearer *phoSys, BaseDiPhotonSmearer * diPhoSys)
{
    // book vertex variables
    static std::vector<HggVertexAnalyzer::getter_t> varMeths_(0);
    if( varMeths_.empty() ) {
	for( size_t iv=0; iv<vtxVarNames_.size(); ++iv ){
	    varMeths_.push_back(HggVertexAnalyzer::dictionary()[vtxVarNames_[iv]].first);
	}
    }

    int cur_type = l.itype[l.current];
    itype_ = cur_type;
    float sampleweight = l.sampleContainer[l.current_sample_index].weight();
    evweight = sampleweight;


    // --- Find reco photons matched to gen photons coming from Higgs
    diphoton_id = -1;    
    
    int i1=l.gh_gen2reco1;
    int i2=l.gh_gen2reco2;
    
    cout << i1 << "  " << i2 << endl;
    
    if (i1>-1 && i2>-1){
	for(int idipho = 0; idipho < l.dipho_n; ++idipho ) {
	    cout << idipho << "  " << l.dipho_leadind[idipho] << "  " << l.dipho_subleadind[idipho] << endl;
	    if(  (l.dipho_leadind[idipho] == i1 && l.dipho_subleadind[idipho]==i2) ||
		 (l.dipho_leadind[idipho] == i2 && l.dipho_subleadind[idipho]==i1)  )
		diphoton_id = idipho;
	}
    }
    cout << "dipho_id = " << diphoton_id <<endl;

    // --- Find reco vtx closest to the true vertex
    TVector3 * genVtx = (TVector3*)l.gv_pos->At(0);
    int closest_id = -1;
    float minDist = 999.;
    for(int vi=0; vi<l.vtx_std_n; ++vi) {
	float dist =  fabs( ( *((TVector3*)l.vtx_std_xyz->At(vi)) - *genVtx ).Z() );
	if( dist < minDist ) {
	    closest_id = vi;
	    minDist = dist;
	}
    }
    	

    TLorentzVector photon1;
    TLorentzVector photon2;
    TLorentzVector Higgs;


    // --- Fill trees
    if (diphoton_id>-1){
	vtxAna_.setPairID(diphoton_id);
	nVert_    = l.vtx_std_n;
	nPU_      = l.pu_n;
	evWeight_ = evweight;
	
	
	// vertex tree
	for(int vi=0; vi<l.vtx_std_n; ++vi) {
	    // recompute p4 for this vertex
	    photon1 = l.get_pho_p4(i1, vi, 0); 
	    photon2 = l.get_pho_p4(i2, vi, 0); 
	    Higgs   = photon1+photon2;

	    pho1_     = &photon1;
	    pho2_     = &photon2;
	    dipho_    = &Higgs;

	    isClosestToGen_ = (vi == closest_id);

	    TH1 * h = ( hMinBiasRef_!=0 ? (TH1*)hMinBiasRef_->Clone("h") : 0);
	    if( h ) { h->Reset("ICE"); }
	    
	    for( size_t ivar=0; ivar<vtxVarNames_.size(); ++ivar ) {
		vtxVars_[ivar] = (vtxAna_.*(varMeths_[ivar]))(vi);
	    }
	    
	    int ntks = l.vtx_std_ntks[vi];
	    for(int ti=0; ti<ntks; ++ti) {
		int tkind = (*l.vtx_std_tkind)[vi][ti];
		if( tkind >= l.tk_n ) { continue; }
		TLorentzVector* tkp4 = (TLorentzVector*)l.tk_p4->At(tkind);
		if( ! isClosestToGen_ ) {
		    hMinBiasSpecturm_->Fill(tkp4->Pt());
		} else {
		    hHiggsSpecturm_->Fill(tkp4->Pt());
		}
		if( h ) { h->Fill(tkp4->Pt()); }
	    }
	    
	    if( h )  {
		ksprob_ = hMinBiasRef_->KolmogorovTest(h);
		delete h;
	    } else {
		ksprob_ = 0.;
	    }
	    
	    l.FillTreeContainer("vtxTree");
	}


	// per event tree (with info from the first N ranked vertices
	vector<int> & rankedVtxs = (*l.vtx_std_ranked_list)[diphoton_id];
	MVA_.assign(storeNVert,-10);
	dZ_.assign(storeNVert,-100);
	diphoM_.assign(storeNVert,-2);
	diphoPt_.assign(storeNVert,-2);

	TLorentzVector truevtx_lead_pho = l.get_pho_p4( l.dipho_leadind[diphoton_id], genVtx, 0);
	TLorentzVector truevtx_sublead_pho = l.get_pho_p4( l.dipho_subleadind[diphoton_id], genVtx, 0);
	TLorentzVector truevtx_dipho = truevtx_lead_pho+truevtx_sublead_pho;
	mTrue_ = gP4.M(); 
	mTrueVtx_ = truevtx_dipho.M(); 
	zTrue_ = genVtx->Z();
	category_ = category;
	nConv_ = vtxAna_.nconv(0);
	    
	dZTrue_ = ( *(TVector3*)l.vtx_std_xyz->At(rankedVtxs[0]) - *genVtx).Z();
	Double_t vtxProbInputs[2];
	float wtot = 0.;
	zRMS_ = 0.;
	for (size_t vi=0;vi<rankedVtxs.size();vi++) {
	    if(vi>=storeNVert) break;
	    MVA_[vi] = vtxAna_.mva(rankedVtxs[vi]);
	    dZ_[vi] = ( *(TVector3*)l.vtx_std_xyz->At(rankedVtxs[vi]) - *(TVector3*)l.vtx_std_xyz->At(rankedVtxs[0])).Z();
	    
	    TLorentzVector lead_pho = l.get_pho_p4( l.dipho_leadind[diphoton_id], rankedVtxs[vi], 0);
	    TLorentzVector sublead_pho = l.get_pho_p4( l.dipho_subleadind[diphoton_id], rankedVtxs[vi], 0);
	    TLorentzVector dipho = lead_pho+sublead_pho;
	    diphoM_[vi]= dipho.M();
	    diphoPt_[vi]= dipho.Pt();
	}
	l.FillTreeContainer("evtTree");
		
	
    }
    
    return false;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
