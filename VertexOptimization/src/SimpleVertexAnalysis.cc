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

    minBiasRefName = "aux/hMinBiasAlpha.root";
    storeNVert = 10;
    runOnReducedNtuples=true;
    
    sigmaT = 0.010;

}

// ----------------------------------------------------------------------------------------------------
SimpleVertexAnalysis::~SimpleVertexAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void SimpleVertexAnalysis::Term(LoopAll& l) 
{
    //    l.outputFile->cd();
    //hMinBias_->Write();
    //hHiggs_->Write();
    //hMinBias_->SetDirectory(0);
    //hHiggs_->SetDirectory(0);
}

// ----------------------------------------------------------------------------------------------------
void SimpleVertexAnalysis::Init(LoopAll& l) 
{
    doSystematics = false;

    // per vertex tree
    pho1_=0, pho2_=0, dipho_=0;
    l.InitTrees("vtxTree");
    l.BookExternalTreeBranch("nVert",   &nVert_   , "vtxTree");
    l.BookExternalTreeBranch("nPU",     &nPU_     , "vtxTree");
    l.BookExternalTreeBranch("itype",&itype_, "vtxTree");
    l.BookExternalTreeBranch("evweight",&evWeight_, "vtxTree");
    l.BookExternalTreeBranch("isClosestToGen",&isClosestToGen_, "vtxTree");
    l.BookExternalTreeBranch("passCiC",&passCiC_, "vtxTree");
    l.BookExternalTreeBranch("pho1",&pho1_,32000,0, "vtxTree");
    l.BookExternalTreeBranch("pho2",&pho2_,32000,0, "vtxTree");
    l.BookExternalTreeBranch("dipho",&dipho_,32000,0, "vtxTree");
    l.BookExternalTreeBranch("alpha02",&alpha02_, "vtxTree");
    l.BookExternalTreeBranch("alpha03",&alpha03_, "vtxTree");
    l.BookExternalTreeBranch("alpha04",&alpha04_, "vtxTree");
    l.BookExternalTreeBranch("alpha05",&alpha05_, "vtxTree");
    l.BookExternalTreeBranch("ksprob02",&ksprob_[0], "vtxTree");
    l.BookExternalTreeBranch("ksprob03",&ksprob_[1], "vtxTree");
    l.BookExternalTreeBranch("ksprob04",&ksprob_[2], "vtxTree");
    l.BookExternalTreeBranch("ksprob05",&ksprob_[3], "vtxTree");
    l.BookExternalTreeBranch("dZTrue",&dZTrue_, "vtxTree");
    l.BookExternalTreeBranch("vtxZ",&vtxZ_, "vtxTree");
    // time related vars.
    l.BookExternalTreeBranch("dtof1",&dtof1_, "vtxTree");
    l.BookExternalTreeBranch("dtof2",&dtof2_, "vtxTree");
    l.BookExternalTreeBranch("tof1",&tof1_, "vtxTree");
    l.BookExternalTreeBranch("tof2",&tof2_, "vtxTree");
    // emulated reco time 
    l.BookExternalTreeBranch("tReco1",&tReco1_, "vtxTree");
    l.BookExternalTreeBranch("tReco2",&tReco2_, "vtxTree");

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




    // min bias ref histograms
    if( minBiasRefName != "" ) {
	TFile * mbRef = TFile::Open(minBiasRefName);
	hMinBiasRef_[0] = (TH1*)mbRef->Get("alpha02_bkg")->Clone("hMinBiasRef02");
	hMinBiasRef_[1] = (TH1*)mbRef->Get("alpha03_bkg")->Clone("hMinBiasRef03");
	hMinBiasRef_[2] = (TH1*)mbRef->Get("alpha04_bkg")->Clone("hMinBiasRef04");
	hMinBiasRef_[3] = (TH1*)mbRef->Get("alpha05_bkg")->Clone("hMinBiasRef05");
	hMinBiasRef_[0]->SetDirectory(0);
	hMinBiasRef_[1]->SetDirectory(0);
	hMinBiasRef_[2]->SetDirectory(0);
	hMinBiasRef_[3]->SetDirectory(0);
	mbRef->Close();
    } else {
	hMinBiasRef_[0] = 0;
	hMinBiasRef_[1] = 0;
	hMinBiasRef_[2] = 0;
	hMinBiasRef_[3] = 0;
    }
        
    hMinBias_ = new TH1F("hMinBias","hMinBias",200,0,20);
    hMinBias_->SetDirectory(0);
    hHiggs_   = new TH1F("hHiggs","hHiggs",200,0,20);
    hHiggs_->SetDirectory(0); 
        
    StatAnalysis::Init(l);
    //PhotonAnalysis::Init(l);

    rnd = new TRandom();


}


// ----------------------------------------------------------------------------------------------------                                                                                                                                     
void SimpleVertexAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
        PhotonAnalysis::ReducedOutputTree(l,0);
}

// ----------------------------------------------------------------------------------------------------                                                                                                                                     
void SimpleVertexAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s )
{
    if (runOnReducedNtuples){
	vtxAna_.setBranchAdresses(t,"vtx_std_");
	vtxAna_.getBranches(t,"vtx_std_",s);
    }
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
    
    // --- Find reco photons matched to gen photons coming from Higgs
    int higgsind=-1;
    int mc1=-1;
    int mc2=-1;
    int i1=-1;
    int i2=-1;
    diphoton_id = -1;    


    vector<int> rankedVtxs;
    rankedVtxs.clear();

    if (runOnReducedNtuples){
	// if running on reduced ntuple can use directly this
	i1=l.gh_gen2reco1;
	i2=l.gh_gen2reco2;
	
	if (i1>-1 && i2>-1){
	    for(int idipho = 0; idipho < l.dipho_n; ++idipho ) {
		//		cout << idipho << "  " << l.dipho_leadind[idipho] << "  " << l.dipho_subleadind[idipho] << endl;
		if(  (l.dipho_leadind[idipho] == i1 && l.dipho_subleadind[idipho]==i2) ||
		     (l.dipho_leadind[idipho] == i2 && l.dipho_subleadind[idipho]==i1)  )
		    diphoton_id = idipho;
	    }
	    
	}
	if (diphoton_id>-1)  
	    rankedVtxs = (*l.vtx_std_ranked_list)[diphoton_id];
    }
    
    else {
	
	l.FindMCHiggsPhotons( higgsind,  mc1,  mc2,  i1,  i2 );
        
	if (i1>-1 && i2>-1){
	    diphoton_id = 0;
	    
	    // --- photons 
	    PhotonInfo pho1=l.fillPhotonInfos(i1,vtxAlgoParams.useAllConversions,0); // 0 = use default energy stored in pho_p4->Energy()
	    PhotonInfo pho2=l.fillPhotonInfos(i2,vtxAlgoParams.useAllConversions,0);
	    
	    
	    // --- Build vtx quantities
	    vtxAna_.clear();
	    l.vertexAnalysis(vtxAna_, pho1, pho2 );
	    rankedVtxs = l.vertexSelection(vtxAna_, vtxConv_, pho1, pho2, vtxVarNames_, mvaVertexSelection, tmvaPerVtxReader_, tmvaPerVtxMethod);
	    
	}
    }
    
    if ( diphoton_id > -1) {

	// --- Check leading/subleading photon	
	TLorentzVector lead_p4    = l.get_pho_p4( i1, rankedVtxs[0]); // use default photon energy (could also pass for example get_pho_p4( i1, rankedVtxs[0],&smeared_pho_energy[0]), if I have smeared energies.
	TLorentzVector sublead_p4 = l.get_pho_p4( i2, rankedVtxs[0]);
	
	int leadind    = i1;
	int subleadind = i2;
	
	if(sublead_p4.Pt()  > lead_p4.Pt() ) {
	    leadind    = i2;
	    subleadind = i1;
	}
	
	
	// --- Fill trees
	vtxAna_.setPairID(diphoton_id);
	nVert_    = l.vtx_std_n;
	nPU_      = l.pu_n;
	evWeight_ = evweight;
	
    
	// vertex tree
	for(int vi=0; vi<l.vtx_std_n; ++vi) {

	    // recompute p4 for this vertex
	    TLorentzVector photon1 = l.get_pho_p4(leadind, vi); 
	    TLorentzVector photon2 = l.get_pho_p4(subleadind, vi); 
	    TLorentzVector Higgs   = photon1+photon2;
		
	    pho1_     = &photon1;
	    pho2_     = &photon2;
	    dipho_    = &Higgs;

	    // vertex vars	
	    isClosestToGen_ = (vi == closest_id);
	    dZTrue_   = ( *(TVector3*)l.vtx_std_xyz->At(vi) - *genVtx).Z();

	    for( size_t ivar=0; ivar<vtxVarNames_.size(); ++ivar ) {
		vtxVars_[ivar] = (vtxAna_.*(varMeths_[ivar]))(vi);
	    }
	

	    // t.o.f
	    float c = 29.9792458; // cm/ns
	    TVector3 *sc1 = (TVector3*)l.sc_xyz->At(l.pho_scind[leadind]);
	    TVector3 *sc2 = (TVector3*)l.sc_xyz->At(l.pho_scind[subleadind]);
	    TVector3 *vtx = (TVector3*)l.vtx_std_xyz->At(vi);
	    float tof01 = (1./c)*sqrt(pow(sc1->X()-0.,2) + pow(sc1->Y()-0.,2) + pow(sc1->Z()-0.,2)); // wrt nominal IP (0,0,0)
	    float tof02 = (1./c)*sqrt(pow(sc2->X()-0.,2) + pow(sc2->Y()-0.,2) + pow(sc2->Z()-0.,2)); // wrt nominal IP (0,0,0)
	    tof1_  = (1./c)*sqrt(pow(sc1->X()-vtx->X(),2) + pow(sc1->Y()-vtx->Y(),2) + pow(sc1->Z()-vtx->Z(),2));
	    tof2_  = (1./c)*sqrt(pow(sc2->X()-vtx->X(),2) + pow(sc2->Y()-vtx->Y(),2) + pow(sc2->Z()-vtx->Z(),2));

	    // fill dtof for each vertex hypothesis
	    dtof1_ = tof1_-tof01; 
	    dtof2_ = tof2_-tof02;
	    // emulation of Treco : Treco = Tcoll + dtof (I don't care about Tcoll as I always take time differences between two objects )
	    float tof1_truevtx  = (1./c)*sqrt(pow(sc1->X()-genVtx->X(),2) + pow(sc1->Y()-genVtx->Y(),2) + pow(sc1->Z()-genVtx->Z(),2));
	    float tof2_truevtx  = (1./c)*sqrt(pow(sc2->X()-genVtx->X(),2) + pow(sc2->Y()-genVtx->Y(),2) + pow(sc2->Z()-genVtx->Z(),2));
	    tReco1_ = (tof1_truevtx-tof01) + rnd->Gaus(0.,sigmaT);
	    tReco2_ = (tof2_truevtx-tof02) + rnd->Gaus(0.,sigmaT);
	    
	    vtxZ_  = vtx->Z(); 

	    if (!runOnReducedNtuples){

		TH1 * h[4];
		for (int ih = 0 ; ih<4; ih++){
		    h[ih]= ( hMinBiasRef_[ih]!=0 ? (TH1*)hMinBiasRef_[ih]->Clone("h") : 0);
		    if( h[ih] ) { h[ih]->Reset("ICE"); }
		}
		
		alpha02_.clear();
		alpha03_.clear();
		alpha04_.clear();
		alpha05_.clear();
		// loop over tracks
		int ntks = l.vtx_std_ntks[vi];
		for(int ti=0; ti<ntks; ++ti) {
		    float tmp02 =0;
		    float tmp03 =0;
		    float tmp04 =0;
		    float tmp05 =0;
		    
		    int tkind = (*l.vtx_std_tkind)[vi][ti];
		    if( tkind >= l.tk_n ) { continue; }
		    TLorentzVector* tk1p4 = (TLorentzVector*)l.tk_p4->At(tkind);
		    float eta1 = tk1p4->Eta();
		    float phi1 = tk1p4->Phi();
		    float pt1  = tk1p4->Pt();
		    // loop over tracks i+1...N
		    for(int tj=ti+1; tj<ntks; ++tj) {
			int tkind2 = (*l.vtx_std_tkind)[vi][tj];
			if( tkind2 >= l.tk_n ) { continue; }
			TLorentzVector* tk2p4 = (TLorentzVector*)l.tk_p4->At(tkind2);
			float eta2 = tk2p4->Eta();
			float phi2 = tk2p4->Phi();
			float pt2  = tk2p4->Pt();
			float deltaR = sqrt( pow(eta1-eta2,2) + pow(phi1-phi2,2) );
			//cout << "deltaR = " <<deltaR <<endl;
			if (deltaR<0.2)
			    tmp02+=pt2/deltaR;
			if (deltaR<0.3)
			    tmp03+=pt2/deltaR;
			if (deltaR<0.4)
			    tmp04+=pt2/deltaR;
			if (deltaR<0.5)
			    tmp05+=pt2/deltaR;
		    }
	    
		    // fill vector of alpha values for this vtx (one value for each track in the vtx)
		    alpha02_.push_back(tmp02);
		    alpha03_.push_back(tmp03);
		    alpha04_.push_back(tmp04);
		    alpha05_.push_back(tmp05);
		    
		    if( ! isClosestToGen_ ) {
			hMinBias_->Fill(tmp03);
		    } else {
			hHiggs_->Fill(tmp03);
		    }
		    
		    if( h[0] )  h[0]->Fill(tmp02);
		    if( h[1] )  h[1]->Fill(tmp03);
		    if( h[2] )  h[2]->Fill(tmp04);
		    if( h[3] )  h[3]->Fill(tmp05);
		    
		}// end loop over vertex tracks
	    

		// compute Kolmogorov probability
		for (int ih=0; ih<4;ih++){
		    if( h[ih] )  {
			ksprob_[ih] = hMinBiasRef_[ih]->KolmogorovTest(h[ih]);
			delete h[ih];
		    } else {
			ksprob_[ih] = 0.;
		    }
		}
	    }
	    
	    l.FillTreeContainer("vtxTree");
	}// end loop over vertices
	
    

	// per event tree (with info from the first N ranked vertices)
	//	vector<int> & rankedVtxs = (*l.vtx_std_ranked_list)[diphoton_id];
	MVA_.assign(storeNVert,-10);
	dZ_.assign(storeNVert,-100);
	diphoM_.assign(storeNVert,-2);
	diphoPt_.assign(storeNVert,-2);
	
	TLorentzVector truevtx_lead_pho = l.get_pho_p4( leadind, genVtx);
	TLorentzVector truevtx_sublead_pho = l.get_pho_p4( subleadind, genVtx);
	TLorentzVector truevtx_dipho = truevtx_lead_pho+truevtx_sublead_pho;
	mTrue_    = gP4.M(); 
	mTrueVtx_ = truevtx_dipho.M(); 
	category_ = category;
	nConv_    = vtxAna_.nconv(0);
	
	zTrue_    = genVtx->Z();	    
	dZTrue_   = ( *(TVector3*)l.vtx_std_xyz->At(rankedVtxs[0]) - *genVtx).Z();
	Double_t vtxProbInputs[2];
	float wtot = 0.;
	zRMS_ = 0.;
	for (size_t vi=0;vi<rankedVtxs.size();vi++) {
	    if(vi>=storeNVert) break;
	    MVA_[vi] = vtxAna_.mva(rankedVtxs[vi]);
	    dZ_[vi]  = ( *(TVector3*)l.vtx_std_xyz->At(rankedVtxs[vi]) - *(TVector3*)l.vtx_std_xyz->At(rankedVtxs[0])).Z();
	    
	    TLorentzVector lead_pho = l.get_pho_p4( leadind, rankedVtxs[vi]);
	    TLorentzVector sublead_pho = l.get_pho_p4( subleadind, rankedVtxs[vi]);
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
