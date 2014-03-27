import ROOT
import sys
import os

from ROOT import *

# ROOT Setup
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetOptTitle(0)

outdir = sys.argv[1]
os.mkdir(outdir)

files = ["~/eos/cms/store/cmst3/user/malberti/HIGGS/Upgrade/vtxPU140_testTiming_10ps_v1/histograms_vertexOpt.root",
         "~/eos/cms/store/cmst3/user/malberti/HIGGS/Upgrade/vtxPU2012_testTiming_10ps_v1/histograms_vertexOpt.root"]
sqrts = ["14TeV","8TeV"]

procs = ["ggh_m125","vbf_m125","wzh_m125"]

vars = ["logsumpt2",
        "ptbal",
        "ptasym",
        "limpulltoconv",
        "nconv",
        "dt",
#        "vtx_mva"
        ]

hsig = {}
hbkg = {}
canvas = {}
for proc in procs:
    hsig[proc]={}
    hbkg[proc]={}
    canvas[proc]={}

canvasEffPt = {}
canvasEffNvtx = {}


    
hpt_mva        = {}
hpt_mva_rv     = {}
hpt_sumpt2     = {}
hpt_sumpt2_rv  = {}

hnvtx_mva        = {}
hnvtx_mva_rv     = {}
hnvtx_sumpt2     = {}
hnvtx_sumpt2_rv  = {}


legend = TLegend(0.52,0.68,0.89,0.89)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetTextFont(42)

legend2 = TLegend(0.50,0.12,0.89,0.35)
legend2.SetFillStyle(0)
legend2.SetBorderSize(0)
legend2.SetTextFont(42)


for i,file in enumerate(files):
    f=TFile.Open(file) 
    for p,proc in enumerate(procs):
        #
        for v,var in enumerate(vars):
            # open new canvas
            if (i==0):
                canvas[proc][var] = TCanvas(proc.replace("_m125","")+"_"+var, proc.replace("_m125","")+"_"+var,500,500)
            # variables histograms
            hsig[proc][var] = f.Get(var+"_rv_cat0_"+proc+"_"+sqrts[i])
            hbkg[proc][var] = f.Get(var+"_wv_cat0_"+proc+"_"+sqrts[i])

            hsig[proc][var].SetLineColor(ROOT.kBlue)
            hbkg[proc][var].SetLineColor(ROOT.kRed)
            ymax = max(hsig[proc][var].GetMaximum()/hsig[proc][var].GetEntries(), hbkg[proc][var].GetMaximum()/hbkg[proc][var].GetEntries())*hsig[proc][var].GetEntries()
            hsig[proc][var].SetMaximum(ymax*1.2)
            
            canvas[proc][var].cd()
            if (i==0):
                hsig[proc][var].SetTitle(proc.replace("_m125","")+"_"+var)
                hsig[proc][var].GetXaxis().SetTitle(var)
                hsig[proc][var].GetYaxis().SetTitle("a.u.")
                hsig[proc][var].DrawNormalized("histo")
                hsig[proc][var].DrawNormalized("histo")
                hbkg[proc][var].DrawNormalized("histo same")
                if v==0 and p==0:
                    legend.AddEntry(hsig[proc][var],"primary vertex, <PU>=140","L")
                    legend.AddEntry(hbkg[proc][var],"pile-up vertex, <PU>=140","L")
            else:
                hsig[proc][var].SetLineStyle(2)
                hbkg[proc][var].SetLineStyle(2)
                hsig[proc][var].DrawNormalized("histo same")
                hbkg[proc][var].DrawNormalized("histo same")
                if v==0 and p==0:
                    legend.AddEntry(hsig[proc][var],"primary vertex, <PU>=20","L")
                    legend.AddEntry(hbkg[proc][var],"pile-up vertex, <PU>=20","L")
            legend.Draw("same")        

            
    # eff plots vs pt
    
        hpt_mva[proc]       = f.Get("pt_cat0_"+proc+"_"+sqrts[i])
        hpt_mva_rv[proc]    = f.Get("pt_rv_cat0_"+proc+"_"+sqrts[i])
        hpt_sumpt2[proc]    = f.Get("pt_highestsumpt2_cat0_"+proc+"_"+sqrts[i])
        hpt_sumpt2_rv[proc] = f.Get("pt_highestsumpt2_rv_cat0_"+proc+"_"+sqrts[i])
        hpt_mva[proc].Sumw2()
        hpt_mva_rv[proc].Sumw2()
        hpt_sumpt2[proc].Sumw2()
        hpt_sumpt2_rv[proc].Sumw2()
        hpt_mva_rv[proc].Divide(hpt_mva_rv[proc],hpt_mva[proc],1,1,"B")
        hpt_sumpt2_rv[proc].Divide(hpt_sumpt2_rv[proc],hpt_sumpt2[proc],1,1,"B")
        # eff plots vs nvtx
        hnvtx_mva[proc]       = f.Get("nvtx_cat0_"+proc+"_"+sqrts[i])
        hnvtx_mva_rv[proc]    = f.Get("nvtx_rv_cat0_"+proc+"_"+sqrts[i])
        hnvtx_sumpt2[proc]    = f.Get("nvtx_highestsumpt2_cat0_"+proc+"_"+sqrts[i])
        hnvtx_sumpt2_rv[proc] = f.Get("nvtx_highestsumpt2_rv_cat0_"+proc+"_"+sqrts[i])
        hnvtx_mva[proc].Sumw2()
        hnvtx_mva_rv[proc].Sumw2()
        hnvtx_sumpt2[proc].Sumw2()
        hnvtx_sumpt2_rv[proc].Sumw2()
        hnvtx_mva_rv[proc].Divide(hnvtx_mva_rv[proc],hnvtx_mva[proc],1,1,"B")
        hnvtx_sumpt2_rv[proc].Divide(hnvtx_sumpt2_rv[proc],hnvtx_sumpt2[proc],1,1,"B")
        
        if i==0:
            canvasEffPt[proc] = TCanvas("efficiency_vs_pt_"+proc.replace("_m125",""), "efficiency_vs_pt_"+proc.replace("_m125",""),500,500)    
            hpt_mva_rv[proc].SetMarkerStyle(20)
            hpt_sumpt2_rv[proc].SetMarkerStyle(24)
            hpt_mva_rv[proc].SetMarkerColor(ROOT.kRed)
            hpt_sumpt2_rv[proc].SetMarkerColor(ROOT.kRed)
            hpt_mva_rv[proc].SetLineColor(ROOT.kRed)
            hpt_sumpt2_rv[proc].SetLineColor(ROOT.kRed)
            hpt_mva_rv[proc].SetTitle("efficiency_vs_pt_"+proc)
            hpt_mva_rv[proc].GetYaxis().SetRangeUser(0,1.1)
            hpt_mva_rv[proc].GetYaxis().SetTitle("fraction |z_{reco}-z_{true}| < 1 cm")
            hpt_mva_rv[proc].GetXaxis().SetTitle("pT (GeV)")
            hpt_mva_rv[proc].Draw()
            hpt_sumpt2_rv[proc].Draw("same")

            canvasEffNvtx[proc] = TCanvas("efficiency_vs_nvtx_"+proc.replace("_m125",""), "efficiency_vs_nvtx_"+proc.replace("_m125",""),500,500)    
            hnvtx_mva_rv[proc].SetMarkerStyle(20)
            hnvtx_sumpt2_rv[proc].SetMarkerStyle(24)
            hnvtx_mva_rv[proc].SetMarkerColor(ROOT.kRed)
            hnvtx_sumpt2_rv[proc].SetMarkerColor(ROOT.kRed)
            hnvtx_mva_rv[proc].SetLineColor(ROOT.kRed)
            hnvtx_sumpt2_rv[proc].SetLineColor(ROOT.kRed)
            hnvtx_mva_rv[proc].SetTitle("efficiency_vs_pt_"+proc)
            hnvtx_mva_rv[proc].GetYaxis().SetRangeUser(0,1.1)
            hnvtx_mva_rv[proc].GetYaxis().SetTitle("fraction |z_{reco}-z_{true}| < 1 cm")
            hnvtx_mva_rv[proc].GetXaxis().SetTitle("number of reconstructed vertices")
            hnvtx_mva_rv[proc].Draw()
            hnvtx_sumpt2_rv[proc].Draw("same")

            if p==0:
                legend2.AddEntry(hpt_mva_rv[proc],"<PU>=140, 14 TeV, vtxMVA","P")
                legend2.AddEntry(hpt_sumpt2_rv[proc],"<PU>=140, 14 TeV, max(sumpt2)","P")
        else:
            canvasEffPt[proc].cd() 
            hpt_mva_rv[proc].SetMarkerStyle(20)
            hpt_sumpt2_rv[proc].SetMarkerStyle(24)
            hpt_mva_rv[proc].SetMarkerColor(ROOT.kBlue)
            hpt_sumpt2_rv[proc].SetMarkerColor(ROOT.kBlue)
            hpt_mva_rv[proc].Draw("same")
            hpt_sumpt2_rv[proc].Draw("same")

            canvasEffNvtx[proc].cd() 
            hnvtx_mva_rv[proc].SetMarkerStyle(20)
            hnvtx_sumpt2_rv[proc].SetMarkerStyle(24)
            hnvtx_mva_rv[proc].SetMarkerColor(ROOT.kBlue)
            hnvtx_sumpt2_rv[proc].SetMarkerColor(ROOT.kBlue)
            hnvtx_mva_rv[proc].Draw("same")
            hnvtx_sumpt2_rv[proc].Draw("same")

            if p==0:
                legend2.AddEntry(hpt_mva_rv[proc],"<PU>=20, 8 TeV, vtxMVA","P")
                legend2.AddEntry(hpt_sumpt2_rv[proc],"<PU>=20, 8 TeV, max(sumpt2)","P")

        canvasEffPt[proc].cd()        
        legend2.Draw("same")
        canvasEffNvtx[proc].cd()        
        legend2.Draw("same")        


for proc in procs:
    for typ in ".pdf",".png",".C":
        canvasEffPt[proc].SaveAs(outdir+"/"+canvasEffPt[proc].GetName()+typ)
        canvasEffNvtx[proc].SaveAs(outdir+"/"+canvasEffNvtx[proc].GetName()+typ)
    for var in vars:
        for typ in ".pdf",".png",".C":
            canvas[proc][var].SaveAs(outdir+"/"+canvas[proc][var].GetName()+typ)

raw_input("ok?")
