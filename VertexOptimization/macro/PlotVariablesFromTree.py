import ROOT
import sys
import os

from ROOT import *

# ROOT Setup
ROOT.gROOT.SetStyle("Plain")
#ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetOptTitle(0)


def bookHisto(name, nbins, min, max, col):
    h = TH1F(name,name,nbins,min, max)
    h.SetLineColor(col)
    h.SetLineWidth(2)
    return h

outdir = sys.argv[1]
os.mkdir(outdir)

ptmin=0
ptmax=20

file = TFile.Open("~/eos/cms/store/cmst3/user/malberti/HIGGS/Upgrade/vtxPU140_testTiming_10ps/histograms_vertexOpt.root")
#file = TFile.Open("/afs/cern.ch/work/m/malberti/HGG/Upgrade/CMSSW_6_1_2/src/h2gglobe/AnalysisScripts/vertexOpt.root")

procs = ["ggh_m125_14TeV","vbf_m125_14TeV","wzh_m125_14TeV"]
#procs = ["ggh_m125_14TeV"]

vars = {"nVert":[150,50,200],
        "logsumpt2":[150,-10,20],
        "ptbal":[240,-20,100],
        "ptasym":[100,-1,1],
        "limpulltoconv":[96,0,12],
        "(tReco1-dtof1-(tReco2-dtof2))":[200,-1,1],
        "mva":[200,-1,1],
#        "alpha02":[600,0,300],
#        "alpha03":[600,0,300],
#        "alpha04":[600,0,300],
#        "alpha05":[600,0,300]
        }

hsig = {}
hbkg = {}
hpt  = {}
hptgood  = {}
hnvtx = {}
hnvtxgood = {}

for proc in procs:
    tree = file.Get("vtxTree/%s"%proc)    
    maxentries = tree.GetEntries()
    #maxentries = 10000
    for v in vars:
        sel_sig="isClosestToGen==1 && dipho.Pt() < %d && dipho.Pt() > %d"%(ptmax,ptmin)     
        sel_bkg="isClosestToGen==0 && dipho.Pt() < %d && dipho.Pt() > %d"%(ptmax,ptmin)     
        toplot = v
        if v == "limpulltoconv":
            sel_sig=sel_sig+" && nconv>0"
            sel_bkg=sel_bkg+" && nconv>0"
        #sig
        hname = "%s_sig_%s"%(v,proc)
        if '(' in hname:
            hname = "dtof_sig_%s"%(proc)
        h = bookHisto(hname, vars[v][0],vars[v][1],vars[v][2],ROOT.kGreen+1)
        if proc not in hsig:
            hsig[proc] = [h]
        else:
            hsig[proc].append(h)
            
        tree.Project(hname,toplot,sel_sig,"",maxentries)
        #bkg
        hname = "%s_bkg_%s"%(v,proc)
        if '(' in hname:
            hname = "dtof_bkg_%s"%(proc)
        h = bookHisto(hname, vars[v][0],vars[v][1],vars[v][2],ROOT.kRed)
        if proc not in hbkg:
            hbkg[proc] = [h]
        else:
            hbkg[proc].append(h)
        tree.Project(hname,toplot,sel_bkg,"",maxentries)
    # eff plots vs pt
    tree = file.Get("evtTree/%s"%proc)    
    hname = "pt_good_%s"%proc
    hptgood[proc] = TH1F(hname, hname, 40,0,200)
    tree.Project(hname,"diphoPt0","abs(dZTrue)<1","",maxentries)
    hname = "pt_%s"%proc
    hpt[proc] = TH1F(hname, hname, 40,0,200)
    tree.Project(hname,"diphoPt0","","",maxentries)
    hptgood[proc].Sumw2()
    hpt[proc].Sumw2()    
    hptgood[proc].Divide(hptgood[proc],hpt[proc],1,1,"B")
    # eff plots vs nvtx 
    hname = "nvtx_good_%s"%proc
    hnvtxgood[proc] = TH1F(hname, hname, 200,0,200)
    tree.Project(hname,"nVert","abs(dZTrue)<1","",maxentries)
    hname = "nvtx_%s"%proc
    hnvtx[proc] = TH1F(hname, hname, 200,0,200)
    tree.Project(hname,"nVert","","",maxentries)
    hnvtxgood[proc].Sumw2()
    hnvtx[proc].Sumw2()    
    hnvtxgood[proc].Divide(hnvtxgood[proc],hnvtx[proc],1,1,"B")



#plots
for proc in procs:
    cname = "efficiency_vs_pt_%s"%proc
    ceffPt = TCanvas(cname,cname,500,500)
    hptgood[proc].GetYaxis().SetRangeUser(0,1.1)
    hptgood[proc].SetMarkerStyle(20)
    hptgood[proc].Draw("e")
    for typ in ".pdf",".png",".C":
        ceffPt.SaveAs(outdir+"/"+cname+typ)

    cname = "efficiency_vs_nvtx_%s"%proc
    ceffNvtx = TCanvas(cname,cname,500,500)
    hnvtxgood[proc].GetYaxis().SetRangeUser(0,1.1)
    hnvtxgood[proc].SetMarkerStyle(20)
    hnvtxgood[proc].Draw("e")
    for typ in ".pdf",".png",".C":
        ceffNvtx.SaveAs(outdir+"/"+cname+typ)

    for i,v in enumerate(vars):
        if 'tReco' in v:
            v = 'dt'
        cname = "%s_%s"%(v,proc)
        c = TCanvas(cname,cname,500,500)
        c.cd()
        gStyle.SetOptTitle(0)
        hbkg[proc][i].SetXTitle(v)
        ymax = max(hsig[proc][i].GetMaximum()/hsig[proc][i].GetEntries(), hbkg[proc][i].GetMaximum()/hbkg[proc][i].GetEntries())*hbkg[proc][i].GetEntries()
        hbkg[proc][i].SetMaximum(ymax*1.2)
        hbkg[proc][i].DrawNormalized()
        hsig[proc][i].DrawNormalized("same")
        
        for typ in ".pdf",".png",".C":
            c.SaveAs(outdir+"/"+cname+typ)

raw_input("ok?")
