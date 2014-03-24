#!/bin/env python

from ROOT import TFile, TF1, TCanvas, kBlack, kBlue, kRed, kGreen, kCyan, kMagenta, kWhite, TGraph, TF2, TF12, kFullCircle, gROOT, TLegend, TPaveText, TGraphAsymmErrors
from sys import argv
import math

def eff(fin,name,var="pt",rebin=5):

    num=fin.Get("%s_rv_cat0_tot" % var )
    den=fin.Get("%s_cat0_tot"    % var )

    if rebin > 0:
        num.Rebin(rebin)
        den.Rebin(rebin)


    eff = num.Clone(name)
    eff.Divide( num, den, 1., 1., "B" )

    num.Rebin( num.GetNbinsX() )
    den.Rebin( den.GetNbinsX() )
    avg_eff = num.Clone()
    avg_eff.Divide( num, den, 1., 1., "B" )
    print "%1.3f +- %1.3f" % ( avg_eff.GetBinContent(1), avg_eff.GetBinError(1) )
    return eff

def prob(fin,name,var="pt",profile=True,rebin=5):

    prb = fin.Get("vtxprob_%s_cat0_tot" % var)

    if rebin > 0:
        prb.RebinX(rebin)
    
    if profile:
        return prb.ProfileX(name,1,-1,"w")

    return prb.Clone(name)

gROOT.LoadMacro("../Macros/rootglobestyle.C")
from ROOT import setTDRStyle, gStyle
setTDRStyle()
gStyle.SetOptTitle(0)


file=argv[1]

fileUp=argv[2]
fileDown=argv[3]

nvtxrebin = 0

fin=TFile.Open(file)
eff1 = eff(fin,"eff1",rebin=10)
eff2 = eff(fin,"eff2","nvtx",rebin=nvtxrebin)
print eff2

prob1 = prob(fin,"prob1")
prob2 = prob(fin,"prob2","nvtx",True,nvtxrebin)

#slope up
finUp=TFile.Open(fileUp)
prob1up = prob(finUp,"prob1up")
prob2up = prob(finUp,"prob2up","nvtx",True,nvtxrebin)

prob1up.SetLineColor(2)
prob2up.SetLineColor(2)

#slope down
finDown=TFile.Open(fileDown)
prob1down = prob(finDown,"prob1down")
prob2down = prob(finDown,"prob2down","nvtx",True,nvtxrebin)
prob1down.SetLineColor(6)
prob2down.SetLineColor(6)

#build error band
ge1 = TGraphAsymmErrors();
for bin in range(1,prob1.GetNbinsX()+1):
    ge1.SetPoint(bin-1,prob1.GetBinCenter(bin),prob1.GetBinContent(bin))
    errdown  = abs(prob1up.GetBinContent(bin)-prob1.GetBinContent(bin))
    errup    = abs(prob1down.GetBinContent(bin)-prob1.GetBinContent(bin))
    errstat  = prob1.GetBinError(bin)
    ge1.SetPointEYlow(bin-1,math.sqrt(errdown*errdown+errstat*errstat))
    ge1.SetPointEYhigh(bin-1,math.sqrt(errup*errup+errstat*errstat))
    #ge1.SetPointEYlow(bin-1,abs(prob1up.GetBinContent(bin)-prob1.GetBinContent(bin)))
    #ge1.SetPointEYhigh(bin-1,abs(prob1down.GetBinContent(bin)-prob1.GetBinContent(bin)))    
ge1.SetFillColor(4);
ge1.SetLineColor(4);
ge1.SetMarkerColor(4);

ge2 = TGraphAsymmErrors();
for bin in range(1,prob2.GetNbinsX()+1):
    ge2.SetPoint(bin-1,prob2.GetBinCenter(bin),prob2.GetBinContent(bin))
    errdown  = abs(prob2up.GetBinContent(bin)-prob2.GetBinContent(bin))
    errup    = abs(prob2down.GetBinContent(bin)-prob2.GetBinContent(bin))
    errstat  = prob2.GetBinError(bin)
    ge2.SetPointEYlow(bin-1,math.sqrt(errdown*errdown+errstat*errstat))
    ge2.SetPointEYhigh(bin-1,math.sqrt(errup*errup+errstat*errstat))
    #ge2.SetPointEYlow(bin-1,abs(prob2up.GetBinContent(bin)-prob2.GetBinContent(bin)))
    #ge2.SetPointEYhigh(bin-1,abs(prob2down.GetBinContent(bin)-prob2.GetBinContent(bin)))
ge2.SetFillColor(4);
ge2.SetLineColor(4);
ge2.SetMarkerColor(4);

eff1.SetMarkerStyle(21);
eff2.SetMarkerStyle(21);

eff1.SetLineColor(kBlack);
eff2.SetLineColor(kBlack);

prob1.SetFillColor(4);
prob1.SetLineColor(4);
prob1.SetMarkerColor(4);

prob2.SetFillColor(4);
prob2.SetLineColor(4);
prob2.SetMarkerColor(4);

#pt = TPaveText(51.37856,0.5654967,190.8217,0.6623778,"br")
pt = TPaveText(0.,1.11,240.,1.15,"br")
pt.SetFillColor(0)
pt.SetLineColor(0)
pt.SetTextAlign(13)
pt.SetTextFont(42)
pt.SetTextSize(0.035)
pt.AddText("CMS Preliminary Simulation, #sqrt{s} = 8 TeV")
pt.Draw("same")

pt2 = TPaveText(5.,1.02,100,1.095,"br")
pt2.SetFillColor(0)
pt2.SetLineColor(0)
pt2.SetTextAlign(13)
pt2.SetTextFont(42)
pt2.AddText("H#rightarrow#gamma#gamma (m_{H} = 125 GeV)")
pt2.AddText("<PU> = 19.9")
pt2.Draw("same")

c1 = TCanvas("vtxProbPt","vtxProbPt")
eff1.SetTitle("Vertex efficiency;p_{T}^{#gamma#gamma} (GeV);Fraction | z_{reco} - z_{true} | < 10 mm")
prob1.SetTitle("Average vertex probability;p_{T}^{#gamma#gamma} (GeV);Fraction | z_{reco} - z_{true} | < 10 mm")
eff1.GetYaxis().SetRangeUser(0.5, 1.1)
eff1.GetXaxis().SetTitleOffset(1.1);
eff1.GetYaxis().SetTitleOffset(1.2);

eff1.Draw("e0")
#prob1.Draw("e3 same")
ge1.Draw("e3 same")
eff1.Draw("e0 same")
#prob1up.Draw("same")
#prob1down.Draw("same")

leg1 = TLegend(0.40,0.24,0.89,0.40)
leg1.SetShadowColor(kWhite), leg1.SetLineColor(kWhite), leg1.SetFillColor(kWhite), leg1.SetTextFont(60)

leg1.AddEntry(eff1,"","pe")
leg1.AddEntry(prob1,"","f")
leg1.Draw("same")
pt.Draw("same")
pt2.Draw("same")


c2 = TCanvas("vtxProbNvtx","vtxProbNvtx")
eff2.SetTitle("True vertex eff.;N_{vtx};Fraction | z_{reco} - z_{true} | < 10 mm")
prob2.SetTitle("Aveage vertex prob.;N_{vtx};Fraction | z_{reco} - z_{true} | < 10 mm")

eff2.GetYaxis().SetRangeUser(0.5, 1.1)
eff2.GetXaxis().SetRangeUser(0, 35)
eff2.GetXaxis().SetTitleOffset(1.1);
eff2.GetYaxis().SetTitleOffset(1.2);
eff2.Draw("e0")
prob2.Draw("e3 same")
#ge2.Draw("e3 same")
eff2.Draw("e0 same")
prob2up.Draw("same")
prob2down.Draw("same")


leg2 = TLegend(0.52,0.62,0.82,0.82)
leg2.SetShadowColor(kWhite), leg2.SetLineColor(kWhite), leg2.SetFillColor(kWhite), leg2.SetTextFont(60)
leg2.AddEntry(eff2,"","pe")
leg2.AddEntry(prob2,"","f")
leg2.Draw("same")
pt.Draw("same")
pt2.Draw("same")


c0 = TCanvas("vtxEffPt","vtxEffPt")
eff0 = eff1.Clone()
eff0.SetTitle("True vertex eff.;p_{T}(#gamma #gamma) (GeV);Fraction of events")
eff0.Draw("e1")
leg0=TLegend(0.332215, 0.357143, 0.60906, 0.423345)
leg0.AddEntry(eff0,"| z_{reco} - z_{true} | < 10 mm","lpf")
leg0.SetShadowColor(kWhite), leg0.SetLineColor(kWhite), leg0.SetFillColor(kWhite), leg0.SetTextFont(60)
leg0.Draw("same")
pt.Draw("same")
pt2.Draw("same")

c3 = TCanvas("vtxEffNvtx","vtxEffNvtx")
eff3 = eff2.Clone()
eff3.SetTitle("True vertex eff.;N_{vtx};Fraction of events")
eff3.Draw("e1")
leg3=TLegend(0.434564,0.684669,0.711409,0.750871)
leg3.SetShadowColor(kWhite), leg3.SetLineColor(kWhite), leg3.SetFillColor(kWhite), leg3.SetTextFont(60)
leg3.AddEntry(eff3,"| z_{reco} - z_{true} | < 10 mm","lpf")
leg3.Draw("same")
pt.Draw("same")
pt2.Draw("same")

for c in c0, c1, c2, c3:
    for fmt in "C", "png", "pdf", "root":
        c.SaveAs( "%s.%s" % (c.GetName(), fmt) )


raw_input("ok?")
