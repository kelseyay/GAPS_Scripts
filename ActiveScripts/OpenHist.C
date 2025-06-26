#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TColor.h>
#include <TProfile.h>
#include <TMath.h>
#include <TChain.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TMinuit.h>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

//FIXME: does this work on mac?
#include <sys/stat.h>

//#include "CRawTrk.hh"

#include "CEventMc.hh"
#include "CAnalysisManager.hh"
#include "GAnalysisIdentification.hh"
#include "GBasicTrigger.hh"
#include "GSimulationParameter.hh"
#include "GPreselection.hh"
#include "CraneConstants.hh"
#include "CraneLogging.hh"
#include "GPlottingTools.hh"
#include "CNet.hh"
#include "CBackpropagation.hh"

#include "GGeometry.hh"

#ifdef USE_BOOST_PROGRAM_OPTIONS
#include "GOptionParser.hh"
#include "GFileIO.hh"
#endif

//Namespaces help you call functions without having to explicitly write out std, Crane, etc... each time.
using namespace std;

void  HistOpen(){

TFile *f = TFile::Open("212hstrip.root");
TFile *mcmf = TFile::Open("p4hstrip.root");

double xlow = 0.3; //Low range for histogram MeV
double xhigh = 4; //High range for histogram MeV
double fitlow = 0.5;
double fithigh = 3.5;
const Int_t NBins = 50;

//Just need one histogram and one fit for one strip

TH1F * hpstrip = new TH1F ("hpstrip","hpstrip", NBins, xlow,xhigh);
hpstrip = (TH1F*)mcmf->Get("h0_l0r5m3s17");
TF1 * gpstrip = new TF1("gp1", "landau", fitlow, fithigh);

TH1F * h212strip = new TH1F ("h212strip","h212strip", NBins, xlow,xhigh);	
h212strip = (TH1F*)f->Get("h0_l0r5m3s17");
TF1 * gstrip = new TF1("g1", "landau", fitlow, fithigh);

TCanvas * c1 = new TCanvas("c1", "c1", 200, 10, 900, 900);
c1->SetLeftMargin(0.12);
c1->SetRightMargin(0.1);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.1);

hpstrip->SetLineColor(1);
gpstrip->SetLineColor(1);
gstrip->SetLineColor(2);
h212strip->SetLineColor(2);

hpstrip->GetXaxis()->SetRangeUser(xlow, xhigh);
hpstrip->GetXaxis()->SetTitle("Energy Deposition of Hit #times  Cos(#theta) [MeV]");
hpstrip->GetYaxis()->SetTitle("Number of Hits on Track");


gStyle->SetOptStat(0);
gPad->SetLogy(1);
hpstrip->Draw();
h212strip->Draw("same");

h212strip->Fit(gstrip,"R");
hpstrip->Fit(gpstrip,"R");

stringstream mcms;
mcms << "MPV:      " << std::fixed << std::setprecision(3) << gpstrip->GetParameter(1) << " #pm   " << std::setprecision(3) << gpstrip->GetParError(1)  << " MeV ";

stringstream s212;
s212 << "MPV:      " << std::fixed << std::setprecision(3) << gstrip->GetParameter(1) << " #pm   " << std::setprecision(3) << gstrip->GetParError(1)  << " MeV ";

hpstrip->SetTitle("");

auto legend = new TLegend(0.4,0.6,0.9,0.9); //TLegend (Double_t x1, Double_t y1, Double_t x2, Double_t y2, const char *header="", Option_t *option="brNDC")
legend->SetTextSize(0.025);
legend->AddEntry(hpstrip,"McMurdo data: Run 9125","l");
legend->AddEntry((TObject*)0, (mcms.str()).c_str()  , "");
legend->AddEntry(h212strip,"Muon MC","l");
legend->AddEntry((TObject*)0, (s212.str()).c_str()  , "");

//legend->AddEntry((TObject*)0, TString::Format("Sigma: \t \t \t \t \t %f #pm %f ",gpstrip->GetParameter(2),gpstrip->GetParError(2)   ), "");
//legend->AddEntry((TObject*)0, TString::Format("#chi^{2}/Ndof: \t \t %f/%d ",gpstrip->GetChisquare(), gpstrip->GetNDF()   ), "");
//legend->AddEntry((TObject*)0, TString::Format("MPV: \t \t \t \t \t \t \t %s #pm %f ", ss.str(),gstrip->GetParError(1)   ), "");
//legend->AddEntry((TObject*)0, TString::Format("Sigma: \t \t \t \t \t %f #pm %f ",gstrip->GetParameter(2),gstrip->GetParError(2)   ), "");
//legend->AddEntry((TObject*)0, TString::Format("#chi^{2}/Ndof: \t \t %f/%d ",gstrip->GetChisquare(), gstrip->GetNDF()   ), "");

legend->Draw();


TCanvas * c2 = new TCanvas("c2", "c2", 200, 10, 900, 900);
c2->SetLeftMargin(0.12);
c2->SetRightMargin(0.1);
c2->SetTopMargin(0.1);
c2->SetBottomMargin(0.1);

gStyle->SetOptStat(0);
gPad->SetLogy(1);
hpstrip->Draw();

auto legend2 = new TLegend(0.4,0.75,0.9,0.9); //TLegend (Double_t x1, Double_t y1, Double_t x2, Double_t y2, const char *header="", Option_t *option="brNDC")

legend2->SetTextSize(0.025);
//legend->SetLegendTextSize(0.5);
legend2->AddEntry(hpstrip,"McMurdo data: Run 9125","l");
legend2->AddEntry((TObject*)0, (mcms.str()).c_str()  , "");
legend2->Draw();

c1->SaveAs("MCM_MC.png");
c2->SaveAs("MCM.png");
}
