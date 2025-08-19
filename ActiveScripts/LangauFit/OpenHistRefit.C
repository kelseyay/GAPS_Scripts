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
#include "/home/kelsey/Downloads/langaufun.C"

#include "GGeometry.hh"

#ifdef USE_BOOST_PROGRAM_OPTIONS
#include "GOptionParser.hh"
#include "GFileIO.hh"
#endif

//Namespaces help you call functions without having to explicitly write out std, Crane, etc... each time.
using namespace std;

void  HistOpen(){

TFile *mcmf = TFile::Open("p4hstrip.root");

double xlow = 0.3; //Low range for histogram MeV
double xhigh = 4; //High range for histogram MeV
double fitlow = 0.5;
double fithigh = 3.5;
const Int_t NBins = 50;

//Just need one histogram and one fit for one strip

TH1F * hpstrip = new TH1F ("p4hstrip","p4hstrip", NBins, xlow,xhigh);
hpstrip = (TH1F*)mcmf->Get("h0_l0r5m3s17");

double min = 0, max = 10; //Change this to your preferences
int numParameter = 4;

TF1*  f_landau_mod = new TF1("f_landau_gauss",langaufun,min,max,numParameter);

f_landau_mod->SetParameter(0, 1);
f_landau_mod->SetParameter(1, 1);
f_landau_mod->SetParameter(2, 1);
f_landau_mod->SetParameter(3, 1);

hpstrip->Fit(f_landau_mod,"R");

hpstrip->Draw();


TCanvas * c1 = new TCanvas("c1", "c1", 200, 10, 900, 900);
c1->SetLeftMargin(0.12);
c1->SetRightMargin(0.1);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.1);

hpstrip->SetLineColor(1);

hpstrip->GetXaxis()->SetRangeUser(xlow, xhigh);
hpstrip->GetXaxis()->SetTitle("Energy Deposition of Hit #times  Cos(#theta) [MeV]");
hpstrip->GetYaxis()->SetTitle("Number of Hits on Track");


gStyle->SetOptStat(0);
gPad->SetLogy(1);
hpstrip->Draw();
hpstrip->SetTitle("");
hpstrip->SaveAs("p4hLGFit.root");

}
