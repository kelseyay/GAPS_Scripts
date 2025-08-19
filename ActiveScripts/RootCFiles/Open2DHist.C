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

//I really need to change this to taking user input, just see DQO

string Filepath = "/home/kelsey/simulations/uhcrapull/LRMSPlots/NewFit241312BetaCut0.2/";
//string Filepath = "~/simulations/uhcrapull/LRMSPlots/NewFitPass4241312BetaCut0.2NoCBEbot/";
//string Filepath = "/home/kelsey/GitScripts/ActiveScripts/EdepTest/";
//string Filepath = "/home/kelsey/GitScripts/ActiveScripts/212RecoBetaCut0.8NoCBEbot3/";
//string Filepath = "/home/kelsey/GitScripts/ActiveScripts/210RecoBetaCut0.8NoCBEbot/";
//string Filepath = "/home/kelsey/simulations/uhcrapull/LRMSPlots/241312Beta0.2NoCBEbot/";
//string Filepath = "/home/kelsey/simulations/uhcrapull/LRMSPlots/Pass4241312BetaCut0.2NoCBEbot/";
char Filename[400];
sprintf(Filename, "%s%s" ,Filepath.c_str(),"hcol21.root" );

//const char Filename[nlen] = (Filepath + "hcol21.root").c_str();

//Fit Range
double mpvmin = 0.5;
double mpvmax = 0.75;

const int nstrips = 32;
const int nmods = 6;
const int nrows = 6;
const int nlayers = 7;
int dt[4] = {3,4,1,2};

auto hcol21 = new TH2F("hcol21","MPV Full Tracker",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);
TFile *f = TFile::Open(Filename);
hcol21 = (TH2F*)f->Get("hcol21");

TCanvas * c1 = new TCanvas("c1", "c1", 200, 10, 900, 900);
c1->SetLeftMargin(0.1);
c1->SetRightMargin(0.15);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.1);
//hcol21->SetBit(TH1::kNoStats);
hcol21->Draw("COLZ");
hcol21->SetStats(false);
hcol21->GetZaxis()->SetTitle("Energy Deposition MPV * Cos(#theta) (MeV)");
hcol21->GetXaxis()->SetTitle("row(0-6)*32 + det(0-3)*8 + strp (0-7)");
hcol21->GetYaxis()->SetTitle("layer(0-7)*6 + mod(0-6)");
hcol21->SetMaximum(mpvmax);
hcol21->SetMinimum(mpvmin);

char Savename[400];
sprintf(Savename, "%s%s" ,Filepath.c_str(),"OpenHistMPV.png" );
c1->SaveAs(Savename);

}
