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
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TMinuit.h>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include "TH1.h"
#include "TCanvas.h"

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

using namespace std;
using namespace Crane::Analysis;
namespace ca = Crane::Analysis;
namespace cl = Crane::Common;
//using Crane::Calibration;

//This function will take in a VolumeID number and then two other numbers that specify which part of the VolumeID you want to interrogate
//so ideally volspec(VolumeID,0,3) will give you the first three numbers of VolumeID which can tell you which CBE part it is. 
//How to use, volspec(12345,1,4) outputs 234 as an integer for your comparing needs
int volspec(int volnum,int a, int b){
	stringstream ss;
	ss << volnum;
        return atoi(ss.str().substr(a, b).c_str());
	}


void  BetaTest(){

int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TrackerCut = 0.4; //Threshold for an energy deposition to be considered a hit
double TofCutLow = 0;

double xlow = 0.3; //Hist range
double xhigh = 6; //Hist range

double coshigh = 0; //0.92; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
double eventbetacut = 0.2; //Cut that is applied to all events
double betacut = 0.9; //Separation of quandrants in the beta plot
const Int_t NBins = 50;

//Create Counters for each Beta quadrant
int BetaGenCts = 0;
int TrueBetaOvrCutRecBOvrCutCts = 0; //True Beta > betacut, and Reco Beta > betacut
int TrueBetaSubCutRecBSubCutCts = 0; //True Beta < betacut, and Reco beta > betacut
int TrueBetaOvrCutRecBSubCutCts = 0; //True Beta > betacut, and Reco beta < betacut
int TrueBetaSubCutRecBOvrCutCts = 0; //True Beta < betacut, and Reco beta < betacut

//double TofCut = 0;//300 // Currently no TofCut 


//Filename
char FilenameRoot[400];
//The filename is hard coded here. I'll fix that later.
//sprintf(FilenameRoot,"/data1/nextcloud/cra_data/data/2024/reconstructed/pre-launch/bfsw241210_tof241201_sd241210/runs/gse5/91229125/ethernet241213_1*.root");
//sprintf(FilenameRoot,"/home/kelsey/simulations/test/ethernet241213_145/mvtest/ethernet241213_145*.root"); //Personal Computer
//sprintf(FilenameRoot,"/home/kelsey/simulations/test/ethernet241213_145/ethernet241213_1451_rec.root"); //Personal Computer
//sprintf(FilenameRoot,"/data1/nextcloud/cra_data/nobackup/MC/v2.1.0/full/mu-/mu-_gaps_triggerlevel1_FTFP_BERT_HP_1721258929_rec.root"); //UHCRA simulated
//sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/simrec/mu-_gaps_triggerlevel1_FTFP_BERT_HP_1721258929_rec.root");  //Simu data on my computer!
//sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/simnew/*.root"); //New Sim on my computer!
sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/simnew/mu-_gaps_triggerlevel1_FTFP_BERT_1744*_rec.root"); //New Sim Beta down to 0.2?
//sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/simnew/mu-_gaps_triggerlevel1_FTFP_BERT_1746072*rec.root"); //New Sim Beta 0.9 and up?
//sprintf(FilenameRoot,"/data1/nextcloud/cra_data/nobackup/MC/v2.1.2/full/mu-/mu-_gaps_triggerlevel1_FTFP_BERT_*_rec.root"); //New UHCRA simulated!

//Prepare reconstronstructed event
CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

CAnalysisManager AnalysisManagerRec;
AnalysisManagerRec.SetEvent(Event);

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

//All of the plots are declared here
//("Title",Number of bins,xmin,xmax,"xlabel","ylabel",ymin,ymax)

//Beta Section
TH1D * HBetaProxy = Plotting.DefineTH1D("HBetaProxy", 50, 0, 2, "1/#sqrt{E_{dep,mean} cos(zenith)}", "events", 0.5, 1e7);
TH1D * HBeta = Plotting.DefineTH1D("HBeta", 101, -5, 5, "#beta_{rec}", "events", 0.5, 1e7);
TH1D * HPriBeta = Plotting.DefineTH1D("HPriBeta", 101, -0.1, 1.2, "#beta_{rec}", "events", 0.5, 1e7);
TH1D * HPriGenBeta = Plotting.DefineTH1D("HPriGenBeta", 101, -0.1, 1.2, "#beta", "events", 0.5, 1e7);

//TH2D * HBetaBetaProxy = Plotting.DefineTH2D("HBetaBetaProxy", 30, 0, 1, 60, 0, 3, "|#beta_{rec}|", "1/#sqrt{E_{dep,mean} cos(zenith)}", "events", 0.5, 1e5);
TH2D* HBetaBetaProxy = Plotting.DefineTH2D("HBetaBetaProxy", 30, 0, 1, 60, 0, 3, "|#beta_{rec}|", "1/#sqrt{E_{dep,mean} cos(zenith)}", "events", 0.5, 1e5);

//Zenith Cos
TH2D* HTrackerCosZenithEnergyDepositionAll = Plotting.DefineTH2D("HTrackerCosZenithEnergyDepositionAll", 20, -1, 0, 100, 0, 10, "cos(zenith)", "energy deposition [MeV]", "events [normalized]", 1e-5, 1);

//Cos
TH1D * HPriCos = Plotting.DefineTH1D("HPriCos", 101, -1.1, 0.1, "#cos_{rec}(theta)", "events", 0.5, 1e7);
TH1D * HPriGenCos = Plotting.DefineTH1D("HPriGenCos", 101, -1.1, 0.1, "#cos_{truth}(theta)", "events", 0.5, 1e7);

//2D Histos
auto hbetacomp = new TH2F("hbetacomp","Beta Generated vs Beta Reconstructed",50,eventbetacut,1,50,0,1.6);
auto hcoscomp = new TH2F("hcoscomp","|Cos(theta) Generated| vs |Cos(theta) Reconstructed|",50,0,1,50,0,1);


//Now we can go over the loop

int PercentageStep = 5;
TreeRec->GetEntry(0);

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
//for(unsigned int i = 5360; i < 5380; i+=MainLoopScaleFactor) //Let's try a smaller maximum number of events
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
//for(unsigned int i = TreeRec->GetEntries()-1000; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor)
	TreeRec->GetEntry(i);

	//Cuts are implemented in this chunk:
	if(Event->GetNTracks() == 1){  //First select the single track event

		//cout << "Single Track Event " << endl;
		bool Umbflag = 0;
		bool CBEtopflag = 0;
		bool CBEbotflag = 0;

		CTrackRec* pt = Event->GetPrimaryTrack();
		uint pt_index = 0;
   	        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

		//if(pt != nullptr) //You need some kind of cut to make sure primary track exists!
		//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
		if(-fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBetaGenerated()) > eventbetacut){		  	

			HPriGenBeta->Fill(Event->GetPrimaryBetaGenerated());
			HPriBeta->Fill(fabs(Event->GetPrimaryBeta()));
			HPriGenCos->Fill(Event->GetPrimaryMomentumDirectionGenerated().CosTheta());
			HPriCos->Fill(Event->GetPrimaryMomentumDirection().CosTheta());
			hbetacomp->Fill(Event->GetPrimaryBetaGenerated(),fabs(Event->GetPrimaryBeta()));
			hcoscomp->Fill(fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()),fabs(Event->GetPrimaryMomentumDirection().CosTheta()));
			
			BetaGenCts++;
			if(fabs(Event->GetPrimaryBetaGenerated()) > betacut && fabs(Event->GetPrimaryBeta()) > betacut){ TrueBetaOvrCutRecBOvrCutCts++; }
			if(fabs(Event->GetPrimaryBetaGenerated()) > betacut && fabs(Event->GetPrimaryBeta()) < betacut){ TrueBetaOvrCutRecBSubCutCts++; }
			if(fabs(Event->GetPrimaryBetaGenerated()) < betacut && fabs(Event->GetPrimaryBeta()) > betacut){ TrueBetaSubCutRecBOvrCutCts++; }
			if(fabs(Event->GetPrimaryBetaGenerated()) < betacut && fabs(Event->GetPrimaryBeta()) < betacut){ TrueBetaSubCutRecBSubCutCts++; }

			//First loop over the event for your needed flags and variables
			for(unsigned int isig = 0; isig < Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
				unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Event->GetVolumeId().at(isig); //Check the VolumeId of the event
				if(volspec(VolumeId,0,3) == 100){ Umbflag = 1;} // cout << "UMB hit!" <<endl ;
				if(volspec(VolumeId,0,3) == 110) {CBEtopflag = 1;}// cout << "CBE top hit!" << endl;
				if(volspec(VolumeId,0,3) == 111) {CBEbotflag = 1;}// cout << "CBE bot hit!" << endl;
			}


			//If the loop was successful, proceed to fill the relevant histograms.
			if(Umbflag && CBEtopflag /*&& CBEbotflag*/ && (pt->GetChi2()/pt->GetNdof()) < 3.2 /*true*/){					
				
				/*	
				HPriGenBeta->Fill(Event->GetPrimaryBetaGenerated());
				HPriBeta->Fill(fabs(Event->GetPrimaryBeta()));
				HPriGenCos->Fill(Event->GetPrimaryMomentumDirectionGenerated().CosTheta());
				HPriCos->Fill(Event->GetPrimaryMomentumDirection().CosTheta());
				hbetacomp->Fill(Event->GetPrimaryBetaGenerated(),fabs(Event->GetPrimaryBeta()));
				hcoscomp->Fill(fabs(Event->GetPrimaryMomentumDirectionGenerated().CosTheta()),fabs(Event->GetPrimaryMomentumDirection().CosTheta()));
				
				BetaGenCts++;
				if(fabs(Event->GetPrimaryBetaGenerated()) > betacut && fabs(Event->GetPrimaryBeta()) > betacut){ TrueBetaOvrCutRecBOvrCutCts++; }
				if(fabs(Event->GetPrimaryBetaGenerated()) > betacut && fabs(Event->GetPrimaryBeta()) < betacut){ TrueBetaOvrCutRecBSubCutCts++; }
				if(fabs(Event->GetPrimaryBetaGenerated()) < betacut && fabs(Event->GetPrimaryBeta()) > betacut){ TrueBetaSubCutRecBOvrCutCts++; }
				if(fabs(Event->GetPrimaryBetaGenerated()) < betacut && fabs(Event->GetPrimaryBeta()) < betacut){ TrueBetaSubCutRecBSubCutCts++; }
				*/
			}
			



		}
	}
}

cout << "Beta Cut " << betacut << endl;
cout << "Total Generated Beta " << BetaGenCts << endl;
cout << "True β > " << betacut << ", Reco β > " << betacut << " NEvents = " << TrueBetaOvrCutRecBOvrCutCts << endl;
cout << "True β > " << betacut << ", Reco β < " << betacut <<  " NEvents = " << TrueBetaOvrCutRecBSubCutCts << endl;
cout << endl;
cout << "True β < " << betacut << ", Reco β < " << betacut << " NEvents = " << TrueBetaSubCutRecBSubCutCts << endl;
cout << "True β < " << betacut << ", Reco β > " << betacut << " NEvents = " << TrueBetaSubCutRecBOvrCutCts << endl;
cout << endl;
cout << "You would lose " << 1.0*TrueBetaOvrCutRecBSubCutCts / (1.0*TrueBetaOvrCutRecBSubCutCts + 1.0*TrueBetaOvrCutRecBOvrCutCts)  << "  of events that actually have β > 0.9 with a > 0.9 Reconstructed β cut" << endl;

//Histogram section
//Generated Beta and Theta
//2D Histo 

TCanvas * c5 = new TCanvas("c5", "c5", 200, 10, 900, 900);
c5->SetLeftMargin(0.1);
c5->SetRightMargin(0.16);
c5->SetTopMargin(0.1);
c5->SetBottomMargin(0.1);
hbetacomp->SetBit(TH1::kNoStats);
hbetacomp->GetXaxis()->SetTitle("Generated Beta");
hbetacomp->GetYaxis()->SetTitle("Reconstructed Beta");
hbetacomp->GetZaxis()->SetTitle("Number of Entries");
hbetacomp->Draw("COLZ");
gPad->SetLogz();
//hbetacomp->SetMaximum(1.3*mpvave);
//hbetacomp->SetMinimum(0.7*mpvave);

string bettitle = "BetaGen_vs_Reco";
char histname[400];
sprintf(histname, "%s.root",bettitle.c_str());
c5->SaveAs(histname);
sprintf(histname, "%s.png",bettitle.c_str());
c5->SaveAs(histname);
sprintf(histname, "%s.pdf",bettitle.c_str());
c5->SaveAs(histname);


TCanvas * c6 = new TCanvas("c6", "c6", 200, 10, 900, 900);
c6->SetLeftMargin(0.1);
c6->SetRightMargin(0.16);
c6->SetTopMargin(0.1);
c6->SetBottomMargin(0.1);
hcoscomp->SetBit(TH1::kNoStats);
hcoscomp->GetXaxis()->SetTitle("Generated Cos(theta)");
hcoscomp->GetYaxis()->SetTitle("Reconstructed Cos(theta)");
hcoscomp->GetZaxis()->SetTitle("Number of Entries");
hcoscomp->Draw("COLZ");
gPad->SetLogz();
//hbetacomp->SetMaximum(1.3*mpvave);
//hbetacomp->SetMinimum(0.7*mpvave);

string costitle = "CosGen_vs_Reco";
sprintf(histname, "%s.root",costitle.c_str());
c6->SaveAs(histname);
sprintf(histname, "%s.png",costitle.c_str());
c6->SaveAs(histname);
sprintf(histname, "%s.pdf",costitle.c_str());
c6->SaveAs(histname);



TCanvas * CPriGenBeta = new TCanvas("CPriGenBeta", "CPriGenBeta", 200, 10, 900, 900);

CPriGenBeta->SetLeftMargin(0.11);
CPriGenBeta->SetRightMargin(0.04);
CPriGenBeta->SetTopMargin(0.04);
//CPriGenBeta->SetTitle("Generated vs Reconstructed Beta Histogram");

HPriGenBeta->SetLineColor(2);
HPriGenBeta->SetTitle((to_string(betacut) + " Beta Cut Generated vs Reconstructed Beta Histogram").c_str());
HPriBeta->SetLineColor(1);
HPriGenBeta->Draw("hist");
HPriBeta->Draw("same");

gPad->SetGridx(1);
gPad->SetGridy(1);
gPad->SetLogy(1);

auto legend = new TLegend(0.2,0.7,0.5,0.9); //(x1, y1, x2, y2) it will draw a rectangle with those coordinates
//legend->SetHeader("Legend"); // option "C" allows to center the header
legend->AddEntry(HPriGenBeta,"Generated Beta");
legend->AddEntry(HPriBeta,"Reconstructed Beta");
legend->Draw();


sprintf(text, "PriGenBeta.root");
CPriGenBeta->SaveAs(text);
sprintf(text, "PriGenBeta.pdf");
CPriGenBeta->SaveAs(text);
sprintf(text, "PriGenBeta.png");
CPriGenBeta->SaveAs(text);


//--------------------------------------

TCanvas * CPriGenCos = new TCanvas("CPriGenCos", "CPriGenCos", 200, 10, 900, 900);

CPriGenCos->SetLeftMargin(0.11);
CPriGenCos->SetRightMargin(0.04);
CPriGenCos->SetTopMargin(0.04);
//CPriGenCos->SetTitle("Generated vs Reconstructed Cos Histogram");

HPriGenCos->SetLineColor(2);
HPriGenCos->SetTitle( ((to_string(betacut) +  " Beta Cut Generated vs Reconstructed Cos Histogram" )).c_str());
HPriCos->SetLineColor(1);
HPriGenCos->Draw("hist");
HPriCos->Draw("same");

gPad->SetGridx(1);
gPad->SetGridy(1);
gPad->SetLogy(1);

auto legend2 = new TLegend(0.2,0.7,0.5,0.9);
//legend->SetHeader("Legend"); // option "C" allows to center the header
legend2->AddEntry(HPriGenCos,"Generated Cos");
legend2->AddEntry(HPriCos,"Reconstructed Cos");
legend2->Draw();


sprintf(text, "PriGenCos.root");
CPriGenCos->SaveAs(text);
sprintf(text, "PriGenCos.pdf");
CPriGenCos->SaveAs(text);
sprintf(text, "PriGenCos.png");
CPriGenCos->SaveAs(text);

cout << endl << "I am done" << endl; 

}


