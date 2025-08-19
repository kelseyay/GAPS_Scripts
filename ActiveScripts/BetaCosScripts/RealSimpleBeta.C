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

/*
string directory = "dirtest";
char SaveDir[600];
sprintf(SaveDir, "mkdir %s", directory.c_str());
*/

int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TrackerCut = 0.4; //Threshold for an energy deposition to be considered a hit
double TofCutLow = 0;

double xlow = 0.3; //Hist range
double xhigh = 6; //Hist range


double coshigh = 0; //0.92; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
double eventbetacut = 0; //Cut that is applied to all events
const Int_t NBins = 50;

//Filename
char FilenameRoot[400];
//The filename is hard coded here. I'll fix that later.
//sprintf(FilenameRoot,"/data1/nextcloud/cra_data/data/2024/reconstructed/pre-launch/bfsw241210_tof241201_sd241210/runs/gse5/91229125/ethernet241213_1*.root");
sprintf(FilenameRoot,"/home/kelsey/simulations/test/ethernet241213_145/mvtest/ethernet241213_145*.root"); //Personal Computer
//sprintf(FilenameRoot,"/home/kelsey/simulations/test/ethernet241213_145/ethernet241213_1451_rec.root"); //Personal Computer
//sprintf(FilenameRoot,"/data1/nextcloud/cra_data/nobackup/MC/v2.1.0/full/mu-/mu-_gaps_triggerlevel1_FTFP_BERT_HP_1721258929_rec.root"); //UHCRA simulated
//sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/simrec/mu-_gaps_triggerlevel1_FTFP_BERT_HP_1721258929_rec.root");  //Simu data on my computer!
//sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/simnew/*.root"); //New Sim on my computer!
//sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/simnew/mu-_gaps_triggerlevel1_FTFP_BERT_1744*_rec.root"); //New Sim Beta down to 0.2?
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
TH1D * HPrBeta = Plotting.DefineTH1D("HPrBeta", 101, -5, 5, "#beta_{rec}", "events", 0.5, 1e7); //For now I've removed Beta Proxy

//Cos
TH1D * HPriCos = Plotting.DefineTH1D("HPriCos", 101, -1.1, 0.1, "Cos_{rec}(theta)", "events", 0.5, 1e7);

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

		//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
		if(pt != nullptr && fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) > eventbetacut){		  	

			//First loop over the event for your needed flags and variables
			for(unsigned int isig = 0; isig < Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
				unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Event->GetVolumeId().at(isig); //Check the VolumeId of the event
				if(volspec(VolumeId,0,3) == 100){ Umbflag = 1;} // cout << "UMB hit!" <<endl ;
				if(volspec(VolumeId,0,3) == 110) {CBEtopflag = 1;}// cout << "CBE top hit!" << endl;
				if(volspec(VolumeId,0,3) == 111) {CBEbotflag = 1;}// cout << "CBE bot hit!" << endl;
			}

			//If the flag checks are successful, proceed to fill the relevant histograms.
			if(Umbflag && CBEtopflag && CBEbotflag && (pt->GetChi2()/pt->GetNdof()) < 3.2 /*true*/){					
					
				HPrBeta->Fill(fabs(Event->GetPrimaryBeta()));
				HPriCos->Fill(Event->GetPrimaryMomentumDirection().CosTheta());

			}
			
		}

	}
}

//Histgrams
//-------------------------------------

TCanvas * CPrRecoBeta = new TCanvas("CPrRecoBeta", "CPrRecoBeta", 200, 10, 900, 900);

CPrRecoBeta->SetLeftMargin(0.11);
CPrRecoBeta->SetRightMargin(0.04);
CPrRecoBeta->SetTopMargin(0.04);
CPrRecoBeta->SetTitle("Reconstructed Beta Histogram");

HPrBeta->SetTitle("Reconstructed Beta Histogram");
HPrBeta->SetLineColor(2);
HPrBeta->Draw("hist");

gPad->SetGridx(1);
gPad->SetGridy(1);
gPad->SetLogy(1);

sprintf(text, "PrRecoBeta.root");
CPrRecoBeta->SaveAs(text);
sprintf(text, "PrRecoBeta.pdf");
CPrRecoBeta->SaveAs(text);
sprintf(text, "PrRecoBeta.png");
CPrRecoBeta->SaveAs(text);


//--------------------------------------

TCanvas * CPriRecoCos = new TCanvas("CPriRecoCos", "CPriRecoCos", 200, 10, 900, 900);

CPriRecoCos->SetLeftMargin(0.11);
CPriRecoCos->SetRightMargin(0.04);
CPriRecoCos->SetTopMargin(0.04);
CPriRecoCos->SetTitle("Reconstructed Cos Histogram");

HPriCos->SetTitle("Reconstructed Cos Histogram");
HPriCos->SetLineColor(2);
HPriCos->Draw("hist");

gPad->SetGridx(1);
gPad->SetGridy(1);
gPad->SetLogy(1);

sprintf(text, "PriRecoCos.root");
CPriRecoCos->SaveAs(text);
sprintf(text, "PriRecoCos.pdf");
CPriRecoCos->SaveAs(text);
sprintf(text, "PriRecoCos.png");
CPriRecoCos->SaveAs(text);

cout << endl << "I am done" << endl; 

}


