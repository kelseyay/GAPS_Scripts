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


void  OccuTest(){

int MainLoopScaleFactor = 10; //Set this number to scale the step size. Larger means runs faster and fewer events
double TrackerCut = 0.4; //Threshold for an energy deposition to be considered a hit
double TofCutLow = 0;

double xlow = 0.3; //Hist range
double xhigh = 6; //Hist range

double coshigh = 0; //0.92; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
double eventbetacut = 0.2; //Cut that is applied to all events
double betacut = 0.9; //Separation of quandrants in the beta plot
const Int_t NBins = 50;

//Filename
char FilenameRoot[400];
//The filename is hard coded here. I'll fix that later.
//sprintf(FilenameRoot,"/data1/nextcloud/cra_data/data/2024/reconstructed/pre-launch/bfsw241210_tof241201_sd241210/runs/gse5/91229125/ethernet241213_1*.root");
//sprintf(FilenameRoot,"/home/kelsey/simulations/test/ethernet241213_145/mvtest/ethernet241213_145*.root"); //Personal Computer
//sprintf(FilenameRoot,"/home/kelsey/simulations/test/ethernet241213_145/ethernet241213_1451_rec.root"); //Personal Computer
//sprintf(FilenameRoot,"/data1/nextcloud/cra_data/nobackup/MC/v2.1.0/full/mu-/mu-_gaps_triggerlevel1_FTFP_BERT_HP_1721258929_rec.root"); //UHCRA simulated
sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/simrec/mu-_gaps_triggerlevel1_FTFP_BERT_HP_1721258929_rec.root");  //Simu data on my computer!
//sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/simnew/*.root"); //New Sim on my computer!
//sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/simnew/mu-_gaps_triggerlevel1_FTFP_BERT_1744*_rec.root"); //New Sim Beta down to 0.2?
//sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/simnew/mu-_gaps_triggerlevel1_FTFP_BERT_1746072*rec.root"); //New Sim Beta 0.9 and up?
//sprintf(FilenameRoot,"/data1/nextcloud/cra_data/nobackup/MC/v2.1.2/full/mu-/mu-_gaps_triggerlevel1_FTFP_BERT_*_rec.root"); //New UHCRA simulated!

//Prepare reconstronstructed event
CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

//CAnalysisManager AnalysisManagerRec;
//AnalysisManagerRec.SetEvent(Event);
//These two lines are currently causing a segfault at the very end of the code, I'm confused. I think these are a relic of Philip's BetaProx whatever
//I conclude these are no longer useful so just move on lol. 

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

//All of the plots are declared here
//("Title",Number of bins,xmin,xmax,"xlabel","ylabel",ymin,ymax)

//2D Histos
TH2D* HTofOccu = Plotting.DefineTH2D("HTofOccu", 30, -2000, 2000, 30, -2000, 2000, "rec. hit position x [mm]", "rec. hit position y [mm]", "events", 0.5, 1e3);
TH2D* HTofCBEtopOccu = Plotting.DefineTH2D("HTofCBEtopOccu", 30, -2000, 2000, 30, -2000, 2000, "rec. hit position x [mm]", "rec. hit position y [mm]", "events", 0.5, 1e3);

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
		if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBetaGenerated()) > eventbetacut){		  	
	
			//cout << endl << "Event is " << i << endl;
			for(unsigned int k = 0; k < Event->GetTrack(0)->GetEnergyDeposition().size(); k++){	
				unsigned int VolumeId = Event->GetTrack(0)->GetVolumeId(k);
				if(GGeometryObject::IsTofVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(k) > TofCutLow){
					if(volspec(VolumeId,0,3) == 100){
                                                HTofOccu->Fill(Event->GetTrack(0)->GetPosition(k).X()+Event->GetTrack(0)->GetPositionResidual(k).X(), Event->GetTrack(0)->GetPosition(k).Y()+Event->GetTrack(0)->GetPositionResidual(k).Y());
					}
					if(volspec(VolumeId,0,3) == 110){
						HTofCBEtopOccu->Fill(Event->GetTrack(0)->GetPosition(k).X()+Event->GetTrack(0)->GetPositionResidual(k).X(), Event->GetTrack(0)->GetPosition(k).Y()+Event->GetTrack(0)->GetPositionResidual(k).Y());
					}
				}
			}

		/*
			//First loop over the event for your needed flags and variables
			for(unsigned int isig = 0; isig < Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
				unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Event->GetVolumeId().at(isig); //Check the VolumeId of the event
				if(volspec(VolumeId,0,3) == 100){ Umbflag = 1;} // cout << "UMB hit!" <<endl ;
				if(volspec(VolumeId,0,3) == 110) {CBEtopflag = 1;}// cout << "CBE top hit!" << endl;
				if(volspec(VolumeId,0,3) == 111) {CBEbotflag = 1;}// cout << "CBE bot hit!" << endl;
			}


			//If the loop was successful, proceed to fill the relevant histograms.
			if(Umbflag && CBEtopflag && CBEbotflag && (pt->GetChi2()/pt->GetNdof()) < 3.2){					
				
			}
			
			//cout << endl;
			//cout << "Low Angle " << Event->GetTrack(0)->GetEnergyDeposition().size() << endl;
		*/
		}

	}
}


//Histogram section


TCanvas * c1 = new TCanvas("c1", "c1", 200, 10, 900, 900);
c1->SetLeftMargin(0.1);
c1->SetRightMargin(0.16);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.1);
HTofOccu->SetBit(TH1::kNoStats);
HTofOccu->GetXaxis()->SetTitle("Generated Cos(theta)");
HTofOccu->GetYaxis()->SetTitle("Reconstructed Cos(theta)");
HTofOccu->GetZaxis()->SetTitle("Number of Entries");
HTofOccu->Draw("COLZ");
gPad->SetLogz();

char histname[400];
string TofUmbTitle = "TofUmbOccu";
sprintf(histname, "%s.root",TofUmbTitle.c_str());
c1->SaveAs(histname);
sprintf(histname, "%s.png",TofUmbTitle.c_str());
c1->SaveAs(histname);
sprintf(histname, "%s.pdf",TofUmbTitle.c_str());
c1->SaveAs(histname);



TCanvas * c2 = new TCanvas("c2", "c2", 200, 10, 900, 900);
c2->SetLeftMargin(0.1);
c2->SetRightMargin(0.16);
c2->SetTopMargin(0.1);
c2->SetBottomMargin(0.1);
HTofCBEtopOccu->SetBit(TH1::kNoStats);
HTofCBEtopOccu->GetXaxis()->SetTitle("Generated Cos(theta)");
HTofCBEtopOccu->GetYaxis()->SetTitle("Reconstructed Cos(theta)");
HTofCBEtopOccu->GetZaxis()->SetTitle("Number of Entries");
HTofCBEtopOccu->Draw("COLZ");
gPad->SetLogz();

string TofCBEtopTitle = "TofCBEtopOccu";
sprintf(histname, "%s.root",TofCBEtopTitle.c_str());
c2->SaveAs(histname);
sprintf(histname, "%s.png",TofCBEtopTitle.c_str());
c2->SaveAs(histname);
sprintf(histname, "%s.pdf",TofCBEtopTitle.c_str());
c2->SaveAs(histname);



cout << endl << "I am done" << endl; 

}


