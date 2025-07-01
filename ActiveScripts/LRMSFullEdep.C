//Declarations of headers for classes, functions, etc...
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
using namespace Crane::Analysis;
namespace ca = Crane::Analysis;
namespace cl = Crane::Common;
//using Crane::Calibration;


//Here are some functions that I write to help me in my analysis.

//This function will take in a VolumeID number and then two other numbers that specify which part of the VolumeID you want to interrogate
//so ideally volspec(VolumeID,0,3) will give you the first three numbers of VolumeID which can tell you which CBE part it is. 
//How to use, volspec(12345,1,4) outputs 234 as an integer for your comparing needs
int volspec(int volnum,int a, int b){
	stringstream ss;
	ss << volnum;
        return atoi(ss.str().substr(a, b).c_str());
}


//The two functions below aim to take the SimpleDet mod from 0-35 and translate that into a row and module
int getmod(int sdlayer, int sdmod){
        if(sdlayer % 2 == 0){
                return 5 - (sdmod % 6);
        }else{
                return 5 - (floor(sdmod/6));
        }
}

int getrow(int sdlayer, int sdmod){
        if(sdlayer % 2 == 0){
                return floor(sdmod/6);
        }else{
                return sdmod % 6;
        }
}

//This function aims to take layer, sddet, and sdstrip and translate that into a channel 0-31
int getch(int sdlayer, int sddet, int sdstrp){
	if(sdlayer % 2 == 0){
		if(sddet == 0){return sdstrp + 24;}
		if(sddet == 1){return 23-sdstrp;}
		if(sddet == 2){return sdstrp;}
		if(sddet == 3){return 15-sdstrp;}
		else{return -1;}
	}
	if(sdlayer % 2 != 0){
                if(sddet == 0){return 7 - sdstrp;}
                if(sddet == 1){return 31 - sdstrp;}
                if(sddet == 2){return sdstrp + 8;}
                if(sddet == 3){return sdstrp + 16;}
		else{return -1;}
	}
	else{return -1;}
}




//This is the actual function that you call when you want to run the analysis.
void  FitMe(){

//Filename
char FilenameRoot[400];
//The filename is hard coded here. I'll fix that later.
//sprintf(FilenameRoot,"/data1/nextcloud/cra_data/data/2024/reconstructed/pre-launch/bfsw241210_tof241201_sd241210/runs/gse5/91229125/ethernet241213_*.root"); //UHCRA Run
//sprintf(FilenameRoot,"/data1/nextcloud/cra_data/nobackup/MC/v2.1.2/full/mu-/*.root"); //212 New UHCRA simulated
//sprintf(FilenameRoot,"/data1/nextcloud/cra_data/nobackup/MC/v2.1.0/full/mu-/mu-_gaps_triggerlevel1_FTFP_BERT_HP_1721258929_rec.root"); //210 UHCRA simulated

//sprintf(FilenameRoot,"/home/kelsey/simulations/test/ethernet241213_145/mvtest/ethernet241213_145*.root"); //Personal Computer
//sprintf(FilenameRoot,"/home/kelsey/simulations/test/ethernet241213_145/ethernet241213_1451_rec.root"); //Personal Computer
//sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/simrec/mu-_gaps_triggerlevel1_FTFP_BERT_HP_1721258929_rec.root"); //Old simu data on my computer!
sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/simnew/*.root"); //New sim personal computer

int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TrackerCut = 0.3; //Threshold for an energy deposition to be considered a hit

double xlow = 0.3; //Low range for histogram MeV
double xhigh = 5; //High range for histogram MeV

double coshigh = 0.54; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
double fitlow = 0.5;
double fithigh = 3.5;
const Int_t NBins = 50;
double betacut = 0.8; //Currently we're only doing a beta > 0 cutoff for real data. Beta > 0.8 recommended for sim

//Full tracker histogram range
double mpvmin = 0.66;
double mpvmax = 0.75;

//Number of layers and strips
const int nstrips = 32;
const int nmods = 6;
const int nrows = 6;
const int nlayers = 7;
int dt[4] = {3,4,1,2};

int lyr[nlayers];
int rw[nrows];
int md[nmods];
int strps[nstrips];

for(int l = 0; l < nlayers; l++){lyr[l] = l;}
for(int r = 0; r < nrows; r++){rw[r] = r;}
for(int m = 0; m < nmods; m++){md[m] = m;}
for(int s = 0; s < nstrips; s++){strps[s] = s;}


//Need a histogram and fitting function for every strip
TH1F * h[nlayers][nrows][nmods][nstrips];
TF1 * g1[nlayers][nrows][nmods][nstrips];

auto hcol21 = new TH2F("hcol21","MPV Full Tracker",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);
auto hnentries = new TH2F("hnentries","Full Tracker Strip-Level NHits",nrows*nstrips,0,nrows*nstrips,nlayers*nmods,0,nlayers*nmods);

for(int l = 0;l<nlayers;l++){
        for(int r = 0;r<nrows;r++){
                for(int m = 0; m < nmods; m++){
                        for(int s = 0; s < nstrips; s++){
                                h[l][r][m][s] = new TH1F (TString::Format("h0_l%ir%im%is%i",l,r,m,s), ("Edep l" + to_string(lyr[l]) + "r" + to_string(rw[r]) + "m" + to_string(md[m]) + "s" + to_string(strps[s])).c_str(), NBins, xlow,xhigh);
                                g1[l][r][m][s] = new TF1("g1", "landau", fitlow, fithigh);
                        }
                }
        }
}


//Prepare reconstronstructed event
CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);


//Prepare textile for saving values
std::ofstream myfile;
myfile.open("EdepList.txt");
myfile << TString::Format( "Filename : %s", FilenameRoot )  << endl;
myfile << TString::Format( "Beta Cut : %f", betacut) << endl;
myfile << TString::Format( "Total Entries : %i", static_cast<int>(TreeRec->GetEntries()/MainLoopScaleFactor)) << endl;
myfile << "Layer \t Row \t Mod \t Det \t Strip \t MPV \t Sigma \t NEntries" << endl;
myfile.close();

//For mapping:
map<int, unsigned int> TofIndexVolumeIdMap;

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta cut " << betacut << endl;
cout << "MPV range on 2D Full Tracker will be " << mpvmin << " - " << mpvmax << endl; 

//Using i to loop over every event in the tree
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
	TreeRec->GetEntry(i);

	//Cuts are implemented in this chunk:
	if(Event->GetNTracks() == 1){  //First select the single track event
		bool Umbflag = 0;
		bool CBEtopflag = 0;
		bool CBEbotflag = 0;

		CTrackRec* pt = Event->GetPrimaryTrack();
		uint pt_index = 0;
   	        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

		//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
		if(pt != nullptr && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && /*fabs(Event->GetPrimaryBetaGenerated())*/ fabs(Event->GetPrimaryBeta()) >  betacut){		  	
			//cout << "Event is " << i << endl;

			//-----------EVENT LEVEL CUT APPLIED

			for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                                unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event
                                if(volspec(VolumeId,0,3) == 100){ Umbflag = 1;} // cout << "UMB hit!" <<endl ;
                                if(volspec(VolumeId,0,3) == 110) {CBEtopflag = 1;}// cout << "CBE top hit!" << endl;
                                if(volspec(VolumeId,0,3) == 111) {CBEbotflag = 1;}// cout << "CBE bot hit!" << endl;
                        }

			if(Umbflag && CBEtopflag /*&& CBEbotflag*/ && (pt->GetChi2()/pt->GetNdof()) < 3.2 ){
				//cout << "Event number " << i << " passes the cuts!" << endl;
				for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
					unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig);		
					if(GGeometryObject::IsTrackerVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){

						int layer = GGeometryObject::GetTrackerLayer(VolumeId);
						
						int sdmod = GGeometryObject::GetLayerModule(VolumeId);
						int det = GGeometryObject::GetModuleDetector(VolumeId);
						int sdstrip = GGeometryObject::GetDetectorStrip(VolumeId);
						
						int row = getrow(layer,sdmod);
						int mod = getmod(layer,sdmod);
						int strip = getch(layer, det, sdstrip);

						//cout << "Edep * Cos(theta) " << Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta() << endl;
						h[layer][row][mod][strip]->Fill(  (Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()))   );
						hnentries->Fill(row*32+strip,layer*6+mod);

					} //Closed bracket for Tracker volume and tracker cutoff

				} //Closed bracket for iteration over event with TOF cuts

			} //Closed bracket for if statement for cuts

			//-----------EVENT LEVEL CUTS END


		} //Closed bracket for event level cut

	} //Closed bracket for single track cut

} //Closed bracket for iteration through tree events, move on to the next event i



//Histogram section
//--------------------------------------

myfile.open("EdepList.txt",std::ios::app);

for(int l = 0; l<nlayers;l++){
        for(int r = 0; r<nrows;r++){
                for(int m = 0; m < nmods; m++){
                        for(int j = 0; j < ceil(nstrips/8); j++){
                                cout << "j is " << j << endl;
                                TCanvas * EdepCompare = new TCanvas("EdepCompare", "EdepCompare", 200, 10, 1800, 900);
                                EdepCompare->SetLeftMargin(0.11);
                                EdepCompare->SetRightMargin(0.04);
                                EdepCompare->SetTopMargin(0.04);
                                TLegend* LegEdepCompare = new TLegend(0.5, 0.75, 0.95, 0.95);
                                LegEdepCompare->SetFillColor(0);

                                EdepCompare->Divide(4,2);

                                for(int s = 0 + 8*j; s < 8 + 8*j;s++){
                                        h[l][r][m][s]->SetLineColor(1);
                                        h[l][r][m][s]->GetXaxis()->SetTitle("Energy Deposition of Hit * Cos(#theta) (MeV)");
                                        h[l][r][m][s]->GetYaxis()->SetTitle("Number of Events");
                                        EdepCompare->cd((s % 8)+1);
                                        h[l][r][m][s]->Fit(g1[l][r][m][s],"R");
                                        gPad->SetGridx(1);
                                        gPad->SetGridy(1);
                                        gPad->SetLogy(1);
                                        gStyle->SetTitleW(0.9);
                                        gStyle->SetOptFit();
                                        h[l][r][m][s]->Draw();

					hcol21->Fill(r*32 + s,l*6 + m,g1[l][r][m][s]->GetParameter(1));
					//mpv[r*32+s][l*6+k] = g1[l][r][k][s]->GetParameter(1); //Save the calculated MPV, it will be used for the histogram
                                	myfile << (TString::Format(    "%i \t %i \t %i \t %i \t %f \t %f \t %i \n",l,r,m,s,g1[l][r][m][s]->GetParameter(1),g1[l][r][m][s]->GetParameter(2), static_cast<int>(h[l][r][m][s]->GetEntries())   ));
				
				}
	
                                string title = "FullEdepl" + to_string(lyr[l]) + "r" + to_string(rw[r]) + "m" + to_string(md[m]) + "d" + to_string(dt[j]);
                                char name[400];
                                sprintf(name, "%s.root",title.c_str());
                                EdepCompare->SaveAs(name);
                                sprintf(name, "%s.png",title.c_str());
                                EdepCompare->SaveAs(name);
                                delete EdepCompare;
                        }
                }
        }
}

myfile.close();

TCanvas * c1 = new TCanvas("c1", "c1", 200, 10, 900, 900);
c1->SetLeftMargin(0.1);
c1->SetRightMargin(0.16);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.1);
hcol21->SetBit(TH1::kNoStats);
hcol21->GetXaxis()->SetTitle("row(0-6)*32 + det(0-3)*8 + strp (0-7)");
hcol21->GetYaxis()->SetTitle("layer(0-7)*6 + mod(0-6)");
hcol21->GetZaxis()->SetTitle("Energy Deposition MPV * Cos(#theta) (MeV)");
hcol21->Draw("COLZ");
hcol21->SetMaximum(mpvmax);
hcol21->SetMinimum(mpvmin);

hcol21->SaveAs("hcol21.root");
hcol21->SaveAs("hcol21");

string title = "HistFullTrackerMPV";
char histname[400];
sprintf(histname, "%s.root",title.c_str());
c1->SaveAs(histname);
sprintf(histname, "%s.png",title.c_str());
c1->SaveAs(histname);

//Histogram for NEntries at a strip level
TCanvas * c2 = new TCanvas("c2", "c2", 200, 10, 900, 900);
c2->SetLeftMargin(0.1);
c2->SetRightMargin(0.16);
c2->SetTopMargin(0.1);
c2->SetBottomMargin(0.1);
hnentries->SetBit(TH1::kNoStats);
hnentries->GetXaxis()->SetTitle("row(0-6)*32 + strp (0-31)");
hnentries->GetYaxis()->SetTitle("layer(0-7)*6 + mod(0-6)");
hnentries->GetZaxis()->SetTitle("Number of Hits");
hnentries->Draw("COLZ");

hnentries->SaveAs("hnentries.root");
hnentries->SaveAs("hnentries");

title = "HistFullTrackerNEntries";
sprintf(histname, "%s.root",title.c_str());
c2->SaveAs(histname);
sprintf(histname, "%s.png",title.c_str());
c2->SaveAs(histname);

//--------------------------------------

cout << endl << "I am done" << endl; 

}

