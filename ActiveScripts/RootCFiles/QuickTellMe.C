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


void  TellMe(){

int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
//double TrackerCut = 0.4; //Threshold for an energy deposition to be considered a hit
//double TofCut = 0;//300 // Currently no TofCut 


//Filename
char FilenameRoot[400];
//The filename is hard coded here. I'll fix that later.
//sprintf(FilenameRoot,"/data1/nextcloud/cra_data/data/2024/reconstructed/pre-launch/bfsw241210_tof241201_sd241210/runs/91229125/ethernet241213_145*.root");
//sprintf(FilenameRoot,"/home/kelsey/simulations/test/ethernet241213_145/mvtest/ethernet241213_145*.root"); //Personal Computer
//sprintf(FilenameRoot,"/home/kelsey/simulations/test/ethernet241213_145/ethernet241213_1451_rec.root"); //Personal Computer
//sprintf(FilenameRoot,"mu_gaps_FTFP_BERT_HP_Trigger1_09_1_beta_Flat_FTFP_BERT_HP_1708424520.root"); //Personal Computer
//sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/simrec/mu-_gaps_triggerlevel1_FTFP_BERT_HP_1721258929_rec.root");  //Ols Simu data on my computer!
sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/simnew/mu-_gaps_triggerlevel1_FTFP_BERT_1744342800_rec.root");
//TFile *FilenameRoot = new TFile("./mu_gaps_FTFP_BERT_HP_Trigger1_09_1_beta_Flat_FTFP_BERT_HP_1708424520.root","READ");

//Prepare reconstronstructed event
CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
//TTree * TreeRec = (TTree*)FilenameRoot->Get("TreeRec");
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

//TGeoManager* geoMan = LoadGeometryRoot("/home/kelsey/simulations/test/ethernet241213_145/ethernet241213_1451_rec.root");
//TGeoManager* geoMan = LoadGeometryRoot("/home/kelsey/simulations/simdat/simrec/mu-_gaps_triggerlevel1_FTFP_BERT_HP_1721258929_rec.root");
TGeoManager* geoMan = LoadGeometryRoot("/home/kelsey/simulations/simdat/simnew/mu-_gaps_triggerlevel1_FTFP_BERT_1744342800_rec.root");
GGeometryMapPtr mcGeo = UnpackGeometry(geoMan, true, true);


//  CEventRec* reco_event_simu_ = new CEventRec(); //New CEventRec (SimpleDet class reference) 
//  TTree * _recoTree_Simu_ = (TTree*)file_simu_mu->Get("TreeRec"); //New TTree in muon simulated file. TreeRec
//  _recoTree_Simu_->SetBranchAddress("Rec",&reco_event_simu_);

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;
int singevents = 0;
int pritracks = 0;
int neventspass = 0;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

//All of the plots are declared here
//("Title",Number of bins,xmin,xmax,"xlabel","ylabel",ymin,ymax)
TH1D * HNTracks = Plotting.DefineTH1D("HNTracks", 51, -0.5, 50.5, "number of tracks", "events", 0.5, 1e8);

//Since number of hits is integer, then just use 15 bins
TH1D * HCosLowAng = Plotting.DefineTH1D("HCosLowAng", 15, 0, 15, "Num Hits on Track", "events", 0.5, 1e6);
TH1D * HEdepLowAng = Plotting.DefineTH1D("HEdepLowAng", 80, 0, 20, "Energy Deposition of Hit (MeV)", "Number of Events", 0.5, 1e6);

TH1D * HChi2 = Plotting.DefineTH1D("HChi2", 51, -1, 50, "Chi2", "events", 0.5, 1e4);
TH1D * HNdof = Plotting.DefineTH1D("HNdof", 50, 0, 10, "Ndof", "events", 0.5, 1e4);
TH1D * HChi2_Ndof = Plotting.DefineTH1D("HChi2_Ndof", 50, -5, 10, "Ndof", "events", 0.5, 1e4);

//Now we can go over the loop

int PercentageStep = 5;
TreeRec->GetEntry(0);

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
for(unsigned int i = 0; i < 500; i+=MainLoopScaleFactor){ //Let's try a smaller maximum number of events
//for(int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
        TreeRec->GetEntry(i);
	if( TreeRec->GetEntries() % (i+1) == 0){cout << "Time at Event " << i << " = " << Event->GetEventTime() << endl;}

	CTrackRec* pt = Event->GetPrimaryTrack();
	uint pt_index = 0;
	for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

	bool Umbflag = 0;
	bool CBEtopflag = 0;
	bool CBEbotflag = 0;

	if(Event->GetNTracks() == 1){ 
		cout << "Event is " << i << endl; 
		if(pt != nullptr){
			for(unsigned int isig = 0; isig < Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){ //Loop over
				unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Event->GetVolumeId().at(isig); //Check the VolumeId of the event
                                //unsigned int VolumeId = Event->GetVolumeId().at(isig);
                                //cout << "First three Volume ID Check" << volspec(VolumeId,0,3)  << endl;
                                //cout << "VolumeId is " << VolumeId <<endl;
                                if(GGeometryObject::IsUmbrellaVolume(VolumeId)){ Umbflag = 1; cout << "UMB hit!" <<endl; }
                                //if(volspec(VolumeId,0,3) == 110) {CBEtopflag = 1; cout << "CBE top hit!" << endl; }
                                //if(volspec(VolumeId,0,3) == 111) {CBEbotflag = 1; cout << "CBE bot hit!" << endl; }
				cout << "Hit is " << isig << " Energy Dep is " << Event->GetTrack(0)->GetEnergyDeposition(isig) << endl;
				//cout << GGeometryObject::GetVolumeIdFromRaw(1, 1, 1, 1) << endl;
				//cout << GetVolumeIdFromRaw(uint layer, uint row, uint module, uint channel) << endl;
				//cout << Crane::Calibration::CTrkCalib::GetVolumeIdFromRaw(uint layer, uint row, uint module, uint channel);
				//cout << Event->GetVolumeIdFromRaw(uint layer, uint row, uint module, uint channel);
				//cout << "layer? " << layer << endl;
			}
		}		
	}
	//cout << endl;
	//cout << "New Event " << i << endl;
	/*
	//Cuts are implemented in this chunk:
	if(Event->GetNTracks() == 1){ //First select the single track event
		singevents++;
		//cout << "Single Track Event " << endl;
		bool Umbflag = 0;
		bool CBEtopflag = 0;
		bool CBEbotflag = 0;

		CTrackRec* pt = Event->GetPrimaryTrack();
		uint pt_index = 0;
		if(pt == nullptr){cout << "Null pointer event number " << i << endl; }
		if(pt != nullptr){
   	        	for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;
                	pritracks++;	
		} 
	}
	*/
}

cout << "Number of single track events " << singevents << endl;
cout << "Number of primary track events " << pritracks << endl;
//cout << "Number of events that pass UMB, CBETop, CBEbot flags " << neventspass << endl; 
cout << endl << "I am done" << endl; 

}

