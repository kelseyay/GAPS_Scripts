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
//sprintf(FilenameRoot,"/home/kelsey/simulations/data/MC/v2.1.0/GAPSSim_proton_TL1_1000det_FTFP_BERT_17_4242_FindPrimaryStarIterative.root");
sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/v.2.1.2/mu-_gaps_triggerlevel1_FTFP_BERT_1744342800_rec.root");

//Prepare MC event
CEventMc* Event = new CEventMc(); //New reconstructed event
TChain * TreeMC = new TChain("TreeMc"); //New TreeMC Tchain object (this is new to me)
//TTree * TreeMC = (TTree*)FilenameRoot->Get("TreeMC");
TreeMC->SetBranchAddress("Mc", &Event); //Set the branch address using Event (defined above)
TreeMC->Add(FilenameRoot);

//Maybe also need ot open the reconstructed event in "tandem" with the reconstructed event.
//Prepare reconstronstructed event
CEventRec* RecEvent = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
//TTree * TreeRec = (TTree*)FilenameRoot->Get("TreeRec");
TreeRec->SetBranchAddress("Rec", &RecEvent); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

//TGeoManager* geoMan = LoadGeometryRoot("/home/kelsey/simulations/test/ethernet241213_145/ethernet241213_1451_rec.root");
//TGeoManager* geoMan = LoadGeometryRoot("/home/kelsey/simulations/simdat/simrec/mu-_gaps_triggerlevel1_FTFP_BERT_HP_1721258929_rec.root");
//TGeoManager* geoMan = LoadGeometryRoot("/home/kelsey/simulations/data/MC/v3.0.0/anti_proton_gaps_triggerlevel2_FTFP_BERT_1759260144_rec.root");
//GGeometryMapPtr mcGeo = UnpackGeometry(geoMan, true, true);


//  CEventRec* reco_event_simu_ = new CEventRec(); //New CEventRec (SimpleDet class reference)
//  TTree * _recoTree_Simu_ = (TTree*)file_simu_mu->Get("TreeMC"); //New TTree in muon simulated file. TreeMC
//  _recoTree_Simu_->SetBranchAddress("Rec",&reco_event_simu_);

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;
int singevents = 0;
int pritracks = 0;
int neventspass = 0;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

//Now we can go over the loop

int PercentageStep = 5;
TreeMC->GetEntry(0);
TreeRec->GetEntry(0);

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeMC->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
for(unsigned int i = 0; i < 10; i+=MainLoopScaleFactor){ //Let's try a smaller maximum number of events
//for(int i = 0; i < TreeMC->GetEntries(); i+=MainLoopScaleFactor){
//
    TreeRec->GetEntry(i);
    TreeMC->GetEntry(i);

    cout << "Event is " << i << endl;
    cout << "Event ID? " << Event->GetEventId() << endl;
    cout << "Event Number? " << Event->GetEventNumber() << endl; //Gviz2D is for sure pulling Event number!
	//if( TreeMC->GetEntries() % (i+1) == 0){cout << "Time at Event " << i << " = " << Event->GetEventTime() << endl;}

	CTrackMc* pt = Event->GetPrimaryTrack();
	//CTrackRec* pt = Event->GetPrimaryTrack();
	uint pt_index = 0;
	for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

	bool Umbflag = 0;
	bool CBEtopflag = 0;
	bool CBEbotflag = 0;

	bool cool = 1;
	int inttrk = 0;

	cout << "Number of tracks " << Event->GetNTracks() << endl;

	if(Event->GetNTracks() > 0){ //Let's try picking at the track information: I am curious about energy deposits and DEF IN MC SPECIES!!!
	    //Mmk, got it! Of course the Gviz and Gviewer are code and code is readable! In the SimpleDet tools, the .cc files have the info :3
	    //If GetNTracks > 1, then there are tracks to iterate over!
		for(uint t = 0; t < Event->GetNTracks();t++){
		    //What hits needed to call this track worthwhile? More than one non-zero energy deposition, right?
		    for(unsigned int isig = 0; isig < Event->GetTrack(t)->GetEnergyDeposition().size(); isig++){
				if(Event->GetTrack(t)->GetEnergyDeposition(isig) > 0){ inttrk++; }
			}
			if(inttrk > 1){cool = 1;}
			inttrk = 0;

		    if(cool){
				cout << "Track is " << t << endl;
				cout << "GetTrackId()? " << Event->GetTrack(t)->GetTrackId()<< endl; //t is not the same as GetTrackId, interesting!
				    //Try to pick at some information?
				cout << "IsXray? " << Event->GetTrack(t)->IsXray() << endl;
				cout << "GetParentId()? " << Event->GetTrack(t)->GetParentId()<< endl;
				cout << "GetProcessType()? " << Event->GetTrack(t)->GetProcessType()<< endl;
				cout << "GetPdg()? " << Event->GetTrack(t)->GetPdg()<< endl;
				//cout << "IsMc()? " << Event->GetTrack(t)->IsMc()<< endl;
				//cout << "IsXray? " << Event->GetTrack(t)->IsXray() << endl;

				for(unsigned int isig = 0; isig < Event->GetTrack(t)->GetEnergyDeposition().size(); isig++){
				    //Loop over each energy deposition on each track
				    cout << "Energy deposition " << Event->GetTrack(t)->GetEnergyDeposition(isig) << " at Volume " << Event->GetTrack(t)->GetVolumeId(isig) << endl;
								//The volume ID is weird? Multiple different energy deposits at the same volumeID? Maybe that's a simulations thing...?
					//So there's a ton of "tracks" but so many just have one or a bunch of Energy depositions = 0...
					//How can there be a track also with only one point too haha?
					cool = 1;
				}
			}
		}
	}



	/* //Single track analysis!
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
				cout << "Hit is " << isig << " Energy Dep is " << Event->GetTrack(0)->GetVolumeId(isig); << endl;
				//cout << GGeometryObject::GetVolumeIdFromRaw(1, 1, 1, 1) << endl;
				//cout << GetVolumeIdFromRaw(uint layer, uint row, uint module, uint channel) << endl;
				//cout << Crane::Calibration::CTrkCalib::GetVolumeIdFromRaw(uint layer, uint row, uint module, uint channel);
				//cout << Event->GetVolumeIdFromRaw(uint layer, uint row, uint module, uint channel);
				//cout << "layer? " << layer << endl;
			}
		}
	}
	*/


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


	cout << endl; //End of event
}

//cout << "Number of events that pass UMB, CBETop, CBEbot flags " << neventspass << endl;
cout << endl << "I am done" << endl;

}
