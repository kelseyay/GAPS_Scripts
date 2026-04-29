//I have been convinced not to do this

#include "KYtools.C"

void  Brec_Compare(){

int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events

//Filename
char FilenameRoot[400];
sprintf(FilenameRoot,"/home/kelsey/simulations/simdat/flight/251221/FPSI/starlink251221_0003_FindPrimaryStarIterative_rec.root");

//Prepare reconstronstructed event
CEventRec* FPSIEvent = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &FPSIEvent); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

TreeRec->GetEntry(0);

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
for(unsigned int i = 0; i < 10; i+=MainLoopScaleFactor){ //Let's try a smaller maximum number of events
//for(int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);
	//if( TreeRec->GetEntries() % (i+1) == 0){cout << "Time at Event " << i << " = " << FPSIEvent->GetEventTime() << endl;}

	CTrackRec* pt = FPSIEvent->GetPrimaryTrack();
	uint pt_index = 0;
	for( ; pt_index < FPSIEvent->GetNTracks(); pt_index++) if( FPSIEvent->GetTrack(pt_index)->IsPrimary() ) break;

	cout << "Event is " << i << " Rec Beta " << FPSIEvent->GetPrimaryBeta() << endl;
}

cout << endl << "I am done" << endl;

}
