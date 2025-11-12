#if !defined(__CINT__) || defined(__MAKECINT__)


#include <TSystem.h>
#include <TStyle.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TStopwatch.h>
#include <TCut.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TF1.h>
#include <TProfile.h>

#include <TSpline.h>

#include <TMatrixD.h>  // Includi la libreria per le matrici TMatrixD


#include "Math/GenVector/EulerAngles.h"
#include "Math/Vector3D.h"
using namespace ROOT::Math;

#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>


#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <bitset>
#include <iomanip>

using namespace std;

#include <GFileIO.hh>

#include "CCalib.hh"
#include <CEventRec.hh>
#include "CRawTrk.hh"
#include "GGeometryTools.hh"
using namespace Crane;
using namespace Reconstruction;
using namespace TrackFit;
using namespace Common;
using namespace Calibration;


#include "GetChain.C"


#define NLAYERS    7
#define NROWS      6
#define NMODULES   6
#define NDETECTORS 4
#define NSTRIPS    8
#define NCHANNELS  NSTRIPS*NDETECTORS

//  evaluate the common noise dividing the strips in group of module/NCN
//#define NCN 4 //detector
#define NCN 2 //half module
//#define NCN 1 //module




#endif

///////////////////////////////////////////
// Pedestals, Sigmas and Common noise
//
// flist = ascii file containing a list of names
// The script creates a TChain with files having path  ddir/name1+suffix , ddir/name2+suffix , etc...
// subcn = true(false) to evaluate and subtract the common noise
// nev = event to be used for pedestal evaluation
//
// The script write two files
// 1) flist-pedestals or flist-pedestals-cn
// 2) flsit-pedestals.root or flsit-pedestals-cn.root
// The former contains pedestals and sigmas to be used by DataCalibration
// The latter contains histograms (including covariance matrix)
//
void EvaluatePedestals(TString flist, int cnmod=0,Long64_t nev=9999999999,Long64_t evidmin=0, Long64_t evidmax=9999999999, TString ddir="",TString suffix="");//**
//
// some scripts to plot pedestal sigmas and differences between different files
//
void CheckPedestalFile(std::string file, double max=2050);//**
void CheckSigmaFile(std::string file, double max=100);//**
void DiffFiles(std::string file1,std::string file2, bool getsigmas=false);//**

void CheckMaskFile(std::string file, double max=2);//**
void CheckGainFile(std::string file, double max=100);//**


void CheckStripRate(TString flist, int nev=999999999, float fraction_max=1,  TString ddir="",TString suffix="",string inputfile="");
//
// Script to show covariance matrix for a given module
// indicated as module = 100*layer+10*row+module
//
void ShowRho(TFile* f,int module);//**

map<uint,double> GetPedestals(std::string file, bool getsigmas=false); // key = layer*10000 + row*1000 + module*100 + channel ;/


///////////////////////////////////////////
// Transfer functions

map<uint,TGraph*> GetTransferFunctions(std::string file, int channel = -1, std::string maskfile="", bool TESTLIB=true);

void CheckTFFile(std::string file, int detector = -1, std::string maskfile="");
void DiffTFFile(std::string file1,std::string file2, int detector = -1, bool TESTLIB=true);


///////////////////////////////////////////
// baseline

TTree* EvaluateBaseline(TString flist, TString fped , long nev = 999999999999, TString ddir="",TString suffix="");
void ShowBaseline(TFile* file, int detector);

///////////////////////////////////////////

TTree* ExtractRow(TString flist, int layer, int row, std::string file = "pedestals-sigmas-20240212-l6-l5-l4-l3-l2.txt", TString ddir="cal_data_level0/",TString suffix="-calib.root");


void DrawTitle(TVirtualPad *c1, TString title, double x=0.07, double y=0.93){
  if(!c1)return;
  c1->cd();  // c1 is the TCanvas
  TLatex *lat = new TLatex();
  lat->SetTextSize(0.03);
  lat->DrawLatexNDC(x,y,title.Data());
}
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////

TTree* ExtractRow(TString flist, int ll, int rr, std::string file , TString ddir,TString suffix){

  TChain *craw = GetChain(flist.Data(),ddir,"TreeRaw",suffix);
  Crane::Calibration::CRawTrk *rawtrk = new Crane::Calibration::CRawTrk();
  craw->SetBranchAddress("Trk", &rawtrk);

  float raw[NMODULES][NCHANNELS];
  float adc[NMODULES][NCHANNELS];
  float cn[NMODULES][NCHANNELS];
  // float adc[NMODULES*NCHANNELS];
  //float adc[192];
  int nm = NMODULES;
  int nc = NCHANNELS;
  float adc00;

  TFile* fi = new TFile(Form("out_%i.root",10*ll+rr),"recreate");
  TTree *tree = new TTree(Form("tree_%i",10*ll+rr),"tree");
  tree->Branch("raw",raw,Form("raw[%i][%i]/F",nm,nc));
  tree->Branch("adc",adc,Form("adc[%i][%i]/F",nm,nc));
  tree->Branch("cn",cn,Form("cn[%i][2]/F",nm));

  cout << endl << " N.entries "<<craw->GetEntries();


  map<uint,double> pedmap = GetPedestals(file);
  cout << endl << " map size "<<pedmap.size();


  float cov[NMODULES][NCHANNELS][NCHANNELS];
  float ave[NMODULES][NCHANNELS];
  float n[NMODULES][NCHANNELS];
  fill_n(&ave[0][0],NMODULES*NCHANNELS,0);
  fill_n(&cov[0][0][0],NMODULES*NCHANNELS*NCHANNELS,0);
  fill_n(&n[0][0],NMODULES*NCHANNELS,0);


  ////////////////////////////////////////////////////
  cout << endl << "OOOOOOO Loop over the events... "<<endl;
  ////////////////////////////////////////////////////
  for (uint ie=0; ie<craw->GetEntries(); ie++){
    craw->GetEntry(ie);

    fill_n(&raw[0][0],NMODULES*NCHANNELS,0);
    fill_n(&adc[0][0],NMODULES*NCHANNELS,0);
    fill_n(&cn[0][0],NMODULES*2,0);
    //fill_n(adc,NMODULES*NCHANNELS,0);

    int nval=0;
    for(uint i=0; i<(uint)rawtrk->adcdata.size(); i++){
      int layer  = rawtrk->layer[i];
      int row    = rawtrk->row[i];
      int module = rawtrk->module[i];
      int channel = rawtrk->channel[i];
      if(layer!=ll)continue;
      if(row!=rr)continue;
      nval++;
    }
    if(nval!=NMODULES*NCHANNELS){
      cout << endl << ie<<" n.val "<<nval<<" - NOT FULL - skip";
      continue;
    }
    for(uint i=0; i<(uint)rawtrk->adcdata.size(); i++){
      int layer  = rawtrk->layer[i];
      int row    = rawtrk->row[i];
      int module = rawtrk->module[i];
      int channel = rawtrk->channel[i];
      int half = (channel < (int)NCHANNELS/2 ? 0 : 1 );
      double val = (double)(rawtrk->adcdata[i]);
      if(layer!=ll)continue;
      if(row!=rr)continue;


      uint key = layer*10000+row*1000+module*100+channel;
      std::map<uint,double>::iterator it = pedmap.find(key);
      if (it != pedmap.end()){
	adc[module][channel] = (val-it->second);
	//	adc[module*NCHANNELS+channel] = (val-it->second);
      }else{
	cout << endl << "missing ped key "<<key;
      }
      if(module==0&&channel==0)adc00=val-it->second;


    }

    for(int im = 0; im < NMODULES; im++){
      for(int ic = 0; ic < NCHANNELS; ic++){

	for(int jc = 0; jc < NCHANNELS; jc++){
	}
      }
    }

    tree->Fill();
  }

  fi->cd();
  tree->Write();
  fi->Close();
  return tree;
};
/////////////////////////////////////
void EvaluatePedestals(TString flist, int cnmod ,Long64_t nev,Long64_t evidmin, Long64_t evidmax, TString ddir,TString suffix, std::string maskfile){


  bool subcn = cnmod>0;

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TChain *craw = GetChain(flist.Data(),ddir,"TreeRaw",suffix);
  TChain *ccal = GetChain(flist.Data(),ddir,"TreeRec",suffix);
  craw->AddFriend(ccal);

  //----------------------------------------- reconstructed event
  CEventRec *evrec= new CEventRec();
  craw->SetBranchAddress("Rec", &evrec);
  //----------------------------------------- raw event
  Crane::Calibration::CRawTrk *rawtrk = new Crane::Calibration::CRawTrk();
  craw->SetBranchAddress("Trk", &rawtrk);
  //----------------------------------------- Prepare mask file
  CCalib *data_calib = NULL;
  data_calib = new  Crane::Calibration::CCalib();
  data_calib->GetTrkCalib().SetMasks(maskfile);

  //-----------------------------------------
  // define histograms
  //
  TH2F *hped = (TH2F*)gROOT->FindObject("hped");
  if(hped)delete hped;
  hped = new TH2F("hped","",NCHANNELS*NROWS,0,NCHANNELS*NROWS,NLAYERS*NMODULES,0,NLAYERS*NMODULES);
  hped->GetXaxis()->SetTitle(Form("row*%i+channel",NCHANNELS));
  hped->GetYaxis()->SetTitle(Form("layer*%i+module",NMODULES));
  hped->SetStats(0);

  TH1F* hpeddist = new TH1F("hpeddist","PED",300,0,2050);

  TH2F *hsig = (TH2F*)gROOT->FindObject("hsig");
  if(hsig)delete hsig;
  hsig = new TH2F("hsig","",NCHANNELS*NROWS,0,NCHANNELS*NROWS,NLAYERS*NMODULES,0,NLAYERS*NMODULES);
  hsig->GetXaxis()->SetTitle(Form("row*%i+channel",NCHANNELS));
  hsig->GetYaxis()->SetTitle(Form("layer*%i+module",NMODULES));
  hsig->SetStats(0);

  TH1F* hsigdist = new TH1F("hsigdist","SIG",500,0,100);

  TH2F *hped_0 = (TH2F*)gROOT->FindObject("hped_0");
  if(hped_0)delete hped_0;
  hped_0 = new TH2F("hped_0","",NCHANNELS*NROWS,0,NCHANNELS*NROWS,NLAYERS*NMODULES,0,NLAYERS*NMODULES);
  hped_0->GetXaxis()->SetTitle(Form("row*%i+channel",NCHANNELS));
  hped_0->GetYaxis()->SetTitle(Form("layer*%i+module",NMODULES));
  hped_0->SetStats(0);

  TH1F* hped_0dist = new TH1F("hped_0dist","PED",300,0,2050);

  TH2F *hsig_0 = (TH2F*)gROOT->FindObject("hsig_0");
  if(hsig_0)delete hsig_0;
  hsig_0 = new TH2F("hsig_0","",NCHANNELS*NROWS,0,NCHANNELS*NROWS,NLAYERS*NMODULES,0,NLAYERS*NMODULES);
  hsig_0->GetXaxis()->SetTitle(Form("row*%i+channel",NCHANNELS));
  hsig_0->GetYaxis()->SetTitle(Form("layer*%i+module",NMODULES));
  hsig_0->SetStats(0);

  TH1F* hsig_0dist = new TH1F("hsig_0dist","SIG",300,0,100);

  TGraph *gsig0 = new TGraph(); gsig0->SetName("gsig0");
  TGraph *gsig  = new TGraph(); gsig->SetName("gsig");
  TGraph *gsig3_min  = new TGraph(); gsig3_min->SetName("gsig3_min");
  TGraph *gsig3_best  = new TGraph(); gsig3_best->SetName("gsig3_best");


  TH2F *hgain = (TH2F*)gROOT->FindObject("hgain");
  if(hgain)delete hgain;
  hgain = new TH2F("hgain","",NCHANNELS*NROWS,0,NCHANNELS*NROWS,NLAYERS*NMODULES,0,NLAYERS*NMODULES);
  hgain->GetXaxis()->SetTitle(Form("row*%i+channel",NCHANNELS));
  hgain->GetYaxis()->SetTitle(Form("layer*%i+module",NMODULES));
  hgain->SetStats(0);

  TH1F* hgaindist = new TH1F("hgaindist","GAIN",300,-15,15);




  TH2F* hfull = new TH2F("hfull","NZS modules",NROWS,0,NROWS,NLAYERS*NMODULES,0,NLAYERS*NMODULES);
  hfull->SetStats(0);


  TH2F *hpulser = (TH2F*)gROOT->FindObject("hpulser");
  if(hpulser)delete hpulser;
  hpulser = new TH2F("hpulser","",NCHANNELS*NROWS,0,NCHANNELS*NROWS,NLAYERS*NMODULES,0,NLAYERS*NMODULES);
  hpulser->GetXaxis()->SetTitle(Form("row*%i+channel",NCHANNELS));
  hpulser->GetYaxis()->SetTitle(Form("layer*%i+module",NMODULES));
  hpulser->SetStats(0);




  vector<double> sigrank;
  vector<TH1F*> hsok;
  vector<TH1F*> hsnotok;
  sigrank.push_back(2.);
  sigrank.push_back(8.);
  sigrank.push_back(100.);
  for(int i=0; i<=sigrank.size();i++){
    float rank = (i<sigrank.size() ? sigrank.at(i) : 9999);
    hsok.push_back( new TH1F(Form("hsok_%i",i),Form("ADC-PED ok (SIG<%f)",rank),500,-200,1000) );
    hsnotok.push_back( new TH1F(Form("hsnotok_%i",i),Form("ADC-PED !ok  (SIG<%f)",rank),500,-200,1000) );
  }

  TH1F* hsall0[NLAYERS][NROWS][NMODULES]; /// all, before cn correction
  TH1F* hsped[NLAYERS][NROWS][NMODULES]; ///used for ped evaluation
  TH1F* hsmip[NLAYERS][NROWS][NMODULES]; ///excluded from ped evaluation
  for(int il=0; il<NLAYERS; il++){
    for(int ir=0; ir<NROWS; ir++){
      for(int im=0; im<NMODULES; im++){
        hsall0[il][ir][im] =  new TH1F(Form("hsall0_%i%i%i",il,ir,im),Form("Module %i%i%i",il,ir,im),600,-200,1000);
        hsped[il][ir][im] =  new TH1F(Form("hsped_%i%i%i",il,ir,im),Form("Module %i%i%i",il,ir,im),600,-200,1000);
        hsmip[il][ir][im] =  new TH1F(Form("hsmip_%i%i%i",il,ir,im),Form("Module %i%i%i",il,ir,im),600,-200,1000);
      }
    }
  }

  // -----------------------------------------------------
  // define vectors for pedestal evaluation
  //
  double  adc[NLAYERS][NROWS][NMODULES][NCHANNELS];
  fill_n( &adc[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);

  double  sum[NLAYERS][NROWS][NMODULES][NCHANNELS];
  double sum2[NLAYERS][NROWS][NMODULES][NCHANNELS];
  int       n[NLAYERS][NROWS][NMODULES][NCHANNELS];
  fill_n( &sum[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);
  fill_n(&sum2[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);
  fill_n(   &n[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);

  double  ped[NLAYERS][NROWS][NMODULES][NCHANNELS];
  double  sig[NLAYERS][NROWS][NMODULES][NCHANNELS];
  fill_n( &ped[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);
  fill_n( &sig[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);

  bool disconnected[NLAYERS][NROWS][NMODULES][NCHANNELS];
  fill_n( &disconnected[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);

  double cov[NLAYERS][NROWS][NMODULES][NCHANNELS][NCHANNELS];
  fill_n( &cov[0][0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS*NCHANNELS,0);

  // double  sumra[NLAYERS][NROWS][NMODULES][NCHANNELS];
  // int     nsumra[NLAYERS][NROWS][NMODULES][NCHANNELS];
  // fill_n( &sumra[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,1);
  // fill_n( &nsumra[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,1);
  double  gain[NLAYERS][NROWS][NMODULES][NCHANNELS];
  fill_n( &gain[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,1);

  double  sig0[NLAYERS][NROWS][NMODULES][NCHANNELS];
  fill_n( &sig0[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);

  vector<int> chmin[NLAYERS][NROWS][NMODULES];

  double  sumcn[NLAYERS][NROWS][NMODULES][NCN];
  double  sumcn2[NLAYERS][NROWS][NMODULES][NCN];
  int     nsumcn[NLAYERS][NROWS][NMODULES][NCN];
  fill_n( &sumcn[0][0][0][0],NLAYERS*NROWS*NMODULES*NCN,0);
  fill_n( &sumcn2[0][0][0][0],NLAYERS*NROWS*NMODULES*NCN,0);
  fill_n( &nsumcn[0][0][0][0],NLAYERS*NROWS*NMODULES*NCN,0);

  bool module_ok[NLAYERS][NROWS][NMODULES][2];
  int module_n[NLAYERS][NROWS][NMODULES][2];


  double cn[NLAYERS][NROWS][NMODULES][NCN];
  int    ncn[NLAYERS][NROWS][NMODULES][NCN];
  // double min[NLAYERS][NROWS][NMODULES][NCN];

  // -----------------------------------------------------
  // define other histograms
  //

  TH1F *hhh = new TH1F("hhh","n*sig",300,-10,100);

  double sigcn[NLAYERS][NROWS][NMODULES][NCN];
  TH1F* hcn[NLAYERS][NROWS][NMODULES][NCN];
  TH1F* hs[NLAYERS][NROWS][NMODULES][NCHANNELS];
  // TH1F* hmin[NLAYERS][NROWS][NMODULES][NCN];
  for(int il=0; il<NLAYERS; il++){
    for(int ir=0; ir<NROWS; ir++){
      for(int im=0; im<NMODULES; im++){
  	for(int ih=0; ih<NCN; ih++){
  	  hcn[il][ir][im][ih] = NULL;
	  sigcn[il][ir][im][ih] = 0;
  	  // hmin[il][ir][im][ih] = NULL;
  	}
  	for(int ic=0; ic<NCHANNELS; ic++){
  	  hs[il][ir][im][ic] = NULL;
	}
      }
    }
  }

  TH2F *hsigcn = (TH2F*)gROOT->FindObject("hsigcn");
  if(hsigcn)delete hsigcn;
  hsigcn = new TH2F("hsigcn","",NCN*NROWS,0,NCN*NROWS,NLAYERS*NMODULES,0,NLAYERS*NMODULES);
  hsigcn->GetXaxis()->SetTitle(Form("row*%i+detector",NCN));
  hsigcn->GetYaxis()->SetTitle(Form("layer*%i+module",NMODULES));
  hsigcn->SetStats(0);

  TH1F* hsigcndist = new TH1F("hsigcndist","SIG CN",300,0,100);

  ////////////////////////////////////////////////////
  cout << endl << "OOOOOOO Loop over the events... (0)"<<endl;
  ////////////////////////////////////////////////////

  cout << endl << " N.tot entries "<<craw->GetEntries();
  nev = TMath::Min((Long64_t)nev,(Long64_t)craw->GetEntries());
  //  cout << endl << " using "<<nev<<" events";
  int cnt=0;
  for (uint ie=0; ie<craw->GetEntries(); ie++){

    if(cnt==nev)break;

    craw->GetEntry(ie);

    if(rawtrk->eventid < evidmin )continue;
    if(rawtrk->eventid > evidmax )continue;
    cnt++;

    // reset vectors
    fill_n(&module_n[0][0][0][0],NLAYERS*NROWS*NMODULES*2,0);
    fill_n(&module_ok[0][0][0][0],NLAYERS*NROWS*NMODULES*2,true);
    fill_n( &adc[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);


    //-----------------------------------------
    // first evaluation of pedestal and sigmans
    // without the exlcusion of particle signals
    // without common noise subtraction
    //-----------------------------------------

    //
    // counts transmitted channels for each half module
    // to identify events with zero suppression
    //
    for(uint i=0; i<(uint)rawtrk->adcdata.size(); i++){
      int layer  = rawtrk->layer[i];
      int row    = rawtrk->row[i];
      int module = rawtrk->module[i];
      int channel = rawtrk->channel[i];
      int half = (channel < (int)NCHANNELS/2 ? 0 : 1 );
      module_n[layer][row][module][half]++;
    }
    //
    // fill adc vector for this event
    // if the event is acquired with zero suppression, it is exlcuded
    //
    for(uint i=0; i<(uint)rawtrk->adcdata.size(); i++){
      int layer  = rawtrk->layer[i];
      int row    = rawtrk->row[i];
      int module = rawtrk->module[i];
      int channel = rawtrk->channel[i];
      int half = (channel < (int)NCHANNELS/2 ? 0 : 1 );

      bool FULL = module_n[layer][row][module][half] == (int)NCHANNELS/2;
      if(!FULL)continue;//go to next event

      double val = (double)(rawtrk->adcdata[i]);
      adc[layer][row][module][channel] = val;

    }

    // -------------------
    // increment counters
    // ------------------
    for(int il = 0 ; il<NLAYERS; il++){
      for(int ir = 0 ; ir<NROWS; ir++){
	for(int im = 0 ; im<NMODULES; im++){

	  // check if the module is transmitted full
	  if(module_n[il][ir][im][0] != (int)NCHANNELS/2)continue;
	  if(module_n[il][ir][im][1] != (int)NCHANNELS/2)continue;
          hfull->Fill(ir,il*NMODULES+im);
	  //ok, is full. increment counters
	  for(int ic=0; ic<NCHANNELS;ic++){
	    sum[il][ir][im][ic]+=adc[il][ir][im][ic];
	    sum2[il][ir][im][ic]+=adc[il][ir][im][ic]*adc[il][ir][im][ic];
	    n[il][ir][im][ic]++;
	  }
	}
      }
    }


  }//loop over events
  std::cout<<std::endl<<" selected events "<<cnt;

  ////////////////////////////////////////////////////
  cout << endl << "OOOOOOO Evaluate PED SIG (0)";
  ////////////////////////////////////////////////////
  // evaluate pedestals and sigmas
  for(int il = 0 ; il<NLAYERS; il++){
    for(int ir = 0 ; ir<NROWS; ir++){
      for(int im = 0 ; im<NMODULES; im++){
	bool OK = false;
	for(int ic = 0 ; ic<NCHANNELS; ic++){
	  if( n[il][ir][im][ic] == 0  )continue;
	  ped[il][ir][im][ic] = sum[il][ir][im][ic]/n[il][ir][im][ic];
	  sig[il][ir][im][ic] = sum2[il][ir][im][ic]/n[il][ir][im][ic];
	  sig[il][ir][im][ic] -= ped[il][ir][im][ic]*ped[il][ir][im][ic];
	  sig[il][ir][im][ic] = TMath::Sqrt( sig[il][ir][im][ic] );
	  OK=true;
	}
      }
    }
  }
  // fill the histogramns
  for(int il = 0 ; il<NLAYERS; il++){
    for(int ir = 0 ; ir<NROWS; ir++){
      for(int im = 0 ; im<NMODULES; im++){
	for(int ic = 0 ; ic<NCHANNELS; ic++){
	  if( n[il][ir][im][ic] == 0  )continue;
	  hped_0->Fill(ir*NCHANNELS+ic,il*NMODULES+im,ped[il][ir][im][ic]);
	  hped_0dist->Fill(ped[il][ir][im][ic]);
	  hsig_0->Fill(ir*NCHANNELS+ic,il*NMODULES+im,sig[il][ir][im][ic]);
	  hsig_0dist->Fill(sig[il][ir][im][ic]);

	}
      }
    }
  }

  //-----------------------------------------
  // iterative evaluation of pedestal and sigmans
  // with exclusion of particle signals
  // and common noise subtraction (if required)
  //-----------------------------------------


  vector<TH2F*> vhrho;
  vector<TH1F*> vhsig;
  vector<TH2F*> vhrho_ex;
  vector<TH1F*> vhsig_ex;

  //
  // define some parameters to control the iteration
  //
  // double sigcut = 5;
  // double sigcut = 4;
  // double sigcut = 10;
  double sigcut = 10;//5;
  double mipcut = 200;//cut on signal to identify mips
  double sigmin = 1.7;//cut on sig to tag disconnected channels // 1.7 ??
  double sigcncut = sigcut;//4; // cut to exlude event with anomalous cn fluctuations
  int niter = 4;//tot n.iterations
  int niter_cn = 3;//subtract common noise starting from this iteration
  int ncnmin = 4;//minimum number of channels good to evaluate baseline

  for (uint it=1 ; it<=niter; it++){ // iterate

    bool DOCNEVAL = subcn&(it>=niter_cn);

    ////////////////////////////////////////////////////
    cout << endl << "OOOOOOO Loop over the events... ("<<it<<")";
    if(DOCNEVAL)cout<<" --subtract CN--";
    cout << endl;
    ////////////////////////////////////////////////////
    // reset vectors
    fill_n( &sum[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);
    fill_n(&sum2[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);
    fill_n(   &n[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);

    fill_n(&module_n[0][0][0][0],NLAYERS*NROWS*NMODULES*2,0);
    fill_n(&module_ok[0][0][0][0],NLAYERS*NROWS*NMODULES*2,true);

    fill_n( &cov[0][0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS*NCHANNELS,0);

    ////////////////////////////////////////////////////////
    // loop over the events
    ////////////////////////////////////////////////////////
    int cnt=0;
    for (uint ie=0; ie<craw->GetEntries(); ie++){

      if(cnt==nev)break;
      craw->GetEntry(ie);
      if(rawtrk->eventid < evidmin )continue;
      if(rawtrk->eventid > evidmax )continue;
      cnt++;

      fill_n(&module_n[0][0][0][0],NLAYERS*NROWS*NMODULES*2,0);
      fill_n(&module_ok[0][0][0][0],NLAYERS*NROWS*NMODULES*2,true);
      fill_n(&cn[0][0][0][0],NLAYERS*NROWS*NMODULES*NCN,0);
      fill_n(&ncn[0][0][0][0],NLAYERS*NROWS*NMODULES*NCN,0);

      fill_n( &adc[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);

      ////////////////////////////////////////////////////////
      // Tag modules with particle signals
      // count transmitted channels to find full modules
      // increment common noise counter
      ////////////////////////////////////////////////////////
      for(uint i=0; i<(uint)rawtrk->adcdata.size(); i++){

	int layer  = rawtrk->layer[i];
	int row    = rawtrk->row[i];
	int module = rawtrk->module[i];
	int channel = rawtrk->channel[i];
	int half = (channel < (int)NCHANNELS/2 ? 0 : 1 );
	int detector = (int)(channel/8);
	int cnindex = 0;//(NCN==2 ? half : detector);
        if(NCN==2) cnindex = half;
        if(NCN==4) cnindex = detector;
 	double val = (double)(rawtrk->adcdata[i]) - ped[layer][row][module][channel];
	double baseline = 0; //do not refer tominimum
	//	dounle baseline = min[layer][row][module][cnindex];
	double nsig =  (val - baseline )/sig[layer][row][module][channel];

	module_n[layer][row][module][half]++;

        //
        // condition to identify particle signals
        //
	bool ISHIT = false;
	if(TMath::Abs(nsig)>sigcut) ISHIT=true;
	if(val > mipcut) ISHIT=true;
	if(ISHIT) module_ok[layer][row][module][0]=false;
	if(ISHIT) module_ok[layer][row][module][1]=false;
        //
	//increment baseline counter
        //
	if( DOCNEVAL &&
            //            !ISHIT && //ATTENZIONE !!!
	    !disconnected[layer][row][module][channel] &&
	    true){
	  cn[layer][row][module][cnindex] += val ;
	  ncn[layer][row][module][cnindex]++;
	}
	if(it==niter && !disconnected[layer][row][module][channel])hhh->Fill(nsig);

      }

      ////////////////////////////////////////////////////////
      // Evaluate baseline
      ////////////////////////////////////////////////////////
      for(int il=0; il<NLAYERS; il++){
	for(int ir=0; ir<NROWS; ir++){
	  for(int im=0; im<NMODULES; im++){
	    for(int ih=0; ih<NCN; ih++){
	      if( ncn[il][ir][im][ih]>=ncnmin )
		cn[il][ir][im][ih] /= ncn[il][ir][im][ih];
	      else
		cn[il][ir][im][ih] = 0;
	      if( ncn[il][ir][im][ih] == 0 ) continue;
	      //
	      // check baseline shift
              // (exclude events with too large baseline fluctuations...)
	      //
	      if(sigcn[il][ir][im][ih]==0)continue;
	      if( TMath::Abs(cn[il][ir][im][ih]) > sigcncut*sigcn[il][ir][im][ih]  ){
	      	module_ok[il][ir][im][0]=false;
	      	module_ok[il][ir][im][1]=false;
	      }
	      //cout << endl << TMath::Abs(cn[il][ir][im][ih])<<" "<<sigcn[il][ir][im][ih];
	      //
	    }
	  }
	}
      }


      ////////////////////////////////////////////////////////
      // Increment subtract baseline
      ////////////////////////////////////////////////////////
      for(uint i=0; i<(uint)rawtrk->adcdata.size(); i++){
	int layer  = rawtrk->layer[i];
	int row    = rawtrk->row[i];
	int module = rawtrk->module[i];
	int channel = rawtrk->channel[i];
	int half = (channel < (int)NCHANNELS/2 ? 0 : 1 );
	int detector = (int)(channel/8);
	double val = (double)(rawtrk->adcdata[i]);
	bool FULL = module_n[layer][row][module][half] == (int)NCHANNELS/2;
	if(!FULL)continue;

	int irank = sigrank.size();
	for(int ir=sigrank.size()-1; ir>=0; ir--){
	  if(sig[layer][row][module][channel]<sigrank[ir]){
	    irank = ir;
	  }
	}

	///////////////////////////////////////////
	int cnindex = 0;//(NCN==2 ? half : detector);
        if(NCN==2) cnindex = half;
        if(NCN==4) cnindex = detector;
	double cncn = 0;
	if(
	   !disconnected[layer][row][module][channel] &&
	   true){
	  cncn = cn[layer][row][module][cnindex]*gain[layer][row][module][channel]  ;
	}
	///////////////////////////////////////////
	adc[layer][row][module][channel] = val-cncn;
	///////////////////////////////////////////

	if(it==niter){ //@last iteration

	  if( module_ok[layer][row][module][half] ){

            hsped[layer][row][module]->Fill(adc[layer][row][module][channel]-ped[layer][row][module][channel]);
	    hsok[irank]->Fill(adc[layer][row][module][channel]-ped[layer][row][module][channel]);
	    if(!hs[layer][row][module][channel]) hs[layer][row][module][channel] = new TH1F(Form("hs_%i_%i",100*layer+10*row+module,channel),"R-PED-B",100,-50,50);
	    hs[layer][row][module][channel]->Fill(adc[layer][row][module][channel]-ped[layer][row][module][channel] );

	  }else{

            hsmip[layer][row][module]->Fill(adc[layer][row][module][channel]-ped[layer][row][module][channel]);
	    hsnotok[irank]->Fill(adc[layer][row][module][channel]-ped[layer][row][module][channel]);

	  }
	}//end last iteration condition

        if(it==niter_cn-1){ //@ last iteration before cn correction
          hsall0[layer][row][module]->Fill(adc[layer][row][module][channel]-ped[layer][row][module][channel]);
          // cout << endl << adc[layer][row][module][channel]-ped[layer][row][module][channel];
        }


      } // end loop over hits

      ////////////////////////////////////////////////////////
      // increment PED SIG COV cunters
      ////////////////////////////////////////////////////////
       for(int il = 0 ; il<NLAYERS; il++){
	for(int ir = 0 ; ir<NROWS; ir++){
	  for(int im = 0 ; im<NMODULES; im++){
	    // check if the module is transmitted full
	    if(module_n[il][ir][im][0] != (int)NCHANNELS/2)continue;
	    if(module_n[il][ir][im][1] != (int)NCHANNELS/2)continue;
	    //ok, is full. increment counters

	    if( module_ok[il][ir][im][0] &&
		module_ok[il][ir][im][1] &&
		true){
              //-----------------------
	      //--> GOOD for pedestals
              //-----------------------

	      for(int ic=0; ic<NCHANNELS;ic++){
		for(int jc=0; jc<NCHANNELS;jc++){
		  cov[il][ir][im][ic][jc] += adc[il][ir][im][ic]*adc[il][ir][im][jc];//
		}
		sum[il][ir][im][ic]+=adc[il][ir][im][ic];
		sum2[il][ir][im][ic]+=adc[il][ir][im][ic]*adc[il][ir][im][ic];
		n[il][ir][im][ic]++;

	      }


	      if(DOCNEVAL){
		for(int icn = 0; icn<NCN; icn++){
		  //--------------------------- fill CN histos
		  if(!hcn[il][ir][im][icn]) hcn[il][ir][im][icn]= new TH1F(Form("hcn_%i_%i",100*il+10*ir+im,icn),"Baseline (pedestals events)",1000,-50,50);
		  hcn[il][ir][im][icn]->Fill( cn[il][ir][im][icn] );
                  sumcn[il][ir][im][icn]  += cn[il][ir][im][icn];
                  sumcn2[il][ir][im][icn] += TMath::Power(cn[il][ir][im][icn],2);;
                  nsumcn[il][ir][im][icn]++;

		}
	      }



	    }else{

              //-----------------------
	      //--> NOT GOOD for pedestals
              //-----------------------

	    }

	  }
	}
      }

    }//end loop over events

    ////////////////////////////////////////////////////
    cout << endl << "OOOOOOO Evaluate PED SIG ("<<it<<")";
    ////////////////////////////////////////////////////
    fill_n( &ped[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);
    fill_n( &sig[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);
    fill_n( &disconnected[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);

    for(int il = 0 ; il<NLAYERS; il++){
      for(int ir = 0 ; ir<NROWS; ir++){
	for(int im = 0 ; im<NMODULES; im++){

	  ////////////////// CN
          for(int ih=0; ih<NCN; ih++){
            sigcn[il][ir][im][ih] = 0;//init
            if( !hcn[il][ir][im][ih] )continue;
            //	    sigcn[il][ir][im][ih]=TMath::Abs(hcn[il][ir][im][ih]->GetRMS());
            double ave = sumcn[il][ir][im][ih]/nsumcn[il][ir][im][ih];
	    sigcn[il][ir][im][ih] = sumcn2[il][ir][im][ih]/nsumcn[il][ir][im][ih];
            sigcn[il][ir][im][ih]-= ave*ave;
            sigcn[il][ir][im][ih] = TMath::Sqrt( sigcn[il][ir][im][ih] );
	    if(it<niter)hcn[il][ir][im][ih]->Reset();
	  }

	  ////////////////// PED & SIG
	  for(int ic = 0 ; ic<NCHANNELS; ic++){


	    if( n[il][ir][im][ic] == 0  )continue;


	    ped[il][ir][im][ic] = sum[il][ir][im][ic]/n[il][ir][im][ic];
	    sig[il][ir][im][ic] = sum2[il][ir][im][ic]/n[il][ir][im][ic];
	    sig[il][ir][im][ic] -=  ped[il][ir][im][ic]*ped[il][ir][im][ic];
	    sig[il][ir][im][ic] = TMath::Sqrt( sig[il][ir][im][ic] );
            disconnected[il][ir][im][ic] = (sig[il][ir][im][ic]<sigmin);
	    if(it==niter_cn-1){
              sig0[il][ir][im][ic] = sig[il][ir][im][ic];//copy
            }
            // if(disconnected[il][ir][im][ic]){
            //   cout << endl << "disconnected "<<il<<ir<<im<<" - "<<ic;
            // }
	    if( n[il][ir][im][ic] < 100 )
              cout << endl <<" MODULE "<<100*il+10*ir+im<<" ch "<<ic<<"  --- PED evaluated with n.events "<<n[il][ir][im][ic];
	  }

        }
      }
    }


    //---------------------------------------------------
    // Do things to be done aftefirst baseline evaluation
    //---------------------------------------------------
    fill_n( &gain[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,1);
    if(it == niter_cn && cnmod == 2){
      ////////////////////////////////////////////////////
      cout << endl << "OOOOOOO Evaluate GAIN ("<<it<<")";
      ////////////////////////////////////////////////////
      double avsig2[NCN];
      int    navsig2[NCN];
      fill_n( &avsig2[0],NCN,0);
      fill_n( &navsig2[0],NCN,0);
      for(int il = 0 ; il<NLAYERS; il++){
        for(int ir = 0 ; ir<NROWS; ir++){
          for(int im = 0 ; im<NMODULES; im++){
            // for(int ic = 0 ; ic<NCHANNELS; ic++){
            //   int half = (ic < (int)NCHANNELS/2 ? 0 : 1 );
            //   int detector = (int)(ic/8);
            //   int cnindex = 0;
            //   if(NCN==2) cnindex = half;
            //   if(NCN==4) cnindex = detector;
            //   if(sig[il][ir][im][ic]==0)continue;
            //   if(disconnected[il][ir][im][ic])continue;
            //   avsig2[cnindex] += sig[il][ir][im][ic]*sig[il][ir][im][ic];
            //   navsig2[cnindex]++;
            // }
            // for(int ih=0; ih<NCN;ih++)if(navsig2[ih]>0)avsig2[ih] /= navsig2[ih]*navsig2[ih];
            //            for(int ih=0; ih<NCN;ih++)cout<<endl<<"**** "<<avsig2[ih] ;
            for(int ic = 0 ; ic<NCHANNELS; ic++){
              int half = (ic < (int)NCHANNELS/2 ? 0 : 1 );
              int detector = (int)(ic/8);
              int cnindex = 0;
              if(NCN==2) cnindex = half;
              if(NCN==4) cnindex = detector;
              if(sig[il][ir][im][ic]==0)continue;
              if(sigcn[il][ir][im][cnindex]==0)continue;
              if(disconnected[il][ir][im][ic])continue;
              double varbase = sigcn[il][ir][im][cnindex]*sigcn[il][ir][im][cnindex] - avsig2[cnindex];
              double G = 0; //Gain part right here!
              //double GG = 0; //Testing the Gigi Gain

	      //Previous Elena method!
	      //G = sig0[il][ir][im][ic]*sig0[il][ir][im][ic] - sig[il][ir][im][ic]*sig[il][ir][im][ic];
              //G /= (2 *  varbase);
              //G += 0.5;

	      G =  ((NCHANNELS/NCN) - 1)*sig0[il][ir][im][ic]*sig0[il][ir][im][ic] - (NCHANNELS/NCN)* sig[il][ir][im][ic]*sig[il][ir][im][ic];
	      G /= (2 * (NCHANNELS/NCN) * varbase);
	      G += 0.5;

	      //cout << "Elena method G: " << G << endl;
	      //cout << "Gianluigi method GG: " << GG << endl << endl;

              gain[il][ir][im][ic] = G;
              hgaindist->Fill(gain[il][ir][im][ic]);
              hgain->Fill(ir*NCHANNELS+ic,il*NMODULES+im,gain[il][ir][im][ic]);
              // cout << endl <<" var(x)   "<< sig0[il][ir][im][ic]*sig0[il][ir][im][ic];
              // cout <<        " var(x-b) "<<sig[il][ir][im][ic]*sig[il][ir][im][ic];
              // cout <<        " var(b)   "<<sigcn[il][ir][im][cnindex]*sigcn[il][ir][im][cnindex];
              // cout <<        " gain     "<< gain[il][ir][im][ic] ;
              if(G<0){
                cout << endl<< "negative gain     "<< gain[il][ir][im][ic]<<" "<<il<<ir<<im<<" - "<<ic<<" - sig "<<sig[il][ir][im][ic] ;
              }

            }
          }
        }
      }
    }


    //---------------------------------------------
    // Do things to be done @ last iteration only
    //---------------------------------------------
    if(it==niter){

      ////////////////////////////////////////////////////
      cout << endl << "OOOOOOO Evaluate covariance matrix ("<<it<<")";
      ////////////////////////////////////////////////////


      for(int il = 0 ; il<NLAYERS; il++){
        for(int ir = 0 ; ir<NROWS; ir++){
          for(int im = 0 ; im<NMODULES; im++){

            int module = 100*il+10*ir+im;

            ////////////////// correlation factor GOOD events

            TH2F* h = new TH2F(Form("hrho_%i",module),Form("RHO module %i",module),NCHANNELS,0,NCHANNELS,NCHANNELS,0,NCHANNELS);
            h->SetStats(0);
            TH1F* hh = new TH1F(Form("hsig_%i",module),Form("SIG module %i",module),NCHANNELS,0,NCHANNELS);
            for(int ic = 0 ; ic<NCHANNELS; ic++){
              if( n[il][ir][im][ic] == 0  )continue;
              hh->Fill(ic,sig[il][ir][im][ic]);
              for(int jc = 0 ; jc<NCHANNELS; jc++){
                if( n[il][ir][im][jc] == 0  )continue;
                cov[il][ir][im][ic][jc] /= n[il][ir][im][ic];
                cov[il][ir][im][ic][jc] -= ped[il][ir][im][ic]*ped[il][ir][im][jc];
                h->Fill(ic,jc,cov[il][ir][im][ic][jc]/(sig[il][ir][im][ic]*sig[il][ir][im][jc]));
              }
            }
            vhrho.push_back(h);
            vhsig.push_back(hh);
          }
        }
      }

      if(cnmod == 2){
        ////////////////////////////////////////////////////
        cout << endl << "OOOOOOO Evaluate best channel to be pulsed ("<<it<<")";
        ////////////////////////////////////////////////////

        for(int il = 0 ; il<NLAYERS; il++){
          for(int ir = 0 ; ir<NROWS; ir++){
            for(int im = 0 ; im<NMODULES; im++){
              /////////////////////////////////////////// gains - intrinsic - c
              const int N = NCHANNELS/NCN;
              double sigma_n[N];
              double cgain[N];
              double sigma_c;

              for(int ih=0; ih<NCN; ih++){	    //1-2-4
                int ic_best = -1;
                int i_best  = -1;
                int ic_min  = -1;
                int i_min   = -1;
                double val_best = 9999999999;
                sigma_c = sigcn[il][ir][im][ih]; //RMS(baseline)
                //------------------------------------- search minimum/best
                for(int i=0; i<N; i++){
                  uint ic = ih*N+i;
                  if(disconnected[il][ir][im][ic])continue;
                  if(!data_calib->GetTrkCalib().IsGood(il,ir,im,ic))continue;
                  if(data_calib->GetTrkCalib().IsGood(il,ir,im,ic)) cout << "Good channel" << il << ir << im << ic << endl;
                  sigma_n[i] = sig[il][ir][im][ic];//intrinsic noise
                  cgain[i]   = gain[il][ir][im][ic];//gain
                  if( sigma_n[i]/cgain[i] < val_best ){
                    ic_min   = ic;
                    i_min    = i;
                    val_best = sigma_n[i]/cgain[i];
                    if( sigma_n[i]/cgain[i] < sigma_c ){
                      ic_best = ic;
                      i_best  = i;
                    }
                  }
                }
                // fill histo
                for(int i=0; i<N; i++){
                  uint ic = ih*N+i;
                  if(ped[il][ir][im][ic]==0)continue;
                  if(sig[il][ir][im][ic]==0)continue;
                  if(disconnected[il][ir][im][ic])continue;
                  if(!data_calib->GetTrkCalib().IsGood(il,ir,im,ic))continue;
                  hpulser->Fill(ir*NCHANNELS+ic,il*NMODULES+im);//1
                  if(ic_best<0)continue;
                  hpulser->Fill(ir*NCHANNELS+ic,il*NMODULES+im);//1+1
                }
                if(ic_best>=0){
                  hpulser->Fill(ir*NCHANNELS+ic_best,il*NMODULES+im);//1+1+1
                  //                  std::cout<<std::endl<<il<<ir<<im<<" - "<<ih<<" - channel best "<<ic_best;
                }
                if(ic_min>=0)chmin[il][ir][im].push_back(ic_min);

                ///////////////////////////////
                for(int i=0; i<N; i++){
                  uint ic = ih*N+i;
                  uint ch = il*NROWS*NMODULES*NCHANNELS + ir*NMODULES*NCHANNELS + im*NCHANNELS + ic;
                  double sig3_min = TMath::Sqrt( sigma_n[i]*sigma_n[i]+sigma_n[i_min]*sigma_n[i_min]*cgain[i]*cgain[i]/cgain[i_min]/cgain[i_min]);
                  double sig3_best = TMath::Sqrt(sigma_n[i]*sigma_n[i]+sigma_c*sigma_c*cgain[i]*cgain[i]);
                  if(ic_best>0)sig3_best = sig3_min;
                  gsig3_best->AddPoint(ch,sig3_best);
                  gsig3_min->AddPoint(ch,sig3_min);
                }

              }

            }//modules
          }//rows
        }//layers
      }//end subcn condition
      //////////////////////////////////////////////////////



    }//endl last iteration condition

  }//end iterations





  for(int il=0; il<NLAYERS; il++){
    for(int ir=0; ir<NROWS; ir++){
      for(int im=0; im<NMODULES; im++){
  	for(int ih=0; ih<NCN; ih++){
  	  if( !hcn[il][ir][im][ih]  )continue;
	  //	  if(hcn[il][ir][im][ih]->GetRMS()<1 )continue;
          //	  cout << endl << il*100+ir*10+im<<" RMS(baseline) "<<sigcn[il][ir][im][ih];
	  hsigcn->Fill(ir*NCN+ih,il*NMODULES+im,sigcn[il][ir][im][ih]);
	  hsigcndist->Fill(sigcn[il][ir][im][ih]);
  	}
      }
    }
  }




  for(int il = 0 ; il<NLAYERS; il++){
    for(int ir = 0 ; ir<NROWS; ir++){
      for(int im = 0 ; im<NMODULES; im++){
	for(int ic = 0 ; ic<NCHANNELS; ic++){
	  if( n[il][ir][im][ic] == 0  )continue;
	  hped->Fill(ir*NCHANNELS+ic,il*NMODULES+im,ped[il][ir][im][ic]);
	  hpeddist->Fill(ped[il][ir][im][ic]);
	  hsig->Fill(ir*NCHANNELS+ic,il*NMODULES+im,sig[il][ir][im][ic]);
	  hsigdist->Fill(sig[il][ir][im][ic]);
          uint ch = il*NROWS*NMODULES*NCHANNELS + ir*NMODULES*NCHANNELS + im*NCHANNELS + ic;
          gsig0->AddPoint(ch,sig0[il][ir][im][ic]);
          gsig->AddPoint(ch,sig[il][ir][im][ic]);
	}

      }
    }
  }

  TString oofile =flist+"-pedestals";
  if(subcn)oofile += Form("-cn%i-mod%i",NCN,cnmod);
  oofile += ".txt";
  fstream fs (oofile.Data(), std::fstream::out);
  fs << "Layer	Row	Module	Channel	Pedestal";
  for(int il = 0 ; il<NLAYERS; il++){
    for(int ir = 0 ; ir<NROWS; ir++){
      for(int im = 0 ; im<NMODULES; im++){
	for(int ic = 0 ; ic<NCHANNELS; ic++){
	  if( n[il][ir][im][ic] == 0  )continue;
	  fs << endl;
	  fs << setw(10) <<il;
	  fs << setw(10) <<ir;
	  fs << setw(10) <<im;
	  fs << setw(10) <<ic;
	  fs << setw(10) <<ped[il][ir][im][ic];
	  fs << setw(10) <<sig[il][ir][im][ic];
	}
      }
    }
  }
  fs.close();
  cout << endl <<"PED written ti file >> "<< oofile;

  if(subcn){
    TString oofile =flist+"-impulse";
    if(subcn)oofile += Form("-cn%i-mod%i",NCN,cnmod);
    oofile += ".txt";
    fstream fs (oofile.Data(), std::fstream::out);
    fs << "L	R	M	C";
    for(int il = 0 ; il<NLAYERS; il++){
      for(int ir = 0 ; ir<NROWS; ir++){
        for(int im = 0 ; im<NMODULES; im++){
          fs << endl;
          fs << setw(4) <<il;
          fs << setw(4) <<ir;
          fs << setw(4) <<im;
          for(auto ic:chmin[il][ir][im])fs << setw(4) << ic;
        }
      }
    }
    fs.close();
    cout << endl <<"Channel to impulse written ti file >> "<< oofile;
  }

  TString oofileroot = oofile;
  oofileroot += ".root";

  TFile* fi = new TFile(oofileroot.Data(),"recreate");
  hped->Write();
  hpeddist->Write();
  hsig->Write();
  hsigdist->Write();
  hsigcn->Write();
  hsigcndist->Write();
  for(auto o :hsok)o->Write();
  for(auto o :hsnotok)o->Write();
  for(auto o :vhrho)o->Write();
  for(auto o :vhsig)o->Write();
  for(auto o :vhrho_ex)o->Write();
  for(auto o :vhsig_ex)o->Write();

  for(int il=0; il<NLAYERS; il++){
    for(int ir=0; ir<NROWS; ir++){
      for(int im=0; im<NMODULES; im++){
  	for(int ih=0; ih<NCN; ih++){
  	  if(hcn[il][ir][im][ih]) hcn[il][ir][im][ih]->Write();
  	}
  	for(int ic=0; ic<NCHANNELS; ic++){
  	  if(hs[il][ir][im][ic]) hs[il][ir][im][ic]->Write();
	}
      }
    }
  }
  gsig0->Write();
  gsig->Write();
  hfull->Write();

  fi->Close();
  cout << endl << "Histo written to file "<<oofileroot;

  vector<TLine*> lines;
  for (uint i=0; i<NROWS; i++)  lines.push_back( new TLine((i+1)*NCHANNELS,            0,(i+1)*NCHANNELS,NLAYERS*NMODULES));
  for (uint i=0; i<NLAYERS; i++)lines.push_back( new TLine(             0,(i+1)*NMODULES, NROWS*NCHANNELS, (i+1)*NMODULES));

  TString oofilepdf =oofile;
  oofilepdf += ".pdf";

  vector<TCanvas*> vca;

  //==============================================================
  TCanvas *ccr3 = new TCanvas("ccr3","Trk Raw",900,600);
  ccr3->Divide(2,1);

  /////////////////////////////////////////////// figura 1
  ccr3->GetPad(1)->cd();

  if(hped)hped->Draw("colz");

  for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");
  /////////////////////////////////////////////// figura 2
  ccr3->GetPad(2)->cd();

  if(hpeddist)hpeddist->Draw("");


  vca.push_back(ccr3);


  //==============================================================
  TCanvas *ccr4 = new TCanvas("ccr4","Trk Raw",900,600);
  ccr4->Divide(2,1);

  /////////////////////////////////////////////// figura 1
  ccr4->GetPad(1)->cd();

  if(hsig)hsig->Draw("colz");

  for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");
  /////////////////////////////////////////////// figura 2
  ccr4->GetPad(2)->cd();

  if(hsigdist)hsigdist->Draw("");


  vca.push_back(ccr4);

  //==============================================================

  TCanvas *ccr4p5 = new TCanvas("ccr4p5","Trk Raw",900,600);
  ccr4p5->Divide(2,1);

  /////////////////////////////////////////////// figura 1
  ccr4p5->GetPad(1)->cd();

  if(hsig_0)hsig_0->Draw("colz");

  for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");
  /////////////////////////////////////////////// figura 2
  ccr4p5->GetPad(2)->cd();

  if(hsig_0dist)hsig_0dist->Draw("");


  vca.push_back(ccr4p5);

  cout << endl ;

  //===============================================================

  // TCanvas *cc = new TCanvas("cc","",900,600);
  // cc->Divide(3,1);
  // for(int ih=0; ih<hsok.size(); ih++){

  //   cc->GetPad(ih+1)->cd();
  //   cc->GetPad(ih+1)->SetLogy();
  //   hsok[ih]->Draw();
  //   hsnotok[ih]->SetLineColor(kRed);
  //   hsnotok[ih]->Draw("same");

  // }

  // return;
  // vca.push_back(cc);


 //==============================================================
  TCanvas *ccr5 = new TCanvas("ccr5","Trk Raw",900,600);
  ccr5->Divide(2,1);

  /////////////////////////////////////////////// figura 1
  ccr5->GetPad(1)->cd();

  if(hgain)hgain->Draw("colz");

  for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");
  /////////////////////////////////////////////// figura 2
  ccr5->GetPad(2)->cd();
  ccr5->GetPad(2)->SetLogy();

  if(hgaindist)hgaindist->Draw("");

  vca.push_back(ccr5);



  //==============================================================
  TCanvas *cccn = new TCanvas("cccn","Trk Raw",900,600);
  cccn->Divide(2,1);

  /////////////////////////////////////////////// figura 1
  cccn->GetPad(1)->cd();

  if(hsigcn)hsigcn->Draw("colz");

  //  for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");
  /////////////////////////////////////////////// figura 2
  cccn->GetPad(2)->cd();

  if(hsigcndist)hsigcndist->Draw("");

  vca.push_back(cccn);


  //==============================================================
  TCanvas *ccr6 = new TCanvas("ccr6","Trk Raw",900,600);
  ccr6->Divide(2,1);

  /////////////////////////////////////////////// figura 1
  ccr6->GetPad(1)->cd();

  if(hpulser)hpulser->Draw("colz");

  for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");
  /////////////////////////////////////////////// figura 2
  ccr6->GetPad(2)->cd();



  vca.push_back(ccr6);




  //==============================================================
  TCanvas *css = new TCanvas("css","",900,600);
  css->SetLogy();

  TH1F* hss = new TH1F("hss","",10,0,NLAYERS*NROWS*NMODULES*NCHANNELS);
  hss->SetStats(0);
  hss->SetMaximum(100);
  hss->SetMinimum(0.1);
  hss->Draw();



  gsig0->SetMarkerStyle(7);
  gsig0->SetMarkerColor(kBlack);
  gsig0->SetLineColor(kBlack);
  gsig0->Draw("p");

  gsig3_min->SetMarkerStyle(7);
  gsig3_min->SetMarkerColor(kBlue);
  gsig3_min->SetLineColor(kBlue);
  gsig3_min->Draw("p");

  gsig3_best->SetMarkerStyle(7);
  gsig3_best->SetMarkerColor(kCyan);
  gsig3_best->SetLineColor(kCyan);
  gsig3_best->Draw("p");


  gsig->SetMarkerStyle(7);
  gsig->SetMarkerColor(kRed);
  gsig->SetLineColor(kRed);
  gsig->Draw("p");

  TLine *ldis = new TLine(gsig0->GetPointX(0),sigmin,gsig0->GetPointX(gsig0->GetN()-1),sigmin);
  ldis->SetLineStyle(2);
  ldis->Draw("same");

  for(int il = 0 ; il<NLAYERS; il++){
    uint ch = il*NROWS*NMODULES*NCHANNELS ;
    TLine *l = new TLine((float)ch,0,(float)ch,100);
    l->Draw("same");
  }



  vca.push_back(css);

  TCanvas *csss = new TCanvas("csss","",900,600);
  //  csss->SetLogy();
  hfull->Draw("colz");
  vca.push_back(csss);

  for(int ic=0;ic<vca.size(); ic++){
    if     (ic==0)             {
      vca[ic]->Print(Form("%s[",oofilepdf.Data()));
      vca[ic]->Print(Form("%s",oofilepdf.Data()));
    }else if(ic==vca.size()-1 ){
      vca[ic]->Print(Form("%s",oofilepdf.Data()));
      vca[ic]->Print(Form("%s]",oofilepdf.Data()));
    }else
      vca[ic]->Print(Form("%s",oofilepdf.Data()));
  }


  gROOT->SetBatch(kTRUE);

  vca.clear();
  /////////////////////////////////////////////////////////////////////
  for(int il = 0 ; il<NLAYERS; il++){
    for(int ir = 0 ; ir<NROWS; ir++){

      ////////////////////////////////////////////////////////////
      TCanvas *c = new TCanvas(Form("c_%i%i",il,ir),"",900,600);
      c->Divide(3,2);
      for(int im = 0 ; im<NMODULES; im++){
        c->GetPad(im+1)->cd();
        c->GetPad(im+1)->SetLogy();

        if( hsall0[il][ir][im] )hsall0[il][ir][im]->SetLineColor(38);
        if( hsall0[il][ir][im] )hsall0[il][ir][im]->Draw();
        //if( hsmip[il][ir][im] )hsmip[il][ir][im]->SetLineColor(kRed);
        //if( hsmip[il][ir][im] )hsmip[il][ir][im]->Draw("same");
        if( hsped[il][ir][im] )hsped[il][ir][im]->SetLineColor(kRed);
        if( hsped[il][ir][im] )hsped[il][ir][im]->Draw("same");
        else cout << endl <<"L"<<il<<" R"<<ir<<" M"<<im<<" empty  ";
      }
      c->Update();
      vca.push_back( c );
      ////////////////////////////////////////////////////////////

    }
  }
  TString oo = oofilepdf+"-dist.pdf";
  for(int ic=0;ic<vca.size(); ic++){
    if     (ic==0)             {
      vca[ic]->Print(Form("%s[",oo.Data()));
      vca[ic]->Print(Form("%s",oo.Data()));
    }else if(ic==vca.size()-1 ){
      vca[ic]->Print(Form("%s",oo.Data()));
      vca[ic]->Print(Form("%s]",oo.Data()));
    }else
      vca[ic]->Print(Form("%s",oo.Data()));
  }


    gROOT->SetBatch(kFALSE);

}

void ShowRho(TFile* f,int detector){
  if(!f)return;
  TH2F* hrho = (TH2F*)f->Get(Form("hrho_%i",detector));
  if(!hrho)return;
  TH1F* hsig = (TH1F*)f->Get(Form("hsig_%i",detector));
  if(!hsig)return;
  hrho->SetStats(0);
  hsig->SetStats(0);

  TCanvas *ccr3 = new TCanvas("ccr3","Trk Raw",900,600);
  ccr3->Divide(2,1);
  ccr3->cd(1);
  hrho->GetZaxis()->SetRangeUser(-1,1);
  hrho->Draw("colz");
  ccr3->cd(2);
  ccr3->GetPad(2)->SetLogy();
  hsig->SetMinimum(0.1);
  hsig->SetMaximum(500);
  hsig->Draw("hist");

}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void ShowBaseline(TFile* file, int mo, int ncn){
  if(!file)return;

  TTree *tree = (TTree*)file->Get("tree");
  if(!tree)return;

  int layer = (int)(mo/100);
  int row = (int)((mo-layer*100)/10);
  int module = (int)((mo-layer*100-row*10));


  TH2F *h[ncn];
  TProfile *hp[ncn];

  TH1F *hs[ncn];
  TH1F *hb[ncn];

  for(int id=0; id<ncn; id++){


    TString hsname = Form("hs%i",id);
    TH1F* hhs = (TH1F*)gDirectory->FindObject(hsname.Data());
    if(hhs)hhs->Delete();
    tree->Draw(
	       Form("signal[%i][%i][%i][%i]>>%s(500,-100,1000)",layer,row,module,id,hsname.Data()),
	       "",
	       "");

    TString hbname = Form("hb%i",id);
    TH1F* hhb = (TH1F*)gDirectory->FindObject(hbname.Data());
    if(hhb)hhb->Delete();
    tree->Draw(
	       Form("base[%i][%i][%i][%i]>>%s(500,-100,1000)",layer,row,module,id,hbname.Data()),
	       "",
	       "");

    hhs = (TH1F*)gDirectory->FindObject(hsname.Data());
    hhb = (TH1F*)gDirectory->FindObject(hbname.Data());
    hhs->SetLineColor(kRed);

    hs[id]=hhs;
    hb[id]=hhb;

    // ------------------------------------------------
    TCut cut = Form("signal[%i][%i][%i][%i]>200",layer,row,module,id);
    cut = cut && Form("nbase[%i][%i][%i][]>=7",layer,row,module);
    // ------------------------------------------------

    TString hname = Form("h%i",id);
    TH2F* hh = (TH2F*)gDirectory->FindObject(hname.Data());
    if(hh)hh->Delete();
    tree->Draw(
	       Form("base[%i][%i][%i][]:Iteration$>>%s",layer,row,module,hname.Data()),
	       cut,
	       "box");
    //
    TString hpname = Form("hp%i",id);
    TProfile* hhp = (TProfile*)gDirectory->FindObject(hpname.Data());
    if(hhp)hhp->Delete();
    tree->Draw(
	       Form("base[%i][%i][%i][]:Iteration$>>%s",layer,row,module,hpname.Data()),
	       cut,
	       "prof");


    hh = (TH2F*)gDirectory->FindObject(hname.Data());
    hhp = (TProfile*)gDirectory->FindObject(hpname.Data());
    hhp->SetLineColor(kRed);
    hhp->SetLineWidth(2);
    h[id]=hh;
    hp[id]=hhp;

  }

  TCanvas *ccc = new TCanvas("ccc","ccc",900,300);
  ccc->Divide(ncn,1);

  for(int id=0; id<ncn; id++){

    ccc->GetPad(id+1)->SetLogy();
    ccc->cd(id+1);
    hs[id]->Draw();
    hb[id]->Draw("same");

  }
  TCanvas *cc = new TCanvas("cc","cc",900,300);
  cc->Divide(ncn,1);

  for(int id=0; id<ncn; id++){

    cc->cd(id+1);
    h[id]->Draw();
    hp[id]->Draw("same");

  }

}

TTree* EvaluateBaseline(TString flist, TString fped , long nev, TString ddir,TString suffix){

  TChain *craw = GetChain(flist.Data(),ddir,"TreeRaw",suffix);
  TChain *ccal = GetChain(flist.Data(),ddir,"TreeRec",suffix);
  craw->AddFriend(ccal);

  //----------------------------------------- reconstructed event
  CEventRec *evrec= new CEventRec();
  craw->SetBranchAddress("Rec", &evrec);
  //  evrec->ChooseReconstruction("FindHough3D");
  //----------------------------------------- raw event
  Crane::Calibration::CRawTrk *rawtrk = new Crane::Calibration::CRawTrk();
  craw->SetBranchAddress("Trk", &rawtrk);
  //-----------------------------------------

  map<uint,double> pedmap = GetPedestals(fped.Data());
  map<uint,double> sigmap = GetPedestals(fped.Data(),true);
  cout << endl << " ped map size "<<pedmap.size();
  double sigmin = 1.7;//cut on sig to tag disconnected channels
  double sigmax = 100;//8;//cut on sig to tag disconnected channels

  bool good[NLAYERS][NROWS][NMODULES][NCHANNELS];
  fill_n( &good[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);
  int conta=0;
  for (std::map<uint,double>::iterator it=sigmap.begin(); it!=sigmap.end(); ++it){
    int key    = it->first;
    int il = (int)key/10000;
    int ir = (int)(key-il*10000)/1000;
    int im = (int)(key-il*10000-ir*1000)/100;
    int ic = (int)(key-il*10000-ir*1000-im*100);
    double sig = it->second;
    if( sig > sigmin && sig < sigmax )good[il][ir][im][ic]=1;
    else conta++;
  }
  cout << endl << "Disconnected channels "<<100*conta/sigmap.size()<<"%";

  int    module_n[NLAYERS][NROWS][NMODULES];
  float  adc[NLAYERS][NROWS][NMODULES][NCHANNELS];
  float  adc_corr[NLAYERS][NROWS][NMODULES][NCHANNELS];

  float  su[NLAYERS][NROWS][NMODULES][NCN];
  float  vu[NLAYERS][NROWS][NMODULES][NCN];
  int    nsu[NLAYERS][NROWS][NMODULES][NCN];
  float  min[NLAYERS][NROWS][NMODULES][NCN];
  float  smin[NLAYERS][NROWS][NMODULES][NCN];

  float b[NLAYERS][NROWS][NMODULES][NCN];
  float s[NLAYERS][NROWS][NMODULES][NCN];


  cout << endl << "Creating file "<<Form("%s-baseline.root",flist.Data());
  TFile* fi = new TFile(Form("%s-baseline.root",flist.Data()),"recreate");
  TTree *tree = new TTree("tree","tree");
  tree->Branch("adc",adc,Form("adc[%i][%i][%i][%i]/F",NLAYERS,NROWS,NMODULES,NCHANNELS));
  tree->Branch("adc_corr",adc_corr,Form("adc_corr[%i][%i][%i][%i]/F",NLAYERS,NROWS,NMODULES,NCHANNELS));

  tree->Branch("base",b,Form("base[%i][%i][%i][%i]/F",NLAYERS,NROWS,NMODULES,NCN));
  tree->Branch("signal",s,Form("signal[%i][%i][%i][%i]/F",NLAYERS,NROWS,NMODULES,NCN));
  tree->Branch("vbase",vu,Form("vbase[%i][%i][%i][%i]/F",NLAYERS,NROWS,NMODULES,NCN));
  tree->Branch("nbase",nsu,Form("nbase[%i][%i][%i][%i]/I",NLAYERS,NROWS,NMODULES,NCN));

  // double  sum[NLAYERS][NROWS][NMODULES][NCHANNELS];
  // double sum2[NLAYERS][NROWS][NMODULES][NCHANNELS];
  // int       n[NLAYERS][NROWS][NMODULES][NCHANNELS];
  // fill_n( &sum[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);
  // fill_n(&sum2[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);
  // fill_n(   &n[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);

  // double cov[NLAYERS][NROWS][NMODULES][NCHANNELS][NCHANNELS];
  // fill_n( &cov[0][0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS*NCHANNELS,0);

  // double  ave[NLAYERS][NROWS][NMODULES][NCHANNELS];
  // double  rms[NLAYERS][NROWS][NMODULES][NCHANNELS];
  // fill_n( &ave[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);
  // fill_n( &rms[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);




  cout << endl << " N.entries "<<craw->GetEntries();

  for (uint ie=0; ie< TMath::Min((long)nev,(long)craw->GetEntries()); ie++){
  //  for (uint ie=0; ie<100; ie++){

    craw->GetEntry(ie);

    fill_n(&module_n[0][0][0],NLAYERS*NROWS*NMODULES,0);
    fill_n( &adc[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);

    fill_n(&min[0][0][0][0],NLAYERS*NROWS*NMODULES*NCN,99999999);
    fill_n(&smin[0][0][0][0],NLAYERS*NROWS*NMODULES*NCN,0);


    // fill_n( &cov[0][0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS*NCHANNELS,0);
    // fill_n( &sum[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);
    // fill_n(&sum2[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);
    // fill_n(   &n[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);


    // ----------------------------
    // conta le hit per ogni modulo
    // ----------------------------
    for(uint i=0; i<(uint)rawtrk->adcdata.size(); i++){
      int layer  = rawtrk->layer[i];
      int row    = rawtrk->row[i];
      int module = rawtrk->module[i];
      int channel = rawtrk->channel[i];
      module_n[layer][row][module]++;
    }
    // ----------------------------
    // riempie la matrice adc, sottraendo il piedistallo
    // cerca il minimo
    // ----------------------------
    for(uint i=0; i<(uint)rawtrk->adcdata.size(); i++){
      int layer  = rawtrk->layer[i];
      int row    = rawtrk->row[i];
      int module = rawtrk->module[i];
      int channel = rawtrk->channel[i];

      int half = (channel < (int)NCHANNELS/2 ? 0 : 1 );
      int detector = (int)(channel/8);
      int ih = 0;//(NCN==2 ? half : detector);
      if(NCN==2) ih = half;
      if(NCN==4) ih = detector;


      bool FULL = module_n[layer][row][module] == NCHANNELS;
      if(!FULL)continue;

      uint key = layer*10000 + row*1000 + module*100 + channel ;
      std::map<uint,double>::iterator itp = pedmap.find(key);
      std::map<uint,double>::iterator its = sigmap.find(key);

      double val = (double)(rawtrk->adcdata[i]);
      if( itp == pedmap.end() || its == sigmap.end() ){
	cout << endl << "MISSING PED-SIG - skip ";
	continue;
      }
      adc[layer][row][module][channel] = val-itp->second;

      if(  adc[layer][row][module][channel] < min[layer][row][module][ih]){
	min[layer][row][module][ih]  = adc[layer][row][module][channel];
	smin[layer][row][module][ih] = its->second;
      }
    }
    // -----------------------------------
    // calcola la baseline, iterativamente
    // -----------------------------------

    fill_n( &adc_corr[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);
    fill_n( &b[0][0][0][0],NLAYERS*NROWS*NMODULES*NCN,0);//baseline

    int nit = 3;
    double sigcut = 10;//4;
    double mipcut = 200;//cut on signal to identify mips
    int nsumin = 4;
    for(int it=0; it<nit;it++){ //iterations


      fill_n(&su[0][0][0][0],NLAYERS*NROWS*NMODULES*NCN,0);
      fill_n(&nsu[0][0][0][0],NLAYERS*NROWS*NMODULES*NCN,0);
      fill_n(&vu[0][0][0][0],NLAYERS*NROWS*NMODULES*NCN,0);
      fill_n(&s[0][0][0][0],NLAYERS*NROWS*NMODULES*NCN,0);
      for(int il=0; il<NLAYERS; il++){
	for(int ir=0; ir<NROWS; ir++){
	  for(int im=0; im<NMODULES; im++){

	    bool FULL = module_n[il][ir][im] == NCHANNELS;
	    if(!FULL)continue;

	    // --------------------------
	    // increment counters
	    // --------------------------
	    for(int ic = 0 ; ic<NCHANNELS; ic++){

	      if(!good[il][ir][im][ic])continue;//skip disconnected channels

              int half = (ic < (int)NCHANNELS/2 ? 0 : 1 );
              int detector = (int)(ic/8);
              int ih = 0;//(NCN==2 ? half : detector);
              if(NCN==2) ih = half;
              if(NCN==4) ih = detector;

	      uint key = il*10000 + ir*1000 + im*100 + ic ;
	      std::map<uint,double>::iterator its = sigmap.find(key);
	      if(its == sigmap.end())continue;
	      double sig = its->second;

	      bool ISHIT=false;
	      if(it==0){//first iteration
		sig = TMath::Sqrt(sig*sig+smin[il][ir][im][ih]*smin[il][ir][im][ih]);
		if( adc[il][ir][im][ic] - min[il][ir][im][ih] > sigcut * sig)ISHIT=true;
		if( adc[il][ir][im][ic] - min[il][ir][im][ih]   > mipcut)ISHIT=true;
	      }else{
		if( adc[il][ir][im][ic] - b[il][ir][im][ih]   > sigcut * sig)ISHIT=true;
		if( adc[il][ir][im][ic] - b[il][ir][im][ih]   > mipcut)ISHIT=true;
	      }
	      // segnale
	      if(ISHIT){
		s[il][ir][im][ih]+=adc[il][ir][im][ic];
		continue;
	      }
	      // --------
	      // baseline
	      // --------
	      nsu[il][ir][im][ih]++;
	      su[il][ir][im][ih] += adc[il][ir][im][ic] ;
	      vu[il][ir][im][ih] += adc[il][ir][im][ic]*adc[il][ir][im][ic] ;
	    }//end loop over channels

	    // --------------------------
	    // baseline evaluation
	    // --------------------------
	    for(int ih=0; ih<NCN; ih++){
	      b[il][ir][im][ih] = 0;
	      if( nsu[il][ir][im][ih] >= nsumin ) {
		b[il][ir][im][ih]=su[il][ir][im][ih]/nsu[il][ir][im][ih];
		vu[il][ir][im][ih] = vu[il][ir][im][ih]/nsu[il][ir][im][ih];
		vu[il][ir][im][ih] -= b[il][ir][im][ih]*b[il][ir][im][ih];
		vu[il][ir][im][ih] = TMath::Sqrt( vu[il][ir][im][ih] );
		//		cout << endl << b[il][ir][im][ih] << " "<<vu[il][ir][im][ih];
	      }
	    }

	    // --------------------------
	    // baseline subtraction
	    // --------------------------

	    for(int ic = 0 ; ic<NCHANNELS; ic++){
	      int ih = (int)(ic/NCN);
	      adc_corr[il][ir][im][ic] = adc[il][ir][im][ic] - b[il][ir][im][ih];
	    }

	  }//modules
	}//rows
      }//layers
    }//end iterations

    tree->Fill();

  }//endl loop over events



  fi->cd();
  tree->Write();
  // for(auto o :vhrho)o->Write();

  fi->Close();
  return tree;

}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void CheckMaskFile(std::string file, double max){

  std::fstream filestr (file, std::fstream::in );
  if(!filestr.is_open()) {
    cout << endl <<" File "<<file<<" not found ";
    return;
  }
  CTrkCalib calib;
  calib.SetMasks(file);


  TH2F *hmask = (TH2F*)gROOT->FindObject("hmask");
  if(hmask)delete hmask;
  hmask = new TH2F("hmask","",NCHANNELS*NROWS,0,NCHANNELS*NROWS,NLAYERS*NMODULES,0,NLAYERS*NMODULES);
  hmask->GetXaxis()->SetTitle(Form("row*%i+channel",NCHANNELS));
  hmask->GetYaxis()->SetTitle(Form("layer*%i+module",NMODULES));
  hmask->SetStats(0);

  TH1F* hmaskdist = new TH1F("hmaskdist","PED",2,0,2);


  // map<uint,double> pedmap = GetPedestals(file);
  map<uint,bool> maskmap = calib.GetMaskMap();

  cout << endl << " map size "<<maskmap.size();


  map<uint,pair<int,int>> mcmap = CTrkCalib::GetModuleChannelMap();


  for (std::map<uint,bool>::iterator it=maskmap.begin(); it!=maskmap.end(); ++it){
    int key    = it->first;
    int val = (int)(it->second);
    // int il = (int)key/10000;
    // int ir = (int)(key-il*10000)/1000;
    // int im = (int)(key-il*10000-ir*1000)/100;
    // int ic = (int)(key-il*10000-ir*1000-im*100);
    // int il = GGeometryObject::GetTrackerLayer(key);
    // int ir = GGeometryObject::GetLayerRow(key);
    // int im = GGeometryObject::GetRowModule(key);
    // int ic = GGeometryObject::GetModuleChannel(key);

    auto itt = mcmap.find(key);
    if( itt != mcmap.end() ){
      int il =  itt->second.first/100;//layer
      int ir = (itt->second.first-il*100)/10;//row
      int im = (itt->second.first-il*100-ir*10);//module
      int ic =  itt->second.second; //channel

      int valcheck = calib.IsGood(il,ir,im,ic);

      hmask->Fill(ir*NCHANNELS+ic,il*NMODULES+im,valcheck);
      hmaskdist->Fill(valcheck);

    }

  }

  vector<TLine*> lines;
  for (uint i=0; i<NROWS; i++)  lines.push_back( new TLine((i+1)*NCHANNELS,            0,(i+1)*NCHANNELS,NLAYERS*NMODULES));
  for (uint i=0; i<NLAYERS; i++)lines.push_back( new TLine(             0,(i+1)*NMODULES, NROWS*NCHANNELS, (i+1)*NMODULES));


  TCanvas *ccr3 = new TCanvas("ccr3","Trk Raw",900,600);
  ccr3->Divide(2,1);

  /////////////////////////////////////////////// figura 1
  ccr3->GetPad(1)->cd();

  if(hmask)hmask->Draw("colz");

  for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");
  /////////////////////////////////////////////// figura 2
  ccr3->GetPad(2)->cd();

  if(hmaskdist)hmaskdist->DrawNormalized("");


  DrawTitle(ccr3->GetPad(1),Form("%s",file.data()));


  return;

}

void CheckPedestalFile(std::string file, double max){

  std::fstream filestr (file, std::fstream::in );
  if(!filestr.is_open()) {
    cout << endl <<" Tracker pedestal file "<<file<<" not found ";
    return;
  }

  TH2F *hped = (TH2F*)gROOT->FindObject("hped");
  if(hped)delete hped;
  hped = new TH2F("hped","",NCHANNELS*NROWS,0,NCHANNELS*NROWS,NLAYERS*NMODULES,0,NLAYERS*NMODULES);
  hped->GetXaxis()->SetTitle(Form("row*%i+channel",NCHANNELS));
  hped->GetYaxis()->SetTitle(Form("layer*%i+module",NMODULES));
  hped->SetStats(0);

  TH1F* hpeddist = new TH1F("hpeddist","PED",300,0,max);


  map<uint,double> pedmap = GetPedestals(file);
  cout << endl << " map size "<<pedmap.size();


  for (std::map<uint,double>::iterator it=pedmap.begin(); it!=pedmap.end(); ++it){
    int key    = it->first;
    double ped = it->second;
    if(ped==0)cout << endl << "skip PED="<<ped<<" "<<it->first;
    if(ped==0)continue;
    int il = (int)key/10000;
    int ir = (int)(key-il*10000)/1000;
    int im = (int)(key-il*10000-ir*1000)/100;
    int ic = (int)(key-il*10000-ir*1000-im*100);
    hped->Fill(ir*NCHANNELS+ic,il*NMODULES+im,ped);
    hpeddist->Fill(ped);
  }

  vector<TLine*> lines;
  for (uint i=0; i<NROWS; i++)  lines.push_back( new TLine((i+1)*NCHANNELS,            0,(i+1)*NCHANNELS,NLAYERS*NMODULES));
  for (uint i=0; i<NLAYERS; i++)lines.push_back( new TLine(             0,(i+1)*NMODULES, NROWS*NCHANNELS, (i+1)*NMODULES));


  TCanvas *ccr3 = new TCanvas("ccr3","Trk Raw",900,600);
  ccr3->Divide(2,1);

  /////////////////////////////////////////////// figura 1
  ccr3->GetPad(1)->cd();

  if(hped)hped->Draw("colz");

  for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");
  /////////////////////////////////////////////// figura 2
  ccr3->GetPad(2)->cd();

  if(hpeddist)hpeddist->Draw("");


  DrawTitle(ccr3->GetPad(1),Form("%s",file.data()));


  return;

}

void CheckSigmaFile(std::string file, double max){


  std::fstream filestr (file, std::fstream::in );
  if(!filestr.is_open()) {
    cout << endl <<" Tracker pedestal file "<<file<<" not found ";
    return;
  }

  TH2F *hsig = (TH2F*)gROOT->FindObject("hsig");
  if(hsig)delete hsig;
  hsig = new TH2F("hsig","",NCHANNELS*NROWS,0,NCHANNELS*NROWS,NLAYERS*NMODULES,0,NLAYERS*NMODULES);
  hsig->GetXaxis()->SetTitle(Form("row*%i+channel",NCHANNELS));
  hsig->GetYaxis()->SetTitle(Form("layer*%i+module",NMODULES));
  hsig->SetStats(0);

  TH1F* hsigdist = new TH1F("hsigdist","SIG",300,0,max);


  map<uint,double> pedmap = GetPedestals(file,true);//read sigmas instead
  cout << endl << " map size "<<pedmap.size();


  for (std::map<uint,double>::iterator it=pedmap.begin(); it!=pedmap.end(); ++it){
    int key    = it->first;
    double ped = it->second;
    int il = (int)key/10000;
    int ir = (int)(key-il*10000)/1000;
    int im = (int)(key-il*10000-ir*1000)/100;
    int ic = (int)(key-il*10000-ir*1000-im*100);
    hsig->Fill(ir*NCHANNELS+ic,il*NMODULES+im,ped);
    hsigdist->Fill(ped);

  }

  vector<TLine*> lines;
  for (uint i=0; i<NROWS; i++)  lines.push_back( new TLine((i+1)*NCHANNELS,            0,(i+1)*NCHANNELS,NLAYERS*NMODULES));
  for (uint i=0; i<NLAYERS; i++)lines.push_back( new TLine(             0,(i+1)*NMODULES, NROWS*NCHANNELS, (i+1)*NMODULES));


  TCanvas *ccr3 = new TCanvas("ccr3","Trk Raw",900,600);
  ccr3->Divide(2,1);

  /////////////////////////////////////////////// figura 1
  ccr3->GetPad(1)->cd();

  hsig->GetZaxis()->SetRangeUser(0,max);
  if(hsig)hsig->Draw("colz");

  for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");
  /////////////////////////////////////////////// figura 2
  ccr3->GetPad(2)->cd();

  if(hsigdist)hsigdist->Draw("");


  DrawTitle(ccr3->GetPad(1),Form("%s",file.data()));


  return;


}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
map<uint,double> GetPedestals(std::string file, bool getsigmas ){

  map<uint,double> pedmap;

  //  std::string file = std::getenv("GAPS") + std::string("/resources/calibration/SiLi-pedestals.txt");
  cout << endl <<"Loading tracker pedestals from file: "<<file<<" ";
  std::fstream filestr (file, std::fstream::in );
  if(!filestr.is_open()) {
    cout << endl <<" Tracker pedestal file "<<file<<" not found ";
    return pedmap;
  }

  double ped = 0;
  int lineskip=1; //NB  !!! skip 1 line

  while(filestr.is_open()){
    char line[200];
    for(uint i=0 ; i<(uint)lineskip; i++)filestr.getline(line,200); //skip first lines
    for(;;){
      uint layer,row,module,channel;
      float val1=0;
      float val2=0;
      filestr >> layer;
      filestr >> row;
      filestr >> module;
      filestr >> channel;
      filestr >> val1;
      filestr >> val2;
      uint key = layer*10000 + row*1000 + module*100 + channel ;//GetVolumeId(layer,row,module,channel);

      //     if( !GoodChannel(layer,row,module,channel) )continue;
      if(!getsigmas)      pedmap.insert( pair<uint,double>(key,val1));
      else                pedmap.insert( pair<uint,double>(key,val2));
      if(!filestr.good())break;

    };
    filestr.close();
  }

  return pedmap;

}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void DiffFiles(std::string file1,std::string file2, bool getsigmas){

  map<uint,double> ped1 = GetPedestals(file1,getsigmas);
  cout << endl << "map size "<<ped1.size();

  map<uint,double> ped2 = GetPedestals(file2,getsigmas);
  cout << endl << "map size "<<ped2.size();

  TH2F *hdiff = (TH2F*)gROOT->FindObject("hdiff");
  if(hdiff)delete hdiff;
  hdiff = new TH2F("hdiff","",NCHANNELS*NROWS,0,NCHANNELS*NROWS,NLAYERS*NMODULES,0,NLAYERS*NMODULES);
  hdiff->GetXaxis()->SetTitle(Form("row*%i+channel",NCHANNELS));
  hdiff->GetYaxis()->SetTitle(Form("layer*%i+module",NMODULES));
  hdiff->SetStats(0);
  vector<TLine*> lines;
  for (uint i=0; i<NROWS; i++)  lines.push_back( new TLine((i+1)*NCHANNELS,            0,(i+1)*NCHANNELS,NLAYERS*NMODULES));
  for (uint i=0; i<NLAYERS; i++)lines.push_back( new TLine(              0,(i+1)*NMODULES, NROWS*NCHANNELS, (i+1)*NMODULES));

  TH1F* hdiffdist = new TH1F("hdiffdist","diff",300,-50,50);

  for (std::map<uint,double>::iterator it1=ped1.begin(); it1!=ped1.end(); ++it1){
    int key1  = it1->first;
    double v1 = it1->second;
    std::map<uint,double>::iterator it2 = ped2.find(key1);
    if (it2 != ped2.end()){
      double v2 = it2->second;
      int il = (int)key1/10000;
      int ir = (int)(key1-il*10000)/1000;
      int im = (int)(key1-il*10000-ir*1000)/100;
      int ic = (int)(key1-il*10000-ir*1000-im*100);
      hdiff->Fill(ir*NCHANNELS+ic,il*NMODULES+im,v1-v2);
      hdiffdist->Fill(v1-v2);
    }
  }

  TCanvas *ccr3 = new TCanvas("ccr3","Trk Raw",900,600);
  ccr3->Divide(2,1);

  /////////////////////////////////////////////// figura 1
  ccr3->GetPad(1)->cd();

  if(hdiff)hdiff->Draw("colz");
  for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");

  DrawTitle(ccr3->GetPad(1),file1.data(),0.07,0.95);
  DrawTitle(ccr3->GetPad(1),file2.data(),0.07,0.95-0.04);

  /////////////////////////////////////////////// figura 2
  ccr3->GetPad(2)->cd();

  if(hdiffdist)hdiffdist->Draw("");

};
///////////////////////////////////////////////////////////////////////////
void CheckTFFile(std::string file, int detector, std::string maskfile ){

  TH2F *htf = (TH2F*)gROOT->FindObject("htf");
  if(htf)delete htf;
  htf = new TH2F("htf","",NCHANNELS*NROWS,0,NCHANNELS*NROWS,NLAYERS*NMODULES,0,NLAYERS*NMODULES);
  htf->GetXaxis()->SetTitle(Form("row*%i+channel",NCHANNELS));
  htf->GetYaxis()->SetTitle(Form("layer*%i+module",NMODULES));
  htf->SetStats(0);

  TH1F *h001 = new TH1F("h001","0.01 MeV",100,0,1000);
  TH1F *h005 = new TH1F("h005","0.05 MeV",100,0,1000);
  TH1F *h01  = new TH1F("h01" ,"0.1 MeV",100,0,1000);
  TH1F *h1   = new TH1F("h1"  ,"1 MeV",100,0,1000);

  double mevmin = 1e-3;
  double mevmax = 100;
  double adcmax = 2047;

  TH1F* h2 = new TH1F("h2","",10,0.001,mevmax);
  h2->GetXaxis()->SetTitle("MeV");
  h2->GetYaxis()->SetTitle("S (ADC counts) ");
  h2->GetYaxis()->SetTitleOffset(1.5);
  h2->GetYaxis()->SetRangeUser(mevmin,mevmax);
  h2->SetStats(0);
  h2->SetMaximum(adcmax);

  TH2F* hgtf = new TH2F("hgtf","",500,-3,2,500,0,adcmax);


  vector<TGraph*> g;

  map<uint,TGraph*> gmap = GetTransferFunctions(file,detector,maskfile,true);


  for (std::map<uint,TGraph*>::iterator it=gmap.begin(); it!=gmap.end(); ++it){
    g.push_back( it->second );
    int key    = it->first;
    int il = (int)key/10000;
    int ir = (int)(key-il*10000)/1000;
    int im = (int)(key-il*10000-ir*1000)/100;
    int ic = (int)(key-il*10000-ir*1000-im*100);
    htf->Fill(ir*NCHANNELS+ic,il*NMODULES+im);

    cout << endl <<"get key "<<key;

    TGraph *gra = it->second;
    h001->Fill( gra->Eval(0.01) );
    h005->Fill( gra->Eval(0.05) );
    h01->Fill( gra->Eval(0.1) );
    h1->Fill( gra->Eval(1) );


    // for (int i = 1; i <= hgtf->GetYaxis()->GetNbins(); i++) {
    //   double y = hgtf->GetYaxis()->GetBinCenter(i);
    //   double x = TMath::Log10(gra->Eval(y));
    //   // cout<<endl<<x<<" "<<y;
    //   hgtf->Fill(x,y);
    // }

  }

  vector<TLine*> lines;
  for (uint i=0; i<NROWS; i++)  lines.push_back( new TLine((i+1)*NCHANNELS,            0,(i+1)*NCHANNELS,NLAYERS*NMODULES));
  for (uint i=0; i<NLAYERS; i++)lines.push_back( new TLine(             0,(i+1)*NMODULES, NROWS*NCHANNELS, (i+1)*NMODULES));


  TCanvas *ccr3 = new TCanvas("ccr3","Trk Raw",900,600);
  ccr3->Divide(2,1);

  /////////////////////////////////////////////// figura 1
  ccr3->GetPad(1)->cd();

  if(htf)htf->Draw("colz");

  for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");
  /////////////////////////////////////////////// figura 2
  ccr3->GetPad(2)->cd();
  ccr3->GetPad(2)->SetLogx();



  h2->Draw();
  for(auto gg: g)gg->Draw("pl");

  DrawTitle(ccr3->GetPad(1),Form("%s %i",file.data(),detector));

  TCanvas *cv = new TCanvas("cv","Trk Raw",900,600);
  cv->Divide(2,1);

  cv->GetPad(1)->cd();

  //  hgtf->Draw("colz");


  cv->GetPad(2)->cd();
  cv->GetPad(2)->SetLogy();

  h001->SetLineColor(kGreen);
  h001->Draw();

  h005->SetLineColor(kBlue);
  h005->Draw("same");

  h01->SetLineColor(kBlack);
  h01->Draw("same");

  h1->SetLineColor(kRed);
  h1->Draw("same");


}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void CheckGainFile(std::string file, double max){

  std::fstream filestr (file, std::fstream::in );
  if(!filestr.is_open()) {
    cout << endl <<" File "<<file<<" not found ";
    return;
  }
  CTrkCalib calib;
  calib.SetMasks(file);


  TH2F *hmask = (TH2F*)gROOT->FindObject("hmask");
  if(hmask)delete hmask;
  hmask = new TH2F("hmask","",NCHANNELS*NROWS,0,NCHANNELS*NROWS,NLAYERS*NMODULES,0,NLAYERS*NMODULES);
  hmask->GetXaxis()->SetTitle(Form("row*%i+channel",NCHANNELS));
  hmask->GetYaxis()->SetTitle(Form("layer*%i+module",NMODULES));
  hmask->SetStats(0);

  TH1F* hmaskdist = new TH1F("hmaskdist","Gain",300,0,max);


  // map<uint,double> pedmap = GetPedestals(file);
  //map<uint,double> maskmap = calib.GetGainMap();

  //cout << endl << " map size "<<maskmap.size();


  map<uint,pair<int,int>> mcmap = CTrkCalib::GetModuleChannelMap();

/*
  for (auto it=maskmap.begin(); it!=maskmap.end(); ++it){
    int key    = it->first;
    double val = (double)(it->second);
    // int il = (int)key/10000;
    // int ir = (int)(key-il*10000)/1000;
    // int im = (int)(key-il*10000-ir*1000)/100;
    // int ic = (int)(key-il*10000-ir*1000-im*100);
    // int il = GGeometryObject::GetTrackerLayer(key);
    // int ir = GGeometryObject::GetLayerRow(key);
    // int im = GGeometryObject::GetRowModule(key);
    // int ic = GGeometryObject::GetModuleChannel(key);
    cout << endl << val;
    auto itt = mcmap.find(key);
    if( itt != mcmap.end() ){
      int il =  itt->second.first/100;//layer
      int ir = (itt->second.first-il*100)/10;//row
      int im = (itt->second.first-il*100-ir*10);//module
      int ic =  itt->second.second; //channel

      //int valcheck = calib.IsGood(il,ir,im,ic);

      hmask->Fill(ir*NCHANNELS+ic,il*NMODULES+im,val);
      hmaskdist->Fill(val);

    }

  }
*/

  vector<TLine*> lines;
  for (uint i=0; i<NROWS; i++)  lines.push_back( new TLine((i+1)*NCHANNELS,            0,(i+1)*NCHANNELS,NLAYERS*NMODULES));
  for (uint i=0; i<NLAYERS; i++)lines.push_back( new TLine(             0,(i+1)*NMODULES, NROWS*NCHANNELS, (i+1)*NMODULES));


  TCanvas *ccr3 = new TCanvas("ccr3","Trk Raw",900,600);
  ccr3->Divide(2,1);

  /////////////////////////////////////////////// figura 1
  ccr3->GetPad(1)->cd();

  if(hmask)hmask->Draw("colz");

  for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");
  /////////////////////////////////////////////// figura 2
  ccr3->GetPad(2)->cd();

  if(hmaskdist)hmaskdist->Draw("");


  DrawTitle(ccr3->GetPad(1),Form("%s",file.data()));


  return;

}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



map<uint,TGraph*> GetTransferFunctions(std::string file, int detector, std::string maskfile,  bool TESTLIB){


  std::map<uint,TF1>    fmap ;
  std::map<uint,TGraph> gmap ;

  CCalib *data_calib = NULL;
  if(detector>=0)cout << endl <<"<<<<< read detector "<<detector<<" >>>>> ";

  if(TESTLIB){

    cout << endl ;
    cout << endl <<" --------------------- ";


    data_calib = new  Crane::Calibration::CCalib() ;

    // data_calib->GetTrkCalib().SetTransferFunctions(file,1,detector);
    // data_calib->GetTrkCalib().SetMasks(maskfile,detector);
    //    data_calib->GetTrkCalib().SetTransferFunctions(file,1);
    data_calib->GetTrkCalib().SetTransferFunctions(file);
    data_calib->GetTrkCalib().SetMasks(maskfile);



    fmap = data_calib->GetTrkCalib().GetTransferFunctionMap_TF1();
    gmap = data_calib->GetTrkCalib().GetTransferFunctionMap_TGraph();


    cout << endl <<" --------------------- ";
    cout << endl ;


  }else{
    std::string tf_filename = file;
    TFile* tf_file = new TFile(tf_filename.c_str(), "read");
    if (!(tf_file)) cout << endl << "Error when opening file " << tf_filename;
    if (!(tf_file->IsOpen())) cout << endl <<"problems opening file " << tf_filename;
    for(uint il=0; il<NLAYERS; il++){
      for(uint ir=0; ir<NROWS; ir++){
	for(uint im=0; im<NMODULES; im++){
	  for (uint ch=0; ch<NCHANNELS; ch++)  {

	    if( detector >=0 &&
		detector != 100*il+10*ir+im &&
		true ) continue;

	    // uint key = Crane::Calibration::CTrkCalib::GetVolumeId(il,ir,im,ch);

	    uint key = CTrkCalib::GetVolumeIdFromRaw(il,ir,im,ch);
            //	    uint key = GGeometryObject::GetTrkVolumeIdFromRaw(il,ir,im,ch);


	    TString graph_ch = Form("Layer%iRow%iModule%iCh%i",il,ir,im,ch);
	    TString func_ch = Form("TF1_Layer%iRow%iModule%iCh%i",il,ir,im,ch);

	    TGraph* gtf = (TGraph*)tf_file->Get(graph_ch);
	    TF1* func_tf = (TF1*)tf_file->Get(func_ch);

	    TGraph* anti_gtf = new TGraph();
	    if(gtf)
	      for(Int_t i=0; i<gtf->GetN(); i++)
		{anti_gtf->SetPoint(anti_gtf->GetN(), gtf->GetY()[i], gtf->GetX()[i]);}


	    if(func_tf)fmap.insert( std::pair<uint,TF1> ( key, TF1( *func_tf )  ));
	    if(gtf)gmap.insert( std::pair<uint,TGraph> ( key, TGraph( *anti_gtf )));
	    if(gtf && func_tf){ // if both are present...

	    }else{    //... if missing ...

	    };
	  }
	}
      }
    }

  }


  cout << endl <<" TF file "<<file;
  cout << endl <<" TF1 map size    "<<fmap.size();
  cout << endl <<" TGraph map size "<<gmap.size();



  map<uint,TGraph*> g;

  map<uint,pair<int,int>> mcmap = CTrkCalib::GetModuleChannelMap();

  for (std::map<uint,TGraph>::iterator it=gmap.begin(); it!=gmap.end(); ++it){

    uint volId = it->first;
    cout << endl << volId;

    auto itt = mcmap.find(volId);
    if( itt == mcmap.end() )continue;
    int l =  itt->second.first/100;//layer
    int r = (itt->second.first-l*100)/10;//row
    int m = (itt->second.first-l*100-r*10);//module
    int c =  itt->second.second; //channel

    uint key = itt->second.first;


    if( detector >=0 &&
	detector != key &&
	true ) continue;

    cout <<" - is detector "<<detector;



    bool OK = data_calib->GetTrkCalib().IsGood(l,r,m,c);
    if(!OK){
      if(detector>=0)cout << endl <<"BAD - channel "<<c<<" module "<<l*100+r*10+m;
      continue;
    }

    cout <<" - is OK ";

    if( TESTLIB &&
	data_calib->GetTrkCalib().TrackerEnergyResponseFunction(200,l,r,m,c) > 1 )
      cout << endl <<"ANOMALOUS TF - module "<<l*100+r*10+m<<" channel "<<c;


    TGraph *ggg = new TGraph();
    double max = 2048;//2047;
    int nstep = 512;//2047;
    for(uint i=0; i<(uint)nstep; i++){

      if(TESTLIB){
	double x = max/nstep*i;
	double yyy = data_calib->GetTrkCalib().TrackerEnergyResponseFunction(x,l,r,m,c);
	ggg->SetPoint(ggg->GetN(), yyy, x );
      }else{
	double mV2MeV=0.841/1000;
	double x = max/nstep*i;
	double yyy = 0; //DA CALCOLARE
	ggg->SetPoint(ggg->GetN(), yyy, x );
	cout << endl << "TESTLIB=false non implementato!!!";
      }
    }

    g.insert( pair<uint,TGraph*> (detector*100+c,ggg) );
    cout <<" - is inserted "<<detector*100+c;

  }

  cout << endl ;
  cout << endl ;
  cout << endl ;
  cout << endl<<" n.tf found "<<g.size() ;

  return g;

};

void DiffTFFile(std::string file1,std::string file2, int detector, bool TESTLIB){

  map<uint,TGraph*> tf1 = GetTransferFunctions(file1,detector,"",TESTLIB);
  cout << endl << "map size "<<tf1.size();

  map<uint,TGraph*> tf2 = GetTransferFunctions(file2,detector,"",TESTLIB);
  cout << endl << "map size "<<tf2.size();


  TH2F *htf = (TH2F*)gROOT->FindObject("htf");
  if(htf)delete htf;
  htf = new TH2F("htf","",NCHANNELS*NROWS,0,NCHANNELS*NROWS,NLAYERS*NMODULES,0,NLAYERS*NMODULES);
  htf->GetXaxis()->SetTitle(Form("row*%i+channel",NCHANNELS));
  htf->GetYaxis()->SetTitle(Form("layer*%i+module",NMODULES));
  htf->SetStats(0);

  vector<TGraph*> g;

 for (std::map<uint,TGraph*>::iterator it1=tf1.begin(); it1!=tf1.end(); ++it1){
    int key1  = it1->first;
    int il = (int)key1/10000;
    int ir = (int)(key1-il*10000)/1000;
    int im = (int)(key1-il*10000-ir*1000)/100;
    int ic = (int)(key1-il*10000-ir*1000-im*100);

    if( detector >=0 &&
	detector != 100*il+10*ir+im &&
	true ) continue;


    TGraph* g1 = it1->second;
    std::map<uint,TGraph*>::iterator it2 = tf2.find(key1);
    if (it2 != tf2.end()){


      TGraph* g2 = it2->second;

      //      cout << endl << key1<<" "<<g1<<" "<<g2;


      htf->Fill(ir*NCHANNELS+ic,il*NMODULES+im);
      //      cout << endl << ir*NCHANNELS+ic << " "<<il*NMODULES+im;

      TGraph *ggg = new TGraph();
      for(int ip=0; ip<g2->GetN(); ip++){
	double x1,y1;
	g1->GetPoint(ip,x1,y1);// x=MeV y=ADC
	double x2,y2;
	g2->GetPoint(ip,x2,y2);// x=MeV y=ADC
	// la y e` la stessa, per costruzione
	ggg->AddPoint(y1,2*(x1-x2)/(x2+x1));
	//ggg->AddPoint(y1,(x1/x2));
	//	cout << endl << y1 <<" "<<y2<<" "<<x1<<" "<<x2<<" "<<x1/x2;
      }
      g.push_back(ggg);
      //g.push_back(g1);
      //g.push_back(g2);
    }
 }




  vector<TLine*> lines;
  for (uint i=0; i<NROWS; i++)  lines.push_back( new TLine((i+1)*NCHANNELS,            0,(i+1)*NCHANNELS,NLAYERS*NMODULES));
  for (uint i=0; i<NLAYERS; i++)lines.push_back( new TLine(             0,(i+1)*NMODULES, NROWS*NCHANNELS, (i+1)*NMODULES));

  TCanvas *ccr3 = new TCanvas("ccr3","Trk Raw",900,600);
  ccr3->Divide(2,1);

  /////////////////////////////////////////////// figura 1
  ccr3->GetPad(1)->cd();

  if(htf)htf->Draw("colz");

  for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");

  DrawTitle(ccr3->GetPad(1),file1.data(),0.07,0.95);
  DrawTitle(ccr3->GetPad(1),file2.data(),0.07,0.95-0.04);


  /////////////////////////////////////////////// figura 2
  ccr3->GetPad(2)->cd();
  //  ccr3->GetPad(2)->SetLogy();

  double mevmin = 1e-2;
  double mevmax = 200;
  double adcmax = 2047;
  // TH1F* h2 = new TH1F("h2","",10,0.001,mevmax);
  TH1F* h2 = new TH1F("h2","",10,0.,adcmax);

  h2->GetYaxis()->SetTitle("dE/E");
  h2->GetXaxis()->SetTitle("S (ADC counts) ");
  h2->GetYaxis()->SetTitleOffset(1.5);
  // h2->GetYaxis()->SetRangeUser(0.9,1.1);
  h2->GetYaxis()->SetRangeUser(-0.2,0.2);
  h2->SetStats(0);
  // h2->SetMaximum(adcmax);
  h2->Draw();
  for(auto gg: g)gg->Draw("pl");

};


void CheckTFDerivative(std::string file, int detector, int channel   , bool TESTLIB     ){


  map<uint,TGraph*> gmap = GetTransferFunctions(file,detector,"",TESTLIB);
  vector<TGraph*> g;
  vector<TGraph*> gde;

  for (std::map<uint,TGraph*>::iterator it=gmap.begin(); it!=gmap.end(); ++it){

    TGraph* g1 = it->second;
    ///    TGraph* g1inv = new TGraph();

    //    int ch = GGeometryObject::GetModuleChannel(it->first);


    map<uint,pair<int,int>> mcmap = CTrkCalib::GetModuleChannelMap();
    auto itt = mcmap.find(it->first);
    if( itt == mcmap.end() )continue;

    int ch =  itt->second.second; //channel


    if(ch!=channel)continue;


    g.push_back( g1 );



    TGraph *gd = new TGraph();
    for(int ip=1; ip<g1->GetN(); ip++){
	double x1,y1;
	g1->GetPoint(ip-1,x1,y1);// x=MeV y=ADC
	double x2,y2;
	g1->GetPoint(ip,x2,y2);// x=MeV y=ADC
	double der = (x2-x1)/(y2-y1);
	gd->AddPoint((x1+x2)/2,der);
    }

    //   TSpline3* spline = new TSpline3("spline", graph);


    gde.push_back( gd );

  }

  TCanvas *ccr3 = new TCanvas("ccr3","Trk Raw",900,600);
  ccr3->Divide(2,1);


  ccr3->GetPad(1)->cd();
  ccr3->GetPad(1)->SetLogx();


  double mevmin = 1e-2;
  double mevmax = 200;
  double adcmax = 2047;
  TH1F* h2 = new TH1F("h2","",10,0.001,mevmax);

  h2->GetXaxis()->SetTitle("MeV");
  h2->GetYaxis()->SetTitle("S (ADC counts) ");
  h2->GetYaxis()->SetTitleOffset(1.5);
  h2->GetYaxis()->SetRangeUser(0.001,100);
  h2->SetStats(0);
  h2->SetMaximum(adcmax);
  h2->Draw();
  for(auto gg: g)gg->Draw("pl");

  ccr3->GetPad(2)->cd();
  ccr3->GetPad(2)->SetLogx();
  ccr3->GetPad(2)->SetLogy();

  TH1F* hh2 = new TH1F("hh2","",10,0.001,mevmax);

  hh2->GetXaxis()->SetTitle("MeV");
  hh2->GetYaxis()->SetTitle("dE/dS (MeV/ADC counts) ");
  hh2->GetYaxis()->SetTitleOffset(1.5);
  hh2->GetYaxis()->SetRangeUser(0.0001,1);
  hh2->SetStats(0);
  hh2->SetMaximum(adcmax);
  hh2->Draw();
  for(auto gg: gde)gg->Draw("pl");

  DrawTitle(ccr3->GetPad(1),Form("%s %i",file.data(),detector));


}


/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////



void writeMaskToFile(bool mask[NLAYERS][NROWS][NMODULES][NCHANNELS], const std::string& filename) {
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    for (int i = 0; i < NLAYERS; ++i) {
        for (int j = 0; j < NROWS; ++j) {
            for (int k = 0; k < NMODULES; ++k) {
                // Calcola il valore per la prima colonna
                int index = i * 100 + j * 10 + k;
                // outFile << std::setw(4) <<index << " ";
                outFile << i<<j<<k << " ";

                // Crea il campo di bit per mask[i][j][k][c] in esadecimale
                std::bitset<NCHANNELS> bitField;
                for (int c = 0; c < NCHANNELS; ++c) {
                    bitField[c] = mask[i][j][k][c];
                }
		outFile << "0x" << std::setw(8) << std::setfill('0') << std::hex << std::uppercase << bitField << std::endl;

		//                outFile << "0x" << std::hex << std::setw(8) << bitField.to_ulong() << std::dec<< std::endl;
            }
        }
    }

    outFile.close();
}




void CheckStripRate(TString flist, int nev, float fraction_max ,TString ddir,TString suffix , string inputfile ){

  TChain *craw = GetChain(flist.Data(),ddir,"TreeRaw",suffix);
  if(!craw)return;

  Crane::Calibration::CRawTrk *rawtrk = new Crane::Calibration::CRawTrk();
  craw->SetBranchAddress("Trk", &rawtrk);


  CTrkCalib calib;
  calib.SetMasks(inputfile);
  map<uint,bool> mmap = calib.GetMaskMap();
  cout << endl << " Mask file "<<inputfile;
  cout << endl << " map size  "<<mmap.size();


  //////////////////////////
  // histograms
  //////////////////////////

  vector<TLine*> lines;
  for (uint i=0; i<NROWS; i++)  lines.push_back( new TLine((i+1)*NCHANNELS,            0,(i+1)*NCHANNELS,NLAYERS*NMODULES));
  for (uint i=0; i<NLAYERS; i++)lines.push_back( new TLine(             0,(i+1)*NMODULES, NROWS*NCHANNELS, (i+1)*NMODULES));

  TH2F *hmap = (TH2F*)gROOT->FindObject("hmap");
  if(hmap)delete hmap;
  hmap = new TH2F("hmap","Hit fraction",NCHANNELS*NROWS,0,NCHANNELS*NROWS,NLAYERS*NMODULES,0,NLAYERS*NMODULES);
  hmap->GetXaxis()->SetTitle(Form("row*%i+channel",NCHANNELS));
  hmap->GetYaxis()->SetTitle(Form("layer*%i+module",NMODULES));
  hmap->SetStats(0);

  TH1F *hmapdist = (TH1F*)gROOT->FindObject("hmapdist");
  if(hmapdist)delete hmapdist;
  // hmapdist = new TH1F("hmapdist","Hit fraction",1000,0,0.1);
  hmapdist = new TH1F("hmapdist","Hit fraction",1000,0,0.01);

  TH1F *hitdist = (TH1F*)gROOT->FindObject("hitdist");
  if(hitdist)delete hitdist;
  // hitdist = new TH1F("hitdist","#hits per event",12000,0,12000);
  hitdist = new TH1F("hitdist","#hit per event",1000,0,1000);

  TH1F *hadc = (TH1F*)gROOT->FindObject("hadc");
  if(hadc)delete hadc;
  hadc = new TH1F("hadc","Raw signal (ADC counts)",500,0,2050);


  //////////////////////////
  // counters
  //////////////////////////

  int  nhit[NLAYERS][NROWS][NMODULES][NCHANNELS];
  fill_n( &nhit[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,0);

  bool  mask[NLAYERS][NROWS][NMODULES][NCHANNELS];
  fill_n( &mask[0][0][0][0],NLAYERS*NROWS*NMODULES*NCHANNELS,true);

  //////////////////////////
  // loop
  //////////////////////////
  std::cout<< std::endl <<" Hit fraction cut "<<fraction_max;

  cout << endl << "OOOOOOO Loop over the events...";
  cout << endl ;
  int ntot = 0;
  int ntrk = 0;
  for (uint i=0; i<craw->GetEntries(); i++){

    if( i > nev )break;

    if( i%1000==0 )std::cout<<"*";

    craw->GetEntry(i);

    // ----------------------------
    // selection
    // ----------------------------

    if( rawtrk->flag!=7 )continue;
    if( rawtrk->eventid==0 )continue;

    ntot++;

    // cout << endl << rawtrk->adcdata.size();
    // ----------------------------
    // conta le hit per ogni modulo
    // ----------------------------
    int ngoodhit =0;
    for(uint i=0; i<(uint)rawtrk->adcdata.size(); i++){
      int layer  = rawtrk->layer[i];
      int row    = rawtrk->row[i];
      int module = rawtrk->module[i];
      int channel = rawtrk->channel[i];

      //-------------------------------------
      bool good = true;
      //      uint key = GGeometryObject::GetTrkVolumeIdFromRaw(layer,row,module,channel);
      uint key = CTrkCalib::GetVolumeIdFromRaw(layer,row,module,channel);
      std::map<uint,bool>::iterator it = mmap.find(key);
      if (it != mmap.end()){
	good = (it->second);
      }else{
	// cout << endl << "missing mask key "<<key;
      }
      //-------------------------------------

      if(!good)continue;

      ngoodhit++;
      nhit[layer][row][module][channel]++;
      hadc->Fill(rawtrk->adcdata[i]);
    }

    hitdist->Fill(ngoodhit);

  }
  std::cout << std::endl << " Ntot "<<ntot;
  if(ntot==0)return ;
  //////////////////////////
  // end loop
  //////////////////////////


  for(int il=0; il<NLAYERS; il++){
    for(int ir=0; ir<NROWS; ir++){
      for(int im=0; im<NMODULES; im++){
  	for(int ic=0; ic<NCHANNELS; ic++){
  	  if( nhit[il][ir][im][ic] > 0 ){
	    double fraction = (double)nhit[il][ir][im][ic]/(double)ntot;
	    // cout << endl << nhit[il][ir][im][ic];
	    hmap->Fill( ir*NCHANNELS+ic,il*NMODULES+im, fraction );
	    hmapdist->Fill(fraction);

	    //============================
	    // mask condition
	    //============================
	    if( fraction > fraction_max )mask[il][ir][im][ic]=false;

	  }
	}
      }
    }
  }

  std::string filename = "mask_output.txt";
  std::cout<< std::endl <<" Write mask to file "<<filename;
  writeMaskToFile(mask, filename);

  /////////////////////////////////////////////// figura 1
  TCanvas *ccr3 = new TCanvas("ccr3","Trk Raw",900,600);
  ccr3->Divide(2,1);
  ccr3->GetPad(1)->cd();
  ccr3->GetPad(1)->SetLogz();

  if(hmap)hmap->GetZaxis()->SetRangeUser(0.0001,1);
  if(hmap)hmap->Draw("colz");

  for( uint i=0; i<lines.size(); i++)lines.at(i)->Draw("same");
  ccr3->GetPad(2)->cd();
  ccr3->GetPad(2)->SetLogy();

  if(hmapdist)hmapdist->Draw("");


  DrawTitle(ccr3->GetPad(1),Form("%s + %s",flist.Data(),inputfile.c_str()));


  /////////////////////////////////////////////// figura 2
  TCanvas *ccr = new TCanvas("ccr","Trk Raw",900,600);
  ccr->Divide(2,1);
  ccr->GetPad(1)->cd();
  ccr->GetPad(1)->SetLogy();

  if(hadc)hadc->Draw("");

  ccr->GetPad(2)->cd();
  ccr->GetPad(2)->SetLogy();

  if(hitdist)hitdist->Draw("");


  DrawTitle(ccr->GetPad(1),Form("%s + %s",flist.Data(),inputfile.c_str()));


}
