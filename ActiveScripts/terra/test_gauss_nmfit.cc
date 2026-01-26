#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <vector>
#include <iterator>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TH1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TString.h"
#include "TRegexp.h"
#include <TLeaf.h>
#include <TApplication.h>
#include <Rtypes.h>
#include <TSystem.h>
#include <TLine.h>
#include <TStyle.h>
#include <TF1.h>
#include <TText.h>
#include <TLatex.h>
#include <TMinuit.h>

using namespace std;

//Might as well start off ambitious, let's try to make a function for any number of strips?
//It might be best to just have that me this function...? There's a fit involved. Don't want to have a variable for every strip. Can do what was done for LRMS then
//

void test_gauss_nmfit()
{
    //There are the variables that change every time!!
    int detnum = 136;
    int run = 001;
    int temp = 37; //Temperature in negative c, sorry, change me if needed
    string strips[] = {"A2","E","F1","F2","G"};

    //string strips[] = {"D","C","B","A"};
    //Variable strip and strip names. It's okay if some are missing and it's okay if some are not found. It is also okay if the order is different than the .dat file.
    //Colors will be based on the order of the the letters in the strips variable.
    //string strips[] = {"A","B","C","potato"};

    //You can change axis ranges if you want:
    float xmin = 0.5;
    float xmax = 20;
    float ymin = 1;
    float ymax = 20;

    //End change variables

    const int nstrips = sizeof(strips) / sizeof(strips[0]); //This is the number of strips that will be fitted based on the user-provided list
    int npstrip[nstrips] = { 0 };
    double I[nstrips];
    double eI[nstrips];
    double Af[nstrips];
    double eAf[nstrips];
    double Rs[nstrips];
    double eRs[nstrips];

    double T = static_cast<double>(-1*temp) + 273.;

    //This could be an input in an executable. It might help with human error actually to leave this as a .cc file though given all the changes.
    TFile *fin = new TFile(TString::Format("summary-%i.root",detnum),"read");


    int np =0;
    /*double T = -43+273.;    // temperature for testing, -37C */
    double q = 1.6e-19; // electron charge
    double k = 1.38e-23; // Boltzmann constant
    double eps = 3.6;    // ionization energy of silicon, eV
    double Rp = 100e6;   // parallel resistance of preamp, 100 MOhm
    double gm = 18e-3;   // transconductance in FET, 18 ms
    double Bita = 1;
    double factor = (2.355*eps*1e-3/q)*(2.355*eps*1e-3/q);

    double Af_sim = 7.3e-14; //simulated ASIC Af
    double Ieff =2.5*1e-9; //A
    double Sw = 0.54*1e-18;
    double Fi = 0.45;   // noise form factor, gauss
    double Fv = 1.02;    // noise form factor, gauss
    double Fvf = 0.52;  // noise form factor, gauss
    double pi = 3.1415926;
    double enc2 = 0.;
    double fwhm = 0.;

    double Ctot = 140e-12; //pF double the value of 8-strip C = eA/d

    TTree *tin = (TTree*)fin->Get("tree");

    np=tin->GetEntries(); //Number of entries for the tree.
    string sn;
    string strip;
    float tpeak;
    float res;
    float res_err;
    std::string *m_sn = new std::string;
    std::string *m_strip = new std::string;
    tin->SetBranchAddress("sn", &m_sn);
    tin->SetBranchAddress("strip", &m_strip);
    tin->SetBranchAddress("tpeak", &tpeak);
    tin->SetBranchAddress("res", &res);
    tin->SetBranchAddress("res_err", &res_err);

    //Success!
    TGraphErrors *g[nstrips];
    TF1 *f[nstrips];
    for(int j = 0;j<nstrips;j++){
        f[j] = new TF1(TString::Format("f%s",strips[j].c_str()),"sqrt([0]*x*1e-6+[1]/(x*1e-6)+[2])",xmin,xmax);
        f[j]->SetLineColor(j+1);
        f[j]->SetParameters(5e5, 1e-5, 1);
        g[j] = new TGraphErrors(36);
    }

    //Iterate over the summary.root data
    for(int i=0; i<np; i++){
        tin->GetEntry(i);

        //for every entry cross-check the name
        for(int j=0;j<nstrips;j++){
            //if(m_strip->compare("stripA")==0){
            if(m_strip->compare("strip"+strips[j])==0){
                npstrip[j]++;
                g[j]->SetPoint(npstrip[j]-1, tpeak, res);
                g[j]->SetPointError(npstrip[j]-1, 0, res_err);
            }
        }
    }

    //Make a canvas for plotting
    TCanvas *myc =new TCanvas("myc","",1000,800);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    myc->SetTitle(0);
    myc->cd();
    myc->SetLogx();
    myc->SetLogy();
    myc->SetGridx();
    myc->SetGridy();

    auto legend = new TLegend(0.15,0.7,0.88,0.89);

    //Set all of the graph parameters on the first plot and the others are drawn on the same plot

    g[0]->Draw();
    g[0]->SetLineWidth(0);
    g[0]->SetMarkerColor(1);
    g[0]->SetMarkerStyle(3);
    g[0]->SetTitle(0);
    g[0]->GetXaxis()->CenterTitle();
    g[0]->GetXaxis()->SetTitle("Peaking time [#mus]");
    g[0]->GetXaxis()->SetTitleOffset(1.3);
    g[0]->GetYaxis()->CenterTitle();
    g[0]->GetYaxis()->SetTitle("X-ray FWHM [keV]");
    g[0]->Fit(f[0],"","same", xmin,xmax);

    auto xaxis = g[0]->GetXaxis();
    auto yaxis = g[0]->GetYaxis();
    xaxis->SetMoreLogLabels();
    yaxis->SetMoreLogLabels();
    xaxis->SetLimits(xmin,xmax);
    g[0]->GetHistogram()->SetMinimum(2.);
    g[0]->GetHistogram()->SetMaximum(20.);

    if(nstrips > 1){ //No error even if nstrips = 1.
        for(int j=1;j<nstrips;j++){
            g[j]->Draw("sameP");
            g[j]->SetLineWidth(0);
            g[j]->SetMarkerColor(j+1);
            g[j]->SetMarkerStyle(3);
            g[j]->Fit(f[j],"","same", xmin,xmax);
        }
    }


    ofstream fout;
    //fout.open;
    fout.open ("test");

    for(int j = 0;j<nstrips;j++){
        I[j] = (f[j]->GetParameter(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
        eI[j] = (f[j]->GetParError(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
        Af[j] = f[j]->GetParameter(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
        eAf[j] = f[j]->GetParError(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
        Rs[j] = f[j]->GetParameter(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;
        eRs[j] = f[j]->GetParError(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;
        cout << strips[j] <<fixed<<setprecision(2)<< " Current " << I[j]<< " Af " << Af[j] <<" Rs "<<Rs[j]<<" "<<endl;
        fout<< "Strip " << strips[j] << "\t" << fixed<<setprecision(2)<<I[j]<<" "<<Af[j]<<" "<<Rs[j]<<" "<<endl;
        legend->AddEntry(g[j], TString::Format("Sh%04i_%s: Ileak =%10.2f nA, Af =%10.2f#times10^{-13} V^{-2}, Rs =%10.2f #Omega ",detnum,(("strip"+strips[j]).c_str()),I[j],Af[j],Rs[j]),"lep");
    }

    fout.close();
    legend->Draw();

    myc->SaveAs(TString::Format("det%i-run%i-%iC-nmfit.pdf",detnum,run,temp));
    myc->SaveAs(TString::Format("det%i-run%i-%iC-nmfit.png",detnum,run,temp));

    //cout << "Test of format " << TString::Format("strip%s",strips[0]) << endl;
    //cout << "Test of format " << "strip"+strips[0] << endl;
    //cout << "Nstrips = " << nstrips << endl;
    //cout << "Strip[0] is " << strips[0] << endl;
    //cout << "npstrips[0] = " << npstrip[0] << endl;
    //cout << "Hello World!!" << endl;

}
