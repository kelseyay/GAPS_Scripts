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

using namespace std;
//double data[100000][7];

int main(int argc, char** argv)
{
    if (argc!=3)
    {
	cout<<"Syntax: "<<argv[0]<<" "<<"[input file] [output file]"<<endl;
	exit(1);
    }
    int np =0.;
    char c;
    ifstream myfile1(argv[1]);
    if (!myfile1)
    {
        cerr<<" cannot open the data file" <<endl;
        return -1;
    }
    while (myfile1.get(c))
    {
        if (c=='\n') np++;
    }
    cout<<"***total data points: "<<np<<"****"<<endl;

    ifstream myfile(argv[1]);
    vector<string> data;
    string line;
    while (getline(myfile, line))
    {
        data.push_back(line);
    }
    vector<string> col0;
    vector<string> col1;
    vector<float> col2;
    vector<float> col3;
    vector<float> col4;
    for (auto it=data.begin(); it!=data.end(); it++)
    {
//        cout<< *it<<endl;
        istringstream is(*it);
        string s;
        int pam=0;
        while (is>>s)
        {
            if(pam ==0)
            {
                col0.push_back(s);
            }
            if(pam==1)
            {
                col1.push_back(s);
            }
            if(pam==2)
            {
                float a2 = atof(s.c_str());
                col2.push_back(a2);
            }
            if(pam==3)
            {
                float a3 = atof(s.c_str());
                col3.push_back(a3);
            }
            if(pam==4)
            {
                float a4 = atof(s.c_str());
                col4.push_back(a4);
            }
            pam ++;
        }
    }
    double Af_sim = 7.3e-14; //simulated ASIC Af
    double q = 1.6e-19; // electron charge
    double eps = 3.6; //eV
    double Ieff =2.5*1e-9; //A
    double Sw = 0.54*1e-18;
    double Fi = 0.64;   // noise form factor, ASIC
    double Fv = 0.853;    // noise form factor, ASIC
    double Fvf = 0.543;  // noise form factor, ASIC
    double pi = 3.1415926;
    double enc2 = 0.;
    double fwhm = 0.;

    TFile *fout = new TFile(argv[2], "recreate");
    TTree *tree = new TTree("tree", "");
    string sn;
    string strip;
    float tpeak;
    float res;
    float res_err;
    tree->Branch("sn", &sn);
    tree->Branch("strip", &strip);
    tree->Branch("tpeak", &tpeak, "tpeak/F");
    tree->Branch("res", &res, "res/F");
    tree->Branch("res_err", &res_err, "res_err/F");
    for (int i=0; i<np; i++)
    {
        sn = col0[i];
        strip = col1[i];
        tpeak = col2[i]/1000.; //ns-->us
//        cout<<"strip: peaking time: "<<strip<<" "<<tpeak<<endl;
        res = col3[i];
        res_err = col4[i];
        tree->Fill();
    }
    tree->Write();
    tree->Delete();
    fout->Close();

    return 0;
}
