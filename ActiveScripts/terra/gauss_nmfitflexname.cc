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
//double data[100000][7];

//int main(int argc, char** argv)
void gauss_nmfitflexname()
//int main()
{
//    if (argc!=3)
//    {
//	cout<<"Syntax: "<<argv[0]<<" "<<"[input file] [output file]"<<endl;
//	exit(1);
//    }
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
    int detnum = 788;
    int run = 001;
    int temp = 42; //Temperature, change me if needed
    double T = static_cast<double>(-1*temp) + 273.;

    //This should be an input in an executable
    //TFile *fin = new TFile(TString::Format("/home/kelsey/fittest/det%i-gauss-run%i-%iC/summary-%i.root",detnum,run,temp,detnum),"read");
    TFile *fin = new TFile(TString::Format("summary-%i.root",detnum),"read");

    TTree *tin = (TTree*)fin->Get("tree");
    np=tin->GetEntries();
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

    auto ga = new TGraphErrors(36);
    auto gb = new TGraphErrors(36);
    auto gc = new TGraphErrors(36);
    auto gd = new TGraphErrors(36);
    //auto ge = new TGraphErrors(36);
    //auto gf = new TGraphErrors(36);
    //auto gg = new TGraphErrors(36);
    //auto gh = new TGraphErrors(36);

    string aname = "stripA";
    string bname = "stripB";
    string cname = "stripC";
    string dname = "stripD";

    //string ename = "stripE";
    //string fname = "stripF";
    //string gname = "stripG";
    //string hname = "stripH";

    int npa = 0;
    int npb = 0;
    int npc = 0;
    int npd = 0;
    //int npe = 0;

    //int npf = 0;
    //int npg = 0;
    //int nph = 0;

    for(int i=0; i<np; i++)
    {
        tin->GetEntry(i);


        if(m_strip->compare(aname)==0)
        {
           npa++;
           ga->SetPoint(npa-1, tpeak, res);
           ga->SetPointError(npa-1, 0, res_err);
        }

        if(m_strip->compare(bname)==0)
        {
           npb++;
           gb->SetPoint(npb-1, tpeak, res);
           gb->SetPointError(npb-1, 0, res_err);
        }
        if(m_strip->compare(cname)==0)
        {
           npc++;
           gc->SetPoint(npc-1, tpeak, res);
           gc->SetPointError(npc-1, 0, res_err);
        }
        if(m_strip->compare(dname)==0)
        {
           npd++;
           gd->SetPoint(npd-1, tpeak, res);
           gd->SetPointError(npd-1, 0, res_err);
        }

	/*
        if(m_strip->compare(ename)==0)
        {
           npe++;
           ge->SetPoint(npe-1, tpeak, res);
           ge->SetPointError(npe-1, 0, res_err);
        }

        if(m_strip->compare(fname)==0)
        {
           npf++;
           gf->SetPoint(npf-1, tpeak, res);
           gf->SetPointError(npf-1, 0, res_err);
        }

        if(m_strip->compare(gname)==0)
        {
           npg++;
           gg->SetPoint(npg-1, tpeak, res);
           gg->SetPointError(npg-1, 0, res_err);
        }

        if(m_strip->compare(hname)==0)
        {
           nph++;
           gh->SetPoint(nph-1, tpeak, res);
           gh->SetPointError(nph-1, 0, res_err);
        }
	*/
    }

    TCanvas *myc =new TCanvas("myc","",1000,800);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    myc->SetTitle(0);
    myc->cd();
    myc->SetLogx();
    myc->SetLogy();
    myc->SetGridx();
    myc->SetGridy();



    ga->Draw();
    ga->SetLineWidth(0);
    ga->SetMarkerColor(1);
    ga->SetMarkerStyle(3);


    //gb->Draw();


    gb->Draw("sameP");
    gb->SetLineWidth(0);
    gb->SetMarkerColor(2);
    gb->SetMarkerStyle(3);


    gc->Draw("sameP");
    gc->SetLineWidth(0);
    gc->SetMarkerColor(3);
    gc->SetMarkerStyle(3);


    gd->Draw("sameP");
    gd->SetLineWidth(0);
    gd->SetMarkerColor(4);
    gd->SetMarkerStyle(3);

	/*
    ge->Draw("sameP");
    ge->SetLineWidth(0);
    ge->SetMarkerColor(kCyan+2);
    ge->SetMarkerStyle(3);

    gf->Draw("sameP");
    gf->SetLineWidth(0);
    gf->SetMarkerColor(6);
    gf->SetMarkerStyle(3);

    gg->Draw("sameP");
    gg->SetLineWidth(0);
    gg->SetMarkerColor(7);
    gg->SetMarkerStyle(3);

    gh->Draw("sameP");
    gh->SetLineWidth(0);
    gh->SetMarkerColor(8);
    gh->SetMarkerStyle(3);
    */

    auto xaxis = ga->GetXaxis(); // If strip A is off: CHANGE ga TO gb in this line and the one below -IAN
    auto yaxis = ga->GetYaxis();
    xaxis->SetMoreLogLabels();
    yaxis->SetMoreLogLabels();
    xaxis->SetLimits(0.4,31);


    ga->SetTitle(0); // Also change here -IAN
    ga->GetHistogram()->SetMinimum(2.);
    ga->GetHistogram()->SetMaximum(20.);
    ga->GetXaxis()->CenterTitle();
    ga->GetXaxis()->SetTitle("Peaking time [#mus]");
    ga->GetXaxis()->SetTitleOffset(1.3);
    ga->GetYaxis()->CenterTitle();
    ga->GetYaxis()->SetTitle("X-ray FWHM [keV]");


    TF1 *fa=new TF1("fa","sqrt([0]*x*1e-6+[1]/(x*1e-6)+[2])",0,40);
    fa->SetLineColor(1);
    fa->SetParameters(5e5, 1e-5, 1);
    ga->Fit("fa","","same", 0.5,30);


    /*TCanvas *c2 = new TCanvas("c2","contours",10,10,600,800);
    c2->cd();

    TGraph *gr12 = (TGraph*)gMinuit->Contour(40,0,1);
    gr12->Draw("alp");
    Int_t n = gr12->GetN();
    Double_t ax[n],ay[n];
    for(Int_t i=0; i<n; i++){
      gr12->GetPoint(i,ax[i],ay[i]);
      std::cout<<ax[i]<<" "<<ay[i]<<std::endl;
    };
      TGraph* gr13 = new TGraph(n,ax[i],ay[i]);
	TCanvas *c2 = new TCanvas("c2","contours",10,10,600,800);
      c2->cd();


    myc->cd();
    */



    TF1 *fb=new TF1("fb","sqrt([0]*x*1e-6+[1]/(x*1e-6)+[2])",0,40);
    fb->SetLineColor(2);
    fb->SetParameters(5e5, 1e-5, 1);
    gb->Fit("fb","","same", 0.5,30);



    TF1 *fc=new TF1("fc","sqrt([0]*x*1e-6+[1]/(x*1e-6)+[2])",0,40);
    fc->SetLineColor(3);
    fc->SetParameters(5e5, 1e-5, 1);
    gc->Fit("fc","","same", 0.5,30);


    TF1 *fd=new TF1("fd","sqrt([0]*x*1e-6+[1]/(x*1e-6)+[2])",0,40);
    fd->SetLineColor(4);
    fd->SetParameters(5e5, 1e-5, 1);
    gd->Fit("fd","","same", 0.5,30);

	/*
    TF1 *fe=new TF1("fe","sqrt([0]*x*1e-6+[1]/(x*1e-6)+[2])",0,40);
    fe->SetLineColor(kCyan+2);
    fe->SetParameters(5e5, 1e-5, 1);
    ge->Fit("fe","","same", 0.5,30);

    TF1 *ff=new TF1("ff","sqrt([0]*x*1e-6+[1]/(x*1e-6)+[2])",0,40);
    ff->SetLineColor(6);
    ff->SetParameters(5e5, 1e-5, 1);
    gf->Fit("ff","","same", 0.5,30);

    TF1 *fg=new TF1("fg","sqrt([0]*x*1e-6+[1]/(x*1e-6)+[2])",0,40);
    fg->SetLineColor(7);
    fg->SetParameters(5e5, 1e-5, 1);
    gg->Fit("fg","","same", 0.5,30);

    TF1 *fh=new TF1("fh","sqrt([0]*x*1e-6+[1]/(x*1e-6)+[2])",0,40);
    fh->SetLineColor(8);
    fh->SetParameters(5e5, 1e-5, 1);
    gh->Fit("fh","","same", 0.5,30);
    */
    double Ctot = 140e-12; //pF double the value of 8-strip C = eA/d
    double Ia, Ib, Ic, Id;//, Ie, If, Ig, Ih;
    double eIa, eIb, eIc, eId;//, eIe, eIf, eIg, eIh;
    double Afa, Afb, Afc, Afd;//, Afe, Aff, Afg, Afh;
    double eAfa, eAfb, eAfc, eAfd;//, eAfe, eAff, eAfg, eAfh;
    double Rsa, Rsb, Rsc, Rsd;//, Rse, Rsf, Rsg, Rsh;
    double eRsa, eRsb, eRsc, eRsd;//, eRse, eRsf, eRsg, eRsh;


    Ia = (ga->GetFunction("fa")->GetParameter(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
    eIa = (ga->GetFunction("fa")->GetParError(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
    Afa = ga->GetFunction("fa")->GetParameter(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
    eAfa = ga->GetFunction("fa")->GetParError(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
    Rsa = ga->GetFunction("fa")->GetParameter(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;
    eRsa = ga->GetFunction("fa")->GetParError(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;


    Ib = (gb->GetFunction("fb")->GetParameter(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
    eIb = (gb->GetFunction("fb")->GetParError(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
    Afb = gb->GetFunction("fb")->GetParameter(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
    eAfb = gb->GetFunction("fb")->GetParError(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
    Rsb = gb->GetFunction("fb")->GetParameter(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;
    eRsb = gb->GetFunction("fb")->GetParError(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;


    Ic = (gc->GetFunction("fc")->GetParameter(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
    eIc = (gc->GetFunction("fc")->GetParError(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
    Afc = gc->GetFunction("fc")->GetParameter(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
    eAfc = gc->GetFunction("fc")->GetParError(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
    Rsc = gc->GetFunction("fc")->GetParameter(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;
    eRsc = gc->GetFunction("fc")->GetParError(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;


    Id = (gd->GetFunction("fd")->GetParameter(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
    eId = (gd->GetFunction("fd")->GetParError(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
    Afd = gd->GetFunction("fd")->GetParameter(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
    eAfd = gd->GetFunction("fd")->GetParError(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
    Rsd = gd->GetFunction("fd")->GetParameter(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;
    eRsd = gd->GetFunction("fd")->GetParError(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;

    /*
    Ie = (ge->GetFunction("fe")->GetParameter(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
    eIe = (ge->GetFunction("fe")->GetParError(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
    Afe = ge->GetFunction("fe")->GetParameter(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
    eAfe = ge->GetFunction("fe")->GetParError(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
    Rse = ge->GetFunction("fe")->GetParameter(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;
    eRse = ge->GetFunction("fe")->GetParError(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;

    If = (gf->GetFunction("ff")->GetParameter(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
    eIf = (gf->GetFunction("ff")->GetParError(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
    Aff = gf->GetFunction("ff")->GetParameter(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
    eAff = gf->GetFunction("ff")->GetParError(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
    Rsf = gf->GetFunction("ff")->GetParameter(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;
    eRsf = gf->GetFunction("ff")->GetParError(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;

    Ig = (gg->GetFunction("fg")->GetParameter(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
    eIg = (gg->GetFunction("fg")->GetParError(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
    Afg = gg->GetFunction("fg")->GetParameter(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
    eAfg = gg->GetFunction("fg")->GetParError(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
    Rsg = gg->GetFunction("fg")->GetParameter(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;
    eRsg = gg->GetFunction("fg")->GetParError(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;

    Ih = (gh->GetFunction("fh")->GetParameter(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
    eIh = (gh->GetFunction("fh")->GetParError(0)/factor/Fi-4*k*T/Rp)/2./q*1e9;
    Afh = gh->GetFunction("fh")->GetParameter(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
    eAfh = gh->GetFunction("fh")->GetParError(2)/factor/Ctot/Ctot/Fvf/2./pi*1e13;
    Rsh = gh->GetFunction("fh")->GetParameter(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;
    eRsh = gh->GetFunction("fh")->GetParError(1)/factor/Fv/Ctot/Ctot/(4.*k*T)-Bita/gm;
    */
    cout<<fixed<<setprecision(2)<<Ia<<" "<<Afa<<" "<<Rsa<<" "<<endl;
    cout<<fixed<<setprecision(2)<<Ib<<" "<<Afb<<" "<<Rsb<<" "<<endl;
    cout<<fixed<<setprecision(2)<<Ic<<" "<<Afc<<" "<<Rsc<<" "<<endl;
    cout<<fixed<<setprecision(2)<<Id<<" "<<Afd<<" "<<Rsd<<" "<<endl;
    // cout<<fixed<<setprecision(2)<<Ie<<" "<<Afe<<" "<<Rse<<" "<<endl;
    //cout<<fixed<<setprecision(2)<<If<<" "<<Aff<<" "<<Rsf<<" "<<endl;
    //cout<<fixed<<setprecision(2)<<Ig<<" "<<Afg<<" "<<Rsg<<" "<<endl;
    //cout<<fixed<<setprecision(2)<<Ih<<" "<<Afh<<" "<<Rsh<<" "<<endl;

    ofstream fout;
    fout.open (TString::Format("/home/kelsey/calib_results/waveform/det%i-gauss-run%i-%iC/noisepars-gauss-det%i-run%i.txt",detnum,run,temp,detnum,run));

    fout<<"stripA  " <<fixed<<setprecision(2)<<Ia<<" "<<Afa<<" "<<Rsa<<" "<<endl;
    fout<<"stripB  "<<fixed<<setprecision(2)<<Ib<<" "<<Afb<<" "<<Rsb<<" "<<endl;
    fout<<"stripC  "<<fixed<<setprecision(2)<<Ic<<" "<<Afc<<" "<<Rsc<<" "<<endl;
    fout<<"stripD  "<<fixed<<setprecision(2)<<Id<<" "<<Afd<<" "<<Rsd<<" "<<endl;
    //fout<<"stripE  "<<fixed<<setprecision(2)<<Ie<<" "<<Afe<<" "<<Rse<<" "<<endl;

    /*
    fout<<"stripF  "<<fixed<<setprecision(2)<<If<<" "<<Aff<<" "<<Rsf<<" "<<endl;
    fout<<"stripG  "<<fixed<<setprecision(2)<<Ig<<" "<<Afg<<" "<<Rsg<<" "<<endl;
    fout<<"stripH  "<<fixed<<setprecision(2)<<Ih<<" "<<Afh<<" "<<Rsh<<" "<<endl;
    */
    fout.close();


    auto legend = new TLegend(0.15,0.7,0.88,0.89);
//    legend->AddEntry(ga, TString::Format("Sh0643_stripA: Ileak=(%10.2f+/-%10.2f)nA, Af=(%10.2f+/-%10.2f)#times10^{-13} V^{-2}, Rs=(%10.2f+/-%10.2f)#Omega",Ia,eIa,Afa,eAfa,Rsa,eRsa),"lep");
    legend->AddEntry(ga, TString::Format("Sh%04i_%s: Ileak =%10.2f nA, Af =%10.2f#times10^{-13} V^{-2}, Rs =%10.2f #Omega ",detnum,aname.c_str(),Ia,Afa,Rsa),"lep");
    legend->AddEntry(gb, TString::Format("Sh%04i_B: Ileak =%10.2f nA, Af =%10.2f#times10^{-13} V^{-2}, Rs =%10.2f #Omega ",detnum,Ib,Afb,Rsb),"lep");
    legend->AddEntry(gc, TString::Format("Sh%04i_C: Ileak =%10.2f nA, Af =%10.2f#times10^{-13} V^{-2}, Rs =%10.2f #Omega ",detnum,Ic,Afc,Rsc),"lep");
    legend->AddEntry(gd, TString::Format("Sh%04i_D: Ileak =%10.2f nA, Af =%10.2f#times10^{-13} V^{-2}, Rs =%10.2f #Omega ",detnum,Id,Afd,Rsd),"lep");
    //legend->AddEntry(ge, TString::Format("Sh%04i_E: Ileak =%10.2f nA, Af =%10.2f#times10^{-13} V^{-2}, Rs =%10.2f #Omega ",detnum,Ie,Afe,Rse),"lep");
    /*legend->AddEntry(gf, TString::Format("Sh%04i_F: Ileak =%10.2f nA, Af =%10.2f#times10^{-13} V^{-2}, Rs =%10.2f #Omega ",detnum,If,Aff,Rsf),"lep");
    legend->AddEntry(gg, TString::Format("Sh%04i_G: Ileak =%10.2f nA, Af =%10.2f#times10^{-13} V^{-2}, Rs =%10.2f #Omega ",detnum,Ig,Afg,Rsg),"lep");
    legend->AddEntry(gh, TString::Format("Sh%04i_H: Ileak =%10.2f nA, Af =%10.2f#times10^{-13} V^{-2}, Rs =%10.2f #Omega ",detnum,Ih,Afh,Rsh),"lep");*/

    legend->Draw();

    TLatex latex;
    myc->SaveAs(TString::Format("det%i-run%i-%iC-nmfit.pdf",detnum,run,temp));
    myc->SaveAs(TString::Format("det%i-run%i-%iC-nmfit.png",detnum,run,temp));
    //return 0;
}
