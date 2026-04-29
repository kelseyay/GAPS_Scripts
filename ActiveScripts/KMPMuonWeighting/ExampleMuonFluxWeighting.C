//Commenting this right now to make sure I understand it!

#include <TFile.h>
#include <TStyle.h>
#include <TTree.h>
#include <TROOT.h>
#include <TMath.h>
#include <TChain.h>
#include <vector>
#include <string>

#include "CAnalysisManager.hh"
#include "GSimulationParameter.hh"
#include "GPlottingTools.hh"

#include "GGeometry.hh"
#include <CEventRec.hh>

using namespace ROOT::Math;

using namespace std;
using namespace Crane::Analysis;
namespace ca = Crane::Analysis;

// example for making a weighted histogram of beta values for simulation v300

void ExampleWeighting(){

    /////////////////////////////////////////
    // v.3.0.0 Simulation
    /////////////////////////////////////////

    cout<<" Simulation v.3.0.0"<<endl;

    //step size to sample over the files, helpful for testing to not run over the full files, which is potentially slow
    int MainLoopScaleFactor = 1;

    string Directory;
    Directory = "/home/kmp2215/scripts/";

    char placeholder[400];
    sprintf(placeholder,"%s/%s","the code segfaults without this but completes the script", "unclear why");

    CAnalysisManager AnalysisManagerRec;


    //--------------------------------------

    //This object is managing the plotting
    ca::GPlottingTools Plotting;

    int xlim = 2;
    int xbin = 100;

    TH1D * HBeta = Plotting.DefineTH1D("HBeta", xbin, 0, xlim, "Beta", "Muon Rate", 0.5, 1e3);
    TH1D * HCosTheta = Plotting.DefineTH1D("HCosTheta", xbin, 0, xlim, "cos(theta)", "Muon Rate", 0.5, 100);


    /////////////////////////////////////////
    // Setting up muon weighting
    /////////////////////////////////////////

    double FluxScaleFactor = 71.1552/32.058750*0.25;

    vector<pair<double, double> > CosZenithCut;
    CosZenithCut.push_back(make_pair(-0.75, -1));
    CosZenithCut.push_back(make_pair(-0.5, -0.75));
    CosZenithCut.push_back(make_pair(-0.25, -0.5));
    CosZenithCut.push_back(make_pair(0, -0.25));

    vector<TGraph*> GMuonTotalFluxUnscaled;
    GMuonTotalFluxUnscaled.push_back(Plotting.ConvertEnergyFluxToBetaFlux(Plotting.GetTH1D(getenv("GAPS") + string("/resources/fluxes/total_fluxes_coszenith_100_m_-13_antarctica.root"), "c6a", "p_total_altitude_zenith_energy_0.875"), 0.1057));
    GMuonTotalFluxUnscaled.push_back(Plotting.ConvertEnergyFluxToBetaFlux(Plotting.GetTH1D(getenv("GAPS") + string("/resources/fluxes/total_fluxes_coszenith_100_m_-13_antarctica.root"), "c6a", "p_total_altitude_zenith_energy_0.625"), 0.1057));
    GMuonTotalFluxUnscaled.push_back(Plotting.ConvertEnergyFluxToBetaFlux(Plotting.GetTH1D(getenv("GAPS") + string("/resources/fluxes/total_fluxes_coszenith_100_m_-13_antarctica.root"), "c6a", "p_total_altitude_zenith_energy_0.375"), 0.1057));
    GMuonTotalFluxUnscaled.push_back(Plotting.ConvertEnergyFluxToBetaFlux(Plotting.GetTH1D(getenv("GAPS") + string("/resources/fluxes/total_fluxes_coszenith_100_m_-13_antarctica.root"), "c6a", "p_total_altitude_zenith_energy_0.125"), 0.1057));

    //currently not working, will be updated for newer simulations
    //std::string recoName = "FindHough3D";
    //-----------------------------------------------------------------------

    char FilenameRoot[400];
    sprintf(FilenameRoot,"%s/%s.root","/data1/nextcloud/cra_data/nobackup/MC/v3.0.0/full/mu-/", "mu-_gaps_triggerlevel1_FTFP_BERT_1757442*");

    //prepare reconstronstructed event
    CEventRec* Event = new CEventRec;
    TChain * TreeRec = new TChain("TreeRec");
    TreeRec->SetBranchAddress("Rec", &Event);
    TreeRec->Add(FilenameRoot);
    TreeRec->GetEntry(0);
    Event->SetEventTime(double(Event->GetEventTime())/(1000./64.)+1631030675);//placeholder for fc conversion to unix time

    //These object are managing the analysis functionality
    AnalysisManagerRec.SetEvent(Event);

    //------------------------------------------------------------------------
    //plots are shown as a function of beta
    //this is the number of bins
    int BetaBins = 25;
    double StartingPlaneAcceptance = 1;
    TH1D* HPrimaryBeta = nullptr;
    double BinWidthFactor = 1;
    std::vector<double> PrimaryBetaLowHigh;

    TChain*TreeSimulationParameter = new TChain("SimulationParameterTree");
    TreeSimulationParameter->Add(FilenameRoot);
    GSimulationParameter * Parameter = new GSimulationParameter;
    TreeSimulationParameter->SetBranchAddress("SimulationParameter", &Parameter);
    TreeSimulationParameter->GetEntry(0);

    AnalysisManagerRec.SetGSimulationParameterTChain(TreeSimulationParameter);
    StartingPlaneAcceptance = AnalysisManagerRec.GetStartingPlaneAcceptance();

    //Find low and high range of beta from simulation parameters
    PrimaryBetaLowHigh = AnalysisManagerRec.GetPrimaryBetaLowHigh();

    HPrimaryBeta = AnalysisManagerRec.GetHPrimaryBeta();
    //HPrimaryBeta = new TH1D("HPrimaryBeta", "", BetaBins, PrimaryBetaLowHigh.at(0), PrimaryBetaLowHigh.at(1));

    BinWidthFactor = (PrimaryBetaLowHigh.at(1)-PrimaryBetaLowHigh.at(0))/double(BetaBins) / HPrimaryBeta->GetBinWidth(1);

    //------------------------------------

    int counts = 0;

    //loop over simulated events
    for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor)
    {
            //------------------------------------
            AnalysisManagerRec.EventReset();
            TreeRec->GetEntry(i);

            //choose a reconstruction for beta, currently contains bugs, v300 only contains FindPrimaryStarIterative so this remains commented out
            //if (AnalysisManagerRec.IsRec())
            //        {
            //        AnalysisManagerRec.GetCEventRec()->ChooseReconstruction(recoName);
            //        AnalysisManagerRec.GetCEventRec()->ListAvailableReconstructions();
            //        }

            double AcceptanceScale;
            if (TreeSimulationParameter != nullptr)
                    {
                    //acceptance scaling factor based on beta of the primary
                    AcceptanceScale = MainLoopScaleFactor*StartingPlaneAcceptance/(BinWidthFactor*HPrimaryBeta->GetBinContent(HPrimaryBeta->FindBin(Event->GetPrimaryBetaGenerated())));

                    if(HPrimaryBeta->GetBinContent(HPrimaryBeta->FindBin(Event->GetPrimaryBetaGenerated())) == 0) AcceptanceScale = 0;
                    }
            else AcceptanceScale = 1;

	    int AngularRegion = -1;
            for(unsigned int a = 0; a < CosZenithCut.size(); a++) if(Event->GetPrimaryMomentumDirectionGenerated().CosTheta() < CosZenithCut.at(a).first && Event->GetPrimaryMomentumDirectionGenerated().CosTheta() > CosZenithCut.at(a).second) AngularRegion = a;
	    if(AngularRegion < 0) continue;
            double RateScale = FluxScaleFactor*AcceptanceScale*GMuonTotalFluxUnscaled.at(AngularRegion)->Eval(Event->GetPrimaryBetaGenerated());


            //------------------------------------
	    double BetaTruth = Event->GetPrimaryBetaGenerated();
	    double costheta = Event->GetPrimaryMomentumDirectionGenerated().CosTheta();

	    HBeta->Fill(BetaTruth,RateScale);
	    HCosTheta->Fill(-costheta,RateScale);
	    counts++;

    }

    cout<<"number of events in v.3.0.0"<<endl;
    cout<<counts<<endl;

    gROOT->Reset();
    TStyle * plain = new TStyle("plain","plain");
    plain->SetCanvasBorderMode(0);
    plain->SetPadBorderMode(0);
    plain->SetPadColor(0);
    plain->SetCanvasColor(0);
    plain->SetTitleColor(1);
    plain->SetStatColor(0);
    plain->SetTitleFillColor(0);
    plain->SetLineWidth(2);
    plain->SetHistLineWidth(4);

    gROOT->SetStyle("plain");

    TCanvas * CBeta = new TCanvas("CBeta", "CBeta", 200, 10, 900, 900);
    CBeta->SetLeftMargin(0.11);
    CBeta->SetRightMargin(0.04);
    CBeta->SetTopMargin(0.04);
    HBeta->GetXaxis()->SetRangeUser(0, 1.25);

    HBeta->SetLineColor(6);
    HBeta->Draw("hist");


    gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->SetLogy(1);

    char text[400];
    sprintf(text, "%sweightedbetav300.png", Directory.c_str());

    cout<<"muon rate: "<<HBeta->Integral(0,xbin)<<"Hz"<<endl;

    CBeta->SaveAs(text);

    TCanvas * CCosTheta = new TCanvas("CCosTheta", "CCosTheta", 200, 10, 900, 900);
    CCosTheta->SetLeftMargin(0.11);
    CCosTheta->SetRightMargin(0.04);
    CCosTheta->SetTopMargin(0.04);
    HCosTheta->GetXaxis()->SetRangeUser(0, 1.25);

    HCosTheta->SetLineColor(1);
    HCosTheta->Draw("hist");

    cout<<"muon rate: "<<HCosTheta->Integral(0,xbin)<<"Hz"<<endl;

    gPad->SetGridx(1);
    gPad->SetGridy(1);
    //gPad->SetLogy(1);

    sprintf(text, "%sweightedcosthetav300.png", Directory.c_str());
    CCosTheta->SaveAs(text);
}
