//Output of this

#include <iostream>
#include <Math/RootFinderAlgorithms.h>
#include <TF1.h>
#include <Math/RootFinder.h>
#include <Math/Functor.h>
using namespace ROOT::Math;
using namespace std;


//Hmm so this seems to work lol. How SLOOWWWWW will it be when I implement the complicated version of this?
float Edep = 0.77; //0.679 (beta = 0.965) is the lower limit right now for the silicon!
//Let's just go to 0.95 for now. Close to the minimum, which is ~0.96. Then you get redundant solutions which I don't want to deal with right this second.

float z = 1;
float Zeff = 14;
float Aeff = 28;
float rho = 2.33; //Density g/cm^2
float L = 0.25;

float ion = (0.000016 * pow(Zeff,0.9));
float C_1 = 0.3071/2; //MeV/ g/cm^2 #2*pi*constants not 4*pi*constants for MPV
float me = 0.511; //mass of electron * c^2


//MPV = (C1*z**2 *(rho*L*Zeff/Aeff)*(1/beta**2))*(math.log( (1.022 * gamma**2 * beta**2)/ion ) + math.log( (C1*z**2 *(rho*L*Zeff/Aeff)* (1/beta**2)) /ion) + 0.2 - beta**2)

double myfunc(double x)
{
    //float toybeta = 0.9;

    //cout << "toy MPV " << C_1*pow(z,2)*(rho*L*Zeff/Aeff)*(1/pow(toybeta,2))* ( log( (1.022 * pow(toybeta,2)/(1 - pow(toybeta,2)) )/ion ) + log(C_1*pow(z,2)*(rho*L*Zeff/Aeff)*(1/pow(toybeta,2))/ion)  + 0.2 - pow(toybeta,2) ) << endl;

    return -Edep + C_1*pow(z,2)*(rho*L*Zeff/Aeff)*(1/pow(x,2))* ( log( (1.022 * pow(x,2)/(1 - pow(x,2)) )/ion ) + log(C_1*pow(z,2)*(rho*L*Zeff/Aeff)*(1/pow(x,2))/ion)  + 0.2 - pow(x,2) );
    //return -Edep + log(x * sqrt(1/x));
}

void find_root()
{
    RootFinder *k = new RootFinder();
    k->SetMethod(RootFinder::kGSL_BISECTION);
    ROOT::Math::Functor1D f(&myfunc);
    //ROOT::Math::Functor1D f(&myfunc2);
    k->SetFunction(f, 0.001, 0.96); //Surely this is the range of values that get tested?
    //How to make it so I can change the MPV in the function without having to re-initialize the function or anything?
    //So maybe let's say only use this technique up to like 0.94 to make this cool and good nice
    //YEAHHHHHH I THINK I GOT IT!!!!! WOOOOOO!!!!!!
    k->Solve();

    double c = k->Root();
    cout << "Edep is " << Edep << "Calculated MPV is " <<  c << endl;

}



double myfunc2(double x)
{
    return C_1*pow(z,2)*(rho*L*Zeff/Aeff)*(1/pow(x,2))* ( log( (1.022 * pow(x,2)/(1 - pow(x,2)) )/ion ) + log(C_1*pow(z,2)*(rho*L*Zeff/Aeff)*(1/pow(x,2))/ion)  + 0.2 - pow(x,2) );
}


//YAYYY!!
void find_root_2()
{
    auto fa3 = new TF1("fa3","myfunc2(x)",0.01,0.95);
    double root = fa3->GetX(0.69);
    cout << "Edep PLEASE " << root << endl;
}

double ZOne_Tkr(double x)
{
    return C_1*pow(z,2)*(rho*L*Zeff/Aeff)*(1/pow(x,2))* ( log( (1.022 * pow(x,2)/(1 - pow(x,2)) )/ion ) + log(C_1*pow(z,2)*(rho*L*Zeff/Aeff)*(1/pow(x,2))/ion)  + 0.2 - pow(x,2) );
}

//Try figuring out why this worked: https://root.cern/doc/v620/classTF1.html
TF1 *Z1_tkr_Solve = new TF1("Z1_tkr_Solve", [](double *x, double *p){ return ZOne_Tkr(x[0]); }, 0.01, 0.95, 0);
double solve_beta(double mpv){
    double root = Z1_tkr_Solve->GetX(mpv);
    return root;
}
