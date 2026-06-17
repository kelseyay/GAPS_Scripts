//Output of this

#include <iostream>
#include <Math/RootFinderAlgorithms.h>
#include <TF1.h>
#include <Math/RootFinder.h>
#include <Math/Functor.h>
using namespace ROOT::Math;
using namespace std;


//Hmm so this seems to work lol. How SLOOWWWWW will it be when I implement the complicated version of this?
float Edep = 0.7372;

float z = 1;
float Zeff = 14;
float Aeff = 28;
float rho = 2.33; //Density g/cm^2
float L = 0.22;

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
    k->SetFunction(f, 0.001, 0.95); //Surely this is the range of values that get tested?
    //So maybe let's say only use this technique up to like 0.94 to make this cool and good nice
    //YEAHHHHHH I THINK I GOT IT!!!!! WOOOOOO!!!!!!
    k->Solve();

    double c = k->Root();
    cout << c << endl;
}
