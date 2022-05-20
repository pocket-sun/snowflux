/* ---------------------------------------------
   Calculating the event spectrum in the energy range [Emin,Emax] for a given neutrino flux
 
   function:     void getSpec(double (*f)(float,float,fluxpara), float (*spec)[2], float Eth, float binStep, int nbin, fluxpara flux)
   parameters:   count function, energy threshold (keV), bin step (keV), number of bins, neutrino flux
   count function: ResNovaCounts, ArgoCounts, DarwinCounts 

   return:       spec[nbin][2], {Er, N} for each bin with the central energy Er (keV)
   -----------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <time.h>
#include "detectors.hpp"
using namespace std;
using std::setw;

extern "C" {
void pygetSpec(double inputs[], double outputs[]);
}
void getSpec(double (*fun)(float,float,fluxpara), float (*spec)[2], float Eth, float binStep, int nbin, fluxpara flux);

int main()
{
    // spectral parameters
    float Eth = 0.5;
    float binStep = 1;
    int nbin = 20;

    // neutrino flux for tests
    fluxpara f;        // data from [Horiuchi et al, PRD 2018]
    f.distance = 10;
    f.ve[0] = 7.7;
    f.vae[0] = 6.44;
    f.vx[0] = 5.88;
    f.ve[1] = 14.1;
    f.vae[1] = 16.3;
    f.vx[1] = 17.2;
    f.ve[2] = 2.67;
    f.vae[2] = 3.28;
    f.vx[2] = 2.20;
    
    fluxpara f0;         // largest fluence
    f0.distance = 10;
    f0.ve[0] = 10;
    f0.vae[0] = 10;
    f0.vx[0] = 10;
    f0.ve[1] = 30;
    f0.vae[1] = 30;
    f0.vx[1] = 30;
    f0.ve[2] = 3.5;
    f0.vae[2] = 3.5;
    f0.vx[2] = 3.5;

    fluxpara f1;         // smallest fluence
    f1.distance = 10;
    f1.ve[0] = 2;
    f1.vae[0] = 2;
    f1.vx[0] = 2;
    f1.ve[1] = 5;
    f1.vae[1] = 5;
    f1.vx[1] = 5;
    f1.ve[2] = 1.5;
    f1.vae[2] = 1.5;
    f1.vx[2] = 1.5;

    // compute the spectrum for ResNova
    float spec[50][2];

    // record the run time
    clock_t start,end;
    start = clock();

    getSpec(ResNovaCounts,spec,Eth,binStep,nbin,f);

    end = clock();
    cout<<"time = "<<double(end-start)/CLOCKS_PER_SEC<<"s"<<endl; // (s)

    // output the results
    cout << setw(10) << "Er (keV)" << setw(15) << "Counts" << endl;
    for (int i = 0; i < nbin; i++)
    {
        cout << setw(10) << spec[i][0] << setw(15) << spec[i][1] << endl;
    }

    return 0;
}

/*------------------------------------------------*/
// calculate the expected spectrum
void getSpec(double (*fun)(float,float,fluxpara), float (*spec)[2], float Eth, float binStep, int nbin, fluxpara flux)
{
    float x=0, L=0, H=0;
    for (int i = 0; i < nbin; i++)
    {
        x = Eth + (i+0.5)*binStep;
        L = x - 0.5*binStep;
        H = x + 0.5*binStep;
        spec[i][0] = x;
        spec[i][1] = fun(L,H,flux);
    }
}
/*------------------------------------------------*/
// interface to ctypes
#define pyEth 0.5
#define pybinStep 1
#define pynbin 20
static float pyspec[50][2];
// length(inputs)=9, length(outputs)=pynbin
// which should be guarteed by user and no safty check is warranted
void pygetSpec(double inputs[], double outputs[]) {

    // spectral parameters
//    static const float Eth = 0.5;
//    static const float binStep = 1;
//    static const int nbin = 20;

    // neutrino flux for tests
    static fluxpara f;     
    f.distance = 10;
    for(size_t k = 0; k != 3; ++k) {
        f.ve[k] = inputs[3*k];
        f.vae[k] = inputs[3*k+1];
        f.vx[k] = inputs[3*k+2];
    }
    getSpec(ResNovaCounts,pyspec,pyEth,pybinStep,pynbin,f);
    for(size_t k = 0; k != pynbin; ++k) {
        outputs[k] = pyspec[k][1];
    }
}
