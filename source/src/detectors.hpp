#ifndef DETECTORSHPP
#define DETECTORSHPP
/* ---------------------------------------------
   Calculating the event counts in the energy bin [Emin,Emax] for a given neutrino flux
 
   function:     double ResNovaCounts(float Emin, float Emax, fluxpara flux)
                 double ArgoCounts(float Emin, float Emax, fluxpara flux)
                 double DarwinCounts(float Emin, float Emax, fluxpara flux)
   parameters:   recoil energy bin [Emin,Emax] (keV), neutrino flux
   
   return:       N_i
   -----------------------------------------------*/

/* ---------------------------------------------
   Define a detector
   -----------------------------------------------*/
#include "constants.hpp"
#include "xSection.hpp"
#include "flux.hpp"

double ResNovaCounts(float Emin, float Emax, fluxpara flux);
double ArgoCounts(float Emin, float Emax, fluxpara flux);
double DarwinCounts(float Emin, float Emax, fluxpara flux);
/*-------------------------------------------------------*/
// Define integration functions
double IntegrandGSL(double Enu, void *params);
double dRateGSL(double Er, void *params);
/* --------------------------------------------------------*/

/* --------------------------------------------------------*/
/* Define detectors. 
   The default construction function will create the Pb208-part of Res-Nova, 
   and a SN at 10 kpc. */
class Detectors
{
private:
    double Natom;                 // number of the target nucleus
    float Enumax;                 // top limit for the neutrino energy in MeV
    static const float step;            // step for the integration on Enu
    static const int nstep;             // number of steps for the integration on Er
public:
    float mass;          // detector mass (ton)
    float molarmass;     // molarmass for the detector material (kg/mol)
    Nucleus nucleus;     // target nucleus
    fluxpara flux;       // neutrino fluence
    Detectors();
    Detectors(float m, float molarm, Nucleus N, fluxpara f);
    ~Detectors(){}

    /* Reduce the top limit for the neutrino energy and nstep if you try to accelerate it
       Caution: you may lose the information from the high energy tail of neutrino spectrum !
       The default value is 150.0 MeV */
    void setEnumax(float Enu){Enumax = Enu;}
    float getEnumax(){return Enumax;}

    // Calculating the scattering rate and counts
    double dRate(double Er);
    double counts(float Emin, float Emax);
};
/*----------------------------------------------------------------*/
#endif
