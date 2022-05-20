/* ---------------------------------------------
   Calculating the differential neutrino fluence on the Earth 
   for a given neutrino energy
 
   function:     double dFluence(double ene, fluxpara para)
   parameters:   E_\nu in MeV, 9 spectrum parameters, distance in kpc

   return:       \Phi(E_\nu) for all flavor, in MeV^-1 cm^-2
   -----------------------------------------------*/

#include <cmath>
//#include "constants.hpp"

//double dFluence(double ene, fluxpara para);

/*----------------------------------------------------------------*/
// Define the data struct type for the flux parameters
typedef struct fluxpara{

/* data struct for each flavor: {Etot,Eave,alp}
   ----------------------------
   Etot: total energy in 10^52 erg;
   Eave: averaged v energy in MeV;
   alp: spectral pinching;

   distance in kpc;
   ----------------------------*/
    float ve[3];
    float vae[3];
    float vx[3];   //the averaged luminosity for one flavor
    float distance;
}fluxpara;

/*----------------------------------------------------------------*/
// Define the f(E_\nu) function
double spectrum(double ene, float Eave, float alp){
    double AA = pow(alp+1,alp+1)/(Eave*tgamma(alp+1));
    return AA * pow(ene/Eave, alp) * exp(-(alp+1)*ene/Eave);
}

double dFluence(double ene, fluxpara para){
    /*-----------------------------------------------------------
    parameters
    ----------------------
    ene : neutrino energy in MeV
    para : spectral parameters in the form of fluxpara

    returns
    --------------------
    double     total differential fluence on the Earth (in MeV^-1 cm^-2)
    ----------------------------------------------------------------*/
    double dd = para.distance * kpctocm;
    double Lve = para.ve[0]*spectrum(ene,para.ve[1],para.ve[2])/para.ve[1];
    double Lvae = para.vae[0]*spectrum(ene,para.vae[1],para.vae[2])/para.vae[1];
    double Lvx = para.vx[0]*spectrum(ene,para.vx[1],para.vx[2])/para.vx[1];
    return (Lve+Lvae+4.0*Lvx)*1.0e52/(4*pi*dd*dd*MeVtoerg);
}

/*----------------------------------------------------------------*/