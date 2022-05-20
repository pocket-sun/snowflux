#ifndef FLUXHPP
#define FLUXHPP
/* ---------------------------------------------
   Calculating the differential neutrino fluence on the Earth 
   for a given neutrino energy
 
   function:     double dFluence(double ene, fluxpara para)
   parameters:   E_\nu in MeV, 9 spectrum parameters, distance in kpc

   return:       \Phi(E_\nu) for all flavor, in MeV^-1 cm^-2
   -----------------------------------------------*/

#include <cmath>
//#include "constants.hpp"


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

double spectrum(double ene, float Eave, float alp);
double dFluence(double ene, fluxpara para);

#endif
