#ifndef XSECHPP
#define XSECHPP
/* ------------------------------------------------
   Calculating the differential cross section for the coherent scattering

   function:    double dxSection(double Enu, double Er, Nucleus nucleus)
   parameters:  neutrino energy in MeV, recoil energy in keV, nucleus

   return:      differential cross section in cm^2/keV
   ------------------------------------------------*/

#include <cmath>
//#include "constants.hpp"

class Nucleus;
double dxSection(double Enu, double Er, Nucleus nucleus);

/* ------------------------------------------------
   Define a nucleus.
   ------------------------------------------------*/

// Define the nucleus, the default construction function will create a Pb208.
class Nucleus
{
private:
    float sp, sn;                 // surface thickness in fm
    double Mass;                  // the nucleus mass in GeV
public:
    // Parameters
    int A;                        // number of nucleons
    int Z;                        // number of protons
    double Rp;                    // rms radius of proton in fm
    double Rpn;                   // neutron skin in fm
    double BEperA;                // binding energy per nucleon in MeV
    double abundance;             // element abundance

    // Basic functions
    Nucleus();
    Nucleus(int aa, int zz, double rp, double rpn, double be, double abun);
    ~Nucleus(){}

    // Accesses to the surface thickness
    void setSp(float Sp){sp = Sp;}
    float getSp(){return sp;}
    void setSn(float Sn){sn = Sn;}
    float getSn(){return sn;}

    // Get the nucleus mass in GeV
    double getMass(){return Mass;}

    // Form factor (Helm) & weak charge
    double HelmFF(double q, double R, float s);
    double Gweak(double q);

    // Get the minimun neutrino energy (MeV) for a given recoil energy (keV).
    double getEmin(double Er){return sqrt(getMass()*Er/2);}
};

/*----------------------------------------------------------------*/
#endif
