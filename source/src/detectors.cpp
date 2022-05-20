#include "constants.hpp"
#include "xSection.hpp"
#include "flux.hpp"
#include "detectors.hpp"
#include <gsl/gsl_integration.h>


// the main approach to adjust the speed of computation
// larger step & smaller nstep -> higher speed but poorer precision
const float Detectors::step = 0.1;  // (0.25) step for the integration on Enu,
const int Detectors::nstep = 20;    // (25) number of steps for the integration on Er

Detectors::Detectors()
{
    Nucleus Pb208;
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

    mass = 465;
    molarmass = 0.2072;
    nucleus = Pb208;
    flux = f;
    Natom = Pb208.abundance*mass*1e3*NA/molarmass;
    Enumax = 150.0;
}

Detectors::Detectors(float m, float molarm, Nucleus N, fluxpara f)
{
    mass = m;
    molarmass = molarm;
    nucleus = N;
    flux = f;
    Natom = N.abundance*mass*1e3*NA/molarmass;
    Enumax = 150.0;
}

// Differential scattering rate (in events/keV)

double Detectors::dRate(double Er)
{
    // integrate Natom*dxSection(Enu,Er,nucleus)*dFluence(Enu,flux) in [Enumin, Enumax]
    double x,sum,H,L;
    int NN;

    // integration range
    H = Enumax;
    L = nucleus.getEmin(Er);
    NN = floor((H-L)/step);

    x = L + 0.5*step;
    sum = 0;
    for(int i=0; i<NN; i++)
    {
        sum += dxSection(x,Er,nucleus)*dFluence(x,flux)*step;  // integrand
        x += step;
    }
    return Natom*sum;
}

// Differential scattering rate with GSL (in events/keV)
struct IntGSL_par
{
    double Er;
    Nucleus nucleus;
    fluxpara flux;
};

double IntegrandGSL(double Enu, void *params)
{
    IntGSL_par *mp = (IntGSL_par *)params;
    return dxSection(Enu,mp->Er,mp->nucleus)*dFluence(Enu,mp->flux);
}

double dRateGSL(double Er, void *params)
{
    Detectors *dec = (Detectors *)params;
    IntGSL_par mp = {Er,dec->nucleus,dec->flux};

    gsl_function f;
    f.function = IntegrandGSL;
    f.params = &mp;

    double result, err;
    size_t n;
    // integration range
    double H = dec->getEnumax();
    double L = dec->nucleus.getEmin(Er);

    gsl_integration_qng(&f,L,H,1e-32,1e-4,&result,&err,&n);
    //results for resnova4: (1,30)->(1e-27,1e-30)

    return result;
}

// Counts in the energy bin [Emin, Emax] (keV)
double Detectors::counts(float Emin, float Emax)
{
    gsl_function f;
    f.function = dRateGSL;
    f.params = this;

    double result, err;
    size_t n;

    gsl_integration_qng(&f,Emin,Emax,1e-32,1e-4,&result,&err,&n);
    //results for resnova4: 1keV bins in(1,30)->(1e-27,1e-30)

    return result*Natom;
}

/*
double Detector::counts(float Emin, float Emax)
{
    // integrate dRate(Er) in[Emin, Emax]
    double x,dx,sum;
    double H,L;

    // integration range
    H = Emax;
    L = Emin;

    dx = (H-L)/nstep;
    x = L + 0.5*dx;
    sum = 0;
    for(int i=0; i<nstep; i++)
    {
        sum += dRate(x)*dx; // integrand
        x += dx;
    }

    return sum/Natom;
}*/

/*----------------------------------------------------------------*/
// Detector objects

// Res-Nova
double ResNovaCounts(float Emin, float Emax, fluxpara flux)
{
    float mass = 465;
    float molmass = 0.2072;

    Nucleus Pb204(204,82,5.4803,0.283,7.87993,0.014);
    Nucleus Pb206(206,82,5.4902,0.283,7.87536,0.241);
    Nucleus Pb207(207,82,5.4943,0.283,7.86987,0.221);
    Nucleus Pb208(208,82,5.5012,0.283,7.86745,0.524);

    Detectors resnova1(mass,molmass,Pb204,flux);
    Detectors resnova2(mass,molmass,Pb206,flux);
    Detectors resnova3(mass,molmass,Pb207,flux);
    Detectors resnova4(mass,molmass,Pb208,flux);

    return resnova1.counts(Emin,Emax)+resnova2.counts(Emin,Emax)
           +resnova3.counts(Emin,Emax)+resnova4.counts(Emin,Emax);
}

// Argo
double ArgoCounts(float Emin, float Emax, fluxpara flux)
{
    float mass = 300;
    float molmass = 0.03995;

    Nucleus Ar40(40,18,3.4274,0.1,8.5953,1.0);

    Detectors argo(mass,molmass,Ar40,flux);

    return argo.counts(Emin,Emax);
}

// Darwin
double DarwinCounts(float Emin, float Emax, fluxpara flux)
{
    float mass = 50;
    float molmass = 0.1313;

    // Data for Xe, to be updated
    Nucleus Xe132(132,54,4.7859,0.2,8.4275,1.0);

    Detectors darwin(mass,molmass,Xe132,flux);

    return darwin.counts(Emin,Emax);
}

