#include <cmath>
#include "flux.hpp"
#include "constants.hpp"
#include "xSection.hpp"

Nucleus::Nucleus()
{
    A = 208;
    Z = 82;
    Rp = 5.5012;
    Rpn = 0.283;
    BEperA = 7.86745;
    abundance = 0.524;
    sp = 0.91;
    sn = 1.02;
    Mass = (A-Z)*MassN + Z*MassP - A*BEperA*1e-3;
}

Nucleus::Nucleus(int aa, int zz, double rp, double rpn, double be, double abun)
{
    A = aa;
    Z = zz;
    Rp = rp;
    Rpn = rpn;
    BEperA = be;
    abundance = abun;

// Averaged surface thickness for Pb, in fm, see Ref. Huang(2022)
    sp = 0.91;
    sn = 1.02;
    Mass = (A-Z)*MassN + Z*MassP - A*BEperA*1e-3;
}

// Calculate the form factor
double Nucleus::HelmFF(double q, double R, float s)
{
    double R0 = sqrt((R*R - 3*s*s)*5.0/3.0);
    double x = q*R0/HbarC;
    double j1 = sin(x)/(x*x) - cos(x)/x;

    return (3*j1/x)*exp(-(q*s/HbarC)*(q*s/HbarC)/2);
}

// Calculate the weak charge
double Nucleus::Gweak(double q)
{return gvp*Z*HelmFF(q,Rp,sp) + gvn*(A-Z)*HelmFF(q,Rp+Rpn,sn);}

/*----------------------------------------------------------------*/

// Calculating the cross section
double dxSection(double Enu, double Er, Nucleus nucleus)
{
    double cons = 1e-38*pow(HbarC*GF,2)/(2*pi);
    double M = nucleus.getMass();
    double qq = sqrt(2*M*Er);
    double GW = nucleus.Gweak(qq);

    return cons*M*(1 - M*Er/(Enu*Enu) + pow((1-Er*1e-3/Enu),2))*GW*GW;
}
