/* ------------------------------------------------
   Define constants.
   ------------------------------------------------*/

const double HbarC = 197.327;         // hc/2pi in MeV fm
const double MeVtoerg = 1.602177e-6;           // MeV -> erg
const double kpctocm  = 3.0857e21;             // kpc -> cm

const double MassP = 0.938272;        // proton mass in GeV
const double MassN = 0.939565;        // neutron mass in GeV
const double NA = 6.02214e23;         // Avogadro constant in mol-1

const double GF = 1.16638e-5;         // Fermi constant in GeV-2
const double wma = 0.23857;           // weak mixing angle at low energy
const double gvp = 0.5 - 2*wma;       // coupling constant for weak interaction
const float gvn = -0.5;

const double pi = 3.1415926;

/*----------------------------------------------------------------*/

/*----------------------------------------------------------------
double integral(double (*fun)(double), double L, double H, long n = 10000);
double integral(double (*fun)(double), double L, double H, long n)
{
    double x,dx,sum;
    dx = (H-L)/n;
    x = L + 0.5*dx;
    sum = 0;
    for(int i=0; i<n; i++)
    {
        sum += fun(x)*dx;
        x += dx;
    }
    return sum;
}
----------------------------------------------------------------*/