#include "pysnow.h"
#include "pySNvC.h"
#include <cmath>

using std::log;

extern "C" {
    double log_likelihood(double x[]);
}

double log_prior(double x[], size_t x_size=9) {
    static double default_paras[] = {2.5, 2.5, 2.5,
                                     9.5, 12., 15.6, 
                                     5.,  5.,  5.};
    static double ratio[] = {.33333333, 2.};
    for(size_t i = 0; i != x_size; ++i) {
        if(x[i]<default_paras[i]*ratio[0] || x[i]>default_paras[i]*ratio[1]) {
            return -1e150;
        }
    }
    return 0.;
}


static double lambHK[NBinsHK]; // output rates per bin
static double lambDUNE[NBinsDUNE]; 
static double lambRes[NBinsRes];
double log_likelihood(double x[]) {

    static const double N_HK[NBinsHK] = {0.};
    static const double N_DUNE[NBinsDUNE] = {0.};
    static const double N_RESNOV[NBinsRes] = {0.};
    
    size_t x_size = 9;
    double res = log_prior(x);
    if(res < -1e149) return -1e150;

    for(size_t i = 0; i != x_size; ++i) {
        if(x[i] < 0.) return -1e150;
    }
    
    // HyperK
    rateGen(expr_hk, x, 10., lambHK, NBinsHK);
    for(size_t i = 0; i != NBinsHK-1; ++i) {
        res += -lambHK[i] + N_HK[i] * log(lambHK[i]);
    }

    // DUNE
    rateGen(expr_dune, x, 10., lambDUNE, NBinsDUNE);
    for(size_t i = 0; i != NBinsDUNE; ++i) {
        res += -lambDUNE[i] + N_DUNE[i] * log(lambDUNE[i]);
    }

    // ResNova
    pygetSpec(x, lambRes); // 10kpc by default
    for(size_t i = 0; i != NBinsRes; ++i) {
        res += -lambRes[i] + N_RESNOV[i] * log(lambRes[i]);
    }

    return res;
}
