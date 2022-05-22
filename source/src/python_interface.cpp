#include "pysnow.h"
#include "pySNvC.h"
#include <cmath>

using std::log;

extern "C" {
    double log_likelihood(double x[]);
}

#ifdef DEBUG
void myprt(const double *a, const char* info, size_t asz=10) {
    printf("%s:\n", info);
    for(size_t k = 0; k != asz; ++k) {
        printf("%7.2lf  ", a[k]);
    }
    printf("\n");
}
#endif

double log_prior(double x[], size_t x_size=9) {
    static double default_paras[] = {2.5, 2.5, 2.5,
                                     9.5, 12., 15.6, 
                                     5.,  5.,  5.};
    static double ratio[] = {.33333333, 2.};
    for(size_t i = 0; i != 3; ++i) {
        if(x[i]<1.5 || x[i]>3.5) 
            return -1e150;
    }
    for(size_t i = 3; i != x_size; ++i) {
        if(x[i]<default_paras[i]*ratio[0] || x[i]>default_paras[i]*ratio[1]) {
            return -1e150;
        }
    }
    return 0.;
}


//static double lambHK[NBinsHK]; // output rates per bin
static double lambHKibd[NBinsHKibd];
static double lambHKes[NBinsHKes]; 
static double lambDUNE[NBinsDUNE]; 
static double lambRes[NBinsRes];
double log_likelihood(double x[]) {

    /* data input */
//  static const double N_HK[NBinsHK] = {0.};
    static const double N_HKibd[NBinsHKibd] = {1170,3632,4314,3428,2250,1354,693,347,160,91};
    static const double N_HKes[NBinsHKes] = {164,214,133,90,45,19,9,6,2,2};
    static const double N_DUNE[NBinsDUNE] = {341,496,390,211,110,42,19,9,5,1};
    static const double N_RESNOV[NBinsRes] = {6318,2780,1457,831,466,310,205,119,70,43};
    
    static const size_t x_size = 9;
    double res = log_prior(x);
    if(res < -1e149) return -1e150;

    for(size_t i = 0; i != x_size; ++i) {
        if(x[i] < 0.) return -1e150;
    }
    
    // HyperK
//    rateGen(expr_hk, x, 10., lambHK, NBinsHK);
//    for(size_t i = 0; i != NBinsHK-1; ++i) {
//        res += -lambHK[i] + N_HK[i] * log(lambHK[i]);
//    }
    
    rateGen(expr_hkibd, x, 10., lambHKibd, NBinsHKibd);
    for(size_t i = 0; i != NBinsHKibd-1; ++i) {
        res += -lambHKibd[i] + N_HKibd[i] * log(lambHKibd[i]);
    }
    rateGen(expr_hkes, x, 10., lambHKes, NBinsHKes);
    for(size_t i = 0; i != NBinsHKes-1; ++i) {
        res += -lambHKes[i] + N_HKes[i] * log(lambHKes[i]);
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

#ifdef DEBUG
    myprt(lambHKibd, "lambHKibd");
    myprt(N_HKibd, "N_HKibd");
    myprt(lambHKes, "lambHKes");
    myprt(N_HKes, "N_HKes");
    myprt(lambDUNE, "lambDUNE");
    myprt(N_DUNE, "N_DUNE");
    myprt(lambRes, "lambRes");
    myprt(N_RESNOV, "N_RESNOV");
#endif

#ifdef DEBUG
    printf("res=%lf\n", res);
#endif
    return res;
}

