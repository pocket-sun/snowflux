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
static double lambJUNOpes[NBinsJUNO];
static double lambJUNOes[NBinsJUNO];
static double lambJUNOibd[NBinsJUNO];
static double lambArgo[NBinsArgo];
double log_likelihood(double x[]) {


    /* data input */
//  static const double N_HK[NBinsHK] = {0.};
    static const double N_HKibd[NBinsHKibd] = {3698,8090,9497,7572,4970,2965,1543,767,357,192};
    static const double N_HKes[NBinsHKes] = {965,476,295,189,97,46,16,12,7,4};
    static const double N_DUNE[NBinsDUNE] = {493,555,413,221,115,43,19,9,5,1};
    static const double N_RESNOV[NBinsRes] = {6318,2780,1457,831,466,310,205,119,70,43};
    static const double N_JUNOpes[NBinsJUNO] = {83};
//    static const double N_JUNOes[NBinsJUNO] = {149};
//    static const double N_JUNOibd[NBinsJUNO] = {51};
    static const double N_ARGO[NBinsArgo] = {282,227,164,140,128,104,101,77,73,66};
    
    static const size_t x_size = 9;
    double res = log_prior(x);
    if(res < -1e149) return -1e150;

    for(size_t i = 0; i != x_size; ++i) {
        if(x[i] < 0.) return -1e150;
    }
    
    // HyperK
    /*
    rateGen(expr_hk, x, 10., lambHK, NBinsHK);
    for(size_t i = 0; i != NBinsHK-1; ++i) {
        res += -lambHK[i] + N_HK[i] * log(lambHK[i]);
    }
    */
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
    pygetSpecResNova(x, lambRes); // 10kpc by default
    for(size_t i = 0; i != NBinsRes; ++i) {
        res += -lambRes[i] + N_RESNOV[i] * log(lambRes[i]);
    }

    // JUNO
    rateGen(expr_junopes, x, 10., lambJUNOpes, NBinsJUNO);
    for(size_t i = 0; i != NBinsJUNO; ++i) {
        res += -lambJUNOpes[i] + N_JUNOpes[i] * log(lambJUNOpes[i]);
    }
    /*
    rateGen(expr_junoes, x, 10., lambJUNOes, NBinsJUNO);
    for(size_t i = 0; i != NBinsJUNO; ++i) {
        res += -lambJUNOes[i] + N_JUNOes[i] * log(lambJUNOes[i]);
    }
    rateGen(expr_junoibd, x, 10., lambJUNOibd, NBinsJUNO);
    for(size_t i = 0; i != NBinsJUNO; ++i) {
        res += -lambJUNOibd[i] + N_JUNOibd[i] * log(lambJUNOibd[i]);
    }
    */

    // Argo
    pygetSpecArgo(x, lambArgo); // 10kpc by default
    for(size_t i = 0; i != NBinsArgo; ++i) {
        res += -lambArgo[i] + N_ARGO[i] * log(lambArgo[i]);
    }
    
#ifdef DEBUG
    myprt(lambHKibd, "lambHKibd", NBinsHKibd);
    myprt(N_HKibd, "N_HKibd", NBinsHKibd);
    myprt(lambHKes, "lambHKes", NBinsHKes);
    myprt(N_HKes, "N_HKes", NBinsHKes);
    myprt(lambDUNE, "lambDUNE", NBinsDUNE);
    myprt(N_DUNE, "N_DUNE", NBinsDUNE);
    myprt(lambRes, "lambRes", NBinsRes);
    myprt(N_RESNOV, "N_RESNOV", NBinsRes);
    myprt(lambJUNOpes, "lambJUNOpes", NBinsJUNO);
    myprt(N_JUNOpes, "N_JUNOpes", NBinsJUNO);
//    myprt(lambJUNOes, "lambJUNOes", NBinsJUNO);
//    myprt(N_JUNOes, "N_JUNOes", NBinsJUNO);
//    myprt(lambJUNOibd, "lambJUNOibd", NBinsJUNO);
//    myprt(N_JUNOibd, "N_JUNOibd", NBinsJUNO);
    myprt(lambArgo, "lambArgo", NBinsArgo);
    myprt(N_ARGO, "N_ARGO", NBinsArgo);
#endif

#ifdef DEBUG
    printf("res=%lf\n", res);
#endif
    return res;
}

