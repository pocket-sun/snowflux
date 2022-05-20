#include <cstdio>
#include <cstring>
#include <iostream>
#include "detector.h"
#include "pinched.h"
#include "pysnow.h"
#include <globes/globes.h>   /* GLoBES library */

using namespace std;

extern "C" {
    Detector expr_hk(HyperK); // definition
    Detector expr_dune(DUNE); // definition
}

// energy bin array should be entailed program lifetime
double hk_energies[NEHK], dune_energies[NEDUNE]; // GeV
void init_two_exps() { 

#ifdef DEBUG
    cout << "initialization finished" << endl;
#endif
    // energy settings
    for(size_t i = 0; i != NEHK; ++i) {
        // 5MeV~50MeV, 10bins
        hk_energies[i] = 5e-3 + (45e-3*i)/(NEHK-1);
    }
    for(size_t i = 0; i != NEDUNE; ++i) {
        // 6MeV~60MeV, 10bins
        dune_energies[i] = 6e-3 + (54e-3*i)/(NEDUNE-1);
    }
    expr_hk.setEnergyBins(hk_energies, NEHK);
    expr_dune.setEnergyBins(dune_energies, NEDUNE);

    // globe intialization 
    expr_hk.glbinit(); 
    expr_dune.glbinit(); 

}

void rateGen(Detector &expr, double paras[], double dist, double res[], size_t res_sz) {

//    double alpha[3] = {2.5, 2.5, 2.5};
//    double E0[3] = {9.5, 12, 15.6}; // MeV
//    double L[3] = {5e52, 5e52, 5e52}; // erg
//    double dist = 10.; // kpc
//    double selected_energies[15];  // GeV

    static bool not_init = 1; // definition excuted only once
    if(not_init) {init_two_exps(); not_init=0;}
#ifdef ERROR
    if(res_sz != expr.getNumberEnergies()-1) {
        fprintf(stderr, "%s: res_sz mismatch, res_sz=%ld,NumberEnergies-1=%ld",
            expr.getGlbFileName(), res_sz, expr.getNumberEnergies()-1);
        exit(-3);
    }
#endif

    double alpha[3] = {paras[0], paras[1], paras[2]};
    double E0[3] = {paras[3], paras[4], paras[5]};
    double L[3] = {paras[6]*1e52, paras[7]*1e52, paras[8]*1e52};
    if(GarchingFluence(alpha, E0, L, dist) == 0) {
        expr.glbreload();
        expr.generateRates(res, res_sz); // NumberEnergies-1
    }
}


// for test only
int pytest(Detector &expr) {
                  
    // user defined parameters 
    double alpha[3] = {2.5, 2.5, 2.5};
    double E0[3] = {9.5, 12, 15.6}; // MeV
    double L[3] = {5e52, 5e52, 5e52}; // erg
    double dist = 10.; // kpc

    double selected_energies[15];  // GeV
    for(size_t k = 0; k != 15; ++k) {
        selected_energies[k] = (k+1)*50./16*1e-3;
    }

    expr.setEnergyBins(selected_energies, 15);
    expr.glbinit();

    double tmp[14];
    for(size_t k = 0; k != 2; ++k) {
    E0[0] = E0[0] * (1.+k/4.);
    if(GarchingFluence(alpha, E0, L, dist) == 0) {
        expr.glbreload();
        expr.generateRates(tmp, 14);
        expr.printEnergyBins();
        expr.printEnergyIntervals();
    }
    }

    return 0;

}
