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

void init(Detector &expr) { expr.glbinit(); }

void rateGen(Detector &expr, double paras[], double dist, double res[]) {

//    double alpha[3] = {2.5, 2.5, 2.5};
//    double E0[3] = {9.5, 12, 15.6}; // MeV
//    double L[3] = {5e52, 5e52, 5e52}; // erg
//    double dist = 10.; // kpc

    double selected_energies[15];  // GeV
    for(size_t k = 0; k != 15; ++k) {
        selected_energies[k] = (k+1)*50./16*1e-3;
    }
    expr.setEnergyBins(selected_energies, 15);
    
    double alpha[3] = {paras[0], paras[1], paras[2]};
    double E0[3] = {paras[3], paras[4], paras[5]};
    double L[3] = {paras[6]*1e52, paras[7]*1e52, paras[8]*1e52};
    if(GarchingFluence(alpha, E0, L, dist) == 0) {
        expr.glbreload();
        expr.generateRates(res, 42); // 14*3
    }
}


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

    double tmp[45];
    for(size_t k = 0; k != 4; ++k) {
    E0[0] = E0[0] * (1.+k/4.);
    if(GarchingFluence(alpha, E0, L, dist) == 0) {
        expr.glbreload();
        expr.generateRates(tmp, 42);
        expr.printEnergyBins();
    }
    }

    return 0;

}
