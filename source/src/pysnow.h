#ifndef PYSNOW
#define PYSNOW

#include "detector.h"
#include "pysnow.h"

enum { 
    NEHK = 11, // NumberEnergies
    NBinsHK = 10, // NumberEnergies-1
    NEDUNE = 11,
    NBinsDUNE = 10
};

extern "C" {
    extern Detector expr_hk; // declaration
    extern Detector expr_dune; // declaration
}

extern "C" {
    void rateGen(Detector &expr, double paras[], double dist, double res[], size_t res_sz);
    void init_two_exps();
    int pytest(Detector&);
}


#endif
