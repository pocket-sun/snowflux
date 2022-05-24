#ifndef PYSNOW
#define PYSNOW

#include "detector.h"
#include "pysnow.h" 

enum { 
    NEHK = 11, // NumberEnergies
    NBinsHK = 10, // NumberEnergies-1
    NEHKibd = 11, 
    NBinsHKibd = 10,
    NEHKes = 11,
    NBinsHKes = 10,
    NEDUNE = 11,
    NBinsDUNE = 10,
    NEJUNO = 2,
    NBinsJUNO = 1
};

extern "C" {
    extern Detector expr_hk; // declaration
    extern Detector expr_hkibd; // declaration
    extern Detector expr_hkes; // declaration
    extern Detector expr_dune; // declaration
    extern Detector expr_junopes; // definition
    extern Detector expr_junoes; // definition
    extern Detector expr_junoibd; // definition
}

extern "C" {
    void rateGen(Detector &expr, double paras[], double dist, double res[], size_t res_sz);
    void init_two_exps();
    int pytest(Detector&);
}


#endif
