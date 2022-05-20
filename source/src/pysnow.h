#ifndef PYSNOW
#define PYSNOW

#include "detector.h"
#include "pysnow.h"

extern "C" {
    extern Detector expr_hk; // declaration
    extern Detector expr_dune; // declaration
}

extern "C" {
    void rateGen(Detector &expr, double paras[], double dist, double[]);
    void init(Detector&);
    int pytest(Detector&);
}


#endif
