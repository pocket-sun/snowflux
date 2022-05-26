#ifndef PYSNVC
#define PYSNVC

#include "detectors.hpp"
#include "pySNvC.h"

enum { 
    NBinsRes = 10,
    NBinsArgo = 10
};

void getSpec(double (*fun)(float,float,fluxpara), float (*spec)[2], float Eth, 
    float binStep, int nbin, fluxpara flux);
void pygetSpecResNova(double inputs[], double outputs[]);
void pygetSpecArgo(double inputs[], double outputs[]);

extern "C"{
    int resnova_test();
}

#endif
