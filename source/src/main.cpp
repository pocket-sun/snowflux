#include <cstdio>
#include <cstring>
#include "detector.h"
#include "pinched.h"
#include <globes/globes.h>   /* GLoBES library */

using namespace std;

int main() {
                  
    Detector expr(DUNE);

    /* Initialize libglobes */
    char name[10]; strcpy(name, "snowglb");
    glbInit(name); // this function will invoke glb_init(char*) and set glbSetChannelPrintFunction
                   // passed by glb_builtin_channel_printf
                  
    /* Initialize experiment NFstandard.glb */
    // must be called only once!!!! too much time spent in this call
    char glbfilename_tmp[64]; strcpy(glbfilename_tmp, expr.getGlbFileName());
    // this return -1 not 0
    glbInitExperiment(glbfilename_tmp,&glb_experiment_list[0],&glb_num_of_exps);

    /* Zero oscillation parameters */
    double theta12 = 0.;
    double theta13 = 0.;
    double theta23 = 0.;
    double deltacp = 0.;
    double sdm = 0.;
    double ldm = 0;

    /* Initialize parameter vector(s) */
    glb_params true_values = glbAllocParams();
    glb_params test_values = glbAllocParams();

    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbSetDensityParams(true_values,1.0,GLB_ALL);
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);  
    glbSetDensityParams(test_values,1.0,GLB_ALL);

    /* The simulated data are computed */
    glbSetOscillationParameters(true_values);
    glbSetRates();

    // user defined parameters 
    
    double alpha[3] = {2.5, 2.5, 2.5};
    double E0[3] = {9.5, 12, 15.6}; // MeV
    double L[3] = {5e52, 5e52, 5e52}; // erg
    double dist = 10.; // kpc

    double selected_energies[15];  // GeV
    for(size_t k = 0; k != 15; ++k) {
        selected_energies[k] = (k+1)*50./16*1e-3;
    }

    for(size_t k = 0; k != 4; ++k) {
    E0[0] = E0[0] * (1.+k/4.);
    if(GarchingFluence(alpha, E0, L, dist) == 0) {
        // clear Globe settings
        glbFreeParams(true_values);
        glbFreeParams(test_values);
        glbClearExperimentList();

        // reload new experiment
        glbInitExperiment(glbfilename_tmp,&glb_experiment_list[0],&glb_num_of_exps);
        true_values = glbAllocParams();
        test_values = glbAllocParams();
        glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
        glbSetDensityParams(true_values,1.0,GLB_ALL);
        glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);  
        glbSetDensityParams(test_values,1.0,GLB_ALL);
        glbSetOscillationParameters(true_values);
        glbSetRates();

        expr.setEnergyBins(selected_energies, 15);
        expr.generateRates();
        expr.printEnergyBins();
    }
    }

    return 0;

}
