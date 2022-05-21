#include "detector.h"
#include "pinched.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>
#include <globes/globes.h>   /* GLoBES library */
#include <sys/types.h>
#include <unistd.h>

using namespace std;

constexpr double Detector::default_energies[200]; // declaration


Detector::Detector(ExpName DetectorName) {

    size_t chcnt = 0;
    switch(DetectorName) {
        case HyperK:
            strcpy(this->detconfigname, "hyperk");
            strcpy(this->channelname, "water");
            for(size_t k = 0; k != 7; ++k) {
                ++chcnt;
                this->SelectedChannels[k] = k;
            }
            this->NumberChannels = chcnt;
            break;
        case HyperKibd:
            strcpy(this->detconfigname, "hyperk");
            strcpy(this->channelname, "water");
            for(size_t k = 0; k != 1; ++k) {
                ++chcnt;
                this->SelectedChannels[k] = k;
            }
            this->NumberChannels = chcnt;
            break;
        case HyperKES:
            strcpy(this->detconfigname, "hyperk");
            strcpy(this->channelname, "water");
            for(size_t k = 1; k != 7; ++k) {
                ++chcnt;
                this->SelectedChannels[k-1] = k;
            }
            this->NumberChannels = chcnt;
            break;
        case DUNE:
            strcpy(this->detconfigname, "ar40kt");
            strcpy(this->channelname, "argon");
            for(size_t k = 0; k != 8; ++k) {
                ++chcnt;
                this->SelectedChannels[k] = k;
            }
            this->NumberChannels = chcnt;
            break;
        default:
            fprintf(stderr, "wrong detector name");
            exit(-1);
            break;
    }

    getChainFile();

}

void Detector::createGLBFile(unsigned flux_num) {

#ifdef PARADEBUG
    cout << "createGLB:" << getFnum() << endl;
#endif
    char sflux_num[64]; sprintf(sflux_num, "%u", flux_num);
    char systemcall[256], glbfilename[64];
    // glb file name
    strcpy(glbfilename, sflux_num);
    strcat(glbfilename, this->detconfigname);
    strcat(glbfilename, this->channelname);
    char nchannel[4];
    sprintf(nchannel, "%lu", this->NumberChannels);
    strcat(glbfilename, nchannel);
    strcat(glbfilename, ".glb");
    // system call
    char *ppath;
    ppath = getenv("SNOWGLOBEPATH"); // "/" should be included in
    if(ppath == NULL) {fprintf(stderr, "SNOWGLOBEPATH should be set"); exit(-4);}
    strcpy(systemcall, "cd ");
    strcat(systemcall, ppath);
    strcat(systemcall, "snowglobe && ./supernova.pl ");
    strcat(systemcall, this->channelname);
    strcat(systemcall, " ");
    strcat(systemcall, this->detconfigname);
    strcat(systemcall, " ");
    strcat(systemcall, glbfilename);
    strcat(systemcall, " ");
    strcat(systemcall, sflux_num);
    strcat(systemcall, " && cd ..");
    if(system(systemcall) != 0) {
        fprintf(stderr, "fail to invoke supernova.pl");
        exit(-5);
    } else {
        strcpy(this->glb_file_name, ppath);
        strcat(this->glb_file_name, "snowglobe/");
        strcat(this->glb_file_name, glbfilename);
    }

}

void Detector::getChainFile() {

    char channel_file_name[128];
    /* Read the channels to process from the file */
    char *ppath;
    ppath = getenv("SNOWGLOBEPATH"); // no "/" should included in
    if(ppath == NULL) {fprintf(stderr, "SNOWGLOBEPATH should be set"); exit(-4);}
    strcpy(channel_file_name, ppath);
    strcat(channel_file_name, "snowglobe/channels/channels_");
    strcat(channel_file_name, channelname);
    strcat(channel_file_name, ".dat");

    FILE* fp_chans = fopen(channel_file_name,"r");
    if (fp_chans == NULL) {
        printf ("Cannot open file \n");
        exit(-5);
    }

    char chbuf[1000];
    size_t nchans = 0; // type int will lead arcane error!!!!!!! CLS
    char cp, flav;
    int skipline = 0;
    char* eofcheck; /* need to check return of fgets for whether line was read successfully */
    while (!feof(fp_chans)) {
        skipline = 0;
        eofcheck = fgets(chbuf,1000,fp_chans);
        for (int ichar = 0; ichar != 5; ichar++){
            if (chbuf[ichar] == '%') skipline = 1;
        }
        if (skipline) continue;
        if (eofcheck != NULL) {
            sscanf(chbuf,"%s %i %s %s %i",chan_name[nchans], &chan_num[nchans],
                &cp, &flav, &num_target_factors[nchans]);
            nchans++;
        }
    }
    fclose(fp_chans);
}

// res[] used to store rates between energy interval central points
// so res_size should fit with NumberEnergies-1
void Detector::generateRates(double res[], size_t res_size) { 

    checkinit();
    checkene();
#ifdef ERROR
    if(res_size != NumberEnergies-1) {
        fprintf(stderr, "res size error");
        exit(-9);
    }
#endif
    size_t index;
  
    double energy, countings; // GeV, numbers/0.5MeV
    double countings_per_interval[NumberEnergies-1] = {0.}; // all initialized to 0.

    char line[256];
    map<double, double> mcountings_per_interval;

    for(index = 0; index != NumberChannels; ++index) {
        // create a FILE in memory
        char *membuf; // dynamically allocated, need free at last
        size_t buflen;
        FILE *fmem = open_memstream(&membuf, &buflen);
        if(fmem == NULL) {
            fprintf(stderr, "fail to open_memstream");
            exit(-1);
        }
        // smeared, weighted 
        glbShowChannelRates(fmem,0,SelectedChannels[index],GLB_POST,GLB_W_EFF,GLB_W_BG);
        fflush(fmem);
        for(size_t l, k = 0; k != buflen;) {
            l = 0;
            while((line[l++]=membuf[k++]) != '\n');
            line[l] = '\0';
            //printf("%s", line);
            int check = sscanf(line, "%lf %lf", &energy, &countings);
            if(check == 2)
                mcountings_per_interval[energy] += countings 
                    * num_target_factors[SelectedChannels[index]];
        }
        fclose(fmem);
        free(membuf);
    }
#ifdef ERROR
    index = 0;
    for(const auto &k: mcountings_per_interval) {
        //printf("%lf, %lf\n", k.first, k.second);
        if(fabs(k.first-default_energies[index++])>1e-7) {
            printf("%lf, %lf\n", k.first, default_energies[index-1]);
            fprintf(stderr, "snowglobe interval error");
            exit(-10);
        }
    }
    if(index != 200) {
        fprintf(stderr, "default energies error");
        exit(-12);
    }
#endif
    index = 0;
    // include [a0, a1), [a1, a2), ..., [an-1, an]
    // ai itself is interval
    for(const auto &k: mcountings_per_interval) {
        if(k.first>=SelectedEnergies[index] && k.first<SelectedEnergies[index+1]) {
            countings_per_interval[index] += k.second;
        } else if(k.first >= SelectedEnergies[index+1]) {
            if(index != NumberEnergies-2) {
                countings_per_interval[++index] += k.second;
            } else {
                countings_per_interval[index] += k.second;
                break;
            }
        }
    }


#ifdef DEBUG    
    printf("%s:", glb_file_name);
    printf("\n------ original spectral -------\n");
    printf("Ev (MeV)                countings\n");
    for(auto &k: mcountings_per_interval) {
        printf("%-10.4lf              %-10.4lf\n", k.first*1e3, k.second);
    }
    printf("\n-------- final spectral --------\n");
    printf("Ev (MeV)                countings\n");
    for(index = 0; index != NumberEnergies-2; ++index) {
        printf("[%-8.5lf, %-8.5lf]    %-11.1lf\n",
            SelectedEnergies[index]*1e3-0.25, SelectedEnergies[index+1]*1e3-0.25,
                countings_per_interval[index]);
    }
    printf("[%-8.5lf, %-8.5lf]    %-11.1lf\n",
            SelectedEnergies[index]*1e3-0.25, SelectedEnergies[index+1]*1e3+0.25,
                countings_per_interval[index]);
#endif
    for(index = 0; index != NumberEnergies-1; ++index) {
        res[index] = countings_per_interval[index];
    }
            
    // for now, no background is assumed
    /*
    sprintf(bgfile,"backgrounds/bg_chan_%s.dat",expt_config_name);

    if (stat(bgfile,&buf) == 0) {
        printf("background iterface needed\n");
    } else {
        //printf("No background file\n");
    }
    */
}

// GeV
void Detector::setEnergyBins(double *ebins, size_t binsNumber) {
    SelectedEnergies = ebins;
    NumberEnergies = binsNumber;
    size_t index = 0, k = 0;
    while(index != binsNumber && k != 199) {
        if(ebins[index] >= default_energies[k] 
            && ebins[index] < default_energies[k+1]) {
            ebins[index++] = default_energies[k];
        }
        ++k;
    }
}

void Detector::printEnergyBins() const {
    checkene();
    printf("\nbin central energies, unit: MeV\n");
    for(size_t k = 0; k != NumberEnergies; ++k) {
        printf("%-8.3lf ", SelectedEnergies[k]*1e3);
        if((k+1) % 5 == 0) printf("\n");
    }
    printf("\n");
}

void Detector::printEnergyIntervals() const {
    checkene();
    size_t index;
    printf("\nenergy intervals, unit: MeV\n");
    for(index = 0; index != NumberEnergies-2; ++index) {
        printf("[%-8.5lf, %-8.5lf]\n", SelectedEnergies[index]*1e3-0.25
            , SelectedEnergies[index+1]*1e3-0.25);
    }
    printf("[%-8.5lf, %-8.5lf]\n", SelectedEnergies[index]*1e3-0.25
        , SelectedEnergies[index+1]*1e3+0.25);
}

void Detector::glbinit() {

    if(init_toggle == 0) {

    static double alpha[3] = {2.5, 2.5, 2.5};
    static double E0[3] = {9.5, 12, 15.6}; // MeV
    static double L[3] = {5e52, 5e52, 5e52}; // erg
    static double dist = 10.; // kpc
#ifdef PARADEBUG
    cout << "init:" << getFnum() << endl;
#endif
    GarchingFluence(alpha, E0, L, dist, getFnum());
    createGLBFile(getFnum());

    /* Initialize libglobes */
    char name[10]; strcpy(name, "snowglb");
    glbInit(name); // this function will invoke glb_init(char*) and set glbSetChannelPrintFunction
                   // passed by glb_builtin_channel_printf
                  
    /* Initialize experiment NFstandard.glb */
    // must be called only once!!!! too much time spent in this call
    char glbfilename_tmp[64]; 
    strcpy(glbfilename_tmp, this->glb_file_name);
    // this return -1 not 0
    glbInitExperiment(glbfilename_tmp,&glb_experiment_list[0],&glb_num_of_exps);
    true_values = glbAllocParams();
    test_values = glbAllocParams();

    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbSetDensityParams(true_values,1.0,GLB_ALL);
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);  
    glbSetDensityParams(test_values,1.0,GLB_ALL);

    /* The simulated data are computed */
    glbSetOscillationParameters(true_values);
    glbSetRates();

    init_toggle = 1;

    }

}

// every time load new flux, or modify experimental parameters
// ,or change glbfile, glbreload() should be called
// in principle, before the call to generate event rates, glbreload()
// should be invoked
void Detector::glbreload() {

    checkinit();
    // clear Globe settings
    glbFreeParams(true_values);
    glbFreeParams(test_values);
    glbClearExperimentList();

    // reload new experiment
    char glbfilename_tmp[64]; 
    strcpy(glbfilename_tmp, this->glb_file_name);
#ifdef PARADEBUG
    cout << "reload:" << glb_file_name << endl;
#endif
    glbInitExperiment(glbfilename_tmp,&glb_experiment_list[0],&glb_num_of_exps);
    true_values = glbAllocParams();
    test_values = glbAllocParams();
    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbSetDensityParams(true_values,1.0,GLB_ALL);
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);  
    glbSetDensityParams(test_values,1.0,GLB_ALL);
    glbSetOscillationParameters(true_values);
    glbSetRates();

}

unsigned Detector::getFnum() const {
    return getpid();
}
