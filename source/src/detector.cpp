#include "detector.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>
#include <globes/globes.h>   /* GLoBES library */

using namespace std;

constexpr double Detector::default_energies[200];

Detector::Detector(ExpName DetectorName) {

    switch(DetectorName) {
        case HyperK:
            strcpy(this->detconfigname, "wc100kt30prct");
            strcpy(this->channelname, "water");
            for(size_t k = 0; k != 7; ++k) {
                this->SelectedChannels[k] = k;
            }
            this->NumberChannels = 7;
            break;
        case DUNE:
            strcpy(this->detconfigname, "ar40kt");
            strcpy(this->channelname, "argon");
            for(size_t k = 6; k != 8; ++k) {
                this->SelectedChannels[k-6] = k;
            }
            this->NumberChannels = 2;
            break;
        default:
            fprintf(stderr, "wrong detector name");
            exit(-1);
            break;
    }

    char systemcall[256], glbfilename[64];
    // glb file name
    strcpy(glbfilename, this->detconfigname);
    strcat(glbfilename, this->channelname);
    strcat(glbfilename, ".glb");
    // system call
    strcpy(systemcall, "cd snowglobe && ./supernova.pl ");
    strcat(systemcall, this->channelname);
    strcat(systemcall, " ");
    strcat(systemcall, this->detconfigname);
    strcat(systemcall, " ");
    strcat(systemcall, glbfilename);
    strcat(systemcall, " && cd ..");
    if(system(systemcall) != 0) {
        fprintf(stderr, "fail to invoke supernova.pl");
        exit(-5);
    } else {
        strcpy(this->glb_file_name, "./snowglobe/");
        strcat(this->glb_file_name, glbfilename);
    }

    getChainFile();

}

void Detector::getChainFile() {

    char channel_file_name[128];
    /* Read the channels to process from the file */
    strcpy(channel_file_name, "./snowglobe/channels/channels_");
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

void Detector::generateRates() { 

    check();
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
    check();
    printf("\nenergy unit: MeV\n");
    for(size_t k = 0; k != NumberEnergies; ++k) {
        printf("%.3lf ", SelectedEnergies[k]*1e3);
    }
    printf("\n");
}

