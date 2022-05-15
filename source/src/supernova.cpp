#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/stat.h>
#include <map>

#include <globes/globes.h>   /* GLoBES library */
#include "pinched.h"

using namespace std;

int ret;
char channel_file_name[256];
static int maxchans = 32;

int main(int argc, char *argv[])
{ 

    // flux parameters 
    double alpha[3] = {2.5, 2.5, 2.5};
    double E0[3] = {9.5, 12, 15.6}; // MeV
    double L[3] = {5e52, 5e52, 5e52}; // erg
    double dist = 10.; // kpc
    if(GarchingFluence(alpha, E0, L, dist) != 0) {
        fprintf(stderr, "error in generating flux");
        exit(-2);
    }
    // channels and energy binning
    // energy difference must greater than 0.5 MeV
    int selected_channels[] = {0, 1}; // refer to snowglobe/channels/xx.dat second cloumn
    double selected_energies[15];
    for(size_t k = 0; k != 15; ++k) {
        selected_energies[k] = (k+1)*50./16*1e-3;
        printf("%.3lf ", selected_energies[k]);
    }
    printf("\n");
    // double selected_energies[] = {3.5e-3, 10.e-3, 50.e-3}; // GeV, (a0~a1, ... , an-1~an)


    /* Initialize libglobes */
    glbInit(argv[0]); // this function will invoke glb_init(char*) and set glbSetChannelPrintFunction
                    // passed by glb_builtin_channel_printf

    if (argc<2) {
    printf("Arguments required: channels filename\n");
    exit(-3);
    }

		
    /* Read the channels to process from the file */
    strcat(channel_file_name, "snowglobe/channels/channels_");
    strcat(channel_file_name, argv[1]);
    strcat(channel_file_name, ".dat");
    //strncpy(channel_file_name,argv[1],256);
    //printf("Channels from %s\n",channel_file_name);

    FILE* fp_chans = fopen(channel_file_name,"r");

    if (fp_chans == NULL) {
    printf ("Cannot open file \n");
    exit(-5);
    }

    int chan_num[maxchans];
    char chan_name[maxchans][256];
    char chbuf[1000];

    int num_chans=0;

    char cp;
    char flav;

    int num_target_factors[maxchans];

    int skipline = 0;
    char* eofcheck; /* need to check return of fgets for whether line was read successfully */
    while (!feof(fp_chans)) {
    skipline = 0;
    eofcheck = fgets(chbuf,1000,fp_chans);
    for (int ichar = 0; ichar != 5; ichar++){
      if (chbuf[ichar] == '%'){
        skipline = 1;
      }
    }
    if (skipline) {continue;}
    if (eofcheck != NULL) {
      ret = sscanf(chbuf,"%s %i %s %s %i",chan_name[num_chans],
        &chan_num[num_chans],&cp, &flav, &num_target_factors[num_chans]);
      num_chans++;
    }
    }

    fclose(fp_chans);


    //printf("Number of channels found: %i\n",num_chans);


    /* Initialize experiment NFstandard.glb */
    char glbfilename[64] = "./snowglobe/supernova.glb"; // CLS, avoid C++ warning 
    glbInitExperiment(glbfilename,&glb_experiment_list[0],&glb_num_of_exps); 
 

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


    size_t ifile;
  
    double energy, countings; // GeV, numbers/0.5MeV
    //
    size_t num_channels = sizeof(selected_channels) / sizeof(selected_channels[0]);
    size_t num_energies = sizeof(selected_energies) / sizeof(selected_energies[0]);
    double countings_per_interval[num_energies-1] = {0.}; // all initialized to 0.

    char line[256];
    map<double, double> mcountings_per_interval;
    for(ifile=0; ifile != num_channels; ++ifile) {
        // create a FILE in memory
        char *membuf; // dynamically allocated, need free at last
        size_t buflen;
        FILE *fmem = open_memstream(&membuf, &buflen);
        if(fmem == NULL) {
            fprintf(stderr, "fail to open_memstream");
            exit(-1);
        }
        // smeared, weighted 
        ret = glbShowChannelRates(fmem,0,chan_num[selected_channels[ifile]],GLB_POST,GLB_W_EFF,GLB_W_BG);
        fflush(fmem);
        for(size_t l, k = 0; k != buflen;) {
            l = 0;
            while((line[l++]=membuf[k++]) != '\n');
            line[l] = '\0';
            //printf("%s", line);
            int check = sscanf(line, "%lf %lf", &energy, &countings);
            if(check == 2)
                mcountings_per_interval[energy] += countings 
                    * num_target_factors[selected_channels[ifile]];
        }
        fclose(fmem);
        free(membuf);
    }
    const double de = 0.25e-3; // 0.5e-3/2
    size_t cnt = 0, index = 0;
    double elast; // energy last
    for(const auto &k: mcountings_per_interval) {
        cnt++;
        if(cnt>1 && fabs(k.first-elast-2*de) > 5e-4) {
            printf("%lf, %lf\n", k.first, elast);
            fprintf(stderr, "snowglobe interval error");
            exit(-10);
        }
        elast = k.first;
        if(selected_energies[index] >= elast-de 
            && selected_energies[index] <= elast+de)
            selected_energies[index++] = elast-de;
    }
    if(index != num_energies) {
        fprintf(stderr, "energy limit overflow");
        exit(-3);
    }
    index = 0;
    bool hold = 0;
    const double ethre = 1e-4; // 1e-5
    for(const auto &k: mcountings_per_interval) {
        if(hold || (index!=num_energies-1 
            && fabs(k.first-de-selected_energies[index])<ethre)) {
            hold = 1;
            countings_per_interval[index] += k.second;
            if(fabs(k.first+de-selected_energies[index+1])<ethre) {
                ++index;
                hold = 0;
            }
        }
    }
    if(index != num_energies-1) {
        fprintf(stderr, "sth went wrong");
        exit(-23);
    }
    
    printf("\n------ original spectral -------\n");
    printf("Ev (MeV)                countings\n");
    for(auto &k: mcountings_per_interval) {
        printf("%-10.4lf              %-10.4lf\n", k.first*1e3, k.second);
    }
    printf("\n-------- final spectral --------\n");
    printf("Ev (MeV)                countings\n");
    for(size_t k = 0; k != num_energies-1; ++k) {
        printf("[%-8.5lf, %-8.5lf]    %-11.1lf\n",
            selected_energies[k]*1e3, selected_energies[k+1]*1e3, countings_per_interval[k]);
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
  
    /* Destroy parameter vector(s) */
    glbFreeParams(true_values);
    glbFreeParams(test_values); 

  
    exit(0);
}
