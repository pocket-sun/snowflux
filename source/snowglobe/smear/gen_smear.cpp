#include <stdio.h>
#include <string.h>
#include <unistd.h> // Linux

int main(int argc, char *argv[]) {

    char *channelname = argv[1];
    char *detname = argv[2];
    
    char filename[64];
    strcpy(filename, "smear_");
    strcat(filename, channelname);
    strcat(filename, "_");
    strcat(filename, detname);
    strcat(filename, ".dat");
    if(access(filename, F_OK) == 0) {
        fprintf(stderr, "file: %s, already exsit!", filename);
        return -1;
    }
    FILE *pf = fopen(filename, "w");
    if(pf == NULL) {
        fprintf(stderr, "fail to create file: %s\n", filename); 
        return -2;
    }

    fprintf(pf, "energy(#%s_smear)<\n", channelname);
    fprintf(pf, "@energy = "); 
    for(size_t row = 0; row != 200; ++row) {
        fprintf(pf, "{0,199,"); // always assume 200 energies
        for(size_t col = 0; col != 200; ++col) {
            // operation region ----------------+
            if(col == row) {                 // |
                fprintf(pf, "1.");           // |
            } else {                         // |
                fprintf(pf, "0");            // |
            }                                // |
            // operation region ----------------+

            // end call of each column ---------+
            if(col != 199) {                 // |
                fprintf(pf, ",");            // |
            } else if(row != 199) {          // |
                fprintf(pf, "}:\n");         // |
            } else {                         // |
                fprintf(pf, "};\n>");        // |
            }                                // |
            // end call, empty below -----------+
        }
    }

    fclose(pf);

    return 0;
}

