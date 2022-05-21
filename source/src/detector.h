#ifndef DETECTOR
#define DETECTOR

#include <cstdio>
#include <cstdlib>
#include "detector.h"
#include <globes/globes.h>   /* GLoBES library */

#define MAXCHANS 32

// to add an experiment
// case in Detector(ExpName) in detector.h should be added as well
enum ExpName {
    HyperK = 1024,
    HyperKibd = 1025,
    HyperKES = 1026,
    DUNE = 2049
};
        
class Detector {
public:
    // copy control
    Detector()=delete;
    Detector(const Detector&)=delete;
    Detector &operator=(const Detector&)=delete;
    Detector(ExpName DetectorName);

    // method
    void setEnergyBins(double ebins[], size_t binsNumber); // GeV
    void generateRates(double res[], size_t res_size);
    void printEnergyBins() const;
    void printEnergyIntervals() const;
    char *getGlbFileName() { return glb_file_name; }
    size_t getNumberEnergies() const { return NumberEnergies; }
    unsigned getFnum() const;
    void glbreload();
    void glbinit();

    // destructor
    ~Detector() {
        if(init_toggle == 1) {
            // clear Globe settings
            glbFreeParams(true_values);
            glbFreeParams(test_values);
        }
    }

private:
    void checkinit() const;
    void checkene() const;
    void getChainFile();
    void createGLBFile(unsigned);

    double *SelectedEnergies = NULL;// GeV
    size_t NumberEnergies = 0;
    size_t SelectedChannels[MAXCHANS];
    size_t NumberChannels;

    double theta12 = 0., theta13 = 0., theta23 = 0.;
    double deltacp = 0., sdm = 0., ldm = 0;
    /* Initialize parameter vector(s) */
    glb_params true_values;
    glb_params test_values;
    bool init_toggle = 0; // declaration

    char glb_file_name[128], detconfigname[64], channelname[64];
    int chan_num[MAXCHANS];
    int num_target_factors[MAXCHANS];
    char chan_name[MAXCHANS][64];
    static constexpr double default_energies[200] = { // decl && initialization
0.7488e-3,
1.2460e-3,
1.7440e-3,
2.2410e-3,
2.7390e-3,
3.2360e-3,
3.7340e-3,
4.2310e-3,
4.7290e-3,
5.2260e-3,
5.7240e-3,
6.2210e-3,
6.7190e-3,
7.2160e-3,
7.7140e-3,
8.2110e-3,
8.7090e-3,
9.2060e-3,
9.7040e-3,
10.2000e-3,
10.7000e-3,
11.2000e-3,
11.6900e-3,
12.1900e-3,
12.6900e-3,
13.1900e-3,
13.6800e-3,
14.1800e-3,
14.6800e-3,
15.1800e-3,
15.6700e-3,
16.1700e-3,
16.6700e-3,
17.1700e-3,
17.6600e-3,
18.1600e-3,
18.6600e-3,
19.1600e-3,
19.6500e-3,
20.1500e-3,
20.6500e-3,
21.1500e-3,
21.6400e-3,
22.1400e-3,
22.6400e-3,
23.1400e-3,
23.6300e-3,
24.1300e-3,
24.6300e-3,
25.1300e-3,
25.6200e-3,
26.1200e-3,
26.6200e-3,
27.1200e-3,
27.6100e-3,
28.1100e-3,
28.6100e-3,
29.1100e-3,
29.6000e-3,
30.1000e-3,
30.6000e-3,
31.1000e-3,
31.5900e-3,
32.0900e-3,
32.5900e-3,
33.0900e-3,
33.5800e-3,
34.0800e-3,
34.5800e-3,
35.0800e-3,
35.5700e-3,
36.0700e-3,
36.5700e-3,
37.0700e-3,
37.5600e-3,
38.0600e-3,
38.5600e-3,
39.0600e-3,
39.5500e-3,
40.0500e-3,
40.5500e-3,
41.0500e-3,
41.5400e-3,
42.0400e-3,
42.5400e-3,
43.0400e-3,
43.5300e-3,
44.0300e-3,
44.5300e-3,
45.0300e-3,
45.5200e-3,
46.0200e-3,
46.5200e-3,
47.0200e-3,
47.5100e-3,
48.0100e-3,
48.5100e-3,
49.0100e-3,
49.5000e-3,
50.0000e-3,
50.5000e-3,
51.0000e-3,
51.4900e-3,
51.9900e-3,
52.4900e-3,
52.9900e-3,
53.4800e-3,
53.9800e-3,
54.4800e-3,
54.9800e-3,
55.4700e-3,
55.9700e-3,
56.4700e-3,
56.9700e-3,
57.4600e-3,
57.9600e-3,
58.4600e-3,
58.9600e-3,
59.4500e-3,
59.9500e-3,
60.4500e-3,
60.9500e-3,
61.4400e-3,
61.9400e-3,
62.4400e-3,
62.9400e-3,
63.4300e-3,
63.9300e-3,
64.4300e-3,
64.9300e-3,
65.4200e-3,
65.9200e-3,
66.4200e-3,
66.9200e-3,
67.4100e-3,
67.9100e-3,
68.4100e-3,
68.9100e-3,
69.4000e-3,
69.9000e-3,
70.4000e-3,
70.9000e-3,
71.3900e-3,
71.8900e-3,
72.3900e-3,
72.8900e-3,
73.3800e-3,
73.8800e-3,
74.3800e-3,
74.8800e-3,
75.3700e-3,
75.8700e-3,
76.3700e-3,
76.8700e-3,
77.3600e-3,
77.8600e-3,
78.3600e-3,
78.8600e-3,
79.3500e-3,
79.8500e-3,
80.3500e-3,
80.8500e-3,
81.3400e-3,
81.8400e-3,
82.3400e-3,
82.8400e-3,
83.3300e-3,
83.8300e-3,
84.3300e-3,
84.8300e-3,
85.3200e-3,
85.8200e-3,
86.3200e-3,
86.8200e-3,
87.3100e-3,
87.8100e-3,
88.3100e-3,
88.8100e-3,
89.3000e-3,
89.8000e-3,
90.3000e-3,
90.8000e-3,
91.2900e-3,
91.7900e-3,
92.2900e-3,
92.7900e-3,
93.2800e-3,
93.7800e-3,
94.2800e-3,
94.7800e-3,
95.2700e-3,
95.7700e-3,
96.2700e-3,
96.7700e-3,
97.2600e-3,
97.7600e-3,
98.2600e-3,
98.7600e-3,
99.2500e-3,
99.7500e-3}; // GeV

};

// check glbinit() and setEnergyBins
inline
void Detector::checkinit() const {
    if(init_toggle == 0) {
        std::fprintf(stderr, "you must initialize the glb by Detector::glbinit()");
        exit(-23);
    }
}
inline
void Detector::checkene() const {
    if(SelectedEnergies == NULL) {
        std::fprintf(stderr, "energies muse be set by Detector::setEnergyBins()");
        exit(-4);
    }
}

#endif
