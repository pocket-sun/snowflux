#include <random>
#include <cstdio>

using namespace std;

void possion_dis(double x[], size_t x_sz) {

    static default_random_engine e(42);
    for(size_t k = 0; k != x_sz; ++k) {
        if(k != x_sz-1) {
            printf("%.2lf,", x[k]);
        } else {
            printf("%.2lf\n", x[k]);
        }
    }
    for(size_t k = 0; k != x_sz; ++k) {
        poisson_distribution<unsigned long> pos(x[k]);
        if(k != x_sz-1) {
            printf("%lu,", pos(e));
        } else {
            printf("%lu\n", pos(e));
        }
    }

}
    
int main() {

// flux parameters:
// 2.5  2.5  2.5  9.5  12  15.6  5  5  5  

    double hkibd[10] = {
1166.2,
3731.8,
4322.3,
3470.1,
2275.8,
1337.1,
715.8 ,
351.0 ,
166.5 ,
81.8};

    double hkes[10] = {
153.2,
211.0,
135.2,
78.6 ,
43.4 ,
23.6 ,
13.0 ,
7.3  ,
4.1  ,
2.5};

    double dune[10] = {
310.3,
491.9,
383.3,
194.9,
102.3,
43.0 ,
17.1 ,
6.6  ,
2.5  ,
1.0};

    double resnov[10] = {
6248.01,
2833.22,
1470.94,
821.207,
480.601,
290.871,
180.55,
114.308,
73.5238,
47.9048};
    
   possion_dis(hkibd, 10); 
   possion_dis(hkes, 10); 
   possion_dis(dune, 10); 
   possion_dis(resnov, 10);

   return 0;
}
