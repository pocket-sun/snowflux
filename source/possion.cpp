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
3693.83,8237.44,9509.19,7634.23,5006.77,2941.70,1574.82,772.17,366.34,179.95};

    double hkes[10] = {
942.02,472.84,297.35,172.97,95.50,51.84,28.55,16.12,9.06,5.52};

    double dune[10] = {
456.43,550.78,406.51,203.65,106.49,44.73,17.85,6.93,2.66,1.07};

    double resnov[10] = {
6248.01,2833.22,1470.94,821.21,480.60,290.87,180.55,114.31,73.52,47.90};

    double junopes[1] = {95.04};
    double junoes[1] = {135.25};
    double junoibd[1] = {63.40};
    
   possion_dis(hkibd, 10); 
   possion_dis(hkes, 10); 
   possion_dis(dune, 10); 
   possion_dis(resnov, 10);
   possion_dis(junopes, 1);
   possion_dis(junoes, 1);
   possion_dis(junoibd, 1);

   return 0;
}
