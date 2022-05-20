#include "pysnow.h"
#include <iostream>

using namespace std;

int main() {

    cout << "==============HyperK==============" << endl;
    pytest(expr_hk);
    cout << "===============DUNE===============" << endl;
    pytest(expr_dune);
    return 0;

}
