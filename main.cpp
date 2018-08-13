#include <iostream>
#include <ctime>
#include "armadillo"
#include "propagators.h"

using namespace arma;
using namespace std;


int main(int argc, char** argv){
    cout << "Armadillo version: " << arma_version::as_string() << endl;

    clock_t begin,end;
    double norm2;
    int nt = 1025;

    Col<int> ml(1);
    ml(0)     = 4;
    Col<double> lambda(2);
    Col<int> Nl(2);
    Nl(0)     = nt;
    lambda(0) = 0.9;
    std::cout << "begin" << std::endl;
    for(int i = 1; i < ml.n_elem+1; i++){
        lambda(i) = std::pow(0.9, i+1);
        Nl(i)     = (Nl(i-1) - 1) / ml(i-1) + 1;
    }
    std::cout << "done" << std::endl;
    begin = clock();
    sp_mat E_F = get_E_F(lambda, Nl, ml);
    end = clock();
    cout << double(end-begin)/CLOCKS_PER_SEC << endl;
    begin = clock();
    norm2 = norm(E_F, 2);
    end = clock();
    cout << double(end-begin)/CLOCKS_PER_SEC << endl;
    mat fullE_F = mat(E_F);
    begin = clock();
    norm2 = norm(fullE_F, 2);
    end = clock();
    cout << double(end-begin)/CLOCKS_PER_SEC << endl;
    vec singVals(1);
    begin = clock();
    svds(singVals, E_F, 1);
    end = clock();
    cout << double(end-begin)/CLOCKS_PER_SEC << endl;

    // save to disk
    mat(E_F).save("E_F.txt", raw_ascii);

    // load from disk
    mat C;
    C.load("E_F.txt");

    return 0;
}
