#include <mpi.h>
#include <iostream>
#include <ctime>
#include "armadillo"
#include "operators.hpp"
#include "propagators.hpp"
#include "io_routines.hpp"
#include "rk_routines.hpp"

using namespace arma;
using namespace std;


int main(int argc, char** argv){
    int world_size;
    int world_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    cout << "Armadillo version: " << arma_version::as_string() << endl;
    cout << "MPI rank: " << world_rank << " / " << world_size << endl;

    clock_t begin,end;
    double norm2;
    int nt = 1025;
    sp_mat *E_F = new sp_mat();
    sp_mat *E_FCF = new sp_mat();
    mat pinvE_F;
    mat pinvE_FCF;

    Col<cx_double>  *din = new Col<cx_double>(3);
    Col<cx_double> *dout = new Col<cx_double>(3);
    (*din).fill({0.6,0.3});
    cout << rkconst::L_stable_SDIRK4 << endl;
    stability_function(rkconst::L_stable_SDIRK1, din, dout);
    stability_function(rkconst::L_stable_SDIRK2, din, dout);
    stability_function(rkconst::L_stable_SDIRK3, din, dout);
    stability_function(rkconst::L_stable_SDIRK4, din, dout);

    Col<int> ml(1);
    ml(0)     = 4;
    Col<double> lambda(2);
    Col<int> Nl(2);
    Nl(0)     = nt;
    lambda(0) = 0.9;
    for(int i = 1; i < ml.n_elem+1; i++){
        lambda(i) = std::pow(0.9, i+1);
        Nl(i)     = (Nl(i-1) - 1) / ml(i-1) + 1;
    }

    /* F-relaxation */
    begin = clock();
    get_E_F(E_F, lambda, Nl, ml);
    end = clock();
    cout << double(end-begin)/CLOCKS_PER_SEC << endl;
    // compute largest singular value
    begin = clock();
    // note: Armadillo 6.500.5 throws an error for matrices of size 10241^2 if norm() is used for sparse matrices
    // todo: check with later versions, if this is still an issue
    norm2 = norm(mat(*E_F), 2);
//    norm2 = norm(cx_mat(*E_F), 2);
    end = clock();
    cout << double(end-begin)/CLOCKS_PER_SEC << " " << norm2 << endl;
    // do the same for the pseudo-inverse
    begin = clock();
    bool pinvSuccess = true;
    pinvSuccess = pinv(pinvE_F, mat(*E_F));
    if(!pinvSuccess){
        cout << ">>>ERROR: Computing pseudo-inverse of error propagator failed.";
        throw;
    }
    norm2 = norm(pinvE_F, 2);
    end = clock();
    cout << double(end-begin)/CLOCKS_PER_SEC << " " << 1.0/norm2 << endl;

    /* FCF-relaxation */
    begin = clock();
    get_E_FCF(E_FCF, lambda, Nl, ml);
    end = clock();
    cout << double(end-begin)/CLOCKS_PER_SEC << endl;
    // compute largest singular value
    begin = clock();
    // note: Armadillo 6.500.5 throws an error for matrices of size 10241^2 if norm() is used for sparse matrices
    // todo: check with later versions, if this is still an issue
    norm2 = norm(mat(*E_FCF), 2);
//    norm2 = norm(cx_mat(*E_FCF), 2);
    end = clock();
    cout << double(end-begin)/CLOCKS_PER_SEC << " " << norm2 << endl;
    // do the same for the pseudo-inverse
    begin = clock();
    bool pinvSuccess2 = true;
    pinvSuccess2 = pinv(pinvE_FCF, mat(*E_FCF));
    if(!pinvSuccess2){
        cout << ">>>ERROR: Computing pseudo-inverse of error propagator failed.";
        throw;
    }
    norm2 = norm(pinvE_FCF, 2);
    end = clock();
    cout << double(end-begin)/CLOCKS_PER_SEC << " " << 1.0/norm2 << endl;

    // save to disk
    export_matrix(E_F, "E_F", raw_ascii);
    export_matrix(E_FCF, "E_FCF", raw_ascii);

    // load from disk
//    mat C;
//    C.load("E_F_real.txt");

    MPI_Finalize();
}
