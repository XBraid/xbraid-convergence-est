#include "rk_routines.hpp"

// evaluate stability function
void stability_function(const int method, arma::Col<double> *z, arma::Col<double> *rz){
    // get Butcher tableau for requested method
    arma::mat           A;
    arma::Col<double>   b;
    arma::Col<double>   c;
    get_butcher_tableau(method, A, b, c);
    for(int evalIdx = 0; evalIdx < (*z).n_elem; evalIdx++){
        (*rz)(evalIdx) = arma::det(arma::eye(b.n_elem,b.n_elem) - (*z)(evalIdx) * A + (*z)(evalIdx) * arma::ones(b.n_elem, 1) * b.t())
                            / arma::det(arma::eye(b.n_elem,b.n_elem) - (*z)(evalIdx) * A);
    }
}

// evaluate stability function (complex version)
void stability_function(const int method, arma::Col<arma::cx_double> *z, arma::Col<arma::cx_double> *rz){
    // get Butcher tableau for requested method
    arma::mat           A;
    arma::Col<double>   b;
    arma::Col<double>   c;
    get_butcher_tableau(method, A, b, c);
    for(int evalIdx = 0; evalIdx < (*z).n_elem; evalIdx++){
        (*rz)(evalIdx) = arma::det(arma::eye(b.n_elem,b.n_elem) - (*z)(evalIdx) * A + (*z)(evalIdx) * arma::ones(b.n_elem, 1) * b.t())
                            / arma::det(arma::eye(b.n_elem,b.n_elem) - (*z)(evalIdx) * A);
    }
}

void get_butcher_tableau(const int method, arma::mat &A, arma::Col<double> &b, arma::Col<double> &c){
    switch(method){
        // implicit 1st-order method, SDIRK1
        // see [Dobrev et al. (2017)]
        case rkconst::L_stable_SDIRK1:
        {
            A.set_size(1,1);
            b.set_size(1);
            c.set_size(1);
            A(0,0)  = 1.0;
            b(0)    = 1.0;
            c(0)    = 1.0;
            break;
        }
        // implicit 2nd-order method, SDIRK2
        // see [Dobrev et al. (2017)]
        case rkconst::L_stable_SDIRK2:
        {
            A.set_size(2,2);
            b.set_size(2);
            c.set_size(2);
            double  gamma = 1.0 - 1.0 / sqrt(2.0);
            A       = {{gamma,      0.0},
                       {1.0-gamma,  gamma}};
            b       = {1.0-gamma,   gamma};
            c       = {gamma,       1.0};
            break;
        }
        // implicit 3rd-order method, SDIRK3
        // see [Dobrev et al. (2017)]
        // see talk by Butcher: https://www.math.auckland.ac.nz/~butcher/CONFERENCES/TRONDHEIM/trondheim.pdf
        // see MFEM, http://mfem.github.io/doxygen/html/ode_8cpp_source.html
        case rkconst::L_stable_SDIRK3:
        {
            A.set_size(3,3);
            b.set_size(3);
            c.set_size(3);
            double  q   = 0.435866521508458999416019;
            double  r   = 1.20849664917601007033648;
            double  s   = 0.717933260754229499708010;
            A       = {{q,      0.0,        0.0},
                       {s-q,    q,          0.0},
                       {r,      1.0-q-r,    q}};
            b       = {r,       1.0-q-r,    q};
            c       = {q,       s,          r};
            break;
        }
        // implicit 4th-order method, SDIRK4
        // see HairerWanner1996, Table IV.6.5
        // see DuarteDobbinsSmooke2016, Appendix C
        case rkconst::L_stable_SDIRK4:
        {
            A.set_size(5,5);
            b.set_size(5);
            c.set_size(5);
            A       = {{0.25,           0.0,            0.0,        0.0,        0.0},
                       {0.5,            0.25,           0.0,        0.0,        0.0},
                       {17.0/50.0,      -1.0/25.0,      0.25,       0.0,        0.0},
                       {371.0/1360.0,   -137.0/2720.0,  15.0/544.0, 0.25,       0.0},
                       {25.0/24.0,      -49.0/48.0,     125.0/16.0, -85.0/12.0, 0.25}};
            b       = {25.0/24.0,       -49.0/48.0,     125.0/16.0, -85.0/12.0, 0.25};
            c       = {0.25,            0.75,           11.0/20.0,  0.5,        1.0};
            break;
        }
    }
}