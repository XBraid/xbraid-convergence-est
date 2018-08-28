#include "rk_routines.hpp"

/**
    evaluate stability function for given Runge-Kutta method
*/
void stability_function(const int method,       ///< Runge-Kutta method, see constants.hpp
                        arma::Col<double> *z,   ///< z = dt*eta, with time step size dt and spatial eigenvalue eta
                        arma::Col<double> *rz   ///< on return, the value of the stability function
                        ){
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

/**
    evaluate stability function for given Runge-Kutta method (complex version)
*/
void stability_function(const int method,               ///< Runge-Kutta method, see constants.hpp
                        arma::Col<arma::cx_double> *z,  ///< z = dt*eta, with time step size dt and spatial eigenvalue eta
                        arma::Col<arma::cx_double> *rz  ///< on return, the value of the stability function
                        ){
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

/**
    evaluate stability function for Runge-Kutta method with provided Butcher tableau
*/
void stability_function(arma::mat A,            ///< Butcher coefficients, matrix A
                        arma::Col<double> b,    ///< Butcher coefficients, vector b
                        arma::Col<double> c,    ///< Butcher coefficients, vector c
                        arma::Col<double> *z,   ///< z = dt*eta, with time step size dt and spatial eigenvalue eta
                        arma::Col<double> *rz   ///< on return, the value of the stability function
                        ){
    for(int evalIdx = 0; evalIdx < (*z).n_elem; evalIdx++){
        (*rz)(evalIdx) = arma::det(arma::eye(b.n_elem,b.n_elem) - (*z)(evalIdx) * A + (*z)(evalIdx) * arma::ones(b.n_elem, 1) * b.t())
                            / arma::det(arma::eye(b.n_elem,b.n_elem) - (*z)(evalIdx) * A);
    }
}

/**
    evaluate stability function for Runge-Kutta method with provided Butcher tableau (complex version)
*/
void stability_function(arma::mat A,                    ///< Butcher coefficients, matrix A
                        arma::Col<double> b,            ///< Butcher coefficients, vector b
                        arma::Col<double> c,            ///< Butcher coefficients, vector c
                        arma::Col<arma::cx_double> *z,  ///< z = dt*eta, with time step size dt and complex spatial eigenvalue eta
                        arma::Col<arma::cx_double> *rz  ///< on return, the value of the stability function
                        ){
    for(int evalIdx = 0; evalIdx < (*z).n_elem; evalIdx++){
        (*rz)(evalIdx) = arma::det(arma::eye(b.n_elem,b.n_elem) - (*z)(evalIdx) * A + (*z)(evalIdx) * arma::ones(b.n_elem, 1) * b.t())
                            / arma::det(arma::eye(b.n_elem,b.n_elem) - (*z)(evalIdx) * A);
    }
}

/** 
    sets the coefficients of a given Runge-Kutta method corresponding to Butcher tableau of form
        <table>
            <tr><td>c<td>A
            <tr><td> <td>b^T
        </table>
*/
void get_butcher_tableau(const int method,      ///< Runge-Kutta method, see constants.hpp
                         arma::mat &A,          ///< matrix A of Butcher tableau
                         arma::Col<double> &b,  ///< vector b of Butcher tableau
                         arma::Col<double> &c   ///< vector c of Butcher tableau
                         ){
    switch(method){
        case rkconst::A_stable_SDIRK2:
        {
            A.set_size(2,2);
            b.set_size(2);
            c.set_size(2);
            A       = {{0.25,   0.0},
                       {0.5,    0.25}};
            b       = {0.5,     0.5};
            c       = {0.25,    0.75};
        }
        case rkconst::A_stable_SDIRK3:
        {
            A.set_size(2,2);
            b.set_size(2);
            c.set_size(2);
            double gamma = (3.0 + sqrt(3.0)) / 6.0;
            A       = {{gamma,          0.0},
                       {1.0-2.0*gamma,  gamma}};
            b       = {0.5,             0.5};
            c       = {gamma,           1.0-gamma};
        }
        case rkconst::A_stable_SDIRK4:
        {
            A.set_size(3,3);
            b.set_size(3);
            c.set_size(3);
            double  q = 1.0 / sqrt(3.0) * cos(rkconst::pi / 18.0) + 0.5;
            double  r = 1.0 / (6.0 *(2.0 * q - 1.0) * (2.0 * q - 1.0));
            A       = {{q,      0.0,        0.0},
                       {0.5-q,  q,          0.0},
                       {2.0*q,  1.0-4.0*q,  q}};
            b       = {r,       1.0-2.0*r,  r};
            c       = {q,       0.5,        1.0-q};
        }
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
        default:
        {
            std::cerr << "Runge-Kutta method " << method << " not implemented." << std::endl;
            throw;
        }
    }
}