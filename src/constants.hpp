#ifndef CONSTANTS_H
#define CONSTANTS_H

// constants referring to a specific Runge-Kutta method
namespace rkconst{
    // A-stable SDIRK methods
    const int A_stable_SDIRK1 = 101;
    const int A_stable_SDIRK2 = 102;
    const int A_stable_SDIRK3 = 103;
    const int A_stable_SDIRK4 = 104;
    // L-stable SDIRK methods
    const int L_stable_SDIRK1 = 201;
    const int L_stable_SDIRK2 = 202;
    const int L_stable_SDIRK3 = 203;
    const int L_stable_SDIRK4 = 204;
}

#endif
