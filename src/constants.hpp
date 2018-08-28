#ifndef CONSTANTS_H
#define CONSTANTS_H

// constants referring to a specific Runge-Kutta method
namespace rkconst{
    // A-stable SDIRK methods
    const int A_stable_SDIRK2 = 102;    ///< 2nd-order A-stable SDIRK method
    const int A_stable_SDIRK3 = 103;    ///< 3rd-order A-stable SDIRK method
    const int A_stable_SDIRK4 = 104;    ///< 4th-order A-stable SDIRK method
    // L-stable SDIRK methods
    /** 1st-order L-stable SDIRK method, Butcher tableau:
        <table>
            <tr>    <td>    \f$1\f$ <td>    \f$1\f$
            <tr>    <td>            <td>    \f$1\f$
        </table>
        See [Dobrev et al. (2017)].
    */
    const int L_stable_SDIRK1 = 201;
    /** 2nd-order L-stable SDIRK method, Butcher tableau:
        <table>
            <tr>    <td>    \f$\gamma\f$    <td>    \f$\gamma\f$   <td>    \f$0\f$
            <tr>    <td>    \f$1\f$         <td>    \f$1-\gamma\f$ <td>    \f$\gamma\f$
            <tr>    <td>                    <td>    \f$1-\gamma\f$ <td>    \f$\gamma\f$
        </table>
        with \f$\gamma = 1.0 - 1.0 / \sqrt{2.0}\f$.
        See [Dobrev et al. (2017)].
    */
    const int L_stable_SDIRK2 = 202;
    /** 3rd-order L-stable SDIRK method, Butcher tableau:
        <table>
            <tr>    <td>    \f$q\f$     <td>    \f$q\f$     <td>    \f$0\f$     <td>    \f$0\f$
            <tr>    <td>    \f$s\f$     <td>    \f$s-q\f$   <td>    \f$q\f$     <td>    \f$0\f$
            <tr>    <td>    \f$r\f$     <td>    \f$r\f$     <td>    \f$1-q-r\f$ <td>    \f$q\f$
            <tr>    <td>                <td>    \f$r\f$     <td>    \f$1-q-r\f$ <td>    \f$q\f$
        </table>
        with \f$q = 0.435866521508458999416019\f$, \f$r = 1.20849664917601007033648\f$ and \f$s = 0.717933260754229499708010\f$.
        See [Dobrev et al. (2017)],
        talk by <a href="https://www.math.auckland.ac.nz/~butcher/CONFERENCES/TRONDHEIM/trondheim.pdf" target="_blank">talk by Butcher</a> 
        and <a href="http://mfem.github.io/doxygen/html/ode_8cpp_source.html" target="_blank">MFEM</a>.
    */
    const int L_stable_SDIRK3 = 203;
    /** 4th-order L-stable SDIRK method, Butcher tableau:
        <table>
            <tr>    <td>    \f$0.25\f$  <td>    \f$0.25\f$      <td>    \f$0\f$         <td>    \f$0\f$         <td>    \f$0\f$         <td>    \f$0\f$
            <tr>    <td>    \f$0.75\f$  <td>    \f$0.5\f$       <td>    \f$0.25\f$      <td>    \f$0\f$         <td>    \f$0\f$         <td>    \f$0\f$
            <tr>    <td>    \f$11/20\f$ <td>    \f$17/50\f$     <td>    \f$-1/25\f$     <td>    \f$0.25\f$      <td>    \f$0\f$         <td>    \f$0\f$
            <tr>    <td>    \f$0.5\f$   <td>    \f$371/1360\f$  <td>    \f$-137/2720\f$ <td>    \f$15/544\f$    <td>    \f$0.25\f$      <td>    \f$0\f$
            <tr>    <td>    \f$1\f$     <td>    \f$25/24\f$     <td>    \f$-49/48\f$    <td>    \f$125/16\f$    <td>    \f$-85/12\f$    <td>    \f$0.25\f$
            <tr>    <td>                <td>    \f$25/24\f$     <td>    \f$-49/48\f$    <td>    \f$125/16\f$    <td>    \f$-85/12\f$    <td>    \f$0.25\f$
        </table>
        See [Hairer, Wanner (1996)], Table IV.6.5 and [Duarte, Dobbins, Smooke (2016)], Appendix C.
    */
    const int L_stable_SDIRK4 = 204;
}

#endif
