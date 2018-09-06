#ifndef CONSTANTS_H
#define CONSTANTS_H

// general constants
namespace constants{
    /**
     * multilevel theory assumes stable time steppers, i.e., that \f$eig(\Phi_l) < 1\f$.
     * note: we allow for some numerical error here.
     */
    const double time_stepper_stability_limit = 1.0 + 1.0e-6;
    const double pi = 3.141592653589793238462643383279502884197169399375105820974;
}

// constants defining type of estimate/bound
namespace mgritestimate{
    //----------------------------------------------------------------------------------------------------------------//
    /**
     *  F-relaxation on all levels
     */
    const int F_relaxation = 0;
    /**
     *  FCF-relaxation on all levels
     */
    const int FCF_relaxation = 1;
    //----------------------------------------------------------------------------------------------------------------//
    /**
     *  \f$\ell^2\f$-norm of numerically constructed error propagator
     */
    const int error_l2_upper_bound = 1;
    /**
     *  \f$\ell^2\f$-norm of numerically constructed error propagator is bounded by \f$\sqrt{ \| E \|_1 \| E \|_\infty }\f$
     */
    const int error_l2_sqrt_upper_bound = 2;
    /**
     *  \f$\ell^2\f$-norm of error propagator is bounded by \f$\sqrt{ \| E \|_1 \| E \|_\infty }\f$.
     *
     *  Note: Expression-like bound, faster than sqrt_upper_bound.
     */
    const int error_l2_sqrt_expression_upper_bound = 3;
    /**
     *  Tight expression-like upper bound for \f$\ell^2\f$-norm of two-grid error propagator (asymptotically exact bound)
     */
    const int error_l2_tight_twogrid_upper_bound = 4;
    /**
     *  \f$\ell^2\f$-norm of pseudo-inverse of numerically constructed error propagator
     */
    const int error_l2_lower_bound = -1;
    /**
     *  \f$\ell^2\f$-norm of numerically constructed error propagator is bounded using triangle inequality and \f$\sqrt{ \| \cdot \|_1 \| \cdot \|_\infty }\f$
     */
    const int error_l2_sqrt_lower_bound = -2;
    /**
     *  Tight expression-like lower bound for \f$\ell^2\f$-norm of two-grid error propagator (asymptotically exact bound)
     */
    const int error_l2_tight_twogrid_lower_bound = -4;
    //----------------------------------------------------------------------------------------------------------------//
    /**
     *  \f$\ell^2\f$-norm of numerically constructed residual propagator
     */
    const int residual_l2_upper_bound = 101;
    /**
     *  \f$\ell^2\f$-norm of numerically constructed residual propagator is bounded by \f$\sqrt{ \| R \|_1 \| R \|_\infty }\f$
     */
    const int residual_l2_sqrt_upper_bound = 102;
    /**
     *  \f$\ell^2\f$-norm of pseudo-inverse of numerically constructed error propagator
     */
    const int residual_l2_lower_bound = -101;
    //----------------------------------------------------------------------------------------------------------------//
}

/**
 *  Constants that refer to a specific Runge-Kutta method.
 *
 *  Note: L-stable is also A-stable.
 */
namespace rkconst{
    //-------------------------------------- A-stable methods --------------------------------------------------------//
    /** 2nd-order A-stable Lobatto IIIA method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$0\f$ <td>    \f$0\f$             <td>    \f$0\f$
     *      <tr>    <td>    \f$1\f$ <td>    \f$\frac{1}{2}\f$   <td>    \f$\frac{1}{2}\f$
     *      <tr>    <td>            <td>    \f$\frac{1}{2}\f$   <td>    \f$\frac{1}{2}\f$
     *  </table>
     *  See [Hairer, Wanner (1996)], Table 5.7 and [Jay (2015), pp. 817-826], Table 1
     */
    const int A_stable_LobattoIIIA_order2 = 402;
    /** 4th-order A-stable Lobatto IIIA method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$0\f$             <td>    \f$0\f$             <td>    \f$0\f$             <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{1}{2}\f$   <td>    \f$\frac{5}{24}\f$  <td>    \f$\frac{1}{3}\f$   <td>    \f$-\frac{1}{24}\f$
     *      <tr>    <td>    \f$1\f$             <td>    \f$\frac{1}{6}\f$   <td>    \f$\frac{2}{3}\f$   <td>    \f$\frac{1}{6}\f$
     *      <tr>    <td>                        <td>    \f$\frac{1}{6}\f$   <td>    \f$\frac{2}{3}\f$   <td>    \f$\frac{1}{6}\f$
     *  </table>
     *  See [Hairer, Wanner (1996)], Table 5.7 and [Jay (2015), pp. 817-826], Table 1
     */
    const int A_stable_LobattoIIIA_order4 = 404;
    /** 6th-order A-stable Lobatto IIIA method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$0\f$                     <td>    \f$0\f$                         <td>    \f$0\f$                         <td>    \f$0\f$                         <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{5-\sqrt{5}}{10}\f$ <td>    \f$\frac{11+\sqrt{5}}{120}\f$   <td>    \f$\frac{25-\sqrt{5}}{120}\f$   <td>    \f$\frac{25-13\sqrt{5}}{120}\f$ <td>    \f$\frac{-1+\sqrt{5}}{120}\f$
     *      <tr>    <td>    \f$\frac{5+\sqrt{5}}{10}\f$ <td>    \f$\frac{11-\sqrt{5}}{120}\f$   <td>    \f$\frac{25+13\sqrt{5}}{120}\f$ <td>    \f$\frac{25+\sqrt{5}}{120}\f$   <td>    \f$\frac{-1-\sqrt{5}}{120}\f$
     *      <tr>    <td>    \f$1\f$                     <td>    \f$\frac{1}{12}\f$              <td>    \f$\frac{5}{12}\f$              <td>    \f$\frac{5}{12}\f$              <td>    \f$\frac{1}{12}\f$
     *      <tr>    <td>                                <td>    \f$\frac{1}{12}\f$              <td>    \f$\frac{5}{12}\f$              <td>    \f$\frac{5}{12}\f$              <td>    \f$\frac{1}{12}\f$
     *  </table>
     *  See [Hairer, Wanner (1996)], Table 5.8 and [Jay (2015), pp. 817-826], Table 1
     */
    const int A_stable_LobattoIIIA_order6 = 406;
    /** 2nd-order A-stable Lobatto IIIB method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$0\f$ <td>    \f$\frac{1}{2}\f$   <td>    \f$0\f$
     *      <tr>    <td>    \f$1\f$ <td>    \f$\frac{1}{2}\f$   <td>    \f$0\f$
     *      <tr>    <td>            <td>    \f$\frac{1}{2}\f$   <td>    \f$\frac{1}{2}\f$
     *  </table>
     *  See [Jay (2015), pp. 817-826], Table 2
     */
    const int A_stable_LobattoIIIB_order2 = 502;
    /** 4th-order A-stable Lobatto IIIB method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$0\f$             <td>    \f$\frac{1}{6}\f$  <td>    \f$-\frac{1}{6}\f$ <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{1}{2}\f$   <td>    \f$\frac{1}{6}\f$  <td>    \f$\frac{1}{3}\f$  <td>    \f$0\f$
     *      <tr>    <td>    \f$1\f$             <td>    \f$\frac{1}{6}\f$  <td>    \f$\frac{5}{6}\f$  <td>    \f$0\f$
     *      <tr>    <td>                        <td>    \f$\frac{1}{6}\f$  <td>    \f$\frac{2}{3}\f$  <td>    \f$\frac{1}{6}\f$
     *  </table>
     *  See [Jay (2015), pp. 817-826], Table 2
     */
    const int A_stable_LobattoIIIB_order4 = 504;
    /** 6th-order A-stable Lobatto IIIB method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$0\f$                     <td>    \f$\frac{1}{12}\f$  <td>    \f$\frac{-1-\sqrt{5}}{24}\f$    <td>    \f$\frac{-1+\sqrt{5}}{24}\f$    <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{5-\sqrt{5}}{10}\f$ <td>    \f$\frac{1}{12}\f$  <td>    \f$\frac{25+\sqrt{5}}{120}\f$   <td>    \f$\frac{25-13\sqrt{5}}{120}\f$ <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{5+\sqrt{5}}{10}\f$ <td>    \f$\frac{1}{12}\f$  <td>    \f$\frac{25+13\sqrt{5}}{120}\f$ <td>    \f$\frac{25-\sqrt{5}}{120}\f$   <td>    \f$0\f$
     *      <tr>    <td>    \f$1\f$                     <td>    \f$\frac{1}{12}\f$  <td>    \f$\frac{11-\sqrt{5}}{24}\f$    <td>    \f$\frac{11+\sqrt{5}}{24}\f$    <td>    \f$0\f$
     *      <tr>    <td>                                <td>    \f$\frac{1}{12}\f$  <td>    \f$\frac{5}{12}\f$              <td>    \f$\frac{5}{12}\f$              <td>    \f$\frac{1}{12}\f$
     *  </table>
     *  See [Jay (2015), pp. 817-826], Table 2
     */
    const int A_stable_LobattoIIIB_order6 = 506;
    /** 8th-order A-stable Lobatto IIIB method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$0\f$                         <td>    \f$\frac{1}{20}\f$  <td>    \f$\frac{-7-\sqrt{21}}{120}\f$      <td>    \f$\frac{1}{15}\f$                  <td>    \f$\frac{-7+\sqrt{21}}{120}\f$      <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{7-\sqrt{21}}{14}\f$    <td>    \f$\frac{1}{20}\f$  <td>    \f$\frac{343+9\sqrt{21}}{2520}\f$   <td>    \f$\frac{56-15\sqrt{21}}{315}\f$    <td>    \f$\frac{343-69\sqrt{21}}{2520}\f$  <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{1}{2}\f$               <td>    \f$\frac{1}{20}\f$  <td>    \f$\frac{49+12\sqrt{21}}{360}\f$    <td>    \f$\frac{8}{45}\f$                  <td>    \f$\frac{49-12\sqrt{21}}{360}\f$    <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{7+\sqrt{21}}{14}\f$    <td>    \f$\frac{1}{20}\f$  <td>    \f$\frac{343+69\sqrt{21}}{2520}\f$  <td>    \f$\frac{56+15\sqrt{21}}{315}\f$    <td>    \f$\frac{343-9\sqrt{21}}{2520}\f$   <td>    \f$0\f$
     *      <tr>    <td>    \f$1\f$                         <td>    \f$\frac{1}{20}\f$  <td>    \f$\frac{119-3\sqrt{21}}{360}\f$    <td>    \f$\frac{13}{45}\f$                 <td>    \f$\frac{119+3\sqrt{21}}{360}\f$    <td>    \f$0\f$
     *      <tr>    <td>                                    <td>    \f$\frac{1}{20}\f$  <td>    \f$\frac{49}{180}\f$                <td>    \f$\frac{16}{45}\f$                 <td>    \f$\frac{49}{180}\f$                <td>    \f$\frac{1}{20}\f$
     * </table>
     *  See [Jay (2015), pp. 817-826], Table 2
     */
    const int A_stable_LobattoIIIB_order8 = 508;
    /** 2nd-order A-stable SDIRK method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$\frac{1}{4}\f$   <td>    \f$\frac{1}{4}\f$   <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{3}{4}\f$   <td>    \f$\frac{1}{2}\f$   <td>    \f$\frac{1}{4}\f$
     *      <tr>    <td>                        <td>    \f$\frac{1}{2}\f$   <td>    \f$\frac{1}{2}\f$
     *  </table>
     *  See [Bonaventura, Della Rocca (2015)], Section 5.1.
     */
    const int A_stable_SDIRK2 = 102;
    /** 3rd-order A-stable SDIRK method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$\gamma\f$    <td>    \f$\gamma\f$        <td>    \f$0\f$
     *      <tr>    <td>    \f$1-\gamma\f$  <td>    \f$1-2\gamma\f$     <td>    \f$\gamma\f$
     *      <tr>    <td>                    <td>    \f$\frac{1}{2}\f$   <td>    \f$\frac{1}{2}\f$
     *  </table>
     *  See <a href="http://mfem.github.io/doxygen/html/ode_8cpp_source.html" target="_blank">MFEM</a>.
     */
    const int A_stable_SDIRK3 = 103;
    /** 4th-order A-stable SDIRK method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$q\f$             <td>    \f$q\f$             <td>    \f$0\f$     <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{1}{2}\f$   <td>    \f$\frac{1}{2}-q\f$ <td>    \f$q\f$     <td>    \f$0\f$
     *      <tr>    <td>    \f$1-q\f$           <td>    \f$2q\f$            <td>    \f$1-4q\f$  <td>    \f$q\f$
     *      <tr>    <td>                        <td>    \f$2q\f$            <td>    \f$1-4q\f$  <td>    \f$q\f$
     *  </table>
     *  See <a href="http://mfem.github.io/doxygen/html/ode_8cpp_source.html" target="_blank">MFEM</a>.
     */
    const int A_stable_SDIRK4 = 104;
    //-------------------------------------- L-stable methods --------------------------------------------------------//
    /** 1st-order L-stable SDIRK method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$1\f$ <td>    \f$1\f$
     *      <tr>    <td>            <td>    \f$1\f$
     *  </table>
     *  See [Dobrev et al. (2017)], Table 3.1.
     */
    const int L_stable_SDIRK1 = 201;
    /** 2nd-order L-stable SDIRK method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$\gamma\f$    <td>    \f$\gamma\f$   <td>    \f$0\f$
     *      <tr>    <td>    \f$1\f$         <td>    \f$1-\gamma\f$ <td>    \f$\gamma\f$
     *      <tr>    <td>                    <td>    \f$1-\gamma\f$ <td>    \f$\gamma\f$
     *  </table>
     *  with \f$\gamma = 1.0 - 1.0 / \sqrt{2.0}\f$.
     *  See [Dobrev et al. (2017)], Table 3.1.
     */
    const int L_stable_SDIRK2 = 202;
    /** 3rd-order L-stable SDIRK method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$q\f$     <td>    \f$q\f$     <td>    \f$0\f$     <td>    \f$0\f$
     *      <tr>    <td>    \f$s\f$     <td>    \f$s-q\f$   <td>    \f$q\f$     <td>    \f$0\f$
     *      <tr>    <td>    \f$r\f$     <td>    \f$r\f$     <td>    \f$1-q-r\f$ <td>    \f$q\f$
     *      <tr>    <td>                <td>    \f$r\f$     <td>    \f$1-q-r\f$ <td>    \f$q\f$
     *  </table>
     *  with \f$q = 0.435866521508458999416019\f$, \f$r = 1.20849664917601007033648\f$ and \f$s = 0.717933260754229499708010\f$.
     *  See [Dobrev et al. (2017)], Table 3.1,
     *  talk by <a href="https://www.math.auckland.ac.nz/~butcher/CONFERENCES/TRONDHEIM/trondheim.pdf" target="_blank">talk by Butcher</a> 
     *  and <a href="http://mfem.github.io/doxygen/html/ode_8cpp_source.html" target="_blank">MFEM</a>.
     */
    const int L_stable_SDIRK3 = 203;
    /** 4th-order L-stable SDIRK method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$0.25\f$  <td>    \f$0.25\f$      <td>    \f$0\f$         <td>    \f$0\f$         <td>    \f$0\f$         <td>    \f$0\f$
     *      <tr>    <td>    \f$0.75\f$  <td>    \f$0.5\f$       <td>    \f$0.25\f$      <td>    \f$0\f$         <td>    \f$0\f$         <td>    \f$0\f$
     *      <tr>    <td>    \f$11/20\f$ <td>    \f$17/50\f$     <td>    \f$-1/25\f$     <td>    \f$0.25\f$      <td>    \f$0\f$         <td>    \f$0\f$
     *      <tr>    <td>    \f$0.5\f$   <td>    \f$371/1360\f$  <td>    \f$-137/2720\f$ <td>    \f$15/544\f$    <td>    \f$0.25\f$      <td>    \f$0\f$
     *      <tr>    <td>    \f$1\f$     <td>    \f$25/24\f$     <td>    \f$-49/48\f$    <td>    \f$125/16\f$    <td>    \f$-85/12\f$    <td>    \f$0.25\f$
     *      <tr>    <td>                <td>    \f$25/24\f$     <td>    \f$-49/48\f$    <td>    \f$125/16\f$    <td>    \f$-85/12\f$    <td>    \f$0.25\f$
     *  </table>
     *  See [Hairer, Wanner (1996)], Table IV.6.5 and [Duarte, Dobbins, Smooke (2016)], Appendix C.
     */
    const int L_stable_SDIRK4 = 204;
    /** 5th-order L-stable SDIRK method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$\frac{4024571134387}{14474071345096}\f$  <td>    \f$\frac{4024571134387}{14474071345096}\f$  <td>    \f$0\f$                                     <td>    \f$0\f$                                     <td>    \f$0\f$                                     <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{5555633399575}{5431021154178}\f$   <td>    \f$\frac{9365021263232}{12572342979331}\f$  <td>    \f$\frac{4024571134387}{14474071345096}\f$  <td>    \f$0\f$                                     <td>    \f$0\f$                                     <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{5255299487392}{12852514622453}\f$  <td>    \f$\frac{2144716224527}{9320917548702}\f$   <td>    \f$-\frac{397905335951}{4008788611757}\f$   <td>    \f$\frac{4024571134387}{14474071345096}\f$  <td>    \f$0\f$                                     <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{3}{20}\f$                          <td>    \f$-\frac{291541413000}{6267936762551}\f$   <td>    \f$\frac{226761949132}{4473940808273}\f$    <td>    \f$-\frac{1282248297070}{9697416712681}\f$  <td>    \f$\frac{4024571134387}{14474071345096}\f$  <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{10449500210709}{14474071345096}\f$ <td>    \f$-\frac{2481679516057}{4626464057815}\f$  <td>    \f$-\frac{197112422687}{6604378783090}\f$   <td>    \f$\frac{3952887910906}{9713059315593}\f$   <td>    \f$\frac{4906835613583}{8134926921134}\f$   <td>    \f$\frac{4024571134387}{14474071345096}\f$
     *      <tr>    <td>                                                <td>    \f$-\frac{2522702558582}{12162329469185}\f$ <td>    \f$\frac{1018267903655}{12907234417901}\f$  <td>    \f$\frac{4542392826351}{13702606430957}\f$  <td>    \f$\frac{5001116467727}{12224457745473}\f$  <td>    \f$\frac{1509636094297}{3891594770934}\f$
     *  </table>
     *  See [Kennedy, Carpenter (2016)], NASA/TM-2016-219173, Table 24.
     */
    const int L_stable_SDIRK5 = 205;
    /** 1st-order L-stable Radau IIA method (same as rkconst::L_stable_SDIRK1), Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$1\f$ <td>    \f$1\f$
     *      <tr>    <td>            <td>    \f$1\f$
     *  </table>
     *  See [Hairer, Wanner (1996)], Table 5.5
     */
    const int L_stable_RadauIIA_order1 = 301;
    /** 3rd-order L-stable Radau IIA method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$\frac{1}{3}\f$   <td>    \f$\frac{5}{12}\f$  <td>    \f$-\frac{1}{12}\f$
     *      <tr>    <td>    \f$1\f$             <td>    \f$\frac{3}{4}\f$   <td>    \f$\frac{1}{4}\f$
     *      <tr>    <td>                        <td>    \f$\frac{3}{4}\f$   <td>    \f$\frac{1}{4}\f$
     *  </table>
     *  See [Hairer, Wanner (1996)], Table 5.5
    */
    const int L_stable_RadauIIA_order3 = 303;
    /** 5th-order L-stable Radau IIA method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$\frac{4-\sqrt{6}}{10}\f$ <td>    \f$\frac{88-7\sqrt{6}}{360}\f$      <td>    \f$\frac{296-169\sqrt{6}}{1800}\f$  <td>    \f$\frac{-2+3\sqrt{6}}{225}\f$
     *      <tr>    <td>    \f$\frac{4+\sqrt{6}}{10}\f$ <td>    \f$\frac{296+169\sqrt{6}}{1800}\f$  <td>    \f$\frac{88+7\sqrt{6}}{360}\f$      <td>    \f$\frac{-2-3\sqrt{6}}{225}\f$
     *      <tr>    <td>    \f$1\f$                     <td>    \f$\frac{16-\sqrt{6}}{36}\f$        <td>    \f$\frac{16+\sqrt{6}}{36}\f$        <td>    \f$\frac{1}{9}\f$
     *      <tr>    <td>                                <td>    \f$\frac{16-\sqrt{6}}{36}\f$        <td>    \f$\frac{16+\sqrt{6}}{36}\f$        <td>    \f$\frac{1}{9}\f$
     *  </table>
     *  See [Hairer, Norsett, Wanner (1987)], Table 7.7 and [Hairer, Wanner (1996)], Table 5.6
    */
    const int L_stable_RadauIIA_order5 = 305;
    /** 2nd-order L- and B-stable Lobatto IIIC method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$0\f$   <td>    \f$\frac{1}{2}\f$    <td>    \f$-\frac{1}{2}\f$
     *      <tr>    <td>    \f$1\f$   <td>    \f$\frac{1}{2}\f$    <td>    \f$\frac{1}{2}\f$
     *      <tr>    <td>              <td>    \f$\frac{1}{2}\f$    <td>    \f$\frac{1}{2}\f$
     *  </table>
     *  See [Jay (2015), pp. 817-826], Table 3
     */
    const int L_stable_LobattoIIIC_order2 = 602;
    /** 4th-order L- and B-stable Lobatto IIIC method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$0\f$             <td>    \f$\frac{1}{6}\f$   <td>    \f$-\frac{1}{3}\f$  <td>    \f$\frac{1}{6}\f$
     *      <tr>    <td>    \f$\frac{1}{2}\f$   <td>    \f$\frac{1}{6}\f$   <td>    \f$\frac{5}{12}\f$  <td>    \f$-\frac{1}{12}\f$
     *      <tr>    <td>    \f$1\f$             <td>    \f$\frac{1}{6}\f$   <td>    \f$\frac{2}{3}\f$   <td>    \f$\frac{1}{6}\f$
     *      <tr>    <td>                        <td>    \f$\frac{1}{6}\f$   <td>    \f$\frac{2}{3}\f$   <td>    \f$\frac{1}{6}\f$
     *  </table>
     *  See [Jay (2015), pp. 817-826], Table 3
     */
    const int L_stable_LobattoIIIC_order4 = 604;
    /** 6th-order L- and B-stable Lobatto IIIC method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$0\f$                     <td>    \f$\frac{1}{12}\f$  <td>    \f$-\frac{\sqrt{5}}{12}\f$      <td>    \f$\frac{\sqrt{5}}{12}\f$       <td>    \f$-\frac{1}{12}\f$
     *      <tr>    <td>    \f$\frac{5-\sqrt{5}}{10}\f$ <td>    \f$\frac{1}{12}\f$  <td>    \f$\frac{1}{4}\f$               <td>    \f$\frac{10-7\sqrt{5}}{60}\f$   <td>    \f$\frac{\sqrt{5}}{60}\f$
     *      <tr>    <td>    \f$\frac{5+\sqrt{5}}{10}\f$ <td>    \f$\frac{1}{12}\f$  <td>    \f$\frac{10+7\sqrt{5}}{60}\f$   <td>    \f$\frac{1}{4}\f$               <td>    \f$-\frac{\sqrt{5}}{60}\f$
     *      <tr>    <td>    \f$1\f$                     <td>    \f$\frac{1}{12}\f$  <td>    \f$\frac{5}{12}\f$              <td>    \f$\frac{5}{12}\f$              <td>    \f$\frac{1}{12}\f$
     *      <tr>    <td>                                <td>    \f$\frac{1}{12}\f$  <td>    \f$\frac{5}{12}\f$              <td>    \f$\frac{5}{12}\f$              <td>    \f$\frac{1}{12}\f$
     *  </table>
     *  See [Jay (2015), pp. 817-826], Table 3
     */
    const int L_stable_LobattoIIIC_order6 = 606;
    /** 8th-order L- and B-stable Lobatto IIIC method, Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$0\f$                         <td>    \f$\frac{1}{20}\f$    <td>    \f$-\frac{7}{60}\f$                   <td>    \f$\frac{2}{15}\f$                  <td>    \f$-\frac{7}{60}\f$                 <td>    \f$\frac{1}{20}\f$
     *      <tr>    <td>    \f$\frac{7-\sqrt{21}}{14}\f$    <td>    \f$\frac{1}{20}\f$    <td>    \f$\frac{29}{180}\f$                  <td>    \f$\frac{47-15\sqrt{21}}{315}\f$    <td>    \f$\frac{203-30\sqrt{21}}{1260}\f$  <td>    \f$-\frac{3}{140}\f$
     *      <tr>    <td>    \f$\frac{1}{2}\f$               <td>    \f$\frac{1}{20}\f$    <td>    \f$\frac{329+105\sqrt{21}}{2880}\f$   <td>    \f$\frac{73}{360}\f$                <td>    \f$\frac{329-105\sqrt{21}}{2880}\f$ <td>    \f$\frac{3}{160}\f$
     *      <tr>    <td>    \f$\frac{7+\sqrt{21}}{14}\f$    <td>    \f$\frac{1}{20}\f$    <td>    \f$\frac{203+30\sqrt{21}}{1260}\f$    <td>    \f$\frac{47+15\sqrt{21}}{315}\f$    <td>    \f$\frac{29}{180}\f$                <td>    \f$-\frac{3}{140}\f$
     *      <tr>    <td>    \f$1\f$                         <td>    \f$\frac{1}{20}\f$    <td>    \f$\frac{49}{180}\f$                  <td>    \f$\frac{16}{45}\f$                 <td>    \f$\frac{49}{180}\f$                <td>    \f$\frac{1}{20}\f$
     *      <tr>    <td>                                    <td>    \f$\frac{1}{20}\f$    <td>    \f$\frac{49}{180}\f$                  <td>    \f$\frac{16}{45}\f$                 <td>    \f$\frac{49}{180}\f$                <td>    \f$\frac{1}{20}\f$
     *  </table>
     *  See [Jay (2015), pp. 817-826], Table 3
     */
    const int L_stable_LobattoIIIC_order8 = 608;
    //-------------------------------------- not A-stable methods ----------------------------------------------------//
    /** 2nd-order Lobatto IIIC\* method (not A-stable, not B-stable), Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$0\f$ <td>    \f$0\f$             <td>    \f$0\f$
     *      <tr>    <td>    \f$1\f$ <td>    \f$1\f$             <td>    \f$0\f$
     *      <tr>    <td>            <td>    \f$\frac{1}{2}\f$   <td>    \f$\frac{1}{2}\f$
     *  </table>
     *  See [Jay (2015), pp. 817-826], Table 4
     */
    const int LobattoIIICast_order2 = 702;
    /** 4th-order Lobatto IIIC\* method (not A-stable, not B-stable), Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$0\f$             <td>    \f$0\f$             <td>    \f$0\f$             <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{1}{2}\f$   <td>    \f$\frac{1}{4}\f$   <td>    \f$\frac{1}{4}\f$   <td>    \f$0\f$
     *      <tr>    <td>    \f$1\f$             <td>    \f$0\f$             <td>    \f$1\f$             <td>    \f$0\f$
     *      <tr>    <td>                        <td>    \f$\frac{1}{6}\f$   <td>    \f$\frac{2}{3}\f$   <td>    \f$\frac{1}{6}\f$
     *  </table>
     *  See [Jay (2015), pp. 817-826], Table 4
     */
    const int LobattoIIICast_order4 = 704;
    /** 6th-order Lobatto IIIC\* method (not A-stable, not B-stable), Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$0\f$                     <td>    \f$0\f$                     <td>    \f$0\f$                         <td>    \f$0\f$                         <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{5-\sqrt{5}}{10}\f$ <td>    \f$\frac{5+\sqrt{5}}{60}\f$ <td>    \f$\frac{1}{6}\f$               <td>    \f$\frac{15-7\sqrt{5}}{60}\f$   <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{5+\sqrt{5}}{10}\f$ <td>    \f$\frac{5-\sqrt{5}}{60}\f$ <td>    \f$\frac{15+7\sqrt{5}}{60}\f$   <td>    \f$\frac{1}{6}\f$               <td>    \f$0\f$
     *      <tr>    <td>    \f$1\f$                     <td>    \f$\frac{1}{6}\f$           <td>    \f$\frac{5-\sqrt{5}}{12}\f$     <td>    \f$\frac{5+\sqrt{5}}{12}\f$     <td>    \f$0\f$
     *      <tr>    <td>                                <td>    \f$\frac{1}{12}\f$          <td>    \f$\frac{5}{12}\f$              <td>    \f$\frac{5}{12}\f$              <td>    \f$\frac{1}{12}\f$
     *  </table>
     *  See [Jay (2015), pp. 817-826], Table 4
     */
    const int LobattoIIICast_order6 = 706;
    /** 8th-order Lobatto IIIC\* method (not A-stable, not B-stable), Butcher tableau:
     *  <table>
     *      <tr>    <td>    \f$0\f$                         <td>    \f$0\f$             <td>    \f$0\f$                             <td>    \f$0\f$                         <td>    \f$0\f$                             <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{7-\sqrt{21}}{14}\f$    <td>    \f$\frac{1}{14}\f$  <td>    \f$\frac{1}{9}\f$                   <td>    \f$\frac{13-3\sqrt{21}}{63}\f$  <td>    \f$\frac{14-3\sqrt{21}}{126}\f$     <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{1}{2}\f$               <td>    \f$\frac{1}{32}\f$  <td>    \f$\frac{91+21\sqrt{21}}{576}\f$    <td>    \f$\frac{11}{72}\f$             <td>    \f$\frac{91-21\sqrt{21}}{576}\f$    <td>    \f$0\f$
     *      <tr>    <td>    \f$\frac{7+\sqrt{21}}{14}\f$    <td>    \f$\frac{1}{14}\f$  <td>    \f$\frac{14+3\sqrt{21}}{126}\f$     <td>    \f$\frac{13+3\sqrt{21}}{63}\f$  <td>    \f$\frac{1}{9}\f$                   <td>    \f$0\f$
     *      <tr>    <td>    \f$1\f$                         <td>    \f$0\f$             <td>    \f$\frac{7}{18}\f$                  <td>    \f$\frac{2}{9}\f$               <td>    \f$\frac{7}{18}\f$                  <td>    \f$0\f$
     *      <tr>    <td>                                    <td>    \f$\frac{1}{20}\f$  <td>    \f$\frac{49}{180}\f$                <td>    \f$\frac{16}{45}\f$             <td>    \f$\frac{49}{180}\f$                <td>    \f$\frac{1}{20}\f$
     *  </table>
     *  See [Jay (2015), pp. 817-826], Table 4
     */
    const int LobattoIIICast_order8 = 708;
}

#endif
