#ifndef SPT_H
#define SPT_H

#include "Common.h"


/******************************************************************************
 * SPT
 *
 * 1-loop standard perturbation theory.  Indices (a,b) are used to indicate
 * either the density-density, density-velocity, or velocity-velocity spectra,
 *   P_{11} <--> P_{\delta\delta}
 *   P_{12} <--> P_{\delta\theta}
 *   P_{22} <--> P_{\theta\theta}
 * whereas the labels "22" or "13" refer to different terms in the perturbation
 * series.
 ******************************************************************************/
class SPT {
public:
    SPT(const Cosmology& C, const PowerSpectrum& P_L, real epsrel = 1e-5);

    /* 1-loop power spectrum: $P_{ab}(k) = P_L(k) + P_{ab}^{(1,3)}(k) + P_{ab}^{(2,2)}(k)$ */
    real P(real k, int a = 1, int b = 1) const;

    /* $P_{ab}^{(2,2)}(k)$ */
    real P22(real k, int a = 1, int b = 1) const;

    /* $P_{ab}^{(1,3)}(k)$ */
    real P13(real k, int a = 1, int b = 1) const;

    /* $P_{\delta\delta}^{(2,2)}(k)$ */
    real P22_dd(real k) const;

    /* $P_{\delta\theta}^{(2,2)}(k)$ */
    real P22_dt(real k) const;

    /* $P_{\theta\theta}^{(2,2)}(k)$ */
    real P22_tt(real k) const;

    /* $P_{\delta\delta}^{(1,3)}(k)$ */
    real P13_dd(real k) const;

    /* $P_{\delta\theta}^{(1,3)}(k)$ */
    real P13_dt(real k) const;

    /* $P_{\theta\theta}^{(1,3)}(k)$ */
    real P13_tt(real k) const;

    /* 1-loop propagator (normalized to 1 at k = 0) */
    real G(real k) const;

private:
    const Cosmology& C;
    const PowerSpectrum& P_L;
    real epsrel;
};

#endif // SPT_H
