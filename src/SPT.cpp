#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/bind.hpp>
using boost::cref;

#include "Cosmology.h"
#include "PowerSpectrum.h"
#include "Quadrature.h"
#include "SPT.h"

/* Limits of integration for second-order power spectrum */
const real QMIN = 1e-5;
const real QMAX = 1e2;


SPT::SPT(const Cosmology& C_, const PowerSpectrum& P_L_, real epsrel_)
    : C(C_), P_L(P_L_)
{
    epsrel = epsrel_;
}

/* P_{ab}(k) */
real SPT::P(real k, int a, int b) const {
    switch(a*b) {
        case 1:
            return P_L(k) + P13_dd(k) + P22_dd(k);
        case 2:
            return P_L(k) + P13_dt(k) + P22_dt(k);
        case 4:
            return P_L(k) + P13_tt(k) + P22_tt(k);
        default:
            warning("SPT: invalid indices, a = %d, b = %d\n", a, b);
            return 0;
    }
}

/* P_{ab}^{(22)} */
real SPT::P22(real k, int a, int b) const {
    switch(a*b) {
        case 1:
            return P22_dd(k);
        case 2:
            return P22_dt(k);
        case 4:
            return P22_tt(k);
        default:
            warning("SPT: invalid indices, a = %d, b = %d\n", a, b);
            return 0;
    }
}

/* P_{ab}^{(13)} */
real SPT::P13(real k, int a, int b) const {
    switch(a*b) {
        case 1: // 1,1
            return P13_dd(k);
        case 2: // 1,2 or 2,1
            return P13_dt(k);
        case 4: // 2,2
            return P13_tt(k);
        default:
            warning("SPT: invalid indices, a = %d, b = %d\n", a, b);
            return 0;
    }
}

static real f22_dd(const PowerSpectrum& P_L, real k, real logq, real x) {
    double q = exp(logq), r = q/k, r2 = r*r, rx = r*x, x2 = x*x, d = 1 + r2 - 2*rx;
    if(d < 1e-10)
        return 0;
    else
        return q * P_L(q) * P_L(k*sqrt(d)) * pow2(3*r + 7*x - 10*r*x2) / pow2(d);
}

/* P_{\delta\delta}^{(22)} */
real SPT::P22_dd(real k) const {
    if(k <= 0)
        return 0;

    real a[2] = { log(QMIN), -1 };
    real b[2] = { log(QMAX), +1 };
    real V = k*k/(98*4*M_PI*M_PI);
    return V * Integrate<2>(bind(f22_dd, cref(P_L), k, _1, _2), a, b, epsrel, epsrel*P_L(k)/V);
}

static real f13_dd(const PowerSpectrum& P_L, real k, real logq) {
    real q = exp(logq), r = q/k;
    real s;
    if(r < 1e-2)
        s = -168 + (928./5.)*pow2(r) - (4512./35.)*pow4(r) + (416./21.)*pow6(r);
    else if(fabs(r-1) < 1e-10)
        s = -88 + 8*(r-1);
    else if(r > 100)
        s = -488./5. + (96./5.)/pow2(r) - (160./21.)/pow4(r) - (1376./1155.)/pow6(r);
    else
        s = 12/pow2(r) - 158 + 100*pow2(r) - 42*pow4(r) + 3/pow3(r) * pow3(r*r - 1) * (7*r*r + 2) * log((1+r)/fabs(1-r));

    return q * P_L(q) * s;
}

/* P_{\delta\delta}^{(13)} */
real SPT::P13_dd(real k) const {
    real a = log(QMIN);
    real b = log(QMAX);
    real V = k*k/(252*4*M_PI*M_PI) * P_L(k);
    return V * Integrate(bind(f13_dd, cref(P_L), k, _1), a, b, epsrel, epsrel*P_L(k)/V);
}

static real f22_dt(const PowerSpectrum& P_L, real k, real logq, real x) {
    double q = exp(logq), r = q/k, r2 = r*r, rx = r*x, x2 = x*x, d = 1 + r2 - 2*rx;
    if(d < 1e-10)
        return 0;
    else
        return q * P_L(q) * P_L(k*sqrt(d)) * (3*r + 7*x - 10*r*x*x)*(7*x - r - 6*r*x*x) / pow2(d);
}

/* P_{\delta\theta}^{(22)} */
real SPT::P22_dt(real k) const {
    if(k <= 0)
        return 0;

    real a[2] = { log(QMIN), -1 };
    real b[2] = { log(QMAX), +1 };
    real V = k*k/(98*4*M_PI*M_PI);
    return V * Integrate<2>(bind(f22_dt, cref(P_L), k, _1, _2), a, b, epsrel, epsrel*P_L(k)/V);
}

static real f13_dt(const PowerSpectrum& P_L, real k, real logq) {
    real q = exp(logq), r = q/k;
    real s;
    if(r < 1e-2)
        s = -168 + (416./5.)*pow2(r) - (2976./35.)*pow4(r) + (224./15.)*pow6(r);
    else if(fabs(r-1) < 1e-10)
        s = -152 - 56*(r-1);
    else if(r > 100)
        s = -200 + (2208./35.)/pow2(r) - (1312./105.)/pow4(r) - (1888./1155.)/pow6(r);
    else
        s = 24/pow2(r) - 202 + 56*pow2(r) - 30*pow4(r) + 3/pow3(r) * pow3(r*r - 1) * (5*r*r + 4) * log((1+r)/fabs(1-r));

    return q * P_L(q) * s;
}

/* P_{\delta\theta}^{(13)} */
real SPT::P13_dt(real k) const {
    real a = log(QMIN);
    real b = log(QMAX);
    real V = k*k/(252*4*M_PI*M_PI) * P_L(k);
    return V * Integrate(bind(f13_dt, cref(P_L), k, _1), a, b, epsrel, epsrel*P_L(k)/V);
}

static real f22_tt(const PowerSpectrum& P_L, real k, real logq, real x) {
    double q = exp(logq), r = q/k, r2 = r*r, rx = r*x, x2 = x*x, d = 1 + r2 - 2*rx;
    if(d < 1e-10)
        return 0;
    else
        return q * P_L(q) * P_L(k*sqrt(d)) * pow2(7*x - r - 6*r*x*x) / pow2(d);
}

/* P_{\theta\theta}^{(22)} */
real SPT::P22_tt(real k) const {
    if(k <= 0)
        return 0;

    real a[2] = { log(QMIN), -1 };
    real b[2] = { log(QMAX), +1 };
    real V = k*k/(98*4*M_PI*M_PI);
    return V * Integrate<2>(bind(f22_tt, cref(P_L), k, _1, _2), a, b, epsrel, epsrel*P_L(k)/V);
}

static real f13_tt(const PowerSpectrum& P_L, real k, real logq) {
    real q = exp(logq), r = q/k;
    real s;
    if(r < 1e-2)
        s = -56 - (32./5.)*pow2(r) - (96./7.)*pow4(r) + (352./105.)*pow6(r);
    else if(fabs(r-1) < 1e-10)
        s = -72 - 40*(r-1);
    else if(r > 100)
        s = -504./5. + (1248./35.)/pow2(r) - (608./105.)/pow4(r) - (160./231.)/pow6(r);
    else
        s = 12/pow2(r) - 82 + 4*pow2(r) - 6*pow4(r) + 3/pow3(r) * pow3(r*r - 1) * (r*r + 2) * log((1+r)/fabs(1-r));

    return q * P_L(q) * s;
}

/* P_{\theta\theta}^{(13)} */
real SPT::P13_tt(real k) const {
    real a = log(QMIN);
    real b = log(QMAX);
    real V = k*k/(84*4*M_PI*M_PI) * P_L(k);
    return V * Integrate(bind(f13_tt, cref(P_L), k, _1), a, b, epsrel, epsrel*P_L(k)/V);
}

real SPT::G(real k) const {
    return 1 + 0.5*P13_dd(k)/P_L(k);
}
