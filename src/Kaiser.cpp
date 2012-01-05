#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "Cosmology.h"
#include "Kaiser.h"
#include "PowerSpectrum.h"
#include "Spline.h"
#include "fftlog.h"


Kaiser::Kaiser(const PowerSpectrum& P_, real f_)
    : P(P_), f(f_)
{
}

real Kaiser::P_s(real k, real mu) {
    return pow2(1 + f*mu*mu) * P(k);
}

Spline Kaiser::ComputeP_ell(int ell, real kmin, real kmax, int Nk) {
    real prefactor;
    if(ell == 0)
        prefactor = 1. + (2./3.)*f + (1./5.)*f*f;
    else if(ell == 2)
        prefactor = (4./3.)*f + (4./7.)*f*f;
    else if(ell == 4)
        prefactor = (8./35.)*f*f;
    else
        prefactor = 0.;

    array k = array::logspace(kmin, kmax, Nk);
    array pk = prefactor * P(k);
    return CubicSpline(k, pk);
}

Spline Kaiser::ComputeXi_ell(int ell, real kmin, real kmax, int N) {
    real sign = (ell % 2 == 1) ? 0. : ((ell % 4 == 0) ? +1. : -1.);

    real prefactor = sign;
    if(ell == 0)
        prefactor *= 1. + (2./3.)*f + (1./5.)*f*f;
    else if(ell == 2)
        prefactor *= (4./3.)*f + (4./7.)*f*f;
    else if(ell == 4)
        prefactor *= (8./35.)*f*f;

    array k = array::logspace(kmin, kmax, N);
    array pk = prefactor * P(k);
    array r(N), xi(N);
    ComputeXiLM(ell, 2, N, k, pk, r, xi);
    return CubicSpline(r, xi);
}

void Kaiser::ComputeXi_ell(int ell, real f, const array& k, const array& pk, array& r, array& xi) {
    int N = (int)k.size();
    real sign = (ell % 2 == 1) ? 0. : ((ell % 4 == 0) ? +1. : -1.);
    real prefactor = sign;
    if(ell == 0)
        prefactor *= 1. + (2./3.)*f + (1./5.)*f*f;
    else if(ell == 2)
        prefactor *= (4./3.)*f + (4./7.)*f*f;
    else if(ell == 4)
        prefactor *= (8./35.)*f*f;

    array pks = prefactor * pk;
    r.resize(N);
    xi.resize(N);
    ComputeXiLM(ell, 2, N, k, pks, r, xi);
}
