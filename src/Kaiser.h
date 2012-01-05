#ifndef KAISER_H
#define KAISER_H

#include "Common.h"

/* Formulas for Kaiser's linear theory redshift-space power spectrum, and
 * Hamilton's corresponding formulas for the redshift-space correlation
 * function. */
class Kaiser {
public:
    Kaiser(const PowerSpectrum& P_L, real f);

    real P_s(real k, real mu);

    Spline ComputeP_ell(int ell, real kmin = 1e-4, real kmax = 1e1, int Nk = 1024);
    Spline ComputeXi_ell(int ell, real kmin = 1e-4, real kmax = 1e1, int Nk = 1024);

    static void ComputeXi_ell(int ell, real f, const array& k, const array& pk, array& r, array& xi);

protected:
    const PowerSpectrum& P;
    real f;
};

#endif // KAISER_H
