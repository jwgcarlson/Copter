#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <cassert>
#include <cmath>
#include <cstdio>

#include <boost/bind.hpp>

#include "CorrelationFunction.h"
#include "PowerSpectrum.h"
#include "Quadrature.h"
#include "array.h"


CorrelationFunction::CorrelationFunction(const PowerSpectrum& P_, real kmin_, real kmax_)
    : P(P_), kmin(kmin_), kmax(kmax_)
{ }

static real f(const PowerSpectrum& P, real r, real k) {
    return k*sin(k*r)*P(k);
}

real CorrelationFunction::Evaluate(real r) const {
    return 1/(2*M_PI*M_PI*r) * Integrate<ExpSub>(bind(f, boost::cref(P), r, _1), kmin, kmax);
}

array CorrelationFunction::EvaluateMany(const array& r) const {
    int i, n = (int)r.size();
    array xi(n);
    #pragma omp parallel default(shared) private(i)
    #pragma omp for schedule(dynamic)
    for(i = 0; i < n; i++)
        xi[i] = Evaluate(r[i]);

    return xi;
}
