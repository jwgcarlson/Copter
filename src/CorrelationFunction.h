#ifndef CORRELATION_FUNCTION_H
#define CORRELATION_FUNCTION_H

#include "Common.h"
#include "array.h"

class CorrelationFunction {
public:
    CorrelationFunction(const PowerSpectrum& P, real kmin = 1e-4, real kmax = 1e1);

    real Evaluate(real r) const;
    real operator()(real r) const { return Evaluate(r); }

    array EvaluateMany(const array& r) const;
    array operator()(const array& r) const { return EvaluateMany(r); }

    const PowerSpectrum& GetPowerSpectrum() const { return P; }

protected:
    const PowerSpectrum& P;
    real kmin, kmax;
};

#endif // CORRELATION_FUNCTION_H
