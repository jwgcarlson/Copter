#include <cmath>
#include <cstdio>

#include <boost/bind.hpp>
using boost::bind;

#include <Copter/Quadrature.h>

const double infty = 1e100;


double f(double x) {
    return x*exp(-x);
}

double g(double s, double t) {
    return pow(t, s-1) * exp(-t);
}

double gamma(double s, double x) {
    return Integrate(bind(g, s, _1), 0, x);
}

double Gamma(double s, double x) {
    if(x > 0)
        return Integrate<InverseSub>(bind(g, s, _1), x, infty);
    else
        return 0;
}

int main() {
    printf("I = %g\n", Integrate(f, 0, 1));
    printf("gamma(4, 10) = %g\n", gamma(4, 10));
    printf("Gamma(4, 10) = %g\n", Gamma(4, 10));
    printf("sum = %g\n", gamma(4, 10) + Gamma(4, 10));
}
