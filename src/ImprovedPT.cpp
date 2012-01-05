#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/bind.hpp>
using boost::cref;

#include "ImprovedPT.h"
#include "Quadrature.h"
#include "SpecialFunctions.h"

/* Limits and tolerance for scale-dependent factors (f,g,h,i) in nonlinear propagator */
const real QMIN = 1e-4;
const real QMAX = 1e1;
const real EPS = 1e-5;


ImprovedPT::ImprovedPT(const Cosmology& C_, const PowerSpectrum& P_i_, real z_i_, real z_, int Neta_, int Nk_, real kcut_)
    : C(C_), P_i(P_i_)
{
    z_i = z_i_;
    z = z_;
    Neta = Neta_;
    Nk = Nk_;
    kcut = kcut_;

    D = GrowthFunction(C);
    deta = log(D(z)/D(z_i))/(Neta-1);

    PrecomputePropagator();
}

ImprovedPT::~ImprovedPT() {
}

real ImprovedPT::U(real X) {
    real x = sqrt(2*fabs(X));
    if(x < 0.01)
        return 1 - pow2(x)/2 + pow4(x)/12 - pow6(x)/144;
    else
        return BesselJ1(2*x)/x;    // cylindrical Bessel function
}


/***************************
 * Renormalized propagator *
 ***************************/

/* Time-dependent factors */
real ImprovedPT::G_alpha(real eta, real etap) const {
    return exp(2*etap) * (exp(2*(eta-etap)) - 1.4*exp(eta-etap) + 0.4*exp(-1.5*(eta-etap)));
}
real ImprovedPT::G_beta_g(real eta, real etap) const {
    return exp(2*etap) * (0.6*exp(eta-etap) - 1 + 0.4*exp(-1.5*(eta-etap)));
}
real ImprovedPT::G_beta_d(real eta, real etap) const {
    return exp(2*etap) * (0.6*exp(3.5*(eta-etap)) - exp(2.5*(eta-etap)) + 0.4*exp(eta-etap));
}
real ImprovedPT::G_gamma_g(real eta, real etap) const {
    return exp(2*etap) * (0.4*exp(eta-etap) - exp(-0.5*(eta-etap)) + 0.6*exp(-1.5*(eta-etap)));
}
real ImprovedPT::G_gamma_d(real eta, real etap) const {
    return exp(2*etap) * (0.4*exp(3.5*(eta-etap)) - exp(2*(eta-etap)) + 0.6*exp(eta-etap));
}
real ImprovedPT::G_delta(real eta, real etap) const {
    return exp(2*etap) * (0.4*exp(3.5*(eta-etap)) - 1.4*exp(eta-etap) + 1);
}

/* k-dependent integrands */
static real B_f(const PowerSpectrum& P, real k, real q) {
    real r = q/k;
    real r2 = pow2(r), r4 = pow4(r);
    if(r < 1e-4)
        return P(q) * (-84 + 464./5.*r2 - 2256./35.*r4);
    else if(fabs(1 - r) < 1e-8)
        return P(q) * (-44 - 4*(1-r));
    else if(r > 10)
        return P(q) * (-244./5. + 48./5/r2 - 80./21./r4);
    else
        return P(q) * (6/r2 - 79 + 50*r2 - 21*r4 + 0.75/pow3(r)*pow3(1-r2)*(2+7*r2)*log(pow2((1-r)/(1+r))));
}
static real B_g(const PowerSpectrum& P, real k, real q) {
    real r = q/k;
    real r2 = pow2(r), r4 = pow4(r);
    if(r < 1e-4)
        return P(q) * (-28 - 16./5.*r2 - 48./7.*r4);
    else if(fabs(1 - r) < 1e-8)
        return P(q) * (-36 + 20*(1-r));
    else if(r > 10)
        return P(q) * (-252./5. + 624./35/r2 - 304./105./r4);
    else
        return P(q) * (6/r2 - 41 + 2*r2 - 3*r4 + 0.75/pow3(r)*pow3(1-r2)*(2+r2)*log(pow2((1-r)/(1+r))));
}
static real B_h(const PowerSpectrum& P, real k, real q) {
    real r = q/k;
    real r2 = pow2(r), r4 = pow4(r);
    if(r < 1e-4)
        return P(q) * (-4 + 64./5.*r2 + 544./35.*r4);
    else if(fabs(1 - r) < 1e-8)
        return P(q) * (16 - 24*(1-r));
    else if(r > 10)
        return P(q) * (76./5. + 256./35/r2 - 32./7./r4);
    else
        return P(q) * (6/r2 + 1 + 9*r4 + 0.75/pow3(r)*pow2(1-r2)*(2+5*r2+3*r4)*log(pow2((1-r)/(1+r))));
}
static real B_i(const PowerSpectrum& P, real k, real q) {
    real r = q/k;
    real r2 = pow2(r), r4 = pow4(r);
    if(r < 1e-4)
        return P(q) * (12 - 16./5.*r2 + 400./7.*r4);
    else if(fabs(1 - r) < 1e-8)
        return P(q) * (44 - 60*(1-r));
    else if(r > 10)
        return P(q) * (268./5. - 16./35/r2 - 208./35./r4);
    else
        return P(q) * (6/r2 + 29 - 18*r2 + 27*r4 + 0.75/pow3(r)*pow2(1-r2)*(2+9*r2+9*r4)*log(pow2((1-r)/(1+r))));
}

/* k-dependent factors */
real ImprovedPT::G_f(real k) const {
    return k*k/(4*M_PI*M_PI) * 1/252. * Integrate<ExpSub>(bind(B_f, cref(P_i), k, _1), QMIN, QMAX, EPS);
}
real ImprovedPT::G_g(real k) const {
    return k*k/(4*M_PI*M_PI) * 1/84. * Integrate<ExpSub>(bind(B_g, cref(P_i), k, _1), QMIN, QMAX, EPS);
}
real ImprovedPT::G_h(real k) const {
    return k*k/(4*M_PI*M_PI) * 1/12. * Integrate<ExpSub>(bind(B_h, cref(P_i), k, _1), QMIN, QMAX, EPS);
}
real ImprovedPT::G_i(real k) const {
    return -k*k/(4*M_PI*M_PI) * 1/36. * Integrate<ExpSub>(bind(B_i, cref(P_i), k, _1), QMIN, QMAX, EPS);
}

void ImprovedPT::ComputeG_ab(real k, real eta1, real eta2, real& g11, real& g12, real& g21, real& g22) const {
    real eta = eta1 - eta2;
    real a = G_alpha(eta1, eta2);
    real bg = G_beta_g(eta1, eta2);
    real bd = G_beta_d(eta1, eta2);
    real gg = G_gamma_g(eta1, eta2);
    real gd = G_gamma_d(eta1, eta2);
    real d = G_delta(eta1, eta2);
    real f = G_f(k);
    real g = G_g(k);
    real h = G_h(k);
    real i = G_i(k);

    g11 = 0.6*exp(eta)*U(a*f - bg*i) + 0.4*exp(-1.5*eta)*U(d*g - gd*h);
    g12 = 0.4*exp(eta)*U(a*f - bg*h) - 0.4*exp(-1.5*eta)*U(d*f - gd*h);
    g21 = 0.6*exp(eta)*U(a*g + gg*h) - 0.6*exp(-1.5*eta)*U(d*g + bd*i);
    g22 = 0.4*exp(eta)*U(a*g - 3/2.*gg*i) + 0.6*exp(-1.5*eta)*U(d*f - 2/3.*bd*h);
}

void ImprovedPT::PrecomputePropagator() {
    g11.resize(Neta*Neta);
    g12.resize(Neta*Neta);
    g21.resize(Neta*Neta);
    g22.resize(Neta*Neta);

    /* Compute scale-dependent functions */
    array k(Nk), f(Nk), g(Nk), h(Nk), i(Nk);
    int n;
    #pragma omp parallel default(shared) private(n)
    #pragma omp for schedule(dynamic)
    for(n = 0; n < Nk; n++) {
        k[n] = n*kcut/(Nk-1);
        f[n] = G_f(k[n]);
        g[n] = G_g(k[n]);
        h[n] = G_h(k[n]);
        i[n] = G_i(k[n]);
    }

    /* Spline the propagator at each point in time */
    array g11tmp(Nk), g12tmp(Nk), g21tmp(Nk), g22tmp(Nk);
    for(int t = 0; t < Neta; t++) {
        real eta = t*deta;
        for(int tp = 0; tp < Neta; tp++) {
            real etap = tp*deta;
            real a = G_alpha(eta, etap);
            real bg = G_beta_g(eta, etap);
            real bd = G_beta_d(eta, etap);
            real gg = G_gamma_g(eta, etap);
            real gd = G_gamma_d(eta, etap);
            real d = G_delta(eta, etap);

            for(int n = 0; n < Nk; n++) {
                g11tmp[n] = 0.6*exp(eta-etap)*U(a*f[n] - bg*i[n]) + 0.4*exp(-1.5*(eta-etap))*U(d*g[n] - gd*h[n]);
                g12tmp[n] = 0.4*exp(eta-etap)*U(a*f[n] - bg*h[n]) - 0.4*exp(-1.5*(eta-etap))*U(d*f[n] - gd*h[n]);
                g21tmp[n] = 0.6*exp(eta-etap)*U(a*g[n] + gg*i[n]) - 0.6*exp(-1.5*(eta-etap))*U(d*g[n] + bd*i[n]);
                g22tmp[n] = 0.4*exp(eta-etap)*U(a*g[n] - 3/2.*gg*i[n]) + 0.6*exp(-1.5*(eta-etap))*U(d*f[n] - 2/3.*bd*h[n]);
            }
            g11[Neta*t + tp] = CubicSpline(Nk, k, g11tmp);
            g12[Neta*t + tp] = CubicSpline(Nk, k, g12tmp);
            g21[Neta*t + tp] = CubicSpline(Nk, k, g21tmp);
            g22[Neta*t + tp] = CubicSpline(Nk, k, g22tmp);
        }
    }
}

real ImprovedPT::G_11(real k, int t, int tp) const {
    if(t == -1)
        t += Neta;
    return g11[Neta*t + tp](k);
}

real ImprovedPT::G_12(real k, int t, int tp) const {
    if(t == -1)
        t += Neta;
    return g12[Neta*t + tp](k);
}

real ImprovedPT::G_21(real k, int t, int tp) const {
    if(t == -1)
        t += Neta;
    return g21[Neta*t + tp](k);
}

real ImprovedPT::G_22(real k, int t, int tp) const {
    if(t == -1)
        t += Neta;
    return g22[Neta*t + tp](k);
}

real ImprovedPT::G_1(real k, int t) const {
    return G_11(k, t) + G_12(k, t);
}

real ImprovedPT::G_2(real k, int t) const {
    return G_21(k, t) + G_22(k, t);
}


/************
 * Vertices *
 ************/

/* $\alpha(\vec{q}, \vec{k}-\vec{q}) = \frac{\vec{k}\cdot\vec{q}}{q^2}$ */
real ImprovedPT::alpha(real k, real q, real r) {
    real k2 = k*k, q2 = q*q, r2 = r*r;
    return (k2 + q2 - r2)/(2*q2);
}
/* $\beta(\vec{q}, \vec{k}-\vec{q}) = \frac{k^2 \vec{q}\cdot(\vec{k}-\vec{q})}{2 q^2 |\vec{k}-\vec{q}|^2}$ */
real ImprovedPT::beta(real k, real q, real r) {
    real k2 = k*k, q2 = q*q, r2 = r*r;
    return k2*(k2 - q2 - r2)/(4*q2*r2);
}


/***************************
 * Mode-coupling integrals *
 ***************************/

/* $I_a(\vec{k},\vec{q};\eta) = \int_0^\eta d\eta_1 ~ G_{ab}(k;\eta,\eta_1) \gamma_{bcd}(\vec{q},\vec{k}-\vec{q}) \tilde G_c(q;\eta_1) \tilde G_d(|\vec{k}-\vec{q}|;\eta_1)$ */
real ImprovedPT::I_1(real k, real q, real kq, int t) const {
    if(t == -1)
        t += Neta;
    real gamma112 = alpha(k, kq, q)/2;
    real gamma121 = alpha(k, q, kq)/2;
    real gamma222 = beta(k, q, kq);
    real i1[t+1];
    for(int t1 = 0; t1 <= t; t1++)
        i1[t1] = G_11(k,t,t1)*(gamma112*G_1(q,t1)*G_2(kq,t1) + gamma121*G_2(q,t1)*G_1(kq,t1))
               + G_12(k,t,t1)*gamma222*G_2(q,t1)*G_2(kq,t1);
    return DiscreteIntegrate(t+1, i1, deta);
}

real ImprovedPT::I_2(real k, real q, real kq, int t) const {
    if(t == -1)
        t += Neta;
    real gamma112 = alpha(k, kq, q)/2;
    real gamma121 = alpha(k, q, kq)/2;
    real gamma222 = beta(k, q, kq);
    real i2[t+1];
    for(int t1 = 0; t1 <= t; t1++)
        i2[t1] = G_21(k,t,t1)*(gamma112*G_1(q,t1)*G_2(kq,t1) + gamma121*G_2(q,t1)*G_1(kq,t1))
               + G_22(k,t,t1)*gamma222*G_2(q,t1)*G_2(kq,t1);
    return DiscreteIntegrate(t+1, i2, deta);
}

/* $J_a(\vec{k},\vec{p},\vec{q};\eta) = \int_0^\eta d\eta_1 ~ G_{ab}(k;\eta,\eta_1) \gamma_{bcd}(\vec{p},\vec{k}-\vec{p}) I_c(\vec{p},\vec{q};\eta_1) \tilde G_d(|\vec{k}-\vec{p}|;\eta_1)$ */
real ImprovedPT::J_1(real k, real p, real kp, real q, real pq, int t) const {
    if(t == -1)
        t += Neta;
    real gamma112 = alpha(k, kp, p)/2;
    real gamma121 = alpha(k, p, kp)/2;
    real gamma222 = beta(k, p, kp);
    real j1[t+1];
    for(int t1 = 0; t1 <= t; t1++) {
        real i1 = I_1(p, q, pq, t1);
        real i2 = I_2(p, q, pq, t1);
        j1[t1] = G_11(k,t,t1)*(gamma112*i1*G_2(kp,t1) + gamma121*i2*G_1(kp,t1))
               + G_12(k,t,t1)*gamma222*i2*G_2(kp,t1);
    }
    return DiscreteIntegrate(t+1, j1, deta);
}

real ImprovedPT::J_2(real k, real p, real kp, real q, real pq, int t) const {
    if(t == -1)
        t += Neta;
    real gamma112 = alpha(k, kp, p)/2;
    real gamma121 = alpha(k, p, kp)/2;
    real gamma222 = beta(k, p, kp);
    real j2[t+1];
    for(int t1 = 0; t1 <= t; t1++) {
        real i1 = I_1(p, q, pq, t1);
        real i2 = I_2(p, q, pq, t1);
        j2[t1] = G_21(k,t,t1)*(gamma112*i1*G_2(kp,t1) + gamma121*i2*G_1(kp,t1))
               + G_22(k,t,t1)*gamma222*i2*G_2(kp,t1);
    }
    return DiscreteIntegrate(t+1, j2, deta);
}

/* Integrand for 2-dimensional integral in $A(k,\mu;z)$, from Eq. (A3) of
 * Taruya, Nishimichi, and Saito:
 *   f_A2(k,f,\mu,q,x) = \sum_{m,n=1}^3 \mu^{2m} f^n \frac{k^3}{(2\pi)^2}
 *                       \times \{ A_{mn}(r,x) P_L(k;z) + \tilde{A}_{mn}(r,x) P_L(kr;z) \}
 *                       \times \frac{P_L(k\sqrt{1 + r^2 - 2rx};z)}{(1 + r^2 - 2rx)^2}
 * with q = kr. */
static real f_A2(const PowerSpectrum& P_L, real k, real f, real mu, real logq, real x) {
    real q = exp(logq);
    real r = q/k, r2 = r*r, r3 = r2*r, r4 = r2*r2;
    real x2 = x*x, x3 = x2*x, x4 = x2*x2;
    real rx = r*x, rx2 = r*x2, r2x = r2*x;
    real d = 1 + r2 - 2*rx;

    /* The integrand is smooth over r = x = 1, but the way we compute it leads
     * to problems.  We cheat to avoid this trouble spot. */
    if(fabs(d) < 1e-5)
        return f_A(P_L, k, f, mu, q, x - 0.0001);

    /* $A_{mn}(r,x)$ */
    real A12 = (r4/14.)*(x2 - 1)*(-1 + 7*rx - 6*x2);
    real A22 = (r3/14.)*(r2x*(13 - 41*x2) - 4*(x + 6*x3) + r*(5 + 9*x2 + 42*x4));
    real A23 = A12;
    real A33 = (r3/14.)*(1 - 7*rx + 6*x2)*(-2*x - r + 3*rx2);

    /* $\tilde{A}_{mn}(r,x)$ */
    real At11 = (1/7.)*(x + r - 2*r*x2)*(3*r + 7*x - 10*rx2);
    real At12 = (r/14.)*(x2 - 1)*(3*r + 7*x - 10*rx2);
    real At22 = (1/14.)*(28*x2 + rx*(25 - 81*x2) + r2*(1 - 27*x2 + 54*x4));
    real At23 = (r/14.)*(1 - x2)*(r - 7*x + 6*rx2);
    real At33 = (1/14.)*(r - 7*x + 6*rx2)*(-2*x - r + 3*rx2);

    real mu2 = mu*mu, mu4 = mu2*mu2, mu6 = mu4*mu2;
    real f2 = f*f, f3 = f2*f;
    real pk = P_L(k), pq = P_L(q), pkq = P_L(k*sqrt(d));
    real S = mu2*f*At11*pq + mu2*f2*(A12*pk + At12*pq) + mu4*f2*(A22*pk + At22*pq)
             + mu4*f3*(A23*pk + At23*pq) + mu6*f3*(A33*pk + At33*pq);
    S *= pkq / pow2(d);
    S *= q * pow2(k)/pow2(2*M_PI);      // multiply by an extra factor q/k compared to Eq. (A3), because we integrate in d(log q) = (k/q) dr
}

static real f_A1(const PowerSpectrum& P_L, real k, real f, real mu, real logq) {
    real q = exp(logq);
    real r = q/k, r2 = r*r, r3 = r2*r, r4 = r2*r2, r6 = r4*r2;
    real mu2 = mu*mu, mu4 = mu2*mu2, mu6 = mu4*mu2;
    real f2 = f*f, f3 = f2*f;
    real pk = P_L(k), pq = P_L(q);

    real a11, a12, a22, a33;
    if(r < 1e-2) {
        a11 = -(2/3.) + (8/7.)*r2 - (24/35.)*r4 + (24/245.)*r6;
        a12 = -(16/35.)*r2 + (48/248.)*r4 - (16/735.)*r6;
        a22 = -(4/3.) + (64/35.)*r2 - (288/245.)*r4 + (128/735.)*r6;
        a33 = -(2/3.) + (24/35.)*r2 - (24/49.)*r4 + (8/105.)*r6;
    }
    else if(fabs(r - 1) < 1e-4) {
        real s = r-1, s2 = s*s;
        a11 = -(2/21.) + (2/7.)*s - (5/7.)*s2;
        a12 = -(2/7.) - (2/7.)*s + (2/7.)*s2;
        a22 = -(10/21.) + (2/7.)*s - (8/7.)*s2;
        a33 = -(8/21.) - (3/7.)*s2;
    }
    else if(r > 100.) {
        a11 = (2/105.) - (24/245.)/r2 - (8/735.)/r4 - (8/2695.)/r6;
        a12 = -(16/35.) + (48/245.)/r2 - (16/735.)/r4 - (16/8085.)/r6;
        a22 = -(44/105.) - (32/735.)/r4 - (64/8085.)/r6;
        a33 = -(46/105.) + (24/245.)/r2 - (8/245.)/r4 - (8/1617.)/r6;
    }
    else {
        real t = r2 - 1, t2 = t*t, t3 = t2*t, t4 = t2*t2;
        real u = log(fabs((r+1)/(r-1)));
        a11 = -1/(84*r) * (2*r(19 - 24*r2 + 9*r4) - 9*t3*u);
        a12 = 1/(112*r3) * (2*r*(r2 + 1)*(3 - 14*r2 + 3*r4) - 3*t4*u);
        a22 = 1/(336*r3) * (2*r*(9 - 185*r2 + 159*r4 - 63*r6) + 9*t3*(7*r2 + 1)*u);
        a33 = 1/(336*r3) * (2*r*(9 - 109*r2 + 63*r4 - 27*r6) + 9*t3*(3*r2 + 1)*u);
    }

    return q*pow2(k/(2*M_PI)) * pk * (a11*mu2*f + a12*mu2*f2 + a22*mu4*f2 + a33*mu6*f3) * pq;
}

static real A(const PowerSpectrum& P_L, real k, real f, real mu) {
    real a[2] = { logqmin, -1 };
    real b[2] = { logqmax, +1 };
    real epsabs = 1e-3 * P_L(k);
    real A1 = Integrate(bind(f_A1, cref(P_L), k, f, mu, _1), logqmin, logqmax, epsrel, epsabs);
    real A2 = Integrate<2>(bind(f_A2, cref(P_L), k, f, mu, _1, _2), a, b, epsrel, epsabs);
    return A1 + A2;
}


/*****************************
 * Tree-level power spectrum *
 *****************************/

/* $P^{(1)}_{ab}(k;\eta) = \tilde G_a(k;\eta) \tilde G_b(k;\eta) P_i(k)$ */
real ImprovedPT::P1_11(real k) const {
    return (k == 0) ? 0 : pow2(G_1(k)) * P_i(k);
}
real ImprovedPT::P1_12(real k) const {
    return (k == 0) ? 0 : G_1(k)*G_2(k) * P_i(k);
}
real ImprovedPT::P1_22(real k) const {
    return (k == 0) ? 0 : pow2(G_2(k)) * P_i(k);
}

real ImprovedPT::P1(real k, int a, int b) const {
    switch(a*b) {
        case 1:
            return P1_11(k);
        case 2:
            return P1_12(k);
        case 4:
            return P1_22(k);
        default:
            warning("ImprovedPT: invalid indices, a = %d, b = %d\n", a, b);
            return 0;
    }
}

/*****************************
 * 1-loop mode-coupling term *
 *****************************/

/* $P^{(2)}_{ab}(k;\eta) = 2 \int \frac{d^3q}{(2\pi)^3} I_a(\vec{k},\vec{q};\eta,\eta_i) I_b(\vec{k},\vec{q};\eta,\eta_i) P_i(q) P_i(|\vec{k}-\vec{q}|)$ */
real ImprovedPT::P2_11(real k) const {
#if 0
    real a[2] = { QMIN, -1 };
    real b[2] = { kcut, +1 };
    return (k == 0) ? 0 : 2*Integrate<2>(bind(&ImprovedPT::f2_11, cref(*this), k, _1, _2), a, b, 1e-3, 1e-3*P_L(k));
#endif
    real a[2] = { k/M_SQRT2, 0 };
    real b[2] = { kcut, k/M_SQRT2 };
    return (k == 0) ? 0 : 2*Integrate<2>(bind(&ImprovedPT::F2_11, cref(*this), k, _1, _2), a, b, 1e-3, 1e-3*P_L(k)/2);
}
real ImprovedPT::P2_12(real k) const {
#if 0
    real a[2] = { QMIN, -1 };
    real b[2] = { kcut, +1 };
    return (k == 0) ? 0 : 2*Integrate<2>(bind(&ImprovedPT::f2_12, cref(*this), k, _1, _2), a, b, 1e-4, 1e-4*P_L(k));
#endif
    real a[2] = { k/M_SQRT2, 0 };
    real b[2] = { kcut, k/M_SQRT2 };
    return (k == 0) ? 0 : 2*Integrate<2>(bind(&ImprovedPT::F2_12, cref(*this), k, _1, _2), a, b, 1e-3, 1e-3*P_L(k)/2);
}
real ImprovedPT::P2_22(real k) const {
#if 0
    real a[2] = { QMIN, -1 };
    real b[2] = { kcut, +1 };
    return (k == 0) ? 0 : 2*Integrate<2>(bind(&ImprovedPT::f2_22, cref(*this), k, _1, _2), a, b, 1e-4, 1e-4*P_L(k));
#endif
    real a[2] = { k/M_SQRT2, 0 };
    real b[2] = { kcut, k/M_SQRT2 };
    return (k == 0) ? 0 : 2*Integrate<2>(bind(&ImprovedPT::F2_22, cref(*this), k, _1, _2), a, b, 1e-3, 1e-3*P_L(k)/2);
}

real ImprovedPT::P2(real k, int a, int b) const {
    switch(a*b) {
        case 1:
            return P2_11(k);
        case 2:
            return P2_12(k);
        case 4:
            return P2_22(k);
        default:
            warning("ImprovedPT: invalid indices, a = %d, b = %d\n", a, b);
            return 0;
    }
}

#if 0
real ImprovedPT::f2_11(real k, real q, real x) const {
    real kq = sqrt(k*k + q*q - 2*k*q*x);
    return q*q/pow2(2*M_PI) * pow2(I_1(k, q, kq)) * P_i(q)*P_i(kq);
}
real ImprovedPT::f2_12(real k, real q, real x) const {
    real kq = sqrt(k*k + q*q - 2*k*q*x);
    return q*q/pow2(2*M_PI) * I_1(k, q, kq)*I_2(k, q, kq) * P_i(q)*P_i(kq);
}
real ImprovedPT::f2_22(real k, real q, real x) const {
    real kq = sqrt(k*k + q*q - 2*k*q*x);
    return q*q/pow2(2*M_PI) * pow2(I_2(k, q, kq)) * P_i(q)*P_i(kq);
}
#endif
real ImprovedPT::F2_11(real k, real x, real y) const {
    real q = (x+y)/M_SQRT2;
    real r = (x-y)/M_SQRT2;
    return 1/pow2(2*M_PI) * q*r/k * 2*pow2(I_1(k,q,r)) * P_i(q)*P_i(r);
}
real ImprovedPT::F2_12(real k, real x, real y) const {
    real q = (x+y)/M_SQRT2;
    real r = (x-y)/M_SQRT2;
    return 1/pow2(2*M_PI) * q*r/k * 2*I_1(k,q,r)*I_2(k,q,r) * P_i(q)*P_i(r);
}
real ImprovedPT::F2_22(real k, real x, real y) const {
    real q = (x+y)/M_SQRT2;
    real r = (x-y)/M_SQRT2;
    return 1/pow2(2*M_PI) * q*r/k * 2*pow2(I_2(k,q,r)) * P_i(q)*P_i(r);
}

real ImprovedPT::P_L(real k) const {
    return pow2(D(z)/D(z_i)) * P_i(k);
}
