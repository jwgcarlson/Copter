#include <cstdlib>

#include <Copter/Common.h>
#include <Copter/Cosmology.h>
#include <Copter/ODE.h>
#include <Copter/Spline.h>
#include <Copter/array.h>

using Constants::c;
using Constants::Mpc;

Cosmology C;

real drdz(real z, real r) {
    real a = 1/(1 + z);
    return c/C.H(a) / (Mpc/C.h);
}

int main(int argc, char* argv[]) {
    if(argc < 3)
        error("Usage: %s cosmology z1 [z2 ...]\n", argv[0]);

    C = Cosmology(argv[1]);

    array z;
    real zmax = 0;
    for(int i = 2; i < argc; i++) {
        z.push_back(atof(argv[i]));
        if(z.back() > zmax)
            zmax = z.back();
    }

    vector<real> ztmp, rtmp;
    RKDP(0, zmax+1e-3, 0, &drdz, ztmp, rtmp, 1e-5);
    Spline r = CubicSpline(ztmp, rtmp);

    for(int i = 0; i < (int)z.size(); i++) {
        info("z = %g :\n", z[i]);
        info("  D_C = %g Mpc/h\n", r(z[i]));
        info("  D_A = %g Mpc/h\n", r(z[i])/(1+z[i]));
        info("  D_L = %g Mpc/h\n", r(z[i])*(1+z[i]));
    }

    return 0;
}
