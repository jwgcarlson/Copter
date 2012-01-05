#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>

#include <Copter/Cosmology.h>
#include <Copter/LinearPS.h>
#include <Copter/SPT.h>


int main(int argc, char* argv[]) {
    /* Set parameters */
    const char* cstr = "example";
    real z = 0;
    real epsrel = 1e-4;
    int Nk = 1000;
    real kmin = 1e-4;
    real kmax = 10;
    const char* output = "spt2.dat";

    /* Open output file */
    FILE* fp = fopen(output, "w");
    fprintf(fp, "# 1-loop SPT for %s cosmology at z = %g\n", cstr, z);
    fprintf(fp, "# k       -- P_11      -- P_12      -- P_22\n");

    /* Calculate SPT power spectrum */
    Cosmology C(cstr);
    LinearPS P_L(C, z);
    SPT spt(C, P_L, epsrel);
    real p11, p12, p22;
    for(int i = 0; i < Nk; i++) {
        real k = kmin * exp(i*log(kmax/kmin)/(Nk-1));
        p11 = spt.P(k, 1,1);    // density-density power spectrum
        p12 = spt.P(k, 1,2);    // density-velocity cross spectrum
        p22 = spt.P(k, 2,2);    // velocity-velocity power spectrum
        printf("%d/%d: %e %e %e %e\n", i+1, Nk, k, p11, p12, p22);
        fprintf(fp, "%e %e %e %e\n", k, p11, p12, p22);
    }

    fclose(fp);
    return 0;
}
