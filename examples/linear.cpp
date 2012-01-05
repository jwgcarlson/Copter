#include <cstdio>

#include <Copter/Cosmology.h>
#include <Copter/LinearPS.h>


int main(int argc, char* argv[]) {
    /* Set parameters */
    const char* cstr = "example";
    real z = 0.5;
    int Nk = 5000;
    real kmin = 1e-5;
    real kmax = 1e2;
    const char* output = "linear.dat";

    /* Open output file */
    FILE* fp = fopen(output, "w");
    fprintf(fp, "# Linear theory for %s cosmology at z = %g\n", cstr, z);
    fprintf(fp, "# k       -- P_L\n");

    /* Calculate SPT power spectrum */
    Cosmology C(cstr);
    LinearPS P_L(C, z);
    array k = array::logspace(kmin, kmax, Nk);
    array pk = P_L(k);

    for(int i = 0; i < Nk; i++)
        fprintf(fp, "%e %e\n", k[i], pk[i]);

    fclose(fp);
    return 0;
}
