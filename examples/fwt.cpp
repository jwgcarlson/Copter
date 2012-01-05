#include <cstdio>

#include <Copter/Cosmology.h>
#include <Copter/LinearPS.h>
#include <Copter/FlowingWithTime.h>


int main(int argc, char* argv[]) {
    /* Set parameters */
    const char* cstr = "example";       // cosmology
    const real z_i = 100;               // starting redshift
    const int Nk = 100;                 // number of points in interpolated power spectrum
    const real kmin = 1e-3;
    const real kmax = 10;
    const int Nz = 2;                   // number of redshift outputs
    const real z[Nz] = { 0, 1 };
    const char* output = "fwt.dat";     // output file

    /* Initialize */
    Cosmology C(cstr);
    LinearPS P_i(C, z_i);
    array k(Nk);
    for(int i = 0; i < Nk; i++)
        k[i] = kmin*exp(i*log(kmax/kmin)/(Nk-1));
    FlowingWithTime fwt(C, z_i, P_i, k);

    /* Calculate power spectrum */
    vector<InterpolatedP_ab> P = fwt.CalculateP_ab(Nz, z);

    /* Write data to file */
    FILE* fp = fopen(output, "w");
    fprintf(fp, "# FlowingWithTime for %s cosmology\n", cstr);
    fprintf(fp, "# Starting redshift: %g", z_i);
    fprintf(fp, "# Output redshifts:");
    fprintf(fp, "# Columns: k   P_11(k; z[0])   P_12(k; z[0])   P_22(k; z[0])   P_11(k; z[1]) ...\n");
    for(int j = 0; j < Nz; j++)
        fprintf(fp, " %g", z[j]);
    fprintf(fp, "\n");
    for(int i = 0; i < Nk; i++) {
        fprintf(fp, "%e", k[i]);
        for(int j = 0; j < Nz; j++)
            fprintf(fp, "  %e %e %e", P[j](1,1, k[i]), P[j](1,2, k[i]), P[j](2,2, k[i]));
        fprintf(fp, "\n");
    }
    fclose(fp);

    return 0;
}
