#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>

#include <Copter/Cosmology.h>
#include <Copter/LinearPS.h>
#include <Copter/RPT.h>


int mkdirp(const char* path) {
    int status = mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(status == -1) {
        if(errno == EEXIST)
            status = errno = 0; // ignore this error
        else
            perror("mkdirp(): ");

    }
    return status;
}


int main(int argc, char* argv[]) {
    if(argc < 3)
        error("Usage: %s cosmology z [z_i Neta kcut Nk kmin kmax]\n", argv[0]);

    const char* cosmology = argv[1];
    real z = atof(argv[2]);
    real z_i = (argc >= 4) ? atof(argv[3]) : 100;
    int Neta = (argc >= 5) ? atoi(argv[4]) : 50;
    real kcut = (argc >= 6) ? atof(argv[5]) : 10;
    int Nk = (argc >= 7) ? atoi(argv[6]) : 500;
    real kmin = (argc >= 8) ? atof(argv[7]) : 1e-3;
    real kmax = (argc >= 9) ? atof(argv[8]) : 0.5;

    /* Open output file */
    char filename[256];
    mkdirp(cosmology);
    snprintf(filename, 256, "%s/z%g", cosmology, z);
    mkdirp(filename);
    snprintf(filename, 256, "%s/z%g/rpt2.dat", cosmology, z);
    FILE* f = fopen(filename, "w");
    fprintf(f, "# 1-loop RPT for %s cosmology at z = %g\n", cosmology, z);
    fprintf(f, "# k       -- P_11      -- P_12      -- P_22\n");

    Cosmology C(cosmology);
    LinearPS P_i(C, z_i);
    RPT rpt(C, P_i, z_i, z, Neta, 1000, kcut);
    for(int i = 0; i < Nk; i++) {
        info("%g%%\n", 100*(i+1.)/Nk);
        real k = kmin * exp(i*log(kmax/kmin)/(Nk-1));
        fprintf(f, "%e %e %e %e\n", k, rpt.P1(k, 1,1) + rpt.P2(k, 1,1),
                                       rpt.P1(k, 1,2) + rpt.P2(k, 1,2),
                                       rpt.P1(k, 2,2) + rpt.P2(k, 2,2));
    }

    fclose(f);
    return 0;
}
