#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint.h"
#include "gto/nr_ecp.h"

int ECPscalar_sph(double*, int*, int*, int*, int, int*, int, double*, ECPOpt*, double*);
int ECPscalar_ipnuc_sph(double*, int*, int*, int*, int, int*, int, double*, ECPOpt*, double*);
int ECPscalar_ipipnuc_sph(double*, int*, int*, int*, int, int*, int, double*, ECPOpt*, double*);
int ECPscalar_ipnucip_sph(double*, int*, int*, int*, int, int*, int, double*, ECPOpt*, double*);
int ECPscalar_iprinv_sph(double*, int*, int*, int*, int, int*, int, double*, ECPOpt*, double*);
int ECPscalar_ipiprinv_sph(double*, int*, int*, int*, int, int*, int, double*, ECPOpt*, double*);
int ECPscalar_iprinvip_sph(double*, int*, int*, int*, int, int*, int, double*, ECPOpt*, double*);
int ECPscalar_ignuc_sph(double*, int*, int*, int*, int, int*, int, double*, ECPOpt*, double*);
void ECPscalar_optimizer(ECPOpt**, int*, int, int*, int, double*);
void ECPdel_optimizer(ECPOpt**);

typedef int (*IntorFn)(double*, int*, int*, int*, int, int*, int, double*, ECPOpt*, double*);

static IntorFn get_fn(const char *name, int *comp, int *use_opt) {
    *use_opt = 1;
    if (!strcmp(name, "scalar"))   { *comp = 1; return ECPscalar_sph; }
    if (!strcmp(name, "ipnuc"))    { *comp = 3; return ECPscalar_ipnuc_sph; }
    if (!strcmp(name, "ipipnuc"))  { *comp = 9; return ECPscalar_ipipnuc_sph; }
    if (!strcmp(name, "ipnucip"))  { *comp = 9; return ECPscalar_ipnucip_sph; }
    if (!strcmp(name, "iprinv"))   { *comp = 3; *use_opt = 0; return ECPscalar_iprinv_sph; }
    if (!strcmp(name, "ipiprinv")) { *comp = 9; *use_opt = 0; return ECPscalar_ipiprinv_sph; }
    if (!strcmp(name, "iprinvip")) { *comp = 9; *use_opt = 0; return ECPscalar_iprinvip_sph; }
    if (!strcmp(name, "ignuc"))    { *comp = 3; return ECPscalar_ignuc_sph; }
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 3) return 2;
    const char *path = argv[1];
    const char *intor = argv[2];
    const char *dump_out = argc > 3 ? argv[3] : NULL;
    int comp = 1, use_opt = 1;
    IntorFn fn = get_fn(intor, &comp, &use_opt);
    if (!fn) return 2;

    FILE *f = fopen(path, "rb");
    if (!f) { perror("fopen"); return 2; }
    int hdr[4];
    fread(hdr, sizeof(int), 4, f);
    int natm = hdr[0], nbas_tot = hdr[1], nbas_ao = hdr[2], nenv = hdr[3];
    int *atm = malloc(sizeof(int) * natm * ATM_SLOTS);
    int *bas = malloc(sizeof(int) * nbas_tot * BAS_SLOTS);
    double *env = malloc(sizeof(double) * nenv);
    fread(atm, sizeof(int), natm * ATM_SLOTS, f);
    fread(bas, sizeof(int), nbas_tot * BAS_SLOTS, f);
    fread(env, sizeof(double), nenv, f);
    fclose(f);

    ECPOpt *opt = NULL;
    if (use_opt) ECPscalar_optimizer(&opt, atm, natm, bas, nbas_ao, env);

    FILE *fo = dump_out ? fopen(dump_out, "wb") : NULL;
    double global_max = 0.0;
    for (int i = 0; i < nbas_ao; i++) {
        for (int j = 0; j < nbas_ao; j++) {
            int shls[2] = {i, j};
            int li = bas[i*BAS_SLOTS + ANG_OF], lj = bas[j*BAS_SLOTS + ANG_OF];
            int di = (2*li + 1) * bas[i*BAS_SLOTS + NCTR_OF];
            int dj = (2*lj + 1) * bas[j*BAS_SLOTS + NCTR_OF];
            int dims[2] = {di, dj};
            size_t nout = (size_t)comp * di * dj;
            double *out = calloc(nout, sizeof(double));
            int csz = fn(NULL, dims, shls, atm, natm, bas, nbas_ao, env, opt, NULL);
            double *cache = malloc(sizeof(double) * csz);
            fn(out, dims, shls, atm, natm, bas, nbas_ao, env, opt, cache);
            for (size_t k = 0; k < nout; k++) {
                if (fabs(out[k]) > global_max) global_max = fabs(out[k]);
            }
            if (fo) fwrite(out, sizeof(double), nout, fo);
            free(cache);
            free(out);
        }
    }
    if (fo) fclose(fo);
    printf("PASS %s max=%.12e\n", intor, global_max);
    ECPdel_optimizer(&opt);
    free(atm); free(bas); free(env);
    return 0;
}
