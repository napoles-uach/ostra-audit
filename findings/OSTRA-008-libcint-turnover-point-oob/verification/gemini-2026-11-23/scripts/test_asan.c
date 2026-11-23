#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint.h"
#include "cint_funcs.h"
#include "rys_roots.h"

static void run_eri_2e(int l0, int l1, int l2, int l3, double omega)
{
    int natm = 2;
    int nbas = 4;
    int atm[2 * ATM_SLOTS];
    int bas[4 * BAS_SLOTS];
    double env[100];
    memset(atm, 0, sizeof(atm));
    memset(bas, 0, sizeof(bas));
    memset(env, 0, sizeof(env));

    int off = PTR_ENV_START;
    env[PTR_RANGE_OMEGA] = omega;

    // Atom 0 at (0, 0, 0), Atom 1 at (0, 0, 1.5)
    atm[0 * ATM_SLOTS + CHARGE_OF] = 1;
    atm[0 * ATM_SLOTS + PTR_COORD] = off;
    env[off + 0] = 0.0; env[off + 1] = 0.0; env[off + 2] = 0.0;
    off += 3;

    atm[1 * ATM_SLOTS + CHARGE_OF] = 1;
    atm[1 * ATM_SLOTS + PTR_COORD] = off;
    env[off + 0] = 0.0; env[off + 1] = 0.0; env[off + 2] = 1.5;
    off += 3;

    int ls[4] = {l0, l1, l2, l3};
    int atms[4] = {0, 0, 1, 1};
    double exps[4] = {1.2, 0.9, 1.1, 0.8};

    for (int i = 0; i < 4; i++) {
        bas[i * BAS_SLOTS + ATOM_OF] = atms[i];
        bas[i * BAS_SLOTS + ANG_OF] = ls[i];
        bas[i * BAS_SLOTS + NPRIM_OF] = 1;
        bas[i * BAS_SLOTS + NCTR_OF] = 1;
        bas[i * BAS_SLOTS + PTR_EXP] = off;
        env[off++] = exps[i];
        bas[i * BAS_SLOTS + PTR_COEFF] = off;
        env[off++] = CINTgto_norm(ls[i], exps[i]);
    }

    int shls[4] = {0, 1, 2, 3};
    int di = 2 * l0 + 1;
    int dj = 2 * l1 + 1;
    int dk = 2 * l2 + 1;
    int dl = 2 * l3 + 1;
    size_t n = (size_t)di * dj * dk * dl;
    double *buf = (double *)calloc(n, sizeof(double));
    CACHE_SIZE_T cache_size = int2e_sph(NULL, NULL, shls, atm, natm, bas, nbas, env, NULL, NULL);
    double *cache = (double *)malloc((size_t)cache_size * sizeof(double));
    int ret = int2e_sph(buf, NULL, shls, atm, natm, bas, nbas, env, NULL, cache);
    double s = 0.0;
    for (size_t i = 0; i < n; i++) s += fabs(buf[i]);
    printf("int2e_sph l=(%d,%d|%d,%d) omega=%g -> ret=%d sum_abs=%.12e buf[0]=%.12e\n",
           l0, l1, l2, l3, omega, ret, s, buf[0]);
    free(buf);
    free(cache);
}

static void run_eri_3c2e(int l0, int l1, int l2, double omega)
{
    int natm = 2;
    int nbas = 3;
    int atm[2 * ATM_SLOTS];
    int bas[3 * BAS_SLOTS];
    double env[100];
    memset(atm, 0, sizeof(atm));
    memset(bas, 0, sizeof(bas));
    memset(env, 0, sizeof(env));

    int off = PTR_ENV_START;
    env[PTR_RANGE_OMEGA] = omega;

    atm[0 * ATM_SLOTS + CHARGE_OF] = 1;
    atm[0 * ATM_SLOTS + PTR_COORD] = off;
    env[off + 0] = 0.0; env[off + 1] = 0.0; env[off + 2] = 0.0;
    off += 3;

    atm[1 * ATM_SLOTS + CHARGE_OF] = 1;
    atm[1 * ATM_SLOTS + PTR_COORD] = off;
    env[off + 0] = 0.0; env[off + 1] = 0.0; env[off + 2] = 1.5;
    off += 3;

    int ls[3] = {l0, l1, l2};
    int atms[3] = {0, 0, 1};
    double exps[3] = {1.2, 0.9, 1.1};

    for (int i = 0; i < 3; i++) {
        bas[i * BAS_SLOTS + ATOM_OF] = atms[i];
        bas[i * BAS_SLOTS + ANG_OF] = ls[i];
        bas[i * BAS_SLOTS + NPRIM_OF] = 1;
        bas[i * BAS_SLOTS + NCTR_OF] = 1;
        bas[i * BAS_SLOTS + PTR_EXP] = off;
        env[off++] = exps[i];
        bas[i * BAS_SLOTS + PTR_COEFF] = off;
        env[off++] = CINTgto_norm(ls[i], exps[i]);
    }

    int shls[3] = {0, 1, 2};
    size_t n = (size_t)(2*l0+1) * (2*l1+1) * (2*l2+1);
    double *buf = (double *)calloc(n, sizeof(double));
    CACHE_SIZE_T cache_size = int3c2e_sph(NULL, NULL, shls, atm, natm, bas, nbas, env, NULL, NULL);
    double *cache = (double *)malloc((size_t)cache_size * sizeof(double));
    int ret = int3c2e_sph(buf, NULL, shls, atm, natm, bas, nbas, env, NULL, cache);
    double s = 0.0;
    for (size_t i = 0; i < n; i++) s += fabs(buf[i]);
    printf("int3c2e_sph l=(%d,%d|%d) omega=%g -> ret=%d sum_abs=%.12e\n",
           l0, l1, l2, omega, ret, s);
    free(buf);
    free(cache);
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <case>\n", argv[0]);
        return 2;
    }
    const char *mode = argv[1];
    double f[128];
    long double lf[128];

    if (strcmp(mode, "boys_ok") == 0) {
        gamma_inc_like(f, 10.0, 39);
        fmt_erfc_like(f, 10.0, 0.3, 39);
        printf("boys_ok m=39: f[0]=%.12e f[39]=%.12e\n", f[0], f[39]);
    } else if (strcmp(mode, "boys_gamma_40") == 0) {
        gamma_inc_like(f, 10.0, 40);
        printf("boys_gamma_40: f[40]=%.12e\n", f[40]);
    } else if (strcmp(mode, "boys_lgamma_40") == 0) {
        lgamma_inc_like(lf, 10.0L, 40);
        printf("boys_lgamma_40: f[40]=%.12Le\n", lf[40]);
    } else if (strcmp(mode, "boys_erfc_40") == 0) {
        fmt_erfc_like(f, 10.0, 0.3, 40);
        printf("boys_erfc_40: f[40]=%.12e\n", f[40]);
    } else if (strcmp(mode, "boys_lerfc_41") == 0) {
        fmt_lerfc_like(lf, 10.0L, 0.3L, 41);
        printf("boys_lerfc_41: f[41]=%.12Le\n", lf[41]);
    } else if (strcmp(mode, "eri_sr_l9") == 0) {
        run_eri_2e(10, 10, 10, 9, -0.5);
    } else if (strcmp(mode, "eri_sr_l10") == 0) {
        run_eri_2e(10, 10, 10, 10, -0.5);
    } else if (strcmp(mode, "eri_full_l10") == 0) {
        run_eri_2e(10, 10, 10, 10, 0.0);
    } else if (strcmp(mode, "eri_lr_l10") == 0) {
        run_eri_2e(10, 10, 10, 10, 0.5);
    } else if (strcmp(mode, "eri_sr_3c2e") == 0) {
        run_eri_3c2e(14, 14, 12, -0.5);
    } else {
        fprintf(stderr, "Unknown mode: %s\n", mode);
        return 2;
    }
    return 0;
}
