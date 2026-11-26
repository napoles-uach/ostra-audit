#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "dftd3.h"

static void check_err(dftd3_error err, const char *ctx) {
    if (dftd3_check_error(err)) {
        char buf[512];
        int len = sizeof(buf);
        dftd3_get_error(err, buf, &len);
        fprintf(stderr, "ERROR in %s: %s\n", ctx, buf);
        exit(2);
    }
}

int main(void) {
    const int nat = 3;
    const int num[3] = {8, 1, 1};
    const double xyz[9] = {
        0.0, 0.0, 0.0,
        1.80, 0.35, 0.15,
       -0.45, 1.75, -0.25
    };
    const double lat[9] = {
        8.0, 0.2, 0.1,
        0.3, 8.2, 0.15,
        0.1, 0.25, 8.5
    };
    const bool pbc[3] = {true, true, true};

    dftd3_error err = dftd3_new_error();
    dftd3_structure mol = dftd3_new_structure(err, nat, num, xyz, lat, pbc);
    check_err(err, "new_structure");
    dftd3_model disp = dftd3_new_d3_model(err, mol);
    check_err(err, "new_d3_model");
    dftd3_param param = dftd3_load_rational_damping(err, "pbe", true);
    check_err(err, "load_rational_damping");

    double e_ref = 999.0, g_ref[9], s_ref[9];
    for (int i = 0; i < 9; i++) { g_ref[i] = 999.0; s_ref[i] = 999.0; }
    dftd3_get_dispersion(err, mol, disp, param, &e_ref, g_ref, s_ref);
    check_err(err, "ref");

    double e = 999.0;
    dftd3_get_dispersion(err, mol, disp, param, &e, NULL, NULL);
    check_err(err, "(NULL,NULL)");
    printf("C API dftd3_get_dispersion (NULL,NULL):  |dE|=%.3e\n", fabs(e - e_ref));

    double g[9];
    e = 999.0;
    for (int i = 0; i < 9; i++) g[i] = 999.0;
    dftd3_get_dispersion(err, mol, disp, param, &e, g, NULL);
    check_err(err, "(grad,NULL)");
    double mg = 0.0;
    for (int i = 0; i < 9; i++) mg = fmax(mg, fabs(g[i] - g_ref[i]));
    printf("C API dftd3_get_dispersion (grad,NULL):  |dE|=%.3e max|dg|=%.3e\n", fabs(e - e_ref), mg);

    double s[9];
    e = 999.0;
    for (int i = 0; i < 9; i++) s[i] = 999.0;
    dftd3_get_dispersion(err, mol, disp, param, &e, NULL, s);
    check_err(err, "(NULL,sigma)");
    double ms = 0.0;
    for (int i = 0; i < 9; i++) ms = fmax(ms, fabs(s[i] - s_ref[i]));
    printf("C API dftd3_get_dispersion (NULL,sigma): |dE|=%.3e max|ds|=%.3e\n", fabs(e - e_ref), ms);

    dftd3_delete_param(&param);
    dftd3_delete_model(&disp);
    dftd3_delete_structure(&mol);
    dftd3_delete_error(&err);
    return 0;
}
