#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "dftd4.h"

static void check_err(dftd4_error err, const char *ctx) {
    if (dftd4_check_error(err)) {
        char buf[512];
        int len = sizeof(buf);
        dftd4_get_error(err, buf, &len);
        fprintf(stderr, "ERROR in %s: %s\n", ctx, buf);
        exit(2);
    }
}

int main(int argc, char **argv) {
    const char *mode = argc > 1 ? argv[1] : "all_safe";
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
    const double charge = 0.0;

    dftd4_error err = dftd4_new_error();
    dftd4_structure mol = dftd4_new_structure(err, nat, num, xyz, &charge, lat, pbc);
    check_err(err, "new_structure");
    dftd4_model disp = dftd4_new_d4_model(err, mol);
    check_err(err, "new_d4_model");
    dftd4_param param = dftd4_load_rational_damping(err, "pbe", true);
    check_err(err, "load_rational_damping");

    double e_ref = 999.0, g_ref[9], s_ref[9];
    for (int i = 0; i < 9; i++) { g_ref[i] = 999.0; s_ref[i] = 999.0; }
    dftd4_get_dispersion(err, mol, disp, param, &e_ref, g_ref, s_ref);
    check_err(err, "get_dispersion(T,T)");
    printf("REF (grad,sigma): E=%.12e max|g|=%.6e max|s|=%.6e\n", e_ref,
           fmax(fabs(g_ref[0]), fabs(g_ref[3])), fmax(fabs(s_ref[0]), fabs(s_ref[4])));

    if (strcmp(mode, "all_safe") == 0) {
        // (F, F)
        double e = 999.0;
        dftd4_get_dispersion(err, mol, disp, param, &e, NULL, NULL);
        check_err(err, "get_dispersion(NULL,NULL)");
        printf("C API (NULL,NULL): |dE|=%.3e\n", fabs(e - e_ref));

        // (T, F): gradient != NULL, sigma == NULL
        double g[9];
        e = 999.0;
        for (int i = 0; i < 9; i++) g[i] = 999.0;
        dftd4_get_dispersion(err, mol, disp, param, &e, g, NULL);
        check_err(err, "get_dispersion(grad,NULL)");
        double mg = 0.0;
        for (int i = 0; i < 9; i++) mg = fmax(mg, fabs(g[i] - g_ref[i]));
        printf("C API (grad,NULL): |dE|=%.3e max|dg|=%.3e\n", fabs(e - e_ref), mg);

        // Properties: all 16 combinations!
        double cn_ref[3], q_ref[3], c6_ref[9], a_ref[3];
        dftd4_get_properties(err, mol, disp, cn_ref, q_ref, c6_ref, a_ref);
        check_err(err, "get_properties(1,1,1,1)");
        double max_prop_err = 0.0;
        for (int mask = 0; mask < 16; mask++) {
            double cn[3], q[3], c6[9], a[3];
            for (int i = 0; i < 3; i++) { cn[i]=999.0; q[i]=999.0; a[i]=999.0; }
            for (int i = 0; i < 9; i++) c6[i]=999.0;
            dftd4_get_properties(err, mol, disp,
                (mask & 1) ? cn : NULL,
                (mask & 2) ? q : NULL,
                (mask & 4) ? c6 : NULL,
                (mask & 8) ? a : NULL);
            check_err(err, "get_properties(mask)");
            if (mask & 1) for (int i=0;i<3;i++) max_prop_err = fmax(max_prop_err, fabs(cn[i]-cn_ref[i]));
            if (mask & 2) for (int i=0;i<3;i++) max_prop_err = fmax(max_prop_err, fabs(q[i]-q_ref[i]));
            if (mask & 4) for (int i=0;i<9;i++) max_prop_err = fmax(max_prop_err, fabs(c6[i]-c6_ref[i]));
            if (mask & 8) for (int i=0;i<3;i++) max_prop_err = fmax(max_prop_err, fabs(a[i]-a_ref[i]));
        }
        printf("C API dftd4_get_properties (all 16 masks): max_err=%.3e\n", max_prop_err);
    } else if (strcmp(mode, "sigma_only") == 0) {
        double e = 999.0, s[9];
        for (int i = 0; i < 9; i++) s[i] = 999.0;
        printf("Calling dftd4_get_dispersion(NULL, sigma)...\n");
        fflush(stdout);
        dftd4_get_dispersion(err, mol, disp, param, &e, NULL, s);
        check_err(err, "get_dispersion(NULL,sigma)");
        double ms = 0.0;
        for (int i = 0; i < 9; i++) ms = fmax(ms, fabs(s[i] - s_ref[i]));
        printf("C API (NULL,sigma): |dE|=%.3e max|ds|=%.3e\n", fabs(e - e_ref), ms);
    }

    dftd4_delete_param(&param);
    dftd4_delete_model(&disp);
    dftd4_delete_structure(&mol);
    dftd4_delete_error(&err);
    return 0;
}
