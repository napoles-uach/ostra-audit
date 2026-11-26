#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tblite.h"

int main(void) {
    const int nat = 3;
    const int num[3] = {8, 1, 1};
    const double xyz[9] = {
        0.0, 0.0, 0.0,
        1.80, 0.35, 0.15,
       -0.45, 1.75, -0.25
    };
    tblite_error err = tblite_new_error();
    tblite_context ctx = tblite_new_context();
    tblite_set_context_verbosity(ctx, 0);
    tblite_structure mol = tblite_new_structure(err, nat, num, xyz, NULL, NULL, NULL, NULL);
    tblite_calculator calc = tblite_new_gfn2_calculator(ctx, mol, NULL);
    tblite_result res = tblite_new_result();

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check_context(ctx)) {
        char buf[512];
        int bufsz = (int)sizeof(buf);
        tblite_get_context_error(ctx, buf, &bufsz);
        fprintf(stderr, "C API ctx error: %s\n", buf);
        return 2;
    }

    double e = 999.0, g[9], s[9], q[3];
    for (int i=0;i<9;i++) { g[i]=999.0; s[i]=999.0; }
    for (int i=0;i<3;i++) q[i]=999.0;
    tblite_get_result_energy(err, res, &e);
    tblite_get_result_gradient(err, res, g);
    tblite_get_result_virial(err, res, s);
    tblite_get_result_charges(err, res, q);
    printf("tblite C API: E=%.12e max|g|=%.6e max|s|=%.6e sum(q)=%.3e err=%d\n",
           e, fmax(fabs(g[0]), fabs(g[3])), fmax(fabs(s[0]), fabs(s[4])),
           q[0]+q[1]+q[2], tblite_check_error(err));

    tblite_delete_result(&res);
    tblite_delete_calculator(&calc);
    tblite_delete_structure(&mol);
    tblite_delete_context(&ctx);
    tblite_delete_error(&err);
    return 0;
}
