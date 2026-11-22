/* Reproducer: LDA_X_1D_SOFT exchange energy collapses for beta*k_F > ~38.6.
   Expected: eps_x -> -1/(2 beta) = -0.5 (beta = 1) monotonically at high density.
   Build: cc repro_1d_soft.c -I<libxc>/include -L<libxc>/lib -lxc -lm */
#include <stdio.h>
#include <xc.h>
#include <xc_funcs.h>

int main(void) {
  double rho[] = {10.0, 24.5, 24.7, 30.0, 100.0, 1000.0};
  double zk[6], vrho[6];
  xc_func_type f;
  if (xc_func_init(&f, XC_LDA_X_1D_SOFT, XC_UNPOLARIZED) != 0) return 1;
  xc_lda_exc_vxc(&f, 6, rho, zk, vrho);
  printf("libxc %s\n", xc_version_string());
  for (int i = 0; i < 6; i++)
    printf("rho = %8.2f   eps_x = % .10f   v_x = % .10f\n", rho[i], zk[i], vrho[i]);
  xc_func_end(&f);
  return 0;
}
