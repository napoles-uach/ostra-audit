#include <cmath>
#include <cstdio>
#include <vector>
#include "libecpint/bessel.hpp"
#include "libecpint/mathutil.hpp"

int main() {
    libecpint::initFactorials();
    const int lMax = 24;
    libecpint::BesselFunction bf(lMax, 1600, 200, 1e-15);

    int total_violations = 0;
    for (int L : {0, 1, 2, 3, 5, 6, 8, 10, 12, 16, 24}) {
        int viol_L = 0;
        double worst_ratio = 1e9, worst_z = 0.0, worst_ub = 0.0, worst_val = 0.0;
        for (double z = 0.025; z <= 35.0; z += 0.05) {
            double val = bf.calculate(z, L);
            double ub  = bf.upper_bound(z, L);
            if (ub < val * (1.0 - 1e-12)) {
                viol_L++;
                total_violations++;
                double r = ub / val;
                if (r < worst_ratio) {
                    worst_ratio = r; worst_z = z; worst_ub = ub; worst_val = val;
                }
            }
        }
        if (viol_L > 0) {
            std::printf("L=%2d: VIOLATIONS=%3d  worst at z=%5.2f: upper_bound=%.6e < K_L(z)=%.6e (ratio=%.6e)\n",
                        L, viol_L, worst_z, worst_ub, worst_val, worst_ratio);
        } else {
            std::printf("L=%2d: VIOLATIONS=  0  (valid upper bound across [0.02, 35.0])\n", L);
        }
    }
    std::printf("TOTAL VIOLATIONS: %d\n", total_violations);
    return 0;
}
