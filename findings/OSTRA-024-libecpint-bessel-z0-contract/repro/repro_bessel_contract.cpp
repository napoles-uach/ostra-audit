#include <cstdio>
#include <vector>
#include "libecpint/bessel.hpp"
#include "libecpint/mathutil.hpp"

int main() {
    libecpint::initFactorials();
    const int lMax = 6;
    libecpint::BesselFunction bf(lMax, 1600, 200, 1e-15);

    // 1) Vector overload: reuse buffer populated by z = 5.0, then evaluate at z = 0.0
    std::vector<double> vals(lMax + 1, 0.0);
    bf.calculate(5.0, lMax, vals);
    std::vector<double> prev = vals;
    bf.calculate(0.0, lMax, vals);
    std::printf("=== Vector calculate(0.0, %d, vals) on reused buffer ===\n", lMax);
    for (int l = 0; l <= lMax; l++) {
        double exact = (l == 0) ? 1.0 : 0.0;
        std::printf("  l=%d: got=%+.10e  prev=%+.10e  exact=%.1f  stale=%s\n",
                    l, vals[l], prev[l], exact, (l > 0 && vals[l] != 0.0) ? "YES" : "no");
    }

    // 2) Scalar overload: calculate(z <= 0, L) for L = 0..6
    std::printf("=== Scalar calculate(z, L) for z = 0.0 and z = -1.0 ===\n");
    for (int l = 0; l <= lMax; l++) {
        double v0 = bf.calculate(0.0, l);
        double vm = bf.calculate(-1.0, l);
        double vp = bf.calculate(1e-12, l);
        double exact = (l == 0) ? 1.0 : 0.0;
        std::printf("  l=%d: calc(0.0)=%+.6e  calc(-1.0)=%+.6e  calc(1e-12)=%+.6e  exact(0)=%.1f\n",
                    l, v0, vm, vp, exact);
    }
    return 0;
}
