"""B_k(x) = int_{-1}^{1} t^k e^{-xt} dt: 12-term series (gcp.f90 bint) vs closed-form upward recursion, against mpmath."""
import mpmath as mp, math
mp.mp.dps = 30
exact = lambda x, k: mp.quad(lambda t: t**k * mp.exp(-x * t), [-1, 1])
def series(x, k):
    return sum((1 - (-1)**(k + i + 1)) / (math.factorial(i) * (k + i + 1)) * (-x)**i for i in range(13))
def closed(x, k):
    b = (math.exp(x) - math.exp(-x)) / x
    for j in range(1, k + 1): b = ((-1)**j * math.exp(x) - math.exp(-x)) / x + j / x * b
    return b
for x in (0.05, 0.1, 1.0, 3.0, 5.0):
    print("x=%-5g" % x, " ".join("k%d: series %.1e closed %.1e" % (k, abs(series(x, k) / exact(x, k) - 1), abs(closed(x, k) / exact(x, k) - 1)) for k in (0, 2, 5, 7, 8)))
