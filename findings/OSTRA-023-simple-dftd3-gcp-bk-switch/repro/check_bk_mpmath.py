import math, mpmath as mp
mp.mp.dps = 50

def B_exact(k, x):
    x = mp.mpf(x)
    s = mp.mpf(0)
    for i in range(0, 120):
        if (k + i) % 2 == 0:
            s += 2 * ((-x)**i) / (mp.factorial(i) * (k + i + 1))
    return float(s)

def bint_12(k, x):
    if abs(x) < 1e-6:
        return (1.0 + (-1.0)**k) / (k + 1.0)
    s = 0.0
    for i in range(0, 13):
        xx = 1.0 - ((-1.0)**(k + i + 1))
        yy = math.factorial(i) * (k + i + 1)
        s += xx / yy * ((-x)**i)
    return s

def baux_all(x, kmax=8):
    ep = math.exp(x); em = math.exp(-x); rx = 1.0/x
    b = [0.0]*(kmax+1)
    for k in range(kmax+1):
        term = rx
        sgn = 1.0 if (k % 2 == 0) else -1.0
        sp = sgn*term; sm = term
        for j in range(1, k+1):
            term = term * (k - j + 1) * rx
            sgn = -sgn
            sp += sgn*term
            sm += term
        b[k] = ep*sp - em*sm
    return b

print("=== 1) bint(x, k) 12-term truncation & |x|<1e-6 jump vs 50-digit mpmath ===")
for x in [9.9e-7, 1.0e-6, 1.0, 1.49, 2.0, 3.0, 5.0]:
    parts = []
    for k in [0, 1, 5, 6]:
        ex = B_exact(k, x); got = bint_12(k, x)
        rel = abs(got - ex) / abs(ex) if ex != 0 else abs(got - ex)
        parts.append(f"B{k}:{rel:.2e}")
    print(f"  x={x:8.2e}: " + "  ".join(parts))

print("\n=== 2) baux_all(x, b(0:8)) closed-form cancellation vs 50-digit mpmath ===")
for x in [0.05, 0.10, 0.1218, 0.15, 0.20, 0.50]:
    b = baux_all(x, 8)
    parts = []
    for k in [6, 7, 8]:
        ex = B_exact(k, x)
        rel = abs(b[k] - ex) / abs(ex)
        parts.append(f"B{k}:{b[k]:+.4e}(rel {rel:.2e})")
    print(f"  x={x:6.4f}: " + "  ".join(parts))
