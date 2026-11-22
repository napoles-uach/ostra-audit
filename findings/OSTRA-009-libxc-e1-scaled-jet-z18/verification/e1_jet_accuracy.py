"""Python transcription of xc_E1_scaled_jet's asymptotic branch (src/expint_e1.c, z > 18) vs mpmath."""
import mpmath as mp
mp.mp.dps = 40
def asym(z, n):
    out = []
    for k in range(n + 1):
        ssum, prev, term = 0.0, float("inf"), 1.0
        for m in range(1, k + 1): term *= m
        for m in range(k + 1): term /= z
        for m in range(60):
            if term > prev: break
            ssum += -term if m & 1 else term; prev = term; term *= (m + k + 1) / z
        out.append(-ssum if k & 1 else ssum)
    return out
ref = lambda z, k: mp.diff(lambda x: mp.exp(x) * mp.e1(x), z, k)
for z in (18.000001, 20, 25, 30, 40):
    a = asym(z, 3)
    print("z=%-10g " % z + "  ".join("k%d rel %.1e" % (k, abs((a[k] - ref(mp.mpf(z), k)) / ref(mp.mpf(z), k))) for k in range(4)))
