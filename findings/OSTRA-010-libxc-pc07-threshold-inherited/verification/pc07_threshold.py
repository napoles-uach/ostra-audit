"""PC07 switching function (maple/mgga_exc/mgga_k_pc07.mpl, eq. 15) at the clamping points z = a/40 and 39a/40."""
import mpmath as mp
mp.mp.dps = 30
fab0 = lambda z, a, b: mp.exp(-a*b/z) * (1 + mp.exp(-a/(a - z)))**b / (mp.exp(-a/z) + mp.exp(-a/(a - z)))**b
for a, b, lab in ((0.5389, 3, "PC07"), (1.784720, 0.258304, "PC07-opt (SCAN-L, r2SCAN-L)")):
    print("%-28s f(a/40)=%.2e   1-f(39a/40)=%.2e" % (lab, fab0(a/40, a, b), 1 - fab0(a*39/40, a, b)))
