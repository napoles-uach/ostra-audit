# OSTRA-001: Tanh-sinh quadrature tail check in xc_integrate truncates integrals with mass away from the midpoint

| Field | Value |
|---|---|
| library | libxc |
| location | src/integrate.c (xc_integrate) |
| versions | 7.1.2 and devel before a19650f34 |
| detected | 2026-11-19 |
| detected_by | Explorer 1 (Claude, main session) |
| oracle | Analytic reference (closed-form integrals of 2K0, 2tK0) |
| libraries | libxc |
| pattern | convergence-fixed |
| oracle_class | analytic-reference |
| class | New numerical bug |
| impact_level | III (niche: lda_x_1d_soft for beta*kF >~ 38.6) |
| status | Resolved upstream |
| report | Fixed upstream: [libxc#646](https://gitlab.com/libxc/libxc/-/issues/646) (commit a19650f34) |

## Symptom
The exchange energy per particle of `lda_x_1d_soft` collapses and changes sign for beta*kF >~ 38.6.

## Cause
The tail check of the tanh-sinh quadrature stops when the weight falls below tolerance even if the integrand
still carries mass far from the interval midpoint. Fix (adopted upstream):
`if(w <= tol && w * (fabs(pts[0]) + fabs(pts[1])) <= tol * (1.0 + fabs(T)))`.

## Verification
Integrals checked against closed forms (int 2K0 -> pi; int 2tK0 = 2(1 - R K1(R))); fix verified on devel after merge.
