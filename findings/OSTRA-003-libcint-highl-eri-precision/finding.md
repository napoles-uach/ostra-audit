# OSTRA-003: Four-centre ERIs lose precision at high angular momentum (Rys roots)

| Field | Value |
|---|---|
| library | libcint |
| location | src/rys_roots.c |
| versions | PySCF 2.15.1 bundled libcint |
| detected | 2026-11-20 |
| detected_by | Explorer 1 |
| oracle | Invariances and analytic references |
| class | Documented limitation |
| impact_level | IV (Delta E_HF = 1e-13 Eh with cc-pV6Z) |
| status | Documented; not reported |
| report | none (documented by libcint) |

Absolute errors ~1e-6 (l=6) and ~1e-4 (l=7) on small elements; propagated to HF energies the effect is 1e-13 Eh.
