# Respuesta al pendiente de la Tanda 1: banda de `NaN` en `lda_c_pk09` cerca de $r_s \approx 22.739$

## Por qué la verificación independiente no lo reprodujo

La verificación independiente probó `lda_c_pk09` **sin espín (`XC_UNPOLARIZED`)** y **polarizado con $\zeta = 0$**. En ambos casos la expresión que desborda **no se activa**, porque en `/Users/napolesd/Documents/libxc-devel/maple/lda_exc/lda_c_pk09.mpl:36`:

```maple
alpha_z := (rs, z) -> 2/((1 + z)^ss(rs) + (1 - z)^ss(rs)):
```

Cuando $\zeta = 0$ (o en la rama no polarizada, donde $\zeta \equiv 0$), $(1 + 0)^{s_s} = 1^{s_s} = 1$ para cualquier $s_s$ finito — incluso cuando $|s_s| \sim 10^4$ o $10^{10}$! En cambio, cuando **`XC_POLARIZED` (`spin = 2`) y $0 < |\zeta| < 1$** (por ejemplo $\zeta = 0.4$), $|s_s| \gtrsim 709 / \ln(1 + |\zeta|)$ provoca el desbordamiento a $+\infty$ de $(1 + |\zeta|)^{s_s}$ (o $(1 - |\zeta|)^{s_s}$ si $s_s < 0$), lo que genera `inf / inf` o `0 * inf` $\to$ `NaN` tanto en `zk` ($\varepsilon_c$) como en `vrho`.

## Condiciones exactas de reproducción

- **Versión y commit probados:** `libxc` rama `devel` (reporta versión `7.1.2`, commit `3c7794c50`), compilado en `/Users/napolesd/Documents/Antigravity/PhysAud/build-libxc/libxc.dylib`.
- **Funcional:** `lda_c_pk09` (`XC_LDA_C_PK09`, id `554`).
- **Polarización de espín:** **`XC_POLARIZED` (`spin = 2`)** con **$0 < |\zeta| < 1$** (p. ej. $\zeta = 0.4$, $\rho_\uparrow = 0.7\,\rho_{\text{tot}}$, $\rho_\downarrow = 0.3\,\rho_{\text{tot}}$).
- **Causa exacta en código fuente:** `maple/lda_exc/lda_c_pk09.mpl:22-25` protege el denominador `f_r_den(rs)` de `ss(rs)` con `DETOL = 1e-12`, pero el polo real $r_{s,\text{polo}} \approx 22.739412132480812$ ($\rho_{\text{tot}} = 3 / (4\pi r_s^3) \approx 2.028219 \times 10^{-5}$) tiene una derivada enorme $f_{\text{den}}'(r_{s,\text{polo}}) \approx 2.74 \times 10^3$. Por tanto, `|f_r_den(rs)| < 1e-12` solo se dispara para $|r_s - r_{s,\text{polo}}| < 3.6 \times 10^{-16}$. En toda una banda de anchura $\Delta r_s \sim 10^{-5}$ alrededor del polo, `f_r_den` no activa `DETOL`, $|s_s|$ supera $10^3$–$10^4$, y $(1 \pm \zeta)^{s_s}$ desborda a `Inf` $\to$ `NaN`.

## Salida reproducida (`reproducers/repro_06_secondary_lak_pk09_lb94_br89.py`)

Con `LIBXC_PATH=~/Documents/Antigravity/PhysAud/build-libxc/libxc.dylib`, `spin = 2`, $\zeta = 0.4$, $r_{s,\text{polo}} = 22.739412132480812$:

| $\Delta r_s = r_s - r_{s,\text{polo}}$ | $\rho_{\text{tot}}$ | $\rho_\uparrow$ | $\rho_\downarrow$ | `zk` ($\varepsilon_c$) | `vrho[0]` | `vrho[1]` |
|---:|---:|---:|---:|---:|---:|---:|
| `-1.0e-4` | `2.028246e-05` | `1.419772e-05` | `6.084737e-06` | `-1.287066e+03` | `+1.402124e+09` | `+1.402124e+09` |
| **`-1.0e-6`** | `2.028219e-05` | `1.419753e-05` | `6.084658e-06` | **`nan`** | **`nan`** | **`nan`** |
| **`+0.0e+0`** | `2.028219e-05` | `1.419753e-05` | `6.084657e-06` | **`nan`** | **`nan`** | **`nan`** |
| **`+1.0e-6`** | `2.028219e-05` | `1.419753e-05` | `6.084656e-06` | **`nan`** | **`nan`** | **`nan`** |
| `+1.0e-4` | `2.028192e-05` | `1.419735e-05` | `6.084577e-06` | `-2.681227e-02` | `-1.664180e-02` | `-1.923443e-02` |

*(Nótese además la divergencia física no-NaN justo a la izquierda del polo en $\Delta r_s = -10^{-4}$: `zk = -1.287e+03`, `vrho = +1.402e+09`, causada por el cambio de signo de `ss(rs)` a negativo al cruzar el polo.)*
