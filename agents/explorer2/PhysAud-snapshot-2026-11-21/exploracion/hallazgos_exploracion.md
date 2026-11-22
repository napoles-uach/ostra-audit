# Exploración adicional multi-librería (XCFun 2.1.1, PySCF 2.14.0, libcint 6, Libxc 7.0.0)

Auditoría física exploratoria realizada bajo el protocolo obligatorio de seis pasos (Pasos 1–5 completados; Paso 6 pendiente de aprobación humana).

## Entorno y versiones medidas
- **Python**: 3.13 (`.venv` en `/Users/napolesd/Documents/Antigravity/PhysAud/.venv`)
- **PySCF**: `2.14.0` (incluye `libxc.15.dylib` = Libxc `7.0.0`, `libxcfun.2.dylib` = XCFun `2.1.1`, `libcint.6.dylib`, `libcgto.dylib`, `libpbc.dylib`, `libdft.dylib`)
- **Libxc devel** (para verificación de `ext_params`): compilación local `build-libxc/libxc.dylib` (commit `3c7794c50`)

---

## Resumen de hallazgos verificados (Pasos 1–5)

| # | Librería | Hallazgo | Clase | Oráculo | Impacto medido (Paso 5) | Reproductor |
|---|---|---|---|---|---|---|
| 1 | **XCFun 2.1.1** (`libxcfun.2.dylib`) | **`REVTPSSX` (`id=44`)**: triple error (registrado como GGA, denominador sin cuadrado, término `(p·σ)²` dimensional) | Error numérico nuevo | Escalamiento uniforme de intercambio + diferencial vs Libxc + contrato API | `NaN` si se evalúa solo; viola escalamiento uniforme y cota Lieb–Oxford; polo en `x = -κ = -0.804`; en H₂O (`def2-TZVP`): `ΔE_x = -3.786 Eh` (`-2,376 kcal/mol`, 42 %); SCF diverge a `+13.84 Eh` (`+90.3 Eh`) | [repro_01_xcfun_revtpssx.py](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion/reproducers/repro_01_xcfun_revtpssx.py) |
| 2 | **XCFun 2.1.1** (`libxcfun.2.dylib`) | **`REVTPSSC` (`id=43`)**: exponente `z²` en lugar de `z³` en Ec. (24) + `NaN` en `(s=0, α=0)` | Error numérico nuevo | Diferencial vs Libxc + análisis funcional en `z = τ_W/τ` | Sigiloso: pasa UEG (`z=0`), 1-electrón/He (`z=1`) y FD4 interno. Desplaza SCF (`TPSSx+revTPSSc`, `def2-TZVP`, grid 4): **H₂O `+1.939 kcal/mol`**, **N₂ `+2.136 kcal/mol`**, **CH₄ `+1.477 kcal/mol`** | [repro_02_xcfun_revtpssc.py](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion/reproducers/repro_02_xcfun_revtpssc.py) |
| 3 | **XCFun 2.1.1** (`libxcfun.2.dylib`) | **`R4SCANX` (`id=54`)**: parámetros `dp4 = da4 = 0.232` en lugar de `dp4 = 0.802, da4 = 0.178` | Error de parámetros / versión preprint | Diferencial vs Libxc + verificación con `xc_func_set_ext_params` (`3.7×10⁻¹⁵`) | Apaga el término de 4.º orden en `5×10⁻¹⁰` ya en `s = 0.5`, colapsando `r4SCAN` hacia `r2SCAN`. Desplaza SCF (`r4SCAN`, `def2-TZVP`, grid 4): **H₂O `+0.455 kcal/mol`**, **N₂ `+0.879 kcal/mol`**, **CH₄ `+0.929 kcal/mol`** | [repro_03_xcfun_r4scanx.py](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion/reproducers/repro_03_xcfun_r4scanx.py) |
| 4 | **XCFun 2.1.1** (`libxcfun.2.dylib`) | **`SCANX`/`SCANC` (`ids=46,45`)**: `NaN`/`-inf` en `α = 1`, `v_σ = NaN` en `s = 0` (`SCANC`) y salto del 49 % en `v_σ` en `s = 0` (`SCANX`) | Error numérico en la orilla | Frontera (`α = 1`, `s = 0`) + diferencial vs Libxc | En `α = 1` exacto (cualquier `s`): `SCANX → NaN`, `SCANC → -inf`; en `s = 0`: `SCANC v_σ = NaN` y `SCANX v_σ` salta 49 %. En mallas moleculares Becke (H₂O, Ne, H₂) se evita el punto exacto (`diff ~ 4×10⁻⁷ Eh`) | [repro_04_xcfun_scan_and_tfk.py](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion/reproducers/repro_04_xcfun_scan_and_tfk.py) |
| 5 | **XCFun 2.1.1** (`libxcfun.2.dylib`) | **`TFK` (`id=24`, Thomas–Fermi)**: ignora polarización de espín | Error físico / implementación | Relación exacta de escalamiento de espín $T_s[\rho_\alpha,\rho_\beta] = \frac{1}{2}(T_s[2\rho_\alpha]+T_s[2\rho_\beta])$ | Evalúa $C_F(\rho_\alpha+\rho_\beta)^{5/3}$ en lugar de $2^{2/3} C_F(\rho_\alpha^{5/3}+\rho_\beta^{5/3})$: error de **8.66 %** en $\zeta=0.41$ y **37.0 %** ($1-2^{-2/3}$) en $\zeta=1$ | [repro_04_xcfun_scan_and_tfk.py](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion/reproducers/repro_04_xcfun_scan_and_tfk.py) |
| 6 | **PySCF 2.14.0** (`pyscf/dft/xcfun.py`) | **`SIGSEGV` (`-11`) en los 19 funcionales VV10/rVV10** con backend `xcfun` | Contrato de API | Contrato de API (ejecución de todos los alias de `xcfun.XC_CODES`) | `xcfun.py:368` inyecta IDs ficticios `5000..5018` para `VV10_XC`, pero `eval_xc1` (`xcfun.py:1037-1061`) los envía sin filtrar a C `_itrf.XCFUN_eval_xc` → muerte instantánea por `SIGSEGV` en `wb97m_v`, `b97m_v`, `scan_rvv10`, `pbe+vv10`, etc. | [repro_05_pyscf_xcfun_vv10_segfault.py](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion/reproducers/repro_05_pyscf_xcfun_vv10_segfault.py) |
| 7 | **Libxc 7.0.0 × PySCF 2.14.0 TDDFT** | **Impacto observable de `s_scaling_2` (`s=1`) en gradientes analíticos TDDFT** (`HSE06`, `LC-ωPBE`, `GGA_X_WPBEH`) | Impacto Paso 5 de discontinuidad en `s = 1` | Consistencia fuerza analítica vs FD4 (escaneo de `h` y nivel de malla 4, 6, 8) | En CO (`def2-SVP`), `CAM-B3LYP` converge a `2.77×10⁻⁹ Eh/Bohr` (malla 8) y `GGA_X_HJS_PBE` a `6.20×10⁻⁷` (malla 6), pero **`HSE06` (`3.38×10⁻⁵`)**, **`LC-ωPBE` (`3.94×10⁻⁵`)** y **`GGA_X_WPBEH` (`2.32×10⁻⁵`)** no convergen al refinar la malla (error **12,000–14,000× mayor** que `CAM-B3LYP`) | [repro_06_tddft_wpbeh_grad.py](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion/reproducers/repro_06_tddft_wpbeh_grad.py) |

---

## Detalle técnico y verificación de causas raíz (Pasos 3 y 4)

### 1. XCFun 2.1.1 — `REVTPSSX` (`id = 44`)
- **Evidencia medida:**
  1. `xcfun.xc_type("REVTPSSX") == "GGA"` (el registro C devuelve `1` en lugar de `2` = MGGA). Al llamarse solo, devuelve `exc = NaN, vrho = NaN, vsigma = NaN`.
  2. Al emparejarlo con `REVTPSSC` (forzando modo MGGA) y restar `REVTPSSC`, el factor de realce $F_x(s, \alpha)$ depende fuertemente de $\rho$, violando el escalamiento uniforme exacto del intercambio: para $s=0.5, \alpha=1.0$, $F_x = 1.0410$ ($\rho=0.1$), $1.0086$ ($\rho=1.0$), $-2.3831$ ($\rho=3.0$) y $+5.0486$ ($\rho=3.5$), atravesando un **polo** en $x = -\kappa = -0.804$ ($\rho \approx 3.23$) y violando la cota de Lieb–Oxford ($1+\kappa = 1.804$).
  3. Inversión algebraica exacta contra la Ec. (10) de Perdew et al., *PRL* **103**, 026403 (2009) (con $c = 2.35204, e = 2.1677, b = 0.40, \kappa = 0.804, \mu = 0.14$): la evaluación de XCFun coincide a **`4.44×10⁻¹⁵`** en 11 puntos de prueba cuando se introducen exactamente dos errores simultáneos en $x(p, z)$:
     - Denominador `1 + sqrt(e) * p` en lugar de `(1 + sqrt(e) * p)^2`.
     - Dentro de la raíz $-\frac{73}{405}\tilde{q}_b \sqrt{\frac{1}{2}(\frac{9}{25} z^2 + p^2)}$, $p^2$ fue sustituido por `(p * sigma)^2`, donde $\sigma = |\nabla\rho|^2 = 4 g_{\alpha\alpha} \propto \rho^{8/3}$ es dimensional.
- **Verificación upstream:** No se localizaron issues abiertos/cerrados en `dftlibs/xcfun` que reporten este error en búsqueda web (nota: no se dispone de clon local de `dftlibs/xcfun`, por lo que conviene verificar el historial del repositorio antes de abrir un issue).

### 2. XCFun 2.1.1 — `REVTPSSC` (`id = 43`)
- **Evidencia medida:**
  - Cumple consistencia interna FD4, coincide en UEG ($z = 0$) y en el límite de 1 electrón ($z = 1$, donde $z^2 = z^3 = 1$).
  - Sin embargo, en 120 puntos con $0 < z < 1$ (no polarizados y polarizados), la discrepancia frente a Libxc `MGGA_C_REVTPSS` cae de **`2.44×10⁻³`** a **`1.37×10⁻⁷`** (atribuible al redondeo `beta_a = 0.066725` del artículo frente a `0.0667245506` en Libxc) al sustituir $z^3 \to z^2$ en la Ec. (24): $\varepsilon_c = \varepsilon_c^{\text{PKZB}}(1 + d\,\varepsilon_c^{\text{PKZB}} z^3)$ con $d = 2.8$.
  - Además devuelve `NaN` en `(s = 0, α = 0)` (división por $\tau = 0$).
  - Por su propiedad de anularse en $z = 0$ y $z = 1$, pasa pruebas estándar de UEG, átomo de H y átomo de He, pero desplaza silenciosamente la energía de moléculas reales en **`+1.48` a `+2.14 kcal/mol`**.

### 3. XCFun 2.1.1 — `R4SCANX` (`id = 54`)
- **Evidencia medida:**
  - En Furness et al., *J. Chem. Phys.* **156**, 034109 (2022), los parámetros de amortiguamiento de cuarto orden son $d_{p4} = 0.802, d_{\alpha 4} = 0.178$ (verificados también como valores por defecto en Libxc `mgga_x_r4scan`).
  - Ajustando `_dp4 = 0.232` y `_da4 = 0.232` en Libxc mediante `xc_func_set_ext_params`, Libxc coincide con XCFun `R4SCANX` a **`3.72×10⁻¹⁵`** en $\varepsilon_x$ y **`1.4–8.6×10⁻¹⁴`** en $v_\rho, v_\sigma, v_\tau$ sobre 10 puntos diversos.
  - Debido a que $0.232^4 = 2.90\times 10^{-3}$ frente a $0.802^4 = 0.4137$, el factor $\exp(-s^4 / d_{p4}^4)$ cae a $5\times 10^{-10}$ ya en $s = 0.5$, suprimiendo casi por completo la corrección de cuarto orden en regiones moleculares típicas (`+0.46` a `+0.93 kcal/mol` en H₂O, N₂, CH₄).
  - Nota (Paso 3): No se pudo verificar directamente el texto completo del preprint inicial (arXiv v1) ni los commits de `dftlibs/xcfun` para confirmar si `0.232` provino de una versión preliminar del manuscrito o de un error de transcripción; en cualquier caso, discrepa silenciosamente de los parámetros publicados de `r4SCAN`.

### 4. XCFun 2.1.1 — `SCANX` / `SCANC` y `TFK`
- **`SCANX`/`SCANC`**: En $\alpha = 1.0$ exacto (cualquier $s$), la función de conmutación $f_{x,c}(\alpha)$ contiene $1/(1-\alpha)$ sin guarda, devolviendo `NaN` (`SCANX`) y `-inf` (`SCANC`). En $\alpha = 1 + 10^{-10}$, `SCANX v_sigma` colapsa a `+2.02×10⁻¹²` frente a `-5.13×10⁻²` en Libxc. En $s = 0$ ($\sigma = 0$), `SCANC` devuelve `v_sigma = NaN` para todo $\alpha$, y `SCANX v_sigma` cae un 49 % (`-5.08×10⁻²` vs `-9.93×10⁻²` en $\alpha = 1.2$).
- **`TFK` (`id = 24`)**: Coincide con `Libxc LDA_K_TF` a `5.5×10⁻¹⁶` sin polarizar, pero al pasar densidades polarizadas $(\rho_\alpha, \rho_\beta)$ devuelve exactamente el valor no polarizado evaluado en $\rho_\alpha + \rho_\beta$, violando $T_s[\rho_\alpha,\rho_\beta] = \frac{1}{2}(T_s[2\rho_\alpha]+T_s[2\rho_\beta])$ hasta en un **37.0 %** en $\zeta = 1$. (`VWK`, `PW91K` y `BTK` en el mismo XCFun sí cumplen el escalamiento de espín a precisión de máquina.)

### 5. PySCF 2.14.0 — `pyscf/dft/xcfun.py` `SIGSEGV` en VV10/rVV10
- **Líneas exactas:**
  - [xcfun.py](file:///Users/napolesd/Documents/Antigravity/PhysAud/.venv/lib/python3.13/site-packages/pyscf/dft/xcfun.py#L368): `XC_CODES.update([(key, 5000+i) for i, key in enumerate(VV10_XC)])`
  - [xcfun.py](file:///Users/napolesd/Documents/Antigravity/PhysAud/.venv/lib/python3.13/site-packages/pyscf/dft/xcfun.py#L1037-L1061): `eval_xc1` construye `fn_ids = [x[0] for x in fn_facs]` y los pasa sin filtrar `xid >= 5000` directamente a `_itrf.XCFUN_eval_xc(...)`, causando `SIGSEGV` (`returncode = -11`) en los 19 alias VV10/rVV10.

### 6. Libxc 7.0.0 × PySCF 2.14.0 — Impacto observable de `s_scaling_2` (`s = 1`) en gradientes TDDFT
- En CO (`def2-SVP`, `z = 2.10 Bohr`, `h = 10⁻³ Bohr`), comparando gradiente analítico TDDFT (`grid_response=True`, estado 1) contra diferencias finitas de 4 puntos (FD4):
  - `CAM-B3LYP`: `7.40×10⁻⁸` (malla 4) → `3.68×10⁻⁸` (malla 6) → **`2.77×10⁻⁹ Eh/Bohr`** (malla 8).
  - `GGA_X_PBE,GGA_C_PBE`: `5.45×10⁻⁷` (malla 6); `GGA_X_HJS_PBE,GGA_C_PBE` (otro PBE de corto alcance, sin `s_scaling_2`): `6.20×10⁻⁷` (malla 6).
  - `GGA_X_WPBEH,GGA_C_PBE` (GGA puro con `s_scaling_2`): **`2.32×10⁻⁵ Eh/Bohr`** (malla 6).
  - `HSE06`: `3.74×10⁻⁵` (malla 4) → `1.74×10⁻⁵` (malla 6) → **`3.38×10⁻⁵ Eh/Bohr`** (malla 8; `g_fd4` converge a `4.7×10⁻⁷` entre mallas 6 y 8, pero `g_ana` oscila y se queda desplazado).
  - `LC-ωPBE`: `1.85×10⁻⁵` (malla 4) → `1.43×10⁻⁵` (malla 6) → **`3.94×10⁻⁵ Eh/Bohr`** (malla 8; `g_fd4` converge a `4.1×10⁻⁸`).
- **Conclusión Paso 5:** La discontinuidad de `s_scaling_2` en `s = 1` (presente en `gga_x_wpbeh`, compartido por `HSE03/06/12` y `LC-ωPBE`) degrada en **4 órdenes de magnitud** la consistencia de gradientes analíticos TDDFT en estados excitados.

---

## Controles negativos y señales descartadas por el protocolo (Paso 2 / Paso 5)
- **PySCF PBC (GDF y RSDF) en celdas unimodulares sesgadas (`Si` primitivo, `gth-szv`, `gth-pbe`)**: Ewald (`0.0`), `RHF(GDF)` (`5.3×10⁻¹⁰ Eh`), `E_coul` y `E_1` son completamente invariantes. El desvío inicial de `-6.06×10⁻⁵ Eh` en `RKS` (`skew2`) proviene únicamente de `BeckeGrids` (desciende a `-7.76×10⁻⁶ Eh` al subir a `Becke level 5`). A diferencia de `multicharge`, PySCF PBC GDF/RSDF **no** presenta error de imágenes fijas.
- **libcint 6 `int2e_erfc` (integrales de corto alcance, `l = 5`)**: Para exponentes sintéticos (`α = 0.5, R = 0.7, ω = 20`), `lrys_wheeler` imprime `singular value` y `int2e_erfc` explota a `2.63×10⁶`. Sin embargo, en pruebas con base real (`Ne2 cc-pV5Z`, capas `H` en dos centros, `R = 2.5 Bohr`, `ω = 0.11..20.0`), `LR + SR = Full` se mantiene exacto a `≤ 8.6×10⁻¹¹`.
- **PySCF VV10 (`r2SCAN+VV10`, `ωB97M-V`, `B97M-V`) y solvente implícito (`ddCOSMO`, `IEF-PCM`, `C-PCM`, `COSMO`, `SS(V)PE`)**: gradientes analíticos vs FD4 limpios (`6×10⁻¹⁰` a `1.5×10⁻⁹ Eh/Bohr`).
- **PySCF `libcgto` ECP (`Au-Cl`, `def2-TZVP`)**: diferencia `ECPscalar_ipnuc` vs FD4 de `1.8×10⁻⁶` en `h = 10⁻⁴` que baja a `1.5×10⁻⁹` en `h = 10⁻⁵` (cambio discreto de orden en la cuadratura radial ECP a nivel `~10⁻¹⁰`).
- **XCFun `M05X` / `M052XX` vs Libxc**: cocientes exactos `0.72` y `0.44` (`1 - a_HF`) por convención de normalización de híbridos (no es un error).
