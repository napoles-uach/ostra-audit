# Verificación independiente (a ciegas) de OSTRA-027 (`TW`) y OSTRA-028 (`SPBEC`) en XCFun

**Fecha:** 2026-11-26  
**Versión auditada:** `github.com/dftlibs/xcfun`, commit `a89b783` (`master`, v2.1.1; misma base empaquetada en PySCF 2.14.0/2.15.1)  
**Entorno:** macOS arm64, `Apple clang++ 21.0.0` (`-std=c++11 -O2 -DXCFUN_MAX_ORDER=3`), Python 3.13, PySCF 2.14.0 (XCFun 2.1.1 + libxc 7.0.0)  
**Metodología:** Fase 1 ejecutada **completamente a ciegas** sin leer `ostra-audit` (lectura de código fuente, derivación desde primeros principios, oráculos físicos/analíticos, compilación de parches mínimos y verificación bit a bit de los otros 80 funcionales de XCFun). Fase 2 realizada solo tras congelar los resultados de Fase 1 para contrastar con `OSTRA-027` y `OSTRA-028`.

---

## Resumen ejecutivo de veredictos (tres niveles de evidencia)

| Hallazgo | Anomalía reproducida | Causa identificada | Defecto confirmado (parche mínimo probado) | Veredicto global frente a `ostra-audit` |
|---|---|---|---|---|
| **A — OSTRA-027** (`TW`, `src/functionals/tw.cpp:21`) | **Confirmada** (H 1s exacto: `0.007505 Eh` vs `0.5 Eh`; H$_2$O: `1.318e7 Eh` vs `57.283 Eh`, 230,145×) | **Confirmada** (`pow(d.gaa + d.gbb, 2.0)` eleva al cuadrado gradientes ya al cuadrado y omite `2*d.gab`) | **Confirmado** (parche compilado + identidad bit a bit 80/80 otros funcionales) | **Confirmado con corrección importante al parche de OSTRA-027**: el parche `1/8*(gaa+2*gab+gbb)/n` en `OSTRA-027` **no** coincide con `libxc GGA_K_VW` en sistemas abiertos parcialmente polarizados (difiere en **-28.71 kcal/mol** en el átomo de Li) porque incumple el escalamiento de espín de Oliver–Perdew; además XCFun ya tiene un funcional duplicado y correcto (`VWK` en `vonw.cpp`). |
| **B — OSTRA-028** (`SPBEC`, `src/functionals/spbec.cpp:22-27`) | **Confirmada** (H$_2$O: `-0.595478 Eh` vs `-0.362040 Eh`, error `-146.48 kcal/mol` / `-64.48%`; viola el límite exacto $t \to \infty$) | **Confirmada** (`gamm` y `beta` intercambiados en L22–23 **y** `beta_gamma` en lugar de `gamm` dentro de `expm1` en L27) | **Confirmado** (parche compilado + identidad bit a bit 80/80 otros funcionales) | **Confirmado con matiz sobre VWN5 vs PW92 y evidencia adicional**: 99.0% del error (`-144.99 kcal/mol`) y la violación del límite $t \to \infty$ vienen de L22–27; el 1.0% restante (`-1.49 kcal/mol`) es VWN5 (convención ADF, mantenida en `OSTRA-028/report/spbec.patch`) vs PW92 (`libxc GGA_C_SPBE`, igualado a `1.67e-15 Eh` en nuestro parche). |

---

## Fase 1 (a ciegas): análisis, oráculos y parches mínimos

### Caso A: XCFun `TW` (`src/functionals/tw.cpp`)

#### 1. Lectura del código y derivación desde primeros principios
En `src/functionals/tw.cpp:20-22` (commit `a89b783`):
```cpp
template <typename num> static num tw(const densvars<num> & d) {
  return 1. / 8. * pow(d.gaa + d.gbb, 2.0) / d.n;
}
```
Según `src/densvars.hpp:59-62, 158-160`, las variables de `densvars` ya representan productos escalares de gradientes:
- `d.gaa` $= \sigma_{\alpha\alpha} = |\nabla n_\alpha|^2$,
- `d.gab` $= \sigma_{\alpha\beta} = \nabla n_\alpha \cdot \nabla n_\beta$,
- `d.gbb` $= \sigma_{\beta\beta} = |\nabla n_\beta|^2$,
- `d.gnn` $= |\nabla n|^2 = \sigma_{\alpha\alpha} + 2\sigma_{\alpha\beta} + \sigma_{\beta\beta}$.

A partir de la definición de la energía cinética de von Weizsäcker existen dos expresiones relevantes:
1. **Forma de densidad total (sin escalamiento de espín):**
   $$T_W[n] = \frac{1}{8} \int \frac{|\nabla n|^2}{n}\,d^3r = \frac{1}{8}\int \frac{\text{d.gnn}}{d.n}\,d^3r$$
2. **Forma con escalamiento exacto de espín** (Oliver–Perdew 1979, $T_s[n_\alpha, n_\beta] = \frac{1}{2}T_s[2n_\alpha] + \frac{1}{2}T_s[2n_\beta]$, usada por `libxc` `GGA_K_VW` y por el otro funcional de von Weizsäcker presente en XCFun, `src/functionals/vonw.cpp` / `XC_VWK`):
   $$T_W[n_\alpha, n_\beta] = \frac{1}{8}\int \left(\frac{|\nabla n_\alpha|^2}{n_\alpha} + \frac{|\nabla n_\beta|^2}{n_\beta}\right)d^3r = \frac{1}{8}\int\left(\frac{\text{d.gaa}}{d.a} + \frac{\text{d.gbb}}{d.b}\right)d^3r$$

Ambas coinciden exactamente en sistemas no polarizados ($n_\alpha = n_\beta$) y en sistemas de un solo electrón ($n_\beta = 0$), pero difieren en densidades parcialmente polarizadas con $\nabla n_\alpha/n_\alpha \neq \nabla n_\beta/n_\beta$ (por Cauchy–Schwarz, $T_W[n] \le T_W[n_\alpha, n_\beta]$).

**Errores en `src/functionals/tw.cpp:21`:**
- **Error 1 (catastrófico, dimensionalmente incorrecto):** `pow(d.gaa + d.gbb, 2.0)` eleva al cuadrado cantidades que ya están al cuadrado, evaluando $\frac{1}{8}(|\nabla n_\alpha|^2 + |\nabla n_\beta|^2)^2 / n$ (dimensión $L^{-11}$ en vez de $L^{-5}$).
- **Error 2:** omite el término cruzado `2 * d.gab` (incluso eliminando `pow(..., 2.0)`, `(d.gaa + d.gbb)/(8*d.n)` sería exactamente la mitad del valor correcto en cualquier densidad cerrada no polarizada).
- **Error 3 (escalamiento de espín):** no aplica el escalamiento de espín de Oliver–Perdew, a diferencia de `libxc GGA_K_VW` y de `XC_VWK` (`src/functionals/vonw.cpp:17-23`) dentro de la propia XCFun.
- Adicionalmente, `tw.cpp:29` no tiene datos de referencia de prueba (`xcfun_test()` reporta `TW has no test`).

#### 2. Oráculos exactos y resultados numéricos (baseline)
Script: [verify_tw.py](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/scripts/verify_tw.py) | Salidas: [verify_tw_baseline.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/outputs/verify_tw_baseline.log), [verify_tw_li_openshell.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/outputs/verify_tw_li_openshell.log).

Para la densidad exacta $1s$ del átomo de hidrógeno ($\rho_\alpha(r) = \frac{1}{\pi} e^{-2r}$, $\rho_\beta = 0$, $|\nabla\rho_\alpha| = 2\rho_\alpha$), se tienen dos integrales analíticas exactas:
- Exacta: $T_W = \frac{1}{8}\int \frac{4\rho^2}{\rho}\,d^3r = 0.5\text{ E}_h$.
- Integrando erróneo de `tw.cpp`: $\frac{1}{8}\int \frac{(4\rho^2)^2}{\rho}\,d^3r = 2\int \rho^3\,d^3r = \frac{2}{27\pi^2} \approx 0.007505272862395\text{ E}_h$ (66.6× demasiado pequeño).

| Sistema / Punto | Oráculo exacto / `libxc GGA_K_VW` | XCFun `'VWK'` (`vonw.cpp`) | XCFun `'TW'` (`tw.cpp`, baseline) | Error / Factor de `'TW'` |
|---|---|---|---|---|
| H $1s$ exacto (malla nivel 5, $T_W = 0.5\text{ E}_h$) | `0.499999999983 Eh` (`err -1.73e-11`) | `0.499999999990 Eh` | **`0.007505272862 Eh`** ($= \frac{2}{27\pi^2}$ a `8.7e-18`) | `-0.492495 Eh` (**66.6× menor**) |
| H UHF/cc-pVQZ en PySCF ($\langle T \rangle = 0.499944788\text{ E}_h$) | `0.499944788003 Eh` (`diff -8.5e-12`) | `0.499944788012 Eh` | **`0.007481252349 Eh`** | `-0.492464 Eh` (**66.8× menor**) |
| H$_2$O RKS PBE/cc-pVTZ (cerrada) | `57.282882450107 Eh` | `57.282882450106 Eh` | **`1.318335556586e+07 Eh`** | **230,145× mayor (>5 órdenes)** |
| Punto no polarizado ($n=0.5, |\nabla n|^2=1.0$) | `0.250000000000` | `0.250000000000` | **`0.062500000000`** | 4× menor (sin `pow`: 2× menor) |
| Punto polarizado ($n_\alpha=0.35, n_\beta=0.15$) | `0.177380952381` (spin-scaled) | `0.177380952381` | **`0.024025000000`** (nota: `gnn/(8n)=0.072500`) | 7.4× menor |
| Átomo abierto Li (UKS PBE/cc-pVTZ) | `7.178192879667 Eh` | `7.178192879659 Eh` | **`4.141229471299e+02 Eh`** (nota: `gnn/(8n)=7.132449 Eh`) | **57.7× mayor** |

---

### Caso B: XCFun `SPBEC` (`src/functionals/spbec.cpp`)

#### 1. Lectura del código y fórmula publicada (Swart, Solà y Bickelhaupt, *J. Chem. Phys.* **131**, 094103, 2009)
En `src/functionals/spbec.cpp:21-45`:
```cpp
/* These are the constants in the paper, not very high precision.. */
static const parameter gamm = 0.066725;
static const parameter beta = 0.031091;
static const parameter beta_gamma = beta / gamm;

template <typename num> static num G(const num & eps, const num & phi3) {
  return beta_gamma / expm1(-eps / (beta_gamma * phi3));
}

template <typename num>
static num H_spbe(const num & t2, const num & eps, const num & phi3) {
  return gamm * phi3 * log(1 + beta_gamma * t2 / (1 + t2 * G(eps, phi3)));
}
...
template <typename num> static num spbec(const densvars<num> & d) {
  num eps = vwn::vwn5_eps(d);
  num p = phi(d);
  num t2 = (cbrt(M_PI / 3) / 16) * d.gnn * d.n_m13 / pow2(p * d.n);
  return d.n * (eps + H_spbe(t2, eps, pow3(p)));
}
```
La correlación sPBE publicada (Swart, Solà y Bickelhaupt 2009, basada en PBE 1996 eliminando el término cuadrático $A t^2$ del numerador y $A^2 t^4$ del denominador) es:
$$E_c^{\text{sPBE}} = \int n \left[\varepsilon_c^{\text{unif}}(n,\zeta) + H_{\text{sPBE}}(t^2, \varepsilon_c^{\text{unif}}, \phi)\right]\,d^3r$$
$$H_{\text{sPBE}} = \gamma \phi^3 \ln\!\left(1 + \frac{\beta}{\gamma}\,\frac{t^2}{1 + A t^2}\right), \qquad A = \frac{\beta}{\gamma}\left[\exp\!\left(-\frac{\varepsilon_c^{\text{unif}}}{\gamma \phi^3}\right) - 1\right]^{-1}$$
con las constantes de PBE: $\gamma = \frac{1-\ln 2}{\pi^2} \approx \mathbf{0.03109069087}$ (redondeado `0.031091`) y $\beta \approx \mathbf{0.06672455060}$ (redondeado `0.066725`), y $t^2 = \frac{1}{16}(\pi/3)^{1/3}\frac{|\nabla n|^2}{\phi^2 n^{7/3}}$ (verificado algebraicamente: L43 y `phi(d)` en L37 son correctas).

**Discrepancias identificadas en `src/functionals/spbec.cpp`:**
1. **Líneas 22–23 (constantes intercambiadas — defecto inequívoco):** `gamm = 0.066725` contiene $\beta$ y `beta = 0.031091` contiene $\gamma$. Esto invierte $\beta/\gamma$ ($0.4660$ en vez de $2.1461$, factor 4.61×) y multiplica por $2.15\times$ el prefactor $\gamma\phi^3$.
2. **Línea 27 (error de fórmula en el argumento de `expm1` — defecto inequívoco):** `expm1(-eps / (beta_gamma * phi3))` usa `beta_gamma` en vez de `gamm`. Físicamente esto **viola la restricción exacta de PBE/sPBE en el límite de variación rápida ($t \to \infty$)**:
   $$\lim_{t\to\infty} H_{\text{sPBE}} = \gamma\phi^3 \ln\!\left(1 + \frac{\beta/\gamma}{A}\right) = -\varepsilon_c^{\text{unif}} \implies \lim_{t\to\infty} \varepsilon_c^{\text{sPBE}} = 0.$$
   Con `beta_gamma` dentro de `expm1`, $\lim_{t\to\infty} H_{\text{sPBE}} = -\frac{\gamma^2}{\beta}\varepsilon_c^{\text{unif}} \neq -\varepsilon_c^{\text{unif}}$, por lo que la energía de correlación nunca se anula al crecer $|\nabla n|$.
3. **Línea 41 y precisión de constantes (convención LDA VWN5 vs PW92):** `spbec.cpp:41` emplea `vwn::vwn5_eps(d)` (convención ADF, análoga a `XC_VWN_PBEC` en `pbec.cpp`), mientras que `libxc GGA_C_SPBE` y el estándar PBE emplean `pw92eps::pw92eps(d)` y constantes de alta precisión (`xcfun_constants::param_gamma`, `param_beta_accurate`).

#### 2. Oráculo independiente y ablación de errores (baseline)
Script: [verify_spbec.py](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/scripts/verify_spbec.py) | Salida: [verify_spbec_baseline.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/outputs/verify_spbec_baseline.log).

| Evaluación / Variante | H$_2$O (RKS PBE/cc-pVTZ) | Diferencia vs `libxc GGA_C_SPBE` | Punto polarizado ($n_\alpha=0.35, n_\beta=0.15, \sigma_{nn}=0.29$) |
|---|---|---|---|
| **XCFun `'SPBEC'` (baseline)** | **`-0.595477895157 Eh`** | **`-0.233437 Eh` (`-146.48 kcal/mol`, `-64.48%`)** | **`-0.060204412236421`** (`-6.35%`; $v_\sigma$: `0.00281` vs `0.00824`) |
| Nuestra réplica exacta del bug (`spbec.cpp`) | `-0.595477895157 Eh` | (diff vs XCFun: `+8.06e-14 Eh`) | `-0.060204412236421` (diff vs XCFun: `-1.39e-17`) |
| sPBE propia (VWN5, constantes 6 dígitos, L22–27 corregidas) | `-0.364418585787 Eh` | `-0.002378 Eh` (`-1.49 kcal/mol`) | `-0.056969382478512` (`-3.59e-4`) |
| sPBE propia (PW92, constantes 6 dígitos, L22–27 corregidas) | `-0.362039114221 Eh` | `+1.378e-06 Eh` (`+0.00086 kcal/mol`) | — |
| sPBE propia (PW92, constantes exactas PBE, L22–27 corregidas) | **`-0.362040491755 Eh`** | **`+0.000e+00 Eh`** | **`-0.056610716012294`** (diff vs libxc: `-6.94e-18`) |
| **`libxc 'GGA_C_SPBE'`** | **`-0.362040491755 Eh`** | referencia (`0`) | **`-0.056610716012294`** |

**Violación del límite físico exacto $t \to \infty$ ($n = 0.1$, no polarizado):**

| $|\nabla n|$ | XCFun `'SPBEC'` baseline ($\varepsilon_c$) | `libxc GGA_C_SPBE` ($\varepsilon_c$) | XCFun `'SPBEC'` parcheado ($\varepsilon_c$) |
|---:|---:|---:|---:|
| `1.0` | `-4.5886050273e-02` | `-3.6197316246e-03` | `-3.6197316246e-03` |
| `10.0` | `-4.5752183512e-02` | `-3.9424840917e-05` | `-3.9424840917e-05` |
| `100.0` | `-4.5750820717e-02` | `-3.9460505804e-07` | `-3.9460505804e-07` |
| `10000.0` | **`-4.5750806951e-02`** (estancado, no se anula) | `-3.9460865131e-11` ($\propto t^{-2} \to 0$) | **`-3.9460865131e-11`** ($\propto t^{-2} \to 0$) |

---

### Parches mínimos compilados y verificación bit a bit (Fase 1, paso 4)

Se compilaron dos bibliotecas compartidas independientes desde `a89b783` (`build-base/libxcfun.dylib` sin modificar y `build-patched/libxcfun.dylib` con los parches mínimos [tw.patch](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/patches/tw.patch) y [spbec.patch](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/patches/spbec.patch), más el binario `testall`).  
Script: [verify_patched_and_bitid.py](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/scripts/verify_patched_and_bitid.py) | Salidas: [verify_patched_and_bitid.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/outputs/verify_patched_and_bitid.log), [testall_base.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/outputs/testall_base.log), [testall_patched.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/outputs/testall_patched.log).

1. **Parche probado en `tw.cpp`** (`0.125 * (d.gaa / d.a + d.gbb / d.b)`):
   - H $1s$ exacto: `0.499999999990097 Eh` (err `-9.90e-12`, limitado por la malla).
   - H$_2$O RKS PBE/cc-pVTZ: `57.282882450105 Eh` vs `libxc GGA_K_VW` (`diff = -2.26e-12 Eh`).
   - Punto polarizado `[e, va, vb, vgaa, vgab, vgbb]` vs `libxc GGA_K_VW`: **`max|diff| = 1.110e-15`**; consistencia FD de 4 puntos (`h = 1e-4`): **`7.797e-13`**.
2. **Parche probado en `spbec.cpp`** (`pw92eps`, `xcfun_constants::param_gamma/param_beta_gamma` y `param_gamma * phi3` dentro de `expm1`):
   - H$_2$O RKS PBE/cc-pVTZ: `-0.362040491755 Eh` vs `libxc GGA_C_SPBE` (**`diff = -1.665e-15 Eh`**).
   - Punto polarizado `[e, va, vb, vgaa, vgab, vgbb]` vs `libxc GGA_C_SPBE`: **`max|diff| = 4.163e-17`**; consistencia FD de 4 puntos (`h = 1e-4`): **`7.247e-14`**.
   - Restaura exactamente el límite $t \to \infty$ (`-3.946e-11` en $|\nabla n| = 10^4$).
3. **Regresión y verificación bit a bit (0 ULP):**
   - `xcfun_test()` y `testall`: **0 fallos** (`All tests ok`) en baseline y parcheado.
   - Comparación exacta de bytes IEEE-754 (órdenes 0, 1 y 2) sobre **todos los demás funcionales de XCFun**: los **80/80 funcionales restantes** (76 funcionales estándar + 4 dependientes de corriente `BRX`, `BRC`, `BRXC`, `CSC`, totalizando **4,492 valores `double`**) son **100% idénticos bit a bit (0 ULP, 0 bytes de diferencia)**. Los únicos funcionales que cambian son `TW` y `SPBEC`.

---

## Fase 2: comparación con `OSTRA-027` y `OSTRA-028`

Solo después de concluir y congelar la Fase 1 se leyeron `findings/OSTRA-027-xcfun-tw-squared-gradient/finding.md` (+ `report/tw.patch`) y `findings/OSTRA-028-xcfun-spbec-swapped-constants/finding.md` (+ `report/spbec.patch`).

### 1. Comparación con OSTRA-027 (`TW`, `src/functionals/tw.cpp`)

- **Lo confirmado (plena coincidencia):**
  - Ubicación exacta (`src/functionals/tw.cpp:21`), mecanismo (`pow(d.gaa + d.gbb, 2.0)` eleva al cuadrado gradientes ya al cuadrado y omite `2*d.gab`), ausencia de pruebas unitarias (`TW has no test`), valor en el átomo de H (`0.0075 Eh` vs `0.5 Eh`) y explosión por órdenes de magnitud en moléculas (`1.318e7 Eh` en H$_2$O).
  - Nivel de evidencia 3 (*defecto confirmado*) e impacto Nivel III.
- **Lo que contradice / corrige a `OSTRA-027` (hallazgo clave de la verificación):**
  - `OSTRA-027` propone en `report/tw.patch` la línea `return 1. / 8. * (d.gaa + 2.0 * d.gab + d.gbb) / d.n;` (`d.gnn / (8 * d.n)`) y concluye que iguala a `libxc GGA_K_VW` porque solo lo probó en el átomo de H ($n_\beta = 0$).
  - **Corrección:** esa expresión corresponde únicamente a la forma **sin escalamiento de espín** $T_W[n] = \int \frac{|\nabla n|^2}{8n}$, que coincide con `libxc GGA_K_VW` solo cuando $n_\alpha = n_\beta$ o $n_\beta = 0$. En cualquier sistema abierto parcialmente polarizado, **contradice a `libxc GGA_K_VW` y viola la relación exacta de escalamiento de espín de Oliver–Perdew** (el mismo defecto de espín documentado en `OSTRA-016` para `TFK`):
    - En nuestro punto polarizado ($n_\alpha=0.35, n_\beta=0.15$): `OSTRA-027` da `0.072500` frente a `0.177381` de `libxc GGA_K_VW` (error `-59.1%`).
    - En el **átomo de Li abierto (UKS PBE/cc-pVTZ)**: el parche de `OSTRA-027` da `7.132449 Eh` frente a `7.178193 Eh` de `libxc GGA_K_VW` (**diferencia de `-0.045744 Eh` = `-28.71 kcal/mol`**).
  - **Observación adicional ausente en `OSTRA-027`:** XCFun **ya incluye un segundo funcional de von Weizsäcker correcto y con escalamiento de espín**, `XC_VWK` en `src/functionals/vonw.cpp:17-23` (`d.gaa/(8*d.a) + d.gbb/(8*d.b)`), expuesto en PySCF como `'VWK'`, que ya iguala a `libxc GGA_K_VW` en H, H$_2$O y Li (`diff = -8e-12 Eh`). Por tanto, antes de enviar un PR upstream conviene decidir explícitamente si `TW` debe corregirse a la forma con escalamiento de espín (`0.125 * (d.gaa/d.a + d.gbb/d.b)`, equivalente a `VWK` y `GGA_K_VW`), o si se mantiene `0.125 * d.gnn / d.n` documentándola explícitamente como variante sin escalamiento de espín (no comparable a `GGA_K_VW` en capas abiertas).

### 2. Comparación con OSTRA-028 (`SPBEC`, `src/functionals/spbec.cpp`)

- **Lo confirmado (plena coincidencia):**
  - Ubicación exacta (`src/functionals/spbec.cpp:22-27`), ambos errores (`gamm` y `beta` intercambiados + `beta_gamma * phi3` en vez de `gamm * phi3` dentro de `expm1`), magnitud del error (~63–64.5% demasiado negativo en H$_2$O) y ausencia de referencia de pruebas (`SPBEC has no test`).
  - Nivel de evidencia 3 (*defecto confirmado*).
- **Diferencia de criterio en el parche mínimo y evidencia nueva aportada:**
  - `OSTRA-028/report/spbec.patch` corrige solo las líneas 22–27 y conserva `vwn::vwn5_eps(d)` + constantes de 6 dígitos (`0.031091`, `0.066725`), en consonancia con la convención ADF (y con `XC_VWN_PBEC` en XCFun). Eso elimina el **99.0%** del error en H$_2$O (`-144.99` de `-146.48 kcal/mol`) y deja una diferencia de `-2.378 mEh` (`-1.49 kcal/mol`) frente a `libxc GGA_C_SPBE` (que usa PW92).
  - Nuestro parche de Fase 1 ([spbec.patch](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/patches/spbec.patch)) además cambió `vwn5_eps` a `pw92eps` y usó `xcfun_constants::param_gamma/param_beta_gamma`, igualando a `libxc GGA_C_SPBE` a **`-1.67e-15 Eh`** en H$_2$O y **`4.16e-17`** en `[e, vrho, vsigma]` polarizado.
  - **Recomendación:** mantener el parche conservador de `OSTRA-028` (solo L22–27) como corrección inequívoca si XCFun desea preservar la convención VWN5 de ADF, y señalar explícitamente a los mantenedores la diferencia VWN5 vs PW92 (`-1.49 kcal/mol` en H$_2$O frente a `libxc`), del mismo modo que `pbec.cpp` distingue `PBEC` y `VWN_PBEC`.
  - **Evidencia física adicional para enriquecer `OSTRA-028`:**
    1. Demostración de que la línea 27 (`beta_gamma * phi3` en `expm1`) viola la restricción exacta $\lim_{t\to\infty}\varepsilon_c = 0$ (baseline estancado en `-4.575e-2` vs `-3.95e-11` al corregirlo).
    2. Verificación de derivadas en punto polarizado ($v_\sigma$ errado por ~3× en baseline; `max|diff| = 4.16e-17` y FD de 4 puntos `7.25e-14` tras el parche).
    3. Matiz de clasificación en la ficha: `class` indica `"New bug (parameters)"`, pero la línea 27 (`beta_gamma` en lugar de `gamm` dentro de `expm1`) es un **error de fórmula** adicional al intercambio de parámetros (`"New bug (parameters + formula)"`).

---

## Inventario de scripts, parches y salidas generadas

- **Scripts:** [verify_tw.py](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/scripts/verify_tw.py), [verify_spbec.py](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/scripts/verify_spbec.py), [verify_patched_and_bitid.py](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/scripts/verify_patched_and_bitid.py)
- **Parches probados (`git diff` sobre `a89b783`):** [tw.patch](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/patches/tw.patch), [spbec.patch](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/patches/spbec.patch)
- **Salidas crudas:** [verify_tw_baseline.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/outputs/verify_tw_baseline.log), [verify_tw_li_openshell.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/outputs/verify_tw_li_openshell.log), [verify_spbec_baseline.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/outputs/verify_spbec_baseline.log), [verify_patched_and_bitid.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/outputs/verify_patched_and_bitid.log), [testall_base.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/outputs/testall_base.log), [testall_patched.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028/outputs/testall_patched.log)
