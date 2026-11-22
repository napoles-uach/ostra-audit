# Cand-04: `libecpint` `BesselFunction::upper_bound(z, L)` viola la cota superior para todo `L >= 1`

- **Librería:** `libecpint` (`v1.1`, commit `d7966f196a5276ac0203c8ff0cca547d840ef51d`, 2026-11-18)
- **Archivos fuente:** `include/libecpint/bessel.hpp` (líneas 102–107), `src/lib/bessel.cpp` (líneas 144–152) y prueba unitaria `tests/units/bessel/bessel.cpp` (líneas 89–120)
- **Patrón del catálogo:** **Patrón 2 — Comentarios/documentación que afirman un hecho numérico no verificado** (segundo caso independiente respecto a `OSTRA-010` en `libxc`)
- **Oráculo:** Desigualdad analítica de cota superior $\text{upper\_bound}(z, L) \ge K_L(z)$ donde $K_L(z) = e^{-z} i_L(z)$
- **Clase:** Error numérico nuevo / contrato de cota violado
- **Nivel de evidencia:** **Defecto confirmado** (3,250+ violaciones reproducidas en código sin modificar; 0 violaciones tras parche mínimo en copia aislada)
- **Estado upstream:** Presente en `v1.1` / `HEAD` (`d7966f19`); no consta en `ostra-audit` ni en issues conocidos de `robashaw/libecpint`. Sin reportar (Paso 6).

---

## 1. Detección (Paso 1)

El encabezado público `include/libecpint/bessel.hpp:102-107` documenta:
```cpp
/// Calculates an upper bound to the Bessel function value at a given point
/// for a given angular momentum
double upper_bound(double z, int L) const;
```
Sin embargo, la implementación (`src/lib/bessel.cpp:145-152`) devuelve el punto tabulado **a la izquierda** en `[0, 16]`, acotado a `[minix, N]`:
```cpp
int ix = std::floor(N * z / 16.0);
int minix = L > 0 ? 1 : 0;
ix = std::min(N, std::max(minix, ix));
int lx = std::min(L, lMax);
return K[ix * K_stride + lx];
```
Esto asume implícitamente que $K_L(z) = e^{-z} i_L(z)$ es monótonamente decreciente. Ese supuesto solo es cierto para $L = 0$. Para **todo $L \ge 1$**, $K_L(0) = 0$ y $K_L(z)$ es **estrictamente creciente** hasta un máximo $z_{\max}(L) > 0$ (con $z_{\max}(1) \approx 1.61$, $z_{\max}(2) \approx 3.83$, y **$z_{\max}(L) > 16$ para todo $L \ge 6$**!). En consecuencia:
1. En todo el intervalo creciente $(0.01, z_{\max}(L))$, el punto tabulado a la izquierda `K[ix, L]` es una **cota inferior**, no una cota superior.
2. Para $L \ge 6$ y $z > 16$, el código fija `ix = N` (`z = 16.0`), subestimando $K_L(z)$ en **factores de hasta ~2,000×** para $L = 24$ (el máximo momento angular `2*LIBECPINT_MAX_L + 8 = 24` que usa la propia prueba `BesselAccuracy` de `libecpint`).

### Evidencia observada (sin parchear)

Compilando `src/lib/bessel.cpp` y `src/lib/mathutil.cpp` sin modificar con `clang++ -fsanitize=address,undefined` sobre una malla `z = 0.025, 0.075, ..., 34.975` (`evidence/unpatched.txt`):

- **`L = 0`:** 0 violaciones (único caso monótonamente decreciente).
- **`L = 1..5`:** violaciones sistemáticas en toda la rama creciente $z < z_{\max}(L)$ (p. ej. `L = 1` en `z = 0.025`: `ub < K_1(z)`; `L = 5` en `z = 2.475`: `ub/K_5 = 0.988`).
- **`L = 6..24`:** violaciones en toda la región $z \in (0.025, 35.0]$:
  - `L = 8, z = 34.975`: `ub = 3.240e-03 < K_8 = 5.058e-03` (`ratio = 0.641`)
  - `L = 12, z = 34.975`: `ub = 2.604e-04 < K_12 = 1.525e-03` (`ratio = 0.171`, subestima `5.9×`)
  - `L = 16, z = 34.975`: `ub = 9.475e-06 < K_16 = 2.974e-04` (`ratio = 0.0319`, subestima `31×`)
  - **`L = 24, z = 34.975`: `ub = 1.688e-09 < K_24 = 3.305e-06` (`ratio = 5.11e-04`, subestima `1,958×`!)**

## 2. Filtrado de artefactos (Paso 2)

- No hay errores de redondeo involucrados: la violación es estructural y alcanza casi tres órdenes de magnitud (`ratio = 5.11e-4`), excediendo por ~12 órdenes de magnitud la tolerancia de comparación `1e-12`.
- **Por qué las pruebas unitarias de `libecpint` no lo detectaron (`tests/units/bessel/bessel.cpp:89-120`):** `UpperBoundSmallZ`, `UpperBoundMidZ` y `UpperBoundBigZ` solo probaron `L = 0, 1, 2` en tres puntos: `z = SMALL/2` (donde `minix = 1` salta hacia adelante a `z = 0.01`), `z = 5.054` (mayor que $z_{\max}(2) \approx 3.83$) y `z = 17.0` (también mayor que $z_{\max}(2)$). Ningún caso de prueba evaluó un punto en la región creciente $0.01 < z < z_{\max}(L)$ ni ningún $L \ge 6$.

## 3. Cruce con documentación y código fuente (Paso 3)

- `upper_bound` es método `public` de `libecpint::BesselFunction` documentado en Doxygen (`bessel.hpp:102-107`) y cubierto por tres pruebas explícitas `EXPECT_TRUE(ub >= values[l])` en `tests/units/bessel/bessel.cpp`.
- En el estado actual del repositorio (`v1.1`), `grep -rn upper_bound` muestra que ningún núcleo interno de `libecpint` lo invoca todavía (solo la API pública y su suite de pruebas). Por ello el impacto actual en integrales internas de `libecpint` es nulo, mientras que cualquier cliente externo o algoritmo de cribado (*screening*) que use `BesselFunction::upper_bound` recibirá cotas inválidas que descartan integrales no despreciables.

## 4. Aislamiento de la causa y parche mínimo (Paso 4)

- **Líneas exactas:** `src/lib/bessel.cpp:145-152`.
- **Parche mínimo verificado en copia aislada (`repro/minimal.patch`):**
  - Para `z <= 0`: devolver `L == 0 ? 1.0 : 0.0`.
  - Para `z >= 16`: si `K[N, lx] <= K[N-1, lx]` ($K_L$ ya decrece en $z = 16$), `K[N, lx]` es cota superior válida; en caso contrario ($L \ge 6$, cuyo pico está por encima de `16`), evaluar `calculate(z, lx)`.
  - Para `0 < z < 16`: tomar `max(K[ix, lx], K[ix+1, lx])` (y evaluar `calculate(z, lx)` solo en el único subintervalo que contiene el máximo).
- **Resultado con parche (`evidence/patched.txt`):** **0 violaciones** en toda la malla para todos los momentos angulares `L = 0, 1, 2, 3, 5, 6, 8, 10, 12, 16, 24`.

## 5. Impacto en observables (Paso 5)

- **Nulo en las integrales ECP internas de `libecpint v1.1`** (no se llama internamente en `radial.cpp` ni `screening.cpp`).
- **Alto en API pública / cribado externo:** un cribado basado en `upper_bound(z, L)` para $L \ge 6$ subestima $K_L(z)$ hasta por un factor de `1,958×` (`L = 24`), pudiendo descartar falsamente bloques de integrales por encima del umbral.

## 6. Reproducibilidad

Ver `repro/repro_upper_bound.cpp`, `repro/minimal.patch`, `evidence/unpatched.txt` y `evidence/patched.txt`.
