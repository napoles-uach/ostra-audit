# Cand-06: Contrato implícito de `libecpint` `BesselFunction::calculate` en `z <= 0`: deja `values[1..maxL]` sin escribir y la versión escalar devuelve `1.0` para `L >= 1`

- **Librería:** `libecpint` (`v1.1`, commit `d7966f196a5276ac0203c8ff0cca547d840ef51d`, 2026-11-18)
- **Archivos fuente:** `src/lib/bessel.cpp` (líneas 156–164 y 220–225), `include/libecpint/bessel.hpp` (líneas 87–99), `src/lib/radial.cpp` (línea 252), `tests/units/bessel/bessel.cpp` (líneas 29–34)
- **Patrón del catálogo:** **Patrón 3 — Contratos de API implícitos (salidas no escritas / postcondiciones rotas)** (segundo caso independiente, análogo directo de `OSTRA-004` en `PySCF libcgto`)
- **Oráculo:** Contrato de escritura completa de buffer reutilizado + límite analítico $K_L(0) = \delta_{L,0}$ y continuidad $\lim_{z \to 0^+} K_L(z) = 0$ para $L \ge 1$
- **Clase:** Contrato de API + error numérico nuevo (en la sobrecarga escalar)
- **Nivel de evidencia:** **Defecto confirmado** (ambos fallos reproducidos en código original y corregidos con un parche de 2 líneas en copia aislada)
- **Estado upstream:** Presente en `v1.1` / `HEAD` (`d7966f19`); no consta en `ostra-audit` ni en issues conocidos. Sin reportar (Paso 6).

---

## 1. Detección (Paso 1)

`libecpint::BesselFunction` expone dos sobrecargas públicas para evaluar $K_L(z) = e^{-z} i_L(z)$:

1. **Sobrecarga vectorial `calculate(double z, int maxL, std::vector<double>& values)` (`src/lib/bessel.cpp:162-163`):**
   ```cpp
   // Set K_0(z) = 1.0, and K_l(z) = 0.0 (for l != 0) if z <= 0
   if (z <= 0) values[0] = 1.0;
   ```
   A pesar de que el comentario inmediatamente anterior declara que establece $K_l(z) = 0.0$ para $l \ne 0$, el código **solo escribe `values[0]`** y deja `values[1..maxL]` completamente sin escribir. Cuando el llamador reutiliza un vector `values` entre llamadas (práctica habitual en lazos de cuadratura), `values[1..maxL]` conserva silenciosamente los valores residuales de la evaluación anterior.

2. **Sobrecarga escalar `calculate(double z, int L)` (`src/lib/bessel.cpp:223-224`):**
   ```cpp
   if (z <= 0)
     value = 1.0;
   ```
   Devuelve **`1.0` para todo `L >= 1`** (el valor analítico exacto es $K_L(0) = 0.0$), produciendo un salto discontinuo de $O(z^L / (2L+1)!!)$ a $1.0$ en $z = 0$.

### Evidencia observada (sin parchear)

Compilando `src/lib/bessel.cpp` y `src/lib/mathutil.cpp` sin modificar (`evidence/unpatched.txt`):

- **Sobrecarga vectorial** (`bf.calculate(5.0, 6, vals)` seguido de `bf.calculate(0.0, 6, vals)`):
  - `l = 0`: `1.000000e+00` (exacto `1.0`)
  - `l = 1`: `8.000545e-02` (**dato obsoleto de `z = 5.0`**, exacto `0.0`)
  - `l = 2`: `5.199219e-02` (**dato obsoleto**, exacto `0.0`)
  - `l = 3..6`: `2.801e-02, 1.277e-02, 5.021e-03, 1.728e-03` (**todos obsoletos**, exacto `0.0`)
- **Sobrecarga escalar** (`bf.calculate(0.0, L)` y `bf.calculate(-1.0, L)`):
  - `l = 1`: `calc(0.0) = +1.000000e+00`, `calc(-1.0) = +1.000000e+00`, mientras `calc(1e-12) = +3.333333e-13` (exacto en `0`: `0.0`).
  - `l = 2..6`: `calc(0.0) = +1.000000e+00` para todos, frente a `calc(1e-12) = 6.67e-26 .. 7.40e-78` (exacto en `0`: `0.0`).

## 2. Filtrado de artefactos (Paso 2)

- Verificado con `clang++ -fsanitize=address,undefined` directamente sobre el código fuente sin modificar.
- **Por qué no lo detectó la prueba unitaria (`tests/units/bessel/bessel.cpp:29-34`):** `CalculateSmallZ` inicializa `std::vector<double> values(3, 0.0)` justo antes de llamar `bessie.calculate(-1.0, 2, values)` y solo comprueba `EXPECT_DOUBLE_EQ(values[0], 1.0)`. Nunca reutilizó un vector con datos previos ni invocó la sobrecarga escalar `calculate(0.0, L)` para `L >= 1`.

## 3. Cruce con documentación y código fuente (Paso 3)

- La documentación pública (`bessel.hpp:87-99`) no advierte en ningún lugar que el llamador deba poner `values` a cero antes de invocar `calculate`.
- Confirmación interna del carácter implícito del contrato: en `src/lib/radial.cpp:252`, un desarrollador de `libecpint` tuvo que insertar manualmente un rodeo con el comentario explícito:
  `// calculate only writes K_0 then, so the other orders are set here`
  En lugar de corregir `bessel.cpp`, el parche local en `radial.cpp` dejó expuestos tanto a la API pública como a cualquier otro llamador interno o externo.

## 4. Aislamiento de la causa y parche mínimo (Paso 4)

- **Líneas exactas:** `src/lib/bessel.cpp:163` y `src/lib/bessel.cpp:224`.
- **Parche mínimo verificado en copia aislada (`repro/minimal.patch`):**
  1. En la sobrecarga vectorial (`bessel.cpp:163`): `for (int l = 1; l <= maxL; l++) values[l] = 0.0;`.
  2. En la sobrecarga escalar (`bessel.cpp:224`): `value = (L == 0) ? 1.0 : 0.0;`.
- **Resultado con parche (`evidence/patched.txt`):**
  - Sobrecarga vectorial: `values[0] = 1.0`, `values[1..6] = 0.0` incluso en buffers reutilizados.
  - Sobrecarga escalar: `calc(0.0, 0) = 1.0`, `calc(0.0, L) = 0.0` para todo `L >= 1`, restaurando continuidad con `z -> 0+`.

## 5. Impacto en observables (Paso 5)

- En el propio `libecpint v1.1`, `radial.cpp:252` mitigó manualmente la sobrecarga vectorial para su ruta concreta; sin embargo, `BesselFunction::calculate` (ambas sobrecargas) forma parte de la API pública instalada (`libecpint/bessel.hpp`). Cualquier código consumidor que reutilice `values` o invoque `calculate(0.0, L)` con `L >= 1` (p. ej. integrales con un centro coincidente con el núcleo ECP, donde $z = 2\alpha |\mathbf{R}_A - \mathbf{R}_C| r = 0$) obtiene datos obsoletos o `1.0` en lugar de `0.0` (error absoluto `+1.0`, infinito relativo).

## 6. Reproducibilidad

Ver `repro/repro_bessel_contract.cpp`, `repro/minimal.patch`, `evidence/unpatched.txt` y `evidence/patched.txt`.
