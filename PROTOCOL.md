# Audit protocol

## Oracles
1. **Consistency** — the potential must be the derivative of the energy (v = ∂(ρε)/∂ρ, σ, τ), checked with
   4-point finite differences and a step-size scan (keep the best step; a deviation that moves with the step is an artifact).
2. **Invariances** — translation, rotation, permutation of atoms, equivalent unit cells, supercells.
3. **Exact identities and constraints** — spin-polarized evaluation at ρ↑ = ρ↓ equals unpolarized; exchange spin
   scaling; erfc + erf = full Coulomb; uniform-gas limits.
4. **Analytic references** — closed forms or arbitrary-precision evaluation (mpmath).
5. **Independent implementations** — the same functional in two libraries (e.g. XCFun vs libxc), with exact parameter attribution.
6. **API contracts** — postconditions a user would expect (outputs written, setters effective on composite objects, no crash).
7. **Branch boundaries** — value and derivatives on both sides of every hard switch or threshold found in the source.

## Six steps (no shortcut)
1. **Detect** with an oracle.
2. **Filter artifacts** — FD step, internal density/σ/τ thresholds (lower them and re-test), badly scaled metrics, zeros by symmetry.
3. **Cross-check documentation and the original paper** — documented limitation or deliberate design?
4. **Isolate the cause** — minimal reproducer and exact source line; attribution by changing only the suspected
   part (patch, parameter, post-hoc subtraction) and showing the discrepancy disappears.
5. **Measure the impact** on real observables (energies, forces, charges), not only on point values.
6. **Report** — only with explicit human approval, after checking it is not already known. See [`REPORTING.md`](REPORTING.md).

## Classes
New numerical bug · numerical precision · design (deliberate, possibly known in the literature) · documented
limitation · API contract · memory safety · convergence · rediscovery (already known/reported) · no impact.

## Impact levels
- **I** changes observables at a publishable level.
- **II** breaks numerical consistency (energy vs force, energy vs potential) without changing typical values.
- **III** niche, but severe or silent when it occurs.
- **IV** no measurable impact.

## Status workflow
`candidate → verified → attributed → impact measured → reported → resolved`, plus `reclassified` / `refuted`.
Refutations and reclassifications are kept: they are part of the measurement of the method.

## Verification by a second agent
Findings proposed by one agent are re-derived by another that did not write the original harness: it reads the
cited source, writes its own test and reports confirmed / corrected / not reproduced. See `agents/LOG.md`.
