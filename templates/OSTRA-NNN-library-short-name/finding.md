# OSTRA-NNN: <one-line title: what fails, where, under which condition>

| Field | Value |
|---|---|
| library | <name (component, if relevant)> |
| location | <path/to/file.c (function), lines> |
| versions | <release(s) and commit hash(es) where it was observed> |
| detected | <YYYY-MM-DD> |
| detected_by | <person or agent (model, tool), batch; verified by ...> |
| oracle | <which oracle from PROTOCOL.md detected it, in words> |
| libraries | <where the defect lives, comma-separated; names from tools/graph_config.json> |
| pattern | <fragile-pattern slugs from tools/graph_config.json, comma-separated, or none> |
| oracle_class | <oracle slugs from tools/graph_config.json, comma-separated> |
| class | <class from PROTOCOL.md>, new |
| impact_level | <I / II / III / IV> (<one-line justification>) |
| status | candidate |
| report | Pending |

## Observation
What was seen, with the exact inputs and the reference used. Numbers go in a short table; raw outputs in `evidence/`.

## Artifact filter
What was done to rule out numerical artifacts (FD step scan, internal thresholds lowered, other build, etc.).

## Documentation and paper
What the documentation and the original publication (including errata) say. Design or defect?

## Cause
Exact source line(s) and the mechanism. Attribution: the minimal change tested and its effect (`repro/`).

## Impact
Effect on real observables (energies, forces, charges...). What is affected and what is not.

## Verification
Independent check by a second agent or person (`verification/`): confirmed / corrected / not reproduced.

## Corrections
Later changes to this record, dated. Never overwrite the original observation.
