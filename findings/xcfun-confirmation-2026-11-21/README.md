# XCFun confirmation phase (2026-11-21)

Minimal isolated fixes to XCFun master (a89b783) for OSTRA-012 (REVTPSSX), OSTRA-013 (REVTPSSC), OSTRA-014 (R4SCANX) and
OSTRA-016 (TFK), built with XCFUN_MAX_ORDER=3 and loaded into PySCF 2.15.1 with LD_PRELOAD; compared with libxc 7.1.2 on
400 random meta-GGA points (`xcfun_fix_validate.py`). Four builds in `validation_results.txt`: PySCF-bundled, vanilla
source, minimal fixes, minimal fixes + libxc constant conventions (more digits for revTPSS c and e; exact PBE beta).
Publications checked: arXiv:0906.5029 (revTPSS), PRL 106, 179902 (erratum: tables only), arXiv:2110.00647v1 (r4SCAN).
