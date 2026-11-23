# Re-check on 2026-11-23 (explorer 1), PySCF 2.15.1 wheel, Linux x86_64; dft/xcfun.py, numint.py, libxc.py and
# lib/dft/xcfun_itrf.c unchanged on master 6b288fd (2026-11-23).

## Crash location (python -X faulthandler)
xcfun.py:338 xc_type -> _itrf.XCFUN_xc_type(xid) with xid = 5000+ (placeholder) -> xcfun_enumerate_parameters(5000)
returns NULL (XCFun src/XCFunctional.cpp:302-310) -> xcfun_set(fun, NULL, 1.) -> SIGSEGV. Call stack:
rks.get_veff -> numint.nr_rks:1117 -> numint._xc_type:2815 -> xcfun.xc_type:338.
Same ids are also passed unfiltered to XCFUN_eval_xc in xcfun.eval_xc (xcfun.py:1043-1060).

## Which xc strings are affected (xcfun.parse_xc)
wb97m_v -> ((5001,1),); b97m_v -> ((5000,1),); wb97x_v -> ((5002,1),); lc_vv10 -> ((5004,1),);
scan_vv10 -> ((5007,1),); b3lyp+vv10 -> B3LYP ids + (5003,1). All 19 VV10_XC keys are mapped to 5000..5018.
Named VV10 functionals map to the placeholder ONLY: their semilocal part is not represented in the XCFun mapping.

## Supported route works (mf.nlc), H2/cc-pVDZ
b3lyp + nlc=vv10: libxc -1.16416231, xcfun -1.16416231
scan  + nlc=vv10: libxc -1.15674014, xcfun -1.15674014

## A filter-only fix would be silently wrong (repro/vv10_filter_trap.py)
scan_vv10 libxc -1.16323561 / xcfun filtered -0.52497908
lc_vv10   libxc -1.16626798 / xcfun filtered -0.52497908
b97m_v    libxc -1.16513535 / xcfun filtered -0.52497908
wb97m_v   libxc -1.15338673 / xcfun filtered -0.52497908
(-0.525 Eh = no exchange-correlation at all.) So the crash currently prevents a silent error of ~0.6 Eh; a fix
must either map the semilocal part or raise a clear error.

## Not yet done
Search of existing PySCF issues (open and closed): web search found nothing; GitHub issue search to be done by the
human reviewer.
