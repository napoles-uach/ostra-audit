#!/usr/bin/env python3
"""Generates the complete, line-verified CSV inventory of numerical thresholds,
branch switches, clamps, fixed series/iteration orders, and hardcoded constants
across maple/ and src/ of libxc devel (commit 3c7794c50), merging measured
boundary jumps from results_branch_boundaries.csv and results_c_helpers_boundaries.csv.
"""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SG = json.loads((ROOT / "scripts" / "source_graph.json").read_text())
FUNCS = {int(k): v for k, v in SG["funcs"].items()}
MPL_TO_FIDS = SG["mpl_to_fids"]
C_TO_FIDS = SG["c_to_all_fids"]

PRIORITY_NAMES = [
    ("gga_x_pbe", "PBE"),
    ("gga_c_pbe", "PBE"),
    ("hyb_gga_xc_pbeh", "PBE0"),
    ("hyb_gga_xc_hse06", "HSE06"),
    ("mgga_x_scan", "SCAN"),
    ("mgga_c_scan", "SCAN"),
    ("mgga_x_r2scan", "r2SCAN"),
    ("mgga_c_r2scan", "r2SCAN"),
    ("mgga_x_r2scanl", "r2SCAN-L"),
    ("mgga_c_r2scanl", "r2SCAN-L"),
    ("mgga_x_tpss", "TPSS"),
    ("mgga_c_tpss", "TPSS"),
    ("hyb_mgga_xc_tpssh", "TPSSh"),
    ("hyb_gga_xc_b3lyp", "B3LYP"),
    ("hyb_gga_xc_wb97x", "wB97X"),
]
NAME_TO_ID = {v["name"]: fid for fid, v in FUNCS.items()}


def load_boundary_measurements():
    by_key = {}
    f_csv = ROOT / "tests" / "results_branch_boundaries.csv"
    with f_csv.open() as fh:
        for r in csv.DictReader(fh):
            key = (r["archivo_linea"], r["funcional"], r["cambio_de_rama"])
            by_key[key] = {
                "probado": "Si",
                "val_1e9": f"{float(r['zk_rel_jump_wide']):.3e}",
                "val_1e14": f"{float(r['zk_rel_jump_tight']):.3e}",
                "der_1e9": f"{float(r['der_rel_jump_wide']):.3e}",
                "der_1e14": f"{float(r['der_rel_jump_tight']):.3e}",
                "canal": r["der_canal_peor"],
            }
            # Also store by (archivo_linea, round(xc,6)) as fallback when unique
            by_key.setdefault((r["archivo_linea"], round(float(r["x_c"]), 6)), by_key[key])

    c_csv = ROOT / "tests" / "results_c_helpers_boundaries.csv"
    with c_csv.open() as fh:
        for r in csv.DictReader(fh):
            key = (r["file_line"], r["helper"])
            by_key[key] = {
                "probado": "Si (C)",
                "val_1e9": f"{float(r['val_rel_jump_1e9']):.3e}",
                "val_1e14": f"{float(r['val_rel_jump_1e14']):.3e}",
                "der_1e9": f"{float(r['der_rel_jump_1e9']):.3e}",
                "der_1e14": f"{float(r['der_rel_jump_1e14']):.3e}",
                "canal": "f'(x)",
            }
            by_key.setdefault((r["file_line"], round(float(r["xc"]), 6)), by_key[key])
    return by_key


def resolve_fids(spec):
    """Resolve list of functional IDs from a file path or explicit list of functional names."""
    if isinstance(spec, list):
        out = set()
        for item in spec:
            if item in NAME_TO_ID:
                out.add(NAME_TO_ID[item])
            elif item in MPL_TO_FIDS:
                out.update(MPL_TO_FIDS[item])
            elif item in C_TO_FIDS:
                out.update(C_TO_FIDS[item])
        return sorted(out)
    if spec == "ALL":
        return sorted(FUNCS.keys())
    if spec == "ALL_GGA":
        return sorted(fid for fid, v in FUNCS.items() if v["family"] in (2, 32))
    if spec == "ALL_MGGA":
        return sorted(fid for fid, v in FUNCS.items() if v["family"] in (4, 64))
    if spec == "ALL_FHC":
        return sorted(fid for fid, v in FUNCS.items() if v["flags"] & (1 << 17))
    if spec == "ALL_MIX_FHC":
        return sorted(fid for fid, v in FUNCS.items() if (v["flags"] & (1 << 17)) and v["aux_ids"])
    if spec in MPL_TO_FIDS:
        return MPL_TO_FIDS[spec]
    if spec in C_TO_FIDS:
        return C_TO_FIDS[spec]
    return []


def summarize_fids(fids):
    fset = set(fids)
    prios = []
    for nm, tag in PRIORITY_NAMES:
        fid = NAME_TO_ID.get(nm)
        if fid in fset and tag not in prios:
            prios.append(tag)
    ex = [FUNCS[fid]["name"] for fid in fids[:6]]
    if len(fids) > 6:
        ex_str = ", ".join(ex) + f" (+{len(fids)-6} mas)"
    else:
        ex_str = ", ".join(ex)
    return len(fids), (", ".join(prios) if prios else "Ninguno"), ex_str


def main():
    meas = load_boundary_measurements()

    # Tuple schema:
    # (prio, cat, file, line, sym, var, cut_str, fids_spec, justif, ref, meas_lookup_key, notes)
    raw = [
        # =====================================================================
        # 1. ALTA PRIORIDAD (USO REAL: PBE, PBE0, HSE06, SCAN, r2SCAN, r2SCAN-L, TPSS, B3LYP, wB97X)
        # =====================================================================
        ("1_alta_uso_real", "cambio_de_rama", "maple/scaling.mpl", 17, "s_scaling_2 (corte inferior)",
         "s", "s = 1.0", "maple/scaling.mpl", "Parcial / ambiguo (atribuido a 1a version Henderson/Gaussian; CP2K usa corte continuo s=8.3)",
         "Henderson et al., JCP 131, 044108 (2009) [texto tras pago no verificado directamente]",
         ("maple/scaling.mpl:17", "HSE06 (hyb_gga_xc_hse06)", "s_scaling_2 (corte inferior)"),
         "CANDIDATO #2a: Salto C0 persistente en s_bar (-5.14e-4) -> salto 9.84e-5 en zk y 6.80e-4 en vsigma de HSE06 (1.51e-4 y 1.53e-3 en LC-wPBE)."),

        ("1_alta_uso_real", "cambio_de_rama", "maple/scaling.mpl", 18, "s_scaling_2 (corte superior + smax)",
         "s", "s = 15.0 (smax = 8.572844)", "maple/scaling.mpl", "No (recorte numerico para evitar desborde de exp(s - smax))",
         "Henderson et al., JCP 131, 044108 (2009) [smax=8.572844 fisico; corte duro s=15 numerico]",
         ("maple/scaling.mpl:18", "HSE06 (hyb_gga_xc_hse06)", "s_scaling_2 (corte superior)"),
         "CANDIDATO #2b (NUEVO): Salto +1.62e-3 en s_bar (3x mayor que en s=1) y caida del 100% en ds_bar/ds -> salto 1.58 (158%) en vsigma de HSE06."),

        ("1_alta_uso_real", "cambio_de_rama", "src/expint_e1.c", 140, "xc_E1_scaled_jet",
         "z (aux5 en wpbeh)", "z = 18.0", ["maple/gga_exc/gga_x_wpbeh.mpl", "maple/gga_exc/gga_c_ft97.mpl", "maple/lda_exc/lda_x_1d_exponential.mpl"],
         "No (artefacto numerico de implementacion en refactor de jets 7.1+)",
         "SLATEC / serie asintotica divergente e^z E1(z)",
         ("src/expint_e1.c:140", "LC-wPBE (hyb_gga_xc_lc_wpbe)", "xc_E1_scaled_jet (z=18)"),
         "CANDIDATO #1 (REGRESION NUEVA en devel vs 7.0.0): Resto superasintotico ~e^-18 = 1.5e-8 trunca E1s'(z) con salto 1.64e-6; rompe consistencia vrho vs d(rho*eps)/drho en LC-wPBE (3.90e-6) y HSE06 (2.66e-7)."),

        ("1_alta_uso_real", "cambio_de_rama", "maple/gga_exc/gga_x_wpbeh.mpl", 66, "eb1 / term1 (wcutoff)",
         "w = omega / k_F", "w = 14.0", "maple/gga_exc/gga_x_wpbeh.mpl",
         "No (corte entre ajuste 8-gaussiano y 1-gaussiano + asimetria < vs > en w==14)",
         "Heyd, Scuseria & Ernzerhof, JCP 118, 8207 (2003); Henderson 2009",
         ("maple/gga_exc/gga_x_wpbeh.mpl:294", "LC-wPBE (hyb_gga_xc_lc_wpbe)", "term1 / eb1 (wcutoff=14)"),
         "CANDIDATO #2c: Salto de rama 8.69e-8 en zk (LC-wPBE) y 7.39e-9 (HSE06), MAS pico puntual de 0.58% (5.8e-3) exactamente en el float w == 14.0 porque eb1 usa w<14 y term1 usa w>14."),

        ("1_alta_uso_real", "cambio_de_rama", "maple/gga_exc/gga_x_wpbeh.mpl", 121, "wpbeh_EG (EGscut)",
         "s", "s = 0.08 (sombreado por wpbeh_scut = 0.13)", "maple/gga_exc/gga_x_wpbeh.mpl",
         "No (expansion Taylor orden 22 para evitar 0/0; codigo muerto bajo s=0.13)",
         "Heyd, Scuseria & Ernzerhof, JCP 118, 8207 (2003)",
         ("maple/gga_exc/gga_x_wpbeh.mpl:121", "HSE06 (hyb_gga_xc_hse06)", "wpbeh_EG (sombra bajo 0.13)"),
         "Suave a nivel ULP (3.28e-15); codigo inactivo porque term2345_smooth sustituye toda la rama para s < 0.13."),

        ("1_alta_uso_real", "cambio_de_rama", "maple/gga_exc/gga_x_wpbeh.mpl", 208, "term2345_smooth / f_wpbeh_quad",
         "s", "s = 0.13 (serie orden 23 en s)", "maple/gga_exc/gga_x_wpbeh.mpl",
         "No (estabilizacion numerica para eliminar ruido de vsigma en s pequeno)",
         "Heyd, Scuseria & Ernzerhof, JCP 118, 8207 (2003)",
         ("maple/gga_exc/gga_x_wpbeh.mpl:208", "HSE06 (hyb_gga_xc_hse06)", "term2345_smooth"),
         "Suave a precision de maquina (zk: 3.28e-15, vsigma: 6.40e-13); resuelve el ruido historico de HSE para s pequeno."),

        ("1_alta_uso_real", "cambio_de_rama", "maple/gga_exc/gga_x_wpbeh.mpl", 232, "t10 (log vs xc_log1p)",
         "aux4(w, s)", "aux4 = wpbeh_D = 0.57786348", "maple/gga_exc/gga_x_wpbeh.mpl",
         "No (identidad algebraica libre de cancelacion)",
         "Heyd, Scuseria & Ernzerhof, JCP 118, 8207 (2003)",
         ("maple/gga_exc/gga_x_wpbeh.mpl:232", "HSE06 (hyb_gga_xc_hse06)", "t10 (log vs log1p)"),
         "Continuo a precision de maquina (1.46e-13 en zk, 6.41e-13 en vsigma)."),

        ("1_alta_uso_real", "umbral_o_recorte", "maple/gga_exc/gga_x_wpbeh.mpl", 78, "aux4 / f_wpbeh0 pisos",
         "aux4, s_bar", "aux4 >= eps^(2/3) (~3.69e-11); s_bar >= 1e-15", "maple/gga_exc/gga_x_wpbeh.mpl",
         "No (evita inf - inf entre E1_scaled(aux5) y log(aux4/aux6) cuando w=s=0)",
         "Heyd, Scuseria & Ernzerhof, JCP 118, 8207 (2003)",
         None, "Preserva la cancelacion logaritmica exacta en el limite PBE0 (omega=0, s->0)."),

        ("1_alta_uso_real", "cambio_de_rama", "maple/mgga_exc/mgga_k_pc07.mpl", 27, "pc07_fab (corte inferior zlo en r2SCAN-L / SCAN-L)",
         "z_pc07", "zlo = a/40 = 0.044618 (para a=1.784720, b=0.258304)",
         ["mgga_x_r2scanl", "mgga_c_r2scanl", "mgga_x_scanl", "mgga_c_scanl", "mgga_x_revscanl", "mgga_c_revscanl", "mgga_k_pc07_opt"],
         "No (umbral 1/40 calibrado para b=3 donde da 4e-51, heredado sin escalar cuando b bajo a 0.258304)",
         "Perdew & Constantin, PRB 75, 155109 (2007); Mejia-Rodriguez & Trickey, PRB 102, 121109 (2020)",
         ("maple/mgga_exc/mgga_k_pc07.mpl:27", "r2SCAN-L_C (mgga_c_r2scanl)", "pc07_fab (zlo = a/40, b=0.258)"),
         "CANDIDATO #3 (ERROR NUMERICO NUEVO): pc07_fab0(a/40) vale 4.59e-5 (no 1.6e-17!); salto persistente 1.79e-6 en zk (r2SCAN-L_C), 2.03e-7 (r2SCAN-L_X) y 100% en vlapl."),

        ("1_alta_uso_real", "cambio_de_rama", "maple/mgga_exc/mgga_k_pc07.mpl", 27, "pc07_fab (corte superior zhi en r2SCAN-L / SCAN-L)",
         "z_pc07", "zhi = 39*a/40 = 1.740102",
         ["mgga_x_r2scanl", "mgga_c_r2scanl", "mgga_x_scanl", "mgga_c_scanl", "mgga_k_pc07", "mgga_k_pc07_opt"],
         "No (corte numerico contra 0/0 en z=a; 1 - fab0(39a/40) ~ 1.1e-18)",
         "Perdew & Constantin, PRB 75, 155109 (2007)",
         ("maple/mgga_exc/mgga_k_pc07.mpl:27", "r2SCAN-L_X (mgga_x_r2scanl)", "pc07_fab (zhi = 39*a/40)"),
         "Suave a nivel ULP (1.61e-14 en zk, 1.35e-13 en vlapl) porque exp(-40/39)/exp(-40) sigue siendo ~1e17."),

        ("1_alta_uso_real", "cambio_de_rama", "maple/tpss_c.mpl", 34, "tpss_C0 (corte ferromagnetico)",
         "1 - |zeta|", "1 - |zeta| = 1e-12", "maple/tpss_c.mpl",
         "Limitacion analitica del funcional + parche duro de libxc",
         "Tao, Perdew, Staroverov & Scuseria, PRL 91, 146401 (2003), Ec. (34)",
         ("maple/tpss_c.mpl:34", "TPSS_C (mgga_c_tpss, grad_zeta=0)", "tpss_C0 (1-|z|=1e-12, grad_z=0)"),
         "CANDIDATO #4: Ec. (34) diverge como (1-|z|)^(-1/3) en derivadas; el corte duro 1e-12 hace saltar vsigma_bb de +873.33 a +0.060 Eh (100% rel) y vrho_b en 0.025 Eh."),

        ("1_alta_uso_real", "cambio_de_rama", "maple/tpss_c.mpl", 40, "tpss_aux = m_min(tau_W/tau, 1)",
         "tau_W / tau", "tau_W / tau = 1.0 (alpha = 0)", "maple/tpss_c.mpl",
         "Si (definicion max(tau_W, tau) / m_min en regiones de 1 electron)",
         "Tao, Perdew, Staroverov & Scuseria, PRL 91, 146401 (2003), Ec. (11)-(12)",
         ("maple/tpss_c.mpl:40", "TPSS_C (mgga_c_tpss)", "tpss_aux m_min(..., 1)"),
         "Quiebre C1 de diseno cuando se desactiva FHC (zk continuo 4.01e-14; vtau salta a 0 para tau < tau_W)."),

        ("1_alta_uso_real", "umbral_o_recorte", "src/functionals.c", 511, "xc_func_set_fhc_enforcement (contrato API)",
         "XC_FLAGS_ENFORCE_FHC", "No propaga a p->func_aux[ii] (38 funcionales mixtos afectados)",
         "ALL_MIX_FHC", "No (omision de bucle sobre func_aux respecto a lineas 463-509)",
         "Lehtola & Marques, JCP 159, 114116 (2023) / Libxc 7.0 API",
         None,
         "CANDIDATO #5 (CONTRATO DE API NUEVO): En los 38 meta-GGA mixtos/hibridos con FHC (TPSSh, r2SCAN0, SCAN0, SCAN-rVV10, PW6B95...), desactivar FHC es un no-op silencioso."),

        ("1_alta_uso_real", "cambio_de_rama", "maple/mgga_exc/mgga_x_scan.mpl", 28, "scan_f_alpha_left / right / a=1",
         "alpha", "alpha = 0.98180, 1.0, 1.02165", "maple/mgga_exc/mgga_x_scan.mpl",
         "Si (a=1 en Ec. 9 de Sun 2015; cortes 0.9818 y 1.0217 donde exp < DBL_EPSILON)",
         "Sun, Ruzsinszky & Perdew, PRL 115, 036402 (2015), Ec. (9)",
         ("maple/mgga_exc/mgga_x_scan.mpl:33", "SCAN_X (mgga_x_scan)", "scan_f_alpha (a=1)"),
         "C2 continuo en el papel y suave a precision de maquina (zk: 6.37e-16, vtau: 9.37e-13 en a=1)."),

        ("1_alta_uso_real", "cambio_de_rama", "maple/mgga_exc/mgga_x_scan.mpl", 63, "scan_gx / m_one_minus_recexp",
         "s", "p^(1/4)/4.9479 = -1/ln(eps) (s = 0.01884)", "maple/mgga_exc/mgga_x_scan.mpl",
         "No (guarda contra subdesborde de exp(-1/x) < DBL_EPSILON)",
         "Sun, Ruzsinszky & Perdew, PRL 115, 036402 (2015)",
         ("maple/mgga_exc/mgga_x_scan.mpl:63", "SCAN_X (mgga_x_scan)", "scan_gx / m_one_minus_recexp"),
         "Continuo a precision de maquina (zk: 2.97e-16, vsigma: 1.34e-11)."),

        ("1_alta_uso_real", "cambio_de_rama", "maple/scan_f_alpha.mpl", 29, "r2scan_f_alpha (empalmes en a=0 y a=2.5)",
         "alpha_reg", "alpha_reg = 0.0 y alpha_reg = 2.5 (polinomio grado 7)",
         ["maple/mgga_exc/mgga_x_r2scan.mpl", "maple/mgga_exc/mgga_c_r2scan.mpl"],
         "Si (interpolacion polinomica C3 en [0, 2.5] definida en Furness et al. 2020)",
         "Furness et al., JPCL 11, 8208 (2020), Ec. (7)-(8)",
         ("maple/scan_f_alpha.mpl:29", "r2SCAN_X (mgga_x_r2scan)", "r2scan_f_alpha (a=0)"),
         "Empalme C3 exacto en a=0 (3.67e-14) y suave a ~3e-9 en vtau en a=2.5 (coeficientes publicados redondeados a ~9 cifras)."),

        ("1_alta_uso_real", "cambio_de_rama", "maple/util.mpl", 192, "x_asinh_cx (serie Taylor orden 16 vs asinh)",
         "c^2 * xs^2", "c^2 * p = 1/16 (xs = 0.25 en B88 / B3LYP)",
         ["maple/gga_exc/gga_x_b88.mpl", "maple/gga_exc/gga_x_pw91.mpl", "maple/gga_exc/gga_x_2d_b88.mpl", "maple/gga_exc/gga_x_hcth_a.mpl"],
         "No (evita 1/sqrt(p) en derivadas pares para gradiente pequeno)",
         "Becke, PRA 38, 3098 (1988)",
         ("maple/util.mpl:192", "B3LYP (hyb_gga_xc_b3lyp)", "x_asinh_cx (p = 1/16)"),
         "Grupo de control (B3LYP / PW91): continuo a nivel de ULP (zk: 5.27e-16, vsigma: 1.12e-14)."),

        ("1_alta_uso_real", "cambio_de_rama", "maple/util.mpl", 163, "m_saturate(y)",
         "y", "y = 1.0 (y/(1+y) vs 1 - 1/(1+y))",
         ["maple/gga_exc/gga_c_pbe.mpl", "maple/gga_exc/gga_x_optx.mpl", "maple/gga_exc/gga_x_s12.mpl", "maple/mgga_exc/mgga_x_pkzb.mpl"],
         "No (reescritura algebraica libre de cancelacion en cola)",
         "Perdew, Burke & Ernzerhof, PRL 77, 3865 (1996)",
         ("maple/gga_exc/gga_c_pbe.mpl:73", "M08-HX_C (mgga_c_m08 via gga_c_pbe)", "m_saturate(A*f1)"),
         "Grupo de control (PBE / PBE0 / OPTX): identidad algebraica exacta a nivel de ULP (4.01e-14 en zk, 2.76e-13 en vrho)."),

        ("1_alta_uso_real", "cambio_de_rama", "maple/attenuation.mpl", 51, "attenuation_erf / one_minus_attenuation_erf",
         "a = omega / (2*k_F)", "a = 1/16 (SR ord 12), 1/2 (mid), 2^24 (LR ord 4)",
         "maple/attenuation.mpl",
         "No (empalmes enforce_smooth entre serie de corto alcance, forma exacta, cola xc_erf_tail y serie asintotica)",
         "Gill, Adamson & Pople, Mol. Phys. 88, 1005 (1996); Chai & Head-Gordon, JCP 128, 084106 (2008)",
         ("maple/enforce_smooth_cutoffs.mpl:21", "wB97X (hyb_gga_xc_wb97x)", "attenuation_erf (mid_cutoff=1/2)"),
         "Grupo de control (wB97X / CAM-B3LYP / LRC): continuo hasta ~1e-12..1e-13 en los tres cortes."),

        ("1_alta_uso_real", "terminos_o_iteraciones_fijos", "maple/b97.mpl", 18, "b97_g / b97mv_g (series de potencia u, w)",
         "i", "n = 5 terminos (wB97X: c_0..c_4); B97M-V/wB97M-V: combinaciones fijas",
         ["maple/b97.mpl", "maple/b97mv.mpl"],
         "Si (definicion variacional del ansatz de Becke 1997 / Mardirossian & Head-Gordon)",
         "Becke, JCP 107, 8554 (1997); Chai & Head-Gordon, JCP 128, 084106 (2008)",
         None, "Grado polinomico fijo definido en el articulo original."),

        # =====================================================================
        # 2. INFRAESTRUCTURA COMPARTIDA (src/ y maple/util.mpl, attenuation.mpl)
        # =====================================================================
        ("2_infraestructura_compartida", "umbral_o_recorte", "src/functionals.c", 360, "dens_threshold / sigma_threshold / zeta_threshold / tau_threshold",
         "rho, sigma, zeta, tau", "dens_thr = 1e-15 (GGA/MGGA std), 1e-12 (PBE_C y 54 mas), 1e-10..1e-6 (M06/HJS); sigma_thr = dens_thr^(4/3); zeta_thr = eps; tau_thr = 1e-20",
         "ALL", "No (umbrales de cribado numerico de libxc)",
         "Marques, Oliveira & Burnus, CPC 183, 2272 (2012); Lehtola et al., SoftwareX 7, 1 (2018)",
         None, "En work_{lda,gga,mgga}_inc.c y screen_dens en util.mpl; causan discrepancias aparentes pol/unpol si no se escalan por canal."),

        ("2_infraestructura_compartida", "umbral_o_recorte", "src/work_gga_inc.c", 87, "Recorte Cauchy-Schwarz de sigma_ab",
         "sigma[1]", "|sigma_ab| <= sqrt(sigma_aa*sigma_bb)*(1 - 8*eps), piso s_min = 16*eps*s_ave",
         ["ALL_GGA", "ALL_MGGA"], "Si (desigualdad fisica |grad n_a . grad n_b| <= |grad n_a||grad n_b| + margen 8*eps)",
         "Lehtola & Marques, JCP 159, 114116 (2023)",
         None, "Aplicado identicamente en work_gga_inc.c:87-91 y work_mgga_inc.c:95-99."),

        ("2_infraestructura_compartida", "umbral_o_recorte", "src/work_mgga_inc.c", 62, "Recorte de curvatura del hueco de Fermi (FHC)",
         "sigma_ss", "sigma_ss <= 8 * rho_s * tau_s (cuando XC_FLAGS_ENFORCE_FHC activo)",
         "ALL_FHC", "Si (cota fisica tau >= tau_W; opcional en tiempo de ejecucion para PAW)",
         "Lehtola & Marques, JCP 159, 114116 (2023)",
         None, "Controlado por xc_func_set_fhc_enforcement (ver Candidato #5 para funcionales mixtos)."),

        ("2_infraestructura_compartida", "umbral_o_recorte", "maple/util.mpl", 652, "xc_reduced_floor y mgga_alpha_safe_p",
         "xt, xs, (8t - p/X2S^2)/(8t)", "xc_reduced_floor = eps^(2/3) (~3.69e-11)",
         "maple/util.mpl", "No (evita 0/0 y singularidades de derivadas en gradiente nulo o tau=tau_W)",
         "Infraestructura interna de libxc", None, "Usado transversalmente en util.mpl, wpbeh, scan, gaploc, q2d, pjs18, tb09."),

        ("2_infraestructura_compartida", "cambio_de_rama", "src/util_inline.h", 15, "xc_x_minus_log1p / xc_atan_tail3 / xc_log1p_tail3 / xc_exp_tail / xc_erf_tail",
         "x", "|x|=0.25 (x-log1p), |x|=0.6 (atan_tail3), |x|=0.5 (log1p_tail3), |x|=4 (exp_tail), |z|=2 (erf_tail)",
         "maple/util.mpl", "No (colas de Taylor con criterios de parada |term| <= eps*|sum|/2, max 100-200 terminos)",
         "Infraestructura C de cancelacion nula en libxc",
         ("src/util_inline.h:15", "xc_x_minus_log1p(+0.25)"),
         "Todos los 7 cortes verificados en test_c_helpers_boundaries.py: saltos 1e-13..1e-15 en valor y derivada."),

        ("2_infraestructura_compartida", "cambio_de_rama", "src/faddeeva.c", 518, "xc_erfcx y xc_erfcx_jet",
         "x", "x = -6.1, 0.0, 4.0 (100 bins en [0,4]), 8.0 (jet), 50.0, 5e7",
         "src/faddeeva.c", "No (algoritmo de Johnson / Faddeeva + serie asintotica donde e^(-x^2) < 10^-28)",
         "Johnson, Faddeeva package (MIT)",
         ("src/faddeeva.c:518", "xc_erfcx_jet(8.0)"),
         "A diferencia de xc_E1_scaled_jet(18), en x=8 el resto no perturbativo es e^(-64) ~ 1.6e-28 << eps: salto medido 3.27e-13."),

        ("2_infraestructura_compartida", "cambio_de_rama", "src/expint_e1.c", 100, "xc_expint_e1_impl (ajustes Chebyshev SLATEC)",
         "x", "x = -10, -4, -1, +1, +4 (grados Chebyshev 39, 25, 19, 16, 25, 26)",
         "src/expint_e1.c", "No (aproximacion por tramos Chebyshev de Fullerton / SLATEC)",
         "Fullerton, SLATEC library",
         ("src/expint_e1.c:116", "xc_expint_e1_scaled(+4)"),
         "La funcion valor xc_E1_scaled(x) es suave a 1e-15 en todos sus cortes; solo el jet xc_E1_scaled_jet en z=18 falla."),

        ("2_infraestructura_compartida", "cambio_de_rama", "src/bessel.c", 54, "xc_bessel_{I0,I1,K0,K1}[_scaled]",
         "x", "x = 2*sqrt(eps), 2*sqrt(2*eps), 2.0, 3.0, 8.0 (Chebyshev grados 11-22)",
         "src/bessel.c", "No (rutinas SLATEC de Fullerton)",
         "Fullerton, SLATEC library",
         ("src/bessel.c:80", "xc_bessel_I0(3.0)"),
         "Todos los 16 cortes de Bessel continuos a ~1e-14..1e-15 en valor (en 2*sqrt(eps) xc_bessel_I0 fija derivada nula con salto O(sqrt(eps))=9.9e-8)."),

        ("2_infraestructura_compartida", "cambio_de_rama", "src/special_functions.c", 45, "LambertW y xc_lambertw_jet (iteracion de Halley)",
         "z", "|z| = eps^(1/3), z = -0.31408624, +1.14987649; max 20 iteraciones Halley, tol 2*eps",
         "src/special_functions.c", "No (inicializador por tramos + refinamiento de Halley)",
         "Corless et al., Adv. Comput. Math. 5, 329 (1996)",
         ("src/special_functions.c:49", "LambertW(-0.31408624)"),
         "Continuo a precision de maquina (1.38e-13 en valor, 2.34e-13 en derivada)."),

        ("2_infraestructura_compartida", "cambio_de_rama", "src/special_functions.c", 193, "xc_dilogarithm (serie Chebyshev grado 38)",
         "x", "x = -1, 0, 0.5, 1.0, 2.0, +-2/eps",
         ["gga_x_gg99"], "No (transformaciones de Euler de Li2 + Chebyshev 38 terminos)",
         "SLATEC / Morris, Math. Comp. (dspenc)",
         ("src/special_functions.c:220", "xc_dilogarithm(-1.0)"),
         "Continuo a 1.69e-14 en (-inf, 1); en x=1 Li2 tiene punto de ramificacion logaritmico (gga_x_gg99 solo evalua en (-1, 0))."),

        ("2_infraestructura_compartida", "terminos_o_iteraciones_fijos", "src/integrate.c", 53, "xc_integrate (cuadratura tanh-sinh)",
         "level", "XC_TANH_SINH_MAX_LEVELS = 12, tol = 2*XC_EPSILON, parada nivel >= 3",
         ["lda_x_1d_exponential", "lda_x_1d_soft", "hyb_gga_xc_edf2"],
         "No (cuadratura adaptativa numerica)",
         "Takahasi & Mori, Publ. RIMS 9, 721 (1974); corregido en commit a19650f34 (#646)",
         None, "Ya auditado y corregido upstream en #646."),

        ("2_infraestructura_compartida", "terminos_o_iteraciones_fijos", "src/mgga_x_br89_get_x_inline.h", 31, "xc_mgga_x_br89_get_x / mbrxc_get_x / math_brent",
         "rhs, iter", "Corte serie rhs = 2*eps; Brent tol = 2*eps, MAX_ITER = 500 (br89) / 100 (mbrxc)",
         ["maple/mgga_exc/mgga_x_br89.mpl", "maple/mgga_exc/mgga_x_mbrxc_bg.mpl"],
         "No (inversion numerica de la ecuacion no lineal de Becke-Roussel)",
         "Becke & Roussel, PRA 39, 3761 (1989)",
         ("src/mgga_x_br89_get_x_inline.h:32", "br89_get_x(+2*eps)"),
         "El inversor C es exacto a 2*eps (salto 0.00e+00), pero ver Candidato #8 sobre el recorte br89_cQ en Maple."),

        ("2_infraestructura_compartida", "cambio_de_rama", "maple/util.mpl", 41, "e1_scaled_G(y)",
         "y", "y = E1_SCALED_G_CROSSOVER = 100 (serie asintotica 16 terminos)",
         ["gga_c_ft97"], "No (evita cancelacion 1 - y*E1_scaled(y) para y grande; 16!/100^16 ~ 2e-19)",
         "Filatov & Thiel, IJQC 62, 603 (1997)",
         ("maple/util.mpl:41", "FT97_C (gga_c_ft97)", "e1_scaled_G (mu=100)"),
         "Continuo a precision de maquina (zk: 2.00e-13, vsigma: 8.01e-13)."),

        ("2_infraestructura_compartida", "cambio_de_rama", "maple/util.mpl", 134, "m_recexp(x) / m_one_minus_recexp(x)",
         "x", "x = -1/ln(XC_EPSILON) (~0.02774)",
         ["maple/mgga_exc/mgga_x_scan.mpl", "maple/mgga_exc/mgga_x_lak.mpl", "maple/mgga_exc/mgga_x_task.mpl", "maple/mgga_exc/hyb_mgga_xc_coach.mpl"],
         "No (corte exacto donde exp(-1/x) = DBL_EPSILON)",
         "Infraestructura interna de libxc",
         ("maple/mgga_exc/mgga_x_scan.mpl:63", "SCAN_X (mgga_x_scan)", "scan_gx / m_one_minus_recexp"),
         "Continuo a precision de maquina."),

        ("2_infraestructura_compartida", "cambio_de_rama", "maple/util.mpl", 175, "m_saturate2(y)",
         "y", "y = 1.0", ["gga_x_rge2"],
         "No (identidad algebraica libre de cancelacion)", "Ruzsinszky, Csonka & Scuseria, JCTC 5, 763 (2009)",
         ("maple/util.mpl:175", "RGE2_X (gga_x_rge2)", "m_saturate2 (y=1)"),
         "Continuo a nivel de ULP (zk: 2.72e-14, vsigma: 2.09e-13)."),

        ("2_infraestructura_compartida", "cambio_de_rama", "maple/util.mpl", 203, "m_expm1_over(y)",
         "y", "y = 1/16 (serie Taylor orden 13 vs -expm1(-y)/y)",
         ["gga_x_vmt_pbe", "gga_x_vmt_ge", "gga_x_vmt84_pbe", "gga_x_vmt84_ge", "gga_k_vt84f", "hyb_gga_x_pbe_erf_laikov"],
         "No (elimina 0/0 de -expm1(-y)/y en y -> 0)",
         "Vela, Medel & Trickey, JCP 130, 244103 (2009); Laikov, JCP 151, 094106 (2019)",
         ("maple/util.mpl:203", "Laikov_SR_X (hyb_gga_x_pbe_erf_laikov)", "m_expm1_over (gamma*s^2=1/16)"),
         "Continuo a nivel de ULP (zk: 1.37e-14, vsigma: 3.73e-14)."),

        ("2_infraestructura_compartida", "cambio_de_rama", "maple/util.mpl", 219, "m_bound(y)",
         "y", "|y| = 1.0 (y/sqrt(1+y^2) vs sgn(y)*sqrt(1 - 1/(1+y^2)))",
         ["hyb_mgga_xc_b98"], "No (evita cancelacion 1e26x en v2sigma2 en la cola)",
         "Becke, JCP 108, 2092 (1998)",
         ("maple/util.mpl:219", "B98 (hyb_mgga_xc_b98)", "m_bound (gamma_ss*q = +1)"),
         "Continuo a nivel de ULP (zk: 8.62e-15, vsigma: 7.67e-14)."),

        ("2_infraestructura_compartida", "cambio_de_rama", "maple/util.mpl", 636, "mgga_series_w(a, n, t)",
         "1 + w = 2*K/(K + t)", "1 + w = 1/5 (re-expansion binomial alrededor de w = -1)",
         ["mgga_x_m06_l", "hyb_mgga_x_m05", "mgga_x_m08_hx", "mgga_c_m08_hx", "hyb_mgga_x_dldf", "mgga_x_tau_hcth", "mgga_x_br89"],
         "No (evita cancelacion catastrofica de sum a_i (-1)^i cuando tau >> K)",
         "Zhao & Truhlar, JCP 125, 194101 (2006) (serie M05/M06/M08)",
         ("maple/util.mpl:636", "M06-L_X (mgga_x_m06_l)", "mgga_series_w (1+w = 1/5)"),
         "Identidad polinomial exacta (zk: 4.84e-13, vtau: 6.14e-12)."),

        ("2_infraestructura_compartida", "cambio_de_rama", "maple/attenuation.mpl", 67, "attenuation_erf_f2 y attenuation_erf_f3",
         "a", "a = 1/16 (SR), 1/4 (cuadratura Gauss 32 pts -> cola 16 pts), 2^24 / 2^12 (LR)",
         ["hyb_mgga_x_js18", "hyb_mgga_x_pjs18"],
         "No (cuadratura y series de alta precision generadas por check_attenuation.mpl)",
         "Patra, Jana & Samal, JCP 148, 134117 (2018)",
         ("maple/attenuation.mpl:81", "PJS18_X (hyb_mgga_x_pjs18)", "attenuation_erf_f2/f3 (mid=1/4)"),
         "Continuo a ~1.7e-12 en valor y derivada."),

        ("2_infraestructura_compartida", "cambio_de_rama", "maple/attenuation.mpl", 94, "attenuation_gau (sin consumidores activos)",
         "a", "a = 1/16, 31/64, 2^24", [],
         "No (infraestructura preparada para atenuacion gaussiana erf_gau; actualmente sin consumidores en src/)",
         "Song et al., JCP 127, 154109 (2007)", None, "Definida en attenuation.mpl y enforce_smooth_cutoffs.mpl pero ningun funcional actual la invoca."),

        ("2_infraestructura_compartida", "cambio_de_rama", "maple/attenuation.mpl", 134, "attenuation_yukawa (y piso ATT_YUK_SMALL_FLOOR)",
         "a", "a = 1/16 (SR ord 7, piso 1e-100), 1/2 (mid), 2^27 (LR ord 4)",
         ["lda_x_yukawa", "gga_x_sfat", "gga_x_sfat_pbe", "gga_x_ityh_opx", "hyb_gga_xc_camy_b3lyp", "hyb_gga_xc_camy_blyp", "hyb_gga_xc_lcy_blyp", "hyb_gga_xc_lcy_pbe"],
         "No (evita 0*inf en a->0 y cancelacion 1 - O(1/a^2) en a->inf)",
         "Savin, Recent Dev. Mod. DFT (1996); Akinaga & Ten-no, CPL 462, 348 (2008)",
         ("maple/attenuation.mpl:134", "Yukawa_LDA_X (lda_x_yukawa)", "attenuation_yukawa (mid=1/2)"),
         "Continuo a ~8e-13 en valor y derivada."),

        # =====================================================================
        # 3. FUNCIONALES SECUNDARIOS / ESPECIFICOS (LDA, GGA, MGGA)
        # =====================================================================
        ("3_secundaria", "cambio_de_rama", "maple/mgga_exc/mgga_c_lak.mpl", 61, "lak_fc y lak_at (extension a alpha < 0)",
         "alpha", "alpha = 0.0 (con m_max(rs*a, XC_EPSILON) en linea 64)",
         ["mgga_c_lak"],
         "No (error algebraico en linea 64: m_max(rs*a, eps) en vez de m_max(rs*m_abs(a), eps))",
         "Lebeda, Aschebrock & Kummel, PRL 133, 136402 (2024), SI Ec. (19)-(20)",
         ("maple/mgga_exc/mgga_c_lak.mpl:61", "LAK_C (mgga_c_lak)", "lak_fc (alpha=0)"),
         "CANDIDATO #6 (ERROR NUMERICO NUEVO cuando FHC desactivado): Para todo alpha < 0, m_max(rs*a, eps) evalua a eps, congelando lak_fc ~ 1 y causando salto 100% en vtau en alpha=0."),

        ("3_secundaria", "cambio_de_rama", "maple/mgga_exc/mgga_c_lak.mpl", 59, "lak_fc_piecewise y lak_H1",
         "x, lak_beta_t(rs)", "x = 1/eps (serie en inf orden 6); lak_beta_t(rs) = 0 (rs = 55.17)",
         ["mgga_c_lak"],
         "Si para lak_beta_t=0 (SI Ec. 25); numerico para x=1/eps",
         "Lebeda, Aschebrock & Kummel, PRL 133, 136402 (2024)",
         ("maple/mgga_exc/mgga_c_lak.mpl:79", "LAK_C (mgga_c_lak)", "lak_H1 (lak_beta_t(rs)=0)"),
         "Empalme lak_H1p / lak_H1m en rs = 55.17 es suave a 1.62e-13 en zk y 8.69e-13 en vsigma."),

        ("3_secundaria", "cambio_de_rama", "maple/mgga_exc/mgga_x_lak.mpl", 46, "lak_fx (serie en alpha = 0)",
         "alpha", "alpha = XC_EPSILON (serie orden 6)",
         ["mgga_x_lak"], "No (evita 0/0 en alpha = 0)",
         "Lebeda, Aschebrock & Kummel, PRL 133, 136402 (2024)",
         ("maple/mgga_exc/mgga_x_lak.mpl:46", "LAK_X (mgga_x_lak)", "lak_fx (alpha = eps)"),
         "Continuo a precision de maquina (zk: 0.00e+00, vtau: 4.91e-16)."),

        ("3_secundaria", "cambio_de_rama", "maple/lda_exc/lda_c_pk09.mpl", 25, "f_r (criba de polo DETOL = 1e-12 y desborde en alpha_z)",
         "rs", "|f_r_den(rs)| <= 1e-12 alrededor del polo rs = 22.73941213 (y pisos K_TOL = 1e-11)",
         ["lda_c_pk09"],
         "Parcial (polo real esta en Ec. 24 de PK09; criba DETOL=1e-12 es insuficiente para zeta != 0)",
         "Proynov & Kong, JCTC 5, 430 (2009); Erratum JCTC 17, 6632 (2021)",
         ("maple/lda_exc/lda_c_pk09.mpl:25", "PK09 (lda_c_pk09, polo f_r)", "f_r (|den|<=1e-12 cerca rs=22.74)"),
         "CANDIDATO #7 (ERROR / INESTABILIDAD GRAVE): Fuera de |den|<=1e-12, f_r crece hasta 7e14 y (1+-z)^ss desborda a inf, devolviendo NaN en zk y vrho para |rs - 22.739412| < ~3e-5 y |vrho| > 1e7 Eh en |drs| = 1e-4."),

        ("3_secundaria", "cambio_de_rama", "maple/gga_vxc/gga_x_lb.mpl", 31, "lb_f0 (corte asintotico p < 300^2)",
         "xs = |grad rho_s|/rho_s^(4/3)", "xs = 300.0 (p = 90000)",
         ["gga_x_lb", "gga_x_lbm"],
         "No (omite el 1 + del denominador para xs >= 300 cuando 3*beta*300*asinh(300) ~ 305)",
         "van Leeuwen & Baerends, PRA 49, 2421 (1994)",
         ("maple/gga_vxc/gga_x_lb.mpl:31", "LB94 (gga_x_lb)", "lb_f0 (asintota p < 300^2)"),
         "CANDIDATO #8: Corte demasiado temprano; introduce un salto persistente de 3.21e-3 (0.32%) en el potencial vrho en xs = 300."),

        ("3_secundaria", "cambio_de_rama", "maple/mgga_exc/mgga_x_br89.mpl", 36, "br89_cQ y br89_mx (banda muerta |Q| < 5e-13)",
         "Q", "|Q| = 5.0e-13 y salto de signo en Q = 0",
         ["mgga_x_br89", "mgga_x_br89_1", "mgga_x_b00", "mgga_x_br89_explicit", "mgga_x_br89_explicit_1", "mgga_c_b94", "hyb_mgga_xc_b94_hyb", "mgga_x_tb09", "mgga_x_rpp09"],
         "No (evita Q = 0 en d(br89_x)/dQ donde Q^2 aparece en denominador, aunque x(Q) es analitica en Q = 0)",
         "Becke & Roussel, PRA 39, 3761 (1989); Proynov, Gan & Kong, CPL 455, 103 (2008)",
         ("maple/mgga_exc/mgga_x_br89.mpl:36", "mgga_x_br89", "br89_cQ (banda |Q|=5e-13)"),
         "CANDIDATO #9: En |Q| < 5e-13 congela Q = +-5e-13 (vtau=vlapl de Q se anulan en banda muerta, salto relativo 1.77e-3..4.41e-3 en |Q|=5e-13; en br89_explicit salto 3.10e-7 en Q=0 por asimetria Pade)."),

        ("3_secundaria", "cambio_de_rama", "maple/lda_exc/lda_c_pz.mpl", 35, "ec (empalme rs = 1 de Perdew-Zunger 1981)",
         "rs", "rs = 1.0",
         ["lda_c_pz", "lda_c_pz_mod", "lda_c_ob_pz"],
         "Si (constantes redondeadas a 4 cifras en Tabla V de PZ81; lda_c_pz_mod usa empalme continuo)",
         "Perdew & Zunger, PRB 23, 5048 (1981)",
         ("maple/lda_exc/lda_c_pz.mpl:35", "PZ81 (lda_c_pz)", "ec (rs = 1)"),
         "Limitacion historica conocida del funcional original PZ81 (zk: 5.38e-4, vrho: 4.16e-4; corregido en lda_c_pz_mod)."),

        ("3_secundaria", "cambio_de_rama", "maple/lda_exc/lda_c_gk72.mpl", 22, "f (empalmes rs = 0.7 y rs = 10.0 de Gordon-Kim 1972)",
         "rs", "rs = 0.7 y rs = 10.0",
         ["lda_c_gk72"],
         "Si (interpolacion por tramos publicada en Ec. 20-23 de Gordon & Kim 1972)",
         "Gordon & Kim, JCP 56, 3122 (1972)",
         ("maple/lda_exc/lda_c_gk72.mpl:22", "GK72 (lda_c_gk72)", "f (rs = 10.0)"),
         "Discontinuidad propia de la definicion original de 1972 (zk: 4.36e-4, vrho: 7.50e-2 en rs=10)."),

        ("3_secundaria", "cambio_de_rama", "maple/lda_exc/lda_c_hl.mpl", 46, "hl_f0 (serie de 24 terminos para rs/r > 4)",
         "rs", "rs / r = 4 (rs = 84 en HL)",
         ["lda_c_hl", "lda_c_gl", "lda_c_vbh"],
         "No (evita cancelacion de log(1 + r/rs) en baja densidad; 24 terminos)",
         "Hedin & Lundqvist, J. Phys. C 4, 2064 (1971)",
         ("maple/lda_exc/lda_c_hl.mpl:46", "HL (lda_c_hl)", "hl_f0 (rs/r = 4)"),
         "Continuo a precision de maquina (zk: 1.73e-13, vrho: 1.71e-13)."),

        ("3_secundaria", "cambio_de_rama", "maple/lda_exc/lda_c_pw_erf.mpl", 83, "lda_c_pw_erf_g0c",
         "rs", "0.7524 * rs = 1 (rs = 1.32908)",
         ["lda_c_pmgb06"], "No (reescritura estable del exponente)",
         "Paziani et al., PRB 73, 155111 (2006)",
         ("maple/lda_exc/lda_c_pw_erf.mpl:83", "PMGB06 (lda_c_pmgb06)", "lda_c_pw_erf_g0c (0.7524*rs=1)"),
         "Continuo a precision de maquina (zk: 8.95e-14, vrho: 1.05e-13)."),

        ("3_secundaria", "cambio_de_rama", "maple/lda_exc/lda_c_2d_amgb.mpl", 48, "amgb_ex6_residual (serie orden 24)",
         "zeta", "|zeta| = 0.5",
         ["lda_c_2d_amgb"], "No (evita cancelacion O(z^6) para polarizacion pequena)",
         "Attaccalite et al., PRL 88, 256601 (2002)",
         ("maple/lda_exc/lda_c_2d_amgb.mpl:48", "2D-AMGB (lda_c_2d_amgb)", "amgb_ex6_residual (|z|=0.5)"),
         "Continuo a precision de maquina (zk: 2.58e-13, vrho: 4.22e-12)."),

        ("3_secundaria", "cambio_de_rama", "maple/lda_exc/lda_x_rel.mpl", 34, "lda_x_rel_bracket (Taylor orden 58)",
         "beta = k_F / c", "beta = 0.5 (rs = 0.02801)",
         ["lda_x_rel"], "No (evita cancelacion b*sqrt(1+b^2) - asinh(b) ~ 2*b^3/3)",
         "Rajagopal, J. Phys. C 11, L943 (1978); MacDonald & Vosko, J. Phys. C 12, 2977 (1979)",
         ("maple/lda_exc/lda_x_rel.mpl:34", "LDA_X_REL (lda_x_rel)", "lda_x_rel_bracket (b=0.5)"),
         "Continuo a ~5e-12 en valor y derivada."),

        ("3_secundaria", "cambio_de_rama", "maple/gga_exc/gga_c_revtca.mpl", 24, "msinc (Taylor orden 7)",
         "x = aa*s/rs", "x = eps^(1/4) (~1.221e-4)",
         ["gga_c_revtca"], "No (evita 0/0 en sin(x)/x)",
         "Tognetti, Cortona & Adamo, CPL 460, 536 (2008)",
         ("maple/gga_exc/gga_c_revtca.mpl:24", "revTCA (gga_c_revtca)", "msinc (x=eps^(1/4))"),
         "Continuo en zk (0.00e+00) y ~8e-9 en vsigma_bb."),

        ("3_secundaria", "cambio_de_rama", "maple/gga_exc/gga_c_w94.mpl", 30, "w94_log_az (log1p vs log)",
         "1 - |zeta|", "1 - |zeta| = 0.5",
         ["gga_c_w94"], "No (identidad algebraica log(|z|) = log1p(-(1-|z|)))",
         "Wilson, J. Mol. Struct. (Theochem) 313, 359 (1994)",
         ("maple/gga_exc/gga_c_w94.mpl:30", "W94 (gga_c_w94)", "w94_log_az (1-|z|=0.5)"),
         "Continuo a precision de maquina (1.53e-13)."),

        ("3_secundaria", "cambio_de_rama", "maple/gga_exc/gga_c_zpbeint.mpl", 55, "f (zpbeint_u = -0.5)",
         "u = alpha*t^3*ln(phi)", "u = -0.5",
         ["gga_c_zpbeint", "gga_c_zpbesol"], "No (reescritura estable de potencias de phi)",
         "Constantin, Fabiano & Della Sala, PRB 84, 233103 (2011)",
         ("maple/gga_exc/gga_c_zpbeint.mpl:56", "zPBEint (gga_c_zpbeint)", "f (zpbeint_u = -0.5)"),
         "Continuo a precision de maquina (1.03e-13 en zk, 7.87e-13 en vrho)."),

        ("3_secundaria", "cambio_de_rama", "maple/gga_exc/gga_k_lkt.mpl", 53, "lkt_sech (Taylor orden 20)",
         "p = s^2", "a * s = 1/4 (s = 0.1923)",
         ["gga_k_lkt"], "No (evita 1/s en derivadas pares de sech(a*s))",
         "Luo, Karasiev & Trickey, PRB 98, 041111 (2018)",
         ("maple/gga_exc/gga_k_lkt.mpl:53", "LKT (gga_k_lkt)", "lkt_sech (a*s = 1/4)"),
         "Continuo a nivel de ULP (6.46e-14)."),

        ("3_secundaria", "cambio_de_rama", "maple/gga_exc/gga_k_meyer.mpl", 43, "meyer_one_minus_feta_over_w / meyer_lambda / meyer_f",
         "w = (s/6)^2", "w = 1/16 (s=1.5), w = 1 (s=6.0), v = 1/w = 1e-3 (s=189.7)",
         ["gga_k_meyer"], "w=1 es singularidad logaritmica analitica de Lindhard (q=2*k_F); w=1/16 y v=1e-3 son numericos",
         "Meyer, Wang & Young, Z. Naturforsch. A 31, 898 (1976)",
         ("maple/gga_exc/gga_k_meyer.mpl:66", "Meyer (gga_k_meyer)", "meyer_lambda (w=1)"),
         "Descartado como artefacto/diseno en w=1: d(lambda)/dw diverge como ln|1-sqrt(w)| (anomalia de Kohn); w=1/16 y v=1e-3 suaves a 1e-13."),

        ("3_secundaria", "cambio_de_rama", "maple/gga_exc/gga_x_bpccac.mpl", 40, "bpccac_f / bpccac_fab_c",
         "xs", "xs = 19.0 (transicion Fermi-Dirac)",
         ["gga_x_bpccac"], "No (corte numerico donde exp(-c*(xs-19)) se satura)",
         "Bremond et al., JCP 140, 244109 (2014)",
         ("maple/gga_exc/gga_x_bpccac.mpl:40", "BPCCAC (gga_x_bpccac)", "bpccac_f (xs=19)"),
         "Continuo a precision de maquina (zk: 0.00e+00, vsigma: 9.72e-13)."),

        ("3_secundaria", "cambio_de_rama", "maple/gga_exc/gga_x_chachiyo.mpl", 54, "chachiyo_R (Taylor orden 16)",
         "q = cx * xs", "q = 1/16 (xs = 0.3489)",
         ["gga_x_chachiyo", "hyb_gga_xc_apbe0"], "No (evita cancelacion en log(x+1) para gradiente pequeno)",
         "Chachiyo & Chachiyo, Comput. Theor. Chem. 1172, 112669 (2020)",
         ("maple/gga_exc/gga_x_chachiyo.mpl:54", "Chachiyo_X (gga_x_chachiyo)", "chachiyo_R (q=1/16)"),
         "Continuo a nivel de ULP (6.47e-16)."),

        ("3_secundaria", "cambio_de_rama", "maple/gga_exc/gga_x_fd_lb94.mpl", 64, "fd_f (serie orden 24)",
         "X = xs / 2^(1/3)", "X = 1/4 (xs = 0.3150)",
         ["gga_x_fd_lb94", "gga_x_fd_revlb94", "hyb_gga_xc_fd_lb94"],
         "No (resuelve integral eliptica / logaritmica en X -> 0)",
         "Gaiduk, Firaha & Staroverov, PRL 108, 253005 (2012)",
         ("maple/gga_exc/gga_x_fd_lb94.mpl:64", "FD-LB94 (gga_x_fd_lb94)", "fd_f (X=1/4)"),
         "Continuo a nivel de ULP (4.84e-15)."),

        ("3_secundaria", "cambio_de_rama", "maple/gga_exc/gga_x_gg99.mpl", 74, "gg99_r2 / gg99_asinh_sq / gg99_f0",
         "p = xs^2, v, w", "p = b^(2/3) = 35.87 (xs = 5.989); v = 1/16; w = 1/16",
         ["gga_x_gg99", "gga_x_kgg99"],
         "Si para p = b^(2/3) (cambio de rama Cardano real <-> trigonometrico); numerico para v,w = 1/16",
         "Gilbert & Gill, CPL 312, 511 (1999)",
         ("maple/gga_exc/gga_x_gg99.mpl:74", "GG99 (gga_x_gg99)", "gg99_r2 (p = b^(2/3))"),
         "Continuo (zk: 7.43e-15, vsigma: 1.05e-10 por raiz cuadrada en discriminante cubico p^3 = b^2)."),

        ("3_secundaria", "cambio_de_rama", "maple/gga_exc/gga_x_htbs.mpl", 40, "htbs_f0 (interpolacion esplines s = 0.6 y 2.6)",
         "s", "s = 0.6 y s = 2.6",
         ["gga_x_htbs"], "Si (definido por tramos en Ec. 7-9 de Haas et al. 2011)",
         "Haas, Tran, Blaha & Schwarz, PRB 83, 205117 (2011)",
         ("maple/gga_exc/gga_x_htbs.mpl:40", "HTBS (gga_x_htbs)", "htbs_f0 (s=2.6)"),
         "Continuo C1 por construccion (zk: 2.12e-14, vsigma: 2.98e-13)."),

        ("3_secundaria", "cambio_de_rama", "maple/gga_exc/gga_x_lv_rpw86.mpl", 33, "enforce_smooth_lr(lv_A/B, t, 100^2, 10)",
         "s", "s = 100.0 (expansion asintotica orden 10)",
         ["gga_x_lv_rpw86"], "No (evita cancelacion a s grande)",
         "Berland & Hyldgaard, PRB 89, 035412 (2014)",
         ("maple/gga_exc/gga_x_lv_rpw86.mpl:33", "LV-rPW86 (gga_x_lv_rpw86)", "enforce_smooth_lr (s=100)"),
         "Continuo a precision de maquina (7.99e-14 en zk, 3.19e-13 en vsigma)."),

        ("3_secundaria", "cambio_de_rama", "maple/gga_exc/gga_x_ncap.mpl", 34, "ncap_tanh_asinh (Taylor orden 16)",
         "s^2", "s^2 = 1/16 (s = 0.25)",
         ["gga_x_ncap", "hyb_gga_xc_ncap"], "No (evita 1/s en derivadas pares)",
         "Carmona-Espindola et al., JCP 150, 054105 (2019)",
         ("maple/gga_exc/gga_x_ncap.mpl:34", "NCAP (gga_x_ncap)", "ncap_tanh_asinh (s^2=1/16)"),
         "Continuo a nivel de ULP (2.17e-14)."),

        ("3_secundaria", "cambio_de_rama", "maple/gga_exc/gga_x_pbe_erf_gws.mpl", 64, "pbe_x_erf_gws_num_smooth",
         "mu_t = nu/2", "mu_t = 1/16, 1/2, 2^24",
         ["gga_x_pbe_erf_gws", "hyb_gga_xc_pbe_erf_gws"], "No (empalmes enforce_smooth)",
         "Goll, Werner & Stoll, PCCP 7, 3917 (2005)",
         ("maple/gga_exc/gga_x_pbe_erf_gws.mpl:64", "PBE-erf-GWS (gga_x_pbe_erf_gws)", "gws_num_smooth (mid=1/2)"),
         "Continuo a ~1e-12..3e-12."),

        ("3_secundaria", "cambio_de_rama", "maple/gga_exc/gga_x_pbepow.mpl", 35, "pbepow_f0 (suma finita 100 terminos, s^2 = 1)",
         "s^2", "s^2 = 1.0 (m = 100 terminos)",
         ["gga_x_pbepow"], "Si para m = 100; numerico para la reescritura en s^2 = 1",
         "Bremond, Tognetti, Chermette & Adamo, JCP 145, 224102 (2016)",
         ("maple/gga_exc/gga_x_pbepow.mpl:35", "PBEpow (gga_x_pbepow)", "pbepow_f0 (s^2=1)"),
         "Continuo a precision de maquina (6.33e-14 en zk, 1.56e-13 en vsigma)."),

        ("3_secundaria", "cambio_de_rama", "maple/gga_exc/gga_x_pbetrans.mpl", 50, "pbetrans_f0 (s^2 = 1)",
         "s^2", "s^2 = 1.0",
         ["gga_x_pbetrans"], "No (reescritura estable del sigmoide Fermi-Dirac)",
         "Bremond et al., JCTC 12, 459 (2016)",
         ("maple/gga_exc/gga_x_pbetrans.mpl:50", "PBEtrans (gga_x_pbetrans)", "pbetrans_f0 (s^2=1)"),
         "Continuo a precision de maquina (5.32e-14)."),

        ("3_secundaria", "cambio_de_rama", "src/hyb_gga_xc_case21.c", 53, "xc_bspline (10 B-splines cubicos uniformes)",
         "u = gammax*s^2/(1+gammax*s^2)", "Nudos u = m/7 (k=3, N=10, PMAX=8)",
         ["hyb_gga_xc_case21", "gga_xc_ Adamson?"], "Si (definicion B-spline cubica C2 del funcional CASE21)",
         "Sparrow et al., JPCL 13, 6896 (2022)",
         ("src/hyb_gga_xc_case21.c:53", "CASE21 (hyb_gga_xc_case21)", "xc_bspline (nudo u=3/7)"),
         "Continuo C2 en los nudos (zk: 2.35e-14, vsigma: 2.16e-13)."),

        ("3_secundaria", "cambio_de_rama", "maple/mgga_exc/mgga_x_gx.mpl", 37, "gx_f_a (empalme en alpha = 1 de PBE-GX)",
         "alpha", "alpha = 1.0",
         ["mgga_x_gx", "mgga_x_pbe_gx"], "Si (quiebre C1 de primera derivada presente en la definicion de Loos 2017)",
         "Loos, JCP 146, 114108 (2017)",
         ("maple/mgga_exc/mgga_x_gx.mpl:37", "PBE-GX (mgga_x_gx)", "gx_f_a (alpha=1, diseno)"),
         "Control conocido: zk continuo (5.28e-14), vtau salta 8.36e-1 por diseno del funcional."),

        ("3_secundaria", "cambio_de_rama", "maple/mgga_exc/mgga_x_regtpss.mpl", 37, "reg_one_minus_f_a y reg_one_minus_A",
         "alpha, s", "|alpha| = 0.5 y 3*s^2 = 1 (s = 0.57735)",
         ["mgga_x_regtpss"], "No (reescritura libre de cancelacion cerca de alpha=1 y s=0)",
         "Ruzsinszky et al., JCTC 8, 2078 (2012)",
         ("maple/mgga_exc/mgga_x_regtpss.mpl:37", "regTPSS_X (mgga_x_regtpss)", "reg_one_minus_f_a (|a|=0.5)"),
         "Continuo a precision de maquina (1.93e-14 y 2.26e-14)."),

        ("3_secundaria", "cambio_de_rama", "maple/mgga_exc/mgga_x_edmgga.mpl", 50, "edmgga_x (serie asintotica Qb < -1e4/a)",
         "Qb", "Qb = -1e4 / 0.704 = -14204.5 (nota: errata inocua -1e-4/a en m_min de linea 51)",
         ["mgga_x_edmgga", "hyb_mgga_xc_edmggah"], "No (evita cancelacion y + sqrt(1+y^2) en y -> -inf; orden 8)",
         "Tao, JCP 115, 3519 (2001)",
         ("maple/mgga_exc/mgga_x_edmgga.mpl:50", "EDMGGA (mgga_x_edmgga)", "edmgga_x (Qb = -1e4/a)"),
         "zk continuo (1.85e-13); pequeno salto residual 2.36e-8 en vsigma en Qb = -14204.5 por truncamiento/redondeo."),

        ("3_secundaria", "cambio_de_rama", "maple/mgga_exc/mgga_x_eel.mpl", 24, "eel_G (corte plano x <= x0)",
         "x = s^2", "x = x0 = 0.14118848 (s = 0.37575)",
         ["mgga_x_eel"], "Si (transicion C^inf tipo exp(-c*(x-x0)^(-1/4)) en Ec. 14 de Aschebrock 2023)",
         "Aschebrock, Lebeda & Kummel, JCP 159, 234107 (2023)",
         ("maple/mgga_exc/mgga_x_eel.mpl:24", "EEL (mgga_x_eel)", "eel_G (s^2 = x0)"),
         "C^inf exacto a precision de maquina (zk: 1.13e-13, vtau: 3.06e-14)."),

        ("3_secundaria", "cambio_de_rama", "maple/mgga_exc/mgga_x_ft98.mpl", 58, "ft98_q20 (series en q3=0 y -inf, orden 9)",
         "q3 = xs^4 - u^2 - b2", "|q3| = eps^(1/4) (~1.221e-4) y q3 = -eps^(-1/4) (-8192)",
         ["mgga_x_ft98"], "No (estabilizacion numerica de 1/(q3 + sqrt(1 + q3^2)))",
         "Filatov & Thiel, Mol. Phys. 95, 189 (1998)",
         ("maple/mgga_exc/mgga_x_ft98.mpl:60", "FT98 (mgga_x_ft98)", "ft98_q20 (|q3| = eps^(1/4))"),
         "Continuo a precision de maquina (0.00e+00)."),

        ("3_secundaria", "cambio_de_rama", "maple/mgga_exc/mgga_x_mbrxc_bg.mpl", 56, "mbrxc_v (Taylor orden 7)",
         "x(Q)", "x = sqrt(eps) (~1.49e-8)",
         ["mgga_x_mbrxc_bg", "mgga_x_mbrxh_bg"], "No (evita 0/0 en x = 0; inalcanzable para tau >= 0 donde x >= 1.15)",
         "Patra, Jana & Samal, JCP 153, 184112 (2020)",
         ("maple/mgga_exc/mgga_x_mbrxc_bg.mpl:56", "MBRXC-BG (mgga_x_mbrxc_bg)", "mbrxc_v (x = sqrt(eps))"),
         "Truncamiento de 4a derivada ~4e-21 << eps; continuo a nivel de ULP."),

        ("3_secundaria", "umbral_o_recorte", "maple/mgga_exc/mgga_x_ms.mpl", 29, "ms_fa (recorte |alpha| <= 1e6)",
         "alpha", "|alpha| = ms_amax = 1e6",
         ["mgga_x_ms0", "mgga_x_ms1", "mgga_x_ms2", "hyb_mgga_x_ms2h", "mgga_x_mvs", "hyb_mgga_x_mvsh"],
         "No (evita desborde de alpha^8 en el cociente racional de Sun et al.)",
         "Sun, Haunschild et al., JCP 137, 051101 (2012)",
         ("maple/mgga_exc/mgga_x_ms.mpl:29", "MS0 (mgga_x_ms0)", "ms_fa clamp (alpha=1e6)"),
         "En alpha = 1e6 zk ya esta saturado a precision de maquina (0.00e+00); vtau ~ 1e-18 cae a 0."),

        ("3_secundaria", "cambio_de_rama", "maple/mbeef.mpl", 38, "mbeef_one_minus_xj (enforce_smooth_lr en alpha = 1e4)",
         "alpha", "alpha = 1e4 (orden 4) y expansion Chebyshev 8x8 (64 terminos) / 5x5 en mbeefvdw",
         ["mgga_x_mbeef", "mgga_x_mbeefvdw"], "No para el corte 1e4; Si para base polinomial 8x8 / 5x5",
         "Wellendorff et al., JCP 140, 144107 (2014); Lundgaard et al., PRB 93, 235162 (2016)",
         ("maple/mbeef.mpl:38", "mBEEF (mgga_x_mbeef)", "mbeef_one_minus_xj (alpha=1e4)"),
         "zk continuo (0.00e+00); salto absoluto despreciable en vtau (~1e-20, rel 2e-8 sobre cola nula)."),

        ("3_secundaria", "cambio_de_rama", "maple/mgga_exc/mgga_k_csk.mpl", 40, "csk_zI",
         "z_csk", "z = -(-ln(eps))^(-1/a) y z = -(-ln(1-eps))^(-1/a)",
         ["mgga_k_csk1", "mgga_k_csk4", "mgga_k_csk_loc1", "mgga_k_csk_loc4"],
         "No (cortes donde exp(-1/|z|^a) cae bajo eps o se satura a 1)",
         "Cancio, Stewart & Kuna, JCP 144, 084107 (2016)",
         ("maple/mgga_exc/mgga_k_csk.mpl:40", "CSK1 (mgga_k_csk1)", "csk_zI (z = -1/(-log(eps)))"),
         "Continuo a precision de maquina (zk: 2.06e-13, vsigma: 6.47e-14)."),

        ("3_secundaria", "cambio_de_rama", "maple/mgga_exc/mgga_k_rpp.mpl", 65, "rpp_ap (regularizacion por tramos x in [0, x0])",
         "x", "x = 0 y x = x0 = 0.819411 (polinomio grado 5 C3)",
         ["mgga_k_rpp", "mgga_x_ofr2"], "Si (construccion C3 del funcional RPP / OFR2)",
         "Kaplan & Perdew, PRR 4, 083803 (2022)",
         ("maple/mgga_exc/mgga_k_rpp.mpl:65", "RPP (mgga_k_rpp)", "rpp_ap (x = x0 = 0.819411)"),
         "Polinomio de grado 5 ajustado con continuidad C3 exacta en x = 0 y x = x0."),

        ("3_secundaria", "umbral_o_recorte", "maple/op.mpl", 18, "op_qab (guarda |b88_zab| < XC_EPSILON)",
         "b88_zab", "|zab| = XC_EPSILON (y 1 - |z| <= zeta_threshold)",
         "maple/op.mpl", "No (evita division por cero en op_qab cuando zab -> 0)",
         "Tsuneda, Suzumura & Hirao, JCP 110, 10664 (1999)",
         None, "Guarda puntual en |zab| < 2.22e-16."),

        ("3_secundaria", "umbral_o_recorte", "src/lda_xc_ksdt.c", 127, "Recorte de temperatura minima T >= 1e-8",
         "T", "params->T = max(T, 1e-8)",
         ["lda_xc_ksdt", "lda_xc_corrksdt", "gga_xc_kt1"], "No (evita 0/0 en t = T/T_F -> 0)",
         "Karasiev, Sjostrom, Dufty & Trickey, PRL 112, 076403 (2014)",
         None, "Piso fijo sobre el parametro externo de temperatura reducida."),

        ("3_secundaria", "umbral_o_recorte", "maple/mgga_exc/mgga_x_2d_prhg07.mpl", 19, "prhg07_y (argumento de LambertW)",
         "x / Pi", "x / Pi >= -0.9999999999 (cerca del punto de ramificacion -1/e)",
         ["mgga_x_2d_prhg07", "mgga_x_2d_prhg07_prp10"],
         "No (evita NaN en LambertW(z) para z < -1/e; truncado a 10 nueves)",
         "Pittalis et al., PRB 76, 235314 (2007)",
         None, "Constante cableada -0.9999999999 * exp(-1) protege contra argumento complejo en W_0."),

        ("3_secundaria", "constante_cableada", "maple/mgga_exc/mgga_c_rscan.mpl", 28, "Regularizadores taur y alphar de rSCAN",
         "taur, alphar", "taur = 1e-4, alphar = 1e-3",
         ["mgga_x_rscan", "mgga_c_rscan", "hyb_mgga_x_wr2scan"],
         "Si (parametros de regularizacion publicados por Bartok & Yates 2019)",
         "Bartok & Yates, JCP 150, 161101 (2019)",
         None, "Parte de la definicion de rSCAN (conocido por romper invariancia de escalamiento uniforme, motivando r2SCAN)."),

        ("3_secundaria", "terminos_o_iteraciones_fijos", "maple/gga_exc/gga_x_beefvdw.mpl", 26, "Expansiones Legendre / Chebyshev de funcionales ML/empiricos",
         "grado", "BEEF-vdW: 30 Legendre; mBEEF / MCML / VCML: 8x8=64 Chebyshev; B97/SOGGA11/N12: 4-6 terminos",
         ["gga_x_beefvdw", "mgga_x_mcml", "mgga_x_vcml", "gga_x_sogga11", "gga_x_n12"],
         "Si (definicion parametrica de cada funcional en su publicacion)",
         "Wellendorff et al., PRB 85, 235149 (2012)",
         None, "Ordenes polinomiales fijos propios del modelo fisico."),
    ]

    out_rows = []
    for idx, item in enumerate(raw, 1):
        prio, cat, fpath, line, sym, var, cut_val, fspec, justif, ref, mkey, notes = item
        fids = resolve_fids(fspec)
        n_fids, prio_str, ex_str = summarize_fids(fids)
        m = meas.get(mkey) if mkey else None
        out_rows.append({
            "id_entrada": f"INV-{idx:02d}",
            "prioridad": prio,
            "categoria": cat,
            "archivo": fpath,
            "linea": line,
            "simbolo_o_funcion": sym,
            "variable_corte": var,
            "valor_corte_o_parametro": cut_val,
            "n_funcionales_afectados": n_fids,
            "funcionales_prioritarios": prio_str,
            "funcionales_ejemplo": ex_str,
            "justificado_en_articulo": justif,
            "referencia_original": ref,
            "probado_en_frontera": m["probado"] if m else "N/A (umbral/orden)",
            "salto_rel_valor_1e9": m["val_1e9"] if m else "",
            "salto_rel_valor_1e14": m["val_1e14"] if m else "",
            "salto_rel_deriv_1e9": m["der_1e9"] if m else "",
            "salto_rel_deriv_1e14": m["der_1e14"] if m else "",
            "canal_peor_deriv": m["canal"] if m else "",
            "observaciones_y_estado": notes,
        })

    out_path = ROOT / "inventario_umbrales_libxc.csv"
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"Wrote {len(out_rows)} inventory entries to {out_path}")


if __name__ == "__main__":
    main()
