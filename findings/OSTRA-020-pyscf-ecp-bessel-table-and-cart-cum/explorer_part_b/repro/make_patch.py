import pathlib, sys

src_dir = pathlib.Path(sys.argv[1])
dst_dir = pathlib.Path(sys.argv[2])

# 1. Patch nr_ecp.h: CART_CUM 456 (l<=12) -> 816 (l<=15, matching ANG_MAX=15 in cint.h)
sh = (src_dir / "nr_ecp.h").read_text()
old_h = "#define CART_CUM        (455+1) // upto l = 12"
new_h = "#define CART_CUM        816     // upto l = 15 (matching ANG_MAX in cint.h)"
assert old_h in sh
(dst_dir / "nr_ecp.h").write_text(sh.replace(old_h, new_h, 1))

# 2. Patch nr_ecp.c
sc = (src_dir / "nr_ecp.c").read_text()

old_py = "        4, 3, 2, 1, 0,14,13,12,11,10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,\n};"
new_py = "        4, 3, 2, 1, 0,14,13,12,11,10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,\n       15,14,13,12,11,10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,\n};"
assert sc.count(old_py) == 1
sc = sc.replace(old_py, new_py, 1)

old_pz = "        9,10,11,12,13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,\n};"
new_pz = "        9,10,11,12,13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,\n        0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,\n};"
assert sc.count(old_pz) == 1
sc = sc.replace(old_pz, new_pz, 1)

old_off = """static int _offset_cart[] = {0, 1, 4, 10, 20, 35, 56, 84, 120,
                             165, 220, 286, 364, 455, 560};"""
new_off = """static int _offset_cart[] = {0, 1, 4, 10, 20, 35, 56, 84, 120,
                             165, 220, 286, 364, 455, 560, 680, 816};"""
assert old_off in sc
sc = sc.replace(old_off, new_off, 1)

old_t2 = """        double *rur, *radi, *radj, *angi, *angj, *rad_all, *buf;
        double *prad, *plast;
        MALLOC_INSTACK(rad_all, nci*ncj*(li+lj+1)*(li+ECP_LMAX+1)*(lj+ECP_LMAX+1));
        MARK_STACK;
        MALLOC_INSTACK(angi, (li+1)*nfi*(ECP_LMAX*2+1)*(li+ECP_LMAX+1));
        MALLOC_INSTACK(angj, (lj+1)*nfj*(ECP_LMAX*2+1)*(lj+ECP_LMAX+1));
        MALLOC_INSTACK(buf, nfi*(ECP_LMAX*2+1)*(lj+ECP_LMAX+1));

        RESTORE_STACK;
        MALLOC_INSTACK(rur, nrs*(li+lj+1));
        MALLOC_INSTACK(radi, nci*(li+ECP_LMAX+1)*nrs);
        MALLOC_INSTACK(radj, ncj*(lj+ECP_LMAX+1)*nrs);
        int lilj1 = li + lj + 1;
        int lilc1 = li + ECP_LMAX + 1;
        int ljlc1 = lj + ECP_LMAX + 1;"""
new_t2 = """        int ecp_lmax = ECP_LMAX;
        for (iloc = 0; iloc < necpbas; iloc++) {
                int lc0 = ecpbas[ANG_OF+iloc*BAS_SLOTS];
                if (lc0 > ecp_lmax) { ecp_lmax = lc0; }
        }
        double *rur, *radi, *radj, *angi, *angj, *rad_all, *buf;
        double *prad, *plast;
        MALLOC_INSTACK(rad_all, nci*ncj*(li+lj+1)*(li+ecp_lmax+1)*(lj+ecp_lmax+1));
        MARK_STACK;
        MALLOC_INSTACK(angi, (li+1)*nfi*(ecp_lmax*2+1)*(li+ecp_lmax+1));
        MALLOC_INSTACK(angj, (lj+1)*nfj*(ecp_lmax*2+1)*(lj+ecp_lmax+1));
        MALLOC_INSTACK(buf, nfi*(ecp_lmax*2+1)*(lj+ecp_lmax+1));
        double *cache_ang = cache;

        RESTORE_STACK;
        MALLOC_INSTACK(rur, nrs*(li+lj+1));
        MALLOC_INSTACK(radi, nci*(li+ecp_lmax+1)*nrs);
        MALLOC_INSTACK(radj, ncj*(lj+ecp_lmax+1)*nrs);
        int lilj1 = li + lj + 1;
        int lilc1 = li + ecp_lmax + 1;
        int ljlc1 = lj + ecp_lmax + 1;"""
assert old_t2 in sc
sc = sc.replace(old_t2, new_t2, 1)

old_t2_ang = """        type2_facs_ang(angi, li, lc, rca, cache);
        type2_facs_ang(angj, lj, lc, rcb, cache);"""
new_t2_ang = """        double *ang_work = (cache > cache_ang) ? cache : cache_ang;
        type2_facs_ang(angi, li, lc, rca, ang_work);
        type2_facs_ang(angj, lj, lc, rcb, ang_work);"""
assert old_t2_ang in sc
sc = sc.replace(old_t2_ang, new_t2_ang, 1)

old_t1_z = "        for (i = 0; i < d3; i++) { rad_all[i] = 0; }"
new_t1_z = "        for (i = 0; i < npi*npj*d2; i++) { rad_all[i] = 0; }"
assert sc.count(old_t1_z) == 1
sc = sc.replace(old_t1_z, new_t1_z, 1)

old_csz = """int ECPscalar_cache_size(int comp, int *shls,
                         int *atm, int natm, int *bas, int nbas, double *env)
{
        const int ish = shls[0];
        const int jsh = shls[1];
        const int li = bas[ANG_OF+ish*BAS_SLOTS];
        const int lj = bas[ANG_OF+jsh*BAS_SLOTS];
        const int nfi = (li+1) * (li+2) / 2;
        const int nfj = (lj+1) * (lj+2) / 2;
        const int npi = bas[NPRIM_OF+ish*BAS_SLOTS];
        const int npj = bas[NPRIM_OF+jsh*BAS_SLOTS];
        const int nci = bas[NCTR_OF+ish*BAS_SLOTS];
        const int ncj = bas[NCTR_OF+jsh*BAS_SLOTS];
        const int lilj1 = li + lj + 1;
        const int lilc1 = li + ECP_LMAX + 1;
        const int ljlc1 = lj + ECP_LMAX + 1;
        const int d1 = lilj1;
        const int d2 = d1 * d1;
        const int d3 = d2 * d1;
        const int di1 = li + 1;
        const int di2 = di1 * di1;
        const int di3 = di2 * di1;
        const int dj1 = lj + 1;
        const int dj2 = dj1 * dj1;
        const int dj3 = dj2 * dj1;
        const int nrs =  1 << LEVEL_MAX;
        int size1 = nci*ncj*(li+lj+1);
        size1 += ((li+1)*nfi*(ECP_LMAX*2+1)*(li+ECP_LMAX+1) +
                  (lj+1)*nfj*(ECP_LMAX*2+1)*(lj+ECP_LMAX+1));
        size1 += nfi*(ECP_LMAX*2+1)*(lj+ECP_LMAX+1);
        size1 += (lj+1)*nfj*(ECP_LMAX*2+1)*(lj+ECP_LMAX+1)*3;
        int size2 = nrs * (li+lj+1 + 1 +
                           nci*(li+ECP_LMAX+1) +
                           ncj*(lj+ECP_LMAX+1) +
                           MAX(npi*lilc1, npj*ljlc1));
        size2 += lilc1 * ljlc1;
        size2 += MAX(di1*di1*di1*lilc1, dj1*dj1*dj1*ljlc1) * (ECP_LMAX*2 + 1);
        int size = nfi*nfj*(nci*ncj+2) * comp;
        size += nci*ncj*(li+lj+1)*(li+ECP_LMAX+1)*(lj+ECP_LMAX+1);
        //size += MAX(size1, size2);  bugs in bufsize estimation, not sure where's the error
        size += size1 + size2 + 120;
        size += nfi*(ECP_LMAX*2+1)*(lj+ECP_LMAX+1);
        size += npi*npj*d2;
        size += d3;
        size += nci*ncj*d3;
        size += nfi*di3;
        size += nfj*dj3;
        size += npi*nci;
        size += npj*ncj;
        size += natm;
        return size;
}"""
new_csz = """int ECPscalar_cache_size_deriv(int comp, int dli, int dlj, int *shls,
                               int *atm, int natm, int *bas, int nbas, double *env)
{
        const int ish = shls[0];
        const int jsh = shls[1];
        const int li = bas[ANG_OF+ish*BAS_SLOTS] + dli;
        const int lj = bas[ANG_OF+jsh*BAS_SLOTS] + dlj;
        const int nfi = (li+1) * (li+2) / 2;
        const int nfj = (lj+1) * (lj+2) / 2;
        const int npi = bas[NPRIM_OF+ish*BAS_SLOTS];
        const int npj = bas[NPRIM_OF+jsh*BAS_SLOTS];
        const int nci = bas[NCTR_OF+ish*BAS_SLOTS];
        const int ncj = bas[NCTR_OF+jsh*BAS_SLOTS];
        const int necpbas = (int)env[AS_NECPBAS];
        const int *ecpbas = bas + ((int)env[AS_ECPBAS_OFFSET])*BAS_SLOTS;
        int ecp_lmax = ECP_LMAX;
        int i;
        for (i = 0; i < necpbas; i++) {
                int lc0 = ecpbas[ANG_OF+i*BAS_SLOTS];
                if (lc0 > ecp_lmax) { ecp_lmax = lc0; }
        }
        const int lilj1 = li + lj + 1;
        const int lilc1 = li + ecp_lmax + 1;
        const int ljlc1 = lj + ecp_lmax + 1;
        const int d1 = lilj1;
        const int d2 = d1 * d1;
        const int d3 = d2 * d1;
        const int di1 = li + 1;
        const int di2 = di1 * di1;
        const int di3 = di2 * di1;
        const int dj1 = lj + 1;
        const int dj2 = dj1 * dj1;
        const int dj3 = dj2 * dj1;
        const int nrs =  1 << LEVEL_MAX;
        int size1 = nci*ncj*(li+lj+1);
        size1 += ((li+1)*nfi*(ecp_lmax*2+1)*(li+ecp_lmax+1) +
                  (lj+1)*nfj*(ecp_lmax*2+1)*(lj+ecp_lmax+1));
        size1 += nfi*(ecp_lmax*2+1)*(lj+ecp_lmax+1);
        size1 += (lj+1)*nfj*(ecp_lmax*2+1)*(lj+ecp_lmax+1)*3;
        int size2 = nrs * (li+lj+1 + 1 +
                           nci*(li+ecp_lmax+1) +
                           ncj*(lj+ecp_lmax+1) +
                           MAX(npi*lilc1, npj*ljlc1));
        size2 += lilc1 * ljlc1;
        size2 += MAX(di1*di1*di1*lilc1, dj1*dj1*dj1*ljlc1) * (ecp_lmax*2 + 1);
        size2 += 3 * MAX(di2, dj2);
        int size = nfi*nfj*(nci*ncj+2) * (comp + 4);
        size += nci*ncj*(li+lj+1)*(li+ecp_lmax+1)*(lj+ecp_lmax+1);
        size += size1 + size2 + 120;
        size += nfi*(ecp_lmax*2+1)*(lj+ecp_lmax+1);
        size += npi*npj*d2 + d2;
        size += d3;
        size += nci*ncj*d3;
        size += nfi*di3;
        size += nfj*dj3;
        size += npi*nci + npj*ncj + npi*npj;
        size += (npi + npj) * BAS_SLOTS + necpbas + 1;
        size += natm;
        return size;
}

int ECPscalar_cache_size(int comp, int *shls,
                         int *atm, int natm, int *bas, int nbas, double *env)
{
        return ECPscalar_cache_size_deriv(comp, 0, 0, shls, atm, natm, bas, nbas, env);
}"""
assert old_csz in sc
sc = sc.replace(old_csz, new_csz, 1)
(dst_dir / "nr_ecp.c").write_text(sc)

# 3. Patch nr_ecp_deriv.c
sd = (src_dir / "nr_ecp_deriv.c").read_text()
old_decl = """int ECPscalar_cache_size(int comp, int *shls,
                         int *atm, int natm, int *bas, int nbas, double *env);"""
new_decl = """int ECPscalar_cache_size(int comp, int *shls,
                         int *atm, int natm, int *bas, int nbas, double *env);
int ECPscalar_cache_size_deriv(int comp, int dli, int dlj, int *shls,
                               int *atm, int natm, int *bas, int nbas, double *env);"""
assert old_decl in sd
sd = sd.replace(old_decl, new_decl, 1)
sd = sd.replace(
"""static int _cart_factory(Function_cart intor_cart, double *out, int comp,
                         int *dims, int *shls, int *ecpbas, int necpbas,""",
"""static int _cart_factory(Function_cart intor_cart, double *out, int comp,
                         int dli, int dlj, int *dims, int *shls, int *ecpbas, int necpbas,"""
)
sd = sd.replace(
"""        if (out == NULL) {
                int cache_size = ECPscalar_cache_size(comp*2, shls,
                                                      atm, natm, bas, nbas, env);
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                int cache_size = ECPscalar_cache_size(comp*2, shls,
                                                      atm, natm, bas, nbas, env);""",
"""        if (out == NULL) {
                int cache_size = ECPscalar_cache_size_deriv(comp*2, dli, dlj, shls,
                                                            atm, natm, bas, nbas, env);
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                int cache_size = ECPscalar_cache_size_deriv(comp*2, dli, dlj, shls,
                                                            atm, natm, bas, nbas, env);"""
)
sd = sd.replace(
"""static int _sph_factory(Function_cart intor_cart, double *out, int comp,
                        int *dims, int *shls, int *ecpbas, int necpbas,""",
"""static int _sph_factory(Function_cart intor_cart, double *out, int comp,
                        int dli, int dlj, int *dims, int *shls, int *ecpbas, int necpbas,"""
)
sd = sd.replace(
"""        if (out == NULL) {
                int cache_size = ECPscalar_cache_size(comp*2+2, shls,
                                                      atm, natm, bas, nbas, env);
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                int cache_size = ECPscalar_cache_size(comp*2+2, shls,
                                                      atm, natm, bas, nbas, env);""",
"""        if (out == NULL) {
                int cache_size = ECPscalar_cache_size_deriv(comp*2+2, dli, dlj, shls,
                                                            atm, natm, bas, nbas, env);
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                int cache_size = ECPscalar_cache_size_deriv(comp*2+2, dli, dlj, shls,
                                                            atm, natm, bas, nbas, env);"""
)
sd = sd.replace("_cart_factory(_deriv1_cart, out, 3,", "_cart_factory(_deriv1_cart, out, 3, 1, 0,")
sd = sd.replace("_sph_factory(_deriv1_cart, out, 3,", "_sph_factory(_deriv1_cart, out, 3, 1, 0,")
sd = sd.replace("_cart_factory(_ipipv_cart, out, 9,", "_cart_factory(_ipipv_cart, out, 9, 2, 0,")
sd = sd.replace("_sph_factory(_ipipv_cart, out, 9,", "_sph_factory(_ipipv_cart, out, 9, 2, 0,")
sd = sd.replace("_cart_factory(_ipvip_cart, out, 9,", "_cart_factory(_ipvip_cart, out, 9, 1, 1,")
sd = sd.replace("_sph_factory(_ipvip_cart, out, 9,", "_sph_factory(_ipvip_cart, out, 9, 1, 1,")
sd = sd.replace("_cart_factory(_igv_cart, out, comp,", "_cart_factory(_igv_cart, out, comp, 1, 0,")
sd = sd.replace("_sph_factory(_igv_cart, out, comp,", "_sph_factory(_igv_cart, out, comp, 1, 0,")
(dst_dir / "nr_ecp_deriv.c").write_text(sd)
print("Patch applied cleanly!")
