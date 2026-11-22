"""Thin ctypes harness over the raw libxc C API (bundled with PySCF)."""
import ctypes, glob, os
import numpy as np

_libpath = os.environ.get("LIBXC_PATH") or glob.glob(os.path.join(os.path.dirname(__import__("pyscf").__file__), "lib", "deps", "lib*", "libxc.so*"))[0]
L = ctypes.CDLL(_libpath)

dp = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
NULL = None

L.xc_func_alloc.restype = ctypes.c_void_p
L.xc_func_init.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
L.xc_func_end.argtypes = [ctypes.c_void_p]
L.xc_func_free.argtypes = [ctypes.c_void_p]
L.xc_func_get_info.argtypes = [ctypes.c_void_p]
L.xc_func_get_info.restype = ctypes.c_void_p
L.xc_func_info_get_family.argtypes = [ctypes.c_void_p]
L.xc_func_info_get_flags.argtypes = [ctypes.c_void_p]
L.xc_func_info_get_kind.argtypes = [ctypes.c_void_p]
L.xc_functional_get_name.restype = ctypes.c_char_p
L.xc_functional_get_name.argtypes = [ctypes.c_int]
L.xc_version_string.restype = ctypes.c_char_p
L.xc_func_set_dens_threshold.argtypes = [ctypes.c_void_p, ctypes.c_double]

FAM = {1: "LDA", 2: "GGA", 4: "MGGA", 32: "GGA", 64: "MGGA", 128: "LDA"}
FLAG_EXC, FLAG_VXC, FLAG_FXC = 1, 2, 4
FLAG_NEEDS_LAPL = 1 << 15


def all_ids():
    n = L.xc_number_of_functionals()
    arr = (ctypes.c_int * n)()
    L.xc_available_functional_numbers(arr)
    return list(arr)


def name(i):
    return L.xc_functional_get_name(i).decode()


class Func:
    def __init__(self, fid, nspin):
        self.p = L.xc_func_alloc()
        if L.xc_func_init(self.p, fid, nspin) != 0:
            L.xc_func_free(self.p)
            self.p = None
            raise RuntimeError("init failed")
        info = L.xc_func_get_info(self.p)
        self.family = FAM.get(L.xc_func_info_get_family(info), "OTHER")
        self.flags = L.xc_func_info_get_flags(info)
        self.kind = L.xc_func_info_get_kind(info)
        self.nspin = nspin
        self.fid = fid

    def close(self):
        if self.p:
            L.xc_func_end(self.p); L.xc_func_free(self.p); self.p = None

    # x: dict with rho (np,ns), sigma (np,nsig), lapl (np,ns), tau (np,ns)
    def exc_vxc(self, x):
        rho = np.ascontiguousarray(x["rho"]); n = rho.shape[0]
        ns = self.nspin; nsig = 1 if ns == 1 else 3
        zk = np.zeros(n); vrho = np.zeros((n, ns))
        out = {"zk": zk, "vrho": vrho}
        if self.family == "LDA":
            L.xc_lda_exc_vxc(ctypes.c_void_p(self.p), ctypes.c_size_t(n), _p(rho), _p(zk), _p(vrho))
        elif self.family == "GGA":
            sig = np.ascontiguousarray(x["sigma"]); vs = np.zeros((n, nsig))
            L.xc_gga_exc_vxc(ctypes.c_void_p(self.p), ctypes.c_size_t(n), _p(rho), _p(sig), _p(zk), _p(vrho), _p(vs))
            out["vsigma"] = vs
        elif self.family == "MGGA":
            sig = np.ascontiguousarray(x["sigma"]); lap = np.ascontiguousarray(x["lapl"]); tau = np.ascontiguousarray(x["tau"])
            vs = np.zeros((n, nsig)); vl = np.zeros((n, ns)); vt = np.zeros((n, ns))
            L.xc_mgga_exc_vxc(ctypes.c_void_p(self.p), ctypes.c_size_t(n), _p(rho), _p(sig), _p(lap), _p(tau),
                              _p(zk), _p(vrho), _p(vs), _p(vl), _p(vt))
            out.update(vsigma=vs, vlapl=vl, vtau=vt)
        else:
            raise RuntimeError("family")
        return out


def _p(a):
    return a.ctypes.data_as(ctypes.c_void_p)
