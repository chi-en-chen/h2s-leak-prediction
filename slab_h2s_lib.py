# ============================================================
#  slab_h2s_lib.py
#  基于 SLAB 重气扩散模型的硫化氢(H2S)扩散分析工具库
#  ------------------------------------------------------------
#  功能:
#   1) 生成 SLAB 输入文件(按 H2S 物性)
#   2) 调用 Slab.exe 计算
#   3) 解析 PREDICT 输出中的"浓度等值线参数"表
#   4) 依据 SLAB 文档公式重构浓度场 C(x,y,z,t)
#      c(x,y,z,t) = cc(x) * (erf(xa)-erf(xb)) * (erf(ya)-erf(yb))
#                    * (exp(-za^2)+exp(-zb^2))
#      xa=(x-xc(t)+bx(t))/(sqrt(2)*betax(t))   ya=(y+b(x))/(sqrt(2)*betac(x))
#      xb=(x-xc(t)-bx(t))/(sqrt(2)*betax(t))   yb=(y-b(x))/(sqrt(2)*betac(x))
#      za=(z-zc(x))/(sqrt(2)*sig(x))           zb=(z+zc(x))/(sqrt(2)*sig(x))
#  单位: 浓度为体积分数(无量纲), 1 ppm = 1e-6 体积分数
# ============================================================
import os
import re
import shutil
import subprocess
import tempfile

import numpy as np
from scipy.special import erf

# ------------------------------------------------------------
# H2S 物性参数 (kg, m, s, K, J 单位制)
# ------------------------------------------------------------
H2S_PROPS = {
    "wms":   0.034081,   # 分子量 (kg/mol)
    "cps":   1010.0,     # 气相定压比热 (J/kg/K)
    "tbp":   213.0,      # 沸点 (K), -60.2 C
    "cmed0": 0.0,        # 初始液相质量分数
    "dhe":   548000.0,   # 汽化潜热 (J/kg)
    "cpsl":  1800.0,     # 液相比热 (J/kg/K)
    "rhosl": 949.0,      # 液相密度 (kg/m3)
    "spb":   -1.0,       # 饱和蒸汽压常数, -1=按沸点自动推算
    "spc":   0.0,
}

# ------------------------------------------------------------
# 输入文件生成
# ------------------------------------------------------------
def _fmt10(v):
    """把数值格式化成不超过10列的字符串, 保留尽量多有效数字.
    FORTRAN 按 F10.3 读取: 字段内含小数点时按实际小数读取."""
    for nd in (6, 5, 4, 3, 2, 1):
        s = f"{v:>10.{nd}f}"
        if len(s) <= 10:
            return s
    return f"{v:>10.1f}"


def write_slab_input(params, path, props=None):
    """按 SLAB 输入顺序写出输入文件.

    params: dict, 至少包含
        idspl   源类型: 1=蒸发池 2=水平射流 3=垂直喷口 4=瞬时蒸发池
        ncalc   计算子步乘数(默认1)
        ts      源温度 K
        qs      泄露速率 kg/s
        as      源面积 m2 (蒸发池面积/喷口截面积)
        tsd     连续泄露持续时间 s
        qtis    瞬时源质量 kg
        hs      源高度 m
        us      水平射流速度 m/s (idspl=2)
        ws      垂直喷口速度 m/s (idspl=3)
        tav     浓度平均时间 s
        xffm    最大下风向距离 m
        zp      浓度计算高度 [4个]
        z0      地面粗糙度 m
        za      风速测量高度 m
        ua      环境风速 m/s
        ta      环境温度 K
        rh      相对湿度 %
        stab    大气稳定度: 0=按ala, 1..6=A..F
        ala     逆蒙宁-奥布霍夫长度倒数 1/m (仅 stab=0 时需要)
    props: 气体物性, 默认 H2S
    """
    p = dict(props or {})
    p.update(params)
    v = lambda k: p[k]
    lines = []
    lines.append(f"{int(p['idspl']):5d}")
    lines.append(f"{int(p.get('ncalc', 1)):5d}")
    for k in ("wms", "cps", "tbp", "cmed0", "dhe", "cpsl", "rhosl", "spb", "spc"):
        lines.append(_fmt10(v(k)))
    for k in ("ts", "qs", "as", "tsd", "qtis", "hs"):
        lines.append(_fmt10(v(k)))
    for k in ("tav", "xffm", "zp1", "zp2", "zp3", "zp4"):
        lines.append(_fmt10(v(k)))
    lines.append(_fmt10(v("z0")))
    for k in ("za", "ua", "ta", "rh", "stab"):
        lines.append(_fmt10(v(k)))
    if float(p.get("stab", 0.0)) == 0.0:
        lines.append(_fmt10(v("ala")))
    lines.append(_fmt10(-1.0))   # 结束标记(与自带算例一致, 实际不会被读取)
    with open(path, "w", encoding="ascii") as f:
        f.write("\n".join(lines) + "\n")
    return path


def default_h2s_params(**overrides):
    """H2S 检测井泄露默认工况, 可覆盖任意参数."""
    p = {
        # 源项(按检测井实际情况修改)
        "idspl": 1,          # 1=地面蒸发池(泄漏液体积聚蒸发)
        "ncalc": 1,
        "ts":    213.0,      # 液池温度=沸点
        "qs":    1.0,        # 泄露速率 kg/s
        "as":    100.0,      # 液池面积 m2
        "tsd":   600.0,      # 持续泄露时间 s
        "qtis":  0.0,
        "hs":    0.0,
        "us":    0.0,
        "ws":    0.0,
        # 计算域
        "tav":   10.0,       # 平均时间 s
        "xffm":  3000.0,     # 最大下风向距离 m
        "zp1":   1.5,        # 浓度计算高度(呼吸带)
        "zp2":   0.0,
        "zp3":   0.0,
        "zp4":   0.0,
        # 气象(现场条件)
        "z0":    0.1,        # 地面粗糙度 m (田野)
        "za":    10.0,       # 风速计高度 m
        "ua":    2.0,        # 10m 高度风速 m/s
        "ta":    300.0,      # 环境温度 K
        "rh":    50.0,       # 相对湿度 %
        "stab":  4.0,        # 稳定度 D(中性), 1..6=A..F
        "ala":   0.0,
    }
    p.update(overrides)
    return p


# ------------------------------------------------------------
# 运行 SLAB
# ------------------------------------------------------------
def run_slab(exe_path="Slab.exe", input_path=None, output_path=None,
             workdir=None):
    """复制输入文件为 input, 调用 Slab.exe, 将结果 PREDICT 改名.

    返回输出文件路径.
    """
    if workdir is None:
        workdir = os.getcwd()
    exe_abs = os.path.abspath(exe_path) if os.path.isabs(exe_path) \
        else os.path.abspath(os.path.join(workdir, exe_path))
    if input_path is None:
        input_path = os.path.join(workdir, "input")
    if output_path is None:
        output_path = os.path.join(workdir, "predict_h2s.txt")

    workdir = os.path.abspath(workdir)
    inp_dst = os.path.join(workdir, "input")
    predict = os.path.join(workdir, "PREDICT")

    shutil.copyfile(os.path.abspath(input_path), inp_dst)
    if os.path.exists(predict):
        os.remove(predict)     # SLAB 硬编码输出名, 运行前必须不存在
    r = subprocess.run([exe_abs], cwd=workdir, capture_output=True, text=True,
                       timeout=300)
    if r.returncode != 0 or not os.path.exists(predict):
        raise RuntimeError(f"SLAB 运行失败 rc={r.returncode}\n{r.stdout}\n{r.stderr}")
    os.replace(predict, os.path.abspath(output_path))
    return os.path.abspath(output_path)


# ------------------------------------------------------------
# 解析 SLAB 输出
# ------------------------------------------------------------
_NUM = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"


def _parse_echo(text):
    """解析文件开头的 problem input 回显, 返回参数字典."""
    meta = {}
    pats = {
        "idspl": r"idspl\s*=\s*(\d+)",
        "ncalc": r"ncalc\s*=\s*(\d+)",
        "wms":   r"wms\s*=\s*(" + _NUM + r")",
        "cps":   r"cps\s*=\s*(" + _NUM + r")",
        "tbp":   r"tbp\s*=\s*(" + _NUM + r")",
        "cmed0": r"cmed0\s*=\s*(" + _NUM + r")",
        "dhe":   r"dhe\s*=\s*(" + _NUM + r")",
        "cpsl":  r"cpsl\s*=\s*(" + _NUM + r")",
        "rhosl": r"rhosl\s*=\s*(" + _NUM + r")",
        "spb":   r"spb\s*=\s*(" + _NUM + r")",
        "spc":   r"spc\s*=\s*(" + _NUM + r")",
        "ts":    r"\bts\s*=\s*(" + _NUM + r")",
        "qs":    r"\bqs\s*=\s*(" + _NUM + r")",
        "as":    r"\bas\s*=\s*(" + _NUM + r")",
        "tsd":   r"tsd\s*=\s*(" + _NUM + r")",
        "qtis":  r"qtis\s*=\s*(" + _NUM + r")",
        "hs":    r"\bhs\s*=\s*(" + _NUM + r")",
        "tav":   r"tav\s*=\s*(" + _NUM + r")",
        "xffm":  r"xffm\s*=\s*(" + _NUM + r")",
        "z0":    r"z0\s*=\s*(" + _NUM + r")",
        "za":    r"za\s*=\s*(" + _NUM + r")",
        "ua":    r"ua\s*=\s*(" + _NUM + r")",
        "ta":    r"ta\s*=\s*(" + _NUM + r")",
        "rh":    r"rh\s*=\s*(" + _NUM + r")",
        "stab":  r"stab\s*=\s*(" + _NUM + r")",
        "ala":   r"ala\s*=\s*(" + _NUM + r")",
    }
    for k, pat in pats.items():
        m = re.search(pat, text)
        if m:
            meta[k] = float(m.group(1))
    if "idspl" in meta:
        meta["idspl"] = int(meta["idspl"])
        meta["ncalc"] = int(meta.get("ncalc", 1))
    return meta


def parse_predict(path):
    """解析 SLAB 输出文件.

    返回 dict:
        meta       输入参数回显
        mode       'continuous'(idspl=1/2/3) 或 'puff'(idspl=4)
        x          下风向位置表 (61,)
        cc,b,betac,zc,sig   随 x (continuous) 或随 t (puff) 的表 (61,)
        t, xc, bx, betax    随时间 t 的表 (61,)
        zplane     dict: x, timp, tcld, bbc, cv (6个横向位置, 61,)
        tav        平均时间 s
    """
    text = open(path, encoding="latin-1").read()
    lines = text.splitlines()
    meta = _parse_echo(text)

    # ---- 定位"浓度等值线参数"表 ----
    sec = None
    for i, ln in enumerate(lines):
        if "concentration contour parameters" in ln:
            sec = i
            break
    if sec is None:
        raise ValueError("输出文件中未找到 concentration contour parameters 段")

    header_i = None
    for j in range(sec, min(sec + 40, len(lines))):
        lj = lines[j]
        if ("cc(x)" in lj and "xc(t)" in lj and "betac" in lj) or \
           ("cc(t)" in lj and "xc(t)" in lj and "betac" in lj):
            header_i = j
            break
    if header_i is None:
        raise ValueError("未找到浓度等值线参数表头")

    mode = "puff" if "cc(t)" in lines[header_i] else "continuous"

    rows = []
    k = header_i + 1
    while k < len(lines) and len(rows) < 61:
        nums = re.findall(_NUM, lines[k])
        if len(nums) >= 9:
            rows.append([float(v) for v in nums[:10]])
        elif len(nums) < 9 and len(nums) > 0 and "1" in lines[k].strip():
            pass
        k += 1
    rows = np.array(rows, dtype=float)

    if mode == "continuous":
        x, cc, b, betac, zc, sig, t, xc, bx, betax = rows.T
    else:
        t, cc, b, betac, zc, sig, xc, bx, betax = rows.T[:9]
        x = t

    # ---- z 平面表(用于校验) ----
    zplane = None
    for i, ln in enumerate(lines):
        mz = re.search(r"concentration in the z =\s*(" + _NUM + r")\s*plane", ln)
        if mz and "plane" in ln:
            zpl = float(mz.group(1))
            for j in range(i, min(i + 20, len(lines))):
                if "x (m)" in lines[j]:
                    k = j + 1
                    zr = []
                    while k < len(lines) and len(zr) < 61:
                        nums = re.findall(_NUM, lines[k])
                        if len(nums) >= 10:
                            zr.append([float(v) for v in nums[:10]])
                        k += 1
                    zr = np.array(zr)
                    zplane = dict(z=zpl, x=zr[:, 0], timp=zr[:, 1],
                                  tcld=zr[:, 2], bbc=zr[:, 3], cv=zr[:, 4:10])
                    break
            break

    return dict(meta=meta, mode=mode, x=x, cc=cc, b=b, betac=betac, zc=zc,
                sig=sig, t=t, xc=xc, bx=bx, betax=betax,
                zplane=zplane, tav=meta.get("tav", 10.0))


# ------------------------------------------------------------
# 浓度场重构
# ------------------------------------------------------------
class SlabField:
    """由 SLAB 输出参数重构的浓度场 C(x, y, z, t)."""

    def __init__(self, parsed, name=""):
        self.mode = parsed["mode"]
        self.meta = parsed.get("meta", {})
        self.name = name
        self.x  = np.asarray(parsed["x"],  float)
        self.cc = np.asarray(parsed["cc"], float)
        self.b  = np.asarray(parsed["b"],  float)
        self.betac = np.asarray(parsed["betac"], float)
        self.zc = np.asarray(parsed["zc"], float)
        self.sig = np.asarray(parsed["sig"], float)
        self.t  = np.asarray(parsed["t"],  float)
        self.xc = np.asarray(parsed["xc"], float)
        self.bx = np.asarray(parsed["bx"], float)
        self.betax = np.asarray(parsed["betax"], float)
        self.zplane = parsed.get("zplane")
        self.tav = parsed.get("tav", 10.0)
        # 消除 x 表第一行的 0 值导致的插值问题
        for arr in (self.cc, self.b, self.betac, self.zc, self.sig):
            if arr[0] == 0 and len(arr) > 1:
                arr[0] = arr[1]
        # 保证 t 轴单调不减
        idx = np.argsort(self.t)
        self.t = self.t[idx]
        for arr in (self.xc, self.bx, self.betax):
            arr[:] = arr[idx]

    # ----- 基础插值 -----
    def _interp(self, tab_x, tab_y, q):
        q = np.asarray(q, dtype=float)
        return np.interp(q, tab_x, tab_y)

    # ----- 浓度计算 (体积分数) -----
    def c_volume(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            *(np.asarray(v, dtype=float) for v in (x, y, z, t)))
        sr2 = np.sqrt(2.0)
        if self.mode == "continuous":
            cc  = self._interp(self.x, self.cc, x)
            b   = self._interp(self.x, self.b,  x)
            bc  = self._interp(self.x, self.betac, x)
            zc  = self._interp(self.x, self.zc, x)
            sg  = self._interp(self.x, self.sig, x)
        else:
            cc  = self._interp(self.t, self.cc, t)
            b   = self._interp(self.t, self.b,  t)
            bc  = self._interp(self.t, self.betac, t)
            zc  = self._interp(self.t, self.zc, t)
            sg  = self._interp(self.t, self.sig, t)
        xc  = self._interp(self.t, self.xc, t)
        bx  = self._interp(self.t, self.bx, t)
        btx = self._interp(self.t, self.betax, t)

        with np.errstate(divide="ignore", invalid="ignore"):
            xa = (x - xc + bx) / (sr2 * btx)
            xb = (x - xc - bx) / (sr2 * btx)
            ya = (y + b) / (sr2 * bc)
            yb = (y - b) / (sr2 * bc)
            za = (z - zc) / (sr2 * sg)
            zb = (z + zc) / (sr2 * sg)
        c = cc * (erf(xa) - erf(xb)) * (erf(ya) - erf(yb)) * \
            (np.exp(-za * za) + np.exp(-zb * zb))
        c = np.nan_to_num(c, nan=0.0, posinf=0.0, neginf=0.0)
        # ??(??/??????)????, ???????????
        if self.mode == "continuous":
            valid = (x >= 0.0) & (x <= self.x[-1]) & (t >= 0.0) & (t <= self.t[-1])
        else:
            valid = (t >= 0.0) & (t <= self.t[-1])
        c = np.where(valid, c, 0.0)
        return np.clip(c, 0.0, 1.0)

    def ppm(self, x, y, z, t):
        return self.c_volume(x, y, z, t) * 1e6

    # ----- 固定点时间序列 -----
    def time_series(self, x, y, z, times):
        times = np.asarray(times, dtype=float)
        return self.ppm(np.full_like(times, x),
                        np.full_like(times, y),
                        np.full_like(times, z), times)

    # ----- 时间扫描网格 -----
    def time_grid(self, n=800, margin=1.0):
        t0, t1 = self.t[0], self.t[-1]
        return np.linspace(t0, t1 + margin * (t1 - t0), n)

    # ----- 下风向中心线最大浓度剖面 -----
    def max_centerline(self, z, xs, n_t=800):
        """返回 (xs, cmax_ppm, t_of_max). cmax 为各 x 处沿中心线(y=0)随时间最大值."""
        xs = np.asarray(xs, dtype=float)
        ts = self.time_grid(n_t)
        C = self.ppm(xs[:, None], 0.0, z, ts[None, :])      # (nx, nt)
        i = np.argmax(C, axis=1)
        cmax = C[np.arange(len(xs)), i]
        tmax = ts[i]
        return xs, cmax, tmax

    # ----- 危害距离: 最大浓度首次低于阈值的下风向距离 -----
    def hazard_distance(self, threshold_ppm, z=1.5, x_max=None, n_x=400,
                        n_t=800):
        x_lim = x_max if x_max else float(self.x[-1])
        xs = np.linspace(5.0, x_lim, n_x)
        _, cmax, _ = self.max_centerline(z, xs, n_t)
        above = cmax >= threshold_ppm
        if not above.any():
            return None, xs, cmax
        return float(xs[above][-1]), xs, cmax

    # ----- 到达/离开/持续时间 -----
    def exceedance_times(self, x, y, z, threshold_ppm, times=None, n=1200):
        """返回 (t_first, t_last, duration_s). C>=阈值的首/末时刻.
        未达到阈值时返回 (None, None, 0)."""
        if times is None:
            times = self.time_grid(n)
        c = self.time_series(x, y, z, times)
        hit = c >= threshold_ppm
        if not hit.any():
            return None, None, 0.0
        t_first = float(times[hit][0])
        t_last = float(times[hit][-1])
        return t_first, t_last, t_last - t_first

    # ----- 与输出文件 z 平面表校验 -----
    def validate(self):
        """与输出文件中的 z=0 平面表对比, 返回 (max_abs_err, max_rel_err)."""
        if self.zplane is None:
            return None, None
        zp = self.zplane
        c0 = self.c_volume(zp["x"], 0.0, zp.get("z", 0.0), zp["timp"])
        ref = zp["cv"][:, 0]
        mask = (zp["x"] >= 0.0) & (ref > 0)
        if not mask.any():
            return 0.0, 0.0
        abs_err = np.max(np.abs(c0[mask] - ref[mask]))
        rel_err = np.max(np.abs(c0[mask] - ref[mask]) / ref[mask])
        return float(abs_err), float(rel_err)
