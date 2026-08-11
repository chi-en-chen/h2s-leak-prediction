# -*- coding: utf-8 -*-
# ============================================================
#  h2s_sensor_analysis.py
#  油井 H2S 泄漏扩散分析与四传感器验证 (基于美国 SLAB 模型)
#  ------------------------------------------------------------
#  流程:
#   1) 按微气象要素(风速/气温/气压/风向/H2S特征)与源参数生成 SLAB 输入
#   2) 运行 Slab.exe 计算扩散场
#   3) 按风向旋转坐标, 计算四个传感器位置的浓度随时间变化
#   4) 绘制 1 / 3 / 5 / 15 分钟四个时刻的地面浓度足迹(地理坐标)
#   5) 导入传感器实测读数(CSV), 与模型预测对比验证
#
#  用法:
#    python h2s_sensor_analysis.py                          # 完整计算(无实测时仅预测)
#    python h2s_sensor_analysis.py --demo-data              # 生成示例传感器读数并验证
#    python h2s_sensor_analysis.py --sensor-data 传感器读数.csv
#    python h2s_sensor_analysis.py --no-run --outdir results_sensors_v1
# ============================================================
import argparse
import csv
import os
import sys

import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator

from slab_h2s_lib import (H2S_PROPS, SlabField, default_h2s_params,
                          parse_predict, run_slab, write_slab_input)

# ------------------------------------------------------------
# 中文字体
# ------------------------------------------------------------
def _setup_font():
    installed = {f.name for f in font_manager.fontManager.ttflist}
    for f in ("Microsoft YaHei", "Noto Sans SC", "SimHei", "DengXian"):
        if f in installed:
            plt.rcParams["font.sans-serif"] = [f]
            break
    plt.rcParams["axes.unicode_minus"] = False

_setup_font()

# ------------------------------------------------------------
# 工况配置: 油井泄漏 + 微气象 + 四个传感器 (按现场实际修改)
# ------------------------------------------------------------
CONFIG = dict(
    # -- 源项(油井泄漏) ---------------------------------------
    idspl=1,            # 1=地面蒸发池(液体积聚蒸发); 3=垂直喷口(井口高速喷出)
    qs=1.0,             # H2S 泄漏速率 kg/s
    tsd=900.0,          # 持续泄漏时间 s (10 min)
    as_=100.0,          # 液池面积 m2 (idspl=1)
    ts=213.0,           # 源温度 K (=H2S 沸点, 液池工况)
    hs=0.0,             # 源高度 m
    us=0.0,             # 水平射流速度 m/s (idspl=2)
    ws=0.0,             # 垂直喷口速度 m/s (idspl=3)
    jet_d=0.25,         # 喷口直径 m (idspl=3 垂直喷口, 换算截面积)
    jet_h=2.0,          # 喷口离地高度 m (idspl=3 垂直喷口)
    # -- 微气象要素 --------------------------------------------
    za=10.0,            # 风速测量高度 m
    ua=3.0,             # 环境风速 m/s (现场实测)
    ta=300.0,           # 环境气温 K (27 C)
    press=1013.25,      # 环境气压 hPa (用于 ppm -> mg/m3 换算)
    rh=60.0,            # 相对湿度 %
    stab=4.0,           # 大气稳定度 1..6 = A..F (D=中性)
    z0=0.1,             # 地面粗糙度 m
    wind_dir=225.0,     # 风向(气象来向): 0=北 90=东 顺时针; 225=西南风(吹向东北)
    # -- 计算域 ------------------------------------------------
    tav=10.0,           # 浓度平均时间 s
    xffm=3000.0,        # 最大下风向距离 m
    zp=1.5,             # 分析高度 m (呼吸带)
    # -- 泄露感知 ----------------------------------------------
    t_leak_start="14:00:00",   # 井内检测到泄漏的时刻(记录用, 模型 t=0 即此 刻)
    # -- 四个传感器: 钻井台正方形四角(边长60m), 油井在正方形中心 --
    #    正北S1 / 正西S2 / 正南S3 / 正东S4, 距井心 = 60/sqrt(2) ≈ 42.4 m
    #    坐标: E=东为正, N=北为正, z=高度 m
    sensors=[
        dict(name="S1", E=0.0,   N=42.4,  z=1.5, desc="钻井台正北 (距井42.4m)"),
        dict(name="S2", E=-42.4, N=0.0,   z=1.5, desc="钻井台正西 (距井42.4m)"),
        dict(name="S3", E=0.0,   N=-42.4, z=1.5, desc="钻井台正南 (距井42.4m)"),
        dict(name="S4", E=42.4,  N=0.0,   z=1.5, desc="钻井台正东 (距井42.4m)"),
    ],
    # -- 足迹图出图时刻: 1/3/5/15 分钟 -------------------------
    map_times=[60.0, 180.0, 300.0, 900.0],
    map_xmax=1200.0,    # 足迹图下风向范围 m
    map_ymax=400.0,     # 足迹图横风向范围 m
    # -- 危害阈值 ppm ------------------------------------------
    thresholds={
        100.0:  "IDLH / ERPG-3 立即威胁生命与健康",
        30.0:   "ERPG-2 不可逆或严重健康影响",
        10.0:   "短时接触限值 STEL(参考)",
        1.0:    "长期接触限值 TLV-TWA(参考)",
    },
)

# H2S 分子量 g/mol (mg/m3 换算用)
H2S_MW = H2S_PROPS["wms"] * 1000.0


# ------------------------------------------------------------
# 工具
# ------------------------------------------------------------
def wind_to_downwind(params, E, N):
    """地理坐标(E东,N北) -> 下风向坐标(x_down, y_cross).
    wind_dir 为气象来向, 下风向方位角 = wind_dir + 180."""
    beta = np.radians((params["wind_dir"] + 180.0) % 360.0)
    E, N = np.asarray(E, float), np.asarray(N, float)
    xd = E * np.sin(beta) + N * np.cos(beta)
    yc = E * np.cos(beta) - N * np.sin(beta)
    return xd, yc


def bearing(E, N):
    """油井看向传感器的地理方位角, 0=北, 90=东, 顺时针."""
    return float(np.degrees(np.arctan2(E, N)) % 360.0)


def ppm_to_mgm3(ppm, params):
    """体积浓度 ppm -> 质量浓度 mg/m3, 使用现场气压与气温."""
    R = 8.314462618           # J/(mol*K)
    P = params["press"] * 100.0   # hPa -> Pa
    return np.asarray(ppm, float) * H2S_MW * P / (R * params["ta"])


def source_geometry(params):
    """按源类型给出 SLAB 所需的源面积与源高度.

    idspl=1/4 蒸发池: 使用 as_(液池面积), 源高 0
    idspl=2   水平射流: 使用 as_(喷口截面积), 源高 hs
    idspl=3   垂直喷口: 由喷口直径 jet_d 换算截面积, 源高 jet_h
    """
    idspl = int(params.get("idspl", 1))
    if idspl == 3:
        d = float(params.get("jet_d", 0.25))
        h = float(params.get("jet_h", 2.0))
        return math.pi * (d / 2.0) ** 2, h
    return float(params.get("as_", 100.0)), float(params.get("hs", 0.0))


def build_input(params, path="h2s_sensor_input.txt"):
    p = default_h2s_params()
    as_src, hs_src = source_geometry(params)
    p["as"] = as_src
    p["hs"] = hs_src
    for k in ("idspl", "qs", "tsd", "ts", "us", "ws", "za", "ua",
              "ta", "rh", "stab", "z0", "tav", "xffm"):
        p[k] = params[k]
    p["zp1"] = params["zp"]
    gas_props = dict(H2S_PROPS)
    for k in ("wms", "tbp", "rhosl", "dhe", "cps", "cpsl"):
        if k in params and params[k] is not None:
            gas_props[k] = params[k]
    if "tbp" in params and params["tbp"] is not None:
        p["ts"] = params["tbp"]
    write_slab_input(p, path, props=gas_props)
    return path


def load_field(args, params):
    out = os.path.join(args.outdir, "predict_h2s.txt")
    if not args.no_run:
        build_input(params, os.path.join(args.outdir, "h2s_sensor_input.txt"))
        out = run_slab(exe_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "Slab.exe"),
                       input_path=os.path.join(args.outdir, "h2s_sensor_input.txt"),
                       output_path=out, workdir=args.outdir)
    else:
        out = os.path.abspath(args.input) if args.input else out
    return SlabField(parse_predict(out), name=os.path.basename(out))


# ------------------------------------------------------------
# 传感器读数读取 / 示例数据
# ------------------------------------------------------------
def sensor_names(params):
    return [s["name"] for s in params["sensors"]]


def load_obs(path, params):
    """读取传感器读数 CSV: 首列 t_s(相对泄漏检测触发的秒数), 其余列为传感器名."""
    with open(path, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError("传感器数据文件为空: " + path)
    names = sensor_names(params)
    t = np.array([float(r["t_s"]) for r in rows])
    data = {}
    for n in names:
        if n in rows[0]:
            data[n] = np.array([float(r[n]) for r in rows])
    if not data:
        raise ValueError(f"CSV 中未找到传感器列 {names}, 请按模板提供")
    return t, data


def make_demo_obs(field, params, outdir):
    """用模型预测叠加约20%随机误差, 生成示例传感器读数(用于演示验证流程)."""
    rng = np.random.default_rng(20260810)
    t = np.arange(0.0, 1200.0 + 1e-6, 10.0)
    out = {"t_s": t}
    for s in params["sensors"]:
        xd, yc = wind_to_downwind(params, s["E"], s["N"])
        pred = np.maximum(field.ppm(xd, yc, s["z"], t), 0.0)
        obs = np.where(pred > 0.05, pred * 10 ** rng.normal(0.0, 0.08, t.size), 0.0)
        out[s["name"]] = np.round(obs, 3)
    path = os.path.join(outdir, "sensor_data_demo.csv")
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["t_s"] + sensor_names(params))
        for i in range(t.size):
            w.writerow([out["t_s"][i]] + [out[n][i] for n in sensor_names(params)])
    return path


def write_template(params, outdir):
    path = os.path.join(outdir, "sensor_data_template.csv")
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["t_s", "S1", "S2", "S3", "S4"])
        w.writerow(["0", "0.0", "0.0", "0.0", "0.0"])
        w.writerow(["60", "1234.5", "56.7", "8.9", "0.1"])
        w.writerow(["120", "999.9", "88.8", "12.3", "0.2"])
    with open(os.path.join(outdir, "传感器数据格式说明.txt"), "w",
              encoding="utf-8-sig") as f:
        f.write("传感器读数文件要求:\n")
        f.write("1) CSV 格式, 首行表头, 编码 UTF-8;\n")
        f.write("2) 首列 t_s = 相对'井内检测到泄漏'时刻的秒数(模型 t=0 即泄漏检测触发);\n")
        f.write(f"3) 其余列为各传感器浓度 ppm, 列名必须为: {sensor_names(params)};\n")
        f.write("4) 未触发读数为 0 即可; 单位统一为 ppm.\n")
    return path


# ------------------------------------------------------------
# 统计计算
# ------------------------------------------------------------
def sensor_stats(field, params):
    """每个传感器的峰值/时刻/到达时刻(1/10/100 ppm)/通过时间."""
    rows = []
    for s in params["sensors"]:
        xd, yc = wind_to_downwind(params, s["E"], s["N"])
        ts = np.linspace(0.0, 1.2 * field.t[-1], 1500)
        c = field.time_series(xd, yc, s["z"], ts)
        imax = int(np.argmax(c))
        info = dict(name=s["name"], E=s["E"], N=s["N"], z=s["z"], desc=s["desc"],
                    dist=float(np.hypot(s["E"], s["N"])),
                    bearing=bearing(s["E"], s["N"]),
                    xd=float(xd), yc=float(yc),
                    peak=float(c[imax]), tpeak=float(ts[imax]))
        for thr in (100.0, 10.0, 1.0):
            t0, t1, dur = field.exceedance_times(xd, yc, s["z"], thr)
            info[f"t_first_{thr:g}"] = t0
            info[f"t_last_{thr:g}"] = t1
            info[f"dur_{thr:g}"] = dur
        rows.append(info)
    return rows


def matrix_at_times(field, params):
    """四个时间点各传感器的模型浓度 (ppm 与 mg/m3)."""
    out = []
    for s in params["sensors"]:
        xd, yc = wind_to_downwind(params, s["E"], s["N"])
        ppm = [float(field.ppm(xd, yc, s["z"], t)) for t in params["map_times"]]
        out.append((s["name"], ppm,
                    [float(ppm_to_mgm3(v, params)) for v in ppm]))
    return out


def calc_validation(field, params, obs_t, obs):
    """逐点对比: 模型预测 vs 实测, 计算验证指标."""
    names = sensor_names(params)
    pts = []
    for s in params["sensors"]:
        xd, yc = wind_to_downwind(params, s["E"], s["N"])
        pred = field.ppm(xd, yc, s["z"], obs_t)
        o = obs[s["name"]]
        m = (o > 0.05) | (pred > 0.05)
        for tt, pp, oo in zip(obs_t[m], pred[m], o[m]):
            pts.append((float(tt), s["name"], float(oo), float(pp)))
    metrics = dict(n=len(pts))
    if pts:
        t = np.array([p[0] for p in pts])
        nm = np.array([p[1] for p in pts])
        oo = np.array([p[2] for p in pts])
        pp = np.array([p[3] for p in pts])
        both = (oo > 0.05) & (pp > 0.05)
        lo, lp = np.log10(np.maximum(oo, 1e-3)), np.log10(np.maximum(pp, 1e-3))
        metrics["R_log"] = float(np.corrcoef(lo, lp)[0, 1])
        metrics["RMSE_log"] = float(np.sqrt(np.mean((lp - lo) ** 2)))
        gb = 10.0 ** np.mean(lp[both] - lo[both])
        metrics["geo_bias"] = float(gb)
        ratio = pp[both] / np.maximum(oo[both], 1e-12)
        metrics["FAC2"] = float(np.mean((ratio >= 0.5) & (ratio <= 2.0)))
        # 峰值比对(按传感器)
        pk = []
        for n in names:
            m2 = nm == n
            if m2.any() and (oo[m2] > 0.05).any():
                pk.append((n, float(oo[m2].max()), float(pp[m2].max())))
        metrics["peaks"] = pk
    return metrics, pts


# ------------------------------------------------------------
# 绘图
# ------------------------------------------------------------
def _geo_grid(params):
    xd = np.linspace(0.0, params["map_xmax"], 300)
    yc = np.linspace(-params["map_ymax"], params["map_ymax"], 220)
    Xd, Yc = np.meshgrid(xd, yc)
    beta = np.radians((params["wind_dir"] + 180.0) % 360.0)
    E = Xd * np.sin(beta) + Yc * np.cos(beta)
    N = Xd * np.cos(beta) - Yc * np.sin(beta)
    return Xd, Yc, E, N


def _wind_label(params):
    b = (params["wind_dir"] + 180.0) % 360.0
    def name(a):
        a = a % 360.0
        dirs = ["北", "东北", "东", "东南", "南", "西南", "西", "西北"]
        return dirs[int((a + 22.5) // 45) % 8]
    return (f"风向: {name(params['wind_dir'])}风({params['wind_dir']:.0f}°) "
            f"-> 下风向 {name(b)} ({b:.0f}°)")


def fig_sensor_map(params, path):
    beta = np.radians((params["wind_dir"] + 180.0) % 360.0)
    fig, ax = plt.subplots(figsize=(7.6, 6.8))
    for s in params["sensors"]:
        ax.plot(s["E"], s["N"], "o", ms=9, color="#1f77b4", mec="k", zorder=5)
        ax.annotate(f"{s['name']}\n距离{np.hypot(s['E'], s['N']):.0f}m 方位{bearing(s['E'], s['N']):.0f}°",
                    (s["E"], s["N"]), textcoords="offset points",
                    xytext=(10, 10), fontsize=8.5,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#999", alpha=0.9))
    ax.plot(0, 0, marker="*", ms=18, color="k", zorder=6)
    ax.annotate("油井(泄漏源)", (0, 0), textcoords="offset points",
                xytext=(-46, -26), fontsize=9, fontweight="bold")
    ax.quiver(0, 0, np.sin(beta), np.cos(beta), scale=1.6,
              scale_units="xy", width=0.02, color="darkred", angles="xy", zorder=4)
    ax.text(0.02, 0.94, _wind_label(params), transform=ax.transAxes,
            fontsize=10, color="darkred",
            bbox=dict(boxstyle="round,pad=0.3", fc="#fff5f0", ec="darkred", alpha=0.9))
    ax.set_xlabel("东向距离 E (m)")
    ax.set_ylabel("北向距离 N (m)")
    ax.set_title("油井 H2S 泄漏: 泄漏源与四个传感器布局")
    ax.grid(alpha=0.3)
    ax.set_aspect("equal")
    lim = 1.15 * max(np.max(np.abs([s["E"] for s in params["sensors"]] + [0])),
                     np.max(np.abs([s["N"] for s in params["sensors"]] + [0])))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    # ?????
    plot_range = 2 * lim
    raw_step = plot_range / 8.0
    mag = 10 ** int(np.floor(np.log10(max(raw_step, 1.0))))
    frac = raw_step / mag
    if frac <= 1.5:
        major = mag
    elif frac <= 3.5:
        major = 2 * mag
    elif frac <= 7.5:
        major = 5 * mag
    else:
        major = 10 * mag
    ax.xaxis.set_major_locator(MultipleLocator(major))
    ax.yaxis.set_major_locator(MultipleLocator(major))
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_footprints(field, params, path):
    Xd, Yc, E, N = _geo_grid(params)
    z = params["zp"]
    times = [t for t in params["map_times"] if t <= field.t[-1] * 1.05]
    labels = {60.0: "1 分钟", 180.0: "3 分钟", 300.0: "5 分钟", 900.0: "15 分钟"}
    n = len(times)
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 11))
    axes = np.atleast_1d(axes).ravel()
    levels = [0.1, 1.0, 10.0, 30.0, 100.0, 1000.0, 1e4, 1e5]
    cmap = plt.get_cmap("turbo").copy()
    cmap.set_under("#f0f0f0")
    norm = mcolors.LogNorm(vmin=levels[0], vmax=levels[-1])
    beta = np.radians((params["wind_dir"] + 180.0) % 360.0)
    cf = None
    for ax, t in zip(axes, times):
        C = np.maximum(field.ppm(Xd, Yc, z, t), levels[0] * 0.02)
        cf = ax.contourf(E, N, C, levels=levels, cmap=cmap, norm=norm,
                         extend="both")
        ax.contour(E, N, C, levels=levels, colors="k", linewidths=0.4, alpha=0.4)
        for s in params["sensors"]:
            ax.plot(s["E"], s["N"], "o", ms=7, color="lime", mec="k", zorder=5)
            ax.annotate(s["name"], (s["E"], s["N"]),
                        textcoords="offset points", xytext=(5, 5),
                        fontsize=8, fontweight="bold", zorder=6)
        ax.plot(0, 0, marker="*", ms=15, color="k", zorder=6)
        ax.quiver(0, 0, np.sin(beta), np.cos(beta), scale=1.6,
                  scale_units="xy", width=0.02, color="darkred",
                  angles="xy", zorder=6)
        # 中心线峰值与100ppm距离
        xf = np.linspace(0.0, params["map_xmax"], 800)
        ccen = field.ppm(xf, 0.0, z, t)
        txt = f"峰值约 {ccen.max():,.0f} ppm"
        ab = ccen >= 100.0
        if ab.any():
            txt += f"\n100ppm 达 {xf[ab][-1]:.0f} m"
        ax.text(0.985, 0.96, txt, transform=ax.transAxes, ha="right", va="top",
                fontsize=9,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#999", alpha=0.9))
        ax.text(0.015, 0.04, f"t = {t:g} s", transform=ax.transAxes,
                fontsize=9, color="#555")
        ax.set_title(labels.get(t, f"t={t:g} s"), fontsize=13, fontweight="bold")
        e0, e1 = float(E.min()), float(E.max())
        n0, n1 = float(N.min()), float(N.max())
        ax.set_xlim(e0 - 0.05 * (e1 - e0), e1 + 0.05 * (e1 - e0))
        ax.set_ylim(n0 - 0.05 * (n1 - n0), n1 + 0.05 * (n1 - n0))
        ax.set_aspect("equal")
        ax.set_xlabel("东向距离 E (m)")
        ax.set_ylabel("北向距离 N (m)")
        ax.grid(alpha=0.25)
    for ax in axes[n:]:
        ax.axis("off")
    cb = fig.colorbar(cf, ax=axes.tolist(), shrink=0.9, pad=0.02)
    cb.set_label("H2S 浓度 (ppm, 对数刻度)")
    fig.suptitle(f"油井 H2S 泄漏扩散足迹(呼吸带 z={z:g} m, {_wind_label(params)})\n"
                 f"源强 {params['qs']:g} kg/s, 持续 {params['tsd']/60:g} min, "
                 f"风速 {params['ua']:g} m/s, 稳定度 {params['stab']:g}",
                 fontsize=13)
    fig.subplots_adjust(left=0.05, right=0.96, top=0.90, bottom=0.06,
                     wspace=0.20, hspace=0.28)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_time_series(field, params, obs_t, obs, path):
    names = sensor_names(params)
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), sharex=True)
    axes = np.atleast_1d(axes).ravel()
    thr = sorted(params["thresholds"])
    tgrid = np.linspace(0.0, 1.2 * field.t[-1], 1500)
    for ax, s in zip(axes, params["sensors"]):
        xd, yc = wind_to_downwind(params, s["E"], s["N"])
        c = field.time_series(xd, yc, s["z"], tgrid)
        ax.plot(tgrid / 60.0, c, lw=1.8, color="#1f77b4",
                label="模型预测 (SLAB)")
        if obs is not None and s["name"] in obs:
            ax.scatter(obs_t / 60.0, obs[s["name"]], s=12, color="darkorange",
                       zorder=5, label="传感器实测")
        for t in thr:
            ax.axhline(t, color="#d62728", ls="--", lw=0.8, alpha=0.55)
        ax.text(1.002, 100.0, "100", transform=ax.get_yaxis_transform(),
                color="#d62728", fontsize=7, va="center")
        ax.text(1.002, 10.0, "10", transform=ax.get_yaxis_transform(),
                color="#d62728", fontsize=7, va="center")
        ax.set_yscale("log")
        ax.set_ylim(1e-3, 2e6)
        ax.set_title(f"{s['name']}  {s['desc']}\n"
                     f"下风向 {xd:.0f} m, 横风 {yc:+.0f} m",
                     fontsize=10.5)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="upper right")
        ax.set_ylabel("H2S 浓度 (ppm)")
    for ax in axes[len(params["sensors"]):]:
        ax.axis("off")
    for ax in axes:
        ax.set_xlabel("时间 (min, 相对井内泄漏检测触发)")
    fig.suptitle(f"四传感器 H2S 浓度时间序列 (源强 {params['qs']:g} kg/s, "
                 f"持续 {params['tsd']/60:g} min, 风速 {params['ua']:g} m/s, "
                 f"{_wind_label(params)})", fontsize=13)
    fig.tight_layout(rect=(0, 0, 0.98, 0.94))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_validation(metrics, pts, names, path):
    fig, ax = plt.subplots(figsize=(8, 7.2))
    if len(pts):
        colors = plt.cm.Set2(np.linspace(0, 1, len(names)))
        for i, n in enumerate(names):
            m = np.array([p[1] == n for p in pts])
            if m.any():
                po = np.array([p[2] for p in pts])[m]
                pp = np.array([p[3] for p in pts])[m]
                ax.scatter(np.maximum(po, 0.01), np.maximum(pp, 0.01),
                           s=22, color=colors[i], label=n, alpha=0.85)
    lo = np.logspace(-2, 6, 100)
    ax.plot(lo, lo, "k-", lw=1.2)
    ax.plot(lo, 2 * lo, "k--", lw=0.8, alpha=0.5)
    ax.plot(lo, 0.5 * lo, "k--", lw=0.8, alpha=0.5)
    ax.text(0.98, 0.06, "1:1", transform=ax.transAxes, ha="right", fontsize=8)
    ax.text(0.98, 0.115, "FAC2 边界", transform=ax.transAxes, ha="right", fontsize=7, color="gray")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.01, 1e6); ax.set_ylim(0.01, 1e6)
    ax.set_xlabel("传感器实测浓度 (ppm)")
    ax.set_ylabel("模型预测浓度 (ppm)")
    ax.set_title("SLAB 模型验证: 预测 vs 实测")
    if len(pts):
        txt = (f"有效点数 n = {metrics['n']}\n"
               f"对数空间相关系数 R = {metrics['R_log']:.3f}\n"
               f"对数均方根误差 = {metrics['RMSE_log']:.3f}\n"
               f"几何平均偏差 = {metrics['geo_bias']:.2f}\n"
               f"FAC2 (0.5~2倍) = {metrics['FAC2']*100:.0f}%")
        ax.text(0.04, 0.96, txt, transform=ax.transAxes, va="top", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", fc="#f4f4f4", ec="#888", alpha=0.95))
    else:
        ax.text(0.5, 0.5, "未提供传感器读数, 无验证数据", transform=ax.transAxes,
                ha="center", fontsize=12, color="#888")
    if len(pts):
        ax.legend(loc="lower right", fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ------------------------------------------------------------
# 汇总输出
# ------------------------------------------------------------
def write_summary(field, params, stats, mat, metrics, pts, obs_src, args, path):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["油井 H2S 泄漏扩散分析与传感器验证 (SLAB 模型)"])
        w.writerow([])
        w.writerow(["一、工况与微气象参数"])
        w.writerow(["泄漏速率 qs (kg/s)", params["qs"]])
        w.writerow(["持续泄漏 tsd (s)", params["tsd"]])
        w.writerow(["液池面积 as (m2)", params["as_"]])
        w.writerow(["源温度 ts (K)", params["ts"]])
        w.writerow(["源类型 idspl", params["idspl"]])
        w.writerow(["风速 ua (m/s)", params["ua"]])
        w.writerow(["气温 ta (K)", params["ta"]])
        w.writerow(["气压 press (hPa)", params["press"]])
        w.writerow(["风向 wind_dir (来向, 度)", params["wind_dir"]])
        w.writerow(["大气稳定度 stab", params["stab"]])
        w.writerow(["相对湿度 rh (%)", params["rh"]])
        w.writerow(["地面粗糙度 z0 (m)", params["z0"]])
        w.writerow(["分析高度 zp (m)", params["zp"]])
        w.writerow(["平均时间 tav (s)", params["tav"]])
        w.writerow(["井内检测到泄漏时刻", params["t_leak_start"]])
        w.writerow([])
        w.writerow(["说明: SLAB 采用标准大气压计算体积浓度(ppm); "
                    "气压/气温用于 ppm 换算 mg/m3."])
        w.writerow([])
        w.writerow(["二、传感器布局(油井为原点, E东/N北)"])
        w.writerow(["名称", "E (m)", "N (m)", "z (m)", "方位(度)", "距离 (m)",
                    "下风向 x (m)", "横风 y (m)", "位置说明"])
        for s in stats:
            w.writerow([s["name"], s["E"], s["N"], s["z"],
                        f"{s['bearing']:.0f}", f"{s['dist']:.0f}",
                        f"{s['xd']:.0f}", f"{s['yc']:.0f}", s["desc"]])
        w.writerow([])
        w.writerow(["三、传感器预测统计(模型)"])
        w.writerow(["名称", "峰值 ppm", "峰值 mg/m3", "峰值时刻 (min)",
                    "到达1ppm (s)", "离开1ppm (s)",
                    "到达10ppm (s)", "离开10ppm (s)",
                    "到达100ppm (s)", "离开100ppm (s)"])
        for s in stats:
            def g(k):
                v = s[k]
                return f"{v:.0f}" if v else "-"
            w.writerow([s["name"], f"{s['peak']:.1f}",
                        f"{ppm_to_mgm3(s['peak'], params):.1f}",
                        f"{s['tpeak']/60:.1f}",
                        g("t_first_1"), g("t_last_1"),
                        g("t_first_10"), g("t_last_10"),
                        g("t_first_100"), g("t_last_100")])
        w.writerow([])
        w.writerow(["四、1/3/5/15 分钟时刻浓度预测 (ppm / mg/m3)"])
        w.writerow(["传感器"] + [f"{t/60:g}min_ppm" for t in params["map_times"]]
                   + [f"{t/60:g}min_mgm3" for t in params["map_times"]])
        for name, ppm, mg in mat:
            w.writerow([name] + [f"{v:.3f}" for v in ppm]
                       + [f"{v:.3f}" for v in mg])
        w.writerow([])
        w.writerow(["五、模型验证(预测 vs 传感器实测)"])
        w.writerow(["实测数据来源", obs_src])
        if len(pts):
            w.writerow(["有效对比点数", metrics["n"]])
            w.writerow(["对数空间相关系数 R", f"{metrics['R_log']:.3f}"])
            w.writerow(["对数空间 RMSE", f"{metrics['RMSE_log']:.3f}"])
            w.writerow(["几何平均偏差", f"{metrics['geo_bias']:.2f}"])
            w.writerow(["FAC2 (预测/实测 0.5~2 比例)", f"{metrics['FAC2']*100:.0f}%"])
            w.writerow([])
            w.writerow(["峰值比对", "实测峰值 ppm", "预测峰值 ppm",
                        "预测/实测"])
            for n, po, pp in metrics["peaks"]:
                w.writerow([n, f"{po:.1f}", f"{pp:.1f}",
                            f"{pp/max(po,1e-9):.2f}"])
            w.writerow([])
            w.writerow(["逐点对照", "时间 (s)", "传感器", "实测 ppm",
                        "预测 ppm", "预测/实测"])
            for t, n, oo, pp in pts:
                w.writerow(["", f"{t:.0f}", n, f"{oo:.3f}", f"{pp:.3f}",
                            f"{pp/max(oo,1e-9):.2f}"])
        else:
            w.writerow(["未提供传感器读数文件, 仅模型预测. "
                        "可用 --demo-data 生成示例数据演示验证流程."])
        w.writerow([])
        w.writerow(["六、图件清单"])
        for fn in ("fig1_传感器布局.png", "fig2_浓度足迹_1_3_5_15分钟.png",
                   "fig3_四传感器时间序列.png", "fig4_模型验证.png"):
            w.writerow([fn, os.path.join(os.path.abspath(args.outdir), fn)])


def print_summary(field, params, stats, metrics, obs_src, args):
    print("=" * 70)
    print("  油井 H2S 泄漏扩散分析与传感器验证 (SLAB)")
    print("=" * 70)
    wd = (params["wind_dir"] + 180.0) % 360.0
    print(f"  源强 {params['qs']:g} kg/s, 持续 {params['tsd']/60:g} min, "
          f"风速 {params['ua']:g} m/s, 气温 {params['ta']-273.15:.0f} C, "
          f"气压 {params['press']:.0f} hPa, 稳定度 {params['stab']:g}")
    print(f"  风向 {params['wind_dir']:.0f} 度(来向) -> 下风向 {wd:.0f} 度")
    ab, rel = field.validate()
    if ab is not None:
        print(f"  浓度重构与 SLAB z=0 表对比: 绝对误差 {ab:.2e}, 相对 {rel*100:.1f}%")
    print("-" * 70)
    print("  传感器预测峰值 (ppm):")
    for s in stats:
        print(f"    {s['name']}: {s['peak']:10.1f} ppm  "
              f"(t={s['tpeak']/60:5.1f} min, 下风向 {s['xd']:.0f} m, "
              f"横风 {s['yc']:+.0f} m)")
    print("-" * 70)
    print("  1/3/5/15 min 时刻中心线 100ppm 危害距离:")
    for t in params["map_times"]:
        xs = np.linspace(5.0, params["map_xmax"], 600)
        ct = field.ppm(xs, 0.0, params["zp"], t)
        ab100 = ct >= 100.0
        d = xs[ab100][-1] if ab100.any() else None
        print(f"    t={t/60:g} min: " + (f"{d:.0f} m" if d else "未达到100ppm"))
    if len(metrics.get("peaks", [])):
        print("-" * 70)
        print("  模型验证 (示例/实测数据源: " + obs_src + "):")
        print(f"    n={metrics['n']}, R(log)={metrics['R_log']:.3f}, "
              f"FAC2={metrics['FAC2']*100:.0f}%")
        for n, po, pp in metrics["peaks"]:
            print(f"    {n}: 实测峰值 {po:.1f} ppm vs 预测 {pp:.1f} ppm "
                  f"(比 {pp/max(po,1e-9):.2f})")
    print("-" * 70)
    print(f"  输出目录: {args.outdir}")
    print("=" * 70)


# ------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="油井 H2S 泄漏扩散分析与传感器验证 (SLAB)")
    ap.add_argument("--no-run", action="store_true",
                    help="不重跑 Slab.exe, 解析已有输出文件")
    ap.add_argument("--input", default=None,
                    help="已有 SLAB 输出文件路径 (配合 --no-run)")
    ap.add_argument("--outdir", default="results_sensors_v1",
                    help="输出目录")
    ap.add_argument("--sensor-data", default=None,
                    help="传感器实测读数 CSV 路径 (首列 t_s, 其余列为传感器名)")
    ap.add_argument("--demo-data", action="store_true",
                    help="用模型预测生成示例传感器读数, 演示验证流程")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    params = dict(CONFIG)

    print(">> 运行 SLAB 计算扩散场 ...")
    field = load_field(args, params)

    # 传感器读数
    obs_t, obs = None, None
    obs_src = "未提供"
    if args.sensor_data:
        obs_t, obs = load_obs(args.sensor_data, params)
        obs_src = args.sensor_data
    elif args.demo_data:
        obs_t, obs = load_obs(make_demo_obs(field, params, args.outdir), params)
        obs_src = "sensor_data_demo.csv (示例数据, 请替换为真实读数)"
        write_template(params, args.outdir)
    else:
        write_template(params, args.outdir)

    print(">> 计算传感器统计 ...")
    stats = sensor_stats(field, params)
    mat = matrix_at_times(field, params)
    metrics, pts = calc_validation(field, params, obs_t, obs) if obs is not None \
        else (dict(n=0, peaks=[]), np.empty((0, 4)))

    print(">> 绘图 ...")
    fig_sensor_map(params, os.path.join(args.outdir, "fig1_传感器布局.png"))
    fig_footprints(field, params,
                   os.path.join(args.outdir, "fig2_浓度足迹_1_3_5_15分钟.png"))
    fig_time_series(field, params, obs_t, obs,
                    os.path.join(args.outdir, "fig3_四传感器时间序列.png"))
    fig_validation(metrics, pts, sensor_names(params),
                   os.path.join(args.outdir, "fig4_模型验证.png"))

    write_summary(field, params, stats, mat, metrics, pts, obs_src, args,
                  os.path.join(args.outdir, "h2s_sensor_summary.csv"))
    print_summary(field, params, stats, metrics, obs_src, args)


if __name__ == "__main__":
    main()