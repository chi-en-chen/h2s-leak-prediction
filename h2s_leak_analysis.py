# -*- coding: utf-8 -*-
# ============================================================
#  h2s_leak_analysis.py
#  检测井 H2S 泄露后随时间扩散情况分析
#  ------------------------------------------------------------
#  基于当前目录下的 SLAB 重气扩散模型(Slab.exe + 输出解析):
#   1) 生成 H2S 工况输入文件
#   2) 运行 Slab.exe 得到预测结果
#   3) 重构浓度场 C(x,y,z,t), 分析:
#       - 各监测点浓度随时间变化曲线
#       - 下风向中心线最大浓度剖面与危害距离
#       - 多个时刻的地面浓度足迹(扩散过程)
#       - 各距离到达时刻 / 云团通过持续时间
#   4) 输出 PNG 图 + CSV 汇总
#  ------------------------------------------------------------
#  用法:
#    python h2s_leak_analysis.py                 # 生成工况并运行 SLAB 后分析
#    python h2s_leak_analysis.py --no-run        # 不重跑 SLAB, 解析现有 predict_h2s.txt
#    python h2s_leak_analysis.py --input 某输出文件 --no-run
#    python h2s_leak_analysis.py --outdir results
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
# 工况配置(检测井 H2S 泄露, 按实际情况修改)
# ------------------------------------------------------------
CONFIG = dict(
    # -- 源项 --------------------------------------------------
    idspl=1,            # 1=地面蒸发池; 2=水平射流; 3=垂直喷口; 4=瞬时蒸发池
    qs=1.0,             # H2S 泄露速率 kg/s
    tsd=600.0,          # 连续泄露持续时间 s (10 min)
    as_=100.0,          # 液池面积 m2
    ts=213.0,           # 源温度 K (=H2S 沸点, 液池工况)
    hs=0.0,             # 源高度 m (井口/池面)
    us=0.0,             # 水平射流速度 m/s (idspl=2)
    ws=0.0,             # 垂直喷口速度 m/s (idspl=3)
    # -- 气象 --------------------------------------------------
    za=10.0,            # 风速测量高度 m
    ua=2.0,             # 环境风速 m/s
    ta=300.0,           # 环境温度 K
    rh=50.0,            # 相对湿度 %
    stab=4.0,           # 大气稳定度 1..6 = A..F (D=中性)
    z0=0.1,             # 地面粗糙度 m
    # -- 计算域 ------------------------------------------------
    tav=10.0,           # 浓度平均时间 s
    xffm=3000.0,        # 最大下风向距离 m
    zp=1.5,             # 浓度分析高度 m (呼吸带)
    # -- 监测点 (x, y, z) m ------------------------------------
    monitors=[(50.0, 0.0, 1.5), (100.0, 0.0, 1.5), (200.0, 0.0, 1.5),
              (500.0, 0.0, 1.5), (1000.0, 0.0, 1.5), (2000.0, 0.0, 1.5)],
    # -- 足迹图出图时刻 s --------------------------------------
    map_times=[60.0, 180.0, 360.0, 600.0, 900.0, 1200.0],
    map_xmax=2000.0,    # 足迹图下风向范围 m
    map_ymax=350.0,     # 足迹图横风向范围 m
    # -- 危害阈值 ppm ------------------------------------------
    thresholds={
        100.0:  "IDLH / ERPG-3 立即威胁生命与健康",
        30.0:   "ERPG-2 不可逆或严重健康影响",
        10.0:   "短时接触限值 STEL(参考)",
        1.0:    "长期接触限值 TLV-TWA(参考)",
        0.1:    "ERPG-1 轻微可逆影响",
    },
)

# 垂直喷口(井口泄漏)备选工况: idspl=3 时使用
JET_PRESET = dict(
    idspl=3, qs=1.0, tsd=600.0, as_=0.008, ts=300.0, hs=2.0,
    ws=30.0, us=0.0,
)

# ------------------------------------------------------------
def build_input(params, path="h2s_input.txt"):
    p = default_h2s_params()
    if int(params.get("idspl", 1)) == 3:
        # 垂直喷口: 由喷口直径换算截面积, 源高取喷口离地高度
        d = float(params.get("jet_d", 0.25))
        p["as"] = math.pi * (d / 2.0) ** 2
        p["hs"] = float(params.get("jet_h", 2.0))
    else:
        p["as"] = params["as_"]
        p["hs"] = params.get("hs", 0.0)
    for k in ("idspl", "qs", "tsd", "ts", "us", "ws", "za", "ua",
              "ta", "rh", "stab", "z0", "tav", "xffm"):
        p[k] = params[k]
    p["zp1"] = params["zp"]
    write_slab_input(p, path, props=H2S_PROPS)
    return path


def load_field(args, params):
    """生成输入 + 运行 SLAB(或解析已有文件) -> SlabField"""
    out = os.path.join(args.outdir, "predict_h2s.txt")
    if not args.no_run:
        build_input(params, os.path.join(args.outdir, "h2s_input.txt"))
        out = run_slab(exe_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "Slab.exe"),
                       input_path=os.path.join(args.outdir, "h2s_input.txt"),
                       output_path=out,
                       workdir=args.outdir)
    else:
        out = os.path.abspath(args.input) if args.input else out
    field = SlabField(parse_predict(out), name=os.path.basename(out))
    return field


# ------------------------------------------------------------
def analyze_time_series(field, params):
    z = params["zp"]
    t1 = field.t[-1]
    times = np.linspace(0.0, 1.2 * t1, 1500)
    rows = []
    for (x, y, hz) in params["monitors"]:
        c = field.time_series(x, y, hz, times)
        imax = int(np.argmax(c))
        rows.append(dict(x=x, y=y, z=hz, cmax_ppm=float(c[imax]),
                         tmax_s=float(times[imax])))
    return times, rows


def analyze_hazard(field, params):
    z = params["zp"]
    xs = np.logspace(np.log10(5.0), np.log10(field.x[-1]), 300)
    _, cmax, tmax = field.max_centerline(z, xs, n_t=1000)
    hz = {}
    for thr in sorted(params["thresholds"]):
        above = cmax >= thr
        d = float(xs[above][-1]) if above.any() else None
        tm = float(tmax[above][-1]) if above.any() else None
        hz[thr] = (d, tm)
    return xs, cmax, tmax, hz


def analyze_exceedance(field, params):
    """到达时刻/离开时刻/持续时间 随下风向距离"""
    z = params["zp"]
    xs = np.linspace(20.0, min(2000.0, field.x[-1] - 10), 80)
    ts = np.linspace(0.0, 1.2 * field.t[-1], 2500)
    out = {}
    for thr in (10.0, 100.0):
        if thr not in params["thresholds"]:
            continue
        C = field.ppm(xs[:, None], 0.0, z, ts[None, :])
        hit = C >= thr
        has = hit.any(axis=1)
        idx_first = np.argmax(hit, axis=1)
        idx_last = hit.shape[1] - 1 - np.argmax(hit[:, ::-1], axis=1)
        out[thr] = dict(xs=xs,
                        t_first=np.where(has, ts[idx_first], np.nan),
                        t_last=np.where(has, ts[idx_last], np.nan))
    return out


# ------------------------------------------------------------
# 绘图
# ------------------------------------------------------------
def fig_time_series(field, times, rows, params, path):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for r in rows:
        x, y, hz = r["x"], r["y"], r["z"]
        c = field.time_series(x, y, hz, times)
        ax.plot(times / 60.0, c, lw=1.6, label=f"x={x:.0f} m (y=0, z={hz:.0f} m)")
    thr = sorted(params["thresholds"])
    colors = plt.cm.Reds(np.linspace(0.4, 0.95, len(thr)))
    for k, t in enumerate(thr):
        ax.axhline(t, color=colors[k], ls="--", lw=1.0)
        ax.text(1.02, t, f"{t:g} ppm", color=colors[k], va="center",
                transform=ax.get_yaxis_transform(), fontsize=8)
    ax.set_xlabel("时间 (min)")
    ax.set_ylabel("H2S 体积浓度 (ppm)")
    ax.set_yscale("log")
    ax.set_ylim(0.01, 2e6)
    ax.set_title("检测井 H2S 泄露后各监测点浓度随时间变化\n"
                 f"(源强 {params['qs']:g} kg/s, 持续 {params['tsd']/60:g} min, "
                 f"风速 {params['ua']:g} m/s, 稳定度 {params['stab']:g})")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_footprints(field, params, path):
    xmax, ymax = params["map_xmax"], params["map_ymax"]
    nx, ny = 320, 260
    xs = np.linspace(0.0, xmax, nx)
    ys = np.linspace(-ymax, ymax, ny)
    X, Y = np.meshgrid(xs, ys)
    z = params["zp"]
    times = [t for t in params["map_times"] if t <= field.t[-1] * 1.2]
    ncol = 3
    nrow = int(np.ceil(len(times) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.4 * ncol, 3.6 * nrow))
    axes = np.atleast_1d(axes).ravel()

    levels = [0.1, 1.0, 10.0, 30.0, 100.0, 1000.0, 1e4, 1e5]
    cmap = plt.get_cmap("turbo").copy()
    cmap.set_under("#f0f0f0")
    norm = mcolors.LogNorm(vmin=levels[0], vmax=levels[-1])

    for ax, t in zip(axes, times):
        C = np.maximum(field.ppm(X, Y, z, t), levels[0] * 0.01)
        cf = ax.contourf(X, Y, C, levels=levels, cmap=cmap, norm=norm,
                         extend="both")
        ax.contour(X, Y, C, levels=levels, colors="k", linewidths=0.5,
                   alpha=0.45)
        xf = np.linspace(0.0, xmax, 1000)
        ccen = field.ppm(xf, 0.0, z, t)
        txt = f"峰值约 {ccen.max():,.0f} ppm"
        above = ccen >= 100.0
        if above.any():
            txt += f"\n100 ppm 达 {xf[above][-1]:.0f} m"
        ax.text(0.985, 0.96, txt, transform=ax.transAxes, ha="right",
                va="top", fontsize=8.5,
                bbox=dict(boxstyle="round,pad=0.25", fc="white",
                          ec="#999999", alpha=0.85))
        ax.plot([0], [0], marker="*", ms=12, color="k", zorder=5)
        ax.set_title(f"t = {t:g} s")
        ax.set_xlim(0, xmax)
        ax.set_ylim(-ymax, ymax)
        ax.set_aspect("equal")
        ax.set_xlabel("下风向距离 x (m)")
        ax.set_ylabel("横向距离 y (m)")
    for ax in axes[len(times):]:
        ax.axis("off")
    cb = fig.colorbar(cf, ax=axes.tolist(), shrink=0.85, ticks=levels)
    cb.set_label("H2S 浓度 (ppm, 对数)")
    cb.set_ticklabels(["0.1", "1", "10", "30", "100", "1e3", "1e4", "1e5"])
    fig.suptitle("检测井 H2S 泄露后地面(呼吸带)浓度足迹随时间演变", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_max_centerline(field, xs, cmax, hz, params, path):
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.plot(xs, cmax, lw=2, color="tab:red", label="最大浓度(随时间)剖面")
    thr = sorted(params["thresholds"])
    for t in thr:
        d, _ = hz[t]
        if d is None:
            continue
        ax.axvline(d, color="grey", ls=":", lw=1.0)
        ax.axhline(t, color="grey", ls=":", lw=1.0)
        ax.annotate(f"{t:g} ppm -> {d:.0f} m",
                    xy=(d, t), xytext=(d * 0.85, t * 1.8),
                    fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(xs[0], xs[-1])
    ax.set_ylim(0.05, 2e6)
    ax.set_xlabel("下风向距离 x (m)")
    ax.set_ylabel("呼吸带高度最大 H2S 浓度 (ppm)")
    ax.set_title("下风向中心线最大 H2S 浓度与危害距离\n"
                 f"(z={params['zp']:g} m, 源强 {params['qs']:g} kg/s, "
                 f"风速 {params['ua']:g} m/s)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_exceedance(ex, params, path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, key, ylab, ttl in [
        (axes[0], "t_first", "到达时刻 (min)", "浓度首次达到阈值的时刻"),
        (axes[1], "t_last", "离开时刻 (min)", "浓度最后低于阈值的时刻")]:
        for thr in (10.0, 100.0):
            d = ex[thr]
            m = ~np.isnan(d[key])
            ax.plot(d["xs"][m] / 1000.0, d[key][m] / 60.0,
                    marker="o", ms=3, lw=1.2, label=f"{thr:g} ppm")
        ax.set_xlabel("下风向距离 x (km)")
        ax.set_ylabel(ylab)
        ax.set_title(ttl)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
    fig.suptitle("H2S 云团沿下风向的到达/通过时间 (呼吸带高度)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ------------------------------------------------------------
def write_summary(field, params, hz, rows, args, path):
    thr = sorted(params["thresholds"])
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["检测井 H2S 泄露扩散分析汇总 (SLAB 模型)"])
        w.writerow([])
        w.writerow(["一、工况参数", "", "", ""])
        w.writerow(["泄露速率 qs (kg/s)", params["qs"]])
        w.writerow(["泄露持续时间 tsd (s)", params["tsd"]])
        w.writerow(["液池面积 as (m2)", params["as_"]])
        w.writerow(["源高度 hs (m)", params["hs"]])
        w.writerow(["环境风速 ua (m/s)", params["ua"]])
        w.writerow(["大气稳定度", params["stab"]])
        w.writerow(["环境温度 ta (K)", params["ta"]])
        w.writerow(["平均时间 tav (s)", params["tav"]])
        w.writerow([])
        w.writerow(["二、危害距离(呼吸带高度, 下风向中心线最大浓度)", "", "", ""])
        w.writerow(["阈值 ppm", "危害距离 m", "达到时刻 s", "说明"])
        for t in thr:
            d, tm = hz[t]
            w.writerow([t, f"{d:.0f}" if d else "未达到",
                        f"{tm:.0f}" if tm else "-",
                        params["thresholds"][t]])
        w.writerow([])
        w.writerow(["三、监测点浓度统计", "", "", ""])
        w.writerow(["x m", "y m", "z m", "最大浓度 ppm", "出现时刻 s"])
        for r in rows:
            w.writerow([r["x"], r["y"], r["z"], f"{r['cmax_ppm']:.1f}",
                        f"{r['tmax_s']:.0f}"])
        w.writerow([])
        w.writerow(["四、校验: 与 SLAB 输出 z=0 表对比最大误差", "", "", ""])
        ab, rel = field.validate()
        if ab is not None:
            w.writerow(["绝对误差", f"{ab:.3e}"])
            w.writerow(["相对误差", f"{rel*100:.2f}%"])
        else:
            w.writerow(["无 z 平面表, 未校验"])


def print_summary(field, params, hz, rows, args):
    print("=" * 64)
    print("  检测井 H2S 泄露扩散分析结果 (SLAB 模型)")
    print("=" * 64)
    print(f"  源强 qs = {params['qs']:g} kg/s, 持续 {params['tsd']/60:g} min, "
          f"风速 {params['ua']:g} m/s, 稳定度 {params['stab']:g}")
    ab, rel = field.validate()
    if ab is not None:
        print(f"  与 SLAB 输出 z=0 表对比: 最大绝对误差 {ab:.2e}, "
              f"最大相对误差 {rel*100:.1f}%")
    print("-" * 64)
    print("  危害距离 (下风向中心线, 呼吸带高度):")
    for t in sorted(params["thresholds"]):
        d, tm = hz[t]
        if d:
            print(f"    {t:6g} ppm: {d:7.0f} m   (约 {tm/60:.1f} min 后达到)")
        else:
            print(f"    {t:6g} ppm: 未达到")
    print("-" * 64)
    print("  监测点 (x, y=0, z=1.5 m) 最大浓度:")
    for r in rows:
        print(f"    x={r['x']:6.0f} m: {r['cmax_ppm']:12.1f} ppm  "
              f"(t={r['tmax_s']/60:.1f} min)")
    print("-" * 64)
    print(f"  图件已保存: {args.outdir}")
    print("=" * 64)


# ------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="检测井 H2S 泄露扩散分析 (SLAB)")
    ap.add_argument("--no-run", action="store_true",
                    help="不重跑 Slab.exe, 解析已有输出文件")
    ap.add_argument("--input", default=None,
                    help="已有 SLAB 输出文件路径 (配合 --no-run)")
    ap.add_argument("--outdir", default="results",
                    help="输出目录 (默认 results)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    params = dict(CONFIG)

    print(">> 准备工况并运行 SLAB ...")
    field = load_field(args, params)

    print(">> 计算监测点时间序列 ...")
    times, rows = analyze_time_series(field, params)

    print(">> 计算危害距离 ...")
    xs, cmax, tmax, hz = analyze_hazard(field, params)

    print(">> 计算到达/通过时间 ...")
    ex = analyze_exceedance(field, params)

    print(">> 绘图 ...")
    fig_time_series(field, times, rows, params,
                    os.path.join(args.outdir, "fig1_监测点浓度时间序列.png"))
    fig_footprints(field, params,
                   os.path.join(args.outdir, "fig2_浓度足迹演变.png"))
    fig_max_centerline(field, xs, cmax, hz, params,
                       os.path.join(args.outdir, "fig3_中心线最大浓度与危害距离.png"))
    fig_exceedance(ex, params,
                   os.path.join(args.outdir, "fig4_到达与通过时间.png"))

    write_summary(field, params, hz, rows, args,
                  os.path.join(args.outdir, "h2s_summary.csv"))
    print_summary(field, params, hz, rows, args)


if __name__ == "__main__":
    main()
