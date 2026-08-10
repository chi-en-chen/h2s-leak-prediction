# -*- coding: utf-8 -*-
"""
油井 H2S 泄漏扩散 - 核心可视化元素 (SLAB 重气模型输出)
 1) 传感器点位: 4 个传感器位于油井正北/正西/正南/正东, 距井口 60 m
    统一图标 + 编号 + 实时浓度读数(模型示意值), 绿=正常 红=告警
 2) H2S 重气扩散云团: 半透明浓度分级渲染
    低浓度浅黄(1-10ppm) / 中浓度橙(10-100ppm) / 高浓度深红(>100ppm)
    叠加等值线: 1ppm 警戒区边界 / 10ppm 危险区边界 / 100ppm 毒性阈值边界
 3) 应急风险分区: 红色实线=核心应急隔离区(≥10ppm外扩) / 黄色虚线=警戒疏散区(≥1ppm外扩)
    带箭头人员疏散路线: 从危险区指向外围安全区域
"""
import argparse
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from matplotlib.transforms import Affine2D

import h2s_sensor_analysis as base
import h2s_time_slices as tsl
from slab_h2s_lib import parse_predict

# ------------------------------------------------------------
# 传感器布局: 正北/正西/正南/正东, 距井口 60 m (最新要求)
# ------------------------------------------------------------
SENSOR_LAYOUT = [
    dict(name="S1", E=0.0,  N=60.0,  desc="正北"),
    dict(name="S2", E=-60.0, N=0.0,  desc="正西"),
    dict(name="S3", E=0.0,  N=-60.0, desc="正南"),
    dict(name="S4", E=60.0, N=0.0,   desc="正东"),
]
ALARM_PPM = 10.0      # 传感器告警阈值 ppm
ZP = 1.5              # 传感器/呼吸带高度 m

# 浓度分级与等值线 (ppm)
LEVELS    = [1.0, 10.0, 100.0, 1000.0]
ZONE_META = [
    (1.0,   1.0,  "#E6B800", "1 ppm 警戒区边界"),
    (10.0,  10.0, "#FF8C1A", "10 ppm 危险区边界"),
    (100.0, 100.0, "#C00000", "100 ppm 毒性阈值边界"),
]
CMAP_CLOUD = ListedColormap(["#FFF7B0", "#FFD84D", "#FF8C1A", "#B30000"])
ALPHA_CLOUD = 0.45

# 分区外扩 (m)
CORE_EXPAND = 60.0
WARN_EXPAND = 100.0


def load_field(outdir, input_path=None):
    import argparse as _ap
    args = _ap.Namespace(outdir=outdir, no_run=True, input=input_path)
    return tsl.load_field(args, dict(base.CONFIG))


def cloud_geo_grid(params, nx=320, ny=301, xmax=3400.0, ymax=420.0):
    """下风向坐标网格上的地面浓度 -> 地理坐标(E,N)网格."""
    beta = np.radians((params["wind_dir"] + 180.0) % 360.0)
    xd = np.linspace(0.0, xmax, nx)
    yc = np.linspace(-ymax, ymax, ny)
    C = np.empty((nx, ny))
    field = None
    return field  # placeholder (实现在 draw_panel 内)

def sensor_readings(field, params, t):
    """传感器位置呼吸带浓度与状态."""
    out = []
    for s in SENSOR_LAYOUT:
        xd, yc = base.wind_to_downwind(params, s["E"], s["N"])
        ppm = float(field.ppm(np.array([xd]), np.array([yc]), ZP, t)[0])
        out.append(dict(name=s["name"], desc=s["desc"], E=s["E"], N=s["N"],
                        ppm=ppm, alarm=ppm >= ALARM_PPM))
    return out


def bbox_rot(field, params, t, thr, expand):
    """浓度>=thr 云团在旋转坐标下的外接矩形(带外扩), 返回 (cx, cy, w, h)."""
    beta = np.radians((params["wind_dir"] + 180.0) % 360.0)
    nx, ny, xmax, ymax = 360, 361, 3600.0, 500.0
    xd = np.linspace(0.0, xmax, nx)
    yc = np.linspace(-ymax, ymax, ny)
    C = field.ppm(xd[:, None], yc[None, :], ZP, t)
    mask = C >= thr
    if not mask.any():
        return 0.0, 0.0, 2 * expand, 2 * expand
    xs, ys = np.where(mask)
    x0, x1 = xd[xs.min()], xd[xs.max()]
    y0, y1 = yc[ys.min()], yc[ys.max()]
    cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0
    return cx, cy, (x1 - x0) + 2 * expand, (y1 - y0) + 2 * expand


def draw_panel(ax, field, params, t, beta_deg):
    """绘制单个时刻的核心可视化面板 (地理坐标 E-N)."""
    beta = np.radians(beta_deg)
    nx, ny, xmax, ymax = 320, 301, 3600.0, 500.0
    xd = np.linspace(0.0, xmax, nx)
    yc = np.linspace(-ymax, ymax, ny)
    C = field.ppm(xd[:, None], yc[None, :], ZP, t)      # (nx, ny) 旋转坐标
    Xe = xd[:, None] * np.sin(beta) + yc[None, :] * np.cos(beta)  # (nx, ny)
    Ne = xd[:, None] * np.cos(beta) - yc[None, :] * np.sin(beta)

    # ---- 2) 重气云团: 半透明浓度分级 ----
    Cm = np.ma.masked_less(C.T, LEVELS[0])
    pcm = ax.pcolormesh(Xe.T, Ne.T, Cm, cmap=CMAP_CLOUD, alpha=ALPHA_CLOUD,
                        norm=BoundaryNorm(LEVELS + [1e5], CMAP_CLOUD.N),
                        shading="auto", zorder=2)
    # ---- 等值线: 1/10/100 ppm ----
    for thr, _, color, _ in ZONE_META:
        ls = "-" if thr >= 100 else "--"
        lw = 2.0 if thr >= 100 else 1.4
        cs = ax.contour(Xe.T, Ne.T, C.T, levels=[thr], colors=[color],
                        linewidths=[lw], linestyles=[ls], zorder=3)
        if thr == 100.0:
            try:
                ax.clabel(cs, fmt="%d ppm", fontsize=9, inline=True, colors=[color])
            except Exception:
                pass

    # ---- 3) 应急风险分区框 (旋转坐标外接矩形) ----
    rot = Affine2D().rotate_deg(beta_deg) + ax.transData
    for thr, expand, edge, ls, lw, label, fc in (
            (10.0, CORE_EXPAND, "#D50000", "-",  2.4, "核心应急隔离区", "#D50000"),
            (1.0,  WARN_EXPAND, "#B58900", (0, (6, 4)), 1.8, "警戒疏散区域", "#B58900")):
        cx, cy, w, h = bbox_rot(field, params, t, thr, expand)
        rect = Rectangle((cx - w / 2.0, cy - h / 2.0), w, h,
                         fill=False, edgecolor=edge, linestyle=ls, lw=lw,
                         zorder=4, transform=rot)
        ax.add_patch(rect)
        ax.text(cx, cy + h / 2.0 + 12, label, color=edge, fontsize=10.5,
                fontweight="bold", ha="center", rotation=beta_deg,
                transform=rot, zorder=5)

    # ---- 泄漏井口 ----
    ax.plot(0.0, 0.0, "*", ms=17, color="#D50000", mec="white", mew=1.2,
            zorder=6)
    ax.text(0.0, -18, "泄漏井口", ha="center", fontsize=9, color="#D50000",
            fontweight="bold", zorder=6)

    # ---- 1) 传感器点位 ----
    for r in sensor_readings(field, params, t):
        fc = "#D50000" if r["alarm"] else "#2E9E44"
        ax.plot(r["E"], r["N"], "D", ms=13, color=fc, mec="white", mew=1.4,
                zorder=7)
        txt = f"{r['name']} {r['desc']} {r['ppm']:.2f} ppm"
        ax.text(r["E"], r["N"] + 26, txt, ha="center", fontsize=8.6,
                color=fc, fontweight="bold", zorder=8,
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec=fc,
                          alpha=0.88, lw=0.8))

    # ---- 人员疏散路线 (危险区 -> 外围安全区) ----
    front = float(np.max(np.where(C >= 1.0, xd[:, None], 0.0)))
    half = float(np.max(np.abs(np.where(C >= 1.0, yc[None, :], 0.0))))
    ext = max(420.0, (front + half) * 0.75)
    for dE, dN in ((0.0, 1.0), (-0.72, 0.72), (-1.0, 0.0)):
        p0 = np.array([0.0, 0.0])
        p1 = np.array([dE, dN]) * ext
        arr = FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=26,
                              lw=2.6, color="#1F6FB2", alpha=0.95, zorder=5)
        ax.add_patch(arr)
    ax.text(8.0, ext * 0.94, "人员疏散路线", color="#1F6FB2", fontsize=10,
            fontweight="bold", zorder=6)

    # ---- 风向标注 (地理坐标) ----
    m = ext
    ax.annotate("", xy=(m * 0.35, 0.0), xytext=(0.0, 0.0),
                arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1.6))
    ax.text(m * 0.37, m * 0.02, f"下风向 {beta_deg:.0f}°", fontsize=9,
            color="#555555", va="center")

    # ---- 图幅自适应 ----
    lim = max(500.0, (front + half) * 0.707 * 1.15 + 260.0)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.grid(True, ls=":", color="#999999", alpha=0.6, zorder=1)
    ax.set_xlabel("东向距离 E (m)")
    ax.set_ylabel("北向距离 N (m)")
    ax.set_title(f"t = {t/60:.0f} 分钟   峰值 {np.max(C):.3g} ppm   1ppm前锋 {front:.0f} m",
                 fontsize=12, fontweight="bold")
    return pcm


def legend_handles():
    h = [
        Patch(facecolor="#FFF7B0", edgecolor="#999999", alpha=0.9,
              label="警戒区 1–10 ppm (浅黄)"),
        Patch(facecolor="#FF8C1A", edgecolor="#999999", alpha=0.9,
              label="危险区 10–100 ppm (橙)"),
        Patch(facecolor="#B30000", edgecolor="#999999", alpha=0.9,
              label="毒性区 >100 ppm (深红)"),
        Line2D([0], [0], color="#C00000", lw=2.0,
               label="100 ppm 毒性阈值边界"),
        Line2D([0], [0], color="#FF8C1A", lw=1.4, ls="--",
               label="10 ppm 危险区边界"),
        Line2D([0], [0], color="#E6B800", lw=1.4, ls="--",
               label="1 ppm 警戒区边界"),
        Line2D([0], [0], color="#D50000", lw=2.4,
               label="核心应急隔离区 (红实线)"),
        Line2D([0], [0], color="#B58900", lw=1.8, ls=(0, (6, 4)),
               label="警戒疏散区域 (黄虚线)"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor="#2E9E44",
               markeredgecolor="white", ms=9, label="传感器 正常 (绿)"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor="#D50000",
               markeredgecolor="white", ms=9, label="传感器 告警 (红, ≥10 ppm)"),
        Line2D([0], [0], color="#1F6FB2", lw=2.6,
               label="人员疏散路线 (危险区→安全区)"),
    ]
    return h


def write_sensor_csv(field, params, times, path):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["油井 H2S 泄漏 - 传感器实时读数(模型示意值, SLAB)"])
        w.writerow(["传感器距井口 60 m: S1正北 S2正西 S3正南 S4正东; "
                    "告警阈值 10 ppm"])
        w.writerow([])
        head = ["时刻(min)", "传感器", "方位", "E(m)", "N(m)", "浓度(ppm)", "状态"]
        w.writerow(head)
        for t in times:
            for r in sensor_readings(field, params, t):
                w.writerow([f"{t/60:g}", r["name"], r["desc"], r["E"], r["N"],
                            f"{r['ppm']:.3g}", "告警" if r["alarm"] else "正常"])


def main():
    ap = argparse.ArgumentParser(description="H2S 泄漏扩散核心可视化 (SLAB)")
    ap.add_argument("--input", default=None,
                    help="已有 SLAB 输出文件(如 results_time_slices_v1/predict_h2s.txt)")
    ap.add_argument("--outdir", default="results_risk_v1", help="输出目录")
    ap.add_argument("--times", default="60,180,300,900",
                    help="时刻(秒,逗号分隔)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    params = dict(base.CONFIG)
    times = [float(v) for v in args.times.split(",")]
    field = load_field(args.outdir, args.input)

    beta_deg = (params["wind_dir"] + 180.0) % 360.0

    # ---- 单时刻大图 (以 15 分钟为例) ----
    fig, ax = plt.subplots(figsize=(10.5, 9.5))
    draw_panel(ax, field, params, times[-1], beta_deg)
    fig.legend(handles=legend_handles(), loc="lower center", ncol=3,
               fontsize=9.5, frameon=True, columnspacing=1.4)
    fig.suptitle(f"油井 H2S 泄漏核心可视化 (t={times[-1]/60:.0f} 分钟, "
                 f"源强 {params['qs']:g} kg/s, 风速 {params['ua']:g} m/s, "
                 f"{base._wind_label(params)})",
                 fontsize=14, fontweight="bold")
    fig.subplots_adjust(left=0.07, right=0.98, top=0.90, bottom=0.14)
    p_hero = os.path.join(args.outdir, f"fig6_核心可视化_{times[-1]/60:.0f}min.png")
    fig.savefig(p_hero, dpi=150)
    plt.close(fig)

    # ---- 四时刻合并图 (2x2) ----
    fig, axes = plt.subplots(2, 2, figsize=(20.5, 18.5))
    pcm = None
    for k, (ax, t) in enumerate(zip(axes.ravel(), times)):
        pcm = draw_panel(ax, field, params, t, beta_deg)
        if k % 2 == 1:
            ax.set_yticklabels([])
        if k < 2:
            ax.set_xticklabels([])
    fig.legend(handles=legend_handles(), loc="lower center", ncol=3,
               fontsize=10.5, frameon=True, columnspacing=1.6)
    fig.suptitle(f"油井 H2S 泄漏核心可视化: 传感器 + 重气云团 + 应急分区 "
                 f"(源强 {params['qs']:g} kg/s, 风速 {params['ua']:g} m/s, "
                 f"{base._wind_label(params)}, 传感器距井口 60 m)",
                 fontsize=16, fontweight="bold")
    fig.subplots_adjust(left=0.06, right=0.99, top=0.93, bottom=0.10,
                        wspace=0.10, hspace=0.16)
    p_all = os.path.join(args.outdir, "fig6_核心可视化_四时刻.png")
    fig.savefig(p_all, dpi=150)
    plt.close(fig)

    write_sensor_csv(field, params, times,
                     os.path.join(args.outdir, "传感器读数汇总.csv"))

    print("输出目录:", args.outdir)
    print("  ", os.path.basename(p_hero))
    print("  ", os.path.basename(p_all))
    print("   传感器读数汇总.csv")


if __name__ == "__main__":
    main()