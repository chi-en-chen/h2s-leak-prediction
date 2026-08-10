# -*- coding: utf-8 -*-
"""
油井 H2S 泄漏扩散 - 分时段(1/3/5/15 min)云团特征分析 (基于 SLAB 模型)
输出: 云团厚度 / 横向铺展范围 / 空间浓度分布 / 浓度热力云图 / 扩散动态动画
"""
import argparse, csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
from matplotlib import font_manager

from slab_h2s_lib import (H2S_PROPS, SlabField, default_h2s_params,
                          parse_predict, run_slab, write_slab_input)
import h2s_sensor_analysis as base

CONFIG = dict(base.CONFIG)

# ---- 计算网格 / 出图范围 (m) -----------------------------------
HORIZ_XMAX = 3000.0   # 水平足迹图下风向范围
VERT_XMAX  = 2500.0   # 垂直剖面图下风向范围
YMAX       = 400.0    # 横风向半范围
ZMAX       = 250.0    # 垂直方向最大高度
THR_MAIN   = 1.0      # 云团边界阈值 ppm (TLV 参考)
THRESHOLDS = (1.0, 10.0, 100.0)
CMAP       = "jet"
VMAX_PPM   = 1.0e5
VMIN_PPM   = 0.01


def build_input(params, path):
    """生成 SLAB 输入文件(复用 h2s_sensor_analysis 的源类型几何换算)."""
    return base.build_input(params, path)


def load_field(args, params):
    out = os.path.join(args.outdir, "predict_h2s.txt")
    if not args.no_run:
        build_input(params, os.path.join(args.outdir, "h2s_time_slices_input.txt"))
        out = run_slab(exe_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "Slab.exe"),
                       input_path=os.path.join(args.outdir, "h2s_time_slices_input.txt"),
                       output_path=out, workdir=args.outdir)
    else:
        out = os.path.abspath(args.input) if args.input else out
    return SlabField(parse_predict(out), name=os.path.basename(out))


# ------------------------------------------------------------
# 单个时刻的云团特征
# ------------------------------------------------------------
def slice_metrics(field, params, t, nx=380, ny=381, nz=120):
    zp = params["zp"]
    xs = np.linspace(0.0, HORIZ_XMAX, nx)
    ys = np.linspace(-YMAX, YMAX, ny)
    zs = np.linspace(0.0, ZMAX, nz)
    ccl = field.ppm(xs, 0.0, zp, t)                 # 中心线(呼吸带)浓度
    cg0 = field.ppm(xs, 0.0, 0.0, t)                # 中心线地面浓度
    Cxy = field.ppm(xs[:, None], ys[None, :], zp, t)  # 水平足迹 (nx, ny)
    Cxz = field.ppm(xs[:, None], 0.0, zs[None, :], t) # 中心线垂直剖面 (nx, nz)

    i_peak = int(np.argmax(ccl))
    x_peak = float(xs[i_peak])
    peak   = float(ccl[i_peak])

    res = dict(t=t, xs=xs, ys=ys, zs=zs, ccl=ccl, cg0=cg0,
               Cxy=Cxy, Cxz=Cxz, i_peak=i_peak, x_peak=x_peak, peak=peak)
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]

    for thr in THRESHOLDS:
        ab = ccl >= thr
        front = float(xs[ab][-1]) if ab.any() else 0.0

        # 云团厚度: 中心线各下风向位置, 浓度>=thr 的最高高度
        maskz = Cxz >= thr
        has = maskz.any(axis=1)
        thick_x = np.zeros(nx)
        if has.any():
            last = nz - 1 - np.argmax(maskz[has, ::-1], axis=1)
            thick_x[has] = zs[last]
        thick_at_peak = float(thick_x[i_peak])
        thick_max = float(thick_x.max())

        # 横向铺展: 各下风向位置, 浓度>=thr 的最大 |y|
        half = np.max(np.abs(np.where(Cxy >= thr, ys, 0.0)), axis=1)
        half_at_peak = float(half[i_peak])
        half_max = float(half.max())

        # 覆盖面积 (呼吸带水平面)
        area = float(np.count_nonzero(Cxy >= thr) * dx * dy)

        res[thr] = dict(front=front, thick_x=thick_x,
                        thick_at_peak=thick_at_peak, thick_max=thick_max,
                        half=half, half_at_peak=half_at_peak,
                        half_max=half_max, width_max=2.0 * half_max,
                        area=area)
    return res


def wind_vec(params):
    beta = np.radians((params["wind_dir"] + 180.0) % 360.0)
    return np.sin(beta), np.cos(beta)


# ------------------------------------------------------------
# 剖面 CSV (空间浓度分布)
# ------------------------------------------------------------
def write_profiles(s, params, path):
    thr = THR_MAIN
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow([f"t={s['t']/60:g} min 中心线浓度剖面 (y=0, z={params['zp']:g} m)"])
        w.writerow(["下风向x(m)", "呼吸带ppm", "地面ppm"])
        for i in range(len(s["xs"])):
            w.writerow([f"{s['xs'][i]:.1f}", f"{s['ccl'][i]:.4g}", f"{s['cg0'][i]:.4g}"])
        w.writerow([])
        w.writerow([f"峰值位置 x={s['x_peak']:.0f} m 处横向剖面 (z={params['zp']:g} m)"])
        w.writerow(["横风y(m)", "ppm"])
        row = s["Cxy"][s["i_peak"], :]
        for i in range(len(s["ys"])):
            w.writerow([f"{s['ys'][i]:.1f}", f"{row[i]:.4g}"])
        w.writerow([])
        w.writerow([f"峰值位置 x={s['x_peak']:.0f} m 处垂向剖面 (y=0)"])
        w.writerow(["高度z(m)", "ppm"])
        col = s["Cxz"][s["i_peak"], :]
        for i in range(len(s["zs"])):
            w.writerow([f"{s['zs'][i]:.1f}", f"{col[i]:.4g}"])


# ------------------------------------------------------------
# 浓度热力云图: 上排=水平足迹, 下排=中心线垂直剖面
# ------------------------------------------------------------
def fig_heatmaps(field, params, slices, path):
    times = params["times"]
    fig, axes = plt.subplots(2, len(times), figsize=(21, 9.5))
    pcm_h = pcm_v = None
    ex, ey, ez = np.array([0]), np.array([0]), np.array([0])
    ux, uy = wind_vec(params)
    for j, s in enumerate(slices):
        xs, ys, zs, Cxy, Cxz = s["xs"], s["ys"], s["zs"], s["Cxy"], s["Cxz"]
        Cm_h = np.ma.masked_less(Cxy.T, VMIN_PPM)
        Cm_v = np.ma.masked_less(Cxz.T, VMIN_PPM)

        ax = axes[0][j]
        pcm_h = ax.pcolormesh(xs, ys, Cm_h, norm=mcolors.LogNorm(VMIN_PPM, VMAX_PPM),
                              cmap=CMAP, shading="auto")
        ax.contour(xs, ys, Cxy.T, levels=[1.0, 10.0, 100.0, 1000.0],
                   colors="k", linewidths=0.5, alpha=0.6)
        ax.quiver(0.02 * HORIZ_XMAX, -0.86 * YMAX, ux, uy, scale=14,
                  color="white", width=0.02)
        ax.text(0.02 * HORIZ_XMAX + 40, -0.86 * YMAX + 40, "风",
                color="white", fontsize=10, fontweight="bold")
        # 传感器位置(下风向坐标)
        for sn in params["sensors"]:
            xd, yc = base.wind_to_downwind(params, sn["E"], sn["N"])
            ax.plot(xd, yc, "w^", ms=6)
            ax.text(xd + 12, yc + 12, sn["name"], color="white", fontsize=7,
                    fontweight="bold")
        m1 = s[1.0]
        m100 = s[100.0]
        txt = (f"峰值 {s['peak']:.3g} ppm\n"
               f"1ppm前锋 {m1['front']:.0f} m | 100ppm {m100['front']:.0f} m\n"
               f"最大宽度 {m1['width_max']:.0f} m | 面积 {m1['area']/1e4:.1f} 万m2")
        ax.text(0.02, 0.97, txt, transform=ax.transAxes, va="top", fontsize=8.5,
                color="white",
                bbox=dict(boxstyle="round,pad=0.3", fc="#222", ec="none", alpha=0.72))
        ax.set_title(f"t = {times[j]/60:g} 分钟", fontsize=13, fontweight="bold")
        ax.set_xlim(0, HORIZ_XMAX); ax.set_ylim(-YMAX, YMAX)
        ax.set_xlabel("下风向距离 (m)"); ax.set_ylabel("横风向 (m)")
        ax.set_aspect("auto")
        if j > 0:
            ax.set_yticklabels([])

        ax = axes[1][j]
        pcm_v = ax.pcolormesh(xs, zs, Cm_v, norm=mcolors.LogNorm(VMIN_PPM, VMAX_PPM),
                              cmap=CMAP, shading="auto")
        ax.contour(xs, zs, Cxz.T, levels=[1.0, 10.0, 100.0],
                   colors="k", linewidths=0.5, alpha=0.6)
        ax.axhline(params["zp"], color="w", lw=0.8, ls="--", alpha=0.8)
        m1 = s[1.0]; m10 = s[10.0]; m100 = s[100.0]
        txt = (f"1ppm厚度 {m1['thick_max']:.0f} m (峰值处 {m1['thick_at_peak']:.0f} m)\n"
               f"10ppm厚度 {m10['thick_max']:.0f} m | 100ppm厚度 {m100['thick_max']:.0f} m")
        ax.text(0.02, 0.97, txt, transform=ax.transAxes, va="top", fontsize=8.5,
                color="white",
                bbox=dict(boxstyle="round,pad=0.3", fc="#222", ec="none", alpha=0.72))
        ax.set_xlim(0, VERT_XMAX); ax.set_ylim(0, ZMAX)
        ax.set_xlabel("下风向距离 (m)"); ax.set_ylabel("高度 (m)")
        if j > 0:
            ax.set_yticklabels([])

    cb = fig.colorbar(pcm_h, ax=list(axes[0]), shrink=0.85, pad=0.015)
    cb.set_label("浓度 (ppm, 对数刻度)")
    cb2 = fig.colorbar(pcm_v, ax=list(axes[1]), shrink=0.85, pad=0.015)
    cb2.set_label("浓度 (ppm, 对数刻度)")
    fig.suptitle(f"H2S 云团浓度热力云图 (源强 {params['qs']:g} kg/s, 风速 "
                 f"{params['ua']:g} m/s, {base._wind_label(params)})",
                 fontsize=15, fontweight="bold")
    fig.subplots_adjust(left=0.05, right=0.90, top=0.90, bottom=0.08,
                        wspace=0.14, hspace=0.32)
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ------------------------------------------------------------
# 扩散动态动画 (水平足迹 + 垂直剖面)
# ------------------------------------------------------------
def make_animation(field, params, path_gif, path_mp4, t_end=900.0, step=5.0,
                   fps=10, nx=300, ny=181, nz=100):
    times = np.arange(0.0, t_end + 1e-6, step)
    zp = params["zp"]
    xs = np.linspace(0.0, HORIZ_XMAX, nx)
    ys = np.linspace(-YMAX, YMAX, ny)
    zs = np.linspace(0.0, ZMAX, nz)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.6, 5.8))
    ux, uy = wind_vec(params)

    C0h = np.ma.masked_less(np.zeros((ny, nx)) + VMIN_PPM, VMIN_PPM)
    C0v = np.ma.masked_less(np.zeros((nz, nx)) + VMIN_PPM, VMIN_PPM)
    pcm_h = ax1.pcolormesh(xs, ys, C0h, norm=mcolors.LogNorm(VMIN_PPM, VMAX_PPM),
                           cmap=CMAP, shading="auto")
    pcm_v = ax2.pcolormesh(xs, zs, C0v, norm=mcolors.LogNorm(VMIN_PPM, VMAX_PPM),
                           cmap=CMAP, shading="auto")
    ax1.quiver(0.02 * HORIZ_XMAX, -0.86 * YMAX, ux, uy, scale=14,
               color="white", width=0.02)
    ax1.text(0.02 * HORIZ_XMAX + 40, -0.86 * YMAX + 40, "风", color="white",
             fontsize=10, fontweight="bold")
    ax2.axhline(zp, color="w", lw=0.8, ls="--", alpha=0.8)
    ax2.text(0.99 * VERT_XMAX, zp + 8, f"呼吸带 {zp:g} m", color="white",
             fontsize=8, ha="right")
    for sn in params["sensors"]:
        xd, yc = base.wind_to_downwind(params, sn["E"], sn["N"])
        ax1.plot(xd, yc, "w^", ms=5)
    ax1.set_xlim(0, HORIZ_XMAX); ax1.set_ylim(-YMAX, YMAX)
    ax1.set_xlabel("下风向距离 (m)"); ax1.set_ylabel("横风向 (m)")
    ax1.set_title("水平足迹 (呼吸带 z=%.0f m)" % zp)
    ax2.set_xlim(0, VERT_XMAX); ax2.set_ylim(0, ZMAX)
    ax2.set_xlabel("下风向距离 (m)"); ax2.set_ylabel("高度 (m)")
    ax2.set_title("中心线垂直剖面")
    cb = fig.colorbar(pcm_h, ax=[ax1, ax2], shrink=0.9, pad=0.02)
    cb.set_label("浓度 (ppm, 对数刻度)")
    tinfo = ax1.text(0.02, 0.96, "", transform=ax1.transAxes, va="top",
                     fontsize=11, fontweight="bold", color="white",
                     bbox=dict(boxstyle="round,pad=0.3", fc="#222", ec="none",
                               alpha=0.75))
    stat = ax1.text(0.02, 0.74, "", transform=ax1.transAxes, va="top",
                    fontsize=8.5, color="white",
                    bbox=dict(boxstyle="round,pad=0.3", fc="#222", ec="none",
                              alpha=0.6))
    fig.subplots_adjust(left=0.06, right=0.90, top=0.94, bottom=0.10, wspace=0.18)

    def update(i):
        t = times[i]
        Cxy = field.ppm(xs[:, None], ys[None, :], zp, t)
        Cxz = field.ppm(xs[:, None], 0.0, zs[None, :], t)
        pcm_h.set_array(np.ma.masked_less(Cxy.T, VMIN_PPM).ravel())
        pcm_v.set_array(np.ma.masked_less(Cxz.T, VMIN_PPM).ravel())
        ccl = field.ppm(xs, 0.0, zp, t)
        ab = ccl >= THR_MAIN
        front = xs[ab][-1] if ab.any() else 0.0
        maskz = Cxz >= THR_MAIN
        thick = 0.0
        has = maskz.any(axis=1)
        if has.any():
            last = nz - 1 - np.argmax(maskz[has, ::-1], axis=1)
            thick = float(np.max(zs[last]))
        half = np.max(np.abs(np.where(Cxy >= THR_MAIN, ys, 0.0)), axis=1)
        width = 2.0 * float(half.max()) if half.size else 0.0
        tinfo.set_text(f"t = {t/60:5.1f} min ({t:.0f} s)")
        stat.set_text(f"1ppm 前锋 {front:.0f} m | 云团厚度 {thick:.0f} m\n"
                      f"最大横向宽度 {width:.0f} m")
        return pcm_h, pcm_v, tinfo, stat

    anim = FuncAnimation(fig, update, frames=len(times), blit=False)
    if path_gif:
        anim.save(path_gif, writer=PillowWriter(fps=fps))
    if path_mp4:
        try:
            anim.save(path_mp4, writer=FFMpegWriter(fps=fps, codec="libx264",
                                                    extra_args=["-pix_fmt", "yuv420p"]))
        except Exception as e:
            print("  MP4 保存失败(已跳过):", e)
    plt.close(fig)
    return len(times)


# ------------------------------------------------------------
# 汇总输出
# ------------------------------------------------------------
def write_summary(field, params, slices, args, path):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["油井 H2S 泄漏扩散 - 分时段云团特征 (SLAB 模型)"])
        w.writerow([])
        w.writerow(["一、工况参数"])
        w.writerow(["泄漏速率 qs (kg/s)", params["qs"]])
        w.writerow(["持续泄漏 tsd (s)", params["tsd"]])
        w.writerow(["风速 ua (m/s)", params["ua"]])
        w.writerow(["气温 ta (K)", params["ta"]])
        w.writerow(["气压 press (hPa)", params["press"]])
        w.writerow(["风向(来向) wind_dir (度)", params["wind_dir"]])
        w.writerow(["大气稳定度 stab", params["stab"]])
        w.writerow(["分析高度 zp (m)", params["zp"]])
        w.writerow(["说明: 四组工况为同一瞬态扩散过程在 1/3/5/15 分钟时刻的切片; "
                    "云团边界阈值取 1 ppm (TLV 参考)."])
        w.writerow([])
        w.writerow(["二、各时刻云团特征"])
        w.writerow(["时刻(min)", "中心线峰值(ppm)", "峰值位置(m)",
                    "1ppm前锋(m)", "10ppm前锋(m)", "100ppm前锋(m)",
                    "1ppm厚度(m)", "10ppm厚度(m)", "100ppm厚度(m)",
                    "1ppm最大半宽(m)", "1ppm最大宽度(m)", "100ppm最大宽度(m)",
                    "1ppm覆盖面积(m2)"])
        for s in slices:
            w.writerow([f"{s['t']/60:g}", f"{s['peak']:.3g}",
                        f"{s['x_peak']:.0f}",
                        f"{s[1.0]['front']:.0f}", f"{s[10.0]['front']:.0f}",
                        f"{s[100.0]['front']:.0f}",
                        f"{s[1.0]['thick_max']:.0f}", f"{s[10.0]['thick_max']:.0f}",
                        f"{s[100.0]['thick_max']:.0f}",
                        f"{s[1.0]['half_max']:.0f}", f"{s[1.0]['width_max']:.0f}",
                        f"{s[100.0]['width_max']:.0f}",
                        f"{s[1.0]['area']:.0f}"])
        w.writerow([])
        w.writerow(["三、空间浓度分布文件"])
        for t in params["times"]:
            w.writerow([f"t={t/60:g}min", f"剖面_{t/60:g}min.csv"])
        w.writerow([])
        ab, rel = field.validate()
        if ab is not None:
            w.writerow(["模型校验", f"与 SLAB z=0 表最大绝对误差 {ab:.2e}, "
                        f"最大相对误差 {rel*100:.1f}%"])


def print_summary(params, slices, args):
    print("=" * 70)
    print("  油井 H2S 泄漏扩散 - 分时段云团特征 (SLAB)")
    print("=" * 70)
    wd = (params["wind_dir"] + 180.0) % 360.0
    print(f"  源强 {params['qs']:g} kg/s, 持续 {params['tsd']/60:g} min, "
          f"风速 {params['ua']:g} m/s, 风向 {params['wind_dir']:.0f} 度(来向) -> 下风向 {wd:.0f} 度")
    print("-" * 70)
    print("  时刻   峰值ppm  峰值x   1ppm前锋  1ppm厚度  最大宽度  1ppm面积")
    for s in slices:
        m1 = s[1.0]
        print(f"  {s['t']/60:4.0f}min {s['peak']:9.3g} {s['x_peak']:6.0f}m "
              f"{m1['front']:8.0f}m {m1['thick_max']:7.0f}m {m1['width_max']:8.0f}m "
              f"{m1['area']/1e4:8.1f}万m2")
    print("-" * 70)
    print(f"  输出目录: {args.outdir}")
    print("=" * 70)


def main():
    ap = argparse.ArgumentParser(description="油井 H2S 分时段扩散分析 (SLAB)")
    ap.add_argument("--no-run", action="store_true", help="不重跑 Slab.exe, 复用已有输出")
    ap.add_argument("--input", default=None, help="已有 SLAB 输出文件(配合 --no-run)")
    ap.add_argument("--outdir", default="results_time_slices_v1", help="输出目录")
    ap.add_argument("--times", default="60,180,300,900",
                    help="分析时刻(秒, 逗号分隔), 默认 1/3/5/15 分钟")
    ap.add_argument("--t-end", type=float, default=900.0, help="动画结束时刻(秒)")
    ap.add_argument("--step", type=float, default=5.0, help="动画时间步长(秒)")
    ap.add_argument("--no-gif", action="store_true", help="不生成 GIF 动画")
    ap.add_argument("--no-mp4", action="store_true", help="不生成 MP4 动画")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    params = dict(CONFIG)
    params["times"] = [float(v) for v in args.times.split(",")]

    print(">> 运行 SLAB 计算扩散场 ...")
    field = load_field(args, params)

    print(">> 计算 1/3/5/15 分钟云团特征 ...")
    slices = [slice_metrics(field, params, t) for t in params["times"]]

    print(">> 输出空间浓度分布(剖面 CSV) ...")
    for t, s in zip(params["times"], slices):
        write_profiles(s, params, os.path.join(args.outdir, f"剖面_{t/60:g}min.csv"))

    print(">> 绘制浓度热力云图 ...")
    fig_heatmaps(field, params, slices,
                 os.path.join(args.outdir, "fig5_浓度热力云图_四时刻.png"))

    print(">> 制作扩散动态动画 ...")
    n_frames = make_animation(
        field, params,
        None if args.no_gif else os.path.join(args.outdir, "h2s_扩散动画.gif"),
        None if args.no_mp4 else os.path.join(args.outdir, "h2s_扩散动画.mp4"),
        t_end=args.t_end, step=args.step)
    print(f"    动画帧数: {n_frames}")

    write_summary(field, params, slices, args,
                  os.path.join(args.outdir, "h2s_time_slices_summary.csv"))
    print_summary(params, slices, args)


if __name__ == "__main__":
    main()