# -*- coding: utf-8 -*-
"""h2s_web runner - 网页平台计算任务调度(封装已有 SLAB 分析脚本)."""

import os

import sys

import time

import traceback

from types import SimpleNamespace

import numpy as np

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

if ROOT not in sys.path:

    sys.path.insert(0, ROOT)

import h2s_sensor_analysis as base

import h2s_time_slices as tsl

import h2s_risk_visual as rv

JOBS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "jobs")

def _make_sensor_layout(dist):
    """根据距井口距离生成四传感器布局 (正北/正西/正南/正东)."""

    return [

        dict(name="S1", E=0.0,    N=dist,   z=1.5, desc="正北"),
        dict(name="S2", E=-dist,  N=0.0,    z=1.5, desc="正西"),
        dict(name="S3", E=0.0,    N=-dist,  z=1.5, desc="正南"),
        dict(name="S4", E=dist,   N=0.0,    z=1.5, desc="正东"),
    ]

# 传感器布局: 正北/正西/正南/正东, 距井口可调(默认 60 m)

SENSOR_LAYOUT = _make_sensor_layout(60.0)

STAB_TO_CODE = {"A": 1, "B": 2, "C": 3, "D": 4, "E": 5, "F": 6}

CODE_TO_STAB = {v: k for k, v in STAB_TO_CODE.items()}

def _to_num(form, key, default):
    """数值参数: 空白回退默认值, 非法值给出清晰提示."""

    raw = form.get(key)

    if raw is None or str(raw).strip() == "":

        return default

    try:

        return float(str(raw).strip())

    except Exception:

        raise ValueError("字段 %s 不是有效数字: %r，请检查后重试" % (key, raw))

def _parse_times(raw):
    """解析分析时刻: 兼容中英文逗号、顿号、分号与空格; 默认 1,3,5,15 分钟."""

    s = str(raw).replace("，", ",").replace("、", ",").replace("；", ";").replace(";", ",")

    parts = [x.strip() for x in s.split(",") if x.strip()]

    times_min = []

    for x in parts:

        try:

            v = float(x)

        except Exception:

            raise ValueError("分析时刻包含无法识别的时间: %r，请用英文逗号分隔" % x)

        if v <= 0:

            raise ValueError("分析时刻必须为正数: %r" % x)

        if v not in times_min:

            times_min.append(v)

    if not times_min:

        times_min = [1.0, 3.0, 5.0, 15.0]

    times_min.sort()

    if times_min[0] < 0.1 or times_min[-1] > 180:

        raise ValueError("分析时刻请在 0.1~180 分钟之间选择")

    return [t * 60.0 for t in times_min]

def make_params(form):
    """网页表单 -> SLAB 工况参数 dict (时刻按分钟输入, 内部换算为秒)."""

    times = _parse_times(form.get("times", "1,3,5,15"))

    qs = _to_num(form, "qs", 1.0)

    if qs <= 0:

        raise ValueError("泄漏速率 qs 必须大于 0")

    ua = _to_num(form, "ua", 3.0)

    if ua <= 0:

        raise ValueError("风速 ua 必须大于 0")

    idspl = int(float(_to_num(form, "idspl", 1)))

    jet_d = _to_num(form, "jet_d", 0.25)

    jet_h = _to_num(form, "jet_h", 2.0)

    if idspl == 3 and jet_d <= 0:

        raise ValueError("喷口直径必须大于 0")

    if idspl == 3 and jet_h < 0:

        raise ValueError("喷口高度不能为负")

    p = dict(base.CONFIG)

    p.update(
        idspl=idspl,
        jet_d=jet_d,
        jet_h=jet_h,
        qs=qs,
        tsd=_to_num(form, "tsd_min", 15.0) * 60.0,
        as_=_to_num(form, "as_", 100.0),
        ua=ua,
        ta=_to_num(form, "ta_c", 27.0) + 273.15,
        press=_to_num(form, "press", 1013.25),
        rh=_to_num(form, "rh", 60.0),
        stab=STAB_TO_CODE.get(str(form.get("stab", "D")).upper(), 4),
        wind_dir=_to_num(form, "wind_dir", 225.0) % 360.0,
        zp=_to_num(form, "zp", 1.5),
        times=times,
        map_times=[t for t in times],
        sensor_dist=_to_num(form, "sensor_dist", 60.0),
        sensors=[dict(s) for s in _make_sensor_layout(_to_num(form, "sensor_dist", 60.0))],
        make_anim=bool(form.get("make_anim", True)),
        anim_min=_to_num(form, "anim_min", 15.0),
        wms=_to_num(form, "wms", 0.034081),
        tbp=_to_num(form, "tbp", 213.0),
        rhosl=_to_num(form, "rhosl", 949.0),
        dhe=_to_num(form, "dhe_kj", 548.0) * 1000.0,
        cps=_to_num(form, "cps", 1010.0),
        cpsl=_to_num(form, "cpsl", 1800.0),
    )

    return p

class Job:
    """后台计算任务."""

    def __init__(self, job_id, params):

        self.id = job_id

        self.params = params

        self.outdir = os.path.join(JOBS_DIR, job_id)

        self.status = "queued"      # queued / running / done / failed

        self.logs = []

        self.progress = 0

        self.error = None

        self.results = dict(files=[], metrics=[], sensors=[], summary={})

    def log(self, msg, progress=None):

        self.logs.append("[%s] %s" % (time.strftime("%H:%M:%S"), msg))

        if progress is not None:

            self.progress = max(self.progress, int(progress))

    def result_file(self, name, category, typ="image"):

        url = "/api/files/%s/%s" % (self.id, name.replace(" ", "%20"))

        return dict(name=name, url=url, category=category, type=typ)

def run_job(job):
    """执行完整计算流程: SLAB -> 分时段特征 -> 热力图 -> 动画 -> 核心可视化 -> 传感器分析."""

    job.status = "running"

    try:

        p = job.params

        os.makedirs(job.outdir, exist_ok=True)

        args = SimpleNamespace(outdir=job.outdir, no_run=False, input=None)

        job.log("正在生成 SLAB 输入文件并启动计算 ...", 5)

        field = tsl.load_field(args, p)

        ab, rel = field.validate()

        if ab is not None:

            job.log("SLAB 计算完成; 浓度重构校验 绝对误差 %.2e, 相对误差 %.1f%%" % (ab, rel * 100.0), 18)

        else:

            job.log("SLAB 计算完成", 18)

        job.log("正在计算 1/3/5/15 分钟云团特征(厚度/铺展/浓度分布) ...", 28)

        slices = [tsl.slice_metrics(field, p, t) for t in p["times"]]

        for t, s in zip(p["times"], slices):

            tsl.write_profiles(s, p, os.path.join(job.outdir, "剖面_%gmin.csv" % (t / 60.0)))

        job.log("分时段云团特征计算完成", 42)

        job.log("正在绘制浓度热力云图 ...", 46)

        tsl.fig_heatmaps(field, p, slices,
                         os.path.join(job.outdir, "fig5_浓度热力云图_四时刻.png"))

        job.log("浓度热力云图完成", 58)

        if p.get("make_anim", True):

            job.log("正在制作扩散动态动画(GIF, 约需 1-2 分钟) ...", 62)

            n = tsl.make_animation(field, p,
                                   os.path.join(job.outdir, "h2s_扩散动画.gif"),
                                   None,
                                   t_end=p.get("anim_min", 15.0) * 60.0, step=5.0)

            job.log("扩散动态动画完成(%d 帧)" % n, 74)

        else:

            job.log("已跳过动画生成(设置中关闭)", 74)

        job.log("正在绘制核心可视化图(传感器+云团+应急分区) ...", 78)

        rv.SENSOR_LAYOUT[:] = [dict(s) for s in p["sensors"]]

        beta_deg = (p["wind_dir"] + 180.0) % 360.0

        hero_path = os.path.join(job.outdir, "fig6_核心可视化_%gmin.png" % (p["times"][-1] / 60.0))

        fig, ax = plt.subplots(figsize=(10.5, 9.5))

        rv.draw_panel(ax, field, p, p["times"][-1], beta_deg)

        fig.legend(handles=rv.legend_handles(), loc="lower center", ncol=3,
                   fontsize=9.5, frameon=True, columnspacing=1.4)

        fig.suptitle("油井 H2S 泄漏核心可视化 (t=%g 分钟, 源强 %g kg/s, 风速 %g m/s, %s)"

                     % (p["times"][-1] / 60.0, p["qs"], p["ua"], base._wind_label(p)),
                     fontsize=14, fontweight="bold")

        fig.subplots_adjust(left=0.07, right=0.98, top=0.90, bottom=0.14)

        fig.savefig(hero_path, dpi=150)

        plt.close(fig)

        nt = len(p["times"])

        if nt >= 4:

            fig, axes = plt.subplots(2, 2, figsize=(20.5, 18.5))

            for k, (ax, t) in enumerate(zip(axes.ravel(), p["times"])):

                rv.draw_panel(ax, field, p, t, beta_deg)

                if k % 2 == 1:

                    ax.set_yticklabels([])

                if k < 2:

                    ax.set_xticklabels([])

            fig.legend(handles=rv.legend_handles(), loc="lower center", ncol=3,
                       fontsize=10.5, frameon=True, columnspacing=1.6)

            fig.suptitle("油井 H2S 泄漏核心可视化: 传感器 + 重气云团 + 应急分区 "
                         "(源强 %g kg/s, 风速 %g m/s, %s, 传感器距井口 %g m)"

                         % (p["qs"], p["ua"], base._wind_label(p),
                            p.get("sensor_dist", 60.0)),
                         fontsize=16, fontweight="bold")

            fig.subplots_adjust(left=0.06, right=0.99, top=0.93, bottom=0.10,
                                wspace=0.10, hspace=0.16)

        else:

            fig, ax = plt.subplots(figsize=(10.5, 9.5))

            rv.draw_panel(ax, field, p, p["times"][0], beta_deg)

        all_path = os.path.join(job.outdir, "fig6_核心可视化_四时刻.png")

        fig.savefig(all_path, dpi=150)

        plt.close(fig)

        job.log("核心可视化图完成", 88)

        job.log("正在绘制传感器分析图(布局/足迹/时间序列) ...", 90)

        base.fig_sensor_map(p, os.path.join(job.outdir, "fig1_传感器布局.png"))

        base.fig_footprints(field, p, os.path.join(job.outdir, "fig2_浓度足迹_多时刻.png"))

        base.fig_time_series(field, p, None, None,
                             os.path.join(job.outdir, "fig3_四传感器时间序列.png"))

        job.log("传感器分析图完成", 95)

        job.log("正在生成汇总数据 ...", 96)

        tsl.write_summary(field, p, slices, args,
                          os.path.join(job.outdir, "h2s_time_slices_summary.csv"))

        stats = base.sensor_stats(field, p)

        mat = base.matrix_at_times(field, p)

        base.write_summary(field, p, stats, mat, dict(n=0, peaks=[]),
                           np.empty((0, 4)), "网页平台预测(无实测)",
                           args, os.path.join(job.outdir, "h2s_sensor_summary.csv"))

        rv.write_sensor_csv(field, p, p["times"],
                            os.path.join(job.outdir, "传感器读数汇总.csv"))

        job.results = build_results(job, field, p, slices, stats)

        job.status = "done"

        job.log("全部计算完成", 100)

    except Exception:

        job.error = traceback.format_exc()

        job.status = "failed"

        job.log("计算失败: %s" % (job.error.splitlines()[-1] if job.error else ""))

def build_results(job, field, p, slices, stats):

    files = []

    files.append(job.result_file("fig6_核心可视化_四时刻.png", "核心可视化"))

    files.append(job.result_file("fig6_核心可视化_%gmin.png" % (p["times"][-1] / 60.0),
                                 "核心可视化"))

    files.append(job.result_file("fig5_浓度热力云图_四时刻.png", "浓度热力云图"))

    if p.get("make_anim", True):

        files.append(job.result_file("h2s_扩散动画.gif", "扩散动画"))

    files.append(job.result_file("fig1_传感器布局.png", "传感器分析"))

    files.append(job.result_file("fig2_浓度足迹_多时刻.png", "传感器分析"))

    files.append(job.result_file("fig3_四传感器时间序列.png", "传感器分析"))

    files.append(job.result_file("传感器读数汇总.csv", "数据文件", "csv"))

    files.append(job.result_file("h2s_sensor_summary.csv", "数据文件", "csv"))

    files.append(job.result_file("h2s_time_slices_summary.csv", "数据文件", "csv"))

    for t in p["times"]:

        files.append(job.result_file("剖面_%gmin.csv" % (t / 60.0), "数据文件", "csv"))

    metrics = []

    for s in slices:

        metrics.append(dict(
            t_min="%g" % (s["t"] / 60.0),
            peak_ppm="%.3g" % s["peak"],
            x_peak="%.0f" % s["x_peak"],
            front_1="%.0f" % s[1.0]["front"],
            front_10="%.0f" % s[10.0]["front"],
            front_100="%.0f" % s[100.0]["front"],
            thick_1="%.0f" % s[1.0]["thick_max"],
            thick_10="%.0f" % s[10.0]["thick_max"],
            width_1="%.0f" % s[1.0]["width_max"],
            area_1="%.1f" % (s[1.0]["area"] / 1e4),
        ))

    sensors = []

    for t in p["times"]:

        for r in rv.sensor_readings(field, p, t):

            sensors.append(dict(t_min="%g" % (t / 60.0), name=r["name"], desc=r["desc"],
                                ppm="%.4g" % r["ppm"],
                                mgm3="%.4g" % base.ppm_to_mgm3(r["ppm"], p),
                                status="告警" if r["alarm"] else "正常"))

    src_name = "地面蒸发池" if p["idspl"] == 1 else "垂直喷口"

    geom = "液池面积 %.3g m²" % p["as_"] if p["idspl"] == 1 else \
        "喷口直径 %.2g m, 离地高度 %.2g m" % (p["jet_d"], p["jet_h"])

    summary = dict(
        src_type="%s（%s）" % (src_name, geom),
        qs="%.3g kg/s" % p["qs"],
        tsd_min="%g min" % (p["tsd"] / 60.0),
        ua="%.2g m/s" % p["ua"],
        ta="%.1f C" % (p["ta"] - 273.15),
        press="%.1f hPa" % p["press"],
        rh="%.0f %%" % p["rh"],
        wind_dir="%.0f 度(来向) -> 下风向 %.0f 度" % (p["wind_dir"], (p["wind_dir"] + 180.0) % 360.0),
        stab=CODE_TO_STAB.get(p["stab"], str(p["stab"])),
        zp="%.1f m" % p["zp"],
        times=" / ".join("%g" % (t / 60.0) for t in p["times"]) + " 分钟",
        sensor_layout="S1正北 S2正西 S3正南 S4正东, 距井口 %g m" % p.get("sensor_dist", 60.0),
    )

    return dict(files=files, metrics=metrics, sensors=sensors, summary=summary)

if __name__ == "__main__":

    # 命令行自测: python runner.py --qs 1 --ua 3 ...

    import argparse

    ap = argparse.ArgumentParser()

    ap.add_argument("--qs", default="1.0"); ap.add_argument("--tsd_min", default="15")

    ap.add_argument("--ua", default="3.0"); ap.add_argument("--ta_c", default="27")

    ap.add_argument("--press", default="1013.25"); ap.add_argument("--rh", default="60")

    ap.add_argument("--wind_dir", default="225"); ap.add_argument("--stab", default="D")

    ap.add_argument("--zp", default="1.5"); ap.add_argument("--times", default="1,3,5,15")

    ap.add_argument("--anim", type=int, default=1)

    ap.add_argument("--anim_min", default="15")

    args = ap.parse_args().__dict__

    form = {k: str(v) for k, v in args.items() if k in ("qs", "tsd_min", "ua", "ta_c",
            "press", "rh", "wind_dir", "stab", "zp", "times", "anim_min")}

    form["make_anim"] = bool(args["anim"])

    params = make_params(form)

    job_id = time.strftime("selftest_%Y%m%d_%H%M%S")

    job = Job(job_id, params)

    run_job(job)

    print("STATUS:", job.status)

    for ln in job.logs:

        print(ln)

    if job.error:

        print(job.error)

    print("FILES:", [f["name"] for f in job.results["files"]])