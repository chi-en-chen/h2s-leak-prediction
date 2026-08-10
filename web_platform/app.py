# -*- coding: utf-8 -*-
"""h2s_web 网页平台 - Flask 服务入口.

启动:  python web_platform/app.py
访问:  http://127.0.0.1:5000
"""
import os
import socket
import threading
import time
import uuid

from flask import Flask, abort, jsonify, render_template, request, send_from_directory

import runner

app = Flask(__name__)

JOBS = {}
LOCK = threading.Lock()


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/run", methods=["POST"])
def api_run():
    form = request.get_json(force=True, silent=True) or {}
    job_id = time.strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex[:6]
    try:
        params = runner.make_params(form)
    except Exception as e:
        return jsonify(ok=False, error="参数解析失败: %s" % e), 400
    job = runner.Job(job_id, params)
    with LOCK:
        JOBS[job_id] = job
    threading.Thread(target=runner.run_job, args=(job,), daemon=True).start()
    return jsonify(ok=True, job_id=job_id)


@app.route("/api/status/<job_id>")
def api_status(job_id):
    job = JOBS.get(job_id)
    if job is None:
        return jsonify(ok=False, error="任务不存在"), 404
    return jsonify(ok=True, job_id=job_id, status=job.status, progress=job.progress,
                   logs=job.logs[-100:], results=job.results, error=job.error)


@app.route("/api/files/<job_id>/<path:name>")
def api_files(job_id, name):
    job = JOBS.get(job_id)
    if job is None or not os.path.isdir(job.outdir):
        abort(404)
    return send_from_directory(job.outdir, name)


def _lan_ip():
    """获取本机局域网 IP (用于分享给其他电脑访问)."""
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("8.8.8.8", 80))
        ip = s.getsockname()[0]
        s.close()
        return ip
    except Exception:
        return "127.0.0.1"


def _write_access_tip(ip):
    try:
        lines = [
            "==================================================",
            "  油井硫化氢(H2S)泄漏扩散预测平台 - 访问地址",
            "--------------------------------------------------",
            "  本机访问:   http://127.0.0.1:5000",
            "  局域网访问: http://%s:5000" % ip,
            "",
            "  提示: 请确保本机防火墙已放行 5000 端口",
            "        (运行 开启局域网访问.bat 可自动放行)",
            "==================================================",
        ]
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "访问地址.txt")
        with open(path, "w", encoding="gbk") as f:
            f.write("\r\n".join(lines) + "\r\n")
    except Exception:
        pass


if __name__ == "__main__":
    ip = _lan_ip()
    _write_access_tip(ip)
    print("H2S 泄漏扩散预测平台已启动")
    print("本机访问:   http://127.0.0.1:5000")
    print("局域网访问: http://%s:5000" % ip)
    app.run(host="0.0.0.0", port=5000, debug=False, threaded=True)