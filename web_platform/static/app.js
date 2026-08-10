/* h2s_web 前端逻辑 */
"use strict";

const $ = (id) => document.getElementById(id);
const defaults = { jet_d: "0.25", jet_h: "2.0",
  qs: "1.0", tsd_min: "15", as_: "100", idspl: "1",
  ua: "3.0", wind_dir_sel: "225", wind_dir: "225",
  ta_c: "27", press: "1013.25", rh: "60", stab: "D",
  zp: "1.5", times: "1,3,5,15", make_anim: true, anim_min: "15"
};

let currentJob = null;
let pollTimer = null;

function resetForm() {
  for (const [k, v] of Object.entries(defaults)) {
    const el = $(k);
    if (!el) continue;
    if (el.type === "checkbox") el.checked = v;
    else el.value = v;
  }
  syncWindDir();
  syncSourceType();
}

function syncWindDir() {
  const sel = $("wind_dir_sel").value;
  const custom = sel === "custom";
  $("wind_dir_custom").classList.toggle("hidden", !custom);
  if (!custom) $("wind_dir").value = sel;
}

function syncSourceType() {
  const jet = $("idspl").value === "3";
  $("as_group").classList.toggle("hidden", jet);
  $("jet_group").classList.toggle("hidden", !jet);
}

function collectParams() {
  const windSel = $("wind_dir_sel").value;
  return {
    qs: $("qs").value, tsd_min: $("tsd_min").value, as_: $("as_").value,
    idspl: $("idspl").value, jet_d: $("jet_d").value, jet_h: $("jet_h").value, ua: $("ua").value,
    wind_dir: windSel === "custom" ? $("wind_dir").value : windSel,
    ta_c: $("ta_c").value, press: $("press").value, rh: $("rh").value,
    stab: $("stab").value, zp: $("zp").value, times: $("times").value.replace(/[，、；;]/g, ",").replace(/\s+/g, ""),
    make_anim: $("make_anim").checked, anim_min: $("anim_min").value
  };
}

function switchTab(name) {
  document.querySelectorAll(".tab").forEach((t) => {
    t.classList.toggle("active", t.dataset.tab === name);
  });
  document.querySelectorAll(".tab-pane").forEach((p) => {
    p.classList.toggle("hidden", p.id !== "tab-" + name);
  });
}

function setStatus(visible) {
  $("status-box").classList.toggle("hidden", !visible);
  $("result-box").classList.toggle("hidden", visible);
  $("empty-box").classList.toggle("hidden", visible);
}

async function runJob() {
  const btn = $("btn-run");
  btn.disabled = true;
  btn.textContent = "计算中…";
  $("log-box").textContent = "";
  $("progress-bar").style.width = "0%";
  $("progress-text").textContent = "正在提交任务 …";
  setStatus(true);
  try {
    const resp = await fetch("/api/run", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(collectParams())
    });
    let data;
    try { data = await resp.json(); } catch (_) { throw new Error("服务响应异常，请确认服务窗口仍在运行"); }
    if (!resp.ok || !data.ok) throw new Error(data.error || "提交失败（HTTP " + resp.status + "）");
    currentJob = data.job_id;
    pollTimer = setInterval(pollStatus, 1500);
  } catch (e) {
    $("progress-text").textContent = "提交失败：" + e.message + "。请检查输入参数后重试";
    btn.disabled = false;
    btn.textContent = "开始计算";
  }
}

async function pollStatus() {
  if (!currentJob) return;
  try {
    const resp = await fetch("/api/status/" + currentJob);
    const data = await resp.json();
    if (!data.ok) throw new Error(data.error || "状态获取失败");
    $("progress-bar").style.width = (data.progress || 0) + "%";
    $("progress-text").textContent = statusText(data);
    if (data.logs && data.logs.length) {
      $("log-box").textContent = data.logs.join("\n");
      $("log-box").scrollTop = $("log-box").scrollHeight;
    }
    if (data.status === "done") {
      clearInterval(pollTimer);
      renderResults(data.results);
      setStatus(false);
      $("btn-run").disabled = false;
      $("btn-run").textContent = "重新计算";
      switchTab("core");
    } else if (data.status === "failed") {
      clearInterval(pollTimer);
      $("progress-text").textContent = "计算失败，请检查参数后重试";
      $("log-box").textContent += "\n" + (data.error || "");
      $("btn-run").disabled = false;
      $("btn-run").textContent = "开始计算";
    }
  } catch (e) {
    clearInterval(pollTimer);
    $("progress-text").textContent = "状态获取失败：" + e.message;
    $("btn-run").disabled = false;
    $("btn-run").textContent = "开始计算";
  }
}

function statusText(data) {
  const map = {
    queued: "任务排队中 …", running: "正在计算 …", done: "计算完成 ✓", failed: "计算失败"
  };
  return (map[data.status] || data.status) + "（" + (data.progress || 0) + "%）";
}

function renderResults(res) {
  if (!res) return;
  const s = res.summary || {};
  const lines = [
    "工况：" + (s.src_type || "") + "；泄漏率 " + (s.qs || "-") + "，持续 " + (s.tsd_min || "-") + "，" +
    "风速 " + (s.ua || "-") + "，气温 " + (s.ta || "-") + "，气压 " + (s.press || "-") + "，" +
    "风向 " + (s.wind_dir || "-") + "，稳定度 " + (s.stab || "-") + "，分析高度 " + (s.zp || "-")
  ];
  $("summary-line").textContent = lines[0];

  // 图片
  const byName = {};
  (res.files || []).forEach((f) => { byName[f.name] = f.url; });
  const setImg = (id, key) => {
    if (byName[key]) $("img-" + id).src = byName[key];
  };
  const heroKey = Object.keys(byName).find((n) => /核心可视化_\d/.test(n));
  setImg("core-all", "fig6_核心可视化_四时刻.png");
  if (heroKey) setImg("core-hero", heroKey);
  setImg("heatmap", "fig5_浓度热力云图_四时刻.png");
  setImg("anim", "h2s_扩散动画.gif");

  // 传感器表格
  renderSensorTable(res.sensors || []);
  // 云团特征表格
  renderMetricsTable(res.metrics || []);
  // 文件下载
  renderFiles(res.files || []);
}

function renderSensorTable(rows) {
  const tbl = $("tbl-sensors");
  if (!rows.length) { tbl.innerHTML = ""; return; }
  const heads = ["时刻(min)", "传感器", "方位", "浓度(ppm)", "浓度(mg/m³)", "状态"];
  let h = "<tr>" + heads.map((x) => "<th>" + x + "</th>").join("") + "</tr>";
  const body = rows.map((r) => {
    const cls = r.status === "告警" ? "alarm" : "ok";
    return "<tr><td>" + r.t_min + "</td><td><b>" + r.name + "</b></td><td>" + r.desc +
      "</td><td>" + r.ppm + "</td><td>" + r.mgm3 + "</td><td class=\"" + cls + "\">" +
      r.status + "</td></tr>";
  }).join("");
  tbl.innerHTML = h + body;
}

function renderMetricsTable(rows) {
  const tbl = $("tbl-metrics");
  if (!rows.length) { tbl.innerHTML = ""; return; }
  const heads = [
    ["t_min", "时刻(min)"], ["peak_ppm", "峰值(ppm)"], ["x_peak", "峰值位置(m)"],
    ["front_1", "1ppm前锋(m)"], ["front_10", "10ppm前锋(m)"], ["front_100", "100ppm前锋(m)"],
    ["thick_1", "1ppm厚度(m)"], ["width_1", "1ppm最大宽度(m)"], ["area_1", "1ppm面积(万m²)"]
  ];
  let h = "<tr>" + heads.map(([, t]) => "<th>" + t + "</th>").join("") + "</tr>";
  const body = rows.map((r) =>
    "<tr>" + heads.map(([k]) => "<td>" + (r[k] ?? "-") + "</td>").join("") + "</tr>"
  ).join("");
  tbl.innerHTML = h + body;
}

function renderFiles(files) {
  const ul = $("file-list");
  const groups = {
    "核心可视化": "可视化图片", "浓度热力云图": "可视化图片", "扩散动画": "动态动画",
    "传感器分析": "可视化图片", "数据文件": "数据文件"
  };
  const order = ["核心可视化", "浓度热力云图", "扩散动画", "传感器分析", "数据文件"];
  ul.innerHTML = "";
  for (const cat of order) {
    const items = files.filter((f) => f.category === cat);
    if (!items.length) continue;
    const h = document.createElement("li");
    h.innerHTML = "<b>" + groups[cat] + "</b>";
    ul.appendChild(h);
    const sub = document.createElement("ul");
    items.forEach((f) => {
      const li = document.createElement("li");
      const a = document.createElement("a");
      a.href = f.url; a.textContent = f.name + (f.type === "csv" ? "（CSV）" : "");
      a.download = f.name;
      li.appendChild(a);
      sub.appendChild(li);
    });
    ul.appendChild(sub);
  }
}

document.addEventListener("DOMContentLoaded", () => {
  resetForm();
  syncSourceType();
  $("idspl").addEventListener("change", syncSourceType);
  $("wind_dir_sel").addEventListener("change", syncWindDir);
  $("btn-run").addEventListener("click", runJob);
  $("btn-reset").addEventListener("click", () => { resetForm(); setStatus(false); });
  document.querySelectorAll(".tab").forEach((t) => {
    t.addEventListener("click", () => switchTab(t.dataset.tab));
  });
});