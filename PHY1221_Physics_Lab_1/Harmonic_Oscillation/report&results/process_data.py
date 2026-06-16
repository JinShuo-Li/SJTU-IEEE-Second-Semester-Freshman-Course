#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data processing for the experiment report:

Run from the directory containing exp1.csv ... exp5.csv. If those exact names
are absent, the script also accepts files whose names start with exp1 ... exp5.
It generates all required fit figures plus summary.md and results.json.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# ----------------------------- constants ---------------------------------
G = 9.794  # m/s^2
M_T_G = 120.57
M_U_G = 129.53
M0_G = 21.74
X0_EXP3_CM = 62.55
X0_EXP4_CM = 62.55
D_M = 0.01010
X_SUPPORT1_CM = 31.80
X_SUPPORT2_CM = 117.98
H_CM = 2.102

M_T_KG = M_T_G / 1000.0
M_U_KG = M_U_G / 1000.0
M0_KG = M0_G / 1000.0

# --------------------------- plotting config ------------------------------
def setup_matplotlib() -> None:
    candidates = [
        "Noto Sans CJK SC", "Noto Serif CJK SC", "AR PL UMing CN",
        "AR PL SungtiL GB", "WenQuanYi Zen Hei", "SimHei",
    ]
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in candidates:
        if name in available:
            plt.rcParams["font.sans-serif"] = [name]
            plt.rcParams["font.family"] = "sans-serif"
            break
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["figure.dpi"] = 150


def find_csv(stem: str, directory: Path) -> Path:
    exact = directory / f"{stem}.csv"
    if exact.exists():
        return exact
    matches = sorted(directory.glob(f"{stem}*.csv"))
    if not matches:
        raise FileNotFoundError(f"Cannot find {stem}.csv or {stem}*.csv in {directory}")
    return matches[0]


def read_csv(stem: str, directory: Path) -> pd.DataFrame:
    path = find_csv(stem, directory)
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    return df


def linfit(x, y) -> Tuple[float, float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    yhat = slope * x + intercept
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot != 0 else float("nan")
    return float(slope), float(intercept), float(r2)


def make_fit_plot(x, y, slope, intercept, r2, xlabel, ylabel, title, out_path, eq_prefix="y") -> None:
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.scatter(x, y, marker="o", s=34, color="black", label="实验数据")
    xr = np.linspace(float(np.min(x)), float(np.max(x)), 200)
    ax.plot(xr, slope * xr + intercept, linestyle="--", linewidth=1.5, color="black", label="线性拟合")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, linestyle=":", linewidth=0.7, alpha=0.85)
    eq = f"{eq_prefix} = {slope:.6g} x {intercept:+.6g}\n$R^2$ = {r2:.6f}"
    ax.text(0.05, 0.95, eq, transform=ax.transAxes, va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="0.5", alpha=0.9))
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def make_period_amplitude_plot(x, y, slope, intercept, r2, mean_T, out_path) -> None:
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.scatter(x, y, marker="o", s=34, color="black", label="实验数据")
    xr = np.linspace(float(np.min(x)), float(np.max(x)), 200)
    ax.plot(xr, slope * xr + intercept, linestyle="--", linewidth=1.5, color="black", label="线性拟合")
    ax.axhline(mean_T, linestyle=":", linewidth=1.4, color="black", label="平均周期")
    ax.set_xlabel("振幅 A / m")
    ax.set_ylabel("周期 T / s")
    ax.set_title("实验4：周期 T 与振幅 A 的关系")
    ax.grid(True, linestyle=":", linewidth=0.7, alpha=0.85)
    eq = f"T = {slope:.6g} A {intercept:+.6g}\n$R^2$ = {r2:.6f}\n$\\bar{{T}}$ = {mean_T:.6f} s"
    ax.text(0.05, 0.95, eq, transform=ax.transAxes, va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="0.5", alpha=0.9))
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def compute_all(directory: Path) -> Tuple[Dict, Dict[str, pd.DataFrame]]:
    exp1 = read_csv("exp1", directory)
    exp2 = read_csv("exp2", directory)
    exp3 = read_csv("exp3", directory)
    exp4 = read_csv("exp4", directory)
    exp5 = read_csv("exp5", directory)

    # Experiment 1
    L10 = float(exp1.loc[exp1["N"] == 0, "L_1(cm)"].iloc[0])
    L20 = float(exp1.loc[exp1["N"] == 0, "L_2(cm)"].iloc[0])
    exp1["x1_m"] = (exp1["L_1(cm)"] - L10) / 100.0
    exp1["x2_m"] = (exp1["L_2(cm)"] - L20) / 100.0
    exp1["F_N"] = exp1["M_i(g)"] / 1000.0 * G
    k1, b1, r2_1 = linfit(exp1["x1_m"], exp1["F_N"])
    k2, b2, r2_2 = linfit(exp1["x2_m"], exp1["F_N"])
    k_static = k1 + k2

    # Experiment 2
    exp2["T_s"] = exp2["10T(ms)"] / 1000.0 / 10.0
    exp2["T2_s2"] = exp2["T_s"] ** 2
    exp2["M_kg"] = M_T_KG + exp2["M_i(g)"] / 1000.0
    a2, b2_fit, r2_2_period = linfit(exp2["T2_s2"], exp2["M_kg"])
    k_period = 4.0 * math.pi ** 2 * a2
    m0_prime = -b2_fit

    # Experiment 3
    exp3["A_m"] = np.abs(exp3["x'(cm)"] - X0_EXP3_CM) / 100.0
    exp3["t_s"] = exp3["t(ms)"] / 1000.0
    exp3["v_m_s"] = D_M / exp3["t_s"]
    exp3["A2_m2"] = exp3["A_m"] ** 2
    exp3["v2_m2_s2"] = exp3["v_m_s"] ** 2
    c3, d3_fit, r2_3 = linfit(exp3["A2_m2"], exp3["v2_m2_s2"])
    k_velocity = c3 * (M_U_KG + m0_prime)

    # Experiment 4
    exp4["A_m"] = np.abs(exp4["x'(cm)"] - X0_EXP4_CM) / 100.0
    exp4["T_s"] = exp4["10T(ms)"] / 1000.0 / 10.0
    s4, b4_fit, r2_4 = linfit(exp4["A_m"], exp4["T_s"])
    mean_T4 = float(exp4["T_s"].mean())
    std_T4 = float(exp4["T_s"].std(ddof=1))
    range_percent_T4 = float((exp4["T_s"].max() - exp4["T_s"].min()) / mean_T4 * 100.0)

    # Experiment 5
    exp5["T_s"] = exp5["10T(ms)"] / 1000.0 / 10.0
    exp5["T2_s2"] = exp5["T_s"] ** 2
    exp5["M_kg"] = M_T_KG + exp5["M_i(g)"] / 1000.0
    a5, b5_fit, r2_5 = linfit(exp5["T2_s2"], exp5["M_kg"])
    k_incline = 4.0 * math.pi ** 2 * a5
    m0_prime_incline = -b5_fit

    support_distance_m = (X_SUPPORT2_CM - X_SUPPORT1_CM) / 100.0
    h_m = H_CM / 100.0
    sin_theta = h_m / support_distance_m
    theta_rad = math.asin(sin_theta)
    theta_deg = theta_rad * 180.0 / math.pi

    results = {
        "constants": {
            "g": G,
            "M_T_g": M_T_G,
            "M_U_g": M_U_G,
            "m0_g": M0_G,
            "D_m": D_M,
            "x0_exp3_cm": X0_EXP3_CM,
            "x0_exp4_cm": X0_EXP4_CM,
            "support_distance_m": support_distance_m,
            "h_m": h_m,
            "sin_theta": sin_theta,
            "theta_rad": theta_rad,
            "theta_deg": theta_deg,
        },
        "exp1": {
            "L10_cm": L10, "L20_cm": L20,
            "k1": k1, "b1": b1, "r2_spring1": r2_1,
            "k2": k2, "b2": b2, "r2_spring2": r2_2,
            "k_static": k_static,
        },
        "exp2": {
            "slope": a2, "intercept": b2_fit, "r2": r2_2_period,
            "k_period": k_period, "m0_prime_kg": m0_prime, "m0_prime_g": m0_prime * 1000.0,
        },
        "exp3": {
            "slope": c3, "intercept": d3_fit, "r2": r2_3,
            "k_velocity": k_velocity,
        },
        "exp4": {
            "slope": s4, "intercept": b4_fit, "r2": r2_4,
            "mean_T": mean_T4, "std_T": std_T4, "range_percent": range_percent_T4,
        },
        "exp5": {
            "slope": a5, "intercept": b5_fit, "r2": r2_5,
            "k_incline": k_incline,
            "m0_prime_incline_kg": m0_prime_incline,
            "m0_prime_incline_g": m0_prime_incline * 1000.0,
        },
    }

    results["errors"] = {
        "delta_period_percent": abs(k_period - k_static) / k_static * 100.0,
        "delta_velocity_percent": abs(k_velocity - k_static) / k_static * 100.0,
        "delta_incline_percent": abs(k_incline - k_static) / k_static * 100.0,
        "m0_prime_to_m0_percent": m0_prime / M0_KG * 100.0,
        "incline_vs_period_percent": abs(k_incline - k_period) / k_period * 100.0,
    }

    frames = {"exp1": exp1, "exp2": exp2, "exp3": exp3, "exp4": exp4, "exp5": exp5}
    return results, frames


def save_figures(directory: Path, results: Dict, frames: Dict[str, pd.DataFrame]) -> None:
    exp1 = frames["exp1"]
    make_fit_plot(exp1["x1_m"], exp1["F_N"], results["exp1"]["k1"], results["exp1"]["b1"],
                  results["exp1"]["r2_spring1"], "伸长量 x1 / m", "拉力 F / N",
                  "实验1：弹簧1 静态拉伸拟合", directory / "exp1_spring1_fit.pdf", eq_prefix="F")
    make_fit_plot(exp1["x2_m"], exp1["F_N"], results["exp1"]["k2"], results["exp1"]["b2"],
                  results["exp1"]["r2_spring2"], "伸长量 x2 / m", "拉力 F / N",
                  "实验1：弹簧2 静态拉伸拟合", directory / "exp1_spring2_fit.pdf", eq_prefix="F")

    exp2 = frames["exp2"]
    make_fit_plot(exp2["T2_s2"], exp2["M_kg"], results["exp2"]["slope"], results["exp2"]["intercept"],
                  results["exp2"]["r2"], r"周期平方 $T^2$ / s$^2$", r"质量 $M_T+M_i$ / kg",
                  "实验2：周期平方与质量关系拟合", directory / "exp2_period_mass_fit.pdf", eq_prefix="M")

    exp3 = frames["exp3"]
    make_fit_plot(exp3["A2_m2"], exp3["v2_m2_s2"], results["exp3"]["slope"], results["exp3"]["intercept"],
                  results["exp3"]["r2"], r"振幅平方 $A^2$ / m$^2$", r"最大速度平方 $v_{max}^2$ / m$^2$ s$^{-2}$",
                  "实验3：最大速度平方与振幅平方拟合", directory / "exp3_velocity_amplitude_fit.pdf", eq_prefix=r"v^2")

    exp4 = frames["exp4"]
    make_period_amplitude_plot(exp4["A_m"], exp4["T_s"], results["exp4"]["slope"], results["exp4"]["intercept"],
                               results["exp4"]["r2"], results["exp4"]["mean_T"], directory / "exp4_period_amplitude.pdf")

    exp5 = frames["exp5"]
    make_fit_plot(exp5["T2_s2"], exp5["M_kg"], results["exp5"]["slope"], results["exp5"]["intercept"],
                  results["exp5"]["r2"], r"周期平方 $T^2$ / s$^2$", r"质量 $M_T+M_i$ / kg",
                  "实验5：有倾角时周期平方与质量关系拟合", directory / "exp5_incline_period_mass_fit.pdf", eq_prefix="M")


def write_summary(directory: Path, results: Dict) -> None:
    md = f"""# 简谐振动实验结果汇总

## 主要常量

| 物理量 | 数值 |
|---|---:|
| 当地重力加速度 $g$ | {G:.3f} m/s² |
| $M_T$ | {M_T_G:.2f} g |
| $M_U$ | {M_U_G:.2f} g |
| 两根弹簧总质量 $m_0$ | {M0_G:.2f} g |
| U 型挡光片有效宽度 $d$ | {D_M:.5f} m |
| 倾角 $\\theta$ | {results['constants']['theta_rad']:.6f} rad = {results['constants']['theta_deg']:.3f}° |

## 拟合与计算结果

| 项目 | 拟合式 | $R^2$ | 结果 |
|---|---|---:|---:|
| 弹簧 1 静态法 | $F={results['exp1']['k1']:.6f}x{results['exp1']['b1']:+.6e}$ | {results['exp1']['r2_spring1']:.6f} | $k_1={results['exp1']['k1']:.6f}$ N/m |
| 弹簧 2 静态法 | $F={results['exp1']['k2']:.6f}x{results['exp1']['b2']:+.6e}$ | {results['exp1']['r2_spring2']:.6f} | $k_2={results['exp1']['k2']:.6f}$ N/m |
| 静态法等效劲度系数 | - | - | $k_{{static}}={results['exp1']['k_static']:.6f}$ N/m |
| 周期-质量法 | $M={results['exp2']['slope']:.6f}T^2{results['exp2']['intercept']:+.6f}$ | {results['exp2']['r2']:.6f} | $k_{{period}}={results['exp2']['k_period']:.6f}$ N/m, $m_0'={results['exp2']['m0_prime_g']:.3f}$ g |
| 速度-振幅法 | $v_{{max}}^2={results['exp3']['slope']:.6f}A^2{results['exp3']['intercept']:+.6f}$ | {results['exp3']['r2']:.6f} | $k_{{velocity}}={results['exp3']['k_velocity']:.6f}$ N/m |
| 周期-振幅关系 | $T={results['exp4']['slope']:.6f}A{results['exp4']['intercept']:+.6f}$ | {results['exp4']['r2']:.6f} | $\\bar T={results['exp4']['mean_T']:.6f}$ s, 极差占比 {results['exp4']['range_percent']:.4f}% |
| 有倾角周期-质量法 | $M={results['exp5']['slope']:.6f}T^2{results['exp5']['intercept']:+.6f}$ | {results['exp5']['r2']:.6f} | $k_{{incline}}={results['exp5']['k_incline']:.6f}$ N/m, $m_{{0,incline}}'={results['exp5']['m0_prime_incline_g']:.3f}$ g |

## 相对误差（以静态法为参考）

| 比较对象 | 相对误差 |
|---|---:|
| $k_{{period}}$ 相对 $k_{{static}}$ | {results['errors']['delta_period_percent']:.3f}% |
| $k_{{velocity}}$ 相对 $k_{{static}}$ | {results['errors']['delta_velocity_percent']:.3f}% |
| $k_{{incline}}$ 相对 $k_{{static}}$ | {results['errors']['delta_incline_percent']:.3f}% |
| $k_{{incline}}$ 相对 $k_{{period}}$ | {results['errors']['incline_vs_period_percent']:.3f}% |

## 简要结论

静态拉伸法、周期-质量法和速度-振幅法得到的等效劲度系数均约为 $4.8\\sim4.9$ N/m，量级一致。实验 4 中周期随振幅变化的极差仅为 {results['exp4']['range_percent']:.4f}%，说明在实验振幅范围内周期基本不随振幅变化。实验 5 中倾角约为 {results['constants']['theta_deg']:.3f}°，倾角主要改变平衡位置，对周期-质量关系斜率影响很小。
"""
    (directory / "summary.md").write_text(md, encoding="utf-8")


def main() -> None:
    directory = Path(__file__).resolve().parent
    setup_matplotlib()
    results, frames = compute_all(directory)
    save_figures(directory, results, frames)
    write_summary(directory, results)
    # Save processed tables and complete numeric results.
    with (directory / "results.json").open("w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    for name, frame in frames.items():
        frame.to_csv(directory / f"{name}_processed.csv", index=False, encoding="utf-8-sig")
    print(json.dumps({
        "k_static": results["exp1"]["k_static"],
        "k_period": results["exp2"]["k_period"],
        "k_velocity": results["exp3"]["k_velocity"],
        "k_incline": results["exp5"]["k_incline"],
        "m0_prime_g": results["exp2"]["m0_prime_g"],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
