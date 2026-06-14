#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process the I-V data of the solar-cell experiment.

Outputs:
1. Summary table of photovoltaic parameters.
2. Grouped I-V and P-V figures.
3. Individual figures for each complex connection.
4. LaTeX table files containing the original data.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

CELL_LENGTH_CM = 5.45
CELL_WIDTH_CM = 5.95
CELL_AREA_M2 = CELL_LENGTH_CM * CELL_WIDTH_CM * 1e-4

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
FIGURE_DIR = BASE_DIR / "figures"
FIGURE_DIR.mkdir(exist_ok=True)

FILE_ORDER = [
    "single-0365W-per-m^2.csv",
    "single-0697W-per-m^2.csv",
    "single-1173W-per-m^2.csv",
    "serial-749W+653W-per-m^2.csv",
    "parallel-757W+653W-per-m^2.csv",
    "complex1s3p-352W+644W+268W+755W-per-m^2.csv",
    "complex2p2s-752W+266W+652W+356W-per-m^2.csv",
]

DISPLAY_NAMES = {
    "single-0365W-per-m^2.csv": "single, 365 W/m²",
    "single-0697W-per-m^2.csv": "single, 697 W/m²",
    "single-1173W-per-m^2.csv": "single, 1173 W/m²",
    "serial-749W+653W-per-m^2.csv": "series, 749+653 W/m²",
    "parallel-757W+653W-per-m^2.csv": "parallel, 757+653 W/m²",
    "complex1s3p-352W+644W+268W+755W-per-m^2.csv": "complex 1S+3P",
    "complex2p2s-752W+266W+652W+356W-per-m^2.csv": "complex 2S2P",
}

DISPLAY_NAMES_CN = {
    "single-0365W-per-m^2.csv": "单块电池, 365 W/m²",
    "single-0697W-per-m^2.csv": "单块电池, 697 W/m²",
    "single-1173W-per-m^2.csv": "单块电池, 1173 W/m²",
    "serial-749W+653W-per-m^2.csv": "两块串联, 749+653 W/m²",
    "parallel-757W+653W-per-m^2.csv": "两块并联, 757+653 W/m²",
    "complex1s3p-352W+644W+268W+755W-per-m^2.csv": "复杂连接 1S+3P, 352+644+268+755 W/m²",
    "complex2p2s-752W+266W+652W+356W-per-m^2.csv": "复杂连接 2S2P, 752+266+652+356 W/m²",
}

GROUP_NON_COMPLEX = [
    "single-0365W-per-m^2.csv",
    "single-0697W-per-m^2.csv",
    "single-1173W-per-m^2.csv",
    "serial-749W+653W-per-m^2.csv",
    "parallel-757W+653W-per-m^2.csv",
]

GROUP_COMPLEX = [
    "complex1s3p-352W+644W+268W+755W-per-m^2.csv",
    "complex2p2s-752W+266W+652W+356W-per-m^2.csv",
]


def latex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def parse_irradiances(filename: str) -> list[float]:
    return [float(x) for x in re.findall(r"(\d+)W", filename)]


def load_curve(filename: str) -> pd.DataFrame:
    df = pd.read_csv(DATA_DIR / filename)
    df = df.rename(columns={"U(V)": "U_V", "I(mA)": "I_mA"})
    df["P_mW"] = df["U_V"] * df["I_mA"]
    return df


def sorted_unique_curve(df: pd.DataFrame) -> pd.DataFrame:
    out = df.groupby("U_V", as_index=False)["I_mA"].mean().sort_values("U_V")
    out["P_mW"] = out["U_V"] * out["I_mA"]
    return out


def smooth_curve(df: pd.DataFrame, points: int = 500) -> pd.DataFrame:
    clean = sorted_unique_curve(df)
    if len(clean) < 3:
        return clean
    u = clean["U_V"].to_numpy()
    i = clean["I_mA"].to_numpy()
    grid = np.linspace(float(u.min()), float(u.max()), points)
    interpolator = PchipInterpolator(u, i)
    smooth_i = np.maximum(interpolator(grid), 0.0)
    return pd.DataFrame({"U_V": grid, "I_mA": smooth_i, "P_mW": grid * smooth_i})


def analyze_file(filename: str) -> dict[str, float | str | int]:
    df = load_curve(filename)
    irradiances = parse_irradiances(filename)
    input_power_w = CELL_AREA_M2 * sum(irradiances)

    isc_row = df.loc[df["U_V"].abs().idxmin()]
    voc_row = df.loc[df["I_mA"].abs().idxmin()]
    max_row = df.loc[df["P_mW"].idxmax()]

    voc_v = float(voc_row["U_V"])
    isc_mA = float(isc_row["I_mA"])
    vm_v = float(max_row["U_V"])
    im_mA = float(max_row["I_mA"])
    pm_mW = float(max_row["P_mW"])
    ff = pm_mW / (voc_v * isc_mA) if voc_v * isc_mA != 0 else np.nan
    eta_percent = (pm_mW * 1e-3) / input_power_w * 100 if input_power_w != 0 else np.nan

    return {
        "file": filename,
        "condition": DISPLAY_NAMES_CN[filename],
        "irradiance_sum_W_m2": sum(irradiances),
        "cell_count": len(irradiances),
        "input_power_W": input_power_w,
        "Voc_V": voc_v,
        "Isc_mA": isc_mA,
        "Vm_V": vm_v,
        "Im_mA": im_mA,
        "Pm_mW": pm_mW,
        "FF": ff,
        "eta_percent": eta_percent,
        "points": len(df),
    }


def plot_group(filenames: Iterable[str], kind: str, output_name: str, title: str) -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Tinos", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
    })
    plt.figure(figsize=(7.4, 4.8))
    for filename in filenames:
        df = load_curve(filename)
        smooth = smooth_curve(df)
        if kind == "iv":
            raw_y = df["I_mA"]
            smooth_y = smooth["I_mA"]
            ylabel = "Current I / mA"
        else:
            raw_y = df["P_mW"]
            smooth_y = smooth["P_mW"]
            ylabel = "Power P / mW"
        plt.plot(smooth["U_V"], smooth_y, linewidth=1.7, label=DISPLAY_NAMES[filename])
        plt.scatter(df["U_V"], raw_y, s=14, alpha=0.78)
    plt.title(title)
    plt.xlabel("Voltage U / V")
    plt.ylabel(ylabel)
    plt.grid(True, linewidth=0.45, alpha=0.55)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / f"{output_name}.pdf")
    plt.savefig(FIGURE_DIR / f"{output_name}.png", dpi=260)
    plt.close()


def save_summary_tables(summary: pd.DataFrame) -> None:
    summary.to_csv(BASE_DIR / "summary_results.csv", index=False)

    table = summary[["condition", "Voc_V", "Isc_mA", "Vm_V", "Im_mA", "Pm_mW", "FF", "eta_percent"]].copy()
    for col in ["Voc_V", "Isc_mA", "Vm_V", "Im_mA", "Pm_mW", "eta_percent"]:
        table[col] = table[col].map(lambda x: f"{x:.2f}")
    table["FF"] = table["FF"].map(lambda x: f"{x:.3f}")
    table.columns = ["实验条件", "$V_{oc}$/V", "$I_{sc}$/mA", "$V_m$/V", "$I_m$/mA", "$P_m$/mW", "$FF$", "$\\eta$/\\%"]
    latex = table.to_latex(index=False, escape=False, column_format="lrrrrrrr")
    (BASE_DIR / "summary_table.tex").write_text(latex, encoding="utf-8")


def raw_data_table_for_file(filename: str, standalone: bool = False) -> str:
    df = load_curve(filename)
    lines = []
    caption = f"{DISPLAY_NAMES_CN[filename]} 的原始测量数据"
    filename_text = latex_escape(filename)
    lines.append(r"\begin{longtable}{rrrr}")
    lines.append(rf"\caption{{{caption}. 数据文件: \texttt{{{filename_text}}}.}}\\")
    lines.append(r"\toprule")
    lines.append(r"序号 & $U$/V & $I$/mA & $P=UI$/mW \\")
    lines.append(r"\midrule")
    lines.append(r"\endfirsthead")
    lines.append(r"\toprule")
    lines.append(r"序号 & $U$/V & $I$/mA & $P=UI$/mW \\")
    lines.append(r"\midrule")
    lines.append(r"\endhead")
    lines.append(r"\midrule")
    lines.append(r"\multicolumn{4}{r}{续下页} \\")
    lines.append(r"\endfoot")
    lines.append(r"\bottomrule")
    lines.append(r"\endlastfoot")
    for idx, row in df.reset_index(drop=True).iterrows():
        lines.append(f"{idx + 1} & {row['U_V']:.2f} & {row['I_mA']:.2f} & {row['P_mW']:.2f} \\\\")
    lines.append(r"\end{longtable}")
    if not standalone:
        lines.append("")
    else:
        lines.append(r"\newpage")
    return "\n".join(lines)


def save_raw_data_tables() -> None:
    tables = [raw_data_table_for_file(filename, standalone=False) for filename in FILE_ORDER]
    (BASE_DIR / "raw_data_tables.tex").write_text("\n\n".join(tables), encoding="utf-8")

    standalone_tables = [raw_data_table_for_file(filename, standalone=True) for filename in FILE_ORDER]
    # The final \newpage is harmless but visually creates no extra page in most engines when followed by \end{document}.
    body = "\n\n".join(standalone_tables)
    raw_tex = rf"""
\documentclass[UTF8,a4paper,11pt]{{ctexart}}
\usepackage{{geometry}}
\geometry{{left=2.1cm,right=2.1cm,top=2.1cm,bottom=2.1cm}}
\usepackage{{fontspec}}
\setmainfont{{Tinos}}
\setCJKmainfont{{Noto Serif CJK SC}}
\usepackage{{booktabs,longtable,array}}
\usepackage{{siunitx}}
\usepackage[hidelinks]{{hyperref}}
\sisetup{{detect-all=true}}
\title{{太阳电池伏安特性曲线测量实验原始数据表}}
\author{{Li Jinshuo}}
\date{{\today}}
\begin{{document}}
\maketitle
\tableofcontents
\newpage
\section*{{说明}}
本文件整理所有 CSV 文件中的原始 $U-I$ 测量数据. 表中 $P=UI$ 一列由原始电压和电流逐点计算得到, 用于作图和寻找最大功率点; 其中电压单位为 V, 电流单位为 mA, 因而功率单位为 mW.
\addcontentsline{{toc}}{{section}}{{说明}}
\newpage
{body}
\end{{document}}
"""
    (BASE_DIR / "raw_data.tex").write_text(raw_tex, encoding="utf-8")


def make_figures() -> None:
    plot_group(GROUP_NON_COMPLEX, "iv", "non_complex_iv", "I-V curves: single cells, series and parallel")
    plot_group(GROUP_NON_COMPLEX, "pv", "non_complex_pv", "P-V curves: single cells, series and parallel")
    plot_group(GROUP_COMPLEX, "iv", "complex_iv", "I-V curves: complex connections")
    plot_group(GROUP_COMPLEX, "pv", "complex_pv", "P-V curves: complex connections")

    # Extra individual figures for the two complex CSV files. These are included to make the two
    # complex cases printable and reviewable independently.
    for filename in GROUP_COMPLEX:
        stem = "complex1s3p" if filename.startswith("complex1") else "complex2p2s"
        plot_group([filename], "iv", f"{stem}_iv", f"I-V curve: {DISPLAY_NAMES[filename]}")
        plot_group([filename], "pv", f"{stem}_pv", f"P-V curve: {DISPLAY_NAMES[filename]}")


def main() -> None:
    summary = pd.DataFrame([analyze_file(filename) for filename in FILE_ORDER])
    save_summary_tables(summary)
    save_raw_data_tables()
    make_figures()
    print(summary.to_string(index=False))
    print(f"\nCell area = {CELL_AREA_M2:.8f} m^2")


if __name__ == "__main__":
    main()
