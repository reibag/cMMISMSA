#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MMISMSA Analyzer - Complete Single-Script PyMOL Plugin
======================================================
An advanced, production-ready PyMOL plugin for parsing, analyzing, 
visualizing, and reporting MMISMSA molecular dynamics output files.
DATE: 2026-08-05
AUTHOR: Javier García Marín (assisted by Gemini Flash)
AFFILIATION: University of Alcalá (Department of Organic and Inorganic Chemistry)
LICENSE: GNU GENERAL PUBLIC LICENSE (GPLv3)
==============================================================================
"""

import base64
import csv
import logging
import queue
import shlex
import shutil
import subprocess
import sys
import threading
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

# PyMOL API imports
from pymol import cmd, plugins

# Configure professional logging
logger = logging.getLogger("MMISMSAAnalyzer")
logger.setLevel(logging.INFO)
if not logger.handlers:
    # Direct log output to standard stdout for PyMOL terminal inspection
    handler = logging.StreamHandler(sys.stdout)
    # Define a clean, standardized format: Timestamp [Level] LoggerName: Message
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s"))
    logger.addHandler(handler)

# ==============================================================================
# 1. DYNAMIC QT COMPATIBILITY BINDING (PySide6 / PyQt5)
# ==============================================================================
QT_BINDING = None
try:
    from PySide6 import QtCore, QtWidgets
    QT_BINDING = "PySide6"
    Slot = QtCore.Slot
except ImportError:
    try:
        from PyQt5 import QtCore, QtWidgets
        QT_BINDING = "PyQt5"
        Slot = QtCore.pyqtSlot
    except ImportError as err:
        logger.critical("Neither PySide6 nor PyQt5 could be found!")
        raise ImportError("MMISMSA Analyzer requires PySide6 or PyQt5.") from err

# Matplotlib Qt Setup (optional)
MPL_AVAILABLE = True
MPL_IMPORT_ERROR = ""
try:
    import matplotlib
    matplotlib.use("Qt5Agg")
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
except ImportError as err:
    MPL_AVAILABLE = False
    MPL_IMPORT_ERROR = str(err)
    logger.warning("matplotlib is not installed. Plotting features will be disabled.")

# ==============================================================================
# 2. CONSTANTS, LOGO & DATA MODELS
# ==============================================================================
VERSION = "1.2.0"

# Neutral UI palette to keep the interface modern and visually consistent
UI_COLORS: Dict[str, str] = {
    "bg": "#f4f6f8",
    "panel": "#ffffff",
    "panel_alt": "#eef1f4",
    "border": "#cfd6de",
    "text": "#2f3a45",
    "muted": "#6f7c89",
    "accent": "#5f7384",
    "accent_hover": "#6e8293",
    "accent_pressed": "#536576",
}

CHART_COLORS: Dict[str, str] = {
    "bar": "#5f7384",
    "bar_edge": "#4a5764",
    "error": "#8b7b6f",
    "line": "#6b7c8c",
    "neg": "#6f8f7a",
    "pos": "#9a7a7a",
    "grid": "#cfd6de",
    "axis": "#5f6b77",
    "figure_bg": "#ffffff",
    "axes_bg": "#f8fafb",
}

# Placeholder Base64 image (Transparent 1x1 PNG or embed your custom logo string here)
# Personal Logo Data (encoded in Base64)
PERSONAL_LOGO_BASE64 = b"""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
"""
ENERGY_TERMS: List[str] = [
    "Vdw", 
    "Colomb", 
    "DesolvLigand", 
    "DesolvReceptor", 
    "Apolar", 
    "Total"
]

RESIDUE_TERMS: List[str] = [
    "Total",
    "Colomb",
    "Vdw",
    "Desolv"
]

ENERGY_TERMS2: List[str] = [
    "Van der Waals", 
    "Coulombic", 
    "Ligand desolvation", 
    "Receptor desolvation", 
    "Apolar", 
    "Total"
]

RESIDUE_TERMS2: List[str] = [
    "Total",
    "Coulombic",
    "Van der Waals",
    "Desolvation"
]

RESIDUE_TERM_KEY_ALIASES: Dict[str, str] = {
    "Total": "Total energy",
    "Coulombic": "Coulombic energy",
    "Van der Waals": "Van der Waals energy",
    "Desolvation": "Desolvation energy",
}

# Bi-directional term translation maps
TERM_MAPPING: Dict[str, str] = dict(zip(ENERGY_TERMS, ENERGY_TERMS2))

# Mapping palette names to PyMOL ramp color definitions
PALETTE_RAMP_MAP: Dict[str, Union[str, List[str]]] = {
    "bwr": ["blue", "white", "red"],
    "rainbow": "rainbow",
    "red_white_blue": ["red", "white", "blue"],
    "coolwarm": ["blue", "white", "red"]
}

@dataclass(frozen=True)
class EnergyTermStats:
    """Dataclass storing summary metrics for an energy term."""
    term: str
    mean: float
    std: float


class SortableResidueTable(QtWidgets.QTableWidget):
    """QTableWidget that sorts numeric data correctly when a header is clicked."""

    def sortItems(self, column: int, order: QtCore.Qt.SortOrder = QtCore.Qt.AscendingOrder) -> None:
        rows = self.rowCount()
        items = []
        for row in range(rows):
            item = self.item(row, column)
            if item is None:
                continue
            raw_value = item.data(QtCore.Qt.UserRole)
            if raw_value is None:
                try:
                    raw_value = float(item.text())
                except (TypeError, ValueError):
                    raw_value = item.text()
            items.append((raw_value, row, item))

        def compare_values(left, right):
            left_value, _, _ = left
            right_value, _, _ = right
            if isinstance(left_value, (int, float)) and isinstance(right_value, (int, float)):
                return -1 if left_value < right_value else 1 if left_value > right_value else 0
            if isinstance(left_value, str) and isinstance(right_value, str):
                return -1 if left_value.lower() < right_value.lower() else 1 if left_value.lower() > right_value.lower() else 0
            return 0

        items.sort(key=lambda entry: entry[0], reverse=(order == QtCore.Qt.DescendingOrder))
        if not all(isinstance(entry[0], (int, float)) for entry in items):
            items.sort(key=lambda entry: str(entry[0]).lower(), reverse=(order == QtCore.Qt.DescendingOrder))

        for new_row, (_, old_row, _) in enumerate(items):
            for col in range(self.columnCount()):
                source_item = self.takeItem(old_row, col)
                if source_item is not None:
                    self.setItem(new_row, col, source_item)

        self.sortByColumn(column, order)


class MMISMSADataSet:
    """Stores trajectory arrays and calculates NumPy summary statistics."""

    def __init__(self, frames: np.ndarray, data: Dict[str, np.ndarray]) -> None:
        self.frames = frames
        # Map raw column names to descriptive ENERGY_TERMS2 labels
        self.data: Dict[str, np.ndarray] = {}
        for raw_k, arr in data.items():
            descriptive_k = TERM_MAPPING.get(raw_k, raw_k)
            self.data[descriptive_k] = arr

        self._stats: Dict[str, EnergyTermStats] = {}
        self._compute_all_statistics()

    def _compute_all_statistics(self) -> None:
        """Vectorized computation of mean and sample standard deviation."""
        for term in ENERGY_TERMS2:
            if term in self.data:
                arr = self.data[term]
                mean_val = float(np.mean(arr))
                std_val = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
                self._stats[term] = EnergyTermStats(term=term, mean=mean_val, std=std_val)

    def get_stats(self, term: str) -> EnergyTermStats:
        if term not in self._stats:
            raise KeyError(f"Term '{term}' not found in calculated statistics.")
        return self._stats[term]

    def get_all_stats(self) -> Dict[str, EnergyTermStats]:
        return self._stats

    @property
    def frame_count(self) -> int:
        return len(self.frames)


def resolve_residue_term_key(term: str, data_dict: Dict[str, np.ndarray]) -> str:
    """Return a valid residue term key from UI label or stored descriptive key."""
    if term in data_dict:
        return term

    mapped = RESIDUE_TERM_KEY_ALIASES.get(term)
    if mapped and mapped in data_dict:
        return mapped

    raise KeyError(term)


def ensure_matplotlib_available(parent: QtWidgets.QWidget, action: str) -> bool:
    """Show a user-facing error if matplotlib is missing and return availability."""
    if MPL_AVAILABLE:
        return True

    details = f"\n\nTechnical detail: {MPL_IMPORT_ERROR}" if MPL_IMPORT_ERROR else ""
    QtWidgets.QMessageBox.critical(
        parent,
        "Matplotlib not installed",
        f"Cannot {action} because matplotlib is not installed.{details}",
    )
    return False

# ==============================================================================
# 3. CSV PARSER & EXPORTERS
# ==============================================================================
class MMISMSAParserError(Exception):
    """Custom parser exception."""
    pass


class MMISMSACSVReader:
    """Robust parser capable of handling global CSVs and per-residue CSV files."""

    @classmethod
    def load_file(cls, filepath: Path) -> MMISMSADataSet:
        path = Path(filepath)
        if not path.exists():
            raise MMISMSAParserError(f"File not found: {path}")

        frames: List[int] = []
        raw_columns: Dict[str, List[float]] = {term: [] for term in ENERGY_TERMS}

        try:
            with open(path, "r", encoding="utf-8") as f:
                lines = [line.strip() for line in f if line.strip()]

            if not lines:
                raise MMISMSAParserError("Empty CSV file encountered.")

            cls._validate_headers(lines[0])

            for line_idx, line in enumerate(lines[1:], start=2):
                parts = [p.strip() for p in line.split(";") if p.strip()]
                
                if len(parts) > 0 and "," in parts[0]:
                    sub_parts = parts[0].split(",")
                    parts = sub_parts + parts[1:]

                if len(parts) < 7:
                    logger.warning("Line %d has insufficient columns (%d), skipping.", line_idx, len(parts))
                    continue

                try:
                    frame_num = int(parts[0])
                    vdw_val = float(parts[1])
                    col_val = float(parts[2])
                    dl_val = float(parts[3])
                    dr_val = float(parts[4])
                    apol_val = float(parts[5])
                    tot_val = float(parts[6])

                    frames.append(frame_num)
                    raw_columns["Vdw"].append(vdw_val)
                    raw_columns["Colomb"].append(col_val)
                    raw_columns["DesolvLigand"].append(dl_val)
                    raw_columns["DesolvReceptor"].append(dr_val)
                    raw_columns["Apolar"].append(apol_val)
                    raw_columns["Total"].append(tot_val)
                except ValueError as ve:
                    logger.warning("Line %d contains invalid numbers (%s), skipping.", line_idx, ve)

        except Exception as err:
            raise MMISMSAParserError(f"Failed parsing file: {err}") from err

        if not frames:
            raise MMISMSAParserError("No valid data rows parsed from CSV.")

        return MMISMSADataSet(
            frames=np.array(frames, dtype=int),
            data={k: np.array(v, dtype=float) for k, v in raw_columns.items()}
        )

    @classmethod
    def load_residue_file(cls, filepath: Path) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
        """
        Parses per-residue CSV file returning dictionaries of means and stds for each descriptive term in RESIDUE_TERMS2.
        """
        path = Path(filepath)
        if not path.exists():
            raise MMISMSAParserError(f"Residue file not found: {path}")

        data_matrix = []
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line: continue
                parts = [float(x) for x in line.split(";") if x.strip()]
                if parts:
                    data_matrix.append(parts)

        if not data_matrix:
            raise MMISMSAParserError("No numeric rows found in residue file.")

        arr = np.array(data_matrix)
        n_cols = arr.shape[1]
        n_residues = n_cols // 3
        
        reshaped = arr[:, :n_residues*3].reshape(-1, n_residues, 3)
        
        coulomb_matrix = reshaped[:, :, 0]
        vdw_matrix = reshaped[:, :, 1]
        desolv_matrix = reshaped[:, :, 2]
        total_matrix = coulomb_matrix + vdw_matrix + desolv_matrix

        ddof = 1 if len(arr) > 1 else 0

        res_means = {
            "Total energy": np.mean(total_matrix, axis=0),
            "Coulombic energy": np.mean(coulomb_matrix, axis=0),
            "Van der Waals energy": np.mean(vdw_matrix, axis=0),
            "Desolvation energy": np.mean(desolv_matrix, axis=0)
        }

        res_stds = {
            "Total energy": np.std(total_matrix, axis=0, ddof=ddof),
            "Coulombic energy": np.std(coulomb_matrix, axis=0, ddof=ddof),
            "Van der Waals energy": np.std(vdw_matrix, axis=0, ddof=ddof),
            "Desolvation energy": np.std(desolv_matrix, axis=0, ddof=ddof)
        }
        
        return res_means, res_stds

    @staticmethod
    def _validate_headers(header_line: str) -> None:
        for term in ENERGY_TERMS:
            if term not in header_line:
                raise MMISMSAParserError(f"Missing mandatory energy term header: '{term}'")


class MMISMSAExporter:
    """Exports structured data files, figures, and HTML reports."""

    @staticmethod
    def export_statistics_csv(dataset: MMISMSADataSet, output_path: Path) -> None:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["Term", "Mean", "Std"])
            for term in ENERGY_TERMS2:
                st = dataset.get_stats(term)
                writer.writerow([st.term, f"{st.mean:.6f}", f"{st.std:.6f}"])

    @staticmethod
    def export_full_analysis(dataset: MMISMSADataSet, output_path: Path) -> None:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["Frame"] + ENERGY_TERMS2)
            for idx, frame in enumerate(dataset.frames):
                writer.writerow([frame] + [dataset.data[term][idx] for term in ENERGY_TERMS2])

            writer.writerow([])
            writer.writerow(["# SUMMARY"])
            writer.writerow(["Term", "Mean", "Std"])
            for term in ENERGY_TERMS2:
                st = dataset.get_stats(term)
                writer.writerow([st.term, f"{st.mean:.6f}", f"{st.std:.6f}"])

    @staticmethod
    def export_residue_statistics_csv(
        res_means: Dict[str, np.ndarray], 
        res_stds: Dict[str, np.ndarray], 
        output_path: Path,
        start_res: int = 1
    ) -> None:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        n_res = len(res_means["Total energy"])

        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow([
                "Residue", 
                "Total_Energy_Mean", "Total_Energy_Std",
                "Coulombic_Energy_Mean", "Coulombic_Energy_Std",
                "Vdw_Energy_Mean", "Vdw_Energy_Std",
                "Desolvation_Energy_Mean", "Desolvation_Energy_Std"
            ])
            for idx in range(n_res):
                res_num = start_res + idx
                writer.writerow([
                    f"Residue_{res_num}",
                    f"{res_means['Total energy'][idx]:.6f}", f"{res_stds['Total energy'][idx]:.6f}",
                    f"{res_means['Coulombic energy'][idx]:.6f}", f"{res_stds['Coulombic energy'][idx]:.6f}",
                    f"{res_means['Van der Waals energy'][idx]:.6f}", f"{res_stds['Van der Waals energy'][idx]:.6f}",
                    f"{res_means['Desolvation energy'][idx]:.6f}", f"{res_stds['Desolvation energy'][idx]:.6f}"
                ])

    @staticmethod
    def generate_html_report(
        dataset: MMISMSADataSet,
        source_filename: str,
        output_dir: Path,
        chart_bytes: Optional[bytes] = None
    ) -> Path:
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        MMISMSAExporter.export_statistics_csv(dataset, out_dir / "statistics.csv")
        MMISMSAExporter.export_full_analysis(dataset, out_dir / "full_analysis.csv")

        def _to_b64(b: Optional[bytes]) -> str:
            return f"data:image/png;base64,{base64.b64encode(b).decode('utf-8')}" if b else ""

        rows = "".join([
            f"<tr><td><b>{t}</b></td><td>{dataset.get_stats(t).mean:.4f}</td><td>{dataset.get_stats(t).std:.4f}</td></tr>\n"
            for t in ENERGY_TERMS2
        ])

        html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>MMISMSA Analysis Report</title>
    <style>
        :root {{
            --bg: {UI_COLORS['bg']};
            --panel: {UI_COLORS['panel']};
            --panel-alt: {UI_COLORS['panel_alt']};
            --border: {UI_COLORS['border']};
            --text: {UI_COLORS['text']};
            --muted: {UI_COLORS['muted']};
            --accent: {UI_COLORS['accent']};
            --chart-bg: {CHART_COLORS['axes_bg']};
        }}
        body {{
            font-family: "Segoe UI", "Noto Sans", "Helvetica Neue", Arial, sans-serif;
            margin: 25px;
            background: var(--bg);
            color: var(--text);
            line-height: 1.45;
        }}
        h1, h2 {{
            color: var(--text);
            border-bottom: 1px solid var(--border);
            padding-bottom: 7px;
            margin-bottom: 14px;
        }}
        .box {{
            background: var(--panel-alt);
            border: 1px solid var(--border);
            padding: 14px;
            border-radius: 10px;
            margin-bottom: 20px;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            background: var(--panel);
            border: 1px solid var(--border);
            border-radius: 10px;
            overflow: hidden;
        }}
        th, td {{
            padding: 10px;
            border: 1px solid #dde3e8;
            text-align: left;
        }}
        th {{
            background: var(--panel-alt);
            color: var(--text);
            font-weight: 600;
        }}
        tr:nth-child(even) td {{
            background: #f8fafb;
        }}
        .flex {{ display: flex; gap: 20px; flex-wrap: wrap; margin-top: 15px; }}
        .card {{
            background: var(--panel);
            padding: 15px;
            border: 1px solid var(--border);
            border-radius: 10px;
            flex: 1;
            min-width: 450px;
            text-align: center;
        }}
        .card h3 {{
            color: var(--muted);
            margin-top: 2px;
            margin-bottom: 12px;
            font-weight: 600;
        }}
        .chart-wrap {{
            background: var(--chart-bg);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 8px;
        }}
        img {{ max-width: 100%; height: auto; }}
    </style>
</head>
<body>
    <h1>MMISMSA Execution Report</h1>
    <div class="box">
        <p><strong>File:</strong> {source_filename}</p>
        <p><strong>Timestamp:</strong> {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
        <p><strong>Frames Analyzed:</strong> {dataset.frame_count}</p>
    </div>
    <h2>Energy Statistics</h2>
    <table>
        <thead><tr><th>Term</th><th>Mean (kcal/mol)</th><th>Std Dev (kcal/mol)</th></tr></thead>
        <tbody>{rows}</tbody>
    </table>
    <h2>Global Energy Means & Error Bars</h2>
    <div class="flex">
        {f'<div class="card"><h3>Energy Profile (Mean ± SD)</h3><div class="chart-wrap"><img src="{_to_b64(chart_bytes)}"/></div></div>' if chart_bytes else ''}
    </div>
</body>
</html>"""
        report_path = out_dir / "report.html"
        report_path.write_text(html, encoding="utf-8")
        return report_path

# ==============================================================================
# 4. EMBEDDED WIDGETS
# ==============================================================================
class MatplotlibWidget(QtWidgets.QWidget):
    """Matplotlib embedded canvas for bar charts with error bars."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.figure = None
        self.canvas = None
        self.ax = None

        if MPL_AVAILABLE:
            self.figure = Figure(figsize=(6, 3.2), dpi=100)
            self.canvas = FigureCanvas(self.figure)
            self.ax = self.figure.add_subplot(111)
            self._apply_chart_style()
            layout.addWidget(self.canvas)
        else:
            info = QtWidgets.QLabel(
                "Matplotlib is not installed. Charts are disabled in this environment."
            )
            info.setAlignment(QtCore.Qt.AlignCenter)
            info.setWordWrap(True)
            layout.addWidget(info)

    def _ensure_available(self) -> None:
        if not MPL_AVAILABLE or self.figure is None or self.canvas is None or self.ax is None:
            raise RuntimeError("matplotlib is not installed")

    def _apply_chart_style(self) -> None:
        """Uses an explicit chart theme so screen and exports always match."""
        self._ensure_available()
        self.figure.patch.set_facecolor(CHART_COLORS["figure_bg"])
        self.ax.set_facecolor(CHART_COLORS["axes_bg"])
        self.ax.tick_params(colors=CHART_COLORS["axis"])
        for spine in self.ax.spines.values():
            spine.set_color("#bac4cd")

    def plot_bar_chart_with_errors(
        self, 
        labels: List[str], 
        means: List[float], 
        stds: List[float], 
        title: str = "MMISMSA Energy Components (Mean ± SD)", 
        ylabel: str = "Energy (kcal/mol)"
    ) -> None:
        self._ensure_available()
        self.ax.clear()
        self._apply_chart_style()
        x = np.arange(len(labels))

        bars = self.ax.bar(
            x, means, yerr=stds, capsize=5, color=CHART_COLORS["bar"],
            edgecolor=CHART_COLORS["bar_edge"], alpha=0.9,
            error_kw={"ecolor": CHART_COLORS["error"], "linewidth": 1.2}
        )

        self.ax.set_title(title, fontsize=11, fontweight="bold")
        self.ax.set_ylabel(ylabel, fontsize=9)
        self.ax.set_xticks(x)
        self.ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
        self.ax.grid(axis="y", linestyle="--", alpha=0.7, color=CHART_COLORS["grid"])

        # Mean and std labels in barplot
        #
        #for bar, m, s in zip(bars, means, stds):
        #    h = bar.get_height()
        #    va = 'bottom' if h >= 0 else 'top'
        #    self.ax.annotate(
        #        f'{m:.2f} ± {s:.2f}', 
        #        xy=(bar.get_x() + bar.get_width() / 2, h + (s if h >= 0 else -s)),
        #        xytext=(0, 3 if h >= 0 else -3), 
        #        textcoords="offset points",
        #        ha='center', va=va, fontsize=7.5, fontweight='bold'
        #    )

        self.figure.tight_layout()
        self.canvas.draw()

    def plot_residue_profile(self, means: np.ndarray, stds: np.ndarray, term_name: str = "Total energy", start_res: int = 1) -> None:
        self._ensure_available()
        self.ax.clear()
        self._apply_chart_style()
        x = np.arange(start_res, start_res + len(means))

        colors = [CHART_COLORS["neg"] if m < 0 else CHART_COLORS["pos"] for m in means]

        self.ax.bar(x, means, yerr=stds, capsize=2, color=colors, alpha=0.9, ecolor='#8a9299', linewidth=0.5)
        self.ax.set_title(f"Per-Residue Energy Profile ({term_name}) [Mean ± SD]", fontsize=11, fontweight="bold")
        self.ax.set_xlabel("Residue Index", fontsize=9)
        self.ax.set_ylabel("Energy (kcal/mol)", fontsize=9)
        self.ax.grid(axis="y", linestyle="--", alpha=0.7, color=CHART_COLORS["grid"])
        self.ax.axhline(0, color='black', linewidth=0.8, linestyle='-')

        self.figure.tight_layout()
        self.canvas.draw()

    def get_png_bytes(self, dpi: int = 300) -> bytes:
        self._ensure_available()
        import io
        buf = io.BytesIO()
        self.figure.savefig(
            buf,
            format="png",
            dpi=dpi,
            facecolor=self.figure.get_facecolor(),
            edgecolor="none",
        )
        buf.seek(0)
        return buf.read()


class TrajectoryWidget(QtWidgets.QWidget):
    """Interactive trajectory plotter using Matplotlib when available."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        if MPL_AVAILABLE:
            self.mpl = MatplotlibWidget()
            layout.addWidget(self.mpl)
        else:
            self.mpl = None
            info = QtWidgets.QLabel(
                "Matplotlib is not installed. Trajectory plots are unavailable."
            )
            info.setAlignment(QtCore.Qt.AlignCenter)
            info.setWordWrap(True)
            layout.addWidget(info)

    def plot_trajectory(self, frames: np.ndarray, energies: np.ndarray, term: str) -> None:
        if not MPL_AVAILABLE or self.mpl is None:
            raise RuntimeError("matplotlib is not installed")

        self.mpl.ax.clear()
        self.mpl._apply_chart_style()
        self.mpl.ax.plot(frames, energies, color=CHART_COLORS["line"], linewidth=1.5)
        self.mpl.ax.set_title(f"Energy Trajectory: {term}")
        self.mpl.ax.set_xlabel("Frame")
        self.mpl.ax.set_ylabel("Energy (kcal/mol)")
        self.mpl.ax.grid(True, linestyle="--", alpha=0.7, color=CHART_COLORS["grid"])
        self.mpl.figure.tight_layout()
        self.mpl.canvas.draw()

# ==============================================================================
# 5. UI TABS & CONTROLS
# ==============================================================================
class TrajectoryTab(QtWidgets.QWidget):
    """Tab housing frame trajectory plots."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._dataset: Optional[MMISMSADataSet] = None
        self._setup_ui()

    def _setup_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        ctrls = QtWidgets.QHBoxLayout()
        ctrls.addWidget(QtWidgets.QLabel("Energy Term:"))
        self.combo = QtWidgets.QComboBox()
        self.combo.addItems(ENERGY_TERMS2)
        self.combo.setCurrentText("Total energy")
        ctrls.addWidget(self.combo)

        self.btn_plot = QtWidgets.QPushButton("Plot Trajectory")
        self.btn_plot.clicked.connect(self._on_plot)
        ctrls.addWidget(self.btn_plot)
        layout.addLayout(ctrls)

        self.traj_widget = TrajectoryWidget()
        layout.addWidget(self.traj_widget)

    def set_dataset(self, dataset: MMISMSADataSet) -> None:
        self._dataset = dataset
        if MPL_AVAILABLE:
            self._on_plot()

    @Slot()
    def _on_plot(self) -> None:
        if not self._dataset: return
        if not ensure_matplotlib_available(self, "plot trajectory"):
            return
        term = self.combo.currentText()
        self.traj_widget.plot_trajectory(self._dataset.frames, self._dataset.data[term], term)


class CombinedGlobalAnalysisTab(QtWidgets.QWidget):
    """
    Fused Global Analysis Tab: combines energy bar charts with error bars,
    a summary statistics table, and export options.
    """

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._dataset: Optional[MMISMSADataSet] = None
        self._source_filepath: Optional[Path] = None
        self._load_handler = None
        self._setup_ui()

    def _setup_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        # Local load controls (kept only inside this tab by design)
        load_row = QtWidgets.QHBoxLayout()
        self.btn_load_global = QtWidgets.QPushButton("Load Global CSV")
        self.btn_load_global.setStyleSheet("font-weight: 600;")
        self.btn_load_global.clicked.connect(self._on_load_clicked)

        self.lbl_file_info = QtWidgets.QLabel("No file loaded.")
        self.lbl_file_info.setStyleSheet("color: #7f8c8d; font-style: italic;")

        load_row.addWidget(self.btn_load_global)
        load_row.addWidget(self.lbl_file_info)
        load_row.addStretch()
        layout.addLayout(load_row)

        # Matplotlib Plot Canvas
        self.canvas_widget = MatplotlibWidget()
        layout.addWidget(self.canvas_widget, stretch=3)

        # Statistics Table
        self.table = QtWidgets.QTableWidget()
        self.table.setColumnCount(3)
        self.table.setHorizontalHeaderLabels(["Energy Component", "Mean (kcal/mol)", "Std Dev (kcal/mol)"])
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        layout.addWidget(self.table, stretch=2)

        # Integrated Export Section
        grp_exports = QtWidgets.QGroupBox("Global Data & Report Exports")
        lay_exports = QtWidgets.QVBoxLayout(grp_exports)

        row1 = QtWidgets.QHBoxLayout()
        self.btn_exp_stats = QtWidgets.QPushButton("Export Statistics CSV")
        self.btn_exp_stats.clicked.connect(self._exp_stats)
        self.btn_exp_all = QtWidgets.QPushButton("Export All Results")
        self.btn_exp_all.setStyleSheet("font-weight: 600;")
        self.btn_exp_all.clicked.connect(self._exp_all)

        self.btn_gen_report = QtWidgets.QPushButton("Generate Report Bundle")
        self.btn_gen_report.setStyleSheet("font-weight: 600;")
        self.btn_gen_report.clicked.connect(self._gen_report)

        row1.addWidget(self.btn_exp_stats)
        row1.addWidget(self.btn_exp_all)
        row1.addWidget(self.btn_gen_report)

        lay_exports.addLayout(row1)

        layout.addWidget(grp_exports)

    def set_load_handler(self, handler) -> None:
        self._load_handler = handler

    def set_file_info(self, text: str) -> None:
        self.lbl_file_info.setText(text)

    @Slot()
    def _on_load_clicked(self) -> None:
        if self._load_handler is not None:
            self._load_handler()

    def set_dataset(self, dataset: MMISMSADataSet, filepath: Path) -> None:
        self._dataset = dataset
        self._source_filepath = filepath
        self.refresh_display()

    def refresh_display(self) -> None:
        if not self._dataset: return

        stats = self._dataset.get_all_stats()

        # 1. Update Chart (if matplotlib is available)
        if MPL_AVAILABLE:
            means = [stats[t].mean for t in ENERGY_TERMS2]
            stds = [stats[t].std for t in ENERGY_TERMS2]
            self.canvas_widget.plot_bar_chart_with_errors(ENERGY_TERMS2, means, stds)

        # 2. Update Table
        self.table.setRowCount(len(ENERGY_TERMS2))
        for row, term in enumerate(ENERGY_TERMS2):
            st = stats[term]
            self.table.setItem(row, 0, QtWidgets.QTableWidgetItem(st.term))
            self.table.setItem(row, 1, QtWidgets.QTableWidgetItem(f"{st.mean:.6f}"))
            self.table.setItem(row, 2, QtWidgets.QTableWidgetItem(f"{st.std:.6f}"))

    @Slot()
    def _exp_stats(self) -> None:
        if not self._check_dataset(): return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save Stats CSV", "statistics.csv", "CSV (*.csv)")
        if path: MMISMSAExporter.export_statistics_csv(self._dataset, Path(path))

    @Slot()
    def _exp_all(self) -> None:
        if not self._check_dataset(): return
        if not ensure_matplotlib_available(self, "export chart images"):
            return
        dir_p = QtWidgets.QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir_p:
            out = Path(dir_p) / "MMISMSA_Exports"
            out.mkdir(parents=True, exist_ok=True)
            self.refresh_display()
            MMISMSAExporter.export_statistics_csv(self._dataset, out / "statistics.csv")
            MMISMSAExporter.export_full_analysis(self._dataset, out / "full_analysis.csv")
            (out / "global_energy_profile.png").write_bytes(self.canvas_widget.get_png_bytes())
            QtWidgets.QMessageBox.information(
                self,
                "Success",
                f"Exported files to:\n{out}\n\nIncluded chart: global_energy_profile.png",
            )

    @Slot()
    def _gen_report(self) -> None:
        if not self._check_dataset(): return
        if not ensure_matplotlib_available(self, "generate report charts"):
            return
        dir_p = QtWidgets.QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir_p:
            out = Path(dir_p) / "MMISMSA_Report"
            self.refresh_display()
            chart_bytes = self.canvas_widget.get_png_bytes()

            source_name = self._source_filepath.name if self._source_filepath else "global.csv"
            rep = MMISMSAExporter.generate_html_report(self._dataset, source_name, out, chart_bytes)
            QtWidgets.QMessageBox.information(self, "Report Created", f"Report generated:\n{rep}")

    def _check_dataset(self) -> bool:
        if not self._dataset:
            QtWidgets.QMessageBox.warning(self, "No Data", "Please load a global CSV file first.")
            return False
        return True


class ResidueAnalysisTab(QtWidgets.QWidget):
    """Tab dedicated to per-residue energy statistics, PDB loading, custom renumbering, CSV export, and PyMOL 3D mapping with scalebar."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._res_means: Optional[Dict[str, np.ndarray]] = None
        self._res_stds: Optional[Dict[str, np.ndarray]] = None
        self._loaded_pdb_obj: Optional[str] = None
        self._setup_ui()

    def _setup_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        # File Management Group
        ctrl_box = QtWidgets.QGroupBox("Per-Residue Files & Renumbering")
        ctrl_lay = QtWidgets.QVBoxLayout(ctrl_box)

        r1 = QtWidgets.QHBoxLayout()
        self.btn_load_res = QtWidgets.QPushButton("Load Per-Residue CSV")
        self.btn_load_res.setStyleSheet("font-weight: 600;")
        self.btn_load_res.clicked.connect(self._load_res_file)

        self.btn_load_pdb = QtWidgets.QPushButton("Load MMISMSA PDB File")
        self.btn_load_pdb.setStyleSheet("font-weight: 600;")
        self.btn_load_pdb.clicked.connect(self._load_pdb_file)

        self.btn_exp_res_csv = QtWidgets.QPushButton("Export Residue Statistics CSV")
        self.btn_exp_res_csv.clicked.connect(self._export_res_csv)

        r1.addWidget(self.btn_load_res)
        r1.addWidget(self.btn_load_pdb)
        r1.addWidget(self.btn_exp_res_csv)

        r2 = QtWidgets.QHBoxLayout()
        r2.addWidget(QtWidgets.QLabel("Set Start Residue ID:"))
        self.spin_start_res = QtWidgets.QSpinBox()
        self.spin_start_res.setRange(-9999, 99999)
        self.spin_start_res.setValue(1)
        self.spin_start_res.setToolTip("AMBER/ISMSA starts residue indexing at 1. Set this to match your target PDB numbering.")
        self.spin_start_res.valueChanged.connect(self._on_start_residue_changed)
        r2.addWidget(self.spin_start_res)

        self.lbl_res_info = QtWidgets.QLabel("No residue CSV loaded.")
        self.lbl_res_info.setStyleSheet("color: #7f8c8d; font-style: italic;")
        self.lbl_pdb_info = QtWidgets.QLabel("No PDB loaded.")
        self.lbl_pdb_info.setStyleSheet("color: #7f8c8d; font-style: italic;")

        r2.addWidget(self.lbl_res_info)
        r2.addWidget(self.lbl_pdb_info)
        r2.addStretch()

        ctrl_lay.addLayout(r1)
        ctrl_lay.addLayout(r2)
        layout.addWidget(ctrl_box)

        # Table
        self.table = SortableResidueTable()
        self.table.setColumnCount(9)
        self.table.setHorizontalHeaderLabels([
            "Residue", 
            "Total Mean", "Total Std",
            "Coulombic Mean", "Coulombic Std",
            "Vdw Mean", "Vdw Std",
            "Desolvation Mean", "Desolvation Std"
        ])
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.table.setSortingEnabled(True)
        self.table.itemSelectionChanged.connect(self._on_residue_table_selection_changed)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        layout.addWidget(self.table)

        # PyMOL Mapping Controls
        map_box = QtWidgets.QGroupBox("PyMOL Structure Mapping & Scalebar")
        map_lay = QtWidgets.QHBoxLayout(map_box)

        map_lay.addWidget(QtWidgets.QLabel("Term Component:"))
        self.combo_term = QtWidgets.QComboBox()
        self.combo_term.addItems(RESIDUE_TERMS2)
        map_lay.addWidget(self.combo_term)

        map_lay.addWidget(QtWidgets.QLabel("Metric:"))
        self.combo_metric = QtWidgets.QComboBox()
        self.combo_metric.addItems(["mean", "std"])
        map_lay.addWidget(self.combo_metric)

        map_lay.addWidget(QtWidgets.QLabel("Selection Target:"))
        self.txt_sel = QtWidgets.QLineEdit("all")
        map_lay.addWidget(self.txt_sel)

        map_lay.addWidget(QtWidgets.QLabel("Palette:"))
        self.combo_pal = QtWidgets.QComboBox()
        self.combo_pal.addItems(["bwr", "rainbow", "red_white_blue", "coolwarm"])
        map_lay.addWidget(self.combo_pal)

        self.btn_map = QtWidgets.QPushButton("Map Statistics & Show Scalebar")
        self.btn_map.setStyleSheet("font-weight: 600;")
        self.btn_map.clicked.connect(self._map_pymol_residues)
        map_lay.addWidget(self.btn_map)

        layout.addWidget(map_box)

    def get_start_residue(self) -> int:
        return self.spin_start_res.value()

    @Slot()
    def _load_res_file(self) -> None:
        f_str, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open MMISMSA Per-Residue CSV File", "", "CSV Files (*.csv *.txt);;All Files (*)"
        )
        if not f_str: return

        try:
            p = Path(f_str)
            means_dict, stds_dict = MMISMSACSVReader.load_residue_file(p)
            self._res_means = means_dict
            self._res_stds = stds_dict
            n_res = len(means_dict["Total energy"])

            self.lbl_res_info.setText(f"CSV: {p.name} ({n_res} residues)")
            self._update_table_display()

            QtWidgets.QMessageBox.information(
                self, "Loaded", f"Parsed residue decomposition for {n_res} residues across all terms."
            )
        except Exception as err:
            logger.exception("Failed loading per-residue CSV")
            QtWidgets.QMessageBox.critical(self, "Error", f"Failed loading residue file:\n{err}")

    @Slot()
    def _on_start_residue_changed(self) -> None:
        if self._res_means is not None:
            self._update_table_display()

    def _update_table_display(self) -> None:
        if self._res_means is None or self._res_stds is None: return
        n_res = len(self._res_means["Total energy"])
        start_res = self.get_start_residue()

        self.table.setRowCount(n_res)
        for idx in range(n_res):
            display_res_num = start_res + idx
            original_res_num = idx + 1
            row_items = []
            residue_item = QtWidgets.QTableWidgetItem(f"Residue {display_res_num}")
            residue_item.setData(QtCore.Qt.UserRole, display_res_num)
            residue_item.setData(QtCore.Qt.UserRole + 1, original_res_num)
            row_items.append(residue_item)
            row_items.append(QtWidgets.QTableWidgetItem(f"{self._res_means['Total energy'][idx]:.4f}"))
            row_items.append(QtWidgets.QTableWidgetItem(f"{self._res_stds['Total energy'][idx]:.4f}"))
            row_items.append(QtWidgets.QTableWidgetItem(f"{self._res_means['Coulombic energy'][idx]:.4f}"))
            row_items.append(QtWidgets.QTableWidgetItem(f"{self._res_stds['Coulombic energy'][idx]:.4f}"))
            row_items.append(QtWidgets.QTableWidgetItem(f"{self._res_means['Van der Waals energy'][idx]:.4f}"))
            row_items.append(QtWidgets.QTableWidgetItem(f"{self._res_stds['Van der Waals energy'][idx]:.4f}"))
            row_items.append(QtWidgets.QTableWidgetItem(f"{self._res_means['Desolvation energy'][idx]:.4f}"))
            row_items.append(QtWidgets.QTableWidgetItem(f"{self._res_stds['Desolvation energy'][idx]:.4f}"))

            for col_idx, item in enumerate(row_items):
                if col_idx == 0:
                    item.setData(QtCore.Qt.UserRole, display_res_num)
                else:
                    try:
                        value = float(item.text())
                    except ValueError:
                        value = item.text()
                    item.setData(QtCore.Qt.UserRole, value)
                self.table.setItem(idx, col_idx, item)

    def _parse_selected_residue_number(self) -> Optional[int]:
        if self.table.currentRow() < 0:
            return None

        item = self.table.item(self.table.currentRow(), 0)
        if item is None:
            return None

        original_residue = item.data(QtCore.Qt.UserRole + 1)
        if isinstance(original_residue, (int, float)):
            return int(original_residue)

        text = item.text().strip()
        try:
            return int(text.split()[-1])
        except ValueError:
            return None

    @Slot()
    def _on_residue_table_selection_changed(self) -> None:
        if not self._loaded_pdb_obj:
            return

        res_num = self._parse_selected_residue_number()
        if res_num is None:
            return

        selection_target = (self.txt_sel.text().strip() or self._loaded_pdb_obj).strip()
        pymol_sel = f"({selection_target}) and resi {res_num}"
        cmd.select("mmismsa_selected_residue", pymol_sel)
        cmd.zoom("mmismsa_selected_residue", 5)
        cmd.center("mmismsa_selected_residue")

    @Slot()
    def _load_pdb_file(self) -> None:
        f_str, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open MMISMSA PDB Structure File", "", "PDB Files (*.pdb);;All Files (*)"
        )
        if not f_str: return

        try:
            p = Path(f_str)
            obj_name = "mmismsa_complex"
            cmd.load(str(p), obj_name)
            self._loaded_pdb_obj = obj_name
            self.txt_sel.setText(obj_name)
            self.lbl_pdb_info.setText(f"PDB: {p.name} (Object: '{obj_name}')")
            QtWidgets.QMessageBox.information(
                self, "PDB Loaded", f"Successfully loaded PDB structure '{p.name}' into PyMOL as '{obj_name}'."
            )
        except Exception as err:
            logger.exception("Failed loading PDB file into PyMOL")
            QtWidgets.QMessageBox.critical(self, "Error", f"Failed loading PDB file:\n{err}")

    @Slot()
    def _export_res_csv(self) -> None:
        if self._res_means is None or self._res_stds is None:
            QtWidgets.QMessageBox.warning(self, "No Data", "Please load a per-residue CSV file first.")
            return

        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save Residue Statistics CSV", "residue_statistics.csv", "CSV Files (*.csv)"
        )
        if path:
            start_res = self.get_start_residue()
            MMISMSAExporter.export_residue_statistics_csv(self._res_means, self._res_stds, Path(path), start_res=start_res)
            QtWidgets.QMessageBox.information(self, "Success", f"Exported residue statistics to:\n{path}")

    @Slot()
    def _map_pymol_residues(self) -> None:
        if self._res_means is None or self._res_stds is None:
            QtWidgets.QMessageBox.warning(self, "No Data", "Please load a per-residue CSV file first.")
            return

        try:
            term = self.combo_term.currentText()
            metric = self.combo_metric.currentText()
            term_key = resolve_residue_term_key(term, self._res_means)
            vals = self._res_means[term_key] if metric == "mean" else self._res_stds[term_key]
            sel = self.txt_sel.text().strip() or "all"
            pal = self.combo_pal.currentText()
            start_res = self.get_start_residue()

            # Clean previous scalebar ramp object if present
            if "mmismsa_scalebar" in cmd.get_names():
                cmd.delete("mmismsa_scalebar")

            # Map values to atom b-factors using custom residue numbering offset
            for idx, val in enumerate(vals):
                res_num = start_res + idx
                target = f"({sel}) and resi {res_num}"
                cmd.alter(target, f"b={val}")

            cmd.sort(sel)
            min_val = float(np.min(vals))
            max_val = float(np.max(vals))

            # Avoid division by zero if all values are identical
            if min_val == max_val:
                max_val = min_val + 0.001

            # Map palette names to explicit PyMOL spectrum color strings
            SPECTRUM_COLOR_MAP = {
                "bwr": "blue white red",
                "rainbow": "blue green yellow red",
                "red_white_blue": "red white blue",
                "coolwarm": "blue white red"
            }
            spec_colors = SPECTRUM_COLOR_MAP.get(pal, "blue white red")

            # Apply PyMOL spectrum with explicit color bounds
            cmd.spectrum("b", spec_colors, sel, minimum=min_val, maximum=max_val)

            # Retrieve robust ramp color definition for scalebar
            ramp_colors = PALETTE_RAMP_MAP.get(pal, ["blue", "white", "red"])
            cmd.ramp_new("mmismsa_scalebar", sel, range=[min_val, max_val], color=ramp_colors)
            cmd.rebuild()

            QtWidgets.QMessageBox.information(
                self, "Success", 
                f"Projected per-residue {term} ({metric}) values across {len(vals)} residues in PyMOL!\n"
                f"Residue range: {start_res} to {start_res + len(vals) - 1}\n"
                f"Scalebar created: [{min_val:.2f} to {max_val:.2f}] kcal/mol"
            )
        except Exception as err:
            logger.error("Residue mapping failed: %s", err)
            QtWidgets.QMessageBox.critical(self, "Error", f"Residue mapping failed:\n{err}")


class ResiduePlotsTab(QtWidgets.QWidget):
    """Tab dedicated cleanly to plotting per-residue energy profiles."""

    def __init__(self, res_stats_tab: ResidueAnalysisTab, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._res_stats_tab = res_stats_tab
        self._has_plotted_profile = False
        self._setup_ui()
        self._res_stats_tab.spin_start_res.valueChanged.connect(self._on_start_residue_changed)

    def _setup_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        ctrls = QtWidgets.QHBoxLayout()
        ctrls.addWidget(QtWidgets.QLabel("Select Term to Plot:"))
        self.combo_term = QtWidgets.QComboBox()
        self.combo_term.addItems(RESIDUE_TERMS2)
        ctrls.addWidget(self.combo_term)

        self.btn_plot = QtWidgets.QPushButton("Plot Residue Profile")
        self.btn_plot.setStyleSheet("font-weight: 600;")
        self.btn_plot.clicked.connect(self.plot_profile)
        ctrls.addWidget(self.btn_plot)

        self.btn_save_plot = QtWidgets.QPushButton("Save Plot")
        self.btn_save_plot.clicked.connect(self._save_plot)
        ctrls.addWidget(self.btn_save_plot)

        ctrls.addStretch()

        layout.addLayout(ctrls)

        self.canvas_widget = MatplotlibWidget()
        layout.addWidget(self.canvas_widget)

    def _plot_profile_impl(self, show_warnings: bool = True) -> bool:
        means_dict = self._res_stats_tab._res_means
        stds_dict = self._res_stats_tab._res_stds

        if means_dict is None or stds_dict is None:
            if show_warnings:
                QtWidgets.QMessageBox.warning(
                    self, "No Residue Data", "Please load a per-residue CSV file in the 'Residue Analysis' tab first."
                )
            return False
        if not ensure_matplotlib_available(self, "plot residue profile"):
            return False

        term = self.combo_term.currentText()
        term_key = resolve_residue_term_key(term, means_dict)
        means = means_dict[term_key]
        stds = stds_dict[term_key]
        start_res = self._res_stats_tab.get_start_residue()
        self.canvas_widget.plot_residue_profile(means, stds, term_name=term_key, start_res=start_res)
        return True

    def _on_start_residue_changed(self, _value: int) -> None:
        # Refresh an existing plot so x-axis residue numbering follows the new start residue.
        if self._has_plotted_profile:
            self._plot_profile_impl(show_warnings=False)

    @Slot()
    def _save_plot(self) -> None:
        if not ensure_matplotlib_available(self, "save residue plot"):
            return
        if not self._has_plotted_profile:
            QtWidgets.QMessageBox.warning(
                self,
                "No Plot",
                "Please generate a residue profile plot before saving.",
            )
            return
        if self.canvas_widget.figure is None:
            QtWidgets.QMessageBox.critical(
                self,
                "Save Error",
                "Plot backend is not available.",
            )
            return

        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save Residue Plot",
            "residue_profile.png",
            "PNG Files (*.png);;PDF Files (*.pdf);;SVG Files (*.svg)",
        )
        if not path:
            return

        out_path = Path(path)
        if out_path.suffix == "":
            out_path = out_path.with_suffix(".png")

        try:
            self.canvas_widget.figure.savefig(
                str(out_path),
                dpi=300,
                facecolor=self.canvas_widget.figure.get_facecolor(),
                edgecolor="none",
            )
            QtWidgets.QMessageBox.information(
                self,
                "Saved",
                f"Residue plot saved to:\n{out_path}",
            )
        except Exception as err:
            QtWidgets.QMessageBox.critical(
                self,
                "Save Error",
                f"Failed to save residue plot:\n{err}",
            )

    @Slot()
    def plot_profile(self) -> None:
        self._has_plotted_profile = self._plot_profile_impl(show_warnings=True)


class cMMISMSAExecutionTab(QtWidgets.QWidget):
    """Tab to launch the cMMISMSA binary from the plugin and capture its output."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._process: Optional[subprocess.Popen] = None
        self._reader_thread: Optional[threading.Thread] = None
        self._output_queue: Optional[queue.Queue] = None
        self._setup_ui()

    def _setup_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        info = QtWidgets.QLabel(
            "Use this tab to run the cMMISMSA executable directly from PyMOL. Select the binary, add CLI arguments, and review the console output."
        )
        info.setWordWrap(True)
        layout.addWidget(info)

        row_exec = QtWidgets.QHBoxLayout()
        row_exec.addWidget(QtWidgets.QLabel("Executable:"))
        self.txt_binary = QtWidgets.QLineEdit()
        self.txt_binary.setPlaceholderText("Path to cMMISMSA binary")
        self.txt_binary.setText(self._find_default_binary())
        row_exec.addWidget(self.txt_binary, stretch=1)
        self.btn_browse = QtWidgets.QPushButton("Browse")
        self.btn_browse.clicked.connect(self._browse_binary)
        row_exec.addWidget(self.btn_browse)
        layout.addLayout(row_exec)

        row_args = QtWidgets.QHBoxLayout()
        row_args.addWidget(QtWidgets.QLabel("Arguments:"))
        self.txt_args = QtWidgets.QLineEdit()
        self.txt_args.setPlaceholderText("Example: --topology step3_input.parm7 -xtc trajectory.xtc --mask 123 --output energy")
        row_args.addWidget(self.txt_args, stretch=1)
        layout.addLayout(row_args)

        row_work = QtWidgets.QHBoxLayout()
        row_work.addWidget(QtWidgets.QLabel("Working directory:"))
        self.txt_workdir = QtWidgets.QLineEdit()
        self.txt_workdir.setText(str(Path.cwd()))
        row_work.addWidget(self.txt_workdir, stretch=1)
        layout.addLayout(row_work)

        row_buttons = QtWidgets.QHBoxLayout()
        self.btn_run = QtWidgets.QPushButton("Run cMMISMSA")
        self.btn_run.setStyleSheet("font-weight: 600;")
        self.btn_run.clicked.connect(self._run_binary)
        row_buttons.addWidget(self.btn_run)

        self.btn_stop = QtWidgets.QPushButton("Stop")
        self.btn_stop.clicked.connect(self._stop_binary)
        self.btn_stop.setEnabled(False)
        row_buttons.addWidget(self.btn_stop)
        layout.addLayout(row_buttons)

        self.lbl_status = QtWidgets.QLabel("Status: idle")
        layout.addWidget(self.lbl_status)

        self.log_output = QtWidgets.QPlainTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setPlaceholderText("Process output will appear here...")
        layout.addWidget(self.log_output, stretch=1)

        self._output_queue = queue.Queue()
        self._timer = QtCore.QTimer(self)
        self._timer.timeout.connect(self._drain_output)
        self._timer.start(100)

    def _find_default_binary(self) -> str:
        candidates = ["cMMISMSA", "cMMISMSA.exe", "cmmismsa", "cmmismsa.exe"]
        for name in candidates:
            found = shutil.which(name)
            if found:
                return found
        return ""

    @Slot()
    def _browse_binary(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Select cMMISMSA binary", "", "Executables (*)"
        )
        if path:
            self.txt_binary.setText(path)

    @Slot()
    def _run_binary(self) -> None:
        binary = self.txt_binary.text().strip()
        if not binary:
            QtWidgets.QMessageBox.warning(self, "Missing binary", "Please select the cMMISMSA executable first.")
            return

        binary_path = Path(binary)
        if not binary_path.exists() and not shutil.which(binary):
            QtWidgets.QMessageBox.warning(self, "Executable not found", f"The selected path does not exist: {binary}")
            self.lbl_status.setText("Status: error (executable not found)")
            return

        args_text = self.txt_args.text().strip()
        args = shlex.split(args_text) if args_text else []
        workdir = self.txt_workdir.text().strip() or str(Path.cwd())

        self.log_output.clear()
        self.log_output.appendPlainText(f"> {binary} {' '.join(args) if args else ''}".strip())
        self.lbl_status.setText("Status: launching...")
        self.btn_run.setEnabled(False)
        self.btn_stop.setEnabled(True)

        self._output_queue = queue.Queue()
        program, proc_args = self._build_command(binary, args)
        try:
            self._process = subprocess.Popen(
                [program] + proc_args,
                cwd=workdir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL,
                text=True,
                bufsize=1,
                universal_newlines=True,
            )
        except Exception as err:
            self._process = None
            self.lbl_status.setText(f"Status: error ({err})")
            QtWidgets.QMessageBox.critical(self, "Execution Error", str(err))
            self.btn_run.setEnabled(True)
            self.btn_stop.setEnabled(False)
            return

        self._reader_thread = threading.Thread(target=self._pump_output, daemon=True)
        self._reader_thread.start()

    @Slot()
    def _stop_binary(self) -> None:
        if self._process is not None:
            try:
                self._process.terminate()
                self._process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                self._process.kill()
                self._process.wait(timeout=2)
            except Exception:
                pass
            self._finish_process()

    def _build_command(self, binary: str, args: List[str]) -> Tuple[str, List[str]]:
        if sys.platform.startswith("win") and Path(binary).suffix.lower() in {".bat", ".cmd"}:
            return "cmd.exe", ["/c", binary] + args
        return binary, args

    def _pump_output(self) -> None:
        if not self._process or not self._process.stdout:
            self._output_queue.put("Unable to capture process output.\n")
            self._output_queue.put(None)
            return

        try:
            for line in self._process.stdout:
                if line:
                    self._output_queue.put(line)
        finally:
            if self._process.stdout:
                self._process.stdout.close()
            self._output_queue.put(None)

    def _drain_output(self) -> None:
        if not self._output_queue:
            return

        while True:
            try:
                item = self._output_queue.get_nowait()
            except queue.Empty:
                break

            if item is None:
                self._finish_process()
                break

            self.log_output.insertPlainText(item)
            self.log_output.ensureCursorVisible()

    def _finish_process(self) -> None:
        if self._process is None:
            return

        exit_code = self._process.poll()
        if exit_code is None:
            return

        self.btn_run.setEnabled(True)
        self.btn_stop.setEnabled(False)
        self.lbl_status.setText(f"Status: finished (exit code {exit_code})")
        self._process = None
        self._reader_thread = None

# ==============================================================================
# 6. MAIN WINDOW & PYMOL ENTRYPOINT
# ==============================================================================
class MMISMSAAnalyzerDock(QtWidgets.QDockWidget):
    """Main Dockable Panel inside PyMOL."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__("MMISMSA Analyzer", parent)
        self.setObjectName("MMISMSAAnalyzerDock")
        self._dataset: Optional[MMISMSADataSet] = None
        self._setup_ui()

    def _setup_ui(self) -> None:
        self._apply_modern_theme()

        container = QtWidgets.QWidget()
        main = QtWidgets.QVBoxLayout(container)

        # Fused Tabs
        self.tabs = QtWidgets.QTabWidget()
        self.tab_run_binary = cMMISMSAExecutionTab()
        self.tab_combined_global = CombinedGlobalAnalysisTab()
        self.tab_combined_global.set_load_handler(self._load_csv)
        self.tab_traj = TrajectoryTab()
        self.tab_res_analysis = ResidueAnalysisTab()
        self.tab_res_plots = ResiduePlotsTab(res_stats_tab=self.tab_res_analysis)

        # Keep execution workflow as the first tab for quicker access.
        self.tabs.addTab(self.tab_run_binary, "Run cMMISMSA")
        self.tabs.addTab(self.tab_combined_global, "Global Analysis")
        self.tabs.addTab(self.tab_traj, "Trajectory")
        self.tabs.addTab(self.tab_res_analysis, "Residue Analysis")
        self.tabs.addTab(self.tab_res_plots, "Residue Plots")
        
        # Developer & Institutional Acknowledgements Tab
        self.init_acknowledgements_tab()

        main.addWidget(self.tabs)
        self.setWidget(container)
        self.resize(880, 680)

    def _apply_modern_theme(self) -> None:
        """Applies a neutral, cohesive Qt stylesheet to the full dock UI."""
        self.setStyleSheet(f"""
            QDockWidget {{
                background: {UI_COLORS['bg']};
                color: {UI_COLORS['text']};
                border: 1px solid {UI_COLORS['border']};
            }}
            QDockWidget::title {{
                background: {UI_COLORS['panel_alt']};
                color: {UI_COLORS['text']};
                padding: 6px 10px;
                border-bottom: 1px solid {UI_COLORS['border']};
                text-align: left;
            }}
            QWidget {{
                background: {UI_COLORS['bg']};
                color: {UI_COLORS['text']};
                font-size: 10pt;
            }}
            QGroupBox {{
                background: {UI_COLORS['panel']};
                border: 1px solid {UI_COLORS['border']};
                border-radius: 8px;
                margin-top: 10px;
                padding: 12px;
                font-weight: 600;
            }}
            QGroupBox::title {{
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 6px;
                color: {UI_COLORS['muted']};
            }}
            QPushButton {{
                background: {UI_COLORS['accent']};
                color: #ffffff;
                border: 1px solid {UI_COLORS['accent_pressed']};
                border-radius: 6px;
                padding: 6px 10px;
                font-weight: 600;
            }}
            QPushButton:hover {{
                background: {UI_COLORS['accent_hover']};
            }}
            QPushButton:pressed {{
                background: {UI_COLORS['accent_pressed']};
            }}
            QPushButton:disabled {{
                background: #aeb7bf;
                border-color: #9ca7b1;
                color: #edf1f4;
            }}
            QLineEdit, QComboBox, QSpinBox {{
                background: {UI_COLORS['panel']};
                border: 1px solid {UI_COLORS['border']};
                border-radius: 6px;
                padding: 5px 7px;
                selection-background-color: #9eacb8;
            }}
            QTableWidget {{
                background: {UI_COLORS['panel']};
                alternate-background-color: #f6f8fa;
                border: 1px solid {UI_COLORS['border']};
                gridline-color: #dde3e8;
            }}
            QHeaderView::section {{
                background: {UI_COLORS['panel_alt']};
                color: {UI_COLORS['text']};
                border: 1px solid {UI_COLORS['border']};
                padding: 6px;
                font-weight: 600;
            }}
            QTabWidget::pane {{
                border: 1px solid {UI_COLORS['border']};
                background: {UI_COLORS['panel']};
                border-radius: 8px;
                top: -1px;
            }}
            QTabBar::tab {{
                background: {UI_COLORS['panel_alt']};
                color: {UI_COLORS['muted']};
                border: 1px solid {UI_COLORS['border']};
                border-bottom: none;
                border-top-left-radius: 7px;
                border-top-right-radius: 7px;
                padding: 8px 14px;
                margin-right: 2px;
            }}
            QTabBar::tab:selected {{
                background: {UI_COLORS['panel']};
                color: {UI_COLORS['text']};
            }}
            QLabel {{
                color: {UI_COLORS['text']};
            }}
            QPlainTextEdit {{
                background: {UI_COLORS['panel']};
                border: 1px solid {UI_COLORS['border']};
                border-radius: 8px;
                color: {UI_COLORS['text']};
            }}
        """)

    def init_acknowledgements_tab(self) -> None:
        """Renders developer metadata, laboratory affiliations, institutional logos, and software license states."""
        tab = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout()

        title = QtWidgets.QLabel(f"<h1 align='center'>MMISMSA Analyzer v{VERSION}</h1>")
        title.setAlignment(QtCore.Qt.AlignCenter)

        info = QtWidgets.QLabel(
            "<p align='center' style='font-size:11pt; line-height:140%;'>"
            "<b>Developer:</b> Javier García Marín<br>"
            "<i>Department of Organic and Inorganic Chemistry</i><br>"
            "<b>University of Alcalá</b><br>"
            "<a href='mailto:javier.garciamarin@uah.es'>javier.garciamarin@uah.es</a><br><br>"
            "<font color='#880088'><b>License: GPLv3</b></font><br>"
            "</p>"
        )
        info.setAlignment(QtCore.Qt.AlignCenter)
        info.setOpenExternalLinks(True)

        logo_label = QtWidgets.QLabel()
        logo_label.setText(
            f"<p align='center'><img src='data:image/png;base64,{PERSONAL_LOGO_BASE64.decode('utf-8')}' "
            "width='100' height='100'></p>"
        )
        logo_label.setAlignment(QtCore.Qt.AlignCenter)

        layout.addStretch()
        layout.addWidget(title)
        layout.addWidget(info)
        layout.addWidget(logo_label)
        layout.addStretch()

        tab.setLayout(layout)
        self.tabs.addTab(tab, "About")

    @Slot()
    def _load_csv(self) -> None:
        f_str, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open MMISMSA global.csv File", "", "CSV Files (*.csv *.txt);;All Files (*)"
        )
        if not f_str: return

        p = Path(f_str)
        try:
            ds = MMISMSACSVReader.load_file(p)
            self._dataset = ds
            self.tab_combined_global.set_dataset(ds, p)
            self.tab_traj.set_dataset(ds)

            self.tab_combined_global.set_file_info(f"Loaded: {p.name} ({ds.frame_count} frames)")
            QtWidgets.QMessageBox.information(self, "Loaded", f"Parsed {ds.frame_count} trajectory frames.")
        except Exception as err:
            logger.exception("Error loading CSV")
            QtWidgets.QMessageBox.critical(self, "Error", f"Failed loading CSV:\n{err}")


# Global reference preventing GC
_plugin_window: Optional[MMISMSAAnalyzerDock] = None

def __init_plugin__(app=None) -> None:
    """PyMOL plugin manager entry point."""
    def launch_plugin():
        global _plugin_window
        parent = plugins.get_qt_parent() if hasattr(plugins, "get_qt_parent") else None
        if _plugin_window is None:
            _plugin_window = MMISMSAAnalyzerDock(parent=parent)
        _plugin_window.show()
        _plugin_window.raise_()
        _plugin_window.activateWindow()

    # Prefer Qt menu registration so the plugin appears in the main Plugins menu.
    if hasattr(plugins, "addmenuitemqt"):
        plugins.addmenuitemqt("MMISMSA Analyzer", launch_plugin)
        logger.info("MMISMSA Analyzer registered in main Plugins menu (Qt).")
    else:
        plugins.addmenuitem("MMISMSA Analyzer", launch_plugin)
        logger.info("MMISMSA Analyzer registered with legacy plugin menu API.")
