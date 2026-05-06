#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RIS-Assisted CubeSat-to-Ground Communication Simulator
Ports MATLAB simulation suite (main.m, pso_convergence.m, elevation_sweep.m)
to an interactive Python/PyQt6 application.

Author: APache | MCTE | 2026
"""

import sys, subprocess

def _ensure(pkg):
    try:
        __import__(pkg.split("[")[0].replace("-", "_"))
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg, "-q"])

for _p in ["PyQt6", "numpy", "matplotlib", "scipy"]:
    _ensure(_p)

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import matplotlib
matplotlib.use("QtAgg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QTabWidget, QPushButton, QLabel, QSpinBox, QProgressBar,
    QFrame, QSizePolicy, QFileDialog,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt6.QtGui import QFont

# ─────────────────────────────────────────────────────────────────────────────
#  Design tokens
# ─────────────────────────────────────────────────────────────────────────────
_BG    = "#07101f"
_CARD  = "#0d1e35"
_PANEL = "#0a1828"
_LINE  = "#162a45"
_CYAN  = "#00c8ff"
_TEXT  = "#c4d4e8"
_DIM   = "#3d5470"
_PLT   = "#0d1e35"
_GRID  = "#111f30"

_RED   = "#f87171"
_AMB   = "#fbbf24"
_BLUE  = "#60a5fa"
_GRN   = "#4ade80"

CSS = f"""
QMainWindow, QWidget {{
    background:{_BG}; color:{_TEXT};
    font-family:'Segoe UI',Arial,sans-serif; font-size:13px;
}}
QTabWidget::pane {{ border:1px solid {_LINE}; background:{_CARD}; }}
QTabBar::tab {{
    background:{_PANEL}; color:{_DIM};
    padding:10px 22px; border:none; border-bottom:3px solid transparent;
    font-size:13px; font-weight:600;
}}
QTabBar::tab:selected {{ color:{_CYAN}; border-bottom:3px solid {_CYAN}; background:{_CARD}; }}
QTabBar::tab:hover:!selected {{ color:{_TEXT}; background:{_CARD}; }}
QPushButton {{
    background:{_CYAN}; color:#00060f; border:none;
    padding:9px 0; border-radius:6px; font-weight:700; font-size:13px;
}}
QPushButton:hover {{ background:#22d3ee; }}
QPushButton:pressed {{ background:#0891b2; }}
QPushButton:disabled {{ background:#152030; color:#2a4060; }}
QPushButton[flat="true"] {{
    background:transparent; color:{_DIM}; border:1px solid {_LINE};
}}
QPushButton[flat="true"]:hover {{ color:{_CYAN}; border-color:{_CYAN}; }}
QLabel {{ color:{_TEXT}; background:transparent; }}
QSpinBox {{
    background:{_PANEL}; color:{_CYAN}; border:1px solid {_LINE};
    border-radius:5px; padding:4px 8px; font-size:13px; font-weight:700;
    min-width:72px; max-width:90px;
}}
QSpinBox:focus {{ border-color:{_CYAN}; }}
QSpinBox::up-button, QSpinBox::down-button {{
    background:{_LINE}; border:none; width:18px; border-radius:3px;
}}
QSpinBox::up-button:hover, QSpinBox::down-button:hover {{
    background:#1e3a5f;
}}
QSpinBox::up-arrow  {{ image:none; width:0; height:0;
    border-left:4px solid transparent; border-right:4px solid transparent;
    border-bottom:5px solid {_CYAN}; margin:auto; }}
QSpinBox::down-arrow {{ image:none; width:0; height:0;
    border-left:4px solid transparent; border-right:4px solid transparent;
    border-top:5px solid {_CYAN}; margin:auto; }}
QProgressBar {{
    border:none; background:{_PANEL}; border-radius:3px; height:6px;
    text-align:center; color:transparent;
}}
QProgressBar::chunk {{
    background:qlineargradient(x1:0,y1:0,x2:1,y2:0,
        stop:0 #0891b2,stop:1 {_CYAN});
    border-radius:3px;
}}
"""

# ─────────────────────────────────────────────────────────────────────────────
#  Simulation engine  (faithful vectorised port of MATLAB files)
# ─────────────────────────────────────────────────────────────────────────────

def rician_channel(n: int, k_db: float) -> np.ndarray:
    k    = 10.0 ** (k_db / 10.0)
    los  = np.sqrt(k / (k + 1)) * np.exp(1j * 2 * np.pi * np.random.rand(n))
    nlos = np.sqrt(1.0 / (2*(k+1))) * (np.random.randn(n) + 1j*np.random.randn(n))
    return los + nlos


def k_db_elev(elev_deg: float) -> float:
    return -5.0 + 17.0 * np.sin(np.deg2rad(elev_deg))


def pso_optimize(h_d, cascaded, swarm=20, n_iter=30):
    n = len(cascaded)
    w, c1, c2, vmax = 0.7, 1.5, 1.5, np.pi / 4

    pos = 2*np.pi * np.random.rand(swarm, n)
    vel = np.zeros((swarm, n))
    fit = np.abs(h_d + np.exp(1j*pos) @ cascaded) ** 2

    pb_pos = pos.copy();  pb_fit = fit.copy()
    gi = int(np.argmax(pb_fit))
    gb_fit = float(pb_fit[gi]);  gb_pos = pb_pos[gi].copy()
    hist = np.zeros(n_iter)

    for it in range(n_iter):
        r1, r2 = np.random.rand(swarm, n), np.random.rand(swarm, n)
        vel = w*vel + c1*r1*(pb_pos-pos) + c2*r2*(gb_pos[np.newaxis]-pos)
        vel = np.clip(vel, -vmax, vmax)
        pos = np.mod(pos + vel, 2*np.pi)
        fit = np.abs(h_d + np.exp(1j*pos) @ cascaded) ** 2
        up  = fit > pb_fit;  pb_pos[up] = pos[up];  pb_fit[up] = fit[up]
        ni  = int(np.argmax(pb_fit))
        if pb_fit[ni] > gb_fit:
            gb_fit = float(pb_fit[ni]);  gb_pos = pb_pos[ni].copy()
        hist[it] = gb_fit

    return np.exp(1j*gb_pos), hist


def _bpsk_detect(h_eff, x, noise):
    y = h_eff * x + noise
    return (np.real(y * np.conj(h_eff)) > 0).astype(int)


def sim_ber(snr_arr, elev, N, n_mc, n_bits, sw, ni, cb=None):
    k_d, k_r = k_db_elev(elev), 10.0
    m = len(snr_arr)
    out = np.zeros((4, m))
    for i, snr in enumerate(snr_arr):
        nv = 1.0 / (10**(snr/10))
        errs = np.zeros(4, dtype=int);  total = 0
        for _ in range(n_mc):
            h_d  = rician_channel(1, k_d)[0]
            casc = rician_channel(N, k_r) * rician_channel(N, k_r)
            phi_r = np.exp(1j * 2*np.pi * np.random.rand(N))
            phi_o = np.exp(1j * (np.angle(h_d) - np.angle(casc)))
            phi_p, _ = pso_optimize(h_d, casc, sw, ni)
            Hs = [h_d,
                  h_d + casc@phi_r,
                  h_d + casc@phi_p,
                  h_d + casc@phi_o]
            bits  = np.random.randint(0, 2, n_bits)
            x     = 2*bits - 1
            noise = np.sqrt(nv/2)*(np.random.randn(n_bits)+1j*np.random.randn(n_bits))
            for j, h in enumerate(Hs):
                errs[j] += np.sum(_bpsk_detect(h, x, noise) != bits)
            total += n_bits
        out[:, i] = np.maximum(errs/total, 1e-6)
        if cb: cb(int((i+1)/m*100))
    return out


def sim_pso_conv(N, k_d, k_r, sw, ni, trials):
    h_d  = rician_channel(1, k_d)[0]
    casc = rician_channel(N, k_r) * rician_channel(N, k_r)
    phi_o   = np.exp(1j*(np.angle(h_d) - np.angle(casc)))
    opt_fit = float(np.abs(h_d + casc@phi_o)**2)
    hists   = [pso_optimize(h_d, casc, sw, ni)[1] for _ in range(trials)]
    return np.array(hists), opt_fit


def sim_elev(elev_arr, snr, N, n_mc, n_bits, sw, ni, cb=None):
    k_r = 10.0;  nv = 1.0/(10**(snr/10))
    out = np.zeros((2, len(elev_arr)))
    for i, elev in enumerate(elev_arr):
        k_d = k_db_elev(elev)
        errs = np.zeros(2, dtype=int);  total = 0
        for _ in range(n_mc):
            h_d  = rician_channel(1, k_d)[0]
            casc = rician_channel(N, k_r) * rician_channel(N, k_r)
            phi_p, _ = pso_optimize(h_d, casc, sw, ni)
            Hs = [h_d, h_d + casc@phi_p]
            bits  = np.random.randint(0, 2, n_bits)
            x     = 2*bits - 1
            noise = np.sqrt(nv/2)*(np.random.randn(n_bits)+1j*np.random.randn(n_bits))
            for j, h in enumerate(Hs):
                errs[j] += np.sum(_bpsk_detect(h, x, noise) != bits)
            total += n_bits
        out[:, i] = np.maximum(errs/total, 1e-6)
        if cb: cb(int((i+1)/len(elev_arr)*100))
    return out


# ─────────────────────────────────────────────────────────────────────────────
#  Worker threads
# ─────────────────────────────────────────────────────────────────────────────

class _Worker(QThread):
    progress = pyqtSignal(int)
    done     = pyqtSignal(object)
    def __init__(self, fn, kw): super().__init__(); self._fn = fn; self._kw = kw
    def run(self): self.done.emit(self._fn(**self._kw))


# ─────────────────────────────────────────────────────────────────────────────
#  UI helpers
# ─────────────────────────────────────────────────────────────────────────────

def _sep():
    f = QFrame()
    f.setFrameShape(QFrame.Shape.HLine)
    f.setStyleSheet(f"background:{_LINE}; max-height:1px; border:none;")
    return f


def _label(text, size=13, color=_TEXT, bold=False, dim=False):
    w = QLabel(text)
    c = _DIM if dim else color
    w.setStyleSheet(
        f"color:{c}; font-size:{size}px;"
        f"font-weight:{'700' if bold else '400'}; background:transparent;")
    return w


class _SpinField(QWidget):
    """Label + QSpinBox row — drop-in replacement for the old slider widget."""
    valueChanged = pyqtSignal(int)

    def __init__(self, label, lo, hi, val, unit="", parent=None):
        super().__init__(parent)
        self.setFixedHeight(34)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(8)

        lbl = _label(label, 12, _DIM)
        row.addWidget(lbl)
        row.addStretch()

        self._box = QSpinBox()
        self._box.setRange(lo, hi)
        self._box.setValue(val)
        self._box.setSuffix(unit)
        self._box.setAlignment(Qt.AlignmentFlag.AlignRight)
        self._box.setButtonSymbols(QSpinBox.ButtonSymbols.UpDownArrows)
        self._box.valueChanged.connect(self.valueChanged)
        row.addWidget(self._box)

    def value(self):       return self._box.value()
    def setValue(self, v): self._box.setValue(v)


class _Card(QWidget):
    def __init__(self, label, value="--", color=_CYAN, parent=None):
        super().__init__(parent)
        self.setFixedHeight(62)
        self.setStyleSheet(
            f"background:{_PANEL}; border-radius:6px;"
            f"border-left:3px solid {color}; border:1px solid {_LINE};"
            f"border-left-color:{color};")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 6, 12, 6)
        lay.setSpacing(1)
        self._v = QLabel(value)
        self._v.setStyleSheet(
            f"color:{color}; font-size:17px; font-weight:700;"
            f"background:transparent; border:none;")
        self._l = QLabel(label)
        self._l.setStyleSheet(
            f"color:{_DIM}; font-size:10px; background:transparent; border:none;")
        self._l.setWordWrap(True)
        lay.addWidget(self._v)
        lay.addWidget(self._l)

    def set(self, v): self._v.setText(str(v))


class _Canvas(FigureCanvas):
    def __init__(self, fig, parent=None):
        super().__init__(fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.updateGeometry()


def _ax_style(ax):
    ax.set_facecolor(_PLT)
    ax.tick_params(colors=_TEXT, labelsize=10, length=4, width=0.6)
    ax.xaxis.label.set_color(_TEXT)
    ax.yaxis.label.set_color(_TEXT)
    ax.title.set_color("#dde8f8")
    ax.title.set_fontsize(13)
    ax.title.set_fontweight("bold")
    for sp in ax.spines.values():
        sp.set_color(_GRID)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, color=_GRID, linewidth=0.8, linestyle="--", alpha=1.0)
    ax.set_axisbelow(True)
    return ax


def _side_panel(parent, width=268):
    w = QWidget(parent)
    w.setFixedWidth(width)
    w.setStyleSheet(f"background:{_PANEL}; border-right:1px solid {_LINE};")
    lay = QVBoxLayout(w)
    lay.setContentsMargins(18, 18, 18, 18)
    lay.setSpacing(8)
    return w, lay


# ─────────────────────────────────────────────────────────────────────────────
#  Tab 1 — BER vs SNR
# ─────────────────────────────────────────────────────────────────────────────

class BERTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._worker = None
        self._setup()

    def _setup(self):
        root = QHBoxLayout(self)
        root.setSpacing(0)
        root.setContentsMargins(0, 0, 0, 0)

        # ── control panel ─────────────────────────────────────────────────
        side, sl = _side_panel(self)

        sl.addWidget(_label("Parameters", 11, _DIM))
        sl.addWidget(_sep())

        self.s_N    = _SpinField("RIS Elements  N",  4, 128,  64)
        self.s_elev = _SpinField("Elevation angle",  5,  85,  10, "°")
        self.s_mc   = _SpinField("Monte Carlo trials", 10, 200, 40)
        self.s_bits = _SpinField("Bits per trial",   50, 500, 200)
        self.s_sw   = _SpinField("PSO swarm size",    5,  50,  20)
        self.s_it   = _SpinField("PSO iterations",    5,  80,  30)
        for s in [self.s_N, self.s_elev, self.s_mc, self.s_bits, self.s_sw, self.s_it]:
            sl.addWidget(s)

        self._k_lbl = _label("", 10, _DIM)
        self._refresh_k()
        self.s_elev.valueChanged.connect(lambda _: self._refresh_k())
        sl.addWidget(self._k_lbl)

        sl.addWidget(_sep())

        self._btn = QPushButton("▶   Run Simulation")
        self._btn.clicked.connect(self._run)
        self._exp = QPushButton("Export PNG")
        self._exp.setProperty("flat", "true")
        self._exp.clicked.connect(self._export)
        sl.addWidget(self._btn)
        sl.addWidget(self._exp)

        sl.addWidget(_sep())
        self._c_gain = _Card("PSO gain vs No-RIS  (BER 10⁻³)", "--", _BLUE)
        self._c_gap  = _Card("Gap to optimal", "--", _GRN)
        sl.addWidget(self._c_gain)
        sl.addWidget(self._c_gap)

        sl.addStretch()
        self._pbar = QProgressBar()
        self._pbar.setValue(0)
        self._stat = _label("Ready", 11, _DIM)
        sl.addWidget(self._pbar)
        sl.addWidget(self._stat)

        root.addWidget(side)

        # ── plot ──────────────────────────────────────────────────────────
        pw = QWidget()
        pw.setStyleSheet(f"background:{_CARD};")
        pl = QVBoxLayout(pw)
        pl.setContentsMargins(20, 16, 20, 16)

        self._fig = Figure(figsize=(9, 6), facecolor=_PLT)
        self._ax  = self._fig.add_subplot(111)
        _ax_style(self._ax)
        self._cv  = _Canvas(self._fig, pw)
        pl.addWidget(self._cv)
        self._idle("BER vs SNR  |  CubeSat Link Performance")

        root.addWidget(pw)

    def _refresh_k(self):
        e = self.s_elev.value()
        k = 10 ** (k_db_elev(e) / 10)
        self._k_lbl.setText(f"  Direct K-factor: {k_db_elev(e):.1f} dB  (K ≈ {k:.2f})")

    def _idle(self, title):
        self._ax.clear(); _ax_style(self._ax)
        self._ax.text(0.5, 0.5, "Press  ▶  Run Simulation  to begin",
                      transform=self._ax.transAxes, ha="center", va="center",
                      color=_DIM, fontsize=15, style="italic")
        self._ax.set_title(title)
        self._cv.draw()

    def _run(self):
        if self._worker and self._worker.isRunning(): return
        self._btn.setEnabled(False)
        self._pbar.setValue(0)
        self._stat.setText("Running...")
        self._c_gain.set("--"); self._c_gap.set("--")

        self._p = dict(
            snr_arr=np.arange(-10, 22, 2, dtype=float),
            elev=self.s_elev.value(), N=self.s_N.value(),
            n_mc=self.s_mc.value(), n_bits=self.s_bits.value(),
            sw=self.s_sw.value(), ni=self.s_it.value())

        def _fn(snr_arr, elev, N, n_mc, n_bits, sw, ni):
            return sim_ber(snr_arr, elev, N, n_mc, n_bits, sw, ni,
                           cb=self._worker.progress.emit)

        self._worker = _Worker(_fn, self._p)
        self._worker.progress.connect(self._pbar.setValue)
        self._worker.done.connect(self._done)
        self._worker.start()

    def _done(self, out):
        snr = self._p["snr_arr"]
        N = self._p["N"]; elv = self._p["elev"]
        ber_no, ber_rand, ber_pso, ber_opt = out

        self._ax.clear(); _ax_style(self._ax)
        kw = dict(linewidth=2.2, markevery=2, markersize=7)
        self._ax.semilogy(snr, ber_no,   color=_RED,  marker="o",  label="No RIS", **kw)
        self._ax.semilogy(snr, ber_rand, color=_AMB,  marker="s",  linestyle="--",
                          label="Random RIS", **kw)
        self._ax.semilogy(snr, ber_pso,  color=_BLUE, marker="^",  label="PSO-Optimised RIS", **kw)
        self._ax.semilogy(snr, ber_opt,  color=_GRN,  marker="d",  label="Closed-form Optimal", **kw)

        tgt = 1e-3
        self._ax.axhline(tgt, color=_DIM, linestyle=":", linewidth=1.1)
        self._ax.text(snr[-1]-0.3, tgt*1.8, "BER = 10⁻³",
                      color=_DIM, fontsize=8, ha="right")

        def snr_at(ber, t=1e-3):
            v = (ber > 0) & (ber < 1)
            if v.sum() < 2: return np.nan
            lb = np.log10(ber[v]); sv = snr[v]
            if t < 10**lb.min() or t > 10**lb.max(): return np.nan
            return float(np.interp(np.log10(t), lb[::-1], sv[::-1]))

        s_no = snr_at(ber_no); s_pso = snr_at(ber_pso); s_opt = snr_at(ber_opt)

        if not (np.isnan(s_no) or np.isnan(s_pso)):
            gain = s_no - s_pso
            self._ax.annotate("",
                xy=(s_pso, tgt*0.55), xytext=(s_no, tgt*0.55),
                arrowprops=dict(arrowstyle="<->", color=_CYAN, lw=2.0))
            self._ax.text((s_no+s_pso)/2, tgt*0.25, f"{gain:.1f} dB gain",
                          color=_CYAN, fontsize=10, ha="center", fontweight="bold")
            self._c_gain.set(f"{gain:.1f} dB")

        if not (np.isnan(s_pso) or np.isnan(s_opt)):
            self._c_gap.set(f"{s_pso - s_opt:.1f} dB")

        self._ax.set_xlabel("SNR (dB)", fontsize=12)
        self._ax.set_ylabel("Bit Error Rate (BER)", fontsize=12)
        self._ax.set_title(
            f"BER vs SNR  |  {elv}° Elevation  |  N = {N} RIS Elements",
            fontsize=13, fontweight="bold", pad=12)
        self._ax.set_ylim([1e-5, 1])
        self._ax.set_xlim([snr[0], snr[-1]])
        self._ax.legend(fontsize=10, facecolor=_PANEL, edgecolor=_LINE,
                        labelcolor=_TEXT, loc="lower left", framealpha=0.95)
        self._fig.tight_layout()
        self._cv.draw()

        self._pbar.setValue(100)
        self._stat.setText("Complete")
        self._btn.setEnabled(True)

    def _export(self):
        p, _ = QFileDialog.getSaveFileName(self, "Save", "BER_vs_SNR.png", "PNG (*.png)")
        if p: self._fig.savefig(p, dpi=150, bbox_inches="tight", facecolor=_PLT)


# ─────────────────────────────────────────────────────────────────────────────
#  Tab 2 — PSO Convergence
# ─────────────────────────────────────────────────────────────────────────────

class PSOTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._worker = None
        self._setup()

    def _setup(self):
        root = QHBoxLayout(self)
        root.setSpacing(0)
        root.setContentsMargins(0, 0, 0, 0)

        side, sl = _side_panel(self)

        sl.addWidget(_label("Parameters", 11, _DIM))
        sl.addWidget(_sep())

        self.s_N  = _SpinField("RIS Elements  N",    4, 128,  64)
        self.s_kd = _SpinField("Direct K-factor",   -10,  15,  -3, " dB")
        self.s_sw = _SpinField("Swarm size",          5,  50,  20)
        self.s_it = _SpinField("Max iterations",     10, 200, 100)
        self.s_tr = _SpinField("Independent trials",  1,  10,   5)
        for s in [self.s_N, self.s_kd, self.s_sw, self.s_it, self.s_tr]:
            sl.addWidget(s)

        sl.addWidget(_sep())

        self._btn = QPushButton("▶   Run Analysis")
        self._btn.clicked.connect(self._run)
        self._exp = QPushButton("Export PNG")
        self._exp.setProperty("flat", "true")
        self._exp.clicked.connect(self._export)
        sl.addWidget(self._btn)
        sl.addWidget(self._exp)

        sl.addWidget(_sep())
        self._c_opt = _Card("Optimal fitness",         "--", _GRN)
        self._c_pso = _Card("PSO mean final fitness",   "--", _BLUE)
        self._c_pct = _Card("Achieved (% of optimal)", "--", "#c084fc")
        for c in [self._c_opt, self._c_pso, self._c_pct]:
            sl.addWidget(c)

        sl.addStretch()
        self._stat = _label("Ready", 11, _DIM)
        sl.addWidget(self._stat)

        root.addWidget(side)

        pw = QWidget()
        pw.setStyleSheet(f"background:{_CARD};")
        pl = QVBoxLayout(pw)
        pl.setContentsMargins(20, 16, 20, 16)
        self._fig = Figure(figsize=(9, 6), facecolor=_PLT)
        self._ax  = self._fig.add_subplot(111)
        _ax_style(self._ax)
        self._cv  = _Canvas(self._fig, pw)
        pl.addWidget(self._cv)
        self._idle("PSO Convergence  |  Phase Optimisation")
        root.addWidget(pw)

    def _idle(self, title):
        self._ax.clear(); _ax_style(self._ax)
        self._ax.text(0.5, 0.5, "Press  ▶  Run Analysis  to begin",
                      transform=self._ax.transAxes, ha="center", va="center",
                      color=_DIM, fontsize=15, style="italic")
        self._ax.set_title(title)
        self._cv.draw()

    def _run(self):
        if self._worker and self._worker.isRunning(): return
        self._btn.setEnabled(False)
        self._stat.setText("Running...")

        self._p = dict(
            N=self.s_N.value(), k_d=float(self.s_kd.value()),
            k_r=10.0, sw=self.s_sw.value(),
            ni=self.s_it.value(), trials=self.s_tr.value())

        def _fn(N, k_d, k_r, sw, ni, trials):
            return sim_pso_conv(N, k_d, k_r, sw, ni, trials)

        self._worker = _Worker(_fn, self._p)
        self._worker.done.connect(self._done)
        self._worker.start()

    def _done(self, out):
        hists, opt = out
        N = self._p["N"]; ni = self._p["ni"]; tr = self._p["trials"]
        it = np.arange(1, ni+1)

        self._ax.clear(); _ax_style(self._ax)

        clrs = plt.cm.cool(np.linspace(0.2, 0.9, tr))
        for i, h in enumerate(hists):
            self._ax.plot(it, h, color=clrs[i], lw=1.5, alpha=0.65, label=f"Trial {i+1}")

        mean_h = hists.mean(0)
        self._ax.plot(it, mean_h, color=_BLUE, lw=3.0, label="Mean", zorder=5)
        self._ax.axhline(opt, color=_GRN, linestyle="--", lw=2.2,
                         label="Closed-form Optimal", zorder=6)
        self._ax.fill_between(it, hists.min(0), hists.max(0), color=_BLUE, alpha=0.07)
        self._ax.axhline(0.95*opt, color=_DIM, lw=0.9, linestyle=":")
        self._ax.text(it[-1], 0.95*opt*1.005, "  95 %",
                      color=_DIM, fontsize=8, va="bottom")

        self._ax.set_xlabel("Iteration", fontsize=12)
        self._ax.set_ylabel("|h_eff|²  (fitness)", fontsize=12)
        self._ax.set_title(
            f"PSO Convergence  |  N = {N}  |  {tr} Independent Trials",
            fontsize=13, fontweight="bold")
        self._ax.set_xlim([1, ni])
        self._ax.legend(fontsize=9, facecolor=_PANEL, edgecolor=_LINE,
                        labelcolor=_TEXT, loc="lower right", framealpha=0.95)
        self._fig.tight_layout(); self._cv.draw()

        pf = float(hists[:, -1].mean())
        self._c_opt.set(f"{opt:.3f}")
        self._c_pso.set(f"{pf:.3f}")
        self._c_pct.set(f"{100*pf/opt:.1f} %")
        self._stat.setText("Complete"); self._btn.setEnabled(True)

    def _export(self):
        p, _ = QFileDialog.getSaveFileName(self, "Save", "PSO_Convergence.png", "PNG (*.png)")
        if p: self._fig.savefig(p, dpi=150, bbox_inches="tight", facecolor=_PLT)


# ─────────────────────────────────────────────────────────────────────────────
#  Tab 3 — Elevation Sweep
# ─────────────────────────────────────────────────────────────────────────────

class ElevTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._worker = None
        self._setup()

    def _setup(self):
        root = QHBoxLayout(self)
        root.setSpacing(0)
        root.setContentsMargins(0, 0, 0, 0)

        side, sl = _side_panel(self)

        sl.addWidget(_label("Parameters", 11, _DIM))
        sl.addWidget(_sep())

        self.s_N    = _SpinField("RIS Elements  N",   4, 128,  64)
        self.s_snr  = _SpinField("Fixed SNR",        -10,  30,   5, " dB")
        self.s_mc   = _SpinField("Monte Carlo trials", 10, 200,  40)
        self.s_bits = _SpinField("Bits per trial",    50, 500, 200)
        self.s_sw   = _SpinField("PSO swarm size",     5,  50,  20)
        self.s_it   = _SpinField("PSO iterations",     5,  80,  30)
        for s in [self.s_N, self.s_snr, self.s_mc, self.s_bits, self.s_sw, self.s_it]:
            sl.addWidget(s)

        sl.addWidget(_sep())
        self._btn = QPushButton("▶   Run Sweep")
        self._btn.clicked.connect(self._run)
        self._exp = QPushButton("Export PNG")
        self._exp.setProperty("flat", "true")
        self._exp.clicked.connect(self._export)
        sl.addWidget(self._btn)
        sl.addWidget(self._exp)

        sl.addWidget(_sep())
        self._c_no  = _Card("BER at 10°  — No RIS",  "--", _RED)
        self._c_pso = _Card("BER at 10°  — PSO-RIS", "--", _BLUE)
        self._c_win = _Card("Extra usable window",    "--", _CYAN)
        for c in [self._c_no, self._c_pso, self._c_win]:
            sl.addWidget(c)

        sl.addStretch()
        self._pbar = QProgressBar(); self._pbar.setValue(0)
        self._stat = _label("Ready", 11, _DIM)
        sl.addWidget(self._pbar); sl.addWidget(self._stat)

        root.addWidget(side)

        pw = QWidget()
        pw.setStyleSheet(f"background:{_CARD};")
        pl = QVBoxLayout(pw)
        pl.setContentsMargins(20, 16, 20, 16)
        self._fig = Figure(figsize=(9, 6), facecolor=_PLT)
        self._ax  = self._fig.add_subplot(111)
        _ax_style(self._ax)
        self._cv  = _Canvas(self._fig, pw)
        pl.addWidget(self._cv)
        self._idle("BER vs Elevation  |  Satellite Pass Analysis")
        root.addWidget(pw)

    def _idle(self, title):
        self._ax.clear(); _ax_style(self._ax)
        self._ax.text(0.5, 0.5, "Press  ▶  Run Sweep  to begin",
                      transform=self._ax.transAxes, ha="center", va="center",
                      color=_DIM, fontsize=15, style="italic")
        self._ax.set_title(title); self._cv.draw()

    def _run(self):
        if self._worker and self._worker.isRunning(): return
        self._btn.setEnabled(False)
        self._pbar.setValue(0); self._stat.setText("Running...")

        self._p = dict(
            elev_arr=np.arange(5, 90, 5, dtype=float),
            snr=float(self.s_snr.value()), N=self.s_N.value(),
            n_mc=self.s_mc.value(), n_bits=self.s_bits.value(),
            sw=self.s_sw.value(), ni=self.s_it.value())

        def _fn(elev_arr, snr, N, n_mc, n_bits, sw, ni):
            return sim_elev(elev_arr, snr, N, n_mc, n_bits, sw, ni,
                            cb=self._worker.progress.emit)

        self._worker = _Worker(_fn, self._p)
        self._worker.progress.connect(self._pbar.setValue)
        self._worker.done.connect(self._done)
        self._worker.start()

    def _done(self, out):
        elev = self._p["elev_arr"]; N = self._p["N"]; snr = self._p["snr"]
        ber_no, ber_pso = out

        self._ax.clear(); _ax_style(self._ax)
        kw = dict(linewidth=2.5, markersize=7)
        self._ax.semilogy(elev, ber_no,  color=_RED,  marker="o", label="No RIS", **kw)
        self._ax.semilogy(elev, ber_pso, color=_BLUE, marker="^", label="PSO-Optimised RIS", **kw)

        thr = 1e-2
        self._ax.axhline(thr, color=_DIM, linestyle=":", lw=1.1)
        self._ax.text(elev[-1], thr*1.6, "  Usability threshold",
                      color=_DIM, fontsize=8, ha="right")
        self._ax.fill_between(elev, ber_no, ber_pso,
                              where=ber_no > ber_pso,
                              color=_BLUE, alpha=0.08, label="Improvement region")

        self._ax.set_xlabel("Elevation Angle (°)", fontsize=12)
        self._ax.set_ylabel("Bit Error Rate (BER)", fontsize=12)
        self._ax.set_title(
            f"BER vs Elevation  |  SNR = {snr:.0f} dB  |  N = {N}",
            fontsize=13, fontweight="bold")
        self._ax.set_ylim([1e-5, 1])
        self._ax.legend(fontsize=10, facecolor=_PANEL, edgecolor=_LINE,
                        labelcolor=_TEXT, framealpha=0.95)
        self._fig.tight_layout(); self._cv.draw()

        i10 = int(np.argmin(np.abs(elev - 10)))
        self._c_no.set(f"{ber_no[i10]:.3f}")
        self._c_pso.set(f"{ber_pso[i10]:.3f}")
        u_no  = float(np.sum(ber_no  < thr)) / len(elev) * (elev[-1] - elev[0])
        u_pso = float(np.sum(ber_pso < thr)) / len(elev) * (elev[-1] - elev[0])
        self._c_win.set(f"+{max(u_pso - u_no, 0):.0f}°")

        self._pbar.setValue(100); self._stat.setText("Complete")
        self._btn.setEnabled(True)

    def _export(self):
        p, _ = QFileDialog.getSaveFileName(self, "Save", "BER_vs_Elevation.png", "PNG (*.png)")
        if p: self._fig.savefig(p, dpi=150, bbox_inches="tight", facecolor=_PLT)


# ─────────────────────────────────────────────────────────────────────────────
#  Tab 4 — System diagram (animated)
# ─────────────────────────────────────────────────────────────────────────────

class DiagramTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._angle = 50.0
        self._timer = None
        self._setup()

    def _setup(self):
        lay = QVBoxLayout(self)
        lay.setContentsMargins(32, 20, 32, 16)
        lay.setSpacing(10)

        hdr = _label("System Architecture  —  RIS-Assisted Satellite Link",
                     15, "#dde8f8", bold=True)
        hdr.setAlignment(Qt.AlignmentFlag.AlignCenter)
        sub = _label(
            "Direct path  h_d  (weak at low elevation)  +"
            "  Sat→RIS  h_sr  →  PSO phase optimisation  →  RIS→GS  h_rg",
            11, _DIM)
        sub.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lay.addWidget(hdr)
        lay.addWidget(sub)

        self._fig = Figure(figsize=(12, 5.2), facecolor=_PLT)
        self._ax  = self._fig.add_subplot(111)
        self._cv  = _Canvas(self._fig, self)
        lay.addWidget(self._cv)

        leg = QHBoxLayout()
        leg.setSpacing(28)
        for txt, clr in [("Direct  h_d", _RED), ("Sat→RIS  h_sr", _BLUE),
                         ("RIS→GS  h_rg", _GRN), ("RIS panel", _CYAN),
                         ("Orbit track", _DIM)]:
            l = QLabel(f"●  {txt}")
            l.setStyleSheet(f"color:{clr}; font-size:11px; background:transparent;")
            leg.addWidget(l)
        leg.addStretch()
        lay.addLayout(leg)

    def showEvent(self, e):
        super().showEvent(e)
        if self._timer is None:
            self._timer = QTimer(self)
            self._timer.timeout.connect(self._tick)
            self._timer.start(180)

    def _tick(self):
        self._angle = (self._angle + 2.5) % 360
        self._draw()

    def _draw(self):
        ax = self._ax
        ax.clear()
        ax.set_facecolor(_PLT)
        ax.set_xlim(0, 10); ax.set_ylim(0, 6)
        ax.axis("off")
        ax.set_title("CubeSat  →  RIS Panel  →  Ground Station  |  Live Signal Flow",
                     color="#dde8f8", fontsize=12, fontweight="bold", pad=8)

        # earth
        ax.add_patch(mpatches.FancyBboxPatch(
            (0, 0), 10, 0.6, boxstyle="square",
            facecolor="#0a1e0e", edgecolor="#183020", lw=1))
        ax.text(5, 0.3, "Earth Surface", color="#2d5a35",
                ha="center", va="center", fontsize=8.5)

        # ground station
        gx, gy = 1.3, 0.6
        ax.add_patch(mpatches.Arc((gx, gy+0.30), 0.72, 0.36, 0, 0, 180,
                                   color="#7080a0", lw=2))
        ax.plot([gx, gx],           [gy, gy+0.30], color="#7080a0", lw=3)
        ax.plot([gx-0.26, gx+0.26], [gy, gy],      color="#7080a0", lw=3)
        ax.text(gx, 0.08, "Ground Station", color="#7080a0", ha="center", fontsize=8.5)

        # RIS panel
        rx, ry = 4.85, 1.65
        for row in range(7):
            for col in range(4):
                ax.add_patch(mpatches.Rectangle(
                    (rx - 0.12 + col*0.095, ry + row*0.11),
                    0.075, 0.088,
                    facecolor=_CYAN, edgecolor=_PLT, alpha=0.87))
        ax.add_patch(mpatches.Ellipse(
            (rx+0.13, ry+0.36), 0.68, 0.88,
            edgecolor=_CYAN, facecolor="none", lw=1.1, alpha=0.3, linestyle="--"))
        ax.text(rx+0.13, ry-0.25, "RIS Panel  (N elements)",
                color=_CYAN, ha="center", fontsize=8.5, fontweight="bold")

        # satellite
        cx, cy, orx, ory = 5.0, 4.55, 3.35, 1.32
        th_arr = np.linspace(0, 2*np.pi, 300)
        ax.plot(cx + orx*np.cos(th_arr), cy + ory*np.sin(th_arr),
                "--", color=_DIM, lw=0.7, alpha=0.4)

        th = np.deg2rad(self._angle)
        sx, sy = cx + orx*np.cos(th), cy + ory*np.sin(th)

        ax.add_patch(mpatches.FancyBboxPatch(
            (sx-0.15, sy-0.085), 0.30, 0.17,
            boxstyle="round,pad=0.02",
            facecolor="#b4c5dc", edgecolor="#6888a8", lw=1.4))
        for dx in [(-0.36, -0.15), (0.15, 0.36)]:
            ax.plot([sx+dx[0], sx+dx[1]], [sy, sy],
                    color=_BLUE, lw=4.5, solid_capstyle="round")
        ax.text(sx, sy+0.24, "CubeSat", color="#b4c5dc", ha="center", fontsize=8)

        # signal paths
        ax.annotate("", xy=(gx+0.28, gy+0.47), xytext=(sx, sy),
            arrowprops=dict(arrowstyle="->", color=_RED, lw=1.8,
                            connectionstyle="arc3,rad=-0.12", linestyle="dashed"))
        ax.text((gx+sx)/2 - 0.5, (gy+sy)/2 + 0.5,
                "h_d  (weak)", color=_RED, fontsize=8.5,
                ha="center", style="italic", rotation=-22)

        ax.annotate("", xy=(rx+0.13, ry+0.78), xytext=(sx, sy-0.09),
            arrowprops=dict(arrowstyle="->", color=_BLUE, lw=2.1,
                            connectionstyle="arc3,rad=0.08"))
        ax.text((rx+sx)/2+0.18, (ry+sy)/2+0.1,
                "h_sr", color=_BLUE, fontsize=11, fontweight="bold")

        ax.annotate("", xy=(gx+0.20, gy+0.50), xytext=(rx-0.13, ry+0.34),
            arrowprops=dict(arrowstyle="->", color=_GRN, lw=2.1,
                            connectionstyle="arc3,rad=-0.08"))
        ax.text((rx+gx)/2, (ry+gy)/2+0.52,
                "h_rg", color=_GRN, fontsize=11, fontweight="bold")

        # formula
        ax.text(7.9, 5.2,
                "θ*ₙ = ∠h_d − ∠(h_sr,ₙ · h_rg,ₙ)",
                color=_CYAN, fontsize=9.5, ha="center", va="center",
                fontfamily="monospace",
                bbox=dict(boxstyle="round,pad=0.45",
                          facecolor=_PANEL, edgecolor=_LINE, alpha=0.92))
        ax.text(7.9, 4.72, "Optimal Phase Formula",
                color=_DIM, fontsize=8, ha="center")

        # coherent note
        ax.annotate("Coherent\nCombining",
            xy=(gx+0.33, gy+0.62), xytext=(2.9, 2.65),
            fontsize=8, color="#a78bfa", ha="center",
            arrowprops=dict(arrowstyle="->", color="#a78bfa", lw=1.1))

        # elevation arc
        raw = float(np.rad2deg(np.arctan2(max(sy-gy, 0.01), max(sx-gx, 0.01))))
        ea  = max(5.0, min(85.0, raw))
        ax.add_patch(mpatches.Arc(
            (gx, gy+0.28), 0.95, 0.95, 0, 0, ea,
            color=_AMB, lw=1.5, alpha=0.7))
        ax.text(gx+0.68, gy+0.55, f"{ea:.0f}°",
                color=_AMB, fontsize=8.5, fontstyle="italic")

        self._fig.tight_layout()
        self._cv.draw()


# ─────────────────────────────────────────────────────────────────────────────
#  Main window
# ─────────────────────────────────────────────────────────────────────────────

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("RIS-Assisted CubeSat Simulator  |  APache  |  MCTE")
        self.setMinimumSize(1280, 800)
        self.setStyleSheet(CSS)
        self._build()

    def _build(self):
        cw = QWidget()
        self.setCentralWidget(cw)
        root = QVBoxLayout(cw)
        root.setSpacing(0)
        root.setContentsMargins(0, 0, 0, 0)

        # header
        hdr = QWidget()
        hdr.setFixedHeight(68)
        hdr.setStyleSheet(
            f"background:qlineargradient(x1:0,y1:0,x2:1,y2:0,"
            f"stop:0 {_BG},stop:0.5 #0c1e38,stop:1 {_BG});"
            f"border-bottom:1px solid {_LINE};")
        hl = QHBoxLayout(hdr)
        hl.setContentsMargins(28, 0, 28, 0)
        hl.setSpacing(16)

        icon = QLabel("⬡")
        icon.setStyleSheet(f"color:{_CYAN}; font-size:34px; background:transparent;")
        hl.addWidget(icon)

        tb = QVBoxLayout(); tb.setSpacing(2)
        t1 = QLabel("RIS-Assisted CubeSat-to-Ground Communication")
        t1.setStyleSheet("color:#eef4ff; font-size:16px; font-weight:700; background:transparent;")
        t2 = QLabel("Particle Swarm Optimisation  ·  Rician Fading  ·  Monte Carlo BER  |  MCTE / APache")
        t2.setStyleSheet(f"color:{_DIM}; font-size:11px; background:transparent;")
        tb.addWidget(t1); tb.addWidget(t2)
        hl.addLayout(tb)
        hl.addStretch()

        for lbl, val in [("Algorithm", "PSO + Optimal"), ("Channel", "Rician"), ("Modulation", "BPSK")]:
            pill = QWidget()
            pill.setStyleSheet(
                f"background:{_PANEL}; border:1px solid {_LINE}; border-radius:6px;")
            pl2 = QVBoxLayout(pill); pl2.setContentsMargins(12, 5, 12, 5); pl2.setSpacing(0)
            v = QLabel(val); l = QLabel(lbl)
            v.setStyleSheet(f"color:{_CYAN}; font-size:11px; font-weight:700; background:transparent;")
            l.setStyleSheet(f"color:{_DIM}; font-size:10px; background:transparent;")
            pl2.addWidget(v); pl2.addWidget(l)
            hl.addWidget(pill)

        root.addWidget(hdr)

        tabs = QTabWidget()
        tabs.addTab(BERTab(),     "  BER vs SNR  ")
        tabs.addTab(PSOTab(),     "  PSO Convergence  ")
        tabs.addTab(ElevTab(),    "  Elevation Sweep  ")
        tabs.addTab(DiagramTab(), "  System Diagram  ")
        root.addWidget(tabs)


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    win = MainWindow()
    win.show()
    app.exec()


if __name__ == "__main__":
    main()
