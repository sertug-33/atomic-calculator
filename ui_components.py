import math

from PySide6.QtCore import (Property, QPropertyAnimation, QEasingCurve, Qt,
                            QPointF, QRectF, QEvent, Signal)
from PySide6.QtWidgets import (QWidget, QSizePolicy, QComboBox, QCompleter,
                               QVBoxLayout, QLabel, QPushButton, QScrollArea,
                               QButtonGroup, QFrame, QApplication)
from PySide6.QtGui import (QPainter, QPen, QBrush, QColor, QFont, QFontMetrics,
                           QPolygonF, QRadialGradient, QLinearGradient)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.ticker import MultipleLocator

FONT_DIN   = "'DIN Alternate', 'DIN Condensed', 'Helvetica Neue', sans-serif"
FONT_PANEL = "'Helvetica Neue', 'Helvetica', 'Arial', sans-serif"
FONT_DATA  = "'Andale Mono', 'PT Mono', 'Menlo', monospace"
FONT_PRINT = "'Courier New', 'Courier', monospace"
FONT_BOOK  = "'Palatino', 'Palatino Linotype', 'Iowan Old Style', serif"

PALETTE_DARK = {
    "chassis":     "#141715",
    "panel":       "#1d2120",
    "panel_alt":   "#242927",
    "tab_bg":      "#191d1c",
    "seam":        "#0b0e0d",
    "bevel":       "#3d453f",
    "inset":       "#07090a",
    "text":        "#d6ddd4",
    "label":       "#8e9a92",
    "faint":       "#5d6762",
    "accent":      "#f2a33c",
    "readout":     "#f2a33c",
    "accent_hi":   "#ffcd87",
    "accent_dim":  "#a86f22",
    "green":       "#86c98f",
    "red":         "#e2614b",
    "blue":        "#4fa8b2",
    "sel_bg":      "#3a4440",
    "sel_fg":      "#ffcd87",
    "grid_major":  "#2b3431",
    "grid_minor":  "#1e2523",
    "frame":       "#414a45",
}

PALETTE_LIGHT = {
    "chassis":     "#e4ddcc",
    "panel":       "#f4efe3",
    "panel_alt":   "#eae3d3",
    "tab_bg":      "#dcd4c1",
    "seam":        "#b6ab94",
    "bevel":       "#fffdf6",
    "inset":       "#fbf8f0",
    "text":        "#1f2528",
    "label":       "#6b6456",
    "faint":       "#928a79",
    "accent":      "#b2341f",
    "readout":     "#1b4965",
    "accent_hi":   "#8d2413",
    "accent_dim":  "#c8846f",
    "green":       "#2d6a4f",
    "red":         "#a8321c",
    "blue":        "#1b4965",
    "sel_bg":      "#e0d6bf",
    "sel_fg":      "#8d2413",
    "grid_major":  "#c98f6d",
    "grid_minor":  "#e0b79c",
    "frame":       "#8d8372",
}

def palette(is_dark: bool) -> dict:
    return PALETTE_DARK if is_dark else PALETTE_LIGHT

_QSS = """
QMainWindow, QDialog { background-color: %(chassis)s; }

QTabWidget::pane {
    background: %(panel)s;
    border: 1px solid %(seam)s;
    border-top: 2px solid %(accent_dim)s;
    border-radius: 0px;
    top: -1px;
}
QTabBar { qproperty-drawBase: 0; }
QTabBar::tab {
    background: %(tab_bg)s;
    color: %(label)s;
    min-width: 108px;
    padding: 8px 12px;
    margin-right: 2px;
    border: 1px solid %(seam)s;
    border-bottom: none;
    border-radius: 0px;
    font-family: %(font_din)s;
    font-weight: bold;
    font-size: 11px;
    letter-spacing: 1.6px;
}
QTabBar::tab:selected {
    background: %(panel)s;
    color: %(accent)s;
    border-top: 2px solid %(accent)s;
    border-left: 1px solid %(bevel)s;
    border-right: 1px solid %(seam)s;
}
QTabBar::tab:hover:!selected { background: %(panel_alt)s; color: %(text)s; }

QLabel {
    color: %(text)s;
    font-family: %(font_panel)s;
    font-size: 12px;
}
QLabel#Engraved {
    color: %(label)s;
    font-family: %(font_din)s;
    font-weight: bold;
    font-size: 10px;
    letter-spacing: 2.0px;
}
QLabel#Readout {
    color: %(readout)s;
    font-family: %(font_data)s;
    font-size: 17px;
    font-weight: bold;
}
QLabel#Unit {
    color: %(faint)s;
    font-family: %(font_data)s;
    font-size: 10px;
    letter-spacing: 0.6px;
}

QLineEdit {
    background-color: %(inset)s;
    color: %(readout)s;
    selection-background-color: %(sel_bg)s;
    selection-color: %(sel_fg)s;
    border: 1px solid %(seam)s;
    border-bottom: 1px solid %(bevel)s;
    border-right: 1px solid %(bevel)s;
    border-radius: 0px;
    padding: 8px 10px;
    font-family: %(font_data)s;
    font-size: 14px;
    letter-spacing: 0.4px;
}
QLineEdit:focus { border: 1px solid %(accent_dim)s; background-color: %(inset)s; }

QComboBox {
    background-color: %(panel_alt)s;
    color: %(text)s;
    border-top: 1px solid %(bevel)s;
    border-left: 1px solid %(bevel)s;
    border-right: 1px solid %(seam)s;
    border-bottom: 1px solid %(seam)s;
    border-radius: 0px;
    padding: 7px 10px;
    font-family: %(font_panel)s;
    font-size: 12px;
}
QComboBox:hover { background-color: %(panel)s; }
QComboBox::drop-down {
    subcontrol-origin: padding;
    subcontrol-position: top right;
    width: 24px;
    border-left: 1px solid %(seam)s;
    background: %(tab_bg)s;
}
QComboBox::down-arrow {
    image: none;
    width: 0px; height: 0px;
    border: none;
    background: transparent;
}
QComboBox QAbstractItemView {
    background-color: %(panel)s;
    color: %(text)s;
    border: 1px solid %(seam)s;
    border-radius: 0px;
    selection-background-color: %(sel_bg)s;
    selection-color: %(sel_fg)s;
    outline: none;
    padding: 2px;
    font-family: %(font_data)s;
    font-size: 12px;
}
QComboBox QAbstractItemView::item { min-height: 26px; padding: 3px 8px; }

QPushButton {
    background-color: %(panel_alt)s;
    color: %(text)s;
    border-top: 1px solid %(bevel)s;
    border-left: 1px solid %(bevel)s;
    border-right: 1px solid %(seam)s;
    border-bottom: 1px solid %(seam)s;
    border-radius: 0px;
    padding: 9px 18px;
    font-family: %(font_din)s;
    font-weight: bold;
    font-size: 12px;
    letter-spacing: 1.4px;
}
QPushButton:hover { background-color: %(panel)s; color: %(text)s; }
QPushButton:pressed {
    background-color: %(chassis)s;
    border-top: 1px solid %(seam)s;
    border-left: 1px solid %(seam)s;
    border-right: 1px solid %(bevel)s;
    border-bottom: 1px solid %(bevel)s;
    padding: 10px 17px 8px 19px;
}
QPushButton#PrimaryLever {
    background-color: %(panel_alt)s;
    color: %(accent)s;
    border-left: 3px solid %(accent_dim)s;
    font-size: 13px;
    padding: 12px 18px;
}
QPushButton#PrimaryLever:hover { border-left: 3px solid %(accent)s; }

QTextEdit {
    background-color: %(inset)s;
    color: %(readout)s;
    border: 1px solid %(seam)s;
    border-bottom: 1px solid %(bevel)s;
    border-right: 1px solid %(bevel)s;
    border-radius: 0px;
    padding: 14px 16px;
    font-family: %(font_print)s;
    font-size: 12px;
    selection-background-color: %(sel_bg)s;
    selection-color: %(sel_fg)s;
}

QTableWidget {
    background-color: %(inset)s;
    alternate-background-color: %(panel)s;
    color: %(text)s;
    gridline-color: %(grid_minor)s;
    border: 1px solid %(seam)s;
    border-radius: 0px;
    font-family: %(font_data)s;
    font-size: 12px;
}
QTableWidget::item { padding: 3px 6px; }
QTableWidget::item:selected { background: %(sel_bg)s; color: %(sel_fg)s; }
QHeaderView::section {
    background-color: %(panel_alt)s;
    color: %(label)s;
    padding: 7px 6px;
    border: none;
    border-right: 1px solid %(seam)s;
    border-bottom: 2px solid %(accent_dim)s;
    font-family: %(font_din)s;
    font-weight: bold;
    font-size: 10px;
    letter-spacing: 1.4px;
}
QTableCornerButton::section { background: %(panel_alt)s; border: none; }

QSlider::groove:horizontal {
    height: 4px;
    background: %(inset)s;
    border-top: 1px solid %(seam)s;
    border-bottom: 1px solid %(bevel)s;
    border-radius: 0px;
}
QSlider::sub-page:horizontal { background: %(accent_dim)s; }
QSlider::handle:horizontal {
    background: %(panel_alt)s;
    border-top: 1px solid %(bevel)s;
    border-left: 1px solid %(bevel)s;
    border-right: 1px solid %(seam)s;
    border-bottom: 1px solid %(seam)s;
    width: 12px;
    height: 26px;
    margin: -12px 0;
    border-radius: 0px;
}
QSlider::handle:horizontal:hover { background: %(panel)s; }

QStatusBar {
    background-color: %(chassis)s;
    color: %(label)s;
    border-top: 1px solid %(seam)s;
    padding: 3px 12px;
    font-family: %(font_data)s;
    font-size: 11px;
    letter-spacing: 0.6px;
}
QStatusBar::item { border: none; }

QScrollBar:vertical   { background: %(chassis)s; width: 11px;  margin: 0; }
QScrollBar:horizontal { background: %(chassis)s; height: 11px; margin: 0; }
QScrollBar::handle:vertical, QScrollBar::handle:horizontal {
    background: %(frame)s;
    min-height: 24px; min-width: 24px;
    border-radius: 0px;
    border: 1px solid %(seam)s;
}
QScrollBar::handle:vertical:hover, QScrollBar::handle:horizontal:hover {
    background: %(label)s;
}
QScrollBar::add-line, QScrollBar::sub-line { width: 0px; height: 0px; border: none; }
QScrollBar::add-page, QScrollBar::sub-page { background: none; }
"""

def build_qss(is_dark: bool) -> str:
    p = dict(palette(is_dark))
    p.update({
        "font_din": FONT_DIN, "font_panel": FONT_PANEL,
        "font_data": FONT_DATA, "font_print": FONT_PRINT, "font_book": FONT_BOOK,
    })
    return _QSS % p

DARK_THEME_QSS = build_qss(True)
LIGHT_THEME_QSS = build_qss(False)

def panel_qss(is_dark: bool, accent_rule: bool = True) -> str:
    p = palette(is_dark)
    top = f"border-top: 2px solid {p['accent_dim']};" if accent_rule \
        else f"border-top: 1px solid {p['bevel']};"
    return (
        "QFrame {"
        f" background-color: {p['panel']};"
        f" border-left: 1px solid {p['bevel']};"
        f" border-right: 1px solid {p['seam']};"
        f" border-bottom: 1px solid {p['seam']};"
        f" {top}"
        " border-radius: 0px;"
        " padding: 16px;"
        "}"
    )

class PanelComboBox(QComboBox):
    def __init__(self, parent=None, is_dark=True):
        super().__init__(parent)
        self.is_dark = is_dark

    def apply_theme(self, is_dark: bool):
        self.is_dark = is_dark
        self.update()

    def paintEvent(self, ev):
        super().paintEvent(ev)
        p = palette(self.is_dark)
        q = QPainter(self)
        q.setRenderHint(QPainter.Antialiasing)
        cx = self.width() - 12
        cy = self.height() / 2
        tri = QPolygonF([QPointF(cx - 4.5, cy - 2.2),
                         QPointF(cx + 4.5, cy - 2.2),
                         QPointF(cx, cy + 3.2)])
        q.setPen(Qt.NoPen)
        q.setBrush(QBrush(QColor(p["accent"])))
        q.drawPolygon(tri)

class PanelDivider(QWidget):
    def __init__(self, parent=None, is_dark=True):
        super().__init__(parent)
        self.is_dark = is_dark
        self.setFixedHeight(2)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)

    def apply_theme(self, is_dark: bool):
        self.is_dark = is_dark
        self.update()

    def paintEvent(self, _):
        p = palette(self.is_dark)
        q = QPainter(self)
        q.setPen(QPen(QColor(p["seam"]), 1))
        q.drawLine(0, 0, self.width(), 0)
        q.setPen(QPen(QColor(p["bevel"]), 1))
        q.drawLine(0, 1, self.width(), 1)

class IndicatorLamp(QWidget):
    LAMP_D = 15
    TEXT_X = 33

    def __init__(self, caption="", parent=None, is_dark=True):
        super().__init__(parent)
        self.caption = caption
        self.state = "idle"
        self.is_dark = is_dark
        self.setFixedHeight(30)
        self._refit()

    def _refit(self):
        f = QFont("DIN Alternate", 10, QFont.Bold)
        f.setLetterSpacing(QFont.AbsoluteSpacing, 1.5)
        self.setMinimumWidth(self.TEXT_X + QFontMetrics(f).horizontalAdvance(self.caption) + 10)
        self.updateGeometry()

    def apply_theme(self, is_dark: bool):
        self.is_dark = is_dark
        self.update()

    def set_state(self, state: str, caption: str | None = None):
        self.state = state
        if caption is not None and caption != self.caption:
            self.caption = caption
            self._refit()
        self.update()

    def paintEvent(self, _):
        p = palette(self.is_dark)
        q = QPainter(self)
        q.setRenderHint(QPainter.Antialiasing)

        d = self.LAMP_D
        cy = self.height() / 2
        cx = 4 + d / 2
        rect = QRectF(4, cy - d / 2, d, d)

        base = {"ok": p["green"], "fail": p["red"]}.get(self.state, p["faint"])
        col = QColor(base)
        lit = self.state in ("ok", "fail")

        grad = QRadialGradient(QPointF(cx - 2, cy - 2), d * 0.85)
        if lit:
            grad.setColorAt(0.0, col.lighter(165))
            grad.setColorAt(0.55, col)
            grad.setColorAt(1.0, col.darker(230))
        else:
            grad.setColorAt(0.0, QColor(p["panel_alt"]))
            grad.setColorAt(1.0, QColor(p["seam"]))
        q.setBrush(QBrush(grad))
        q.setPen(Qt.NoPen)
        q.drawEllipse(rect)

        q.setBrush(Qt.NoBrush)
        q.setPen(QPen(QColor(p["seam"]), 1.4))
        q.drawArc(rect.adjusted(-1, -1, 1, 1), 180 * 16, 180 * 16)
        q.setPen(QPen(QColor(p["bevel"]), 1.4))
        q.drawArc(rect.adjusted(-1, -1, 1, 1), 0, 180 * 16)

        f = QFont("DIN Alternate", 10, QFont.Bold)
        f.setLetterSpacing(QFont.AbsoluteSpacing, 1.5)
        q.setFont(f)
        q.setPen(QColor(col if lit else QColor(p["label"])))
        q.drawText(QRectF(self.TEXT_X, 0, self.width() - self.TEXT_X - 4, self.height()),
                   Qt.AlignVCenter | Qt.AlignLeft, self.caption)

class AnalogGauge(QWidget):
    def __init__(self, caption="ÖLÇÜM", unit="%", vmax=100.0, warn=None, parent=None, is_dark=True):
        super().__init__(parent)
        self.caption, self.unit, self.vmax = caption, unit, float(vmax)
        self.warn = warn
        self.is_dark = is_dark
        self._value = 0.0
        self._target = 0.0
        self.setMinimumSize(210, 150)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.setFixedHeight(172)

        self.anim = QPropertyAnimation(self, b"needle")
        self.anim.setDuration(620)
        self.anim.setEasingCurve(QEasingCurve.OutBack)

    def get_needle(self): return self._value
    def set_needle(self, v): self._value = v; self.update()
    needle = Property(float, get_needle, set_needle)

    def apply_theme(self, is_dark: bool): self.is_dark = is_dark; self.update()

    def set_value(self, v: float, animate: bool = True):
        self._target = max(0.0, min(self.vmax, float(v)))
        if animate:
            self.anim.stop()
            self.anim.setStartValue(self._value)
            self.anim.setEndValue(self._target)
            self.anim.start()
        else:
            self.set_needle(self._target)

    def _angle(self, v):
        frac = 0.0 if self.vmax == 0 else max(0.0, min(1.0, v / self.vmax))
        return math.radians(-60 + 120 * frac)

    def paintEvent(self, _):
        p = palette(self.is_dark)
        q = QPainter(self)
        q.setRenderHint(QPainter.Antialiasing)

        w, h = self.width(), self.height()
        q.fillRect(self.rect(), QColor(p["inset"]))
        q.setPen(QPen(QColor(p["seam"]), 1))
        q.drawLine(0, 0, w, 0); q.drawLine(0, 0, 0, h)
        q.setPen(QPen(QColor(p["bevel"]), 1))
        q.drawLine(w - 1, 0, w - 1, h); q.drawLine(0, h - 1, w, h - 1)

        cx, cy = w / 2, h - 50
        R = min(w * 0.42, h * 0.62)

        q.setPen(QPen(QColor(p["grid_minor"]), 6, Qt.SolidLine, Qt.FlatCap))
        q.drawArc(QRectF(cx - R * 0.86, cy - R * 0.86, R * 1.72, R * 1.72), int(30 * 16), int(120 * 16))

        q.setPen(QPen(QColor(p["frame"]), 1.3, Qt.SolidLine, Qt.FlatCap))
        q.drawArc(QRectF(cx - R, cy - R, R * 2, R * 2), int(30 * 16), int(120 * 16))

        if self.warn is not None and 0 <= self.warn < self.vmax:
            warn_deg = 150.0 - 120.0 * (self.warn / self.vmax)
            q.setPen(QPen(QColor(p["red"]), 4, Qt.SolidLine, Qt.FlatCap))
            q.drawArc(QRectF(cx - R, cy - R, R * 2, R * 2), int(30 * 16), int((warn_deg - 30.0) * 16))

        fnum = QFont("DIN Alternate", 9, QFont.Bold)
        q.setFont(fnum)
        for i in range(21):
            v = self.vmax * i / 20.0
            a = self._angle(v)
            major = (i % 5 == 0)
            ln = R * (0.16 if major else 0.08)
            sx, sy = cx + math.sin(a) * R, cy - math.cos(a) * R
            ex, ey = cx + math.sin(a) * (R - ln), cy - math.cos(a) * (R - ln)
            q.setPen(QPen(QColor(p["text"] if major else p["label"]), 1.6 if major else 0.9))
            q.drawLine(QPointF(sx, sy), QPointF(ex, ey))
            if major:
                tx = cx + math.sin(a) * (R - ln - 13)
                ty = cy - math.cos(a) * (R - ln - 13)
                q.setPen(QColor(p["label"]))
                q.drawText(QRectF(tx - 16, ty - 8, 32, 16), Qt.AlignCenter, f"{v:g}")

        a = self._angle(self._value)
        tipx, tipy = cx + math.sin(a) * (R - 6), cy - math.cos(a) * (R - 6)
        tailx, taily = cx - math.sin(a) * (R * 0.16), cy + math.cos(a) * (R * 0.16)
        needle_col = QColor(p["red"] if (self.warn is not None and self._value >= self.warn) else p["accent"])
        q.setPen(QPen(needle_col, 2.0, Qt.SolidLine, Qt.RoundCap))
        q.drawLine(QPointF(tailx, taily), QPointF(tipx, tipy))
        q.setBrush(QBrush(QColor(p["frame"])))
        q.setPen(QPen(QColor(p["seam"]), 1))
        q.drawEllipse(QPointF(cx, cy), 5, 5)

        fc = QFont("DIN Alternate", 10, QFont.Bold)
        fc.setLetterSpacing(QFont.AbsoluteSpacing, 1.8)
        q.setFont(fc)
        q.setPen(QColor(p["label"]))
        q.drawText(QRectF(0, 6, w, 14), Qt.AlignCenter, self.caption)

        q.setFont(QFont("Andale Mono", 13, QFont.Bold))
        q.setPen(needle_col)
        q.drawText(QRectF(0, h - 28, w, 20), Qt.AlignCenter, f"{self._target:.2f} {self.unit}")

class SemfBarMeter(QWidget):
    GUTTER = 178
    VALCOL = 84
    ROW = 26

    def __init__(self, parent=None, is_dark=True):
        super().__init__(parent)
        self.is_dark = is_dark
        self.terms = []
        self.total = 0.0
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.setMinimumWidth(460)
        self.setFixedHeight(self.ROW * 5 + 64)

    def apply_theme(self, is_dark: bool):
        self.is_dark = is_dark
        self.update()

    def set_terms(self, terms, total=0.0):
        self.terms = list(terms)
        self.total = total
        self.update()

    @staticmethod
    def _nice_step(span, target=8):
        if span <= 0: return 1.0
        raw = span / target
        mag = 10 ** math.floor(math.log10(raw))
        for m in (1, 2, 5, 10):
            if raw <= m * mag: return m * mag
        return 10 * mag

    def paintEvent(self, _):
        if not self.terms: return
        p = palette(self.is_dark)
        q = QPainter(self)
        q.setRenderHint(QPainter.Antialiasing)

        w, h = self.width(), self.height()
        q.fillRect(self.rect(), QColor(p["inset"]))
        q.setPen(QPen(QColor(p["seam"]), 1))
        q.drawRect(0, 0, w - 1, h - 1)

        plot_x = self.GUTTER
        plot_w = max(60, w - plot_x - self.VALCOL - 14)

        vals = [v for _, v in self.terms]
        maxpos = max([v for v in vals if v > 0] or [0.0])
        maxneg = max([-v for v in vals if v < 0] or [0.0])
        span = (maxpos + maxneg) or 1.0
        px_per_mev = plot_w / span
        zero_x = plot_x + maxneg * px_per_mev

        top, bot = 20, h - 24

        q.setFont(QFont("Andale Mono", 8))
        step = self._nice_step(span)
        k = -int(maxneg // step) - 1
        while True:
            v = k * step
            k += 1
            if v > maxpos + step: break
            x = zero_x + v * px_per_mev
            if not (plot_x - 1 <= x <= plot_x + plot_w + 1): continue
            q.setPen(QPen(QColor(p["grid_minor"]), 1))
            q.drawLine(QPointF(x, top), QPointF(x, bot))
            q.setPen(QColor(p["faint"]))
            q.drawText(QRectF(x - 26, top - 15, 52, 13), Qt.AlignCenter, f"{v:+.0f}")

        q.setPen(QPen(QColor(p["frame"]), 1.6))
        q.drawLine(QPointF(zero_x, top), QPointF(zero_x, bot))

        fl = QFont("DIN Alternate", 10, QFont.Bold)
        fl.setLetterSpacing(QFont.AbsoluteSpacing, 1.0)
        fv = QFont("Andale Mono", 11, QFont.Bold)
        fm = QFontMetrics(fl)

        y = top + 2
        for caption, val in self.terms:
            col = QColor(p["green"] if val >= 0 else p["red"])
            length = val * px_per_mev
            bar = QRectF(min(zero_x, zero_x + length), y + 5, max(abs(length), 1.0), self.ROW - 13)

            q.setPen(Qt.NoPen)
            q.setBrush(QBrush(col.darker(155) if self.is_dark else col.lighter(150)))
            q.drawRect(bar)
            q.setPen(QPen(col, 1.2))
            q.setBrush(Qt.NoBrush)
            q.drawRect(bar)

            q.setFont(fl)
            q.setPen(QColor(p["label"]))
            avail = self.GUTTER - 18
            q.drawText(QRectF(8, y, avail, self.ROW), Qt.AlignVCenter | Qt.AlignRight,
                       fm.elidedText(caption, Qt.ElideMiddle, int(avail)))

            q.setFont(fv)
            q.setPen(col)
            q.drawText(QRectF(w - self.VALCOL - 8, y, self.VALCOL, self.ROW),
                       Qt.AlignVCenter | Qt.AlignRight, f"{val:+.2f}")
            y += self.ROW

        q.setPen(QPen(QColor(p["frame"]), 1))
        q.drawLine(8, h - 22, w - 8, h - 22)
        f = QFont("DIN Alternate", 10, QFont.Bold)
        f.setLetterSpacing(QFont.AbsoluteSpacing, 1.8)
        q.setFont(f)
        q.setPen(QColor(p["label"]))
        q.drawText(QRectF(8, h - 21, self.GUTTER - 18, 20), Qt.AlignVCenter | Qt.AlignRight, "TOPLAM Eb")
        q.setFont(QFont("Andale Mono", 12, QFont.Bold))
        q.setPen(QColor(p["accent"]))
        q.drawText(QRectF(self.GUTTER, h - 21, w - self.GUTTER - 8, 20), Qt.AlignVCenter | Qt.AlignLeft, f"{self.total:.3f} MeV")

class EngravedScale(QWidget):
    def __init__(self, vmin=0, vmax=100, divisions=10, unit="", parent=None, is_dark=True):
        super().__init__(parent)
        self.vmin, self.vmax, self.div, self.unit = vmin, vmax, divisions, unit
        self.is_dark = is_dark
        self.setFixedHeight(22)

    def apply_theme(self, is_dark: bool):
        self.is_dark = is_dark
        self.update()

    def paintEvent(self, _):
        p = palette(self.is_dark)
        q = QPainter(self)
        w = self.width()
        pad = 7
        q.setFont(QFont("Andale Mono", 8))
        for i in range(self.div + 1):
            frac = i / self.div
            x = pad + frac * (w - 2 * pad)
            major = (i % 2 == 0)
            q.setPen(QPen(QColor(p["frame"] if major else p["grid_minor"]), 1.2 if major else 0.8))
            q.drawLine(QPointF(x, 0), QPointF(x, 6 if major else 3))
            if major:
                v = self.vmin + (self.vmax - self.vmin) * frac
                box = QRectF(x - 22, 7, 44, 13)
                align = Qt.AlignCenter
                if box.left() < 0:
                    box = QRectF(0, 7, 44, 13); align = Qt.AlignLeft | Qt.AlignVCenter
                elif box.right() > w:
                    box = QRectF(w - 44, 7, 44, 13); align = Qt.AlignRight | Qt.AlignVCenter
                q.setPen(QColor(p["faint"]))
                q.drawText(box, align, f"{v:g}")

class ScanlineOverlay(QWidget):
    def __init__(self, parent=None, enabled=True):
        super().__init__(parent)
        self.setAttribute(Qt.WA_TransparentForMouseEvents, True)
        self.setAttribute(Qt.WA_NoSystemBackground, True)
        self.enabled = enabled
        if parent is not None:
            parent.installEventFilter(self)
            self.setGeometry(parent.rect())
            self.raise_()

    def eventFilter(self, obj, ev):
        if obj is self.parent() and ev.type() in (QEvent.Type.Resize, QEvent.Type.Show):
            self.setGeometry(self.parent().rect())
            self.raise_()
        return False

    def set_enabled(self, on: bool):
        self.enabled = on
        self.setVisible(on)
        self.update()

    def paintEvent(self, _):
        if not self.enabled: return
        q = QPainter(self)
        q.setPen(QPen(QColor(0, 0, 0, 26), 1))
        for y in range(0, self.height(), 3):
            q.drawLine(0, y, self.width(), y)
        g = QLinearGradient(0, 0, 0, self.height())
        g.setColorAt(0.0, QColor(0, 0, 0, 34))
        g.setColorAt(0.14, QColor(0, 0, 0, 0))
        g.setColorAt(0.86, QColor(0, 0, 0, 0))
        g.setColorAt(1.0, QColor(0, 0, 0, 34))
        q.fillRect(self.rect(), QBrush(g))

class ThemeTransitionOverlay(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAttribute(Qt.WA_TransparentForMouseEvents, True)
        self.setAttribute(Qt.WA_NoSystemBackground, True)
        self._pixmap = None
        self._opacity = 1.0
        self.hide()

        self.anim = QPropertyAnimation(self, b"overlayOpacity")
        self.anim.setDuration(460)
        self.anim.setEasingCurve(QEasingCurve.InOutQuart)
        self.anim.finished.connect(self._on_finished)

    def _on_finished(self):
        self._pixmap = None
        self.hide()

    def get_opacity(self): return self._opacity
    def set_opacity(self, val): self._opacity = val; self.update()
    overlayOpacity = Property(float, get_opacity, set_opacity)

    def start_fade(self, pixmap):
        if not pixmap or pixmap.isNull(): return
        self._pixmap = pixmap
        self.setGeometry(self.parent().rect())
        self.show()
        self.raise_()
        self.anim.stop()
        self.anim.setStartValue(1.0)
        self.anim.setEndValue(0.0)
        self.anim.start()

    def paintEvent(self, event):
        if self._pixmap and not self._pixmap.isNull():
            painter = QPainter(self)
            painter.setOpacity(self._opacity)
            painter.drawPixmap(self.rect(), self._pixmap)

class SplashScreen(QWidget):
    """Uygulama açılışında kısaca gösterilen 'güç açılışı' ekranı — logo,
    başlık ve gerçek kuruluş adımlarını izleyen kısa bir telemetri şeridi,
    ardından ana pencereye yumuşak geçiş (fade). Her zaman koyu temada
    çizilir (kullanıcının tema tercihinden bağımsız — klasik bir cihaz
    açılış ekranı hissi için)."""

    finished = Signal()

    STEPS = [
        "NUBASE2020 / AME2020 VERİTABANI YÜKLENİYOR...",
        "FİZİK MOTORU BAŞLATILIYOR...",
        "ARAYÜZ MODÜLLERİ İNŞA EDİLİYOR...",
        "SİSTEM HAZIR.",
    ]

    def __init__(self, pixmap=None):
        super().__init__(None, Qt.FramelessWindowHint | Qt.Tool)
        self.setFixedSize(480, 320)
        self._pixmap = pixmap
        self._step_index = 0

        screen = QApplication.primaryScreen()
        if screen is not None:
            geo = screen.geometry()
            self.move(geo.center().x() - self.width() // 2,
                      geo.center().y() - self.height() // 2)

        self._fade_anim = QPropertyAnimation(self, b"windowOpacity")
        self._fade_anim.setDuration(420)
        self._fade_anim.setEasingCurve(QEasingCurve.InOutQuart)
        self._fade_anim.finished.connect(self._on_fade_finished)

    def set_step(self, index: int):
        """Gerçek kuruluş ilerlemesine göre çağrılır (ana pencere __init__
        içindeki kontrol noktalarından) — sahte bir zamanlayıcı yerine
        gerçekten yapılan işi yansıtır. Hemen (senkron) yeniden çizer."""
        self._step_index = max(0, min(index, len(self.STEPS) - 1))
        self.repaint()

    def start_fade_out(self):
        self._fade_anim.stop()
        self._fade_anim.setStartValue(1.0)
        self._fade_anim.setEndValue(0.0)
        self._fade_anim.start()

    def _on_fade_finished(self):
        # ÖNEMLİ SIRALAMA: önce 'finished' yayınlanır (bu, dinleyicide ana
        # pencereyi gösterir), ANCAK SONRA kapatılır. Tersi sırayla
        # (önce close()) QApplication'ın varsayılan
        # quitOnLastWindowClosed=True davranışı, o anda görünür tek pencere
        # olan splash kapanır kapanmaz bir çıkış olayı kuyruğa alabilir —
        # ana pencere hemen ardından gösterilse bile uygulama sonradan
        # sessizce kapanabilir (gözlemlendi, düzeltildi).
        self.finished.emit()
        self.close()

    def paintEvent(self, _):
        p = PALETTE_DARK
        q = QPainter(self)
        q.setRenderHint(QPainter.Antialiasing)

        q.fillRect(self.rect(), QColor(p["chassis"]))
        q.setPen(QPen(QColor(p["frame"]), 1.2))
        q.drawRect(self.rect().adjusted(0, 0, -1, -1))
        q.setPen(QPen(QColor(p["accent_dim"]), 2.4))
        q.drawLine(0, 1, self.width(), 1)

        if self._pixmap and not self._pixmap.isNull():
            px = self._pixmap
            x = (self.width() - px.width()) // 2
            q.drawPixmap(x, 46, px)

        f_title = QFont("DIN Alternate", 16, QFont.Bold)
        f_title.setLetterSpacing(QFont.AbsoluteSpacing, 2.4)
        q.setFont(f_title)
        q.setPen(QColor(p["accent"]))
        q.drawText(QRectF(0, 188, self.width(), 28), Qt.AlignCenter, "SİNOP ÜNİVERSİTESİ")

        f_sub = QFont("Andale Mono", 9)
        q.setFont(f_sub)
        q.setPen(QColor(p["label"]))
        q.drawText(QRectF(20, 216, self.width() - 40, 18), Qt.AlignCenter,
                   "NÜKLEER ENERJİ MÜHENDİSLİĞİ ENSTRÜMANTASYON KONSOLU")

        bar_rect = QRectF(40, 260, self.width() - 80, 6)
        q.setPen(Qt.NoPen)
        q.setBrush(QColor(p["inset"]))
        q.drawRect(bar_rect)
        frac = self._step_index / max(1, len(self.STEPS) - 1)
        fill_rect = QRectF(bar_rect.x(), bar_rect.y(), bar_rect.width() * frac, bar_rect.height())
        q.setBrush(QColor(p["accent"]))
        q.drawRect(fill_rect)
        q.setPen(QPen(QColor(p["frame"]), 1))
        q.setBrush(Qt.NoBrush)
        q.drawRect(bar_rect)

        f_step = QFont("Andale Mono", 8)
        q.setFont(f_step)
        q.setPen(QColor(p["faint"]))
        q.drawText(QRectF(40, 271, self.width() - 80, 16), Qt.AlignLeft | Qt.AlignVCenter,
                   self.STEPS[self._step_index])

        f_model = QFont("Andale Mono", 8)
        q.setFont(f_model)
        q.setPen(QColor(p["faint"]))
        q.drawText(QRectF(0, 296, self.width(), 16), Qt.AlignCenter, "MODEL SNU-1956-ATOM-IV")

class SmoothLogoWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._blend = 0.0
        self.pix_dark = None
        self.pix_light = None
        self.setFixedHeight(120)
        self.setMinimumWidth(140)

        self.anim = QPropertyAnimation(self, b"blendFactor")
        self.anim.setDuration(460)
        self.anim.setEasingCurve(QEasingCurve.InOutQuart)

    def get_blend(self): return self._blend
    def set_blend(self, val): self._blend = val; self.update()
    blendFactor = Property(float, get_blend, set_blend)

    def set_pixmaps(self, pix_dark, pix_light):
        self.pix_dark = pix_dark
        self.pix_light = pix_light
        self.update()

    def transition_to(self, to_light: bool):
        self.anim.stop()
        self.anim.setStartValue(self._blend)
        self.anim.setEndValue(1.0 if to_light else 0.0)
        self.anim.start()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        painter.setRenderHint(QPainter.SmoothPixmapTransform)
        w, h = self.width(), self.height()
        if self.pix_dark and not self.pix_dark.isNull() and self._blend < 1.0:
            x = (w - self.pix_dark.width()) // 2
            y = (h - self.pix_dark.height()) // 2
            painter.setOpacity(1.0 - self._blend)
            painter.drawPixmap(x, y, self.pix_dark)
        if self.pix_light and not self.pix_light.isNull() and self._blend > 0.0:
            x = (w - self.pix_light.width()) // 2
            y = (h - self.pix_light.height()) // 2
            painter.setOpacity(self._blend)
            painter.drawPixmap(x, y, self.pix_light)

class MplCanvas(FigureCanvas):
    def __init__(self):
        self.fig = Figure(figsize=(6, 4), facecolor=PALETTE_DARK["chassis"])
        self.axes = self.fig.add_subplot(111)
        self.axes.set_facecolor(PALETTE_DARK["inset"])
        super().__init__(self.fig)
        self._is_dark = True

    def apply_theme(self, is_dark: bool):
        self._is_dark = is_dark
        p = palette(is_dark)
        self.fig.set_facecolor(p["chassis"])
        self.axes.set_facecolor(p["inset"])

    def trace(self, x, y, color, label, lw=1.7, bloom=True):
        if bloom and self._is_dark:
            self.axes.plot(x, y, color=color, lw=lw * 3.6, alpha=0.10, solid_capstyle="round", zorder=2)
        self.axes.plot(x, y, color=color, lw=lw, label=label, solid_capstyle="round", zorder=3)

    def apply_graticule(self, x_div=10, y_div=8, subdiv=5):
        p = palette(self._is_dark)
        ax = self.axes
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        ax.xaxis.set_major_locator(MultipleLocator((x1 - x0) / x_div))
        ax.yaxis.set_major_locator(MultipleLocator((y1 - y0) / y_div))
        ax.xaxis.set_minor_locator(MultipleLocator((x1 - x0) / (x_div * subdiv)))
        ax.yaxis.set_minor_locator(MultipleLocator((y1 - y0) / (y_div * subdiv)))

        if self._is_dark:
            ax.grid(True, which="major", color=p["grid_major"], lw=0.8, ls="-", alpha=1.0)
            ax.grid(True, which="minor", color=p["grid_minor"], lw=0.5, ls="-", alpha=1.0)
        else:
            ax.grid(True, which="minor", color=p["grid_minor"], lw=0.45, ls="-", alpha=0.85)
            ax.grid(True, which="major", color=p["grid_major"], lw=0.8, ls="-", alpha=0.85)

        ax.set_axisbelow(True)
        for s in ax.spines.values():
            s.set_color(p["frame"])
            s.set_linewidth(1.1)

        ax.tick_params(which="both", direction="in", top=True, right=True, colors=p["label"], labelsize=8.5)
        ax.tick_params(which="major", length=5, width=1.0)
        ax.tick_params(which="minor", length=2.5, width=0.6)
        for lbl in ax.get_xticklabels() + ax.get_yticklabels():
            lbl.set_fontfamily("monospace")

    def readout_box(self, text: str, loc="lower right"):
        p = palette(self._is_dark)
        xa, ha = (0.985, "right") if "right" in loc else (0.015, "left")
        ya, va = (0.03, "bottom") if "lower" in loc else (0.97, "top")
        self.axes.text(
            xa, ya, text, transform=self.axes.transAxes,
            ha=ha, va=va, fontfamily="monospace", fontsize=8.5,
            color=p["readout"],
            bbox=dict(boxstyle="square,pad=0.45",
                      facecolor=p["chassis"] if self._is_dark else p["panel"],
                      edgecolor=p["frame"], linewidth=0.9),
            zorder=6,
        )

    def style_legend(self):
        p = palette(self._is_dark)
        leg = self.axes.legend(
            loc="upper right", framealpha=1.0, fancybox=False, borderpad=0.7,
            facecolor=p["chassis"] if self._is_dark else p["panel"],
            edgecolor=p["frame"], labelcolor=p["text"],
            prop={"family": "monospace", "size": 8.5},
        )
        leg.get_frame().set_linewidth(0.9)
        return leg


class SearchableIsotopeCombo(PanelComboBox):
    """PanelComboBox + yazarak filtreleme (contains-eşleşme, prefix değil).
    3000+ nüklitlik listelerde (Reaksiyon hedef/mermi, Bağlanma izotop seçimi)
    belirli bir izotopu hızlıca bulmak için. userData tabanlı seçim aynen
    PanelComboBox gibi çalışmaya devam eder; sadece yazarken filtre gevşer."""

    def __init__(self, parent=None, is_dark=True, placeholder="İzotop ara… (örn. Fe, 56Fe, Demir)"):
        super().__init__(parent, is_dark=is_dark)
        self.setEditable(True)
        self.setInsertPolicy(QComboBox.NoInsert)
        # setEditable(True) bir varsayılan QCompleter kurar (kendi modeliyle,
        # userData/seçim bağlantısı korunur); sadece eşleşme modunu gevşetiyoruz.
        comp = self.completer()
        if comp is not None:
            comp.setCaseSensitivity(Qt.CaseInsensitive)
            comp.setFilterMode(Qt.MatchContains)
            comp.setCompletionMode(QCompleter.PopupCompletion)
        le = self.lineEdit()
        if le is not None:
            le.setPlaceholderText(placeholder)


class SidebarNav(QWidget):
    """Kategorilere ayrılmış dikey modül menüsü (retro enstrümantasyon paneli
    kanal listesi görünümünde). Modül sayısı arttıkça yatay sekme şeridi yerine
    bu liste büyür/kaydırılır — üst sınır yok."""

    itemSelected = Signal(int)

    def __init__(self, parent=None, is_dark=True, width=228):
        super().__init__(parent)
        self.is_dark = is_dark
        self.setFixedWidth(width)
        self._buttons = []
        self._group = QButtonGroup(self)
        self._group.setExclusive(True)

        self._scroll = QScrollArea(self)
        self._scroll.setWidgetResizable(True)
        self._scroll.setFrameShape(QFrame.NoFrame)
        self._scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self._inner = QWidget()
        self._layout = QVBoxLayout(self._inner)
        self._layout.setContentsMargins(0, 10, 0, 16)
        self._layout.setSpacing(1)
        self._layout.addStretch(1)   # her zaman en altta kalır; öğeler öncesine eklenir
        self._insert_pos = 0

        self._scroll.setWidget(self._inner)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(self._scroll)

        self._header_lbl = None
        self._footer_lbl = None

        self.setObjectName("SidebarNav")
        self.apply_theme(is_dark)

    def set_header(self, title: str, subtitle: str = ""):
        """En üste, panel kimliğini vurgulayan küçük bir 'künye plakası' ekler
        (yalnız bir kere çağrılmalı — sections/items eklenmeden ÖNCE)."""
        lbl = QLabel()
        lbl.setTextFormat(Qt.RichText)
        lbl.setWordWrap(True)
        self._header_lbl = lbl
        self._header_title = title
        self._header_subtitle = subtitle
        self._layout.insertWidget(self._insert_pos, lbl)
        self._insert_pos += 1
        self._style_header()

    def set_footer(self, text: str):
        """En alta (kayan listenin dibine, her zaman görünür) küçük bir
        telemetri/özet satırı ekler."""
        lbl = QLabel()
        lbl.setTextFormat(Qt.RichText)
        lbl.setWordWrap(True)
        self._footer_lbl = lbl
        self._footer_text = text
        self._layout.addWidget(lbl)   # stretch'ten SONRA -> daima dipte kalır
        self._style_footer()

    def set_footer_text(self, text: str):
        self._footer_text = text
        if self._footer_lbl is not None:
            self._style_footer()

    def _style_header(self):
        if self._header_lbl is None:
            return
        p = palette(self.is_dark)
        self._header_lbl.setText(
            f"<div style='padding:14px 16px 10px 16px;'>"
            f"<span style='color:{p['accent']}; font-family:{FONT_DIN}; font-weight:bold;"
            f" font-size:13px; letter-spacing:1.6px;'>{self._header_title}</span><br>"
            f"<span style='color:{p['faint']}; font-family:{FONT_DATA}; font-size:9px;"
            f" letter-spacing:0.8px;'>{self._header_subtitle}</span></div>"
        )
        self._header_lbl.setStyleSheet(
            f"QLabel {{ background: transparent; border-bottom: 1px solid {p['seam']}; }}"
        )

    def _style_footer(self):
        if self._footer_lbl is None:
            return
        p = palette(self.is_dark)
        self._footer_lbl.setText(
            f"<div style='padding:9px 16px; color:{p['faint']}; font-family:{FONT_DATA};"
            f" font-size:9px; letter-spacing:0.8px;'>{self._footer_text}</div>"
        )
        self._footer_lbl.setStyleSheet(
            f"QLabel {{ background: transparent; border-top: 1px solid {p['seam']}; }}"
        )

    def add_section(self, title: str):
        lbl = QLabel(title)
        lbl.setObjectName("SidebarSection")
        self._layout.insertWidget(self._insert_pos, lbl)
        self._insert_pos += 1
        self._style_section_label(lbl)

    def add_item(self, label: str) -> int:
        btn = QPushButton(label)
        btn.setCheckable(True)
        btn.setCursor(Qt.PointingHandCursor)
        idx = len(self._buttons)
        btn.clicked.connect(lambda _checked=False, i=idx: self.itemSelected.emit(i))
        self._group.addButton(btn)
        self._buttons.append(btn)
        self._layout.insertWidget(self._insert_pos, btn)
        self._insert_pos += 1
        self._style_button(btn)
        if idx == 0:
            btn.setChecked(True)
        return idx

    def set_current(self, index: int):
        if 0 <= index < len(self._buttons):
            self._buttons[index].setChecked(True)

    def apply_theme(self, is_dark: bool):
        self.is_dark = is_dark
        p = palette(is_dark)
        self._scroll.setStyleSheet(
            f"QScrollArea {{ background: {p['tab_bg']}; border: none;"
            f" border-right: 1px solid {p['seam']}; }}"
            f"QScrollArea > QWidget > QWidget {{ background: {p['tab_bg']}; }}"
        )
        for i in range(self._layout.count()):
            item = self._layout.itemAt(i)
            w = item.widget() if item else None
            if w is self._header_lbl:
                continue
            elif w is self._footer_lbl:
                continue
            elif isinstance(w, QLabel):
                self._style_section_label(w)
            elif isinstance(w, QPushButton):
                self._style_button(w)
        self._style_header()
        self._style_footer()

    def _style_section_label(self, lbl: QLabel):
        p = palette(self.is_dark)
        lbl.setStyleSheet(
            f"QLabel {{ color: {p['faint']}; font-family: {FONT_DIN}; font-weight: bold;"
            f" font-size: 10px; letter-spacing: 1.6px;"
            f" padding: 16px 16px 5px 16px; background: transparent; }}"
        )

    def _style_button(self, btn: QPushButton):
        p = palette(self.is_dark)
        btn.setStyleSheet(f"""
            QPushButton {{
                text-align: left; border: none; border-left: 3px solid transparent;
                background: transparent; color: {p['label']};
                padding: 9px 14px 9px 13px; margin: 0 6px 0 0;
                font-family: {FONT_PANEL}; font-size: 12px;
            }}
            QPushButton:hover {{ background: {p['panel_alt']}; color: {p['text']}; }}
            QPushButton:checked {{
                background: {p['panel']}; color: {p['accent']};
                border-left: 3px solid {p['accent']}; font-weight: bold;
            }}
        """)
