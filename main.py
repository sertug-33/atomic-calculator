import sys
import os
import re
import io
import math
import numpy as np
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QStackedWidget, QVBoxLayout,
    QHBoxLayout, QLabel, QLineEdit, QPushButton, QTextEdit,
    QTableWidget, QTableWidgetItem, QHeaderView, QSlider, QFrame, QGridLayout,
    QSizePolicy, QFileDialog, QScrollArea
)
from PySide6.QtCore import Qt, QTimer, QDateTime, QMarginsF
from PySide6.QtGui import QPixmap, QIcon, QTextDocument, QFont, QPageLayout
from PySide6.QtPrintSupport import QPrinter
from matplotlib.ticker import LogLocator, NullFormatter

from database import ISOTOPES
from physics_engine import (
    semf_binding_energy, semf_breakdown, separation_energies,
    calculate_reaction, calculate_bateman, solve_channel,
    calculate_moderation, MODERATORS,
    binding_energy, semf_deviation,
    decay_channels, specific_activity,
    decay_series_chain, DECAY_SERIES,
    gamma_attenuation, gamma_thickness_for_factor, SHIELDING_MATERIALS,
    gamma_attenuation_multiline, RADIOACTIVE_SOURCES, buildup_factor,
    dose_rate_point_source, dose_rate_multiline, air_mu_en_rho,
    eta_factor, thermal_utilization_factor, resonance_escape_probability, k_infinity
)
from isotope_explain import explain_isotope
from density import calculate_material, MATERIALS, parse_formula
from ui_components import (
    MplCanvas, SmoothLogoWidget, ThemeTransitionOverlay, SplashScreen,
    build_qss, palette, panel_qss,
    IndicatorLamp, AnalogGauge, SemfBarMeter, EngravedScale,
    ScanlineOverlay, PanelDivider, PanelComboBox,
    SearchableIsotopeCombo, SidebarNav,
    DARK_THEME_QSS, LIGHT_THEME_QSS
)

ELEMENT_NAMES = {
    "H": "Hidrojen", "He": "Helyum", "Li": "Lityum", "Be": "Berilyum", "B": "Bor",
    "C": "Karbon", "N": "Azot", "O": "Oksijen", "F": "Flor", "Ne": "Neon",
    "Na": "Sodyum", "Mg": "Magnezyum", "Al": "Alüminyum", "Si": "Silisyum",
    "P": "Fosfor", "S": "Kükürt", "Cl": "Klor", "Ar": "Argon", "K": "Potasyum",
    "Ca": "Kalsiyum", "Sc": "Skandiyum", "Ti": "Titanyum", "V": "Vanadyum",
    "Cr": "Krom", "Mn": "Mangan", "Fe": "Demir", "Co": "Kobalt", "Ni": "Nikel",
    "Cu": "Bakır", "Zn": "Çinko", "Kr": "Kripton", "Rb": "Rubidyum", "Sr": "Stronsiyum",
    "Y": "İtriyum", "Zr": "Zirkonyum", "Mo": "Molibden", "Tc": "Teknesyum",
    "Ag": "Gümüş", "Cd": "Kadmiyum", "I": "İyot", "Xe": "Ksenon", "Cs": "Sezyum",
    "Ba": "Baryum", "Nd": "Neodim", "Sm": "Samaryum", "Gd": "Gadolinyum",
    "Au": "Altın", "Pb": "Kurşun", "Bi": "Bismut", "Po": "Polonyum", "Rn": "Radon",
    "Ra": "Radyum", "Th": "Toryum", "U": "Uranyum", "Np": "Neptünyum",
    "Pu": "Plütonyum", "Am": "Amerikyum", "Cf": "Kaliforniyum"
}

# test_physics.py'deki toplam doğrulama testi sayısı (sidebar footer'ında
# gösterilir) — yeni bir test_* fonksiyonu eklendiğinde burası da güncellenmeli.
APP_TEST_COUNT = 20

DOSE_ACTIVITY_UNITS = {
    "Bq": 1.0, "kBq": 1.0e3, "MBq": 1.0e6, "GBq": 1.0e9,
    "µCi": 3.7e4, "mCi": 3.7e7, "Ci": 3.7e10,
}
DOSE_DISTANCE_UNITS = {"cm": 1.0, "m": 100.0}

PRESET_REACTIONS = [
    ("—  KATALOG REFERANSI SEÇİNİZ",                              ""),
    ("REF-01   ENDO   Rutherford yapay transmutasyon · 1919",     "14N + alpha -> 17O + p"),
    ("REF-02   EKZO   Chadwick nötron keşfi · 1932",              "9Be + alpha -> 12C + n"),
    ("REF-03   ENDO   PET F-18 medikal siklotron sentezi",        "18O + p -> 18F + n"),
    ("REF-04   ENDO   Monokromatik nötron kaynağı · Li-p",        "7Li + p -> 7Be + n"),
    ("REF-05   EKZO   D-T termonükleer füzyon odası",             "d + t -> alpha + n"),
    ("REF-06   EKZO   U-235 termal nötron fisyon kanalı",         "235U + n -> 92Kr + 141Ba + 3n"),
    ("REF-07   EKZO   BNCT bor nötron yakalama terapisi",         "10B + n -> 7Li + alpha"),
    ("REF-08   EKZO   Am-241 iyonizasyon dedektörü · alfa",       "241Am -> 237Np + alpha"),
]

LIGHT_PARTICLES = [
    ("p  (Proton)", "p"),
    ("n  (Nötron)", "n"),
    ("alpha  (Alfa · 4He)", "alpha"),
    ("d  (Döteron · 2H)", "d"),
    ("t  (Triton · 3H)", "t")
]

W_FIS = 58

def _fis(label, value, unit=""):
    sag = f"{value}  {unit}".rstrip()
    nokta = "." * max(2, W_FIS - len(label) - len(sag) - 2)
    return f" {label} {nokta} {sag}\n"

def get_iso_display_name(iso_key):
    m = re.match(r"^(\d+)([a-zA-Z]+)$", iso_key)
    if m:
        a_str, sym = m.groups()
        el_name = ELEMENT_NAMES.get(sym, sym)
        return f"{iso_key} — {el_name}-{a_str}"
    return iso_key

_STABLE_EBA_CURVE = None

def get_stable_eba_curve():
    """Kararlı nüklitlerin (A, Eb/A) noktalarını bir kez hesaplayıp önbelleğe alır.
    Bağlanma enerjisi eğrisi grafiği için kullanılır (klasik 'Eb/A vs A' eğrisi)."""
    global _STABLE_EBA_CURVE
    if _STABLE_EBA_CURVE is None:
        pts = []
        for k, v in ISOTOPES.items():
            if k in ("n", "p"):
                continue
            if v.get("half_life") != "Kararlı":
                continue
            try:
                eb = binding_energy(k)
            except Exception:
                continue
            if eb["Eb_per_A"] and eb["Eb_per_A"] > 0:
                pts.append((v["A"], eb["Eb_per_A"]))
        _STABLE_EBA_CURVE = sorted(pts)
    return _STABLE_EBA_CURVE

_SEGRE_CHART_DATA = None

def _decay_category(mode: str, is_stable: bool) -> str:
    """NUBASE2020 bozunma modu etiketini kaba bir kategoriye indirger
    (Segre / N-Z haritasında renklendirme için)."""
    if is_stable:
        return "stable"
    primary = mode.split("/")[0] if mode else "-"
    if primary in ("β-", "2β-"):
        return "beta_minus"
    if primary in ("β+", "EC"):
        return "beta_plus"
    if primary == "α":
        return "alpha"
    return "other"

def get_segre_chart_data():
    """Tüm nüklitlerin (N, Z) koordinatlarını bozunma kategorisine göre
    gruplandırıp bir kez önbelleğe alır (~3558 nokta, N-Z 'kararlılık vadisi' haritası)."""
    global _SEGRE_CHART_DATA
    if _SEGRE_CHART_DATA is None:
        groups = {"stable": ([], []), "beta_minus": ([], []), "beta_plus": ([], []),
                  "alpha": ([], []), "other": ([], [])}
        for k, v in ISOTOPES.items():
            if k in ("n", "p"):
                continue
            Z, A = v["Z"], v["A"]
            N = A - Z
            is_stable = v.get("half_life") == "Kararlı"
            cat = _decay_category(v.get("mode", "-"), is_stable)
            xs, ys = groups[cat]
            xs.append(N); ys.append(Z)
        _SEGRE_CHART_DATA = groups
    return _SEGRE_CHART_DATA

_NZ_LOOKUP = None

def get_nz_lookup():
    """(N, Z) tam sayı çiftinden izotop anahtarına eşleme — Nüklit Haritası'ndaki
    tıklama olaylarında en yakın noktayı bulmak için bir kez inşa edilip önbelleğe
    alınır."""
    global _NZ_LOOKUP
    if _NZ_LOOKUP is None:
        d = {}
        for k, v in ISOTOPES.items():
            if k in ("n", "p"):
                continue
            d[(v["A"] - v["Z"], v["Z"])] = k
        _NZ_LOOKUP = d
    return _NZ_LOOKUP

def load_smart_pixmap(image_path: str, target_height: int = 115) -> QPixmap:
    if not os.path.exists(image_path):
        return None
    try:
        from PIL import Image
        with Image.open(image_path) as im:
            bbox = im.getbbox()
            if bbox:
                im = im.crop(bbox)
            buf = io.BytesIO()
            im.save(buf, format="PNG")
            pix = QPixmap()
            pix.loadFromData(buf.getvalue())
            return pix.scaledToHeight(target_height, Qt.SmoothTransformation)
    except Exception:
        pix = QPixmap(image_path)
        if not pix.isNull():
            return pix.scaledToHeight(target_height, Qt.SmoothTransformation)
        return None

class NuclearCalculatorApp(QMainWindow):
    def __init__(self, splash: SplashScreen | None = None):
        super().__init__()
        self.setWindowTitle("Sinop Üniversitesi • Nükleer Enerji Mühendisliği Enstrümantasyon Konsolu [SNU-1956-ATOM-IV]")
        self.resize(1360, 860)
        self.setMinimumSize(1300, 760)

        self.is_dark = True
        self.setStyleSheet(build_qss(self.is_dark))

        base_dir = os.path.dirname(os.path.abspath(__file__))
        logo_light_path = os.path.join(base_dir, "logo.png")
        logo_dark_path = os.path.join(base_dir, "lacivertlogo.png")

        self.cached_pix_light = load_smart_pixmap(logo_light_path, 115)
        self.cached_pix_dark = load_smart_pixmap(logo_dark_path, 115) or self.cached_pix_light

        if self.cached_pix_light:
            self.setWindowIcon(QIcon(self.cached_pix_light))

        if splash is not None:
            splash.set_step(1)   # veritabanı zaten import anında yüklendi -> "motor başlatılıyor"
            QApplication.processEvents()

        # --- Sol kenar modül menüsü + sağda tek sayfalık içerik yığını -----
        # (Eskiden QTabWidget ile yatay sekmelerdi; modül sayısı arttıkça
        # sıkışmaması için kategorilere ayrılmış dikey listeye geçildi.)
        self.stack = QStackedWidget()
        self.sidebar = SidebarNav(is_dark=self.is_dark)
        self.sidebar.itemSelected.connect(self._on_nav_selected)
        self.sidebar.set_header("⚛ ATOMİK HESAPLAYICI", "MODEL SNU-1956-ATOM-IV")

        if splash is not None:
            splash.set_step(2)
            QApplication.processEvents()

        MODULES = [
            ("BÖLÜM I · TEMEL",         [("Künye",          self.create_about_tab),
                                          ("Nüklit Fişleri",  self.create_db_tab),
                                          ("Nüklit Haritası", self.create_chart_tab)]),
            ("BÖLÜM II · REAKSİYON",   [("Reaksiyon",       self.create_q_tab),
                                          ("Bağlanma",        self.create_binding_tab)]),
            ("BÖLÜM III · BOZUNMA",    [("Bozunma",         self.create_decay_tab),
                                          ("Bozunma Modları", self.create_decay_modes_tab),
                                          ("Bozunma Serisi",  self.create_decay_series_tab)]),
            ("BÖLÜM IV · NÖTRON FİZİĞİ",[("Moderasyon",     self.create_moderation_tab),
                                          ("Kritiklik",      self.create_criticality_tab)]),
            ("BÖLÜM V · MALZEME",      [("Yoğunluk",        self.create_density_tab)]),
            ("BÖLÜM VI · RADYASYON",   [("Zırhlama",        self.create_shielding_tab),
                                          ("Doz Hızı",        self.create_dose_tab)]),
        ]
        self._breadcrumbs = []
        for section_title, items in MODULES:
            self.sidebar.add_section(section_title)
            for label, factory in items:
                self.sidebar.add_item(label)
                self.stack.addWidget(factory())
                self._breadcrumbs.append((section_title, label))
        self.sidebar.set_footer(f"{len(self._breadcrumbs)} MODÜL   ·   {APP_TEST_COUNT} DOĞRULAMA TESTİ   ·   ÇEVRİMDIŞI")
        self.update_about_display()   # Künye ilk oluşturulduğunda modül sayısı henüz tam değildi

        self.theme_btn = QPushButton()
        self.theme_btn.setObjectName("ModeSwitch")
        self.theme_btn.setCursor(Qt.PointingHandCursor)
        self.theme_btn.clicked.connect(self.toggle_theme)
        self._sync_theme_btn()

        self.breadcrumb_lbl = QLabel()
        self.breadcrumb_lbl.setObjectName("Breadcrumb")

        self.topbar = QWidget()
        self.topbar.setObjectName("TopBar")
        topbar_l = QHBoxLayout(self.topbar)
        topbar_l.setContentsMargins(14, 6, 10, 6)
        topbar_l.addWidget(self.breadcrumb_lbl)
        topbar_l.addStretch(1)
        topbar_l.addWidget(self.theme_btn)
        self._sync_topbar()
        self._sync_breadcrumb(0)

        body = QWidget()
        body_l = QHBoxLayout(body)
        body_l.setContentsMargins(0, 0, 0, 0)
        body_l.setSpacing(0)
        body_l.addWidget(self.sidebar)
        body_l.addWidget(self.stack, 1)

        central = QWidget()
        outer_l = QVBoxLayout(central)
        outer_l.setContentsMargins(0, 0, 0, 0)
        outer_l.setSpacing(0)
        outer_l.addWidget(self.topbar)
        outer_l.addWidget(body, 1)

        self.stack.setCurrentIndex(0)
        self.setCentralWidget(central)

        self.info_status = QLabel()
        self.statusBar().addWidget(self.info_status)

        self.clock_lbl = QLabel()
        self.clock_lbl.setStyleSheet("font-family: 'Andale Mono', 'Menlo', monospace; font-weight: bold; padding-right: 12px;")
        self.statusBar().addPermanentWidget(self.clock_lbl)

        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_clock)
        self.timer.start(1000)
        self.update_clock()

        self.transition_overlay = ThemeTransitionOverlay(self)

        if splash is not None:
            splash.set_step(3)
            QApplication.processEvents()

    def _sync_theme_btn(self):
        p = palette(self.is_dark)
        self.theme_btn.setText("PANEL  ▮ KATOT" if self.is_dark else "PANEL  ▮ BÜLTEN")
        self.theme_btn.setStyleSheet(
            f"QPushButton#ModeSwitch {{"
            f" background:{p['chassis']}; color:{p['label']};"
            f" border:1px solid {p['seam']}; border-radius:0px;"
            f" padding:5px 12px; margin:0 8px 2px 0;"
            f" font-family:'DIN Alternate', 'Helvetica Neue', sans-serif; font-size:10px;"
            f" font-weight:bold; letter-spacing:1.6px; }}"
            f"QPushButton#ModeSwitch:hover {{ color:{p['accent']};"
            f" border-color:{p['frame']}; }}"
        )

    def _sync_topbar(self):
        p = palette(self.is_dark)
        self.topbar.setStyleSheet(
            f"QWidget#TopBar {{ background: {p['chassis']}; border-bottom: 1px solid {p['seam']}; }}"
        )
        self.breadcrumb_lbl.setStyleSheet(
            f"QLabel#Breadcrumb {{ color: {p['label']}; font-family: 'DIN Alternate', 'Helvetica Neue', sans-serif;"
            f" font-size: 11px; font-weight: bold; letter-spacing: 1.4px; }}"
        )

    def _sync_breadcrumb(self, index: int):
        if 0 <= index < len(self._breadcrumbs):
            section, label = self._breadcrumbs[index]
            p = palette(self.is_dark)
            self.breadcrumb_lbl.setText(
                f"{section}   <span style='color:{p['accent']};'>›</span>   {label}"
            )
            self.breadcrumb_lbl.setTextFormat(Qt.RichText)

    def _on_nav_selected(self, index: int):
        self.stack.setCurrentIndex(index)
        self._sync_breadcrumb(index)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if hasattr(self, "transition_overlay"):
            self.transition_overlay.setGeometry(self.rect())

    def update_clock(self):
        now = QDateTime.currentDateTime()
        self.clock_lbl.setText("TELEMETRİ  " + now.toString("dd.MM.yyyy   HH:mm:ss"))
        nabiz = "●" if now.time().second() % 2 == 0 else "○"
        self.info_status.setText(
            f"{nabiz}  ÇALIŞMA MODU AKTİF   |   SİNOP ÜNİVERSİTESİ NÜKLEER ENERJİ MÜHENDİSLİĞİ"
        )

    def toggle_theme(self):
        old_screen = self.grab()
        self.is_dark = not self.is_dark
        self.setStyleSheet(build_qss(self.is_dark))
        self._sync_theme_btn()
        self._sync_topbar()
        self._sync_breadcrumb(self.stack.currentIndex())

        widgets = [
            self.sidebar,
            self.lamp_z, self.lamp_a, self.lamp_thr, self.semf_meter,
            self.dev_gauge, self.decay_scale, self.preset_combo,
            self.iso_combo, self.divider_left,
            self.combo_target, self.combo_proj, self.combo_emit,
            self.mod_combo, self.divider_mod, self.db_lamp,
            self.dens_combo, self.dens_frame, self.chart_combo, self.chart_frame,
            self.dm_combo, self.dm_card_left, self.dm_card_right,
            self.ds_combo, self.ds_frame_left, self.ds_frame_right,
            self.shield_combo, self.shield_frame, self.shield_source_combo,
            self.dose_source_combo, self.dose_frame, self.dose_activity_unit,
            self.dose_distance_unit, self.dose_shield_combo,
            self.crit_frame, self.crit_mod_combo,
            self.shield_lamp, self.dose_lamp, self.crit_lamp,
            self.nuclide_divider1, self.nuclide_divider2,
        ]
        for wgt in widgets:
            if hasattr(wgt, "apply_theme"):
                wgt.apply_theme(self.is_dark)
            
        self.crt_overlay.set_enabled(self.is_dark)

        self.smooth_logo.transition_to(to_light=not self.is_dark)
        self.update_about_display()
        self.update_binding_display()
        self.update_decay_plot()
        self._sync_channel_ui()
        self.update_moderation_display()
        self.dens_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        self.chan_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        self.search_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        self.chart_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        self.chart_hint_lbl.setStyleSheet(f"color: {palette(self.is_dark)['faint']}; font-family: monospace; font-size: 9.5px;")
        self.nuclide_panel_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        self.handle_q_calc()
        self.handle_density_calc()
        self._draw_segre_chart()
        iso_for_panel = self.chart_combo.currentData() if hasattr(self, "chart_combo") else None
        self._update_nuclide_info_panel(iso_for_panel or "56Fe")
        self.update_decay_modes_display()
        self.ds_frame_left.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        self.ds_frame_right.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        self.update_decay_series_display()
        self.shield_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        self.handle_shielding_calc()
        self.dose_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        self.handle_dose_calc()
        self.crit_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        self.handle_criticality_calc()

        self.transition_overlay.start_fade(old_screen)

    def _export_report_pdf(self, get_text_fn, default_name: str):
        """Verilen metni (ölçüm fişi) PDF olarak diske kaydeder — QTextDocument +
        QPrinter kullanır, harici bir bağımlılık gerekmez. get_text_fn: () -> str."""
        text = get_text_fn()
        if not text or not text.strip():
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Raporu PDF Olarak Kaydet", default_name, "PDF Dosyası (*.pdf)"
        )
        if not path:
            return
        if not path.lower().endswith(".pdf"):
            path += ".pdf"

        doc = QTextDocument()
        doc.setDefaultFont(QFont("Andale Mono", 9))
        doc.setPlainText(text)

        printer = QPrinter(QPrinter.HighResolution)
        printer.setOutputFormat(QPrinter.PdfFormat)
        printer.setOutputFileName(path)
        printer.setPageMargins(QMarginsF(14, 14, 14, 14), QPageLayout.Millimeter)
        doc.print_(printer)

        self.info_status.setText(
            f"○  RAPOR KAYDEDİLDİ: {os.path.basename(path)}   |   SİNOP ÜNİVERSİTESİ NÜKLEER ENERJİ MÜHENDİSLİĞİ"
        )

    def _make_export_button(self, on_click):
        btn = QPushButton("RAPORU PDF OLARAK KAYDET")
        btn.setObjectName("PrimaryLever")
        btn.clicked.connect(on_click)
        return btn

    def create_about_tab(self):
        w = QWidget()
        outer_layout = QVBoxLayout(w)
        outer_layout.setContentsMargins(0, 0, 0, 0)
        outer_layout.setSpacing(0)

        # İçerik (lambalar/istatistikler eklendikçe) küçük pencerelerde panel
        # yüksekliğini aşabilir; QVBoxLayout bu durumda alt widget'ları
        # sıkıştırıp üst üste bindirebildiğinden (gözlemlendi, düzeltildi),
        # kart bir QScrollArea içine alınır — taşarsa kırpma/örtüşme yerine
        # kaydırma olur.
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setStyleSheet("QScrollArea { background: transparent; border: none; }"
                              "QScrollArea > QWidget > QWidget { background: transparent; }")
        self.about_scroll = scroll

        content = QWidget()
        main_layout = QVBoxLayout(content)
        main_layout.setAlignment(Qt.AlignHCenter | Qt.AlignTop)
        main_layout.setContentsMargins(20, 20, 20, 20)
        main_layout.setSpacing(12)

        self.about_card = QFrame()
        self.about_card.setMaximumWidth(740)
        card_layout = QVBoxLayout(self.about_card)
        card_layout.setSpacing(12)
        card_layout.setContentsMargins(16, 16, 16, 16)

        self.smooth_logo = SmoothLogoWidget()
        self.smooth_logo.set_pixmaps(self.cached_pix_dark, self.cached_pix_light)
        card_layout.addWidget(self.smooth_logo, 0, Qt.AlignCenter)

        self.uni_title_lbl = QLabel()
        self.uni_title_lbl.setAlignment(Qt.AlignCenter)
        card_layout.addWidget(self.uni_title_lbl)

        self.proj_title_lbl = QLabel()
        self.proj_title_lbl.setAlignment(Qt.AlignCenter)
        card_layout.addWidget(self.proj_title_lbl)

        self.about_details_lbl = QLabel()
        self.about_details_lbl.setTextFormat(Qt.RichText)
        self.about_details_lbl.setWordWrap(True)
        card_layout.addWidget(self.about_details_lbl)

        self.about_divider = PanelDivider(is_dark=self.is_dark)
        card_layout.addWidget(self.about_divider)

        status_cap = QLabel("SİSTEM DURUMU"); status_cap.setObjectName("Engraved")
        card_layout.addWidget(status_cap)

        lamps_row1 = QHBoxLayout(); lamps_row1.setSpacing(22)
        lamps_row2 = QHBoxLayout(); lamps_row2.setSpacing(22)
        self.about_lamp_db     = IndicatorLamp("NUBASE2020 VERİTABANI", is_dark=self.is_dark)
        self.about_lamp_engine = IndicatorLamp("FİZİK MOTORU DOĞRULANDI", is_dark=self.is_dark)
        self.about_lamp_pdf    = IndicatorLamp("PDF DIŞA AKTARMA HAZIR", is_dark=self.is_dark)
        self.about_lamp_offline = IndicatorLamp("ÇEVRİMDIŞI ÇALIŞMA", is_dark=self.is_dark)
        for l in (self.about_lamp_db, self.about_lamp_engine):
            l.set_state("ok")
            lamps_row1.addWidget(l)
        lamps_row1.addStretch(1)
        for l in (self.about_lamp_pdf, self.about_lamp_offline):
            l.set_state("ok")
            lamps_row2.addWidget(l)
        lamps_row2.addStretch(1)
        card_layout.addLayout(lamps_row1)
        card_layout.addLayout(lamps_row2)

        self.about_stats_lbl = QLabel()
        self.about_stats_lbl.setTextFormat(Qt.RichText)
        self.about_stats_lbl.setAlignment(Qt.AlignCenter)
        card_layout.addWidget(self.about_stats_lbl)

        main_layout.addWidget(self.about_card)
        scroll.setWidget(content)
        outer_layout.addWidget(scroll)
        self.update_about_display()
        return w

    def update_about_display(self):
        p = palette(self.is_dark)
        self.about_card.setStyleSheet(panel_qss(self.is_dark))
        self.uni_title_lbl.setText(f"<h2 style='color: {p['accent']}; margin: 0; text-align: center; letter-spacing: 2.2px; font-family: sans-serif;'>SİNOP ÜNİVERSİTESİ</h2>")
        self.proj_title_lbl.setText(f"<p style='color: {p['label']}; margin: 2px 0 0 0; text-align: center; font-family: monospace; font-size: 11px; letter-spacing: 1.2px;'>[ NÜKLEER ENERJİ MÜHENDİSLİĞİ ENSTRÜMANTASYON KONSOLU • MODEL: SNU-1956-ATOM-IV ]</p>")

        details_html = f"""
        <table width='100%' cellpadding='8' cellspacing='0'
               style='background-color:{p['inset']}; color:{p['text']};
                      border: 1px solid {p['seam']}; font-family:monospace; font-size:12px; line-height: 1.8;'>
            <tr>
                <td width='35%' style='color:{p['label']};'>ÖĞRENCİ</td>
                <td style='color:{p['accent']}; font-weight: bold;'>Ali Sertuğ Çetin</td>
            </tr>
            <tr>
                <td style='color:{p['label']};'>VERİTABANI MİMARİSİ</td>
                <td>AME2020 / NUBASE Çift Hassasiyetli Nükleer Veri Kütüphanesi</td>
            </tr>
        </table>
        
        <div style='color:{p['label']}; font-size: 11px; margin-top: 10px; line-height: 1.7; font-family: sans-serif;'>
            <b>TEORİK HESAPLAMA METODOLOJİSİ (ATOMIC AGE SERİSİ):</b><br>
            • <b>Kinematik Q-Değeri & Eşik Enerjisi (Eth):</b> Kütle merkezi momentum korunumu.<br>
            • <b>SEMF & Ayrıştırma Enerjileri (Sn, Sp):</b> Bethe-Weizsäcker sıvı damlası modeli.<br>
            • <b>Moderasyon & Kritiklik:</b> Letarji/ξ analitiği, dört-faktör formülü (k∞=η·f·p·ε).<br>
            • <b>Zırhlama, Doz Hızı & Bozunma Serisi:</b> NIST/IAEA referanslı gama fiziği.
        </div>
        """
        self.about_details_lbl.setText(details_html)

        n_modules = len(getattr(self, "_breadcrumbs", None) or []) or 13
        n_nuclides = len(ISOTOPES)
        stats_html = (
            f"<div style='margin-top:6px; font-family:monospace; font-size:11px;"
            f" letter-spacing:0.6px; color:{p['accent']};'>"
            f"{n_modules} MODÜL&nbsp;&nbsp;·&nbsp;&nbsp;{APP_TEST_COUNT} DOĞRULAMA TESTİ&nbsp;&nbsp;·&nbsp;&nbsp;"
            f"{n_nuclides} NÜKLİT&nbsp;&nbsp;·&nbsp;&nbsp;4 DOĞAL BOZUNMA SERİSİ</div>"
        )
        self.about_stats_lbl.setText(stats_html)
        if hasattr(self, "about_divider"):
            self.about_divider.apply_theme(self.is_dark)
        for l in (getattr(self, "about_lamp_db", None), getattr(self, "about_lamp_engine", None),
                  getattr(self, "about_lamp_pdf", None), getattr(self, "about_lamp_offline", None)):
            if l is not None:
                l.apply_theme(self.is_dark)

    def create_q_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(10)

        chan_frame = self.chan_frame = QFrame()
        chan_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        chan_layout = QVBoxLayout(chan_frame)
        chan_layout.setContentsMargins(12, 10, 12, 12)
        chan_layout.setSpacing(8)

        chan_title = QLabel("HIZLANDIRICI REAKSİYON KANAL SEÇİCİ   ·   X (a, b) Y")
        chan_title.setObjectName("Engraved")
        chan_layout.addWidget(chan_title)

        selectors_layout = QHBoxLayout()
        selectors_layout.setSpacing(8)

        vbox_x = QVBoxLayout()
        lbl_x = QLabel("HEDEF (X)"); lbl_x.setObjectName("Engraved")
        vbox_x.addWidget(lbl_x)
        self.combo_target = SearchableIsotopeCombo(is_dark=self.is_dark)
        for k in sorted([iso for iso in ISOTOPES.keys() if iso not in ["n", "p"]]):
            self.combo_target.addItem(get_iso_display_name(k), userData=k)
        idx_14n = self.combo_target.findData("14N")
        if idx_14n >= 0: self.combo_target.setCurrentIndex(idx_14n)
        vbox_x.addWidget(self.combo_target)
        selectors_layout.addLayout(vbox_x, 3)

        lbl_plus1 = QLabel("+"); lbl_plus1.setObjectName("Readout")
        selectors_layout.addWidget(lbl_plus1, 0, Qt.AlignCenter)

        vbox_a = QVBoxLayout()
        lbl_a = QLabel("MERMİ (a)"); lbl_a.setObjectName("Engraved")
        vbox_a.addWidget(lbl_a)
        self.combo_proj = PanelComboBox(is_dark=self.is_dark)
        for label, val in LIGHT_PARTICLES:
            self.combo_proj.addItem(label, userData=val)
        idx_alpha = self.combo_proj.findData("alpha")
        if idx_alpha >= 0: self.combo_proj.setCurrentIndex(idx_alpha)
        vbox_a.addWidget(self.combo_proj)
        selectors_layout.addLayout(vbox_a, 2)

        lbl_arrow = QLabel("➔"); lbl_arrow.setObjectName("Readout")
        selectors_layout.addWidget(lbl_arrow, 0, Qt.AlignCenter)

        vbox_b = QVBoxLayout()
        lbl_b = QLabel("EMİSYON (b)"); lbl_b.setObjectName("Engraved")
        vbox_b.addWidget(lbl_b)
        self.combo_emit = PanelComboBox(is_dark=self.is_dark)
        for label, val in LIGHT_PARTICLES:
            self.combo_emit.addItem(label, userData=val)
        idx_p = self.combo_emit.findData("p")
        if idx_p >= 0: self.combo_emit.setCurrentIndex(idx_p)
        vbox_b.addWidget(self.combo_emit)
        selectors_layout.addLayout(vbox_b, 2)

        lbl_plus2 = QLabel("+"); lbl_plus2.setObjectName("Readout")
        selectors_layout.addWidget(lbl_plus2, 0, Qt.AlignCenter)

        vbox_y = QVBoxLayout()
        lbl_y = QLabel("ARTIK ÇEKİRDEK (Y)"); lbl_y.setObjectName("Engraved")
        vbox_y.addWidget(lbl_y)
        self.lbl_residual = QLabel("17O — Oksijen-17")
        self.lbl_residual.setObjectName("Readout")
        self.lbl_residual.setStyleSheet("padding: 5px 8px; border: 1px solid #323943; background: #07090a;")
        vbox_y.addWidget(self.lbl_residual)
        selectors_layout.addLayout(vbox_y, 3)

        chan_layout.addLayout(selectors_layout)
        layout.addWidget(chan_frame)

        self.combo_target.currentIndexChanged.connect(self.handle_channel_changed)
        self.combo_proj.currentIndexChanged.connect(self.handle_channel_changed)
        self.combo_emit.currentIndexChanged.connect(self.handle_channel_changed)

        row = QHBoxLayout(); row.setSpacing(8)
        cap = QLabel("KATALOG REFERANSI"); cap.setObjectName("Engraved")
        row.addWidget(cap)
        self.preset_combo = PanelComboBox(is_dark=self.is_dark)
        for label, val in PRESET_REACTIONS:
            self.preset_combo.addItem(label, userData=val)
        self.preset_combo.currentIndexChanged.connect(self.handle_preset_selected)
        row.addWidget(self.preset_combo, 1)
        layout.addLayout(row)

        cap2 = QLabel("REAKSİYON DENKLEMİ   ·   GİRENLER →  ÜRÜNLER")
        cap2.setObjectName("Engraved")
        layout.addWidget(cap2)

        self.rxn_input = QLineEdit("14N + alpha -> 17O + p")
        self.rxn_input.returnPressed.connect(self.handle_q_calc)
        layout.addWidget(self.rxn_input)

        btn_row = QHBoxLayout()
        calc_btn = QPushButton("HESAPLA  ·  Q-DEĞERİ VE EŞİK KİNEMATİĞİ")
        calc_btn.setObjectName("PrimaryLever")
        calc_btn.clicked.connect(self.handle_q_calc)
        btn_row.addWidget(calc_btn)
        btn_row.addStretch(1)
        btn_row.addWidget(self._make_export_button(
            lambda: self._export_report_pdf(self.rxn_output.toPlainText, "reaksiyon_raporu.pdf")))
        layout.addLayout(btn_row)

        lamps = QHBoxLayout(); lamps.setSpacing(24)
        self.lamp_z   = IndicatorLamp("YÜK  Z  KORUNDU", is_dark=self.is_dark)
        self.lamp_a   = IndicatorLamp("NÜKLEON  A  KORUNDU", is_dark=self.is_dark)
        self.lamp_thr = IndicatorLamp("EŞİK ENERJİSİ GEREKLİ", is_dark=self.is_dark)
        for l in (self.lamp_z, self.lamp_a, self.lamp_thr):
            lamps.addWidget(l)
        lamps.addStretch(1)
        layout.addLayout(lamps)

        result_row = QHBoxLayout(); result_row.setSpacing(10)
        self.rxn_output = QTextEdit()
        self.rxn_output.setReadOnly(True)
        result_row.addWidget(self.rxn_output, 1)

        self.rxn_canvas = MplCanvas()
        result_row.addWidget(self.rxn_canvas, 1)
        layout.addLayout(result_row, 1)

        self.crt_overlay = ScanlineOverlay(self.rxn_output.viewport(), enabled=self.is_dark)
        self.handle_channel_changed()
        return w

    def _draw_energy_diagram(self, r):
        """Girenler/Çıkanlar için basit bir 'reaksiyon koordinatı' enerji seviyesi
        şeması. Q>0: ürünler alt seviyede (enerji açığa çıkar) — Q<0: üst seviyede."""
        p = palette(self.is_dark)
        self.rxn_canvas.apply_theme(self.is_dark)
        ax = self.rxn_canvas.axes
        ax.clear()
        q = r["q_val"]

        if abs(q) < 1e-9:
            ax.hlines(0, 0.5, 9.5, colors=p["accent"], linewidth=3)
            ax.text(5.0, 0.18, "Q = 0   (ELASTİK SAÇILMA)", ha="center",
                    color=p["accent"], fontfamily="monospace", fontsize=9.5)
            ax.set_ylim(-1.2, 1.2)
        else:
            left_y, right_y = 0.0, -q
            up_color = p["green"] if q > 0 else p["red"]
            ax.hlines(left_y, 0.5, 3.5, colors=p["accent"], linewidth=3.4)
            ax.hlines(right_y, 6.5, 9.5, colors=up_color, linewidth=3.4)
            ax.plot([3.5, 6.5], [left_y, right_y], color=p["faint"],
                    linestyle=(0, (4, 3)), linewidth=1.2)
            ax.annotate("", xy=(5.0, right_y), xytext=(5.0, left_y),
                        arrowprops=dict(arrowstyle="-|>", color=up_color, lw=1.6))
            mid_y = (left_y + right_y) / 2.0
            ax.text(5.25, mid_y, f"Q = {q:+.4f} MeV", color=up_color,
                    fontfamily="monospace", fontsize=9, va="center")
            span = max(abs(left_y), abs(right_y), 0.5)
            ax.set_ylim(-span * 1.35 - 0.4, span * 1.35 + 0.4)
            if r.get("e_threshold") is not None:
                loc = "upper left" if right_y > left_y else "lower left"
                self.rxn_canvas.readout_box(
                    f"E_th (lab.) = {r['e_threshold']:.4f} MeV", loc=loc)

        ax.axhline(0, color=p["grid_minor"], lw=0.6, ls=":")
        ax.set_xlim(0, 10)
        ax.set_xticks([2.0, 8.0])
        ax.set_xticklabels(["GİRENLER", "ÇIKANLAR"], fontfamily="monospace", fontsize=9.5)
        ax.set_ylabel("BAĞIL ENERJİ  [MeV]", color=p["label"], fontfamily="monospace", fontsize=9)
        ax.set_facecolor(p["inset"])
        for s in ax.spines.values():
            s.set_color(p["frame"]); s.set_linewidth(1.1)
        ax.tick_params(colors=p["label"], labelsize=8.5, length=0)
        ax.grid(True, axis="y", color=p["grid_minor"], lw=0.5, ls=":")
        self.rxn_canvas.fig.tight_layout()
        self.rxn_canvas.draw()

    def _sync_channel_ui(self):
        p = palette(self.is_dark)
        if hasattr(self, "lbl_residual"):
            self.lbl_residual.setStyleSheet(f"padding: 5px 8px; border: 1px solid {p['seam']}; background: {p['inset']}; color: {p['readout']};")

    def handle_channel_changed(self):
        t = self.combo_target.currentData()
        a = self.combo_proj.currentData()
        b = self.combo_emit.currentData()
        if not t or not a or not b: return

        try:
            res_key, zy, ay, equation = solve_channel(t, a, b)
            # Kanal seçiciden gelen denklemde mermi bellidir -> eşik hesabına aktarılır
            self._channel_eq, self._channel_proj = equation, a
            if res_key:
                disp_name = get_iso_display_name(res_key)
                self.lbl_residual.setText(disp_name)
                self.rxn_input.setText(equation)
                self.handle_q_calc()
            else:
                self.lbl_residual.setText(f"Z={zy}, A={ay} (Bilinmiyor)")
                self.rxn_input.setText(equation)
                self.handle_q_calc()
        except Exception as e:
            self.lbl_residual.setText("İmkânsız")
            self.rxn_output.setPlainText(f" [KANAL HATASI]: {str(e)}")

    def handle_preset_selected(self):
        formula = self.preset_combo.currentData()
        if formula:
            self.rxn_input.setText(formula)
            self.handle_q_calc()

    def handle_q_calc(self):
        text = self.rxn_input.text().strip()
        proj = self._channel_proj if text == getattr(self, "_channel_eq", None) else None
        try:
            r = calculate_reaction(text, projectile=proj)
        except Exception as e:
            self.rxn_output.setPlainText(f" ARIZA · GİRDİ ÇÖZÜMLENEMEDİ\n {e}")
            for l in (self.lamp_z, self.lamp_a, self.lamp_thr):
                l.set_state("idle")
            self.rxn_canvas.apply_theme(self.is_dark)
            self.rxn_canvas.axes.clear()
            self.rxn_canvas.draw()
            return

        self.lamp_z.set_state("ok" if r["z_conserved"] else "fail")
        self.lamp_a.set_state("ok" if r["a_conserved"] else "fail")

        faz_map = {
            "EKZOTERMİK": "EKZOTERMİK · NET ENERJİ SALINIMI",
            "ENDOTERMİK": "ENDOTERMİK · DIŞ KİNETİK ENERJİ GEREKİR",
            "ELASTİK":    "ELASTİK SAÇILMA · Q = 0",
        }
        faz = faz_map[r["phase"]]

        # Lamba: tek çekirdek -> bozunma; iki çekirdek -> X(a,b)Y reaksiyonu.
        # Not: Q>0 bir reaksiyon "kendiliğinden" olmaz; mermi yine de hedefe ulaşmalı
        # (yüklü parçacıkta Coulomb bariyeri). Sadece kinematik eşik yoktur.
        if r["n_in"] == 1:
            if r["phase"] == "EKZOTERMİK":
                self.lamp_thr.set_state("ok", "KENDİLİĞİNDEN BOZUNABİLİR (Q>0)")
            else:
                self.lamp_thr.set_state("fail", "BOZUNMA ENERJETİK OLARAK YASAK (Q≤0)")
        elif r["phase"] == "ENDOTERMİK":
            self.lamp_thr.set_state("fail", "EŞİK ENERJİSİ GEREKLİ")
        elif r["phase"] == "ELASTİK":
            self.lamp_thr.set_state("ok", "ELASTİK · KİNEMATİK EŞİK YOK")
        else:
            self.lamp_thr.set_state("ok", "KİNEMATİK EŞİK YOK (Q>0)")

        rap  = " SİNOP ÜNİVERSİTESİ · NÜKLEER ENERJİ MÜHENDİSLİĞİ\n"
        rap += " REAKSİYON KİNEMATİĞİ ÖLÇÜM FİŞİ" + " " * 11 + "FORM SNU-2/56\n"
        rap += "=" * (W_FIS + 2) + "\n"
        rap += _fis("GİREN KÜTLE TOPLAMI",  f"{r['m_in']:12.6f}",  "u")
        rap += _fis("ÇIKAN KÜTLE TOPLAMI",  f"{r['m_out']:12.6f}", "u")
        rap += _fis("KÜTLE KUSURU  Δm",     f"{r['delta_m']:12.6f}", "u")
        rap += "-" * (W_FIS + 2) + "\n"
        rap += _fis("REAKSİYON Q-DEĞERİ",   f"{r['q_val']:12.4f}", "MeV")
        rap += _fis("TERMODİNAMİK FAZ",     faz)

        if r["e_threshold"] is not None:
            rap += "\n KİNEMATİK EŞİK ANALİZİ\n" + "-" * (W_FIS + 2) + "\n"
            rap += _fis("LABORATUVAR EŞİĞİ  E_th", f"{r['e_threshold']:12.4f}", "MeV")
            rap += _fis("YAKLAŞIK  -Q(1+m_a/m_X)", f"{r['e_threshold_approx']:12.4f}", "MeV")
            rap += ("\n NOT: Hedef çekirdeğin durağan kabul edildiği laboratuvar\n"
                    " sisteminde, gelen merminin kinetik enerjisi bu değerin\n"
                    " altındayken reaksiyon kanalı kapalıdır.\n")

        rap += "\n" + "=" * (W_FIS + 2) + "\n"
        rap += f" ÖLÇÜM {QDateTime.currentDateTime().toString('dd.MM.yyyy HH:mm:ss')}"
        rap += "   ·   OPERATÖR: A.S.ÇETİN\n"
        self.rxn_output.setPlainText(rap)
        self._draw_energy_diagram(r)

    def create_binding_tab(self):
        w = QWidget()
        main_layout = QVBoxLayout(w)
        main_layout.setContentsMargins(20, 20, 20, 20)
        main_layout.setSpacing(12)

        top_layout = QHBoxLayout()
        cap = QLabel("ANALİZ EDİLECEK İZOTOP"); cap.setObjectName("Engraved")
        top_layout.addWidget(cap)

        self.iso_combo = SearchableIsotopeCombo(is_dark=self.is_dark)
        self.iso_combo.setMaxVisibleItems(14)
        
        iso_keys = [k for k in ISOTOPES.keys() if k not in ["n", "p"]]
        for k in iso_keys:
            display_str = get_iso_display_name(k)
            data = ISOTOPES[k]
            mode_info = f" ({data['mode']})" if data['mode'] != '-' else ""
            full_label = f"{display_str}  [Z={data['Z']}, A={data['A']}{mode_info}]"
            self.iso_combo.addItem(full_label, userData=k)

        fe_idx = iso_keys.index("56Fe") if "56Fe" in iso_keys else 0
        self.iso_combo.setCurrentIndex(fe_idx)
        self.iso_combo.currentIndexChanged.connect(self.update_binding_display)
        top_layout.addWidget(self.iso_combo, 1)
        main_layout.addLayout(top_layout)

        content_layout = QHBoxLayout()
        content_layout.setSpacing(12)

        self.card_left = QFrame()
        self.layout_left = QVBoxLayout(self.card_left)
        self.layout_left.setSpacing(10)
        self.layout_left.setContentsMargins(16, 16, 16, 16)
        
        self.iso_info_lbl = QLabel()
        self.iso_info_lbl.setTextFormat(Qt.RichText)
        self.layout_left.addWidget(self.iso_info_lbl)

        self.divider_left = PanelDivider(is_dark=self.is_dark)
        self.layout_left.addWidget(self.divider_left)

        self.iso_table_lbl = QLabel()
        self.iso_table_lbl.setTextFormat(Qt.RichText)
        self.layout_left.addWidget(self.iso_table_lbl)
        self.layout_left.addStretch(1)

        content_layout.addWidget(self.card_left, 1)

        self.card_right = QFrame()
        self.layout_right = QVBoxLayout(self.card_right)
        self.layout_right.setSpacing(12)
        self.layout_right.setContentsMargins(16, 16, 16, 16)

        cap_curve = QLabel("BAĞLANMA ENERJİSİ EĞRİSİ   ·   Eb/A vs A  (KARARLI NÜKLİTLER)")
        cap_curve.setObjectName("Engraved")
        cap_curve.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        self.layout_right.addWidget(cap_curve)

        self.eba_canvas = MplCanvas()
        self.layout_right.addWidget(self.eba_canvas, 1)

        cap_semf = QLabel("SEMF 5-TERİMLİ ENERJİ DENGESİ (MEV)")
        cap_semf.setObjectName("Engraved")
        cap_semf.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        self.layout_right.addWidget(cap_semf)

        self.semf_meter = SemfBarMeter(is_dark=self.is_dark)
        self.layout_right.addWidget(self.semf_meter)

        self.dev_gauge = AnalogGauge("SEMF MODEL SAPMASI", "%", vmax=10.0, warn=5.0, is_dark=self.is_dark)
        self.layout_right.addWidget(self.dev_gauge)

        content_layout.addWidget(self.card_right, 1)
        main_layout.addLayout(content_layout)

        exp_row = QHBoxLayout()
        exp_row.addStretch(1)
        exp_row.addWidget(self._make_export_button(
            lambda: self._export_report_pdf(self._build_binding_report_text, "baglanma_raporu.pdf")))
        main_layout.addLayout(exp_row)

        self.update_binding_display()
        return w

    def update_binding_display(self):
        iso = self.iso_combo.currentData()
        if not iso or iso not in ISOTOPES: return
        data = ISOTOPES[iso]
        Z, A, m_exp = data["Z"], data["A"], data["mass"]
        N = A - Z
        half_life = data["half_life"]
        mode = data["mode"]

        # Bağlanma enerjisi fizik motorundan (Z·M(1H) + N·m_n − M_atom; elektron dengeli)
        be = binding_energy(iso)
        delta_m = be["delta_m"]
        eb_exp = be["Eb"]
        eb_exp_per_a = be["Eb_per_A"]

        sep = separation_energies(iso)
        sn_str = f"{sep['sn']:.3f} MeV ({sep['sn_type']})" if sep and sep['sn'] is not None else "Tanımsız (A=1)"
        sp_str = f"{sep['sp']:.3f} MeV ({sep['sp_type']})" if sep and sep['sp'] is not None else "Tanımsız (Z=0)"

        semf = semf_breakdown(A, Z)
        eb_semf = semf["total"]
        eb_semf_per_a = eb_semf / A if A > 0 else 0
        dev = semf_deviation(iso)
        error = abs(dev) if dev is not None else 0.0

        n_z_ratio = N / Z if Z > 0 else 0
        is_stable = half_life == "Kararlı"
        
        p = palette(self.is_dark)
        status_color = p["green"] if is_stable else p["red"]
        status_text = "KARARLI NÜKLİT" if is_stable else f"RADYOAKTİF ({mode} BOZUNMASI)"
        display_name = get_iso_display_name(iso)

        self.card_left.setStyleSheet(panel_qss(self.is_dark))
        self.card_right.setStyleSheet(panel_qss(self.is_dark))

        self.iso_info_lbl.setText(f"""
        <h3 style='color: {p['accent']}; margin: 0; font-family: monospace; letter-spacing: 1px;'>{display_name}</h3>
        <p style='color: {status_color}; font-weight: bold; margin: 4px 0 0 0; font-family: monospace; font-size: 11px;'>[●] {status_text}</p>
        """)

        self.iso_table_lbl.setText(f"""
        <table style='width: 100%; color: {p['text']}; font-family: monospace; font-size: 12px; line-height: 2.0;'>
            <tr><td style='color: {p['label']};'>Proton Sayısı (Z):</td><td><b>{Z}</b></td></tr>
            <tr><td style='color: {p['label']};'>Nötron Sayısı (N):</td><td><b>{N}</b></td></tr>
            <tr><td style='color: {p['label']};'>Kütle Numarası (A):</td><td><b>{A}</b></td></tr>
            <tr><td style='color: {p['label']};'>Nötron / Proton (N/Z):</td><td><b>{n_z_ratio:.3f}</b></td></tr>
            <tr><td style='color: {p['label']};'>Deneysel Kütle:</td><td>{m_exp:.6f} u</td></tr>
            <tr><td style='color: {p['label']};'>Kütle Kusuru (Δm):</td><td>{delta_m:.6f} u</td></tr>
            <tr><td style='color: {p['label']};'>Yarılanma Ömrü (T½):</td><td><b>{half_life}</b></td></tr>
            <tr><td style='color: {p['label']};'>Deneysel Eb/A:</td><td><b>{eb_exp_per_a:.3f} MeV/nükleon</b></td></tr>
            <tr><td style='color: {p['accent']};'>Nötron Ayrıştırma (Sn):</td><td><b>{sn_str}</b></td></tr>
            <tr><td style='color: {p['accent']};'>Proton Ayrıştırma (Sp):</td><td><b>{sp_str}</b></td></tr>
        </table>
        """)

        self.semf_meter.set_terms([
            ("HACİM  a_v · A",        semf["vol"]),
            ("YÜZEY  a_s · A^2/3",    semf["surf"]),
            ("COULOMB  Z(Z−1)/A^1/3", semf["coul"]),
            ("ASİMETRİ  (A−2Z)²/A",   semf["asym"]),
            ("ÇİFTLENİM  δ",          semf["pair"]),
        ], total=eb_semf)
        self.dev_gauge.set_value(error)
        self._draw_binding_curve(iso, A, eb_exp_per_a)

    def _build_binding_report_text(self) -> str:
        iso = self.iso_combo.currentData()
        if not iso or iso not in ISOTOPES:
            return ""
        data = ISOTOPES[iso]
        Z, A, m_exp = data["Z"], data["A"], data["mass"]
        N = A - Z
        half_life = data["half_life"]
        be = binding_energy(iso)
        sep = separation_energies(iso)
        sn_str = f"{sep['sn']:.3f} MeV ({sep['sn_type']})" if sep and sep['sn'] is not None else "Tanımsız (A=1)"
        sp_str = f"{sep['sp']:.3f} MeV ({sep['sp_type']})" if sep and sep['sp'] is not None else "Tanımsız (Z=0)"
        semf = semf_breakdown(A, Z)
        dev = semf_deviation(iso)

        rap  = " SİNOP ÜNİVERSİTESİ · NÜKLEER ENERJİ MÜHENDİSLİĞİ\n"
        rap += " BAĞLANMA ENERJİSİ ÖLÇÜM FİŞİ" + " " * 12 + "FORM SNU-3/56\n"
        rap += "=" * (W_FIS + 2) + "\n"
        rap += _fis("NÜKLİT", get_iso_display_name(iso))
        rap += _fis("PROTON SAYISI  Z", str(Z))
        rap += _fis("NÖTRON SAYISI  N", str(N))
        rap += _fis("KÜTLE NUMARASI  A", str(A))
        rap += _fis("YARILANMA ÖMRÜ  T½", half_life)
        rap += "-" * (W_FIS + 2) + "\n"
        rap += _fis("KÜTLE KUSURU  Δm", f"{be['delta_m']:12.6f}", "u")
        rap += _fis("DENEYSEL Eb", f"{be['Eb']:12.4f}", "MeV")
        rap += _fis("DENEYSEL Eb/A", f"{be['Eb_per_A']:12.4f}", "MeV/nükleon")
        rap += _fis("NÖTRON AYRIŞTIRMA  Sn", sn_str)
        rap += _fis("PROTON AYRIŞTIRMA  Sp", sp_str)
        rap += "-" * (W_FIS + 2) + "\n"
        rap += _fis("SEMF TOPLAM Eb", f"{semf['total']:12.4f}", "MeV")
        rap += _fis("SEMF SAPMASI", f"{dev:12.3f}", "%") if dev is not None else ""
        rap += "\n" + "=" * (W_FIS + 2) + "\n"
        rap += f" ÖLÇÜM {QDateTime.currentDateTime().toString('dd.MM.yyyy HH:mm:ss')}"
        rap += "   ·   OPERATÖR: A.S.ÇETİN\n"
        return rap

    def _draw_binding_curve(self, iso, A, eb_per_a):
        """Klasik 'Bağlanma Enerjisi Eğrisi': kararlı nüklitler için Eb/A vs A,
        seçili izotop ve ⁵⁶Fe zirvesi vurgulanır."""
        p = palette(self.is_dark)
        self.eba_canvas.apply_theme(self.is_dark)
        ax = self.eba_canvas.axes
        ax.clear()

        curve = get_stable_eba_curve()
        xs = [a for a, _ in curve]
        ys = [e for _, e in curve]
        ax.plot(xs, ys, color=p["faint"], lw=1.1, marker="o", markersize=2.2,
                alpha=0.65, zorder=2, solid_capstyle="round")

        fe_eb = binding_energy("56Fe")["Eb_per_A"]
        ax.scatter([56], [fe_eb], marker="*", s=90, color=p["green"], zorder=4,
                   label="⁵⁶Fe  (en kararlı bölge)")

        ax.scatter([A], [eb_per_a], marker="o", s=80, color=p["accent"],
                   edgecolors=p["accent_hi"], linewidths=1.4, zorder=6)
        ax.annotate(f"{get_iso_display_name(iso)}\n{eb_per_a:.3f} MeV/nükleon",
                    xy=(A, eb_per_a), xytext=(0, 14 if eb_per_a < 8.4 else -22),
                    textcoords="offset points", ha="center",
                    color=p["accent"], fontfamily="monospace", fontsize=8.5)

        ax.set_xlabel("KÜTLE NUMARASI  A", color=p["label"], fontfamily="monospace", fontsize=9)
        ax.set_ylabel("Eb/A  [MeV/nükleon]", color=p["label"], fontfamily="monospace", fontsize=9)
        ax.set_ylim(0, 9.5)
        ax.set_xlim(0, 245)
        ax.set_facecolor(p["inset"])
        for s in ax.spines.values():
            s.set_color(p["frame"]); s.set_linewidth(1.1)
        ax.tick_params(colors=p["label"], labelsize=8.5, length=0)
        ax.grid(True, color=p["grid_minor"], lw=0.5, ls=":")
        ax.legend(loc="lower right", fontsize=7.5, frameon=True,
                        facecolor=p["chassis"] if self.is_dark else p["panel"],
                        edgecolor=p["frame"], labelcolor=p["text"])
        self.eba_canvas.fig.tight_layout()
        self.eba_canvas.draw()

    def create_decay_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(10)

        cap = QLabel("BATEMAN KİNETİK BOZUNMA ANALİZİ   ·   NÜKLİT A ➔ B ➔ C (KARARLI)")
        cap.setObjectName("Engraved")
        layout.addWidget(cap)

        self.canvas = MplCanvas()
        layout.addWidget(self.canvas, 1)

        slider_box = QVBoxLayout()
        slider_box.setSpacing(2)
        
        lbl_box = QHBoxLayout()
        lbl_s = QLabel("NÜKLİT A YARILANMA ÖMRÜ (T½)"); lbl_s.setObjectName("Engraved")
        lbl_box.addWidget(lbl_s)
        self.val_slider_lbl = QLabel("10.0 s")
        self.val_slider_lbl.setObjectName("Readout")
        lbl_box.addStretch(1)
        lbl_box.addWidget(self.val_slider_lbl)
        slider_box.addLayout(lbl_box)

        self.t_slider = QSlider(Qt.Horizontal)
        self.t_slider.setRange(1, 100)
        self.t_slider.setValue(20)
        self.t_slider.valueChanged.connect(self.update_decay_plot)
        slider_box.addWidget(self.t_slider)

        self.decay_scale = EngravedScale(0.5, 50.0, divisions=10, is_dark=self.is_dark)
        slider_box.addWidget(self.decay_scale)
        layout.addLayout(slider_box)

        self.telemetry_frame = QFrame()
        self.telemetry_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=False))
        t_layout = QGridLayout(self.telemetry_frame)
        t_layout.setContentsMargins(12, 8, 12, 8)
        t_layout.setHorizontalSpacing(20)

        c1 = QLabel("BOZUNMA SABİTİ (λ)"); c1.setObjectName("Engraved")
        c2 = QLabel("ORTALAMA ÖMÜR (τ)"); c2.setObjectName("Engraved")
        c3 = QLabel("BAŞLANGIÇ AKTİVİTESİ (A₀)"); c3.setObjectName("Engraved")
        c4 = QLabel("A₀ · CURIE CİNSİNDEN"); c4.setObjectName("Engraved")

        self.lbl_lam = QLabel("0.0693 s⁻¹"); self.lbl_lam.setObjectName("Readout")
        self.lbl_tau = QLabel("14.43 s"); self.lbl_tau.setObjectName("Readout")
        self.lbl_act_bq = QLabel("69.31 Bq"); self.lbl_act_bq.setObjectName("Readout")
        self.lbl_act_ci = QLabel("1.87 nCi"); self.lbl_act_ci.setObjectName("Readout")

        t_layout.addWidget(c1, 0, 0); t_layout.addWidget(self.lbl_lam, 1, 0)
        t_layout.addWidget(c2, 0, 1); t_layout.addWidget(self.lbl_tau, 1, 1)
        t_layout.addWidget(c3, 0, 2); t_layout.addWidget(self.lbl_act_bq, 1, 2)
        t_layout.addWidget(c4, 0, 3); t_layout.addWidget(self.lbl_act_ci, 1, 3)

        layout.addWidget(self.telemetry_frame)

        exp_row = QHBoxLayout()
        exp_row.addStretch(1)
        exp_row.addWidget(self._make_export_button(
            lambda: self._export_report_pdf(self._build_bateman_report_text, "bozunma_kinetigi_raporu.pdf")))
        layout.addLayout(exp_row)

        self.update_decay_plot()
        return w

    def _build_bateman_report_text(self) -> str:
        rap  = " SİNOP ÜNİVERSİTESİ · NÜKLEER ENERJİ MÜHENDİSLİĞİ\n"
        rap += " BATEMAN KİNETİK BOZUNMA ÖLÇÜM FİŞİ" + " " * 6 + "FORM SNU-6/56\n"
        rap += "=" * (W_FIS + 2) + "\n"
        rap += _fis("MODEL", "NÜKLİT A ➔ B ➔ C  (C kararlı ürün)")
        rap += _fis("NÜKLİT A YARILANMA ÖMRÜ  T½", self.val_slider_lbl.text())
        rap += _fis("NÜKLİT B YARILANMA ÖMRÜ  T½", "25.0 s  (sabit)")
        rap += "-" * (W_FIS + 2) + "\n"
        rap += _fis("BOZUNMA SABİTİ  λ", self.lbl_lam.text())
        rap += _fis("ORTALAMA ÖMÜR  τ", self.lbl_tau.text())
        rap += _fis("BAŞLANGIÇ AKTİVİTESİ  A₀", self.lbl_act_bq.text())
        rap += _fis("A₀ · CURIE CİNSİNDEN", self.lbl_act_ci.text())
        rap += "\n" + "=" * (W_FIS + 2) + "\n"
        rap += f" ÖLÇÜM {QDateTime.currentDateTime().toString('dd.MM.yyyy HH:mm:ss')}"
        rap += "   ·   OPERATÖR: A.S.ÇETİN\n"
        return rap

    def update_decay_plot(self):
        t_half_A = self.t_slider.value() / 2.0
        if hasattr(self, "val_slider_lbl"):
            self.val_slider_lbl.setText(f"{t_half_A:.1f} s")
        t, N_A, N_B, N_C = calculate_bateman(t_half_A)

        lam_A = math.log(2) / max(0.001, t_half_A)
        tau_A = 1.0 / lam_A
        a0_bq = lam_A * 1000.0
        a0_nci = (a0_bq / 3.7e10) * 1e9

        if hasattr(self, "lbl_lam"):
            self.lbl_lam.setText(f"{lam_A:.4f} s⁻¹")
            self.lbl_tau.setText(f"{tau_A:.2f} s")
            self.lbl_act_bq.setText(f"{a0_bq:.2f} Bq")
            self.lbl_act_ci.setText(f"{a0_nci:.3f} nCi")

        self.canvas.apply_theme(self.is_dark)
        p = palette(self.is_dark)
        ax = self.canvas.axes
        ax.clear()

        self.canvas.trace(t, N_A, p["accent"], f"NÜKLİT A   T½ ={t_half_A:6.1f} s")
        self.canvas.trace(t, N_B, p["blue"],    "NÜKLİT B   T½ =  25.0 s")
        self.canvas.trace(t, N_C, p["green"],   "NÜKLİT C   kararlı ürün")

        ax.set_xlim(float(t[0]), float(t[-1]))
        ax.set_ylim(0, None)
        ax.set_xlabel("ZAMAN  [s]", color=p["label"], fontfamily="monospace", fontsize=9)
        ax.set_ylabel("POPÜLASYON  [çekirdek]", color=p["label"], fontfamily="monospace", fontsize=9)

        self.canvas.apply_graticule()
        self.canvas.style_legend()
        x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
        self.canvas.readout_box(f"H: {(x1-x0)/10:5.1f} s/div\nV: {(y1-y0)/8:5.0f} /div")
        self.canvas.fig.tight_layout()
        self.canvas.draw()

    def create_decay_modes_tab(self):
        w = QWidget()
        main_layout = QVBoxLayout(w)
        main_layout.setContentsMargins(20, 20, 20, 20)
        main_layout.setSpacing(12)

        top_layout = QHBoxLayout()
        cap = QLabel("BOZUNMA MODLARI ANALİZ EDİLECEK İZOTOP"); cap.setObjectName("Engraved")
        top_layout.addWidget(cap)

        self.dm_combo = SearchableIsotopeCombo(is_dark=self.is_dark)
        iso_keys = [k for k in ISOTOPES.keys() if k not in ("n", "p")]
        for k in iso_keys:
            display_str = get_iso_display_name(k)
            data = ISOTOPES[k]
            mode_info = f" ({data['mode']})" if data['mode'] != '-' else ""
            self.dm_combo.addItem(f"{display_str}  [Z={data['Z']}, A={data['A']}{mode_info}]", userData=k)
        co60_idx = iso_keys.index("60Co") if "60Co" in iso_keys else 0
        self.dm_combo.setCurrentIndex(co60_idx)
        self.dm_combo.currentIndexChanged.connect(self.update_decay_modes_display)
        top_layout.addWidget(self.dm_combo, 1)
        main_layout.addLayout(top_layout)

        content_layout = QHBoxLayout()
        content_layout.setSpacing(12)

        self.dm_card_left = QFrame()
        self.dm_layout_left = QVBoxLayout(self.dm_card_left)
        self.dm_layout_left.setSpacing(10)
        self.dm_layout_left.setContentsMargins(16, 16, 16, 16)

        self.dm_info_lbl = QLabel()
        self.dm_info_lbl.setTextFormat(Qt.RichText)
        self.dm_layout_left.addWidget(self.dm_info_lbl)

        self.dm_divider = PanelDivider(is_dark=self.is_dark)
        self.dm_layout_left.addWidget(self.dm_divider)

        self.dm_table_lbl = QLabel()
        self.dm_table_lbl.setTextFormat(Qt.RichText)
        self.dm_layout_left.addWidget(self.dm_table_lbl)
        self.dm_layout_left.addStretch(1)

        content_layout.addWidget(self.dm_card_left, 1)

        self.dm_card_right = QFrame()
        self.dm_layout_right = QVBoxLayout(self.dm_card_right)
        self.dm_layout_right.setSpacing(12)
        self.dm_layout_right.setContentsMargins(16, 16, 16, 16)

        cap_chart = QLabel("KANAL BAZLI Q-DEĞERLERİ   ·   YEŞİL = ENERJİK OLARAK MÜMKÜN (Q>0)")
        cap_chart.setObjectName("Engraved")
        cap_chart.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        self.dm_layout_right.addWidget(cap_chart)

        self.dm_canvas = MplCanvas()
        self.dm_layout_right.addWidget(self.dm_canvas, 1)

        content_layout.addWidget(self.dm_card_right, 1)
        main_layout.addLayout(content_layout)

        exp_row = QHBoxLayout()
        exp_row.addStretch(1)
        exp_row.addWidget(self._make_export_button(
            lambda: self._export_report_pdf(self._build_decay_modes_report_text, "bozunma_modlari_raporu.pdf")))
        main_layout.addLayout(exp_row)

        self.update_decay_modes_display()
        return w

    def _build_decay_modes_report_text(self) -> str:
        iso = self.dm_combo.currentData()
        if not iso or iso not in ISOTOPES:
            return ""
        data = ISOTOPES[iso]
        Z, A = data["Z"], data["A"]
        N = A - Z
        half_life = data["half_life"]

        rap  = " SİNOP ÜNİVERSİTESİ · NÜKLEER ENERJİ MÜHENDİSLİĞİ\n"
        rap += " BOZUNMA MODLARI ÖLÇÜM FİŞİ" + " " * 14 + "FORM SNU-5/56\n"
        rap += "=" * (W_FIS + 2) + "\n"
        rap += _fis("NÜKLİT", get_iso_display_name(iso))
        rap += _fis("PROTON SAYISI  Z", str(Z))
        rap += _fis("NÖTRON SAYISI  N", str(N))
        rap += _fis("KÜTLE NUMARASI  A", str(A))
        rap += _fis("YARILANMA ÖMRÜ  T½", half_life)
        rap += _fis("NUBASE MODU", data.get("mode", "-"))
        rap += "-" * (W_FIS + 2) + "\n"
        rap += " KANAL BAZLI Q-DEĞERLERİ (AME2020 ATOMİK KÜTLELERİ)\n"
        rap += "-" * (W_FIS + 2) + "\n"
        channels = decay_channels(iso) or {}
        for key in ("alpha", "beta_minus", "beta_plus", "ec"):
            ch = channels.get(key)
            if not ch:
                continue
            durum = "MÜMKÜN" if ch["q"] > 0 else "YASAK"
            rap += _fis(f"{ch['label']} → {get_iso_display_name(ch['daughter'])}",
                        f"{ch['q']:+.4f} MeV  ({durum})")
        sa = specific_activity(iso)
        rap += "-" * (W_FIS + 2) + "\n"
        if sa:
            rap += _fis("BOZUNMA SABİTİ  λ", f"{sa['lambda']:.4e}", "s⁻¹")
            rap += _fis("SPESİFİK AKTİVİTE", f"{sa['bq_per_g']:.4e}", "Bq/g")
            rap += _fis("SPESİFİK AKTİVİTE", f"{sa['ci_per_g']:.4e}", "Ci/g")
        else:
            rap += _fis("SPESİFİK AKTİVİTE", "Tanımsız (kararlı nüklit)")
        rap += "\n" + "=" * (W_FIS + 2) + "\n"
        rap += f" ÖLÇÜM {QDateTime.currentDateTime().toString('dd.MM.yyyy HH:mm:ss')}"
        rap += "   ·   OPERATÖR: A.S.ÇETİN\n"
        return rap

    def update_decay_modes_display(self):
        iso = self.dm_combo.currentData()
        if not iso or iso not in ISOTOPES:
            return
        data = ISOTOPES[iso]
        Z, A = data["Z"], data["A"]
        N = A - Z
        half_life = data["half_life"]
        mode = data["mode"]
        is_stable = half_life == "Kararlı"

        p = palette(self.is_dark)
        status_color = p["green"] if is_stable else p["red"]
        status_text = "KARARLI NÜKLİT" if is_stable else f"RADYOAKTİF · NUBASE MODU: {mode}"
        display_name = get_iso_display_name(iso)

        self.dm_card_left.setStyleSheet(panel_qss(self.is_dark))
        self.dm_card_right.setStyleSheet(panel_qss(self.is_dark))

        self.dm_info_lbl.setText(f"""
        <h3 style='color: {p['accent']}; margin: 0; font-family: monospace; letter-spacing: 1px;'>{display_name}</h3>
        <p style='color: {status_color}; font-weight: bold; margin: 4px 0 0 0; font-family: monospace; font-size: 11px;'>[●] {status_text}</p>
        """)

        sa = specific_activity(iso)
        if sa:
            sa_rows = f"""
            <tr><td style='color: {p['accent']};'>Bozunma Sabiti (λ):</td><td><b>{sa['lambda']:.4e} s⁻¹</b></td></tr>
            <tr><td style='color: {p['accent']};'>Spesifik Aktivite:</td><td><b>{sa['bq_per_g']:.4e} Bq/g</b></td></tr>
            <tr><td style='color: {p['accent']};'>Spesifik Aktivite:</td><td><b>{sa['ci_per_g']:.4e} Ci/g</b></td></tr>
            """
        else:
            sa_rows = f"<tr><td style='color: {p['label']};'>Spesifik Aktivite:</td><td>Tanımsız (kararlı nüklit)</td></tr>"

        self.dm_table_lbl.setText(f"""
        <table style='width: 100%; color: {p['text']}; font-family: monospace; font-size: 12px; line-height: 2.0;'>
            <tr><td style='color: {p['label']};'>Proton Sayısı (Z):</td><td><b>{Z}</b></td></tr>
            <tr><td style='color: {p['label']};'>Nötron Sayısı (N):</td><td><b>{N}</b></td></tr>
            <tr><td style='color: {p['label']};'>Kütle Numarası (A):</td><td><b>{A}</b></td></tr>
            <tr><td style='color: {p['label']};'>Yarılanma Ömrü (T½):</td><td><b>{half_life}</b></td></tr>
            {sa_rows}
        </table>
        """)

        self._draw_decay_modes_chart(iso)

    def _draw_decay_modes_chart(self, iso):
        p = palette(self.is_dark)
        self.dm_canvas.apply_theme(self.is_dark)
        ax = self.dm_canvas.axes
        ax.clear()

        channels = decay_channels(iso) or {}
        if not channels:
            ax.text(0.5, 0.5, "BU NÜKLİT İÇİN TANIMLI KANAL YOK\n(Z≤2 VEYA ÜRÜN VERİTABANINDA YOK)",
                    ha="center", va="center", transform=ax.transAxes,
                    color=p["faint"], fontfamily="monospace", fontsize=9)
        else:
            order = [k for k in ("alpha", "beta_minus", "beta_plus", "ec") if k in channels]
            labels = [channels[k]["label"] + f"\n→ {get_iso_display_name(channels[k]['daughter'])}" for k in order]
            qs = [channels[k]["q"] for k in order]
            colors = [p["green"] if q > 0 else p["red"] for q in qs]
            y_pos = list(range(len(order)))

            ax.barh(y_pos, qs, color=colors, height=0.55, zorder=3)
            ax.set_yticks(y_pos)
            ax.set_yticklabels(labels, fontfamily="monospace", fontsize=8.5)
            ax.invert_yaxis()
            ax.axvline(0, color=p["frame"], lw=1.2)

            max_abs = max(abs(q) for q in qs) or 1.0
            for y, q in zip(y_pos, qs):
                if abs(q) > 0.22 * max_abs:
                    # Çubuk yeterince geniş: etiketi çubuğun içine, koyu renkle yerleştir.
                    ax.text(q / 2.0, y, f"{q:+.4f} MeV", va="center", ha="center",
                            color=p["chassis"], fontfamily="monospace", fontsize=8.5,
                            fontweight="bold", zorder=4)
                else:
                    # Çubuk çok kısa: etiketi ucunun dışına, panel metin rengiyle yerleştir.
                    ax.text(q, y, f"  {q:+.4f} MeV  ", va="center",
                            ha="left" if q >= 0 else "right",
                            color=p["text"], fontfamily="monospace", fontsize=8.5, zorder=4)

            ax.set_xlabel("Q-DEĞERİ  [MeV]", color=p["label"], fontfamily="monospace", fontsize=9)

        ax.set_facecolor(p["inset"])
        for s in ax.spines.values():
            s.set_color(p["frame"]); s.set_linewidth(1.1)
        ax.tick_params(colors=p["label"], labelsize=8.5, length=0)
        ax.grid(True, axis="x", color=p["grid_minor"], lw=0.5, ls=":")
        # tight_layout(), en uzun kanal etiketi ("ELEKTRON YAKALAMA (EC)\n→ ...")
        # bazı pencere genişliklerinde soldan kırpılmasına yol açabiliyor;
        # sabit, cömert bir sol pay ayırmak bunu pencere boyutundan bağımsız çözer.
        self.dm_canvas.fig.subplots_adjust(left=0.34, right=0.97, top=0.90, bottom=0.16)
        self.dm_canvas.draw()

    def create_decay_series_tab(self):
        w = QWidget()
        main_layout = QVBoxLayout(w)
        main_layout.setContentsMargins(20, 20, 20, 20)
        main_layout.setSpacing(12)

        top_layout = QHBoxLayout()
        cap = QLabel("BOZUNMA SERİSİ (DOĞAL ZİNCİR) SEÇİMİ"); cap.setObjectName("Engraved")
        top_layout.addWidget(cap)

        self.ds_combo = PanelComboBox(is_dark=self.is_dark)
        for start_key, info in DECAY_SERIES.items():
            self.ds_combo.addItem(
                f"{info['name']}  ·  {get_iso_display_name(start_key)} → {get_iso_display_name(info['final'])}",
                userData=start_key
            )
        self.ds_combo.currentIndexChanged.connect(self.update_decay_series_display)
        top_layout.addWidget(self.ds_combo, 1)
        main_layout.addLayout(top_layout)

        content_layout = QHBoxLayout()
        content_layout.setSpacing(12)

        self.ds_frame_left = QFrame()
        left_l = QVBoxLayout(self.ds_frame_left)
        left_l.setContentsMargins(12, 10, 12, 12)
        left_l.setSpacing(8)

        cap_left = QLabel("ZİNCİR ADIMLARI  (BASKIN NUBASE MODUNA GÖRE)")
        cap_left.setObjectName("Engraved")
        cap_left.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        left_l.addWidget(cap_left)

        self.ds_output = QTextEdit()
        self.ds_output.setReadOnly(True)
        left_l.addWidget(self.ds_output, 1)

        content_layout.addWidget(self.ds_frame_left, 1)

        self.ds_frame_right = QFrame()
        right_l = QVBoxLayout(self.ds_frame_right)
        right_l.setContentsMargins(12, 10, 12, 12)
        right_l.setSpacing(8)

        cap_right = QLabel("ZİNCİR HARİTASI  ·  A (KÜTLE NO) vs Z (PROTON NO)")
        cap_right.setObjectName("Engraved")
        cap_right.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        right_l.addWidget(cap_right)

        self.ds_canvas = MplCanvas()
        right_l.addWidget(self.ds_canvas, 1)

        content_layout.addWidget(self.ds_frame_right, 1)
        main_layout.addLayout(content_layout, 1)

        exp_row = QHBoxLayout()
        exp_row.addStretch(1)
        exp_row.addWidget(self._make_export_button(
            lambda: self._export_report_pdf(self._build_decay_series_report_text, "bozunma_serisi_raporu.pdf")))
        main_layout.addLayout(exp_row)

        self._ds_chain = []
        self.update_decay_series_display()
        return w

    def update_decay_series_display(self):
        start_key = self.ds_combo.currentData()
        if not start_key:
            return
        info = DECAY_SERIES.get(start_key, {})
        chain = decay_series_chain(start_key)
        self._ds_chain = chain

        p = palette(self.is_dark)
        self.ds_frame_left.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        self.ds_frame_right.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))

        html = [f"<div style='font-family:monospace; font-size:11.5px; color:{p['text']}; line-height:1.55;'>"]
        html.append(f"<p style='color:{p['accent']}; font-weight:bold; font-size:13px; margin-bottom:6px;'>{info.get('name', '')}</p>")
        mode_color = {"α": p["red"], "β-": p["blue"], "β+": p["accent"], "EC": p["accent"]}
        for i, n in enumerate(chain):
            if i > 0:
                mode_disp = chain[i - 1]["step_mode"] or "?"
                col = mode_color.get(mode_disp, p["label"])
                q_disp = f"Q={n['q']:+.3f} MeV" if n['q'] is not None else ""
                html.append(f"<p style='color:{col}; margin:2px 0;'>&nbsp;&nbsp;↓&nbsp; {mode_disp} &nbsp; {q_disp}</p>")
            status_color = p["green"] if n["stable"] else p["text"]
            status_tag = "  [KARARLI — ZİNCİR SONU]" if n["stable"] else ""
            html.append(
                f"<p style='margin:1px 0;'><b style='color:{status_color};'>{n['name']}</b>"
                f"&nbsp;&nbsp;Z={n['Z']} &nbsp;A={n['A']} &nbsp;T½={n['half_life']}{status_tag}</p>"
            )
            if n.get("terminal_reason"):
                html.append(f"<p style='color:{p['red']};'>&nbsp;&nbsp;⚠ {n['terminal_reason']}</p>")
        html.append(f"<p style='color:{p['label']}; margin-top:8px;'>TOPLAM ADIM: {max(0, len(chain) - 1)}</p>")
        html.append("</div>")
        self.ds_output.setHtml("".join(html))

        self._draw_decay_series_chart(chain)

    def _draw_decay_series_chart(self, chain):
        p = palette(self.is_dark)
        self.ds_canvas.apply_theme(self.is_dark)
        ax = self.ds_canvas.axes
        ax.clear()

        if not chain:
            ax.text(0.5, 0.5, "VERİ YOK", ha="center", va="center",
                     transform=ax.transAxes, color=p["faint"], fontfamily="monospace")
            self.ds_canvas.draw()
            return

        mode_colors = {"α": p["red"], "β-": p["blue"], "β+": p["accent"], "EC": p["accent"]}

        for i in range(1, len(chain)):
            prev, cur = chain[i - 1], chain[i]
            col = mode_colors.get(prev["step_mode"], p["label"])
            ax.plot([prev["Z"], cur["Z"]], [prev["A"], cur["A"]], color=col, lw=1.8, zorder=2)

        # Ardışık düğümler sık sık aynı A satırında (β⁻: Z+1, A=) yan yana düşer;
        # etiketleri her düğümde farklı bir yöne kaydırarak çakışmayı azaltıyoruz.
        offsets = [(7, 7), (7, -11), (-9, 7), (-9, -11)]
        for i, n in enumerate(chain):
            col = p["green"] if n["stable"] else p["faint"]
            ax.scatter([n["Z"]], [n["A"]], s=70 if n["stable"] else 42,
                        color=col, edgecolors=p["frame"], linewidths=0.9, zorder=3)
            ox, oy = offsets[i % len(offsets)]
            ha = "left" if ox >= 0 else "right"
            ax.annotate(n["name"], (n["Z"], n["A"]), textcoords="offset points",
                         xytext=(ox, oy), fontsize=7.3, fontfamily="monospace",
                         color=p["text"], ha=ha, zorder=4)

        from matplotlib.lines import Line2D
        legend_handles = [
            Line2D([0], [0], color=p["red"], lw=2, label="α  (Z−2, A−4)"),
            Line2D([0], [0], color=p["blue"], lw=2, label="β⁻  (Z+1, A=)"),
        ]
        leg = ax.legend(handles=legend_handles, loc="lower left", fontsize=7.5,
                          framealpha=0.9, facecolor=p["inset"], edgecolor=p["frame"])
        for txt in leg.get_texts():
            txt.set_color(p["text"])
            txt.set_fontfamily("monospace")

        ax.set_xlabel("PROTON SAYISI  Z", color=p["label"], fontfamily="monospace", fontsize=9)
        ax.set_ylabel("KÜTLE NUMARASI  A", color=p["label"], fontfamily="monospace", fontsize=9)
        ax.set_facecolor(p["inset"])
        for s in ax.spines.values():
            s.set_color(p["frame"]); s.set_linewidth(1.1)
        ax.tick_params(colors=p["label"], labelsize=8)
        ax.grid(True, color=p["grid_minor"], lw=0.5, ls=":")
        ax.margins(x=0.16, y=0.09)

        self.ds_canvas.fig.subplots_adjust(left=0.14, right=0.96, top=0.94, bottom=0.14)
        self.ds_canvas.draw()

    def _build_decay_series_report_text(self) -> str:
        start_key = self.ds_combo.currentData()
        if not start_key:
            return ""
        info = DECAY_SERIES.get(start_key, {})
        chain = self._ds_chain or decay_series_chain(start_key)

        rap  = " SİNOP ÜNİVERSİTESİ · NÜKLEER ENERJİ MÜHENDİSLİĞİ\n"
        rap += " BOZUNMA SERİSİ ÖLÇÜM FİŞİ" + " " * 15 + "FORM SNU-8/56\n"
        rap += "=" * (W_FIS + 2) + "\n"
        rap += _fis("SERİ", info.get("name", "-"))
        rap += _fis("BAŞLANGIÇ NÜKLİDİ", get_iso_display_name(start_key))
        rap += _fis("BEKLENEN SON ÜRÜN", get_iso_display_name(info.get("final", "-")))
        rap += _fis("TOPLAM ADIM SAYISI", str(max(0, len(chain) - 1)))
        rap += "-" * (W_FIS + 2) + "\n"
        for i, n in enumerate(chain):
            rap += _fis(f"{i+1:02d}. {n['name']}", f"Z={n['Z']} A={n['A']}")
            rap += _fis("   T½", n["half_life"])
            if i > 0:
                prev = chain[i - 1]
                q_disp = f"{prev['q']:+.4f} MeV" if prev["q"] is not None else "—"
                rap += _fis(f"   ← {prev['step_mode'] or '?'} BOZUNMASI  Q", q_disp)
            if n.get("terminal_reason"):
                rap += _fis("   NOT", n["terminal_reason"])
        rap += "\n" + "=" * (W_FIS + 2) + "\n"
        rap += f" ÖLÇÜM {QDateTime.currentDateTime().toString('dd.MM.yyyy HH:mm:ss')}"
        rap += "   ·   OPERATÖR: A.S.ÇETİN\n"
        return rap

    def create_db_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(10)

        search_frame = self.search_frame = QFrame()
        search_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        s_layout = QVBoxLayout(search_frame)
        s_layout.setContentsMargins(12, 10, 12, 10)
        s_layout.setSpacing(8)

        s_title = QLabel("NÜKLİT ARŞİVİ VE İZOTOPİK ARAMA PANELİ  (AME2020 / NUBASE ÇEVRİMDIŞI KÜTÜPHANESİ)")
        s_title.setObjectName("Engraved")
        s_layout.addWidget(s_title)

        s_bar_layout = QHBoxLayout()
        s_bar_layout.setSpacing(10)

        self.db_search_input = QLineEdit()
        self.db_search_input.setPlaceholderText("İzotop, Element veya Özellik Ara (Örn: 131I, 137Cs, U, Pu, 60Co, Kararlı, alfa)...")
        self.db_search_input.textChanged.connect(self.filter_db_table)
        s_bar_layout.addWidget(self.db_search_input, 1)

        self.db_lamp = IndicatorLamp("YEREL AME2020 AKTİF", is_dark=self.is_dark)
        self.db_lamp.set_state("ok")
        s_bar_layout.addWidget(self.db_lamp)

        s_layout.addLayout(s_bar_layout)

        self.lbl_db_telemetry = QLabel(f"SİSTEMDE {len(ISOTOPES)} NÜKLİT KAYITLI (TAM ÇEVRİMDIŞI ÇALIŞMA)")
        self.lbl_db_telemetry.setObjectName("Engraved")
        s_layout.addWidget(self.lbl_db_telemetry)

        layout.addWidget(search_frame)

        self.db_table = QTableWidget(len(ISOTOPES), 5)
        self.db_table.setHorizontalHeaderLabels(["NÜKLİT KODU", "Z (PROTON)", "A (KÜTLE)", "KÜTLE [u]", "YARILANMA ÖMRÜ"])
        self.db_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.populate_db_table()
        layout.addWidget(self.db_table, 1)

        return w

    def populate_db_table(self):
        self.db_table.setRowCount(len(ISOTOPES))
        for row, (k, v) in enumerate(ISOTOPES.items()):
            self.db_table.setItem(row, 0, QTableWidgetItem(k))
            self.db_table.setItem(row, 1, QTableWidgetItem(str(v["Z"])))
            self.db_table.setItem(row, 2, QTableWidgetItem(str(v["A"])))
            self.db_table.setItem(row, 3, QTableWidgetItem(f"{v['mass']:.6f}"))
            self.db_table.setItem(row, 4, QTableWidgetItem(v["half_life"]))
        if hasattr(self, "lbl_db_telemetry"):
            self.lbl_db_telemetry.setText(f"SİSTEMDE {len(ISOTOPES)} NÜKLİT KAYITLI (TAM ÇEVRİMDIŞI ÇALIŞMA)")

    def filter_db_table(self):
        text = self.db_search_input.text().strip().lower()
        visible_count = 0
        for row in range(self.db_table.rowCount()):
            match = False
            for col in range(self.db_table.columnCount()):
                item = self.db_table.item(row, col)
                if item and text in item.text().lower():
                    match = True
                    break
            self.db_table.setRowHidden(row, not match)
            if match: visible_count += 1
        
        if text:
            self.lbl_db_telemetry.setText(f"FİLTRELENEN: {visible_count} / {len(ISOTOPES)} NÜKLİT")
        else:
            self.lbl_db_telemetry.setText(f"SİSTEMDE {len(ISOTOPES)} NÜKLİT KAYITLI (TAM ÇEVRİMDIŞI ÇALIŞMA)")

    def create_chart_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(10)

        self.chart_frame = QFrame()
        self.chart_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        f_layout = QVBoxLayout(self.chart_frame)
        f_layout.setContentsMargins(12, 10, 12, 12)
        f_layout.setSpacing(8)

        title = QLabel("NÜKLİT HARİTASI   ·   N-Z KARARLILIK VADİSİ  (SEGRE ÇİZELGESİ)")
        title.setObjectName("Engraved")
        f_layout.addWidget(title)

        row = QHBoxLayout(); row.setSpacing(8)
        lbl = QLabel("HARİTADA VURGULA"); lbl.setObjectName("Engraved")
        row.addWidget(lbl)
        self.chart_combo = SearchableIsotopeCombo(is_dark=self.is_dark)
        self.chart_combo.addItem("— Seçilmedi —", userData=None)
        for k in sorted([iso for iso in ISOTOPES.keys() if iso not in ("n", "p")]):
            self.chart_combo.addItem(get_iso_display_name(k), userData=k)
        idx_fe = self.chart_combo.findData("56Fe")
        if idx_fe >= 0: self.chart_combo.setCurrentIndex(idx_fe)
        self.chart_combo.currentIndexChanged.connect(self.handle_chart_selection)
        row.addWidget(self.chart_combo, 1)

        self.chart_reset_btn = QPushButton("SIFIRLA GÖRÜNÜM")
        self.chart_reset_btn.setObjectName("ModeSwitch")
        self.chart_reset_btn.setCursor(Qt.PointingHandCursor)
        self.chart_reset_btn.clicked.connect(self._reset_segre_view)
        row.addWidget(self.chart_reset_btn)
        f_layout.addLayout(row)

        hint = QLabel("İPUCU: Bir noktaya tıklayarak künyesini açın  ·  kaydırma tekerleği ile yakınlaştırın  ·  sürükleyerek kaydırın.")
        hint.setStyleSheet(f"color: {palette(self.is_dark)['faint']}; font-family: monospace; font-size: 9.5px;")
        f_layout.addWidget(hint)
        self.chart_hint_lbl = hint
        layout.addWidget(self.chart_frame)

        content_row = QHBoxLayout()
        content_row.setSpacing(12)

        self.segre_canvas = MplCanvas()
        # Sürükleme/yakınlaştırma durumu (izleyici -> _draw_segre_chart bunu
        # koruyarak yeniden çizer, aksi hâlde her seçimde/tema değişiminde
        # görünüm sıfırlanırdı).
        self._segre_xlim = None
        self._segre_ylim = None
        self._segre_press_px = None
        self._segre_last_px = None
        self._segre_press_data = None
        self._segre_dragging = False
        self.segre_canvas.mpl_connect("button_press_event", self._on_segre_press)
        self.segre_canvas.mpl_connect("motion_notify_event", self._on_segre_motion)
        self.segre_canvas.mpl_connect("button_release_event", self._on_segre_release)
        self.segre_canvas.mpl_connect("scroll_event", self._on_segre_scroll)
        content_row.addWidget(self.segre_canvas, 3)

        content_row.addWidget(self._build_nuclide_info_panel(), 2)

        layout.addLayout(content_row, 1)

        self._draw_segre_chart()
        self._update_nuclide_info_panel("56Fe")
        return w

    def _reset_segre_view(self):
        self._segre_xlim = None
        self._segre_ylim = None
        self._draw_segre_chart()

    def _select_from_data_coords(self, n0, z0):
        """(N, Z) veri koordinatına en yakın izotopu bulup künye panelini
        günceller — hem doğrudan tıklama hem de sürükleme-olmayan basma/bırakma
        çiftinden çağrılır."""
        if n0 is None or z0 is None:
            return
        lookup = get_nz_lookup()
        best_key, best_dist = None, None
        for dn in range(-2, 3):
            for dz in range(-2, 3):
                cand = (round(n0) + dn, round(z0) + dz)
                key = lookup.get(cand)
                if key is None:
                    continue
                dist = (cand[0] - n0) ** 2 + (cand[1] - z0) ** 2
                if best_dist is None or dist < best_dist:
                    best_dist, best_key = dist, key
        if best_key is None or best_dist > 4.0:
            return
        idx = self.chart_combo.findData(best_key)
        if idx >= 0:
            self.chart_combo.setCurrentIndex(idx)   # -> handle_chart_selection tetiklenir
        else:
            self._update_nuclide_info_panel(best_key)

    def _on_segre_press(self, event):
        if event.button != 1 or event.inaxes != self.segre_canvas.axes:
            self._segre_press_px = None
            return
        self._segre_press_px = (event.x, event.y)
        self._segre_last_px = (event.x, event.y)
        self._segre_press_data = (event.xdata, event.ydata)
        self._segre_dragging = False

    def _on_segre_motion(self, event):
        if self._segre_press_px is None or event.x is None or event.y is None:
            return
        dx_px = event.x - self._segre_press_px[0]
        dy_px = event.y - self._segre_press_px[1]
        if not self._segre_dragging and (dx_px * dx_px + dy_px * dy_px) > 16:
            self._segre_dragging = True
        if not self._segre_dragging:
            return
        ax = self.segre_canvas.axes
        inv = ax.transData.inverted()
        xd0, yd0 = inv.transform(self._segre_last_px)
        xd1, yd1 = inv.transform((event.x, event.y))
        dx, dy = xd0 - xd1, yd0 - yd1
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        ax.set_xlim(xlim[0] + dx, xlim[1] + dx)
        ax.set_ylim(ylim[0] + dy, ylim[1] + dy)
        self._segre_xlim, self._segre_ylim = ax.get_xlim(), ax.get_ylim()
        self.segre_canvas.draw_idle()
        self._segre_last_px = (event.x, event.y)

    def _on_segre_release(self, event):
        if self._segre_press_px is not None and not self._segre_dragging:
            self._select_from_data_coords(*self._segre_press_data)
        self._segre_press_px = None
        self._segre_last_px = None
        self._segre_press_data = None
        self._segre_dragging = False

    def _on_segre_scroll(self, event):
        ax = self.segre_canvas.axes
        if event.inaxes != ax or event.xdata is None or event.ydata is None:
            return
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        zoom_in = event.button == "up"
        scale = 0.85 if zoom_in else (1.0 / 0.85)
        new_w = (xlim[1] - xlim[0]) * scale
        new_h = (ylim[1] - ylim[0]) * scale
        new_w = max(6.0, min(170.0, new_w))
        new_h = max(4.0, min(110.0, new_h))
        relx = (event.xdata - xlim[0]) / (xlim[1] - xlim[0])
        rely = (event.ydata - ylim[0]) / (ylim[1] - ylim[0])
        new_xlim = (event.xdata - relx * new_w, event.xdata + (1 - relx) * new_w)
        new_ylim = (event.ydata - rely * new_h, event.ydata + (1 - rely) * new_h)
        ax.set_xlim(new_xlim)
        ax.set_ylim(new_ylim)
        self._segre_xlim, self._segre_ylim = new_xlim, new_ylim
        self.segre_canvas.draw_idle()

    def _build_nuclide_info_panel(self):
        """Nüklit Haritası'nda tıklanan izotop için künye + fiziksel/genel
        açıklama paneli (sağ sabit yan panel). Taşma ihtimaline karşı
        QScrollArea içine alınır (Künye sekmesindeki aynı önlem)."""
        self.nuclide_panel_frame = QFrame()
        self.nuclide_panel_frame.setMinimumWidth(320)
        self.nuclide_panel_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        outer = QVBoxLayout(self.nuclide_panel_frame)
        outer.setContentsMargins(12, 10, 12, 12)
        outer.setSpacing(8)

        title = QLabel("İZOTOP KÜNYESİ"); title.setObjectName("Engraved")
        outer.addWidget(title)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setStyleSheet("QScrollArea { background: transparent; border: none; }"
                              "QScrollArea > QWidget > QWidget { background: transparent; }")
        self.nuclide_panel_scroll = scroll

        content = QWidget()
        card_layout = QVBoxLayout(content)
        card_layout.setContentsMargins(2, 2, 2, 2)
        card_layout.setSpacing(10)

        self.nuclide_name_lbl = QLabel()
        self.nuclide_name_lbl.setTextFormat(Qt.RichText)
        card_layout.addWidget(self.nuclide_name_lbl)

        self.nuclide_facts_lbl = QLabel()
        self.nuclide_facts_lbl.setTextFormat(Qt.RichText)
        card_layout.addWidget(self.nuclide_facts_lbl)

        self.nuclide_divider1 = PanelDivider(is_dark=self.is_dark)
        card_layout.addWidget(self.nuclide_divider1)

        phys_title = QLabel("FİZİKSEL AÇIKLAMA"); phys_title.setObjectName("Engraved")
        card_layout.addWidget(phys_title)
        self.nuclide_phys_title_lbl = phys_title

        self.nuclide_phys_lbl = QLabel()
        self.nuclide_phys_lbl.setWordWrap(True)
        self.nuclide_phys_lbl.setTextFormat(Qt.RichText)
        card_layout.addWidget(self.nuclide_phys_lbl)

        self.nuclide_divider2 = PanelDivider(is_dark=self.is_dark)
        card_layout.addWidget(self.nuclide_divider2)

        gen_title = QLabel("GENEL BİLGİ"); gen_title.setObjectName("Engraved")
        card_layout.addWidget(gen_title)
        self.nuclide_gen_title_lbl = gen_title

        self.nuclide_general_lbl = QLabel()
        self.nuclide_general_lbl.setWordWrap(True)
        self.nuclide_general_lbl.setTextFormat(Qt.RichText)
        card_layout.addWidget(self.nuclide_general_lbl)

        card_layout.addStretch(1)
        scroll.setWidget(content)
        outer.addWidget(scroll, 1)
        return self.nuclide_panel_frame

    def handle_chart_selection(self):
        self._draw_segre_chart()
        iso = self.chart_combo.currentData() if hasattr(self, "chart_combo") else None
        if iso:
            self._update_nuclide_info_panel(iso)

    def _update_nuclide_info_panel(self, iso_key: str):
        if not hasattr(self, "nuclide_name_lbl"):
            return
        info = explain_isotope(iso_key)
        p = palette(self.is_dark)
        if not info:
            self.nuclide_name_lbl.setText(f"<h3 style='color:{p['accent']};margin:0;font-family:monospace;'>—</h3>")
            self.nuclide_facts_lbl.setText("")
            self.nuclide_phys_lbl.setText("")
            self.nuclide_general_lbl.setText("")
            return

        status_color = p["green"] if info["is_stable"] else p["red"]
        status_text = "KARARLI" if info["is_stable"] else f"KARARSIZ · {info['mode']}"

        self.nuclide_name_lbl.setText(f"""
        <h3 style='color:{p['accent']}; margin:0; font-family:monospace; letter-spacing:1px;'>{info['display_name']}</h3>
        <p style='color:{status_color}; font-weight:bold; margin:4px 0 0 0; font-family:monospace; font-size:11px;'>[●] {status_text}</p>
        """)

        self.nuclide_facts_lbl.setText(f"""
        <table style='width:100%; color:{p['text']}; font-family:monospace; font-size:11.5px; line-height:1.9;'>
            <tr><td style='color:{p['label']};'>Proton (Z):</td><td><b>{info['Z']}</b></td></tr>
            <tr><td style='color:{p['label']};'>Nötron (N):</td><td><b>{info['N']}</b></td></tr>
            <tr><td style='color:{p['label']};'>Kütle No (A):</td><td><b>{info['A']}</b></td></tr>
            <tr><td style='color:{p['label']};'>Yarı Ömür:</td><td><b>{info['half_life']}</b></td></tr>
        </table>
        """)

        self.nuclide_phys_lbl.setText(
            f"<p style='color:{p['text']}; font-family:monospace; font-size:10.5px; line-height:1.55;'>{info['physical_text']}</p>"
        )
        self.nuclide_general_lbl.setText(
            f"<p style='color:{p['text']}; font-family:monospace; font-size:10.5px; line-height:1.55;'>{info['general_text']}</p>"
        )

    def _draw_segre_chart(self):
        p = palette(self.is_dark)
        self.segre_canvas.apply_theme(self.is_dark)
        ax = self.segre_canvas.axes
        ax.clear()

        groups = get_segre_chart_data()
        style = {
            "stable":     (p["green"], "KARARLI",              10),
            "beta_minus": (p["blue"],  "β⁻ BOZUNUMU",           4),
            "beta_plus":  (p["accent"], "β⁺ / EC BOZUNUMU",     4),
            "alpha":      (p["red"],   "α BOZUNUMU (+ BİLEŞİK)", 4),
            "other":      (p["faint"], "DİĞER (SF, p, n, ...)",  4),
        }
        # Kararlılık vadisi en son çizilsin ki üstte, belirgin kalsın.
        for cat in ("other", "alpha", "beta_plus", "beta_minus", "stable"):
            xs, ys = groups[cat]
            color, label, size = style[cat]
            ax.scatter(xs, ys, s=size, color=color, alpha=0.85 if cat == "stable" else 0.55,
                       linewidths=0, label=label, zorder=3 if cat == "stable" else 2)

        ax.plot([0, 160], [0, 160], color=p["faint"], lw=0.8, ls=":", alpha=0.5, zorder=1, label="N = Z")

        iso = self.chart_combo.currentData() if hasattr(self, "chart_combo") else None
        if iso and iso in ISOTOPES:
            v = ISOTOPES[iso]
            N = v["A"] - v["Z"]
            ax.scatter([N], [v["Z"]], s=140, facecolors="none",
                       edgecolors=p["accent_hi"], linewidths=2.0, zorder=6)
            ax.annotate(get_iso_display_name(iso), xy=(N, v["Z"]), xytext=(8, 8),
                        textcoords="offset points", color=p["accent_hi"],
                        fontfamily="monospace", fontsize=8.5, zorder=7)

        ax.set_xlabel("NÖTRON SAYISI  N", color=p["label"], fontfamily="monospace", fontsize=9)
        ax.set_ylabel("PROTON SAYISI  Z", color=p["label"], fontfamily="monospace", fontsize=9)
        if getattr(self, "_segre_xlim", None) and getattr(self, "_segre_ylim", None):
            ax.set_xlim(self._segre_xlim)
            ax.set_ylim(self._segre_ylim)
        else:
            ax.set_xlim(-2, 160)
            ax.set_ylim(-2, 100)
        ax.set_facecolor(p["inset"])
        for s in ax.spines.values():
            s.set_color(p["frame"]); s.set_linewidth(1.1)
        ax.tick_params(colors=p["label"], labelsize=8.5, length=0)
        ax.grid(True, color=p["grid_minor"], lw=0.5, ls=":")
        leg = ax.legend(loc="lower right", fontsize=7.2, frameon=True, markerscale=1.8,
                        facecolor=p["chassis"] if self.is_dark else p["panel"],
                        edgecolor=p["frame"], labelcolor=p["text"])
        self.segre_canvas.fig.tight_layout()
        self.segre_canvas.draw()

    def create_moderation_tab(self):
        w = QWidget()
        main_layout = QVBoxLayout(w)
        main_layout.setContentsMargins(20, 20, 20, 20)
        main_layout.setSpacing(12)

        top_layout = QHBoxLayout()
        cap = QLabel("MODERATÖR MALZEME KANALI"); cap.setObjectName("Engraved")
        top_layout.addWidget(cap)

        self.mod_combo = PanelComboBox(is_dark=self.is_dark)
        for k, v in MODERATORS.items():
            self.mod_combo.addItem(f"{v['name']}  [MR ≈ {v['mr']:.0f}]", userData=k)
        self.mod_combo.currentIndexChanged.connect(self.update_moderation_display)
        top_layout.addWidget(self.mod_combo, 1)
        main_layout.addLayout(top_layout)

        content_layout = QHBoxLayout()
        content_layout.setSpacing(12)

        self.card_mod_left = QFrame()
        self.layout_mod_left = QVBoxLayout(self.card_mod_left)
        self.layout_mod_left.setSpacing(10)
        self.layout_mod_left.setContentsMargins(16, 16, 16, 16)

        self.mod_info_lbl = QLabel()
        self.mod_info_lbl.setTextFormat(Qt.RichText)
        self.layout_mod_left.addWidget(self.mod_info_lbl)

        self.divider_mod = PanelDivider(is_dark=self.is_dark)
        self.layout_mod_left.addWidget(self.divider_mod)

        self.mod_table_lbl = QLabel()
        self.mod_table_lbl.setTextFormat(Qt.RichText)
        self.layout_mod_left.addWidget(self.mod_table_lbl)
        self.layout_mod_left.addStretch(1)

        content_layout.addWidget(self.card_mod_left, 1)

        self.card_mod_right = QFrame()
        self.layout_mod_right = QVBoxLayout(self.card_mod_right)
        self.layout_mod_right.setSpacing(10)
        self.layout_mod_right.setContentsMargins(16, 16, 16, 16)

        cap_chart = QLabel("NÖTRON ENERJİ YAVAŞLAMA KİNETİĞİ   ·   E(k) = E₀ · e^(−k·ξ)")
        cap_chart.setObjectName("Engraved")
        self.layout_mod_right.addWidget(cap_chart)

        self.mod_canvas = MplCanvas()
        self.layout_mod_right.addWidget(self.mod_canvas, 1)

        content_layout.addWidget(self.card_mod_right, 1)
        main_layout.addLayout(content_layout)

        exp_row = QHBoxLayout()
        exp_row.addStretch(1)
        exp_row.addWidget(self._make_export_button(
            lambda: self._export_report_pdf(self._build_moderation_report_text, "moderasyon_raporu.pdf")))
        main_layout.addLayout(exp_row)

        self.update_moderation_display()
        return w

    def _build_moderation_report_text(self) -> str:
        mod_k = self.mod_combo.currentData() if hasattr(self, "mod_combo") else "H2O"
        if not mod_k: mod_k = "H2O"
        res = calculate_moderation(mod_k)
        d = res["data"]

        rap  = " SİNOP ÜNİVERSİTESİ · NÜKLEER ENERJİ MÜHENDİSLİĞİ\n"
        rap += " NÖTRON MODERASYON ÖLÇÜM FİŞİ" + " " * 12 + "FORM SNU-4/56\n"
        rap += "=" * (W_FIS + 2) + "\n"
        rap += _fis("MODERATÖR", d["name"])
        rap += _fis("REAKTÖR KULLANIMI", d["reactor_type"])
        rap += _fis("EFEKTİF KÜTLE NO  A", f"{d['A_eff']:.1f}")
        rap += "-" * (W_FIS + 2) + "\n"
        rap += _fis("KİNEMATİK ÇARPAN  α", f"{res['alpha']:12.4f}")
        rap += _fis("ORT. LOG ENERJİ KAYBI  ξ", f"{d['xi']:12.3f}")
        rap += _fis("TEK VURUŞTA MAX KAYIP", f"{res['max_loss']:12.1f}", "%")
        rap += _fis("LOG-ORT. ENERJİ KAYBI", f"{res['avg_loss']:12.1f}", "% / çarpışma")
        rap += _fis("ARİTMETİK ORT. KAYIP", f"{res['avg_loss_arith']:12.1f}", "% / çarpışma")
        rap += _fis("TOPLAM LETARJİ ARTIŞI  Δu", f"{res['delta_u']:12.2f}", "(2 MeV ➔ 0.025 eV)")
        rap += "-" * (W_FIS + 2) + "\n"
        rap += _fis("TERMALLEŞME ÇARPIŞMA SAYISI", f"{res['n_coll']:12.1f}", "çarpışma")
        rap += _fis("MAKROSKOPİK YAVAŞLAMA", f"{d['msdp']:12.3f}", "cm⁻¹")
        rap += _fis("MAKROSKOPİK YUTULMA  Σa", f"{d['Sigma_a']:12.6f}", "cm⁻¹")
        rap += _fis("MODERASYON ORANI  MR", f"{d['mr']:12.1f}")
        rap += "\n" + "=" * (W_FIS + 2) + "\n"
        rap += f" ÖLÇÜM {QDateTime.currentDateTime().toString('dd.MM.yyyy HH:mm:ss')}"
        rap += "   ·   OPERATÖR: A.S.ÇETİN\n"
        return rap

    def update_moderation_display(self):
        mod_k = self.mod_combo.currentData() if hasattr(self, "mod_combo") else "H2O"
        if not mod_k: mod_k = "H2O"
        res = calculate_moderation(mod_k)
        d = res["data"]
        p = palette(self.is_dark)

        self.card_mod_left.setStyleSheet(panel_qss(self.is_dark))
        self.card_mod_right.setStyleSheet(panel_qss(self.is_dark))

        self.mod_info_lbl.setText(f"""
        <h3 style='color: {p['accent']}; margin: 0; font-family: monospace; letter-spacing: 1px;'>{d['name']}</h3>
        <p style='color: {p['green']}; font-weight: bold; margin: 4px 0 0 0; font-family: monospace; font-size: 11px;'>[●] REAKTÖR KULLANIMI: {d['reactor_type']}</p>
        """)

        self.mod_table_lbl.setText(f"""
        <table style='width: 100%; color: {p['text']}; font-family: monospace; font-size: 12px; line-height: 2.1;'>
            <tr><td style='color: {p['label']};'>Efektif Kütle No (A):</td><td><b>{d['A_eff']:.1f}</b></td></tr>
            <tr><td style='color: {p['label']};'>Kinematik Çarpan (α):</td><td><b>{res['alpha']:.4f}</b></td></tr>
            <tr><td style='color: {p['accent']};'>Ort. Log Enerji Kaybı (ξ):</td><td><b>{d['xi']:.3f}</b></td></tr>
            <tr><td style='color: {p['label']};'>Tek Vuruşta Max Kayıp:</td><td><b>%{res['max_loss']:.1f}</b></td></tr>
            <tr><td style='color: {p['label']};'>Log-Ort. Enerji Kaybı (1−e^−ξ):</td><td><b>%{res['avg_loss']:.1f} / çarpışma</b></td></tr>
            <tr><td style='color: {p['label']};'>Aritmetik Ort. Kayıp ((1−α)/2):</td><td><b>%{res['avg_loss_arith']:.1f} / çarpışma</b></td></tr>
            <tr><td style='color: {p['label']};'>Toplam Letarji Artışı (Δu):</td><td><b>{res['delta_u']:.2f}</b> (2 MeV ➔ 0.025 eV)</td></tr>
            <tr><td style='color: {p['green']}; font-weight:bold;'>Termalleşme Çarpışma Sayısı (n):</td><td><b>~{res['n_coll']:.1f} çarpışma</b></td></tr>
            <tr><td style='color: {p['label']};'>Makroskopik Yavaşlama (MSDP):</td><td><b>{d['msdp']:.3f} cm⁻¹</b></td></tr>
            <tr><td style='color: {p['label']};'>Makroskopik Yutulma (Σa):</td><td><b>{d['Sigma_a']:.6f} cm⁻¹</b></td></tr>
            <tr><td style='color: {p['accent']}; font-weight:bold;'>Moderasyon Oranı (MR):</td><td><b>{d['mr']:.1f}</b></td></tr>
        </table>
        """)

        self.mod_canvas.apply_theme(self.is_dark)
        ax = self.mod_canvas.axes
        ax.clear()

        colors_ref = {"H2O": "#ffb454", "D2O": "#4fa8b2", "Be": "#86c98f", "C": "#e2614b"}
        for k_other, d_other in MODERATORS.items():
            if k_other != mod_k:
                res_other = calculate_moderation(k_other)
                ax.plot(res_other["k_arr"], res_other["E_arr"], color=colors_ref.get(k_other, p["label"]),
                        lw=1.2, ls="--", alpha=0.55, label=f"{k_other} (n={res_other['n_coll']:.0f})")

        self.mod_canvas.trace(res["k_arr"], res["E_arr"], d["color"],
                              f"{mod_k} SEÇİLİ (n={res['n_coll']:.0f})", lw=2.4)

        ax.axhline(0.0253, color=p["red"], lw=1.3, ls=":", label="Termal Eşik (0.0253 eV)")

        ax.set_yscale("log")
        ax.set_xlim(0, 130)
        ax.set_ylim(0.01, 3.0e6)

        ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=8))
        ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=(0.2, 0.4, 0.6, 0.8), numticks=8))
        ax.yaxis.set_minor_formatter(NullFormatter())

        ax.set_xlabel("ÇARPIŞMA SAYISI  [k]", color=p["label"], fontfamily="monospace", fontsize=9)
        ax.set_ylabel("NÖTRON ENERJİSİ  [eV]", color=p["label"], fontfamily="monospace", fontsize=9)

        ax.grid(True, which="major", color=p["grid_major"], lw=0.8, ls="-")
        ax.grid(True, which="minor", color=p["grid_minor"], lw=0.4, ls=":")

        for s in ax.spines.values():
            s.set_color(p["frame"])
            s.set_linewidth(1.1)

        ax.tick_params(which="both", direction="in", top=True, right=True, colors=p["label"], labelsize=8.5)
        for lbl in ax.get_xticklabels() + ax.get_yticklabels():
            lbl.set_fontfamily("monospace")

        self.mod_canvas.style_legend()
        self.mod_canvas.readout_box(f"E₀ = 2.0 MeV\nEth = 0.025 eV\nξ = {d['xi']:.3f}")
        self.mod_canvas.fig.tight_layout()
        self.mod_canvas.draw()

    def create_density_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(10)

        self.dens_frame = QFrame()
        self.dens_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        f_layout = QVBoxLayout(self.dens_frame)
        f_layout.setContentsMargins(12, 10, 12, 12)
        f_layout.setSpacing(10)

        title = QLabel("MALZEME KANALI   ·   ATOM YOĞUNLUĞU N = ρ·N_A / M")
        title.setObjectName("Engraved")
        f_layout.addWidget(title)

        top_row = QHBoxLayout(); top_row.setSpacing(10)

        vbox_mat = QVBoxLayout()
        lbl_mat = QLabel("MALZEME"); lbl_mat.setObjectName("Engraved")
        vbox_mat.addWidget(lbl_mat)
        self.dens_combo = PanelComboBox(is_dark=self.is_dark)
        for k, v in MATERIALS.items():
            self.dens_combo.addItem(v["name"], userData=k)
        self.dens_combo.currentIndexChanged.connect(self.handle_material_changed)
        vbox_mat.addWidget(self.dens_combo)
        top_row.addLayout(vbox_mat, 3)

        vbox_formula = QVBoxLayout()
        lbl_formula = QLabel("KİMYASAL FORMÜL"); lbl_formula.setObjectName("Engraved")
        vbox_formula.addWidget(lbl_formula)
        self.dens_formula = QLineEdit("U")
        self.dens_formula.editingFinished.connect(self.handle_density_calc)
        vbox_formula.addWidget(self.dens_formula)
        top_row.addLayout(vbox_formula, 2)

        vbox_rho = QVBoxLayout()
        lbl_rho = QLabel("YOĞUNLUK  ρ  [g/cm³]"); lbl_rho.setObjectName("Engraved")
        vbox_rho.addWidget(lbl_rho)
        self.dens_rho = QLineEdit("19.10")
        self.dens_rho.editingFinished.connect(self.handle_density_calc)
        vbox_rho.addWidget(self.dens_rho)
        top_row.addLayout(vbox_rho, 2)

        vbox_enr = QVBoxLayout()
        self.lbl_enr = QLabel("ZENGİNLEŞTİRME  w(²³⁵U)  [ağırlıkça %]")
        self.lbl_enr.setObjectName("Engraved")
        vbox_enr.addWidget(self.lbl_enr)
        self.dens_enrich = QLineEdit("0.711")
        self.dens_enrich.editingFinished.connect(self.handle_density_calc)
        vbox_enr.addWidget(self.dens_enrich)
        top_row.addLayout(vbox_enr, 2)

        f_layout.addLayout(top_row)
        layout.addWidget(self.dens_frame)

        btn_row = QHBoxLayout()
        calc_btn = QPushButton("HESAPLA  ·  ATOM YOĞUNLUĞU")
        calc_btn.setObjectName("PrimaryLever")
        calc_btn.clicked.connect(self.handle_density_calc)
        btn_row.addWidget(calc_btn)
        btn_row.addStretch(1)
        btn_row.addWidget(self._make_export_button(
            lambda: self._export_report_pdf(self.dens_output.toPlainText, "yogunluk_raporu.pdf")))
        layout.addLayout(btn_row)

        result_row = QHBoxLayout(); result_row.setSpacing(10)
        self.dens_output = QTextEdit()
        self.dens_output.setReadOnly(True)
        result_row.addWidget(self.dens_output, 1)

        self.dens_canvas = MplCanvas()
        result_row.addWidget(self.dens_canvas, 1)
        layout.addLayout(result_row, 1)

        self.handle_material_changed()
        return w

    def handle_material_changed(self):
        key = self.dens_combo.currentData()
        mat = MATERIALS.get(key, MATERIALS["U_METAL"])
        self.dens_formula.setText(mat["formula"])
        self.dens_rho.setText(f"{mat['rho']:.3f}")
        self.handle_density_calc()

    def _draw_density_chart(self, r):
        """Formüldeki her element/izotop için atom yoğunluğunu (log ölçek)
        yatay çubuk grafiği olarak çizer. U varsa ²³⁵U/²³⁸U ayrı çubuklar olur."""
        p = palette(self.is_dark)
        self.dens_canvas.apply_theme(self.is_dark)
        ax = self.dens_canvas.axes
        ax.clear()

        labels, values, colors = [], [], []
        for key, e in r["elements"].items():
            if key == "U" and "N_235" in e:
                labels.append("²³⁵U"); values.append(e["N_235"]); colors.append(p["accent"])
                labels.append("²³⁸U"); values.append(e["N_238"]); colors.append(p["green"])
            else:
                labels.append(key); values.append(e["N"]); colors.append(p["accent"])

        y_pos = list(range(len(labels)))
        ax.barh(y_pos, values, color=colors, height=0.55, zorder=3)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontfamily="monospace", fontsize=9.5)
        ax.invert_yaxis()
        ax.set_xscale("log")
        ax.set_xlabel("ATOM YOĞUNLUĞU  N  [atom/cm³]", color=p["label"],
                       fontfamily="monospace", fontsize=9)

        for y, v in zip(y_pos, values):
            ax.text(v, y, f"  {v:.3e}", va="center", ha="left",
                    color=p["readout"], fontfamily="monospace", fontsize=8)

        ax.set_facecolor(p["inset"])
        for s in ax.spines.values():
            s.set_color(p["frame"]); s.set_linewidth(1.1)
        ax.tick_params(colors=p["label"], labelsize=8.5, length=0)
        ax.grid(True, axis="x", which="both", color=p["grid_minor"], lw=0.5, ls=":")
        self.dens_canvas.fig.tight_layout()
        self.dens_canvas.draw()

    def _formula_has_uranium(self, formula: str) -> bool:
        try:
            return any(sym == "U" for sym, _ in parse_formula(formula))
        except Exception:
            return False

    def handle_density_calc(self):
        formula = self.dens_formula.text().strip()
        has_u = self._formula_has_uranium(formula)
        self.lbl_enr.setVisible(has_u)
        self.dens_enrich.setVisible(has_u)

        try:
            rho = float(self.dens_rho.text().replace(",", "."))
        except ValueError:
            self.dens_output.setPlainText(" ARIZA · YOĞUNLUK (ρ) SAYISAL DEĞİL")
            self._clear_density_chart()
            return

        w235 = None
        if has_u:
            try:
                w235 = float(self.dens_enrich.text().replace(",", "."))
            except ValueError:
                self.dens_output.setPlainText(" ARIZA · ZENGİNLEŞTİRME YÜZDESİ SAYISAL DEĞİL")
                self._clear_density_chart()
                return

        try:
            r = calculate_material(rho, formula, w235)
        except Exception as e:
            self.dens_output.setPlainText(f" ARIZA · FORMÜL ÇÖZÜMLENEMEDİ\n {e}")
            self._clear_density_chart()
            return

        rap  = " SİNOP ÜNİVERSİTESİ · NÜKLEER ENERJİ MÜHENDİSLİĞİ\n"
        rap += " ATOM YOĞUNLUĞU ÖLÇÜM FİŞİ" + " " * 15 + "FORM SNU-7/56\n"
        rap += "=" * (W_FIS + 2) + "\n"
        rap += _fis("FORMÜL", formula)
        rap += _fis("YOĞUNLUK  ρ", f"{rho:12.4f}", "g/cm³")
        rap += _fis("FORMÜL AĞIRLIĞI  M", f"{r['M_formula']:12.4f}", "g/mol")
        rap += _fis("MOLEKÜL/FORMÜL YOĞUNLUĞU", f"{r['N_molecule']:12.4e}", "molekül/cm³")
        rap += "-" * (W_FIS + 2) + "\n"

        for key, e in r["elements"].items():
            rap += _fis(f"{key}  (adet/formül={e['count']})  M={e['M']:.3f}",
                        f"{e['N']:12.4e}", "atom/cm³")
            if key == "U" and "N_235" in e:
                mix = r["mix"]
                rap += _fis(f"  ↳ ²³⁵U  ({mix['atom_pct_235']:.4f} atom%)",
                            f"{e['N_235']:12.4e}", "atom/cm³")
                rap += _fis(f"  ↳ ²³⁸U  ({mix['atom_pct_238']:.4f} atom%)",
                            f"{e['N_238']:12.4e}", "atom/cm³")

        rap += "\n" + "=" * (W_FIS + 2) + "\n"
        rap += (" NOT: Atomik ağırlıklar elle girilmez; doğal izotopik bolluğa\n"
                " (NUBASE2020) göre ağırlıklı ortalama olarak hesaplanır.\n")
        rap += f" ÖLÇÜM {QDateTime.currentDateTime().toString('dd.MM.yyyy HH:mm:ss')}"
        rap += "   ·   OPERATÖR: A.S.ÇETİN\n"
        self.dens_output.setPlainText(rap)
        self._draw_density_chart(r)

    def _clear_density_chart(self):
        self.dens_canvas.apply_theme(self.is_dark)
        self.dens_canvas.axes.clear()
        self.dens_canvas.draw()

    def create_shielding_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(10)

        self.shield_frame = QFrame()
        self.shield_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        f_layout = QVBoxLayout(self.shield_frame)
        f_layout.setContentsMargins(12, 10, 12, 12)
        f_layout.setSpacing(10)

        title = QLabel("ZIRHLAMA KANALI   ·   GAMA ZAYIFLAMASI  I = I₀ · e^(−μx)  ·  BUILDUP  B(μx)")
        title.setObjectName("Engraved")
        f_layout.addWidget(title)

        top_row = QHBoxLayout(); top_row.setSpacing(10)

        vbox_src = QVBoxLayout()
        lbl_src = QLabel("GELEN IŞIN  (KAYNAK)"); lbl_src.setObjectName("Engraved")
        vbox_src.addWidget(lbl_src)
        self.shield_source_combo = PanelComboBox(is_dark=self.is_dark)
        self.shield_source_combo.addItem("Serbest Enerji (MeV Gir)", userData="FREE")
        for k, v in RADIOACTIVE_SOURCES.items():
            n_lines = len(v["lines"])
            self.shield_source_combo.addItem(f"{v['name']}  [{n_lines} çizgi]", userData=k)
        self.shield_source_combo.currentIndexChanged.connect(self.handle_shielding_source_changed)
        vbox_src.addWidget(self.shield_source_combo)
        top_row.addLayout(vbox_src, 3)

        vbox_mat = QVBoxLayout()
        lbl_mat = QLabel("ZIRH MALZEMESİ"); lbl_mat.setObjectName("Engraved")
        vbox_mat.addWidget(lbl_mat)
        self.shield_combo = PanelComboBox(is_dark=self.is_dark)
        for k, v in SHIELDING_MATERIALS.items():
            self.shield_combo.addItem(f"{v['name']}  [ρ={v['density']:.2f} g/cm³]", userData=k)
        self.shield_combo.currentIndexChanged.connect(self.handle_shielding_calc)
        vbox_mat.addWidget(self.shield_combo)
        top_row.addLayout(vbox_mat, 3)

        vbox_e = QVBoxLayout()
        lbl_e = QLabel("FOTON ENERJİSİ  [MeV]"); lbl_e.setObjectName("Engraved")
        vbox_e.addWidget(lbl_e)
        self.shield_energy = QLineEdit("1.0")
        self.shield_energy.editingFinished.connect(self.handle_shielding_calc)
        vbox_e.addWidget(self.shield_energy)
        top_row.addLayout(vbox_e, 2)

        vbox_x = QVBoxLayout()
        lbl_x = QLabel("ZIRH KALINLIĞI  x  [cm]"); lbl_x.setObjectName("Engraved")
        vbox_x.addWidget(lbl_x)
        self.shield_thickness = QLineEdit("5.0")
        self.shield_thickness.editingFinished.connect(self.handle_shielding_calc)
        vbox_x.addWidget(self.shield_thickness)
        top_row.addLayout(vbox_x, 2)

        f_layout.addLayout(top_row)

        self.shield_source_info = QLabel()
        self.shield_source_info.setObjectName("Engraved")
        self.shield_source_info.setWordWrap(True)
        self.shield_source_info.setVisible(False)
        f_layout.addWidget(self.shield_source_info)

        layout.addWidget(self.shield_frame)

        btn_row = QHBoxLayout()
        calc_btn = QPushButton("HESAPLA  ·  ZAYIFLAMA")
        calc_btn.setObjectName("PrimaryLever")
        calc_btn.clicked.connect(self.handle_shielding_calc)
        btn_row.addWidget(calc_btn)
        self.shield_lamp = IndicatorLamp("BUILDUP VERİSİ MEVCUT", is_dark=self.is_dark)
        btn_row.addSpacing(16)
        btn_row.addWidget(self.shield_lamp)
        btn_row.addStretch(1)
        btn_row.addWidget(self._make_export_button(
            lambda: self._export_report_pdf(self.shield_output.toPlainText, "zirhlama_raporu.pdf")))
        layout.addLayout(btn_row)

        result_row = QHBoxLayout(); result_row.setSpacing(10)
        self.shield_output = QTextEdit()
        self.shield_output.setReadOnly(True)
        result_row.addWidget(self.shield_output, 1)

        self.shield_canvas = MplCanvas()
        result_row.addWidget(self.shield_canvas, 1)
        layout.addLayout(result_row, 1)

        self.handle_shielding_source_changed()
        return w

    def handle_shielding_source_changed(self):
        source_key = self.shield_source_combo.currentData() or "FREE"
        is_free = (source_key == "FREE")
        self.shield_energy.setEnabled(is_free)

        if is_free:
            self.shield_source_info.setVisible(False)
        else:
            src = RADIOACTIVE_SOURCES[source_key]
            parts = []
            for l in src["lines"]:
                e_disp = f"{l['energy']*1000:.1f} keV" if l["energy"] < 1.0 else f"{l['energy']:.3f} MeV"
                parts.append(f"{e_disp} (%{l['intensity']*100:.2f})")
            self.shield_source_info.setText(
                f"KAYNAK ÇİZGİLERİ (IAEA referans şiddetleri): {'   ·   '.join(parts)}"
            )
            self.shield_source_info.setVisible(True)

        self.handle_shielding_calc()

    def handle_shielding_calc(self):
        source_key = self.shield_source_combo.currentData() or "FREE" if hasattr(self, "shield_source_combo") else "FREE"
        material_key = self.shield_combo.currentData() or "Pb"

        try:
            thickness = float(self.shield_thickness.text().replace(",", "."))
            if thickness < 0:
                raise ValueError
        except ValueError:
            self.shield_output.setPlainText(" ARIZA · KALINLIK DEĞERİ GEÇERSİZ (negatif olamaz)")
            self._clear_shielding_chart()
            self.shield_lamp.set_state("idle")
            return

        if source_key == "FREE":
            try:
                energy = float(self.shield_energy.text().replace(",", "."))
                if energy <= 0:
                    raise ValueError
            except ValueError:
                self.shield_output.setPlainText(" ARIZA · FOTON ENERJİSİ GEÇERSİZ (0'dan büyük olmalı)")
                self._clear_shielding_chart()
                self.shield_lamp.set_state("idle")
                return
            self._calc_shielding_single(material_key, energy, thickness)
        else:
            self._calc_shielding_multiline(material_key, source_key, thickness)

    def _calc_shielding_single(self, material_key, energy, thickness):
        r = gamma_attenuation(material_key, energy, thickness, I0=100.0)
        x10 = gamma_thickness_for_factor(material_key, energy, 10.0)
        x100 = gamma_thickness_for_factor(material_key, energy, 100.0)
        x1000 = gamma_thickness_for_factor(material_key, energy, 1000.0)

        rap  = " SİNOP ÜNİVERSİTESİ · NÜKLEER ENERJİ MÜHENDİSLİĞİ\n"
        rap += " ZIRHLAMA / GAMA ZAYIFLAMASI ÖLÇÜM FİŞİ" + " " * 3 + "FORM SNU-9/56\n"
        rap += "=" * (W_FIS + 2) + "\n"
        rap += _fis("GELEN IŞIN", "Serbest enerji (kaynak seçilmedi)")
        rap += _fis("ZIRH MALZEMESİ", r["material"])
        rap += _fis("YOĞUNLUK  ρ", f"{r['rho']:12.4f}", "g/cm³")
        rap += _fis("FOTON ENERJİSİ", f"{energy:12.4f}", "MeV")
        rap += _fis("ZIRH KALINLIĞI  x", f"{thickness:12.4f}", "cm")
        rap += "-" * (W_FIS + 2) + "\n"
        rap += _fis("KÜTLE ZAYIFLAMA KATSAYISI  μ/ρ", f"{r['mu_rho']:12.5f}", "cm²/g")
        rap += _fis("LİNEER ZAYIFLAMA KATSAYISI  μ", f"{r['mu']:12.5f}", "cm⁻¹")
        rap += _fis("KALINLIK  (MFP CİNSİNDEN)  μx", f"{r['mu_x']:12.4f}", "mfp")
        rap += "-" * (W_FIS + 2) + "\n"
        rap += " DAR-DEMET (NARROW-BEAM, SAÇILAN FOTON YOK)\n"
        rap += _fis("GİRİŞ ŞİDDETİ  I₀", f"{r['I0']:12.4f}", "(keyfi birim)")
        rap += _fis("ÇIKIŞ ŞİDDETİ  I", f"{r['I']:12.4f}", "(keyfi birim)")
        rap += _fis("GEÇİRGENLİK  I/I₀", f"{r['transmission']*100:12.4f}", "%")
        rap += _fis("ZAYIFLAMA", f"{r['attenuation_pct']:12.4f}", "%")
        rap += "-" * (W_FIS + 2) + "\n"
        rap += self._buildup_report_block(r)
        rap += "-" * (W_FIS + 2) + "\n"
        rap += _fis("YARI-DEĞER KALINLIĞI  HVL", f"{r['hvl']:12.4f}", "cm")
        rap += _fis("ONDA-BİR DEĞER KALINLIĞI  TVL", f"{r['tvl']:12.4f}", "cm")
        rap += _fis("SEÇİLEN KALINLIK  (HVL CİNSİNDEN)", f"{r['n_hvl']:12.3f}", "× HVL")
        rap += "-" * (W_FIS + 2) + "\n"
        rap += " GEREKLİ KALINLIK (DAR-DEMET, AZALTMA FAKTÖRÜNE GÖRE)\n"
        if x10 is not None:
            rap += _fis("  10× AZALTMA İÇİN", f"{x10:12.4f}", "cm")
        if x100 is not None:
            rap += _fis("  100× AZALTMA İÇİN", f"{x100:12.4f}", "cm")
        if x1000 is not None:
            rap += _fis("  1000× AZALTMA İÇİN", f"{x1000:12.4f}", "cm")
        rap += "\n" + "=" * (W_FIS + 2) + "\n"
        rap += self._shielding_report_footer()
        self.shield_output.setPlainText(rap)
        self.shield_lamp.set_state("ok" if r["B"] is not None else "fail")
        self._draw_shielding_chart_single(r)

    def _calc_shielding_multiline(self, material_key, source_key, thickness):
        m = gamma_attenuation_multiline(material_key, source_key, thickness, I0_total=100.0)

        rap  = " SİNOP ÜNİVERSİTESİ · NÜKLEER ENERJİ MÜHENDİSLİĞİ\n"
        rap += " ZIRHLAMA / GAMA ZAYIFLAMASI ÖLÇÜM FİŞİ" + " " * 3 + "FORM SNU-9/56\n"
        rap += "=" * (W_FIS + 2) + "\n"
        rap += _fis("GELEN IŞIN  (KAYNAK)", m["source"])
        rap += _fis("ZIRH MALZEMESİ", m["material"])
        rap += _fis("ZIRH KALINLIĞI  x", f"{thickness:12.4f}", "cm")
        rap += "-" * (W_FIS + 2) + "\n"
        rap += " ÇİZGİ BAZLI ZAYIFLAMA (HER ENERJİ AYRI μ İLE HESAPLANIR)\n"
        rap += "-" * (W_FIS + 2) + "\n"
        for l in m["lines"]:
            e_disp = f"{l['energy']*1000:.1f} keV" if l["energy"] < 1.0 else f"{l['energy']:.4f} MeV"
            rap += _fis(f"  {e_disp}  (şiddet payı %{l['intensity_frac']*100:.2f})",
                        f"{l['I0']:10.3f} → {l['I']:.4f}", "(keyfi birim)")
            rap += _fis(f"    μ = {l['mu']:.5f} cm⁻¹ , μx = {l['mu_x']:.3f} mfp",
                        f"HVL={l['hvl']:.3f} cm")
        rap += "-" * (W_FIS + 2) + "\n"
        rap += " BİRLEŞİK (AĞIRLIKLI TOPLAM) SONUÇ — DAR-DEMET\n"
        rap += _fis("TOPLAM GİRİŞ ŞİDDETİ  I₀", f"{m['I0']:12.4f}", "(keyfi birim)")
        rap += _fis("TOPLAM ÇIKIŞ ŞİDDETİ  I", f"{m['I']:12.4f}", "(keyfi birim)")
        rap += _fis("BİRLEŞİK GEÇİRGENLİK  I/I₀", f"{m['transmission']*100:12.4f}", "%")
        rap += _fis("BİRLEŞİK ZAYIFLAMA", f"{m['attenuation_pct']:12.4f}", "%")
        rap += "-" * (W_FIS + 2) + "\n"
        if m["has_buildup"]:
            rap += _fis("GENİŞ-DEMET ÇIKIŞ ŞİDDETİ  I·B", f"{m['I_broad']:12.4f}", "(keyfi birim)")
            rap += _fis("GENİŞ-DEMET GEÇİRGENLİK", f"{m['transmission_broad']*100:12.4f}", "%")
        else:
            rap += " NOT: Bu kaynağın çizgilerinden en az biri doğrulanmış buildup\n"
            rap += " tablosu aralığı (1-10 MeV, Su/Beton/Demir) dışında kaldığından\n"
            rap += " birleşik geniş-demet sonucu hesaplanmadı — sonuç dar-demettir.\n"
        rap += "\n" + "=" * (W_FIS + 2) + "\n"
        rap += self._shielding_report_footer()
        self.shield_output.setPlainText(rap)
        self.shield_lamp.set_state("ok" if m["has_buildup"] else "fail")
        self._draw_shielding_chart_multi(m)

    def _buildup_report_block(self, r) -> str:
        rap = " GENİŞ-DEMET (BROAD-BEAM, TAYLOR BUILDUP FAKTÖRÜ İLE)\n"
        if r["B"] is not None:
            rap += _fis("BUILDUP FAKTÖRÜ  B(μx)", f"{r['B']:12.4f}")
            rap += _fis("GENİŞ-DEMET ÇIKIŞ ŞİDDETİ  I·B", f"{r['I_broad']:12.4f}", "(keyfi birim)")
            rap += _fis("GENİŞ-DEMET GEÇİRGENLİK", f"{r['transmission_broad']*100:12.4f}", "%")
        else:
            rap += " NOT: Bu enerji/malzeme kombinasyonu doğrulanmış buildup tablosu\n"
            rap += " aralığının (1-10 MeV; yalnız Su/Beton/Demir) dışında kaldığından\n"
            rap += " buildup hesaplanmadı — sonuç dar-demettir (gerçek doz bundan\n"
            rap += " YÜKSEK olabilir).\n"
        return rap

    def _shielding_report_footer(self) -> str:
        return (
            " NOT: μ/ρ değerleri NIST X-Ray/Gamma-Ray Mass Attenuation Coefficients\n"
            " referans tablolarından (log-log enterpolasyon) alınmıştır. Buildup\n"
            " faktörü B(μx), Taylor'ın iki-üstel formuyle (A1,α1,α2 katsayıları\n"
            " Univ. of Illinois NPRE 441 ders notlarından doğrulanmıştır) hesaplanır\n"
            " ve yalnız Su/Beton/Demir + 1-10 MeV aralığında geçerlidir. Kaynak gama\n"
            " çizgileri IAEA referans şiddet veritabanındandır (≥%1 şiddetli çizgiler).\n"
            f" ÖLÇÜM {QDateTime.currentDateTime().toString('dd.MM.yyyy HH:mm:ss')}"
            "   ·   OPERATÖR: A.S.ÇETİN\n"
        )

    def _draw_shielding_chart_single(self, r):
        p = palette(self.is_dark)
        self.shield_canvas.apply_theme(self.is_dark)
        ax = self.shield_canvas.axes
        ax.clear()

        mu, hvl = r["mu"], r["hvl"]
        xmax = max(r["thickness"] * 1.3, (hvl * 6 if hvl not in (0.0, float("inf")) else 1.0), 1.0)
        xs = np.linspace(0, xmax, 220)
        ys_narrow = r["I0"] * np.exp(-mu * xs)

        ax.plot(xs, ys_narrow, color=p["accent"], lw=2.2, zorder=3, label="Dar-demet  I₀·e^(−μx)")

        if r["B"] is not None:
            ys_broad = np.array([r["I0"] * math.exp(-mu * x) *
                                 (buildup_factor(r["material_key"], r["energy"], mu * x) or 1.0)
                                 for x in xs])
            ax.plot(xs, ys_broad, color=p["red"], lw=1.8, ls="--", zorder=3,
                    label="Geniş-demet  I₀·e^(−μx)·B(μx)")
            ax.scatter([r["thickness"]], [max(r["I_broad"], 1e-6)], color=p["red"], s=40, zorder=5)

        ax.axvline(r["thickness"], color=p["faint"], lw=1.3, ls="--", zorder=2)
        ax.scatter([r["thickness"]], [max(r["I"], 1e-6)], color=p["accent"], s=55, zorder=4)
        ax.annotate(f"x={r['thickness']:.2f} cm", (r["thickness"], max(r["I"], 1e-6)),
                    textcoords="offset points", xytext=(8, 10), color=p["text"],
                    fontfamily="monospace", fontsize=8, zorder=5)

        ax.set_yscale("log")
        ax.set_ylim(max(r["I0"] * 1e-4, 1e-3), r["I0"] * 12.0)
        ax.set_xlim(0, xmax)
        ax.set_xlabel("ZIRH KALINLIĞI  x  [cm]", color=p["label"], fontfamily="monospace", fontsize=9)
        ax.set_ylabel("ŞİDDET  I  (I₀=100)", color=p["label"], fontfamily="monospace", fontsize=9)
        ax.set_facecolor(p["inset"])
        for s in ax.spines.values():
            s.set_color(p["frame"]); s.set_linewidth(1.1)
        ax.tick_params(colors=p["label"], labelsize=8.5, length=0)
        ax.grid(True, which="major", color=p["grid_minor"], lw=0.5, ls=":")
        self.shield_canvas.style_legend()

        b_txt = f"B = {r['B']:.3f}" if r["B"] is not None else "B = n/a"
        self.shield_canvas.readout_box(f"μ = {mu:.4f} cm⁻¹\nHVL = {hvl:.3f} cm\n{b_txt}")
        self.shield_canvas.fig.tight_layout()
        self.shield_canvas.draw()

    def _draw_shielding_chart_multi(self, m):
        p = palette(self.is_dark)
        self.shield_canvas.apply_theme(self.is_dark)
        ax = self.shield_canvas.axes
        ax.clear()

        line_colors = [p["accent"], p["red"], p["blue"], p["green"]]
        max_hvl = max((l["hvl"] for l in m["lines"] if l["hvl"] not in (0.0, float("inf"))), default=1.0)
        xmax = max(m["thickness"] * 1.3, max_hvl * 6, 1.0)
        xs = np.linspace(0, xmax, 220)

        total_narrow = np.zeros_like(xs)
        total_broad = np.zeros_like(xs)
        for i, l in enumerate(m["lines"]):
            ys = l["I0"] * np.exp(-l["mu"] * xs)
            total_narrow += ys
            if l["B"] is not None:
                ys_b = np.array([l["I0"] * math.exp(-l["mu"] * x) *
                                 (buildup_factor(m["material_key"], l["energy"], l["mu"] * x) or 1.0)
                                 for x in xs])
            else:
                ys_b = ys
            total_broad += ys_b
            e_disp = f"{l['energy']*1000:.0f} keV" if l["energy"] < 1.0 else f"{l['energy']:.3f} MeV"
            ax.plot(xs, ys, color=line_colors[i % len(line_colors)], lw=1.0, ls=":", alpha=0.65,
                    zorder=2, label=f"{e_disp} (dar-demet)")

        ax.plot(xs, total_narrow, color=p["text"], lw=2.4, zorder=3, label="Birleşik (dar-demet)")
        if m["has_buildup"]:
            ax.plot(xs, total_broad, color=p["red"], lw=1.8, ls="--", zorder=3,
                    label="Birleşik (geniş-demet, B dahil)")

        ax.axvline(m["thickness"], color=p["faint"], lw=1.3, ls="--", zorder=1)
        ax.scatter([m["thickness"]], [max(m["I"], 1e-6)], color=p["text"], s=55, zorder=4)

        ax.set_yscale("log")
        ax.set_ylim(max(m["I0"] * 1e-4, 1e-3), m["I0"] * 12.0)
        ax.set_xlim(0, xmax)
        ax.set_xlabel("ZIRH KALINLIĞI  x  [cm]", color=p["label"], fontfamily="monospace", fontsize=9)
        ax.set_ylabel("ŞİDDET  I  (I₀=100)", color=p["label"], fontfamily="monospace", fontsize=9)
        ax.set_facecolor(p["inset"])
        for s in ax.spines.values():
            s.set_color(p["frame"]); s.set_linewidth(1.1)
        ax.tick_params(colors=p["label"], labelsize=8.5, length=0)
        ax.grid(True, which="major", color=p["grid_minor"], lw=0.5, ls=":")
        self.shield_canvas.style_legend()
        self.shield_canvas.fig.tight_layout()
        self.shield_canvas.draw()

    def _clear_shielding_chart(self):
        self.shield_canvas.apply_theme(self.is_dark)
        self.shield_canvas.axes.clear()
        self.shield_canvas.draw()

    # --- Doz Hızı (Nokta-Çekirdek Yöntemi) ------------------------------------

    def create_dose_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(10)

        self.dose_frame = QFrame()
        self.dose_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        f_layout = QVBoxLayout(self.dose_frame)
        f_layout.setContentsMargins(12, 10, 12, 12)
        f_layout.setSpacing(10)

        title = QLabel("DOZ HIZI KANALI   ·   NOKTA-ÇEKİRDEK YÖNTEMİ  Ḋ = φ·E·(μₑₙ/ρ)ₕₐᵥₐ")
        title.setObjectName("Engraved")
        f_layout.addWidget(title)

        row1 = QHBoxLayout(); row1.setSpacing(10)

        vbox_src = QVBoxLayout()
        lbl_src = QLabel("KAYNAK"); lbl_src.setObjectName("Engraved")
        vbox_src.addWidget(lbl_src)
        self.dose_source_combo = PanelComboBox(is_dark=self.is_dark)
        self.dose_source_combo.addItem("Serbest Enerji (MeV Gir)", userData="FREE")
        for k, v in RADIOACTIVE_SOURCES.items():
            n_lines = len(v["lines"])
            self.dose_source_combo.addItem(f"{v['name']}  [{n_lines} çizgi]", userData=k)
        self.dose_source_combo.currentIndexChanged.connect(self.handle_dose_source_changed)
        vbox_src.addWidget(self.dose_source_combo)
        row1.addLayout(vbox_src, 3)

        vbox_e = QVBoxLayout()
        lbl_e = QLabel("FOTON ENERJİSİ  [MeV]"); lbl_e.setObjectName("Engraved")
        vbox_e.addWidget(lbl_e)
        self.dose_energy = QLineEdit("1.25")
        self.dose_energy.editingFinished.connect(self.handle_dose_calc)
        vbox_e.addWidget(self.dose_energy)
        row1.addLayout(vbox_e, 2)

        vbox_p = QVBoxLayout()
        lbl_p = QLabel("EMİSYON OLASILIĞI  [foton/bozunma]"); lbl_p.setObjectName("Engraved")
        vbox_p.addWidget(lbl_p)
        self.dose_emission_prob = QLineEdit("1.0")
        self.dose_emission_prob.editingFinished.connect(self.handle_dose_calc)
        vbox_p.addWidget(self.dose_emission_prob)
        row1.addLayout(vbox_p, 2)

        f_layout.addLayout(row1)

        row2 = QHBoxLayout(); row2.setSpacing(10)

        vbox_a = QVBoxLayout()
        lbl_a = QLabel("AKTİVİTE  A"); lbl_a.setObjectName("Engraved")
        vbox_a.addWidget(lbl_a)
        a_row = QHBoxLayout(); a_row.setSpacing(4)
        self.dose_activity = QLineEdit("1.0")
        self.dose_activity.editingFinished.connect(self.handle_dose_calc)
        a_row.addWidget(self.dose_activity, 2)
        self.dose_activity_unit = PanelComboBox(is_dark=self.is_dark)
        for u in DOSE_ACTIVITY_UNITS:
            self.dose_activity_unit.addItem(u, userData=u)
        self.dose_activity_unit.setCurrentIndex(list(DOSE_ACTIVITY_UNITS).index("Ci"))
        self.dose_activity_unit.currentIndexChanged.connect(self.handle_dose_calc)
        a_row.addWidget(self.dose_activity_unit, 1)
        vbox_a.addLayout(a_row)
        row2.addLayout(vbox_a, 3)

        vbox_d = QVBoxLayout()
        lbl_d = QLabel("MESAFE  d"); lbl_d.setObjectName("Engraved")
        vbox_d.addWidget(lbl_d)
        d_row = QHBoxLayout(); d_row.setSpacing(4)
        self.dose_distance = QLineEdit("1.0")
        self.dose_distance.editingFinished.connect(self.handle_dose_calc)
        d_row.addWidget(self.dose_distance, 2)
        self.dose_distance_unit = PanelComboBox(is_dark=self.is_dark)
        for u in DOSE_DISTANCE_UNITS:
            self.dose_distance_unit.addItem(u, userData=u)
        self.dose_distance_unit.setCurrentIndex(list(DOSE_DISTANCE_UNITS).index("m"))
        self.dose_distance_unit.currentIndexChanged.connect(self.handle_dose_calc)
        d_row.addWidget(self.dose_distance_unit, 1)
        vbox_d.addLayout(d_row)
        row2.addLayout(vbox_d, 3)

        vbox_sm = QVBoxLayout()
        lbl_sm = QLabel("ZIRH MALZEMESİ  (OPSİYONEL)"); lbl_sm.setObjectName("Engraved")
        vbox_sm.addWidget(lbl_sm)
        self.dose_shield_combo = PanelComboBox(is_dark=self.is_dark)
        self.dose_shield_combo.addItem("— Zırhsız —", userData="NONE")
        for k, v in SHIELDING_MATERIALS.items():
            self.dose_shield_combo.addItem(f"{v['name']}  [ρ={v['density']:.2f} g/cm³]", userData=k)
        self.dose_shield_combo.currentIndexChanged.connect(self.handle_dose_calc)
        vbox_sm.addWidget(self.dose_shield_combo)
        row2.addLayout(vbox_sm, 3)

        vbox_sx = QVBoxLayout()
        lbl_sx = QLabel("ZIRH KALINLIĞI  x  [cm]"); lbl_sx.setObjectName("Engraved")
        vbox_sx.addWidget(lbl_sx)
        self.dose_shield_thickness = QLineEdit("0.0")
        self.dose_shield_thickness.editingFinished.connect(self.handle_dose_calc)
        vbox_sx.addWidget(self.dose_shield_thickness)
        row2.addLayout(vbox_sx, 2)

        f_layout.addLayout(row2)

        self.dose_source_info = QLabel()
        self.dose_source_info.setObjectName("Engraved")
        self.dose_source_info.setWordWrap(True)
        self.dose_source_info.setVisible(False)
        f_layout.addWidget(self.dose_source_info)

        layout.addWidget(self.dose_frame)

        btn_row = QHBoxLayout()
        calc_btn = QPushButton("HESAPLA  ·  DOZ HIZI")
        calc_btn.setObjectName("PrimaryLever")
        calc_btn.clicked.connect(self.handle_dose_calc)
        btn_row.addWidget(calc_btn)
        self.dose_lamp = IndicatorLamp("ZIRH UYGULANDI", is_dark=self.is_dark)
        btn_row.addSpacing(16)
        btn_row.addWidget(self.dose_lamp)
        btn_row.addStretch(1)
        btn_row.addWidget(self._make_export_button(
            lambda: self._export_report_pdf(self.dose_output.toPlainText, "doz_hizi_raporu.pdf")))
        layout.addLayout(btn_row)

        result_row = QHBoxLayout(); result_row.setSpacing(10)
        self.dose_output = QTextEdit()
        self.dose_output.setReadOnly(True)
        result_row.addWidget(self.dose_output, 1)

        self.dose_canvas = MplCanvas()
        result_row.addWidget(self.dose_canvas, 1)
        layout.addLayout(result_row, 1)

        self.handle_dose_source_changed()
        return w

    def handle_dose_source_changed(self):
        source_key = self.dose_source_combo.currentData() or "FREE"
        is_free = (source_key == "FREE")
        self.dose_energy.setEnabled(is_free)
        self.dose_emission_prob.setEnabled(is_free)

        if is_free:
            self.dose_source_info.setVisible(False)
        else:
            src = RADIOACTIVE_SOURCES[source_key]
            parts = []
            for l in src["lines"]:
                e_disp = f"{l['energy']*1000:.1f} keV" if l["energy"] < 1.0 else f"{l['energy']:.3f} MeV"
                parts.append(f"{e_disp} (%{l['intensity']*100:.2f})")
            self.dose_source_info.setText(
                f"KAYNAK ÇİZGİLERİ (IAEA referans şiddetleri): {'   ·   '.join(parts)}"
            )
            self.dose_source_info.setVisible(True)

        self.handle_dose_calc()

    def handle_dose_calc(self):
        source_key = self.dose_source_combo.currentData() or "FREE" if hasattr(self, "dose_source_combo") else "FREE"
        shield_key = self.dose_shield_combo.currentData() or "NONE"
        activity_unit = self.dose_activity_unit.currentData() or "Bq"
        distance_unit = self.dose_distance_unit.currentData() or "cm"

        try:
            activity_val = float(self.dose_activity.text().replace(",", "."))
            distance_val = float(self.dose_distance.text().replace(",", "."))
            shield_x = float(self.dose_shield_thickness.text().replace(",", "."))
            if activity_val <= 0 or distance_val <= 0 or shield_x < 0:
                raise ValueError
        except ValueError:
            self.dose_output.setPlainText(" ARIZA · AKTİVİTE / MESAFE / KALINLIK DEĞERİ GEÇERSİZ")
            self._clear_dose_chart()
            self.dose_lamp.set_state("idle")
            return

        activity_Bq = activity_val * DOSE_ACTIVITY_UNITS[activity_unit]
        distance_cm = distance_val * DOSE_DISTANCE_UNITS[distance_unit]
        shield_mat = None if shield_key == "NONE" else shield_key

        if source_key == "FREE":
            try:
                energy = float(self.dose_energy.text().replace(",", "."))
                prob = float(self.dose_emission_prob.text().replace(",", "."))
                if energy <= 0 or prob <= 0:
                    raise ValueError
            except ValueError:
                self.dose_output.setPlainText(" ARIZA · FOTON ENERJİSİ / EMİSYON OLASILIĞI GEÇERSİZ")
                self._clear_dose_chart()
                self.dose_lamp.set_state("idle")
                return
            self._calc_dose_single(activity_Bq, activity_val, activity_unit, energy, prob,
                                    distance_cm, distance_val, distance_unit, shield_mat, shield_x)
        else:
            self._calc_dose_multiline(activity_Bq, activity_val, activity_unit, source_key,
                                       distance_cm, distance_val, distance_unit, shield_mat, shield_x)

    def _calc_dose_single(self, activity_Bq, activity_val, activity_unit, energy, prob,
                           distance_cm, distance_val, distance_unit, shield_mat, shield_x):
        r = dose_rate_point_source(activity_Bq, energy, prob, distance_cm, shield_mat, shield_x)

        rap  = " SİNOP ÜNİVERSİTESİ · NÜKLEER ENERJİ MÜHENDİSLİĞİ\n"
        rap += " DOZ HIZI ÖLÇÜM FİŞİ  (NOKTA-ÇEKİRDEK)" + " " * 3 + "FORM SNU-10/56\n"
        rap += "=" * (W_FIS + 2) + "\n"
        rap += _fis("KAYNAK", "Serbest enerji (kaynak seçilmedi)")
        rap += _fis("AKTİVİTE  A", f"{activity_val:12.4f}", activity_unit)
        rap += _fis("FOTON ENERJİSİ  E", f"{energy:12.4f}", "MeV")
        rap += _fis("EMİSYON OLASILIĞI  p", f"{prob:12.4f}", "foton/bozunma")
        rap += _fis("MESAFE  d", f"{distance_val:12.4f}", distance_unit)
        rap += "-" * (W_FIS + 2) + "\n"
        rap += self._dose_shield_report_block(r, shield_mat, shield_x)
        rap += "-" * (W_FIS + 2) + "\n"
        rap += _fis("FOTON AKI YOĞUNLUĞU  φ", f"{r['phi']:12.4e}", "foton/(cm²·s)")
        rap += _fis("HAVA μₑₙ/ρ", f"{r['mu_en_rho']:12.5f}", "cm²/g")
        rap += "-" * (W_FIS + 2) + "\n"
        rap += " SONUÇ (HAVA KERMA HIZI ≈ EŞDEĞER DOZ HIZI, w_R=1)\n"
        rap += _fis("DOZ HIZI", f"{r['dose_rate_uSv_h']:12.4f}", "µSv/sa")
        rap += _fis("DOZ HIZI", f"{r['dose_rate_mGy_h']:12.5f}", "mGy/sa")
        rap += _fis("DOZ HIZI", f"{r['dose_rate_mrem_h']:12.4f}", "mrem/sa")
        rap += "\n" + "=" * (W_FIS + 2) + "\n"
        rap += self._dose_report_footer()
        self.dose_output.setPlainText(rap)
        self._set_dose_lamp(shield_mat, shield_x)
        self._draw_dose_chart([r["energy"]], [r["emission_prob"]], activity_Bq, distance_cm,
                               shield_mat, shield_x, r["dose_rate_uSv_h"])

    def _set_dose_lamp(self, shield_mat, shield_x):
        if shield_mat is None or shield_x <= 0:
            self.dose_lamp.set_state("idle", "ZIRH UYGULANMADI")
        else:
            self.dose_lamp.set_state("ok", "ZIRH UYGULANDI")

    def _calc_dose_multiline(self, activity_Bq, activity_val, activity_unit, source_key,
                              distance_cm, distance_val, distance_unit, shield_mat, shield_x):
        m = dose_rate_multiline(activity_Bq, source_key, distance_cm, shield_mat, shield_x)
        src = RADIOACTIVE_SOURCES[source_key]

        rap  = " SİNOP ÜNİVERSİTESİ · NÜKLEER ENERJİ MÜHENDİSLİĞİ\n"
        rap += " DOZ HIZI ÖLÇÜM FİŞİ  (NOKTA-ÇEKİRDEK)" + " " * 3 + "FORM SNU-10/56\n"
        rap += "=" * (W_FIS + 2) + "\n"
        rap += _fis("KAYNAK", m["source"])
        rap += _fis("AKTİVİTE  A", f"{activity_val:12.4f}", activity_unit)
        rap += _fis("MESAFE  d", f"{distance_val:12.4f}", distance_unit)
        rap += "-" * (W_FIS + 2) + "\n"
        rap += self._dose_shield_report_block(None, shield_mat, shield_x)
        rap += "-" * (W_FIS + 2) + "\n"
        rap += " ÇİZGİ BAZLI DOZ KATKISI (HER ENERJİ AYRI HESAPLANIR)\n"
        rap += "-" * (W_FIS + 2) + "\n"
        for l in m["lines"]:
            e_disp = f"{l['energy']*1000:.1f} keV" if l["energy"] < 1.0 else f"{l['energy']:.4f} MeV"
            rap += _fis(f"  {e_disp}  (p=%{l['emission_prob']*100:.2f})",
                        f"{l['dose_rate_uSv_h']:12.4f}", "µSv/sa")
        rap += "-" * (W_FIS + 2) + "\n"
        rap += " TOPLAM SONUÇ  (TÜM ÇİZGİLERİN DOZ KATKISI TOPLANIR)\n"
        rap += _fis("DOZ HIZI", f"{m['dose_rate_uSv_h']:12.4f}", "µSv/sa")
        rap += _fis("DOZ HIZI", f"{m['dose_rate_mGy_h']:12.5f}", "mGy/sa")
        rap += _fis("DOZ HIZI", f"{m['dose_rate_mrem_h']:12.4f}", "mrem/sa")
        rap += "\n" + "=" * (W_FIS + 2) + "\n"
        rap += self._dose_report_footer()
        self.dose_output.setPlainText(rap)
        self._set_dose_lamp(shield_mat, shield_x)
        energies = [l["energy"] for l in src["lines"]]
        probs = [l["intensity"] for l in src["lines"]]
        self._draw_dose_chart(energies, probs, activity_Bq, distance_cm,
                               shield_mat, shield_x, m["dose_rate_uSv_h"])

    def _dose_shield_report_block(self, r_single, shield_mat, shield_x) -> str:
        if shield_mat is None or shield_x <= 0:
            return " ZIRH: Yok (serbest alanda / havada doz hızı)\n"
        rap = f" ZIRH: {SHIELDING_MATERIALS[shield_mat]['name']}, x = {shield_x:.3f} cm\n"
        if r_single is not None and r_single["shield"] is not None:
            sh = r_single["shield"]
            mode = "geniş-demet (buildup dahil)" if sh["I_broad"] is not None else "dar-demet (buildup verisi yok)"
            rap += f" ZIRH GEÇİRGENLİĞİ ({mode}): %{r_single['transmission']*100:.4f}\n"
        return rap

    def _dose_report_footer(self) -> str:
        return (
            " NOT: μₑₙ/ρ (havanın kütle enerji-soğurma katsayısı) NIST X-Ray/Gamma-Ray\n"
            " Mass Attenuation Coefficients (dry air) referans tablosundan log-log\n"
            " enterpolasyonla alınmıştır — ZIRHLAMA modülündeki μ/ρ (zayıflama)\n"
            " katsayısından FARKLI bir büyüklüktür. Nokta-çekirdek yöntemi: φ=A·p/(4πd²),\n"
            " Ḋ=φ·E·(μₑₙ/ρ)·1.602176634e-10 Gy/s. Fotonlar için w_R=1 kabulüyle eşdeğer\n"
            " doz (Sv) sayısal olarak hava kermaya (Gy) eşit alınmıştır. Kaynak gama\n"
            " çizgileri IAEA referans şiddet veritabanındandır (≥%1 şiddetli çizgiler).\n"
            f" ÖLÇÜM {QDateTime.currentDateTime().toString('dd.MM.yyyy HH:mm:ss')}"
            "   ·   OPERATÖR: A.S.ÇETİN\n"
        )

    def _draw_dose_chart(self, energies, probs, activity_Bq, distance_cm, shield_mat, shield_x, current_uSv_h):
        p = palette(self.is_dark)
        self.dose_canvas.apply_theme(self.is_dark)
        ax = self.dose_canvas.axes
        ax.clear()

        dmax = max(distance_cm * 6.0, 300.0)
        dmin = max(distance_cm * 0.1, 1.0)
        ds = np.geomspace(dmin, dmax, 200)

        def total_dose_at(d):
            total = 0.0
            for e, pr in zip(energies, probs):
                rr = dose_rate_point_source(activity_Bq, e, pr, d, shield_mat, shield_x)
                total += rr["dose_rate_uSv_h"]
            return total

        ys = np.array([total_dose_at(d) for d in ds])
        ax.plot(ds / 100.0, ys, color=p["accent"], lw=2.2, zorder=3, label="Doz Hızı  ∝ 1/d²")
        ax.scatter([distance_cm / 100.0], [max(current_uSv_h, 1e-9)], color=p["accent"], s=55, zorder=5)
        ax.annotate(f"d={distance_cm/100.0:.2f} m", (distance_cm / 100.0, max(current_uSv_h, 1e-9)),
                    textcoords="offset points", xytext=(8, 10), color=p["text"],
                    fontfamily="monospace", fontsize=8, zorder=5)

        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlabel("MESAFE  d  [m]", color=p["label"], fontfamily="monospace", fontsize=9)
        ax.set_ylabel("DOZ HIZI  [µSv/sa]", color=p["label"], fontfamily="monospace", fontsize=9)
        ax.set_facecolor(p["inset"])
        for s in ax.spines.values():
            s.set_color(p["frame"]); s.set_linewidth(1.1)
        ax.tick_params(colors=p["label"], labelsize=8.5, length=0)
        ax.grid(True, which="major", color=p["grid_minor"], lw=0.5, ls=":")
        self.dose_canvas.style_legend()
        self.dose_canvas.readout_box(f"Ḋ = {current_uSv_h:.3f} µSv/sa\nd = {distance_cm/100.0:.3f} m")
        self.dose_canvas.fig.tight_layout()
        self.dose_canvas.draw()

    def _clear_dose_chart(self):
        self.dose_canvas.apply_theme(self.is_dark)
        self.dose_canvas.axes.clear()
        self.dose_canvas.draw()

    # --- Kritiklik / Dört-Faktör Formülü (k∞ = η·f·p·ε) -----------------------

    def create_criticality_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(10)

        self.crit_frame = QFrame()
        self.crit_frame.setStyleSheet(panel_qss(self.is_dark, accent_rule=True))
        f_layout = QVBoxLayout(self.crit_frame)
        f_layout.setContentsMargins(12, 10, 12, 12)
        f_layout.setSpacing(10)

        title = QLabel("KRİTİKLİK KANALI   ·   DÖRT-FAKTÖR FORMÜLÜ  k∞ = η·f·p·ε")
        title.setObjectName("Engraved")
        f_layout.addWidget(title)

        row1 = QHBoxLayout(); row1.setSpacing(10)

        vbox_mod = QVBoxLayout()
        lbl_mod = QLabel("MODERATÖR"); lbl_mod.setObjectName("Engraved")
        vbox_mod.addWidget(lbl_mod)
        self.crit_mod_combo = PanelComboBox(is_dark=self.is_dark)
        for k, v in MODERATORS.items():
            self.crit_mod_combo.addItem(v["name"], userData=k)
        self.crit_mod_combo.currentIndexChanged.connect(self.handle_criticality_calc)
        vbox_mod.addWidget(self.crit_mod_combo)
        row1.addLayout(vbox_mod, 3)

        vbox_e = QVBoxLayout()
        lbl_e = QLabel("YAKIT ZENGİNLEŞTİRMESİ  [U-235 atom %]"); lbl_e.setObjectName("Engraved")
        vbox_e.addWidget(lbl_e)
        self.crit_enrichment = QLineEdit("3.2")
        self.crit_enrichment.editingFinished.connect(self.handle_criticality_calc)
        vbox_e.addWidget(self.crit_enrichment)
        row1.addLayout(vbox_e, 2)

        vbox_r = QVBoxLayout()
        lbl_r = QLabel("YAKIT/MODERATÖR ATOM ORANI  N_F/N_M"); lbl_r.setObjectName("Engraved")
        vbox_r.addWidget(lbl_r)
        self.crit_nf_nm = QLineEdit("0.02")
        self.crit_nf_nm.editingFinished.connect(self.handle_criticality_calc)
        vbox_r.addWidget(self.crit_nf_nm)
        row1.addLayout(vbox_r, 2)

        f_layout.addLayout(row1)

        row2 = QHBoxLayout(); row2.setSpacing(10)

        vbox_i = QVBoxLayout()
        lbl_i = QLabel("EFEKTİF REZONANS İNTEGRALİ  I_eff  [barn]"); lbl_i.setObjectName("Engraved")
        vbox_i.addWidget(lbl_i)
        self.crit_ieff = QLineEdit("10.0")
        self.crit_ieff.editingFinished.connect(self.handle_criticality_calc)
        vbox_i.addWidget(self.crit_ieff)
        row2.addLayout(vbox_i, 2)

        vbox_eps = QVBoxLayout()
        lbl_eps = QLabel("HIZLI FİSYON FAKTÖRÜ  ε"); lbl_eps.setObjectName("Engraved")
        vbox_eps.addWidget(lbl_eps)
        self.crit_epsilon = QLineEdit("1.00")
        self.crit_epsilon.editingFinished.connect(self.handle_criticality_calc)
        vbox_eps.addWidget(self.crit_epsilon)
        row2.addLayout(vbox_eps, 2)

        row2.addStretch(3)
        f_layout.addLayout(row2)

        note = QLabel(
            "NOT: I_eff ders kitabınızdan/verilen değerden girilmelidir — homojen karışımlar için "
            "güvenilir bir I_eff(N_M/N₂₈) ampirik korelasyonu doğrulanamadı (uydurulmadı). "
            "ε, homojen sistemlerde tipik olarak ≈1.00'dır (heterojen/topaklı yakıtta ~1.02-1.03)."
        )
        note.setObjectName("Engraved")
        note.setWordWrap(True)
        f_layout.addWidget(note)

        layout.addWidget(self.crit_frame)

        btn_row = QHBoxLayout()
        calc_btn = QPushButton("HESAPLA  ·  k∞")
        calc_btn.setObjectName("PrimaryLever")
        calc_btn.clicked.connect(self.handle_criticality_calc)
        btn_row.addWidget(calc_btn)
        self.crit_lamp = IndicatorLamp("k∞ > 1  ·  KRİTİK OLABİLİR", is_dark=self.is_dark)
        btn_row.addSpacing(16)
        btn_row.addWidget(self.crit_lamp)
        btn_row.addStretch(1)
        btn_row.addWidget(self._make_export_button(
            lambda: self._export_report_pdf(self.crit_output.toPlainText, "kritiklik_raporu.pdf")))
        layout.addLayout(btn_row)

        result_row = QHBoxLayout(); result_row.setSpacing(10)
        self.crit_output = QTextEdit()
        self.crit_output.setReadOnly(True)
        result_row.addWidget(self.crit_output, 1)

        self.crit_canvas = MplCanvas()
        result_row.addWidget(self.crit_canvas, 1)
        layout.addLayout(result_row, 1)

        self.handle_criticality_calc()
        return w

    def handle_criticality_calc(self):
        mod_key = self.crit_mod_combo.currentData() or "H2O" if hasattr(self, "crit_mod_combo") else "H2O"
        try:
            enrichment = float(self.crit_enrichment.text().replace(",", "."))
            nf_nm = float(self.crit_nf_nm.text().replace(",", "."))
            i_eff = float(self.crit_ieff.text().replace(",", "."))
            epsilon = float(self.crit_epsilon.text().replace(",", "."))
            if not (0.0 < enrichment <= 100.0) or nf_nm <= 0 or i_eff < 0 or epsilon <= 0:
                raise ValueError
        except ValueError:
            self.crit_output.setPlainText(" ARIZA · GİRDİ DEĞERLERİ GEÇERSİZ (zenginleştirme 0-100%, oran/I_eff/ε > 0 olmalı)")
            self._clear_criticality_chart()
            return

        k = k_infinity(mod_key, enrichment, nf_nm, i_eff, epsilon)
        if k is None:
            self.crit_output.setPlainText(" ARIZA · HESAPLANAMADI (girdi kombinasyonunu kontrol edin)")
            self._clear_criticality_chart()
            self.crit_lamp.set_state("idle")
            return
        self.crit_lamp.set_state("ok" if k["k_inf"] > 1.0 else "fail")

        mod_name = MODERATORS[mod_key]["name"]
        rap  = " SİNOP ÜNİVERSİTESİ · NÜKLEER ENERJİ MÜHENDİSLİĞİ\n"
        rap += " KRİTİKLİK / DÖRT-FAKTÖR FORMÜLÜ FİŞİ" + " " * 3 + "FORM SNU-11/56\n"
        rap += "=" * (W_FIS + 2) + "\n"
        rap += _fis("MODERATÖR", mod_name)
        rap += _fis("YAKIT ZENGİNLEŞTİRMESİ", f"{enrichment:12.4f}", "% (U-235 atom)")
        rap += _fis("YAKIT/MODERATÖR ORANI  N_F/N_M", f"{nf_nm:12.5f}")
        rap += _fis("EFEKTİF REZONANS İNTEGRALİ  I_eff", f"{i_eff:12.4f}", "barn")
        rap += "-" * (W_FIS + 2) + "\n"
        rap += " ARA DEĞERLER\n"
        rap += _fis("YAKIT  σa (ORTALAMA)", f"{k['sigma_a_fuel']:12.4f}", "barn")
        rap += _fis("YAKIT  σf (ORTALAMA)", f"{k['sigma_f_fuel']:12.4f}", "barn")
        rap += _fis("MODERATÖR  σa (TÜRETİLMİŞ)", f"{k['sigma_a_M']:12.5f}", "barn")
        rap += _fis("YAKITTA  N₂₈ (U-238 ATOM YOĞ.)", f"{k['N_28']:12.4e}", "1/cm³")
        rap += _fis("MODERATÖR  ξΣs", f"{k['xi_sigma_s']:12.4f}", "cm⁻¹")
        rap += "-" * (W_FIS + 2) + "\n"
        rap += " DÖRT-FAKTÖR FORMÜLÜ\n"
        rap += _fis("η  ÜRETİM (REPRODÜKSİYON) FAKTÖRÜ", f"{k['eta']:12.5f}")
        rap += _fis("f  TERMAL KULLANIM FAKTÖRÜ", f"{k['f']:12.5f}")
        rap += _fis("p  REZONANS KAÇMA OLASILIĞI", f"{k['p']:12.5f}")
        rap += _fis("ε  HIZLI FİSYON FAKTÖRÜ", f"{k['epsilon']:12.5f}")
        rap += "-" * (W_FIS + 2) + "\n"
        rap += _fis("k∞  SONSUZ ORTAM ÇOĞALMA FAKTÖRÜ", f"{k['k_inf']:12.5f}")
        durum = "KRİTİK OLABİLİR (k∞>1)" if k["k_inf"] > 1.0 else "KRİTİK OLAMAZ (k∞≤1, sızıntısız halde bile)"
        rap += _fis("DEĞERLENDİRME", durum)
        rap += "\n" + "=" * (W_FIS + 2) + "\n"
        rap += self._criticality_report_footer()
        self.crit_output.setPlainText(rap)
        self._draw_criticality_chart(k)

    def _criticality_report_footer(self) -> str:
        return (
            " NOT: η/f formülleri ve U-235/U-238 2200 m/s tesir kesitleri (σa=680.8/2.70 b,\n"
            " σf=582.2 b, ν=2.42) standart literatür değerleridir (çapraz doğrulandı). p için\n"
            " genel üstel form p=exp(−N₂₈·I_eff/ξΣs) doğrulanmıştır; I_eff KULLANICI GİRDİSİDİR\n"
            " (homojen karışım için güvenilir ampirik korelasyon bulunamadığından uydurulmadı).\n"
            " k_eff = k∞·(sızıntısız olasılık) BU MODÜLDE HESAPLANMAZ — geometri/B² gerektirir,\n"
            " ayrı bir modül olarak planlanabilir. Bu sonuç yalnız SONSUZ ORTAM (k∞) içindir.\n"
            f" ÖLÇÜM {QDateTime.currentDateTime().toString('dd.MM.yyyy HH:mm:ss')}"
            "   ·   OPERATÖR: A.S.ÇETİN\n"
        )

    def _draw_criticality_chart(self, k):
        p = palette(self.is_dark)
        self.crit_canvas.apply_theme(self.is_dark)
        ax = self.crit_canvas.axes
        ax.clear()

        labels = ["η", "f", "p", "ε", "k∞"]
        values = [k["eta"], k["f"], k["p"], k["epsilon"], k["k_inf"]]
        colors = [p["accent"], p["blue"], p["green"], p["red"], p["text"]]

        bars = ax.bar(labels, values, color=colors, zorder=3, width=0.6)
        ax.axhline(1.0, color=p["faint"], lw=1.3, ls="--", zorder=2)
        ax.annotate("k=1 (KRİTİK EŞİK)", (0, 1.0), textcoords="offset points",
                    xytext=(4, 6), color=p["label"], fontfamily="monospace", fontsize=8)

        for b, v in zip(bars, values):
            ax.annotate(f"{v:.4f}", (b.get_x() + b.get_width() / 2, v),
                        textcoords="offset points", xytext=(0, 6), ha="center",
                        color=p["text"], fontfamily="monospace", fontsize=8.5, zorder=5)

        ymax = max(values + [1.0]) * 1.25
        ax.set_ylim(0, ymax)
        ax.set_ylabel("DEĞER", color=p["label"], fontfamily="monospace", fontsize=9)
        ax.set_facecolor(p["inset"])
        for s in ax.spines.values():
            s.set_color(p["frame"]); s.set_linewidth(1.1)
        ax.tick_params(colors=p["label"], labelsize=9.5, length=0)
        ax.grid(True, which="major", axis="y", color=p["grid_minor"], lw=0.5, ls=":")

        self.crit_canvas.readout_box(f"k∞ = {k['k_inf']:.4f}")
        self.crit_canvas.fig.tight_layout()
        self.crit_canvas.draw()

    def _clear_criticality_chart(self):
        self.crit_canvas.apply_theme(self.is_dark)
        self.crit_canvas.axes.clear()
        self.crit_canvas.draw()


if __name__ == "__main__":
    import time

    app = QApplication(sys.argv)
    # Splash -> ana pencere geçişi sırasında, splash kapanırken görünür başka
    # pencere olmayabileceği an için Qt'nin varsayılan "son pencere kapanınca
    # çık" davranışını geçici olarak kapatıyoruz (SplashScreen zaten ana
    # pencereyi önce gösterip sonra kendini kapatacak şekilde sıralanmıştır —
    # bu ikinci bir güvenlik katmanıdır). Ana pencere kapatılınca uygulamanın
    # yine de kapanması için tekrar açılır.
    app.setQuitOnLastWindowClosed(False)

    # --- Açılış ekranı: logo + kısa "sistem başlatılıyor" telemetrisi, sonra
    # ana pencereye yumuşak geçiş. Adımlar sahte bir zamanlayıcıyla değil,
    # ana pencerenin __init__'indeki gerçek kuruluş noktalarında ilerler;
    # kuruluş çok hızlı biterse bile (küçük veritabanı, hızlı makine) açılış
    # ekranı en az MIN_SPLASH_MS boyunca ekranda kalır ki göze çarpsın.
    base_dir = os.path.dirname(os.path.abspath(__file__))
    splash_pix = (load_smart_pixmap(os.path.join(base_dir, "lacivertlogo.png"), 92)
                  or load_smart_pixmap(os.path.join(base_dir, "logo.png"), 92))
    splash = SplashScreen(splash_pix)
    splash.show()
    splash.set_step(0)
    app.processEvents()

    t0 = time.monotonic()
    window = NuclearCalculatorApp(splash=splash)
    elapsed_ms = (time.monotonic() - t0) * 1000.0

    MIN_SPLASH_MS = 900
    remaining_ms = max(0, int(MIN_SPLASH_MS - elapsed_ms))

    splash.finished.connect(window.show)
    splash.finished.connect(lambda: app.setQuitOnLastWindowClosed(True))
    QTimer.singleShot(remaining_ms, splash.start_fade_out)

    sys.exit(app.exec())
