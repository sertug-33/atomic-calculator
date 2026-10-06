# -*- coding: utf-8 -*-
"""
Atom Yoğunluğu Modülü  (syllabus: "Nükleer Reaksiyonlar I", madde 4)

Üç alt başlık:
  a) Elementler için:              N = rho * N_A / M
  b) Bileşikler için:              N_molekul = rho * N_A / M_bilesik ; N_i = n_i * N_molekul
  c) Zenginleştirilmiş bileşikler: U'nun izotopik karışımı ağırlık yüzdesinden
                                    efektif atomik ağırlığa çevrilir, sonra (b) uygulanır.

Atomik ağırlıklar ELLE girilmez: ABUNDANCE (doğal bolluk, NUBASE2020) ve ISOTOPES
(izotop kütleleri, AME2020) üzerinden bolluk-ağırlıklı ortalama olarak hesaplanır.
Yoğunluk (rho) fiziksel bir malzeme özelliğidir, nükleer veriden çıkmaz — kullanıcı
(ödev/problem) tarafından verilir; MATERIALS kataloğu sadece tipik/ders kitabı
değerleriyle bir başlangıç noktası sunar, kullanıcı istediği rho'yu yazabilir.
"""
import re
from database import ISOTOPES, ABUNDANCE

N_A = 6.02214076e23  # Avogadro sayısı (1/mol) — CODATA 2018

# --- UI'daki "MALZEME KANALI" açılır menüsü için hazır katalog ---------------
# rho: tipik/teorik yoğunluk (g/cm^3, ders kitabı değeri) — kullanıcı değiştirebilir.
# "U" formülünde geçen her malzeme zenginleştirme kutusunu otomatik aktif eder.
MATERIALS = {
    "U_METAL": {"name": "Uranyum Metal",              "formula": "U",    "rho": 19.10},
    "UO2":     {"name": "Uranyum Dioksit (UO₂)",       "formula": "UO2",  "rho": 10.97},
    "U3O8":    {"name": "Triuranyum Oktaoksit (U₃O₈)", "formula": "U3O8", "rho": 8.30},
    "UC":      {"name": "Uranyum Karbür (UC)",         "formula": "UC",  "rho": 13.63},
    "H2O":     {"name": "Hafif Su (H₂O)",               "formula": "H2O", "rho": 1.000},
    "D2O":     {"name": "Ağır Su (D₂O)",                "formula": "D2O", "rho": 1.104},
    "C":       {"name": "Reaktör Grafiti (C)",          "formula": "C",   "rho": 1.70},
    "BE":      {"name": "Berilyum (Be)",                "formula": "Be",  "rho": 1.85},
    "ZR":      {"name": "Zirkonyum (Zircaloy ≈ Zr)",    "formula": "Zr",  "rho": 6.56},
    "AL":      {"name": "Alüminyum (Al)",               "formula": "Al",  "rho": 2.70},
    "FE":      {"name": "Demir / Çelik (Fe)",           "formula": "Fe",  "rho": 7.87},
}

# --- Element sembolü <-> Z eşlemesi (NUBASE2020'den türetilir) --------------
_SYMBOL_TO_Z = None


def _symbol_to_z():
    global _SYMBOL_TO_Z
    if _SYMBOL_TO_Z is None:
        _SYMBOL_TO_Z = {}
        for k, v in ISOTOPES.items():
            if k in ("n", "p"):
                continue
            sym = v["name"].split("-")[0]
            _SYMBOL_TO_Z.setdefault(sym, v["Z"])
    return _SYMBOL_TO_Z


def _mass_by_za():
    d = {}
    for k, v in ISOTOPES.items():
        if k in ("n", "p"):
            continue
        d[(v["Z"], v["A"])] = v["mass"]
    return d


_MASS_BY_ZA = _mass_by_za()

_Z_TO_SYMBOL = {}
for _k, _v in ISOTOPES.items():
    if _k in ("n", "p"):
        continue
    _Z_TO_SYMBOL.setdefault(_v["Z"], _v["name"].split("-")[0])

# D (döteryum) ve T (trityum): doğal element karışımı DEĞİL, saf tek izotop.
# Ağır su gibi formüllerde "D2O" yazılabilsin diye özel olarak çözülür.
_SPECIAL_ISOTOPES = {"D": "2H", "T": "3H"}


def atomic_weight(Z: int):
    """Doğal izotopik bolluğa göre ağırlıklı ortalama atomik ağırlık (g/mol).
    NUBASE2020'nin IS= (isotopic abundance) alanlarından hesaplanır."""
    isos = ABUNDANCE.get(Z)
    if not isos:
        return None
    tot_ab = sum(ab for _, ab in isos)
    if tot_ab <= 0:
        return None
    num = sum(ab * _MASS_BY_ZA[(Z, A)] for A, ab in isos if (Z, A) in _MASS_BY_ZA)
    return num / tot_ab


def element_symbol(Z: int) -> str:
    return _Z_TO_SYMBOL.get(Z, f"Z{Z}")


def element_atom_density(rho: float, Z: int):
    """N = rho * N_A / M   (atom/cm^3). rho: g/cm^3."""
    M = atomic_weight(Z)
    if M is None:
        raise ValueError(f"Z={Z} için doğal bolluk verisi yok.")
    return {"M": M, "N": rho * N_A / M, "symbol": element_symbol(Z)}


def _resolve_token(sym: str):
    """Formüldeki bir sembolü (Z, atomik_ağırlık, döngü_anahtarı) üçlüsüne çözer.
    'D'/'T' -> saf 2H/3H izotop kütlesi (doğal H karışımı DEĞİL)."""
    if sym in _SPECIAL_ISOTOPES:
        iso_key = _SPECIAL_ISOTOPES[sym]
        return 1, ISOTOPES[iso_key]["mass"], sym
    Z = _symbol_to_z().get(sym)
    if Z is None:
        raise ValueError(f"Element/izotop tanınmadı: '{sym}'")
    M = atomic_weight(Z)
    if M is None:
        raise ValueError(f"'{sym}' için doğal bolluk verisi yok.")
    return Z, M, sym


# --- Kimyasal formül ayrıştırma ---------------------------------------------
_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(formula: str):
    """'UO2' -> [('U',1), ('O',2)] ; 'D2O' -> [('D',2), ('O',1)]"""
    formula = formula.strip()
    out, pos = [], 0
    while pos < len(formula):
        m = _FORMULA_TOKEN.match(formula, pos)
        if not m or not m.group(1):
            raise ValueError(f"Formül çözümlenemedi: '{formula}' (konum {pos})")
        out.append((m.group(1), int(m.group(2)) if m.group(2) else 1))
        pos = m.end()
    if not out:
        raise ValueError(f"Boş formül: '{formula}'")
    return out


def compound_atom_density(rho: float, formula: str):
    """Doğal (zenginleştirilmemiş) bileşik/element için atom başına yoğunluk.
    Tek elementli formül ('Fe' gibi) element_atom_density ile aynı sonucu verir."""
    parts = parse_formula(formula)
    M_formula, resolved = 0.0, []
    for sym, count in parts:
        Z, M, key = _resolve_token(sym)
        M_formula += count * M
        resolved.append((key, Z, count, M))

    N_molecule = rho * N_A / M_formula
    elements = {key: {"Z": Z, "count": count, "M": M, "N": N_molecule * count}
                for key, Z, count, M in resolved}
    return {"M_formula": M_formula, "N_molecule": N_molecule, "elements": elements}


def enriched_uranium_weight(w235_pct: float):
    """U-235 AĞIRLIK yüzdesinden karışımın efektif atomik ağırlığı ve
    U-235/U-238 ATOM (mol) yüzdeleri.  1/M_mix = w235/M235 + w238/M238"""
    M235 = _MASS_BY_ZA[(92, 235)]
    M238 = _MASS_BY_ZA[(92, 238)]
    w235 = w235_pct / 100.0
    w238 = 1.0 - w235
    M_mix = 1.0 / (w235 / M235 + w238 / M238)
    n235, n238 = w235 / M235, w238 / M238
    tot_n = n235 + n238
    return {"M_mix": M_mix, "atom_pct_235": 100.0 * n235 / tot_n,
            "atom_pct_238": 100.0 * n238 / tot_n, "M235": M235, "M238": M238}


def enriched_compound_atom_density(rho: float, formula: str, w235_pct: float):
    """Formülde 'U' geçen HERHANGİ bir bileşik/element için (UO2, U3O8, UC, U metal...)
    zenginleştirmeyi hesaba katan tam çözüm. U atom yoğunluğu 235U/238U'ya bölünür."""
    mix = enriched_uranium_weight(w235_pct)
    parts = parse_formula(formula)
    if not any(sym == "U" for sym, _ in parts):
        raise ValueError("Formülde 'U' yok; zenginleştirme uygulanamaz.")

    M_formula, resolved = 0.0, []
    for sym, count in parts:
        if sym == "U":
            Z, M, key = 92, mix["M_mix"], "U"
        else:
            Z, M, key = _resolve_token(sym)
        M_formula += count * M
        resolved.append((key, Z, count, M))

    N_molecule = rho * N_A / M_formula
    elements = {}
    for key, Z, count, M in resolved:
        N_total = N_molecule * count
        if key == "U":
            elements["U"] = {
                "Z": 92, "count": count, "M": M, "N": N_total,
                "N_235": N_total * mix["atom_pct_235"] / 100.0,
                "N_238": N_total * mix["atom_pct_238"] / 100.0,
            }
        else:
            elements[key] = {"Z": Z, "count": count, "M": M, "N": N_total}
    return {"M_formula": M_formula, "N_molecule": N_molecule, "elements": elements, "mix": mix}


def calculate_material(rho: float, formula: str, w235_pct: float | None = None):
    """UI'nin çağıracağı tek giriş noktası. Formülde 'U' varsa ve w235_pct
    verilmişse zenginleştirilmiş yol, aksi halde doğal bolluk yolu kullanılır."""
    has_u = any(sym == "U" for sym, _ in parse_formula(formula))
    if has_u and w235_pct is not None:
        return enriched_compound_atom_density(rho, formula, w235_pct)
    return compound_atom_density(rho, formula)


if __name__ == "__main__":
    r = element_atom_density(19.1, 92)
    print(f"Doğal U metal: M={r['M']:.3f} g/mol, N={r['N']:.4e} atom/cm^3")

    r = compound_atom_density(10.97, "UO2")
    print(f"Doğal UO2: M={r['M_formula']:.3f} g/mol, "
          f"N_U={r['elements']['U']['N']:.4e}, N_O={r['elements']['O']['N']:.4e}")

    r = enriched_compound_atom_density(10.5, "UO2", 3.2)
    u = r["elements"]["U"]
    print(f"%3.2 zg. UO2: M_U={u['M']:.3f}, N_U235={u['N_235']:.4e}, "
          f"N_U238={u['N_238']:.4e}, N_O={r['elements']['O']['N']:.4e}")

    r = compound_atom_density(1.104, "D2O")
    print(f"Ağır su D2O: M={r['M_formula']:.3f}, "
          f"N_D={r['elements']['D']['N']:.4e}, N_O={r['elements']['O']['N']:.4e}")
