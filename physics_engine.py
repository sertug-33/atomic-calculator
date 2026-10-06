import re
import math
import numpy as np
from database import ISOTOPES, ALIASES, U_TO_MEV, M_P, M_N, M_H1, M_E

# |Q| bu değerin altındaysa reaksiyon elastik kabul edilir (MeV)
Q_ELASTIC_TOL = 1e-6


def atomic_mass(iso_key: str) -> float:
    """Q-değeri hesabında kullanılacak ATOMİK kütle.
    Veritabanındaki kütleler atomik (elektronlu) olduğu için çıplak proton ("p")
    yerine 1H atom kütlesi kullanılır; iki taraftaki elektronlar böylece dengelenir.
    (Nötronun elektronu yoktur, kütlesi olduğu gibi kullanılır.)"""
    if iso_key == "p":
        return M_H1
    return ISOTOPES[iso_key]["mass"]


def binding_energy(iso_key: str):
    """Deneysel bağlanma enerjisi (AME2020 atomik kütlelerinden).
    B = [Z·M(1H) + N·m_n − M_atom] · 931.494 MeV
    Dikkat: Z·M_P (çıplak proton) KULLANILMAZ — Z elektron eksik kalır (~0.511·Z MeV hata)."""
    if iso_key not in ISOTOPES:
        return None
    iso = ISOTOPES[iso_key]
    Z, A, m = iso["Z"], iso["A"], iso["mass"]
    N = A - Z
    if A <= 1:
        return {"delta_m": 0.0, "Eb": 0.0, "Eb_per_A": 0.0}
    delta_m = Z * M_H1 + N * M_N - m
    Eb = delta_m * U_TO_MEV
    return {"delta_m": delta_m, "Eb": Eb, "Eb_per_A": Eb / A}


def semf_deviation(iso_key: str):
    """SEMF'in deneysel bağlanma enerjisinden yüzde sapması (işaretli)."""
    exp = binding_energy(iso_key)
    if not exp or exp["Eb"] == 0:
        return None
    iso = ISOTOPES[iso_key]
    semf = semf_binding_energy(iso["A"], iso["Z"])
    return (semf - exp["Eb"]) / exp["Eb"] * 100.0

# STANDART LİTERATÜR MODERATÖR VERİLERİ (Lamarsh & Baratta, Table 7.1 & 7.2)
MODERATORS = {
    "H2O": {
        "name": "Hafif Su (H₂O)",
        "formula": "H₂O",
        "A_eff": 1.0,
        "xi": 0.920,
        "sigma_s": 44.0,       # Epithermal / thermal scattering (barn)
        "Sigma_s": 1.48,       # cm^-1
        "Sigma_a": 0.022,      # cm^-1 (0.0253 eV termal yutulma)
        "msdp": 1.36,          # xi * Sigma_s (cm^-1)
        "mr": 62.0,            # Moderating Ratio
        "color": "#ffb454",
        "reactor_type": "PWR / BWR (Zenginleştirilmiş U şarttır; MR düşüktür)"
    },
    "D2O": {
        "name": "Ağır Su (D₂O)",
        "formula": "D₂O",
        "A_eff": 2.0,
        "xi": 0.509,
        "sigma_s": 10.6,
        "Sigma_s": 0.35,
        "Sigma_a": 0.000033,   # Neredeyse sıfır nötron yutumu
        "msdp": 0.178,
        "mr": 5400.0,
        "color": "#4fa8b2",
        "reactor_type": "CANDU / PHWR (Doğal Uranyum ile kritikliği sağlar)"
    },
    "Be": {
        "name": "Berilyum (Be-9)",
        "formula": "⁹Be",
        "A_eff": 9.0,
        "xi": 0.207,
        "sigma_s": 6.1,
        "Sigma_s": 0.65,
        "Sigma_a": 0.0012,
        "msdp": 0.135,
        "mr": 112.0,
        "color": "#86c98f",
        "reactor_type": "Araştırma ve Test Reaktörleri (Nötron yansıtıcı kalkan)"
    },
    "C": {
        "name": "Reaktör Grafiti (C-12)",
        "formula": "¹²C",
        "A_eff": 12.0,
        "xi": 0.158,
        "sigma_s": 4.8,
        "Sigma_s": 0.38,
        "Sigma_a": 0.00038,
        "msdp": 0.060,
        "mr": 158.0,
        "color": "#d6ddd4",
        "reactor_type": "RBMK / HTGR / Magnox (Yüksek sıcaklık gaz reaktörleri)"
    }
}

def semf_breakdown(A: int, Z: int):
    if A <= 0:
        return {"total": 0.0, "vol": 0.0, "surf": 0.0, "coul": 0.0, "asym": 0.0, "pair": 0.0}
    N = A - Z
    av, as_, ac, aa, ap = 15.8, 18.3, 0.714, 23.2, 12.0
    vol = av * A
    surf = -as_ * (A ** (2/3))
    coul = -ac * (Z * (Z - 1)) / (A ** (1/3)) if A > 1 else 0.0
    asym = -aa * ((A - 2 * Z) ** 2) / A
    
    if A % 2 != 0:
        pair = 0.0
    elif Z % 2 == 0 and N % 2 == 0:
        pair = ap / (A ** 0.5)
    else:
        pair = -ap / (A ** 0.5)
        
    total = max(0.0, vol + surf + coul + asym + pair)
    return {
        "total": total, "vol": vol, "surf": surf,
        "coul": coul, "asym": asym, "pair": pair
    }

def semf_binding_energy(A: int, Z: int) -> float:
    return semf_breakdown(A, Z)["total"]

def separation_energies(iso_key: str):
    if iso_key not in ISOTOPES:
        return None
    iso = ISOTOPES[iso_key]
    z, a, m = iso["Z"], iso["A"], iso["mass"]
    
    sn, sn_type = None, "Deneysel"
    if a > 1:
        daughter_n_mass = None
        for k, v in ISOTOPES.items():
            if v["Z"] == z and v["A"] == a - 1:
                daughter_n_mass = v["mass"]
                break
        if daughter_n_mass is not None:
            sn = (daughter_n_mass + M_N - m) * U_TO_MEV
        else:
            sn = semf_binding_energy(a, z) - semf_binding_energy(a - 1, z)
            sn_type = "SEMF"

    sp, sp_type = None, "Deneysel"
    if a > 1 and z >= 1:
        daughter_p_mass = None
        for k, v in ISOTOPES.items():
            if v["Z"] == z - 1 and v["A"] == a - 1:
                daughter_p_mass = v["mass"]
                break
        if daughter_p_mass is not None:
            # Atomik kütlelerle çalışıldığı için proton yerine 1H atomu (p + e)
            sp = (daughter_p_mass + M_H1 - m) * U_TO_MEV
        else:
            sp = semf_binding_energy(a, z) - semf_binding_energy(a - 1, z - 1)
            sp_type = "SEMF"

    return {"sn": sn, "sn_type": sn_type, "sp": sp, "sp_type": sp_type}

def resolve_alias(symbol: str) -> str:
    return ALIASES.get(symbol, symbol)

def parse_nuclear_item(item_str: str):
    item_str = item_str.strip()
    norm = resolve_alias(item_str)
    if norm in ISOTOPES:
        return 1, norm
    m = re.match(r"^(\d+)\s*\*?\s*([a-zA-Z0-9α\-_]+)$", item_str)
    if m:
        coeff_str, sym = m.groups()
        sym_norm = resolve_alias(sym)
        if sym_norm in ISOTOPES:
            return int(coeff_str), sym_norm
    raise ValueError(f"Bilinmeyen bileşen: '{item_str}'")

def calculate_reaction(text: str, projectile: str | None = None):
    """X + a -> Y + b reaksiyonu için Q-değeri ve eşik enerjisi.
    projectile: mermi çekirdeğin anahtarı (örn. "4He"). Verilmezse soldaki en hafif
    parçacık mermi kabul edilir (ters kinematikte yanlış olur — kanal seçiciden geçirin)."""
    if "->" not in text:
        raise ValueError("Denklemde '->' işareti bulunmalıdır.")

    def parse_side(side_str):
        parts = [p.strip() for p in side_str.split("+") if p.strip()]
        mass, tot_z, tot_a = 0.0, 0, 0
        items = []
        for item in parts:
            coeff, iso = parse_nuclear_item(item)
            m_iso = atomic_mass(iso)          # p -> 1H (elektron dengesi)
            mass += coeff * m_iso
            tot_z += coeff * ISOTOPES[iso]["Z"]
            tot_a += coeff * ISOTOPES[iso]["A"]
            items.extend([(iso, m_iso)] * coeff)
        return mass, tot_z, tot_a, items

    left, right = text.split("->")
    m_in, z_in, a_in, in_items = parse_side(left)
    m_out, z_out, a_out, _ = parse_side(right)
    delta_m = m_in - m_out
    q_val = delta_m * U_TO_MEV

    is_elastic = abs(q_val) < Q_ELASTIC_TOL
    if is_elastic:
        q_val = 0.0
        phase = "ELASTİK"
    elif q_val > 0:
        phase = "EKZOTERMİK"
    else:
        phase = "ENDOTERMİK"

    e_threshold = None          # tam (rölativistik) eşik, MeV
    e_threshold_approx = None   # ders kitabı yaklaşımı: -Q(1 + m_a/m_X)
    if phase == "ENDOTERMİK" and len(in_items) == 2:
        proj_key = resolve_alias(projectile) if projectile else None
        if proj_key is not None and proj_key in [k for k, _ in in_items]:
            idx = [k for k, _ in in_items].index(proj_key)
        else:
            idx = 0 if in_items[0][1] <= in_items[1][1] else 1
        m_a = in_items[idx][1]
        m_X = in_items[1 - idx][1]
        # Rölativistik tam eşik: T_th = -Q · (Σm_giren + Σm_çıkan) / (2·m_X)
        # (düşük enerjide ders kitabındaki -Q(1 + m_a/m_X) ile 4. haneye kadar aynıdır)
        e_threshold = -q_val * (m_in + m_out) / (2.0 * m_X)
        e_threshold_approx = -q_val * (1.0 + m_a / m_X)

    return {
        "m_in": m_in, "m_out": m_out, "delta_m": delta_m,
        "q_val": q_val, "z_in": z_in, "z_out": z_out,
        "a_in": a_in, "a_out": a_out,
        "z_conserved": z_in == z_out,
        "a_conserved": a_in == a_out,
        "is_exothermic": phase == "EKZOTERMİK",
        "is_elastic": is_elastic,
        "phase": phase,                 # "EKZOTERMİK" / "ENDOTERMİK" / "ELASTİK"
        "n_in": len(in_items),          # 1 ise bozunma, 2 ise X(a,b)Y reaksiyonu
        "e_threshold": e_threshold,
        "e_threshold_approx": e_threshold_approx,
    }

def solve_channel(target: str, projectile: str, emitted: str):
    t_norm = resolve_alias(target)
    p_norm = resolve_alias(projectile)
    e_norm = resolve_alias(emitted)

    if t_norm not in ISOTOPES or p_norm not in ISOTOPES or e_norm not in ISOTOPES:
        raise ValueError("Seçilen parçacıklardan biri veritabanında bulunamadı.")

    zt = ISOTOPES[t_norm]["Z"]
    at = ISOTOPES[t_norm]["A"]
    zp = ISOTOPES[p_norm]["Z"]
    ap = ISOTOPES[p_norm]["A"]
    ze = ISOTOPES[e_norm]["Z"]
    ae = ISOTOPES[e_norm]["A"]

    zy = zt + zp - ze
    ay = at + ap - ae

    if zy < 0 or ay < 0:
        return None, zy, ay, f"İmkânsız Kanal (Z={zy}, A={ay})"

    residual = None
    for k, v in ISOTOPES.items():
        if v["Z"] == zy and v["A"] == ay:
            residual = k
            break

    equation = f"{t_norm} + {p_norm} -> {residual if residual else 'X'} + {e_norm}"
    return residual, zy, ay, equation

def calculate_bateman(t_half_A: float, t_half_B: float = 25.0, N0: float = 1000.0):
    t = np.linspace(0, 100, 300)
    lam_A = np.log(2) / max(0.01, t_half_A)
    lam_B = np.log(2) / max(0.01, t_half_B)
    N_A = N0 * np.exp(-lam_A * t)
    if abs(lam_A - lam_B) > 1e-5:
        N_B = N0 * (lam_A / (lam_B - lam_A)) * (np.exp(-lam_A * t) - np.exp(-lam_B * t))
    else:
        N_B = N0 * lam_A * t * np.exp(-lam_A * t)
    N_C = np.maximum(0, N0 - N_A - N_B)
    return t, N_A, N_B, N_C

def calculate_moderation(mod_key: str, E0: float = 2.0e6, E_th: float = 0.0253):
    """
    Nötron yavaşlama kinetiği: E(k) = E0 * exp(-k * xi)
    k: Çarpışma indeksi
    E0: Başlangıç fisyon nötron enerjisi (2.0 MeV)
    E_th: Hedef termal enerji (0.0253 eV)
    """
    data = MODERATORS.get(mod_key, MODERATORS["H2O"])
    xi = data["xi"]
    delta_u = math.log(E0 / E_th)
    n_coll = delta_u / xi
    n_steps = int(math.ceil(n_coll)) + 1
    
    k_arr = np.linspace(0, n_coll, 150)
    E_arr = E0 * np.exp(-k_arr * xi)
    
    # Kinematik alfa parametresi
    A = data["A_eff"]
    alpha = ((A - 1.0) / (A + 1.0)) ** 2 if A > 1.0 else 0.0
    max_loss = (1.0 - alpha) * 100.0
    # Logaritmik ortalama kayıp: E'nin ortalama ln'i ξ kadar azalır -> 1 - e^(-ξ)
    avg_loss = (1.0 - math.exp(-xi)) * 100.0
    # Aritmetik ortalama kesirsel kayıp (tek çekirdek, izotropik CM saçılması): (1-α)/2
    avg_loss_arith = (1.0 - alpha) / 2.0 * 100.0

    return {
        "data": data,
        "delta_u": delta_u,
        "n_coll": n_coll,
        "alpha": alpha,
        "max_loss": max_loss,
        "avg_loss": avg_loss,            # etiket: "Log-Ortalama Enerji Kaybı"
        "avg_loss_arith": avg_loss_arith,  # etiket: "Aritmetik Ort. Enerji Kaybı"
        "k_arr": k_arr,
        "E_arr": E_arr
    }


# --- Bozunma Modları: kanal-bazlı Q-değeri ve spesifik aktivite --------------
# Element sembolü <-> Z eşlemesi (ISOTOPES'in "name" alanından türetilir,
# density.py'deki _Z_TO_SYMBOL ile aynı yöntem — modüller birbirinden bağımsız
# kalsın diye burada ayrıca inşa edilir).
_Z_TO_SYMBOL = {}
for _k, _v in ISOTOPES.items():
    if _k in ("n", "p"):
        continue
    _Z_TO_SYMBOL.setdefault(_v["Z"], _v["name"].split("-")[0])


def _iso_key_for(Z: int, A: int):
    sym = _Z_TO_SYMBOL.get(Z)
    return f"{A}{sym}" if sym else None


def decay_channels(iso_key: str):
    """Verilen nüklit için enerjik/geometrik olarak tanımlı bozunma kanallarının
    (α, β⁻, β⁺, EC) Q-değerlerini döndürür. AME2020 ATOMİK kütleleri kullanılır:
      Q_α  = [M(Z,A) − M(Z-2,A-4) − M(4He)]·931.494
      Q_β- = [M(Z,A) − M(Z+1,A)]·931.494                (elektronlar kendiliğinden dengeli)
      Q_β+ = [M(Z,A) − M(Z-1,A) − 2m_e]·931.494          (pozitron + fazla atomik elektron)
      Q_EC = [M(Z,A) − M(Z-1,A)]·931.494                 (yakalanan elektron zaten atomik kütlede)
    Kanal yalnızca ürün nüklit veritabanında mevcutsa hesaba katılır; Q>0 ise
    o kanal enerjik olarak MÜMKÜNDÜR (gerçekte gözlenip gözlenmediği ayrı konudur)."""
    if iso_key not in ISOTOPES:
        return None
    iso = ISOTOPES[iso_key]
    Z, A, M = iso["Z"], iso["A"], iso["mass"]
    out = {}

    if Z > 2 and A > 4:
        d_key = _iso_key_for(Z - 2, A - 4)
        if d_key and d_key in ISOTOPES:
            q = (M - ISOTOPES[d_key]["mass"] - ISOTOPES["4He"]["mass"]) * U_TO_MEV
            out["alpha"] = {"q": q, "daughter": d_key, "label": "α BOZUNMASI"}

    d_key = _iso_key_for(Z + 1, A)
    if d_key and d_key in ISOTOPES:
        q = (M - ISOTOPES[d_key]["mass"]) * U_TO_MEV
        out["beta_minus"] = {"q": q, "daughter": d_key, "label": "β⁻ BOZUNMASI"}

    if Z > 1:
        d_key = _iso_key_for(Z - 1, A)
        if d_key and d_key in ISOTOPES:
            q = (M - ISOTOPES[d_key]["mass"] - 2.0 * M_E) * U_TO_MEV
            out["beta_plus"] = {"q": q, "daughter": d_key, "label": "β⁺ BOZUNMASI"}

            q_ec = (M - ISOTOPES[d_key]["mass"]) * U_TO_MEV
            out["ec"] = {"q": q_ec, "daughter": d_key, "label": "ELEKTRON YAKALAMA (EC)"}

    return out


# --- Bozunma Serisi Zinciri: baskın moda göre kararlıya kadar izleme ---------
# Baskın bozunma modu -> (ΔZ, ΔA) net dönüşüm. Gecikmeli parçacık emisyonlu
# bileşik modlar (β⁻n, β⁻p, β⁻α, β⁺p, β⁺α) iki adımın net etkisi olarak girildi.
DECAY_TRANSITIONS = {
    "α": (-2, -4), "β-": (1, 0), "β+": (-1, 0), "EC": (-1, 0),
    "p": (-1, -1), "n": (0, -1),
    "β-n": (1, -1), "β-p": (0, -1), "β-α": (-1, -4),
    "β+p": (-2, -1), "β+α": (-3, -4),
    "2β-": (2, 0), "2β+": (-2, 0), "2EC": (-2, 0),
}

# Doğada bulunan / tarihsel olarak önemli dört bozunma serisi (kütle numarası mod 4'e göre).
DECAY_SERIES = {
    "238U":  {"name": "Uranyum Serisi (4n+2)",                 "final": "206Pb"},
    "235U":  {"name": "Aktinyum Serisi (4n+3)",                "final": "207Pb"},
    "232Th": {"name": "Toryum Serisi (4n)",                    "final": "208Pb"},
    "237Np": {"name": "Neptünyum Serisi (4n+1, soyu tükenmiş)", "final": "205Tl"},
}


def decay_series_chain(start_key: str, max_steps: int = 40):
    """Verilen nüklitten başlayıp NUBASE'deki BASKIN bozunma moduna göre (mode
    alanındaki '/' ile ayrılmış listenin ilk elemanı) kararlı bir nüklide
    ulaşana kadar zinciri izler. Her adımda o adımın Q-değerini de hesaplar.
    Dönüş: düğüm sözlüklerinden oluşan liste — her biri
      {key, name, Z, A, half_life, mode, stable, step_mode, q, terminal_reason}
    """
    if start_key not in ISOTOPES:
        return []
    chain = []
    visited = set()
    key = start_key
    for _ in range(max_steps):
        if key not in ISOTOPES or key in visited:
            break
        visited.add(key)
        iso = ISOTOPES[key]
        Z, A, M = iso["Z"], iso["A"], iso["mass"]
        stable = (iso["half_life"] == "Kararlı")
        mode_str = iso["mode"]
        dominant = mode_str.split("/")[0] if mode_str and mode_str != "-" else None
        node = {
            "key": key, "name": iso["name"], "Z": Z, "A": A,
            "half_life": iso["half_life"], "mode": mode_str,
            "stable": stable, "step_mode": dominant, "q": None,
            "terminal_reason": None,
        }
        chain.append(node)
        if stable:
            break
        if dominant not in DECAY_TRANSITIONS:
            node["terminal_reason"] = "Bilinmeyen / karmaşık bozunma kanalı (zincir burada kesiliyor)"
            break
        dZ, dA = DECAY_TRANSITIONS[dominant]
        new_Z, new_A = Z + dZ, A + dA
        if new_Z < 0 or new_A < 1:
            node["terminal_reason"] = "Geçersiz ürün (Z<0 veya A<1)"
            break
        new_key = _iso_key_for(new_Z, new_A)
        if not new_key or new_key not in ISOTOPES:
            node["terminal_reason"] = "Ürün nüklit veritabanında yok"
            break
        if dominant == "α":
            node["q"] = (M - ISOTOPES[new_key]["mass"] - ISOTOPES["4He"]["mass"]) * U_TO_MEV
        elif dominant == "β-":
            node["q"] = (M - ISOTOPES[new_key]["mass"]) * U_TO_MEV
        elif dominant in ("β+",):
            node["q"] = (M - ISOTOPES[new_key]["mass"] - 2.0 * M_E) * U_TO_MEV
        elif dominant == "EC":
            node["q"] = (M - ISOTOPES[new_key]["mass"]) * U_TO_MEV
        else:
            node["q"] = None  # bileşik/nadir kanallar için Q burada hesaplanmaz
        key = new_key
    return chain


# --- Gama Zayıflama / Zırhlama -----------------------------------------------
# Kütle zayıflama katsayıları (μ/ρ, cm²/g) NIST XCOM/X-Ray Mass Attenuation
# Coefficients referans tablolarından alınmıştır (fotoelektrik+Compton+çift
# oluşum toplamı, koherent saçılma dahil). Yoğunluklar (g/cm³) standart
# literatür değerleridir. Beton için NIST "Ordinary Concrete" bileşimi kullanıldı.
SHIELDING_MATERIALS = {
    "Pb": {
        "name": "Kurşun (Pb)", "density": 11.35, "color": "#8a8f98",
        "mu_rho": {0.1: 5.549, 0.15: 2.014, 0.2: 0.9985, 0.3: 0.4031, 0.4: 0.2323,
                   0.5: 0.1614, 0.6: 0.1248, 0.8: 0.0887, 1.0: 0.07102, 1.5: 0.05222,
                   2.0: 0.04606, 3.0: 0.04234, 4.0: 0.04197, 5.0: 0.04272, 6.0: 0.04391,
                   8.0: 0.04675, 10.0: 0.04972},
    },
    "Fe": {
        "name": "Demir / Çelik (Fe)", "density": 7.874, "color": "#b3752c",
        "mu_rho": {0.1: 0.3717, 0.15: 0.1964, 0.2: 0.1460, 0.3: 0.1099, 0.4: 0.0940,
                   0.5: 0.08414, 0.6: 0.07704, 0.8: 0.06699, 1.0: 0.05995, 1.5: 0.04883,
                   2.0: 0.04265, 3.0: 0.03621, 4.0: 0.03312, 5.0: 0.03146, 6.0: 0.03057,
                   8.0: 0.02991, 10.0: 0.02994},
    },
    "Concrete": {
        "name": "Beton (Ordinary Concrete)", "density": 2.3, "color": "#9a9484",
        "mu_rho": {0.1: 0.1738, 0.15: 0.1436, 0.2: 0.1282, 0.3: 0.1097, 0.4: 0.09783,
                   0.5: 0.08915, 0.6: 0.08236, 0.8: 0.07227, 1.0: 0.06495, 1.5: 0.05288,
                   2.0: 0.04557, 3.0: 0.03701, 4.0: 0.03217, 5.0: 0.02908, 6.0: 0.02697,
                   8.0: 0.02432, 10.0: 0.02278},
    },
    "Al": {
        "name": "Alüminyum (Al)", "density": 2.699, "color": "#c9ccd1",
        "mu_rho": {0.1: 0.1704, 0.15: 0.1378, 0.2: 0.1223, 0.3: 0.1042, 0.4: 0.09276,
                   0.5: 0.08445, 0.6: 0.07802, 0.8: 0.06841, 1.0: 0.06146, 1.5: 0.05006,
                   2.0: 0.04324, 3.0: 0.03541, 4.0: 0.03106, 5.0: 0.02836, 6.0: 0.02655,
                   8.0: 0.02437, 10.0: 0.02318},
    },
    "H2O": {
        "name": "Su (H₂O)", "density": 1.0, "color": "#4fa8b2",
        "mu_rho": {0.1: 0.1707, 0.15: 0.1505, 0.2: 0.1370, 0.3: 0.1186, 0.4: 0.1061,
                   0.5: 0.09687, 0.6: 0.08956, 0.8: 0.07865, 1.0: 0.07072, 1.5: 0.05754,
                   2.0: 0.04942, 3.0: 0.03969, 4.0: 0.03403, 5.0: 0.03031, 6.0: 0.02770,
                   8.0: 0.02429, 10.0: 0.02219},
    },
}


def gamma_mu_rho(material_key: str, energy_MeV: float) -> float:
    """Verilen enerjide kütle zayıflama katsayısını (cm²/g) tablo noktaları
    arasında LOG-LOG enterpolasyonla döndürür (μ/ρ birkaç büyüklük mertebesi
    değiştiği için log-log, doğrusal enterpolasyondan çok daha isabetlidir).
    Tablo aralığının dışında en yakın uç değere sabitlenir (ekstrapolasyon yok)."""
    mat = SHIELDING_MATERIALS[material_key]
    energies = sorted(mat["mu_rho"].keys())
    if energy_MeV <= energies[0]:
        return mat["mu_rho"][energies[0]]
    if energy_MeV >= energies[-1]:
        return mat["mu_rho"][energies[-1]]
    log_e = np.log(energies)
    log_mu = np.log([mat["mu_rho"][e] for e in energies])
    return float(np.exp(np.interp(np.log(energy_MeV), log_e, log_mu)))


# Taylor iki-üstel buildup faktörü katsayıları: B(μx) = A1·e^(α1·μx) + (1−A1)·e^(α2·μx)
# (μx = kalınlığın ortalama serbest yol (mfp) cinsinden değeri, boyutsuz).
# Kaynak: Univ. of Illinois NPRE 441 "Principles of Radiation Protection" ders notları
# (M. Ragheb), Table III — nokta izotropik kaynak, üç bağımsız sorguda birebir tutarlı
# doğrulandı. YALNIZCA bu üç malzeme ve 1-10 MeV aralığı için buildup uygulanır; tablo
# aralığı dışında veya Pb/Al için güvenilir/doğrulanmış katsayı bulunamadığından buildup
# hesaba katılmaz (dar-demet sonucu döner) — uydurma katsayıyla yanlış sonuç vermektense
# "veri yok" demeyi tercih ediyoruz.
BUILDUP_TAYLOR = {
    "H2O": {1.0: (11.0, -0.104, 0.030), 2.0: (6.4, -0.076, 0.092), 3.0: (5.2, -0.062, 0.110),
            4.0: (4.5, -0.055, 0.117), 6.0: (3.55, -0.050, 0.124), 8.0: (3.05, -0.045, 0.128),
            10.0: (2.7, -0.042, 0.130)},
    "Concrete": {1.0: (10.0, -0.088, 0.029), 2.0: (6.3, -0.069, 0.058), 3.0: (4.7, -0.062, 0.073),
                 4.0: (3.9, -0.059, 0.079), 6.0: (3.1, -0.059, 0.083), 8.0: (2.7, -0.056, 0.086),
                 10.0: (2.6, -0.050, 0.084)},
    "Fe": {1.0: (8.0, -0.089, 0.040), 2.0: (5.5, -0.079, 0.070), 3.0: (4.4, -0.077, 0.075),
           4.0: (3.75, -0.075, 0.082), 6.0: (2.9, -0.082, 0.075), 8.0: (2.35, -0.083, 0.055),
           10.0: (2.0, -0.095, 0.012)},
}


def buildup_factor(material_key: str, energy_MeV: float, mu_x: float):
    """Taylor formuyla geniş-demet (broad-beam) buildup faktörü B(μx):
        B(μx) = (1−A1)·e^(α1·μx) + A1·e^(α2·μx)
    (A1 büyük katsayı YAVAŞ BÜYÜYEN üstel terimle (α2>0) eşleşir, küçük/negatif
    katsayı HIZLI SÖNEN terimle (α1<0) eşleşir — bu eşleşme B(0)=1 sınır koşulunu
    sağlayan VE dört olası işaret/eşleşme kombinasyonu arasında tekil biçimde hem
    pozitif hem monoton artan hem de mfp arttıkça makul mertebede kalan sonucu
    veren kombinasyondur; ters eşleşme (A1 ↔ α1) negatif/anlamsız B üretir ve
    elendi. Referans: Univ. of Illinois NPRE 441 ders notları, Tablo III.)
    Doğrulanmış katsayı tablosu olmayan malzeme/enerji kombinasyonları için None
    döner (çağıran taraf bunu 'buildup hesaplanamadı, dar-demet sonucu' olarak
    yorumlamalı — B=1 varsaymak SESSİZCE yanlış sonuç üretir, bu yüzden None tercih
    edildi)."""
    table = BUILDUP_TAYLOR.get(material_key)
    if not table:
        return None
    energies = sorted(table.keys())
    if energy_MeV < energies[0] or energy_MeV > energies[-1]:
        return None
    A1 = np.interp(energy_MeV, energies, [table[e][0] for e in energies])
    a1 = np.interp(energy_MeV, energies, [table[e][1] for e in energies])
    a2 = np.interp(energy_MeV, energies, [table[e][2] for e in energies])
    B = (1.0 - A1) * math.exp(a1 * mu_x) + A1 * math.exp(a2 * mu_x)
    return max(float(B), 1.0)


def gamma_attenuation(material_key: str, energy_MeV: float, thickness_cm: float, I0: float = 100.0):
    """Lambert-Beer zayıflama yasası: I = I0 · e^(−μx) (dar-demet/narrow-beam).
    Dönüş: μ/ρ, μ (lineer, cm⁻¹), I, geçirgenlik oranı, HVL (yarı-değer kalınlığı,
    ln2/μ), TVL (onda-bir değer kalınlığı, ln10/μ) ve — veri varsa — Taylor buildup
    faktörüyle düzeltilmiş geniş-demet (broad-beam) şiddeti I_broad = I·B(μx)."""
    mat = SHIELDING_MATERIALS[material_key]
    mu_rho = gamma_mu_rho(material_key, energy_MeV)
    rho = mat["density"]
    mu = mu_rho * rho
    mu_x = mu * thickness_cm
    I = I0 * math.exp(-mu_x)
    hvl = math.log(2) / mu if mu > 0 else float("inf")
    tvl = math.log(10) / mu if mu > 0 else float("inf")
    B = buildup_factor(material_key, energy_MeV, mu_x)
    I_broad = I * B if B is not None else None
    return {
        "material": mat["name"], "material_key": material_key,
        "energy": energy_MeV, "thickness": thickness_cm,
        "mu_rho": mu_rho, "mu": mu, "rho": rho, "mu_x": mu_x,
        "I0": I0, "I": I,
        "transmission": (I / I0) if I0 else 0.0,
        "attenuation_pct": (1.0 - I / I0) * 100.0 if I0 else 0.0,
        "hvl": hvl, "tvl": tvl,
        "n_hvl": (thickness_cm / hvl) if hvl not in (0.0, float("inf")) else 0.0,
        "B": B, "I_broad": I_broad,
        "transmission_broad": (I_broad / I0) if (I0 and I_broad is not None) else None,
    }


# --- Radyoaktif Kaynak Spektrumları (çoklu gama çizgisi) ---------------------
# Enerjiler (MeV) ve emisyon olasılıkları (bozunma başına foton, kesir) IAEA
# "Xgamma Reference Gamma-Ray Standards" veri tabanından alınmıştır
# (https://www-nds.iaea.org/xgamma_standards/). Sadece ≥%1 şiddetli ana çizgiler
# dahil edilmiştir. Na-22'de 0.511 MeV çizgisinin olasılığı 1'i aşar (pozitron
# başına 2 foton yayınlanır) — bu fiziksel olarak doğrudur, hata değildir.
RADIOACTIVE_SOURCES = {
    "Co60": {"name": "Kobalt-60 (⁶⁰Co)", "lines": [
        {"energy": 1.173228, "intensity": 0.9985},
        {"energy": 1.332492, "intensity": 0.99983},
    ]},
    "Cs137": {"name": "Sezyum-137 (¹³⁷Cs)", "lines": [
        {"energy": 0.661657, "intensity": 0.8499},
    ]},
    "Na22": {"name": "Sodyum-22 (²²Na)", "lines": [
        {"energy": 0.511000, "intensity": 1.798},
        {"energy": 1.274537, "intensity": 0.99940},
    ]},
    "I131": {"name": "İyot-131 (¹³¹I)", "lines": [
        {"energy": 0.080185, "intensity": 0.02607},
        {"energy": 0.284305, "intensity": 0.0606},
        {"energy": 0.364489, "intensity": 0.812},
        {"energy": 0.636989, "intensity": 0.0726},
    ]},
    "Am241": {"name": "Amerikyum-241 (²⁴¹Am)", "lines": [
        {"energy": 0.026345, "intensity": 0.024},
        {"energy": 0.059541, "intensity": 0.3578},
    ]},
    "Ba133": {"name": "Baryum-133 (¹³³Ba)", "lines": [
        {"energy": 0.080998, "intensity": 0.329},
        {"energy": 0.302851, "intensity": 0.1834},
        {"energy": 0.356013, "intensity": 0.6205},
    ]},
}


def gamma_attenuation_multiline(material_key: str, source_key: str, thickness_cm: float, I0_total: float = 100.0):
    """Çok çizgili bir radyoaktif kaynağın her gama çizgisi için AYRI AYRI
    zayıflama hesaplanır (her çizginin μ'sü farklıdır), sonra çizgiler bozunum
    başına foton olasılıklarıyla ağırlıklandırılarak BİRLEŞİK bir geçirgenlik
    elde edilir. Bu, gerçek bir dedektörün/dozimetrenin ölçtüğü toplam şiddete
    karşılık gelir (basit toplamsal model; çizgiler arası girişim/kaskad açısal
    korelasyonu ihmal edilir)."""
    source = RADIOACTIVE_SOURCES[source_key]
    total_intensity = sum(l["intensity"] for l in source["lines"])
    lines_out = []
    I0_sum = I_sum = 0.0
    I_broad_sum = 0.0
    any_buildup = False
    for line in source["lines"]:
        frac = line["intensity"] / total_intensity
        r = gamma_attenuation(material_key, line["energy"], thickness_cm, I0=I0_total * frac)
        r["intensity_frac"] = frac
        lines_out.append(r)
        I0_sum += r["I0"]
        I_sum += r["I"]
        if r["I_broad"] is not None:
            I_broad_sum += r["I_broad"]
            any_buildup = True
        else:
            I_broad_sum += r["I"]

    return {
        "material": SHIELDING_MATERIALS[material_key]["name"], "material_key": material_key,
        "source": source["name"], "source_key": source_key, "thickness": thickness_cm,
        "lines": lines_out,
        "I0": I0_sum, "I": I_sum,
        "transmission": (I_sum / I0_sum) if I0_sum else 0.0,
        "attenuation_pct": (1.0 - I_sum / I0_sum) * 100.0 if I0_sum else 0.0,
        "I_broad": I_broad_sum if any_buildup else None,
        "transmission_broad": (I_broad_sum / I0_sum) if (I0_sum and any_buildup) else None,
        "has_buildup": any_buildup,
    }


def gamma_thickness_for_factor(material_key: str, energy_MeV: float, reduction_factor: float):
    """Verilen bir azaltma faktörüne (I0/I) ulaşmak için gereken kalınlığı (cm)
    döndürür. reduction_factor <= 1 için tanımsızdır (azaltma yok/artış)."""
    if reduction_factor <= 1:
        return None
    mat = SHIELDING_MATERIALS[material_key]
    mu = gamma_mu_rho(material_key, energy_MeV) * mat["density"]
    if mu <= 0:
        return None
    return math.log(reduction_factor) / mu


# --- Doz Hızı Hesaplayıcısı (Nokta-Çekirdek Yöntemi) -------------------------
# Havanın kütle enerji-soğurma katsayısı (μₑₙ/ρ, cm²/g) — NIST "X-Ray Mass
# Attenuation Coefficients" tablosundan (dry air, near sea level composition):
# https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/air.html
# DİKKAT: Bu, zırhlama modülünde kullanılan μ/ρ (kütle ZAYIFLAMA katsayısı)
# İLE AYNI ŞEY DEĞİLDİR — μₑₙ/ρ, foton enerjisinin ortama ne kadarının
# GERÇEKTEN SOĞURULUP dozu oluşturduğunu verir (Compton'da saçılan foton
# enerjisi soğurulmaz, bu yüzden μₑₙ/ρ < μ/ρ).
AIR_MU_EN_RHO = {
    0.02: 0.5389, 0.03: 0.1537, 0.04: 0.06833, 0.05: 0.04098, 0.06: 0.03041,
    0.08: 0.02407, 0.1: 0.02325, 0.15: 0.02496, 0.2: 0.02672, 0.3: 0.02872,
    0.4: 0.02949, 0.5: 0.02966, 0.6: 0.02953, 0.8: 0.02882, 1.0: 0.02789,
    1.25: 0.02666, 1.5: 0.02547, 2.0: 0.02345, 3.0: 0.02057, 4.0: 0.01870,
    5.0: 0.01740, 6.0: 0.01647, 8.0: 0.01525, 10.0: 0.01450,
}

# MeV/g'dan Gy/kg'a (= Gy) dönüşüm sabiti: 1 MeV = 1.602176634e-13 J,
# 1 g = 1e-3 kg  =>  1 MeV/g = 1.602176634e-13 / 1e-3 J/kg = 1.602176634e-10 Gy
MEV_PER_G_TO_GY = 1.602176634e-10


def air_mu_en_rho(energy_MeV: float) -> float:
    """Verilen enerjide havanın kütle enerji-soğurma katsayısını (cm²/g)
    LOG-LOG enterpolasyonla döndürür (gamma_mu_rho ile aynı yöntem)."""
    energies = sorted(AIR_MU_EN_RHO.keys())
    if energy_MeV <= energies[0]:
        return AIR_MU_EN_RHO[energies[0]]
    if energy_MeV >= energies[-1]:
        return AIR_MU_EN_RHO[energies[-1]]
    log_e = np.log(energies)
    log_mu = np.log([AIR_MU_EN_RHO[e] for e in energies])
    return float(np.exp(np.interp(np.log(energy_MeV), log_e, log_mu)))


def dose_rate_point_source(activity_Bq: float, energy_MeV: float, emission_prob: float,
                            distance_cm: float, shield_material: str | None = None,
                            shield_thickness_cm: float = 0.0):
    """Nokta-çekirdek (point-kernel) yöntemiyle serbest-alan (havada) hava
    kerma hızı. Foton eşdeğeri: fluence oranı φ = (A·p)/(4πd²)  [foton/(cm²·s)]
    Kerma hızı: K̇ = φ · E · (μₑₙ/ρ)_hava · 1.602176634e-10   [Gy/s]
    Fotonlar için w_R=1 kabul edilerek eşdeğer doz hızı (Sv/s) sayısal olarak
    aynı alınır (ICRP basitleştirmesi — nötron/alfa için GEÇERLİ DEĞİLDİR).
    Opsiyonel zırhlama verilirse önce gamma_attenuation ile (buildup varsa
    geniş-demet, yoksa dar-demet) zayıflatılmış fluence kullanılır."""
    if distance_cm <= 0:
        return None
    phi0 = (activity_Bq * emission_prob) / (4.0 * math.pi * distance_cm ** 2)  # foton/(cm^2 s), zırhsız

    shield_info = None
    transmission = 1.0
    if shield_material and shield_thickness_cm > 0:
        shield_info = gamma_attenuation(shield_material, energy_MeV, shield_thickness_cm, I0=100.0)
        if shield_info["I_broad"] is not None:
            transmission = shield_info["transmission_broad"]
        else:
            transmission = shield_info["transmission"]

    phi = phi0 * transmission
    mu_en_rho = air_mu_en_rho(energy_MeV)
    dose_rate_Gy_s = phi * energy_MeV * mu_en_rho * MEV_PER_G_TO_GY
    dose_rate_Gy_h = dose_rate_Gy_s * 3600.0
    return {
        "energy": energy_MeV, "emission_prob": emission_prob, "distance": distance_cm,
        "phi0": phi0, "phi": phi, "mu_en_rho": mu_en_rho,
        "transmission": transmission, "shield": shield_info,
        "dose_rate_Gy_s": dose_rate_Gy_s, "dose_rate_Gy_h": dose_rate_Gy_h,
        "dose_rate_mGy_h": dose_rate_Gy_h * 1000.0,
        "dose_rate_uSv_h": dose_rate_Gy_h * 1e6,   # foton için Sv sayısal olarak Gy'ye eşit
        "dose_rate_mrem_h": dose_rate_Gy_h * 100000.0,
    }


def dose_rate_multiline(activity_Bq: float, source_key: str, distance_cm: float,
                         shield_material: str | None = None, shield_thickness_cm: float = 0.0):
    """Çok çizgili bir radyoaktif kaynağın (RADIOACTIVE_SOURCES) her gama
    çizgisi için AYRI AYRI doz hızı hesaplanır (her çizginin hem μₑₙ/ρ'si hem
    zırh varsa zayıflaması farklıdır), sonra toplam doz hızı elde etmek için
    TOPLANIR (dozlar toplamsaldır — attenuation'ın aksine ağırlıklı ortalama
    değil, doğrudan toplam)."""
    source = RADIOACTIVE_SOURCES[source_key]
    lines_out = []
    total_Gy_s = 0.0
    for line in source["lines"]:
        r = dose_rate_point_source(activity_Bq, line["energy"], line["intensity"],
                                    distance_cm, shield_material, shield_thickness_cm)
        lines_out.append(r)
        total_Gy_s += r["dose_rate_Gy_s"]
    total_Gy_h = total_Gy_s * 3600.0
    return {
        "source": source["name"], "source_key": source_key, "distance": distance_cm,
        "lines": lines_out,
        "dose_rate_Gy_s": total_Gy_s, "dose_rate_Gy_h": total_Gy_h,
        "dose_rate_mGy_h": total_Gy_h * 1000.0,
        "dose_rate_uSv_h": total_Gy_h * 1e6,
        "dose_rate_mrem_h": total_Gy_h * 100000.0,
    }


def specific_activity(iso_key: str):
    """Spesifik aktivite: A_s = λ·N_A / M   (Bq/g ve Ci/g).
    Kararlı nüklitler ve yarı ömrü tanımsız kayıtlar için None döner."""
    if iso_key not in ISOTOPES:
        return None
    iso = ISOTOPES[iso_key]
    t_half_s = iso.get("half_life_s")
    M = iso["mass"]
    if not t_half_s or t_half_s == float("inf") or t_half_s <= 0:
        return None
    N_A = 6.02214076e23
    lam = math.log(2) / t_half_s
    bq_per_g = lam * N_A / M
    return {
        "lambda": lam, "t_half_s": t_half_s,
        "bq_per_g": bq_per_g, "ci_per_g": bq_per_g / 3.7e10,
    }


# --- Kritiklik / Dört-Faktör Formülü (k∞ = η·f·p·ε) --------------------------
# 2200 m/s (0.0253 eV) termal mikroskopik tesir kesitleri — iki bağımsız
# kaynakla (nuclear-power.com; Lamarsh-türevi ders materyali) çapraz
# doğrulandı, iç tutarlılık kontrolü: σ_f + σ_γ ≈ σ_a (582.2+99≈681 ✓).
U235_SIGMA_A = 680.8   # barn, toplam termal soğurma (fisyon+radyatif yakalama)
U235_SIGMA_F = 582.2   # barn, termal fisyon
U235_NU = 2.42          # fisyon başına ortalama yayınlanan nötron sayısı
U238_SIGMA_A = 2.70     # barn, termal soğurma (termal fisyon ihmal edilebilir düzeyde küçük)

# Saf moderatör için molar kütle (g/mol) ve yoğunluk (g/cm³) — moderatörün
# mikroskopik soğurma tesir kesitini (σa, barn), MODERATORS sözlüğündeki
# ZATEN DOĞRULANMIŞ makroskopik Σa (Lamarsh & Baratta Tablo 7.1/7.2)
# değerinden N=ρ·N_A/M yoluyla TÜRETMEK için kullanılır — yeni/doğrulanmamış
# bir literatür sayısı eklemek yerine mevcut, test edilmiş veriyle iç
# tutarlılık tercih edildi.
_MODERATOR_MOLAR_MASS = {"H2O": 18.015, "D2O": 20.028, "Be": 9.012, "C": 12.011}
_MODERATOR_DENSITY = {"H2O": 1.0, "D2O": 1.104, "Be": 1.848, "C": 1.60}
_N_A = 6.02214076e23


def eta_factor(enrichment_pct: float):
    """Üretim (reprodüksiyon) faktörü η = ν·Σf,yakıt / Σa,yakıt.
    enrichment_pct: yakıttaki U-235 ATOM yüzdesi (doğal U için ≈0.711).
    U-238'in termal fisyonu ihmal edilir (fiziksel olarak doğru yaklaşım).
    Doğrulama: e=%100 için η≈2.08, e=%0.711 (doğal U) için η≈1.34 —
    literatürdeki (nuclear-power.com) bilinen asimptotik/doğal-U değerleriyle
    birebir örtüşüyor."""
    if enrichment_pct is None or not (0.0 < enrichment_pct <= 100.0):
        return None
    e = enrichment_pct / 100.0
    sigma_a_fuel = e * U235_SIGMA_A + (1.0 - e) * U238_SIGMA_A
    sigma_f_fuel = e * U235_SIGMA_F
    eta = U235_NU * sigma_f_fuel / sigma_a_fuel
    return {"eta": eta, "sigma_a_fuel": sigma_a_fuel, "sigma_f_fuel": sigma_f_fuel}


def moderator_sigma_a_micro_barn(mod_key: str):
    """Moderatörün mikroskopik soğurma tesir kesiti (barn) — bkz. yukarıdaki
    modül notu (MODERATORS'taki doğrulanmış Σa'dan türetilir)."""
    if mod_key not in _MODERATOR_MOLAR_MASS:
        return None
    mat = MODERATORS[mod_key]
    M = _MODERATOR_MOLAR_MASS[mod_key]
    rho = _MODERATOR_DENSITY[mod_key]
    N = rho * _N_A / M                       # molekül veya atom / cm³
    sigma_a_cm2 = mat["Sigma_a"] / N
    return sigma_a_cm2 * 1e24                # cm² -> barn


def moderator_atom_density(mod_key: str):
    """Saf moderatörün (normal yoğunlukta) atom/molekül yoğunluğu (1/cm³)."""
    if mod_key not in _MODERATOR_MOLAR_MASS:
        return None
    M = _MODERATOR_MOLAR_MASS[mod_key]
    rho = _MODERATOR_DENSITY[mod_key]
    return rho * _N_A / M


def thermal_utilization_factor(mod_key: str, enrichment_pct: float, NF_over_NM: float):
    """Termal kullanım faktörü:  f = 1 / [1 + (N_M/N_F)·(σa,M/σa,F)]
    NF_over_NM: karışımdaki yakıt/moderatör ATOM yoğunluğu oranı (N_F/N_M,
    boyutsuz) — reaktör tasarımına bağlı bir girdi parametresidir."""
    eta_data = eta_factor(enrichment_pct)
    sigma_a_M = moderator_sigma_a_micro_barn(mod_key)
    if eta_data is None or sigma_a_M is None or not NF_over_NM or NF_over_NM <= 0:
        return None
    sigma_a_F = eta_data["sigma_a_fuel"]
    NM_over_NF = 1.0 / NF_over_NM
    f = 1.0 / (1.0 + NM_over_NF * (sigma_a_M / sigma_a_F))
    return {"f": f, "sigma_a_M": sigma_a_M, "sigma_a_F": sigma_a_F}


def resonance_escape_probability(mod_key: str, enrichment_pct: float,
                                  NF_over_NM: float, I_eff_barn: float):
    """Rezonans kaçma olasılığı, homojen karışım için GENEL/doğrulanmış
    üstel form ile:  p = exp[ −N₂₈·I_eff / (ξΣs)_moderatör ]
    (Lamarsh Böl. 6; Univ. of Illinois NPRE ders notları [mragheb.com];
    CANDU Course 227 ile bağımsız kaynaklarda tutarlı — genel form).

    ÖNEMLİ DÜRÜSTLÜK NOTU: Homojen karışımlar için I_eff'i doğrudan
    N_moderatör/N_U238 oranından üreten, İKİ BAĞIMSIZ GÜVENİLİR KAYNAKLA
    doğrulanmış bir ampirik korelasyon bulunamadı (yalnız çubuk/lump
    geometrisi için farklı bir korelasyon — I_eff=a+b/d — bulundu, o da bu
    homojen modele uygulanamaz). Bu yüzden I_eff BURADA UYDURULMAZ; ders
    kitabınızdan veya verilen bir değerden I_EFF_BARN PARAMETRESİ OLARAK
    GİRİLMELİDİR. Üstel formülün kendisi doğrulanmıştır, yalnızca I_eff'in
    kaynağı kullanıcıya aittir."""
    N_M = moderator_atom_density(mod_key)
    if N_M is None or not NF_over_NM or NF_over_NM <= 0 or I_eff_barn is None or I_eff_barn < 0:
        return None
    if enrichment_pct is None or not (0.0 < enrichment_pct <= 100.0):
        return None
    e = enrichment_pct / 100.0
    N_F = NF_over_NM * N_M
    N_28 = N_F * (1.0 - e)                    # yakıttaki U-238 atom yoğunluğu
    xi_sigma_s = MODERATORS[mod_key]["msdp"]   # ξΣs (cm⁻¹), zaten doğrulanmış tabloda mevcut
    if xi_sigma_s <= 0:
        return None
    exponent = -(N_28 * I_eff_barn * 1e-24) / xi_sigma_s
    p = math.exp(exponent)
    return {"p": p, "N_28": N_28, "xi_sigma_s": xi_sigma_s}


def k_infinity(mod_key: str, enrichment_pct: float, NF_over_NM: float,
               I_eff_barn: float, epsilon: float = 1.0):
    """Sonsuz ortam çoğalma faktörü:  k∞ = η · f · p · ε
    ε (hızlı fisyon faktörü): homojen karışımlarda ε≈1.00 (literatür:
    nuclear-power.com) — heterojen/topaklı yakıt tasarımlarında (~1.02-1.03)
    fiziksel mekanizma (hızlı nötronların moderasyondan önce U-238'e
    çarpması) farklı çalışır. Bu yüzden genel/homojen bir kapalı-form
    formülü yerine kullanıcı girdisi olarak bırakıldı (varsayılan 1.00).
    NOT: k_eff = k∞·P_sızıntısız hesaplanmaz (v1 kapsamı dışı) — kaçak
    olasılığı geometriye (B², M²=L²+τ) bağlıdır ve ayrı bir modül gerektirir."""
    eta_data = eta_factor(enrichment_pct)
    f_data = thermal_utilization_factor(mod_key, enrichment_pct, NF_over_NM)
    p_data = resonance_escape_probability(mod_key, enrichment_pct, NF_over_NM, I_eff_barn)
    if eta_data is None or f_data is None or p_data is None or epsilon is None or epsilon <= 0:
        return None
    eta, f, p = eta_data["eta"], f_data["f"], p_data["p"]
    k_inf = eta * f * p * epsilon
    return {
        "eta": eta, "f": f, "p": p, "epsilon": epsilon, "k_inf": k_inf,
        "sigma_a_fuel": eta_data["sigma_a_fuel"], "sigma_f_fuel": eta_data["sigma_f_fuel"],
        "sigma_a_M": f_data["sigma_a_M"], "N_28": p_data["N_28"], "xi_sigma_s": p_data["xi_sigma_s"],
    }
