"""
NUBASE2020 / AME2020 tam veritabanı yükleyicisi.

Kaynak: F.G. Kondev, M. Wang, W.J. Huang, S. Naimi, G. Audi,
        "The NUBASE2020 evaluation of nuclear physics properties",
        Chin. Phys. C45, 030001 (2021).

`nubase2020.txt` dosyasını (bu modülle aynı klasörde) satır satır okuyup
database.py'nin beklediği ISOTOPES şemasına (Z, A, mass, half_life, mode, name)
dönüştürür. Ayrıca element başına doğal izotopik bolluk tablosunu (ABUNDANCE)
üretir — atom yoğunluğu / doğal bulunurluk modülü bunu kullanacak.

Sadece taban durumları (isomer değil) okunur; NUBASE'de i=0 taban durumunu
gösterir (i=1,2,... izomerler bu sürümde ihmal edilir).
"""
import os
import re

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_DATA_FILE = os.path.join(_THIS_DIR, "nubase2020.txt")

U_TO_MEV = 931.49410242  # 1 u = 931.494 MeV/c^2 (CODATA 2018)

# NUBASE yarı ömür birim kodu -> (çarpan[saniye], Türkçe birim)
_UNIT_SECONDS = {
    "ys": 1e-24, "zs": 1e-21, "as": 1e-18, "fs": 1e-15, "ps": 1e-12,
    "ns": 1e-9, "us": 1e-6, "ms": 1e-3, "s": 1.0,
    "m": 60.0, "h": 3600.0, "d": 86400.0,
    "y": 3.1556952e7, "ky": 3.1556952e10, "My": 3.1556952e13,
    "Gy": 3.1556952e16, "Ty": 3.1556952e19, "Py": 3.1556952e22,
    "Ey": 3.1556952e25, "Zy": 3.1556952e28, "Yy": 3.1556952e31,
}

# Baskın bozunma etiketleri -> okunabilir Türkçe kısaltma
_MODE_LABELS = [
    ("B-n", "β-n"), ("B-A", "β-α"), ("B-p", "β-p"),
    ("B-", "β-"), ("B+p", "β+p"), ("B+A", "β+α"), ("B+", "β+"),
    ("EC", "EC"), ("2EC", "2EC"), ("2B-", "2β-"), ("2B+", "2β+"),
    ("IT", "IT"), ("A", "α"), ("SF", "SF"), ("p", "p"), ("n", "n"),
]


def _parse_mode(br_field: str) -> str:
    """'B-=100' ; 'B+=61.52 26;B-=38.48 26' ; 'A=100;SF=1.9e-7 1' -> 'β+/β-' vb."""
    if not br_field:
        return "-"
    found = []
    for token in br_field.split(";"):
        token = token.strip()
        if not token or token.startswith("IS="):
            continue
        m = re.match(r"^([A-Za-z0-9+\-]+)[=~<]", token)
        if not m:
            continue
        tag = m.group(1)
        for pat, label in _MODE_LABELS:
            if tag == pat and label not in found:
                found.append(label)
                break
    return "/".join(found) if found else "-"


def _parse_abundance(br_field: str):
    """BR alanından 'IS=91.754 106' -> 91.754 (yüzde) döndürür, yoksa None."""
    if not br_field:
        return None
    m = re.search(r"IS=([\d.]+)", br_field)
    return float(m.group(1)) if m else None


def _half_life_str(val_str: str, unit: str):
    """NUBASE (değer, birim) -> (saniye_float veya None, Türkçe görüntü metni)."""
    val_str = val_str.strip()
    unit = unit.strip()
    if val_str in ("stbl", ""):
        return (float("inf"), "Kararlı")
    if val_str in ("p-unst",):
        return (0.0, "Kararsız (parçacık)")
    try:
        val = float(val_str.replace("#", ""))
    except ValueError:
        return (None, "Bilinmiyor")
    sec = val * _UNIT_SECONDS.get(unit, float("nan"))
    if sec != sec:  # NaN -> tanınmayan birim
        return (None, "Bilinmiyor")

    # Türkçe okunabilir gösterim: küçük ölçekte doğal birim, büyükte yıl (bilimsel gösterim)
    if sec < 1e-6:
        disp = f"{sec*1e9:.3g} ns"
    elif sec < 1e-3:
        disp = f"{sec*1e6:.3g} µs"
    elif sec < 1.0:
        disp = f"{sec*1e3:.3g} ms"
    elif sec < 60.0:
        disp = f"{sec:.3g} s"
    elif sec < 3600.0:
        disp = f"{sec/60.0:.3g} dk"
    elif sec < 86400.0:
        disp = f"{sec/3600.0:.3g} saat"
    elif sec < 3.1556952e7 * 2:
        disp = f"{sec/86400.0:.3g} gün"
    else:
        years = sec / 3.1556952e7
        disp = f"{years:.3g} yıl" if years < 1e4 else f"{years:.4e} yıl".replace("e+0", "e").replace("e+", "e")
    return (sec, disp)


def load_nubase(path: str = _DATA_FILE):
    """Dönüş: (isotopes: dict, abundance: dict[Z] -> list[(A, yüzde)])
    isotopes[key] = {"Z", "A", "mass", "half_life", "half_life_s", "mode", "name"}
    """
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"NUBASE2020 veri dosyası bulunamadı: {path}\n"
            "Bu dosyayı proje klasörüne 'nubase2020.txt' adıyla koyun."
        )

    isotopes = {}
    abundance = {}

    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or len(line) < 90:
                continue
            try:
                A = int(line[0:3])
                Z = int(line[4:7])
                isomer_flag = line[7]
            except ValueError:
                continue
            if isomer_flag != "0":
                continue  # sadece taban durumu

            symbol = line[11:17].strip()
            if not symbol:
                continue
            key = "n" if symbol == "1n" else symbol

            mass_excess_str = line[18:31].strip().replace("#", "")
            if not mass_excess_str:
                continue
            try:
                mass_excess_keV = float(mass_excess_str)
            except ValueError:
                continue
            mass = A + (mass_excess_keV / 1000.0) / U_TO_MEV

            hl_val = line[69:78]
            hl_unit = line[78:80]
            hl_sec, hl_disp = _half_life_str(hl_val, hl_unit)

            br_field = line[119:220] if len(line) > 119 else ""
            mode = _parse_mode(br_field) if hl_disp != "Kararlı" else "-"

            el_name_match = re.match(r"^\d+([A-Za-z]+)$", symbol) or re.match(r"^([A-Za-z]+)$", symbol)
            el_symbol = el_name_match.group(1) if el_name_match else symbol

            isotopes[key] = {
                "Z": Z, "A": A, "mass": mass,
                "half_life": hl_disp, "half_life_s": hl_sec,
                "mode": mode, "name": f"{el_symbol}-{A}",
            }

            ab = _parse_abundance(br_field)
            if ab is not None:
                abundance.setdefault(Z, []).append((A, ab))

    return isotopes, abundance


if __name__ == "__main__":
    iso, ab = load_nubase()
    print(f"{len(iso)} taban-durum nüklit yüklendi.")
    for k in ("n", "1H", "2H", "4He", "56Fe", "40K", "235U"):
        print(k, iso.get(k))
    print("Fe (Z=26) doğal bolluk:", ab.get(26))
