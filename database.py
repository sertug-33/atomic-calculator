# -*- coding: utf-8 -*-
"""
Nükleer veri katmanı.

Bu modül artık AME2020/NUBASE2020'nin TAMAMINI (3558 taban-durum nüklit)
`nubase.py` üzerinden yükler — elle girilmiş 158 satırlık sabit liste yerine.
Aynı klasörde `nubase2020.txt` (NUBASE2020 ham veri dosyası) bulunmalıdır.

Geriye dönük uyumluluk: ISOTOPES ve ALIASES'ın şeması (anahtar/alan adları)
önceki elle-yazılmış sürümle birebir aynıdır, bu yüzden physics_engine.py,
main.py ve ui_components.py hiçbir değişiklik yapılmadan çalışmaya devam eder.
"""
from nubase import load_nubase

# CODATA 2018 sabitleri
U_TO_MEV = 931.49410242        # 1 u = 931.494 MeV/c^2
M_P  = 1.007276466621          # ÇIPLAK proton kütlesi (u) — atomik kütlelerle KARIŞTIRMAYIN
M_N  = 1.00866491595           # nötron kütlesi (u)
M_E  = 0.000548579909065       # elektron kütlesi (u)
M_H1 = 1.00782503223           # ATOMİK 1H kütlesi (p + e), AME2020

# NOT: ISOTOPES içindeki kütleler ATOMİK kütlelerdir (elektronlar dahil, AME2020).
# Bağlanma enerjisi, Q-değeri ve Sp hesaplarında protonlar için M_P değil M_H1
# kullanılmalıdır; böylece iki taraftaki elektron sayıları kendiliğinden dengelenir.
# ("p" kaydı çıplak proton olarak kalır; physics_engine reaksiyonlarda onu 1H'ye çevirir.)

ALIASES = {
    "n": "n", "p": "p", "d": "2H", "t": "3H",
    "alpha": "4He", "alfa": "4He", "α": "4He", "He-4": "4He"
}

# --- NUBASE2020'yi yükle -----------------------------------------------------
ISOTOPES, ABUNDANCE = load_nubase()

# "p" (çıplak proton) NUBASE'de yok (orada sadece atomik 1H var) — elle ekleniyor.
# Reaksiyon/bozunma hesaplarında ("p" seçildiğinde) physics_engine bunun yerine
# atomik 1H kütlesini (M_H1) kullanır; bu kayıt yalnızca arayüzde "Proton" olarak
# ayrı bir mermi/emisyon seçeneği görünmesi için tutuluyor.
ISOTOPES["p"] = {
    "Z": 1, "A": 1, "mass": M_P, "half_life": "Kararlı",
    "half_life_s": float("inf"), "mode": "-", "name": "Proton",
}

# ---------------------------------------------------------------------------
# Doğal bulunurluk / atom yoğunluğu modülleri için: element sembolü -> Z eşlemesi.
ELEMENT_Z = {}
for _key, _v in ISOTOPES.items():
    if _key in ("n", "p"):
        continue
    _sym = _v["name"].split("-")[0]
    ELEMENT_Z.setdefault(_sym, _v["Z"])

if __name__ == "__main__":
    print(f"{len(ISOTOPES)} nüklit yüklü (NUBASE2020).")
    print("56Fe:", ISOTOPES["56Fe"])
    print("40K :", ISOTOPES["40K"])
    print("Fe doğal bolluk (ABUNDANCE[26]):", ABUNDANCE[26])
