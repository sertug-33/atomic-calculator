# -*- coding: utf-8 -*-
"""
Nüklit Haritası künye paneli için açıklama motoru.

Bu modül, tıklanan bir nüklit için İKİ tür açıklama üretir:
  1) FİZİKSEL açıklama: kararlılık vadisine göre konumu (SEMF'ten türetilen
     aynı katsayılarla — physics_engine.semf_breakdown ile TUTARLI), sihirli
     sayılara yakınlığı, çift-çift/tek-tek yapısı ve — kararsızsa —
     decay_channels() ile AME2020 kütlelerinden HESAPLANAN gerçek Q-değerine
     dayanan bozunma gerekçesi.
  2) GENEL/SÖZEL açıklama: elements_info.py'deki ansiklopedik element bilgisi
     + bu izotopa özgü bağlam (doğal bolluk, yarı ömür, kütle numarası notu).

Disiplin: hiçbir sayısal değer uydurulmaz. Fiziksel açıklamadaki tüm sayılar
ya doğrudan veri tabanından (NUBASE2020/AME2020) ya da mevcut, zaten
doğrulanmış physics_engine fonksiyonlarından (decay_channels, semf_breakdown)
gelir. SEMF kararlılık-vadisi tahmini yalnızca NİTELİKSEL bir konumlandırma
(proton-zengin / nötron-zengin) için kullanılır, kesin bir kehanet olarak
sunulmaz.
"""
from database import ISOTOPES, ABUNDANCE
from physics_engine import decay_channels, atomic_mass
from elements_info import get_element_info

# Sihirli sayılar (kabuk modeli, kapalı proton/nötron kabukları)
MAGIC_NUMBERS = (2, 8, 20, 28, 50, 82, 126)

# semf_breakdown ile AYNI katsayılar (physics_engine.py) — tutarlılık için
# buraya da kopyalanır (dairesel import'tan kaçınmak için sabit olarak).
_AC = 0.714   # Coulomb katsayısı (MeV)
_AA = 23.2    # Asimetri katsayısı (MeV)


def stability_valley_Z(A: int) -> float:
    """SEMF'in Z'ye göre extremum koşulundan (dB/dZ=0) türetilen, verilen A
    için 'en kararlı' proton sayısının yaklaşık tahmini:
        Z_A ≈ A / (2 + (a_c / (2 a_a)) · A^(2/3))
    Bu, ders kitaplarında (örn. Krane, Lamarsh) verilen Z ≈ A/(1.98+0.0155·A^(2/3))
    yaklaşımıyla aynı formun bir türevidir; burada projede zaten kullanılan
    a_c=0.714, a_a=23.2 katsayılarıyla iç tutarlılık sağlanır."""
    if A <= 0:
        return 0.0
    return A / (2.0 + (_AC / (2.0 * _AA)) * A ** (2.0 / 3.0))


def _parity_label(Z: int, N: int) -> str:
    z_even, n_even = (Z % 2 == 0), (N % 2 == 0)
    if z_even and n_even:
        return "çift-çift (Z çift, N çift) — eşleşme enerjisi ekstra kararlılık sağlar"
    if not z_even and not n_even:
        return "tek-tek (Z tek, N tek) — eşleşme enerjisi kaybı nedeniyle en az kararlı grup"
    return "tek-çift (Z ve N'den biri tek)"


def _magic_label(Z: int, N: int) -> str | None:
    z_magic = Z in MAGIC_NUMBERS
    n_magic = N in MAGIC_NUMBERS
    if z_magic and n_magic:
        return f"ÇİFT SİHİRLİ çekirdek (Z={Z} ve N={N} ikisi de kapalı kabuk) — olağanüstü kararlı"
    if z_magic:
        return f"Z={Z} sihirli sayı (kapalı proton kabuğu) — ekstra kararlılık katkısı"
    if n_magic:
        return f"N={N} sihirli sayı (kapalı nötron kabuğu) — ekstra kararlılık katkısı"
    return None


_MODE_REASON = {
    "β-": ("nötron fazlası", "beta_minus"),
    "2β-": ("nötron fazlası (çift β⁻)", "beta_minus"),
    "β+": ("proton fazlası", "beta_plus"),
    "EC": ("proton fazlası", "ec"),
    "2β+": ("proton fazlası (çift β⁺)", "beta_plus"),
    "α": ("ağır/proton bakımından zengin çekirdeklerde Coulomb itmesinin baskınlığı", "alpha"),
}


def _format_q(q: float) -> str:
    return f"{q:+.3f} MeV"


# decay_channels()'ın ürettiği etiketlerin doğru Türkçe küçük harfli biçimleri
# — Python'un varsayılan (Türkçe olmayan) str.lower()'ı "I" harfini yanlışlıkla
# noktalı "i"ye çevirip "bozunması" yerine "bozunmasi" yazdığı için elle eşlenir.
_LABEL_LOWER = {
    "α BOZUNMASI": "α bozunması",
    "β⁻ BOZUNMASI": "β⁻ bozunması",
    "β⁺ BOZUNMASI": "β⁺ bozunması",
    "ELEKTRON YAKALAMA (EC)": "elektron yakalama (EC)",
}


def _label_lower(label: str) -> str:
    return _LABEL_LOWER.get(label, label.lower())


def explain_isotope(iso_key: str) -> dict | None:
    """Verilen nüklit anahtarı için tam açıklama paketini döndürür, yoksa None."""
    if iso_key not in ISOTOPES or iso_key in ("n", "p"):
        return None
    iso = ISOTOPES[iso_key]
    Z, A = iso["Z"], iso["A"]
    N = A - Z
    is_stable = iso.get("half_life") == "Kararlı"
    mode = iso.get("mode", "-")
    primary_mode = mode.split("/")[0] if mode and mode != "-" else None

    el = get_element_info(Z)
    el_symbol = el["symbol"] if el else iso["name"].split("-")[0]
    el_tr_name = el["tr_name"] if el else el_symbol
    display_name = f"{el_tr_name}-{A}  ({el_symbol}-{A})"

    Z_est = stability_valley_Z(A)
    delta_z = Z - Z_est
    # A arttıkça vadi genişliği de artar; kaba bir tolerans olarak ~%6·Z_est kullanılır
    # (yalnızca "belirgin fazlalık" var mı yok mu demek için, kesin sınır değildir).
    tol = max(0.6, 0.06 * Z_est)
    if delta_z > tol:
        richness = "proton_zengin"
        richness_text = (f"Z={Z}, bu kütle numarası (A={A}) için SEMF'ten tahmin edilen "
                          f"kararlılık-vadisi merkezinden (Z≈{Z_est:.1f}) belirgin biçimde "
                          f"YÜKSEK — yani çekirdek görece PROTON BAKIMINDAN ZENGİN.")
    elif delta_z < -tol:
        richness = "notron_zengin"
        richness_text = (f"Z={Z}, bu kütle numarası (A={A}) için SEMF'ten tahmin edilen "
                          f"kararlılık-vadisi merkezinden (Z≈{Z_est:.1f}) belirgin biçimde "
                          f"DÜŞÜK — yani çekirdek görece NÖTRON BAKIMINDAN ZENGİN.")
    else:
        richness = "dengeli"
        richness_text = (f"Z={Z}, bu kütle numarası (A={A}) için SEMF'ten tahmin edilen "
                          f"kararlılık-vadisi merkezine (Z≈{Z_est:.1f}) yakın — proton/nötron "
                          f"dengesi kararlılık vadisinin ortasına oturuyor.")

    parity_text = _parity_label(Z, N)
    magic_text = _magic_label(Z, N)

    channels = decay_channels(iso_key) or {}

    # --- Fiziksel açıklama metni ------------------------------------------------
    phys_lines = [richness_text, f"Çekirdek yapısı: {parity_text}."]
    if magic_text:
        phys_lines.append(magic_text + ".")

    if is_stable:
        phys_lines.append(
            "NUBASE2020 verisine göre bu nüklit KARARLIDIR (gözlenen bir bozunum modu yok). "
            "Kararlılığı; yukarıdaki N/Z dengesi, çift-çift/tek-tek yapı ve (varsa) kapalı "
            "kabuk etkilerinin bir arada, çekirdeği enerjik olarak mümkün hiçbir bozunum "
            "kanalına göre daha düşük (ya da onlarla kıyaslanabilir/imkânsız) bir toplam "
            "enerjide tutmasıyla açıklanır."
        )
        # Enerjik olarak "mümkün" ama gözlenmeyen kanal varsa (Q>0 fakat kararlı
        # işaretli — çok uzun ömürlü/yasaklı geçiş olabilir) bunu da şeffafça belirt.
        possible = [(k, v) for k, v in channels.items() if v["q"] > 0]
        if possible:
            best = max(possible, key=lambda kv: kv[1]["q"])
            phys_lines.append(
                f"Not: AME2020 kütlelerine göre {_label_lower(best[1]['label'])} kanalı enerjik "
                f"olarak mümkün görünüyor (Q={_format_q(best[1]['q'])}) ama NUBASE2020 bunu "
                f"gözlenmemiş/ihmal edilebilir olarak işaretliyor — bu genellikle spin/parite "
                f"seçim kurallarının geçişi son derece yavaşlattığı (pratikte 'kararlı' "
                f"sayılabilecek kadar uzun ömürlü) durumlarda görülür."
            )
    else:
        reason, chan_key = _MODE_REASON.get(primary_mode, (None, None))
        q_info = channels.get(chan_key) if chan_key else None
        if reason and q_info:
            phys_lines.append(
                f"Baskın bozunum modu: {primary_mode} — {reason} nedeniyle çekirdek "
                f"kararlılık vadisine yaklaşmak üzere {_label_lower(q_info['label'])} yapıyor "
                f"(→ {q_info['daughter']}). AME2020 kütle farkından hesaplanan "
                f"Q-değeri: {_format_q(q_info['q'])} (Q>0, yani enerjik olarak mümkün ve "
                f"kendiliğinden gerçekleşebilir)."
            )
        elif primary_mode:
            phys_lines.append(
                f"Baskın bozunum modu: {primary_mode}. Bu kanal için Q-değeri, mevcut "
                f"kütle tablosunda ürün nüklit eksikliği nedeniyle bu sürümde hesaplanamadı."
            )
        else:
            phys_lines.append("Bozunum modu NUBASE2020 kaydında belirtilmemiş.")

        hl = iso.get("half_life", "?")
        phys_lines.append(f"Yarı ömür: {hl}.")

    physical_text = " ".join(phys_lines)

    # --- Genel/sözel açıklama metni ---------------------------------------------
    general_lines = []
    if el:
        general_lines.append(f"{el_tr_name} ({el_symbol}), periyodik tabloda Z={Z} sıra "
                              f"numarasıyla yer alan bir {el['category'].lower()}dir. {el['desc']}")
    else:
        general_lines.append(f"{el_symbol} elementi hakkında genel bilgi bu sürümde mevcut değil.")

    ab_list = ABUNDANCE.get(Z, [])
    ab_match = next((pct for a, pct in ab_list if a == A), None)
    if is_stable and ab_match is not None:
        general_lines.append(
            f"Bu izotop ({el_symbol}-{A}), doğal {el_tr_name.lower()} içinde yaklaşık "
            f"%{ab_match:.2f} oranında bulunur."
        )
    elif is_stable and ab_list:
        general_lines.append(
            f"Bu izotop kararlı olsa da doğal {el_tr_name.lower()} karışımında ölçülebilir "
            f"bir bolluk kaydı bulunmuyor (çok düşük/ihmal edilebilir olabilir)."
        )
    elif not is_stable:
        general_lines.append(
            f"Bu izotop doğada ya bulunmaz ya da yalnızca başka bir nüklidin bozunumu "
            f"sonucu geçici olarak ortaya çıkar; kararlı değildir."
        )

    general_text = " ".join(general_lines)

    return {
        "key": iso_key,
        "display_name": display_name,
        "symbol": el_symbol,
        "tr_name": el_tr_name,
        "Z": Z, "A": A, "N": N,
        "is_stable": is_stable,
        "mode": mode,
        "half_life": iso.get("half_life", "?"),
        "Z_stable_estimate": Z_est,
        "richness": richness,
        "magic": magic_text is not None,
        "physical_text": physical_text,
        "general_text": general_text,
        "channels": channels,
    }
