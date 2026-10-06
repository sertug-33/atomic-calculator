"""Fizik motoru referans testleri (AME2020 / NUBASE2020 değerleri).
Çalıştırma:  python -m pytest test_physics.py -q     (veya: python test_physics.py)
"""
from physics_engine import (binding_energy, semf_deviation, separation_energies,
                            calculate_reaction, calculate_moderation, calculate_bateman,
                            decay_channels, specific_activity,
                            decay_series_chain, DECAY_SERIES,
                            gamma_attenuation, gamma_thickness_for_factor,
                            buildup_factor, gamma_attenuation_multiline, RADIOACTIVE_SOURCES,
                            air_mu_en_rho, dose_rate_point_source, dose_rate_multiline,
                            eta_factor, thermal_utilization_factor,
                            resonance_escape_probability, k_infinity)
from density import (atomic_weight, element_atom_density, compound_atom_density,
                     enriched_uranium_weight, enriched_compound_atom_density)
from isotope_explain import explain_isotope, stability_valley_Z
from elements_info import get_element_info, ELEMENT_INFO
import numpy as np

def near(x, ref, tol): assert abs(x - ref) <= tol, f"{x} != {ref} (±{tol})"

def test_baglanma():
    near(binding_energy("56Fe")["Eb_per_A"], 8.79036, 2e-4)   # ders kitabının klasiği
    near(binding_energy("238U")["Eb_per_A"], 7.57013, 2e-4)
    near(binding_energy("4He")["Eb"], 28.29566, 2e-4)
    near(binding_energy("2H")["Eb"], 2.22457, 2e-4)
    near(semf_deviation("56Fe"), -0.34, 0.02)

def test_ayristirma():
    near(separation_energies("56Fe")["sp"], 10.1837, 2e-3)
    near(separation_energies("57Fe")["sn"], 7.6464, 2e-3)

def test_q_degerleri():
    for eq, ref in [("14N + 4He -> 17O + p", -1.19185),     # Rutherford 1919
                    ("3H + p -> 3He + n", -0.76376),
                    ("226Ra -> 222Rn + 4He", 4.87062),
                    ("2H + 3H -> 4He + n", 17.5893),
                    ("6Li + n -> 3H + 4He", 4.78347),
                    ("10B + n -> 7Li + 4He", 2.78991)]:
        near(calculate_reaction(eq)["q_val"], ref, 2e-4)

def test_elastik_ve_esik():
    d = calculate_reaction("14N + 4He -> 14N + 4He")
    assert d["phase"] == "ELASTİK" and not d["is_exothermic"] and d["e_threshold"] is None
    d = calculate_reaction("14N + 4He -> 17O + p")
    near(d["e_threshold"], 1.5326, 2e-3)

def test_moderasyon():
    m = calculate_moderation("H2O")
    near(m["n_coll"], 19.77, 0.02)
    near(m["avg_loss_arith"], 50.0, 1e-9)

def test_bateman_korunum():
    t, A, B, C = calculate_bateman(10.0, 25.0, 1000.0)
    assert np.allclose(A + B + C, 1000.0)
    near(t[np.argmax(B)], np.log(2.5) / (np.log(2)/10 - np.log(2)/25), 0.4)  # B'nin tepe zamanı

def test_atomik_agirlik():
    near(atomic_weight(92), 238.029, 0.01)   # doğal U (periyodik tablo: 238.03)
    near(atomic_weight(8), 15.999, 0.01)     # doğal O
    near(atomic_weight(26), 55.845, 0.01)    # doğal Fe

def test_yogunluk_dogal_uranyum_metal():
    r = element_atom_density(19.1, 92)       # Lamarsh & Baratta Örnek 2.2 tarzı
    near(r["N"], 4.83e22, 0.02e22)

def test_yogunluk_dogal_uo2():
    r = compound_atom_density(10.97, "UO2")
    near(r["elements"]["U"]["N"], 2.446e22, 0.01e22)
    near(r["elements"]["O"]["N"], 4.893e22, 0.01e22)

def test_yogunluk_zenginlestirilmis_uo2():
    mix = enriched_uranium_weight(3.2)
    near(mix["M_mix"], 237.95, 0.05)
    assert 3.2 < mix["atom_pct_235"] < 3.3       # atom% > kütle% (235 daha hafif)
    r = enriched_compound_atom_density(10.5, "UO2", 3.2)
    u = r["elements"]["U"]
    near(u["N_235"] + u["N_238"], u["N"], 1e17)  # tutarlılık: parçalar toplamı

def test_yogunluk_agir_su():
    r = compound_atom_density(1.104, "D2O")
    near(r["M_formula"], 20.028, 0.01)
    near(r["elements"]["D"]["N"], 6.64e22, 0.02e22)

def test_bozunma_kanallari():
    # 226Ra alfa Q-değeri: reaksiyon motoruyla (calculate_reaction) çapraz kontrol
    near(decay_channels("226Ra")["alpha"]["q"], calculate_reaction("226Ra -> 222Rn + 4He")["q_val"], 1e-6)
    # 14C -> 14N + beta- (bilinen değer: 156.5 keV)
    near(decay_channels("14C")["beta_minus"]["q"], 0.1565, 2e-3)
    # 18F: sadece beta+/EC enerjik olarak mümkün (bilinen: 634 keV / 1656 keV), alfa/beta- YASAK olmalı
    d18f = decay_channels("18F")
    near(d18f["beta_plus"]["q"], 0.6339, 2e-3)
    near(d18f["ec"]["q"], 1.6559, 2e-3)
    assert d18f["alpha"]["q"] < 0 and d18f["beta_minus"]["q"] < 0
    # 60Co -> 60Ni + beta- (bilinen: ~2.824 MeV)
    near(decay_channels("60Co")["beta_minus"]["q"], 2.8228, 2e-3)

def test_spesifik_aktivite():
    # 226Ra: Curie birimi tarihsel olarak 1 g 226Ra'nın aktivitesi olarak tanımlanmıştır (~1 Ci/g)
    near(specific_activity("226Ra")["ci_per_g"], 0.989, 0.01)
    # 60Co: bilinen spesifik aktivite ~1131 Ci/g
    near(specific_activity("60Co")["ci_per_g"], 1131.6, 1.0)
    # 137Cs: bilinen spesifik aktivite ~86.9 Ci/g
    near(specific_activity("137Cs")["ci_per_g"], 86.9, 0.5)
    # Kararlı nüklit için tanımsız (None) dönmeli
    assert specific_activity("56Fe") is None

def test_bozunma_serisi_zinciri():
    # Dört klasik doğal seri de literatürdeki bilinen kararlı son ürüne ulaşmalı.
    for start, ref_final in [("238U", "206Pb"), ("235U", "207Pb"),
                              ("232Th", "208Pb"), ("237Np", "205Tl")]:
        chain = decay_series_chain(start)
        assert chain, f"{start} zinciri boş döndü"
        last = chain[-1]
        assert last["stable"], f"{start} zinciri kararlı bir nüklitte bitmedi"
        assert last["key"] == ref_final, f"{start} -> {last['key']} != {ref_final}"
    # U-238 alfa Q-değeri: reaksiyon motoruyla çapraz kontrol (zincirin ilk adımı)
    chain_u238 = decay_series_chain("238U")
    near(chain_u238[0]["q"], calculate_reaction("238U -> 234Th + 4He")["q_val"], 1e-6)

def test_zirhlama():
    # Kurşunda 1 MeV foton için HVL literatürde ~0.85-0.90 cm civarındadır (NIST μ/ρ=0.07102 cm²/g).
    r = gamma_attenuation("Pb", 1.0, 5.0, I0=100.0)
    near(r["hvl"], 0.860, 0.02)
    # Lambert-Beer tutarlılığı: I = I0 * e^(-mu*x) formülü elle de aynı sonucu vermeli
    near(r["I"], 100.0 * np.exp(-r["mu"] * 5.0), 1e-9)
    # Tam olarak 1 HVL kalınlıkta geçirgenlik %50 olmalı
    r_hvl = gamma_attenuation("Pb", 1.0, r["hvl"], I0=100.0)
    near(r_hvl["transmission"], 0.5, 1e-6)
    # gamma_thickness_for_factor tersine tutarlı olmalı: o kalınlıkta gerçekten 10x azalma olmalı
    x10 = gamma_thickness_for_factor("Pb", 1.0, 10.0)
    r10 = gamma_attenuation("Pb", 1.0, x10, I0=100.0)
    near(r10["I"], 10.0, 1e-6)

def test_buildup_faktoru():
    # Sınır koşulu: kalınlık 0 iken buildup 1 olmalı (hiç saçılma yok)
    near(buildup_factor("H2O", 1.0, 0.0), 1.0, 1e-9)
    # Kalınlık arttıkça buildup artmalı (monotonluk) — tüm doğrulanmış malzemelerde
    for mat in ("H2O", "Concrete", "Fe"):
        b_low = buildup_factor(mat, 2.0, 2.0)
        b_high = buildup_factor(mat, 2.0, 8.0)
        assert b_high > b_low > 1.0, f"{mat}: buildup monoton artmıyor"
    # Katsayı tablosu olmayan malzeme (Pb) ve tablo dışı enerji (0.5 MeV) için
    # None dönmeli — uydurma B=1 varsayımı YAPILMAMALI.
    assert buildup_factor("Pb", 1.0, 5.0) is None
    assert buildup_factor("H2O", 0.5, 5.0) is None

def test_coklu_cizgi_kaynak_zayiflamasi():
    # Co-60'ın iki çizgisi (1.173 / 1.333 MeV) ayrı ayrı zayıflatılıp
    # bozunum-başına-foton olasılıklarıyla ağırlıklandırılarak birleştirilmeli.
    m = gamma_attenuation_multiline("Pb", "Co60", 5.0)
    assert len(m["lines"]) == 2
    near(sum(l["intensity_frac"] for l in m["lines"]), 1.0, 1e-6)
    near(m["I0"], 100.0, 1e-6)
    assert 0.0 < m["transmission"] < 1.0
    # Kurşun için buildup verisi olmadığından (has_buildup=False) narrow-beam sonuç dönmeli
    assert m["has_buildup"] is False
    # Cs-137 tek çizgili — sonuç, tek-enerji gamma_attenuation ile birebir aynı olmalı
    single = gamma_attenuation("Concrete", RADIOACTIVE_SOURCES["Cs137"]["lines"][0]["energy"], 10.0, I0=100.0)
    multi = gamma_attenuation_multiline("Concrete", "Cs137", 10.0)
    near(multi["I"], single["I"], 1e-6)
    # Beton + Cs137, 10 mfp'lik kalın bir zırhta buildup verisi (1-10 MeV aralığında
    # DEĞİL, çünkü Cs-137 0.662 MeV) olmadığından has_buildup False olmalı
    assert multi["has_buildup"] is False

def test_doz_hizi():
    # Co-60, 1 Ci (=3.7e10 Bq) @ 1 m, zirhsiz: literatur kurali ~1.32 R/h ~= 11.57 mGy/h
    # (saglik fizigi ders kitaplarinda yaygin kullanilan referans deger).
    r = dose_rate_multiline(3.7e10, "Co60", 100.0)
    near(r["dose_rate_mGy_h"], 11.57, 1.0)
    # Mesafe 2 katina cikarsa doz hizi ters-kare yasasiyla 1/4'e dusmeli.
    r2 = dose_rate_multiline(3.7e10, "Co60", 200.0)
    near(r2["dose_rate_mGy_h"], r["dose_rate_mGy_h"] / 4.0, 0.05)
    # Tek enerjili nokta-kaynak fonksiyonu, coklu-hat sarmalayicisinin tek satirlik
    # toplamiyla tutarli olmali (Cs-137 tek cizgili).
    single = dose_rate_point_source(3.7e10, RADIOACTIVE_SOURCES["Cs137"]["lines"][0]["energy"],
                                     RADIOACTIVE_SOURCES["Cs137"]["lines"][0]["intensity"], 100.0)
    multi = dose_rate_multiline(3.7e10, "Cs137", 100.0)
    near(multi["dose_rate_Gy_s"], single["dose_rate_Gy_s"], 1e-12)
    # Zirhlama eklenince doz hizi azalmali (Pb 5 cm, Co-60).
    r_shielded = dose_rate_multiline(3.7e10, "Co60", 100.0, shield_material="Pb", shield_thickness_cm=5.0)
    assert r_shielded["dose_rate_mGy_h"] < r["dose_rate_mGy_h"]
    # air_mu_en_rho tablo ucu disina cikmamali (sinir degerlere kilitlenmeli).
    near(air_mu_en_rho(0.01), air_mu_en_rho(0.02), 1e-9)
    near(air_mu_en_rho(50.0), air_mu_en_rho(10.0), 1e-9)


def test_kritiklik_dort_faktor():
    # Saf U-235 (e=%100): literatur (nuclear-power.com) asimptotik eta ~= 2.08
    near(eta_factor(100.0)["eta"], 2.08, 0.02)
    # Dogal uranyum (e=%0.711): literatur eta ~= 1.34
    near(eta_factor(0.711)["eta"], 1.34, 0.02)
    # Gecersiz zenginlestirme (0 veya >100) None donmeli
    assert eta_factor(0.0) is None
    assert eta_factor(150.0) is None

    # f: N_F/N_M arttikca (yakit yogunlasinca) f 1'e yaklasmali (monotonluk)
    f_low = thermal_utilization_factor("H2O", 3.2, 0.001)["f"]
    f_high = thermal_utilization_factor("H2O", 3.2, 1.0)["f"]
    assert 0.0 < f_low < f_high < 1.0

    # p: I_eff=0 iken rezonansta hic kayip yok -> p=1 sinir kosulu
    p0 = resonance_escape_probability("H2O", 3.2, 0.01, 0.0)["p"]
    near(p0, 1.0, 1e-9)
    # I_eff arttikca p azalmali (monoton)
    p_low_I = resonance_escape_probability("H2O", 3.2, 0.01, 5.0)["p"]
    p_high_I = resonance_escape_probability("H2O", 3.2, 0.01, 50.0)["p"]
    assert 1.0 > p_low_I > p_high_I > 0.0

    # k_inf = eta*f*p*epsilon tutarliligi (elle carpimla ayni sonucu vermeli)
    k = k_infinity("H2O", 3.2, 0.05, 10.0, epsilon=1.0)
    near(k["k_inf"], k["eta"] * k["f"] * k["p"] * k["epsilon"], 1e-9)
    # Gecersiz girdi (eksik/negatif) icin None donmeli
    assert k_infinity("H2O", 3.2, -1.0, 10.0) is None
    assert k_infinity("H2O", 3.2, 0.05, -5.0) is None


def test_izotop_aciklama():
    # Element bilgi tablosu Z=1..118 tam olmali, bosluk/duplike olmamali.
    assert len(ELEMENT_INFO) == 118
    assert all(z in ELEMENT_INFO for z in range(1, 119))
    fe = get_element_info(26)
    assert fe["symbol"] == "Fe" and fe["tr_name"] == "Demir"
    assert get_element_info(0) is None and get_element_info(200) is None

    # Kararlilik-vadisi Z tahmini: Fe-56 icin gercek Z=26'ya yakin olmali (SEMF yaklasimi).
    near(stability_valley_Z(56), 26.0, 1.5)

    # Kararli bir cift-sihirli cekirdek (208Pb): stable=True, mode aciklamasinda
    # "KARARLI" gecmeli, fiziksel metinde sihirli sayi notu bulunmali.
    pb = explain_isotope("208Pb")
    assert pb["is_stable"] is True
    assert "KARARLI" in pb["physical_text"] or "kararlıdır" in pb["physical_text"].lower()
    assert "SİHİRLİ" in pb["physical_text"].upper()
    assert pb["Z"] == 82 and pb["N"] == 126

    # Notron-zengin bir beta- yayicisi (14C): dogru zenginlik yonu ve pozitif Q ile
    # tutarli aciklama uretmeli (Q>0 gercek AME2020 kutle farkindan gelir).
    c14 = explain_isotope("14C")
    assert c14["is_stable"] is False
    assert c14["richness"] == "notron_zengin"
    assert "β-" in c14["mode"] or "β⁻" in c14["physical_text"]
    assert "+0." in c14["physical_text"] or "Q=" in c14["physical_text"]

    # Proton-zengin bir β+/EC yayicisi yonu dogru tespit edilmeli.
    na22 = explain_isotope("22Na")
    assert na22["richness"] in ("proton_zengin", "dengeli")

    # Bilinmeyen/gecersiz anahtar icin None donmeli (uydurma veri uretilmemeli).
    assert explain_isotope("999Xx") is None
    assert explain_isotope("n") is None and explain_isotope("p") is None

    # Turkce kucuk harfe cevirme hatasi (dotless-I) duzeltmesi: "bozunmasi" DEGIL
    # dogru Turkce "bozunması" gecmeli.
    u238 = explain_isotope("238U")
    assert "bozunması" in u238["physical_text"]
    assert "bozunmasi" not in u238["physical_text"]


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"): fn(); print("✓", name)
