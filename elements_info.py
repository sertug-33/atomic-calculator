# -*- coding: utf-8 -*-
"""
Elementlerin genel (kimyasal/ansiklopedik) bilgi tablosu — Z=1..118.

Bu veri, izotop bazlı nükleer verilerden (NUBASE2020/AME2020) ayrı, genel
kimya/element bilgisidir: sembol, Türkçe adı, periyodik tablo kategorisi ve
kısa bir tanıtım cümlesi (öne çıkan kullanım alanı veya tarihsel not).

Kapsam kasıtlı olarak GENEL ve iyi bilinen ansiklopedik gerçeklerle
sınırlıdır (element kategorisi, başlıca kullanım alanı, keşif/adlandırma
notu) — spesifik sayısal nükleer verilerle (yarı ömür, Q-değeri, bolluk vb.)
KARIŞTIRILMAMALIDIR; onlar için nubase2020.txt / database.py / physics_engine.py
kullanılır. Nüklit Haritası künye panelinde, seçilen izotopun bağlı olduğu
elementin genel tanıtımı olarak gösterilir.
"""

# category kodları: alkali, alkali_toprak, gecis, post_gecis, metaloid,
# ametal, halojen, soy_gaz, lantanit, aktinit, super_agir

ELEMENT_INFO = {
    1:   ("H",  "Hidrojen",      "ametal",       "Evrenin en bol ve en hafif elementi; su ve tüm organik moleküllerin temel bileşeni, yıldızların füzyon yakıtıdır."),
    2:   ("He", "Helyum",        "soy_gaz",      "Havadan hafif soy gaz; balonlarda ve süper iletken mıknatısların kriyojenik soğutmasında kullanılır."),
    3:   ("Li", "Lityum",        "alkali",       "En hafif metal; lityum-iyon pillerin ve bazı ruh hâli dengeleyici ilaçların temel bileşenidir."),
    4:   ("Be", "Berilyum",      "alkali_toprak","Hafif ve sert; Chadwick'in 1932'de nötronu keşfettiği Be(α,n) reaksiyonunda kullanıldı, nötron kaynaklarında hâlâ kullanılır."),
    5:   ("B",  "Bor",           "metaloid",     "Borosilikat cam ve seramiklerde, yüksek nötron yakalama kesiti nedeniyle BNCT kanser tedavisinde ve reaktör kontrol çubuklarında kullanılır."),
    6:   ("C",  "Karbon",        "ametal",       "Tüm organik yaşamın omurgası; elmas, grafit ve grafen gibi çok farklı allotropları vardır, reaktörlerde moderatör olarak da kullanılır."),
    7:   ("N",  "Azot",          "ametal",       "Atmosferin yaklaşık %78'ini oluşturur; gübre ve patlayıcı üretiminde temel hammaddedir."),
    8:   ("O",  "Oksijen",       "ametal",       "Solunum ve yanmanın temel gazı; yer kabuğunda kütlece en bol bulunan elementtir."),
    9:   ("F",  "Flor",          "halojen",      "En elektronegatif element; diş macunlarında (florür) ve teflon (PTFE) üretiminde kullanılır."),
    10:  ("Ne", "Neon",          "soy_gaz",      "Kimyasal olarak neredeyse tamamen inert; karakteristik kırmızı-turuncu ışığıyla neon tabelalarda tanınır."),
    11:  ("Na", "Sodyum",        "alkali",       "Sofra tuzunun (NaCl) bileşeni; vücutta sinir iletiminde ve bazı hızlı nötron reaktörlerinde soğutucu olarak kullanılır."),
    12:  ("Mg", "Magnezyum",     "alkali_toprak","Hafif yapı alaşımlarında kullanılır; klorofil molekülünün merkez atomudur."),
    13:  ("Al", "Alüminyum",     "post_gecis",   "Hafifliği ve korozyon direnci nedeniyle uçak, ambalaj ve inşaat sektöründe yaygın kullanılan bir metaldir."),
    14:  ("Si", "Silisyum",      "metaloid",     "Modern elektroniğin (bilgisayar çipleri, güneş panelleri) temel yarı iletken malzemesidir."),
    15:  ("P",  "Fosfor",        "ametal",       "DNA, RNA ve ATP'nin yapı taşıdır; gübre üretiminde temel hammaddedir."),
    16:  ("S",  "Kükürt",        "ametal",       "Sülfürik asit üretiminde ve bazı amino asitlerin (sistein, metiyonin) yapısında bulunur."),
    17:  ("Cl", "Klor",          "halojen",      "İçme suyu dezenfeksiyonunda ve PVC plastiği üretiminde yaygın olarak kullanılır."),
    18:  ("Ar", "Argon",         "soy_gaz",      "Havada bolca bulunan inert bir gaz; kaynak işlemlerinde koruyucu atmosfer olarak kullanılır."),
    19:  ("K",  "Potasyum",      "alkali",       "Hücre içi sinir ve kas fonksiyonlarında kritik rol oynar; doğal K-40 izotopu insan vücudundaki başlıca doğal radyoaktivite kaynağıdır."),
    20:  ("Ca", "Kalsiyum",      "alkali_toprak","Kemik ve diş yapısının temel mineral bileşenidir."),
    21:  ("Sc", "Skandiyum",     "gecis",        "Hafif alüminyum-skandiyum alaşımlarında (havacılık, spor ekipmanları) kullanılır."),
    22:  ("Ti", "Titanyum",      "gecis",        "Yüksek mukavemet/düşük ağırlık oranı ve biyouyumluluğu nedeniyle havacılık ve ortopedik implantlarda kullanılır."),
    23:  ("V",  "Vanadyum",      "gecis",        "Çelik alaşımlarına sertlik ve yorulma direnci katmak için kullanılır."),
    24:  ("Cr", "Krom",          "gecis",        "Paslanmaz çeliğe karakteristik korozyon direncini kazandıran alaşım elementidir."),
    25:  ("Mn", "Mangan",        "gecis",        "Çelik üretiminde sertleştirici ve oksijen giderici olarak yaygın kullanılan bir alaşım elementidir."),
    26:  ("Fe", "Demir",         "gecis",        "Çelik üretiminin temelidir; hemoglobinde oksijen taşınmasında ve nükleer bağlanma enerjisi eğrisinin zirvesinde (en kararlı bölge) yer alır."),
    27:  ("Co", "Kobalt",        "gecis",        "Lityum-iyon pil katotlarında kullanılır; Co-60 izotopu radyoterapi ve endüstriyel radyografide yaygın bir gama kaynağıdır."),
    28:  ("Ni", "Nikel",         "gecis",        "Paslanmaz çelik alaşımlarında ve madeni paralarda yaygın olarak kullanılır."),
    29:  ("Cu", "Bakır",         "gecis",        "Gümüşten sonra en yüksek elektrik iletkenliğine sahip metal; kablo ve elektronik devrelerde temel malzemedir."),
    30:  ("Zn", "Çinko",         "gecis",        "Çeliği paslanmaya karşı kaplamada (galvanizleme) ve pirinç alaşımında kullanılır."),
    31:  ("Ga", "Galyum",        "post_gecis",   "Oda sıcaklığına yakın erime noktasıyla bilinir; GaAs ve GaN yarı iletkenleri LED ve yüksek frekans elektroniğinde kullanılır."),
    32:  ("Ge", "Germanyum",     "metaloid",     "İlk transistörlerde kullanıldı; günümüzde fiber optik kablolarda ve kızılötesi optikte kullanılır."),
    33:  ("As", "Arsenik",       "metaloid",     "GaAs yarı iletken bileşiklerinde kullanılır; tarihsel olarak zehirliliğiyle de bilinir."),
    34:  ("Se", "Selenyum",      "ametal",       "Fotokopi makinelerinde ve cam üretiminde kullanılır; vücutta iz miktarda gerekli bir antioksidan elementtir."),
    35:  ("Br", "Brom",          "halojen",      "Oda sıcaklığında sıvı hâlde bulunan tek ametal element; alev geciktirici kimyasallarda kullanılır."),
    36:  ("Kr", "Kripton",       "soy_gaz",      "Yüksek performanslı aydınlatma ve eski nesil fotoğraf flaşlarında kullanılan inert bir gazdır."),
    37:  ("Rb", "Rubidyum",      "alkali",       "Atomik saatlerde frekans referansı olarak kullanılır; oda sıcaklığına yakın erime noktasına sahiptir."),
    38:  ("Sr", "Stronsiyum",    "alkali_toprak","Havai fişeklerde kırmızı renk verir; Sr-90 izotopu ise fisyon ürünü olarak nükleer atıkta bulunan önemli bir radyoizotoptur."),
    39:  ("Y",  "İtriyum",       "gecis",        "LED fosforlarında, lazer kristallerinde (YAG) ve süper iletken seramiklerde kullanılır."),
    40:  ("Zr", "Zirkonyum",     "gecis",        "Düşük nötron yutma kesiti ve korozyon direnci nedeniyle nükleer reaktör yakıt çubuğu kaplamalarında (Zircaloy) kullanılır."),
    41:  ("Nb", "Niyobyum",      "gecis",        "Süper iletken mıknatıslarda (MRI, parçacık hızlandırıcıları) ve çelik alaşımlarında kullanılır."),
    42:  ("Mo", "Molibden",      "gecis",        "Yüksek sıcaklığa dayanıklı çelik alaşımlarında ve Mo-99 üzerinden tıbbi Tc-99m üretiminde kullanılır."),
    43:  ("Tc", "Teknesyum",     "gecis",        "Doğada kararlı izotopu bulunmayan ilk yapay element; Tc-99m tıbbi sintigrafi görüntülemesinde dünya çapında en yaygın kullanılan radyoizotoptur."),
    44:  ("Ru", "Rutenyum",      "gecis",        "Platin grubu metal; elektrik kontaklarında ve bazı kimyasal katalizörlerde kullanılır."),
    45:  ("Rh", "Rodyum",        "gecis",        "Platin grubu metal; otomobil katalitik konvertörlerinde azot oksitleri indirgemek için kullanılır."),
    46:  ("Pd", "Palladyum",     "gecis",        "Katalitik konvertörlerde ve hidrojen gazını depolama/arıtma özelliğiyle bilinir."),
    47:  ("Ag", "Gümüş",         "gecis",        "Bilinen en yüksek elektrik ve ısı iletkenliğine sahip metal; takı, ayna kaplama ve elektronik kontaklarda kullanılır."),
    48:  ("Cd", "Kadmiyum",      "gecis",        "Nikel-kadmiyum pillerde ve yüksek nötron yutma kesiti nedeniyle nükleer reaktör kontrol çubuklarında kullanılır."),
    49:  ("In", "İndiyum",      "post_gecis",   "Dokunmatik ekran ve LCD panellerde şeffaf iletken kaplama (indiyum kalay oksit, ITO) olarak kullanılır."),
    50:  ("Sn", "Kalay",         "post_gecis",   "Lehim alaşımlarında ve bronzun (bakır-kalay alaşımı) bileşeninde kullanılır; kararlı izotop sayısı en fazla olan elementtir (Z=50, sihirli sayı)."),
    51:  ("Sb", "Antimon",       "metaloid",     "Alev geciktirici bileşiklerde ve kurşun akü alaşımlarında kullanılır."),
    52:  ("Te", "Tellür",        "metaloid",     "CdTe ince film güneş panellerinde ve bazı metal alaşımlarının işlenebilirliğini artırmada kullanılır."),
    53:  ("I",  "İyot",          "halojen",      "Tiroid hormonu üretimi için gereklidir; tıbbi dezenfektan ve I-131 ile tiroid tedavisinde kullanılır."),
    54:  ("Xe", "Ksenon",        "soy_gaz",      "Yüksek yoğunluklu araç farlarında ve uydu/uzay araçlarının iyon itki motorlarında kullanılan inert bir gazdır."),
    55:  ("Cs", "Sezyum",        "alkali",       "Cs-133 geçişi saniyenin uluslararası tanımının temelidir (atomik saatler); Cs-137 ise önemli bir fisyon ürünü radyoizotoptur."),
    56:  ("Ba", "Baryum",        "alkali_toprak","Radyolojide sindirim sistemi görüntülemesinde kontrast madde ('baryum lahmacunu') olarak kullanılır."),
    57:  ("La", "Lantan",        "lantanit",     "Nikel-metal hidrit pillerde ve yüksek kaliteli kamera/teleskop lenslerinde kullanılır."),
    58:  ("Ce", "Seryum",        "lantanit",     "Çakmak taşlarında ve cam parlatma tozlarında kullanılan en bol lantanit elementtir."),
    59:  ("Pr", "Praseodim",     "lantanit",     "Güçlü kalıcı mıknatıslarda ve camlara/seramiklere sarı-yeşil renk vermede kullanılır."),
    60:  ("Nd", "Neodim",        "lantanit",     "Bilinen en güçlü kalıcı mıknatısların (Nd-Fe-B) ve lazer kristallerinin temel bileşenidir."),
    61:  ("Pm", "Prometyum",     "lantanit",     "Doğada kararlı izotopu bulunmayan nadir lantanit; radyoizotop pillerde ve ışıldayan boyalarda kullanılmıştır."),
    62:  ("Sm", "Samaryum",      "lantanit",     "Sm-Co kalıcı mıknatıslarında ve yüksek nötron yutma kesitiyle reaktör kontrol çubuklarında kullanılır."),
    63:  ("Eu", "Evropyum",      "lantanit",     "TV ve LED ekranlarda kırmızı ve mavi fosfor bileşeni olarak kullanılır."),
    64:  ("Gd", "Gadolinyum",    "lantanit",     "MRI kontrast maddelerinde ve bilinen en yüksek termal nötron yutma kesitine sahip elementlerden biri olarak reaktör kontrolünde kullanılır."),
    65:  ("Tb", "Terbiyum",      "lantanit",     "Yeşil fosforlarda ve manyeto-optik veri depolama malzemelerinde kullanılır."),
    66:  ("Dy", "Disprozyum",    "lantanit",     "Yüksek sıcaklıkta manyetik özelliğini koruyan kalıcı mıknatıs alaşımlarında kullanılır."),
    67:  ("Ho", "Holmiyum",      "lantanit",     "Elementler arasında en yüksek manyetik momente sahiptir; lazer ve özel mıknatıslarda kullanılır."),
    68:  ("Er", "Erbiyum",       "lantanit",     "Fiber optik iletişim ağlarındaki erbiyum katkılı fiber amplifikatörlerde (EDFA) kullanılır."),
    69:  ("Tm", "Tulyum",        "lantanit",     "Taşınabilir X-ışını cihazlarında kompakt bir radyasyon kaynağı olarak kullanılır."),
    70:  ("Yb", "İterbiyum",     "lantanit",     "Atomik saatlerde ve fiber lazer sistemlerinde kullanılır."),
    71:  ("Lu", "Lutesyum",      "lantanit",     "PET tarayıcılardaki dedektör kristallerinin (LSO) temel bileşenidir."),
    72:  ("Hf", "Hafniyum",      "gecis",        "Yüksek nötron yutma kesiti nedeniyle nükleer reaktör kontrol çubuklarında ve bilgisayar çiplerinde yalıtkan katman olarak kullanılır."),
    73:  ("Ta", "Tantal",        "gecis",        "Elektronik kondansatörlerde ve yüksek biyouyumluluğu nedeniyle vücut implantlarında kullanılır."),
    74:  ("W",  "Volfram (Tungsten)", "gecis",   "Bilinen en yüksek erime noktasına sahip metal; ampul filamanlarında ve ağır/sert alaşımlarda kullanılır."),
    75:  ("Re", "Renyum",        "gecis",        "Jet motoru süper alaşımlarında yüksek sıcaklık dayanımı sağlamak için kullanılır."),
    76:  ("Os", "Osmiyum",       "gecis",        "Bilinen en yoğun doğal element; dolma kalem uçları ve sert alaşımlarda kullanılır."),
    77:  ("Ir", "İridyum",       "gecis",        "En korozyona dayanıklı metallerden biri; buji elektrotlarında ve tarihsel olarak metre/kilogram prototip standartlarında kullanılmıştır."),
    78:  ("Pt", "Platin",        "gecis",        "Katalitik konvertörlerde, laboratuvar ekipmanlarında ve takıda kullanılan değerli bir metaldir."),
    79:  ("Au", "Altın",         "gecis",        "Yüksek korozyon direnci ve iletkenliği nedeniyle takı, elektronik kontak ve rezerv para birimi olarak kullanılır."),
    80:  ("Hg", "Civa",          "gecis",        "Oda sıcaklığında sıvı hâlde bulunan tek metal; geleneksel termometre ve barometrelerde kullanılırdı."),
    81:  ("Tl", "Talyum",        "post_gecis",   "Oldukça toksik bir elementtir; Tl-201 izotopu kalp kası sintigrafisinde tıbbi görüntüleme amacıyla kullanılır."),
    82:  ("Pb", "Kurşun",        "post_gecis",   "Yüksek yoğunluğu nedeniyle radyasyon zırhlamada ve akülerde kullanılır; birçok doğal bozunma serisinin kararlı son ürünüdür."),
    83:  ("Bi", "Bizmut",        "post_gecis",   "Uzun süre 'en ağır kararlı element' sanılmıştır — aslında son derece uzun ömürlü (yarı ömrü evrenin yaşının milyarlarca katı) bir α-yayıcıdır; mide ilaçlarında kullanılır."),
    84:  ("Po", "Polonyum",      "metaloid",     "Marie Curie tarafından 1898'de keşfedilmiştir (Polonya'ya ithafen adlandırılmıştır); yoğun α-radyoaktivitesiyle bilinir."),
    85:  ("At", "Astatin",       "halojen",      "Doğada en nadir bulunan elementlerden biridir; tüm izotopları radyoaktif ve kısa ömürlüdür."),
    86:  ("Rn", "Radon",         "soy_gaz",      "Radyoaktif bir soy gazdır; toprak ve kayalardan sızarak kapalı mekanlarda birikebilir ve akciğer kanseri riski oluşturur."),
    87:  ("Fr", "Fransiyum",     "alkali",       "Doğada bulunan en kararsız elementlerden biridir; en uzun ömürlü izotopunun yarı ömrü bile dakikalar mertebesindedir."),
    88:  ("Ra", "Radyum",        "alkali_toprak","Marie ve Pierre Curie tarafından keşfedilmiştir; eskiden ışıldayan boyalarda kullanılmış, günümüzde radyasyon güvenliği açısından dikkatle ele alınır."),
    89:  ("Ac", "Aktinyum",      "aktinit",      "Aktinit serisine adını veren element; güçlü α-yayıcı olması nedeniyle hedefe yönelik radyoterapide (Ac-225) araştırılmaktadır."),
    90:  ("Th", "Toryum",        "aktinit",      "Th-232, potansiyel bir alternatif nükleer yakıt döngüsünün (toryum-uranyum-233 döngüsü) başlangıç malzemesidir."),
    91:  ("Pa", "Protaktinyum",  "aktinit",      "Doğal uranyum-235 bozunma serisinde ara ürün olarak çok az miktarda doğada bulunur."),
    92:  ("U",  "Uranyum",       "aktinit",      "Nükleer reaktörlerin ve bazı silahların temel yakıt elementidir; U-235 termal nötronlarla fisyona uğrayabilen (fisil) doğal izotoptur."),
    93:  ("Np", "Neptünyum",     "aktinit",      "Uranyumdan sonra keşfedilen ilk transuranyum element; plütonyum üretim zincirinde ara ürün olarak oluşur."),
    94:  ("Pu", "Plütonyum",     "aktinit",      "Pu-239, hem nükleer reaktör yakıtı hem de nükleer silahlarda kullanılan fisil bir izotoptur."),
    95:  ("Am", "Amerikyum",     "aktinit",      "Am-241 izotopu, evlerde kullanılan iyonizasyon tipi duman dedektörlerinin radyasyon kaynağıdır."),
    96:  ("Cm", "Küriyum",       "aktinit",      "Marie ve Pierre Curie onuruna adlandırılmıştır; uzay araçlarında radyoizotop güç kaynağı olarak kullanılmıştır."),
    97:  ("Bk", "Berkelyum",     "aktinit",      "Berkeley, Kaliforniya'da sentezlenmiştir; yalnızca araştırma amaçlı iz miktarlarda üretilir."),
    98:  ("Cf", "Kaliforniyum",  "aktinit",      "Güçlü bir nötron kaynağıdır; endüstriyel nötron radyografisinde ve reaktör başlatma kaynaklarında kullanılır."),
    99:  ("Es", "Aynştaynyum",   "aktinit",      "Albert Einstein onuruna adlandırılmıştır; ilk kez 1952'deki termonükleer patlama kalıntılarında tespit edilmiştir."),
    100: ("Fm", "Fermiyum",      "aktinit",      "Enrico Fermi onuruna adlandırılmıştır; yalnızca parçacık hızlandırıcılarında iz miktarlarda üretilebilir."),
    101: ("Md", "Mendelevyum",   "aktinit",      "Periyodik tablonun yaratıcısı Dmitri Mendeleyev onuruna adlandırılmıştır."),
    102: ("No", "Nobelyum",      "aktinit",      "Alfred Nobel onuruna adlandırılmıştır; tüm izotopları kısa ömürlü ve yalnızca laboratuvarda üretilir."),
    103: ("Lr", "Lawrensiyum",   "aktinit",      "Parçacık hızlandırıcısının (siklotron) mucidi Ernest Lawrence onuruna adlandırılmıştır."),
    104: ("Rf", "Rutherfordiyum","super_agir",   "Atom çekirdeğini keşfeden Ernest Rutherford onuruna adlandırılmıştır; ağır iyon hızlandırıcılarında sentezlenir."),
    105: ("Db", "Dubniyum",      "super_agir",   "Rusya'daki Dubna Birleşik Nükleer Araştırmalar Enstitüsü onuruna adlandırılmıştır."),
    106: ("Sg", "Seaborgiyum",   "super_agir",   "Glenn Seaborg onuruna adlandırılmıştır — hayattayken bir elementin kendi adıyla anıldığı nadir örneklerden biridir."),
    107: ("Bh", "Bohriyum",      "super_agir",   "Fizikçi Niels Bohr onuruna adlandırılmıştır."),
    108: ("Hs", "Hassiyum",      "super_agir",   "Almanya'nın Hessen eyaleti onuruna adlandırılmıştır; GSI Darmstadt'ta keşfedilmiştir."),
    109: ("Mt", "Meitneryum",    "super_agir",   "Nükleer fizyonun teorik açıklamasına katkı sağlayan Lise Meitner onuruna adlandırılmıştır."),
    110: ("Ds", "Darmstadtiyum", "super_agir",   "Keşfedildiği şehir olan Almanya'nın Darmstadt kentine ithafen adlandırılmıştır."),
    111: ("Rg", "Röntgenyum",    "super_agir",   "X-ışınlarını keşfeden Wilhelm Röntgen onuruna adlandırılmıştır."),
    112: ("Cn", "Kopernikyum",   "super_agir",   "Gök bilimci Nicolaus Copernicus onuruna adlandırılmıştır."),
    113: ("Nh", "Nihonyum",      "super_agir",   "Japonya'da (Nihon) keşfedilmiştir; Asya'da keşfedilip adlandırılan ilk elementtir."),
    114: ("Fl", "Flerovyum",     "super_agir",   "Rusya'daki Flerov Nükleer Reaksiyonlar Laboratuvarı onuruna adlandırılmıştır."),
    115: ("Mc", "Moskovyum",     "super_agir",   "Rusya'nın Moskova bölgesi onuruna adlandırılmıştır; Dubna'da sentezlenmiştir."),
    116: ("Lv", "Livermoryum",   "super_agir",   "ABD'deki Lawrence Livermore Ulusal Laboratuvarı onuruna adlandırılmıştır."),
    117: ("Ts", "Tennessin",     "super_agir",   "ABD'nin Tennessee eyaleti onuruna adlandırılmıştır; periyodik tabloda keşfedilen en son halojendir."),
    118: ("Og", "Oganesson",     "super_agir",   "Fizikçi Yuri Oganessian onuruna adlandırılmıştır; bilinen en ağır elementtir."),
}

CATEGORY_LABELS = {
    "alkali": "Alkali Metal",
    "alkali_toprak": "Toprak Alkali Metal",
    "gecis": "Geçiş Metali",
    "post_gecis": "Geçiş Sonrası Metal",
    "metaloid": "Metaloid",
    "ametal": "Ametal",
    "halojen": "Halojen",
    "soy_gaz": "Soy Gaz",
    "lantanit": "Lantanit",
    "aktinit": "Aktinit",
    "super_agir": "Süper Ağır Element",
}


def get_element_info(Z: int):
    """Z proton sayısı için (sembol, Türkçe ad, kategori_etiketi, açıklama) döndürür.
    Bilinmeyen Z için None döner (mevcut kapsam 1..118)."""
    entry = ELEMENT_INFO.get(Z)
    if not entry:
        return None
    sym, tr_name, cat_code, desc = entry
    return {
        "symbol": sym,
        "tr_name": tr_name,
        "category": CATEGORY_LABELS.get(cat_code, cat_code),
        "desc": desc,
    }
