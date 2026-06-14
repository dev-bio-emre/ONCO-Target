# ==============================================================================
# SCRIPT: PyMOL ile Protein-Ligand Kenetlenme (Docking) Görselleştirmesi
# ==============================================================================

# 1. Dosya yolları ana dizine (çalıştırma noktasına) göre ayarlandı
load data/pdb/6OM3_clean.pdb, receptor
load data/ligands/Docking_Result.pdbqt, ligand

# 2. Arka planı beyaz yap ve genel görünümü temizle
bg_color white
hide everything

# 3. Proteini 3 boyutlu "kurdela" (cartoon) modelinde göster ve rengini gri yap
show cartoon, receptor
color gray70, receptor

# 4. İlacımızı (Ligand) "çubuk" (sticks) modelinde, belirgin bir yeşil renkle göster
show sticks, ligand
color green, ligand

# 5. Protein cebindeki atomları (liganda 5 Angstrom yakınlıktakileri) göster
select pocket, byres receptor within 5.0 of ligand
show sticks, pocket
color cyan, pocket

# Karbon, Azot ve Oksijenleri standart renklerine boya (Yorum satırı alt satıra alındı)
util.cnc pocket  

# 6. Ligand ve Protein arasındaki Hidrojen bağlarını (H-bonds) hesapla ve sarı kesik çizgilerle çiz
distance h_bonds, receptor, ligand, 3.2, mode=2
set dash_color, yellow
set dash_gap, 0.3
set dash_radius, 0.05

# 7. Kamerayı tam olarak savaş alanına (ligandın olduğu cebe) odakla ve yakınlaş
center ligand
zoom ligand, 12

# 8. Yüksek kaliteli görüntü için ışık ve gölgeleri aç
