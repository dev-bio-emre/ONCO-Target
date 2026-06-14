# ==============================================================================
# SCRIPT: AutoDock Vina ile Moleküler Kenetlenme (Docking) Simülasyonu
# ==============================================================================

import os
import subprocess
from Bio.PDB import PDBParser

# Dosya yolları
base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
pdb_dir = os.path.join(base_dir, "data", "pdb")
lig_dir = os.path.join(base_dir, "data", "ligands")

protein_pdb = os.path.join(pdb_dir, "6OM3_clean.pdb")
protein_pdbqt = os.path.join(pdb_dir, "6OM3_clean.pdbqt")

ligand_sdf = os.path.join(lig_dir, "Benzamidine_Analog.sdf")
ligand_pdbqt = os.path.join(lig_dir, "Benzamidine_Analog.pdbqt")

output_pdbqt = os.path.join(lig_dir, "Docking_Result.pdbqt")

print("1. PDBQT Format Dönüşümleri Yapılıyor (OpenBabel ile)...")
# Protein için: -xr (esnek yan zincirleri sabitle)
os.system(f"obabel {protein_pdb} -O {protein_pdbqt} -xr 2>/dev/null")
# Ligand için
os.system(f"obabel {ligand_sdf} -O {ligand_pdbqt} 2>/dev/null")

print("2. Arama Kutusu (Grid Box) Merkez Koordinatları Hesaplanıyor...")
# Proteinin merkezini buluyoruz ki Vina nereye çarpacağını bilsin (Blind Docking)
parser = PDBParser(QUIET=True)
structure = parser.get_structure("6OM3", protein_pdb)

x_sum = y_sum = z_sum = 0
count = 0
for atom in structure.get_atoms():
    coord = atom.get_coord()
    x_sum += coord[0]
    y_sum += coord[1]
    z_sum += coord[2]
    count += 1

center_x = x_sum / count
center_y = y_sum / count
center_z = z_sum / count

print(f"Kutu Merkezi -> X: {center_x:.2f}, Y: {center_y:.2f}, Z: {center_z:.2f}")

print("3. AutoDock Vina Motoru Ateşleniyor (Bu işlem işlemcine bağlı olarak 15-30 saniye sürebilir)...")
# 30x30x30 Angstromluk devasa bir arama kutusu
vina_cmd = (
    f"vina --receptor {protein_pdbqt} --ligand {ligand_pdbqt} "
    f"--center_x {center_x:.2f} --center_y {center_y:.2f} --center_z {center_z:.2f} "
    f"--size_x 30 --size_y 30 --size_z 30 --out {output_pdbqt}"
)

# Vina'yı çalıştır
os.system(vina_cmd)

print("=========================================================")
print("SİMÜLASYON TAMAMLANDI! Sonuçlar ekrana yazdırıldı.")
print("=========================================================")