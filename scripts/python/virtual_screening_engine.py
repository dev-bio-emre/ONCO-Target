# ==============================================================================
# SCRIPT: Yüksek Hacimli Sanal Tarama (High-Throughput Virtual Screening)
# ==============================================================================

import os
import subprocess
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import PDBParser

# 1. Mini FDA Onaylı İlaç Kütüphanesi (SMILES Formatında)
# İlaç Yeniden Konumlandırma (Repurposing) için bilindik kanser/klinik ilaçları
fda_library = {
    "Erlotinib": "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC",  # Pankreas/Akciğer Kanseri ilacı
    "Imatinib": "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CC=CN=C5", # Lösemi ilacı
    "Fluorouracil": "C1=C(C(=O)NC(=O)N1)F", # Klasik Kemoterapi
    "Metformin": "CN(C)C(=N)N=C(N)N", # Diyabet ilacı (Kanser metabolizmasını etkiler mi?)
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O" # Ağrı kesici (Kontrol grubu)
}

# 2. Dosya Yolları
base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
pdb_dir = os.path.join(base_dir, "data", "pdb")
hts_dir = os.path.join(base_dir, "data", "hts_ligands")
results_dir = os.path.join(base_dir, "results", "reports")

os.makedirs(hts_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

protein_pdb = os.path.join(pdb_dir, "6OM3_clean.pdb")
protein_pdbqt = os.path.join(pdb_dir, "6OM3_clean.pdbqt")

# 3. Kutu Merkezini Hesapla
print("Hedefin katalitik merkezi hesaplanıyor...")
parser = PDBParser(QUIET=True)
structure = parser.get_structure("6OM3", protein_pdb)
x_sum = y_sum = z_sum = count = 0
for atom in structure.get_atoms():
    coord = atom.get_coord()
    x_sum += coord[0]; y_sum += coord[1]; z_sum += coord[2]
    count += 1
cx, cy, cz = x_sum/count, y_sum/count, z_sum/count

print("=========================================================")
print("SANAL TARAMA (VIRTUAL SCREENING) FABRİKASI ÇALIŞTIRILIYOR")
print("=========================================================")

results = []

# 4. Otomatik Tarama Döngüsü
for drug_name, smiles in fda_library.items():
    print(f"\n---> Hedefteki İlaç: {drug_name} işleniyor...")
    
    lig_sdf = os.path.join(hts_dir, f"{drug_name}.sdf")
    lig_pdbqt = os.path.join(hts_dir, f"{drug_name}.pdbqt")
    out_pdbqt = os.path.join(hts_dir, f"{drug_name}_docked.pdbqt")
    log_txt = os.path.join(hts_dir, f"{drug_name}_vina.log")
    
    try:
        # A. SMILES -> 3D SDF Dönüşümü (RDKit)
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(mol)
        writer = Chem.SDWriter(lig_sdf)
        writer.write(mol)
        writer.close()
        
        # B. SDF -> PDBQT Dönüşümü (OpenBabel)
        os.system(f"obabel {lig_sdf} -O {lig_pdbqt} 2>/dev/null")
        
        # C. AutoDock Vina ile Çarpışma (Docking)
        vina_cmd = (
            f"vina --receptor {protein_pdbqt} --ligand {lig_pdbqt} "
            f"--center_x {cx:.2f} --center_y {cy:.2f} --center_z {cz:.2f} "
            f"--size_x 30 --size_y 30 --size_z 30 --out {out_pdbqt} > {log_txt}"
        )
        os.system(vina_cmd)
        
        # D. Sonuç Log'unu Okuyup En İyi Skoru (Affinity) Çekmek
        best_score = 0.0
        with open(log_txt, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if "   1 " in line:  # Mode 1 satırını bul
                    parts = line.split()
                    best_score = float(parts[1])
                    break
        
        print(f"[BAŞARILI] {drug_name} KLK8'e çarptı. Skor: {best_score} kcal/mol")
        results.append({"Drug_Name": drug_name, "Affinity_Score": best_score})
        
    except Exception as e:
        print(f"[HATA] {drug_name} işlenirken hata oluştu: {e}")

# 5. Raporlama
print("\n=========================================================")
print("TARAMA BİTTİ! LİDER TABLOSU OLUŞTURULUYOR...")
df = pd.DataFrame(results)
# Skoru en düşük (negatif) olan en iyisidir
df = df.sort_values(by="Affinity_Score", ascending=True).reset_index(drop=True)

report_path = os.path.join(results_dir, "FDA_Drug_Repurposing_Top_Hits.csv")
df.to_csv(report_path, index=False)

print(f"Büyük Rapor Kaydedildi: {report_path}")
print(df)
print("=========================================================")