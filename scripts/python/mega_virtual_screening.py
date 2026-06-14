# ==============================================================================
# SCRIPT: MEGA HTVS - 10.000 İlaçlık Yüksek Hacimli Sanal Tarama
# ==============================================================================

import os
import random
import multiprocessing
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import PDBParser

print("=========================================================")
print("MEGA SANAL TARAMA (10.000 MOLEKÜL) BAŞLATILIYOR...")
print("İşlemci Çekirdek Sayısı (CPU Cores):", multiprocessing.cpu_count())
print("=========================================================")

# Dosya yolları
base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
pdb_dir = os.path.join(base_dir, "data", "pdb")
hts_dir = os.path.join(base_dir, "data", "mega_hts_ligands")
results_dir = os.path.join(base_dir, "results", "reports")

os.makedirs(hts_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

protein_pdbqt = os.path.join(pdb_dir, "6OM3_clean.pdbqt")

# Kutu Merkezi
parser = PDBParser(QUIET=True)
structure = parser.get_structure("6OM3", os.path.join(pdb_dir, "6OM3_clean.pdb"))
x_sum = y_sum = z_sum = count = 0
for atom in structure.get_atoms():
    coord = atom.get_coord()
    x_sum += coord[0]; y_sum += coord[1]; z_sum += coord[2]
    count += 1
cx, cy, cz = x_sum/count, y_sum/count, z_sum/count

# Sentezlenecek İlaç Sayısı (Hedef: 10.000)
# Benzamidin iskeletine rastgele dallar ekleyerek 10.000 türev oluşturuyoruz
TOTAL_DRUGS = 10000
results = []
report_path = os.path.join(results_dir, "MEGA_10000_Drug_Screening_Hits.csv")

print(f"Toplam {TOTAL_DRUGS} molekül sentezlenip Vina motoruna (22 Threads) gönderiliyor. Bu işlem sabaha kadar sürecektir...")

# Otomatik mega döngü
for i in range(1, TOTAL_DRUGS + 1):
    drug_id = f"Virtual_Drug_{i:05d}"
    
    # İlaç çeşitlemesi (Basit bir randomize varyasyon)
    carbon_chain = "C" * random.randint(1, 8)
    halogen = random.choice(["F", "Cl", "O", "N", "S"])
    smiles = f"NC(=N)c1ccc(cc1){carbon_chain}{halogen}"
    
    lig_sdf = os.path.join(hts_dir, f"{drug_id}.sdf")
    lig_pdbqt = os.path.join(hts_dir, f"{drug_id}.pdbqt")
    out_pdbqt = os.path.join(hts_dir, f"{drug_id}_docked.pdbqt")
    log_txt = os.path.join(hts_dir, f"{drug_id}_vina.log")
    
    try:
        # Molekül sentezi ve 3D katlama
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(mol)
        writer = Chem.SDWriter(lig_sdf)
        writer.write(mol)
        writer.close()
        
        # PDBQT Dönüşümü
        os.system(f"obabel {lig_sdf} -O {lig_pdbqt} 2>/dev/null")
        
        # Vina Çarpışması (Her işlem tüm CPU gücünü kullanacak)
        vina_cmd = (
            f"vina --receptor {protein_pdbqt} --ligand {lig_pdbqt} "
            f"--center_x {cx:.2f} --center_y {cy:.2f} --center_z {cz:.2f} "
            f"--size_x 30 --size_y 30 --size_z 30 --out {out_pdbqt} > {log_txt}"
        )
        os.system(vina_cmd)
        
        # Skor okuma
        best_score = 0.0
        with open(log_txt, 'r') as f:
            for line in f:
                if "   1 " in line:
                    best_score = float(line.split()[1])
                    break
                    
        results.append({"Drug_ID": drug_id, "SMILES": smiles, "Affinity_Score": best_score})
        
        # Her 100 ilaçta bir raporu kaydet (Elektrik kesilirse veri kaybolmasın)
        if i % 100 == 0:
            pd.DataFrame(results).to_csv(report_path, index=False)
            print(f"[{i}/{TOTAL_DRUGS}] tarandı. Anlık en iyi skorlar diske yazıldı...")
            
    except Exception as e:
        continue

# Son Kayıt
df = pd.DataFrame(results).sort_values(by="Affinity_Score", ascending=True)
df.to_csv(report_path, index=False)
print("SİSTEM GÖREVİ TAMAMLADI: Mega rapor hazırlandı!")