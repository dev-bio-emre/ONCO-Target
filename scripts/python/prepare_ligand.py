# ==============================================================================
# SCRIPT: SMILES Formatından 3D Ligand (İlaç Molekülü) Sentezi
# ==============================================================================

import os
from rdkit import Chem
from rdkit.Chem import AllChem

# Benzamidin iskeleti (Serin proteazlar için klasik bir inhibitör başlangıcı)
smiles_string = "NC(=N)c1ccc(cc1)"
ligand_name = "Benzamidine_Analog"

# Kayıt yolları
base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ligand_dir = os.path.join(base_dir, "data", "ligands")
os.makedirs(ligand_dir, exist_ok=True)

out_sdf = os.path.join(ligand_dir, f"{ligand_name}.sdf")
out_pdb = os.path.join(ligand_dir, f"{ligand_name}.pdb")

print(f"SMILES şifresi çözülüyor: {smiles_string}")

try:
    # 1. SMILES'tan 2D molekül objesi oluştur ve hidrojenleri ekle
    mol = Chem.MolFromSmiles(smiles_string)
    mol = Chem.AddHs(mol)
    
    # 2. Molekülü 3 Boyutlu uzaya katla (Konformasyon optimizasyonu)
    print("Molekül 3D uzayda katlanıyor ve enerji optimizasyonu yapılıyor (MMFF94)...")
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)
    
    # 3. İki farklı formatta kaydet (SDF ve PDB)
    writer_sdf = Chem.SDWriter(out_sdf)
    writer_sdf.write(mol)
    writer_sdf.close()
    
    Chem.MolToPDBFile(mol, out_pdb)
    
    print("-" * 50)
    print("LİGAND BAŞARIYLA SENTEZLENDİ!")
    print(f"Toplam Atom Sayısı (Hidrojenler dahil): {mol.GetNumAtoms()}")
    print(f"SDF Konumu: {out_sdf}")
    print(f"PDB Konumu: {out_pdb}")
    print("-" * 50)

except Exception as e:
    print(f"HATA: Molekül sentezlenirken bir sorun oluştu! Detay: {e}")