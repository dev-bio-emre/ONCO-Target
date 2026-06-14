# ==============================================================================
# SCRIPT: Ham PDB Yapısını Kenetlenme (Docking) İçin Temizleme
# ==============================================================================

import os
from Bio.PDB import PDBParser, PDBIO, Select

target_pdb_id = "6OM3"

# Dosya yolları
base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
raw_pdb = os.path.join(base_dir, "data", "pdb", f"{target_pdb_id}.pdb")
clean_pdb = os.path.join(base_dir, "data", "pdb", f"{target_pdb_id}_clean.pdb")

print(f"[{target_pdb_id}] Ham protein yapısı temizleme motoruna alınıyor...")

# Biyoinformatik Filtresi: Sadece 'A' zincirini tut, su (HOH) ve yabancı iyonları at!
class ReceptorSelect(Select):
    def accept_residue(self, residue):
        # Eğer kalıntı (residue) 'A' zincirindeyse VE bir heteroatom (su, ligand vb.) değilse
        if residue.get_parent().id == 'A' and residue.id[0] == ' ':
            return True
        return False

try:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(target_pdb_id, raw_pdb)
    
    # Yeni ve temiz yapıyı diske yazıcı (IO) objesi
    io = PDBIO()
    io.set_structure(structure)
    
    # Filtremizi (ReceptorSelect) uygulayarak kaydet
    io.save(clean_pdb, ReceptorSelect())
    
    # Temizlenmiş yapıyı test edip atom sayısını okuyalım
    clean_structure = parser.get_structure(f"{target_pdb_id}_clean", clean_pdb)
    atom_count = sum(1 for _ in clean_structure.get_atoms())
    
    print("-" * 50)
    print("TEMİZLİK BAŞARILI: Su molekülleri ve kopya zincirler buharlaştırıldı!")
    print(f"Orijinal Atom Sayısı : 28906")
    print(f"Kalan Saf Atom Sayısı: {atom_count}")
    print(f"Temiz Dosya Konumu   : {clean_pdb}")
    print("-" * 50)

except Exception as e:
    print(f"HATA: Temizleme işlemi sırasında bir pürüz çıktı! Detay: {e}")