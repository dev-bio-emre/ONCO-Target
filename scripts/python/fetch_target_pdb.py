# ==============================================================================
# SCRIPT: RCSB PDB Veritabanından Hedef Proteinin 3 Boyutlu Yapısını İndirme
# ==============================================================================

import os
import urllib.request
from Bio.PDB import PDBParser

# Hedef proteinimiz: KLK8 (PDB ID: 6OM3)
target_pdb_id = "6OM3"

# Dosyaların kaydedileceği klasör yapısını kur
base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
pdb_dir = os.path.join(base_dir, "data", "pdb")
os.makedirs(pdb_dir, exist_ok=True)

pdb_file_path = os.path.join(pdb_dir, f"{target_pdb_id}.pdb")

print(f"[{target_pdb_id}] kristal yapısı PDB sunucularından aranıyor...")

# PDB dosyasını doğrudan indir
url = f"https://files.rcsb.org/download/{target_pdb_id}.pdb"

try:
    urllib.request.urlretrieve(url, pdb_file_path)
    print(f"BAŞARILI: {target_pdb_id}.pdb dosyası diske kaydedildi -> {pdb_file_path}")
    
    # Biopython ile yapıyı test edip atom sayılarını okuyalım
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(target_pdb_id, pdb_file_path)
    
    atom_count = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_count += 1
                    
    print("-" * 50)
    print(f"PROTEİN YAPISI DOĞRULANDI!")
    print(f"Toplam tespit edilen atom sayısı: {atom_count}")
    print("-" * 50)
    print("Sistem moleküler kenetlenme (Molecular Docking) için hazır.")

except Exception as e:
    print(f"HATA: Dosya indirilirken bir sorun oluştu! Detay: {e}")