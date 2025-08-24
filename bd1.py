import requests
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.Seq import Seq
from Bio import pairwise2
import os

# === 配置 ===
pdb_ids = ['3MXF', '6SWN', '2YEL', '3U5L', '3ZYU', '4BW1', '4CFL', '5BT4', '5CQT', '5F62', '6CD4', '6E4A', '6ZB3']  # 示例，你可以替换
bd1_sequence = 'WPFQQPVDAVKLNLPDYYKIIKTPMDMGTIKKRLENNYYWNAQECIQDFNTMFTNCYIYNKPGDDIVLMAEA'  # 示例序列（你提供完整的）

matched = []

def download_cif(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    r = requests.get(url)
    if r.status_code == 200:
        with open(f"{pdb_id}.cif", "wb") as f:
            f.write(r.content)
        return True
    return False

def extract_sequences(pdb_id):
    parser = MMCIFParser()
    structure = parser.get_structure(pdb_id, f"{pdb_id}.cif")
    sequences = {}
    for model in structure:
        for chain in model:
            seq = ''
            for res in chain:
                if res.get_id()[0] == ' ':
                    if 'CA' in res:
                        resname = res.get_resname()
                        try:
                            aa = Seq('').three_letter_to_one_letter(resname)
                        except:
                            from Bio.Data.IUPACData import protein_letters_3to1
                            aa = protein_letters_3to1.get(resname, 'X')
                        seq += aa
            if len(seq) > 0:
                sequences[chain.id] = seq
    return sequences

def is_match(chain_seq, bd1_seq):
    alignments = pairwise2.align.localms(chain_seq, bd1_seq, 2, -1, -0.5, -0.1)
    if alignments:
        top = alignments[0]
        identity = top[2] / len(bd1_seq)
        return identity > 0.9  # 可调节匹配阈值
    return False

# === 主程序 ===
for pdb in pdb_ids:
    print(f"Processing {pdb}...")
    if not os.path.exists(f"{pdb}.cif"):
        success = download_cif(pdb)
        if not success:
            print(f"Failed to download {pdb}")
            continue
    seqs = extract_sequences(pdb)
    for chain_id, seq in seqs.items():
        if is_match(seq, bd1_sequence):
            print(f"Match found in {pdb} (Chain {chain_id})")
            matched.append(pdb)
            break

print("\n✅ Matched complexes:")
print(matched)
