import os
import io
import sys
import gzip
import requests

# ------- 可调参数 -------
PDBS = ['3MXF', '6SWN', '2YEL', '3U5L', '3ZYU', '4BW1', '4CFL', '5BT4', '5CQT', '5F62', '6CD4', '6E4A', '6ZB3']
OUT_DIR = "filtered_pdbs"
EXCLUDE_WATER = True
KEEP_IONS = True  # 若想排除离子，把它设为 False 并补充 ION_SET
ION_SET = {"NA", "K", "CL", "MG", "CA", "ZN", "MN", "FE", "CU", "CO", "CD", "NI", "BR", "I", "F", "SR", "BA"}

# ------- 依赖检查 -------
try:
    from Bio.PDB import PDBParser, PDBIO, Select, is_aa
except Exception as e:
    print("需要安装 Biopython： pip install biopython")
    raise

os.makedirs(OUT_DIR, exist_ok=True)

def download_pdb(pdb_id: str) -> str:
    """
    下载 PDB 文本（优先 .pdb.gz，不行则 .pdb），返回字符串内容。
    """
    pdb_id = pdb_id.upper()
    urls = [
        f"https://files.rcsb.org/download/{pdb_id}.pdb.gz",
        f"https://files.rcsb.org/download/{pdb_id}.pdb",
    ]
    last_err = None
    for url in urls:
        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            data = r.content
            if url.endswith(".gz"):
                try:
                    data = gzip.decompress(data)
                except OSError:
                    # 不是有效 gzip，直接当作文本
                    pass
            txt = data.decode("utf-8", errors="replace")
            return txt
        except Exception as e:
            last_err = e
    raise RuntimeError(f"下载 {pdb_id} 失败：{last_err}")

def find_protein_chains(structure):
    """
    返回蛋白质链 ID 列表（按出现顺序）。
    判定标准：链中含有至少一个标准氨基酸残基。
    """
    protein_chains = []
    for model in structure:
        for chain in model:
            # 只要该链中出现至少一个标准氨基酸，就视为蛋白链
            if any(is_aa(res, standard=True) for res in chain):
                if chain.id not in protein_chains:
                    protein_chains.append(chain.id)
        break  # 只看第一个 model 即可
    return protein_chains

class ProteinAPlusLigandSelect(Select):
    """
    仅保留：
      - 选中的蛋白质链（仅氨基酸残基）；
      - 所有 HETATM 小分子（可选排除水/离子）。
    """
    def __init__(self, protein_chain_id, exclude_water=True, keep_ions=True, ion_set=None):
        super().__init__()
        self.protein_chain_id = protein_chain_id
        self.exclude_water = exclude_water
        self.keep_ions = keep_ions
        self.ion_set = ion_set or set()

    def accept_model(self, model):
        # 保留所有 model（通常只有 model 0）
        return 1

    def accept_chain(self, chain):
        # 不在链层面丢弃；在残基层面过滤
        return 1

    def accept_residue(self, residue):
        hetflag, resseq, icode = residue.id  # hetflag: ' ' (标准) 或 'H_' (HETATM)
        resname = residue.get_resname().strip().upper()

        # 氨基酸：只保留在目标蛋白链
        if hetflag == " ":
            if is_aa(residue, standard=True) and residue.get_parent().id == self.protein_chain_id:
                return 1
            else:
                return 0

        # HETATM（小分子/辅因子/离子/水）
        # 排除水
        if self.exclude_water and resname in {"HOH", "WAT"}:
            return 0
        # 若不保留离子，则用一个简单的离子集合进行排除
        if not self.keep_ions and resname in self.ion_set:
            return 0
        return 1

    def accept_atom(self, atom):
        return 1

def filter_and_save(pdb_id: str):
    pdb_text = download_pdb(pdb_id)
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    structure = parser.get_structure(pdb_id.upper(), io.StringIO(pdb_text))

    protein_chains = find_protein_chains(structure)
    if not protein_chains:
        print(f"[{pdb_id}] 未识别到蛋白质链，跳过。")
        return

    # 选 A 链，若不存在则取第一个蛋白链
    if "A" in protein_chains:
        protein_chain_id = "A"
    else:
        protein_chain_id = protein_chains[0]

    selector = ProteinAPlusLigandSelect(
        protein_chain_id=protein_chain_id,
        exclude_water=EXCLUDE_WATER,
        keep_ions=KEEP_IONS,
        ion_set=ION_SET
    )
    io_writer = PDBIO()
    io_writer.set_structure(structure)

    out_path = os.path.join(OUT_DIR, f"{pdb_id.upper()}_cleaned.pdb")
    io_writer.save(out_path, select=selector)
    print(f"[{pdb_id}] 蛋白链选择：{protein_chain_id}  → 已保存：{out_path}")

if __name__ == "__main__":
    for pid in PDBS:
        try:
            filter_and_save(pid)
        except Exception as e:
            print(f"[{pid}] 处理失败：{e}", file=sys.stderr)

