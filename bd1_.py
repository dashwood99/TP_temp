import requests
from collections import OrderedDict

# 需要检查的 PDB 条目
PDBS = ['3MXF', '6SWN', '2YEL', '3U5L', '3ZYU', '4BW1', '4CFL', '5BT4', '5CQT', '5F62', '6CD4', '6E4A', '6ZB3']

# 你的目标序列（去掉换行空格）
target = ("WPFQQPVDAVKLNLPDYYKIIKTPMDMGTIKKRLENNYYWNAQECIQDFNTMFTNCYIYNKPGDDIVLMAEA"
          ).strip().upper()

def fetch_entry_fasta(pdb_id: str) -> str:
    # RCSB 文档说明的新版下载端点（per entry）
    # 也可改用 /fasta/entity/{PDB}_{entity_id}/download 或 /fasta/chain/{PDB}.{asym_id}/download
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/download"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.text

def parse_only_proteins(fasta_text: str):
    """
    只保留蛋白（polypeptide）实体。
    RCSB 的 FASTA header 里会带有类似：
      >2YEM_1|Chains A,B|BROMODOMAIN-CONTAINING PROTEIN 4 ...
    我们将每个条目展开到链层级（A、B…）存成字典。
    """
    out = OrderedDict()
    header, seq = None, []
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            if header and seq:
                out[header] = "".join(seq).replace(" ", "").upper()
            # 解析 header，抽出 PDB、实体与链
            h = line[1:].strip()
            # 形如 "2YEM_1|Chains A,B|BROMODOMAIN-CONTAINING PROTEIN 4"
            parts = h.split("|")
            id_part = parts[0]            # e.g. 2YEM_1
            chains = []
            for p in parts[1:]:
                if p.strip().lower().startswith("chains"):
                    # "Chains A,B" → ["A","B"]
                    chains = [c.strip() for c in p.split(" ",1)[1].split(",")]
                    break
            # 没明确写 Chains 的（有些条目可能写 Instance/Label），就只用实体级键
            if chains:
                # 用实体级临时键，等拿到序列后再复制到每条链
                header = (id_part, tuple(chains))  # 先暂存
            else:
                header = (id_part, tuple())
            seq = []
        else:
            seq.append(line.strip())

    if header and seq:
        out[header] = "".join(seq).replace(" ", "").upper()

    # 展开：一个实体对应多条链 → 为每条链复制同一条蛋白序列
    expanded = OrderedDict()
    for (id_part, chains), sequence in out.items():
        pdb, entity = id_part.split("_", 1) if "_" in id_part else (id_part, "")
        if chains:
            for ch in chains:
                expanded[f"{pdb}_{entity}_chain_{ch}"] = sequence
        else:
            expanded[f"{pdb}_{entity}"] = sequence
    return expanded

def build_seq_dict(pdb_list):
    seq_dict = OrderedDict()
    for pdb in pdb_list:
        fa = fetch_entry_fasta(pdb)
        seqs = parse_only_proteins(fa)
        seq_dict.update(seqs)
    return seq_dict

def find_subseq(seq, subseq):
    """返回所有出现的位置（1-based 起止）。"""
    hits = []
    i = 0
    while True:
        j = seq.find(subseq, i)
        if j == -1: break
        hits.append((j+1, j+len(subseq)))  # 1-based
        i = j + 1
    return hits

# 1) 拉取并构建 {键: 序列} 的字典
seq_dict = build_seq_dict(PDBS)

# 2) 在每条序列里搜索你的目标序列
results = {}
for key, seq in seq_dict.items():
    hits = find_subseq(seq, target)
    if hits:
        results[key] = hits

print("共获取到的蛋白序列条目：", len(seq_dict))
for k in list(seq_dict.keys())[:5]:
    print("示例键：", k, "长度=", len(seq_dict[k]))

print("\n比对结果（出现则列出起止位置，1-based）:")
if results:
    for k, hitlist in results.items():
        print(k, "=>", hitlist)
else:
    print("未在任何序列中找到完整子序列。")

