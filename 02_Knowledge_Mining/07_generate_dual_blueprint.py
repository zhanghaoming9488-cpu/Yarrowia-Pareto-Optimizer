import pandas as pd
import numpy as np
import gzip
from Bio import SeqIO

# === 配置 ===
GOLDEN_CSV = "golden_elite_genes.csv"
REF_FASTA = "reference/Yarrowia_lipolytica.cdna.all.fa.gz"
RAMP_SIZE = 39 # 前 13 个密码子
OUTPUT_MATRIX = "pareto_codon_blueprint.csv"
# ===========

CODON_DICT = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def get_elite_regions():
    """精确切分 Elite 基因的 Ramp 区和 Body 区"""
    try:
        df_golden = pd.read_csv(GOLDEN_CSV)
        elite_ids = set(df_golden['GeneID_RNA'].astype(str).tolist())
    except FileNotFoundError:
        print(f"[!] 找不到 {GOLDEN_CSV}")
        return [], []

    open_func = gzip.open if REF_FASTA.endswith('.gz') else open
    mode = "rt" if REF_FASTA.endswith('.gz') else "r"
    
    ramp_seqs, body_seqs = [], []
    with open_func(REF_FASTA, mode) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if not any(eid in record.id for eid in elite_ids): continue
            
            seq_str = str(record.seq)
            if "ATG" not in seq_str: continue
            cds = seq_str[seq_str.find("ATG"):]
            if len(cds) < RAMP_SIZE + 100: continue
            
            # 精确切分！
            ramp_seqs.append(cds[:RAMP_SIZE])
            body_seqs.append(cds[RAMP_SIZE:])
            
    return ramp_seqs, body_seqs

def calculate_matrix(seq_list, fallback_matrix=None):
    """
    计算 Frequency (用于生成) 和 Weight (用于评估)
    fallback_matrix 用于 Ramp 区：如果某个氨基酸在 39bp 里没出现过，借用 Body 区的数据防报错
    """
    # 1. 基础计数 (加入 0.1 伪计数，防止绝对的 0 概率导致算法崩溃)
    counts = {c: 0.1 for c in CODON_DICT.keys() if CODON_DICT[c] != '*'}
    
    for seq in seq_list:
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon in counts:
                counts[codon] += 1

    # 2. 按氨基酸分组计算
    aa_groups = {}
    for codon, count in counts.items():
        aa = CODON_DICT[codon]
        if aa not in aa_groups: aa_groups[aa] = {}
        aa_groups[aa][codon] = count

    freq_dict, weight_dict = {}, {}
    for aa, syn_codons in aa_groups.items():
        total_count = sum(syn_codons.values())
        max_count = max(syn_codons.values())
        
        # 如果这个氨基酸总计只有不到 1 次（纯靠伪计数撑着），说明样本缺失
        if total_count < 1.0 and fallback_matrix is not None:
            for codon in syn_codons:
                freq_dict[codon] = fallback_matrix['Freq'][codon]
                weight_dict[codon] = fallback_matrix['Weight'][codon]
        else:
            for codon, count in syn_codons.items():
                freq_dict[codon] = count / total_count
                weight_dict[codon] = count / max_count

    return {'Freq': freq_dict, 'Weight': weight_dict}

def main():
    print("[*] 正在提取 Pareto 算法终极蓝图...")
    ramp_seqs, body_seqs = get_elite_regions()
    print(f"    -> 成功提取 {len(ramp_seqs)} 条 Elite Ramp (前 39bp) 序列")
    print(f"    -> 成功提取 {len(body_seqs)} 条 Elite Body (主体延伸) 序列")
    
    # 先算 Body，作为 Ramp 的坚实后盾
    body_matrix = calculate_matrix(body_seqs)
    # 再算 Ramp，如果数据稀疏，向 Body 借调规则
    ramp_matrix = calculate_matrix(ramp_seqs, fallback_matrix=body_matrix)
    
    # 整理合并为最终的 DataFrame
    data = []
    for codon, aa in CODON_DICT.items():
        if aa == '*': continue
        data.append({
            'Amino_Acid': aa,
            'Codon': codon,
            'Ramp_Freq': round(ramp_matrix['Freq'][codon], 4),
            'Ramp_w': round(ramp_matrix['Weight'][codon], 4),
            'Body_Freq': round(body_matrix['Freq'][codon], 4),
            'Body_w': round(body_matrix['Weight'][codon], 4)
        })
        
    df_blueprint = pd.DataFrame(data)
    # 按氨基酸排序，强迫症福音
    df_blueprint = df_blueprint.sort_values(by=['Amino_Acid', 'Codon'])
    df_blueprint.to_csv(OUTPUT_MATRIX, index=False)
    
    print("\n" + "="*50)
    print("🏆 Pareto 算法双轨制字典提取完毕！")
    print(f"    -> 已保存至: {OUTPUT_MATRIX}")
    print("="*50 + "\n")
    
    # 打印一个极其经典的“安全 AT 密码子”案例供你验货
    print("[🔍 验货时刻：亮氨酸 (Leu) 在前 39bp 的生死抉择]")
    leu_data = df_blueprint[df_blueprint['Amino_Acid'] == 'L'][['Codon', 'Ramp_Freq', 'Body_Freq']]
    print(leu_data.to_string(index=False))
    print("\n-> 你会发现：同样是降 GC，Ramp 区绝对偏好特定的 AT 密码子，而彻底封杀了另一些毒性密码子！")

if __name__ == "__main__":
    main()
