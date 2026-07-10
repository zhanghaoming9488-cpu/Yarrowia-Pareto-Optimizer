import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import gzip
from Bio import SeqIO
from collections import defaultdict

# === 配置 ===
REF_FASTA = "reference/Yarrowia_lipolytica.cdna.all.fa.gz" # 确保路径正确
OUTPUT_HEATMAP = "Yali_CPS_Heatmap.png"
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

def get_valid_cds_from_genome():
    open_func = gzip.open if REF_FASTA.endswith('.gz') else open
    mode = "rt" if REF_FASTA.endswith('.gz') else "r"
    seqs = []
    with open_func(REF_FASTA, mode) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_str = str(record.seq)
            if "ATG" not in seq_str: continue
            cds = seq_str[seq_str.find("ATG"):]
            cds = cds[:len(cds) - (len(cds) % 3)]
            if len(cds) > 300:
                seqs.append(cds)
    return seqs

def mine_all_cps(seq_list):
    """计算并返回全量 3721 种密码子对的 CPS"""
    codon_count = defaultdict(float)
    aa_count = defaultdict(float)
    codon_pair_count = defaultdict(float)
    aa_pair_count = defaultdict(float)
    
    for seq in seq_list:
        codons = [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3]
        for i in range(len(codons) - 1):
            c1, c2 = codons[i], codons[i+1]
            if c1 not in CODON_DICT or c2 not in CODON_DICT: continue
            a1, a2 = CODON_DICT[c1], CODON_DICT[c2]
            if a1 == '*' or a2 == '*': continue
            
            codon_count[c1] += 1
            codon_count[c2] += 1 
            aa_count[a1] += 1
            aa_count[a2] += 1
            codon_pair_count[f"{c1}-{c2}"] += 1
            aa_pair_count[f"{a1}-{a2}"] += 1

    results = []
    pseudo = 1.0 
    
    for c1 in CODON_DICT:
        for c2 in CODON_DICT:
            a1, a2 = CODON_DICT[c1], CODON_DICT[c2]
            if a1 == '*' or a2 == '*': continue
            
            pair = f"{c1}-{c2}"
            aa_pair = f"{a1}-{a2}"
            O_AB = codon_pair_count[pair]
            
            # 使用 pseudo 防止没出现过的组合报错
            F_A = codon_count[c1] + pseudo
            F_B = codon_count[c2] + pseudo
            F_AA_A = aa_count[a1] + pseudo
            F_AA_B = aa_count[a2] + pseudo
            F_AA_pair = aa_pair_count[aa_pair] + pseudo
            
            expected = ((F_A * F_B) / (F_AA_A * F_AA_B)) * F_AA_pair
            cps = math.log((O_AB + pseudo) / (expected + pseudo))
            
            results.append({
                'Codon_1': c1,
                'Codon_2': c2,
                'CPS': round(cps, 4)
            })
    return pd.DataFrame(results)

def main():
    print("[*] 正在加载全基因组序列并计算全量 CPS...")
    seqs = get_valid_cds_from_genome()
    df_all_cps = mine_all_cps(seqs)
    
    print("[*] 正在绘制 61x61 密码子对适应度热图...")
    # 将扁平的数据表转换成 61x61 的二维矩阵
    cps_matrix = df_all_cps.pivot(index='Codon_1', columns='Codon_2', values='CPS')

    # 设置画图参数
    plt.figure(figsize=(16, 14))
    
    # 核心绘图：RdBu 为红蓝渐变。center=0 表示 0 分处是白色。
    # cmap 解释：负分(有毒)会变红，正分(偏好)会变蓝
    sns.heatmap(cps_matrix, cmap='RdBu', center=0, 
                cbar_kws={'label': 'Codon Pair Score (CPS)'})
                
    plt.title('Global Codon Pair Score (CPS) Landscape of Yarrowia lipolytica\n(Dark Red indicates Toxic Pairs)', fontsize=20, pad=20)
    plt.xlabel('Second Codon (Position i+1)', fontsize=16)
    plt.ylabel('First Codon (Position i)', fontsize=16)
    
    plt.xticks(fontsize=8, rotation=90)
    plt.yticks(fontsize=8, rotation=0)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_HEATMAP, dpi=300)
    print(f"\n[+] 完美！高分辨率热图已保存至: {OUTPUT_HEATMAP}")

if __name__ == "__main__":
    main()
