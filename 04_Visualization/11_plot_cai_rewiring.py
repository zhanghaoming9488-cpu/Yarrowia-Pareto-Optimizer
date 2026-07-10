import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
import gzip

# === 配置 ===
GOLDEN_CSV = "golden_elite_genes.csv"
REF_FASTA = "reference/Yarrowia_lipolytica.cdna.all.fa.gz"
RAMP_SIZE = 39 # 斩首行动：严格切掉前 13 个密码子
# ===========

# 标准密码子字典
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

def get_sequences():
    """解析全集，区分 Elite (去头) 和 Standard (全长)"""
    try:
        df_golden = pd.read_csv(GOLDEN_CSV)
        elite_ids = set(df_golden['GeneID_RNA'].astype(str).tolist())
    except FileNotFoundError:
        print(f"[!] 找不到 {GOLDEN_CSV}")
        return [], []

    open_func = gzip.open if REF_FASTA.endswith('.gz') else open
    mode = "rt" if REF_FASTA.endswith('.gz') else "r"
    
    bg_seqs = []
    elite_seqs = []
    
    with open_func(REF_FASTA, mode) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_str = str(record.seq)
            if "ATG" not in seq_str: continue
            cds = seq_str[seq_str.find("ATG"):]
            if len(cds) < 300: continue
            
            # 全基因组全部作为 Standard Background
            bg_seqs.append(cds)
            
            # Elite 基因切掉前 39 bp (13个密码子)
            if any(eid in record.id for eid in elite_ids):
                if len(cds) > RAMP_SIZE + 100:
                    elite_seqs.append(cds[RAMP_SIZE:])
                    
    return bg_seqs, elite_seqs

def count_and_calculate_weights(seq_list):
    """统计密码子频率并计算 CAI 权重 w"""
    counts = {c: 0 for c in CODON_DICT.keys()}
    for seq in seq_list:
        # 确保只读完整的密码子
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon in counts:
                counts[codon] += 1
                
    # 按氨基酸分组计算权重 w = count / max_count_for_this_AA
    aa_counts = {}
    for codon, count in counts.items():
        aa = CODON_DICT[codon]
        if aa == '*': continue
        if aa not in aa_counts: aa_counts[aa] = {}
        aa_counts[aa][codon] = count
        
    weights = {}
    for aa, syn_codons in aa_counts.items():
        max_count = max(syn_codons.values())
        for codon, count in syn_codons.items():
            weights[codon] = count / max_count if max_count > 0 else 0
            
    return weights, aa_counts

def main():
    print("[*] 正在加载序列并提取纯净 CDS...")
    bg_seqs, elite_seqs = get_sequences()
    print(f"    -> Standard 全量基因: {len(bg_seqs)} 条")
    print(f"    -> Elite 主体序列 (切除前 39bp): {len(elite_seqs)} 条")
    
    print("[*] 正在计算 CAI 权重 (w) 矩阵...")
    bg_w, bg_aa_counts = count_and_calculate_weights(bg_seqs)
    elite_w, _ = count_and_calculate_weights(elite_seqs)
    
    # 将字典转为 DataFrame 方便导出和绘图
    df_compare = pd.DataFrame({
        'Amino_Acid': [CODON_DICT[c] for c in bg_w.keys()],
        'Codon': list(bg_w.keys()),
        'Standard_w': list(bg_w.values()),
        'Elite_w': list(elite_w.values())
    })
    # 算一下偏好性漂移差值
    df_compare['Shift'] = df_compare['Elite_w'] - df_compare['Standard_w']
    df_compare.to_csv("cai_matrix_comparison.csv", index=False)
    print("    [+] 新旧 CAI 对照矩阵已保存为 cai_matrix_comparison.csv")
    
    # ================= 极具冲击力的可视化设计 =================
    print("[*] 正在生成双图联动直观对比...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={'width_ratios': [1, 1.2]})
    
    # --- 图 A: 散点对角线图 (Identity Plot) ---
    ax1.plot([-0.1, 1.1], [-0.1, 1.1], color='gray', linestyle='--', alpha=0.6, label='y=x (No Change)')
    
    # 用颜色标记 GC 含量 (高 GC 标红，低 GC 标蓝)
    for i, row in df_compare.iterrows():
        gc = (row['Codon'].count('G') + row['Codon'].count('C')) / 3
        color = '#e74c3c' if gc > 0.5 else ('#3498db' if gc < 0.5 else '#2ecc71')
        ax1.scatter(row['Standard_w'], row['Elite_w'], color=color, s=50, alpha=0.8, edgecolors='white')
        
        # 标记出偏好性发生巨大倒转的密码子 (变化率 > 0.4)
        if abs(row['Shift']) > 0.4:
            ax1.text(row['Standard_w'] + 0.02, row['Elite_w'] - 0.02, row['Codon'], fontsize=9)
            
    ax1.set_title("Genome vs. Elite Codon Preferences ($w$)", fontsize=13, fontweight='bold')
    ax1.set_xlabel("Standard $w$ (All Genes)", fontsize=11)
    ax1.set_ylabel("Elite $w$ (Top Genes, Body Only)", fontsize=11)
    ax1.grid(True, linestyle=':', alpha=0.5)
    
    # --- 图 B: 赢家与输家 (Top Shifted Amino Acids) ---
    # 找到偏好性发生最大重组的前 4 个氨基酸
    shift_sum = df_compare.groupby('Amino_Acid')['Shift'].apply(lambda x: x.abs().sum()).sort_values(ascending=False)
    # 排除 M 和 W (因为它们只有一个密码子)
    top_aas = [aa for aa in shift_sum.index if aa not in ['M', 'W']][:4]
    
    plot_data = df_compare[df_compare['Amino_Acid'].isin(top_aas)].copy()
    
    x_pos = np.arange(len(plot_data))
    width = 0.35
    
    ax2.bar(x_pos - width/2, plot_data['Standard_w'], width, label='Standard', color='#bdc3c7')
    ax2.bar(x_pos + width/2, plot_data['Elite_w'], width, label='Elite', color='#e67e22')
    
    ax2.set_xticks(x_pos)
    # 标签格式为 "氨基酸_密码子"
    ax2.set_xticklabels(plot_data['Amino_Acid'] + "\n" + plot_data['Codon'], fontsize=9)
    ax2.set_title("Codon Usage 'Re-wiring' in Elite Genes", fontsize=13, fontweight='bold')
    ax2.set_ylabel("CAI Weight ($w$)", fontsize=11)
    ax2.legend()
    ax2.grid(axis='y', linestyle='--', alpha=0.4)
    
    plt.tight_layout()
    plt.savefig("cai_comparison_dashboard.png", dpi=300)
    print("    [+] 顶级对比图表已保存为 cai_comparison_dashboard.png！")

if __name__ == "__main__":
    main()
