import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

# === 配置 ===
GOLDEN_CSV = "golden_elite_genes.csv"
RNA_FASTA = "high_expression_genes.fasta" # 原始提取的序列库
# ===========

def get_sequences_from_golden_set():
    """提取 143 个黄金基因的完整序列"""
    df = pd.read_csv(GOLDEN_CSV)
    # 获取精英 RNA 的 ID 列表
    elite_ids = set(df['GeneID_RNA'].tolist())
    
    seqs = []
    for record in SeqIO.parse(RNA_FASTA, "fasta"):
        if record.id in elite_ids:
            # 找到起始密码子 ATG 并截取后面的 CDS
            seq_str = str(record.seq)
            if "ATG" in seq_str:
                cds = seq_str[seq_str.find("ATG"):]
                if len(cds) > 300: # 确保基因足够长
                    seqs.append(cds)
    return seqs

def calculate_at_content(seq_list, window_size):
    """
    计算特定窗口内的平均 AT 含量。
    在生信中，5'端局部极高的 AT 含量是降低 mRNA 二级结构 (减小 \Delta G)、
    促进核糖体结合的最核心标志。
    """
    at_contents = []
    for seq in seq_list:
        if len(seq) < window_size: continue
        segment = seq[:window_size]
        at_count = segment.count('A') + segment.count('T')
        at_contents.append(at_count / window_size)
    return np.mean(at_contents)

def main():
    print("[*] 启动 5' 端动态窗口扫描...")
    elite_seqs = get_sequences_from_golden_set()
    print(f"    -> 成功载入 {len(elite_seqs)} 条黄金序列 CDS。")
    
    # 我们测试从 15 bp 到 150 bp (步长为 3 bp，即 1 个密码子)
    window_sizes = np.arange(15, 153, 3)
    at_means = []
    
    print("    正在计算各窗口大小下的局部序列特征...")
    for w in window_sizes:
        at_mean = calculate_at_content(elite_seqs, w)
        at_means.append(at_mean)
        
    # 计算全局 AT 含量作为基线
    global_at = calculate_at_content(elite_seqs, 1000) # 取足够长代表全局
    
    # === 寻找拐点 (一阶导数法) ===
    # 当局部 AT 含量开始急剧向全局 AT 含量逼近并趋于平缓时，该区域即为 5' 调控边界
    gradients = np.gradient(at_means)
    
    # 绘图展示
    plt.figure(figsize=(10, 6))
    plt.plot(window_sizes, at_means, marker='o', color='#2c3e50', linewidth=2, label="Local AT Content (5' Ramp)")
    plt.axhline(y=global_at, color='#e74c3c', linestyle='--', label=f"Global AT Baseline ({global_at:.3f})")
    
    plt.title("Y. lipolytica Elite Genes: 5' Translational Ramp Analysis", fontsize=14, fontweight='bold')
    plt.xlabel("Distance from Start Codon (bp)", fontsize=12)
    plt.ylabel("Mean AT Content", fontsize=12)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    
    plt.savefig("5prime_window_analysis.png", dpi=300)
    print("\n[+] 扫描完成！分析图表已保存为 5prime_window_analysis.png")
    print("    -> 观察图表：局部特征在哪一个 bp 数值处开始平缓，那里就是真正的 5' 调控边界！")

if __name__ == "__main__":
    main()
