import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
import gzip

# === 配置 ===
GOLDEN_CSV = "golden_elite_genes.csv"
# 这里务必填入你包含全基因组 6000+ 个基因的 Fasta 文件路径！
REF_FASTA = "reference/Yarrowia_lipolytica.cdna.all.fa.gz" 
# ===========

def get_cds(seq_str):
    """寻找最靠近 5' 端的 ATG 并截取下游序列"""
    if "ATG" in seq_str:
        cds = seq_str[seq_str.find("ATG"):]
        if len(cds) > 300: # 只取有效长度的基因
            return cds
    return None

def calculate_at_content(seq_list, window_size):
    """计算特定窗口内的平均 AT 含量"""
    at_contents = []
    for seq in seq_list:
        if len(seq) < window_size: continue
        segment = seq[:window_size]
        at_count = segment.count('A') + segment.count('T')
        at_contents.append(at_count / window_size)
    return np.mean(at_contents) if at_contents else 0

def main():
    print("[*] 启动顶刊级对照扫描：Elite vs. All Genes...")
    
    # 1. 读取黄金基因的 ID
    try:
        df_golden = pd.read_csv(GOLDEN_CSV)
        # 把 ID 转为字符串集合，方便快速匹配
        elite_ids = set(df_golden['GeneID_RNA'].astype(str).tolist())
    except FileNotFoundError:
        print(f"[!] 找不到 {GOLDEN_CSV}")
        return

    # 2. 从全集 Fasta 中分离出 Elite 和 Background
    elite_seqs = []
    all_seqs = []
    
    print(f"    正在解析全基因组库 ({REF_FASTA})，这可能需要几秒钟...")
    
    # 兼容 .gz 压缩包或普通 .fasta 文件
    open_func = gzip.open if REF_FASTA.endswith('.gz') else open
    mode = "rt" if REF_FASTA.endswith('.gz') else "r"
    
    try:
        with open_func(REF_FASTA, mode) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                cds = get_cds(str(record.seq))
                if not cds: continue
                
                # 记录所有有效的基因作为 Background
                all_seqs.append(cds)
                
                # 判断是否是 Elite 基因 (模糊匹配，因为 FASTA 表头可能很长)
                is_elite = any(eid in record.id for eid in elite_ids)
                if is_elite:
                    elite_seqs.append(cds)
    except FileNotFoundError:
        print(f"[!] 找不到基因组文件 {REF_FASTA}，请检查路径。")
        return

    print(f"    -> 成功提取 Elite 精英序列: {len(elite_seqs)} 条")
    print(f"    -> 成功提取 全局背景序列 (All Genes): {len(all_seqs)} 条")

    # 3. 滑动窗口计算 (从 15bp 到 150bp，步长 3bp)
    window_sizes = np.arange(15, 153, 3)
    elite_at_means = []
    all_at_means = []
    
    print("    正在计算各窗口大小下的局部 AT 含量特征...")
    for w in window_sizes:
        elite_at_means.append(calculate_at_content(elite_seqs, w))
        all_at_means.append(calculate_at_content(all_seqs, w))
        
    # 计算全局基准线 (取全基因组前 1000bp 的平均值)
    global_at_baseline = calculate_at_content(all_seqs, 1000)

    # 4. 绘制神仙打架的对比图
    plt.figure(figsize=(10, 6))
    
    # 画精英基因 (醒目红色)
    plt.plot(window_sizes, elite_at_means, marker='o', color='#e74c3c', linewidth=2.5, 
             label="Elite Genes (Top 143)")
             
    # 画全量基因 (低调蓝色虚线)
    plt.plot(window_sizes, all_at_means, marker='s', color='#3498db', linewidth=2, linestyle='--',
             label=f"All Genes Background (n={len(all_seqs)})")
             
    # 画全局基准线
    plt.axhline(y=global_at_baseline, color='#7f8c8d', linestyle=':', linewidth=1.5,
                label=f"Genome-wide AT Baseline ({global_at_baseline:.3f})")
    
    plt.title("Translational Ramp Comparison: Elite vs. All Genes", fontsize=15, fontweight='bold')
    plt.xlabel("Distance from Start Codon (bp)", fontsize=13)
    plt.ylabel("Mean AT Content", fontsize=13)
    plt.legend(fontsize=11, loc='upper right')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    
    plot_filename = "5prime_ramp_comparison.png"
    plt.savefig(plot_filename, dpi=300)
    print(f"\n[+] 完美！高分文章级对比图表已保存为 {plot_filename}")

if __name__ == "__main__":
    main()