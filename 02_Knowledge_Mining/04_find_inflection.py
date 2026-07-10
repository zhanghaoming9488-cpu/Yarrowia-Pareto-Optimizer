import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

FILE_PATH = "high_expression_stats.csv"

def find_log_knee_point(x, y):
    """
    在对数空间 (Log-space) 寻找几何拐点。
    这是单细胞测序和转录组学中最常用的阈值划分方法。
    """
    # 1. 对 TPM 进行 Log10 变换 (+1 防止 log(0))
    y_log = np.log10(y + 1)
    
    # 2. 将 X 和 Log(Y) 缩放到 0~1 之间
    x_norm = (x - np.min(x)) / (np.max(x) - np.min(x))
    y_norm = (y_log - np.min(y_log)) / (np.max(y_log) - np.min(y_log))
    
    # 3. 计算到对角线的距离
    # 直线方程 x + y = 1 (因为是从 0,1 到 1,0)
    # 距离正比于 x_norm + y_norm (找最小值点，即曲线向内凹得最深的点)
    distances = x_norm + y_norm
    return np.argmin(distances)

def main():
    print("[*] Loading RNA-seq data for Log-space analysis...")
    df = pd.read_csv(FILE_PATH)
    
    # 提取 TPM > 1 的基因并排序 (过滤掉完全不表达的噪声)
    df = df[df['mean_tpm'] > 1.0].sort_values(by='mean_tpm', ascending=False).reset_index(drop=True)
    df['rank'] = df.index + 1
    
    x = df['rank'].values
    y = df['mean_tpm'].values
    
    knee_idx = find_log_knee_point(x, y)
    knee_rank = df['rank'].iloc[knee_idx]
    knee_tpm = df['mean_tpm'].iloc[knee_idx]
    
    print("\n" + "="*50)
    print("🎯 Data-Driven Inflection Point (Log-Space) Found!")
    print(f"    -> Gene Rank (排名): {knee_rank}")
    print(f"    -> TPM Threshold (原始阈值): {knee_tpm:.2f}")
    print("="*50 + "\n")

    # 绘图展示
    plt.figure(figsize=(10, 6))
    
    # 注意这里 Y 轴用对数坐标
    plt.plot(x, y, color='#2980b9', linewidth=2.5, label='Sorted TPM (Log Scale)')
    plt.yscale('log') 
    
    plt.axvline(x=knee_rank, color='#e74c3c', linestyle='--', alpha=0.8)
    plt.axhline(y=knee_tpm, color='#e74c3c', linestyle='--', alpha=0.8)
    plt.scatter(knee_rank, knee_tpm, color='#e74c3c', s=100, zorder=5, 
                label=f'Log-Knee Point\n(Rank: {knee_rank}, TPM: {knee_tpm:.0f})')
    
    plt.title('RNA-seq TPM Distribution (Log Scale) & Optimal Cutoff', fontsize=14, fontweight='bold')
    plt.xlabel('Gene Rank', fontsize=12)
    plt.ylabel('Mean TPM (Log Scale)', fontsize=12)
    plt.legend(fontsize=11)
    plt.grid(True, which="both", ls="--", alpha=0.4)
    plt.tight_layout()
    
    plt.savefig('tpm_log_inflection.png', dpi=300)
    print("[*] Plot saved as 'tpm_log_inflection.png'.")

if __name__ == "__main__":
    main()
