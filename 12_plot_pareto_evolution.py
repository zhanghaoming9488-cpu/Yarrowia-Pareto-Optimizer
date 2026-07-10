import pandas as pd
import numpy as np
import random
import math
import os
import matplotlib.pyplot as plt
import seaborn as plt_sns # 避免与 plt 冲突

# === 最终配置 (Final Configuration for mCherry) ===
BLUEPRINT_CSV = "pareto_codon_blueprint.csv"
TOXIC_CSV = "yali_toxic_pairs.csv"
OUTPUT_FILE = "Final_mCherry_Optimal_Sequences.csv"

# 真实的 mCherry 蛋白质序列 (236 AA)
AA_SEQ = "MASVEEIRNAQRAKGPATILAIGTATPDHCVYQSDYADFYFRVTKSEHMTALKKKFNRICDKSMIKKRYIHLTEEMLEEHPNIGAYMAPSLNIRQEIITAEVPKLGKEAALKALKEWGQPKSKITHLVFCTTSGVEMPGADYKLANLLGLEPSVRRVMLYHQGCYAGGTVLRTAKDLAENNAGARVLVVCSEITVVTFRGPSEDALDSLVGQALFGDGSAAVIVGSDPDISIERPLFQLVSAAQTFIPNSAGAIAGNLREVGLTFHLWPNVPTLISENIEKCLTQAFDPLGISDWNSLFWIAHPGGPAILDAVEAKLNLDKKKLEATRHVLSEYGNMSSACVLFILDEMRKKSLKGERATTGEGLDWGVLFGFGPGLTIETIVLHSIPMVTN"

POPULATION_SIZE = 2000   
GENERATIONS = 1000       
MUTATION_RATE = 0.005     
RAMP_SIZE = 13          
# ====================================

def load_knowledge_base():
    if not os.path.exists(BLUEPRINT_CSV) or not os.path.exists(TOXIC_CSV):
        print(f"[!] 找不到蓝图或黑名单 CSV 文件，请检查路径。")
        return None, None, None
        
    df_bp = pd.read_csv(BLUEPRINT_CSV)
    gen_dict, eval_dict = {}, {}
    for aa in df_bp['Amino_Acid'].unique():
        sub = df_bp[df_bp['Amino_Acid'] == aa]
        gen_dict[aa] = {
            'Ramp_Codons': sub['Codon'].tolist(), 'Ramp_Probs': sub['Ramp_Freq'].tolist(),
            'Body_Codons': sub['Codon'].tolist(), 'Body_Probs': sub['Body_Freq'].tolist()
        }
        for _, row in sub.iterrows():
            eval_dict[row['Codon']] = {'Ramp_w': row['Ramp_w'], 'Body_w': row['Body_w']}
            
    df_toxic = pd.read_csv(TOXIC_CSV)
    toxic_dict = {row['Codon_Pair']: abs(row['CPS']) for _, row in df_toxic.iterrows()}
    return gen_dict, eval_dict, toxic_dict

class DNAIndividual:
    def __init__(self, seq=None):
        self.seq = seq
        self.ramp_cai = 0.0  
        self.body_cai = 0.0  
        self.toxic = 0.0     
        self.rank = 0
        self.domination_count = 0
        self.dominated_solutions = []

def random_cds(aa_seq, gen_dict):
    codons = []
    for i, aa in enumerate(aa_seq):
        if i < RAMP_SIZE:
            c, p = gen_dict[aa]['Ramp_Codons'], gen_dict[aa]['Ramp_Probs']
        else:
            c, p = gen_dict[aa]['Body_Codons'], gen_dict[aa]['Body_Probs']
        p = np.array(p) / sum(p)
        codons.append(np.random.choice(c, p=p))
    return "".join(codons)

def calculate_objectives(seq, eval_dict, toxic_dict):
    ramp_w, body_w, toxic_penalty = [], [], 0.0
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if i // 3 < RAMP_SIZE:
            ramp_w.append(eval_dict.get(codon, {}).get('Ramp_w', 0.01))
        else:
            body_w.append(eval_dict.get(codon, {}).get('Body_w', 0.01))
            
    ramp_cai = math.exp(np.mean(np.log([w for w in ramp_w if w > 0]))) if ramp_w else 1.0
    body_cai = math.exp(np.mean(np.log([w for w in body_w if w > 0]))) if body_w else 1.0
    
    for i in range(0, len(seq) - 3, 3):
        pair = f"{seq[i:i+3]}-{seq[i+3:i+6]}"
        if pair in toxic_dict:
            toxic_penalty += toxic_dict[pair]
            
    return ramp_cai, body_cai, toxic_penalty

def mutate(cds, aa_seq, gen_dict):
    cds_list = [cds[i:i+3] for i in range(0, len(cds), 3)]
    idx = random.randint(0, len(cds_list)-1)
    aa = aa_seq[idx]
    
    if idx < RAMP_SIZE:
        codons, probs = gen_dict[aa]['Ramp_Codons'], gen_dict[aa]['Ramp_Probs']
    else:
        codons, probs = gen_dict[aa]['Body_Codons'], gen_dict[aa]['Body_Probs']
        
    probs = np.array(probs) / sum(probs)
    if len(codons) > 1:
        cds_list[idx] = np.random.choice(codons, p=probs)
    return "".join(cds_list)

def dominates(p, q):
    better_or_eq = (p.ramp_cai >= q.ramp_cai and p.body_cai >= q.body_cai and p.toxic <= q.toxic)
    strictly_better = (p.ramp_cai > q.ramp_cai or p.body_cai > q.body_cai or p.toxic < q.toxic)
    return better_or_eq and strictly_better

def assign_ranks_flat(population):
    for p in population:
        p.domination_count = 0
        p.dominated_solutions = []
        
    current_front = []
    for p in population:
        for q in population:
            if dominates(p, q):
                p.dominated_solutions.append(q)
            elif dominates(q, p):
                p.domination_count += 1
        if p.domination_count == 0:
            p.rank = 1
            current_front.append(p)
            
    rank = 1
    while current_front:
        next_front = []
        for p in current_front:
            for q in p.dominated_solutions:
                q.domination_count -= 1
                if q.domination_count == 0:
                    q.rank = rank + 1
                    next_front.append(q)
        rank += 1
        current_front = next_front

def main():
    gen_dict, eval_dict, toxic_dict = load_knowledge_base()
    if not gen_dict: return
    
    population = []
    for _ in range(POPULATION_SIZE):
        seq = random_cds(AA_SEQ, gen_dict)
        ind = DNAIndividual(seq)
        ind.ramp_cai, ind.body_cai, ind.toxic = calculate_objectives(seq, eval_dict, toxic_dict)
        population.append(ind)
        
    print(f"[*] 启动 mCherry 帕累托优化引擎 (Gen 0-{GENERATIONS})...")
    
    # === 用于绘图的数据记录 ===
    history = {'gen': [], 'max_body': [], 'max_ramp': [], 'min_toxic': []}
    snapshots = {}  # 记录第 1, 250, 500 代的所有个体状态
    snapshot_gens = [1, GENERATIONS // 20, GENERATIONS]
    
    for gen in range(1, GENERATIONS + 1):
        offspring = []
        for parent in population:
            new_seq = parent.seq
            if random.random() < MUTATION_RATE:
                new_seq = mutate(new_seq, AA_SEQ, gen_dict)
            child = DNAIndividual(new_seq)
            child.ramp_cai, child.body_cai, child.toxic = calculate_objectives(new_seq, eval_dict, toxic_dict)
            offspring.append(child)
            
        combined = population + offspring
        assign_ranks_flat(combined)
        combined.sort(key=lambda x: (x.rank, x.toxic, -x.body_cai))
        population = combined[:POPULATION_SIZE]
        
        # --- 记录数据 ---
        best_body = max(p.body_cai for p in population)
        best_ramp = max(p.ramp_cai for p in population)
        min_toxic = min(p.toxic for p in population)
        
        history['gen'].append(gen)
        history['max_body'].append(best_body)
        history['max_ramp'].append(best_ramp)
        history['min_toxic'].append(min_toxic)
        
        if gen in snapshot_gens:
            snapshots[gen] = {
                'ramp': [p.ramp_cai for p in population],
                'body': [p.body_cai for p in population],
                'toxic': [p.toxic for p in population]
            }
        
        if gen % 20 == 0:
            print(f"    Gen {gen}: Max Body={best_body:.3f}, Max Ramp={best_ramp:.3f}, Min Toxic={min_toxic:.4f}")

    # ================= 绘图代码部分 =================
    print("\n[*] 优化完成！正在生成科研绘图...")
    
    # 图 1：帕累托前沿散点图
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True, sharey=True)
    cmap = 'RdYlGn_r'
    
    for i, g in enumerate(snapshot_gens):
        sc = axes[i].scatter(snapshots[g]['ramp'], snapshots[g]['body'], 
                             c=snapshots[g]['toxic'], cmap=cmap, vmin=0, vmax=5, 
                             alpha=0.7, edgecolors='none', s=20)
        axes[i].set_title(f"Generation {g}")
        axes[i].set_xlabel("Ramp CAI")
        if i == 0: axes[i].set_ylabel("Body CAI")
        axes[i].grid(True, linestyle='--', alpha=0.5)

    cbar = fig.colorbar(sc, ax=axes, orientation='vertical', fraction=0.02, pad=0.02)
    cbar.set_label("Toxic Penalty (Red = High, Green = Low)")
    plt.suptitle("Pareto Front Evolution (mCherry)", fontsize=16)
    plt.savefig("RealData_Pareto_Evolution.png", dpi=300, bbox_inches='tight')
    
    # 图 2：多指标折线图
    fig2, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(history['gen'], history['max_body'], label="Max Body CAI", color='#1f77b4', lw=2)
    ax1.plot(history['gen'], history['max_ramp'], label="Max Ramp CAI", color='#ff7f0e', lw=2)
    ax1.set_xlabel("Generations")
    ax1.set_ylabel("CAI Score")
    ax1.grid(True, linestyle='--', alpha=0.5)
    
    ax2 = ax1.twinx()
    ax2.plot(history['gen'], history['min_toxic'], label="Min Toxic Penalty", color='#d62728', lw=2, linestyle='-.')
    ax2.set_ylabel("Toxic Penalty Score")
    
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='center right')
    plt.title("Multi-Objective Optimization Trends")
    plt.savefig("RealData_Pareto_Trends.png", dpi=300, bbox_inches='tight')
    
    print("[+] 绘图完毕！已在当前目录生成 png 文件。")

if __name__ == "__main__":
    main()