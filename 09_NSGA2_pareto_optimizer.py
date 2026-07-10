import pandas as pd
import numpy as np
import random
import math
import os

# === 最终配置 (Final Configuration) ===
BLUEPRINT_CSV = "pareto_codon_blueprint.csv"
TOXIC_CSV = "yali_toxic_pairs.csv"
OUTPUT_FILE = "Final_Yali_Optimal_Sequences.csv"

AA_SEQ = "MKLSTILFTACATLALALAGRRRSVQWCAVSQPEATKCFQWQRNMRKVRGPPVSCIKRDSPIQCIQAIAENRADAVTLDGGFIYEAGLAPYKLRPVAAEVYGTERQPRTHYYAVAVVKKGGSFQLNELQGLKSCHTGLRRTAGWNVPIGTLRPFLNWTGPPEPIEAAVARFFSASCVPGADKGQFPNLCRLCAGTGENKCAFSSQEPYFSYSGAFKCLRDGAGDVAFIRESTVFEDLSDEAERDEYELLCPDNTRKPVDKFKDCHLARVPSHAVVARSVNGKEDAIWNLLRQAQEKFGKDKSPKFQLFGSPSGQKDLLFKDSAIGFSRVPPRIDSGLYLGSGYFTAIQNLRKSEEEVAARRARVVWCAVGEQELRKCNQWSGLSEGSVTCSSASTTEDCIALVLKGEADAMSLDGGYVYTAGKCGLVPVLAENYKSQQSSDPDPNCVDRPVEGYLAVAVVRRSDTSLTWNSVKGKKSCHTAVDRTAGWNIPMGLLFNQTGSCKFDEYFSQSCAPGSDPRSNLCALCIGDEQGENKCVPNSNERYYGYTGAFRCLAENAGDVAFVKDVTVLQNTDGNNNEAWAKDLKLADFALLCLDGKRKPVTEARSCHLAMAPNHAVVSRMDKVERLKQVLLHQQAKFGRNGSDCPDKFCLFQSETKNLLFNDNTECLARLHGKTTYEKYLGPQYVAGITNLKKCSTSPLLEACEFLRK"

POPULATION_SIZE = 3000   
GENERATIONS = 2000       
MUTATION_RATE = 0.002     
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
        # 扁平化参数
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
    """三维空间支配判断"""
    better_or_eq = (p.ramp_cai >= q.ramp_cai and p.body_cai >= q.body_cai and p.toxic <= q.toxic)
    strictly_better = (p.ramp_cai > q.ramp_cai or p.body_cai > q.body_cai or p.toxic < q.toxic)
    return better_or_eq and strictly_better

def assign_ranks_flat(population):
    """
    终极扁平化排序法：不再生成嵌套列表，直接修改每个个体的 self.rank！
    彻底杜绝了 len(front) 的崩溃 bug。
    """
    for p in population:
        p.domination_count = 0
        p.dominated_solutions = []
        
    current_front = []
    # 找出所有第一梯队 (Rank 1)
    for p in population:
        for q in population:
            if dominates(p, q):
                p.dominated_solutions.append(q)
            elif dominates(q, p):
                p.domination_count += 1
        if p.domination_count == 0:
            p.rank = 1
            current_front.append(p)
            
    # 逐层推导后续梯队的 Rank
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
        
    print(f"[*] 启动终极扁平架构达尔文引擎 (Gen 0-{GENERATIONS})...")
    
    for gen in range(GENERATIONS):
        offspring = []
        for parent in population:
            new_seq = parent.seq
            if random.random() < MUTATION_RATE:
                new_seq = mutate(new_seq, AA_SEQ, gen_dict)
            child = DNAIndividual(new_seq)
            child.ramp_cai, child.body_cai, child.toxic = calculate_objectives(new_seq, eval_dict, toxic_dict)
            offspring.append(child)
            
        combined = population + offspring
        
        # 1. 盖戳：给所有人打上 Rank 排名
        assign_ranks_flat(combined)
        
        # 2. 暴力大排序：优先排 Rank (小越好)，然后比毒性 (小越好)，最后比极速 (大越好)
        combined.sort(key=lambda x: (x.rank, x.toxic, -x.body_cai))
        
        # 3. 粗暴截断：直接砍掉后 500 名，保留前 500 名！
        population = combined[:POPULATION_SIZE]
        
        if (gen + 1) % 30 == 0:
            best_body = max(p.body_cai for p in population)
            best_ramp = max(p.ramp_cai for p in population)
            min_toxic = min(p.toxic for p in population)
            print(f"    Gen {gen+1}: Max Body CAI={best_body:.3f}, Max Ramp CAI={best_ramp:.3f}, Min Toxic={min_toxic:.4f}")

    print("\n[Done] 正在执行“三位一体”科研对比筛选 (抓取 A/B/C 种子)...")
    
    # 1. 基础过滤：首先必须是零毒性 (Toxic Penalty = 0) 的精英
    safe_pop = [p for p in population if p.toxic == 0]
    
    # 防御性编程：如果没人能做到绝对零毒性，则选毒性最低的一批人
    if not safe_pop:
        min_toxic = min(p.toxic for p in population)
        safe_pop = [p for p in population if p.toxic <= min_toxic]

    # 2. 定向抓取三名性格鲜明的选手
    # 选手 A: 极速先锋 (Max Body CAI)
    design_A = max(safe_pop, key=lambda x: x.body_cai)
    
    # 选手 B: 起步之王 (Max Ramp CAI)
    design_B = max(safe_pop, key=lambda x: x.ramp_cai)
    
    # 选手 C: 均衡大师 (Balanced - 离理想点 [1.0, 1.0] 最近)
    def get_dist(p):
        return math.sqrt((1.0 - p.ramp_cai)**2 + (1.0 - p.body_cai)**2)
    design_C = min(safe_pop, key=get_dist)

    # 3. 汇总去重并导出
    final_candidates = [
        {'ID': 'Design_A_MaxBody', 'Desc': 'Priority: Translation Speed', 'Obj': design_A},
        {'ID': 'Design_B_MaxRamp', 'Desc': 'Priority: Start Efficiency', 'Obj': design_B},
        {'ID': 'Design_C_Balanced', 'Desc': 'Priority: Multi-objective Balance', 'Obj': design_C}
    ]

    output = []
    seen_seqs = set()
    for item in final_candidates:
        ind = item['Obj']
        if ind.seq not in seen_seqs:
            output.append({
                'Design_ID': item['ID'],
                'Strategy': item['Desc'],
                'Ramp_CAI': ind.ramp_cai,
                'Body_CAI': ind.body_cai,
                'Toxic_Penalty': ind.toxic,
                'Sequence': ind.seq
            })
            seen_seqs.add(ind.seq)

    # 4. 打印对照报告
    print("-" * 70)
    for res in output:
        print(f"[{res['Design_ID']}] -> {res['Strategy']}")
        print(f"   指标: Body_CAI={res['Body_CAI']:.3f} | Ramp_CAI={res['Ramp_CAI']:.3f} | Toxic={res['Toxic_Penalty']}")
        print(f"   Seq (前30bp): {res['Sequence'][:30]}...")
    print("-" * 70)

    pd.DataFrame(output).to_csv(OUTPUT_FILE, index=False)
    print(f"\n[+] 科研对照组已生成！请查看: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()