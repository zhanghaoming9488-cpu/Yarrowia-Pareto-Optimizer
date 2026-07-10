import os
import pandas as pd
import gzip
import re
from Bio import SeqIO
from collections import defaultdict

# === 配置 (扩大范围!) ===
QUANT_DIR = "quant_results"           # RNA-seq 结果目录
REF_GENOME = "reference/Yarrowia_lipolytica.cdna.all.fa.gz"
PROTEIN_FILE = "proteinGroups.txt"
TOP_N_SELECTION = 2000                # 扩大到前 2000 名
OUTPUT_CODON = "codon_stats_golden.csv" # 生成新的黄金标准
OUTPUT_RAMP = "ramp_stats_golden.csv"
# =======================

# --- 1. RNA 处理部分 ---
def load_rna_top_genes():
    print(f"[*] Loading RNA-seq data (Targeting Top {TOP_N_SELECTION})...")
    tpm_dict = {}
    for sample_dir in os.listdir(QUANT_DIR):
        f = os.path.join(QUANT_DIR, sample_dir, "abundance.tsv")
        if os.path.exists(f):
            df = pd.read_csv(f, sep='\t')
            tpm_dict[sample_dir] = dict(zip(df['target_id'], df['tpm']))
            
    tpm_df = pd.DataFrame(tpm_dict)
    tpm_df['mean_tpm'] = tpm_df.mean(axis=1)
    # 筛选 Top N
    top_genes = tpm_df.sort_values(by='mean_tpm', ascending=False).head(TOP_N_SELECTION)
    print(f"    Selected {len(top_genes)} high mRNA genes (Min TPM: {top_genes['mean_tpm'].min():.2f})")
    return top_genes.index.tolist()

def extract_rna_sequences(gene_list):
    print("[*] Extracting sequences & descriptions...")
    rna_data = []
    with gzip.open(REF_GENOME, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in gene_list:
                # 清洗描述，用于匹配
                desc = record.description.replace(record.id, "").lower()
                # 简单的清洗正则
                desc = re.sub(r"os=.*", "", desc) 
                desc = re.sub(r"[^a-z\s]", "", desc)
                rna_data.append({
                    "GeneID": record.id,
                    "Description": " ".join(desc.split()), # 去重空格
                    "Sequence": str(record.seq).upper()
                })
    return pd.DataFrame(rna_data)

# --- 2. Protein 处理部分 ---
def load_protein_top_genes():
    print(f"[*] Loading Proteomics data (Targeting Top {TOP_N_SELECTION})...")
    try:
        df = pd.read_csv(PROTEIN_FILE, sep='\t')
        # 找 Intensity 列
        if 'Intensity' in df.columns:
            df = df.sort_values(by='Intensity', ascending=False)
        elif 'iBAQ' in df.columns:
            df = df.sort_values(by='iBAQ', ascending=False)
        
        # 筛选 Top N
        df = df.head(TOP_N_SELECTION)
        
        # 提取描述
        prot_data = []
        for _, row in df.iterrows():
            header = str(row['Fasta headers']).lower()
            # 清洗
            desc = re.sub(r"os=.*", "", header)
            desc = re.sub(r"sp\|.*?\|", "", desc) # 去掉 sp|ID|
            desc = re.sub(r"tr\|.*?\|", "", desc)
            desc = re.sub(r"[^a-z\s]", "", desc)
            prot_data.append({
                "ProteinID": row.get('Majority protein IDs'),
                "Description": " ".join(desc.split())
            })
        print(f"    Selected {len(prot_data)} high abundance proteins.")
        return pd.DataFrame(prot_data)
    except Exception as e:
        print(f"[!] Error loading proteins: {e}")
        return pd.DataFrame()

# --- 3. 匹配与计算 ---
def calculate_rscu_and_ramp(sequences):
    # (这是之前的逻辑，整合到这里)
    print(f"[*] Calculating RSCU & Ramp for {len(sequences)} intersection genes...")
    
    # RSCU
    codon_counts = defaultdict(int)
    aa_counts = defaultdict(int)
    genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    
    # Ramp Stats
    head_gcs = []
    body_gcs = []

    valid_seqs = 0
    for seq in sequences:
        if len(seq) < 100 or len(seq) % 3 != 0: continue
        valid_seqs += 1
        
        # RSCU loop
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            if codon in genetic_code:
                codon_counts[codon] += 1
                aa_counts[genetic_code[codon]] += 1
        
        # Ramp loop (Head 40bp)
        head = seq[:40]
        body = seq[40:]
        head_gc = (head.count('G') + head.count('C')) / len(head)
        body_gc = (body.count('G') + body.count('C')) / len(body)
        head_gcs.append(head_gc)
        body_gcs.append(body_gc)

    # 导出 RSCU
    rscu_data = []
    for codon, aa in genetic_code.items():
        total_aa = aa_counts[aa]
        if total_aa == 0: continue
        
        # Count synonymous codons
        syn_codons = [k for k,v in genetic_code.items() if v == aa]
        expected = total_aa / len(syn_codons)
        
        observed = codon_counts[codon]
        rscu = observed / expected if expected > 0 else 0
        rscu_data.append({'AminoAcid': aa, 'Codon': codon, 'RSCU': round(rscu, 4)})
        
    pd.DataFrame(rscu_data).to_csv(OUTPUT_CODON, index=False)
    
    # 打印 Ramp 结果
    avg_head = sum(head_gcs)/len(head_gcs)
    avg_body = sum(body_gcs)/len(body_gcs)
    print(f"\n=== Golden Set Stats (N={valid_seqs}) ===")
    print(f"Avg Head GC: {avg_head:.2%}")
    print(f"Avg Body GC: {avg_body:.2%}")
    print(f"Difference:  {avg_body - avg_head:.2%}")
    if (avg_body - avg_head) > 0.02: # 大于2%就是显著
        print("[success] Strong ramp signal detected in Dual-High data!")

def main():
    # 1. Get Top 2000 RNA
    gene_ids = load_rna_top_genes()
    df_rna = extract_rna_sequences(gene_ids)
    
    # 2. Get Top 2000 Protein
    df_prot = load_protein_top_genes()
    
    if df_prot.empty: 
        print("Protein data failed. Stop.")
        return

    # 3. Match
    print("[*] Matching intersection...")
    intersection_seqs = []
    matched_count = 0
    
    # 使用 Set 加速查找
    # 将蛋白描述拆成关键词集合，比如 "elongation", "factor", "1-alpha"
    prot_keywords = []
    for idx, row in df_prot.iterrows():
        # 只有长度大于3的词才算关键词
        words = set([w for w in row['Description'].split() if len(w) > 3])
        prot_keywords.append(words)

    for idx, row in df_rna.iterrows():
        rna_desc = row['Description']
        rna_words = set([w for w in rna_desc.split() if len(w) > 3])
        if not rna_words: continue
        
        # 只要有 2 个以上的关键词重叠，就算匹配 (比较宽松，但因为都在 Top 列表里，风险可控)
        is_match = False
        for prot_words in prot_keywords:
            common = rna_words & prot_words
            if len(common) >= 2: # 至少重叠2个关键词
                is_match = True
                break
        
        if is_match:
            intersection_seqs.append(row['Sequence'])
            matched_count += 1
            
    print(f"    Found {matched_count} genes in the intersection (Dual High).")
    
    if matched_count < 50:
        print("[!] Warning: Still low. But let's see the stats.")
    else:
        print("[*] Great sample size!")

    # 4. Calculate Stats
    if intersection_seqs:
        calculate_rscu_and_ramp(intersection_seqs)
        print(f"\n[Done] New golden rules saved to {OUTPUT_CODON}")

if __name__ == "__main__":
    main()
