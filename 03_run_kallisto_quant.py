import os
import subprocess
from tqdm import tqdm

# === 配置 ===
KALLISTO_PATH = "./kallisto.exe"
INDEX_PATH = "Yarrowia.idx"
REF_PATH = "reference/Yarrowia_lipolytica.cdna.all.fa.gz"
RAW_DATA_DIR = "raw_data"
OUTPUT_DIR = "quant_results"
# ===========

def run_quantification():
    # 1. 检查索引是否存在
    if not os.path.exists(INDEX_PATH):
        print("[!] Index not found. Please rebuild index first!")
        return

    # 2. 确保输出主目录存在 (关键修复！)
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print(f"[*] Created output directory: {OUTPUT_DIR}")

    # 3. 扫描文件
    files = [f for f in os.listdir(RAW_DATA_DIR) if f.endswith(".gz")]
    samples = {}
    for f in files:
        sample_id = f.split("_")[0] if "_" in f else f.split(".")[0]
        if sample_id not in samples: samples[sample_id] = []
        samples[sample_id].append(os.path.join(RAW_DATA_DIR, f))
    
    print(f"[*] Found {len(samples)} samples. Starting quantification...")

    # 4. 运行
    for sample_id, input_files in samples.items():
        # 构造输出路径
        output_subdir = os.path.join(OUTPUT_DIR, sample_id)
        
        # === 关键修复：Windows 路径兼容 ===
        # 把反斜杠 \ 强制换成正斜杠 /，防止 Kallisto 报错
        output_subdir = output_subdir.replace("\\", "/")
        
        print(f"\n--- Processing {sample_id} ---")
        
        # 构造命令
        cmd = [KALLISTO_PATH, "quant", "-i", INDEX_PATH, "-o", output_subdir, "-b", "0"]
        
        # 处理输入文件路径 (也换成正斜杠)
        input_files = [f.replace("\\", "/") for f in input_files]
        
        if len(input_files) == 1:
            cmd.extend(["--single", "-l", "200", "-s", "20", input_files[0]])
        else:
            input_files.sort()
            cmd.extend(input_files)
            
        try:
            # 这里的 stderr=None 让报错显示出来，以便调试
            subprocess.run(cmd, check=True)
            print(f"    [OK] Done.")
        except subprocess.CalledProcessError as e:
            print(f"    [FAIL] Kallisto crashed for {sample_id}")

if __name__ == "__main__":
    run_quantification()