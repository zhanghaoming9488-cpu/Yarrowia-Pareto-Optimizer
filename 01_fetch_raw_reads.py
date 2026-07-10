import os
import subprocess
import pandas as pd
import sys

# === 你指定的补漏名单 ===
TARGET_IDS = ["SRR9021311", "SRR9021310"]

# === 配置 ===
OUTPUT_DIR = "raw_data"
ARIA2_PATH = "./aria2c.exe"  # 确保 aria2c.exe 在这里
METADATA_FILE = os.path.join(OUTPUT_DIR, "metadata.csv")
# ============

def download_with_aria2(url, dest_dir):
    cmd = [ARIA2_PATH, "-x", "16", "-s", "16", "-d", dest_dir, url]
    try:
        print(f"[*] 正在调用 aria2c 下载: {os.path.basename(url)}")
        subprocess.run(cmd, check=True)
    except Exception as e:
        print(f"[!] Aria2 下载出错: {e}")

def main():
    # 1. 检查 metadata 是否存在
    if not os.path.exists(METADATA_FILE):
        print(f"[!] 错误: 找不到 {METADATA_FILE}。请先运行一次之前的脚本以生成表格，或者手动下载 metadata。")
        return

    # 2. 读取表格
    try:
        df = pd.read_csv(METADATA_FILE)
    except Exception as e:
        print(f"[!] 读取 metadata.csv 失败: {e}")
        return

    print(f"[*] 正在从表格中检索目标文件: {TARGET_IDS} ...")
    
    found_count = 0
    
    # 3. 遍历目标ID
    for target_id in TARGET_IDS:
        # 在表格中查找 run_accession 等于 target_id 的行
        row = df[df['run_accession'] == target_id]
        
        if row.empty:
            print(f"[!] 警告: 在表格中未找到 ID {target_id}，跳过。")
            continue
            
        found_count += 1
        fastq_urls = row.iloc[0]['fastq_ftp']
        
        if pd.isna(fastq_urls):
            print(f"[!] 错误: {target_id} 没有下载链接。")
            continue

        # 处理链接 (ENA 可能有多个文件，用分号隔开)
        urls = str(fastq_urls).split(';')
        for u in urls:
            if not u.strip(): continue
            
            # 补全 http
            if not u.startswith('http'):
                download_url = 'http://' + u
            else:
                download_url = u
            
            filename = os.path.basename(download_url)
            local_path = os.path.join(OUTPUT_DIR, filename)
            
            # 检查文件是否已存在
            if os.path.exists(local_path) and os.path.getsize(local_path) > 0:
                print(f"    [√] {filename} 已存在，跳过。")
            else:
                download_with_aria2(download_url, OUTPUT_DIR)

    print("-" * 50)
    print(f"任务完成。共处理 {found_count}/{len(TARGET_IDS)} 个目标 ID。")

if __name__ == "__main__":
    main()