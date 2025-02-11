import pandas as pd
import os

# 定义标志物列表和结局类型列表
markers = ["miRNA", "mRNA", "lncRNA", "protein", "CNV", "methylation", "APA", "AS", "SNP"]
types = ["OS", "DSS", "DFI", "PFI"]


# 定义输入和输出的根目录
base_input_dir = "/home/wuzj/survival/03.result"
base_output_dir = "/home/wuzj/survival/03.result"

# 设置分块读取大小
chunksize = 100000  # 每次读取10万行

# 遍历每种标志物和结局类型
for marker in markers:
    print(f"{marker} start.")
    for outcome_type in types:
        # 构建输入文件路径
        input_file = os.path.join(base_input_dir, f"res_{outcome_type}", marker, f"{marker}_all_survival_results.csv")
        output_file = os.path.join(base_output_dir, f"res_{outcome_type}", marker, f"filtered_{marker}_all_survival_results.csv")

        # 检查输入文件是否存在
        if not os.path.exists(input_file):
            print(f"输入文件 {input_file} 不存在，跳过处理...")
            continue

        print(f"{marker} {outcome_type} start.")

        # 初始化表头保存标志
        header_saved = False

        # 循环逐块读取文件
        for chunk in pd.read_csv(input_file, sep='\t', chunksize=chunksize, low_memory=False):
            # 将列强制转换为数值类型，无法转换的值变为 NaN
            for col in ['KMp', 'UniCoxp', 'MultiCoxp']:
                chunk[col] = pd.to_numeric(chunk[col], errors='coerce')

            # 筛选条件：KMp、UniCoxp 和 MultiCoxp 均小于0.05
            filtered_chunk = chunk[
                (chunk['KMp'] < 0.05) & (chunk['UniCoxp'] < 0.05) & (chunk['MultiCoxp'] < 0.05)
            ]

            # 将筛选结果追加写入文件
            if not header_saved:
                # 第一次写入包括表头
                filtered_chunk.to_csv(output_file, sep='\t', mode='w', index=False)
                header_saved = True
            else:
                # 之后不再写入表头
                filtered_chunk.to_csv(output_file, sep='\t', mode='a', index=False, header=False)

        print(f"{marker} {outcome_type} done.")
        
    print(f"{marker} done.")

print("all marker done.")
