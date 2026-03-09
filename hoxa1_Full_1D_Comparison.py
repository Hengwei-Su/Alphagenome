from alphagenome.data import genome
from alphagenome.models import dna_client
import matplotlib.pyplot as plt
import numpy as np
import os

# 1. 建立连接
API_KEY = os.environ.get("ALPHAGENOME_API_KEY", "AIzaSyB4I97Pl77lK1TO-0PNVnJ6q8Jo0w-W1Hg")
model = dna_client.create(API_KEY)

# 2. 定义 Hoxa1 位点 (1MB 窗口)
variant = genome.Variant(chromosome="chr6", position=52155000, reference_bases="C", alternate_bases="T")
interval = genome.Interval(chromosome="chr6", start=52155000-524288, end=52155000+524288)

# 3. 使用你刚才探测出的精确名称
target_modalities = [
    dna_client.OutputType.RNA_SEQ,
    dna_client.OutputType.ATAC,
    dna_client.OutputType.DNASE,
    dna_client.OutputType.CAGE,
    dna_client.OutputType.PROCAP,
    dna_client.OutputType.CHIP_HISTONE,
    dna_client.OutputType.CHIP_TF,
    dna_client.OutputType.SPLICE_SITES,
    dna_client.OutputType.SPLICE_JUNCTIONS,
    dna_client.OutputType.SPLICE_SITE_USAGE,
    # 注意：CONTACT_MAPS 返回的是二维矩阵，暂时不放进这个一维绘图清单里
]

print("🚀 正在请求全模态 10 种一维数据预测 (由于数据量大，请耐心等待)...")

pred = model.predict_variant(
    interval=interval,
    variant=variant,
    requested_outputs=target_modalities,
    ontology_terms=["CL:0000127"], # 使用最通用的肺组织
)

# 4. 绘图：将所有返回的模态堆叠显示
n = len(target_modalities)
fig, axes = plt.subplots(n, 1, figsize=(15, 2.5 * n), sharex=True)

for i, mod_type in enumerate(target_modalities):
    ax = axes[i]
    # 根据枚举值获取数据
    attr_name = mod_type.name.lower()
    ref_data = getattr(pred.reference, attr_name).values.squeeze()
    alt_data = getattr(pred.alternate, attr_name).values.squeeze()
    
    # 绘图
    ax.plot(ref_data, label="WT (REF)", color='royalblue', alpha=0.6)
    ax.plot(alt_data, label="Mutant (ALT)", color='crimson', linestyle='--', alpha=0.6)
    ax.set_ylabel(attr_name.upper())
    if i == 0: ax.legend()

plt.suptitle("AlphaGenome Full 1D Modalities: Hoxa1 Case Study", fontsize=16)
plt.tight_layout()
plt.savefig("Hoxa1_Full_1D_Comparison.png", dpi=300)
print("✅ 成功生成 10 种模态对比图：Hoxa1_Full_1D_Comparison.png")