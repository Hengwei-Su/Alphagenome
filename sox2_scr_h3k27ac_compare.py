import os
import matplotlib.pyplot as plt
from alphagenome.data import genome
from alphagenome.visualization import plot_components

# 注意：这里改回 dna_client，而不是 dna_model
from alphagenome.models import dna_client

# ─── 1. 扔掉 Kaggle，使用云端 API 模式 ────────────────
# 填入你之前的 API KEY
os.environ["ALPHAGENOME_API_KEY"] = "AIzaSyA_DRr-O_t6d39cotLBTSvn2H7hKG-UeWU"
API_KEY = os.environ["ALPHAGENOME_API_KEY"]

# 创建客户端，后续所有的 1Mb 序列都会发送到云端计算，不消耗你本地的显卡
model = dna_client.create(API_KEY)

# ─── 2. 定义预测区间 (mm10 坐标系) ────────────────
# Sox2 新坐标起于 34649995。
# 1Mb 窗口 (1,048,576 bp) 居中计算：
# 新 Start = 34126939, 新 End = 35175515
interval = genome.Interval(
    chromosome='chr3', 
    start=34126939, 
    end=35175515
)

# ─── 3. 定义 Sox2 敲除变异 (mm10 序列) ────────────
# 你提供的 2466bp 真实 mm10 序列 (去除了换行符)
actual_ref_seq = (
    "GGATGGTTGTCTATTAACTTGTTCAAAAAAGTATCAGGAGTTGTCAAGGC"
    "AGAGAAGAGAGTGTTTGCAAAAAGGGAAAAGTACTTTGCTGCCTCTTTAA"
    "GACTAGGGCTGGGAGAAAGAAGAGGAGAGAGAAAGAAAGGAGAGAAGTTT"
    "GGAGCCCGAGGCTTAAGCCTTTCCAAAAACTAATCACAACAATCGCGGCG"
    "GCCCGAGGAGGAGAGCGCCTGTTTTTTCATCCCAATTGCACTTCGCCCGT"
    "CTCGAGCTCCGCTTCCCCCCAACTATTCTCCGCCAGATCTCCGCGCAGGG"
    "CCGTGCACGCCGAGGCCCCCGCCCGCGGCCCCTGCATCCCGGCCCCCGAG"
    "CGCGGCCCCCACAGTCCCGGCCGGGCCGAGGGTTGGCGGCCGCCGGCGGG"
    "CCGCGCCCGCCCAGCGCCCGCATGTATAACATGATGGAGACGGAGCTGAA"
    "GCCGCCGGGCCCGCAGCAAGCTTCGGGGGGCGGCGGCGGAGGAGGCAACG"
    "CCACGGCGGCGGCGACCGGCGGCAACCAGAAGAACAGCCCGGACCGCGTC"
    "AAGAGGCCCATGAACGCCTTCATGGTATGGTCCCGGGGGCAGCGGCGTAA"
    "GATGGCCCAGGAGAACCCCAAGATGCACAACTCGGAGATCAGCAAGCGCC"
    "TGGGCGCGGAGTGGAAACTTTTGTCCGAGACCGAGAAGCGGCCGTTCATC"
    "GACGAGGCCAAGCGGCTGCGCGCTCTGCACATGAAGGAGCACCCGGATTA"
    "TAAATACCGGCCGCGGCGGAAAACCAAGACGCTCATGAAGAAGGATAAGT"
    "ACACGCTTCCCGGAGGCTTGCTGGCCCCCGGCGGGAACAGCATGGCGAGC"
    "GGGGTTGGGGTGGGCGCCGGCCTGGGTGCGGGCGTGAACCAGCGCATGGA"
    "CAGCTACGCGCACATGAACGGCTGGAGCAACGGCAGCTACAGCATGATGC"
    "AGGAGCAGCTGGGCTACCCGCAGCACCCGGGCCTCAACGCTCACGGCGCG"
    "GCACAGATGCAACCGATGCACCGCTACGACGTCAGCGCCCTGCAGTACAA"
    "CTCCATGACCAGCTCGCAGACCTACATGAACGGCTCGCCCACCTACAGCA"
    "TGTCCTACTCGCAGCAGGGCACCCCCGGTATGGCGCTGGGCTCCATGGGC"
    "TCTGTGGTCAAGTCCGAGGCCAGCTCCAGCCCCCCCGTGGTTACCTCTTC"
    "CTCCCACTCCAGGGCGCCCTGCCAGGCCGGGGACCTCCGGGACATGATCA"
    "GCATGTACCTCCCCGGCGCCGAGGTGCCGGAGCCCGCTGCGCCCAGTAGA"
    "CTGCACATGGCCCAGCACTACCAGAGCGGCCCGGTGCCCGGCACGGCCAT"
    "TAACGGCACACTGCCCCTGTCGCACATGTGAGGGCTGGACTGCGAACTGG"
    "AGAAGGGGAGAGATTTTCAAAGAGATACAAGGGAATTGGGAGGGGTGCAA"
    "AAAGAGGAGAGTAGGAAAAATCTGATAATGCTCAAAAGGAAAAAAAATCT"
    "CCGCAGCGAAACGACAGCTGCGGAAAAAAACCACCAATCCCATCCAAATT"
    "AACGCAAAAACCGTGATGCCGACTAGAAAACTTTTATGAGAGATCTTGGG"
    "ACTTCTTTTTGGGGGACTATTTTTGTACAGAGAAAACCTGAGGGCGGCGG"
    "GGAGGGCGGGGGAATCGGACCATGTATAGATCTGGAGGAAAAAAACTACG"
    "CAAAACTTTTTTTTAAAGTTCTAGTGGTACGTTAGGCGCTTCGCAGGGAG"
    "TTCGCAAAAGTCTTTACCAGTAATATTTAGAGCTAGACTCCGGGCGATGA"
    "AAAAAAAGTTTTAATATTTGCAAGCAACTTTTGTACAGTATTTATCGAGA"
    "TAAACATGGCAATCAAATGTCCATTGTTTATAAGCTGAGAATTTGCCAAT"
    "ATTTTTCGAGGAAAGGGTTCTTGCTGGGTTTTGATTCTGCAGCTTAAATT"
    "TAGGACCGTTACAAACAAGGAAGGAGTTTATTCGGATTTGAACATTTTAG"
    "TTTTAAAATTGTACAAAAGGAAAACATGAGAGCAAGTACTGGCAAGACCG"
    "TTTTCGTGGTCTTGTTTAAGGCAAACGTTCTAGATTGTACTAAATTTTTA"
    "ACTTACTGTTAAAGGCAAAAAAAAAATGTCCATGCAGGTTGATATCGTTG"
    "GTAATTTATAATAGCTTTTGTTCAATCCTACCCTTTCATTTTGTTCACAT"
    "AAAAAATATGGAATTACTGTGTTTGAAATATTTTCTTATGGTTTGTAATA"
    "TTTCTGTAAATTGTGATATTTTAAGGTTTTTCCCCCCTTTTATTTTCCGT"
    "AGTTGTATTTTAAAAGATTCGGCTCTGTTATTGGAATCAGGCTGCCGAGA"
    "ATCCATGTATATATTTGAACTAATACCATCCTTATAACAGCTACATTTTC"
    "AACTTAAGTTTTTACTCCATTATGCACAGTTTGAGATAAATAAATTTTTG"
    "AAATATGGACACTGAA"
)

variant = genome.Variant(
    chromosome='chr3',
    position=34649995,               # 更新为 mm10 的起始坐标
    reference_bases=actual_ref_seq,  # 传入真实的参考序列
    alternate_bases='',              # 空字符串表示 Deletion
)

# --- 3. 运行多模态预测 ---
outputs = model.predict_variant(
    interval=interval,
    variant=variant,
    organism=dna_client.Organism.MUS_MUSCULUS,
    ontology_terms=[],  # 依然保持为空，拿回全量数据
    requested_outputs=[
        dna_client.OutputType.CHIP_HISTONE,
        dna_client.OutputType.RNA_SEQ  # <--- 新增 RNA 表达预测
    ]
)

# ─── 5. 分别过滤 H3K27ac 和 RNA-seq 的 Track ───────────────
# 5.1 过滤 H3K27ac (使用之前成功的 ES-Bruce4)
meta_histone = outputs.reference.chip_histone.metadata
mask_h3k27ac = (meta_histone['histone_mark'] == 'H3K27ac') & (meta_histone['biosample_name'] == 'ES-Bruce4')
wt_chip_histone = outputs.reference.chip_histone.filter_tracks(mask_h3k27ac.values)
ko_chip_histone = outputs.alternate.chip_histone.filter_tracks(mask_h3k27ac.values)

# 5.2 过滤 RNA-seq 
# RNA-seq 没有 histone_mark 列，所以只用细胞系过滤
meta_rna = outputs.reference.rna_seq.metadata
# 这里使用模糊匹配，因为 RNA-seq 的细胞系命名可能和 ChIP-seq 略有不同 (比如叫 ES-E14 或 Bruce4)
mask_rna = meta_rna.astype(str).apply(lambda x: x.str.contains('ES-Bruce4|ES-E14|ESC', case=False, na=False)).any(axis=1)
wt_rna = outputs.reference.rna_seq.filter_tracks(mask_rna.values)
ko_rna = outputs.alternate.rna_seq.filter_tracks(mask_rna.values)

# ─── 6. 多模态可视化 ──────────────────────────────
plot_components.plot(
    [
        # 第一张子图：H3K27ac
        plot_components.OverlaidTracks(
            tdata={
                'WT H3K27ac': wt_chip_histone,
                'Sox2-KO H3K27ac': ko_chip_histone,
            },
            colors={'WT H3K27ac': 'dimgrey', 'Sox2-KO H3K27ac': 'red'},
        ),
        # 第二张子图：RNA-seq
        plot_components.OverlaidTracks(
            tdata={
                'WT RNA-seq': wt_rna,
                'Sox2-KO RNA-seq': ko_rna,
            },
            colors={'WT RNA-seq': 'blue', 'Sox2-KO RNA-seq': 'orange'},
        ),
    ],
    interval=outputs.reference.chip_histone.interval.resize(200000),
    annotations=[plot_components.VariantAnnotation([variant], alpha=0.8)],
)

fig = plt.gcf()
fig.suptitle("mESC Multi-modal: Sox2 KO vs Distal SCR (mm10)", y=0.95)
plt.savefig("sox2_ko_multimodal_prediction.png", dpi=300, bbox_inches='tight')
print("多模态实验图表已生成！")