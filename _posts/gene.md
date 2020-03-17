


# TSS处理

识别TSS就是要找到5'端的位置。目前公开5'端数据集包括

已有的<kbd>参考基因注释数据集</kbd>包括[NCBI基因](https://www.ncbi.nlm.nih.gov/genome/guide/human/)，[GENCODE](https://www.gencodegenes.org/)和 [ENSEMBL基因](https://asia.ensembl.org/Homo_sapiens/Info/Index)。通过在各种实验条件下将**DNA甲基化区域定位到参考基因组**（如[GRCh39](https://www.ncbi.nlm.nih.gov/grc)）可以分析特定启动子在哪些组织和阶段被甲基化。 通过基于参考基因集中的基因注释对**表达水平进行定量分析**可以轻松获得分配给这些参考基因的其他信息（例如相关疾病），这些可以看成是<kbd>精准医学</kbd>的背景。利用通过测序获得的大规模数据集的研究领域之一是转录调控领域。 **RNA-seq**可以测量转录本和基因的表达水平。 **CAGE（Cap Analysis of Gene Expression, 基因表达的上限（帽）分析）**可以识别哪些启动子和增强子区域是活跃的，并且可以从转录开始的位置定义确切的基因组位置[^5]。



# 计算TSS的Enrichment score


# 参考文献:

[^1]: [转录起始位点是启动子吗？](https://www.zhihu.com/question/290262821)
[^2]: [Kapranov, Philipp. "From transcription start site to cell biology." Genome biology 10.4 (2009): 217.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2688922/)
[^3]: [Promoters](https://www.addgene.org/mol-bio-reference/promoters/)
[^4]: [Valen, Eivind, and Albin Sandelin. "Genomic and chromatin signals underlying transcription start-site selection." Trends in genetics 27.11 (2011): 475-485.](https://www.sciencedirect.com/science/article/abs/pii/S016895251100134X)
[^5]: [Abugessaisa, Imad, et al. "refTSS: a reference data set for human and mouse transcription start sites." Journal of molecular biology 431.13 (2019): 2407-2422.](https://www.sciencedirect.com/science/article/pii/S0022283619302530#bb0010)






