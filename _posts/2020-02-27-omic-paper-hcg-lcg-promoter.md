---
title: "组学之转录起始位点综述论文之：转录起始位点选择的基因组和染色质信号"
date: 2020-02-27 12:14:00 +0800
categories: [Genomics, Paper]
tags: [Transcription]
---

参考文献 [Valen Eivind, and Albin Sandelin. "Genomic and chromatin signals underlying transcription start-site selection." Trends in genetics 27.11 (2011): 475-485.](https://www.sciencedirect.com/science/article/abs/pii/S016895251100134X)


通过对全基因组测序数据分析不仅可以确定TSS的位置，而且可以确定蛋白质与DNA的相互作用以及整个基因组中的染色质特征， 从而可以构建基础生物学的定量和预测模型。

在任何给定时间, 特定细胞类型或组织中仅表达有限数量的基因, 这些以多种方式表达的一组基因定义了该特定细胞或组织的功能。

基因转录调控着基因表达，该调控是从TSS[^1]处的转录开始的。在TSS位点处，RNA聚合酶II（RNAPII）[^2]及其辅因子共同形成起始前复合物（Pre-initiation complex, PIC[^3]）并与核心启动子[^4]绑定后开始转录下游基因。

影响转录过程的因素：
- 通过转录因子蛋白(TF[^5])与DNA中特定位点的结合（TF绑定位点，TFBS[^6]）来稳定或抑制转录过程。
- DNA本身的化学修饰(modification)也可以影响转录。
- DNA包裹在组蛋白周围，于是组蛋白包裹DNA的位置及组蛋白与DNA的相互关联水平(the level of association with DNA)都会影响PIC和TF绑定位点的可及性(accessibility)。
- 组蛋白不仅是静态实体，还可以进行化学修饰（染色质信号）从而影响其与DNA和其他蛋白质的结合（图1）。

![Signals important for transcription activation]({{ "/assets/img/blogs/tss-signals-and-features-2011-Valen-Eivind.png" | relative_url }})

启动子：与TSS相邻的区域中的*DNA序列及染色质中存在的调节信号*通常被称为“启动子”，启动子在调节转录中的起着重要作用。

TATA框（TATA box）：TSS周围的区域包含相似的序列信号，这些信号形成一个固定模板，该模板由TSS上游富含TA的模式（约30 nt）主导（称为TATA框）。含有TATA框的启动子较少（只有大约10–25％的启动子具有清晰的TATA模式）。

CpG岛[^7]：与基因组平均值相比，脊椎动物（如人类）的启动子通常具有更大的整体序列偏向性：CG二核苷酸(CG dinucleotides)的过量表达。这种富含CG的局部片段被称为CpG岛（CGI）。
- 40-70％的人类启动子带有CGI（这说明带有CGI的人类启动子还是比较少），这是因为脊椎动物的CG二核苷酸是甲基化事件的靶标，随着时间的推移由于脱氨作用会将Cs突变为Ts,导致CG二核苷酸在脊椎动物基因组中相对较少。但也有例外，比如与大脑和发育相关的特异性启动子通常富含CpG。
- 与TATA框不同，CGI通常与普遍表达的基因相关联(CGIs are often associated with ubiquitously expressed genes)。
- 启动子中TATA框和CG含量的差异表明了不同的调节机制。
  - 从人到酵母，TATA框是保守的，
  - 很大程度上，CGI是脊椎动物特异的，
  - 与含有CGIs的脊椎动物启动子相关的特征（包括染色质和TF信号）在真核生物（如酵母）中都是保守的。

![Sequence patterns commonly found around the transcription start site]({{ "/assets/img/blogs/tss-box1-2011-Valen-Eivind.png" | relative_url }})

哺乳动物基因组测序技术促进两种方法的开发:
- 基于染色质免疫沉淀（ChIP）[^8] 的全基因组方法（框2）：主要目的是检测结合DNA的蛋白质(DNA-bound proteins)，包括PIC[^3]的成分和具有特定化学修饰的组蛋白，
- 基于捕获RNA / cDNA的方法：主要目的是检测RNA的5'端边缘。

![Genome-wide methods for finding regulatory signals and promoters]({{ "/assets/img/blogs/tss-box2-2011-Valen-Eivind.png" | relative_url }})

启动子与TSS：
- 给定的启动子可能含有单个TSS，也可能包含分布在覆盖30-100 nt或更多区域的多个TSS，
- 给定的启动子的TSS分布通常在物种之间是保守的。

启动子分类（比较不同基因的TSS分布表明，启动子可以根据其TSS的分布大致分为“sharp”和“borad”两类启动子，图2）：
- sharp类启动子通常对应于具有TATA boxes的tissue-specific promoters （~25–30％的case有TATA盒）；
- borad类启动子具有CGI的过量表达，此类启动子通常在许多组织(tissue)中有活性。

这种基于TSS分布的二分法也反映在其他的的一些关于标准化CG含量遵循双峰分布的脊椎动物的研究中：<kbd>大多数启动子的标准化CG含量高或低，但很少有中间值</kbd>。基于CG含量的启动子分类还揭示了每个类别特有的其他特征：
- 高CG含量（HCG）[^9]的启动子对应于board类的启动子，因为它们通常无处不在的表达，而且经常与housekeeping基因相关（are often associated with housekeeping genes）。而在没有CGI的生物中board类启动子与特定motif(基序)的富集相关，而不是与一般CpG富集相关。
- 低CG含量（LCG）[^10]的启动子对应于sharp类的启动子，因为它们的表达更可能受到时间或空间的限制，并具有TATA框。与非TATA启动子相比，*它们往往在表达上具有更大的变异性*，非TATA启动子往往表现出稳定的表达。

不同启动子类别中影响TSS位置的信号主要是TFBS和染色质信号的不同（表1）。

![Summary of general properties of HCG/LCG promoters]({{ "/assets/img/blogs/tss-table1-2011-Valen-Eivind.png" | relative_url }})

TF结合位点
- TF通常分为通用TF[^11]和特定TF。通用TF与PIC相关，而特定TF通常与更远离TSS的位点结合。 特定TF在不同细胞之间的表达可能会有所不同（图1a，方框3）。 **在大多数文献中，“转录因子”是指后者的特定类别**。
- 因为TF的表达是与给定细胞中特定结合位点结合的前提，所以各个TF的表达水平和细胞特异性对于基因调节非常重要。
- HCG和LCG启动子的TFBS landscape明显不同：HCG减少了组织特异性TF的位点，而LCG却比较丰富。

![TFs: proteins binding DNA with a role in gene regulation]({{ "/assets/img/blogs/tss-box3-2011-Valen-Eivind.png" | relative_url }})

染色质信号Chromatin signals(核小体的数量)
DNA的可及性(accessibility)是转录调控的重要特征。 可及性取决于染色质状态以及DNA及其相互作用的蛋白的结构特性（框4）。早期关于人类，果蝇和酵母中核小体定位的研究表明，平均而言活跃转录的基因的TSS区耗尽了核小体，形成了一个单核宽的无核小体区域（NFR）。这个术语有点误导，因为后来证明NFR被不稳定结合的核小体占据，从而使调节机制更容易接近DNA。 在苍蝇和人类中，这种不稳定的核小体包含替代性的组蛋白变体H3.3和H2A.Z（代替了规范的H3和H2A组蛋白）。 同样，最近在酵母中也观察到不稳定的核小体，表明该特征对真核生物是保守的。

![Chromatin signals: chemical modification of DNA and chromatin]({{ "/assets/img/blogs/tss-box4-2011-Valen-Eivind.png" | relative_url }})

除了改变组蛋白类型，细胞还可以化学修饰组蛋白以改变其功能。 特别是在TSS周围的组蛋白中观察到了几种修饰，这些修饰与活跃或抑制基因(repressed)特别相关。 这些修饰通常称为<kbd>“表观遗传标记”</kbd>，尽管该术语也可以用于其他染色质特征。 与转录活性TSS相关的标记研究最多的是<kbd>H3K4me3和几个乙酰化标记</kbd>，这些标记主要集中在+1和-1核小体上。但是，这些标记的存在与转录活性并不完全相关，因为许多带有H3K4me3标记的启动子不会产生mRNA，而是将RNAPII保持在“平衡”状态。 这增加了第二层调节控制，许多小组已对伸长率的调节进行了深入研究。

<kbd>基于HCG / LCG或宽域（broad）/敏锐(sharp)的特性对启动子进行分类显示，在序列水平上观察到的二分法在染色质水平上也很明显。 在HCG /宽域中，NFR区域平均比在LCG /敏锐启动子中更为明显，并且通常与CGI一致。</kbd> 根据其核小体的占有率，这种分裂被称为“占用occupied”或“耗尽depleted”启动子。 然而，在没有CGI的情况下，均聚物homopolymeric poly（dA：dT）片段会形成刚性的DNA片段，不利于包裹组蛋白，这表明不同物种已开发出建立NFR的不同方法。

宽域启动子的下游区域通常显示规则间隔regular intervals（即定相phasing）排列的核小体这种典型模式(canonical pattern)，其中<kbd>TSS下游的第一个核小体位点很强（定位固定），随后的核小体位置逐渐变化。在NFR上游区域，人们也可以观察到核小体不太明显的分阶段定位</kbd>。即使在没有基因表达的情况下，这些启动子中的NFR侧翼的核小体也富含H​​2A.Z和H3K4me3和乙酰化标记，表明这些信号可能不是指示延伸，而是保持平衡。

<kbd>与HCG /宽域启动子相反，LCG /敏锐启动子通常缺乏NFR</kbd>。LCG /敏锐启动子具有占据该区域的相对稳定的核小体。同样，平均下游核小体的分布似乎比HCG的有序。核小体对NFR稳定占据的合理解释是，默认情况下，这些启动子中的大多数都不具有转录活性（这是因为许多启动子仅用于特定细胞），而这些沉默的启动子含有覆盖TSS的核小体的区域，这使得PIC无法访问。激活这种启动子必须涉及去除核小体。

在局部水平上，核心启动子区域内的序列信号，尤其是+/-1二核苷酸和–30附近的TATA样图谱，可以将高使用率的TSS与低使用率的TSS区分开。此发现基于以下假设：相同序列的信号在LCG和HCG启动子中具有同等的重要性。最近的一项研究表明，对于敏锐启动子，可以定义一组高度可预测TSS位置的TFBS，包括TATA等一般模式以及CREB和AP2等其他TF。然而，当尝试使用相同的框架分析广泛的启动子时，这更具挑战性，因为组蛋白修饰数据比序列特征更具信息性。

最近的两项研究表明，核心启动子区域周围的少数染色质标记可预测几种物种的启动子活性。研究中组蛋白标记最重要的排名最高的是H4K20me1，H3K27ac，H3K79me1和H2BK5ac ，而另一项发现H3K79me2，H3K79me3和H3K4me2标记是最预测的 。研究之间的差异可能归因于所使用的不同数据集，但也可能是由于标记中所含信息的大量冗余所致，因此可以从分析中删除多个标记，并且仍然可以预测高表达准确性。这可能是由于特定标记之间的高度相关性。回到上面讨论的启动子二分法，H3K27ac和H4K20me1是HCG中最预测的表达信号，而H3K4me3和H3K79me1在LCG中最预测。

总之，这些研究表明，可以使用序列和染色质信号来预测TSS的位置和活性，并查明每种类型中最重要的信号。

结束语
大量采用不同类型数据的研究至少支持两种主要类型的脊椎动物启动子，而这些启动子在选择TSS时所用的信号却大不相同。 HCG启动子通常具有许多TSS，表达广泛，并且较少依赖TF结合位点表达。这些启动子的高CG含量本身可以通过染色质信号引发启动子进行转录激活。相反，LCG启动子似乎遵循类似于转录起始的经典观点的调节逻辑，该启动子默认情况下是无活性的，但是可以被募集染色质重塑复合物的特定TFs激活，而后者又去除了核小体以暴露出更多的TFBS或TSS。

由于HCG / LCG启动子具有与每个类别相关的许多功能，因此可以使用广泛不同的数据集对启动子进行分类：从50种关注TSS分布的标签数据到DNA序列和染色质修饰。

[^1]: 转录起始位点（TSS）：在mRNA中转录的第一个基因组核苷酸。
[^2]: RNA聚合酶II（RNAPII）：读取DNA产生mRNA的酶。
[^3]: 起始复合物（PIC）：一种蛋白质复合物，在TSS上组装以启动转录，由RNAPII，通用TF，TAF以及许多其他蛋白质组成。
[^4]: 核心启动子：TSS周围的DNA区域，带有调节信号以结合PIC。
[^5]: 转录因子（TF）：与DNA绑定（或结合）并调节转录的蛋白质
[^6]: 转录因子结合位点（TFBS）：由TF绑定（或结合）的短DNA序列(通常为5-15 nt)。
[^7]: CpG是指胞嘧啶C后接鸟嘌呤G（“p”是指核苷酸之间的磷酸酯连接子），通常用CpG岛（CGI）来讨论。CpG岛（CGI）是指CpG过多的区域。由于这是一种统计属性，因此存在许多定义，通常基于标准化的CpG含量。
[^8]: 染色质免疫沉淀（ChIP）是一种针对具有特异性抗体的特定蛋白质提取特定DNA-蛋白质复合物的技术。可以与测序（ChIP-seq）或杂交技术（ChIP-chip）结合使用以研究绑定DNA(bound DNA)。
[^9]: 高CG（HCG）启动子：基于单个C和G核苷酸的数量，CpG含量高于预期的启动子区域。
[^10]: 低CG（LCG）低启动子：基于单个C和G核苷酸的数量，CpG含量低于或类似于基因组背景的启动子区域。
[^11]: 通用转录因子：无处不在的TF，通常是PIC的一部分。
[^12]: 归一化CG含量：根据给定的DNA序列中Cs和Gs的比例，将区域归一化到预期的区域的CpG含量，通常用于CGI的定义。
[^13]: 无核小体区域（NFR）：TSS上游紧邻的区域，平均而言，该区域耗尽了典型的核小体。

