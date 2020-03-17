---
title: "组学之工具"
date: 2020-03-04 10:14:00 +0800
categories: [Genomics, Tools]
tags: [Omics tools]
---

# 比对 bwa
常见的是使用 `bam` 软件进行比对。根据数据情况的不同， `bam` 的参数设置也不同。

* 如果测序长度 >70bp, 则使用 `bwa mem`；如果测序长度较短（如35bp，NIPT常用的读长）， 则使用 `bwa aln`；

* 如果是单端比对，则使用 `bwa samse`； 如果是双端比对，则使用 `bwa sampe`。

例如对于NIPT短读长(35bp)，单端比对，使用如下比对命令（命令的顺序固定）：
```shell 
bwa aln -n 0 -e 0 -k 0 -t 16 hg19.fa s1.fq.gz > s1.sai
bwa samse -n -1 hg19.fa s1.sai s1.fq.gz > s1.sam
```

# 生成bam，排序并建索引
一旦比对生成 `.sam`文件，需要使用 `samtools` 对比对数据进行排序，并生成 `.bam` 文件，同时建立索引 `index`。
```shell
samtools view s1.sam -bSh > s1.bam
samtools sort -@ 16 s1.bam -o s1.sorted.bam
samtools index s1.sorted.bam
```

# 重复reads过滤
```shell
samtools rmdup -s s1.sorted.bam s1.rmdups.bam
samtools view -F 4 s1.rmdups.bam -bSh > s1.final.bam
samtools index s1.final.bam
```

deeptools-3.3.2 deeptoolsintervals-0.1.9 matplotlib-3.1.3 plotly-4.5.2 retrying-1.3.3
pysam-0.15.2(pysam.__samtools_version__=1.9)

# 术语
* read: 测序仪返回的原始序列.一个read可以包括多个segment。read之间的先后顺序表示被测序仪读到的时间前后关系.
* segment: 一段连续的序列或子序列
linear alignment: 线性联配表示一个read比对到单个参考序列，可以存在插入，缺失，跳过(skip),剪切(clip), 但是不存在方向改变的情况(比如说一部分和正向链联配，另一个位置则是和负向链联配）。最简单的判断的方式就是，一个linear alignment只用一行记录。
* chimeric alignment: 嵌合联配需要多行记录。比如说r003第一个记录是后6个匹配，第二个记录则是反向序列的后5个匹配。第一个被称之为"representative",其他都是"supplementary"
read alignment: 无论是linear alignment, 还是chimeric alignment, 只要能完整表示一个read，都成为是read alignment
* multiple mapping: 由于存在重复区，一个read 可能比对到参考基因组的不同区域。其中一个被认为是primary，其他都是secondary.
* 两个系统|1-based coordinate system（SAM,VCF,GFF,wiggle)和0-based coordinate system(BAM, BCFv2, BED, PSL).自行用R和Python感受一下两者的不同。
> chimeric alignment 可能是结构变异，基因融合，参考序列误组装，RNA-Seq，实验protocol等因素造成。对于chimeric alignment的里面每一个linear alignment而言，由于相互之前不存在重叠，故而联配质量较高，适合用于SNP/INDEL calling.相反, multiple mapping则是因为重复造成(read越长出现的概率越低), 相互之间存在重叠，仅有其中一条有最优的匹配，其他联配质量过低会被SNP/INDEL caller忽略。


# 参考

1. [SAMtools: SAM格式的处理利器](https://www.jianshu.com/p/8d01019f33f2)

http://www.bio-info-trainee.com/4452.html


https://zhuanlan.zhihu.com/p/29456819

https://www.melbournebioinformatics.org.au/tutorials/tutorials/var_detect_advanced/var_detect_advanced_background/