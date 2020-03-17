






http://wiki2.gene-go.com/bbs/thread-508-1-1.html

重点：
https://pzweuj.github.io/2019/03/20/NIPT.html

http://blog.sina.com.cn/s/blog_83f77c940102valj.html

http://bioinformaticsinstitute.ru/sites/default/files/kozyulina_20170916.pdf


重点：
https://blog.csdn.net/Cassiel60/article/details/88742058
https://cran.r-project.org/web/packages/NIPTeR/vignettes/NIPTeR.html

https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2557-8

https://deeptools.readthedocs.io/en/develop/content/tools/computeGCBias.html

Sensitivity of Noninvasive Prenatal Detection of Fetal Aneuploidy from Maternal Plasma Using Shotgun Sequencing Is Limited Only by Counting Statistics

# bias
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3906086/

Investigating and Correcting Plasma DNA Sequencing Coverage Bias to Enhance Aneuploidy Discovery

> One such bias is the Guanine-Cytosine content (GC) bias which is uni-modal in nature as the read coverage is maximized for genomic regions with 40–50% GC with coverage decreasing at the extreme values. The current tests predominantly incorporate a locally weighted scatterplot smoothing (LOESS) correction step as described by Alkan et al. [34] to correct for GC. This involves binning the read counts into non-overlapping windows and calculating the GC content of these windows (a 50 Kb bin size is in wide use). A LOESS curve is fitted to the plot of bin counts vs. bin GC to obtain predicted values that are used to correct the raw counts.

> A second source of bias stemming from the alignment of short read data is the inability of sequence reads to map unambiguously or uniquely to highly repetitive regions of the genome. This mappability bias is generally dealt with by removal of such regions (annotated in the RepeatMasker database [35]) rather than quantifying the mappability and correcting for it.

> We also show an additional source of bias due to fragmentation that is inflated in plasma cfDNA when compared to genomic DNA. The fragmentation effect is the position-specific pattern of nucleotides around DNA fragment ends. Plasma DNA shows different and stronger patterns than genomic DNA due to the fragments originating from a biological process and possible nuclease activity rather than a random shearing during the DNA fragmentation step of library preparation. This bias is investigated in depth to illustrate its relationship with the GC bias. We then extend the improved GC correction to incorporate the fragmentation bias.

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3378858/
Summarizing and correcting the GC content bias in high-throughput sequencing