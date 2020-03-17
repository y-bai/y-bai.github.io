---
title: "组学之deeptools之一: computeGCBias"
date: 2020-03-04 10:14:00 +0800
categories: [Genomics, Tools]
tags: [Omics tools]
---

# deeptools
这里使用的版本为deeptools = 3.3.2, pysam = 0.15.2(使用`pip list`查看,)
# computeGCBias

参考 https://github.com/deeptools/deepTools/blob/3.3.2/deeptools/computeGCBias.py

## 参数
1. 必须参数
    
    * `-b`: 待处理的`.bam`或`.cram`文件。
    * `--effectiveGenomeSize`: 参考基因组去的有效大小（在[deeptools官网](https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html)上可以查到）。
    * `-g`: 参考基因组路径(`.2bit`格式）,可以通过如下进行转换：
    ```shell
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
    chmod +x faToTwoBit
    faToTwoBit hg38.fasta hg38.2bit
    ```
    * `-freq`：保存结果的文件路径。
2. 可选参数
    * `-l`: 测序模板长度(`pysam`中`AlignedSegment`对象的`template_length`值)，**它不是read_length(不是pysam中得到query_length或read_length)**。这个值可以通过deeptools中的`get_read_and_fragment_length()`方法可以得到。但是在实际计算过程中没有得到`template_length`，只得到了`read_length`（这是因为`.cram`没有存这个值，即在`.cram`中设置为0）。这种情况下需要直接显式指定。关于`template`和`read`可以参看https://www.jianshu.com/p/9c99e09630da 和 https://blog.csdn.net/xcaryyz/article/details/79257604。

        参看https://pysam.readthedocs.io/en/latest/usage.html#usage 发现 **35bp**的测序长度的模板长度是**167**。因为NIPT测序是用cfDNA,所以一般测序长度是160bp-180bp,而一般DNA的建库长度为300bp。

    * `--sampleSize`：采样大小，默认值是`5e7`.
    * `--extraSampling`: BED 文件，指明了需要额外采样的区域（这些区域在基因组中可能被过低表征），默认值是`None`。
    * `--biasPlot`：保存GC bias图像路径（.png格式）.默认是`None`,即不存GC bias的图像。
    * `--regionSize`: 区域大小（在计算每个GC对应的reads时的bin大小, 在computeGCBias里是画图的时候使用）。默认是300bp。
    * `-r`：指定染色体号，默认是`None`,即使用全部23条染色体。
    * `-p`：并行计算的CPU核的数量，默认是1.
    * `-v`： 参看处理过程信息

[注意]几个影响结果的参数：
* `-l`参数：由于项目使用的`.cram`中没有指定测序建库的长度，所以在计算中必须要显式指定，不能使用读长(35bp)，否则计算结果是错的。cfDNA的建库长度（即template length）一般是160bp-180bp，这里设置为180bp。

* `--effectiveGenomeSize`: 这是参考基因组的有效大小,可以使用如下软件计算得出(只计算一个DNA的一个链)
```shell
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faCount
chmod +x faCount
./faCount hg38.fasta -summary
```

结果如下：

```
#seq    len     A       C       G       T       N       cpg
total   3217346917      896424917       622999561       625196255       898832735       173893449       30914151
prcnt   1.0     0.2786  0.1936  0.1943  0.2794  0.0540  0.0096
```

由上面的结果可知，hg38的原始大小是3217346917bp，这个大小是等于`A+C+G+T+N`的结果。

根据https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html 上 option1的说明
> 1. The number of non-N bases in the genome.
> 2. The number of regions (of some size) in the genome that are uniquely mappable (possibly given some maximal edit distance).

那么这个时候`--effectiveGenomeSize`应该是`3217346917-173893449=3043453468`。 但是官网用的是`2913022398`,不知道还去除了什么点位（可能是由于原始参考基因组中出了1-22、X、Y、M染色体，还有其他染色体，比如`chr1_KI270706v1_random=175055bp`等,他们都有长度.）。使用下面的代码可以输入每个染色体的长度

```python
# chromSizes: list of tuples
chromSizes = [(bam.references[i], bam.lengths[i])
                  for i in range(len(bam.references))]
chromSizes = [x for x in chromSizes if x[0] in tbit.chroms()]
```

结果如下(部分)：

```
[('chr1', 248956422), ('chr2', 242193529), ('chr3', 198295559), ('chr4', 190214555), ('chr5', 181538259), ('chr6', 170805979), ('chr7', 159345973), ('chr8', 145138636), ('chr9', 138394717), ('chr10', 133797422), ('chr11', 135086622), ('chr12', 133275309), ('chr13', 114364328), ('chr14', 107043718), ('chr15', 101991189), ('chr16', 90338345), ('chr17', 83257441), ('chr18', 80373285), ('chr19', 58617616), ('chr20', 64444167), ('chr21', 46709983), ('chr22', 50818468), ('chrX', 156040895), ('chrY', 57227415), ('chrM', 16569), ('chr1_KI270706v1_random', 175055), ('chr1_KI270707v1_random', 32032), ('chr1_KI270708v1_random', 127682),...
```
从chr1到chr22, chrX, chrY, chrM的总长度是3088286401.

这里`--effectiveGenomeSize`使用官网提供的2913022398。

* `--sampleSize`: 这里默认的是`5e7`。从源代码看这个参数主要有三个地方用到了：

    1. 用`possion`估计最大和最小read数量的时候的p-value,代码如下：

        ```python
            if args.fragmentLength:
            fragment_len_dict = \
                {'median': args.fragmentLength}

        else:
            fragment_len_dict, __ = \
                get_read_and_fragment_length(args.bamfile, None,
                                             numberOfProcessors=args.numberOfProcessors,
                                             verbose=args.verbose)
            if not fragment_len_dict:
                print("\nPlease provide the fragment length used for the "
                      "sample preparation.\n")
                exit(1)

            fragment_len_dict = {'median': int(fragment_len_dict['median'])}

        chrNameBitToBam = tbitToBamChrName(list(tbit.chroms().keys()), bam.references)

        bam, mapped, unmapped, stats = bamHandler.openBam(global_vars['bam'], returnStats=True, nThreads=args.numberOfProcessors)

        global_vars['total_reads'] = mapped

        global_vars['reads_per_bp'] = \
            float(global_vars['total_reads']) / args.effectiveGenomeSize

        confidence_p_value = float(1) / args.sampleSize

        # chromSizes: list of tuples
        chromSizes = [(bam.references[i], bam.lengths[i])
                     for i in range(len(bam.references))]
        chromSizes = [x for x in chromSizes if x[0] in tbit.chroms()]

        # use poisson distribution to identify peaks that should be discarted.
        # I multiply by 4, because the real distribution of reads
        # vary depending on the gc content
        # and the global number of reads per bp may a be too low.
        # empirically, a value of at least 4 times as big as the
        # reads_per_bp was found.
        # Similarly for the min value, I divide by 4.
        global_vars['max_reads'] = poisson(4 * global_vars['reads_per_bp'] * fragment_len_dict['median']).isf(confidence_p_value)
        # this may be of not use, unless the depth of sequencing is really high
        # as this value is close to 0
        global_vars['min_reads'] = poisson(0.25 * global_vars['reads_per_bp'] * fragment_len_dict['median']).ppf(confidence_p_value)

        ```

    2. 计算`stepSize`的时候, `stepSize`作为参数计算`GC Content`

        ```python
        # the GC of the genome is sampled each stepSize bp.
        stepSize = max(int(global_vars['genome_size'] / args.sampleSize), 1) 

        data = tabulateGCcontent(fragment_len_dict,
                                 chrNameBitToBam, stepSize,
                                 chromSizes,
                                 numberOfProcessors=args.numberOfProcessors,
                                 verbose=args.verbose,
                                 region=args.region)

        np.savetxt(args.GCbiasFrequenciesFile.name, data)
        ```

        这里的`global_vars['genome_size']`就是参考基因组的原始大小，即`'genome_size': 3217346917`。上面的结果是`stepSize=64`。


    3. 计算GC含量含量的时候
    
        ```python
        if sum(F_gc) == 0:
            sys.exit("No fragments included in the sampling! Consider decreasing (or maybe increasing) the --sampleSize parameter")
        scaling = float(sum(N_gc)) / float(sum(F_gc))
        ```

        从上面2和3可以看出，计算GC含量是根据`stepSize`计算的。但是如果参数`--sampleSize`使用默认的`5e7`的话，`stepSize`为64。


## 如何计算GC含量

根据deepTools提供的`Tester`例子分析, 例子中使用的参数如下：

```python

fragmentLength = {'median': 10}, 
chrNameBitToBam = {'chr2L': '2L'} #这个对应参考基因组和bam中染色体的名字对应关系，例子中得到参考基因组只有一个染色体
stepSize = 2
chromSizes = [('2L', 23011544)], 1) # 因为bam对应的reference只有一个染色体2L， 对应代码如下：
###
bam = bamHandler.openBam(global_vars['bam'])
chromSizes = [(bam.references[i], bam.lengths[i])
                    for i in range(len(bam.references))]
###

# 其他的参数设置如下
global_vars = {'2bit': self.tbitFile,
                       'bam': self.bamFile,
                       'filter_out': None,
                       'mappability': self.mappability,
                       'extra_sampling_file': None,
                       'max_reads': 5,
                       'min_reads': 0,
                       'min_reads': 0,
                       'reads_per_bp': 0.3,
                       'total_reads': mapped, # 这里的计算结果是 201
                       'genome_size': sum(tbit.chroms().values()), # 这里的计算结果是 1000
                       'verbose': True
                       }
```

根据以上参数调用`tabulateGCcontent()`函数， 代码如下：

```python 
def tabulateGCcontent(fragmentLength, chrNameBitToBam, stepSize,
                      chromSizes, numberOfProcessors=None, verbose=False,
                      region=None):
    r"""
    Subdivides the genome or the reads into chunks to be analyzed in parallel
    using several processors. This codes handles the creation of
    workers that tabulate the GC content for small regions and then
    collects and integrates the results
    >>> test = Tester()
    >>> arg = test.testTabulateGCcontent()
    >>> res = tabulateGCcontent(*arg)
    >>> res
    array([[  0.        ,  18.        ,   1.        ],
           [  3.        ,  63.        ,   0.45815996],
           [  7.        , 159.        ,   0.42358185],
           [ 25.        , 192.        ,   1.25278115],
           [ 28.        , 215.        ,   1.25301422],
           [ 16.        , 214.        ,   0.71935396],
           [ 12.        ,  95.        ,   1.21532959],
           [  9.        ,  24.        ,   3.60800971],
           [  3.        ,  11.        ,   2.62400706],
           [  0.        ,   0.        ,   1.        ],
           [  0.        ,   0.        ,   1.        ]])
    """
    global global_vars
    chrNameBamToBit = dict([(v, k) for k, v in chrNameBitToBam.items()])
    # 这一步的结果是 chrNameBamToBit: {'2L': 'chr2L'}， 就是对传入的chrNameBitToBam做了一下转换
    
    chunkSize = int(min(2e6, 4e5 / global_vars['reads_per_bp']))
    # 这一步就是把整个基因组切分成一个个的chunk, 后面对每个chunk进行mapReduce计算。这里chunkSize是根据`reads_per_bp`来计算的，最小也是2e6. 这里的结果是 chunkSize=1333333

    chromSizes = [(k, v) for k, v in chromSizes if k in list(chrNameBamToBit.keys())]
    # 这里输出染色体的大小。 结果是chromSizes:=[('2L', 23011544)]

    # 调用 mapReduce
    imap_res = mapReduce.mapReduce((stepSize,
                                    fragmentLength, chrNameBamToBit,
                                    verbose),
                                   tabulateGCcontent_wrapper,
                                   chromSizes,
                                   genomeChunkLength=chunkSize,
                                   numberOfProcessors=numberOfProcessors,
                                   region=region)

    for subN_gc, subF_gc in imap_res:
        try:
            F_gc += subF_gc
            N_gc += subN_gc
        except NameError:
            F_gc = subF_gc
            N_gc = subN_gc

    if sum(F_gc) == 0:
        sys.exit("No fragments included in the sampling! Consider decreasing (or maybe increasing) the --sampleSize parameter")
    
    # 这部分就是论文https://academic.oup.com/nar/article/40/10/e72/2411059里的公式了
    scaling = float(sum(N_gc)) / float(sum(F_gc))

    R_gc = np.array([float(F_gc[x]) / N_gc[x] * scaling
                     if N_gc[x] and F_gc[x] > 0 else 1
                     for x in range(len(F_gc))])

    data = np.transpose(np.vstack((F_gc, N_gc, R_gc)))
    return data

```

下面来分析`mapReduce`里的`tabulateGCcontent_worker()`函数。对于`tabulateGCcontent()`中计算出的每一个`chunk`， `mapReduce`都会调用一次`tabulateGCcontent_worker()`：
```python
def tabulateGCcontent_worker(chromNameBam, start, end, stepSize,
                             fragmentLength,
                             chrNameBamToBit, verbose=False):
    r""" given genome regions, the GC content of the genome is tabulated for
    fragments of length 'fragmentLength' each 'stepSize' positions.

    >>> test = Tester()
    >>> args = test.testTabulateGCcontentWorker()
    >>> N_gc, F_gc = tabulateGCcontent_worker(*args)

    The forward read positions are:
    [1,  4,  10, 10, 16, 18]
    which correspond to a GC of
    [1,  1,  1,  1,  2,  1]
    The evaluated position are
    [0,  2,  4,  6,  8, 10, 12, 14, 16, 18]
    the corresponding GC is
    [2,  1,  1,  2,  2,  1,  2,  3,  2,  1]

    >>> print(N_gc)
    [0 4 5 1]
    >>> print(F_gc)
    [0 4 1 0]
    >>> test.set_filter_out_file()
    >>> chrNameBam2bit =  {'2L': 'chr2L'}

    Test for the filter out option
    >>> N_gc, F_gc = tabulateGCcontent_worker('2L', 0, 20, 2,
    ... {'median': 3}, chrNameBam2bit)
    >>> test.unset_filter_out_file()

    The evaluated positions are
    [ 0  2  8 10 12 14 16 18]
    >>> print(N_gc)
    [0 3 4 1]
    >>> print(F_gc)
    [0 3 1 0]

    Test for extra_sampling option
    >>> test.set_extra_sampling_file()
    >>> chrNameBam2bit =  {'2L': 'chr2L'}
    >>> res = tabulateGCcontent_worker('2L', 0, 20, 2,
    ... {'median': 3}, chrNameBam2bit)

    The new positions evaluated are
    [0, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18]
    and the GC is
    [2, 1, 1, 0, 1, 2, 2, 1,  2,  3,  2,  1]
    >>> print(res[0])
    [1 5 5 1]
    >>> print(res[1])
    [0 5 1 0]

    """
    if start > end:
        raise NameError("start %d bigger that end %d" % (start, end))

    chromNameBit = chrNameBamToBit[chromNameBam]

    # 对于每一次调用，该函数的输入参数是：
    # 2L 0 1333333 2 {'median': 10} {'2L': 'chr2L'} False
    # 2L 1333333 2666666 2 {'median': 10} {'2L': 'chr2L'} False
    # ...

    # array to keep track of the GC from regions of length 'fragmentLength'
    # from the genome. The index of the array is used to
    # indicate the gc content. The values inside the
    # array are counts. Thus, if N_gc[10] = 3, that means
    # that 3 regions have a gc_content of 10.
    # 通过以下的分析得知，这里的regions的大小是fragmentLength，按照这个大小从以stepSize在窗口上滑动(这里有overlap)
    subN_gc = np.zeros(fragmentLength['median'] + 1, dtype='int')
    subF_gc = np.zeros(fragmentLength['median'] + 1, dtype='int')

    tbit = py2bit.open(global_vars['2bit'])
    bam = bamHandler.openBam(global_vars['bam'])
    peak = 0
    startTime = time.time()

    if verbose:
        print("[{:.3f}] computing positions to "
              "sample".format(time.time() - startTime))

    positions_to_sample = getPositionsToSample(chromNameBit,
                                               start, end, stepSize)
    # 这里计算的positions_to_sample的结果是（对应每一个chunk）
    # positions_to_sample:[      0       2       4 ... 1333328 1333330 1333332]
    # positions_to_sample:[1333333 1333335 1333337 ... 2666661 2666663 2666665]
    # ...
    # 由此可以看出stepSize实际是用在这里了

    read_counts = []
    # Optimize IO.
    # if the sample regions are far apart from each
    # other is faster to go to each location and fetch
    # the reads found there.
    # Otherwise, if the regions to sample are close to
    # each other, is faster to load all the reads in
    # a large region into memory and consider only
    # those falling into the positions to sample.
    # The following code gets the reads
    # that are at sampling positions that lie close together

    # 这里就说明了stepSize不能设置的过大，最大是1000
    if np.mean(np.diff(positions_to_sample)) < 1000:
        start_pos = min(positions_to_sample)
        end_pos = max(positions_to_sample)
        if verbose:
            print("[{:.3f}] caching reads".format(time.time() - startTime))

        # 这里使用pysam 的fetch()函数获取给定染色体的给定区域比对有的reads.
        # start and stop denote an interval in 0-based, half-open coordinates (like BED files and Python slices).
        # 而pysam是以位置0开始计算，而且依据python的语法，末尾是不计算在内的。
        # pysam必须要手动使得其包含末尾的位点，这就是为什么这里的end_pos要+1的原因。

        # np.bincount的结果是计算每个pos的的数量，counts的结果意义是:index是相对position, 对应的值是数量,如count[0]=1,表示pos为0的read有1个。这里的pos是每个read的起始位置。这里read_count的计算其实和stepSize是没有关系的。
        counts = np.bincount([r.pos - start_pos
                              for r in bam.fetch(chromNameBam, start_pos,
                                                 end_pos + 1)
                              if not r.is_reverse and not r.is_unmapped and r.pos >= start_pos],
                             minlength=end_pos - start_pos + 2)

        read_counts = counts[positions_to_sample - min(positions_to_sample)]
        if verbose:
            print("[{:.3f}] finish caching reads.".format(
                time.time() - startTime))

    countTime = time.time()

    # 这里就和stepSize有关了
    # positions_to_sample 是genemo的绝对位置
    c = 1
    for index in range(len(positions_to_sample)):
        i = positions_to_sample[index]
        # stop if the end of the chromosome is reached
        # tbit.chroms(chromNameBit) 是参考基因组对应染色体长度，根据参数的设置，这里是1000
        # 所以例子中只计算前1000的位点。
        if i + fragmentLength['median'] > tbit.chroms(chromNameBit):
            break

        try:
            # 从这里可以看出GC是计算自参考基因的量，而且在计算参考基因组上的GC时用的窗口大小是fragmentLength的大小, 如果fraction=False，那么GC是绝对数量
            gc = getGC_content(tbit, chromNameBit, int(i), int(i + fragmentLength['median']), fraction=False)
        except Exception as detail:
            if verbose:
                print(detail)
            continue

        subN_gc[gc] += 1 # 这个是存储GC数量为gc的fragmentLength有多少个

        # count all reads at position 'i'
        if len(read_counts) == 0:  # case when no cache was done
            num_reads = len([x.pos for x in bam.fetch(chromNameBam, i, i + 1)
                             if x.is_reverse is False and x.pos == i])
        else:
            num_reads = read_counts[index]
        # 计算当前位置的read数量，并判断是否是达到了peak
        if num_reads >= global_vars['max_reads']:
            peak += 1
            continue

        # 如果没有达到peak
        subF_gc[gc] += num_reads # 这里存储的是GC数量为gc的reads数有多少。
        if verbose:
            if index % 50000 == 0:
                endTime = time.time()
                print("%s processing %d (%.1f per sec) @ %s:%s-%s %s" %
                      (multiprocessing.current_process().name,
                       index, index / (endTime - countTime),
                       chromNameBit, start, end, stepSize))
        c += 1

    if verbose:
        endTime = time.time()
        print("%s processing %d (%.1f per sec) @ %s:%s-%s %s" %
              (multiprocessing.current_process().name,
               index, index / (endTime - countTime),
               chromNameBit, start, end, stepSize))
        print("%s total time %.1f @ %s:%s-%s %s" % (multiprocessing.current_process().name,
                                                    (endTime - startTime), chromNameBit, start, end, stepSize))

    return(subN_gc, subF_gc)
```


以上分析完成。总结就是：
1. 按照fregmentLenght有over-lap的滑动计算GC

现在来看GCbias画图

```python
if args.biasPlot:
    reads_per_gc = countReadsPerGC(args.regionSize,
                                       chrNameBitToBam, stepSize * 10,
                                       chromSizes,
                                       numberOfProcessors=args.numberOfProcessors,
                                       verbose=args.verbose,
                                       region=args.region)
    if args.plotFileFormat == "plotly":
        plotlyGCbias(args.biasPlot, data, reads_per_gc, args.regionSize)
    else:
        plotGCbias(args.biasPlot, data, reads_per_gc, args.regionSize, image_format=args.plotFileFormat)

```
首先调用`countReadsPerGC()`计算每个GC的reads:

```python
def countReadsPerGC(regionSize, chrNameBitToBam, stepSize,
                    chromSizes, numberOfProcessors=None, verbose=False,
                    region=None):
    r"""
    Computes for a region of size regionSize, the GC of the region
    and the number of reads that overlap it.
    >>> test = Tester()
    >>> arg = test.testCountReadsPerGC()
    >>> reads_per_gc = countReadsPerGC(*arg)
    >>> reads_per_gc[0:5,:]
    array([[132.        ,   0.44      ],
           [132.        ,   0.44      ],
           [133.        ,   0.44      ],
           [134.        ,   0.43666667],
           [134.        ,   0.44      ]])
    """
    global global_vars

    chrNameBamToBit = dict([(v, k) for k, v in chrNameBitToBam.items()])
    chunkSize = int(min(2e6, 4e5 / global_vars['reads_per_bp']))

    imap_res = mapReduce.mapReduce((stepSize,
                                    regionSize, chrNameBamToBit,
                                    verbose),
                                   countReadsPerGC_wrapper,
                                   chromSizes,
                                   genomeChunkLength=chunkSize,
                                   numberOfProcessors=numberOfProcessors,
                                   region=region)

    reads_per_gc = []
    for sub_reads_per_gc in imap_res:
        reads_per_gc += sub_reads_per_gc

    reads_per_gc = np.asarray(reads_per_gc)
    return reads_per_gc
```
`regionSize`使用默认值300， 在计算GC bias中，只用于画图，且作为binsize使用。


```python
def countReadsPerGC_worker(chromNameBam,
                           start, end, stepSize, regionSize,
                           chrNameBamToBit, verbose=False):
    """given a genome region defined by
    (start, end), the GC content is quantified for
    regions of size regionSize that are contiguous
    """

    chromNameBit = chrNameBamToBit[chromNameBam]
    tbit = py2bit.open(global_vars['2bit'])
    bam = bamHandler.openBam(global_vars['bam'])
    c = 1
    sub_reads_per_gc = []
    positions_to_sample = getPositionsToSample(chromNameBit,
                                               start, end, stepSize)

    for index in range(len(positions_to_sample)):
        i = positions_to_sample[index]
        # stop if region extends over the chromosome end
        if tbit.chroms(chromNameBit) < i + regionSize:
            break

        try:
            # 这里计算的是Fraction
            # regionSize被当成binSize来使用 
            gc = getGC_content(tbit, chromNameBit, int(i), int(i + regionSize))
        except Exception as detail:
            if verbose:
                print("{}:{}-{}".format(chromNameBit, i, i + regionSize))
                print(detail)
            continue
        numberReads = bam.count(chromNameBam, i, i + regionSize)
        sub_reads_per_gc.append((numberReads, gc))
        c += 1

    return sub_reads_per_gc
```

```python
def plotGCbias(file_name, frequencies, reads_per_gc, region_size, image_format=None):
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['svg.fonttype'] = 'none'
    import matplotlib.pyplot as plt

    # prepare data for boxplot
    reads, GC = reads_per_gc.T
    reads_per_gc, bin_labels = bin_by(reads, GC, nbins=100)
    # 这里源代码是 0.2 <= x <= 0.7， 在这里改成0.3 <= x <= 0.6
    # 当 -regionSize 过大（比如50k）， 或--sampleSize过小导致stepSize>1000时，0.2 <= x <= 0.7都会导致报错
    # 当 -regionSize=50k，0.3 <= x <= 0.6运行正常
    to_keep = [idx for idx, x in enumerate(bin_labels) if 0.3 <= x <= 0.6]
    reads_per_gc = [reads_per_gc[x] for x in to_keep]
    bin_labels = [bin_labels[x] for x in to_keep]

    title = "reads per regions of {} bp".format(region_size)
    fig = plt.figure(figsize=(6, 8))
    ax1 = fig.add_subplot(211, title=title)
    ax2 = fig.add_subplot(212,
                          title='normalized observed/expected read counts')

    # make boxplot

    bp = ax1.boxplot(reads_per_gc, notch=0, patch_artist=True)
    plt.setp(bp['boxes'], color='black', facecolor='LightGreen')
    plt.setp(bp['medians'], color='black')
    plt.setp(bp['whiskers'], color='black', linestyle='dashed')
    plt.setp(bp['fliers'], marker='None')
    # get the whisker that spands the most
    y_max = max([x.get_data()[1][1] for x in bp['whiskers']])
    ax1.set_ylim(0 - (y_max * 0.05), y_max * 1.05)
    ax1.set_ylabel('Number of reads')
    ax1.set_xlabel('GC fraction')

    xticks = [idx for idx, x in enumerate(bin_labels) if int(x * 100) % 10 == 0]

    ax1.set_xticks(xticks)
    ax1.set_xticklabels(["{:.1f}".format(bin_labels[x]) for x in xticks])

    x = np.linspace(0, 1, frequencies.shape[0])
    y = np.log2(frequencies[:, 2])
    ax2.plot(x, y, color='#8c96f0')
    ax2.set_xlabel('GC fraction')
    ax2.set_ylabel('log2ratio observed/expected')
    ax2.set_xlim(0.2, 0.7)
    y_max = max(y[np.where(x >= 0.2)[0][0]:np.where(x <= 0.7)[0][-1] + 1])
    y_min = min(y[np.where(x >= 0.2)[0][0]:np.where(x <= 0.7)[0][-1] + 1])
    if y_max > 0:
        y_max *= 1.1
    else:
        y_max *= 0.9
    if y_min < 0:
        y_min *= 1.1
    else:
        y_min *= 0.9
    ax2.set_ylim(y_min, y_max)
    plt.tight_layout()
    plt.savefig(file_name, bbox_inches='tight', dpi=100, format=image_format)
    plt.close()

```





