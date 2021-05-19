---
title: "TermProject_mission01"
author: "박소영"
date: "5/11/2021"
mainfont: NanumGothic
output:
  html_document:
    keep_md: yes
    latex_engine: xelatex
    toc: true
    toc_float: true
    df_print: paged
---




# 0. Import Library and Data

```r
# libraries
library(dplyr)
library(ggplot2)
library(knitr)

# data
cnts <- read.table("read-counts.txt", sep="\t", header=T)
cnts$clip_enrichment <- cnts$CLIP.35L33G.bam/cnts$RNA.control.bam
cnts$rden_change <- cnts$RPF.siLin28a.bam/cnts$RNA.siLin28a.bam/
  (cnts$RPF.siLuc.bam/cnts$RNA.siLuc.bam)
mouselocal <- read.table("mouselocal.txt", sep="\t", header = T)
```


```r
# count data
kable(head(cnts))
```



|Geneid               |Chr                                |Start                                                   |End                                                     |Strand        | Length| CLIP.35L33G.bam| RNA.control.bam| RNA.siLin28a.bam| RNA.siLuc.bam| RPF.siLin28a.bam| RPF.siLuc.bam| clip_enrichment| rden_change|
|:--------------------|:----------------------------------|:-------------------------------------------------------|:-------------------------------------------------------|:-------------|------:|---------------:|---------------:|----------------:|-------------:|----------------:|-------------:|---------------:|-----------:|
|ENSMUSG00000102693.2 |chr1                               |3143476                                                 |3144545                                                 |+             |   1070|               0|               0|                0|             0|                0|             0|             NaN|         NaN|
|ENSMUSG00000064842.3 |chr1                               |3172239                                                 |3172348                                                 |+             |    110|               0|               0|                0|             0|                0|             0|             NaN|         NaN|
|ENSMUSG00000051951.6 |chr1;chr1;chr1;chr1;chr1;chr1;chr1 |3276124;3276746;3283662;3283832;3284705;3491925;3740775 |3277540;3277540;3285855;3286567;3287191;3492124;3741721 |-;-;-;-;-;-;- |   6094|               4|               1|                1|             1|                0|             0|               4|         NaN|
|ENSMUSG00000102851.2 |chr1                               |3322980                                                 |3323459                                                 |+             |    480|               3|               0|                0|             0|                0|             0|             Inf|         NaN|
|ENSMUSG00000103377.2 |chr1                               |3435954                                                 |3438772                                                 |-             |   2819|               0|               0|                0|             0|                0|             0|             NaN|         NaN|
|ENSMUSG00000104017.2 |chr1                               |3445779                                                 |3448011                                                 |-             |   2233|               0|               0|                0|             0|                0|             0|             NaN|         NaN|


```r
# protein localization data
kable(head(mouselocal))
```



|gene_id            |Gene.names           |type      |
|:------------------|:--------------------|:---------|
|ENSMUSG00000000001 |Gnai3                |cytoplasm |
|ENSMUSG00000000028 |Cdc45 Cdc45l Cdc45l2 |nucleus   |
|ENSMUSG00000000049 |Apoh B2gp1           |cytoplasm |
|ENSMUSG00000000058 |Cav2                 |cytoplasm |
|ENSMUSG00000000085 |Scmh1                |nucleus   |
|ENSMUSG00000000093 |Tbx2                 |nucleus   |
# 1. The Problem to Solve out
## 1.1. The raw graph
* 그냥 그리면 논문과 아주 많이 다른 그래프를 어떻게 더 프로세스할 것인가?

```r
#without any other process
ggplot(cnts, aes(x = log2(clip_enrichment), y = log2(rden_change))) +
  geom_point(size = 0.5, alpha=0.5)
```

```
## Warning: Removed 34128 rows containing missing values (geom_point).
```

![](mission01_files/figure-html/unnamed-chunk-4-1.png)<!-- -->


## 1.2 Explore more on the read count data

```r
colnames(cnts)
```

```
##  [1] "Geneid"           "Chr"              "Start"            "End"             
##  [5] "Strand"           "Length"           "CLIP.35L33G.bam"  "RNA.control.bam" 
##  [9] "RNA.siLin28a.bam" "RNA.siLuc.bam"    "RPF.siLin28a.bam" "RPF.siLuc.bam"   
## [13] "clip_enrichment"  "rden_change"
```

```r
par(mfrow = c(3,2))
hist(cnts$CLIP.35L33G.bam)
hist(cnts$RNA.control.bam)
hist(cnts$RNA.siLin28a.bam)
hist(cnts$RNA.siLuc.bam)
hist(cnts$RPF.siLin28a.bam)
hist(cnts$RPF.siLuc.bam)
```

![](mission01_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

* 히스토그램으로는 보기 어려우니 statistic summary로 보자

```r
colnames(cnts)
```

```
##  [1] "Geneid"           "Chr"              "Start"            "End"             
##  [5] "Strand"           "Length"           "CLIP.35L33G.bam"  "RNA.control.bam" 
##  [9] "RNA.siLin28a.bam" "RNA.siLuc.bam"    "RPF.siLin28a.bam" "RPF.siLuc.bam"   
## [13] "clip_enrichment"  "rden_change"
```

```r
summary(cnts[7:12])
```

```
##  CLIP.35L33G.bam     RNA.control.bam     RNA.siLin28a.bam   RNA.siLuc.bam     
##  Min.   :      0.0   Min.   :     0.00   Min.   :     0.0   Min.   :     0.0  
##  1st Qu.:      0.0   1st Qu.:     0.00   1st Qu.:     0.0   1st Qu.:     0.0  
##  Median :      0.0   Median :     0.00   Median :     0.0   Median :     0.0  
##  Mean   :    246.2   Mean   :    93.95   Mean   :   222.9   Mean   :   176.1  
##  3rd Qu.:     25.0   3rd Qu.:    16.00   3rd Qu.:    19.0   3rd Qu.:    15.0  
##  Max.   :2451339.0   Max.   :119734.00   Max.   :214518.0   Max.   :237513.0  
##  RPF.siLin28a.bam  RPF.siLuc.bam    
##  Min.   :      0   Min.   :      0  
##  1st Qu.:      0   1st Qu.:      0  
##  Median :      0   Median :      0  
##  Mean   :    193   Mean   :    234  
##  3rd Qu.:      5   3rd Qu.:      7  
##  Max.   :4698218   Max.   :5531458
```
* 모든 read count의 median이 다 0이다. 사분위수 아닌 십분위수로도 summary 해보자

```r
for (i in 7:12){
  print(paste0("========", colnames(cnts)[i], "========="))
  print(quantile(cnts[,i], seq(0,1,0.1)))
}
```

```
## [1] "========CLIP.35L33G.bam========="
##      0%     10%     20%     30%     40%     50%     60%     70%     80%     90% 
##       0       0       0       0       0       0       5      14      55     338 
##    100% 
## 2451339 
## [1] "========RNA.control.bam========="
##     0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
##      0      0      0      0      0      0      2      8     40    201 119734 
## [1] "========RNA.siLin28a.bam========="
##       0%      10%      20%      30%      40%      50%      60%      70% 
##      0.0      0.0      0.0      0.0      0.0      0.0      3.0      9.0 
##      80%      90%     100% 
##     58.0    468.2 214518.0 
## [1] "========RNA.siLuc.bam========="
##     0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
##      0      0      0      0      0      0      2      7     43    361 237513 
## [1] "========RPF.siLin28a.bam========="
##        0%       10%       20%       30%       40%       50%       60%       70% 
##       0.0       0.0       0.0       0.0       0.0       0.0       0.0       2.0 
##       80%       90%      100% 
##      18.0     151.2 4698218.0 
## [1] "========RPF.siLuc.bam========="
##      0%     10%     20%     30%     40%     50%     60%     70%     80%     90% 
##       0       0       0       0       0       0       1       3      22     210 
##    100% 
## 5531458
```

```r
print("================")
```

```
## [1] "================"
```

```r
sprintf("the number of current genes : %d", nrow(cnts))
```

```
## [1] "the number of current genes : 55359"
```

* 전체 55359개의 유전자 중 약 5000개 내외의 유전자를 남기고자하며, 따라서 cut-off를 read count 100정도롤 잡아본다.

```r
colnames(cnts)
```

```
##  [1] "Geneid"           "Chr"              "Start"            "End"             
##  [5] "Strand"           "Length"           "CLIP.35L33G.bam"  "RNA.control.bam" 
##  [9] "RNA.siLin28a.bam" "RNA.siLuc.bam"    "RPF.siLin28a.bam" "RPF.siLuc.bam"   
## [13] "clip_enrichment"  "rden_change"
```

```r
cnts2 <- cnts[apply(cnts[,7:12], 1, function (x) all(x >=100)),]

ggplot(cnts2, aes(x = log2(clip_enrichment), y = log2(rden_change))) +
  geom_point(size = 0.5, alpha=0.3)
```

![](mission01_files/figure-html/unnamed-chunk-8-1.png)<!-- -->
  
* 훨씬 나아 보인다.

# 2. Merge Subcellular Location
read count data의 Geneid를 protein localization데이터셋의 gene_id와 맞춰서 데이터셋을 합쳐준다.

```r
#extract gene_id form Geneid
cnts2 <- cnts2 %>% mutate(geneid2 = sub("(*)[.]\\d+", "\\1", Geneid))
cnts2$geneid2[1:5]
```

```
## [1] "ENSMUSG00000033845" "ENSMUSG00000033813" "ENSMUSG00000033793"
## [4] "ENSMUSG00000025907" "ENSMUSG00000051285"
```

```r
#merge two dataframe
colnames(cnts2)
```

```
##  [1] "Geneid"           "Chr"              "Start"            "End"             
##  [5] "Strand"           "Length"           "CLIP.35L33G.bam"  "RNA.control.bam" 
##  [9] "RNA.siLin28a.bam" "RNA.siLuc.bam"    "RPF.siLin28a.bam" "RPF.siLuc.bam"   
## [13] "clip_enrichment"  "rden_change"      "geneid2"
```

```r
df <- merge(mouselocal, cnts2[13:15], by.x = "gene_id", by.y = "geneid2")
dim(df)
```

```
## [1] 3287    5
```
* 두개의 데이터 세트에서 정보가 모두 있는 유전자들만 합쳤더니 데이터 사이즈가 3000개 조금 넘도록 또 많이 줄었다.  
한번 그래프가 어떻게 보이는 지 보자.

```r
ggplot(df, aes(x = log2(clip_enrichment), y = log2(rden_change))) +
  geom_point(size = 0.5, alpha=0.3)
```

![](mission01_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

```r
ggplot(df, aes(x = log2(clip_enrichment), y = log2(rden_change), color = type)) +
  geom_point(size = 0.5, alpha=0.7)
```

![](mission01_files/figure-html/unnamed-chunk-10-2.png)<!-- -->
  
* 여전히 나쁘지 않게 보인다. 이제 df를 가지고 그래프 최적화 추가작업을 진행한다. 

# 3. Optimize Visualization
## 3.1. Random sampling for each localization

```r
localtype <- unique(df$type)
df2 <- data.frame()
for (t in localtype) {
  temp <- df %>% filter(type == t)
  set.seed(123)
  temp <- temp %>% filter(gene_id %in%
                            sample(temp$gene_id, 300))
  df2 <- rbind(df2, temp)
}

ggplot(df2, aes(x = log2(clip_enrichment), y = log2(rden_change), color = type)) +
  geom_point(size = 1, alpha=0.8)
```

![](mission01_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

## 3.2. Customizing the graph

```r
ggplot(df2, aes(x = log2(clip_enrichment), y = log2(rden_change), color = type)) +
  geom_point(size = 1, alpha=0.8) +  lims(y=c(-3,2),x=c(-4,4.5)) +
  labs(x = expression("LIN28A CLIP enrichment (log"[2]*")"),
       y = expression("Ribosome density change upon "*italic(Lin28a)*" knockdown (log"[2]*")")) +
  theme_bw() + theme(legend.position = c(0.25,0.87), legend.title = element_blank(),
                     legend.background = element_rect(colour = "#888888", size =0.5),
                     legend.text = element_text(size = 10)) + 
  scale_color_manual(values = c("#008000", "#DC143C", "#1E90FF")) +
  coord_fixed(ratio = 2)
```

![](mission01_files/figure-html/unnamed-chunk-12-1.png)<!-- -->
