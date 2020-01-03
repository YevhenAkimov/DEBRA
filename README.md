DEBRA demo
================
Yevhen Akimov
3/1/2020

DEBRA - DESeq-based Barcode Representation Analysis
===================================================

DEBRA algorithm is developed for the identification of differentially represented clones (clones that differentially respond to the assay conditions) in clone tracing experiments. The clone tracing experiments are often associated with treatment-induced reduction in sample size which has been shown to affect the statistical properties of barcode read counts. Traditional algorithms (DESeq,DESeq2,edgeR) used for the inference of the differentially represented sequencing tags fail to account for these discrepancies and produce high number of false discoveries when the difference in size between control and condition samples is high. DEBRA utilizes the DESeq/DESeq2 functionality to reliably call differentially represented barcodes (DRBs) in clone tracing experiments.

Start by installing the package
-------------------------------

``` r
library(devtools)
install_github("YevhenAkimov/DEBRA")
```

Let's attach the DEBRA library and load example dataset
-------------------------------------------------------

``` r
library(DEBRA)
```

    ## Loading required package: locfit

    ## locfit 1.5-9.1    2013-03-22

``` r
data(bar_counts)
head(bar_counts)
```

    ##                                           null_660.1 null_660.2 null_330.1
    ## AAAAACAGACGCGTATGGGGCGGTCTAGCCCAATGAATGGA          0          3         13
    ## AAAACATGATCTCCGAGCGATGGCCAATGATAAAGGGAGGA          3          5          0
    ## AAAACATGATTGGAGAAGAATGGTCCGTTCTGTAGAATGGT         87        119        110
    ## AAAATGGCCCGTCCCCGAGACGGCCACAGAAGTCGAGAGGT          2          6          6
    ## AAACAGTAAGTGAATAGAGAAGGACCTCGGTGATCAGCGGA          6          4          5
    ## AAACCCTTAGGATCATTGGGGGGATGACGATTCAGGAAGGC        485        501        445
    ##                                           null_330.2 m_null_40.p35.1
    ## AAAAACAGACGCGTATGGGGCGGTCTAGCCCAATGAATGGA          1               2
    ## AAAACATGATCTCCGAGCGATGGCCAATGATAAAGGGAGGA          0               0
    ## AAAACATGATTGGAGAAGAATGGTCCGTTCTGTAGAATGGT         96              63
    ## AAAATGGCCCGTCCCCGAGACGGCCACAGAAGTCGAGAGGT          0               0
    ## AAACAGTAAGTGAATAGAGAAGGACCTCGGTGATCAGCGGA          5               0
    ## AAACCCTTAGGATCATTGGGGGGATGACGATTCAGGAAGGC        462             319
    ##                                           m_null_40.p35.2
    ## AAAAACAGACGCGTATGGGGCGGTCTAGCCCAATGAATGGA               0
    ## AAAACATGATCTCCGAGCGATGGCCAATGATAAAGGGAGGA               4
    ## AAAACATGATTGGAGAAGAATGGTCCGTTCTGTAGAATGGT              65
    ## AAAATGGCCCGTCCCCGAGACGGCCACAGAAGTCGAGAGGT               0
    ## AAACAGTAAGTGAATAGAGAAGGACCTCGGTGATCAGCGGA               4
    ## AAACCCTTAGGATCATTGGGGGGATGACGATTCAGGAAGGC             293

**Define control and condition sample names and run DEBRA script**
------------------------------------------------------------------

``` r
control_names = c("null_660.1","null_660.2","null_330.1","null_330.2" )
condition_names = c("m_null_40.p35.1","m_null_40.p35.2")

drb=DEBRA(bar_counts, control_names=control_names, condition_names=condition_names, method="DESeq2(Wald)")
```

    ## [1] KS tests completed
    ## [1] "Overlap estimation completed"

**Extracting results**
----------------------

Sorted results

``` r
results=resultsDRB(drb)
head(results)
```

    ##                                            baseMean log2FoldChange
    ## TGTGACACCTTATTCCGGGGAGGCAACATTCAAGGGGCGGA 3013.6896     -0.6348834
    ## CGTTACTCGGGGGATTCCGGGGGACATGCGTGAGGCGTGGA  970.3050     -0.7324366
    ## TTGCACACGATATCGCAAGGAGGAATTTGTTGTTGGGAGGT  950.7817      0.7207876
    ## CCTGGGCTCGGCGTCTGAGATGGTAAGTGATTCCAGGTGGC 1614.1302     -0.6200769
    ## ATTGTCGGCAGGAAGAAGCGTGGTGCCTGCCTTGGGAGGGT 1186.2894     -0.6587834
    ## TCTTTTTCGCGAAAGTAGGGCGGTGAACACTTTCCAGAGGC  199.3760     -1.4956945
    ##                                               lfcSE      stat       pvalue
    ## TGTGACACCTTATTCCGGGGAGGCAACATTCAAGGGGCGGA 0.1108983 -5.724918 1.034837e-08
    ## CGTTACTCGGGGGATTCCGGGGGACATGCGTGAGGCGTGGA 0.1425014 -5.139854 2.749521e-07
    ## TTGCACACGATATCGCAAGGAGGAATTTGTTGTTGGGAGGT 0.1402106  5.140750 2.736445e-07
    ## CCTGGGCTCGGCGTCTGAGATGGTAAGTGATTCCAGGTGGC 0.1232661 -5.030393 4.894754e-07
    ## ATTGTCGGCAGGAAGAAGCGTGGTGCCTGCCTTGGGAGGGT 0.1335535 -4.932731 8.108773e-07
    ## TCTTTTTCGCGAAAGTAGGGCGGTGAACACTTTCCAGAGGC 0.3018721 -4.954729 7.243127e-07
    ##                                                   padj
    ## TGTGACACCTTATTCCGGGGAGGCAACATTCAAGGGGCGGA 8.992732e-06
    ## CGTTACTCGGGGGATTCCGGGGGACATGCGTGAGGCGTGGA 7.964445e-05
    ## TTGCACACGATATCGCAAGGAGGAATTTGTTGTTGGGAGGT 7.964445e-05
    ## CCTGGGCTCGGCGTCTGAGATGGTAAGTGATTCCAGGTGGC 1.063385e-04
    ## ATTGTCGGCAGGAAGAAGCGTGGTGCCTGCCTTGGGAGGGT 1.174421e-04
    ## TCTTTTTCGCGAAAGTAGGGCGGTGAACACTTTCCAGAGGC 1.174421e-04

Sorted results with no independent filtering applied

``` r
results=resultsDRB(drb,filtered=F)
head(results)
```

    ##                                            baseMean log2FoldChange
    ## TGTGACACCTTATTCCGGGGAGGCAACATTCAAGGGGCGGA 3013.6896     -0.6348834
    ## CGTTACTCGGGGGATTCCGGGGGACATGCGTGAGGCGTGGA  970.3050     -0.7324366
    ## TTGCACACGATATCGCAAGGAGGAATTTGTTGTTGGGAGGT  950.7817      0.7207876
    ## CCTGGGCTCGGCGTCTGAGATGGTAAGTGATTCCAGGTGGC 1614.1302     -0.6200769
    ## ATTGTCGGCAGGAAGAAGCGTGGTGCCTGCCTTGGGAGGGT 1186.2894     -0.6587834
    ## TCTTTTTCGCGAAAGTAGGGCGGTGAACACTTTCCAGAGGC  199.3760     -1.4956945
    ##                                               lfcSE      stat       pvalue
    ## TGTGACACCTTATTCCGGGGAGGCAACATTCAAGGGGCGGA 0.1108983 -5.724918 1.034837e-08
    ## CGTTACTCGGGGGATTCCGGGGGACATGCGTGAGGCGTGGA 0.1425014 -5.139854 2.749521e-07
    ## TTGCACACGATATCGCAAGGAGGAATTTGTTGTTGGGAGGT 0.1402106  5.140750 2.736445e-07
    ## CCTGGGCTCGGCGTCTGAGATGGTAAGTGATTCCAGGTGGC 0.1232661 -5.030393 4.894754e-07
    ## ATTGTCGGCAGGAAGAAGCGTGGTGCCTGCCTTGGGAGGGT 0.1335535 -4.932731 8.108773e-07
    ## TCTTTTTCGCGAAAGTAGGGCGGTGAACACTTTCCAGAGGC 0.3018721 -4.954729 7.243127e-07
    ##                                                   padj
    ## TGTGACACCTTATTCCGGGGAGGCAACATTCAAGGGGCGGA 3.511201e-05
    ## CGTTACTCGGGGGATTCCGGGGGACATGCGTGAGGCGTGGA 3.109708e-04
    ## TTGCACACGATATCGCAAGGAGGAATTTGTTGTTGGGAGGT 3.109708e-04
    ## CCTGGGCTCGGCGTCTGAGATGGTAAGTGATTCCAGGTGGC 4.151975e-04
    ## ATTGTCGGCAGGAAGAAGCGTGGTGCCTGCCTTGGGAGGGT 4.585511e-04
    ## TCTTTTTCGCGAAAGTAGGGCGGTGAACACTTTCCAGAGGC 4.585511e-04

**Now we can display analytical plots**
---------------------------------------

MA-plot

``` r
MAplot(drb)
```

![](DEBRA_markdown_files/figure-markdown_github/plot-1.png)

``` r
MAplot(drb,FDR=0.1)
```

![](DEBRA_markdown_files/figure-markdown_github/plot-2.png)

MA-plot for non-filtered data

``` r
MAplot(drb, filtered =F)
```

![](DEBRA_markdown_files/figure-markdown_github/unnamed-chunk-6-1.png)

Visualize independet filtering results

``` r
filterPlot(drb)
```

![](DEBRA_markdown_files/figure-markdown_github/unnamed-chunk-7-1.png)

**Visualized data related to beta threshold estimation**
--------------------------------------------------------

Plot theoretical and empirical Kolmogorov-Smirnov (KS) statistic values

``` r
KS_plot(drb)
```

![](DEBRA_markdown_files/figure-markdown_github/unnamed-chunk-8-1.png)

Plot local overlap between Gamma-modelled empirical and theoretical Kolmogorov-Smirnov statistic values

``` r
overlapFit(drb)
```

![](DEBRA_markdown_files/figure-markdown_github/unnamed-chunk-9-1.png)

**Advanced DEBRA workflow**
---------------------------

We start by creating the DEBRADataSet. In this example we will use DESeq2(LRT) test.

``` r
drb_2=DEBRADataSet( counts =  bar_counts,
              control_names  = control_names,
              condition_names = condition_names,
              method=c("DESeq2(LRT)"),
              trended=T)
```

Next, we estimate the beta threshold with default parameters.

``` r
drb_2= estimateBeta(drb_2)
```

    ## [1] KS tests completed
    ## [1] "Overlap estimation completed"

``` r
KS_plot(drb_2)
```

![](DEBRA_markdown_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
overlapFit(drb_2)
```

![](DEBRA_markdown_files/figure-markdown_github/unnamed-chunk-11-2.png)

The estimation was successfull, but what if it is not? Let's change the "s" parameter to 200, so that there is not enough of the data points for a proper fitting of the overlap between theoretical and empirical KS test statistics.

``` r
drb_2= estimateBeta(drb_2,s = 200)
```

    ## [1] "found beta value, replacing it"
    ## [1] KS tests completed
    ## [1] "Overlap estimation completed"
    ## [1] "Failed fitting with LL.4() sigmoid; trying G.3() (see drc::drm for details)"
    ## [1] "Failed to fit mean-KS overlap with sigmoid; run KS_plot(drb) to assess the NB fit quality and manually set beta threshold. Setting beta to default_beta value. "

Let's see the KS test plot

``` r
KS_plot(drb_2)
```

![](DEBRA_markdown_files/figure-markdown_github/unnamed-chunk-13-1.png)

Now we can visually assess the mean value at which two KS test statstics start to ovrlap and set the beta value manually.

``` r
drb_2=DEBRADataSet( counts =  bar_counts,
              control_names  = control_names,
              condition_names = condition_names,
              method=c("DESeq2(LRT)"),
              trended=T,
              beta=15)
```

As we set the beta threshold manually, we can skip estimatBeta() command and run the testDRB() right away followed by independent filtering

``` r
drb_2=testDRB(drb_2)
```

    ## converting counts to integer mode

``` r
drb_2=independentFilteringDRB(drb_2)
```

Let's now call DRBs using original DESeq2(LRT) and compare these to the DEBRA results

``` r
require(ggplot2)
```

    ## Loading required package: ggplot2

``` r
drb_original=DEBRA(bar_counts, control_names=control_names, condition_names=condition_names, method="DESeq2(LRT)",modified = F)
```

    ## [1] " no beta value found, using beta = 0 "

``` r
res_orig=resultsDRB(drb_original)
res_debra=resultsDRB(drb_2)


x=rbind(cbind.data.frame(padj=res_orig$padj,method="DESeq2",baseMean=res_orig$baseMean),
        cbind.data.frame(padj=res_debra$padj,method="DEBRA",baseMean=res_debra$baseMean))

ggplot(x[x$padj<0.25,], aes(x=log10(baseMean),fill=method))+geom_histogram( position="dodge",bins=10)+theme_minimal()
```

![](DEBRA_markdown_files/figure-markdown_github/unnamed-chunk-16-1.png)

We observe the large number of DRBs comes from the low count region when assesed with the original DESeq2
