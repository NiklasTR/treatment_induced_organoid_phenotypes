---
title: "EDA organoid partition"
author: "Niklas Rindtorff"
output:
  html_document:
    keep_md: true
  pdf_document: default
params:
  data: "data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_sampled.Rds"
  echo: FALSE
---




Loading packages


```
## [1] "parameter input:"
```

```
## [1] "data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_sampled.Rds"
```

loading input data and annotation. Note that on the central cluster, with access to the complete data table, the definition of the input can easily be changed. For remote work, the subsampled dataset "umap_drugs_sampled.Rds" is the default choice.



# Partition inspection

We are able to observe 4 partitions in our data. 
After manual inspection, it becomes cleat that the two smallest partitions are mostly consisting of 



![](nr_morphology_partition_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

![](nr_morphology_partition_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


## drug overrepresentation

![](nr_morphology_partition_files/figure-html/unnamed-chunk-6-1.png)<!-- -->



![](nr_morphology_partition_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

## chi-square


I wonder wether certain batches or organoid lines are overrepresented in each section. 



![](nr_morphology_partition_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

## line overrepresentation




![](nr_morphology_partition_files/figure-html/unnamed-chunk-11-1.png)<!-- -->


I plot chisq residuals for each plate

I recognize no difference between reimaged plates (leading digit is "9", plates were reimaged due to errors during the first pass) and plates that were not reimaged.

![](nr_morphology_partition_files/figure-html/unnamed-chunk-12-1.png)<!-- -->


## multinomial regression

We run a multinomial regression using the *nnet* package.


```
## .
##      1      2      3      4 
## 283583  21320   3385   1228
```

```
## # weights:  8 (3 variable)
## initial  value 429080.285479 
## iter  10 value 103934.587764
## final  value 103929.549407 
## converged
```

```
## # weights:  56 (39 variable)
## initial  value 429080.285479 
## iter  10 value 93273.146824
## iter  20 value 89261.723946
## iter  30 value 87909.198169
## iter  40 value 87104.566970
## iter  50 value 86086.654631
## iter  60 value 86028.760916
## iter  70 value 86010.747848
## iter  80 value 86003.923852
## iter  90 value 86002.864912
## iter 100 value 85998.378411
## final  value 85998.378411 
## stopped after 100 iterations
```

```
## # weights:  316 (234 variable)
## initial  value 429080.285479 
## iter  10 value 82518.666121
## iter  20 value 80085.517701
## iter  30 value 79545.739683
## iter  40 value 79080.968936
## iter  50 value 78664.558187
## iter  60 value 78273.598036
## iter  70 value 77701.258495
## iter  80 value 76890.920329
## iter  90 value 76383.229317
## iter 100 value 76309.712980
## final  value 76309.712980 
## stopped after 100 iterations
```

```
## # weights:  36 (24 variable)
## initial  value 429080.285479 
## iter  10 value 96298.512754
## iter  20 value 92401.865374
## iter  30 value 91225.990894
## iter  40 value 90530.988437
## iter  50 value 90490.948091
## iter  60 value 90480.224419
## iter  70 value 90470.557343
## final  value 90470.537303 
## converged
```



![](nr_morphology_partition_files/figure-html/unnamed-chunk-14-1.png)<!-- -->



![](nr_morphology_partition_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

# observation
* after processing, organoids organize in 4 distinct phenotype partitions
* the distribution of organoids across these partitions is non-random
  * the screening plate influences the distribution of organoids across partitions, 3 plates show strong deviation from the expected distribution, both in a chi-square test **and** in a multinomial regression. The feature **plate** is more predictive than the organoid **line** or **screen_id**
    * D027T01P906L03
    * D020T01P906L03
    * D013T01P001L02
  * drug treatment influences the distribution of organoids across the partitions
    * DMSO control treatment are depleted in sector 2, sector 2 has previously been shown to contain dead and small organoids
    * SN38 and bortezomib are enriched in sector 2
  
# conclusion
* given the patterns above, I believe it is most likely we are seeing systematic errors (not dependent on the biological axes of line and drug) in plates:
  * D027T01P906L03
  * D020T01P906L03
  * D013T01P001L02
* in addition we are observing deviations that are drug dependent, that I consider signal

# next steps
* inspect the three plates manually
* 
    

# figure






```
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: macOS  10.16
## 
## Matrix products: default
## BLAS/LAPACK: /Users/rindtorf/github/promise/env/lib/R/lib/libRblas.dylib
## 
## locale:
## [1] C
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] nnet_7.3-12     gridExtra_2.3   ggridges_0.5.3  scico_1.2.0    
##  [5] princurve_2.1.4 cowplot_1.1.1   ggrastr_0.2.3   here_0.1       
##  [9] forcats_0.4.0   stringr_1.4.0   dplyr_0.8.0.1   purrr_0.3.2    
## [13] readr_1.3.1     tidyr_0.8.3     tibble_2.1.1    ggplot2_3.1.1  
## [17] tidyverse_1.2.1
## 
## loaded via a namespace (and not attached):
##  [1] beeswarm_0.3.1     tidyselect_0.2.5   xfun_0.6          
##  [4] haven_2.1.0        lattice_0.20-38    colorspace_1.4-1  
##  [7] generics_0.0.2     htmltools_0.3.6    yaml_2.2.0        
## [10] rlang_0.3.4        pillar_1.3.1       glue_1.3.1        
## [13] withr_2.1.2        RColorBrewer_1.1-2 modelr_0.1.4      
## [16] readxl_1.3.1       plyr_1.8.4         munsell_0.5.0     
## [19] gtable_0.3.0       cellranger_1.1.0   rvest_0.3.3       
## [22] codetools_0.2-16   evaluate_0.13      labeling_0.3      
## [25] knitr_1.22         vipor_0.4.5        broom_0.5.2       
## [28] Rcpp_1.0.1         scales_1.0.0       backports_1.1.4   
## [31] jsonlite_1.6       hms_0.4.2          digest_0.6.18     
## [34] stringi_1.4.3      grid_3.6.1         rprojroot_1.3-2   
## [37] cli_1.1.0          tools_3.6.1        magrittr_1.5      
## [40] lazyeval_0.2.2     crayon_1.3.4       pkgconfig_2.0.2   
## [43] pheatmap_1.0.12    xml2_1.2.0         ggbeeswarm_0.6.0  
## [46] lubridate_1.7.4    assertthat_0.2.1   rmarkdown_1.12    
## [49] httr_1.4.0         rstudioapi_0.10    R6_2.4.0          
## [52] nlme_3.1-139       compiler_3.6.1
```


