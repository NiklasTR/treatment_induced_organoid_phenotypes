---
title: "Organoid Unsupervised Exploration"
author: Niklas Rindtorff
output:
  html_document:
    keep_md: true
  pdf_document: default
params:
  data: "data/processed/PhenotypeSpectrum/filtered_radical/umap_absolute_all_drugs_sampled.Rds"
  remote: FALSE
  cache: FALSE
---



Loading packages


```r
library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(readr)
library(here)
library(ggrastr)
library(cowplot)
library(princurve)
library(scico)
library(ggridges)

# parameter
print("parameter input:")
```

```
## [1] "parameter input:"
```

```r
print(params$data)
```

```
## [1] "data/processed/PhenotypeSpectrum/filtered_radical/umap_absolute_all_drugs_sampled.Rds"
```

loading input data and annotation. Note that on the central cluster, with access to the complete data table, the definition of the input can easily be changed. For remote work, the subsampled dataset "umap_drugs_sampled.Rds" is the default choice.


```r
# I wish I could solve my path problems with the here() package, but experienced unreliable behavior 
# PATH = "/dkfz/groups/shared/OE0049/B110-Isilon2/promise/"
PATH = paste0(here::here(), "/")

#umap_df <- read_rds(paste0(PATH, "data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_tidy.Rds"))
#umap_df <- read_rds(paste0(PATH, "data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_sampled.Rds"))
umap_df <- read_rds(here::here(params$data))

organoid_morphology <- read_delim(here::here("references/imaging/visual_classification_organoids.csv"), ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  dplyr::select(line = organoid, morphology = visual_inspection_v2)
```

# Partition inspection

We are able to observe 4 partitions in our data. 
After manual inspection, it becomes cleat that the two smallest partitions are mostly consisting of 


```r
umap_df %>% 
  ggplot(aes(v1, v2, color = factor(partition))) + 
  geom_point_rast(alpha = 0.5, size = 0.35) + 
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme_cowplot() +
  labs(x = "UMAP 1",
       y = "UMAP 2",
       color = "partition") + 
  theme(legend.position = "bottom") + 
    coord_fixed()
```

![](organoid_unsupervised_exploration_files/figure-html/unnamed-chunk-3-1.png)<!-- -->


```r
umap_df %>% 
  dplyr::count(partition) %>% 
  mutate(ratio = n/sum(n)) %>% 
  arrange(desc(ratio))
```

```
## # A tibble: 2 x 3
##   partition      n  ratio
##   <fct>      <int>  <dbl>
## 1 1         213667 0.914 
## 2 2          19979 0.0855
```

I remove 2 partitions from all main figures for ease of reading. Below, it is easy to toggle the removal of partitions on and off to make sure this filtering step is robust


```r
gg_cluster <- umap_df %>%
  filter(partition %in% c(1,2)) %>%
  ggplot(aes(v1, v2, color = factor(cluster))) + 
  ggrastr::geom_point_rast(alpha = 0.5, size = 0.35) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(12, "Set3"), "#fb9a99")) +
  cowplot::theme_cowplot() +
  labs(x = "UMAP 1",
       y = "UMAP 2",
       color = "partition") + 
  theme(legend.position = "bottom") + 
    coord_fixed()

gg_cluster + ggsave(here("reports/figures/gg_cluster.pdf"), width = 4, height = 4)
```

![](organoid_unsupervised_exploration_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


# Organoid Size Distributions

I plot a size-distribution. 


```r
gg_size_dist <- umap_df %>% 
  filter(partition %in% c(1,2)) %>%
  ggplot(aes(size)) + 
  geom_histogram() + 
  theme_cowplot()

gg_size_dist_log <- gg_size_dist + 
  scale_x_log10() 

gg_size_dist_log
```

![](organoid_unsupervised_exploration_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

I add the eCDF.


```r
df <- umap_df %>% filter(partition %in% c(1,2))

gg_ecdf <- ggplot(df %>% filter(drug == "DMSO")) +
  stat_ecdf(aes(x = size, group = line), 
              geom = "step", size = 1) +
  #scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  labs(y = "f(size)",
       x = "organoid size [pixels]") + 
  theme_cowplot()

gg_ecdf
```

![](organoid_unsupervised_exploration_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

For more details about distributions, please refer to *reports/Phenotypespectrum/'xyz'_dist.pdf*.


```r
line_param <- umap_df %>% filter(partition %in% c(1,2)) %>% 
  nest(-line, -replicate) %>% 
  mutate(fit = map(data, ~ fitdistrplus::fitdist(.x$size, "lnorm")),
         param = map(fit, ~ .x$estimate %>% broom::tidy()))

df <- line_param %>% unnest(param) %>% 
  filter(names == "meanlog") %>% 
  group_by(line) %>% 
  mutate(mean_meanlog = mean(x)) %>% 
  arrange(mean_meanlog) %>% 
  ungroup() %>%
  mutate(line = factor(line, levels = .$line %>% unique()))

organoid_size_factor <- df$line %>% levels()

df <- df %>% 
    dplyr::select(line, replicate, x) %>%
    # tidyr::pivot_wider(names_from = replicate, 
    #             values_from = x)
    tidyr::spread(key = replicate, value = x)

r_size = df %>% ungroup() %>% dplyr::select(-line) %>% as.matrix %>% cor() %>% min()

gg_size_replicate <- df %>% 
  ggplot(aes(`1`, `2`)) +
  geom_smooth(method = "lm", se = FALSE, color = "grey") +
  geom_point() + 
  theme_cowplot() + 
  labs(x = "Replicate 1, mu [ln size]",
       y = "Replicate 2, mu [ln size]",
       caption = paste0("2 replicates, r= ", round(r_size, 2))) + 
  coord_fixed(ratio = 1) 
  #geom_abline(slope = 1, color = "grey")
  

gg_size_replicate
```

![](organoid_unsupervised_exploration_files/figure-html/unnamed-chunk-8-1.png)<!-- -->


```r
organoid_size_factor_09 <- umap_df %>% filter(partition %in% c(1,2)) %>% group_by(line) %>% 
  summarise(x = quantile(size_log, 0.9)) %>% 
  #summarise(x = mean(size_log)) %>% 
  arrange(x) %>% .$line


gg_size_dist_morph_ridge <- umap_df %>% filter(partition %in% c(1,2)) %>% filter(drug == "DMSO") %>% 
  mutate(line = factor(line, levels = organoid_size_factor_09)) %>% 
  ggplot() +
  geom_density_ridges_gradient(aes(y = line, x = size_log, fill = stat(x)), scale = 1) +
  #geom_density(aes(x = size_log, group = replicate, color = morphological_class)) + 
  #facet_wrap(~ line) + 
  scale_fill_viridis_c() +
  labs(caption = "DMSO treated organoids",
       x = "ln(size)",
       fill = "size") + 
  theme(legend.position = "bottom") +
  theme_cowplot() 

gg_size_dist_morph_ridge
```

![](organoid_unsupervised_exploration_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


```r
umap_size <- function(umap){
  umap %>%
  #filter(Size < 1000) %>%
  ggplot(aes(v1, v2, color = size_log)) + 
  geom_point_rast(alpha = 0.5, size = 0.35) + 
  scale_color_viridis_c() +
  theme_cowplot() +
  labs(x = "UMAP 1",
       y = "UMAP 2",
       color = "ln(size)") + 
  theme(legend.position = "bottom") + 
    coord_fixed()
}

gg_size <- umap_size(umap_df %>% filter(partition %in% c(1,2)))
```



# Organoid Size during Drug Treatment



```r
drug_size <- umap_df %>% filter(partition %in% c(1,2)) %>% filter(drug == "DMSO" | drug == "Paclitaxel") %>% 
  mutate(concentration = ifelse(drug == "DMSO", 0, concentration)) %>%
  #filter(morphological_class == "disorganized") %>% 
  #filter(morphological_class != "other") %>% 
  mutate(concentration = factor(concentration, levels = c("0", "0.0016", "0.008", "0.04", "0.2", "1.0")))

ggdrug_size <- ggplot(drug_size) +
  geom_density(aes(x = log(size), group = concentration, color = concentration), size = 1.5) + 
  scico::scale_color_scico_d() +
  theme_cowplot() + 
  #scale_x_continuous(limits = c(0, 15000)) +
  theme(legend.position = "bottom") + 
  labs(color = "Paclitaxel Concentration Factor",
       title = "Organoid size distribution",
       x = "ln(size)")

ggdrug_size 
```

![](organoid_unsupervised_exploration_files/figure-html/unnamed-chunk-11-1.png)<!-- -->


```r
drug_count <- drug_size %>%
  dplyr::count(concentration, line, replicate, well)

ggdrug_count <- ggplot(drug_count) +
  geom_density(aes(x = log(n), group = concentration, color = concentration), size = 1.5) + 
  scico::scale_color_scico_d() +
  theme_cowplot() + 
  theme(legend.position = "bottom") + 
  labs(color = "Paclitaxel Concentration Factor",
       title = "Organoid count distribution",
       x = "ln(n)")

ggdrug_count
```

![](organoid_unsupervised_exploration_files/figure-html/unnamed-chunk-12-1.png)<!-- -->



```r
set.seed(234)
loi = c("D022T01", "D046T01")


df <- umap_df  %>%
  filter(partition %in% c(1,2)) %>%
  filter(drug == "DMSO" | drug == "Paclitaxel") %>% 
  mutate(concentration = ifelse(drug == "DMSO", 0, concentration)) %>% 
  filter(line == loi)

gg_drug <- umap_df %>% filter(partition %in% c(1,2)) %>%
  dplyr::select(-line, -concentration) %>%
  ggplot(aes(v1, v2)) + 
  geom_point_rast(alpha = 1, size = 0.35, color = "#f1f1f1") + 
  geom_point_rast(data = df  %>%
    group_by(concentration) %>% 
    sample_n(1000, replace = TRUE),
  aes(color = concentration),alpha = 1, size = 1.5, shape=16) + 
  #facet_wrap( ~ concentration, ncol = 1) + 
  #scale_color_brewer(type = "seq", palette = "YlOrRd") + 
  #geom_density2d(color = "black") + 
  theme_classic() +
  labs(x = "UMAP 1",
       y = "UMAP 2",
       caption = paste0(paste(loi, collapse=" "), ", Paclitaxel"),
       color = "Concentration Factor") + 
  scico::scale_color_scico_d() + 
  facet_wrap(~ line, ncol = 2) +
  theme(legend.position = "bottom") +
  #theme_cowplot(font_size = 8) + 
  theme(legend.position = "bottom")
 
gg_drug
```

![](organoid_unsupervised_exploration_files/figure-html/unnamed-chunk-13-1.png)<!-- -->




```r
loi <- c("D022T01", "D055T01")

drug_size_param <- drug_size %>% 
  nest(-concentration, -line) %>% 
  mutate(fit = map(data, ~ fitdistrplus::fitdist(.x$size, "lnorm")),
         param = map(fit, ~ .x$estimate %>% broom::tidy()))

df <- drug_size_param %>% unnest(param) %>% 
  filter(names == "meanlog") %>%
  mutate(concentration = factor(concentration, levels = c("0", "0.0016", "0.008", "0.04", "0.2", "1.0"))) 

gg_size_drug <- df %>%
  filter(line %in% loi) %>%
  #mutate(concentration = as.numeric(as.character(concentration))) %>%
  ggplot(aes(concentration, x)) + 
    geom_point(color = "grey") + 
  geom_line(data = df %>% dplyr::rename(line_h = line) , aes(group = line_h), color = "grey") +
  geom_point(data = df %>% dplyr::rename(line_h = line), color = "grey") +
  geom_point(color = "black") +
  geom_line(aes(group = line), color = "black") +
  labs(y = 'mu ln(size)' ,
       x = "Concentration Factor",
       caption = paste0(paste(loi, collapse=" "), ", Paclitaxel")) +
  facet_wrap(~ line)+
  theme_cowplot()

gg_size_drug
```

![](organoid_unsupervised_exploration_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

# Organoid Heterogeneity

I plot 2 organoid lines treated with DMSO control


```r
set.seed(234)


df <- umap_df  %>% filter(partition %in% c(1,2)) %>%
  mutate(cystic = if_else(line == "D013T01" &  well == "D24" & plate == "D013T01P001L02", TRUE, FALSE)) %>%
  mutate(compact = if_else(line == "D055T01" &  well == "D24" & plate == "D055T01P007L02", TRUE, FALSE))

gg_cys_comp <- df %>% 
  sample_frac(0.01) %>%
  ggplot(aes(v1, v2, color = size_log)) + 
  #scale_color_brewer(type = "qual", palette = 2) +
  geom_point_rast(alpha = 0.1, size = 0.35) + 
  geom_point_rast(color = "#F4B400", alpha = 1, size = 0.5, data = df %>% filter(cystic == TRUE)) +
  geom_point_rast(color = "#DB4437", alpha = 1, size = 0.5, data = df %>% filter(compact == TRUE)) + 
  scale_color_viridis_c() +
  labs(color = "size",
       caption = "yellow: D013T01, red: D055T01",
       x = "UMAP 1",
       y = "UMAP 2") +
  theme_cowplot() + 
  coord_fixed()

gg_cys_comp
```

![](organoid_unsupervised_exploration_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

In general, DMSO treated organoid lines cover the same latent space than drug treated organoids. This is likely influenced by the large number of untreated organoids in the dataset.


```r
set.seed(123)

df <- umap_df %>% filter(partition %in% c(1,2))


gg_size_supp <- df %>%
  mutate(drug = if_else(drug == "DMSO", "DMSO", "other")) %>%
  ggplot(aes(v1, v2, color = size_log)) + 
  geom_point_rast(alpha = 0.5, size = 0.35) + 
  scale_color_viridis_c() +
  theme_cowplot() +
  labs(x = "UMAP 1",
       y = "UMAP 2") + 
  theme(legend.position = "bottom") +
  facet_wrap(~ drug)

gg_size_supp + ggsave(paste0(PATH, "reports/figures/gg_size_all.pdf"))
```

![](organoid_unsupervised_exploration_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

# Organoid line differences

I create a single plot showing the two extreme organoid lines and their distribution within the embedding. 


```r
set.seed(123)

loi <- c("D046T01", "D007T01",  "D027T01", "D018T01") #c("D055T01", "D007T01",  "D021T01", "D019T01", "D027T01")
#loi <- umap_df$line %>% unique()

df <- umap_df %>%
  filter(drug == "DMSO") %>% 
  filter(partition %in% c(1,2))


gg_line <- df %>% dplyr::select(-line) %>%
  ggplot(aes(v1, v2)) + 
  geom_point_rast(alpha = 1, size = 0.35, color = "#f1f1f1") + 
  geom_point_rast(data = umap_df %>%
                    filter(drug == "DMSO") %>% 
                   # filter(line %in% c("D021T01")) %>%
    filter(line %in% loi) %>% 
    mutate(line = factor(line, levels = loi)), #%>% 
      #sample_frac(0.1),
    #mutate(line = factor(line, levels = c("D021T01"))),
  aes(color = line),alpha = .4, size = 0.35, shape=16) + 
  facet_wrap( ~ line, ncol =2) +
  scale_color_brewer(type = "qual", palette = "Set2") +
  #scale_color_manual(values = c(c("#D80D12", "#461C01", "#9a4c91", "#70BE6F", "#24345E"))) + 	
  #geom_density2d(color = "black") + 
  theme_classic() +
  labs(x = "UMAP 1",
       y = "UMAP 2")+
       #caption = "control treated organoids") + 
  theme_cowplot(font_size = 8) + 
  theme(legend.position = "nothing")

gg_line + ggsave(paste0(PATH, "reports/figures/gg_size_all.pdf"), width = 4, height = 4)
```

![](organoid_unsupervised_exploration_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

I am focusing on cystic vs solid organoid lines


```r
#UMAP Cystic (Lines 18, 13, 27, 30) vs. Solid (others) treated with DMSO, for Figure 1 / matching expression analysis done for cystic vs. rest

set.seed(123)

cystic_l <- organoid_morphology %>% filter(morphology == "cystic") %>%.$line %>% paste0(., "01")
dense_l <- organoid_morphology %>% filter(morphology == "solid") %>%.$line %>% paste0(., "01")

df <- umap_df %>%
  filter(drug == "DMSO") %>% 
  filter(partition %in% c(1,2)) %>%
  mutate(morphology = case_when(line %in% cystic_l ~ "cystic",
                                line %in% dense_l ~ "solid",
                                 TRUE ~ "other"))

gg_cystic <- umap_df %>% 
  ggplot(aes(v1, v2)) + 
  geom_point_rast(alpha = 1, size = 0.35, color = "#f1f1f1") + 
  
  # geom_point_rast(data = df %>%
  #   filter(morphology == "cystic") %>% 
  #     sample_frac(0.05),
  # aes(color = morphology),alpha = .4, size = 0.35, shape=16) + 
  
  geom_density_2d(data = df %>% #geom_density_2d_filled
    filter(morphology != "other"), # %>% 
    #  sample_frac(0.05)
    aes(fill = morphology), size = 1.5) +

  #scale_color_brewer(type = "qual", palette = "Set2") +
  #scale_fill_manual(values = c("#0571b0", "#ca0020")) + 	
  scale_color_manual(values = c("#0571b0", "#ca0020")) + 	
  #geom_density2d(color = "black") + 
  theme_classic() +
  labs(x = "UMAP 1",
       y = "UMAP 2")+
       #caption = "control treated organoids") + 
  theme_cowplot(font_size = 8) + 
  #theme(legend.position = "nothing")  + 
  coord_fixed() + 
  scale_x_reverse()


gg_cystic
```

![](organoid_unsupervised_exploration_files/figure-html/unnamed-chunk-18-1.png)<!-- -->


# Plot Export


```r
#gg_size_dist + ggsave(paste0(PATH, "reports/figures/gg_size_dist.pdf"))
#gg_size_dist_log +  ggsave(paste0(PATH, "reports/figures/gg_size_dist_log.pdf"))
#ggdrug_size + ggsave(paste0(PATH, "reports/figures/gg_drug_dist_size_log.pdf"))
#gg_ecdf + ggsave(paste0(PATH, "reports/figures/gg_size_dist_ecdf.pdf"))
#gg_size_dist_morph_ridge +  ggsave(paste0(PATH, "reports/figures/gg_size_dist_morph_ridge.pdf"), width = 4, height = 4)

# ggdrug_count + ggsave(paste0(PATH, "reports/figures/gg_drug_dist_n_log.pdf"))
# gg_drug + ggsave(paste0(PATH, "reports/figures/gg_drug.pdf"), width = 8, height = 4)
# gg_size_drug + ggsave(paste0(PATH, "reports/figures/gg_trametinib_size_dose.pdf"), width = 3.65, height = 3.65)
gg_cystic + ggsave(paste0(PATH, "reports/figures/gg_cystic.pdf"), width = 4, height = 4)
```

![](organoid_unsupervised_exploration_files/figure-html/unnamed-chunk-19-1.png)<!-- -->


```r
plot_grid(plot_grid(gg_size_dist_morph_ridge, gg_size_replicate, labels = c('A', 'B'), label_size = 12, ncol = 2),
          gg_size,
          plot_grid(gg_line, gg_cystic, labels = c('D', 'E'), label_size = 12, ncol = 2),
          labels = c('', 'C', ''), label_size = 12, ncol = 1) +
  ggsave(paste0(PATH, "reports/panels/panel_size_dist.pdf"), width = 8, height = 16)

plot_grid(plot_grid(ggdrug_size, ggdrug_count, labels = c('A', 'B'), label_size = 12),
          gg_drug,
          gg_size_drug,
          labels = c('', 'C', 'D'), label_size = 12, ncol = 1) +
  ggsave(paste0(PATH, "reports/panels/panel_size_drug.pdf"), width = 8, height = 12)

gg_size_supp
```



```r
sessionInfo()
```

```
## R version 4.0.0 (2020-04-24)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.2 LTS
## 
## Matrix products: default
## BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] ggridges_0.5.2  scico_1.1.0     princurve_2.1.4 cowplot_1.0.0  
##  [5] ggrastr_0.2.3   here_0.1        readr_1.3.1     purrr_0.3.4    
##  [9] magrittr_1.5    tidyr_1.1.0     dplyr_1.0.0     ggplot2_3.3.1  
## 
## loaded via a namespace (and not attached):
##  [1] beeswarm_0.2.3     tidyselect_1.1.0   xfun_0.14          lattice_0.20-41   
##  [5] splines_4.0.0      colorspace_1.4-1   vctrs_0.3.1        generics_0.0.2    
##  [9] viridisLite_0.3.0  htmltools_0.4.0    mgcv_1.8-31        yaml_2.2.1        
## [13] utf8_1.1.4         survival_3.1-12    rlang_0.4.6        isoband_0.2.1     
## [17] pillar_1.4.4       glue_1.4.1         withr_2.2.0        fitdistrplus_1.1-1
## [21] RColorBrewer_1.1-2 lifecycle_0.2.0    plyr_1.8.6         stringr_1.4.0     
## [25] munsell_0.5.0      gtable_0.3.0       evaluate_0.14      labeling_0.3      
## [29] knitr_1.28         Cairo_1.5-12       vipor_0.4.5        fansi_0.4.1       
## [33] broom_0.5.6        Rcpp_1.0.4.6       scales_1.1.1       backports_1.1.7   
## [37] farver_2.0.3       hms_0.5.3          digest_0.6.25      stringi_1.4.6     
## [41] grid_4.0.0         rprojroot_1.3-2    cli_2.0.2          tools_4.0.0       
## [45] tibble_3.0.1       crayon_1.3.4       pkgconfig_2.0.3    MASS_7.3-51.5     
## [49] Matrix_1.2-18      ellipsis_0.3.1     ggbeeswarm_0.6.0   assertthat_0.2.1  
## [53] rmarkdown_2.2      R6_2.4.1           nlme_3.1-147       compiler_4.0.0
```

