library(ggplot2)

setwd("/Users/jansauer/Thesis/Projects/PROMISE/FeatureAnalysis/python_package/dead_organoid_classifier")
fname = "roc_data_D004T01_on_D021T01.csv"
clf_cl = strsplit(fname, "_")[[1]][3]
data_cl = strsplit(strsplit(fname, "_")[[1]][5], ".", fixed = TRUE)[[1]][1]

# Set up data
dat = read.csv(fname, header = TRUE, row.names = 1, comment.char = "#")
dat$Type = "ROC"
dat = rbind(dat, data.frame(
  "FalsePosRate" = c(0, 1), "TruePosRate" = c(0, 1), "Type" = "Diag"))

# Load AUC
auc = readLines(fname)[2]
auc = trimws(strsplit(auc, ":")[[1]][2])

# Plot
ggplot(data = dat) + geom_line(
    aes(x = FalsePosRate, y = TruePosRate, color = Type, linetype = Type), 
    size = 1.5) + 
  xlab("False Positive Rate") + ylab("True Positive Rate") + 
  theme(legend.position = "None") + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  ggtitle(
    label = sprintf("ROC for %s Classifier Applied to %s Data", clf_cl, data_cl), 
    subtitle = sprintf("AUC = %0.5s", auc))

out_fname = gsub("\\.csv$", "\\.png", fname)
ggsave(filename = out_fname, width = 6, height = 6)
