library(ggplot2)
library(reshape2)

dat = read.csv(
  "diagnostic_matrix.csv", header = TRUE, comment.char = "#", 
  row.names = "CLASSIFIER_DATA")
# rownames(dat) = paste0(rownames(dat), "_clf")
# colnames(dat) = paste0(colnames(dat), "_data")

dat = as.data.frame(as.table(as.matrix(dat)))
colnames(dat) = c("Classifiers", "Val.Data", "Accuracy")

p = ggplot(data = dat) + geom_tile(aes(x = Classifiers, y = Val.Data, fill = Accuracy)) + 
  scale_fill_gradientn(colors=c("skyblue","yellow","red")) + 
  xlab("Classifiers") + ylab("Validation Data") + 
  ggtitle(
    "Accuracies of Dead Organoid Classifiers", 
    subtitle = "Classifiers were trained on one cell line and applied to datasets of all cell lines") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot(p)
ggsave("diagnostic_matrix.png", width = 6, height = 6)
