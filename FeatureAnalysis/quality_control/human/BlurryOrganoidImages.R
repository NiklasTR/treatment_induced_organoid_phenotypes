# This script generates the blurry organoid images
library(ggplot2)
library(pheatmap)

# Set up data
in_dir = "/Users/jansauer/Thesis/Projects/PROMISE/FeatureAnalysis/python_package/blurry_organoid_statistics"

all_plate_files = list.files(in_dir, pattern = "^M001|^D0")
worst_plates = list()
i = 1
for(plate in all_plate_files) {
  blurry_data = read.csv(file.path(in_dir, plate), header = TRUE, row.names = 1)
  focused = blurry_data$Focused_Organoids + blurry_data$Focused_Shrapnel
  blurry = blurry_data$Blurry_Organoids + blurry_data$Blurry_Shrapnel
  total_num = matrix(focused + blurry, nrow=16, ncol=24, byrow = TRUE)
  ratio = matrix(focused / (focused + blurry), nrow=16, ncol=24, byrow = TRUE)
  worst_plates[[i]] = data.frame(plate = strsplit(plate, "_")[[1]][1], min_num = min(total_num), min_ratio = min(ratio))
  i = i + 1
  # rownames(ratio) = LETTERS[1:16]
  # colnames(ratio) = as.character(1:24)
  # pheatmap(ratio, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = total_num)
}
worst_plates = do.call(rbind, worst_plates)
