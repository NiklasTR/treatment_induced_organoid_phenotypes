basedir = "/Volumes/collab-bernd-fischer/PROMISE/layouts"
layout_files = list.files(file.path(basedir, "processed_layout_files"), pattern = ".csv")

layouts = vector("list", length(layout_files))
names(layouts) = layout_files
for(layout_file in layout_files) {
  layout = read.csv2(file.path(basedir, "processed_layout_files", layout_file), stringsAsFactors = FALSE, row.names = 1)
  # Ensure that the absolute minimum of headers are present:
  # Product.Name;Well_ID_384;Row_ID_384;Col_ID_384;Library_ID
  if(!("Product.Name" %in% colnames(layout) & 
     "Well_ID_384" %in% colnames(layout) & 
     "Row_ID_384" %in% colnames(layout) & 
     "Col_ID_384" %in% colnames(layout) & 
     "Library_ID" %in% colnames(layout))) {
    stop(sprintf("Missing essential headers in file %s", layout_file))
  }
  layouts[[layout_file]] = layout
}

# # Combine layouts into one data frame and save it
# layout_headers = lapply(layouts, colnames)
# common_layout_headers = Reduce(intersect, layout_headers)
# layouts = lapply(layouts, function(x) x[,common_layout_headers])
# layouts = do.call(rbind, layouts)
# write.csv2(x = layouts, file = file.path(basedir, "layouts.dat"))

# Save layouts as individual files (for backwards compatibility of code)
for(layout in layouts) {
  lib_ids = unique(layout$Library_ID)
  for(lib_id in lib_ids) {
    liblayout = layout[layout$Library_ID == lib_id,]
    if(nrow(liblayout) != 384) stop(sprintf("Incomplete library for ID %s", lib_id))
    liblayout$Library_ID = sprintf("%02d", lib_id)
    libname = paste0("L", sprintf("%02d", lib_id), ".txt")
    write.csv2(x = liblayout, file = file.path(basedir, libname), row.names = FALSE)
  }
}