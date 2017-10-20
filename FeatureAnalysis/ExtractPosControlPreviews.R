# This script retrieves the HTML previews for the positive controls as a visual verification
library(PROMISE)
library(hwriter)
configdir = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/configdir"
setwd("~/Thesis/Projects/PROMISE/FeatureAnalysis/")

# Set up controls
pos_ctrls = list(
  "Bortezomib", "Irinotecan / SN-38", "Volasertib", 
  "Methotrexate", "Staurosporine_500nM")
pos_ctrls_conc = list(
  c(0.2, 1), c(0.2, 1), c(0.2, 1), c(0.2, 1), c(0.0016, 0.008, 0.04, 0.2, 1))

# Load layout
lib = loadLibrary(library = "L08", configdir = configdir)
lib$Product.Name = as.character(lib$Product.Name)
lib$Well_ID_384 = as.character(lib$Well_ID_384)
lib$concentration = as.numeric(lib$concentration)

# Load wells
pos_ctrls_wells = c()
pos_ctrls_wells_id = c()
for(i in seq_along(pos_ctrls)) {
  ctrl = pos_ctrls[[i]]
  conc = pos_ctrls_conc[[i]]
  relevant_entries = lib[intersect(
    grep(ctrl, lib$Product.Name, ignore.case = TRUE), 
    which(lib$concentration %in% conc)),]
  pos_ctrls_wells = c(pos_ctrls_wells, relevant_entries$Well_ID_384)
  pos_ctrls_wells_id = c(pos_ctrls_wells_id, paste0(
    relevant_entries$Product.Name, " @ ", relevant_entries$concentration))
}

# Load plates
all_l08_plates = list.files(featuresdir, pattern = "^D.*L08$")

# Remove plates that were imaged in duplicate (9XX vs 0XX)
cell_lines = substr(all_l08_plates, 1, 7)
plates_to_use = c()
for(cl in unique(cell_lines)) {
  cl_plates = grep(cl, all_l08_plates, value = TRUE)
  plate_nums = substr(cl_plates, 9, 11)
  for(i in seq_along(cl_plates)) {
    if(substr(plate_nums[[i]], 1, 1) == "9") {
      plates_to_use = c(plates_to_use, TRUE)
      next
    }
    reimaged = paste0("9", substr(plate_nums[[i]], 2, 3))
    if(reimaged %in% plate_nums) {
      plates_to_use = c(plates_to_use, FALSE)
    } else {
      plates_to_use = c(plates_to_use, TRUE)
    }
  }
}
all_l08_plates = all_l08_plates[plates_to_use]

# Copy preview files into local directory and load them into matrices
preview_files = matrix(
  "", nrow = length(all_l08_plates), ncol = length(pos_ctrls_wells), 
  dimnames = list(all_l08_plates, pos_ctrls_wells_id))
small_preview_files = matrix(
  "", nrow = length(all_l08_plates), ncol = length(pos_ctrls_wells), 
  dimnames = list(all_l08_plates, pos_ctrls_wells_id))
for(r in seq_along(all_l08_plates)) {
  plate = all_l08_plates[[r]]
  for(c in seq_along(pos_ctrls_wells)) {
    well = pos_ctrls_wells[[c]]
    img_fn = file.path(htmldir, plate, sprintf("%s_%s_%s_1.jpeg", plate, substr(well, 1, 1), substr(well, 2, 3)))
    small_img_fn = file.path(htmldir, plate, sprintf("%s_%s_%s_2.jpeg", plate, substr(well, 1, 1), substr(well, 2, 3)))
    # file.copy(img_fn, "pos_ctrl_previews/")
    # file.copy(small_img_fn, "pos_ctrl_previews/")
    # img_fn = file.path("pos_ctrl_previews", basename(img_fn))
    # small_img_fn = file.path("pos_ctrl_previews", basename(small_img_fn))
    preview_files[r, c] = img_fn
    small_preview_files[r, c] = small_img_fn
  }
}

# Create html file
page = openPage(file.path("PosControlPreviews.html"), link.css="hwriter.css")
preview_files = hwriteImage(preview_files, table=FALSE)
row.names(preview_files) = all_l08_plates
colnames(preview_files) = paste0(pos_ctrls_wells_id, " (", pos_ctrls_wells, ")")
hwrite(preview_files, page=page, br=TRUE)
closePage(page, splash=FALSE)
