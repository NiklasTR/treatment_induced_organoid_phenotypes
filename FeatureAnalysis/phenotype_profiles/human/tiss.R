library(PROMISE)

calc_tiss_vector = function(
    drug, cell_line, out_fn,
    feature_dir = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/features", 
    config_dir = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/configdir") {
  
  source(file.path(config_dir, "watchdogConfig.R"))
  
  all_plates = list.files(feature_dir, pattern = paste0(cell_line, "P[0-9]{3}L08$"))
  
  # Use only the reimaged plates (9XX vs 0XX)
  plate_ids = substr(all_plates, 9, 11)
  use_plate = c()
  for(plate_id in plate_ids) {
    if(substr(plate_id, 0, 0) == "9") {
      use_plate = c(use_plate, TRUE)
      next
    }
    reimaged = paste0("9", substr(plate_id, 2, 3))
    if(reimaged %in% plate_ids) {
      use_plate = c(use_plate, FALSE)
    } else {
      use_plate = c(use_plate, TRUE)
    }
  }
  
  all_plates = all_plates[use_plate]
  all_plates = all_plates[order(substr(all_plates, 10, 14))]
  
  # Load image features  
  lib = loadLibrary("L08", config_dir)
  drug_wells = as.character(lib[lib$Product.Name == drug, "Well_ID_384"])
  dmso_wells = as.character(lib[lib$Product.Name == "DMSO", "Well_ID_384"])
  drug_features = list()
  drug_otypes = c()
  dmso_features = list()
  dmso_otypes = c()
  for(plate in all_plates) {
    for(well in drug_wells) {
      feat_fn = file.path(
        feature_dir, plate, "wells_normalized", 
        sprintf("%s_%s_%s_features_normalized.h5", 
                plate, substr(well, 1, 1), 
                substr(well, 2, 3)))
      feats = h5read(feat_fn, "features_organoids")
      colnames(feats) = h5read(feat_fn, "feature_names_organoids")
      well_id = paste0(plate, "_", substr(well, 1, 1), "_", substr(well, 2, 3))
      rownames(feats) = paste0(well_id, "__", seq_len(nrow(feats)))
      otype = h5read(feat_fn, "object_type_organoids")
      feats = feats[otype != "BLURRY", ]
      drug_otypes = c(drug_otypes, otype[otype != "BLURRY"])
      drug_features[[well_id]] = feats
    }
    for(well in dmso_wells) {
      feat_fn = file.path(
        feature_dir, plate, "wells_normalized", 
        sprintf("%s_%s_%s_features_normalized.h5", 
                plate, substr(well, 1, 1), 
                substr(well, 2, 3)))
      feats = h5read(feat_fn, "features_organoids")
      colnames(feats) = h5read(feat_fn, "feature_names_organoids")
      well_id = paste0(plate, "_", substr(well, 1, 1), "_", substr(well, 2, 3))
      rownames(feats) = paste0(well_id, "__", seq_len(nrow(feats)))
      otype = h5read(feat_fn, "object_type_organoids")
      feats = feats[otype != "BLURRY", ]
      dmso_otypes = c(dmso_otypes, otype[otype != "BLURRY"])
      dmso_features[[well_id]] = feats
    }
  }
  drug_features = do.call(rbind, drug_features)
  drug_metadata = data.frame(
    "CONCENTRATION" = lib[
      sapply(
        strsplit(rownames(drug_features), "_"), 
        function(x) paste0(x[2], x[3])), 
      "concentration"], 
    "TYPE" = drug_otypes)
  dmso_features = do.call(rbind, dmso_features)
  dmso_metadata = data.frame("TYPE" = dmso_otypes)
  
  # Split into replicate plates
  # Drug
  plates = substr(rownames(drug_features), 1, 14)
  plate_ids = substr(plates, 12, 14)
  plate_nums = substr(plates, 10, 11)
  replicates = rep(0, length(plate))
  for(plate_id in unique(plate_ids)) {
    layout_plates = plate_nums[plate_ids == plate_id]
    replicates[plate_ids == plate_id & plate_nums == min(layout_plates)] = 1
    replicates[plate_ids == plate_id & plate_nums == max(layout_plates)] = 2
  }
  drug_features_rep1 = drug_features[replicates == 1, , drop=FALSE]
  drug_metadata_rep1 = drug_metadata[replicates == 1, , drop=FALSE]
  drug_features_rep2 = drug_features[replicates == 2, , drop=FALSE]
  drug_metadata_rep2 = drug_metadata[replicates == 2, , drop=FALSE]
  
  # DMSO
  plates = substr(rownames(dmso_features), 1, 14)
  plate_ids = substr(plates, 12, 14)
  plate_nums = substr(plates, 10, 11)
  replicates = rep(0, length(plate))
  for(plate_id in unique(plate_ids)) {
    layout_plates = plate_nums[plate_ids == plate_id]
    replicates[plate_ids == plate_id & plate_nums == min(layout_plates)] = 1
    replicates[plate_ids == plate_id & plate_nums == max(layout_plates)] = 2
  }
  dmso_features_rep1 = dmso_features[replicates == 1, , drop=FALSE]
  dmso_metadata_rep1 = dmso_metadata[replicates == 1, , drop=FALSE]
  dmso_features_rep2 = dmso_features[replicates == 2, , drop=FALSE]
  dmso_metadata_rep2 = dmso_metadata[replicates == 2, , drop=FALSE]
  
  # Calculate KS vector
  drug_stats_rep1 = list()
  drug_stats_rep2 = list()
  drug_stats_pvals_rep1 = list()
  drug_stats_pvals_rep2 = list()
  for(concentration in sort(unique(drug_metadata$CONCENTRATION))) {
    conc_features_rep1 = drug_features_rep1[drug_metadata_rep1$CONCENTRATION == concentration, ]
    conc_metadata_rep1 = drug_metadata_rep1[drug_metadata_rep1$CONCENTRATION == concentration, ]
    conc_features_rep2 = drug_features_rep2[drug_metadata_rep2$CONCENTRATION == concentration, ]
    conc_metadata_rep2 = drug_metadata_rep2[drug_metadata_rep2$CONCENTRATION == concentration, ]
    
    feature_stats_rep1 = c()
    feature_stats_rep2 = c()
    feature_stats_pvals_rep1 = c()
    feature_stats_pvals_rep2 = c()
    for(feature in colnames(conc_features_rep1)) {
      drug_vals_rep1 = conc_features_rep1[, feature]
      drug_vals_rep2 = conc_features_rep2[, feature]
      dmso_vals_rep1 = dmso_features_rep1[, feature]
      dmso_vals_rep2 = dmso_features_rep2[, feature]
      
      # Perform both KS tests to determine the statistic.
      success = tryCatch({
        dmso_less_test_rep1 = ks.test(dmso_vals_rep1, drug_vals_rep1, alternative = "less")
        dmso_less_val_rep1 = dmso_less_test_rep1$statistic
        dmso_less_pval_rep1 = dmso_less_test_rep1$p.value
        
        dmso_less_test_rep2 = ks.test(dmso_vals_rep2, drug_vals_rep2, alternative = "less")
        dmso_less_val_rep2 = dmso_less_test_rep2$statistic
        dmso_less_pval_rep2 = dmso_less_test_rep2$p.value
        
        dmso_greater_test_rep1 = ks.test(dmso_vals_rep1, drug_vals_rep1, alternative = "greater")
        dmso_greater_val_rep1 = dmso_greater_test_rep1$statistic
        dmso_greater_pval_rep1 = dmso_greater_test_rep1$p.value
        
        dmso_greater_test_rep2 = ks.test(dmso_vals_rep2, drug_vals_rep2, alternative = "greater")
        dmso_greater_val_rep2 = dmso_greater_test_rep2$statistic
        dmso_greater_pval_rep2 = dmso_greater_test_rep2$p.value
        TRUE}, 
        error = function(e) FALSE)
      
      if(!success) {
        feature_stats_rep1[[feature]] = NA
        feature_stats_pvals_rep1[[feature]] = NA
        feature_stats_rep2[[feature]] = NA
        feature_stats_pvals_rep2[[feature]] = NA
        next
      }
      
      statistic_rep1 = NA
      statistic_rep2 = NA
      pval_rep1 = 1
      pval_rep2 = 1
      if(dmso_greater_val_rep1 > dmso_less_val_rep1) {
        statistic_rep1 = -dmso_greater_val_rep1
        pval_rep1 = dmso_greater_pval_rep1
      }
      if(dmso_greater_val_rep1 < dmso_less_val_rep1) {
        statistic_rep1 = dmso_less_val_rep1
        pval_rep1 = dmso_less_pval_rep1
      }
      feature_stats_rep1[[feature]] = statistic_rep1
      feature_stats_pvals_rep1[[feature]] = pval_rep1
      if(dmso_greater_val_rep2 > dmso_less_val_rep2) {
        statistic_rep2 = -dmso_greater_val_rep2
        pval_rep2 = dmso_greater_pval_rep2
      }
      if(dmso_greater_val_rep2 < dmso_less_val_rep2) {
        statistic_rep2 = dmso_less_val_rep2
        pval_rep2 = dmso_less_pval_rep2
      }
      feature_stats_rep2[[feature]] = statistic_rep2
      feature_stats_pvals_rep2[[feature]] = pval_rep2
    }
    
    drug_stats_rep1[[as.character(concentration)]] = feature_stats_rep1
    drug_stats_rep2[[as.character(concentration)]] = feature_stats_rep2
    drug_stats_pvals_rep1[[as.character(concentration)]] = feature_stats_pvals_rep1
    drug_stats_pvals_rep2[[as.character(concentration)]] = feature_stats_pvals_rep2
  }
  
  drug_stats_rep1 = unlist(drug_stats_rep1)
  drug_stats_rep2 = unlist(drug_stats_rep2)
  drug_stats_pvals_rep1 = unlist(drug_stats_pvals_rep1)
  drug_stats_pvals_rep2 = unlist(drug_stats_pvals_rep2)
  
  saveRDS(object = list(
    "rep1" = drug_stats_rep1, "rep2" = drug_stats_rep2, 
    "rep1_pval" = drug_stats_pvals_rep1, "rep2_pval" = drug_stats_pvals_rep2), 
    file = out_fn)
}