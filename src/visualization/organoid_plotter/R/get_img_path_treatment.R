anno_df <- read_rds(here::here("../data/compound_libraries/annotation.Rds"))
path_to_images_dir <- here::here("../data/profiling/htmldir")


get_img_path_treatment <- function(path_to_images, treatment_of_interest, anno){
  list.files(path_to_images) %>% 
    tibble(file = ., path = paste0(path_to_images, "/", file)) %>% 
    # removing files that I do not want to work with
    filter(! grepl(pattern = ".md5s", x = file)) %>% 
    filter(! grepl(pattern = ".html", x = file)) %>% 
    filter(! grepl(pattern = "copylog", x = file)) %>% 
    filter(! grepl(pattern = "D0", x = file)) %>% 
    mutate(library_id = substr(file, 13, 14) %>% as.numeric()) %>% 
    left_join(anno, by = "library_id") %>%
    drop_na() %>% 
    filter(product_name %in% treatment_of_interest) %>% 
    mutate(path_img = paste0(path, "/", file, "_", row_id_384, "_", 
                             stringr::str_pad(col_id_384, width = 2, side = "left", pad = "0"),
                             "_1.jpeg")) %>% 
    dplyr::select(file, path, path_img, product_name, concentration) %>% 
    mutate(line = substr(file, 1, 7)) %>%
    return()
}

get_img_path_well <- function(path_to_images, library, col_id, row_id, anno){
  list.files(path_to_images) %>% 
    tibble(file = ., path = paste0(path_to_images, "/", file)) %>% 
    # removing files that I do not want to work with
    filter(! grepl(pattern = ".md5s", x = file)) %>% 
    filter(! grepl(pattern = ".html", x = file)) %>% 
    filter(! grepl(pattern = "copylog", x = file)) %>% 
    filter(! grepl(pattern = "D0", x = file)) %>% 
    mutate(library_id = substr(file, 13, 14) %>% as.numeric()) %>% 
    left_join(anno, by = "library_id") %>%
    drop_na() %>% 
    filter(col_id_384 == col_id,
           row_id_384 == row_id,
           library_id == library) %>% 
    mutate(path_img = paste0(path, "/", file, "_", row_id_384, "_", 
                             stringr::str_pad(col_id_384, width = 2, side = "left", pad = "0"),
                             "_1.jpeg")) %>% 
    dplyr::select(file, path, path_img, product_name, concentration) %>% 
    mutate(line = substr(file, 1, 7)) %>%
    return()
}

get_img_path_exact <- function(path_to_images, plate, col_id, row_id, anno){
  list.files(path_to_images) %>% 
    tibble(file = ., path = paste0(path_to_images, "/", file)) %>% 
    # removing files that I do not want to work with
    filter(! grepl(pattern = ".md5s", x = file)) %>% 
    filter(! grepl(pattern = ".html", x = file)) %>% 
    filter(! grepl(pattern = "copylog", x = file)) %>% 
    filter(! grepl(pattern = "D0", x = file)) %>% 
    mutate(library_id = substr(file, 13, 14) %>% as.numeric()) %>% 
    left_join(anno, by = "library_id") %>%
    drop_na() %>% 
    filter(col_id_384 == col_id,
           row_id_384 == row_id,
           file == plate) %>% 
    mutate(path_img = paste0(path, "/", file, "_", row_id_384, "_", 
                             stringr::str_pad(col_id_384, width = 2, side = "left", pad = "0"),
                             "_1.jpeg")) %>% 
    dplyr::select(file, path, path_img, product_name, concentration) %>% 
    mutate(line = substr(file, 1, 7)) %>%
    return()
}

get_img_data <- function(path){
  set.seed(312)
  
  img_df <- path %>%
    group_by(concentration, product_name, line) %>% 
    sample_n(1) %>%
    filter(!(product_name == "DMSO" & concentration != 1)) %>%
    mutate(img_raw = purrr::map(path_img, ~ imager::load.image(.x))) %>%
    mutate(img = purrr::map(img_raw, ~ .x %>% as.data.frame(wide = "c") %>% 
                              mutate(rgb.val = rgb(c.1, c.2, c.3))))
  
  return(img_df)
}

plot_img <- function(data, treatment_of_interest){
  img_df_plot <- data %>% 
    unnest(img) 
  
  img_df_plot %>%  
    ggplot(aes(x,y)) + 
    geom_raster(aes(fill = rgb.val)) + 
    scale_fill_identity() + 
    scale_y_reverse() + 
    facet_grid(line ~ product_name + concentration) + 
    theme_bw() + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ggsave(file = paste0(treatment_of_interest[1], "_mco.pdf"), 
           height = 3*(img_df_plot$line %>% unique() %>% length()),
           width = 3*(img_df_plot %>% ungroup %>% dplyr::select(product_name, concentration) %>% 
                        distinct() %>% nrow()))
}