
#' @title Get Replicate IDs for Plates
#'
#' @description Get hard-coded replicate IDs for all the plates.
#'
#' @usage get_replicate_ids_for_plates()
#'
#' @return A named vector containing the replicate IDs
#'
#' @author Jan Sauer
#'
#' @examples print(get_replicate_ids_for_plates)
#' @export
get_replicate_ids_for_plates = function(plates = NULL) {
  rep_ids = c(
    "M001A03P003L02"= 1, "M001A03P004L03"= 1, "M001A03P005L04"= 1,
    "M001A03P006L05"= 1, "M001A03P007L06"= 1, "M001A03P008L07"= 1,
    "M001A03P009L02"= 2, "M001A03P010L03"= 2, "M001A03P011L04"= 2,
    "M001A03P012L05"= 2, "M001A03P013L06"= 2, "M001A03P014L07"= 2,
    "M001B04P003L02"= 1, "M001B04P004L03"= 1, "M001B04P005L04"= 1,
    "M001B04P006L05"= 1, "M001B04P007L06"= 1, "M001B04P008L07"= 1,
    "M001B04P009L02"= 2, "M001B04P010L03"= 2, "M001B04P011L04"= 2,
    "M001B04P012L05"= 2, "M001B04P013L06"= 2, "M001B04P014L07"= 2,
    "M001K02P003L02"= 1, "M001K02P004L03"= 1, "M001K02P005L04"= 1,
    "M001K02P006L05"= 1, "M001K02P007L06"= 1, "M001K02P008L07"= 1,
    "M001K02P009L02"= 2, "M001K02P010L03"= 2, "M001K02P011L04"= 2,
    "M001K02P012L05"= 2, "M001K02P013L06"= 2, "M001K02P014L07"= 2,
    "M001W01P003L02"= 1, "M001W01P004L03"= 1, "M001W01P005L04"= 1,
    "M001W01P006L05"= 1, "M001W01P007L06"= 1, "M001W01P908L07"= 1,
    "M001W01P009L02"= 2, "M001W01P010L03"= 2, "M001W01P011L04"= 2,
    "M001W01P012L05"= 2, "M001W01P013L06"= 2, "M001W01P014L07"= 2)
  if(is.null(plates)) {
    return(rep_ids)
  } else {
    return(rep_ids[plates])
  }
}

#' @title Get All Lines
#'
#' @description Get the names of all lines
#'
#' @usage get_lines()
#'
#' @return A vector of line IDs
#'
#' @author Jan Sauer
#'
#' @examples print(get_lines)
#' @export
get_lines = function() {
  return(unique(substr(names(get_replicate_ids_for_plates()), 1, 7)))
}

#' @title Get All Plates
#'
#' @description Get the names of all plates
#'
#' @usage get_plates()
#'
#' @return A vector of plate IDs
#'
#' @author Jan Sauer
#'
#' @examples print(get_plates)
#' @export
get_plates = function() {
  return(unique(names(get_replicate_ids_for_plates())))
}
