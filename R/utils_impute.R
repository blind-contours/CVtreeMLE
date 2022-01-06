#' @title Impute NA values and create indicator variables
#' @param data Input data
impute_NA_vals <- function(data) {
  check_NA_in_cols <- sapply(data, function(x)any(is.na(x)))
  cols <- names(check_NA_in_cols[check_NA_in_cols == TRUE])

  append_these <-  data.frame( is.na(data[, cols] ))
  append_these  <- sapply(append_these, as.numeric)
  append_these <- as.data.frame(append_these)
  names(append_these) <- paste(cols, "NA", sep = "_")
  data = cbind(data, append_these)
  data[,cols]  <- sapply(data[,cols], as.numeric)

  for(i in 1:length(cols)){
    col <- cols[i]
    data[is.na(data[,col]), col] <- mean(data[,col], na.rm = TRUE)
  }

  return(list("data" = data, "impute cols" = names(append_these)))
}
