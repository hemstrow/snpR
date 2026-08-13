factor_DT <- setClass(Class = 'factor_DT',
                      slots = c(data = "data.table",
                                cols = "character"), 
                      prototype = methods::prototype(data = "data.table",
                                                     cols = character(0)))

.as_factor_DT <- function(DT, factor_cols){
  if(any(duplicated(factor_cols))){
    stop("Duplicated column names not permitted.\n")
  }
  factor_cols <- factor_cols[which(factor_cols %in% colnames(DT))]
  
  DT[, `f&` := as.factor(do.call(paste, c(.SD, sep = "&"))), .SDcols = factor_cols]
  DT[,c(factor_cols) := NULL]
  setcolorder(DT, c("f&", setdiff(names(DT), "f&")))
  
  
  return(new("factor_DT", data = DT, cols = factor_cols))
}



setMethod("[", "factor_DT", function(x, i, j, ..., drop = TRUE){
  cat('hello')
  sub_i <- substitute(i)
  sub_j <- substitute(j)
})