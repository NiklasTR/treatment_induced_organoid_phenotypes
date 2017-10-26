selectByStability <- function(
  features, metadata, preselect, Rdim = 40, verbose = TRUE) {
  
  # Set up features for stability selection
  Drep1 = as.matrix(features[metadata$REPLICATE == 1, ])
  Drep2 = as.matrix(features[metadata$REPLICATE == 2, ])
  D = array(
    NA_real_, dim = c(nrow(Drep1), 2, ncol(Drep1)), 
    dimnames = list(NULL, NULL, colnames(features)))
  D[,1,] = Drep1
  D[,2,] = Drep2
  
  correlation = rep(NA_real_, Rdim)
  correlationAll = list()
  ratioPositive = rep(NA_real_, Rdim)
  selected = rep("", Rdim)
  
  DSel = D[,1,]
  DSel[] = NA_real_
  colnames(DSel) = NULL
  for (k in 1:Rdim) {
    if (k > length(preselect)) {
      for (i in length(preselect):dim(D)[3]) {
        if(verbose) cat("k=",k," i=",i,"\r")
        model = lm(as.vector(D[,1,i]) ~ DSel[,1:(k-1),drop=FALSE]+0)
        D[,1,i] = model$residuals
        model = lm(as.vector(D[,2,i]) ~ DSel[,1:(k-1),drop=FALSE]+0)
        D[,2,i] = model$residuals
      }
    }
    C = apply(D, 3, function(x) {
      a = x[,1]
      b = x[,2]
      I = which(is.finite(a) & is.finite(b))
      cor(a[I],b[I]) 
    })
    if (k > length(preselect)) {
      I = names(C)[which.max(C)]
    } else {
      I = preselect[k]
    }
    correlation[k] = C[I]
    ratioPositive[k] = sum(C > 0, na.rm=TRUE) / length(C)
    selected[k] = I
    if(verbose) cat(
      "k=",k," selected = ",selected[k], " cor = ", correlation[k], 
      " r = ", ratioPositive[k], "\n")
    correlationAll[[k]] = C
    DSel[,k] = apply(D[,,I,drop=FALSE],1,mean,na.rm=TRUE)
    D = D[,,dimnames(D)[[3]] != I,drop=FALSE]
    if(dim(D)[3] < length(preselect)) break
  }
  
  res = list(selected = selected, correlation = correlation, 
             ratioPositive = ratioPositive, correlationAll = correlationAll)
  return(res)
}