summary.FstEnv <-
function(fstenv){
  num.model <- length(fstenv)
  models <- NULL
  TIC <- NULL
  R2 <- NULL
  coeffs <- list()
  for(i in 1:11){
    models <- c(models,fstenv[[i]]$model) 
    TIC <- c(TIC,fstenv[[i]]$TIC) 
    R2 <- c(R2,fstenv[[i]]$R2)
    coeffs[[i]] <- fstenv[[i]]$coefficients 
  }

  models.tab <- data.frame(models,TIC,R2)
  row.names(models.tab) <- names(coeffs) <- paste0("model",1:num.model)

  bestmodel <- rownames(models.tab)[order(models.tab[,"TIC"])[1]]
  bestmodel.coeff <- coeffs[[bestmodel]]
  bestmodel.model <- models.tab[bestmodel,]

  outlist <- list(
    models=models.tab,
    best.model=as.character(bestmodel.model[,"models"]),
    best.coefficients=bestmodel.coeff,
    best.TIC=bestmodel.model[,"TIC"],
    best.R2=bestmodel.model[,"R2"]
  )
  class(outlist) <- "summary.FstEnv"
  return(outlist)
}
