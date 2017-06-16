print.summary.FstEnv <-
function(sfe){
  print(sfe$models)
  message("\nBest model: ", sfe$best.model)
  message("Coefficients:")
  print(sfe$best.coefficients)
  message("TIC= ", round(sfe$best.TIC,4))
  message("R2= ", round(sfe$best.R2,4))
}
