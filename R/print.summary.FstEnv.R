print.summary.FstEnv <-
function(x, ...){
  LFx <- intToUtf8(0x0A)
  print(x$models)
  message(LFx, "Best model: ", x$best.model)
  message("Coefficients:")
  print(x$best.coefficients)
  message("TIC= ", round(x$best.TIC,4))
  message("R2= ", round(x$best.R2,4))
}
