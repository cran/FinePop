calc.PopSt <-
function(popdata, method){
  numPop <- popdata$npops
  numMarker <- popdata$nloci
  numInd <- popdata$pop_sizes
  numAllele <- max(popdata$nalleles)

  af <- array(0, c(numPop,numMarker,numAllele))
  for(cmak in 1:numMarker){
    af[,cmak,1:popdata$nalleles[cmak]] <- t(popdata$allele_freq[[cmak]])
  }
  af2 <- af^2

  PopSt <- array(0, c(numPop,numPop))
  dimnames(PopSt) <- list(popdata$pop_names,popdata$pop_names)
  message("Calculating population ", appendLF=F)
  cprogressPop <- ""
  spop <- 2
  for(i in 1:(numPop-1)){
  for(j in (i+1):numPop){
    message(paste0(rep("\b", nchar(cprogressPop)), collapse=""), appendLF=F)
    cprogressPop <- paste0(i, ":", j, " ")
    message(cprogressPop, appendLF=F); flush.console()

    hs <- 1 - colMeans(apply(af2[c(i,j),,], c(1,2), sum))
    ht <- 1 - rowSums(apply(af[c(i,j),,], c(2,3), mean)^2)
    cn <- 1 / mean(1/numInd[c(i,j)])
    Hs <- 2*cn/(2*cn-1) * hs
    Ht <- ht + Hs/(2*cn*spop)
    mHs <- mean(Hs,na.rm=T)

    PopSt[i,j] <- PopSt[j,i] <- 
      switch(method,
        "GstN" = 1 - mean(hs,na.rm=T)/mean(ht,na.rm=T),
        "GstNC" = 1 - mHs/mean(Ht,na.rm=T),
        "GstH" = (1 - mHs/mean(Ht,na.rm=T)) * (1+mHs)/(1-mHs),
        "DJ" = (mean(Ht,na.rm=T)-mHs)/(1-mHs) * spop/(spop-1),
        NA
      )
  }}
  message(" done.")
  return(PopSt)
}
