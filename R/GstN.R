GstN <-
function(popdata){
numPop <- popdata$npops
numMarker <- popdata$nloci
numInd <- popdata$pop_sizes
numAllele <- max(popdata$nalleles)

af <- array(0, c(numPop,numMarker,numAllele))
for(cmak in 1:numMarker){
  af[,cmak,1:popdata$nalleles[cmak]] <- t(popdata$allele_freq[[cmak]])
}
af2 <- af^2

# Gst (Nei 1973)
message("Calculating Gst(Nei1973)... ", appendLF=F); flush.console()
gstN<-array(0, c(numPop,numPop))
dimnames(gstN) <- list(popdata$pop_names,popdata$pop_names)
message("Populations ", appendLF=F)
cprogressPop <- ""
for(i in 1:(numPop-1)){
for(j in (i+1):numPop){
  message(paste0(rep("\b", nchar(cprogressPop)), collapse=""), appendLF=F)
  cprogressPop <- paste0(i, ":", j, " ")
  message(cprogressPop, appendLF=F); flush.console()
  hs <- 1 - colMeans(apply(af2[c(i,j),,], c(1,2), sum))
  ht <- 1 - rowSums(apply(af[c(i,j),,], c(2,3), mean)^2)
  gstN[i,j] <- 1 - mean(hs,na.rm=T)/mean(ht,na.rm=T)
}}
gstN <- t(gstN)
gstN[upper.tri(gstN,diag=T)] <- NA
gstN <- as.matrix(as.dist(gstN))
#gstN <- gstN[-1,-numPop]
message(" done.")

return(gstN)
}
