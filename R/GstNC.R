GstNC <-
function(popdata){
message("Calculating Gst (Nei&Chesser 1983)... ", appendLF=F); flush.console()
numPop <- popdata$npops
numMarker <- popdata$nloci
numInd <- popdata$pop_sizes
numAllele <- max(popdata$nalleles)
#nadigit <- paste(rep("0",nchar(popdata$pop_allele[[1]][[1]][1])),collapse="") ### eq9$11

af <- array(0, c(numPop,numMarker,numAllele))
for(cmak in 1:numMarker){
  af[,cmak,1:popdata$nalleles[cmak]] <- t(popdata$allele_freq[[cmak]])
}
af2 <- af^2

#sumX <- array(0,c(numPop,numMarker)) ### eq9&11
#for(i in 1:numPop){
#  cgt1 <- popdata$pop_allele[[i]][[1]]; cgt1 <- gsub(nadigit,NA,cgt1)
#  cgt2 <- popdata$pop_allele[[i]][[2]]; cgt2 <- gsub(nadigit,NA,cgt2)
#  sumX[i,] <- colMeans(cgt1==cgt2, na.rm=T)
#}

gstNC <- array(0, c(numPop,numPop))
dimnames(gstNC) <- list(popdata$pop_names,popdata$pop_names)
message("Populations ", appendLF=F)
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

  #H0 <- 1 - colMeans(sumX[c(i,j),])    ### eq9&11 
  #Hs <- cn/(cn-1) * (hs - H0/(2*cn))
  #Ht <- ht + Hs/(cn*2) - H0/(2*cn*spop)

  Hs <- 2*cn/(2*cn-1) * hs              ### eq15&16
  Ht <- ht + Hs/(2*cn*spop)

  gstNC[i,j] <- 1 - mean(Hs,na.rm=T)/mean(Ht,na.rm=T)
  #gstNC[i,j] <- 1 - mean(Hs/Ht,na.rm=T)
}}
gstNC <- t(gstNC)
gstNC[upper.tri(gstNC,diag=T)] <- NA
gstNC <- gstNC[-1,-numPop]
message(" done.")

return(gstNC)
}
