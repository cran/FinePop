thetaWC.pair <-
function(popdata){
numpop <- popdata$npops
numlocus <- popdata$nloci
rpop <- 2

cfstmat <- matrix(0, nrow=numpop, ncol=numpop)
cat("Populations ")
cprogressPop <- ""
for(cpop1 in 1:(numpop-1)){
for(cpop2 in (cpop1+1):numpop){
cat(paste0(rep("\b", nchar(cprogressPop)), collapse=""))
cprogressPop <- paste0(cpop1, ":", cpop2, " ")
cat(cprogressPop); flush.console()

sum_a <- 0
sum_abc <- 0
bucket <- numlocus%/%10; if(bucket<10){bucket <- numlocus}
for(cloc in 1:numlocus){
  nsamples12 <- popdata$indtyp[[cloc]][c(cpop1, cpop2)]

  genotype1A <- popdata$pop_allele[[cpop1]][[1]][,cloc]
  genotype1a <- popdata$pop_allele[[cpop1]][[2]][,cloc]
  genotype2A <- popdata$pop_allele[[cpop2]][[1]][,cloc]
  genotype2a <- popdata$pop_allele[[cpop2]][[2]][,cloc]
  allele_names <- rownames(popdata$allele_freq[[cloc]])

  for(callele in 1:((popdata$nalleles)[cloc])){
    n_bar <- mean(nsamples12)
    n_c <- (n_bar * rpop - sum(nsamples12^2)/(n_bar*rpop)) / (rpop-1)

    freqA <- popdata$allele_freq[[cloc]][callele,c(cpop1,cpop2)]
    p_bar <- sum(freqA * nsamples12) / (rpop*n_bar)

    s2 <- sum(nsamples12 * (freqA - p_bar)^2) / ((rpop-1)*n_bar)
    #s2 <- sum(nsamples12 * (freqA - p_bar)^2) / ((rpop-1)*n_bar) * (rpop-1)/rpop

    callelename <- allele_names[callele]
    hetero1 <- ((genotype1A==callelename)&(genotype1a!=callelename)) |
               ((genotype1A!=callelename)&(genotype1a==callelename))
    hetero1 <- hetero1[!is.na(hetero1)]
    hetero2 <- ((genotype2A==callelename)&(genotype2a!=callelename)) |
               ((genotype2A!=callelename)&(genotype2a==callelename))
    hetero2 <- hetero2[!is.na(hetero2)]
    freqAH <- c(sum(hetero1),sum(hetero2))/nsamples12
    #freqAH <- 2 * freqA * (1-freqA)
    #freqAH <- 1 - sum(freqA^2)
    h_bar <- sum(freqAH * nsamples12) / (rpop*n_bar)
             #0

    #WCa <- n_bar/n_c * (s2 - 1/(n_bar-1) * (p_bar*(1-p_bar) - (rpop-1)/rpop * s2 - h_bar/4 ))
    WCa <- (rpop-1)/rpop * n_bar/n_c * (s2 - 1/(n_bar-1) * (p_bar*(1-p_bar) - (rpop-1)/rpop * s2 - h_bar/4 ))
    #WCa <- n_bar/n_c * (s2 - 1/(n_bar-1) * (p_bar*(1-p_bar) - s2 - h_bar/4 ))
    #WCa <- n_bar/n_c * (s2*(rpop-1)/rpop - 1/(n_bar-1) * (p_bar*(1-p_bar) - (rpop-1)/rpop * s2 - h_bar/4 ))

    WCb <- n_bar/(n_bar-1) * (p_bar*(1-p_bar) - (rpop-1)/rpop * s2 - (2*n_bar-1)/(4*n_bar) * h_bar)
    #WCb <- n_bar/(n_bar-1) * (p_bar*(1-p_bar) - s2 - (2*n_bar-1)/(4*n_bar) * h_bar)
    #WCb <- p_bar*(1-p_bar)# - (rpop-1)/rpop * s2 - (2*n_bar-1)/(4*n_bar) * h_bar)

    WCc <- h_bar/2
    #WCc <- 0

    if(is.finite(WCa)){
       sum_a <- sum_a + WCa;
       sum_abc <- sum_abc + WCa + WCb + WCc
    }
  }#allele
}#loc

cfstmat[cpop1, cpop2] <- (sum_a / sum_abc)
}}#pop1 pop2
cat("\n")

dimnames(cfstmat) <- list(popdata$pop_names,popdata$pop_names)
cfstmat <- t(cfstmat)
cfstmat[upper.tri(cfstmat,diag=T)] <- NA
cfstmat <- cfstmat[-1,-numpop]
return(cfstmat)
}
