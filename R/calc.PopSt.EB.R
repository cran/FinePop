calc.PopSt.EB <-
function(popdata, method){
########### STEP 1: Preprocess ##########
nPop <- popdata$npops             # number of sub pops
nLoci <- popdata$nloci            # number of markers
LocusNames <- popdata$loci_names
PopNames <- popdata$pop_names

nAlleles <- popdata$nalleles
nAllelesMax <- max(nAlleles)
n.lpa <- array(NA, c(nLoci,nPop,nAllelesMax))
for (cl in 1:nLoci){
  obs.an <- t(popdata$obs_allele_num[[cl]])
  check0 <- colSums(obs.an)==0        # remove 0 freq allele
  obs.an <- obs.an[,!check0, drop=F]
  cna <- ncol(obs.an)
  nAlleles[cl] <- cna
  n.lpa[cl,,1:cna] <- obs.an
}
nAllelesMax <- max(nAlleles)
n.lpa <- n.lpa[,,1:nAllelesMax, drop=F]

n.la <- apply(n.lpa, c(1,3), sum, na.rm=T)
n.l <- rowSums(n.la, na.rm=T)
n.lp <- apply(n.lpa, c(1,2), sum, na.rm=T)
n.lp.array <- array(rep(n.lp, nAllelesMax), c(nLoci,nPop,nAllelesMax))

# Excluding one-sided markers (major alelle freq > 99%)
maf_check <- apply(n.la/n.l, 1, max, na.rm=T) > 0.99
if(sum(maf_check)>0){
  message(
    "WARNING:\n",
    " Detected inappropriate markers (major allele frequency > 99%).\n",
    " Locus Names: ", paste(LocusNames[maf_check],collapse=", "))
}


########## STEP 2: Estimate global Fst ########## 
message("Computing global differenciation with ML method... ",appendLF=F);flush.console()
 
p.lpa <- n.lpa / n.lp.array
p.ave.la <- apply(p.lpa, c(1,3), mean, na.rm=T)
p.ave.la.array <- array(rep(p.ave.la, nPop), c(nLoci,nAllelesMax,nPop))
p.ave.la.array <- aperm(p.ave.la.array, c(1,3,2))

sigma2 <- apply((p.lpa - p.ave.la.array)^2, c(1,3), sum, na.rm=T) / (nPop-1)
theta_init <- mean(p.ave.la*(1-p.ave.la)/sigma2, na.rm=T) - 1

nlogL <- function(log.theta){
  theta <- exp(log.theta)
  alpha <- theta * p.ave.la
  alpha.array <- array(rep(alpha,nPop), c(nLoci,nAllelesMax,nPop))
  lgtheta <- lgamma(theta)
  lgalpha <- lgamma(alpha.array)
  clogL <- lgtheta * nLoci * nPop - sum(lgamma(n.lp+theta)) +
           sum(lgamma(n.lap+alpha.array)-lgalpha, na.rm=T)
  return(-clogL)
}
n.lap <- aperm(n.lpa, c(1,3,2))
opt_result <- optim(par=log(theta_init), fn=nlogL, method="L-BFGS-B",
                    lower=log(theta_init/1000), upper=log(theta_init*1000),
                    hessian=T, control=list(maxit=1000))
rm(n.lap)

theta_est <- exp(opt_result$par)

globalFst_est <- 1/(theta_est+1)

message("done.")

########## STEP 3: Estimate pairwise Fst ########## 
message("Computing pairwise differenciation with EB method... ",appendLF=F);flush.console()
Nit <- 1000

GstVec<-function(an1,an2,Nit){
 lam1 <- rgamma(length(an1)*Nit, rep(an1,Nit), scale=1)
 lam1 <- matrix(lam1, nrow=Nit, byrow=T)
 p1 <- lam1 / rowSums(lam1)

 lam2 <- rgamma(length(an2)*Nit, rep(an2,Nit), scale=1)
 lam2 <- matrix(lam2, nrow=Nit, byrow=T)
 p2 <- lam2 / rowSums(lam2)

 mHt <- mean(1-rowSums((0.5*(p1+p2))^2))
 mHs <- mHt - mean(0.25 * rowSums((p1-p2)^2))
 return(c(mHs,mHt))
}

calcHtHsMat <- function(AN){
  pair_Hs <- array(0, c(nPop,nPop,nLoci))
  pair_Ht <- array(0, c(nPop,nPop,nLoci))
  message("pop ",appendLF=F)
  cstep <- ""
  for(cp1 in 1:(nPop-1)){
  for(cp2 in (cp1+1):nPop){
    message(rep("\b",nchar(cstep)),appendLF=F)
    cstep <- paste0(cp1, ":", cp2)
    message(cstep,appendLF=F);flush.console() 
    for(cl in 1:nLoci){
      cna <- nAlleles[cl]
      AN.l1 <- AN[cl,cp1,1:cna]
      AN.l2 <- AN[cl,cp2,1:cna]
      cGstVec <- GstVec(AN.l1, AN.l2, Nit)
      pair_Hs[cp1,cp2,cl] <- pair_Hs[cp2,cp1,cl] <- cGstVec[1]
      pair_Ht[cp1,cp2,cl] <- pair_Ht[cp2,cp1,cl] <- cGstVec[2]
    }
  }}
  message(rep("\b",nchar(cstep)+4),"done.")
  return(list(Hs=pair_Hs,Ht=pair_Ht))
}

AN <- theta_est * p.ave.la.array + n.lpa
HtHsMat <- calcHtHsMat(AN)

ebhs <- apply(HtHsMat$Hs,c(1,2),mean)
ebht <- apply(HtHsMat$Ht,c(1,2),mean)

pairGstH <- (1-ebhs/ebht)*(1+ebhs)/(1-ebhs)
diag(pairGstH) <- 0
pairD <- (ebht-ebhs)/(1-ebhs)*2/(2-1)
pairFst <- 1-ebhs/ebht
diag(pairFst) <- 0

dimnames(pairFst) <- list(PopNames, PopNames)
dimnames(pairGstH) <- list(PopNames, PopNames)
dimnames(pairD) <- list(PopNames, PopNames)

if(method=="Fst"){
########## STEP 4: Variance of Fst ########## 
message("Estimating variance of differenciation with BS method... ",appendLF=F);flush.console()

calcFst <- function(HtHsMat){
  1 - apply(HtHsMat$Hs,c(1,2),mean) / apply(HtHsMat$Ht,c(1,2),mean)
}
bootFst <- function(HtHsMat){
  locus_boot <- sample(1:nLoci, nLoci, replace=T)
  cHs <- HtHsMat$Hs
  cHt <- HtHsMat$Ht
  Hs_boot <- cHs[,,locus_boot]
  Ht_boot <- cHt[,,locus_boot]
  return(calcFst(list(Hs=Hs_boot,Ht=Ht_boot)))
}

nboot <- 100
result_boot <- replicate(nboot,bootFst(HtHsMat))
pairFst_boot_est <- apply(result_boot, c(1,2), mean)
pairFst_boot_sd <- apply(result_boot, c(1,2), sd)

diag(pairFst_boot_est) <- 0
diag(pairFst_boot_sd) <- 0
dimnames(pairFst_boot_est) <- list(PopNames, PopNames)
dimnames(pairFst_boot_sd) <- list(PopNames, PopNames)

message("done.")
}

message("!!! FINISH !!!")

########## STEP 5: Output ##########
result <-
  switch(method,
     "Fst" = list(
               global=list(
                 theta=theta_est,
                 fst=globalFst_est
               ),
               pairwise=list(
                 fst=pairFst,
                 fst.boot=pairFst_boot_est,
                 fst.boot.sd=pairFst_boot_sd
               )
             ),
     "GstH" = pairGstH,
     "DJ" = pairD,
     NA
   )
return(result)
}
