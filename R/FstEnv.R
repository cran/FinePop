FstEnv <-
function(fst.bs, environment, distance=NULL){
  bsrep <- length(fst.bs$bs.pop.list)
  numpop <- length(fst.bs$bs.pop.list[[1]])
  distance <- as.matrix(as.dist(distance))

### model list ###
  gelist <- colnames(environment)
  if(!is.null(distance)){gelist <- c("distance", gelist)}
  funclist <- NULL
  for(cnf in 1:length(gelist)){
    fcomb <- combn(gelist, cnf)
    #fcombI <- apply(fcomb, 2, paste, collapse="*")
    if(cnf==1){
      fcombI <- NULL
    }else{
      fcombI <- apply(fcomb, 2, paste, collapse="+")
      fcombI <- paste0("(",fcombI,")^2")
    }
    fcombS <- apply(fcomb, 2, paste, collapse="+")
    fcomb <- unique(as.vector(rbind(fcombI,fcombS)))
    funclist <- c(funclist, fcomb)
  }
  funclist <- paste0("fst~-1+", funclist)

### model fitting and TIC ###
  out.list <- list() 
  for(cmodel in 1:length(funclist)){
    infunc <- funclist[cmodel]
    ## lm - bootstrap ##
    efflist <- NULL; r2list <- NULL
    for(crep in bsrep:0){
      cfstmat <- NULL; cenv <- NULL; cgeo <- NULL

      if(crep==0){
        cfstmat <- fst.bs$org.fst
        bspops <- 1:numpop
      }else{
        cfstmat <- fst.bs$bs.fst.list[[crep]]
        bspops <- fst.bs$bs.pop.list[[crep]]
      }
      cgeo <- matrix(0, nrow=numpop, ncol=numpop)
      for(cpop1 in 1:(numpop-1)){
      for(cpop2 in (cpop1+1):numpop){
        cbspop1 <- min(bspops[cpop1], bspops[cpop2])
        cbspop2 <- max(bspops[cpop1], bspops[cpop2])
        cgeo[cpop1,cpop2] <- distance[cbspop1, cbspop2]
      }}
      cenv <- environment[bspops,]
      num.env <- ncol(environment)

      clmmat <- NULL
      for(cpop1 in 1:(numpop-1)){
      for(cpop2 in (cpop1+1):numpop){
        cenvvec <- NULL
        for(ce in 1:num.env){cenvvec <- c(cenvvec, abs(cenv[cpop1,ce]-cenv[cpop2,ce]))}
        clmmat <- rbind(clmmat, c(cfstmat[cpop1,cpop2], cgeo[cpop1,cpop2], cenvvec))
      }}#pops
      colnames(clmmat) <- c("fst",gelist)
      clmmat.scale <- data.frame(scale(data.frame(clmmat)))
      lmresult <- lm(infunc, data=clmmat.scale)
      efflist <- rbind(efflist, lmresult$coeff)
      r2list <- c(r2list, summary(lmresult)$r.sq)
    }#rep
    efflist <- efflist[-nrow(efflist),]
    r2list <- r2list[-length(r2list)]

    Scoeff <- if(is.vector(efflist)){T}else{F}

    lmresult0 <- summary(lmresult)$coeff
    colnames(lmresult0) <- c("coeff0","SE","t","p0")
    varnames <- rownames(lmresult0)

    r2mean <- mean(r2list)
    coeffSD <- if(Scoeff){sd(efflist)}else{apply(efflist, 2, sd)}
    coeffM <- if(Scoeff){mean(efflist)}else{apply(efflist, 2, mean)}
    coeffZ0 <- lmresult0[,"coeff0"] / coeffSD
    coeffZb <- coeffM / coeffSD
    cpv0 <- pnorm(abs(coeffZ0), lower.tail=F) * 2
    cpvb <- pnorm(abs(coeffZb), lower.tail=F) * 2
    lm_summary <- cbind(coeffB=coeffM, SD=coeffSD, ZB=coeffZb, pB=cpvb, ZB0=coeffZ0, pB0=cpv0)
    rownames(lm_summary) <- varnames

    coeffSE <- lmresult0[,"SE"]
    scoeffSESDsq <- sum((coeffSD/coeffSE)^2)
    cTICsdse <- -2 * as.numeric(logLik(lmresult)) + 2 * scoeffSESDsq

    cTICtrcov <- NA
    if(!Scoeff){
      covA <- vcov(lmresult)
      covB <- cov(efflist)
      covAB <- solve(covA) %*% covB
      cTICtrcov <- -2 * as.numeric(logLik(lmresult)) + 2 * sum(diag(covAB),na.rm=T)
    }

    cTICout <- if(Scoeff){cTICsdse}else{cTICtrcov}

    all.result <- cbind(lmresult0, lm_summary)
    all.result <- all.result[,c("coeff0","SD", "ZB0", "pB0"), drop=F]
    colnames(all.result) <- c("Estimate","SD","Z","p")

    out.list[[cmodel]] <- list(
      model=infunc,
      coefficients=all.result,
      TIC=cTICout,
      R2=summary(lmresult)$r.sq
    ) 
  }# model fitting, TIC

  class(out.list) <- "FstEnv"
  return(out.list)
}
