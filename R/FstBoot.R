FstBoot <-
function(popdata, fst.method="EBFST", bsrep=100, log.bs=F, locus=F){
  calcFst <- eval(parse(text=fst.method))
### rawdata ###
  gpdata0 <- popdata
  numpop <- gpdata0$npops
  numloci <- gpdata0$nloci
  if(fst.method=="EBFST"){
    fst0 <- calcFst(gpdata0, locus=locus)
    fst0.l <- fst0$pairwise$fst.locus
    fst0 <- fst0$pairwise$fst
  }else{
    fst0 <- calcFst(gpdata0)
  }

### loc & ind resample ###
gp_locind_resample <- function(popdata, crep){
  gpdat <- popdata$pop_list
  indnames <- popdata$ind_names

  gpout <- c("bs", popdata$loci_names)
  select_pop <- sample(1:numpop, numpop, replace=T)
  for(cpop in 1:numpop){
    cbspop <- select_pop[cpop]
    cgp <- gpdat[[cbspop]]
    cnind <- nrow(cgp)
    select_ind <- sample(1:cnind, cnind, replace=T)
    cgp <- cbind(indnames[[cbspop]][select_ind], ",", cgp[select_ind,])
    cgp <- apply(cgp, 1, paste, collapse=" ")
    gpout <- c(gpout, "POP", cgp)
  }
  cat(gpout, file="finepopbs.tmp", sep="\n")
  gpdat <- read.genepop("finepopbs.tmp")
  if(log.bs){file.rename("finepopbs.tmp", paste0("gtdata_bs",crep,".txt"))
  }else{file.remove("finepopbs.tmp")}
  return(list(gp=gpdat, pop=select_pop))
}# function loc & ind resample

### bootstrap cal fst ###
  bs.poplist <- list()
  bs.fstlist <- list()
  bs.fstlist.l <- list()
  crep <- 1
  trial.success <- ""
  while(crep <= bsrep){
    message("Bootstrapping ", crep, "/", bsrep, " ", trial.success)
    gc()
    bs.gpdata <- gp_locind_resample(gpdata0, crep)
    cgp <- bs.gpdata$gp
    cpop <- bs.gpdata$pop

    cfst <- cfst.l <- NULL
    try({
      if(fst.method=="EBFST"){
        cfst <- calcFst(cgp, locus=locus)
        cfst.l <- cfst$pairwise$fst.locus
        cfst <- cfst$pairwise$fst
      }else{
        cfst <- calcFst(cgp)
      }
    })

    if(sum(!is.finite(cfst))==0 & !is.null(cfst)){
      bs.poplist[[crep]] <- cpop
      bs.fstlist[[crep]] <- cfst
      bs.fstlist.l[[crep]] <- cfst.l
      crep <- crep+1
      trial.success <- ""
    }else{
      trial.success <- "retry"
    }
  }#bs calc fst

  result <- list(
     bs.pop.list=bs.poplist,
     bs.fst.list=bs.fstlist,
     org.fst=fst0
  )
  if(fst.method=="EBFST"&locus==T){
   result$bs.fst.list.locus <- bs.fstlist.l
   result$org.fst.locus <- fst0.l
  }

  return(result)
}
