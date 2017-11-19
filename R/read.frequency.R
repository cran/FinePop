read.frequency <-
function(frequency, popname=NULL){
# read frequency file
LFx <- intToUtf8(0x0A)
HTx <- intToUtf8(0x09)
all_lines <- scan(frequency, what=character(), quiet=T, sep=LFx, blank.lines.skip=F)
all_lines <- tolower(all_lines)
all_lines <- gsub(HTx, " ", all_lines)
all_lines <- gsub(" {1,}"," ", all_lines)
all_lines <- gsub(" $", "", all_lines)

# num. subpopulation
line_numpop <- grep("#the number of subpopulations", all_lines) + 1
numPop <- as.integer(all_lines[line_numpop])

# pop name
if(is.null(popname)){
  pop.names <- paste0("pop",1:numPop)
}else{
  pop.names <- scan(popname, what=character(), quiet=T, blank.lines.skip=T)
}

# num. locus
line_numlocus <- grep("#the number of loci", all_lines) + 1
numMarker <- as.integer(all_lines[line_numlocus])

# allele freq
lines_af <- grep("#locus", all_lines)
lines_af <- cbind(start=lines_af+1,end=lines_af+numPop)

AlleleCount <- list()
AlleleFreq <- list()
numAlleles <- rep(0,numMarker)
AlleleList <- list()
IndObs <- list()
for(cmak in 1:numMarker){
  cac <- all_lines[lines_af[cmak,"start"]:lines_af[cmak,"end"]]
  cac <- matrix(as.integer(unlist(strsplit(cac, " "))), nrow=numPop, byrow=T)
  cac <- cac[,!is.na(cac[1,])] 
  cac <- t(cac)
  cni <- colSums(cac)
  caf <- t(t(cac)/cni)
  cna <- nrow(cac)
  ca <- paste0("000", 1:cna)
  ca <- substr(ca, nchar(ca)-2, nchar(ca))
  dimnames(cac) <- dimnames(caf) <- list(ca,pop.names)

  AlleleCount[[cmak]] <- cac
  AlleleFreq[[cmak]] <- caf
  numAlleles[cmak] <- cna
  AlleleList[[cmak]] <- ca
  IndObs[[cmak]] <- cni
}

MarkerList <- paste0("Locus", 1:numMarker)
numCall <- matrix(unlist(IndObs),nrow=numMarker,ncol=numPop,byrow=T)
numInd <- apply(numCall, 2, max)
CallRate <- rowSums(numCall)/sum(numInd)

return(list(
  obs_allele_num=AlleleCount,
  allele_freq=AlleleFreq,
  indtyp=IndObs,
  npops=numPop,
  pop_sizes=numInd,
  pop_names=pop.names,
  nloci=numMarker,
  loci_names=MarkerList,
  nalleles=numAlleles,
  call_rate=CallRate,
  all_alleles=AlleleList
))
}
