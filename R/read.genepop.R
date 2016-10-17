read.genepop <-
function(genepop, popname=NULL){
# read genepop file
all_lines <- scan(genepop, what=character(), quiet=T, sep="\n", blank.lines.skip=F)

# title 
gp_title <- all_lines[1]

# locus pop sample count
cline <- gsub(" ", "", all_lines)
cline <- gsub("\t", "", cline)
poploc <- which(toupper(cline)=="POP")
MarkerList <- cline[2:(poploc[1]-1)]
MarkerList <- gsub(",", "\n", MarkerList)
MarkerList <- unlist(strsplit(MarkerList, "\n"))
numMarker <- length(MarkerList)
numPop <- length(poploc)
popstart <- poploc + 1
popend <- c(poploc[-1]-1, length(all_lines))
numInd <- popend - popstart + 1
numIndAll <- sum(numInd)
PopID <- rep(1:numPop, numInd)
rm(cline);gc()

if(is.null(popname)){
  pop.names <- paste0("pop",1:numPop)
}else{
  pop.names <- scan(popname, what=character(), quiet=T, blank.lines.skip=T)
}

# marker genotype
gtdata <- all_lines[-poploc]
gtdata <- gtdata[-(1:(poploc[1]-1))]
gtdata <- unlist(strsplit(gtdata, ","))
IndID <- gtdata[c(T,F)]
IndID <- gsub(" ", "", IndID)
IndID <- gsub("\t", "", IndID)
gtdata <- gsub(" ", "\t", gtdata[c(F,T)])
gtdata <- matrix(unlist(strsplit(gtdata, "\t")), nrow=numIndAll, byrow=T)
gtdata <- gtdata[,-which(colSums(gtdata=="")!=0)]
gp_digit <- as.integer(nchar(as.character(gtdata[1,1]))/2)
gp_na <- paste(rep("0", gp_digit), collapse="")
rm(all_lines);gc()

# gpdata
htdata1 <- substr(gtdata,1,gp_digit)
htdata2 <- substr(gtdata,gp_digit+1,gp_digit*2)

haplo <- list()
diplo <- list()
ind_names <- list()
for(cpop in 1:numPop){
  cpopind <- PopID==cpop
  haplo[[cpop]] <- list(htdata1[cpopind,],htdata2[cpopind,])
  diplo[[cpop]] <- gtdata[cpopind,]
  ind_names[[cpop]] <- IndID[cpopind]
}

htdata <- rbind(htdata1, htdata2)
rm(htdata1, htdata2, gtdata);gc()

AlleleCount <- list()
AlleleFreq <- list()
IndObs <- list()
numAlleles <- rep(0,numMarker)
AlleleList <- list()
CallRate <- rep(0, numMarker)
for(cm in 1:numMarker){
  cgt <- table(htdata[,cm], c(PopID, PopID), exclude=gp_na)
  colnames(cgt) <- NULL
  numAlleles[cm] <- nrow(cgt)
  AlleleList[[cm]] <- rownames(cgt)
  cgtnum <- colSums(cgt)
  AlleleCount[[cm]] <- cgt
  numcall <- as.integer(cgtnum/2)
  IndObs[[cm]] <- numcall
  CallRate[cm] <- sum(numcall)/numIndAll
  AlleleFreq[[cm]] <- t(t(cgt) / cgtnum) 
}

IndNames <- list()
for(cpop in 1:numPop){IndNames[[cpop]] <- IndID[PopID==cpop]}

return(list(pop_allele=haplo,
            pop_list=diplo,
            obs_allele_num=AlleleCount,
            allele_freq=AlleleFreq,
            indtyp=IndObs,
            npops=numPop,
            pop_sizes=numInd,
            pop_names=pop.names,
            ind_names=IndNames,
            nloci=numMarker,
            loci_names=MarkerList,
            nalleles=numAlleles,
            call_rate=CallRate,
            all_alleles=AlleleList
            ))
}
