clip.genepop.name <-
function(infile, outfile, remove.list){
####
cat("Reading GENEPOP file... ");flush.console()
all_lines <- scan(infile, what=character(), quiet=T, sep="\n", blank.lines.skip=F)
cat("done.\n"); flush.console()

####
cat("Processing  GENEPOP data... "); flush.console()
# title 
genepop_title <- all_lines[1]

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
PopID <- rep(1:numPop, numInd)

rm(cline);gc()

# marker genotype
gtdata <- all_lines[-poploc]
gtdata <- gtdata[-(1:(poploc[1]-1))]
gtdata <- unlist(strsplit(gtdata, ","))

IndID <- gtdata[c(T,F)]
gtdata <- gsub(" ", "\t", gtdata[c(F,T)])
gtdata <- matrix(unlist(strsplit(gtdata, "\t")), nrow=sum(numInd), byrow=T)
gtdata <- gtdata[,-which(colSums(gtdata=="")!=0)]

genepop_digit <- as.integer(nchar(as.character(gtdata[1,1]))/2)

if(genepop_digit==2){genepop_na <- "00"
}else if(genepop_digit==3){genepop_na <- "000"
}else{cat("!!! Disorder in reading GENEPOP file. !!!\n");return(0)}

rm(all_lines); gc()
cat("done.\n"); flush.console()

# remove
remloci <- match(remove.list, MarkerList)
numNG <- length(remove.list)
numOK <- numMarker - numNG

####
cat("Writing GENEPOP file without low MAF markers... "); flush.console()

gtdata <- cbind(paste0(IndID, ","), gtdata[,-remloci])
gtdata <- apply(gtdata, 1, paste, collapse="\t")

outline <- rep("", 1 + numOK + numPop + sum(numInd))
outline[1] <- paste0("CLIPPED: ", genepop_title)
outline[2:(numOK+1)] <- MarkerList[-remloci]
outline[poploc-numNG] <- "POP"
for(cp in 1:numPop){outline[popstart[cp]:popend[cp]-numNG] <- gtdata[PopID==cp]}
cat(outline, file=outfile, sep="\n")

cat("done.\n")
cat("Number of populations = ", numPop, "\n",
    "Number of markers = ", numMarker, " -> ", numOK, "\n", sep="")
cat(numNG, "markers removed.\n")
}
