library(LOLA)
library(ggplot2)
library(GenomicRanges)

source("themes.R")


dbPath = system.file("extdata", "hg19", package="LOLA")
regionDB = loadRegionDB(dbLocation=dbPath)
data("sample_universe", package="LOLA")
data("sample_input", package="LOLA")

userSetsRedefined =	redefineUserSets(userSets, userUniverse)

resRedefined = runLOLA(userSetsRedefined, userUniverse, regionDB, cores=1)

hist(resRedefined$pValueLog)
  
# userSetA = reduce(do.call(c, (LOLA::sampleGRL(regionDB$regionGRL,
#                                         prop=c(.1,.25,.05,.05,0)))))
# userSetB = reduce(do.call(c, (sampleGRL(regionDB$regionGRL,
#                                         prop=c(.2,.05,.05,.05,0)))))

userSetB = read.table(file = "lola_vignette_data/setB_100.bed", header = F) 
colnames(userSetB) <- c('chr','start','end','id','score','strand')
userSetB <- with(userSetB, GRanges(chr, IRanges(start+1, end), strand, score, id=id))

userSetC = read.table(file = "lola_vignette_data/setC_100.bed", header = F) 
colnames(userSetC) <- c('chr','start','end','id','score','strand')
userSetC <- with(userSetC, GRanges(chr, IRanges(start+1, end), strand, score, id=id))

userSets = GRangesList(list(setB=userSetB, setC=userSetC))
userSets = GRangesList(setB=userSetB)

userUniverse = read.table(file = "lola_vignette_data/activeDHS_universe.bed", header = F) 
colnames(userUniverse) <- c('chr','start','end','id','score','strand')
userUniverse <- with(userUniverse, GRanges(chr, IRanges(start+1, end), strand, score, id=id))

userSetsRedefined =	redefineUserSets(userSets, userUniverse)

resRedefined = runLOLA(userSetsRedefined, userUniverse, regionDB, cores=1)

ggplot(resRedefined, aes(logOddsRatio)) +
  geom_histogram() +
  facet_wrap(~userSet, ncol = 1) +
  theme_ns()
