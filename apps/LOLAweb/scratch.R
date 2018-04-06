library(LOLA)
library(ggplot2)
library(GenomicRanges)

source("themes.R")


dbPath = system.file("extdata", "hg19", package="LOLA")
regionDB = loadRegionDB(dbLocation=dbPath)

data("sample_universe", package="LOLA")
# data("sample_input", package="LOLA")

userSet = read.table(file = "lola_vignette_data/setB_100.bed", header = F) 
colnames(userSet) <- c('chr','start','end','id','score','strand')
userSet <- with(userSet, GRanges(chr, IRanges(start+1, end), strand, score, id=id))
userSetsRedefined =	redefineUserSets(userSet, userUniverse)

resRedefined = runLOLA(userSetsRedefined, userUniverse, regionDB, cores=1)

resRedefined = subset(as.data.frame(resRedefined),
                      oddsRatio > 2,
                      pValueLog > 2,
                      support > 2
)

summary(resRedefined$pValueLog)

hist(resRedefined$pValueLog)

View(resRedefined)
View(resRedefined$pValueLog)
resRedefined[order(pValueLog)]$pValueLog[5]
resRedefined
dplyr::arrange(resRedefined, userSet)
sort(resRedefined, collection)
resRedefined[order(resRedefined[,"collection"])]
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

?sort

resRedefined
resRedefined[order(meanRnk, maxRnk)]

system.time({
  dbPath = system.file("extdata", "hg19", package="LOLA")
  regionDB = loadRegionDB(dbLocation="~/Downloads/scratch/ns5bc/resources/regions/LOLACore/hg19/")
})

system.time({
  dbPath = system.file("extdata", "hg38", package="LOLA")
  regionDB = loadRegionDB(dbLocation="~/Downloads/scratch/ns5bc/resources/regions/LOLACore/hg38/")
})

system.time({
  dbPath = system.file("extdata", "mm10", package="LOLA")
  regionDB = loadRegionDB(dbLocation="~/Downloads/scratch/ns5bc/resources/regions/LOLACore/mm10/")
})



ggplot(resRedefined, aes(reorder(description, sort(b, decreasing = F)), oddsRatio, fill = userSet)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  theme_ns()


regionDB = loadRegionDB("reference/LOLAJaspar/hg19/")
# regionDB = loadRegionDB("reference/LOLARoadmap/hg19/")
# regionDB = loadRegionDB("reference/CORE/hg19/")

userSet = read.table(file = "lola_vignette_data/setB_100.bed", header = F, stringsAsFactors = F) 
colnames(userSet) <- c('chr','start','end','id','score','strand')
userSet <- with(userSet, GRanges(chr, IRanges(start+1, end), strand, score, id=id))
userSet = GRangesList(userSet)



userUniverse = read.table(file = "universes/hg19/tiles1000.hg19.bed", header = F, stringsAsFactors = F)
colnames(userUniverse) <- c('chr','start','end','id','score','strand')
userUniverse <- with(userUniverse, GRanges(chr, IRanges(start+1, end), strand, score, id=id))


userSetsRedefined =	redefineUserSets(userSet, userUniverse)

resRedefined = runLOLA(userSetsRedefined,
                       userUniverse,
                       regionDB,
                       cores=1)


# GenomicDistributions plots

library(GenomicDistributions)
library(GenomicRanges)

query = read.table(file = "lola_vignette_data/setC_complete.bed", header = F)
#query = read.table(file = "lola_vignette_data/setB_100.bed", header = F) 
colnames(query) = c('chr','start','end','id','score','strand')
query = with(query, GRanges(chr, IRanges(start+1, end), strand, score, id=id))

x = genomicDistribution(query, "hg19")

# Then, plot the result:
plotGenomicDist(x)



EnsDb = loadEnsDb("hg19")
featsWide = ensembldb::genes(EnsDb, columns=NULL)

# Now, grab just a single base pair at the TSS
feats = promoters(featsWide, 1, 1)

# Change from ensembl-style chrom annotation to UCSC_style
seqlevels(feats) = paste0("chr", seqlevels(feats))

system.time({
  featureDistance = featureDistribution(query, feats)
})

# Then plot the result:
plotFeatureDist(featureDistance)


featureDistance2 = featureDistribution(queryList, feats)
plotFeatureDist(featureDistance2)



keyphrase <- "486TQ3MHFJPU2D9"
loadCaches(keyphrase, assignToVariable = "cipher", loadEnvir = globalenv(), cacheDir = "cache")

cipher <- get("cipher", envir = globalenv())

# keyphrase
key <- hash(charToRaw(keyphrase))

dat <- data_decrypt(cipher, key)

res <- unserialize(dat)

res$resRedefined$id <- paste(res$resRedefined$description, res$resRedefined$dbSet, sep = "_")
res$resRedefined$axis_label <- strtrim(res$resRedefined$description, 50)

res$resRedefined %>%
  filter(maxRnk < 90) %>%
  arrange(desc(meanRnk)) %>%
  # ggplot(aes(rev(reorder(rev(axis_label), eval(parse(text = "meanRnk")))), oddsRatio, fill = userSet, group = id)) +
  ggplot(aes(reorder(axis_label, rev(eval(parse(text = "meanRnk")))), oddsRatio, fill = userSet, group = id)) +
  # ggplot(aes(reorder(axis_label, eval(parse(text = "meanRnk"))), oddsRatio, fill = userSet, group = id)) +
  # ggplot(aes(reorder(axis_label, eval(parse(text = "maxRnk"))), oddsRatio, fill = userSet, group = id)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Description") +
  ylab("Odds Ratio") +
  coord_flip() +
  theme_ns()

tst <-
  res$resRedefined %>%
  filter(maxRnk < 90) %>%
  mutate(invRnk = 1 - (maxRnk / nrow(res$resRedefined)))


tst_reordered <- tst[order(tst$meanRnk),]


as.character(tst_reorder)
rev(tst_reorder)


keyphrase <- "GN2ZKFXCHW7TRPA"
env <- new.env()

simpleCache(keyphrase, assignToVariable = "cipher", loadEnvir = environment(), cacheDir = "cache")

cipher <- get("cipher", envir = env)

# keyphrase
key <- hash(charToRaw(keyphrase))

dat <- data_decrypt(cipher, key)

res <- unserialize(dat)



plot_input <- function(res, metric, ylabel, sortcol) {
  
  # conditional for inverting rank sorting
  
  if(grepl("rnk", sortcol, ignore.case = TRUE)) {
    
    # need to order data frame by sort col if it's a rank
    dat <- res[order(as.data.frame(res)[,sortcol]), ]
    
    # now construt base layer for plot with reverse on the sort
    p <- ggplot(dat, aes(reorder(axis_label, rev(eval(parse(text = sortcol)))), eval(parse(text = metric)), fill = userSet, group = id))
    
  } else {
    
    p <- ggplot(res, aes(reorder(axis_label, eval(parse(text = sortcol))), eval(parse(text = metric)), fill = userSet, group = id))
    
  }
  
  p +
    geom_bar(stat = "identity", position = "dodge") +
    xlab("Description") +
    ylab(ylabel) +
    coord_flip() +
    theme_ns()
  
}

library(ggplot2)

y <- "mpg"

ggplot(data = mtcars, aes(wt,eval(parse(text = y)))) +
  geom_point()


simpleCache()

dbPath = system.file("extdata", "hg19", package="LOLA")

LOLARoadmap_hg38 = loadRegionDB(dbLocation="reference/LOLARoadmap/hg38")

Core_hg19 = loadRegionDB(dbLocation="reference/Core/hg19")
object.size(Core_hg19)

LOLARoadmap_hg19 = loadRegionDB(dbLocation="reference/LOLARoadmap/hg19")
object.size(LOLARoadmap_hg19)


x <- system.time({
  
  rnorm(1000000)
})



##### ggplotly

# userSet = readBed(file = "lola_vignette_data/setB_100.bed") 
# dbPath = "reference/Core/hg19/"
# regionDB = loadRegionDB(dbLocation=dbPath)
# 
# userUniverse = readBed("universes/hg19/tiles1000.hg19.bed")
# # userUniverse = read.table(file = "universes/hg38/tiles.hg38.5000.bed", header = F) 
# # colnames(userUniverse) <- c('chr','start','end','id','score','strand')
# # userUniverse <- with(userUniverse, GRanges(chr, IRanges(start+1, end), strand, score, id=id))
# 
# userSetsRedefined =	redefineUserSets(userSet, userUniverse)
# 
# 
# resRedefined = runLOLA(userSetsRedefined, userUniverse, regionDB, cores=1)

library(sodium)

keyphrase <- "486TQ3MHFJPU2D9"
loadCaches(keyphrase, assignToVariable = "cipher", loadEnvir = globalenv(), cacheDir = "cache")

cipher <- get("cipher", envir = globalenv())

# keyphrase
key <- hash(charToRaw(keyphrase))

dat <- data_decrypt(cipher, key)

res <- unserialize(dat)

res$resRedefined$id <- paste(res$resRedefined$description, res$resRedefined$dbSet, sep = "_")
res$resRedefined$axis_label <- strtrim(res$resRedefined$description, 50)

library(magrittr)

p <-
  res$resRedefined %>%
  filter(maxRnk < 90) %>%
  arrange(desc(meanRnk)) %>%
  # ggplot(aes(description, eval(parse(text = "meanRnk")), oddsRatio, fill = userSet, group = id)) +
  ggplot(aes(reorder(axis_label, eval(parse(text = "meanRnk"))), rev(oddsRatio), fill = userSet, group = id)) +
  # ggplot(aes(reorder(axis_label, eval(parse(text = "maxRnk"))), oddsRatio, fill = userSet, group = id)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Description") +
  ylab("Odds Ratio") +
  coord_flip() +
  theme_ns()

library(plotly)
  
ggplotly(p)


q <- 
  res$resRedefined %>%
  filter(maxRnk < 100) %>%
  arrange(desc(meanRnk)) %>%
  ggplot(aes(pValueLog, oddsRatio, size = log(support), 
             text = paste0("Collection: ", 
                           collection, 
                           "\n",
                           "Description: ",
                           axis_label))) +
  geom_point() +
  xlab("P Value Log") +
  ylab("Odds Ratio") +
  # coord_flip() +
  geom_blank(aes(text = collection)) +
  theme_ns() +
  theme(legend.position = "bottom")

ggplotly(q, tooltip = c("x", "y", "size", "text"))

ggplotly(q, tooltip = c("x", "y", "size")) %>%
  layout(showlegend = TRUE, legendgroup = "size")
