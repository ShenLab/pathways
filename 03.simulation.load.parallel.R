library(doParallel)
simulation_time <- 20000
cl <- makeCluster(48)
registerDoParallel(cl)
load("data/simulation.Rdata")
source("pathway_burden_analysis.R")
# change pathway gene IDs to symbols
pathways.data$Symbol.cls3 = rep(NA, dim(pathways.data)[1])
Symbol.cls3 <- 
foreach(i=1:dim(pathways.data)[1], .combine = rbind) %dopar% {
  geneID <- as.character(pathways.data$genes.cls3[i])
  LGD.geneID <- strsplit(geneID, split = '\\|')[[1]][1]
  mis.geneID <- strsplit(geneID, split = '\\|')[[1]][2]
  LGD.genes <- strsplit(LGD.geneID, split = ", ")[[1]]
  mis.genes <- strsplit(mis.geneID, split = ", ")[[1]]
  LGD.genes <- as.character(rescaled.reference$HGNC[match(LGD.genes, rescaled.reference$GeneID)])
  mis.genes <- as.character(rescaled.reference$HGNC[match(mis.genes, rescaled.reference$GeneID)])
  outID <- paste(toString(LGD.genes), toString(mis.genes), sep = "|")
  outID
}
pathways.data$Symbol.cls3 <- Symbol.cls3
pathwayname = 'c5.all'
wanted.pathways <- pathways
wanted.pathways.data <- pathways.data
# filter to have at least 2 variants
wanted.pathways.names <- row.names(wanted.pathways.data)
pois.p.value.matrix = foreach(i=1:simulation_time, .combine=cbind) %dopar% {
  simulation.i <- readRDS(paste0("simulation.result/random.seed.", i, ".RDS"))
  simulation.i$pois.p.cls3[match(wanted.pathways.names, row.names(simulation.i))]
}
saveRDS(pois.p.value.matrix, file = "data/pois.p.value.matrix.RDS")
pois.p.value.matrix <- readRDS("data/pois.p.value.matrix.RDS")
corrected.pois.p.value.matrix <- foreach(i=1:dim(wanted.pathways.data)[1], .combine=rbind) %dopar% {
  min(1, sum(pois.p.value.matrix<=wanted.pathways.data$pois.p.cls3[i])/simulation_time)
}
saveRDS(corrected.pois.p.value.matrix, file = "data/corrected.pois.p.value.matrix.RDS")
stopCluster(cl)
## load
corrected.pois.p.value.matrix <- readRDS("data/corrected.pois.p.value.matrix.RDS")
wanted.pathways.data$pois.FWER.cls3 <- corrected.pois.p.value.matrix
wanted.pathways.data <- wanted.pathways.data[order(wanted.pathways.data$pois.FWER.cls3),]
# write output
dir.create('pathway.enrichment.results/')
write.table(wanted.pathways.data, file = "pathway.enrichment.results/wanted.pathway.enrichment.tsv", sep = "\t", quote = F)
wanted.pathways.data.significant <- wanted.pathways.data[wanted.pathways.data$pois.FWER.cls3<=0.05,]
write.table(wanted.pathways.data[wanted.pathways.data$pois.FWER.cls3<=0.05,],
            file = "pathway.enrichment.results/wanted.pathway.enrichment.significant.tsv", sep = "\t", quote = F)

# load data
wanted.pathways.data.significant <- read.table("pathway.enrichment.results/wanted.pathway.enrichment.significant.tsv", sep = "\t", header = T)
# output gtf file
pathwayname = 'c5.all'
pathways.tmp = as.list(GSA::GSA.read.gmt(paste("data/genesets/", pathwayname, ".v7.2.symbols.gmt", sep = "")))
wanted.pathways.data.significant.gtf <- pathways.tmp
wanted.pathways.data.significant.gtf$geneset.names <- pathways.tmp$geneset.names[pathways.tmp$geneset.names %in% rownames(wanted.pathways.data.significant)]
wanted.pathways.data.significant.gtf$genesets <- pathways.tmp$genesets[pathways.tmp$geneset.names %in% rownames(wanted.pathways.data.significant)]
wanted.pathways.data.significant.gtf$geneset.descriptions <- pathways.tmp$geneset.descriptions[pathways.tmp$geneset.names %in% rownames(wanted.pathways.data.significant)]
to.gtf <- data.frame(names=wanted.pathways.data.significant.gtf$geneset.names,
                     descriptions=wanted.pathways.data.significant.gtf$geneset.descriptions,
                     symbols=rep(NA,length(wanted.pathways.data.significant.gtf$geneset.names)))
for (i in 1:length(wanted.pathways.data.significant.gtf$genesets)) {
  to.gtf$symbols[i] <- paste(wanted.pathways.data.significant.gtf$genesets[[i]],collapse = " ")
}
write.table(to.gtf, file = "data/genesets/significant.pathways.symbols.gmt", quote = F, row.names = F, col.names = F, sep = "\t")
# do a clustering analysis and visualization don't need it for now.
get.pathway.corr.by.mutrate <- function(pathway.genes, reference) {
  result <- matrix(0, length(pathway.genes), length(pathway.genes))
  for (i in 1:(length(pathway.genes)-1)) {
    for (j in (i+1):length(pathway.genes)) {
      i.genes <- pathway.genes[[i]]
      j.genes <- pathway.genes[[j]]
      all.genes <- as.character(unique(c(i.genes, j.genes)))
      all.genes <- all.genes[all.genes%in%reference$EntrezID]
      overlap.genes <- i.genes[i.genes %in% j.genes]
      overlap.genes <- overlap.genes[overlap.genes %in% reference$EntrezID]
      all.mutrate <- sum(reference$Mu_Missense[match(all.genes, reference$EntrezID)]) +sum(reference$Mu_LoF[match(all.genes, reference$EntrezID)])
      overlap.mutrate <- sum(reference$Mu_Missense[match(overlap.genes, reference$EntrezID)]) +sum(reference$Mu_LoF[match(overlap.genes, reference$EntrezID)])
      result[i, j] <- overlap.mutrate / all.mutrate
      result[j, i] <- overlap.mutrate / all.mutrate
      }
  }
  result
}
get.pathway.corr.by.simulation <- function(pathway.name, simulation.matrix) {
  significant.pathways.corr <- cor(simulation.matrix)
  diag(significant.pathways.corr) <- 0
  colnames(significant.pathways.corr) <- pathway.name
  rownames(significant.pathways.corr) <- pathway.name
  p.matrix <- significant.pathways.corr
  validation <- significant.pathways.corr
  for (i in 1:(dim(simulation.matrix)[2]-1)) {
    for (j in (i+1):dim(simulation.matrix)[2]) {
      cortest <- cor.test(simulation.matrix[,i], simulation.matrix[,j])
      p.matrix[i, j] = cortest$p.value
      p.matrix[j, i] = cortest$p.value
      validation[i, j] = cortest$estimate
      validation[j, i] = cortest$estimate
    }
  }
  result <- list(corr=significant.pathways.corr, p=p.matrix, validation=validation)
  result
}

get.pathway.size.mutrate <- function(pathway.genes, reference) {
  result <- rep(0, length(pathway.genes))
  for (i in 1:(length(pathway.genes))) {
      i.genes <- pathway.genes[[i]]
      i.genes <- i.genes[i.genes %in% reference$EntrezID]
      i.mutrate <- sum(reference$Mu_Missense[match(i.genes, reference$EntrezID)]) +sum(reference$Mu_LoF[match(i.genes, reference$EntrezID)])
      result[i] <- i.mutrate
  }
  result
}

get.pathway.size.dnvs <- function(pathway.genes, pathways.data) {
  result <- rep(0, length(pathway.genes))
  for (i in 1:(length(pathway.genes))) {
    result[i] <- pathways.data$dn.cls3[i]
  }
  result
}

get.pathway.size.genes <- function(pathway.genes, pathways.data) {
  result <- rep(0, length(pathway.genes))
  for (i in 1:(length(pathway.genes))) {
    result[i] <- length(unique(c(strsplit(pathways.data$genes.cls1[i], split = ", ")[[1]],
                                 strsplit(pathways.data$genes.cls2[i], split = ", ")[[1]])))
  }
  result
}

significant.pathways.simulation <- pois.p.value.matrix[match(row.names(wanted.pathways.data.significant),
                                                             wanted.pathways.names),]
significant.pathways.genes <- wanted.pathways$genesets[match(row.names(wanted.pathways.data.significant),
                                                             wanted.pathways$geneset.names)]
reference = read.table("~/Data/EA.TEF.mutrate.3mer.txt", sep = "\t", header = TRUE, na.strings = ".")

significant.pathways.corr.simulation <- get.pathway.corr.by.simulation(row.names(wanted.pathways.data.significant), t(significant.pathways.simulation))
significant.pathways.corr.mutrate <- get.pathway.corr.by.mutrate(significant.pathways.genes, reference)
significant.pathways.size.mutrate <- get.pathway.size.mutrate(significant.pathways.genes, reference)
significant.pathways.size.dnvs <- get.pathway.size.dnvs(significant.pathways.genes, wanted.pathways.data.significant)
significant.pathways.size.genes <- get.pathway.size.genes(significant.pathways.genes, wanted.pathways.data.significant)


significant.pathways.names <- row.names(wanted.pathways.data.significant)
significant.pathways <- wanted.pathways
significant.pathways$genesets <- wanted.pathways$genesets[match(significant.pathways.names, wanted.pathways$geneset.names)]
significant.pathways$geneset.names <- wanted.pathways$geneset.names[match(significant.pathways.names, wanted.pathways$geneset.names)]
significant.pathways$geneset.descriptions <- wanted.pathways$geneset.descriptions[match(significant.pathways.names, wanted.pathways$geneset.names)]


to.cytoscape.edge <- matrix(NA, nrow = length(significant.pathways.names)*(length(significant.pathways.names)-1)/2, ncol=4)
colnames(to.cytoscape.edge) <- c("source", "target", "Jaccard Index", "correlation")
significant.pathways.names <- tolower(significant.pathways.names)
significant.pathways.names <- gsub("^go_", "GO: ", significant.pathways.names)
significant.pathways.names <- gsub("^hp_", "HPO: ", significant.pathways.names)
significant.pathways.names <- gsub("_", " ", significant.pathways.names)
k=1
for (i in 1:(length(significant.pathways.names)-1)) {
  for (j in (i+1):length(significant.pathways.names)) {
    to.cytoscape.edge[k, 1] = significant.pathways.names[i]
    to.cytoscape.edge[k, 2] = significant.pathways.names[j]
    to.cytoscape.edge[k, 4] = significant.pathways.corr.simulation[i, j]
    to.cytoscape.edge[k, 3] = significant.pathways.corr.mutrate[i, j]
    k=k+1
  }
}
to.cytoscape.node <- data.frame(name=significant.pathways.names,
                                mutrate=significant.pathways.size.mutrate,
                                DNVs=significant.pathways.size.dnvs,
                                genes=significant.pathways.size.genes,
                                FWER=wanted.pathways.data.significant$pois.FWER.cls3)
write.table(to.cytoscape.edge, file = "to.cytoscape.edge.tsv", sep = "\t", row.names = F, quote = F)
write.table(to.cytoscape.node, file = "to.cytoscape.node.tsv", sep = "\t", row.names = F, quote = F)



