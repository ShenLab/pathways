# load cases
mydata <- read.table("data/06.denovo.annotated.txt", sep = "\t", header = T, na.strings = ".") # your data here
samplenumber=185
# load controls
NDD_control <- read.table("./data/NDD_control_DNVs_Anno_wSpliceAI.b37.txt",
                          sep = "\t", na.strings = ".", header = T)
for (i in 1:dim(NDD_control)[1]) {
  NDD_control$DS[i] <- max(as.numeric(varhandle::unfactor(NDD_control[i,c("DS_AG","DS_AL","DS_DG","DS_DL")])))
}
NDD_control$spliceAI_Damage <- NDD_control$DS >= 0.8
NDD_control_trios <- read.table("./data/SPARK_SSC.control.trios.txt", sep = "\t", header = T)
NDD_control_trios <- length(unique(NDD_control_trios$IID))
# NDD vclass
source("./annotate.var.class.R")
NDD_control <- annotate.var.class(NDD_control, annotate.splicing = T)
NDD_control <- annotate.Dmis(NDD_control, "VEST3", 0.85)
mydata <- annotate.Dmis(mydata, "VEST3", 0.85)
# set cases
cases <- mydata[!is.na(mydata$vclass),]
NDD_control <- NDD_control[!is.na(NDD_control$vclass),]
reference = read.table("./data/EA.TEF.mutrate.3mer.txt", sep = "\t", header = TRUE, na.strings = ".")
blacklist <- read.table("./data/GENCODEV19_blacklist.txt", sep = "\t", header = F)
reference <- reference[!reference$GeneID %in% blacklist$V1,]
geneset = as.character(unique(reference$GeneID))
# pathways are downloaded from gsea gene sets: http://www.gsea-msigdb.org/gsea/msigdb/index.jsp
# I will try C3 first
# for future flexibility, the input should be list
source("pathway_burden_analysis.R")
# corrected_control_samplenumber <- floor(sum(NDD_control$vclass=="syn"&NDD_control$Symbol%in%geneset)
#                                         /sum(cases$vclass=="syn"&cases$Symbol%in%geneset)*samplenumber)
sequenced.pheno <- read.table("data/sequenced.pheno.txt", sep = "\t", header = T)
sequenced.pheno <- sequenced.pheno[!sequenced.pheno$Study.ID.on.list %in% c("CARE4-11", "R0035450", "R0035951", "R0045232"),]
complex.IID <- as.character(sequenced.pheno$Study.ID.on.list[sequenced.pheno$Syndromic.Non.syndromic=="Syndromic"])
samplenumber_complex <- length(complex.IID)
NDD_control_mis <- NDD_control
cases_comp_mis <- varhandle::unfactor(cases[cases$Syndromic.Non.syndromic=="Syndromic",])
cases_comp_mis$vclass[cases_comp_mis$vclass=="mis"]="Dmis"
NDD_control_mis$vclass[NDD_control_mis$vclass=="mis"]="Dmis"
# corrected_control_samplenumber <- floor(sum(NDD_control_mis$vclass=="syn"&NDD_control_mis$Symbol%in%geneset)
#                                         /sum(cases_comp_mis$vclass=="syn"&cases_comp_mis$Symbol%in%geneset)*samplenumber_complex)
samplenumber_control <- NDD_control_trios
pathwayname = 'c5.all'
pathways = as.list(GSA::GSA.read.gmt(paste("data/genesets/", pathwayname, ".v7.2.entrez.gmt", sep = "")))

rescale.factor <- get.mytada.data.case.control.mutrate(geneset, cases, NDD_control_mis, reference)
rescaled.reference <- reference
rescaled.reference[,5:dim(reference)[2]] = reference[,5:dim(reference)[2]] * sum(rescale.factor$dn.cls0) / 185 / sum(rescale.factor$mut.cls0) / 2

mytada.data <- get.mytada.data.case.control.mutrate(geneset, cases_comp_mis, NDD_control_mis, rescaled.reference)

pathways.data <- get.pathway.data.case.control.mutrate(mytada.data, pathways, samplenumber_complex, NDD_control_trios, rescaled.reference)

# filter pathways
min_expect <- 2
filtered.pathways <- pathways
filtered.pathways$genesets <- pathways$genesets[pathways.data$mut.cls3 >= min_expect/samplenumber_complex/2]
filtered.pathways$geneset.names <- pathways$geneset.names[pathways.data$mut.cls3 >= min_expect/samplenumber_complex/2]
filtered.pathways$geneset.descriptions <- pathways$geneset.descriptions[pathways.data$mut.cls3 >= min_expect/samplenumber_complex/2]
filtered.pathways.data <- pathways.data[pathways.data$mut.cls3 >= min_expect/samplenumber_complex/2,]
pathways <- filtered.pathways
pathways.data <- filtered.pathways.data
# save interm data to do simulation
save(mytada.data, pathways, samplenumber_complex, samplenumber_control, pathways.data, rescaled.reference, file = "data/simulation.Rdata")
# prepare for parallel simulation results folder
dir.create('simulation.result/')

