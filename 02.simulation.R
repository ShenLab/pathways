load("data/simulation.Rdata")
source("pathway_burden_analysis.R")
# do simulations to identify effect test sizes
# first generate enough bootstraping genes
# set seed just incase further replication
randomize.get.pathway.data <- function(i, mytada.data, pathways, samplenumber_case, samplenumber_control, reference) {
  set.seed(i)
  mytada.data.rand <- mytada.data
  mytada.data.rand$dn.cls0 <- rpois(dim(mytada.data)[1], mytada.data$mut.cls0*samplenumber_case*2)
  mytada.data.rand$dn.cls1 <- rpois(dim(mytada.data)[1], mytada.data$mut.cls1*samplenumber_case*2)
  mytada.data.rand$dn.cls2 <- rpois(dim(mytada.data)[1], mytada.data$mut.cls2*samplenumber_case*2)
  result <- get.pathway.data.case.control.mutrate(mytada.data.rand,
                                                  pathways,
                                                  samplenumber_case, samplenumber_control,
                                                  reference)
  result
}
args = commandArgs(trailingOnly=TRUE)
i <- as.integer(args[1])
simulation <- randomize.get.pathway.data(i, mytada.data, pathways, samplenumber_complex, samplenumber_control, rescaled.reference)
saveRDS(simulation, file=paste0("simulation.result/random.seed.", i, ".RDS"))

