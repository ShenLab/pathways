pathway_burden_analysis_mutrate <- function(geneset, pathways, cases, samplenumber, reference) {
  # form a proper input for TADA
  mytada.data = data.frame(gene.id=geneset,
                            mut.cls0=numeric(length(geneset)),
                            mut.cls1=numeric(length(geneset)),
                            mut.cls2=numeric(length(geneset)),
                            dn.cls0=numeric(length(geneset)),
                            dn.cls1=numeric(length(geneset)),
                            dn.cls2=numeric(length(geneset)))

  for (i in 1:length(geneset)) {
    index = match(as.character(geneset[i]), reference$HGNC)
    # check whether we found the gene
    if (!is.na(index)) {
      mytada.data$mut.cls0[i] = reference$Mu_Silent[index]
      mytada.data$mut.cls1[i] = reference$Mu_LoF[index]
      mytada.data$mut.cls2[i] = reference$Mu_Dmis_REVEL0.5[index]
      
    } else {
      mytada.data$mut.cls0[i] = .Machine$double.eps
      mytada.data$mut.cls1[i] = .Machine$double.eps
      mytada.data$mut.cls2[i] = .Machine$double.eps
      
    }
    if (mytada.data$mut.cls0[i] <= 0) {
      mytada.data$mut.cls0[i] = .Machine$double.eps
    }
    if (mytada.data$mut.cls1[i] <= 0) {
      mytada.data$mut.cls1[i] = .Machine$double.eps
    }
    if (mytada.data$mut.cls2[i] <= 0) {
      mytada.data$mut.cls2[i] = .Machine$double.eps
    }
    gene_cases = cases[cases$genes==as.character(geneset[i]),]
    
    mytada.data$dn.cls0[i] = length(unique(gene_cases$IID[(gene_cases$vclass=="syn")]))
    mytada.data$dn.cls1[i] = length(unique(gene_cases$IID[(gene_cases$vclass=="LGD")]))
    mytada.data$dn.cls2[i] = length(unique(gene_cases$IID[gene_cases$vclass=="mis"
                                                            & gene_cases$REVEL>=0.5
                                                            & !is.na(gene_cases$REVEL)]))
  }
  message("sanity check before correction")
  syn_oe = sum(mytada.data$dn.cls0) / (sum(mytada.data$mut.cls0) * samplenumber * 2)
  message(paste0("syn_oe = ", syn_oe))
  message("adjusting mutation rate according to syn_oe")
  # mytada.data[,2:4]=mytada.data[,2:4]*syn_oe
  
  syn_oe = sum(mytada.data$dn.cls0) / (sum(mytada.data$mut.cls0) * samplenumber * 2)
  LGD_oe = sum(mytada.data$dn.cls1) / (sum(mytada.data$mut.cls1) * samplenumber * 2)
  Dmis_oe = sum(mytada.data$dn.cls2) / (sum(mytada.data$mut.cls2) * samplenumber * 2)
  
  overall.burden = data.frame(obs_syn = sum(mytada.data$dn.cls0),
                      exp_syn = (sum(mytada.data$mut.cls0) * samplenumber * 2),
                      burden_syn = syn_oe,
                      obs_LGD = sum(mytada.data$dn.cls1),
                      exp_LGD = (sum(mytada.data$mut.cls1) * samplenumber * 2),
                      burden_LGD = LGD_oe,
                      obs_Dmis = sum(mytada.data$dn.cls2),
                      exp_Dmis = (sum(mytada.data$mut.cls2) * samplenumber * 2),
                      burden_Dmis = Dmis_oe
                      )
  
  pathways.data <- data.frame(dn.cls0 = rep(0, length(pathways$geneset.names)),
                              mut.cls0 = rep(0, length(pathways$geneset.names)),
                              burden.cls0 = rep(0, length(pathways$geneset.names)),
                              poisson.p.cls0 = rep(0, length(pathways$geneset.names)),
                              
                              dn.cls1 = rep(0, length(pathways$geneset.names)),
                              mut.cls1 = rep(0, length(pathways$geneset.names)),
                              burden.cls1 = rep(0, length(pathways$geneset.names)),
                              poisson.p.cls1 = rep(0, length(pathways$geneset.names)),
                              
                              dn.cls2 = rep(0, length(pathways$geneset.names)),
                              mut.cls2 = rep(0, length(pathways$geneset.names)),
                              burden.cls2 = rep(0, length(pathways$geneset.names)),
                              poisson.p.cls2 = rep(0, length(pathways$geneset.names)))
  for (i in 1:dim(pathways.data)[1]) {
    pathways.data$dn.cls0[i] = sum(mytada.data$dn.cls0[mytada.data$gene.id %in% unlist(pathways$genesets[i])])
    pathways.data$dn.cls1[i] = sum(mytada.data$dn.cls1[mytada.data$gene.id %in% unlist(pathways$genesets[i])])
    pathways.data$dn.cls2[i] = sum(mytada.data$dn.cls2[mytada.data$gene.id %in% unlist(pathways$genesets[i])])
    pathways.data$mut.cls0[i] = sum(mytada.data$mut.cls0[mytada.data$gene.id %in% unlist(pathways$genesets[i])])
    pathways.data$mut.cls1[i] = sum(mytada.data$mut.cls1[mytada.data$gene.id %in% unlist(pathways$genesets[i])])
    pathways.data$mut.cls2[i] = sum(mytada.data$mut.cls2[mytada.data$gene.id %in% unlist(pathways$genesets[i])])
  }
  pathways.data$burden.cls0 <- pathways.data$dn.cls0 / pathways.data$mut.cls0 / samplenumber / 2
  pathways.data$burden.cls1 <- pathways.data$dn.cls1 / pathways.data$mut.cls1 / samplenumber / 2
  pathways.data$burden.cls2 <- pathways.data$dn.cls2 / pathways.data$mut.cls2 / samplenumber / 2
  
  pathways.data$poisson.p.cls0 <- ppois(pathways.data$dn.cls0, pathways.data$mut.cls0*samplenumber*2, lower.tail = F)
  pathways.data$poisson.p.cls1 <- ppois(pathways.data$dn.cls1, pathways.data$mut.cls1*samplenumber*2, lower.tail = F)
  pathways.data$poisson.p.cls2 <- ppois(pathways.data$dn.cls2, pathways.data$mut.cls2*samplenumber*2, lower.tail = F)
  
  row.names(pathways.data) <- pathways$geneset.names
  
  pathways.data$description <- pathways$geneset.descriptions
  
  result <- list(overall.burden, pathways.data)
  result
}

get.mytada.data.case.control.mutrate <- function(geneset, cases, controls, reference) {
  mytada.data = data.frame(gene.id=geneset,
                           mut.cls0=numeric(length(geneset)),
                           mut.cls1=numeric(length(geneset)),
                           mut.cls2=numeric(length(geneset)),
                           ctl.cls0=numeric(length(geneset)),
                           ctl.cls1=numeric(length(geneset)),
                           ctl.cls2=numeric(length(geneset)),
                           dn.cls0=numeric(length(geneset)),
                           dn.cls1=numeric(length(geneset)),
                           dn.cls2=numeric(length(geneset)),
                           entrez.id=reference$EntrezID[match(geneset, reference$GeneID)])
  
  for (i in 1:length(geneset)) {
    # check whether we found the gene
    index = match(as.character(geneset[i]), reference$GeneID)
    # get mut rate
    if (!is.na(index)) {
      mytada.data$mut.cls0[i] = reference$Mu_Silent[index]
      mytada.data$mut.cls1[i] = reference$Mu_LoF[index]
      mytada.data$mut.cls2[i] = reference$Mu_Missense[index]
    } else {
      mytada.data$mut.cls0[i] = .Machine$double.eps
      mytada.data$mut.cls1[i] = .Machine$double.eps
      mytada.data$mut.cls2[i] = .Machine$double.eps
    }
    # get control, note here Dmis just mean mis
    gene_controls = controls[grep(as.character(geneset[i]), controls$GeneID),]
    mytada.data$ctl.cls0[i] = sum(gene_controls$vclass=="syn")
    mytada.data$ctl.cls1[i] = sum(gene_controls$vclass=="LGD")
    mytada.data$ctl.cls2[i] = sum(gene_controls$vclass=="Dmis")
    # get case, note here Dmis just mean mis
    gene_cases = cases[grep(as.character(geneset[i]), cases$GeneID),]
    mytada.data$dn.cls0[i] = sum(gene_cases$vclass=="syn")
    mytada.data$dn.cls1[i] = sum(gene_cases$vclass=="LGD")
    mytada.data$dn.cls2[i] = sum(gene_cases$vclass=="Dmis")
  }
  mytada.data
}

get.pathway.data.case.control.mutrate <- function(mytada.data, pathways, samplenumber_case, samplenumber_control, reference) {
  pathways.data <- data.frame(dn.cls0 = rep(0, length(pathways$geneset.names)),
                              ctl.cls0 = rep(0, length(pathways$geneset.names)),
                              mut.cls0 = rep(0, length(pathways$geneset.names)),
                              burden.cls0 = rep(0, length(pathways$geneset.names)),
                              binom.p.cls0 = rep(0, length(pathways$geneset.names)),
                              pois.p.cls0 = rep(0, length(pathways$geneset.names)),
                              genes.cls0 = rep(NA, length(pathways$geneset.names)),
                              
                              dn.cls1 = rep(0, length(pathways$geneset.names)),
                              ctl.cls1 = rep(0, length(pathways$geneset.names)),
                              mut.cls1 = rep(0, length(pathways$geneset.names)),
                              burden.cls1 = rep(0, length(pathways$geneset.names)),
                              binom.p.cls1 = rep(0, length(pathways$geneset.names)),
                              pois.p.cls1 = rep(0, length(pathways$geneset.names)),
                              genes.cls1 = rep(NA, length(pathways$geneset.names)),
                              
                              dn.cls2 = rep(0, length(pathways$geneset.names)),
                              ctl.cls2 = rep(0, length(pathways$geneset.names)),
                              mut.cls2 = rep(0, length(pathways$geneset.names)),
                              burden.cls2 = rep(0, length(pathways$geneset.names)),
                              binom.p.cls2 = rep(0, length(pathways$geneset.names)),
                              pois.p.cls2 = rep(0, length(pathways$geneset.names)),
                              genes.cls2 = rep(NA, length(pathways$geneset.names)),
                              
                              dn.cls3 = rep(0, length(pathways$geneset.names)),
                              ctl.cls3 = rep(0, length(pathways$geneset.names)),
                              mut.cls3 = rep(0, length(pathways$geneset.names)),
                              burden.cls3 = rep(0, length(pathways$geneset.names)),
                              binom.p.cls3 = rep(0, length(pathways$geneset.names)),
                              pois.p.cls3 = rep(0, length(pathways$geneset.names)),
                              genes.cls3 = rep(NA, length(pathways$geneset.names)))
  for (i in 1:dim(pathways.data)[1]) {
    for (j in 1:3) {
      pathways.data[i, j*7-6] = sum(mytada.data[,7+j][mytada.data$entrez.id %in% unlist(pathways$genesets[i])])
      pathways.data[i, j*7-5] = sum(mytada.data[,4+j][mytada.data$entrez.id %in% unlist(pathways$genesets[i])])
      pathways.data[i, j*7-4] = sum(mytada.data[,1+j][mytada.data$entrez.id %in% unlist(pathways$genesets[i])])
      
      dnv.genes <- as.character(mytada.data$gene.id[mytada.data[,7+j]!=0])
      dnv.genes.entrez.id <- as.character(mytada.data$entrez.id[mytada.data[,7+j]!=0])
      dnv.genes <- dnv.genes[dnv.genes.entrez.id %in% unlist(pathways$genesets[i])]
      if (length(dnv.genes)!=0) {
        pathways.data[i, j*7] = toString(as.character(dnv.genes))
      }
    }
  }
  pathways.data$dn.cls3 = pathways.data$dn.cls1 + pathways.data$dn.cls2
  pathways.data$mut.cls3 = pathways.data$mut.cls1 + pathways.data$mut.cls2
  pathways.data$ctl.cls3 = pathways.data$ctl.cls1 + pathways.data$ctl.cls2
  
  pathways.data$genes.cls3 = paste(pathways.data$genes.cls1, pathways.data$genes.cls2, sep = "|")
  pathways.data$burden.cls0 <- pathways.data$dn.cls0 / pathways.data$mut.cls0 / samplenumber_case / 2
  pathways.data$burden.cls1 <- pathways.data$dn.cls1 / pathways.data$mut.cls1 / samplenumber_case / 2
  pathways.data$burden.cls2 <- pathways.data$dn.cls2 / pathways.data$mut.cls2 / samplenumber_case / 2
  pathways.data$burden.cls3 <- pathways.data$dn.cls3 / pathways.data$mut.cls3 / samplenumber_case / 2
  binom.p.case.control <- function (dn_case, dn_ctl, samplenumber_case, samplenumber_control) {
    rho <- samplenumber_case / (samplenumber_case + samplenumber_control)
    if (dn_case==0 & dn_ctl==0) {
      p <- 1
    } else {
      p <- 1-pbinom(dn_case-1,
                    dn_case+dn_ctl,
                    rho)
    }
    p
  }
  pois.p.case.control <- function (dn_case, mut_rate, samplenumber_case) {
    lambda <- mut_rate
    if (lambda==0) {
      lambda = .Machine$double.eps
    }
    if (dn_case==0) {
      p <- 1
    } else {
      p <- 1-ppois(dn_case-1,
                   samplenumber_case*2*lambda)
    }
  }
  for (i in 1:dim(pathways.data)[1]) {
    for (j in 1:4) {
      pathways.data[i, j*7-2] <- binom.p.case.control(pathways.data[i, j*7-6], pathways.data[i, j*7-5], samplenumber_case, samplenumber_control)
      pathways.data[i, j*7-1] <- pois.p.case.control(pathways.data[i, j*7-6], pathways.data[i, j*7-4], samplenumber_case)
    }
  }
  row.names(pathways.data) <- pathways$geneset.names
  pathways.data$description <- pathways$geneset.descriptions
  message("finished pathway data table")
  gc()
  pathways.data
}

get.mytada.data.case.control <- function(geneset, cases, controls) {
  mytada.data = data.frame(gene.id=geneset,
                           mut.cls0=numeric(length(geneset)),
                           mut.cls1=numeric(length(geneset)),
                           mut.cls2=numeric(length(geneset)),
                           dn.cls0=numeric(length(geneset)),
                           dn.cls1=numeric(length(geneset)),
                           dn.cls2=numeric(length(geneset)))
  
  for (i in 1:length(geneset)) {
    # check whether we found the gene
    gene_controls = controls[controls$Symbol==as.character(geneset[i]),]
    mytada.data$mut.cls0[i] = sum(gene_controls$vclass=="syn")
    mytada.data$mut.cls1[i] = sum(gene_controls$vclass=="LGD")
    mytada.data$mut.cls2[i] = sum(gene_controls$vclass=="Dmis")
    gene_cases = cases[cases$Symbol==as.character(geneset[i]),]
    mytada.data$dn.cls0[i] = sum(gene_cases$vclass=="syn")
    mytada.data$dn.cls1[i] = sum(gene_cases$vclass=="LGD")
    mytada.data$dn.cls2[i] = sum(gene_cases$vclass=="Dmis")
  }
  mytada.data
}

get.pathway.data.case.control <- function(mytada.data, pathways, samplenumber_case, samplenumber_control) {
  pathways.data <- data.frame(dn.cls0 = rep(0, length(pathways$geneset.names)),
                              mut.cls0 = rep(0, length(pathways$geneset.names)),
                              burden.cls0 = rep(0, length(pathways$geneset.names)),
                              binom.p.cls0 = rep(0, length(pathways$geneset.names)),
                              pois.p.cls0 = rep(0, length(pathways$geneset.names)),
                              genes.cls0 = rep(NA, length(pathways$geneset.names)),
                              
                              dn.cls1 = rep(0, length(pathways$geneset.names)),
                              mut.cls1 = rep(0, length(pathways$geneset.names)),
                              burden.cls1 = rep(0, length(pathways$geneset.names)),
                              binom.p.cls1 = rep(0, length(pathways$geneset.names)),
                              pois.p.cls1 = rep(0, length(pathways$geneset.names)),
                              genes.cls1 = rep(NA, length(pathways$geneset.names)),
                              
                              dn.cls2 = rep(0, length(pathways$geneset.names)),
                              mut.cls2 = rep(0, length(pathways$geneset.names)),
                              burden.cls2 = rep(0, length(pathways$geneset.names)),
                              binom.p.cls2 = rep(0, length(pathways$geneset.names)),
                              pois.p.cls2 = rep(0, length(pathways$geneset.names)),
                              genes.cls2 = rep(NA, length(pathways$geneset.names)),
                              
                              dn.cls3 = rep(0, length(pathways$geneset.names)),
                              mut.cls3 = rep(0, length(pathways$geneset.names)),
                              burden.cls3 = rep(0, length(pathways$geneset.names)),
                              binom.p.cls3 = rep(0, length(pathways$geneset.names)),
                              pois.p.cls3 = rep(0, length(pathways$geneset.names)),
                              genes.cls3 = rep(NA, length(pathways$geneset.names)))
  for (i in 1:dim(pathways.data)[1]) {
    for (j in 1:3) {
      pathways.data[i, j*6-5] = sum(mytada.data[,4+j][mytada.data$gene.id %in% unlist(pathways$genesets[i])])
      pathways.data[i, j*6-4] = sum(mytada.data[,1+j][mytada.data$gene.id %in% unlist(pathways$genesets[i])])
      dnv.genes <- as.character(mytada.data$gene.id[mytada.data[,4+j]!=0])
      dnv.genes <- dnv.genes[dnv.genes %in% unlist(pathways$genesets[i])]
      if (length(dnv.genes)!=0) {
        pathways.data[i, j*6] = toString(as.character(dnv.genes))
      }
    }
  }
  pathways.data$dn.cls3 = pathways.data$dn.cls1 + pathways.data$dn.cls2
  pathways.data$mut.cls3 = pathways.data$mut.cls1 + pathways.data$mut.cls2
  pathways.data$genes.cls3 = paste(pathways.data$genes.cls1, pathways.data$genes.cls2, sep = "|")
  pathways.data$burden.cls0 <- pathways.data$dn.cls0 / pathways.data$mut.cls0 / samplenumber_case * samplenumber_control
  pathways.data$burden.cls1 <- pathways.data$dn.cls1 / pathways.data$mut.cls1 / samplenumber_case * samplenumber_control
  pathways.data$burden.cls2 <- pathways.data$dn.cls2 / pathways.data$mut.cls2 / samplenumber_case * samplenumber_control
  pathways.data$burden.cls3 <- pathways.data$dn.cls3 / pathways.data$mut.cls3 / samplenumber_case * samplenumber_control
  binom.p.case.control <- function (dn_case, dn_ctl) {
    rho <- samplenumber_case / (samplenumber_case + samplenumber_control)
    if (dn_case==0 & dn_ctl==0) {
      p <- 1
    } else {
      p <- 1-pbinom(dn_case-1,
                    dn_case+dn_ctl,
                    rho)
    }
    p
  }
  pois.p.case.control <- function (dn_case, dn_ctl) {
    lambda <- dn_ctl/samplenumber_control/2
    if (lambda==0) {
      lambda = .Machine$double.eps
    }
    if (dn_case==0 & dn_ctl==0) {
      p <- 1
    } else {
      p <- 1-ppois(dn_case-1,
                   samplenumber_case*2*lambda)
    }
  }
  for (i in 1:dim(pathways.data)[1]) {
    for (j in 1:4) {
      pathways.data[i, j*6-2] <- binom.p.case.control(pathways.data[i, j*6-5], pathways.data[i, j*6-4])
      pathways.data[i, j*6-1] <- pois.p.case.control(pathways.data[i, j*6-5], pathways.data[i, j*6-4])
    }
  }
  row.names(pathways.data) <- pathways$geneset.names
  pathways.data$description <- pathways$geneset.descriptions
  message("finished pathway data table")
  gc()
  pathways.data
}

pathway_burden_analysis_casecontrol <- function(geneset, pathways, cases, controls, samplenumber_case, samplenumber_control, simulation_time = 10000) {
  # form a proper input, here mut does not mean true mut but just occurance in controls
  mydata.data <- get.mytada.data.case.control(geneset, cases, controls)
  message("sanity check before correction")
  syn_oe = sum(mytada.data$dn.cls0) / (sum(mytada.data$mut.cls0) * samplenumber_case / samplenumber_control)
  message(paste0("syn_oe = ", syn_oe))
  message("adjusting mutation rate according to syn_oe")
  # mytada.data[,2:4]=mytada.data[,2:4]*syn_oe
  
  syn_oe = sum(mytada.data$dn.cls0) / (sum(mytada.data$mut.cls0) * samplenumber_case / samplenumber_control)
  LGD_oe = sum(mytada.data$dn.cls1) / (sum(mytada.data$mut.cls1) * samplenumber_case / samplenumber_control)
  Dmis_oe = sum(mytada.data$dn.cls2) / (sum(mytada.data$mut.cls2) * samplenumber_case / samplenumber_control)
  
  overall.burden = data.frame(obs_syn = sum(mytada.data$dn.cls0),
                              exp_syn = (sum(mytada.data$mut.cls0) * samplenumber_case / samplenumber_control),
                              burden_syn = syn_oe,
                              obs_LGD = sum(mytada.data$dn.cls1),
                              exp_LGD = (sum(mytada.data$mut.cls1) * samplenumber_case / samplenumber_control),
                              burden_LGD = LGD_oe,
                              obs_Dmis = sum(mytada.data$dn.cls2),
                              exp_Dmis = (sum(mytada.data$mut.cls2) * samplenumber_case / samplenumber_control),
                              burden_Dmis = Dmis_oe)
  # real pathways.data
  pathways.data <- get.pathway.data.case.control(mytada.data, pathways, samplenumber_case, samplenumber_control)
  # filter pathways and keep only pathways with >= 2 variants in control
  filtered.pathways <- pathways
  filtered.pathways$genesets <- pathways$genesets[pathways.data$mut.cls3 >= 3]
  filtered.pathways$geneset.names <- pathways$geneset.names[pathways.data$mut.cls3 >= 3]
  filtered.pathways$geneset.descriptions <- pathways$geneset.descriptions[pathways.data$mut.cls3 >= 3]
  pathways.data <- pathways.data[pathways.data$mut.cls3 >= 3,]
  result <- list(overall.burden, pathways.data)
  gc()
  result
}

load_pathway_output_burden_file <- function(pathwayname) {
  pathways = as.list(GSA::GSA.read.gmt(paste("genesets/", pathwayname, ".v7.2.symbols.gmt", sep = "")))
  message(paste0("pathway ", pathwayname, " has ", length(pathways$geneset.names), " genesets"))
  message(paste0("begin pathway ", pathwayname, " all Dmis"))
  pathways_burden <- pathway_burden_analysis_casecontrol(geneset, pathways,
                                                         cases, NDD_control,
                                                         samplenumber,
                                                         corrected_control_samplenumber)
  saveRDS(pathways_burden, file = paste("pathway.enrichment.results.Dmis/", pathwayname,
                                        ".pathway.result.RDS", sep = ""))
  message(paste0("finished pathway ", pathwayname, " all Dmis"))
  pathway.data <- pathways_burden[[2]]
  pathway.data <- pathway.data[order(pathway.data$binom.p.cls3),]
  write.table(pathway.data, file = paste("pathway.enrichment.results.Dmis/", pathwayname,
                                         ".pathway.enrichment.tsv", sep = ""),
              quote = F, row.names = T, col.names = T, sep = "\t", na = ".")
  rm(pathways_burden, pathway.data)
  gc()
  message(paste0("begin pathway ", pathwayname, " complex mis"))
  pathways_burden <- pathway_burden_analysis_casecontrol(geneset, pathways,
                                                         cases_mis, NDD_control_mis,
                                                         samplenumber_complex,
                                                         corrected_control_samplenumber)
  saveRDS(pathways_burden, file = paste("pathway.enrichment.results/", pathwayname,
                                        ".pathway.result.RDS", sep = ""))
  message(paste0("finished pathway ", pathwayname, " complex mis"))
  pathway.data <- pathways_burden[[2]]
  pathway.data <- pathway.data[order(pathway.data$binom.p.cls3),]
  write.table(pathway.data, file = paste("pathway.enrichment.results/", pathwayname,
                                         ".pathway.enrichment.tsv", sep = ""),
              quote = F, row.names = T, col.names = T, sep = "\t", na = ".")
  rm(pathways_burden, pathway.data)
  gc()
}
