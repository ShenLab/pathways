annotate.var.class <- function(DNV_table, annotate.splicing = FALSE) {
  effs = as.character(DNV_table$GeneEff)
  standards = c("frameshift"="LGD",
                "stop_gained"="LGD",
                "missense"="mis",
                "splice_region"=NA,
                "synonymous"="syn",
                "splice_acceptor"="LGD",
                "splice_donor"="LGD",
                "start_lost"="LGD",
                "inframe_deletion"="inframe_deletion",
                "inframe_insertion"="inframe_insertion",
                "synonymous;missense"="mis",
                "missense;missense"="mis",
                "stop_lost"="LGD",
                "protein_altering"="LGD",
                "stop_retained"="syn",
                "intronic"=NA)
  standard_var = effs
  if (!annotate.splicing) {
    for (i in 1:length(effs)) {
      standard_var[i] = standards[strsplit(effs[i], split = ";")[[1]][1]]
    }
  } else {
    for (i in 1:length(effs)) {
      if (DNV_table$spliceAI_Damage[i]&!is.na(DNV_table$spliceAI_Damage[i])) {
        standard_var[i] = "LGD"
      } else {
        standard_var[i] = standards[strsplit(effs[i], split = ";")[[1]][1]]
      }
    }
  }
  DNV_table$vclass <- standard_var
  DNV_table
}

annotate.Dmis <- function(DNV_table, score_name, cutoff) {
  original_vclass <- DNV_table$vclass
  # wipe out last annotation
  original_vclass[original_vclass=="Dmis"]="mis"
  # use new Dmis annotation
  original_vclass[!is.na(DNV_table$vclass) & DNV_table$vclass %in% c("mis","Dmis") & !is.na(DNV_table[,score_name]) & DNV_table[, score_name]>=cutoff] = "Dmis"
  DNV_table$vclass <- original_vclass
  DNV_table
}

extract.scores.only <- function(DNV_table) {
  topick <- !is.na(DNV_table$vclass)
  result = data.frame(genes = as.character(DNV_table$Symbol[topick]),
                      vclass = as.character(DNV_table$vclass[topick]),
                      REVEL = as.numeric(paste(DNV_table$REVEL[topick])),
                      MPC = as.numeric(paste(DNV_table$MPC[topick])),
                      MVP2 = as.numeric(paste(DNV_table$MVP2[topick])),
                      CADD = as.numeric(paste(DNV_table$CADD13[topick])),
                      PrimateAI = as.numeric(paste(DNV_table$PrimateAI[topick])),
                      Polyphen2 = as.numeric(paste(DNV_table$Polyphen2[topick])),
                      # SIFT = as.numeric(paste(DNV_table$SIFT[topick])),
                      VEST = as.numeric(paste(DNV_table$VEST3[topick])),
                      MCAP = as.numeric(paste(DNV_table$MCAP[topick])),
                      EIGEN = as.numeric(paste(DNV_table$EIGEN[topick])),
                      GERP = as.numeric(paste(DNV_table$GERP[topick])),
                      gnomAD_oe = as.numeric(paste(DNV_table$LoFOvsE[topick]))
                      # gnomAD_pLI = as.numeric(paste(DNV_table$gnomADpLI[topick]))
                      )
  result
}