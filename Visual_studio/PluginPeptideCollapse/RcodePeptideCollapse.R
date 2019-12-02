##################################################
# Project: DIA consolidation
# Author(s): Alexander Hogrebe
# Date: 2019-04-03
# Version: 1.3
# Script purpose: Perseus plugin scripting
##
##
##################################################

#read Perseus command-line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("Perseus should provide three arguments: parameters inFile outFile", call.=FALSE)
}
paramFile <- args[1]
inFile <- args[2]
outFile <- args[3]

# test settings
# paramFile <- "C:/Users/sxk864/AppData/Local/Temp/tmp95A9.tmp"
# inFile <- "C:/Users/sxk864/AppData/Local/Temp/tmp87D1.tmp"
# outFile <- "C:/Users/sxk864/AppData/Local/Temp/tmp88EB.tmp"
# rm(list=(grep("^par\\..*$", ls(), value=T, perl=T)))

# test settings
# paramFile <- "C:/Users/sxk864/AppData/Local/Temp/tmp5919.tmp"
# inFile <- "C:/Users/sxk864/AppData/Local/Temp/tmp5436.tmp"
# outFile <- "C:/Users/sxk864/AppData/Local/Temp/tmp5966.tmp"
# rm(list=(grep("^par\\..*$", ls(), value=T, perl=T)))

#installing packages if not installed yet
options(repos='http://cran.rstudio.com/') #without this option, packages cannot be installed on a fresh R installation
list.of.packages <- c("PerseusR", "data.table", "XML", "doParallel", "stringr", "MASS", "BiocManager", "Rcpp")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#load packages
suppressMessages(library(data.table))
suppressMessages(library(XML))
suppressMessages(library(doParallel))
suppressMessages(library(stringr))
suppressMessages(library(MASS))
suppressMessages(library(BiocManager))

#install Biostrings for FASTA readin
list.of.packages <- c("Biostrings", "Biobase")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages, ask = F)

#load Biostrings
suppressMessages(library(Biostrings))
suppressMessages(library(PerseusR)) #requires Biobase to work properly


#str(parameters)
#Perseus Data Read
mdata <- read.perseus(inFile)
parameters <- xmlToList(paramFile)
#extract parameters from xml file
#condition column renaming; 0 = grouping by condition column, 1 = wide-format no grouping
par.cond <- as.integer(parameters[[1]][[1]][[1]])
if(par.cond == 0){
  #read out column name
  par.cond.col <- unlist(parameters[[1]][[1]][[3]][[1]][[1]][[1]][[2]])[as.integer(parameters[[1]][[1]][[3]][[1]][[1]][[1]][[1]])+1]
  
  #rewrite column name using . for spaces and brackets
  par.cond.col <- gsub("( |\\(|\\)|\\[|\\])", ".", par.cond.col, perl=T)
}
# if(par.cond == 1){
#   #no additional input read out
# }

#collapse level; 0 = localized site level, 1 = localized peptide level (stoich), 2 = ModSpec peptide level
par.level <- as.integer(parameters[[1]][[2]][[1]])
#read out info for 0 = localized site level
if(par.level == 0){
  #read out probability column type; 0 = EG.PTMLocalizationProbabilities, 1 = EG.PTMAssayProbability, 2 = No probability column
  par.level.col <- as.integer(parameters[[1]][[2]][[3]][[1]][[1]][[1]][[1]])
  
  if(par.level.col == 0){
    par.level.cutoff <- as.numeric(parameters[[1]][[2]][[3]][[1]][[1]][[1]][[3]][[1]][[1]][[1]][[1]])
  }
  if(par.level.col == 1){
    par.level.cutoff <- as.numeric(parameters[[1]][[2]][[3]][[1]][[1]][[1]][[3]][[2]][[1]][[1]][[1]])
  }
  if(par.level.col == 2){
    # par.level.col.nocut <- unlist(parameters[[1]][[2]][[3]][[1]][[1]][[1]][[3]][[3]][[1]][[1]][[2]])[
    #   as.integer(parameters[[1]][[2]][[3]][[1]][[1]][[1]][[3]][[3]][[1]][[1]][[1]])+1]
    par.level.col.nocut <- as.integer(parameters[[1]][[2]][[3]][[1]][[1]][[1]][[3]][[3]][[1]][[1]][[1]])
  }
  
  #read out site level column type; 0 = PG.Genes, 1 = PG.ProteinGroups
  par.level.col.genprot <- unlist(parameters[[1]][[2]][[3]][[1]][[1]][[2]][[2]])[as.integer(parameters[[1]][[2]][[3]][[1]][[1]][[2]][[1]])+1]
  
  #read out FASTA file path
  par.FASTA.file <- gsub("\\", "/", parameters[[1]][[2]][[3]][[1]][[1]][[3]][[1]], fixed=T)
  
  #read out FASTA parsing rule
  par.FASTA.parse <- parameters[[1]][[2]][[3]][[1]][[1]][[4]][[1]]
}

#read out info for 1 = localized peptide level (stoich)
if(par.level == 1){
  #read out probability column type; 0 = EG.PTMLocalizationProbabilities, 1 = EG.PTMAssayProbability, 2 = No probability column
  par.level.col <- as.integer(parameters[[1]][[2]][[3]][[2]][[1]][[1]][[1]])
  
  if(par.level.col == 0){
    par.level.cutoff <- as.numeric(parameters[[1]][[2]][[3]][[2]][[1]][[1]][[3]][[1]][[1]][[1]][[1]])
  }
  if(par.level.col == 1){
    par.level.cutoff <- as.numeric(parameters[[1]][[2]][[3]][[2]][[1]][[1]][[3]][[2]][[1]][[1]][[1]])
  }
  if(par.level.col == 2){
    # par.level.col.nocut <- unlist(parameters[[1]][[2]][[3]][[2]][[1]][[1]][[3]][[3]][[1]][[1]][[2]])[
    #   as.integer(parameters[[1]][[2]][[3]][[2]][[1]][[1]][[3]][[3]][[1]][[1]][[1]])+1]
    par.level.col.nocut <- as.integer(parameters[[1]][[2]][[3]][[2]][[1]][[1]][[3]][[3]][[1]][[1]][[1]])
  }
  if(par.level.col == 3){
    #MaxQuant modified sequence
    par.level.cutoff <- as.numeric(parameters[[1]][[2]][[3]][[2]][[1]][[1]][[3]][[4]][[1]][[1]][[1]])
    
    #read out MQ probability column, only used for this specific setup
    par.level.MQ.col <- unlist(parameters[[1]][[2]][[3]][[2]][[1]][[1]][[3]][[4]][[1]][[2]][[2]])[
      as.integer(parameters[[1]][[2]][[3]][[2]][[1]][[1]][[3]][[4]][[1]][[2]][[1]])+1]
    #rewrite column name using . for spaces and brackets
    par.level.MQ.col <- gsub("( |\\(|\\)|\\[|\\])", ".", par.level.MQ.col, perl=T)
  }
  
  #read out stoichiometry settings; 0 = Calculate stoichiometries, 1 = Skip
  par.level.stoich <- as.integer(parameters[[1]][[2]][[3]][[2]][[1]][[2]][[1]])
}

#read out info for 2 = ModSpec peptide level
if(par.level == 2){
  #read out probability column type; 0 = EG.PTMLocalizationProbabilities, 1 = EG.PTMAssayProbability, 2 = No probability column
  par.level.col <- as.integer(parameters[[1]][[2]][[3]][[3]][[1]][[1]][[1]])
  
  if(par.level.col == 0){
    par.level.cutoff <- as.numeric(parameters[[1]][[2]][[3]][[3]][[1]][[1]][[3]][[1]][[1]][[1]][[1]])
  }
  if(par.level.col == 1){
    par.level.cutoff <- as.numeric(parameters[[1]][[2]][[3]][[3]][[1]][[1]][[3]][[2]][[1]][[1]][[1]])
  }
  if(par.level.col == 2){
    # par.level.col.nocut <- unlist(parameters[[1]][[2]][[3]][[3]][[1]][[1]][[3]][[3]][[1]][[1]][[2]])[
    #   as.integer(parameters[[1]][[2]][[3]][[3]][[1]][[1]][[3]][[3]][[1]][[1]][[1]])+1]
    par.level.col.nocut <- as.integer(parameters[[1]][[2]][[3]][[3]][[1]][[1]][[3]][[3]][[1]][[1]][[1]])
  }
}

#read out PTM strings
par.PTM <- unlist(strsplit(parameters[[1]][[3]][[1]], split=";"))

#read out aggregation type; 0 = Linear modeling based, 1 = summing
par.agg <- as.integer(parameters[[1]][[4]][[1]])

#read out number of CPU threads
par.CPU <- as.integer(parameters[[1]][[5]][[1]])



## START DATA PROCESSING

#transform main columns and annotation column into data.table object
data <- data.table(mdata@main, mdata@annotCols)
#data <- data[, .SD[1:1000,], by=par.cond.col]


#check if condition column contains at least 2 distinct conditions
if(par.cond == 0L){
  if(data[, length(unique(get(par.cond.col)))<2]){
    stop("Condition column contains less than 2 distinct conditions, please select a different condition column or renaming parameters")
  }
}

#check if necessary columns present
if(par.level==0){
  #need PEP.PeptidePosition
  if(!"PEP.PeptidePosition" %in% names(data)){
    stop("Missing: site-level consolidation requires PEP.PeptidePosition column to work")
  }

  #need PG.Genes
  if(par.level.col.genprot==0 & !"PG.Genes" %in% names(data)){
    stop("Missing: site-level consolidation with chosen settings requires PG.Genes column to work")
  }

  #need PG.ProteinGroups
  if(par.level.col.genprot==1 & !"PG.ProteinGroups" %in% names(data)){
    stop("Missing: site-level consolidation with chosen settings requires PG.ProteinGroups column to work")
  }
}

#if MQ modified sequence cutoff set, check if at least 1 "Probability" column provided
if(par.level == 1L & par.level.col == 3L){
  if( length(grep("Probabilities$", names(data), perl=T)) == 0L){
    stop("When probability filtering with MaxQuant modified sequence input, please provide at least 1 ... Probabilities column")
  }
}


# setnames(data, "EG.TotalQuantity..Settings.", "Int1")
# data[, Int2 := Int1*2]
# main.name <- c("Int1", "Int2")
# data[, R.FileName := gsub("20171125_QE7_nLC14_DBJ_SA_DIAphos_RPE1_pilot2_", "DIA_", R.FileName, fixed=T)]


#extract main column names
main.name <- names(mdata@main)
else.name <- names(mdata@annotCols)

# #stop function if more than 1 main column provided
# if(length(main.name)>1){
#   stop("Only 1 main column may be provided for long-to-wide format transformation")
# }

#delete "Filtered", "0" and "1" entries, then reload as numeric columns
for(col in main.name){
  data[get(col)=="Filtered", paste(col) := NA]
  data[get(col)==1, paste(col) := NA]
  data[get(col)==0, paste(col) := NA]
  data[get(col)=="NaN", paste(col) := NA]
  data[, paste(col) := as.numeric(get(col))]
}
rm(col)



#correct EG.PTMLocProb and/or extract PTM_0 pos/probs from EG.PTMLocProb/EG.PrecId
#initiate multiple cores
cl <- makeCluster(par.CPU)
registerDoParallel(cl)

#for data testing only
#temp <- copy(data)
#PTM <- par.PTM[1]

if(par.level == 0){
  
  if(any(par.level.col == 0, if(exists("par.level.col.nocut")){par.level.col.nocut == 0}, na.rm = T)){
    
    #if par.level = 0 & (par.level.col = 0 OR par.level.col.nocut = 0) -> EG.PTMLocProb-based: correct EG.PrecId & extract PTM_0 positions/probs
    
    #define function
    PTM.valid.extract.site.string <- function(temp, par.PTM){
      
      PTM.count <- 0
      for(PTM in par.PTM){
        #remove square brackets from PTM and mark round brackets as escaped, then insert into pattern
        pat <- gsub("[^:]*", gsub(")", "\\)", gsub("(", "\\(", gsub("(\\[|\\])", "", PTM, perl=T), fixed=T), fixed=T), "\\[[^:]*\\: (.{1,5})\\%\\]", fixed=T)
        pat.del <- gsub("[^:]*", gsub(")", "\\)", gsub("(", "\\(", gsub("(\\[|\\])", "", PTM, perl=T), fixed=T), fixed=T), "\\[[^:]*\\: .{1,5}\\%\\]", fixed=T)
        
        #create PTM valid row vector
        PTM.exist <- paste("PTM_", PTM.count, "_exist", sep="")
        temp[, paste(PTM.exist) :=F]
        temp[like(EG.PTMLocalizationProbabilities, pat), paste(PTM.exist) :=T]
        
        #use pattern to extract total number of respective PTM probabilities summed up
        temp[, paste("PTM_", PTM.count, "_num", sep="") :=
               sapply(EG.PTMLocalizationProbabilities, function(x){ round(sum(as.numeric(gsub(pat, "\\1", unlist(str_extract_all(x, pat)), perl=T)))/100, 0) })]
        
        #use pattern to extract respective PTM probabilities
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_prob", sep="") :=
               sapply(EG.PTMLocalizationProbabilities, function(x){ as.numeric(gsub(pat, "\\1", unlist(str_extract_all(x, pat)), perl=T))/100 })]
        
        #create EG.PTMLocalizationProbabilities with other PTMs deleted
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_seq", sep="") := EG.PTMLocalizationProbabilities]
        for(PTM.other in setdiff(par.PTM, PTM)){
          pat.other <- gsub("[^:]*", gsub(")", "\\)", gsub("(", "\\(", gsub("(\\[|\\])", "", PTM.other, perl=T), fixed=T), fixed=T), "\\[[^:]*\\: .{1,5}\\%\\]", fixed=T)
          temp[, paste(paste("PTM_", PTM.count, "_seq", sep="")) := gsub(pat.other, "", get(paste("PTM_", PTM.count, "_seq", sep="")), perl=T)]
        }
        
        #extract positions of respective PTMs
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_pos", sep="") := sapply(get(paste("PTM_", PTM.count, "_seq", sep="")), function(x){
          #problem: if PTM on last/first position can create exemption if "_" removed -> use as anchor to mark start/stop of sequence and then substract
          temp <- nchar(unlist(strsplit(x, split=pat.del)))
          temp <- temp-c(1, rep(0, length(temp)-2), 1)
          return(cumsum(temp[-length(temp)]))
        })]
        #extract PTM valid localization vector: eg T, T, F, T, F
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_log", sep="") :=
               mapply(function(x, y){ rank(-x, ties.method = "first")<=y },
                      x=get(paste("PTM_", PTM.count, "_prob", sep="")), y=get(paste("PTM_", PTM.count, "_num", sep="")))]
        
        #extract PTM valid positions: eg 3,8,15
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_pos_val", sep="") :=
               mapply(function(x, y){ list(x[y]) }, x=get(paste("PTM_", PTM.count, "_pos", sep="")), y=get(paste("PTM_", PTM.count, "_log", sep="")))]
        
        #extract PTM valid probabilities: eg 1,1,1
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_prob_val", sep="") :=
               mapply(function(x, y){ list(x[y]) }, x=get(paste("PTM_", PTM.count, "_prob", sep="")), y=get(paste("PTM_", PTM.count, "_log", sep="")))]
        
        #increase PTM counter
        PTM.count <- PTM.count+1
      }
      
      #create base sequence column
      temp[, PTM_base_seq := gsub("\\[[^:]*\\: .{1,5}\\%\\]", "", gsub("_", "", EG.PTMLocalizationProbabilities, fixed=T), perl=T)]
      #create summed PTM num column
      temp[, PTM_num := apply(.SD, 1, function(x){ sum(x, na.rm=T) }), .SDcols = grep("^PTM_.*_num$", names(temp), value=T, perl=T)]
      #create combined PTM position column
      temp[PTM_num>0, PTM_pos_val := apply(.SD, 1, function(x){ unlist(x) }), .SDcols = grep("^PTM_.*_pos_val$", names(temp), value=T, perl=T)]
      #create combined PTM type column
      temp[PTM_num>0, PTM_type_val := list(apply(.SD, 1, function(x){ unlist(mapply(function(y, z){rep(z, y)}, y=x, z=par.PTM)) })),
           .SDcols = grep("^PTM_.*_num$", names(temp), value=T, perl=T)]
      #create sorted combined PTM position column
      temp[PTM_num>0, PTM_pos_val_sorted := lapply(PTM_pos_val, sort)]
      #create sorted combined PTM type column
      temp[PTM_num>0, PTM_type_val_sorted := mapply(function(x, y){ x[order(y)] }, x=PTM_type_val, y=PTM_pos_val)]
      
      #recreate EG.PrecursorId column with valid types only
      temp[PTM_num==0, EG.PrecursorId.PTM.val := PTM_base_seq]
      temp[PTM_num>0, EG.PrecursorId.PTM.val := mapply(function(posit, types, baseseq){
        #step by step create res character string
        res <- character()
        pos.save <- 0
        #iterate for loop over each individual pos/type
        for(i in 1:(length(posit)) ){
          #build res from base sequence until pos
          res <- paste(res, substr(baseseq, start=pos.save, stop=posit[i]), sep="")
          #insert PTM type
          res <- paste(res, types[i], sep="")
          #update pos.save
          pos.save <- posit[i]+1
        }
        #finish by attaching last sequence substring
        return(paste(res, substr(baseseq, start=pos.save, stop=nchar(baseseq)), sep=""))
      }, posit=PTM_pos_val_sorted, types=PTM_type_val_sorted, baseseq=PTM_base_seq)]
      
      #create column combined with charge state for dcast
      temp[, PTM_group := paste("_", EG.PrecursorId.PTM.val, "_", gsub("_.*_(\\..)$", "\\1", EG.PrecursorId, perl=T), sep="")]
      
      #delete all rows except PTM_0 and PTM_base_seq (needed for aa type extraction)
      for(col in setdiff(grep("^PTM_[1-9][0-9]*_.*$", names(temp), value=T, perl=T), c("PTM_base_seq", "PTM_group"))){
        temp[, paste(col) := NULL]
      }
      
      #rename probability column
      setnames(temp, "PTM_0_prob_val", "PTM_localization")
      
      #return data table to top level
      return(temp)
    }
    
    #filter out all rows NOT containing target PTM (-> useless for site-level collapsing)
    pat <- gsub("[^:]*", gsub(")", "\\)", gsub("(", "\\(", gsub("(\\[|\\])", "", par.PTM[1], perl=T), fixed=T), fixed=T), "\\[[^:]*\\: (.{1,5})\\%\\]", fixed=T)
    
    if(data[!like(EG.PTMLocalizationProbabilities, pat), .N]>0){
      print(paste("Cave: Removed ", data[!like(EG.PTMLocalizationProbabilities, pat), .N],
                  " rows not containin target PTM in EG.PTMLocalizationProbabilities", sep=""))
    }
    data <- copy(data[like(EG.PTMLocalizationProbabilities, pat), ])
    rm(pat)
    
    #divide dataset into sub-datasets for parallelization
    data[, subd := sort(c(rep((1:(par.CPU-1)), round((.N/par.CPU), 0)), rep(par.CPU, .N-(round((.N/par.CPU), 0)*(par.CPU-1)) )))]
    data.sub <- split(data, by="subd")
    
    data <- copy(rbindlist(foreach(testp=data.sub, .packages=c("data.table", "stringr")) %dopar% {
      PTM.valid.extract.site.string(temp=testp, par.PTM=par.PTM)
    }))
    print("Done: EG.PTMLocalizationProbabilities-based EG.PrecursorId correction and PTM position/probability extraction")
    #temp[PTM_0_num>0, .(PTM_0_pos_val, PTM_localization, EG.PrecursorId.PTM.val, EG.PrecursorId.PTM.val.charge)]
    
  } else if(any(par.level.col == 1, if(exists("par.level.col.nocut")){par.level.col.nocut == 1}, na.rm = T)){
    #if par.level = 0 & (par.level.col = 1 OR par.level.col.nocut = 1) -> EG.PrecId-based: extract PTM_0 positions
    
    #define function
    PTM.valid.extract.site.string <- function(temp, par.PTM){
      
      #extract num of PTM 0 from EG.PrecId
      temp[, PTM_0_num := sapply(EG.PrecursorId, function(x){ length(unlist(strsplit(x, split=par.PTM[1], fixed=T)))-1 })]
      
      #create base sequence column, and PTM 0 only seq column
      temp[, PTM_base_seq := gsub("^(.*)\\..$", "\\1", gsub("\\[[^[]*\\]", "", gsub("_", "", EG.PrecursorId, fixed=T), perl=T), perl=T)]
      #define deletion pattern for all PTMs except target PTM
      pat.del.all <- paste("\\[(?!", gsub(")", "\\)", gsub("(", "\\(", substr(par.PTM[1], start = 2, stop = nchar(par.PTM[1])-1), fixed=T), fixed=T),
                           ")[^[]*\\]", sep="")
      temp[, PTM_0_seq := gsub(pat.del.all, "", gsub("^(.*)\\..$", "\\1", EG.PrecursorId, perl=T), perl=T)]
      
      #create PTM_0_pos_val
      temp[PTM_0_num>0, PTM_0_pos_val := sapply(PTM_0_seq, function(x){
        #problem: if PTM on last/first position can create exemption if "_" removed -> use as anchor to mark start/stop of sequence and then substract
        tempo <- nchar(unlist(strsplit(x, split=par.PTM[1], fixed=T)))
        tempo <- tempo-c(1, rep(0, length(tempo)-2), 1)
        return(cumsum(tempo[-length(tempo)]))
      })]
      
      #delete PTM_0_seq
      temp[, PTM_0_seq := NULL]
      
      #create PTM group
      temp[, PTM_group := EG.PrecursorId]
      
      #return data table to top level
      return(temp)
    }
    
    #filter out all rows NOT containing target PTM (-> useless for site-level collapsing)
    pat <- gsub("(", "\\(", gsub(")", "\\)", gsub("[", "\\[", gsub("]", "\\]", par.PTM[1], fixed=T), fixed=T), fixed=T), fixed=T)
    if(data[!like(EG.PrecursorId, pat), .N]>0){
      print(paste("Cave: Removed ", data[!like(EG.PrecursorId, pat), .N],
                  " rows not containin target PTM in EG.PTMLocalizationProbabilities", sep=""))
    }
    data <- data[like(EG.PrecursorId, pat), ]
    rm(pat)
    
    #divide dataset into sub-datasets for parallelization
    data[, subd := sort(c(rep((1:(par.CPU-1)), round((.N/par.CPU), 0)), rep(par.CPU, .N-(round((.N/par.CPU), 0)*(par.CPU-1)) )))]
    data.sub <- split(data, by="subd")
    
    data <- copy(rbindlist(foreach(testp=data.sub, .packages=c("data.table", "stringr")) %dopar% {
      PTM.valid.extract.site.string(temp=testp, par.PTM=par.PTM)
    }))
    print("Done: EG.PrecursorId-based PTM position extraction")
    #temp[PTM_0_num>0, .(PTM_0_pos_val)]
    
  }
} else
  if(par.level == 1){
  
  if(any(par.level.col == 0, if(exists("par.level.col.nocut")){par.level.col.nocut == 0}, na.rm = T)){
    #if par.level = 1 & (par.level.col = 0 OR par.level.col.nocut = 0) -> EG.PTMLocProb-based: correct EG.PrecId
    
    #define function
    PTM.valid.extract.site.string <- function(temp, par.PTM){
      
      PTM.count <- 0
      for(PTM in par.PTM){
        #remove square brackets from PTM and mark round brackets as escaped, then insert into pattern
        pat <- gsub("[^:]*", gsub(")", "\\)", gsub("(", "\\(", gsub("(\\[|\\])", "", PTM, perl=T), fixed=T), fixed=T), "\\[[^:]*\\: (.{1,5})\\%\\]", fixed=T)
        pat.del <- gsub("[^:]*", gsub(")", "\\)", gsub("(", "\\(", gsub("(\\[|\\])", "", PTM, perl=T), fixed=T), fixed=T), "\\[[^:]*\\: .{1,5}\\%\\]", fixed=T)
        
        #create PTM valid row vector
        PTM.exist <- paste("PTM_", PTM.count, "_exist", sep="")
        temp[, paste(PTM.exist) :=F]
        temp[like(EG.PTMLocalizationProbabilities, pat), paste(PTM.exist) :=T]
        
        #use pattern to extract total number of respective PTM probabilities summed up
        temp[, paste("PTM_", PTM.count, "_num", sep="") :=
               sapply(EG.PTMLocalizationProbabilities, function(x){ round(sum(as.numeric(gsub(pat, "\\1", unlist(str_extract_all(x, pat)), perl=T)))/100, 0) })]
        
        #use pattern to extract respective PTM probabilities
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_prob", sep="") :=
               sapply(EG.PTMLocalizationProbabilities, function(x){ as.numeric(gsub(pat, "\\1", unlist(str_extract_all(x, pat)), perl=T))/100 })]
        
        #create EG.PTMLocalizationProbabilities with other PTMs deleted
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_seq", sep="") := EG.PTMLocalizationProbabilities]
        for(PTM.other in setdiff(par.PTM, PTM)){
          pat.other <- gsub("[^:]*", gsub(")", "\\)", gsub("(", "\\(", gsub("(\\[|\\])", "", PTM.other, perl=T), fixed=T), fixed=T), "\\[[^:]*\\: .{1,5}\\%\\]", fixed=T)
          temp[, paste(paste("PTM_", PTM.count, "_seq", sep="")) := gsub(pat.other, "", get(paste("PTM_", PTM.count, "_seq", sep="")), perl=T)]
        }
        
        #extract positions of respective PTMs
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_pos", sep="") := sapply(get(paste("PTM_", PTM.count, "_seq", sep="")), function(x){
          #problem: if PTM on last/first position can create exemption if "_" removed -> use as anchor to mark start/stop of sequence and then substract
          temp <- nchar(unlist(strsplit(x, split=pat.del)))
          temp <- temp-c(1, rep(0, length(temp)-2), 1)
          return(cumsum(temp[-length(temp)]))
        })]
        #extract PTM valid localization vector: eg T, T, F, T, F
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_log", sep="") :=
               mapply(function(x, y){ rank(-x, ties.method = "first")<=y },
                      x=get(paste("PTM_", PTM.count, "_prob", sep="")), y=get(paste("PTM_", PTM.count, "_num", sep="")))]
        
        #extract PTM valid positions: eg 3,8,15
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_pos_val", sep="") :=
               mapply(function(x, y){ list(x[y]) }, x=get(paste("PTM_", PTM.count, "_pos", sep="")), y=get(paste("PTM_", PTM.count, "_log", sep="")))]
        
        #extract PTM valid probabilities: eg 1,1,1
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_prob_val", sep="") :=
               mapply(function(x, y){ list(x[y]) }, x=get(paste("PTM_", PTM.count, "_prob", sep="")), y=get(paste("PTM_", PTM.count, "_log", sep="")))]
        
        #increase PTM counter
        PTM.count <- PTM.count+1
      }
      
      
      #create base sequence column
      temp[, PTM_base_seq := gsub("\\[[^:]*\\: .{1,5}\\%\\]", "", gsub("_", "", EG.PTMLocalizationProbabilities, fixed=T), perl=T)]
      #create summed PTM num column
      temp[, PTM_num := apply(.SD, 1, function(x){ sum(x, na.rm=T) }), .SDcols = grep("^PTM_.*_num$", names(temp), value=T, perl=T)]
      #create combined PTM position column
      temp[PTM_num>0, PTM_pos_val := apply(.SD, 1, function(x){ unlist(x) }), .SDcols = grep("^PTM_.*_pos_val$", names(temp), value=T, perl=T)]
      #create combined PTM type column
      temp[PTM_num>0, PTM_type_val := list(apply(.SD, 1, function(x){ unlist(mapply(function(y, z){rep(z, y)}, y=x, z=par.PTM)) })),
           .SDcols = grep("^PTM_.*_num$", names(temp), value=T, perl=T)]
      #create sorted combined PTM position column
      temp[PTM_num>0, PTM_pos_val_sorted := lapply(PTM_pos_val, sort)]
      #create sorted combined PTM type column
      temp[PTM_num>0, PTM_type_val_sorted := mapply(function(x, y){ x[order(y)] }, x=PTM_type_val, y=PTM_pos_val)]
      
      #recreate EG.PrecursorId column with valid types only
      temp[PTM_num==0, EG.PrecursorId.PTM.val := PTM_base_seq]
      temp[PTM_num>0, EG.PrecursorId.PTM.val := mapply(function(posit, types, baseseq){
        #step by step create res character string
        res <- character()
        pos.save <- 0
        #iterate for loop over each individual pos/type
        for(i in 1:(length(posit)) ){
          #build res from base sequence until pos
          res <- paste(res, substr(baseseq, start=pos.save, stop=posit[i]), sep="")
          #insert PTM type
          res <- paste(res, types[i], sep="")
          #update pos.save
          pos.save <- posit[i]+1
        }
        #finish by attaching last sequence substring
        return(paste(res, substr(baseseq, start=pos.save, stop=nchar(baseseq)), sep=""))
      }, posit=PTM_pos_val_sorted, types=PTM_type_val_sorted, baseseq=PTM_base_seq)]
      
      #create column combined with charge state for dcast
      temp[, PTM_group := paste("_", EG.PrecursorId.PTM.val, "_", gsub("_.*_(\\..)$", "\\1", EG.PrecursorId, perl=T), sep="")]
      
      #calculate mean peptide probability column -> target PTM probabilities only
      temp[, PTM_localization := NA_real_]
      #temp[PTM_num>0, PTM_localization := apply(.SD, 1, function(x){mean(as.numeric(unlist(x)), na.rm=T)}), .SDcols = grep("^PTM_.*_prob_val$", names(temp), value=T)]
      temp[PTM_0_num>0, PTM_localization := sapply(PTM_0_prob_val, function(x){mean(as.numeric(unlist(x)), na.rm=T)})]
      
      #delete all rows except PTM_0 and PTM_base_seq (needed for aa type extraction)
      for(col in setdiff(grep("^PTM_[1-9][0-9]*_.*$", names(temp), value=T, perl=T), c("PTM_base_seq", "PTM_localization", "PTM_group"))){
        temp[, paste(col) := NULL]
      }
      
      #return data table to top level
      return(temp)
    }
    
    #divide dataset into sub-datasets for parallelization
    data[, subd := sort(c(rep((1:(par.CPU-1)), round((.N/par.CPU), 0)), rep(par.CPU, .N-(round((.N/par.CPU), 0)*(par.CPU-1)) )))]
    data.sub <- split(data, by="subd")
    
    data <- copy(rbindlist(foreach(testp=data.sub, .packages=c("data.table", "stringr")) %dopar% {
      PTM.valid.extract.site.string(temp=testp, par.PTM=par.PTM)
    }))
    print("Done: EG.PTMLocalizationProbabilities-based EG.PrecursorId correction")
    #temp[PTM_0_num>0, .(PTM_0_pos_val, PTM_localization, EG.PrecursorId.PTM.val, EG.PrecursorId.PTM.val.charge)]
    
  } else if(any(par.level.col == 1, if(exists("par.level.col.nocut")){par.level.col.nocut == 1}, na.rm = T)){
    #if par.level = 1 & (par.level.col = 1 OR par.level.col.nocut = 1) -> nothing (use raw EG.PrecId later)
    
    #define function
    PTM.valid.extract.site.string <- function(temp, par.PTM){
      
      #extract number of PTM
      temp[, paste("PTM_0_num", sep="") := sapply(EG.PrecursorId, function(x){ length(unlist(strsplit(x, split=par.PTM[1], fixed=T)))-1 })]
      
      #define PTM_group
      temp[, PTM_group := EG.PrecursorId]
      
      #return data table to top level
      return(temp)
    }
    
    #divide dataset into sub-datasets for parallelization
    data[, subd := sort(c(rep((1:(par.CPU-1)), round((.N/par.CPU), 0)), rep(par.CPU, .N-(round((.N/par.CPU), 0)*(par.CPU-1)) )))]
    data.sub <- split(data, by="subd")
    
    data <- copy(rbindlist(foreach(testp=data.sub, .packages=c("data.table", "stringr")) %dopar% {
      PTM.valid.extract.site.string(temp=testp, par.PTM=par.PTM)
    }))
    print("Done: EG.PrecursorId-based PTM number extraction")
    
  } else if(if(exists("par.level.col.nocut")){par.level.col.nocut == 2}else{F}){
    #if MQ Modified.sequence without probability column
    
    #define function
    PTM.valid.extract.site.string <- function(temp, par.PTM){
      
      #extract number of PTM
      temp[, paste("PTM_0_num", sep="") := sapply(Modified.sequence, function(x){ length(unlist(strsplit(x, split=par.PTM[1], fixed=T)))-1 })]
      
      #define PTM_group
      temp[, PTM_group := Modified.sequence]
      
      #return data table to top level
      return(temp)
    }
    
    #divide dataset into sub-datasets for parallelization
    data[, subd := sort(c(rep((1:(par.CPU-1)), round((.N/par.CPU), 0)), rep(par.CPU, .N-(round((.N/par.CPU), 0)*(par.CPU-1)) )))]
    data.sub <- split(data, by="subd")
    
    data <- copy(rbindlist(foreach(testp=data.sub, .packages=c("data.table", "stringr")) %dopar% {
      PTM.valid.extract.site.string(temp=testp, par.PTM=par.PTM)
    }))
    print("Done: EG.PrecursorId-based PTM number extraction")
    
  } else if(par.level.col == 3L){
    #if MQ Modified.sequence AND probability extraction
    
    #define function
    PTM.valid.extract.site.string <- function(temp, par.PTM){
      
      #extract number of PTM
      temp[, paste("PTM_0_num", sep="") := sapply(Modified.sequence, function(x){ length(unlist(strsplit(x, split=par.PTM[1], fixed=T)))-1 })]
      
      #define PTM_group
      temp[, PTM_group := Modified.sequence]
      
      #create PTM valid row vector
      temp[, PTM_0_exist := F]
      temp[like(get(par.level.MQ.col), "\\("), PTM_0_exist :=T]
      
      #use pattern to extract respective PTM probabilities
      temp[PTM_0_exist==T, PTM_0_prob :=
             sapply(get(par.level.MQ.col), function(x){ as.numeric(gsub("\\((.*)\\)", "\\1", unlist(str_extract_all(x, "\\([^(]*\\)")), perl=T)) })]
      
      #extract positions of respective PTMs
      temp[PTM_0_exist==T, PTM_0_pos := sapply(get(par.level.MQ.col), function(x){
        #problem: if PTM on last/first position can create exemption if "_" removed -> use as anchor to mark start/stop of sequence and then substract
        temp <- nchar(unlist(strsplit(paste("_", x, "_", sep=""), split="\\([^(]*\\)")))
        temp <- temp-c(1, rep(0, length(temp)-2), 1)
        return(cumsum(temp[-length(temp)]))
      })]
      
      #extract PTM valid localization vector: eg T, T, F, T, F
      temp[PTM_0_exist==T, PTM_0_log := mapply(function(x, y){ rank(-x, ties.method = "first")<=y }, x=PTM_0_prob, y=PTM_0_num)]
      
      #extract PTM valid positions: eg 3,8,15
      temp[PTM_0_exist==T, PTM_0_pos_val := mapply(function(x, y){ list(x[y]) }, x=PTM_0_pos, y=PTM_0_log)]
      
      #extract PTM valid probabilities: eg 1,1,1
      temp[PTM_0_exist==T, PTM_0_prob_val := mapply(function(x, y){ list(x[y]) }, x=PTM_0_prob, y=PTM_0_log)]
      
      #calculate mean peptide probability column -> target PTM probabilities only
      temp[, PTM_localization := NA_real_]
      #temp[PTM_num>0, PTM_localization := apply(.SD, 1, function(x){mean(as.numeric(unlist(x)), na.rm=T)}), .SDcols = grep("^PTM_.*_prob_val$", names(temp), value=T)]
      temp[PTM_0_num>0, PTM_localization := sapply(PTM_0_prob_val, function(x){mean(as.numeric(unlist(x)), na.rm=T)})]
      
      #return data table to top level
      return(temp)
    }
    
    #divide dataset into sub-datasets for parallelization
    data[, subd := sort(c(rep((1:(par.CPU-1)), round((.N/par.CPU), 0)), rep(par.CPU, .N-(round((.N/par.CPU), 0)*(par.CPU-1)) )))]
    data.sub <- split(data, by="subd")
    
    data <- copy(rbindlist(foreach(testp=data.sub, .packages=c("data.table", "stringr")) %dopar% {
      PTM.valid.extract.site.string(temp=testp, par.PTM=par.PTM)
    }))
    print("Done: EG.PrecursorId-based PTM number extraction")
    
  }
  
} else
  if(par.level == 2){
  
  if(any(par.level.col == 0, if(exists("par.level.col.nocut")){par.level.col.nocut == 0}, na.rm = T)){
    #if par.level = 2 & (par.level.col = 0 OR par.level.col.nocut = 0) -> extract num, type and prob from EG.PTMLocProb
    
    #define function
    PTM.valid.extract.site.string <- function(temp, par.PTM){
      
      PTM.count <- 0
      for(PTM in par.PTM){
        
        #remove square brackets from PTM and mark round brackets as escaped, then insert into pattern
        pat <- gsub("[^:]*", gsub(")", "\\)", gsub("(", "\\(", gsub("(\\[|\\])", "", PTM, perl=T), fixed=T), fixed=T), "\\[[^:]*\\: (.{1,5})\\%\\]", fixed=T)
        pat.del <- gsub("[^:]*", gsub(")", "\\)", gsub("(", "\\(", gsub("(\\[|\\])", "", PTM, perl=T), fixed=T), fixed=T), "\\[[^:]*\\: .{1,5}\\%\\]", fixed=T)
        
        #create PTM valid row vector
        PTM.exist <- paste("PTM_", PTM.count, "_exist", sep="")
        temp[, paste(PTM.exist) :=F]
        temp[like(EG.PTMLocalizationProbabilities, pat), paste(PTM.exist) :=T]
        
        #use pattern to extract total number of respective PTM probabilities summed up
        temp[, paste("PTM_", PTM.count, "_num", sep="") :=
               sapply(EG.PTMLocalizationProbabilities, function(x){ round(sum(as.numeric(gsub(pat, "\\1", unlist(str_extract_all(x, pat)), perl=T)))/100, 0) })]
        
        #use pattern to extract respective PTM probabilities
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_prob", sep="") :=
               sapply(EG.PTMLocalizationProbabilities, function(x){ as.numeric(gsub(pat, "\\1", unlist(str_extract_all(x, pat)), perl=T))/100 })]
        
        #create EG.PTMLocalizationProbabilities with other PTMs deleted
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_seq", sep="") := EG.PTMLocalizationProbabilities]
        for(PTM.other in setdiff(par.PTM, PTM)){
          pat.other <- gsub("[^:]*", gsub(")", "\\)", gsub("(", "\\(", gsub("(\\[|\\])", "", PTM.other, perl=T), fixed=T), fixed=T), "\\[[^:]*\\: .{1,5}\\%\\]", fixed=T)
          temp[, paste(paste("PTM_", PTM.count, "_seq", sep="")) := gsub(pat.other, "", get(paste("PTM_", PTM.count, "_seq", sep="")), perl=T)]
        }
        
        #extract positions of respective PTMs
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_pos", sep="") := sapply(get(paste("PTM_", PTM.count, "_seq", sep="")), function(x){
          #problem: if PTM on last/first position can create exemption if "_" removed -> use as anchor to mark start/stop of sequence and then substract
          temp <- nchar(unlist(strsplit(x, split=pat.del)))
          temp <- temp-c(1, rep(0, length(temp)-2), 1)
          return(cumsum(temp[-length(temp)]))
        })]
        #extract PTM valid localization vector: eg T, T, F, T, F
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_log", sep="") :=
               mapply(function(x, y){ rank(-x, ties.method = "first")<=y },
                      x=get(paste("PTM_", PTM.count, "_prob", sep="")), y=get(paste("PTM_", PTM.count, "_num", sep="")))]
        
        #extract PTM valid positions: eg 3,8,15
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_pos_val", sep="") :=
               mapply(function(x, y){ list(x[y]) }, x=get(paste("PTM_", PTM.count, "_pos", sep="")), y=get(paste("PTM_", PTM.count, "_log", sep="")))]
        
        #extract PTM valid probabilities: eg 1,1,1
        temp[get(PTM.exist), paste("PTM_", PTM.count, "_prob_val", sep="") :=
               mapply(function(x, y){ list(x[y]) }, x=get(paste("PTM_", PTM.count, "_prob", sep="")), y=get(paste("PTM_", PTM.count, "_log", sep="")))]
        
        
        #increase PTM counter
        PTM.count <- PTM.count+1
      }
      
      #create base sequence column from EG.PTMLocProb
      temp[, PTM_base_seq := gsub("\\[[^:]*\\: .{1,5}\\%\\]", "", gsub("_", "", EG.PTMLocalizationProbabilities, fixed=T), perl=T)]
      
      #create summed PTM num column
      temp[, PTM_num := apply(.SD, 1, function(x){ sum(x, na.rm=T) }), .SDcols = grep("^PTM_.*_num$", names(temp), value=T, perl=T)]
      
      #create MQ modification specific peptides-like sequence modification column
      temp[, PTM_collapse_key := paste(PTM_base_seq, apply(.SD, 1, function(x){ paste(x[x>0], par.PTM[x>0], collapse=",", sep="") }), sep="_"),
           .SDcols = grep("^PTM_.*_num$", names(temp), value=T)]
      
      #calculate mean peptide probability column
      temp[, PTM_localization := NA_real_]
      #temp[PTM_num>0, PTM_localization := apply(.SD, 1, function(x){mean(as.numeric(unlist(x)), na.rm=T)}), .SDcols = grep("^PTM_.*_prob_val$", names(temp), value=T)]
      temp[PTM_0_num>0, PTM_localization := sapply(PTM_0_prob_val, function(x){mean(as.numeric(unlist(x)), na.rm=T)})]
      
      
      #create combined PTM position column
      temp[PTM_num>0, PTM_pos_val := apply(.SD, 1, function(x){ unlist(x) }), .SDcols = grep("^PTM_.*_pos_val$", names(temp), value=T, perl=T)]
      #create combined PTM type column
      temp[PTM_num>0, PTM_type_val := list(apply(.SD, 1, function(x){ unlist(mapply(function(y, z){rep(z, y)}, y=x, z=par.PTM)) })),
           .SDcols = grep("^PTM_.*_num$", names(temp), value=T, perl=T)]
      #create sorted combined PTM position column
      temp[PTM_num>0, PTM_pos_val_sorted := lapply(PTM_pos_val, sort)]
      #create sorted combined PTM type column
      temp[PTM_num>0, PTM_type_val_sorted := mapply(function(x, y){ x[order(y)] }, x=PTM_type_val, y=PTM_pos_val)]
      
      #recreate EG.PrecursorId column with valid types only
      temp[PTM_num==0, EG.PrecursorId.PTM.val := PTM_base_seq]
      temp[PTM_num>0, EG.PrecursorId.PTM.val := mapply(function(posit, types, baseseq){
        #step by step create res character string
        res <- character()
        pos.save <- 0
        #iterate for loop over each individual pos/type
        for(i in 1:(length(posit)) ){
          #build res from base sequence until pos
          res <- paste(res, substr(baseseq, start=pos.save, stop=posit[i]), sep="")
          #insert PTM type
          res <- paste(res, types[i], sep="")
          #update pos.save
          pos.save <- posit[i]+1
        }
        #finish by attaching last sequence substring
        return(paste(res, substr(baseseq, start=pos.save, stop=nchar(baseseq)), sep=""))
      }, posit=PTM_pos_val_sorted, types=PTM_type_val_sorted, baseseq=PTM_base_seq)]
      
      #create column combined with charge state for dcast
      temp[, PTM_group := paste("_", EG.PrecursorId.PTM.val, "_", gsub("_.*_(\\..)$", "\\1", EG.PrecursorId, perl=T), sep="")]
      
      
      #delete all columns except PTM_collapse_key, PTM_group, PTM_localization and PTM_0_prob_val
      for(col in setdiff(grep("^PTM_.*$", names(temp), value=T, perl=T), c("PTM_group", "PTM_collapse_key", "PTM_localization", "PTM_0_num"))){
        temp[, paste(col) := NULL]
      }
      
      #return data table to top level
      return(temp)
    }
    
    #divide dataset into sub-datasets for parallelization
    data[, subd := sort(c(rep((1:(par.CPU-1)), round((.N/par.CPU), 0)), rep(par.CPU, .N-(round((.N/par.CPU), 0)*(par.CPU-1)) )))]
    data.sub <- split(data, by="subd")
    
    data <- copy(rbindlist(foreach(testp=data.sub, .packages=c("data.table", "stringr")) %dopar% {
      PTM.valid.extract.site.string(temp=testp, par.PTM=par.PTM)
    }))
    print("Done: EG.PTMLocalizationProbabilities-based PTM position/probability extraction")
    #temp[PTM_0_num>0, .(PTM_group, PTM_collapse_key, PTM_localization)]
    
  } else if(par.level == 2 & any(par.level.col == 1, if(exists("par.level.col.nocut")){par.level.col.nocut == 1}, na.rm = T)){
    #if par.level = 2 & (par.level.col = 1 OR par.level.col.nocut = 1) -> extract num and type from EG.PrecId
    
    #define function
    PTM.valid.extract.site.string <- function(temp, par.PTM){
      
      PTM.count <- 0
      for(PTM in par.PTM){
        
        #extract number of PTM
        temp[, paste("PTM_", PTM.count, "_num", sep="") := sapply(EG.PrecursorId, function(x){ length(unlist(strsplit(x, split=PTM, fixed=T)))-1 })]
        
        #increase PTM counter
        PTM.count <- PTM.count+1
      }
      
      #create PTM base seq from EG.PrecId
      temp[, PTM_base_seq := gsub("^(.*)\\..$", "\\1", gsub("\\[[^[]*\\]", "", gsub("_", "", EG.PrecursorId, fixed=T), perl=T), perl=T)]
      
      #create MQ modification specific peptides-like sequence modification column
      temp[, PTM_collapse_key := paste(PTM_base_seq, apply(.SD, 1, function(x){ paste(x[x>0], par.PTM[x>0], collapse=",", sep="") }), sep="_"),
           .SDcols = grep("^PTM_.*_num$", names(temp), value=T)]
      
      #delete all rows except PTM_collapse_key
      for(col in setdiff(grep("^PTM_.*$", names(temp), value=T, perl=T), c("PTM_collapse_key", "PTM_0_num"))){
        temp[, paste(col) := NULL]
      }
      
      #define PTM_group
      temp[, PTM_group := EG.PrecursorId]
      
      #return data table to top level
      return(temp)
    }
    
    #divide dataset into sub-datasets for parallelization
    data[, subd := sort(c(rep((1:(par.CPU-1)), round((.N/par.CPU), 0)), rep(par.CPU, .N-(round((.N/par.CPU), 0)*(par.CPU-1)) )))]
    data.sub <- split(data, by="subd")
    
    data <- copy(rbindlist(foreach(testp=data.sub, .packages=c("data.table", "stringr")) %dopar% {
      PTM.valid.extract.site.string(temp=testp, par.PTM=par.PTM)
    }))
    print("Done: EG.PrecursorId-based PTM position extraction")
    #temp[PTM_0_num>0, .(PTM_group, PTM_collapse_key)]
    
  }
  }

#finish parallel processing
stopCluster(cl)
data[, subd := NULL]




#site-level only: row expansion
if(par.level == 0){
  
  #expand dataset row-wise
  data[, rown := .I]
  
  if(any(par.level.col == 0, if(exists("par.level.col.nocut")){par.level.col.nocut == 0}, na.rm = T)){
    #localization column extracted from EG.PTMLocProb
    data.exp <- copy(rbind(data[PTM_0_num>0, .(PTM_localization = unlist(PTM_localization), PTM_0_pos_val = unlist(PTM_0_pos_val)), by=rown],
                           data[PTM_0_num==0, .(rown, PTM_localization = NA, PTM_0_pos_val = NA_integer_)])[order(rown)])
    #merge missing data onto dataset, excluding PTM_localization column
    data <- merge(x=data.exp, y=data[, .SD, .SDcols = setdiff(names(data), c("PTM_localization", "PTM_0_pos_val"))],
                  by.x="rown", by.y="rown", all.x=T, all.y=F)
    
  } else {
    #localization column is EG.PTMAssayProbability
    data.exp <- copy(rbind(data[PTM_0_num>0, .(PTM_0_pos_val = unlist(PTM_0_pos_val)), by=rown],
                           data[PTM_0_num==0, .(rown, PTM_0_pos_val = NA_integer_)])[order(rown)])
    #merge missing data onto dataset, excluding pSTY.loc.val column
    data <- merge(x=data.exp, y=data[, .SD, .SDcols = setdiff(names(data), c("PTM_0_pos_val"))],
                  by.x="rown", by.y="rown", all.x=T, all.y=F)
  }
  
  data[, rown := NULL]
  
  #create PTM_0_aa column
  data[, PTM_0_aa := substr(PTM_base_seq, start=PTM_0_pos_val, stop=PTM_0_pos_val)]
  
  print("Done: site-level expansion")
}


#PTM-localization formatting
if(par.level.col %in% c(0, 3)){
  data[, PTM_localization := as.numeric(PTM_localization)]
}
if(par.level.col == 1){
  data[, PTM_localization := as.numeric(EG.PTMAssayProbability)]
}



#transform data from long into wide format based on condition column
if(par.cond==0L){
  
  #create cond.names column -> if more than 1 main column, equals combination of main.name and cond.col arguments
  if(length(main.name)>1){
    cond.names <- as.vector(outer(main.name, data[, sort(unique(get(par.cond.col)))], paste, sep="_"))
  } else {
    #if only 1 column, equals cond.col arguments
    cond.names <- data[, sort(unique(get(par.cond.col)))]
  }
  
  #site-level and probability included
  if(par.level == 0 & par.level.col %in% c(0,1,3)){
    
    #create vector of all user defined columns to keep during transformation -> PTM_group and co excluded; PTM_0_num and PTM_0_aa included for multiplicity tag
    #keep <- c(setdiff(names(data), c(grep("^PTM_.*$", names(data), value=T, perl=T), par.cond.col, main.name) ), "PTM_0_num", "PTM_0_aa")
    keep <- c(setdiff(else.name, par.cond.col), "PTM_0_num", "PTM_0_aa")
    
    #perform transformation via dcast
    data <- cbind(dcast(data, PTM_group + PTM_0_pos_val ~ get(par.cond.col), value.var = main.name,
                        fun.aggregate = function(x) if(sum(is.na(x)) == length(x)){NA_real_}else{sum(x, na.rm = T)}), #sum only if >0 non-NA
                  dcast(data, PTM_group + PTM_0_pos_val ~ ., value.var = keep,
                        fun.aggregate = function(x) paste(unique(unlist(strsplit(as.character(x), split=";", fixed=T))), collapse=";"))[, -c(1,2)],
                  dcast(data, PTM_group + PTM_0_pos_val ~ ., value.var = c("PTM_localization"),
                        fun.aggregate = function(x) if(sum(is.na(x)) == length(x)){NA_real_}else{max(x, na.rm=T)})[, -c(1,2)]) #max only if >0 non-NA
    setnames(data, ".", "PTM_localization")
    
    print("Done: Long- to wide-format transformation")
  }
  #site-level and probability ignored
  if(par.level == 0 & par.level.col == 2){
    
    #create vector of all user defined columns to keep during transformation -> PTM_group and co excluded; PTM_0_num and PTM_0_aa included for multiplicity tag
    #keep <- c(setdiff(names(data), c(grep("^PTM_.*$", names(data), value=T, perl=T), par.cond.col, main.name) ), "PTM_0_num", "PTM_0_aa")
    keep <- c(setdiff(else.name, par.cond.col), "PTM_0_num", "PTM_0_aa")
    
    #perform transformation via dcast
    data <- cbind(dcast(data, PTM_group + PTM_0_pos_val ~ get(par.cond.col), value.var = main.name,
                        fun.aggregate = function(x) if(sum(is.na(x)) == length(x)){NA_real_}else{sum(x, na.rm = T)}), #sum only if >0 non-NA
                  dcast(data, PTM_group + PTM_0_pos_val ~ ., value.var = keep,
                        fun.aggregate = function(x) paste(unique(unlist(strsplit(as.character(x), split=";", fixed=T))), collapse=";"))[, -c(1,2)])
    
    print("Done: Long- to wide-format transformation")
  }
  
  #localized/ModSpec peptide-level and probability included
  if(par.level %in% c(1,2) & par.level.col %in% c(0,1,3)){
    
    #for ModSpec only, include PTM_collapse_key in keep -> should always be unique
    if(par.level==1){
      #create vector of all user defined columns to keep during transformation -> PTM_group and co excluded; PTM_0_num included
      #keep <- c(setdiff(names(data), c(grep("^PTM_.*$", names(data), value=T, perl=T), par.cond.col, main.name) ))
      keep <- c(setdiff(else.name, par.cond.col), "PTM_0_num")
    }
    if(par.level==2){
      #create vector of all user defined columns to keep during transformation -> PTM_group and co excluded; PTM_0_num included
      #keep <- c(setdiff(names(data), c(grep("^PTM_.*$", names(data), value=T, perl=T), par.cond.col, main.name) ), "PTM_collapse_key")
      keep <- c(setdiff(else.name, par.cond.col), "PTM_0_num", "PTM_collapse_key")
    }
    
    #perform transformation via dcast
    data <- cbind(dcast(data, PTM_group ~ get(par.cond.col), value.var = main.name,
                        fun.aggregate = function(x) if(sum(is.na(x)) == length(x)){NA_real_}else{sum(x, na.rm = T)}), #sum only if >0 non-NA
                  dcast(data, PTM_group ~ ., value.var = keep,
                        fun.aggregate = function(x) paste(unique(unlist(strsplit(as.character(x), split=";", fixed=T))), collapse=";"))[, -c(1)],
                  dcast(data, PTM_group ~ ., value.var = c("PTM_localization"),
                        fun.aggregate = function(x) if(sum(is.na(x)) == length(x)){NA_real_}else{max(x, na.rm=T)})[, -c(1)]) #max only if >0 non-NA
    setnames(data, ".", "PTM_localization")
    
    print("Done: Long- to wide-format transformation")
  }
  #localized/ModSpec peptide-level and probability ignored
  if(par.level %in% c(1,2) & par.level.col == 2){
    
    #for ModSpec only, include PTM_collapse_key in keep -> should always be unique
    if(par.level==1){
      #create vector of all user defined columns to keep during transformation -> PTM_group and co excluded; PTM_0_num and PTM_0_aa included for multiplicity tag
      #keep <- c(setdiff(names(data), c(grep("^PTM_.*$", names(data), value=T, perl=T), par.cond.col, main.name) ))
      keep <- c(setdiff(else.name, par.cond.col), "PTM_0_num")
    }
    if(par.level==2){
      #create vector of all user defined columns to keep during transformation -> PTM_group and co excluded; PTM_0_num and PTM_0_aa included for multiplicity tag
      #keep <- c(setdiff(names(data), c(grep("^PTM_.*$", names(data), value=T, perl=T), par.cond.col, main.name) ), "PTM_collapse_key")
      keep <- c(setdiff(else.name, par.cond.col), "PTM_0_num", "PTM_collapse_key")
    }
    
    #perform transformation via dcast
    data <- cbind(dcast(data, PTM_group ~ get(par.cond.col), value.var = main.name,
                        fun.aggregate = function(x) if(sum(is.na(x)) == length(x)){NA_real_}else{sum(x, na.rm = T)}), #sum only if >0 non-NA
                  dcast(data, PTM_group ~ ., value.var = keep,
                        fun.aggregate = function(x) paste(unique(unlist(strsplit(as.character(x), split=";", fixed=T))), collapse=";"))[, -c(1)])
    
    print("Done: Long- to wide-format transformation")
  }
} else {
  #create cond.names column which equals main.name only
  cond.names <- main.name
}




#prepare consolidation
#delete all non-target PTM rows? -> leave that to manual Perseus pre-filtering; already done for site-level only
#site-level preparation -> key on PTM_mult123_tag
if(par.level==0){
  #create gene/protein column based on user defined gene/protein column -> ignore first entry if ""
  data[, PTM_genprot := sapply(get(par.level.col.genprot), function(x){
    strsplit(x, split="(;|,)", perl=T)[[1]][!strsplit(x, split="(;|,)", perl=T)[[1]] == ""][1]
    })]
  
  #remove rows that are empty in genprot
  if(data[(is.na(PTM_genprot) | PTM_genprot==""| PTM_genprot=='""'), .N]>0){
    print(paste("Cave: Removed ", data[(is.na(PTM_genprot) | PTM_genprot==""| PTM_genprot=='""'), .N], " rows due to missing gene/prot information", sep=""))
  }
  data <- data[!(is.na(PTM_genprot) | PTM_genprot==""| PTM_genprot=='""'), ]
  
  #create position column based on position column -> ignore first entry if ""
  data[, PTM_pep_pos := sapply(PEP.PeptidePosition, function(x){
    strsplit(x, split="(;|,)", perl=T)[[1]][!strsplit(x, split="(;|,)", perl=T)[[1]] == ""][1]
    })]
  
  #create multiplicity vector, reading out number of target PTM per peptides
  data[, PTM_mult123 := as.integer(PTM_0_num)]
  data[PTM_mult123>=3, PTM_mult123 := 3]
  #PTM_mult123_tag, but renamed to PTM_collapse_key
  data[, PTM_collapse_key := paste(PTM_genprot, "_", PTM_0_aa, (PTM_0_pos_val+as.integer(PTM_pep_pos)-1), "_M", PTM_mult123, sep="")]
  setkey(data, "PTM_collapse_key")
}
#localized peptide-level preparation -> key on localized EG.PrecId without charge
if(par.level==1){
  #create PTM_group without charge state
  data[, PTM_collapse_key := gsub("^(.*)\\..$", "\\1", PTM_group, perl=T)]
  
  #delete all PTM info except target PTM from collapse key
  for(PTM.other in par.PTM[-1]){
    data[, PTM_collapse_key := gsub(PTM.other, "", PTM_collapse_key, fixed=T)]
  }
  
  #key on PTM_collapse_key
  setkey(data, "PTM_collapse_key")
}
#ModSpec peptide-level preparation -> key on PTM_ModSpec
if(par.level==2){
  #key on PTM_collapse_key
  setkey(data, "PTM_collapse_key")
}




#initiate consolidation
#prepare list of subset data dt, which are passed to the cores individually
if(par.cond == 0L){
  if(par.level.col %in% c(0,1,3)){
    keep <- c(cond.names, setdiff(else.name, par.cond.col), "PTM_0_num", "PTM_group", "PTM_collapse_key", "PTM_localization")
  }
  if(par.level.col == 2){
    keep <- c(cond.names, setdiff(else.name, par.cond.col), "PTM_0_num", "PTM_group", "PTM_collapse_key")
  }
} else if(par.cond == 1L){
  if(par.level.col %in% c(0,1,3)){
    keep <- c(cond.names, else.name, "PTM_0_num", "PTM_group", "PTM_collapse_key", "PTM_localization")
  }
  if(par.level.col == 2){
    keep <- c(cond.names, else.name, "PTM_0_num", "PTM_group", "PTM_collapse_key")
  }
}


#create vector with PTM_collapse_key group intervals
sub_vec <- data[, round(quantile(1:length(unique(PTM_collapse_key)), probs = seq(0,1,length.out = par.CPU+1)[-(par.CPU+1)]),0)]

#create PTM_subgrp column based on sub_vec, which marks groups of rows to be subset into lists for parallel processing
data[, PTM_subgrp := data[, .(PTM_grp_nr = rep(.GRP, .N)), by=PTM_collapse_key #create PTM_collapse_key group numbers
                          ][, .(PTM_grp_test = c(.BY %in% sub_vec, rep(F, .N-1))), by=PTM_grp_nr #mark each first row of PTM_grp_nr that matches vec as T
                            ][, cumsum(PTM_grp_test)] ] #create PTM_subgrp column

#create subdt: 8 (=number of cores) sub-data tables within list to be passed into function separately
subdt <- split(data[, .SD, .SDcols = c(keep, "PTM_subgrp")], by="PTM_subgrp")

# #check if total number of rows in subdt fits those in test
# sum(unlist(lapply(subdt, nrow)))
# #check if first and last entry different
# for(i in 1:(par.CPU-1)){
#   print(paste(subdt[[i]][.N, .(PTM_collapse_key)]==subdt[[i+1]][1, .(PTM_collapse_key)],
#               subdt[[i]][.N, .(PTM_collapse_key)], subdt[[i+1]][1, .(PTM_collapse_key)], sep=" ; "))
# }

#load function: consolidation based on equal division (2018-03-21)
consolidate <- function(cons, cond){
  #check if only NA -> return data table with 1 row filled only with NA immediately
  if(cons[, all(is.na(.SD))]){
    return(as.data.table(matrix(rep(NA, length(cond)), nrow = 1, dimnames = list(NULL, cond))))
  }
  
  #sort rows from lowest median intensity to highest, to later make sure that (if needed), lowest median rows are kicked out first
  cons <- cons[order(cons[, apply(.SD, 1, function(x) median(x, na.rm=T)), .SDcols = cond])]
  
  #check if consolidation neccessary at all -> check if all conditions have no missing values
  while(cons[, any(!sapply(.SD, function(x) length(which(!is.na(x))))[cons[, sapply(.SD, function(x) length(which(!is.na(x))))]>0]==.N)]){
    #create loop variables for first iteration
    tempc <- NULL
    
    #start for loop iteration
    for(num in (1:cons[, .N])){
      #check if any NA in this row; if yes, execute normalization, if not simply copy values
      if(any(is.na(unlist(transpose(cons)[, num, with=F])))){
        #create logical vector reading out NA entries
        logv <- is.na(unlist(transpose(cons)[, num, with=F]))
        
        #this line calculates the missing value substitution values, and then reorders them together with the non-missing values
        #finally, it adds them into the tempc matrix via rbind
        #need to separate: if >1 NA value, can use apply; but if less then directly use median instead of apply
        if(length(which(logv))==1){
          #trigger this function if there is exactly 1 NA in this row
          tempc <- rbind(tempc, c(median(sapply(transpose(cons)[, -num, with=F], function(x) median(unlist(transpose(cons)[, num, with=F])/x, na.rm=T)*x[logv]), na.rm=T),
                                  unlist(transpose(cons)[, num, with=F], use.names = F)[!logv])[order(c(which(logv), which(!logv)))])
        } else {
          #trigger this function if there are more then 1 NA in this row
          tempc <- rbind(tempc, c(apply(sapply(transpose(cons)[, -num, with=F], function(x) median(unlist(transpose(cons)[, num, with=F])/x, na.rm=T)*x[logv]), 1,
                                        function(x) median(x, na.rm=T)),
                                  unlist(transpose(cons)[, num, with=F], use.names = F)[!logv])[order(c(which(logv), which(!logv)))])
        }
        
      } else {
        #copy values instead
        tempc <- rbind(tempc, unlist(transpose(cons)[, num, with=F]))
      }
    }
    
    #when for loop done, check if tempc matrix equals starting matrix
    if(identical(cons, setNames(as.data.table(tempc), cond))){
      #if identical, means that consolidation does not work anymore -> remove cons row with lowest number of
      cons <- cons[-(grep(min(apply(cons, 1, function(x) length(which(!is.na(x))))), apply(cons, 1, function(x) length(which(!is.na(x)))))[1]), ]
      
    } else {
      #if consolidation was triggered, then make tempc new cons
      cons <- setNames(as.data.table(tempc), cond)
    }
    #now, check while condition again and if neccessary, rerurn while loop
  }
  
  #whether or not while loop was triggered, should now have dt without NA except for completely empty conditions
  #thus, calculate col medians and return (CAVE: we are in intensity space; CAVE: transposition)
  #option 1: center all rows around 0 in log space, then med, then exponential transformation; preferrable because no row is over-/under-weighted?
  #return(setNames(cons[, lapply(transpose(lapply(transpose(.SD), function(x) log(x)-median(log(x), na.rm=T))), function(x) exp(median(x)))], cond))
  #option 2: sum up all intensities in log space -> high intensity rows gain more weight?!?
  #return(setNames(cons[, lapply(.SD, function(x) exp(sum(log(x))))], cond)) #problem in stoich -> exponents too high!
  #return(setNames(cons[, lapply(.SD, function(x) exp(sum(log(x))/length(x)))], cond))
  #option 2b: sum all intensities in normal intensity space -> put higher weight on high intensities
  return(setNames(cons[, lapply(.SD, function(x) sum(x))], cond))
  #option 3: mean or median of all intensities in log space -> high intensity rows gain more weight?!?
  #return(setNames(cons[, lapply(.SD, function(x) exp(mean(log(x))))], cond))
  #return(setNames(cons[, lapply(.SD, function(x) exp(median(log(x))))], cond))
}

if(par.level==0){ print("Done: Preparation for site-level collapse") }
if(par.level==1){ print("Done: Preparation for localized peptide-level collapse") }
if(par.level==2){ print("Done: Preparation for modification specific peptide-level collapse") }


#initiate parallel processing
cl <- makeCluster(par.CPU)
registerDoParallel(cl)

#start collapse: include localization
if(par.level.col %in% c(0,1,3)){
  #load function
  consdt.pept.num <- function(temp, cond, type){
    #remove all entries without or with zero-length keycol entries
    temp <- temp[!is.na(PTM_collapse_key) & nchar(PTM_collapse_key)>0,]
    
    #need at least two valid entries per row to consider further
    temp <- temp[temp[, apply(.SD, 1, function(x) length(which(!is.na(x)))) >1, .SDcols = cond],]
    
    #write number of entries per keycol group
    temp[, PTM_collapse_key_num := .N, by = PTM_collapse_key]
    
    #perform consolidation in chunks devided by PTM_collapse_key, rbind with PTM_collapse_key_entries equal to one, and overwrite temp
    if(type==0){
      #perform linear modelling
      temp <- copy(cbind(temp[, consolidate(cons = .SD, cond = cond), .SDcols = cond, by=PTM_collapse_key],
                         temp[, .(PTM_localization = max(PTM_localization, na.rm=T)), by=PTM_collapse_key][,-1],
                         temp[, lapply(.SD, function(x) paste(unique(x), collapse=";")),
                              .SDcols = setdiff(names(temp), c("PTM_collapse_key", "PTM_localization", cond)), by=PTM_collapse_key][,-1]))
    }
    if(type==1){
      #simply add intensities
      temp <- copy(cbind(temp[, lapply(.SD, function(x){ sum(x, na.rm = T) }), .SDcols = cond, by=PTM_collapse_key],
                         temp[, .(PTM_localization = max(PTM_localization, na.rm=T)), by=PTM_collapse_key][,-1],
                         temp[, lapply(.SD, function(x) paste(unique(x), collapse=";")),
                              .SDcols = setdiff(names(temp), c("PTM_collapse_key", "PTM_localization", cond)), by=PTM_collapse_key][,-1]))
      
      #re-delete all 0s in columns that were generated by summing function
      for(col in cond.names){
        temp[get(col)==0, paste(col) := NA]
      }
    }
    
    #export temp into global environment
    return(temp)
  }
  
  #parallel collapse
  data <- copy(rbindlist(foreach(testp=subdt,.packages=c("data.table")) %dopar% {
    consdt.pept.num(temp=testp, cond=cond.names, type=par.agg)
  }))
  
  #overwrite all -Inf PTM_localization that occured from rows without non-missing localizations
  data[is.infinite(PTM_localization), PTM_localization := NA]
  
  print("Done: Collapse including maximum localization probability")
}

#start collapse: ignore localization
if(par.level.col == 2){
  #load function
  consdt.pept.num <- function(temp, cond, type){
    #remove all entries without or with zero-length keycol entries
    temp <- temp[!is.na(PTM_collapse_key) & nchar(PTM_collapse_key)>0,]
    
    #need at least two valid entries per row to consider further
    temp <- temp[temp[, apply(.SD, 1, function(x) length(which(!is.na(x)))) >1, .SDcols = cond],]
    
    #write number of entries per keycol group
    temp[, PTM_collapse_key_num := .N, by = PTM_collapse_key]
    
    #perform consolidation in chunks devided by PTM_collapse_key, rbind with PTM_collapse_key_entries equal to one, and overwrite temp
    if(type==0){
      #perform linear modelling
      temp <- copy(cbind(temp[, consolidate(cons = .SD, cond = cond), .SDcols = cond, by=PTM_collapse_key],
                         temp[, lapply(.SD, function(x) paste(unique(x), collapse=";")),
                              .SDcols = setdiff(names(temp), c("PTM_collapse_key", cond)), by=PTM_collapse_key][,-1]))
    }
    if(type==1){
      #simply add intensities
      temp <- copy(cbind(temp[, lapply(.SD, function(x){ sum(x, na.rm = T) }), .SDcols = cond, by=PTM_collapse_key],
                         temp[, lapply(.SD, function(x) paste(unique(x), collapse=";")),
                              .SDcols = setdiff(names(temp), c("PTM_collapse_key", cond)), by=PTM_collapse_key][,-1]))
      
      #re-delete all 0s in columns that were generated by summing function
      for(col in cond.names){
        temp[get(col)==0, paste(col) := NA]
      }
    }
    
    #export temp into global environment
    return(temp)
  }
  
  #parallel collapse
  data <- copy(rbindlist(foreach(testp=subdt,.packages=c("data.table")) %dopar% {
    consdt.pept.num(temp=testp, cond=cond.names, type=par.agg)
  }))
  
  print("Done: Collapse ignoring localization probability")
}

#close parallel processing
stopCluster(cl)






#PTM localization probability filtering
if(par.level.col %in% c(0,1,3)){
  if(par.level.cutoff>0){
    print(paste("Cave: Filtered ", data[!(PTM_localization>=par.level.cutoff | is.na(PTM_localization)), .N], " rows below cutoff of ", par.level.cutoff, sep=""))
    
    #filter all rows below cutoff
    data <- data[PTM_localization>=par.level.cutoff | is.na(PTM_localization), ]
  }
}


#amino acid type extrection if site-level
if( par.level == 0){
  data[, PTM_0_aa := gsub("^.*_(.).*_M[1-3]$", "\\1", PTM_collapse_key, perl=T)]
}


#PTM motif sequence matching only if site-level & FASTA file path defined
if( par.level == 0 & if(exists("par.FASTA.file")){!identical(par.FASTA.file, character(0))}else{F} ){
  #par.FASTA.parse <- ".*\\|(.*)\\|.*"
  #par.FASTA.parse <- ".*GN=([^ ]*) .*"
  
  #read fasta file into R
  fasta <- readAAStringSet(par.FASTA.file)
  fasta <- data.table(header = names(fasta), ident_seq = paste(fasta))
  fasta[, ident := gsub(par.FASTA.parse, "\\1", header, perl=T)]
  fasta[, pat_match := grepl(par.FASTA.parse, header)]
  #keep only rows with parse pattern
  fasta <- fasta[pat_match==T, ]
  fasta[, header := NULL]
  fasta[, pat_match := NULL]
  
  #need to first check if par.level.col.genprot are part of FASTA file
  data[, ident.first := gsub("^;(.*)$", "\\1", get(par.level.col.genprot), perl=T)]
  data[, ident.first := gsub("^([^;]+);.*", "\\1", ident.first, perl=T)]
  data[, ident.fasta.match := get(par.level.col.genprot) %in% fasta$ident]
  print(paste("Cave: could not map ", data[ident.fasta.match==F, .N], " rows out of ", data[, .N], " for FASTA info", sep=""))
  
  #map fasta information onto data
  data <- merge(x=data, y=fasta, by.x="ident.first", by.y="ident", all.x=T, all.y=F)
  
  #re-extract PTM position from PTM_collapse_key
  data[, PTM_pos := as.integer(gsub("^.*_.(.*)_M.$", "\\1", PTM_collapse_key, perl=T))]
  
  #create PTM sequence window, starting by determining sub-cases
  data[, seq_end_length := nchar(ident_seq)-PTM_pos]
  
  #if PTM_pos >=16 AND seq_end_length >=15
  data[PTM_pos>=16 & seq_end_length>=15, PTM_seq := substr(ident_seq, (as.integer(PTM_pos)-15), (as.integer(PTM_pos)+15))]
  
  #if PTM_pos <16 AND seq_end_length >=15
  data[PTM_pos<16 & seq_end_length>=15, PTM_seq := 
               mapply(function(x,y){paste(paste(rep("_", (16-x)), collapse=""), substr(y, 1, (as.integer(x)+15)), sep="")}, x=PTM_pos, y=ident_seq)]
  
  #if PTM_pos >=16 AND seq_end_length <15
  data[PTM_pos>=16 & (seq_end_length<15 & seq_end_length>=0), PTM_seq :=
               mapply(function(x,y){paste(substr(y, (as.integer(x)-15), (nchar(y))), paste(rep("_", (15-(nchar(y)-x))), collapse=""), sep="")}, x=PTM_pos, y=ident_seq)]
  
  #if PTM_pos <16 AND seq_end_length <15
  data[PTM_pos<16 & (seq_end_length<15 & seq_end_length>=0), PTM_seq :=
               mapply(function(x,y){paste(paste(rep("_", (16-x)), collapse=""),
                                          substr(y, (as.integer(x)-15), (nchar(y))), paste(rep("_", (15-(nchar(y)-x))), collapse=""), sep="")}, x=PTM_pos, y=ident_seq)]
  
  data[, seq_end_length := NULL]
  data[, ident.fasta.match := NULL]
  data[, ident_seq := NULL]
  data[, PTM_pos := NULL]
  
  print("Done: PTM sequence motif extraction")
}



#initiate stoichiometry calculation
if(par.level==1 & if(exists("par.level.stoich")){par.level.stoich == 0}else{F}){
  
  #create PTM_stoich_key sequence, deleting target PTM from PTM_collapse_key
  data[, PTM_stoich_key := gsub(par.PTM[1], "", PTM_collapse_key, fixed=T)]
  setkey(data, PTM_stoich_key, PTM_0_num, PTM_collapse_key)
  
  #prepare list of subset data dt, which are passed to the cores individually
  if(par.cond == 0L){
    if(par.level.col %in% c(0,1,3)){
      keep <- c(cond.names, setdiff(else.name, par.cond.col), "PTM_0_num", "PTM_group", "PTM_collapse_key", "PTM_collapse_key_num", "PTM_stoich_key",
                "PTM_localization")
    }
    if(par.level.col == 2){
      keep <- c(cond.names, setdiff(else.name, par.cond.col), "PTM_0_num", "PTM_group", "PTM_collapse_key", "PTM_collapse_key_num", "PTM_stoich_key")
    }
  } else if(par.cond == 1L){
    if(par.level.col %in% c(0,1,3)){
      keep <- c(cond.names, else.name, "PTM_0_num", "PTM_group", "PTM_collapse_key", "PTM_collapse_key_num", "PTM_stoich_key",
                "PTM_localization")
    }
    if(par.level.col == 2){
      keep <- c(cond.names, else.name, "PTM_0_num", "PTM_group", "PTM_collapse_key", "PTM_collapse_key_num", "PTM_stoich_key")
    }
  }
  
  
  
  #create vector describing split points for dataset between collapsing units
  RAWseq_vec <- data[, unique(PTM_stoich_key)]
  vec <- 1L
  for(col in RAWseq_vec[round(quantile(1:length(RAWseq_vec), probs = seq(0,1,length.out = par.CPU+1)[-1]),0)]){
    vec <- c(vec, grep(col, data$PTM_stoich_key, fixed = T)[length(grep(col, data$PTM_stoich_key, fixed=T))])
  }
  
  #create subdt: 8 (=number of cores) sub-data tables within list to be passed into function separately
  subdt <- list(data[1:vec[2], .SD, .SDcols=keep])
  for(col in 2:par.CPU){
    subdt <- c(subdt, list(data[(vec[col]+1):vec[col+1], .SD, .SDcols=keep]))
  }
  
  print("Done: Stoichiometry preparation")
  
  # #check if total number of rows in subdt fits those in test
  # sum(unlist(lapply(subdt, nrow)))
  # #check if first and last entry different
  # for(i in 1:(par.CPU-1)){
  #   print(paste(subdt[[i]][.N, .(PTM_stoich_key)]==subdt[[i+1]][1, .(PTM_stoich_key)],
  #               subdt[[i]][.N, .(PTM_stoich_key)], subdt[[i+1]][1, .(PTM_stoich_key)], sep=" ; "))
  # }
  
  #load functions
  rlmfun <- function(x, skiperror=T, multiphos=2, PTMnum = PTMnum){
    #to avoid singularity, need to have overfit -> at least as many valid value matches over all peptides as conditions
    #eg, if we have three occseq peptide variants, we need at least 3 conditions with complete matches over all occseq peptide variants!
    #initiate logical selv vector which will define number of rows that can be used for rlm modeling via TRUE, and excludes duplicate peptide rows
    selv <- !duplicated(x)
    #tried to solve problem with duplicate rows leading to singularity, but also kicks out working rows -> disregarded for now, using error trap instead
    # if(length(which(apply(x, 2, function(y) !any(duplicated(y)) )))>=nrow(x)){
    #   selv <- !duplicated(x)
    # } else {
    #   selv <- c(T, rep(F, (nrow(x)-1)))
    # }
    
    #multiphos: if 1, should only consider occgrp with PTMnum equaling 2 different unique PTMnum; if 0, should only equal 2 different total PTMnum; if 2 use all
    if(multiphos==1){
      if(length(unique(PTMnum))!=2){
        x[!(PTMnum == sort(unique(PTMnum))[1]| PTMnum == sort(unique(PTMnum))[2])] <- NA
      }
    }
    if(multiphos==0){
      if(length(unique(PTMnum))!=2 | length(PTMnum)!=2){
        #calculate lowest two PTMnum values and for each individually pick the matrix row with the highest medium intensity -> overwrite all other entries with NA
        x[!(rowMedians(x, na.rm=T)==sort(rowMedians(x[PTMnum == sort(unique(PTMnum))[1],,drop=F], na.rm=T), decreasing = T)[1] |
              rowMedians(x, na.rm=T)==sort(rowMedians(x[PTMnum == sort(unique(PTMnum))[2],,drop=F], na.rm=T), decreasing = T)[1])] <- NA
      }
    }
    
    # #check if there are duplicate values between rows; if yes, add rnorm onto rows to enable calculation
    # if(any(apply(as.matrix(bla[, .SD, .SDcols=condi]), 2, function(x) any(duplicated(x))))){}
    # duplicated(as.matrix(bla[, .SD, .SDcols=condi]))
    # t(apply(as.matrix(bla[, .SD, .SDcols=condi]), 1, function(x) x+rnorm(n=length(x), mean=(median(x, na.rm=T)/100))))
    
    #count number of conditions with full matches and check if at least equal to number of occseq
    if(!(length(which(apply(x, 2, function(x) length(which(!is.na(x))))==nrow(x))) >= nrow(x))){
      #need to remove x-row with lowest number of non-NA-values and check again
      #if all values removed, then no calculation possible
      #at same time remember their original position, so that reconstruction of slope for different peptide rows is possible
      nanv <- apply(x, 1, function(x) length(which(!is.na(x))))
      
      #create permanent selection vector for formula creation subsetting
      selv <- nanv>min(nanv[selv], na.rm=T) & selv
      
      #check if at least 2 positive lines remaining in selv, because otherwise rlm will not be possible
      if(length(which(selv))>1){
        y <- x[(1:nrow(x))[selv],]
        
        #repeat reducing selv vector if matching conditions still not met
        while(!(length(which(apply(y, 2, function(x) length(which(!is.na(x))))==nrow(y))) >= nrow(y))){
          #create second selection vector, which is then compared to initial x-based selection vector and overwrite selv
          nanv2 <- apply(y, 1, function(x) length(which(!is.na(x))))
          selv <- nanv>min(nanv2, na.rm=T) & selv
          y <- x[(1:nrow(x))[selv],]
          
          #check again if at least 2 positive lines remaining in selv, because otherwise rlm will not be possible
          if(length(which(selv))<2){break}
        }
      }
    }
    
    #check if more than one row existent in selv, since at least two matching occseq entries are necessary for rlm
    if(length(which(selv))>1){
      #create function and named vector selector as string
      
      f <- paste("x[", which(selv)[1], ",]~x[", which(selv)[2], ",]", sep="")
      if(length(which(selv)) > 2){
        for(repn in which(selv)[c(-1,-2)]){
          f <- paste(f, "+x[", repn ,",]", sep="")
          #g <- c(g, paste("x[",repn,", ]", sep=""))
        }
      }
      #should NOT add plus 0, as entries cannot be expected to correlate through 0
      #f <- paste(f, "+0", sep="")
      
      g <- c("x[1, ]", "x[2, ]")
      if(nrow(x)>2){
        for(repn in 3:nrow(x)){
          g <- c(g, paste("x[",repn,", ]", sep=""))
        }
      }
      
      #perform rlm calculation and return data table
      #importantly, cannot call value by name if rlm models only two entries, then have to construct vector manyally with slope in between
      if(length(which(selv))==2){
        if(skiperror){
          return(data.table(slope = as.numeric(c(rep(NA, (which(selv)[2]-1)),
                                                 tryCatch(coefficients(summary(do.call("rlm", list(as.formula(f), maxit=1000))))[-1,"Value"],
                                                          error=function(e){print("RLM error occured, returning NA; see FAQ");NA}),
                                                 rep(NA, (length(selv)-which(selv)[2])))), calc = selv))
        } else {
          return(data.table(slope = as.numeric(c(rep(NA, (which(selv)[2]-1)), coefficients(summary(do.call("rlm", list(as.formula(f), maxit=1000))))[-1,"Value"],
                                                 rep(NA, (length(selv)-which(selv)[2])))), calc = selv))
        }
        
      } else {
        if(skiperror){
          return(data.table(slope = tryCatch(coefficients(summary(do.call("rlm", list(as.formula(f), maxit=1000))))[,"Value"][g],
                                             error=function(e){print("RLM error occured, returning NA; see FAQ");as.numeric(rep(NA, nrow(x)))}), calc = selv))
        } else {
          return(data.table(slope = as.numeric(coefficients(summary(do.call("rlm", list(as.formula(f), maxit=1000))))[,"Value"][g]), calc = selv))
        }
        
      }
      
    } else {
      #if only one row, return NA data table
      return(data.table(slope = as.numeric(rep(NA, nrow(x))), calc = rep(F, nrow(x))))
    }
    
  }
  
  stoichiometry <- function(pepi, calc, cond){
    #check if enough calc entries present, if yes execute occ calculation
    stoi <- setNames(as.data.table(rep(list(rep(as.numeric(NA), length(calc))), length(cond))), gsub("^(.*)$", "Occ_\\1", cond))
    
    if(length(which(calc))>1){
      #first entry, calculate occ and write into first medlist evd group
      #for the first entry: except target entry, times own slope; two slopes plus, one slope minus
      stoi[which(calc)[1], ] <- pepi[, .SD[which(calc)[1],-1]-(.SD[which(calc)[2],-1]*as.numeric(.SD[which(calc)[2],1]))]
      if(length(which(calc))>2){
        for(repn in which(calc)[-c(1,2)]){stoi[which(calc)[1], ] <- stoi[which(calc)[1], ] - pepi[, .SD[repn,-1]*as.numeric(.SD[repn,1])]}
      }
      stoi[which(calc)[1], gsub("^(.*)$", "Occ_\\1", cond) := pepi[, .SD[which(calc)[1],-1]] / .SD]
      
      #for all subsequent entries, calculate occ and write into subsequent evd groups
      #for other entries: except target entry, times own slope, divided by target slope (y entry exception: only divided by target); two slopes plus, one slope minus
      for(repn2 in which(calc)[-1] ){
        stoi[repn2, ] <- pepi[, .SD[repn2,-1]-(.SD[which(calc)[1],-1]/as.numeric(.SD[repn2,1]))]
        
        if(length(which(calc))>2){
          for(repn in which(calc)[-c(1,grep(repn2, which(calc)))]){stoi[repn2, ] <- stoi[repn2, ] + pepi[, .SD[repn,-1]*as.numeric(.SD[repn,1])/as.numeric(.SD[repn2,1])]}
        }
        
        #calculate occupancies and write results into medlist for subsequent entries
        #multiply with if testing for slope being NA even though calc is TRUE, which then needs to set all stoi to NA for this entry
        stoi[repn2, gsub("^(.*)$", "Occ_\\1", cond) := pepi[, .SD[repn2,-1]*if(is.na(as.numeric(.SD[repn2,1]))){NA}else{1}] / .SD]
      }
      return(stoi)
    } else {
      #if not enough calc entries present, return NA data table
      return(stoi)
    }
  }
  
  stoichdt <- function(temp, cond.names){
    #calculate rlm and store as model, then report back occupancies
    temp <- copy(cbind(temp[, .SD, by=PTM_stoich_key], temp[, rlmfun(as.matrix(.SD), skiperror=T, multiphos=2, PTMnum=as.numeric(PTM_0_num)),
                                                            .SDcols = cond.names, by=PTM_stoich_key][,-1]))
    
    #calculate stoichiometry values for each occgrp
    temp <- copy(cbind(temp[, .SD, by=PTM_stoich_key], temp[, stoichiometry(.SD, calc, cond.names), .SDcols=c("slope", cond.names), by=PTM_stoich_key][,-1]))
    
    #delete illegal
    for (col in gsub("^(.*)$", "Occ_\\1", cond.names)){ temp[get(col)<0 | get(col)>1, (col) := NA] }
    
    #calculate number of PTM_stoich_key sequences per PTM_stoich_key grp and number of PTMs per PTM_stoich_key grp
    temp[calc==T, PTM_stoich_key_num := .N, by=PTM_stoich_key]
    temp[calc==T, PTM_stoich_key_PTMnum := paste(sort(unique(PTM_0_num)), collapse=";"), by=PTM_stoich_key]
    
    #return data table to top level
    return(temp)
  }
  
  #initiate parallel processing
  cl <- makeCluster(par.CPU)
  registerDoParallel(cl)
  
  #parallel stoichiometry calculation
  data <- copy(rbindlist(foreach(testp=subdt,.packages=c("data.table", "MASS")) %dopar% {
    stoichdt(temp=testp, cond.names=cond.names)
  }))
  
  #close parallel processing
  stopCluster(cl)
  
  cond.names <- c(cond.names, gsub("^(.*)$", "Occ_\\1", cond.names))
  
  print("Done: Stoichiometry calculation")
}





#write res into main including row names
maind <- as.data.frame(data[, .SD, .SDcols = cond.names])
rownames(maind) <- paste("Row.", seq(1,nrow(maind),1), sep="")
mdata@main <- maind


#prepare writing additional columns into annotCols
keep <- character()

#main columns
if(par.level.col %in% c(0,1,3)){
  keep <- c(keep, else.name, "PTM_0_num", "PTM_group", "PTM_collapse_key", "PTM_collapse_key_num", "PTM_localization")
} else if(par.level.col == 2){
  keep <- c(keep, else.name, "PTM_0_num", "PTM_group", "PTM_collapse_key", "PTM_collapse_key_num")
}


if(if(exists("par.level.stoich")){par.level.stoich == 0}else{F}){
  #collapse and stoichiometry
  keep <- c(keep, "PTM_stoich_key", "PTM_stoich_key_num", "PTM_stoich_key_PTMnum")
  
} else if( par.level == 0 & if(exists("par.FASTA.file")){!identical(par.FASTA.file, character(0))}else{F} ){
  #collapse and PTM motif sequence
  keep <- c(keep, "PTM_seq", "PTM_0_aa")
  
} else if(par.level == 0){
  #collapse site-level only
  keep <- c(keep, "PTM_0_aa")
  
} else {
  #collapse only
  
}

#if grouping performed, delete par.cond.col from keep list, since does not exist anymore
if(par.cond == 0){
  keep <- setdiff(keep, par.cond.col)
}

#write additional columns into annotCols
c.annot <- as.data.frame(data[, .SD, .SDcols = keep])
rownames(c.annot) <- paste("Row.", seq(1,nrow(c.annot),1), sep="")
mdata@annotCols <- c.annot

# if(par.level.col %in% c(0,1)){
#   
#   if(if(exists("par.level.stoich")){par.level.stoich == 0}else{F}){
#     #collapse and stoichiometry
#     c.annot <- as.data.frame(data[, .SD, .SDcols = setdiff(c(else.name, "PTM_0_num", "PTM_group", "PTM_collapse_key", "PTM_collapse_key_num", "PTM_localization",
#                                                              "PTM_stoich_key", "PTM_stoich_key_num", "PTM_stoich_key_PTMnum"), par.cond.col)])
#     rownames(c.annot) <- paste("Row.", seq(1,nrow(c.annot),1), sep="")
#     mdata@annotCols <- c.annot
#   } else if( par.level == 0 & if(exists("par.FASTA.file")){!identical(par.FASTA.file, character(0))}else{F} ){
#     #collapse and PTM motif sequence
#     c.annot <- as.data.frame(data[, .SD, .SDcols = setdiff(c(else.name, "PTM_0_num", "PTM_group", "PTM_collapse_key", "PTM_collapse_key_num", "PTM_localization",
#                                                              "PTM_seq", "PTM_0_aa"), par.cond.col)])
#     rownames(c.annot) <- paste("Row.", seq(1,nrow(c.annot),1), sep="")
#     mdata@annotCols <- c.annot
#   } else if(par.level == 0){
#     #collapse site-level only
#     c.annot <- as.data.frame(data[, .SD, .SDcols = setdiff(c(else.name, "PTM_0_num", "PTM_group", "PTM_collapse_key", "PTM_collapse_key_num", "PTM_localization",
#                                                              "PTM_0_aa"), par.cond.col)])
#     rownames(c.annot) <- paste("Row.", seq(1,nrow(c.annot),1), sep="")
#     mdata@annotCols <- c.annot
#   } else {
#     #collapse only
#     c.annot <- as.data.frame(data[, .SD, .SDcols = setdiff(c(else.name, "PTM_0_num", "PTM_group", "PTM_collapse_key", "PTM_collapse_key_num", "PTM_localization"),
#                                                            par.cond.col)])
#     rownames(c.annot) <- paste("Row.", seq(1,nrow(c.annot),1), sep="")
#     mdata@annotCols <- c.annot
#   }
# }
# if(par.level.col == 2){
#   
#   if(if(exists("par.level.stoich")){par.level.stoich == 0}else{F}){
#     #collapse and stoichiometry
#     c.annot <- as.data.frame(data[, .SD, .SDcols = setdiff(c(else.name, "PTM_0_num", "PTM_group", "PTM_collapse_key", "PTM_collapse_key_num",
#                                                              "PTM_stoich_key", "PTM_stoich_key_num", "PTM_stoich_key_PTMnum"), par.cond.col)])
#     rownames(c.annot) <- paste("Row.", seq(1,nrow(c.annot),1), sep="")
#     mdata@annotCols <- c.annot
#   } else if( par.level == 0 & if(exists("par.FASTA.file")){!identical(par.FASTA.file, character(0))}else{F} ){
#     #collapse and PTM motif sequence
#     c.annot <- as.data.frame(data[, .SD, .SDcols = setdiff(c(else.name, "PTM_0_num", "PTM_group", "PTM_collapse_key", "PTM_collapse_key_num", "PTM_seq",
#                                                              "PTM_0_aa"), par.cond.col)])
#     rownames(c.annot) <- paste("Row.", seq(1,nrow(c.annot),1), sep="")
#     mdata@annotCols <- c.annot
#   } else if(par.level == 0){
#     #collapse only
#     c.annot <- as.data.frame(data[, .SD, .SDcols = setdiff(c(else.name, "PTM_0_num", "PTM_group", "PTM_collapse_key", "PTM_collapse_key_num", "PTM_0_aa"),
#                                                            par.cond.col)])
#     rownames(c.annot) <- paste("Row.", seq(1,nrow(c.annot),1), sep="")
#     mdata@annotCols <- c.annot
#   } else {
#     #collapse site-level only
#     c.annot <- as.data.frame(data[, .SD, .SDcols = setdiff(c(else.name, "PTM_0_num", "PTM_group", "PTM_collapse_key", "PTM_collapse_key_num"),
#                                                            par.cond.col)])
#     rownames(c.annot) <- paste("Row.", seq(1,nrow(c.annot),1), sep="")
#     mdata@annotCols <- c.annot
#   }
# }


#write descriptions
#mdata@description <- c(mname, setdiff(names(res), mname))

#write output file
write.perseus(mdata, outFile)
