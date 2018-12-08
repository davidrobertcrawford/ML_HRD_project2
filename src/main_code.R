#!/usr/bin/env Rscript

# "main_code.R"



# Preamble ----------------------------------------------------------------

rm(list=ls())
gc()

library(stringr)
library(caret)
library(glmnet)
library(pROC)
library(randomForest)
library(ggplot2)
library(gridExtra)
library(data.table)


# Set number of runs for run3ML
num_runs = 20



args = commandArgs(trailingOnly=TRUE)

add_feat = fread(args[1])
#add_feat = fread("additional_features.ICGC-BRCA-EU_BRCA_Davies2017.tsv")
counts_brca = fread(args[2])
#counts_brca = fread("counts.ICGC-BRCA-EU_BRCA_22.SBS-96.tsv")
sample_key = fread(args[3])
#sample_key = fread("samples.ICGC-BRCA-EU_BRCA_22.tsv")
raw_counts = fread(args[4])
#raw_counts = fread("counts.ICGC-BRCA-EU_BRCA_22.Letouze2017.tsv")







set.seed(2048)

source("run3ML.R")


# Functions ---------------------------------------------------------------

log_col = function(the_vector){
    the_vector1 = log(the_vector+1)
    return(the_vector1)
}

norm_col = function(the_vector){
    the_mean = mean(the_vector)
    the_sd = sd(the_vector)
    the_vector1 = the_vector-the_mean
    the_vector2 = the_vector1/the_sd
    return(the_vector2)
}


# Load files --------------------------------------------------------------

# Harmonize sample ids ----------------------------------------------------

#length(unique(sample_key$Patient))
#length(unique(sample_key$Sample))
#length(unique(sample_key$`Alternative Sample Name`))

#colnames(sample_key)
setnames(sample_key, "Sample", "Sample_id_1")
setnames(sample_key, "Alternative Sample Name", "Sample_id_2")

#colnames(add_feat)
setnames(add_feat, "Sample", "Sample_id_2")

#colnames(counts_brca)
setnames(counts_brca, "V1", "Sample_id_1")

#colnames(clin_feat)
#setnames(clin_feat, "Sample", "Sample_id_2")

#colnames(raw_counts)
setnames(raw_counts, "Sample", "Sample_id_1")



# Check for data completeness ---------------------------------------------

#sum(is.na(add_feat))
#sum(is.na(clin_feat))
#sum(is.na(counts_brca))
#sum(is.na(raw_counts))

#summary(add_feat)

#Remove rows w/ NA
add_feat = add_feat[!(is.na(hrd))]
#sum(is.na(add_feat))


# Set uniform ids & set status column -------------------------------------

add_feat = merge(add_feat, sample_key, by="Sample_id_2")
add_feat[,c("Sample_id_2", "Patient"):=NULL]

#Remove samples without methylation data (see Davies et al. online methods)
add_feat = add_feat[add_feat$isUsedForEvaluation==TRUE]

#Subset genome status (should have 371 unique ids)
status_dt = add_feat[,.SD,.SDcols=c("Sample_id_1", "isQuiescentGenomeControl")]

#Add genome status into two other datasets
counts_brca = merge(counts_brca, status_dt, by="Sample_id_1")

#Merge the raw data
raw_dt = merge(counts_brca, raw_counts, by="Sample_id_1")


# Prep data tables --------------------------------------------------------

#dim(raw_dt)
proc_dt = data.table(add_feat)
proc_dt[,c("isUsedForEvaluation"):=NULL]


# Prepare test/train sets -------------------------------------------------

cases = unique(status_dt$Sample_id_1)

cases_false = unique(status_dt[status_dt$isQuiescentGenomeControl==TRUE,
                                Sample_id_1])
cases_true = unique(status_dt[status_dt$isQuiescentGenomeControl==FALSE,
                               Sample_id_1])

count_all = length(cases)
count_false = length(cases_false)
count_true = length(cases_true)

#Get size for holdout
holdout_size = round(count_all*.1)

prepare_tt = function(cases){
    holdout = sample(cases, size = holdout_size, replace=FALSE)
    
    holdout_false = status_dt[Sample_id_1 %in% holdout & 
                                  isQuiescentGenomeControl==FALSE, 
                              Sample_id_1]
    
    holdout_true = status_dt[Sample_id_1 %in% holdout &
                                 isQuiescentGenomeControl==TRUE,
                             Sample_id_1]
   
    holdout_size_each = min(length(holdout_false), length(holdout_true))
    
    holdout_false = sample(holdout_false, size=holdout_size_each,
                           replace=FALSE)
    holdout_true = sample(holdout_true, size=holdout_size_each,
                          replace=FALSE)

    holdout_final = c(holdout_false, holdout_true)
    
    tt_list = status_dt[!(Sample_id_1 %in% holdout), Sample_id_1]
    
    tt_false = status_dt[Sample_id_1 %in% tt_list &
                             isQuiescentGenomeControl==FALSE,
                         Sample_id_1]
    
    tt_true = status_dt[Sample_id_1 %in% tt_list &
                            isQuiescentGenomeControl==TRUE,
                        Sample_id_1]
    
    tt_size_each = min(length(tt_false), length(tt_true))
    
    tt_false = sample(tt_false, size=tt_size_each,
                      replace=FALSE)
    
    tt_true = sample(tt_true, size=tt_size_each,
                     replace=FALSE)
    
    tt_final = c(tt_true, tt_false)
    
    all_final = c(tt_final, holdout_final)
    
    sample_round_data = data.table(status_dt)
    sample_round_data = sample_round_data[Sample_id_1 %in% all_final]
    
    sample_round_data[,isHoldout:=ifelse(Sample_id_1 %in% holdout_final,
                                         TRUE, FALSE)]
    
    sample_round_data[,Sample_id_1:=as.factor(Sample_id_1)]
    
    return(sample_round_data)
}



# Run trials --------------------------------------------------------------



auc_dt = data.table(Set_num = integer(),
                    Raw_or_proc = character(),
                    Normalize = logical(),
                    Log_transform = logical(),
                    AUC_glm = numeric(),
                    AUC_svm = numeric(),
                    AUC_RF = numeric(),
                    Rank_glm = integer(),
                    Rank_svm = integer(),
                    Rank_RF = integer())

for(num in seq(1,num_runs,1)){
    set_1 = prepare_tt(cases)
    results_1rff = run3ML(sample_data=set_1, 
                          raw_or_proc="raw",
                          normalize=FALSE,
                          log_transform=FALSE)
    rank_1rff = rank(-results_1rff)
    
    results_1pff = run3ML(sample_data=set_1, 
                          raw_or_proc="proc",
                          normalize=FALSE,
                          log_transform=FALSE)
    rank_1pff = rank(-results_1pff)
    
    results_1bff = run3ML(sample_data=set_1, 
                          raw_or_proc="both",
                          normalize=FALSE,
                          log_transform=FALSE)
    rank_1bff = rank(-results_1bff)
    
    results_1rtt = run3ML(sample_data=set_1, 
                          raw_or_proc="raw",
                          normalize=TRUE,
                          log_transform=TRUE)
    rank_1rtt = rank(-results_1rtt)
    
    results_1ptt = run3ML(sample_data=set_1, 
                          raw_or_proc="proc",
                          normalize=TRUE,
                          log_transform=TRUE)
    rank_1ptt = rank(-results_1ptt)
    
    results_1btt = run3ML(sample_data=set_1, 
                          raw_or_proc="both",
                          normalize=TRUE,
                          log_transform=TRUE)
    rank_1btt = rank(-results_1btt)
    
    results_1rtf = run3ML(sample_data=set_1, 
                          raw_or_proc="raw",
                          normalize=TRUE,
                          log_transform=FALSE)
    rank_1rtf = rank(-results_1rtf)
    
    results_1ptf = run3ML(sample_data=set_1, 
                          raw_or_proc="proc",
                          normalize=TRUE,
                          log_transform=FALSE)
    rank_1ptf = rank(-results_1ptf)
    
    results_1btf = run3ML(sample_data=set_1, 
                          raw_or_proc="both",
                          normalize=TRUE,
                          log_transform=FALSE)
    rank_1btf = rank(-results_1btf)
    
    results_1rft = run3ML(sample_data=set_1, 
                          raw_or_proc="raw",
                          normalize=FALSE,
                          log_transform=TRUE)
    rank_1rft = rank(-results_1rft)
    
    results_1pft = run3ML(sample_data=set_1, 
                          raw_or_proc="proc",
                          normalize=FALSE,
                          log_transform=TRUE)
    rank_1pft = rank(-results_1pft)
    
    results_1bft = run3ML(sample_data=set_1, 
                          raw_or_proc="both",
                          normalize=FALSE,
                          log_transform=TRUE)
    rank_1bft = rank(-results_1bft)
    
    
    row_rff = list(num, "raw", FALSE, FALSE, results_1rff[1], 
                   results_1rff[2],
                   results_1rff[3],
                   rank_1rff[1],
                   rank_1rff[2],
                   rank_1rff[3])
    row_pff= list(num, "proc", FALSE, FALSE, results_1pff[1], 
                  results_1pff[2],
                  results_1pff[3],
                  rank_1pff[1],
                  rank_1pff[2],
                  rank_1pff[3])
    row_bff= list(num, "both", FALSE, FALSE, results_1bff[1], 
                  results_1bff[2],
                  results_1bff[3],
                  rank_1bff[1],
                  rank_1bff[2],
                  rank_1bff[3])
    
    row_rtt= list(num, "raw", TRUE, TRUE, results_1rtt[1], 
                  results_1rtt[2],
                  results_1rtt[3],
                  rank_1rtt[1],
                  rank_1rtt[2],
                  rank_1rtt[3])
    row_ptt= list(num, "proc", TRUE, TRUE, results_1ptt[1], 
                  results_1ptt[2],
                  results_1ptt[3],
                  rank_1ptt[1],
                  rank_1ptt[2],
                  rank_1ptt[3])
    row_btt= list(num, "both", TRUE, TRUE, results_1btt[1], 
                  results_1btt[2],
                  results_1btt[3],
                  rank_1btt[1],
                  rank_1btt[2],
                  rank_1btt[3])
    
    row_rft= list(num, "raw", FALSE, TRUE, results_1rft[1], 
                  results_1rft[2],
                  results_1rft[3],
                  rank_1rft[1],
                  rank_1rft[2],
                  rank_1rft[3])
    row_pft= list(num, "proc", FALSE, TRUE, results_1pft[1], 
                  results_1pft[2],
                  results_1pft[3],
                  rank_1pft[1],
                  rank_1pft[2],
                  rank_1pft[3])
    row_bft= list(num, "both", FALSE, TRUE, results_1bft[1], 
                  results_1bft[2],
                  results_1bft[3],
                  rank_1bft[1],
                  rank_1bft[2],
                  rank_1bft[3])
    
    
    row_rtf= list(num, "raw", TRUE, FALSE, results_1rtf[1], 
                  results_1rtf[2],
                  results_1rtf[3],
                  rank_1rtf[1],
                  rank_1rtf[2],
                  rank_1rtf[3])
    row_ptf= list(num, "proc", TRUE, FALSE, results_1ptf[1], 
                  results_1ptf[2],
                  results_1ptf[3],
                  rank_1ptf[1],
                  rank_1ptf[2],
                  rank_1ptf[3])
    row_btf= list(num, "both", TRUE, FALSE, results_1btf[1], 
                  results_1btf[2],
                  results_1btf[3],
                  rank_1btf[1],
                  rank_1btf[2],
                  rank_1btf[3])
    
    
    auc_dt = rbind(auc_dt, row_rff, row_pff, row_bff, 
                   row_rtt, row_ptt, row_btt,
                   row_rft, row_pft, row_bft,
                   row_rtf, row_ptf, row_btf)
    
    gc()
}

auc_dt[,Rank_glm:=as.integer(Rank_glm)]
auc_dt[,Rank_svm:=as.integer(Rank_svm)]
auc_dt[,Rank_RF:=as.integer(Rank_RF)]

#file_name = "auc_dt_run3ML_2018_12_06.RData"

#save(auc_dt, file=file_name)


rank_all = round(c(mean(auc_dt$Rank_RF),
                   mean(auc_dt$Rank_svm),
                   mean(auc_dt$Rank_glm)), 2)
auc_all = round(c(mean(auc_dt$AUC_RF), mean(auc_dt$AUC_svm), 
                  mean(auc_dt$AUC_glm)), 3)

# Subset means ------------------------------------------------------------

mean(auc_dt$AUC_glm)
mean(auc_dt$AUC_svm)
mean(auc_dt$AUC_RF)

median(auc_dt$AUC_glm)
median(auc_dt$AUC_svm)
median(auc_dt$AUC_RF)

# Log t/f
auc_dt_LOG_t = auc_dt[Log_transform==TRUE]
auc_dt_LOG_f = auc_dt[Log_transform==FALSE]

rank_log_t = round(c(mean(auc_dt_LOG_t$Rank_RF),
                   mean(auc_dt_LOG_t$Rank_svm),
                   mean(auc_dt_LOG_t$Rank_glm)), 2)
auc_log_t = round(c(mean(auc_dt_LOG_t$AUC_RF), 
                    mean(auc_dt_LOG_t$AUC_svm), 
                  mean(auc_dt_LOG_t$AUC_glm)), 3)

rank_log_f = round(c(mean(auc_dt_LOG_f$Rank_RF),
                    mean(auc_dt_LOG_f$Rank_svm),
                    mean(auc_dt_LOG_f$Rank_glm)), 2)
auc_log_f = round(c(mean(auc_dt_LOG_f$AUC_RF), 
                   mean(auc_dt_LOG_f$AUC_svm), 
                   mean(auc_dt_LOG_f$AUC_glm)), 3)


auc_dt_n_t = auc_dt[Normalize==TRUE]
auc_dt_n_f = auc_dt[Normalize==FALSE]

rank_n_t = round(c(mean(auc_dt_n_t$Rank_RF),
                   mean(auc_dt_n_t$Rank_svm),
                   mean(auc_dt_n_t$Rank_glm)), 2)
auc_n_t = round(c(mean(auc_dt_n_t$AUC_RF), 
                  mean(auc_dt_n_t$AUC_svm), 
                  mean(auc_dt_n_t$AUC_glm)), 3)

rank_n_f = round(c(mean(auc_dt_n_f$Rank_RF),
                   mean(auc_dt_n_f$Rank_svm),
                   mean(auc_dt_n_f$Rank_glm)), 2)
auc_n_f = round(c(mean(auc_dt_n_f$AUC_RF), 
                  mean(auc_dt_n_f$AUC_svm), 
                  mean(auc_dt_n_f$AUC_glm)), 3)





auc_dt_both = auc_dt[Raw_or_proc=="both"]

rank_both = round(c(mean(auc_dt_both$Rank_RF),
                   mean(auc_dt_both$Rank_svm),
                   mean(auc_dt_both$Rank_glm)), 2)
auc_both = round(c(mean(auc_dt_both$AUC_RF), 
                  mean(auc_dt_both$AUC_svm), 
                  mean(auc_dt_both$AUC_glm)), 3)


auc_dt_raw = auc_dt[Raw_or_proc=="raw"]

rank_raw = round(c(mean(auc_dt_raw$Rank_RF),
                    mean(auc_dt_raw$Rank_svm),
                    mean(auc_dt_raw$Rank_glm)), 2)
auc_raw = round(c(mean(auc_dt_raw$AUC_RF), 
                   mean(auc_dt_raw$AUC_svm), 
                   mean(auc_dt_raw$AUC_glm)), 3)



auc_dt_proc = auc_dt[Raw_or_proc=="proc"]

rank_proc = round(c(mean(auc_dt_proc$Rank_RF),
                    mean(auc_dt_proc$Rank_svm),
                    mean(auc_dt_proc$Rank_glm)), 2)
auc_proc = round(c(mean(auc_dt_proc$AUC_RF), 
                   mean(auc_dt_proc$AUC_svm), 
                   mean(auc_dt_proc$AUC_glm)), 3)




auc_dt_ff = auc_dt[Normalize==FALSE & Log_transform==FALSE]

rank_ff = round(c(mean(auc_dt_ff$Rank_RF),
                    mean(auc_dt_ff$Rank_svm),
                    mean(auc_dt_ff$Rank_glm)), 2)
auc_ff = round(c(mean(auc_dt_ff$AUC_RF), 
                   mean(auc_dt_ff$AUC_svm), 
                   mean(auc_dt_ff$AUC_glm)), 3)

auc_dt_tt = auc_dt[Normalize==TRUE & Log_transform==TRUE]

rank_tt = round(c(mean(auc_dt_tt$Rank_RF),
                    mean(auc_dt_tt$Rank_svm),
                    mean(auc_dt_tt$Rank_glm)), 2)
auc_tt = round(c(mean(auc_dt_tt$AUC_RF), 
                   mean(auc_dt_tt$AUC_svm), 
                   mean(auc_dt_tt$AUC_glm)), 3)

auc_dt_tf = auc_dt[Normalize==TRUE & Log_transform==FALSE]

rank_tf = round(c(mean(auc_dt_tf$Rank_RF),
                    mean(auc_dt_tf$Rank_svm),
                    mean(auc_dt_tf$Rank_glm)), 2)
auc_tf = round(c(mean(auc_dt_tf$AUC_RF), 
                   mean(auc_dt_tf$AUC_svm), 
                   mean(auc_dt_tf$AUC_glm)), 3)

auc_dt_ft = auc_dt[Normalize==FALSE & Log_transform==TRUE]

rank_ft = round(c(mean(auc_dt_ft$Rank_RF),
                    mean(auc_dt_ft$Rank_svm),
                    mean(auc_dt_ft$Rank_glm)), 2)
auc_ft = round(c(mean(auc_dt_ft$AUC_RF), 
                   mean(auc_dt_ft$AUC_svm), 
                   mean(auc_dt_ft$AUC_glm)), 3)



rank_table = data.table(Model = c("RF", "SVM", "GLM"),
                     All = rank_all,
                     Raw = rank_raw,
                     Proc = rank_proc,
                     Both = rank_both,
                     Log_Norm_TT = rank_tt,
                     Log_Norm_TF = rank_tf,
                     Log_Norm_FT = rank_ft,
                     Log_Norm_FF = rank_ff)

auc_table = data.table(Model = c("RF", "SVM", "GLM"),
                     All = auc_all,
                     Raw = auc_raw,
                     Proc = auc_proc,
                     Both = auc_both,
                     Log_Norm_TT = auc_tt,
                     Log_Norm_TF = auc_tf,
                     Log_Norm_FT = auc_ft,
                     Log_Norm_FF = auc_ff)


# GGPLOT ------------------------------------------------------------------
auc_cols = colnames(auc_dt)

auc_auc = auc_dt[,.SD,.SDcols=auc_cols[1:7]]

auc_melt = melt(auc_auc, id.vars=auc_cols[1:4])

setnames(auc_melt, "variable", "ML_method")
setnames(auc_melt, "value", "AUC")

the_cols = c(auc_cols[1:4], auc_cols[8:10])

auc_rank = data.table(auc_dt)
auc_rank[,AUC_glm:=NULL]
auc_rank[,AUC_RF:=NULL]
auc_rank[,AUC_svm:=NULL]

#auc_rank = auc_dt[,.SD,.SDcols=the_cols]
rank_melt = melt(auc_rank, id.vars=c("Raw_or_proc", "Log_transform", 
                                     "Normalize", "Set_num"))
setnames(rank_melt, "variable", "ML_method")
setnames(rank_melt, "value", "Rank")
rank_melt[,Rank:=as.integer(Rank)]
#rank_melt[,Raw_or_proc:=as.factor("Raw_or_proc")]
#rank_melt[,Log_transform:=as.factor("Log_transform")]
#rank_melt[,Normalize:=as.factor("Normalize")]
#rank_melt[,ML_method:=as.factor("ML_method")]

#colnames(auc_melt)

#auc_melt[,Raw_or_proc:=as.factor("Raw_or_proc")]
#auc_melt[,Log_transform:=as.factor("Log_transform")]
#auc_melt[,Normalize:=as.factor("Normalize")]
#auc_melt[,ML_method:=as.factor("ML_method")]

plot_1 = ggplot(data=auc_melt)+
    geom_violin(aes(x=ML_method, y=AUC, color=ML_method, fill=ML_method))+
    theme_bw()+
    ggtitle("AUC by Machine Learning Method")+
    theme(plot.title = element_text(hjust=.5))

plot_2 = ggplot(data=auc_melt)+
    geom_violin(aes(x=ML_method, y=AUC, color=ML_method, fill=ML_method))+
    theme_bw()+
    ggtitle("AUC by Machine Learning Method & Feature Set (both, proc, raw)")+
    theme(plot.title = element_text(hjust=.5))+
    facet_grid(.~Raw_or_proc) #facet_wrap(~Raw_or_proc, ncol=3)#

plot_3 = ggplot(data=auc_melt)+
    geom_violin(aes(x=ML_method, y=AUC, color=ML_method, fill=ML_method))+
    theme_bw()+
    ggtitle("AUC by Machine Learning Method & Normalize (T/F) & Log_transform (T/F)")+
    theme(plot.title = element_text(hjust=.5))+
    facet_grid(.~Log_transform+Normalize)




plot_4 = ggplot(data=rank_melt)+
    geom_histogram(aes(x=Rank, color=ML_method, fill=ML_method), bins=3)+
    theme_bw()+
    ggtitle("Rank by Machine Learning Method & Feature Set (both, proc, raw)")+
    theme(plot.title = element_text(hjust=.5),
          legend.position="none")+
    facet_grid(Raw_or_proc~ML_method) # 


jpeg(filename=args[5], width=600, height=400, units="px", 
     pointsize=36, quality=100)
plot_1
dev.off()

jpeg(filename=args[6], width=600, height=400, units="px", 
     pointsize=36, quality=100)
plot_2
dev.off()

jpeg(filename=args[7], width=600, height=400, units="px", 
     pointsize=36, quality=100)
plot_3
dev.off()

jpeg(filename=args[8], width=600, height=400, units="px", 
     pointsize=36, quality=100)
plot_4
dev.off()



# Tables ------------------------------------------------------------------





# Tables with AUROC mean, median, etc.

jpeg(filename=args[9], width=600, height=400, units="px", 
     pointsize=36, quality=100)
grid.table(auc_table)
dev.off()


jpeg(filename=args[10], width=600, height=400, units="px", 
     pointsize=36, quality=100)
grid.table(rank_table)
dev.off()