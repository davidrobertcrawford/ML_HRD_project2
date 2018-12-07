
#!/usr/bin/env Rscript

### "addeR.R"
# Produce multiple ROC AUC plot from Davies et al. 2017 (Fig4a)

args = commandArgs(trailingOnly=TRUE)

set.seed(3141)

library(tidyverse)
library(data.table)
library(glmnet)
library(pROC)
library(gridExtra)



add_feat = fread(args[1])
clin_brca = fread(args[2])





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

#Source outer_fold function
source("outer_fold.R")

#Source simple_train function
source("simple_train.R")

#Source run_training_cv function
source("run_training_cv.R")



# Explore Nonzero_counts --------------------------------------------------
# Run 20 x cv_glmnet (100 outer, 10 inner folds), once for each combination of arguments
r1 = run_training_cv(the_measure="mae", to_balance=TRUE, to_normalize=TRUE)
r2 = run_training_cv(the_measure="mae", to_balance=FALSE, to_normalize=TRUE)
r3 = run_training_cv(the_measure="mae", to_balance=TRUE, to_normalize=FALSE)
r4 = run_training_cv(the_measure="mae", to_balance=FALSE, to_normalize=FALSE)

r5 = run_training_cv(the_measure="mse", to_balance=TRUE, to_normalize=TRUE)
r6 = run_training_cv(the_measure="mse", to_balance=FALSE, to_normalize=TRUE)
r7 = run_training_cv(the_measure="mse", to_balance=TRUE, to_normalize=FALSE)
r8 = run_training_cv(the_measure="mse", to_balance=FALSE, to_normalize=FALSE)

r9 = run_training_cv(the_measure="auc", to_balance=TRUE, to_normalize=TRUE)
r10 = run_training_cv(the_measure="auc", to_balance=FALSE, to_normalize=TRUE)
r11 = run_training_cv(the_measure="auc", to_balance=TRUE, to_normalize=FALSE)
r12 = run_training_cv(the_measure="auc", to_balance=FALSE, to_normalize=FALSE)

r13 = run_training_cv(the_measure="deviance", to_balance=TRUE, to_normalize=TRUE)
r14 = run_training_cv(the_measure="deviance", to_balance=FALSE, to_normalize=TRUE)
r15 = run_training_cv(the_measure="deviance", to_balance=TRUE, to_normalize=FALSE)
r16 = run_training_cv(the_measure="deviance", to_balance=FALSE, to_normalize=FALSE)

r17 = run_training_cv(the_measure="class", to_balance=TRUE, to_normalize=TRUE)
r18 = run_training_cv(the_measure="class", to_balance=FALSE, to_normalize=TRUE)
r19 = run_training_cv(the_measure="class", to_balance=TRUE, to_normalize=FALSE)
r20 = run_training_cv(the_measure="class", to_balance=FALSE, to_normalize=FALSE)

# Combine results
results_20 = rbind(r1,r2,r3,r4,
                   r5,r6,r7,r8,
                   r9,r10,r11,r12,
                   r13,r14,r15,r16,
                   r17,r18,r19,r20)

# Plot the nonzero counts
results_nonzero = unique(results_20[,.SD,.SDcols=c("Variables", "Nonzero_count", "Settings")])
results_nonzero = results_nonzero[!(Variables=="(Intercept)")]

nonzero_plot = ggplot(data=results_nonzero)+
    geom_jitter(aes(x=Variables, y=Nonzero_count, color=Variables),
                size=3,
                alpha=.75,
                width=.25,
                show.legend=FALSE)+
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust=0.4, hjust=1))

# Save plot to jpg
jpeg(filename=args[4], width=600, height=400, units="px", pointsize=36, quality=100)
nonzero_plot
dev.off()

# Build table of nonzero counts
results_nonzero[,NZ_Max:=max(Nonzero_count), by=Variables]
results_nonzero = results_nonzero[NZ_Max>0]
results_nonzero[,NZ_Max:=NULL]

results_cast = dcast(results_nonzero, Settings~Variables, value.var="Nonzero_count")

# Save table to jpg
jpeg(filename=args[5], width=800, height=500, units="px", pointsize=36, quality=100)
grid.table(results_cast)
dev.off()


# Plot Values to Find Settings --------------------------------------------
# Make deep copy of results table
results_20_slim = data.table(results_20)

results_20_slim[,Value_Max:=max(value), by=Variables]

results_20_slim[,value:=NULL]
results_20_slim[,variable:=NULL]
results_20_slim = unique(results_20_slim)

# Remove variables (genomic features) with no nonzero counts
results_20_slim=results_20_slim[Value_Max>0 | Variables=="(Intercept)"]

# Build plot of features with and without internal standardization
standardize_plot = ggplot(data = results_20_slim)+
    geom_histogram(aes(x=Value_Avg), bins=100, alpha=.75)+
    facet_grid(Variables~Standardize, scales="free")+
    theme_bw()+
    theme(strip.text.y = element_text(angle=0))+
    theme(axis.text.y=element_text(angle=0, vjust=0.4, hjust=1))+
    geom_vline(xintercept=0, lty=3)

# Save plot to jpg
jpeg(filename=args[6], width=600, height=800, units="px", pointsize=36, quality=100)
standardize_plot
dev.off()

# Build table of counts with standardization==FALSE
results_20_sFalse = results_20_slim[Standardize==FALSE]
results_20_sFalse = results_20_sFalse[Value_Max>0 | Variables=="(Intercept)"]

balance_plot = ggplot(data = results_20_sFalse)+
    geom_histogram(aes(x=Value_Avg), bins=100, alpha=.75)+
    facet_grid(Variables~Balance)+
    theme_bw()+
    theme(strip.text.y = element_text(angle=0))+
    theme(axis.text.y=element_text(angle=0, vjust=0.4, hjust=1))+
    geom_vline(xintercept=0, lty=3)

# Save plot to jpg
jpeg(filename=args[7], width=600, height=800, units="px", pointsize=36, quality=100)
balance_plot
dev.off()


# Build plot comparing with respect to measures
measure_plot = ggplot(data = results_20_sFalse)+
    geom_histogram(aes(x=Value_Avg), bins=100)+
    facet_grid(Variables~Measure, scales="free")+
    theme_bw()+
    theme(strip.text.y = element_text(angle=0))+
    theme(axis.text.y=element_text(angle=0, vjust=0.4, hjust=1))+
    geom_vline(xintercept=0, lty=3)

# Save plot to jpg
jpeg(filename=args[8], width=600, height=800, units="px", pointsize=36, quality=100)
measure_plot
dev.off()

# Build table plot of nonzero counts
results_nonzero = unique(results_20_sFalse[,.SD,.SDcols=c("Variables", "Nonzero_count", "Settings")])
results_nonzero = results_nonzero[!(Variables=="(Intercept)")]
results_nonzero[,NZ_Max:=max(Nonzero_count), by=Variables]
results_nonzero = results_nonzero[NZ_Max>0]
results_nonzero[,NZ_Max:=NULL]

results_cast = dcast(results_nonzero, Settings~Variables, value.var="Nonzero_count")

# Save table plot to jpg
jpeg(filename=args[9], width=700, height=250, units="px", pointsize=36, quality=100)
grid.table(results_cast)
dev.off()

# Build comparison table for their supp.12 and 14
# For table 12 (nonzero counts)
comp_table = fread("davies_fig_12_nonzero.txt")
results_nonzero_final = results_nonzero[Settings=="mae_sFALSE_bFALSE"]
results_nonzero_final[,Settings:=NULL]
comp_table = merge(comp_table, results_nonzero_final, by.x="Genomic Feature", by.y="Variables")

# Set to descending order according to counts from Davies et al.
setorderv(comp_table, cols="Nonzero Count (Davies)", order=-1)

setnames(comp_table, "Nonzero_count", "Nonzero Count (Crawford)")

# Build table and save to jpg
jpeg(filename=args[10], width=500, height=250, units="px", pointsize=36, quality=100)
grid.table(comp_table)
dev.off()

# For table 14 (mean weights)
weight_table = fread("davies_fig_14_weights.txt")
our_weights = results_20_sFalse[Settings=="mae_sFALSE_bFALSE"]
our_weights = our_weights[,.SD,.SDcols=c("Variables", "Value_Avg")]
setnames(our_weights, "Variables", "Genomic Feature")
setnames(our_weights, "Value_Avg", "Mean Weight (Crawford)")
weight_combined = merge(weight_table, our_weights, by="Genomic Feature")

#Set to descending order according to weights from Davies et al.
setorderv(weight_combined, cols="Mean Weight (Davies)", order=-1)

# Build table and save to jpg
jpeg(filename=args[11], width=450, height=200, units="px", pointsize=36, quality=100)
grid.table(weight_combined)
dev.off()









# AUC Plots ---------------------------------------------------------------

# Train on 311 cases to get model for final AUC curves

# Build training dataset
glm_input = add_feat[isUsedForEvaluation==TRUE]
input_positive = unique(clin_brca[Loss_of_Alternative_Allele=="Yes",Sample])
glm_input = glm_input[Sample %in% input_positive | isQuiescentGenomeControl==TRUE]
glm_input[,Case:=as.factor(ifelse(isQuiescentGenomeControl==TRUE, 0, 1))]
glm_response = glm_input$Case

glm_input_trim = data.table(glm_input)
glm_input_trim[,isQuiescentGenomeControl:=NULL]
glm_input_trim[,isUsedForEvaluation:=NULL]
glm_input_trim[,Sample:=NULL]
glm_input_trim[,Case:=NULL]

the_names = colnames(glm_input_trim)
glm_input_trim[,(the_names):=lapply(.SD,log_col),.SDcols=the_names]
glm_input_trim[,(the_names):=lapply(.SD,norm_col),.SDcols=the_names]

glm_mat = as.matrix(glm_input_trim)
glm_input_trim[,Case:=glm_response]

# Train model
the_dt = glm_input_trim
#Divide dataset into cases
the_dt_0 = the_dt[Case==0]
the_dt_1 = the_dt[Case==1]

#Count rows in each case dt
the_rows_0 = nrow(the_dt_0)
the_rows_1 = nrow(the_dt_1)

#Assign to train or test (holdout)
#Get size for holdout
holdout_0 = round(the_rows_0*.1)
holdout_1 = round(the_rows_1*.1)

#Select holdout rows
holdout_rows_0 = sample.int(the_rows_0, size=holdout_0, replace=FALSE)
holdout_rows_1 = sample.int(the_rows_1, size=holdout_1, replace=FALSE)

#Assign row number column
the_dt_0[,the_row:=1]
the_dt_0[,the_row:=cumsum(the_row)]
the_dt_1[,the_row:=1]
the_dt_1[,the_row:=cumsum(the_row)]

#Assign Shuffle according to holdouts
the_dt_0[,Shuffle:=ifelse(the_row %in% holdout_rows_0, 1, 0)]
the_dt_1[,Shuffle:=ifelse(the_row %in% holdout_rows_1, 1, 0)]

#Get rid of the_row
the_dt_0[,the_row:=NULL]
the_dt_1[,the_row:=NULL]

the_train_0 = the_dt_0[Shuffle==0]
the_train_1 = the_dt_1[Shuffle==0]

the_test_0 = the_dt_0[Shuffle==1]
the_test_1 = the_dt_1[Shuffle==1]

#Recombine dts
the_train = rbind(the_train_0, the_train_1)
the_test = rbind(the_test_0, the_test_1)

the_train[,Shuffle:=NULL]
the_test[,Shuffle:=NULL]

train_response = the_train$Case
test_response = the_test$Case

the_train[,Case:=NULL]
the_test[,Case:=NULL]

#Convert NaN to 0 [this results from subsamples w/ no nonzero for variable]
for (j in seq_len(ncol(the_train))){
    set(the_train,which(is.nan(the_train[[j]])),j,0)
}

for (j in seq_len(ncol(the_test))){
    set(the_test,which(is.nan(the_test[[j]])),j,0)
}

train_mat = as.matrix(the_train)
test_mat = as.matrix(the_test)

sum(is.na(train_mat))
sum(is.na(test_mat))

train_result = cv.glmnet(x = train_mat,
                         y = train_response,
                         nfolds=10,
                         lower.limits=0,
                         type.measure = "mae",
                         family="binomial",
                         standardize=FALSE)

the_fit = train_result$glmnet.fit

#Build test set for final AUC curves
glm_input = add_feat[isUsedForEvaluation=="TRUE"]
glm_input[,Case:=as.factor(ifelse(isQuiescentGenomeControl==TRUE, 0, 1))]
glm_response = glm_input$Case

glm_input_trim = data.table(glm_input)
glm_input_trim[,isQuiescentGenomeControl:=NULL]
glm_input_trim[,isUsedForEvaluation:=NULL]
glm_input_trim[,Sample:=NULL]
glm_input_trim[,Case:=NULL]

# Standardize
the_names = colnames(glm_input_trim)
glm_input_trim[,(the_names):=lapply(.SD,log_col),.SDcols=the_names]
glm_input_trim[,(the_names):=lapply(.SD,norm_col),.SDcols=the_names]

glm_input_trim[,c("SV4", "SV2", "SV1", "SV6", "e.1",
                  "e.2", "e.5", "e.6", "e.13", "e.17",
                  "e.18", "e.20", "e.26", "del.rep.prop",
                  "del.none.prop"):=0]

test_mat = as.matrix(glm_input_trim)

test_371 = predict(the_fit, type="response", newx=test_mat)

temp_dt = data.table(test_371)
the_cols = colnames(temp_dt)
col_num = length(the_cols)
the_pick = the_cols[col_num]
the_result = temp_dt[,.SD,.SDcols=the_pick]
the_result = unlist(the_result)

roc_obj = roc(glm_response, the_result)

# Plot for each of six features considered singly

# e.8
glm_input_trim_e8 = data.table(glm_input_trim)
glm_input_trim_e8[,c("del.mh.prop","hrd","e.3","SV5","SV3",
                     "SV4", "SV2", "SV1", "SV6", "e.1",
                     "e.2", "e.5", "e.6", "e.13", "e.17",
                     "e.18", "e.20", "e.26", "del.rep.prop",
                     "del.none.prop"):=0]

test_mat = as.matrix(glm_input_trim_e8)
test_371 = predict(the_fit, type="response", newx=test_mat)

temp_dt = data.table(test_371)
the_cols = colnames(temp_dt)
col_num = length(the_cols)
the_pick = the_cols[col_num]
the_result = temp_dt[,.SD,.SDcols=the_pick]
the_result = unlist(the_result)

roc_obj_e8= roc(glm_response, the_result)

# hrd
glm_input_trim_hrd = data.table(glm_input_trim)
glm_input_trim_hrd[,c("del.mh.prop","e.8","e.3","SV5","SV3",
                      "SV4", "SV2", "SV1", "SV6", "e.1",
                      "e.2", "e.5", "e.6", "e.13", "e.17",
                      "e.18", "e.20", "e.26", "del.rep.prop",
                      "del.none.prop"):=0]

test_mat = as.matrix(glm_input_trim_hrd)
test_371 = predict(the_fit, type="response", newx=test_mat)

temp_dt = data.table(test_371)
the_cols = colnames(temp_dt)
col_num = length(the_cols)
the_pick = the_cols[col_num]
the_result = temp_dt[,.SD,.SDcols=the_pick]
the_result = unlist(the_result)

roc_obj_hrd= roc(glm_response, the_result)


# e.3
glm_input_trim_e3 = data.table(glm_input_trim)
glm_input_trim_e3[,c("del.mh.prop","hrd","e.8","SV5","SV3",
                     "SV4", "SV2", "SV1", "SV6", "e.1",
                     "e.2", "e.5", "e.6", "e.13", "e.17",
                     "e.18", "e.20", "e.26", "del.rep.prop",
                     "del.none.prop"):=0]

test_mat = as.matrix(glm_input_trim_e3)
test_371 = predict(the_fit, type="response", newx=test_mat)

temp_dt = data.table(test_371)
the_cols = colnames(temp_dt)
col_num = length(the_cols)
the_pick = the_cols[col_num]
the_result = temp_dt[,.SD,.SDcols=the_pick]
the_result = unlist(the_result)

roc_obj_e3= roc(glm_response, the_result)




# SV5
glm_input_trim_SV5 = data.table(glm_input_trim)
glm_input_trim_SV5[,c("del.mh.prop","hrd","e.3","e.8","SV3",
                      "SV4", "SV2", "SV1", "SV6", "e.1",
                      "e.2", "e.5", "e.6", "e.13", "e.17",
                      "e.18", "e.20", "e.26", "del.rep.prop",
                      "del.none.prop"):=0]

test_mat = as.matrix(glm_input_trim_SV5)
test_371 = predict(the_fit, type="response", newx=test_mat)

temp_dt = data.table(test_371)
the_cols = colnames(temp_dt)
col_num = length(the_cols)
the_pick = the_cols[col_num]
the_result = temp_dt[,.SD,.SDcols=the_pick]
the_result = unlist(the_result)

roc_obj_SV5= roc(glm_response, the_result)



# SV3
glm_input_trim_SV3 = data.table(glm_input_trim)
glm_input_trim_SV3[,c("e.8","hrd","e.3","SV5","del.mh.prop",
                      "SV4", "SV2", "SV1", "SV6", "e.1",
                      "e.2", "e.5", "e.6", "e.13", "e.17",
                      "e.18", "e.20", "e.26", "del.rep.prop",
                      "del.none.prop"):=0]

test_mat = as.matrix(glm_input_trim_SV3)
test_371 = predict(the_fit, type="response", newx=test_mat)

temp_dt = data.table(test_371)
the_cols = colnames(temp_dt)
col_num = length(the_cols)
the_pick = the_cols[col_num]
the_result = temp_dt[,.SD,.SDcols=the_pick]
the_result = unlist(the_result)

roc_obj_SV3= roc(glm_response, the_result)


# DEL
glm_input_trim_DEL = data.table(glm_input_trim)
glm_input_trim_DEL[,c("e.8","hrd","e.3","SV5","SV3",
                      "SV4", "SV2", "SV1", "SV6", "e.1",
                      "e.2", "e.5", "e.6", "e.13", "e.17",
                      "e.18", "e.20", "e.26", "del.rep.prop",
                      "del.none.prop"):=0]

test_mat = as.matrix(glm_input_trim_DEL)
test_371 = predict(the_fit, type="response", newx=test_mat)

temp_dt = data.table(test_371)
the_cols = colnames(temp_dt)
col_num = length(the_cols)
the_pick = the_cols[col_num]
the_result = temp_dt[,.SD,.SDcols=the_pick]
the_result = unlist(the_result)

roc_obj_DEL= roc(glm_response, the_result)


all_the_roc = list(del.mh.prop = roc_obj_DEL,
                   e.3 = roc_obj_e3,
                   e.8 = roc_obj_e8,
                   ALL = roc_obj,
                   hrd = roc_obj_hrd,
                   SV5 = roc_obj_SV5,
                   SV3 = roc_obj_SV3)

the_auc_plot = ggroc(all_the_roc)

jpeg(filename=args[3], width=600,
    height=600, units="px",
    pointsize=36, quality=100)
the_auc_plot
dev.off()


#Build table of AUC values:
the_auc_full = auc(roc_obj)[1]
the_auc_e8 = auc(roc_obj_e8)[1]
the_auc_hrd = auc(roc_obj_hrd)[1]
the_auc_SV5 = auc(roc_obj_SV5)[1]
the_auc_SV3 = auc(roc_obj_SV3)[1]
the_auc_del = auc(roc_obj_DEL)[1]
the_auc_e3 = auc(roc_obj_e3)[1]

genomic_features = c("ALL", "e.8", "hrd", "e.3", "SV5", "SV3", "del.mh.prop")
auc_values = c(round(the_auc_full,2), round(the_auc_e8,2), round(the_auc_hrd,2), round(the_auc_e3,2),
               round(the_auc_SV5,2), round(the_auc_SV3,2), round(the_auc_del,2))

auc_dt = data.table("Genomic Feature" = genomic_features,
                    "ROC AUC" = auc_values)

# Build table and save to jpg
jpeg(filename=args[12], width=250, height=200, units="px", pointsize=36, quality=100)
grid.table(auc_dt)
dev.off()
