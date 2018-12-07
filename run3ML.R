# "run3ML.R"
# Run GLM, RF, SVM @ input data subset, return respective AUCs on holdout test

run3ML = function(sample_data, raw_or_proc, normalize, log_transform){
    
    ### Prep data
    
    if(raw_or_proc=="raw"){
        the_data = raw_dt[Sample_id_1 %in% sample_data$Sample_id_1]
    }
    
    if(raw_or_proc=="proc"){
        the_data = proc_dt[Sample_id_1 %in% sample_data$Sample_id_1]
    }
    
    if(raw_or_proc=="both"){
        the_data_r = raw_dt[Sample_id_1 %in% sample_data$Sample_id_1]
        the_data_p = proc_dt[Sample_id_1 %in% sample_data$Sample_id_1]
        the_data = merge(the_data_r, the_data_p, 
                         by=c("Sample_id_1", "isQuiescentGenomeControl"))
    }
    
    the_data[,isControl:=ifelse(isQuiescentGenomeControl==TRUE, "ctrl_0", 
                                "ctrl_1")]
    the_data[,isControl:=as.factor(isControl)]
    the_data[,isQuiescentGenomeControl:=NULL]
    
    the_holdout_ids = sample_data[isHoldout==TRUE, Sample_id_1]
    the_tt_ids = sample_data[isHoldout==FALSE, Sample_id_1]
    
    the_holdout = the_data[Sample_id_1 %in% the_holdout_ids]
    the_tt = the_data[Sample_id_1 %in% the_tt_ids]
   
    the_holdout[,Sample_id_1:=NULL]
    the_tt[,Sample_id_1:=NULL]
    
    if(normalize==TRUE | log_transform==TRUE){
        holdout_case = the_holdout$isControl
        tt_case = the_tt$isControl
        
        the_holdout[,isControl:=NULL]
        the_tt[,isControl:=NULL]
        
        holdout_names = colnames(the_holdout)
        tt_names = colnames(the_tt)
        
        if(log_transform==TRUE){
            the_holdout[,(holdout_names):=lapply(.SD,log_col),
                        .SDcols=holdout_names]
            the_tt[,(tt_names):=lapply(.SD,log_col),
                   .SDcols=tt_names]
        }
        
        if(normalize==TRUE){
            the_holdout[,(holdout_names):=lapply(.SD,norm_col),
                        .SDcols=holdout_names]
            the_tt[,(tt_names):=lapply(.SD,norm_col),
                   .SDcols=tt_names]
        }
  
        the_holdout[,isControl:=holdout_case]
        the_tt[,isControl:=tt_case]
    }
    
    
    #Convert NaN to 0 [this results from subsamples w/ no nonzero for variable]
    for (j in seq_len(ncol(the_holdout))){
        set(the_holdout,which(is.nan(the_holdout[[j]])),j,0)
    }
    
    for (j in seq_len(ncol(the_tt))){
        set(the_tt,which(is.nan(the_tt[[j]])),j,0)
    }
    
    
    
    cols_tt = colnames(the_tt)
    new_cols_tt = str_replace_all(cols_tt, "[^a-zA-Z0-9]", "")
    colnames(the_tt) = new_cols_tt
    
    cols_holdout = colnames(the_holdout)
    new_cols_holdout = str_replace_all(cols_holdout, "[^a-zA-Z0-9]","")
    colnames(the_holdout) = new_cols_holdout
    
    ### GLM
    
    tt_noCase = data.table(the_tt)
    tt_train_control = tt_noCase$isControl
    tt_noCase[,isControl:=NULL]
    
    tt_noCase = as.matrix(tt_noCase)
    
    glm_train = cv.glmnet(x=tt_noCase, y=tt_train_control,
                          nfolds=10, lower.limits=0,
                          type.measure="auc", family="binomial",
                          standardize=FALSE)
    
    glm_fit = glm_train$glmnet.fit
    
    holdout_noCase = data.table(the_holdout)
    holdout_test_control = holdout_noCase$isControl
    holdout_noCase[,isControl:=NULL]    
    
    holdout_noCase = as.matrix(holdout_noCase)    
    
    glm_test = predict(glm_fit, type="response", newx = holdout_noCase)
    
    glm_test = as.data.table(glm_test)
    the_cols = colnames(glm_test)
    col_num = length(the_cols)
    the_pick = the_cols[col_num]
    glm_result = glm_test[,.SD,.SDcols=the_pick]
    glm_result = unlist(glm_result)
    
    roc_obj_glm = roc(predictor=glm_result, response=holdout_test_control)
    
    auc_glm = roc_obj_glm$auc
    
    
    ### SVM
    
    svm_ctrl = trainControl(method="repeatedcv",
                            number=5,
                            repeats=5,
                            classProbs=TRUE,
                            summaryFunction=twoClassSummary,
                            savePredictions=TRUE)
    
    svm_train = train(isControl~., the_tt,
                      method="svmRadial",
                      tuneLength=5,
                      trControl=svm_ctrl,
                      metric="ROC",
                      verbose=FALSE)
    
    svm_probs = predict(svm_train, newdata=the_holdout, type="prob")
    
    svm_pred = svm_probs$ctrl_1
    svm_resp = the_holdout$isControl
    
    roc_obj_svm = roc(predictor=svm_pred, response=svm_resp)
    
    auc_svm = roc_obj_svm$auc
    
    
    ### RF
    
    model_rf = randomForest(isControl~., data=the_tt)
    
    pred_rf = predict(object=model_rf, newdata=the_holdout,
                      type="prob", norm.votes=TRUE)
     
    pred_rf = as.data.table(pred_rf)                   
    for_auc = cbind(pred_rf, the_holdout$isControl)
    for_auc[,isControl:=as.factor(V2)]
    for_auc[,Prob_1:=ctrl_1]
    
    roc_obj_RF = roc(predictor=for_auc$Prob_1, response=for_auc$isControl)
    auc_RF = roc_obj_RF$auc
    
    
    ### Combine (return dt w/ classification probabilities @ holdout)
    
    results = c(auc_glm, auc_svm, auc_RF)
    
    return(results)
}
