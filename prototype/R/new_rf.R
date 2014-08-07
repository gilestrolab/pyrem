rm(list=ls())
DBG = 1
printdbg <- function(){
    print(DBG)
    assign("DBG", DBG+1, envir = .GlobalEnv)
    
}
set.seed(2)

library("randomForest")
library("caret")
library("abind")
library("parallel")


OUT_CSV = "/data/pyrem/Ellys/all_features.csv"
       
exclude_coeffs <-function(name_list, df){
    non_feature_cols <- grep("\\.", colnames(df), invert=T)
    l <- lapply(name_list, grep, colnames(df))
    
    col_idx <- c(non_feature_cols, unique(do.call("c", l)))
    return(df[,col_idx])

}

curate_df <- function(dfo){
    df <- dfo
    # remove cols containing NAs
    s = sapply(df, function(col)sum(is.na(col)))
    df <- df[,s ==0]
    
    
    df$t <- df$X
    df$X <- NULL
    df$animal <- sprintf("%s_%s", df$treatment, df$animal)
    df$animal  <- as.factor(df$animal)
    
    
    # fixme Hjorth activity == mean power
    df  <- df[ ,grep("activity", colnames(df),invert=T)]
    
    # our Y variable is a factor (vigilance state)
    df$y <- as.factor(df[,"vigilance_state.vigil.value"])
    df$vigilance_state.vigil.value <- NULL
    df$py <- df[,"vigilance_state.vigil.proba"]
    df$vigilance_state.vigil.proba <- NULL
    
    # remove ambiguous vigilance states ?
#~     df <- subset( df, vigilance_state.vigil.proba == 1)
#~     df$vigilance_state.vigil.proba <- NULL 
    
    
    return(df)
}

#################################
make_time_features <- function(lag, d, min_max_lag){
    
    # padding with vectro of NAs
    min_lag <-  min_max_lag[1]
    max_lag <- min_max_lag[2]
    
    range <- max_lag -  min_lag
    
    n_na_before <- range - (lag - min_lag) 
    n_na_after <- range - n_na_before
    
    to_concatenate = list(NULL, NULL, NULL)
    
    if(n_na_before > 0)
        to_concatenate[[1]] <- matrix(rep(NA, ncol(d)*n_na_before), nrow=n_na_before, dimnames=list(NULL, colnames(d)))
    if(n_na_after > 0)
        to_concatenate[[3]] <- matrix(rep(NA, ncol(d)*n_na_after), nrow=n_na_after, dimnames=list(NULL, colnames(d)))
    
    to_concatenate[[2]] <- d
    out <- do.call("rbind", to_concatenate)
    colnames(out) <- sprintf("lag_%s.%s", ifelse(lag<0, paste("m",abs(lag),sep=""), as.character(lag)), colnames(out))

    return(out)
    
    }
entropy <- function(v){
    v[v==0] <- 1e-100
    h <- 1 + sum(v * log2(v)) /log2(length(v))
    return(h)
}

add_lagged_time_features <- function(d, pattern_valid_cols = "\\.", min_max_lag=c(-2,2)){
    if (min_max_lag[2] == 0  & min_max_lag[1] ==0) 
        return (d)
    # we ensure data is sorted by time
    d <- d[order(d$t),]
    # we ensure we have an homogenous time series:
    stopifnot(
         length(unique(diff(d$t))) == 1
    )
    
    cols_to_lag_idx <- grep(pattern_valid_cols, colnames(d))
    
    cols_to_lag <- d[,cols_to_lag_idx ]
    
    other_cols <- subset(d,select = -cols_to_lag_idx)
    
    min_lag <-  min_max_lag[1]
    max_lag <- min_max_lag[2]
    
    stopifnot(min_lag < max_lag)
    
    
    other_cols <- other_cols[(abs(min_lag) + 1): (nrow(other_cols) - abs(max_lag)),, drop=F]
    
    
    dfs <- lapply(min_lag : max_lag, make_time_features, d=cols_to_lag, min_max_lag )
    
    all_lagged_cols <- na.omit(do.call("cbind",dfs))
    
    return(cbind(other_cols, all_lagged_cols))
}


########################################################

#~ 



crossval_test <- function(out_level, original_df){
        print(out_level)
        train_df <- subset(original_df, animal !=out_level)
        test_df <- subset(original_df, animal ==out_level)
        train_df$animal <- NULL
        test_df$animal <- NULL
        
        rf <- randomForest(y ~ ., train_df, ntree=40, sampsize=c(200,200,200))        
        
        preds <- predict(rf, test_df, type="prob")
        pred_values <- apply(preds,1,function(v){colnames(preds)[which.max(v)]})
        pred_values <- factor(pred_values, levels=levels(test_df$y))
        #out <- data.frame(real = test_df$y, preds = preds)
        out <- sum(test_df$y == pred_values) / length(pred_values)
        plot(rf)
        strt <- 500
        stp<- 2000
        l = stp - strt +1
        h_rows <- apply(preds, 1, entropy)        

        confidence <- 0:99/100
        CV_accuracy <- sapply(confidence, function(r){
                        valid <- (h_rows > r & h_rows <= r+0.1)
                        sum_posit <- sum((test_df$y == pred_values)[valid])
                        if (sum(valid) == 0)
                            return(NA)
                        return(sum_posit/sum(valid))
                    }
                    
                )
        cm <- t(confusionMatrix (pred_values, test_df$y)$byClass)
        d <- data.frame(cm)
        d <- reshape(d, varying = colnames(d) ,v.names = "y", direction = "long",ids=rownames(d))        
        d$animal <- out_level
        return(d)
        
}

select_variables_pca <- function(df, t=0.99){

#~ 
    pca <- princomp(~. , df[,grep("\\.", colnames(df))])
    
    components <- pca$scores
    cumpca <- cumsum(pca$sdev)
    cumpca <- cumpca - min(cumpca)
    cumpca <- cumpca / max(cumpca)
    components <- components[,1:which(cumpca > t)[1]]
    
    d_out <- df[, grep("\\.", colnames(df), invert=T)]
    d_out <- cbind(data.frame(components), d_out)
    d <- subset( d_out, py == 1)
    d$py <- NULL 
    d$treatment <- NULL 
    d$t <- NULL 
    d$animal <- NULL 
    d$y <- droplevels(d$y)
    rf <- randomForest(y ~ ., d, ntree=100, sampsize=c(1000,1000,1000))
    print(rf)
#~     #d <- subset( d, py == 1)
    return(d_out)
#~     
    #remove ammbiguous vigilance states    
    
}
select_variables <- function(df, t=0.50){
    
    d <- subset( df, py == 1)
    d$py <- NULL 
    d$treatment <- NULL 
    d$t <- NULL 
    d$animal <- NULL 
    d$y <- droplevels(d$y)
    
#~     result <- rfcv(subset(d, select=-y),d$y, cv.fold=10, ntree=50, sampsize=c(1000,1000,1000))
    
    
#~     with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))
#~     return(out)
    
    rf <- randomForest(y ~ ., d, ntree=50, sampsize=c(1000,1000,1000))
#~     return(rf)
    print(rf)
    print(importance(rf))
    imp_df <- data.frame(importance(rf))
    bad_vars <- subset(imp_df, MeanDecreaseGini < quantile(rf$importance,t))
    match_bad_name = match(rownames(bad_vars), colnames(df))
    return(subset(df, select = -match_bad_name))
#~ 
#~     
    
    }
important_subbands <- function(df){
    d <- subset( df, py == 1)
    d$py <- NULL 
    d$treatment <- NULL 
    d$t <- NULL 
    d$animal <- NULL 
    d$y <- droplevels(d$y)
    
    rf <- randomForest(y ~ ., d, ntree=50, sampsize=c(1000,1000,1000))
    imp_df <- data.frame(do.call("rbind",strsplit(rownames(rf$importance), "\\.")))
    imp_df$y <- rf$importance[,1]

    
    }

test_different_lags <- function(tau, df){
#~     #                0,1,2,3,4,5,6,7,8,9,0
#~     clustersize <- c(4,4,4,3,3,3,2,2,2,2,2) - 1
    
    print("==========================")
    print(paste("lag =", tau))
    anim_df = split(df,df$animal)
    lagged_anim_dfs <- lapply(anim_df, add_lagged_time_features,  min_max_lag=c(-tau,+tau))
    df <- do.call("rbind",lagged_anim_dfs)

    #remove ambiguous vigilance states
    df <- subset( df, py == 1)
    df$py <- NULL 
    df$y <- droplevels(df$y)
    df$treatment <- NULL 
    
    print(dim(df))
#~     l = sapply(levels(df$animal), crossval_test, original_df=df)
#~     print(clustersize[tau+1])
#~     cl <- makeCluster(clustersize[tau+1])
#~     clusterExport(cl, "randomForest")
#~     clusterExport(cl, "entropy")
#~     clusterExport(cl, "confusionMatrix")
#~     l = parLapply(cl, levels(df$animal), crossval_test, original_df=df)
    l = lapply(levels(df$animal), crossval_test, original_df=df)
#~     stopCluster(cl)
    d <- do.call("rbind",l)
    d$lag <- tau
    toprint <-subset(d,id == "Sensitivity" |id == "Specificity" |id == "Pos Pred Value"   )
    
    print (aggregate(y ~ id * time , toprint, mean))
    return(d)
}

"
Starts here
"


dfo <- read.csv(OUT_CSV , na.string="NaN")
#dfo <- readRDS(file="/data/pyrem/Ellys/all_features.rds")

df <- curate_df(dfo)


print("Selecting variables")
#~ df <- select_variables(df)


#~ good_coeffs <-c(
#~ 
#~                 "EMG_1_cD_1",
#~                 "EMG_1_cD_2",
#~                 "EMG_1_cD_3",
#~                 "EEG_parietal_cereb_cD_3",
#~                 "EEG_parietal_cereb_cD_4",
#~                 "EEG_parietal_cereb_cD_5",
#~                 "EEG_parietal_cereb_cD_6"
#~                 )
#~ 
#~ good_coeffs <-c("EMG_1_cD_4",
#~                 "EMG_1_cD_5",
#~                 "EMG_1_cD_6")
#~          


#df <- exclude_coeffs(good_coeffs,df)
#~ 
#~ rf <- select_variables(df)
#~ df <- select_variables_pca(df)
#df <- select_variables(df)

l_dfs <- lapply(rep(0:3, 3), test_different_lags, df)

stop()
############################################################################
#~
#df <- subset(df, df$vigil_value==87)
#df$vigil_value <- NULL

#~
#~ crossval <- function(out_level, original_df){
#~         train_df <- subset(original_df, animal !=out_level)
#~         #train_df <- original_df
#~         # train_df <- original_df
#~         test_df <- subset(original_df, animal ==out_level)
#~         print(out_level)
#~
#~          # print(head(test_df))
#~          # print(unique(train_df$treatment))
#~          # print(unique(test_df$treatment))
#~
#~          train_df$animal <- NULL
#~          test_df$animal <- NULL
#~          rf <- randomForest(treatment ~ ., train_df, ntree=50)
#~          varImpPlot(rf)
#~
#~          votes <- predict(rf, test_df,type="prob")
#~
#~          votes <- matrix(votes, nrow=length(rownames(votes)),dimnames=list(NULL,colnames(votes)))
#~
#~          h_rows <- apply(votes, 1, entropy)
#~          weighted_votes <- colSums(votes * h_rows )
#~
#~
#~          #
#~          # weighted_votes = weighted_votes /sum(weighted_votes)
#~          #
#~          # # tab <- matrix(weighted_votes, nrow=1,dimnames=list(NULL,colnames(votes)))
#~          p = max(weighted_votes) / sum(weighted_votes)
#~          out <- data.frame(pred = colnames(votes)[which.max(weighted_votes)], real=unique(test_df$treatment), p=p)
#~
#~          print(weighted_votes)
#~          return(out)
#~ }

#~
#~ entropy <- function(v){
#~     v[v==0] <- 1e-100
#~     h <- 1 + sum(v * log2(v)) /log2(length(v))
#~     return(h)
#~ }
#~ set.seed(1)
#~ l = lapply(levels(df$animal), crossval, original_df=df)
#~
#~ d = do.call("rbind",l)
#~ print(d)

# predictions = do.call("rbind", l)
# predictions <- within(predictions, pred <- ifelse(GFP > TelC, "GFP", "TelC"))

#
######################################################################################################


#~
#~

#
#
# dpdf("/tmp/test.pdf")
# for (a in levels(df$animal)){
#   sub_df <- subset(df, animal == a)
#   plot( log10(power_median) ~ t,sub_df, pch=20, col=treatment,ylim=c(6,10),main=a)
#   }
# dev.off()




#############
#~ set.seed(1)
#~ out_level <- "GFP_D"
#~ train_df <- subset(df, animal !=out_level)
#~ test_df <- subset(df, animal ==out_level)
#~
#~ train_df$animal <- NULL
#~ test_df$animal <- NULL
#~
#~ ytrain <- train_df$treatment
#~ ytest <- test_df$treatment
#~
#~ xtrain<- train_df
#~ xtrain$treatment <- NULL
#~ xtest <- test_df
#~ xtest$treatment <- NULL
#~
#~ rf <- randomForest(x = xtrain, y= ytrain, ntree=100)
#~
#~ table(predict(rf, xtest))


'
Class: 78  Class: 82 Class: 87
0.9094319 0.84323282 0.8801143  # "coeffs eeg and emg, D_2:6 (no A)"
0.8943559 0.82718649 0.8678393 # "coeffs eeg D_2:6 , emg_5_6 (no A)"
0.9014322 0.84034602 0.8681586 # "coeffs eMG D_2:6 , eEG_5_6 (no A)"
0.8803510 0.83611262 0.8641646 # "coeffs eeg and emg, D_4:6 (no A)"
0.9141507 0.81931153 0.8861858 # all data
0.8872067 0.84670554 0.8764197 #  "coeffs eeg and emg, D_3:6 (no A)"
0.8777605 0.82924432 0.8722724 # "coeffs eMG D_1:2 , eEG_4_6 (no A)"

'
