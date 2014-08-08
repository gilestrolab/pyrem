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
library("MASS")


OUT_CSV = "/data/pyrem/Ellys/all_features.csv"
       
keep_only_grep <-function(name_list, df){
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
        
        #rf <- randomForest(y ~ ., train_df, ntree=40, sampsize=c(200,200,200))        
        rf <- randomForest(y ~ ., train_df, ntree=100, sampsize=c(5000,1000,5000))        
        
        
        preds <- predict(rf, test_df, type="prob")
        pred_values <- apply(preds,1,function(v){colnames(preds)[which.max(v)]})
        pred_values <- factor(pred_values, levels=levels(test_df$y))
        #out <- data.frame(real = test_df$y, preds = preds)
        out <- sum(test_df$y == pred_values) / length(pred_values)
        
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
        plot(rf, ylim=c(0,1),lwd=2)
        print(rf)
        d <- data.frame(cm)
        d <- reshape(d, varying = colnames(d) ,v.names = "y", direction = "long",ids=rownames(d))        
        d$animal <- out_level
        return(d)
#~         
}
#~ 
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
    
    rf <- randomForest(y ~ ., d, ntree=100, sampsize=c(1000,1000,1000))
    print(rf)
    plot(rf)
    imp_df <- data.frame(importance(rf))
    bad_vars <- subset(imp_df, MeanDecreaseGini < mean(rf$importance) * t)
    print(paste("excluding", nrow(bad_vars), "variables"))
    match_bad_name = match(rownames(bad_vars), colnames(df))
    return(subset(df, select = -match_bad_name))
#~ 
#~     
    
    }

make_definitive_rf <- function (df){
    anim_df = split(df,df$animal)
    lagged_anim_dfs <- lapply(anim_df, add_lagged_time_features,  min_max_lag=c(-4,+4))
    df <- do.call("rbind",lagged_anim_dfs)
    d <- subset( df, py == 1)
    d$py <- NULL 
    d$treatment <- NULL 
    d$t <- NULL 
    d$animal <- NULL 
    d$y <- droplevels(d$y)

    rf <- randomForest(y ~ ., d, ntree=100, sampsize=c(5000,2000,5000))
    return(rf)
}
important_subbands <- function(df){
    rf <- make_definitive_rf(df)
    
    print(rf)
    plot(rf)
    imp_df <- data.frame(do.call("rbind",strsplit(rownames(rf$importance), "\\.")))
    imp_df$y <- rf$importance[,1]
    return(imp_df)
    
    }
#ddd <- important_subbands(df)
#~  
 
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
#exclude animal:
df <- subset(df, animal != "GFP_3")
df$animal <- droplevels(df$animal)
df <- keep_only_grep(c("power"),df)
df <- keep_only_grep(c("mean"),df)

#~ df <- keep_only_grep(good_coeffs,df)
#~ 
#~ rf <- select_variables(df)
#~ df <- select_variables_pca(df)
#~ df <- select_variables(df)

l_dfs <- lapply(rep(0:3, 1), test_different_lags, df)
cv_df <- do.call("rbind", l_dfs)
boxplot(y ~ lag * time, subset(cv_df, id == "Sensitivity"), log="y")
boxplot(y ~ lag * time, subset(cv_df, id == "Specificity"), log="y")
boxplot(y ~ lag * time, subset(cv_df, id == "Pos Pred Value"), log="y")

stop()


show_2d_hist <- function(x, z, y, levels= c(.25,0.75,0.5)){
    
    mat <- cbind(x, z)
    contour(kde2d(x, z,n=100,  h=c(0.15,.15)), levels=c(0), main= "Distribution of actual labels")
    
    sub_mats <- split(data.frame(mat), y)
    
    kerns = lapply(sub_mats,function(m)kde2d(m[,"x"], m[,"z"]))
    cols =as.numeric(as.factor(names(sub_mats)))
    
    for (i in 1:length(kerns))
        contour(kerns[[i]],levels=levels, add=T, col=cols[i],n=100,  h=c(0.3,.3))
    legend("bottomleft", legend=names(sub_mats), lwd=2,col=cols,title="78 -> N; 83 -> R; 88->W")
    
    #sub sample d to keep 2000 of each class:
    l <- split(data.frame(x,z), y)
    l <- lapply(l, function(xx){
        idx = 1:nrow(xx)
        return(xx[sample(idx, 5000, replace=T),])
        })
    dsub <- do.call("rbind", l)
    kern <- kde2d(dsub[,"x"],dsub[,"z"],n=100,  h=c(0.3,0.3))
    contour(kern, levels=0:9/10,main="Distribution of a balanced sample")
#~     
#~     if (! is.null(rf))
#~     contour(kde2d(x, z,n=100,  h=c(0.15,.15)), levels=c(0))
#~         new_y <- predict(rf, dsub)
#~         sub_mats <- split(dsub, new_y)
#~         
#~         kerns = lapply(sub_mats,function(m)kde2d(df$lag_0.EEG_parietal_frontal_cD_6.power.mean / df$lag_0.EEG_parietal_frontal_cD_1.power.mean, m$lag_0.EMG_REF_cD_3.power.mean))
#~         cols =as.numeric(as.factor(names(sub_mats)))
#~         
#~         for (i in 1:length(kerns))
#~             contour(kerns[[i]],levels=levels, add=T, col=cols[i],n=100,  h=c(0.3,.3))
#~         legend("bottomleft", legend=names(sub_mats), lwd=2,col=cols)
    }
#show_2d_hist(log10(d$EEG_parietal_frontal_cD_5.power.median), log10(d$EMG_REF_cD_3.power.median), d$y)

df <- curate_df(dfo)
#exclude animal:
df <- subset(df, animal != "GFP_3")
df$animal <- droplevels(df$animal)
df <- keep_only_grep(c("power"),df)
df <- keep_only_grep(c("mean"),df)
anim_df = split(df,df$animal)
lagged_anim_dfs <- lapply(anim_df, add_lagged_time_features,  min_max_lag=c(-4,+4))
df <- do.call("rbind",lagged_anim_dfs)
d <- subset( df, py == 1)
d$py <- NULL 
d$treatment <- NULL 
d$t <- NULL 
d$animal <- NULL 
d$y <- droplevels(d$y)

pdf("/tmp/todel.pdf", h=9, w=16)
par(mfrow=c(1,3))

#~ rf <- make_definitive_rf(df)

dsub <- d
dsub$y <- predict(rf, dsub)

l <- split(dsub, new_y)
l <- lapply(l, function(xx){
    idx = 1:nrow(xx)
    return(xx[sample(idx, 5000, replace=T),])
    })
dsub <- do.call("rbind", l)


plot(log10(dsub$lag_0.EEG_parietal_frontal_cD_6.power.mean / dsub$lag_0.EEG_parietal_frontal_cD_1.power.mean), log10(dsub$lag_0.EMG_REF_cD_3.power.mean), col=dsub$y,
 pch=20,
 xlab="cD_6 / cD_1 (~ Theta:Delta)",
 ylab="EMG's c_D3 mean power",
 main="Model predictions")
 
 show_2d_hist(log10(d$lag_0.EEG_parietal_frontal_cD_6.power.mean / d$lag_0.EEG_parietal_frontal_cD_1.power.mean), log10(d$lag_0.EMG_REF_cD_3.power.mean), d$y, levels= 2:9 /10)
dev.off()
