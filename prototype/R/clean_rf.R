


#~ 
#~ In order to assess performance of our predictors, we compared cross validation misscalssification rate.
#~ two important precaussions were taken to make this assesment fair:
#~ 1) Since features are actually time series, we can suspect an significant statistical dependence between sucessive states.
#~  Conventionnal cross validation would test a predictor trained with a random sample of the whole data against another random part of the data.
#~  This would provide a fair assessment of a classifier and allow to determine if it is overfit.
#~  In the case of temporal, or spatial data, a "blindly" random sample would not account for the temporal structure.
#~  Such a procedure would fail to detect an overfit model because most of the data that were removed is very corelated to some of the data remaining in the training set.
#~  For this reason, a stratified crossvalidation procedure was performed to asses the classifiers; one entiere time as used as for cross validation while all the other time series constituted the training set.
#~ 2) The different sleep stages have unequel prevalence, so, testting missclassification on en entiere time series would not account for misscalssification of minority class.
#~  for instance, an error rate of only 10% would be observed if all REM sleep epoch were misclassified.



######################################################
rm(list = ls())
set.seed(12134)
library(markovchain)
library(randomForest)
library("abind")
######################################################
OUT_PREFIX<- "/data/pyrem/Ellys/all_features"



curate_df <- function(df){
    
    # remove cols containing NAs
    s = sapply(df, function(col)sum(is.na(col)))
    if(sum(s) != 0)
		print(paste("removing columns: ", colnames(df)[which(s >0)]))
	
    df <- df[,s ==0]
    
    
    df$t <- df$X
    df$X <- NULL
    df$animal <- sprintf("%s_%s", df$treatment, df$animal)
    df$animal  <- as.factor(df$animal)

    # our Y variable is a factor (vigilance state)
    df$y <- as.factor(df[,"vigilance_state.vigil.value"])
    df$vigilance_state.vigil.value <- NULL
    df$py <- df[,"vigilance_state.vigil.proba"]
    df$vigilance_state.vigil.proba <- NULL
       
    return(df)
}




cache_load_file <- function(prefix){
	out_csv = paste(prefix, ".csv",sep="")
	out_rds= paste(prefix, ".rds",sep="")
	if(file.exists(out_rds)){
		return(readRDS(file=out_rds))
	}
	else{
		df <- read.csv(out_csv, na.string="NaN")
		saveRDS(df, file=out_rds)
		return(df)
	}
}


to_homogenous_univariate_ts <- function(df){
	

	df <- df[order(df$t),]
	# we ensure we have an homogenous time series:
	stopifnot(
		 length(unique(diff(df$t))) == 1
	)
	df <- subset(df, py ==1)
	return(na.omit(df)$y)
}

#######
#' A simple description of the time series in terms of prevalence and
describe_time_series <- function(df){
	anim_dfs = split(df,df$animal)
	time_series <- lapply(anim_dfs, to_homogenous_univariate_ts)
	time_series <- lapply(time_series, droplevels)
	prevalences <- sapply(time_series, function(x){out <- table(x); return(out/sum(out))})
	
	print("Prevalence means:")
	print(round(apply(prevalences, 1, mean),3))
	print("Prevalence sds:")
	print(round(apply(prevalences, 1, sd),3))
	transit_mats <-lapply(time_series, function(x)markovchainFit(data=x)$estimate[,])
	trans_array <- abind(transit_mats, along=3)
	print("Empirical average transition matrix:")
	print(round(apply(trans_array, c(1,2), mean),3))
	print("SDs:")
	print(round(apply(trans_array, c(1,2), sd),3))
	
}

strip_df_for_ml <- function(df){
	d <- subset( df, py == 1)
	d$py <- NULL 
	d$treatment <- NULL 
	d$t <- NULL 
	d$animal <- NULL 
	d$y <- droplevels(d$y)
	return(d)
	}
	
stratified_cv_one_level <- function(out_level, df, rf){
	test_df <- subset(df, animal ==out_level)
	sample_size <- min(table(test_df$y))
	test_classes <- split( test_df, test_df$y)
	test_classes_l <- lapply(test_classes, function(d){
		idxs <- sample(1: nrow(d))
		return(d[idxs,])
		})
	test_df <- do.call("rbind", test_classes_l)
	test_df <- strip_df_for_ml(test_df)
	
	preds <- predict(rf, test_df)
	error <- sum(test_df$y != preds) / length(preds)
	return(error)
}

stratified_cv <- function(df, rf, nprocesses=1){
	if(nprocesses < 2){
		return(mean(sapply(levels(df$animal), stratified_cv_one_level, df, rf)))
	}
}
select_variables <- function(df, scale=1.2){
	d <- strip_df_for_ml(df)
	# importantly!! importance has to be computed on a balanced sample
    rf <- randomForest(y ~ ., d, ntree=50, sampsize=c(1000,1000,1000))
    row_names <- gsub("sample_2_0\\.", "sample_2_0d",rownames(rf$importance))
    cv_error <- stratified_cv(df, rf)
    imp_df <- data.frame(do.call("rbind",strsplit(row_names, "\\.")))
    imp_df$y <- rf$importance[,1]
    imp_df$var_name <- rownames(rf$importance)
	t <- quantile(rf$importance, 1-1/scale)
	
    var_to_remove <- names(rf$importance[rf$importance < t,])
	bad_idxs = match(var_to_remove, colnames(df))
	if(length(bad_idxs) == 0)
		return(NULL)
    out <- df[,-bad_idxs]
    return(list(data=out, cv_error=cv_error, n_vars= ncol(d) - 1))

}

recursive_variable_elimination <- function (df, scale=1.2,target=2){
	

	error = numeric()
	n_vars = numeric()
	variables = list()
	i=1
	while(TRUE){
		selection <- select_variables(df, scale=1.3)
		if (is.null(selection)){
			break
		}
		error <- c(error, selection$cv_error)
		n_vars <- c(n_vars, selection$n_vars)
		cn <- colnames(df)
		df <- selection$data
		variables[[i]] <- cn
		i <- i + 1
		print(i)
		if(min(n_vars) <= target){
			break
		}
	}
	return(list(data = data.frame(n_vars, error), variables))
}
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
add_lagged_time_features <- function(d, min_max_lag=c(-2,2),pattern_valid_cols = "\\."){
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

add_lagged_time_features_to_whole <- function(df,min_max_lag=c(-2,2),pattern_valid_cols = "\\."){
	l_dfs <- split(df, df$animal)
	l_dfs <- lapply(l_dfs,add_lagged_time_features, min_max_lag=min_max_lag,pattern_valid_cols=pattern_valid_cols)
	return(do.call("rbind", l_dfs))
	}

analize_error_vs_lag <- function(df, max_tau, replicates=3){
	taus = numeric()
	errors = numeric()
	for (i in 1:replicates){
		for (tau in 0:max_tau){
			dd <- add_lagged_time_features_to_whole(df, min_max_lag=c(-tau, +tau))
			d <- strip_df_for_ml(dd)
		# importantly!! importance has to be computed on a balanced sample
			rf <- randomForest(y ~ ., d, ntree=50, sampsize=c(1000,1000,1000))
			print(ncol(d) - 1)
			taus <- c(taus,tau)
			errors <- c(errors, stratified_cv(dd, rf))
		}

	}
	error_lag_df <- data.frame(error=errors, lag=taus)
	plot(error ~ lag, error_lag_df, pch=20)
	return(error_lag_df)
}
	
analize_error_vs_nvar <- function(df, replicates=3){
	l <- lapply(1:replicates, function(dummy, df)recursive_variable_elimination(df)$data, df)
	error_nvar_df <- do.call("rbind", l)
	plot(error ~ n_vars, error_nvar_df, pch=20,log="x")
	return(error_nvar_df)
	}
	

#~ main <- function(){
	df <- cache_load_file(OUT_PREFIX)
	dfo <- curate_df(df)
	
	describe_time_series(dfo)
	error_nvar_df <- analize_error_vs_nvar(dfo,replicates=10)
 	
	out <- recursive_variable_elimination(dfo, target=25)
	
    df <- subset(dfo, select = out[[2]][[length(out[[2]])]])

    error_lag_df <- analize_error_vs_lag(df, 4, replicates=10)

	#todo:
#~ 	cv_transition matrices and prevalence...
###finalrf:
#~ 	dd <- add_lagged_time_features_to_whole(df, min_max_lag=c(-2, +2))
#~ 	d <- strip_df_for_ml(dd)
#~ 	# importantly!! importance has to be computed on a balanced sample
#~ 	rf <- randomForest(y ~ ., d, ntree=50)
#~ 	
#~ }

