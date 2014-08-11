


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
graphics.off()
set.seed(12134)
library(markovchain)
library(randomForest)
library("abind")
library("parallel")
######################################################
OUT_PREFIX<- "/data/pyrem/Ellys/all_features"
VARIABLE_TO_USE <-c(
	"EEG_parietal_frontal_cD_1.power.mean",        
	"EEG_parietal_frontal_cD_1.power.min",        
	
	"EEG_parietal_frontal_cD_6.power.mean",        
	"EEG_parietal_frontal_cD_6.power.min",        
	
	"EEG_parietal_frontal_cD_5.hjorth.complexity",   
	"EEG_parietal_frontal_cD_5.hjorth.morbidity",
	
	"EEG_parietal_frontal_cD_5.power.mean",
	"EEG_parietal_frontal_cD_5.power.min",
	
	"EMG_REF_cD_1.power.min",                                          
	"EMG_REF_cD_1.power.mean",                                          
	"EMG_REF_cD_2.power.min",                       
	"EMG_REF_cD_2.power.mean",                       
	"EMG_REF_cD_3.power.min",
	"EMG_REF_cD_3.power.mean"
)

entropy_conf <- function(v){
    v[v==0] <- 1e-100
    h <- 1 + sum(v * log2(v)) /log2(length(v))
    return(h)
}

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
	
stratified_cv_one_level <- function(out_level, df, rf_fun, return_preds = FALSE){
	test_df <- subset(df, animal ==out_level)
	train_df <- subset(df, animal !=out_level)
	test_df <- strip_df_for_ml(test_df)
	sample_size <- 2000
	
	test_classes <- split( test_df, test_df$y)
	test_classes_l <- lapply(test_classes, function(d){
		idxs <- sample(1: nrow(d), sample_size, replace=T)
		return(d[idxs,])
		})

	test_df <- do.call("rbind", test_classes_l)
	rf <- rf_fun(train_df)

	
	if(return_preds){
		preds_p <- predict(rf, test_df, type="prob")
		conf <- apply(preds_p, 1, entropy_conf)
		preds <- apply(preds_p, 1, which.max)
		preds <- levels(test_df$y)[preds]
		out_df <- data.frame(y= test_df$y, preds = preds, conf, animal=out_level)
		return(out_df)
	}
	
	preds <- predict(rf, test_df)
	error <- sum(test_df$y != preds) / length(preds)
	
	
	return(list(error=error, rf=rf))
}



stratified_cv <- function(df, rf_fun, return_preds=FALSE, nprocesses=4){
	
	if(nprocesses < 2){
#~ 		if(! return_rfs){
		errors <- sapply(levels(df$animal), function(x)stratified_cv_one_level(x, df, rf_fun)$error)
#~ 		}
#~ 		else{
#~ 			errors_rfs <- lapply(levels(df$animal), stratified_cv_one_level,df, rf_fun)
#~ 			mean_error <- mean(sapply(errors_rfs,function(l)l$error))
#~ 			rfs <- lapply(errors_rfs,function(l)l$rf)
#~ 		}
	}
	else{
		cl <- makeCluster(nprocesses)
		clusterExport(cl, "randomForest")
		clusterExport(cl, "stratified_cv_one_level")
		clusterExport(cl, "strip_df_for_ml")
		clusterExport(cl, "entropy_conf")
		clusterExport(cl, "rf_fun")
		if(! return_preds)
			out <- parSapply(cl, levels(df$animal), function(x)stratified_cv_one_level(x, df, rf_fun)$error)
		else{
			dfs <- parLapply(cl, levels(df$animal), function(x)stratified_cv_one_level(x, df, rf_fun,return_preds))
			out <- do.call("rbind", dfs)
		}
			
		stopCluster(cl)
		
	}
	
	return(out)
}
select_variables <- function(df, scale=1.2){
	
	# importantly!! importance has to be computed on a balanced sample
    rf_fun <- function(dd){
				d <- strip_df_for_ml(dd);
				return(randomForest(y ~ ., d, ntree=50, sampsize=c(1000,1000,1000)))
				}

    cv_error_vec <- stratified_cv(df, rf_fun)
    
    
    rf <- rf_fun(df)

	row_names <- gsub("sample_2_0\\.", "sample_2_0d",rownames(rf$importance))
	t <- quantile(rf$importance, 1-1/scale)
	
    var_to_remove <- names(rf$importance[rf$importance < t,])
	bad_idxs = match(var_to_remove, colnames(df))
	if(length(bad_idxs) == 0)
		return(NULL)
    out <- df[,-bad_idxs]
    
    # hacky way to get the real number of variable used
    d <- strip_df_for_ml(df)
    
    return(list(data=out, cv_error_vec=cv_error_vec, n_vars= ncol(d) - 1))
}

recursive_variable_elimination <- function (df, scale=1.2,target=2){
	
	variables = list()
	i=1
	while(TRUE){
		
		selection <- select_variables(df, scale=1.3)
		if (is.null(selection)){
			break
		}
		selection$select_variables
		n_vars <- selection$n_vars
		out_tmp <- data.frame(error = selection$cv_error_vec, 
							n_vars = n_vars,
							animal=names(selection$cv_error_vec))
		
		if(exists("out_df"))
			out_df <- rbind(out_df, out_tmp)
		else
			out_df <- out_tmp
		
		cn <- colnames(df)
		df <- selection$data
		variables[[i]] <- cn
		
		i <- i + 1
		
		if(min(n_vars) <= target){
			break
		}
		print(out_df)
	}
	return(list(data = out_df, variables))
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

	for (i in 1:replicates){
		for (tau in 0:max_tau){
			dd <- add_lagged_time_features_to_whole(df, min_max_lag=c(-tau, +tau))
			
		# importantly!! importance has to be computed on a balanced sample
			rf_fun <- function(dd){
				d <- strip_df_for_ml(dd);
				return(randomForest(y ~ ., d, ntree=50, sampsize=c(1000,1000,1000)))
				}
			
			errors <- stratified_cv(dd, rf_fun)
			if(exists("error_lag_df"))
				error_lag_df <- rbind(error_lag_df, data.frame(error=errors, lag=tau))
			else
				error_lag_df <- data.frame(error=errors, lag=tau)
			
		}

	}
	 
	boxplot(error ~ lag, error_lag_df, pch=20, xlab="tau", ylab=c("Average cross validation class error"),ylab=c(0,0.15))
	return(error_lag_df)
}
	
analize_error_vs_nvar <- function(df, replicates=3){
	l <- lapply(1:replicates, function(dummy, df)recursive_variable_elimination(df)$data, df)
	error_nvar_df <- do.call("rbind", l)
	boxplot(error ~ n_vars, error_nvar_df, pch=20, xlab=expression(italic(p)), ylab=c("Average cross validation class error"),ylab=c(0,0.15))
	return(error_nvar_df)
	}
	
# makes a  new df where y is predicted y. prediction is realised through CV

analize_error_vs_confidence <- function(df){
	rf_fun <- function(dd){
		d <- strip_df_for_ml(dd);
		# 200 trees for accuracy of confidence
		return(randomForest(y ~ ., d, ntree=200, sampsize=c(1000,1000,1000)))
		}
	dd <- stratified_cv(df, rf_fun, TRUE)
	
	conf_bins <- 0:9/10
	error_vs_conf <- sapply(0:9/10, function(c, d){
		sub_d <- subset(d, conf>c &  conf <= c+0.1)
		errors <- sub_d$y != sub_d$preds
		return(sum(errors) / nrow(sub_d))
		}, dd)
	barplot(error_vs_conf, width=1,space=0, names.arg=conf_bins, ylim=c(0,0.5), ylab= "Actual error", xlab="Prediction confidence")
	hist(dd$conf, nclass=20, col="grey", freq=F, xlab="Prediction confidence", main="")
}
make_cv_predictions <- function(df){
	
	}

REPLICATES <- 1
#~ main <- function(){
	pdf("/tmp/todel.pdf")
	df <- cache_load_file(OUT_PREFIX)
	dfo <- curate_df(df)
	
#~ 	describe_time_series(dfo)
#~ 	error_nvar_df <- analize_error_vs_nvar(dfo,replicates=REPLICATES)
#~  	
#~ 	out <- recursive_variable_elimination(dfo, target=25)
#~ 	
#~ 	
	non_predictor_cols <- grep("\\.", colnames(dfo), invert=T, value=TRUE)
	df <- dfo[,c(VARIABLE_TO_USE, non_predictor_cols)]
#~     error_lag_df <- analize_error_vs_lag(df, 2, replicates=REPLICATES)
#~     

	analize_error_vs_confidence(df)

	#todo:
#~ 	CV transition matrices and prevalence...

#~ 	dev.off()

###finalrf:
#~ 	dd <- add_lagged_time_features_to_whole(df, min_max_lag=c(-2, +2))
#~ 	d <- strip_df_for_ml(dd)
#~ 	# importantly!! importance has to be computed on a balanced sample
#~ 	rf <- randomForest(y ~ ., d, ntree=50)
#~ 	
#~ }




#sleep stages? within non-rem?
#~ plot3d(x=log10(df$EMG_REF_cD_1.power.mean), y=log10(df$EEG_parietal_frontal_cD_6.power.mean), z=log10(df$EEG_parietal_frontal_cD_1.power.mean), col=df$y)
