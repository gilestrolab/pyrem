


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
library(lmerTest)
library("abind")
library("parallel")
library("betareg")
library("betareg")
library(MASS)

source("./pair_of_classes_importance.R")
######################################################


show_2d_hist <- function(x, z, y, levels= c(.25,0.75,0.5),...){
		mat <- cbind(x, z)
		contour(kde2d(x, z,n=50,  h=c(0.15,.15)), levels=c(10), ...)
		sub_mats <- split(data.frame(mat), y)
		kerns = lapply(sub_mats,function(m)kde2d(m[,"x"], m[,"z"]))
		cols =as.numeric(as.factor(names(sub_mats)))
		for (i in 1:length(kerns))
			contour(kerns[[i]],levels=levels, add=T, col=cols[i],n=100,  h=c(0.3,.3))
		legend("bottomleft", legend=names(sub_mats), lwd=2,col=cols,title="78 -> N; 83 -> R; 88->W")
		
		
#~ 		errors <- preds != y
#~ 		hcols <- heat.colors(10)
#~  		points(x[errors], z[errors], col=hcols[round(conf[errors]*10)], pch=".",cex=3) 
#~  		points(x[errors], z[errors], col=hcols[round(conf[errors]*10)], pch=".",cex=3) 
#~ 		contour(kde2d(x[errors], z[errors],n=50,  h=c(0.15,.15)), levels=c(0:4/4), add=T)#. col="grey")
#~  		hist(conf[errors])
 		
	}
#~ show_2d_hist(log10(dd$win_7.EEG_parietal_frontal_cD_6.power.mean), log10(dd$win_7.EMG_REF_cD_1.power.mean), preds$true_y, preds$y, preds$conf)




make_length_per_class <- function(series){
	
	diff <-  c(FALSE, series[1:length(series)-1] != series[-1])
	trans_idx <- which(diff) 
	starts <- c(1, trans_idx)
	ends <-  c(trans_idx-1, length(series))
	start_end <- cbind(starts, ends)
	length_per_class <- apply(start_end, 1, function(x){
		sub_series <- series[x[1]: x[2]]
		class <- unique(sub_series)
		stopifnot(length(class) == 1)
		data.frame(class = class, length=length(sub_series))
		}
		)
	length_per_class <- do.call("rbind",length_per_class)
	
	}
epoch_lengths_and_number <- function(series){
	length_per_class <- make_length_per_class(series)
	mean_length_per_class <- aggregate(length ~ class, length_per_class,mean)
	n_epoch_per_class <- table(length_per_class$class)
	return(list(length=mean_length_per_class, number=n_epoch_per_class))
	}
	
entropy_conf <- function(v){
    v[v==0] <- 1e-100
    h <- 1 + sum(v * log2(v)) /log2(length(v))
    return(h)
}

normalise_cols <- function(d, pattern_valid_cols = "\\.power.m"){
			d <- d[order(d$t),]
			cols_to_norm_idx <- grep(pattern_valid_cols, colnames(d))

			cols_to_norm <- d[,cols_to_norm_idx]
			
			cols_to_norm <- log10(cols_to_norm)

			cols_to_norm <- apply(cols_to_norm, 2, function(x){
				return((x - mean(x)) / sd(x))
				})
			other_cols <- subset(d,select = -cols_to_norm_idx)
			
			out <- cbind(cols_to_norm,other_cols)
				
			}
			
curate_df <- function(df,log_norm_all_Xs = T,pattern_valid_cols = "\\.power.m"){
    
    # remove cols containing NAs
    s = sapply(df, function(col){
		if( sum(is.na(col)) + sum(!is.finite(col)) > 0)
			return(T)
		return(F)
		})
    
    if(sum(s) != 0)
		print(paste("removing columns: ", colnames(df)[s]))
	
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
       
    if(log_norm_all_Xs){
		l_dfs <- split(df, df$animal)
		
		l_dfs <- lapply(l_dfs,normalise_cols, pattern_valid_cols)
		df <- do.call("rbind", l_dfs)
		}
    return(df)
}

cache_load_file <- function(prefix){
	out_csv = paste(prefix, ".csv",sep="")
	out_rds= paste(prefix, ".rds",sep="")
	if(file.exists(out_rds)){
		out <- readRDS(file=out_rds)
		
		return(out)
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

strip_df_for_ml <- function(df, keep_animal=F){
	d <- subset( df, py == 1)
	d$py <- NULL 
	d$treatment <- NULL 
	d$t <- NULL 
	
	if(!keep_animal)
		d$animal <- NULL 
		
	d$y <- droplevels(d$y)
	return(d)
	}
	
stratified_cv_one_level <- function(out_level, df, rf_fun,return_preds=FALSE, resample_test=TRUE, sample_size=750){
	test_df <- subset(df, animal ==out_level)
	train_df <- subset(df, animal !=out_level)
	test_df <- strip_df_for_ml(test_df)
	
	if(resample_test){
		test_classes <- split( test_df, test_df$y)
		test_classes_l <- lapply(test_classes, function(d){
			idxs <- sample(1: nrow(d), sample_size, replace=T)
			return(d[idxs,])
			})

		test_df <- do.call("rbind", test_classes_l)
	}
	
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

stratified_cv <- function(df, rf_fun, return_preds=FALSE, resample_test=TRUE, nprocesses=N_THREADS){
	
	if(nprocesses < 2){
		
		if(! return_preds){
			out <- sapply(levels(df$animal), function(x)stratified_cv_one_level(x, df, rf_fun,return_preds, resample_test)$error)
		}
		else{
			dfs <- lapply(levels(df$animal), function(x)stratified_cv_one_level(x, df, rf_fun,return_preds, resample_test))
			out <- do.call("rbind", dfs)
		}
		}
	else{
		cl <- makeCluster(nprocesses)
		clusterExport(cl, "randomForest")
		clusterExport(cl, "stratified_cv_one_level")
		clusterExport(cl, "strip_df_for_ml")
		clusterExport(cl, "entropy_conf")
#~ 		
#~ 		clusterExport(cl, "rf_fun")
		tryCatch({clusterExport(cl, "rf_fun")},
		error = function(e){
			return(NA)
			})
		
		if(! return_preds){
			out <- parSapply(cl, levels(df$animal), function(x)stratified_cv_one_level(x, df, rf_fun,return_preds, resample_test)$error)
		}
		else{
			dfs <- parLapply(cl, levels(df$animal), function(x)stratified_cv_one_level(x, df, rf_fun,return_preds, resample_test))
			out <- do.call("rbind", dfs)
		}
			
		stopCluster(cl)
		
		}
	
	return(out)
	}
select_variables <- function(df, scale=1.2){
	
	# importantly!!, the importance has to be computed on a balanced sample
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
make_window_features <- function(winsize, d, dummy){
	kernel = rep(1/winsize, winsize)
	
#~     out <- apply(d, 2, 
#~ 		function(v,k){
#~ 			runmed(v, k)
#~ 			}, winsize)
    out <- apply(d, 2, 
		function(v,k){
			runmed(v, k)
			}, winsize)
    colnames(out) <- sprintf("win_%i.%s", winsize, colnames(out))
    return(out)
    
    }
    
    
add_lagged_time_features <- function(d, min_max_lag=c(-2,2),pattern_valid_cols = "\\.", use_window=FALSE){
    if(! use_window)
		if (min_max_lag[2] == 0  & min_max_lag[1] ==0) 
			return (d)
    # we ensure data is sorted by time
    d <- d[order(d$t),]
    # we ensure we have an homogenous time series:
    if(!use_window)
		stopifnot(length(unique(diff(d$t))) == 1)
	
    cols_to_lag_idx <- grep(pattern_valid_cols, colnames(d))
    
    cols_to_lag <- d[,cols_to_lag_idx ]
    
    other_cols <- subset(d,select = -cols_to_lag_idx)
    
    min_lag <-  min_max_lag[1]
    max_lag <- min_max_lag[2]
    
    
    if(!use_window)
		stopifnot(min_lag < max_lag)
    
    
    if(use_window){  
		dfs <- lapply(min_max_lag, make_window_features, d=cols_to_lag, min_max_lag )
		all_lagged_cols <- do.call("cbind",dfs)
		out <- na.omit(cbind(other_cols, all_lagged_cols))
		}
	else{
		other_cols <- other_cols[(abs(min_lag) + 1): (nrow(other_cols) - abs(max_lag)),, drop=F]
    
		dfs <- lapply(min_lag : max_lag, make_time_features, d=cols_to_lag, min_max_lag )
		all_lagged_cols <- na.omit(do.call("cbind",dfs))
		out <- cbind(other_cols, all_lagged_cols)
    }
    
    
    return(out)
}

add_mean_features <- function(d, pattern_valid_cols = "\\."){
    
    cols_to_use_idx <- grep(pattern_valid_cols, colnames(d))
    
    cols_to_use <- d[,cols_to_use_idx ]
    
    other_cols <- subset(d,select = -cols_to_use_idx)
    

	medians <- apply(cols_to_use,2, mean)
	medians <- matrix(medians,nrow=nrow(cols_to_use), ncol=length(medians) ,byrow=T)
	colnames(medians) <- sprintf("_median.%s", colnames(cols_to_use))
	
	out <- cbind(d, medians)
    return(out)
}
####   df <- add_lagged_time_features_to_whole(df, min_max_lag=c(1,3))
add_mean_features_to_whole <- function(df, pattern_valid_cols = "\\."){
	l_dfs <- split(df, df$animal)
	l_dfs <- lapply(l_dfs,add_mean_features, pattern_valid_cols=pattern_valid_cols)
	return(do.call("rbind", l_dfs))
	}


add_lagged_time_features_to_whole <- function(df,min_max_lag=c(-2,2),pattern_valid_cols = "\\.", use_window=FALSE){
	l_dfs <- split(df, df$animal)
	l_dfs <- lapply(l_dfs,add_lagged_time_features, min_max_lag=min_max_lag,pattern_valid_cols=pattern_valid_cols,use_window=use_window)
	return(do.call("rbind", l_dfs))
	}

analize_error_vs_lag <- function(df, max_tau, replicates=3, use_window =FALSE){
	
	if(use_window == TRUE){
		tau_list <- max_tau
		dd_fun <- function(df, tau){
			add_lagged_time_features_to_whole(df, min_max_lag=tau, use_window=TRUE)
			}
		}
	if(use_window == FALSE){
		 tau_list <- 0:max_tau
		dd_fun <- function(df, tau){
			add_lagged_time_features_to_whole(df, min_max_lag=c(-tau, +tau), use_window=FALSE)
			}
		}
		
	for (i in 1:replicates){
		for (tau in tau_list){
			
		dd <- dd_fun(df,tau)
		# importantly!! importance has to be computed on a balanced sample
			rf_fun <- function(dd){
				d <- strip_df_for_ml(dd);
				return(randomForest(y ~ ., d, ntree=50, sampsize=c(1000,1000,1000)))
				}
			
			errors <- stratified_cv(dd, rf_fun)
			if(exists("error_lag_df"))
				error_lag_df <- rbind(error_lag_df, data.frame(error=errors, lag=paste0(tau, collapse="_")))
			else
				error_lag_df <- data.frame(error=errors, lag=paste0(tau, collapse="_"))
			
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
	
analize_error_vs_nvar_comb <- function(df, replicates=3, targets=c(1,2,3,4,5,10,15,20)){
	l <- lapply(1:replicates, function(dummy)select_comb_vars(df,targets))
	error_nvar_df <- do.call("rbind", l)
	boxplot(error ~ n_vars, error_nvar_df, pch=20, xlab=expression(italic(p)), ylab=c("Average cross validation class error"),ylab=c(0,0.15))
	return(error_nvar_df)
	}
	
analize_error_vs_confidence <- function(df){
	rf_fun <- function(dd){
		d <- strip_df_for_ml(dd)
		# 200 trees for accuracy of confidence
		rf <-randomForest(y ~ ., d, ntree=200, sampsize=c(1000,1000,1000))
		return(rf)
		}
	dd <- stratified_cv(df, rf_fun, TRUE)
	
	conf_bins <- 0:9/10
	error_vs_conf <- sapply(conf_bins, function(c, d){
		sub_d <- subset(d, conf>c &  conf <= c+0.1)
		errors <- sub_d$y != sub_d$preds
		return(sum(errors) / nrow(sub_d))
		}, dd)
	barplot(error_vs_conf, width=1,space=0, names.arg=conf_bins, ylim=c(0,0.5), ylab= "Actual (cross validation) error", xlab="Prediction confidence")
	hist(dd$conf, nclass=20, col="grey", freq=F, xlab="Prediction confidence", main="")
	return(list(error_vs_conf=error_vs_conf, data=dd))
}
analize_predicted_ts_features <- function(df){

	# 3 strip series
	# 4 predict
	rf_fun <- function(dd){
		d <- strip_df_for_ml(dd);
		return(randomForest(y ~ ., d, ntree=100, sampsize=c(5000,500,5000)))
		}
	preds <- stratified_cv(df, rf_fun, resample_test=F, return_preds=T)
	preds$y <- preds$preds
	preds$preds <- NULL
	dd <- strip_df_for_ml(df) 
	preds$true_y <- dd$y
	

#~ 	for(a in levels(preds$animal)){
#~ 		subdd <- dd[preds$animal == a,]
#~ 		subpreds <- subset(preds, animal == a)
#~ 
#~ 		show_2d_hist(log10(subdd$win_3.EEG_parietal_frontal_cD_6.power.mean / subdd$win_3.EEG_parietal_frontal_cD_1.power.mean), log10(subdd$win_3.EMG_REF_cD_1.power.mean), subpreds$y,main=a)
#~ 		show_2d_hist(log10(subdd$win_3.EEG_parietal_frontal_cD_6.power.mean / subdd$win_3.EEG_parietal_frontal_cD_1.power.mean), log10(subdd$win_3.EMG_REF_cD_1.power.mean), subpreds$true_y,main=a)
#~ 		
#~ 	}

	anim_dfs = split(preds,preds$animal)
	pred_time_series <- sapply(anim_dfs, function(x)x$y)
	true_time_series <- sapply(anim_dfs, function(x)x$true_y)
	
	pred_prevalences <- sapply(pred_time_series, function(x){out <- table(x); return(out/sum(out))})	
	true_prevalences <- sapply(true_time_series, function(x){out <- table(x); return(out/sum(out))})
	
	stat_df <- expand.grid(class=rownames(true_prevalences), animal=colnames(true_prevalences), prediction=c(T,F))
	stat_df <- cbind(stat_df,prevalence=c(as.numeric(pred_prevalences), as.numeric(true_prevalences)))
	boxplot( prevalence ~ prediction * class, stat_df, col=c("grey", "red"),  ylab="Prevalence of sleep stages")
	# no effect of prediction & no interaction => no evidence of difference in prevalence
	mod <- betareg( prevalence ~ prediction * class, stat_df)
	print(summary(mod))

	
	stat_df$number <- c(
		sapply( pred_time_series, function(x)epoch_lengths_and_number(x)$number),
		sapply( true_time_series, function(x)epoch_lengths_and_number(x)$number))
		
	# effect of interaction  => not as many event -> some events are more frequents
	boxplot( number ~  prediction * class, stat_df, col=c("grey", "red"), ylab="Number of events")
	mod <- glmer( number ~  prediction *class + (1|animal), stat_df, family="poisson")
	print(summary(mod))
	
	
	stat_df$length <- c(
		sapply( pred_time_series, function(x)epoch_lengths_and_number(x)$length$length),
		sapply( true_time_series, function(x)epoch_lengths_and_number(x)$length$length))
	
	stat_df$number <- c(
		sapply( pred_time_series, function(x)epoch_lengths_and_number(x)$number),
		sapply( true_time_series, function(x)epoch_lengths_and_number(x)$number))
	# no effect of interaction nor prediction approximatively the same length
	boxplot( length ~  prediction * class, stat_df, col=c("grey", "red"), ylab="Length of events")
	mod <- lmer( length ~  prediction *class + (1|animal), stat_df)
	print(summary(mod))
	
	#todo use caret to compute average sensitivity...
	
#~ 	
	pred_transit_mats <-lapply(pred_time_series, function(x)markovchainFit(data=x)$estimate[,])
	pred_trans_array <- abind(pred_transit_mats, along=3)
	print("Predicted average transition matrix:")
	print(round(apply(pred_trans_array, c(1,2), mean),3))
	print(round(apply(pred_trans_array, c(1,2), sd),3))
	
		
	true_transit_mats <-lapply(true_time_series, function(x)markovchainFit(data=x)$estimate[,])
	true_trans_array <- abind(true_transit_mats, along=3)
	print("Empirical average transition matrix:")
	print(round(apply(true_trans_array, c(1,2), mean),3))
	print(round(apply(true_trans_array, c(1,2), sd),3))
	
	
	###
#~ 	dd$preds <-preds$y
	return(stat_df)
	}





main <- function(input_file_prefix){
	print(input_file_prefix)
	prefix <- format(Sys.time(), "%Y_%m_%d-%X")
	pdf(sprintf("/tmp/Rout_%s.pdf",prefix)); sink(sprintf("/tmp/Rout_%s.txt",prefix))	
	dfo <- cache_load_file(input_file_prefix)
	dfo <- curate_df(dfo, F)

#~ 	dfo <- curate_df(dfo,F)
	
	
#~ 	
	describe_time_series(dfo)
	error_nvar_df <- analize_error_vs_nvar(dfo,replicates=REPLICATES)
#~ 	error_nvar_df$method = "recursive"
#~ 	error_nvar_df_comb <- analize_error_vs_nvar_comb(dfo,replicates=REPLICATES, COMB_NVARS)
#~  	error_nvar_df_comb$method = "comb"
#~  	
#~  	error_nvar_df <- rbind(error_nvar_df, error_nvar_df_comb)
#~  	
 	
#~ 	out <- recursive_variable_elimination(dfo, target=25)
#~ 	
	
	
	
	
	non_predictor_cols <- grep("\\.", colnames(dfo), invert=T, value=TRUE)
	df <- dfo[,c(VARIABLE_TO_USE, non_predictor_cols)]
	
	
	
	
    error_lag_df_window <- analize_error_vs_lag(df, list(
									c(1),
									c(1,3),
									c(1,3,7),
									c(1,3,7, 15),
									c(1,3,7, 15, 31),
									c(1,3,7, 15, 31, 61),
									c(1,3,5, 9, 17, 33, 65,129)),
									 replicates=REPLICATES, use_window=TRUE)
	
#~     error_lag_df_l <- analize_error_vs_lag(df, 3, replicates=REPLICATES)
	
	
#~ 	df <- add_lagged_time_features_to_whole(df, min_max_lag=c(-TAU_TO_USE, +TAU_TO_USE))
#~ 	df <- add_lagged_time_features_to_whole(df, min_max_lag=c(-TAU_TO_USE, +TAU_TO_USE))
	df <- add_lagged_time_features_to_whole(df, min_max_lag=WINDOW_SIZE_TO_USE, use_window=TRUE)
#~ 	dfm <- add_mean_features_to_whole(df)
	err_conf_data <- analize_error_vs_confidence(df)
	
	stat_df <- analize_predicted_ts_features(df)
	
	#todo:
#~ 	CV transition matrices and prevalence...

	dev.off();sink()
	return(list(error_nvar_df, err_conf_data,error_lag_df_window, stat_df))
}



VARIABLE_TO_USE <- {c(
#~ 	"EEG_parietal_frontal_cD_1.power.mean",        
#~ 	"EEG_parietal_frontal_cD_1.power.min",        
#~ 	
#~ 	"EEG_parietal_frontal_cD_6.power.mean",        
#~ 	"EEG_parietal_frontal_cD_6.power.min",        
#~ 	
#~ 	"EEG_parietal_frontal_cD_5.hjorth.complexity",   
#~ 	"EEG_parietal_frontal_cD_5.hjorth.morbidity",
#~ 	
#~ 	"EEG_parietal_frontal_cD_5.power.mean",
#~ 	"EEG_parietal_frontal_cD_5.power.min",
#~ 	
#~ 	"EMG_REF_cD_1.power.min",                                          
#~ 	"EMG_REF_cD_1.power.mean",                                          
#~ 	"EMG_REF_cD_2.power.min",                       
#~ 	"EMG_REF_cD_2.power.mean",                       
#~ 	"EMG_REF_cD_3.power.min",
#~ 	"EMG_REF_cD_3.power.mean"
#~ 	
#~ 	
#~ 	
)}


VARIABLE_TO_USE <- {c(
"EEG_parietal_frontal_cA_6.power.mean",
"EEG_parietal_frontal_cA_6.power.min",
"EEG_parietal_frontal_cD_1.power.mean",     
"EEG_parietal_frontal_cD_1.power.min",     
"EEG_parietal_frontal_cD_2.hjorth.complexity",
"EEG_parietal_frontal_cD_2.hjorth.morbidity", 
"EEG_parietal_frontal_cD_5.hjorth.complexity",
"EEG_parietal_frontal_cD_5.hjorth.morbidity", 
"EEG_parietal_frontal_cD_6.power.mean", 
"EEG_parietal_frontal_cD_6.power.min",        
"EEG_parietal_frontal_cD_6.power.sd",         
"EMG_REF_cD_1.power.mean",                    
"EMG_REF_cD_1.power.min",                     
"EMG_REF_cD_2.power.mean",                    
"EMG_REF_cD_2.power.min",                    
"EMG_REF_cD_3.power.mean",                    
"EMG_REF_cD_3.power.min")}

TAU_TO_USE <- 1
WINDOW_SIZE_TO_USE <- c(1,3,7, 15, 31, 61)
COMB_NVARS <- c(1,2,3,4,5,10,15,20)
REPLICATES <- 10
N_THREADS <- 3
#~ args <- commandArgs(trailingOnly = TRUE)
#~ main(args[1])
out <- main("/data/pyrem/Ellys/all_features_e5s")

############################################################################

###finalrf:
	
#~ 	dd <- add_lagged_time_features_to_whole(df, min_max_lag=c(-2, +2))
#~ 	d <- strip_df_for_ml(dd)
#~ 	# importantly!! importance has to be computed on a balanced sample
#~ 	rf <- randomForest(y ~ ., d, ntree=50)
#~ 	
#~ }

#~ plot3d(x$EEG_parietal_frontal_cD_2.hjorth.morbidity ,log10(x$EMG_REF_cD_1.power.mean), log10(x$EEG_parietal_frontal_cD_5.power.mean),  col=y, size=5)
#~ 

#sleep stages? within non-rem?
#~ plot3d(x=log10(df$EMG_REF_cD_1.power.mean), y=log10(df$EEG_parietal_frontal_cD_6.power.mean), z=log10(df$EEG_parietal_frontal_cD_1.power.mean), col=df$y)

#~ xp <- as.numeric(pred_time_series[[2]])[1000:1100]
#~ xt <- as.numeric(true_time_series[[2]])[1000:1100]
#~ plot(1:length(xt), rep(1, length(xt)), col=xt, pch=20,cex=1)
#~ points(1:length(xp), rep(1.2, length(xt)), col=xp, pch=20,cex=1)
#~ points(1:length(xp), rep(0.8, length(xt)), col=ifelse(xt==xp,0,1), pch=20,cex=1)


#~ 

#~ df2 <- add_lagged_time_features_to_whole(df, min_max_lag=WINDOW_SIZE_TO_USE, use_window=TRUE)
#~ 
#~ df2 <- strip_df_for_ml(df2, keep_animal=T)
#~ 
#~ anim <- "TelC_3"
#~ train <- df2[df2$animal == anim,]
#~ test =df2[df2$animal == anim,]
#~ df2$animal <- NULL
#~ test$animal <- NULL
#~ rf <- randomForest(subset(df2, select=-y),df2$y, 
#~ 		sampsize=c(5000,500,5000),
#~ 		ytest=test$y, xtest=subset(test, select=-y), 
#~ 		 ntree=100, do.trace=T )
		 
#~ 
#~ rf <- randomForest(subset(train, select=-y),train$y, 
#~ 		ytest=test$y, xtest=subset(test, select=-y), 
#~ 		keep.forest=T,
#~ 		sampsize=c(5000,500,5000),
#~ 		 ntree=100, do.trace=T )
#~ 
#~ ani_df <- subset(df, animal==anim)
#~ new_y <- predict(rf, ani_df)
#~ dfout <- cbind(ani_df, new_y)
#~ towrite <- dfout[, c("t","y", "new_y")]
#~ towrite$eq <- ifelse(as.character(towrite$y) == as.character(towrite$new_y), "", "!")
#~ write.table(towrite, "/tmp/out.tsv", quote=F, row.names=F)

#~ show_2d_hist(log10(dfout$win_7.EEG_parietal_frontal_cD_6.power.mean), log10(dfout$win_7.EMG_REF_cD_1.power.mean), dfout$new_y, dfout$new_y, NULL)
#~ show_2d_hist(log10(dfout$win_7.EEG_parietal_frontal_cD_6.power.mean), log10(dfout$win_7.EMG_REF_cD_1.power.mean), dfout$y, dfout$y, NULL)
#~ plot3d(log10(dfout$win_1.EEG_parietal_frontal_cD_6.power.mean), log10(dfout$win_1.EEG_parietal_frontal_cD_1.power.mean), log10(dfout$win_1.EMG_REF_cD_1.power.mean), col=dfout$y)
#~ plot3d(log10(dfout$win_7.EEG_parietal_frontal_cD_6.power.mean), log10(dfout$win_7.EEG_parietal_frontal_cD_2.hjorth.complexity), log10(dfout$win_7.EMG_REF_cD_1.power.mean), col=dfout$new_y)
