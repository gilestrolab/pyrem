rm(list=ls())
library(randomForest)


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

strip_df_for_ml <- function(df){
	d <- subset( df, py == 1)
	d$py <- NULL 
	d$treatment <- NULL 
	d$t <- NULL 
	d$animal <- NULL 
	d$y <- droplevels(d$y)
	return(d)
	}
	
load_data <- function(input){
	df <- readRDS(file=input)
	df <- curate_df(df)
#~ 	df <-subset(df, animal == "GFP_2")
	df <- strip_df_for_ml(df)
	return(df)
	}

one_comb_rf_importance <- function(comb, x,y,...){
	p_one <- comb[1]
	p_two <- comb[2]
	
	new_y <- droplevels(y[ y == p_one | y == p_two])
	new_x <-  x[ y == p_one | y == p_two,]
	rf <- randomForest(new_x,new_y,...)

#~ 	print(rf)
	return(rf$importance)
}

one_exlusive_rf_importance <- function(level, x,y,...){

	new_y <- as.factor(ifelse(y == level, 0,1))
	rf <- randomForest(x,new_y,...)
	return(rf$importance)
}



comb_rf_importance <- function(x,y, ...){
	level_combinations <- combn(levels(y),2)
	m <- apply(level_combinations,2, one_comb_rf_importance, x,y, ...)
	cn <- apply(level_combinations,2, paste0, collapse="_")
	rownames(m) <-colnames(x)
	colnames(m) <- cn
	rank_m <- apply(-m,2,rank,ties.method="first")
	return(rank_m)
}


exclusive_rf_importance <- function(x,y, ...){
	
	m <- sapply(levels(y), one_exlusive_rf_importance, x,y, ...)
	rownames(m) <-colnames(x)
	colnames(m) <- levels(y)
	rank_m <- apply(-m,2,rank,ties.method="first")
	return(rank_m)
}

#~ combinations <- split(level_combinations, 1:nrow(level_combinations))

select_comb_variables_name <- function(comb_importance_ranks, n_uniq_vars){
	
	ordered_vars <- apply(comb_importance_ranks, 2,order)
	cns <- colnames(comb_importance_ranks)
	used_vars <- numeric()
	per_comb_vars <- list()
	for(i in 1:nrow(ordered_vars)){
		v <- ordered_vars[i,]
		
		used_vars <- c(used_vars, v)
		for(w in 1:length(v)){
			value <- v[w]
			# iif the value was not used
			if(sum(value == used_vars) ==1){
				# add new value to the var list
				per_comb_vars[[cns[w]]] <- c(per_comb_vars[[cns[w]]], value)
				
				}
			else{
				# only if this new value has not been recorded yet
				if (sum(sapply(per_comb_vars, function(vec, p)p %in% vec, value)) ==0)
					per_comb_vars[["shared"]] <- c(per_comb_vars[["shared"]], value)
				}
			}
		}

	vars <- do.call("c", lapply(per_comb_vars,function(x, max)x[1:max], n_uniq_vars))
	names(vars) <- rep(c(cns,"shared"), each=n_uniq_vars)
	vars <- na.omit(vars)
	
	out <- rownames(comb_importance_ranks)[vars]
	names(out) <- names(vars)
	
	return(out)
	}

select_comb_vars <- function(df,targets){
	d <- strip_df_for_ml(df)
	x <- subset(d, select=-y)
	y <- d$y
	comb_imp <- comb_rf_importance(x,y, ntree=50, sampsize=c(5000,5000))
	
	# importantly!!, the importance has to be computed on a balanced sample
	rf_fun <- function(dd){
				d <- strip_df_for_ml(dd);
				return(randomForest(y ~ ., d, ntree=50, sampsize=c(1000,1000,1000)))
				}
	out <- lapply(targets, 
			function(t, df, comb_imp){
				good_names <- select_comb_variables_name(comb_imp, t)
				all_extra_vars <- grep("\\.",colnames(df), invert=T, value=T)
				df2 <- df[,c(good_names,all_extra_vars)]
				cv_error_vec <- stratified_cv(df2, rf_fun)
				return(data.frame(error=cv_error_vec, n_vars=length(good_names), animal=names(cv_error_vec)))
			}, 
			df, comb_imp)
			
   
    out <- do.call("rbind", out)
    return(out)
    
	}

#~ select_comb_vars(dfo, 10)

#~ 
#~ dd <- load_data("/data/pyrem/Ellys/all_features_e5s.rds")
#~ 
#~ x <- subset(dd, select=-y)
#~ y <- dd$y
#~ 
#~ comb_imp <- comb_rf_importance(x,y, ntree=50, sampsize=c(1000,1000))
#~ 
#~ 
#~ selection <- select_comb_variables(comb_imp, 3)
#~ selected_vars <- colnames(x)[selection]
#~ names(selected_vars) <- names(selection)
#~ 
#~ n <- 5000
#~ mx <- do.call("rbind", lapply(split(x,y), function(x)x[sample(1:nrow(x),n),])); 
#~ my <- as.factor(rep(levels(y), each=n))

