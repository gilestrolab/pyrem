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
	df <-subset(df, animal == "GFP_2")
	df <- strip_df_for_ml(df)
	return(df)
	}

one_comb_rf_importance <- function(comb, x,y,...){
	p_one <- comb[1]
	p_two <- comb[2]
	
	new_y <- droplevels(y[ y == p_one | y == p_two])
	new_x <-  x[ y == p_one | y == p_two,]
	rf <- randomForest(new_x,new_y,...)
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


dd <- load_data("/data/pyrem/Ellys/all_features_e5s.rds")

x <- subset(dd, select=-y)
y <- dd$y
#~ 
#~ good_vars <- apply(rank_m, 2, function(x,n){
#~ 	names(x)[match(1:n, x)]}, 5)

com_imp <- comb_rf_importance(x,y, ntree=100)



select_comb_variables <- function(comb_importance_ranks, n_uniq_vars){
	
	ordered_vars <- apply(comb_importance_ranks, 2,order)
	cns <- colnames(com_imp)
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
	
	return(na.omit(vars))
	}
select_comb_variables(com_imp)

