rm(list=ls())
DBG = 1
printdbg <- function(){
    print(DBG)
    assign("DBG", DBG+1, envir = .GlobalEnv)
    
}

library("randomForest")


OUT_CSV = "/data/pyrem/Ellys/all_features.csv"

curate_df <- function(dfo){
    df <- dfo
    # remove erroneous kurtosis
    #df <- subset(df, power_kurtosis != Inf)
    #df <- subset(df, power_kurtosis != Inf)

    # remove ambiguous vigilance states
    df <- subset( df, vigilance_state.vigil.proba == 1)

    df$vigilance_state.vigil.proba <- NULL
    df$t <- df$X
    df$X <- NULL
    df$animal <- sprintf("%s_%s", df$treatment, df$animal)
    df <- subset(df, animal != "TelC_E")
    df$animal  <- as.factor(df$animal)

    df$y <- as.factor(df[,"vigilance_state.vigil.value"])
    df$vigilance_state.vigil.value <- NULL
    s = sapply(df, function(col)sum(is.na(col)))
    df <- df[,s ==0]
    return(df)
}

########################################################

dfo <- read.csv(OUT_CSV , na.string="NaN")
df <- curate_df(dfo)



crossval_test <- function(out_level, original_df){
        printdbg()
        train_df <- subset(original_df, animal !=out_level)
        test_df <- subset(original_df, animal ==out_level)
        train_df$animal <- NULL
        test_df$animal <- NULL
        printdbg()
        rf <- randomForest(y ~ ., train_df, ntree=50)
        printdbg()
        print(rf)
        varImpPlot(rf)
        preds <- predict(rf, test_df)

        #out <- data.frame(real = test_df$y, preds = preds)
        out <- sum(test_df$y == preds) / length(preds)
        return(out)
}

l = sapply(levels(df$animal), crossval_test, original_df=df)
print(l)
exit()
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

