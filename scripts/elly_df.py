__author__ = 'quentin'

import glob
import os
import pandas as pd
import pyrem as pr
from multiprocessing import Pool
import numpy as np


DATA_FILE_PATTERN= "/data/pyrem/Ellys/pkls/*.pkl"
SAMPLING_RATE = 200
OUT_CSV = "/data/pyrem/Ellys/all_features.csv"
LAG_WINDOW = 100
N_PROCESSES = 8

# DATA_FILE_PATTERN= "/data/pyrem/Ellys/pkls/GFP_*A*.pkl"
# OUT_CSV = "/tmp/all_features.csv"




def features_one_file(f):
    file_name = os.path.basename(f).split(".")[0]
    treatment, animal = file_name.split("_")

    pol = pr.polygraph_from_pkl(f)
    pol = pol.normalise()
    pol = pr.preprocess_eegs(pol)
    print "processing " + f

    tmp_df = feature_factory.make_features_for_epochs(pol,10,LAG_WINDOW, add_major_annotations=True)

    tmp_df["animal"] = animal
    tmp_df["treatment"] = treatment


    return tmp_df

if __name__ == "__main__":

    files = glob.glob(DATA_FILE_PATTERN)

    feature_factory = pr.features.FeatureFactory([
        pr.features.PeriodFeatures(),
        pr.features.PowerFeatures(),
        pr.features.NonLinearFeatures(),
        # pr.features.EntropyFeatures(),
        pr.features.HjorthFeatures(),
        pr.features.WaveletsFeaturesDB4(),
        pr.features.MSEFeatures(),


    ])
    if N_PROCESSES > 1 :
        p = Pool(N_PROCESSES)
        dfs = p.map(features_one_file, sorted(files))
    else:

        dfs = map(features_one_file, sorted(files))

    out_df = pd.concat(dfs)

    out_df.to_csv(OUT_CSV, float_format="%e")

""" R>

library("randomForest")
# library("foreach")
# library("doParallel")


######################################################
OUT_PREFIX<- "/data/pyrem/Ellys/all_features"

filename <- paste(OUT_PREFIX,"rds",sep=".")

tryCatch({

    dfo <- readRDS(file=filename)

}, error = function(e) {
    csv_filename <- paste(OUT_PREFIX,"csv",sep=".")
    dfo <- read.csv(csv_filename, na.string="NaN")
    saveRDS(dfo, file=filename)
}
)

df <- dfo

# remove erroneous kurtosis
df <- subset(df, power_kurtosis != Inf)
# remove ambiguous vigilance states
df <- subset(df, vigil_prob == 1)
df$vigil_prob <- NULL



# use only eeg 2
df <- subset(df, channel=="EEG_parietal_frontal")



df$animal <- sprintf("%s_%s", df$treatment, df$animal)
df$animal  <- as.factor(df$animal)



df$t <- df$X
df$X <- NULL
df$channel <- NULL;

#df <- subset(df, df$vigil_value==87)
#df$vigil_value <- NULL
df$vigil_value <- as.factor(df$vigil_value)

# crossval <- function(out_level, original_df){
#         train_df <- subset(original_df, animal !=out_level)
#         # train_df <- original_df
#
#         test_df <- subset(original_df, animal ==out_level)
#
#         # print(head(train_df))
#         # print(head(test_df))
#         # print(unique(train_df$treatment))
#         # print(unique(test_df$treatment))
#
#         train_df$animal <- NULL
#         test_df$animal <- NULL
#         rf <- randomForest(treatment ~ ., train_df, ntree=50)
#         varImpPlot(rf)
#         votes <- predict(rf, test_df,type="prob")
#
#         votes <- matrix(votes, nrow=length(rownames(votes)),dimnames=list(NULL,colnames(votes)))
#
#         h_rows <- apply(votes, 1, entropy)
#         weighted_votes <- colSums(votes * h_rows )
#
#
#         #
#         # weighted_votes = weighted_votes /sum(weighted_votes)
#         #
#         # # tab <- matrix(weighted_votes, nrow=1,dimnames=list(NULL,colnames(votes)))
#         #
#         out <- data.frame(pred = colnames(votes)[which.max(weighted_votes)], real=unique(test_df$treatment))
#         # out$real <- unique(test_df$treatment)
#         # out$animal <- out_level
#         print(out)
#         return(out)
# }
#
# entropy <- function(v){
#     v[v==0] <- 1e-100
#     h <- 1 + sum(v * log2(v)) /log2(length(v))
#     return(h)
# }
# set.seed(1)
# l = lapply(levels(df$animal), crossval, original_df=df)
# d = do.call("rbind",l)
# predictions = do.call("rbind", l)
# predictions <- within(predictions, pred <- ifelse(GFP > TelC, "GFP", "TelC"))
#
######################################################################################################


crossval_test <- function(out_level, original_df){
        train_df <- subset(original_df, animal !=out_level)
        test_df <- subset(original_df, animal ==out_level)
        train_df$animal <- NULL
        test_df$animal <- NULL
        rf <- randomForest(vigil_value ~ ., train_df, ntree=50)
        print(rf)
        varImpPlot(rf)
        preds <- predict(rf, test_df)

        out <- data.frame(real = test_df$vigil_value, preds = preds)
        out <- sum(test_df$vigil_value == preds) / length(preds)
        return(out)
}

l = sapply(levels(df$animal), crossval_test, original_df=df)



#
#
# dpdf("/tmp/test.pdf")
# for (a in levels(df$animal)){
#   sub_df <- subset(df, animal == a)
#   plot( log10(power_median) ~ t,sub_df, pch=20, col=treatment,ylim=c(6,10),main=a)
#   }
# dev.off()




#############
set.seed(1)
out_level <- "GFP_D"
train_df <- subset(df, animal !=out_level)
test_df <- subset(df, animal ==out_level)

train_df$animal <- NULL
test_df$animal <- NULL

ytrain <- train_df$treatment
ytest <- test_df$treatment

xtrain<- train_df
xtrain$treatment <- NULL
xtest <- test_df
xtest$treatment <- NULL

rf <- randomForest(x = xtrain, y= ytrain, ntree=100)

table(predict(rf, xtest))


"""