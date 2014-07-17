__author__ = 'quentin'

import glob
import os
import pandas as pd
import pyrem as pr
from multiprocessing import Pool

DATA_FILE_PATTERN= "/data/pyrem/eds/all/*.csv"
SAMPLING_RATE = 200
OUT_CSV = "/tmp/test.csv"


def features_one_file(f):
    file_name = os.path.basename(f).split(".")[0]
    animal, dose = file_name.split("_")

    print animal, dose

    pol = pr.polygraph_from_csv(f,SAMPLING_RATE)

    tmp_df = feature_feactory.make_features_for_epochs(pol,5,1)

    tmp_df["animal"] = animal
    try:
        tmp_df["dose"] = float(dose)
    except:
        tmp_df["dose"] = "NA"

    return tmp_df

if __name__ == "__main__":

    files = glob.glob(DATA_FILE_PATTERN)

    feature_feactory = pr.features.FeatureFactory([
        pr.features.PeriodFeatures(),
        pr.features.PowerFeatures(),
        pr.features.NonLinearFeatures(),
        pr.features.EntropyFeatures(),
        pr.features.HjorthFeatures(),
        pr.features.WaveletsFeaturesDB4(),


    ])

    p = Pool(8)
    dfs = p.map(features_one_file, sorted(files))
    out_df = pd.concat(dfs)

    out_df.to_csv(OUT_CSV, float_format="%e")

""" R>
library("randomForest")
df = read.csv("/tmp/test.csv")
df <- subset(df, power_kurtosis != Inf)
df$X <- NULL; df$channel <- NULL
#df$rand <- rnorm(nrow(df))
#df$dose <- as.factor(df$dose)

dfC <- subset(df, animal=="C")
dfA <- subset(df, animal=="A")
dfB <- subset(df, animal=="B")


rfA = randomForest(dose ~ ., dfA, ntree=1000)

rfB = randomForest(dose ~ .,dfB,
                #xtest = subset(dfA, select=-c(dose)), ytest=dfA$dose,
                ntree=1000)

rf = randomForest(dose ~ .,na.omit(df),
                ntree=1000)


pdf("/tmp/out.pdf", w=12, h=9)

# plot(rfA$test$predicted ~ dfB$dose, pch=20)
# plot(rfB$test$predicted ~ dfA$dose, pch=20)


dfC <- subset(df, animal=="C" & is.na(dose), select=-dose)

varImpPlot(rf)


dff = subset(df, animal != "C")
dff$animal <- as.character(dff$animal)
par(mfrow=c(2,1))
boxplot(power_mean ~ dose, dff,pch=20, ylab="mean power", xlab="dose o/oo")
boxplot(power_mean ~ animal*dose, dff,pch=20, ylab="mean power", xlab="dose o/oo", col=c(grey(0.3),grey(0.7)))

boxplot(welch_median ~ dose, dff,pch=20, ylab="median frequency", xlab="dose o/oo")
boxplot(welch_median ~ animal*dose, dff,pch=20, ylab="median frequency", xlab="dose o/oo", col=c(grey(0.3),grey(0.7)))


boxplot(entropy_sample_2_1000 ~  dose,dff ,pch=20, ylab="Sample entropy (m=2, tau=1, r = 1000)", xlab="dose o/oo")
boxplot(entropy_sample_2_1000 ~ animal*dose, dff,pch=20, ylab="Sample entropy (m=2, tau=1, r = 1000)", xlab="dose o/oo", col=c(grey(0.3),grey(0.7)))


boxplot(power_sd ~  dose,dff ,pch=20, ylab="sd of power", xlab="dose o/oo")
boxplot(power_sd~ animal*dose, dff,pch=20, ylab="sd of power", xlab="dose o/oo", col=c(grey(0.3),grey(0.7)))
dev.off()

#plot(entropy_sample_2_1000 ~ dose, df)
#plot (predict(rfB, dfA) ~ dfA$dose)
"""