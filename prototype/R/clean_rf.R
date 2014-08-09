


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
