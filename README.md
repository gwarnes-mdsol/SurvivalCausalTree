# SurvivalCausalTree

This repository implemented the survival causal tree (SCT) method proposed in Bioinformatics paper. The main computation part is written in C to improve efficiency.For replication purpose, please see http://nugget.unisa.edu.au/Thus/SCT.zip

This code is based on rpart package from CRAN, and Susan Athey's causalTree package at https://github.com/susanathey/causalTree.



tree <- causalTree(as.formula(paste("y~",paste(f))),
                          data=dataTrain, treatment=dataTrain$w,
                          split.Rule="survival", split.Honest=F, cv.option=cv.option.temp, minsize = minsize.temp,
                          split.alpha = 1, cv.alpha = 1, xval=0, cp=0,propensity = rep(.5,999),completeCase = rep(1,999))
                          
where completeCase is a vector storing the completeCase indicator for each sample,  propensity is a vector storing the propensity score.
