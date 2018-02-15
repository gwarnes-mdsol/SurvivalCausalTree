# SurvivalCausalTree

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Linux Build Status](https://travis-ci.org/gwarnes-mdsol/SurvivalCausalTree.svg)](https://travis-ci.org/gwarnes-mdsol/SurvivalCausalTree)
[![Coverage Status](https://codecov.io/gh/gwarnes-mdsol/SurvivalCausalTree/branch/master/graph/badge.svg)](https://codecov.io/gh/gwarnes-mdsol/SurvivalCausalTree)


This repository implemented the survival causal tree (SCT) method proposed in Bioinformatics paper. The main computation part is written in C which significantly improves efficiency. The old code is at http://nugget.unisa.edu.au/Thus/SCT.zip

This code is based on rpart package from CRAN, and Susan Athey's causalTree package at https://github.com/susanathey/causalTree.

In addition, the propensity score in causalTree package has to be the same for all samples, in this package we fixed this problem so that it can be specified as a vector.

Use the script in the "tests" folder for a demo.
                          
