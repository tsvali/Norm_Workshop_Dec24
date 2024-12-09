### R code from vignette source 'ROTS.Rnw'

###################################################
### code chunk number 1: load
###################################################
library(ROTS)
data(upsSpikeIn)
input = upsSpikeIn


###################################################
### code chunk number 2: ROTS.Rnw:118-120
###################################################
groups = c(rep(0,3), rep(1,3))
groups


###################################################
### code chunk number 3: ROTS
###################################################
results = ROTS(data = input, groups = groups , B = 100 , K = 500 , seed = 1234)
names(results) 


###################################################
### code chunk number 4: summary
###################################################
summary(results, fdr = 0.05)


###################################################
### code chunk number 5: volcano
###################################################
plot(results, fdr = 0.05, type = "volcano")


###################################################
### code chunk number 6: heatmap
###################################################
plot(results, fdr = 0.05, type = "heatmap")


