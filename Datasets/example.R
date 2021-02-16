setwd("Datasets")

y <- read.csv("AMI.csv",header=T)
post_bnb(y[,1],K=1000)

post_unif(y[,1],K=1000)

post_pois(y[,1],K=1000)
