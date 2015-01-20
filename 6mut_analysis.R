###Before this analysis, convert the output of 


require(plyr)
setwd("/Volumes/mattrich/DUNHAM/SUL1/pMR002-library/6mut-library/20141118/")

seqlist <- data.frame(read.table("../6mut_order2.txt", colClasses=c("character"))) 
names(seqlist) <- "seqID"

d1 <- data.frame(read.table("5-1.counts.txt", sep="\t", colClasses=c("character", "numeric")))  
d2 <- data.frame(read.table("5-3.counts.txt", sep="\t", colClasses=c("character", "numeric")))  
d3 <- data.frame(read.table("5-5.counts.txt", sep="\t", colClasses=c("character", "numeric")))  
d4 <- data.frame(read.table("5-7.counts.txt", sep="\t", colClasses=c("character", "numeric")))  

names(d1) <- c("seqID", "count.1") 
names(d2) <- c("seqID", "count.3") 
names(d3) <- c("seqID", "count.5") 
names(d4) <- c("seqID", "count.7") 

dat <- join_all(list(seqlist,d1,d2,d3,d4), by="seqID")
time <- c(5.5, 19.3, 30, 41.7)

dat$freq.1 <- dat$count.1/sum(dat$count.1)
dat$freq.3 <- dat$count.3/sum(dat$count.3)
dat$freq.5 <- dat$count.5/sum(dat$count.5)
dat$freq.7 <- dat$count.7/sum(dat$count.7)
dat$ratio.1 <- log2(dat$freq.1/dat$freq.1)
dat$ratio.3 <- log2(dat$freq.3/dat$freq.1)
dat$ratio.5 <- log2(dat$freq.5/dat$freq.1)
dat$ratio.7 <- log2(dat$freq.7/dat$freq.1)

fitness <- c()
intercepts <- c()
for( i in seq(32) ){
  fits <- lm(c(dat$ratio.1[i], dat$ratio.3[i], dat$ratio.5[i], dat$ratio.7[i]) ~ time)
  fitness <- c(fitness, fits$coefficients[2])
  intercepts <- c(intercepts, fits$coefficients[1])
}

dat$intercept <- intercepts
dat$slope <- fitness
dat$norm_slope <- dat$slope - dat$slope[1]

par(mfrow=c(1,1))
b <- barplot(dat$norm_slope, names.arg=dat$seqID, las=2, ylim=c(-0.005,.12))
