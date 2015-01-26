#setwd("Desktop/SUL1/")
setwd("/Volumes/mattrich/DUNHAM/SUL1/pMR002-library/20140708/enrich/20140729/log/FY3sul1/")
pdiff <- read.table("fimo.25.pdiff.filt.txt", sep="\t", header=T)
singles <- read.table("singles.pos.class.tsv", header=T, sep="\t")
s <- singles[singles$count.0 >= 50,]
c(mean(s$slope) - 2*sd(s$slope), mean(s$slope)+2*sd(s$slope))

#v <- read.table("variants.tsv", header=T, sep="\t")
#bottom = quantile(s$slope, .05, na.rm=T)
#top = quantile(s$slope, .95, na.rm=T)
top = .05
bottom = -.05
###
plot(log(pdiff$p_wt), log(pdiff$p_diff))
###

m <- merge(s, pdiff, by="sequence")

m <- m[m$p < .0001 | m$p_wt < 0.0001,] #only p-values less than 0.0001 kept

median(m[log2(m$p_diff) > 0,]$p_diff)
plot(m$slope ~ log2(m$p_diff), type="n", xaxt="n", xlab='log2(p_mut/p_wt)', ylab="Fitness")
#plot(m$slope ~ log2(m$p_diff), type="n", xaxt="n", xlab='log2(p_mut/p_wt)', ylab="Fitness", xlim=c(-7,-2), ylim=c(.04,.1))
#plot(m$slope ~ log2(m$p_diff), type="n", xaxt="n", xlab='log2(p_mut/p_wt)', ylab="Fitness", xlim=c(2,7), ylim=c(.04,.1))
abline(v=log2(10), lty=2)
abline(v=log2(.1), lty=2)
abline(h=top, lty=2)
abline(h=bottom, lty=2)
axis(1, at=seq(-10,10,by=1))

cex2 = 1
points(m$slope ~ log2(m$p_diff), cex=cex2*.8, pch=16, col="grey")
points(m[m$class == "A",]$slope ~ log2(m[m$class == "A",]$p_diff), pch=16, col="blue", cex=cex2)
points(m[m$class == "B",]$slope ~ log2(m[m$class == "B",]$p_diff), pch=16, col="green", cex=cex2)
points(m[m$class == "C",]$slope ~ log2(m[m$class == "C",]$p_diff), pch=16, col="red", cex=cex2)
points(m[m$class == "D",]$slope ~ log2(m[m$class == "D",]$p_diff), pch=16, col="cyan", cex=cex2)
points(m[m$class == "TATA",]$slope ~ log2(m[m$class == "TATA",]$p_diff), pch=16, col="orange", cex=cex2)
points(m[m$class == "UORF",]$slope ~ log2(m[m$class == "UORF",]$p_diff), pch=16, col="black", cex=cex2)
points(m[m$class == "N" & m$slope >= top,]$slope ~ m[m$class == "N" & m$slope >= top,]$p_diff, cex=cex2)
legend("bottomright", legend=c("A", "B", "C", "TATA", "uORF", "Unannotated"), pch=16, col=c("blue", "green", "red", "orange", "black", "grey"))




log2(median(m[log2(m$p_diff) < 0,]$p_diff))
log2(median(m[log2(m$p_diff) > 0,]$p_diff))

f <- m[(m$p_diff <= -3.3 | m$p_diff >= 3.3) & (m$slope >= top | m$slope <= bottom),]) 
m[m$p_diff < .1 & m$slope > top,]
m[m$p_diff > 10 & m$slope > top,]
m[m$p_diff > 10 & m$slope < bottom,]
m[m$p_diff < .1 & m$slope < bottom,]

plot(f$slope ~ log2(f$p_diff))


