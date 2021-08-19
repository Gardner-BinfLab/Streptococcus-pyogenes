install.packages("remotes")
remotes::install_github("choisy/cutoff")

library(cutoff)
d <- read.csv('~/Documents/winter_projects/spyogenes/doc/protein.99.eggnog.txt', sep='\t', header=FALSE)
names(d) <- c('RefSeq','eggNOG','Evalue','Bitscore')


hist(log(d$Bitscore), 100, F, xlab="Bitscore", ylab="Density", main=NULL, col="grey")
lines(density(log(d$Bitscore)), lwd=1.5, col="blue")

(bitscore <- em(d$Bitscore,"lognormal","lognormal"))
