#Estimated Ancestral Character States Using ape and REML
library(ape)

##data accession
#internal node averages
ina <- as.matrix(read.csv("ina.csv", row.names=(1)))
inaMat <- as.matrix(read.csv("ina.csv", row.names=(1)))[,1]
inaLong <- as.matrix(read.csv("ina.csv", row.names=(1)))[,2]
inaRat <- as.matrix(read.csv("ina.csv", row.names=(1)))[,3]
inaDNDS <- as.matrix(read.csv("ina.csv", row.names=(1)))[,4]

#tree and continuous data
X <- as.matrix(read.csv("sml.csv", row.names=(1)))
X <- X[-13,]
X <- X[-length(X)]
tree <- read.tree("final_newick")
tree$edge.length[tree$edge.length == 0.000] <- 0.00001
text <- "((Rattus_norvegicus:0.0,(Urocitellus_parryii:0.06975,(((Canis_lupus_familiaris:0.04300999999999999,Neovison_vison:0.03615):0.06042,(Bos_taurus:0.07979,Monodon_monoceros:0.04198):0.013509999999999994):0.013020000000000004,(Mus_musculus:0.14374,Loxodonta_africana:0.12806):0.007280000000000009):0.0045499999999999985):0.03386):0.01,((Homo_sapiens:0.0,Pan_troglodytes:0.0):0.015329999999999996,(Macaca_mulatta:0.00519,Theropithecus_gelada:0.0):0.02538):0.010000000000000004);"
sexualMaturity <- as.matrix(read.csv("sml.csv", row.names=1))[,1]
longevity <- as.matrix(read.csv("sml.csv", row.names=1))[,2]
ratio <- as.matrix(read.csv("sml.csv", row.names=1))[,3]
dnds <- as.matrix(read.csv("sml.csv", row.names=1))[,4]

#ancestral character state estimation with REML
#"It could be shown that, with a continuous character, REML results in unbiased estimates of the variance of the Brownian motion process while ML gives a downward bias. Therefore the former is recommanded."
matREML <- ace(sexualMaturity, tree, type="continuous", method="REML",CI=TRUE)
longREML <- ace(longevity, tree, type="continuous", method="REML",CI=TRUE)
ratREML <- ace(ratio, tree, type="continuous", method="REML",CI=TRUE)
dndsREML <- ace(dnds, tree, type="continuous", method="REML", CI=TRUE)

#plotting each plot (for lines and text on points, run the plot, then the lines, then the text)
#ex: plot(), then abline(), then abline(), then text()
plot(inaLong, longREML$ace, xlim=c(10,95), ylim=c(10,95), pch=19,cex=3.5,xlab="Longevity 'Truth' in Years",ylab="REML Longevity in Years",main="Longevity Estimated Ancestral Character States")
plot(inaMat, matREML$ace, xlim=c(450,4250), ylim=c(450,4250), pch=19,cex=3.5,xlab="Sexual Maturity 'Truth' in Days",ylab="Sexual Maturity REML in Days",main="Sexual Maturity Estimated Ancestral Character States")
plot(inaRat, ratREML$ace, xlim=c(8,21), ylim=c(8,21), pch=19,cex=3.5,xlab="Ratio (SM/L) 'Truth'",ylab="Ratio (SM/L) REML",main="Ratio (SM/L) Estimated Ancestral Character States")
plot(inaDNDS, dndsREML$ace, xlim=c(0,.8), ylim=c(0,.8), pch=19,cex=3.5,xlab="dN/dS 'Truth'",ylab="dN/dS REML",main="dN/dS Estimated Ancestral Character States")

#plot straight line (what it should look like, if REML is perfect)
abline(0,1,col="Orange",lwd=2, lty=2)

#plot linear regression
REML_long <- lm(longREML$ace ~ inaLong)
abline(REML_long, lwd=2, col="Blue")
REML_mat <- lm(matREML$ace ~ inaMat)
abline(REML_mat, lwd=2, col="Blue")
REML_rat <- lm(ratREML$ace ~ inaRat)
abline(REML_rat, lwd=2, col="Blue")
REML_dnds <- lm(dndsREML$ace ~ inaDNDS)
abline(REML_dnds, lwd=2, col="Blue")

#put text on the points
text(longREML$ace ~inaLong, labels=names(inaLong), col="White")
text(matREML$ace ~inaMat, labels=names(inaMat), col="White")
text(ratREML$ace ~inaRat, labels=names(inaRat), col="White")
text(dndsREML$ace ~inaDNDS, labels=names(inaDNDS), col="white")

#print the 95 confidence interval
longREML$CI95
matREML$CI95
ratREML$CI95
dndsREML$CI95

#linear regression summaries
long_summary <- summary(REML_long)
mat_summary <- summary(REML_mat)
rat_summary <- summary(REML_rat)
dnds_summary <- summary(REML_dnds)
