#looking for phylogenetic correlation
library(phytools)

#tree and continuous data
X <- as.matrix(read.csv("sml.csv", row.names=(1)))
X <- X[-13,]
tree <- read.tree("final_newick")
tree$edge.length[tree$edge.length == 0.000] <- 0.00001
text <- "((Rattus_norvegicus:0.0,(Urocitellus_parryii:0.06975,(((Canis_lupus_familiaris:0.04300999999999999,Neovison_vison:0.03615):0.06042,(Bos_taurus:0.07979,Monodon_monoceros:0.04198):0.013509999999999994):0.013020000000000004,(Mus_musculus:0.14374,Loxodonta_africana:0.12806):0.007280000000000009):0.0045499999999999985):0.03386):0.01,((Homo_sapiens:0.0,Pan_troglodytes:0.0):0.015329999999999996,(Macaca_mulatta:0.00519,Theropithecus_gelada:0.0):0.02538):0.010000000000000004);"
sexualMaturity <- as.matrix(read.csv("sml.csv", row.names=1))[,1]
longevity <- as.matrix(read.csv("sml.csv", row.names=1))[,2]
ratio <- as.matrix(read.csv("sml.csv", row.names=1))[,3]
dnds <- as.matrix(read.csv("sml.csv", row.names=1))[,4]

obj <- phyl.vcv(X,vcv.phylo(tree),1)

r.matlong <- cov2cor(obj$R)["dNdS", "Longevity"]
t.matlong <- r.matlong*sqrt((Ntip(tree)-2)/(1-r.matlong^2))
P.matlong <- 2*pt(abs(t.matlong),df=Ntip(tree)-2,lower.tail=FALSE)
P.matlong

r.ratdnds <- cov2cor(obj$R)["dNdS", "Ratio"]
t.ratdnds <- r.ratdnds*sqrt((Ntip(tree)-2)/(1-r.ratdnds^2))
P.ratdnds <- 2*pt(abs(t.ratdnds),df=Ntip(tree)-2,lower.tail=FALSE)
P.ratdnds

r.matdnds <- cov2cor(obj$R)["dNdS", "Maturity"]
t.matdnds <- r.matdnds*sqrt((Ntip(tree)-2)/(1-r.matdnds^2))
P.matdnds <- 2*pt(abs(t.matdnds),df=Ntip(tree)-2,lower.tail=FALSE)
P.matdnds

ratio_dnds <- X[,-1]
ratio_dnds <- ratio_dnds[,-1]

mat_dnds <-X[,-2]
mat_dnds <-mat_dnds[,-2]

long_dnds <- X[,-1]
long_dnds <- long_dnds[,-2]


#Getting the phylogenetic generalized least squares for dnds and sexual matrurity, longevity, and their ratio
pglsRatio <- gls(dnds~ratio, correlation = corBrownian(phy = tree), method = "ML")
summary(pglsRatio)
phylomorphospace(tree, ratio_dnds)
abline(a = coef(pglsRatio)[1], b = coef(pglsRatio)[2], lwd=2, col="Blue")

pglsLongevity <- gls(dnds~longevity, correlation = corBrownian(phy = tree), method = "ML")
summary(pglsLongevity)
phylomorphospace(tree, long_dnds)
abline(a = coef(pglsLongevity)[1], b = coef(pglsLongevity)[2], lwd=2, col="Blue")

pglsMaturity <- gls(dnds~sexualMaturity, correlation = corBrownian(phy = tree), method = "ML")
summary(pglsMaturity)
phylomorphospace(tree, mat_dnds)
abline(a = coef(pglsMaturity)[1], b = coef(pglsMaturity)[2], lwd=2, col="Blue")


