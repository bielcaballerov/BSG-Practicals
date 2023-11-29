
# ###########################################################################
# Practical 4: population substructure MDS
# ###########################################################################

rm(list=ls())
library(data.table)
library(MASS)

filename <- "/Users/work/lectures_UPC 2023_2024/MDS-BSG/lecture_practicals/04_substructure/Chr21.dat"
raw_data <- fread(filename, header = T)

SNPdata <- raw_data[,7:ncol(raw_data)]

n <- nrow(SNPdata)
p <- ncol(SNPdata)
n   # individuals
p   # SNPs variants

# 1. How many variants are there in this database? 
p

# 1. What percentage of the data is missing? 
perc.mis <- 100*sum(is.na(SNPdata))/(n*p)
perc.mis

# 2. Compute manhattan distance
D <- dist(SNPdata, method = "manhattan")
class(D)
D

D.man <- as.matrix(D)
class(D.man)

D.man[1:5,1:5]
dim(D.man)

# 4. Apply MDS with 2 dimensions
mds.out <- cmdscale(D.man,k=2,eig=TRUE)
mds.points <- mds.out$points[,1:2]

plot(mds.points[,1], mds.points[,2],asp=1,xlab="First principal axis",
     ylab="Second principal axis", main = "MDS map of individuals")
text(mds.points[,1], mds.points[,2], rownames(SNPdata), cex=0.5,pos=1)

first_PA <- mds.points[,1]
print(length(first_PA[first_PA<0]))

print(length(first_PA[first_PA>0]))

# 5. GOF
mds.gof = mds.out$GOF
mds.gof

# ALTERNATIVE: check if matrix is positive semi definite
mds.eig = mds.out$eig 

if  (sum(mds.eig<0)>=1) {
  gof_nonEuclidean <- sum(mds.eig[1:2])/sum(abs(mds.eig[1:n-1]))
} else {
  gof_euclidean <- sum(mds.eig[1:2])/sum(mds.eig[1:n-1])
}

gof_nonEuclidean
gof_euclidean

# Does the perfect representation of this distance matrix exist?
gof <- c()
for (i in 1:(n-1)){
  out <- cmdscale(D.man,k=i,eig=TRUE)
  gof <- c(gof, out$GOF[2])
}
plot(gof)


# 6. Estimated vs observed distances
estimatedD <- as.matrix(dist(mds.points))#, method = "manhattan"))
dim(estimatedD)

Dobs.vec <- D.man[lower.tri(D.man)]
Dest.vec <- estimatedD[lower.tri(estimatedD)]

fitted_line = lm(Dest.vec~Dobs.vec)

cor(Dobs.vec,Dest.vec)

plot(Dobs.vec,Dest.vec,xlab="Observed",ylab="Fitted",main='Observed vs Fitted Distances')
abline(fitted_line,col=2,lwd=2)

summary(fitted_line)$r.squared



# 7. non-metric MDS k=2
set.seed(12345)

init <- scale(matrix(runif(2*n),ncol=2),scale=FALSE)
nmds <- isoMDS(D.man,y=init,k=2, trace=FALSE)

nmds.points <- nmds$points
plot(nmds.points[,1], nmds.points[,2],pch=19, xlab='MDS-dim1', ylab='MDS-dim2')
text(nmds.points[,1], nmds.points[,2], rownames(SNPdata), cex=0.5,pos=1)


# 8. non-metric MDS with different random initializations
par(mfcol = c(3, 2))
stress.ini <- c()
for (i in 1:6) {
  init <- scale(matrix(runif(2*n),ncol=2),scale=FALSE)
  nmds <- isoMDS(D.man,y=init,k=2, trace=FALSE)
  stress.ini <- c(stress.ini, nmds$stress)
  Y <- nmds$points
  plot(Y[,1],Y[,2],pch=19, xlab='MDS-dim1', ylab='MDS-dim2')
  text(Y[,1], Y[,2], rownames(SNPdata), cex=0.5,pos=1)
}

stress.ini 

# 9. non-metric MDS for n-dimensions
par(mfcol = c(1, 1))

stress.dim <- c()

for (i in 1:50) {
  init <- scale(matrix(runif(i*n),ncol=i),scale=FALSE)
  out <- isoMDS(D.man,y=init, k=i, trace=FALSE)
  stress.dim <- c(stress.dim, out$stress)
}

plot(stress.dim,main="Stress against number of dimensions")
abline(h=10)

min(which(stress.dim<10))


# 10. Compare classic MDS to non-metric MDS
best_stress <- Inf 
best_Y <- Inf
worst_stress <- 0 
worst_Y <- Inf

for (i in 1:50) {
  init <- scale(matrix(runif(2*n),ncol=2),scale=FALSE) 
  out <- isoMDS(D.man,y=init, k=2, trace=FALSE)
  Y <- out$points
  stress <- out$stress
  if (stress > worst_stress) {
    worst_stress <- stress
    worst_Y <- Y }
  if (stress < best_stress) {
    best_stress <- stress
    best_Y <- Y }
}

par(mfcol = c(3, 1))
colors = (best_Y[,1]>0)*1+1
plot(best_Y[,1],best_Y[,2],col=colors,asp=1,xlab="First principal axis",
     ylab="Second principal axis", main = "best isoMDS map out of 100")


colors = (mds.points[,1]>0)*1+1
plot(mds.points[,1],mds.points[,2],col=colors,asp=1,xlab="First principal axis",
     ylab="Second principal axis", main = "MDS map of individuals")
text(mds.points[,1], mds.points[,2], rownames(SNPdata), cex=0.5,pos=1)



colors = (worst_Y[,1]>0)*1+1
plot(worst_Y[,1],worst_Y[,2],asp=1,xlab="First principal axis",
     ylab="Second principal axis", main = "worse isoMDS map out of 100")
text(worst_Y[,1], worst_Y[,2], rownames(SNPdata), cex=0.5,pos=1)


first_bestPA <- best_Y[,1]
print(length(first_bestPA[first_bestPA<0]))
print(length(first_bestPA[first_bestPA>0]))

first_PA <- mds.points[,1]
print(length(first_PA[first_PA<0]))
print(length(first_PA[first_PA>0]))
 
# mds
gr1.mds <- which(first_PA<0)
gr2.mds <- which(first_PA>0)

# best nmds
gr1.nmds <- which(first_bestPA<0)
gr2.nmds <- which(first_bestPA>0)

# worse nmds
first_worstPA <- worst_Y[,1]
gr1.bad.nmds <- which(first_worstPA<0)
gr2.bad.nmds <- which(first_worstPA>0)

gr1.mds
gr1.nmds
gr1.bad.nmds

gr2.mds
gr2.nmds
gr2.bad.nmds

# 11. correlation

R <- cbind(mds.points, best_Y)
colnames(R) <- c("MDS-1","MDS-2","NMDS-1","NMDS-2")


round(cor(R),digits=2)



# FINISH HERE
# ###########################################################################
# ###########################################################################



# 
image(D.man)


# HWE
genotype.counts <- MakeCounts(SNPdata)[,1:3]

gntype.gr1 <- MakeCounts(SNPdata[gr1.mds,])[,1:3]
gntype.gr2 <- MakeCounts(SNPdata[gr2.mds,])[,1:3]

# 
alpha <- 0.05
pvalues_exact <- HWExactStats(genotype.counts)
pvalues_significant_exact <- pvalues_exact[pvalues_exact <= alpha]

pvalues_exact.gr1 <- HWExactStats(gntype.gr1)
pvalues_significant_exact.gr1 <- pvalues_exact.gr1[pvalues_exact.gr1 <= alpha]

pvalues_exact.gr2 <- HWExactStats(gntype.gr2)
pvalues_significant_exact.gr2 <- pvalues_exact.gr2[pvalues_exact.gr2 <= alpha]

length(pvalues_significant_exact)/nrow(genotype.counts)*100
length(pvalues_significant_exact.gr1)/nrow(gntype.gr1)*100
length(pvalues_significant_exact.gr2)/nrow(gntype.gr2)*100

a <- 1:120
b <- a[seq(1, length(a), 6)]

sequence <-as.integer(seq(1, dim(SNPdata)[2], 200))
length(sequence)

col.names <- colnames(SNPdata[,..sequence])
length(col.names)
all_genotypes <- makeGenotypes(SNPdata, convert=col.names, method=as.genotype.allele.count)

LD_pairs <- LD(all_genotypes)

Dm <- LD_pairs$D
Dp <- LD_pairs$"D'"
R2 <- LD_pairs$"R^2"
X2 <- LD_pairs$"X^2"

image(R2, col=rainbow(25), axes=F)


