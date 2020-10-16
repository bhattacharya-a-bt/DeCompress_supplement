library(rafalib)
library(pheatmap)

set.seed(3)

target <- matrix(
  c(1,1,0,0,
    1,1,0,0,
    1,0,0,0,
    0,0,1,0),
  ncol=4, byrow=TRUE)

makeTargetSamples <- matrix(
  rnorm(160,
        mean=c(rep(c(1,0),c(10,30)),
             rep(c(0,1,0),c(10,10,20)),
             rep(c(0,1,0),c(20,10,10)),
             rep(c(0,1),c(30,10))), sd=.1),
  ncol=4)

targetSamples <- t(makeTargetSamples %*% t(target))
#imagemat(targetSamples)

colnames(targetSamples) <- 1:40

ref <- matrix(
  c(0,0,0,
    0,0,0,
    1,0,0,
    0,1,0,
    1,1,0,
    1,1,0,
    1,1,0,    
    0,0,1,
    0,0,1,
    0,0,1),
  ncol=3, byrow=TRUE)

makeRefSamples <- matrix(
  rnorm(90,
        mean=c(rep(c(1,0),c(10,20)),
               rep(c(0,1,0),c(10,10,10)),
               rep(c(0,1),c(20,10))), sd=.1),
  ncol=3)

refSamples <- t(makeRefSamples %*% t(ref))
#imagemat(refSamples)

combined <- cbind(rbind(targetSamples,matrix(0,ncol=40,nrow=6)), 
                  refSamples)
colnames(combined) <- 1:70

src <- rep(c("target","reference"),c(40,30))
annoCombined <- data.frame(
  source=factor(src, levels=unique(src)),
  row.names=1:70)

pheatmap(combined, annotation_col = annoCombined,
         cluster_rows = FALSE, cluster_cols = FALSE, border=NA,
         show_colnames = FALSE)

project <- matrix(
  c(1,0,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,0,0,
    0,0,1,0,1,1,1,0,0,0,
    0,0,0,1,1,1,1,0,0,0) + rnorm(40,0,.05),
  ncol=10, byrow=TRUE) 

projectedSamples <- t(t(targetSamples) %*% project)
colnames(projectedSamples) <- 1:40

#imagemat(cbind(projectedSamples, refSamples))

colAnno <- data.frame(group1 = factor(rep(c(1,2,1,2), each=10)),
                      group2 = factor(rep(1:2,each=20)),
                      row.names=1:40)
                                      
pheatmap(targetSamples, 
         cluster_rows=FALSE, 
         show_colnames=FALSE,
         annotation_col=colAnno,
         border_color=NA)

pheatmap(projectedSamples, 
         cluster_rows=FALSE, 
         show_colnames=FALSE,
         annotation_col=colAnno,
         border_color=NA)
