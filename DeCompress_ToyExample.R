library(rafalib)
library(pheatmap)

set.seed(3)

target <- matrix(
  c(1,1,0,0,
    1,1,0,0,
    1,0,0,0,
    0,0,1,0,
    0,1,0,1),
  ncol=4, byrow=TRUE)

makeTargetSamples <- matrix(
  rnorm(160,
        mean=c(rep(c(1,0),c(10,30)),
             rep(c(0,1,0),c(10,10,20)),
             rep(c(0,1,0),c(20,10,10)),
             rep(c(0,1),c(30,10))), sd=.1),
  ncol=4)

targetSamples <- t(makeTargetSamples %*% t(target))
rownames(targetSamples) <- 1:5
#imagemat(targetSamples)

colnames(targetSamples) <- 1:40

ref <- matrix(
  c(0,0,0,
    0,0,0,
    1,0,0,
    0,1,0,
    0,0,1,
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
rownames(combined) <- 1:11

src <- rep(c("target","reference"),c(40,30))
annoCombined <- data.frame(
  source=factor(src, levels=unique(src)),
  row.names=1:70)

geneCombined <- data.frame(
  geneGroup=factor(rep(LETTERS[1:6], c(2,1,1,1,3,3))),
  row.names=1:11)

ann_colors = list(
  geneGroup = RColorBrewer::brewer.pal(6, "Dark2"),
  source = c(target="slateblue",reference="grey"),
  cluster1 = c("1"="grey","2"="tomato2"),
  cluster2 = c("1"="grey","2"="mediumorchid2")
)
names(ann_colors$geneGroup) <- LETTERS[1:6]

png(file="suppfig1-1.png", width=600, height=900, res=125)
pheatmap(combined, 
         annotation_col = annoCombined,
         annotation_row = geneCombined,
         cluster_rows = FALSE, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         border=NA, show_colnames = FALSE, show_rownames=FALSE)
dev.off()

project <- matrix(
  c(1,0,0,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,0,0,0,
    0,0,1,0,0,1,1,1,0,0,0,
    0,0,0,1,0,1,1,1,0,0,0,
    0,0,0,0,1,0,0,0,1,1,1) + rnorm(55,0,.05),
  ncol=11, byrow=TRUE) 

projectedSamples <- t(t(targetSamples) %*% project)
colnames(projectedSamples) <- 1:40
rownames(projectedSamples) <- 1:11

#imagemat(cbind(projectedSamples, refSamples))

colAnno <- data.frame(cluster1 = factor(rep(c(1,2,1,2), each=10)),
                      cluster2 = factor(rep(1:2,each=20)),
                      row.names=1:40)

geneTarget <- data.frame(
  geneGroup=factor(rep(LETTERS[1:4], c(2,1,1,1))),
  row.names=1:5)

png(file="suppfig1-2.png", width=600, height=900, res=125)
pheatmap(targetSamples, 
         cluster_rows=FALSE, 
         show_colnames=FALSE,
         show_rownames=FALSE,
         annotation_col=colAnno,
         annotation_row=geneTarget,
         annotation_colors = ann_colors,
         border_color=NA)
dev.off()

png(file="suppfig1-3.png", width=600, height=900, res=125)
pheatmap(projectedSamples, 
         cluster_rows=FALSE, 
         show_colnames=FALSE,
         show_rownames=FALSE,
         annotation_col=colAnno,
         annotation_row=geneCombined,
         annotation_colors = ann_colors,
         border_color=NA)
dev.off()
