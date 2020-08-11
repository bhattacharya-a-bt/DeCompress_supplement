### THIS PROVIDES AN EXAMPLE OF GENERATING PSEUDOTARGETED PANELS
### FROM GEO DATASETS. THIS EXAMPLES USES GSE19830.
### OTHER DATASETS USED: GSE97284, GSE123604, GSE64098

source('deconvolution/benchmarking_functions.R')

#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE19830", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL1355", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

exp = exprs(gset)
pData = pData(gset)

extractPercentages <- function(str){

    split = strsplit(str,'[%/]')
    return(as.numeric(split[[1]][c(1,3,5)]))

}

truth = t(sapply(as.character(pData$`tissue:ch1`),extractPercentages))/100
colnames(truth) = c('Liver','Brain','Lung')
rownames(truth) = colnames(exp)


### MEAN REFERENCES
mean.reference =
    data.frame(Liver =
                   as.numeric(rowMeans(exp[,colnames(exp) %in%
                                               rownames(truth)[1:3]])),
               Brain =
                   as.numeric(rowMeans(exp[,
                                           colnames(exp) %in%
                                               rownames(truth)[4:6]])),
               Lung =
                   as.numeric(rowMeans(exp[,
                                           colnames(exp) %in%
                                               rownames(truth)[7:9]])))
rownames(mean.reference) = rownames(exp)

getGenesbyCV = function(exp){

    require(ggplot2)
    means = rowMeans(exp)
    sds = apply(exp,1,sd)
    df = data.frame(Mean = means, SD = sds)
    plot = ggplot(data = df, aes(Mean,SD)) + geom_point() + geom_hline(yintercept = median(df$SD),color='red',linetype=2) +
        geom_vline(xintercept = median(df$Mean),color='red',linetype=2)
    index = which(means >= median(means) & sds > median(sds))
    return(list(Plot = plot,Index = as.numeric(index)))

}


for (s in c(200, 500, 800)){

    lll = list()
    for (i in 1:25){
        iii = getGenesbyCV(exp)$Index
        samp = sample(iii,s)
        o.mix = exp[samp,]
        t.mix = exp[iii,]
        t.prop = truth
        lll = rlist::list.append(lll,
                                 list(observed.mixed = o.mix,
                                      true.mixed.expression = t.mix,
                                      true.proportions = truth))}

    saveRDS(lll,
            paste0('K',s,'_GSE19830.RDS'))


}
