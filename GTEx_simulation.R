require(data.table)

### DOWNLOAD GTEX MEDIAN PROFILES
gtex_counts =
    fread('https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz')
gtex_counts$Ensembl =
    sapply(strsplit(gtex_counts$gene_id,'[.]'),
           function(vec) vec[1])


### FUNCTIONS TO GENERATE SIMS

### create a proportion matrix
runif.unit <- function (n, min = 0, max = 1){

    x <- runif(n = n, min = min, max = max)
    return(x/sum(x))

}

### generate proportion matrix, pure profiles, mixed expression reference,
### and target datasets for parameters: sample size (N),
### compartment profiles (pure), number of genes in the target (n.obs.genes),
### and multiplicative noise (noise.sd)
gen.mix.exp <- function (N, pure, n.obs.genes, noise.sd) {

    pure = as.matrix(pure)
    props = t(replicate(N,runif.unit(ncol(pure),0,1)))
    true.mixed = pure %*% t(props)
    colnames(true.mixed) = paste0('Samp',1:ncol(true.mixed))
    noise = abs(matrix( rnorm(nrow(true.mixed)*ncol(true.mixed),
                              mean=0,sd=noise.sd),
                        nrow(true.mixed),
                        ncol(true.mixed) ))
    true.mixed = true.mixed * noise
    s = sample(nrow(true.mixed),n.obs.genes)
    obs.mixed = true.mixed[s,]
    obs.pure = pure[s,]
    return(list(true.proportions = props,
                true.pure.profiles = pure,
                true.mixed.expression = true.mixed,
                observed.mixed = obs.mixed))

}

generate_all_sim <- function(pure, noise.sd, ng, N,cell.num){
    require(rlist)
    comp = list()

    for (i in 1:25){
        comp = list.append(comp,
                           gen.mix.exp(N = N,
                                       pure = as.matrix(pure),
                                       n.obs.genes = ng,
                                       noise.sd = noise.sd))
    }
    outName = paste0('NEW_SD',noise.sd,'_K',ng,'_',cell.num,'.RDS')
    saveRDS(comp,outName)
}

### pick the full pure matrices that simulate normal breast
full.pure.2 = as.data.frame(gtex_counts[,c('Breast - Mammary Tissue',
                                           'Cells - EBV-transformed lymphocytes')])
full.pure.3 = as.data.frame(gtex_counts[,c('Breast - Mammary Tissue',
                                           'Cells - Transformed fibroblasts',
                                           'Cells - EBV-transformed lymphocytes')])
full.pure.4 = as.data.frame(gtex_counts[,c('Breast - Mammary Tissue',
                                           'Cells - Transformed fibroblasts',
                                           'Cells - EBV-transformed lymphocytes',
                                           'Adipose - Subcutaneous')])
full.pure.5 = as.data.frame(gtex_counts[,c('Breast - Mammary Tissue',
                                           'Cells - Transformed fibroblasts',
                                           'Cells - EBV-transformed lymphocytes',
                                           'Adipose - Subcutaneous',
                                           'Whole Blood')])
rownames(full.pure.2) = rownames(full.pure.3) =
    rownames(full.pure.4) = rownames(full.pure.5) = gtex_counts$Ensembl

### pure matrix for dissimilar tissues
full.pure.4.bad = as.data.frame(gtex_counts[,c('Breast - Mammary Tissue',
                                               'Whole Blood',
                                               'Pancreas',
                                               'Pituitary')])



for (noise.sd in c(8,4)){
    for (K in c(200,500,800)){
        generate_all_sim(full.pure.2[rowSums(full.pure.2)!=0,],
                         noise.sd = noise.sd,N = 200,ng=K,
                         cell.num=2)
        generate_all_sim(full.pure.3[rowSums(full.pure.3)!=0,],
                         noise.sd = noise.sd,N = 200,ng=K,
                         cell.num=3)
        generate_all_sim(full.pure.4[rowSums(full.pure.4)!=0,],
                         noise.sd = noise.sd,N = 200,ng=K,
                         cell.num=4)
    }}

for (noise.sd in c(8,4)){
    for (K in c(200,500,800)){
        generate_all_sim(full.pure.4.bad[rowSums(full.pure.4.bad)!=0,],
                         noise.sd = noise.sd,N = 200,ng=K,cell.num=4)
        }
    }
