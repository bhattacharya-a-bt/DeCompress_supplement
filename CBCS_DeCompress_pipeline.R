### SAMPLE CODE FOR CBCS ANALYSIS - CBCS EXPRESSION AVAILABLE AT GSE148426

require(DeCompress)

### STEP 1: feature selection of reference expression matrix
inform.set = findInformSet(tcga,
                           method = 'variance',
                           n_genes = 1000,
                           n.types = 5,
                           scree = 'cumvar')

### STEP 2: train compression matrix for projection of
###         target feature space --> inform.set
cs.cbcs = trainCS(tcga[match(rownames(cbcs),rownames(tcga)),],
                  inform.set,
                  method = c('lasso',
                             'enet',
                             'ridge'),
                    par = T,
                    n.cores = 10,
                    lambda = .1)

### STEP 3: DeCompress the target panel
decompressed.cbcs = expandTarget(cbcs,
                                 cs.cbcs$compression_mat)

### STEP 4: run ensemble deconvolution
decompress.res = bestDeconvolution(decompressed.cbcs,
                                   n.types = 5,
                                   scree = 'cumvar',
                                   logTransform = F,
                                   known.props = NULL,
                                   methods = c('TOAST',
                                          'linseed',
                                          'celldistinguisher'))
