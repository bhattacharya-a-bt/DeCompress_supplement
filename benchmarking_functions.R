### functions used in generating pseudo-targeted panels

regBetweenRows <- function(shift,mat){
  require(vsn)
  expressions <- log(mat+shift)
  means <- rowMeans(expressions)
  sds <- rowSdDiffs(expressions)
  return(as.numeric(abs(coef(lm(sds~means))[2])))
}

mse.matrices = function(list.mat){
  
  require(DescTools)
  col.mse <- function(col,mat1,mat2){
    v1 = mat1[,col]
    v2 = mat2[,col]
    return(mean((v1-v2)^2))
  }
  
  mat1 = list.mat[[1]]
  mat2 = list.mat[[2]]
  return(mean(sapply(1:ncol(mat1),col.mse,mat1,mat2)))
  
}

get_num_components = function(q){
  sapply(q,function(list) ncol(list[[1]]))
}

sing.val.dist <- function(l){

	return(abs(sum(svd(l[[2]])$d)^.5 - sum(svd(l[[1]])$d)^.5))

}

mse.error <- function(list.mat,rand=F){

	mat1 = list.mat[[1]]
	mat2 = list.mat[[2]]

	require(hydroGOF)
	mse.row = function(i,mat1,mat2){
		return(mse(mat1[i,],mat2[i,]))
	}

	mse.all = function(mat1,mat2){
		return(mean(sapply(1:nrow(mat1),mse.row,mat1=mat1,mat2=mat2)))
	}

	if (rand == F){
		return(mse.all(mat1,mat2))
	}

	if (rand == T){

		require(gtools)
		p = permutations(n = ncol(mat1),r = ncol(mat1),v = 1:ncol(mat1),repeats.allowed=F)
		error = 100
		for (j in 1:nrow(p)){
			error = min(error,mse.all(mat1,mat2[,p[j,]]))
		}
		return(error)

	}

}