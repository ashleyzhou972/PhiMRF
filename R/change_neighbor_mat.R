change_mat<-function(sub_neighbor){
	#sub_neighbor should be a adjacency matrix with
	#1 indicating edge and 0 indicating no edge
	dim = dim(sub_neighbor)[1]
	rowSum_mat = matrix(rep(rowSums(sub_neighbor),dim), byrow=F, nrow=dim)
	colSum_mat = matrix(rep(colSums(sub_neighbor),dim), byrow=T, nrow=dim)
	sum_mat = rowSum_mat + colSum_mat
	sum_mat[sum_mat==0]<--Inf
	new_neighbor = sub_neighbor/sum_mat
	return(new_neighbor)
}


subset_mat<-function(datafolder, key_net, key_y, trans_already, genename_vector,
                     outputfolder, output_key_net, output_key_y){
  #The exisiting neighborhood should be located in the datafolder
  #with name [key_net]_neighbors_trans.rds if [trans_already]=T
  #or with name [key_net]_neighbors.rds if [trans_already]=F
  #y should be located in [datafolder]/y/
  #with name [key_y]_y.rds
  #colnames of [key_net]_neighbors{_trans}.rds should agree with y
  #[outputfolder] can be NULL, if NULL, not save
  #output subsetted neighbor in [outputfolder]
  #with name [output_key_net]_neighbors_trans.rds
  #output y in [outputfolder]/y/
  #with name [output_key_y]_y.rds
  if (trans_already){
    ori_neighbor = readRDS(paste0(datafolder, '/', key_net, "_neighbors_trans.rds"))
  } else{
    ori_neighbor = readRDS(paste0(datafolder, '/', key_net, "_neighbors.rds"))
  }
  if (class(ori_neighbor)!='dgCMatrix'){
    stop("Original neighbor must of the class dgCMatrix in the Matrix package\n")
  }
  names = colnames(ori_neighbor)
  #print(length(names))
  ori_y = readRDS(paste0(datafolder, '/y/', key_y, '_y.rds'))
  # TODO 20190917: if length(genename_vector)==1
  ind = match(genename_vector, names)
  #print(genename_vector[1:3])
  #print(names[1:3])
  ind = ind[!is.na(ind)]
  new_neighbor = ori_neighbor[ind,ind]
  if (length(ind)==1){
    # new_neighbor = methods::as(as.matrix(new_neighbor), "dgCMatrix")
    attributes(new_neighbor)$Dimnames[[1]]<-attributes(ori_neighbor)$Dimnames[[1]][ind]
    attributes(new_neighbor)$Dimnames[[2]]<-attributes(ori_neighbor)$Dimnames[[2]][ind]
    new_y = matrix(ori_y[ind, ], ncol = 2)
    colnames(new_y)<-colnames(ori_y)
  } else {
    new_y = ori_y[ind, ]
  }
  attributes(new_neighbor)$x<-rep(1, length(attributes(new_neighbor)$x))
  if (dim(new_neighbor)[1]>5000){
    stop("Subset neighbor too large, use preprocess\n")
  }
  new_neighbor_trans = change_mat(new_neighbor)
  rownames(new_y) = names[ind]
  if (!is.null(outputfolder) && !is.null(output_key_net) && !is.null(output_key_y)){
    saveRDS(new_neighbor_trans, paste0(outputfolder, "/", output_key_net, "_neighbors_trans.rds"))
    saveRDS(new_y, paste0(outputfolder, "/y/", output_key_y, "_y.rds"))
  }
  return(list(net = new_neighbor_trans, y = new_y))
}
