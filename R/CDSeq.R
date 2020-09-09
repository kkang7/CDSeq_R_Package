#' Complete deconvolution using sequencing data.
#' 
#' \code{CDSeq} takes bulk RNA-seq data as input and simultaneously returns estimates of both cell-type-specific gene expression profiles and sample-specific cell-type proportions (versioin 1.0.6).
#' 
#' @param bulk_data RNA-Seq read counts matrix. Columns represent samples and rows represent genes.
#' 
#' @param beta beta is a scalar or a vector of length G where G is the number of genes; default value for beta is 0.5; When beta=Null, CDSeq uses reference_gep to estimate beta.
#' 
#' @param alpha alpha is a scalar or a vector of length cell_type_number where cell_type_number is the number of cell type; default value for alpha is 5.
#' 
#' @param cell_type_number number of cell types. cell_type_number can be an integer or a vector of different integers. To estimate the number of cell types, please provide a vector for cell_type_number,
#' e.g. cell_type_number <- 2:30, then CDSeq will estimate the number of cell types.
#' 
#' @param mcmc_iterations number of iterations for the Gibbs sampler; default value is 700.
#' 
#' @param dilution_factor a scalar to dilute the read counts for speeding up; default value is 1. CDSeq will use bulk_data/dilution_factor. 
#' 
#' @param gene_subset_size number of genes randomly sampled for each block. Default is NULL.
#' 
#' @param block_number number of blocks. Each block contains gene_subset_size genes. Default is 1.
#' 
#' @param cpu_number number of cpu cores that can be used for parellel computing; 
#' Default is NULL and CDSeq will detect the available number of cores on the device and use number of all cores - 1 for parallel computing.
#' 
#' @param gene_length a vector of the effective lenght (gene length - read length + 1) of each gene; Default is NULL.
#' 
#' @param reference_gep a reference gene expression profile can be used to determine the cell type and/or estimate beta; Default is NULL.
#' 
#' @param print_progress_msg_to_file print progress message to a text file. Set 1 if need to print progress msg to a file and set 0 if no printing. Default is 0;
#' 
#' @importFrom stats cor
#' 
#' @examples 
#' result<-CDSeq(bulk_data =  mixtureGEP, cell_type_number = 6, mcmc_iterations = 10, 
#'        dilution_factor = 1, block_number = 1, gene_length = as.vector(gene_length), 
#'        reference_gep = refGEP, cpu_number=1, print_progress_msg_to_file=1)
#' 
#' result<-CDSeq(bulk_data =  mixtureGEP, cell_type_number=2:8, mcmc_iterations = 10, 
#'        dilution_factor = 1, block_number = 1, gene_length = as.vector(gene_length), 
#'        reference_gep = refGEP, cpu_number=2, print_progress_msg_to_file=0)
#' @export
#'  
#' @return CDSeq returns estimates of both cell-type-specific gene expression profiles and sample-specific cell-type proportions. CDSeq will also return estimated number of cell types.
#'  and the log posterior values for different number of cell types. 


#---------
# Feature
#---------
#1. Use reduce-recovery method and CPU parallel computing to speed up the deconvolution.
#2. Estimate hyperparameter for cell-type-specific GEPs (i.e. beta) using reference profile when cell_type_number is scalar.
#3. Estimate number of cell types when cell_type_number is a vector of integers.
#4. When block_number (number of partition on the bulk RNASeq data) is 1, whole bulk_data will be used. GEP is not from reduce-recovery.



#------------
# Parameters
#------------
#1. bulk_data:          RNA-Seq read counts matrix. Columns represent samples and rows represent genes. 
#2. beta:               beta is a scalar or a vector of length G where G is the number of genes; default value for beta is 0.5; When beta=Null, CDSeq uses reference_gep to estimate beta.
#3. alpha:              alpha is a scalar or a vector of length cell_type_number where cell_type_number is the number of cell type; default value for alpha is 5.
#4. cell_type_number:   number of cell types. cell_type_number can be an integer or a vector of different integers. To estimate the number of cell types, please provide a vector for cell_type_number, 
#                       e.g. cell_type_number <- 2:30, then CDSeq will estimate the number of cell types.
#5. mcmc_iterations:    number of iterations for the Gibbs sampler; default value is 700.
#6. dilution_factor:    a scalar to dilute the read counts for speeding up; default value is 1. CDSeq will use bulk_data/dilution_factor. 
#7. gene_subset_size:   number of genes randomly sampled for each block. Default is NULL.
#8. block_number:       number of genes randomly sampled for each block. Default is 1.
#9. cpu_number:         number of cpu cores that can be used for parellel computing; Default is NULL and CDSeq will detect the available number of cores on the device and use number of all cores - 1 for parallel computing.
#10.gene_length:        a vector of the effective lenght (gene length - read length + 1) of each gene; Default is NULL.
#11.reference_gep:      a reference gene expression profile can be used to determine the cell type and/or estimate beta; Default is NULL.




CDSeq <- function( bulk_data,  
                   beta = 0.5, 
                   alpha = 5, 
                   cell_type_number = NULL, 
                   mcmc_iterations = 700, 
                   dilution_factor = 1,
                   gene_subset_size = NULL,
                   block_number = 1, 
                   cpu_number = NULL,
                   gene_length = NULL,
                   reference_gep = NULL,
                   print_progress_msg_to_file = 0) {
  
  cat("\n============================================ CDSeq R package version 1.0.7 ==================================================\n")
  cat("|  Coders     : Kai Kang, David Huang                                                                                         |\n")
  cat("|  Reference 1: CDSeq: A novel complete deconvolution method for dissecting heterogeneous samples using gene expression data  |\n")
  cat("|  Reference 2: CDSeq: An R package for fast complete deconvolution using gene expression data                                |\n")
  cat("|  Maintainer : kangkai0714@gmail.com                                                                                         |\n")
  cat("===============================================================================================================================\n")
  
  start_time <- Sys.time()
  ########################
  # check input arguments
  ########################
  #check number of arguments
  if(nargs()>11){stop("CDSeq function requires at most 11 arguments.")}
  
  #check bulk_data
  if(is.null(bulk_data)){stop("bulk_data is NULL. Please provide valid bulk RNA-Seq counts data.")}
  if(!is.matrix(bulk_data) && !is.data.frame(bulk_data)){stop("The input bulk_data has to be either a matrix or dataframe.")}
  if(is.data.frame(bulk_data)){bulk_data <- as.matrix(bulk_data)}
  
  #check block_number
  if(is.null(block_number)){block_number<-1;warning("block_number is NULL. CDSeq assigned block_number to be 1.")}
  if(block_number < 1){stop("block_number has to be greater than or equal to 1.")}
  if(block_number%%1!=0){block_number<-ceiling(block_number);warning("block_number has to be an integer. CDSeq has rounded it up to an integer.")}
  if(block_number>1){cat("CDSeq is running in Reduce-Recover mode which breaks up the whole data into blocks. To run CDSeq on the whole data, please set block_number = 1.\n")
    }else{cat("CDSeq is running in non Reduce-Recover mode. To use Reduce-Recover mode, assign a value to block_number that is greater than 1.\n")}
  
  #check gene_subset_size and bulk_data
  if(is.null(gene_subset_size) & block_number>1){stop("gene_subset_size is null. Please assign a positive value to gene_subset_size, e.g. gene_subset_size = 1000, or set block_number to be 1")}
  if(!is.null(gene_subset_size)){
    if(!is.numeric(gene_subset_size)){stop("gene_subset_size has to be an integer.")}else if( gene_subset_size%%1!=0){ gene_subset_size <- round(gene_subset_size) }
    if(gene_subset_size>nrow(bulk_data)){stop("gene_subset_size is greater than nrows(bulk_data). Please assign a smaller value to gene_subset_size.")}
    if(gene_subset_size<nrow(bulk_data) & block_number==1){warning("gene_subset_size is smaller than nrows(bulk_data) and block_number is 1. This may result in inaccurate estimations. We suggest to increase block_number or set gene_subset_size <- nrows(bulk_data).")}
    if(nrow(bulk_data)<=1000 & block_number>1){message("We suggest to run the whole data instead of using Reduce-Recovery when the total number of genes in bulk_data is not that large.")}
  }
  
  #check beta
  if(is.null(beta)){
    if(is.null(reference_gep)){stop("Please provide a value for beta. beta can be a scalar or a vector of length equals to the number of genes. Alternatively, provide a reference_gep if beta is null.")}}
  else{
    if(!is.numeric(beta) | is.matrix(beta) | is.data.frame(beta) ){stop("beta has to be a real number or a vector of real numbers.")}
    #if(length(cell_type_number)>1 & is.null(beta)){stop("cell_type_number should be a scalar if beta is Null")}
    if(length(beta)>1 & block_number==1 & length(beta)!=nrow(bulk_data)){stop("Length of beta should be equal to total number of genes if beta is a vector")}
    if(length(beta)>1 & block_number>1 & !is.null(gene_subset_size)){ if(length(beta)!=gene_subset_size) stop("When block_number is greater than 1, beta should be either a scalar or a vector of length gene_subset_size.")}
    if(length(beta)==1){
      if(block_number>1){
        beta <- rep(beta,gene_subset_size)
      }else if(block_number==1){
        if(!is.null(gene_subset_size)){
          if(gene_subset_size<=nrow(bulk_data)) { beta <- rep(beta,gene_subset_size) }
        }else{
          beta <- rep(beta,nrow(bulk_data))
        }
      }else{stop("block_number has to be a positive integer.")}
      
    } # Gibbs sampler requires beta to be a vector for computation
  }
  
  #check alpha
  if(is.null(alpha)){stop(" Input alpha is missing. alpha should be a positive real number.")}
  if(length(alpha)>1){stop(" alpha should be a positive real number, not a vector.")}
  
  #Check cell_type_number
  if(is.null(cell_type_number)){stop(" cell_type_number, the number of cell types, is missing. cell_type_number can be a scalar value or a vector.")}
  if(sum(cell_type_number<2)>0){stop("cell_type_number, the number of cell types, has to be greater than 2.")}
  
  #Check mcmc_iterations
  if(is.null(mcmc_iterations)){stop(" mcmc_iterations, the number of MCMC iterations, is missing.")}
  if(mcmc_iterations<1){stop(" mcmc_iterations, the number of MCMC iterations, has to be greater than 1.")}
  
  #Check dilution_factor
  #if(is.null(dilution_factor) || dilution_factor==1){dilution_factor<-1; warning("dilution_factor is NOT used for speeding up ")}
  if(length(dilution_factor)>1){stop("dilution_factor has to be an positive integer, NOT a vector.")}
  if(dilution_factor<0){stop("dilution_factor has to be positive.")}

  #Check Gene length
  #if not given, gl = 0, estimate RNA proportions, not cell proportions 
  if(is.null(gene_length)){
    warning("gene_length is NOT provided. CDSeq will estiamte read rate not gene rate. Please provide gene length if you are interested in GEP estimation.")
    gl<-0
  }else{
    if(is.matrix(gene_length)){stop("gene_length should be a vector, NOT a matrix.")}
    if(length(gene_length)!=nrow(bulk_data)){stop("length(gene_length) should be equal to nrow(bulk_data)")}
    if(any(gene_length < 1)){stop("gene_length denotes effective legnth of genes and should be a vector of positive values that are greater than 1.")}
    #if(floor(sum(gene_length))!=sum(gene_length)){gene_length = round(gene_length);warning("The provided gene_length seems to contain non-integers. It has been round to integer.")}
    if(any(gene_length%%1!=0)){gene_length = round(gene_length);warning("The provided gene_length seems to contain non-integers. CDSeq has rounded up to integers.")}
    gl <- 1
  }

  #Check reference_gep
  #if NOT given, then rpkm = ref = 0, and cell type assignment and RPKM normalization will NOT be performed
  if(is.null(reference_gep)){
    rpkm<-0;ref<-0;estimate_beta<-0
    warning("Reference gene expression profile is missing. Cell type identification and RPKM normalization will NOT be performed by CDSeq. Users can identify CDSeq-identified cell-types using marker genes or reference gene expression profiles.")
  }else{ 
    if(nrow(reference_gep)!=nrow(bulk_data)){stop("nrow(reference_gep) should be equal to nrow(bulk_data)")}
    if(!is.matrix(reference_gep) & !is.data.frame(reference_gep)){stop("The reference GEPs should be either a matrix or data frame")}
    rpkm<-1;ref<-1
    #cls<-sum(reference_gep[,1]) #check if reference profile is raw read counts data
    if(any(reference_gep%%1!=0)){ #check if reference profile is raw read counts data
      if(is.null(beta)){
        stop("reference_gep should be read counts data for estimating beta. Please provide read counts reference_gep or provide a value (either a scalar or vector of length nrow(bulk_data)) for beta.")
      }else{warning("reference_gep is NOT read counts data, RPKM normalization is NOT performed.");rpkm<-0;estimate_beta<-0}
    }else{
      if(is.null(beta) & block_number == 1 ){estimate_beta<-1}else{estimate_beta<-0}
    }
    if(gl!=1){
      rpkm<-0
      warning("Gene length is NOT provided, RPKM normalization is NOT performed.")
     }
  }
  
  # check the cell_type_number and ncol(refGEP), 
  if(!is.null(reference_gep)){
    if(max(cell_type_number)>ncol(reference_gep)){
      rpkm <- 0
      ref <- 0
    }
  }
  
  
  ###################################
  # Save sample names and gene names
  ###################################
  
  gene_names<-row.names(bulk_data)
  sample_names<-colnames(bulk_data)
  cell_types <- NULL  

  #Check if gene name is consistent in bulk_data and reference_gep
  if(!is.null(reference_gep)){
    if(!is.null(gene_names) && !is.null(row.names(reference_gep))){
      if( any(gene_names!=row.names(reference_gep)) ){
        warning("gene names in bulk_data may be different from the gene names in reference profile. Please make sure. CDSeq will use gene names from bulk_data.")
      }
    }
    cell_types<-colnames(reference_gep)
    reference_gep<-as.matrix(reference_gep)
  }
  
  if(is.null(gene_names)){gene_names<-paste("gene",1:nrow(bulk_data),sep = "_")}
  if(is.null(sample_names)){sample_names<-paste("sample",1:ncol(bulk_data),sep = "_")}

  
  ###############
  # Data Dilution
  ###############
  # check if bulk_data is read counts data
  if(any(bulk_data%%1!=0)){ warning(" bulk_data is NOT read count data. Please provide read counts data if possible for potentially better estimations.");rawcount<-0 }else{rawcount<-1}
  bulk_data_diluted<-ceiling(bulk_data/dilution_factor) # total read counts of bulk_data is reduced by data dilution 

  
  #############################
  # Break up bulk_data into blocks 
  #############################
  bulk_data_blocks <- c() # is this a good way of initialize the list?
  refGEPlist <- c()
  if(block_number==1){
    bulk_data_blocks[[1]] <- bulk_data_diluted
    if(ref==1){refGEPlist[[1]] <- reference_gep}
  }else{
    #total number of genes
    totalg <- nrow(bulk_data_diluted)
    
    #save the gene index for each block
    groupmatrix <- matrix(0,nrow = block_number ,ncol = gene_subset_size )
    for (i in 1:block_number) {groupmatrix[i,] <- sample(totalg,gene_subset_size,replace = TRUE)}
    
    #GEP for each block
    for(i in 1:block_number){bulk_data_blocks[[i]] <-  bulk_data_diluted[groupmatrix[i,],]}
    if(estimate_beta==1){for(i in 1:block_number){ refGEPlist[[i]] <-  reference_gep[groupmatrix[i,],]}}
  }
  
  #####################################################
  # Estimate beta from ref(only use when cell_type_number is a scalar)
  #####################################################
  beta_est <- NULL
  #if(is.null(beta))
  if(estimate_beta==1){
    beta_est <- matrix(0,nrow = block_number,ncol = gene_subset_size)
    for(i in 1:block_number){
      fit <- dirmult(t(refGEPlist[[i]]))
      beta_est[i,] <- fit$gamma 
    }
  }
  
  ###########################
  # parallel computing setup
  ###########################
  if(is.null(cpu_number)){
    registerDoParallel(detectCores()-1)
    message("cpu_number is not provided. CDSeq uses detectCores()-1 number of cpu cores for parallel computing.\n")
  }else{
    if(!is.numeric(cpu_number) || !(is.atomic(cpu_number) && length(cpu_number)==1) || is.matrix(cpu_number)){stop("cpu_number has to be an integer that is greater than or equal to 1.")}
    if(cpu_number > detectCores()-1 | cpu_number < 1){stop("cpu_number is invalid. To detect how many cores are available, use detectCores().")}
    if(cpu_number%%1!=0){
      cpu_number <- ceiling(cpu_number)
      warning("cpu_number is not an integer. CDSeq uses ceiling(cpu_integer).")
    }
    registerDoParallel(cpu_number)
  }
  
  celltype_assignment <- NULL
  cellTypeAssignSplit <- NULL
  # keep track of the worker process ID for parallel progress printing
  processIDs <- matrix(0,length(cell_type_number),block_number) 
  #=====================================================================================================
  # if cell_type_number is a scalar
  #=====================================================================================================
  
  if(length(cell_type_number)==1){
    CDSeq_tmp_log <- tempfile(pattern = "CDSeq_tmp_log_",fileext = ".txt")
    #Check if rpkm normalization should be performed
    # if the reference profile is not read counts data, then RPKM will not be performed
    # if the number of cell types in reference profile is less than cell_type_number, then only part of CDSeq-identified cell types will be associated with reference profiles
    if(ref==1 && cell_type_number>ncol(reference_gep)){
      warning("The number of cell types cell_type_number is greater than the number of cell types in reference profile, so not all CDSeq-identified cell types can be associated with reference cell types.")
      rpkm<-0
    }
    
    if(block_number>1){
      if(print_progress_msg_to_file==1){cat(sprintf("CDSeq is running in parallel, it may take some time...\n(Progress is being printed to %s. You may delete the file if will not use it.)\n",CDSeq_tmp_log))}
      else{cat("CDSeq is running in parallel, it may take some time...\n")}
      printout <- 0}
    else{printout <- 1}
    
    #Gibbs sampler and store its running time
    gibbsRunningTime<-system.time(
      {
        
        allresult <- foreach(i=1:block_number, .inorder = FALSE) %dopar% {
          
          if(!is.null(beta_est)){beta <- beta_est[i,]}

          result <- gibbsSampler(alpha,beta,bulk_data_blocks[[i]],cell_type_number,mcmc_iterations, printout, Sys.getpid(), i, CDSeq_tmp_log, print_progress_msg_to_file)
          
          #outputs are two vectors. estGEP_vec is gene x cell type; estSSp_vec is sample x cell type
          estGEP_vec<-result$csGEP_vec
          estSSp_vec<-result$SSP_vec
          
          #vector to matrix
          samplesize<-ncol(bulk_data)
          estGEP_mat<-t(matrix(estGEP_vec,nrow = cell_type_number,ncol = nrow(bulk_data_blocks[[i]]))) # gene by cell_type
          estSSp_mat<-matrix(estSSp_vec,nrow = cell_type_number,ncol = samplesize) #cell_type by sample_size
          
          #estimated proportions and GEPs using bulk_data blocks
          estProp<-t(t(estSSp_mat+alpha)/colSums(estSSp_mat+alpha))  #cell_type by sample_size
          estGEP_read<-t(t(estGEP_mat+beta)/colSums(estGEP_mat+beta))#gene by cell_type
          
          if(block_number==1){
            if(gl==1){estGEP<-read2gene(estGEP_read,gene_length)}else{estGEP<-estGEP_read} # read to gene
            if(ref==1){
              if(rawcount==1){
                corr_GEP<-cor(estGEP_read,refGEPlist[[i]])
              }else{
                corr_GEP<-cor(estGEP,refGEPlist[[i]])
                if(gl!=1){warning("Gene length is NOT provided and the reference GEP is NOT read counts data, the cell type association may be inaccurate.")}
              }
              hungarian_result<-hungarian_Rcpp(1-corr_GEP)
              celltype_assignment<-hungarian_result$cost_assignment+1
              
              if(rawcount==1){estProp <- RNA2Cell(colSums(refGEPlist[[i]][,celltype_assignment]),estProp)}
              if(rpkm==1){estGEP <- gene2rpkm(estGEP,gene_length,refGEPlist[[i]][,celltype_assignment])}
            }
            return(list(estProp = estProp, estGEP = estGEP, celltype_assignment = celltype_assignment,cellTypeAssignSplit = result$cellTypeAssignSplit, processID = Sys.getpid()))
          }else{
            #return(estProp)
            return(list(estProp = estProp, processID = Sys.getpid()))
          }
        }#foreach end
      }
    )#gibbs time
    
    #stop parallel
    stopImplicitCluster()
    
    if(block_number>1){
      #reorder the estProp for each block
      estproplist <- c()
      for(i in 1:block_number){
        corr <- cor(t(allresult[[1]]$estProp), t(allresult[[i]]$estProp))
        hungarian_result<-hungarian_Rcpp(1-corr)
        ordervalue <- hungarian_result$cost_assignment+1
        estproplist[[i]] <- as.matrix(allresult[[i]]$estProp)[ordervalue,]
        processIDs[i] <-allresult[[i]]$processID
      }
      
      #get the average of the prop
      nsample <- ncol(allresult[[1]]$estProp)
      sumestprop <- matrix(0,nrow = cell_type_number,ncol = nsample)
      for(i in 1:block_number){sumestprop <- estproplist[[i]]+ sumestprop}
      averestprop <- sumestprop/block_number
      averestprop <- sweep(averestprop,MARGIN=2,FUN="/",STATS=colSums(averestprop))
      
      #reduce recover method
      bulk_data_norm<-t(t(bulk_data)/colSums(bulk_data))
      estProp_right_inv<-t(averestprop)%*%ginv(crossprod(t(averestprop)))
      estGEP_read<-abs(bulk_data_norm%*%estProp_right_inv)
      estGEP_read<-t(t(estGEP_read)/colSums(estGEP_read)) # renormalize
      estProp<-averestprop
      if(gl==1){estGEP<-read2gene(estGEP_read,gene_length)}else{estGEP<-estGEP_read} # read to gene
      
      # cell type association
      if(ref==1){
        if(rawcount==1){
          #corr_GEP<-cor(estGEP_read,refGEPlist[[i]])
          corr_GEP<-cor(estGEP_read,reference_gep)
        }else{
          if(gl!=1){corr_GEP<-cor(estGEP,reference_gep);warning("Gene length is NOT provided and the reference GEP is NOT read counts data, the cell type association may be inaccurate.")}
        }
        hungarian_result<-hungarian_Rcpp(1-corr_GEP)
        celltype_assignment<-hungarian_result$cost_assignment+1
        
        if(rawcount==1){estProp <- RNA2Cell(colSums(reference_gep[,celltype_assignment]),estProp)}
        if(rpkm==1){estGEP <- gene2rpkm(estGEP,gene_length,reference_gep[,celltype_assignment])}
      }
    }else{#for block_number=1
      estProp <- allresult[[1]]$estProp
      estGEP <- allresult[[1]]$estGEP
      celltype_assignment<-allresult[[1]]$celltype_assignment
      cellTypeAssignSplit <- allresult[[1]]$cellTypeAssignSplit
      processIDs[1] <-allresult[[1]]$processID
    } 
    # keep all the parameters
    parameters <- list( beta = beta, 
                        alpha = alpha, 
                        cell_type_number = cell_type_number, 
                        mcmc_iterations = mcmc_iterations, 
                        dilution_factor = dilution_factor,
                        gene_subset_size = gene_subset_size,
                        block_number = block_number, 
                        cpu_number = cpu_number,
                        gene_length = gene_length,
                        reference_gep = reference_gep,
                        print_progress_msg_to_file = print_progress_msg_to_file)
    #Final output
    CDSeq_result<-list(estProp=estProp,estGEP=estGEP,gibbsRunningTime = gibbsRunningTime, cell_type_assignment = celltype_assignment, cellTypeAssignSplit = cellTypeAssignSplit,processIDs = processIDs, parameters = parameters)
    
    if(ref==0){
      cell_types<-paste("unknown_cell_type",1:cell_type_number,sep = "_")
      rownames(CDSeq_result$estGEP)<-gene_names
      colnames(CDSeq_result$estGEP)<-cell_types
      
      rownames(CDSeq_result$estProp)<-cell_types
      colnames(CDSeq_result$estProp)<-sample_names
      
      dimnames(CDSeq_result$cellTypeAssignSplit)[[1]] <- gene_names
      dimnames(CDSeq_result$cellTypeAssignSplit)[[2]] <- sample_names
      dimnames(CDSeq_result$cellTypeAssignSplit)[[3]] <- cell_types
      
    }else{
      refcol<-ncol(reference_gep)
      if(is.null(cell_types)){cell_types<-paste("ref_cell_type",1:refcol,sep = "_")}
      if(refcol>=cell_type_number){cell_types_tmp<-cell_types[celltype_assignment]}
      if(refcol<cell_type_number){
        cell_types_tmp<-paste("unknow_cell_type",1:cell_type_number)
        cell_types_idx<-which(celltype_assignment!=0)
        cell_types_idx_2<-which(celltype_assignment==0)
        cell_types_tmp[cell_types_idx]<-cell_types[celltype_assignment[cell_types_idx]]
        cell_types_tmp[cell_types_idx_2]<-paste("unknow_cell_type",1:length(cell_types_idx_2))
      }
      
      colnames(CDSeq_result$estGEP)<-cell_types_tmp
      rownames(CDSeq_result$estGEP)<-gene_names
      
      rownames(CDSeq_result$estProp)<-cell_types_tmp
      colnames(CDSeq_result$estProp)<-sample_names
    }
    end_time <- Sys.time()
    CDSeq_running_time <- as.numeric( end_time - start_time, units = "hours")
    cat(sprintf("CDSeq completed successfully using %.4f hours\n",CDSeq_running_time))
    return(CDSeq_result)
  }
  
  #=====================================================================================================
  # if cell_type_number is a vector
  #=====================================================================================================
  if(length(cell_type_number)>1){
    #two parallel loop
    printout=0
    CDSeq_tmp_log <- tempfile(pattern = "CDSeq_tmp_log_",fileext = ".txt")
    if(print_progress_msg_to_file==1){
      cat(sprintf("CDSeq is running in parallel, it may take some time...\n(Progress is being printed to %s. You may delete the file if will not use it.)\n",CDSeq_tmp_log))
    }
    #cat("CDSeq is running in parallel, it may take some time...\n(Progress is being printed to CDSeq_logfile.txt in working directory. You may delete the file if will not use it.)\n")
    est_all<-foreach(j=1:length(cell_type_number), .combine = 'rbind') %:% 
      foreach(i=1:block_number, .combine = 'c') %dopar% {
      if(!is.null(beta_est)){beta <- beta_est[i,]}
      processIDs[j,i] <- Sys.getpid()
      result <- gibbsSampler(alpha,beta,bulk_data_blocks[[i]],cell_type_number[j],mcmc_iterations, printout, processIDs[j,i], i, CDSeq_tmp_log, print_progress_msg_to_file)

      #output two vectors. estGEP_vec is gene x cell type; estSSp_vec is sample x cell type
      estGEP_vec<-result$csGEP_vec
      estSSp_vec<-result$SSP_vec
      cellTypeAssignSplit <- result$cellTypeAssignSplit
      
      #vector to matrix
      samplesize<-ncol(bulk_data)
      estGEP_mat<-t(matrix(estGEP_vec,nrow = cell_type_number[j],ncol = nrow(bulk_data_blocks[[i]]))) # gene by cell_type
      estSSp_mat<-matrix(estSSp_vec,nrow = cell_type_number[j],ncol = samplesize) #  cell_type by sample_size
      
      #estimated proportions and GEPs on reduced data set
      estProp<-t(t(estSSp_mat+alpha)/colSums(estSSp_mat+alpha))#  cell_type by sample_size
      estGEP_read<-t(t(estGEP_mat+beta)/colSums(estGEP_mat+beta))# gene by cell_type
      
      #calculate logposterior 
      lgpst<-logpost(estProp, estGEP_read, bulk_data_blocks[[i]], alpha ,beta)
      
      if ( block_number == 1){
        if ( gl==1 ){estGEP<-read2gene(estGEP_read,gene_length)}else{estGEP<-estGEP_read} # read to gene
        if ( ref==1 ){
          if ( rawcount==1 & cell_type_number[j]<=ncol(refGEPlist[[i]])){
            corr_GEP<-cor(estGEP_read,refGEPlist[[i]])
          }else{
            if(gl!=1){warning("Gene length is NOT provided and the reference GEP is NOT read counts data, the cell type association may be inaccurate.")}
            corr_GEP<-cor(estGEP,refGEPlist[[i]])
          }
          hungarian_result<-hungarian_Rcpp(1-corr_GEP)
          celltype_assignment<-hungarian_result$cost_assignment+1
          if(rawcount==1 & nrow(estProp)<=ncol(refGEPlist[[i]])){estProp <- RNA2Cell(colSums(refGEPlist[[i]][,celltype_assignment]),estProp)}
          if(rpkm==1){estGEP <- gene2rpkm(estGEP,gene_length,refGEPlist[[i]][,celltype_assignment])}
        }
      }# save cellTypeAssignSplit only when block number is 1, since for block number greater than 1, cellTypeAssignSplit info is not useful. 
      if(block_number>1){return(list(estProp,lgpst))}else{return(list(estProp,estGEP,lgpst,celltype_assignment, cellTypeAssignSplit))}
    }#foreach loop end
    
    #stop parallel
    stopImplicitCluster()
    #use this to save the result for all cell_type_number value
    CDseq_all <- c()
    
    #when block_number > 1 or = 1
    if(block_number>1){
      estPropall <- est_all[,seq(1,2*block_number,by=2)]
      lgpstall <- est_all[,seq(2,2*block_number,by=2)]
    
      for(j in 1:length(cell_type_number)){
        #reorder the estProp
        estproplist <- c()
        for(i in 1:block_number){
          corr <- cor(t(estPropall[j,1][[1]]), t(estPropall[j,i][[1]]))
          hungarian_result<-hungarian_Rcpp(1-corr)
          ordervalue <- hungarian_result$cost_assignment+1
          estproplist[[i]] <- as.matrix(estPropall[j,i][[1]])[ordervalue,]
        }
        
        #get the average for each block
        nsample <- ncol(estPropall[j,1][[1]])
        sumestprop <- matrix(0,nrow = cell_type_number[j],ncol = nsample)
        for(i in 1:block_number){sumestprop <- estproplist[[i]]+ sumestprop}
        averestprop <- sumestprop/block_number
        averestprop <- sweep(averestprop,MARGIN=2,FUN="/",STATS=colSums(averestprop))
        
        #reduce recover
        bulk_data_norm<-t(t(bulk_data)/colSums(bulk_data))
        estProp_right_inv<-t(averestprop)%*%ginv(crossprod(t(averestprop)))
        estGEP_read<-abs(bulk_data_norm%*%estProp_right_inv)
        estGEP_read<-t(t(estGEP_read)/colSums(estGEP_read)) # renormalize
        
        # using the gene length information us convert read to gene
        if(gl==1){estGEP<-read2gene(estGEP_read,gene_length)}else{estGEP<-estGEP_read} # read to gene
        
        mlgpst <- mean(unlist(lgpstall[j,]))
        
        # add colnames and rownames
        cell_types<-paste("CDSeq_estimated_cell_type",1:ncol(estGEP),sep = "_")
        rownames(estGEP)<-gene_names
        colnames(estGEP)<-cell_types
        rownames(averestprop)<-cell_types
        colnames(averestprop)<-sample_names
        
        CDseq_all[[j]]<-list(estProp=averestprop,estGEP=estGEP,lgpst=mlgpst)
        
      }
    }else{#for block_number=1 
      estPropall <- est_all[,1]
      estGEPall <- est_all[,2]
      lgpstall <- est_all[,3]
      celltype_assignment_all<-est_all[,4]
      cellTypeAssignSplit_all <- est_all[,5]
      for(j in 1:length(cell_type_number)){
        averestprop <- estPropall[[j]]
        estGEP <- estGEPall[[j]]
        mlgpst <- lgpstall[[j]]
        celltype_assignment<-celltype_assignment_all[[j]]
        cellTypeAssignSplit <- cellTypeAssignSplit_all[[j]]
        
        # add colnames and rownames
        cell_types<-paste("CDSeq_estimated_cell_type",1:ncol(estGEP),sep = "_")
        rownames(estGEP)<-gene_names
        colnames(estGEP)<-cell_types
        rownames(averestprop)<-cell_types
        colnames(averestprop)<-sample_names
        
        dimnames(cellTypeAssignSplit)[[1]] <- gene_names
        dimnames(cellTypeAssignSplit)[[2]] <- sample_names
        dimnames(cellTypeAssignSplit)[[3]] <- cell_types
        
        CDseq_all[[j]]<-list(estProp=averestprop,estGEP=estGEP,lgpst=mlgpst, cell_type_assignment = celltype_assignment, cellTypeAssignSplit = cellTypeAssignSplit)
      }
    }
  
  
  #get the max of mlgpst and it's estProp and estGEP
  alllgpst <- c()
  for(i in 1:length(cell_type_number)){alllgpst[i] <- CDseq_all[[i]]$lgpst}
  maxTindex <- which(alllgpst==max(alllgpst))
  maxT <- cell_type_number[maxTindex]
  
  maxaverestprop <- CDseq_all[[maxTindex]]$estProp
  maxGEP <- CDseq_all[[maxTindex]]$estGEP
  maxlgpst <- max(alllgpst)
  
  maxcellTypeAssignSplit <- CDseq_all[[maxTindex]]$cellTypeAssignSplit
  
  # cell type association
  if(ref==1 & block_number>1){
    if(rawcount==1){
      corr_GEP<-cor(maxGEP,reference_gep)
    }else{
      if(gl!=1){corr_GEP<-cor(maxGEP,reference_gep);warning("Gene length is NOT provided and the reference GEP is NOT read counts data, the cell type association may be inaccurate.")}
    }
    hungarian_result<-hungarian_Rcpp(1-corr_GEP)
    celltype_assignment<-hungarian_result$cost_assignment+1
    
    if(rawcount==1 & nrow(maxaverestprop)<=ncol(reference_gep)){maxaverestprop <- RNA2Cell(colSums(reference_gep[,celltype_assignment]),maxaverestprop)}
    if(rpkm==1 & ncol(maxGEP)<=ncol(reference_gep)){maxGEP <- gene2rpkm(maxGEP,gene_length,reference_gep[,celltype_assignment])}
  }
  if(ref==1 & block_number==1){celltype_assignment<-CDseq_all[[maxTindex]]$cell_type_assignment}
  # keep all the parameters
  parameters <- list( beta = beta, 
                      alpha = alpha, 
                      cell_type_number = cell_type_number, 
                      mcmc_iterations = mcmc_iterations, 
                      dilution_factor = dilution_factor,
                      gene_subset_size = gene_subset_size,
                      block_number = block_number, 
                      cpu_number = cpu_number,
                      gene_length = gene_length,
                      reference_gep = reference_gep,
                      print_progress_msg_to_file = print_progress_msg_to_file)
  #Final output
  CDSeq_result_max<-list(estProp=maxaverestprop, estGEP=maxGEP, cell_type_assignment = celltype_assignment,cellTypeAssignSplit = maxcellTypeAssignSplit, lgpst=maxlgpst, estT=maxT, est_all = CDseq_all, parameters = parameters)
  
  if(ref==0){
    cell_types<-paste("unknown_cell_type",1:CDSeq_result_max$estT,sep = "_")
    rownames(CDSeq_result_max$estGEP)<-gene_names
    colnames(CDSeq_result_max$estGEP)<-cell_types
    
    rownames(CDSeq_result_max$estProp)<-cell_types
    colnames(CDSeq_result_max$estProp)<-sample_names
  }else{
    if(is.null(celltype_assignment)){cat("\n no cell type assignment\n")}
    refcol<-ncol(reference_gep)
    if(is.null(cell_types)){cell_types<-paste("ref_cell_type",1:refcol,sep = "_")}
    if(refcol>=CDSeq_result_max$estT){cell_types_tmp<-cell_types[celltype_assignment]}
    if(refcol<CDSeq_result_max$estT){
      cell_types_tmp<-paste("unknow_cell_type",1:CDSeq_result_max$estT)
      cell_types_idx<-which(celltype_assignment!=0)
      cell_types_idx_2<-which(celltype_assignment==0)
      cell_types_tmp[cell_types_idx]<-cell_types[celltype_assignment[cell_types_idx]]
      cell_types_tmp[cell_types_idx_2]<-paste("unknow_cell_type",1:length(cell_types_idx_2))
    }
    
    colnames(CDSeq_result_max$estGEP)<-cell_types_tmp
    rownames(CDSeq_result_max$estGEP)<-gene_names
    
    rownames(CDSeq_result_max$estProp)<-cell_types_tmp
    colnames(CDSeq_result_max$estProp)<-sample_names
  }
  
  end_time <- Sys.time()
  CDSeq_running_time <- as.numeric( end_time - start_time, units = "hours")
  #if(CDSeq_running_time<1){CDSeq_running_time < - as.numeric( end_time - start_time, units = "mins")}
  cat(sprintf("CDSeq completed successfully using %.4f hours\n",CDSeq_running_time))
  return(CDSeq_result_max)
  }#end for cell_type_number is a vector
  
}#end function

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  