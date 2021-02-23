#' \code{cellTypeAssignMarkerGenes} assigns CDSeq-identified cell types using user-provided marker gene list and plots heatmap.
#' @param cell_gep gene expression profile matrix with G rows (genes) and M columns (cell types).
#' @param marker_gene_list a G (genes) by C (cell types with known identities) matrix or dataframe that contains the marker genes for each cell type. Column names must be CellType and GeneName.
#' @param threshold a numeric value that provides the threshold of whether a known cell type in the marker gene list can be identified. 
#' @param fig_path the location where the heatmap figure is saved. 
#' @param rowlabelsize row label size
#' @param collabelsize column label size
#' @param margins a vector of length 2 indicates row and column label margins
#' @param fig_width figure width for pdf figure
#' @param fig_height figure height for pdf figure
#' @param keysize color key size for heatmap
#' @param srtcol column label angle
#' @param keypar color key layout
#' @param heatmap_name the name of heatmap figure of one-to-one assignment.
#' @param heatmap_name_fuzzy_assign the name of heatmap figure of fuzzy assignment. 
#' @param verbose if TRUE, some information will be printed.
#' @importFrom grDevices pdf dev.off 
#' @importFrom gplots heatmap.2 bluered
#' @importFrom clue solve_LSAP
#' @importFrom Biobase rowMax
#' @return cellTypeAssignMarkerGenes returns a list containing: 
#' GEP_markerSum (a A by B matrix where A is nrow(marker_gene_list), B is ncol(cell_gep)), 
#' 
#' GEP_markerSum_zscore (row-wise z score of GEP_markerSum), 
#' 
#' GEP_matched is cell_gep(,cell_type_idx),
#' 
#' cell_type_idx (column index of cell_gep that are considered matching with cell types in marker_gene_list), 
#' 
#' cell_type_matched stores the cell types in marker_gene_list that are considered to be matched with cell_gep, 
#' 
#' GEP_markerSum_zscore_matched contains only the rows of GEP_markerSum_zscore that are considered to be matched with some cell types in cell_gep. GEP_markerSum_zscore_matched and GEP_markerSum_zscore have same columns. 
#' 
#' cell_type_matched_fuzzy is a zero-one matrix that has the same size as GEP_markerSum_zscore_matched. If (i,j) element is one, means ith cell type in marker_gene_list is assigned to jth element in cell_gep.

cellTypeAssignMarkerGenes <- function(cell_gep = NULL,
                                     marker_gene_list = NULL,
                                     threshold = 2,
                                     fig_path = getwd(),
                                     rowlabelsize = 1,
                                     collabelsize = 1,
                                     margins = c(3,0),
                                     fig_width = 100,
                                     fig_height = 100,
                                     keysize = 1,
                                     srtcol = 45,
                                     keypar = c(3.5,0,3,0),
                                     heatmap_name = "cellTypeAssign_heatmap.pdf",# may remove .pdf extension  heatmap_name_twoway = NULL
                                     heatmap_name_fuzzy_assign = "cellTypeAssign_heatmap_fuzzy.pdf",
                                     verbose=FALSE){
  if(verbose){
    cat("-----------------in cell_type_assign_heatmap-----------------\n")
    cat("ncol(cell_gep) = ", ncol(cell_gep)," nrow = (cell_gep) = ", nrow(cell_gep),"\n")
  }
  
  # check input data
  if(is.null(rownames(cell_gep))){stop("cell_gep has no row name. row name should be gene names.")}
  
  # sum over marker gene expressions of cell_gep
  #colnames(marker_gene_list) <- c('CellType','GeneName')
  if(is.null(colnames(marker_gene_list))){stop("marker_gene_list needs to have column names: CellType and GeneName")}
  if(length(setdiff(colnames(marker_gene_list), c("CellType","GeneName")))!=0){stop("marker_gene_list column names have to be: CellType and GeneName")}
  
  celltypes <- unique(marker_gene_list$CellType)
  GEP_markerSum <- matrix(0, nrow = length(celltypes), ncol = ncol(cell_gep))
  rownames(GEP_markerSum) <- celltypes
  if(!is.null(colnames(cell_gep))){
    colnames(GEP_markerSum) <- colnames(cell_gep)
  }
  
  GEP_matched <- cell_gep
  cell_type_idx <- 1:ncol(cell_gep)
  cell_type_matched <- NULL
  GEP_markerSum_zscore <- NULL
  GEP_markerSum_zscore_matched <- NULL
  cell_type_matched_fuzzy <- NULL
  cell_type_assignment_fuzzy <- NULL
  cell_type_assignment <- NULL
  
  for (i in 1:length(celltypes)) {
    celltype_idx <- which(marker_gene_list$CellType %in% celltypes[i])
    if(length(celltype_idx)==0){stop("cell type ---",celltypes[i],"--- is missing. Something wrong with marker_gene_list. ")}
    markers <- marker_gene_list$GeneName[celltype_idx]
    markers_idx <- which(rownames(cell_gep) %in% markers)
    if(is.null(markers_idx)){
      warning("Markers for cell type --", celltypes[i] ,"-- not found in the cell_gep genes.")
    }else if(length(markers_idx)>1){
      GEP_markerSum[i,] <- colSums(cell_gep[markers_idx,])
    }else if(length(markers_idx)==1){
      GEP_markerSum[i,] <- cell_gep[markers_idx,]
    }else if(length(markers_idx)==0){
      warning("Markers for cell type --", celltypes[i] ,"-- not found in the cell_gep genes.")
    }
  }
  
  # compute zscore of the cell_gep
  
  # this computes the z score of the rows of GEP_markerSum
  GEP_markerSum_zscore <- t(scale(t(GEP_markerSum))) 
  GEP_markerSum_zscore[is.nan(GEP_markerSum_zscore)] <- 0
  
  # this chooses the rows whose max values are greater than threshold. 
  # this is saying if the max value of that row is greater than threshold
  # then we consider there is a match of that cell type
  cell_type_matched_idx <- which(Biobase::rowMax(GEP_markerSum_zscore) > threshold) 
  if(verbose){
    cat("length(cell_type_matched_idx) = ", length(cell_type_matched_idx), "threshold = ",threshold,"\n")
  }
  if(length(cell_type_matched_idx)==0){stop("no cell type assigned. all z scores are smaller than the threshold. Try lower the threshold value.")}
  cell_type_matched <- celltypes[cell_type_matched_idx]
  GEP_markerSum_zscore_matched <- GEP_markerSum_zscore[cell_type_matched_idx,]
  
  
  if(nrow(GEP_markerSum_zscore_matched) > ncol(GEP_markerSum_zscore_matched)){
    # if known cell types are greater than the cell types need to be assigned, then do transpose and do Munkres.
    cell_assign <- solve_LSAP(exp(t(GEP_markerSum_zscore_matched)),maximum = TRUE)
    marker_idx <- cell_assign[1:length(cell_assign)]
    
    cell_type_assignment <- data.frame(cbind(cell.type.in.marker.list = rownames(GEP_markerSum_zscore_matched[marker_idx,]), 
                                             cell.type.in.cell.gep = colnames(GEP_markerSum_zscore_matched[marker_idx,])))
    
    filename <- paste0(fig_path,heatmap_name)
    pdf(file = filename, width = fig_width, height = fig_height)
    heatmap.2(GEP_markerSum_zscore_matched[marker_idx,],
              col = bluered(100),
              dendrogram='none',
              Rowv=FALSE, 
              Colv=FALSE,
              cexRow = rowlabelsize,#0.7, 
              cexCol=collabelsize,#0.8,
              trace = "none",
              density.info = "none",
              labRow = cell_type_matched,
              margins=margins,#c(3,0), 
              keysize=keysize,
              key.par=list(mar=keypar),
              srtCol=srtcol,
              lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(1, 5), lwid=c(1, 10, 1))
    dev.off()
  }else{
    # plot heatmap this is one to one assignment 
    cell_assign <- solve_LSAP(exp(GEP_markerSum_zscore_matched),maximum = TRUE)
    cell_type_idx <- cell_assign[1:length(cell_assign)]
    
    GEP_matched <- cell_gep[,cell_type_idx]
    colnames(GEP_matched) <- cell_type_matched
    
    cell_type_assignment <- data.frame(cbind(cell.type.in.marker.list = rownames(GEP_markerSum_zscore_matched[,cell_type_idx]), 
                                             cell.type.in.cell.gep = colnames(GEP_markerSum_zscore_matched[,cell_type_idx])))
    
    filename <- paste0(fig_path,heatmap_name)
    pdf(file = filename, width = fig_width, height = fig_height)
    heatmap.2(GEP_markerSum_zscore_matched[,cell_type_idx],
              col = bluered(100),
              dendrogram='none',
              Rowv=FALSE, 
              Colv=FALSE,
              cexRow = rowlabelsize,#0.7, 
              cexCol=collabelsize,#0.8,
              trace = "none",
              density.info = "none",
              labRow = cell_type_matched,
              margins=margins,#c(3,0), 
              keysize=keysize,
              key.par=list(mar=keypar),
              srtCol=srtcol,
              lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(1, 5), lwid=c(1, 10, 1))
    dev.off()
    
  }
  
  # in the case where number of CDSeq estimated cell types is smaller than the cell types for matching
  # we do fuzzy assignment, meaning multiple CDSeq estimated cell types could match to one known cell types using marker genes
  # this is a slightly variant version of Munkes assignment problem, a possible to way of doing this the following
  # idea from : https://stackoverflow.com/questions/48108496/hungarian-algorithm-multiple-jobs-per-worker is interesting but cannot be applied in my case
  # I propose the following interative procedure for fuzzy assignment
  # when number of known cell type, i.e. ncol(GEP_markerSum_zscore_matched), is smaller than the number of CDSeq estimated cell type
  # then we could run Munkres assignment recursivly until all CDSeq-estimated cell type has been assigned.
  
  if(ncol(GEP_markerSum_zscore_matched) > nrow(GEP_markerSum_zscore_matched)){
    cell_type_matched_fuzzy <- matrix(0,ncol = ncol(GEP_markerSum_zscore_matched), nrow = nrow(GEP_markerSum_zscore_matched))
    colnames(cell_type_matched_fuzzy) <- colnames(GEP_markerSum_zscore_matched)
    rownames(cell_type_matched_fuzzy) <- rownames(GEP_markerSum_zscore_matched)
    
    # current assignment from previous Munkres solution
    if(length(cell_type_idx)!=nrow(GEP_markerSum_zscore_matched)){
      warning("There is an error in fuzzy assign. no figure plot for fuzzy assign.")
    }else{
      for (i in 1:length(cell_type_idx)) {
        cell_type_matched_fuzzy[i,cell_type_idx[i]] <- 1
      }
      # remove assigned CDSeq estimated cell types
      tmp_mat <- as.matrix(GEP_markerSum_zscore_matched[,-cell_type_idx])
      # get the index for unassigned cell types
      unassigned_idx <- c(1:ncol(GEP_markerSum_zscore_matched))[-cell_type_idx]
      
      # run Munkres recursively
      while (any( colSums(cell_type_matched_fuzzy)==0 )) {
        if(verbose){
          cat("ncol(tmp_mat) = ",ncol(tmp_mat)," nrow(tmp_mat) = ", nrow(tmp_mat),"length(unassigned_idx) = ",length(unassigned_idx),"\n")
        }
        if(ncol(tmp_mat)>=nrow(tmp_mat)){
          if(verbose){cat("in if...\n")}
          
          # case 1: number of un-assigned CDSeq estGEP is greater than cell types in marker gene list
          cell_assign <- solve_LSAP(exp(tmp_mat),maximum = TRUE)
          tmp_idx <- cell_assign[1:length(cell_assign)]
          tmp_mat <- as.matrix(tmp_mat[,-tmp_idx])# use as.matrix is because if tmp_mat left with one column, the ncol(tmp_mat) is NULL
          
          newly_assigned <- unassigned_idx[tmp_idx]
          unassigned_idx <- unassigned_idx[-tmp_idx]
          for (i in 1:length(tmp_idx)) {
            cell_type_matched_fuzzy[i,newly_assigned[i]] <- 1
          }
        }else{
          if(verbose) {cat("in else...\n")}
          # case 2: number of un-assigned CDSeq estGEP is smaller than cell types in marker gene list
          if(length(unassigned_idx)>0){
            cell_assign <- solve_LSAP(exp(t(tmp_mat)),maximum = TRUE)
            marker_idx <- cell_assign[1:length(cell_assign)]
            #newly_assigned <- unassigned_idx[tmp_idx]
            #unassigned_idx <- unassigned_idx[-tmp_idx]
            for (i in 1:length(marker_idx)) {
              cell_type_matched_fuzzy[marker_idx[i],unassigned_idx[i]] <- 1
              if(verbose){cat("marker_idx[i]=",marker_idx[i],"unassigned_idx[i]=",unassigned_idx[i],"\n")}
            }
          }
          # if(length(unassigned_idx)>0){
          #   rem_mat <- GEP_markerSum[,unassigned_idx]
          #   rem_idx <- max.col(t(rem_mat))
          #   tmp_mat <- rem_mat[rem_idx,]
          #   
          #   cell_assign <- solve_LSAP(exp(tmp_mat),maximum = TRUE)
          #   tmp_idx <- cell_assign[1:length(cell_assign)]
          #   #tmp_mat <- tmp_mat[,-tmp_idx]
          #   
          #   newly_assigned <- unassigned_idx[tmp_idx]
          #   unassigned_idx <- unassigned_idx[-tmp_idx]
          #   for (i in 1:length(tmp_idx)) {
          #     cell_type_matched_fuzzy[i,newly_assigned[i]] <- 1
          #   }
          #   
          # }
        }
      }
      
      if(sum(cell_type_matched_fuzzy)!=ncol(cell_type_matched_fuzzy)){warning("something wrong with fuzzy assignment.no figure plotted for fuzzy assign--1")}
      filename <- paste0(fig_path,heatmap_name_fuzzy_assign)
      pdf(file = filename, width = fig_width, height = fig_height)
      cell_type_assing_fuzzy_idx <- numeric(0)
      for (i in 1:nrow(cell_type_matched_fuzzy)) {
        cell_type_assing_fuzzy_idx <- c(cell_type_assing_fuzzy_idx, which(cell_type_matched_fuzzy[i,] == 1))
      }
      if(length(cell_type_assing_fuzzy_idx)!=ncol(cell_type_matched_fuzzy)){warning("something wrong with fuzzy assignment.no figure plotted for fuzzy assign--2")}
      heatmap.2(GEP_markerSum_zscore_matched[,cell_type_assing_fuzzy_idx],
                col = bluered(100),
                dendrogram='none',
                Rowv=FALSE, 
                Colv=FALSE,
                cexRow = rowlabelsize,#0.7, 
                cexCol=collabelsize,#0.8,
                trace = "none",
                density.info = "none",
                labRow = cell_type_matched,
                margins=margins,#c(3,0), 
                keysize=keysize,
                key.par=list(mar=keypar),
                srtCol=srtcol,
                lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(1, 5), lwid=c(1, 10, 1))
      dev.off()
    }
    
    tmp <- rep('a',ncol(cell_type_matched_fuzzy))
    for (i in 1:length(tmp)) {
      tmp[i] <- rownames(cell_type_matched_fuzzy)[ which(cell_type_matched_fuzzy[,i]==1) ]
    }
    cell_type_assignment_fuzzy <- data.frame(cbind(cell.type.in.marker.list = tmp, 
                                             cell.type.in.cell.gep = colnames(cell_type_matched_fuzzy)))
  }
  
  
  # if(!is.null(heatmap_name_twoway)){
  #   filename2 <- paste0(fig_path,heatmap_name_twoway)
  #   pdf(file = filename2, width = 60,height = 20,paper = "a4")
  #   lmat = rbind(c(0,3),c(2,1),c(0,4))
  #   lwid = c(1.5,4)
  #   lhei = c(1.5,4,1)
  #   distance.row <- dist(as.matrix(GEP_markerSum_zscore), method = "euclidean")
  #   cluster.row <- hclust(distance.row, method = "complete")
  #   distance.col <- dist(t(as.matrix(GEP_markerSum_zscore)),method = "euclidean")
  #   cluster.col <- hclust(distance.col, method = "complete")
  #   heatmap.2(cell_gep,scale="row",col = bluered(100),dendrogram='both',Rowv=TRUE, Colv = TRUE,#Colv=as.dendrogram(cluster.col),
  #             cexRow = 0.7, cexCol=0.8,trace = "none",density.info = "none",lmat = lmat, lwid = lwid, lhei = lhei)
  #   dev.off()
  # }
  
  output <- list()
  output$GEP_markerSum <- GEP_markerSum
  output$GEP_markerSum_zscore <- GEP_markerSum_zscore
  output$GEP_markerSum_zscore_matched <- GEP_markerSum_zscore_matched
  output$GEP_matched <- GEP_matched
  output$cell_type_matched_idx <- cell_type_idx
  output$cell_type_matched_fuzzy_idx <- cell_type_assing_fuzzy_idx
  output$cell_type_matched <- cell_type_matched
  output$cell_type_assign_fuzzy <- cell_type_matched_fuzzy
  output$cell_type_assignment_fuzzy <- cell_type_assignment_fuzzy
  output$cell_type_assignment <- cell_type_assignment
  #if(fig_marker){output$marker.plots <- marker.plots}
  return(output)
}