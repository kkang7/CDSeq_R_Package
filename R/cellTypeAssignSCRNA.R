#' \code{cellTypeAssignSCRNA} assigns CDSeq-identified cell types using single cell RNAseq data.
#' @param cdseq_gep CDSeq-estimated gene expression profile matrix with G rows (genes) and T columns (cell types).
#' @param cdseq_prop CDSeq-estimated sample-specific cell-type proportion, a matrix with T rows (cell type) and M (sample size).
#' @param cdseq_gep_sample_specific CDSeq-estimated sample-specific cell type gene expression, in the form of read counts. It is a 3 dimension array, i.e. gene by sample by cell type. The element cdseq_gep_sample_specific[i,j,k] represents the reads mapped to gene i from cell type k in sample j.
#' @param sc_gep a G (genes) by N (cell) matrix or dataframe that contains the gene expression profile for N single cells.
#' @param sc_annotation a dataframe contains two columns "cell_id"  and "cell_type". cell_id needs to match with the cell_id in sc_gep but not required to have the same size. cell_type is the cell type annotation for the single cells.
#' @param sc_batch a vector contains batch information of single cell data, i.e. sc_gep, and length(sc_batch) = ncol(sc_gep). 
#' @param batch_correction perform Harmony batch correction if it is 1.
#' @param harmony_iter Maximum number of rounds to run Harmony. One round of Harmony involves one clustering and one correction step.
#' @param harmony_cluster Maximum number of rounds to run clustering at each round of Harmony.
#' @param nb_size size parameter for negative binomial distribution, check rnbinom for details.
#' @param nb_mu mu parameter for negative binomial distribution, check rnbinom for details.
#' @param breaksList parameter for pheatmap controling the color scale. See pheatmap function for details.
#' @param corr_threshold if the correlation between CDSeq-estimated GEPs and the scRNAseq GEP is below this value, then it is considered the two cell types are not matching. 
#' @param pseudo_cell_count an integer indicating how many pseudo cells will be generated from CDSeq-estimated cell-type-specific gene expression profiles. Default values is 1. 
#' @param seurat_count_threshold this parameter will be passed to Seurat subset function (subset = nCount_RNA > seurat_count_threshold) for filtering out single cells whose total counts is less this threshold. 
#' @param seurat_scale_factor this parameter will be passed to scale.factor in Seurat function NormalizeData.
#' @param seurat_norm_method this parameter will be passed to normalization.method in Seurat function NormalizeData.
#' @param seurat_select_method this parameter will be passed to selection.method in Seurat function FindVariableFeatures
#' @param seurat_nfeatures this parameter will be passed to nfeatures in Seurat function FindVariableFeatures.
#' @param seurat_npcs this parameter will be passed to npcs in Seurat function RunPCA.
#' @param seurat_dims this parameter will be passed to dims in Seurat function FindNeighbors.
#' @param seurat_reduction this parameter will be passed to reduction in Seurat function FindNeighbors.
#' @param seurat_resolution this parameter will be passed to resolution in Seurat function FindClusters.
#' @param seurat_find_marker this parameter controls if run seurat FindMarker function, default is FALSE.
#' @param seurat_DE_test this parameter will be passed to test.use in Seurat function FindAllMarkers.
#' @param seurat_DE_logfc this parameter will be passed to logfc.threshold in Seurat function FindAllMarkers.
#' @param seurat_top_n_markers the number of top DE markers saved from Seurat output. 
#' @param sc_pt_size point size of single cell data in umap and tsne plots
#' @param cdseq_pt_size point size of CDSeq-estimated cell types in umap and tsne plots
#' @param plot_umap set 1 to plot umap figure of scRNAseq and CDSeq-estimated cell types, 0 otherwise.
#' @param plot_tsne set 1 to plot tsne figure of scRNAseq and CDSeq-estimated cell types, 0 otherwise.
#' @param plot_per_sample currently disabled for debugging
#' @param fig_save 1 or 0. 1 means save figures to local and 0 means do not save figures to local.
#' @param fig_path the location where the heatmap figure is saved. 
#' @param fig_name the name of umap and tsne figures. Umap figure will have the name of fig_name_umap_date and tsne figure will be named fig_name_tsne_date.
#' @param fig_format "pdf", "jpeg", or "png".
#' @param corr_heatmap_fontsize font size of the correlation heatmap between scRNAseq GEP and CDSeq-estimated GEPs.
#' @param fig_dpi figure dpi
#' @param verbose if TRUE, some calculation information will be print.
#' @importFrom grDevices rainbow
#' @importFrom stats rmultinom rnbinom var
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures ScaleData FindNeighbors FindClusters FindAllMarkers RunUMAP RunTSNE RunPCA 
#' @importFrom ggplot2 guide_legend guides aes scale_size_manual scale_shape_manual scale_fill_manual xlab ylab theme ggsave ggplot ggtitle geom_point element_text scale_colour_manual 
#' @importFrom dplyr top_n group_by %>%
#' @importFrom pheatmap pheatmap
#' @importFrom rlang .data
#' @importFrom Matrix colSums
#' @importFrom harmony RunHarmony
#' @export
#' @return cellTypeAssignSCRNA returns a list containing following fields: 
#' fig_path: same as the input fig_path
#' 
#' fig_name: same as the input fig_name
#' 
#' cdseq_synth_scRNA: synthetic scRNAseq data generated using CDSeq-estiamted GEPs
#' 
#' cdseq_scRNA_umap: ggplot figure of the umap outcome
#' 
#' cdseq_scRNA_tsne: ggplot figure of the tsne outcome
#' 
#' cdseq_synth_scRNA_seurat: Seurat object containing the scRNAseq combined with CDSeq-estimated cell types. Cell id for CDSeq-estimated cell types start with "CDSeq". 
#' 
#' seurat_cluster_purity: for all cells in a Seurat cluster i,  the ith value in seurat_cluster_purity is the proportion of the mostly repeated cell annotation from sc_annotation. 
#' For example, after Seurat clustering, suppose there are 100 cells in cluster 1, out of these 100 cells, 90 cells' annotation in sc_annotation is cell type A, then 
#' the fist value in seurat_cluster_purity is 0.9. This output can be used to assess the agreement between Seurat clustering and the given sc_annotation. 
#' 
#' seurat_unique_clusters: Unique Seurat cluster numbering. This can be used together with seurat_cluster_gold_label to match the Seurat clusters with given annotations.
#' 
#' seurat_cluster_gold_label: The cell type annotations for each unique Seurat cluster based on sc_annotation.
#' 
#' seurat_markers: DE genes for each Seurat cluster.
#' 
#' seurat_top_markers: Top seurat_top_n_markers DE genes for each Seurat cluster.
#' 
#' CDSeq_cell_type_assignment_df: cell type assignment for CDSeq-estimated cell types.
#' 
#' CDSeq_cell_type_assignment_confidence: cell type assignment confidence matrix, only available when pseudo_cell_count > 1.
#' 
#' CDSeq_cell_type_assignment_df_all: cell type assignment for CDSeq-estimated cell types, only available when pseudo_cell_count > 1.
#' 
#' cdseq_prop_merged: CDSeq-estimated cell type proportions with cell type annotations (annotated using clustering with scRNAseq).
#' 
#' cdseq_gep_sample_specific_merged: sample-specific cell-type read counts. It is a 3d array with dimensions: gene, sample, cell type. 
#' 
#' input_list: values for input parameters
#' 
#' cdseq_sc_comb_umap_df: dataframe for umap plot
#' 
#' cdseq_sc_comb_tsne_df: dataframe for tsne plot
#' 
#' cdseq_prop_merged_byCorr: CDSeq-estimated cell type proportions with cell type annotations (annotated using correlation with scRNAseq).
#' 
#' cdseq_gep_merged_byCorr: CDSeq-estimated cell-type-specific GEPs with cell type annotations (annotated using correlation with scRNAseq).
#' 
#' cdseq_annotation_byCorr: CDSeq-estimated cell type annotations (annotated using correlation with scRNAseq)


cellTypeAssignSCRNA <- function(cdseq_gep = NULL,
                                cdseq_prop = NULL,
                                cdseq_gep_sample_specific = NULL,
                                sc_gep = NULL,
                                sc_annotation = NULL,
                                sc_batch = NULL,
                                batch_correction = 1,
                                harmony_iter = 10,
                                harmony_cluster = 20,
                                nb_size = NULL,
                                nb_mu = NULL,
                                corr_threshold = 0,
                                breaksList = seq(0,1,0.01),
                                pseudo_cell_count = 1,
                                seurat_count_threshold = 0,
                                seurat_scale_factor = 10000,
                                seurat_norm_method = "LogNormalize",
                                seurat_select_method = 'vst',
                                seurat_nfeatures = 1000,
                                seurat_npcs = 30,
                                seurat_dims = 1:30,
                                seurat_reduction = 'pca',
                                seurat_resolution = 0.8,
                                seurat_find_marker = FALSE,
                                seurat_DE_test = "wilcox",
                                seurat_DE_logfc = 0.25,
                                seurat_top_n_markers = 10,
                                sc_pt_size = 1,
                                cdseq_pt_size = 3, 
                                plot_umap = 1,
                                plot_tsne = 1,
                                plot_per_sample = 0,
                                fig_save = 0,
                                fig_path = getwd(),
                                fig_name = "CDSeqCellTypeAssignSCRNA",
                                fig_format = "jpeg",
                                fig_dpi = 100,
                                corr_heatmap_fontsize = 10,
                                verbose=TRUE){
  #################################################################
  ##                         check input                         ##
  #################################################################
  # save all the inputs
  input_list <- list(cdseq_gep = cdseq_gep,
                     cdseq_prop = cdseq_prop,
                     cdseq_gep_sample_specific = cdseq_gep_sample_specific,
                     sc_gep = sc_gep,
                     sc_annotation = sc_annotation,
                     sc_batch = sc_batch,
                     batch_correction = batch_correction,
                     pseudo_cell_count = pseudo_cell_count,
                     nb_size = nb_size,
                     nb_mu = nb_mu,
                     corr_threshold = corr_threshold,
                     seurat_count_threshold = seurat_count_threshold,
                     seurat_scale_factor = seurat_scale_factor,
                     seurat_norm_method = seurat_norm_method,
                     seurat_select_method = seurat_select_method,
                     seurat_nfeatures = seurat_nfeatures,
                     seurat_npcs = seurat_npcs,
                     seurat_dims = seurat_dims,
                     seurat_reduction = seurat_reduction,
                     seurat_resolution = seurat_resolution,
                     seurat_DE_test = seurat_DE_test,
                     seurat_DE_logfc = seurat_DE_logfc,
                     seurat_top_n_markers = seurat_top_n_markers,
                     plot_umap = plot_umap,
                     plot_tsne = plot_tsne,
                     plot_per_sample = plot_per_sample,
                     fig_path = fig_path,
                     fig_name = fig_name,
                     fig_format = fig_format,
                     corr_heatmap_fontsize = 10,
                     fig_dpi = fig_dpi)
  
  
  #cat("CDSeq version 1.0.7\n")
  if(verbose){cat("checking input ... \n")}
  # check dimension
  nr_cdseq <- nrow(cdseq_gep)
  nc_cdseq <- ncol(cdseq_gep)
  if(!is.null(sc_gep)){
    nr_sc <- nrow(sc_gep)
    nc_sc <- ncol(sc_gep)
  }
  
  # check batch 
  if(!is.null(sc_batch)){
    if(length(sc_batch) != ncol(sc_gep)){stop("length(sc_batch) has to be equal to ncol(sc_gep).")}
  }else{
    if(!is.null(sc_gep)){
      warning("sc_bath is NOT given. Please provide batch information for scRNAseq if there is any.")
    }
  }
  
  # check input pseudo_cell_count
  if(is.null(pseudo_cell_count)){stop("pseudo_cell_count has to be an integer.")}
  if(!is.numeric(pseudo_cell_count)){stop("pseudo_cell_count has to an integer.")}
  if(length(pseudo_cell_count)>1){stop("pseudo_cell_count has to be a scalar.")}
  if(pseudo_cell_count<0){stop("pseudo_cell_count has to be greater than or equal to 1.")}
  if(is.null(sc_gep) && pseudo_cell_count < 100){
    stop("Increase pseudo_cell_count to be greater than 100 at least.")
  }
  # check gene names
  if(!is.null(sc_gep)){
    gene_cdseq <- rownames(cdseq_gep)
    gene_sc <- rownames(sc_gep)
    if(is.null(gene_cdseq) || is.null(gene_sc)){stop("cdseq_gep and sc_gep should have gene names as row names.")}
    #if(sum(gene_cdseq == gene_sc) != nr_cdseq ){stop("cdseq_gep and sc_gep should have the same gene names and in the same order.")}
    if(!identical(gene_cdseq,gene_sc)){stop("cdseq_gep and sc_gep should have the same gene names and in the same order.")}
  }

  # check cell type proportions
  if(!is.null(cdseq_prop)){
    nr_prop <- nrow(cdseq_prop)
    nc_prop <- ncol(cdseq_prop)
    if(!is.null(sc_gep)){
      if(nr_cdseq!=nr_sc){stop("cdseq_gep and sc_gep should have the same number of rows (genes).")}
    }
    if(nr_prop!=nc_cdseq){stop("column number of cdseq_gep should equal to row number of cdseq_prop.")}
    
    if(is.null(colnames(cdseq_prop))){stop("cdseq_prop column name is missing. Column names are sample ids.")}
    if(is.null(rownames(cdseq_prop))){stop("cdseq_prop row name is missing. Row names are the cell type names.")}
    if(sum(rownames(cdseq_prop) == colnames(cdseq_gep)) != nr_prop){stop("Row names of cdseq_prop and column names of cdseq_gep should be the same and in same order. They both denote cell types.")}
  }else{
    warning("CDSeq-estimated cell type proportions is not provided.")
  }
  
  # check cell ids
  if(!is.null(sc_gep)){
    sc_id <- colnames(sc_gep)
    if(is.null(sc_id)){stop("Please provide cell id in sc_gep.")}
  }
  
  # check annotation
  if(is.null(sc_annotation) && is.null(sc_gep)){
    warning("sc_annotation (single cell annotation) and single cell data are not provided. We will perform cluster and find cluster markers on CDSeq-estimated cell types.")
    seurat_find_marker <- TRUE
  }else if( !is.null(sc_annotation) && !is.null(sc_gep)){
    if(is.null(dim(sc_annotation))){stop("sc_annotation has to be a two dimensional data frame with column names: cell_id and cell_type")}
    if(nrow(sc_annotation)!=nc_sc){stop("nrow(sc_annotation) should be equal to ncol(sc_gep).")}
    #if(sum(colnames(sc_annotation) == c("cell_id","cell_type"))!=2){stop("sc_annotation has to be a two dimensional data frame with column names: cell_id and cell_type")}
    if(!identical(colnames(sc_annotation), c("cell_id","cell_type"))){stop("sc_annotation has to be a two dimensional data frame with column names: cell_id and cell_type")}
    #cat("sum(sc_id == sc_annotation$cell_id) = ",sum(sc_id == sc_annotation$cell_id),"length(sc_id)=",length(sc_id),"\n")
    #if(sum(sc_id == sc_annotation$cell_id) != length(sc_id)){stop("cell id in sc_gep has to be the same as the cell id in sc_annotation.")}
    if(!identical(sc_id,sc_annotation$cell_id)){stop("cell id in sc_gep has to be the same as the cell id in sc_annotation.")}
  }
  
  # construct figure names
  current_time <- format(Sys.time(),"%Y_%m_%d_%H_%M_%S")
  
  if(identical(fig_name,"CDSeqCellTypeAssignSCRNA")){
    # if user does not provide names, then use default names appended with current time
    fig_name_corr <- paste0(fig_name,"_corr_", current_time, ".",fig_format)
    fig_name_umap <- paste0(fig_name,"_UMAP_", current_time, ".",fig_format)
    fig_name_umap_perSample <- paste0(fig_name,"_UMAP_perSample_", current_time,".",fig_format)
    fig_name_tsne <- paste0(fig_name,"_TSNE_", current_time, ".",fig_format)
    fig_name_tsne_perSample <- paste0(fig_name,"_TSNE_perSample_", current_time,".", fig_format)
  }else{
    # if user provided names, then do not append current time
    fig_name_corr <- paste0(fig_name,"_corr.",fig_format)
    fig_name_umap <- paste0(fig_name,"_UMAP.",fig_format)
    fig_name_umap_perSample <- paste0(fig_name,"_UMAP_perSample." ,fig_format)
    fig_name_tsne <- paste0(fig_name,"_TSNE",fig_format)  
    fig_name_tsne_perSample <- paste0(fig_name,"_TSNE_perSample.",fig_format)
    
  }
  
  
  #####################################################################################################
  ##  cell type assignment using correlation with pseudo-bulk by merging scRNAseq of the same types  ##
  #####################################################################################################
  if(!is.null(sc_gep) && !is.null(sc_annotation)){
    if(verbose){cat("Annotating CDSeq-estimated cell type using correlation with scRNAseq ...\n")}
    sc_tmp <- sc_gep
    colnames(sc_tmp) <- sc_annotation[["cell_type"]]
    sc_merge <- t(rowsum(t(as.matrix(sc_tmp)),colnames(sc_tmp)))
    
    cdseq_correlation <- cor(cdseq_gep,sc_merge)
    cdseq_annotation_byCorr <- rep("unknown",ncol(cdseq_gep))
    for (i in 1:ncol(cdseq_gep)) {
      if(max(cdseq_correlation[i,])<corr_threshold){
        cdseq_annotation_byCorr[i] <- "Unknown"
      }else{
        cdseq_annotation_byCorr[i] <- names(which.max(cdseq_correlation[i,]))
      }
    }
    #=======================  merge gep ===================
    cdseq_gep_tmp <- cdseq_gep
    colnames(cdseq_gep_tmp) <- cdseq_annotation_byCorr
    
    cdseq_gep_merged_byCorr <- t(rowsum(t(as.matrix(cdseq_gep_tmp)),colnames(cdseq_gep_tmp)))
    comm_ct <- intersection(list.vector = list(colnames(sc_merge),colnames(cdseq_gep_merged_byCorr)), order = "stable")
    cdseq_gep_merged_byCorr <- cdseq_gep_merged_byCorr[,comm_ct$index[[2]],  drop = FALSE]
    colnames(cdseq_gep_merged_byCorr) <- paste0("CDSeq_",colnames(cdseq_gep_merged_byCorr))
    
    #breaksList <- seq(0,1,0.01)
    pheatmap(cor(cdseq_gep_merged_byCorr, sc_merge),
             fontsize = corr_heatmap_fontsize,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             breaks = breaksList,
             filename = file.path(fig_path,fig_name_corr) )
    if(verbose){cat("Heatmap of correlation between CDSeq-estimated GEP and scRNAseq GEP is saved to ", file.path(fig_path,fig_name_corr),"\n")}
    
    #======================= merge proportion ====================
    cdseq_prop_tmp <- cdseq_prop
    rownames(cdseq_prop_tmp) <- cdseq_annotation_byCorr
    
    cdseq_prop_merged_byCorr <- rowsum(t(as.matrix(t(cdseq_prop_tmp))),colnames(t(cdseq_prop_tmp)))
    rownames(cdseq_prop_merged_byCorr) <- paste0("CDSeq_",rownames(cdseq_prop_merged_byCorr))

  }

  ###################################################################
  ##  Generate synthetic scRNAseq data from CDSeq estimated GEPs   ##
  ##                 and combine with scRNAseq data                ##
  ###################################################################
  # declare variable to avoid R CMD check notes for no visible binding for global variable
  avg_logFC <- cluster <- nCount_RNA <- NULL
  
  if(verbose){cat("generating synthetic scRNAseq data using CDSeq estimated GEPs ...\n")}
  # generate synthetic scRNAseq 
  #ncell <-  1
  ncell <- round(pseudo_cell_count)
  cdseq_synth_scRNA <- matrix(0,nrow = nr_cdseq, ncol = nc_cdseq*ncell)
  
  # estimate mean and variance of scRNAseq and use those values for nb_mu and nb_size
  if(!is.null(sc_gep) && is.null(nb_mu) && is.null(nb_size)){
    #nb_mu <- max(Matrix::colSums(sc_gep)); nb_size <- nb_mu^2
    nb_mu <- mean(Matrix::colSums(sc_gep))# mean
    nb_size <- nb_mu^2/( stats::var(Matrix::colSums(sc_gep)) - nb_mu) # dispersion
  }
  
  if(is.null(sc_gep) && is.null(nb_mu) && is.null(nb_size)){
    nb_mu <- 5*nrow(cdseq_gep)
    nb_size <- nb_mu
  }
  
  for (i in 1:nc_cdseq) {
    nreads <- rnbinom(n = ncell,size = nb_size, mu = nb_mu)
    for (j in 1:ncell) {
      cdseq_synth_scRNA[,(i-1)*ncell + j] <-  rmultinom(1,nreads[j],cdseq_gep[,i])
    }
  }
  #colnames(cdseq_synth_scRNA) <- paste("CDSeq_SynthCell",rep(1:nc_cdseq,each = ncell), rep(1:ncell), sep = ".")  
  colnames(cdseq_synth_scRNA) <- paste(rep(colnames(cdseq_gep),each = ncell), "CDSeq" ,rep(1:ncell), sep = ".")  
  rownames(cdseq_synth_scRNA) <- rownames(cdseq_gep)
  
  if(!is.null(sc_gep)){
    cdseq_sc_comb <- cbind(sc_gep,cdseq_synth_scRNA)
  }else{
    if(verbose){cat("sc_gep is not provided. Please check cdseq_synth_scRNA_seurat_markers to further identify CDSeq-estimated cell types. \n")}
    cdseq_sc_comb <- cdseq_synth_scRNA
    seurat_find_marker <- TRUE
  }
  
  ##################################################################
  ##                          Run seurat                          ##
  ##################################################################
  if(verbose){cat("Calling Seurat pipeline ... \n")}
  cdseq_synth_scRNA_seurat <- Seurat::CreateSeuratObject(counts = cdseq_sc_comb, project = "CDSeqSyntheticData")
  
  # add batch information
  if(is.null(sc_batch) && !is.null(sc_gep)){
    cdseq_synth_scRNA_seurat$batch <-c(rep("scRNA",ncol(sc_gep)), rep("CDSeq",ncol(cdseq_synth_scRNA)))# train_batch#
  }else if(!is.null(sc_batch) && !is.null(sc_gep)){
    cdseq_synth_scRNA_seurat$batch <-c(sc_batch, rep("CDSeq",ncol(cdseq_synth_scRNA)))# train_batch#
  }else if(is.null(sc_batch) && is.null(sc_gep)){
    cdseq_synth_scRNA_seurat$batch <- rep("CDSeq",ncol(cdseq_synth_scRNA))# train_batch#
  }
  # filter cells
  cdseq_synth_scRNA_seurat <- subset(cdseq_synth_scRNA_seurat, subset = nCount_RNA > seurat_count_threshold )#nFeature_RNA > 20 & nFeature_RNA < 2500)
  # normalize
  cdseq_synth_scRNA_seurat <- Seurat::NormalizeData(cdseq_synth_scRNA_seurat, normalization.method = seurat_norm_method , scale.factor = seurat_scale_factor, verbose = FALSE)
  # select genes
  cdseq_synth_scRNA_seurat <- Seurat::FindVariableFeatures(cdseq_synth_scRNA_seurat, selection.method = seurat_select_method, nfeatures = seurat_nfeatures, verbose = FALSE)

  cdseq_synth_scRNA_seurat <- Seurat::ScaleData(cdseq_synth_scRNA_seurat,verbose = FALSE)
  
  cdseq_synth_scRNA_seurat <- Seurat::RunPCA(cdseq_synth_scRNA_seurat, npcs = seurat_npcs,verbose = FALSE) #features = VariableFeatures(object = cdseq_synth_scRNA_seurat), 
  
  # batch correction
  if(batch_correction){
    if(verbose) {cat("running Harmony..\n")}
    cdseq_synth_scRNA_seurat <- RunHarmony(cdseq_synth_scRNA_seurat,"batch", max.iter.harmony = harmony_iter, max.iter.cluster=harmony_cluster,verbose = FALSE)
    seurat_reduction <- "harmony"
  }
  
  cdseq_synth_scRNA_seurat <- Seurat::FindNeighbors(cdseq_synth_scRNA_seurat, dims = seurat_dims, reduction = seurat_reduction, verbose = FALSE)
  
  cdseq_synth_scRNA_seurat <- Seurat::FindClusters(cdseq_synth_scRNA_seurat, resolution = seurat_resolution, verbose = FALSE)
  
  ##################################################################
  ##               use seurat cluster and DE genes                ##
  ##            for CDSeq-estimatd cell type annotation           ##
  ##################################################################
  # find makers
  if(seurat_find_marker){
    if(verbose){cat("Finding markers for clusters...\n")}
    cdseq_synth_scRNA_seurat_markers <- Seurat::FindAllMarkers(cdseq_synth_scRNA_seurat, logfc.threshold = seurat_DE_logfc, test.use = seurat_DE_test,verbose = FALSE)
    #cdseq_synth_scRNA_seurat_markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 2, wt = avg_logFC)
    seurat_top_markers <- cdseq_synth_scRNA_seurat_markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = seurat_top_n_markers, wt = avg_logFC)
  }else{
    cdseq_synth_scRNA_seurat_markers <- NULL
    seurat_top_markers <- NULL
  }
 
  #seurat_top_markers_df <- as.data.frame(cbind(seurat_top_markers, singleCellAnnotation = rep("NA",nrow(seurat_top_markers))))

  if(!is.null(sc_annotation) && !is.null(sc_gep)){
    ##################################################################
    ##                  use scRNAseq annotation to                  ##
    ##              annotate CDSeq estimated cell types             ##
    ##################################################################
    if(seurat_find_marker){
      cdseq_synth_scRNA_seurat_markers_df <- cbind(cdseq_synth_scRNA_seurat_markers, 
                                                   singleCellAnnotation = rep("Unknown",nrow(cdseq_synth_scRNA_seurat_markers)),
                                                   stringsAsFactors = FALSE)
    }else{
      cdseq_synth_scRNA_seurat_markers_df <- NULL
    }
    
    # grep cell id and gene names in seurat output
    scRNAseq_cell_gene_name <- cdseq_synth_scRNA_seurat@assays$RNA@counts@Dimnames[[1]]
    scRNAseq_cell_id_seurat <- cdseq_synth_scRNA_seurat@assays$RNA@counts@Dimnames[[2]]

    # construct gold standard cluster label using sc_annotation
    scRNA_seurat_comm_cell <- intersection(list(scRNAseq_cell_id_seurat,sc_annotation$cell_id),order = 'stable')
    scRNA_seurat_comm_cell_id <- scRNA_seurat_comm_cell$comm.value
    scRNA_seurat_comm_idx <- scRNA_seurat_comm_cell$index[[2]]
    seurat_scRNA_comm_idx <- scRNA_seurat_comm_cell$index[[1]]
    
    scRNAseq_seurat_goldstandard_annotation <- sc_annotation$cell_type[scRNA_seurat_comm_idx]
    scRNAseq_seurat_goldstandard_annotation_cell_id <- scRNAseq_cell_id_seurat[seurat_scRNA_comm_idx]
    if(sum(scRNAseq_seurat_goldstandard_annotation_cell_id == sc_annotation$cell_id[scRNA_seurat_comm_idx]) != length(scRNAseq_seurat_goldstandard_annotation_cell_id)){
      stop("Something is wrong when taking intersection between sc_annotation cell id and seurat cell id.")
    }
    
    # extract seurat cluster label for the single cells with cell type annotations
    # NOTICE THAT [scRNAseq_cell_cluster_label_seurat] and [scRNAseq_seurat_goldstandard_annotation] should be in the same order in term of cell ids
    scRNAseq_cell_cluster_label_seurat <- cdseq_synth_scRNA_seurat@meta.data$seurat_clusters[seurat_scRNA_comm_idx]
    
    cdseq_cell_idx <- grep("CDSeq",scRNAseq_cell_id_seurat)
    cdseq_cell_seurat_cluster_id <- cdseq_synth_scRNA_seurat@meta.data$seurat_clusters[cdseq_cell_idx]
    cdseq_cell_id <- scRNAseq_cell_id_seurat[cdseq_cell_idx]
    
    
    CDSeq_cell_type_assignment_df <- data.frame(cdseq_cell_id = cdseq_cell_id, 
                                                seurat_cluster = cdseq_cell_seurat_cluster_id, 
                                                cdseq_cell_type_assignment = rep("unknown",length(cdseq_cell_id)),
                                                stringsAsFactors = FALSE)

    
    # compute the cluster score (there is a bug due to the large number of cells )
    #seurat_cluster_scores <- cluster_scores(c = as.numeric(scRNAseq_cell_cluster_label_seurat) , k = as.numeric(scRNAseq_seurat_goldstandard_annotation) ,beta = 1)
    
    # assign CDSeq-estimated cell types based on single cell annotation
    seurat_unique_clusters <- sort(as.numeric(unique(scRNAseq_cell_cluster_label_seurat)))
    seurat_cluster_purity <- rep(0, length(seurat_unique_clusters))
    seurat_cluster_gold_label <- rep("unknown",length(seurat_unique_clusters))
    for (i in 1:length(seurat_unique_clusters)) {
      seurat_cluster_members <- which(as.numeric(scRNAseq_cell_cluster_label_seurat) == seurat_unique_clusters[i])
      seurat_cluster_gold_standard_label <- as.character(scRNAseq_seurat_goldstandard_annotation[seurat_cluster_members])
      seurat_cluster_assess <- max_rep(seurat_cluster_gold_standard_label)
      if(length(seurat_cluster_assess$max_element)>1){
        warning("In cluster ", i, ", there s a tie between two or more possible cell type annotations. A random one is chosen. Try to reduce pseudo_cell_count or adjust resolution parameter.")
        seurat_cluster_purity[i] <- seurat_cluster_assess$max_element_proportion
        seurat_cluster_gold_label[i] <- seurat_cluster_assess$max_element[1]
      }else{
        seurat_cluster_purity[i] <- seurat_cluster_assess$max_element_proportion
        seurat_cluster_gold_label[i] <- seurat_cluster_assess$max_element
      }
      # add sc_annotation to DE marker data frame
      if(seurat_find_marker){
        tmp_idx <- which(as.numeric(cdseq_synth_scRNA_seurat_markers_df$cluster) == seurat_unique_clusters[i])
        cdseq_synth_scRNA_seurat_markers_df$singleCellAnnotation[tmp_idx] <- seurat_cluster_gold_label[i]
      }
      cdseq_tmp_idx <- which(as.numeric(CDSeq_cell_type_assignment_df$seurat_cluster) == seurat_unique_clusters[i])
      CDSeq_cell_type_assignment_df$cdseq_cell_type_assignment[cdseq_tmp_idx] <- seurat_cluster_gold_label[i]
    }
    
    # if pseudo_cell_count > 1, construct a cell type assignment confidence matrix
    # each row of this matrix denotes a CDSeq-estimated cell type and each column of this matrix denotes a cell type in sc_annotation
    # entry (i, j) is a value in [0, 1] denoting the proportion of the pseudo_cell_count cells for the ith CDSeq-identified cell types 
    # annotated as the jth cell types in sc_annotation 
    if(pseudo_cell_count>1){
      # create confidence matrix, each column denotes a CDSeq-identified cell type, each row denotes a cell type in sc_annotation
      CDSeq_cell_type_assignment_confidence <- matrix(0, ncol = ncol(cdseq_gep), nrow = length(unique(sc_annotation$cell_type)))
      colnames(CDSeq_cell_type_assignment_confidence) <- colnames(cdseq_gep)
      rownames(CDSeq_cell_type_assignment_confidence) <- unique(sc_annotation$cell_type)
      
      # reduce CDSeq_cell_type_assignment_df to one CDSeq-estimated cell type matching one cell type assignment
      CDSeq_cell_type_assignment_df_all <- CDSeq_cell_type_assignment_df
      
      CDSeq_cell_type_assignment_df <- data.frame(cdseq_cell_id = colnames(cdseq_gep), 
                                                  #seurat_cluster = rep("not_assigned",ncol(cdseq_gep)), 
                                                  cdseq_cell_type_assignment = rep("unknown",ncol(cdseq_gep)),
                                                  stringsAsFactors = FALSE)
      # merge annotations
      if(verbose) {cat("merging annotations...\n")}
      for (i in 1:ncol(CDSeq_cell_type_assignment_confidence)) {
        # get the idx for each CDSeq estimated cell types
        temp_idx <- grep(colnames(cdseq_gep)[i], CDSeq_cell_type_assignment_df_all$cdseq_cell_id)
        # get the corresponding cell types
        temp_ct_df <- as.data.frame(base::table(CDSeq_cell_type_assignment_df_all$cdseq_cell_type_assignment[temp_idx]))
        temp_comm <- intersection(list.vector = list(rownames(CDSeq_cell_type_assignment_confidence), as.character(temp_ct_df[[1]])), order="stable")
        temp_mat_idx_1 <- temp_comm$index[[1]]#which( rownames(CDSeq_cell_type_assignment_confidence) %in% as.character(temp_ct_df[[1]]))
        temp_mat_idx_2 <- temp_comm$index[[2]]#match(as.character(temp_ct_df[[1]]), rownames(CDSeq_cell_type_assignment_confidence) )
        CDSeq_cell_type_assignment_confidence[temp_mat_idx_1 , i] <- temp_ct_df[[2]][temp_mat_idx_2] #
        # refill CDSeq_cell_type_assignment_df
        #temp_cluster_max <- max_rep(CDSeq_cell_type_assignment_df_all$seurat_cluster[temp_idx])
        # cat("i = ", i, "---",colnames(CDSeq_cell_type_assignment_confidence)[i],
        #     " ---- length(temp_cluster_max$max_element) = ", length(temp_cluster_max$max_element),
        #     "\n")
        # for (aa in 1:length(temp_cluster_max$max_element)) {
        #   cat("aa = ",aa,"temp_cluster_max$max_element[[",aa,"]]=",temp_cluster_max$max_element[[aa]],"\n ---------  \n")
        # }
        # cat("======================================\n")
        #CDSeq_cell_type_assignment_df$seurat_cluster[i] <- temp_cluster_max$max_element
        temp_ct_max <- max_rep(CDSeq_cell_type_assignment_df_all$cdseq_cell_type_assignment[temp_idx])
        temp_idx_2 <- grep(colnames(cdseq_gep)[i],CDSeq_cell_type_assignment_df$cdseq_cell_id)
        if(length(temp_ct_max$max_element)>1){
          warning("There is a tie between two or more cell types. A random cell type is chosen. Try to reduce pseudo_cell_count or adjust resolution parameter.")
          CDSeq_cell_type_assignment_df$cdseq_cell_type_assignment[temp_idx_2] <- temp_ct_max$max_element[1]
        }else{
          CDSeq_cell_type_assignment_df$cdseq_cell_type_assignment[temp_idx_2] <- temp_ct_max$max_element
        }
      }
      # normalize 
      CDSeq_cell_type_assignment_confidence <- t(t(CDSeq_cell_type_assignment_confidence)/colSums(CDSeq_cell_type_assignment_confidence))
      
      
    }else{
      CDSeq_cell_type_assignment_df_all <- NULL
      CDSeq_cell_type_assignment_confidence <- NULL
    }
    
    if(seurat_find_marker){
      seurat_top_markers_df <- cdseq_synth_scRNA_seurat_markers_df %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = seurat_top_n_markers, wt = avg_logFC)
    }else{
      seurat_top_markers_df <- NULL
    }
    
    if(!is.null(cdseq_prop)){
       # current version assumes each CDSeq-estGEP only generate 1 synthetic single cell. May extend to more general case later
       # the gsub function extract patterns between two dots.
       #cdseq_est_cell_ids <- as.numeric(unlist(lapply(CDSeq_cell_type_assignment_df$cdseq_cell_id, function(x){gsub("^.+?\\.(.+?)\\..*$", "\\1",x)})))
      
       # the sub function extract pattern before a dot
       #cdseq_est_cell_ids <- unlist(lapply(CDSeq_cell_type_assignment_df$cdseq_cell_id, function(x){sub("\\..*","",x)}))
       #if(sum(cdseq_est_cell_ids == rownames(cdseq_prop)) != nr_prop){stop("Please make sure the rownames(cdseq_prop) and colnames(cdseq_gep) are the same and in same order.")}
       if(!identical(rownames(cdseq_prop), sub("\\..*","",CDSeq_cell_type_assignment_df[["cdseq_cell_id"]]))){
         cat("rownames(cdseq_prop)",rownames(cdseq_prop),"\n")
         cat(sub("\\..*","",CDSeq_cell_type_assignment_df$cdseq_cell_id),"\n")
         stop("Please make sure the rownames(cdseq_prop) and colnames(cdseq_gep) are the same and in same order. Or lower seurat_count_threshold to make sure some of the CDSeq estimated cell types are NOT filtered out.")
       } 
      rownames(cdseq_prop) <- CDSeq_cell_type_assignment_df$cdseq_cell_type_assignment
       cdseq_prop_merged <- rowsum(cdseq_prop, group = rownames(cdseq_prop))
    }
    
    if(!is.null(cdseq_gep_sample_specific)){
      cdseq_gep_sample_specific_merged <- array(numeric(),c(dim(cdseq_gep_sample_specific)[1], dim(cdseq_gep_sample_specific)[2],length(unique(CDSeq_cell_type_assignment_df$cdseq_cell_type_assignment)) ))
      if(is.null(dimnames(cdseq_gep_sample_specific))){
        if(is.null(cdseq_prop)){stop("Need to provide cdseq_prop or provide dimnames for cdseq_gep_sample_specific.")}
        message("dimnames(cdseq_gep_sample_specific) is NULL. cellTypeAssignSCRNA assumes its dimnames are same as: rownames(cdseq_gep), colnames(cdseq_prop), rownames(cdseq_prop)")
        dimnames(cdseq_gep_sample_specific)[[1]] <- rownames(cdseq_gep)
        dimnames(cdseq_gep_sample_specific)[[2]] <- colnames(cdseq_prop)
        dimnames(cdseq_gep_sample_specific)[[3]] <- rownames(cdseq_prop) 
        for (i in 1:dim(cdseq_gep_sample_specific)[2]) {
          cdseq_gep_sample_specific_merged[,i,] <- t(rowsum(t(cdseq_gep_sample_specific[,i,]), group = rownames(t(cdseq_gep_sample_specific[,i,]))))
        }
        dimnames(cdseq_gep_sample_specific_merged)[[1]] <- dimnames(cdseq_gep_sample_specific)[[1]]
        dimnames(cdseq_gep_sample_specific_merged)[[2]] <- dimnames(cdseq_gep_sample_specific)[[2]]
        # accorading to rowsum funtion, if reorder = TRUE (default), then the result will be in order of sort(unique(group))
        dimnames(cdseq_gep_sample_specific_merged)[[3]] <- sort(unique(dimnames(cdseq_gep_sample_specific)[[3]]))
      }else{
        cdseq_est_cell_ids <- unlist(lapply(CDSeq_cell_type_assignment_df$cdseq_cell_id, function(x){sub("\\..*","",x)}))
        if(sum(cdseq_est_cell_ids == dimnames(cdseq_gep_sample_specific)[[3]]) != nr_prop){stop("Please make sure the dimnames(cdseq_gep_sample_specific)[[3]] and colnames(cdseq_gep) are the same and in same order.")}
        dimnames(cdseq_gep_sample_specific)[[3]] <- CDSeq_cell_type_assignment_df$cdseq_cell_type_assignment
        for (i in 1:dim(cdseq_gep_sample_specific)[2]) {
          cdseq_gep_sample_specific_merged[,i,] <- t(rowsum(t(cdseq_gep_sample_specific[,i,]), group = rownames(t(cdseq_gep_sample_specific[,i,]))))
        }
        dimnames(cdseq_gep_sample_specific_merged)[[1]] <- dimnames(cdseq_gep_sample_specific)[[1]]
        dimnames(cdseq_gep_sample_specific_merged)[[2]] <- dimnames(cdseq_gep_sample_specific)[[2]]
        # accorading to rowsum funtion, if reorder = TRUE (default), then the result will be in order of sort(unique(group))
        dimnames(cdseq_gep_sample_specific_merged)[[3]] <- sort(unique(dimnames(cdseq_gep_sample_specific)[[3]]))
      }
      
      ####################################################################
      ##  use cell-type-specific-per-sample read counts for clustering  ##
      ####################################################################
      cdseq_gep_sample_specific_mat <- matrix(cdseq_gep_sample_specific,nrow = dim(cdseq_gep_sample_specific)[1],ncol = prod(dim(cdseq_gep_sample_specific)[2:3]))
      rownames(cdseq_gep_sample_specific_mat) <- rownames(cdseq_gep)
      colnames(cdseq_gep_sample_specific_mat) <- paste0("CDSeq_",
                                                        dimnames(cdseq_gep_sample_specific)[[2]],
                                                        "_",
                                                        rep(dimnames(cdseq_gep_sample_specific)[[3]],each=length(dimnames(cdseq_gep_sample_specific)[[2]])))
      cdseq_gep_sample_specific_mat_sample_ids <- rep(dimnames(cdseq_gep_sample_specific)[[2]], each = length(dimnames(cdseq_gep_sample_specific)[[3]]) )
      cdseq_persample_sc_comb <- cbind(sc_gep,cdseq_gep_sample_specific_mat)
      if(verbose){cat("colnames(cdseq_persample_sc_comb)[ncol(sc_gep) + 1] = ", colnames(cdseq_persample_sc_comb)[ncol(sc_gep) + 1],"\n")}
      ##################################################################
      ##                          Run seurat                          ##
      ##################################################################
      if(verbose){cat("Calling Seurat pipeline: per-sample-cell-type-specific ... \n")}
      cdseq_synth_scRNA_seurat_persample <- Seurat::CreateSeuratObject(counts = cdseq_persample_sc_comb, project = "cdseq_synth_scRNAseq_per_sample")
      # filter cells
      cdseq_synth_scRNA_seurat_persample <- subset(cdseq_synth_scRNA_seurat_persample, subset = nCount_RNA > seurat_count_threshold,verbose = FALSE )#nFeature_RNA > 20 & nFeature_RNA < 2500)
      # normalize
      cdseq_synth_scRNA_seurat_persample <- Seurat::NormalizeData(cdseq_synth_scRNA_seurat_persample, normalization.method = seurat_norm_method , scale.factor = seurat_scale_factor,verbose = FALSE)
      # select genes
      cdseq_synth_scRNA_seurat_persample <- Seurat::FindVariableFeatures(cdseq_synth_scRNA_seurat_persample, selection.method = seurat_select_method, nfeatures = seurat_nfeatures,verbose = FALSE)
      
      cdseq_synth_scRNA_seurat_persample <- Seurat::ScaleData(cdseq_synth_scRNA_seurat_persample,verbose = FALSE)
      
      cdseq_synth_scRNA_seurat_persample <- Seurat::RunPCA(cdseq_synth_scRNA_seurat_persample, npcs = seurat_npcs,verbose = FALSE) #features = VariableFeatures(object = cdseq_synth_scRNA_seurat), 
      
      cdseq_synth_scRNA_seurat_persample <- Seurat::FindNeighbors(cdseq_synth_scRNA_seurat_persample, dims = seurat_dims, reduction = seurat_reduction,verbose = FALSE)
      
      cdseq_synth_scRNA_seurat_persample <- Seurat::FindClusters(cdseq_synth_scRNA_seurat_persample, resolution = seurat_resolution,verbose = FALSE)
      
      ##################################################################
      ##   CDSeq-estimated per-sample-cell-type-specific annotation   ##
      ##################################################################
      # grep single cell id in seurat output
      scRNAseq_cell_gene_name <- cdseq_synth_scRNA_seurat_persample@assays$RNA@counts@Dimnames[[1]]
      scRNAseq_cell_id_seurat <- cdseq_synth_scRNA_seurat_persample@assays$RNA@counts@Dimnames[[2]]
      if(verbose){cat("length(scRNAseq_cell_id_seurat) = ", length(scRNAseq_cell_id_seurat),"\n" )}
      
      # construct gold standard cluster label using sc_annotation
      scRNA_seurat_comm_cell <- intersection(list(scRNAseq_cell_id_seurat,sc_annotation$cell_id),order = 'stable')
      scRNA_seurat_comm_cell_id <- scRNA_seurat_comm_cell$comm.value
      scRNA_seurat_comm_idx <- scRNA_seurat_comm_cell$index[[2]]
      seurat_scRNA_comm_idx <- scRNA_seurat_comm_cell$index[[1]]
      
      scRNAseq_seurat_goldstandard_annotation <- sc_annotation$cell_type[scRNA_seurat_comm_idx]
      scRNAseq_seurat_goldstandard_annotation_cell_id <- scRNAseq_cell_id_seurat[seurat_scRNA_comm_idx]
      if(sum(scRNAseq_seurat_goldstandard_annotation_cell_id == sc_annotation$cell_id[scRNA_seurat_comm_idx]) != length(scRNAseq_seurat_goldstandard_annotation_cell_id)){
        stop("Something is wrong when taking intersection between sc_annotation cell id and seurat cell id.")
      }
      
      # construct seurat cluster label 
      # NOTICE THAT [scRNAseq_cell_cluster_label_seurat] and [scRNAseq_seurat_goldstandard_annotation] should be in the same order in term of cell ids
      scRNAseq_cell_cluster_label_seurat_persample <- cdseq_synth_scRNA_seurat_persample@meta.data$seurat_clusters[seurat_scRNA_comm_idx]
      
      cdseq_cell_idx <- grep("CDSeq",scRNAseq_cell_id_seurat)
      cdseq_cell_seurat_cluster_id <- cdseq_synth_scRNA_seurat_persample@meta.data$seurat_clusters[cdseq_cell_idx]
      cdseq_cell_id <- scRNAseq_cell_id_seurat[cdseq_cell_idx]
      
      if(verbose){
        cat("length(cdseq_cell_id) = ", length(cdseq_cell_id),
            " length(cdseq_gep_sample_specific_mat_sample_ids) = ",length(cdseq_gep_sample_specific_mat_sample_ids),
            " length(cdseq_cell_seurat_cluster_id) = ",length(cdseq_cell_seurat_cluster_id),"\n")
      }
      
      CDSeq_cell_type_assignment_df_persample <- data.frame(cdseq_cell_id = cdseq_cell_id, 
                                                            #sample_id = cdseq_gep_sample_specific_mat_sample_ids, 
                                                            seurat_cluster = cdseq_cell_seurat_cluster_id, 
                                                            cdseq_cell_type_assignment = rep("not_assigned",length(cdseq_cell_id)),stringsAsFactors = FALSE)
      
      # compute the cluster score (there is a bug due to the large number of cells )
      #seurat_cluster_scores <- cluster_scores(c = as.numeric(scRNAseq_cell_cluster_label_seurat) , k = as.numeric(scRNAseq_seurat_goldstandard_annotation) ,beta = 1)
      
      # assign CDSeq-estimated cell types based on single cell annotation
      seurat_unique_clusters_persample <- sort(as.numeric(unique(scRNAseq_cell_cluster_label_seurat_persample)))
      seurat_cluster_purity <- rep(0, length(seurat_unique_clusters_persample))
      seurat_cluster_gold_label_persample <- rep("no_label",length(seurat_unique_clusters_persample))
      for (i in 1:length(seurat_unique_clusters_persample)) {
        seurat_cluster_members <- which(as.numeric(scRNAseq_cell_cluster_label_seurat_persample) == seurat_unique_clusters_persample[i])
        seurat_cluster_gold_standard_label <- as.character(scRNAseq_seurat_goldstandard_annotation[seurat_cluster_members])
        seurat_cluster_assess <- max_rep(seurat_cluster_gold_standard_label)
        seurat_cluster_purity[i] <- seurat_cluster_assess$max_element_proportion
        seurat_cluster_gold_label_persample[i] <- seurat_cluster_assess$max_element
        # add sc_annotation to DE marker data frame
        #tmp_idx <- which(as.numeric(cdseq_synth_scRNA_seurat_markers_df$cluster) == seurat_unique_clusters[i])
        #cdseq_synth_scRNA_seurat_markers_df$singleCellAnnotation[tmp_idx] <- seurat_cluster_gold_label[i]
        
        cdseq_tmp_idx <- which(as.numeric(CDSeq_cell_type_assignment_df_persample$seurat_cluster) == seurat_unique_clusters_persample[i])
        CDSeq_cell_type_assignment_df_persample$cdseq_cell_type_assignment[cdseq_tmp_idx] <- seurat_cluster_gold_label_persample[i]
      }
      
      # annotation cdseq-estimate proportions
      
      
    }
  }

  
  ##################################################################
  ##                             UMAP                             ##
  ##################################################################
  if(plot_umap){
    if(verbose){cat("Running UMAP: plot synthetic CDSeq-scRNAseq ... \n")}
    #################################################################
    ##                plot synthetic CDSeq-scRNAseq                ##
    #################################################################
    cdseq_synth_scRNA_seurat <- Seurat::RunUMAP(cdseq_synth_scRNA_seurat, dims = seurat_dims,reduction = seurat_reduction ,verbose = FALSE)
    
    # plot scRNAseq and cdseq estimates
    if(is.null(sc_annotation)){
      cluster_label <- paste0("cluster_",cdseq_synth_scRNA_seurat$seurat_clusters)#sub_grp#clusterlouvain$membership
    }else if(!is.null(sc_annotation) && !is.null(sc_gep)){
      tmp_label <- as.numeric(cdseq_synth_scRNA_seurat$seurat_clusters)
      cluster_label <- rep("unknown",length(tmp_label))
      for (i in 1:length(seurat_unique_clusters)) {
        cluster_label[which(tmp_label==seurat_unique_clusters[i])] <- seurat_cluster_gold_label[i]
      }
    }
    

    cluster_label <- factor(cluster_label, levels = unique(cluster_label))
    cdseq_sc_comb_umap <- cdseq_synth_scRNA_seurat@reductions$umap@cell.embeddings
    
    cell_sources <- rownames(cdseq_sc_comb_umap)
    CDSeq_idx <- grep("CDSeq",rownames(cdseq_sc_comb_umap))
    cell_sources[CDSeq_idx] <- "CDSeq"
    point_stroke <- c(rep(1,length(CDSeq_idx)))
    
    if(!is.null(sc_gep)){
      umap_tot <- 1:nrow(cdseq_sc_comb_umap)
      scRNA_idx <- umap_tot[-CDSeq_idx]
      cell_sources[scRNA_idx] <- "scRNA"
      cell_sources <- factor(cell_sources,levels = c('scRNA','CDSeq') )
      point_stroke <- c(rep(0,length(scRNA_idx)), rep(2,length(CDSeq_idx)))
    }
    
    #cell_sources <- factor(c(rep('scRNA',length(scRNA_idx)), rep('CDSeq',length(CDSeq_idx))), levels = c('scRNA','CDSeq'))
    cdseq_sc_comb_umap_df <- data.frame(V1 = cdseq_sc_comb_umap[,1], V2 = cdseq_sc_comb_umap[,2], cell_sources= cell_sources , cluster_label = cluster_label)
    
    #cat("nrow(cdseq_sc_comb_umap_df) = ", nrow(cdseq_sc_comb_umap_df), "length(point_stroke) = ", length(point_stroke),"...\n")

    if(!is.null(sc_gep)){
      cdseq_scRNA_umap<- ggplot(cdseq_sc_comb_umap_df,aes(x=.data$V1, y=.data$V2, colour=as.factor(cell_sources), size=as.factor(cell_sources), shape=as.factor(cell_sources), fill=as.factor(cluster_label), stroke=point_stroke)) + 
        ggtitle(paste0('CDSeq-estimated cell types and scRNAseq')) + 
        xlab("UMAP 1") + ylab("UMAP 2") +
        geom_point() + # this is the edge color
        scale_colour_manual(name="cell source",values=c("gray","black"),na.value = NA,label=levels(cell_sources)) +
        scale_fill_manual(name="clusters",values=c(rainbow(length(unique(cluster_label)))),na.value = NA,label=levels(cluster_label)) +
        scale_size_manual(name="cell source", values = c(sc_pt_size,cdseq_pt_size), na.value = NA,label=levels(cell_sources)) +
        scale_shape_manual(name="cell source", values = c(21,23), na.value = NA,label=levels(cell_sources)) +
        theme(plot.title = element_text(size = 45, face = "bold"), 
              axis.title.x = element_text(color="black", size=35,face="bold"),
              axis.title.y = element_text(color="black", size=35,face="bold"),
              legend.title = element_text(color = "blue", face = "bold", size = 30),
              legend.text = element_text(color = "blue", face = "bold", size = 30)) + 
        guides(color = guide_legend(override.aes = list(size=5)), fill = guide_legend(override.aes = list(shape=21,size=5)))
    }else{
      cdseq_scRNA_umap<- ggplot(cdseq_sc_comb_umap_df,aes(x=.data$V1, y=.data$V2, fill=as.factor(cluster_label), stroke=point_stroke)) + 
        ggtitle(paste0('CDSeq-estimated cell types')) + 
        xlab("UMAP 1") + ylab("UMAP 2") +
        geom_point(size = cdseq_pt_size, shape=21) + # this is the edge color
        #scale_colour_manual(name="cell source",values=c("gray","black"),na.value = NA,label=levels(cell_sources)) +
        scale_fill_manual(name="clusters",values=c(rainbow(length(unique(cluster_label)))),na.value = NA,label=levels(cluster_label)) +
        #scale_size_manual(name="cell source", values = c(sc_pt_size,cdseq_pt_size), na.value = NA,label=levels(cell_sources)) +
        #scale_shape_manual(name="cell source", values = c(21,23), na.value = NA,label=levels(cell_sources)) +
        theme(plot.title = element_text(size = 45, face = "bold"), 
              axis.title.x = element_text(color="black", size=35,face="bold"),
              axis.title.y = element_text(color="black", size=35,face="bold"),
              legend.title = element_text(color = "blue", face = "bold", size = 30),
              legend.text = element_text(color = "blue", face = "bold", size = 30)) + 
        guides(color = guide_legend(override.aes = list(size=5)), fill = guide_legend(override.aes = list(shape=21,size=5)))
      
    }
    
    if(fig_save){
      #fig_tmp_name <- paste0(fig_path,fig_name,"_umap_",gsub("-|\\s|:","_",Sys.time()),".",fig_format)
      #if(verbose){cat("save umap at ",fig_tmp_name ,"...\n")}
      if(verbose){cat("save umap at ",file.path(fig_path,fig_name_umap) ,"...\n")}
      ggsave(filename = fig_name_umap,
             plot = cdseq_scRNA_umap,
             path = fig_path,
             width = 25,
             height = 20,
             dpi = fig_dpi)
    }
    
    ##################################################################
    ##    plot CDSeq-estimated cell-type-specific gep per sample    ##
    ##################################################################
    if(plot_per_sample){
      warning("plot_per_sample option is not working properly. I'm working on it.")
      if(verbose){cat("Running UMAP: plot CDSeq-estimated cell-type-specific gep per sample ... \n")}
      
      cdseq_synth_scRNA_seurat_persample <- Seurat::RunUMAP(cdseq_synth_scRNA_seurat_persample, reduction = seurat_reduction ,dims = seurat_dims,verbose = FALSE)
      
      # plot scRNAseq and cdseq estimates
      if(is.null(sc_annotation)){
        cluster_label <- paste0("cluster_",cdseq_synth_scRNA_seurat_persample$seurat_clusters)#sub_grp#clusterlouvain$membership
      }else{
        tmp_label <- as.numeric(cdseq_synth_scRNA_seurat_persample$seurat_clusters)
        cluster_label <- rep("no_label",length(tmp_label))
        for (i in 1:length(seurat_unique_clusters_persample)) {
          cluster_label[which(tmp_label==seurat_unique_clusters_persample[i])] <- seurat_cluster_gold_label_persample[i]
        }
      }
      
      cluster_label <- factor(cluster_label, levels = unique(cluster_label))
      cdseq_sc_comb_umap <- cdseq_synth_scRNA_seurat_persample@reductions$umap@cell.embeddings
      
      cell_sources <- rownames(cdseq_sc_comb_umap)
      umap_tot <- 1:nrow(cdseq_sc_comb_umap)
      CDSeq_idx <- grep("CDSeq",rownames(cdseq_sc_comb_umap))
      scRNA_idx <- umap_tot[-CDSeq_idx]
      
      if(verbose){
        cat("length(CDSeq_idx) = ", length(CDSeq_idx), 
            " length(scRNA_idx) = ", length(scRNA_idx),
            " length(cell_sources) = ",length(cell_sources),
            " length(cluster_label) = ",length(cluster_label),
            " nrow(cdseq_sc_comb_umap) = ",nrow(cdseq_sc_comb_umap),
            " ncol(cdseq_sc_comb_umap) = ",ncol(cdseq_sc_comb_umap),"\n")
      }
      
      cell_sources[CDSeq_idx] <- "CDSeq"
      cell_sources[scRNA_idx] <- "scRNA"
      cell_sources <- factor(cell_sources,levels = c('scRNA','CDSeq') )
      #cell_sources <- factor(c(rep('scRNA',length(scRNA_idx)), rep('CDSeq',length(CDSeq_idx))), levels = c('scRNA','CDSeq'))
      cdseq_sc_comb_umap_per_sample_df <- data.frame(V1 = cdseq_sc_comb_umap[,1], V2 = cdseq_sc_comb_umap[,2], cell_sources= cell_sources , cluster_label = cluster_label)
      point_stroke <- c(rep(0,length(scRNA_idx)), rep(2,length(CDSeq_idx)))
      if(verbose){
        cat("nrow(cdseq_sc_comb_umap_per_sample_df) = ", nrow(cdseq_sc_comb_umap_per_sample_df), "length(point_stroke) = ", length(point_stroke),"...\n")
      }
      
      cdseq_scRNA_umap_per_sample<- ggplot(cdseq_sc_comb_umap_per_sample_df,aes(x=.data$V1, y=.data$V2, colour=as.factor(cell_sources), size=as.factor(cell_sources), shape=as.factor(cell_sources), fill=as.factor(cluster_label), stroke=point_stroke)) + 
        ggtitle(paste0('CDSeq estimated per-sample-cell-type-specific GEP and scRNAseq')) + 
        xlab("UMAP 1") + ylab("UMAP 2") +
        geom_point() + # this is the edge color
        scale_colour_manual(name="cell source",values=c("gray","black"),na.value = NA,label=levels(cell_sources)) +
        scale_fill_manual(name="clusters",values=c(rainbow(length(unique(cluster_label)))),na.value = NA,label=levels(cluster_label)) +
        scale_size_manual(name="cell source", values = c(sc_pt_size,cdseq_pt_size), na.value = NA,label=levels(cell_sources)) +
        scale_shape_manual(name="cell source", values = c(21,23), na.value = NA,label=levels(cell_sources)) +
        theme(plot.title = element_text(size = 45, face = "bold"), 
              axis.title.x = element_text(color="black", size=35,face="bold"),
              axis.title.y = element_text(color="black", size=35,face="bold"),
              legend.title = element_text(color = "blue", face = "bold", size = 30),
              legend.text = element_text(color = "blue", face = "bold", size = 30)) + 
        guides(color = guide_legend(override.aes = list(size=5)), fill = guide_legend(override.aes = list(shape=21,size=5)))
      if(fig_save){
        #fig_tmp_name <- paste0(fig_path,fig_name,"_umap_per_sample_cts_",gsub("-|\\s|:","_",Sys.time()),".",fig_format)
        #if(verbose){cat("save umap at ",fig_tmp_name ,"...\n")}
        if(verbose){cat("save umap at ",file.path(fig_path,fig_name_umap_perSample) ,"...\n")}
        ggsave(filename = fig_name_umap_perSample,#fig_tmp_name,#paste0(fig_path,fig_name,"_umap.",fig_format),
               path = fig_path,
               plot = cdseq_scRNA_umap_per_sample,
               width = 25,
               height = 20,
               dpi = fig_dpi)
      }
    }
    
    
    
  }
  ##################################################################
  ##                             TSNE                             ##
  ##################################################################
  if(plot_tsne){
    if(verbose){cat("Running TSNE: synthetic CDSeq-estimated cell types ... \n")}
    cdseq_synth_scRNA_seurat <- Seurat::RunTSNE(cdseq_synth_scRNA_seurat, dims = seurat_dims, reduction = seurat_reduction ,check_duplicates = FALSE, verbose = FALSE)
    
    # plot scRNAseq and cdseq estimates
    if(is.null(sc_annotation)){
      cluster_label <- paste0("cluster_",cdseq_synth_scRNA_seurat$seurat_clusters)#sub_grp#clusterlouvain$membership
    }else{
      tmp_label <- as.numeric(cdseq_synth_scRNA_seurat$seurat_clusters)
      cluster_label <- rep("a",length(tmp_label))
      for (i in 1:length(seurat_unique_clusters)) {
        cluster_label[which(tmp_label==seurat_unique_clusters[i])] <- seurat_cluster_gold_label[i]
      }
    }
    
    cluster_label <- factor(cluster_label, levels = unique(cluster_label))
    cdseq_sc_comb_tsne <- cdseq_synth_scRNA_seurat@reductions$tsne@cell.embeddings
    
    cell_sources <- rownames(cdseq_sc_comb_tsne)
    CDSeq_idx <- grep("CDSeq",rownames(cdseq_sc_comb_tsne))
    cell_sources[CDSeq_idx] <- "CDSeq"
    point_stroke <- c(rep(1,length(CDSeq_idx)))
    
    if(!is.null(sc_gep)){
      tsne_tot <- 1:nrow(cdseq_sc_comb_tsne)
      scRNA_idx <- tsne_tot[-CDSeq_idx]
      cell_sources[scRNA_idx] <- "scRNA"
      cell_sources <- factor(cell_sources,levels = c('scRNA','CDSeq') )
      point_stroke <- c(rep(0,length(scRNA_idx)), rep(2,length(CDSeq_idx)))
    }
    
    
    cdseq_sc_comb_tsne_df <- data.frame(V1 = cdseq_sc_comb_tsne[,1], V2 = cdseq_sc_comb_tsne[,2], cell_sources= cell_sources , cluster_label = cluster_label)
    
    if(!is.null(sc_gep)){
      cdseq_scRNA_tsne<- ggplot(cdseq_sc_comb_tsne_df,aes(x=.data$V1, y=.data$V2, colour=as.factor(cell_sources), size=as.factor(cell_sources), shape=as.factor(cell_sources), fill=as.factor(cluster_label), stroke=point_stroke)) +
        ggtitle(paste0('CDSeq estimated cell types and scRNAseq')) + 
        xlab("TSNE 1") + ylab("TSNE 2") +
        geom_point() + # this is the edge color
        scale_colour_manual(name="cell source",values=c("gray","black"),na.value = NA,label=levels(cell_sources)) +
        scale_fill_manual(name="clusters",values=c(rainbow(length(unique(cluster_label)))),na.value = NA,label=levels(cluster_label)) +
        scale_size_manual(name="cell source", values = c(sc_pt_size,cdseq_pt_size), na.value = NA,label=levels(cell_sources)) +
        scale_shape_manual(name="cell source", values = c(21,23), na.value = NA,label=levels(cell_sources)) +
        theme(plot.title = element_text(size = 45, face = "bold"), 
              axis.title.x = element_text(color="black", size=35,face="bold"),
              axis.title.y = element_text(color="black", size=35,face="bold"),
              legend.title = element_text(color = "blue", face = "bold", size = 30),
              legend.text = element_text(color = "blue", face = "bold", size = 30)) + 
        guides(color = guide_legend(override.aes = list(size=5)), fill = guide_legend(override.aes = list(shape=21,size=5)))
    }else{
      cdseq_scRNA_tsne<- ggplot(cdseq_sc_comb_tsne_df,aes(x=.data$V1, y=.data$V2, fill=as.factor(cluster_label), stroke=point_stroke)) +
        ggtitle(paste0('CDSeq estimated cell types')) + 
        xlab("TSNE 1") + ylab("TSNE 2") +
        geom_point(size = cdseq_pt_size, shape=21) + # this is the edge color
        #scale_colour_manual(name="cell source",values=c("gray","black"),na.value = NA,label=levels(cell_sources)) +
        scale_fill_manual(name="clusters",values=c(rainbow(length(unique(cluster_label)))),na.value = NA,label=levels(cluster_label)) +
        #scale_size_manual(name="cell source", values = c(sc_pt_size,cdseq_pt_size), na.value = NA,label=levels(cell_sources)) +
        #scale_shape_manual(name="cell source", values = c(21,23), na.value = NA,label=levels(cell_sources)) +
        theme(plot.title = element_text(size = 45, face = "bold"), 
              axis.title.x = element_text(color="black", size=35,face="bold"),
              axis.title.y = element_text(color="black", size=35,face="bold"),
              legend.title = element_text(color = "blue", face = "bold", size = 30),
              legend.text = element_text(color = "blue", face = "bold", size = 30)) + 
        guides(color = guide_legend(override.aes = list(size=5)), fill = guide_legend(override.aes = list(shape=21,size=5)))
    }
    
    
    if(fig_save){
      #fig_tmp_name <- paste0(fig_path,fig_name,"_tsne_",gsub("-|\\s|:","_",Sys.time()),".",fig_format)
      #if(verbose){cat("save umap at ", fig_tmp_name ,"...\n")}
      if(verbose){cat("save umap at ", file.path(fig_path,fig_name_tsne) ,"...\n")}
      ggsave(filename = fig_name_tsne,
             path = fig_path,
             plot = cdseq_scRNA_tsne,
             width = 25,
             height = 20,
             dpi = fig_dpi)
    }

    
    ##################################################################
    ##    plot CDSeq-estimated cell-type-specific gep per sample    ##
    ##################################################################
    if(plot_per_sample){
      warning("plot_per_sample option is not working properly. I'm working on it.")
      if(verbose){cat("Running TSNE: per sample CDSeq-estimated cell types ... \n")}
      cdseq_synth_scRNA_seurat_persample <- Seurat::RunTSNE(cdseq_synth_scRNA_seurat_persample, dims = seurat_dims, reduction = seurat_reduction ,check_duplicates = FALSE,verbose = FALSE)
      
      # plot scRNAseq and cdseq estimates
      if(is.null(sc_annotation)){
        cluster_label <- paste0("cluster_",cdseq_synth_scRNA_seurat_persample$seurat_clusters)#sub_grp#clusterlouvain$membership
      }else{
        tmp_label <- as.numeric(cdseq_synth_scRNA_seurat_persample$seurat_clusters)
        cluster_label <- rep("a",length(tmp_label))
        for (i in 1:length(seurat_unique_clusters_persample)) {
          cluster_label[which(tmp_label==seurat_unique_clusters_persample[i])] <- seurat_cluster_gold_label_persample[i]
        }
      }
      
      cluster_label <- factor(cluster_label, levels = unique(cluster_label))
      cdseq_sc_comb_tsne <- cdseq_synth_scRNA_seurat_persample@reductions$tsne@cell.embeddings
      
      cell_sources <- rownames(cdseq_sc_comb_tsne)
      tsne_tot <- 1:nrow(cdseq_sc_comb_tsne)
      CDSeq_idx <- grep("CDSeq",rownames(cdseq_sc_comb_tsne))
      scRNA_idx <- tsne_tot[-CDSeq_idx]
      
      cell_sources[CDSeq_idx] <- "CDSeq"
      cell_sources[scRNA_idx] <- "scRNA"
      cell_sources <- factor(cell_sources,levels = c('scRNA','CDSeq') )
      
      cdseq_sc_comb_tsne_df_persample <- data.frame(V1 = cdseq_sc_comb_tsne[,1], V2 = cdseq_sc_comb_tsne[,2], cell_sources= cell_sources , cluster_label = cluster_label)
      point_stroke <- c(rep(0,length(scRNA_idx)), rep(2,length(CDSeq_idx)))
      
      cdseq_scRNA_tsne_persample <- ggplot(cdseq_sc_comb_tsne_df_persample,aes(x=.data$V1, y=.data$V2, colour=as.factor(cell_sources), size=as.factor(cell_sources), shape=as.factor(cell_sources), fill=as.factor(cluster_label), stroke=point_stroke)) +
        ggtitle(paste0('CDSeq estimated cell types and scRNAseq')) + 
        xlab("TSNE 1") + ylab("TSNE 2") +
        geom_point() + # this is the edge color
        scale_colour_manual(name="cell source",values=c("gray","black"),na.value = NA,label=levels(cell_sources)) +
        scale_fill_manual(name="clusters",values=c(rainbow(length(unique(cluster_label)))),na.value = NA,label=levels(cluster_label)) +
        scale_size_manual(name="cell source", values = c(sc_pt_size,cdseq_pt_size), na.value = NA,label=levels(cell_sources)) +
        scale_shape_manual(name="cell source", values = c(21,23), na.value = NA,label=levels(cell_sources)) +
        theme(plot.title = element_text(size = 45, face = "bold"), 
              axis.title.x = element_text(color="black", size=35,face="bold"),
              axis.title.y = element_text(color="black", size=35,face="bold"),
              legend.title = element_text(color = "blue", face = "bold", size = 30),
              legend.text = element_text(color = "blue", face = "bold", size = 30)) + 
        guides(color = guide_legend(override.aes = list(size=5)), fill = guide_legend(override.aes = list(shape=21,size=5)))
      
      if(fig_save){
        #fig_tmp_name <- paste0(fig_path,fig_name,"_tsne_per_sample_cts_",gsub("-|\\s|:","_",Sys.time()),".",fig_format)
        #if(verbose){cat("save umap at ", fig_tmp_name ,"...\n")}
        if(verbose){cat("save umap at ", file.path(fig_path,fig_name_tsne_perSample),"...\n")}
        ggsave(filename = fig_name_tsne_perSample,#paste0(fig_path,fig_name,"_tsne.",fig_format),
               path = fig_path,
               plot = cdseq_scRNA_tsne_persample,
               width = 25,
               height = 20,
               dpi = fig_dpi)
      }
      
    }
    
  }
  

  output <- list()
  #output$fig_path <- fig_path
  #output$fig_name <- fig_name
  output$cdseq_synth_scRNA <- cdseq_synth_scRNA
  if(plot_umap){output$cdseq_scRNA_umap <- cdseq_scRNA_umap}
  if(plot_tsne){output$cdseq_scRNA_tsne <- cdseq_scRNA_tsne}
  output$cdseq_synth_scRNA_seurat <- cdseq_synth_scRNA_seurat
  
  if(!is.null(sc_annotation) && !is.null(sc_gep)){
    output$seurat_cluster_purity <- seurat_cluster_purity
    output$seurat_cluster_gold_label <- seurat_cluster_gold_label
    output$seurat_unique_clusters <- seurat_unique_clusters
    output$seurat_markers <- cdseq_synth_scRNA_seurat_markers_df
    output$CDSeq_cell_type_assignment_df <- CDSeq_cell_type_assignment_df
    output$CDSeq_cell_type_assignment_confidence <- CDSeq_cell_type_assignment_confidence
    output$CDSeq_cell_type_assignment_df_all <- CDSeq_cell_type_assignment_df_all
    
    output$cdseq_gep_merged_byCorr <- cdseq_gep_merged_byCorr
    output$cdseq_annotation_byCorr <- cdseq_annotation_byCorr
  }else{
    output$seurat_markers <- cdseq_synth_scRNA_seurat_markers
  }
  
  if(!is.null(cdseq_prop) && !is.null(sc_annotation) && !is.null(sc_gep)){
    output$cdseq_prop_merged <- cdseq_prop_merged
    output$cdseq_prop_merged_byCorr <- cdseq_prop_merged_byCorr
    # cdseq_prop_merged_byCorr
    # cdseq_gep_merged_byCorr
    # cdseq_annotation_byCorr
    }
  if(!is.null(cdseq_gep_sample_specific)){output$cdseq_gep_sample_specific_merged <- cdseq_gep_sample_specific_merged}
  if(seurat_find_marker){
    output$cdseq_synth_scRNA_seurat_markers <- cdseq_synth_scRNA_seurat_markers
    output$seurat_top_markers <- seurat_top_markers#seurat_top_markers_df
  }
  output$input_list <- input_list
  if(plot_umap){output$cdseq_sc_comb_umap_df <- cdseq_sc_comb_umap_df}
  if(plot_tsne){output$cdseq_sc_comb_tsne_df <- cdseq_sc_comb_tsne_df}
  return(output)
  
  
  
  
  
  
  
  
  
}