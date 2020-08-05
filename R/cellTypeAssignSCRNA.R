#' \code{cellTypeAssignSCRNA} assigns CDSeq-identified cell types using single cell RNAseq data.
#' @param cdseq_gep CDSeq-estimated gene expression profile matrix with G rows (genes) and T columns (cell types).
#' @param cdseq_prop CDSeq-estimated sample-specific cell-type proportion, a matrix with T rows (cell type) and M (sample size).
#' @param cdseq_gep_sample_specific CDSeq-estimated sample-specific cell type gene expression, in the form of read counts. It is a 3 dimension array, i.e. gene by sample by cell type. The element cdseq_gep_sample_specific[i,j,k] represents the reads mapped to gene i from cell type k in sample j.
#' @param sc_gep a G (genes) by N (cell) matrix or dataframe that contains the gene expression profile for N single cells.
#' @param sc_annotation a dataframe contains two columns "cell_id"  and "cell_type". cell_id needs to match with the cell_id in sc_gep but not required to have the same size. cell_type is the cell type annotation for the single cells.
#' @param nb_size size parameter for negative binomial distribution, check rnbinom for details.
#' @param nb_mu mu parameter for negative binomial distribution, check rnbinom for details.
#' @param seurat_count_threshold this parameter will be passed to Seurat subset function (subset = nCount_RNA > seurat_count_threshold) for filtering out single cells whose total counts is less this threshold. 
#' @param seurat_scale_factor this parameter will be passed to scale.factor in Seurat function NormalizeData.
#' @param seurat_norm_method this parameter will be passed to normalization.method in Seurat function NormalizeData.
#' @param seurat_select_method this parameter will be passed to selection.method in Seurat function FindVariableFeatures
#' @param seurat_nfeatures this parameter will be passed to nfeatures in Seurat function FindVariableFeatures.
#' @param seurat_npcs this parameter will be passed to npcs in Seurat function RunPCA.
#' @param seurat_dims this parameter will be passed to dims in Seurat function FindNeighbors.
#' @param seurat_reduction this parameter will be passed to reduction in Seurat function FindNeighbors.
#' @param seurat_resolution ths parameter will be passed to resolution in Seurat function FindClusters.
#' @param seurat_DE_test this parameter will be passed to test.use in Seurat function FindAllMarkers.
#' @param seurat_DE_logfc this parameter will be passed to logfc.threshold in Seurat function FindAllMarkers.
#' @param seurat_top_n_markers the number of top DE markers saved from Seurat output. 
#' @param plot_umap set 1 to plot umap figure of scRNAseq and CDSeq-estimated cell types, 0 otherwise.
#' @param plot_tsne set 1 to plot tsne figure of scRNAseq and CDSeq-estimated cell types, 0 otherwise.
#' @param fig_path the location where the heatmap figure is saved. 
#' @param fig_name the name of umap and tsne figures. Umap figure will have the name of fig_name_umap_date and tsne figure will be named fig_name_tsne_date.
#' @param fig_format "pdf", "jpeg", or "png".
#' @param fig_dpi figure dpi
#' @importFrom grDevices rainbow
#' @importFrom stats rmultinom rnbinom
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures ScaleData FindNeighbors FindClusters FindAllMarkers RunUMAP RunTSNE RunPCA 
#' @importFrom ggplot2 guide_legend guides aes scale_size_manual scale_shape_manual scale_fill_manual xlab ylab theme ggsave ggplot ggtitle geom_point element_text scale_colour_manual 
#' @importFrom dplyr top_n group_by 
#' @importFrom magrittr %>%
#' @importFrom rlang .data
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
#' cdseq_prop_merged: CDSeq-estimated cell type proportions with cell type annotations.
#' 
#' cdseq_gep_sample_specific_merged: sample-specific cell-type read counts. It is a 3d array with dimensions: gene, sample, cell type. 


cellTypeAssignSCRNA <- function(cdseq_gep = NULL,
                                cdseq_prop = NULL,
                                cdseq_gep_sample_specific = NULL,
                                sc_gep = NULL,
                                sc_annotation = NULL,
                                nb_size = 100,
                                nb_mu = 1000,
                                seurat_count_threshold = 100,
                                seurat_scale_factor = 10000,
                                seurat_norm_method = "LogNormalize",
                                seurat_select_method = 'vst',
                                seurat_nfeatures = 100,
                                seurat_npcs = 50,
                                seurat_dims = 1:10,
                                seurat_reduction = 'pca',
                                seurat_resolution = 0.8,
                                seurat_DE_test = "wilcox",
                                seurat_DE_logfc = 0.25,
                                seurat_top_n_markers = 10,
                                plot_umap = 1,
                                plot_tsne = 1,
                                fig_path = getwd(),
                                fig_name = "cellTypeAssignSCRNA",
                                fig_format = "pdf",
                                fig_dpi = 300){
  #################################################################
  ##                         check input                         ##
  #################################################################
  cat("checking input ... \n")
  # check dimension
  nr_cdseq <- nrow(cdseq_gep)
  nc_cdseq <- ncol(cdseq_gep)
  nr_sc <- nrow(sc_gep)
  nc_sc <- ncol(sc_gep)
  
  
  # check gene names
  gene_cdseq <- rownames(cdseq_gep)
  gene_sc <- rownames(sc_gep)
  if(is.null(gene_cdseq) || is.null(gene_sc)){stop("cdseq_gep and sc_gep should have gene names as row names.")}
  if(sum(gene_cdseq == gene_sc) != nr_cdseq ){stop("cdseq_gep and sc_gep should have the same gene names and in the same order.")}
  
  # check cell type proportions
  if(!is.null(cdseq_prop)){
    nr_prop <- nrow(cdseq_prop)
    nc_prop <- ncol(cdseq_prop)
    if(nr_cdseq!=nr_sc){stop("cdseq_gep and sc_gep should have the same number of rows (genes).")}
    if(nr_prop!=nc_cdseq){stop("column number of cdseq_gep should equal to row number of cdseq_prop.")}
    
    if(is.null(colnames(cdseq_prop))){stop("cdseq_prop column name is missing. Column names are sample ids.")}
    if(is.null(rownames(cdseq_prop))){stop("cdseq_prop row name is missing. Row names are the cell type names.")}
    if(sum(rownames(cdseq_prop) == colnames(cdseq_gep)) != nr_prop){stop("Row names of cdseq_prop and column names of cdseq_gep should be the same and in same order. They both denote cell types.")}
  }
  
  # check cell ids
  sc_id <- colnames(sc_gep)
  if(is.null(sc_id)){stop("Please provide cell id in sc_gep.")}
  
  # check annotation
  if(is.null(sc_annotation)){
    warning("sc_annotation (single cell annotation) is not provided. Cluster labels and DE genes will be generated. It is recommanded to use annotated scRNAseq in this function.")
  }else{
    if(is.null(dim(sc_annotation))){stop("sc_annotation has to be a two dimensional data frame with column names: cell_id and cell_type")}
    if(nrow(sc_annotation)!=nc_sc){stop("nrow(sc_annotation) should be equal to ncol(sc_gep).")}
    if(sum(colnames(sc_annotation) == c("cell_id","cell_type"))!=2){stop("sc_annotation has to be a two dimensional data frame with column names: cell_id and cell_type")}
    #cat("sum(sc_id == sc_annotation$cell_id) = ",sum(sc_id == sc_annotation$cell_id),"length(sc_id)=",length(sc_id),"\n")
    if(sum(sc_id == sc_annotation$cell_id) != length(sc_id)){stop("cell id in sc_gep has to be the same as the cell id in sc_annotation.")}
  }
  
  
  # declare variable to avoid R CMD check notes for no visible binding for global variable
  avg_logFC <- cluster <- nCount_RNA <- NULL
  ###################################################################
  ##  Generate synthetic scRNAseq data from CDSeq estimated GEPs   ##
  ##                 and combine with scRNAseq data                ##
  ###################################################################
  cat("generating synthetic scRNAseq data using CDSeq estimated GEPs ...\n")
  # generate synthetic scRNAseq 
  ncell <-  1
  cdseq_synth_scRNA <- matrix(0,nrow = nr_cdseq, ncol = nc_cdseq*ncell)
  for (i in 1:nc_cdseq) {
    nreads <- rnbinom(n = ncell,size = nb_size, mu = nb_mu)
    for (j in 1:ncell) {
      cdseq_synth_scRNA[,(i-1)*ncell + j] <-  rmultinom(1,nreads[j],cdseq_gep[,i])
    }
  }
  #colnames(cdseq_synth_scRNA) <- paste("CDSeq_SynthCell",rep(1:nc_cdseq,each = ncell), rep(1:ncell), sep = ".")  
  colnames(cdseq_synth_scRNA) <- paste(rep(colnames(cdseq_gep),each = ncell), "CDSeq" ,rep(1:ncell), sep = ".")  
  rownames(cdseq_synth_scRNA) <- rownames(cdseq_gep)
  
  cdseq_sc_comb <- cbind(sc_gep,cdseq_synth_scRNA)
  
  ##################################################################
  ##                          Run seurat                          ##
  ##################################################################
  cat("Calling Seurat pipeline ... \n")
  cdseq_synth_scRNA_seurat <- Seurat::CreateSeuratObject(counts = cdseq_sc_comb, project = "cdseq_synth_scRNAseq")
  # filter cells
  cdseq_synth_scRNA_seurat <- subset(cdseq_synth_scRNA_seurat, subset = nCount_RNA > seurat_count_threshold )#nFeature_RNA > 20 & nFeature_RNA < 2500)
  # normalize
  cdseq_synth_scRNA_seurat <- Seurat::NormalizeData(cdseq_synth_scRNA_seurat, normalization.method = seurat_norm_method , scale.factor = seurat_scale_factor)
  # select genes
  cdseq_synth_scRNA_seurat <- Seurat::FindVariableFeatures(cdseq_synth_scRNA_seurat, selection.method = seurat_select_method, nfeatures = seurat_nfeatures)

  cdseq_synth_scRNA_seurat <- Seurat::ScaleData(cdseq_synth_scRNA_seurat)
  
  cdseq_synth_scRNA_seurat <- Seurat::RunPCA(cdseq_synth_scRNA_seurat, npcs = seurat_npcs) #features = VariableFeatures(object = cdseq_synth_scRNA_seurat), 
  
  cdseq_synth_scRNA_seurat <- Seurat::FindNeighbors(cdseq_synth_scRNA_seurat, dims = seurat_dims, reduction = seurat_reduction)
  
  cdseq_synth_scRNA_seurat <- Seurat::FindClusters(cdseq_synth_scRNA_seurat, resolution = seurat_resolution)
  
  ##################################################################
  ##               use seurat cluster and DE genes                ##
  ##            for CDSeq-estimatd cell type annotation           ##
  ##################################################################
  # find makers
  cdseq_synth_scRNA_seurat_markers <- Seurat::FindAllMarkers(cdseq_synth_scRNA_seurat, min.pct = 0.25, logfc.threshold = seurat_DE_logfc, test.use = seurat_DE_test)
  cdseq_synth_scRNA_seurat_markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 2, wt = avg_logFC)
  seurat_top_markers <- cdseq_synth_scRNA_seurat_markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = seurat_top_n_markers, wt = avg_logFC)
  
  #seurat_top_markers_df <- as.data.frame(cbind(seurat_top_markers, singleCellAnnotation = rep("NA",nrow(seurat_top_markers))))
  
  if(!is.null(sc_annotation)){
    ##################################################################
    ##                  use scRNAseq annotation to                  ##
    ##              annotate CDSeq estimated cell types             ##
    ##################################################################
    cdseq_synth_scRNA_seurat_markers_df <- data.frame(cdseq_synth_scRNA_seurat_markers, singleCellAnnotation = rep("NA",nrow(cdseq_synth_scRNA_seurat_markers)),stringsAsFactors = FALSE)
    
    # grep single cell id in seurat output
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
    
    # construct seurat cluster label 
    # NOTICE THAT [scRNAseq_cell_cluster_label_seurat] and [scRNAseq_seurat_goldstandard_annotation] should be in the same order in term of cell ids
    scRNAseq_cell_cluster_label_seurat <- cdseq_synth_scRNA_seurat@meta.data$seurat_clusters[seurat_scRNA_comm_idx]
    
    cdseq_cell_idx <- grep("CDSeq",scRNAseq_cell_id_seurat)
    cdseq_cell_seurat_cluster_id <- cdseq_synth_scRNA_seurat@meta.data$seurat_clusters[cdseq_cell_idx]
    cdseq_cell_id <- scRNAseq_cell_id_seurat[cdseq_cell_idx]
    
    CDSeq_cell_type_assignment_df <- data.frame(cdseq_cell_id = cdseq_cell_id, seurat_cluster = cdseq_cell_seurat_cluster_id, cdseq_cell_type_assignment = rep("NA",length(cdseq_cell_id)),stringsAsFactors = FALSE)
    
    # compute the cluster score (there is a bug due to the large number of cells )
    #seurat_cluster_scores <- cluster_scores(c = as.numeric(scRNAseq_cell_cluster_label_seurat) , k = as.numeric(scRNAseq_seurat_goldstandard_annotation) ,beta = 1)
    
    # assign CDSeq-estimated cell types based on single cell annotation
    seurat_unique_clusters <- sort(as.numeric(unique(scRNAseq_cell_cluster_label_seurat)))
    seurat_cluster_purity <- rep(0, length(seurat_unique_clusters))
    seurat_cluster_gold_label <- rep("a",length(seurat_unique_clusters))
    for (i in 1:length(seurat_unique_clusters)) {
      seurat_cluster_members <- which(as.numeric(scRNAseq_cell_cluster_label_seurat) == seurat_unique_clusters[i])
      seurat_cluster_gold_standard_label <- as.character(scRNAseq_seurat_goldstandard_annotation[seurat_cluster_members])
      seurat_cluster_assess <- max_rep(seurat_cluster_gold_standard_label)
      seurat_cluster_purity[i] <- seurat_cluster_assess$max_element_proportion
      seurat_cluster_gold_label[i] <- seurat_cluster_assess$max_element
      # add sc_annotation to DE marker data frame
      tmp_idx <- which(as.numeric(cdseq_synth_scRNA_seurat_markers_df$cluster) == seurat_unique_clusters[i])
      cdseq_synth_scRNA_seurat_markers_df$singleCellAnnotation[tmp_idx] <- seurat_cluster_gold_label[i]
      
      cdseq_tmp_idx <- which(as.numeric(CDSeq_cell_type_assignment_df$seurat_cluster) == seurat_unique_clusters[i])
      CDSeq_cell_type_assignment_df$cdseq_cell_type_assignment[cdseq_tmp_idx] <- seurat_cluster_gold_label[i]
    }
    seurat_top_markers_df <- cdseq_synth_scRNA_seurat_markers_df %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = seurat_top_n_markers, wt = avg_logFC)
    
    if(!is.null(cdseq_prop)){
       # current version assumes each CDSeq-estGEP only generate 1 synthetic single cell. May extend to more general case later
       # the gsub function extract patterns between two dots.
       #cdseq_est_cell_ids <- as.numeric(unlist(lapply(CDSeq_cell_type_assignment_df$cdseq_cell_id, function(x){gsub("^.+?\\.(.+?)\\..*$", "\\1",x)})))
      
       # the sub function extract pattern before a dot
       cdseq_est_cell_ids <- unlist(lapply(CDSeq_cell_type_assignment_df$cdseq_cell_id, function(x){sub("\\..*","",x)}))
       if(sum(cdseq_est_cell_ids == rownames(cdseq_prop)) != nr_prop){stop("Please make sure the rownames(cdseq_prop) and colnames(cdseq_gep) are the same and in same order.")}
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
      }else{
        cdseq_est_cell_ids <- unlist(lapply(CDSeq_cell_type_assignment_df$cdseq_cell_id, function(x){sub("\\..*","",x)}))
        if(sum(cdseq_est_cell_ids == dimnames(cdseq_gep_sample_specific)[[3]]) != nr_prop){stop("Please make sure the dimnames(cdseq_gep_sample_specific)[[3]] and colnames(cdseq_gep) are the same and in same order.")}
        dimnames(cdseq_gep_sample_specific)[[3]] <- CDSeq_cell_type_assignment_df$cdseq_cell_type_assignment
        for (i in 1:dim(cdseq_gep_sample_specific)[2]) {
          cdseq_gep_sample_specific_merged[,i,] <- t(rowsum(t(cdseq_gep_sample_specific[,i,]), group = rownames(t(cdseq_gep_sample_specific[,i,]))))
        }
      }
    }
  }

  
  ##################################################################
  ##                             UMAP                             ##
  ##################################################################
  if(plot_umap){
    cat("Running UMAP ... \n")
    cdseq_synth_scRNA_seurat <- RunUMAP(cdseq_synth_scRNA_seurat, dims = seurat_dims)
    
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
    cdseq_sc_comb_umap <- cdseq_synth_scRNA_seurat@reductions$umap@cell.embeddings
    
    cell_sources <- rownames(cdseq_sc_comb_umap)
    umap_tot <- 1:nrow(cdseq_sc_comb_umap)
    CDSeq_idx <- grep("CDSeq",rownames(cdseq_sc_comb_umap))
    scRNA_idx <- umap_tot[-CDSeq_idx]
    
    cell_sources[CDSeq_idx] <- "CDSeq"
    cell_sources[scRNA_idx] <- "scRNA"
    cell_sources <- factor(cell_sources,levels = c('scRNA','CDSeq') )
    #cell_sources <- factor(c(rep('scRNA',length(scRNA_idx)), rep('CDSeq',length(CDSeq_idx))), levels = c('scRNA','CDSeq'))
    cdseq_sc_comb_umap_df <- data.frame(V1 = cdseq_sc_comb_umap[,1], V2 = cdseq_sc_comb_umap[,2], cell_sources= cell_sources , cluster_label = cluster_label)
    point_stroke <- c(rep(0,length(scRNA_idx)), rep(2,length(CDSeq_idx)))
    
    cat("nrow(cdseq_sc_comb_umap_df) = ", nrow(cdseq_sc_comb_umap_df), "length(point_stroke) = ", length(point_stroke),"...\n")
    
    cdseq_scRNA_umap<- ggplot(cdseq_sc_comb_umap_df,aes(x=.data$V1, y=.data$V2, colour=as.factor(cell_sources), size=as.factor(cell_sources), shape=as.factor(cell_sources), fill=as.factor(cluster_label), stroke=point_stroke)) + 
      ggtitle(paste0('CDSeq estimated cell types and scRNAseq')) + 
      xlab("UMAP 1") + ylab("UMAP 2") +
      geom_point() + # this is the edge color
      scale_colour_manual(name="cell source",values=c("gray","black"),na.value = NA,label=levels(cell_sources)) +
      scale_fill_manual(name="clusters",values=c(rainbow(length(unique(cluster_label)))),na.value = NA,label=levels(cluster_label)) +
      scale_size_manual(name="cell source", values = c(1,3), na.value = NA,label=levels(cell_sources)) +
      scale_shape_manual(name="cell source", values = c(21,23), na.value = NA,label=levels(cell_sources)) +
      theme(plot.title = element_text(size = 45, face = "bold"), 
            axis.title.x = element_text(color="black", size=35,face="bold"),
            axis.title.y = element_text(color="black", size=35,face="bold"),
            legend.title = element_text(color = "blue", face = "bold", size = 30),
            legend.text = element_text(color = "blue", face = "bold", size = 30)) + 
      guides(color = guide_legend(override.aes = list(size=5)), fill = guide_legend(override.aes = list(shape=21,size=5)))
    
    fig_tmp_name <- paste0(fig_path,fig_name,"_umap_",gsub("-|\\s|:","_",Sys.time()),".",fig_format)
    cat("save umap at ",fig_tmp_name ,"...\n")
    ggsave(filename = fig_tmp_name,#paste0(fig_path,fig_name,"_umap.",fig_format),
           plot = cdseq_scRNA_umap,
           width = 25,
           height = 20,
           dpi = fig_dpi)
  }
  ##################################################################
  ##                             TSNE                             ##
  ##################################################################
  if(plot_tsne){
    cat("Running TSNE ... \n")
    cdseq_synth_scRNA_seurat <- RunTSNE(cdseq_synth_scRNA_seurat, dims = seurat_dims, check_duplicates = FALSE)
    
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
    tsne_tot <- 1:nrow(cdseq_sc_comb_tsne)
    CDSeq_idx <- grep("CDSeq",rownames(cdseq_sc_comb_tsne))
    scRNA_idx <- tsne_tot[-CDSeq_idx]
    
    cell_sources[CDSeq_idx] <- "CDSeq"
    cell_sources[scRNA_idx] <- "scRNA"
    cell_sources <- factor(cell_sources,levels = c('scRNA','CDSeq') )
    
    cdseq_sc_comb_tsne_df <- data.frame(V1 = cdseq_sc_comb_tsne[,1], V2 = cdseq_sc_comb_tsne[,2], cell_sources= cell_sources , cluster_label = cluster_label)
    point_stroke <- c(rep(0,length(scRNA_idx)), rep(2,length(CDSeq_idx)))
    
    cdseq_scRNA_tsne<- ggplot(cdseq_sc_comb_tsne_df,aes(x=.data$V1, y=.data$V2, colour=as.factor(cell_sources), size=as.factor(cell_sources), shape=as.factor(cell_sources), fill=as.factor(cluster_label), stroke=point_stroke)) +
      ggtitle(paste0('CDSeq estimated cell types and scRNAseq')) + 
      xlab("TSNE 1") + ylab("TSNE 2") +
      geom_point() + # this is the edge color
      scale_colour_manual(name="cell source",values=c("gray","black"),na.value = NA,label=levels(cell_sources)) +
      scale_fill_manual(name="clusters",values=c(rainbow(length(unique(cluster_label)))),na.value = NA,label=levels(cluster_label)) +
      scale_size_manual(name="cell source", values = c(1,3), na.value = NA,label=levels(cell_sources)) +
      scale_shape_manual(name="cell source", values = c(21,23), na.value = NA,label=levels(cell_sources)) +
      theme(plot.title = element_text(size = 45, face = "bold"), 
            axis.title.x = element_text(color="black", size=35,face="bold"),
            axis.title.y = element_text(color="black", size=35,face="bold"),
            legend.title = element_text(color = "blue", face = "bold", size = 30),
            legend.text = element_text(color = "blue", face = "bold", size = 30)) + 
      guides(color = guide_legend(override.aes = list(size=5)), fill = guide_legend(override.aes = list(shape=21,size=5)))
    
    fig_tmp_name <- paste0(fig_path,fig_name,"_tsne_",gsub("-|\\s|:","_",Sys.time()),".",fig_format)
    cat("save umap at ", fig_tmp_name ,"...\n")
    ggsave(filename = fig_tmp_name,#paste0(fig_path,fig_name,"_tsne.",fig_format),
           plot = cdseq_scRNA_tsne,
           width = 25,
           height = 20,
           dpi = fig_dpi)
  }

  output <- list()
  output$fig_path <- fig_path
  output$fig_name <- fig_name
  output$cdseq_synth_scRNA <- cdseq_synth_scRNA
  if(plot_umap){output$cdseq_scRNA_umap <- cdseq_scRNA_umap}
  if(plot_tsne){output$cdseq_scRNA_tsne <- cdseq_scRNA_tsne}
  output$cdseq_synth_scRNA_seurat <- cdseq_synth_scRNA_seurat
  if(!is.null(sc_annotation)){
    output$seurat_cluster_purity <- seurat_cluster_purity
    output$seurat_cluster_gold_label <- seurat_cluster_gold_label
    output$seurat_unique_clusters <- seurat_unique_clusters
    output$seurat_markers <- cdseq_synth_scRNA_seurat_markers_df
    output$CDSeq_cell_type_assignment_df <- CDSeq_cell_type_assignment_df
  }else{
    output$seurat_markers <- cdseq_synth_scRNA_seurat_markers
  }
  if(!is.null(cdseq_prop)){output$cdseq_prop_merged <- cdseq_prop_merged}
  if(!is.null(cdseq_gep_sample_specific)){output$cdseq_gep_sample_specific_merged <- cdseq_gep_sample_specific_merged}
  output$seurat_top_markers <- seurat_top_markers_df
  return(output)
  
  
  
  
  
  
  
  
  
}