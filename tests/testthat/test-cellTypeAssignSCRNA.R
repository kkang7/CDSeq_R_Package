cdseq.result <- CDSeq::CDSeq(bulk_data = pbmc_mix,
                             cell_type_number = 3,
                             mcmc_iterations = 5,
                             cpu_number = 1,
                             gene_length = rep(1,100),
                             dilution_factor = 50,
                             reference_gep = sc_gep)

cdseq.result.celltypeassign <- cellTypeAssignSCRNA(cdseq_gep = cdseq.result$estGEP, # CDSeq-estimated cell-type-specific GEPs
                                                   cdseq_prop = cdseq.result$estProp, # CDSeq-estimated cell type proportions
                                                   sc_gep = sc_gep,         # PBMC single cell data
                                                   sc_annotation = sc_annotation,# PBMC single data annotations
                                                   sc_pt_size = 3,
                                                   cdseq_pt_size = 6,
                                                   seurat_nfeatures = 10,
                                                   seurat_npcs = 3,
                                                   seurat_dims=1:3,
                                                   plot_umap = 0,
                                                   plot_tsne = 0)

test_that("GEP and proportion assignments are returned", {
  expect_true("CDSeq_cell_type_assignment_df" %in% names(cdseq.result.celltypeassign))
  expect_true("cdseq_prop_merged" %in% names(cdseq.result.celltypeassign))
})
