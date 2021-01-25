
result1<-CDSeq(bulk_data =  mixtureGEP, cell_type_number = 6, mcmc_iterations = 5, 
        dilution_factor = 50, block_number = 1, gene_length = as.vector(gene_length), 
        reference_gep = refGEP, cpu_number=1, print_progress_msg_to_file=0)

test_that("GEP and proportion estimates are returned", {
  expect_true("estGEP" %in% names(result1))
  expect_true("estProp" %in% names(result1))
})
#> Test passed 

test_that("GEP and proportion are in the right dimension", {
  expect_equal(nrow(result1[["estGEP"]]), nrow(mixtureGEP))
  expect_equal(ncol(result1[["estGEP"]]),6)
  expect_equal(nrow(result1[["estProp"]]),6)
  expect_equal(ncol(result1[["estProp"]]), ncol(mixtureGEP))
})
#> Test passed 

result2<-CDSeq(bulk_data =  mixtureGEP, cell_type_number = 2:3, mcmc_iterations = 5, 
               dilution_factor = 50, block_number = 1, gene_length = as.vector(gene_length), 
               reference_gep = refGEP, cpu_number=1, print_progress_msg_to_file=0)

test_that("Check est_all are returned when cell_type_number is a vector", {
  expect_true("est_all" %in% names(result2))
})
#> Test passed 