#' \code{max_rep} Find the element that repeats the most in a given vector and calculate its proportion.
#' @param v a vector 
#' @return max_rep_value contains two elements: max_element and max_element_proportion. max_element is the element that repeats the most in v, and max_element_proportion
#' is its proportion.

# find the number of maximum repititions

max_rep <- function(v){
  unique_v <- unique(v)
  n_reps <- rep(0,length(unique_v))
  max_rep_value <- list()
  for (i in 1:length(unique_v)) {
    idx <- which(v == unique_v[i])
    n_reps[i] <- length(idx)
  }
  max_element_idx <- which(n_reps==max(n_reps))
  max_element <- unique_v[max_element_idx]
  max_rep_value$max_element <- max_element
  max_rep_value$max_element_proportion <- max(n_reps)/sum(n_reps)
  return(max_rep_value)
}