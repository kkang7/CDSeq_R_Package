#' \code{intersection} take intersection of multiple lists and return the common set and index
#' @param list.vector this is a list of list contain all the data.
#' @param order this is either sort or stable. If choose sort, the output common value will be sorted. If choose stable, the output common value will be in the same order as appear in the first element in list.vector.
#' @return The common values among lists and their indices.





#' intersection function: 
#' input: list.vector is a list of list contain all the data
#' for example, if we need to find the common elements of a, b, c, then input should be list(a,b,c)

intersection <- function(list.vector,
                         order='sort'
                         ){
  if(order!='sort' && order!='stable'){stop('order has to be either \'sort\' or \'stable\'.')}
  list.length <- length(list.vector)
  if(list.length==1){stop('list.vector has to be a list of list, for example, list(a,b,c)')}
  
  idx <- vector("list",list.length)
  list.comm <- Reduce(intersect,list.vector)
  

  
  if(order=='sort'){
    list.comm <- sort(list.comm)
    for (i in 1:list.length) {
      idx[[i]] <- match(list.comm,list.vector[[i]])
    }
  }
  
  if(order=='stable'){
    for (i in 1:list.length) {
      idx[[i]] <- match(list.comm,list.vector[[i]])
    }
  }
  
  interaction.return <- list(comm.value=list.comm, index=idx)
  return(interaction.return)
}