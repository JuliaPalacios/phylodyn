
#' Generation of A-matrix
#' 
#' @description This function generates the A-matrix that encodes the labeling information of an input tree.
#' @param t A \code{phylo} object.
#' @return A-matrix of the input tree.
#' @examples 
#' t <- ape::rcoal(5)
#' Amat <- gen_Amat(t)
#' @export
gen_Amat <- function(t) { 
  n_tip = t$Nnode + 1 
  edges = t$edge
  edge_lengths = t$edge.length
  
  # initialize
  edge_mat <- as.data.frame(cbind(edges, edge_lengths))
  colnames(edge_mat) <- c("parent", "child", "edge_length")
  leaves <- seq(1,n_tip)
  heights <- rep(0, n_tip)
  leaf_mat <- as.data.frame(cbind(leaves, heights))
  A <- matrix(0, n_tip-1, n_tip) 
  internal_mat <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(internal_mat) <- c("parent", "height_parent", "child1", "child2")
  
  while(nrow(edge_mat) > 0){
    # find the parents whose 2 children are both leaves
    filter_mat <- edge_mat[edge_mat$child %in% leaf_mat$leaves,]
    for(p in unique(filter_mat$parent)){
      if(nrow(filter_mat[filter_mat$parent == p,]) < 2){
        filter_mat <- filter_mat[filter_mat$parent != p,]
      }
    }
    filter_mat$height_child <- rep(0, nrow(filter_mat))
    filter_mat$height_parent <- rep(0, nrow(filter_mat))
    for(i in seq(1:nrow(filter_mat))){
      filter_mat$height_child[i] <- leaf_mat[leaf_mat$leaves == filter_mat$child[i],]$heights
    }
    filter_mat$height_parent <- filter_mat$height_child + filter_mat$edge_length
    
    # update leaf_mat, internal_mat, edge_mat (remove edges)
    for(p in unique(filter_mat$parent)){
      temp_mat <- filter_mat[filter_mat$parent == p,]
      leaf_mat <- rbind(leaf_mat, c(temp_mat$parent[1], temp_mat$height_parent[1]))
      internal_mat <- rbind(internal_mat, c(temp_mat$parent[1], temp_mat$height_parent[1], temp_mat$child[1], temp_mat$child[2]))
      edge_mat <- edge_mat[edge_mat$parent != p,]
    }
  }
  
  colnames(internal_mat) <- c("parent", "height_parent", "child1", "child2")
  internal_mat <- internal_mat[order(internal_mat$height_parent),] # order by height_parent
  internal_mat$row <- seq(n_tip-1, 1, -1)
  
  # update A matrix
  for(i in 1:nrow(internal_mat)){
    children <- c(internal_mat$child1[i], internal_mat$child2[i])
    cur_row <- internal_mat$row[i]
    for(c in children){
      if(c %in% seq(1:n_tip)){
        A[cur_row, c] <- 1
      }else{
        A[cur_row,] <- A[cur_row,] + A[internal_mat[internal_mat$parent == c,]$row,]
      }
    }
  }
  return(A)
}

#' Distance between labeled trees
#' 
#' @description This function generates the distance between two labeled trees with the same number of tips. 
#' It is calculated as a weighted average of F-matrix norm and A-matrix norm.
#' @param t1 A \code{phylo} object.
#' @param t2 A \code{phylo} object, with the same number of tips as t1.
#' @param alpha A number between 0 and 1. The weight that we put on F-matrix.
#' @return The distance between two labeled trees.
#' @examples 
#' t1 <- ape::rcoal(5)
#' t2 <- ape::rcoal(5)
#' dist <- label_dist(t1, t2)
#' @export
label_dist <- function(t1, t2, alpha) {
  F1 <- gen_Fmat(t1)
  F2 <- gen_Fmat(t2)
  A1 <- gen_Amat(t1)
  A2 <- gen_Amat(t2)
  return(alpha * norm(F1-F2, type = "F") + (1-alpha) * norm(A1-A2, type = "F"))
}


