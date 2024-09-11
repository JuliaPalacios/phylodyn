
# Functions to calculate coalescent rates under Beta and Lambda coalescents 

# Calculate (b choose k) lambda_{b,k}(alpha) under a Beta-coalescent 
beta_coal_factor <- function(b, k, alpha) {
  b <- as.integer(b)
  k <- as.integer(k)
  if (b == 1) {
    return(1)
  }
  if (alpha < 0 | alpha >2 ) {
    return(NA)
  } else if (alpha == 2) {
    if (k > 2) {
      return(0)
    } else if (k == 2) {
      return(1)
    }
  }
  return(choose(b, k)*beta(k - alpha, b-k+alpha)/beta(2-alpha, alpha))
}

# Calculate (b choose k) lambda_{b,k}(alpha) for k=2,...,b 
beta_coal_factor_rates <- function(b, alpha) {
  b <- as.integer(b)
  if (b==1) {
    return(1)
  }
  return(sapply(2:b, beta_coal_factor, b=b, alpha=alpha))
}

# Calculate lambda_b(alpha)
beta_coal_factor_total <- function(b, alpha) {
  b <- as.integer(b)
  if (b==1) {
    return(1)
  }
  return(sum(sapply(2:b, beta_coal_factor, b=b, alpha=alpha)))
}

# Calculate ratio of (b choose k) lambda_{b,k}/lambda_b
beta_coal_factor_prop <- function(b, alpha) {
  b <- as.integer(b)
  if (b==1) {
    return(1)
  }
  all_rates <- beta_coal_factor_rates(b,alpha)
  return(all_rates/sum(all_rates))
}

# check if given function is density on [0,1]
is_density <- function(lambda_meas) {
  integral <- integrate(lambda_meas, lower =0, upper = 1)$value
  diff <- abs(integral-1)
  return( diff < 1e-5)
}

# check if given discrete distribution is pmf on [0,1]
is_pmf <- function(lambda_pmf) {
  if (any( names(lambda_pmf) != c('values', 'probs'))) {
    return(FALSE)
  }
  
  if (any(lambda_pmf$values <0) | any(lambda_pmf$values >1 )) {
    return(FALSE)
  }
  
  if (abs(sum(lambda_pmf$probs) - 1) > 1e-5) {
    return(FALSE)
  }
  return(TRUE)
}

# Calculate (b choose k) lambda_b,k for some Lambda measure
lambda_coal_factor <- function(b, k, lambda_meas) {
  b <- as.integer(b)
  k <- as.integer(k)
  if (b == 1) {
    return(1)
  }
  
  
  if (class(lambda_meas) == 'list') {
    if (!is_pmf(lambda_meas)) {
      stop('Invalid discrete Lambda-measure')
    }
    coal_factor_values <- lambda_meas$values^(k-2) * (1-lambda_meas$values)^(b-k)
    coal_factor <- sum(lambda_meas$probs * coal_factor_values)
    return(coal_factor * choose(b,k))
  }
  
  if (class(lambda_meas) == 'function') {
    if (!is_density(lambda_meas)) {
      stop('Invalid continuous Lambda-measure')
    }
    coal_factor_integrand <- function(x) { x^(k-2) * (1-x)^(b-k) * lambda_meas(x)}
    coal_factor <- integrate(coal_factor_integrand, lower =0, upper = 1)$value
    return(coal_factor * choose(b,k))
  }
  
  stop('Invalid Lambda-measure')
}

# Calculate (b choose k) lambda_b,k, k=2,...,b for some Lambda measure
lambda_coal_factor_rates <- function(b, lambda_meas) {
  b <- as.integer(b)
  if (b==1) {
    return(1)
  }
  return(sapply(2:b, lambda_coal_factor, b=b, lambda_meas=lambda_meas))
}

# Calculate lambda_b for some Lambda measure
lambda_coal_factor_total <- function(b, lambda_meas) {
  b <- as.integer(b)
  if (b==1) {
    return(1)
  }
  return(sum(sapply(2:b, lambda_coal_factor, b=b, lambda_meas = lambda_meas)))
}

# Calculate (b choose k) lambda_b,k/lambda_b for some Lambda measure
lambda_coal_factor_prop <- function(b, lambda_meas) {
  b <- as.integer(b)
  if (b==1) {
    return(1)
  }
  all_rates <- lambda_coal_factor_rates(b,lambda_meas)
  return(all_rates/sum(all_rates))
}

# Block-size likelihood under Beta-coalescent 
# Wrapper function 
beta_coal_block_size_log_lik <- function(alpha, block_sizes, lineages) {
  return(sum(mapply(function(n,k) {log(beta_coal_factor(n,k, alpha)) - log(beta_coal_factor_total(n, alpha)) }, 
                    lineages, block_sizes)))
}


#' Estimate of alpha parameter of a multifurcating tree under
#' the Beta-coalescent model using only topology information
#'
#' @param tree A \code{phylo} object 
#'
#' @return The block-size MLE for estimating alpha 
#' @export
beta_coal_block_size_alpha_mle <- function(tree) {
  summarize_tree <- summarize_multif_phylo(tree)
  block_sizes <- summarize_tree$block_sizes
  
  if (all(block_sizes == 2)) {
    return(2)
  }
  
  N <- sum(summarize_tree$n_sampled)
  coal_times <- summarize_tree$coal_times
  K <- length(block_sizes)
  lineages <- c(N, N-cumsum(block_sizes-1))
  lineages <- lineages[-length(lineages)]
  
  optimum <- optimize(beta_coal_block_size_log_lik, interval=c(0,1.99), 
                      block_sizes = block_sizes, lineages=lineages, maximum=TRUE)
  return(optimum$maximum)
}

# Many functions to compute log-likelihood of tree under Beta-coalescent model
beta_coal_lik_init <- function(alpha_multif, block_sizes, samp_times, n_sampled, coal_times, grid) {
  ns = length(samp_times)
  nc = length(coal_times)
  ng = length(grid)-1
  
  if (length(samp_times) != length(n_sampled)) {
    stop("samp_times vector of differing length than n_sampled vector.")
  }
  
  if (max(samp_times, coal_times) > max(grid)) {
    stop("Grid does not envelop all sampling and/or coalescent times.")
  } 
  
  t = sort(unique(c(samp_times, coal_times, grid)))
  l = rep(0, length(t))
  
  for (i in 1:ns) {
    l[t >= samp_times[i]] = l[t >= samp_times[i]] + n_sampled[i]
  }
  
  for (i in 1:nc) {
    l[t >= coal_times[i]] = l[t >= coal_times[i]] - block_sizes[i] + 1
  }
  
  l_diff <- -diff(l)
  l_block <- c(ifelse(l_diff<=0, 2, l_diff+1),2)
  
  if (sum((l < 1) & (t >= min(samp_times))) > 0) {
    stop("Number of active lineages falls below 1 after the first sampling point.")
  }
  
  mask = l > 0
  t = t[mask]
  l = utils::head(l[mask], -1)
  l_block = utils::head(l_block[mask], -1)
  
  gridrep = rep(0, ng)
  
  for (i in 1:ng) {
    gridrep[i] = sum(t > grid[i] & t <= grid[i+1])
  }
  
  C_block = mapply(beta_coal_factor, b=l, k = l_block, alpha= alpha_multif)
  C_total = sapply(l, beta_coal_factor_total, alpha= alpha_multif)
  D = diff(t)
  
  y = rep(0, length(D))
  y[t[-1] %in% coal_times] = 1
  
  bins = cut(x = samp_times, breaks = t, include.lowest = TRUE)
  tab <- stats::aggregate(n_sampled ~ bins, FUN = sum, labels = FALSE)
  count <- rep(0, length(D))
  count[as.numeric(tab$bins)] <- tab$n_sampled
  count[utils::head(t, -1) >= max(samp_times)] <- NA
  
  rep_idx = cumsum(gridrep)
  rep_idx = cbind(rep_idx-gridrep+1,rep_idx)
  
  return(list(t=t, l=l, l_block=l_block, C_block=C_block, C_total=C_total, D=D, y=y, count=count, gridrep=gridrep, ns=sum(n_sampled), nc=nc, ng=ng, rep_idx=rep_idx, 
              args=list(alpha_multif= alpha_multif, block_sizes = block_sizes, samp_times=samp_times, n_sampled=n_sampled, coal_times=coal_times, grid=grid)))
}

beta_coal_loglik <- function(init, f, grad=FALSE) {
  if (init$ng != length(f))
    stop(paste("Incorrect length for f; should be", init$ng))
  
  f = rep(f, init$gridrep)
  
  llnocoal = init$D * init$C_total * exp(-f)
  
  if (!grad) {
    lls = init$y * log(init$C_block) -init$y * f - llnocoal
    ll = sum(lls[!is.nan(lls)])
    return(ll)
  } else {
    dll = apply(init$rep_idx,1,function(idx)sum(-init$y[idx[1]:idx[2]]+llnocoal[idx[1]:idx[2]]))
    return(dll)
  }
}

beta_coal_loglik_wrapper <- function(alpha, block_sizes, samp_times, n_sampled, coal_times, f, grid) {
  tree_init <- beta_coal_lik_init(alpha, block_sizes, samp_times, n_sampled, coal_times, grid)
  return(-beta_coal_loglik(tree_init, f))
}

#' Simulate a multifurcating tree under the Beta-coalescent 
#'
#' @param N An integer giving the number of tips in the tree
#' @param alpha A number between 0 and 2 giving the parameter for the Beta coalescent  
#' @param tip.label a character vector giving the tip labels; if not specified, the tips \code{t1, t2}, ..., are
#' given 
#'
#' @return An object of class \code{phylo}
#' @export
rbeta_coal <- function (N, alpha, tip.label = NULL) {
  N <- as.integer(N)
  
  if (alpha <= 0 | alpha > 2) {
    stop('Invalid alpha parameter for Beta-coalescent')
  }
  
  if (alpha == 2) {
    return(rcoal(N))
  } 
  edges <- matrix(nrow=0, ncol=2) 
  if (N == 2) {
    edge <- matrix(c(3,1,3,2), ncol=2, byrow=TRUE)
    Nnode <- 1 
    if (is.null(tip.label)) {
      tip.label <- c('t1', 't2')
    }
    tip.label <- sample(tip.label)
    edge.length <- rep(rexp(1), 2)
    
    tree <- list(edge = edge, tip.label = tip.label, Nnode = Nnode, edge.length=edge.length)
    class(tree) <- "phylo"
    return(tree)
  }
  
  curr_lineages <- 1:N 
  curr_n_lineages <- length(curr_lineages)
  curr_internal_node <- N 
  branch_lengths <- c()
  all_edge_lengths <- c()
  
  while (curr_n_lineages >1) {
    
    curr_internal_node <- curr_internal_node + 1 
    curr_total_coal_rate <- beta_coal_factor_total(curr_n_lineages, alpha)
    curr_br_length <- rexp(1, rate = curr_total_coal_rate)
    branch_lengths <- c(curr_br_length, branch_lengths)
    
    curr_block_size <- sample(curr_n_lineages - 1, 1, prob = beta_coal_factor_prop(curr_n_lineages, alpha))+ 1
    curr_merge <- sample(curr_lineages, curr_block_size)
    
    new_edges <- matrix(c(rep(curr_internal_node, curr_block_size), curr_merge), ncol=2)
    edges <- rbind(edges, new_edges)
    
    new_edges_node_depth <- new_edges 
    new_edges_node_depth[new_edges[,2] <= N,2] <- N
    new_node_depth <- new_edges_node_depth[,1] - new_edges_node_depth[,2]
    
    new_branch_lengths <- cumsum(branch_lengths)[new_node_depth]
    all_edge_lengths <- c(all_edge_lengths, new_branch_lengths)
    
    curr_lineages <- setdiff(curr_lineages, curr_merge)
    curr_lineages <- c(curr_lineages, curr_internal_node)
    curr_n_lineages <- length(curr_lineages)
  }
  
  if (is.null(tip.label)) {
    tip.label <- paste("t", 1:N, sep = "")
  }
  tip.label <- sample(tip.label)
  
  edges <- edges[nrow(edges):1, ]
  edges[edges[, 2] <= N, 2] <- 1:N
  all_edge_lengths <- rev(all_edge_lengths)
  edges[edges[, 1] > N, 1] <- curr_internal_node - edges[edges[, 1] > N, 1] + N+1
  edges[edges[, 2] > N, 2] <- curr_internal_node - edges[edges[, 2] > N, 2] + N+1
  
  tree <- list(edge = edges, tip.label = tip.label, Nnode = curr_internal_node - N, edge.length=all_edge_lengths)
  
  class(tree) <- "phylo"
  
  tree <- read.tree(text=write.tree(tree))
  
  return(tree)
  
}

#' Simulate from inhomogeneous, heterochronous Beta-coalescent
#' 
#' @param samp_times numeric vector of sampling times.
#' @param n_sampled numeric vector of samples taken per sampling time.
#' @param alpha number between 0 and 2 giving the parameter for the Beta distribution 
#' @param traj function that returns effective population size at time t.
#' @param method which sampling method to use. "tt" invoke time-transformation
#'   method, "thin" invokes thinning method.
#' @param val_upper numeric used by time-transformation method to set a starting
#'   point for its dynamic numerical integration upper bound.
#' @param lower_bound numeric lower limit of \code{traj} function on its
#'   support.  Used only by thinning method.
#' @param ... additional arguments to be passed to \code{traj} function.
#'   
#' @return A list containing vectors of coalescent times \code{coal_times}, 
#'   intercoalescent times \code{intercoal_times}, number of active lineages
#'   \code{lineages}, block-size at each coalescent event \code{block_sizes},
#'   as well as passing along \code{samp_times} and \code{n_sampled}.
#' @export
#' 
#' @examples
#' coalsim(0:2, 3:1, alpha =1, unif_traj)
beta_coalsim <- function (samp_times, n_sampled, alpha, traj, method = "tt", val_upper = 10, lower_bound = 1, ...) {
  
  if (length(n_sampled) != length(samp_times)) {
    stop('Unequal number of sampling times and length of number sampled')
  }
  
  if (method == "tt") {
    result = beta_coalsimtt(samp_times, n_sampled, alpha, traj, val_upper, 
                            ...)
  } else if (method == "thin") {
    result = beta_coalsim_thin(samp_times, n_sampled, alpha, traj, lower_bound, 
                               ...)
  } else {
    stop("Argument method not recognized.")
  }
  return(result)
}

# internal functions used for simulating multifurcating tree 
beta_hazard_uniroot_stepfun <- function (traj_inv_stepfun, alpha, lineages, start, target) {
  knots = knots(traj_inv_stepfun)
  lin_factor = beta_coal_factor_total(lineages, alpha)
  t = start
  while (target > 0 && sum(knots > t) > 0) {
    next_knot = min(knots[knots > t])
    if ((next_knot - t) * lin_factor * traj_inv_stepfun(mean(c(t, next_knot))) > target) {
      result = t + target/(lin_factor * traj_inv_stepfun(mean(c(t, next_knot))))
      target = 0
    }
    else {
      target = target - (next_knot - t) * lin_factor * traj_inv_stepfun(mean(c(t, next_knot)))
      t = next_knot
    }
  }
  if (sum(knots > t) < 1) 
    result = t + target/(lin_factor * traj_inv_stepfun(t + 1))
  return(result - start)
}

beta_coalsimtt <- function (samp_times, n_sampled, alpha, traj, val_upper = 10, ...) {
  if (stats::is.stepfun(traj)) {
    knots = knots(traj)
    midpts = c(min(knots) - 1, knots[-1] - diff(knots)/2, max(knots) + 1)
    traj_inv <- stats::stepfun(x = knots, y = 1/traj(midpts))
    hazard <- function(t, lins, alpha, start, target) { beta_coal_factor_total(lins, alpha) * 
        integrate_step_fun(traj_inv, start,  start + t) - target }
    is_stepfun = TRUE
  } else {
    traj_inv <- function(t) 1/traj(t, ...)
    hazard <- function(t, lins, alpha, start, target){ beta_coal_factor_total(lins, alpha) *
        stats::integrate(traj_inv, start, start + t)$value - target} 
    is_stepfun = FALSE
  }
  
  coal_times = NULL
  lineages = NULL
  block_sizes = NULL
  
  curr = 1
  active_lineages = n_sampled[curr]
  time = samp_times[curr]
  
  while (time <= max(samp_times) || active_lineages > 1) {
    if (active_lineages == 1) {
      curr <- curr + 1
      active_lineages <- active_lineages + n_sampled[curr]
      time <- samp_times[curr]
    }
    target <- stats::rexp(1)
    if (is_stepfun) {
      y <- beta_hazard_uniroot_stepfun(traj_inv_stepfun = traj_inv, alpha = alpha,
                                       lineages = active_lineages, start = time, target = target)
    }
    else {
      y <- stats::uniroot(hazard, lins = active_lineages, alpha = alpha, start = time, target = target, 
                          lower = 0, upper = val_upper, extendInt = "upX")$root
    }
    
    while (curr < length(samp_times) && time + y >= samp_times[curr + 1]) {
      target <- -hazard(t = samp_times[curr + 1] - time,lins = active_lineages, alpha=alpha,
                        start = time, target = target)
      curr <- curr + 1
      active_lineages <- active_lineages + n_sampled[curr]
      time <- samp_times[curr]
      if (is_stepfun) {
        y <- beta_hazard_uniroot_stepfun(traj_inv_stepfun = traj_inv, alpha = alpha,
                                         lineages = active_lineages, start = time, target = target)
      }
      else {
        y <- stats::uniroot(hazard, lins = active_lineages, alpha=alpha, start = time, target = target, 
                            lower = 0, upper = val_upper, extendInt = "upX")$root
      }
    }
    time <- time + y
    coal_times = c(coal_times, time)
    lineages = c(lineages, active_lineages)
    
    curr_block_size = sample(active_lineages - 1, 1, prob = beta_coal_factor_prop(active_lineages, alpha))+ 1
    
    block_sizes = c(block_sizes, curr_block_size)
    active_lineages = active_lineages - curr_block_size + 1
  }
  
  return(list(coal_times = coal_times, lineages = lineages, block_sizes = block_sizes,
              intercoal_times = c(coal_times[1], diff(coal_times)), 
              samp_times = samp_times, n_sampled = n_sampled, alpha=alpha))
}

beta_coalsim_thin <- function (samp_times, n_sampled, alpha, traj, lower_bound, ...) {
  coal_times = NULL
  lineages = NULL
  block_sizes = NULL
  
  curr = 1
  active_lineages = n_sampled[curr]
  time = samp_times[curr]
  while (time <= max(samp_times) || active_lineages > 1) {
    
    if (active_lineages == 1) {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    
    time = time + stats::rexp(1, beta_coal_factor_total(active_lineages, alpha)/lower_bound)
    if (curr < length(samp_times) && time >= samp_times[curr + 1]) {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    } else if (stats::runif(1) <= lower_bound/traj(time, ...)) {
      coal_times = c(coal_times, time)
      lineages = c(lineages, active_lineages)
      
      curr_block_size = sample(active_lineages - 1, 1, prob = beta_coal_factor_prop(active_lineages, alpha))+ 1
      block_sizes = c(block_sizes, curr_block_size)
      active_lineages = active_lineages - curr_block_size + 1
    }
  }
  return(list(coal_times = coal_times, lineages = lineages, block_sizes = block_sizes, 
              intercoal_times = c(coal_times[1], diff(coal_times)), 
              samp_times = samp_times, n_sampled = n_sampled, alpha = alpha))
}

#' Returns a phylo object from the argumentes generated with beta_coalsim
#' 
#' @param args is a list containing vectors of coalescent times \code{coal_times}, sampling times \code{samp_times}, 
#'  block sizes \code{block_sizes}, and number sampled per sampling time \code{n_sampled}, etc. This list is the output of coalsim
#' @return A list with two elements \code{newikck} contains the tree in phylo format, \code{lables} a vector with tip labels 
#' @export
#' 
#' @examples
#' simulation1<-beta_coalsim(samp_times=c(0,1), n_sampled = c(25,25), alpha=1.5, traj=exp_traj)
#' tree1<- generate_newick_multif(simulation1)$newick
#' plot(tree1$newick)
generate_newick_multif <- function(args) {
  n <- sum(args$n_sampled) 
  if (length(args$lineages) == n-1) {
    return(phylodyn::generate_newick_rand(args))
  }
  
  labels <- paste(rep("t", n), seq(1, n, 1), rep("_", n), rep(args$samp_times[1], args$n_sampled[1]), sep = "")
  tb <- args$n_sampled[1]
  s <- 0
  temp_labels <- labels[1:tb]
  temp_times <- rep(args$samp_times[1], args$n_sampled[1])
  initial.row <- 2
  args2 <- gen_INLA_args_multif(args$block_sizes, args$samp_times, args$coal_times, args$n_sampled, args$alpha)
  args2$lineages <- c(args2$lineages ,1)
  for (j in 2:length(args2$event)) {
    if (args2$event[j] == 1) { 
      s <- args2$s[j] 
      n_merge <- args2$lineages[j-1] - args2$lineages[j] + 1
      ra <- sample(tb, n_merge)
      ra <- sort(ra)
      if (tb>1) {
        new_label <- paste("(", temp_labels[ra[1]], ":", s - temp_times[ra[1]], ",", sep='')
        ind_remove <- ra[-1]
        for (ind in ind_remove) {
          new_label <- paste(new_label, temp_labels[ind], ":", s - temp_times[ind], ',', sep='')
        }
        new_label <- substr(new_label, 1, nchar(new_label)-1)
        new_label <- paste(new_label, ')', sep='')
        ind_remove <- ra[-1]
        temp_labels[ra[1]] <- new_label
        temp_labels <- temp_labels[-ind_remove]
        temp_times[ra] <- s
        temp_times <- temp_times[-ind_remove]
      }
      tb <- tb - n_merge+1
    } else {
      s <- args2$s[j]
      if (args$n_sample[initial.row] == 1) {
        temp_labels <- c(temp_labels, labels[cumsum(args$n_sampled)[initial.row]])
        initial.row <- initial.row + 1
        tb <- tb + 1
        temp_times <- c(temp_times, s)
      }
      else {
        end <- cumsum(args$n_sampled)[initial.row]
        ini <- cumsum(args$n_sampled)[initial.row - 1] + 1
        for (k in ini:end) {
          temp_labels <- c(temp_labels, labels[k])
          tb <- tb + 1
          temp_times <- c(temp_times, s)
        }
        initial.row <- initial.row + 1
      }
    }
  }
  
  out.tree <- ape::read.tree(text = paste(temp_labels, ";", sep = ""))
  return(list(newick = out.tree, labels = labels))
}


#' Summarize a multifurcating phylogeny.
#' 
#' @param phy a \code{phylo} object containing a phylogeny.
#'   
#' @return A list containing vectors of sampling times \code{samp_times}, number 
#'   sampled per sampling time \code{n_sampled}, coalescent times \code{coal_times},
#'   and block sizes per coalescent event \code{block_sizes}. 
#' @export
summarize_multif_phylo<-function(tree){
  hgpstat <- heterochronous_gp_stat_multif(tree)
  return(list(block_sizes=hgpstat$block_sizes,samp_times = hgpstat$samp_times, n_sampled = hgpstat$n_sampled, 
              coal_times = hgpstat$coal_times))
}

# Internal function 
branching_sampling_times_multif <- function(tree) {
  if (class(tree) != "phylo") 
    stop("object \"tree\" is not of class \"phylo\"")
  edge.mat <- tree$edge
  
  K <- tree$Nnode
  N <- max(tree$edge) - tree$Nnode
  t.tot <- max(ape::node.depth.edgelength(tree))
  n.t <- t.tot - ape::node.depth.edgelength(tree)
  xx <- as.numeric(rep(NA, N+K))
  
  names(xx) <- as.character(c(-(1:K), 1:N))
  
  xx[1:K] <- sort(n.t[(N+1):length(n.t)], decreasing = TRUE)
  xx[(K+1):length(xx)] <- sort(n.t[1:(N)], decreasing = TRUE)
  return(xx)
}

# Internal function 
heterochronous_gp_stat_multif<-function(phy,tol = 10^(-6)){ 
  b.s.times = branching_sampling_times_multif(phy)
  int.ind = which(as.numeric(names(b.s.times)) < 0)
  tip.ind = which(as.numeric(names(b.s.times)) > 0)
  num.tips = length(tip.ind)
  num.coal.events = length(int.ind)
  sampl.suf.stat = rep(NA, num.coal.events)
  sorted.coal.times = sort(b.s.times[int.ind])
  names(sorted.coal.times) = NULL
  sampling.times = sort((b.s.times[tip.ind]))
  
  for (i in 2:length(sampling.times)) {
    if ((sampling.times[i] - sampling.times[i - 1]) < tol) {
      sampling.times[i] <- sampling.times[i - 1]
    }
  }
  unique.sampling.times <- unique(sampling.times)
  sampled.lineages = NULL
  for (sample.time in unique.sampling.times) {
    sampled.lineages = c(sampled.lineages, sum(sampling.times == sample.time))
  }
  
  edges <-phy$edge
  K <- phy$Nnode
  N <- max(phy$edge) - K
  
  values1<-max(ape::node.depth.edgelength(phy))-ape::node.depth.edgelength(phy)
  values1<-values1[(N+1):length(values1)]
  values<-sort(values1,decreasing=T)
  correct.label<-seq(N+1,N+K)
  reference<-matrix(c(values1,seq(N+1,N+K),rep(0,length(values1))),ncol=3)
  
  for (j in 1:length(values1)){
    reference[reference[,1]==values[j],3]<-correct.label[j]
  }
  
  newedges<-matrix(0,nrow=nrow(edges),ncol=4)
  newedges[,1:2]<-edges
  newedges[,3:4]<-edges
  
  for (j in (N+1):(N+K)){
    newedges[newedges[,1]==j,3]<-reference[reference[,2]==j,3]
    newedges[newedges[,2]==j,4]<-reference[reference[,2]==j,3]
  }
  
  edges<-newedges[,3:4]
  ordered.block.sizes<- sapply(1:K, function(i) {length(which(edges[,1] == N+i))})
  ordered.block.sizes <- rev(ordered.block.sizes)
  
  return(list(block_sizes = ordered.block.sizes, coal_times = sorted.coal.times, samp_times = unique.sampling.times, 
              n_sampled = sampled.lineages))
}

# Internal function 
gen_INLA_args_multif<-function(block_sizes, samp_times, coal_times, n_sampled, alpha_multif, lambda_meas=NULL) {
  if (length(intersect(coal_times, samp_times)) > 0) {
    warning("Coincident sampling event and coalescent event: results may be unpredictable.")
  }
  
  if (is.null(alpha_multif) & is.null(lambda_meas)) {
    stop('Invalid Lambda measure provided')
  } else if (!is.null(alpha_multif)) {
    if (alpha_multif > 2 | alpha_multif <= 0) {
      stop('Invalid alpha parameter for Beta-coalescent')
    }
  }
  
  l <- length(samp_times)
  m <- length(coal_times)
  sorting <- sort(c(samp_times, coal_times), index.return = TRUE)
  lineage_change <- c(n_sampled, -block_sizes+1)[sorting$ix]
  lineages <- utils::head(cumsum(lineage_change), -1)
  
  if (!is.null(alpha_multif)) {
    if (alpha_multif==2) {
      coal_factor <-choose(lineages,block_sizes)
    } else{
      coal_factor <- sapply(lineages, beta_coal_factor_rates, alpha=alpha_multif)
      coal_factor <- unlist(lapply(coal_factor, sum))
    }
  } else {
    coal_factor <- sapply(lineages, lambda_coal_factor_rates, lambda_meas=lambda_meas)
    coal_factor <- unlist(lapply(coal_factor, sum))
  }
  event <- c(rep(0, l), rep(1, m))[sorting$ix]
  return(list(coal_factor = coal_factor, s = sorting$x, event = event,lineages = lineages))
}

# Internal function 
coal_stats_multif<-function (grid, block_sizes, samp_times, coal_times, n_sampled, alpha_multif, lambda_meas=NULL, log_zero = -100) {
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2
  if (is.null(n_sampled)) {
    n_sampled <- rep(1, length(samp_times))
  }
  
  args <- gen_INLA_args_multif(block_sizes=block_sizes,samp_times = samp_times, 
                               coal_times = coal_times, n_sampled = n_sampled, 
                               alpha_multif = alpha_multif, lambda_meas=lambda_meas)
  
  coal_factor <- args$coal_factor
  s <- args$s
  event <- args$event
  grid_trimmed <- setdiff(x = grid, y = s)
  sorting <- sort(c(grid_trimmed, s), index.return = TRUE)
  sgrid <- sorting$x
  ordering <- sorting$ix
  time_index <- cut(x = sgrid[-1], breaks = grid, labels = FALSE)
  time <- field[time_index]
  event_out <- c(rep(0, length(grid_trimmed)), event)[ordering]
  Cfun <- stats::stepfun(x = s, y = c(0, coal_factor, 0), right = TRUE)
  Cvec <- Cfun(sgrid[-1])
  E <- diff(sgrid) * Cvec
  E_log = log(E)
  E_log[E == 0] = log_zero
  return(data.frame(time = time, event = event_out[-1], E = E, E_log = E_log))
}

# Internal function 
infer_coal_samp_multif <-function(block_sizes, samp_times, coal_times,  n_sampled, alpha_multif, lambda_meas=NULL,
                                  lengthout = 100, fns = NULL, prec_alpha = 0.01, prec_beta = 0.01, beta1_prec = 0.001, use_samp = FALSE,
                                  log_fns = TRUE, simplify = TRUE, events_only = FALSE, derivative = FALSE, grid = NULL)  {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("INLA needed for this function to work. Use install.packages(\"INLA\", repos=c(getOption(\"repos\"), INLA=\"https://inla.r-inla-download.org/R/stable\"), dep=TRUE).", 
         call. = FALSE)
  }
  if (min(coal_times) < min(samp_times)) {
    stop("First coalescent time occurs before first sampling time")
  }
  if (max(samp_times) > max(coal_times)) {
    stop("Last sampling time occurs after last coalescent time")
  }
  
  if (is.null(grid)) {
    grid <- seq(min(samp_times), max(coal_times), length.out = lengthout +1)
  }
  
  if (is.null(n_sampled)) {
    n_sampled <- rep(1, length(samp_times))
  }
  
  coal_data <- coal_stats_multif(grid = grid, block_sizes=block_sizes, samp_times = samp_times, coal_times = coal_times, 
                                 n_sampled = n_sampled, alpha_multif = alpha_multif,lambda_meas=lambda_meas)
  
  if (simplify) {
    coal_data <- with(coal_data, phylodyn:::condense_stats(time = time, event = event, E = E))
  }
  
  
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  
  # if (!use_samp & !multia) {
  #     data <- with(coal_data, data.frame(y = event, time = time, E_log = E_log))
  #     formula <- y ~ -1 + f(time, model = "rw1", hyper = hyper, 
  #                           constr = FALSE, scale.model = TRUE)
  #     family <- "poisson"
  # } else if (!use_samp & multia) {
  #   sorting <- sort(c(samp_times, coal_times), index.return = TRUE)
  #   lineage_change <- c(n_sampled, -block_sizes+1)[sorting$ix]
  #   lineages <- utils::head(cumsum(lineage_change), -1)
  #   s <-sorting$x 
  # 
  #   b_lineages <- lineages[ sapply(grid, function(x) { which(s == min(s[s >= x])) }) - 1]
  #   # b_lineages <- stats::stepfun(x = s, y = c(0, lineages, 0), right = TRUE)
  # 
  #   # E_log_new <- log(b_lineages-1) +coal_data$E_log
  #   E_log_new <- log(b_lineages-1) + log(diff(grid))
  #   
  #   data <- with(coal_data, data.frame(y = event, time = time, E_log = E_log_new))
  #   data$size <- log(b_lineages/2)
  #   formula <- y ~ -1 + size+ f(time, model = "rw1", hyper = hyper, 
  #                         constr = FALSE, scale.model = TRUE)
  #   family <- "poisson"
  
  if (!use_samp) {
    data <- with(coal_data, data.frame(y = event, time = time, E_log = E_log))
    formula <- y ~ -1 + f(time, model = "rw1", hyper = hyper,
                          constr = FALSE, scale.model = TRUE)
    family <- "poisson"
  } else if (use_samp) {
    if (events_only) {
      samp_data <- phylodyn:::samp_stats(grid = grid, samp_times = samp_times)
    } else {
      samp_data <- phylodyn:::samp_stats(grid = grid, samp_times = samp_times, n_sampled = n_sampled)
    }
    
    data <- phylodyn:::joint_stats(coal_data = coal_data, samp_data = samp_data)
    
    if (is.null(fns)) {
      formula <- Y ~ -1 + beta0 + f(time, model = "rw1",  hyper = hyper, constr = FALSE) + 
        f(time2, w, copy = "time", fixed = FALSE, param = c(0, beta1_prec))
    } else {
      vals <- NULL
      bins <- sum(data$beta0 == 0)
      for (fni in fns) {
        if (log_fns) 
          vals <- cbind(vals, c(rep(0, bins), log(fni(samp_data$time))))
        else vals <- cbind(vals, c(rep(0, bins), fni(samp_data$time)))
      }
      data$fn <- vals
      formula <- Y ~ -1 + beta0 + fn + f(time, model = "rw1", hyper = hyper, constr = FALSE) + 
        f(time2, w, copy = "time", fixed = FALSE, param = c(0, beta1_prec))
    }
    
    family <- c("poisson", "poisson")
  } else {
    stop("Invalid use_samp value, should be boolean.")
  }
  
  if (derivative) {
    Imat <- diag(lengthout)
    A <- utils::head(Imat, -1) - utils::tail(Imat, -1)
    field <- grid[-1] - diff(grid)/2
    A <- diag(1/diff(field)) %*% A
    A[A == 0] <- NA
    lc_many <- INLA::inla.make.lincombs(time = A)
  } else {
    lc_many <- NULL
  }
  
  mod <- INLA::inla(formula, family = family, data = data,lincomb = lc_many, 
                    offset = data$E_log, control.predictor = list(compute = TRUE), 
                    control.compute = list(config = TRUE), control.inla(strategy = "adaptive", int.strategy = "eb"),num.threads=1)
  return(list(result = mod, data = data, grid = grid, x = coal_data$time))
}

#' Bayesian nonparametric phylodynamic reconstruction of multifurcating tree 
#' under the Lambda-coalescent 
#' 
#' @param data \code{phylo} object or list containing vectors of coalescent 
#'   times \code{coal_times}, sampling times \code{samp_times}, block sizes
#'   \code{block_sizes}, and number sampled per sampling time \code{n_sampled}.
#' @param alpha_multif number between 0 and 2 giving the parameter for the Beta coalescent
#' @param lambda_meas A probability density function on the unit interval, or a named list 
#'   of \code{values} and \code{probs} representing a probaiblity mass function on the 
#'   unit interval. Not used if \code{alpha_multif} is specified.  
#' @param lengthout numeric specifying number of grid points.
#' @param pref logical. Should the preferential sampling model be used? Default is \code{FALSE}. 
#' @param prec_alpha,prec_beta numerics specifying gamma prior for precision \eqn{\kappa}.
#' @param beta1_prec numeric specifying precision for normal prior on \eqn{\beta_1}.
#' @param fns list containing functions of covariates.
#' @param log_fns logical whether or not to to apply a log-transformation to
#'   the output of the functions in \code{fns}.
#' @param simplify logical whether to fully bucket all Poisson points.
#' @param derivative logical whether to calculate estimates of the 
#'   log-derivative.
#' @param forward logical whether to use the finite difference approximations of
#'   the log-derivative as a forward or backward derivative.
#' @param grid it gives the opportunity to evaluate BNPR on a prespecified grid
#'   
#' @return Phylodynamic reconstruction of effective population size at grid 
#'   points. \code{result} contains the INLA output, \code{data} contains the 
#'   information passed to INLA, \code{grid} contains the grid end points, 
#'   \code{x} contains the grid point centers, \code{effpop} contains a vector 
#'   of the posterior median effective population size estimates, 
#'   \code{effpop025} and \code{effpop975} contain the 2.5th and 97.5th 
#'   posterior percentiles, \code{summary} contains a data.frame of the 
#'   estimates, and \code{derivative} (if \code{derivative = TRUE}) contains a
#'   data.frame summarizing the log-derivative.
#' @export
BNPR_Lambda <- function(data, alpha_multif, lambda_meas=NULL, lengthout = 100, pref=FALSE, 
                        prec_alpha=0.01, prec_beta=0.01, beta1_prec = 0.001, fns = NULL, log_fns = TRUE,
                        simplify = TRUE, derivative = FALSE, forward = TRUE,grid=NULL) {
  if (is.null(alpha_multif) & is.null(lambda_meas)) {
    stop('Invalid Lambda-measure provided')
  } else if (!is.null(alpha_multif)) {
    if (alpha_multif > 2 | alpha_multif <= 0) {
      stop('Invalid alpha parameter for Beta-coalescent')
    }
  }
  
  if (class(data) == "phylo") {
    phy <- summarize_multif_phylo(data)
  } else if (all(c("block_sizes","coal_times", "samp_times", "n_sampled") %in% names(data))) {
    phy <- with(data, list(block_sizes=block_sizes, samp_times = samp_times, 
                           coal_times = coal_times, n_sampled = n_sampled))
  } else {
    stop('Invalid input data type')
  }
  
  if (!is.null(alpha_multif)) {
    if (alpha_multif == 2 & !all(phy$block_sizes == 2)) {
      stop('Attempting to infer from binary tree model given multifurcating tree')
    }
  }
  
  result <- infer_coal_samp_multif(block_sizes=phy$block_sizes,samp_times = phy$samp_times, coal_times = phy$coal_times,
                                   alpha_multif, lambda_meas, 
                                   n_sampled = phy$n_sampled, fns = fns, lengthout = lengthout, prec_alpha =prec_alpha,
                                   prec_beta = prec_beta, beta1_prec = beta1_prec, use_samp = pref, log_fns = log_fns,
                                   derivative = derivative, grid = grid, simplify = simplify)
  
  result$samp_times <- phy$samp_times
  result$n_sampled  <- phy$n_sampled
  result$coal_times <- phy$coal_times
  result$block_sizes <- phy$block_sizes
  
  result$effpop     <- exp(-result$result$summary.random$time$`0.5quant`)
  result$effpopmean <- exp(-result$result$summary.random$time$mean)
  result$effpop975  <- exp(-result$result$summary.random$time$`0.025quant`)
  result$effpop025  <- exp(-result$result$summary.random$time$`0.975quant`)
  
  result$summary <- with(result$result$summary.random$time,
                         data.frame(time = ID, mean = exp(-mean),
                                    sd = sd * exp(-mean),
                                    quant0.025 = exp(-`0.975quant`),
                                    quant0.5 = exp(-`0.5quant`),
                                    quant0.975 = exp(-`0.025quant`)))
  
  if (derivative) {
    if (forward) {
      ind <- c(1:(lengthout-1), (lengthout-1))
    } else {
      ind <- c(1, 1:(lengthout-1))
    }
    
    result$derivative <- with(result$result$summary.lincomb,
                              data.frame(time = result$x, mean = -mean[ind], sd = sd[ind],
                                         quant0.025 = -`0.975quant`[ind],
                                         quant0.5   = -`0.5quant`[ind],
                                         quant0.975 = -`0.025quant`[ind]))
  }
  
  if (pref){
    result$beta0     <- result$result$summary.fixed["beta0","0.5quant"]
    result$beta0summ <- result$result$summary.fixed["beta0",]
    rownames(result$beta0summ) <- "Beta0"
    result$beta1     <- result$result$summary.hyperpar[2,"0.5quant"]
    result$beta1summ <- result$result$summary.hyperpar[2,]
    rownames(result$beta1summ) <- "Beta1"
  }
  
  return(result)
}

#' Joint inference of alpha that characterizes the Beta-coalescent and 
#' the effective population size trajectory via hybrid method 
#'
#' @param tree A \code{phylo} object. 
#' @param tol tolerance for stopping the iteration.
#' @param lengthout numeric specifying number of grid points.
#' @param maxiter maximum number of iterations before termination if not converged. 
#'
#' @return \code{alpha_opt} optimal \eqn{\alpha} value found by hybrid method, 
#'   \code{alpha_block_size_MLE} block-size MLE of \eqn{\alpha}, \code{niter}
#'   the number of iterations it took to converge, and \code{BNPR} the output of 
#'   BNPR using the \code{alpha_opt} value. 
#' @export
posterior_beta_coal_hybrid <- function(tree, tol=1e-3, lengthout=100, maxiter=50) {
  phy <- summarize_multif_phylo(tree)
  alpha_block_size_MLE <- beta_coal_block_size_alpha_mle(tree)
  alpha_MLE_curr <- alpha_block_size_MLE
  alpha_MLE_diff <- 10
  alpha_niters <- 0
  while(alpha_MLE_diff > tol & alpha_niters <= maxiter) {
    result_INLA <- infer_coal_samp_multif(block_sizes=phy$block_sizes,samp_times = phy$samp_times, 
                                          coal_times = phy$coal_times, n_sampled = phy$n_sampled, alpha_multif=alpha_MLE_curr, lengthout = lengthout)
    opt_alpha_INLA <- optim(alpha_MLE_curr, fn=beta_coal_loglik_wrapper, lower= 1e-6, upper = 2-1e-6, method='L-BFGS-B',
                            block_sizes=phy$block_sizes, samp_times=phy$samp_times, n_sampled = phy$n_sampled, coal_times=phy$coal_times,
                            f = -result_INLA$result$summary.random$time$`0.5quant`, grid=result_INLA$grid)
    alpha_MLE_Ne_INLA <- opt_alpha_INLA$par
    alpha_MLE_diff <- abs(alpha_MLE_curr - alpha_MLE_Ne_INLA)
    alpha_MLE_curr <- alpha_MLE_Ne_INLA
    alpha_niters <- alpha_niters + 1
  }
  
  BNPR_final <- BNPR_Lambda(tree, alpha_multif = alpha_MLE_curr, lengthout=lengthout)
  print(paste('# iterations', alpha_niters))
  result <- list(alpha_opt = alpha_MLE_curr, alpha_block_size_MLE = alpha_block_size_MLE,  niter = alpha_niters, BNPR = BNPR_final)
  return(result)
}

# Internal functions for sHMC
beta_U_split <- function (theta, init, invC, alpha, beta, grad = FALSE) {
  D = length(theta)
  f = theta[-D]
  tau = theta[D]
  invCf = invC %*% f
  if (!grad) {
    loglik = beta_coal_loglik(init, f)
    logpri = ((D - 1)/2 + alpha) * tau - (t(f) %*% invCf/2 + beta) * exp(tau)
    return(list(loglik = -loglik, logpri = -logpri, logpos = -(loglik + logpri)))
  } else {
    dU_res = -c(beta_coal_loglik(init, f, grad), ((D - 1)/2 + alpha) - beta * exp(tau))
    return(dU_res)
  }
}

beta_coal_sampling <- function (data, para, alg, setting, init, verbose = TRUE, printevery = 100) {
  lik_init = data$lik_init
  Ngrid = lik_init$ng + 1
  alpha = para$alpha
  beta = para$beta
  invC = para$invC
  rtEV = para$rtEV
  EVC = para$EVC
  cholC = para$cholC
  stepsz = setting$stepsz
  Nleap = setting$Nleap
  
  rand_leap = setting$rand_leap
  
  NSAMP = setting$NSAMP
  NBURNIN = setting$NBURNIN
  NSUBSAMP = setting$NSUBSAMP
  recorded_iters = seq.int(from = NBURNIN + 1, to = NSAMP, 
                           by = NSUBSAMP)
  samp = matrix(NA, length(recorded_iters), Ngrid)
  acpi = 0
  acpt = rep(NA, length(recorded_iters))
  logpri = rep(NA, length(recorded_iters))
  loglik = rep(NA, length(recorded_iters))
  logpos = rep(NA, length(recorded_iters))
  theta = init$theta
  u = init$u
  du = init$du
  
  Ufun = function(theta, grad = FALSE) beta_U_split(theta = theta,init = lik_init, invC = invC, 
                                                    alpha = alpha, beta = beta, grad = grad)
  colnames(samp) = c(paste("f", 1:(Ngrid - 1), sep = ""), "tau")
  
  start_time = Sys.time()
  # cat("Running ", alg, " sampling...\n")
  for (iter in 1:NSAMP) {
    res = phylodyn:::splitHMC(theta, u, du, Ufun, rtEV, EVC, stepsz,Nleap, rand_leap)
    acpi <- acpi + res$Ind
    theta[1:Ngrid] <- res$q
    u <- res$u
    du <- res$du
    output_index = match(x = iter, table = recorded_iters)
    if (!is.na(output_index)) {
      samp[output_index, ] <- theta
      acpt[output_index] <- res$Ind
      logpri[output_index] <- res$pos_summ$logpri
      loglik[output_index] <- res$pos_summ$loglik
      logpos[output_index] <- res$pos_summ$logpos
    }
    if (verbose && iter%%printevery == 0) {
      cat(iter, " iterations have been finished!\n")
      cat("Online acceptance rate is ", acpi/printevery, 
          "\n")
      acpi = 0
    }
  }
  #stop_time <- Sys.time()
  #time <- stop_time - start_time
  #cat("\nTime consumed : ", time, units(time))
  #cat("\nFinal Acceptance Rate: ", sum(acpt)/(NSAMP - NBURNIN), "\n")
  pos_summ = data.frame(acpt = acpt, logpri = logpri, loglik = loglik, 
                        logpos = logpos)
  result = list(samp = samp, alg = alg, time = time, pos_summ = pos_summ)
  return(result)
}

beta_coal_mcmc_sampling = function(dataset, alpha_multif, alg='splitHMC', nsamp, nburnin=0, nsubsamp=1, ngrid=100,
                                   nugget="1,1", prec_alpha = 1e-2, prec_beta = 1e-2,
                                   TrjL=NULL, Nleap=NULL, szkappa=NULL, rand_leap=NULL,
                                   f_init = rep(1, ngrid-1), kappa = 1,
                                   covariates=NULL, betas=rep(0, 2+length(covariates)),
                                   samp_alg = "none", kappa_alg = "gibbs",
                                   beta_vars = rep(100, length(betas)), printevery=100) {
  if (alg != 'splitHMC') {
    stop('Algorithm not supported')
  }
  
  if (class(dataset) == "phylo") {
    phy <- summarize_multif_phylo(dataset)
  } else if (all(c("coal_times", "samp_times", "n_sampled", 'block_sizes') %in% names(dataset))) {
    phy <- with(dataset, list(samp_times = samp_times, coal_times = coal_times,
                              n_sampled = n_sampled, block_sizes= block_sizes))
  }
  
  samp_times = phy$samp_times
  n_sampled  = phy$n_sampled
  coal_times = phy$coal_times
  block_sizes = phy$block_sizes
  
  # Jump tuning parameters--should probably have an option to change in the arguments
  if (is.null(TrjL)){
    TrjL = 3 
  }
  if (is.null(Nleap)) {
    Nleap =15
  }
  
  rand_leap = TRUE
  
  stepsz = TrjL/Nleap
  
  grid_bds = range(c(coal_times,samp_times))
  #Ngrid = 100
  
  grid = seq(grid_bds[1],grid_bds[2],length.out=ngrid)
  intl = grid[2]-grid[1]
  midpts = grid[-1]-intl/2
  
  covar_vals = NULL
  if (!is.null(covariates)) {
    for (fcn in covariates) {
      covar_vals = cbind(covar_vals, log(fcn(midpts)), deparse.level = 0)
    }
  }
  
  # initialize likelihood calculation
  lik_init = beta_coal_lik_init(alpha_multif=alpha_multif, block_sizes = block_sizes, samp_times=samp_times, n_sampled=n_sampled, coal_times=coal_times, grid=grid)
  
  # calculate intrinsic precision matrix
  invC <- phylodyn:::Q_matrix(midpts,0,1)
  
  # fudge to be able to compute the cholC
  if (nugget == "1,1")
    invC[1,1] <- invC[1,1]+.0001 # nugget at (1,1)
  else if (nugget == "diag")
    #Julia: I had a warning when using this--needs to be corrected
    diag(invC)<-diag(invC)+.0001 # nugget for the whole diagonal
  else if (nugget == "none")
    warning("No nugget may result in a non-full-rank matrix.")
  else
    stop(paste("Unrecognized argument nugget = '", nugget, "', please use '1,1', 'diag', or 'none'.", sep = ""))
  
  eig  = spam::eigen.spam(invC, TRUE)
  # eig = eigen(invC, TRUE)
  rtEV = sqrt(eig$values)
  EVC  = eig$vectors
  
  #C = spam::solve.spam(invC)
  C = solve(invC) 
  cholC = chol(C)
  
  # initializations
  theta = c(f_init, kappa)
  
  u  = beta_U_split(theta,lik_init,invC,prec_alpha,prec_beta)$logpos
  du = beta_U_split(theta,lik_init,invC,prec_alpha,prec_beta, TRUE)
  
  # MCMC sampling preparation
  dataset = list(lik_init = lik_init, covar_vals = covar_vals)
  para = list(alpha = prec_alpha, beta = prec_beta, invC = invC, rtEV = rtEV,
              EVC = EVC, cholC = cholC, betas = betas, beta_vars = beta_vars)
  setting = list(stepsz = stepsz, Nleap = Nleap,
                 NSAMP = nsamp, NBURNIN = nburnin, NSUBSAMP = nsubsamp,
                 szkappa = szkappa, rand_leap=rand_leap,
                 proposal_sds = rep(0.3, length(betas)))
  init = list(theta = theta, u = u, du = du, betas = betas)
  
  # Run MCMC sampler
  res_MCMC = beta_coal_sampling(data = dataset, para = para, alg = alg, setting = setting,
                                init = init, verbose=FALSE, printevery = printevery)
  
  
  res_MCMC$alg = alg
  res_MCMC$samp_alg = samp_alg
  res_MCMC$kappa_alg = kappa_alg
  res_MCMC$Ngrid = ngrid
  
  logfmat = res_MCMC$samp[,1:(ngrid-1)]
  params = matrix(res_MCMC$samp[,ngrid])
  
  res_MCMC$grid = grid
  res_MCMC$x = midpts
  res_MCMC$samp_times = samp_times
  res_MCMC$n_sampled = n_sampled
  res_MCMC$coal_times = coal_times
  res_MCMC$block_sizes = block_sizes
  
  return(res_MCMC)
}

#' Joint inference of alpha that characterizes the Beta-coalescent and 
#' the effective population size trajectory via MCMC 
#'
#' @param dataset \code{phylo} object or list containing vectors of coalescent 
#'   times \code{coal_times}, sampling times \code{samp_times}, block sizes
#'   \code{block_sizes}, and number sampled per sampling time \code{n_sampled}.
#' @param nsamp_total Integer representing the number of MCMC samples desired 
#' @param lengthout numeric specifying number of grid points.
#' @param burn-in numeric between 0 and 1 indicating proportion of samples to toss 
#'   at the start of the MCMC chain 
#'   
#' @return MCMC results of alpha and EPS sampled via sHMC: \code{\alpha} contains the 
#'   posterior sample of \eqn{\alpha},  \code{effpop} contains a vector of the posterior median 
#'   effective population size estimates, \code{grid} contains the grid end points, 
#'   \code{x} contains the grid point centers, \code{effpop025} and \code{effpop975} contain the 
#'   2.5th and 97.5th posterior percentiles, as well as the original list containing coalescent times,
#'   sampling times, block sizes, and number sampled per sampling time 
#' @export
posterior_beta_coal_MCMC_w_sHMC <- function(tree, nsamp_total=20000, lengthout=100, burn_in=0.1) {
  phy <- summarize_multif_phylo(tree)
  print(paste('# MCMC samples', nsamp_total))
  # initialize alpha and f 
  alpha_curr <- runif(1,0,2)
  # need to store f and kappa current 
  sHMC_curr <- beta_coal_mcmc_sampling(phy, alpha_curr, 'splitHMC', nsamp=1)
  f_curr <- sHMC_curr$samp[1, ][-lengthout]
  kappa_curr <- sHMC_curr$samp[1, ][lengthout]
  
  grid <- sHMC_curr$grid 
  
  # set the grid for discretization 
  alpha_grid <- 0.0025+ seq(0, 1.995, by=0.005)
  phy_lik_init <- lapply(alpha_grid, function(x) { beta_coal_lik_init(x, phy$block_sizes, phy$samp_times, phy$n_sampled, phy$coal_times, grid) })
  
  post_sample_alpha <- vector('double', length=nsamp_total)
  post_sample_f <- vector('list', length=nsamp_total)
  post_sample_kappa <- vector('double', length=nsamp_total)
  
  for (i in 1:nsamp_total) {
    discrete_lik_curr <- sapply(phy_lik_init, function(x) { beta_coal_loglik(x, f_curr)})
    discrete_prob_curr <- exp(discrete_lik_curr - max(discrete_lik_curr))
    discrete_prob_curr <- discrete_prob_curr/sum(discrete_prob_curr)
    if(! any(is.na(discrete_prob_curr))) {
      alpha_curr <- sample(alpha_grid, 1, prob=discrete_prob_curr)
      alpha_curr <- runif(1, min=alpha_curr-0.0025, max = alpha_curr + 0.0025)
    } else {
      print('subsample')
      discrete_lik_curr_subsample_ind <- order(discrete_lik_curr, decreasing=TRUE)[1:400]
      discrete_prob_curr <- exp(discrete_lik_curr[discrete_lik_curr_subsample_ind] - 
                                  max(discrete_lik_curr[discrete_lik_curr_subsample_ind]))
      print(discrete_lik_curr[discrete_lik_curr_subsample_ind] + 
              max(discrete_lik_curr[discrete_lik_curr_subsample_ind]))
      discrete_prob_curr <- discrete_prob_curr/sum(discrete_prob_curr)
      alpha_grid_curr <- alpha_grid[discrete_lik_curr_subsample_ind]
      alpha_curr <- sample(alpha_grid_curr, 1, prob=discrete_prob_curr)
      alpha_curr <- runif(1, min=alpha_curr-0.0025, max = alpha_curr + 0.0025)
    }
    if (i %% 1000 == 0) {
      print(paste('iter', i, ':', alpha_curr))
    }
    
    sHMC_curr <- beta_coal_mcmc_sampling(phy, alpha_curr, 'splitHMC', nsamp=1, f_init = f_curr, kappa = kappa_curr, prec_alpha = 1e-3)
    f_curr <- sHMC_curr$samp[1, ][-lengthout]
    kappa_curr <- sHMC_curr$samp[1, ][lengthout]
    
    post_sample_alpha[i] <- alpha_curr
    post_sample_f[[i]] <- as.vector(f_curr)
    post_sample_kappa[i] <- kappa_curr
  }
  
  post_sample_f <- do.call(cbind, post_sample_f)
  
  # add summary with burn-in  
  if (burn_in >0) {
    n_remove_burnin <- floor(burn_in * nsamp_total)
    post_sample_f <- post_sample_f[,-(1:n_remove_burnin)]
    post_sample_alpha <- post_sample_alpha[-(1:n_remove_burnin)]
    post_sample_kappa <- post_sample_kappa[-(1:n_remove_burnin)]
  }
  
  post_sample_results <- list(alpha= post_sample_alpha, eps = exp(-post_sample_f), 
                              f= -post_sample_f, grid=grid,x = sHMC_curr$x, 
                              block_sizes=phy$block_sizes,samp_times = phy$samp_times, 
                              coal_times = phy$coal_times, n_sampled = phy$n_sampled, kappa = post_sample_kappa)
  
  f_quantiles <- rowQuantiles(post_sample_results$f, probs=c(0.025, 0.5, 0.975))
  post_sample_results$effpop <- exp(-f_quantiles[,2])
  post_sample_results$effpopmean <- exp(rowMeans(-post_sample_results$f))
  post_sample_results$effpop975  <- exp(-f_quantiles[,1])
  post_sample_results$effpop025  <- exp(-f_quantiles[,3])
  
  return(post_sample_results)
}

# Performance measures that incorporate the structure of the output of MCMC 

#' Computes the coverage of the 95% credible interval of a reconstructed EPS 
#' trajectory given a reference trajectory 
#'
#' @param INLA_out output directly from \code{BNPR_Lambda} or \code{posterior_beta_coal_MCMC_w_sHMC},
#'  or the \code{BNPR} argument from the output of \code{posterior_beta_coal_hybrid}
#' @param traj the reference trajectory to compare to 
#'   
#' @return List of number of grid points covered \code{total} or percentage covered \code{avg}.
#' @export
envelope_new <- function(INLA_out, traj,  hilim=Inf, lolim=0, yhilim=Inf, ylolim=0, ...) {
  if (is.null(INLA_out$result)) { 
    grid <- INLA_out$grid 
    n <- length(grid) -1
    grid_pts <- (grid[2:(n+1)] + grid[1:n])/2
    
    mask = grid_pts <= hilim & grid_pts >= lolim & traj(grid_pts, ...) <= yhilim & traj(grid_pts, ...) >= lolim
    grid_pts = grid_pts[mask]
    n <- length(grid_pts)
    
    lo <- INLA_out$effpop025[mask]
    hi <- INLA_out$effpop975[mask]
  } else {
    mod = INLA_out$result$summary.random$time
    
    mask = mod$ID <= hilim & mod$ID >= lolim & traj(mod$ID, ...) <= yhilim & traj(mod$ID, ...) >= lolim
    
    grid_pts = mod$ID[mask]
    n = length(grid_pts)
    lo = exp(-mod$"0.975quant"[mask])
    hi = exp(-mod$"0.025quant"[mask])
  }
  
  
  truth = traj(grid_pts, ...) 
  result = sum(truth < hi & truth > lo)
  return(list(tot = result, avg = result / n))
}

#' Computes the MSE of a reconstructed EPS trajectory given a reference trajectory 
#'
#' @param INLA_out output directly from \code{BNPR_Lambda} or \code{posterior_beta_coal_MCMC_w_sHMC},
#'  or the \code{BNPR} argument from the output of \code{posterior_beta_coal_hybrid}
#' @param traj the reference trajectory to compare to 
#'   
#' @return List of total MSE \code{total} or averaged MSE by number of grid points \code{avg}.
#' @export
mse <- function(INLA_out, traj, hilim=Inf, lolim=0, yhilim=Inf, ylolim=0, ...) {
  if (is.null(INLA_out$result)) { 
    grid <- INLA_out$grid 
    n <- length(grid) -1
    grid_pts <- (grid[2:(n+1)] + grid[1:n])/2
    mask = grid_pts <= hilim & grid_pts >= lolim & traj(grid_pts, ...) <= yhilim & traj(grid_pts, ...) >= lolim
    grid_pts = grid_pts[mask]
    n <- length(grid_pts)
    med <- INLA_out$effpop[mask]
  } else { 
    mod = INLA_out$result$summary.random$time
    mask = mod$ID <= hilim & mod$ID >= lolim & traj(mod$ID, ...) <= yhilim & traj(mod$ID, ...) >= lolim
    
    grid_pts = mod$ID[mask]
    n = length(grid_pts)
    med = exp(-mod$"0.5quant"[mask])
  }
  
  truth = traj(grid_pts, ...)
  result = sum( (med - truth)^2/truth )
  return(list(tot = result, avg = result / n ))
}

#' Computes the biase of a reconstructed EPS trajectory given a reference trajectory 
#'
#' @param INLA_out output directly from \code{BNPR_Lambda} or \code{posterior_beta_coal_MCMC_w_sHMC},
#'  or the \code{BNPR} argument from the output of \code{posterior_beta_coal_hybrid}
#' @param traj the reference trajectory to compare to 
#'   
#' @return List of total biase \code{total} or averaged bias by number of grid points \code{avg}.
#' @export
bias_new <- function(INLA_out, traj, hilim=Inf, lolim=0, yhilim=Inf, ylolim=0, ...) {
  if (is.null(INLA_out$result)) { 
    grid <- INLA_out$grid 
    n <- length(grid) -1
    grid_pts <- (grid[2:(n+1)] + grid[1:n])/2
    mask = grid_pts <= hilim & grid_pts >= lolim & traj(grid_pts, ...) <= yhilim & traj(grid_pts, ...) >= lolim
    grid_pts = grid_pts[mask]
    n <- length(grid_pts)
    med <- INLA_out$effpop[mask]
  } else { 
    mod = INLA_out$result$summary.random$time
    mask = mod$ID <= hilim & mod$ID >= lolim & traj(mod$ID, ...) <= yhilim & traj(mod$ID, ...) >= lolim
    grid_pts = mod$ID[mask]
    n = length(grid_pts)
    med = exp(-mod$"0.5quant"[mask])
  }
  
  truth = traj(grid_pts, ...)
  
  result = sum( (med - truth)/truth )
  return(list(tot = result, avg = result / n ))
}

#' Computes the deviance of a reconstructed EPS trajectory given a reference trajectory 
#'
#' @param INLA_out output directly from \code{BNPR_Lambda} or \code{posterior_beta_coal_MCMC_w_sHMC},
#'  or the \code{BNPR} argument from the output of \code{posterior_beta_coal_hybrid}
#' @param traj the reference trajectory to compare to 
#'   
#' @return List of total deviance \code{total} or averaged deviance by number of grid points \code{avg}.
#' @export
dev_new <- function(INLA_out, traj, hilim=Inf, lolim=0, yhilim=Inf, ylolim=0, ...) {
  if (is.null(INLA_out$result)) { 
    grid <- INLA_out$grid
    n <- length(grid) -1
    grid_pts <- (grid[2:(n+1)] + grid[1:n])/2
    
    mask = grid_pts <= hilim & grid_pts >= lolim & traj(grid_pts, ...) <= yhilim & traj(grid_pts, ...) >= lolim
    grid_pts = grid_pts[mask]
    n <- length(grid_pts)
    med <- INLA_out$effpop[mask]
  } else { 
    mod = INLA_out$result$summary.random$time
    mask = mod$ID <= hilim & mod$ID >= lolim & traj(mod$ID, ...) <= yhilim & traj(mod$ID, ...) >= lolim
    
    grid_pts = mod$ID[mask]
    n = length(grid_pts)
    med = exp(-mod$"0.5quant"[mask])
  }
  
  truth = traj(grid_pts, ...)
  result = sum( abs(med - truth)/truth )
  return(list(tot = result, avg = result / n ))
}

