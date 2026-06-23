#### Data Simulation
generate_true_M_and_data <- function(num_tips, traj = "exp_traj", rate = 1, seq_len = 1000,
                                     write_files = FALSE,
                                     name_fasta = "sequences.fasta",
                                     name_tree  = "true_tree.newick") {
  samp_times <- 0
  M_true_tree <- generate_newick(
    coalsim(
      samp_times = samp_times,
      n_sampled  = num_tips,
      traj       = traj,
      method     = "tt",
      val_upper  = 11
    )
  )
  M_true_tree$newick$tip.label <- sapply(M_true_tree$newick$tip.label,
                                         function(x) sub("_0", "", x))
  M_true_tree$labels <- sapply(M_true_tree$labels,
                               function(x) sub("_0", "", x))
  M_true_newick_string <- write.tree(M_true_tree$newick)

  temp_fasta <- tempfile(fileext = ".fasta")
  fasta_path <- if (write_files) name_fasta else tempfile(fileext = ".fasta")
  on.exit({
    unlink(temp_fasta)
    if (!write_files) unlink(fasta_path)
  }, add = TRUE)

  opts <- paste0("-mHKY -t2.0 -f0.25,0.25,0.25,0.25 -l", seq_len, " -s", rate)
  seqgen(opts = opts, newick.tree = M_true_newick_string, temp.file = temp_fasta)
  data2 <- read.phylip(temp_fasta)
  dat2fasta(data2, outfile = fasta_path)
  mydna <- as.phyDat(read.FASTA(fasta_path))

  if (write_files) {
    write.tree(M_true_tree$newick, file = name_tree)
  }

  list(M_true_tree = M_true_tree$newick, sequences = mydna)
}

#### Normalizing Constant Estimation
generate_sample_uniform <- function(num_samps, num_tips = 10) {
  encoding_chain <- rEncod(m = num_samps, n = num_tips, distr = "uniform")
  samples <- vector("list", num_samps)
  for (i in 1:num_samps) {
    samples[[i]] <- Fmat_from_myencod(encoding_chain[[i]])
  }
  samples
}

precompute_tree_chain_distance_cache <- function(tree_chain) {
  valid_tree_chain <- Filter(Negate(is.null), tree_chain)
  if (length(valid_tree_chain) == 0) {
    stop("tree_chain must contain at least one matrix.")
  }

  template <- valid_tree_chain[[1]]
  template_dim <- dim(template)
  if (length(template_dim) != 2 || template_dim[1] != template_dim[2]) {
    stop("tree_chain entries must be square matrices.")
  }

  lower_idx <- lower.tri(template, diag = TRUE)
  feature_count <- sum(lower_idx)

  tree_matrix <- t(vapply(
    valid_tree_chain,
    function(tree) {
      if (!identical(dim(tree), template_dim)) {
        stop("All tree matrices must have the same dimensions.")
      }
      tree[lower_idx]
    },
    numeric(feature_count)
  ))

  list(
    tree_matrix = tree_matrix,
    tree_ss = rowSums(tree_matrix * tree_matrix),
    lower_idx = lower_idx,
    template_dim = template_dim,
    n_trees = nrow(tree_matrix)
  )
}

tree_chain_squared_distances <- function(cache, M) {
  if (is.null(cache$tree_matrix) || is.null(cache$tree_ss) || is.null(cache$lower_idx)) {
    stop("cache must come from precompute_tree_chain_distance_cache().")
  }

  if (!identical(dim(M), cache$template_dim)) {
    stop("M must have the same dimensions as the cached tree matrices.")
  }

  M_vec <- M[cache$lower_idx]
  M_ss <- sum(M_vec * M_vec)

  as.numeric(cache$tree_ss - 2 * as.vector(cache$tree_matrix %*% M_vec) + M_ss)
}

zigzag <- function(n) {
  fact<-list()
  zig<-list()
  fact[[1]] = 1;
  for (i in (1:n)){
    fact[[i+1]] = fact[[i]] *i
  }
  zig[[1]] = 1;
  zig[[2]] = 1;
  for (i in (2:(n-1))){
    sum = 0
    for (k in (0:(i-1))){
      sum = sum+(fact[[i]]/(fact[[i- k]]*fact[[k+1]]))*zig[[k+1]]*zig[[i-k]]
    }
    zig[[i+1]] = sum / 2
  }
  return(zig[[length(zig)]])
}


compute_log_Z_est <- function(beta, M, cache, diam = 1) {
  num_tips <- dim(M)[1] + 1
  distances <- tree_chain_squared_distances(cache, M)


  log_terms <- -(beta / diam) * distances
  if (any(!is.finite(log_terms))) {
    stop(
      "compute_log_Z_est() produced non-finite values: ",
      "beta=", beta,
      ", diam=", diam

    )
  }

  m <- max(log_terms)

  log_Z_est <- log(zigzag(num_tips - 1)) + m + log(mean(exp(log_terms - m)))
  if (!is.finite(log_Z_est)) {
    stop(
      "compute_log_Z_est() returned a non-finite log_Z estimate: ",
      "beta=", beta,
      ", diam=", diam
    )
  }

  log_Z_est
}


###Likelihood Function

### REPLACE WITH NEW FUNCTION USING UPGMA STYLE THING
log_likelihood_given_tree <- function(tree_fmat, coal_times, sequences, mode = "max", R = 1) {
  best_ll <- -Inf
  best_tree <- NULL
  ll_sum <- 0
  for (r in seq_len(R)) {
    rooted_tree <- mytree_from_F(tree_fmat, coal_times)
    ll <- pml(
      rooted_tree,
      sequences,
      bf = c(0.25, 0.25, 0.25, 0.25),
      Q  = c(1, 2, 1, 1, 2, 1)
    )$log
    ll_sum <- ll_sum + ll
    if (ll > best_ll) {
      best_ll <- ll
      best_tree <- rooted_tree
    }
  }
  log_likelihood <- if (mode == "max") best_ll else ll_sum / R
  list(log_likelihood = log_likelihood, rooted_tree = best_tree)
}


### GRADIENT ESTIMATORS

estimate_grad_g2 <- function(data,
                             coal_times,
                             samples_gibbs,
                             cache_gibbs,
                             M_est,
                             b_est,
                             log_Z,
                             num_tip_label_iters) {
  num_tips <- dim(M_est)[1] + 1
  expectation <- 0
  elbo <- 0
  expected_sq_dist <- mean(tree_chain_squared_distances(cache_gibbs, M_est))
  for (mat in samples_gibbs) {
    l_lik <- log_likelihood_given_tree(sequences = data,
                                       tree_fmat = mat,
                                       coal_times = coal_times,
                                       mode = "average",
                                       R = num_tip_label_iters)$log_likelihood
    l_var <- -(num_tips - 1 - num_cherries(mat)) * log(2) +
      lfactorial(num_tips - 1) - log_Z - b_est * norm(M_est - mat, type = "F")^2
    grad_g_log_q <- b_est * (expected_sq_dist - norm(mat - M_est, type = "F")^2)
    expectation <- expectation + grad_g_log_q * (l_lik - l_var)
    elbo <- elbo + (l_lik - l_var)
  }
  list(grad = expectation, elbo = elbo)
}

estimate_grad_M2 <- function(data,
                             coal_times,
                             samples_gibbs,
                             M_est,
                             b_est,
                             log_Z,
                             num_tip_label_iters) {
  num_tips <- dim(M_est)[1] + 1
  expectation <- 0
  elbo <- 0
  mean_f <- Reduce('+', samples_gibbs) / length(samples_gibbs)
  for (mat in samples_gibbs) {
    l_lik <- log_likelihood_given_tree(sequences = data,
                                       tree_fmat = mat,
                                       coal_times = coal_times,
                                       mode = "average",
                                       R = num_tip_label_iters)$log_likelihood
    grad_M_log_q <- 2 * b_est * (mat - mean_f)
    l_var <- -(num_tips - 1 - num_cherries(mat)) * log(2) +
      lfactorial(num_tips - 1) - log_Z - b_est * norm(M_est - mat, type = "F")^2
    expectation <- expectation + grad_M_log_q * (l_lik - l_var)
    elbo <- elbo + (l_lik - l_var)
  }
  list(grad = expectation, elbo = elbo)
}



estimate_grad_M_g <- function(data,
                              coal_times,
                              samples_gibbs,
                              cache_gibbs,
                              M_est,
                              b_est,
                              log_Z,
                              num_tip_label_iters) {
  num_tips <- dim(M_est)[1] + 1
  expectation_M <- 0
  expectation_g <- 0
  elbo <- 0
  mean_f <- Reduce('+', samples_gibbs) / length(samples_gibbs)
  expected_sq_dist <- mean(tree_chain_squared_distances(cache_gibbs, M_est))
  for (mat in samples_gibbs) {
    l_lik <- log_likelihood_given_tree(sequences = data,
                                       tree_fmat = mat,
                                       coal_times = coal_times,
                                       mode = "average",
                                       R = num_tip_label_iters)$log_likelihood
    grad_M_log_q <- 2 * b_est * (mat - mean_f)
    l_var <- -(num_tips - 1 - num_cherries(mat)) * log(2) +
      lfactorial(num_tips - 1) - log_Z - b_est * norm(M_est - mat, type = "F")^2
    expectation_M <- expectation_M + grad_M_log_q * (l_lik - l_var)
    grad_g_log_q <- b_est * (expected_sq_dist - norm(mat - M_est, type = "F")^2)
    expectation_g <- expectation_g + grad_g_log_q * (l_lik - l_var)
    elbo <- elbo + (l_lik - l_var)
  }
  list(grad_g = expectation_g, grad_M = expectation_M, elbo = elbo)
}

#### Sampling Utilities

Logdistance_prob <- function(fmat, M, beta=1, diam=1, N=1){
  n <- dim(fmat)[1] + 1
  -(beta/diam) * distance_Fmat(fmat, M, dist="l2")^2/N^4
}


sampleF <- function(M, beta, iter, diam, startF = NULL, take_every = 1){
  acceptance_rate <- 0
  n <- nrow(M) + 1
  
  chainF <- vector("list", iter/take_every)
  
  if (is.null(startF)) {
    current_tree <- rcoal(n)
    chainF[[1]]  <- gen_Fmat(current_tree, tol=13)
  } else {
    chainF[[1]] <- startF
  }
  
  current_Encod <- my_encod(chainF[[1]])
  f_current     <- Fmat_from_myencod(current_Encod)
  Logdist_curr  <- Logdistance_prob(f_current, M, beta=beta, diam=diam)
  
  MeanM <- matrix(0, nrow=n-1, ncol=n-1)
  
  
  for (j in 1:iter){
    
    proposed    <- proposal_myencod(current_Encod)
    f_proposed  <- Fmat_from_myencod(proposed)
    Logdist_prop <- Logdistance_prob(f_proposed, M, beta=beta, diam=diam)
    
    prob <- exp(Logdist_prop - Logdist_curr)
    
    if (runif(1) < prob){
      if(j%%take_every == 0){
        chainF[[j/take_every]] <- f_proposed
      }
      
      current_Encod   <- proposed
      f_current       <- f_proposed
      Logdist_curr    <- Logdist_prop
      acceptance_rate <- acceptance_rate + 1
    } else {
      if(j%%take_every == 0){
        chainF[[j/take_every]] <- f_current
      }
    }
    if(j%%take_every == 0){
      MeanM <- MeanM + chainF[[j/take_every]]
    }

  }

  list(
    chainF = chainF,
    MeanM = MeanM / length(chainF),
    acceptance_rate = acceptance_rate / iter
  )
}

