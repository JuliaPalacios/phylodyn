#library(ape)
#library(phylodyn)
#library(phyclust)
#library(phangorn)
#library(phylotools)
#library(fmatrix)

#source("R/VB_Utils.R")



# ---- Config ----
num_tips            <- 10
step_size           <- 0.1
num_samps           <- 200
num_unif_samples    <- 10000
num_grad_desc_steps <- 50
num_tip_label_iters <- 30
take_every          <- 20
burnin              <- ceiling(round(num_samps * 0.1) / take_every) * take_every
joint               <- TRUE
init_method         <- "upgma"  # "caterpillar" or "upgma"

# ---- Generate data + initialize ----
init_results <- generate_true_M_and_data(num_tips)

inter_coal_times <- coalescent.intervals(init_results$M_true_tree)$interval.length
inter_coal_times[inter_coal_times <= 0.001] <- 0.01
coal_times <- cumsum(inter_coal_times)

if (init_method == "upgma") {
  upgma_tree <- upgma(dist.hamming(init_results$sequences))
  M_est_tree <- update_time(upgma_tree, coal_times)
} else {
  init_Fmat  <- gen_caterpillar(num_tips)
  M_est_tree <- mytree_from_F(init_Fmat, coal_times)
}
M_est <- round(gen_Fmat(M_est_tree, tol = 8), 0)
M_init     <- M_est
M_true     <- round(gen_Fmat(init_results$M_true_tree, tol = 8), 0)
g_est      <- log(0.01)

# Uniform proposal samples + cache for log-Z — depends only on num_tips, build once.
# For num_tips in 5..11 the fmatrix package ships F.list<n> with every F-matrix
# enumerated, so use that exact list instead of sampling.
flist_name <- paste0("F.list", num_tips)
samples_uniform <- if (exists(flist_name, mode = "list")) {
  get(flist_name)
} else {
  generate_sample_uniform(num_samps = num_unif_samples, num_tips = num_tips)
}
uniform_cache <- precompute_tree_chain_distance_cache(samples_uniform)

# ---- AdaGrad-style gradient descent ----
M_grads <- list(matrix(step_size, nrow = nrow(M_est), ncol = nrow(M_est)))
g_grads <- list(step_size)
elbos2  <- c()

start_time <- Sys.time()
for (i in 1:num_grad_desc_steps) {
  print(paste0("Gradient Descent Iteration ", i))
  print(M_est)
  print(paste0("Estimated log beta: ", g_est))

  samples <- generate_sample_markov_chain(
    M = M_est, b = exp(g_est),
    num_samps = num_samps, burnin = burnin, take_every = take_every
  )
  log_Z <- compute_log_Z_est(beta = exp(g_est), M = M_est, cache = uniform_cache)

  if (!joint) {
    cache_gibbs <- precompute_tree_chain_distance_cache(samples)
    g_up <- estimate_grad_g2(
      data                = init_results$sequences,
      coal_times          = coal_times,
      samples_gibbs       = samples,
      cache_gibbs         = cache_gibbs,
      M_est               = M_est,
      b_est               = exp(g_est),
      log_Z               = log_Z,
      num_tip_label_iters = num_tip_label_iters
    )
    g_grad <- g_up$grad
    g_grads[[length(g_grads) + 1]] <- g_grad^2
    g_est <- g_est + step_size / sqrt(Reduce('+', g_grads)) * g_grad

    samples <- generate_sample_markov_chain(
      M = M_est, b = exp(g_est),
      num_samps = num_samps, burnin = burnin
    )
    log_Z <- compute_log_Z_est(beta = exp(g_est), M = M_est, cache = uniform_cache)

    M_up <- estimate_grad_M2(
      data                = init_results$sequences,
      coal_times          = coal_times,
      samples_gibbs       = samples,
      M_est               = M_est,
      b_est               = exp(g_est),
      log_Z               = log_Z,
      num_tip_label_iters = num_tip_label_iters
    )
    M_grad <- M_up$grad
    M_grads[[length(M_grads) + 1]] <- M_grad^2
    M_est <- M_est + step_size / sqrt(Reduce('+', M_grads)) * M_grad
    elbos2 <- c(elbos2, M_up$elbo)
  } else {
    cache_gibbs <- precompute_tree_chain_distance_cache(samples)
    J_up <- estimate_grad_M_g(
      data                = init_results$sequences,
      coal_times          = coal_times,
      samples_gibbs       = samples,
      cache_gibbs         = cache_gibbs,
      M_est               = M_est,
      b_est               = exp(g_est),
      log_Z               = log_Z,
      num_tip_label_iters = num_tip_label_iters
    )
    M_grad <- J_up$grad_M
    M_grads[[length(M_grads) + 1]] <- M_grad^2
    M_est <- M_est + step_size / sqrt(Reduce('+', M_grads)) * M_grad

    g_grad <- J_up$grad_g
    g_grads[[length(g_grads) + 1]] <- g_grad^2
    g_est <- g_est + step_size / sqrt(Reduce('+', g_grads)) * g_grad

    elbos2 <- c(elbos2, J_up$elbo)
  }
}
end_time <- Sys.time()
print(end_time - start_time)

plot(elbos2)

M_estimated_tree <- mytree_from_F(nearby_Fmat(M_est), coal_times)
plot(M_estimated_tree)
plot(init_results$M_true_tree)

print(paste0("L2 distance to true M:    ", distance_Fmat(M_est, M_true, dist = "l2")))
print(paste0("L2 distance to initial M: ", distance_Fmat(M_est, M_init, dist = "l2")))
