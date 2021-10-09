gen_INLA_args <- function(samp_times, n_sampled, coal_times)
{
  if (sum(n_sampled) != length(coal_times) + 1)
    stop("Number sampled not equal to number of coalescent events + 1.")
  
  if (length(intersect(coal_times, samp_times)) > 0)
    warning("Coincident sampling event and coalescent event: results may be unpredictable.")
  
  l <- length(samp_times)
  m <- length(coal_times)
  sorting <- sort(c(samp_times, coal_times), index.return=TRUE)
  
  lineage_change <- c(n_sampled, rep(-1, m))[sorting$ix]
  lineages <- utils::head(cumsum(lineage_change), -1) # remove entry for the post-final-coalescent-event open interval
  coal_factor <- lineages*(lineages-1)/2
  
  event <- c(rep(0, l), rep(1, m))[sorting$ix]
  
  return(list(coal_factor=coal_factor, s=sorting$x, event=event, lineages=lineages))
}

gen_summary = function(coal_times, samp_times, n_sampled)
{
  args = gen_INLA_args(coal_times, samp_times, n_sampled)
  n = length(args$s)
  return(data.frame(cbind(lineages=args$lineages, start_time=args$s[1:(n-1)], stop_time=args$s[2:n], end_event=args$event[2:n], change=diff(c(args$indicator,1)))))
}

#' Bayesian nonparametric phylodynamic reconstruction.
#' 
#' @param data \code{phylo} object or list containing vectors of coalescent 
#'   times \code{coal_times}, sampling times \code{samp_times}, and number 
#'   sampled per sampling time \code{n_sampled}.
#' @param lengthout numeric specifying number of grid points.
#' @param pref logical. Should the preferential sampling model be used?
#' @param prec_alpha,prec_beta numerics specifying gamma prior for precision 
#'   \eqn{\kappa}.
#' @param beta1_prec numeric specifying precision for normal prior on 
#'   \eqn{\beta_1}.
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
#' 
#' @examples
#' data("NY_flu")
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'  res = BNPR(NY_flu)
#'  plot_BNPR(res)
#' }
BNPR <- function(data, lengthout = 100, pref=FALSE, prec_alpha=0.01,
                 prec_beta=0.01, beta1_prec = 0.001, fns = NULL, log_fns = TRUE,
                 simplify = TRUE, derivative = FALSE, forward = TRUE,grid=NULL)
{
  if (class(data) == "phylo")
  {
    phy <- summarize_phylo(data)
  }
  else if (all(c("coal_times", "samp_times", "n_sampled") %in% names(data)))
  {
    phy <- with(data, list(samp_times = samp_times, coal_times = coal_times,
                           n_sampled = n_sampled))
  }
  
  result <- infer_coal_samp(samp_times = phy$samp_times, coal_times = phy$coal_times,
                            n_sampled = phy$n_sampled, fns = fns, lengthout = lengthout,
                            prec_alpha = prec_alpha, prec_beta = prec_beta,
                            beta1_prec = beta1_prec, use_samp = pref, log_fns = log_fns,
                            simplify = simplify, derivative = derivative,grid=grid)
  
  result$samp_times <- phy$samp_times
  result$n_sampled  <- phy$n_sampled
  result$coal_times <- phy$coal_times
  
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
 
  
  if (derivative)
  {
    if (forward)
      ind <- c(1:(lengthout-1), (lengthout-1))
    else
      ind <- c(1, 1:(lengthout-1))
    
    result$derivative <- with(result$result$summary.lincomb,
                              data.frame(time = result$x, mean = -mean[ind], sd = sd[ind],
                                         quant0.025 = -`0.975quant`[ind],
                                         quant0.5   = -`0.5quant`[ind],
                                         quant0.975 = -`0.025quant`[ind]))
  }
  
  if (pref)
  {
    result$beta0     <- result$result$summary.fixed["beta0","0.5quant"]
    result$beta0summ <- result$result$summary.fixed["beta0",]
    rownames(result$beta0summ) <- "Beta0"
    result$beta1     <- result$result$summary.hyperpar[2,"0.5quant"]
    result$beta1summ <- result$result$summary.hyperpar[2,]
    rownames(result$beta1summ) <- "Beta1"
  }
  
  return(result)
}

#' Bayesian nonparametric phylodynamic reconstruction from multiple correlated trees
#' 
#' @param data \code{multiPhylo} multiPhylo object  
#' @param lengthout numeric specifying number of grid points.
#' @param prec_alpha,prec_beta numerics specifying gamma prior for precision 
#'   \eqn{\tau}.
#' @param zero_dates, a vector with actual dates of the most recent sample per phylogeny 
#' @return Phylodynamic reconstruction of effective population size at grid 
#'   points. \code{data} contains the information passed to INLA, \code{grid} contains the grid end points, 
#'   \code{x} contains the grid point centers, \code{effpop} contains a vector 
#'   of the posterior median effective population size estimates, 
#'   \code{effpop025} and \code{effpop975} contain the 2.5th and 97.5th 
#'   posterior percentiles, \code{summary} contains a data.frame of the 
#'   estimates, and \code{derivative} (if \code{derivative = TRUE}) contains a
#'   data.frame summarizing the log-derivative.
#' @export
#' 
#' @examples
#' data("NY_flu")
#' res = BNPR(NY_flu)
#' plot_BNPR(res)
BNPR_multiple <- function(data, lengthout = 100, prec_alpha=0.01,
                 prec_beta=0.01, zero_dates=NA)
{
  if (sum(abs(ape::coalescent.intervals(data[[1]])$interval.length-ape::coalescent.intervals(data[[2]])$interval.length))==0)
    
    
  
    stop("The realizations are the same. Code needs to be fixed")
  
  if (class(data) != "multiPhylo")
  {
    result<-0
  }
  else{
    if ((is.na(zero_dates))[1]) {zero_dates<-rep(0,length(data))}
    maxdate<-max(unlist(zero_dates))
    maxcoaltimes<-0
    summarylist<-vector("list",length(data))
    for (j in 1:length(data)){
      summarylist[[j]]<-summarize_phylo(data[[j]])
      summarylist[[j]]$samp_times<-summarylist[[j]]$samp_times+maxdate-zero_dates[[j]]
      summarylist[[j]]$coal_times<-summarylist[[j]]$coal_times+maxdate-zero_dates[[j]]
      maxcoaltimes<-max(maxcoaltimes,summarylist[[j]]$coal_times)
    }
    grid <- seq(0, maxcoaltimes, length.out = lengthout+1)
    coal_data<-vector("list",length(data))
    indicators<-matrix(NA,nrow=lengthout*length(data),ncol=length(data))
    for (j in 1:length(data)){
    coal_data[[j]] <- coal_stats(grid = grid,summarylist[[j]]$samp_times, summarylist[[j]]$n_sampled,
                            coal_times = summarylist[[j]]$coal_times)
    coal_data[[j]] <- with(coal_data[[j]], condense_stats(time = time, event = event, E=E))
    
    if (j==1){y<-coal_data[[j]]$event;
              E_log<-coal_data[[j]]$E_log
              indicators[1:lengthout,1]<-coal_data[[j]]$time}else{
                y<-c(y,coal_data[[j]]$event)
                E_log<-c(E_log,coal_data[[j]]$E_log)
                indicators[(lengthout*(j-1)+1):(lengthout*j),j]<-coal_data[[j]]$time
                }
    }
   i<-indicators[,1]
   j<-indicators[,2]
   k<-indicators[,2]
   hyper.rw1 = list(prec = list(param = c(prec_alpha,prec_beta)))
   formula = y ~ f(i,model="rw1",hyper = hyper.rw1,constr = FALSE)+ f(j,copy="i",
                        hyper = list(beta=list(fixed=FALSE)))+f(k, model="iid") -1
   
   datause<-data.frame(y=y,i=i,j=j,k=k,E_log=E_log)
    
   result = inla(formula, family="Poisson",data = datause,offset=datause$E_log,control.predictor = list(compute=TRUE))
    
  #result$samp_times <- phy$samp_times
  #result$n_sampled  <- phy$n_sampled
  #result$coal_times <- phy$coal_times
  
  summarylist$effpop_1     <- exp(-result$summary.random$i$`0.5quant`)
  summarylist$effpopmean_1 <- exp(-result$summary.random$i$mean)
  summarylist$effpop975_1  <- exp(-result$summary.random$i$`0.025quant`)
  summarylist$effpop025_1  <- exp(-result$summary.random$i$`0.975quant`)
  
  summarylist$effpop_2    <- exp(-result$summary.random$j$`0.5quant`-result$summary.random$k$`0.5quant`)
  summarylist$effpopmean_2 <- exp(-result$summary.random$j$mean-result$summary.random$k$mean)
  summarylist$effpop975_2  <- exp(-result$summary.random$j$`0.025quant`-result$summary.random$k$`0.025quant`)
  summarylist$effpop025_2  <- exp(-result$summary.random$j$`0.975quant`-result$summary.random$k$`0.975quant`)
  summarylist$hypers<-result$summary.hyperpar
  summarylist$grid<-grid
  
  result<-summarylist
  
  }
  return(result)
}

#' @describeIn BNPR Uses preferential sampling model.
#' @export
BNPR_PS <- function(data, lengthout = 100, prec_alpha=0.01, prec_beta=0.01,
                    beta1_prec = 0.001, fns = NULL, log_fns = TRUE,
                    simplify = TRUE, derivative = FALSE, forward = TRUE)
{
  return(BNPR(data = data, lengthout = lengthout, pref = TRUE,
              prec_alpha = prec_alpha, prec_beta = prec_beta,
              beta1_prec = beta1_prec, fns = fns, log_fns = log_fns,
              simplify = simplify, derivative = derivative, forward = forward))
}

coal_stats <- function(grid, samp_times, coal_times, n_sampled = NULL,
                       log_zero = -100)
{
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2
  
  if (is.null(n_sampled))
    n_sampled <- rep(1, length(samp_times))
  args <- gen_INLA_args(samp_times = samp_times, n_sampled = n_sampled,
                        coal_times = coal_times)
  
  coal_factor <- args$coal_factor
  s <- args$s
  event <- args$event
  
  grid_trimmed <- setdiff(x = grid, y = s)
  sorting <- sort(c(grid_trimmed, s), index.return=TRUE)
  sgrid <- sorting$x
  ordering <- sorting$ix
  
  time_index <- cut(x = sgrid[-1], breaks = grid, labels = FALSE)
  time <- field[time_index]
  
  event_out <- c(rep(0, length(grid_trimmed)), event)[ordering]
  
  Cfun <- stats::stepfun(x = s, y = c(0, coal_factor, 0), right = TRUE)
  Cvec <- Cfun(sgrid[-1])
  E <- diff(sgrid)*Cvec
  
  E_log = log(E)
  E_log[E == 0] = log_zero
  
  return(data.frame(time = time, event = event_out[-1], E = E, E_log = E_log))
}

condense_stats <- function(time, event, E, log_zero = -100)
{
  result <- stats::aggregate(event ~ time, FUN = sum)
  result$E <- stats::aggregate(E ~ time, FUN = sum)$E
  
  E_log = log(result$E)
  E_log[result$E == 0] = log_zero
  result$E_log <- E_log
  
  return(result)
}

infer_coal <- function(samp_times, coal_times, n_sampled = NULL, lengthout = 100,
                       prec_alpha = 0.01, prec_beta = 0.01, simplify = FALSE,
                       derivative = FALSE)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE).',
         call. = FALSE)
  }
  
  if (min(coal_times) < min(samp_times))
    stop("First coalescent time occurs before first sampling time")
  
  if (max(samp_times) > max(coal_times))
    stop("Last sampling time occurs after last coalescent time")
  
  grid <- seq(min(samp_times), max(coal_times), length.out = lengthout+1)
  
  if (is.null(n_sampled))
    n_sampled <- rep(1, length(samp_times))
  
  coal_data <- coal_stats(grid = grid, samp_times = samp_times, n_sampled = n_sampled,
                          coal_times = coal_times)
  
  if (simplify)
    coal_data <- with(coal_data, condense_stats(time = time, event = event, E=E))
  
  data <- with(coal_data, data.frame(y = event, time = time, E_log = E_log))
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  formula <- y ~ -1 + f(time, model="rw1", hyper = hyper, constr = FALSE)
  
  if (derivative)
  {
    Imat <- diag(lengthout)
    A <- utils::head(Imat, -1) - utils::tail(Imat, -1)
    field <- grid[-1] - diff(grid)/2
    A <- diag(1/diff(field)) %*% A
    A[A == 0] <- NA
    
    lc_many <- INLA::inla.make.lincombs(time = A)
    
    #mod <- INLA::inla(formula, family = "poisson", data = data, lincomb = lc_many,
     #                 control.predictor = list(compute=TRUE))  
    mod <- INLA::inla(formula, family = "poisson", data = data, lincomb = lc_many,
                      control.predictor = list(compute=TRUE),
                      control.inla = list(lincomb.derived.only=FALSE))
  }
  else
  {
    mod <- INLA::inla(formula, family = "poisson", data = data, offset = data$E_log,
                      control.predictor = list(compute=TRUE))
  }
  
  return(list(result = mod, data = data, grid = grid, x = coal_data$time))
}

samp_stats <- function(grid, samp_times, n_sampled = NULL, trim_end = FALSE)
{
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2
  E <- diff(grid)
  
  bins <- cut(x = samp_times, breaks = grid, include.lowest = TRUE)
  
  if (is.null(n_sampled))
    count <- as.vector(table(bins))
  else
  {
    tab <- stats::aggregate(n_sampled ~ bins, FUN = sum, labels = FALSE)
    count <- rep(0, lengthout)
    count[as.numeric(tab$bins)] <- tab$n_sampled
  }
  
  count[utils::head(grid, -1) >= max(samp_times)] <- NA
  result <- data.frame(time = field, count = count, E = E, E_log = log(E))
  
  if (trim_end)
    result <- result[stats::complete.cases(result),]
  
  return(result)
}

joint_stats <- function(coal_data, samp_data)
{
  n1 <- length(coal_data$time)
  n2 <- length(samp_data$time)
  beta0 <- c(rep(0, n1), rep(1, n2))
  E_log <- c(coal_data$E_log, samp_data$E_log)
  Y <- matrix(c(coal_data$event, rep(NA, n2), rep(NA, n1), samp_data$count),
              nrow = n1 + n2, byrow = FALSE)
  w <- c(rep(1, n1), rep(-1, n2))
  time  <- c(coal_data$time, rep(NA, n2))
  time2 <- c(rep(NA, n1), samp_data$time)
  
  return(list(Y = Y, beta0 = beta0, time = time, time2 = time2, w = w, E_log = E_log))
}

infer_coal_samp <- function(samp_times, coal_times, n_sampled=NULL, fns = NULL,
                            lengthout=100, prec_alpha=0.01, prec_beta=0.01,
                            beta1_prec=0.001, use_samp = FALSE, log_fns = TRUE,
                            simplify = FALSE, events_only = FALSE,
                            derivative = FALSE,grid=NULL)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE).',
         call. = FALSE)
  }
  
  if (min(coal_times) < min(samp_times))
    stop("First coalescent time occurs before first sampling time")
  
  if (max(samp_times) > max(coal_times))
    stop("Last sampling time occurs after last coalescent time")
  
  if (is.null(grid))
    grid <- seq(min(samp_times), max(coal_times), length.out = lengthout+1)
  
  if (is.null(n_sampled))
    n_sampled <- rep(1, length(samp_times))
  
  coal_data <- coal_stats(grid = grid, samp_times = samp_times, n_sampled = n_sampled,
                          coal_times = coal_times)
  
  if (simplify)
    coal_data <- with(coal_data, condense_stats(time=time, event=event, E=E))
  
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  
  if (!use_samp)
  {
    data <- with(coal_data, data.frame(y = event, time = time, E_log = E_log))
    
    formula <- y ~ -1 + f(time, model="rw1", hyper = hyper, constr = FALSE, scale.model = TRUE)
    family <- "poisson"
  }
  else if (use_samp)
  {
    if (events_only)
      samp_data <- samp_stats(grid = grid, samp_times = samp_times)
    else
      samp_data <- samp_stats(grid = grid, samp_times = samp_times,
                              n_sampled = n_sampled)
    
    data <- joint_stats(coal_data = coal_data, samp_data = samp_data)
    
    if (is.null(fns))
    {
      formula <- Y ~ -1 + beta0 +
        f(time, model="rw1", hyper = hyper, constr = FALSE) +
        f(time2, w, copy="time", fixed=FALSE, param=c(0, beta1_prec))
    }
    else
    {
      vals <- NULL
      bins <- sum(data$beta0 == 0)
      for (fni in fns)
      {
        if (log_fns)
          vals <- cbind(vals, c(rep(0, bins), log(fni(samp_data$time))))
        else
          vals <- cbind(vals, c(rep(0, bins), fni(samp_data$time)))
      }
      data$fn <- vals
      
      formula <- Y ~ -1 + beta0 + fn +
        f(time, model="rw1", hyper = hyper, constr = FALSE) +
        f(time2, w, copy="time", fixed=FALSE, param=c(0, beta1_prec))
    }
    
    family <- c("poisson", "poisson")
  }
  else
    stop("Invalid use_samp value, should be boolean.")
  
  if (derivative)
  {
    Imat <- diag(lengthout)
    A <- utils::head(Imat, -1) - utils::tail(Imat, -1)
    field <- grid[-1] - diff(grid)/2
    A <- diag(1/diff(field)) %*% A
    A[A == 0] <- NA
    
    lc_many <- INLA::inla.make.lincombs(time = A)
  }
  else
  {
    lc_many <- NULL
  }
  
  mod <- INLA::inla(formula, family = family, data = data,
                    lincomb = lc_many, offset = data$E_log,
                    control.predictor = list(compute=TRUE),
                    control.compute=list(config = TRUE))
                    
  #mod <- INLA::inla(formula, family = family, data = data,
   #                 lincomb = lc_many, offset = data$E_log,
    #                control.predictor = list(compute=TRUE),
     #               control.inla = list(lincomb.derived.only=FALSE))
  
  return(list(result = mod, data = data, grid = grid, x = coal_data$time))
}

infer_coal_deriv <- function(samp_times, coal_times, n_sampled = NULL, lengthout = 100,
                             prec_alpha = 0.01, prec_beta = 0.01, simplify = FALSE)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE).',
         call. = FALSE)
  }
  
  if (min(coal_times) < min(samp_times))
    stop("First coalescent time occurs before first sampling time")
  
  if (max(samp_times) > max(coal_times))
    stop("Last sampling time occurs after last coalescent time")
  
  grid <- seq(min(samp_times), max(coal_times), length.out = lengthout+1)
  
  if (is.null(n_sampled))
    n_sampled <- rep(1, length(samp_times))
  
  coal_data <- coal_stats(grid = grid, samp_times = samp_times, n_sampled = n_sampled,
                          coal_times = coal_times)
  
  if (simplify)
    coal_data <- with(coal_data,
                      condense_stats(time = time, event = event, E=E))
  
  data <- with(coal_data, data.frame(y = event, time = time, E_log = E_log))
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  formula <- y ~ -1 + f(time, model="rw1", hyper = hyper, constr = FALSE) + offset(data$E_log)
  
  Imat <- diag(lengthout)
  A <- utils::head(Imat, -1) - utils::tail(Imat, -1)
  A[A == 0] <- NA
  
  #lcmat <- rbind(c(1, -1, rep(NA, lengthout - 2)), c(NA, 1, -1, rep(NA, lengthout - 3)))
  lc_many <- INLA::inla.make.lincombs(time = A)
  #names(lc1) = "lc1"
  
  mod <- INLA::inla(formula, family = "poisson", data = data, lincomb = lc_many,
                    control.predictor = list(compute=TRUE),
                    control.inla = list(lincomb.derived.only=FALSE))
  
  return(list(result = mod, data = data, grid = grid, x = coal_data$time))
}

infer_samp_exper <- function(samp_times, fns, n_sampled = NULL, lengthout = 100,
                             prec_alpha = 0.01, prec_beta = 0.01)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE).',
         call. = FALSE)
  }
  
  grid <- seq(min(samp_times),max(samp_times),length.out=lengthout+1)
  
  samp_data <- samp_stats(grid = grid, samp_times = samp_times,
                          n_sampled = n_sampled)
  
  vals = NULL
  for (fni in fns)
  {
    vals = cbind(vals, log(fni(samp_data$time)))
  }
  
  data <- with(samp_data, list(y = count, time = time, E_log = E_log))
  data$fn=vals
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  formula_sampling <- y ~ 1 + fn + f(time, model="rw1", hyper = hyper, constr=FALSE)
  
  mod <- INLA::inla(formula_sampling, family="poisson", data=data,
                    offset=data$E_log, control.predictor=list(compute=TRUE))
  
  return(list(result = mod, data = data, grid = grid, x = samp_data$time))
}

infer_coal_samp_exper <- function(samp_times, coal_times, n_sampled=NULL, fns = NULL,
                                  lengthout=100, prec_alpha=0.01, prec_beta=0.01,
                                  beta1_prec=0.001, use_samp = FALSE, log_fns = TRUE,
                                  simplify = FALSE, events_only = FALSE)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE).',
         call. = FALSE)
  }
  
  if (min(coal_times) < min(samp_times))
    stop("First coalescent time occurs before first sampling time")
  
  if (max(samp_times) > max(coal_times))
    stop("Last sampling time occurs after last coalescent time")
  
  grid <- seq(min(samp_times), max(coal_times), length.out = lengthout+1)
  
  if (is.null(n_sampled))
    n_sampled <- rep(1, length(samp_times))
  
  coal_data <- coal_stats(grid = grid, samp_times = samp_times, n_sampled = n_sampled,
                          coal_times = coal_times)
  
  if (simplify)
    coal_data <- with(coal_data, condense_stats(time=time, event=event, E=E))
  
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  
  if (events_only)
    samp_data <- samp_stats(grid = grid, samp_times = samp_times)
  else
    samp_data <- samp_stats(grid = grid, samp_times = samp_times,
                            n_sampled = n_sampled)
  
  data <- joint_stats(coal_data = coal_data, samp_data = samp_data)
  
  vals <- NULL
  bins <- sum(data$beta0 == 0)
  for (fni in fns)
  {
    vals <- cbind(vals, c(rep(0, bins), log(fni(samp_data$time))))
  }
  data$fn <- vals
  
  formula <- Y ~ -1 + beta0 + fn +
    f(time, model="rw1", hyper = hyper, constr = FALSE) +
    f(time2, w, copy="time", fixed=FALSE, param=c(0, beta1_prec))
  
  family <- c("poisson", "poisson")
  
  mod <- INLA::inla(formula, family = family, data = data, offset = data$E_log,
                    control.predictor = list(compute=TRUE),
                    control.inla = list(lincomb.derived.only=FALSE))
  
  return(list(result = mod, data = data, grid = grid, x = coal_data$time))
}

infer_samp <- function(samp_times, n_sampled = NULL, lengthout = 100, grid = NULL,
                       prec_alpha = 0.01, prec_beta = 0.01)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE).',
         call. = FALSE)
  }
  
  if (is.null(grid)) {
    grid <- seq(min(samp_times),max(samp_times),length.out=lengthout+1)
  }
  
  samp_data <- samp_stats(grid = grid, samp_times = samp_times,
                          n_sampled = n_sampled)
  
  data <- with(samp_data, data.frame(y = count, time = time, E_log = E_log))
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  formula_sampling <- y ~ 1 + f(time, model="rw1", hyper = hyper, constr=FALSE)
  
  mod <- INLA::inla(formula_sampling, family="poisson", data=data,
                    offset=data$E_log, control.predictor=list(compute=TRUE))
  
  return(list(result = mod, data = data, grid = grid, x = samp_data$time))
}

#' Nonparametric estimate of sampling intensity
#' 
#' @param data \code{phylo} object or list containing vectors of sampling times
#'   \code{samp_times} and number sampled per sampling time \code{n_sampled}.
#' @param lengthout numeric specifying number of grid points.
#' @param grid numeric vector of endpoints of the latent field.
#' @param prec_alpha,prec_beta numerics specifying gamma prior for precision 
#'   \eqn{\kappa}.
#'   
#' @return
#' @export
#' 
#' @examples
BNPR_samp_only <- function(data, lengthout = 100, grid = NULL, tmrca = NULL,
                           prec_alpha = 0.01, prec_beta = 0.01)
{
  if (class(data) == "phylo")
  {
    phy <- summarize_phylo(data)
  }
  else if (all(c("coal_times", "samp_times", "n_sampled") %in% names(data)))
  {
    phy <- with(data, list(samp_times = samp_times, coal_times = coal_times,
                           n_sampled = n_sampled))
  }
  
  if (is.null(grid)) {
    if (is.null(tmrca)) {
      grid <- seq(min(phy$samp_times), max(phy$samp_times), length.out=lengthout+1)
    } else {
      grid <- seq(min(phy$samp_times), tmrca, length.out=lengthout+1)
    }
  }
  midpts <- grid[-1] - diff(grid)/2
  
  foo <- infer_samp(samp_times = phy$samp_times, n_sampled = phy$n_sampled, 
                       lengthout = lengthout, grid = grid,
                       prec_alpha = prec_alpha, prec_beta = prec_beta)
  
  quants <- foo$result$summary.random$time
  beta0 <- foo$result$summary.fixed$`0.5quant`
  
  result <- list(q025Samp = exp(quants$`0.025quant` + beta0),
                medianSamp = exp(quants$`0.5quant` + beta0),
                q975Samp = exp(quants$`0.975quant` + beta0),
                beta0 = beta0, grid = grid, midpts = midpts,
                samp_times = gene$samp_times, n_sampled = gene$n_sampled, internals = foo)
  
  return(result)
}

log_samp_int = function(logpop, betas,
                        covariates = NULL, cov_betas = NULL,
                        pow_covariates = NULL, pow_cov_betas = NULL)
{
  result = betas[1] + logpop * betas[2]
  
  if (!is.null(covariates))
  {
    for (i in 1:length(covariates))
    {
      result = result + covariates[[i]] * cov_betas[i]
    }
  }
  
  if (!is.null(pow_covariates))
  {
    for (i in 1:length(pow_covariates))
    {
      result = result + pow_covariates[[i]] * logpop * pow_cov_betas[i]
    }
  }
  
  return(result)
}

#' Calculate posterior sampling intensities
#' 
#' @param popTbl table of log-effective population sizes.
#' @param betaTbl table of log-linear sampling model coefficients.
#' @param covariates list of covariates
#' @param powBetaTbl table of log-linear sampling model coefficients for terms 
#'   crossed with the log-effective population size.
#' @param pow_covariates list of covariates to be crossed with the log-effective
#'   population size.
#'   
#' @return a matrix of posterior sampling intensities.
#' @export
#' 
#' @examples 2 + 2
log_samp_mat <- function(popTbl, betaTbl, covariates = NULL, 
                         powBetaTbl = NULL, pow_covariates = NULL)
{
  n <- nrow(popTbl)
  m <- ncol(popTbl)
  #m <- length(popTbl[1,])
  #grid <- epoch_width * 0:m
  
  samp_mat<- matrix(0, nrow = n, ncol = m)
  for (i in 1:n) {
    logpop <- as.numeric(popTbl[i, ])
    betas <- as.numeric(betaTbl[i, ])
    
    cov_betas <- NULL
    if (length(betas) > 2) {
      cov_betas <- tail(betas, -2)
      betas <- betas[1:2]
    }
    
    pow_cov_betas <- NULL
    if (ncol(powBetaTbl) > 0) {
      pow_cov_betas <- as.numeric(powBetaTbl[i, ])
    }
    
    samp_mat[i, ] <- log_samp_int(logpop = logpop, betas = betas,
                                  covariates = covariates, cov_betas = cov_betas,
                                  pow_covariates = pow_covariates, pow_cov_betas = pow_cov_betas)
  }
  
  return(samp_mat)
}


##Functions for infering ``selection" coefficient 

#' Bayesian nonparametric phylodynamic reconstruction of two EPSs, the two trees share baseline Ne
#' 
#' @param \code{tree1} and \code{tree2} are two phylo objects 
#'   sampling times \code{samp_times1} and  \code{samp_times2} encodes sampling information
#' @param lengthout numeric specifying number of grid points.
#' @param prec_alpha,prec_beta numerics specifying gamma prior for precision 
#'   \eqn{\kappa}.
#' @param beta1_prec numeric specifying precision for normal prior on 
#'   \eqn{\beta_1}.
#' @param beta0_remove wheter to have an intercept 
#'   \eqn{\beta_1}.
#' @param simplify logical whether to fully bucket all Poisson points.
#' @param derivative logical whether to calculate estimates of the 
#'   log-derivative.
#' @param forward logical whether to use the finite difference approximations of
#'   the log-derivative as a forward or backward derivative.
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
#' 
#' @examples
#' data("NY_flu")
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'  res = BNPR(NY_flu)
#'  plot_BNPR(res)
#' }
BNPR_sel <- function (tree1,tree2, samp_times1,samp_times2, lengthout = 100, prec_alpha = 0.01, 
                      prec_beta = 0.01, beta1_prec = 0.001, beta0_remove=FALSE,
                      simplify = TRUE, derivative = FALSE, forward = TRUE) 
{
  
  phy1 <- summarize_phylo(tree1)
  phy1$samp_times <- phy1$samp_times + min(samp_times1)
  phy1$coal_times <- phy1$coal_times + min(samp_times1)
  phy2 <- summarize_phylo(tree2)
  phy2$samp_times <- phy2$samp_times + min(samp_times2)
  phy2$coal_times <- phy2$coal_times + min(samp_times2)
  
  
  result <- infer_coal_samp_selection(phy1,phy2,lengthout = lengthout, 
                                      prec_alpha = prec_alpha, prec_beta = prec_beta, beta1_prec = beta1_prec, 
                                      simplify,beta0_remove)
  
  #result$samp_times <- phy$samp_times
  #result$n_sampled <- phy$n_sampled
  #result$coal_times <- phy$coal_times
  result$effpop <- exp(-result$result$summary.random$time$`0.5quant`)
  result$effpopmean <- exp(-result$result$summary.random$time$mean)
  result$effpop975 <- exp(-result$result$summary.random$time$`0.025quant`)
  result$effpop025 <- exp(-result$result$summary.random$time$`0.975quant`)
  result$summary <- with(result$result$summary.random$time, 
                         data.frame(time = ID, mean = exp(-mean), sd = sd * exp(-mean), 
                                    quant0.025 = exp(-`0.975quant`), quant0.5 = exp(-`0.5quant`), 
                                    quant0.975 = exp(-`0.025quant`)))
  
  if (beta0_remove==FALSE){
  result$beta0 <- result$result$summary.fixed["beta0", "0.5quant"]
  result$beta0summ <- result$result$summary.fixed["beta0",]
  rownames(result$beta0summ) <- "Beta0"
  result$beta0post <- result$result$marginals.fixed$beta0
  }
  result$beta1 <- result$result$summary.hyperpar[2, "0.5quant"]
  result$beta1summ <- result$result$summary.hyperpar[2,]
  rownames(result$beta1summ) <- "Beta1"
  result$beta1post <- result$result$marginals.hyperpar$`Beta for time2`
  
  return(result)
}





infer_coal_samp_selection <- function(phy1,phy2, lengthout=100, prec_alpha=0.01, prec_beta=0.01, 
                                      beta1_prec=0.001, simplify = TRUE, beta0_remove=FALSE,
                                      events_only = FALSE)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE).',
         call. = FALSE)
  }
  
  
  coal_times1 <- phy1$coal_times
  coal_times2 <- phy2$coal_times
  n_sampled1 <- phy1$n_sampled
  n_sampled2 <- phy2$n_sampled
  samp_times1 <- phy1$samp_times
  samp_times2 <- phy2$samp_times
  
  grid <- seq(min(c(samp_times1,samp_times2)), max(c(coal_times1,coal_times2)), length.out = lengthout+1)
  
  #if (is.null(n_sampled))
  #  n_sampled <- rep(1, length(samp_times))
  
  coal_data1 <- coal_stats(grid = grid, samp_times = samp_times1, n_sampled = n_sampled1,
                           coal_times = coal_times1)
  if (simplify){coal_data1 <- with(coal_data1, condense_stats(time=time, event=event, E=E))}
  coal_data1 <- identify_off_grid(coal_data1,samp_times1,coal_times1)
  
  coal_data2 <- coal_stats(grid = grid, samp_times = samp_times2, n_sampled = n_sampled2,
                           coal_times = coal_times2)
  if (simplify){coal_data2 <- with(coal_data2, condense_stats(time=time, event=event, E=E))}
  coal_data2 <- identify_off_grid(coal_data2,samp_times2,coal_times2)
  
  
  
  #simplify is used to avoid duplicate in time 
  
  
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  
  
  
  data <- joint_coal_stats(coal_data1 = coal_data1, coal_data2 = coal_data2)
  
  
  if (beta0_remove){
    formula <- Y ~ -1 + 
      f(time, model="rw1", hyper = hyper, constr = FALSE) +
      f(time2,copy="time", fixed=FALSE, param=c(0, beta1_prec))
  } else {
    formula <- Y ~ -1 + beta0 +
      f(time, model="rw1", hyper = hyper, constr = FALSE) +
      f(time2,copy="time", fixed=FALSE, param=c(0, beta1_prec))
  }
  
  family <- c("poisson", "poisson")
  
  lc_many <- NULL
  
  mod <- INLA::inla(formula, family = family, data = data,
                    lincomb = lc_many, offset = data$E_log,
                    control.predictor = list(compute=TRUE),
                    control.inla = list(lincomb.derived.only=FALSE))
  
  return(list(result = mod, data = data, grid = grid, x = coal_data1$time))
}



identify_off_grid <- function(coal_data,samp_times,coal_times)
{ 
  if (samp_times[1]>0){
    id <- which.min(abs(coal_data$time-samp_times[1]))
    coal_data$event[1:(id-1)] <- NA
  }
  if (coal_times[length(coal_times)]<coal_data$time[length(coal_data$time)]){
    id <- which.min(abs(coal_data$time-coal_times[length(coal_times)])) #Note that there are double points in $time, that's why id+1
    coal_data$event[(id+1):length(coal_data$event)] <- NA
  }
  return(coal_data)
}


joint_coal_stats <- function(coal_data1, coal_data2)
{
  n1 <- length(coal_data1$time) #Here are the midpts 
  n2 <- length(coal_data2$time) #
  beta0 <- c(rep(0, n1), rep(1, n2))
  E_log <- c(coal_data1$E_log, coal_data2$E_log) #samp_data$E_log does not include NA in the pref_samp, so I am not modifying it here either
  Y <- matrix(c(coal_data1$event, rep(NA, n2), rep(NA, n1), coal_data2$event),
              nrow = n1 + n2, byrow = FALSE)
  #Need to correct Y for the parts of the grid that need to be eccluded
  #Note, there are some parts of Y where they are both NAs
  w <- c(rep(1, n1), rep(1, n2))  #what does the w do? I think that takes into account for the fact that one is at the numerator,
  #the other one at the numerator. 
  
  #we are not correcting the time, i.e. there are times also to part that are not defined
  time  <- c(coal_data1$time, rep(NA, n2))
  time2 <- c(rep(NA, n1), coal_data2$time)
  
  return(list(Y = Y, beta0 = beta0, time = time, time2 = time2, w = w, E_log = E_log))
}


