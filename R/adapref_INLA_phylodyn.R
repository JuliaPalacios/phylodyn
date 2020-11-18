#########Function to Adaptive preferential sampling with INLA



BNPR_PS_ada<-function (data, lengthout = 100, zeta1 = 0.01, zeta2 = 0.01, simplify = TRUE,
                       derivative = FALSE, forward = TRUE)
{
  return(BNPR_ada(data = data, lengthout = lengthout, pref = TRUE,
                  zeta1 = zeta1, zeta2 = zeta2,  simplify = simplify, derivative = derivative,
                  forward = forward))
}


#' Adaptive preferential sampling method_INLA approximation
#'
#' @param data \code{phylo} object or list containing vectors of coalescent
#'   times \code{coal_times}, sampling times \code{samp_times}, and number
#'   sampled per sampling time \code{n_sampled}.
#' @param lengthout number of grid points.
#' @param zeta1 precision parameter field on effective pop size
#' @param zeta2 precision parameter field on beta
#' @param simplify logical whether to fully bucket all Poisson points.
#' @param derivative logical whether to calculate estimates of the
#'   log-derivative.
#' @param forward logical whether to use the finite difference approximations of
#'   the log-derivative as a forward or backward derivative.
#'
#' @return Effective pop size and beta posterior,
#'  \code{result} contains the INLA output, \code{data} contains the
#'   information passed to INLA, \code{grid}  grid boundaries,
#'   \code{x}  grid mid points \code{effpop} median effective population size estimates,
#'   \code{effpop025} and \code{effpop975}  2.5th and 97.5th
#'   posterior percentiles, \code{summary} contains a data.frame of the
#'   estimates, and \code{derivative} (if \code{derivative = TRUE})  log-derivative,
#'   \code{sampInt} median beta, \code{sampIntmean} mean beta,
#'   \code{sampInt975} and \code{sampInt025}97.5th and 2.5th posterior percentile beta
#' @export

BNPR_PS_ada<-function (data, lengthout = 100,  zeta1 = 0.01,
                    zeta2 = 0.01,
                    simplify = TRUE, derivative = FALSE, forward = TRUE)
{
  if (class(data) == "phylo") {
    phy <- summarize_phylo(data)
  }
  else if (all(c("coal_times", "samp_times", "n_sampled") %in%
               names(data))) {
    phy <- with(data, list(samp_times = samp_times, coal_times = coal_times,
                           n_sampled = n_sampled))
  }
  result <- infer_coal_samp_ada(samp_times = phy$samp_times, coal_times = phy$coal_times,
                                n_sampled = phy$n_sampled,  lengthout = lengthout,
                                zeta1 = zeta1, zeta2 = zeta2,simplify = simplify,
                                derivative = derivative)
  result$samp_times <- phy$samp_times
  result$n_sampled <- phy$n_sampled
  result$coal_times <- phy$coal_times
  result$effpop <- exp(-result$result$summary.random$time$`0.5quant`)
  result$effpopmean <- exp(-result$result$summary.random$time$mean)
  result$effpop975 <- exp(-result$result$summary.random$time$`0.025quant`)
  result$effpop025 <- exp(-result$result$summary.random$time$`0.975quant`)
  result$summary <- with(result$result$summary.random$time,
                         data.frame(time = ID, mean = exp(-mean), sd = sd * exp(-mean),
                                    quant0.025 = exp(-`0.975quant`), quant0.5 = exp(-`0.5quant`),
                                    quant0.975 = exp(-`0.025quant`)))
  if (derivative) {
    if (forward)
      ind <- c(1:(lengthout - 1), (lengthout - 1))
    else ind <- c(1, 1:(lengthout - 1))
    result$derivative <- with(result$result$summary.lincomb,
                              data.frame(time = result$x, mean = -mean[ind], sd = sd[ind],
                                         quant0.025 = -`0.975quant`[ind], quant0.5 = -`0.5quant`[ind],
                                         quant0.975 = -`0.025quant`[ind]))
  }

    result$sampInt <- exp(result$result$summary.random$time3$`0.5quant`)
    result$sampIntmean <- exp(result$result$summary.random$time3$mean)
    result$sampInt975 <- exp(result$result$summary.random$time3$`0.975quant`)
    result$sampInt025 <- exp(result$result$summary.random$time3$`0.025quant`)
    result$sampIntsummary <- with(result$result$summary.random$time3,
                                  data.frame(time = ID, mean = exp(mean), sd = sd * exp(mean),
                                             quant0.025 = exp(`0.025quant`), quant0.5 = exp(`0.5quant`),
                                             quant0.975 = exp(`0.975quant`)))

  return(result)
}


infer_coal_samp_ada <- function(samp_times, coal_times, n_sampled=NULL,
                                lengthout=100, zeta1=0.01, zeta2=0.01,
                                simplify = FALSE, events_only = FALSE,
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

  coal_data <- coal_stats_ada(grid = grid, samp_times = samp_times, n_sampled = n_sampled,
                              coal_times = coal_times)
  #practically above we are compunting for all sampling time, coalescent event and grid point the corresponding

  if (simplify)
    coal_data <- with(coal_data, condense_stats_ada(time=time, event=event, E=E)) #not sure what's happening here
  #not sure what's above

  hyper <- list(prec = list(param = c(zeta1,  zeta2)))

  if (events_only)
    samp_data <- samp_stats_new(grid = grid, samp_times = samp_times)
  else
    samp_data <- samp_stats_new(grid = grid, samp_times = samp_times,
                                n_sampled = n_sampled)

  data <- joint_stats_ada(coal_data = coal_data, samp_data = samp_data)


    data$beta0<-NULL
    data$time3=data$time2
    formula <- Y ~ -1 +
      f(time, model="rw1", hyper = hyper, constr = FALSE) +
      f(time2, w, copy="time") + f(time3, model="rw1",hyper = hyper, constr=FALSE)



  family <- c("poisson", "poisson")



  lc_many <- NULL


  mod <- INLA::inla(formula, family = family, data = data,
                    lincomb = lc_many, offset = data$E_log,
                    control.predictor = list(compute=TRUE),
                    control.inla = list(lincomb.derived.only=FALSE))

  return(list(result = mod, data = data, grid = grid, x = coal_data$time))
}


coal_stats_ada <- function(grid, samp_times, coal_times, n_sampled = NULL,
                           log_zero = -100)
{
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2 ##these are the midpoints of the grid.

  if (is.null(n_sampled))
    n_sampled <- rep(1, length(samp_times))
  args <- gen_INLA_args_ada(samp_times = samp_times, n_sampled = n_sampled,
                            coal_times = coal_times)

  coal_factor <- args$coal_factor
  s <- args$s
  event <- args$event

  grid_trimmed <- setdiff(x = grid, y = s)
  sorting <- sort(c(grid_trimmed, s), index.return=TRUE)
  sgrid <- sorting$x
  ordering <- sorting$ix

  time_index <- cut(x = sgrid[-1], breaks = grid, labels = FALSE) #this tells where grid points, coalescent times and sampling times are, in which grid box
  time <- field[time_index] #convert smapling times,coalescent time etc into mid grid points

  event_out <- c(rep(0, length(grid_trimmed)), event)[ordering] #determines what's what.

  Cfun <- stats::stepfun(x = s, y = c(0, coal_factor, 0), right = TRUE) #this is an interpolating function for the coalescent factors
  Cvec <- Cfun(sgrid[-1]) #this is the coalescent factor at each grid point, coalescent times and sampling times
  E <- diff(sgrid)*Cvec #delta of each interval I, times the corresponding grid point

  E_log = log(E)  #and log of it.
  E_log[E == 0] = log_zero

  return(data.frame(time = time, event = event_out[-1], E = E, E_log = E_log))
}


gen_INLA_args_ada <- function(samp_times, n_sampled, coal_times)
{
  if (sum(n_sampled) != length(coal_times) + 1)
    stop("Number sampled not equal to number of coalescent events + 1.")

  if (length(intersect(coal_times, samp_times)) > 0)
    warning("Coincident sampling event and coalescent event: results may be unpredictable.")

  l <- length(samp_times)
  m <- length(coal_times)
  sorting <- sort(c(samp_times, coal_times), index.return=TRUE)

  lineage_change <- c(n_sampled, rep(-1, m))[sorting$ix] ##this has the number of linearges
  lineages <- utils::head(cumsum(lineage_change), -1) # remove entry for the post-final-coalescent-event open interval
  #lineages above has the number of activ lineages
  coal_factor <- lineages*(lineages-1)/2

  event <- c(rep(0, l), rep(1, m))[sorting$ix] ### 0 is sampling event and 1 is a coaleascent event

  return(list(coal_factor=coal_factor, s=sorting$x, event=event, lineages=lineages))
}

condense_stats_ada <- function(time, event, E, log_zero = -100)
{
  result <- stats::aggregate(event ~ time, FUN = sum)
  result$E <- stats::aggregate(E ~ time, FUN = sum)$E

  E_log = log(result$E)
  E_log[result$E == 0] = log_zero
  result$E_log <- E_log

  return(result)
}

samp_stats_new <- function(grid, samp_times, n_sampled = NULL, trim_end = FALSE)
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

joint_stats_ada <- function(coal_data, samp_data)
{
  id_max_samp<-max(which(samp_data$count>=1))
  samp_data<-samp_data[1:id_max_samp,]
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



#Function taken from spmrf


varRef <- function(nn, kap, omg2, order=1){
  #nn is number of thetas, kap = 1/gam^2 (precision), omg2=omega^2 (variance of theta1)
  #returns vector of marginal variances
  nnv <- 1:nn
  if (order==1) out <- omg2 + (nnv-1)*(1/kap)
  if (order==2) out <- omg2 + (1/kap)*nnv*(nnv-1)*(2*nnv-1)/6
  if (order >= 3) out <- NULL
  return(out)
}


