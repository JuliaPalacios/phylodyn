############################################
######## PLOTTING TOOLS ####################
############################################

#' Plot the posterior of beta
#' it is the equivalent of plot_BNPR in phylodyn
#'
#' @param BNPR_out output of BNPR_PS_ada
#' @param traj function summarizing the true effective population size
#'   trajectory.
#' @param xlim numeric x-axis interval.
#' @param ylim numeric y-axis interval.
#' @param nbreaks integer number of bins for sampling heatmap.
#' @param lty line type for estimated trajectory.
#' @param lwd line width for estimated trajectory.
#' @param col color for estimated trajectory.
#' @param main character main plot title.
#' @param log character which axes to plot log-scale. Defaults to "y".
#' @param ylab character y-axis label.
#' @param xlab character x-axis label.
#' @param xmarline numeric if not using default x-axis labels, how far to put
#'   the labels from the axis.
#' @param axlabs character vector x-axis labels.
#' @param traj_lty,traj_lwd,traj_col line type, line width, and line color for
#'   the true trajectory.
#' @param newplot boolean whether to create a new plot or superimpose over a
#'   previously open plot.
#' @param credible_region logical whether to display pointwise credible region.
#' @param heatmaps boolean whether to display sampling and coalescent heatmaps.
#' @param heatmap_labels boolean whether to display labels on heatmaps.
#' @param heatmap_labels_side string which side of plot to display heatmaps.
#' @param heatmap_labels_cex numeric scaling factor for heatmap labels.
#' @param heatmap_width numeric how wide heatmaps should be.
#' @param yscale numeric scaling applied to all effective population
#'   calculations.
#' @param ... additional arguments to be passed onto plot().
#'
#' @export



plot_beta<-function (BNPR_out, traj = NULL, xlim = NULL, ylim = NULL, nbreaks = 40,
                   lty = 1, lwd = 2, col = "black", main = "", log = "y", ylab = "Sampling Intensity",
                   xlab = "Time", xmarline = 3, axlabs = NULL, traj_lty = 2,
                   traj_lwd = 2, traj_col = col, newplot = TRUE, credible_region = TRUE,
                   heatmaps = TRUE, heatmap_labels = TRUE, heatmap_labels_side = "right",
                   heatmap_labels_cex = 0.7, heatmap_width = 7, yscale = 1, max_samp,
                   ...)
{
  grid = BNPR_out$grid
  if (is.null(xlim)) {
    xlim = c(max(grid), min(grid))
  }
  mask = BNPR_out$x >= min(xlim) & BNPR_out$x <= max(xlim)
  t = BNPR_out$x[mask]
  #t=decimal_date(max_samp)-t
  xlim_dec=xlim
  #xlim=(c(decimal_date(max_samp)-xlim[1],decimal_date(max_samp)-xlim[2]))
  y = BNPR_out$sampInt[mask] * yscale
  yhi = BNPR_out$sampInt975[mask] * yscale
  ylo = BNPR_out$sampInt025[mask] * yscale
  if (newplot) {
    if (is.null(ylim)) {
      ymax = max(yhi)
      ymin = min(ylo)
    }
    else {
      ymin = min(ylim)
      ymax = max(ylim)
    }
    if (heatmaps) {
      yspan = ymax/ymin
      yextra = yspan^(1/10)
      ylim = c(ymin/(yextra^1.35), ymax)
    }
    else {
      ylim = c(ymin, ymax)
    }
    if (is.null(axlabs)) {
      graphics::plot(1, 1, type = "n", log = log, xlab = xlab,
                     ylab = ylab, main = main, xlim = xlim, ylim = ylim,
                     ...)
    }
    else {
      graphics::plot(1, 1, type = "n", log = log, xlab = "",
                     ylab = ylab, main = main, xlim = xlim, ylim = ylim,
                     xaxt = "n", ...)
      graphics::axis(1, at = axlabs$x, labels = axlabs$labs, cex.axis=1,
                     las = 1)
      graphics::mtext(text = xlab, side = 1, line = xmarline)
    }
  }
  if (credible_region) {
    shade_band(x = t, ylo = ylo, yhi = yhi, col = "lightgray")
  }
  if (!is.null(traj)) {
    graphics::lines(t, traj(t), lwd = traj_lwd, lty = traj_lty,
                    col = traj_col)
  }
  if (newplot) {
    if (heatmaps) {
      samps = rep(BNPR_out$samp_times, BNPR_out$n_sampled)
      #samps = decimal_date(max_samp)-samps[samps <= max(xlim_dec) & samps >= min(xlim_dec)]
      samps = samps[samps <= max(xlim_dec) & samps >= min(xlim_dec)]
      coals = BNPR_out$coal_times
      coals = coals[coals <= max(xlim) & coals >= min(xlim)]
      breaks = seq(min(xlim_dec), max(xlim_dec), length.out = nbreaks)
      h_samp = graphics::hist(samps, breaks = breaks, plot = FALSE)
      #h_coal = graphics::hist(coals, breaks = breaks, plot = FALSE)
      hist2heat(h_samp, y = ymin/yextra^0.5, wd = heatmap_width)
      #hist2heat(h_coal, y = ymin/yextra, wd = heatmap_width)
      if (heatmap_labels) {
        if (heatmap_labels_side == "left") {
          lab_x = max(xlim)
          lab_adj = 0
        }
        else if (heatmap_labels_side == "right") {
          lab_x = min(xlim)
          lab_adj = 1
        }
        else {
          warning("heatmap_labels_side not \"left\" or \"right\", defaulting to right")
          lab_x = min(xlim)
          lab_adj = 1
        }
        graphics::text(x = lab_x, y = ymin/(yextra^0.2)+0.0075,
                       labels = "Sampling events", adj = c(lab_adj,
                                                           0), cex = heatmap_labels_cex)
        #graphics::text(x = lab_x, y = ymin/(yextra^1.25),
        #   labels = "Coalescent events", adj = c(lab_adj,
        # 1), cex = heatmap_labels_cex)
      }
    }
  }
  graphics::lines(t, y, lwd = lwd, col = col, lty = lty)
}



############################################
######## SIMULATION SAMPLING TIMES #########
############################################


#' Sample sampling times from a Poisson process
#'
#' @param max.n number of samples requires
#' @param traj function describe the rate of Poisson proces
#' @param xlim time interval in which to do sampling
#' @param ... additional arguments to be passed to the traj function
#'
#' @export


sampsim_thin <- function(max.n, traj,xlim, ...)
{

  #curr = 1
  #active_lineages = n_sampled[curr]
  #time = samp_times[curr]
  time=0
  samp_times<-c()
  x<-seq(xlim[1],xlim[2],by=0.1)
  lambda<-max(traj(x))
  while (length(samp_times) <= max.n) {
    time=time-log(runif(1))/lambda
    if (runif(1)<=traj(time)/lambda) {samp_times<-c(samp_times,time)}
  }

  return(samp_times)
}




#skyline and rough estimates of the sampling intensity
skyLine_and_samp_int <- function(phy){
  # classic skyline for isochronous or heterochronous data.
  # phy is a list from summarize_phylo() with
  # coal_times, samp_times, and n_sampled
  # outputs vector of Ne estimates sorted from most recent.
  #The whole function is identical to the skyLine function in sprmf packaged except for the last few lines of code
  st <- phy$samp_times
  ns <- phy$n_sampled
  ct <- phy$coal_times
  sm <- data.frame(type='s', time=st, nadd=ns)
  cm <- data.frame(type='c', time=ct, nadd=-1)
  zm <- merge(sm, cm, all=T)
  zm <- zm[order(zm$time), ]
  zm$ncount <- cumsum(zm$nadd) #number of lineages at a the time
  zn <- nrow(zm) #number sampling events + number of coalescent event
  zuid <- numeric(zn)
  zuid[zm$type=='c'] <- 1:nrow(cm)
  for (j in (zn-1):1) {
    if (zuid[j]==0) zuid[j] <- zuid[j+1]
  } #it is the coalescent event time in which we are in
  zm$zuid <- zuid
  wk <- diff(zm$time)
  dm <- data.frame(uid=zm$zuid[-1], ncount=zm$ncount[-zn], wk=round(wk,8) ) #uid: which coalescent we are at
  #ncount: number of lineages, wk: length of the interval, thsub: binom coeff * lenght of interval
  dm$thsub <- dm$ncount*(dm$ncount-1)*dm$wk/2 #binomial coefficient
  adm <- aggregate(dm[,"thsub"], by=list(uid=dm$uid), sum) #sum contributions in the same coalescent interval
  zc1 <- c(0,cm$time[-nrow(cm)]) #coalescent times including 0, excluding TMRCA
  zc2 <- cm$time #coalescent times
  ctmid <- zc1 + (zc2-zc1)/2 #midpts between coalescent time
  #numerical correction skyle
  adm$x[adm$x<.000000001]<-.000000001
  #rough estimates of the sampling intensity
  last.samp<-min(which(zc2>max(st)))
  zm[zm[,"nadd"]==-1,"nadd"] <- 0
  samp.count <- as.vector(aggregate(zm[,"nadd"],by=list(zuid=zm$zuid),sum))
  beta <- samp.count[1:last.samp,2]/adm$x[1:last.samp]
  #beta[beta<.000000001]<-.000000001
  beta<-beta[beta>0]
  #adm$uid is the list of coalescent event, ctmid is the mid length point, adm$x is the estimate of the parameter
  return(list(k=adm$uid, mid_ctime=ctmid, theta=adm$x,beta=beta))
}



#' Precision parameters empirica proposals for the two fields
#'
#' @param phylo \code{phylo} object or list containing vectors of coalescent
#'   times \code{coal_times}, sampling times \code{samp_times}, and number
#'   sampled per sampling time \code{n_sampled}.
#' @param ncell number of grid points.
#' @param upBound
#' @param alpha quantile to consider
#' @param order.c order field on effective population size
#' @param order.s order field on beta
#'
#' @return Precision parameters,
#'  \code{zeta1}  and \code{zeta2}
#' @export


set_prec_parameters <- function(phylo, ncell, upBound=NULL, alpha=0.05, order.c=1,order.s=1){
  zsky <- skyLine_and_samp_int(phylo)
  vld <- var(log(zsky$theta))
  if (order.c==1)	cmr <- varRef(ncell, 1,  vld, order=1)
  if (order.c==2)	cmr <- varRef(ncell, 1,  vld, order=2)
  sref <- exp(mean(0.5*log(cmr)))
  if (is.null(upBound)) upBound <- sqrt(vld)
  zeta1 <- upBound/(sref*(tan((pi/2)*(1-alpha))))

  vlb <- var(log(zsky$beta))
  if (order.s==1)	cmr <- varRef(ncell, 1,  vlb, order=1)
  if (order.s==2)	cmr <- varRef(ncell, 1,  vlb, order=2)
  sref <- exp(mean(0.5*log(cmr)))
  if (is.null(upBound)) upBound <- sqrt(vlb)
  zeta2 <- upBound/(sref*(tan((pi/2)*(1-alpha))))

  return(zetas=list(zeta1=zeta1,zeta2=zeta2))
}




#' Extract posterior beta and effective population size from a stan object and puts them
#' into the phylodyn format in order to plot
#'
#' @param mfit stan fit object (output of adapref_sampling)
#' @param grid grid boundaries
#' @param midpts grid mid points
#' @param samp_times
#' @param coal_times
#' @param n_samples
#' @param alpha which BCI to keep
#'
#' @return BNPR style output
#' @export


extract_out <- function(mfit,grid,midpts,samp_times,coal_times,n_sampled, alpha=0.05){

  ##Assumes that you are extracting an output of a stan fit model
  tmp.th <- rstan::extract(mfit, "theta")[[1]]
  tmp.alp <- rstan::extract(mfit, "alpha")[[1]]

  plow <- alpha/2
  phigh <- 1 - alpha/2


  res<-list()
  res$grid <- grid
  res$x <- midpts
  res$n_sampled <- n_sampled
  res$samp_times <- samp_times
  res$coal_times <- coal_times

  res$effpop <- exp(apply(tmp.th, 2, median))
  res$effpopmean <- exp(apply(tmp.th, 2, mean))
  res$effpop025 <- exp(apply(tmp.th, 2, quantile, probs=plow))
  res$effpop975 <- exp(apply(tmp.th, 2, quantile, probs=phigh))
  res$summary<-data.frame(time=res$x,mean=res$effpopmean,quant0.025 =res$effpop025 , quant0.5 = res$effpop,
                          quant0.975 = res$effpop975)

  res$sampInt <- exp(apply(tmp.alp, 2, median))
  res$sampIntmean <- exp(apply(tmp.alp, 2, mean))
  res$sampInt025 <- exp(apply(tmp.alp, 2, quantile, probs=plow))
  res$sampInt975 <- exp(apply(tmp.alp, 2, quantile, probs=phigh))
  res$sampIntsummary<-data.frame(time=res$x,mean=c(res$sampIntmean,rep(NA,length(res$x)-length(res$sampIntmean))),quant0.025 =c(res$sampInt025,rep(NA,length(res$x)-length(res$sampIntmean))) ,
                                 quant0.5 = c(res$sampInt,rep(NA,length(res$x)-length(res$sampIntmean))), quant0.975 = c(res$sampInt975,rep(NA,length(res$x)-length(res$sampIntmean))))


  return(res)
}


#' Extract  effective population size from a stan object obtained with the sprmf package
#' and puts them into the phylodyn format in order to plot
#'
#' @param mfit stan fit object (output of adapref_sampling)
#' @param grid grid boundaries
#' @param midpts grid mid points
#' @param samp_times
#' @param coal_times
#' @param n_samples
#' @param alpha which BCI to keep
#'
#' @return BNPR style output
#' @export

extract_spmrf <- function(mfit,grid,midpts,samp_times,coal_times,n_sampled, alpha=0.05){

  ##Assumes that you are extracting an output of a stan fit model
  tmp.th <- rstan::extract(mfit, "theta")[[1]]

  plow <- alpha/2
  phigh <- 1 - alpha/2


  res<-list()
  res$grid <- grid
  res$x <- midpts
  res$n_sampled <- n_sampled
  res$samp_times <- samp_times
  res$coal_times <- coal_times

  res$effpop <- exp(apply(tmp.th, 2, median))
  res$effpopmean <- exp(apply(tmp.th, 2, mean))
  res$effpop025 <- exp(apply(tmp.th, 2, quantile, probs=plow))
  res$effpop975 <- exp(apply(tmp.th, 2, quantile, probs=phigh))
  res$summary<-data.frame(time=res$x,mean=res$effpopmean,quant0.025 =res$effpop025 , quant0.5 = res$effpop,
                          quant0.975 = res$effpop975)


  return(res)
}

