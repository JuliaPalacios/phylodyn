###############################
###Trial analysis function ####
###############################
###Author: Lorenzo Cappello


#correct distance matrix if we are in the heterochrnous case. 
correct_distance_het2<-function(n1,samplingdates,dm1){
  c1=matrix(rep(samplingdates,n1),n1,n1,byrow=FALSE)
  c2=matrix(rep(samplingdates,n1),n1,n1,byrow=TRUE)
  correction=c1+c2-2*diag(samplingdates)
  CorrDistance=dm1+correction
  return(CorrDistance)
} 

#NEW DOUBLE CHECK
upgma_tree_correction<-function(CorrDistance,n_sample,samp_times,name_samp){
  #Create UPGMA tree
  treeUPGMA1<-upgma(CorrDistance)
  tmrca<-sum(coalescent.intervals(treeUPGMA1)$interval.length)
  if(length(samp_times)>1){#heterochronous samples
    treeUPGMA1_het<-trimBranches2(treeUPGMA1,n_sample=dim(name_samp)[1],"plus",name_samp)
    treeUPGMA1_het$edge.length[treeUPGMA1_het$edge.length<.000000000001]<-max(min(treeUPGMA1_het$edge.length/10),.000000000001)
    #correct for possibly negative edges
    while(sum(treeUPGMA1_het$edge.length<=0*1)>=1){
      MinEdge=min(abs(treeUPGMA1_het$edge.length))/10 #MinEdge length in the UPGMA tree: we assign it to the negative edge
      idn=min(which(treeUPGMA1_het$edge.length<=0))
      parNode=treeUPGMA1_het$edge[idn,1]
      idPar=which(treeUPGMA1_het$edge[,2]==parNode)
      idChild=which(treeUPGMA1_het$edge[,1]==parNode) #Look for the child node complementary to the negative one
      chilNode=treeUPGMA1_het$edge[idChild[idChild!=idn],2]
      fill=abs(treeUPGMA1_het$edge.length[idn])+MinEdge #Edge length to compensate
      treeUPGMA1_het$edge.length[idn]=MinEdge #correct the negative edge
      treeUPGMA1_het$edge.length[idChild[idChild!=idn]]=treeUPGMA1_het$edge.length[idChild[idChild!=idn]]+fill #Compensate 
      treeUPGMA1_het$edge.length[idPar]=treeUPGMA1_het$edge.length[idPar]-fill #Compensate 
    }
  } else {treeUPGMA1_het=treeUPGMA1}
  treeUPGMA1_het$edge.length[treeUPGMA1_het$edge.length<.000000000001]<-max(min(treeUPGMA1_het$edge.length)/10,.000000000001)
  treeUPGMA1_het$tmrca<-tmrca
  return(tree=treeUPGMA1_het)
}



# function to trim branches from UPGMA tree
trimBranches2 <- function(tree, n_sample, sign,name_samp) {
  #Define the corr vector
  #Define samp_times of 
  idx<-c()
  for (i in 1:length(tree$tip.label)){
    idx<-c(idx,which(name_samp[,2]==tree$tip.label[i]))
  }
  corr<-as.numeric(name_samp[idx,1])
  if (sign=="minus"){
    corr<-corr*-1
  }
  idLeafEDGE<-which(tree$edge[,2]<=n_sample)
  tree$edge.length[idLeafEDGE]<-tree$edge.length[idLeafEDGE]-corr[tree$edge[idLeafEDGE,2]]
  return(tree)
}  


serial_upgma<-function(fastafile,mu, samp_times, name_samp, model="JC69"){
  dm <- as.matrix(dist.ml(fastafile,model))
  dm1=as.matrix(dm/mu)
  corr_distance<-correct_distance_het2(dim(name_samp)[1],samp_times,dm1)
  tree<-upgma_tree_correction(corr_distance,n_sample=dim(name_samp)[1],samp_times,name_samp)
  return(tree)
}

serial_upgma_inputDist<-function(distList,mu, samp_times, name_samp){
  dm <- as.matrix(distList$distGen)
  dm1=as.matrix(dm/mu)
  corr_distance<-correct_distance_het2(dim(name_samp)[1],samp_times,dm1)
  tree<-upgma_tree_correction(corr_distance,n_sample=dim(name_samp)[1],samp_times,name_samp)
  return(tree)
}


mu_linear_reg<-function(alig){
  n = length(alig)
  ids = names(alig)
  dates_samp<-ymd(substr(ids, start = nchar(ids) - 9, stop = nchar(ids)))
  ##A rough estimate of mutation rate
  x<-dates_samp-dates_samp[length(dates_samp)]
  hamming<-as.matrix(dist.hamming(as.phyDat(alig),ratio=T))
  reg<-lm(hamming[nrow(hamming),-nrow(hamming)]~-1+x[-nrow(hamming)])
  #Convert data to 0s and 1s - Assumes the reference sequence is the last row
  data.matrix<-as.character(as.matrix(alig))
  mu<-reg$coefficients[[1]]*365
  return(mu)
}
#this is the original function as a backup
# mu_linear_reg_inputDist<-function(distList){
#   n = distList$n
#   ids = distList$seq_names
#   dates_samp<-ymd(substr(ids, start = nchar(ids) - 9, stop = nchar(ids)))
#   ##A rough estimate of mutation rate
#   x<-dates_samp-dates_samp[length(dates_samp)]
#   hamming<-as.matrix(distList$hamming)
#   reg<-lm(hamming[nrow(hamming),-nrow(hamming)]~-1+x[-nrow(hamming)])
#   #Convert data to 0s and 1s - Assumes the reference sequence is the last row
#   mu<-reg$coefficients[[1]]*365
#   return(mu)
# }

mu_linear_reg_inputDist<-function(distList){
  n = distList$n
  ids = distList$seq_names
  dates_samp<-ymd(substr(ids, start = nchar(ids) - 9, stop = nchar(ids)))
  anc_date<-ymd("2019-12-31")
  ##A rough estimate of mutation rate
 # x<-dates_samp-dates_samp[length(dates_samp)]
  x<-dates_samp-anc_date
  hamming<-as.matrix(distList$hamming)
  reg<-lm(hamming~-1+x)
  #Convert data to 0s and 1s - Assumes the reference sequence is the last row
  mu<-reg$coefficients[[1]]*365
  return(mu)
}






axis_label<-function(bnp,lastdate,byy=5/365){
  bnp1<-bnp
  bnp1$grid <- lastdate - bnp1$grid
  bnp1$x <- lastdate - bnp1$x
  bnp1$samp_times <- lastdate - bnp1$samp_times
  bnp1$coal_times <- lastdate - bnp1$coal_times
  labs<-format(date_decimal(bnp1$x), "%d-%m-%Y")
  mindate<-min(bnp1$x)
  maxdate<-max(bnp1$x)
  axlabs = list(x =seq(mindate,maxdate,by=byy), labs = format(date_decimal(seq(mindate,maxdate,by=byy)), "%b-%d"),cexlab=.1)
  axlabs$x<-lastdate-axlabs$x
  return(axlabs)
}


plot_BNPR2<-function (BNPR_out, traj = NULL, xlim = NULL, ylim = NULL, nbreaks = 40,
                      lty = 1, lwd = 2, col = "black", main = "", log = "y", ylab = "Effective Population Size",
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
  y = BNPR_out$effpop[mask] * yscale
  yhi = BNPR_out$effpop975[mask] * yscale
  ylo = BNPR_out$effpop025[mask] * yscale
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


##the functions below do not really need to be added if we add all the above to phylodyn. 

shade_band = function(x, ylo, yhi, xlim=NULL, col="gray")
{
  if (is.null(xlim))
    xlim = c(0, Inf)
  mask = x >= min(xlim) & x <= max(xlim)
  
  x = x[mask]
  ylo = ylo[mask]
  yhi = yhi[mask]
  
  graphics::polygon(c(x, rev(x)), c(yhi, rev(ylo)), col=col, border=NA)
}

float2gray = function(f)
{
  return(sprintf("#%02x%02x%02x", floor((1-f) * 255), floor((1-f) * 255), floor((1-f) * 255)))
}

#' Plot heatmap of a histogram
#' 
#' @param hist \code{histogram} object to be displayed.
#'
#' @param y numeric y-coordinate to display heatmap.
#' @param wd numeric width of heatmap in y-units.
#'
#' @export
hist2heat = function(hist, y, wd)
{
  breaks = hist$breaks
  counts = hist$counts
  upper  = max(counts)
  n = length(counts)
  cols = float2gray(counts / upper)
  graphics::segments(x0 = breaks[1:n], y0=y, x1 = breaks[2:(n+1)], y1=y, lwd=wd, col=cols, lend=1)
}

