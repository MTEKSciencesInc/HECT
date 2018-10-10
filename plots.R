#' Probability of superiority plot
#'
#' @param trial An object of class 'trial'
#'
#' @return Probabilities of superiority for each arm evolving through the course of the trial
#' @export


psup_plot = function(trial, upper) {
  psup = trial$psup
  df = melt(psup)
  names(df) = c('treatment', 'interim_look', 'p.best')
  df$treatment = as.factor(df$treatment)
  df$interim_look = as.integer(df$interim_look) - 1
  p = ggplot(df, aes(x = interim_look, y = p.best, color = treatment)) +
    geom_line(aes(linetype = treatment),size = 1) + geom_hline(yintercept = upper, color = 'darkgrey') +
    scale_color_manual(values = cbPalette) +
    ylab('probability of superiority') + xlab('interim look') +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          strip.text.x = element_text(size = 8, face = "bold"))
  return(p)
}

#' Posterior density plot
#'
#' @param trial An object of class 'trial'
#' @param type a character representing the type of the parameter; "absolute", "OR" (odds ratio),
#' or "RR" (relative risk)
#' @return 5 snapshots of the posterior density plots for the effect sizes through teh trial
#' @export

post_plot = function(trial, type) {
  if (type == "absolute") {
    theta = trial$theta
    np = dim(theta)[3]
    select = floor(seq(1, np, length = 5))
    df = melt(theta)[,-1]
    names(df) = c('treatment', 'interim_look', 'effect')
    df0 = df[df$interim_look %in% select,]
    df0$interim_look = as.factor(df0$interim_look)
    #df1 = data.frame(theta0 = theta0, treatment = as.factor(sort(unique(df0$treatment))))
    p = ggplot(df0, aes(x = effect, fill = interim_look)) +
      geom_density(alpha = .5, color = 'grey') +
      #geom_vline(data = df1, aes(xintercept = theta0)) +
      facet_grid(treatment ~ ., labeller = label_both) +
      scale_fill_brewer() + xlim(-3, 3) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"),
            strip.text.y = element_text(size = 12, face = "bold"))
  }
  if (type == "count") {
    theta = exp(trial$theta)
    np = dim(theta)[3]
    select = floor(seq(1, np, length = 5))
    df = melt(theta)[,-1]
    names(df) = c('treatment', 'interim_look', 'intensity')
    df0 = df[df$interim_look %in% select,]
    df0$interim_look = as.factor(df0$interim_look)
    #df1 = data.frame(theta0 = theta0, treatment = as.factor(sort(unique(df0$treatment))))
    p = ggplot(df0, aes(x = intensity, fill = interim_look)) +
      geom_density(alpha = .5, color = 'grey') +
      #geom_vline(data = df1, aes(xintercept = theta0)) +
      facet_grid(treatment ~ ., labeller = label_both) +
      scale_fill_brewer() + xlim(0, 10) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"),
            strip.text.y = element_text(size = 12, face = "bold"))
  }
  if (type == "rate") {
    theta = exp(trial$theta)/(1 + exp(trial$theta))
    np = dim(theta)[3]
    select = floor(seq(1, np, length = 5))
    df = melt(theta)[,-1]
    names(df) = c('treatment', 'interim_look', 'rate')
    df0 = df[df$interim_look %in% select,]
    df0$interim_look = as.factor(df0$interim_look)
    #df1 = data.frame(theta0 = theta0, treatment = as.factor(sort(unique(df0$treatment))))
    p = ggplot(df0, aes(x = rate, fill = interim_look)) +
      geom_density(alpha = .5, color = 'grey') +
      #geom_vline(data = df1, aes(xintercept = theta0)) +
      facet_grid(treatment ~ ., labeller = label_both) +
      scale_fill_brewer() + xlim(0, 1) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"),
            strip.text.y = element_text(size = 12, face = "bold"))
  }
  return(p)
}

#' Data plot
#'
#' @param trial An object of class 'trial'
#'
#' @return A visual summary of the trial data
#' @export

data_plot = function(trial) {
  y = trial$y
  x = trial$x
  nt = nrow(x)
  treat = 1:nt
  df = data.frame(cbind(y, t(x)%*%treat, 1:length(y)))
  names(df) = c('response', 'treatment', 'patient')
  df$treatment = as.factor(df$treatment)
  if (all(y %in% c(0,1))) {
    df$response = as.factor(df$response)
    p = ggplot(df, aes(x = patient, y = response, color = response)) +
      geom_point(size = 3) + facet_grid(treatment ~ ., labeller = label_both) +
      scale_color_manual(values = cbPalette[2:3]) + 
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"),
            strip.text.y = element_text(size = 12, face = "bold"))
  } else {
    p = ggplot(df, aes(x = patient, y = response, color = response)) +
      geom_point(size = 3) + facet_grid(treatment ~ ., labeller = label_both) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"),
            strip.text.y = element_text(size = 12, face = "bold"))
  }
  return(p)
}


#' Estimates plot
#'
#' @param trial An object of class 'trial'
#'
#' @return Graph of estimates and 95% credible intervals
#' @export
#' 
estPlot = function(trial) {
  est = trial$est
  nt = nrow(est)
  treat = 1:nt
  est = data.frame(est, treatment = treat)
  est = est[rev(order(est$treatment)),]
  est$treatment = factor(est$treatment, levels = est$treatment[rev(order(est$treatment))])
  p = ggplot() + geom_linerange(data=est, mapping=aes(x=as.factor(treatment), ymin=low, ymax=up), size=1, color="darkblue") +
    geom_point(data=est, mapping=aes(x=as.factor(treatment), y=p.est), size=4, shape=21, fill="#33FFFF", color = 'darkgrey') +
    coord_flip() + ylab('estimate') + xlab('treatment') + theme(axis.text=element_text(size=12),
                                                              axis.title=element_text(size=14,face="bold"),
                                                              strip.text.y = element_text(size = 12, face = "bold"))
  return(p)
}

#' Design plot
#'
#' @param trial An object of class 'trial'
#'
#' @return A visual representation of the trial design
#' @export

designPlot = function(trial) {
  x = trial$x
  nt = nrow(x)
  nint = ncol(trial$psup) - 1
  intS = ncol(x)/nint
  start = c()
  end = c()
  for (i in 1:nt) {
    c0 = c()
    for (j in 1:nint) c0[j] = 1 %in% x[i, ((j-1)*intS +1):(j*intS)]
    start[i] = min(which(c0 == TRUE)) - 1
    end[i] = max(which(c0 == TRUE))
  }
  df = data.frame(arm = paste('Treatment', 1:nt), start = as.integer(start), end = as.integer(end))
  #df$arm = factor(df$arm, levels = levels(df$arm)[rev(order(df$arm))])
  p = ggplot(df, aes(colour=arm, x=start, y = arm, label = arm)) + 
    geom_segment(aes(x=start, xend=end, y=arm, yend=arm), size=15, alpha = .75) +
    xlab("Interim look") + geom_text(aes(x = start + nint/10), fontface = "bold", color = 'black') + 
    theme_classic() + 
    scale_color_manual(values = cbPalette) + 
    theme(axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x=element_text(size=12),
          axis.title.x=element_text(size=14,face="bold"),
          legend.position="none")
  return(p)
}


#' HECTvsBRCT plot
#'
#' @param a0 alpha for HECT
#' @param a1 alpha for BRCT
#' @param p0 power for HECT
#' @param p1 power for BRCT
#' @param c0 cost for HECT
#' @param c1 cost for BRCT
#'
#' @return Comparison plots
#' @export

HECTvsBRCTPlot = function(a0, a1, p0, p1, c0, c1, compCon) {
  dfc = data.frame(cost = c(c0, c1), design = c('HECT', 'RCT'))
  if (compCon == T) {
    nt0 = length(a0)
    df = data.frame(alpha = c(a0, a1), power = c(p0, p1), 
                    design = c(rep('HECT', nt0) , rep('RCT', nt0)), 
                    arm = rep(paste('Treatment', 2:(nt0+1)), 2))
    p1 = ggplot(df, aes(fill=design, x=design, y = alpha)) + geom_bar(stat = 'identity', alpha = .75) +
      theme_classic() + facet_grid(arm ~.) + 
      scale_fill_manual(values=cbPalette) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"),
            strip.text = element_text(face="bold", size=12),
            legend.position="none")
    p2 = ggplot(df, aes(fill=design, x=design, y = power)) + geom_bar(stat = 'identity', alpha = .75) +
      theme_classic() + facet_grid(arm ~.) + 
      scale_fill_manual(values=cbPalette[3:4]) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"),
            strip.text = element_text(face="bold", size=12),
            legend.position="none")
  } else {
    df = data.frame(alpha = c(a0, a1), power = c(p0, p1), design = c('HECT', 'RCT'))
    p1 = ggplot(df, aes(fill=design, x=design, y = alpha)) + geom_bar(stat = 'identity', alpha = .75) +
      theme_classic() + ylab('type I error rate') +
      scale_fill_manual(values=cbPalette) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"),
            strip.text = element_text(face="bold", size=12),
            legend.position="none")
    p2 = ggplot(df, aes(fill=design, x=design, y = power)) + geom_bar(stat = 'identity', alpha = .75) +
      theme_classic() + 
      scale_fill_manual(values=cbPalette[3:4]) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"),
            strip.text = element_text(face="bold", size=12),
            legend.position="none")
  }
  p3 = ggplot(dfc, aes(fill=design, x=design, y = cost)) + geom_bar(stat = 'identity', alpha = .75) +
    theme_classic() + 
    scale_fill_manual(values=cbPalette[5:6]) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          strip.text = element_text(face="bold", size=12),
          legend.position="none")
  return(grid.arrange(p1, p2, p3, nrow = 1, ncol = 3))
}
  


