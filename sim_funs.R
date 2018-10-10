#' @importFrom stats pnorm
#' @importFrom stats rbinom
#' @importFrom stats rmultinom
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @importFrom stats var
#' @importFrom stats dnorm
#' @importFrom stats dpois
#' @importFrom stats rpois
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom stats filter
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 scale_color_brewer
#' @importFrom ggplot2 scale_fill_brewer
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom reshape2 melt
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom abind abind
NULL

#' Closed form Gaussian Posterior
#'
#' Returns updated Gaussian posterior parameters for the effect sizes
#'
#' @param yb vector of new responses, length equals number of patients in the most recent batch
#' @param xb design matrix for the most recent batch
#' @param prior matrix of Gaussian parameters up to the most recent batch, number of rows is the
#' number of arms, first column contains the means and second column the sd's.
#' @param ysd observation variance.
#' @return matrix of means and sd's for the Gaussian posterior
#'
#' @export

post_Gauss = function(yb, xb, prior, ysd) {
  if (length(yb) > 0) {
    nb = rowSums(xb)
    mp = prior[,1]
    sdp = prior[,2]
    ybar = (xb%*%yb)/rowSums(xb)
    ybar[is.na(ybar)] = 0
    var = 1 / (nb/ysd^2 + 1/sdp^2)
    post_mean = var * (mp/sdp^2 + nb*ybar/ysd^2)
    out = cbind(post_mean, post.sd = sqrt(var))
  }
  else out = cbind(rep(0, nrow(prior)), rep(0, nrow(prior)))
  return(out)
}

#' Closed form Beta Posterior
#'
#' Returns updated Beta posterior parameters for the rate of success
#'
#' @param yb vector of new responses, length equals number of patients in the most recent batch
#' @param xb design matrix for the most recent batch
#' @param prior matrix of Beta parameters up to the most recent batch, number of rows is the
#' number of arms.
#' @return matrix of parameters of the Beta posterior
#'
#' @export

post_beta = function(yb, xb, prior) {
  if (length(yb) > 0) {
    nb = rowSums(xb)
    alphap = prior[,1]
    betap = prior[,2]
    y = (xb%*%yb)
    #ybar[is.na(ybar)] = 0
    alpha = y + alphap
    beta = nb - y + betap
    out = cbind(alpha, beta)
  }
  else out = cbind(rep(0, nrow(prior)), rep(0, nrow(prior)))
  return(out)
}

#' Closed form Gamma Posterior
#'
#' Returns updated Beta posterior parameters for the rate of success
#'
#' @param yb vector of new responses, length equals number of patients in the most recent batch
#' @param xb design matrix for the most recent batch
#' @param prior matrix of Gamma parameters up to the most recent batch, number of rows is the
#' number of arms.
#' @return matrix of parameters of the Gamma posterior
#'
#' @export

post_gamma = function(yb, xb, prior) {
  if (length(yb) > 0) {
    nb = rowSums(xb)
    alphap = prior[,1]
    betap = prior[,2]
    y = (xb%*%yb)
    #ybar[is.na(ybar)] = 0
    alpha = y + alphap
    beta = nb + betap
    out = cbind(alpha, beta)
  }
  else out = cbind(rep(0, nrow(prior)), rep(0, nrow(prior)))
  return(out)
}

#' Bernoulli likelihood
#'
#' Returns the log-likelihood for a vector of parameter values, rate responses and design matrix
#'
#' @param theta vector of effect sizes, length equals number of arms
#' @param y vector of count responses, length equals number of patients
#' @param x design matrix
#' @param ... covariate related arguments, \code{beta}: vector of coefficients, length equals number of covariates,
#' \code{z}: covariate matrix
#' @return log-likelihood value
#'
#' @export

llik_bernoulli = function(theta, y, x, ...) {
  more = list(...)
  if (!is.null(more$beta)){
    beta = more$beta
    z = more$beta
    add = t(z)%*%beta
  } else add = 0
  if (!all(y %in% c(0,1))) stop('non-binary response, y takes 0 or 1 values')
  c = c(t(x)%*%theta + add)
  out = c(mapply(function(y, c) log(y * (1 - pnorm(0, c, 1)) +
                                      (1 - y) * pnorm(0, c, 1)), y, c))
  out = sum(out)
  if (out == -Inf) out = -1e6
  return(out)
}

#' Poisson likelihood
#'
#' Returns the log-likelihood for a vector of parameter values, count responses and design matrix
#'
#' @param theta vector of effect sizes, length equals number of arms
#' @param y vector of count responses, length equals number of patients
#' @param x design matrix
#' @param ... covariate related arguments, \code{beta}: vector of coefficients, length equals number of covariates,
#' \code{z}: covariate matrix
#' @return log-likelihood value
#'
#' @export

llik_poisson = function(theta, y, x, ...) {
  more = list(...)
  if (!is.null(more$beta)){
    beta = more$beta
    z = more$beta
    add = t(z)%*%beta
  } else add = 0
  if (!all(y >= 0)) stop('negative response, y takes non-negative')
  #if (all(y %in% c(0,1))) warning('response values are 0 or 1; response.type = rate may be a better option')
  c = exp(c(t(x)%*%theta + add))
  out = c(mapply(function(y, c) dpois(y, c, log = T), y, c))
  out = sum(out)
  if (out == -Inf) out = -1e6
  return(out)
}

#' Prior
#'
#' Returns log-prior for a vector of parameter values - this is the default normal prior function in case the user does not want to
#' bother to decide on a prior distribution
#'
#' @param par vector of parameter values
#' @param ... the normal prior hyper-parameters (optional), can be specified as vectors of the same size as \code{par} or
#' scalars. The default is \code{pmean = 0} and \code{psd = 10}.
#' @return log-prior value
#'
#' @export

lp_default = function(par, ...) {
  args = list(...)
  pmean = if(is.null(args$pmean)) 0 else args$pmean
  psd = if(is.null(args$psd)) 10 else args$psd
  out = if(is.null(par)) 0 else sum(mapply(dnorm, x = par, mean = pmean, sd = psd, MoreArgs = list(log = T)))
  return(out)
}

#' Posterior
#'
#' Returns log-posterior for a vector of parameter values, count responses and design matrix
#'
#' @param theta vector of effect sizes, length equals number of arms
#' @param y vector of responses, length equals number of patients
#' @param x design matrix
#' @param response.type response type, either 'rate' or 'count'
#' @param lp_theta a function returning log-prior for effect sizes; the default is lp_default which is a diffuse normal prior,
#' but can be replaced with a user-specified functin
#' @param ... covariate related arguments, \code{beta}: vector of coefficients, length equals number of covariates,
#' \code{z}: covariate matrix, \code{lp_beta}: a function returning log-prior for coefficients; the default is lp_default which is a diffuse normal prior,
#' but can be replaced with a user-specified functin
#' @return log-posterior value
#'
#' @export

lpost = function(theta, y, x, response.type, lp_theta = lp_default, ...) {
  more = list(...)
  if (!is.null(more$beta)) {
    beta = more$beta
    z = more$z
    lp_beta = more$lp_beta
    bprior = lp_beta(beta)
  } else bprior = 0
  if (response.type == 'rate') ll = llik_bernoulli(theta, y, x, ...)
  if (response.type == 'count') ll = llik_poisson(theta, y, x, ...)
  lp = lp_theta(theta) + bprior
  return(ll + lp)
}

#' MCMC step
#'
#' This function perfomrs the nested MCMC steps required in SMC
#'
#' @param theta vector of effect sizes, length equals number of arms
#' @param y vector of responses, length equals number of patients
#' @param x design matrix
#' @param qt proposal standard deviation
#' @param response.type response type, either 'rate' or 'count'
#' @param ng number of MCMC steps for each particle. It can be as small as 1, but need to be larger if batch updates are made
#' @param ... covariate related arguments, \code{beta}: vector of coefficients, length equals number of covariates,
#' \code{z}: covariate matrix, \code{lp_beta}: a function returning log-prior for coefficients; the default is lp_default which is a diffuse normal prior,
#' but can be replaced with a user-specified functin, \code{qb}: proposal variance for \code{beta}
#' @return vector of updated parameters
#'
#' @export

Gibbs = function(theta, y, x, qt, response.type, ng, ...) {
  more = list(...)
  nt = length(qt)
  new = theta + qt*rnorm(nt)
  for (j  in 1:ng) {
    for (i in 1:nt) {
      new_th = theta
      new_th[i] = new[i]
      lp = lpost(theta, y, x, response.type = response.type, ...)
      lpnew = lpost(new_th, y, x, response.type = response.type, ...)
      p = min(1, exp(lpnew - lp))
      if (runif(1) <= p) theta = new_th
    }
    if (!is.null(more$beta)) {
      beta = more$beta
      z = more$z
      qb = more$qb
      if (!is.null(more$lp_beta)) lp_beta = more$lp_beta else lp_beta = lp_default
      nb = length(qb)
      newb = beta + qb*rnorm(nb)
      for (i in 1:nb) {
        new_b = beta
        new_b[i] = newb[i]
        lp = lpost(theta, y, x, response.type = response.type, beta = beta, z = z, lp_beta = lp_beta)
        lpnew = lpost(theta, y, x, response.type = response.type, beta = newb, z = z, lp_beta = lp_beta)
        p = min(1, exp(lpnew - lp))
        if (runif(1) <= p) beta = new_b
      }
    } else beta = NULL
  }
  out = list(theta = theta, beta = beta)
  return(out)
}

#' SMC update
#'
#' This function updates a posterior sample of effect sizes according to newly observed data
#'
#' @param theta matrix of posterior samples for effect sizes, columns: arms, rows: samples
#' @param y vector of recent responses
#' @param x design matrix for the recent observations
#' @param N number of particles
#' @param response.type response type, either 'rate' or 'count'
#' @param ... covariate related arguments, \code{beta}: vector of coefficients, length equals number of covariates,
#' \code{z}: covariate matrix
#' @return updated posterior samples
#'
#' @export

SMC_step = function(theta, y, x, N, response.type = 'rate', ...) {
  more = list(...)
  smA = split(theta, row(theta))
  if (!is.null(more$beta)) {
    beta = more$beta
    z = more$z
    smB = split(beta, row(beta))
    if (response.type == 'rate') ll = mapply(llik_bernoulli, smA, beta = smB, MoreArgs = list(x = x, y = y, z = z))
    if (response.type == 'count') ll = mapply(llik_poisson, smA, beta = smB,  MoreArgs = list(x = x, y = y, z = z))
  }
  else {
    if (response.type == 'rate') ll = apply(theta, 1, llik_bernoulli, x = x, y = y)
    if (response.type == 'count') ll = apply(theta, 1, llik_poisson, x = x, y = y)
    beta = NULL
  }
  w = exp(ll)/sum(exp(ll))
  index = sample(1:N, N, prob = w, replace = T)
  theta = theta[index,]
  out = list(theta = theta, beta = beta)
  return(out)
}

#' Probabilities of superiority
#'
#' This function estimates the probability of seperiority for each treatment arm from a posterior sample of effect sizes
#'
#' @param theta matrix of posterior samples for effect sizes, columns: arms, rows: samples
#' @return vector of probabilities of superiority
#'
#' @export

sup_check = function(theta) {
  n = length(theta)
  check = c()
  for (i in 1:n) check[i] = (theta[i] == max(theta))
  return(check)
}

#' futility check
#'
#' This function is used to estimate the probability of futility of each treatment arm compare to control
#'
#' @param theta matrix of posterior samples for effect sizes, columns: arms, rows: samples
#' @param MID minimally important difference
#' @return vector of probabilities of superiority
#'
#' @export

con_fut_check = function(theta, MID) {
  n = length(theta)
  check = c()
  for (i in 2:n) check = c(check, (theta[i] - theta[1])<MID)
  return(check)
}

#' Superiority to control check
#'
#' This function estimates the probability of seperiority for each treatment arm from a posterior sample of effect sizes
#'
#' @param theta matrix of posterior samples for effect sizes, columns: arms, rows: samples
#' @return vector of probabilities of superiority
#'
#' @export

con_sup_check = function(theta) {
  n = length(theta)
  check = c()
  for (i in 2:n) check = c(check, theta[i] > theta[1])
  return(check)
}

#' RAR simulator
#'
#' This function simulates a RAR trial
#'
#' @param nt Number of treatment arms
#' @param theta0 vector of effect sizes
#' @param nb batch size for each update
#' @param maxN maximum number of patients that the trial can run for; the trial is stopped when \code{maxN} is achieved
#' @param N Number of particles
#' @param upper upper threshold for probability of superiority; the trial is stopped when one of the trial arms reaches
#' this value
#' @param lower lower threshold for probability of superiority; an arm is stopped when its corresponding probability of
#' superiority falls below this value
#' @param burn the size of pre-adaptation sample; the default is 10 times number of treatments
#' @param response.type response type, either 'rate', 'count' or 'absolute'
#' @param conjugate_prior logical argument, default is True
#' @param padhere vector of adherence rates for each treatment
#' @param compCon logical, compare all arms with control (first arm)
#' @param adapt logical, adapt randomization probabilities, used only when \code{compCon = F}
#' @param platf logical, use a platform design, i.e., add the last arm only when at least one arm is dropped
#' @param MID minimally important difference for futility test of paiwise comparisons with control; used only when \code{compCon = T}
#' @return updated posterior samples
#'
#' @export

RAR_sim = function(nt, theta0, good.out = T, nb = 1, maxN = 500, N = 1000, upper = 0.975, uppfut = 0.95, lower = .01,
                   burn = 10*nt, response.type, conjugate_prior = T, padhere = rep(1,nt), adapt = T,
                   platf = T, compCon = F, MID = 0) {
  if (good.out == F) {
    if (response.type == 'absolute') theta0 = - theta0
    if (response.type == 'rate') theta0 = 1 - theta0
  }
  ng = nb
  addarmlater = rep(0, nt)
  if (platf == T) addarmlater[nt] = Inf
  j = 0
  x = array(0, dim = c(nt, 1))
  x0 = array(0, dim = c(nt, 1))
  y = NULL
  theta = array(rnorm(N*nt, 0, 10), dim = c(N, nt, 1))
  check = array(0, dim = c(N, nt))
  fcheck = array(0, dim = c(N, nt))
  if (conjugate_prior != T | response.type == "absolute") post0 = cbind(rep(0, nt), rep(10, nt))
  if (conjugate_prior == T & response.type == "rate") {
    post0 = cbind(rep(1, nt), rep(1, nt))
    p_new = apply(post0, 1, function(x) rbeta(N, x[1], x[2]))
    theta = array(- log((1 - p_new) / p_new), dim = c(N, nt, 1))
  }
  if (conjugate_prior == T & response.type == "count") {
    post0 = cbind(rep(.1, nt), rep(.001, nt))
    lambda_new = apply(post0, 1, function(x) rgamma(N, x[1], x[2]))
    theta = array(log(lambda_new), dim = c(N, nt, 1))
  }
  psup = array(rep(1/nt, nt), dim = c(nt, 1))
  nt0 = length(which(addarmlater == 0))
  psup_out = array(rep(1/nt0, nt), dim = c(nt, 1))
  psup_out[which(addarmlater>0),1] = 0
  prand = rep(1/nt, nt)
  prand[which(addarmlater>0)] = 0
  fut.stop = NULL
  repeat {
    j = j + 1
    # if (is.function(updateProgress)) {
    #   text <- paste0("updating results based on batch/patient:", j)
    #   updateProgress(detail = text)
    # }
    if (j*nb>maxN) nb = maxN - (j-1)*nb
    xb = rmultinom(nb, 1, prob = prand)
    xb0 = xb
    for (k in 1:nb) {
      xb0[which(xb0[,k] == 1), k] = rbinom(1, 1, padhere[which(xb0[,k] == 1)])
    }
    if (response.type == 'rate') {
      yb = apply(t(xb0)%*%theta0, 1, function(z) rbinom(1, 1, prob = z))
    }
    if (response.type == 'count') {
      yb = apply(t(xb0)%*%theta0, 1, function(z) rpois(1, exp(z)))
    }
    if (response.type == 'absolute') {
      yb = apply(t(xb0)%*%theta0, 1, function(z) rnorm(1, z))
    }
    x = abind(x, xb, along = 2)
    x0 = abind(x0, xb0, along = 2)
    y = c(y, yb)
    if ((response.type == 'rate' | response.type == 'count') & conjugate_prior != T) {
      smc_out = SMC_step(theta[,,j], x = xb0, y = yb, N = N, response.type)
      theta_new = smc_out$theta
      qt = apply(theta_new, 2, sd)/3
      Gibbs_out = apply(theta_new, 1, Gibbs, y = y, x = x[,-1], qt = qt,
                        response.type = response.type, ng)
      theta_new = t(sapply(Gibbs_out, function(l) return(l$theta)))
    }
    if (response.type == 'rate' & conjugate_prior == T) {
      post0 = post_beta(yb, xb, post0)
      p_new = apply(post0, 1, function(x) rbeta(N, x[1], x[2]))
      theta_new = - log((1 - p_new) / p_new)
    }
    if (response.type == 'count' & conjugate_prior == T) {
      post0 = post_gamma(yb, xb, post0)
      lambda_new = apply(post0, 1, function(x) rgamma(N, x[1], x[2]))
      theta_new = log(lambda_new)
    }
    if (response.type == 'absolute') {
      ysd = ifelse(length(y) == 1, 1, sd(y))
      post0 = post_Gauss(yb, xb, post0, ysd)
      theta_new = apply(post0, 1, function(x) rnorm(N, x[1], x[2]))
    }
    #theta_new[,which(addarmlater>j)] = -Inf
    theta = abind(theta, theta_new, along = 3)
    
    if (compCon == F) {
      check[,which(prand>0)] = t(apply(theta[,which(prand>0),j+1], 1, sup_check))
      check[,which(prand==0)] = 0
      psup_out = abind(psup_out, apply(check, 2, mean), along = 2)
      psup_out[which(addarmlater>j), j+1] = 0
      if (length(y) < burn) {
        psup = abind(psup, psup[,j], along = 2)
      } else psup = abind(psup, apply(check, 2, mean), along = 2)
      ll = NULL
      if (length(which(psup[,j+1] < lower))>0) {
        ll = which(psup[,j+1] < lower)
        }
      if (length(ll[addarmlater[ll]<j]) != 0) {
        fut.stop = c(fut.stop, ll[addarmlater[ll]<j])
      }
      if (!is.null(fut.stop)) {
        psup[fut.stop,j+1] = 0
        if (platf == T) addarmlater[nt] = j
      }
      ntj = length(which(addarmlater<=j & psup[,j+1]>0)) 
      if (sum(addarmlater==j)>0) ntj = ntj + sum(addarmlater==j)
      prand = rep(0, nt)
      if (adapt == T & length(y)>=burn)  {
        prand[which(addarmlater<=j & psup[,j+1]>0)] = ((ntj - 1)/ntj)*sqrt(psup[which(addarmlater<=j & psup[,j+1]>0),j+1])
      } else prand[which(addarmlater<=j & psup[,j+1]>0)] = 1/ntj
      prand[which(addarmlater==j)] = 1/ntj
      
      
    } else {
      mat = apply(theta[,which(prand>0),j+1], 1, con_sup_check)
      if (is.null(dim(mat))) mat = matrix(mat, N, length(which(prand>0)) - 1) else mat = t(mat)
      check[,which(prand>0)] = cbind(rep(0, N), mat)
      check[,which(prand==0)] = 0
      psup_out = abind(psup_out, apply(check, 2, mean), along = 2)
      psup_out[which(addarmlater>j), j+1] = 0
      if (length(y) < burn) {
        psup = abind(psup, psup[,j], along = 2)
      } else psup = abind(psup, apply(check, 2, mean), along = 2)
      if (response.type == 'absolute') fmat = apply(theta[,which(prand>0),j+1], 1, con_fut_check, MID = MID)
      if (response.type == 'rate') fmat = apply(p_new[,which(prand>0)], 1, con_fut_check, MID = MID)
      if (is.null(dim(fmat))) fmat = matrix(fmat, N, length(which(prand>0)) - 1) else fmat = t(fmat)
      fcheck[,which(prand>0)] = cbind(rep(0, N), fmat)
      fcheck[,which(prand==0)] = 0
      pfut = apply(fcheck, 2, mean)
      pfut[which(addarmlater>j)] = 0
      ff = NULL
      if (length(which(pfut > uppfut))>0) ff = which(pfut > uppfut)
      if (!is.null(ff)) {
        fut.stop = c(fut.stop, ff[addarmlater[ff]<j])
        
      }
      if (!is.null(fut.stop)) {
        psup[fut.stop,j+1] = 0
        if (platf == T) addarmlater[nt] = j
      }
      psup_out[,1] = c(0,rep(1/(nt0-1), (nt - 1)))
      psup_out[which(addarmlater>0),1] = 0
      ntj = length(which(addarmlater<=j & psup[,j+1]>0)) + 1
      if (sum(addarmlater==j)>0) ntj = ntj + sum(addarmlater==j)
      prand = rep(0, nt)
      prand[which(addarmlater<=j & psup[,j+1]>0)]  = 1
      prand[1] = 1
      prand[which(addarmlater==j)] = 1
    }
    condition_A = (platf == F | nt > 3 | j>addarmlater[nt]) & (max(psup[,j+1]) > upper)
    if (condition_A | sum(prand>0) <= 1  | length(y) >= maxN) break
    #if (sum(prand>0) <= 1  | max(psup[,j+1]) > upper | length(y) >= maxN) break
  }
  if (conjugate_prior == T) {
    if (response.type == 'absolute') {
      p.est = post0[,1]
      low = apply(post0, 1, function(z) qnorm(.025, z[1], z[2]))
      up = apply(post0, 1, function(z) qnorm(.975, z[1], z[2]))
    }
    if (response.type == "rate") {
      p.est = post0[,1]/apply(post0,1,sum)
      low = apply(post0, 1, function(z) qbeta(.025, z[1], z[2]))
      up = apply(post0, 1, function(z) qbeta(.975, z[1], z[2]))
    }
    if (response.type == "count") {
      p.est = post0[,1]/post0[,2]
      low = apply(post0, 1, function(z) qgamma(.025, z[1], z[2]))
      up = apply(post0, 1, function(z) qgamma(.975, z[1], z[2]))
    }
    est = data.frame(p.est = p.est, low = low, up = up)
  } else est = NULL
  if (good.out == F) {
    if (response.type == 'absolute') {
      theta = - theta
      est = -est
      low = est$low
      est$low = est$up
      est$up = low
    }
    if (response.type == 'rate') {
      theta = 1 - theta
      est = 1 - est
      low = est$low
      est$low = est$up
      est$up = low
    }
  }
  out = list(psup0 = psup, psup = psup_out, theta = theta, est = est, y = y, x = x[,-1], x0 = x0[,-1])
  class(out) = 'trial'
  return(out)
}


#' RAR simulator wrapper
#'
#'
#' @param i counter
#' @inheritParams RAR_sim
#' @return power and alpha indicators for one trial
#'
#' @export

sim_wrapper = function(i, nt, theta0, good.out = T, nb = 1, maxN = 500, N = 1000, upper = 0.975, uppfut = 0.95, lower = .05,
                       burn = 10*nt, response.type, conjugate_prior = T, padhere = rep(0.95,nt), adapt = T,
                       platf = T, compCon = F, MID = 0) {
  if (good.out == F) {
    if (response.type == 'absolute') theta0 = - theta0
    if (response.type == 'rate') theta0 = 1 - theta0
  }
  ng = nb
  addarmlater = rep(0, nt)
  if (platf == T) addarmlater[nt] = Inf
  j = 0
  x = array(0, dim = c(nt, 1))
  x0 = array(0, dim = c(nt, 1))
  y = NULL
  theta = array(rnorm(N*nt, 0, 10), dim = c(N, nt, 1))
  check = array(0, dim = c(N, nt))
  fcheck = array(0, dim = c(N, nt))
  if (conjugate_prior != T | response.type == "absolute") post0 = cbind(rep(0, nt), rep(10, nt))
  if (conjugate_prior == T & response.type == "rate") {
    post0 = cbind(rep(1, nt), rep(1, nt))
    p_new = apply(post0, 1, function(x) rbeta(N, x[1], x[2]))
    theta = array(- log((1 - p_new) / p_new), dim = c(N, nt, 1))
  }
  if (conjugate_prior == T & response.type == "count") {
    post0 = cbind(rep(.1, nt), rep(.001, nt))
    lambda_new = apply(post0, 1, function(x) rgamma(N, x[1], x[2]))
    theta = array(log(lambda_new), dim = c(N, nt, 1))
  }
  psup = array(rep(1/nt, nt), dim = c(nt, 1))
  nt0 = length(which(addarmlater == 0))
  psup_out = array(rep(1/nt0, nt), dim = c(nt, 1))
  psup_out[which(addarmlater>0),1] = 0
  prand = rep(1/nt, nt)
  prand[which(addarmlater>0)] = 0
  fut.stop = NULL
  repeat {
    j = j + 1
    # if (is.function(updateProgress)) {
    #   text <- paste0("updating results based on batch/patient:", j)
    #   updateProgress(detail = text)
    # }
    if (j*nb>maxN) nb = maxN - (j-1)*nb
    xb = rmultinom(nb, 1, prob = prand)
    xb0 = xb
    for (k in 1:nb) {
      xb0[which(xb0[,k] == 1), k] = rbinom(1, 1, padhere[which(xb0[,k] == 1)])
    }
    if (response.type == 'rate') {
      yb = apply(t(xb0)%*%theta0, 1, function(z) rbinom(1, 1, prob = z))
    }
    if (response.type == 'count') {
      yb = apply(t(xb0)%*%theta0, 1, function(z) rpois(1, exp(z)))
    }
    if (response.type == 'absolute') {
      yb = apply(t(xb0)%*%theta0, 1, function(z) rnorm(1, z))
    }
    x = abind(x, xb, along = 2)
    x0 = abind(x0, xb0, along = 2)
    y = c(y, yb)
    if ((response.type == 'rate' | response.type == 'count') & conjugate_prior != T) {
      smc_out = SMC_step(theta[,,j], x = xb0, y = yb, N = N, response.type)
      theta_new = smc_out$theta
      qt = apply(theta_new, 2, sd)/3
      Gibbs_out = apply(theta_new, 1, Gibbs, y = y, x = x[,-1], qt = qt,
                        response.type = response.type, ng)
      theta_new = t(sapply(Gibbs_out, function(l) return(l$theta)))
    }
    if (response.type == 'rate' & conjugate_prior == T) {
      post0 = post_beta(yb, xb, post0)
      p_new = apply(post0, 1, function(x) rbeta(N, x[1], x[2]))
      theta_new = - log((1 - p_new) / p_new)
    }
    if (response.type == 'count' & conjugate_prior == T) {
      post0 = post_gamma(yb, xb, post0)
      lambda_new = apply(post0, 1, function(x) rgamma(N, x[1], x[2]))
      theta_new = log(lambda_new)
    }
    if (response.type == 'absolute') {
      ysd = ifelse(length(y) == 1, 1, sd(y))
      post0 = post_Gauss(yb, xb, post0, ysd)
      theta_new = apply(post0, 1, function(x) rnorm(N, x[1], x[2]))
    }
    #theta_new[,which(addarmlater>j)] = -Inf
    #theta = abind(theta, theta_new, along = 3)
    
    if (compCon == F) {
      check[,which(prand>0)] = t(apply(theta_new[,which(prand>0)], 1, sup_check))
      check[,which(prand==0)] = 0
      psup_out = abind(psup_out, apply(check, 2, mean), along = 2)
      psup_out[which(addarmlater>j), j+1] = 0
      if (length(y) < burn) {
        psup = abind(psup, psup[,j], along = 2)
      } else psup = abind(psup, apply(check, 2, mean), along = 2)
      ll = NULL
      if (length(which(psup[,j+1] < lower))>0) ll = which(psup[,j+1] < lower)
      if (length(ll[addarmlater[ll]<j]) != 0) {
        fut.stop = c(fut.stop, ll[addarmlater[ll]<j])
      }
      if (!is.null(fut.stop)) {
        psup[fut.stop,j+1] = 0
        if (platf == T) addarmlater[nt] = j
      }
      ntj = length(which(addarmlater<=j & psup[,j+1]>0)) 
      if (sum(addarmlater==j)>0) ntj = ntj + sum(addarmlater==j)
      if (adapt == T & length(y)>=burn)  {
        prand = rep(0, nt)
        prand[which(addarmlater<=j & psup[,j+1]>0)] = ((ntj - 1)/ntj)*sqrt(psup[which(addarmlater<=j & psup[,j+1]>0),j+1])
        prand[which(addarmlater==j)] = 1/ntj
      } else {
        prand = rep(0, nt)
        prand[which(addarmlater<=j & psup[,j+1]>0)] = 1/ntj
        prand[which(addarmlater==j)] = 1/ntj
      }
    } else {
      mat = apply(theta_new[,which(prand>0)], 1, con_sup_check)
      if (is.null(dim(mat))) mat = matrix(mat, N, length(which(prand>0)) - 1) else mat = t(mat)
      check[,which(prand>0)] = cbind(rep(0, N), mat)
      check[,which(prand==0)] = 0
      psup_out = abind(psup_out, apply(check, 2, mean), along = 2)
      psup_out[which(addarmlater>j), j+1] = 0
      if (length(y) < burn) {
        psup = abind(psup, psup[,j], along = 2)
      } else psup = abind(psup, apply(check, 2, mean), along = 2)
      if (response.type == 'absolute') fmat = apply(theta_new[,which(prand>0)], 1, con_fut_check, MID = MID)
      if (response.type == 'rate') fmat = apply(p_new[,which(prand>0)], 1, con_fut_check, MID = MID)
      if (is.null(dim(fmat))) fmat = matrix(fmat, N, length(which(prand>0)) - 1) else fmat = t(fmat)
      fcheck[,which(prand>0)] = cbind(rep(0, N), fmat)
      fcheck[,which(prand==0)] = 0
      pfut = apply(fcheck, 2, mean)
      pfut[which(addarmlater>j)] = 0
      ff = NULL
      if (length(which(pfut > uppfut))>0) ff = which(pfut > uppfut)
      if (!is.null(ff)) {
        fut.stop = c(fut.stop, ff[addarmlater[ff]<j])
        
      }
      if (!is.null(fut.stop)) {
        psup[fut.stop,j+1] = 0
        if (platf == T) addarmlater[nt] = j
      }
      psup_out[,1] = c(0,rep(1/(nt0-1), (nt - 1)))
      psup_out[which(addarmlater>0),1] = 0
      ntj = length(which(addarmlater<=j & psup[,j+1]>0)) + 1
      if (sum(addarmlater==j)>0) ntj = ntj + sum(addarmlater==j)
      prand = rep(0, nt)
      prand[which(addarmlater<=j & psup[,j+1]>0)]  = 1
      prand[1] = 1
      prand[which(addarmlater==j)] = 1
    }
    condition_A = (platf == F | nt > 3 | j>addarmlater[nt]) & (max(psup[,j+1]) > upper)
    if (condition_A | sum(prand>0) <= 1  | length(y) >= maxN) break
  }
  psup_last = psup[,j+1]
  pbin = ifelse(psup_last > upper, 1, 0)
  pow0 = sum(pbin[which.max(theta0)] == 1)
  alpha = sum(sum(pbin) == 1)
  pow = pbin[-1]
  out = list(power0 = pow0, power = pow, alpha = alpha, N_terminate = length(y))
  
  incProgress(1)#, detail = paste("Finished simulation", i))
  
  return(out)
}


#' RAR simulator simulator
#'
#' This function simulates multiple RAR trials to calculate power and type I error rate
#' @inheritParams RAR_sim
#' @return estimated power and type I error rate
#'
#' @export
power_compute = function(nt, theta0, good.out = T, nb = 1, maxN = 500, N = 1000, upper = 0.975, uppfut = 0.95, lower = .05,
                         burn = 10*nt, response.type, conjugate_prior = T, padhere = rep(0.95,nt), adapt = T,
                         platf = T, compCon = F, MID = 0, M = 100) {
  df = data.frame(t(sapply(1:M, sim_wrapper, nt, theta0 = theta0, good.out = good.out, nb = nb, maxN = maxN, N = N, upper = upper, uppfut = uppfut, lower = lower,
                           burn = burn, response.type, conjugate_prior = conjugate_prior, padhere = padhere, adapt = adapt,
                           platf = platf, compCon = compCon, MID = MID, simplify = T)))
  if (compCon == T) power = apply(do.call(rbind, df$power), 2, mean) else power = mean(unlist(df$power0))
  out = list(power = power, Nt = as.numeric(df$N_terminate))
  return(out)
}

#' RAR simulator simulator
#'
#' This function simulates multiple RAR trials to calculate power and type I error rate
#' @inheritParams RAR_sim
#' @return estimated power and type I error rate
#'
#' @export
alpha_compute = function(nt, theta0, good.out = T, nb = 1, maxN = 500, N = 1000, upper = 0.975, uppfut = 0.95, lower = .05,
                         burn = 10*nt, response.type, conjugate_prior = T, padhere = rep(0.95,nt), adapt = T,
                         platf = T, compCon = F, MID = 0, M = 100) {
  if (response.type == 'absolute' | response.type == 'count') theta0 = rep(0, nt)
  if (response.type == 'rate') theta0 = rep(mean(theta0), nt)
  df = data.frame(t(sapply(1:M, sim_wrapper, nt, theta0 = theta0, good.out = good.out, nb = nb, maxN = maxN, N = N, upper = upper, uppfut = uppfut, lower = lower,
                           burn = burn, response.type, conjugate_prior = conjugate_prior, padhere = padhere, adapt = adapt,
                           platf = platf, compCon = compCon, MID = MID, simplify = T)))
  if (compCon == T) alpha = apply(do.call(rbind, df$power), 2, mean) else alpha = mean(unlist(df$alpha))
  out = list(alpha = alpha)
  return(out)
}


#' Sample size
#'
#' This function calculates sample size for a conventional multi-arm RCT
#'
#' @param pars vector of effect sizes, length equals number of arms
#' @param power required statistical power
#' @param alpha required type I error rate
#' @param dropout drop-out rate
#' @param type effect type "absolute" or "rate"
#' @param sigmas variances
#' @return sample size per arm \code{n} and total sample size \code{N}
#'
#' @export
sampsize = function(pars, power = .8, alpha = .05, dropout = 0.2, type, sigmas) {
  if (missing(sigmas)) {
    if (type == 'absolute') sigmas = rep(1, length(pars))
  }
  if (type == 'rate') sigmas = pars * (1 - pars)
  i = which.max(pars)
  diffs = pars[i] - pars[-i]
  diff0 = min(diffs)
  sig1 = sigmas[i]
  sig2 = sigmas[-i][which.max(pars[-i])]
  c = qnorm(1-alpha/2) + qnorm(power)
  n = ceiling((((sig1 + sig2) * c^2)/diff0^2)/(1-dropout))
  N = n * length(pars)
  out = list(n = n, N = N)
  return(out)
}

#' Frequentist pairwise power
#'
#' This function calculates power given sample size and effect size for a pairwise comparison 
#'
#' @param delta effect size, mean or proportion difference
#' @param n sample size per arm
#' @param alpha required type I error rate
#' @param type effect type "absolute" or "rate"
#' @param sigmas variances
#' @return power
#'
#' @export
power_freq = function(delta, n, alpha = 0.05, type, sigma) {
  c = qnorm(1 - alpha/2)
  z0 = delta * sqrt(n)/sqrt(sigma)
  pnorm(z0 - c) + pnorm(- z0 - c)
}

#' Frequentist multi-arm power
#'
#' This function calculates power given sample size and effect size for a multi-arm trial 
#'
#' @param n sample size per arm
#' @param pars effect sizes
#' @param alpha required type I error rate
#' @param type effect type "absolute" or "rate"
#' @param sigmas variances
#' @param target used for optimization, default is zero
#' @return power
#'
#' @export
power_freqsup = function(n, pars, alpha = 0.05, type, sigmas, target = 0){
  if (missing(sigmas)) {
    if (type == 'absolute') sigmas = rep(1, length(pars))
  }
  if (type == 'rate') sigmas = pars * (1 - pars)
  nc = choose(length(pars), 2)
  alphai = alpha/nc
  j = which.max(pars)
  delta = pars[j] - pars[-j]
  sigma = sigmas[j] + sigmas[-j]
  pows = power_freq(delta, n, alpha = alphai, type = type, sigma = sigma)
  prod(pows) - target
}

#' Frequentist compare-to-control power
#'
#' This function calculates power for comparing each arm to control given sample size and effect size for a multi-arm trial 
#'
#' @param n sample size per arm
#' @param pars effect sizes
#' @param alpha required type I error rate
#' @param type effect type "absolute" or "rate"
#' @param sigmas variances
#' @param target used for optimization, default is zero
#' @return power
#'
#' @export
power_freqcon = function(n, pars, alpha = 0.05, type, sigmas, target = 0){
  if (missing(sigmas)) {
    if (type == 'absolute') sigmas = rep(1, length(pars))
  }
  if (type == 'rate') sigmas = pars * (1 - pars)
  nc = length(pars) - 1
  alphai = alpha/nc
  j = 1
  delta = pars[j] - pars[-j]
  sigma = sigmas[j] + sigmas[-j]
  pows = power_freq(delta, n, alpha = alphai, type = type, sigma = sigma)
  pows - target
}

#' Multi-arm frequentist sample size 
#'
#' This function calculates sample size per arm for a given power effect sizes for a multi-arm trial 
#'
#' @param power target power, default is 0.8
#' @param pars effect sizes
#' @param alpha required type I error rate
#' @param type effect type "absolute" or "rate"
#' @param sigmas variances
#' @return sample size per arm
#'
#' @export
sampsize_freqsup = function(power = .8, pars, alpha = 0.05, type, sigmas) {
  if (missing(sigmas)) {
    if (type == 'absolute') sigmas = rep(1, length(pars))
  }
  if (type == 'rate') sigmas = pars * (1 - pars)
  ceiling(uniroot(power_freqsup, c(1, 100000), pars = pars, alpha = 0.05, type = type, sigmas = sigmas, 
                  target = 0.8)$root)
}


#' RAR simulator
#'
#' This function simulates a RAR trial
#'
#' @param trial object of class \code{trial}
#' @param nt number of treatments - should be teh same as in \code{trial}
#' @param theta_ae effect on secondary response
#' @param response.type response type, either 'rate', 'count' or 'absolute'
#' @return updated posterior samples
#'
#' @export

AE_sim = function(trial,nt, theta_ae, response.type, good.out = T) {
  if (good.out == F) {
    if (response.type == 'absolute') theta_ae = - theta_ae
    if (response.type == 'rate') theta_ae = 1 - theta_ae
  }
  x0 = trial$x0
  x = trial$x
  y = trial$y
  if (response.type == 'absolute') {
    post0 = cbind(rep(0, nt), rep(10, nt))
    y_ae = apply(t(x0)%*%theta_ae, 1, function(z) rnorm(1, z))
    ysd = sd(y)
    post = post_Gauss(y_ae, x, post0, ysd)
    p.est = post[,1]
    low = apply(post, 1, function(z) qnorm(.025, z[1], z[2]))
    up = apply(post, 1, function(z) qnorm(.975, z[1], z[2]))
  }
  if (response.type == "rate") {
    post0 = cbind(rep(1, nt), rep(1, nt))
    y_ae = apply(t(x0)%*%theta_ae, 1, function(z) rbinom(1, 1, prob = z))
    post = post_beta(y_ae, x, post0)
    p.est = post[,1]/apply(post,1,sum)
    low = apply(post, 1, function(z) qbeta(.025, z[1], z[2]))
    up = apply(post, 1, function(z) qbeta(.975, z[1], z[2]))
  }
  if (response.type == "count") {
    post0 = cbind(rep(.1, nt), rep(.001, nt))
    y_ae = apply(t(x0)%*%theta_ae, 1, function(z) rpois(1, z))
    post = post_gamma(y_ae, x, post0)
    p.est = post[,1]/post[,2]
    low = apply(post, 1, function(z) qgamma(.025, z[1], z[2]))
    up = apply(post, 1, function(z) qgamma(.975, z[1], z[2]))
  }
  est = data.frame(p.est = p.est, low = low, up = up)
  if (good.out == F) {
    if (response.type == 'absolute') est = -est
    if (response.type == 'rate') est = 1 - est
  }
  out = list(y = y_ae, x = x, est = est)
  class(out) = 'trial'
  return(out)
}

#' RCT simulator
#'
#' This function simulates a RCT trial
#'
#' @param nt Number of treatment arms
#' @param theta0 vector of effect sizes
#' @param maxN maximum number of patients that the trial can run for; the trial is stopped when \code{maxN} is achieved
#' @param N Number of particles
#' @param upper upper threshold for probability of superiority; the trial is stopped when one of the trial arms reaches
#' this value
#' @param response.type response type, either 'rate', 'count' or 'absolute'
#' @param conjugate_prior logical argument, default is True
#' @param padhere vector of adherence rates for each treatment
#' @return updated posterior samples
#'
#' @export
sim_wrapper_RCT = function(i, nt, theta0, maxN = 500, N = 1000, upper = 0.95, good.out = T,
                           response.type, conjugate_prior = T, padhere = rep(1,nt), compCon = F) {
  if (good.out == F) {
    if (response.type == 'absolute') theta0 = - theta0
    if (response.type == 'rate') theta0 = 1 - theta0
  }
  ng = maxN
  nb = maxN
  j = 0
  x = array(0, dim = c(nt, 1))
  y = NULL
  if (conjugate_prior != T | response.type == 'absolute') post0 = cbind(rep(0, nt), rep(10, nt))
  if (conjugate_prior == T & response.type == "rate") post0 = cbind(rep(1, nt), rep(1, nt))
  if (conjugate_prior == T & response.type == "count") post0 = cbind(rep(.1, nt), rep(.001, nt))
  psup = rep(1/nt, nt)
  xb = rmultinom(nb, 1, prob = sqrt(psup))
  xb0 = xb
  for (k in 1:nb) {
    xb0[which(xb0[,k] == 1), k] = rbinom(1, 1, padhere[which(xb0[,k] == 1)])
  }
  if (response.type == 'rate') {
    yb = apply(t(xb0)%*%theta0, 1, function(z) rbinom(1, 1, prob = z))
  }
  if (response.type == 'count') {
    yb = apply(t(xb0)%*%theta0, 1, function(z) rpois(1, exp(z)))
  }
  if (response.type == 'absolute') {
    yb = apply(t(xb0)%*%theta0, 1, function(z) rnorm(1, z))
  }
  x = xb
  y = yb
  if (response.type == 'rate' & conjugate_prior == T) {
    post0 = post_beta(yb, xb, post0)
    p_new = apply(post0, 1, function(x) rbeta(N, x[1], x[2]))
    theta_new = - log((1 - p_new) / p_new)
  }
  if (response.type == 'count' & conjugate_prior == T) {
    post0 = post_gamma(yb, xb, post0)
    lambda_new = apply(post0, 1, function(x) rgamma(N, x[1], x[2]))
    theta_new = log(lambda_new)
  }
  if (response.type == 'absolute') {
    ysd = ifelse(length(y) == 1, 1, sd(y))
    post0 = post_Gauss(yb, xb, post0, ysd)
    theta_new = apply(post0, 1, function(x) rnorm(N, x[1], x[2]))
  }
  check = t(apply(theta_new, 1, sup_check))
  ccheck = t(apply(theta_new, 1, con_sup_check))
  psupp = apply(check, 2, mean)
  cpsupp = apply(ccheck, 2, mean)
  pbin = ifelse(psupp > upper, 1, 0)
  pow = sum(pbin[which.max(theta0)] == 1)
  pow1 = ifelse(cpsupp > upper, 1, 0)
  alpha = sum(sum(pbin) == 1)
  if (compCon == T) power = pow1 else power = pow
  out = list(power = power, alpha = alpha)
  
  incProgress(1)#, detail = paste("Finished simulation", i))
  
  return(out)
}

#' RCT power estimation
#'
#' This function simulates multiple RCT trials to estimate power
#' @inheritParams sim_wrapper_RCT
#' @return estimated power
#'
#' @export
power_compute_RCT = function(nt, theta0, maxN = 500, N = 1000, upper = 0.95, good.out = T, response.type, 
                             conjugate_prior = T, padhere = rep(1, nt), compCon = F, M = 500) {
  df = data.frame(t(sapply(1:M, sim_wrapper_RCT, nt = nt, theta0 = theta0, maxN = maxN,
                           N = N, upper = upper, response.type = response.type, good.out = good.out,
                           conjugate_prior = conjugate_prior, padhere = padhere, compCon = compCon,
                           simplify = T)))
  if (compCon == T & length(df$power[[1]])>1) power = apply(do.call(rbind, df$power), 2, mean) else power = mean(unlist(df$power))
  return(power)
}

#' RCT type I error rate estimation
#'
#' This function simulates multiple RCT trials to estimate type I error rate
#' @inheritParams sim_wrapper_RCT
#' @return estimated type I error rate
#'
#' @export
alpha_compute_RCT = function(nt, theta0, maxN = 500, N = 1000, upper = 0.95, response.type, good.out = T,
                             conjugate_prior = T, compCon = F, M = 500) {
  if (response.type == 'absolute' | response.type == 'count') theta0 = rep(0, nt)
  if (response.type == 'rate') theta0 = rep(mean(theta0), nt)
  df = data.frame(t(sapply(1:M, sim_wrapper_RCT, nt = nt, theta0 = theta0, maxN = maxN,
                           N = N, upper = upper, response.type = response.type, good.out = good.out,
                           conjugate_prior = conjugate_prior, compCon = compCon, simplify = T)))
  if (compCon == T & length(df$power[[1]])>1) 
    alpha = apply(do.call(rbind, df$power), 2, mean) else alpha = mean(unlist(df$alpha))
  return(alpha)
}
