
#' Posterior estimates table
#'
#' @param trial An object of class 'trial'
#' @param credibility Level of credibility for credible intervals, default is 0.95
#' @param type a character representing the type of the parameter; "absolute", "rate",
#' or "count"
#' @return A data frame of effect size estimates and credible intervals
#' @export

est = function(trial, credibility = 0.95, type) {
  a = (1 - credibility) / 2
  cr = credibility*100
  nt = nrow(trial$x)
  if (!is.null(trial$est)) {
    es = trial$est
    p.est = es$p.est
    lower = es$low
    upper = es$up
    CI = paste("(", round(lower,2), ",", "", round(upper,2), ")", sep = ' ')
    out = data.frame(round(p.est,2), CI)
    names(out) = c("Estimate", paste(cr, "%", "CI", sep = ""))
    row.names(out) = paste("treatment", 1:nt)
  } else {
    if (type == "absolute") {
      theta = trial$theta
      np = dim(theta)[3]
      nt = dim(theta)[2]
      theta = theta[,,np]
      est = apply(theta, 2, mean)
      lower = apply(theta, 2, quantile, p = a)
      upper = apply(theta, 2, quantile, p = (1-a))
      CI = paste("(", round(lower,2), ",", "", round(upper,2), ")", sep = ' ')
      out = data.frame(round(est,2), CI)
      names(out) = c("Estimate", paste(cr, "%", "CI", sep = ""))
      row.names(out) = paste("treatment", 1:nt)
    }
    if (type == "rate") {
      theta = exp(trial$theta)/(1 + exp(trial$theta))
      np = dim(theta)[3]
      nt = dim(theta)[2]
      theta = theta[,,np]
      est = apply(theta, 2, mean)
      lower = apply(theta, 2, quantile, p = a)
      upper = apply(theta, 2, quantile, p = (1-a))
      CI = paste("(", round(lower,2), ",", "", round(upper,2), ")", sep = ' ')
      out = data.frame(round(est,2), CI)
      names(out) = c("Estimate (Rate)", paste(cr, "%", "CI", sep = ""))
      row.names(out) = paste("treatment", 1:nt)
    }
    if (type == "count") {
      theta = exp(trial$theta)
      np = dim(theta)[3]
      nt = dim(theta)[2]
      theta = theta[,,np]
      est = apply(theta, 2, mean)
      lower = apply(theta, 2, quantile, p = a)
      upper = apply(theta, 2, quantile, p = (1-a))
      CI = paste("(", round(lower,2), ",", "", round(upper,2), ")", sep = ' ')
      out = data.frame(round(est,2), CI)
      names(out) = c("Estimate (Intensity)", paste(cr, "%", "CI", sep = ""))
      row.names(out) = paste("treatment", 1:nt)
    }
  }
  return(out)
}


#' @param trial An object of class 'trial'
#'
#' @return A data frame of probabilities of superiority
#' @export

psup = function(trial) {
  psup = trial$psup
  nt = dim(psup)[1]
  np = dim(psup)[2]
  psup = psup[,np]
  out = data.frame(psup)
  names(out) = c("Probability of Superiosity")
  row.names(out) = paste("treatment", 1:nt)
  return(out)
}


#' @param trial An object of class 'trial'
#'
#' @return A data frame of probabilities of superiority
#' @export

datat = function(trial) {
  y = trial$y
  x = trial$x
  nt = nrow(x)
  Narm = rowSums(x)
  Parm = round(rowSums(x)/length(y), 2)
  smean = round((x %*% y) / Narm, 2)
  out = data.frame(smean, Narm, Parm)
  names(out) = c("sample mean", "number of assignments", "proportion of assignments")
  row.names(out) = paste("treatment", 1:nt)
  return(out)
}


