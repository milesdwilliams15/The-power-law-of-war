#############################################################
#Reference Class definition for inverse Burr
#############################################################
library(poweRlaw)
library(actuar)
coninvburr =
  setRefClass("coninvburr",
              contains = "ctn_distribution",
              fields = list(
                dat = function(x) {
                  if (!missing(x) && !is.null(x)) {
                    #check_ctn_data(x)
                    d = sort(x)
                    internal[["cum_n"]] <<- rev(seq_along(d))
                    internal[["dat"]] <<- d
                    xmin <<- d[1]
                  } else internal[["dat"]]
                },
                xmin = function(x) {
                  if (!missing(x) && !is.null(x)) {
                    if ("estimate_xmin" %in% class(x)) {
                      pars <<- x$pars
                      x = x$xmin
                    }
                    internal[["xmin"]] <<- x
                    if (length(internal[["dat"]])) {
                      selection = min(which(internal[["dat"]] >= (x - .Machine$double.eps ^ 0.5)))
                      internal[["n"]] <<- internal[["cum_n"]][selection]
                    }
                  } else  internal[["xmin"]]
                },
                pars = function(x) {
                  if (!missing(x) && !is.null(x)) {
                    if ("estimate_pars" %in% class(x)) x = x$pars
                    internal[["pars"]] <<- x
                  } else internal[["pars"]]
                }
              )
  )
#############################################################
#Initialisation
#############################################################
coninvburr$methods(
  list(
    initialize = function(dat) {
      no_pars <<- 3
      ##Use the internal attribute for copying
      if (!missing(dat)) {
        #check_ctn_data(dat)
        d = sort(dat)
        internal[["cum_n"]] <<- rev(seq_along(d))
        internal[["dat"]] <<- d
        xmin <<- d[1]
      }
    }
  )
)

#############################################################
#PDF method
#############################################################

setMethod("dist_pdf",
          signature = signature(m = "coninvburr"),
          definition = function(m, q = NULL, log = FALSE) {
            xmin = m$getXmin(); pars = m$getPars()
            if (is.null(q)) q = m$dat
            
            pdf = log(dinvburr(q, pars[1], pars[2], scale = pars[3])) -
              log(pinvburr(xmin, pars[1], pars[2], scale = pars[3], lower.tail = FALSE))
            if (!log) {
              pdf = exp(pdf)
              pdf[q < xmin] = 0
            } else {
              pdf[q < xmin] = -Inf
            }
            pdf
          }
)

#############################################################
#CDF method
#############################################################
setMethod("dist_cdf",
          signature = signature(m = "coninvburr"),
          definition = function(m, q = NULL, lower_tail = TRUE) {
            pars = m$pars; xmin = m$xmin
            if (is.null(pars)) stop("Model parameters not set.")
            if (is.null(q)) q = m$dat
            
            if (lower_tail) {
              p = pinvburr(q, pars[1], pars[2], scale = pars[3], lower.tail = lower_tail)
              C = pinvburr(xmin, pars[1], pars[2], scale = pars[3], lower.tail = FALSE)
              pdf = (p / C - 1 / C + 1)
            } else {
              log_p = pinvburr(q, pars[1], pars[2], scale = pars[3], lower.tail = FALSE, log.p = TRUE)
              log_C = pinvburr(xmin, pars[1], pars[2], scale = pars[3], lower.tail = FALSE, log.p = TRUE)
              pdf = exp(log_p - log_C)
            }
            pdf[q < xmin] = 0
            pdf
          }
)

setMethod("dist_all_cdf",
          signature = signature(m = "coninvburr"),
          definition = function(m, lower_tail = TRUE, xmax = 1e5) {
            xmin = m$getXmin()
            xmax = min(max(m$dat), xmax)
            dist_cdf(m, q = xmin:xmax, lower_tail = lower_tail)
          }
)

#############################################################
#ll method
#############################################################
setMethod("dist_ll",
          signature = signature(m = "coninvburr"),
          definition = function(m) {
            q = m$dat
            n = m$internal[["n"]]; N = length(q)
            q = q[(N - n + 1):N]
            
            coninvburr_tail_ll(q, m$getPars(), m$getXmin())
          }
)

########################################################
#Log-likelihood
########################################################
coninvburr_tail_ll = function(x, pars, xmin) {
  if (is.vector(pars)) pars = t(as.matrix(pars))
  n = length(x)
  joint_prob = colSums(apply(exp(pars), 1,
                             function(i) dinvburr(x, i[1], i[2], scale = i[3], log = TRUE)))
  
  prob_over = apply(exp(pars), 1, function(i)
    pinvburr(xmin, i[1], i[2], scale = i[3], log.p = TRUE, lower.tail = FALSE))
  joint_prob - n * prob_over
}

########################################################
#Rand number generator
########################################################
setMethod("dist_rand",
          signature = signature(m = "coninvburr"),
          definition = function(m, n = "numeric") {
            xmin = m$getXmin(); pars = m$getPars()
            rns = numeric(n)
            i = 0; N = 0
            tail_prob = pinvburr(xmin, pars[1L], pars[2L], scale = pars[3L], lower.tail = FALSE)
            if ((1 / tail_prob) > 10e10) {
              stop("It appears that your parameters put in you in the __very__ extreme tail
                     of the inverse Burr distribution. This means it is impossible to generate
                     random numbers in a finite amount of time.")
            }
            ## n-0.5 to avoid floating point silliness.
            while (i < (n - 0.5)) {
              ## Since we reject RNs less than xmin we should simulate N > n rns
              ## If we simulate N Rns (below), we will keep n-i (or reject N-(n-i))
              N = ceiling((n - i) / tail_prob)
              
              ## Simple rejection sampler
              x = rinvburr(N, pars[1L], pars[2L], scale = pars[3L])
              x = x[x > xmin]
              if (length(x)) {
                x = x[1:min(length(x), n - i)]
                rns[(i + 1L):(i + length(x))] = x
                i = i + length(x)
              }
            }
            rns
          }
)

#############################################################
#MLE method
#############################################################
coninvburr$methods(
  mle = function(set = TRUE, initialise = NULL) {
    x = dat
    x = x[x > xmin]
    if (is.null(initialise)) theta_0 = c(0, 0, 0) else
      theta_0 = initialise
    # Chop off values below
    # negloglike = function(par) {
    #   r = -coninvburr_tail_ll(x, par, xmin)
    #   if (!is.finite(r)) r = 1e12
    #   r
    # }
    negloglike = function(x, pars) {
      ps = exp(pars)
      alpha = ps[1]
      theta = ps[2]
      mu <- ps[3]
      sum(-log(dinvburr(x = x, shape1 = alpha,
                        shape2 = theta, 
                        scale = mu)))
    }
    mle = suppressWarnings(optim(par = theta_0,
                                 x = x,
                                 fn = negloglike))
    mle$par = exp(mle$par)
                                 # method = "L-BFGS-B",
                                 # lower = c(-Inf, .Machine$double.eps)))
    if (set)
      pars <<- exp(mle$par)
    class(mle) = "estimate_pars"
    names(mle)[1L] = "pars"
    #mle$pars = pars
    mle
    
  }
)



# need to create a new bootstrap_p function --------------------------------

sample_p_helper = function(i, m, x_lower) {
  ## Total sample size
  N = get_n(m)
  ntail_prop = get_ntail(m, prop = TRUE)
  
  ## Proportion to sample
  n1 = sum(runif(N) > ntail_prop) # less than xmin
  
  # q should be of length N
  c(sample(x_lower, n1, replace = TRUE), #less than xmin
    dist_rand(m, N - n1))
}

bootstrap_p_helper = function(i, m, x_lower, xmins, pars, xmax, distance) {
  
  q = sample_p_helper(i, m, x_lower)
  m_cpy = m$getRefClass()$new(q)
  
  est = estimate_xmin(m_cpy, xmins = xmins, pars = pars, xmax = xmax, distance = distance)
  ## Remove the character now, since we will change to data frame.
  est["distance"] = NULL
  unlist(est)
}


get_bootstrap_p_sims = function(m, no_of_sims, seed, threads = 1) {
  
  if (is.null(m$getPars())) {
    stop("Parameters need to be set. See ?estimate_xmin")
  }
  # cl = parallel::makeCluster(threads)
  # on.exit(parallel::stopCluster(cl))
  
  x = m$dat
  x_lower = x[x < m$xmin]
  ## Set cluster seed
  # parallel::clusterSetRNGStream(cl, seed)
  # parallel::clusterExport(cl, "get_n")
  sapply(1:no_of_sims, sample_p_helper, m = m, x_lower = x_lower)
}


ib_bootstrap_p = function(m, xmins = NULL, pars = NULL, xmax = 1e9,
                       no_of_sims = 100, threads = 1,
                       seed = NULL, distance = "ks") {
  
  if (is.null(m$getPars())) {
    message("Parameters will be initially estimated via estimate_xmin")
  }
  
  m_cpy = m$copy()
  # time = timer()
  # time$start()
  # gof_v = estimate_xmin(m_cpy, xmins = xmins, pars = pars,
  #                       xmax = xmax, distance = distance)
  gof_v <- list(gof = NULL)
  gof_v$gof = suppressWarnings(ks.test(
    m_cpy$dat,
    "pinvburr",
    shape1 = m_cpy$pars[1],
    shape2 = m_cpy$pars[2],
    scale = m_cpy$pars[3]
  )$statistic)
  #time$stop()
  if (is.na(gof_v$gof) || is.infinite(gof_v$gof)) {
    stop("Unable to estimate initial xmin using estimate_xmin(), so we can't bootstrap.")
  }
  
  if (min(m_cpy$dat) > xmax) {
    stop("The smallest value in your data set is larger than xmax. The xmax
         parameter is the upper limit of the xmin search space.")
  }
  
  if (max(m_cpy$dat) > xmax) {
    message("Some of your data is larger than xmax. The xmax parameter is
            the upper bound of the xmin search space. You could try increasing
            it. If the estimated values are below xmax, it's probably OK not to
            worry about this.")
  }
  
  # message("Expected total run time for ", no_of_sims,
  #         " sims, using ", threads, " threads is ",
  #         signif(time$get() * no_of_sims / threads, 3), " seconds.")
  # 
  if (is.null(m$getPars()) || is.null(m$getXmin())) m_cpy$setXmin(gof_v)
  
  x = m_cpy$dat
  x_lower = x[x < m_cpy$xmin]
  
  ## Start clock and parallel bootstrap
  # time$start()
  # cl = parallel::makeCluster(threads)
  # on.exit(stopCluster(cl))
  
  ## Set cluster seed
  #if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
  
  #parallel::clusterExport(cl, c("dist_rand", "estimate_xmin"))
  # nof = sapply(1:no_of_sims,
  #                 bootstrap_p_helper, m_cpy,
  #                 x_lower, xmins, pars, xmax, distance)
  nof = tibble(
    its = 1:no_of_sims,
    nof = map(
      its, ~ {
        simx <- dist_rand(m_cpy, n = length(m_cpy$dat))
        simm <- coninvburr$new(simx)
        simf <- suppressWarnings(estimate_pars(simm))
        simD <- ks.test(
          simx,
          "pinvburr",
          shape1 = simf$pars[1],
          shape2 = simf$pars[2],
          scale = simf$pars[3]
        )$statistic
        tibble(
          gof = simD,
          xmin = simm$xmin,
          pars1 = simf$pars[1],
          pars2 = simf$pars[2],
          pars3 = simf$pars[3],
          ntail = length(simm$dat)
        )
      }
    )
  )
  ## Stop clock and cluster
  #total_time = time$get(stop = TRUE) * threads
  
  #bootstraps = as.data.frame(t(nof))
  bootstraps = nof$nof |> bind_rows() |>
    as.data.frame()
  l = list(p = sum(bootstraps$gof >= gof_v[["gof"]], na.rm = TRUE) / no_of_sims,
           gof = gof_v[["gof"]],
           bootstraps = bootstraps,
           #sim_time = total_time[[1]] / no_of_sims,
           seed = seed,
           package_version = packageVersion("poweRlaw"),
           distance = distance)
  class(l) = "bs_p_xmin"
  l
}

my_bootstrap_p <- function(m, xmins = NULL, pars = NULL,
               xmax = 1e09, no_of_sims = 100,
               threads = 1, seed = NULL, distance = "ks") {
  if(class(m)%in%"coninvburr") {
    message("Running for 'coninvburr' class, so no parallel bootstrapping.")
    ib_bootstrap_p(m, xmins,pars, xmax, no_of_sims,
                   threads, seed, distance)
  } else {
    bootstrap_p(m, xmins, pars, xmax, no_of_sims,
                threads, seed, distance)
  }
}

