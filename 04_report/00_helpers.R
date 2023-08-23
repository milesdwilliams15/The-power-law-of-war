#############################################
## Some helper functions to make reporting ##
## results in the report/paper easier      ##
#############################################


# model fitting -----------------------------------------------------------

fit_inburr <- function(dat) {
  inburr <- function(dat, pars) {
    pars <- exp(pars)
    alpha <- pars[1]
    theta <- pars[2]
    mu <- pars[3]
    sum(-log(dinvburr(x = dat, shape1 = alpha,
                      shape2 = theta, 
                      scale = mu)))
  }
  optim(
    par = c(0, 0, 0),
    fn = inburr,
    dat = dat
  ) -> opt_out
  opt_out$par <- exp(opt_out$par)
  opt_out
}
boot_inburr <- function(dat, its = 200) {
  the_fit <- fit_inburr(dat)
  Ddat <- rank(dat) / max(rank(dat))
  Sdat <- pinvburr(
    dat, shape1 = the_fit$par[1],
    shape2 = the_fit$par[2],
    scale = the_fit$par[3]
  )
  the_D <- ks.test(
    dat,
    "pinvburr",
    shape1 = the_fit$par[1],
    shape2 = the_fit$par[2],
    scale = the_fit$par[3]
  )$statistic
  the_res <- Ddat - Sdat
  tibble(
    its = 1:its,
    bdat = map(
      its, ~ 
        tibble(
          bdat = rinvburr(
            n = length(dat),
            shape1 = the_fit$par[1],
            shape2 = the_fit$par[2],
            scale = the_fit$par[3]
          ))
    ),
    bootstrap = map(
      bdat, ~ {
        bfit <- fit_inburr(.x$bdat)
        tibble(
          alpha = bfit$par[1],
          theta = bfit$par[2],
          mu = bfit$par[3]
        )
      }
    ),
    D = map2(
      bdat, bootstrap, ~ {
        ks.test(
          .x$bdat,
          "pinvburr",
          shape1 = .y$alpha,
          shape2 = .y$theta,
          scale = .y$mu
        )$statistic
      }
    ) |> unlist()
  ) -> boot_out
  pval <- mean(boot_out$D > the_D)
  pars <- the_fit$par
  names(pars) <- c("alpha", "theta", "mu")
  list(
    dat = dat,
    pars = pars,
    bootstrap = boot_out$bootstrap |> 
      bind_rows() |>
      mutate(D = boot_out$D),
    the_D = the_D,
    pval = pval,
    the_res = the_res,
    logLik = -the_fit$value
  )
}
fit_lnorm  <- function(dat) {
  lnorm <- function(dat, pars) {
    sum(-log(
      dlnorm(dat, mean = pars[1],
            sd = exp(pars[2]))
    ))
  }
  optim(
    par = c(0, 0),
    fn = lnorm,
    dat = dat
  ) -> opt_out
  opt_out$par[2] <- exp(opt_out$par[2])
  opt_out
}
boot_lnorm <- function(dat, its = 200) {
  the_fit <- fit_lnorm(dat)
  Ddat <- rank(dat) / max(rank(dat))
  Sdat <- pnorm(
    dat, mean = the_fit$par[1],
    sd = the_fit$par[2]
  )
  the_D <- ks.test(
    log(dat), 
    "plnorm",
    meanlog = the_fit$par[1],
    sdlog = the_fit$par[2])$statistic
  the_res <- Ddat - Sdat
  tibble(
    its = 1:its,
    bdat = map(
      its, ~ 
        tibble(
          bdat = rnorm(
            n = length(dat),
            mean = the_fit$par[1],
            sd = the_fit$par[2]
          ) |> exp())
    ),
    bootstrap = map(
      bdat, ~ {
        bfit <- fit_lnorm(.x$bdat)
        tibble(
          meanlog = bfit$par[1],
          sdlog = bfit$par[2]
        )
      }
    ),
    D = map2(
      bdat, bootstrap, ~ {
        ks.test(
          log(.x$bdat), 
          "plnorm",
          meanlog = .y$meanlog,
          sdlog = .y$sdlog)$statistic
      }
    ) |> unlist()
  ) -> boot_out
  pval <- mean(boot_out$D > the_D)
  pars <- the_fit$par
  names(pars) <- c("meanlog", "sdlog")
  list(
    dat = dat,
    pars = pars,
    bootstrap = boot_out$bootstrap |> 
      bind_rows() |>
      mutate(D = boot_out$D),
    the_D = the_D,
    pval = pval,
    the_res = the_res,
    logLik = -the_fit$value
  )
}
fit_pl <- function(dat, force_xmin = NULL) {
  if(is.integer(dat)) {
    the_fit <- displ$new(dat)
    if(is.null(force_xmin)) {
      xmin <- estimate_xmin(the_fit,
                    xmins = unique(dat))
    } else {
      xmin <- force_xmin
    }
    the_fit$setXmin(xmin)
    the_fit$setPars(estimate_pars(the_fit))
    Ddat <- rank(dat) / max(rank(dat))
    Ddat <- Ddat[dat >= the_fit$xmin]
    list(
      par = c(the_fit$xmin, the_fit$pars, max(1 - Ddat)),
      value = dist_ll(the_fit)
    )
  } else {
    if(is.null(force_xmin)) {
      xmin <- estimate_xmin(the_fit,
                            xmins = unique(dat))
    } else {
      xmin <- force_xmin
    }
    the_fit <- conpl$new(dat)
    the_fit$setXmin(xmin)
    the_fit$setPars(estimate_pars(the_fit))
    Ddat <- rank(dat) / max(rank(dat))
    Ddat <- Ddat[dat >= the_fit$xmin]
    list(
      par = c(the_fit$xmin, the_fit$pars, max(1 - Ddat)),
      value = dist_ll(the_fit)
    )
  }
}
boot_pl <- function(dat, its = 200, force_xmin = NULL) {
  the_fit <- fit_pl(dat, force_xmin)
  Ddat <- rank(dat) / max(rank(dat))
  Ddat <- Ddat[dat >= the_fit$par[1]]
  Sdat <- pplcon(
    dat[dat >= the_fit$par[1]], xmin = the_fit$par[1],
    alpha = the_fit$par[2]
  )
  the_D <- ks.test(
    dat[dat >= the_fit$par[1]], 
    "pplcon",
    xmin = the_fit$par[1],
    alpha = the_fit$par[2])$statistic
  the_res <- Ddat - Sdat
  tibble(
    its = 1:its,
    bdat = map(
      its, ~ 
        tibble(
          bdat = rplcon(
            n = length(dat),
            xmin = the_fit$par[1],
            alpha = the_fit$par[2]
          ))
    ),
    bootstrap = map(
      bdat, ~ {
        bfit <- SimDesign::quiet(fit_pl(.x$bdat))
        tibble(
          xmin = bfit$par[1],
          alpha = bfit$par[2],
          pmin = bfit$par[3]
        )
      }
    ),
    D = map2(
      bdat, bootstrap, ~
        ifelse(
          sum(.x$bdat >= .y$xmin) == 1,
          NA_real_,
          ks.test(
            .x$bdat[.x$bdat >= .y$xmin], 
            "pplcon",
            xmin = .y$xmin,
            alpha = .y$alpha)$statistic
        )
    ) |> unlist()
  ) -> boot_out
  pval <- mean(boot_out$D > the_D, na.rm=T)
  pars <- the_fit$par
  names(pars) <- c("xmin", "alpha", "pmin")
  list(
    dat = dat[dat > the_fit$par[1]],
    pars = pars,
    bootstrap = boot_out$bootstrap |> 
      bind_rows() |>
      mutate(D = boot_out$D),
    the_D = the_D,
    pval = pval,
    the_res = the_res,
    logLik = the_fit$value
  )
}


# model inference ---------------------------------------------------------

inf_pl <- function(dat, its = 100) {
  fit <- fit_pl(dat)
  names(fit$par) <- c("xmin", "alpha", "pmin")
  bfit <- tibble(
    its = 1:its,
    map(
      .x = its,
      .f = ~ {
        fit <- fit_pl(
          sample(dat, length(dat), T)
        )
        tibble(
          xmin = fit$par[1],
          alpha = fit$par[2],
          logLik = fit$value
        )
      }
    )
  ) |> unnest()
  list(
    par = fit$par,
    logLik = fit$value,
    bpar = bfit
  )
}
inf_inburr <- function(dat, its = 100) {
  fit <- fit_inburr(dat)
  names(fit$par) <- c("alpha", "theta", "mu")
  bfit <- tibble(
    its = 1:its,
    map(
      .x = its,
      .f = ~ {
        fit <- fit_inburr(
          sample(dat, length(dat), T)
        )
        tibble(
          alpha = fit$par[1],
          theta = fit$par[2],
          mu = fit$par[3],
          logLik = fit$value
        )
      }
    )
  ) |> unnest()
  list(
    par = fit$par,
    logLik = fit$value,
    bpar = bfit
  )
}
inf_lnorm <- function(dat, its = 100) {
  fit <- fit_lnorm(dat)
  names(fit$par) <- c("mean", "sd")
  bfit <- tibble(
    its = 1:its,
    map(
      .x = its,
      .f = ~ {
        fit <- fit_lnorm(
          sample(dat, length(dat), T)
        )
        tibble(
          mean = fit$par[1],
          sd = fit$par[2],
          logLik = fit$value
        )
      }
    )
  ) |> unnest()
  list(
    par = fit$par,
    logLik = fit$value,
    bpar = bfit
  )
}

# compare_models ----------------------------------------------------------

## This function takes one or more of the estimated
## models and returns AIC and BIC values for each.
express_fit_pl <- function(dat, force_xmin = NULL) {
  ## fit the model
  m <- fit_pl(dat, force_xmin)
  
  ## name the fitted values
  names(m$par) <- c("xmin", "alpha", "pmin")
  
  ## calculate the D statistic
  D <- ks.test(
    x = dat[dat>=m$par["xmin"]],
    "pplcon",
    xmin = m$par["xmin"],
    alpha = m$par["alpha"]
  )$statistic
  
  ## return the data, parameters, and D
  list(
    dat = dat,
    pars = m$par,
    logLik = m$value,
    D = D
  )
}
express_fit_ln <- function(dat) {
  ## fit the model
  m <- fit_lnorm(dat)
  
  ## name the parameters
  names(m$par) <- c("meanlog", "sdlog")
  
  ## calculate the D statistic
  D <- ks.test(
    x = dat,
    "plnorm",
    meanlog = m$par["meanlog"],
    sdlog = m$par["sdlog"]
  )$statistic
  
  ## return the data, model parameters, and D
  list(
    dat = dat,
    pars = m$par,
    logLik = -m$value,
    D = D
  )
}
express_fit_ib <- function(dat) {
  ## fit the model
  m <- fit_inburr(dat)
  
  ## name the parameters
  names(m$par) <- c("alpha", "theta", "mu")
  
  ## calculate the D statistic
  D <- ks.test(
    x = dat,
    "pinvburr",
    shape1 = m$par["alpha"],
    shape2 = m$par["theta"],
    scale = m$par["mu"]
  )$statistic
  
  ## return the data, the parameters, and D
  list(
    dat = dat,
    pars = m$par,
    logLik = - m$value,
    D = D
  )
}

compare_models <- function(dat, its = 200, force_xmin = NULL) {
  ## Fit each of the models to the data
  ## (only fit lnorm and invburr for x >= xmin)
  pl_fit <- express_fit_pl(dat, force_xmin)
  xmin <- pl_fit$pars[1]
  ib_fit <- express_fit_ib(dat[dat>=xmin])
  ln_fit <- express_fit_ln(dat[dat>=xmin])
  out <- tibble(
    Model = c("Power-law", "Inverse Burr", "Log-normal"),
    D = c(pl_fit$D, ib_fit$D, ln_fit$D)
  )

  ## Now bootstrap
  boot_out <- tibble(
    its = 1:its,
    out = map(
      its, ~ {
        bdat <- sample(dat, length(dat), T)
        pl_fit <- suppressWarnings( 
          express_fit_pl(bdat, force_xmin))
        xmin <- pl_fit$pars[1]
        ib_fit <- suppressWarnings(
          express_fit_ib(bdat[bdat>=xmin]))
        ln_fit <- suppressWarnings(
          express_fit_ln(bdat[bdat>=xmin]))
        the_Ds <- c(
          pl_fit$D, ib_fit$D, ln_fit$D
        )
        tibble(Model = c(
          "Power-law", "Inverse Burr", "Log-normal"
        ),
        bootD = the_Ds
        )
      }
    ))
  
  boot_out <- boot_out$out |>
    bind_rows() |>
    group_by(Model) |>
    summarize(se = sd(bootD),
              lo = quantile(bootD, 0.09),
              hi = quantile(bootD, 0.91))
    
  ## return
  full_join(out, boot_out, by = "Model")
  
}


# plotting the distributions ----------------------------------------------

## The following plot and line functions make it
## easy to plot Pr(X > x) over x in log-log space
## and overlay the fitted models.

my_plot <- function(dat, pch = 19, col = "gray", 
                    xlab = "x",
                    ylab = "Pr(X > x)",
                    main = NULL) {
  # plot(
  #   x = dat,
  #   y = 1 - rank(dat) / max(rank(dat)),
  #   log = "xy",
  #   pch = pch,
  #   col = col,
  #   ...
  # )
  my_new_plot <<- ggplot() +
    aes(
      x = dat,
      y = 1 - rank(dat) / max(rank(dat))
    ) +
    geom_point(
      color = col,
      shape = pch
    ) +
    scale_x_log10(
      labels = scales::comma
    ) +
    scale_y_log10() +
    labs(
      x = xlab,
      y = ylab,
      title = main
    )
  #print(my_new_plot)
}
my_lines <- function(m, col = "black", lty = 1, lsz = 1) {
  if("mu" %in% names(m$pars)) {
    my_new_plot <<- my_new_plot +
      geom_line(
        aes(
          x = sort(m$dat),
          y = 1 - pinvburr(
            sort(m$dat),
            m$pars[1],
            m$pars[2],
            scale = m$pars[3]
          ),
          color = col
        ),
        linetype = lty,
        linesize = lsz
      )
    # lines(
    #   x = sort(m$dat),
    #   y = 1 - pinvburr(
    #     sort(m$dat),
    #     m$pars[1],
    #     m$pars[2],
    #     scale = m$pars[3]
    #   ),
    #   #log = "xy",
    #   col = col,
    #   lty = lty,
    #   ...
    # )
  } else if("meanlog" %in% names(m$pars)) {
    my_new_plot <<- my_new_plot +
      geom_line(
        aes(
          x = sort(m$dat),
          y = 1 - pnorm(
            log(sort(m$dat)),
            mean = m$pars[1],
            sd = m$pars[2]
          ),
          color = col
        ),
        linetype = lty,
        linesize = lsz
      )
    # lines(
    #   x = sort(m$dat),
    #   y = 1 - pnorm(
    #     log(sort(m$dat)),
    #     mean = m$pars[1],
    #     sd = m$pars[2]
    #   ),
    #   #log = "xy",
    #   col = col,
    #   lty = lty,
    #   ...
    # )
  } else {
    my_new_plot <<- my_new_plot +
      geom_line(
        aes(
          x = sort(m$dat[m$dat >= m$pars[1]]),
          y = (1 - pplcon(
            sort(m$dat[m$dat >= m$pars[1]]),
            xmin = m$pars[1],
            alpha = m$pars[2]
          )) * m$pars[3],
          color = col
        ),
        linetype = lty,
        linesize = lsz
      )
    # lines(
    #   x = sort(m$dat[m$dat >= m$pars[1]]),
    #   y = (1 - pplcon(
    #     sort(m$dat[m$dat >= m$pars[1]]),
    #     xmin = m$pars[1],
    #     alpha = m$pars[2]
    #   )) * m$pars[3],
    #   col = col,
    #   lty = lty,
    #   ...
    # )
  }
  #print(my_new_plot)
}



# inverse burr regression -------------------------------------------------

## This function performs an inverse Burr regression
## where covariates determine the value of mu, while the 
## shape parameters theta and alpha are estimated as constants.

inbur_reg <- function(formula, data = NULL) {
  
  ## The Data
  y  <- model.frame(formula, data)[, 1]
  x  <- model.matrix(formula, data)
  
  ## The likelihood
  inbur_lik <- function(x, y, pars) {
    b <- matrix(
      data = pars[1:ncol(x)],
      nrow = ncol(x),
      ncol = 1
    )
    mu <- exp(x %*% b)
    alpha <- exp(pars[ncol(x) + 1])
    theta <- exp(pars[ncol(x) + 2])
    sum(
      - log(
        dinvburr(
          y, 
          shape1 = alpha,
          shape2 = theta,
          scale = mu
        )
      )
    )
  }
  
  ## Estimation
  optim(
    par = rep(0, len = ncol(x) + 2),
    fn = inbur_lik,
    x = x,
    y = y,
    hessian = F
  ) -> opt_out
  
  ## Bootstrapping
  tibble(
    its = 1:200,
    bout = map(
      its,
      ~ {
        bkeep <- sample(1:nrow(x), nrow(x), T)
        optim(
          par = rep(0, len = ncol(x) + 2),
          fn = inbur_lik,
          x = x[bkeep, ],
          y = y[bkeep],
          hessian = F
        ) -> opt_out
        tibble(
          pars = 1:length(opt_out$par),
          vals = opt_out$par
        )
      }
    )
  ) |>
    unnest(cols = bout) |>
    group_by(pars) |>
    summarize(
      std.error = sd(vals, na.rm=T)
    ) -> boot_se
  
  list(
    out = tibble(
      term = c(colnames(x),
               "log(alpha)","log(theta)"),
      estimate = opt_out$par,
      std.error = boot_se$std.error,
      statistic = estimate / std.error,
      p.value = 1 - pnorm(
        abs(statistic)
      ) |> round(3)
    ),
    logLik = -opt_out$value,
    dat = model.frame(formula, data)
  )
}

# get aic and bic for inverse burr regression models ----------------------

## These functions return the AIC and BIC for inverse
## Burr models.

ib_aic <- function(m) {
  logLik <- m$logLik
  n <- nrow(m$dat)
  k <- nrow(m$out)
  - 2 * logLik  + 2 * k 
}
ib_bic <- function(m) {
  logLik <- m$logLik
  n <- nrow(m$dat)
  k <- nrow(m$out)
  - 2 * logLik + k * log(n)
}



# get predictions for inverse burr model ----------------------------------

ib_predict <- function(m) {
  tibble(
    var = 1:ncol(m$dat[, -1]),
    term = colnames(m$dat[, -1]),
    pred = map(
      var,
      ~ {
        var_names <- colnames(m$dat)[-1]
        the_var <- var_names[.x]
        pd <- m$dat |>
          mutate(
            across(var_names[-.x], mean)
          )
        X <- cbind(1, pd[,-1]) |> as.matrix()
        b <- m$out$estimate[1:ncol(X)]
        mu <- exp(X %*% b)
        alpha <- exp(m$out$estimate[ncol(X)+1])
        theta <- exp(m$out$estimate[ncol(X)+2])
        pred <- qinvburr(0.5, alpha, theta, scale = mu)
        out <- tibble(
          #var = var_names[.x],
          val = X[, 1+.x],
          fit = pred
        )
        ## caculate the SE of the fit
        tibble(
          x = 1:nrow(X),
          map(
            x,
            ~ {
              mus <- 0
              for(i in 1:200) {
                b <- m$out$estimate[1:ncol(X)] +
                  rnorm(ncol(X), 0, sd = m$out$std.error[1:ncol(X)]) *
                  (colnames(X) == the_var)
                mus[i] <- exp(X[.x, ] %*% b)
              }
              pred <- qinvburr(0.5, alpha, theta, scale = mus)
              tibble(
                se = sd(pred, na.rm=T)
              )
            }
          )
        ) |> unnest(2) -> sim_out
        bind_cols(out, sim_out) |>
          mutate(
            lwr = fit - 1.96 * se,
            upr = fit + 1.96 * se
          )
      }
    )
  ) |> unnest(pred)
}
