# ======================================================
# Empirical analysis for section:
# "Implications: Parameterizing Correlates of War Size"
# ======================================================

# Make a function that will estimate an inverse Burr regression
inbur_reg <- function(formula, data = NULL, its = 2000) {
  
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
    its = 1:its,
    bout = future_map(
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
      },
      .options = furrr_options(seed = T)
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

# now estimate models with bootstrapped SEs:
# the right-hand side, which includes:
# - log of world population
# - log of combined population of countries fighting a war
# - log of military personnel
# - worst V-Dem score among fighting countries
# - post 1950 indicator
rhs <- ~ log(wpop) + log(tpop) + 
  log(milper) + min_polyarchy + I(year > 1950)

# estimate inverse Burr models for each outcome
set.seed(1)
ib_fit1 <- inbur_reg(
  update(rhs, batdeath ~ .), data = dt
)$out
ib_fit2 <- inbur_reg(
  update(rhs, 1e06 * batdeath / wpop ~ .), data = dt
)$out
ib_fit3 <- inbur_reg(
  update(rhs, 1e06 * batdeath / tpop ~ .), data = dt
)$out

# estimate the log-linear models via OLS with robust SEs
ln_fit1 <- lm_robust(
  update(rhs, log(batdeath) ~ .), data = dt
) |> tidy()
ln_fit2 <- lm_robust(
  update(rhs, log(1e06 * batdeath / wpop) ~ .), data = dt
) |> tidy()
ln_fit3 <- lm_robust(
  update(rhs, log(1e06 * batdeath / tpop) ~ .), data = dt
) |> tidy()

# report the results as a coefficient plot
set_palette(
  binary = qual[c(1, 4)],
  from_coolors = F
)
bind_rows(
  ib_fit1 |> mutate(model = "Inverse Burr",
                    outcome = "Severity"),
  ib_fit2 |> mutate(model = "Inverse Burr", 
                    outcome = "Prevalence"),
  ib_fit3 |> mutate(model = "Inverse Burr",
                    outcome = "Intensity"),
  ln_fit1 |> mutate(model = "Log-normal", 
                    outcome = "Severity"),
  ln_fit2 |> mutate(model = "Log-normal", 
                    outcome = "Prevalence"),
  ln_fit3 |> mutate(model = "Log-normal", 
                    outcome = "Intensity")
) |> filter(
  !(term %in% c("log(alpha)", "log(theta)", "(Intercept)"))
) |>
  mutate(
    outcome = factor(
      outcome, levels = c("Severity", "Prevalence", "Intensity")
    )
  ) |>
  select(term:p.value, model:outcome) |>
  ggplot() +
  aes(
    x = estimate,
    y = term,
    xmin = estimate - 1.96 * std.error,
    xmax = estimate + 1.96 * std.error,
    color = model
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  geom_pointrange(
    position = ggstance::position_dodgev(-.5),
    size = .3
  ) +
  geom_vline(
    xintercept = -11,
    color = "gray80",
    linewidth = 30
  ) +
  geom_text(
    aes(x = -11,
        label = paste0(round(estimate, 2), 
                       gtools::stars.pval(p.value))),
    position = ggstance::position_dodgev(-.75),
    show.legend = F,
    fontface = "bold",
    hjust = 0.2
  ) +
  facet_wrap(~ outcome, scales = "free_x") +
  ggpal(type = "binary") +
  scale_y_discrete(
    labels = c("Post-1950", "Military Size (ln)",
               "Belligerent Pop. (ln)", "Global Pop. (ln)",
               "Democracy")
  ) +
  labs(
    x = "Coefficient with 95% CIs",
    y = NULL,
    color = NULL,
    caption = "*** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.1"
  ) +
  ggthemes::theme_few() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(
      size = 16
    ),
    axis.text = element_text(
      size = 14
    ),
    axis.title = element_text(
      size = 14
    ),
    legend.text = element_text(
      size = 14
    ),
    panel.grid.major.y = element_line(
      linetype = 3,
      color = "gray30"
    )
  ) 