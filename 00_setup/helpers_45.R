# Helpers for power laws and alternatives #


# Pr(X > x) function ------------------------------------------------------

PXx <- function(x) unlist(
  lapply(x, function(X) mean(X < x, na.rm=T))
)

# The logit function ------------------------------------------------------

logit <- function(x) log(x / (1 - x))
logistic <- function(x) 1 / (1 + exp(-x))

# The power law estimator

pl.est <- function(x) {
  xmins <- unique(x)
  dat   <- numeric(length(xmins))
  z     <- sort(x)
  for(i in 1:length(xmins)) {
    xmin <- xmins[i]
    z1   <- z[z>=xmin]
    n    <- length(z1)
    a    <- 1 + n*(sum(log(z1/xmin)))^-1
    cx   <- (n:1)/n
    cf   <- (z1 / xmin)^(-a+1)
    dat[i] <- max(abs(cf-cx))
  }
  D <- min(dat[dat>0], na.rm=T)
  xmin <- min(xmins[which(dat==D)])
  z <- x[x>=xmin]
  z <- sort(z)
  n <- length(z)
  beta <- 1 + n*(sum(log(z/xmin)))^(-1)
  c(xmin = xmin, beta = beta)
}

pl.boot <- function(x, R) {
  replicate(
    n = R,
    list(
      pl.est(sample(x, length(x), T))
    )
  ) -> out
  bind_rows(out) %>%
    mutate(rep = 1:n())
}

pl <- function(x, R = 200) {
  est  <- pl.est(x)
  boot <- pl.boot(x, R = R)
  list(
    estimate  = est,
    bootstrap = boot
  )
}


# logit-log form ----------------------------------------------------------

ll.est <- function(x) {
  px <- PXx(x)
  data <- data.frame(x = log(x), px = px) %>%
    filter(x != -Inf & x != Inf)
  out <- coef(glm(px ~ x, data, family = quasibinomial))
  out[2] <- - 1 / out[2]
  out[1] <- out[1] * out[2]
  names(out) <- c("mu", "sigma")
  out
}

ll.boot <- function(x, R) {
  replicate(
    n = R,
    list(ll.est(sample(x, length(x), T)))
  ) -> out
  bind_rows(out) %>%
    mutate(rep = 1:n())
}

ll <- function(x, R = 200) {
  est <- ll.est(x)
  boot <- ll.boot(x, R = R)
  list(
    estimate = est,
    bootstrap = boot
  )
}

# log-normal --------------------------------------------------------------

ln.est <- function(x) {
  data <- data.frame(x = log(x)) %>%
    filter(x != -Inf & x != Inf)
  out <- c(mean(data$x, na.rm=T), sd(data$x, na.rm=T))
  names(out) <- c("mu", "sigma")
  out
}

ln.boot <- function(x, R) {
  replicate(
    n = R,
    list(ln.est(sample(x, length(x), T)))
  ) -> out
  bind_rows(out) %>%
    mutate(rep = 1:n())
}

ln <- function(x, R = 200) {
  est <- ln.est(x)
  boot <- ln.boot(x, R = R)
  list(
    estimate = est,
    bootstrap = boot
  )
}

# compare -----------------------------------------------------------------

