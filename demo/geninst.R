library(TLPR)

get_adjustment_weights <- function(env) {
  # Generate powers of nSI from 1 to nI
  si_powers <- with(env, nSI ^ (1L:nI))
  # Calculate the base value (nSI ^ nI)
  base_value <- with(env, nSI ^ nI)
  # Generate products of base_value and powers of nSJ from 0 to (nJ - 1)
  sj_powers <- base_value * with(env, (nSJ ^ (1L:(nJ - 1L))))
  # Combine the results
  c(1L, si_powers, sj_powers)
}

env <- new.env()
generate_cssap(env, rate = 4.0, tau = 4L, nB = 5L, nCS =10L, nCO = 2L, nI = 3L, nJ = 2L)

# apply(env$Cb, 1L, function(x) x, simplify = "array") -> env$Cbx

env$winner    |>
  seq_along() |>
  as.list()   |>
  (\(lst) {names(lst) <- names(env$winner); lst})() |>
  assign(x = "carrierIdx", envir = env)

env$winner |>
  names()  |>
  assign(x = "winnerKey", envir = env)

data.frame(Carrier = names(env$winner)[env$carriers[seq(env$nL_)]], CTb = env$CTb) |>
  (\(df) split(df$CTb, df$Carrier))() |>
  assign(x = "CTb_list", pos = -1L, envir = env)

# data.frame(Carrier = rep(names(env$winner), each = env$tau), Cb = c(env$Cb)) |>
#   (\(df) split(df$Cb, df$Carrier))() |>
#   assign(x = "Cbxx", pos = -1L, envir = env)

with(
  env, {
    max.S <- R <- 10L

    nSI  <- max.S + 1L
    nSJ  <- 2L * max.S + 1L
    nSdx <- (nSI^nI)*(nSJ^nJ)
    nA   <- 4L
    alpha <- c(rep(5.0, env$nI), rep(c(4.0, 8.0), each = env$nJ))
    Q <- list(
      "vals" = c(0L, 4L, 8L),
      "prob" = c(0.4, 0.3, 0.3)
    )
    D <- list(
      "vals" = c(0L, 4L, 8L),
      "prob" = c(0.25, 0.25, 0.5)
    )
    W <- list(
      "vals" = c(7.0, 22.0),
      "prob" = dnorm(c(7.0, 22.0), 18.0, 12.0) / sum(dnorm(c(7.0, 22.0), 18.0, 12.0))
    )
    nQ <- length(Q$vals)
    nD <- length(D$vals)
    nW <- length(W$vals)
    nOmega <- (nQ ^ nI) * (nD ^ nJ) * (nW ^ nCO)
    adjw <- get_adjustment_weights(env)
  }
)

json <- jsonlite::toJSON(sapply(names(env), get, envir = env), pretty = TRUE)
write(json, file = "/Users/gva/drayage/TLPR/src/instance3x2_001.json")
