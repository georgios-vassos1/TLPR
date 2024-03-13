library(TLPR)

env <- new.env()
generate_cssap(env, rate = 10.0, tau = 12L, nCS = 5L, nCO = 1L, nB = 3L, nI = 4L, nJ = 3L)

env$alpha <- c(rep(1.0, env$nI), rep(c(1.0, 2.0), each = env$nJ))

# Choose the effective constraints
#   Constraints 2L and 3L depend on the state of the system.
#   This entails significant computational burden in evaluating the portfolio contract.
Idx    <- c(1L,4L)
rhs_dx <- unlist(list(
  seq(env$nCS+env$nCO), 
  env$nCS+env$nCO+seq(env$nI),
  env$nCS+env$nCO+env$nI+seq(env$nJ),
  env$nCS+env$nCO+env$nI+env$nJ+1L)[Idx])

model <- create_model(env, Idx)

policy <- list(
  "myopic" = optimal_assignment,
  "random" = random_assignment
)
npi <- length(policy)

args <- list(
  model = model,
  obj_ = NULL,
  rhs_ = NULL,
  n = NULL,
  k = env$nvars,
  x = NULL
)

exog <- list(
  "Q" = list("func"=rpois, "params"=list('n'=env$tau*env$nI, 'lambda'=10L)),
  "D" = list("func"=rpois, "params"=list('n'=env$tau*env$nJ, 'lambda'= 5L))
)

get_from_routes(env)
get_to_routes(env)

N <- 50L
costs <- matrix(NA, nrow = env$tau, ncol = npi*N)
i <- 1L
while (i <= N) {
  print(i)
  result <- simulate_system(env, policy, args, exog, rhs_dx, correction = TRUE)
  if (result$status) {
    costs[,(i-1L)*npi+seq(npi)] <- result$cost
    i <- i + 1L
  }
}

stats_ <- do.call(cbind, lapply(seq(npi), function(x) compute_stats(costs[,seq(x, npi*N, by = npi)], N = 1L)))
matplot(stats_, type = 'l', lwd = rep(c(2L,1L,1L), 2L), lty = 1L, pch = 19L, ylab = "Cumulative Cost", col = rep(c(1L,2L,3L), each = 3L))

## Single run plots
matplot(apply(result$cost, 2L, cumsum), type = 'b', lwd = 2L, lty = 1L, pch = 19L, ylab = "Cumulative Cost")

pdx <- 1L
idx <- outer((seq(env$tau)-1L)*env$nI, seq(env$nI), '+')
i <- 2L
plot(result$S.I[idx[,i],pdx], type = 's', ylim = c(min(result$S.I[,pdx],result$X.I[,pdx]), max(result$S.I[,pdx],result$X.I[,pdx])))
lines(result$X.I[idx[,i],pdx], type = 's', col = 2L)
lines(result$Q[idx[,i]], type = 's', col = 3L)
