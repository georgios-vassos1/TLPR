library(TLPR)

env <- new.env()
generate_cssap(env, tau = 12L)
model <- create_model(env)

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
  k = env$n_pairs
)

get_from_routes(env)
get_to_routes(env)

N <- 30L
costs <- matrix(NA, nrow = env$tau, ncol = npi*N)
i <- 1L
while (i <= N) {
  print(i)
  exog <- list(
    "Q" = list("func"=rpois, "params"=list('n'=env$tau*env$nI, 'lambda'=40L)),
    "D" = list("func"=rpois, "params"=list('n'=env$tau*env$nI, 'lambda'=40L))
  )
  result <- simulate_system(env, policy, args, exog)
  if (result$status) {
    costs[,(i-1L)*npi+seq(npi)] <- result$cost
    i <- i + 1L
  }
}

stats_ <- do.call(cbind, lapply(seq(npi), function(x) compute_stats(costs[,seq(x, npi*N, by = npi)], N)))
matplot(stats_, type = 'l', lwd = rep(c(2L,1L,1L), 2L), lty = 1L, pch = 19L, ylab = "Cumulative Cost", col = rep(c(1L,2L), each = 3L))

## Single run plots
matplot(apply(result$cost, 2L, cumsum), type = 'b', lwd = 2L, lty = 1L, pch = 19L, ylab = "Cumulative Cost")

pdx <- 1L
idx <- outer((seq(env$tau)-1L)*env$nI, seq(env$nI), '+')
i <- 1L
plot(result$S.I[idx[,i],pdx], type = 's', ylim = c(min(result$S.I[,pdx],result$X.I[,pdx]), max(result$S.I[,pdx],result$X.I[,pdx])))
lines(result$X.I[idx[,i],pdx], type = 's', col = 2L)
lines(result$Q[idx[,i]], type = 's', col = 3L)
