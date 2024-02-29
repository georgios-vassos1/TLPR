library(TLPR)

## Dynamic program
env <- new.env()
input_configuration(env)
simulate_auction(env)
initialize(env, rate = 0.4)
model <- create_model(env)

get_from_routes(env)
get_to_routes(env)

max.S <- env$R
inc.S <- 1L

get_state_indices(env, max.S, inc.S)

## Exogenous state vector
Q <- list(
  "vals" = seq(0L, env$R, by = 2L),
  "prob" = c(0.2, 0.4, 0.4)
)
D <- list(
  "vals" = seq(0L, env$R, by = 2L),
  "prob" =c(0.2, 0.4, 0.4, 0.15)
)
W <- list(
  "vals" = c(22.-8., 22., 22.+8.),
  "prob" = dnorm(c(22.-8., 22., 22.+8.), 22.0, 8.0) / sum(dnorm(c(22.-8., 22., 22.+8.), 22.0, 8.0))
)

# Homogeneity: all lanes same cost of transportation in the spot
get_scenario_space(Q, D, W)

## Action space
A_ <- seq(0L, env$nI*(max.S+max(Q$vals)), by = 4L)
nA <- length(A_)

## Parallel run
chunkup <- function(n, k) {
  chunk_size <- n %/% k
  chunks <- seq(1L, n, by = chunk_size)
  chunks[k+1L] <- chunks[k+1L] + (n %% k)
  chunks
}

n_chunks <- 48L
chunks   <- chunkup(nSdx, n_chunks)

library(foreach)
library(doParallel)

cores <- detectCores()
cl    <- makeCluster(cores-4L) # not to overload your computer
clusterExport(cl, c("model", "max.S", "SI_", "SJ_", "Sdx", "A_", "nA", "Phi.t"))
registerDoParallel(cl)

adj.w <- c(1L, nSI, nSI^2L, (nSI^2L)*nSJ)

t <- 1L
K <- c(150L,151L)
foreach (idx = seq(n_chunks), .packages = c('gurobi','TLPR'), .combine = 'rbind') %dopar% {
  start <- chunks[idx]
  end   <- chunks[idx+1L]-1L
  run_scenarios(env, t, start, end, K, adj.w)
} -> result

stopCluster(cl)

## Dynamic Programming
for (t in seq(tau,1L)) {
  for (i in seq(nSdx)) {
    for (j in seq(nA)) {
      next_i <- transit[(i-1L)*nOmega*nA+(j-1L)*nOmega+seq(nOmega),]
      Q[t,(i-1L)*nA+j] <- u[t,(i-1L)*nA+j] + V[t+1L,(next_i-1L)*nA+j] %*% PPhi.t
    }
    V[t,i] <- max(Q[t,(i-1L)*nA+seq(nA)])
  }
}
