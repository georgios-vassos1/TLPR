library(TLPR)

## Dynamic program
env <- new.env()
# input_configuration(env)
# simulate_auction(env)
# initialize(env, rate = 0.4)
generate_cssap(env, rate = 10.0, tau = 4L, nCS = 1L, nCO = 1L, nI = 1L, nJ = 1L)
model <- create_model(env)

get_from_routes(env)
get_to_routes(env)

max.S <- env$R
inc.S <- 1L

get_state_indices(env, max.S, inc.S)

## Exogenous state vector
Q <- list(
  "vals" = seq(0L, env$R %/% env$tau, by = 10L),
  "prob" = c(0.2, 0.4, 0.4)
)
D <- list(
  "vals" = seq(0L, env$R %/% env$tau, by = 10L),
  "prob" =c(0.2, 0.4, 0.4, 0.15)
)
W <- list(
  "vals" = c(22.-8., 22., 22.+8.),
  "prob" = dnorm(c(22.-8., 22., 22.+8.), 22.0, 8.0) / sum(dnorm(c(22.-8., 22., 22.+8.), 22.0, 8.0))
)

# Homogeneity: all lanes same cost of transportation in the spot
get_scenario_space(Q, D, W)

## Action space 
# A_ <- seq(0L, env$nI*(max.S+max(Q$vals)), by = 5L)
A_ <- seq(0L, max(env$Cb+env$Co), by = 5L)
nA <- length(A_)

run_scenarios(env, 1L, 1L, nSdx, 1L, adj.w)

adj.w <- c(1L, 1L)
K     <- seq(nOmega)
result <- matrix(NA, nrow = env$tau*nSdx*nA*nOmega, ncol = 5L)
for (t in seq(env$tau)) {
  print(t)
  result[(t-1L)*nSdx*nA*nOmega+seq(nSdx*nA*nOmega),] <- run_scenarios(env, t, 1L, nSdx, K, adj.w)
}

## Parallel run
chunkup <- function(n, k) {
  chunk_size <- n %/% k
  chunks <- seq(1L, n, by = chunk_size)
  chunks[k+1L] <- chunks[k+1L] + (n %% k)
  chunks
}

n_chunks <- 12L
chunks   <- chunkup(nSdx, n_chunks)

library(foreach)
library(doParallel)

# adj.w <- c(1L, nSI, nSI^2L, (nSI^2L)*nSJ)
# K <- which((rowSums(Phi.t[,1L:2L]) > 4L) & (rowSums(Phi.t[,3L:4L]) > 4L))

adj.w <- c(1L, 1L)
K     <- seq(nOmega)
lenK  <- length(K)

cores <- detectCores()
cl    <- makeCluster(cores-4L) # not to overload your computer
clusterExport(cl, c("model", "max.S", "SI_", "SJ_", "Sdx", "A_", "nA", "Phi.t"))
registerDoParallel(cl)

result <- matrix(NA, nrow = env$tau*nSdx*nA*lenK, ncol = 5L)
for (t in seq(env$tau)) {
  print(t)
  foreach (idx = seq(n_chunks), .packages = c('gurobi','TLPR'), .combine = 'rbind') %dopar% {
    start <- chunks[idx]
    end   <- chunks[idx+1L]-1L
    run_scenarios(env, t, start, end, K, adj.w)
  } -> result[(t-1L)*nSdx*nA*lenK+seq(nSdx*nA*lenK),]
}

stopCluster(cl)

# mget(ls(), envir = .GlobalEnv) |> saveRDS("~/Desktop/instance001.R")

active <- which(complete.cases(result))

result <- cbind2(result, NA)
for (t in seq(env$tau)) {
  result[(t-1L)*nSdx*nA*lenK+seq(nSdx*nA*lenK),6L] <- t
}
transit <- result

tilde_phidx <- 12L # Fix scenario

V <- matrix(NA, nrow = env$tau+1L, ncol = nSdx)
V[env$tau+1L,] <- 0.0
Q <- matrix(NA, nrow = env$tau, ncol = nSdx*nA)
## Dynamic Programming
for (t in seq(env$tau,1L)) {
  print(t)
  for (i in seq(nSdx)) {
    for (j in seq(nA)) {
      position <- (((t-1L)*nSdx+(i-1L))*nA+(j-1L))*lenK
      next_s <- transit[position + seq(lenK),1L]
      cost.t <- transit[position + tilde_phidx,2L] # Fix scenario
      scnidx <- transit[position + seq(lenK),5L]
      prob   <- PPhi.t[scnidx] / sum(PPhi.t[K])
      Q[t,(i-1L)*nA+j] <- - cost.t + c(V[t+1L,next_s] %*% prob)
    }
    V[t,i] <- max(Q[t,(i-1L)*nA+seq(nA)], na.rm = TRUE)
  }
}

sum(apply(t(matrix(Q[1L,], nrow = nA)), 1L, function(x) which.max(x)) != 1L)

      