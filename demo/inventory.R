library(gurobi)

C.t <- function(w.t, x.t, z.t, params = list(OutputFlag = 0L), ...) {

  stopifnot((length(w.t) == length(x.t)) && (z.t <= sum(x.t)))

  n.t <- length(w.t)

  model <- list()

  model$A          <- matrix(1.0, ncol = n.t)
  model$obj        <- w.t
  model$modelsense <- 'min'
  model$lb         <- rep(0L, n.t)
  model$ub         <- x.t
  model$rhs        <- z.t
  model$vtype      <- rep('C', n.t)
  model$sense      <- c('=')

  result <- gurobi(model, params = params)
  return(
    list(
      cost = result$objval,
      allocation = result$x,
      status = result$status
    )
  )
}

## Exact dynamic programming
env <- new.env()
# Number of periods
env$tau   <- 4L
# Inventory limit
env$R <- 20L

## Endogenous state vector
env$stateSupport <- 0L:env$R
## State indices
env$Sdx  <- seq_along(env$stateSupport)
env$nSdx <- length(env$Sdx)

## Decision vector
env$actionSupport <- 0L:(env$R %/% 2L)
env$Adx  <- seq_along(env$actionSupport)
env$nAdx <- length(env$Adx)

## Exogenous state vector
discrete_disrt <- function(range, mu, spread) {
  prob <- dnorm(range, mu, spread)
  prob / sum(prob)
}

env$demandRange <- seq(0L, (env$R %/% 2L), by = 1L)

env$D <- list(
  "vals" = env$demandRange,
  "prob" = discrete_disrt(env$demandRange, 0.8 * env$R, env$R %/% 6L)
)

env$Ddx  <- seq_along(env$demandRange)
env$nDdx <- length(env$Ddx)

env$scndx <- seq(env$nDdx)
env$scnpb <- env$D$prob
# plot(env$demandRange, env$scnpb)

with(env, nScen <- nDdx)

## Problem design parameters
env$profit <- rnorm(env$tau, 100.0, 10.0)

# Number of sources
env$nSources <- 8L # Must be less than S.max
# Holding cost multiplier
env$alpha <- 30.0
# Execution cost
env$executionCost <- matrix(abs(rnorm(env$tau*env$nSources, 50.0, 20.0)), nrow = env$tau, ncol = env$nSources)
# Capacity vector
env$capacity <- matrix(rpois(env$tau*env$nSources, env$R %/% 4L), nrow = env$tau, ncol = env$nSources)

# Utility function
env$u <- matrix(NA, nrow = env$tau, ncol = env$nSdx * env$nDdx)
for (t in seq(env$tau)) {
  # for (i in env$Sdx) {
  #   for (k in env$Ddx) {
  #     env$u[t,(i - 1L) * env$nDdx + k] <- env$profit[t] * env$demandRange[k] - env$alpha * env$stateSupport[i]
  #   }
  # }
  env$u[t,] <- c(outer(env$profit[t] * env$demandRange, env$stateSupport * env$alpha, '-'))
}
matplot(env$stateSupport, t(matrix(env$u[1L,], nrow = env$nDdx)), type = 'l', ylab = 'Utility')

## Precompute immediate cost
env$cost <- matrix(NA, nrow = env$tau, ncol = env$nAdx)
for (t in seq(env$tau)) {
  print(t)
  for (j in env$Adx) {
    optx <- C.t(env$executionCost[t,], env$capacity[t,], env$actionSupport[j])
    env$cost[t,j] <- optx$cost
  }
}
matplot(env$actionSupport, t(env$cost), type = 'l', ylab = 'Cost')

# Dynamic programming
env$transit <- matrix(NA, nrow = env$tau * env$nSdx * env$nAdx * env$nDdx, ncol = 6L)
for (t in seq(env$tau)) {
  print(t)
  for (i in env$Sdx) {
    for (j in env$Adx) {
      for (k in env$Ddx) {
        # sdx <- pmin(pmax(env$stateSupport[i] + env$actionSupport[j] - env$demandRange[k], 0L), env$R) + 1L
        sdx <- env$stateSupport[i] + env$actionSupport[j] - env$demandRange[k] + 1L
        if (sdx < 1L || sdx > env$R + 1L) next
        idx <- (((t - 1L) * env$nSdx + (i - 1L)) * env$nAdx + (j - 1L)) * env$nDdx + k
        env$transit[idx, ] <- c(sdx, env$profit[t] * env$demandRange[k] - env$alpha * env$stateSupport[i] - env$cost[t, j], i, j, k, t)
      }
    }
  }
}
# all(env$transit[,1L] == c(pmin(pmax(outer(c(outer(-env$demandRange, env$actionSupport, '+')), env$stateSupport, '+'), 0L), env$R) + 1L), na.rm = T)

# Simulation
system_transition <- function(env, pi, varphidx, init_s = NULL) {
  S    <- matrix(0.0, nrow = env$tau + 1L, ncol = 1L)
  cost <- matrix(NA,  nrow = env$tau + 1L, ncol = 1L)
  q    <- numeric(env$tau)

  if (!is.null(init_s)) {
    S[1L,] <- init_s
  }

  # Simulation loop
  for (t in seq(env$tau)) {
    i <- S[t,] + 1L
    j <- sample(env$Adx, 1L, prob = pi[t, (i - 1L) * env$nAdx + env$Adx])
    k <- varphidx[t]

    p <- (((t - 1L) * env$nSdx + (i - 1L)) * env$nAdx + (j - 1L)) * env$nScen + k

    q[t]    <- env$actionSupport[j]
    cost[t] <- env$transit[p, 2L] # + env$alpha * env$stateSupport[i]
    next_i  <- env$transit[p, 1L]

    S[t + 1L, 1L] <- env$stateSupport[env$Sdx[next_i]]
  }
  cost[t+1L] <- S[t + 1L, 1L] * env$alpha

  list(
    "S" = S,
    "A" = q,
    "D" = env$demandRange[varphidx],
    "cost" = sum(cost)
  )
}

# Sample scenarios
(varphidx <- sample(env$Ddx, env$tau, replace = TRUE)) # , prob = env$D$prob))

V <- matrix(NA, nrow = env$tau + 1L, ncol = env$nSdx)
V[env$tau+1L,] <- env$stateSupport * (env$alpha %/% 1L)
Q <- matrix(0.0, nrow = env$tau, ncol = env$nSdx * env$nAdx)
pi_rand <- matrix(NA, nrow = env$tau, ncol = env$nSdx * env$nAdx)
pi_star <- matrix(NA, nrow = env$tau, ncol = env$nSdx * env$nAdx)
# Exact dynamic programming
for (t in seq(env$tau,1L)) {
  print(t)
  k <- varphidx[t]
  for (i in env$Sdx) {
    p <- (((t - 1L) * env$nSdx + (i - 1L)) * env$nAdx + (env$Adx - 1L)) * env$nDdx + k

    actions <- env$transit[p, 4L] |> na.omit() |> as.vector()
    rewards <- env$transit[p, 2L] |> na.omit() |> as.vector()
    # next_sdx <- env$transit[p, 1L] |> na.omit() |> as.vector()

    for (jdx in seq_along(actions)) {
      # next_i <- next_sdx[jdx]
      # Q[t, (i - 1L) * env$nAdx + actions[jdx]] <- rewards[jdx] + V[t + 1L, next_i]

      next_p <- (((t - 1L) * env$nSdx + (i - 1L)) * env$nAdx + (actions[jdx] - 1L)) * env$nDdx + env$Ddx
      kdx    <- which(!is.na(env$transit[next_p, 1L]))
      next_i <- env$transit[next_p, 1L][kdx]

      prob <- env$D$prob[kdx] / sum(env$D$prob[kdx])
      Q[t, (i - 1L) * env$nAdx + actions[jdx]] <- rewards[jdx] + sum(prob * V[t + 1L, next_i])
    }
    idx <- (i - 1L) * env$nAdx + env$Adx
    # Value function
    V[t,i] <- max(Q[t, idx], na.rm = TRUE)
    # Policy
    Qxs <- Q[t, idx]
    # Stochastic optimal
    pi_rand[t, idx] <- data.table::fcoalesce(Qxs / sum(Qxs, na.rm = T), 0.0)
    # Deterministic optimal
    pi_star[t, idx] <- data.table::fcoalesce(as.numeric(Q[t, idx] == max(Q[t, idx], na.rm = T)), 0.0)
  }
}
t <- 1L
plot(V[t,])
# i <- 1L
# plot(env$actionSupport, Q[t, (i - 1L) * env$nAdx + env$Adx], type = 's')

S0 <- 0L
N   <- 10000L
sim <- numeric(N)
for (idx in 1L:N) {
  sim[idx] <- system_transition(env, pi_rand, varphidx, init_s = S0)$cost
}
summary(sim)
system_transition(env, pi_star, varphidx, init_s = S0)

sum(sim > system_transition(env, pi_star, varphidx, init_s = S0)$cost) / N

# # Value function (exhaustive)
# V <- matrix(NA, nrow = env$tau + 1L, ncol = env$nSdx * env$nDdx)
# V[env$tau+1L,] <- rep(env$stateSupport * (env$alpha %/% 10L), each = env$nScen)
# Q <- matrix(0.0, nrow = env$tau * env$nSdx * env$nDdx, ncol = env$nAdx)
# # Exact dynamic programming
# for (t in seq(env$tau,1L)) {
#   print(t)
#   for (i in env$Sdx) {
#     for (k in env$Ddx) {
#       qdx <- ((t - 1L) * env$nSdx + (i - 1L)) * env$nDdx + k
#       vdx <- (i - 1L) * env$nDdx + k
# 
#       p <- (((t - 1L) * env$nSdx + (i - 1L)) * env$nAdx + env$Adx - 1L) * env$nDdx + k
# 
#       actions <- env$transit[p, 4L] |> na.omit() |> as.vector()
#       next_i  <- env$transit[p, 1L] |> na.omit() |> as.vector()
#       rewards <- env$transit[p, 2L] |> na.omit() |> as.vector()
# 
#       Q[qdx, actions] <- rewards + c(env$D$prob %*% matrix(V[t + 1L, c(outer(env$Ddx, (next_i - 1L) * env$nDdx, '+'))], nrow = env$nDdx))
#       V[t, vdx] <- max(Q[qdx,], na.rm = T)
#     }
#   }
# }
# t <- 2L
# i <- 5L
# plot(env$stateSupport, c(env$D$prob %*% matrix(V[t,], nrow = env$nDdx)), type = 'o', ylab = 'V(s)', xlab = 's')
# matplot(env$stateSupport, t(matrix(V[t,], nrow = env$nDdx)), type = 'l', ylab = 'V(s)', xlab = 's')
# matplot(env$actionSupport, t(Q[((t - 1L) * env$nSdx + (i - 1L)) * env$nDdx + env$Ddx,]), ylab = 'Q(s,x)', xlab = 'x', type = 's')
