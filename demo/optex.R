# library(plot3D)
library(TLPR)

## Dynamic program
env <- new.env()
# generate_cssap(env, rate = 4.0, tau = 4L, nB = 5L, nCS =10L, nCO = 1L, nI = 1L, nJ = 1L)
# env$alpha <- c(rep(5.0, env$nI), rep(c(4.0, 8.0), each = env$nJ))
# env$R <- 10L

json_path <- "~/drayage/TLPR/src/instances/instance1x1_4_001.json"
jsonlite::fromJSON(json_path) |>
  list2env(envir = env)

# Scale up the value of the instance (experiment)
env$alpha <- 3.0 * env$alpha
# Transformation is needed
env$from_i <- list(c(env$from_i))
env$to_j   <- list(c(env$to_j))
# do.call(data.table::CJ, c(replicate(env$nJ, seq(nSJ), simplify = F), replicate(env$nI, seq(nSI), simplify = F)))[,c(2L,1L)]

model <- create_model(env)

dp_config <- function(env) {
  ## Endogenous state vector
  env$stateSupport <- 0L:env$R
  env$extendedStateSupport <- -env$R:env$R

  ## State indices
  env$Sdx <- do.call(TLPR::CartesianProductX, c(
    replicate(env$nI, seq_along(env$stateSupport), simplify = FALSE),
    replicate(env$nJ, seq_along(env$extendedStateSupport), simplify = FALSE)))

  env$nSdx <- nrow(env$Sdx)

  ## Decision vector
  env$actionSupport <- 0L:env$R
  env$Adx  <- seq_along(env$actionSupport)
  env$nAdx <- length(env$Adx)

  ## Exogenous state vector
  env$nQ <- length(env$Q$vals)
  env$nD <- length(env$D$vals)
  env$nW <- length(env$W$vals)

  env$scndx <- do.call(TLPR::CartesianProductX, c(
    replicate(env$nI,  seq(env$nQ), simplify = FALSE),
    replicate(env$nJ,  seq(env$nD), simplify = FALSE),
    replicate(env$nCO, seq(env$nW), simplify = FALSE)))

  env$scnpb <- apply(env$scndx, 1L, function(x) env$Q$prob[x[1L]] * env$D$prob[x[2L]] * env$W$prob[x[3L]])

  env$Qdx <- do.call(TLPR::CartesianProductX, replicate(env$nI,  seq(env$nQ), simplify = FALSE))
  env$Ddx <- do.call(TLPR::CartesianProductX, replicate(env$nJ,  seq(env$nD), simplify = FALSE))
  env$Wdx <- do.call(TLPR::CartesianProductX, replicate(env$nCO, seq(env$nW), simplify = FALSE))

  env$nQdx <- nrow(env$Qdx)
  env$nDdx <- nrow(env$Ddx)
  env$nWdx <- nrow(env$Wdx)

  with(env, nScen <- nQdx * nDdx * nWdx)
}

dp_config(env)

################################################################################################
## Environment
#############################################
n_chunks <- 12L
chunks   <- chunkup(env$nSdx, n_chunks)

library(foreach)
library(doParallel)

cores <- detectCores()
cl    <- makeCluster(cores-2L) # not to overload your computer
clusterExport(cl, c("env", "model"))
registerDoParallel(cl)

# Measure time lapse from here
timex <- Sys.time()

transit <- matrix(NA, nrow = env$tau * env$nSdx * env$nAdx * env$nScen, ncol = 6L)
for (t in seq(env$tau)) {
  print(t)
  foreach (idx = seq(n_chunks), .packages = c("gurobi","TLPR"), .combine = "rbind") %dopar% {
    start <- chunks[idx]
    end   <- chunks[idx+1L]-1L
    computeEnvironmentRx(env, model, t, start, end)
  } -> transit[(t - 1L) * env$nSdx * env$nAdx * env$nScen + seq(env$nSdx * env$nAdx * env$nScen),]
}

timex <- Sys.time() - timex

stopCluster(cl)
################################################################################################

## Optimize for specific scenario index
dynamic_programming <- function(env, transit, varphidx, ...) {
  # Initialize data containers
  V <- matrix(NA, nrow = env$tau+1L, ncol = env$nSdx)
  Q <- matrix(NA, nrow = env$tau, ncol = env$nSdx * env$nAdx)
  pi_rand <- matrix(NA, nrow = env$tau, ncol = env$nSdx * env$nAdx)
  pi_star <- matrix(NA, nrow = env$tau, ncol = env$nSdx * env$nAdx)

  # Terminal state values
  V[env$tau+1L,] <- - c(
    cbind(env$stateSupport[env$Sdx[,1L]], 
          pmax(env$extendedStateSupport[env$Sdx[,2L]], 0L), 
        - pmin(env$extendedStateSupport[env$Sdx[,2L]], 0L)) %*% c(env$alpha))

  ## Dynamic Programming Loop
  for (t in seq(env$tau,1L)) {
    print(t)
    k <- varphidx[t]
    for (i in seq(env$nSdx)) {
      for (j in seq(env$nAdx)) {
        p <- (((t - 1L) * env$nSdx + (i - 1L)) * env$nAdx + (j - 1L)) * env$nScen

        cost.t   <- transit[p + k, 2L]
        hold.t   <- h.t(env, env$stateSupport[env$Sdx[i, env$I_]], env$extendedStateSupport[env$Sdx[i, env$nI + env$J_]], env$alpha)
        reward.t <- - cost.t - hold.t

        # # Integrating stochasticity
        # next_p <- p + seq(env$nScen)
        # kdx    <- which(!is.na(transit[next_p, 1L]))
        # next_i <- transit[next_p, 1L][kdx]

        # prob   <- env$scnpb[kdx] / sum(env$scnpb[kdx])
        # Q[t, (i - 1L) * env$nAdx + j] <- reward.t + sum(V[t + 1L, next_i] * prob)

        # Deterministic solution (works)
        next_i <- transit[p + k, 1L]
        Q[t, (i - 1L) * env$nAdx + j] <- reward.t + V[t + 1L, next_i]
      }
      idx <- (i - 1L) * env$nAdx + env$Adx
      # Value function
      V[t,i] <- max(Q[t, idx], na.rm = TRUE)
      # Policy
      Qxs <- Q[t, idx] - min(Q[t, idx], na.rm = T) + 1.0
      # Stochastic optimal
      pi_rand[t, idx] <- data.table::fcoalesce(Qxs / sum(Qxs, na.rm = T), 0.0)
      # Deterministic optimal
      pi_star[t, idx] <- data.table::fcoalesce(as.numeric(Q[t, idx] == max(Q[t, idx], na.rm = T)), 0.0)
    }
  }

  list(
    "V" = V,
    "Q" = Q,
    "pi_star" = pi_star,
    "pi_rand" = pi_rand
  )
}

## Simulate policies
system_transition <- function(env, transit, pi, varphidx, init_s = NULL) {
  S.I  <- matrix(0.0, nrow = (env$tau + 1L) * env$nI, ncol = 1L)
  S.J  <- matrix(0.0, nrow = (env$tau + 1L) * env$nJ, ncol = 1L)
  X.I  <- matrix(0.0, nrow = env$tau * env$nI, ncol = 1L)
  X.J  <- matrix(0.0, nrow = env$tau * env$nJ, ncol = 1L)
  cost <- matrix(NA,  nrow = env$tau + 1L, ncol = 1L)
  q    <- numeric(env$tau)

  if (!is.null(init_s)) {
    S.I[1L,] <- init_s[env$I_]
    S.J[1L,] <- init_s[env$nI + env$J_]
  }

  # Simulation loop
  for (t in seq(env$tau)) {
    idx <- (t - 1L) * env$nI + env$I_
    jdx <- (t - 1L) * env$nJ + env$J_

    i <- sum(c(S.I[idx, 1L], env$R + S.J[jdx, 1L]) * env$stateKeys) + 1L
    j <- sample(seq(env$nAdx), 1L, prob = pi[t, (i - 1L) * env$nAdx + seq(env$nAdx)])
    k <- varphidx[t]

    p <- (((t - 1L) * env$nSdx + (i - 1L)) * env$nAdx + (j - 1L)) * env$nScen + k

    q[t] <- env$actionSupport[j]
    cost[t] <- transit[p, 2L] + h.t(env, S.I[idx, 1L], S.J[jdx, 1L], env$alpha)
    next_i  <- transit[p, 1L]

    S.I[env$nI + idx, 1L] <- env$stateSupport[env$Sdx[next_i, env$I_]]
    S.J[env$nJ + jdx, 1L] <- env$extendedStateSupport[env$Sdx[next_i, env$nI + env$J_]]

    X.I[idx, 1L] <- S.I[idx, 1L] + env$Q$vals[env$scndx[k, env$I_]] - S.I[idx + env$nI, 1L]
    X.J[jdx, 1L] <- S.J[jdx + env$nJ, 1L] - S.J[jdx, 1L] + env$D$vals[env$scndx[k, env$nI + env$J_]]
  }
  idx <- t * env$nI + env$I_
  jdx <- t * env$nJ + env$J_
  cost[t+1L] <- c(cbind(S.I[idx, 1L], pmax(S.J[jdx, 1L], 0L), - pmin(S.J[jdx, 1L], 0L)) %*% env$alpha)

  list(
    "S.I" = S.I,
    "S.J" = S.J,
    "X.I" = X.I,
    "X.J" = X.J,
    "cost" = sum(cost)
  )
}

# Simulate from transit 
# Fix scenario, e.g., c(17L, 15L, 17L, 14L), c(8L, 16L, 12L, 10L), c(5L,  8L,  8L, 11L), c(18L, 8L, 8L, 1L), c(9L, 18L, 7L, 10L)
(varphidx <- sample(nrow(env$scndx), env$tau, replace = TRUE, prob = env$scnpb))
for (k in varphidx) {
  print(c(env$Q$vals[env$scndx[k,env$I_]], env$D$vals[env$scndx[k,env$nI+env$J_]], env$W$vals[env$scndx[k,env$nI+env$J_+1L]]))
}

dynamic_programming(env, transit, varphidx) |>
  list2env(envir = .GlobalEnv)

S0  <- c(0L, 2L)
N   <- 5000L
sim <- numeric(N)
for (idx in 1L:N) {
  sim[idx] <- system_transition(env, transit, pi_rand, varphidx, init_s = S0)$cost
}
summary(sim)
system_transition(env, transit, pi_star, varphidx, init_s = S0)$cost

sum(sim < system_transition(env, transit, pi_star, varphidx, init_s = S0)$cost) / N

# Illustration of the system evolution
system_transition(env, transit, pi_star, varphidx, init_s = S0) |>
  list2env(envir = .GlobalEnv)

lanes <- c(env$L_, seq(env$nL))
dtx   <- data.table::as.data.table(cbind(
  t       = seq(env$tau),
  origin  = env$L[lanes,1L],
  destination = env$L[lanes,2L]+env$nI,
  assignment  = apply(cbind(X.I, X.J), 1L, min)))

coords <- list(
  "1x1" = list(
    "x" = c(-1.60, 1.60, -1.00, 1.00),
    "y" = c(-1.00,-1.00, -0.40,-0.40)
  ),
  "2x2" = list(
    "x" = c(-1.34,-1.34, 1.34, 1.34, -1.00,-1.00, 1.00, 1.00),
    "y" = c(-1.00, 1.00,-1.00, 1.00, -0.60, 0.60,-0.60, 0.60)
  ),
  "3x2" = list(
    "x" = c(-1.34,-1.34,-1.34, 1.34, 1.34,-1.00, -1.00,-0.70,  1.00, 1.00),
    "y" = c(-1.00, 0.00, 1.00,-1.00, 0.00,-0.60,  0.40, 1.20, -0.65, 0.35)
  )
)

setoo <- paste0(env$nI, "x", env$nJ)

par(mfrow = c(2L,3L), mai = c(0.2, 0.2, 0.3, 0.2))
for (t_ in seq(env$tau + 1L)) {
  # Create an igraph object
  g <- igraph::graph_from_data_frame(dtx[t==t_,-1L], directed = TRUE, vertices = seq(env$nI+env$nJ))
  # Set vertex and edge attributes
  igraph::V(g)$color <- "lightblue"
  igraph::V(g)$size <- 80L  # Set the size of the vertices
  igraph::V(g)$label.cex <- 1.25  # Increase the font size inside the nodes
  igraph::E(g)$arrow.size <- 0.5  # Set the size of the arrows
  # Set edge weights
  igraph::E(g)$weight <- dtx$assignment[(t_-1L)*env$nL+seq(env$nL)]
  # Plot the graph
  plot(
    g,
    edge.label = igraph::E(g)$weight,
    layout = igraph::layout_on_grid(g),
    xlim = c(-1.5, 1.5),
    edge.label.cex = 1.25)

  # Add the title closer to the graph
  title(main = paste0("t = ", t_), line = -4.0, cex.main = 1.5)  # Reduce space between title and graph

  k <- varphidx[t_]
  # Add labels on top left and top right
  text(
    x = coords[[setoo]]$x,
    y = coords[[setoo]]$y,
    labels = c(
      env$Q$vals[env$scndx[k,env$I_]],
      env$D$vals[env$scndx[k,env$nI+env$J_]],
      S.I[(t_-1L)*env$nI+env$I_,1L],
      S.J[(t_-1L)*env$nJ+env$J_,1L]),
    adj = c(.5, .5),
    cex = 1.25)
}
par(mfrow = c(1L,1L))

# Value function surface
t <- 1L

cbind2(TLPR::CartesianProductX(env$stateSupport, env$extendedStateSupport), V[t,]) |>
  as.data.frame() |>
  reshape2::dcast(V2 ~ V1, value.var = "V3") |>
  dplyr::select(-1L) |>
  as.matrix() |>
  `rownames<-`(seq(-env$R, env$R)) -> Vx

Vx |>
  heatmap(
    Rowv=NA, Colv=NA, 
    labRow = seq(-env$R,env$R), labCol = seq(0L,env$R), 
    col = cm.colors(256L), scale="none",
    xlab = "Entry State", ylab = "Exit State", main = paste0("State Value at t = ", t)
  )

plotVsurface(Vx, t)

# Action value surface for fixed exit state
exit_state <- 0L
Qdx <- get_Qdx_fixed_exit(env, exit_state)
# TLPR::CartesianProductX(env$actionSupport, env$stateSupport, env$extendedStateSupport)[Qdx,]

t <- 1L

cbind2(TLPR::CartesianProductX(env$actionSupport, env$stateSupport, env$extendedStateSupport)[Qdx,c(2L,1L)], Q[t, Qdx]) |>
  as.data.frame() |>
  reshape2::dcast(V2 ~ V1, value.var = "V3") |>
  dplyr::select(-1L) |>
  as.matrix() |>
  `rownames<-`(seq(0L, env$R)) -> Qx

Qx |>
  heatmap(Rowv=NA, Colv=NA, labRow = seq(0L,env$R), labCol = seq(0,env$R), col = cm.colors(256L), scale="none")

plotQsurfaceForFixedExit(Qx, exit_state, t)

# Action value surface for fixed entry state
entry_state <- 4L
Qdx <- get_Qdx_fixed_entry(env, entry_state)
# TLPR::CartesianProductX(env$actionSupport, env$stateSupport, env$extendedStateSupport)[Qdx,]

# matrix(Q[t, Qdx], ncol = env$nAdx)
cbind2(TLPR::CartesianProductX(env$actionSupport, env$stateSupport, env$extendedStateSupport)[Qdx,c(3L,1L)], Q[t, Qdx]) |>
  as.data.frame() |>
  reshape2::dcast(V2 ~ V1, value.var = "V3") |>
  dplyr::select(-1L) |>
  as.matrix() |>
  `rownames<-`(seq(0L, env$R)) -> Qx

Qx |>
  heatmap(Rowv=NA, Colv=NA, labRow = seq(0L,env$R), labCol = seq(-env$R,env$R), col = cm.colors(256L), scale="none")

plotQsurfaceForFixedEntry(Qx, entry_state, t, theta = 36.0, phi = 60.0)

################################################################################################
## Useful indexing checks
################################################################################################
# flowKeys <- c(nQ ^ (0L:(env$nI+env$nJ-1L)))
# env$flowKeys <- c(
#    env$nQ ^ (0L:(env$nI-1L)), 
#   (env$nQ ^ env$nI) * (env$nD ^ (0L:(env$nJ-1L))),
#   (env$nQ ^ env$nI) * (env$nD ^ (env$nJ)) * env$nW ^ (0L:(env$nCO-1L)))

# idx <- numeric(env$nQdx * env$nDdx * env$nWdx)
# for (k1 in seq(env$nQdx)) {
#   for (k2 in seq(env$nDdx)) {
#     for(k3 in seq(env$nWdx)) {
#       p <- ((k1 - 1L) * env$nDdx + (k2 - 1L)) * env$nWdx + k3
#       idx[p] <- sum(c(env$Qdx[k1,] - 1L, env$Ddx[k2,] - 1L, env$Wdx[k3,] - 1L) * env$flowKeys) + 1L
#     }
#   }
# }
# all(sort(idx) == seq(env$nQdx * env$nDdx * env$nWdx))

# env$stateKeys <- c(
#   length(env$stateSupport) ^ (0L:(env$nI-1L)), 
#   length(env$stateSupport) ^ env$nI * length(env$extendedStateSupport) ^ (0L:(env$nJ-1L)))
################################################################################################

#########################################################
# ## Check transition correctness manually
#########################################################
# idx <- 97325L
# 
# i <- transit[idx,3L]
# j <- transit[idx,4L]
# k <- transit[idx,5L]
# next_i <- transit[idx,1L]
# 
# env$stateSupport[env$Sdx[i,env$I_]]
# env$extendedStateSupport[env$Sdx[i,env$nI+env$J_]]
# 
# env$actionSupport[j]
# env$Q$vals[env$scndx[k,env$I_]]
# env$D$vals[env$scndx[k,env$nI+env$J_]]
# 
# env$stateSupport[env$Sdx[next_i,env$I_]]
# env$extendedStateSupport[env$Sdx[next_i,env$nI+env$J_]]
#########################################################

#######################################################################################
## OBSOLETE 
########################################################################################
# ## Dynamic program
# env <- new.env()
# # input_configuration(env)
# # simulate_auction(env)
# # initialize(env, rate = 0.4)
# generate_cssap(env, rate = 4.0, tau = 4L, nB = 5L, nCS =10L, nCO = 1L, nI = 1L, nJ = 1L)
# model <- create_model(env)
#
# # TLPR:::get_from_routes(env)
# # TLPR:::get_to_routes(env)
#
# max.S <- 10L
# env$R <- max.S
# inc.S <- 1L
#
# get_state_indices(.GlobalEnv, env$nI, env$nJ, max.S, inc.S)
#
# do.call(data.table::CJ, c(
#   replicate(env$nJ, seq(nSJ), simplify = F), 
#   replicate(env$nI, seq(nSI), simplify = F)))[,c(2L,1L)]
#
# ## Exogenous state vector
# Q <- list(
#   "vals" = c(0L, 4L, 8L),
#   "prob" = c(0.4, 0.3, 0.3)
# )
# D <- list(
#   "vals" = c(0L, 4L, 8L),
#   "prob" = c(0.25, 0.25, 0.5)
# )
# W <- list(
#   "vals" = c(7.0, 22.0),
#   "prob" = dnorm(c(7.0, 22.0), 18.0, 12.0) / sum(dnorm(c(7.0, 22.0), 18.0, 12.0))
# )
#
# # Homogeneity: all lanes same cost of transportation in the spot
# get_scenario_space(.GlobalEnv, env$nI, env$nJ, Q, D, W)
#
# ## Action space 
# # A_ <- seq(0L, env$nI*(max.S+max(Q$vals)), by = 5L)
# A_ <- seq(0L, 10L, by = 1L)
# nA <- length(A_)
#
# adj.w  <- c(1L, max.S+1L)
# K      <- seq(nOmega)
# alpha  <- c(rep(5.0, env$nI), rep(c(4.0, 8.0), each = env$nJ))
# resultog <- compute_environment(env, .GlobalEnv, 1L, 1L, nSdx, K, adj.w, FUN = optimal_assignment, alpha = alpha)
# # all(na.omit(resultog == result[seq(nSdx*nA*nOmega),]))
#
# # active <- which(complete.cases(result))
# # result[active,]
#
# # Serial execution
# # adj.w <- c(1L, max.S+1L)
# # K     <- seq(nOmega)
# # result <- matrix(NA, nrow = env$tau*nSdx*nA*nOmega, ncol = 5L)
# # for (t in seq(env$tau)) {
# #   print(t)
# #   p <- (t-1L)*nSdx*nA*nOmega+seq(nSdx*nA*nOmega)
# #   result[p,] <- compute_environment(env, t, 1L, nSdx, K, adj.w, alpha = alpha)
# # }
#
# ## Parallel run
# chunkup <- function(n, k) {
#   chunk_size <- n %/% k
#   chunks <- seq(1L, n, by = chunk_size)
#   chunks[k+1L] <- chunks[k+1L] + (n %% k)
#   chunks
# }
#
# n_chunks <- 12L
# chunks   <- chunkup(nSdx, n_chunks)
#
# library(foreach)
# library(doParallel)
#
# adj.w <- c(1L, max.S+1L)
# K     <- seq(nOmega)
# lenK  <- length(K)
#
# # adj.w <- c(1L, nSI, nSI^2L, (nSI^2L)*nSJ)
# # K <- which((rowSums(Phi.t[,1L:2L]) > 4L) & (rowSums(Phi.t[,3L:4L]) > 4L))
# # lenK  <- length(K)
#
# cores <- detectCores()
# cl    <- makeCluster(cores-4L) # not to overload your computer
# clusterExport(cl, c("model", "max.S", "SI_", "SJ_", "Sdx", "A_", "nA", "Phi.t"))
# registerDoParallel(cl)
#
# # result <- matrix(NA, nrow = env$tau*nSdx*nA*lenK, ncol = 5L)
# # foreach (t = seq(env$tau), .packages = c('gurobi','TLPR'), .combine = 'rbind') %dopar% {
# #   compute_environment(env, .GlobalEnv, t, 1L, nSdx, K, adj.w, alpha = alpha)
# # } -> result
#
# # Measure time lapse from here
# timex <- Sys.time()
#
# result <- matrix(NA, nrow = env$tau*nSdx*nA*lenK, ncol = 5L)
# for (t in seq(env$tau)) {
#   print(t)
#   foreach (idx = seq(n_chunks), .packages = c('gurobi','TLPR'), .combine = 'rbind') %dopar% {
#     start <- chunks[idx]
#     end   <- chunks[idx+1L]-1L
#     compute_environment(env, .GlobalEnv, t, start, end, K, adj.w, FUN = heuristic_assignment, alpha = alpha)
#   } -> result[(t-1L)*nSdx*nA*lenK+seq(nSdx*nA*lenK),]
# }
#
# timex <- Sys.time() - timex
#
# stopCluster(cl)
#
# # mget(ls(), envir = .GlobalEnv) |> saveRDS("~/Desktop/instance001.R")
#
# # active <- which(complete.cases(result))
# # is_active <- numeric(env$tau*nSdx*nA*lenK)
# # is_active[active] <- 1L
#
# result <- cbind2(result, NA)
# for (t in seq(env$tau)) {
#   result[(t-1L)*nSdx*nA*lenK+seq(nSdx*nA*lenK),6L] <- t
# }
# transit <- result
#
# tilde_phidx <- c(1L, 1L, 7L, 7L) # Fix scenario
#
# V <- matrix(NA, nrow = env$tau+1L, ncol = nSdx)
# V[env$tau+1L,] <- 0.0
# Q <- matrix(NA, nrow = env$tau, ncol = nSdx*nA)
# ## Dynamic Programming
# for (t in seq(env$tau,1L)) {
#   print(t)
#   for (i in seq(nSdx)) {
#     for (j in seq(nA)) {
#       position <- (((t-1L)*nSdx+(i-1L))*nA+(j-1L))*lenK
#       next_s <- transit[position + seq(lenK),1L]
#       cost.t <- transit[position + tilde_phidx[t],2L] # Fix scenario
#       scnidx <- transit[position + seq(lenK),5L]
#       prob   <- PPhi.t[scnidx] / sum(PPhi.t[K])
#       Q[t,(i-1L)*nA+j] <- - cost.t + c(V[t+1L,next_s] %*% prob)
#     }
#     V[t,i] <- max(Q[t,(i-1L)*nA+seq(nA)], na.rm = TRUE)
#   }
# }
#
# cbind2(cbind2(SI_[Sdx[,env$I_]], SJ_[Sdx[,env$nI+env$J_]]), V[2L,]) |>
#   as.data.frame() |>
#   reshape2::dcast(V2 ~ V1, value.var = "V3") |>
#   dplyr::select(-1L) |>
#   as.matrix() |>
#   `rownames<-`(seq(-env$R, env$R)) |>
#   heatmap(Rowv=NA, Colv=NA, labRow = seq(-env$R,env$R), labCol = seq(0L,env$R), col = cm.colors(256), scale="none")
#
# apply(t(matrix(Q[1L,], nrow = nA)), 1L, function(x) which.max(x))
#
# # Construct policy (time, state, action) to {0,1}
# pi_rand <- matrix(NA, nrow = env$tau, ncol = nSdx*nA)
# pi_star <- matrix(NA, nrow = env$tau, ncol = nSdx*nA)
# for (t in seq(env$tau)) {
#   for (i in seq(nSdx)) {
#       tmp <- (Q[t,(i-1L)*nA+seq(nA)] - min(Q[t,(i-1L)*nA+seq(nA)], na.rm = T)+1.0)
#       pi_rand[t,(i-1L)*nA+seq(nA)] <- data.table::fcoalesce(tmp / sum(tmp, na.rm = T), 0.0)
#       pi_star[t,(i-1L)*nA+seq(nA)] <- data.table::fcoalesce(
#         as.numeric(Q[t,(i-1L)*nA+seq(nA)] == max(Q[t,(i-1L)*nA+seq(nA)], na.rm = T)), 0.0)
#   }
# }
#
# ## Time benchmark
# t <- 1L
# i <- sample(nSdx, 1L)
# j <- sample(seq(nA), 1L)
# obj <- c(env$CTb, env$CTo[t,])
# rhs <- c(env$Cb[t,], env$Co[t,], env$R - SJ_[Sdx[i,env$nI+env$J_]], SI_[Sdx[i,env$I_]] + 
#            sample(c(0L, 4L, 8L), 1L, prob = c(0.2, 0.5, 0.3)), A_[j])
# microbenchmark::microbenchmark(
#   optimal_assignment(model, obj, rhs),
#   heuristic_assignment(model, obj, rhs, env$nvars, nrow(model$A)),
#   times = 100L
# )