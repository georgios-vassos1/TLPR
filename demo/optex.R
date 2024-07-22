library(TLPR)

## Dynamic program
env <- new.env()
# input_configuration(env)
# simulate_auction(env)
# initialize(env, rate = 0.4)
generate_cssap(env, rate = 4.0, tau = 4L, nB = 5L, nCS =10L, nCO = 1L, nI = 2L, nJ = 1L)
model <- create_model(env)

# TLPR:::get_from_routes(env)
# TLPR:::get_to_routes(env)

env$alpha <- c(rep(5.0, env$nI), rep(c(4.0, 8.0), each = env$nJ))

max.S <- 10L
env$R <- max.S
inc.S <- 1L

# do.call(data.table::CJ, c(replicate(env$nJ, seq(nSJ), simplify = F), replicate(env$nI, seq(nSI), simplify = F)))[,c(2L,1L)]

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
env$Adx <- seq_along(env$actionSupport)
env$nAdx <- length(env$Adx)

## Exogenous state vector
env$Q <- list(
  "vals" = c(0L, 4L, 8L),
  "prob" = c(0.4, 0.3, 0.3)
)
env$D <- list(
  "vals" = c(0L, 4L, 8L),
  "prob" = c(0.25, 0.25, 0.5)
)
env$W <- list(
  "vals" = c(7.0, 22.0),
  "prob" = dnorm(c(7.0, 22.0), 18.0, 12.0) / sum(dnorm(c(7.0, 22.0), 18.0, 12.0))
)

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

# flowKeys <- c(nQ ^ (0L:(env$nI+env$nJ-1L)))
env$flowKeys <- c(
   env$nQ ^ (0L:(env$nI-1L)), 
  (env$nQ ^ env$nI) * (env$nD ^ (0L:(env$nJ-1L))),
  (env$nQ ^ env$nI) * (env$nD ^ (env$nJ)) * env$nW ^ (0L:(env$nCO-1L)))

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

env$stateKeys <- c(
  length(env$stateSupport) ^ (0L:(env$nI-1L)), 
  length(env$stateSupport) ^ env$nI * length(env$extendedStateSupport) ^ (0L:(env$nJ-1L)))

with(env, nScen <- nQdx * nDdx * nWdx)

computeEnvironmentRx <- function(env, model, t, start, end) {

  transit <- matrix(NA, nrow = (end - start + 1L) * env$nAdx * env$nScen, ncol = 5L)
  for (i in seq(start, end)) {
    for (j in seq(env$nAdx)) {
      for (k1 in seq(env$nQdx)) {
        q <- env$Q$vals[env$Qdx[k1,]]

        rhs <- c(
          # Carrier capacity
          env$Cb[t,], env$Co[t,], 
          # Storage limits
          env$R - env$extendedStateSupport[env$Sdx[i,env$nI+env$J_]], env$stateSupport[env$Sdx[i,env$I_]] + q, 
          # Transport volume
          env$actionSupport[j])

        for (k3 in seq(env$nWdx)) {
          w <- rep(env$W$vals[env$Wdx[k3,]], env$nL)

          optx <- optimal_assignment(
            model, obj_ = c(env$CTb, w), 
            rhs_ = rhs)

          if (optx$status == "INFEASIBLE") next

          # Obtain origin-destination assignment volumes
          xI <- unlist(lapply(env$from_i, function(l) sum(optx$x[l])))
          xJ <- unlist(lapply(env$to_j,   function(l) sum(optx$x[l])))

          for (k2 in seq(env$nDdx)) {
            d <- env$D$vals[env$Ddx[k2,]]

            # Compute the index of the next state of the system
            next_i <- sum(
              c(
                pmax(pmin(env$stateSupport[env$Sdx[i,env$I_]] + q, env$R) - xI, 0L), 
                pmin(pmax(env$extendedStateSupport[env$Sdx[i,env$nI+env$J_]] - d, -env$R) + xJ, env$R) + env$R
              ) * env$stateKeys) + 1L

            kdx <- ((k3 - 1L) * env$nQdx + (k2 - 1L)) * env$nDdx + k1
            # kdx <- sum(c(Qdx[k1,] - 1L, Ddx[k2,] - 1L, Wdx[k3,] - 1L) * flowKeys) + 1L

            # Store into the transition matrix
            transit[((i - start) * env$nAdx + (j - 1L)) * env$nScen + kdx,] <- c(
              next_i,
              h.t(
                env, 
                env$stateSupport[env$Sdx[next_i,env$I_]], 
                env$extendedStateSupport[env$Sdx[next_i,env$nI+env$J_]], 
                alpha = env$alpha) + optx$objval, 
              i, j, kdx)
          }
        }
      }
    }
  }

  transit
}

## Parallel run

## May cause index out-of-bounds errors
#  chunkup <- function(n, k) {
#    chunk_size <- n %/% k
#    chunks <- seq(1L, n, by = chunk_size)
#    chunks[k+1L] <- chunks[k+1L] + (n %% k)
#    chunks
#  }

chunkup <- function(n, k) {
  chunk_size <- n %/% k
  chunks <- seq(1L, n, by = chunk_size)
  if (length(chunks) < k + 1) {
    chunks <- c(chunks, n + 1L)
  } else {
    chunks[k + 1] <- n + 1L
  }
  chunks
}

n_chunks <- 12L
chunks   <- chunkup(env$nSdx, n_chunks)

library(foreach)
library(doParallel)

cores <- detectCores()
cl    <- makeCluster(cores-4L) # not to overload your computer
clusterExport(cl, c("env", "model"))
registerDoParallel(cl)

# Measure time lapse from here
timex <- Sys.time()

result <- matrix(NA, nrow = env$tau * env$nSdx * env$nAdx * env$nScen, ncol = 5L)
for (t in seq(env$tau)) {
  print(t)
  foreach (idx = seq(n_chunks), .packages = c("gurobi","TLPR"), .combine = "rbind") %dopar% {
    start <- chunks[idx]
    end   <- chunks[idx+1L]-1L
    computeEnvironmentRx(env, model, t, start, end)
  } -> result[(t - 1L) * env$nSdx * env$nAdx * env$nScen + seq(env$nSdx * env$nAdx * env$nScen),]
}

timex <- Sys.time() - timex

stopCluster(cl)

result <- cbind2(result, NA)
for (t in seq(env$tau)) {
  result[(t - 1L) * env$nSdx * env$nAdx * env$nScen + seq(env$nSdx * env$nAdx * env$nScen), 6L] <- t
}
transit <- result

varphidx <- c(1L, 1L, 7L, 7L) # Fix scenario

V <- matrix(NA, nrow = env$tau+1L, ncol = env$nSdx)
V[env$tau+1L,] <- 0.0
Q <- matrix(NA, nrow = env$tau, ncol = env$nSdx * env$nAdx)
## Dynamic Programming
for (t in seq(env$tau,1L)) {
  print(t)
  for (i in seq(env$nSdx)) {
    for (j in seq(env$nAdx)) {
      p <- (((t - 1L) * env$nSdx + (i - 1L)) * env$nAdx + (j - 1L)) * env$nScen
      next_s <- transit[p + seq(env$nScen), 1L]
      cost.t <- transit[p + varphidx[t], 2L] # Fix scenario
      kdx    <- transit[p + seq(env$nScen), 5L]
      prob   <- env$scnpb[kdx] # / sum(env$scnpb[K])
      Q[t, (i - 1L) * env$nAdx + j] <- - cost.t + c(V[t + 1L, next_s] %*% prob)
    }
    V[t,i] <- max(Q[t, (i - 1L) * env$nAdx + seq(env$nAdx)], na.rm = TRUE)
  }
}

cbind2(cbind2(env$stateSupport[env$Sdx[,env$I_]], env$extendedStateSupport[env$Sdx[,env$nI+env$J_]]), V[2L,]) |>
  as.data.frame() |>
  reshape2::dcast(V2 ~ V1, value.var = "V3") |>
  dplyr::select(-1L) |>
  as.matrix() |>
  `rownames<-`(seq(-env$R, env$R)) |>
  heatmap(Rowv=NA, Colv=NA, labRow = seq(-env$R,env$R), labCol = seq(0L,env$R), col = cm.colors(256), scale="none")

apply(t(matrix(Q[1L,], nrow = env$nAdx)), 1L, function(x) which.max(x))


########################################################################################
## Dynamic program
env <- new.env()
# input_configuration(env)
# simulate_auction(env)
# initialize(env, rate = 0.4)
generate_cssap(env, rate = 4.0, tau = 4L, nB = 5L, nCS =10L, nCO = 1L, nI = 1L, nJ = 1L)
model <- create_model(env)

# TLPR:::get_from_routes(env)
# TLPR:::get_to_routes(env)

max.S <- 10L
env$R <- max.S
inc.S <- 1L

do.call(data.table::CJ, c(
  replicate(env$nJ, seq(nSJ), simplify = F), 
  replicate(env$nI, seq(nSI), simplify = F)))[,c(2L,1L)]

get_state_indices(.GlobalEnv, env$nI, env$nJ, max.S, inc.S)

## Exogenous state vector
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

# Homogeneity: all lanes same cost of transportation in the spot
get_scenario_space(.GlobalEnv, env$nI, env$nJ, Q, D, W)

## Action space 
# A_ <- seq(0L, env$nI*(max.S+max(Q$vals)), by = 5L)
A_ <- seq(0L, 10L, by = 1L)
nA <- length(A_)

adj.w  <- c(1L, max.S+1L)
K      <- seq(nOmega)
alpha  <- c(rep(5.0, env$nI), rep(c(4.0, 8.0), each = env$nJ))
resultog <- compute_environment(env, .GlobalEnv, 1L, 1L, nSdx, K, adj.w, FUN = optimal_assignment, alpha = alpha)
# na.omit(resultog == result)

# active <- which(complete.cases(result))
# result[active,]

# Serial execution
# adj.w <- c(1L, max.S+1L)
# K     <- seq(nOmega)
# result <- matrix(NA, nrow = env$tau*nSdx*nA*nOmega, ncol = 5L)
# for (t in seq(env$tau)) {
#   print(t)
#   p <- (t-1L)*nSdx*nA*nOmega+seq(nSdx*nA*nOmega)
#   result[p,] <- compute_environment(env, t, 1L, nSdx, K, adj.w, alpha = alpha)
# }

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

adj.w <- c(1L, max.S+1L)
K     <- seq(nOmega)
lenK  <- length(K)

# adj.w <- c(1L, nSI, nSI^2L, (nSI^2L)*nSJ)
# K <- which((rowSums(Phi.t[,1L:2L]) > 4L) & (rowSums(Phi.t[,3L:4L]) > 4L))
# lenK  <- length(K)

cores <- detectCores()
cl    <- makeCluster(cores-4L) # not to overload your computer
clusterExport(cl, c("model", "max.S", "SI_", "SJ_", "Sdx", "A_", "nA", "Phi.t"))
registerDoParallel(cl)

# result <- matrix(NA, nrow = env$tau*nSdx*nA*lenK, ncol = 5L)
# foreach (t = seq(env$tau), .packages = c('gurobi','TLPR'), .combine = 'rbind') %dopar% {
#   compute_environment(env, .GlobalEnv, t, 1L, nSdx, K, adj.w, alpha = alpha)
# } -> result

# Measure time lapse from here
timex <- Sys.time()

result <- matrix(NA, nrow = env$tau*nSdx*nA*lenK, ncol = 5L)
for (t in seq(env$tau)) {
  print(t)
  foreach (idx = seq(n_chunks), .packages = c('gurobi','TLPR'), .combine = 'rbind') %dopar% {
    start <- chunks[idx]
    end   <- chunks[idx+1L]-1L
    compute_environment(env, .GlobalEnv, t, start, end, K, adj.w, FUN = heuristic_assignment, alpha = alpha)
  } -> result[(t-1L)*nSdx*nA*lenK+seq(nSdx*nA*lenK),]
}

timex <- Sys.time() - timex

stopCluster(cl)

# mget(ls(), envir = .GlobalEnv) |> saveRDS("~/Desktop/instance001.R")

# active <- which(complete.cases(result))
# is_active <- numeric(env$tau*nSdx*nA*lenK)
# is_active[active] <- 1L

result <- cbind2(result, NA)
for (t in seq(env$tau)) {
  result[(t-1L)*nSdx*nA*lenK+seq(nSdx*nA*lenK),6L] <- t
}
transit <- result

tilde_phidx <- c(1L, 1L, 7L, 7L) # Fix scenario

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
      cost.t <- transit[position + tilde_phidx[t],2L] # Fix scenario
      scnidx <- transit[position + seq(lenK),5L]
      prob   <- PPhi.t[scnidx] / sum(PPhi.t[K])
      Q[t,(i-1L)*nA+j] <- - cost.t + c(V[t+1L,next_s] %*% prob)
    }
    V[t,i] <- max(Q[t,(i-1L)*nA+seq(nA)], na.rm = TRUE)
  }
}

cbind2(cbind2(SI_[Sdx[,env$I_]], SJ_[Sdx[,env$nI+env$J_]]), V[2L,]) |>
  as.data.frame() |>
  reshape2::dcast(V2 ~ V1, value.var = "V3") |>
  dplyr::select(-1L) |>
  as.matrix() |>
  `rownames<-`(seq(-env$R, env$R)) |>
  heatmap(Rowv=NA, Colv=NA, labRow = seq(-env$R,env$R), labCol = seq(0L,env$R), col = cm.colors(256), scale="none")

apply(t(matrix(Q[1L,], nrow = nA)), 1L, function(x) which.max(x))

# Construct policy (time, state, action) to {0,1}
pi_rand <- matrix(NA, nrow = env$tau, ncol = nSdx*nA)
pi_star <- matrix(NA, nrow = env$tau, ncol = nSdx*nA)
for (t in seq(env$tau)) {
  for (i in seq(nSdx)) {
      tmp <- (Q[t,(i-1L)*nA+seq(nA)] - min(Q[t,(i-1L)*nA+seq(nA)], na.rm = T)+1.0)
      pi_rand[t,(i-1L)*nA+seq(nA)] <- data.table::fcoalesce(tmp / sum(tmp, na.rm = T), 0.0)
      pi_star[t,(i-1L)*nA+seq(nA)] <- data.table::fcoalesce(
        as.numeric(Q[t,(i-1L)*nA+seq(nA)] == max(Q[t,(i-1L)*nA+seq(nA)], na.rm = T)), 0.0)
  }
}

## Time benchmark
t <- 1L
i <- sample(nSdx, 1L)
j <- sample(seq(nA), 1L)
obj <- c(env$CTb, env$CTo[t,])
rhs <- c(env$Cb[t,], env$Co[t,], env$R - SJ_[Sdx[i,env$nI+env$J_]], SI_[Sdx[i,env$I_]] + 
           sample(c(0L, 4L, 8L), 1L, prob = c(0.2, 0.5, 0.3)), A_[j])
microbenchmark::microbenchmark(
  optimal_assignment(model, obj, rhs),
  heuristic_assignment(model, obj, rhs, env$nvars, nrow(model$A)),
  times = 100L
)
