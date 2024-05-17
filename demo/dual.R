library(TLPR)

## Dynamic program
env <- new.env()
# input_configuration(env)
# simulate_auction(env)
# initialize(env, rate = 0.4)
generate_cssap(env, rate = 4.0, tau = 4L, nCS = 1L, nCO = 1L, nI = 2L, nJ = 2L)
model <- create_model(env)
model$vtype <- rep('C', env$nvars)

env$from_i <- vector(mode = "list", length = env$nI)
for (i in seq(env$nI)) {
  idx <- (seq(env$nJ) - 1L) * env$nI
  msk <- which(apply(outer(env$L_, idx + i, "=="), 1L, any))
  idx <- env$nL_ + c(outer(idx, (seq(env$nCO) - 1L) * env$nL, "+")) + i
  env$from_i[[i]] <- c(msk, idx)
}

env$to_j <- vector(mode = "list", length = env$nJ)
for (j in seq(env$nJ)) {
  idx <- (j - 1L) * env$nI + seq(env$nI)
  msk <- which(apply(outer(env$L_, idx, "=="), 1L, any))
  idx <- env$nL_ + c(outer(idx, (seq(env$nCO) - 1L) * env$nL, 
                           "+"))
  env$to_j[[j]] <- c(msk, idx)
}

get_from_routes(env)
get_to_routes(env)

max.S <- 10L
env$R <- max.S
inc.S <- 1L

get_state_indices(env, max.S, inc.S)

## Exogenous state vector
Q <- list(
  "vals" = c(0L, 4L, 8L),
  "prob" = c(0.2, 0.5, 0.3)
)
D <- list(
  "vals" = c(0L, 4L, 8L),
  "prob" = c(0.25, 0.25, 0.5)
)
W <- list(
  "vals" = 15.0,
  "prob" = 1.0
)

# Homogeneity: all lanes same cost of transportation in the spot
get_scenario_space(Q, D, W)

## Action space 
# A_ <- seq(0L, env$nI*(max.S+max(Q$vals)), by = 5L)
A_ <- seq(0L, 15L, by = 5L)
nA <- length(A_)

adj.w  <- c(1L, max.S+1L, (max.S+1L)^2.)
alpha  <- c(rep(2.5, env$nI), rep(c(2.5, 4.5), each = env$nJ))

cbind(SI_[Sdx[,1L]], SI_[Sdx[,2L]], max.S+SJ_[Sdx[,env$nI+env$J_]]) %*% adj.w + 1L

# Serial execution
scenaria <- 8L
nScen <- length(scenaria)
transit <- matrix(NA, nrow = env$tau*nSdx*nA*nScen, ncol = 5L)

stateIdx <- function(env, Si, Sj, max.S, weights, ...) {
  c(c(Si, Sj+max.S) %*% weights + 1L)
}

library(data.table)
cols <- c("t", "state", "action", "scenario", "objval", "pi.1", "pi.2", "pi.3", "pi.4", "pi.5")
# Initialize the data.table with 10 rows of NA values
dt <- data.table(matrix(NA, nrow = env$tau * nSdx * nA * nScen, ncol = length(cols), dimnames = list(NULL, cols)))
dt[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]

dt$t        <- rep(seq(env$tau), each = nSdx*nA*nScen)
dt$state    <- rep(seq(nSdx), each = nA*nScen, times = env$tau)
dt$action   <- A_[rep(rep(seq(nA), each = nScen, times = env$tau*nSdx))]
dt$scenario <- rep(scenaria, each = env$tau*nSdx*nA)
dt$Cb       <- rep(env$Cb, each = nSdx*nA*nScen)
dt$Co       <- rep(env$Co, each = nSdx*nA*nScen)

setkeyv(dt, c("t", "state", "action", "scenario"))

for (t_ in seq(env$tau)) {
  for (j in seq(2L, nA)) {
    for (i in seq(nSdx)) {
      # if(!(i%%100L)) print(i)
      for (k in seq(nScen)) {
        scndx <- scenaria[k]
        q <- Phi.t[scndx,env$I_]
        d <- Phi.t[scndx,env$nI+env$J_]
        w <- rep(Phi.t[scndx,env$nI+env$nJ+1L], env$nL)

        optx <- optimal_assignment(
          model, obj_ = c(env$CTb, w), 
          rhs_ = c(env$Cb[t_,], env$Co[t_,], env$R - SJ_[Sdx[i,env$nI+env$J_]], SI_[Sdx[i,env$I_]] + q, A_[j]))

        if (optx$status == "INFEASIBLE") next
        print(c(optx$objval, optx$x, optx$pi))
        idx <- (((t_-1L)*nSdx*nA*nScen+(i-1L)*nA*nScen+(j-1L)*nScen+k))
        # dt[idx, c("objval", "pi.1", "pi.2", "pi.3", "pi.4", "pi.5") := as.list(c(optx$objval, optx$pi))]

        xI <- unlist(lapply(env$from_i,   function(l) sum(optx$x[l])))
        xJ <- unlist(lapply(env$to_j, function(l) sum(optx$x[l])))

        next_i <- stateIdx(
          env, 
          pmax(pmin(SI_[Sdx[i,env$I_]] + q, env$R) - xI, 0L), 
          pmin(pmax(SJ_[Sdx[i,env$nI+env$J_]] - d, -env$R) + xJ, env$R), 
          max.S, adj.w)
        transit[(((t_-1L)*nSdx+(i-1L))*nA+(j-1L))*nScen+k,] <- c(
          next_i,
          h.t(env, SI_[Sdx[next_i,env$I_]], SJ_[Sdx[next_i,env$nI+env$J_]], alpha = alpha) + optx$objval, 
          i, j, scndx)
      }
    }
  }
}

cbind2(cbind2(SI_[Sdx[,env$I_]], SJ_[Sdx[,env$nI+env$J_]]), reshape2::dcast(dt[t==1L,], state ~ action, value.var = "pi.1")[,-1L]) |>
  setNames(c("SI", "SJ", paste0("A", seq(nA)))) |>
  (\(dtx) lapply(seq(nA), function(x) reshape2::dcast(dtx, SI ~ SJ, value.var = paste0('A', x))[,-1L] |> as.matrix() |> unname()))() -> lst

lst

lst[[2L]] |>
  heatmap(Rowv=NA, Colv=NA, labCol = seq(-env$R,env$R), labRow = seq(0L,env$R), col = cm.colors(256L), scale="none")

for (state in SI_) {
  for (q in Q$vals) {
    for (action in seq(0L, max(A_))) {
      cat(" ")
      cat((action - state - q > 0.0) * 1L)
    }
    cat("\n")
  }
}
flags <- c((outer(seq(0L, max(A_)), c(outer(Q$vals, SI_, '+')), '-') <= 0.0) * 1L)
nAdx <- max(A_) + 1L
t(((outer(seq(0L, max(A_)), c(outer(Q$vals, SI_, '+')), '-') <= 0.0) * 1L))
for (i in seq(nSI)) {
  for (j in seq_along(Q$vals)) {
    idx <- (i-1L)*length(Q$vals)*nAdx + (j-1L)*nAdx + seq(nAdx)
    print(flags[idx])
  }
}
