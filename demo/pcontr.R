library(TLPR)

## 
evaluate_portfolio <- function(env, A_, plot = FALSE) {

  cost <- rep(NaN, env$tau * nA)
  allocation <- rep(NaN, env$tau * nA * env$nvars)
  nA <- length(A_)

  for (t in seq(env$tau)) {
    obj_ <- c(env$CTb, env$CTo[t,])
    rhs_ <- c(env$Cb[t,], env$Co[t,], NA)
    for (j in seq(nA)) {
      rhs_[env$nCS+env$nCO+1L] <- A_[j]
      optx <- optimal_assignment(model, obj_, rhs_)
      if (optx$status == "OPTIMAL") {
        cost[(t-1L)*nA+j] <- optx$objval
        allocation[((t-1L)*nA+(j-1L))*env$nvars+seq(env$nvars)] <- optx$x
      } else {
        print(optx$status)
      }
    }
  }

  if (plot) {
    matplot(A_, matrix(cost, nrow = nA), type = 'l')
  }

  list(
    "cost" = cost,
    "assignment" = allocation
  )
}

library(data.table)
library(igraph)

portfolio_evaluation_data <- function(env, A_, assignments) {
  ## Auxiliary data
  nA <- length(A_)
  lanes    <- c(env$L_, seq(env$nL))
  carriers <- rep(seq_along(env$winner), env$nLc[-1L])
  rbind2(
    cbind2(unlist(carriers), unlist(env$B)),
    cbind2(rep(env$nCS+seq(env$nCO), each = env$nL), seq(env$nL))
  ) -> carriers
  # carriers[,1L] <- rep(seq_along(rle(carriers[,1L])$values), rle(carriers[,1L])$lengths)
  # Create the data table.
  dt <- as.data.table(cbind(
    t       = rep(seq(env$tau), each = nA * env$nvars),
    action  = rep(rep(A_, each = env$nvars), env$tau),
    carrier = carriers[,1L],
    lane    = lanes,
    origin  = env$L[lanes,2L],
    destination = env$L[lanes,1L]+env$nI,
    assignment  = assignments))
  # Add capacity and cost information.
  dt[, capacity := cbind2(env$Cb,env$Co)[t,carrier], by = .(t, carrier)]
  dt[, cost := c(env$CTb,env$CTo[t,]), by = .(t,action)]
  # Assign to global context.
  assign("dt", dt, envir = globalenv())
  cat("A data.table called `dt' has been stored in the global context.")
}

## Demo
# Initiate instance
env <- new.env()
generate_cssap(env, rate = 10.0, tau = 12L, nCS = 10L, nCO = 1L, nB = 3L, nI = 2L, nJ = 2L)
# Generate optimization model
Idx   <- c(1L,4L)
model <- create_model(env, Idx)
# Action space 
A_ <- seq(0L, max(sweep(env$Cb, 1L, env$Co, '+')), by = 5L)
nA <- length(A_)
cpx <- evaluate_portfolio(env, A_, plot = TRUE)

portfolio_evaluation_data(env, A_, cpx$assignment)

## Plot
# action-lane indexes for a single block
idx <- outer((seq(nA)-1L)*env$nvars, seq(env$nvars), '+')
# lane indexes for a given action over time
j   <- nA
jdx <- c(outer(idx[j,], (seq(env$tau)-1L)*nA*env$nvars, '+'))
# Filter the data table
dtx <- dt[jdx, .(assignment = sum(assignment, na.rm = T)), by = .(t, origin, destination)]

n <- 2L
reflex <- diag(n)[,seq(n,1L)]

par(mfrow = c(3L,4L))
for (t_ in seq(env$tau)) {
  # Create an igraph object
  g <- graph.data.frame(dtx[t==t_,-1L], directed = TRUE, vertices = seq(env$nI+env$nJ))
  # Set vertex and edge attributes
  V(g)$color <- "lightblue"
  V(g)$size <- 50L  # Set the size of the vertices
  E(g)$arrow.size <- 0.25  # Set the size of the arrows
  # Set edge weights
  E(g)$weight <- dtx$assignment[(t_-1L)*env$nL+seq(env$nL)]
  # Plot the graph
  plot(g, edge.label = E(g)$weight, layout = igraph::layout_on_grid(g) %*% reflex, main = paste0("t = ", t_))
  # Add labels on top left and top right
  text(x = c(-1.33,-1.33, 1.33, 1.33), 
       y = c( 1.25,-0.75, 1.25,-0.75), 
       labels = c("Q", "Q", "D", "D"),
       adj = c(.5, .5))
}
par(mfrow = c(1L,1L))

# Manual inspection
j <- nA
dt[(t==6L) & (action==A_[j])]
