library(TLPR)

env <- new.env()
generate_cssap(env, rate = 10.0, tau = 12L, nCS = 10L, nCO = 1L, nB = 3L, nI = 2L, nJ = 2L)

env$alpha <- c(rep(1.0, env$nI), rep(c(1.0, 2.0), each = env$nJ))

# Choose the effective constraints
#   Constraints 2L and 3L depend on the state of the system.
#   This entails significant computational burden in evaluating the portfolio contract.
Idx    <- c(1L,2L,3L,4L)
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
  "D" = list("func"=rpois, "params"=list('n'=env$tau*env$nJ, 'lambda'=10L))
)

get_from_routes(env)
get_to_routes(env)

## Single simulation run
result <- simulate_system(env, policy, args, exog, rhs_dx, correction = FALSE)

## Plot solution
library(data.table)
library(igraph)

n <- 2L
reflex <- diag(n)[,seq(n,1L)]

lanes    <- c(env$L_, seq(env$nL))
carriers <- rep(seq_along(env$winner), env$nLc[-1L])
rbind2(
  cbind2(unlist(carriers), unlist(env$B)),
  cbind2(rep(env$nCS+seq(env$nCO), each = env$nL), seq(env$nL))
) -> carriers
# carriers[,1L] <- rep(seq_along(rle(carriers[,1L])$values), rle(carriers[,1L])$lengths)
# Create the data table.
dt <- as.data.table(cbind(
  t       = rep(seq(env$tau), each = env$nvars),
  carrier = carriers[,1L],
  lane    = lanes,
  origin  = env$L[lanes,2L],
  destination = env$L[lanes,1L]+env$nI,
  assignment  = result$allocation[,1L]))
# Add capacity and cost information.
dt[, capacity := cbind2(env$Cb,env$Co)[t,carrier], by = .(t, carrier)]
dt[, cost := c(env$CTb,env$CTo[t,]), by = .(t)]

dtx <- dt[, .(assignment = sum(assignment, na.rm = T)), by = .(t, origin, destination)]

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
  text(x = c(-1.34,-1.34, 1.34, 1.34, -1.00,-1.00, 1.00, 1.00), 
       y = c(-1.00, 1.00,-1.00, 1.00, -0.60, 0.60,-0.60, 0.60), 
       labels = c(
         result$Q[(t_-1L)*env$nI+seq(env$nI)],
         result$D[(t_-1L)*env$nI+seq(env$nI)],
         result$S.I[(t_-1L)*env$nI+seq(env$nI),1L],
         result$S.J[(t_-1L)*env$nJ+seq(env$nJ),1L]),
       adj = c(.5, .5))
}
par(mfrow = c(1L,1L))

## Statistics
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
matplot(stats_, type = 'l', 
        lwd = rep(c(2L,1L,1L), 2L), lty = 1L, pch = 19L, ylab = "Cumulative Cost", col = rep(c(1L,2L,3L), each = 3L))

## Single run plots
matplot(apply(result$cost, 2L, cumsum), type = 'b', lwd = 2L, lty = 1L, pch = 19L, ylab = "Cumulative Cost")

pdx <- 1L
idx <- outer((seq(env$tau)-1L)*env$nI, seq(env$nI), '+')
i <- 2L
plot(result$S.I[idx[,i],pdx], type = 's', ylim = c(min(result$S.I[,pdx],result$X.I[,pdx]), max(result$S.I[,pdx],result$X.I[,pdx])))
lines(result$X.I[idx[,i],pdx], type = 's', col = 2L)
lines(result$Q[idx[,i]], type = 's', col = 3L)
