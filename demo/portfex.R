library(TLPR)

## Demo
# Initiate instance
env <- new.env()
generate_cssap(env, rate = 10.0, tau = 4L, nCS = 10L, nCO = 1L, nB = 3L, nI = 2L, nJ = 2L)

get_from_routes(env)
get_to_routes(env)

# Generate optimization model
Idx   <- c(1L,2L,3L,4L)
model <- create_model(env, Idx)
# Action space 
A_ <- seq(0L, max(sweep(env$Cb, 1L, env$Co, '+')), by = 5L)
nA <- length(A_)

S0 <- t(rpois(env$nI*env$nJ, 10L))
Q  <- rpois(env$tau*env$nI, 10L)
D  <- rpois(env$tau*env$nJ, 10L)

cpx <- eval_portfolio(env, A_, S0, Q, D, plot = TRUE)

# matplot(A_, t(cpx$cost[[2L]]), type = 'l')

# idx <- rep(seq(nA), each = nA)
# Y <- apply(cpx$cost[[3L]], 2L, function(x) tapply(x, idx, mean, na.rm = T))
# matplot(A_, t(Y), type = 'l')
# matplot(A_, t(cpx$cost[[3L]]), type = 'l')

# idx <- rep(rep(seq(nA), each = nA), nA)
# Y <- apply(cpx$cost[[4L]], 2L, function(x) tapply(x, idx, mean, na.rm = T))
# matplot(A_, t(Y), type = 'l')
# matplot(A_, t(cpx$cost[[4L]]), type = 'l')

t <- 2L
j <- 5L

Sx <- c(0L,0L,0L,0L)

xspan <- seq(0L, 30L, by = 5L)
yspan <- seq(0L, 30L, by = 5L)
Qx <- as.matrix(unname(expand.grid(xspan, yspan)))
cost <- rep(NaN, nrow(Qx))
for (k in seq(nrow(Qx))) {
  cost[k] <- eval_portfolio_t(env, A_[j], t, Sx, Qx[k,])$cost
}

cbind2(Qx, cost) |>
  as.data.frame() |>
  reshape2::dcast(V2 ~ V1, value.var = "V3") |>
  dplyr::select(-1L) |>
  as.matrix() |>
  apply(1L, rev) |>
  heatmap(Rowv=NA, Colv=NA, labRow = rev(xspan), labCol = yspan, col = cm.colors(256), scale="none")
# z <- matrix(cost, nrow = length(xspan), ncol = length(yspan))
# plot3D::persp3D(x = xspan, y = yspan, z = z, col = "lightblue", 
#                 theta = 30, phi = 30, xlab = "X", ylab = "Y", zlab = "Z", main = "3D Surface Plot")

library(data.table)
library(igraph)

## Draw the feasible region and show we must only calculate cost on boundary areas
# data.table::CJ(seq(0L,30L),seq(0L,30L),seq(0L,30L),seq(0L,30L))

stateless_portf_eval_data <- function(env, A_, assignments) {
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

stateless_portf_eval_data(env, A_, cpx$assignment)

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
dt[(t==2L) & (action==A_[j])]
