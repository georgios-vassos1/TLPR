library(TLPR)

env <- new.env()
# generate_cssap(env, rate = 10.0, tau = 12L, nB = 8L, nCS = 10L, nCO = 1L, nI = 1L, nJ = 1L)
# env$alpha <- c(rep(0.0, env$nI), rep(c(0.0, 0.0), each = env$nJ))

json_path <- "~/drayage/TLPR/src/instances/instance1x1_12_001.json"
jsonlite::fromJSON(json_path) |>
  list2env(envir = env)
# Necessary transformation when importing from json
env$from_i <- list(c(env$from_i))
env$to_j   <- list(c(env$to_j))

# Create the gurobi model
model <- create_model(env)

# Define the list of allocation policies
policy <- list(
  "myopic" = optimal_assignment,
  "random" = random_assignment,
  "capacitated_random" = capacitated_random_assignment,
  "heuristic" = heuristic_assignment
)
npi <- length(policy)

# Define the exogenous state variables
exog <- list(
  "Q" = list("func"=sample, "params"=list('x' = env$Q$vals, "size" = env$tau*env$nI, "prob" = env$Q$prob, "replace" = TRUE)),
  "D" = list("func"=sample, "params"=list('x' = env$D$vals, "size" = env$tau*env$nI, "prob" = env$D$prob, "replace" = TRUE))
)
# Sample the spot market vector
env$CTo <- matrix(TLPR:::sample_exogenous_state(sample, list('x' = env$W$vals, "size" = env$tau*env$nCO, "prob" = env$W$prob, "replace" = TRUE)), ncol = env$nCO)

# Initialize the arguments for the simulation
args <- list(
  model = model,
  obj_ = NULL,
  rhs_ = NULL,
  n = NULL,
  k = env$nvars,
  edx = nrow(model$A)
)

## Single simulation run
(result <- simulate_system(env, random_volume, policy, args)) # , env$Q$vals[env$scndx[varphidx,1L]], env$D$vals[env$scndx[varphidx,2L]]))

## Plot transition solution
# Configuration for different instance grids
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

library(data.table)
library(igraph)

# Select policy index to visualize
pdx <- 1L

# Transpose grid if necessary by uncommenting the right side
reflex <- diag(2L) # [,seq(n,1L)]

# Create the dataframe for the graph plot
lanes    <- c(env$L_, seq(env$nL))
carriers <- rep(seq_along(env$winner), env$nLc[-1L])
rbind2(
  cbind2(unlist(carriers), unlist(env$B)),
  cbind2(rep(env$nCS+seq(env$nCO), each = env$nL), seq(env$nL))
) -> carriers
# carriers[,1L] <- rep(seq_along(rle(carriers[,1L])$values), rle(carriers[,1L])$lengths)

dt <- as.data.table(cbind(
  t       = rep(seq(env$tau), each = env$nvars),
  carrier = carriers[,1L],
  lane    = lanes,
  origin  = env$L[lanes,1L],
  destination = env$L[lanes,2L]+env$nI,
  assignment  = result$allocation[,pdx]))
# Add capacity and cost information.
dt[, capacity := cbind2(env$Cb,env$Co)[t,carrier], by = .(t, carrier)]
dt[, cost := c(env$CTb,env$CTo[t,]), by = .(t)]

dtx <- dt[, .(assignment = sum(assignment, na.rm = T)), by = .(t, origin, destination)]

# Select configuration
setoo <- paste0(env$nI, "x", env$nJ)

par(mfrow = c(3L,4L), mai = c(0.2, 0.2, 0.3, 0.2))
for (t_ in seq(env$tau)) {
  # Create an igraph object
  g <- graph_from_data_frame(dtx[t==t_,-1L], directed = TRUE, vertices = seq(env$nI+env$nJ))
  # Set vertex and edge attributes
  V(g)$color <- "lightblue"
  V(g)$size <- 80L  # Set the size of the vertices
  V(g)$label.cex <- 1.25  # Increase the font size inside the nodes
  E(g)$arrow.size <- 0.5  # Set the size of the arrows
  # Set edge weights
  E(g)$weight <- dtx$assignment[(t_-1L)*env$nL+seq(env$nL)]
  # Plot the graph
  plot(
    g,
    edge.label = E(g)$weight,
    layout = igraph::layout_on_grid(g),
    xlim = c(-1.5, 1.5),
    edge.label.cex = 1.25)

  # Add the title closer to the graph
  title(main = paste0("t = ", t_), line = -4.0, cex.main = 1.5)  # Reduce space between title and graph

  # Add labels on top left and top right
  text(x = coords[[setoo]]$x,
       y = coords[[setoo]]$y,
       labels = c(
         result$Q[(t_-1L)*env$nI+env$I_],
         result$D[(t_-1L)*env$nJ+env$J_],
         result$S.I[(t_-1L)*env$nI+env$I_,pdx],
         result$S.J[(t_-1L)*env$nJ+env$J_,pdx]),
       adj = c(.5, .5),
       cex = 1.25)
}
par(mfrow = c(1L,1L))

## Statistics
N <- 50L
costs <- matrix(NA, nrow = env$tau, ncol = npi*N)
i <- 1L
while (i <= N) {
  print(i)
  result <- simulate_system(env, random_volume, policy, args)
  if (result$status) {
    costs[,(i-1L)*npi+seq(npi)] <- result$cost
    i <- i + 1L
  }
}

stats_ <- do.call(cbind, lapply(seq(npi), function(x) compute_stats(costs[,seq(x, npi*N, by = npi)], N = 1L)))
matplot(stats_, type = 'l', 
        lwd = rep(c(2L,1L,1L), npi), lty = 1L, pch = 19L, ylab = "Cumulative Cost", col = rep(seq(npi), each = 3L))

## Single run plots
matplot(apply(result$cost, 2L, cumsum), type = 'b', lwd = 2L, lty = 1L, pch = 19L, ylab = "Cumulative Cost")

pdx <- 1L
idx <- outer((seq(env$tau)-1L)*env$nI, env$I_, '+')
i <- 2L
plot(result$S.I[idx[,i],pdx], type = 's', ylim = c(min(result$S.I[,pdx],result$X.I[,pdx]), max(result$S.I[,pdx],result$X.I[,pdx])))
lines(result$X.I[idx[,i],pdx], type = 's', col = 2L)
lines(result$Q[idx[,i]], type = 's', col = 3L)
