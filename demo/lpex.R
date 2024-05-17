library(TLPR)

env <- new.env()
generate_cssap(env, rate = 10.0, tau = 12L, nCS = 10L, nCO = 1L, nB = 3L, nI = 2L, nJ = 2L)
env$alpha <- c(rep(1.0, env$nI), rep(c(1.0, 2.0), each = env$nJ))

exog <- list(
  "Q" = list("func"=rpois, "params"=list('n'=env$tau*env$nI, 'lambda'=10L)),
  "D" = list("func"=rpois, "params"=list('n'=env$tau*env$nJ, 'lambda'=10L))
)

get_from_routes(env)
get_to_routes(env)

Q <- do.call(exog$Q$func, exog$Q$params)
D <- do.call(exog$Q$func, exog$Q$params)

# Strategic carriers
A1 <- Matrix::spMatrix(ncol = env$nvars, nrow = env$nCS)
for (k in env$CS) {
  A1[k, sum(env$nLc[1L:k])+seq(env$nLc[k+1L])] <- 1L
}
rhs1 <- env$Cb[1L,]
sns1 <- rep("<", env$nCS)

# Spot carriers
A2 <- Matrix::spMatrix(ncol = env$nvars, nrow = env$nCO)
for (k in seq(env$nCO)) {
  A2[k, env$nL_+(k-1L)*env$nL+seq(env$nL)] <- 1L
}
rhs2 <- env$Co[1L,]
sns2 <- rep("<", env$nCO)

rbind(A1, A2)

A12 <- cbind(
  Matrix::spMatrix(ncol = env$nI + 2L*env$nJ, nrow = nrow(A1)+nrow(A2)),
  rbind(A1, A2),
  Matrix::spMatrix(ncol = env$nI + 2L*env$nJ, nrow = nrow(A1)+nrow(A2)))

A3 <- Matrix::spMatrix(ncol = env$nI+2L*env$nJ+env$nvars+env$nI+2L*env$nJ, nrow = env$nI+env$nJ)
for (i in seq(env$nI)) {
  A3[i, i] <- -1L
  idx <- env$from_i[[i]]
  A3[i, env$nI+2L*env$nJ+idx] <- 1L
  A3[i, env$nI+2L*env$nJ+env$nvars+i] <- 1L
}

for (j in seq(env$nJ)) {
  A3[env$nI+j, env$nI+(j-1L)*2L+seq(2L)] <- c(1L, -1L)
  jdx <- env$to_j[[j]]
  A3[env$nI+j, env$nI+2L*env$nJ+jdx] <- 1L
  A3[env$nI+j, env$nI+2L*env$nJ+env$nvars+env$nI+(j-1L)*2L+seq(2L)] <- c(-1L, 1L)
}
rhs3 <- c(Q[seq(env$nI)],D[seq(env$nJ)])
sns3 <- rep("=", env$nI+env$nJ)

A3

A4 <- Matrix::spMatrix(ncol = env$nI, nrow = env$nI)
for (i in seq(env$nI)) {
  A4[i, i] <- 1L
}

A4 <- cbind(
  A4, Matrix::spMatrix(ncol = 2L*env$nJ + env$nvars + env$nI + 2L*env$nJ, nrow = env$nI)
)

rhs4 <- rep(env$R, env$nI)
sns4 <- rep("<", env$nI)

A5 <- Matrix::spMatrix(ncol = 2L*env$nJ + env$nvars, nrow = env$nJ)
for (j in seq(env$nJ)) {
  A5[j, (j-1L)*2L+1L] <- 1L
  jdx <- env$to_j[[j]]
  A5[j, 2L*env$nJ+jdx] <- 1L
}

A5 <- cbind(
  Matrix::spMatrix(ncol = env$nI, nrow = env$nJ),
  A5, Matrix::spMatrix(ncol = env$nI + 2L*env$nJ, nrow = env$nJ)
)

rhs5 <- rep(env$R, env$nJ)
sns5 <- rep("<", env$nJ)

CW <- rep(20.0, env$nI)
CD <- rep(20.0, env$nJ)
CB <- rep(50.0, env$nJ)

obj_ <- c(CW, CD, CB, env$CTb, env$CTo[1L,], CW, CD, CB)

A <- rbind(A12, A3, A4, A5)
rhs <- c(rhs1, rhs2, rhs3, rhs4, rhs5)
sns <- c(sns1, sns2, sns3, sns4, sns5)

model <- list()

model$A          <- A
model$obj        <- obj_
model$modelsense <- 'min'
model$rhs        <- rhs
model$sense      <- sns
model$lb         <- rep(0L,  ncol(model$A))
model$vtype      <- rep('C', ncol(model$A))

model

opt <- gurobi::gurobi(model, params = list(OutputFlag = 0L))
opt$x
opt$pi


offset <- env$nI + 2L*env$nJ
A.tau <- Matrix::spMatrix(nrow = env$tau*nrow(A), ncol = env$tau*(ncol(A)-offset)+offset)
A.tau[seq(nrow(A)),seq(ncol(A))] <- A

obj.tau <- numeric(ncol(A.tau))
obj.tau[seq(ncol(A)-offset)] <- obj_[seq(ncol(A)-offset)]
for (t in seq(2L,env$tau)) {
  rdx <- (t-1L)*nrow(A)+seq(nrow(A))
  cdx <- (t-1L)*(ncol(A)-offset)+seq(ncol(A))
  A.tau[rdx,cdx] <- A
  obj.tau[(t-1L)*(ncol(A)-offset)+seq(ncol(A)-offset)] <- c(CW, CD, CB, env$CTb, env$CTo[t,])
}
obj.tau[seq(ncol(A.tau)-offset+1L,ncol(A.tau))] <- c(CW, CD, CB)

rhs.tau <- rep(rhs, env$tau)
sns.tau <- rep(sns, env$tau)

# rhs <- rep(rhs, 2L)
# sns <- rep(sns, 2L)

model <- list()

model$A          <- A.tau
model$obj        <- obj.tau
model$modelsense <- 'min'
model$rhs        <- rhs.tau
model$sense      <- sns.tau
model$lb         <- rep(0L,  ncol(model$A))
model$vtype      <- rep('I', ncol(model$A))

model

opt <- gurobi::gurobi(model, params = list(OutputFlag = 0L))

S.I  <- matrix(0.0, nrow = (env$tau+1L)*env$nI, ncol = 1L)
S.J  <- matrix(0.0, nrow = (env$tau+1L)*env$nJ, ncol = 1L)
X.I  <- matrix(0.0, nrow = env$tau*env$nI, ncol = 1L)
X.J  <- matrix(0.0, nrow = env$tau*env$nJ, ncol = 1L)
allocation <- matrix(NA, env$tau * env$nvars, ncol = 1L)

blk <- offset + env$nvars
for (t in seq(env$tau+1L)) {
  i <- (t-1L)*env$nI+env$I_
  j <- (t-1L)*env$nJ+env$J_
  idx <- (t-1L)*blk + seq(blk)
  S.I[i,1L] <- opt$x[idx][env$I_]
  S.J[j,1L] <- opt$x[idx][env$nI+env$J_] - opt$x[idx][env$nI+env$nJ+env$J_]

  if (t > env$tau) break
  allocation[(t-1L)*env$nvars+seq(env$nvars),1L] <- opt$x[idx][offset+seq(env$nvars)]
}

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
# Create the data table.
dt <- as.data.table(cbind(
  t       = rep(seq(env$tau), each = env$nvars),
  carrier = carriers[,1L],
  lane    = lanes,
  origin  = env$L[lanes,2L],
  destination = env$L[lanes,1L]+env$nI,
  assignment  = allocation[,1L]))

dtx <- dt[, .(assignment = sum(assignment, na.rm = T)), by = .(t, origin, destination)]

par(mfrow = c(3L,4L))
for (t_ in seq(env$tau)) {
  # Create an igraph object
  g <- graph_from_data_frame(dtx[t==t_,-1L], directed = TRUE, vertices = seq(env$nI+env$nJ))
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
         # Q[(t_-1L)*env$nI+seq(env$nI)],
         # D[(t_-1L)*env$nI+seq(env$nI)],
         rep(rhs3[env$I_], env$tau)[(t_-1L)*env$nI+seq(env$nI)],
         rep(rhs3[env$nI+env$J_], env$tau)[(t_-1L)*env$nI+seq(env$nI)],
         S.I[(t_-1L)*env$nI+seq(env$nI),1L],
         S.J[(t_-1L)*env$nJ+seq(env$nJ),1L]),
       adj = c(.5, .5))
}
par(mfrow = c(1L,1L))
