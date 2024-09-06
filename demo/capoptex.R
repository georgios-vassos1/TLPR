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
env$from_i <- apply(t(env$from_i), 2L, as.integer, simplify = F) # list(c(env$from_i))
env$to_j   <- apply(t(env$to_j),   2L, as.integer, simplify = F) # list(c(env$to_j))

# Reservation cost
env$nSources <- length(env$L_) + env$nL * env$nCO
v      <- rep(0.0, env$nSources * env$tau)
vdx    <- c(outer(seq(env$nSources - 1L), seq(0, env$nSources * (env$tau - 1L), by = env$nSources), "+"))
v[vdx] <- c(8.52, 4.46, 4.25, 9.70) # runif(length(vdx), 3.0, 10.0)

# Refine the contract rates (framework agreement vs option contract)
env$CTb <- env$CTb / 5.0

# Create the model
model <- create_model(env)

# Scenario space specification
env$scndx <- do.call(TLPR::CartesianProductX, c(
  replicate(env$nI,  seq(env$nQ), simplify = FALSE),
  replicate(env$nJ,  seq(env$nD), simplify = FALSE),
  replicate(env$nCO, seq(env$nW), simplify = FALSE)))

env$scnpb <- apply(env$scndx, 1L, function(x) prod(env$Q$prob[x[seq(env$nI)]], env$D$prob[x[env$nI + seq(env$nJ)]], env$W$prob[x[env$nI + env$nJ + 1L]]))

# Choose a scenario
varphidx <- c(9L, 18L, 7L, 10L)
Q <- env$Q$vals[env$scndx[varphidx,1L]]
D <- env$D$vals[env$scndx[varphidx,2L]]
env$CTo[] <- env$W$vals[env$scndx[varphidx,3L]]

# Build the multistage model
ccx <- carrier_capacity_padded(env)
tlx <- transition_logic(env, q = Q[seq(env$nI)], d = D[seq(env$nJ)])
slx <- storage_limits(env, q = Q[seq(env$nI)])

obj_ <- c(env$alpha, env$CTb, env$CTo[1L,], env$alpha)

A   <- rbind(ccx$A, tlx$A, slx$A)
rhs <- c(ccx$rhs, tlx$rhs, slx$rhs)
sns <- c(ccx$sense, tlx$sense, slx$sense)

model <- multiperiod_expansion(env, Q, D, A, obj_, rhs, sns)
model$modelsense <- 'min'
model$vtype <- rep('C', ncol(model$A))

## Capacity optimization
f <- function(model, x, v, capdx, ...) {
  model$rhs[capdx] <- x
  opt <- gurobi::gurobi(model, params = list(OutputFlag = 0L))
  opt$objval + sum(v * x)
}

# 1x1 instance
capdx <- c(outer(seq(env$nCS + env$nCO), seq(0L, 6L * (env$tau - 1L), 6L), '+'))
# 2x2 instance
# capdx <- c(outer(1L:4L, seq(0L, 12L * (env$tau - 1L), 12L), '+'))

x0 <- c(t(cbind2(env$Cb, env$Co)))
v0 <- f(model, x0, v, capdx)

# Define control parameters
control_list <- list(
  maxit = 1000,         # Maximum number of iterations
  factr = 1e3,          # Tolerance for convergence
  pgtol = 1e-9,         # Tolerance for the gradient
  parscale = rep(1.0, length(x0))  # Scaling parameters (if necessary)
)

(x <- optim(
  x0, 
  f, 
  model = model, 
  v     = v,
  capdx = capdx, 
  method = 'L-BFGS-B', 
  lower = 0.0, 
  upper = 10.0, 
  control = control_list
))

CostPerTEU <- function(model, x, v, capdx, ...) {
  model$rhs[capdx] <- x
  opt <- gurobi::gurobi(model, params = list(OutputFlag = 0L))
  (opt$objval + sum(v * x)) / sum(opt$x)
}

# Optimal capacity
CostPerTEU(model, x0, v, capdx)
CostPerTEU(model, round(x$par), v, capdx)

# Capacity performance
M <- 1000000L
capacities <- matrix(sample(seq(0L, 10L), ((env$nCS + env$nCO) * env$tau) * M, replace = TRUE), nrow = M)

cl <- parallel::makeCluster(parallel::detectCores() - 2L)
parallel::clusterExport(cl, c('env', 'model', 'x0', 'x', 'v', 'capdx', 'capacities', 'M'))

timex <- Sys.time()
parallel::parApply(cl, capacities, 1L, CostPerTEU, model = model, v = v, capdx = capdx) -> opcosts
timex <- Sys.time() - timex

parallel::stopCluster(cl)

summary(opcosts)
summary(total_costs)

par(mfrow = c(1L, 2L))

hist(total_costs, main = 'Total Cost Density for 1,000,000 Capacity Vectors', xlab = 'Total Cost', freq = F, ylim = c(0.0, 0.0065))
lines(density(total_costs)$x, density(total_costs)$y, col = 'red', lwd = 2L)
# lines(sort(unique(opcosts)), dnorm(sort(unique(opcosts)), mean = mean(opcosts), sd = sd(opcosts)), col = 'green', lwd = 2L, add = TRUE)

hist(opcosts, main = 'Cost-per-TEU Density for 1,000,000 Capacity Vectors', xlab = 'Cost-per-TEU', freq = F)
lines(density(opcosts)$x, density(opcosts)$y, col = 'red', lwd = 2L)

par(mfrow = c(1L, 1L))

## Approximate capacity optimization
EV <- function(cl, env, model, x, v, scnmat, qdx, ddx, wdx, capdx, ...) {
  # scnmat <- unique(matrix(sample(nrow(env$scndx), env$tau * N, replace = TRUE, prob = env$scnpb), nrow = N))
  prob   <- parallel::parApply(cl, scnmat, 1L, function(idx) prod(env$scnpb[idx]))
  N <- nrow(scnmat)

  # Define the function to be applied in parallel
  job <- function(k) {
    varphidx <- scnmat[k,]
    Q <- env$Q$vals[env$scndx[varphidx, 1L]]
    D <- env$D$vals[env$scndx[varphidx, 2L]]
    env$CTo[] <- env$W$vals[env$scndx[varphidx, 3L]]

    model$rhs[qdx] <- rep(Q, each = 2L)
    model$rhs[ddx] <- D
    model$obj[wdx] <- c(env$CTo)
    model$rhs[capdx] <- x

    gurobi::gurobi(model, params = list(OutputFlag = 0L))$objval
  }

  # Use mclapply to run the function in parallel
  optx <- parallel::parLapply(cl, seq(N), job)

  # Convert the list returned by mclapply to a numeric vector
  optx <- unlist(optx)

  # Calculate the mean of the optx vector
  optcx <- sum(v * x) + optx
  # optcx <- (optcx - min(optcx)) / (max(optcx) - min(optcx))
  sum(optcx * (prob / sum(prob)))
}

qdx <- c(outer(c(3L,5L), seq(0L, 6L * (env$tau - 1L), 6L), '+'))
ddx <- seq(4L, 6L* env$tau, by = 6L)
wdx <- seq(5L, 6L* env$tau, by = 5L)
capdx <- c(outer(1L:2L, seq(0L, 6L * (env$tau - 1L), 6L), '+'))

N <- 1000L
scnmat <- unique(matrix(sample(nrow(env$scndx), env$tau * N, replace = TRUE, prob = env$scnpb), nrow = N))
x0 <- c(t(cbind2(env$Cb, env$Co)))

# cl <- parallel::makeCluster(parallel::detectCores() - 2L)
# parallel::clusterExport(cl, c('env', 'model', 'v', 'x0', 'scnmat', 'qdx', 'ddx', 'wdx', 'capdx'))
# EV(cl, env, model, x0, v, scnmat, qdx, ddx, wdx, capdx)
# parallel::stopCluster(cl)

control_list <- list(
  maxit = 1000L,        # Maximum number of iterations
  factr = 1e3,          # Tolerance for convergence
  pgtol = 1e-9          # Tolerance for the gradient
  # parscale = rep(1.0, length(x0))  # Scaling parameters (if necessary)
)

cl <- parallel::makeCluster(parallel::detectCores() - 2L)
parallel::clusterExport(cl, c('env', 'model', 'x0', 'v', 'scnmat', 'qdx', 'ddx', 'wdx', 'capdx'))
(x.ev <- optim(
  x0,
  EV,
  cl  = cl,
  env = env,
  model = model,
  v = v,
  scnmat = scnmat,
  qdx = qdx,
  ddx = ddx,
  wdx = wdx,
  capdx = capdx,
  method = 'L-BFGS-B',
  lower = 0.0,
  upper = Inf,
  control = control_list
))
parallel::stopCluster(cl)

cl <- parallel::makeCluster(parallel::detectCores() - 2L)
parallel::clusterExport(cl, c('env', 'model', 'x', 'x0', 'v', 'scnmat', 'qdx', 'ddx', 'wdx', 'capdx'))
EV(cl, env, model, x$par, v, scnmat, qdx, ddx, wdx, capdx)
parallel::stopCluster(cl)

##
C.opt <- function(env, model, x, v, varphidx, qdx, ddx, wdx, capdx, ...) {

  Q <- env$Q$vals[env$scndx[varphidx, 1L]]
  D <- env$D$vals[env$scndx[varphidx, 2L]]
  env$CTo[] <- env$W$vals[env$scndx[varphidx, 3L]]

  model$rhs[qdx] <- rep(Q, each = 2L)
  model$rhs[ddx] <- D
  model$obj[wdx] <- c(env$CTo)
  model$rhs[capdx] <- x

  optx <- gurobi::gurobi(model, params = list(OutputFlag = 0L))$objval

  # Calculate the mean of the optx vector
  sum(v * x) + optx
}

compute.regret <- function(env, model, x.ev, v, varphidx, qdx, ddx, wdx, capdx, control_list) {
  x <- optim(
    x.ev$par,
    C.opt,
    env   = env,
    model = model,
    v     = v,
    varphidx = varphidx,
    qdx = qdx,
    ddx = ddx,
    wdx = wdx,
    capdx = capdx,
    method = 'L-BFGS-B',
    lower = 0.0,
    upper = 10.0,
    control = control_list
  )

  C.opt(env, model, x.ev$par, v, varphidx, qdx, ddx, wdx, capdx) - 
    C.opt(env, model, round(x$par), v, varphidx, qdx, ddx, wdx, capdx)
}

cl <- parallel::makeCluster(parallel::detectCores() - 2L)
parallel::clusterExport(cl, c(
  'C.opt', 'compute.regret', 
  'env', 'model', 'x.ev', 'v', 'scnmat', 
  'qdx', 'ddx', 'wdx', 'capdx', 'control_list'))

timex <- Sys.time()

parallel::parApply(
  cl, scnmat, 1L, compute.regret, 
  env = env, model = model, x.ev = x.ev, v = v, 
  qdx = qdx, ddx = ddx, wdx = wdx, capdx = capdx, 
  control_list = control_list) -> regrets

timex <- Sys.time() - timex

parallel::stopCluster(cl)

# Out-of-bag scenarios
N   <- 1000L
oob <- unique(matrix(sample(nrow(env$scndx), env$tau * N, replace = TRUE, prob = env$scnpb), nrow = N))
oob

cl <- parallel::makeCluster(parallel::detectCores() - 2L)
parallel::clusterExport(cl, c(
  'C.opt', 'compute.regret', 
  'env', 'model', 'x.ev', 'v', 'oob', 
  'qdx', 'ddx', 'wdx', 'capdx', 'control_list'))

timex <- Sys.time()

parallel::parApply(
  cl, oob, 1L, compute.regret, 
  env = env, model = model, x.ev = x.ev, v = v, 
  qdx = qdx, ddx = ddx, wdx = wdx, capdx = capdx, 
  control_list = control_list) -> oob.regrets

timex <- Sys.time() - timex

parallel::stopCluster(cl)


par(mfrow = c(1L, 2L))

demp <- density(round(regrets, 2L))
hist(round(regrets, 2L), main = 'In-Sample Regret (1,000 Scenarios)', xlab = 'Regret', freq = F, breaks = 15L)
lines(demp$x[demp$x > 0.0], demp$y[demp$x > 0.0], col = 'red', lwd = 2L)

oob.demp <- density(round(oob.regrets, 2L))
hist(round(oob.regrets, 2L), main = 'Out-of-Sample Regret (1,000 Scenarios)', xlab = 'Regret', freq = F, breaks = 15L)
lines(oob.demp$x[oob.demp$x > 0.0], oob.demp$y[oob.demp$x > 0.0], col = 'red', lwd = 2L)

par(mfrow = c(1L, 1L))

nobs <- min(length(regrets), length(oob.regrets))
inregrets  <- sort(sample(regrets, nobs, replace = F))
outregrets <- sort(sample(oob.regrets, nobs, replace = F))
plot(inregrets, outregrets, main = 'In-Sample vs Out-of-Sample Regret', xlab = 'In-Sample Regret', ylab = 'Out-of-Sample Regret')
lines(inregrets, lm(outregrets ~ inregrets)$fitted.values, col = 'red', lwd = 2L)
lines(outregrets, outregrets, col = 'blue', lwd = 2L, lty = 3L)

# Pre-sample the matrices
N <- 1000L
pre_sampled_matrices <- vector("list", 500L)  # Store the matrices
for (idx in seq_along(pre_sampled_matrices)) {
  pre_sampled_matrices[[idx]] <- unique(matrix(sample(nrow(env$scndx), env$tau * N, replace = TRUE, prob = env$scnpb), nrow = N))
}

cl <- parallel::makeCluster(parallel::detectCores() - 2L)
parallel::clusterExport(cl, c('env', 'model', 'x0', 'x', 'v', 'pre_sampled_matrices', 'qdx', 'ddx', 'wdx', 'capdx'))

# Now use the pre-sampled matrices in the main loop
results <- numeric(500L)
for (idx in seq_along(results)) {
  results[idx] <- EV(cl, env, model, x = x0, v = v, scnmat = pre_sampled_matrices[[idx]], qdx, ddx, wdx, capdx)
}

parallel::stopCluster(cl)

hist(results, main = 'Histogram of the Expected Value', xlab = 'Expected Value')
