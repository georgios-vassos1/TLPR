library(TLPR)

devtools::build("~/drayage/TLPR")
devtools::document("~/drayage/TLPR")
devtools::install("~/drayage/TLPR", reload = F)

nSdx <- 100L
nAdx <-  5L
nUdx <-  80L

n    <- nSdx * nAdx * nUdx
p    <- 5L
mu   <- runif(p, 1.0, 10.0)
covm <- generatePositiveDefiniteMatrix(p)

TLPR::rmvnorm(n, p, mu, covm, 8L)

points(X[1L:(n-1L)], X[2L:n])

microbenchmark::microbenchmark(
  TLPR::rmvnorm(n, p, mu, covm, 8L),
  MASS::mvrnorm(n, mu, covm),
  times = 100L
)

## Cartesian product benchmarks
microbenchmark::microbenchmark(
  data.table::CJ(1:50, 1:50, 1:50, 1:50),
  TLPR::CartesianProductX(1:50, 1:50, 1:50, 1:50),
  TLPR::CartesianProductLB(1:50, 1:50, 1:50, 1:50),
  # expand.grid(1:50, 1:50, 1:50, 1:50),
  times = 100L
)

## Dynamic program
env <- new.env()

json_path <- "/Users/gva/drayage/TLPR/src/instance3x2_001.json"
jsonlite::fromJSON(json_path) |>
  list2env(envir = env)

TLPR::computeEnvironmentSTL(json_path, seq(0L, env$nSI), numThreads = 10L)

args <- replicate(env$nI, env$Q$vals, simplify = FALSE)
do.call(TLPR::CartesianProduct, args) |> nrow() -> nInflow
env$nA <- env$max.S

looper <- function(env) {
  Sdx <- matrix(NA, nrow = env$nSdx * env$nA * nInflow, ncol = env$nI + env$nJ)
  for (i in seq(env$nSdx)) {
    for (j in seq(env$nA)) {
      for (k in seq(nInflow)) {
        Sdx[i,] <- c(rep(3L, env$nI), rep(12L, env$nJ))
      }
    }
  }
  Sdx
}

## Measure time of looper
microbenchmark::microbenchmark(
  looper(env),
  times = 50L
)

result <- TLPR::optimizeModelFromJSON(json_path)

model <- create_model(env)
optx <- TLPR::optimal_assignment(
  model, c(env$CTb, env$CTo[1L,]), 
  c(env$Cb[1L,], env$Co[1L,], c(24L, 22L), c(27L, 37L, 29L), 10L)
)

optx$x
optx$objval
result$objval
result$x

microbenchmark::microbenchmark(
  result <- TLPR::optimizeModelFromJSON(json_path),
  TLPR::optimal_assignment(model, c(env$CTb, env$CTo[1L,]), c(env$Cb[1L,], env$Co[1L,], result$limits[env$nI + env$J_], result$limits[env$I_], 40L)),
  times = 100L
)
