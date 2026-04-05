library(TLPR)

# ── Shared test environment ───────────────────────────────────────────────────
# Minimal 1x1 system: nI=1, nJ=1, nCO=1, R=2, tau=2
# stateSupport        = {0,1,2}       → 3 values
# extendedStateSupport = {-2,-1,0,1,2} → 5 values
# nSdx = 3 * 5 = 15
# actionSupport = {0,1,2} → nAdx = 3
# Q,D: 2-point uniform; W: degenerate → nScen = 2*2*1 = 4

make_dp_env <- function() {
  env        <- new.env(parent = emptyenv())
  env$nI     <- 1L
  env$nJ     <- 1L
  env$nCO    <- 1L
  env$R      <- 2L
  env$tau    <- 2L
  env$I_     <- 1L
  env$J_     <- 1L
  env$alpha  <- c(1.0, 1.0, 1.0)   # length nI + 2*nJ = 3
  env$Q      <- list(vals = c(0L, 1L), prob = c(0.5, 0.5))
  env$D      <- list(vals = c(0L, 1L), prob = c(0.5, 0.5))
  env$W      <- list(vals = c(10.0),   prob = c(1.0))
  env
}

# ── dp_config ─────────────────────────────────────────────────────────────────

test_that("dp_config sets state support arrays correctly", {
  env <- make_dp_env()
  dp_config(env)

  expect_equal(env$stateSupport,         0L:2L)
  expect_equal(env$extendedStateSupport, -2L:2L)
})

test_that("dp_config Sdx has correct dimensions", {
  env <- make_dp_env()
  dp_config(env)

  # (R+1) origin states × (2R+1) destination states
  expect_equal(env$nSdx, 15L)
  expect_equal(nrow(env$Sdx), 15L)
  expect_equal(ncol(env$Sdx), 2L)  # nI + nJ
})

test_that("dp_config action space is correct", {
  env <- make_dp_env()
  dp_config(env)

  expect_equal(env$actionSupport, 0L:2L)
  expect_equal(env$nAdx,          3L)
  expect_equal(env$Adx,           1L:3L)
})

test_that("dp_config exogenous cardinalities are correct", {
  env <- make_dp_env()
  dp_config(env)

  expect_equal(env$nQ, 2L)
  expect_equal(env$nD, 2L)
  expect_equal(env$nW, 1L)
})

test_that("dp_config marginal index dimensions are correct", {
  env <- make_dp_env()
  dp_config(env)

  expect_equal(env$nQdx, 2L)
  expect_equal(env$nDdx, 2L)
  expect_equal(env$nWdx, 1L)
  expect_equal(env$nScen, 4L)  # nQdx * nDdx * nWdx
})

test_that("dp_config joint scenario probabilities sum to 1", {
  env <- make_dp_env()
  dp_config(env)

  expect_equal(sum(env$scnpb), 1.0, tolerance = 1e-12)
  expect_true(all(env$scnpb > 0))
  expect_equal(length(env$scnpb), 4L)
})

test_that("dp_config joint scenario table has correct shape", {
  env <- make_dp_env()
  dp_config(env)

  # scndx: nQ*nD*nW rows, (nI + nJ + nCO) cols
  expect_equal(nrow(env$scndx), 4L)
  expect_equal(ncol(env$scndx), 3L)
})

test_that("dp_config uniform Q and D give equal scenario probabilities", {
  env <- make_dp_env()
  dp_config(env)

  # All 4 scenarios should have prob 0.25
  expect_equal(env$scnpb, rep(0.25, 4L), tolerance = 1e-12)
})

test_that("dp_config flowKeys encode all nScen scenarios bijectively", {
  env <- make_dp_env()
  dp_config(env)

  # Compute mixed-radix keys: nQ^(0:(nI-1)), then nQ^nI * nD^(0:(nJ-1)), then ...
  env$flowKeys <- c(
    env$nQ ^ (0L:(env$nI - 1L)),
    (env$nQ ^ env$nI) * (env$nD ^ (0L:(env$nJ - 1L))),
    (env$nQ ^ env$nI) * (env$nD ^ env$nJ) * env$nW ^ (0L:(env$nCO - 1L))
  )

  kdx_all <- integer(env$nScen)
  p <- 0L
  for (k1 in seq(env$nQdx)) {
    for (k2 in seq(env$nDdx)) {
      for (k3 in seq(env$nWdx)) {
        p <- p + 1L
        kdx_all[p] <- sum(c(env$Qdx[k1, ] - 1L,
                             env$Ddx[k2, ] - 1L,
                             env$Wdx[k3, ] - 1L) * env$flowKeys) + 1L
      }
    }
  }

  # Every index in 1:nScen should appear exactly once
  expect_equal(sort(kdx_all), 1L:env$nScen)
})

# ── dynamic_programming ───────────────────────────────────────────────────────

# Wrapper: suppress expected max/min-on-all-NA warnings from dynamic_programming
# when the transit table has no feasible entries.
dp_no_transitions <- function(env) {
  suppressWarnings(dynamic_programming(env, matrix(NA_real_,
    nrow = env$tau * env$nSdx * env$nAdx * env$nScen, ncol = 6L)))
}

make_all_na_transit <- function(env) {
  matrix(NA_real_,
         nrow = env$tau * env$nSdx * env$nAdx * env$nScen,
         ncol = 6L)
}

test_that("dynamic_programming returns correctly named list", {
  env <- make_dp_env()
  dp_config(env)
  result <- dp_no_transitions(env)

  expect_named(result, c("V", "Q", "pi_star", "pi_rand"))
})

test_that("dynamic_programming output matrices have correct dimensions", {
  env <- make_dp_env()
  dp_config(env)
  result <- dp_no_transitions(env)

  expect_equal(dim(result$V),        c(env$tau + 1L, env$nSdx))
  expect_equal(dim(result$Q),        c(env$tau,       env$nSdx * env$nAdx))
  expect_equal(dim(result$pi_star),  c(env$tau,       env$nSdx * env$nAdx))
  expect_equal(dim(result$pi_rand),  c(env$tau,       env$nSdx * env$nAdx))
})

test_that("dynamic_programming terminal values are all finite", {
  env <- make_dp_env()
  dp_config(env)
  result <- dp_no_transitions(env)

  expect_true(all(is.finite(result$V[env$tau + 1L, ])))
})

test_that("dynamic_programming terminal value is 0 at zero state", {
  env <- make_dp_env()
  dp_config(env)
  result <- dp_no_transitions(env)

  # Find the state index where S.I = 0 and S.J = 0
  zero_idx <- which(
    env$stateSupport[env$Sdx[, env$I_]] == 0L &
    env$extendedStateSupport[env$Sdx[, env$nI + env$J_]] == 0L
  )
  expect_equal(result$V[env$tau + 1L, zero_idx], 0.0, tolerance = 1e-12)
})

test_that("dynamic_programming terminal value is negative for positive inventory", {
  # alpha > 0 ⟹ any non-zero inventory incurs cost ⟹ V(terminal) < 0
  env <- make_dp_env()
  env$alpha <- c(5.0, 4.0, 8.0)
  dp_config(env)
  result <- dp_no_transitions(env)

  # S.I = 2, S.J = 0: cost = 2*5 = 10 → V = -10
  max_si_zero_sj <- which(
    env$stateSupport[env$Sdx[, env$I_]] == 2L &
    env$extendedStateSupport[env$Sdx[, env$nI + env$J_]] == 0L
  )
  expect_equal(result$V[env$tau + 1L, max_si_zero_sj], -10.0, tolerance = 1e-12)
})

test_that("dynamic_programming terminal value for backorder state is negative", {
  env <- make_dp_env()
  env$alpha <- c(5.0, 4.0, 8.0)
  dp_config(env)
  result <- dp_no_transitions(env)

  # S.I = 0, S.J = -2: backorder cost = 2*8 = 16 → V = -16
  backorder_idx <- which(
    env$stateSupport[env$Sdx[, env$I_]] == 0L &
    env$extendedStateSupport[env$Sdx[, env$nI + env$J_]] == -2L
  )
  expect_equal(result$V[env$tau + 1L, backorder_idx], -16.0, tolerance = 1e-12)
})

test_that("dynamic_programming with no feasible transitions leaves Q all NA", {
  env <- make_dp_env()
  dp_config(env)
  result <- dp_no_transitions(env)

  # No feasible transitions → Q is all NA
  expect_true(all(is.na(result$Q)))
  # V[1:tau] is -Inf: max(all-NA, na.rm=TRUE) = -Inf
  expect_true(all(!is.finite(result$V[seq(env$tau), ])))
})

test_that("dynamic_programming pi matrices sum to 0 or 1 per state block", {
  env <- make_dp_env()
  dp_config(env)
  # Use all-NA transit: fcoalesce converts NA → 0
  result <- dp_no_transitions(env)

  for (t in seq(env$tau)) {
    for (i in seq(env$nSdx)) {
      idx   <- (i - 1L) * env$nAdx + seq(env$nAdx)
      s_sum <- sum(result$pi_star[t, idx])
      r_sum <- sum(result$pi_rand[t, idx])
      # With all-NA Q: fcoalesce gives 0, so sums are 0
      expect_true(s_sum == 0.0 || abs(s_sum - 1.0) < 1e-12,
                  label = paste0("pi_star row sum t=", t, " i=", i))
      expect_true(r_sum == 0.0 || abs(r_sum - 1.0) < 1e-12,
                  label = paste0("pi_rand row sum t=", t, " i=", i))
    }
  }
})

# ── get_adjustment_weights ────────────────────────────────────────────────────

test_that("get_adjustment_weights with nJ=1 returns c(1, nSI^1, nSI^2)", {
  # parent = baseenv() so that ^ and : are accessible via with()
  env <- new.env(parent = baseenv())
  env$nSI <- 5L
  env$nI  <- 2L
  env$nJ  <- 1L
  env$nSJ <- 3L

  expect_equal(get_adjustment_weights(env), c(1L, 5L, 25L))
})

test_that("get_adjustment_weights with nJ>1 appends sj_powers", {
  env <- new.env(parent = baseenv())
  env$nSI <- 2L
  env$nI  <- 2L
  env$nJ  <- 2L
  env$nSJ <- 3L

  # si_powers = c(2, 4); base = 4; sj_powers = 4*3 = 12
  expect_equal(get_adjustment_weights(env), c(1L, 2L, 4L, 12L))
})

# ── rolling_dp_cx ─────────────────────────────────────────────────────────────

test_that("rolling_dp_cx returns same structure as dynamic_programming", {
  env <- make_dp_env(); dp_config(env)
  local_mocked_bindings(
    bellmanUpdateCx = function(jsonFile, t, stateSupport, flowSupport,
                                scnpb, alpha, V_next, numThreads) {
      n  <- length(V_next)
      na <- env$nAdx
      list(V_t       = rep(0, n),
           Q_t       = rep(0, n * na),
           pi_star_t = rep(0, n * na),
           pi_rand_t = rep(0, n * na))
    }, .package = "TLPR")
  r <- rolling_dp_cx(env, "dummy.json", 1L)
  expect_named(r, c("V", "Q", "pi_star", "pi_rand"))
  expect_equal(dim(r$V),       c(env$tau + 1L, env$nSdx))
  expect_equal(dim(r$Q),       c(env$tau,      env$nSdx * env$nAdx))
  expect_equal(dim(r$pi_star), c(env$tau,      env$nSdx * env$nAdx))
  expect_equal(dim(r$pi_rand), c(env$tau,      env$nSdx * env$nAdx))
})

test_that("rolling_dp_cx terminal V[tau+1,] matches dynamic_programming", {
  env <- make_dp_env(); dp_config(env)
  local_mocked_bindings(
    bellmanUpdateCx = function(jsonFile, t, stateSupport, flowSupport,
                                scnpb, alpha, V_next, numThreads) {
      n  <- length(V_next)
      na <- env$nAdx
      list(V_t       = rep(0, n),
           Q_t       = rep(NA_real_, n * na),
           pi_star_t = rep(0, n * na),
           pi_rand_t = rep(0, n * na))
    }, .package = "TLPR")
  rolling <- rolling_dp_cx(env, "dummy.json", 1L)
  ref     <- suppressWarnings(dynamic_programming(env, make_all_na_transit(env)))
  expect_equal(rolling$V[env$tau + 1L, ], ref$V[env$tau + 1L, ], tolerance = 1e-12)
})

# ── get_adjustment_weights ────────────────────────────────────────────────────

test_that("get_adjustment_weights length is nI + nJ", {
  env <- new.env(parent = baseenv())
  env$nSI <- 3L
  env$nI  <- 3L
  env$nJ  <- 2L
  env$nSJ <- 4L

  w <- get_adjustment_weights(env)
  expect_equal(length(w), env$nI + env$nJ)
})
