library(TLPR)

# ── input_configuration ───────────────────────────────────────────────────────

test_that("input_configuration sets tau, nI, nJ, nCS defaults", {
  env <- new.env(parent = emptyenv())
  input_configuration(env, tau = 4L, nI = 2L, nJ = 3L, nCS = 5L)

  expect_equal(env$tau, 4L)
  expect_equal(env$nI,  2L)
  expect_equal(env$nJ,  3L)
  expect_equal(env$nCS, 5L)
})

test_that("input_configuration sets nL = nI * nJ", {
  env <- new.env(parent = emptyenv())
  input_configuration(env, tau = 1L, nI = 2L, nJ = 3L, nCS = 1L)

  expect_equal(env$nL, 6L)
})

test_that("input_configuration sets I_, J_, CS correctly", {
  env <- new.env(parent = emptyenv())
  input_configuration(env, tau = 1L, nI = 3L, nJ = 2L, nCS = 4L)

  expect_equal(env$I_, 1L:3L)
  expect_equal(env$J_, 1L:2L)
  expect_equal(env$CS, 1L:4L)
})

test_that("input_configuration accepts a config list", {
  env <- new.env(parent = emptyenv())
  cfg <- list(tau = 6L, nI = 2L, nJ = 3L)
  input_configuration(env, config = cfg)

  expect_equal(env$tau, 6L)
  expect_equal(env$nI,  2L)
  expect_equal(env$nJ,  3L)
})

test_that("input_configuration defaults tau=12 and nI=nJ=nCS=2 when no args", {
  env <- new.env(parent = emptyenv())
  # Suppress random nB sampling side-effect
  set.seed(1L)
  input_configuration(env)

  expect_equal(env$tau, 12L)
  expect_equal(env$nI,   2L)
  expect_equal(env$nJ,   2L)
  expect_equal(env$nCS,  2L)
})

# ── simulate_auction + init_env ───────────────────────────────────────────────

test_that("simulate_auction populates L, B, winner", {
  env <- new.env(parent = emptyenv())
  set.seed(42L)
  input_configuration(env, tau = 4L, nI = 2L, nJ = 2L, nCS = 2L, nB = 3L)
  simulate_auction(env)

  expect_true(is.matrix(env$L))
  expect_equal(ncol(env$L), 2L)           # (origin, destination)
  expect_equal(nrow(env$L), env$nL)       # nI * nJ lanes
  expect_equal(length(env$B), env$nB)
  expect_true(is.list(env$winner))
})

test_that("init_env sets R, nvars, from_i, to_j after auction", {
  env <- new.env(parent = emptyenv())
  set.seed(7L)
  generate_cssap(env, tau = 2L, nI = 1L, nJ = 1L, nCS = 2L, nB = 2L, nCO = 1L)

  expect_true(is.numeric(env$R) && env$R > 0L)
  expect_true(is.integer(env$nvars) || is.numeric(env$nvars))
  expect_equal(length(env$from_i), env$nI)
  expect_equal(length(env$to_j),   env$nJ)
})

test_that("from_i and to_j indices are within 1:nvars", {
  env <- new.env(parent = emptyenv())
  set.seed(99L)
  generate_cssap(env, tau = 2L, nI = 2L, nJ = 2L, nCS = 2L, nB = 4L, nCO = 1L)

  all_from <- unlist(env$from_i)
  all_to   <- unlist(env$to_j)

  expect_true(all(all_from >= 1L & all_from <= env$nvars))
  expect_true(all(all_to   >= 1L & all_to   <= env$nvars))
})
