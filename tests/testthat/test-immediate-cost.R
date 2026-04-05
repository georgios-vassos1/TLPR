library(TLPR)

# ── h.t ───────────────────────────────────────────────────────────────────────

test_that("h.t computes origin holding cost correctly", {
  env <- new.env(parent = emptyenv())
  env$nI <- 2L
  env$nJ <- 1L
  alpha  <- c(2.0, 3.0, 4.0, 5.0)   # nI + 2*nJ = 4

  # S.I = c(1, 2), S.J = 0 (no holding/backorder at destination)
  # holding_origin = 1*2 + 2*3 = 8
  expect_equal(h.t(env, c(1, 2), 0, alpha), 8.0)
})

test_that("h.t computes destination holding cost correctly", {
  env <- new.env(parent = emptyenv())
  env$nI <- 2L
  env$nJ <- 1L
  alpha  <- c(2.0, 3.0, 4.0, 5.0)

  # S.I = c(0, 0), S.J = 3 (positive: holding at destination)
  # holding_dest = 3*4 = 12
  expect_equal(h.t(env, c(0, 0), 3, alpha), 12.0)
})

test_that("h.t computes backorder cost correctly", {
  env <- new.env(parent = emptyenv())
  env$nI <- 2L
  env$nJ <- 1L
  alpha  <- c(2.0, 3.0, 4.0, 5.0)

  # S.I = c(0, 1), S.J = -2 (negative: backorder at destination)
  # holding_origin = 0*2 + 1*3 = 3
  # holding_dest   = max(-2,0)*4 = 0
  # backorder      = -min(-2,0)*5 = 2*5 = 10
  # total = 3 + 0 + 10 = 13
  expect_equal(h.t(env, c(0, 1), -2, alpha), 13.0)
})

test_that("h.t combined holding and backorder", {
  env <- new.env(parent = emptyenv())
  env$nI <- 2L
  env$nJ <- 1L
  alpha  <- c(2.0, 3.0, 4.0, 5.0)

  # S.I = c(1, 2), S.J = 3
  # 1*2 + 2*3 + 3*4 - 0*5 = 2 + 6 + 12 = 20
  expect_equal(h.t(env, c(1, 2), 3, alpha), 20.0)
})

test_that("h.t default alpha uses all-ones", {
  env <- new.env(parent = emptyenv())
  env$nI <- 1L
  env$nJ <- 1L

  # alpha = c(1, 1, 1), S.I = 3, S.J = 2
  # 3*1 + 2*1 - 0*1 = 5
  expect_equal(h.t(env, 3, 2), 5.0)
})

test_that("h.t is zero for zero state", {
  env <- new.env(parent = emptyenv())
  env$nI <- 2L
  env$nJ <- 2L
  alpha  <- runif(6L, 1.0, 10.0)
  expect_equal(h.t(env, c(0, 0), c(0, 0), alpha), 0.0)
})

# ── heuristic_assignment ──────────────────────────────────────────────────────

# Minimal model: 3 vars, 2 constraints: per-var cap <= 10, total = volume
make_cap_model <- function() {
  list(
    A     = rbind(diag(3L), rep(1L, 3L)),
    sense = c("<", "<", "<", "=")
  )
}

test_that("heuristic_assignment returns zero vector when volume is 0", {
  result <- heuristic_assignment(make_cap_model(),
                                 obj_ = c(1, 2, 3),
                                 rhs_ = c(10, 10, 10, 0),
                                 k = 3L, edx = 4L)
  expect_equal(result$status, "HEURISTIC")
  expect_equal(result$x,      c(0, 0, 0))
  expect_equal(result$objval, 0.0)
})

test_that("heuristic_assignment allocation sums to volume", {
  result <- heuristic_assignment(make_cap_model(),
                                 obj_ = c(3, 1, 2),
                                 rhs_ = c(10, 10, 10, 7),
                                 k = 3L, edx = 4L)
  expect_equal(result$status, "HEURISTIC")
  expect_equal(sum(result$x), 7L)
  expect_true(all(result$x >= 0L))
  expect_true(all(result$x <= 10L))
})

test_that("heuristic_assignment respects per-variable capacity", {
  result <- heuristic_assignment(make_cap_model(),
                                 obj_ = c(1, 1, 1),
                                 rhs_ = c(2, 3, 4, 5),
                                 k = 3L, edx = 4L)
  expect_equal(sum(result$x), 5L)
  expect_true(result$x[1L] <= 2L)
  expect_true(result$x[2L] <= 3L)
  expect_true(result$x[3L] <= 4L)
})

# ── capacitated_random_assignment ─────────────────────────────────────────────

test_that("capacitated_random_assignment returns zero vector when volume is 0", {
  result <- capacitated_random_assignment(make_cap_model(),
                                          obj_ = c(1, 2, 3),
                                          rhs_ = c(5, 5, 5, 0),
                                          k = 3L, edx = 4L)
  expect_equal(result$status, "CAPACITATED RANDOM")
  expect_equal(result$x,      c(0, 0, 0))
  expect_equal(result$objval, 0.0)
})

test_that("capacitated_random_assignment allocation sums to volume", {
  set.seed(42L)
  result <- capacitated_random_assignment(make_cap_model(),
                                          obj_ = c(1, 2, 3),
                                          rhs_ = c(5, 5, 5, 8),
                                          k = 3L, edx = 4L)
  expect_equal(result$status, "CAPACITATED RANDOM")
  expect_equal(sum(result$x), 8L)
  expect_true(all(result$x >= 0L))
  expect_true(all(result$x <= 5L))
})
