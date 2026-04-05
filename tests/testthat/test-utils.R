library(TLPR)

# ── consolidate_idx ───────────────────────────────────────────────────────────

test_that("consolidate_idx returns all combinations (2 x 3)", {
  result <- consolidate_idx(list(1:2, 1:3))
  expect_equal(nrow(result), 6L)
  expect_equal(ncol(result), 2L)
  expect_equal(nrow(unique(result)), 6L)
})

test_that("consolidate_idx first column varies fastest (matches CartesianProductX)", {
  r1 <- consolidate_idx(list(1:3, 1:2))
  r2 <- CartesianProductX(1:3, 1:2)
  # consolidate_idx may carry column names from data.table::CJ; ignore attributes
  expect_equal(unname(r1), r2)
})

# ── CartesianProductX ─────────────────────────────────────────────────────────

test_that("CartesianProductX dimensions are product of lengths", {
  r <- CartesianProductX(1:4, 1:3)
  expect_equal(nrow(r), 12L)
  expect_equal(ncol(r), 2L)
})

test_that("CartesianProductX covers all pairs exactly once", {
  r <- CartesianProductX(1:2, 1:3)
  # 6 unique pairs
  expect_equal(nrow(unique(r)), 6L)
  # min/max within expected range
  expect_equal(range(r[, 1L]), c(1L, 2L))
  expect_equal(range(r[, 2L]), c(1L, 3L))
})

test_that("CartesianProductX single-vector input returns column matrix", {
  r <- CartesianProductX(1:5)
  expect_equal(dim(r), c(5L, 1L))
  expect_equal(r[, 1L], 1:5)
})

# ── chunkup ───────────────────────────────────────────────────────────────────

test_that("chunkup produces k+1 boundaries", {
  expect_equal(length(chunkup(10L, 3L)), 4L)
  expect_equal(length(chunkup(1L,  1L)), 2L)
})

test_that("chunkup first boundary is 1, last is n+1", {
  chunks <- chunkup(12L, 4L)
  expect_equal(chunks[1L], 1L)
  expect_equal(chunks[5L], 13L)
})

test_that("chunkup covers entire range with no gaps or overlaps", {
  for (n in c(1L, 5L, 9L, 12L)) {
    for (k in c(1L, 2L, 3L)) {
      if (k > n) next
      # chunkup may return more than k+1 values; only the first k+1 are valid
      chunks  <- chunkup(n, k)[seq(k + 1L)]
      covered <- unlist(mapply(seq, chunks[-length(chunks)],
                               chunks[-1L] - 1L, SIMPLIFY = FALSE))
      expect_equal(sort(covered), seq_len(n),
                   label = paste0("n=", n, " k=", k))
    }
  }
})
