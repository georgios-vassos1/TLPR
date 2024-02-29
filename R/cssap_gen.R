#' Get Lanes
#' 
#' Generates a matrix of lanes.
#' 
#' @param I A vector representing the sources.
#' @param J A vector representing the destinations.
#' @return A matrix of lanes.
get_lanes <- function(I, J) {
  as.matrix(unname(expand.grid(I, J)))
}

#' Generate Bid
#' 
#' Generates a bid.
#' 
#' @param lanes A matrix of lanes.
#' @return A bid.
generate_bid <- function(lanes) {
  L_ <- seq(nrow(lanes))
  w_ <- dnorm(L_, mean = mean(L_))
  w_ <- w_ / sum(w_)
  n  <- sample(L_, 1L, prob = w_)
  sample(L_, n)
}

#' Generate Bids
#' 
#' Generates multiple bids.
#' 
#' @param nB The number of bids.
#' @param L A matrix of lanes.
#' @return A list of bids.
generate_bids <- function(nB, L) {
  B  <- vector(mode = 'list', length = nB)
  for (i in seq(nB)) {
    B[[i]] <- generate_bid(L)
  }
  B
}

#' Simulate Winners
#' 
#' Simulates the winners of bids.
#' 
#' @param nB The number of bids.
#' @param nCS The number of strategic carriers.
#' @return A list of winners.
simulate_winners <- function(nB, nCS) {
  B_ <- matrix(c(seq(nB), rep(NA, nB)), ncol = 2L)
  B_[,2L] <- sample(seq(nCS), nB, replace = TRUE)
  split(B_[,1L], B_[,2L])
}

#' Input Configuration
#' 
#' Configures the input environment.
#' 
#' @param env The environment variable.
#' @param ... Additional arguments.
#' @export
input_configuration <- function(env, ...) {
  args <- list(...)
  if (!is.null(args[["config"]])) {
    for (key in names(config)) {
      assign(key, config[[key]], envir = env)
    }
  } else {
    env$tau <- ifelse(is.null(args[["tau"]]), 12L, args[["tau"]])
    env$nI  <- ifelse(is.null(args[["nI"]]),   2L, args[["nI"]])
    env$nJ  <- ifelse(is.null(args[["nJ"]]),   2L, args[["nJ"]])
    env$nB  <- ifelse(is.null(args[["nB"]]),   3L, args[["nB"]])
    env$nP  <-  1L
    env$nCS <- ifelse(is.null(args[["nCS"]]),   2L, args[["nCS"]])
  }

  env$nL <- env$nI * env$nJ
  env$I_ <- seq(env$nI)
  env$J_ <- seq(env$nJ)
  env$CS <- seq(env$nCS)
}

#' Simulate Auction
#' 
#' Simulates an auction.
#' 
#' @param env The environment variable.
#' @param ... Additional arguments.
#' @export
simulate_auction <- function(env, ...) {
  env$L  <- get_lanes(env$I_, env$J_)
  env$B  <- generate_bids(env$nB, env$L)
  env$winner <- simulate_winners(env$nB, env$nCS)
}

#' Initialize
#' 
#' Initializes the environment.
#' 
#' @param env The environment variable.
#' @param ... Additional arguments.
#' @export
initialize <- function(env, ...) {
  args <- list(...)
  # Storage limits
  if (is.null(rate <- args[["rate"]])) {
    rate <- 10L
  }
  env$R  <- 10L*rate
  # Reset strategic carrier information based on bidding results
  env$CS  <- seq_along(env$winner)
  env$nCS <- length(env$CS)
  # (source, lane) indexing
  env$nLc <- c(0L, unlist(lapply(env$CS, function(idx) length(unlist(env$B[env$winner[[idx]]])))))
  env$L_  <- unlist(lapply(env$CS, function(idx) env$B[env$winner[[idx]]]))
  env$nL_ <- length(env$L_)
  # Transportation cost of strategic carriers
  env$CTb <- rnorm(env$nL_, 12.0, 4.0)
  # Transportation cost of spot carriers
  if (is.null(env$nCO <- args[["nCO"]])) {
    env$nCO <- 1L
  }
  if (is.null(env$CTo <- args[["CTo"]])) {
    env$CTo <- pmax(matrix(rnorm(env$tau * env$nCO * env$nL, 20.0, 10.0), nrow = env$tau), 5.0)
  }
  # Capacity limits
  if (is.null(env$Cb <- args[["Cb"]])) {
    env$Cb <- matrix(sample(env$R%/%rate, env$tau * env$nCS, replace = TRUE), nrow = env$tau)
  }
  if (is.null(env$Co <- args[["Co"]])) {
    env$Co <- matrix(rep((env$R * env$nI) %/% env$nCO, env$tau*env$nCO), ncol = env$nCO)
  }
  # Number of decision variables
  env$nvars <- env$nL_ + ncol(env$CTo)
}

#' Generate CSSAP
#' 
#' Generates a CSSAP (Carrier Scheduling and Shipment Assignment Problem).
#' 
#' @param env The environment variable.
#' @param tau The time horizon. Defaults to 12L.
#' @param ... Additional arguments.
#' @export
generate_cssap <- function(env, tau=12L, ...) {
  do.call(input_configuration, list(env, 'tau'=tau, ...))
  simulate_auction(env)
  do.call(initialize, list(env, ...))
}

#' Print Info
#' 
#' Prints information about the auction.
#' 
#' @param env The environment variable.
#' @param ... Additional arguments.
#' @export
print_info <- function(env, ...) {
  for (i in seq_along(env$winner)) {
    ibids <- unlist(env$winner[[i]])
    cat("Carrier", i, "won bids", paste0(ibids, collapse = ", "), "\n")
    for (bid in env$winner[[i]]) {
      cat(" Bid", bid, "involves lanes", paste0(env$B[[bid]], collapse = ", "), "\n")
    }
  }
  residue <- setdiff(seq_along(env$B), unique(unname(unlist(env$winner))))
  if (length(residue) > 0L) {
    cat("Bids", paste0(residue, collapse = ", "), "are only allocated to spot suppliers.", "\n")
    for (bid in residue) {
      cat(" Bid", bid, "involves lanes", paste0(env$B[[bid]], collapse = ", "), "\n")
    }
  }
}
