#' Get Lanes
#' 
#' Generates a matrix of lanes.
#' 
#' @param I A vector representing the sources.
#' @param J A vector representing the destinations.
#' @return A matrix of lanes.
get_lanes <- function(I, J) {
  unname(as.matrix(data.table::CJ(J, I)))[,c(2L,1L),drop = F] # reverse order of columns
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
    cfg <- args[["config"]]
    for (key in names(cfg)) {
      assign(key, cfg[[key]], envir = env)
    }
  } else {
    env$tau <- ifelse(is.null(args[["tau"]]), 12L, args[["tau"]])
    env$nI  <- ifelse(is.null(args[["nI"]]),   2L, args[["nI"]])
    env$nJ  <- ifelse(is.null(args[["nJ"]]),   2L, args[["nJ"]])
    env$nL  <- env$nI * env$nJ
    ncombs  <- sum(choose(env$nL,seq(env$nL)))
    env$nB  <- min(ifelse(is.null(args[["nB"]]), sample(ncombs, 1L), args[["nB"]]), ncombs)
    env$nP  <-  1L
    env$nCS <- ifelse(is.null(args[["nCS"]]), 2L, args[["nCS"]])
  }

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

#' Initialize Environment
#'
#' Sets up the drayage environment: storage limits, carrier information, cost
#' parameters, capacity limits, and route indexing.
#'
#' @param env The environment variable.
#' @param ... Additional arguments.
#' @export
init_env <- function(env, ...) {
  args <- list(...)
  # Storage limits
  if (is.null(rate <- args[["rate"]])) {
    rate <- 10L
  }
  env$R  <- 10L*rate
  # Minimum transportation cost
  if (is.null(min.TC <- args[["min.TC"]])) {
    min.TC <- 1.0
  }
  # Reset strategic carrier information based on bidding results
  env$CS  <- seq_along(env$winner)
  env$nCS <- length(env$CS)
  # (source, lane) indexing
  env$nLc <- c(0L, unlist(lapply(env$CS, function(idx) length(unlist(env$B[env$winner[[idx]]])))))
  env$L_  <- unlist(lapply(env$CS, function(idx) env$B[env$winner[[idx]]]))
  env$nL_ <- length(env$L_)
  # Transportation cost of strategic carriers
  env$CTb <- pmax(rnorm(env$nL_, 12.0, 4.0), min.TC)
  # Transportation cost of spot carriers
  if (is.null(env$nCO <- args[["nCO"]])) {
    env$nCO <- 1L
  }
  # Carrier indexing on lanes (to be bind with lanes)
  env$carriers <- c(rep(seq_along(env$winner), env$nLc[-1L]), rep(env$nCS+seq(env$nCO), each = env$nL))
  # Transportation cost of spot carriers
  if (is.null(env$CTo <- args[["CTo"]])) {
    env$CTo <- pmax(matrix(rnorm(env$tau * env$nCO * env$nL, 20.0, 10.0), nrow = env$tau), 5.0)
  }
  # Capacity limits
  if (is.null(env$Cb <- args[["Cb"]])) {
    env$Cb <- matrix(sample(rate, env$tau * env$nCS, replace = TRUE), nrow = env$tau)
  }
  if (is.null(env$Co <- args[["Co"]])) {
    env$Co <- matrix(rep((rate * env$nI) %/% env$nCO, env$tau*env$nCO), ncol = env$nCO)
  }
  # Number of decision variables
  env$nvars <- env$nL_ + ncol(env$CTo)
}

#' Get From Routes
#' 
#' Retrieves the routes originating from different sources.
#' 
#' @param env The environment variable.
get_from_routes <- function(env, ...) {
  env$from_i <- vector(mode = 'list', length = env$nI)
  for (i in seq(env$nI)) {
    # All possible routes from origin i to the destinations
    idx <- (seq(env$nJ) - 1L) * env$nI
    # All indices in the bids that start from origin i
    msk <- which(apply(outer(env$L_, idx + i, '=='), 1L, any))
    # Adding spot indices that start from origin i
    idx <- env$nL_ + c(outer(idx, (seq(env$nCO) - 1L)*env$nL, '+')) + i
    # Store result
    env$from_i[[i]] <- c(msk, idx)
  }
}

#' Get To Routes
#' 
#' Retrieves the routes leading to different destinations.
#' 
#' @param env The environment variable.
get_to_routes <- function(env, ...) {
  env$to_j <- vector(mode = 'list', length = env$nJ)
  for (j in seq(env$nJ)) {
    # All possible routes to destination j from all origins
    idx <- (j-1L)*env$nI + seq(env$nI)
    # All indices in the bids that do to destination j
    msk <- which(apply(outer(env$L_, idx, '=='), 1L, any))
    # Adding spot indices that go to destination j
    idx <- env$nL_ + c(outer(idx, (seq(env$nCO)-1L)*env$nL, '+'))
    # Store result
    env$to_j[[j]] <- c(msk, idx)
  }
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
  do.call(init_env, list(env, ...))
  get_from_routes(env)
  get_to_routes(env)
}

# ── Regime preset table ───────────────────────────────────────────────────────
# Returns a list of defaults for Q, D, W, alpha, cap_factor, seasonality
# scaled to the given storage limit R.  All values can be overridden by
# the caller via explicit arguments to generate_instance.
.regime_defaults <- function(regime, R, nI, nJ) {
  r2 <- R %/% 2L
  r4 <- R %/% 4L
  switch(regime,
    balanced = list(
      Q           = list(vals = c(0L, r2, R),           prob = c(0.10, 0.70, 0.20)),
      D           = list(vals = c(0L, r2, R),           prob = c(0.10, 0.50, 0.40)),
      W           = list(vals = c(10.0, 15.0, 22.0),    prob = NULL),
      alpha       = c(rep(5.0, nI), rep(c(4.0, 8.0),  each = nJ)),
      cap_factor  = 1.0,
      seasonality = 0.0
    ),
    # Scarce supply, high demand, expensive spot market, tight carrier capacity.
    tight = list(
      Q           = list(vals = c(0L, r4, r2),              prob = c(0.20, 0.50, 0.30)),
      D           = list(vals = c(r2, R, as.integer(R*1.5)), prob = c(0.20, 0.50, 0.30)),
      W           = list(vals = c(15.0, 25.0, 40.0),         prob = NULL),
      alpha       = c(rep(8.0, nI), rep(c(6.0, 12.0), each = nJ)),
      cap_factor  = 0.5,
      seasonality = 0.0
    ),
    # Abundant supply, low demand, cheap spot market, ample carrier capacity.
    loose = list(
      Q           = list(vals = c(r2, R, as.integer(R*1.5)), prob = c(0.30, 0.50, 0.20)),
      D           = list(vals = c(0L, r4, r2),               prob = c(0.30, 0.50, 0.20)),
      W           = list(vals = c(8.0, 12.0, 16.0),          prob = NULL),
      alpha       = c(rep(3.0, nI), rep(c(2.0, 4.0),  each = nJ)),
      cap_factor  = 2.0,
      seasonality = 0.0
    ),
    # Wide Q/D/W distributions and seasonal spot rate pattern.
    volatile = list(
      Q           = list(vals = c(0L, r4, r2, R),              prob = c(0.15, 0.30, 0.35, 0.20)),
      D           = list(vals = c(0L, r4, r2, R),              prob = c(0.10, 0.20, 0.40, 0.30)),
      W           = list(vals = c(5.0, 10.0, 18.0, 30.0, 50.0), prob = NULL),
      alpha       = c(rep(5.0, nI), rep(c(4.0, 8.0),  each = nJ)),
      cap_factor  = 1.0,
      seasonality = 0.3
    ),
    stop("Unknown regime: ", regime,
         ". Choose from: balanced, tight, loose, volatile.")
  )
}

#' Generate Instance
#'
#' Generates a complete TLPR problem instance, including auction structure,
#' routing, stochastic parameters, and derived indices. Optionally writes
#' the instance to a JSON file.
#'
#' @param nI Number of inbound (source) locations. Defaults to 1L.
#' @param nJ Number of outbound (destination) locations. Defaults to 1L.
#' @param tau Time horizon. Defaults to 4L.
#' @param nB Number of bids. Defaults to 5L.
#' @param nCS Number of strategic carriers. Defaults to 10L.
#' @param nCO Number of spot (outside option) carriers. Defaults to 1L.
#' @param rate Storage-limit rate: \code{R = 10 * rate}. Defaults to 1.0.
#' @param regime One of \code{"balanced"} (default), \code{"tight"},
#'   \code{"loose"}, or \code{"volatile"}. Sets default values for Q, D, W,
#'   alpha, cap_factor, and seasonality, all scaled to R.  Any of those
#'   arguments supplied explicitly will override the regime preset.
#' @param Q List with \code{vals} (supply support, scaled to R by default) and
#'   \code{prob}. \code{NULL} uses the regime preset.
#' @param D List with \code{vals} (demand support) and \code{prob}.
#'   \code{NULL} uses the regime preset.
#' @param W List with \code{vals} (spot-rate support, actual cost units) and
#'   \code{prob} (uniform if \code{NULL}). \code{NULL} uses the regime preset.
#' @param alpha Holding-cost penalty vector (length \code{nI + 2*nJ}).
#'   \code{NULL} uses the regime preset.
#' @param cap_factor Multiplier on carrier capacity: \code{< 1} tightens,
#'   \code{> 1} loosens. \code{NULL} uses the regime preset.
#' @param seasonality Amplitude of the sinusoidal seasonal multiplier applied
#'   to spot rates (\code{CTo}) over time: \code{0} = no seasonality.
#'   \code{NULL} uses the regime preset.
#' @param max.S Maximum storage level. \code{NULL} defaults to \code{R}.
#' @param seed Integer seed for reproducibility.
#' @param path If non-NULL, writes the instance as JSON.  \code{TRUE} uses a
#'   default filename \code{instance{nI}x{nJ}_{tau}_001.json}.
#' @return The populated environment (invisibly).
#' @export
generate_instance <- function(
    nI = 1L, nJ = 1L, tau = 4L,
    nB = 5L, nCS = 10L, nCO = 1L, rate = 1.0,
    regime      = c("balanced", "tight", "loose", "volatile"),
    Q           = NULL,
    D           = NULL,
    W           = NULL,
    alpha       = NULL,
    cap_factor  = NULL,
    seasonality = NULL,
    max.S       = NULL,
    seed        = NULL,
    path        = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  regime <- match.arg(regime)
  R      <- as.integer(round(10.0 * rate))

  # ── Resolve regime presets; explicit args override ─────────────────────────
  defs <- .regime_defaults(regime, R, nI, nJ)
  if (is.null(Q))           Q           <- defs$Q
  if (is.null(D))           D           <- defs$D
  if (is.null(W))           W           <- defs$W
  if (is.null(alpha))       alpha       <- defs$alpha
  if (is.null(cap_factor))  cap_factor  <- defs$cap_factor
  if (is.null(seasonality)) seasonality <- defs$seasonality

  # Normalise W probabilities (uniform if prob = NULL)
  if (is.null(W$prob)) {
    W$prob <- rep(1.0 / length(W$vals), length(W$vals))
  }
  W$prob <- W$prob / sum(W$prob)

  # ── Auction structure + routing ────────────────────────────────────────────
  env <- new.env(parent = baseenv())
  generate_cssap(env, tau = tau, nI = nI, nJ = nJ,
                 nB = nB, nCS = nCS, nCO = nCO, rate = rate)

  # ── Carrier post-processing ────────────────────────────────────────────────
  env$carrierIdx <- setNames(as.list(seq_along(env$winner)), names(env$winner))
  env$winnerKey  <- names(env$winner)
  df             <- data.frame(
    Carrier = names(env$winner)[env$carriers[seq(env$nL_)]],
    CTb     = env$CTb
  )
  env$CTb_list <- split(df$CTb, df$Carrier)

  # ── Carrier capacity with cap_factor ──────────────────────────────────────
  # Cb: per-carrier strategic capacity, R × nCS.  Scaled by cap_factor.
  env$Cb <- matrix(
    pmax(1L, as.integer(round(
      matrix(sample(R, tau * env$nCS, replace = TRUE), nrow = tau) * cap_factor
    ))),
    nrow = tau
  )
  # Co: spot-carrier capacity per period (constant over time).
  env$Co <- matrix(
    rep(pmax(1L, as.integer(R * nI * cap_factor / nCO)), tau * nCO),
    ncol = nCO
  )

  # ── Spot rates (CTo) with optional seasonality ────────────────────────────
  # CTo initialises the LP warm-start for each period; W$vals are the realized
  # rates used in the DP.  Seasonality adds a sinusoidal multiplier over time.
  base_spot <- pmax(
    matrix(rnorm(tau * nCO * env$nL, 20.0, 10.0), nrow = tau),
    5.0
  )
  if (seasonality > 0) {
    season_mult <- 1.0 + seasonality * sin(2.0 * pi * seq(0L, tau - 1L) / tau)
    base_spot   <- sweep(base_spot, 1L, season_mult, `*`)
  }
  env$CTo <- base_spot

  # ── Stochastic parameters ──────────────────────────────────────────────────
  env$max.S  <- if (is.null(max.S)) R else as.integer(max.S)
  env$nSI    <- env$max.S + 1L
  env$nSJ    <- 2L * env$max.S + 1L
  env$nSdx   <- (env$nSI ^ nI) * (env$nSJ ^ nJ)
  env$nA     <- env$max.S + 1L
  env$alpha  <- alpha
  env$Q      <- Q
  env$D      <- D
  env$W      <- W
  env$nQ     <- length(Q$vals)
  env$nD     <- length(D$vals)
  env$nW     <- length(W$vals)
  env$nOmega <- (env$nQ ^ nI) * (env$nD ^ nJ) * (env$nW ^ nCO)

  # ── Mixed-radix keys ──────────────────────────────────────────────────────
  env$stateKeys <- get_adjustment_weights(env)
  env$flowKeys  <- c(
    env$nQ ^ (seq(nI) - 1L),
    (env$nQ ^ nI) * env$nD ^ (seq(nJ) - 1L),
    (env$nQ ^ nI) * (env$nD ^ nJ) * env$nW ^ (seq(nCO) - 1L)
  )

  # ── Optional JSON output ───────────────────────────────────────────────────
  if (!is.null(path)) {
    if (isTRUE(path)) {
      path <- sprintf("instance%dx%d_%d_001.json", nI, nJ, tau)
    }
    json <- jsonlite::toJSON(
      sapply(sort(names(env)), get, envir = env, simplify = FALSE),
      pretty = TRUE
    )
    writeLines(json, path)
    message("Instance written to: ", path)
  }

  invisible(env)
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
