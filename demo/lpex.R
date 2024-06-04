library(TLPR)

env <- new.env()
generate_cssap(env, rate = 10.0, tau = 12L, nCS = 10L, nCO = 1L, nB = 3L, nI = 2L, nJ = 2L)

env$alpha <- c(rep(20.0, env$nI), rep(20.0, env$nJ), rep(50.0, env$nJ))

exog <- list(
  "Q" = list("func"=rpois, "params"=list('n'=env$tau*env$nI, 'lambda'=10L)),
  "D" = list("func"=rpois, "params"=list('n'=env$tau*env$nJ, 'lambda'=10L))
)

Q <- do.call(exog$Q$func, exog$Q$params)
D <- do.call(exog$Q$func, exog$Q$params)

ccx <- carrier_capacity_padded(env)
tlx <- transition_logic(env, q = Q[seq(env$nI)], d = D[seq(env$nJ)])
slx <- storage_limits(env, q = Q[seq(env$nI)])

obj_ <- c(env$alpha, env$CTb, env$CTo[1L,], env$alpha)

A   <- rbind(ccx$A, tlx$A, slx$A)
rhs <- c(ccx$rhs, tlx$rhs, slx$rhs)
sns <- c(ccx$sense, tlx$sense, slx$sense)

model <- multiperiod_expansion(env, A, obj_, rhs, sns)
model$modelsense <- 'min'
model$vtype <- rep('C', ncol(model$A))

opt <- gurobi::gurobi(model, params = list(OutputFlag = 0L))

phs <- post_hoc_simulation(env, opt$x)
dtx <- copmute_graph_dt(env, phs$S.I, phs$S.J, phs$allocation)
plot2x2instance(env, dtx, phs$S.I, phs$S.J, Q, D)