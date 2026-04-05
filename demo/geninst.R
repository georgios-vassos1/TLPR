library(TLPR)

env <- generate_instance(
  nI = 2L, nJ = 1L, tau = 4L,
  nB = 5L, nCS = 10L, nCO = 1L, rate = 4.0,
  seed = 42L,
  path = "~/drayage/TLPR/src/instances/instance2x1_4_001.json"
)
