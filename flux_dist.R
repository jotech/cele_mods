library(sybilSBML)
library(sybil)
library(ggplot2)
library(reshape2)

cele_jul <- readSBMLmod("~/uni/dat/mod/Celegans_2016_2.xml")
cele_con <- readSBMLmod("~/uni/dat/mod/Celegans_consensus.xml")
ex_jul <- findExchReact(cele_jul)
ex_con <- findExchReact(cele_con)

sol_jul <- optimizeProb(cele_jul, algorithm = "mtf")
exch_jul <- getFluxDist(sol_jul)
flux_jul <- getNetFlux(exch_jul)
sum(flux_jul@uptake) + sum(flux_jul@product) # active reactions
plot(exch_jul)


sol_con <- optimizeProb(cele_con, algorithm = "mtf")
exch_con <- getFluxDist(sol_con)
flux_con <- getNetFlux(exch_con)
sum(flux_con@uptake) + sum(flux_con@product) # active reactions
plot(exch_con)
