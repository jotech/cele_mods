library(sybilSBML)
library(sybil)
library(ggplot2)
library(reshape2)

cele_jul <- readSBMLmod("~/uni/dat/mod/Celegans_2016_2.xml")
cele_con <- readSBMLmod("~/uni/dat/mod/Celegans_consensus.xml")
ex_jul <- findExchReact(cele_jul)
ex_con <- findExchReact(cele_con)

sol  <- optimizeProb(cele_jul) 
exch <- getFluxDist(sol, ex_jul)
flux_jul <- getNetFlux(exch)
data.frame("uptake"=cele_jul@met_name[ex_jul@met_pos[flux_jul@uptake]], "value"=flux_jul@rate[flux_jul@uptake])
data.frame("product"=cele_jul@met_name[ex_jul@met_pos[flux_jul@product]], "value"=flux_jul@rate[flux_jul@product])


sol  <- optimizeProb(cele_con) 
exch <- getFluxDist(sol, ex_con)
flux_con <- getNetFlux(exch)
data.frame("uptake"=cele_con@met_name[ex_con@met_pos[flux_con@uptake]], "value"=flux_con@rate[flux_con@uptake])
data.frame("product"=cele_con@met_name[ex_con@met_pos[flux_con@product]], "value"=flux_con@rate[flux_con@product])
