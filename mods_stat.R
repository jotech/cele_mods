library(sybilSBML)
library(sybil)
library(ggplot2)
library(reshape2)

cele_jul <- readSBMLmod("~/uni/dat/mod/Celegans_2016_2.xml")
cele_con <- readSBMLmod("~/uni/dat/mod/Celegans_consensus.xml")

cele_stat <- data.frame()

# basic comparison
ex_jul <- findExchReact(cele_jul)
ex_con <- findExchReact(cele_con)
cele_stat <- rbind(cele_stat, 
                   data.frame("model"="juliane", "reactions"=ncol(S(cele_jul)), 
                              "metabolites"=nrow(S(cele_jul)), "genes"=length(cele_jul@allGenes),
                              "exchanges"=ex_jul@react_num))
cele_stat <- rbind(cele_stat, 
                   data.frame("model"="consensus", "reactions"=ncol(S(cele_con)), 
                              "metabolites"=nrow(S(cele_con)), "genes"=length(cele_con@allGenes),
                              "exchanges"=ex_con@react_num))

cele_melt <- melt(cele_stat, id.vars = "model")

ggplot(cele_melt, aes(x=variable, y=value, fill=model)) + 
  geom_bar(stat="identity", position = position_dodge(width = 0.9)) + xlab("") + ylab("")
#ggsave("~/uni/cele_mods/img/basic_cmp.pdf")

blocked_jul <- sum(blockedReact(cele_jul))
blocked_con <- sum(blockedReact(cele_con))

fva_jul <- fluxVar(cele_jul)
fva_con <- fluxVar(cele_con)
length(which(maxSol(fva_jul, "lp_obj") < 1e-06 & maxSol(fva_jul, "lp_obj") > -1e-06 &
               minSol(fva_jul, "lp_obj") < 1e-06 & minSol(fva_jul, "lp_obj") > -1e-06))
length(which(maxSol(fva_con, "lp_obj") < 1e-06 & maxSol(fva_con, "lp_obj") > -1e-06 &
             minSol(fva_con, "lp_obj") < 1e-06 & minSol(fva_con, "lp_obj") > -1e-06))
