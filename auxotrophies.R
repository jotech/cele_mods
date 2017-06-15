library(sybilSBML)
library(sybil)
library(ggplot2)
library(reshape2)

cele_jul <- readSBMLmod("~/uni/dat/mod/Celegans_2016_2.xml")
cele_con <- readSBMLmod("~/uni/dat/mod/Celegans_consensus.xml")
ex_jul <- findExchReact(cele_jul)
ex_con <- findExchReact(cele_con)

key = "b12"
pos <- grep(key, cele_jul@met_name, ignore.case = T)
data.frame("name"=cele_jul@met_name[pos], "id"=cele_jul@met_id[pos])

pos <- grep(key, cele_con@met_name, ignore.case = T)
data.frame("name"=cele_con@met_name[pos], "id"=cele_con@met_id[pos])


aux_db <- read.csv("~/uni/cele_mods/dat/auxotrophy.csv", stringsAsFactors = F)
aux_df <- data.frame()
for(i in 1:nrow(aux_db)){
  met     <- aux_db$Metabolite[i]
  id_jul  <- aux_db$Juliane[i]
  id_con  <- aux_db$Consensus[i]
  aux_jul <- check_aux(cele_jul, met_id = id_jul)
  aux_con <- check_aux(cele_con, met_id = id_con)
  aux_df  <- rbind(aux_df, data.frame("met"=met, "aux_jul"=aux_jul, "aux_con"=aux_con))  
}

check_aux <- function(mod, met_id){
  ex <- findExchReact(mod)
  pos <- grep(met_id, ex@met_id, fixed = T)
  if(length(pos)==0){
    mod <- addExchReact(mod, met_id, lb=0, ub=-1000, tag="EX_")
    ex <- findExchReact(mod)
    new_obj <- ex@react_id[which(ex@met_id==met_id)]
  }else{
    new_obj <- ex@react_id[pos]
    mod_jul <- changeBounds(mod, new_obj, lb = 0, ub = 1000)
  }
  mod <- changeObjFunc(mod, react = new_obj)
  sol <- optimizeProb(mod, retOptSol=F)
  if(sol$stat == 5 & sol$obj>SYBIL_SETTINGS("TOLERANCE")){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

aux_melt <- melt(aux_df, id.vars = "met", variable.name = "mod")
aux_melt$value <- ifelse(aux_melt$value==TRUE, 0,1)
aux_melt$met <- factor(aux_melt$met, levels = levels(aux_melt$met)[order(levels(aux_melt$met), decreasing = T)])
ggplot(aux_melt) + geom_tile(aes(x=mod, y=met, fill=value)) +
  scale_fill_continuous(low = "white", high = "black") + 
  ggtitle("Auxotrophic substances which could be produced") + theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(size=11, hjust = 1, colour="black"),
        axis.text.y = element_text(size=11, hjust = 1, colour="black"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="none")
#ggsave("~/uni/cele_mods/img/auxotrophy.pdf")
