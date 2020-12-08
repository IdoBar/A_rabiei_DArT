pacman::p_load(tidyverse, lme4, lsmeans, gtools)
# set the threshold
high_path_thresh <- 3
cluster_data <- read_csv("https://raw.githubusercontent.com/IdoBar/A_rabiei_DArT/master/data/cluster_data.csv") %>% 
  mutate(Cluster=LETTERS[Cluster], Virulence=ifelse(Path_rating<high_path_thresh,0, 1)) 
data.glmer <- glmer(Virulence ~ Cluster + (1|id), family = "binomial", data = cluster_data, nAGQ = 5)
lsm <- lsmeans(data.glmer, "Cluster", type = "response")
signif_lsm <- pairs(lsm) %>% as.data.frame() %>% mutate(stars=stars.pval(p.value))