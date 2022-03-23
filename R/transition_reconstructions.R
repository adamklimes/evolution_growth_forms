## transition_reconstructions
library(phytools)
library(phangorn)

## data_loading
dat <- read.csv(file = "data/composed_dataset.csv")
ALLOTB <- read.tree("data/ALLOTB.tre")

## data_preparation
ALLOTBangio <- extract.clade(ALLOTB, "Magnoliophyta")
woody <- dat$woody[dat$spec %in% ALLOTBangio$tip.label & !is.na(dat$woody)]
names(woody) <- dat$spec[dat$spec %in% ALLOTBangio$tip.label & !is.na(dat$woody)]
fabids <- extract.clade(ALLOTB, "mrcaott2ott2737")
malvids <- extract.clade(ALLOTB, "mrcaott96ott607")
campanulids <- extract.clade(ALLOTB, "campanulids")
lamiids <- extract.clade(ALLOTB, "lamiids")

## analyses
aux_simmap <- function(phy, trait, nsim = 1){
  trait <- trait[names(trait) %in% phy$tip.label]
  phy <- keep.tip(phy, names(trait))
  set.seed(10)
  make.simmap(phy, trait, model = "ARD", nsim = nsim)
}
aux_parsimony <- function(phy, trait){
  trait <- trait[names(trait) %in% phy$tip.label]
  phy <- keep.tip(phy, names(trait))
  X <- phyDat(as.matrix(trait), type = "USER", levels = 0:1)
  parsimony(phy, X)
}

smap_fabids <- aux_simmap(fabids, woody)
smap_fabids$parsimony <- aux_parsimony(fabids, woody)
# save(smap_fabids, file = "data/analyses/reconstructions/smap_fabids.RData")
smap_malvids <- aux_simmap(malvids, woody)
smap_malvids$parsimony <- aux_parsimony(malvids, woody)
# save(smap_malvids, file = "data/analyses/reconstructions/smap_malvids.RData")
smap_campanulids <- aux_simmap(campanulids, woody)
smap_campanulids$parsimony <- aux_parsimony(campanulids, woody)
# save(smap_campanulids, file = "data/analyses/reconstructions/smap_campanulids.RData")
smap_lamiids <- aux_simmap(lamiids, woody)
smap_lamiids$parsimony <- aux_parsimony(lamiids, woody)
# save(smap_lamiids, file = "data/analyses/reconstructions/smap_lamiids.RData")
load(file = "data/analyses/reconstructions/smap_fabids.RData")
load(file = "data/analyses/reconstructions/smap_malvids.RData")
load(file = "data/analyses/reconstructions/smap_campanulids.RData")
load(file = "data/analyses/reconstructions/smap_lamiids.RData")

## calculate_changes
calc_changes <- function(smap){
  calc_tip_changes <- function(tip, smap){
    np <- nodepath(smap, tip, length(smap$tip.label) + 1)
    changes <- NULL
    for (i in seq_along(np)[-1]){
      ed <- which(rownames(smap$mapped.edge) == 
        paste(np[i], np[i-1], sep = ","))
      changes <- c(changes, rev(smap$maps[[ed]]))
    }    
    changes
  }
  res <- lapply(seq_along(smap$tip.label), calc_tip_changes, smap)
  sapply(res, function(x) sum(abs(diff(as.numeric(factor(names(x)))))))
}

trans_fabids <- calc_changes(smap_fabids)
trans_malvids <- calc_changes(smap_malvids)
trans_campanulids <- calc_changes(smap_campanulids)
trans_lamiids <- calc_changes(smap_lamiids)

sum_changes <- function(trans, smap, woody){
  res <- table(trans, woody = woody[match(smap$tip.label, names(woody))])
  list(trans = res, prim_woody = res[1, 2] / sum(res[, 2]))
}
sum_changes(trans_fabids, smap_fabids, woody)
sum_changes(trans_malvids, smap_malvids, woody)
sum_changes(trans_campanulids, smap_campanulids, woody)
sum_changes(trans_lamiids, smap_lamiids, woody)
#_