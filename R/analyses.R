## analyses
library(ape)
library(phylolm)

## data_loading
dat <- read.csv(file = "data/composed_dataset.csv")
GBOTB <- read.tree("data/GBOTB.tre")
ALLOTB <- read.tree("data/ALLOTB.tre")

## data_preparation
ALLOTBangio <- extract.clade(ALLOTB, "Magnoliophyta")
GBOTBangio <- extract.clade(GBOTB, node = 79883)
tree_fabales <- extract.clade(ALLOTBangio, 
  grep("Fabales", ALLOTBangio$node.label) + length(ALLOTBangio$tip.label))
tree_lamiales <- extract.clade(ALLOTBangio, 
  grep("Lamiales", ALLOTBangio$node.label) + length(ALLOTBangio$tip.label))
dat$herb <- as.numeric(!dat$woody)

aux_sel <- data.frame(
  phyA = dat$spec %in% ALLOTBangio$tip.label,
  phyG = dat$spec %in% GBOTBangio$tip.label,
  perenNE = with(dat, (!annual | is.na(annual) | 
    (annual & woody & !is.na(woody))) & (!epiphyte) | is.na(epiphyte)),
  biome = !is.na(dat$biome)
)
sel <- data.frame(
  default = aux_sel$phyA & aux_sel$perenNE & aux_sel$biome,
  GBOTB = aux_sel$phyG & aux_sel$perenNE & aux_sel$biome,
  X4 = aux_sel$phyA & aux_sel$perenNE & aux_sel$biome & dat$biome == "X4",
  X3 = aux_sel$phyA & aux_sel$perenNE & aux_sel$biome & dat$biome == "X3",
  fabales = dat$spec %in% tree_fabales$tip.label & aux_sel$perenNE,
  lamiales = dat$spec %in% tree_lamiales$tip.label & aux_sel$perenNE,
  annuals = aux_sel$phyA & aux_sel$biome,
  full = aux_sel$phyA & aux_sel$perenNE
)
aux_pgls <- function(vec, dat, phy, sel, biome = TRUE, std = TRUE, ...){
  set.seed(10)
  cols <- c(vec, "herb", "spec")
  if (biome) cols <- c(cols, "biome")
  dat <- dat[sel, cols]
  dat <- dat[rowSums(is.na(dat)) == 0, ]
  # exclusion of biomes with low number of observations (<20) for hebs or woody plants
  form <- as.formula(paste(vec, "~ herb"))
  if (biome) {
    aux_t <- table(dat$biome, dat$herb)
    sel_biomes <- rownames(aux_t)[apply(aux_t > 19, 1, all)]
    dat <- dat[dat$biome %in% sel_biomes, ]
    form <- as.formula(paste(vec, "~ biome * herb"))
    dat$biome <- factor(dat$biome)
    contrasts(dat$biome) <- contr.sum
  }
  if (any(table(dat$herb) < 20)) warning("Low number of observations")
  phy <- keep.tip(phy, dat[, "spec"])
  rownames(dat) <- dat$spec
  data_nst <- dat[, vec]
  arg <- list(formula = form, data = dat, phy = phy, model = "lambda")
  if (all(range(data_nst) == 0:1)) {
    fun <- phyloglm
    arg$model <- NULL
  } else {
    fun <- phylolm
    if (std) arg$data[, vec] <- 
      (arg$data[, vec] - mean(arg$data[, vec])) / sd(arg$data[, vec])
  }
  mod <- do.call(fun, arg)
  mod$varNames <- vec
  if (biome) {
    arg$data$biome <- factor(arg$data$biome) # reset contrasts
    mod$mod_ct <- do.call(fun, arg)
    mod$biomes <- sort(unique(dat$biome))
  }  
  mod$data_nst <- data_nst
  list(mod)
}
vars <- c("freezingFFD", "freezingT", "freezingE", "drought", 
  "droughtP", "fire", "clonal", "shade", "shadeP", "SLA")

## model_runs
mods <- sapply(vars, aux_pgls, dat, ALLOTBangio, sel$default)
# save(mods, file = "data/analyses/models/default.RData")

mods <- sapply(vars, aux_pgls, dat, GBOTBangio, sel$GBOTB)
# save(mods, file = "data/analyses/models/GBOTB.RData")

mods <- sapply(vars, aux_pgls, dat, ALLOTBangio, sel$X4, biome = FALSE)
# save(mods, file = "data/analyses/models/X4.RData")

mods <- sapply(vars, aux_pgls, dat, ALLOTBangio, sel$X3, biome = FALSE)
# save(mods, file = "data/analyses/models/X3.RData")

mods <- sapply(vars, aux_pgls, dat, tree_fabales, sel$fabales, biome = FALSE)
# save(mods, file = "data/analyses/models/fabales.RData")

mods <- sapply(vars, aux_pgls, dat, tree_lamiales, sel$lamiales, biome = FALSE)
# save(mods, file = "data/analyses/models/lamiales.RData")

mods <- sapply(vars, aux_pgls, dat, ALLOTBangio, sel$annuals)
# save(mods, file = "data/analyses/models/annuals.RData")

mods <- sapply(vars, aux_pgls, dat, ALLOTBangio, sel$full, biome = FALSE)
# save(mods, file = "data/analyses/models/full.RData")

## supplementary_analyses
# data_preparation
fabids <- extract.clade(ALLOTBangio, "mrcaott2ott2737")
malvids <- extract.clade(ALLOTBangio, "mrcaott96ott607")
campanulids <- extract.clade(ALLOTBangio, "campanulids")
lamiids <- extract.clade(ALLOTBangio, "lamiids")

# model_runs
mods <- sapply(vars, aux_pgls, dat, fabids, 
  sel$full & dat$spec %in% fabids$tip.label, biome = FALSE)
# save(mods, file = "data/analyses/models/fabids.RData")

mods <- sapply(vars, aux_pgls, dat, malvids, 
  sel$full & dat$spec %in% malvids$tip.label, biome = FALSE)
# save(mods, file = "data/analyses/models/malvids.RData")

mods <- sapply(vars, aux_pgls, dat, campanulids, 
  sel$full & dat$spec %in% campanulids$tip.label, biome = FALSE)
# save(mods, file = "data/analyses/models/campanulids.RData")

mods <- sapply(vars, aux_pgls, dat, lamiids, 
  sel$full & dat$spec %in% lamiids$tip.label, biome = FALSE)
# save(mods, file = "data/analyses/models/lamiids.RData")

#_