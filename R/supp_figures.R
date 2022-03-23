# Supplementary figures

# data_loading
dat <- read.csv(file = "data/composed_dataset.csv")

#___________________________________________________________________
# Figure S1 - proportion of woody plants for fire tolerance
aux_fire <- with(dat, table(woody[!is.na(fire)], biome[!is.na(fire)]))
aux_fire <- aux_fire[, apply(aux_fire > 19, 2, all)]
fire_prop <- cbind(table(dat$woody[!is.na(dat$fire)])/sum(table(dat$woody[!is.na(dat$fire)])), aux_fire/rep(colSums(aux_fire), each = 2))
aux_woody <- table(dat$woody, dat$biome)
aux_woody <- aux_woody[, colnames(aux_woody) %in% colnames(aux_fire)]
woody_prop <- cbind(table(dat$woody)/sum(table(dat$woody)), aux_woody/rep(colSums(aux_woody), each = 2))

# png("figures/S1_fire_woodyprop.png", width = 480*10, height = 480*10, res = 72*10)
par(mai = c(2.1,0.8,0.4,0.4))
barplot(woody_prop[2:1, ], width = 0.5, space = c(0, rep(1.5,5)), col=c("#de2d26", "#fc9272"), ylab = "Proportion of woody species", names.arg = rep("",6))
barplot(fire_prop[2:1, ], add = TRUE, width = 0.5, space = c(1, rep(1.5,5)), names.arg = rep("",6))
legend(3.65,0.9, fill = c("#de2d26", "#fc9272", grey.colors(2)), legend = c("Woody dataset", "Herbaceous dataset", "Woody fire-subset", "Herbaceous fire-subset"), bg = "white")
aux_lab <- c("Without biomes", "Tropical and subtropical\nmoist broadleaf forests", "Tropical and subtropical\ndry broadleaf forests", "Tropical and subtropical\nconiferous forests", "Temperate broadleaf\nand mixed forests", "Temperate\nconiferous forests")
axis(1, labels = aux_lab, at = 1:6*1.25-0.75, tick = FALSE, las = 2)
axis(3, labels = colSums(cbind(table(dat$woody[!is.na(dat$fire)]), aux_fire)), at = 1:6*1.25-0.75, tick = FALSE, line = -0.5)

#___________________________________________________________________
# Figure S2 - results of the default model of shade tolerance
library(phylolm)
library(msm)
load(file = "data/analyses/models/default.RData")
mod <- mods$shade$mod_ct
aux_d <- function(vars, mod){
  cf <- coef(mod)
  pars <- c(1, which(names(cf) %in% vars))
  form <- as.formula(paste("~", paste0("x", pars, collapse = " + ")))
  c(Estimate = sum(cf[pars]), StdError = deltamethod(form, cf, vcov(mod)))
}
aux_ci <- function(pars, mod){
  qt(c(0.025, 0.975), mod$n - 2) * pars["StdError"] + pars["Estimate"] 
}
pars_clon <- rbind(B1 = aux_d(NULL, mod),
  B1w = aux_d("herb", mod),
  B2 = aux_d("biomeX2", mod),
  B2w = aux_d(c("herb", "biomeX2", "biomeX2:herb"), mod),
  B3 = aux_d("biomeX3", mod),
  B3w = aux_d(c("herb", "biomeX3", "biomeX3:herb"), mod),
  B4 = aux_d("biomeX4", mod),
  B4w = aux_d(c("herb", "biomeX4", "biomeX4:herb"), mod),
  B5 = aux_d("biomeX5", mod),
  B5w = aux_d(c("herb", "biomeX5", "biomeX5:herb"), mod)
)
pars <- cbind(pars_clon, t(apply(pars_clon, 1, aux_ci, mod)))
cols <- c("red", "blue")
nobs <- colSums(mod$X)
nobs[1] <- nobs[1] - sum(nobs[2:5]) 

# png("figures/S2_shade_tol.png", width = 480*10, height = 480*10, res = 72*10)
par(mai = c(1.7,0.8,0.4,0.4))
plot(1:nrow(pars), pars[, "Estimate"], axes = FALSE, xlab = "", 
  ylab = "Standardised shade tolerance",
  ylim = c(min(pars[, 3]), max(pars[, 4])),
  col = cols, pch = 16)
box(bty = "l")
axis(2)
arrows(x0 = 1:nrow(pars), y0 = pars[, 3], y1 = pars[, 4], 
  code = 3, length = 0.1, angle = 90, col = cols)
abline(v = 1:4*2+0.5, lwd = 2, lty = 2, col = "grey")
lab <- c("Tropical and subtropical\nmoist broadleaf forests", "Tropical and subtropical\ndry broadleaf forests", "Tropical and subtropical\nconiferous forests", "Temperate broadleaf\nand mixed forests", "Temperate\nconiferous forests")
axis(1, labels = lab, las = 2, at = 1:5*2-0.5, cex.axis = 0.8)
legend(7.5,-0.4, fill = cols, legend = c("Woody", "Herbaceous"), bg = "white")
text(1:5*2-0.5, -0.67, nobs[1:5], cex = 0.8)

#___________________________________________________________________
# Figure S3 - results of the default model of clonality
library(phylolm)
library(msm)
load(file = "data/analyses/models/default.RData")
mod <- mods$clonal$mod_ct
aux_d <- function(vars, mod){
  cf <- coef(mod)
  pars <- c(1, which(names(cf) %in% vars))
  form <- as.formula(paste("~", paste0("x", pars, collapse = " + ")))
  c(Estimate = sum(cf[pars]), StdError = deltamethod(form, cf, vcov(mod)))
}
aux_ci <- function(pars, mod){
  qt(c(0.025, 0.975), mod$n - 2) * pars["StdError"] + pars["Estimate"] 
}
inv_logit <- function(x) exp(x) / (1 + exp(x))
pars_clon <- rbind(B1 = aux_d(NULL, mod),
  B1w = aux_d("herb", mod),
  B2 = aux_d("biomeX2", mod),
  B2w = aux_d(c("herb", "biomeX2", "biomeX2:herb"), mod),
  B3 = aux_d("biomeX3", mod),
  B3w = aux_d(c("herb", "biomeX3", "biomeX3:herb"), mod),
  B4 = aux_d("biomeX4", mod),
  B4w = aux_d(c("herb", "biomeX4", "biomeX4:herb"), mod),
  B5 = aux_d("biomeX5", mod),
  B5w = aux_d(c("herb", "biomeX5", "biomeX5:herb"), mod)
)
pars <- inv_logit(cbind(pars_clon, t(apply(pars_clon, 1, aux_ci, mod))))
cols <- c("red", "blue")
nobs <- colSums(mod$X)
nobs[1] <- nobs[1] - sum(nobs[2:5]) 

# png("figures/S3_clon_res.png", width = 480*10, height = 480*10, res = 72*10)
par(mai = c(1.7,0.8,0.4,0.4))
plot(1:nrow(pars), pars[, "Estimate"], axes = FALSE, xlab = "", 
  ylab = "Probability of being clonal",
  ylim = c(0, 1),
  col = cols, pch = 16, yaxs = "i")
box(bty = "c")
axis(2)
arrows(x0 = 1:nrow(pars), y0 = pars[, 3], y1 = pars[, 4], 
  code = 3, length = 0.1, angle = 90, col = cols)
abline(v = 1:4*2+0.5, lwd = 2, lty = 2, col = "grey")
lab <- c("Tropical and subtropical\nmoist broadleaf forests", "Tropical and subtropical\ndry broadleaf forests", "Tropical and subtropical\nconiferous forests", "Temperate broadleaf\nand mixed forests", "Temperate\nconiferous forests")
axis(1, labels = lab, las = 2, at = 1:5*2-0.5, cex.axis = 0.8)
legend(7.5,0.3, fill = cols, legend = c("Woody", "Herbaceous"), bg = "white")
text(1:5*2-0.5, 0.05, nobs[1:5], cex = 0.8)
#_