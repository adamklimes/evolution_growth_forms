# Figures
#_________________________________________________________________
# Figure 1 - family representation
library(ape)
library(adephylo)

# loading_data
dat <- read.csv(file = "data/composed_dataset.csv")
taxa <- read.csv("data/family.csv")
ftree <- read.tree("data/Gastauer_2017_familytree.txt")

# data_preparation
taxa$family[taxa$family == "Compositae"] <- "Asteraceae"
taxa$family[taxa$family == "Leguminosae"] <- "Fabaceae" 
taxa$family[taxa$family == "Rhipogonaceae"] <- "Ripogonaceae" 
dat$family <- taxa$family[match(dat$spec, sub("_NA", "", taxa$newtaxa))]
ftree$tip.label <- paste0(toupper(substr(ftree$tip.label,1,1)), 
  substr(ftree$tip.label,2,100))
ftree$tip.label <- sub("_NA", "", ftree$tip.label)

dat_fig_aux <- dat[!is.na(dat$woody) & dat$family %in% ftree$tip.label, 
  c(1,16,6:15)]
dat_fig <- data.frame(dat_fig_aux[, -(1:2)], row.names = dat_fig_aux[, 1])
dat_bin <- data.frame((!is.na(dat_fig)) + 0)
dat_f <- data.frame(apply(dat_bin, 2, tapply, dat_fig_aux$family, 
  function(x) any(x == 1)+0))
miss_f <- matrix(0, ncol = ncol(dat_f), 
  nrow = sum(!ftree$tip.label %in% rownames(dat_f)))
rownames(miss_f) <- ftree$tip.label[!ftree$tip.label %in% rownames(dat_f)]
colnames(miss_f) <- colnames(dat_f)
dat_f <- rbind(dat_f, miss_f)
cols <- data.frame(t(data.frame("white", 
  rev(c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f",
    "#ff7f00","#cab2d6","#6a3d9a")))), 
  row.names = 0:1)
ftree$node.label[ftree$node.label == "lamids"] <- "lamiids"
cap_sel <- ftree$node.label %in% c("ranunculales", "saxifragales", 
  "santalales", "caryophyllales", "ericales")
capitalize <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, 100))
ftree$node.label[cap_sel] <- capitalize(ftree$node.label[cap_sel])
colnames(cols) <- colnames(dat_f)
tips_list <- listTips(ftree)
sel_groups <- c("magnoliids", "monocots", "Ranunculales", "fabids", 
  "malvids", "campanulids", "lamiids", "Saxifragales", "Santalales", 
  "Caryophyllales", "Ericales")
pos_aux <- sapply(sel_groups, function(x, tips_list, rn) 
  which(rn %in% names(tips_list[[x]])), tips_list, rownames(dat_f))
class_f <- rep(NA, nrow(dat_f))
for (i in names(pos_aux)){
  if (any(!is.na(class_f[pos_aux[[i]]]))) stop("Overlapping groups!")
  class_f[pos_aux[[i]]] <- i
}
class_f <- class_f[match(ftree$tip.label, rownames(dat_f))]
labs <- c("Frost-free days", "Minimum temperature", "Freezing exposure", 
  "Drought tolerance", "Drought tolerance (P.)", "Fire tolerance", 
  "Clonality", "Shade tolerance", "Shade tolerance (P.)", "Specific leaf area")
class_f[class_f == "fabids"] <- "Fabidae"
class_f[class_f == "malvids"] <- "Malvidae"
class_f[class_f == "campanulids"] <- "Campanulidae"
class_f[class_f == "lamiids"] <- "Lamiidae"
class_f[class_f == "magnoliids"] <- "Magnoliidae"
class_f[class_f == "monocots"] <- "Liliopsida"

# png("figures/Fig1_family_tree.png", width = 480*10, height = 480*10*1.2, res = 72*10)
par(mai = c(0,0,0,0))
plot(0,0, type = "n", axes = FALSE)
legend(-0.6,-0.75, fill = unlist(cols[2, ]), legend = labs, bty = "n", 
  ncol = 2)
par(new = TRUE, mai = c(1.2,0,0,0))
trait.plot2(compute.brlen(ftree), dat_f, cols, legend = FALSE, 
  class = class_f, margin = 0.65, cex.lab = 0.8, g_lwd = 3)

#_________________________________________________________________
# Figure 2 - transitions
library(phytools)
load(file = "data/analyses/reconstructions/smap_fabids.RData")
load(file = "data/analyses/reconstructions/smap_malvids.RData")
load(file = "data/analyses/reconstructions/smap_campanulids.RData")
load(file = "data/analyses/reconstructions/smap_lamiids.RData")
aux_smap_plot <- function(smap, tit){
  plot(smap, fsize = 0, type = "fan", ftype = "off", 
    colors = setNames(cols, c(1,0)), mar = c(0.1,0.1,2,0.1))
  par(new = TRUE, mai = c(0,0,0,0))
  plot(c(0,1), c(0,1), type = "n", axes = FALSE, ann = FALSE)
  text(0.5, 0.95, tit, cex = 1.5)
  text(0.04, 0.83, paste0("Species:\n", length(smap$tip.label)), adj = c(0,0.5))
  smap_tr <- summary(smap)$Tr
  aux_diagram(c(0.72,0.995,0.72,0.995), smap_tr[1, 2], smap_tr[2, 1], 
    smap$parsimony)
}
aux_diagram <- function(fig, to_w, to_h, min_t){
  fg <- par()$fig
  new_fig <- c(fig[1] * (fg[2] - fg[1]) + fg[1],
    fig[2] * (fg[2] - fg[1]) + fg[1],
    fig[3] * (fg[4] - fg[3]) + fg[3],
    fig[4] * (fg[4] - fg[3]) + fg[3])
  par(new = TRUE, fig = new_fig)
  plot(c(0,1), c(0,1), type = "n", axes = FALSE, ann = FALSE)
  plotrix::draw.circle(0.5, 0.8, radius = 0.2, col = cols[1])
  plotrix::draw.circle(0.5, 0.2, radius = 0.2, col = cols[2])
  xx <- seq(0.1, 0.3, by = 0.01)
  r2 <- 0.16
  lines(xx, (1+sqrt(1-4*(-r2+(xx-0.5)^2+0.25))) / 2, lwd = 2)
  lines(-xx+1, (1+sqrt(1-4*(-r2+(xx-0.5)^2+0.25))) / 2, lwd = 2)
  lines(xx, -(1+sqrt(1-4*(-r2+(xx-0.5)^2+0.25))) / 2+1, lwd = 2)
  lines(-xx+1, -(1+sqrt(1-4*(-r2+(xx-0.5)^2+0.25))) / 2+1, lwd = 2)
  lines(c(0.205,0.3,0.26), c(0.857, 0.8464102, 0.75), lwd = 2)
  lines(-c(0.205,0.3,0.26)+1, -c(0.857, 0.8464102, 0.75)+1, lwd = 2)
  cex_set <- 0.75
  text(0.065, 0.5, to_w, cex = cex_set, adj = c(0.5, 0), srt = 90)
  text(0.935, 0.5, to_h, cex = cex_set, adj = c(0.5, 0), srt = -90)
  text(0.5, 0.5, paste("Min.:", min_t), cex = cex_set)
}
cols <- c("red", "blue")
# cols <- c("#e41a1c", "#377eb8")

# png("figures/Fig2_trans.png", width = 480*10, height = 480*10, res = 72*10)
par(fig = c(0,0.5,0.5,1))
aux_smap_plot(smap_fabids, "Fabidae")
par(new = T, fig = c(0.5,1,0.5,1))
aux_smap_plot(smap_malvids, "Malvidae")
par(new = T, fig = c(0,0.5,0,0.5))
aux_smap_plot(smap_campanulids, "Campanulidae")
par(new = T, fig = c(0.5,1,0,0.5))
aux_smap_plot(smap_lamiids, "Lamiidae")
par(new = TRUE, mai = c(0,0,0,0), fig = c(0,1,0,1))
plot(c(0,1), c(0,1), type = "n", axes = FALSE, ann = FALSE)
legend(0.43, 0.05, fill = cols, legend = c("Woody", "Herbaceous"), 
  bty = "n")

#_________________________________________________________________
# Figure 3 - default model results
library(phylolm)
library(msm)
library(colorspace)

load(file = "data/analyses/models/default.RData")
aux_d <- function(vars, mod){
  cf <- coef(mod)
  pars <- which(names(cf) %in% vars)
  form <- as.formula(paste("~", paste0("x", pars, collapse = " + ")))
  stder <- if(all(vars %in% names(cf))) 
    deltamethod(form, cf, vcov(mod)) else NA
  c(Estimate = sum(cf[pars]), StdError = stder)
}
logit <- function(x) log(x/(1 - x))
inv_logit <- function(x) exp(x) / (1 + exp(x))
get_pars <- function(mod){
  pars_aux <- rbind(
    ME = aux_d("herb", mod),
    Bint = aux_d("herb", mod$mod_ct),
    B1 = c(NA, NA),
    B2 = aux_d(c("herb", "biomeX2:herb"), mod$mod_ct),
    B3 = aux_d(c("herb", "biomeX3:herb"), mod$mod_ct),
    B4 = aux_d(c("herb", "biomeX4:herb"), mod$mod_ct),
    B5 = aux_d(c("herb", "biomeX5:herb"), mod$mod_ct),
    B6 = aux_d(c("herb", "biomeX6:herb"), mod$mod_ct),
    B7 = aux_d(c("herb", "biomeX7:herb"), mod$mod_ct),
    B8 = aux_d(c("herb", "biomeX8:herb"), mod$mod_ct),
    B10 = aux_d(c("herb", "biomeX10:herb"), mod$mod_ct),
    B11 = aux_d(c("herb", "biomeX11:herb"), mod$mod_ct),
    B12 = aux_d(c("herb", "biomeX12:herb"), mod$mod_ct),
    B13 = aux_d(c("herb", "biomeX13:herb"), mod$mod_ct))
  pars_aux[sub("X", "B", mod$biomes[1]), ] <- pars_aux["Bint", ]
  pars_aux <- pars_aux[rownames(pars_aux) != "Bint", ]
  q_mat <- matrix(qt(c(0.025, 0.975), mod$n - mod$p), 
    nrow = nrow(pars_aux), ncol = 2, byrow = TRUE)
  ci <- q_mat * pars_aux[, "StdError"] + pars_aux[, "Estimate"]
  colnames(ci) <- c("lwr", "upr")
  res <- cbind(pars_aux, ci)
  if (all(range(mod$y) == 0:1)){
    cf <- coef(mod)
    eh <- sum(cf[c("(Intercept)", "herb")])
    ew <- cf["(Intercept)"]
    est_n <- (inv_logit(eh) - inv_logit(ew)) / sd(mod$y)
    corr <- est_n / pars_aux["ME", "Estimate"]
    res <- res * corr
  }
  res
}
aux_arrows <- function(pars, ypos, sel_col){
  aux_fn <- function(pars, ypos, sel_col, lwd){
    sel_col <- rep(sel_col, nrow(pars))
    sel_col[sign(pars[, "lwr"]) != sign(pars[, "upr"])] <- 3
    sel_col[sign(pars[, "lwr"]) == 1 & sign(pars[, "upr"]) == 1] <- 
      -sel_col[sign(pars[, "lwr"]) == 1 & sign(pars[, "upr"]) == 1] + 3
    cols_aux <- cols[sel_col]
    if (lwd == 1) cols_aux <- lighten(cols_aux, 0.5)
    arrows(pars[, "lwr"], ypos, x1 = pars[, "upr"], 
      col = cols_aux, length = 0, lwd = lwd) 
    points(pars[, "Estimate"], ypos, col = cols_aux, pch = 16, 
      cex = lwd / 3) 
  }
  pars_s <- pars[!is.na(pars[, "StdError"]) & rownames(pars) != "ME", ]
  ydev <- 0.8
  yres <- 0.3
  aux_seq <- seq(yres, ydev, length.out = nrow(pars_s) / 2) 
  y_coor <- ypos + c(-aux_seq, aux_seq)[1:nrow(pars_s)]
  aux_fn(pars_s, y_coor, sel_col, lwd = 1)
  aux_fn(pars["ME", , drop = FALSE], ypos, sel_col, lwd = 3)
}
nobs.phyloglm <- function(x) x$n

pars <- lapply(mods, get_pars)
cols <- c("blue", "red", "black")
labs <- rev(c("Frost-free days", "Minimum temperature", "Freezing exposure", 
  "Drought tolerance", "Drought tolerance\n(PLANTS database)", "Fire tolerance", 
  "Clonality", "Shade tolerance", "Shade tolerance\n(PLANTS database)", "Specific leaf area"))
hyps <- c("Freezing", "Drought", "Fire", "Shade", "Growth\nrate")

# png("figures/Fig3_def_mod.png", width = 480*10, height = 480*10, res = 72*10)
par(mai = c(0.9,2,0.5,0.1))
plot(c(-1.75,1.6), c(0,27), type = "n", axes = FALSE, ylab = "", 
  xlab = "Standardised effect of herbaceous habit")
axis(1)
box(bty="l")
abline(v = 0, lwd = 2, col = "grey", lty = 2)
axis(2, labels = labs, at = 0:9*3, las = 2, cex.axis = 0.8)
invisible(Map(aux_arrows, pars, 9:0*3, c(1,1,2,2,2,2,2,2,2,2)))
text(-1.83, 9:0*3, sapply(mods, nobs), adj = c(0, 0.5), cex = 0.75)
legend("topright", pch = 16, lwd = 3, col = cols, bty = "n",
  legend = c("In agreement", "In contrast", "No difference"),
  title = "Results vs. hypothesis")
par(new = TRUE, mai = c(0,0,0,0))
plot(c(0,1), c(0,27), type = "n", ann = FALSE, axes = FALSE)
text(0.005, c(-0.5,1.5,3.5,5.5,8) * 2.365 + 3.75, rev(hyps), cex = 1.2, srt = 90)
invisible(sapply(c(0.5,2.5,4.5,6.5) * 2.365 + 3.75, 
  function(x) lines(c(-0.025, 1.025), rep(x, 2))))
#_
