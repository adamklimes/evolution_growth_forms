# extraction_results
library(phylolm)
library(aod)

# functions
aux_extract <- function(type, path){
  load(file = paste0(path, type, ".RData"))
  aux_mod <- function(mod){
    cf <- coef(mod)
    biome_differ <- NULL
    test_h <- coef(summary(mod))["herb", "p.value"]
    est_h <- cf["herb"]    
    if (length(cf) > 2){
      coefs_int <- grep(":herb", names(cf))
      test_int <- wald.test(mod$vcov, cf, coefs_int)$result$chi2["P"]
      if (test_int < 0.05) biome_differ <- "*" else biome_differ <- NULL
    }
    eff <- c("-", "¡", "^")[(est_h < 0) + (est_h > 0)*2 + 1]
    sig <- round(test_h, 3)
    if (sig == 0) sig <- "<0.001"
    paste0(eff, " ", sig, biome_differ)
  }
  sapply(mods, aux_mod)
}
aux_ext_obs <- function(type, path){
  load(file = paste0(path, type, ".RData")) 
  sapply(mods, function(x) {
    tab <- table(x$X[, "herb"])
    paste0(nobs(x), " (", tab[2], ", ", tab[1], ")")
  })
}
aux_ext_lam <- function(type, path){
  load(file = paste0(path, type, ".RData"))  
  sapply(mods, function(x) x$optpar)
}
nobs.phyloglm <- function(x) x$n

# extraction
modtypes <- c("default", "GBOTB", "annuals", "full", "X4", "X3", 
  "fabales", "lamiales")
res_tab <- sapply(modtypes, aux_extract, "data/analyses/models/")
res_tab

# extraction_nobs
res_obs <- sapply(modtypes, aux_ext_obs, "data/analyses/models/")
res_obs

# supplementary_clades
modtypes_clades <- c("fabids", "malvids", "campanulids", "lamiids")
res_tab <- sapply(modtypes_clades, aux_extract, "data/analyses/models/")
res_obs <- sapply(modtypes_clades, aux_ext_obs, "data/analyses/models/")
#_