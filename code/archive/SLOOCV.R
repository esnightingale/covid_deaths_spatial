################################################################################
# Description: Function to perform spatial LOOCV with INLA fit
# 
# Author: Emily S Nightingale
# Date created: 30/09/2020
# 
################################################################################
################################################################################

## BYM spatial effect, geography-dependent RW
f <- n ~ 
  IMD + prop_minority + log(pop_dens) + prop_kw + 
  f(w, model = "rw1",
    hyper = list(prec = prior.prec.tp),
    values = seq(min(dat$w),max(dat$w)),   
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    hyper = list(prec = prior.prec.tp),
    replicate = geog,
    values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
    scale.model = T) +
  f(la, model = "bym2", graph = g,
    scale.model = T,
    constr = T,
    hyper=list(
      phi =list(param =c(0.5, 2/3)),
      prec = prior.prec.sp) 
  )

preds <- matrix(ncol = 1000, nrow = n_distinct(dat$la))

for (ltla in unique(regions$lad19cd)) {
  
id <- which(regions$lad19cd == ltla)

# Exclude all values for one LTLA
cv_dat <- dat
cv_dat$n[cv_dat$lad19cd == ltla] <- NA 

# Also exclude values for neighbouring LTLAs
nbs <- unlist(g$nbs[id])
cv_dat$n[cv_dat$lad19cd %in% regions$lad19cd[nbs]] <- NA 

# Fit model with excluded values and draw posterior samples
fit <- fit_mod(f, cv_dat)
samples <- inla.posterior.sample(fit, n = 1000)

# Extract predictions from each sample for excluded ltla
preds[id,] <- unlist(lapply(samples1, function(sample) sample$latent[id]))

}

saveRDS(preds, "SLOOCV_preds.rds")

# Rescale and summarise over samples
dat_pred <- bind_cols(dplyr::select(dat, geography, lad19cd, lad19nm, la, wod, E_wk, n), preds) %>%
  mutate(pred_n = exp(value)*E_wk,
         pred_err = abs(n - pred_n))

dat_pred %>%
  group_by(lad19cd, name) %>%
  # sum over weeks
  summarise(value = sum(pred_n)) %>%
  # average over samples
  group_by(lad19cd) %>%
  summarise(q50 = median(value),
            q25 = quantile(value, 0.25),
            q75 = quantile(value, 0.75)) %>% 
  pivot_longer(cols = contains("q")) %>% 
  left_join(regions) %>%
  basic_map(fill = "value", rate1e5 = TRUE) +
  facet_wrap(~name) +
  labs(title = "Predicted deaths per 100,000: Quantiles of 1000 posterior samples per local authority")

