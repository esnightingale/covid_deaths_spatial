set_priors <- function(){
  
  # Fit fixed-effect glm to obtain reasonable bound for random effect precision
  nullmod <- MASS::glm.nb(n ~ IMD_quint + prop_minority + offset(log(E_wk)), dat)
  summary(nullmod)
  
  res <- residuals(nullmod)
  sd(res)
  
  # Calculare SD of residuals averaged over random effect index
  dat %>%
    filter(!is.na(n)) %>%
    mutate(res = res) %>%
    group_by(wk_since_first) %>%
    summarise(res = mean(res)) %>%
    pull(res) %>%
    var(na.rm = TRUE) %>%
    sqrt() -> sd_time1
  
  dat %>%
    filter(!is.na(n)) %>%
    mutate(res = res) %>%
    group_by(w) %>%
    summarise(res = mean(res)) %>%
    pull(res) %>%
    var(na.rm = TRUE) %>%
    sqrt() -> sd_time2
  
  dat %>%
    filter(!is.na(n)) %>%
    mutate(res = res) %>%
    group_by(lad19cd) %>%
    summarise(res = mean(res)) %>%
    pull(res) %>%
    var(na.rm = TRUE) %>%
    sqrt() -> sd_space
  
  # Define pc priors P(sig > U) = alpha, setting U based on rule of thumb from Simpson (2017)
  prior.prec.tp1 <- list(prec = list(prior = "pc.prec",
                                     param = c(sd_time1/0.31, 0.01)))
  prior.prec.tp2 <- list(prec = list(prior = "pc.prec",
                                     param = c(sd_time2/0.31, 0.01)))
  prior.prec.sp <- list(prec = list(prior = "pc.prec",
                                    param = c(sd_space/0.31, 0.01)))
  
  # Sample and plot priors on SD scale
  prior.samp.tp1 <- inla.pc.rprec(10000, u = sd_time1/0.31, alpha = 0.01)
  hist(1/sqrt(prior.samp.tp1), breaks = 100)
  
  prior.samp.tp2 <- inla.pc.rprec(10000, u = sd_time2/0.31, alpha = 0.01)
  hist(1/sqrt(prior.samp.tp2), breaks = 100)
  
  prior.samp.sp <- inla.pc.rprec(10000,u = sd_space/0.31, alpha = 0.01)
  hist(1/sqrt(prior.samp.sp), breaks = 100)
  
  return(list(time1 = prior.prec.tp1, time2 = prior.prec.tp2, space = prior.prec.sp))
}
