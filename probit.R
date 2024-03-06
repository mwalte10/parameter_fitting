library(data.table)
library(dplyr)
library(ggplot2)
library(eppasm)
library(readxl)

root <- getwd()
pjnz <- paste0(root, '/spectrum_data/Botswana2023v4 WPP 02_03_2023 KOS.PJNZ')
pop_1 <- paste0(root, '/spectrum_data/Botswana2023v4 WPP 02_03_2023 KOS_pop1.xlsx')
unc <- paste0(root, '/spectrum_data/uncertainty_file_300draws.xlsx')
hiv_births <- readxl::read_xlsx(unc, sheet = "101. Mothers needing PMTCT") %>% data.table()
births <- readxl::read_xlsx(unc, sheet = "21. Total Births" ) %>% data.table()

reshape_draws <- function(x){
  colnames(x) <- as.character(x[2,])
  x <- x[-(1:2),]
  x <- melt(x, id.vars = c('Iteration', 'Sex'))
  return(x)
}

hiv_births <- reshape_draws(hiv_births)
births <- reshape_draws(births)

birth_prev <- merge(hiv_births[,.(Iteration, variable, value)],
                    births[Sex == 'Male+Female',.(Iteration, variable, value)], by = c('Iteration', 'variable'))
birth_prev <- birth_prev[,.(draw = Iteration, year = as.integer(as.character(variable)), prev_est = value.x / value.y)]

#LAF is called frr_scalar in the hivp object "<FRRbyLocation MV>"
##Think I can use ll_ancrtcens

fitdata <- eppasm::prepare_spec_fit(pjnz, proj.end = 2023)
ancsitedat <- rbind(data.table(attr(fitdata$Urban, "eppd")$ancsitedat)[,urban := T],
                    data.table(attr(fitdata$Rural, "eppd")$ancsitedat)[,urban := F])
ancsitedat <- ancsitedat[type == 'ancrt']
ancsitedat <- ancsitedat[,.(prev = weighted.mean(prev,n), n), by = 'year']
ancsitedat <- ancsitedat[,.(prev, n = sum(n)), by = 'year']

dt <- merge(birth_prev[year %in% ancsitedat$year], unique(ancsitedat), by = c('year'), allow.cartesian = T)

##Need to figure out where the err.sigma is coming from
####That's the part that I don't uncertand in Jeff's code where its the variance of the probit transformed prevalence
lm.loss <- function(par, dt){
  laf.par <- par[1]
  err.sigma <- par[2]
  true_prev <- dt$prev
  est.prev <- dt$prev_est

  dt[,W.ancrt := qnorm(prev)]
  dt[,v.ancrt := 2 * pi * exp(W.ancrt^2)*prev*(1-prev)/n]
  dt[,v.ancrt := ifelse(v.ancrt < 0 , 0, v.ancrt)]
  dt[,W.spec := qnorm(prev_est)]
  dt[,laf := W.spec / W.ancrt]


  likelihoods <- dnorm(dt$W.ancrt, mean = dt$W.spec * dt$laf, sd =  sqrt(dt$v.ancrt + err.sigma))

  log.likelihoods <- log(likelihoods)

  out <-  sum(log.likelihoods)

  return(out)
}

par_bad <- c( laf = 8, err.sigma = 3)
par_good <- c( laf = 1, err.sigma = 0.15)
par_true <- c( laf =  1.29153691, err.sigma = 0.05316214)

lm.loss(par_bad, dt)
lm.loss(par_good, dt)
lm.loss(par_true, dt)

parameter.fits <- optim(par = par_good, fn = lm.loss, hessian = F, dt = dt, method = "Nelder-Mead")
hessian <- parameter.fits$hessian
hessian.inv <- solve(hessian)
parameter.se <- sqrt(diag(hessian.inv))

parameter.se
out <- data.table(ml = parameter.fits$par,
                  lb = parameter.fits$par - 1.96 * parameter.se,
                  ub = parameter.fits$par + 1.96 * parameter.se)


rownames(out) <- c('LAF', 'Standard deviation of errors')

##Is the posterior distribution just the rnorm with the new mean and sd?
