library(data.table)
library(dplyr)
library(ggplot2)
library(eppasm)
library(readxl)

root <- getwd()

pjnz <- paste0(root, '/spectrum_data/Botswana2023v4 WPP 02_03_2023 KOS.PJNZ')
birth_data <- fread(paste0(root, '/spectrum_data/birth_draws_laf1.csv'))
ancsitedat <- fread(paste0(root, '/spectrum_data/ancrt_data.csv'))

create_ancsitedat <- function(prev, sd){
  ancsitedat <- data.table(year = 2012:2022,
                            rnorm(n = length(2012:2022),
                                        mean = prev, sd = sd))
  colnames(ancsitedat) <- c('year', paste0('sd_', sd))
  return(ancsitedat)
}


var_.01 <- create_ancsitedat(prev = 0.26, sd = 0.01)
var_.05 <- create_ancsitedat(prev = 0.26, sd = 0.05)
var_.1 <- create_ancsitedat(prev = 0.26, sd = 0.1)
var_.15 <- create_ancsitedat(prev = 0.26, sd = 0.15)
ancsitedat <- merge(var_.01, var_.05, by = 'year')
ancsitedat <- merge(ancsitedat, var_.1, by = 'year')
ancsitedat <- merge(ancsitedat, var_.15, by = 'year')

dt <- merge(birth_data[year %in% ancsitedat$year], ancsitedat, by = 'year')

lm.loss <- function(par, dt, true_name){
  laf.par <- par[1]
  err.sigma <- par[2]
  true_prev <- dt[[true_name]]
  est.prev <- dt$prev_est
  ancrt_sd <- dt$sd

  likelihoods <- dnorm(true_prev, mean = est.prev * laf.par, sd = err.sigma)

  log.likelihoods <- log(likelihoods)

  deviance <- -2 * sum(log.likelihoods)

  return(deviance)
}
par_true <- c( laf =  1.29153691, err.sigma = 0.05316214)

get_val_sd <- function(starting = par_true, name = 'prev'){
  fit <- optim(par = starting, fn = lm.loss, hessian = T, dt = dt, true_name = name)
  hessian <- fit$hessian
  hessian.inv <- solve(hessian)
  parameter.se <- sqrt(diag(hessian.inv))

  out <- data.table(laf = fit$par[1],
                    lb = fit$par[1] - 1.96 * parameter.se[1],
                    ub = fit$par[1] + 1.96 * parameter.se[1],
                    sd = parameter.se[1],
                    run = name)
  return(out)
}


out <- lapply(c(setdiff(colnames(ancsitedat), 'year')), get_val_sd, starting = par_true)
out <- rbindlist(out)
