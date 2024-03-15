library(data.table)
library(dplyr)
library(ggplot2)
library(eppasm)
library(readxl)
library(leapfrog)
root <- getwd()

cd4_mortality_model <- function(dt, par){
  af <- dt[['af']]
  cd4 <- dt[['cd4']]
  tt <- dt[['tt']]

  b1 <- par[['mort1']] ##carrying capacity (0-1)
  b2 <- par[['mort2']] ##y intercept (0-1)
  b3 <- par[['mort3']] ##growth coefficient
  b4 <- par[['mort4']]
  # err <- par[['err']]
  cd4_vec = 1:7
  mean = b3 * (b4-1) / (1 + exp((-b1*(b4-1))))  + b2 * af
  sd = sd(b3 * (cd4_vec-1) / (1 + exp((-b1*(cd4_vec-1))))  + b2 * af)
  mean = 0
  sd = 1

  #mort = pnorm(b3 * (cd4-1) / (1 + exp((-b1*(cd4-1)))) + b2, mean = mean, sd = sd)
  mort = pnorm(b3 * (cd4-1) / (1 + exp((-b1*(cd4-1)))) + b2 * af, mean = mean, sd = sd)



  return(mort)
}

get_dt_mort <- function(dt,  par_vec){
  results <- apply(dt, 1, cd4_mortality_model, par = par_vec)
  dt[,mort := results]
  dt[,mort := ifelse(cd4 ==7 & af > 5, 0, mort)]
  dt <- dt[,.(cd4, tt, af, mort)]
  return(dt)
}

fill_cd4_mort <- function(dt){
  cd4_mort_new <- array(0, dim = c(7,4,15))
  for(af.x in 1:15){
    for(tt.x in 1:4){
      cd4_mort_new[,tt.x,af.x] <- dt[af == af.x & tt == tt.x,mort]
    }
  }
  # cd4_mort_new[,1,1:3] <- dt[af == 1 & tt == 1,mort]
  # cd4_mort_new[,2,1:3] <- dt[af == 1 & tt == 2,mort]
  # cd4_mort_new[,3,1:3] <- dt[af == 1 & tt == 3,mort]
  # cd4_mort_new[,4,1:3] <- dt[af == 1 & tt == 4,mort]
  #
  # cd4_mort_new[,1,4:5] <- dt[af == 2 & tt == 1,mort]
  # cd4_mort_new[,2,4:5] <- dt[af == 2 & tt == 2,mort]
  # cd4_mort_new[,3,4:5] <- dt[af == 2 & tt == 3,mort]
  # cd4_mort_new[,4,4:5] <- dt[af == 2 & tt == 4,mort]
  #
  # cd4_mort_new[,1,6:15] <- dt[af == 3 & tt == 1,mort]
  # cd4_mort_new[,2,6:15] <- dt[af == 3 & tt == 2,mort]
  # cd4_mort_new[,3,6:15] <- dt[af == 3 & tt == 3,mort]
  # cd4_mort_new[,4,6:15] <- dt[af == 3 & tt == 4,mort]
  return(cd4_mort_new)
}

new_cd4_mort <- function(dt, par_vec){
  dt <- get_dt_mort(dt, par_vec)
  cd4_mort <- fill_cd4_mort(dt)
  return(cd4_mort)
}
