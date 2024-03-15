library(data.table)
library(dplyr)
library(ggplot2)
library(eppasm)
library(readxl)
library(leapfrog)
root <- getwd()

cd4_progression_model <- function(dt, par){
  af <- dt[['af']]
  cd4 <- dt[['cd4']]

  b1 <- par[['prog1']] ##carrying capacity (0-1)
  b2 <- par[['prog2']] ##y intercept (0-1)
  b3 <- par[['prog3']] ##growth coefficient
  # b4 <- par[['b4']]
  # err <- par[['err']]
  cd4_vec = 1:7
  mean = mean(exp(b1) + (cd4_vec-1) * exp(b2))
  sd = sd(exp(b1) + (cd4_vec-1) * exp(b2))

  prog <- pnorm(exp(b1) + (cd4-1) * exp(b2), mean = mean, sd = sd)

  return(prog)
}

get_dt_prog <- function(dt,  par_vec){
  results <- apply(dt, 1, cd4_progression_model, par = par_vec)
  dt[,prog := results]
  dt <- dt[,.(cd4, af, prog)]
  return(dt)
}

fill_cd4_prog <- function(dt){
  cd4_prog_new <- array(0, dim = c(6,15))
  cd4_prog_new[1:6,1:5] <- dt[af == 1,prog]
  cd4_prog_new[1:5,6:15] <- dt[af == 2 & cd4 %in% 1:5,prog]

  return(cd4_prog_new)
}

new_cd4_prog <- function(dt, par_vec){
  dt <- get_dt_prog(dt, par_vec)
  cd4_prog <- fill_cd4_prog(dt)
  return(cd4_prog)
}
