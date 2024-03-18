library(data.table)
library(dplyr)
library(ggplot2)
library(eppasm)
library(readxl)
library(leapfrog)
root <- getwd()

cd4_progression_model <- function(dt, par, struct_pars){
  af <- dt[['af']]
  cd4 <- dt[['cd4']]

  b1 <- par[['prog1']] ##carrying capacity (0-1)
  b2 <- par[['prog2']] ##y intercept (0-1)
  b3 <- par[['prog3']] ##growth coefficient
  # b4 <- par[['b4']]
  # err <- par[['err']]
  cd4_vec = 1:struct_pars$hDS
  mean = mean(exp(b1) + (cd4_vec-1) * exp(b2))
  sd = sd(exp(b1) + (cd4_vec-1) * exp(b2))

  prog <- pnorm(exp(b1) + (cd4-1) * exp(b2), mean = mean, sd = sd)

  return(prog)
}

get_dt_prog <- function(dt,  par_vec, struct_pars){
  results <- apply(dt, 1, cd4_progression_model, par = par_vec, struct_pars = struct_pars)
  dt[,prog := results]
  dt <- dt[,.(cd4, af, prog)]
  return(dt)
}

fill_cd4_prog <- function(dt, struct_pars){
  cd4_prog_new <- array(0, dim = c((struct_pars$hDS - 1),15))
  cd4_prog_new[,1:5] <- dt[af == 1,prog]
  cd4_prog_new[,6:15] <- dt[af == 2,prog]

  return(cd4_prog_new)
}

new_cd4_prog <- function(dt, par_vec, struct_pars){
  dt <- get_dt_prog(dt, par_vec, struct_pars)
  cd4_prog <- fill_cd4_prog(dt, struct_pars)
  return(cd4_prog)
}
