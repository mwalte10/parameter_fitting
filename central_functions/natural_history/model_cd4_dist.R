library(data.table)
library(dplyr)
library(ggplot2)
library(eppasm)
library(readxl)
library(leapfrog)
root <- getwd()

cd4_dist_model <- function(cd4, par){
  b1 <- par[['dist1']] ##carrying capacity (0-1)
  # b2 <- par[['dist2']] ##y intercept (0-1)
  # b3 <- par[['dist3']] ##growth coefficient
  # b4 <- par[['b4']]
  # err <- par[['err']]


  dist <- exp(-cd4 * b1) / sum(exp(-cd4 * b1))

  return(dist)
}


new_cd4_dist <- function(par_vec, struct_pars){
  cd4_dist <- cd4_dist_model(1:struct_pars$hDS, par_vec)
  return(cd4_dist)
}
