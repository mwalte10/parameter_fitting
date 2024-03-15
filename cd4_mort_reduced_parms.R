
cd4_mortality <- function(dt, par){
  ag <- dt[['ag']]
  cd4 <- dt[['cd4']]
  tt <- dt[['tt']]

  b1 <- par[['b1']]
  b2 <- par[['b2']]
  b3 <- par[['b3']]
  b4 <- par[['b4']]
  err <- par[['err']]

  ##need to constrain mort betwen zero and one
  mort <- pnorm(1/ (1 + exp((b1 * cd4 + b2 * cd4 * ag))) + b3 * ag + b4 * tt + err)

  return(mort)
}

get_dt <- function(dt,  par_vec){
  results <- apply(dt, 1, cd4_mortality, par = par_vec)
  dt[,mort := results]
  dt[,mort := ifelse(cd4 ==7 & ag == 3, 0, mort)]
  dt <- dt[,.(cd4, tt, ag, mort)]
  return(dt)
}

fill_cd4_mort <- function(dt){
  cd4_mort_new <- array(0, dim = c(7,4,15))
  cd4_mort_new[,1,1:3] <- dt[ag == 1 & tt == 1,mort]
  cd4_mort_new[,2,1:3] <- dt[ag == 1 & tt == 2,mort]
  cd4_mort_new[,3,1:3] <- dt[ag == 1 & tt == 3,mort]
  cd4_mort_new[,4,1:3] <- dt[ag == 1 & tt == 4,mort]

  cd4_mort_new[,1,4:5] <- dt[ag == 2 & tt == 1,mort]
  cd4_mort_new[,2,4:5] <- dt[ag == 2 & tt == 2,mort]
  cd4_mort_new[,3,4:5] <- dt[ag == 2 & tt == 3,mort]
  cd4_mort_new[,4,4:5] <- dt[ag == 2 & tt == 4,mort]

  cd4_mort_new[,1,6:15] <- dt[ag == 3 & tt == 1,mort]
  cd4_mort_new[,2,6:15] <- dt[ag == 3 & tt == 2,mort]
  cd4_mort_new[,3,6:15] <- dt[ag == 3 & tt == 3,mort]
  cd4_mort_new[,4,6:15] <- dt[ag == 3 & tt == 4,mort]
  return(cd4_mort_new)
}

new_cd4_mort <- function(dt, par_vec){
  dt <- get_dt(dt, par_vec)
  cd4_mort <- fill_cd4_mort(dt)
  return(cd4_mort)
}

par_vec <- rep(0.01, 5)
names(par_vec) <- c('b1', 'b2', 'b3', 'b4', 'err')
dt <- data.table(expand.grid(list(ag = 1:3,
                                  cd4 = 1:7,
                                  tt = 1:4)))

new_cd4_mort(dt, par_vec)
