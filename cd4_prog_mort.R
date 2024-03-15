rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(eppasm)
library(readxl)
library(leapfrog)
root <- getwd()

lapply(list.files(paste0(root, '/central_functions/natural_history/'),full.names = T), source)

# ts_in <- 30/365
# years_out <- 15
# seed_inf = 1000
# input <- prep_spectrum_data()
#
# out <- get_pct_surviving(ts_in = ts_in, years_out = years_out, seed_inf = seed_inf, input = unlist(input[c('cd4_prog', 'cd4_mort')]),
#                          input_dist = input$dist, input_transition = input$cd4_transition)
#
# #dev.new(width = 7.5, height = 5,noRStudioGD = T)
# dir.create(paste0(root, '/plots'))
# png(filename = 'plots/surv plot.png', width = 7.5, height = 5, units = 'in', res = 300)
# ggplot(out[ts < 25], aes(ts , pct_surviving, col = as.factor(infection_type))) + geom_line(lwd = 2) +
#   scale_x_continuous(breaks = c(0,1,3,6,12,24), labels = c(0, 1, 3, 6, 12, 24)) +
#   labs(x = 'Months since infection', y = 'Percent surviving', col = 'Infection type') +
#   theme_bw(base_size = 12) +
#   theme(
#     strip.text.x = element_text(size = rel(1)), axis.ticks.y= element_blank())
# dev.off()

optimization_function <- function(par = par_vec, dt = surv, input, years_out, flags = c(cd4_mort =  T, cd4_prog = T)){
  err <- par[length(par)]
  pars <- par#[c(1:3)]
  print(pars)
  cd4_mort_flag = flags['cd4_mort']
  cd4_prog_flag = flags['cd4_prog']

  out <- get_pct_surviving(ts_in = 30/365, years_out = years_out, seed_inf = 1000, input = input, par_vec = pars,
                           cd4_mort_flag, cd4_prog_flag)

  dt <- dt[month < (years_out * 12 + 1),]
  dt[,hr := -log10(surv) / (month+1)]
  out[,month := 1:nrow(out)]
  out[,hr := -log10(pct_surviving) / (month)]

  # likelihoods.peri <- dnorm(dt[type == 'Perinatal',surv], mean = out[infection_type == 'Perinatal',pct_surviving], sd = exp(err))
  likelihoods.peri <- dnorm(dt[type == 'Perinatal',hr], mean = out[infection_type == 'Perinatal',hr], sd = exp(err))
  # likelihoods.bf <- dnorm(dt[type == 'BF',surv], mean = out[infection_type == 'BF',pct_surviving], sd = exp(err))

  log.likelihoods.peri <- log(likelihoods.peri)
  #log.likelihoods.bf <- log(likelihoods.bf)

 # deviance <- -2 * (sum(log.likelihoods.peri) +  sum(log.likelihoods.bf))
  deviance <- -2 * (sum(log.likelihoods.peri))


  return(deviance)

}

surv <- qs::qread(paste0(root, '/survival_curves.qs'))
surv <- surv[month != max(month),]
#prep_spectrum_data() ## located in the 00_prep_data file
input <- qs::qread(paste0(root, '/input.qs'))
pars <- c('mort1', 'mort2', 'mort3', 'mort4',
          'prog1', 'prog2', 'prog3',
          'b4', 'err','err_total')
par_vec <- rep(0, length(pars))
names(par_vec) <- pars
par_vec['mort1'] = 1
par_vec['mort2'] = -0.75
par_vec['mort3'] = 0.3
par_vec['mort4'] = 4.5
flags = c(cd4_mort =  T, cd4_prog = T)

years_out = 2
parameter.fits <- optim(par = par_vec, fn = optimization_function, hessian = F, dt = surv, input = input, years_out = years_out, flags = flags,
                        #trace is the level of reporting
                        control=list(trace=4,
                                     #negative value turns it to a maximation problem
                                     fnscale = 1,
                        #maxit is the maximum number of iterations, have this really low rn to finish faster for testing
                                    # maxit = 1,
                        #frequency of reports
                        REPORT = 1), method="Nelder-Mead")


##Fitting to modelled output and not data, so there isn't a model fit in the statistical sense
##this fit is dependent on what point on these are choosing
##lots of combos that could give this survival curve
##start to fit just a couple parameters rather than all of them
##constrain the parameters a bit more (rate parameters between zero and one)
##statistical inference: uncertainty in the parameter conditional on the data

##ideally would be fit to the same data used to make these survival curves
  ## if we had the uncertainty in the parameters that are fit to this data we could try to infer
  ## Could fit to the hazard function rather than the survival function -> this would reduce the dependence on earlier timesteps
surv_calc <- get_pct_surviving(ts_in = 30/365, years_out = years_out, seed_inf = 1000, par_vec = parameter.fits$par,
                  input = input, cd4_mort_flag = flags[1], cd4_prog_flag = flags[2])
surv_calc[,ts := 1:nrow(surv_calc) - 1]


surv_calc_og <- get_pct_surviving(ts_in = 30/365, years_out = years_out, seed_inf = 1000, par_vec = parameter.fits$par,
                               input = input, cd4_mort_flag = F, cd4_prog_flag = F)
surv_calc_og[,ts := 1:nrow(surv_calc) - 1]
setnames(surv_calc_og, 'pct_surviving', 'Spectrum- current parameters')
surv_calc <- merge(surv_calc, surv_calc_og, by = c('infection_type', 'ts', 'age'))
surv <- merge(surv, surv_calc, by.x = c('type', 'month'), by.y = c('infection_type', 'ts'))
surv[,age := NULL]
setnames(surv, c('pct_surviving','surv'), c('Spectrum- fit parameters', 'Survival curves'))
surv <- melt(surv, id.vars = c('type', 'month'))
surv$type <- factor(surv$type, levels = c('Perinatal', 'BF'))
ggplot(surv, aes(month, value, col = as.factor(variable))) + facet_wrap(~type) +
  geom_line() + ylim(0,1)

dt <- data.table(expand.grid(list(af = 1:15,
                                  cd4 = 1:7,
                                  tt = 1:4)))

cd4_mort = new_cd4_mort(dt, parameter.fits$par)
cd4_mort[,1,1:5]

dt <- data.table(expand.grid(list(af = 1:2,
                                  cd4 = 1:6)))
cd4_prog = new_cd4_prog(dt, parameter.fits$par)
cd4_prog[,1:5]
