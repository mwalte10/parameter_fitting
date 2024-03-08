rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(eppasm)
library(readxl)
library(leapfrog)
root <- getwd()

prep_spectrum_data <- function(pjnz = paste0(root, '/spectrum_data/Benin_2023_02_21final2_KOS.PJNZ')){
  pjnz <- paste0(root, '/spectrum_data/Benin_2023_02_21final2_KOS.PJNZ')
  hivp <- leapfrog::prepare_leapfrog_projp(pjnz = pjnz)
  hivp <- leapfrog:::prepare_hc_leapfrog_projp(pjnz, hivp)

  cd4_transition <- hivp$paed_cd4_transition
  dist <- hivp$paed_cd4_dist
  cd4_mort <- array(data = NA, dim = c(7, 4, 15))
  cd4_mort[1:7,1:4,1:5] <- hivp$paed_cd4_mort
  cd4_mort[1:6,1:4,6:15] <- hivp$adol_cd4_mort
  cd4_mort[7,1:4,6:15] <- 0
  cd4_prog <- array(data = NA, dim = c(7,15))
  cd4_prog[1:7,1:5] <- hivp$paed_cd4_prog
  cd4_prog[1:6,6:15] <- hivp$adol_cd4_prog
  cd4_prog[7,6:15] <- 0

  return(list(cd4_mort = cd4_mort, dist = dist, cd4_prog = cd4_prog, cd4_transition = cd4_transition))
}

prep_output_structure <- function(ts_in = 30/365, years_out = 15, seed_inf = 1000){
  hivstrat_paeds_out <- array(0, dim = c(7, 4, 15, 2, (years_out  * (1/ts_in))))
  hivstrat_paeds_out[,1,1,1:2,1] <- seed_inf * dist
  hivstrat_paeds_out[,2:3,1,1:2,1] <- seed_inf * dist
  hivstrat_paeds_out[,4,2,1:2,1] <- seed_inf * dist

  return(hivstrat_paeds_out)
}

cd4_progression <- function(hivstrat_paeds, cd4_mort, cd4_prog, ts){
  deaths_paeds <- array(0, dim = c(7, 4, 15, 2))
  grad_paeds <- array(0, dim = c(7, 4, 15, 2))

  for(g in 1:2){
    for(hm in 1:7){
      for(af in 1:15){
        for(cat in 1:4){
          deaths_paeds[hm, cat, af, g] = deaths_paeds[hm, cat, af, g] + hivstrat_paeds[hm, cat, af, g] -
            hivstrat_paeds[hm, cat, af, g] * cd4_mort[hm, cat, af] * ts
        }
      }
    }
  }

  for(g in 1:2){
    for(hm in 2:7){
      for(af in 1:15){
        for(cat in 1:4){
          grad_paeds[hm - 1, cat, af, g] = - 0.5 * ts * (deaths_paeds[hm-1, cat, af, g] * cd4_prog[hm - 1,af]  +
             hivstrat_paeds[hm-1, cat, af, g] * cd4_prog[hm - 1,af])

          grad_paeds[hm, cat, af, g] =  0.5 * ts * (deaths_paeds[hm-1, cat, af, g] * cd4_prog[hm - 1,af]  +
                                                       hivstrat_paeds[hm-1, cat, af, g] * cd4_prog[hm - 1,af])
        }
      }
    }
  }
  return(grad_paeds)
}

cd4_mortality <- function(hivstrat_paeds, cd4_mort, ts){
  grad_paeds <- array(0, dim = c(7, 4, 15, 2))
  for(g in 1:2){
    for(hm in 1:7){
      for(af in 1:15){
        for(cat in 1:4){
          grad_paeds[hm, cat, af, g] = grad_paeds[hm, cat, af, g] -
            ts * hivstrat_paeds[hm, cat, af, g] * cd4_mort[hm, cat, af]
        }
      }
    }
  }
  return(grad_paeds)
}

new_hivpop <- function(hivstrat_paeds, cd4_mort, cd4_prog,  cd4_transition, ts){
 out <- array(0, dim = c(7, 4, 15, 2))

 out[, , 1, ] <- out[, , 1, ] + hivstrat_paeds[, , 1, ] * (1-ts)

  for(g in 1:2){
    for(hm in 1:7){
      for(af in 2:5){
        for(cat in 1:4){
          ##assuming no background mortality
            out[hm, cat, af, g] <- out[hm, cat, af, g] + hivstrat_paeds[hm, cat, af-1, g] * ts
            out[hm, cat, af, g] <- out[hm, cat, af, g] + hivstrat_paeds[hm, cat, af, g]* (1-ts)

        }
      }
    }
  }

  for(g in 1:2){
    for(hm in 1:7){
      for(hm_alt in 1:6){
        for(cat in 1:4){
          ##assuming no background mortality
          out[hm_alt, cat, 6, g] <- out[hm_alt, cat, 6, g] + hivstrat_paeds[hm, cat, 5, g] * cd4_transition[hm_alt, hm] * ts

        }
      }
    }
  }

 out[, , 6, ] <- out[, , 6, ] + hivstrat_paeds[, , 6, ] * (1-ts)


  for(g in 1:2){
    for(hm in 1:6){
      for(af in 7:15){
        for(cat in 1:4){
          ##assuming no background mortality
          out[hm, cat, af, g] <- out[hm, cat, af, g] + hivstrat_paeds[hm, cat, af-1, g] * ts
          out[hm, cat, af, g] <- out[hm, cat, af, g] + hivstrat_paeds[hm, cat, af, g] * (1-ts)
        }
      }
    }
  }
  #out <- hivstrat_paeds

  grad_cd4_prog <- cd4_progression(out, cd4_mort, cd4_prog, ts)
  grad_cd4_mort <- cd4_mortality(out, cd4_mort, ts)

  for(g in 1:2){
    for(hm in 1:7){
      for(af in 1:15){
        for(cat in 1:4){
          out[hm, cat, af, g] <- out[hm, cat, af, g] + grad_cd4_prog[hm, cat, af, g] + grad_cd4_mort[hm, cat, af, g]
        }
      }
    }
  }
  return(out)

}

get_pct_surviving <- function(ts_in, years_out, seed_inf = 1000, input, input_dist, input_transition){

  dist <- input_dist
  cd4_mort <- array(unlist(input)[grep('cd4_mort', names(unlist(input)))], dim = c(7,4,15))
  cd4_prog <- array(unlist(input)[grep('cd4_prog', names(unlist(input)))], dim = c(7,15))
  cd4_transition <- input_transition

  hivstrat_paeds_out <- prep_output_structure()

  end_sim <- (years_out  * (1/ts_in)) - 1
  for(ts in 1:end_sim){
    hivstrat_paeds_out[, , , , ts+1] <-  new_hivpop(hivstrat_paeds_out[, , , , ts], cd4_mort, cd4_prog,cd4_transition = cd4_transition, ts = ts_in)
  }

  out <- data.table(melt(hivstrat_paeds_out))
  out <- out[,.(value = sum(value)), by = c('Var5', 'Var2')]
  out[Var2 == 1 & value > 0, keep := 'Perinatal']
  out[Var2 == 2 & value > 0, keep := 'BF 0-6 mo.']
  out[Var2 == 3 & value > 0, keep := 'BF 7-12 mo.']
  out[Var2 == 4 & value > 0, keep := 'BF 12+ mo.']
  out <- out[!is.na(keep)]
  out$keep <- factor(out$keep, levels = c('Perinatal', 'BF 0-6 mo.', 'BF 7-12 mo.', 'BF 12+ mo.'))
  out <- out[,.(ts = Var5, infection_type = keep, pct_surviving = value / (seed_inf * 2))]
  return(out)
}

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

optimization_function <- function(par = c(err, input), dt = surv, input_dist, input_transition){
  err <- par[1]
  input <- par[2:length(par)]
  out <- get_pct_surviving(ts_in = 30/365, years_out = 15, seed_inf = 1000, input = input, input_dist, input_transition)

  likelihoods <- dnorm(dt[type == 'Perinatal',surv], mean = out[infection_type == 'Perinatal',pct_surviving], sd = err)

  log.likelihoods <- log(likelihoods)
  #
  deviance <- -2 * sum(log.likelihoods)

  return(deviance)

}

surv <- qs::qread(paste0(root, '/survival_curves.qs'))
input <- prep_spectrum_data()
par_good <- c(err = 0.1, unlist(input[c('cd4_prog', 'cd4_mort')]))

opt_diffstep=1e-3

###I cant get this to converge
start <- Sys.time()
parameter.fits <- optim(par = par_good, fn = optimization_function, hessian = F, dt = surv,
                        #trace is the level of reporting
                        control=list(trace=4,
                                     #negative value turns it to a maximation problem
                                     fnscale = -1,
                        #maxit is the maximum number of iterations, have this really low rn to finish faster for testing
                                    # maxit = 1,
                        #frequency of reports
                        REPORT = 1), method="BFGS", input_dist = input$dist,
                        input_transition = input$cd4_transition)
end <- Sys.time()


##Fitting to modelled output and not data, so there isn't a model fit in the statistical sense
##this fit is dependent on what point on these are choosing
##lots of combos that could give this survival curve
##start to fit just a couple parameters rather than all of them
##constrain the parameters a bit more (rate parameters between zero and one)
##statistical inference: uncertainty in the parameter conditional on the data

##ideally would be fit to the same data used to make these survival curves
  ## if we had the uncertainty in the parameters that are fit to this data we could try to infer
  ## Could fit to the hazard function rather than the survival function -> this would reduce the dependence on earlier

