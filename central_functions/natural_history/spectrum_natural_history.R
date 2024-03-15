#rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(eppasm)
library(readxl)
library(leapfrog)
root <- getwd()

prep_output_structure <- function(ts_in = 30/365, years_out = 15, seed_inf = 1000, dist){
  hivstrat_paeds_out <- array(0, dim = c(7, 4, (years_out  * (1/ts_in))))
  hivstrat_paeds_out[,1,1] <- seed_inf * dist
  # hivstrat_paeds_out[,2,3] <- seed_inf * dist
  # hivstrat_paeds_out[,3,8] <- seed_inf * dist
  # hivstrat_paeds_out[,4,((1/ts_in) + 1)] <- seed_inf * dist

  return(hivstrat_paeds_out)
}

cd4_progression <- function(hivstrat_paeds, cd4_mort, cd4_prog, ts_step, ts_current){
  deaths_paeds <- array(0, dim = c(7, 4))
  grad_paeds <- array(0, dim = c(7, 4))
  age = floor(ts_current*ts_step) + 1

    for(hm in 1:7){
        for(cat in 1:3){
          deaths_paeds[hm, cat] = hivstrat_paeds[hm, cat] - hivstrat_paeds[hm, cat] * cd4_mort[hm, cat, age] * ts_step
        }
      if((age+1) < 15){
        deaths_paeds[hm, 4] =  hivstrat_paeds[hm, 4] - hivstrat_paeds[hm, 4] * cd4_mort[hm, 4, (age + 1)] * ts_step
      }
    }

    for(hm in 2:7){
        for(cat in 1:3){
          grad_paeds[hm - 1, cat] = - 0.5 * ts_step * (deaths_paeds[hm-1, cat] * cd4_prog[hm - 1,age]  +
                                                           hivstrat_paeds[hm-1, cat] * cd4_prog[hm - 1,age])

          grad_paeds[hm, cat] =  0.5 * ts_step * (deaths_paeds[hm-1, cat] * cd4_prog[hm - 1,age]  +
                                                      hivstrat_paeds[hm-1, cat] * cd4_prog[hm - 1,age])
        }
      if((age+1) < 15){
      grad_paeds[hm - 1, 4] = - 0.5 * ts_step * (deaths_paeds[hm-1, 4] * cd4_prog[hm - 1,age+1]  +
                                                hivstrat_paeds[hm-1, 4] * cd4_prog[hm - 1,age+1])

      grad_paeds[hm, 4] =  0.5 * ts_step * (deaths_paeds[hm-1, 4] * cd4_prog[hm - 1,age+1]  +
                                           hivstrat_paeds[hm-1, 4] * cd4_prog[hm - 1,age+1])}
    }
  return(grad_paeds)
}

cd4_mortality <- function(hivstrat_paeds, cd4_mort, ts_step, ts_current){
  grad_paeds <- array(0, dim = c(7, 4))
  age = floor(ts_current*ts_step) + 1
    for(hm in 1:7){
        for(cat in 1:3){
          grad_paeds[hm, cat] = -ts_step * hivstrat_paeds[hm, cat] * cd4_mort[hm, cat, age]
        }
      if((age + 1) < 15){
        grad_paeds[hm, 4] = -ts_step * hivstrat_paeds[hm, 4] * cd4_mort[hm, 4, (age + 1)]
      }
    }
  return(grad_paeds)
}

new_hivpop <- function(hivstrat_paeds, cd4_mort, cd4_prog, cd4_transition, ts_step, ts_current){
  grad_cd4_prog <- cd4_progression(hivstrat_paeds, cd4_mort, cd4_prog, ts_step, ts_current )
  grad_cd4_mort <- cd4_mortality(hivstrat_paeds, cd4_mort, ts_step, ts_current )

    for(hm in 1:7){
        for(cat in 1:4){
          hivstrat_paeds[hm, cat] <- hivstrat_paeds[hm, cat] + grad_cd4_prog[hm, cat] + grad_cd4_mort[hm, cat]
        }
    }

  return(hivstrat_paeds)

}

get_pct_surviving <- function(ts_in, years_out, seed_inf = 1000, input, par_vec, cd4_mort_flag = T, cd4_prog_flag = T){
  if(is.null(par_vec)){
    cd4_mort = input$cd4_mort
  }else{
    if(cd4_mort_flag){
      dt <- data.table(expand.grid(list(af = 1:15,
                                        cd4 = 1:7,
                                        tt = 1:4)))

      cd4_mort = new_cd4_mort(dt, par_vec)
    }else{
      cd4_mort = input$cd4_mort
    }
  }

  if(is.null(par_vec)){
    cd4_prog = input$cd4_prog
  }else{
    if(cd4_prog_flag){
      dt <- data.table(expand.grid(list(af = 1:2,
                                        cd4 = 1:6)))

      cd4_prog = new_cd4_prog(dt, par_vec)
    }else{
      cd4_prog = input$cd4_prog
    }
  }

  dist <- input$dist
  # cd4_prog[1:6,1:5] <- 0.1
  # cd4_prog[1:5,1:5] <- 0.1

  cd4_transition <- input$cd4_transition

  hivstrat_paeds_out <- prep_output_structure(ts_in, years_out, seed_inf, dist = dist)

  end_sim <- (years_out  * (1/ts_in)) - 1
  for(ts in 1:end_sim){
    hivstrat_paeds_out[, , ts+1] <-  new_hivpop(hivstrat_paeds_out[, , ts],
                                                    cd4_mort = cd4_mort,
                                                    cd4_prog = cd4_prog ,
                                                    cd4_transition = cd4_transition,
                                                    ts_step = ts_in,
                                                    ts_current = ts)
  }

  out <- data.table(melt(hivstrat_paeds_out))
  out[,age := Var3 * ts_in]
  out <- out[,.(value = sum(value)), by = c('age', 'Var2')]
  out <- out[Var2 == 1, keep := 'Perinatal']
  out <- out[keep == 'Perinatal']
  out <- out[,.(age, infection_type = keep, pct_surviving = value / 1000)]

  return(out)
}

# input <- qs::qread(paste0(root, '/input.qs'))
# years_out = 15
# ts_in = 0.1
# seed_inf = 1000
# surv_calc <- get_pct_surviving(ts_in = 0.1, years_out = years_out, seed_inf = 1000,
#                                input = input)
#
# ggplot(surv_calc, aes(age, pct_surviving, col = as.factor(infection_type))) + geom_line()
