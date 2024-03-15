library(data.table)
library(dplyr)
library(ggplot2)
library(eppasm)
library(readxl)
library(leapfrog)

root <- getwd()
pjnz <- paste0(root, '/spectrum_data/Botswana2023v4 WPP 02_03_2023 KOS.PJNZ')
pop_1 <- paste0(root, '/spectrum_data/Botswana2023v4 WPP 02_03_2023 KOS_pop1.xlsx')
hivp <- leapfrog::prepare_leapfrog_projp(pjnz = pjnz)
hivp <- leapfrog:::prepare_hc_leapfrog_projp(pjnz, hivp)
demp <- leapfrog::prepare_leapfrog_demp(pjnz)
input <- qs::qread(paste0(root, '/spectrum_data/draws/birth_prev/1.qs'))
hivp$hiv_input <- input$hivstrat_adult
hivp$art_input <- input$artstrat_adult
hivp$totpop_input <- demp$basepop
hivp$laf <- rep(1, 61)
out <- leapfrogR(demp, projp = hivp)

calc_curr_last <- function(input, af = 1, t = 20){
  hivstrat_adult = input$hivstrat_adult
  artstrat_adult = input$artstrat_adult
  totpop <- input$totpop

  nHIVcurr = sum(c(hivstrat_adult[,af,2,t],
                   artstrat_adult[,,af,2,t]));
  nHIVlast =  sum(c(hivstrat_adult[,af,2,t-1],
                    artstrat_adult[,,af,2,t-1]));
  totpop = totpop[af, t];

  prev = nHIVcurr / totpop;

  return(list(nHIVcurr = nHIVcurr,
              nHIVlast = nHIVlast,
              totpop = totpop,
              prev = prev))
}

calc_df_hivpop <- function(input, af, t, laf = 1, hm = 1){
  hivstrat_adult = input$hivstrat_adult
  artstrat_adult = input$artstrat_adult
  fert_mult_by_age = input$fert_mult_by_age
  fert_mult_offart = input$fert_mult_offart

  df <- 0
  df <- laf * fert_mult_by_age[af] * fert_mult_offart[hm] * ((hivstrat_adult[hm, af, 2, t] + hivstrat_adult[hm, af, 2, t-1]) / 2)
  df <- df + laf * fert_mult_by_age[af] * fert_mult_offart[hm] * ((artstrat_adult[1, hm, af, 2, t] + artstrat_adult[1, hm, af, 2, t-1]) / 2)

  return(df)
}

calc_df_artpop <- function(input, af, t, laf, hm, hu){

  hivstrat_adult = input$hivstrat_adult
  artstrat_adult = input$artstrat_adult
  fert_mult_onart = input$fert_mult_onart

  df = laf * fert_mult_onart[af] * ((artstrat_adult[hu, hm, af, 2, t] + artstrat_adult[hu, hm, af, 2, t-1]) / 2)
  return(df)
}


calc_df_byaf <- function(input, af, t, laf = 1){
  hivstrat_adult = input$hivstrat_adult
  artstrat_adult = input$artstrat_adult
  fert_mult_by_age = input$fert_mult_by_age
  fert_mult_onart = input$fert_mult_onart
  fert_mult_offart = input$fert_mult_offart

  df <- sum(unlist(lapply(c(1:7), calc_df_hivpop, t=t, laf = laf, input = input, af = af))) +
    sum(unlist(lapply(c(1:7), FUN = calc_df_artpop, t=t, input = input, laf = laf, hu = 2, af = af))) +
    sum(unlist(lapply(c(1:7), FUN = calc_df_artpop, t=t, input = input, laf = laf, hu = 3, af = af)))

  curr_last <- calc_curr_last(input = input, af = af, t = t)
  nHIVcurr <- curr_last$nHIVcurr
  nHIVlast <- curr_last$nHIVlast

  if(nHIVcurr > 0){
    df = df / ((nHIVcurr + nHIVlast) / 2);
  }else{
    df = 1;
  }

 return(df)
}

calc_hiv_births_age_spec <- function(input, laf, af, t){
  tfr = input$tfr
  asfr = input$asfr
  asfr_sum = input$asfr_sum
  total_births = input$total_births
  totpop <- input$totpop

  curr_last <- calc_curr_last(input = input, af = af, t = t)
  nHIVcurr <- curr_last$nHIVcurr
  nHIVlast <- curr_last$nHIVlast
  totpop <- curr_last$totpop
  prev <- curr_last$prev

  df = calc_df_byaf(input = input, af = af, t = t, laf = laf)

  births <- (nHIVcurr + nHIVlast) / 2 * tfr[t] * df / (df * prev + 1 - prev) *  asfr[af, t] / asfr_sum[t] ;
  return(births)
}

calc_hiv_births_total_year <- function(input, laf, t){
  births <- Reduce('+', lapply(1:35, calc_hiv_births_age_spec, input = input, t = t, laf = laf))
  return(births)
}

calc_hiv_births <- function(draw, laf, year.idx = 1:61){
  input <- qs::qread(paste0(root, '/spectrum_data/draws/birth_prev/', draw, '.qs'))
  total_births = input$total_births[year.idx]
  birth_prev <- unlist(lapply(year.idx, calc_hiv_births_total_year, input = input, laf = laf)) #/ total_births
  # birth_prev[t]
  return(birth_prev)
}

calc_birth_prev <- function(laf = 1){
  birth_prev <- lapply(1:300, calc_hiv_births, laf = laf)
  return(birth_prev)
}

