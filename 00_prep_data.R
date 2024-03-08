library(data.table)
library(dplyr)
library(ggplot2)
library(eppasm)
library(readxl)

root <- getwd()
pjnz <- paste0(root, '/spectrum_data/Botswana2023v4 WPP 02_03_2023 KOS.PJNZ')
pop_1 <- paste0(root, '/spectrum_data/Botswana2023v4 WPP 02_03_2023 KOS_pop1.xlsx')
unc <- paste0(root, '/spectrum_data/uncertainty_file_300draws.xlsx')
hiv_births <- readxl::read_xlsx(unc, sheet = "101. Mothers needing PMTCT") %>% data.table()
births <- readxl::read_xlsx(unc, sheet = "21. Total Births" ) %>% data.table()

reshape_draws <- function(x){
  colnames(x) <- as.character(x[2,])
  x <- x[-(1:2),]
  x <- melt(x, id.vars = c('Iteration', 'Sex'))
  return(x)
}

hiv_births <- reshape_draws(hiv_births)
births <- reshape_draws(births)

birth_prev <- merge(hiv_births[,.(Iteration, variable, value)],
                    births[Sex == 'Male+Female',.(Iteration, variable, value)], by = c('Iteration', 'variable'))
birth_prev <- birth_prev[,.(draw = Iteration, year = as.integer(as.character(variable)), prev_est = value.x / value.y)]
write.csv(birth_prev, paste0(root, '/spectrum_data/birth_draws_laf1.csv'),row.names = F)


fitdata <- eppasm::prepare_spec_fit(pjnz, proj.end = 2023)
ancsitedat <- rbind(data.table(attr(fitdata$Urban, "eppd")$ancsitedat)[,urban := T],
                    data.table(attr(fitdata$Rural, "eppd")$ancsitedat)[,urban := F])
ancsitedat <- ancsitedat[type == 'ancrt',]
write.csv(ancsitedat, paste0(root, '/spectrum_data/ancrt_data.csv'),row.names = F)

unc <- paste0(root, '/spectrum_data/uncertainty_file_300draws_calibratedlaf.xlsx')
hiv_births <- readxl::read_xlsx(unc, sheet = "101. Mothers needing PMTCT") %>% data.table()
births <- readxl::read_xlsx(unc, sheet = "21. Total Births" ) %>% data.table()


hiv_births <- reshape_draws(hiv_births)
births <- reshape_draws(births)

birth_prev <- merge(hiv_births[,.(Iteration, variable, value)],
                    births[Sex == 'Male+Female',.(Iteration, variable, value)], by = c('Iteration', 'variable'))
birth_prev <- birth_prev[,.(draw = Iteration, year = as.integer(as.character(variable)), prev_est = value.x / value.y)]
write.csv(birth_prev, paste0(root, '/spectrum_data/birth_draws_laf_calibrated.csv'),row.names = F)


inc <- readxl::read_xlsx(unc, sheet = "53. Incidence (15-49) (Percent)") %>% data.table()
inc <- reshape_draws(inc)
inc <- inc[Sex == 'Male+Female']

dir.create(paste0(root, '/spectrum_data/draws/'))
fp <- prepare_directincid(pjnz)
for(i in 1:300){
  fp$incidinput <- inc[Iteration == i, value] / 100
  out <- eppasm::simmod(fp)
  save <- list(hivpop = attr(out, 'hivpop'), artpop = attr(out, 'artpop'))
  qs::qsave(save, paste0(root, '/spectrum_data/draws/', i, '.qs'))
  print(i / 300)
}


prep_birth_prev_struc <- function(draw, pjnz){

  out <- qs::qread(paste0(root, '/spectrum_data/draws/', draw, '.qs'))
  hivstrat_adult <- out$hivpop
  artstrat_adult <- out$artpop
  hivp = leapfrog::prepare_leapfrog_projp(pjnz)
  fert_mult_by_age <- hivp$fert_mult_by_age
  fert_mult_offart <- hivp$fert_mult_offart
  fert_mult_onart <- hivp$fert_mult_onart
  demp <- leapfrog::prepare_leapfrog_demp(pjnz)
  tfr <- demp$tfr
  asfr <- unique(demp$asfr)
  totpop <- demp$basepop
  totpop <- data.table(melt(totpop))
  totpop <- totpop[,.(Var1 = floor(Var1/5)  * 5, value), by = c('Var2', 'Var3')]
  totpop <- totpop[,.(value  = sum(value)), by = c('Var1', 'Var2', 'Var3')]
  totpop <- unique(totpop)
  totpop <- dcast(totpop[Var2 == 'Female' & Var1 %in% 15:45,], Var1 ~ Var3, value.var = 'value')
  totpop <- array(unlist(totpop[,2:62]), dim = c(7,61), dimnames = list(age = seq(15,45, by = 5), year = 1970:2030))

  out <- list(hivstrat_adult = hivstrat_adult,
              artstrat_adult = artstrat_adult,
              fert_mult_by_age = fert_mult_by_age[unique(floor(15:49/5)  * 5) - 14],
              fert_mult_onart = fert_mult_onart[unique(floor(15:49/5)  * 5) - 14],
              fert_mult_offart = fert_mult_offart,
              tfr = tfr,
              asfr = asfr,
              asfr_sum = colSums(asfr),
              total_births = demp$births,
              totpop = totpop)

  qs::qsave(out, paste0(root, '/spectrum_data/draws/birth_prev/', draw, '.qs'))
}

dir.create(paste0(root, '/spectrum_data/draws/birth_prev/'))
lapply(1:300, prep_birth_prev_struc, pjnz = pjnz)
