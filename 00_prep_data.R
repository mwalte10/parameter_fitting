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
