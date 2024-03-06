



##See which countries have ancrt data
pjnz_vec <- list.files('C:/Users/mwalters/Imperial College London/HIV Inference Group - WP - 2023 final restricted/Public spectrum files', full.names = T)
pjnz_vec <- c(pjnz_vec[grepl('.pjnz', pjnz_vec)],
              pjnz_vec[grepl('.PJNZ', pjnz_vec)])



dp_read_incidence_option <- function(dp) {

  dp <- leapfrog:::get_dp_data(dp)
  if (leapfrog:::exists_dptag(dp, "<IncidenceOptions MV>") ) {
    incidence_model_int <- as.integer(leapfrog:::dpsub(dp, "<IncidenceOptions MV>", 2, 4))
    if(incidence_model_int == 1 & leapfrog:::dpsub(dp, "<EpidemicTypeFromEPP MV>", 2, 4) == 'GENERALIZED'){
      incidence_model_int = 6
    }
  } else {
    stop("<IncidenceOptions MV> tag not found. Function probably needs update for this .DP file.")
  }

  aim_incidence_options <- c("Direct incidence input" = 0,
                             "EPP, concentrated" = 1,
                             "AEM" = 2,
                             "CSAVR" = 3,
                             "Fit to mortality data" = 4,
                             "ECDC" = 5,
                             'EPP, generalized' = 6)

  names(aim_incidence_options)[match(incidence_model_int, aim_incidence_options)]
}
get_anc_data <- function(pjnz){
  print(pjnz)
  type = dp_read_incidence_option(pjnz)
  if(type %in% c('EPP, generalized')){
    fitdata <- eppasm::prepare_spec_fit(pjnz, proj.end = 2023)
    ancsitedat <- rbind(data.table(attr(fitdata$Urban, "eppd")$ancsitedat)[,urban := T],
                        data.table(attr(fitdata$Rural, "eppd")$ancsitedat)[,urban := F], fill = T)
    ancsitedat[,loc := pjnz]
    return(ancsitedat)
  }

}

pjnz_vec <- setdiff(pjnz_vec,
                    c("C:/Users/mwalters/Imperial College London/HIV Inference Group - WP - 2023 final restricted/Public spectrum files/Haiti_2023_V9_revsh90.pjnz",
                      "C:/Users/mwalters/Imperial College London/HIV Inference Group - WP - 2023 final restricted/Public spectrum files/namibia_15052023.pjnz",
                      "C:/Users/mwalters/Imperial College London/HIV Inference Group - WP - 2023 final restricted/Public spectrum files/Tanzania_National_Zonal_19Apr2023 KOS.pjnz",
                      "C:/Users/mwalters/Imperial College London/HIV Inference Group - WP - 2023 final restricted/Public spectrum files/Zambia_2023_March161 KOS.pjnz",
                      "C:/Users/mwalters/Imperial College London/HIV Inference Group - WP - 2023 final restricted/Public spectrum files/Mozambique_6.29_29.05.2023.PJNZ",
                      "C:/Users/mwalters/Imperial College London/HIV Inference Group - WP - 2023 final restricted/Public spectrum files/PNG_2023_CurrentModel_v11.1.PJNZ",
                      "C:/Users/mwalters/Imperial College London/HIV Inference Group - WP - 2023 final restricted/Public spectrum files/South Sudan 2023 final.PJNZ"))
ancdata <- rbindlist(lapply(pjnz_vec, get_anc_data), fill = T)
ancdata <- ancdata[!is.na(n)]
ancdata <- ancdata[type == 'ancrt']

ancdata <- ancdata[,.(prev = weighted.mean(prev, n), n), by = c('year','loc')]
ancdata <- unique(ancdata[,.(prev, n = sum(n)), by = c('year', 'loc')])
dt <- dcast(ancdata[,.(loc, year)], loc ~ year, fun.aggregate = 'length')

ggplot() + geom_point(ancdata, aes(year, ))
