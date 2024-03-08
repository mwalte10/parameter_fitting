rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(eppasm)
library(readxl)

root <- getwd()
pjnz <- paste0(root, '/spectrum_data/Botswana2023v4 WPP 02_03_2023 KOS.PJNZ')
pop_1 <- paste0(root, '/spectrum_data/Botswana2023v4 WPP 02_03_2023 KOS_pop1.xlsx')
birth_data <- fread(paste0(root, '/spectrum_data/birth_draws_laf1.csv'))
ancsitedat <- fread(paste0(root, '/spectrum_data/ancrt_data.csv'))

ancsitedat <- ancsitedat[type == 'ancrt']
ancsitedat <- ancsitedat[,.(prev, w.prev = weighted.mean(prev, n), n), by = c('year')]
ancsitedat[,diff := (prev - w.prev)^2]
ancsitedat[,diff := sum(diff), by = 'year']
ancsitedat <- ancsitedat[,.(w.prev, n = sum(n)), by = c('year')]
ancsitedat[,var := w.prev / (n-1)]
ancsitedat <- unique(ancsitedat)
dt <- merge(birth_data[year %in% ancsitedat$year], unique(ancsitedat), by = c('year'), allow.cartesian = T)

logitTransform <- function(p) { log(p/(1-p)) }
invlogitTransform <- function(p) { exp(p) / (1+exp(p)) }


# dt[,offset_ancrt := (prev*n+0.5)/(n+1)]
# dt[,W.ancrt := qnorm(offset_ancrt)]
# dt[,v.ancrt := 2*pi*exp(W.ancrt^2)*offset_ancrt*(1-offset_ancrt)/n]
# dt[,W.prev_est := qnorm(prev_est)]

lm.loss <- function(par, dt, draw.x){
  laf.par <- par[1]
  err.laf <- par[2]
  prev <- dt[draw == draw.x,w.prev]
  prev_est <- dt[draw == draw.x,prev_est]
  v.ancrt <- dt[draw == draw.x,var]

  likelihoods <- dnorm(log(prev), mean = log(prev_est) + laf.par, sd = sqrt(err.laf + v.ancrt), log = T)
  deviance <- -2 * sum(likelihoods)

  return(deviance)
}


##this needs to be run separately for each draw
fits <- list()
for(i in 1:300){
  parameter.fits <- optim(par = c(laf = 1, err.laf = 0.05), fn = lm.loss, hessian = F, dt = dt, draw.x = i)
  hessian <- numDeriv::hessian(lm.loss, parameter.fits$par, dt = dt[draw == i,], draw.x = i)
  hessian.inv <- solve(hessian)
  parameter.se <- sqrt(diag(hessian.inv))
  out <- data.table(ml = parameter.fits$par,
                    se = parameter.se,
                    var = c('laf', 'err'),
                    draw = i)
  fits[[i]] <- out
  print(i / 300)
}

fits <- rbindlist(fits)
fits <- melt(fits, id.vars = c('var', 'draw'))
fits <- dcast(fits, draw + variable ~ var, value.var = 'value')

birth_data <- merge(birth_data, fits[variable == 'ml'], by = 'draw', allow.cartesian = T)
birth_data[,prev_est_new := prev_est * exp(laf)]
birth_data[,median := median(prev_est_new), by = 'year']
birth_data[,lower := quantile(prev_est_new, 0.025), by = 'year']
birth_data[,upper := quantile(prev_est_new, 0.975), by = 'year']

fits_mean <- fits[,.(err = median(err)), by = c('variable')]

spec_birth_data <- fread(paste0(root, '/spectrum_data/birth_draws_laf_calibrated.csv'))
spec_birth_data <- unique(spec_birth_data[,.(prev = median(prev_est), lower = quantile(prev_est, 0.025), upper = quantile(prev_est, 0.975), run = 'Current Spectrum'), by  = 'year'])
birth_data <- unique(birth_data[,.(year, prev = median, lower, upper, run = 'LAF fit 1.28')])
birth_data <- rbind(birth_data, spec_birth_data)
mean <- unique(ancsitedat)
mean <- mean[,err := fits_mean[variable == 'ml',err]]
mean <- mean[,.(year, w.prev, lower = w.prev - 1.96 * sqrt(var + err), upper =  w.prev + 1.96 * sqrt(var + err))]

dev.new(width = 7.5, height = 5,noRStudioGD = T)
dir.create(paste0(root, '/plots'))
png(filename = 'plots/LAF plot.png', width = 7.5, height = 5, units = 'in', res = 300)
ggplot() + geom_line(data = birth_data[year > 1999 & year < 2026], aes(year, prev, col = as.factor(run))) +
  geom_ribbon(data = birth_data[year > 1999 & year < 2026], aes(x = year, ymin = lower, ymax = upper, fill = as.factor(run)), alpha = 0.2) +
  geom_point(data = mean, aes(year, w.prev),  shape = 17, show.legend = F) +
 # geom_errorbar(data = mean, aes(year, ymin = lower, ymax = upper), ) +
   labs(size = 'Women tested at ANC-RT site',
       fill = 'LAF', color = 'LAF', x = NULL, y = 'HIV prevalence among pregnant women',
       caption = 'Triangles are national HIV prevalence')  + theme_bw(base_size = 11) +
  theme(
         strip.text.x = element_text(size = rel(0.8)),
         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1),
         axis.text.y = element_blank(), axis.ticks.y= element_blank(),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
