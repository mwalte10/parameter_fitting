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


dt <- merge(birth_data[year %in% ancsitedat$year], ancsitedat, by = c('year'), allow.cartesian = T)
dt <- dt[type == 'ancrt']

##Need to figure out where the err.sigma is coming from
####That's the part that I don't uncertand in Jeff's code where its the variance of the probit transformed prevalence
lm.loss <- function(par, true_prev, est.prev){
  laf.par <- par[1]
  err.laf <- par[2]

  ##I think that this sd also need to include variance from true_prev and est.prev, just not sure how to combine all of those
  likelihoods <- dnorm(true_prev, mean = est.prev * laf.par, sd = err.laf)
#
  log.likelihoods <- log(likelihoods)
#
  deviance <- -2 * sum(log.likelihoods)

  return(deviance)
}

par_bad <- c( laf = 8, err.laf = 3)
par_good <- c( laf = 1, err.laf = 0.2)
par_true <- c( laf =  1.29153691, err.laf = 0.05316214)

lm.loss(par_bad, true_prev = dt$prev, est.prev = dt$prev_est)
lm.loss(par_good, true_prev = dt$prev, est.prev = dt$prev_est)
lm.loss(par_true, true_prev = dt$prev, est.prev = dt$prev_est)

parameter.fits <- optim(par = par_good, fn = lm.loss, hessian = T, true_prev = dt$prev, est.prev = dt$prev_est)
hessian <- parameter.fits$hessian
hessian.inv <- solve(hessian)
parameter.se <- sqrt(diag(hessian.inv))

parameter.se
out <- data.table(ml = parameter.fits$par,
                  lb = parameter.fits$par - 1.96 * parameter.se,
                  ub = parameter.fits$par + 1.96 * parameter.se,
                  se = parameter.se)


rownames(out) <- c('LAF', 'Standard deviation of errors')

##Is the posterior distribution just the rnorm with the new mean and sd?
posterior_laf <- data.table(laf = rnorm(300, out$ml[1], sd = out$se[1]), draw = 1:300)
birth_data <- merge(birth_data, posterior_laf, by = 'draw')
birth_data[,prev_median := median(prev_est), by = 'year']
birth_data[,prev_est_new := prev_median * laf]
birth_data[,upper_og := quantile(prev_est, 0.975), by = 'year']
# birth_data[,upper_new := quantile(prev_est_new, 0.975), by = 'year']
birth_data[,upper_new := prev_est_new + 1.96 * 7.433817e-04, by = 'year']
birth_data[,lower_og := quantile(prev_est, 0.025), by = 'year']
birth_data[,lower_new := prev_est_new -  1.96 * 7.433817e-04, by = 'year']
# birth_data[,lower_new := quantile(prev_est_new , 0.025), by = 'year']

birth_data <- unique(birth_data[,.(prev_est = median(prev_est), upper_og, upper_new, lower_og, lower_new, prev_est_new = median(prev_est_new)), by = 'year'])
birth_data <- rbind(birth_data[,.(year, prev = prev_est, upper = upper_og, lower = lower_og, run = 'LAF = 1')],
                    birth_data[,.(year, prev = prev_est_new, upper = upper_new, lower = lower_new, run = 'LAF = 1.29')])

spec_birth_data <- fread(paste0(root, '/spectrum_data/birth_draws_laf_calibrated.csv'))
spec_birth_data <- unique(spec_birth_data[,.(prev = median(prev_est), lower = quantile(prev_est, 0.025), upper = quantile(prev_est, 0.975), run = 'Current Spectrum'), by  = 'year'])
birth_data <- rbind(birth_data, spec_birth_data)

mean <- ancsitedat[type == 'ancrt']
mean <- mean[,.(prev = weighted.mean(prev,n), n), by = 'year']
mean <- mean[,.(prev, n = sum(n)), by = 'year']
mean <- unique(mean)

dev.new(width = 7.5, height = 5,noRStudioGD = T)
dir.create(paste0(root, '/plots'))
png(filename = 'plots/LAF plot.png', width = 7.5, height = 5, units = 'in', res = 300)
ggplot() + geom_line(data = birth_data[year > 1999 & year < 2026], aes(year, prev, col = as.factor(run))) +
  geom_ribbon(data = birth_data[year > 1999 & year < 2026], aes(x = year, ymin = lower, ymax = upper, fill = as.factor(run)), alpha = 0.2) +
  geom_point(data = ancsitedat[type == 'ancrt'], aes(year, prev, size = n), alpha = 0.5) +
  geom_point(data = mean, aes(year, prev, size = n), col = 'red', shape = 17, show.legend = F) +
  labs(size = 'Women tested at ANC-RT site',
       fill = 'LAF', color = 'LAF', x = NULL, y = 'HIV prevalence among pregnant women',
       caption = 'Red triangles are national HIV prevalence')  + theme_bw(base_size = 11) +
  theme(
         strip.text.x = element_text(size = rel(0.8)),
         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1),
         axis.text.y = element_blank(), axis.ticks.y= element_blank(),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
