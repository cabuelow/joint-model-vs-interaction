# simulate fish biomass data assuming standardised effects for 
# reserve age, size, and other environmental factors
# then fit models to recover parameters

library(tidyverse)
library(rstan)
library(INLA)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = T)

# simulate data ------------------------------
# number of observation in reserves vs. fished locations?
# loop through from 1% reserves to 41% by increments of 10

loop <- c(0.01, 0.11, 0.21, 0.31, 0.41) # loop increments
total_n <- 10000 # total n
out <- list() # list for storing results

for(i in 1:length(loop)){
  set.seed(1)
  n_res <- total_n*loop[i] # reserves
  n_fis <- total_n - n_res # fished
  
  # covariates
  x_mpa <- rep(1, n_res)
  x_restrictions <- rbinom(n_fis, 1, 0.5) # 50% of unfished sites are restricted
  x_age <- round(runif(n_res, min=1, max=80))
  x_sage <- x_age/sd(x_age) # scale to have an sd of 1
  x_size <- rnorm(n_res, mean = 0, sd = 1) 
  x_enviro_r <- rnorm(n_res, mean = 0, sd = 1) # restricted environmental covariate
  x_enviro_f <- rnorm(n_fis, mean = 0, sd = 1) # fished environmental covariate
  
  # set simulation parameters
  # standardised effect sizes (slopes)
  b_restrictions <- 1
  b_age <- 1
  b_size <- 1
  b_enviro <- 1
  # residual standard deviation
  sigma <- 0.5
  # intercept (logged fish biomass)
  alpha <- 5
  # reserve effect (alpha)
  alpha_res <- 1
  
  # simulate log biomass from normal distribution
  # the effect of mpa, age and size on biomass only occurs in reserves, 
  # effect of restrictions only occurs in fished areas
  logbio_reserve <- rnorm(n_res, alpha + alpha_res + b_age*x_sage + b_size*x_size + b_enviro*x_enviro_r, sigma)
  logbio_fished <- rnorm(n_fis, alpha + b_restrictions*x_restrictions + b_enviro*x_enviro_f, sigma)
  
  # model fitting - recover parameters ------------------------------
  ################
  # simple linear model where there are 0s for age and size in reserves
  # create dataframe
  ################
  dat <- data.frame(logbio = c(logbio_reserve, logbio_fished), 
                    mpa = c(x_mpa, rep(0, n_fis)), 
                    restrictions = c(rep(0, n_res), x_restrictions),
                    age = c(x_sage, rep(0, n_fis)), 
                    size = c(x_size, rep(0, n_fis)),
                    enviro = c(x_enviro_r, x_enviro_f)) %>% 
    mutate(across(c(mpa, restrictions), factor))
  hist(dat$logbio)
  # frequentist
  fit_lm <- lm(logbio ~ mpa + restrictions + age + size + enviro, data = dat)
  fit_lm
  # bayesian
  fit_inla <- inla(logbio ~ mpa + restrictions + age + size + enviro,
                   data = dat,
                   family = 'gaussian',
                   control.predictor = list(compute = TRUE),
                   control.compute=list(return.marginals.predictor=TRUE))
  fit_inla$summary.fixed
  
  ################
  # joint model in stan
  ################
  # create list of data
  stanDat <- list(b = logbio_reserve,
                  b2 = logbio_fished,
                  res = length(logbio_reserve),
                  fis = length(logbio_fished),
                  mpa = x_mpa,
                  rest = x_restrictions,
                  ag = x_sage, 
                  si = x_size,
                  env = x_enviro_r,
                  env2 = x_enviro_f)
  
  fit_stan <- stan(file = 'scripts/models/joint-model.stan', data = stanDat, 
                   chains = 4, iter=2000, warmup=1000, control = list(adapt_delta = 0.9))
  fit_stan2 <- data.frame(parameter = row.names(summary(fit_stan)$summary), 
                          summary(fit_stan)$summary)
  
  ################
  # simple linear model where make mpa categories instead of assuming 0s for size and age in fished areas
  ################
  #dat2 <- dat %>% 
  # mutate(big_old = factor(ifelse(age >= quantile(age, 0.9) & size >= quantile(size, 0.9) & restrictions == 0 & mpa == 1, 1, 0)),
  #        average = factor(ifelse(age < quantile(age, 0.9) & age >= quantile(age, 0.7) & size >= quantile(size, 0.9) & size >= quantile(size, 0.7) & restrictions == 0 & mpa == 1, 1, 0)),
  #        small_new = factor(ifelse(age < quantile(age, 0.7) & size < quantile(size, 0.7) & restrictions == 0 & mpa == 1, 1, 0)))
  # mutate(big = factor(ifelse(size >= quantile(size, 0.5) & restrictions == 0 & mpa == 1, 1, 0)),
  #        old = factor(ifelse(age >= quantile(age, 0.9) & restrictions == 0 & mpa == 1, 1, 0)),
  #        small = factor(ifelse(size < quantile(size, 0.5) & restrictions == 0 & mpa == 1, 1, 0)),
  #        young = factor(ifelse(age < quantile(age, 0.9) & restrictions == 0 & mpa == 1, 1, 0)))
  
  #it_lm_categories <- lm(logbio ~ big + old + small + young + restrictions + enviro, dat2)
  #it_lm_categories
  
  out[[i]] <- data.frame(
    percent_reserves = loop[i]*100, 
    coefficient = c('alpha', 'b_mpa', 'b_restrictions', 'b_age', 'b_size', 'b_enviro'),
    truth = c(alpha-4, alpha_res, b_restrictions, b_age, b_size, b_enviro),
    linear_model = c(fit_inla$summary.fixed[1,1]-4, fit_inla$summary.fixed[2,1], fit_inla$summary.fixed[3,1], fit_inla$summary.fixed[4,1], fit_inla$summary.fixed[5,1], fit_inla$summary.fixed[6,1]),
    joint_model = c(fit_stan2[6,2]-4, fit_stan2[1,2], fit_stan2[2,2], fit_stan2[3,2], fit_stan2[4,2], fit_stan2[5,2]),
    upp_joint_model = c(fit_stan2[6,9]-4, fit_stan2[1,9], fit_stan2[2,9], fit_stan2[3,9], fit_stan2[4,9], fit_stan2[5,9]),
    low_joint_model = c(fit_stan2[6,5]-4, fit_stan2[1,5], fit_stan2[2,5], fit_stan2[3,5], fit_stan2[4,5], fit_stan2[5,5]),
    upp_linear_model = c(fit_inla$summary.fixed[1,5]-4, fit_inla$summary.fixed[2,5], fit_inla$summary.fixed[3,5], fit_inla$summary.fixed[4,5], fit_inla$summary.fixed[5,5], fit_inla$summary.fixed[6,5]),
    low_linear_model = c(fit_inla$summary.fixed[1,3]-4, fit_inla$summary.fixed[2,3], fit_inla$summary.fixed[3,3], fit_inla$summary.fixed[4,3], fit_inla$summary.fixed[5,3], fit_inla$summary.fixed[6,3])
  )
}

# plot to compare accuracy and precision of effect sizes ------------------------------

plot_dat <- do.call(rbind, out)

plot_dat %>% 
  ggplot() +
  #geom_point(aes(x = truth, y = coefficient), col = 'black', alpha = 0.5) +
  geom_vline(xintercept = 1, alpha = 0.5, lty = 'dashed') +
  geom_point(aes(x = linear_model, y = coefficient, color = 'Linear model'), size = 2, alpha = 0.5) +
  geom_errorbarh(aes(xmin = low_linear_model, xmax = upp_linear_model, y = coefficient, color = 'Linear model'), height = 0, size = 2, alpha = 0.5) +
  geom_point(aes(x = joint_model, y = coefficient, color = 'Joint model'), size = 0.5) +
  geom_errorbarh(aes(xmin = low_joint_model, xmax = upp_joint_model, y = coefficient, color = 'Joint model'), height = 0) +
  scale_color_manual(values = c("Linear model" = '#458B74', "Joint model" = 'coral3')) +
  xlab('Deviation from the truth (x = 1)') +
  ylab('Coefficient') +
  facet_wrap(~factor(percent_reserves)) +
  theme_classic() +
  theme(legend.title = element_blank())
ggsave('outputs/plot.png', height = 3, width = 7)
