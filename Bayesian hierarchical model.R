######################################################################################################
####################### Hald Model for C. jejuni sources before and after ############################
######################################################################################################


install.packages("sourceR")
library(sourceR)

#setwd("H:/Antoine- Massey/Thesis/C. jejuni/attribution/attribution")

setwd("C:/Users/Antoine/Documents/Antoine- Massey/Thesis/C. jejuni/attribution/attribution")

db <- read.csv("hald_jejuni.csv")

priors <- list(a = 1, r = 1, theta = c(0.01, 0.00001))


res <- saBayes(formula = Human~Supplier.A+Supplier.B+Supplier.Others+Ruminants +Environmental.water+Wild.bird,
               time=~Time, location=~Location, type=~ST,
               data=db, priors = priors,
               alpha_conc = 1,
               likelihood_dist = "pois", n_iter = 5000,
               params_fix = list(type_effects = FALSE, r = FALSE),
               mcmc_params = list(save_lambda = TRUE, burn_in = 1,
                                  thin = 1, n_r = 200))

## End(Not run)

s <- summary(res, alpha = 0.05, burn_in = 1, thin = 1)

# print the median and HPD for the type effects

s$q

# print the median and HPD for the lambda's for each source
# for time 2, location A

s$lj_proportion$time1$locationA
s$lj_proportion$time2$locationA

