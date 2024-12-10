require(cmdstanr)
require(rstan)
require(psych)
require(bayesplot)
require(pROC)
require(irtoys)
################################################################################

# Import data (long format)

d_long <- read.csv('./data/simdata_long.csv')
d_wide <- read.csv('./data/simdata_wide.csv')

################################################################################
# #Step 1: Generate starting parameters

# We fit a Rasch model first to estimate item difficulty parameters 
# I will use these values later to explore label swapping in different chains

require(TAM)

initial_fit <- TAM::tam.mml(resp=d_wide[,1:50])

b_est <- initial_fit$xsi$xsi
b_est <- (b_est-mean(b_est))/sd(b_est)

################################################################################
#Step 2: Fit DG-IRT with a single chain
# The goal of this is to obtain initial start values
# We will provide the values obtained from single chain as starting values
# when re-fitting model again later with multiple chains

# Input Data

data_resp <- list(
  I              = length(unique(d_long$Item)),
  J              = length(unique(d_long$id)),
  n_obs          = nrow(d_long),
  p_loc          = d_long$id,
  i_loc          = d_long$Item,
  Y              = d_long$R
)

# Compile the model syntax

mod <- cmdstan_model('./script/dgirt.stan')

# Fit the model using a single chain with a small number of iterations

fit_init <- mod$sample(
  data            = data_resp,
  seed            = 1234,
  chains          = 1,
  iter_warmup     = 150,
  iter_sampling   = 150,
  refresh         = 10,
  adapt_delta     = 0.99)


fit_init$cmdstan_summary()

stanfit_init <- rstan::read_stan_csv(fit_init$output_files())


b_start <- summary(stanfit_init, pars = c("b"), probs = c(0.025, 0.975))$summary[,1]
plot(b_est,b_start)
################################################################################

b_start <- as.vector(summary(stanfit_init, pars = c("b"), probs = c(0.025, 0.975))$summary[,1])
theta <- data.frame(thetat = summary(stanfit_init, pars = c("person"), probs = c(0.025, 0.975))$summary[seq(1,1000,2),1],
                    thetac = summary(stanfit_init, pars = c("person"), probs = c(0.025, 0.975))$summary[seq(2,1000,2),1])
rownames(theta) <- NULL
colnames(theta) <- NULL
H <- as.numeric(summary(stanfit_init, pars = c("pH"), probs = c(0.025, 0.975))$summary[,1])
C <- as.numeric(summary(stanfit_init, pars = c("pC"), probs = c(0.025, 0.975))$summary[,1])

# Fit the model

start <- list(b      = b_start,
              person = theta,
              pH     = H,
              pC     = C)

fit <- mod$sample(
  data            = data_resp,
  seed            = 1234,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 500,
  iter_sampling   = 1000,
  refresh         = 10,
  init            = list(start,start,start,start),
  adapt_delta     = 0.99)


# Compile the output files into an rstan object

fit$cmdstan_summary()

stanfit <- rstan::read_stan_csv(fit$output_files())








