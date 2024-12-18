require(cmdstanr)
################################################################################

# Import data (long format)

d_long <- read.csv('./data/simdata_long.csv')
d_wide <- read.csv('./data/simdata_wide.csv')

################################################################################
# Fit DG-IRT with a single chain
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

################################################################################

# Extract the estimates from single chain initialization to feed into
# multi-chain estimation as start parameters

b_start <- as.vector(fit_init$summary("b")$mean)
b_start

theta <- data.frame(thetat = fit_init$summary("person")$mean[1:500],
                    thetac = fit_init$summary("person")$mean[501:1000])
rownames(theta) <- NULL
colnames(theta) <- NULL
head(theta)

H <- as.vector(fit_init$summary("pH")$mean)
head(H)

C <- as.vector(fit_init$summary("pC")$mean)
head(C)


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


# Save the model object with all parameters for future use

fit$save_object(file = "./do_no_upload/model_fit.RDS")










