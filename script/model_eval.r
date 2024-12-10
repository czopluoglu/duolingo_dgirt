
# Estimation Time 

get_elapsed_time(stanfit)

(sum(get_elapsed_time(stanfit))/4)/3600

################################################################################
# Analyze the parameter estimates

View(summary(stanfit, pars = c("mu_thetat",
                               "mu_thetac",
                               "sigma_thetat",
                               "sigma_thetac",
                               "mu_b",
                               "sigma_b",
                               "omega_P"), probs = c(0.025, 0.975))$summary)


View(summary(stanfit, pars = c("b"), probs = c(0.025, 0.975))$summary)
View(summary(stanfit, pars = c("person"), probs = c(0.025, 0.975))$summary)
View(summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary)
View(summary(stanfit, pars = c("pH"), probs = c(0.025, 0.975))$summary)

################################################################################
# Model Diagnostics

# Extract model summary with 95% credible intervals
model_summary <- as.data.frame(
  summary(stanfit, probs = c(0.025, 0.975))$summary
)

# Add a column for parameter types extracted from row names
model_summary$type <- gsub("\\[.*$", "", row.names(model_summary))

# Filter for person and item parameters of interest
model_summary <- model_summary[model_summary$type%in%c('pC','pH','person','b'),]

# ESS

  N    <- dim(model_summary)[[1]]
  iter <- dim(extract(stanfit)[[1]])[[1]]
  
  model_summary$n_eff_ratio <- ratio <- model_summary[,'n_eff'] / iter
  
  psych::describeBy(model_summary$n_eff_ratio,
                    model_summary$type,
                    mat=TRUE)[,c('group1','min')]
  
  psych::describeBy(model_summary$n_eff,
                    model_summary$type,
                    mat=TRUE)[,c('group1','mean','min','max')]
  
  ggplot(model_summary, aes(x = n_eff)) +
    geom_histogram(fill = "blue", color = "black", alpha = 0.7) +
    facet_wrap(~type, scales = "free_y",nrow=2) +
    labs(
      title = "Distribution of ESS Values by Parameter Type",
      x = "R-hat",
      y = "Frequency"
    ) +
    theme_minimal()
  
# Split Rhat
  
  ggplot(model_summary, aes(x = Rhat)) +
    geom_histogram(binwidth = 0.002, fill = "blue", color = "black", alpha = 0.7) +
    facet_wrap(~type, scales = "free_y",nrow=2) +
    labs(
      title = "Distribution of R-hat Values by Parameter Type",
      x = "R-hat",
      y = "Frequency"
    ) +
    theme_minimal()
  
  psych::describeBy(model_summary$Rhat,
                    model_summary$type,
                    mat=TRUE)[,c('group1','mean','min','max')]
  
  sum(model_summary$Rhat<1.01)/nrow(model_summary)

# Tree depth
  
  max_depth = 10
  sampler_params <- get_sampler_params(stanfit, inc_warmup=FALSE)
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  table(treedepths)
  n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  N = length(treedepths)
  100 * n / N

# E-BMFI
  
  sampler_params <- get_sampler_params(stanfit, inc_warmup=FALSE)
  e_bfmi <- c()
  for (n in 1:length(sampler_params)) {
    energies = sampler_params[n][[1]][,'energy__']
    numer = sum(diff(energies)**2) / length(energies)
    denom = var(energies)
    e_bfmi[n] = numer / denom
    print(sprintf('Chain %s: E-BFMI = %s', n, numer / denom))
  }
  
# Divergences
  
  sampler_params <- get_sampler_params(stanfit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
  N = length(divergent)
  100 * (n / N)

################################################################################
# Parameter Estimates
  
  # Person parameters
  
  theta <- summary(stanfit, pars = c("person"), probs = c(0.025, 0.975))$summary
  theta <- matrix(theta[,1],ncol=2,byrow=TRUE)
  psych::describe(theta)
  
  # Item parameters
  
  b <- summary(stanfit, pars = c("b"), probs = c(0.025, 0.975))$summary
  psych::describe(b[,1])
  
  
  # Probability of Items being compromised
  
  C_vec <- c()
  
  for(kk in 1:50){
    C_vec[kk] = unique(d_long[d_long$Item==kk,]$compromised)
  }
  
  
  pC <- as.numeric(summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary[,1])
  
  psych::describeBy(pC,C_vec,mat=TRUE)[,c('group1','n','mean','sd','min','max')]
  
  plot(density(pC[C_vec==0]),xlim=c(0,1),main="",ylim = c(0,8))
  points(density(pC[C_vec==1]),lty=2,type='l')
  
  
  require(pROC)
  
  auc(C_vec,pC)
  
  roc_analysis <- roc(response = C_vec,
                      predictor = pC)
  
  plot(1-roc_analysis$specificities,
       roc_analysis$sensitivities,
       xlim = c(0,1),ylim=c(0,1),
       xlab = 'False Positive Rate (1-Specificity)',
       ylab = 'True Positive Rate (Sensitivity)',
       type='l')
  
  my_thresholds <- seq(from=0.5,to=0.8,by=0.01)
  
  coords(roc_analysis, 
         my_thresholds, 
         input="threshold", 
         ret=c("threshold","specificity", "sensitivity"))

  # Probability of examinees having item preknowledge
  
  pH <- as.numeric(summary(stanfit, pars = c("pH"), probs = c(0.025, 0.975))$summary[,1])
  
  psych::describeBy(pH,d_wide$group,mat=TRUE)[,c('group1','n','mean','sd','min','max')]
  
  plot(density(pH[d_wide$group==0]),xlim=c(0,1),main="",ylim = c(0,8))
  points(density(pH[d_wide$group==1]),lty=2,type='l')
  
  auc(d_wide$group,pH)
  
  roc_analysis <- roc(response = d_wide$group,
                      predictor = pH)
  
  plot(1-roc_analysis$specificities,
       roc_analysis$sensitivities,
       xlim = c(0,1),ylim=c(0,1),
       xlab = 'False Positive Rate (1-Specificity)',
       ylab = 'True Positive Rate (Sensitivity)',
       type='l')
  
  my_thresholds <- seq(from=0.5,to=0.8,by=0.01)
  
  coords(roc_analysis, 
         my_thresholds, 
         input="threshold", 
         ret=c("threshold","specificity", "sensitivity"))
  