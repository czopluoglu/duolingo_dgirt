################################################################################
require(MASS)
require(MBESS)
require(here)
require(psych)
require(pROC)
set.seed(7122024)
################################################################################
# Generate model parameters

N  <- 500     # sample size
n  <- 50      # number of items
pe <- 0.10    # proportion of examinees with item preknowledge
pi <- 0.30    # proportion of compromised items

# Generate the binary status of examinee item preknowledge
  # 1: examinee has item preknowledge
  # 0: examinee has item preknowledge

  tmp <- runif(N,0,1)
  H  <- ifelse(tmp<=quantile(tmp,pe),1,0)
  H
  table(H)

# Generate the binary status of item compromise
  # 1: item is compromised
  # 0: item is not compromised
  
  tmp <- runif(n,0,1)
  C  <- ifelse(tmp<=quantile(tmp,pi),1,0)
  C
  table(C)

# Generate item difficulty parameters
  
  b <- rnorm(n,0,1)
  b <- (b-mean(b))/sd(b)
  b
  describe(b)

# Generate person parameters
  
  mu_t    <- 0      # mean of true latent trait parameters
  mu_c    <- 3      # mean of cheating latent trait parameters
  sigma_t <- 1      # standard dev. of true latent trait parameters
  sigma_c <- 1.25   # standard dev. of cheating latent trait parameters
  corr    <- 0.8    # covariance between cheating and true latent trait parameters

        1/(1+exp(-mu_t)) # probability of correct for average item without preknowledge
        1/(1+exp(-mu_c)) # probability of correct for average item with preknowledge
          
        # some rough idea about effect size, item preknowledge effect
        # Odds ratio
          exp(mu_c)/exp(mu_t)
          
      
  th <- mvrnorm(N,
                mu = c(mu_t,mu_c),
                Sigma = matrix(c(sigma_t,corr,corr,sigma_c),2,2))
      
  theta_t <- th[,1] 
  theta_c <- th[,2]
      
  describe(theta_t)
  describe(theta_c)
  describe(theta_c - theta_t)
  cor(theta_t,theta_c)
      
###############################################################################
# Generate observed responses

r <- matrix(nrow=N,ncol=n)

for(j in 1:N){
  for(i in 1:n){
    
    p_t <- exp(theta_t[j] - b[i])/(1+exp(theta_t[j] - b[i]))
    p_c <- exp(theta_c[j] - b[i])/(1+exp(theta_c[j] - b[i]))
    
    if(H[j] == 1 & C[i] == 1){
      r[j,i] = rbinom(1,1,p_c)
    } else {
      r[j,i] = rbinom(1,1,p_t)
    }
    
  }
}

colnames(r) <- paste0("Y",1:n)
###############################################################################

# Convert it to data frame and add group membership and a unique ID

d       <- as.data.frame(r)
d$group <- H
d$id    <- 1:nrow(d)

# Check the data

head(d)

# Reshape it to long format (for plotting purposes)

d.long <- reshape(data        = d,
                  idvar       = 'id',
                  varying     = list(colnames(r)[1:n]),
                  timevar     = "Item",
                  times       = 1:n,
                  v.names     = c("R"),
                  direction   = "long")

d.long <- d.long[!is.na(d.long$R),]

# Add item status

d.long$compromised <- NA

for(i in 1:n){
  
  d.long[d.long$Item==i,]$compromised = C[i]
  
}


d.long <- d.long[order(d.long$Item),]

describeBy(d.long$R,list(d.long$group,d.long$compromised),mat=TRUE)


################################################################################

write.csv(d.long,'./data/simdata_long.csv',
          row.names = FALSE)  

write.csv(d,'./data/simdata_wide.csv',
          row.names = FALSE)  




