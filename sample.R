# Compute the positive part of a real number x, which is $\max(x, 0)$.
positive_part <- function(x) {ifelse(x > 0, x, 0)}

# This function generates n data ponints from some unimodal population.
# Input: ----------------------------------------------------
# n: sample size;
# mu: the mode of the population, default value is 0.
# skewness: the parameter that reflects the skewness of the distribution, note it is not
#           the exact skewness defined in statistics text book, default value is 0.
# tailedness: the parameter that reflects the tailedness of the distribution, note it is
#             not the exact kurtosis defined in text book, default value is 0.

# When all arguments take their default values, the data will be generated from standard 
# normal distribution.

random_sample <- function(n, mu = 0, skewness = 0, tailedness = 0){
  sigma = 1
  
  # The sampling scheme resembles the rejection sampling. For each step, an initial data point
  # was proposed, and it will be rejected or accepted based on the weights determined by the
  # skewness and tailedness of input. 
  reject_skewness <- function(x){
    scale = 1
    # if `skewness` > 0 (means data are right-skewed), then small values of x will be rejected
    # with higher probability.
    l <- exp(-scale * skewness * x)
    l/(1 + l)
  }
  
  reject_tailedness <- function(x){
    scale = 1
    # if `tailedness` < 0 (means data are lightly-tailed), then big values of x will be rejected with
    # higher probability.
    l <- exp(-scale * tailedness * abs(x))
    l/(1 + l)
  }
  
  # w is another layer option to control the tailedness, the higher the w is, the data will be
  # more heavily-tailed. 
  w = positive_part((1 - exp(-0.5 * tailedness)))/(1 + exp(-0.5 * tailedness))
  
  filter <- function(x){
    # The proposed data points will be accepted only if it satified the following condition, 
    # in which way we controlled the skewness and tailedness of data. (For example, the 
    # proposed data point will be rejected more frequently if it has higher skewness or
    # tailedness.)
    accept <- runif(length(x)) > reject_tailedness(x) * reject_skewness(x)
    x[accept]
  }
  
  result <- filter(mu + sigma * ((1 - w) * rnorm(n) + w * rt(n, 5)))
  # Keep generating data points until the length of data vector reaches n.
  while (length(result) < n) {
    result <- c(result, filter(mu + sigma * ((1 - w) * rnorm(n) + w * rt(n, 5))))
  }
  result[1:n]
}

multimodal <- function(n, Mu, skewness = 0, tailedness = 0) {
  # Deal with the bimodal case.
  mumu <- as.numeric(Mu %*% rmultinom(n, 1, rep(1, length(Mu))))
  mumu + random_sample(n, skewness = skewness, tailedness = tailedness)
}