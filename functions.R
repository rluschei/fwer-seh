source("libraries.R")

corr <- function(n, case){
  n_P2T <- sum(n[2:3,"T"])
  n_P2C <- sum(n[2:3,"C"])
  if(case == "nested"){
    n_P1T <- sum(n[,"T"])
    n_P1C <- sum(n[,"C"])
    sum <- (n[2,"T"] + n[3,"T"])/(n_P1T * n_P2T) + (n[2,"C"] + n[3,"C"])/(n_P1C * n_P2C)
  }
  if(case == "overlap"){
    n_P1T <- sum(n[1:2,"T"])
    n_P1C <- sum(n[1:2,"C"])
    sum <- n[2,"T"]/(n_P1T * n_P2T) + n[2,"C"]/(n_P1C * n_P2C)
  }
  sum / (sqrt(V(1,n,case) * V(2,n,case)))
}

V <- function(i, n, case){
  if(i == 1) {
    if(case == "nested") return (1/sum(n[,"T"]) + 1/(sum(n[,"C"])))
    if(case == "overlap") return (1/sum(n[1:2,"T"]) + 1/(sum(n[1:2,"C"])))
  }
  if(i == 2){
    return (1/sum(n[2:3,"T"]) + 1/(sum(n[2:3,"C"])))
  } 
}

sizes <- function(N, pi, delta, alloc) {
  if(alloc %in% c("A", "B", "C")){
    n <- rep(4, 3) + rmultinom(1, N - 12, pi)
    n <- c(n * delta, n * (1-delta))
  }
  else{
    p <- c(delta * pi, (1 - delta) * pi)
    n <- rep(2, 6) + rmultinom(1, N - 12, p)
  }
  mat <- matrix(c(n[1:3], n[4:6]), nrow = 3, ncol = 2)
  dimnames(mat) <- list(c("Pop1", "Pop2", "Pop3"), c("T", "C"))
  return(mat)
}


fwer.anova <- function(c, nu, rho, N){
  Sigma <- matrix(rep(1,4), nr = 2) 
  Sigma[1,2] <- rho
  Sigma[2,1] <- rho
  1 - pmvt(upper = c, delta = nu, corr = Sigma, df = N-6)[1]
}

crit.anova <- function(rho, alpha, N){
  f <- function(crit){fwer.anova(c(crit, crit), c(0, 0), rho, N) - alpha}
  uniroot(f, interval = c(0,25))$root
}

fwer.marg <- function(c, rho, df){
  Sigma <- matrix(rep(1,4), nr = 2) 
  Sigma[1,2] <- rho
  Sigma[2,1] <- rho
  1 - pmvt(upper = c, corr = Sigma, df = df)[1]
}

crit.marg <- function(rho, df, alpha){
  f <- function(crit){fwer.marg(c(crit, crit), rho, df) - alpha}
  uniroot(f, interval = c(0,25))$root
}

gen.Theta <- function(EHF, pi, case, rate) {
  if(rate == "fwer"){
    Theta <- rep(0, 3)
    Theta[2] <- EHF * runif(1, -1, 1)
    Theta[3] <- -pi[2]/pi[3] * Theta[2]
    
    if(case == "nested"){
      Theta[1] <- (-pi[2] * Theta[2] - pi[3] * Theta[3])/pi[1]
    }
    if(case == "overlap"){
      Theta[1] <- -pi[2]/pi[1] * Theta[2]
    }
    return(Theta)
  }
  if(rate == "power"){
    while (TRUE) {
      Theta <- EHF * runif(3, -1, 1)
      
      if (case == "nested") {
        if (sum(pi * Theta) > 0 && sum(pi[2:3] * Theta[2:3]) > 0) {
          return(Theta)
        }
      }
      if (case == "overlap") {
        if (sum(pi[1:2] * Theta[1:2]) > 0 && sum(pi[2:3] * Theta[2:3]) > 0) {
          return(Theta)
        }
      } 
    }
  }
}

perform_test <- function(n, N, case, test, pi, alpha, sd, Nboot, Theta, EHF, CHF, delta, mu){
  if(test == "marg+t") return(marg.test(n, case, alpha, Theta, CHF, sd, N, mu))
  if(test %in% c("anova+boot", "marg+boot", "marg+shr+boot", "strat+boot")) return(boot.test(alpha, Nboot, test, Theta, sd, CHF, n, case, N, delta, mu))
  if(test == "anova+t") return(anova.test(n, case, alpha, sd, Theta, CHF, mu, N))
  if(test == "rem") return(rem(n, case, alpha, Theta, CHF, sd))
  if(test == "marg+shr+t") return(marg.test(n, case, alpha, Theta, CHF, sd, N, mu, shr = T))
}

gen.pops.list <- function(case){
  list(
    if (case == "nested") 1:3 else 1:2,
    2:3
  )
}

mean_diff <- function(pops, n, mu, type = "pw"){
  if(type == "pw"){
    mt <- sum(n[pops,"T"] * mu[pops,"T"]) / sum(n[pops,"T"])
    mc <- sum(n[pops,"C"] * mu[pops,"C"]) / sum(n[pops,"C"])
    return(mt - mc)
  } 
  if(type == "sw"){
    nS <- n[pops,"T"] + n[pops,"C"]
    mdS <- mu[pops,"T"] - mu[pops,"C"]
    return(sum(nS/sum(n[pops,]) * mdS))
  } 
}


t.anova <- function(case, n, sd, mu) {  
  pops_list <- gen.pops.list(case)
  sapply(pops_list, function(pops) {
    num <- mean_diff(pops, n, mu)
    denom <- sd * sqrt(1 / sum(n[pops,"T"]) + 1 / sum(n[pops,"C"]))
    num / denom
  })
}


gen.mu.true <- function(Theta, CHF) {
  mu.true <- matrix(0, nrow = 3, ncol = 2)
  rownames(mu.true) <- c("Pop1", "Pop2", "Pop3")
  colnames(mu.true) <- c("T", "C")
  
  for (p in 1:3) {
    mu.true[p,"T"] <- (p - 1) * CHF + Theta[p]
    mu.true[p,"C"] <- (p - 1) * CHF 
  }
  
  mu.true
}

gen.mu.est <- function(mu.true, n, sd){
  matrix(
    rnorm(length(mu.true), mean = as.vector(mu.true), sd = sd / sqrt(as.vector(n))),
    nrow = nrow(mu.true),
    dimnames = dimnames(mu.true)
  )
}

anova.test <- function(n, case, alpha, sd, Theta, CHF, mu, N){
  rho <- corr(n, case)
  t <- t.anova(case, n, sd, mu)
  c <- rep(crit.anova(rho, alpha, N), 2)
  if(any(t > c)) return(1)
  return(0)
}

pop.var <- function(pops, Tmt, mu, n, sd){
  sum(n[pops,Tmt] * mu[pops,Tmt]^2 / sum(n[pops,Tmt])) - sum(n[pops,Tmt] * mu[pops,Tmt] / sum(n[pops,Tmt]))^2 + sd^2
}

t.marg <- function(case, n, sd, mu, boot = F, shr = F) {
  pops_list <- gen.pops.list(case)
  
  t_vals <- numeric(2)
  df_vals <- numeric(2)
  
  mu.orig <- mu
  for (i in 1:2) {
    pops <- pops_list[[i]]
    num <- mean_diff(pops, n, mu.orig)
    if(shr) mu <- js(mu.orig, n, sd, pops)
    varT <- pop.var(pops, "T", mu, n, sd)
    varC <- pop.var(pops, "C", mu, n, sd)
    nT <- sum(n[pops, "T"])
    nC <- sum(n[pops, "C"])
    denom <- sqrt(varT / nT + varC / nC)
    t_vals[i] <- num / denom
    if(boot == F) df_vals[i] <- floor((varT/nT + varC/nC)^2 / ((varT/nT)^2 / (nT - 1) + (varC/nC)^2 / (nC - 1)))
  }
  
  return(list(t = t_vals, df = df_vals))
}

gen.sd.est <- function(sd, N){
  sqrt(rchisq(1, N - 6) * sd^2 / (N - 6))
}

marg.test <- function(n, case, alpha, Theta, CHF, sd, N, mu, shr = F){
  marg <- t.marg(case, n, sd, mu, shr)
  t <- marg[[1]]
  df <- marg[[2]]
  rho <- corr(n, case)
  c <- c(crit.marg(rho, df[1], alpha), crit.marg(rho, df[2], alpha))
  if (any(t > c)) return(1)
  return(0)
}


t.swdiff <- function(case, n, N, sd, mu, delta) {
  pops_list <- gen.pops.list(case)
  t_vals <- numeric(2)
  
  for (i in 1:2) {
    pops <- pops_list[[i]]
    num <- mean_diff(pops, n, mu, "sw")
    
    d <- mu[pops, "T"] - mu[pops,"C"]
    nP <- sum(n[pops,])
    nS <- rowSums(n[pops, , drop = FALSE])
    piS <- nS/nP
    Sigma <- nP * (diag(piS) - outer(piS, piS))
    V1 <- d %*% Sigma %*% d / nP^2
    
    V2 <- sd^2/nP^2 * sum((1/delta[pops] + 1/(1-delta[pops])) * nS)
    
    denom <- sqrt(V1 + V2)
    
    t_vals[i] <- num / denom
  }
  t_vals
}

r <- function(pops, n) {
  nS <- rowSums(n[pops, , drop = FALSE])
  piS <- nS / sum(nS)
  rval1 <- setNames(rep(0,3), c("Pop1", "Pop2", "Pop3"))
  rval2 <- setNames(rep(0,3), c("Pop1", "Pop2", "Pop3"))
  rval1[names(piS)] <- piS
  rval2[names(piS)] <- -piS
  append(rval1, rval2)
}

t.boot <- function(case, n, sd, Nboot, test, N, delta, mu) {
  pops_list <- gen.pops.list(case)
  df <- data.frame(mu = c(mu[,"T"], mu[,"C"]), r1 = r(pops_list[[1]], n), r2 = r(pops_list[[2]], n))
  mu.null <- lm("mu ~ 0 + r1 + r2", df)$residuals
  
  p <- c((n[,"T"] + n[,"C"] - 4)/sum(n[,"T"] + n[,"C"] - 4))
  
  base::t(replicate(Nboot, {
    n <- sizes(N, p, delta, "A")
    mu.boot <- matrix(
      rnorm(length(n), mean = mu.null, sd = sd / sqrt(as.vector(n))),
      nrow = nrow(n),
      dimnames = dimnames(n)
    )
    if(test == "anova+boot") return(t.anova(case, n, sd, mu.boot))
    if(test == "marg+boot") return(t.marg(case, n, sd, mu.boot, boot = T)[[1]])
    if(test == "marg+shr+boot") return(t.marg(case, n, sd, mu.boot, boot = T, shr = T)[[1]])
    if(test == "strat+boot") return(t.swdiff(case, n, N, sd, mu.boot, delta))
  }))
}

boot.test <- function(alpha, Nboot, test, Theta, sd, CHF, n, case, N, delta, mu){
  if(test == "anova+boot") t <- t.anova(case, n, sd, mu)
  if(test == "marg+boot") t <- t.marg(case, n, sd, mu, boot = T)[[1]]
  if(test == "marg+shr+boot") t <- t.marg(case, n, sd, mu, boot = T, shr = T)[[1]]
  if(test == "strat+boot") t <- t.swdiff(case, n, N, sd, mu, delta)
  t_boot <- t.boot(case, n, sd, Nboot, test, N, delta, mu)
  max <- apply(t_boot, 1, max)
  p1 <- sum(max >= t[1])/Nboot
  p2 <- sum(max >= t[2])/Nboot
  p <- c(p1, p2)
  if(any(p <= alpha, na.rm = T)) return(1)
  return(0)
}

crit.boot <- function(alpha, Nboot, test, n, case, N, delta, sd, mu){
  t_boot <- t.boot(case, n, sd, Nboot, test, N, delta, mu)
  max <- apply(t_boot, 1, max)
  quantile(max, probs = 1 - alpha)
}

js <- function(mu, n, sd, pops) {
  mu.js <- mu
  k <- length(pops)
  
  for (j in 1:ncol(mu)) {
    mu_j <- mu[pops, j]
    mu_bar <- mean(mu_j)
    s2 <- sum(mu_j^2 * n[pops,j])/sd^2
    shrinkage <- if (s2 > 0) (k - 2) / s2 else 0
    mu.js[pops, j] <- 1 - shrinkage * (mu_j - mu_bar) + mu_bar
  }
  
  rownames(mu.js) <- rownames(mu)
  colnames(mu.js) <- colnames(mu)
  mu.js
}

generate_sample <- function(n, mu, sd){
  n <- ceiling(n)
  tibble(
    Pop = factor(c(
      rep(1:3, times = n[,"T"]),  
      rep(1:3, times = n[,"C"])
    )),
    Tmt = factor(c(rep("T", sum(n[,"T"])), rep("C", sum(n[,"C"]))),
                 levels = c("T", "C")),
    Resp = as.vector(unlist(mapply(function(p, t) {
      rnorm(n[p, t], mean = mu[p,t], sd = sd)
    }, rep(1:3, 2), rep(c("T", "C"), each = 3))))
  )
}

rem <- function(n, case, alpha, Theta, CHF, sd){
  mu.true <- gen.mu.true(Theta, CHF)
  sample <- generate_sample(n, mu.true, sd)
  
  if(case == "nested"){
    sample1 <- sample
    sample2 <- sample %>% filter(Pop != 1)
  }
  if(case == "overlap"){
    sample1 <- sample %>% filter(Pop != 3)
    sample2 <- sample %>% filter(Pop != 1)
  }
  sample1$Tmt <- relevel(sample1$Tmt, ref = "C")
  sample2$Tmt <- relevel(sample2$Tmt, ref = "C")
  
  mod1 <- lmer(Resp ~ Tmt + (1 + Tmt | Pop), data = sample1)
  mod2 <- lmer(Resp ~ Tmt + (1 + Tmt | Pop), data = sample2)
  
  # Extract two-tailed p-values
  p_two_tailed <- c(
    summary(mod1)$coefficients["TmtT", "Pr(>|t|)"],
    summary(mod2)$coefficients["TmtT", "Pr(>|t|)"]
  )
  
  # Calculate right-sided p-values
  p <- p_two_tailed / 2
  coef_signs <- c(
    summary(mod1)$coefficients["TmtT", "Estimate"],
    summary(mod2)$coefficients["TmtT", "Estimate"]
  )
  p <- ifelse(coef_signs > 0,  p, 1-p)
  
  if(any(p <= alpha/2)) return(1)
  return(0)
}

