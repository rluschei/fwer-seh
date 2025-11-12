source("libraries.R")
source("functions.R")

tar_option_set(
  controller = crew_controller_local(workers = detectCores() - 5)
)

sim <- function(N, EHF, CHF, test, rate, seed, alloc, case, alpha = 0.025, sd = 0.5, 
                Nsim = 10^3, Nboot = 10^3) {
  if (EHF == 0 & rate == "power") return(NA)
  
  set.seed(seed)
  pi <- runif(3); pi <- pi / sum(pi)

  Theta <- gen.Theta(EHF, pi, case, rate)
  
  delta <- switch(alloc,
                  "A" = rep(0.5, 3),
                  "B" = rep(0.33, 3),
                  "C" = c(0.5, 0.33, 0.5),
                  "D" = rep(0.5, 3)
  )

  results <- map_dbl(1:Nsim, function(i) {
    n <- sizes(N, pi, delta, alloc)
    if(alloc == "D") delta <- n[,1] / rowSums(n)
    
    mu.true <- gen.mu.true(Theta, CHF)
    mu.est <- gen.mu.est(mu.true, n, sd)
    sd.est <- gen.sd.est(sd, N)
    perform_test(n, N, case, test, pi, alpha, sd.est, Nboot, Theta, EHF, CHF, delta, mu.est)
  })
  
  mean(results)
}

create_row <- function(params){
  rval <- sim(params$N, params$EHF, params$CHF, params$test, params$rate, 
              params$seed, params$alloc, params$case)
  tibble(
    N = params$N, EHF = params$EHF, CHF = params$CHF, test = params$test, 
    rate = params$rate, alloc = params$alloc, case = params$case, seed = params$seed,
    value = rval
  )
}

list(
  
  tar_target(
    params,
    expand_grid(
      N = c(250, 500, 1000),
      EHF = c(0, 1, 10),
      CHF = c(0, 1, 10),
      test = c("anova+t", "anova+boot", "marg+t", "marg+boot", "marg+shr+boot", 
               "strat+boot", "rem"),
      rate = c("fwer", "power"),
      alloc = c("A", "B", "C", "D"),
      case = c("nested", "overlap"),
      seed = 1:100
    )
  ),
  
  tar_target(
    results,
    create_row(params),
    pattern = map(params)
  )
  
)
