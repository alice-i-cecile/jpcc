# Demonic ritual ####
demonic_ritual <- function(effects, tra, sparse=TRUE){
  birth_years <- get_birth_years(tra, sparse)
  birth_index <- sapply(birth_years, get_birth_index, tra=tra, sparse=sparse)
  k <- -birth_index
  
  if (sparse){
    years <- sort(as.numeric(levels(tra$t)))
  } else {
    years <- sort(as.numeric(dimnames(tra)[[2]]))
  } 
  year_index <- match(names(effects$T), as.character(years))
  names(year_index) <- names(effects$T)
  
  if (sparse){
    ages <- sort(as.numeric(levels(tra$a)))
  } else {
    ages <- sort(as.numeric(dimnames(tra)[[3]]))
  }  
  age_index <- match(names(effects$A), as.character(ages))
  names(age_index) <- names(effects$A)
  
  k <- k[!is.na(k)]
  year_index <- year_index[!is.na(year_index)]
  age_index <- age_index[!is.na(age_index)]
  
  ridge_ids <- list(I=k, T=year_index, A=age_index)
  
  return(ridge_ids)
  
}


# Add a demonic intrusion ####
add_demons <- function (m=1, effects, tra, sparse=TRUE, form="multiplicative"){
  ridge_ids <- demonic_ritual(effects, tra, sparse)
  
  corrupted_effects <- effects
  
  if(form=="multiplicative")
  {
    corrupted_effects <- lapply(corrupted_effects, log)
  }
  
  
  # I = I +m*k
  corrupted_effects$I <- corrupted_effects$I + m*ridge_ids$I
  
  # T = T + m*C
  corrupted_effects$T <- corrupted_effects$T + m*ridge_ids$T
  
  # A = A - m*R
  corrupted_effects$A <- corrupted_effects$A - m*ridge_ids$A
  
  if(form=="multiplicative")
  {
    corrupted_effects <- lapply(corrupted_effects, exp)
  }
  
  corrupted_effects <- rescale_effects(corrupted_effects, form)
  return(corrupted_effects)
}

# Estimate a demonic intrusion ####
estimate_demons <- function (corrupted_effects, tra, sparse=TRUE, form="multiplicative"){
  ridge_ids <- demonic_ritual(corrupted_effects, tra, sparse)
  
  
  if(form=="multiplicative")
  {
    corrupted_effects <- lapply(corrupted_effects, log)
  }
  
  # Assuming normal random variation around the demonic line
  demonic_llh <- function (purified_effects)
  {
    sd_I <- sd(purified_effects$I)
    sd_T <- sd(purified_effects$T)
    sd_A <- sd(purified_effects$A)
    
    llh_I <- sum(dnorm(purified_effects$I, sd=sd_I, log=TRUE))
    llh_T <- sum(dnorm(purified_effects$T, sd=sd_T, log=TRUE))
    llh_A <- sum(dnorm(purified_effects$A, mean=mean(purified_effects$A), sd=sd_A, log=TRUE))
    
    llh <- sum(llh_I, llh_T, llh_A)
    
    return(llh)
  } 

  purify_probe <- function(m){
    purified_effects <- remove_demons(m, corrupted_effects, tra, sparse, form="additive")
    llh <- demonic_llh(purified_effects)
    return(-llh)
  }
  
  optimal_m <- optimize(f=purify_probe, interval=c(-1, 1))$minimum
  
  return(optimal_m)
}


# Remove a demonic intrusion ####
remove_demons <- function (m=1, corrupted_effects, tra, sparse=TRUE, form="multiplicative"){
  ridge_ids <- demonic_ritual(effects, tra, sparse)
  
  effects <- corrupted_effects
  
  if(form=="multiplicative")
  {
   effects <- lapply(effects, log)
  }
  
  # I = I + m*k
  effects$I <- effects$I - m*ridge_ids$I
  
  # T = T + m*C
  effects$T <- effects$T - m*ridge_ids$T
  
  # A = A - m*R
  effects$A <- effects$A + m*ridge_ids$A
  
  if(form=="multiplicative")
  {
    effects <- lapply(effects, exp)
  }
  
  effects <- rescale_effects(effects, form)
  return(effects)
}

# Full post-hoc correction of demonic intrusion ####
post_hoc_intercession <- function(corrupted_effects, tra, sparse=TRUE, form="multiplicative")
{
  optimal_m <- estimate_demons(corrupted_effects, tra, sparse, form)
  
  # Return suspiciously high trends
  if (abs(optimal_m) > 0.99){
    return(corrupted_effects)
  }
  
  print(paste("m of", optimal_m, "selected in choosing solution on likelihood ridge."))
  purified <- remove_demons(optimal_m, corrupted_effects, tra, sparse, form)
  
  
  return(purified)
}

# Testing ####
#effects <- ta_2$effects
#effects <- ita_2$effects
#tra <- ta_2$tra

#corrupted_effects <- add_demons(m=1, effects, tra, form="additive")
#purified <- remove_demons(m=1, corrupted_effects, tra, form="additive")

#corrupted_effects <- add_demons(m=1, effects, tra, form="multiplicative")
#purified <- remove_demons(m=1, corrupted_effects, tra, form="multiplicative")

#purified <-  remove_demons(m=-0.00633, effects, tra, form="multiplicative")
#purified <-  post_hoc_intercession(effects, tra, form="multiplicative")

#purified <- post_hoc_intercession(corrupted_effects, tra, form="multiplicative")


#plot(ta_2$effects$A)
#plot(purified$T)
#plot(purified$A)