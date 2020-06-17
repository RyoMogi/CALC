# ===========================================================================================
# Code to calculate both period and cohort EYWC taking into account of mortality using HFD and HMD.
# Ryohei Mogi, rmogi@ced.uab.es
# with Jessica Nisen, Vladimir Canudas-Romo
# 2019/5/30
# ===========================================================================================
# upd: 2020/06/04

library(tidyverse)


### ---- Country selecter ----
country  <- c("JPN", "ESP", "NLD", "HUN", "DNK", "SWE",
              "CZE", "EST", "LTU", "USA", "BLR")
country2 <- c("Japan", "Spain", "The Netherlands", "Hungary", "Denmark", "Sweden", 
              "Czechia", "Estonia", "Lithuania", "The US", "Belarus")

outcome <- c()


for(i in 1:length(country)){
  
  Names  <- country[i]
  Names2 <- country2[i]
  
  # Period data
  A1_per  <- read.table(paste("Data/", Names, "pft.txt", sep = ""), header = TRUE, fill = TRUE, skip = 2)
  A1M_per <- read.table(paste("Data/", Names, ".fltper_1x1.txt", sep = ""), header = TRUE, fill = TRUE, skip = 2)
  
  # Cohort data
  A1_bc  <- read.table(paste("Data/", Names, "cft.txt", sep = ""), header = TRUE, fill = TRUE, skip = 2)
  
  Exp    <- read.table(paste("Data/", Names, ".Exposures_1x1.txt", sep = ""), header = TRUE, fill = TRUE, skip = 2)
  Death  <- read.table(paste("Data/", Names, ".Deaths_lexis.txt", sep = ""), header = TRUE, fill = TRUE, skip = 2)
  
  ### ---- Period EYWC ----
  over50 <- c("50", "51", "52", "53", "54", "55+")
  `%out%` = Negate(`%in%`)
  
  # calculate period EYWC 
  A1_per_sel <- A1_per %>% 
    filter(x %out% over50,
           Year == 2015) %>%
    mutate(q1x = as.numeric(as.character(q1x)),
           px = 1 - q1x,
           lx = cumprod(px)) %>% 
    select(Year, x, q1x, px, lx)
  
  A1M_per16 <- A1M_per %>% 
    filter(Age %out% "110+") %>%
    mutate(Age = as.numeric(as.character(Age))) %>% 
    filter(Age %in% 12:49,
           Year == 2015) %>% 
    select(Year, Age, qx)
  
  lx_per <- c(1, A1_per_sel$lx[-nrow(A1_per_sel)])
  
  Lx_per <- (lx_per[1:37] + lx_per[2:38]) * 0.5
  Lx_per <- c(Lx_per, lx_per[38])
  
  EYWC     <- round(sum(Lx_per), 2)
  adj_EYWC <- round(sum(exp(- cumsum(A1_per_sel$q1x + A1M_per16$qx))), 2)
  
  ### ---- Cohort EYWC ----
  
  #### Death rate
  period2cohort <- function(data){
    
    n <- ncol(data)
    
    bc <- c()
    for(i in 2:n){
      row <- c()
      row <- c(data[i, - c(1:i - 1)], rep(NA, i-1))
      bc  <- rbind(bc, row)
    }
    data <- as.matrix(data)
    bc <- rbind(data[1, ], bc)
    
    years <- as.numeric(as.character(colnames(bc)))
    colnames(bc) <- years - 12
    rownames(bc) <- NULL
    
    bc <- as.data.frame(bc)
    bc <- bc[1:38, ]
    return(bc)
  }
  
  ## Convert exposure data from period to cohort
  # convert long to wide
  Exp_mat <- Exp %>% 
    select(Year, Age, Female) %>% 
    filter(Age %out% "110+") %>% 
    mutate(Age = as.numeric(as.character(Age))) %>% 
    filter(Age %in% 12:49) %>% 
    spread(key = Year, value = Female) %>% 
    select(-Age)
  # convert period to cohort
  Exp_matBC <- period2cohort(data = Exp_mat)
  Exp_matBC66 <- Exp_matBC[, "1966"] 
  
  
  ## Convert death data from period to cohort
  Death_matBC <- Death %>% 
    select(Cohort, Age, Female) %>% 
    filter(Age %out% "110+") %>% 
    mutate(Age = as.numeric(as.character(Age))) %>% 
    filter(Age %in% 12:49) %>%
    group_by(Cohort, Age) %>% 
    summarise(sum(Female)) %>% 
    spread(key = Cohort, value = "sum(Female)") %>% 
    select(-Age)
  
  Death_matBC66 <- as.matrix(Death_matBC[, "1966"])
  
  ## Calculate life table
  Mx1966 <- Death_matBC66 / unlist(Exp_matBC66)
  
  # Death Probabilities between ages x and x+n
  n  <- rep(1, 38)
  ax <- rep(0.5, 38)
  qx1966     <- (n * Mx1966) / (1 + (n - ax) * Mx1966)
  
  
  #### Fertility rate
  A1_bc66 <- A1_bc %>% 
    filter(x %out% over50,
           Cohort == 1966) %>%
    mutate(q1x = as.numeric(as.character(q1x)),
           px = 1 - q1x,
           lx = cumprod(px)) %>% 
    select(Cohort, x, q1x, px, lx)
  
  lx_bc <- c(1, A1_bc66$lx[-nrow(A1_bc66)])
  Lx_bc <- (lx_bc[1:37] + lx_bc[2:38]) * 0.5
  Lx_bc <- c(Lx_bc, lx_bc[38])

  #### Calculate cohort EYWC
  EYWC_bc     <- round(sum(Lx_bc), 2)
  adj_EYWC_bc <- round(sum(exp(- cumsum(A1_bc66$q1x - qx1966))), 2)
  
  
  ### ---- CALC ----

  CALCFunc <- function(px){
    
    # px matrix to lx
    lx <- apply(px, 2, cumprod)
    lx <- rbind(rep(1, ncol(px)), lx)
    
    # lx to CAL lx
    CALlx <- c()
    for(i in 1:ncol(px)){
      order <- 38:1
      CALlx[i] <- lx[order[i], i]
    }
    
    #CALlx <- apply(px, 2, prod, na.rm = T)
    
    CALlx <- rev(CALlx)
    
    # CAL lx to CAL Lx
    CALLx <- (CALlx[1:37] + CALlx[2:38]) * 0.5
    CALLx <- c(CALLx, CALlx[38])
    
    # CALC
    CALC <- round(sum(CALLx), 2)
    
    return(CALC)
  }
  
  ## Function to create a matrix of lx or px
  widelxpx <- function(data, target){
    px <- data %>% 
      as.data.frame() %>% 
      #filter(Cohort >= Y1) %>% # select year from the same year
      filter(x %out% over50) %>%              # select age (12- to 49)
      mutate(q1x = as.numeric(as.character(q1x)),
             px = 1 - q1x) %>% 
      select(Cohort, x, px)
    
    # create a matrix of px
    px_wide <- px %>% 
      mutate(px = ifelse(x %in% c("12-", "13", "14", "15") & is.na(px), 1, px)) %>% 
      spread(key = Cohort, value = px) %>%
      select(-x) %>% 
      as.matrix()
    
    lx_wide <- matrix(NA, ncol = ncol(px_wide), nrow = nrow(px_wide))
    lx_wide[1, ] <- 1
    for(i in 1:(nrow(lx_wide)-1)){
      lx_wide[i+1, ] <- lx_wide[i, ] * px_wide[i, ]
    }
    
    colnames(lx_wide) <- colnames(px_wide)
    
    if(target == "lx"){
      outcome <- lx_wide
    } else {
      outcome <- px_wide
    }
    
    return(outcome)
  }
  
  ## For country 1
  px_bc <- widelxpx(data = A1_bc, target = "px")
  
  if(any(colnames(px_bc) == "1966")){
    px_1966 <- px_bc[, "1966"]
  } else {
    px_1966 <- NA
  }
  
  ## the position of the maximum completed birth cohort
  if(length(px_1966[!is.na(px_1966)]) == 38){
    min1 <- which(colnames(px_bc) == "1966")
    
    px_bc <- px_bc[, min1:ncol(px_bc)]
    
    # extract data in a triangle format
    px_triangle <- c()
    for(k in 1:ncol(px_bc)){
      px_triangle <- cbind(px_triangle, c(px_bc[1:(38 - k + 1), k], rep(NA, k - 1)))
    }
    
    colnames(px_triangle) <- colnames(px_bc)
    
  } else {
    min1 <- px_bc[nrow(px_bc), ]
    min1 <- length(min1[!is.na(min1)])
    
    ## select data from the maximum completed birth cohort
    px_triangle <- px_bc[, min1:ncol(px_bc)]
  }
  
  
  ### Prepare period data to create hypothetical data
  
  # create Age:Year matrix contains lx or px
  data_select <- function(Pdata, lf_bc, target){
    
    startY <- as.numeric(as.character(colnames(lf_bc)))[ncol(lf_bc)] + 12
    lastBC <- lf_bc[, ncol(lf_bc)]
    endY   <- startY + length(lastBC[!is.na(lastBC)]) - 1
    
    px_wide <- Pdata %>%
      as.data.frame() %>% 
      filter(x %out% over50) %>%
      mutate(q1x = as.numeric(as.character(q1x)),
             px = 1 - q1x) %>%
      select(Year, x, px) %>%
      filter(Year >= startY & Year <= endY) %>%
      spread(key = Year, value = px) %>%
      select(-x) %>%
      as.matrix()
    
    lx_wide <- matrix(NA, ncol = ncol(px_wide), nrow = nrow(px_wide))
    lx_wide[1, ] <- 1
    for(i in 1:(nrow(lx_wide)-1)){
      lx_wide[i+1, ] <- lx_wide[i, ] * px_wide[i, ]
    }
    
    colnames(lx_wide) <- colnames(px_wide)
    
    if(target == "lx"){
      outcome <- lx_wide
    } else {
      outcome <- px_wide
    }
    
    return(outcome)
  }
  
  ## Data from country A using period fertility table
  A1_matper <- data_select(Pdata = A1_per, lf_bc = px_triangle, target = "px")
  
  # make new data strage: hypthetical cohort
  period2cohort <- function(data){
    
    n <- ncol(data)
    
    bc <- c()
    for(i in 2:n){
      row <- c()
      row <- c(data[i, - c(1:(i - 1))], rep(NA, i-1))
      bc  <- rbind(bc, row)
    }
    bc <- rbind(data[1,], bc)
    
    years <- as.numeric(as.character(colnames(bc)))
    colnames(bc) <- years - 12
    rownames(bc) <- NULL
    
    return(bc)
  }
  
  A1_hypbc <- period2cohort(data = A1_matper)
  
  ## combine cohort data and hypothetical data
  A1_hypbc <- rbind(A1_hypbc, matrix(NA, dim(px_triangle)[1] - dim(A1_hypbc)[1], ncol(A1_hypbc)))
  px_A1_bc <- cbind(px_triangle, A1_hypbc[, -1])
  
  ### Calculate CALC 
  result_CALC <- round(CALCFunc(px_A1_bc), 2)
  
  ### ---- adjusted CALC ----
  
  ### Fertility
  ## Function to create a matrix of qx
  wideadjqx <- function(Fdata){
    outcome <- Fdata %>% 
      #filter(Cohort >= Y1) %>% # select year from the same year
      filter(x %out% over50) %>%              # select age (12- to 49)
      select(Cohort, x, q1x) %>% 
      mutate(q1x = as.numeric(as.character(q1x)),
             q1x = ifelse((x %in% c("12-", "13", "14", "15")) & is.na(q1x), 0, q1x)) %>% 
      na.omit() %>% 
      spread(key = Cohort, value = q1x) %>%
      select(-x)
    
    outcome <- outcome[, which(colnames(outcome) == "1966"):ncol(outcome)]
    
    return(outcome)
  }
  
  Fqx_bc <- wideadjqx(Fdata = A1_bc)
  
  ## hyptohetical cohort data of qx
  Fqx_hypbc = 1 - A1_hypbc
  
  ## combine cohort data and hypothetical data
  Fert_qx   <- cbind(Fqx_bc, Fqx_hypbc[, -1])
  
  ### Mortality
  
  Death_matBC6603 <- Death_matBC[, which(colnames(Death_matBC) == "1966"):which(colnames(Death_matBC) == "2003")]
  Exp_matBC6603   <- Exp_matBC[, which(colnames(Exp_matBC) == "1966"):which(colnames(Exp_matBC) == "2003")]
  Mort_qx <- Death_matBC6603 / unlist(Exp_matBC6603)
  
  # To make clear triangle matrix
  if(length(Mort_qx[, 2]) == 38){
    
    for(r in 2:ncol(Mort_qx)){
      Mort_qx[40 - r, r] <- NA
    }
  } else {
    Mort_qx <- Mort_qx
  }
  
  ### combine both fertility and mortality
  adjqx_bc <- Fert_qx + Mort_qx
  adjpx_bc <- 1 - adjqx_bc
  adjpx_bc <- as.matrix(adjpx_bc)

  result_adjCALC <- round(CALCFunc(adjpx_bc), 2)
  
  
  ### ---- RESULTS ----
  eywc <- c(Names2, result_CALC, result_adjCALC, EYWC, adj_EYWC, EYWC_bc, adj_EYWC_bc)
  
  
  outcome <- rbind(outcome, eywc)
  print(Names)
}

colnames(outcome) <- c("Country", "CALC", "Adjusted CALC", "Period EYWC", "Adjusted EYWC", "Cohort EYWC", "Adjusted cohort EYWC")
rownames(outcome) <- NULL
write.csv(outcome, "Results/Results_perbcadjEYWC_CALC.csv")
