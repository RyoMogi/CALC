# ===========================================================================================
# Code to calculate both period and cohort EYWC using HFD.
# lx is estimated by px and CALC is estimated by the sum of Lx.
# Ryohei Mogi, rmogi@ced.uab.es
# with Jessica Nisen, Vladimir Canudas-Romo
# ===========================================================================================

library(tidyverse)

### ---- Country selecter ----
country  <- c("BLR", "CZE", "DNK", "EST", "HUN", "JPN", "LTU", 
              "NLD", "ESP", "SWE", "USA")
country2 <- c("Belarus", "Czechia", "Denmark", "Estonia", "Hungary", "Japan", "Lithuania", 
              "The Netherlands", "Spain", "Sweden", "The US")

EYWC <- c()
lxLx <- "Lx"

for(i in 1:length(country)){
  Names  <- country[i]
  Names2 <- country2[i]
  
  # Period data
  A1_per <- read.table(paste("Data/", Names, "pft.txt", sep = ""), header = TRUE, fill = TRUE, skip = 2)
  
  # Cohort data
  A1_bc <- read.table(paste("Data/", Names, "cft.txt", sep = ""), header = TRUE, fill = TRUE, skip = 2)
  
  ### ---- Period EYWC ----
  over50 <- c("50", "51", "52", "53", "54", "55+")
  `%out%` = Negate(`%in%`)
  
  # select year (ideally 2015, if they don't have enough data, put the maximum year)
  length_check <- A1_per %>%
    filter(x %out% over50) %>%
    group_by(Year) %>%
    mutate(q1x = as.numeric(as.character(q1x)),
           px = 1 - q1x) %>%
    drop_na(px) %>%
    summarise(length_px = length(px)) %>%
    filter(length_px == 38)
  
  A1_maxY <- max(length_check$Year)
  
  if(A1_maxY >= 2015){
    A1_maxY <- 2015
  } else {
    A1_maxY <- A1_maxY
  }
  
  # calculate period EYWC 
  A1_per <- A1_per %>% 
    filter(Year == A1_maxY) %>%
    filter(x %out% over50) %>%
    mutate(q1x = as.numeric(as.character(q1x)),
           px = 1 - q1x,
           lx = cumprod(px),
           L0x = as.numeric(as.character(L0x)),
           l0x = as.numeric(as.character(l0x)))
  
  lx_per <- c(1, A1_per$lx[-nrow(A1_per)])
  
  Lx_per <- c()
  for(i in 1:37){
    Lx_per[i] <- 0.5 * (lx_per[i] + lx_per[i+1])
  }
  Lx_per[38] <- lx_per[38]
  
  A1_EYWC_per <- round(sum(Lx_per), 2)
  
  # check
  #sum((A1_per$l0x / 10000) -  lx_per)
  #
  #EYWC_Lx <- round(sum(A1_per$L0x / 10000), 2)
  A1_EYWC_per_lx <- round(0.5 + sum(lx_per[-1]), 2)
  
  
  ### ---- Cohort EYWC ----
  # select cohort (ideally 1966 (2015 - 49), if they don't have enough data, put the maximum cohort)
  A1_maxBC <- A1_maxY - 49
  
  # check whether A1_maxBC is within the A1_bc data range
  bc <- unique(A1_bc$Cohort)
  
  if(A1_maxBC %in% bc) {
    
    # calculate period EYWC 
    A1_bc <- A1_bc %>% 
      filter(Cohort == A1_maxBC) %>%
      filter(x %out% over50) %>%
      mutate(q1x = as.numeric(as.character(q1x)),
             px = 1 - q1x,
             lx = cumprod(px))
    
  } else {
    
    # select birth cohort having full age-range (38)
    length_check_bc <- A1_bc %>%
      filter(x %out% over50) %>%
      group_by(Cohort) %>%
      mutate(q1x = as.numeric(as.character(q1x)),
             px = 1 - q1x) %>%
      drop_na(px) %>%
      summarise(length_px = length(px)) %>%
      filter(length_px == 38)
    
    A1_maxBC <- max(length_check_bc$Cohort)
    
    # calculate period EYWC 
    A1_bc <- A1_bc %>% 
      filter(Cohort == A1_maxBC) %>%
      filter(x %out% over50) %>%
      mutate(q1x = as.numeric(as.character(q1x)),
             px = 1 - q1x,
             lx = cumprod(px))
  }
  
  lx_bc <- c(1, A1_bc$lx[-nrow(A1_bc)])
  
  Lx_bc <- c()
  for(i in 1:37){
    Lx_bc[i] <- 0.5 * (lx_bc[i] + lx_bc[i+1])
  }
  Lx_bc[38] <- lx_bc[38]
  
  A1_EYWC_bc <- round(sum(Lx_bc), 2)
  A1_EYWC_bc_lx <- round(0.5 + sum(lx_bc[-1]), 2)
  
  # RESULTS
  if(lxLx == "Lx"){
    eywc <- c(A1_EYWC_per, A1_EYWC_bc, A1_maxY, A1_maxBC)
  } else {
    eywc <- c(A1_EYWC_per_lx, A1_EYWC_bc_lx, A1_maxY, A1_maxBC)
  }
  
  print(eywc)
  
  EYWC <- rbind(EYWC, eywc)
}

EYWC <- cbind(country2, EYWC)
colnames(EYWC) <- c("Country", "Period EYWC", "Cohort EYWC", "Year", "Cohort")
write.csv(EYWC, paste0("Results/Results_EYWC_perbc_Lx.csv"))
