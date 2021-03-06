# ===========================================================================================
# Code to calculate CALC and decomposition
# Ryohei Mogi, rmogi@ced.uab.es
# with Jessica Nisen, Vladimir Canudas-Romo
# ===========================================================================================

library(tidyverse)
over50 <- c("50", "51", "52", "53", "54", "55+")
`%out%` = Negate(`%in%`)
bpy.colors <- function (n = 100, cutoff.tails = 0.1, alpha = 1){
  n <- as.integer(n[1])
  if (n <= 0)
    return(character(0))
  
  if (cutoff.tails >= 1 || cutoff.tails < 0)
    stop("cutoff.tails should be in [0, 1]")
  i = seq(0.5 * cutoff.tails, 1 - 0.5 * cutoff.tails, length = n)
  r = ifelse(i < .25, 0, ifelse(i < .57, i / .32 - .78125, 1))
  g = ifelse(i < .42, 0, ifelse(i < .92, 2 * i - .84, 1))
  b = ifelse(i < .25, 4 * i, ifelse(i < .42, 1,
                                    ifelse(i < .92, -2 * i + 1.84, i / .08 - 11.5)))
  rgb(r, g, b, alpha)
}

library(RColorBrewer)

mypalette  <- rev(brewer.pal(8, "YlGnBu"))
mypalette2 <- rev(brewer.pal(8, "YlOrRd"))

WildColors <- c(mypalette[1:4], "white", "white", mypalette2[c(6, 4, 2, 1)])

Mycolor <- c(mypalette[1:4], "white", "white", "#FED976", "#FFEDA0", "#D4B9DA", "#FD8D3C")

levels <- c(-1, -0.1, -0.01, -0.001, -0.0001, 0, .0001, .001, .01, .1, 1)

## legend bar
customAxis <- function() { 
  n <- length(levels) 
  y <- seq(min(levels), max(levels), length.out = n) 
  rect(0, y[1:(n-1)], 1, y[2:n], col = WildColors) 
  axis(4, at = y, labels = levels) 
} 


country  <- c("BLR", "CZE", "DNK", "EST", "HUN", "JPN", "LTU", 
              "NLD", "ESP", "SWE", "USA")
country2 <- c("Belarus", "Czechia", "Denmark", "Estonia", "Hungary", "Japan", "Lithuania", 
              "The Netherlands", "Spain", "Sweden", "The US")


Names = c(country[10], country[11])
Names2 = c(country2[10], country2[11])

#### ---- function ----
# create a diagonal matrix for CALC choosing lx or px
lxpx <- function(Names, lxpx){
  
  # target life table function: lx or px
  target <- lxpx
  
  ## Data from country 1
  A1 <- read.table(paste("Data/", Names, "cft.txt", sep = ""), header = TRUE, fill = TRUE, skip = 2)
  
  ## Function to create a matrix of lx or px
  widelxpx <- function(data){
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
  lf_bc <- widelxpx(data = A1)
  
  if(any(colnames(lf_bc) == "1966")){
    lf_1966 <- lf_bc[, "1966"]
  } else {
    lf_1966 <- NA
  }
  
  ## the position of the maximum completed birth cohort
  if(length(lf_1966[!is.na(lf_1966)]) == 38){
    min1 <- which(colnames(lf_bc) == "1966")
    
    lf_bc <- lf_bc[, min1:ncol(lf_bc)]
    
    # extract data in a triangle format
    lf_triangle <- c()
    for(k in 1:ncol(lf_bc)){
      lf_triangle <- cbind(lf_triangle, c(lf_bc[1:(38 - k + 1), k], rep(NA, k - 1)))
    }
    
    colnames(lf_triangle) <- colnames(lf_bc)
    
  } else {
    min1 <- lf_bc[nrow(lf_bc), ]
    min1 <- length(min1[!is.na(min1)])
    
    ## select data from the maximum completed birth cohort
    lf_triangle <- lf_bc[, min1:ncol(lf_bc)]
  }
  
  
  ### Prepare period data to create hypothetical data
  
  # create Age:Year matrix contains lx or px
  data_select <- function(Pdata, lf_bc){
    
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
  A1_per <- read.table(paste("Data/", Names, "pft.txt", sep = ""), header = TRUE, fill = TRUE, skip = 2)
  A1_per <- data_select(Pdata = A1_per, lf_bc = lf_triangle)
  
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
  
  A1_hypbc <- period2cohort(data = A1_per)
  
  ## combine cohort data and hypothetical data
  A1_hypbc <- rbind(A1_hypbc, matrix(NA, dim(lf_triangle)[1] - dim(A1_hypbc)[1], ncol(A1_hypbc)))
  lf_A1_bc <- cbind(lf_triangle, A1_hypbc[, -1])
  
  
  return(lf_A1_bc)
}

# calculate CALC and decomposition
CALCDecompFunction  <- function(px1, px2, lxLx, Name1, Name2){
  CALClx  <- c()
  CALClx1 <- c()
  CALClx2 <- c()
  PxCh   <- c()
  
  PxCh <- log(px2 / px1)
  PxCh <- ifelse(is.na(PxCh), 0, PxCh)
  colnames(PxCh) <- rownames(PxCh) <- NULL
  # change the order: 1st column (the youngest cohort) -> the last column (the oldest cohort)
  PxCh <- PxCh[, ncol(PxCh):1]
  
  px2CALlx <- function(px){
    # px matrix to lx
    lx <- apply(px, 2, cumprod)
    lx <- rbind(rep(1, ncol(px)), lx)
    
    # lx to CAL lx
    CALlx <- c()
    for(i in 1:ncol(px)){
      order <- 38:1
      CALlx[i] <- lx[order[i], i]
    }
    CALlx <- rev(CALlx)
    
    return(CALlx)
  }
  
  CALClx1 <- px2CALlx(px1)
  CALClx2 <- px2CALlx(px2)
  
  #CALClx1 <- apply(px1, 2, prod, na.rm = T)
  #CALClx2 <- apply(px2, 2, prod, na.rm = T)
  
  CALClx_mid <- t(matrix(rep((CALClx1 + CALClx2)/2, 38), length(CALClx1)))
  
  # calculate CALC
  CALC <- function(lx, type){
    
    # CALC using lx
    CALC_lx <- sum(lx[-1]) + 0.5
    
    # CALC using Lx
    Lx <- (lx[1:37] + lx[2:38]) / 2
    Lx <- c(Lx, lx[38])
    CALC_Lx <- sum(Lx)
    
    out <- ifelse(type == "lx", CALC_lx, CALC_Lx)
    
    return(out)
  }
  
  # final output
  A1 <- CALC(lx = CALClx1, type = lxLx)
  A2 <- CALC(lx = CALClx2, type = lxLx)
  A3 <- sum(CALClx2) - sum(CALClx1)
  A4 <- sum(PxCh * CALClx_mid)
  
  
  print(rbind(c(paste("CALC-", Name1), paste("CALC-", Name2), "Diff", "est-Diff"), round(c(A1, A2, A3, A4), 2)))
  return(PxCh * CALClx_mid)
}

CALCFunc <- function(px, lxLx){
  
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
  
  if(lxLx == "Lx"){
    CALC <- round(sum(CALLx), 2)
  } else {
    CALC <- round(0.5 + sum(CALlx[-1]), 2)
  }
  
  return(CALC)
}

#### ---- prepare lx and px ----
#lx1 <- lxpx(Names = Names[1], lxpx = "lx")
#lx2 <- lxpx(Names = Names[2], lxpx = "lx")

px1 <- lxpx(Names = Names[1], lxpx = "px")
px2 <- lxpx(Names = Names[2], lxpx = "px")


#### ---- Calculate CALC and decomposition ----
## Experiment to check the function
px1_f0 = lx1_f0 = matrix(1, 38, 38)
px2_f0 = lx2_f0 = matrix(1, 38, 38)

check <- CALCDecompFunction(px1_f0, px2_f0, "Lx", Names[1], Names[2])
# -> We get 38 CALC for each country.

CALlxDecompBC <- CALCDecompFunction(px1, px2, "Lx", Names[1], Names[2])

#### ---- Results table ----
## The correct assignment of contributions and the cummulative changes
CALlxD  <- matrix(0, 38, 38)
CALlxDS <- CALlxD

Age <- c(12:49)

BC <- 1966:2003

for (y in 1:38){
  for (x in 1:y){
    CALlxD[x, (38 - y + x)]  <- CALlxDecompBC[x, y]			
    CALlxDS[x, (38 - y + x)] <- sum(CALlxDecompBC[(1:x), y])
  }
}

# +12 = the period of CALC
colnames(px1)[ncol(px1)]

tab_contr <- CALlxDecompBC
# table: contribution of each age and cohort
for(i in 1:ncol(CALlxDecompBC)){
  tab_contr[, 1 + ncol(CALlxDecompBC) - i] <- round(CALlxDecompBC[, i], 4)
}
colnames(tab_contr) <- BC
rownames(tab_contr) <- Age

write.csv(tab_contr, paste("Results/Tab_AgeCohortContr-", Names[1], "-", Names[2], ".csv", sep = ""))

# table: cumulative sum of contribution of each age and cohort
tab_contr_cumsum <- apply(tab_contr, 2, cumsum)

write.csv(tab_contr_cumsum, paste("Results/Tab_SUM_AgeCohortContr-", Names[1], "-", Names[2], ".csv", sep = ""))

# calculate all cohort
all_contr <- rev(tab_contr_cumsum[nrow(tab_contr_cumsum), ])
sum_all_contr <- cumsum(all_contr)

Tab2 <- cbind(12:49, all_contr, sum_all_contr, tab_contr[, 1], tab_contr_cumsum[, 1],
              tab_contr[, 13], tab_contr_cumsum[, 13], tab_contr[, 26], tab_contr_cumsum[, 26])
colnames(Tab2) <- c("Age", "Allbc-age", "Allbc-sum", "1966-age", "1966-sum",
                    "1978-age", "1978-sum", "1991-age", "1991-sum")

write.csv(Tab2, paste("Results/Tab2_", Names[1], "-", Names[2], ".csv", sep = ""))

#### ---- Plot ----
options(scipen = 10)
Nm <- paste("Graph/Fig", Names[2], "-", Names[1], "_bc.pdf", sep = "")
pdf(Nm)
par(cex.axis = 1)
par(oma = c(1, 0, 0, 0)) #bottom, right, top, left
filled.contour(BC, Age, t(CALlxDS), levels = levels, 
               col = WildColors, key.axes = customAxis(), 
               ylab = "Cumulative age- & cohort-contribution", xlab = "", cex.lab = 1.1,
               plot.axes = {axis(1, at = c(1966, seq(1970, 2000, by = 5), 2003), 
                                 labels = c("", "1970\n(1982)", "1975\n(1987)", "1980\n(1992)", 
                                            "1985\n(1997)", "1990\n(2002)", "1995\n(2007)", 
                                            "2000\n(2012)", "2003\n(2015)"), hadj = 0.6, padj = 0.5, cex.axis = 0.9)
                 axis(2, at = seq(15, 50, by = 5), labels = seq(15, 50, by = 5))})
mtext("Birth cohort\n(Year)", 1, line = 4.5, adj = 0.4) # adj: (-)left-right(+), line: (+)down-up(-)
mtext(Names2[2], 3, 0.5, adj = 0.9, cex = 0.9)
mtext(Names2[1], 1, 0.5, adj = 0.9, cex = 0.9)
dev.off()

# specific contribution
options(scipen = 10)
Nm <- paste("Graph/Fig", Names[2], "-", Names[1], "_bc_specific.pdf", sep = "")
pdf(Nm)
par(cex.axis = 1)
par(oma = c(1, 0, 0, 0))
filled.contour(BC, Age, t(CALlxD), levels = levels, 
               col = WildColors, key.axes = customAxis(), 
               ylab = "Age- & cohort-contribution", xlab = "", cex.lab = 1.1,
               plot.axes = {axis(1, at = c(1966, seq(1970, 2000, by = 5), 2003), 
                                 labels = c("", "1970\n(1982)", "1975\n(1987)", "1980\n(1992)", 
                                            "1985\n(1997)", "1990\n(2002)", "1995\n(2007)", 
                                            "2000\n(2012)", "2003\n(2015)"), hadj = 0.6, padj = 0.5, cex.axis = 0.9)
                 axis(2, at = seq(15, 50, by = 5), labels = seq(15, 50, by = 5))})
mtext("Birth cohort\n(Year)", 1, line = 4.5, adj = 0.4) # adj: (-)left-right(+), line: (+)down-up(-)
mtext(Names2[2], 3, 0.5, adj = 0.9, cex = 0.9)
mtext(Names2[1], 1, 0.5, adj = 0.9, cex = 0.9)
dev.off()


#### ---- Calculate CALC for all countries ----
CALC <- c()
method <- "Lx"
for(i in 1:length(country)){
  Names  <- country[i]
  Names2 <- country2[i]
  
  px <- lxpx(Names = Names, lxpx = "px")
  calc <- CALCFunc(px, method)
  
  calc_mat <- c(Names2, calc)
  
  CALC <- rbind(CALC, calc_mat)
}
colnames(CALC) <- c("Coutnry", "CALC")
rownames(CALC) <- c()
write.csv(CALC, paste0("Results/Tab_CALC_", method, ".csv"))
