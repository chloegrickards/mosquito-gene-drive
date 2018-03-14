
##GENOTYPES##
#Create empty genotype vectors to fill
AA <- rep(0,tspan)
AB <- rep(0,tspan)
BB <- rep(0,tspan)
AG <- rep(0,tspan)
BG <- rep(0,tspan)
GG <- rep(0,tspan)

#Initial population
AA[1] <- 99
AB[1] <- 0
BB[1] <- 0
AG[1] <- 0
BG[1] <- 0
GG[1] <- 1

##GENE DRIVE RATES##
# gamma is the gene drive dominance (relavent to AG and BG)
gamma <- 0.5
# xi is the actual resistance the gene drive gives
xi <- 0.9

##DEATH RATES##
#Natural Death Rate
#life span
# Units = number of generations
lspan <- 3
#death rate
delta <- 1/lspan

#Infection rate
beta <- 0.3

#Infection Death Rate
delta_i <- 0.5

#Death rates for each genotype
d <- delta + beta*delta_i
dxG <- delta + beta*delta_i*(1 - gamma*xi)
dGG <- delta + beta*delta_i*(1-xi)

transition_matrix <- function(args) {
  
  return(tm)
}

# R, first index = 1 because uggghh why
for(i in 2:tspan) {
  
  #Accounting for death rates
  AA[i-1] <- AA[i-1]*(1 - d)
  AB[i-1] <- AB[i-1]*(1 - d)
  BB[i-1] <- BB[i-1]*(1 - d)
  AG[i-1] <- AG[i-1]*(1 - dxG)
  BG[i-1] <- BG[i-1]*(1 - dxG)
  GG[i-1] <- GG[i-1]*(1 - dGG)
  tot <- AA[i-1] + AB[i-1] + BB[i-1] + AG[i-1] + BG[i-1] + GG[i-1]
  
  #Setting up proportions
  pAA <- AA[i-1]/tot
  pAB <- AB[i-1]/tot
  pBB <- BB[i-1]/tot
  pAG <- AG[i-1]/tot
  pBG <- BG[i-1]/tot
  pGG <- GG[i-1]/tot
  
  
  genotypes[i + 1] = transition_matrix*genotypes[i]
  
  

}

