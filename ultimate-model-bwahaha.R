
##GENOTYPES##
#Create empty genotype vectors to fill
AA <- rep(0,tspan)
AB <- rep(0,tspan)
BB <- rep(0,tspan)
AG <- rep(0,tspan)
BG <- rep(0,tspan)
GG <- rep(0,tspan)
#
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
#
# xi is the actual resistance the gene drive gives
xi <- 0.9

##INFECTION RATES##
#beta is infection rate
beta <- 0.3
#
#delta_i is infection Death Rate
delta_i <- 0.5

##POPULATION RATES##
#Death Rates
#lspan is life span
# Units = number of generations
lspan <- 3
#delta is natural death rate
delta <- 1/lspan
#d, dxG, and dGG are death rates for each genotype
# from natural and infection-caused death
d <- delta + beta*delta_i
dxG <- delta + beta*delta_i*(1 - gamma*xi)
dGG <- delta + beta*delta_i*(1-xi)
#
#Growth Rates
# r = intrinsic growth rate
# Use r for genotypes: AA, AB, BB
#Units = number of mosquitoes
r <- 40
#g = fitness cost of the gene drive
#Units = number of mosquitoes
cost_g <- 0.1*r
#r_xG = adjusted growth rate for genotypes: AG, BG
r_xG <- r - cost_g
#r_GG = adjusted growth rate for genotype: GG
r_GG <- r - 2*cost_g
#
#k is the carrying capacity
#Units = number of mosquitoes
k <- 100


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

