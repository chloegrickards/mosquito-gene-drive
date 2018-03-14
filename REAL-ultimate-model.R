#sorry it's long


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
#
#g is gene drive efficiency
#from Gantz et al.
g <- 0.98
#b is gene drive efficiency rate
# chance of becoming resistant to the gene drive
# arbitratily chosen for now
b <- 0.1

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

#Transition Matrix
# Takes place in two parts
#
# Part 1: Genotype Transition matrix
#This function lists all of the genotype changes from all possible reproduction events
# tracks changes from genotype to genotype, dependent on:
# - Mendelian inheritance
# - resistance to the gene drive
# - gene drive
# - intrinsic growth rate
# - fitness cost to the intrinsic growth rate
#
# Part 2: Population Transition Matrix
# 
genotype_transition_matrix <- function(pAA, pAB, pBB, pAG, pBG, pGG, r, r_xG, r_GG, b, g) {
  
  # First calculate each cell of the matrix. Then assemble the matrix
  #Format = currentgenotype_nextgenotype

  #Row 1
  AA_AA <- r*pAA + (1/2)*(r*pAB + ((r+r_xG)/2)*pAG)
  AA_AB <- r*pBB + (1/2)*(r*pAB + ((r+r_xG)/2)*pBG) + (b/2)*(((r+r_xG)/2)*pAG + ((r+r_xG)/2)*pBG) + b*((r+r_GG)/2)*pGG
  AA_BB <- 0
  AA_AG <- (1-b)*(1-g)*(1/2)*((r+r_xG)/2)*(pAG + pBG) + (1-b)*(1-g)*((r+r_GG)/2)*pGG
  AA_BG <- 0
  AA_GG <- (1-b)*g*(1/2)*((r+r_xG)/2)*(pAG + pBG) + (1-b)*g*((r+r_xG)/2)*pGG

  #Row 2
  AB_AA <- (1/2)*r*pAA + (1/4)*(r*pAB + ((r+r_xG)/2)*pAG)
  AB_AB <- (1/2)*r*(pAA + pAB + pBB) + (1/4)*((r+r_xG)/2)*(pAG + pBG) + (1/4)*b*((r+r_xG)/2)*(pAG + pBG) + (1/2)*b*((r+r_GG)/2)*pGG
  AB_BB <- (1/4)*(r*pAB + ((r+r_xG)/2)*pBG) + (1/2)*r*pBB + (1/4)*b*((r+r_xG)/2)*(pAG + pBG) + (1/2)*b*((r+r_GG)/2)*pGG
  AB_AG <- (1/4)*(1-b)*(1-g)*((r+r_xG)/2)*(pAG + pBG) + (1/2)*(1-b)*(1-g)*((r+r_GG)/2)*pGG
  AB_BG <- (1/4)*(1-b)*((r+r_xG)/2)*(pAG + pBG) + (1/2)*(1-b)*((r+r_GG)/2)*pGG
  AB_GG <- (1/4)*(1-b)*g*((r+r_xG)/2)*(pAG + pBG) + (1/2)*(1-b)*g*((r+r_xG)/2)*pGG

  #Row 3
  BB_AA <- 0
  BB_AB <- r*pAA + (1/2)*(r*pAB + ((r+r_xG)/2)*pAG)
  BB_BB <- r*pBB + (1/2)*(r*pAB + ((r+r_xG)/2)*pBG) + (1/2)*b*((r+r_xG)/2)*(pAG + pBG) + b*((r+r_GG)/2)*pGG
  BB_AG <- 0
  BB_BG <- (1/2)*(1-b)*((r+r_xG)/2)*(pAG + pBG) + (1-b)*((r+r_GG)/2)*pGG
  BB_GG <- 0

  #Row 4
  AG_AA <- (1/2)*r*pAA + (1/4)*(r*pAB + ((r+r_xG)/2)*pAG)
  AG_AB <- (1/2)*b*(r*pAA + ((r+r_xG)/2)*pAG + ((r+r_GG)/2)*pGG) + (1/4)*(r*pAB + ((r+r_xG)/2)*pBG) + (1/4)*b*(r*pAB + ((r+r_xG)/2)*pBG) + (1/2)*r*pBB
  AG_BB <- (1/2)*b*r*pBB + (1/4)*(b^2)*((r+r_xG)/2)*(pAG + pBG) + (1/4)*b*(r*pAB + ((r+r_xG)/2)*pBG) + (1/2)*(b^2)*((r+r_GG)/2)*pGG
  AG_AG <- (1/2)*(1-b)*(1-g)*(r*pAA + ((r+r_xG)/2)*pAG + ((r+r_GG)/2)*pGG) + (1/4)*(1-b)*(1-g)*(r*pAB + ((r+r_xG)/2)*pBG)
  AG_BG <- (1/4)*(1-b)*(r*pAB + ((r+r_xG)/2)*pBG) + (1/2)*(1-b)*r*pBB + (1/4)*(2*b*(1-b))*((r+r_xG)/2)*(pAG + pBG) + (1/2)*(2*b*(1-b))*((r+r_GG)/2)*pGG
  AG_GG <- (1/2)*(1-b)*g*(r*pAA + ((r+r_xG)/2)*pAG + ((r+r_GG)/2)*pGG) + (1/4)*(1-b)*g*(r*pAB + ((r+r_xG)/2)*pBG) + (1/4)*((1-b)^2)*((r+r_xG)/2)*(pAG + pBG) + (1/2)*((1-b)^2)*((r+r_GG)/2)*pGG

  #Row 5
  BG_AA <- 0
  BG_AB <- (1/2)*r*pAA + (1/2)*b*r*pAA + (1/4)*(r*pAB + ((r+r_xG)/2)*pAG) + (1/4)*b*(r*pAB + ((r+r_xG)/2)*pAG)
  BG_BB <- (1/4)*(r*pAB + ((r+r_xG)/2)*pBG) + (1/4)*b*(r*pAB + ((r+r_xG)/2)*pAG) + (1/2)*r*pBB + (1/2)*b*(r*pBB + ((r+r_xG)/2)*pBG + ((r+r_GG)/2)*pGG) + (1/4)*(b^2)*((r+r_xG)/2)*(pAG + pBG) + (1/2)*(b^2)*((r+r_GG)/2)*pGG
  BG_AG <- (1/2)*(1-b)*(1-g)*r*pAA + (1/4)*(1-b)*(1-g)*(r*pAB + ((r+r_xG)/2)*pAG)
  BG_BG <- (1/4)*(1-b)*(r*pAB + ((r+r_xG)/2)*pAG) + (1/2)*(1-b)*(r*pBB + ((r+r_xG)/2)*pBG + ((r+r_GG)/2)*pGG) + (1/4)*(2*b*(1-b))*((r+r_xG)/2)*(pAG + pBG) + (1/2)*(2*b*(1-b))*((r+r_GG)/2)*pGG
  BG_GG <- (1/2)*(1-b)*g*r*pAA + (1/4)*(1-b)*g*(r*pAB + ((r+r_xG)/2)*pAG) + (1/4)*((1-b)^2)*((r+r_xG)/2)*(pAG + pBG) + (1/2)*((1-b)^2)*((r+r_GG)/2)*pGG

  #Row 6                        
  GG_AA <- 0
  GG_AB <- b*r*pAA + (1/2)*b*(r*pAB + ((r+r_xG)/2)*pAG)
  GG_BB <- b*r*pBB + (1/2)*b*(r*pAB + ((r+r_xG)/2)*pBG) + (1/2)*(b^2)*((r+r_xG)/2)*(pAG + pBG) + (b^2)*((r+r_GG)/2)*pGG
  GG_AG <- (1-b)*(1-g)*r*pAA + (1/2)*(1-b)*(1-g)*(r*pAB + ((r+r_xG)/2)*pAG)
  GG_BG <- (1-b)*r*pBB + (1/2)*(1-b)*(r*pAB + ((r+r_xG)/2)*pBG) + (1/2)*(2*b*(1-b))*((r+r_xG)/2)*(pAG + pBG) + (2*b*(1-b))*((r+r_GG)/2)*pGG
  GG_GG <- (1-b)*g*r*pAA + (1/2)*(1-b)*g*(r*pAB + ((r+r_xG)/2)*pAG) + (1/2)*((1-b)^2)*((r+r_xG)/2)*(pAG + pBG) + ((1-b)^2)*((r+r_GG)/2)*pGG

  #Concatenate rows
  gen_rowAA <- c(AA_AA, AA_AB, AA_BB, AA_AG, AA_BG, AA_GG)
  gen_rowAB <- c(AB_AA, AB_AB, AB_BB, AB_AG, AB_BG, AB_GG)
  gen_rowBB <- c(BB_AA, BB_AB, BB_BB, BB_AG, BB_BG, BB_GG)
  gen_rowAG <- c(AG_AA, AG_AB, AG_BB, AG_AG, AG_BG, AG_GG)
  gen_rowBG <- c(BG_AA, BG_AB, BG_BB, BG_AG, BG_BG, BG_GG)
  gen_rowGG <- c(GG_AA, GG_AB, GG_BB, GG_AG, GG_BG, GG_GG)

  # Genotype Transition matrix
  #This function lists all of the genotype changes from all possible reproduction events
  # tracks changes from genotype to genotype, dependent on:
  # - Mendelian inheritance
  # - resistance to the gene drive
  # - gene drive
  # - intrinsic growth rate
  # - fitness cost to the intrinsic growth rate
  gtm <- rbind(gen_rowAA, gen_rowAB, gen_rowBB, gen_rowAG, gen_rowBG, gen_rowGG))
  
  # Logistic Growth Equation
  # Integrated version, so that it fits with the time series format
  # from dN/dt = r*N*(K-N)/K
  # Applied to the gene drive equation since density dependence impact developing (larval) mosquitoes
  # And we need a way to control the gene drive proportions anyway so that the gene drive doesn't explode our model :)
  # Equation below does NOT include r, that's determined by each transition matrix cell
  log_growth <- (1/s_r)*(s_death*k/(s_death + (k - s_death)*math.exp(-s_r)) - s_death)
  
  
  pop_rowAA <- c(AA_AA, AA_AB, AA_BB, AA_AG, AA_BG, AA_GG)
  pop_rowAB <- c(AB_AA, AB_AB, AB_BB, AB_AG, AB_BG, AB_GG)
  pop_rowBB <- c(BB_AA, BB_AB, BB_BB, BB_AG, BB_BG, BB_GG)
  pop_rowAG <- c(AG_AA, AG_AB, AG_BB, AG_AG, AG_BG, AG_GG)
  pop_rowBG <- c(BG_AA, BG_AB, BG_BB, BG_AG, BG_BG, BG_GG)
  pop_rowGG <- c(GG_AA, GG_AB, GG_BB, GG_AG, GG_BG, GG_GG)
  
  return(ptm)
}

population_transition_matrix <- function(args) {
  #AA_AA = 1-(sAA[i-1]*(mu + (1-mu)*rho))/sAA[i-1] + sAA_AA_r/(sAA[i-1]*s_r)
  #AA_AA = sAA[i-1] - death*sAA[i-1] + sAA_AA_r*(growth eqn
  #sAA_AA_r/(sAA[i-1]*s_r)*(s_death*k/(s_death + (k - s_death)*math.exp(-s_r)) - s_death)
  
  rowAA <- c(AA_AA, AA_AB, AA_BB, AA_AG, AA_BG, AA_GG)
  rowAB <- c(AB_AA, AB_AB, AB_BB, AB_AG, AB_BG, AB_GG)
  rowBB <- c(BB_AA, BB_AB, BB_BB, BB_AG, BB_BG, BB_GG)
  rowAG <- c(AG_AA, AG_AB, AG_BB, AG_AG, AG_BG, AG_GG)
  rowBG <- c(BG_AA, BG_AB, BG_BB, BG_AG, BG_BG, BG_GG)
  rowGG <- c(GG_AA, GG_AB, GG_BB, GG_AG, GG_BG, GG_GG)
  
  #Transition matrix
  ptm <- rbind(rowAA, rowAB, rowBB, rowAG, rowBG, rowGG)) * 

  return(ptm)
  
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
  population_prev <- c(AA[i-1], AB[i-1], BB[i-1], AG[i-1], BG[i-1], GG[i-1])
  tot <- sum(population_prev)
  
  #Setting up proportions
  pAA <- AA[i-1]/tot
  pAB <- AB[i-1]/tot
  pBB <- BB[i-1]/tot
  pAG <- AG[i-1]/tot
  pBG <- BG[i-1]/tot
  pGG <- GG[i-1]/tot
  genotypes_prev <- c(pAA, pAB, pBB, pAG, pBG, pGG)
  
  genotypes_curr = genotypes_prev*genotype_transition_matrix(pAA, pAB, pBB, pAG, pBG, pGG, r, r_xG, r_GG, b, g)
  
  population_curr = population_prev*population_transition_matrix(genotypes_curr)

  

}

