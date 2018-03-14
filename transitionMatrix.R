#Initial Genotype Population
#Units = number of mosquitoes
AA <- 99
AB <- 0
BB <- 0
AG <- 0
BG <- 0
GG <- 1

#Important Constants
#k is the carrying capacity
#Units = number of mosquitoes
k <- 100

#Reproduction Rates
#r = intrinsic growth rate
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

#Transition Matrix set-up
#This contains all the transition probabilities. We'll fill these into our transition matrix.
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
rowAA <- c(AA_AA, AA_AB, AA_BB, AA_AG, AA_BG, AA_GG)
rowAB <- c(AB_AA, AB_AB, AB_BB, AB_AG, AB_BG, AB_GG)
rowBB <- c(BB_AA, BB_AB, BB_BB, BB_AG, BB_BG, BB_GG)
rowAG <- c(AG_AA, AG_AB, AG_BB, AG_AG, AG_BG, AG_GG)
rowBG <- c(BG_AA, BG_AB, BG_BB, BG_AG, BG_BG, BG_GG)
rowGG <- c(GG_AA, GG_AB, GG_BB, GG_AG, GG_BG, GG_GG)

#Transition matrix
(trans_matrix <- rbind(rowAA, rowAB, rowBB, rowAG, rowBG, rowGG))

#Gene Drive Efficiency
#from Gantz et al.
g <- 0.98

#Gene Drive Resistance Rate
# chance of becoming resistant to the gene drive
# arbitratily chosen for now
b <- 0.1

