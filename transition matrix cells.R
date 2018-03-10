
#This contains all the transition probabilities. We'll fill these into our transition matrix.

#Format = currentgenotype_nextgenotype

#Row 1
AA_AA <- AA + (1/2)*(AB + AG)
AA_AB <- BB + (1/2)*(AB + BG) + (b/2)*(AG + BG) + b*GG
AA_BB <- 0
AA_AG <- (1-b)*(1-g)*(1/2)*(AG + BG) + (1-b)*(1-g)*GG
AA_BG <- 0
AA_GG <- (1-b)*(g)*(1/2)*(AG + BG) + (1-b)*(g)*GG

#Row 2
AB_AA <- (1/2)*(AA) + (1/4)*(AB + AG)
AB_AB <- (1/2)*(AA + AB + BB) + (1/4)*(AG + BG) + (1/4)*b*(AG + BG) + (1/2)*b*GG
AB_BB <- (1/4)*(AB + BG) + (1/2)*BB + (1/4)*b*(AG + BG) + (1/2)*b*GG
AB_AG <- (1/4)*(1-b)*(1-g)*(AG + BG) + (1/2)*(1-b)*(1-g)*GG
AB_BG <- (1/4)*(1-b)*(AG + BG) + (1/2)*(1-b)*GG
AB_GG <- (1/4)*(1-b)*(g)*(AG + BG) + (1/2)*(1-b)*(g)*GG

#Row 3
BB_AA <- 0
BB_AB <- AA + (1/2)*(AB + AG) + (1/2)*b*AG
BB_BB <- BB + (1/2)*(AB + BG) + (1/2)*b*(AG + BG) + b*GG
BB_AG <- 0
BB_BG <- (1/2)*(1-b)*(AG + BG) + (1-b)*GG
BB_GG <- 0

#Row 4
AG_AA <- (1/2)*AA + (1/4)*(AB + AG)
AG_AB <- (1/2)*b*(AA + AG + GG) + (1/4)*(AB + BG) + (1/4)*b*(AB + BG) + (1/2)*BB
AG_BB <- (1/2)*b*(BG) + (1/4)*(b^2)*(AG + BG) + (1/4)*b*(AB + BG) + (1/2)*(b^2)*GG
AG_AG <- (1/2)*(1-b)*(1-g)*(AA + AG + GG) + (1/4)*(1-b)*(1-g)*(AB + BG)
AG_BG <- (1/4)*(1-b)*(AB + BG) + (1/2)*(1-b)*BB + (1/4)*(2*b*(1-b))*(AG + BG) + (1/2)*(2*b*(1-b))*(GG)
AG_GG <- (1/2)*(1-b)*g*(AA + AG + GG) + (1/4)*(1-b)*g*(AB + BG) + (1/4)*((1-b)^2)*(AG + BG) + (1/2)*((1-b)^2)*(GG)

#Row 5
BG_AA <- 0
BG_AB <- (1/2)*AA + (1/2)*b*AA + (1/4)*(AB + AG) + (1/4)*b*(AB + AG)
BG_BB <- (1/4)*(AB + BG) + (1/4)*b*(AB + AG + BG) + (1/2)*BB + (1/2)*b*(BB + GG) + (1/4)*(b^2)*(AG + BG) + (1/2)*(b^2)*(GG)
BG_AG <- (1/2)*(1-b)*(1-g)*AA + (1/4)*(1-b)*(1-g)*(AB + AG)
BG_BG <- (1/4)*(1-b)*(AB + AG) + (1/2)*(1-b)*(BB + BG + GG) + (1/4)*(2*b*(1-b))*(AG + BG) + (1/2)*(2*b*(1-b))*(GG)
BG_GG <- (1/2)*(1-b)*(g)*AA + (1/4)*(1-b)*(g)*(AB + AG) + (1/4)*((1-b)^2)*(AG + BG) + (1/2)*((1-b)^2)*(GG)

#Row 6                        
GG_AA <- 0
GG_AB <- b*AA + (1/2)*b*(AB + AG)
GG_BB <- b*BB + (1/2)*b*(AB + BG) + (1/2)*(b^2)*(AG + BG) + (b^2)*GG
GG_AG <- (1-b)*(1-g)*AA + (1/2)*(1-b)*(1-g)*(AB + AG)
GG_BG <- (1-b)*BB + (1/2)*(1-b)*(AB + BG) + (1/2)*(2*b*(1-b))*(AG + BG) + (2*b*(1-b))*GG
GG_GG <- (1-b)*g*AA + (1/2)*(1-b)*g*(AB + AG) + (1/2)*((1-b)^2)*(AG + BG) + ((1-b)^2)*(GG)

