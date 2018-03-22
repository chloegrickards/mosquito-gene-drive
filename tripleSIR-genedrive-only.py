# sorry it's long
import numpy as np
import math
import scipy
import pylab as pl

#To Do:
# - add malaria equations for mosquito and human
# - if time, make prettier by putting the transition matrix in another function or scrupt or something

##GENERAL CONSTANTS##
# Number of mosquito generations
ngen = 100

##GENOTYPES##
# Create empty genotype vectors to fill
AA = np.zeros(ngen)
AB = np.zeros(ngen)
BB = np.zeros(ngen)
AG = np.zeros(ngen)
BG = np.zeros(ngen)
GG = np.zeros(ngen)
#
# Initial population
AA[0] = 99
AB[0] = 0.000000000001
BB[0] = 0.000000000001
AG[0] = 0.000000000001
BG[0] = 0.000000000001
GG[0] = 1

##GENE DRIVE RATES##
# gamma is the gene drive dominance (relavent to AG and BG)
gamma = 0.5
#
# xi is a measure of how much infection resistance that a gene drive can pass on to a mosquito
xi = 0.9
#
# g is gene drive efficiency
# from Gantz et al.
g = 0.98
#
# b is gene drive resistance rate
# chance of becoming resistant to the gene drive
# arbitratily chosen for now
b = 0.1

##INFECTION RATES##
# beta is infection rate
beta = 0.3
#
# delta_i is infection Death Rate
delta_i = 0.5

##POPULATION RATES##
# Death Rates
# lspan is life span
# Units = number of generations
lspan = 3
# delta is natural death rate
delta = 1 / lspan
# d, dxG, and dGG are death rates for each genotype
# from natural and infection-caused death
d = delta + beta * delta_i
dxG = delta + beta * delta_i * (1 - gamma * xi)
dGG = delta + beta * delta_i * (1 - xi)
#
# Growth Rates
# r = intrinsic growth rate
# Use r for genotypes: AA, AB, BB
# Units = number of mosquitoes
r = 40
# g = fitness cost of the gene drive
# Units = number of mosquitoes
cost_g = 0.1 * r
# r_xG = adjusted growth rate for genotypes: AG, BG
r_xG = r - cost_g
# r_GG = adjusted growth rate for genotype: GG
r_GG = r - 2 * cost_g
#
# k is the carrying capacity
# Units = number of mosquitoes
k = 100

##GENE DRIVE SEIR##
Sg = np.zeros(ngen)
Eg = np.zeros(ngen)
Ig = np.zeros(ngen)
Rg = np.zeros(ngen)

##MOSQUITO MALARIA##

##HUMAN MALARIA##


# Iterate through each generation
for i in range(1,ngen):

    # Previous population
    AA_prev = AA[i - 1]
    AB_prev = AB[i - 1]
    BB_prev = BB[i - 1]
    AG_prev = AG[i - 1]
    BG_prev = BG[i - 1]
    GG_prev = GG[i - 1]
    population_prev = [AA_prev, AB_prev, BB_prev, AG_prev, BG_prev, GG_prev]


    # STEP 1: Calculate genotype changes for each new generation

    # Reproducing population
    # Accounting for natural and infection-induced death rates
    AA_repr = AA_prev * (1 - d)
    AB_repr = AB_prev * (1 - d)
    BB_repr = BB_prev * (1 - d)
    AG_repr = AG_prev * (1 - dxG)
    BG_repr = BG_prev * (1 - dxG)
    GG_repr = GG_prev * (1 - dGG)
    population_repr = [AA_repr, AB_repr, BB_repr, AG_repr, BG_repr, GG_repr]
    pop_tot_repr = np.sum(population_repr)


    # Previous genotypes
    # Setting up genotype frequencies from the reproducing population
    pAA = AA_repr / pop_tot_repr
    pAB = AB_repr / pop_tot_repr
    pBB = BB_repr / pop_tot_repr
    pAG = AG_repr / pop_tot_repr
    pBG = BG_repr / pop_tot_repr
    pGG = GG_repr / pop_tot_repr
    genotypes_prev = [pAA, pAB, pBB, pAG, pBG, pGG]

    # Transition Matrix calculations
    # So fun.

    # Row 1
    AA_AA = r * pAA + (1 / 2) * (r * pAB + ((r + r_xG) / 2) * pAG)
    AA_AB = r * pBB + (1 / 2) * (r * pAB + ((r + r_xG) / 2) * pBG) + (b / 2) * (((r + r_xG) / 2) * pAG + ((r + r_xG) / 2) * pBG) + b * ((r + r_GG) / 2) * pGG
    AA_BB = 0
    AA_AG = (1 - b) * (1 - g) * (1 / 2) * ((r + r_xG) / 2) * (pAG + pBG) + (1 - b) * (1 - g) * ((r + r_GG) / 2) * pGG
    AA_BG = 0
    AA_GG = (1 - b) * g * (1 / 2) * ((r + r_xG) / 2) * (pAG + pBG) + (1 - b) * g * ((r + r_xG) / 2) * pGG

    # Row 2
    AB_AA = (1 / 2) * r * pAA + (1 / 4) * (r * pAB + ((r + r_xG) / 2) * pAG)
    AB_AB = (1 / 2) * r * (pAA + pAB + pBB) + (1 / 4) * ((r + r_xG) / 2) * (pAG + pBG) + (1 / 4) * b * ((r + r_xG) / 2) * (pAG + pBG) + (1 / 2) * b * ((r + r_GG) / 2) * pGG
    AB_BB = (1 / 4) * (r * pAB + ((r + r_xG) / 2) * pBG) + (1 / 2) * r * pBB + (1 / 4) * b * ((r + r_xG) / 2) * (pAG + pBG) + (1 / 2) * b * ((r + r_GG) / 2) * pGG
    AB_AG = (1 / 4) * (1 - b) * (1 - g) * ((r + r_xG) / 2) * (pAG + pBG) + (1 / 2) * (1 - b) * (1 - g) * ((r + r_GG) / 2) * pGG
    AB_BG = (1 / 4) * (1 - b) * ((r + r_xG) / 2) * (pAG + pBG) + (1 / 2) * (1 - b) * ((r + r_GG) / 2) * pGG
    AB_GG = (1 / 4) * (1 - b) * g * ((r + r_xG) / 2) * (pAG + pBG) + (1 / 2) * (1 - b) * g * ((r + r_xG) / 2) * pGG

    # Row 3
    BB_AA = 0
    BB_AB = r * pAA + (1 / 2) * (r * pAB + ((r + r_xG) / 2) * pAG)
    BB_BB = r * pBB + (1 / 2) * (r * pAB + ((r + r_xG) / 2) * pBG) + (1 / 2) * b * ((r + r_xG) / 2) * (pAG + pBG) + b * ((r + r_GG) / 2) * pGG
    BB_AG = 0
    BB_BG = (1 / 2) * (1 - b) * ((r + r_xG) / 2) * (pAG + pBG) + (1 - b) * ((r + r_GG) / 2) * pGG
    BB_GG = 0

    # Row 4
    AG_AA = (1 / 2) * r * pAA + (1 / 4) * (r * pAB + ((r + r_xG) / 2) * pAG)
    AG_AB = (1 / 2) * b * (r * pAA + ((r + r_xG) / 2) * pAG + ((r + r_GG) / 2) * pGG) + (1 / 4) * (r * pAB + ((r + r_xG) / 2) * pBG) + (1 / 4) * b * (r * pAB + ((r + r_xG) / 2) * pBG) + (1 / 2) * r * pBB
    AG_BB = (1 / 2) * b * r * pBB + (1 / 4) * (b ** 2) * ((r + r_xG) / 2) * (pAG + pBG) + (1 / 4) * b * (r * pAB + ((r + r_xG) / 2) * pBG) + (1 / 2) * (b ** 2) * ((r + r_GG) / 2) * pGG
    AG_AG = (1 / 2) * (1 - b) * (1 - g) * (r * pAA + ((r + r_xG) / 2) * pAG + ((r + r_GG) / 2) * pGG) + (1 / 4) * (1 - b) * (1 - g) * (r * pAB + ((r + r_xG) / 2) * pBG)
    AG_BG = (1 / 4) * (1 - b) * (r * pAB + ((r + r_xG) / 2) * pBG) + (1 / 2) * (1 - b) * r * pBB + (1 / 4) * (2 * b * (1 - b)) * ((r + r_xG) / 2) * (pAG + pBG) + (1 / 2) * (2 * b * (1 - b)) * ((r + r_GG) / 2) * pGG
    AG_GG = (1 / 2) * (1 - b) * g * (r * pAA + ((r + r_xG) / 2) * pAG + ((r + r_GG) / 2) * pGG) + (1 / 4) * (1 - b) * g * (r * pAB + ((r + r_xG) / 2) * pBG) + (1 / 4) * ((1 - b) ** 2) * ((r + r_xG) / 2) * (pAG + pBG) + (1 / 2) * ((1 - b) ** 2) * ((r + r_GG) / 2) * pGG

    # Row 5
    BG_AA = 0
    BG_AB = (1 / 2) * r * pAA + (1 / 2) * b * r * pAA + (1 / 4) * (r * pAB + ((r + r_xG) / 2) * pAG) + (1 / 4) * b * (r * pAB + ((r + r_xG) / 2) * pAG)
    BG_BB = (1 / 4) * (r * pAB + ((r + r_xG) / 2) * pBG) + (1 / 4) * b * (r * pAB + ((r + r_xG) / 2) * pAG) + (1 / 2) * r * pBB + (1 / 2) * b * (r * pBB + ((r + r_xG) / 2) * pBG + ((r + r_GG) / 2) * pGG) + (1 / 4) * (b ** 2) * ((r + r_xG) / 2) * (pAG + pBG) + ( 1 / 2) * (b ** 2) * ((r + r_GG) / 2) * pGG
    BG_AG = (1 / 2) * (1 - b) * (1 - g) * r * pAA + (1 / 4) * (1 - b) * (1 - g) * (r * pAB + ((r + r_xG) / 2) * pAG)
    BG_BG = (1 / 4) * (1 - b) * (r * pAB + ((r + r_xG) / 2) * pAG) + (1 / 2) * (1 - b) * (r * pBB + ((r + r_xG) / 2) * pBG + ((r + r_GG) / 2) * pGG) + (1 / 4) * (2 * b * (1 - b)) * ((r + r_xG) / 2) * (pAG + pBG) + (1 / 2) * (2 * b * (1 - b)) * ((r + r_GG) / 2) * pGG
    BG_GG = (1 / 2) * (1 - b) * g * r * pAA + (1 / 4) * (1 - b) * g * (r * pAB + ((r + r_xG) / 2) * pAG) + (1 / 4) * ((1 - b) ** 2) * ((r + r_xG) / 2) * (pAG + pBG) + (1 / 2) * ((1 - b) ** 2) * ((r + r_GG) / 2) * pGG

    # Row 6
    GG_AA = 0
    GG_AB = b * r * pAA + (1 / 2) * b * (r * pAB + ((r + r_xG) / 2) * pAG)
    GG_BB = b * r * pBB + (1 / 2) * b * (r * pAB + ((r + r_xG) / 2) * pBG) + (1 / 2) * (b ** 2) * ((r + r_xG) / 2) * (pAG + pBG) + (b ** 2) * ((r + r_GG) / 2) * pGG
    GG_AG = (1 - b) * (1 - g) * r * pAA + (1 / 2) * (1 - b) * (1 - g) * (r * pAB + ((r + r_xG) / 2) * pAG)
    GG_BG = (1 - b) * r * pBB + (1 / 2) * (1 - b) * (r * pAB + ((r + r_xG) / 2) * pBG) + (1 / 2) * (2 * b * (1 - b)) * ((r + r_xG) / 2) * (pAG + pBG) + (2 * b * (1 - b)) * ((r + r_GG) / 2) * pGG
    GG_GG = (1 - b) * g * r * pAA + (1 / 2) * (1 - b) * g * (r * pAB + ((r + r_xG) / 2) * pAG) + (1 / 2) * ((1 - b) ** 2) * ((r + r_xG) / 2) * (pAG + pBG) + ((1 - b) ** 2) * ((r + r_GG) / 2) * pGG

    # Concatenate rows
    gen_rowAA = [AA_AA, AA_AB, AA_BB, AA_AG, AA_BG, AA_GG]
    gen_rowAB = [AB_AA, AB_AB, AB_BB, AB_AG, AB_BG, AB_GG]
    gen_rowBB = [BB_AA, BB_AB, BB_BB, BB_AG, BB_BG, BB_GG]
    gen_rowAG = [AG_AA, AG_AB, AG_BB, AG_AG, AG_BG, AG_GG]
    gen_rowBG = [BG_AA, BG_AB, BG_BB, BG_AG, BG_BG, BG_GG]
    gen_rowGG = [GG_AA, GG_AB, GG_BB, GG_AG, GG_BG, GG_GG]

    # Genotype Transition matrix
    # This function lists all of the genotype changes from all possible reproduction events
    # Output = genotype frequencies
    # tracks changes from genotype to genotype, dependent on:
    # - Mendelian inheritance
    # - resistance to the gene drive
    # - gene drive
    # - intrinsic growth rate
    # - fitness cost to the intrinsic growth rate
    #gtm = [[gen_rowAA], [gen_rowAB], [gen_rowBB], [gen_rowAG], [gen_rowBG], [gen_rowGG]]
    gtm = np.matrix([gen_rowAA,
                     gen_rowAB,
                     gen_rowBB,
                     gen_rowAG,
                     gen_rowBG,
                     gen_rowGG])
    genotypes_curr = genotypes_prev * gtm
    gen_tot_curr = np.sum(genotypes_curr)


    ##STEP 2: Calculate Population Changes for each new generation

    # Logistic Growth Equation
    # Integrated version, so that it fits with the time series format
    # from dN/dt = r*N*(K-N)/K
    # Applied to the gene drive equation since density dependence impact developing (larval) mosquitoes
    # And we need a way to control the gene drive proportions anyway so that the gene drive doesn't explode our model :)
    # Equation below does NOT include r, that's determined by each transition matrix cell
    log_growth = (1 / gen_tot_curr) * (pop_tot_repr * k / (pop_tot_repr + (k - pop_tot_repr) * math.exp(-gen_tot_curr)) - pop_tot_repr)


    #pop_rowAA = [AA[i - 1] / (pAA * log_growth) + AA_AA, AA_AB, AA_BB, AA_AG, AA_BG, AA_GG] * pAA * log_growth
    #pop_rowAB = [AB_AA, AB[i - 1] / (pAB * log_growth) + AB_AB, AB_BB, AB_AG, AB_BG, AB_GG] * pAB * log_growth
    #pop_rowBB = [BB_AA, BB_AB, BB[i - 1] / (pBB * log_growth) + BB_BB, BB_AG, BB_BG, BB_GG] * pBB * log_growth
    #pop_rowAG = [AG_AA, AG_AB, AG_BB, AG[i - 1] / (pAG * log_growth) + AG_AG, AG_BG, AG_GG] * pAG * log_growth
    #pop_rowBG = [BG_AA, BG_AB, BG_BB, BG_AG, BG[i - 1] / (pBG * log_growth) + BG_BG, BG_GG] * pBG * log_growth
    #pop_rowGG = [GG_AA, GG_AB, GG_BB, GG_AG, GG_BG, GG[i - 1] / (pGG * log_growth) + GG_GG] * pGG * log_growth

    # Population Transition matrix
    # This takes the genotype ratios from Step 1 and applies population dynamics.
    # Output = individuals of a certain genotype
    # Dependent on the following:
    # - logistic growth and related growth rates
    # - death rates
    #ptm = [[pop_rowAA], [pop_rowAB], [pop_rowBB], [pop_rowAG], [pop_rowBG], [pop_rowGG]]


    pop_rowAA = [AA_repr / AA_prev + AA_AA * pAA * log_growth / AA_prev, AA_AB * pAA * log_growth / AA_prev, AA_BB * pAA * log_growth / AA_prev, AA_AG * pAA * log_growth / AA_prev, AA_BG * pAA * log_growth / AA_prev, AA_GG * pAA * log_growth / AA_prev]
    pop_rowAB = [AB_AA * pAB * log_growth / AB_prev, AB_repr / AB_prev + AB_AB * pAB * log_growth / AB_prev, AB_BB * pAB * log_growth / AB_prev, AB_AG * pAB * log_growth / AB_prev, AB_BG * pAB * log_growth / AB_prev, AB_GG * pAB * log_growth / AB_prev]
    pop_rowBB = [BB_AA * pBB * log_growth / BB_prev, BB_AB * pBB * log_growth / BB_prev, BB_repr / BB_prev + BB_BB * pBB * log_growth / BB_prev, BB_AG * pBB * log_growth / BB_prev, BB_BG * pBB * log_growth / BB_prev, BB_GG * pBB * log_growth / BB_prev]
    pop_rowAG = [AG_AA * pAG * log_growth / AG_prev, AG_AB * pAG * log_growth / AG_prev, AG_BB * pAG * log_growth / AG_prev, AG_repr / AG_prev + AG_AG * pAG * log_growth / AG_prev, AG_BG * pAG * log_growth / AG_prev, AG_GG * pAG * log_growth / AG_prev]
    pop_rowBG = [BG_AA * pBG * log_growth / BG_prev, BG_AB * pBG * log_growth / BG_prev, BG_BB * pBG * log_growth / BG_prev, BG_AG * pBG * log_growth / BG_prev, BG_repr / BG_prev + BG_BG * pBG * log_growth / BG_prev, BG_GG * pBG * log_growth / BG_prev]
    pop_rowGG = [GG_AA * pGG * log_growth / GG_prev, GG_AB * pGG * log_growth / GG_prev, GG_BB * pGG * log_growth / GG_prev, GG_AG * pGG * log_growth / GG_prev, GG_BG * pGG * log_growth / GG_prev, GG_repr / GG_prev + GG_GG * pGG * log_growth / GG_prev]


    ptm = np.matrix([pop_rowAA,
                     pop_rowAB,
                     pop_rowBB,
                     pop_rowAG,
                     pop_rowBG,
                     pop_rowGG])

    #print(ptm)

    population_curr = population_prev * ptm

   # print(population_curr)


    AA[i] = population_curr[0, 0]
    AB[i] = population_curr[0, 1]
    BB[i] = population_curr[0, 2]
    AG[i] = population_curr[0, 3]
    BG[i] = population_curr[0, 4]
    GG[i] = population_curr[0, 5]

    # STEP 3: Gene Drive SEIR model

    # Susceptible to a gene drive "infection"
    Sg[i] = AA[i] + AB[i]
    # Exposed to a gene drive "infection"
    Eg[i] = AG[i]
    # "Infected" with a gene drive
    Ig[i] = GG[i]
    # Recovered from a gene drive "infection"
    Rg[i] = BB[i] + BG[i]

    # STEP 4: Malaria (mosquito) model

    # Number of susceptible mosquitoes in a population
    # may not need Sm[i], esp. if we're just using X, which would require X[i]
    Sm = AA[i] + AB[i] + BB[i] + (1 - gamma * xi) * AG[i] + (1 - gamma * xi) * BG[i] + (1 - xi) * GG[i]

    # STEP 5: Malaria (human) model



#PLOT


tspan = np.arange(ngen)

#transpose vectors for plotting and store as S1, S2, ..., S6
S1 = scipy.transpose(AA)
S2 = scipy.transpose(AB)
S3 = scipy.transpose(BB)
S4 = scipy.transpose(AG)
S5 = scipy.transpose(BG)
S6 = scipy.transpose(GG)

#plot transposed vectors with labels and corresponding colored lines
pl.figure(1)
pl.plot(tspan, S1, 'r-', label='AA')
pl.plot(tspan, S2, 'b-', label='AB')
pl.plot(tspan, S3, 'g-', label='BB')
pl.plot(tspan, S4, 'o-', label='AG')
pl.plot(tspan, S5, 'y-', label='BG')
pl.plot(tspan, S6, 'd-', label='GG')

#add grid, legend, and labels
pl.grid()
pl.legend(loc='best')
pl.xlabel('Generation')
pl.ylabel('Population %')
pl.title('Gene Drive Evolution')
pl.xticks(scipy.linspace(0, ngen, 11))  # sets ticks for x-axis

pl.figure(2)

S = scipy.transpose(Sg)
E = scipy.transpose(Eg)
I = scipy.transpose(Ig)
R = scipy.transpose(Rg)

pl.plot(tspan, S, 'r-', label='Sg')
pl.plot(tspan, E, 'b-', label='Eg')
pl.plot(tspan, I, 'g-', label='Ig')
pl.plot(tspan, R, 'o-', label='Rg')

pl.legend(loc='best')
pl.xlabel('Generation')
pl.ylabel('Population %')
pl.title('SEIR model')
pl.xticks(scipy.linspace(0, ngen, 11))  # sets ticks for x-axis


pl.show()
