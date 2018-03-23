# Import statements
import numpy as np
import math
import scipy
import pylab as pl
from scipy.integrate import odeint

##GENERAL CONSTANTS##
# Number of mosquito generations
ngen = 25

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

##GENE DRIVE CONSTANTS##
# gamma is the gene drive dominance (relevant to AG and BG)
# No literature value yet
gamma = 0.5
#
# xi is a measure of how much infection resistance that a gene drive can pass on to a mosquito
# Taken from the anti-effector gene used by Gantz et al. 2015
# Exact number found in Isaacs et al. 2011
xi = 0.84
#
# g is gene drive efficiency
# from Gantz et al.
g = 0.98
#
# b is gene drive resistance rate
# chance of becoming resistant to the gene drive
# arbitratily chosen for now
# Gantz et al.
# b = 0.05 in our modification
b = 0.25

##POPULATION CONSTANTS##
# lspan is life span
# Units = number of generations
# From CDC
lspan = 2
#
# delta is natural death rate
delta = 1 / lspan
#
# delta_i is infection Death Rate
# Cator et al.
delta_i = 0.25
#
# r = intrinsic growth rate
# r = how many mosquitoes are born from a mating event
# adjusted growth rate for genotypes: AA, AB, BB
# Units = number of mosquitoes
# Rossignol et al.
r = 89
#
# cost_g = fitness cost of the gene drive
# cost_g = how many less mosquitoes are born from a mating event
# Units = number of mosquitoes
# Unverified, somewhat supported by Vella et al.
cost_g = 0.3 * r
#
# r_xG = adjusted growth rate for genotypes: AG, BG
r_xG = r - cost_g
#
# r_GG = adjusted growth rate for genotype: GG
r_GG = r - 2 * cost_g
#
# k is the carrying capacity
# Units = number of mosquitoes
k = 200

##GENE DRIVE SEIR##
Sg = np.zeros(ngen)
Eg = np.zeros(ngen)
Ig = np.zeros(ngen)
Rg = np.zeros(ngen)

##MALARIA CONSTANTS
#solution arrays: XY
Xh = np.zeros(ngen)
Xh[0] = 0.1
Ym = np.zeros(ngen)
Ym[0] = 0.2

#solution arrays: mosquito malaria
Nm0 = AA[0] + AB[0] + BB[0] + (1-xi*gamma)*AG[0] + (1-xi*gamma)*BG[0] + (1-xi)*GG[0]
Sm = np.zeros(ngen)
Sm[0] = (1-Ym[0])*Nm0
Im = np.zeros(ngen)
Im[0] = Ym[0]*Nm0

#solution arrays: mosquito malaria
#number of humans
Nh = 500
Sh = np.zeros(ngen)
Sh[0] = (1-Xh[0])*Nh
Ih = np.zeros(ngen)
Ih[0] = Xh[0]*Nh

## MALARIA CONSTANTS ##
#Generation time in days
gen = 30
# Expected number of bites on humans per mosquito
# a = 0.01 in our modification
a = 0.15
# Probability of infection of an uninfected mosquito by biting an infectious human
c = 0.82
# Probability of surviving one day
p = 0.81
# Force of mortality (i.e. instantaneous mortality rate) for mosquitoes
f = -np.log(p)
# Rate of mosquito emergence
ma = 2.06
# Equilibrium mosquito density per human
m = ma/f
# Transmission efficiency from infectious mosquito to susceptible human
# Churcher 2017
tr = 0.092
# Extrinsic cycle duration in days
# Brasil et al 2011
n = 14
#Duration of illness
#takes 3-7 days to recover
r = 1/7

def XYmodel(y, t):
    Xh, Ym = y
    dydt =[((m*(a**2)*tr*c*np.exp(-f*n))/(f + a*c*Xh))*(Xh * (1-Xh)) - r*Xh, (a * c * Xh) / (f + a * c * Xh)]
    return dydt



# Iterate through each generation
for i in range(1,ngen):

    #Constant calculations
    # d, dxG, and dGG are death rates for each genotype
    # from natural and infection-caused death
    d = delta + Ym[i-1] * delta_i
    dxG = delta + Ym[i-1] * delta_i * (1 - gamma * xi)
    dGG = delta + Ym[i-1] * delta_i * (1 - xi)

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
    total_pop_curr = np.sum(population_curr)

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

    #Number of mosquitoes in susceptible + infected population
    Nm = AA[i] + AB[i] + BB[i] + (1-xi*gamma)*AG[i] + (1-xi*gamma)*BG[i] + (1-xi)*GG[i]

    # proportion of infected humans
    # dXhdt = ((m*(a**2)*t*c*np.exp(-f*n))/(f + a*c*Xh))*(Xh * (1-Xh)) - r*X

    # proportion of infected mosquitos
    # dYmdt = (a*c*Xh) / (f + a*c*Xh)

    #Setting up integration
    #initial conditions
    Xh_prev = Xh[i-1]
    Ym_prev = Ym[i-1]
    y0 = [Xh_prev, Ym_prev]
    #time step
    t_int = np.linspace(i-1, i, 100)
    #solver
    sol = odeint(XYmodel, y0, t_int)


    #Solution set
    Xh_solved = sol[:, 0]
    Ym_solved = sol[:, 1]

    Xh[i] = Xh_solved[-1]
    Ym[i] = Ym_solved[-1]

    # Number of susceptible mosquitoes in a population
    Sm[i] = (1 - Ym[i]) * Nm

    # Number of infected mosquitoes in a population
    Im[i] = Nm * Ym[i]

    # STEP 5: Malaria (human) model

    # Number of suceptible humans in a population
    Sh[i] = (1 - Xh[i]) * Nh

    # Number of infected humans in a population
    Ih[i] = Nh * Xh[i]



print(Xh)
print(Ym)

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

Sg = scipy.transpose(Sg)
Eg = scipy.transpose(Eg)
Ig = scipy.transpose(Ig)
Rg = scipy.transpose(Rg)

pl.plot(tspan, Sg, 'r-', label='Sg')
pl.plot(tspan, Eg, 'g-', label='Eg')
pl.plot(tspan, Ig, 'b-', label='Ig')
pl.plot(tspan, Rg, 'o-', label='Rg')

pl.legend(loc='best')
pl.xlabel('Generation')
pl.ylabel('Population %')
pl.title('Gene Drive SEIR model')
pl.xticks(scipy.linspace(0, ngen, 11))  # sets ticks for x-axis

pl.figure(3)
Sm = scipy.transpose(Sm)
Im = scipy.transpose(Im)

pl.plot(tspan, Sm, 'r-', label='Sm')
pl.plot(tspan, Im, 'b-', label='Im')

pl.legend(loc='best')
pl.xlabel('Generation')
pl.ylabel('Number of Mosquitoes')
pl.title('Malaria SI model')
pl.xticks(scipy.linspace(0, ngen, 11))  # sets ticks for x-axis

pl.figure(4)
Sh = scipy.transpose(Sh)
Ih = scipy.transpose(Ih)

pl.plot(tspan, Sh, 'r-', label='Sh')
pl.plot(tspan, Ih, 'b-', label='Ih')

pl.legend(loc='best')
pl.xlabel('Generation')
pl.ylabel('Number of Humans')
pl.title('Human SI model')
pl.xticks(scipy.linspace(0, ngen, 11))  # sets ticks for x-axis


pl.show()
