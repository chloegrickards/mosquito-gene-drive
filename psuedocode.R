
#packages

#constants and parameters

#loop over generations
ngen = 100

for (i in 1:ngen) {
  # calculate genotypes
  # genotype(i) = genotypes(i-1)*[transition matrix]
  
  #Gene Drive SEIR
  #Sg(i) = genotypes of AA, AB -- susceptible to gene drive infection
  #Eg(i) = genotypes of AG -- carry a gene drive infection
  #Ig(i) = genotypes of GG -- completely gene drive infected
  #Rg(i) = genotypes of BB, BG -- cannot
  # possibility of moving BG --> Eg
  
  #Malaria SIR model
  #birth rate of susceptible mosquitoes = AA(i) + AB(i) + BB(i) + (some %)*AG + (some %)*BG(i) + (some %)*GG(i)
  # = r in the logistic growth equations
  
  #Option 2:
  # modify dXdt
  
  #dXdt from Jamie's equations
  #NOT including mu for birth rate because we already did that above with the genotypes transition matrix
  #dXdt = (logistic growth function, with birth rates = r) - beta*x*y - death rate

  #X(i) = integrate(logistic growth function, with birth rates = r) + integrate(-beta*x*y - death rate)
  
  
  
  #X(i) = integrate(dXdt)
  
  
  
  
  #Sm(i) = Susceptble mosquites
  # some% = susceptibility fr being heterozygous
  # depends on these susceptible phenotypes: AA, AB, BB, some% AG, some% BG
  # + birth
  # - death
  # - infections
  
  # Im = Infectious mosquitoes
  # + infections
  # - death
  
  #Rm = "Recovered" but actully Resistant mosquitoes
  # depends on these resistant phenotypes: GG, (1-some%) AG, (1-some%) BG
  # + birth = total mosquito birth - susceptible mosquito birth
  # - death

  # Human R0 model 
  #dX = R0
  #dX = mab(Im)(1-X) - rX
  
}

