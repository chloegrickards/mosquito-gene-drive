
#packages

#constants and parameters

#loop over generations
ngen = 100

for (i in 1:ngen) {
  # calculate genotypes
  # genotype(i) = genotypes(i-1)*[transition matrix]
  
  #Gene Drive SEIR
  #Sg = genotypes of AA, AB -- susceptible to gene drive infection
  #Eg = genotypes of AG -- carry a gene drive infection
  #Ig = genotypes of GG -- completely gene drive infected
  #Rg = genotypes of BB, BG -- cannot
  # possibility of moving BG --> Eg
  
  #Malaria SIR model
  #Sm = Susceptble mosquites
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

