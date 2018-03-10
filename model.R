
#Something to import params

#Create empty vectors to fill
genotypes = zeros(ic.length, tspan)

# I need to check how things are indexed
for(i in tspan) {
  genotypes[i + 1] = transition_matrix*genotypes[i]
  
  

}
