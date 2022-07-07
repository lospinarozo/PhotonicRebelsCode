################################### About run_mcmc ###################################
# Date: 2021 Nov 01
# This function is custom written to run MCMCglmm on all 1300 trees of beetles 
# We include the function here to keep the main RMarkdown tidy.

# This function includes only 1 random factor:
# 1) phylogeny ("phylogeny_name")
# In this function, we first run 1300 iteration (nitt = 1300, burnin = 0, thin = 1) to kick it running. 
# The last iteration is saved as the starting values of the run on the next tree.
# Next, we repeat 1300 times (nitt = 2000, burnin = 1999, thin = 1.) 
# but discard the first 300 times as burn in.
# The method is adapted from Stuart-Fox et al. 2021 (Ecology Letters, 24, 2207– 2218).

######################################################################################


run_mcmc_byspp <- function(trees, model_chosen, phylo_data){
  
  phylo_data$phylogeny_name = factor(phylo_data$phylogeny_name) 
  species = as.data.frame(unique(phylo_data$phylogeny_name))
  row.names(species) = species[,1]
  temp <- name.check(trees[[1]], species) 
  tree.pruned <- trees[[1]] 
  AinvTA <- inverseA(tree.pruned)$Ainv 
  
  # This MCMCglmm() here is jsut to set a start to kick it running
  m.res = MCMCglmm(model_chosen , 
                   random = ~ phylogeny_name, 
                   # here I account for phylogeny info and species as random factors
                   prior = prior, 
                   ginverse = list(animal= AinvTA), 
                   data = phylo_data, 
                   family = 'gaussian', # depends on the distribution of the parameter you want to estimate
                   nitt = 1300, burnin = 0, thin = 1, # just to get the chains start 
                   # Here we set the nitt number as the loop number below.
                   pl = T,
                   pr = F,
                   verbose = F,
                   singular.ok = T)
  #This is just to set the starting frame, parameters for the following chains
  
  #setting the name of the output file (ie results will be stored in m.res)
  m0 = m.res # we set the result from the above MCMCglmm() as the start point of the following loop
  
  for(i in 1:1300){
    
    species = as.data.frame(unique(phylo_data$phylogeny_name))
    row.names(species) = species[,1] 
    temp <- name.check(trees[[i]], species) 
    tree.pruned <- trees[[i]] # put it i because we want to take different trees for each time
    treeok = tree.pruned # Then we save it as treeok which will be used for incerseA() matrix
    
    #treeok=check_and_fix_ultrametric(treeok)
    
    AinvTA <- inverseA(treeok)$Ainv
    
    start <- list(Liab = m0$Liab[1,], 
                  R = m0$VCV[1,3], 
                  G = list(G1 = m0$VCV[1,1], G2 = m0$VCV[1,2]))
    # This start here is to tell the below MCMCglmm() where to start. 
    #This start is derived from the results of the previous round as we use m0.
    
    # Note:
    # VCV = phylogenetic variance-covariance 
    m0 = MCMCglmm(model_chosen,
                  random = ~animal + phylogeny_name,
                  prior = prior, # distribution where i am gonna sample from
                  ginverse = list(animal = AinvTA),
                  data = phylo_data,
                  family = 'gaussian', 
                  nitt = 2000, burnin = 1999, thin = 1, 
                  #run for 2000 iterations and then save the last one and move to next tree
                  start = start, # just tell you where to start, but varies between each round. sample from posterior distribution
                  pl = T,
                  pr = F,
                  verbose = F,
                  singular.ok = T)
    
    #this sets the burnin, in this case the first 300 trees
    if(i > 300){ 
      
      m.res$VCV[i - 300,] <- m0$VCV[1,]
      m.res$Sol[i - 300,] <- m0$Sol[1,]
      m.res$Liab[i - 300,] <- m0$Liab[1,]
    }
    
    print(i) #shows progress  
    
    
  }
  
  return(m.res) #m.res includes posterior distributions of parameters taking phylogenetic uncertainty into account.
  #plot(m.res)
}