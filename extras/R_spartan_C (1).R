#source('functions_pgls.R')
library(dplyr)
library(ape)
library(caper)
library(picante)
library(plyr)
library(phytools)


###Load functions###
pgls_run=function(x, model, dataset)
  
{tryCatch(
  {
    #dataset= env_nest_size_range
    #x=trees[[5]] #To test the function
    tree=check_and_fix_ultrametric(x) #Check to see if ultrametric
    rownames(dataset)=dataset$TipLabel #name of column with spp
    dataset_ord=match.phylo.data(tree, dataset) ### create dataset that matches order of tree tips for each tree
    treeord=dataset_ord$phy ### rename phylogeny
    dataord=as.data.frame(dataset_ord$dat) ### rename ordered dataset
    
    #Specify the types of variables in the dataset
    dataord$PC1_T=as.numeric(as.character(dataord$PC1_T))
    dataord$PC2_T=as.numeric(as.character(dataord$PC2_T))
    dataord$PC1_P=as.numeric(as.character(dataord$PC1_P))
    dataord$PC2_P=as.numeric(as.character(dataord$PC2_P))
    dataord$range.size.m2=as.numeric(as.character(dataord$range.size.m2))
    dataord$BodyMass.Value=as.numeric(as.character(dataord$BodyMass.Value))
    dataord$nest_cat_cup2=as.factor(dataord$nest_cat_cup2)
    dataord$location=as.factor(dataord$location)
    dataord$nest_cat_cup2 = factor(dataord$nest_cat_cup2, levels=c("Open", "Domed","Cavity"))
    
    #dataord$logweight=log10(as.numeric(as.character(dataord$weight)))
    #dataord$PC1_NICHE1000B=(as.numeric(as.character(dataord$PC1_NICHE1000B)))
    
    #Node labels
    #treeord$node.label=seq(1,276)
    
    #Create comparative data for caper
    compdata <- comparative.data(phy = treeord, data = dataord,
                                 names.col = "TipLabel", vcv = TRUE, 
                                 na.omit = FALSE, warn.dropped = TRUE)
    
    #Model used
    #model <- weightedNIR ~ avgYearSol + weightedVIS
    modA <- pgls(model, data=compdata, param.CI = 0.95,lambda = "ML") 
    #model 1 We are using maximum likelihood to estimate the lambda (phylogenetic signal) of the link, and with that info it runs de model 
    
    model.results<-summary(modA)
    
    if (nrow(model.results[["coefficients"]]) == 5) 
      
    {
      
      coef.res1[1,1]<-coefficients(model.results)[2,1]  ### estimates for predictor 1
      coef.res1[1,2]<-coefficients(model.results)[2,3]
      coef.res1[1,3]<-coefficients(model.results)[2,4]
      
      coef.res2[1,1]<-coefficients(model.results)[3,1]  ### estimates for predictor 1
      coef.res2[1,2]<-coefficients(model.results)[3,3]
      coef.res2[1,3]<-coefficients(model.results)[3,4]
      
      coef.res3[1,1]<-coefficients(model.results)[4,1]  ### estimates for predictor 1
      coef.res3[1,2]<-coefficients(model.results)[4,3]
      coef.res3[1,3]<-coefficients(model.results)[4,4]
      # 
      coef.res4[1,1]<-coefficients(model.results)[5,1]  ### estimates for predictor 1
      coef.res4[1,2]<-coefficients(model.results)[5,3]
      coef.res4[1,3]<-coefficients(model.results)[5,4]
      
      final.results1<-data.frame(coef.res1)  #### result for predictor 1
      colnames(final.results1)<-c("estimateP1", "ts1", "ps1")
      
      final.results2<-data.frame(coef.res2)  #### result for predictor 2
      colnames(final.results2)<-c("estimateP2", "tP2", "pP2")
      
      final.results3<-data.frame(coef.res3)  #### result for predictor 3
      colnames(final.results3)<-c("estimateP3", "tP3", "pP3")
      
      final.results4<-data.frame(coef.res4)  #### result for predictor 4
      colnames(final.results4)<-c("estimateP4", "tP4", "pP4")
      
      #total_results=cbind(final.results1, final.results2)
      
       total_results=cbind(final.results1, final.results2, final.results3,final.results4)
      
      return(total_results)
      
    }
    
    if (nrow(model.results[["coefficients"]]) == 7) #Six predictors plus intercept
    {
      coef.res1[1,1]<-coefficients(model.results)[2,1]  ### estimates for predictor 1
      coef.res1[1,2]<-coefficients(model.results)[2,3]
      coef.res1[1,3]<-coefficients(model.results)[2,4]
      
      coef.res2[1,1]<-coefficients(model.results)[3,1] #### estimates for predictor 2
      coef.res2[1,2]<-coefficients(model.results)[3,3]
      coef.res2[1,3]<-coefficients(model.results)[3,4]
      
      coef.res3[1,1]<-coefficients(model.results)[4,1] #### estimates for predictor 2
      coef.res3[1,2]<-coefficients(model.results)[4,3]
      coef.res3[1,3]<-coefficients(model.results)[4,4]
      
      coef.res4[1,1]<-coefficients(model.results)[5,1] #### estimates for predictor 2
      coef.res4[1,2]<-coefficients(model.results)[5,3]
      coef.res4[1,3]<-coefficients(model.results)[5,4]
      
      coef.res5[1,1]<-coefficients(model.results)[6,1] #### estimates for predictor 2
      coef.res5[1,2]<-coefficients(model.results)[6,3]
      coef.res5[1,3]<-coefficients(model.results)[6,4]
      
      coef.res6[1,1]<-coefficients(model.results)[7,1] #### estimates for predictor 2
      coef.res6[1,2]<-coefficients(model.results)[7,3]
      coef.res6[1,3]<-coefficients(model.results)[7,4]
      
      
      
      final.results1<-data.frame(coef.res1)  #### result for predictor 1
      colnames(final.results1)<-c("estimateP1", "ts1", "ps1")
      
      final.results2<-data.frame(coef.res2)  #### result for predictor 2
      colnames(final.results2)<-c("estimateP2", "tP2", "pP2")
      
      final.results3<-data.frame(coef.res3)  #### result for predictor 3
      colnames(final.results3)<-c("estimateP3", "tP3", "pP3")
      
      final.results4<-data.frame(coef.res4)  #### result for predictor 4
      colnames(final.results4)<-c("estimateP4", "tP4", "pP4")
      
      final.results5<-data.frame(coef.res5)  #### result for predictor 5
      colnames(final.results5)<-c("estimateP5", "tP5", "pP5")
      
      final.results6<-data.frame(coef.res6)  #### result for predictor 6
      colnames(final.results6)<-c("estimateP6", "tP6", "pP6")
      
      total_results=cbind(final.results1, final.results2, final.results3,final.results4,final.results5,final.results6)
      
      return(total_results)
      
    }
    
    else { return('Change dimentions, too many')}
    
  }
  
  
  , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}


pgls_run2=function(x, model, dataset)
  
{tryCatch(
  {
    #dataset= env_nest_size_range
    #x=trees[[5]] #To test the function
    tree=check_and_fix_ultrametric(x) #Check to see if ultrametric
    rownames(dataset)=dataset$TipLabel #name of column with spp
    dataset_ord=match.phylo.data(tree, dataset) ### create dataset that matches order of tree tips for each tree
    treeord=dataset_ord$phy ### rename phylogeny
    dataord=as.data.frame(dataset_ord$dat) ### rename ordered dataset
    
    #Specify the types of variables in the dataset
    #dataord$PC1_T=as.numeric(as.character(dataord$PC1_T))
    #dataord$PC2_T=as.numeric(as.character(dataord$PC2_T))
    #dataord$PC1_P=as.numeric(as.character(dataord$PC1_P))
    #dataord$PC2_P=as.numeric(as.character(dataord$PC2_P))
    dataord$range.size.m2=as.numeric(as.character(dataord$range.size.m2))
    dataord$BodyMass.Value=as.numeric(as.character(dataord$BodyMass.Value))
    dataord$nest_cat_cup2=as.factor(dataord$nest_cat_cup2)
    dataord$location=as.factor(dataord$location)
    dataord$nest_cat_cup2 = factor(dataord$nest_cat_cup2, levels=c("Open", "Domed","Cavity"))
    
    #dataord$logweight=log10(as.numeric(as.character(dataord$weight)))
    #dataord$PC1_NICHE1000B=(as.numeric(as.character(dataord$PC1_NICHE1000B)))
    
    #Node labels
    #treeord$node.label=seq(1,276)
    
    #Create comparative data for caper
    compdata <- comparative.data(phy = treeord, data = dataord,
                                 names.col = "TipLabel", vcv = TRUE, 
                                 na.omit = FALSE, warn.dropped = TRUE)
    
    #Model used
    #model <- weightedNIR ~ avgYearSol + weightedVIS
    modA <- pgls(model, data=compdata, param.CI = 0.95,lambda = "ML") 
    #model 1 We are using maximum likelihood to estimate the lambda (phylogenetic signal) of the link, and with that info it runs de model 
    
    model.results<-summary(modA)
    
    if (nrow(model.results[["coefficients"]]) == 5) 
      
    {
      
      coef.res1[1,1]<-coefficients(model.results)[2,1]  ### estimates for predictor 1
      coef.res1[1,2]<-coefficients(model.results)[2,3]
      coef.res1[1,3]<-coefficients(model.results)[2,4]
      
      coef.res2[1,1]<-coefficients(model.results)[3,1]  ### estimates for predictor 1
      coef.res2[1,2]<-coefficients(model.results)[3,3]
      coef.res2[1,3]<-coefficients(model.results)[3,4]
      
      coef.res3[1,1]<-coefficients(model.results)[4,1]  ### estimates for predictor 1
      coef.res3[1,2]<-coefficients(model.results)[4,3]
      coef.res3[1,3]<-coefficients(model.results)[4,4]
      # 
      coef.res4[1,1]<-coefficients(model.results)[5,1]  ### estimates for predictor 1
      coef.res4[1,2]<-coefficients(model.results)[5,3]
      coef.res4[1,3]<-coefficients(model.results)[5,4]
      
      final.results1<-data.frame(coef.res1)  #### result for predictor 1
      colnames(final.results1)<-c("estimateP1", "ts1", "ps1")
      
      final.results2<-data.frame(coef.res2)  #### result for predictor 2
      colnames(final.results2)<-c("estimateP2", "tP2", "pP2")
      
      final.results3<-data.frame(coef.res3)  #### result for predictor 3
      colnames(final.results3)<-c("estimateP3", "tP3", "pP3")
      
      final.results4<-data.frame(coef.res4)  #### result for predictor 4
      colnames(final.results4)<-c("estimateP4", "tP4", "pP4")
      
      #total_results=cbind(final.results1, final.results2)
      
       total_results=cbind(final.results1, final.results2, final.results3,final.results4)
      
      return(total_results)
      
    }
    
    if (nrow(model.results[["coefficients"]]) == 7) #Six predictors plus intercept
    {
      coef.res1[1,1]<-coefficients(model.results)[2,1]  ### estimates for predictor 1
      coef.res1[1,2]<-coefficients(model.results)[2,3]
      coef.res1[1,3]<-coefficients(model.results)[2,4]
      
      coef.res2[1,1]<-coefficients(model.results)[3,1] #### estimates for predictor 2
      coef.res2[1,2]<-coefficients(model.results)[3,3]
      coef.res2[1,3]<-coefficients(model.results)[3,4]
      
      coef.res3[1,1]<-coefficients(model.results)[4,1] #### estimates for predictor 2
      coef.res3[1,2]<-coefficients(model.results)[4,3]
      coef.res3[1,3]<-coefficients(model.results)[4,4]
      
      coef.res4[1,1]<-coefficients(model.results)[5,1] #### estimates for predictor 2
      coef.res4[1,2]<-coefficients(model.results)[5,3]
      coef.res4[1,3]<-coefficients(model.results)[5,4]
      
      coef.res5[1,1]<-coefficients(model.results)[6,1] #### estimates for predictor 2
      coef.res5[1,2]<-coefficients(model.results)[6,3]
      coef.res5[1,3]<-coefficients(model.results)[6,4]
      
      coef.res6[1,1]<-coefficients(model.results)[7,1] #### estimates for predictor 2
      coef.res6[1,2]<-coefficients(model.results)[7,3]
      coef.res6[1,3]<-coefficients(model.results)[7,4]
      
      
      
      final.results1<-data.frame(coef.res1)  #### result for predictor 1
      colnames(final.results1)<-c("estimateP1", "ts1", "ps1")
      
      final.results2<-data.frame(coef.res2)  #### result for predictor 2
      colnames(final.results2)<-c("estimateP2", "tP2", "pP2")
      
      final.results3<-data.frame(coef.res3)  #### result for predictor 3
      colnames(final.results3)<-c("estimateP3", "tP3", "pP3")
      
      final.results4<-data.frame(coef.res4)  #### result for predictor 4
      colnames(final.results4)<-c("estimateP4", "tP4", "pP4")
      
      final.results5<-data.frame(coef.res5)  #### result for predictor 5
      colnames(final.results5)<-c("estimateP5", "tP5", "pP5")
      
      final.results6<-data.frame(coef.res6)  #### result for predictor 6
      colnames(final.results6)<-c("estimateP6", "tP6", "pP6")
      
      total_results=cbind(final.results1, final.results2, final.results3,final.results4,final.results5,final.results6)
      
      return(total_results)
      
    }
    
    else { return('Change dimentions, too many')}
    
  }
  
  
  , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

check_and_fix_ultrametric <- function(phy){
  
  if (!is.ultrametric(phy)){
    
    vv <- vcv.phylo(phy)
    dx <- diag(vv)
    mxx <- max(dx) - dx
    for (i in 1:length(mxx)){
      phy$edge.length[phy$edge[,2] == i] <- phy$edge.length[phy$edge[,2] == i] + mxx[i]
    }
    if (!is.ultrametric(phy)){
      stop("Ultrametric fix failed\n")
    }	
  }
  
  return(phy)
}

######## Analysis ###################
model.results<-list() #the entire model for each run, stored in a list
###list with coefficients and estimates for each of the predictors in the model, create as many as predictors you have ######
coef.res1<-matrix(data=NA, nrow=1, ncol=3) 
coef.res2<-matrix(data=NA, nrow=1, ncol=3)
coef.res3<-matrix(data=NA, nrow=1, ncol=3)
coef.res4<-matrix(data=NA, nrow=1, ncol=3)

env_nest_size_red= readRDS('env_nest_size_red.R')
trees= readRDS('1300_trees.rds')[51:100]
env_nest_size_red$nest_cat_cup2 = factor(env_nest_size_red$nest_cat_cup2, levels=c("Open", "Domed","Cavity"))

#Model 2: Are there differences in the variation in precipitation across the range?
model.results<-list() #the entire model for each run, stored in a list
###list with coefficients and estimates for each of the predictors in the model, create as many as predictors you have ######
coef.res1<-matrix(data=NA, nrow=1, ncol=3) 
coef.res2<-matrix(data=NA, nrow=1, ncol=3)
coef.res3<-matrix(data=NA, nrow=1, ncol=3)
coef.res4<-matrix(data=NA, nrow=1, ncol=3)


#M <- lm(PC1_P ~ nest_cat_cup2 + BodyMass.Value + location,data=env_nest_size_red)
#check_model(M) #Can have all variables in the same models

model2 <- PC1_P ~ nest_cat_cup2 + log(BodyMass.Value)

runs_1000=lapply(trees,pgls_run,model=model2,dataset=env_nest_size_red) ##half the time with 4 cores
df2 <- ldply (runs_1000, data.frame)[1:50,] #This will turn results into a dataframe, currently only 5 trees
#HPDinterval(as.mcmc(df2[2:7])) #We n

saveRDS(df2,'Model_precipitation_50B_log.R')

