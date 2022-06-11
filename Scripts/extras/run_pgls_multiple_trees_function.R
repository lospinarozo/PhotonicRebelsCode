
# title: "PGLS with multiple phylogenetic trees"

# Packages

library(dplyr)
library(ape)
library(caper)
library(picante)
library(plyr)
library(phytools)



# Define functions

# the following function will have three arguments: 
  ## x = a .tree file with multiple trees. 
  ## model = The model we want to fit
  ## dataset =  A data frame with the predictors and response variables per each species


pgls_run=function(x, model, dataset) 
  
{tryCatch(
  {
    tree=check_and_fix_ultrametric(x) #Check if tree is ultrametric 
    dataset$TipLabel=rownames(dataset) #converts the row names into a column
    dataset_ord=match.phylo.data(tree, dataset) ### create dataset that matches order of tree tips for each tree
    treeord=dataset_ord$phy ### rename phylogeny vector as treeord
    dataord=as.data.frame(dataset_ord$dat) ### rename ordered dataset
    
    #Specify the types of variables in the dataset (numeric or character):
    dataord$Reflectance=as.numeric(as.character(dataord$Reflectance))
    
    dataord$PC1=as.numeric(as.character(dataord$PC1))
    dataord$PC2=as.numeric(as.character(dataord$PC2))
    
    dataord$size=as.numeric(as.character(dataord$size))
    rownames(dataord)<-NULL
    dataord<-dataord[,2:6]
    
    #Create comparative data for caper
    compdata <- comparative.data(phy = treeord, # phylogeny
                                 data = dataord, # explanatory and response var
                                 names.col = "TipLabel", 
                                 vcv = TRUE, 
                                 na.omit = FALSE, 
                                 warn.dropped = TRUE)
    
    #Model used
    modA <- pgls(model, data=compdata, param.CI = 0.95,lambda = "ML") 
    # We are using max likelihood to estimate the lambda (phylogenetic signal) of the link, and with that info it runs de model 
    
   
    if (nrow(model.results[["coefficients"]]) == 5) #4 predictors plus intercept
      
    {
      
      coef.res1[1,1]<-coefficients(model.results)[2,1]  ### estimates for predictor 1
      coef.res1[1,2]<-coefficients(model.results)[2,3]
      coef.res1[1,3]<-coefficients(model.results)[2,4]
      
      coef.res2[1,1]<-coefficients(model.results)[3,1]  ### estimates for predictor 2
      coef.res2[1,2]<-coefficients(model.results)[3,3]
      coef.res2[1,3]<-coefficients(model.results)[3,4]
      
      coef.res3[1,1]<-coefficients(model.results)[4,1]  ### estimates for predictor 3
      coef.res3[1,2]<-coefficients(model.results)[4,3]
      coef.res3[1,3]<-coefficients(model.results)[4,4]
      
      coef.res4[1,1]<-coefficients(model.results)[5,1]  ### estimates for predictor 4
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
      
      
      
      total_results=cbind(final.results1, final.results2, final.results3,final.results4)
      
      return(total_results)
      
    }
    
    if (nrow(model.results[["coefficients"]]) == 6) #5 predictors plus intercept
    {
      coef.res1[1,1]<-coefficients(model.results)[2,1]  ### estimates for predictor 1
      coef.res1[1,2]<-coefficients(model.results)[2,3]
      coef.res1[1,3]<-coefficients(model.results)[2,4]
      
      coef.res2[1,1]<-coefficients(model.results)[3,1] #### estimates for predictor 2
      coef.res2[1,2]<-coefficients(model.results)[3,3]
      coef.res2[1,3]<-coefficients(model.results)[3,4]
      
      coef.res3[1,1]<-coefficients(model.results)[4,1] #### estimates for predictor 3
      coef.res3[1,2]<-coefficients(model.results)[4,3]
      coef.res3[1,3]<-coefficients(model.results)[4,4]
      
      coef.res4[1,1]<-coefficients(model.results)[5,1] #### estimates for predictor 4
      coef.res4[1,2]<-coefficients(model.results)[5,3]
      coef.res4[1,3]<-coefficients(model.results)[5,4]
      
      coef.res5[1,1]<-coefficients(model.results)[6,1] #### estimates for predictor 5
      coef.res5[1,2]<-coefficients(model.results)[6,3]
      coef.res5[1,3]<-coefficients(model.results)[6,4]
      
      
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
      
      total_results=cbind(final.results1, final.results2,
                          final.results3,final.results4,
                          final.results5)
      
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
      
      total_results=cbind(final.results1, final.results2,
                          final.results3,final.results4,
                          final.results5,final.results6)
      
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


## Empty vectors

model.results<-list() #the entire model for each run, stored in a list
###list with coefficients and estimates for each of the predictors in the model, create one for each predictor in your model ######
coef.res1<-matrix(data=NA, nrow=1, ncol=3) 
coef.res2<-matrix(data=NA, nrow=1, ncol=3)
coef.res3<-matrix(data=NA, nrow=1, ncol=3)
coef.res4<-matrix(data=NA, nrow=1, ncol=3)
coef.res5<-matrix(data=NA, nrow=1, ncol=3)


