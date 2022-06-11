### Install packages ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree") #need BiocManager to install, so run previous command 

install.packages("phylogram")
install.packages("devtools")
devtools::install_github("eliocamp/ggnewscale") #Needs Devtools to install this, install if you don't have it

#### Visualising data ####
library(phylogram)
library(ggtree)
library(devtools)
library(ggnewscale)
library(ape)
library(phangorn)
library(caper)

############### JUST RUN ##################################
# Generation of data, ignore, just run from line 19 to 30 #
newick <- "(A,((B,(C,(H,(I)))),(D,(E,(F,G)))));"
new_l <- 0.2
tree <- as.phylo(read.dendrogram(text = newick))
ref <- as.data.frame(as.numeric(as.character(c(80,90,75,70,75,35,39,38,40))))
colnames(ref) <- c('ref')
rownames(ref)<- tree$tip.label
temp<- as.data.frame(as.numeric(as.character(c(35,32,28,24,25,14,23,20,20))))
colnames(temp) <- c('temp')
rownames(temp)<- tree$tip.label
data_ordered <- cbind(ref,temp, rownames(temp))

plot(tree) #Your tree is an R object of class phylogeny, with tip labels and branch lengths

# START EXCERCISE HERE, AFTER RUNNING LINES ABOVE ####
#############################
#First step, lets see how data is distributed.
#Plotting data using ggtree, object tree is our tree, ref is reflectance 
#and temp is temperature values. In real life you'll need to make sure tree tip names match rownames of data.
# ggtree has great documentation! https://guangchuangyu.github.io/ggtree-book/short-introduction-to-r.html

phylo <- ggtree(tree) + geom_tiplab() #create phylogeny object for ggtree, if you type phylo you'll see your tree

p1 <- gheatmap(phylo, ref, offset=0.1, width= 0.1, low="gold",high = "blue") #first plot one variable

p2 <- p1 + new_scale_fill() # add new scale for second variable

gheatmap(p2, temp, offset=0.8, width=0.1, low = "gold", high = "red") #plot second variable

##### Run PGLS to test association #######
# generate comparative dataset with tree and ordered data
comp_data <-  comparative.data(phy = tree, data = data_ordered,
                                          names.col = "rownames(temp)", vcv = TRUE, 
                                          na.omit = FALSE, warn.dropped = TRUE)

model1 <- pgls(ref ~ temp, data=comp_data, param.CI = 0.95,lambda = 0.001)
summary(model1)

######### GAME #######
##### Go to: https://b.socrative.com/login/student/ 
# ROOM: MEDINA4417
# Just one person per team is enough, enter the answers to compete with the other team.

# Question 1: If you had to guess, what would be the phylogenetic signal of NIR reflectance?
# A. 0.001 B.0.1 C.0.9

# Question 2: What do you think is the most likely ancestral state for NIR, for the whole clade? 
# A. High reflectance B. Intermediate reflectance C. Low reflectance

# Question 3: Which one is true about model 1?
# A. There is a significant and negative relationship between temperature and NIR
# B. In places where temperatures are close to zero NIR reflectance is expected to be around 6.8%
# C. For every increase in two degrees, reflectance increases by ~4.5%
# D. For every increase in one percent reflectance, temperature increases by 2.29.

# Question 4: 
# Imagine your phylogenetic tree is re-shuffled and
# species are now placed randomly in the tree (species still keep their traits).
# Which are more likely to be your new results if you run model 1 again using the random phylogeny?
# A. B= 2.29, P=0.006 B. B=2.36, P=0.0002 C. B=2.12, P=0.021

#Question 5:
# If you wanted to test whether there were sex differences between males and females in NIR, 
#using a similar dataset, what would you do?
# A. Run a normal linear regression using sex as a predictor
# B. Do a PGLS using sex as a predictor
# C. None of the above

#Question 6:
# Which are more likely to be evil, true elves or dwarves?
http://phylonetworks.blogspot.com/2016/03/the-phylogeny-of-elves-and-other.html 

#Question 7:
# A model with the following command would give me a result closer to a common linear regression:
#model1 <- pgls(ref ~ temp, data=comp_data, param.CI = 0.95,lambda = 0.3)
# True or False?
