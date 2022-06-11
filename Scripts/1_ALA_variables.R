
#### OBTAIN VARIABLES FROM ALA 

# Packages
#install.packages("dplyr")
#install.packages("devtools") # needed to install ALA
#devtools::install_github("AtlasOfLivingAustralia/ALA4R") # install ALA


# Activate packages ----
library(dplyr)
library(devtools)
library(ALA4R)


# Extract ----

points <- read.csv("../Data/5_Locations.csv")
head(points)
points_ala=points[,c(3,4)]
head(points_ala)

layers<-c(
          "el676", #Mean Annual Net Primary Productivity (tonnes/ha/yr) *
          "el1072", #	Fractional Cover - Bare Soil (2012-03-05)
          "el1081", #	Leaf Area Index (LAI) - 2012-03-05
          "el1080", #	Fraction of Photosynthetically Active Radiation (fPAR)
          "el715" , # Mean annual aridity index *
          "el722" , # Mean annual growth index C3 macrotherm *
          "el731" , # 	Mean annual growth index C3  mesotherm 
          "el751" # 	Mean annual growth index C4 megatherm
          )


extracted_vars <- intersect_points(points_ala, layers)
head(extracted_vars)
extracted_points= cbind(points, extracted_vars[,3:10])
head(extracted_points)
length(extracted_points)
write.csv(extracted_points,"../Data/20220107VegetationVariables.csv")
