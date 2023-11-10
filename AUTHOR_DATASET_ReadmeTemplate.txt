PhotonicRebelsreadme.txt file was generated on 2023-11-10 by Laura Bibiana Ospina-Rozo


######## GENERAL INFORMATION ########

1. Title of Dataset: Photonic Rebels Data Set

2. Author Information

	Laura Ospina-Rozo

Corresponding author
Institution: School of Biosciences - University of Melbourne - Parkville Victoria 3010 – Australia
Email: laura.ospinarozo@student.unimelb.edu.au Telephone: +61 0481882671 ORCID: https://orcid.org/0000-0002-1904-202X

	Iliana Medina
Institution: School of Biosciences - University of Melbourne - Parkville Victoria 3010 – Australia
Email: iliana.medina@unimelb.edu.au - ORCID: https://orcid.org/0000-0002-1021-5035 

	Andrew Hugall

Institution: School of Biosciences - University of Melbourne - Parkville Victoria 3010 – Australia & Sciences Department, Museum Victoria, GPO Box 666E, Melbourne, Victoria, 3001, Australia
Email: ahugall@museum.vic.gov.au 

	Katrina J. Rankin

Institution: School of Biosciences - University of Melbourne - Parkville Victoria 3010 – Australia 
Email:  katrina.rankin@unimelb.edu.au 

	Nicholas W. Roberts

Institution: School of Biological Sciences, University of Bristol, Bristol Life Sciences Building, 24 Tyndall Avenue, Bristol, BS8 1TQ, UK
Email: nicholas.roberts@bristol.ac.uk - ORCID: https://orcid.org/0000-0002-4540-6683 

	Ann Roberts

Institution: ARC Centre of Excellence for Transformative Meta-Optical Systems - School of Physics - University of Melbourne - Parkville Victoria 3010 – Australia
Email: ann.roberts@unimelb.edu.au - ORCID https://orcid.org/0000-0003-4295-9730 

	Andrew Mitchell
Institution: Australian Museum Research Institute, Australian Museum, 1 William Street, Sydney NSW 2010, Australia
Email: Andrew.Mitchell@australian.museum - ORCID: https://orcid.org/0000-0001-5022-5898 

	Chris A. M. Reid
Institution: Australian Museum Research Institute, Australian Museum, 1 William Street, Sydney NSW 2010, Australia
Email: Chris.Reid@australian.museum - ORCID: https://orcid.org/0000-0003-1899-9839 

	Adnan Moussalli
Institution: Sciences Department, Museum Victoria, GPO Box 666E, Melbourne, Victoria, 3001, Australia
Email: amoussalli@museum.vic.gov.au 

	Devi Stuart-Fox

Institution: School of Biosciences - University of Melbourne - Parkville Victoria 3010 – Australia
Email: d.stuart-fox@unimelb.edu.au - ORCID: https://orcid.org/0000-0003-3362-1412 


3. Date of data collection (approx): 2019 to 2022 

4. Geographic location of data collection: Melbourne Australia

5. Information about funding sources that supported the collection of the data: 
We also acknowledge the use of the Atlas of Living Australia, (https://ror.org/018n2ja79). This work was supported by funding from the Australian Research Council (grant numbers DP190102203, FT180100216, CE200100010) and a research grant from the Air Force Office of Scientific Research (AFOSR)/European Office of Aerospace Research and Development (EOARD) (grant number FA9550-19-1-7005) to N.W.R.



######## SHARING/ACCESS INFORMATION ########


2. This data is used in the publication: 
"Photonic rebels: reflectivity and polarization of Christmas beetles are not explained by climate"

3. Links to other publicly accessible locations of the data: 
Interactive code: https://lospinarozo.github.io/PhotonicRebelsCode/ 
Original code and data: https://github.com/lospinarozo/PhotonicRebelsCode

4. Recommended citation for this dataset: Please cite this data with the correspondent Dryad DOI




######## DATA & FILE OVERVIEW ########

1. File List: 

FromCode -- This folder contains data produced by our code that needs to be subsequently imported and used by other scripts.

ParametersMatLab -- This folder contains 15 txt files obtained by our matlab linearization and equalization tool. Since they all have the same format, we provide a general description below

RmdImages -- This folder contains 3 images used in the interactive version of the code

1_Reflectance_HRBatch.csv -- Raw reflectance profile of the subset of individuals from the Pretty Cool Beetles manuscript

2_Reflectance_FirstBatch.csv -- Raw reflectance profile from the rest of the individuals

3_Transmittance_HRBatch.csv -- Raw transmittance profile of the subset of individuals from the Pretty Cool Beetles manuscript

4_SunIrradiance.csv -- Spectral profile correspondent to the sun irradiance

5_Locations.csv -- Includes location data from each beetle speciment and their activity periods

6_VegetationVariables.csv -- Vegetation variables extracted by our code from ALA

7_ClimateVariables.csv --  Climate variables extracted from from continent-wide 0.05° grids of interpolated weather data from the Australian Gridded Climate Data (AGCD)/AWAP database (Australian Bureau of Meteorology official dataset for monthly gridded rainfall analysis)

8_Size.csv -- Sice of each specimen in cm

9_CodesAndSpecies.csv -- Species names for the each of the code names used in the other files

10_Polarization.csv -- RGB values extrated from the calibrated photographs under three different filters: VIS, LCP (left polarized) and RCP (right polarized)

12_knownGreyStdReflectance.csv -- Reflectance of the grey standard under each filter provided by the fabricant. We used the same value for the three of them since all photographs were taken with the lens in the visible spectral range

13_SuppTransmittancvsAbsorbance.csv -- Transmitivity and absorptivity in visible light for the subset of beetles studied in "Pretty Cool Beetles"

ConsReflEcolSpp.csv -- Consolidated file with the reflectivity values for each spectral band, the PC components (sumarizing ecological variables) and size

SupplementaryBeetleDataBase.csv -- Table including the species names and ANIC references for the samples used, when available. This is also available in the interactive version of the code

SupplementaryBeetleDataBase.xlsx -- same as the previous file in .xlsx format

XMAS_mat2b_bst2e.xml -- The BEAST xml file for generating phylogenetic trees from the Christmas beetle supermatrix data. 
Includes the sequence data, partition evolution model, and various prior topological and dating age constraints.

xmas_mat2b_bst2ef_set23nn2_pinct_med.tre -- Maximum clade credibility consensus of the tree set, made by TreeAnnotator. 
Includes node posterior support, age confidence interval and branch rate information.

XMAS_mat2b_bst2ef_set23nn2_pinct.nwk -- The set of 2000 trees randomly drawn from the BEAST posterior sample, used for comparative analysis. 
These have been pruned back to the focal taxa, had tip increment added and labels revised to match the trait dataset.



2. Relationship between files, if important: The data in all these files is combined by our phylogenetic comparison analysis. Details in the codmmented code.

3. Additional related data collected that was not included in the current data package: Not applicable, All data included

4. Are there multiple versions of the dataset? Yes, version control in Github link (above)

5. Relationship with the manuscript:  All these files were used to produce the main figures and tables as well as the supplementary materials in the manuscript.




######## METHODOLOGICAL INFORMATION ########

1. Description of methods used for collection/generation of data: 
All methods udes to collect this data are detailed in the correspondent manuscript or the commented code. 
Note that some of the data used in this manuscript was retrieved from a previous manuscript here abbreviated "Pretty cool beetles": Laura Ospina-Rozo, Jegadesan Subbiah, Ainsley Seago, Devi Stuart-Fox, Pretty Cool Beetles: Can Manipulation of Visible and Near-Infrared Sunlight Prevent Overheating?, Integrative Organismal Biology, Volume 4, Issue 1, 2022, obac036, https://doi.org/10.1093/iob/obac036

2. Methods for processing the data: 
To process this data we have used primarily the software R and R studio. The files in the folder ParametersMatLab where produced with a custom made MatLab code. 

3. Instrument- or software-specific information needed to interpret the data: NA

4. Standards and calibration information, if appropriate: Detailed in the manuscript

5. Environmental/experimental conditions: NA

6. Describe any quality-assurance procedures performed on the data: NA

7. People involved with sample collection, processing, analysis and/or submission:

Data on reflectivity, absorptivity and transmissivity was collected by Laura Ospina-Rozo. 
The raw RGB values from the calibrated photographs used to calculate polarization parameters were collected by Katrina Rankin. 
Phylogenetic was data produced by Andrew Hugall. 
Climate and vegetation data collected by Laura Ospina-Rozo with assistance from Iliana Medina. 




######## DATA-SPECIFIC INFORMATION BY FILE #######


__________________________________________________________________________

TITLE: Parameters MatLab\StandardParamsF01LCP.txt (only one example since they are all the same format)

1. Number of columns: 9 

2. Number of cases/rows: 5

3. Variable List: 
a, b, c and d are the parameters needed in the linearization/equalization equations. 
adjR2 is the correlation coeficient for each fitting done by MatLab. This values should always be around 0.99
seq is just a list of numbers
Channel is the chanel of the photograph, red, green, blue and gray average value. 
Photograph contains the codeof the photograph used for this standarization
Mode details in the manuscript.

4. Missing data codes: no

5. Specialized formats or other abbreviations used: 

__________________________________________________________________________

__________________________________________________________________________

TITLE: 1_Reflectance_HRBatch.csv

1. Number of cols: 57

2. Number of cases/rows: 1302

3. Variable List: First column (wl) is the wavelength and all the others columns are the reflectance (% in comparison to a white standard) for each beetle species

4. Missing data codes: no

5. Specialized formats or other abbreviations used: All the species abbreviations is detailed in the file 9_CodesAndSpecies.csv

__________________________________________________________________________

__________________________________________________________________________

TITLE: 2_Reflectance_FirstBatch

1. Number of variables: 237

2. Number of cases/rows: 901

3. Variable List: First column (wl) is the wavelength and all the others columns are the reflectance (% in comparison to a white standard) for each beetle species

4. Missing data codes: no

5. Specialized formats or other abbreviations used: All the species abbreviations is detailed in the file 9_CodesAndSpecies.csv

__________________________________________________________________________


__________________________________________________________________________

TITLE:3_Transmittance_HRBatch

1. Number of variables: 57

2. Number of cases/rows: 701

3. Variable List: First column (wl) is the wavelength and all the others columns are the transmittance (% in comparison to the full beam of light captured by the detector) for each beetle species

4. Missing data codes: no

5. Specialized formats or other abbreviations used: All the species abbreviations is detailed in the file 9_CodesAndSpecies.csv

__________________________________________________________________________


__________________________________________________________________________

TITLE: 4_SunIrradiance

1. Number of variables: 2

2. Number of cases/rows: 852

3. Variable List: First column (wl) is the wavelength and column 2 is the sun irradiance

4. Missing data codes: no

5. Specialized formats or other abbreviations used: none

__________________________________________________________________________


__________________________________________________________________________

TITLE: 5_Locations

1. Number of variables: 9

2. Number of cases/rows: 277

3. Variable List: 
Spp -- Specimen code
Reg -- ANIC reference if available, otherwise the code of the polarization photo for ID purposes. 
Latitude -- location
Longitude -- location
MonthMaxALA -- month with maximum ALA records for each species
MonthCollectionLabel -- month when each specimen was collected
MonthsActivityALA -- for each species we calculated the percentage of records for each month relative to the total, and only considered the months with a number of records equivalent to >10% of total records for the species. 
Batch	-- HR means this beetle was also studied in the manuscript "Pretty cool beetles". beetles labeled "original" in this column were studied only in this manuscript. 
NameinHR -- for specimens studied in "Pretty cool beetles" the codes used in that manuscript are provided for ID purposes. 

4. Missing data codes: no

5. Specialized formats or other abbreviations used: NA

__________________________________________________________________________


__________________________________________________________________________

TITLE: 6_VegetationVariables

1. Number of columns: 18

2. Number of cases/rows: 277

3. Variable List: 
spp -- beetle specimen code
reg -- ANIC reference if available, otherwise the code of the polarization photo for ID purposes. 
Latitude -- location
Longitude -- location
ALA.month -- Month in ALA with most records
picode -- Polarization Photo number for ID purposes
colection.month -- month when each specimen was collected
Batch	-- HR means this beetle was also studied in the manuscript "Pretty cool beetles". beetles labeled "original" in this column were studied only in this manuscript. 
NameinHR -- for specimens studied in "Pretty cool beetles" the codes used in that manuscript are provided for ID purposes. 
NameinTree -- species name as it appears in the phylogenetic tree

The following are the ecological variables related to vegetation cover recovered from ALA: 
NPPMean
fractionalCoverBareSoil20120305	
leafAreaIndexLAI20120305
fractionOfPhotosyntheticallyActiveRadiationFPAR
aridityIndexAnnualMean
growthIndexC3MacrothermPlantsAnnualMean
growthIndexC3MesothermPlantsAnnualMean
growthIndexC4MegathermPlantsAnnualMean


4. Missing data codes: no

5. Specialized formats or other abbreviations used: no

__________________________________________________________________________


__________________________________________________________________________

TITLE: 7_ClimateVariables

1. Number of columns: 19

2. Number of cases/rows: 276

3. Variable List: 
species	-- beetle specimen code
reg -- ANIC reference if available, otherwise the code of the polarization photo for ID purposes. 
lat -- Location
lon -- Location

The following correspond to the ecological variables extracted from bioClim:
avg_temp_over_35 avg_max_temp	avg_min_temp	avg_sol	avg_year_sol	avg_year_vpr	cloud_cover	avg_rr	avg_vpr	avg_temp_over_35_Coll	avg_max_temp_Coll	avg_min_temp_Coll	avg_sol_Coll	avg_rr_Coll	avg_vpr_Coll


4. Missing data codes: no

5. Specialized formats or other abbreviations used: 

__________________________________________________________________________


__________________________________________________________________________

TITLE: 8_Size

1. Number of columns: 2

2. Number of cases/rows: 276

3. Variable List: 
ind -- beetle specimen code
size -- length of the beetle in cm.

4. Missing data codes: no

5. Specialized formats or other abbreviations used: none

__________________________________________________________________________


_________________________________________________________________________

TITLE: 9_CodesAndSpecies

1. Number of columns: 2

2. Number of cases/rows: 276

3. Variable List: 
ind -- beetle specimen code
phylogeny_name -- species name as it is in the phylogenetic tree

4. Missing data codes: no

5. Specialized formats or other abbreviations used: none

__________________________________________________________________________


__________________________________________________________________________

TITLE: 10_Polarization

1. Number of columns: 14

2. Number of cases/rows: 4942

3. Variable List: 

ind  -- beetle specimen code
SpeciesName -- species name as it is in the phylogenetic tree
ANIC	-- ANIC code
Pic_code -- Number of the calibrated photograph
Filter	-- VIS == visible,  LCP == left handed polarized, RCP == right handed polarized 
tr	-- region sampled in the photograph. Either the elytron, pronotum or gray standard. 
label	-- channel R == red, G == green, B == blue. 
area	-- sampled area
mean	-- mean intensity value
min	-- min intensity value
max	-- max intensity value
camera_cat	-- configuration of the camera. This parameter is necessary because the photos were taken in different days, so there is one set of calibration parameters per day. Each photograph is calibrated according to its correspondent parameters. 
location	-- RoI == sampled in the region of interest, grey == sampled in the greay standard.
patch_ID	-- ID of the patch

4. Missing data codes: no

5. Specialized formats or other abbreviations used: none

__________________________________________________________________________


__________________________________________________________________________

TITLE: 12_knownGreyStdReflectance

1. Number of columns: 2

2. Number of cases/rows: 3

3. Variable List: 0.38446 is the reflectance of the grey standard included on each photograph. This value is given by the manufacturer. We used the same value for the three filters since the reflectance of the standards varies mostly according to spectral range and all images were taken in visible light. 

__________________________________________________________________________


__________________________________________________________________________

TITLE: 13_SuppTransmittancvsAbsorbance

1. Number of columns: 3

2. Number of cases/rows: 57

3. Variable List: 
Beetle	--  beetle specimen code
Td_VIS	--  transmissivity
Ab_VIS	--  absorptivity

Original data in "Pretty cool beetles"

4. Missing data codes: no

5. Specialized formats or other abbreviations used: none

__________________________________________________________________________


__________________________________________________________________________

TITLE: ConsReflEcolSpp

1. Number of columns: 8

2. Number of cases/rows: 48

3. Variable List: 

first row is the species name, but it does not contain header since this is the format required in R for the analysis, i.e. specis name has the be the row label and therefore can not have a header. 

TOT	-- total reflectivity (broadband 400 to 1700nm)
VIS	-- visible reflectivity (400 to 700 nm)
NIR	-- near infrared reflectivity (700 to 1700 nm)
Res	-- residuals from the regression between VIS and NIR reflectivity
PC1	-- first principal component sumarizing ecological variables
PC2	-- second principal component sumarizing ecological variables
Size	-- length in cm.


4. Missing data codes: no

5. Specialized formats or other abbreviations used: NA

__________________________________________________________________________


__________________________________________________________________________

TITLE: SupplementaryBeetleDataBase

1. Number of columns: 7

2. Number of cases/rows: 262

3. Variable List: 
Species_name	-- species name as it is in the phylogenetic tree
Code	-- beetle specimen code
ANIC reference	-- ANIC code if available
Latitude -- location for each specimen
Longitude -- location for each specimen
Batch	-- HR means this beetle was also studied in the manuscript "Pretty cool beetles". beetles labeled "original" in this column were studied only in this manuscript. 
NameinHR -- for specimens studied in "Pretty cool beetles" the codes used in that manuscript are provided for ID purposes. 

4. Missing data codes: no

5. Specialized formats or other abbreviations used: NA

__________________________________________________________________________


__________________________________________________________________________

TITLE: Phylogeny Data

Notes from the contents of folder phylogeny (outside of the Data folder, but still part of the R project) https://github.com/lospinarozo/PhotonicRebelsCode/tree/main/Phylogeny 
Note that the names of the files are the same as in this data set. 

XMAS_mat2b_bst2e.xml: 

The BEAST xml file for generating phylogenetic trees from the Christmas beetle supermatrix data. 
Includes the sequence data, partition evolution model, and various prior topological and dating age constraints.


XMAS_mat2b_bst2ef_set23nn2_pinct.nwk: 

The set of 2000 trees randomly drawn from the BEAST posterior sample, used for comparative analysis. 
These have been pruned back to the focal taxa, had tip increment added and labels revised to match the trait dataset.

xmas_mat2b_bst2ef_set23nn2_pinct_med.tre: 

Maximum clade credibility consensus of the tree set, made by TreeAnnotator. 
Includes node posterior support, age confidence interval and branch rate information.


__________________________________________________________________________