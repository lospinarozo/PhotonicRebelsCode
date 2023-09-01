Notes on the contents of folder https://github.com/lospinarozo/PhotonicRebelsCode/tree/main/Phylogeny 

PAT_script_v5.sh: 

Text file containing various scripts of bash and R code for assisting the assembling and checking of genetic sequence data supermatricies. 
The scripts include comment lines detailing use. Some code requires software MAFFT and IQTree.


XMAS_mat2b_bst2e.xml: 

The BEAST xml file for generating phylogenetic trees from the Christmas beetle supermatrix data. 
Includes the sequence data, partition evolution model, and various prior topological and dating age constraints.


XMAS_mat2b_bst2ef_set23nn2_pinct.nwk: 

The set of 2000 trees randomly drawn from the BEAST posterior sample, used for comparative analysis. 
These have been pruned back to the focal taxa, had tip increment added and labels revised to match the trait dataset.

xmas_mat2b_bst2ef_set23nn2_pinct_med.tre: 

Maximum clade credibility consensus of the tree set, made by TreeAnnotator. 
Includes node posterior support, age confidence interval and branch rate information.
