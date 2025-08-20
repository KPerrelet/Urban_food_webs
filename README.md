# Urban_food_webs
Data and code for the article titled "Urbanization drives the decoupling,  simplification, and homogenization of aquatic and terrestrial food webs"

The information contained in this repository is referenced in the manuscript. 


*****************************************************

covariates.csv - Covariate values extracted for all ground-level sites. As many sites are private, their location cannot be disclosed publically.

While this table cannot be directly recreated as doing so would require the coordinates of the sites, the calculation of the covariate can be followed on the script covariates_calculation.R.



foodweb_metrics.csv - food web metrics for the 54 sites. 

This table can be recreated using the foodweb_metrics_calculation.R script.



metaweb.csv - regional metaweb based on species detected in the study. 

metaweb_processed.csv - processed version of the regional metaweb for later handling. 

Rows = interaction, columns = interaction source and target. 

The metaweb_processed.csv table can be recreated using the prepare_metaweb.R script. 



taxamat.csv - taxa table for all sites. 

Rows = sites, columns = species. Values indicate number of reads per OTU (after filtering, contamination removal, and rarefaction). 



taxo.csv - Taxonomy table for all taxa. 

Rows = taxa, columns = phylogeny. 



Raw sequences for the data generated in this study are available on the European Nucleotide Archive under accession project numbers PRJEB88517.

Data is also available on Zenodo at [10.5281/zenodo.15263460](https://doi.org/10.5281/zenodo.15263459)


*****************************************************


metaweb_fig1.R - generates the network (metaweb) in Figure 1. 

fig2.R - generate Figure 2. 

fig3.R - generate Figure 3. 

fig4.R - generate Figure 4. 

figS2.R - generate Figure S2. 

figS4.R - generate Figure S4. 

figS5.R - generate Figure S5. 

figS6.R - generate Figure S6. 

figS7.R - generate Figure S7. 

figS8.R - generate Figure S8. 

figS8.R - generate Figure S9. 

foodweb_metrics_calculation.R - calculate food web metrics for all sites, including those resulting from null models.

prepare_metaweb.R - execute minor modifications to the regional metaweb for smooth processing. 

covariates_calculation.R - calculate covariates for all sites



*****************************************************

Data dictionary: 

site_no                   |  paired site ID. A green roof site has the same ID as the ground-level site only if they are paired. 

natural_banks             |  percentage of natural banks around the pond (%)

marcophyte                |  marcophyte coverage (%)

distance_forest           |  distance to the forest (m) 

distance_pond             |  distance to the closest pond (m) 

terrestrial_connectivity  |  patch cohesion index, calculated using landscapemetrics

aquatic_connectivity      |  mean distance between ponds in a 500 m buffer, calculated using landscapemetrics

green50                   |  number of green pixels in a 50 m buffer around the pond edge

grey500_frac              |  fraction of grey surface in a 500 m buffer

complementarity           |  landscape diversity metric, calculated using the lsm_l_joinent function in landscapemetrics

overwarming               |  mean deviation of nighttime temperature compare to the regional average (K)

no2                       |  annual mean air NO2 concentration 

population                |  population density in a 500 m buffer

vegetation_height_sd_50   |  standard deviation of the vegetation height in a 50 m buffer (m)

vegetation_density_50     |  vegetation density in a 50 m buffer 
