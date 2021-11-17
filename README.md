# SeedBankBias

Principal Investigator:  Julie Larson, University of Colorado Boulder, julie.e.larson@colorado.edu

Date of data collection   2017-2018.

Geographic location of data collection:  Boulder, CO, USA [See associated manuscript for more information on site locations]

Funding sources or sponsorship that supported the collection of the data:  
This work was supported by the City of Boulder Open Space & Mountain Parks Funded Research Program, Univ. of Colorado Undergraduate Research Opportunities Program, and a graduate research grant from the Dept. of Ecology & Evolutionary Biology (Univ. of Colorado)


--------------------------
SHARING/ACCESS INFORMATION
-------------------------- 

Licenses/restrictions placed on the data, or limitations of reuse:  none

Recommended citation for the data:  This data is being archived in anticipation of manuscript publication in the journal: Ecology  
Search for: Larson, J.E & K.N. Suding. Seed bank bias: Differential tracking of functional traits in the seed bank and vegetation across a gradient. Ecology. 



--------------------
PROJECT OVERVIEW
--------------------

This repository contains one R script and four .CSV data files related to the Seed Bank Bias Project - 
an assessment of how the taxonomic and functional composition of vegetation and 
seedbank change across an edaphic gradient (i.e. soil terraces increasing in elevation 
and surface age) in Boulder, CO, USA. 

This environmental gradient consists of 12 sites located across six previously-described 
soil surfaces varying in age from <5000 years to >1-2 million years. 
In 2017-2018, plant community cover estimates were made (visual cover of species) 
and seed banks were sampled and grown out to estimate seed bank composition (count by species).

For a full description of the project and analyses, see the accepted manuscript (citation information above).
 

--------------------
DATA & FILE OVERVIEW
--------------------

Filenames and brief description of all data files:
1.  Larson&Suding_SeedbankBias_Rscript.R   
This file contains the R code used to run analyses associated with the accepted manuscript (see citation above)

2.  env_data.csv  
[Environmental matrix] This file contains data related to the edaphic gradient, including elevation, soil age rank, soil properties, and other features of each site (all properties quantified at the site-level).

3. traits_commonspp.csv    
[Trait matrix] This file contains mean trait values and characterstics for each species. Full species names and additional information can be found in the Supporting Information attached to the accepted manuscript (see citation above).

4. compiled_AVG_commonspp.csv
[Community matrix] This file contains species abundances (columns) in the vegetation (estimated as % aerial cover) and seed bank (estimated as counts) for each plot (rows). Values reflect averages across 2017 and 2018 sampling years. These are raw data that have not yet been relatived or transformed in any way. However, infrequent species (not found in >5% of plots in either the vegetation OR the seed bank) have already been removed in this version.

5. seedbank_veg_all_species_years.csv
[Supplemental community matrix]. Same as no. 4, but contains abundance data for ALL recorded species (including infrequent species) separated out by sampling year (2017 or 2018).


--------------------------
METHODOLOGICAL INFORMATION
--------------------------
Note that additional descriptions of variables within each data file can be found at the beginning of the R script (Larson&Suding_SeedbankBias_Rscript.R), and in the accepted publication (see citation information above).


