# SeedBankBias

Principal Investigator:  Julie Larson, University of Colorado Boulder, julie.e.larson@colorado.edu

Date of data collection   2017-2018.

Geographic location of data collection:  Boulder, CO, USA [See associated manuscript for more information on site locations]

R Software: All analyses were conducted in R (R version 3.5.1).

Funding sources or sponsorship that supported the collection of the data:  
This work was supported by the City of Boulder Open Space & Mountain Parks Funded Research Program, Univ. of Colorado Undergraduate Research Opportunities Program, and a graduate research grant from the Dept. of Ecology & Evolutionary Biology (Univ. of Colorado). J.L. was supported by a USDA NIFA Pre-Doctoral Fellowship (Accession No. 1019166)


--------------------------
SHARING/ACCESS INFORMATION
-------------------------- 

Licenses/restrictions placed on the data, or limitations of reuse: CC0 1.0 Universal (Public Domain, no limits to reuse) https://creativecommons.org/publicdomain/zero/1.0/legalcode

Contact with the PI regarding specific reuse of data is preferred.

Recommended citation for the data:  
This data is being archived in anticipation of manuscript publication in the journal: Ecology  
Search for: Larson, J.E & K.N. Suding. Seed bank bias: Differential tracking of functional traits in the seed bank and vegetation across a gradient. Ecology. 



--------------------
PROJECT OVERVIEW
--------------------

This repository contains one R script and four .CSV data files related to the Seed Bank Bias Project - 
an assessment of how the taxonomic and functional composition of the vegetation and seedbank change across 
an edaphic gradient (i.e. soil terraces increasing in elevation and surface age) in Boulder, CO, USA. 

The bulk of analyses explore patterns among three data matrices: environment, species composition, and traits. 


[Environment] The environmental gradient consists of 12 sites located across six previously-described soil surfaces varying 
in age from <5000 years to >1-2 million years. Older soil surfaces also occur on higher terraces, and ranked
terrace elevation is used as the main proxy for the gradient in analyses.  Soil variables and other properties 
are also included in the site-level environmental data matrix.

[Species] In 2017-2018, plant community cover estimates were made within several plots at each site (visual cover of species), 
and corresponding seed bank samples were grown out to estimate seed bank composition (count by species).

[Traits] We also compiled plant functional traits for common species in these communities, 
including leaf and root traits (collected in a greenhouse growout) and seed mass (collected / compiled). 


Additional descriptions of the data are provided below and at the beginning of the R script associated with this repository. 
For a full description of the project and analyses, see the accepted manuscript (citation information above).
 

--------------------
DATA & FILE OVERVIEW
--------------------

Filenames and brief description of all data files:
1.  Larson&Suding_SeedbankBias_Rscript.R   
This file contains the R code used to run analyses associated with the accepted manuscript (see citation above). 
It also contains metadata for each CSV file, including descriptions of all rows, columns, and units (as applicable) for each dataset used in analyses (see 'Data Types') within the R script.

2.  env_data.csv  
[Environmental matrix] This file contains data related to the edaphic gradient, including elevation, soil age rank, soil properties, and other features of each site (all properties quantified at the site-level).

3. traits_commonspp.csv    
[Trait matrix] This file contains mean trait values and characterstics for each species. Full species names and additional information can be found in the Supporting Information attached to the accepted manuscript (see citation above).

4. compiled_AVG_commonspp.csv
[Community matrix] This file contains species abundances (columns) in the vegetation (estimated as % aerial cover) and seed bank (estimated as counts) for each plot (rows). Values reflect averages across 2017 and 2018 sampling years. These are raw data that have NOT yet been relatived or transformed in any way. However, infrequent species (not found in >5% of plots in either the vegetation OR the seed bank) have already been removed in this version in preparation for most analyses.

5. seedbank_veg_all_species_years.csv
[Supplemental community matrix]. Same as no. 4, but contains abundance data for ALL recorded species (including infrequent species) separated out by sampling year (2017 or 2018). These data were only used by the authors for supplemental analyses of interannual differences and infrequent species removal effects.


---------
METADATA
---------
Information on the types and structure of data included in each .CSV file, including row and column descriptions, units, levels of sampling, etc., can be found in the R Script associated with repository (Larson&Suding_SeedbankBias_Rscript.R), in the section titled 'Data Types'.


--------------------------
METHODOLOGICAL INFORMATION
--------------------------
Detailed descriptions of the structure and variables for each dataset (i.e. CSV file) can be found in he R script (Larson&Suding_SeedbankBias_Rscript.R), and in the accepted publication.

Here, we include a shortened overview of the methods that is adapted from the accepted manuscript. 
This is meant to be a brief overview to guide interpretation of the data and data preparation/analyses carried out in the R Script.

For full details, search for:
J.E & K.N. Suding. Seed bank bias: Differential tracking of functional traits in the seed bank and vegetation across a gradient. Ecology.


Site and Sampling Design.

We established 12 sites across a 10km stretch of xeric grasslands in Boulder, CO, US.  Climate is semi-arid (488mm annual precip.), but grasslands contain a diverse mix of semi-arid shortgrass and mesic tallgrass species.  All sites are in natural areas, and most are grazed by cattle rotationally. At each site, we sampled vegetative composition, seed banks, and soils along a 50m transect.  At every 12.5m along each transect we installed two 1m2 plots (10 per site) for community composition data and some environmental metrics. 

Some datasets are summarized at the site-level (e.g., environmental data) while others occur at the plot-level (e.g., community composition).


Environmental gradient and data. 

Though the 12 sites were in close proximity (<10km), they occurred along a series of progressively older and higher soil terraces (i.e. terrace elevation gradient).  Elevation gain is relatively small (range: 1672m to 1920m), but estimated soil age increases from five thousand years (lower terraces) to two million years (higher terraces) across six pedologic surface types. This likely shapes local vegetation, possibly via impacts on soil water dynamics.
For site-level soil properties, we collected and pooled 9 replicate soil cores along each site’s 50m transect (2cm diam., 10cm depth [as possible]).  We sieved each sample (2mm), analyzing fresh soils for pH and dried soils for texture (Boyoucos Hydrometer method).  We then ground soils for total soil C and N (% by mass), and soil organic matter (loss-on-ignition method).  We also installed soil moisture sensors at 6 of the 12 sites across the gradient at depths of 10cm and 30cm to estimate shallow and deep soil volumetric water content (soil VWC) in May (wettest month; averaged here across 2017 and 2018).

We also quantified percent aerial cover of bare ground, litter (including detached litter and attached standing dead vegetation), cow manure, and rocks environmental variables in 7 out of 10 plots per site (subsampled due to time constraints), which were averaged for site-level estimates of these variables.  


Community composition.

We measured plant community composition as aerial cover of each species in 7 out of 10 plots per site (sub-sampled due to time constraints; averaged across 2017 & 2018 for each plot for most analyses). 

To characterize seed bank composition, we collected and pooled three seed bank samples (5cm diam x ~2cm depth) from just outside of each plot (n= 10 per site) in April 2017 and 2018.  Samples were grown out in the greenhouse in the same year of their collection (100mL of each pooled soil sample per plot spread on top of potting soil).  We counted emerging plants for 4 months, applied a 2 month dormancy-breaking cold period, then continued counting until emergence ceased (up to 4 months).  Plot-level seed bank counts were also averaged across 2017 & 2018 for most analyses. 


Trait sampling. 

We collected vegetative trait data for 54 species (a majority of vegetation and seed bank communities by abundance; see publication cited above).  For most species, we collected and averaged traits from 4 to 5 replicate plants grown in a common greenhouse environment (target age of 16 weeks, or as early as flowering occurred). Traits included: 
-Plant height (measured from soil to the tallest photosynthetic tissue, standardized by the number of -growing days to account for slight differences in final plant age)
-Root mass ratio (RMR, the ratio of total root dry mass relative to whole plant dry mass)
-Root dry matter content (RDMC, the ratio of dry mass to fresh mass in subsampled root tissue)
-Leaf dry matter content (LDMC, the ratio of dry mass to fresh mass in subsampled leaf tissue)
-Specific leaf area (SLA, leaf area per dry mass of subsampled leaves)
-Specific root length (SRL, fine root length per dry mass of subsampled roots) 
-Root diameter (Rdiam, the average diameter across the length of subsampled fine roots)
-Seed mass (the average dry mass per seed, estimated from the same seedlots when possible, or substituted from the Kew Seed Information Database). 


Analysis. 

Data preparation.  

For both seed bank and vegetative datasets, we pooled some morphologically-similar species at the genus or functional group level.  For analyses, we averaged all community cover and seed bank data across 2017-2018 at the plot-level (but the full, separated dataset is also included here). We also removed infrequent species from the combined vegetation and seed bank dataset (i.e. any species found in fewer than 5% of vegetation plots or seed bank samples; the full dataset is also included here).  
[Completed in R script:]  We remove four environmental variables sharing high correlations with others (soil N, soil organic matter, % sand, and cowpie cover) then check for linearity and normality among the reduced set of eight environmental variables (% clay, % silt, soil pH, soil C, % litter, % bare ground, % rock, % total vegetative cover) plus terrace elevation rank and soil age rank.  For abundance-weighted community analyses, we square root-transform vegetation and seed bank composition to moderate the influence of dominant species, then convert raw cover and counts to plot-level relative abundances before estimating site-level abundances (used for most analyses). We also log-transform seed mass and SRL to improve normality, then scale all traits prior to analyses.

Environmental Analysis.

We used constrained ordination (redundancy analysis, RDA) to characterize the portion of variation in environmental space (including soil age rank) that is explained by terrace elevation rank (i.e. the gradient proxy used in all subsequent analyses).  

Taxonomic seed bank bias. 

First, we use a linear model with site-level richness as a function of terrace elevation rank, community type (seed bank, vegetation, or seed bank plus vegetation species), and their interaction.  By including a combined community type (seed bank plus vegetation species), we aim to assess whether the seed bank stored unique species from the vegetation, increasing site-level richness.  We also use perMANOVA to explore how species composition was affected by community type (vegetation or seed bank), terrace elevation rank, and their interaction.  Finally, we use nonmetric multidimensional scaling (NMDS) to visualize species compositional differences between the seed bank and vegetation across the gradient.  

Functional seed bank bias and environmental tracking.  

To explore how seed banks contribute to functional diversity, we start with the trait by species matrix (8 traits) and three site by species matrices with relative abundances: vegetation only (79 species), seed bank only (74 species) and combined vegetation and seed bank (88 species).  For each species matrix, we estimate site-level functional diversity as functional dispersion, but also compare this to functional richness, Rao’s Q, and to estimates based on presence-absence.  We then model site-level functional dispersion as a function of community type (seed bank, vegetation, or seed bank plus vegetation), terrace elevation rank, and their interaction (linear models). 
To understand how the soil age gradient affects distributions of specific traits in the vegetation and seed bank, we also estimate trait community weighted means (CWMs) for each community type within each site based on species’ relative abundances in the vegetation and seed bank. We then model each trait CWM as a function of community type (vegetation or seed bank), terrace elevation rank, and their interaction (linear models). We also use fourth-corner tests to provide additional insight into the strength of trait-environment relationships within community types.  We conduct fourth-corner tests for the correlation between traits and terrace elevation rank within both vegetation and seed bank communities.  Because the function cannot handle incomplete trait data, we imputed missing values with a principal components analysis model.





