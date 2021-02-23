#' ---
#' title: "Chronosequence Analysis"
#' author: "Julie Larson"
#' date: "10 September 2020"
#' output: github_document
#' ---
#' 
#' 

#'  ############
#' 
#'  **Read in data**  
#'    



# Set user's Drive path
gdrive <- "/Users/Julie\ Larson/Google\ Drive"


#' Environmental data 
env_path <-  paste0(gdrive, "/Github/SeedBankBias/env_data.csv")
env <- read.csv(env_path)

# Full veg and sb data
veg_sb_path <-  paste0(gdrive, "/Github/SeedBankBias/compiled_AVG_commonspp.csv")
veg_sb <- read.csv(veg_sb_path)

# Trait data
trait_path <-  paste0(gdrive, "/Github/SeedBankBias/traits_commonspp.csv")
trait <- read.csv(trait_path)


#'  ############
#' 
#'  **Load packages**  
#'    

#+ results=FALSE, message=FALSE, warning=FALSE
library(tidyverse)
library(psych)
library(vegan)
library(geometry)
library(ade4)
library(FD)
library(plotly)
library(gridExtra)
library(Hmisc)
library(RColorBrewer)
library(grid)
library(ggcorrplot)
library(cowplot)
library(lme4)
library(lmerTest)
library(car)
library(questionr)
library(picante)
library(dichromat)
library(GGally)
library(patchwork)
library(ggcorrplot)


# Oridicenter() source code
source ('http://www.davidzeleny.net/anadat-r/doku.php/en:customized_functions:ordicenter?do=export_code&codeblock=0')

# Ordihull() source code
source ('http://www.davidzeleny.net/anadat-r/doku.php/en:customized_functions:orglhull?do=export_code&codeblock=1')

# Set option to not cut off tibbles
options(tibble.width = Inf)
options(tibble.print_max= Inf)  

# Set option for contrasts (Type III sums of squares (lm's))
options(contrasts = c("contr.sum","contr.poly"))

# Load function to calculate standard errors
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
} 





#' 
#'  **Overview:**
#'  This R script describes and runs analyses related to the Chronosequence Project - an assessment of how
#' the taxonomic and functional composition of vegetation and seedbanks change across an edaphic gradient related
#' to soil terraces increasing in eleveation and surface age in Boulder, CO, USA. 
#' This gradient consists of 12 sites across six previously-described soil surfaces
#' varying in age from <5000 years to >1-2 million years. Plant community cover estimates were made and seedbanks 
#' were sampled in 2017 and 2018.
#' 
#' 
#'  ############
#'   **Data types**
#'  ############ 
#'   

#'  ---
#'   **Environmental data** 
#'  ---
#'  
#'  *Data level* - Collected or aggregated at the site level (units and collection level indicated below in parentheses) 
#'  
#'  *Variables* -  
#'  elevation rank (lowest to highest, 1-12), 
#'  age rank (youngest to oldest, 1-6),   
#'  elevation (m, site-level),  
#'  clay (%, site-level, 2017),  
#'  silt (%, site-level, 2017),  
#'  sand (%, site-level, 2017),  
#'  pH (site-level, 2017),  
#'  soil C (% site-level, 2017),  
#'  soil OM (% site-level, from 2017),   
#'  litter (%cover, plot-level avg 2017-18),   
#'  bare ground (%cover, plot-level avg 2017-18),   
#'  cowpie (%cover, plot-level avg 2017-18),  
#'  rock (%cover, plot-level avg 2017-18),  
#'  total plant cover (%cover, plot-level, avg 2017-18),   
#'  soil N (%, site-level, 2017)  
#'  
#'  *Variables for qualitative assessment only* -  
#'  Soil VWC (May deep and shallow [among other possible metrics] at only 6/12 sites)
#'  


#'  ---
#'  **Vegetation and seedbank data**
#'  ---
#'  
#'  *Data Level* - Collected at plot level [seedbank has 10 plots and veg has a subset of 7 plots within each of 12 sites)]
#'  
#'  *Raw variables* -  
#'  Type ['veg' or 'sb']: Indicates whether data from each plot (row) is from the vegetation or seedbank   
#'  Site [A to L]: One of 12 letters is an indicator for site identity  
#'  Remaining columns: Community composition  
#'  
#'  Note on community comp: Percent plant cover (veg) or seedbank counts (sb) for each species or genera, 
#'  averaged across 2017-18 on a per plot/type basis. Note that the level of classification (e.g., species, 
#'  genera, or broader functional group) was determined as the coarsest level of identification across 
#'  veg and seedbank samples. For example, Heterotheca Villosa and H. foliosa were differentiable in the veg
#'  but not in the seedbank. Therefore, these two species are pooled at the genera level for both datasets
#'  for comparability.  
#'   
#'  *Rare species removal*-  
#'  The loaded dataset already has species removed that were rare in both seedbank and vegetation samples. Species
#'  were only retained if they were found in >5% of either seedbank samples or vegetation samples.  
#'    
#'  First, we started with each full dataset and removed any species not found in 5% of plots:  
#'     Seedbank: 94 known species/genera, unknown specimens, and pooled functional groups, reduced to 54 common  
#'     Vegetation: 110 known species/genera, unknown specimens, and pooled functional groups, reduced to 76 common  
#'    
#'  We then created a combined dataset that included any species retained in either the common sb or veg subset of species.  
#'  **This resulted in a final compiled dataset with 90 species common (>5% of plots) in either the seedbank or the veg**  
#'  
#'  *Data transformation*-  
#'  Seedbank data:  Prior to analysis, data will be squareroot-transformed and converted to relative abundances
#'  on a per plot basis. See 'seedbank_exploration_final.R' for reasoning behind this decision. Quick summary
#'  is that this improves evenness (necessary for trait-based inferences).  
#'  
#'  Veg data: Prior to analysis, data will be squareroot-transformed and converted to relative abundances
#'  on a per plot basis. See 'seedbank_exploration_final.R' for reasoning behind this decision. Quick summary
#'  is that this improves evenness (necessary for trait-based inferences).  
#'  

  
#'  ---
#'  **Trait data**
#'  ---
#' 
#' *Data Level* - Compiled at the species level (with up to 5 underlying reps per species).
#'     For genera (which generally only pool 2 or a few known possible species), we use trait data 
#'     from a single representative species. This information, and info on trait data availability
#'     across the full species/genera list, can be found in the raw data table.
#' 
#' *Raw Variables* -  
#'  Height (per day)  
#'  RMR  
#'  RDMC  
#'  SRL  
#'  Rdiam  
#'  LDMC  
#'  SLA  
#'  Seed mass  
#'    
#'  ############
#' 
#' 
#'  **Research Questions**
#'  
#'  **Gradient:** Does the soil chronosequence align with measured soil properties that could explain potential
#'  mechanisms of interaction with plants? (e.g., pH, soil texture, soil fertility)
#'  
#'  **Q1:** How do the taxonomic richness and composition of the vegetation and seedbank compare to one
#'  another and change across the soil chronsequence? 
#' 
#'  **Q2:** How do functional diversity and composition of the vegetation and seedbank compare to one
#'  another and change across the soil chronosequence?
#'  
#'  **Q3:** Are the functional responses of vegetation and seedbanks similar (i.e. the seedbank mirrors the veg, 
#'  and visa versa), or does the veg-seedbank discrepancy shift across the gradient? 
#' 
#' 
#'  ############
#' 






#' ############ 
#'  
#'  *Gradient:* Does the terrace elevation rank align with measured soil properties that could explain potential
#'  mechanisms of community turnover? 
#'  
#' ############  

#
#' **Data prep - Need to reduce environmental variables down to smaller subset** 
#'    Approach is to explore and remove:  
#'    a) highly correlated variables  
#'    b) any exploratory variables that appear non-functional (e.g., if not correlated with others
#'      and no definite ecological function  
#'    c) soil VWC, due to smaller sample size, after exploring correlations with other variables  


#' 
#' *Figure S1: Correlation heat map* showing degree of positive/negative correlations and significance  
#' 

# Env data structure
str(env)
env2 <- env %>%
  select(age_rank:vwc_may_avg_deep)
env2$elevation_rank <- as.numeric(rank(env2$elevation_m))
env2 <- env2 %>%
  select(elevation_rank, age_rank, elevation_m,  everything())

# Create correlation table and p-value table
env_corr <- round(cor(env2, use="pairwise.complete.obs"), 2)
p_mat <- cor_pmat(env2)

# Heat map for ALL correlations
#+ fig.width=8, fig.height=8
env_cor_fig_all <- ggcorrplot(env_corr, type="lower",method="circle", lab=T, lab_size = 2.5, tl.cex=10)
env_cor_fig_all

# Heat map for SIGNIFICANT correlations
env_cor_fig_sig <- ggcorrplot(env_corr, type="lower", method="circle", p.mat=p_mat, insig="blank", lab=T, lab_size = 2.5, tl.cex=10)
env_cor_fig_sig 

# Save full correlation heat map
tiff(filename="env_corr_w_elev.tiff", res=600, width=7, height = 5, units = "in")
env_cor_fig_all
dev.off()



#' 
#' *Env Correlation Assessment:* Note variables that meet conditions for removal  
#' 

# a) Highly correlated variables (r>0.8)
#     - Elevation is highly correlated with soil age rank (r=0.94), due to the terrace
#       formation involved in the different soil surface ages. Unless we think elevation 
#       has some independent important effects (e.g., on temp, precip, etc.), will 
#       probably remove elevation


#' 
#' Terrace elevation and soil age rank
corr.test(env$elevation_m,env$age_rank)
cor_elev <- ggplot( aes(y=elevation_m, x=age_rank), data=env) + 
  geom_point(alpha=0.4, cex=3) +
  #geom_smooth(method="loess", se=F) +
  geom_smooth(method="lm", se=F, lty=2, col="gray40")+
  geom_text(aes(x=1.25, y=1900, label="A) \n          r = 0.94"), col="gray40") +
  labs(x="Estimated Soil Age Rank", y="Terrace Elevation (m)")+
  scale_x_continuous(breaks=1:6) +
  #ylim(1500,2000)+
  theme(text = element_text(size=12)) +
  theme_bw()
cor_elev


#' Terrace Elevation RANK and  soil age rank
env$elev_rank <- as.numeric(rank(env$elevation_m)) 

corr.test(env$elev_rank,env$age_rank)
cor_elevrank <- ggplot( aes(y=elev_rank, x=age_rank), data=env) + 
  geom_point(alpha=0.4, cex=3) +
  #geom_smooth(method="loess", se=F) +
  geom_smooth(method="lm", se=F, lty=2, col="gray40")+
  geom_text(aes(x=1.25, y=11, label="B) \n          r = 0.93"), col="gray40") +
  labs(x="Estimated Soil Age Rank", y="Terrace Elevation (rank)")+
  scale_x_continuous(breaks=1:6) + 
  scale_y_continuous(breaks=seq(1,12,2)) +
  theme(text = element_text(size=12))+
  theme_bw()
cor_elevrank


#' Correlation panel of the three main traits

oldnames = c("age_rank","elevation_m","elev_rank")
newnames = c('Soil Age Rank','Elevation (m)', 'Elevation (rank)')

cor_elev <- env %>% 
  select(age_rank, elevation_m, elev_rank) %>% 
  rename_at(vars(oldnames), ~ newnames)%>% 
  ggpairs(., lower = list(continuous = wrap("smooth", alpha = 0.5, size=1)))
cor_elev



#     - Soil N and Soil C are highly correlated (r=0.95). Because soil N is missing for
#       one sample, we will remove soil N and retain soil C.

cor2 <- ggplot( aes(y=soilC, x=soilN), data=env) + 
  geom_point() +
  geom_smooth(method="lm", se=F)
cor2

#     - Soil OM and Soil C are highly correlated (r=0.85). Soil C has one more sig. 
#       correlation with other variables, so we will remove soil OM and retain soil C.

cor3 <-ggplot( aes(y=soilC, x=om_2017), data=env) + 
  geom_point() +
  geom_smooth(method="lm", se=F)
cor3

#     - Sand is significantly correlated with clay (r=-0.81) and less so with silt(r=-0.72).
#       However, neither clay nor silt are correlated with eachother. Therefore, we will 
#       retain clay and silt, but remove sand.

cor4 <- ggplot( data=env) + 
  geom_point(aes(x=sand, y=clay), col="green") +
  geom_smooth(aes(x=sand, y=clay), method="lm", se=F, col="green") +
  geom_point(aes(x=sand, y=silt), col="blue") +
  geom_smooth(aes(x=sand, y=silt), method="lm", se=F, col="blue") +
  labs(y="% Clay (green) or % Silt (blue)", x="% Sand")
cor4

cor5 <- ggplot( aes(y=silt, x=clay), data=env) + 
  geom_point() +
  geom_smooth(method="lm", se=F)
cor5

# b) any exploratory variables that appear non-functional (e.g., non-correlated with others)  
#     - Cowpie cover only varies from 0 to 1.8%, and it's not clear that it is
#       accurately describing current grazing pressure. If that were the case, we might expect
#       a negative association with litter, plant cover, or a positive association with bare ground
#       and all we see is a positive association with litter (r=0.7); maybe cows are breaking up veg into
#       litter but not hitting hard enough to break it down/remove? Litter does have some weaker
#       correlations with other variables, so we will retain litter here, but remove cowpie.

cor6 <- ggplot( aes(y=litter_cov, x=cowpie_cov), data=env) + 
  geom_point() +
  geom_smooth(method="lm", se=F)
cor6

# Create grid of all plots
high_cor <- grid.arrange(cor1, cor2, cor3, cor4, cor5, cor6, ncol=2)
high_cor

# c) soil VWC, due to smaller sample size, after exploring correlations with other variables  
#     - Shallow May soil VWC is tightly correlated (r=-0.89) with total vegetation cover
#     - Deep May soil VWC is tightly tied to bare ground (r=0.91)

# Create reduced dataset with only variables of interest
vwc_dat <- env %>%
  select(vwc_may_avg_shallow, vwc_may_avg_deep, bare_cov, veg_cov) %>%
  gather(key="depth", value="vwc", vwc_may_avg_shallow, vwc_may_avg_deep)
vwc_dat$depth <-  factor(vwc_dat$depth,
                            levels = c("vwc_may_avg_shallow", "vwc_may_avg_deep"),
                            labels = c("shallow","deep"))
 
# Create bare ground figure
vwc_bare <- ggplot( aes(y=100*vwc, x=bare_cov), data=vwc_dat) + 
  geom_point(aes(col=depth),cex=2.5) +
  geom_smooth(method="lm", se=F, data=subset(vwc_dat, depth=="deep"),col="dodgerblue4") +
  scale_color_manual(values=c("dodgerblue","dodgerblue4")) +
  labs(x="Bare ground (%, site average)", y="Soil volumetric water content (% in May)")+
  annotate(geom="text", x=24, y=36, label="r=0.91", color="dodgerblue4", fontface="bold") + 
  theme(legend.position = "none")
vwc_bare

# Create veg cover figure
vwc_veg <- ggplot( aes(y=100*vwc, x=veg_cov), data=vwc_dat) + 
  geom_point(aes(col=depth), cex=2.5) +
  geom_smooth(method="lm", se=F, data=subset(vwc_dat, depth=="shallow"),col="dodgerblue") + 
  scale_color_manual(values=c("dodgerblue","dodgerblue4")) +
  labs(x="Vegetative cover (%, site average)", y="Soil volumetric water content (% in May)", col="Soil Depth")+
  annotate(geom="text", x=87, y=18, label="r=-0.89", color="dodgerblue", fontface="bold") +
  xlim(40,100)

# View and save figures
grid.arrange(vwc_bare, vwc_veg, ncol=2, widths=c(1,1.4))

tiff(filename="vwc_corr.tiff", res=600, width=9, height = 4.5, units = "in")
grid.arrange(vwc_bare, vwc_veg, ncol=2, widths=c(1,1.4))
dev.off()





###########
#'  
#'  *RDA*
#'  



#'  *Elevation rank* 
#'  #########
#'  
#'  **Test 1B) RDA with elevation rank as constraining variable **  
#'  
#'

#'  
#'  **Create reduced ENV dataset** with redundant and non-functional variables removed  
#'  

# Remove correlated variables
str(env2)
env_red <- env2 %>% 
  select(-elevation_m, -elevation_rank,  -soilN, -om_2017, -sand, -cowpie_cov, -vwc_may_avg_shallow, -vwc_may_avg_deep)

# View retained variables for linearity/distributions
pairs(env_red)

# One high outlier in tot_cov -- log transforming helps some. 
hist(env_red$veg_cov, breaks=20)
env_red$veg_cov <- log(env_red$veg_cov)
hist(env_red$veg_cov, breaks=20)
pairs(env_red)

# Scale all
env_red_s <- scale(env_red) #scaled

# Save elevation as a separate variable
age <- env %>%
  select(age_rank)
elev_rank <- env2 %>%
  select(elevation_rank)

#'   
#'   *Model 1b:*: RDA (env ~ elev_rank)
#'   

#'   *Assessment 1b*: Total RDA R2 for age rank = 0. (env variation constrained by elevation age rank)
env_rda <- rda(env_red_s~., elev_rank, scale=FALSE) 
summary(env_rda)

#' 
#' *Figure 1b:* RDA1 vs PC1

# View rda() version of plot
plot(env_rda)

# Save data for figure
rda_plot <-plot(env_rda)

rda_scores <- data.frame(rda_plot$sites)  #Site scores
rda_scores$elev_rank <- as.numeric(env2$elevation_rank)
rownames(rda_scores) <- env$site
rda_scores

rda_vecs <- data.frame(rda_plot$species)  #Env vectors
rownames(rda_vecs) <- c("Soil age rank","% Clay", "% Silt", "pH", "soil C", "% Litter", "% Bare", "% Rock", "% Veg Cover")
rda_vecs  
rda_vecs$xlab <- c(1.1,-0.22, 0.57, -0.86, 0.84, .27, 0.75, 0.29, -0.68)
rda_vecs$ylab <- c(.1,.67, 0.71, 0.5, 0.63, 0.36, -0.20, -1.1, 0.63)

rda_plot$biplot# Elevation Loading
rda_elev <- data.frame(cbind(1.329, 0)) 
rownames(rda_elev) <- "Elevation rank"
colnames(rda_elev) <- c("RDA1", "PC1")
rda_elev

rda_vwc <- cbind(rda_scores, env$vwc_may_avg_shallow, env$vwc_may_avg_deep) # Soil VWC correlation
colnames(rda_vwc)[4] <- "shallow vwc"
colnames(rda_vwc)[5] <- "deep vwc"
rda_cor <- round(cor(rda_vwc, use="pairwise.complete.obs"), 2)
rda_vwc <- data.frame(rda_cor[c(4:5),c(1:2)])
rda_vwc


#' *Figure 1b:* Biplot= RDA1 vs PC1 with VWC overlaid 
rda_fig2 <- ggplot () +
  geom_point(aes(x=RDA1, y=PC1, color=elev_rank), cex=3, data=rda_scores)+
  scale_color_gradient(low="#BDCEDF", high="navy") +
  geom_segment(mapping=aes(x=0, y=0, xend=RDA1, yend=PC1), size=1.2, data=rda_vecs) +
  geom_segment(mapping=aes(x=0,y=0, xend=RDA1, yend=PC1), arrow=arrow(length = unit(0.4,"cm")), size=1, col="red", data=rda_elev) +
  geom_segment(mapping=aes(x=0,y=0, xend=RDA1, yend=PC1), arrow=arrow(length = unit(0.4,"cm")), size=1, col="blue", data=rda_vwc) +
  labs(x="RDA Axis 1 - 31% variance", y="PC Axis 1 - 26% variance",col="Elevation Rank") +
  geom_text(mapping=aes(x=xlab, y=ylab, label = rownames(rda_vecs)), data=rda_vecs) +
  geom_text(aes(x=RDA1, y= PC1-0.13, label="Elevation \n Rank"), col="red", data=rda_age) +
  geom_text(aes(x=RDA1+0.36, y= PC1, label= rownames(rda_vwc)), col="blue", data=rda_vwc)+
  xlim(-2.25,1.5)+
  theme_classic()
rda_fig2

tiff(filename="soilage_rda_w_elev_rank.tiff", res=600, width=6, height = 5.5, units = "in")
rda_fig2
dev.off()





#'  ############
#'  
#'  **Q:** How do the taxonomic richness and composition of the vegetation and seedbank compare to one
#'  another and change across the soil chronsequence? 
#'  [i.e. Is there a pattern of community response with implications for management?]
#'  
#'  ############



#' 
#' *Preliminary Test 1a: Vegetation richness as a function of soil age rank*  
#' 
#' Note, rare species have already been removed to include only plants found in >5% of veg or sb samples

# Filter veg data
veg <- veg_sb %>% 
  filter(type == "veg")

# Compile veg by site
veg_site <- veg %>%
  select(-type, -age_rank) %>%
  group_by(site) %>%
  summarise_all(mean)
  
# Reduce site cover to presence/absence
veg_site_pa <- veg_site %>%
  replace(((.)>0),1)

# Sum across columns to generate a total number of species (total richness) per site and year
veg_site_pa %>%
  ungroup() %>%
  select(acevul:virfal) %>% 
  rowSums(na.rm=TRUE) ->
  env$veg_rich


#' 
#' *Model 1a:* lm(veg_rich ~ elevation rank)  
#' *Assessment 1a (elev rank):* 
#' 
env$elev_rank <- env2$elevation_rank

lm_veg_rich2 <- lm(veg_rich ~ elev_rank, data=env)
summary(lm_veg_rich2)


#'  
#'  *Supplemental Info 1a:* Veg taxonomic evenness across sites  
#' 

# Sqrt transform veg data (veg_sqrt)
veg_sqrt <- veg %>%
  select(-type, -site, -age_rank)
veg_sqrt <- sqrt(veg_sqrt)

# Save sqrt-transformed data with site variable (veg_sqrt1)
veg_sqrt1 <- cbind(veg[,2],veg_sqrt)
colnames(veg_sqrt1)[1] <- "site"

# Convert to relative abundances
veg_sqrt_rel <- decostand(veg_sqrt, "total")
veg_sqrt_rel1 <- cbind(cbind(veg[,2], veg_sqrt_rel))
colnames(veg_sqrt_rel1)[1] <- "site"

# Calculate site mean relative abundances
site_avg_veg <- veg_sqrt_rel1 %>%
  group_by(site) %>%
  summarise_all(mean) %>%
  select(-site)

# Site avg evenness
env$veg_even <- diversity(site_avg_veg)/log(env$veg_rich)  # Calculate Pielou's evenness



#' 
#' *Preliminary Test 1b: Seedbank richness as a function of elevation rank*  
#' 
#' Note, rare species have already been removed to include only speices found in >5% of veg or sb samples

# Filter sb data
sb <- veg_sb %>% 
  filter(type == "sb")

# Compile sb by site
sb_site <- sb %>%
  select(-type, -age_rank) %>%
  group_by(site) %>%
  summarise_all(mean)

# Sum across columns of full counts to generate a total seedbank size per site
sb_site %>%
  ungroup() %>%
  select(acevul:virfal) %>% 
  rowSums(na.rm=TRUE) -> 
  env$sb_size
env$sb_size <- 200*env$sb_size  #Mulitply times 200 to go from seeds per 100cm3 to seeds per 20000 cm3

# Reduce site cover to presence/absence
sb_site_pa <- sb_site %>%
  replace(((.)>0),1)

# Sum across columns to generate a total number of species (total richness) per site 
sb_site_pa %>%
  ungroup() %>%
  select(acevul:virfal) %>% 
  rowSums(na.rm=TRUE) ->
  env$sb_rich


#' 
#' *Supp Model 1c:* lm(sb_rich ~ elev rank)  
#' *Assessment 1c (elev rank):* 
#' 

lm_sb_rich <- lm(sb_rich ~ elev_rank, data=env)
summary(lm_sb_rich)



#' 
#' *Supplemental Info 2d*: SB size across sites and as a function of elevation rank
#' 

lm_sb_size2 <- lm(sb_size ~ elev_rank, data=env)
summary(lm_sb_size2)

ggplot(data=env) +
  geom_point(aes(x=elev_rank, y=sb_size)) +
  labs(x="Elevation", y="Seedbank size (per m2)")


#'  
#'  *Supplemental Info 2e:* Estimate SB taxonomic evenness across sites  
#' 

# Sqrt transform sb data (sb_sqrt)
sb_sqrt <- sb %>%
  select(-type, -site, -age_rank)
sb_sqrt <- sqrt(sb_sqrt)

# Save sqrt-transformed data with site variable (veg_sqrt1)
sb_sqrt1 <- cbind(sb[,2],sb_sqrt)
colnames(sb_sqrt1)[1] <- "site"

# Convert to relative abundances
sb_sqrt_rel <- decostand(sb_sqrt, "total")
sb_sqrt_rel1 <- cbind(cbind(sb[,2], sb_sqrt_rel))
colnames(sb_sqrt_rel1)[1] <- "site"

# Calculate site mean relative abundances
site_avg_sb <- sb_sqrt_rel1 %>%
  group_by(site) %>%
  summarise_all(mean) %>%
  select(-site)

# Site avg evenness
env$sb_even <- diversity(site_avg_sb)/log(env$sb_rich)  # Calculate Pielou's evenness


#'
#'    **Supplemental Assesssment: Taxonomic evenness as a function of elevation rank and community type**
#'    

#' *FIGURE* 

#' Elevation rank figure
even_fig_dat <- env %>% 
  select(veg_even, sb_even) %>%
  gather(key="type", value="evenness")
even_fig_dat$age_rank <- rep(1:6,each=2)
even_fig_dat$elev <- env$elevation_m
even_fig_dat$elev_rank <- rep(as.numeric(rank(env$elevation_m)),2)
even_fig_dat$type <-  factor(even_fig_dat$type,
                             levels = c("veg_even","sb_even"),
                             labels = c("Vegetation","Seedbank"))

even_fig2 <- ggplot(aes(x=elev_rank, y=evenness, col=type), data=even_fig_dat) + 
  geom_point() + 
  geom_smooth(method="lm",se=F) +
  scale_x_continuous(breaks = seq(1,12,2)) +
  scale_color_manual(values=c("black","#2171b5")) + 
  labs(x="Terrace Elevation (rank)", y="Pielou's Evenness", col="Community type") +
  ylim(0,1)
even_fig2


tiff(filename="evenness_fig_elev_rank.tiff", res=600, width=6, height=4, units = "in")
even_fig2
dev.off()

#'
#' *Model*
#' 

# Try basic linear model, note that there is moderate significance withOUT an interaction, but not with 
even_lm_rank <- lm(evenness ~ type * elev_rank, data=even_fig_dat)
summary(even_lm_rank)
Anova(even_lm_rank, type="III")
Anova(even_lm_rank, type="II")

even_lm2_rank <- lm(evenness ~ type + elev_rank, data=even_fig_dat)
summary(even_lm2_rank)
Anova(even_lm2_rank, type="III")
Anova(even_lm2_rank, type="II")



#'  
#'  **Taxonomic richness as a function of elevation rank and community type**
#'  

#' 
#' * Prep*
#' 
#' We already have separate veg and sb richness, but need a combined richness estimate

# Compile veg/sb by site
veg_sb_site <- veg_sb %>%
  select(-type) %>%
  group_by(site) %>%
  summarise_all(mean)

# Reduce site cover to presence/absence
veg_sb_site_pa <- veg_sb_site %>%
  replace(((.)>0),1)

# Sum across columns to generate a total number of species (total richness) per site and year
veg_sb_site_pa %>%
  ungroup() %>%
  select(acevul:virfal) %>% 
  rowSums(na.rm=TRUE) ->
  env$veg_sb_rich

#' Create figure
# Select and reconfigure data
rich_fig_dat <- env %>%
  select(veg_rich, sb_rich, veg_sb_rich) %>%
  gather(key="type", value="richness")
rich_fig_dat$age_rank <- rep(1:6, each=2)
rich_fig_dat$elevation <- env$elevation_m
rich_fig_dat$elev_rank <- rep(as.numeric(rank(env$elevation_m)),times=3)
rich_fig_dat <- data.frame(rich_fig_dat)
rich_fig_dat$type <-  factor(rich_fig_dat$type,
         levels = c("veg_rich","sb_rich","veg_sb_rich"),
         labels = c("Vegetation","Seedbank","Veg+Seedbank"))

#' Summary of values
veg_sb_compare <- env %>% 
  select(veg_rich:veg_sb_rich) %>%
  summarise_all(mean)
veg_sb_compare

#'  Average richness at the site-level is higher in the veg (38.8 species) than the 
#'  seedbank (27.9 species). When the vegetation and seedbank at a site are considered 
#'  together, richness increases to 49.6 species on average, suggesting that the seedbank 
#'  contains some unique species and turnover (not simply nested from veg) 


#' 
#' *Figure*
#' 
rich_fig <- ggplot (aes(x=elev_rank, y=richness, col=type), data=rich_fig_dat) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  labs(x= "Elevation (m)", y="Species richness", col="Community type") + 
  scale_color_manual(values=c("#636363","#9ECAE1","black")) 
rich_fig

tiff(filename="richness_fig_elev.tiff", res=600, width=6, height = 4, units = "in")
rich_fig
dev.off()


#'  
#'  *Model: Species richness as a function of type (seedbank or veg) and elevation rank*  
#'  


#' ELEVATION RANK
# Run linear model to report
rich_fig_dat$type <- factor(rich_fig_dat$type, levels = c("Veg+Seedbank",  "Seedbank",  "Vegetation"))
rich_lm_rank <- lm(richness ~ type * elev_rank, data=rich_fig_dat)
summary(rich_lm_rank)
Anova(rich_lm_rank, type="III")
Anova(rich_lm_rank, type="II")

# Removing interaction does not change results qualitatively at all!
rich_lm_rank2 <- lm(richness ~ type + elev_rank, data=rich_fig_dat)
summary(rich_lm_rank2)
Anova(rich_lm_rank2, type="III")
Anova(rich_lm_rank2, type="II")



#'  
#'  **PerMANOVA and NMDS:  Species community composition as a function of type (seedbank or veg) and age rank**  
#'  



#'  *PerMANOVA*
#'  permanova(veg_sb_matrix ~ type * age_rank, strata=site, data=veg_sb)  
#'      [NOTE that in both plot- and site-level analyses here]

#'
#'  Set up analysis to include *all plots* as replicates (not accounting for pseudo-replication within sites)
#'  
veg_sb_sqrt <- veg_sb %>%
  select(-type,-site,-age_rank)
veg_sb_sqrt <- sqrt(veg_sb_sqrt)
veg_sb_rel <- decostand(veg_sb_sqrt, "total")

#' Create dataframe of explanatory factors, including elevation
factors <- veg_sb %>%
  select(type, age_rank,site)
levels(factors$site)
factors$site <- fct_recode(factors$site, 
           "1680" = "A",
           "1672" = "B",
           "1677" = "C",
           "1687" = "D",
           "1689" = "E",
           "1704" = "F",
           "1806" = "G",
           "1809" = "H",
           "1908" = "I",
           "1920" = "J",
           "1907" = "K",
           "1909" = "L")
names(factors)[names(factors) == "site"] <- "elev"
factors$elev <- as.character(factors$elev)
factors$elev <- as.numeric(factors$elev)
factors$site <- veg_sb$site
factors$site <- fct_recode(factors$site, 
           "3" = "A",
           "1" = "B",
           "2" = "C",
           "4" = "D",
           "5" = "E",
           "6" = "F",
           "7" = "G",
           "8" = "H",
           "10" = "I",
           "12" = "J",
           "9" = "K",
           "11" = "L")
names(factors)[names(factors) == "site"] <- "elev_rank"
factors$elev_rank <- as.character(factors$elev_rank)
factors$elev_rank <- as.numeric(factors$elev_rank)
str(factors)


#' Run Plot-level PerMANOVA
# Elevation rank, NO blocking factors
spp_perm_plot2 <-adonis(veg_sb_rel ~ type * elev_rank, data=factors, permutations = 999, method="bray")
spp_perm_plot2


#' 
#' Set up analysis to pool at the *site level*, no replication of plots within sites
#' 
site_veg_sb_rel <- rbind(site_avg_veg, site_avg_sb)
type <- as.factor(rep(c("veg","sb"),each=12))
age_rank <- as.numeric(rep(1:6, each=2,2))
factors_site <- data.frame(type, age_rank)
factors_site$site <- env$site
factors_site$elev <- env$elevation_m
factors_site$elev_rank <- env$elev_rank
str(factors_site)


#' 
#' Run Site-level PerMANOVA
#' 

# ELEVATION rank
spp_perm_site2 <-adonis(site_veg_sb_rel ~ type * elev_rank, data=factors_site, permutations = 999, method="bray")
spp_perm_site2


#' 
#' Set up analysis to pool at the site level, AND include only *presence/absence data*
#' 
site_veg_sb_rel <- rbind(site_avg_veg, site_avg_sb)
site_veg_sb_pa <- site_veg_sb_rel %>%
  replace(((.)>0),1)
type <- as.factor(rep(c("veg","sb"),each=12))
age_rank <- as.numeric(rep(1:6, each=2,2))
factors_site <- data.frame(type, age_rank)
factors_site$site <- env$site
factors_site$elev <- env$elevation_m
factors_site$elev_rank <- env$elev_rank


#' 
#' Run Site-leve, presences-absence PerMANOVA
#' 

# ELEVATION rank
spp_perm_site_pa2 <-adonis(site_veg_sb_pa ~ type * elev_rank, data=factors_site, permutations = 999, method="bray")
spp_perm_site_pa2




#'   
#' *Nonmetric multidimensional scaling (NMDS)*
#'
#'  

#' *Site-level NMDS*
site_nmds <- metaMDS(site_veg_sb_rel,k=2,trymax=100)

site_nmds
stressplot(site_nmds)


#' *Create figure*
# Saving NMDS scores and labels
site_species <- site_nmds$species
site_scores <- cbind(factors_site, site_nmds$points)
site_scores$type_age <- paste(site_scores$type, site_scores$age_rank)
site_scores$site_names <- as.character(env$site_names)
site_scores$type <- dplyr::recode(site_scores$type, sb="S", veg="V")
site_scores$site_names2 <- as.character(c("1a","1b","2a","2b","3a","3b","4a","4b","5a","5b","6a","6b"))
site_scores$site_names3 <- paste(site_scores$type, site_scores$site_names2, sep = "")
site_scores$elev_rank <- env2$elevation_rank
site_scores$type_elev <- paste(site_scores$type, site_scores$elev_rank, sep="")
site_scores

# Generate a vector of colors for Elevations
# For plot-level
colfunc_sb <- colorRampPalette(c("white", "dodgerblue4"))
sb_col <- colfunc_sb(12)
colfunc_veg <- colorRampPalette(c("white", "gray16"))
veg_col <- colfunc_veg(12)
nmds_col <- data.frame(cbind(sb_col, veg_col))
rank(env$elevation_m) 
# Create vector that moves colors to appropriate place based on elevation rank
nmds_col$elev_rank <- c(2,3,1,4,5,6,7,8,11,9, 12,10)
nmds_col_list <- nmds_col %>% 
  arrange(elev_rank) %>% 
  select(-elev_rank) %>% 
  gather(key="type", value="col") #%>%
  #arrange(desc(type))
nmds_col_list$col
elev_col <- as.character(nmds_col_list$col)

# For site-level
colfunc_sb2 <- colorRampPalette(c("skyblue1", "navy"))
sb_col2 <- colfunc_sb2(12)
colfunc_veg2 <- colorRampPalette(c("gray76", "black"))
veg_col2 <- colfunc_veg2(12)
nmds_col2 <- data.frame(cbind(sb_col2, veg_col2))
# Create vector that moves colors to appropriate place based on elevation rank
nmds_col2$elev_rank <- c(2,3,1,4,5,6,7,8,11,9, 12,10)
nmds_col_list2 <- nmds_col2 %>% 
  arrange(elev_rank) %>% 
  select(-elev_rank) %>% 
  gather(key="type", value="col") %>%
  arrange(desc(type))
nmds_col_list2$col
elev_col2 <- as.character(nmds_col_list2$col)


# Saving color vectors
colvecage_type <- c("#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c","#d9d9d9", "#bdbdbd", "#969696", "#737373", "#525252", "#252525")
colvecage_type2 <- rep(c("#d9d9d9", "#bdbdbd", "#969696", "#737373", "#525252", "#252525", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c"), each=2)
colvecage_type3 <- rep(c("#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c","#d9d9d9", "#bdbdbd", "#969696", "#737373", "#525252", "#252525"), each=2)
colvecage_type4 <- rep(c("#9ecae1", "#08519c", "#c6dbef", "#6baed6", "#4292c6", "#2171b5", "#bdbdbd", "#525252", "#d9d9d9", "#969696", "#737373",  "#252525"), each=2)
colvecage_type4 <- rep(c("#9ecae1", "#08519c", "#c6dbef", "#6baed6", "#4292c6", "#2171b5", "#bdbdbd", "#525252", "#d9d9d9", "#969696", "#737373",  "#252525"), each=2)
colvec_elev <- c("#08519c","#9ecae1","#9ecae1","#08519c", "#c6dbef", "#c6dbef", "#6baed6", "#6baed6", "#4292c6", "#2171b5" ,"#4292c6", "#2171b5", "#525252","#bdbdbd" ,"#bdbdbd",  "#525252" ,"#d9d9d9", "#d9d9d9", "#969696" ,"#969696" ,"#737373","#252525",  "#737373" ,"#252525")
colvecage <- c("#FED976", "#FEB24C", "#FD8D3D", "#F03B20", "#BD0026","#800026")
colvecage2 <- rep(c("#FED976", "#FEB24C", "#FD8D3D", "#F03B20", "#BD0026","#800026"),each=2,2)
colvectype <- c("#2171b5","#525252")
colvectype <- c("blue","black")


tiff(filename="site_nmds_elev_spp.tiff", res=600, width=6, height = 5, units = "in")

# Plotting NMDS
par(mar=c(4,4,2,2))
ordiplot(site_nmds, display="si", type="n")
orditorp(site_nmds, display="species", col="gray")
with(factors_site, ordiellipse(site_nmds, type, draw="polygon", col=colvectype, alpha=0.1,lwd=0.5, kind="se", conf=0.95), label=TRUE)
orditorp (site_nmds, display="sites", label=site_scores$type_elev, col=elev_col2, air=0.01, cex=0.75)
with(site_scores, ordicenter(site_nmds, type, col="black", cex=1.2))
#with(site_scores, ordiellipse(site_nmds, type_age, col=colvecage_type4, lwd=0.5, draw="polygon", kind="se", conf=0.95),alpha=0.5)
#with(factors_site, ordiellipse(site_nmds, type, col="black", kind="se", conf=0.95, lwd=0.1, label=TRUE))
#with(site_scores, ordicenter(site_nmds, age_rank, col="darkgreen", cex=0.75))
#orditorp (site_nmds, display="sites", label=site_scores$type_age, col=colvecage_type2, air=0.01, cex=0.75)
#with(site_scores, ordicenter(site_nmds, type_age, col=colvecage2))
#with (factors_site, ordihull(site_nmds, age_rank, draw="polygon", col=colvecage, alpha=0.2, lwd=0.11, label=TRUE))

dev.off()


#' 
#' *Plot-level NMDS*
#' 
plot_nmds <- metaMDS(veg_sb_rel,k=2,trymax=100)

plot_nmds
stressplot(plot_nmds)

#' *Create figure*
# Saving NMDS scores and labels
plot_species <- plot_nmds$species
plot_scores <- cbind(factors, plot_nmds$points)
plot_scores$type <- dplyr::recode(plot_scores$type, sb="S", veg="V")
plot_scores$type_age <- paste(plot_scores$type, plot_scores$age_rank, sep = "-")
plot_scores$site <- veg_sb$site
plot_scores$site2 <- dplyr::recode(plot_scores$site, A="PP1", B="PP2", C="L1", D="L2", E="S1", F="S2", G="V1", H="V2", I="YRF1", J="YRF2", K="ORF1", L="ORF2")
plot_scores$sitesubcode <- dplyr::recode(plot_scores$site, A="a", B="b", C="a", D="b", E="a", F="b", G="a", H="b", I="a", J="b", K="a", L="b")
plot_scores$type_site <- paste(plot_scores$type, plot_scores$site2, sep = "-")
plot_scores$type_site2 <- as.factor(paste(plot_scores$type, plot_scores$age_rank, plot_scores$sitesubcode, sep = ""))
plot_scores$site_rank <- dplyr::recode(plot_scores$site, A="3", B="1", C="2", D="4", E="5", F="6", G="7", H="8", I="10", J="12", K="9", L="11")
plot_scores$type_elev <- paste(plot_scores$type, plot_scores$site_rank, sep = "")
head(plot_scores)
str(plot_scores)
levels(plot_scores$type_site2)

tiff(filename="plot_nmds_by_elev_wSPP.tiff", res=600, width=6, height = 6, units = "in")

# Plotting NMDS with ELEVATION color gradient
par(mar=c(4,4,2,2))
ordiplot(plot_nmds, type="n", ylim=c(-1.5,1.5))
orditorp (plot_nmds, display="species", col="darkgray", air=0.01, cex=0.6)
#with(plot_scores, ordiellipse(plot_nmds, type_site, col=colvecage_type4, lwd=0.5, draw="polygon", kind="se", conf=0.95),alpha=0.5)
with(plot_scores, ordiellipse(plot_nmds, type_site2, col=elev_col, lwd=0.5, draw="polygon", kind="se", conf=0.95),alpha=0.5)
with(plot_scores, ordicenter(plot_nmds, type_elev, col="black", cex=0.7))
#with(factors, ordiellipse(plot_nmds, type, draw="polygon", col=colvectype, alpha=0.3,lwd=0.1, kind="se", conf=0.95), label=TRUE)
#with (factors, ordiellipse(plot_nmds, type, col=colvectype, kind="se", conf=0.95, lwd=0, label=TRUE))
#with(plot_scores, ordiellipse(plot_nmds, type_age, col=colvecage_type, draw="polygon", kind="se", conf=0.99),alpha=0.5)
#with(plot_scores, ordicenter(plot_nmds, type_age, col="black", cex=0.72))

dev.off()


  
  
#'  
#'  **Q2: How do functional diversity and composition of the vegetation and seedbank compare to one**
#'  **another and change across the soil chronosequence?**  
#'  [i.e. What are the possible mechanisms of community response and implications for future managment?)  
#'  

#'  
#'  *Data prep*
#'  

# Grab a subset of columns that contain mature trait data

#' **HIGHLIGHT TO RUN**
trait_dat <- trait %>%
  select(species,height_per_day:seed_mass)
trait_dat <- trait_dat[order(trait_dat$species),]
str(trait_dat)

pairs(trait_dat[2:9])

#' 
#' *Estimate trait coverage in the vegetation*  
#' 

#' Grab square-root transformed vegetation cover data, and up totals for each species (across dataset), 
#' and turn into relative abundances
veg_spp_tot <- veg %>%
  select(acevul:virfal) %>%
  summarise_all(sum)
veg_spp_tot_sqrt <- sqrt(veg_spp_tot)
veg_spp_tot_rel <- decostand(veg_spp_tot_sqrt, "total")
veg_spp_all <- data.frame(t(veg_spp_tot_rel))
colnames(veg_spp_all)[1] <- "rel_cover"
veg_spp_all$species <- rownames(veg_spp_all)

#' Reduce veg and trait datasets to list of species found in veg
veg_spp_only <- veg_spp_all %>%
  filter(rel_cover >0)
trait_veg <- trait_dat %>%
  filter(species %in% veg_spp_only$species)
veg_only1 <- veg_sqrt_rel1 %>%
  select(site, one_of(veg_spp_only$species))

#' Trait coverage
# Height
height_spp <- trait_veg %>% 
  filter(height_per_day > 0)
veg_height <- veg_spp_only %>%
  filter(species %in% height_spp$species) %>%
  select(-species) %>%
  summarise_all(sum)

#RMR 
rmr_spp <- trait_veg %>% 
  filter(RMR > 0)
veg_rmr <- veg_spp_only %>%
  filter(species %in% rmr_spp$species) %>%
  select(-species) %>%
  summarise_all(sum)

#RDMC
rdmc_spp <- trait_veg %>% 
  filter(RDMC > 0)
veg_rdmc <- veg_spp_only %>%
  filter(species %in% rdmc_spp$species) %>%
  select(-species) %>%
  summarise_all(sum)

#SRL
srl_spp <- trait_veg %>% 
  filter(SRL > 0)
veg_srl <- veg_spp_only %>%
  filter(species %in% srl_spp$species) %>%
  select(-species) %>%
  summarise_all(sum)

#Rdiam
rdiam_spp <- trait_veg %>% 
  filter(Rdiam > 0)
veg_rdiam <- veg_spp_only %>%
  filter(species %in% rdiam_spp$species) %>%
  select(-species) %>%
  summarise_all(sum)

#LDMC
ldmc_spp <- trait_veg %>% 
  filter(LDMC > 0)
veg_ldmc <- veg_spp_only %>%
  filter(species %in% ldmc_spp$species) %>%
  select(-species) %>%
  summarise_all(sum)

#SLA
sla_spp <- trait_veg %>% 
  filter(SLA > 0)
veg_sla <- veg_spp_only %>%
  filter(species %in% sla_spp$species) %>%
  select(-species) %>%
  summarise_all(sum)

#seed mass
smass_spp <- trait_veg %>% 
  filter(seed_mass > 0)
veg_smass <- veg_spp_only %>%
  filter(species %in% smass_spp$species) %>%
  select(-species) %>%
  summarise_all(sum)

#' *Estimate trait coverage in the seedbank*  
#' 

#' Grab seedbank count data, and up totals for each species (across dataset), squareroot transform, 
#' and turn into relative abundances
sb_spp_tot <- sb %>%
  select(acevul:virfal) %>%
  summarise_all(sum)
sb_spp_tot_sqrt <- sqrt(sb_spp_tot)
sb_spp_tot_rel <- decostand(sb_spp_tot_sqrt, "total")
sb_spp_all <- data.frame(t(sb_spp_tot_rel))
colnames(sb_spp_all)[1] <- "rel_cover"
sb_spp_all$species <- rownames(sb_spp_all)

#' Reduce veg and trait datasets to list of species found in veg
sb_spp_only <- sb_spp_all %>%
  filter(rel_cover >0)
trait_sb <- trait_dat %>%
  filter(species %in% sb_spp_only$species)
sb_only1 <- sb_sqrt_rel1 %>%
  select(site, one_of(sb_spp_only$species))

#' Trait coverage
# Height
height_spp <- trait_sb %>% 
  filter(height_per_day > 0)
sb_height <- sb_spp_only %>%
  filter(species %in% height_spp$species) %>%
  select(-species) %>%
  summarise_all(sum)

#RMR 
rmr_spp <- trait_sb %>% 
  filter(RMR > 0)
sb_rmr <- sb_spp_only %>%
  filter(species %in% rmr_spp$species) %>%
  select(-species) %>%
  summarise_all(sum)

#RDMC
rdmc_spp <- trait_sb %>% 
  filter(RDMC > 0)
sb_rdmc <- sb_spp_only %>%
  filter(species %in% rdmc_spp$species) %>%
  select(-species) %>%
  summarise_all(sum)

#SRL
srl_spp <- trait_sb %>% 
  filter(SRL > 0)
sb_srl <- sb_spp_only %>%
  filter(species %in% srl_spp$species) %>%
  select(-species) %>%
  summarise_all(sum)

#Rdiam
rdiam_spp <- trait_sb %>% 
  filter(Rdiam > 0)
sb_rdiam <- sb_spp_only %>%
  filter(species %in% rdiam_spp$species) %>%
  select(-species) %>%
  summarise_all(sum)

#LDMC
ldmc_spp <- trait_sb %>% 
  filter(LDMC > 0)
sb_ldmc <- sb_spp_only %>%
  filter(species %in% ldmc_spp$species) %>%
  select(-species) %>%
  summarise_all(sum)

#SLA
sla_spp <- trait_sb %>% 
  filter(SLA > 0)
sb_sla <- sb_spp_only %>%
  filter(species %in% sla_spp$species) %>%
  select(-species) %>%
  summarise_all(sum)

#seed mass
smass_spp <- trait_sb %>% 
  filter(seed_mass > 0)
sb_smass <- sb_spp_only %>%
  filter(species %in% smass_spp$species) %>%
  select(-species) %>%
  summarise_all(sum)

#' Vegetation cover captured by trait dataset
rbind(veg_height, veg_rmr, veg_rdmc, veg_srl, veg_rdiam, veg_ldmc, veg_sla, veg_smass)
#' Seedbank cover captured by trait dataset
rbind(sb_height, sb_rmr, sb_rdmc, sb_srl, sb_rdiam, sb_ldmc, sb_sla, sb_smass)


#'  *Reduce veg and trait dataframes* to only contain species with at least one trait available
# Cactus and  ascste has to be removed from each
# Remaining species=79

veg_only_t1 <- veg_sqrt %>%
  select(one_of(veg_spp_only$species)) %>%
  select(-cactus, -ascste) 
veg_only_t <- decostand(veg_only_t1, "total")
veg_only_t$site <- veg$site
veg_only_site <- veg_only_t %>%
  group_by(site) %>%
  summarise_all(mean) %>%
  select(-site)
rowSums(veg_only_site)

veg_only_site <- veg_only_site[,order(names(veg_only_site))]

# Turn veg dataframe to matrix
veg_only_site_m <- as.matrix(veg_only_site)

# Check trait correlations for high correlations and linearity
trait_veg <- trait_veg %>%
  filter(!species %in% c("cactus","ascste"))

trait_veg_s <- trait_veg %>%
  select(-species)
round(cor(trait_veg_s, use="pairwise.complete.obs"), 2)
pairs(trait_veg_s) 

# Log-transform SRL and seed mass
trait_veg_s$seed_mass <- log(trait_veg_s$seed_mass)
trait_veg_s$SRL <- log(trait_veg_s$SRL)

# Scale trait data
trait_veg_s <- data.frame(scale(trait_veg_s))
rownames(trait_veg_s) <- trait_veg$species

# View correlations post-transformation
round(cor(trait_veg_s, use="pairwise.complete.obs"), 2)
pairs(trait_veg_s) 



#'  *Reduce sb and trait dataframes* to contain only species with >=1 trait
# Cactus has to be removed from each
# Remaining species=74

sb_only_t1 <- sb_sqrt %>%
  select(one_of(sb_spp_only$species)) %>%
  select(-cactus) 
sb_only_t <- decostand(sb_only_t1, "total")
sb_only_t$site <- sb$site
sb_only_site <- sb_only_t %>%
  group_by(site) %>%
  summarise_all(mean) %>%
  select(-site)
rowSums(sb_only_site)

sb_only_site <- sb_only_site[,order(names(sb_only_site))]

sb_only_site_m <- as.matrix(sb_only_site)

# Check trait correlations for high correlations and linearity
trait_sb <- trait_sb %>%
  filter(species !="cactus")

trait_sb_s <- trait_sb %>%
  select(-species)
round(cor(trait_sb_s, use="pairwise.complete.obs"), 2)
pairs(trait_sb_s) 

# Log-transform SRL and seed mass
trait_sb_s$seed_mass <- log(trait_sb_s$seed_mass)
trait_sb_s$SRL <- log(trait_sb_s$SRL)

# Scale trait data
trait_sb_s <- data.frame(scale(trait_sb_s))
rownames(trait_sb_s) <- trait_sb$species

# Recheck correlations and scatterplots after transformation
round(cor(trait_sb_s, use="pairwise.complete.obs"), 2)
pairs(trait_sb_s) 


#'  *Reduce combined veg/sb and trait dataframes* to contain only species with >=1 trait
# Ascste/Cactus has to be removed from each
# Remaining species=88

veg_sb_t1 <- veg_sb_sqrt  %>%
  select(-cactus, -ascste) 
veg_sb_t <- decostand(veg_sb_t1, "total")
veg_sb_t$site <- veg_sb$site
veg_sb_t$type <- veg_sb$type
veg_sb_site2 <- veg_sb_t %>%
  group_by(site, type) %>%
  summarise_all(mean) %>%
  select(-type) %>%
  group_by(site) %>%
  summarise_all(mean) %>%
  select(-site)
rowSums(veg_sb_site2)

veg_sb_site_m <- as.matrix(veg_sb_site2)

# Check trait correlations for high correlations and linearity
trait_veg_sb <- trait_dat %>%
  filter(!species %in% c("cactus","ascste"))

trait_veg_sb_s <- trait_veg_sb %>%
  select(-species)
round(cor(trait_veg_sb_s, use="pairwise.complete.obs"), 2)
pairs(trait_veg_sb_s) 

# Log-transform SRL and seed mass
trait_veg_sb_s$seed_mass <- log(trait_veg_sb_s$seed_mass)
trait_veg_sb_s$SRL <- log(trait_veg_sb_s$SRL)

# Scale trait data
trait_veg_sb_s <- data.frame(scale(trait_veg_sb_s))
rownames(trait_veg_sb_s) <- trait_veg_sb$species

# Check on correlations after transformations
str(trait_veg_sb_s)
round(cor(trait_veg_sb_s, use="pairwise.complete.obs"), 2)
pairs(trait_veg_sb_s) 


# Save scatterplot matrix
tiff(filename="species_trait_corr.tiff", res=600, width=6.5, height = 6.5, units = "in")
pairs(trait_veg_sb_s) 
dev.off()


#'  
#'  ** Veg functional diversity as a function of soil age rank**  
#'  
#'  Notes on Diversity metrics:
#'  *FDis:* 
#'  With abundance data, this metric is the weighted average distance of all
#'  species in a community to the centroid (whose location is shifted totwards the 
#'  most abundant species). It has no upper limit. By construction, it is unaffected
#'  by species richness.  It is computed from the uncorrected specis matrix, and
#'  any axes with negative eigenvalues are corrected following Anderson 2006. Main 
#'  citation is Laliberte & Legendere 2010. 
#'   
#'  *FRic:*
#'  Convex hull volume of present species in trait space.
#'  Cannot account for species relative abundances (Villeger et al. 2008)
#'  
#'  *FDiv:*
#'  For all species in a community in trait space, this metric calculates an average
#'  distace from the centroid. It then looks at the distribution of species 
#'  differences away from this mean. Specifically, the more abundant species that
#'  are FURTHER from the centroid than average, the higher FDiv is. 
#'  CAN account for species relative abundances (Villeger et al. 2008)
#'  
#'  *FEve:*
#'  If species are distributed in trait space and the shortest tree was drawn that
#'  linked all species in a community, FEve is the regularity of the links between 
#'  species. CAN account for species relative abundances (Villeger et al. 2008)
#'  

#'  
#'  *Vegetative Functional Diversity*  
#'  

#' Calculate Gower distance for trait matrix
gower_trait_v <- gowdis(trait_veg_s)

#' Calculate diversity indices at the site-level  
#' FDiv and FDis metric for veg only, m=13 is highest working tried, with Fric=24% of variation 
#' HOWEVER, the m value does not affect the calculation of FDis, so is set lower here for speed

# Relative abundances
dbfd1 <- dbFD(gower_trait_v, veg_only_site_m, w.abun=TRUE, corr="cailliez", m="min")
dbfd1$qual.FRic

# Presence-Absence
dbfd1_pa <- dbFD(gower_trait_v, veg_only_pa_m, w.abun=TRUE, corr="cailliez", m="min")
dbfd1$qual.FRic


#' Fric
env$fric_v <- dbfd1$FRic

#' Fdiv
env$fdiv_v <- dbfd1$FDiv

#' Fdis
env$fdis_v <- dbfd1$FDis
env$fdis_v_pa <- dbfd1_pa$FDis

#' RaoQ
env$raoq_v <- dbfd1$RaoQ



#' 
#' *Seedbank functional diversity*
#' 

#'  Functional Diversity (across all traits) for each site, sb only  
# Calculate Gower distance for trait matrix
gower_trait_s <- gowdis(trait_sb_s)


#' Calculate diversity indices at the site-level  
# Relative abundances
dbfd2 <- dbFD(gower_trait_s, sb_only_site_m, w.abun=TRUE, corr="cailliez", m="min")
dbfd2$qual.FRic

# Presence-Absence
dbfd2_pa <- dbFD(gower_trait_s, sb_only_pa_m, w.abun=TRUE, corr="cailliez", m="min")

#' Fric
env$fric_s <- dbfd2$FRic

#' Fdiv
env$fdiv_s <- dbfd2$FDiv

#' Fdis
env$fdis_s <- dbfd2$FDis
env$fdis_s_pa <- dbfd2_pa$FDis

#' RaoQ
env$raoq_s <- dbfd2$RaoQ

#'
#'  *Combined (veg + seedbank) Functional Diversity*
#'  


#'  Functional Diversity (across all traits) for each site, veg and seedbank relative abundances combined  
# Calculate Gower distance for trait matrix
gower_trait_vs <- gowdis(trait_veg_sb_s)

# Diversity metrics for veg and seedbank data
dbfd3 <- dbFD(gower_trait_vs, veg_sb_site_m, w.abun=TRUE, corr="cailliez", m='min')
dbfd3$qual.FRic

# Create pres-abs matrix to swap out if needed
veg_sb_site_pa <- veg_sb_site2 %>%
  replace(((.)>0),1)
veg_sb_pa_m <- as.matrix(veg_sb_site_pa)

# Presence-Absence
dbfd3_pa<- dbFD(gower_trait_vs, veg_sb_pa_m, w.abun=TRUE, corr="cailliez", m='min')


#' Fric
env$fric_vs <- dbfd3$FRic

#' Fdiv
env$fdiv_vs <- dbfd3$FDiv

#' Fdis
env$fdis_vs <- dbfd3$FDis
env$fdis_vs_pa <- dbfd3_pa$FDis

#' RaoQ
env$raoq_vs <- dbfd3$RaoQ



#'
#'  **Model: Functional dispersion as a function of type (seedbank or veg) and age rank**
#'  With figure and linear mixed model
#'      

# Create dataframe for FDis figure and analysis
fdis_fig_dat <- env %>%
  select(fdis_v, fdis_s, fdis_vs) %>%
  gather(key="type", value="fdis")
fdis_fig_dat$age_rank <- rep(1:6, each=2)
fdis_fig_dat$elev <- env$elevation_m
fdis_fig_dat$elev_rank <- rep(as.numeric(rank(env$elevation_m)), times=3)
fdis_fig_dat <- data.frame(fdis_fig_dat)
fdis_fig_dat$type <-  factor(fdis_fig_dat$type,
                             levels = c("fdis_v","fdis_s","fdis_vs"),
                             labels = c("Vegetation","Seedbank","Veg+Seedbank"))

fdis_pa <- env %>% 
  select(fdis_v_pa, fdis_s_pa, fdis_vs_pa) %>%
  gather(key="type", value="fdis_pa") %>%
  select(fdis_pa)
fdis_fig_dat <- data.frame(cbind(fdis_fig_dat,fdis_pa))


#' Summary of diversity values across seedbank, veg and seedbank+veg
veg_sb_compare_fdis <- env %>% 
  select(fdis_v, fdis_s, fdis_vs, fdiv_v, fdiv_s, fdiv_vs, fric_v, fric_s, fric_vs, raoq_v, raoq_s, raoq_vs) %>%
  summarise_all(mean)
veg_sb_compare_fdis


#'   *Fig.*: Scatterplot with functional dispersion (y) as a function of terrace elevation rank (x), with sites as 
#'      points colored by seedbank (light), veg (med), seedbank+veg (dark)# 


# Figure by ELEVATION RANK
fdis_fig2 <- ggplot (aes(x=elev_rank, y=fdis, col=type), data=fdis_fig_dat) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  labs(x= "Elevation Rank", y="Functional Dispersion", col="Community type") + 
  scale_color_manual(values=c("#636363","#9ECAE1","black")) 
fdis_fig2

# Figure by ELEVATION RANK w/ PRES-ABS
fdis_fig2_pa <- ggplot (aes(x=elev_rank, y=fdis_pa, col=type), data=fdis_fig_dat) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  labs(x= "Elevation Rank", y="Functional Dispersion (pres-abs)", col="Community type") + 
  geom_text(aes(x=3, y= .235 , label="type p=0.074")) + 
  geom_text(aes(x=3, y= .23, label="elev rank p=0.005")) + 
  geom_text(aes(x=3, y= .225, label="t*e NS")) + 
  scale_color_manual(values=c("#636363","#9ECAE1","black")) 
fdis_fig2_pa


tiff(filename="fdis_fig_elevrank.tiff", res=600, width=6, height = 4, units = "in")
fdis_fig2
dev.off()



#' 
#' *Combined species richness and functional dispersion figure*
#' 

div_fig_dat <- cbind (fdis_fig_dat,rich_fig_dat$richness)
colnames(div_fig_dat)[6] <- "richness"

div_fig_dat <- div_fig_dat %>%
  gather("metric","value", richness, fdis)
div_fig_dat$metric <- factor(div_fig_dat$metric, levels=c("richness","fdis"))


### ELEVATION RANK FIG
# Create text data.frames for annotating plot
div_text <- data.frame(label = c("A)", "B)"), metric = c("richness", "fdis"), elev=c(1,1), value=c(65, 0.2))
stats_text1 <- data.frame(label = c("type p<0.001", "type p<0.001"), 
                          metric = c("richness", "fdis"), 
                          elev=c(10.5,3), value=c(63, 0.2))
stats_text2 <- data.frame(label = c("elev rank p<0.001", "elev rank p=0.022"), 
                          metric = c("richness", "fdis"), 
                          elev=c(10.5,3), value=c(60, 0.195))
stats_text3 <- data.frame(label = c("type*elev rank NS", "type*elev rank p=0.004"), 
                          metric = c("richness", "fdis"), 
                          elev=c(10.5,3), value=c(57, 0.19))


div_fig2 <- ggplot (aes(x=elev_rank, y=value, col=type), data=div_fig_dat) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  labs(x= "Terrace Elevation (rank)", col="Community type") + 
  scale_x_continuous(breaks= seq(from = 1, to = 12, by = 2)) +
  scale_color_manual(values=c("#636363","#9ECAE1","black")) +
  facet_wrap(. ~metric, ncol=1, scale="free", strip.position = "left", labeller = as_labeller(c(fdis = "Functional Dispersion", richness = "Species Richness"))) +
  ylab(NULL) +
  theme_bw() +
  theme(strip.background = element_blank(),strip.placement = "outside", strip.text.y = element_text(size = 11)) +
  geom_text(data=div_text, aes(x=elev, y=value, label=label), inherit.aes = FALSE) +
  geom_text(data=stats_text1, aes(x=elev, y=value, label=label), size=3.2, inherit.aes=FALSE) +
  geom_text(data=stats_text2, aes(x=elev, y=value, label=label), size=3.2, inherit.aes=FALSE) +
  geom_text(data=stats_text3, aes(x=elev, y=value, label=label), size=3.2, inherit.aes=FALSE) 
div_fig2




#' 
#'  *Functional richness & Rao's Q figures for comparison*
#'  


#' *Functional richness*
fric_fig_dat <- env %>%
  select(fric_v, fric_s, fric_vs) %>%
  gather(key="type", value="fric")
fric_fig_dat$age_rank <- rep(1:6, each=2)
fric_fig_dat$elev <- env$elevation_m
fric_fig_dat$elev_rank <- rep(as.numeric(rank(env$elevation_m)),3)
fric_fig_dat <- data.frame(fric_fig_dat)
fric_fig_dat$type <-  factor(fric_fig_dat$type,
                             levels = c("fric_v","fric_s","fric_vs"),
                             labels = c("Vegetation","Seedbank","Veg+Seedbank"))


### Elevation Rank
fric_fig2 <- ggplot (aes(x=elev_rank, y=fric, col=type), data=fric_fig_dat) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  labs(x= "Terrace Elevation (rank)", y="Functional Richness", col="Community type") + 
  scale_color_manual(values=c("#636363","#9ECAE1","black")) +
  annotate("text", x = 0.5, y = 7, label = "B)", size=5, fontface =2) +
  annotate("text", x = 10, y = 7, label = "type p<0.001", size=3)+
  annotate("text", x = 10, y = 6.8, label = "elev rank p<0.001", size=3)+
  annotate("text", x = 10, y = 6.6, label = "type*elev rank NS", size=3) +
  scale_x_continuous(breaks=seq(1,12,2))
fric_fig2


tiff(filename="fric_fig_elev_rank.tiff", res=600, width=6, height = 4, units = "in")
fric_fig2
dev.off()



#' *Rao's Q*
raoq_fig_dat <- env %>%
  select(raoq_v, raoq_s, raoq_vs) %>%
  gather(key="type", value="raoq")
raoq_fig_dat$age_rank <- rep(1:6, each=2)
raoq_fig_dat$elev <- env$elevation_m
raoq_fig_dat$elev_rank <- rep(as.numeric(rank(env$elevation_m)),3)
raoq_fig_dat <- data.frame(raoq_fig_dat)
raoq_fig_dat$type <-  factor(raoq_fig_dat$type,
                             levels = c("raoq_v","raoq_s","raoq_vs"),
                             labels = c("Vegetation","Seedbank","Veg+Seedbank"))


#' Elevation Rank
raoq_fig2 <- ggplot (aes(x=elev_rank, y=raoq, col=type), data=raoq_fig_dat) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  labs(x= NULL, y="Rao's Q", col="Community type") + 
  scale_color_manual(values=c("#636363","#9ECAE1","black")) +
  annotate("text", x = 0.5, y = 0.05, label = "A)", size=5, fontface =2) +
  annotate("text", x = 4, y = 0.05, label = "type p=0.001", size=3)+
  annotate("text", x = 4, y = 0.049, label = "elev rank NS",size=3)+
  annotate("text", x = 4, y = 0.048, label = "type*elev rank p=0.001", size=3) +
  scale_x_continuous(breaks = seq(1,12,2))+
  theme(axis.text = element_text(size=7)) 
raoq_fig2


tiff(filename="raoq_fig_elev_rank.tiff", res=600, width=6, height = 4, units = "in")
raoq_fig2
dev.off()



#' *Combined FRic & RaoQ*

tiff(filename="fric_raoq_fig_elev_rank.tiff", res=600, width=6, height = 8, units = "in")
grid.arrange(raoq_fig2,fric_fig2)
dev.off()



#'  
#'    *Functional diversity models*
#'  

fdis_fig_dat$type <- factor(fdis_fig_dat$type, levels = c("Veg+Seedbank", "Vegetation",  "Seedbank"))


#' *Functional dispersion:*  Run linear model for elevation rank
fdis_lm2 <- lm(fdis ~ type * elev_rank, data=fdis_fig_dat)
summary(fdis_lm2)
Anova(fdis_lm2, type="III")
Anova(fdis_lm2, type="II")

#' *Functional dispersion:*  Run PRES-ABS linear model for elevation rank
fdis_lm2_pa <- lm(fdis_pa ~ type * elev_rank, data=fdis_fig_dat)
summary(fdis_lm2_pa)
Anova(fdis_lm2_pa, type="III")
Anova(fdis_lm2_pa, type="II")

fdis_fig_dat %>%
  select(fdis,type) %>%
  group_by(type) %>%
  summarise (mean_fids = mean(fdis))

rich_fig_dat %>%
  select(richness,type) %>%
  group_by(type) %>%
  summarise (mean_rich = mean(richness))




#' *Functional richness*
fric_fig_dat$site <- env$site

options(contrasts = c("contr.sum", "contr.poly"))
fric_fig_dat$type <- factor(fric_fig_dat$type, levels = c("Seedbank", "Veg+Seedbank", "Vegetation"))

#' LM: Ranked elevation
fric_lm_ranked <- lm(fric ~ type * elev_rank, data=fric_fig_dat)    
summary(fric_lm_ranked) 
Anova(fric_lm_ranked, type="III")
Anova(fric_lm_ranked, type="II")


#' *Rao's  Q*

raoq_fig_dat$site <- env$site
options(contrasts = c("contr.sum", "contr.poly"))
raoq_fig_dat$type <- factor(raoq_fig_dat$type, levels = c("Vegetation","Seedbank", "Veg+Seedbank"))

raoq_fig_dat %>%
  select(raoq,type) %>%
  group_by(type) %>%
  summarise (mean_raoq = mean(raoq))

#' LM: Ranked elevation
raoq_lm_ranked <- lm(raoq ~ type * elev_rank, data=raoq_fig_dat)    
summary(raoq_lm_ranked) 
Anova(raoq_lm_ranked, type="III")
Anova(raoq_lm_ranked, type="II")





#' 
#' **Q3:  Functional trait composition (CWMs)**
#' 
#' Functional composition is assessed via trait community-weighted means (CWMs) and community-weighted variances
#' for each site and community type (seedbank or veg) 




#' *Calculate Vegetation CWMs*
#'  Calculate and explore vegetation CWM's based on weighted abundances
#'  

#' Veg CWMS
fc_v <- functcomp(trait_veg_s, veg_only_site_m)

#' Compare CWMs to age rank
cwm_v <- cbind(env$age_rank, fc_v)
cwm_v 
colnames(cwm_v)[1] <- "age_rank"
cwm_v$age_rank <- as.numeric(cwm_v$age_rank)
cwm_v <- cwm_v[,c("height_per_day","RMR", "RDMC", "SRL", "Rdiam", "LDMC", "SLA","seed_mass","age_rank")]
corr.test(cwm_v)
pairs(cwm_v)

#' Create correlation heat map
cwm_veg_corr <- round(cor(cwm_v, use="pairwise.complete.obs"), 2)
p_cwm_veg <- cor_pmat(cwm_v)
cwm_veg_corr_sig <- ggcorrplot(cwm_veg_corr, title="Weighted CWM correlations in veg (p<0.1 shown)", type="lower", method="circle", p.mat=p_cwm_veg, sig.level=0.10, insig="blank", lab=T, lab_size = 2.5, tl.cex=10)
cwm_veg_corr_sig

#' 
#' *Checking for dominant species influence*
#' Because there are some weird correlations between trait CWMS, I want to check whether 
#' dominant species might be driving trends.
#' 
#' *Approach 1:* Compare results to CWM's based on pres_absence
#'  

#' Create pres/abs matrix and compute veg CWMS
veg_only_pa <- veg_only_site_m %>%
  replace(((.)>0),1)
veg_only_pa_m <- as.matrix(veg_only_pa)
fc_pa_v <- functcomp(trait_veg_s, veg_only_pa_m)

#' Compare CWMs to age rank, explore trends
cwm_pa_v <- cbind(env$age_rank, fc_pa_v)
cwm_pa_v
colnames(cwm_pa_v)[1] <- "age_rank"
cwm_pa_v$age_rank <- as.numeric(cwm_pa_v$age_rank)
cwm_pa_v <- cwm_pa_v[,c("height_per_day","RMR", "RDMC", "SRL", "Rdiam", "LDMC", "SLA","seed_mass","age_rank")]
corr.test(cwm_v)
corr.test(cwm_pa_v)
pairs(cwm_pa_v)

#' Create correlation heat map
cwm_veg_pa_corr <- round(cor(cwm_pa_v, use="pairwise.complete.obs"), 2)
p_cwm_veg_pa <- cor_pmat(cwm_pa_v)
cwm_veg_pa_corr_sig <- ggcorrplot(cwm_veg_pa_corr, title="Pres/Abs CWM correlations in veg (p<0.1 shown)", type="lower", method="circle", p.mat=p_cwm_veg_pa, sig.level=0.10,insig="blank", lab=T, lab_size = 2.5, tl.cex=10)
cwm_veg_pa_corr_sig

#' 
#' *Figure - Weighted vs Pres Absence*
#' View CWM correlations from weighted data VW pres/abs side by side

grid.arrange(cwm_veg_corr_sig, cwm_veg_pa_corr_sig)
tiff(filename="Veg_CWM_correlations.tiff", res=600, width=6, height = 7, units = "in")
grid.arrange(cwm_veg_corr_sig, cwm_veg_pa_corr_sig)
dev.off()

#' 
#' *Approach 2* - Identify dominant species and explore their potential for influence
#' 

#' Histogram of species relative abundances for each site
hist_dat_v_t <- data.frame(t(veg_only_site))
colnames(hist_dat_v_t) <- c("A","B","C","D","E","F","G","H","I","J","K","L")

hist_dat_v <- hist_dat_v_t %>%
  gather(key="site", value="abundance", A:L)

spp_dist_v <- ggplot(data=hist_dat_v, aes(abundance)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~site, scales = 'free_x')
spp_dist_v

tiff(filename="veg_spp_dist.tiff", res=300, width=8, height = 6, units = "in")
spp_dist_v
dev.off()

#' 
#' At each site, there seems to be 1-3 dominant species up around or above 10% cover. Even
#' though evenness estimates were quite high for vegetation, we need to be careful that CWM's
#' aren't simply reflecing changing abundances of dominant species across the age gradient as
#' opposed to trait shifts reflecting the community as a whole 
#' 
#' Dominant vegetation species = 
#' Andger (5-21% cover across the gradient)
#' Chondrosum spp (1-8% cover across the gradient)
#' Muhmon (0-16%)
#' Poacom (0.001-15%)
#' 
#' To be careful, I will first look at the cover of these dominant species across the
#' soil age gradient to see whether they track it significantly. If they do, then I will
#' need to check whether their traits are driving patterns. 

dom_spp <- veg_only_site %>%
  select(andger,chondrosum_spp,muhmon,poacom) %>%
  gather(key="species", value="abundance", andger:poacom)
dom_spp$age_rank <- env$age_rank

#' Regression results: 
#' andger NS (p=0.17), chondrosum NS (p=0.17), **muhmon p<0.001**, poacom NS (p=0.34)
#' 

andger <-dom_spp %>%
  filter(species=="andger")
summary(lm(abundance~age_rank, data=andger))
chondrosum <-dom_spp %>%
  filter(species=="chondrosum_spp")
summary(lm(abundance~age_rank, data=chondrosum))
muhmon <-dom_spp %>%
  filter(species=="muhmon")
summary(lm(abundance~age_rank, data=muhmon))
poacom <-dom_spp %>%
  filter(species=="poacom")
summary(lm(abundance~age_rank, data=poacom))

#' Figure for Dominant species
#' 

dom_text <- data.frame(label = c("p=0.17","p=0.17","p<0.001","p=0.34"), species=c("andger","chondrosum_spp","muhmon","poacom"), age_rank=c(2,2,2,2), abundance=c(0.2,0.2,0.2,0.2))

dom_dist_v <- ggplot(data=dom_spp, aes(x=age_rank,y=abundance)) +
  geom_point(cex=2, col="gray20", alpha=0.5)+
  geom_smooth(method="lm",se=F, col="gray80") +
  labs(x="Soil Age Rank", y="Relative cover (%)") +
  facet_wrap(~species, ncol=2)+
  geom_text(data=dom_text, aes(x=age_rank, y=abundance, label=label), inherit.aes=FALSE)
dom_dist_v

tiff(filename="dom_dist_veg.tiff", res=300, width=6, height = 5, units = "in")
dom_dist_v
dev.off()

#' 
#' In the weighted abundance CWMs, several weird trait correlations appeared. If they were 
#' driven by these dominant species, we might also expect to see combinations of
#' trait values in the dominant species that match CWM correlations.
#' 

#' Weird CWM trait correlations that need to be explored:  
#' 1) High root diameter with low LDMC (r=-0.96!!!)
#' 2) high height with high LDMC and low SLA
#' 3) High RDMC with high SLA (and low RDMC)
#' 

# Create dataframe with dominant species only
trait_dom <- trait %>%
  filter(species %in% c("andger","chondrosum_spp","muhmon","poacom"))
trait_dom$label <- c("andger","chondrosum_spp","muhmon","poacom")

#' 1) See whether any dominant species have high LDMC w/ low RDiam (or visa versa)
#' ... Muhmon does
ggplot() +
  geom_point(aes(x=LDMC, y=Rdiam), col="gray", cex=2, data=trait) +
  geom_point(aes(x=LDMC, y=Rdiam), col="red", cex=2, data=trait_dom) +
  geom_text(data=trait_dom, aes(x=LDMC+1, y=Rdiam+0.03, label=label), col="red")

#' 2) See whether any dominant species have high LDMC w/ low height (or visa versa)
#' ... Muhmon does
ggplot() +
  geom_point(aes(x=LDMC, y=height_per_day), col="gray", cex=2, data=trait) +
  geom_point(aes(x=LDMC, y=height_per_day), col="red", cex=2, data=trait_dom) +
  geom_text(data=trait_dom, aes(x=LDMC+1, y=height_per_day+0.015, label=label), col="red")

#' 3) See whether any dominant species have low LDMC w/ high RDMC (or visa versa)
#' Not really.. and we don't have and RDMC estimate for Muhmon
ggplot() +
  geom_point(aes(x=LDMC, y=RDMC), col="gray", cex=2, data=trait) +
  geom_point(aes(x=LDMC, y=RDMC), col="red", cex=2, data=trait_dom) +
  geom_text(data=trait_dom, aes(x=LDMC+1, y=RDMC+1, label=label), col="red")







#' 
#' 
#' *Calculate Seedbank CWMs*
#'  Calculate and explore vegetation CWM's based on weighted abundances
#'  
#'  

#' Seedbank CWMS
fc_s <- functcomp(trait_sb_s, sb_only_site_m)

#' Compare CWMs to age rank
cwm_s <- cbind(env$age_rank, fc_s)
cwm_s
colnames(cwm_s)[1] <- "age_rank"
cwm_s$age_rank <- as.numeric(cwm_s$age_rank)
cwm_s <- cwm_s[,c("height_per_day","RMR", "RDMC", "SRL", "Rdiam", "LDMC", "SLA","seed_mass","age_rank")]
corr.test(cwm_s)
pairs(cwm_s)

#' Create correlation heat map
cwm_sb_corr <- round(cor(cwm_s, use="pairwise.complete.obs"), 2)
p_cwm_sb <- cor_pmat(cwm_s)
cwm_sb_corr_sig <- ggcorrplot(cwm_sb_corr, title="Weighted CWM correlations in seedbank (p<0.1 shown)", type="lower", method="circle", p.mat=p_cwm_sb, sig.level=0.10, insig="blank", lab=T, lab_size = 2.5, tl.cex=10)
cwm_sb_corr_sig

#' 
#' *Checking for dominant species influence*
#' Because there are some weird correlations between trait CWMS, I want to check whether 
#' dominant species might be driving trends.
#' 
#' *Approach 1:* Compare results to CWM's based on pres_absence
#'  

#' Create pres/abs matrix and compute veg CWMS
sb_only_pa <- sb_only_site_m %>%
  replace(((.)>0),1)
sb_only_pa_m <- as.matrix(sb_only_pa)
fc_pa_s <- functcomp(trait_sb_s, sb_only_pa_m)

#' Compare CWMs to age rank, explore trends
cwm_pa_s <- cbind(env$age_rank, fc_pa_s)
cwm_pa_s
colnames(cwm_pa_s)[1] <- "age_rank"
cwm_pa_s$age_rank <- as.numeric(cwm_pa_s$age_rank)
cwm_pa_s <- cwm_pa_s[,c("height_per_day","RMR", "RDMC", "SRL", "Rdiam", "LDMC", "SLA","seed_mass","age_rank")]
corr.test(cwm_pa_s)
pairs(cwm_pa_s)

#' Create correlation heat map
cwm_sb_pa_corr <- round(cor(cwm_pa_s, use="pairwise.complete.obs"), 2)
p_cwm_sb_pa <- cor_pmat(cwm_pa_s)
cwm_sb_pa_corr_sig <- ggcorrplot(cwm_sb_pa_corr, title="Pres/Abs CWM correlations in seedbank (p<0.1 shown)", type="lower", method="circle", p.mat=p_cwm_sb_pa, sig.level=0.10,insig="blank", lab=T, lab_size = 2.5, tl.cex=10)
cwm_sb_pa_corr_sig

#' 
#'  *Figure - CWM Correlations with Weighted vs Pres Absence data*
#' View CWM correlations from weighted data VW pres/abs side by side

tiff(filename="SB_CWM_correlations.tiff", res=600, width=6, height = 7, units = "in")
grid.arrange(cwm_sb_corr_sig, cwm_sb_pa_corr_sig)
dev.off()


#'
#'
#' *Combine / Save all CWMs*
#' 
#' 

# Add identifying factors to each cwm dataframe
cwm_v$type <- "Vegetation"
cwm_v$abund <- "Abundance"
cwm_s$type <-"Seedbank"
cwm_s$abund <- "Abundance"

cwm_pa_v$type <- "Vegetation"
cwm_pa_v$abund <- "Pres-Abs"
cwm_pa_s$type <-"Seedbank"
cwm_pa_s$abund <- "Pres-Abs"

# Check that column names align, then bind dataframe
colnames(cwm_v)
colnames(cwm_s)
colnames(cwm_pa_v)
colnames(cwm_pa_s)

cwms <- rbind(cwm_v,cwm_s,cwm_pa_v,cwm_pa_s)
cwms
write.csv(cwms,'site_cwms_11Sep20.csv')


#' 
#' 
#'  *Read in CWMs if skipping ahead*
#'  
#'  

cwms <- read.csv('site_cwms_11Sep20.csv')
cwms <- cwms %>% select (-X)
cwms$elev <- env$elevation_m
cwms$elev_rank <- rep(as.numeric(rank(env$elevation_m)),4)
str(cwms)

#'
#' *Figure - CWM differences between seedbank with Weighted vs Pres Abs data*
#'


#' Create dataframe with means and standard errors
cwm_bar_fig_dat_mean <- cwms %>%
  select(-age_rank) %>%
  group_by(type,abund) %>%
  summarise_all(mean) %>%
  gather(key="trait", value="cwm",height_per_day:seed_mass)
cwm_bar_fig_dat_se <- cwms %>%
  select(-age_rank) %>%
  group_by(type,abund) %>%
  summarise_all(se) %>%
  gather(key="trait", value="se",height_per_day:seed_mass)

cwm_bar_fig_dat <- bind_cols(cwm_bar_fig_dat_mean, cwm_bar_fig_dat_se) %>%
  select(-trait1, -abund1, -type1)
cwm_bar_fig_dat$trait <- plyr::revalue(cwm_bar_fig_dat$trait , c("height_per_day"="Height d-1", "seed_mass"="Seedmass"))
cwm_bar_fig_dat$trait <- factor(cwm_bar_fig_dat$trait, levels=c("Seedmass","Rdiam","SRL","RMR","LDMC","RDMC","SLA","Height d-1"))

#' Figure
cwm_barplot <- ggplot(data=cwm_bar_fig_dat) + 
  geom_bar(aes(x=type, y=cwm, fill=type), stat="identity", position=position_dodge()) +
  labs(x=NULL, y = "Standardized CWM") +
  scale_fill_manual(values=c("dodgerblue4","black"), name="Type")+
  geom_errorbar(mapping=aes(x=type, ymin = (cwm-se), ymax = (cwm+se)),width=0.25)+
  theme(axis.text.x = element_text(face="bold", angle=90))+
  facet_grid(rows=vars(abund), cols=vars(trait))+ 
  theme(axis.text.x = element_blank(), axis.ticks.x= element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)
  #annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
cwm_barplot

tiff(filename="CWM_barplot.tiff", res=300, width=10, height = 6, units = "in")
cwm_barplot
dev.off()

#' 
#' Figure for abundance only
#' 

cwm_bar_fig_abund <- cwm_bar_fig_dat %>%
  filter(abund == "Abundance")

cwm_barplot_abund <- ggplot(data=cwm_bar_fig_abund) + 
  geom_bar(aes(x=type, y=cwm, fill=type), stat="identity", position=position_dodge()) +
  labs(x=NULL, y = "Standardized CWM") +
  scale_fill_manual(values=c("dodgerblue4","black"), name="Type")+
  geom_errorbar(mapping=aes(x=type, ymin = (cwm-se), ymax = (cwm+se)),width=0.25)+
  theme(axis.text.x = element_text(face="bold", angle=90))+
  facet_wrap(~trait,nrow=1) + 
  theme(axis.text.x = element_blank(), axis.ticks.x= element_blank()) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, col="gray80")
cwm_barplot_abund

tiff(filename="CWM_barplot_abund.tiff", res=300, width=10, height =4, units = "in")
cwm_barplot_abund
dev.off()




#' 
#' *CWM Linear models*
#' 

#' *Abundance weighted*
cwms_abund <- cwms %>% 
  filter(abund == "Abundance")
cwms_abund$site <- env$site

head(cwms_abund)


#' Elevation Rank - LMs
options(contrasts = c("contr.sum", "contr.poly"))
cwms_abund$type <- factor(cwms_abund$type, levels = c("Vegetation", "Seedbank"))


height_lm_rank <- lm(height_per_day ~ type * elev_rank, data=cwms_abund)    
RMR_lm_rank <- lm(RMR ~ type * elev_rank , data=cwms_abund)    
RDMC_lm_rank <- lm(RDMC ~ type * elev_rank , data=cwms_abund)    
SRL_lm_rank <- lm(SRL ~ type * elev_rank , data=cwms_abund)    
Rdiam_lm_rank <- lm(Rdiam ~ type * elev_rank , data=cwms_abund)    
LDMC_lm_rank <- lm(LDMC ~ type * elev_rank , data=cwms_abund)    
SLA_lm_rank <- lm(SLA ~ type * elev_rank , data=cwms_abund)    
seedmass_lm_rank <- lm(seed_mass~ type * elev_rank , data=cwms_abund)    

Anova(height_lm_rank, type="III")
Anova(RMR_lm_rank, type="III")
Anova(RDMC_lm_rank, type="III")
Anova(SRL_lm_rank, type="III")
Anova(Rdiam_lm_rank, type="III")
Anova(LDMC_lm_rank, type="III")
Anova(SLA_lm_rank, type="III")
Anova(seedmass_lm_rank, type="III")

Anova(height_lm_rank, type="II")
Anova(RMR_lm_rank, type="II")
Anova(RDMC_lm_rank, type="II")
Anova(SRL_lm_rank, type="II")
Anova(Rdiam_lm_rank, type="II")
Anova(LDMC_lm_rank, type="II")
Anova(SLA_lm_rank, type="II")
Anova(seedmass_lm_rank, type="II")




#' *Presence-Absence*
cwms_pa <- cwms %>% 
  filter(abund == "Pres-Abs")
cwms_pa$site <- env$site

head(cwms_pa)


#' 
#' Elevation Rank - LMs
#' 

options(contrasts = c("contr.sum", "contr.poly"))
#options(contrasts = c("contr.treatment", "contr.poly"))

height_lm_rank_pa  <- lm(height_per_day ~ type * elev_rank, data=cwms_pa )    
RMR_lm_rank_pa  <- lm(RMR ~ type * elev_rank , data=cwms_pa )    
RDMC_lm_rank_pa  <- lm(RDMC ~ type * elev_rank , data=cwms_pa )    
SRL_lm_rank_pa  <- lm(SRL ~ type * elev_rank , data=cwms_pa )    
Rdiam_lm_rank_pa  <- lm(Rdiam ~ type * elev_rank , data=cwms_pa )    
LDMC_lm_rank_pa  <- lm(LDMC ~ type * elev_rank , data=cwms_pa )    
SLA_lm_rank_pa  <- lm(SLA ~ type * elev_rank , data=cwms_pa )    
seedmass_lm_rank_pa  <- lm(seed_mass~ type * elev_rank , data=cwms_pa )    

Anova(height_lm_rank_pa, type="III")
Anova(RMR_lm_rank_pa, type="III")
Anova(RDMC_lm_rank_pa , type="III")
Anova(SRL_lm_rank_pa , type="III")
Anova(Rdiam_lm_rank_pa , type="III")
Anova(LDMC_lm_rank_pa , type="III")
Anova(SLA_lm_rank_pa , type="III")
Anova(seedmass_lm_rank_pa , type="III")





#' 
#' **Figures:  CWM & CWV~ Age rank (seedbank and veg shown)**
#' 
#' ** Abundance Weighted**
#' 
#' 
#' 

#'  
#'  *Separate CWM figures*
#'  



#' *Height*

# Height - Elevation Rank (minor sig, type p=0.058)
Anova(height_lm_rank, type="III")
summary(height_lm_rank)

height_cwm2 <- 
  ggplot(data=subset(cwm_abun_fig_dat,trait=="Height d-1")) +
  geom_point(aes(x=elev_rank, y=cwm, col=type,alpha=0.4)) + 
  #geom_smooth(aes(x=elev_rank, y=cwm, col=type), lty=3, size=0.5, method="lm", se=F) +
  geom_boxplot(aes(x=6.5, y=cwm, fill=type), alpha=0.5, width=3.5) +
  scale_color_manual(values=c("dodgerblue3","black")) +
  scale_fill_manual(values=c("dodgerblue3","black")) +
  labs(x=NULL, y="Height d-1 \n(std. CWM)", col="Community type", fill="Community type") +
  scale_x_continuous(breaks=seq(1,12,2), limits=c(1,16)) +
  geom_text(aes(x=14, y=0.35, label = "type p=0.058"), size=3.4) +
  geom_text(aes(x=14, y=0.17, label = "elev rank NS"), size=3.4) +
  geom_text(aes(x=14, y=-0.01, label = "t*e NS"), size=3.4)+
  geom_text(aes(x=1, y=0.35, label="B)"), size=3.4) +
  theme(text = element_text(size=10), axis.text.x=element_blank(),axis.text.y=element_text(size=10)) 
height_cwm2


#' *RMR*

# Elevation Rank - RMR (type p=0.040, elev rank NS, t*e NS)
Anova(RMR_lm_rank, type="III")

rmr_cwm2 <- 
  ggplot(data=subset(cwm_abun_fig_dat,trait=="RMR")) +
  geom_point(aes(x=elev_rank, y=cwm, col=type, alpha=0.4)) + 
  #geom_smooth(aes( x=elev_rank, y=cwm, col=type), lty=3, size=1, method="lm", se=F) +
  geom_boxplot(aes(x=6.5, y=cwm, fill=type), alpha=0.5, width=3.5) +
  scale_color_manual(values=c("dodgerblue3","black")) +
  scale_fill_manual(values=c("dodgerblue3","black")) +
  labs(x=NULL, y="RMR \n(std. CWM)", col="Community type", fill="Community type") +
  scale_x_continuous(breaks=seq(1,12,2), limits=c(1,16)) +
  geom_text(aes(x=14, y=0.6, label = "type p=0.026"), size=3.4) +
  geom_text(aes(x=14, y=0.41, label = "elev rank NS"), size=3.4) +
  geom_text(aes(x=14, y=0.24, label = "t*e NS"), size=3.4)+
  geom_text(aes(x=1, y=0.6, label="C)"), size=3.4) +
  theme(text = element_text(size=10), axis.text.x=element_blank(),axis.text.y=element_text(size=10))
rmr_cwm2



#' *RDMC*

# Elev Rank - RDMC (type NS, elev rank NS, t*e NS)
Anova(RDMC_lm_rank, type="III")

rdmc_cwm2 <- 
  ggplot(data=subset(cwm_abun_fig_dat,trait=="RDMC")) +
  geom_point(aes(x=elev_rank, y=cwm, col=type)) + 
  #geom_smooth(aes( x=elev_rank, y=cwm, col=type), method="lm", se=F) +
  #geom_boxplot(aes(x=3.5, y=cwm, fill=type), alpha=0.2, width=3) +
  scale_color_manual(values=c("dodgerblue3","black")) +
  #scale_fill_manual(values=c("dodgerblue3","black")) +
  labs(x=NULL, y="RDMC \n(std. CWM)", col="Community type", fill="Community type") +
  scale_x_continuous(breaks=seq(1,12,2), limits=c(1,16)) +
  geom_text(aes(x=14, y=0.6, label = "type NS"), size=3.4) +
  geom_text(aes(x=14, y=0.47, label = "elev rank NS"), size=3.4) +
  geom_text(aes(x=14, y=0.34, label = "t*e NS"), size=3.4)+
  geom_text(aes(x=1, y=0.6, label="A)"), size=3.4) +
  theme(text = element_text(size=10), axis.text.x=element_blank(),axis.text.y=element_text(size=10))
rdmc_cwm2




#' *SRL*

#' Elev Rank - SRL (type p<0.001, elev_rank NS, t*e NS)
Anova(SRL_lm_rank, type="III")

srl_cwm2 <- 
  ggplot(data=subset(cwm_abun_fig_dat,trait=="SRL")) +
  geom_point(aes(x=elev_rank, y=cwm, col=type, alpha=0.4)) + 
  #geom_smooth(aes( x=elev_rank, y=cwm, col=type), lty=2, size=1, method="lm", se=F) +
  geom_boxplot(aes(x=6.5, y=cwm, fill=type), alpha=0.5, width=3.5) +
  scale_color_manual(values=c("dodgerblue3","black")) +
  scale_fill_manual(values=c("dodgerblue3","black")) +
  labs(x=NULL, y="SRL \n(std. CWM)", col="Community type", fill="Community type") +
  scale_x_continuous(breaks=seq(1,12,2), limits=c(1,16)) +
  geom_text(aes(x=14, y=0.9, label = "type p=0.01"), size=3.4) +
  geom_text(aes(x=14, y=0.7, label = "elev rank NS"), size=3.4) +
  geom_text(aes(x=14, y=0.5, label = "t*e NS"), size=3.4)+
  geom_text(aes(x=1, y=0.9, label="D)"), size=3.4) +
  theme(text = element_text(size=10), axis.text.x=element_blank(),axis.text.y=element_text(size=10))
srl_cwm2


#' *Rdiam*

# Elev Rank - Rdiam (type p =0.016, soil age p=0.81, t*s NS)
Anova(Rdiam_lm_rank, type="III")

rdiam_cwm2 <- ggplot(data=subset(cwm_abun_fig_dat,trait=="Rdiam")) +
  geom_point(aes(x=elev_rank, y=cwm, col=type)) + 
  geom_smooth(aes( x=elev_rank, y=cwm, col=type), method="lm", se=F) +
  #geom_boxplot(aes(x=3.5, y=cwm, fill=type), alpha=0.2, width=3) +
  scale_color_manual(values=c("dodgerblue3","black")) +
  #scale_fill_manual(values=c("dodgerblue3","black")) +
  labs(x=NULL, y="RDiam \n(std. CWM)", col="Community type", fill="Community type") +
  scale_x_continuous(breaks=seq(1,12,2), limits=c(1,16)) +
  geom_text(aes(x=14, y=0.6, label = "type p=0.016"), size=3.4) +
  geom_text(aes(x=14, y=0.4, label = "elev rank p=0.002"), size=3.4) +
  geom_text(aes(x=14, y=0.2, label = "t*e NS"), size=3.4)+
  geom_text(aes(x=1, y=0.6, label="E)"), size=3.4) +
  theme(text = element_text(size=10), axis.text.x=element_blank(),axis.text.y=element_text(size=10))
rdiam_cwm2


#' *LDMC*

# Elev Rank - LDMC (type p=0.079, soil age p=.065, t*s p=0.069))
Anova(LDMC_lm_rank, type="III")
Anova(LDMC_lm_rank)

ldmc_cwm2 <- ggplot(data=subset(cwm_abun_fig_dat,trait=="LDMC")) +
  geom_point(aes(x=elev_rank, y=cwm, col=type)) + 
  geom_smooth(aes( x=elev_rank, y=cwm, col=type), method="lm", se=F) +
  #geom_boxplot(aes(x=3.5, y=cwm, fill=type), alpha=0.2, width=3) +
  scale_color_manual(values=c("dodgerblue3","black")) +
  #scale_fill_manual(values=c("dodgerblue3","black")) +
  labs(x=NULL, y="LDMC \n(std. CWM)", col="Community type", fill="Community type") +
  scale_x_continuous(breaks=seq(1,12,2), limits=c(1,16)) +
  geom_text(aes(x=14, y=0.7, label = "type p=0.079"), size=3.4) +
  geom_text(aes(x=14, y=0.5, label = "elev rank p<0.001"), size=3.4) +
  geom_text(aes(x=14, y=0.3, label = "t*e p=0.069"), size=3.4)+
  geom_text(aes(x=1, y=0.7, label="G)"), size=3.4) +
  theme(text = element_text(size=10), axis.text.x=element_blank(),axis.text.y=element_text(size=10))
ldmc_cwm2


#' *SLA*

#  Elev Rank - SLA (type p=.049, elev rank p=.045, t*e p=0.026)
Anova(SLA_lm_rank, type="III")
summary(SLA_lm_rank)
Anova(SLA_lm_rank)

sla_cwm2 <- 
  ggplot(data=subset(cwm_abun_fig_dat,trait=="SLA")) +
  geom_point(aes(x=elev_rank, y=cwm, col=type)) + 
  geom_smooth(aes( x=elev_rank, y=cwm, col=type), method="lm", se=F) +
  #geom_boxplot(aes(x=3.5, y=cwm, fill=type), alpha=0.2, width=3) +
  scale_color_manual(values=c("dodgerblue3","black")) +
  #scale_fill_manual(values=c("dodgerblue3","black")) +
  labs(x="Terrace Elevation (rank)", y="SLA \n (std. CWM)", col="Community type", fill="Community type") +
  scale_x_continuous(breaks=seq(1,12,2), limits=c(1,16)) +
  geom_text(aes(x=14, y=1, label = "type p=0.049"), size=3.4) +
  geom_text(aes(x=14, y=0.82, label = "elev rank NS"), size=3.4) +
  geom_text(aes(x=14, y=0.64, label = "t*e p=0.026"), size=3.4)+
  geom_text(aes(x=1, y=1, label="H)"), size=3.4) +
  theme(text = element_text(size=10), axis.text.x=element_text(size=10),axis.text.y=element_text(size=10))
sla_cwm2



#' *Seedmass*

# Elev Rank - Seedmass (type p<0.001, elev rank p=0.03, t*e NS)
Anova(seedmass_lm_rank, type="III")
Anova(seedmass_lm_rank)
summary(seedmass_lm_rank)

seedmass_cwm2 <- 
  ggplot(data=subset(cwm_abun_fig_dat,trait=="Seedmass")) +
  geom_point(aes(x=elev_rank, y=cwm, col=type)) + 
  geom_smooth(aes( x=elev_rank, y=cwm, col=type), method="lm", se=F) +
  #geom_boxplot(aes(x=3.5, y=cwm, fill=type), alpha=0.2, width=3) +
  scale_color_manual(values=c("dodgerblue3","black")) +
  #scale_fill_manual(values=c("dodgerblue3","black")) +
  labs(x=NULL, y=" Seed mass \n(std. CWM)", col="Community type", fill="Community type") +
  scale_x_continuous(breaks=seq(1,12,2), limits=c(1,16)) +
  geom_text(aes(x=14, y=0.3, label = "type p<0.001"), size=3.4) +
  geom_text(aes(x=14, y=0.07, label = "elev rank p=0.003"), size=3.4) +
  geom_text(aes(x=14, y=-0.16, label = "t*e NS"), size=3.4) +
  geom_text(aes(x=1, y=0.3, label="F)"), size=3.4) +
  theme(text = element_text(size=10), axis.text.x=element_blank(),axis.text.y=element_text(size=10))
seedmass_cwm2



#' *Save plot grids*
elev_rank_cwm <- plot_grid(rdmc_cwm2, height_cwm2, rmr_cwm2, srl_cwm2, rdiam_cwm2, seedmass_cwm2, ldmc_cwm2, sla_cwm2, ncol=1, rel_heights = c(1,1,1,1,1,1,1,1.37))
elev_rank_cwm

tiff(filename="CWM_elev_rank_III.tiff", res=600, width=7, height = 14, units = "in")
elev_rank_cwm
dev.off()



#'
#' *CWM figures for pres-abs*
#' 

cwm_pres_fig_dat <- cwms_pa %>%
  select(-site) %>%
  gather(key="trait", value="cwm", height_per_day:seed_mass)

# Recode traits
cwm_pres_fig_dat$trait <- plyr::revalue(cwm_pres_fig_dat$trait, c("height_per_day"="Height d-1", "seed_mass"="Seedmass"))
str(cwm_pres_fig_dat)


#' *Height*

# Height - Elevation Rank (all NS)
Anova(height_lm_rank_pa, type="III")
height_cwm_pa2 <- 
  ggplot(data=subset(cwm_pres_fig_dat,trait=="Height d-1")) +
  geom_point(aes(x=elev_rank, y=cwm, col=type)) + 
  #geom_smooth(aes( x=age_rank, y=cwm, col=type), lty=2, size=0.5, method="lm", se=F) +
  #geom_boxplot(aes(x=3.5, y=cwm, fill=type), alpha=0.2, width=3) +
  scale_color_manual(values=c("dodgerblue3","black")) +
  scale_fill_manual(values=c("dodgerblue3","black")) +
  labs(x=NULL, y="Height d-1 \n(std. CWM)", col="Community type", fill="Community type") +
  scale_x_continuous(breaks=seq(1,12,2), limits=c(1,16)) +
  geom_text(aes(x=14, y=0.35, label = "type NS"), size=3.4) +
  geom_text(aes(x=14, y=0.17, label = "elev rank NS"), size=3.4) +
  geom_text(aes(x=14, y=-0.01, label = "t*e NS"), size=3.4)+
  geom_text(aes(x=1, y=0.35, label="B)"), size=3.4) +
  theme(text = element_text(size=10), axis.text.x=element_blank(),axis.text.y=element_text(size=10)) 
height_cwm_pa2



#' *RMR*

# Elevation Rank - RMR (type NS, elev rank p=0.03, t*e NS)
Anova(RMR_lm_rank_pa, type="III")
rmr_cwm_pa2 <- 
  ggplot(data=subset(cwm_pres_fig_dat,trait=="RMR")) +
  geom_point(aes(x=elev_rank, y=cwm, col=type)) + 
  geom_smooth(aes( x=elev_rank, y=cwm), lty=1, color="gray60", method="lm", se=F) +
  #geom_boxplot(aes(x=3.5, y=cwm, fill=type), alpha=0.5, width=3) +
  scale_color_manual(values=c("dodgerblue3","black")) +
  scale_fill_manual(values=c("dodgerblue3","black")) +
  labs(x=NULL, y="RMR \n(std. CWM)",  col="Community type", fill="Community type") +
  scale_x_continuous(breaks=seq(1,12,2), limits=c(1,16)) +
  geom_text(aes(x=14, y=0.6, label = "type NS"), size=3.4) +
  geom_text(aes(x=14, y=0.43, label = "elev rank p=0.006"), size=3.4) +
  geom_text(aes(x=14, y=0.27, label = "t*e NS"), size=3.4)+
  geom_text(aes(x=1, y=0.6, label="C)"), size=3.4) +
  theme(text = element_text(size=10), axis.text.x=element_blank(),axis.text.y=element_text(size=10))
rmr_cwm_pa2


#' *RDMC*

# Elev Rank - RDMC (type NS, elev rank NS, t*e NS)
Anova(RDMC_lm_rank_pa, type="III")
rdmc_cwm_pa2 <- 
  ggplot(data=subset(cwm_pres_fig_dat,trait=="RDMC")) +
  geom_point(aes(x=elev_rank, y=cwm, col=type)) + 
  #geom_smooth(aes( x=elev_rank, y=cwm, col=type), method="lm", se=F) +
  #geom_boxplot(aes(x=3.5, y=cwm, fill=type), alpha=0.2, width=3) +
  scale_color_manual(values=c("dodgerblue3","black")) +
  #scale_fill_manual(values=c("dodgerblue3","black")) +
  labs(x=NULL, y="RDMC \n(std. CWM)", col="Community type", fill="Community type") +
  scale_x_continuous(breaks=seq(1,12,2), limits=c(1,16)) +
  geom_text(aes(x=14, y=0.6, label = "type NS"), size=3.4) +
  geom_text(aes(x=14, y=0.47, label = "elev rank NS"), size=3.4) +
  geom_text(aes(x=14, y=0.34, label = "t*e NS"), size=3.4)+
  geom_text(aes(x=1, y=0.6, label="A)"), size=3.4) +
  theme(text = element_text(size=10), axis.text.x=element_blank(),axis.text.y=element_text(size=10))
rdmc_cwm_pa2



#' *SRL*


#' Elev Rank - SRL (type NS, elev_rank p<0.001, t*e p<0.001)
Anova(SRL_lm_rank_pa, type="III")
srl_cwm_pa2 <- 
  ggplot(data=subset(cwm_pres_fig_dat,trait=="SRL")) +
  geom_point(aes(x=elev_rank, y=cwm, col=type)) + 
  geom_smooth(aes( x=elev_rank, y=cwm, col=type), lty=1, method="lm", se=F) +
  #geom_boxplot(aes(x=3.5, y=cwm, fill=type), alpha=0.5, width=3) +
  scale_color_manual(values=c("dodgerblue3","black")) +
  scale_fill_manual(values=c("dodgerblue3","black")) +
  labs(x=NULL, y="SRL \n(std. CWM)", col="Community type", fill="Community type") +
  scale_x_continuous(breaks=seq(1,12,2), limits=c(1,16)) +
  geom_text(aes(x=14, y=0.9, label = "type NS"), size=3.4) +
  geom_text(aes(x=14, y=0.7, label = "elev rank p<0.001"), size=3.4) +
  geom_text(aes(x=14, y=0.5, label = "t*e p<0.001"), size=3.4)+
  geom_text(aes(x=1, y=0.9, label="D)"), size=3.4) +
  theme(text = element_text(size=10), axis.text.x=element_blank(),axis.text.y=element_text(size=10))
srl_cwm_pa2




#' *Rdiam*

# Elev Rank - Rdiam (type p NS, soil age p<0.001, t*s p=0.071)
Anova(Rdiam_lm_rank_pa, type="III")

rdiam_cwm_pa2 <- ggplot(data=subset(cwm_pres_fig_dat,trait=="Rdiam")) +
  geom_point(aes(x=elev_rank, y=cwm, col=type)) + 
  geom_smooth(aes( x=elev_rank, y=cwm, col=type), method="lm", se=F) +
  #geom_boxplot(aes(x=3.5, y=cwm, fill=type), alpha=0.2, width=3) +
  scale_color_manual(values=c("dodgerblue3","black")) +
  #scale_fill_manual(values=c("dodgerblue3","black")) +
  labs(x=NULL, y="RDiam \n(std. CWM)", col="Community type", fill="Community type") +
  scale_x_continuous(breaks=seq(1,12,2), limits=c(1,16)) +
  geom_text(aes(x=14, y=0.6, label = "type NS"), size=3.4) +
  geom_text(aes(x=14, y=0.4, label = "elev rank p<0.001"), size=3.4) +
  geom_text(aes(x=14, y=0.2, label = "t*e p=0.071"), size=3.4)+
  geom_text(aes(x=1, y=0.6, label="E)"), size=3.4) +
  theme(text = element_text(size=10), axis.text.x=element_blank(),axis.text.y=element_text(size=10))
rdiam_cwm_pa2



#' *LDMC*

# Elev Rank - LDMC (type NS, elev rank NS, t*e p=0.008))
Anova(LDMC_lm_rank_pa, type="III")
ldmc_cwm_pa2 <- ggplot(data=subset(cwm_pres_fig_dat,trait=="LDMC")) +
  geom_point(aes(x=elev_rank, y=cwm, col=type)) + 
  geom_smooth(aes( x=elev_rank, y=cwm, col=type), method="lm", se=F) +
  #geom_boxplot(aes(x=3.5, y=cwm, fill=type), alpha=0.2, width=3) +
  scale_color_manual(values=c("dodgerblue3","black")) +
  #scale_fill_manual(values=c("dodgerblue3","black")) +
  labs(x=NULL, y="LDMC \n(std. CWM)", col="Community type", fill="Community type") +
  scale_x_continuous(breaks=seq(1,12,2), limits=c(1,16)) +
  geom_text(aes(x=14, y=0.7, label = "type NS"), size=3.4) +
  geom_text(aes(x=14, y=0.5, label = "elev rank .026"), size=3.4) +
  geom_text(aes(x=14, y=0.3, label = "t*e p=0.008"), size=3.4)+
  geom_text(aes(x=1, y=0.7, label="G)"), size=3.4) +
  theme(text = element_text(size=10), axis.text.x=element_blank(),axis.text.y=element_text(size=10))
ldmc_cwm_pa2




#' *SLA*

#  Elev Rank - SLA (type NS, elev rank NS, t*e NS)
Anova(SLA_lm_rank_pa, type="III")
sla_cwm_pa2 <- 
  ggplot(data=subset(cwm_pres_fig_dat,trait=="SLA")) +
  geom_point(aes(x=elev_rank, y=cwm, col=type)) + 
  #geom_smooth(aes( x=elev_rank, y=cwm, col=type), method="lm", se=F) +
  #geom_boxplot(aes(x=3.5, y=cwm, fill=type), alpha=0.2, width=3) +
  scale_color_manual(values=c("dodgerblue3","black")) +
  #scale_fill_manual(values=c("dodgerblue3","black")) +
  labs(x="Terrace Elevation (rank)", y="SLA \n (std. CWM)", col="Community type", fill="Community type") +
  scale_x_continuous(breaks=seq(1,12,2), limits=c(1,16)) +
  geom_text(aes(x=14, y=.48, label = "type NS"), size=3.4) +
  geom_text(aes(x=14, y=0.35, label = "elev rank NS"), size=3.4) +
  geom_text(aes(x=14, y=0.22, label = "t*e NS"), size=3.4)+
  geom_text(aes(x=1, y=.5, label="H)"), size=3.4) +
  theme(text = element_text(size=10), axis.text.x=element_text(size=10),axis.text.y=element_text(size=10))
sla_cwm_pa2



#' *Seedmass*

# Elev Rank - Seedmass (type p<0.001, elev rank NS, t*e NS)
Anova(seedmass_lm_rank_pa, type="III")
seedmass_cwm_pa2 <- 
  ggplot(data=subset(cwm_pres_fig_dat,trait=="Seedmass")) +
  geom_point(aes(x=elev_rank, y=cwm, col=type)) + 
  geom_smooth(aes( x=elev_rank, y=cwm, col=type), method="lm", se=F) +
  #geom_boxplot(aes(x=6.5, y=cwm, fill=type), alpha=0.5, width=3.5) +
  scale_color_manual(values=c("dodgerblue3","black")) +
  scale_fill_manual(values=c("dodgerblue3","black")) +
  labs(x=NULL, y=" Seed mass \n(std. CWM)", col="Community type", fill="Community type") +
  scale_x_continuous(breaks=seq(1,12,2), limits=c(1,16)) +
  geom_text(aes(x=14, y=0.3, label = "type p=0.024"), size=3.4) +
  geom_text(aes(x=14, y=0.07, label = "elev rank p=.041"), size=3.4) +
  geom_text(aes(x=14, y=-0.16, label = "t*e NS"), size=3.4) +
  geom_text(aes(x=1, y=0.3, label="F)"), size=3.4) +
  theme(text = element_text(size=10), axis.text.x=element_blank(),axis.text.y=element_text(size=10))
seedmass_cwm_pa2



#' *Save plot grids*
elev_rank_cwm_pa <- plot_grid(rdmc_cwm_pa2, height_cwm_pa2,rmr_cwm_pa2, srl_cwm_pa2, rdiam_cwm_pa2, seedmass_cwm_pa2, ldmc_cwm_pa2, sla_cwm_pa2, ncol=1, rel_heights = c(1,1,1,1,1,1,1,1.37))
elev_rank_cwm_pa


tiff(filename="CWM_elev_rank_pa_III.tiff", res=600, width=7, height = 14, units = "in")
elev_rank_cwm_pa
dev.off()

