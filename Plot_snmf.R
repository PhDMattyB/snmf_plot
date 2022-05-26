##################################################################
## snmf and ggplot admixture
##
## Matt Brachmann (PhDMattyB)
##
## 2022-05-26
##
##################################################################

## load in required packages
library(BiocManager)
library(tidyverse)
# BiocManager::install('LEA')
library(LEA)
library(reshape2)
library(viridis)

## set working directory
setwd('~/Charr_Adaptive_Introgression/Charr_Project_1/Snmf/')


# Run snmf ----------------------------------------------------------------

## snmf needs to use a .geno file format
## ped2geno converts the ped file into the geno file
## this function can be a bit of a pain in the ass if
## the ped file is even slightly off
# convert = ped2geno('Charr_Lab_ChrAC08.ped',
#                    'Charr_Lab_AC08.geno')

snmf('Charr_Poly.geno',
     K = 1:4,
     entropy = T,
     repetitions = 10,
     project = 'new')


# Load snmf project -------------------------------------------------------
## use this section after running snmf

## REMEMBER when loading a project that you ran in another place
## you need to update the filepath in the .snmfProjectFile
project = load.snmfProject("Charr_Lab.snmfProject")
summary(project)

## plots the cross entropy coefficient
## want to pick the k value with the lowest 
## cross entropy coefficient
plot(project,
     cex = 1.2,
     col = "black",
     pch = 19)

## ce just specifies the best k value with the lowest ce
ce = cross.entropy(project, K = 4)
best_ce = which.min(ce)

## Gets the matrix of q-values for the best k value
qmatrix = Q(project, 
            K = 4, 
            run = best_ce)


## This will show the clusters each individual was assigned to
## The order is the same as the .ped file!!
## Take the first two columns of the ped file to group by population

apply(qmatrix, 1, which.max) %>%
  as_tibble() %>%
  dplyr::rename(Genetic_group = value) %>%
  write_tsv('Charr_poly_snmf_groups_k4.txt')


## write the qvalues from the snmf output to a csv
## just so we have the data for later
as_tibble(qmatrix) %>%
  dplyr::rename(Q1 = V1,
                Q2 = V2,
                Q3 = V3,
                Q4 = V4) %>% 
  write_csv('Charr_poly_snmf_qvalues_k4.csv')


# Mr. clean that data  ----------------------------------------------------------

## Need to create some identifiers for the plots

## Read in the ped file with the population and inidividual ids
popn_data = read_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/GeneticData/Charr_Poly_03.22.2021.ped', 
                     col_names = F) %>% 
  dplyr::select(X1, 
                X2) %>% 
  rename(FID = X1, 
         IndivID = X2)

## genetic group information from snmf
## the individuals are ordered the same as the ped file
group = read_tsv('Charr_poly_snmf_groups_k4.txt')

## make a df
identifiers = bind_cols(popn_data, 
                        group)

## want that lat and long data
environment_data = read_csv('~/Charr_Adaptive_Introgression/Charr_Project_1/SampleSiteData/SampleSites_Coords_1June2020.csv') %>% 
  dplyr::rename(FID = 1) %>% 
  dplyr:: select(FID, 
                 Lat, 
                 Long)

## make a clean df
identifiers = full_join(identifiers, 
                        environment_data, 
                        by = 'FID')
qvalues = read_csv('Charr_poly_snmf_qvalues_k4.csv')

## add in the qvalues to estimate admixture
## and we've got a dataframe
snmf_data = bind_cols(identifiers, 
                      qvalues) %>% 
  arrange(Lat)


# Plot that shit ----------------------------------------------------------

## Use the melt function from reshape2 to get
## the data into the proper orientation
snmf_melted = melt(snmf_data, 
                   id.vars = c('FID', 
                               'IndivID', 
                               'Genetic_group', 
                               'Lat',
                               'Long')) %>% 
  as_tibble()

## bw is the king of plots
theme_set(theme_bw())

## plot is currently arranged by Latitude 
## See arrange function in snmf_data

## need a colour palette
test_col = c( '#4E9EBF',
              '#4E458C',
              '#F23545', 
              '#F29F05')
## If custom colour palettes aren't your thing
## uncomment the other scale_fill_manual() function
## it plots it with viridis

## snmf latitude plot
snmf_poly_k4 = ggplot(data = snmf_melted, 
                      aes(x = reorder(IndivID, Lat),
                          y = value, 
                          fill = variable, 
                          group = Lat))+
  geom_bar(stat = "identity", 
           width = 1)+
  scale_fill_manual(values = test_col)+
  # scale_fill_manual(values = magma(n = 4))+
  labs(x = 'Individuals', 
       y = 'Ancestry proportion')+
  theme(axis.text.y = element_text(color = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        ## can add xaxis labels if needed
        # axis.text.x = element_text(angle = 90,
        #                            hjust = 1,
        #                            vjust = -0.09,
        #                            size = 6,
        #                            color = 'black'),
        legend.position = 'none')+
  scale_x_discrete(guide = guide_axis(n.dodge = 5))+
  scale_y_continuous(expand = c(0,0))

snmf_poly_k4


# ggsave my game ----------------------------------------------------------

## ggsave call to make the best saved plot possible
## always specify dpi=retina for highest quality
ggsave('snmf_charr_poly_k4_noaxis_glacialcols.tiff',
       path = '~/Charr_Adaptive_Introgression/Charr_Project_1/Figures/',
       plot = snmf_poly_k4, 
       dpi = 'retina', 
       units = 'cm', 
       width = 25, 
       height = 15)


