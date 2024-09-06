#-----------------------#
#title: "Vivli 2024"
#output: .
#date: "August-2024"
#-----------------------#

#------------------------------------------------------------------------#
# Load necessary libraries ######
library(ggplot2)
library(dplyr)
library(readxl)
library(lubridate)
library(tidyr)
library(readr)
library(zoo)
library(changepoint)
library(broom)
library(purrr)
library(writexl)
library(mgcv)
library(gratia)
library(pscl)
library(MASS)
library(gamlss)
library(gamlss.dist)
library(rlang)
library(ggthemes)
library(grid)
library(cowplot)
library(sf)
library(usmap)
library(viridisLite)
library(viridis)
library(reshape2)
library(sp)
library(rgdal)
library(spdep)
library(PROJ)
library(proj4)
library(tigris)
library(pscl)
library(TMB)
library(glmmTMB)
library(rnaturalearth)
library(rnaturalearthdata)
library(broom)
library(forcats)
library(scales) 
library(patchwork)
library(openxlsx)
library(RColorBrewer)
library(ggrepel)
#######

#------------------------------------------------------------------------#
####DEFINE WORKING DESK/PATH##########
#Define base path:
base_path <- "/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_projects/0_UniversityofOxford/Vivli/data/data"
base_pathOut <- "/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_projects/0_UniversityofOxford/Vivli/data/data/Outputs/ATLAS"
######
#------------------------------------------------------------------------#
########Loading datasets:########
data_atlasp <- paste0(base_path, "/atlas/2024_05_28 atlas_antibiotics.csv")
data_atlas<- read_csv(data_atlasp)

#Creating Age groups variable
data_atlas <- data_atlas %>%
  mutate(Agegroup = case_when(
    `Age Group` %in% c('0 to 2 Years', '13 to 18 Years', '3 to 12 Years') ~ 0,
    `Age Group` == '19 to 64 Years' ~ 1,
    `Age Group` %in% c('65 to 84 Years', '85 and Over') ~ 2,
    TRUE ~ NA_real_  # Adds NA for any cases not covered above
  ))

data_atlas_us <- subset(data_atlas, Country == "United States")
data_atlas_us_orig <- subset(data_atlas, Country == "United States")
variable_names <- names(data_atlas_us)


unique_states <- unique(data_atlas_us$State)
data_atlas_us <- data_atlas_us %>%
  mutate(State = ifelse(State == "Washington DC", "District of Columbia", State))
unique(data_atlas_us$State)

european_countries <- c("France", "Spain", "Belgium", "Italy", "Germany", "Ireland", "Portugal", 
                        "Greece", "United Kingdom", "Poland", "Switzerland", "Hungary", "Austria", 
                        "Finland", "Denmark", "Sweden", "Croatia", "Czech Republic", "Netherlands", 
                        "Russia", "Romania", "Turkey", "Latvia", "Lithuania", "Serbia", 
                        "Ukraine", "Slovenia", "Bulgaria", "Norway", "Slovak Republic", "Estonia")
countries_list_eu_includ <- c(
  "Austria", "Belgium", "Croatia", "Czech Republic",
  "Denmark", "Finland", "France", "Germany",
  "Greece", "Hungary", "Ireland", "Italy", "Latvia",
  "Lithuania", "Netherlands", "Poland", "Portugal",
  "Romania",  "Slovak Republic", "Slovenia",
  "Spain", "Sweden", "Switzerland",
  "United Kingdom"
)

# Filter data_atlas to only include rows where 'Country' is in the list (more than 3 years being reported)
data_atlas_eu <- data_atlas %>%
  filter(Country %in% countries_list_eu_includ)

states_list_us_NOTinclud <- c(
  "Arkansas", "Connecticut", "Delaware", "District of Columbia", "Hawaii", "Maine", "Minnesota","New Mexico","South Carolina","Montana", "New Hampshire","Oklahoma","Vermont", "West Virginia"
)

#print(variable_names)
#####

#------------------------------------------------------------------------#
####DATA Management, fixing antibiotics, creating resistance/susceptible, country and infection syndrome selection ##########
#Adjusting name of ATBs and changing antibiotic values to 0/1.

antibiotics_to_modify <- c("Amikacin_I", "Ampicillin_I", "Cefepime_I", "Amoxycillin clavulanate_I", 
                         "Azithromycin_I", "Ceftazidime_I", "Cefoxitin_I", "Ceftriaxone_I", 
                         "Clarithromycin_I", "Clindamycin_I", "Erythromycin_I", "Levofloxacin_I", 
                         "Imipenem_I", "Linezolid_I", "Meropenem_I", "Minocycline_I", 
                         "Metronidazole_I", "Penicillin_I", "Piperacillin tazobactam_I", 
                         "Tigecycline_I", "Vancomycin_I", "Aztreonam_I", "Ampicillin sulbactam_I", 
                         "Aztreonam avibactam_I", "Cefixime_I", "Ceftaroline avibactam_I", 
                         "Ceftaroline_I", "Ceftazidime avibactam_I", "Ciprofloxacin_I", 
                         "Daptomycin_I", "Colistin_I", "Ertapenem_I", "Doripenem_I", 
                         "Gentamicin_I", "Gatifloxacin_I", "Moxifloxacin_I", "Oxacillin_I", 
                         "Quinupristin dalfopristin_I", "Sulbactam_I", "Tetracycline_I", 
                         "Teicoplanin_I", "Trimethoprim sulfa_I", "Ceftolozane tazobactam_I", 
                         "Meropenem vaborbactam_I", "Cefoperazone sulbactam_I", "Cefpodoxime_I", 
                         "Ceftibuten_I", "Tebipenem_I", "Ceftibuten avibactam_I")

# Loop through each variable and replace values for resistance and susceptible (0,1)
for (var in antibiotics_to_modify) {
  data_atlas_us[[var]] <- ifelse(data_atlas_us[[var]] == "Susceptible", 0, 
                                 ifelse(data_atlas_us[[var]] %in% c("Resistant", "Intermediate"), 1, data_atlas_us[[var]])
  )
}
for (var in antibiotics_to_modify) {
  data_atlas_eu[[var]] <- ifelse(data_atlas_eu[[var]] == "Susceptible", 0, 
                                 ifelse(data_atlas_eu[[var]] %in% c("Resistant", "Intermediate"), 1, data_atlas_eu[[var]])
  )
}

#Creating subgroup definitions for carbapenems and third-generation cephalosporins

carbapenems <- c("Imipenem_I", "Meropenem_I", "Ertapenem_I", "Doripenem_I", "Tebipenem_I", "Meropenem vaborbactam_I")
third_gen_cephalosporins <- c("Cefixime_I", "Ceftriaxone_I", "Ceftazidime_I", "Cefpodoxime_I", 
                              "Ceftibuten_I", "Ceftibuten avibactam_I", "Ceftazidime avibactam_I")

first_line_antibiotics <- c("Amikacin_I", "Ampicillin sulbactam_I", 
                            "Cefepime_I", "Piperacillin tazobactam_I",
                            "Ceftriaxone_I", "Ceftazidime avibactam_I", "Ciprofloxacin_I", 
                            "Gentamicin_I", "Ertapenem_I", "Meropenem_I", 
                            "Meropenem vaborbactam_I", "Trimethoprim sulfa_I", "Minocycline_I", "Amoxycillin clavulanate_I")
#https://www.idsociety.org/practice-guideline/amr-guidance/#Table1andTable2
  
#Below I am creating the AWaRe classification depending on antibiotic tested.
# Lists for AWARE classification based on Mike Thorne's list.
access <- c("Amikacin_I", "Ampicillin_I", "Amoxycillin clavulanate_I", "Cefoxitin_I", 
            "Clindamycin_I", "Cefazolin_I", "Cefalexin_I", "Penicillin_I", 
            "Metronidazole_I", "Gentamicin_I", "Trimethoprim sulfa_I", "Tetracycline_I", "Sulbactam_I","Oxacillin_I","Ampicillin sulbactam_I")

watch <- c("Azithromycin_I", "Ceftriaxone_I", "Ciprofloxacin_I", "Levofloxacin_I", 
           "Imipenem_I", "Meropenem_I", "Vancomycin_I", "Piperacillin tazobactam_I", 
           "Cefepime_I", "Aztreonam_I", "Clarithromycin_I", "Erythromycin_I", 
           "Minocycline_I", "Moxifloxacin_I", "Ceftazidime_I", 
           "Teicoplanin_I", "Daptomycin_I", "Ceftaroline_I", "Tebipenem_I", "Ceftibuten_I", "Gatifloxacin_I","Ertapenem_I","Cefixime_I")

reserve <- c("Linezolid_I", "Aztreonam avibactam_I", "Ceftaroline avibactam_I", "Ceftazidime avibactam_I", 
             "Ceftolozane tazobactam_I", "Meropenem vaborbactam_I", "Cefoperazone sulbactam_I", 
             "Doripenem_I", "Quinupristin dalfopristin_I", "Colistin_I", "Ceftibuten avibactam_I", "Tigecycline_I")

antibiotic_classes <- list(
  beta_lactams = c("Ampicillin_I", "Cefepime_I", "Amoxycillin clavulanate_I", "Ceftazidime_I", 
                   "Cefoxitin_I", "Ceftriaxone_I", "Penicillin_I", "Piperacillin tazobactam_I", 
                   "Aztreonam_I", "Ampicillin sulbactam_I", "Aztreonam avibactam_I", "Cefixime_I", 
                   "Ceftaroline avibactam_I", "Ceftaroline_I", "Ceftazidime avibactam_I", 
                   "Ceftolozane tazobactam_I", "Meropenem vaborbactam_I", "Cefoperazone sulbactam_I", 
                   "Cefpodoxime_I", "Ceftibuten_I", "Tebipenem_I", "Ceftibuten avibactam_I"),
  carbapenems = c("Imipenem_I", "Meropenem_I", "Ertapenem_I", "Doripenem_I"),
  aminoglycosides = c("Amikacin_I", "Gentamicin_I"),
  fluoroquinolones = c("Levofloxacin_I", "Ciprofloxacin_I", "Gatifloxacin_I", "Moxifloxacin_I"),
  macrolides = c("Azithromycin_I", "Clarithromycin_I", "Erythromycin_I"),
  glycopeptides = c("Vancomycin_I", "Teicoplanin_I"),
  lipopeptides = c("Daptomycin_I"),
  oxazolidinones = c("Linezolid_I"),
  tetracyclines = c("Minocycline_I", "Tetracycline_I"),
  polymyxins = c("Colistin_I"),
  lincosamides = c("Clindamycin_I"),
  streptogramins = c("Quinupristin dalfopristin_I"),
  nitroimidazoles = c("Metronidazole_I"),
  sulfonamides = c("Trimethoprim sulfa_I"),
  glycylcyclines = c("Tigecycline_I"))

# Function to classify antibiotics based on WHO AWaRe classification
classify_aware <- function(antibiotic) {
  if (antibiotic %in% access) {
    return("Access")
  } else if (antibiotic %in% watch) {
    return("Watch")
  } else if (antibiotic %in% reserve) {
    return("Reserve")
  } else {
    return("Unclassified/Not Recommended")
  }
}
df_AtbAware <- data.frame(Antibiotic = antibiotics_to_modify, stringsAsFactors = FALSE)

# Apply classification to each antibiotic
df_AtbAware$AWaRe_classification <- sapply(df_AtbAware$Antibiotic, classify_aware)
# Check the dataframe
print(df_AtbAware)

data_atlas_us_orig<- data_atlas_us
data_atlas_eu_orig<- data_atlas_eu
#Only keep BSIs and Urine samples
data_atlas_us <- subset(data_atlas_us, Source %in% c("Blood", "Blood Vessels", "Urine"))
data_atlas_us$Source[data_atlas_us$Source == "Blood Vessels"] <- "Blood"

data_atlas_eu <- subset(data_atlas_eu, Source %in% c("Blood", "Blood Vessels", "Urine"))
data_atlas_eu$Source[data_atlas_eu$Source == "Blood Vessels"] <- "Blood"
# 415582 all // 125340 bsi and uti // 64540.
#######

#------------------------------------------------------------------------#
#Creating sub-datasets for pathogens described in WHO list: https://iris.who.int/bitstream/handle/10665/376776/9789240093461-eng.pdf?sequence=1 #######
#Gram-negative
data_atlas_us_enterob <- subset(data_atlas_us, Family %in% c("Enterobacteriaceae", "Enterobacterales"))
data_atlas_us_enterob_orig <- subset(data_atlas_us_orig, Family %in% c("Enterobacteriaceae", "Enterobacterales"))
data_atlas_us_abaumannii <- subset(data_atlas_us, Species == "Acinetobacter baumannii")
data_atlas_us_pseudomonas <- subset(data_atlas_us, Species == "Pseudomonas aeruginosa")

data_atlas_eu_enterob <- subset(data_atlas_eu, Family %in% c("Enterobacteriaceae", "Enterobacterales"))
data_atlas_eu_enterob_orig <- subset(data_atlas_eu_orig, Family %in% c("Enterobacteriaceae", "Enterobacterales"))
data_atlas_eu_abaumannii <- subset(data_atlas_eu, Species == "Acinetobacter baumannii")
data_atlas_eu_pseudomonas <- subset(data_atlas_eu, Species == "Pseudomonas aeruginosa")


#Gram-positive
data_atlas_us_enterococ <- subset(data_atlas_us, Species == "Enterococcus faecium")
data_atlas_us_staphy <- subset(data_atlas_us, Family == "Staphylococcus spp")
data_atlas_us_streptop <- subset(data_atlas_us, Family == "Streptococcus pneumoniae")

data_atlas_eu_enterococ <- subset(data_atlas_eu, Species == "Enterococcus faecium")
data_atlas_eu_staphy <- subset(data_atlas_eu, Family == "Staphylococcus spp")
data_atlas_eu_streptop <- subset(data_atlas_eu, Family == "Streptococcus pneumoniae")


#Group up by Gram-positive and negative.
data_atlas_us_gramneg <- subset(data_atlas_us, Family == "Enterobacterales" | Family == "Enterobacteriaceae"| Species == "Acinetobacter baumannii" | Species == "Pseudomonas aeruginosa")
data_atlas_us_grampos <- subset(data_atlas_us, Family == "Staphylococcus spp"| Species == "Enterococcus faecium" | Species == "Streptococcus pneumoniae")

data_atlas_eu_gramneg <- subset(data_atlas_eu, Family == "Enterobacterales" | Family == "Enterobacteriaceae"| Species == "Acinetobacter baumannii" | Species == "Pseudomonas aeruginosa")
data_atlas_eu_grampos <- subset(data_atlas_eu, Family == "Staphylococcus spp"| Species == "Enterococcus faecium" | Species == "Streptococcus pneumoniae")


#Group it all.
data_atlas_us_allf <- rbind(data_atlas_us_gramneg, data_atlas_us_grampos)
# Create a new column 'Family2' with NAs
data_atlas_us_allf$Family2 <- NA
data_atlas_us_allf$Family2[data_atlas_us_allf$Family == "Staphylococcus spp"] <- "Staphylococcus aureus"
data_atlas_us_allf$Family2[data_atlas_us_allf$Family == "Enterobacterales"] <- "Enterobacterales"
data_atlas_us_allf$Family2[data_atlas_us_allf$Family == "Enterobacteriaceae"] <- "Enterobacterales"
data_atlas_us_allf$Family2[data_atlas_us_allf$Species == "Acinetobacter baumannii"] <- "Acinetobacter baumannii"
data_atlas_us_allf$Family2[data_atlas_us_allf$Species == "Pseudomonas aeruginosa"] <- "Pseudomonas aeruginosa"
data_atlas_us_allf$Family2[data_atlas_us_allf$Species == "Streptococcus pneumoniae"] <- "Streptococcus pneumoniae"
data_atlas_us_allf$Family2[data_atlas_us_allf$Species == "Enterococcus faecium"] <- "Enterococcus faecium"

data_atlas_eu_allf <- rbind(data_atlas_eu_gramneg, data_atlas_eu_grampos)
# Create a new column 'Family2' with NAs
data_atlas_eu_allf$Family2 <- NA
data_atlas_eu_allf$Family2[data_atlas_eu_allf$Family == "Staphylococcus spp"] <- "Staphylococcus aureus"
data_atlas_eu_allf$Family2[data_atlas_eu_allf$Family == "Enterobacterales"] <- "Enterobacterales"
data_atlas_eu_allf$Family2[data_atlas_eu_allf$Family == "Enterobacteriaceae"] <- "Enterobacterales"
data_atlas_eu_allf$Family2[data_atlas_eu_allf$Species == "Acinetobacter baumannii"] <- "Acinetobacter baumannii"
data_atlas_eu_allf$Family2[data_atlas_eu_allf$Species == "Pseudomonas aeruginosa"] <- "Pseudomonas aeruginosa"
data_atlas_eu_allf$Family2[data_atlas_eu_allf$Species == "Streptococcus pneumoniae"] <- "Streptococcus pneumoniae"
data_atlas_eu_allf$Family2[data_atlas_eu_allf$Species == "Enterococcus faecium"] <- "Enterococcus faecium"


#Group up by critical, high-group, or medium-group WHO bacterial priority pathogens
data_atlas_us_critical <- subset(data_atlas_us, Family == "Enterobacterales" | Family == "Enterobacteriaceae"| Species == "Acinetobacter baumannii")
data_atlas_us_high <- subset(data_atlas_us, Family == "Staphylococcus spp" | Species == "Enterococcus faecium"| Species == "Pseudomonas aeruginosa")
data_atlas_us_medium<- subset(data_atlas_us, Family == "Streptococcus pneumoniae")

data_atlas_eu_critical <- subset(data_atlas_eu, Family == "Enterobacterales" | Family == "Enterobacteriaceae"| Species == "Acinetobacter baumannii")
data_atlas_eu_high <- subset(data_atlas_eu, Family == "Staphylococcus spp" | Species == "Enterococcus faecium"| Species == "Pseudomonas aeruginosa")
data_atlas_eu_medium<- subset(data_atlas_eu, Family == "Streptococcus pneumoniae")
#######

####### See how many cultures there are per year in the US and Europe, by source #######
# Ensure the Year column is a factor
data_atlas_us_allf$Year <- as.factor(data_atlas_us_allf$Year)
# Calculate the count of each source per year
source_count_by_year <- data_atlas_us_allf %>%
  group_by(Year, Source) %>%
  summarise(Count = n()) %>%
  ungroup()
# Calculate the proportion of each source per year
source_proportion_by_year <- source_count_by_year %>%
  group_by(Year) %>%
  mutate(Total_Count = sum(Count),
         Proportion = (Count / Total_Count) * 100) %>%
  ungroup()
# Create labels with the year and number of observations
source_proportion_by_year <- source_proportion_by_year %>%
  group_by(Year) %>%
  mutate(Label = paste0(Year, " (N = ", Total_Count[1], ")")) %>%
  ungroup()

#  color palette
lancet_palette <- c("Blood" = "#D55E00", "Urine" = "#F39C12")  # Adjusted for Lancet style

# Create pie charts for each year
USA_samples<-ggplot(source_proportion_by_year, aes(x = "", y = Proportion, fill = Source)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +  # Black contour
  coord_polar(theta = "y") +
  facet_wrap(~ Label) +  # Use the new Label with year and (N = xxx)
  theme_void() +  # Remove x and y axis
  scale_fill_manual(values = lancet_palette) +  # Lancet-style colors
  geom_text(aes(label = paste0(round(Proportion, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 3, color = "white") +  # Add percentages inside the pie
  labs(title = "Proportion of blood and urine samples by year", fill = "Source") +
  theme(
    plot.title = element_text(hjust = 0, size = 14),  # Align title to the left
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.position = "right"
  )
ggsave(filename = "SamplesBloodUrine_US.tiff", plot = USA_samples, device = "tiff", 
       path = base_pathOut,
       width = 11, height = 8, dpi = 500, units = "in") 

data_atlas_eu_allf$Year <- as.factor(data_atlas_eu_allf$Year)
# Calculate the count of each source per year
source_count_by_year <- data_atlas_eu_allf %>%
  group_by(Year, Source) %>%
  summarise(Count = n()) %>%
  ungroup()
# Calculate the proportion of each source per year
source_proportion_by_year <- source_count_by_year %>%
  group_by(Year) %>%
  mutate(Total_Count = sum(Count),
         Proportion = (Count / Total_Count) * 100) %>%
  ungroup()
# Create labels with the year and number of observations
source_proportion_by_year <- source_proportion_by_year %>%
  group_by(Year) %>%
  mutate(Label = paste0(Year, " (N = ", Total_Count[1], ")")) %>%
  ungroup()
# Create pie charts for each year
EU_samples<-ggplot(source_proportion_by_year, aes(x = "", y = Proportion, fill = Source)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +  # Black contour
  coord_polar(theta = "y") +
  facet_wrap(~ Label) +  # Use the new Label with year and (N = xxx)
  theme_void() +  # Remove x and y axis
  scale_fill_manual(values = lancet_palette) +  # Lancet-style colors
  geom_text(aes(label = paste0(round(Proportion, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 3, color = "white") +  # Add percentages inside the pie
  labs(title = "Proportion of blood and urine samples by year", fill = "Source") +
  theme(
    plot.title = element_text(hjust = 0, size = 14),  # Align title to the left
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.position = "right"
  )
ggsave(filename = "SamplesBloodUrine_EU.tiff", plot = EU_samples, device = "tiff", 
       path = base_pathOut,
       width = 11, height = 8, dpi = 500, units = "in") 


#######

######Table of descriptives and image pathogen/quantity over time#######

data_atlas_us_enterob_o <- subset(data_atlas_us_orig, Family %in% c("Enterobacteriaceae", "Enterobacterales"))

family2_year_summary2 <- data_atlas_us_enterob_o %>%
  group_by(State,Year) %>%
  summarise(Observations = n()) %>%
  arrange(State, Year)

family2_year_summary2 <- data_atlas_us_enterob %>%
  group_by(State,Year) %>%
  summarise(Observations = n()) %>%
  arrange(State, Year)

# Create the plot with text labels for observations
Cultur_bsiuti_stat_time<-ggplot(family2_year_summary2, aes(x = Year, y = State, fill = Observations)) +
  geom_tile(color = "white") +  # Tiles with white borders
  geom_text(aes(label = Observations), color = "black", size = 3) +  # Add text labels for observations
  scale_fill_gradient(name = "Observations", low = "#F4CCCC", high = "#990000") +  # Lancet-style gradient
  scale_x_continuous(breaks = seq(2004, 2022, by = 1)) +  # Set x-axis to display each year from 2004 to 2022
  theme_minimal() +  # Minimal theme for a clean look
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.title = element_text(size = 12),  # Adjust axis title size
    axis.text = element_text(size = 10),  # Adjust axis text size
    plot.title = element_text(size = 14, face = "bold")  # Adjust plot title size and bold
  ) +
  labs(
    title = "Number of observations per state by year",
    x = "Year",
    y = "US State"
  )
ggsave(filename = "Cultur_bsiuti_stat_time.tiff", plot = Cultur_bsiuti_stat_time, device = "tiff", 
       path = base_pathOut,
       width = 11, height = 8, dpi = 500, units = "in") 


# Create a summary table of observations by Family and Year
family2_year_summary <- data_atlas_us_allf %>%
  group_by(Family2, Year) %>%
  summarise(Observations = n()) %>%
  arrange(Family2, Year)
# View the summary table
print(family2_year_summary)

# Convert the data to a wide format for the heatmap
data_wide <- dcast(family2_year_summary, Family2 ~ Year, value.var = "Observations", fill = 0)
# Melt the data back into a long format for ggplot
data_long <- melt(data_wide, id.vars = "Family2", variable.name = "Year", value.name = "Observations")
# Convert Year to a factor for proper ordering in the plot
data_long$Year <- factor(data_long$Year, levels = sort(unique(data_long$Year)))
# Define colors 
red_palette <- c("#FFCCCC", "#FF6666", "#FF0000", "#CC0000", "#990000")
# Create the heatmap
PathogenAvaiTime<-ggplot(data_long, aes(x = Year, y = Family2, fill = Observations)) +
  geom_tile(color = "black") +  # Add black borders around the tiles
  scale_fill_gradientn(colors = red_palette, na.value = "grey80") +  # Use Lancet-style colors
  geom_text(aes(label = Observations), color = "white", size = 3) +  # Add the observation numbers on the tiles
  labs(title = "Number of observations by pathogen and year", fill = "Observations") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 360, hjust = 1),
    plot.title = element_text(hjust = 0, face = "bold"),
    legend.position = "right"
  )

ggsave(filename = "PathogensAvailableTime_US.tiff", plot = PathogenAvaiTime, device = "tiff", 
       path = base_pathOut,
       width = 11, height = 8, dpi = 500, units = "in") 

######

#####Descriptives per state, count of isolates per state 'total 2004-2021'#####
# Calculate the counts per state
state_counts <- data_atlas_us_enterob %>%
  group_by(State) %>%
  summarise(Count = n())

# Get US map data
us_map <- map_data("state")
# Convert state names in the counts data to lowercase for merging
state_counts$State <- tolower(state_counts$State)
# Merge the map data with the state counts
us_map_data <- left_join(us_map, state_counts, by = c("region" = "State"))
# Replace NA counts with 0 for states with no data
us_map_data$Count[is.na(us_map_data$Count)] <- 0


# Plot the map
Count_isolatespState<-ggplot(us_map_data, aes(x = long, y = lat, group = group, fill = Count)) +
  geom_polygon(color = "black", size = 0.2) +  # Draw state borders
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey80") +  # Lancet-style colors
  labs(title = "",
       fill = "Total count\n of isolates\n2004-2021") +
  theme_void() +  # Removes axis titles and ticks
  theme(
    plot.title = element_text(hjust = 0, size = 14, face = "bold"),  # Align title to the left
    legend.position = "right",  # Position the legend on the right
    legend.key.height = unit(2, "cm"),  # Adjust the height of the legend bar
    legend.key.width = unit(0.5, "cm"),  # Adjust the width of the legend bar
    legend.title = element_text(size = 12, face = "bold", hjust = 0),  # Adjust legend title and align left
    legend.text = element_text(size = 10),  # Adjust legend text
    legend.background = element_rect(color = "white", fill = "white"),  # Black contour around the legend bar
    legend.margin = margin(0, 20, 0, 0),  # Add more margin between the legend and the right side of the image
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Adjust plot margins
  ) 

ggsave(filename = "Count_isolatespState.tiff", plot = Count_isolatespState, device = "tiff", 
       path = base_pathOut,
       width = 12, height = 8, dpi = 500, units = "in") 


country_counts <- data_atlas_eu_enterob %>%
  group_by(Country) %>%
  summarise(Count = n())

#####


#------------------------------------------------------------------------#
#FOR THE US
#Subgroup analyses for enterobacterales carbapenems ##########
#data_atlas_us_enterob_carbap <- data_atlas_us_enterob_carbap %>%

variables_to_keep <- c("Year", carbapenems)
data_atlas_us_enterob_carbap <- data_atlas_us_enterob[, variables_to_keep]


data_long <- melt(data_atlas_us_enterob_carbap, id.vars = "Year", measure.vars = carbapenems, 
                  variable.name = "Antibiotic", value.name = "Count")
# Summarize the data to get counts of non-missing values per antibiotic per year
data_summary <- dcast(data_long, Antibiotic ~ Year, fun.aggregate = function(x) sum(!is.na(x)))
# View the matrix
data_summary

#Key carbapenem-specific resistance for Erta/Mero/Imi in the US: Imipenem_I Meropenem_I
data_atlas_us_enterob_carbap <- subset(data_atlas_us_enterob, select = c("Year", "State","Gender","Age Group","Species", carbapenems))

# Calculate the number of missing values for each carbapenem
missing_data_summary <- colSums(is.na(data_atlas_us_enterob_carbap[, carbapenems]))
# Calculate the proportion of missing data for each carbapenem
total_rows <- nrow(data_atlas_us_enterob_carbap)
missing_data_proportion <- missing_data_summary / total_rows * 100
# Combine the results into a data frame for easy viewing
missing_data_report <- data.frame(
  Antibiotic = carbapenems,
  Missing_Count = missing_data_summary,
  Missing_Proportion = missing_data_proportion)
# Display table
print(missing_data_report)

#Carbapenems all
data_atlas_us_enterob_carbap[carbapenems] <- lapply(data_atlas_us_enterob_carbap[carbapenems], as.numeric)
data_atlas_us_enterob_carbap$Carbapenem_Sum <- rowSums(data_atlas_us_enterob_carbap[, carbapenems], na.rm = TRUE)
data_atlas_us_enterob_carbap$Carbapenem_res <- ifelse(data_atlas_us_enterob_carbap$Carbapenem_Sum >= 1, 1, 0)
data_atlas_us_enterob_carbap$Carbapenem_res  <- as.numeric(data_atlas_us_enterob_carbap$Carbapenem_res )
# Step 1: Calculate Meropenem resistance by year across the US
carbapenem_resistance_by_year <- aggregate(Carbapenem_res  ~ Year, data = data_atlas_us_enterob_carbap, FUN = mean, na.rm = TRUE)
carbapenem_resistance_by_year$Carbapenem_res <- carbapenem_resistance_by_year$Carbapenem_res *100
carbapenem_resistance_by_year$Year <- as.factor(carbapenem_resistance_by_year$Year)

# Create the plot
Carbapenem_trendUS<-ggplot(carbapenem_resistance_by_year, aes(x = Year, y = Carbapenem_res, group = 1)) +  # Set group = 1 to treat all data as a single group
  geom_line(color = "#D55E00", size = 1, linetype = "solid") +  
  geom_point(color = "black", size = 3, shape = 21, fill = "#D55E60", stroke = 0.5) + 
  labs(title = "",
       x = "Year",
       y = "Carbapenem resistance (%) in Enterobacterales") +
  theme_minimal() +  # Use a minimal theme for clarity
  theme(
    plot.title = element_text(hjust = 0, size = 14, face = "bold"),  # Align title to the left
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),   # Horizontal x-axis labels
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "gray80")   # Light gray grid lines for better contrast
  )

ggsave(filename = "Carbapenem_trendUS.tiff", plot = Carbapenem_trendUS, device = "tiff", 
       path = base_pathOut,
       width = 11, height = 8, dpi = 500, units = "in") 





#######

#Subgroup analyses for enterobacterales first-line antibitiocs ##########
variables_to_keep <- c("Year", first_line_antibiotics)
data_atlas_us_enterob_firstlin <- data_atlas_us_enterob[, variables_to_keep]


data_long <- melt(data_atlas_us_enterob_firstlin, id.vars = "Year", measure.vars = first_line_antibiotics, 
                  variable.name = "Antibiotic", value.name = "Count")
# Summarize the data to get counts of non-missing values per antibiotic per year
data_summary <- dcast(data_long, Antibiotic ~ Year, fun.aggregate = function(x) sum(!is.na(x)))
# View the matrix
data_summary

#Key firstline ATB resistance
data_atlas_us_enterob_firstlin <- subset(data_atlas_us_enterob, select = c("Year", "State","Gender","Age Group","Species", first_line_antibiotics))

# Calculate the number of missing values for each carbapenem
missing_data_summary <- colSums(is.na(data_atlas_us_enterob_firstlin[, first_line_antibiotics]))
# Calculate the proportion of missing data for each carbapenem
total_rows <- nrow(data_atlas_us_enterob_firstlin)
missing_data_proportion <- missing_data_summary / total_rows * 100
# Combine the results into a data frame for easy viewing
missing_data_report <- data.frame(
  Antibiotic = first_line_antibiotics,
  Missing_Count = missing_data_summary,
  Missing_Proportion = missing_data_proportion)
# Display table
print(missing_data_report)

#First_lineAtbs
data_atlas_us_enterob_firstlin[first_line_antibiotics] <- lapply(data_atlas_us_enterob_firstlin[first_line_antibiotics], as.numeric)
data_atlas_us_enterob_firstlin$Firstlin_Sum <- rowSums(data_atlas_us_enterob_firstlin[, first_line_antibiotics], na.rm = TRUE)
data_atlas_us_enterob_firstlin$Firstlin_res <- ifelse(data_atlas_us_enterob_firstlin$Firstlin_Sum >= 1, 1, 0)
data_atlas_us_enterob_firstlin$Firstlin_res  <- as.numeric(data_atlas_us_enterob_firstlin$Firstlin_res )
# Step 1: Calculate Meropenem resistance by year across the US
Firstlin_resistance_by_year <- aggregate(Firstlin_res  ~ Year, data = data_atlas_us_enterob_firstlin, FUN = mean, na.rm = TRUE)
Firstlin_resistance_by_year$Firstlin_res <- Firstlin_resistance_by_year$Firstlin_res *100
Firstlin_resistance_by_year$Year <- as.factor(Firstlin_resistance_by_year$Year)

# Create the plot
Firstlin_trendUS<-ggplot(Firstlin_resistance_by_year, aes(x = Year, y = Firstlin_res, group = 1)) +  # Set group = 1 to treat all data as a single group
  geom_line(color = "#D55E00", size = 1, linetype = "solid") +  
  geom_point(color = "black", size = 3, shape = 21, fill = "#D55E60", stroke = 0.5) + 
  labs(title = "",
       x = "Year",
       y = "First-line antibiotic resistance (%) in Enterobacterales") +
  theme_minimal() +  # Use a minimal theme for clarity
  theme(
    plot.title = element_text(hjust = 0, size = 14, face = "bold"),  # Align title to the left
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),   # Horizontal x-axis labels
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "gray80")   # Light gray grid lines for better contrast
  )

ggsave(filename = "Firstlin_trendUS.tiff", plot = Firstlin_trendUS, device = "tiff", 
       path = base_pathOut,
       width = 11, height = 8, dpi = 500, units = "in") 





#######

#Subgroup analyses Enterobacterales for Cephalosporins########
#Key 3rd Generation Cephalosporins : Cefixime_I Ceftriaxone_I Ceftazidime_I Cefpodoxime_I "Ceftazidime avibactam_I" 
variables_to_keep <- c("Year", third_gen_cephalosporins)
data_atlas_us_enterob_cepha <- data_atlas_us_enterob[, variables_to_keep]


data_long <- melt(data_atlas_us_enterob_cepha, id.vars = "Year", measure.vars = third_gen_cephalosporins, 
                  variable.name = "Antibiotic", value.name = "Count")
# Summarize the data to get counts of non-missing values per antibiotic per year
data_summary <- dcast(data_long, Antibiotic ~ Year, fun.aggregate = function(x) sum(!is.na(x)))
# View the matrix
data_summary


data_atlas_us_enterob_cepha<- subset(data_atlas_us_enterob, select = c("Year", "State","Gender","Age Group","Species", third_gen_cephalosporins))
missing_data_summary <- colSums(is.na(data_atlas_us_enterob_cepha[, third_gen_cephalosporins]))
# Calculate the proportion of missing data for each cephalosporings
total_rows <- nrow(data_atlas_us_enterob_cepha)
missing_data_proportion <- missing_data_summary / total_rows * 100
# Combine the results into a data frame for easy viewing
missing_data_report <- data.frame(
  Antibiotic = third_gen_cephalosporins,
  Missing_Count = missing_data_summary,
  Missing_Proportion = missing_data_proportion)
# Display the report
print(missing_data_report)

#3rd Generation Cephalosporins ALL ---- #
data_atlas_us_enterob_cepha[third_gen_cephalosporins] <- lapply(data_atlas_us_enterob_cepha[third_gen_cephalosporins], as.numeric)
data_atlas_us_enterob_cepha$Cephalosporin_Sum <- rowSums(data_atlas_us_enterob_cepha[, third_gen_cephalosporins], na.rm = TRUE)

data_atlas_us_enterob_cepha$Cephalosporin_res <- ifelse(data_atlas_us_enterob_cepha$Cephalosporin_Sum >= 1, 1, 0)
data_atlas_us_enterob_cepha$Cephalosporin_res  <- as.numeric(data_atlas_us_enterob_cepha$Cephalosporin_res )
# Step 1: Calculate Meropenem resistance by year across the US
cephalosporin_resistance_by_year <- aggregate(Cephalosporin_res  ~ Year, data = data_atlas_us_enterob_cepha, FUN = mean, na.rm = TRUE)
cephalosporin_resistance_by_year$Cephalosporin_res <- cephalosporin_resistance_by_year$Cephalosporin_res *100
cephalosporin_resistance_by_year$Year <- as.factor(cephalosporin_resistance_by_year$Year)

# Create the plot
Cephalosporin_trendUS<-ggplot(cephalosporin_resistance_by_year, aes(x = Year, y = Cephalosporin_res, group = 1)) +  # Set group = 1 to treat all data as a single group
  geom_line(color = "#D55E00", size = 1, linetype = "solid") +  
  geom_point(color = "black", size = 3, shape = 21, fill = "#D55E60", stroke = 0.5) + 
  labs(title = "",
       x = "Year",
       y = "Third generation cephalosporin resistance (%) in Enterobacterales") +
  theme_minimal() +  # Use a minimal theme for clarity
  theme(
    plot.title = element_text(hjust = 0, size = 14, face = "bold"),  # Align title to the left
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),   # Horizontal x-axis labels
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "gray80")   # Light gray grid lines for better contrast
  )

ggsave(filename = "Cephalosporin_trendUS.tiff", plot = Cephalosporin_trendUS, device = "tiff", 
       path = base_pathOut,
       width = 11, height = 8, dpi = 500, units = "in") 

######


##Subgroup analyses Enterobacterales MDR (at least 3 antibiotic classes) ####

# Convert all antibiotic variables to numeric
antibiotic_list <- unlist(antibiotic_classes)  # Flatten the list of antibiotic classes into a vector of antibiotic names
# Convert each antibiotic variable in the dataframe to numeric
for (antibiotic in antibiotic_list) {
  data_atlas_us[[antibiotic]] <- as.numeric(data_atlas_us[[antibiotic]])
}

for (class_name in names(antibiotic_classes)) {
  antibiotic_list <- antibiotic_classes[[class_name]]
  data_atlas_us[[class_name]] <- rowSums(data_atlas_us[, antibiotic_list], na.rm = TRUE)
}
for (class_name in names(antibiotic_classes)) {
  data_atlas_us[[class_name]] <- ifelse(data_atlas_us[[class_name]] >= 1, 1, 0)
}
# View the updated dataframe with the new variables
head(data_atlas_us)

#Creating MDR variable:
# Calculate the sum of the antibiotic class variables
data_atlas_us$Total_Antibiotic_Class_Sum <- rowSums(data_atlas_us[, c("beta_lactams","carbapenems", "aminoglycosides",  "fluoroquinolones",  "macrolides", "glycopeptides", "lipopeptides", "oxazolidinones", "tetracyclines", "polymyxins", "lincosamides" ,"streptogramins", "nitroimidazoles",  "sulfonamides", "glycylcyclines")], na.rm = TRUE)
# Create the MDR variable based on the sum
data_atlas_us$MDR <- ifelse(data_atlas_us$Total_Antibiotic_Class_Sum >= 3, 1, 0)


data_atlas_us$MDR  <- as.numeric(data_atlas_us$MDR )

MDR_resistance_by_year <- aggregate(MDR  ~ Year, data = data_atlas_us, FUN = mean, na.rm = TRUE)
MDR_resistance_by_year$MDR <- MDR_resistance_by_year$MDR *100
MDR_resistance_by_year$Year <- as.factor(MDR_resistance_by_year$Year)

# Create the plot
MDR_trendUS<-ggplot(MDR_resistance_by_year, aes(x = Year, y = MDR, group = 1)) +  # Set group = 1 to treat all data as a single group
  geom_line(color = "#D55E00", size = 1, linetype = "solid") +  
  geom_point(color = "black", size = 3, shape = 21, fill = "#D55E60", stroke = 0.5) + 
  labs(title = "",
       x = "Year",
       y = "Multidrug resistance (%) in Enterobacterales") +
  theme_minimal() +  # Use a minimal theme for clarity
  theme(
    plot.title = element_text(hjust = 0, size = 14, face = "bold"),  # Align title to the left
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),   # Horizontal x-axis labels
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "gray80")   # Light gray grid lines for better contrast
  )

ggsave(filename = "MDR_trendUS.tiff", plot = MDR_trendUS, device = "tiff", 
       path = base_pathOut,
       width = 11, height = 8, dpi = 500, units = "in") 


#####


#4Maps with first-line resistance, carbapenem, cephalosporin and MDR resistance per state#####
MDR_resistance_by_year_state <- aggregate(MDR ~ State, data = data_atlas_us, FUN = mean, na.rm = TRUE)
MDR_resistance_by_year_state$MDR <- MDR_resistance_by_year_state$MDR *100
#MDR_resistance_by_year_state$Year <- as.factor(MDR_resistance_by_year_state$Year)

cephalosporin_resistance_by_year_state <- aggregate(Cephalosporin_res  ~  State, data = data_atlas_us_enterob_cepha, FUN = mean, na.rm = TRUE)
cephalosporin_resistance_by_year_state$Cephalosporin_res <- cephalosporin_resistance_by_year_state$Cephalosporin_res *100
#cephalosporin_resistance_by_year_state$Year <- as.factor(cephalosporin_resistance_by_year_state$Year)

carbapenem_resistance_by_year_state <- aggregate(Carbapenem_res  ~ State, data = data_atlas_us_enterob_carbap, FUN = mean, na.rm = TRUE)
carbapenem_resistance_by_year_state$Carbapenem_res <- carbapenem_resistance_by_year_state$Carbapenem_res *100
#carbapenem_resistance_by_year_state$Year <- as.factor(carbapenem_resistance_by_year_state$Year)

Firstlin_resistance_by_year_state <- aggregate(Firstlin_res  ~ State, data = data_atlas_us_enterob_firstlin, FUN = mean, na.rm = TRUE)
Firstlin_resistance_by_year_state$Firstlin_res <- Firstlin_resistance_by_year_state$Firstlin_res *100


# Load US map data
us_map <- map_data("state")
# Prepare the data: ensure state names match and join with map data
carbapenem_resistance_by_year_state$region <- tolower(carbapenem_resistance_by_year_state$State)
cephalosporin_resistance_by_year_state$region <- tolower(cephalosporin_resistance_by_year_state$State)
MDR_resistance_by_year_state$region <- tolower(MDR_resistance_by_year_state$State)
Firstlin_resistance_by_year_state$region <- tolower(Firstlin_resistance_by_year_state$State)

us_map_carbapenem <- us_map %>%
  left_join(carbapenem_resistance_by_year_state, by = "region")
us_map_cephalosporin <- us_map %>%
  left_join(cephalosporin_resistance_by_year_state, by = "region")
us_map_MDR <- us_map %>%
  left_join(MDR_resistance_by_year_state, by = "region")
us_map_firstlin <- us_map %>%
  left_join(Firstlin_resistance_by_year_state, by = "region")

lancet_theme <- theme(
  text = element_text(family = "serif", color = "black"),
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  panel.background = element_rect(fill = "white", color = "white"),
  legend.position = "bottom",
  legend.title = element_text(size = 10, face = "bold"),
  legend.text = element_text(size = 8)
)

#Individual graphs:

# Plot for Carbapenem Resistance
carbapenem_heatmap <- ggplot(us_map_carbapenem, aes(x = long, y = lat, group = group, fill = Carbapenem_res)) +
  geom_polygon(color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey80") +
  labs(title = "B. Carbapenem resistance by US state", fill = "Carbapenem\nresistance (%)") +
  theme_void() +
  lancet_theme +
  theme(plot.title = element_text(hjust = 0, size = 12, face = "bold"))

# Plot for Cephalosporin Resistance
cephalosporin_heatmap <- ggplot(us_map_cephalosporin, aes(x = long, y = lat, group = group, fill = Cephalosporin_res)) +
  geom_polygon(color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey80") +
  labs(title = "C. Cephalosporin resistance by US state", fill = "Cephalosporin\nresistance (%)") +
  theme_void() +
  lancet_theme +
  theme(plot.title = element_text(hjust = 0, size = 12, face = "bold"))

# Plot for MDR Resistance
MDR_heatmap <- ggplot(us_map_MDR, aes(x = long, y = lat, group = group, fill = MDR)) +
  geom_polygon(color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey80") +
  labs(title = "D. MDR resistance by US state", fill = "MDR (%)") +
  theme_void() +
  lancet_theme +
  theme(plot.title = element_text(hjust = 0, size = 12, face = "bold"))

# Firstline Resistance
Firstlin_heatmap <- ggplot(us_map_firstlin, aes(x = long, y = lat, group = group, fill = Firstlin_res)) +
  geom_polygon(color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey80") +
  labs(title = "A. First-line AMR by US state", fill = "First-line antibiotic\n resistance (%)") +
  theme_void() +
  lancet_theme +
  theme(plot.title = element_text(hjust = 0, size = 12, face = "bold"))




library(gridExtra)
# Combine the three heatmaps
#combined_plot <- grid.arrange(carbapenem_heatmap, cephalosporin_heatmap, MDR_heatmap, ncol = 1)
# Save the combined plot
#ggsave(filename = "resistance_by_state_combined.tiff", plot = combined_plot, device = "tiff", path = base_pathOut,
#       width = 8, height = 13, dpi = 500, units = "in")


######

#--EU-----------------------------------------------------------------------#

#Subgroup analyses for enterobacterales carbapenems ##########

variables_to_keep <- c("Year", carbapenems, "Country", "Gender", "Species")
data_atlas_eu_enterob_carbap <- data_atlas_eu_enterob[, variables_to_keep]


data_long <- melt(data_atlas_eu_enterob_carbap, id.vars = "Year", measure.vars = carbapenems, 
                  variable.name = "Antibiotic", value.name = "Count")
# Summarize the data to get counts of non-missing values per antibiotic per year
data_summary <- dcast(data_long, Antibiotic ~ Year, fun.aggregate = function(x) sum(!is.na(x)))
# View the matrix
data_summary

#Key carbapenem-specific resistance for Erta/Mero/Imi in the US: Imipenem_I Meropenem_I
#data_atlas_eu_enterob_carbap <- subset(data_atlas_eu_enterob, select = c("Year","Country","Gender","Species", carbapenems))

# Calculate the number of missing values for each carbapenem
missing_data_summary <- colSums(is.na(data_atlas_eu_enterob_carbap[, carbapenems]))
# Calculate the proportion of missing data for each carbapenem
total_rows <- nrow(data_atlas_eu_enterob_carbap)
missing_data_proportion <- missing_data_summary / total_rows * 100
# Combine the results into a data frame for easy viewing
missing_data_report <- data.frame(
  Antibiotic = carbapenems,
  Missing_Count = missing_data_summary,
  Missing_Proportion = missing_data_proportion)
# Display table
print(missing_data_report)

#Carbapenems all
data_atlas_eu_enterob_carbap[carbapenems] <- lapply(data_atlas_eu_enterob_carbap[carbapenems], as.numeric)
data_atlas_eu_enterob_carbap$Carbapenem_Sum <- rowSums(data_atlas_eu_enterob_carbap[, carbapenems], na.rm = TRUE)
data_atlas_eu_enterob_carbap$Carbapenem_res <- ifelse(data_atlas_eu_enterob_carbap$Carbapenem_Sum >= 1, 1, 0)
data_atlas_eu_enterob_carbap$Carbapenem_res  <- as.numeric(data_atlas_eu_enterob_carbap$Carbapenem_res )
# Step 1: Calculate Meropenem resistance by year across the US
carbapenem_resistance_by_year <- aggregate(Carbapenem_res  ~ Year, data = data_atlas_eu_enterob_carbap, FUN = mean, na.rm = TRUE)
carbapenem_resistance_by_year$Carbapenem_res <- carbapenem_resistance_by_year$Carbapenem_res *100
carbapenem_resistance_by_year$Year <- as.factor(carbapenem_resistance_by_year$Year)

# Create the plot
Carbapenem_trendEU<-ggplot(carbapenem_resistance_by_year, aes(x = Year, y = Carbapenem_res, group = 1)) +  # Set group = 1 to treat all data as a single group
  geom_line(color = "#D55E00", size = 1, linetype = "solid") +  
  geom_point(color = "black", size = 3, shape = 21, fill = "#D55E60", stroke = 0.5) + 
  labs(title = "",
       x = "Year",
       y = "Carbapenem resistance (%) in Enterobacterales") +
  theme_minimal() +  # Use a minimal theme for clarity
  theme(
    plot.title = element_text(hjust = 0, size = 14, face = "bold"),  # Align title to the left
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),   # Horizontal x-axis labels
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "gray80")   # Light gray grid lines for better contrast
  )

ggsave(filename = "Carbapenem_trendEU.tiff", plot = Carbapenem_trendEU, device = "tiff", 
       path = base_pathOut,
       width = 11, height = 8, dpi = 500, units = "in") 





#######

#Subgroup analyses for enterobacterales first-line antibitiocs ##########
variables_to_keep <- c("Year", first_line_antibiotics)
data_atlas_eu_enterob_firstlin <- data_atlas_eu_enterob[, variables_to_keep]


data_long <- melt(data_atlas_eu_enterob_firstlin, id.vars = "Year", measure.vars = first_line_antibiotics, 
                  variable.name = "Antibiotic", value.name = "Count")
# Summarize the data to get counts of non-missing values per antibiotic per year
data_summary <- dcast(data_long, Antibiotic ~ Year, fun.aggregate = function(x) sum(!is.na(x)))
# View the matrix
data_summary

#Key firstline ATB resistance
data_atlas_eu_enterob_firstlin <- subset(data_atlas_eu_enterob, select = c("Year", "Country","Gender","Age Group","Species", first_line_antibiotics))

# Calculate the number of missing values for each carbapenem
missing_data_summary <- colSums(is.na(data_atlas_eu_enterob_firstlin[, first_line_antibiotics]))
# Calculate the proportion of missing data for each carbapenem
total_rows <- nrow(data_atlas_eu_enterob_firstlin)
missing_data_proportion <- missing_data_summary / total_rows * 100
# Combine the results into a data frame for easy viewing
missing_data_report <- data.frame(
  Antibiotic = first_line_antibiotics,
  Missing_Count = missing_data_summary,
  Missing_Proportion = missing_data_proportion)
# Display table
print(missing_data_report)

#First_lineAtbs
data_atlas_eu_enterob_firstlin[first_line_antibiotics] <- lapply(data_atlas_eu_enterob_firstlin[first_line_antibiotics], as.numeric)
data_atlas_eu_enterob_firstlin$Firstlin_Sum <- rowSums(data_atlas_eu_enterob_firstlin[, first_line_antibiotics], na.rm = TRUE)
data_atlas_eu_enterob_firstlin$Firstlin_res <- ifelse(data_atlas_eu_enterob_firstlin$Firstlin_Sum >= 1, 1, 0)
data_atlas_eu_enterob_firstlin$Firstlin_res  <- as.numeric(data_atlas_eu_enterob_firstlin$Firstlin_res )
# Step 1: Calculate Meropenem resistance by year across the US
Firstlin_resistance_by_year <- aggregate(Firstlin_res  ~ Year, data = data_atlas_eu_enterob_firstlin, FUN = mean, na.rm = TRUE)
Firstlin_resistance_by_year$Firstlin_res <- Firstlin_resistance_by_year$Firstlin_res *100
Firstlin_resistance_by_year$Year <- as.factor(Firstlin_resistance_by_year$Year)

# Create the plot
Firstlin_trendEU<-ggplot(Firstlin_resistance_by_year, aes(x = Year, y = Firstlin_res, group = 1)) +  # Set group = 1 to treat all data as a single group
  geom_line(color = "#D55E00", size = 1, linetype = "solid") +  
  geom_point(color = "black", size = 3, shape = 21, fill = "#D55E60", stroke = 0.5) + 
  labs(title = "",
       x = "Year",
       y = "First-line antibiotic resistance (%) in Enterobacterales") +
  theme_minimal() +  # Use a minimal theme for clarity
  theme(
    plot.title = element_text(hjust = 0, size = 14, face = "bold"),  # Align title to the left
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),   # Horizontal x-axis labels
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "gray80")   # Light gray grid lines for better contrast
  )

ggsave(filename = "Firstlin_trendEU.tiff", plot = Firstlin_trendEU, device = "tiff", 
       path = base_pathOut,
       width = 11, height = 8, dpi = 500, units = "in") 





#######

#Subgroup analyses Enterobacterales for Cephalosporins########
#Key 3rd Generation Cephalosporins : Cefixime_I Ceftriaxone_I Ceftazidime_I Cefpodoxime_I "Ceftazidime avibactam_I" 
variables_to_keep <- c("Year", third_gen_cephalosporins)
data_atlas_eu_enterob_cepha <- data_atlas_eu_enterob[, variables_to_keep]


data_long <- melt(data_atlas_eu_enterob_cepha, id.vars = "Year", measure.vars = third_gen_cephalosporins, 
                  variable.name = "Antibiotic", value.name = "Count")
# Summarize the data to get counts of non-missing values per antibiotic per year
data_summary <- dcast(data_long, Antibiotic ~ Year, fun.aggregate = function(x) sum(!is.na(x)))
# View the matrix
data_summary


data_atlas_eu_enterob_cepha<- subset(data_atlas_eu_enterob, select = c("Year", "Country","Gender","Age Group","Species", third_gen_cephalosporins))
missing_data_summary <- colSums(is.na(data_atlas_eu_enterob_cepha[, third_gen_cephalosporins]))
# Calculate the proportion of missing data for each cephalosporings
total_rows <- nrow(data_atlas_eu_enterob_cepha)
missing_data_proportion <- missing_data_summary / total_rows * 100
# Combine the results into a data frame for easy viewing
missing_data_report <- data.frame(
  Antibiotic = third_gen_cephalosporins,
  Missing_Count = missing_data_summary,
  Missing_Proportion = missing_data_proportion)
# Display the report
print(missing_data_report)

#3rd Generation Cephalosporins ALL ---- #
data_atlas_eu_enterob_cepha[third_gen_cephalosporins] <- lapply(data_atlas_eu_enterob_cepha[third_gen_cephalosporins], as.numeric)
data_atlas_eu_enterob_cepha$Cephalosporin_Sum <- rowSums(data_atlas_eu_enterob_cepha[, third_gen_cephalosporins], na.rm = TRUE)

data_atlas_eu_enterob_cepha$Cephalosporin_res <- ifelse(data_atlas_eu_enterob_cepha$Cephalosporin_Sum >= 1, 1, 0)
data_atlas_eu_enterob_cepha$Cephalosporin_res  <- as.numeric(data_atlas_eu_enterob_cepha$Cephalosporin_res )
# Step 1: Calculate Meropenem resistance by year across the US
cephalosporin_resistance_by_year <- aggregate(Cephalosporin_res  ~ Year, data = data_atlas_eu_enterob_cepha, FUN = mean, na.rm = TRUE)
cephalosporin_resistance_by_year$Cephalosporin_res <- cephalosporin_resistance_by_year$Cephalosporin_res *100
cephalosporin_resistance_by_year$Year <- as.factor(cephalosporin_resistance_by_year$Year)

# Create the plot
Cephalosporin_trendEU<-ggplot(cephalosporin_resistance_by_year, aes(x = Year, y = Cephalosporin_res, group = 1)) +  # Set group = 1 to treat all data as a single group
  geom_line(color = "#D55E00", size = 1, linetype = "solid") +  
  geom_point(color = "black", size = 3, shape = 21, fill = "#D55E60", stroke = 0.5) + 
  labs(title = "",
       x = "Year",
       y = "Third generation cephalosporin resistance (%) in Enterobacterales") +
  theme_minimal() +  # Use a minimal theme for clarity
  theme(
    plot.title = element_text(hjust = 0, size = 14, face = "bold"),  # Align title to the left
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),   # Horizontal x-axis labels
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "gray80")   # Light gray grid lines for better contrast
  )

ggsave(filename = "Cephalosporin_trendEU.tiff", plot = Cephalosporin_trendEU, device = "tiff", 
       path = base_pathOut,
       width = 11, height = 8, dpi = 500, units = "in") 

######

##Subgroup analyses Enterobacterales MDR (at least 3 antibiotic classes) ####

# Convert all antibiotic variables to numeric
antibiotic_list <- unlist(antibiotic_classes)  # Flatten the list of antibiotic classes into a vector of antibiotic names
# Convert each antibiotic variable in the dataframe to numeric
for (antibiotic in antibiotic_list) {
  data_atlas_eu[[antibiotic]] <- as.numeric(data_atlas_eu[[antibiotic]])
}

for (class_name in names(antibiotic_classes)) {
  antibiotic_list <- antibiotic_classes[[class_name]]
  data_atlas_eu[[class_name]] <- rowSums(data_atlas_eu[, antibiotic_list], na.rm = TRUE)
}
for (class_name in names(antibiotic_classes)) {
  data_atlas_eu[[class_name]] <- ifelse(data_atlas_eu[[class_name]] >= 1, 1, 0)
}
# View the updated dataframe with the new variables
head(data_atlas_eu)

#Creating MDR variable:
# Calculate the sum of the antibiotic class variables
data_atlas_eu$Total_Antibiotic_Class_Sum <- rowSums(data_atlas_eu[, c("beta_lactams","carbapenems", "aminoglycosides",  "fluoroquinolones",  "macrolides", "glycopeptides", "lipopeptides", "oxazolidinones", "tetracyclines", "polymyxins", "lincosamides" ,"streptogramins", "nitroimidazoles",  "sulfonamides", "glycylcyclines")], na.rm = TRUE)
# Create the MDR variable based on the sum
data_atlas_eu$MDR <- ifelse(data_atlas_eu$Total_Antibiotic_Class_Sum >= 3, 1, 0)


data_atlas_eu$MDR  <- as.numeric(data_atlas_eu$MDR )

MDR_resistance_by_year <- aggregate(MDR  ~ Year, data = data_atlas_eu, FUN = mean, na.rm = TRUE)
MDR_resistance_by_year$MDR <- MDR_resistance_by_year$MDR *100
MDR_resistance_by_year$Year <- as.factor(MDR_resistance_by_year$Year)

# Create the plot
MDR_trendEU<-ggplot(MDR_resistance_by_year, aes(x = Year, y = MDR, group = 1)) +  # Set group = 1 to treat all data as a single group
  geom_line(color = "#D55E00", size = 1, linetype = "solid") +  
  geom_point(color = "black", size = 3, shape = 21, fill = "#D55E60", stroke = 0.5) + 
  labs(title = "",
       x = "Year",
       y = "Multidrug resistance (%) in Enterobacterales") +
  theme_minimal() +  # Use a minimal theme for clarity
  theme(
    plot.title = element_text(hjust = 0, size = 14, face = "bold"),  # Align title to the left
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),   # Horizontal x-axis labels
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "gray80")   # Light gray grid lines for better contrast
  )

ggsave(filename = "MDR_trendEU.tiff", plot = MDR_trendEU, device = "tiff", 
       path = base_pathOut,
       width = 11, height = 8, dpi = 500, units = "in") 


#####

#4Maps with first-line, carbapenems, cephalosporin, MDR resistance per state#####
MDR_resistance_by_year_state <- aggregate(MDR ~ Country, data = data_atlas_eu, FUN = mean, na.rm = TRUE)
MDR_resistance_by_year_state$MDR <- MDR_resistance_by_year_state$MDR *100

cephalosporin_resistance_by_year_state <- aggregate(Cephalosporin_res  ~  Country, data = data_atlas_eu_enterob_cepha, FUN = mean, na.rm = TRUE)
cephalosporin_resistance_by_year_state$Cephalosporin_res <- cephalosporin_resistance_by_year_state$Cephalosporin_res *100

carbapenem_resistance_by_year_state <- aggregate(Carbapenem_res  ~ Country, data = data_atlas_eu_enterob_carbap, FUN = mean, na.rm = TRUE)
carbapenem_resistance_by_year_state$Carbapenem_res <- carbapenem_resistance_by_year_state$Carbapenem_res *100
#carbapenem_resistance_by_year_state$Year <- as.factor(carbapenem_resistance_by_year_state$Year)

Firstlin_resistance_by_year_state <- aggregate(Firstlin_res  ~ Country, data = data_atlas_eu_enterob_firstlin, FUN = mean, na.rm = TRUE)
Firstlin_resistance_by_year_state$Firstlin_res <- Firstlin_resistance_by_year_state$Firstlin_res *100



# Load EU map data
# Download the Europe shapefile from Natural Earth
europe_shapefile <- ne_countries(continent = "Europe", returnclass = "sf")
# Ensure all geometries are valid
europe_shapefile <- st_make_valid(europe_shapefile)
# Define the bounding box to include the full extent of Italy and mainland Europe
bounding_box <- st_bbox(c(xmin = -10, xmax = 50, ymin = 30, ymax = 70), crs = st_crs(europe_shapefile))
# Crop the shapefile to the defined bounding box
europe_shapefile <- st_crop(europe_shapefile, bounding_box)
plot(europe_shapefile)
# Filter to remove specific non-European territories if needed
europe_shapefile <- europe_shapefile %>%
  filter(!sovereignt %in% c("Canary Islands", "Madeira", "Azores", "Cyprus"))
europe_shapefile <- europe_shapefile %>%
  mutate(sovereignt = ifelse(sovereignt == "Slovakia", "Slovak Republic", sovereignt))
europe_shapefile <- europe_shapefile %>%
  mutate(sovereignt = ifelse(sovereignt == "Czechia", "Czech Republic", sovereignt))
europe_shapefile <- europe_shapefile %>%
  mutate(sovereignt = ifelse(sovereignt == "Republic of Serbia", "Serbia", sovereignt))

# Prepare the data: ensure state names match and join with map data
carbapenem_resistance_by_year_state$sovereignt <- carbapenem_resistance_by_year_state$Country
cephalosporin_resistance_by_year_state$sovereignt <- cephalosporin_resistance_by_year_state$Country
MDR_resistance_by_year_state$sovereignt <- MDR_resistance_by_year_state$Country
Firstlin_resistance_by_year_state$sovereignt <- Firstlin_resistance_by_year_state$Country

eu_map_carbapenem <- europe_shapefile %>%
  left_join(carbapenem_resistance_by_year_state, by = "sovereignt")
eu_map_cephalosporin <- europe_shapefile  %>%
  left_join(cephalosporin_resistance_by_year_state, by = "sovereignt")
eu_map_MDR <- europe_shapefile  %>%
  left_join(MDR_resistance_by_year_state, by = "sovereignt")
eu_map_Firstlin <- europe_shapefile  %>%
  left_join(Firstlin_resistance_by_year_state, by = "sovereignt")


lancet_theme <- theme(
  text = element_text(family = "serif", color = "black"),
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  panel.background = element_rect(fill = "white", color = "white"),
  legend.position = "bottom",
  legend.title = element_text(size = 10, face = "bold"),
  legend.text = element_text(size = 8)
)

#Individual graphs:
#xlims <- range(st_coordinates(europe_shapefile)[,1])
#ylims <- range(st_coordinates(europe_shapefile)[,2])

# Plot for Carbapenem Resistance
carbapenem_heatmap2 <- ggplot(eu_map_carbapenem, aes(fill = Carbapenem_res)) +
  geom_sf(color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey80") +
  labs(title = "F. Carbapenem resistance in Europe", fill = "Carbapenem\nresistance (%)") +
  coord_sf(expand = FALSE) +  # No expansion beyond the limits of the data
  theme_minimal() +
  lancet_theme +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom",
    axis.text = element_blank(),     # Remove axis text
    axis.title = element_blank(),    # Remove axis titles
    axis.ticks = element_blank(),    # Remove axis ticks
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white", colour = NA), # Optional: set panel background to white
    plot.margin = unit(c(0, 0, 0, 0), "cm") # Set all margins to zero
  )


# Plot for Cephalosporin Resistance
cephalosporin_heatmap2 <- ggplot(eu_map_cephalosporin, aes(fill = Cephalosporin_res)) +
  geom_sf(color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey80") +
  labs(title = "G. Cephalosporin resistance in Europe", fill = "Cephalosporin\n resistance (%)") +
  coord_sf(expand = FALSE) + 
  theme_minimal() +
  lancet_theme +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom",
    axis.text = element_blank(),     # Remove axis text
    axis.title = element_blank(),    # Remove axis titles
    axis.ticks = element_blank(),    # Remove axis ticks
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white", colour = NA), # Optional: set panel background to white
    plot.margin = unit(c(1, 1, 1, 1), "mm") # Set all margins to 1mm
  )

# Plot for MDR Resistance
MDR_heatmap2 <- ggplot(eu_map_MDR, aes(fill = MDR)) +
  geom_sf(color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey80") +
  labs(title = "H. MDR in Europe", fill = "MDR (%)") +
  coord_sf(expand = FALSE) + 
  theme_minimal() +
  lancet_theme +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom",
    axis.text = element_blank(),     # Remove axis text
    axis.title = element_blank(),    # Remove axis titles
    axis.ticks = element_blank(),    # Remove axis ticks
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white", colour = NA), # Optional: set panel background to white
    plot.margin = unit(c(1, 1, 1, 1), "mm") # Set all margins to 1mm
  )


# Firstline Resistance
Firstlin_heatmap2 <- ggplot(eu_map_Firstlin, aes(fill = Firstlin_res)) +
  geom_sf(color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey80") +
  labs(title = "E. First-line AMR in Europe", fill = "First-line antibiotic\n resistance (%)") +
  coord_sf(expand = FALSE) + 
  theme_minimal() +
  lancet_theme +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom",
    axis.text = element_blank(),     # Remove axis text
    axis.title = element_blank(),    # Remove axis titles
    axis.ticks = element_blank(),    # Remove axis ticks
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white", colour = NA), # Optional: set panel background to white
    plot.margin = unit(c(1, 1, 1, 1), "mm") # Set all margins to 1mm
  )

######

#--EU + US states ---------#
#Figure 1: COMBINED GRAPHS or HEATmaps: ######
carbapenem_heatmap2 <- carbapenem_heatmap2 +
  theme(plot.title = element_text(hjust = 0, size = 12, face = "bold"))
cephalosporin_heatmap2 <- cephalosporin_heatmap2 +
  theme(plot.title = element_text(hjust = 0, size = 12, face = "bold"))
MDR_heatmap2 <- MDR_heatmap2 +
  theme(plot.title = element_text(hjust = 0, size = 12, face = "bold"))
Firstlin_heatmap2 <- Firstlin_heatmap2 +
  theme(plot.title = element_text(hjust = 0, size = 12, face = "bold"))

library(gridExtra)
# Combine the three heatmaps
#combined_plot_fig1 <- grid.arrange(Firstlin_heatmap2,Firstlin_heatmap,  carbapenem_heatmap, carbapenem_heatmap2, cephalosporin_heatmap, cephalosporin_heatmap2, MDR_heatmap, MDR_heatmap2, ncol = 2)

# Combine the plots using patchwork
combined_plot_fig1 <- (Firstlin_heatmap + carbapenem_heatmap) / 
  (cephalosporin_heatmap + MDR_heatmap ) / 
  (Firstlin_heatmap2 + carbapenem_heatmap2) / 
  (cephalosporin_heatmap2 + MDR_heatmap2)

#Configure orientation
combined_plot_fig1 <- combined_plot_fig1 + 
  plot_layout(ncol = 2, nrow = 2, heights = rep(1, 4))

# Print the combined plot to check the layout
print(combined_plot_fig1)
# Save the combined plot
ggsave(filename = "resistance_overtime_combined.tiff", plot = combined_plot_fig1, device = "tiff", path = base_pathOut,
       width = 13, height = 8, dpi = 500, units = "in")

######

#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#

#Variable to use population per European country 2019:
pop_europe <- europe_shapefile %>%
  dplyr::select(pop_est, sovereignt)

#Variable to use population per US state 2019:


#------------------------------------------------------------------------------------------------------#
#Data for later models, function to generate dataframes:
#------------------------------------------------------------------------------------------------------#
######
compute_resistance_summary <- function(data, antibiotics, level = "country") {
  # Select relevant columns
  variables_to_keep <- c("Year", "State", antibiotics)
  data_selected <- data[, variables_to_keep]
  
  # Melt data to long format
  data_long <- reshape2::melt(data_selected, id.vars = c("Year", "State"), 
                              measure.vars = antibiotics, 
                              variable.name = "Antibiotic", 
                              value.name = "Count")
  
  # Summarize data to get counts of non-missing values per antibiotic per year and state
  data_summary <- reshape2::dcast(data_long, Antibiotic ~ Year + State, 
                                  fun.aggregate = function(x) sum(!is.na(x)))
  # Display summary table
  print(data_summary)
  
  # Select additional columns for more detailed analysis
  data_detailed <- subset(data, select = c("Year", "State", "Gender", "Age Group", "Species", antibiotics))
  
  # Calculate missing data
  missing_data_summary <- colSums(is.na(data_detailed[, antibiotics]))
  total_rows <- nrow(data_detailed)
  missing_data_proportion <- missing_data_summary / total_rows * 100
  missing_data_report <- data.frame(Antibiotic = antibiotics, 
                                    Missing_Count = missing_data_summary, 
                                    Missing_Proportion = missing_data_proportion)
  print(missing_data_report)
  
  # Convert antibiotics to numeric and calculate resistance
  data_detailed[antibiotics] <- lapply(data_detailed[antibiotics], as.numeric)
  data_detailed$Antibiotic_Sum <- rowSums(data_detailed[, antibiotics], na.rm = TRUE)
  data_detailed$Resistance <- ifelse(data_detailed$Antibiotic_Sum >= 1, 1, 0)
  
  # Compute counts of AMR positive isolates and all isolates tested per year and state
  if (level == "state") {
    # Calculate resistance percentage, AMR positives, and total isolates
    resistance_summary <- aggregate(Resistance ~ Year + State, data = data_detailed, FUN = mean, na.rm = TRUE)
    resistance_summary$AMR_Positive <- aggregate(Resistance ~ Year + State, data = data_detailed, FUN = sum, na.rm = TRUE)$Resistance
    resistance_summary$Total_Isolates <- aggregate(Resistance ~ Year + State, data = data_detailed, FUN = length)$Resistance
  } else {
    resistance_summary <- aggregate(Resistance ~ Year, data = data_detailed, FUN = mean, na.rm = TRUE)
    resistance_summary$AMR_Positive <- aggregate(Resistance ~ Year, data = data_detailed, FUN = sum, na.rm = TRUE)$Resistance
    resistance_summary$Total_Isolates <- aggregate(Resistance ~ Year, data = data_detailed, FUN = length)$Resistance
  }
  
  # Convert year to factor
  resistance_summary$Year <- as.factor(resistance_summary$Year)
  
  # Return the resistance summary
  return(resistance_summary)
}
compute_mdr_summary <- function(data, antibiotic_classes, level = "country") {
  # Flatten the list of antibiotic classes into a vector of antibiotic names
  antibiotic_list <- unlist(antibiotic_classes)
  
  # Convert each antibiotic variable in the dataframe to numeric
  data[antibiotic_list] <- lapply(data[antibiotic_list], as.numeric)
  
  # Calculate resistance for each antibiotic class
  for (class_name in names(antibiotic_classes)) {
    class_antibiotics <- antibiotic_classes[[class_name]]
    data[[class_name]] <- rowSums(data[, class_antibiotics], na.rm = TRUE)
    data[[class_name]] <- ifelse(data[[class_name]] >= 1, 1, 0)
  }
  
  # Calculate the sum of the antibiotic class variables for MDR determination
  data$Total_Antibiotic_Class_Sum <- rowSums(data[, names(antibiotic_classes)], na.rm = TRUE)
  data$MDR <- ifelse(data$Total_Antibiotic_Class_Sum >= 3, 1, 0)
  
  # Aggregate MDR data by year or state level
  if (level == "state") {
    # Calculate MDR percentage, MDR positives, and total isolates
    mdr_summary <- aggregate(MDR ~ Year + State, data = data, FUN = mean, na.rm = TRUE)
    mdr_summary$MDR_Positive <- aggregate(MDR ~ Year + State, data = data, FUN = sum, na.rm = TRUE)$MDR
    mdr_summary$Total_Isolates <- aggregate(MDR ~ Year + State, data = data, FUN = length)$MDR
  } else {
    mdr_summary <- aggregate(MDR ~ Year, data = data, FUN = mean, na.rm = TRUE)
    mdr_summary$MDR_Positive <- aggregate(MDR ~ Year, data = data, FUN = sum, na.rm = TRUE)$MDR
    mdr_summary$Total_Isolates <- aggregate(MDR ~ Year, data = data, FUN = length)$MDR
  }
  
  # Convert year to factor
  mdr_summary$Year <- as.factor(mdr_summary$Year)
  
  # Return the MDR summary
  return(mdr_summary)
}

#Building the datasets for analyses:
state_resistance_carbap<-compute_resistance_summary(data_atlas_us_enterob, carbapenems, "state")
state_resistance_cephalos<-compute_resistance_summary(data_atlas_us_enterob, third_gen_cephalosporins, "state")
state_resistance_firstline<-compute_resistance_summary(data_atlas_us_enterob, first_line_antibiotics , "state")
state_resistance_mdr<- compute_mdr_summary(data_atlas_us_enterob, antibiotic_classes, level = "state")

state_resistance_carbap_o<-compute_resistance_summary(data_atlas_us_enterob_orig, carbapenems, "state")
state_resistance_cephalos_o<-compute_resistance_summary(data_atlas_us_enterob_orig, third_gen_cephalosporins, "state")
state_resistance_firstline_o<-compute_resistance_summary(data_atlas_us_enterob_orig, first_line_antibiotics , "state")
state_resistance_mdr_o<- compute_mdr_summary(data_atlas_us_enterob_orig, antibiotic_classes, level = "state")

compute_resistance_summary_eu <- function(data, antibiotics, level = "country") {
  # Select relevant columns
  variables_to_keep <- c("Year", "Country", antibiotics)
  data_selected <- data[, variables_to_keep]
  
  # Melt data to long format
  data_long <- reshape2::melt(data_selected, id.vars = c("Year", "Country"), 
                              measure.vars = antibiotics, 
                              variable.name = "Antibiotic", 
                              value.name = "Count")
  
  # Summarize data to get counts of non-missing values per antibiotic per year and state
  data_summary <- reshape2::dcast(data_long, Antibiotic ~ Year + Country, 
                                  fun.aggregate = function(x) sum(!is.na(x)))
  # Display summary table
  print(data_summary)
  
  # Select additional columns for more detailed analysis
  data_detailed <- subset(data, select = c("Year", "Country", "Gender", "Age Group", "Species", antibiotics))
  
  # Calculate missing data
  missing_data_summary <- colSums(is.na(data_detailed[, antibiotics]))
  total_rows <- nrow(data_detailed)
  missing_data_proportion <- missing_data_summary / total_rows * 100
  missing_data_report <- data.frame(Antibiotic = antibiotics, 
                                    Missing_Count = missing_data_summary, 
                                    Missing_Proportion = missing_data_proportion)
  print(missing_data_report)
  
  # Convert antibiotics to numeric and calculate resistance
  data_detailed[antibiotics] <- lapply(data_detailed[antibiotics], as.numeric)
  data_detailed$Antibiotic_Sum <- rowSums(data_detailed[, antibiotics], na.rm = TRUE)
  data_detailed$Resistance <- ifelse(data_detailed$Antibiotic_Sum >= 1, 1, 0)
  
  # Compute counts of AMR positive isolates and all isolates tested per year and country
  if (level == "country") {
    # Calculate resistance percentage, AMR positives, and total isolates
    resistance_summary <- aggregate(Resistance ~ Year + Country, data = data_detailed, FUN = mean, na.rm = TRUE)
    resistance_summary$AMR_Positive <- aggregate(Resistance ~ Year + Country, data = data_detailed, FUN = sum, na.rm = TRUE)$Resistance
    resistance_summary$Total_Isolates <- aggregate(Resistance ~ Year + Country, data = data_detailed, FUN = length)$Resistance
  } else {
    resistance_summary <- aggregate(Resistance ~ Year, data = data_detailed, FUN = mean, na.rm = TRUE)
    resistance_summary$AMR_Positive <- aggregate(Resistance ~ Year, data = data_detailed, FUN = sum, na.rm = TRUE)$Resistance
    resistance_summary$Total_Isolates <- aggregate(Resistance ~ Year, data = data_detailed, FUN = length)$Resistance
  }
  
  # Convert year to factor
  resistance_summary$Year <- as.factor(resistance_summary$Year)
  
  # Return the resistance summary
  return(resistance_summary)
}
compute_mdr_summary_eu <- function(data, antibiotic_classes, level = "country") {
  # Flatten the list of antibiotic classes into a vector of antibiotic names
  antibiotic_list <- unlist(antibiotic_classes)
  
  # Convert each antibiotic variable in the dataframe to numeric
  data[antibiotic_list] <- lapply(data[antibiotic_list], as.numeric)
  
  # Calculate resistance for each antibiotic class
  for (class_name in names(antibiotic_classes)) {
    class_antibiotics <- antibiotic_classes[[class_name]]
    data[[class_name]] <- rowSums(data[, class_antibiotics], na.rm = TRUE)
    data[[class_name]] <- ifelse(data[[class_name]] >= 1, 1, 0)
  }
  
  # Calculate the sum of the antibiotic class variables for MDR determination
  data$Total_Antibiotic_Class_Sum <- rowSums(data[, names(antibiotic_classes)], na.rm = TRUE)
  data$MDR <- ifelse(data$Total_Antibiotic_Class_Sum >= 3, 1, 0)
  
  # Aggregate MDR data by year or state level
  if (level == "country") {
    # Calculate MDR percentage, MDR positives, and total isolates
    mdr_summary <- aggregate(MDR ~ Year + Country, data = data, FUN = mean, na.rm = TRUE)
    mdr_summary$MDR_Positive <- aggregate(MDR ~ Year + Country, data = data, FUN = sum, na.rm = TRUE)$MDR
    mdr_summary$Total_Isolates <- aggregate(MDR ~ Year + Country, data = data, FUN = length)$MDR
  } else {
    mdr_summary <- aggregate(MDR ~ Year, data = data, FUN = mean, na.rm = TRUE)
    mdr_summary$MDR_Positive <- aggregate(MDR ~ Year, data = data, FUN = sum, na.rm = TRUE)$MDR
    mdr_summary$Total_Isolates <- aggregate(MDR ~ Year, data = data, FUN = length)$MDR
  }
  
  # Convert year to factor
  mdr_summary$Year <- as.factor(mdr_summary$Year)
  
  # Return the MDR summary
  return(mdr_summary)
}

country_resistance_carbap_eu<-compute_resistance_summary_eu(data_atlas_eu_enterob, carbapenems, "country")
country_resistance_cephalos_eu<-compute_resistance_summary_eu(data_atlas_eu_enterob, third_gen_cephalosporins, "country")
country_resistance_firstline_eu<-compute_resistance_summary_eu(data_atlas_eu_enterob, first_line_antibiotics , "country")
country_resistance_mdr_eu<- compute_mdr_summary_eu(data_atlas_eu_enterob, antibiotic_classes, level = "country")

country_resistance_carbap_o_eu<-compute_resistance_summary_eu(data_atlas_eu_enterob_orig, carbapenems, "country")
country_resistance_cephalos_o_eu<-compute_resistance_summary_eu(data_atlas_eu_enterob_orig, third_gen_cephalosporins, "country")
country_resistance_firstline_o_eu<-compute_resistance_summary_eu(data_atlas_eu_enterob_orig, first_line_antibiotics , "country")
country_resistance_mdr_o_eu<- compute_mdr_summary_eu(data_atlas_eu_enterob_orig, antibiotic_classes, level = "country")


#-------------------------------#
#Gender:
compute_resistance_summaryGen <- function(data, antibiotics, level = "country") {
  # Select relevant columns
  variables_to_keep <- c("Year", "State", "Gender", antibiotics)
  data_selected <- data[, variables_to_keep]
  
  # Melt data to long format
  data_long <- reshape2::melt(data_selected, id.vars = c("Year", "State", "Gender"), 
                              measure.vars = antibiotics, 
                              variable.name = "Antibiotic", 
                              value.name = "Count")
  
  # Summarize data to get counts of non-missing values per antibiotic per year, state, and gender
  data_summary <- reshape2::dcast(data_long, Antibiotic ~ Year + State + Gender, 
                                  fun.aggregate = function(x) sum(!is.na(x)))
  # Display summary table
  print(data_summary)
  
  # Select additional columns for more detailed analysis
  data_detailed <- subset(data, select = c("Year", "State", "Gender", "Age Group", "Species", antibiotics))
  
  # Calculate missing data
  missing_data_summary <- colSums(is.na(data_detailed[, antibiotics]))
  total_rows <- nrow(data_detailed)
  missing_data_proportion <- missing_data_summary / total_rows * 100
  missing_data_report <- data.frame(Antibiotic = antibiotics, 
                                    Missing_Count = missing_data_summary, 
                                    Missing_Proportion = missing_data_proportion)
  print(missing_data_report)
  
  # Convert antibiotics to numeric and calculate resistance
  data_detailed[antibiotics] <- lapply(data_detailed[antibiotics], as.numeric)
  data_detailed$Antibiotic_Sum <- rowSums(data_detailed[, antibiotics], na.rm = TRUE)
  data_detailed$Resistance <- ifelse(data_detailed$Antibiotic_Sum >= 1, 1, 0)
  
  # Compute counts of AMR positive isolates and all isolates tested per year, state, and gender
  if (level == "state") {
    # Calculate resistance percentage, AMR positives, and total isolates
    resistance_summary <- aggregate(Resistance ~ Year + State + Gender, data = data_detailed, FUN = mean, na.rm = TRUE)
    resistance_summary$AMR_Positive <- aggregate(Resistance ~ Year + State + Gender, data = data_detailed, FUN = sum, na.rm = TRUE)$Resistance
    resistance_summary$Total_Isolates <- aggregate(Resistance ~ Year + State + Gender, data = data_detailed, FUN = length)$Resistance
  } else {
    resistance_summary <- aggregate(Resistance ~ Year + Gender, data = data_detailed, FUN = mean, na.rm = TRUE)
    resistance_summary$AMR_Positive <- aggregate(Resistance ~ Year + Gender, data = data_detailed, FUN = sum, na.rm = TRUE)$Resistance
    resistance_summary$Total_Isolates <- aggregate(Resistance ~ Year + Gender, data = data_detailed, FUN = length)$Resistance
  }
  
  # Convert year to factor
  resistance_summary$Year <- as.factor(resistance_summary$Year)
  
  # Return the resistance summary
  return(resistance_summary)
}
compute_mdr_summaryGen <- function(data, antibiotic_classes, level = "country") {
  # Flatten the list of antibiotic classes into a vector of antibiotic names
  antibiotic_list <- unlist(antibiotic_classes)
  
  # Convert each antibiotic variable in the dataframe to numeric
  data[antibiotic_list] <- lapply(data[antibiotic_list], as.numeric)
  
  # Calculate resistance for each antibiotic class
  for (class_name in names(antibiotic_classes)) {
    class_antibiotics <- antibiotic_classes[[class_name]]
    data[[class_name]] <- rowSums(data[, class_antibiotics], na.rm = TRUE)
    data[[class_name]] <- ifelse(data[[class_name]] >= 1, 1, 0)
  }
  
  # Calculate the sum of the antibiotic class variables for MDR determination
  data$Total_Antibiotic_Class_Sum <- rowSums(data[, names(antibiotic_classes)], na.rm = TRUE)
  data$MDR <- ifelse(data$Total_Antibiotic_Class_Sum >= 3, 1, 0)
  
  # Aggregate MDR data by year, state, and gender level
  if (level == "state") {
    # Calculate MDR percentage, MDR positives, and total isolates
    mdr_summary <- aggregate(MDR ~ Year + State + Gender, data = data, FUN = mean, na.rm = TRUE)
    mdr_summary$MDR_Positive <- aggregate(MDR ~ Year + State + Gender, data = data, FUN = sum, na.rm = TRUE)$MDR
    mdr_summary$Total_Isolates <- aggregate(MDR ~ Year + State + Gender, data = data, FUN = length)$MDR
  } else {
    mdr_summary <- aggregate(MDR ~ Year + Gender, data = data, FUN = mean, na.rm = TRUE)
    mdr_summary$MDR_Positive <- aggregate(MDR ~ Year + Gender, data = data, FUN = sum, na.rm = TRUE)$MDR
    mdr_summary$Total_Isolates <- aggregate(MDR ~ Year + Gender, data = data, FUN = length)$MDR
  }
  
  # Convert year to factor
  mdr_summary$Year <- as.factor(mdr_summary$Year)
  
  # Return the MDR summary
  return(mdr_summary)
}

state_resistance_carbapGen<-compute_resistance_summaryGen(data_atlas_us_enterob, carbapenems, "state")
state_resistance_cephalosGen<-compute_resistance_summaryGen(data_atlas_us_enterob, third_gen_cephalosporins, "state")
state_resistance_firstlineGen<-compute_resistance_summaryGen(data_atlas_us_enterob, first_line_antibiotics , "state")
state_resistance_mdrGen<- compute_mdr_summaryGen(data_atlas_us_enterob, antibiotic_classes, level = "state")

state_resistance_carbap_oGen<-compute_resistance_summaryGen(data_atlas_us_enterob_orig, carbapenems, "state")
state_resistance_cephalos_oGen<-compute_resistance_summaryGen(data_atlas_us_enterob_orig, third_gen_cephalosporins, "state")
state_resistance_firstline_oGen<-compute_resistance_summaryGen(data_atlas_us_enterob_orig, first_line_antibiotics , "state")
state_resistance_mdr_oGen<- compute_mdr_summaryGen(data_atlas_us_enterob_orig, antibiotic_classes, level = "state")

compute_resistance_summary_euGen <- function(data, antibiotics, level = "country") {
  # Select relevant columns
  variables_to_keep <- c("Year", "Country", "Gender", antibiotics)
  data_selected <- data[, variables_to_keep]
  
  # Melt data to long format
  data_long <- reshape2::melt(data_selected, id.vars = c("Year", "Country", "Gender"), 
                              measure.vars = antibiotics, 
                              variable.name = "Antibiotic", 
                              value.name = "Count")
  
  # Summarize data to get counts of non-missing values per antibiotic per year, country, and gender
  data_summary <- reshape2::dcast(data_long, Antibiotic ~ Year + Country + Gender, 
                                  fun.aggregate = function(x) sum(!is.na(x)))
  # Display summary table
  print(data_summary)
  
  # Select additional columns for more detailed analysis
  data_detailed <- subset(data, select = c("Year", "Country", "Gender", "Age Group", "Species", antibiotics))
  
  # Calculate missing data
  missing_data_summary <- colSums(is.na(data_detailed[, antibiotics]))
  total_rows <- nrow(data_detailed)
  missing_data_proportion <- missing_data_summary / total_rows * 100
  missing_data_report <- data.frame(Antibiotic = antibiotics, 
                                    Missing_Count = missing_data_summary, 
                                    Missing_Proportion = missing_data_proportion)
  print(missing_data_report)
  
  # Convert antibiotics to numeric and calculate resistance
  data_detailed[antibiotics] <- lapply(data_detailed[antibiotics], as.numeric)
  data_detailed$Antibiotic_Sum <- rowSums(data_detailed[, antibiotics], na.rm = TRUE)
  data_detailed$Resistance <- ifelse(data_detailed$Antibiotic_Sum >= 1, 1, 0)
  
  # Compute counts of AMR positive isolates and all isolates tested per year, country, and gender
  if (level == "country") {
    # Calculate resistance percentage, AMR positives, and total isolates
    resistance_summary <- aggregate(Resistance ~ Year + Country + Gender, data = data_detailed, FUN = mean, na.rm = TRUE)
    resistance_summary$AMR_Positive <- aggregate(Resistance ~ Year + Country + Gender, data = data_detailed, FUN = sum, na.rm = TRUE)$Resistance
    resistance_summary$Total_Isolates <- aggregate(Resistance ~ Year + Country + Gender, data = data_detailed, FUN = length)$Resistance
  } else {
    resistance_summary <- aggregate(Resistance ~ Year + Gender, data = data_detailed, FUN = mean, na.rm = TRUE)
    resistance_summary$AMR_Positive <- aggregate(Resistance ~ Year + Gender, data = data_detailed, FUN = sum, na.rm = TRUE)$Resistance
    resistance_summary$Total_Isolates <- aggregate(Resistance ~ Year + Gender, data = data_detailed, FUN = length)$Resistance
  }
  
  # Convert year to factor
  resistance_summary$Year <- as.factor(resistance_summary$Year)
  
  # Return the resistance summary
  return(resistance_summary)
}
compute_mdr_summary_euGen <- function(data, antibiotic_classes, level = "country") {
  # Flatten the list of antibiotic classes into a vector of antibiotic names
  antibiotic_list <- unlist(antibiotic_classes)
  
  # Convert each antibiotic variable in the dataframe to numeric
  data[antibiotic_list] <- lapply(data[antibiotic_list], as.numeric)
  
  # Calculate resistance for each antibiotic class
  for (class_name in names(antibiotic_classes)) {
    class_antibiotics <- antibiotic_classes[[class_name]]
    data[[class_name]] <- rowSums(data[, class_antibiotics], na.rm = TRUE)
    data[[class_name]] <- ifelse(data[[class_name]] >= 1, 1, 0)
  }
  
  # Calculate the sum of the antibiotic class variables for MDR determination
  data$Total_Antibiotic_Class_Sum <- rowSums(data[, names(antibiotic_classes)], na.rm = TRUE)
  data$MDR <- ifelse(data$Total_Antibiotic_Class_Sum >= 3, 1, 0)
  
  # Aggregate MDR data by year or state level
  if (level == "country") {
    # Calculate MDR percentage, MDR positives, and total isolates
    mdr_summary <- aggregate(MDR ~ Year + Country + Gender, data = data, FUN = mean, na.rm = TRUE)
    mdr_summary$MDR_Positive <- aggregate(MDR ~ Year + Country + Gender, data = data, FUN = sum, na.rm = TRUE)$MDR
    mdr_summary$Total_Isolates <- aggregate(MDR ~ Year + Country + Gender, data = data, FUN = length)$MDR
  } else {
    mdr_summary <- aggregate(MDR ~ Year + Gender, data = data, FUN = mean, na.rm = TRUE)
    mdr_summary$MDR_Positive <- aggregate(MDR ~ Year + Gender, data = data, FUN = sum, na.rm = TRUE)$MDR
    mdr_summary$Total_Isolates <- aggregate(MDR ~ Year + Gender, data = data, FUN = length)$MDR
  }
  
  # Convert year to factor
  mdr_summary$Year <- as.factor(mdr_summary$Year)
  
  # Return the MDR summary
  return(mdr_summary)
}

country_resistance_carbap_euGen<-compute_resistance_summary_euGen(data_atlas_eu_enterob, carbapenems, "country")
country_resistance_cephalos_euGen<-compute_resistance_summary_euGen(data_atlas_eu_enterob, third_gen_cephalosporins, "country")
country_resistance_firstline_euGen<-compute_resistance_summary_euGen(data_atlas_eu_enterob, first_line_antibiotics , "country")
country_resistance_mdr_euGen<- compute_mdr_summary_euGen(data_atlas_eu_enterob, antibiotic_classes, level = "country")

country_resistance_carbap_o_euGen<-compute_resistance_summary_euGen(data_atlas_eu_enterob_orig, carbapenems, "country")
country_resistance_cephalos_o_euGen<-compute_resistance_summary_euGen(data_atlas_eu_enterob_orig, third_gen_cephalosporins, "country")
country_resistance_firstline_o_euGen<-compute_resistance_summary_euGen(data_atlas_eu_enterob_orig, first_line_antibiotics , "country")
country_resistance_mdr_o_euGen<- compute_mdr_summary_euGen(data_atlas_eu_enterob_orig, antibiotic_classes, level = "country")



#-------------------------------#
#Age groups:
compute_resistance_summaryAgeg <- function(data, antibiotics, level = "country") {
  # Select relevant columns
  variables_to_keep <- c("Year", "State", "Agegroup", antibiotics)
  data_selected <- data[, variables_to_keep]
  
  # Melt data to long format
  data_long <- reshape2::melt(data_selected, id.vars = c("Year", "State", "Agegroup"), 
                              measure.vars = antibiotics, 
                              variable.name = "Antibiotic", 
                              value.name = "Count")
  
  # Summarize data to get counts of non-missing values per antibiotic per year, state, and gender
  data_summary <- reshape2::dcast(data_long, Antibiotic ~ Year + State + Agegroup, 
                                  fun.aggregate = function(x) sum(!is.na(x)))
  # Display summary table
  print(data_summary)
  
  # Select additional columns for more detailed analysis
  data_detailed <- subset(data, select = c("Year", "State", "Agegroup", "Species", antibiotics))
  
  # Calculate missing data
  missing_data_summary <- colSums(is.na(data_detailed[, antibiotics]))
  total_rows <- nrow(data_detailed)
  missing_data_proportion <- missing_data_summary / total_rows * 100
  missing_data_report <- data.frame(Antibiotic = antibiotics, 
                                    Missing_Count = missing_data_summary, 
                                    Missing_Proportion = missing_data_proportion)
  print(missing_data_report)
  
  # Convert antibiotics to numeric and calculate resistance
  data_detailed[antibiotics] <- lapply(data_detailed[antibiotics], as.numeric)
  data_detailed$Antibiotic_Sum <- rowSums(data_detailed[, antibiotics], na.rm = TRUE)
  data_detailed$Resistance <- ifelse(data_detailed$Antibiotic_Sum >= 1, 1, 0)
  
  # Compute counts of AMR positive isolates and all isolates tested per year, state, and gender
  if (level == "state") {
    # Calculate resistance percentage, AMR positives, and total isolates
    resistance_summary <- aggregate(Resistance ~ Year + State + Agegroup, data = data_detailed, FUN = mean, na.rm = TRUE)
    resistance_summary$AMR_Positive <- aggregate(Resistance ~ Year + State + Agegroup, data = data_detailed, FUN = sum, na.rm = TRUE)$Resistance
    resistance_summary$Total_Isolates <- aggregate(Resistance ~ Year + State + Agegroup, data = data_detailed, FUN = length)$Resistance
  } else {
    resistance_summary <- aggregate(Resistance ~ Year + Agegroup, data = data_detailed, FUN = mean, na.rm = TRUE)
    resistance_summary$AMR_Positive <- aggregate(Resistance ~ Year + Agegroup, data = data_detailed, FUN = sum, na.rm = TRUE)$Resistance
    resistance_summary$Total_Isolates <- aggregate(Resistance ~ Year + Agegroup, data = data_detailed, FUN = length)$Resistance
  }
  
  # Convert year to factor
  resistance_summary$Year <- as.factor(resistance_summary$Year)
  
  # Return the resistance summary
  return(resistance_summary)
}
compute_mdr_summaryAgeg <- function(data, antibiotic_classes, level = "country") {
  # Flatten the list of antibiotic classes into a vector of antibiotic names
  antibiotic_list <- unlist(antibiotic_classes)
  
  # Convert each antibiotic variable in the dataframe to numeric
  data[antibiotic_list] <- lapply(data[antibiotic_list], as.numeric)
  
  # Calculate resistance for each antibiotic class
  for (class_name in names(antibiotic_classes)) {
    class_antibiotics <- antibiotic_classes[[class_name]]
    data[[class_name]] <- rowSums(data[, class_antibiotics], na.rm = TRUE)
    data[[class_name]] <- ifelse(data[[class_name]] >= 1, 1, 0)
  }
  
  # Calculate the sum of the antibiotic class variables for MDR determination
  data$Total_Antibiotic_Class_Sum <- rowSums(data[, names(antibiotic_classes)], na.rm = TRUE)
  data$MDR <- ifelse(data$Total_Antibiotic_Class_Sum >= 3, 1, 0)
  
  # Aggregate MDR data by year, state, and gender level
  if (level == "state") {
    # Calculate MDR percentage, MDR positives, and total isolates
    mdr_summary <- aggregate(MDR ~ Year + State + Agegroup, data = data, FUN = mean, na.rm = TRUE)
    mdr_summary$MDR_Positive <- aggregate(MDR ~ Year + State + Agegroup, data = data, FUN = sum, na.rm = TRUE)$MDR
    mdr_summary$Total_Isolates <- aggregate(MDR ~ Year + State + Agegroup, data = data, FUN = length)$MDR
  } else {
    mdr_summary <- aggregate(MDR ~ Year + Agegroup, data = data, FUN = mean, na.rm = TRUE)
    mdr_summary$MDR_Positive <- aggregate(MDR ~ Year + Agegroup, data = data, FUN = sum, na.rm = TRUE)$MDR
    mdr_summary$Total_Isolates <- aggregate(MDR ~ Year + Agegroup, data = data, FUN = length)$MDR
  }
  
  # Convert year to factor
  mdr_summary$Year <- as.factor(mdr_summary$Year)
  
  # Return the MDR summary
  return(mdr_summary)
}

state_resistance_carbapAgeg<-compute_resistance_summaryAgeg(data_atlas_us_enterob, carbapenems, "state")
state_resistance_cephalosAgeg<-compute_resistance_summaryAgeg(data_atlas_us_enterob, third_gen_cephalosporins, "state")
state_resistance_firstlineAgeg<-compute_resistance_summaryAgeg(data_atlas_us_enterob, first_line_antibiotics , "state")
state_resistance_mdrAgeg<- compute_mdr_summaryAgeg(data_atlas_us_enterob, antibiotic_classes, level = "state")

state_resistance_carbap_oAgeg<-compute_resistance_summaryAgeg(data_atlas_us_enterob_orig, carbapenems, "state")
state_resistance_cephalos_oAgeg<-compute_resistance_summaryAgeg(data_atlas_us_enterob_orig, third_gen_cephalosporins, "state")
state_resistance_firstline_oAgeg<-compute_resistance_summaryAgeg(data_atlas_us_enterob_orig, first_line_antibiotics , "state")
state_resistance_mdr_oAgeg<- compute_mdr_summaryAgeg(data_atlas_us_enterob_orig, antibiotic_classes, level = "state")


compute_resistance_summary_euAgeg <- function(data, antibiotics, level = "country") {
  # Select relevant columns
  variables_to_keep <- c("Year", "Country", "Agegroup", antibiotics)
  data_selected <- data[, variables_to_keep]
  
  # Melt data to long format
  data_long <- reshape2::melt(data_selected, id.vars = c("Year", "Country", "Agegroup"), 
                              measure.vars = antibiotics, 
                              variable.name = "Antibiotic", 
                              value.name = "Count")
  
  # Summarize data to get counts of non-missing values per antibiotic per year, country, and gender
  data_summary <- reshape2::dcast(data_long, Antibiotic ~ Year + Country + Agegroup, 
                                  fun.aggregate = function(x) sum(!is.na(x)))
  # Display summary table
  print(data_summary)
  
  # Select additional columns for more detailed analysis
  data_detailed <- subset(data, select = c("Year", "Country", "Gender", "Agegroup", "Species", antibiotics))
  
  # Calculate missing data
  missing_data_summary <- colSums(is.na(data_detailed[, antibiotics]))
  total_rows <- nrow(data_detailed)
  missing_data_proportion <- missing_data_summary / total_rows * 100
  missing_data_report <- data.frame(Antibiotic = antibiotics, 
                                    Missing_Count = missing_data_summary, 
                                    Missing_Proportion = missing_data_proportion)
  print(missing_data_report)
  
  # Convert antibiotics to numeric and calculate resistance
  data_detailed[antibiotics] <- lapply(data_detailed[antibiotics], as.numeric)
  data_detailed$Antibiotic_Sum <- rowSums(data_detailed[, antibiotics], na.rm = TRUE)
  data_detailed$Resistance <- ifelse(data_detailed$Antibiotic_Sum >= 1, 1, 0)
  
  # Compute counts of AMR positive isolates and all isolates tested per year, country, and gender
  if (level == "country") {
    # Calculate resistance percentage, AMR positives, and total isolates
    resistance_summary <- aggregate(Resistance ~ Year + Country + Agegroup, data = data_detailed, FUN = mean, na.rm = TRUE)
    resistance_summary$AMR_Positive <- aggregate(Resistance ~ Year + Country + Agegroup, data = data_detailed, FUN = sum, na.rm = TRUE)$Resistance
    resistance_summary$Total_Isolates <- aggregate(Resistance ~ Year + Country + Agegroup, data = data_detailed, FUN = length)$Resistance
  } else {
    resistance_summary <- aggregate(Resistance ~ Year + Agegroup, data = data_detailed, FUN = mean, na.rm = TRUE)
    resistance_summary$AMR_Positive <- aggregate(Resistance ~ Year + Agegroup, data = data_detailed, FUN = sum, na.rm = TRUE)$Resistance
    resistance_summary$Total_Isolates <- aggregate(Resistance ~ Year + Agegroup, data = data_detailed, FUN = length)$Resistance
  }
  
  # Convert year to factor
  resistance_summary$Year <- as.factor(resistance_summary$Year)
  
  # Return the resistance summary
  return(resistance_summary)
}
compute_mdr_summary_euAgeg <- function(data, antibiotic_classes, level = "country") {
  # Flatten the list of antibiotic classes into a vector of antibiotic names
  antibiotic_list <- unlist(antibiotic_classes)
  
  # Convert each antibiotic variable in the dataframe to numeric
  data[antibiotic_list] <- lapply(data[antibiotic_list], as.numeric)
  
  # Calculate resistance for each antibiotic class
  for (class_name in names(antibiotic_classes)) {
    class_antibiotics <- antibiotic_classes[[class_name]]
    data[[class_name]] <- rowSums(data[, class_antibiotics], na.rm = TRUE)
    data[[class_name]] <- ifelse(data[[class_name]] >= 1, 1, 0)
  }
  
  # Calculate the sum of the antibiotic class variables for MDR determination
  data$Total_Antibiotic_Class_Sum <- rowSums(data[, names(antibiotic_classes)], na.rm = TRUE)
  data$MDR <- ifelse(data$Total_Antibiotic_Class_Sum >= 3, 1, 0)
  
  # Aggregate MDR data by year or state level
  if (level == "country") {
    # Calculate MDR percentage, MDR positives, and total isolates
    mdr_summary <- aggregate(MDR ~ Year + Country + Agegroup, data = data, FUN = mean, na.rm = TRUE)
    mdr_summary$MDR_Positive <- aggregate(MDR ~ Year + Country + Agegroup, data = data, FUN = sum, na.rm = TRUE)$MDR
    mdr_summary$Total_Isolates <- aggregate(MDR ~ Year + Country + Agegroup, data = data, FUN = length)$MDR
  } else {
    mdr_summary <- aggregate(MDR ~ Year + Agegroup, data = data, FUN = mean, na.rm = TRUE)
    mdr_summary$MDR_Positive <- aggregate(MDR ~ Year + Agegroup, data = data, FUN = sum, na.rm = TRUE)$MDR
    mdr_summary$Total_Isolates <- aggregate(MDR ~ Year + Agegroup, data = data, FUN = length)$MDR
  }
  
  # Convert year to factor
  mdr_summary$Year <- as.factor(mdr_summary$Year)
  
  # Return the MDR summary
  return(mdr_summary)
}

country_resistance_carbap_euAgeg<-compute_resistance_summary_euAgeg(data_atlas_eu_enterob, carbapenems, "country")
country_resistance_cephalos_euAgeg<-compute_resistance_summary_euAgeg(data_atlas_eu_enterob, third_gen_cephalosporins, "country")
country_resistance_firstline_euAgeg<-compute_resistance_summary_euAgeg(data_atlas_eu_enterob, first_line_antibiotics , "country")
country_resistance_mdr_euAgeg<- compute_mdr_summary_euAgeg(data_atlas_eu_enterob, antibiotic_classes, level = "country")

country_resistance_carbap_o_euAgeg<-compute_resistance_summary_euAgeg(data_atlas_eu_enterob_orig, carbapenems, "country")
country_resistance_cephalos_o_euAgeg<-compute_resistance_summary_euAgeg(data_atlas_eu_enterob_orig, third_gen_cephalosporins, "country")
country_resistance_firstline_o_euAgeg<-compute_resistance_summary_euAgeg(data_atlas_eu_enterob_orig, first_line_antibiotics , "country")
country_resistance_mdr_o_euAgeg<- compute_mdr_summary_euAgeg(data_atlas_eu_enterob_orig, antibiotic_classes, level = "country")




######

#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#GRAPHS OF TRENDS PER STATE and PER EUROPEAN COUNTRY, across antibiotic groups.

#####US graphs ######
##CARBAPENEM resistance in the US.
# Aggregate Carbapenem resistance by Year and State
carbapenem_resistance_by_year_state <- aggregate(Carbapenem_res ~ Year + State, data = data_atlas_us_enterob_carbap, FUN = mean, na.rm = TRUE)
carbapenem_resistance_by_year_state$Carbapenem_res <- carbapenem_resistance_by_year_state$Carbapenem_res * 100
carbapenem_resistance_by_year_state$Year <- as.factor(carbapenem_resistance_by_year_state$Year)
carbapenem_resistance_by_year_state <- carbapenem_resistance_by_year_state[order(carbapenem_resistance_by_year_state$Year, carbapenem_resistance_by_year_state$State),]
carbapenem_resistance_by_year_state$Year <- as.factor(carbapenem_resistance_by_year_state$Year)
# Plotting
# Base color (Choose an appropriate color, here's an example with blue)
base_color <- "#FC6C85"
  # Function to darken color
darken_color <- function(color, amount = 0.2) {
  scales::colour_ramp(c(color, "#FC6C85"))(amount)[1]
}
# Apply the function
dot_color <- darken_color(base_color, 0.2)
# Creating the plot
p <- ggplot(carbapenem_resistance_by_year_state, aes(x = Year, y = Carbapenem_res, group = State)) +
  geom_line(aes(color = base_color), size = 1) +  # Line with base color
  geom_point(aes(color = dot_color), size = 1.2, shape = 20, fill = dot_color, stroke = 1) + # Dots, slightly darker, with black contour
  scale_color_identity() + # Use the colors as is
  facet_wrap(~ State, scales = "free_y") + # Facet by State, with independent y scales
  labs(title = "",
       x = "Year",
       y = "Carbapenem-resistance in Enterobacterales (%)") +
  theme_minimal(base_size = 12) +  # Start with a minimal theme
  theme(
    text = element_text(family = "serif"), # Lancet uses serif fonts
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.background = element_blank(), # Remove facet label backgrounds
    strip.text.x = element_text(size = 10, face = "bold"), # Style facet labels
    panel.grid.major = element_blank(), # No major grid lines
    panel.grid.minor = element_blank(), # No minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5), # Black border
    axis.text.x = element_text(angle = 90, hjust = 1, size=7), # Slant x axis texts for readability
    axis.text.y = element_text(size=8), # Slant x axis texts for readability
    legend.position = "none" # No legend needed
  )
ggsave(filename = "DStrebd_state_resistance_carbap.tiff", plot = p, device = "tiff", path = base_pathOut,
       width = 13, height = 7, dpi = 500, units = "in")

#Cephalosporin-resistance in the US states:
# Aggregate Carbapenem resistance by Year and State
cephalosporin_resistance_by_year_state <- aggregate(Cephalosporin_res ~ Year + State, data = data_atlas_us_enterob_cepha, FUN = mean, na.rm = TRUE)
cephalosporin_resistance_by_year_state$Cephalosporin_res <- cephalosporin_resistance_by_year_state$Cephalosporin_res * 100
cephalosporin_resistance_by_year_state$Year <- as.factor(cephalosporin_resistance_by_year_state$Year)
cephalosporin_resistance_by_year_state <- cephalosporin_resistance_by_year_state[order(cephalosporin_resistance_by_year_state$Year, cephalosporin_resistance_by_year_state$State),]
cephalosporin_resistance_by_year_state$Year <- as.factor(cephalosporin_resistance_by_year_state$Year)
# Plotting
# Base color (Choose an appropriate color, here's an example with blue)
base_color <- "#FC6C85"
  # Function to darken color
darken_color <- function(color, amount = 0.2) {
  scales::colour_ramp(c(color, "#FC6C85"))(amount)[1]
}
# Apply the function
dot_color <- darken_color(base_color, 0.2)
# Creating the plot
p <- ggplot(cephalosporin_resistance_by_year_state, aes(x = Year, y = Cephalosporin_res, group = State)) +
  geom_line(aes(color = base_color), size = 1) +  # Line with base color
  geom_point(aes(color = dot_color), size = 1.2, shape = 20, fill = dot_color, stroke = 1) + # Dots, slightly darker, with black contour
  scale_color_identity() + # Use the colors as is
  facet_wrap(~ State, scales = "free_y") + # Facet by State, with independent y scales
  labs(title = "",
       x = "Year",
       y = "3rd generation cephalosporin-resistance in Enterobacterales (%)") +
  theme_minimal(base_size = 12) +  # Start with a minimal theme
  theme(
    text = element_text(family = "serif"), # Lancet uses serif fonts
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.background = element_blank(), # Remove facet label backgrounds
    strip.text.x = element_text(size = 10, face = "bold"), # Style facet labels
    panel.grid.major = element_blank(), # No major grid lines
    panel.grid.minor = element_blank(), # No minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5), # Black border
    axis.text.x = element_text(angle = 90, hjust = 1, size=7), # Slant x axis texts for readability
    axis.text.y = element_text(size=8), # Slant x axis texts for readability
    legend.position = "none" # No legend needed
  )
ggsave(filename = "DStrebd_state_resistance_cephalos.tiff", plot = p, device = "tiff", path = base_pathOut,
       width = 13, height = 7, dpi = 500, units = "in")

#First-line antibiotic -resistance in the US states:
# Aggregate Carbapenem resistance by Year and State
firstlin_resistance_by_year_state <- aggregate(Firstlin_res ~ Year + State, data = data_atlas_us_enterob_firstlin, FUN = mean, na.rm = TRUE)
firstlin_resistance_by_year_state$firstlin_res <- firstlin_resistance_by_year_state$Firstlin_res * 100
firstlin_resistance_by_year_state$Year <- as.factor(firstlin_resistance_by_year_state$Year)
firstlin_resistance_by_year_state <- firstlin_resistance_by_year_state[order(firstlin_resistance_by_year_state$Year, firstlin_resistance_by_year_state$State),]
firstlin_resistance_by_year_state$Year <- as.factor(firstlin_resistance_by_year_state$Year)
# Plotting
# Base color (Choose an appropriate color, here's an example with blue)
base_color <- "#FC6C85"
  # Function to darken color
darken_color <- function(color, amount = 0.2) {
  scales::colour_ramp(c(color, "#FC6C85"))(amount)[1]
}
# Apply the function
dot_color <- darken_color(base_color, 0.2)
# Creating the plot
p <- ggplot(firstlin_resistance_by_year_state, aes(x = Year, y = Firstlin_res, group = State)) +
  geom_line(aes(color = base_color), size = 1) +  # Line with base color
  geom_point(aes(color = dot_color), size = 1.2, shape = 20, fill = dot_color, stroke = 1) + # Dots, slightly darker, with black contour
  scale_color_identity() + # Use the colors as is
  facet_wrap(~ State, scales = "free_y") + # Facet by State, with independent y scales
  labs(title = "",
       x = "Year",
       y = "First line antibiotic-resistance in Enterobacterales (%)") +
  theme_minimal(base_size = 12) +  # Start with a minimal theme
  theme(
    text = element_text(family = "serif"), # Lancet uses serif fonts
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.background = element_blank(), # Remove facet label backgrounds
    strip.text.x = element_text(size = 10, face = "bold"), # Style facet labels
    panel.grid.major = element_blank(), # No major grid lines
    panel.grid.minor = element_blank(), # No minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5), # Black border
    axis.text.x = element_text(angle = 90, hjust = 1, size=7), # Slant x axis texts for readability
    axis.text.y = element_text(size=8), # Slant x axis texts for readability
    legend.position = "none" # No legend needed
  )
ggsave(filename = "DStrebd_state_resistance_firstlin.tiff", plot = p, device = "tiff", path = base_pathOut,
       width = 13, height = 7, dpi = 500, units = "in")



#First-line antibiotic -resistance in the US states:
# Aggregate Carbapenem resistance by Year and State
mdr_resistance_by_year_state <- aggregate(MDR ~ Year + State, data = data_atlas_us, FUN = mean, na.rm = TRUE)
mdr_resistance_by_year_state$MDR <- mdr_resistance_by_year_state$MDR * 100
mdr_resistance_by_year_state$Year <- as.factor(mdr_resistance_by_year_state$Year)
mdr_resistance_by_year_state <- mdr_resistance_by_year_state[order(mdr_resistance_by_year_state$Year, mdr_resistance_by_year_state$State),]
mdr_resistance_by_year_state$Year <- as.factor(mdr_resistance_by_year_state$Year)
# Plotting
# Base color (Choose an appropriate color, here's an example with blue)
base_color <- "#FC6C85"
  # Function to darken color
darken_color <- function(color, amount = 0.2) {
  scales::colour_ramp(c(color, "#FC6C85"))(amount)[1]
}
# Apply the function
dot_color <- darken_color(base_color, 0.2)
# Creating the plot
p <- ggplot(mdr_resistance_by_year_state, aes(x = Year, y = MDR, group = State)) +
  geom_line(aes(color = base_color), size = 1) +  # Line with base color
  geom_point(aes(color = dot_color), size = 1.2, shape = 20, fill = dot_color, stroke = 1) + # Dots, slightly darker, with black contour
  scale_color_identity() + # Use the colors as is
  facet_wrap(~ State, scales = "free_y") + # Facet by State, with independent y scales
  labs(title = "",
       x = "Year",
       y = "MDR in Enterobacterales (%)") +
  theme_minimal(base_size = 12) +  # Start with a minimal theme
  theme(
    text = element_text(family = "serif"), # Lancet uses serif fonts
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.background = element_blank(), # Remove facet label backgrounds
    strip.text.x = element_text(size = 10, face = "bold"), # Style facet labels
    panel.grid.major = element_blank(), # No major grid lines
    panel.grid.minor = element_blank(), # No minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5), # Black border
    axis.text.x = element_text(angle = 90, hjust = 1, size=7), # Slant x axis texts for readability
    axis.text.y = element_text(size=8), # Slant x axis texts for readability
    legend.position = "none" # No legend needed
  )
ggsave(filename = "DStrebd_state_resistance_mdr.tiff", plot = p, device = "tiff", path = base_pathOut,
       width = 13, height = 7, dpi = 500, units = "in")
######

#####EUROPE graphs ########
##CARBAPENEM resistance in EUROPE
# Aggregate Carbapenem resistance by Year and State
carbapenem_resistance_by_year_country_eu <- aggregate(Carbapenem_res ~ Year + Country, data = data_atlas_eu_enterob_carbap, FUN = mean, na.rm = TRUE)
carbapenem_resistance_by_year_country_eu$Carbapenem_res <- carbapenem_resistance_by_year_country_eu$Carbapenem_res * 100
carbapenem_resistance_by_year_country_eu$Year <- as.factor(carbapenem_resistance_by_year_country_eu$Year)
carbapenem_resistance_by_year_country_eu <- carbapenem_resistance_by_year_country_eu[order(carbapenem_resistance_by_year_country_eu$Year, carbapenem_resistance_by_year_country_eu$Country),]
carbapenem_resistance_by_year_country_eu$Year <- as.factor(carbapenem_resistance_by_year_country_eu$Year)
carbapenem_resistance_by_year_country_eu<-carbapenem_resistance_by_year_country_eu[carbapenem_resistance_by_year_country_eu$Country %in% countries_list_eu_includ, ]

# Plotting
# Base color (Choose an appropriate color, here's an example with blue)
base_color <- "#377eb8"
  # Function to darken color
darken_color <- function(color, amount = 0.2) {
  scales::colour_ramp(c(color, "#377eb8"))(amount)[1]
}
# Apply the function
dot_color <- darken_color(base_color, 0.2)
# Creating the plot
p <- ggplot(carbapenem_resistance_by_year_country_eu, aes(x = Year, y = Carbapenem_res, group = Country)) +
  geom_line(aes(color = base_color), size = 1) +  # Line with base color
  geom_point(aes(color = dot_color), size = 1.2, shape = 20, fill = dot_color, stroke = 1) + # Dots, slightly darker, with black contour
  scale_color_identity() + # Use the colors as is
  facet_wrap(~ Country, scales = "fixed") + # Facet by State, with independent y scales
  ylim(0, 40)+
  labs(title = "",
       x = "Year",
       y = "Carbapenem-resistance in Enterobacterales (%)") +
  theme_minimal(base_size = 12) +  # Start with a minimal theme
  theme(
    text = element_text(family = "serif"), # Lancet uses serif fonts
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.background = element_blank(), # Remove facet label backgrounds
    strip.text.x = element_text(size = 10, face = "bold"), # Style facet labels
    panel.grid.major = element_blank(), # No major grid lines
    panel.grid.minor = element_blank(), # No minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5), # Black border
    axis.text.x = element_text(angle = 90, hjust = 1, size=7), # Slant x axis texts for readability
    axis.text.y = element_text(size=8), # Slant x axis texts for readability
    legend.position = "none" # No legend needed
  )
ggsave(filename = "DStrebd_EUcount_resistance_carbap.tiff", plot = p, device = "tiff", path = base_pathOut,
       width = 13, height = 7, dpi = 500, units = "in")

#Cephalosporin-resistance in the US states:
# Aggregate Carbapenem resistance by Year and State
cephalosporin_resistance_by_year_country <- aggregate(Cephalosporin_res ~ Year + Country, data = data_atlas_eu_enterob_cepha, FUN = mean, na.rm = TRUE)
cephalosporin_resistance_by_year_country$Cephalosporin_res <- cephalosporin_resistance_by_year_country$Cephalosporin_res * 100
cephalosporin_resistance_by_year_country$Year <- as.factor(cephalosporin_resistance_by_year_country$Year)
cephalosporin_resistance_by_year_country <- cephalosporin_resistance_by_year_country[order(cephalosporin_resistance_by_year_country$Year, cephalosporin_resistance_by_year_country$Country),]
cephalosporin_resistance_by_year_country$Year <- as.factor(cephalosporin_resistance_by_year_country$Year)
cephalosporin_resistance_by_year_country<-cephalosporin_resistance_by_year_country[cephalosporin_resistance_by_year_country$Country %in% countries_list_eu_includ, ]
# Plotting
# Creating the plot
p <- ggplot(cephalosporin_resistance_by_year_country, aes(x = Year, y = Cephalosporin_res, group = Country)) +
  geom_line(aes(color = base_color), size = 1) +  # Line with base color
  geom_point(aes(color = dot_color), size = 1.2, shape = 20, fill = dot_color, stroke = 1) + # Dots, slightly darker, with black contour
  scale_color_identity() + # Use the colors as is
  facet_wrap(~ Country, scales = "fixed") + # Facet by State, with independent y scales
  labs(title = "",
       x = "Year",
       y = "3rd generation cephalosporin-resistance in Enterobacterales (%)") +
  theme_minimal(base_size = 12) +  # Start with a minimal theme
  theme(
    text = element_text(family = "serif"), # Lancet uses serif fonts
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.background = element_blank(), # Remove facet label backgrounds
    strip.text.x = element_text(size = 10, face = "bold"), # Style facet labels
    panel.grid.major = element_blank(), # No major grid lines
    panel.grid.minor = element_blank(), # No minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5), # Black border
    axis.text.x = element_text(angle = 90, hjust = 1, size=7), # Slant x axis texts for readability
    axis.text.y = element_text(size=8), # Slant x axis texts for readability
    legend.position = "none" # No legend needed
  )
ggsave(filename = "DStrebd_Eucountry_resistance_cephalos.tiff", plot = p, device = "tiff", path = base_pathOut,
       width = 13, height = 7, dpi = 500, units = "in")

#First-line antibiotic -resistance in the Europe:
# Aggregate Carbapenem resistance by Year and State
firstlin_resistance_by_year_country <- aggregate(Firstlin_res ~ Year + Country, data = data_atlas_eu_enterob_firstlin, FUN = mean, na.rm = TRUE)
firstlin_resistance_by_year_country$firstlin_res <- firstlin_resistance_by_year_country$Firstlin_res * 100
firstlin_resistance_by_year_country$Year <- as.factor(firstlin_resistance_by_year_country$Year)
firstlin_resistance_by_year_country <- firstlin_resistance_by_year_country[order(firstlin_resistance_by_year_country$Year, firstlin_resistance_by_year_country$Country),]
firstlin_resistance_by_year_country$Year <- as.factor(firstlin_resistance_by_year_country$Year)
firstlin_resistance_by_year_country<-firstlin_resistance_by_year_country[firstlin_resistance_by_year_country$Country %in% countries_list_eu_includ, ]

# Plotting

# Creating the plot
p <- ggplot(firstlin_resistance_by_year_country, aes(x = Year, y = Firstlin_res, group = Country)) +
  geom_line(aes(color = base_color), size = 1) +  # Line with base color
  geom_point(aes(color = dot_color), size = 1.2, shape = 20, fill = dot_color, stroke = 1) + # Dots, slightly darker, with black contour
  scale_color_identity() + # Use the colors as is
  facet_wrap(~ Country, scales = "fixed") + # Facet by State, with independent y scales
  labs(title = "",
       x = "Year",
       y = "First line antibiotic-resistance in Enterobacterales (%)") +
  theme_minimal(base_size = 12) +  # Start with a minimal theme
  theme(
    text = element_text(family = "serif"), # Lancet uses serif fonts
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.background = element_blank(), # Remove facet label backgrounds
    strip.text.x = element_text(size = 10, face = "bold"), # Style facet labels
    panel.grid.major = element_blank(), # No major grid lines
    panel.grid.minor = element_blank(), # No minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5), # Black border
    axis.text.x = element_text(angle = 90, hjust = 1, size=7), # Slant x axis texts for readability
    axis.text.y = element_text(size=8), # Slant x axis texts for readability
    legend.position = "none" # No legend needed
  )
ggsave(filename = "DStrebd_EUcounty_resistance_firstlin.tiff", plot = p, device = "tiff", path = base_pathOut,
       width = 13, height = 7, dpi = 500, units = "in")



#MDR
mdr_resistance_by_year_country <- aggregate(MDR ~ Year + Country, data = data_atlas_eu, FUN = mean, na.rm = TRUE)
mdr_resistance_by_year_country$MDR <- mdr_resistance_by_year_country$MDR * 100
mdr_resistance_by_year_country$Year <- as.factor(mdr_resistance_by_year_country$Year)
mdr_resistance_by_year_country <- mdr_resistance_by_year_country[order(mdr_resistance_by_year_country$Year, mdr_resistance_by_year_country$Country),]
mdr_resistance_by_year_country$Year <- as.factor(mdr_resistance_by_year_country$Year)
mdr_resistance_by_year_country<-mdr_resistance_by_year_country[mdr_resistance_by_year_country$Country %in% countries_list_eu_includ, ]

# Plotting
# Creating the plot
p <- ggplot(mdr_resistance_by_year_country, aes(x = Year, y = MDR, group = Country)) +
  geom_line(aes(color = base_color), size = 1) +  # Line with base color
  geom_point(aes(color = dot_color), size = 1.2, shape = 20, fill = dot_color, stroke = 1) + # Dots, slightly darker, with black contour
  scale_color_identity() + # Use the colors as is
  facet_wrap(~ Country, scales = "fixed") + # Facet by State, with independent y scales
  labs(title = "",
       x = "Year",
       y = "MDR in Enterobacterales (%)") +
  theme_minimal(base_size = 12) +  # Start with a minimal theme
  theme(
    text = element_text(family = "serif"), # Lancet uses serif fonts
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.background = element_blank(), # Remove facet label backgrounds
    strip.text.x = element_text(size = 10, face = "bold"), # Style facet labels
    panel.grid.major = element_blank(), # No major grid lines
    panel.grid.minor = element_blank(), # No minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5), # Black border
    axis.text.x = element_text(angle = 90, hjust = 1, size=7), # Slant x axis texts for readability
    axis.text.y = element_text(size=8), # Slant x axis texts for readability
    legend.position = "none" # No legend needed
  )
ggsave(filename = "DStrebd_EUcountr_resistance_mdr.tiff", plot = p, device = "tiff", path = base_pathOut,
       width = 13, height = 7, dpi = 500, units = "in")
######

#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#Figure 2 violin plots for US states and EUROPE
#######

# Function to calculate 3-year averages, adjusting the first group to 4 years
calculate_averages <- function(df, location_var, resistance_var) {
  # Ensure Year is properly formatted and numeric
  df <- df %>%
    mutate(Year = as.numeric(as.character(Year)),
           Year_group = cut(Year, 
                            breaks = c(2003, 2007, 2010, 2013, 2016, 2019, 2022), 
                            labels = c("2004-2007", "2008-2010", "2011-2013", "2014-2016", "2017-2019", "2020-2022")))
  # Calculate averages by location and year group
  df %>%
    filter(!is.na(Year_group)) %>%
    group_by_at(c(location_var, "Year_group")) %>%
    summarise(mean_resistance = mean({{resistance_var}}, na.rm = TRUE), .groups = 'drop') %>%
    ungroup()
}

# Apply the function to your datasets
carbapenem_res_us_avg <- calculate_averages(carbapenem_resistance_by_year_state, "State", Carbapenem_res)
cephalosporin_res_us_avg <- calculate_averages(cephalosporin_resistance_by_year_state, "State", Cephalosporin_res)
mdr_res_us_avg <- calculate_averages(mdr_resistance_by_year_state, "State", MDR)
firstlin_res_us_avg <- calculate_averages(firstlin_resistance_by_year_state, "State", Firstlin_res)

carbapenem_res_eu_avg <- calculate_averages(carbapenem_resistance_by_year_country_eu, "Country", Carbapenem_res)
cephalosporin_res_eu_avg <- calculate_averages(cephalosporin_resistance_by_year_country, "Country", Cephalosporin_res)
mdr_res_eu_avg <- calculate_averages(mdr_resistance_by_year_country, "Country", MDR)
firstlin_res_eu_avg <- calculate_averages(firstlin_resistance_by_year_country, "Country", Firstlin_res)

# Function to combine datasets for plotting with specific attention to the names
combine_for_plot <- function(us_data, eu_data, resistance_type) {
  # Ensure column names are set correctly
  us_data$Region <- 'US'
  eu_data$Region <- 'Europe'
  us_data$Resistance_Type <- resistance_type
  eu_data$Resistance_Type <- resistance_type
  
  # Select and order columns explicitly to avoid mismatch error
  required_columns <- c("Region", "Year_group", "mean_resistance", "Resistance_Type")
  us_data <- us_data[, required_columns, drop = FALSE]  # Ensure column presence and order
  eu_data <- eu_data[, required_columns, drop = FALSE]  # Ensure column presence and order
  
  combined <- rbind(us_data, eu_data)
  combined
}

# Combine data for plotting
combined_data <- bind_rows(
  combine_for_plot(carbapenem_res_us_avg, carbapenem_res_eu_avg, "Carbapenem"),
  combine_for_plot(cephalosporin_res_us_avg, cephalosporin_res_eu_avg, "Cephalosporin"),
  combine_for_plot(mdr_res_us_avg, mdr_res_eu_avg, "MDR"),
  combine_for_plot(firstlin_res_us_avg, firstlin_res_eu_avg, "First-line")
)

# Define a color palette
color_palette <- c("US" = "#1b9e77", "Europe" = "#d95f02")
# Function to create a plot for a given resistance type with properly aligned boxplots
create_plotC <- function(data, resistance_type) {
  plot <- ggplot(data[data$Resistance_Type == resistance_type, ], aes(x = Year_group, y = mean_resistance, fill = Region, color = Region)) +
    geom_violin(trim = TRUE, alpha = 0.6, position = position_dodge(width = 0.8)) +
    geom_boxplot(aes(group = interaction(Region, Year_group)), width = 0.2, position = position_dodge(width = 0.8), outlier.shape = NA, color = "black", fill = "white", alpha = 0.5) +
    geom_point(position = position_dodge(width = 0.8), size = 1, alpha = 0.5, aes(color = Region), shape = 16) +  # Add points with lighter colors
    scale_fill_manual(values = color_palette) +
    labs(title = "", y = "Average carbapenem-resistance (%)", x = "") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 360, hjust = 1),
      strip.text.x = element_text(size = 14, face = "bold"),
      legend.position = "none",
      axis.line = element_line(color = "black", size = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", color = "black"),
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_text(face = "bold", size = 12)
    ) +
    scale_y_continuous(limits = c(0, 30),breaks = seq(0, 30, by = 5), oob = scales::oob_squish)  # Adjust upper limit dynamically
  
  return(plot)
}
# Create a plot for "Carbapenem" resistance type
plot_carbapenem <- create_plotC(combined_data, "Carbapenem")

create_plotCe <- function(data, resistance_type) {
  plot <- ggplot(data[data$Resistance_Type == resistance_type, ], aes(x = Year_group, y = mean_resistance, fill = Region, color = Region)) +
    geom_violin(trim = TRUE, alpha = 0.6, position = position_dodge(width = 0.8)) +
    geom_boxplot(aes(group = interaction(Region, Year_group)), width = 0.2, position = position_dodge(width = 0.8), outlier.shape = NA, color = "black", fill = "white", alpha = 0.5) +
    geom_point(position = position_dodge(width = 0.8), size = 1, alpha = 0.5, aes(color = Region), shape = 16) +  # Add points with lighter colors
    scale_fill_manual(values = color_palette) +
    labs(title = "", y = "Average 3rd generation cephalosporin-resistance (%)", x = "Year group") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 360, hjust = 1),
      strip.text.x = element_text(size = 14, face = "bold"),
      legend.position = "none",
      axis.line = element_line(color = "black", size = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", color = "black"),
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_text(face = "bold", size = 12)
    ) +
    scale_y_continuous(limits = c(0, 70),breaks = seq(0, 70, by = 5), oob = scales::oob_squish)  # Adjust upper limit dynamically
  
  return(plot)
}
# Create a plot for "Carbapenem" resistance type
plot_cephalosporin<- create_plotCe(combined_data, "Cephalosporin")


create_plotFL <- function(data, resistance_type) {
  plot <- ggplot(data[data$Resistance_Type == resistance_type, ], aes(x = Year_group, y = mean_resistance*100, fill = Region, color = Region)) +
    geom_violin(trim = TRUE, alpha = 0.6, position = position_dodge(width = 0.8)) +
    geom_boxplot(aes(group = interaction(Region, Year_group)), width = 0.2, position = position_dodge(width = 0.8), outlier.shape = NA, color = "black", fill = "white", alpha = 0.5) +
    geom_point(position = position_dodge(width = 0.8), size = 1, alpha = 0.5, shape = 16) +
    scale_fill_manual(values = color_palette) +
    labs(title = "", y = "Average first-line antibiotic-resistance (%)", x = "") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 360, hjust = 1),
      strip.text.x = element_text(size = 14, face = "bold"),
      legend.position = c(0.95, 0.2),  # Place legend inside the plot area at the top right
      legend.justification = c(1, 1),  # Anchor the legend at the top right
      legend.box.just = "right",  # Justify the legend box at the right
      legend.text = element_text(size = 11),
      #legend.background = element_rect(fill = "none", colour = "none"),  # Optional: Style the legend background
      axis.line = element_line(color = "black", size = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", color = "black"),
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_text(face = "bold", size = 12)
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10), oob = scales::oob_squish)  # Adjust upper limit dynamically
  
  return(plot)
}
# Create a plot for "Carbapenem" resistance type
plot_firstline<- create_plotFL(combined_data, "First-line")

create_plotmdr <- function(data, resistance_type) {
  plot <- ggplot(data[data$Resistance_Type == resistance_type, ], aes(x = Year_group, y = mean_resistance, fill = Region, color = Region)) +
    geom_violin(trim = TRUE, alpha = 0.6, position = position_dodge(width = 0.8)) +
    geom_boxplot(aes(group = interaction(Region, Year_group)), width = 0.2, position = position_dodge(width = 0.8), outlier.shape = NA, color = "black", fill = "white", alpha = 0.5) +
    geom_point(position = position_dodge(width = 0.8), size = 1, alpha = 0.5, aes(color = Region), shape = 16) +  # Add points with lighter colors
    scale_fill_manual(values = color_palette) +
    labs(title = "", y = "Average MDR (%)", x = "Year group") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 360, hjust = 1),
      strip.text.x = element_text(size = 14, face = "bold"),
      legend.position = "none",
      axis.line = element_line(color = "black", size = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", color = "black"),
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_text(face = "bold", size = 12)
    ) +
    scale_y_continuous(limits = c(0, 70),breaks = seq(0, 70, by = 5), oob = scales::oob_squish)  # Adjust upper limit dynamically
  
  return(plot)
}
# Create a plot for "Carbapenem" resistance type
plot_mdr<- create_plotmdr(combined_data, "MDR")

# Combine the plots into a 2x2 layout
combined_plot <- (plot_firstline | plot_carbapenem) / 
  (plot_cephalosporin | plot_mdr)  # Use | for side-by-side and / for above-below

# Use patchwork to arrange the plots and add tagging
combined_plotxoxx <- combined_plot + 
  plot_layout(guides = 'collect') +  # Collects and merges the guides
  plot_annotation(tag_levels = 'A')  # Adds tags A, B, C, D to each plot

combined_plotxoxx <- combined_plotxoxx & theme(
  plot.tag = element_text(face = "bold", size = 14)  # Set tags in bold and adjust size if needed
)

ggsave(filename = "TrendsEuropeandStates_amr.tiff", plot = combined_plotxoxx, device = "tiff", path = base_pathOut,
       width = 16, height = 10, dpi = 800, units = "in") 

#######

#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#Table 1 geographical area statistics ########
#state_resistance_carbap state_resistance_cephalos state_resistance_firstline  country_resistance_carbap_eu country_resistance_cephalos_eu country_resistance_firstline_eu
#country_resistance_mdr_eu state_resistance_mdr
######
process_data <- function(data) {
  data %>%
    group_by(Year) %>%
    summarise(
      Median_Resistance = median(Resistance, na.rm = TRUE),
      P25_Resistance = quantile(Resistance, 0.25, na.rm = TRUE),  # 25th percentile
      P75_Resistance = quantile(Resistance, 0.75, na.rm = TRUE),  # 75th percentile
      Total_Isolates = sum(Total_Isolates, na.rm = TRUE)
    )
}

results_carbap_state <- process_data(state_resistance_carbap)
results_cephalos_state <- process_data(state_resistance_cephalos)
results_firstline_state <- process_data(state_resistance_firstline)
results_carbap_country <- process_data(country_resistance_carbap_eu)
results_cephalos_country <- process_data(country_resistance_cephalos_eu)
results_firstline_country <- process_data(country_resistance_firstline_eu)

# Set row names for merging
row.names(results_carbap_state) <- results_carbap_state$Year
row.names(results_cephalos_state) <- results_cephalos_state$Year
row.names(results_firstline_state) <- results_firstline_state$Year
row.names(results_carbap_country) <- results_carbap_country$Year
row.names(results_cephalos_country) <- results_cephalos_country$Year
row.names(results_firstline_country) <- results_firstline_country$Year

# Combine by binding columns
final_table <- cbind(
  Carbap_State = results_carbap_state,
  Cephalos_State = results_cephalos_state,
  Firstline_State = results_firstline_state,
  Carbap_Country = results_carbap_country,
  Cephalos_Country = results_cephalos_country,
  Firstline_Country = results_firstline_country
)

# Print the final table
file_name <- "resistance_summary_tableAMR.xlsx"
# Create the full path for the file
full_path <- paste0(base_pathOut, file_name)
# Use openxlsx to write the DataFrame to an Excel file
write.xlsx(final_table, file = full_path)


process_data2 <- function(data) {
  data %>%
    group_by(Year) %>%
    summarise(
      Median_Resistance = median(MDR, na.rm = TRUE),
      P25_Resistance = quantile(MDR, 0.25, na.rm = TRUE),  # 25th percentile
      P75_Resistance = quantile(MDR, 0.75, na.rm = TRUE),  # 75th percentile
      Total_Isolates = sum(Total_Isolates, na.rm = TRUE)
    )
}

results_mdr_state <- process_data2(state_resistance_mdr)
results_mdr_country <- process_data2(country_resistance_mdr_eu)

# Set row names for merging
row.names(results_mdr_state) <- results_mdr_state$Year
row.names(results_mdr_country) <- results_mdr_country$Year

# Combine by binding columns
final_table2 <- cbind(
  MDR_State = results_mdr_state,
  MDR_Country = results_mdr_country
)


# Specify the file name
file_name <- "resistance_summary_tableMDR.xlsx"
# Create the full path for the file
full_path <- paste0(base_pathOut, file_name)
# Use openxlsx to write the DataFrame to an Excel file
write.xlsx(final_table2, file = full_path)





######

#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
# Analyses for the US states GAM SPATIOTEMPORAL 
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#

######
#STATE INFO US:  state_resistance_carbap state_resistance_cephalos state_resistance_firstline state_resistance_mdr.   ///#COUNTRY INFO:  country_resistance_carbap country_resistance_cephalos country_resistance_firstline  country_resistance_mdr 
fit_gam_model <- function(state_resistance_carbap) {
  # 1. Load Spatial Data and Resistance Data
  states_shapefile <- states(cb = TRUE)
  # 2. Convert 'states_shapefile' to 'sf' object if not already
  if (!inherits(states_shapefile, "sf")) {
    states_shapefile <- st_as_sf(states_shapefile)
  }
  state_resistance_carbap$NAME<- state_resistance_carbap$State
  # 3. Merge spatial data with resistance data by 'NAME'
  states_merged <- merge(states_shapefile, state_resistance_carbap, by = "NAME", all.x = TRUE)
  # 4. Convert merged data to 'Spatial' object for neighborhood creation
  states_spatial <- as(states_merged, "Spatial")
  # 5. Create unique geometries for neighborhood structure
  states_spatial_unique <- states_spatial[!duplicated(states_spatial$GEOID), ]
  # 6. Remove Hawaii (GEOID "15") and Alaska (GEOID "02") from the spatial data
  states_spatial_filtered <- states_spatial[!states_spatial$GEOID %in% c("15", "02"), ]
  # 7. Convert to 'Spatial' object if necessary
  if (!inherits(states_spatial_filtered, "Spatial")) {
    states_spatial_filtered <- as(states_spatial_filtered, "Spatial")
  }
  # 8. Remove duplicates based on GEOID to create unique state geometries for neighborhood structure
  states_spatial_unique <- states_spatial_filtered[!duplicated(states_spatial_filtered$GEOID), ]
  # 9. Recreate the neighborhood structure using unique state geometries
  nb <- poly2nb(states_spatial_unique, row.names = states_spatial_unique$GEOID)
  # 10. Assign region names to the neighborhood list
  names(nb) <- states_spatial_unique$GEOID
  print(nb)
  # 11. Clean and ensure data alignment with neighborhood structure
  # Convert 'GEOID' to factor and ensure levels match the neighborhood structure
  states_spatial_filtered@data$GEOID <- factor(states_spatial_filtered@data$GEOID, levels = attr(nb, "region.id"))
  # Convert 'Resistance' to numeric
  states_spatial_filtered@data$Resistance <- as.numeric(states_spatial_filtered@data$Resistance)
  states_spatial_filtered@data$AMR_Positive <- as.numeric(states_spatial_filtered@data$AMR_Positive)
  states_spatial_filtered@data$Total_Isolates <- as.numeric(states_spatial_filtered@data$Total_Isolates)
  # Convert 'Year' to numeric and handle missing data
  states_spatial_filtered@data$Year <- as.numeric(as.character(states_spatial_filtered@data$Year))
  states_spatial_filtered@data <- states_spatial_filtered@data[complete.cases(states_spatial_filtered@data), ]
  # 12. Recreate neighborhood structure using the cleaned and filtered data
  states_spatial_unique <- states_spatial_filtered[!duplicated(states_spatial_filtered$GEOID), ]
  nb <- poly2nb(states_spatial_unique, row.names = states_spatial_unique$GEOID)
  names(nb) <- states_spatial_unique$GEOID
  # 13. Set control parameters for GAM
  ctrl <- gam.control(nthreads = 6)  # Set parallel threads for faster computation
  # 14. Fit the MRF model with spatio-temporal smoothing
  model_gamx1 <- gam(Resistance ~ 
                       s(Year, m=3, k=10, bs = "tp") +   # Smooth term for Year
                       s(GEOID, bs = 'mrf', xt = list(nb = nb)) +  # Smooth term for GEOID
                       te(Year, GEOID, bs=c("tp", "mrf"), m=c(3, NA), xt=list(Year=NULL, GEOID=list(nb=nb))), # Interaction term
                     data = states_spatial_filtered@data,
                     family = binomial, 
                     select = TRUE,
                     method = "REML", 
                     weights = Total_Isolates)
  return(model_gamx1)
}

fit_gam_model_mdr <- function(state_resistance_carbap) {
  # 1. Load Spatial Data and Resistance Data
  states_shapefile <- states(cb = TRUE)
  # 2. Convert 'states_shapefile' to 'sf' object if not already
  if (!inherits(states_shapefile, "sf")) {
    states_shapefile <- st_as_sf(states_shapefile)
  }
  state_resistance_carbap$NAME<- state_resistance_carbap$State
  # 3. Merge spatial data with resistance data by 'NAME'
  states_merged <- merge(states_shapefile, state_resistance_carbap, by = "NAME", all.x = TRUE)
  # 4. Convert merged data to 'Spatial' object for neighborhood creation
  states_spatial <- as(states_merged, "Spatial")
  # 5. Create unique geometries for neighborhood structure
  states_spatial_unique <- states_spatial[!duplicated(states_spatial$GEOID), ]
  # 6. Remove Hawaii (GEOID "15") and Alaska (GEOID "02") from the spatial data
  states_spatial_filtered <- states_spatial[!states_spatial$GEOID %in% c("15", "02"), ]
  # 7. Convert to 'Spatial' object if necessary
  if (!inherits(states_spatial_filtered, "Spatial")) {
    states_spatial_filtered <- as(states_spatial_filtered, "Spatial")
  }
  # 8. Remove duplicates based on GEOID to create unique state geometries for neighborhood structure
  states_spatial_unique <- states_spatial_filtered[!duplicated(states_spatial_filtered$GEOID), ]
  # 9. Recreate the neighborhood structure using unique state geometries
  nb <- poly2nb(states_spatial_unique, row.names = states_spatial_unique$GEOID)
  # 10. Assign region names to the neighborhood list
  names(nb) <- states_spatial_unique$GEOID
  print(nb)
  # 11. Clean and ensure data alignment with neighborhood structure
  # Convert 'GEOID' to factor and ensure levels match the neighborhood structure
  states_spatial_filtered@data$GEOID <- factor(states_spatial_filtered@data$GEOID, levels = attr(nb, "region.id"))
  # Convert 'Resistance' to numeric
  states_spatial_filtered@data$MDR <- as.numeric(states_spatial_filtered@data$MDR)
  states_spatial_filtered@data$MDR_Positive <- as.numeric(states_spatial_filtered@data$MDR_Positive)
  states_spatial_filtered@data$Total_Isolates <- as.numeric(states_spatial_filtered@data$Total_Isolates)
  # Convert 'Year' to numeric and handle missing data
  states_spatial_filtered@data$Year <- as.numeric(as.character(states_spatial_filtered@data$Year))
  states_spatial_filtered@data <- states_spatial_filtered@data[complete.cases(states_spatial_filtered@data), ]
  # 12. Recreate neighborhood structure using the cleaned and filtered data
  states_spatial_unique <- states_spatial_filtered[!duplicated(states_spatial_filtered$GEOID), ]
  nb <- poly2nb(states_spatial_unique, row.names = states_spatial_unique$GEOID)
  names(nb) <- states_spatial_unique$GEOID
  # 13. Set control parameters for GAM
  ctrl <- gam.control(nthreads = 6)  # Set parallel threads for faster computation
  # 14. Fit the MRF model with spatio-temporal smoothing
  model_gamx1 <- gam(MDR ~ 
                       s(Year, m=3, k=10, bs = "tp") +   # Smooth term for Year
                       s(GEOID, bs = 'mrf', xt = list(nb = nb)) +  # Smooth term for GEOID
                       te(Year, GEOID, bs=c("tp", "mrf"), m=c(3, NA), xt=list(Year=NULL, GEOID=list(nb=nb))), # Interaction term
                     data = states_spatial_filtered@data,
                     family = binomial, 
                     select = TRUE,
                     method = "REML", 
                     weights = Total_Isolates)
  return(model_gamx1)
}

#OUTPUT OF THE RESULTS SPATIO-TEMPORAL MODELS:
model_output_carbap <- fit_gam_model(state_resistance_carbap)
model_output_cephalos <- fit_gam_model(state_resistance_cephalos)
model_output_firstline <- fit_gam_model(state_resistance_firstline)
model_output_mdr <-fit_gam_model_mdr(state_resistance_mdr)

# Tidying model outputs
extract_model_summary <- function(model) {
  # Extracting smooth terms
  smooth_summary <- tidy(model, parametric = FALSE)
  
  # Extracting parametric coefficients
  parametric_summary <- tidy(model, exponentiate = FALSE, parametric = TRUE)
  
  # Combine both summaries
  combined_summary <- bind_rows(parametric_summary, smooth_summary)
  
  return(combined_summary)
}

# Apply the function to each model
summary_carbap <- extract_model_summary(model_output_carbap)
summary_cephalos <- extract_model_summary(model_output_cephalos)
summary_firstline <- extract_model_summary(model_output_firstline)
summary_mdr <- extract_model_summary(model_output_mdr)

# Adding model identifiers
summary_carbap$model <- "Carbapenem Resistance"
summary_cephalos$model <- "Cephalosporin Resistance"
summary_firstline$model <- "First-line Antibiotic Resistance"
summary_mdr$model <- "MDR Resistance"

# Combining all summaries into one dataframe
combined_results <- bind_rows(summary_carbap, summary_cephalos, summary_firstline, summary_mdr)

# Optionally, select and rename columns for clarity
final_results <- combined_results %>%
  dplyr::select(Model = model, Term = term, Estimate = estimate, EDF = edf, Ref.DF = ref.df, Std.Error = std.error,
                Statistic = statistic, `P.Value` = p.value)
#####

#-------------------------------------------#
###FIRST DERIVATIVE GRAPHS, US states:
#-------------------------------------------#

#NEW, CARBAPENEM RESISTANCE GRAPH: ######
model_output_carbap <- fit_gam_model(state_resistance_carbap)
final_dataset<- model.frame(model_output_carbap)
original_geo_levels <- levels(final_dataset$GEOID)

years_range <- seq(min(state_resistance_carbap$Year), max(state_resistance_carbap$Year), length.out = 100)
# Create prediction data
new_data <- expand.grid(Year = years_range, GEOID = original_geo_levels)
new_data$Year <- as.numeric(as.character(new_data$Year))  # Ensure Year is nume

# Compute predictions and derivatives
preds <- predict(model_output_carbap, newdata = new_data, type = "terms", se.fit = TRUE, deriv = 1)

# Assuming you have already generated `preds` as before
# Now, include the interaction term in the derivative calculation
new_data$Year_deriv <- preds$fit[, "s(Year)"] + preds$fit[, "te(Year,GEOID)"]

# Calculate confidence intervals including both main and interaction effects
new_data$SE <- sqrt(preds$se.fit[, "s(Year)"]^2 + preds$se.fit[, "te(Year,GEOID)"]^2)

# Using a normal approximation for the confidence interval
alpha <- 0.05
z_score <- qnorm(1 - alpha / 2)
new_data$lower_ci <- new_data$Year_deriv - z_score * new_data$SE
new_data$upper_ci <- new_data$Year_deriv + z_score * new_data$SE

# Convert Spatial*DataFrame to a regular data frame for easier manipulation
states_data_df <- as.data.frame(states_spatial_filtered)

# Ensure GEOID is factor in both data frames for correct merging
new_data$GEOID <- as.factor(new_data$GEOID)
states_data_df$GEOID <- as.factor(states_data_df$GEOID)

# Merge to append the NAME corresponding to each GEOID
new_data <- merge(new_data, states_data_df[, c("GEOID", "NAME")], by = "GEOID", all.x = TRUE)

# Calculate zero crossing: Check if zero is inside or outside the interval
new_data <- new_data %>%
  mutate(
    zero_in_prev = ifelse(lag(lower_ci) <= 0 & lag(upper_ci) >= 0, 1, 0),
    zero_in_curr = ifelse(lower_ci <= 0 & upper_ci >= 0, 1, 0),
    sign_change = zero_in_prev != zero_in_curr  # Change in zero inclusion status
  )

# Filter crossings that happen after the year 2004
crossing_points <- new_data %>%
  filter(sign_change, Year > 2004) %>%
  dplyr::select(NAME, Year) %>%
  distinct(NAME, Year)  # Ensure unique crossing points

# Define the Lancet style theme
lancet_theme <- theme_minimal() +
  theme(text = element_text(family = "sans", color = "black"),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7, face = "bold"),
        legend.position = "none",  # Remove legend
        strip.text.x = element_text(size = 7))

# Plotting with NAME as labels and Lancet style, including multiple crossing points
p <- ggplot(new_data, aes(x = Year, y = Year_deriv, group = NAME, color = NAME)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = NAME), alpha = 0.2) +
  geom_vline(data = crossing_points, aes(xintercept = Year), color = "red", linetype = "dashed", size = 1) +
  facet_wrap(~NAME, scales = "free_y") +
  lancet_theme +
  labs(title = "",
       x = "Year",
       y = "Rate of change (first derivative from spatio-temporal GAM model)")
#First Derivative of Year by State in Carbapenem Resistance Model
ggsave(filename = "p_state_resistance_carbap.tiff", plot = p, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
######

#NEW, CEPHALOSPORIN RESISTANCE GRAPH: ######
model_output_cephalos <- fit_gam_model(state_resistance_cephalos)
final_dataset<- model.frame(model_output_cephalos)
original_geo_levels <- levels(final_dataset$GEOID)

years_range <- seq(min(state_resistance_cephalos$Year), max(state_resistance_cephalos$Year), length.out = 100)
# Create prediction data
new_data <- expand.grid(Year = years_range, GEOID = original_geo_levels)
new_data$Year <- as.numeric(as.character(new_data$Year))  # Ensure Year is nume

# Compute predictions and derivatives
preds <- predict(model_output_cephalos, newdata = new_data, type = "terms", se.fit = TRUE, deriv = 1)

# Assuming you have already generated `preds` as before
# Now, include the interaction term in the derivative calculation
new_data$Year_deriv <- preds$fit[, "s(Year)"] + preds$fit[, "te(Year,GEOID)"]

# Calculate confidence intervals including both main and interaction effects
new_data$SE <- sqrt(preds$se.fit[, "s(Year)"]^2 + preds$se.fit[, "te(Year,GEOID)"]^2)

# Using a normal approximation for the confidence interval
alpha <- 0.05
z_score <- qnorm(1 - alpha / 2)
new_data$lower_ci <- new_data$Year_deriv - z_score * new_data$SE
new_data$upper_ci <- new_data$Year_deriv + z_score * new_data$SE

# Convert Spatial*DataFrame to a regular data frame for easier manipulation
states_data_df <- as.data.frame(states_spatial_filtered)

# Ensure GEOID is factor in both data frames for correct merging
new_data$GEOID <- as.factor(new_data$GEOID)
states_data_df$GEOID <- as.factor(states_data_df$GEOID)

# Merge to append the NAME corresponding to each GEOID
new_data <- merge(new_data, states_data_df[, c("GEOID", "NAME")], by = "GEOID", all.x = TRUE)

# Calculate zero crossing: Check if zero is inside or outside the interval
new_data <- new_data %>%
  mutate(
    zero_in_prev = ifelse(lag(lower_ci) <= 0 & lag(upper_ci) >= 0, 1, 0),
    zero_in_curr = ifelse(lower_ci <= 0 & upper_ci >= 0, 1, 0),
    sign_change = zero_in_prev != zero_in_curr  # Change in zero inclusion status
  )

# Filter crossings that happen after the year 2004
crossing_points <- new_data %>%
  filter(sign_change, Year > 2004) %>%
  dplyr::select(NAME, Year) %>%
  distinct(NAME, Year)  # Ensure unique crossing points

# Define the Lancet style theme
lancet_theme <- theme_minimal() +
  theme(text = element_text(family = "sans", color = "black"),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7, face = "bold"),
        legend.position = "none",  # Remove legend
        strip.text.x = element_text(size = 7))

# Plotting with NAME as labels and Lancet style, including multiple crossing points
p <- ggplot(new_data, aes(x = Year, y = Year_deriv, group = NAME, color = NAME)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = NAME), alpha = 0.2) +
  geom_vline(data = crossing_points, aes(xintercept = Year), color = "red", linetype = "dashed", size = 1) +
  facet_wrap(~NAME, scales = "free_y") +
  lancet_theme +
  labs(title = "",
       x = "Year",
       y = "Rate of change (first derivative from spatio-temporal GAM model)")
#First Derivative of Year by State in Carbapenem Resistance Model
ggsave(filename = "p_state_resistance_cephalos.tiff", plot = p, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
######

#NEW, FIRST-LINE RESISTANCE GRAPH: ######
model_output_firstline <- fit_gam_model(state_resistance_firstline)
final_dataset<- model.frame(model_output_firstline)
original_geo_levels <- levels(final_dataset$GEOID)

state_resistance_firstline$Year <- as.numeric(as.character(state_resistance_firstline$Year))
years_range <- seq(min(state_resistance_firstline$Year), max(state_resistance_firstline$Year), length.out = 100)
# Create prediction data
new_data <- expand.grid(Year = years_range, GEOID = original_geo_levels)
new_data$Year <- as.numeric(as.character(new_data$Year))  # Ensure Year is nume

# Compute predictions and derivatives
preds <- predict(model_output_firstline, newdata = new_data, type = "terms", se.fit = TRUE, deriv = 1)

# Assuming you have already generated `preds` as before
# Now, include the interaction term in the derivative calculation
new_data$Year_deriv <- preds$fit[, "s(Year)"] + preds$fit[, "te(Year,GEOID)"]

# Calculate confidence intervals including both main and interaction effects
new_data$SE <- sqrt(preds$se.fit[, "s(Year)"]^2 + preds$se.fit[, "te(Year,GEOID)"]^2)

# Using a normal approximation for the confidence interval
alpha <- 0.05
z_score <- qnorm(1 - alpha / 2)
new_data$lower_ci <- new_data$Year_deriv - z_score * new_data$SE
new_data$upper_ci <- new_data$Year_deriv + z_score * new_data$SE

# Convert Spatial*DataFrame to a regular data frame for easier manipulation
states_data_df <- as.data.frame(states_spatial_filtered)

# Ensure GEOID is factor in both data frames for correct merging
new_data$GEOID <- as.factor(new_data$GEOID)
states_data_df$GEOID <- as.factor(states_data_df$GEOID)

# Merge to append the NAME corresponding to each GEOID
new_data <- merge(new_data, states_data_df[, c("GEOID", "NAME")], by = "GEOID", all.x = TRUE)

# Calculate zero crossing: Check if zero is inside or outside the interval
new_data <- new_data %>%
  mutate(
    zero_in_prev = ifelse(lag(lower_ci) <= 0 & lag(upper_ci) >= 0, 1, 0),
    zero_in_curr = ifelse(lower_ci <= 0 & upper_ci >= 0, 1, 0),
    sign_change = zero_in_prev != zero_in_curr  # Change in zero inclusion status
  )

# Filter crossings that happen after the year 2004
crossing_points <- new_data %>%
  filter(sign_change, Year > 2004) %>%
  dplyr::select(NAME, Year) %>%
  distinct(NAME, Year)  # Ensure unique crossing points

# Define the Lancet style theme
lancet_theme <- theme_minimal() +
  theme(text = element_text(family = "sans", color = "black"),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7, face = "bold"),
        legend.position = "none",  # Remove legend
        strip.text.x = element_text(size = 7))

# Plotting with NAME as labels and Lancet style, including multiple crossing points
p <- ggplot(new_data, aes(x = Year, y = Year_deriv, group = NAME, color = NAME)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = NAME), alpha = 0.2) +
  geom_vline(data = crossing_points, aes(xintercept = Year), color = "red", linetype = "dashed", size = 1) +
  facet_wrap(~NAME, scales = "free_y") +
  lancet_theme +
  labs(title = "",
       x = "Year",
       y = "Rate of change (first derivative from spatio-temporal GAM model)")
#First Derivative of Year by State in Carbapenem Resistance Model
ggsave(filename = "p_state_resistance_firstline.tiff", plot = p, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
######


#NEW, MDR RESISTANCE GRAPH: ######
model_output_mdr <- fit_gam_model_mdr(state_resistance_mdr)
final_dataset<- model.frame(model_output_mdr)
original_geo_levels <- levels(final_dataset$GEOID)

state_resistance_mdr$Year <- as.numeric(as.character(state_resistance_mdr$Year))
years_range <- seq(min(state_resistance_mdr$Year), max(state_resistance_mdr$Year), length.out = 100)
# Create prediction data
new_data <- expand.grid(Year = years_range, GEOID = original_geo_levels)
new_data$Year <- as.numeric(as.character(new_data$Year))  # Ensure Year is nume

# Compute predictions and derivatives
preds <- predict(model_output_mdr, newdata = new_data, type = "terms", se.fit = TRUE, deriv = 1)

# Assuming you have already generated `preds` as before
# Now, include the interaction term in the derivative calculation
new_data$Year_deriv <- preds$fit[, "s(Year)"] + preds$fit[, "te(Year,GEOID)"]

# Calculate confidence intervals including both main and interaction effects
new_data$SE <- sqrt(preds$se.fit[, "s(Year)"]^2 + preds$se.fit[, "te(Year,GEOID)"]^2)

# Using a normal approximation for the confidence interval
alpha <- 0.05
z_score <- qnorm(1 - alpha / 2)
new_data$lower_ci <- new_data$Year_deriv - z_score * new_data$SE
new_data$upper_ci <- new_data$Year_deriv + z_score * new_data$SE

# Convert Spatial*DataFrame to a regular data frame for easier manipulation
states_data_df <- as.data.frame(states_spatial_filtered)

# Ensure GEOID is factor in both data frames for correct merging
new_data$GEOID <- as.factor(new_data$GEOID)
states_data_df$GEOID <- as.factor(states_data_df$GEOID)

# Merge to append the NAME corresponding to each GEOID
new_data <- merge(new_data, states_data_df[, c("GEOID", "NAME")], by = "GEOID", all.x = TRUE)

# Calculate zero crossing: Check if zero is inside or outside the interval
new_data <- new_data %>%
  mutate(
    zero_in_prev = ifelse(lag(lower_ci) <= 0 & lag(upper_ci) >= 0, 1, 0),
    zero_in_curr = ifelse(lower_ci <= 0 & upper_ci >= 0, 1, 0),
    sign_change = zero_in_prev != zero_in_curr  # Change in zero inclusion status
  )

# Filter crossings that happen after the year 2004
crossing_points <- new_data %>%
  filter(sign_change, Year > 2004) %>%
  dplyr::select(NAME, Year) %>%
  distinct(NAME, Year)  # Ensure unique crossing points

# Define the Lancet style theme
lancet_theme <- theme_minimal() +
  theme(text = element_text(family = "sans", color = "black"),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7, face = "bold"),
        legend.position = "none",  # Remove legend
        strip.text.x = element_text(size = 7))

# Plotting with NAME as labels and Lancet style, including multiple crossing points
p <- ggplot(new_data, aes(x = Year, y = Year_deriv, group = NAME, color = NAME)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = NAME), alpha = 0.2) +
  geom_vline(data = crossing_points, aes(xintercept = Year), color = "red", linetype = "dashed", size = 1) +
  facet_wrap(~NAME, scales = "free_y") +
  lancet_theme +
  labs(title = "",
       x = "Year",
       y = "Rate of change (first derivative from spatio-temporal GAM model)")
#First Derivative of Year by State in Carbapenem Resistance Model
ggsave(filename = "p_state_resistance_mdr.tiff", plot = p, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
######


#-------------------------------------------#
###SUBGROUP ANALYSES, US states:
#-------------------------------------------#




#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
# Analyses for EUROPE
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#country_resistance_carbap_eu country_resistance_cephalos_eu country_resistance_firstline_eu country_resistance_mdr_eu

#------------------------------------------------------------------------------------------------------#
# Analyses for Europe GAM SPATIOTEMPORAL 
#------------------------------------------------------------------------------------------------------#
######
# Load Spatial Data and Resistance Data
k=6
fit_gam_model_eu <- function(state_resistance_carbap) {
  europe_shapefile <- ne_countries(continent = "Europe", returnclass = "sf")
  europe_shapefile <- st_make_valid(europe_shapefile)
  # Rename and modify sovereignt to NAME
  europe_shapefile <- mutate(europe_shapefile, NAME = sovereignt)
  europe_shapefile <- mutate(europe_shapefile,
                             NAME = case_when(
                               NAME == "Slovakia" ~ "Slovak Republic",
                               NAME == "Republic of Serbia" ~ "Serbia",
                               NAME == "Czechia" ~ "Czech Republic",
                               TRUE ~ NAME
                             ))
  
  # Ensure all geometries are valid
  states_shapefile<-europe_shapefile
  state_resistance_carbap$NAME<- state_resistance_carbap$Country
  state_resistance_carbap <- state_resistance_carbap %>%
    filter(NAME != "Turkey") 
    # Merge spatial data with resistance data by 'NAME'
  states_merged <- merge(states_shapefile, state_resistance_carbap, by = "NAME", all.x = TRUE)
  
  # Create unique geometries for neighborhood structure
  states_spatial_unique <- states_merged[!duplicated(states_merged$NAME), ]
  #states_spatial_unique$GEOID <- as.numeric(factor(states_spatial_unique$NAME))
  
  # Merge to append the GEOID corresponding to each NAME
  #states_spatial <- st_join(states_merged, states_spatial_unique[, c("NAME", "GEOID")])
  #states_spatial <- states_spatial %>%
  #  rename(NAME = NAME.x) %>%
  #  dplyr::select(-NAME.y)
  # Remove countries not included in countries_list_eu_includ
  states_spatial_filtered <- states_merged %>%
    dplyr::filter(NAME %in% countries_list_eu_includ)
  
  states_spatial_filtered <- states_spatial_filtered %>%
    mutate(
      GEOID = as.numeric(factor(NAME, levels = unique(NAME)))
    )
  
  if (!inherits(states_spatial_filtered, "Spatial")) {
    states_spatial_filtered <- as(states_spatial_filtered, "Spatial")
  }
  
  
  # 8. Remove duplicates based on GEOID to create unique state geometries for neighborhood structure
  states_spatial_unique <- states_spatial_filtered[!duplicated(states_spatial_filtered$GEOID), ]
  # 9. Recreate the neighborhood structure using unique state geometries
  nb <- poly2nb(states_spatial_unique, row.names = states_spatial_unique$GEOID)
  # 10. Assign region names to the neighborhood list
  names(nb) <- states_spatial_unique$GEOID
  
  # Clean and ensure data alignment with neighborhood structure
  states_spatial_filtered@data$GEOID <- factor(states_spatial_filtered@data$GEOID, levels = attr(nb, "region.id"))
  states_spatial_filtered@data$Resistance <- as.numeric(states_spatial_filtered@data$Resistance)
  states_spatial_filtered@data$AMR_Positive <- as.numeric(states_spatial_filtered@data$AMR_Positive)
  states_spatial_filtered@data$Total_Isolates <- as.numeric(states_spatial_filtered@data$Total_Isolates)
  states_spatial_filtered@data$Year <- as.numeric(as.character(states_spatial_filtered@data$Year))
  #states_spatial_filtered@data <- states_spatial_filtered@data[complete.cases(states_spatial_filtered@data), ]
  
  # Set control parameters for GAM
  ctrl <- gam.control(nthreads = 6)  # Set parallel threads for faster computation
  
  # Fit the MRF model with spatio-temporal smoothing
  model_gamx1 <- gam(Resistance ~ 
                       s(Year, m=3, k=k, bs = "tp") +
                       s(GEOID, bs = 'mrf', xt = list(nb = nb)) +
                       te(Year, GEOID, bs=c("tp", "mrf"), m=c(3, NA), xt=list(Year=NULL, GEOID=list(nb=nb))),
                     data = states_spatial_filtered@data,
                     family = binomial, 
                     select = TRUE,
                     method = "REML", 
                     weights = Total_Isolates)
  return(list(fr=model_gamx1, states_spatial_filtered=states_spatial_filtered))
}
fit_gam_model_mdr_eu <- function(state_resistance_carbap) {
  # Load Spatial Data and Resistance Data
  europe_shapefile <- ne_countries(continent = "Europe", returnclass = "sf")
  europe_shapefile <- st_make_valid(europe_shapefile)
  # Rename and modify sovereignt to NAME
  europe_shapefile <- mutate(europe_shapefile, NAME = sovereignt)
  europe_shapefile <- mutate(europe_shapefile,
                             NAME = case_when(
                               NAME == "Slovakia" ~ "Slovak Republic",
                               NAME == "Republic of Serbia" ~ "Serbia",
                               NAME == "Czechia" ~ "Czech Republic",
                               TRUE ~ NAME
                             ))
  
  # Ensure all geometries are valid
  states_shapefile<-europe_shapefile
  state_resistance_carbap$NAME<- state_resistance_carbap$Country
  state_resistance_carbap <- state_resistance_carbap %>%
    filter(NAME != "Turkey") 
  # Merge spatial data with resistance data by 'NAME'
  states_merged <- merge(states_shapefile, state_resistance_carbap, by = "NAME", all.x = TRUE)
  
  # Create unique geometries for neighborhood structure
  states_spatial_unique <- states_merged[!duplicated(states_merged$NAME), ]
  #states_spatial_unique$GEOID <- as.numeric(factor(states_spatial_unique$NAME))
  
  # Merge to append the GEOID corresponding to each NAME
  #states_spatial <- st_join(states_merged, states_spatial_unique[, c("NAME", "GEOID")])
  #states_spatial <- states_spatial %>%
  #  rename(NAME = NAME.x) %>%
  #  dplyr::select(-NAME.y)
  # Remove countries not included in countries_list_eu_includ
  countries_list_eu_includ <- c(
    "Austria", "Belgium", "Bulgaria", "Croatia", "Czech Republic",
    "Denmark","Finland", "France", "Germany",
    "Greece", "Hungary", "Ireland", "Italy", "Latvia",
    "Lithuania", "Netherlands", "Poland", "Portugal",
    "Romania", "Russia", "Slovak Republic", "Slovenia",
    "Spain", "Sweden", "Switzerland", "Ukraine",
    "United Kingdom"
  )
  countries_list_eu_includ <- c(
    "Austria", "Belgium", "Croatia", "Czech Republic",
    "Denmark", "Finland", "France", "Germany",
    "Greece", "Hungary", "Ireland", "Italy", "Latvia",
    "Lithuania", "Netherlands", "Poland", "Portugal",
    "Romania",  "Slovak Republic", "Slovenia",
    "Spain", "Sweden", "Switzerland",
    "United Kingdom"
  )
  states_spatial_filtered <- states_merged %>%
    dplyr::filter(NAME %in% countries_list_eu_includ)
  
  states_spatial_filtered <- states_spatial_filtered %>%
    mutate(
      GEOID = as.numeric(factor(NAME, levels = unique(NAME)))
    )
  
  if (!inherits(states_spatial_filtered, "Spatial")) {
    states_spatial_filtered <- as(states_spatial_filtered, "Spatial")
  }
  
  
  # 8. Remove duplicates based on GEOID to create unique state geometries for neighborhood structure
  states_spatial_unique <- states_spatial_filtered[!duplicated(states_spatial_filtered$GEOID), ]
  # 9. Recreate the neighborhood structure using unique state geometries
  nb <- poly2nb(states_spatial_unique, row.names = states_spatial_unique$GEOID)
  # 10. Assign region names to the neighborhood list
  names(nb) <- states_spatial_unique$GEOID

  # 11. Clean and ensure data alignment with neighborhood structure
  # Convert 'GEOID' to factor and ensure levels match the neighborhood structure
  states_spatial_filtered@data$GEOID <- factor(states_spatial_filtered@data$GEOID, levels = attr(nb, "region.id"))
  # Convert 'Resistance' to numeric
  states_spatial_filtered@data$MDR <- as.numeric(states_spatial_filtered@data$MDR)
  states_spatial_filtered@data$MDR_Positive <- as.numeric(states_spatial_filtered@data$MDR_Positive)
  states_spatial_filtered@data$Total_Isolates <- as.numeric(states_spatial_filtered@data$Total_Isolates)
  # Convert 'Year' to numeric and handle missing data
  states_spatial_filtered@data$Year <- as.numeric(as.character(states_spatial_filtered@data$Year))
  # 12. Recreate neighborhood structure using the cleaned and filtered data
  states_spatial_unique <- states_spatial_filtered[!duplicated(states_spatial_filtered$GEOID), ]
  nb <- poly2nb(states_spatial_unique, row.names = states_spatial_unique$GEOID)
  names(nb) <- states_spatial_unique$GEOID
  
  # 13. Set control parameters for GAM
  ctrl <- gam.control(nthreads = 6)  # Set parallel threads for faster computation
  # 14. Fit the MRF model with spatio-temporal smoothing
  model_gamx1 <- gam(MDR ~ 
                       s(Year, m=3, k=k, bs = "tp") +   # Smooth term for Year
                       s(GEOID, bs = 'mrf', xt = list(nb = nb)) +  # Smooth term for GEOID
                       te(Year, GEOID, bs=c("tp", "mrf"), m=c(3, NA), xt=list(Year=NULL, GEOID=list(nb=nb))), # Interaction term
                     data = states_spatial_filtered@data,
                     family = binomial, 
                     select = TRUE,
                     method = "REML", 
                     weights = Total_Isolates)
  return(list(sf=model_gamx1,  states_spatial_filtered= states_spatial_filtered))
}

#OUTPUT OF THE RESULTS SPATIO-TEMPORAL MODELS:
model_output_carbap_eu <- fit_gam_model_eu(country_resistance_carbap_eu)
model_output_carbap_eu <- model_output_carbap_eu$fr
model_output_cephalos_eu <- fit_gam_model_eu(country_resistance_cephalos_eu)
model_output_cephalos_eu<- model_output_cephalos_eu$fr
model_output_firstline_eu <- fit_gam_model_eu(country_resistance_firstline_eu)
model_output_firstline_eu<- model_output_firstline_eu$fr
model_output_mdr_eu <-fit_gam_model_mdr_eu(country_resistance_mdr_eu)
model_output_mdr_eu <- model_output_mdr_eu$sf

# Tidying model outputs
extract_model_summary_eu <- function(model) {
  # Extracting smooth terms
  smooth_summary <- tidy(model, parametric = FALSE)
  
  # Extracting parametric coefficients
  parametric_summary <- tidy(model, exponentiate = FALSE, parametric = TRUE)
  
  # Combine both summaries
  combined_summary <- bind_rows(parametric_summary, smooth_summary)
  
  return(combined_summary)
}

# Apply the function to each model
summary_carbap_eu <- extract_model_summary_eu(model_output_carbap_eu)
summary_cephalos_eu <- extract_model_summary_eu(model_output_cephalos_eu)
summary_firstline_eu <- extract_model_summary_eu(model_output_firstline_eu)
summary_mdr_eu <- extract_model_summary_eu(model_output_mdr_eu)

# Adding model identifiers
summary_carbap_eu$model <- "Carbapenem Resistance"
summary_cephalos_eu$model <- "Cephalosporin Resistance"
summary_firstline_eu$model <- "First-line Antibiotic Resistance"
summary_mdr_eu$model <- "MDR Resistance"
# Combining all summaries into one dataframe
combined_results_eu <- bind_rows(summary_carbap_eu, summary_cephalos_eu, summary_firstline_eu, summary_mdr_eu)
# Optionally, select and rename columns for clarity
final_results_eu <- combined_results_eu %>%
  dplyr::select(Model = model, Term = term, Estimate = estimate, EDF = edf, Ref.DF = ref.df, Std.Error = std.error,
                Statistic = statistic, `P.Value` = p.value)

######


#-------------------------------------------------------------------------#
###FIRST/SECOND-derivative & GROWTH rate/DOUBLING/HALVING GRAPHS, Europe:
#-------------------------------------------------------------------------#
#FUNCTIONS for final plots///
# - - - - - -#
gam_predictions <- function(gam_model, newdata) {
  prediction <- predict(gam_model, newdata = newdata, type = "response", se = TRUE)
  prediction_df <- data.frame(
    pred = prediction$fit * 100,  # convert to percentage
    pred_lower = (prediction$fit - (2 * prediction$se.fit)) * 100,
    pred_upper = (prediction$fit + (2 * prediction$se.fit)) * 100,
    se = prediction$se.fit
  )
  prediction_df <- cbind(newdata, prediction_df)
  return(prediction_df)
}
derivatives_mh <- function(gam_model, newdata, n = 200, h1 = 1e-7, h2 = 1e-7, startpoint = 0, type = c("forward", "backward", "central")) {
  if (!type %in% c("forward", "backward", "central")) { 
    stop("Type must be either forward, backward, or central")
  }
  
  linv <- gam_model$family$linkinv 
  set.seed(2929)
  br <- gam.mh(gam_model, thin = 2, ns = 2000, rw.scale = 0.4)$bs 
  
  # Split the data by group and calculate derivatives per group
  results <- newdata %>%
    group_by(GEOID) %>%
    do({
      data <- .
      # Function to calculate first and second derivatives
      calc_derivatives <- function(data, h) {
        X0 <- predict(gam_model, type = "lpmatrix", newdata = data)
        data_subh <- mutate(data, Year = Year - h)
        data_addh <- mutate(data, Year = Year + h)
        X_subh <- predict(gam_model, type = "lpmatrix", newdata = data_subh)
        X_addh <- predict(gam_model, type = "lpmatrix", newdata = data_addh)
        ff0 <- linv(X0 %*% t(br))
        ff_subh <- linv(X_subh %*% t(br))
        ff_addh <- linv(X_addh %*% t(br))
        return(list(ff0 = ff0, ff_subh = ff_subh, ff_addh = ff_addh))
      }
      
      # Calculate derivatives
      derivatives1 <- calc_derivatives(data, h1)
      firstdiv <- switch(type,
                         forward = (derivatives1$ff_addh - derivatives1$ff0) / h1, 
                         backward = (derivatives1$ff0 - derivatives1$ff_subh) / h1, 
                         central = (derivatives1$ff_addh - derivatives1$ff_subh) / (2 * h1))
      
      derivatives2 <- calc_derivatives(data, h2)
      seconddiv <- switch(type,
                          forward = (derivatives2$ff_addh - 2 * derivatives2$ff0 + derivatives2$ff_subh) / (h2^2),
                          backward = (derivatives2$ff_subh - 2 * derivatives2$ff0 + derivatives2$ff_addh) / (h2^2),
                          central = (derivatives2$ff_addh - 2 * derivatives2$ff0 + derivatives2$ff_subh) / (h2^2))
      
      # Prepare output
      ffquant1 <- apply(firstdiv, 1, quantile, probs = c(0.025, 0.5, 0.975))
      sdquant2 <- apply(seconddiv, 1, quantile, probs = c(0.025, 0.5, 0.975))
      
      data %>%
        mutate(
          first_derivative = ffquant1[2, ],
          first_lower = ffquant1[1, ],
          first_upper = ffquant1[3, ],
          second_derivative = sdquant2[2, ], 
          second_lower = sdquant2[1, ], 
          second_upper = sdquant2[3, ],
          second_derivative_signif = ifelse(second_lower <= 0 & second_upper >= 0, 0, 1),
          derivative_breakpoint = ifelse(second_derivative_signif == 1 & lag(second_derivative_signif, default = second_derivative_signif[1]) == 0, 1, 0),
          derivative_breakpoint = ifelse(row_number() < startpoint, 0, derivative_breakpoint),
          first_derivative_sign_change = (lag(first_upper < 0, default = first_upper[1] < 0) & first_lower > 0) | (lag(first_lower > 0, default = first_lower[1] > 0) & first_upper < 0),
          direction = case_when(
            derivative_breakpoint == 1 & second_lower > 0 & lag(second_lower, default = second_lower[1]) < 0 ~ "increasing",
            derivative_breakpoint == 1 & second_upper < 0 & lag(second_upper, default = second_upper[1]) > 0 ~ "decreasing",
            TRUE ~ NA_character_
          )
        )
    })
  
  return(results)
}
simulate_growth_rate_ci <- function(first_derivative, second_derivative, n_sim = 1000) {
  set.seed(123)  # for reproducibility
  
  simulations <- replicate(n_sim, {
    # Handle possible NA or zero values
    first_mean <- ifelse(is.na(first_derivative) | first_derivative == 0, 0.001, first_derivative)
    second_mean <- ifelse(is.na(second_derivative), 0, second_derivative)
    
    # Ensure non-negative SD
    first_sd <- max(0.001, abs(first_mean * 0.1))
    second_sd <- max(0.001, abs(second_mean * 0.1))
    
    first_sim <- rnorm(1, mean = first_mean, sd = first_sd)
    second_sim <- rnorm(1, mean = second_mean, sd = second_sd)
    
    # Avoid division by zero
    if (first_sim == 0) first_sim <- 0.001
    first_sim + (second_sim / first_sim)
  })
  
  # Calculate the quantiles for the simulated growth rates
  ci_lower <- quantile(simulations, probs = 0.25)
  ci_upper <- quantile(simulations, probs = 0.75)
  ci_median <- median(simulations)
  
  # Return as a list to fit within do()
  list(median = ci_median, lower_ci = ci_lower, upper_ci = ci_upper)
}
derivatives_mh2 <- function(gam_model, newdata, startpoint = 0) {
  set.seed(2929)
  br <- gam.mh(gam_model, thin = 2, ns = 2000, rw.scale = 0.4) $bs
  linv <- gam_model$family$linkinv 
  
  # Calculate derivatives and apply GAM MH method
  results <- newdata %>%
    group_by(GEOID) %>%
    mutate(
      # Calculate first and second derivatives
      First_Derivative = c(NA, diff(pred)),
      Second_Derivative = c(NA, diff(First_Derivative)),
      First_Derivative = ifelse(is.na(First_Derivative), 0, First_Derivative),
      Second_Derivative = ifelse(is.na(Second_Derivative), 0, Second_Derivative),
      first_lower= c(NA, diff(pred_lower)),
      first_upper= c(NA, diff(pred_upper)),
      second_upper = c(NA, diff(first_upper)),
      second_lower = c(NA, diff(first_lower)),
      first_lower = ifelse(is.na(first_lower), 0, first_lower),
      second_lower = ifelse(is.na(second_lower), 0, second_lower),
      first_upper = ifelse(is.na(first_upper), 0, first_upper),
      second_upper = ifelse(is.na(second_upper), 0, second_upper)
    ) %>%
    ungroup() %>%
    mutate(
      # Calculate growth rate safely
      growth_rate = ifelse(First_Derivative != 0, 
                           First_Derivative + (Second_Derivative / First_Derivative), 
                           0),
      growth_rate = ifelse(is.nan(growth_rate), 0, growth_rate),
      # Statistical significance and breakpoints
      second_derivative_signif = ifelse(second_lower <= 0 & second_upper >= 0, 0, 1),
      derivative_breakpoint = ifelse(second_derivative_signif == 1 & lag(second_derivative_signif, default = second_derivative_signif[1]) == 0, 1, 0),
      derivative_breakpoint = ifelse(row_number() < startpoint, 0, derivative_breakpoint),
      first_derivative_sign_change = (lag(first_upper < 0, default = first_upper[1] < 0) & first_lower > 0) | (lag(first_lower > 0, default = first_lower[1] > 0) & first_upper < 0),
      direction = case_when(
        derivative_breakpoint == 1 & second_lower > 0 & lag(second_lower, default = second_lower[1]) < 0 ~ "increasing",
        derivative_breakpoint == 1 & second_upper < 0 & lag(second_upper, default = second_upper[1]) > 0 ~ "decreasing",
        TRUE ~ NA_character_
      )
    )
  
  return(results)
}

# - - - - - -#

#------------------------------------------------------------------------------#
#NEW, CARBAPENEM RESISTANCE GRAPH: ######
country_resistance_carbap_eu <- country_resistance_carbap_eu%>%
  filter(Country != "Turkey") 
model_output_carbap_eu <- fit_gam_model_eu(country_resistance_carbap_eu)
final_dataset<- model.frame(model_output_carbap_eu$fr)
original_geo_levels <- levels(final_dataset$GEOID)

country_resistance_carbap_eu$Year <- as.numeric(as.character(country_resistance_carbap_eu$Year))
years_range <- seq(min(country_resistance_carbap_eu$Year), max(country_resistance_carbap_eu$Year), length.out = 100)
# Create prediction data
new_data <- expand.grid(Year = years_range, GEOID = original_geo_levels)
new_data$Year <- as.numeric(as.character(new_data$Year))  # Ensure Year is nume
states_spatial_fitered<- model_output_carbap_eu$states_spatial_filtered
states_data_df <- as.data.frame(states_spatial_fitered@data)

original_geo_levels <- unique(new_data$GEOID)
GEOID_nameC <- states_data_df %>%
  dplyr::select(GEOID, NAME)
unique_GEOID_data <- GEOID_nameC %>%
  distinct(GEOID, .keep_all = TRUE)
# Expand pred_data to include every combination of Year and GEOID
pred_data <- expand.grid(Year = seq(min(new_data$Year), max(new_data$Year), by = 1),
                         GEOID = original_geo_levels)

# Predict from GAM using the created function
predictions <- gam_predictions(model_output_carbap_eu$fr, newdata = pred_data)
predictions <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE)
predictions_carb_EU <- gam_predictions(model_output_carbap_eu$fr, newdata = pred_data)
predictions_carb_EU <- merge(predictions_carb_EU, unique_GEOID_data, by = "GEOID", all.x = TRUE)

resultsGrowt <- derivatives_mh2(model_output_carbap_eu$fr, predictions_carb_EU)


# Usage of the function
derivatives_data <- derivatives_mh(model_output_carbap_eu$fr, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001,startpoint = 0)

derivatives_data$growth_rate <- derivatives_data$first_derivative + (derivatives_data$second_derivative / derivatives_data$first_derivative)

merged_data <- merge(derivatives_data, unique_GEOID_data, by = "GEOID", all.x = TRUE)
predictions <- predictions %>%
  group_by(NAME) %>%
  mutate(
    # Using lag to shift the pred values down
    pred_lag = lag(pred, default = NA),
    pred_lagu = lag(pred_upper, default = NA),
    pred_lagl = lag(pred_lower, default = NA),
    
    # Calculate the growth rate
    growth_rate2 = if_else(
      is.na(pred_lag),
      NA_real_,  # Ensures that the first value where lag is NA gets an NA in growth_rate2
      100 * (pred - pred_lag) / pred_lag),  # Calculate percentage change

    growth_rate2_up = if_else(
      is.na(pred_lagu),
      NA_real_,  # Ensures that the first value where lag is NA gets an NA in growth_rate2
      100 * (pred_upper - pred_lagu) / pred_lagu),  # Calculate percentage change

    growth_rate2_lo = if_else(
      is.na(pred_lagl),
      NA_real_,  # Ensures that the first value where lag is NA gets an NA in growth_rate2
      100 * (pred_lower - pred_lagl) / pred_lagl)  # Calculate percentage change
  ) %>%
  dplyr::select(-pred_lag)  # Remove the temporary lag column

Carb_predictions_grat<-predictions 
Carb_changep_EU<- merged_data

# Plotting
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)
# Calculate the maximum and minimum values for setting y-axis limits
y_max <- max(merged_data$first_upper, na.rm = TRUE)
y_min <- min(merged_data$first_lower, na.rm = TRUE)
# If needed, you can add some padding to ensure all data points are within the view
padding <- (y_max - y_min) * 0.05  # 5% padding
y_max <- y_max + padding
y_min <- y_min - padding
# Plotting with adjustments
# Plotting with universal y-axis limits
p <- ggplot(data = merged_data, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022,labels = as.character(2004:2022))+# Line at y=0
  scale_color_manual(values = palette) +  # Assuming palette is predefined
  scale_fill_manual(values = palette) +
  labs(title = "",
       x = "Year",
       y = "First derivative (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size=7)
  ) +
  facet_wrap(~NAME, scales = "free_y")   # Apply the same y-axis limits to all facets

# Display the plot
ggsave(filename = "first_derivat_eu_carb.tiff", plot = p, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

#-----------------------------#
#PREDICTIONS: per country Europe
#-----------------------------#
max_y <- max(predictions$pred_upper, na.rm = TRUE)
# Generate the plot
ppred <- ggplot(data = predictions, aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +  # Use fixed scales for y-axis across all facets
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(0, max_y + 0), breaks = seq(0, max_y + 0, by = 5)) +
  labs(title = "",
       x = "Year",
       y = "Predicted Carbapenem-Resistance (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
# Display the plot
print(ppred)
ggsave(filename = "predictions_breakpoint_carb.tiff", plot = ppred, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


netherlands_predictions <- predictions %>%
  filter(NAME == "Netherlands")
merged_datanetherlands <- merged_data %>%
  filter(NAME == "Netherlands")
netherlands_plot<-ggplot(data = netherlands_predictions, aes(x = Year, y = pred)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper), alpha = 0.2) +
  geom_vline(data = filter(merged_datanetherlands , first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_datanetherlands, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +  # Use fixed scales for y-axis across all facets
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(0, 25 + 0), breaks = seq(0, 25 + 0, by = 5)) +
  labs(title = "",
       x = "Year",
       y = "Predicted carbapenem-resistance (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
netherlands_plot

#manually calcualted derivatives.
netherlands_dataxx <- data.frame(
  Year = 2004:2022,
  Resistance = c(4.9775106, 2.9577511, 2.2512514, 2.1623578, 2.4997834,
                 3.2371703, 4.3422959, 5.6629015, 6.9452633, 8.0076611,
                 8.8846535, 9.7900113, 10.9687793, 12.5527863, 14.4403619,
                 16.2290809, 17.3129993, 17.1747504, 15.6667596)
)
# Calculate the first derivative (Rate of Change)
netherlands_dataxx$First_Derivative <- c(NA, diff(netherlands_dataxx$Resistance))
netherlands_dataxx$Second_Derivative <- c(NA, diff(netherlands_dataxx$First_Derivative))
netherlands_dataxx$First_Derivative[is.na(netherlands_dataxx$First_Derivative)] <- 0
netherlands_dataxx$Second_Derivative[is.na(netherlands_dataxx$Second_Derivative)] <- 0
netherlands_dataxx$growth_rate <- netherlands_dataxx$First_Derivative + (netherlands_dataxx$Second_Derivative/netherlands_dataxx$First_Derivative)
netherlands_dataxx$growth_rate[is.nan(netherlands_dataxx$growth_rate)] <- 0
# Display the data frame with derivatives
print(netherlands_dataxx)



#-----------------------------#
#PREDICTIONS: All together
#-----------------------------#
# Step 1: Prepare the data for labels
label_data <- predictions %>%
  group_by(NAME) %>%
  filter(Year == max(Year)) %>%
  ungroup()

final_palette <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(24)
# Step 2: Generate the plot
ppred_alle <- ggplot(data = predictions, aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_label_repel(data = label_data, aes(label = NAME, y = pred), 
                   point.padding = 0.2, nudge_x = 1, direction = 'y', 
                   size = 3.5, color = "black", fontface = "bold",
                   box.padding = 0.35, segment.color = "grey50",
                   fill = "white") +  # White background for labels
  scale_color_manual(values = final_palette) +
  scale_fill_manual(values = final_palette) +
  scale_x_continuous(breaks = 2004:2022) +
  scale_y_continuous(limits = c(NA, max(predictions$pred_upper) +4)) +
  scale_y_continuous(breaks = seq(0, 40, by = 5))+
  labs(title = "",
       x = "Year",
       y = "Predicted carbapenem-resistance (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 360, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size=13),
    axis.title.y = element_text(size=13)
  )



# Calculate max and min values
extreme_values <- predictions %>%
  group_by(NAME) %>%
  summarize(max_pred = max(pred, na.rm = TRUE), min_pred = min(pred, na.rm = TRUE)) %>%
  ungroup()
# Ranking for maximum and minimum predictions
extreme_values <- extreme_values %>%
  arrange(desc(max_pred)) %>%
  mutate(rank_max = row_number()) %>%
  arrange(min_pred) %>%
  mutate(rank_min = row_number())
# Filter to get only the top 2 and bottom 2
extreme_values <- extreme_values %>%
  filter(rank_max <= 3 | rank_min <= 3)
# Prepare highlighted names list
highlighted_names <- unique(c(extreme_values$NAME[extreme_values$rank_max <= 3], 
                              extreme_values$NAME[extreme_values$rank_min <= 3]))
predictions <- predictions %>%
  mutate(
    color = ifelse(NAME %in% highlighted_names, NAME, NA),
    pred = ifelse(is.na(pred), 0, pred)  # Replace NA predictions with 0 or another suitable default value
  )
# Check range of 'pred' to adjust y-axis limits
range(predictions$pred, na.rm = TRUE)
# Define the color palette
unique_names <- unique(predictions$color, na.rm = TRUE)
num_colors <- length(unique_names)

# Choose a palette, e.g., 'Set1' which is good for categorical data
# Ensure that there are enough colors, if not repeat the palette
palette_colors <- brewer.pal(min(num_colors, 8), "Set3")
if (num_colors > 8) {
  palette_colors <- rep(palette_colors, length.out = num_colors)
}

final_palette <- setNames(palette_colors, unique_names)
# Create the plot
ppred_alle2x <- ggplot(data = predictions, aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = color), size = 1, na.rm = TRUE) +  # Ensure NA values in color are ignored
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  scale_color_manual(values = final_palette, na.translate = FALSE) +
  scale_fill_manual(values = final_palette, guide = "none") +  # Make sure there is a plus sign at the end
  scale_x_continuous(
    breaks = seq(from = 2004, to = 2022, by = 4),
    limits = c(2004, 2022)  # Closing parenthesis for limits
  ) +  # Closing parenthesis was missing after 2022
  scale_y_continuous(
    limits = c(min(predictions$pred_lower, na.rm = TRUE) - 0, 
               max(predictions$pred_upper, na.rm = TRUE) + 0), 
    breaks = seq(0, 35, by = 5)
  ) +
  labs(title = "C. Predicted CRE (%)",
       x = "Year",
       y = "Predicted carbapenem-resistance (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 360, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13)
  )

# Print the plot
print(ppred_alle2x)


ggsave(filename = "pred_allEurop_toget_carb.tiff", plot = ppred_alle, device = "tiff", path = base_pathOut,
       width = 16, height = 10, dpi = 500, units = "in")

#-----------------------------#
# Plotting SECOND DERIVATIVE:
#-----------------------------#
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)
# Calculate the maximum and minimum values for setting y-axis limits
y_max <- max(merged_data$second_upper, na.rm = TRUE)
y_min <- min(merged_data$second_lower, na.rm = TRUE)
# If needed, you can add some padding to ensure all data points are within the view
padding <- (y_max - y_min) * 0.05  # 5% padding
y_max <- y_max + padding
y_min <- y_min - padding
# Plotting with adjustments
# Plotting with universal y-axis limits
p2<- ggplot(data = merged_data, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Line at y=0
  scale_color_manual(values = palette) +  # Assuming palette is predefined
  scale_fill_manual(values = palette) +
  labs(title = "",
       x = "Year",
       y = "Second derivative (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y") #+
#scale_y_continuous(limits = c(y_min, y_max))  # Apply the same y-axis limits to all facets

# Display the plot
ggsave(filename = "second_derivat_eu_carb.tiff", plot = p2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")



p2netherl<- ggplot(data = merged_datanetherlands, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_datanetherlands, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_datanetherlands, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid", color = "red") +  # Line at y=0
  scale_color_manual(values = palette) +  # Assuming palette is predefined
  scale_fill_manual(values = palette) +
  labs(title = "",
       x = "Year",
       y = "Second derivative (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")



#-----------------------------#
#GROWTH RATE: 
#-----------------------------#
# Assuming derivatives_data is already prepared and using mutate for efficiency
growth_rate_ci <- resultsGrowt  %>% #merged_data
  rowwise() %>%
  mutate(
    growth_rate_results = list(simulate_growth_rate_ci(First_Derivative, Second_Derivative, n_sim = 1000))
  ) %>%
  unnest_wider(c(growth_rate_results)) %>%
  ungroup()

#growth_rate_ci$ci_lower <- growth_rate_ci$growth_rate +growth_rate_ci$ci_lower
#growth_rate_ci$ci_upper <- growth_rate_ci$growth_rate +growth_rate_ci$ci_upper
growth_rate_ci <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE)

growth_rate_ci$doubling_times <-  log(2) / growth_rate_ci$growth_rate2 # Calculate doubling time
growth_rate_ci$halving_times <- log(0.5) / growth_rate_ci$growth_rate2   # Calculate halving time
growth_rate_ci$NAME<- growth_rate_ci$NAME.x

growth_rate_ci <- growth_rate_ci %>%
  mutate(
    growth_rate2 = if_else(Year == 2004, NA, growth_rate2),
    growth_rate2_lo = if_else(Year == 2004, NA, growth_rate2_lo),
    growth_rate2_up = if_else(Year == 2004, NA, growth_rate2_up),
    doubling_times = if_else(Year == 2004, NA, doubling_times),
    halving_times = if_else(Year == 2004, NA, halving_times)
  )


p3<- ggplot(data = growth_rate_ci, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  #geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_ribbon(aes(ymin = growth_rate2_lo, ymax = growth_rate2_up), alpha = 0.2) +  # Shaded area for CIs
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Line at y=0
  scale_color_manual(values = palette) +  # Assuming palette is predefined
  scale_fill_manual(values = palette) +
  labs(title = "",
       x = "Year",
       y = "Growth rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    scale_x_continuous(breaks = 2004:2022,labels = as.character(2004:2022)),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  facet_wrap(~NAME, scales = "free_y")  # Apply the same y-axis limits to all facets

# Display the plot
ggsave(filename = "growth_ratio_carbEU.tiff", plot = p3, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")




netherlands_growth_rate_ci <- growth_rate_ci %>%
  filter(NAME == "Netherlands")
merged_datanetherlands <- merged_data %>%
  filter(NAME == "Netherlands")

p3nethelands<- ggplot(data = netherlands_growth_rate_ci, aes(x = Year, y = growth_rate2)) +
  geom_line(aes(color = NAME), size = 1) +
  #geom_ribbon(aes(ymin = growth_rate2l, ymax = growth_rate2u, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_datanetherlands, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_datanetherlands, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_ribbon(aes(ymin = growth_rate2_lo, ymax = growth_rate2_up), alpha = 0.2) +  # Shaded area for CIs
  geom_hline(yintercept = 0, linetype = "solid", color = "red") +  # Line at y=0
  scale_x_continuous(breaks = 2005:2022, labels = as.character(2005:2022)) +
  scale_color_manual(values = palette) +  # Assuming palette is predefined
  scale_fill_manual(values = palette) +
  labs(title = "",
       x = "Year",
       y = "Growth rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)
  ) +
  facet_wrap(~NAME, scales = "free_y")  # Apply the same y-axis limits to all facets



# Plot Doubling Times
p_doubling <- ggplot(data = growth_rate_ci, aes(x = Year, y = doubling_times, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_point(aes(color = NAME), size = 2) +
  labs(title = "",
       x = "Year",
       y = "Doubling Time (in units of time)") +
  theme_minimal() +
  theme(text = element_text(size = 12, family = "Times New Roman"),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  facet_wrap(~NAME, scales = "free_y")

# Plot Halving Times
p_halving <- ggplot(data = growth_rate_ci, aes(x = Year, y = halving_times, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_point(aes(color = NAME), size = 2) +
  labs(title = "",
       x = "Year",
       y = "Halving Time (in units of time)") +
  theme_minimal() +
  theme(text = element_text(size = 12, family = "Times New Roman"),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  facet_wrap(~NAME, scales = "free_y")

# Print the plots
print(p_doubling)
print(p_halving)


ggsave(filename = "p_doublingEU_carb.tiff", plot = p_doubling, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "p_halvingEU_carb.tiff", plot = p_halving, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")









######
#------------------------------------------------------------------------------#
#NEW, CEPHALOSPORIN RESISTANCE GRAPH: ######
model_output_cephalos_eu <- fit_gam_model_eu(country_resistance_cephalos_eu)
final_dataset<- model.frame(model_output_cephalos_eu$fr)
original_geo_levels <- levels(final_dataset$GEOID)

country_resistance_cephalos_eu$Year <- as.numeric(as.character(country_resistance_cephalos_eu$Year))
years_range <- seq(min(country_resistance_cephalos_eu$Year), max(country_resistance_cephalos_eu$Year), length.out = 100)
# Create prediction data
new_data <- expand.grid(Year = years_range, GEOID = original_geo_levels)
new_data$Year <- as.numeric(as.character(new_data$Year))  # Ensure Year is nume
states_spatial_fitered<- model_output_cephalos_eu$states_spatial_filtered
states_data_df <- as.data.frame(states_spatial_fitered@data)

original_geo_levels <- unique(new_data$GEOID)
GEOID_nameC <- states_data_df %>%
  dplyr::select(GEOID, NAME)
unique_GEOID_data <- GEOID_nameC %>%
  distinct(GEOID, .keep_all = TRUE)
# Expand pred_data to include every combination of Year and GEOID
pred_data <- expand.grid(Year = seq(min(new_data$Year), max(new_data$Year), by = 1),
                         GEOID = original_geo_levels)

# Predict from GAM using the created function
predictions <- gam_predictions(model_output_cephalos_eu$fr, newdata = pred_data)
predictions <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE)

# Usage of the function
derivatives_data <- derivatives_mh(model_output_cephalos_eu$fr, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001)
derivatives_data$growth_rate <- derivatives_data$first_derivative + (derivatives_data$second_derivative / derivatives_data$first_derivative)
merged_data <- merge(derivatives_data, unique_GEOID_data, by = "GEOID", all.x = TRUE)
Cephalos_changep_EU<- merged_data
# Plotting
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)
# Calculate the maximum and minimum values for setting y-axis limits
y_max <- max(merged_data$first_upper, na.rm = TRUE)
y_min <- min(merged_data$first_lower, na.rm = TRUE)
# If needed, you can add some padding to ensure all data points are within the view
padding <- (y_max - y_min) * 0.05  # 5% padding
y_max <- y_max + padding
y_min <- y_min - padding
# Plotting with adjustments
# Plotting with universal y-axis limits
p <- ggplot(data = merged_data, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Line at y=0
  scale_color_manual(values = palette) +  # Assuming palette is predefined
  scale_fill_manual(values = palette) +
  labs(title = "",
       x = "Year",
       y = "First derivative (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")   # Apply the same y-axis limits to all facets

# Display the plot
print(p)
ggsave(filename = "first_derivat_eu_cephalos.tiff", plot = p, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

#-----------------------------#
#PREDICTIONS: per country Europe
#-----------------------------#
min_y <- min(predictions$pred_lower, na.rm = TRUE)
max_y <- max(predictions$pred_upper, na.rm = TRUE)
# Generate the plot
ppred <- ggplot(data = predictions, aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +  # Use fixed scales for y-axis across all facets
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(min_y, max_y + 0), breaks = seq(0, max_y + 0, by = 10)) +
  labs(title = "",
       x = "Year",
       y = "Predicted 3rd generation cephalosporin-resistance (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size=8)
  )
# Display the plot
print(ppred)
ggsave(filename = "predictions_breakpoint_cephalos.tiff", plot = ppred, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

#-----------------------------#
#PREDICTIONS: All together
#-----------------------------#
# Step 1: Prepare the data for labels
label_data <- predictions %>%
  group_by(NAME) %>%
  filter(Year == max(Year)) %>%
  ungroup()

final_palette <- colorRampPalette(RColorBrewer::brewer.pal(10, "Set2"))(27)
# Step 2: Generate the plot
ppred_alle <- ggplot(data = predictions, aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_label_repel(data = label_data, aes(label = NAME, y = pred), 
                   point.padding = 0.2, nudge_x = 1, direction = 'y', 
                   size = 3.5, color = "black", fontface = "bold",
                   box.padding = 0.35, segment.color = "grey50",
                   fill = "white") +  # White background for labels
  scale_color_manual(values = final_palette) +
  scale_fill_manual(values = final_palette) +
  scale_x_continuous(breaks = 2004:2022) +
  scale_y_continuous(limits = c(NA, max(predictions$pred_upper) +5)) +
  scale_y_continuous(breaks = seq(0, 45, by = 5))+
  
  labs(title = "",
       x = "Year",
       y = "Predicted 3rd generation cephalosporin-resistance (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 360, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size=13),
    axis.title.y = element_text(size=13)
  )

# Print the plot
print(ppred_alle)
ggsave(filename = "pred_allEurop_toget_cephalos.tiff", plot = ppred_alle, device = "tiff", path = base_pathOut,
       width = 14, height = 8, dpi = 500, units = "in")

#-----------------------------#
# Plotting SECOND DERIVATIVE:
#-----------------------------#
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)
# Calculate the maximum and minimum values for setting y-axis limits
y_max <- max(merged_data$second_upper, na.rm = TRUE)
y_min <- min(merged_data$second_lower, na.rm = TRUE)
# If needed, you can add some padding to ensure all data points are within the view
padding <- (y_max - y_min) * 0.05  # 5% padding
y_max <- y_max + padding
y_min <- y_min - padding
# Plotting with adjustments
# Plotting with universal y-axis limits
p2<- ggplot(data = merged_data, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Line at y=0
  scale_color_manual(values = palette) +  # Assuming palette is predefined
  scale_fill_manual(values = palette) +
  labs(title = "",
       x = "Year",
       y = "Second derivative (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y") #+
#scale_y_continuous(limits = c(y_min, y_max))  # Apply the same y-axis limits to all facets

# Display the plot
print(p2)
ggsave(filename = "second_derivat_eu_cephalos.tiff", plot = p2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


#-----------------------------#
#GROWTH RATE: 
#-----------------------------#
# Assuming derivatives_data is already prepared and using mutate for efficiency
growth_rate_ci <- derivatives_data %>%
  rowwise() %>%
  mutate(
    growth_rate_results = list(simulate_growth_rate_ci(first_derivative, second_derivative, n_sim = 1000))
  ) %>%
  unnest_wider(growth_rate_results)


growth_rate_ci <- merge(growth_rate_ci, unique_GEOID_data, by = "GEOID", all.x = TRUE)

growth_rate_ci$doubling_times <-  log(2) / growth_rate_ci$growth_rate # Calculate doubling time
growth_rate_ci$halving_times =  log(0.5) / growth_rate_ci$growth_rate   # Calculate halving time

doubling_halving_times <- derivatives_data %>%
  mutate(
    Doubling_Time = log(2) / first_derivative,   # Calculate doubling time
    Halving_Time = log(0.5) / first_derivative   # Calculate halving time
  )

# Filter to remove infinite or undefined times (which occur if first_derivative is 0)
doubling_halving_times <- doubling_halving_times %>%
  filter(!is.infinite(Doubling_Time) & !is.na(Doubling_Time) & 
           !is.infinite(Halving_Time) & !is.na(Halving_Time))

p3<- ggplot(data = growth_rate_ci, aes(x = Year, y = growth_rate, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  #geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Line at y=0
  scale_color_manual(values = palette) +  # Assuming palette is predefined
  scale_fill_manual(values = palette) +
  labs(title = "",
       x = "Year",
       y = "Growth rate") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")  # Apply the same y-axis limits to all facets

# Display the plot
print(p3)
ggsave(filename = "growth_ratio_cephalosEU.tiff", plot = p3, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


# Plot Doubling Times
p_doubling <- ggplot(data = growth_rate_ci, aes(x = Year, y = doubling_times, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_point(aes(color = NAME), size = 2) +
  labs(title = "",
       x = "Year",
       y = "Doubling Time (in units of time)") +
  theme_minimal() +
  theme(text = element_text(size = 12, family = "Times New Roman"),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  facet_wrap(~NAME, scales = "free_y")

# Plot Halving Times
p_halving <- ggplot(data = growth_rate_ci, aes(x = Year, y = halving_times, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_point(aes(color = NAME), size = 2) +
  labs(title = "",
       x = "Year",
       y = "Halving Time (in units of time)") +
  theme_minimal() +
  theme(text = element_text(size = 12, family = "Times New Roman"),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  facet_wrap(~NAME, scales = "free_y")

# Print the plots
print(p_doubling)
print(p_halving)

Cephalos_changep_EU<- growth_rate_ci
ggsave(filename = "p_doublingEU_cephalos.tiff", plot = p_doubling, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "p_halvingEU_cephalos.tiff", plot = p_halving, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


######
#------------------------------------------------------------------------------#
#NEW, FIRST-LINE RESISTANCE GRAPH: ######
model_output_firstline_eu <- fit_gam_model_eu(country_resistance_firstline_eu)
final_dataset<- model.frame(model_output_firstline_eu$fr)
original_geo_levels <- levels(final_dataset$GEOID)

country_resistance_firstline_eu$Year <- as.numeric(as.character(country_resistance_firstline_eu$Year))
years_range <- seq(min(country_resistance_firstline_eu$Year), max(country_resistance_firstline_eu$Year), length.out = 100)
# Create prediction data
new_data <- expand.grid(Year = years_range, GEOID = original_geo_levels)
new_data$Year <- as.numeric(as.character(new_data$Year))  # Ensure Year is nume
# Convert Spatial*DataFrame to a regular data frame for easier manipulation
states_spatial_fitered<- model_output_firstline_eu$states_spatial_filtered
states_data_df <- as.data.frame(states_spatial_fitered@data)



original_geo_levels <- unique(new_data$GEOID)
GEOID_nameC <- states_data_df %>%
  dplyr::select(GEOID, NAME)
unique_GEOID_data <- GEOID_nameC %>%
  distinct(GEOID, .keep_all = TRUE)
# Expand pred_data to include every combination of Year and GEOID
pred_data <- expand.grid(Year = seq(min(new_data$Year), max(new_data$Year), by = 1),
                         GEOID = original_geo_levels)

# Predict from GAM using the created function
predictions <- gam_predictions(model_output_firstline_eu$fr, newdata = pred_data)
predictions <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE)

# Usage of the function
derivatives_data <- derivatives_mh(model_output_firstline_eu$fr, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001)
derivatives_data$growth_rate <- derivatives_data$first_derivative + (derivatives_data$second_derivative / derivatives_data$first_derivative)

merged_data <- merge(derivatives_data, unique_GEOID_data, by = "GEOID", all.x = TRUE)
firstline_changep_EU<- merged_data
# Plotting
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)
# Calculate the maximum and minimum values for setting y-axis limits
y_max <- max(merged_data$first_upper, na.rm = TRUE)
y_min <- min(merged_data$first_lower, na.rm = TRUE)
# If needed, you can add some padding to ensure all data points are within the view
padding <- (y_max - y_min) * 0.05  # 5% padding
y_max <- y_max + padding
y_min <- y_min - padding
# Plotting with adjustments
# Plotting with universal y-axis limits
p <- ggplot(data = merged_data, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Line at y=0
  scale_color_manual(values = palette) +  # Assuming palette is predefined
  scale_fill_manual(values = palette) +
  labs(title = "",
       x = "Year",
       y = "First derivative (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")   # Apply the same y-axis limits to all facets

# Display the plot
print(p)
ggsave(filename = "first_derivat_eu_firstline.tiff", plot = p, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

#-----------------------------#
#PREDICTIONS: per country Europe
#-----------------------------#
max_y <- max(predictions$pred_upper, na.rm = TRUE)
min_y <- min(predictions$pred_lower, na.rm = TRUE)
# Generate the plot
ppred <- ggplot(data = predictions, aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +  # Use fixed scales for y-axis across all facets
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(min_y, max_y + 0), breaks = seq(0, max_y + 0, by = 10)) +
  labs(title = "",
       x = "Year",
       y = "First-line antibiotic-resistance (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size=8)
  )
# Display the plot
print(ppred)
ggsave(filename = "predictions_breakpoint_firstline.tiff", plot = ppred, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

#-----------------------------#
#PREDICTIONS: All together
#-----------------------------#
# Step 1: Prepare the data for labels
label_data <- predictions %>%
  group_by(NAME) %>%
  filter(Year == max(Year)) %>%
  ungroup()

final_palette <- colorRampPalette(RColorBrewer::brewer.pal(10, "Set2"))(27)
# Step 2: Generate the plot
ppred_alle <- ggplot(data = predictions, aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_label_repel(data = label_data, aes(label = NAME, y = pred), 
                   point.padding = 0.2, nudge_x = 1, direction = 'y', 
                   size = 3.5, color = "black", fontface = "bold",
                   box.padding = 0.35, segment.color = "grey50",
                   fill = "white") +  # White background for labels
  scale_color_manual(values = final_palette) +
  scale_fill_manual(values = final_palette) +
  scale_x_continuous(breaks = 2004:2022) +
  scale_y_continuous(limits = c(NA, max(predictions$pred_upper) +5)) +
  scale_y_continuous(breaks = seq(0, 45, by = 5))+
  
  labs(title = "",
       x = "Year",
       y = "First-line antibiotic-resistance (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 360, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size=13),
    axis.title.y = element_text(size=13)
  )

# Print the plot
print(ppred_alle)
ggsave(filename = "pred_allEurop_toget_firstline.tiff", plot = ppred_alle, device = "tiff", path = base_pathOut,
       width = 14, height = 8, dpi = 500, units = "in")

#-----------------------------#
# Plotting SECOND DERIVATIVE:
#-----------------------------#
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)
# Calculate the maximum and minimum values for setting y-axis limits
y_max <- max(merged_data$second_upper, na.rm = TRUE)
y_min <- min(merged_data$second_lower, na.rm = TRUE)
# If needed, you can add some padding to ensure all data points are within the view
padding <- (y_max - y_min) * 0.05  # 5% padding
y_max <- y_max + padding
y_min <- y_min - padding
# Plotting with adjustments
# Plotting with universal y-axis limits
p2<- ggplot(data = merged_data, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Line at y=0
  scale_color_manual(values = palette) +  # Assuming palette is predefined
  scale_fill_manual(values = palette) +
  labs(title = "",
       x = "Year",
       y = "Second derivative (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y") #+
#scale_y_continuous(limits = c(y_min, y_max))  # Apply the same y-axis limits to all facets

# Display the plot
print(p2)
ggsave(filename = "second_derivat_eu_firstline.tiff", plot = p2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


#-----------------------------#
#GROWTH RATE: 
#-----------------------------#
# Assuming derivatives_data is already prepared and using mutate for efficiency
growth_rate_ci <- derivatives_data %>%
  rowwise() %>%
  mutate(
    growth_rate_results = list(simulate_growth_rate_ci(first_derivative, second_derivative, n_sim = 1000))
  ) %>%
  unnest_wider(growth_rate_results)


growth_rate_ci <- merge(growth_rate_ci, unique_GEOID_data, by = "GEOID", all.x = TRUE)

growth_rate_ci$doubling_times <-  log(2) / growth_rate_ci$growth_rate # Calculate doubling time
growth_rate_ci$halving_times =  log(0.5) / growth_rate_ci$growth_rate   # Calculate halving time

doubling_halving_times <- derivatives_data %>%
  mutate(
    Doubling_Time = log(2) / first_derivative,   # Calculate doubling time
    Halving_Time = log(0.5) / first_derivative   # Calculate halving time
  )

# Filter to remove infinite or undefined times (which occur if first_derivative is 0)
doubling_halving_times <- doubling_halving_times %>%
  filter(!is.infinite(Doubling_Time) & !is.na(Doubling_Time) & 
           !is.infinite(Halving_Time) & !is.na(Halving_Time))

p3<- ggplot(data = growth_rate_ci, aes(x = Year, y = growth_rate, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  #geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Line at y=0
  scale_color_manual(values = palette) +  # Assuming palette is predefined
  scale_fill_manual(values = palette) +
  labs(title = "",
       x = "Year",
       y = "Growth rate") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")  # Apply the same y-axis limits to all facets

# Display the plot
print(p3)
ggsave(filename = "growth_ratio_firstlineEU.tiff", plot = p3, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


# Plot Doubling Times
p_doubling <- ggplot(data = growth_rate_ci, aes(x = Year, y = doubling_times, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_point(aes(color = NAME), size = 2) +
  labs(title = "",
       x = "Year",
       y = "Doubling Time (in units of time)") +
  theme_minimal() +
  theme(text = element_text(size = 12, family = "Times New Roman"),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  facet_wrap(~NAME, scales = "free_y")

# Plot Halving Times
p_halving <- ggplot(data = growth_rate_ci, aes(x = Year, y = halving_times, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_point(aes(color = NAME), size = 2) +
  labs(title = "",
       x = "Year",
       y = "Halving Time (in units of time)") +
  theme_minimal() +
  theme(text = element_text(size = 12, family = "Times New Roman"),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  facet_wrap(~NAME, scales = "free_y")

# Print the plots
print(p_doubling)
print(p_halving)


ggsave(filename = "p_doublingEU_firstline.tiff", plot = p_doubling, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "p_halvingEU_firstline.tiff", plot = p_halving, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

firstline_changep_EU<- growth_rate_ci

######
#------------------------------------------------------------------------------#
#NEW, MDR RESISTANCE GRAPH: ######
model_output_mdr_eu <- fit_gam_model_mdr_eu(country_resistance_mdr_eu)
final_dataset<- model.frame(model_output_mdr_eu$sf)
original_geo_levels <- levels(final_dataset$GEOID)

country_resistance_mdr_eu$Year <- as.numeric(as.character(country_resistance_mdr_eu$Year))
years_range <- seq(min(country_resistance_mdr_eu$Year), max(country_resistance_mdr_eu$Year), length.out = 100)
# Create prediction data
new_data <- expand.grid(Year = years_range, GEOID = original_geo_levels)
new_data$Year <- as.numeric(as.character(new_data$Year))  # Ensure Year is nume

# Convert Spatial*DataFrame to a regular data frame for easier manipulation
states_spatial_filtered<- model_output_mdr_eu$states_spatial_filtered
states_data_df <- as.data.frame(states_spatial_filtered@data)

original_geo_levels <- unique(new_data$GEOID)
GEOID_nameC <- states_data_df %>%
  dplyr::select(GEOID, NAME)
unique_GEOID_data <- GEOID_nameC %>%
  distinct(GEOID, .keep_all = TRUE)
# Expand pred_data to include every combination of Year and GEOID
pred_data <- expand.grid(Year = seq(min(new_data$Year), max(new_data$Year), by = 1),
                         GEOID = original_geo_levels)

# Predict from GAM using the created function
predictions <- gam_predictions(model_output_mdr_eu$sf, newdata = pred_data)
predictions <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE)

#predictions_cephalos <- gam_predictions(model_output_cephalos_eu$fr, newdata = pred_data)

# Usage of the function
derivatives_data <- derivatives_mh(model_output_mdr_eu$sf, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001)
derivatives_data$growth_rate <- derivatives_data$first_derivative + (derivatives_data$second_derivative / derivatives_data$first_derivative)

merged_data <- merge(derivatives_data, unique_GEOID_data, by = "GEOID", all.x = TRUE)
mdr_changep_EU<- merged_data
#derivatives_data_cepha <- derivatives_mh(model_output_cephalos_eu$fr, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001)


# Plotting
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)
# Calculate the maximum and minimum values for setting y-axis limits
y_max <- max(merged_data$first_upper, na.rm = TRUE)
y_min <- min(merged_data$first_lower, na.rm = TRUE)
# If needed, you can add some padding to ensure all data points are within the view
padding <- (y_max - y_min) * 0.05  # 5% padding
y_max <- y_max + padding
y_min <- y_min - padding
# Plotting with adjustments
# Plotting with universal y-axis limits
p <- ggplot(data = merged_data, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Line at y=0
  scale_color_manual(values = palette) +  # Assuming palette is predefined
  scale_fill_manual(values = palette) +
  labs(title = "",
       x = "Year",
       y = "First derivative (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")   # Apply the same y-axis limits to all facets

# Display the plot
print(p)
ggsave(filename = "first_derivat_eu_mdr.tiff", plot = p, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

#-----------------------------#
#PREDICTIONS: per country Europe
#-----------------------------#
max_y <- max(predictions$pred_upper, na.rm = TRUE)
min_y <- min(predictions$pred_lower, na.rm = TRUE)
# Generate the plot
ppred <- ggplot(data = predictions, aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +  # Use fixed scales for y-axis across all facets
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(min_y, max_y + 0), breaks = seq(0, max_y + 0, by = 10)) +
  labs(title = "",
       x = "Year",
       y = "Predicted multidrug-resistance (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size=8)
  )
# Display the plot
print(ppred)
ggsave(filename = "predictions_breakpoint_mdr.tiff", plot = ppred, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

#-----------------------------#
#PREDICTIONS: All together
#-----------------------------#
# Step 1: Prepare the data for labels
label_data <- predictions %>%
  group_by(NAME) %>%
  filter(Year == max(Year)) %>%
  ungroup()

final_palette <- colorRampPalette(RColorBrewer::brewer.pal(10, "Set2"))(27)
# Step 2: Generate the plot
ppred_alle <- ggplot(data = predictions, aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_label_repel(data = label_data, aes(label = NAME, y = pred), 
                   point.padding = 0.2, nudge_x = 1, direction = 'y', 
                   size = 3.5, color = "black", fontface = "bold",
                   box.padding = 0.35, segment.color = "grey50",
                   fill = "white") +  # White background for labels
  scale_color_manual(values = final_palette) +
  scale_fill_manual(values = final_palette) +
  scale_x_continuous(breaks = 2004:2022) +
  scale_y_continuous(limits = c(NA, max(predictions$pred_upper) +5)) +
  scale_y_continuous(breaks = seq(0, 45, by = 5))+
  
  labs(title = "",
       x = "Year",
       y = "Predicted multidrug-resistance (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 360, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size=13),
    axis.title.y = element_text(size=13)
  )

# Print the plot
print(ppred_alle)
ggsave(filename = "pred_allEurop_toget_mdr.tiff", plot = ppred_alle, device = "tiff", path = base_pathOut,
       width = 14, height = 8, dpi = 500, units = "in")

#-----------------------------#
# Plotting SECOND DERIVATIVE:
#-----------------------------#
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)
# Calculate the maximum and minimum values for setting y-axis limits
y_max <- max(merged_data$second_upper, na.rm = TRUE)
y_min <- min(merged_data$second_lower, na.rm = TRUE)
# If needed, you can add some padding to ensure all data points are within the view
padding <- (y_max - y_min) * 0.05  # 5% padding
y_max <- y_max + padding
y_min <- y_min - padding
# Plotting with adjustments
# Plotting with universal y-axis limits
p2<- ggplot(data = merged_data, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Line at y=0
  scale_color_manual(values = palette) +  # Assuming palette is predefined
  scale_fill_manual(values = palette) +
  labs(title = "",
       x = "Year",
       y = "Second derivative (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y") #+
#scale_y_continuous(limits = c(y_min, y_max))  # Apply the same y-axis limits to all facets

# Display the plot
print(p2)
ggsave(filename = "second_derivat_eu_mdr.tiff", plot = p2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


#-----------------------------#
#GROWTH RATE: 
#-----------------------------#
# Assuming derivatives_data is already prepared and using mutate for efficiency
growth_rate_ci <- derivatives_data %>%
  rowwise() %>%
  mutate(
    growth_rate_results = list(simulate_growth_rate_ci(first_derivative, second_derivative, n_sim = 1000))
  ) %>%
  unnest_wider(growth_rate_results)


growth_rate_ci <- merge(growth_rate_ci, unique_GEOID_data, by = "GEOID", all.x = TRUE)

growth_rate_ci$doubling_times <-  log(2) / growth_rate_ci$growth_rate # Calculate doubling time
growth_rate_ci$halving_times =  log(0.5) / growth_rate_ci$growth_rate   # Calculate halving time

doubling_halving_times <- derivatives_data %>%
  mutate(
    Doubling_Time = log(2) / first_derivative,   # Calculate doubling time
    Halving_Time = log(0.5) / first_derivative   # Calculate halving time
  )

# Filter to remove infinite or undefined times (which occur if first_derivative is 0)
doubling_halving_times <- doubling_halving_times %>%
  filter(!is.infinite(Doubling_Time) & !is.na(Doubling_Time) & 
           !is.infinite(Halving_Time) & !is.na(Halving_Time))

p3<- ggplot(data = growth_rate_ci, aes(x = Year, y = growth_rate, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  #geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Line at y=0
  scale_color_manual(values = palette) +  # Assuming palette is predefined
  scale_fill_manual(values = palette) +
  labs(title = "",
       x = "Year",
       y = "Growth rate") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")  # Apply the same y-axis limits to all facets

# Display the plot
print(p3)
ggsave(filename = "growth_ratio_mdrEU.tiff", plot = p3, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


# Plot Doubling Times
p_doubling <- ggplot(data = growth_rate_ci, aes(x = Year, y = doubling_times, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_point(aes(color = NAME), size = 2) +
  labs(title = "",
       x = "Year",
       y = "Doubling Time (in units of time)") +
  theme_minimal() +
  theme(text = element_text(size = 12, family = "Times New Roman"),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  facet_wrap(~NAME, scales = "free_y")

# Plot Halving Times
p_halving <- ggplot(data = growth_rate_ci, aes(x = Year, y = halving_times, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_point(aes(color = NAME), size = 2) +
  labs(title = "",
       x = "Year",
       y = "Halving Time (in units of time)") +
  theme_minimal() +
  theme(text = element_text(size = 12, family = "Times New Roman"),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  facet_wrap(~NAME, scales = "free_y")

# Print the plots
print(p_doubling)
print(p_halving)


ggsave(filename = "p_doublingEU_mdr.tiff", plot = p_doubling, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "p_halvingEU_mdr.tiff", plot = p_halving, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

mdr_changep_EU<- growth_rate_ci


######
#------------------------------------------------------------------------------#



#-------------------------------------------------------------------#
###SUBGROUP ANALYSES, European countries: Gender and Age groups.
#-------------------------------------------------------------------#

######

fit_gam_model_euGen <- function(country_resistance_carbap_euGen) {
  europe_shapefile <- ne_countries(continent = "Europe", returnclass = "sf")
  europe_shapefile <- st_make_valid(europe_shapefile)
  # Rename and modify sovereignt to NAME
  europe_shapefile <- mutate(europe_shapefile, NAME = sovereignt)
  europe_shapefile <- mutate(europe_shapefile,
                             NAME = case_when(
                               NAME == "Slovakia" ~ "Slovak Republic",
                               NAME == "Republic of Serbia" ~ "Serbia",
                               NAME == "Czechia" ~ "Czech Republic",
                               TRUE ~ NAME
                             ))
  
  # Ensure all geometries are valid
  states_shapefile<-europe_shapefile
  country_resistance_carbap_euGen$NAME<- country_resistance_carbap_euGen$Country
  # Merge spatial data with resistance data by 'NAME'
  states_merged <- merge(states_shapefile, country_resistance_carbap_euGen, by = "NAME", all.x = TRUE)
  states_merged <- states_merged %>%
    dplyr::filter(NAME %in% countries_list_eu_includ)
  
  # Create unique geometries for neighborhood structure
  states_spatial_unique <- states_merged[!duplicated(states_merged$NAME), ]
  states_spatial_unique$GEOID <- as.numeric(factor(states_spatial_unique$NAME))
  
  states_spatial_filtered <- states_merged %>%
    mutate(
      GEOID = as.numeric(factor(NAME, levels = unique(NAME)))
    )
  
  if (!inherits(states_spatial_filtered, "Spatial")) {
    states_spatial_filtered <- as(states_spatial_filtered, "Spatial")
  }
  
  # 8. Remove duplicates based on GEOID to create unique state geometries for neighborhood structure
  states_spatial_unique <- states_spatial_filtered[!duplicated(states_spatial_filtered$GEOID), ]
  # 9. Recreate the neighborhood structure using unique state geometries
  nb <- poly2nb(states_spatial_unique, row.names = states_spatial_unique$GEOID)
  # 10. Assign region names to the neighborhood list
  names(nb) <- states_spatial_unique$GEOID
  
  # Clean and ensure data alignment with neighborhood structure
  states_spatial_filtered@data$GEOID <- factor(states_spatial_filtered@data$GEOID, levels = attr(nb, "region.id"))
  states_spatial_filtered@data$Resistance <- as.numeric(states_spatial_filtered@data$Resistance)
  states_spatial_filtered@data$AMR_Positive <- as.numeric(states_spatial_filtered@data$AMR_Positive)
  states_spatial_filtered@data$Total_Isolates <- as.numeric(states_spatial_filtered@data$Total_Isolates)
  states_spatial_filtered@data$Year <- as.numeric(as.character(states_spatial_filtered@data$Year))
  #states_spatial_filtered@data <- states_spatial_filtered@data[complete.cases(states_spatial_filtered@data), ]
  
  # Set control parameters for GAM
  ctrl <- gam.control(nthreads = 6)  # Set parallel threads for faster computation
  states_spatial_filtered@data$Gender <- as.factor(states_spatial_filtered@data$Gender)
  
  
  # Fit the MRF model with spatio-temporal smoothing
  model_gamx1 <- gam(Resistance ~ 
                       s(Year, m=3, k=k, bs = "tp") +  # Assuming 'Year' can support more knots
                       s(GEOID, bs = 'mrf', xt = list(nb = nb)) +
                       te(Year, GEOID, bs=c("tp", "mrf"), m=c(3, NA), xt=list(Year=NULL, GEOID=list(nb=nb))) +
                       Gender +  # Treat 'Gender' as a factor (categorical effect)
                       Year:Gender,  
                     data = states_spatial_filtered@data,
                     family = binomial, 
                     select = TRUE,
                     method = "REML", 
                     weights = Total_Isolates)
  
  return(list(fr=model_gamx1, states_spatial_filtered=states_spatial_filtered))
}
fit_gam_model_mdr_euGen <- function(state_resistance_carbap) {
  # Load Spatial Data and Resistance Data
  europe_shapefile <- ne_countries(continent = "Europe", returnclass = "sf")
  europe_shapefile <- st_make_valid(europe_shapefile)
  # Rename and modify sovereignt to NAME
  europe_shapefile <- mutate(europe_shapefile, NAME = sovereignt)
  europe_shapefile <- mutate(europe_shapefile,
                             NAME = case_when(
                               NAME == "Slovakia" ~ "Slovak Republic",
                               NAME == "Republic of Serbia" ~ "Serbia",
                               NAME == "Czechia" ~ "Czech Republic",
                               TRUE ~ NAME
                             ))
  
  # Ensure all geometries are valid
  states_shapefile<-europe_shapefile
  state_resistance_carbap$NAME<- state_resistance_carbap$Country
  # Merge spatial data with resistance data by 'NAME'
  states_merged <- merge(states_shapefile, state_resistance_carbap, by = "NAME", all.x = TRUE)
  states_merged <- states_merged %>%
    dplyr::filter(NAME %in% countries_list_eu_includ)
  # Create unique geometries for neighborhood structure
  states_spatial_unique <- states_merged[!duplicated(states_merged$NAME), ]
  #states_spatial_unique$GEOID <- as.numeric(factor(states_spatial_unique$NAME))
  states_spatial_filtered <- states_merged %>%
    dplyr::filter(NAME %in% countries_list_eu_includ)
  
  states_spatial_filtered <- states_spatial_filtered %>%
    mutate(
      GEOID = as.numeric(factor(NAME, levels = unique(NAME)))
    )
  
  if (!inherits(states_spatial_filtered, "Spatial")) {
    states_spatial_filtered <- as(states_spatial_filtered, "Spatial")
  }
  
  
  # 8. Remove duplicates based on GEOID to create unique state geometries for neighborhood structure
  states_spatial_unique <- states_spatial_filtered[!duplicated(states_spatial_filtered$GEOID), ]
  # 9. Recreate the neighborhood structure using unique state geometries
  nb <- poly2nb(states_spatial_unique, row.names = states_spatial_unique$GEOID)
  # 10. Assign region names to the neighborhood list
  names(nb) <- states_spatial_unique$GEOID
  
  # 11. Clean and ensure data alignment with neighborhood structure
  # Convert 'GEOID' to factor and ensure levels match the neighborhood structure
  states_spatial_filtered@data$GEOID <- factor(states_spatial_filtered@data$GEOID, levels = attr(nb, "region.id"))
  # Convert 'Resistance' to numeric
  states_spatial_filtered@data$MDR <- as.numeric(states_spatial_filtered@data$MDR)
  states_spatial_filtered@data$MDR_Positive <- as.numeric(states_spatial_filtered@data$MDR_Positive)
  states_spatial_filtered@data$Total_Isolates <- as.numeric(states_spatial_filtered@data$Total_Isolates)
  # Convert 'Year' to numeric and handle missing data
  states_spatial_filtered@data$Year <- as.numeric(as.character(states_spatial_filtered@data$Year))
  # 12. Recreate neighborhood structure using the cleaned and filtered data
  states_spatial_unique <- states_spatial_filtered[!duplicated(states_spatial_filtered$GEOID), ]
  nb <- poly2nb(states_spatial_unique, row.names = states_spatial_unique$GEOID)
  names(nb) <- states_spatial_unique$GEOID
  
  # 13. Set control parameters for GAM
  ctrl <- gam.control(nthreads = 6)  # Set parallel threads for faster computation
  states_spatial_filtered@data$Gender <- as.factor(states_spatial_filtered@data$Gender)
  
  
  # Fit the MRF model with spatio-temporal smoothing
  model_gamx1 <- gam(MDR ~ 
                       s(Year, m=3, k=k, bs = "tp") +  # Assuming 'Year' can support more knots
                       s(GEOID, bs = 'mrf', xt = list(nb = nb)) +
                       te(Year, GEOID, bs=c("tp", "mrf"), m=c(3, NA), xt=list(Year=NULL, GEOID=list(nb=nb))) +
                       Gender +  # Treat 'Gender' as a factor (categorical effect)
                       Year:Gender,  
                     data = states_spatial_filtered@data,
                     family = binomial, 
                     select = TRUE,
                     method = "REML", 
                     weights = Total_Isolates)
  return(list(sf=model_gamx1,  states_spatial_filtered= states_spatial_filtered))
}

#OUTPUT OF THE RESULTS SPATIO-TEMPORAL MODELS:

model_output_carbap_euGen <- fit_gam_model_euGen(country_resistance_carbap_euGen)
model_output_cephalos_euGen <- fit_gam_model_euGen(country_resistance_cephalos_euGen)
model_output_firstline_euGen <- fit_gam_model_euGen(country_resistance_firstline_euGen)
model_output_mdr_euGen <-fit_gam_model_mdr_euGen(country_resistance_mdr_euGen)

fit_gam_model_euAgeg <- function(country_resistance_carbap_euGen) {
  europe_shapefile <- ne_countries(continent = "Europe", returnclass = "sf")
  europe_shapefile <- st_make_valid(europe_shapefile)
  # Rename and modify sovereignt to NAME
  europe_shapefile <- mutate(europe_shapefile, NAME = sovereignt)
  europe_shapefile <- mutate(europe_shapefile,
                             NAME = case_when(
                               NAME == "Slovakia" ~ "Slovak Republic",
                               NAME == "Republic of Serbia" ~ "Serbia",
                               NAME == "Czechia" ~ "Czech Republic",
                               TRUE ~ NAME
                             ))
  
  # Ensure all geometries are valid
  states_shapefile<-europe_shapefile
  country_resistance_carbap_euGen$NAME<- country_resistance_carbap_euGen$Country
  # Merge spatial data with resistance data by 'NAME'
  states_merged <- merge(states_shapefile, country_resistance_carbap_euGen, by = "NAME", all.x = TRUE)
  states_merged <- states_merged %>%
    dplyr::filter(NAME %in% countries_list_eu_includ)
  # Create unique geometries for neighborhood structure
  states_spatial_unique <- states_merged[!duplicated(states_merged$NAME), ]
  #states_spatial_unique$GEOID <- as.numeric(factor(states_spatial_unique$NAME))

  states_spatial_filtered <- states_merged %>%
    dplyr::filter(NAME %in% countries_list_eu_includ)
  
  states_spatial_filtered <- states_spatial_filtered %>%
    mutate(
      GEOID = as.numeric(factor(NAME, levels = unique(NAME)))
    )
  
  if (!inherits(states_spatial_filtered, "Spatial")) {
    states_spatial_filtered <- as(states_spatial_filtered, "Spatial")
  }
  
  
  # 8. Remove duplicates based on GEOID to create unique state geometries for neighborhood structure
  states_spatial_unique <- states_spatial_filtered[!duplicated(states_spatial_filtered$GEOID), ]
  # 9. Recreate the neighborhood structure using unique state geometries
  nb <- poly2nb(states_spatial_unique, row.names = states_spatial_unique$GEOID)
  # 10. Assign region names to the neighborhood list
  names(nb) <- states_spatial_unique$GEOID
  
  # Clean and ensure data alignment with neighborhood structure
  states_spatial_filtered@data$GEOID <- factor(states_spatial_filtered@data$GEOID, levels = attr(nb, "region.id"))
  states_spatial_filtered@data$Resistance <- as.numeric(states_spatial_filtered@data$Resistance)
  states_spatial_filtered@data$AMR_Positive <- as.numeric(states_spatial_filtered@data$AMR_Positive)
  states_spatial_filtered@data$Total_Isolates <- as.numeric(states_spatial_filtered@data$Total_Isolates)
  states_spatial_filtered@data$Year <- as.numeric(as.character(states_spatial_filtered@data$Year))
  #states_spatial_filtered@data <- states_spatial_filtered@data[complete.cases(states_spatial_filtered@data), ]
  
  # Set control parameters for GAM
  ctrl <- gam.control(nthreads = 6)  # Set parallel threads for faster computation
  states_spatial_filtered@data$Agegroup <- as.factor(states_spatial_filtered@data$Agegroup)
  
  
  # Fit the MRF model with spatio-temporal smoothing
  model_gamx1 <- gam(Resistance ~ 
                       s(Year, m=3, k=10, bs = "tp") +  # Assuming 'Year' can support more knots
                       s(GEOID, bs = 'mrf', xt = list(nb = nb)) +
                       te(Year, GEOID, bs=c("tp", "mrf"), m=c(3, NA), xt=list(Year=NULL, GEOID=list(nb=nb))) +
                       Agegroup +  # Treat 'Gender' as a factor (categorical effect)
                       Year:Agegroup,  
                     data = states_spatial_filtered@data,
                     family = binomial, 
                     select = TRUE,
                     method = "REML", 
                     weights = Total_Isolates)
  
  return(list(fr=model_gamx1, states_spatial_filtered=states_spatial_filtered))
}
fit_gam_model_mdr_euAgeg <- function(state_resistance_carbap) {
  # Load Spatial Data and Resistance Data
  europe_shapefile <- ne_countries(continent = "Europe", returnclass = "sf")
  europe_shapefile <- st_make_valid(europe_shapefile)
  # Rename and modify sovereignt to NAME
  europe_shapefile <- mutate(europe_shapefile, NAME = sovereignt)
  europe_shapefile <- mutate(europe_shapefile,
                             NAME = case_when(
                               NAME == "Slovakia" ~ "Slovak Republic",
                               NAME == "Republic of Serbia" ~ "Serbia",
                               NAME == "Czechia" ~ "Czech Republic",
                               TRUE ~ NAME
                             ))
  
  # Ensure all geometries are valid
  states_shapefile<-europe_shapefile
  state_resistance_carbap$NAME<- state_resistance_carbap$Country
  # Merge spatial data with resistance data by 'NAME'
  states_merged <- merge(states_shapefile, state_resistance_carbap, by = "NAME", all.x = TRUE)
  states_merged <- states_merged %>%
    dplyr::filter(NAME %in% countries_list_eu_includ)
  # Create unique geometries for neighborhood structure
  states_spatial_unique <- states_merged[!duplicated(states_merged$NAME), ]
  #states_spatial_unique$GEOID <- as.numeric(factor(states_spatial_unique$NAME))
  
  states_spatial_filtered <- states_merged %>%
    dplyr::filter(NAME %in% countries_list_eu_includ)
  
  states_spatial_filtered <- states_spatial_filtered %>%
    mutate(
      GEOID = as.numeric(factor(NAME, levels = unique(NAME)))
    )
  
  if (!inherits(states_spatial_filtered, "Spatial")) {
    states_spatial_filtered <- as(states_spatial_filtered, "Spatial")
  }
  
  
  # 8. Remove duplicates based on GEOID to create unique state geometries for neighborhood structure
  states_spatial_unique <- states_spatial_filtered[!duplicated(states_spatial_filtered$GEOID), ]
  # 9. Recreate the neighborhood structure using unique state geometries
  nb <- poly2nb(states_spatial_unique, row.names = states_spatial_unique$GEOID)
  # 10. Assign region names to the neighborhood list
  names(nb) <- states_spatial_unique$GEOID
  
  # 11. Clean and ensure data alignment with neighborhood structure
  # Convert 'GEOID' to factor and ensure levels match the neighborhood structure
  states_spatial_filtered@data$GEOID <- factor(states_spatial_filtered@data$GEOID, levels = attr(nb, "region.id"))
  # Convert 'Resistance' to numeric
  states_spatial_filtered@data$MDR <- as.numeric(states_spatial_filtered@data$MDR)
  states_spatial_filtered@data$MDR_Positive <- as.numeric(states_spatial_filtered@data$MDR_Positive)
  states_spatial_filtered@data$Total_Isolates <- as.numeric(states_spatial_filtered@data$Total_Isolates)
  # Convert 'Year' to numeric and handle missing data
  states_spatial_filtered@data$Year <- as.numeric(as.character(states_spatial_filtered@data$Year))
  # 12. Recreate neighborhood structure using the cleaned and filtered data
  states_spatial_unique <- states_spatial_filtered[!duplicated(states_spatial_filtered$GEOID), ]
  nb <- poly2nb(states_spatial_unique, row.names = states_spatial_unique$GEOID)
  names(nb) <- states_spatial_unique$GEOID
  
  # 13. Set control parameters for GAM
  ctrl <- gam.control(nthreads = 6)  # Set parallel threads for faster computation
  # 14. Fit the MRF model with spatio-temporal smoothing
  states_spatial_filtered@data$Agegroup <- as.factor(states_spatial_filtered@data$Agegroup)
  model_gamx1 <- gam(MDR ~ 
                       s(Year, m=3, k=10, bs = "tp") +   # Smooth term for Year
                       s(GEOID, bs = 'mrf', xt = list(nb = nb)) +  # Smooth term for GEOID
                       te(Year, GEOID, bs=c("tp", "mrf"), m=c(3, NA), xt=list(Year=NULL, GEOID=list(nb=nb))) +
                       Agegroup +  
                       Year:Agegroup, # Interaction term
                     data = states_spatial_filtered@data,
                     family = binomial, 
                     select = TRUE,
                     method = "REML", 
                     weights = Total_Isolates)
  return(list(sf=model_gamx1,  states_spatial_filtered= states_spatial_filtered))
}

#OUTPUT OF THE RESULTS SPATIO-TEMPORAL MODELS:
model_output_carbap_euAgeg <- fit_gam_model_euAgeg(country_resistance_carbap_euAgeg)
model_output_cephalos_euAgeg <- fit_gam_model_euAgeg(country_resistance_cephalos_euAgeg)
model_output_firstline_euAgeg <- fit_gam_model_euAgeg(country_resistance_firstline_euAgeg)
model_output_mdr_euAgeg <-fit_gam_model_mdr_euAgeg(country_resistance_mdr_euAgeg)


#-----Predict for GENDER --------#
#CHECK NOW HOW TO PREDICT below for GENDER:
original_geo_levels <- levels(final_dataset$GEOID)
country_resistance_carbap_eu$Year <- as.numeric(as.character(country_resistance_carbap_eu$Year))
years_range <- seq(min(country_resistance_carbap_eu$Year), max(country_resistance_carbap_eu$Year), length.out = 100)
# Function to process data, predict, and plot results
compute_and_plot_predictionsGen <- function(data, model_output, years_range, original_geo_levels) {
  # Convert Year to numeric
  data$Year <- as.numeric(as.character(data$Year))
  
  # Prepare new data for predictions
  new_data <- expand.grid(Year = years_range, GEOID = original_geo_levels)
  new_data$Year <- as.numeric(as.character(new_data$Year))
  
  # Obtain unique geographic levels and names
  states_data_df <- as.data.frame(model_output$states_spatial_filtered@data)
  GEOID_nameC <- states_data_df %>%
    dplyr::select(GEOID, NAME) %>%
    distinct(GEOID, .keep_all = TRUE)
  
  # Expand prediction data to include every combination of Year and GEOID
  genders <- c("Male", "Female")
  pred_data <- expand.grid(Year = seq(min(new_data$Year), max(new_data$Year), by = 1),
                           GEOID = unique(new_data$GEOID),
                           Gender = genders)
  
  # Generate predictions using the provided GAM model
  predictions <- gam_predictions(model_output$fr, newdata = pred_data)
  predictions <- merge(predictions, GEOID_nameC, by = "GEOID", all.x = TRUE)
  
  # Merge predictions back with the main data to include gender and other variables
  #predictions <- merge(predictions, data, by = c("Year", "GEOID"))
  
  # Calculate average predictions and confidence intervals by country and gender
  predictions_summary <- predictions %>%
    group_by(Year, Gender) %>%
    summarise(avg_pred = mean(pred),
              avg_pred_lower = mean(pred_lower),
              avg_pred_upper = mean(pred_upper), .groups = 'drop')
  
  # Plotting the results
  ggplot(predictions_summary, aes(x = Year, y = avg_pred, group = Gender, color = Gender)) +
    geom_line(size = 1) +  # Make the line a bit thicker
    geom_ribbon(aes(ymin = avg_pred_lower, ymax = avg_pred_upper, fill = Gender), alpha = 0.2, color = NA) +
    scale_color_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e")) +
    scale_fill_manual(values = c("Male" = adjustcolor("#1f77b4", alpha.f = 0.2), 
                                 "Female" = adjustcolor("#ff7f0e", alpha.f = 0.2))) +  # Use the same colors as the lines but more transparent
    labs(title = "",
         subtitle = "",
         x = "Year",
         y = "Predicted resistance (%)") +
    theme_minimal(base_size = 14) +  # Base font size adjustment for better readability
    theme(legend.position = "bottom",  # Adjust legend positioning
          legend.title = element_blank(),  # Remove legend title
          plot.title = element_text(face = "bold", hjust = 0.5),  # Center and bold the plot title
          plot.subtitle = element_text(hjust = 0.5),  # Center the subtitle
          axis.text = element_text(size = 12),  # Adjust axis text size
          axis.title = element_text(size = 11))  # Adjust axis title size
  
}
compute_and_plot_predictionsGenmdr <- function(data, model_output, years_range, original_geo_levels) {
  # Convert Year to numeric
  data$Year <- as.numeric(as.character(data$Year))
  
  # Prepare new data for predictions
  new_data <- expand.grid(Year = years_range, GEOID = original_geo_levels)
  new_data$Year <- as.numeric(as.character(new_data$Year))
  
  # Obtain unique geographic levels and names
  states_data_df <- as.data.frame(model_output$states_spatial_filtered@data)
  GEOID_nameC <- states_data_df %>%
    dplyr::select(GEOID, NAME) %>%
    distinct(GEOID, .keep_all = TRUE)
  
  # Expand prediction data to include every combination of Year and GEOID
  genders <- c("Male", "Female")
  pred_data <- expand.grid(Year = seq(min(new_data$Year), max(new_data$Year), by = 1),
                           GEOID = unique(new_data$GEOID),
                           Gender = genders)
  
  # Generate predictions using the provided GAM model
  predictions <- gam_predictions(model_output$sf, newdata = pred_data)
  predictions <- merge(predictions, GEOID_nameC, by = "GEOID", all.x = TRUE)
  
  # Merge predictions back with the main data to include gender and other variables
  #predictions <- merge(predictions, data, by = c("Year", "GEOID"))
  
  # Calculate average predictions and confidence intervals by country and gender
  predictions_summary <- predictions %>%
    group_by(Year, Gender) %>%
    summarise(avg_pred = mean(pred),
              avg_pred_lower = mean(pred_lower),
              avg_pred_upper = mean(pred_upper), .groups = 'drop')
  
  # Plotting the results
  ggplot(predictions_summary, aes(x = Year, y = avg_pred, group = Gender, color = Gender)) +
    geom_line(size = 1) +  # Make the line a bit thicker
    geom_ribbon(aes(ymin = avg_pred_lower, ymax = avg_pred_upper, fill = Gender), alpha = 0.2, color = NA) +
    scale_color_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e")) +
    scale_fill_manual(values = c("Male" = adjustcolor("#1f77b4", alpha.f = 0.2), 
                                 "Female" = adjustcolor("#ff7f0e", alpha.f = 0.2))) +  # Use the same colors as the lines but more transparent
    labs(title = "",
         subtitle = "",
         x = "Year",
         y = "Predicted resistance (%)") +
    theme_minimal(base_size = 14) +  # Base font size adjustment for better readability
    theme(legend.position = "bottom",  # Adjust legend positioning
          legend.title = element_blank(),  # Remove legend title
          plot.title = element_text(face = "bold", hjust = 0.5),  # Center and bold the plot title
          plot.subtitle = element_text(hjust = 0.5),  # Center the subtitle
          axis.text = element_text(size = 12),  # Adjust axis text size
          axis.title = element_text(size = 11))  # Adjust axis title size
  
}

plot_resultGen_carb_EU <- compute_and_plot_predictionsGen(country_resistance_carbap_euGen, model_output_carbap_euGen, years_range, original_geo_levels)
plot_resultGen_carb_EU <- plot_resultGen_carb_EU + ylab("Carbapenem-resistance (%)")
plot_resultGen_cephalos_EU <- compute_and_plot_predictionsGen(country_resistance_cephalos_euGen, model_output_cephalos_euGen, years_range, original_geo_levels)
plot_resultGen_cephalos_EU <- plot_resultGen_cephalos_EU  + ylab("3rd generation cephalosporin-resistance (%)")
plot_resultGen_firstline_EU <- compute_and_plot_predictionsGen(country_resistance_firstline_euGen, model_output_firstline_euGen, years_range, original_geo_levels)
plot_resultGen_firstline_EU <- plot_resultGen_firstline_EU + ylab("First-line antibiotic-resistance (%)")
plot_resultGen_mdr_EU <- compute_and_plot_predictionsGenmdr(country_resistance_mdr_euGen, model_output_mdr_euGen, years_range, original_geo_levels)
plot_resultGen_mdr_EU <- plot_resultGen_mdr_EU + ylab("Multidrug resistance (%)")

# Arrange the plots
combined_plotki <- (plot_resultGen_mdr_EU + plot_resultGen_firstline_EU) /
  (plot_resultGen_cephalos_EU + plot_resultGen_carb_EU)
# Collect the guides and add tags
combined_plotki <- combined_plotki + 
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A',
                  tag_prefix = "", 
                  tag_suffix = ".")
# Modify legend position and theme settings for clarity
combined_plotki <- combined_plotki & theme(
  legend.position = "bottom",
  legend.justification = "center",
  legend.box.background = element_rect(color = "white", linetype = "solid"),
  legend.background = element_rect(fill = "white", color = NA),
  plot.margin = unit(c(0,0,0,0), "lines")
)
ggsave(filename = "predictions_gender_EU.tiff", plot = combined_plotki, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


#-----Predict for Agegroups --------#
#CHECK NOW HOW TO PREDICT below for Agegroup:
original_geo_levels <- levels(final_dataset$GEOID)
country_resistance_carbap_eu$Year <- as.numeric(as.character(country_resistance_carbap_eu$Year))
years_range <- seq(min(country_resistance_carbap_eu$Year), max(country_resistance_carbap_eu$Year), length.out = 100)
# Function to process data, predict, and plot results
compute_and_plot_predictionsAgeg <- function(data, model_output, years_range, original_geo_levels) {
  # Convert Year to numeric
  data$Year <- as.numeric(as.character(data$Year))
  
  # Prepare new data for predictions
  new_data <- expand.grid(Year = years_range, GEOID = original_geo_levels)
  new_data$Year <- as.numeric(as.character(new_data$Year))
  
  # Obtain unique geographic levels and names
  states_data_df <- as.data.frame(model_output$states_spatial_filtered@data)
  GEOID_nameC <- states_data_df %>%
    dplyr::select(GEOID, NAME) %>%
    distinct(GEOID, .keep_all = TRUE)
  
  # Expand prediction data to include every combination of Year and GEOID
  agegroups <- c(0, 1, 2)
  new_data$Year <- as.numeric(as.character(new_data$Year))
  pred_data <- expand.grid(Year = seq(min(new_data$Year), max(new_data$Year), by = 1),
                           GEOID = unique(new_data$GEOID),
                           Agegroup = agegroups)
  
  # Generate predictions using the provided GAM model
  predictions <- gam_predictions(model_output$fr, newdata = pred_data)
  predictions <- merge(predictions, GEOID_nameC, by = "GEOID", all.x = TRUE)
  
  # Merge predictions back with the main data to include gender and other variables
  #predictions <- merge(predictions, data, by = c("Year", "GEOID"))
  
  # Calculate average predictions and confidence intervals by country and gender
  predictions_summary <- predictions %>%
    group_by(Year, Agegroup) %>%
    summarise(avg_pred = mean(pred),
              avg_pred_lower = mean(pred_lower),
              avg_pred_upper = mean(pred_upper), .groups = 'drop')
  
  predictions_summary$Agegroup <- factor(predictions_summary$Agegroup)
  predictions_summary$Agegroup <- factor(predictions_summary$Agegroup,
                                         levels = c("0", "1", "2"),
                                         labels = c("≤18yo", "19≤ and ≤64", "≥65"))
  # Now create the plot
  ggplot(predictions_summary, aes(x = Year, y = avg_pred, group = Agegroup, color = Agegroup)) +
    geom_line(size = 1) +  # Make the line a bit thicker
    geom_ribbon(aes(ymin = avg_pred_lower, ymax = avg_pred_upper, fill = Agegroup), alpha = 0.2, color = NA) +
    scale_color_manual(values = c("≤18yo" = "#1f77b4", "19≤ and ≤64" = "#ff9f8e", "≥65" = "#FDFD96")) +
    scale_fill_manual(values = c("≤18yo" = adjustcolor("#1f77b4", alpha.f = 0.2), 
                                 "19≤ and ≤64" = adjustcolor("#ff9f8e", alpha.f = 0.2),
                                 "≥65" = adjustcolor("#FDFD96", alpha.f = 0.2))) +
    labs(title = "",
         subtitle = "",
         x = "Year",
         y = "Predicted resistance (%)") +
    theme_minimal(base_size = 14) +  # Base font size adjustment for better readability
    theme(legend.position = "bottom",  # Adjust legend positioning
          legend.title = element_blank(),  # Remove legend title
          plot.title = element_text(face = "bold", hjust = 0.5),  # Center and bold the plot title
          plot.subtitle = element_text(hjust = 0.5),  # Center the subtitle
          axis.text = element_text(size = 12),  # Adjust axis text size
          axis.title = element_text(size = 11))  # Adjust axis title size
  
}
compute_and_plot_predictionsAgegmdr <- function(data, model_output, years_range, original_geo_levels) {
  # Convert Year to numeric
  data$Year <- as.numeric(as.character(data$Year))
  
  # Prepare new data for predictions
  new_data <- expand.grid(Year = years_range, GEOID = original_geo_levels)
  new_data$Year <- as.numeric(as.character(new_data$Year))
  
  # Obtain unique geographic levels and names
  states_data_df <- as.data.frame(model_output$states_spatial_filtered@data)
  GEOID_nameC <- states_data_df %>%
    dplyr::select(GEOID, NAME) %>%
    distinct(GEOID, .keep_all = TRUE)
  
  # Expand prediction data to include every combination of Year and GEOID
  agegroups <- c(0, 1, 2)
  pred_data <- expand.grid(Year = seq(min(new_data$Year), max(new_data$Year), by = 1),
                           GEOID = unique(new_data$GEOID),
                           Agegroup = agegroups)
  
  # Generate predictions using the provided GAM model
  predictions <- gam_predictions(model_output$sf, newdata = pred_data)
  predictions <- merge(predictions, GEOID_nameC, by = "GEOID", all.x = TRUE)
  
  # Merge predictions back with the main data to include gender and other variables
  #predictions <- merge(predictions, data, by = c("Year", "GEOID"))
  
  # Calculate average predictions and confidence intervals by country and gender
  predictions_summary <- predictions %>%
    group_by(Year, Agegroup) %>%
    summarise(avg_pred = mean(pred),
              avg_pred_lower = mean(pred_lower),
              avg_pred_upper = mean(pred_upper), .groups = 'drop')
  predictions_summary$Agegroup <- factor(predictions_summary$Agegroup,
                                         levels = c("0", "1", "2"),
                                         labels = c("≤18yo", "19≤ and ≤64", "≥65"))
  # Plotting the results

  ggplot(predictions_summary, aes(x = Year, y = avg_pred, group = Agegroup, color = Agegroup)) +
    geom_line(size = 1) +  # Make the line a bit thicker
    geom_ribbon(aes(ymin = avg_pred_lower, ymax = avg_pred_upper, fill = Agegroup), alpha = 0.2, color = NA) +
    scale_color_manual(values = c("≤18yo" = "#1f77b4", "19≤ and ≤64" = "#ff9f8e", "≥65" = "#FDFD96")) +
    scale_fill_manual(values = c("≤18yo" = adjustcolor("#1f77b4", alpha.f = 0.2), 
                                 "19≤ and ≤64" = adjustcolor("#ff9f8e", alpha.f = 0.2),
                                 "≥65" = adjustcolor("#FDFD96", alpha.f = 0.2))) +
    labs(title = "",
         subtitle = "",
         x = "Year",
         y = "Predicted resistance (%)") +
    theme_minimal(base_size = 14) +  # Base font size adjustment for better readability
    theme(legend.position = "bottom",  # Adjust legend positioning
          legend.title = element_blank(),  # Remove legend title
          plot.title = element_text(face = "bold", hjust = 0.5),  # Center and bold the plot title
          plot.subtitle = element_text(hjust = 0.5),  # Center the subtitle
          axis.text = element_text(size = 12),  # Adjust axis text size
          axis.title = element_text(size = 11))  # Adjust axis title size
  
  
}

plot_resultAgeg_carb_EU <- compute_and_plot_predictionsAgeg(country_resistance_carbap_euAgeg, model_output_carbap_euAgeg, years_range, original_geo_levels)
plot_resultAgeg_carb_EU <- plot_resultAgeg_carb_EU + ylab("Carbapenem-resistance (%)")
plot_resultAgeg_cephalos_EU <- compute_and_plot_predictionsAgeg(country_resistance_cephalos_euAgeg, model_output_cephalos_euAgeg, years_range, original_geo_levels)
plot_resultAgeg_cephalos_EU <- plot_resultAgeg_cephalos_EU  + ylab("3rd generation cephalosporin-resistance (%)")
plot_resultAgeg_firstline_EU <- compute_and_plot_predictionsAgeg(country_resistance_firstline_euAgeg, model_output_firstline_euAgeg, years_range, original_geo_levels)
plot_resultAgeg_firstline_EU <- plot_resultAgeg_firstline_EU + ylab("First-line antibiotic-resistance (%)")
plot_resultAgeg_mdr_EU <- compute_and_plot_predictionsAgegmdr(country_resistance_mdr_euAgeg, model_output_mdr_euAgeg, years_range, original_geo_levels)
plot_resultAgeg_mdr_EU <- plot_resultAgeg_mdr_EU + ylab("Multidrug resistance (%)")

# Arrange the plots
combined_plotki2 <- (plot_resultAgeg_mdr_EU + plot_resultAgeg_firstline_EU) /
  (plot_resultAgeg_cephalos_EU + plot_resultAgeg_carb_EU)
# Collect the guides and add tags
combined_plotki2 <- combined_plotki2 + 
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A',
                  tag_prefix = "", 
                  tag_suffix = ".")
# Modify legend position and theme settings for clarity
combined_plotki2 <- combined_plotki2 & theme(
  legend.position = "bottom",
  legend.justification = "center",
  legend.box.background = element_rect(color = "white", linetype = "solid"),
  legend.background = element_rect(fill = "white", color = NA),
  plot.margin = unit(c(0,0,0,0), "lines")
)
ggsave(filename = "predictions_Agegroup_EU.tiff", plot = combined_plotki2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")




#----------------------------------------------------------------#
#Table outputs: 

model_output_carbap_euGen2 <- model_output_carbap_euGen$fr
model_output_cephalos_euGen2<- model_output_cephalos_euGen$fr
model_output_firstline_euGen2<- model_output_firstline_euGen$fr
model_output_mdr_euGen2 <- model_output_mdr_euGen$sf

model_output_carbap_euAgeg2 <- model_output_carbap_euAgeg$fr
model_output_cephalos_euAgeg2<- model_output_cephalos_euAgeg$fr
model_output_firstline_euAgeg2<- model_output_firstline_euAgeg$fr
model_output_mdr_euAgeg2 <- model_output_mdr_euAgeg$sf

# Apply the function to each model
summary_carbap_eu <- extract_model_summary_eu(model_output_carbap_euGen2)
summary_cephalos_eu <- extract_model_summary_eu(model_output_cephalos_euGen2)
summary_firstline_eu <- extract_model_summary_eu(model_output_firstline_euGen2)
summary_mdr_eu <- extract_model_summary_eu(model_output_mdr_euGen2)

summary_carbap_eu2 <- extract_model_summary_eu(model_output_carbap_euAgeg2)
summary_cephalos_eu2 <- extract_model_summary_eu(model_output_cephalos_euAgeg2)
summary_firstline_eu2 <- extract_model_summary_eu(model_output_firstline_euAgeg2)
summary_mdr_eu2 <- extract_model_summary_eu(model_output_mdr_euAgeg2)

# Adding model identifiers
summary_carbap_eu$model <- "Carbapenem Resistance"
summary_cephalos_eu$model <- "Cephalosporin Resistance"
summary_firstline_eu$model <- "First-line Antibiotic Resistance"
summary_mdr_eu$model <- "MDR Resistance"
summary_carbap_eu2$model <- "Carbapenem Resistance"
summary_cephalos_eu2$model <- "Cephalosporin Resistance"
summary_firstline_eu2$model <- "First-line Antibiotic Resistance"
summary_mdr_eu2$model <- "MDR Resistance"
# Combining all summaries into one dataframe
combined_results_eu <- bind_rows(summary_carbap_eu, summary_cephalos_eu, summary_firstline_eu, summary_mdr_eu, summary_carbap_eu2, summary_cephalos_eu2, summary_firstline_eu2, summary_mdr_eu2)
# Optionally, select and rename columns for clarity
final_results_eu2 <- combined_results_eu %>%
  dplyr::select(Model = model, Term = term, Estimate = estimate, EDF = edf, Ref.DF = ref.df, Std.Error = std.error,
                Statistic = statistic, `P.Value` = p.value)
######

#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
  
#------------------------------------------------------------------------#
#Covariates; multiple data sources at country-level #####
  # Load the RData file
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/0_UniversityofOxford/Vivli/data/data/Outputs/ATLAS/")
load("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/0_UniversityofOxford/Vivli/data/data/data_covariates/gbd_covars.RDATA")
gbd_covars<-x

load("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/0_UniversityofOxford/Vivli/data/data/data_covariates/cmmid_comorb.RDATA")
cmmid_comorb<-x

load("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/0_UniversityofOxford/Vivli/data/data/data_covariates/dysentery_inc.RDATA")
dysentery_inc<-x

load("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/0_UniversityofOxford/Vivli/data/data/data_covariates/gbd_clin_inc.RDATA")
gbd_clin_inc<-x

load("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/0_UniversityofOxford/Vivli/data/data/data_covariates/incomedata.RDATA")
incomedata<-x

load("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/0_UniversityofOxford/Vivli/data/data/data_covariates/inform.work.RDATA")
inform.work<-x

load("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/0_UniversityofOxford/Vivli/data/data/Outputs/ATLAS/popdata.RDATA")
popdata<-x

load("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/0_UniversityofOxford/Vivli/data/data/data_covariates/vacc_cov_comb.RDATA")
vacc_cov_comb<-x

load("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/0_UniversityofOxford/Vivli/data/data/data_covariates/sepsis_incidence.RDATA")
sepsis_incidence<-x

load("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/0_UniversityofOxford/Vivli/data/data/data_covariates/wb_covars.RDATA")
wb_covars<-x

load("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/0_UniversityofOxford/Vivli/data/data/data_covariates/gram1518.didaware.RDATA")
gramatb<-x


load("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/0_UniversityofOxford/Vivli/data/data/data_covariates/gram.did.RDATA")
gram_did<-x


#######
#Random forest; Lasso/Ridge/Elastic Net regression

###-------------------------------------###
#USE ANTIBIOTIC CONSUMPTION FROM GRAM ###
###-------------------------------------###
#Let's use GRAM antibiotic consumption! // gramatb
#filter gram variables:
#####
gram_did2 <- gram_did %>%
  mutate(ihme_country_name = case_when(
    ihme_country_name == "Czechia" ~ "Czech Republic",
    ihme_country_name == "Slovakia" ~ "Slovak Republic",
    TRUE ~ ihme_country_name  # Leaves other country names unchanged
  ))

gram_did2<- gram_did2 %>%
  filter(ihme_country_name %in% countries_list_eu_includ,  
         year >= 2004 & year <= 2022)

gram_did2neth<- gram_did2 %>%
  filter(ihme_country_name=="Netherlands")
gram_did2$NAME <- gram_did2$ihme_country_name
gram_did2$Year <- gram_did2$year
#Dataframes with the data:
dataframe_list <- list(
  mdr_changep_EU = mdr_changep_EU,
  firstline_changep_EU = firstline_changep_EU,
  Carb_changep_EU = Carb_changep_EU,
  Cephalos_changep_EU = Cephalos_changep_EU
)

modify_data_frames <- function(df) {
  # Create new columns
  df <- df %>%
    mutate(
      ChangePsecond = derivative_breakpoint,
      ChangePsecondfirst = as.numeric(first_derivative_sign_change)  # Convert logical TRUE/FALSE to 1/0
    )
  return(df)
}

#Create main variables per dataframe:
dataframe_list <- lapply(dataframe_list, modify_data_frames)

dataframe_list$mdr_changep_EU$ChangePsecond
dataframe_list$mdr_changep_EU$ChangePsecondfirst

dataframe_list$firstline_changep_EU$ChangePsecond
dataframe_list$firstline_changep_EU$ChangePsecondfirst

dataframe_list$Carb_changep_EU$ChangePsecond
dataframe_list$Carb_changep_EU$ChangePsecondfirst

dataframe_list$Cephalos_changep_EU$ChangePsecond
dataframe_list$Cephalos_changep_EU$ChangePsecondfirst
#######
#average DID per year in Europe.
average_did_per_year <- gram_did2 %>%
  group_by(Year) %>%  # Group data by the Year column
  summarise(Average_DID = mean(did_total, na.rm = TRUE))

average_did_per_year2 <- gram_did2 %>%
  group_by(Year) %>%  
  summarise(Average_beta = mean(j01d_beta, na.rm = TRUE))  # Calculate average, ignoring NA values


average_did_per_country <- gram_did2 %>%
  group_by(NAME) %>%  # Group data by the country name
  summarise(Average_DID = mean(did_total, na.rm = TRUE)) %>%  # Calculate the average DID, ignoring NA values
  ungroup() %>%  # Remove the grouping
  arrange(desc(Average_DID))  # Sort the results in descending order based on average DID

# Display the country with the highest average DID
highest_consumption_country <- head(average_did_per_country, 1)
print(highest_consumption_country)




average_did_per_country <- gram_did2 %>%
  group_by(NAME) %>%  # Group data by the country name
  summarise(Average_DID = mean(j01d_beta, na.rm = TRUE)) %>%  # Calculate the average DID, ignoring NA values
  ungroup() %>%  # Remove the grouping
  arrange(desc(Average_DID))  # Sort the results in descending order based on average DID

# Display the country with the highest average DID
highest_consumption_country <- head(average_did_per_country, 1)
print(highest_consumption_country)



#one of the main variables: growth_rate
#-----------------------------------------------------------#
#REGRESSION ANALYSIS:
#-----------------------------------------------------------#
gram_did2x<- gram_did2
gram_did2x$NAME<- gram_did2$ihme_country_name
gram_did2x$Year<- gram_did2$year
#Outcomes: ChangePsecondfirst  ChangePsecond
predictions_carb_EU <- gam_predictions(model_output_carbap_eu$fr, newdata = pred_data)
predictions_carb_EU <- merge(predictions_carb_EU, unique_GEOID_data, by = "GEOID", all.x = TRUE)
predictions_cephalos_EU <- gam_predictions(model_output_cephalos_eu$fr, newdata = pred_data)
predictions_cephalos_EU <- merge(predictions_cephalos_EU, unique_GEOID_data, by = "GEOID", all.x = TRUE)
predictions_firstline_EU <- gam_predictions(model_output_firstline_eu$fr, newdata = pred_data)
predictions_firstline_EU <- merge(predictions_firstline_EU, unique_GEOID_data, by = "GEOID", all.x = TRUE)
predictions_mdr_EU <- gam_predictions(model_output_mdr_eu$sf, newdata = pred_data)
predictions_mdr_EU <- merge(predictions_mdr_EU, unique_GEOID_data, by = "GEOID", all.x = TRUE)

merged_dataCarb_EU <- merge(gram_did2x, Carb_predictions_grat, by = c("NAME", "Year"))
merged_dataCephalos_EU <- merge(gram_did2x, Ceph_predictions_grat, by = c("NAME", "Year"))
merged_datafirstline_EU <- merge(gram_did2x,Firstl_predictions_grat, by = c("NAME", "Year"))
merged_datamdr_EU <- merge(gram_did2x, MDR_predictions_grat, by = c("NAME", "Year"))

library(glmnet)
library(knitr)
antibiotics_xx <- c("did_total", "j01a_tet", "j01c_pen", "j01d_beta", "j01e_sulpha", "j01f_macro", "j01g_amino", "j01m_quin")
##Generating data for ATBs & controls for GAMs#######
merged_dataCarb_EU <- merged_dataCarb_EU %>% 
  arrange(NAME,Year)
merged_dataCarb_EU <- merged_dataCarb_EU %>%
  group_by(NAME) %>%
  mutate(pct_change_did_total = (did_total - lag(did_total)) / lag(did_total) * 100) %>%
  mutate(pct_change_did_total = replace_na(pct_change_did_total, 0))
merged_dataCarb_EU <- merged_dataCarb_EU %>%
  group_by(NAME) %>%
  mutate(pct_change_blact_total = (j01d_beta - lag(j01d_beta)) / lag(j01d_beta) * 100) %>%
  mutate(pct_change_blact_total = replace_na(pct_change_blact_total, 0))


########

#MODELS splines and all #######
gam_predictions2 <- function(gam_model, newdata) {
  prediction <- predict(gam_model, newdata = newdata, type = "response", se = TRUE)
  prediction_df <- data.frame(
    pred = prediction$fit ,  # Assuming you need predictions as percentages
    pred_lower = (prediction$fit - 1.96 * prediction$se.fit) ,  # 95% CI Lower Bound
    pred_upper = (prediction$fit + 1.96 * prediction$se.fit) ,  # 95% CI Upper Bound
    se = prediction$se.fit
  )
  prediction_df <- cbind(newdata, prediction_df)
  return(prediction_df)
}

#----------------------------------#
#PREDICTIONS FOR DID
#----------------------------------#
######
model_output_merged_dataCarb_EU <- gam(growth_rate2 ~ 
                     s(did_total, m=3, k=k, bs = "tp") +
                     s(pct_change_did_total, m=3, k=k, bs = "tp") +
                     t2(did_total, pct_change_did_total),
                   data = merged_dataCarb_EU,
                   family = gaussian, 
                   select = TRUE,
                   method = "REML")

set.seed(123)  # for reproducibility
summary(merged_dataCarb_EU$j01d_beta)
summary(merged_dataCarb_EU$pct_change_blact_total)
new_data_frame <- data.frame(
  did_total = rnorm(n = 100, mean = mean(merged_dataCarb_EU$did_total), sd = sd(merged_dataCarb_EU$did_total)),  # Normal distribution
  pct_change_did_total = rnorm(n = 1000, mean = mean(merged_dataCarb_EU$pct_change_did_total), sd = sd(merged_dataCarb_EU$pct_change_did_total))  # Normal distribution
)
new_data_frame <- expand.grid(
  #did_total = seq(0, 50, by = 1)  # Sequence from 0 to 40 by 1
  pct_change_did_total = seq(-30, 30, by = 1)  # Sequence from 0 to 30 by 1
)
predictions_carb_EU_did <- gam_predictions2(model_output_merged_dataCarb_EU, new_data_frame)
predictions_carb_EU_did$pred<- predictions_carb_EU_did$pred/100
predictions_carb_EU_did$pred_lower<- predictions_carb_EU_did$pred_lower/100
predictions_carb_EU_did$pred_upper<- predictions_carb_EU_did$pred_upper/100


line_color <- "#FF3D1F"  # Watermelon color for the line
ci_color <- "#FFB2A8"  # Lighter shade of watermelon for the confidence interval shading
  # Assuming your dataframe is named predictions_carb_EU_did and includes the necessary columns
  # Plotting with a professional color palette
loess_pred <- loess(pred ~ pct_change_did_total, data = predictions_carb_EU_did, control = loess.control(surface = "direct"))
loess_lower <- loess(pred_lower ~ pct_change_did_total, data = predictions_carb_EU_did, control = loess.control(surface = "direct"))
loess_upper <- loess(pred_upper ~ pct_change_did_total, data = predictions_carb_EU_did, control = loess.control(surface = "direct"))
did_total_seq <- data.frame(pct_change_did_total = seq(-40, 40, length.out = 400))
smooth_predictions <- data.frame(
  did_total = did_total_seq$pct_change_did_total,
  smooth_pred = predict(loess_pred, newdata = did_total_seq),
  smooth_lower = predict(loess_lower, newdata = did_total_seq),
  smooth_upper = predict(loess_upper, newdata = did_total_seq)
)
library(ggplot2)
# Plotting with shaded area between smoothed confidence intervals
atb_did_total<- ggplot(smooth_predictions, aes(x = did_total)) +
  geom_ribbon(aes(ymin = smooth_lower, ymax = smooth_upper), fill = "#8FCB2F", alpha = 0.3) +
  geom_line(aes(y = smooth_pred), color = "#007100", size = 1.2, alpha = 0.8) +
  labs(title = "H. Antibiotic use & CRE using GAMs",
       x = "Percentage change of overall antibiotic usage (%)",
       y = "Predicted growth rate for CRE (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(color = "black", face = "bold", size=14, hjust = 0.0),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "white", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "#555555"),
    axis.title = element_text(color = "#555555"),
    legend.position = "none",
    text = element_text(size = 12, family = "Times New Roman"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)  ) +
  xlim(-40, 40)

######



#----------------------------------#
#PREDICTIONS FOR Beta-lactams ######
#----------------------------------#
######
model_output_merged_dataCarb_EU2<- gam(growth_rate2 ~ 
                                        s(j01d_beta, m=3, k=k, bs = "tp") +
                                        s(pct_change_blact_total, m=3, k=k, bs = "tp") +
                                        t2(j01d_beta, pct_change_blact_total),
                                      data = merged_dataCarb_EU,
                                      family = gaussian, 
                                      select = TRUE,
                                      method = "REML")

new_data_frame <- data.frame(
  j01d_beta = rnorm(n = 100, mean = mean(merged_dataCarb_EU$j01d_beta), sd = sd(merged_dataCarb_EU$j01d_beta)),  # Normal distribution
  pct_change_blact_total = rnorm(n = 100, mean = mean(merged_dataCarb_EU$pct_change_blact_total), sd = sd(merged_dataCarb_EU$pct_change_blact_total))  # Normal distribution
)
predictions_carb_EU_blac <- gam_predictions2(model_output_merged_dataCarb_EU2, new_data_frame)
predictions_carb_EU_blac$j01d_beta <- ifelse(predictions_carb_EU_blac$j01d_beta < 0, 0, predictions_carb_EU_blac$j01d_beta)


line_color <- "#FF3D1F"  # Watermelon color for the line
  ci_color <- "#FFB2A8"  # Lighter shade of watermelon for the confidence interval shading
    # Assuming your dataframe is named predictions_carb_EU_did and includes the necessary columns
  # Plotting with a professional color palette
loess_pred <- loess(pred ~ pct_change_blact_total, data = predictions_carb_EU_blac, control = loess.control(surface = "direct"))
loess_lower <- loess(pred_lower ~ pct_change_blact_total, data = predictions_carb_EU_blac, control = loess.control(surface = "direct"))
loess_upper <- loess(pred_upper ~ pct_change_blact_total, data = predictions_carb_EU_blac, control = loess.control(surface = "direct"))
j01d_beta_seq <- data.frame(pct_change_blact_total = seq(0, 10, length.out = 300))
smooth_predictions2 <- data.frame(
    j01d_beta = j01d_beta_seq$pct_change_blact_total,
    smooth_pred = predict(loess_pred, newdata = j01d_beta_seq),
    smooth_lower = predict(loess_lower, newdata = j01d_beta_seq),
    smooth_upper = predict(loess_upper, newdata = j01d_beta_seq)
  )
  library(ggplot2)
  # Plotting with shaded area between smoothed confidence intervals
  atb_j01d_beta_total<- ggplot(smooth_predictions2, aes(x = j01d_beta)) +
    geom_ribbon(aes(ymin = smooth_lower, ymax = smooth_upper), fill = "#8FCB2F", alpha = 0.3) +
    geom_line(aes(y = smooth_pred), color = "#007100", size = 1.2, alpha = 0.8) +
    labs(title = "H. β-lactam use & CRE using GAMs",
         x = "β-lactam use (j01d) in DIDs",
         y = "Predicted carbapenem-resistance (%)") +
    theme_minimal() +
    theme(
      plot.title = element_text(color = "#555555", face = "bold", size=14, hjust = 0.0),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "white"),
      panel.grid.major = element_line(color = "white", size = 0.5),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "#555555"),
      axis.title = element_text(color = "#555555"),
      legend.position = "none",
      text = element_text(size = 12, family = "Times New Roman"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)  ) +
    xlim(-10, 10)
  
  

  

######


# Initialize a list to store summaries
fit_and_summarize_models <- function(data, antibiotics) {
  # Initialize an empty list to store model summaries
  model_summaries <- list()
  # Loop through antibiotics and fit models
  for (ab in antibiotics) {
    full_formula <- as.formula(paste("ChangePsecondfirst ~ s(Year, m=3, k=10, bs = 'tp') + s(", ab, ", m=3, bs='tp')"))
    reduced_formula <- as.formula(paste("ChangePsecondfirst ~ s(", ab, ", m=3, bs='tp')"))
    
    # Fit the full and reduced models using GAM
    full_model <- gam(full_formula, data = data, family = binomial(), select = TRUE, method = "REML")
    reduced_model <- gam(reduced_formula, data = data, family = binomial(), select = TRUE, method = "REML")
    
    # Extract and store formatted summaries
    model_summaries[[ab]] <- list(
      full_model_summary = extract_model_summary_eu(full_model),
      reduced_model_summary = extract_model_summary_eu(reduced_model)
    )
  }
  return(model_summaries)
}
print_model_summaries <- function(summaries) {
  for (ab in names(summaries)) {
    cat("\n--- Summary for", ab, "---\n")
    print(summaries[[ab]]$full_model_summary)
    print(summaries[[ab]]$reduced_model_summary)
  }
}
# Call the function to fit models and extract summaries
model_summaries_Carb_EU <- fit_and_summarize_models(merged_dataCarb_EU, antibiotics_xx)
model_summaries_Cephalos_EU <- fit_and_summarize_models(merged_dataCephalos_EU, antibiotics_xx)
model_summaries_firstline_EU <- fit_and_summarize_models(merged_datafirstline_EU, antibiotics_xx)
model_summaries_mdr_EU <- fit_and_summarize_models(merged_datamdr_EU, antibiotics_xx)


#####





#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#

carbapenem_heatmap2ab<- carbapenem_heatmap2 +
  theme(plot.title = element_text(hjust = 0, size = 14, face = "bold")) + labs(title = "A. CRE(%) 2004-22, Europe") 


combined_plot_fig_examplar2 <- (carbapenem_heatmap2ab + plot_carbapenem_eu_fig1) / 
  (ppred_alle2x + plot_resultGen_carb_EU2 ) / 
  (netherlands_plot + p3nethelands) / 
  (  did_total_plot+atb_did_total)  
  


#Configure orientation
combined_plot_fig_examplar <- combined_plot_fig_examplar2 + 
  plot_layout(ncol = 2, nrow = 2, heights = rep(2, 2)) + facet_wrap(~variable, scales = "free")

# Print the combined plot to check the layout
print(combined_plot_fig_examplar)
# Save the combined plot
ggsave(filename = "Figure_examplar.tiff", plot = combined_plot_fig_examplar, device = "tiff", path = base_pathOut,
       width = 17, height = 10, dpi = 500, units = "in")


#PANEL B#######
carbapenem_res_eu_avg <- calculate_averages(carbapenem_resistance_by_year_country_eu, "Country", Carbapenem_res)
# Apply the function to calculate averages specifically for netherlands now
carbapenem_res_eu_avg_ireland <- calculate_averages(
  carbapenem_resistance_by_year_country_eu %>% filter(Country == 'Netherlands'), 
  "Country", 
  Carbapenem_res)
# Define a color palette for Europe
color_palette <- c("Europe" = "#d95f02")

plot_carbapenem_europe <- function(data, data_ireland) {
  ggplot() +
    geom_violin(data = data, aes(x = Year_group, y = mean_resistance, fill = "Europe"), trim = TRUE, alpha = 0.6) +
    geom_boxplot(data = data, aes(x = Year_group, y = mean_resistance, group = Year_group),
                 width = 0.2, alpha = 0.6, outlier.shape = NA, color = "black", fill = "white") +
    geom_point(data = data, aes(x = Year_group, y = mean_resistance, color = "Each European country included"), 
               position = position_dodge(width = 0.8), size = 1.5, alpha = 0.5) +
    scale_fill_manual(values = c("Europe" = "#d95f02")) +
    scale_color_manual(values = c("Each European country included" = "#d95f02"), name = "", labels = "Each European country included") +
    labs(title = "B. CRE trends in Europe",
         subtitle = "",
         x = "Year group", y = "Carbapenem-resistance (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 360, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.position = "bottom",  # Position the legend at the bottom
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12)) +  # Size for the legend title
    scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 5), oob = scales::oob_squish) +
    guides(fill = guide_none(), color = guide_legend(title = ""))  # Customize legends
}

# Generate and print the plot
plot_carbapenem_eu_fig1 <- plot_carbapenem_europe(carbapenem_res_eu_avg, carbapenem_res_eu_avg_ireland)
plot_carbapenem_eu_fig1<- plot_carbapenem_eu_fig1 +  theme(
  text = element_text(size = 12, family = "Times New Roman"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_rect(fill = "white", color = "black"),
  strip.text = element_text(size = 14, face = "bold")
)
plot_carbapenem_eu_fig1<-plot_carbapenem_eu_fig1 + 
  theme(plot.title = element_text(hjust = 0))

######
#PANEL C#####
netherlands_plot <- netherlands_plot +labs(title = "E. Predicted CRE, exemplar") +theme(plot.title = element_text(face = "bold"))
##Growth rate######
p3nethelands<- p3nethelands +labs(title = "F. Growth rate, exemplar") +theme(plot.title = element_text(face = "bold"))
#####
###NethelerlandsDID total####
gram_did2neth<- gram_did2 %>%
  filter(ihme_country_name=="Netherlands")
library(hrbrthemes)  # For theme_ipsum()
# Constants for the plot
didTotalColor <- "#E76C3C"  # A reddish color for DID total
  j01dBetaColor <- "#0073C2FF"  # A bluish color for j01d_beta
    # Compute the coefficient for scaling
  coeff <- max(gram_did2neth$did_total) / max(gram_did2neth$j01d_beta)
  # Create the plot
did_total_plot <- ggplot(gram_did2neth, aes(x = Year)) +
    geom_line(aes(y = did_total), size = 2, color = didTotalColor) +  # Line for did_total
    geom_point(aes(y = did_total), color = "black", size = 3, shape = 21, fill = didTotalColor) +  # Points for did_total
    geom_line(aes(y = j01d_beta * coeff), size = 2, color = j01dBetaColor) +  # Scaled line for j01d_beta
    geom_point(aes(y = j01d_beta * coeff), color = "black", size = 3, shape = 21, fill = j01dBetaColor) +  # Points for did_total
    scale_y_continuous(
      name = "Total antibiotic usage in DIDs",
      limits = c(4, 12),  # Manual limits for the primary axis
      breaks = seq(4, 12, by = 1),  # Adjust breaks as necessary
      sec.axis = sec_axis(
        trans = ~ . / coeff, 
        name = "Other β-lactams (j01d) in DIDs",
        breaks = seq(0.2, 0.6, by = 0.1)  # Adjust breaks as necessary for the secondary axis
      )) +
    #theme_ipsum() +  # Clean theme with good defaults
    theme(
      axis.title.y = element_text(color = didTotalColor, size = 14, hjust = 0.5),
      axis.title.y.right = element_text(color = j01dBetaColor, size = 14,hjust = 0.5),
      plot.title = element_text(hjust = 0.0),
      axis.title.x= element_text(size=13, hjust = 0.5),
      text = element_text(size = 12, family = "Times New Roman")
    ) +
    labs(
      title = "G. Antibiotic use, exemplar",
      x = "Year"
    ) +scale_x_continuous(breaks = seq(2004, 2018, by = 2)) 


did_total_plot<- did_total_plot+ facet_wrap(~ihme_country_name, scales = "free_y") +theme(plot.title = element_text(face = "bold"))
did_total_plot<- did_total_plot + theme(plot.background = element_rect(fill = "white", color = "white"),  # Sets the plot background to white
panel.background = element_rect(fill = "white", color = "white"),    panel.grid.major = element_blank(),  # Removes major grid lines
panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),     strip.background = element_rect(fill = "white", color = "black", size = 0.4),  # White background with black border for facet labels
strip.text = element_text(color = "black", size = 14, face = "bold"), axis.ticks = element_blank() # Black text for facet labels
)

  
did_total_plot

#####
####second derivative: #####
p2netherl <- p2netherl +labs(title = "G. Second derivative, exemplar") +theme(plot.title = element_text(face = "bold"))
######
#a######
ppred_alle2x <- ppred_alle2x +labs(title = "C. Predicted CRE Europe") +theme(plot.title = element_text(face = "bold"), legend.position = "bottom",  # Moves the legend to the bottom
                                                                                                                               legend.box = "horizontal")
ppred_alle2x<- ppred_alle2x+  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(face = "bold", size=14),
    legend.position = "bottom",  # Moves the legend to the bottom
    legend.box = "horizontal",
    legend.title=element_blank())
#GENgraph####

plot_resultGen_carb_EU2<- plot_resultGen_carb_EU+ labs(title = "D. Predicted CRE Europe, by Gender") +theme(plot.title = element_text(face = "bold"))
plot_resultGen_carb_EU2<- plot_resultGen_carb_EU2 +  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(face = "bold", size=14))

plot_resultGen_carb_EU2 <- plot_resultGen_carb_EU2 +
  theme(legend.position = "bottom",  # Moves the legend to the bottom
        legend.box = "horizontal", 
        legend.title = element_blank())
#####

######
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#

#-----------------------------------------------------------#
#-----------------------------------------------------------#
#-----------------------------------------------------------#
#-----------------------------------------------------------#
#-----------------------------------------------------------#
#-----------------------------------------------------------#
#-----------------------------------------------------------#
#-----------------------------------------------------------#
#-----------------------------------------------------------#

#API & State-level data in the US #######
library(tidyverse)
library(tidyquant)
library(censusapi)
apis <- listCensusApis()
# Example: Get population data by state
#Indicators_us_state <- getCensus(
#  name = "timeseries/poverty/saipe",
#  vars = c("NAME", "SAEPOVRTALL_PT"),
#  region = "state:*",
#  time = "from 2004",
#  key= '859acf11fb09589d811b3bb7d39b3bc2b503dfe1')
#######




