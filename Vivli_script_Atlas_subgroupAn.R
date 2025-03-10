
#                                                                              #


#. EUROPE, below.

#                                                                              #


#Subgroup analyses - Gender-specific change points [Europe]

#------------------------------------------------------------------------------#
#NEW, CARBAPENEM RESISTANCE GRAPH: ######
# Step 1: Prepare filtered dataset with Gender trimmed
country_resistance_carbap_eu <- country_resistance_carbap_euGen %>%
  mutate(Gender = str_trim(Gender)) %>%  # Trim whitespace
  filter(Gender == "Male")  # Filter for Male

# Step 2: Set up the model output
model_output_carbap_eu <- model_output_carbap_euGen

# Step 3: Prepare Year data
country_resistance_carbap_eu$Year <- as.numeric(as.character(country_resistance_carbap_eu$Year))
years_range <- seq(min(country_resistance_carbap_eu$Year), max(country_resistance_carbap_eu$Year), length.out = 100)

# Step 4: Prepare spatial data
states_spatial_filtered <- model_output_carbap_eu$states_spatial_filtered
states_data_df <- as.data.frame(states_spatial_filtered@data)
original_geo_levels <- unique(states_data_df$GEOID)

# Step 5: Prepare unique GEOID data
GEOID_nameC <- states_data_df %>%
  dplyr::select(GEOID, NAME) %>%
  distinct(GEOID, .keep_all = TRUE)

# Step 6: Expand prediction data including Gender
# Assuming Gender levels are "Male" and "Female"
gender_levels <- c("Male", "Female")

pred_data <- expand.grid(Year = seq(min(years_range), max(years_range), by = 1),
                         GEOID = original_geo_levels,
                         Gender = gender_levels)  # Include Gender in prediction

# Step 7: Predict from GAM using the created function
predictions <- gam_predictions(model_output_carbap_eu$fr, newdata = pred_data)
predictions <- merge(predictions, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions <- predictions %>% arrange(GEOID, Gender, Year)



# Step 8: Repeat for carbapenem resistance predictions
predictions_carb_EU <- gam_predictions(model_output_carbap_eu$fr, newdata = pred_data)
predictions_carb_EU <- merge(predictions_carb_EU, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions_carb_EU <- predictions_carb_EU %>% arrange(GEOID, Gender, Year)

# Step 9: Calculate growth derivatives
resultsGrowth <- derivatives_mh2(model_output_carbap_eu$fr, predictions_carb_EU)


library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Step 1: Calculate derivatives by Gender
derivatives_data <- derivatives_mh(model_output_carbap_eu$fr, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001, startpoint = 0)

# Step 2: Calculate growth rates
derivatives_data <- derivatives_data %>%
  mutate(growth_rate = first_derivative + (second_derivative / first_derivative))

# Step 3: Merge with GEOID data
merged_data <- merge(derivatives_data, unique_GEOID_data, by = "GEOID", all.x = TRUE)
merged_data$first_derivative_sign_change[merged_data$Year == 2004] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2004] <- 0
merged_data$first_derivative_sign_change[merged_data$Year == 2005] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2005] <- 0
# Step 4: Process predictions by Gender
predictions <- predictions %>%
  group_by(GEOID, Gender) %>%  # Grouping by Gender
  mutate(
    pred_lag = lag(pred, default = NA),
    pred_lagu = lag(pred_upper, default = NA),
    pred_lagl = lag(pred_lower, default = NA),
    growth_rate2 = if_else(is.na(pred_lag), NA_real_, 100 * (pred - pred_lag) / pred_lag),
    growth_rate2_up = if_else(is.na(pred_lagu), NA_real_, 100 * (pred_upper - pred_lagu) / pred_lagu),
    growth_rate2_lo = if_else(is.na(pred_lagl), NA_real_, 100 * (pred_lower - pred_lagl) / pred_lagl)
  ) %>%
  dplyr::select(-pred_lag, -pred_lagu, -pred_lagl)  # Clean up lag columns

# Step 5: Store results
Carb_predictions_grat <- predictions
Carb_changep_EU <- merged_data

# Step 6: Plotting by Gender
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Calculate y-axis limits
y_max <- max(merged_data$first_upper, na.rm = TRUE)
y_min <- min(merged_data$first_lower, na.rm = TRUE)
padding <- (y_max - y_min) * 0.05
y_max <- y_max + padding
y_min <- y_min - padding

# Plotting with facets for Gender
library(ggplot2)
library(dplyr)

# Step 1: Filter the merged data by Gender
merged_data_male <- merged_data %>% filter(Gender == "Male")
merged_data_female <- merged_data %>% filter(Gender == "Female")

# Step 2: Plot for Males
p_male <- ggplot(data = merged_data_male, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_male, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Males",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 3: Plot for Females
p_female <- ggplot(data = merged_data_female, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Females",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Save the plots
#ggsave(filename = "first_derivat_eu_carb_male.tiff", plot = p_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "first_derivat_eu_carb_female.tiff", plot = p_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Step 1: Filter merged data by Gender
merged_data_male <- merged_data %>% filter(Gender == "Male")
merged_data_female <- merged_data %>% filter(Gender == "Female")

# Step 2: Define color palette
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 3: Plot for Males - Second Derivative
p2_male2 <- ggplot(data = merged_data_male, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_male, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Males",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Plot for Females - Second Derivative
p2_female2 <- ggplot(data = merged_data_female, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Females",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 5: Save the plots
#ggsave(filename = "second_derivat_eu_carb_male.tiff", plot = p2_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "second_derivat_eu_carb_female.tiff", plot = p2_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")



###REVIEW GROWTH RATES, I presume logs are not calculated and bringing NA values due to 0-values.


library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Step 1: Calculate Growth Rate with Gender
# Step 2: Merge with Predictions and GEOID Data (including Gender)
growth_rate_ci <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE) %>%
  arrange(GEOID, Year)

# Step 3: Calculate Doubling and Halving Times
growth_rate_ci <- growth_rate_ci %>%
  mutate(
    doubling_times = log(2) / growth_rate2,
    halving_times = log(0.5) / growth_rate2,
    NAME = NAME.x,  # Ensure proper naming
    growth_rate2 = if_else(Year == 2004, NA_real_, growth_rate2),
    growth_rate2_lo = if_else(Year == 2004, NA_real_, growth_rate2_lo),
    growth_rate2_up = if_else(Year == 2004, NA_real_, growth_rate2_up),
    doubling_times = if_else(Year == 2004, NA_real_, doubling_times),
    halving_times = if_else(Year == 2004, NA_real_, halving_times)
  )

# Step 4: Filter Data by Gender
merged_data2 <- growth_rate_ci %>%
  dplyr::left_join(merged_data %>% dplyr::select(GEOID, Year, Gender, first_derivative_sign_change, derivative_breakpoint),
                   by = c("GEOID", "Year", "Gender"))

growth_rate_male <- merged_data2 %>% filter(Gender == "Male")
growth_rate_female <- merged_data2 %>% filter(Gender == "Female")


# Step 5: Define Color Palette
num_colors <- length(unique(growth_rate_ci$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 6: Plot for Males
p3_male <- ggplot(data = growth_rate_male, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_male , derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Males",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 7: Plot for Females
p3_female <- ggplot(data = growth_rate_female, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Females",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 8: Save the Plots
ggsave(filename = "growth_rate_carbEU_male.tiff", plot = p3_male, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_carbEU_female.tiff", plot = p3_female, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")





max_y<-40

ppred_male <- ggplot(data = filter(predictions, Gender == "Male"), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Gender == "Male"), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Gender == "Male"), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 10)) +
  labs(title = "Predicted Carbapenem-Resistance (%) - Male",
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )

# Plot for Female
ppred_female <- ggplot(data = filter(predictions, Gender == "Female"), 
                       aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Gender == "Female"), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Gender == "Female"), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 10)) +
  labs(title = "Predicted Carbapenem-Resistance (%) - Female",
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
# Display the plot
print(ppred_female)
ggsave(filename = "predictions_breakpoint_carb_female.tiff", plot = ppred_female, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_carb_male.tiff", plot = ppred_male, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")



######
#------------------------------------------------------------------------------#
#NEW, CEPHALOSPORIN RESISTANCE GRAPH: ######
# Step 1: Prepare filtered dataset with Gender trimmed
country_resistance_cephalos_eu <- country_resistance_cephalos_euGen %>%
  mutate(Gender = str_trim(Gender)) %>%  # Trim whitespace
  filter(Gender == "Male")  # Filter for Male

# Step 2: Set up the model output
model_output_cephalos_eu <- model_output_cephalos_euGen

# Step 3: Prepare Year data
country_resistance_cephalos_eu$Year <- as.numeric(as.character(country_resistance_cephalos_eu$Year))
years_range <- seq(min(country_resistance_cephalos_eu$Year), max(country_resistance_cephalos_eu$Year), length.out = 100)

# Step 4: Prepare spatial data
states_spatial_filtered <- model_output_cephalos_eu$states_spatial_filtered
states_data_df <- as.data.frame(states_spatial_filtered@data)
original_geo_levels <- unique(states_data_df$GEOID)

# Step 5: Prepare unique GEOID data
GEOID_nameC <- states_data_df %>%
  dplyr::select(GEOID, NAME) %>%
  distinct(GEOID, .keep_all = TRUE)

# Step 6: Expand prediction data including Gender
# Assuming Gender levels are "Male" and "Female"
gender_levels <- c("Male", "Female")

pred_data <- expand.grid(Year = seq(min(years_range), max(years_range), by = 1),
                         GEOID = original_geo_levels,
                         Gender = gender_levels)  # Include Gender in prediction

# Step 7: Predict from GAM using the created function
predictions <- gam_predictions(model_output_cephalos_eu$fr, newdata = pred_data)
predictions <- merge(predictions, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions <- predictions %>% arrange(GEOID, Gender, Year)



# Step 8: Repeat for carbapenem resistance predictions
predictions_cephalos_EU <- gam_predictions(model_output_cephalos_eu$fr, newdata = pred_data)
predictions_cephalos_EU <- merge(predictions_cephalos_EU, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions_cephalos_EU <- predictions_cephalos_EU %>% arrange(GEOID, Gender, Year)

# Step 9: Calculate growth derivatives
resultsGrowth <- derivatives_mh2(model_output_cephalos_eu$fr, predictions_cephalos_EU)


library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Step 1: Calculate derivatives by Gender
derivatives_data <- derivatives_mh(model_output_cephalos_eu$fr, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001, startpoint = 0)

# Step 2: Calculate growth rates
derivatives_data <- derivatives_data %>%
  mutate(growth_rate = first_derivative + (second_derivative / first_derivative))

# Step 3: Merge with GEOID data
merged_data <- merge(derivatives_data, unique_GEOID_data, by = "GEOID", all.x = TRUE)
merged_data$first_derivative_sign_change[merged_data$Year == 2004] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2004] <- 0
merged_data$first_derivative_sign_change[merged_data$Year == 2005] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2005] <- 0

# Step 4: Process predictions by Gender
predictions <- predictions %>%
  group_by(GEOID, Gender) %>%  # Grouping by Gender
  mutate(
    pred_lag = lag(pred, default = NA),
    pred_lagu = lag(pred_upper, default = NA),
    pred_lagl = lag(pred_lower, default = NA),
    growth_rate2 = if_else(is.na(pred_lag), NA_real_, 100 * (pred - pred_lag) / pred_lag),
    growth_rate2_up = if_else(is.na(pred_lagu), NA_real_, 100 * (pred_upper - pred_lagu) / pred_lagu),
    growth_rate2_lo = if_else(is.na(pred_lagl), NA_real_, 100 * (pred_lower - pred_lagl) / pred_lagl)
  ) %>%
  dplyr::select(-pred_lag, -pred_lagu, -pred_lagl)  # Clean up lag columns

# Step 5: Store results
Carb_predictions_grat <- predictions
Carb_changep_EU <- merged_data

# Step 6: Plotting by Gender
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Calculate y-axis limits
y_max <- max(merged_data$first_upper, na.rm = TRUE)
y_min <- min(merged_data$first_lower, na.rm = TRUE)
padding <- (y_max - y_min) * 0.05
y_max <- y_max + padding
y_min <- y_min - padding

# Plotting with facets for Gender
library(ggplot2)
library(dplyr)

# Step 1: Filter the merged data by Gender
merged_data_male <- merged_data %>% filter(Gender == "Male")
merged_data_female <- merged_data %>% filter(Gender == "Female")

# Step 2: Plot for Males
p_male <- ggplot(data = merged_data_male, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_male, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Males",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 3: Plot for Females
p_female <- ggplot(data = merged_data_female, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Females",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Save the plots
#ggsave(filename = "first_derivat_eu_carb_male.tiff", plot = p_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "first_derivat_eu_carb_female.tiff", plot = p_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Step 1: Filter merged data by Gender
merged_data_male <- merged_data %>% filter(Gender == "Male")
merged_data_female <- merged_data %>% filter(Gender == "Female")

# Step 2: Define color palette
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 3: Plot for Males - Second Derivative
p2_male2 <- ggplot(data = merged_data_male, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_male, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Males",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Plot for Females - Second Derivative
p2_female2 <- ggplot(data = merged_data_female, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Females",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 5: Save the plots
#ggsave(filename = "second_derivat_eu_carb_male.tiff", plot = p2_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "second_derivat_eu_carb_female.tiff", plot = p2_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")



###REVIEW GROWTH RATES, I presume logs are not calculated and bringing NA values due to 0-values.


library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Step 1: Calculate Growth Rate with Gender
# Step 2: Merge with Predictions and GEOID Data (including Gender)
growth_rate_ci <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE) %>%
  arrange(GEOID, Year)

# Step 3: Calculate Doubling and Halving Times
growth_rate_ci <- growth_rate_ci %>%
  mutate(
    doubling_times = log(2) / growth_rate2,
    halving_times = log(0.5) / growth_rate2,
    NAME = NAME.x,  # Ensure proper naming
    growth_rate2 = if_else(Year == 2004, NA_real_, growth_rate2),
    growth_rate2_lo = if_else(Year == 2004, NA_real_, growth_rate2_lo),
    growth_rate2_up = if_else(Year == 2004, NA_real_, growth_rate2_up),
    doubling_times = if_else(Year == 2004, NA_real_, doubling_times),
    halving_times = if_else(Year == 2004, NA_real_, halving_times)
  )

# Step 4: Filter Data by Gender
merged_data2 <- growth_rate_ci %>%
  dplyr::left_join(merged_data %>% dplyr::select(GEOID, Year, Gender, first_derivative_sign_change, derivative_breakpoint),
                   by = c("GEOID", "Year", "Gender"))

growth_rate_male <- merged_data2 %>% filter(Gender == "Male")
growth_rate_female <- merged_data2 %>% filter(Gender == "Female")


# Step 5: Define Color Palette
num_colors <- length(unique(growth_rate_ci$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 6: Plot for Males
p3_male <- ggplot(data = growth_rate_male, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_male , derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Males",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 7: Plot for Females
p3_female <- ggplot(data = growth_rate_female, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Females",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 8: Save the Plots
ggsave(filename = "growth_rate_3GCREU_male.tiff", plot = p3_male, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_3GCREU_female.tiff", plot = p3_female, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")





max_y<-70

ppred_male <- ggplot(data = filter(predictions, Gender == "Male"), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Gender == "Male"), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Gender == "Male"), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 10)) +
  labs(title = "Predicted 3GCR (%) - Male",
       x = "Year",
       y = "Predicted 3GCR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )

# Plot for Female
ppred_female <- ggplot(data = filter(predictions, Gender == "Female"), 
                       aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Gender == "Female"), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Gender == "Female"), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 10)) +
  labs(title = "Predicted 3GCR (%) - Female",
       x = "Year",
       y = "Predicted 3GCR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
# Display the plot
print(ppred_female)
ggsave(filename = "predictions_breakpoint_3GCR_female.tiff", plot = ppred_female, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_3GCR_male.tiff", plot = ppred_male, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


######
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#NEW, MDR RESISTANCE GRAPH: ######
# Step 1: Prepare filtered dataset with Gender trimmed
country_resistance_mdr_eu <- country_resistance_mdr_euGen %>%
  mutate(Gender = str_trim(Gender)) %>%  # Trim whitespace
  filter(Gender == "Male")  # Filter for Male

# Step 2: Set up the model output
model_output_mdr_eu <- model_output_mdr_euGen

# Step 3: Prepare Year data
country_resistance_mdr_eu$Year <- as.numeric(as.character(country_resistance_mdr_eu$Year))
years_range <- seq(min(country_resistance_mdr_eu$Year), max(country_resistance_mdr_eu$Year), length.out = 100)

# Step 4: Prepare spatial data
states_spatial_filtered <- model_output_mdr_eu$states_spatial_filtered
states_data_df <- as.data.frame(states_spatial_filtered@data)
original_geo_levels <- unique(states_data_df$GEOID)

# Step 5: Prepare unique GEOID data
GEOID_nameC <- states_data_df %>%
  dplyr::select(GEOID, NAME) %>%
  distinct(GEOID, .keep_all = TRUE)

# Step 6: Expand prediction data including Gender
# Assuming Gender levels are "Male" and "Female"
gender_levels <- c("Male", "Female")

pred_data <- expand.grid(Year = seq(min(years_range), max(years_range), by = 1),
                         GEOID = original_geo_levels,
                         Gender = gender_levels)  # Include Gender in prediction

# Step 7: Predict from GAM using the created function
predictions <- gam_predictions(model_output_mdr_eu$sf, newdata = pred_data)
predictions <- merge(predictions, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions <- predictions %>% arrange(GEOID, Gender, Year)



# Step 8: Repeat for carbapenem resistance predictions
predictions_mdr_EU <- gam_predictions(model_output_mdr_eu$sf, newdata = pred_data)
predictions_mdr_EU <- merge(predictions_mdr_EU, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions_mdr_EU <- predictions_mdr_EU %>% arrange(GEOID, Gender, Year)

# Step 9: Calculate growth derivatives
resultsGrowth <- derivatives_mh2(model_output_mdr_eu$sf, predictions_mdr_EU)


library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Step 1: Calculate derivatives by Gender
derivatives_data <- derivatives_mh(model_output_mdr_eu$sf, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001, startpoint = 0)

# Step 2: Calculate growth rates
derivatives_data <- derivatives_data %>%
  mutate(growth_rate = first_derivative + (second_derivative / first_derivative))

# Step 3: Merge with GEOID data
merged_data <- merge(derivatives_data, unique_GEOID_data, by = "GEOID", all.x = TRUE)
merged_data$first_derivative_sign_change[merged_data$Year == 2004] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2004] <- 0
merged_data$first_derivative_sign_change[merged_data$Year == 2005] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2005] <- 0

# Step 4: Process predictions by Gender
predictions <- predictions %>%
  group_by(GEOID, Gender) %>%  # Grouping by Gender
  mutate(
    pred_lag = lag(pred, default = NA),
    pred_lagu = lag(pred_upper, default = NA),
    pred_lagl = lag(pred_lower, default = NA),
    growth_rate2 = if_else(is.na(pred_lag), NA_real_, 100 * (pred - pred_lag) / pred_lag),
    growth_rate2_up = if_else(is.na(pred_lagu), NA_real_, 100 * (pred_upper - pred_lagu) / pred_lagu),
    growth_rate2_lo = if_else(is.na(pred_lagl), NA_real_, 100 * (pred_lower - pred_lagl) / pred_lagl)
  ) %>%
  dplyr::select(-pred_lag, -pred_lagu, -pred_lagl)  # Clean up lag columns

# Step 5: Store results
Carb_predictions_grat <- predictions
Carb_changep_EU <- merged_data

# Step 6: Plotting by Gender
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Calculate y-axis limits
y_max <- max(merged_data$first_upper, na.rm = TRUE)
y_min <- min(merged_data$first_lower, na.rm = TRUE)
padding <- (y_max - y_min) * 0.05
y_max <- y_max + padding
y_min <- y_min - padding

# Plotting with facets for Gender
library(ggplot2)
library(dplyr)

# Step 1: Filter the merged data by Gender
merged_data_male <- merged_data %>% filter(Gender == "Male")
merged_data_female <- merged_data %>% filter(Gender == "Female")

# Step 2: Plot for Males
p_male <- ggplot(data = merged_data_male, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_male, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Males",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 3: Plot for Females
p_female <- ggplot(data = merged_data_female, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Females",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Save the plots
#ggsave(filename = "first_derivat_eu_carb_male.tiff", plot = p_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "first_derivat_eu_carb_female.tiff", plot = p_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(RColorBrewer)
merged_data$first_derivative_sign_change[merged_data$Year == 2004] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2004] <- 0


# Step 1: Filter merged data by Gender
merged_data_male <- merged_data %>% filter(Gender == "Male")
merged_data_female <- merged_data %>% filter(Gender == "Female")

# Step 2: Define color palette
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 3: Plot for Males - Second Derivative
p2_male2 <- ggplot(data = merged_data_male, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_male, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Males",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Plot for Females - Second Derivative
p2_female2 <- ggplot(data = merged_data_female, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Females",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 5: Save the plots
#ggsave(filename = "second_derivat_eu_carb_male.tiff", plot = p2_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "second_derivat_eu_carb_female.tiff", plot = p2_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")



###REVIEW GROWTH RATES, I presume logs are not calculated and bringing NA values due to 0-values.


library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Step 1: Calculate Growth Rate with Gender
# Step 2: Merge with Predictions and GEOID Data (including Gender)
growth_rate_ci <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE) %>%
  arrange(GEOID, Year)

# Step 3: Calculate Doubling and Halving Times
growth_rate_ci <- growth_rate_ci %>%
  mutate(
    doubling_times = log(2) / growth_rate2,
    halving_times = log(0.5) / growth_rate2,
    NAME = NAME.x,  # Ensure proper naming
    growth_rate2 = if_else(Year == 2004, NA_real_, growth_rate2),
    growth_rate2_lo = if_else(Year == 2004, NA_real_, growth_rate2_lo),
    growth_rate2_up = if_else(Year == 2004, NA_real_, growth_rate2_up),
    doubling_times = if_else(Year == 2004, NA_real_, doubling_times),
    halving_times = if_else(Year == 2004, NA_real_, halving_times)
  )

# Step 4: Filter Data by Gender
merged_data2 <- growth_rate_ci %>%
  dplyr::left_join(merged_data %>% dplyr::select(GEOID, Year, Gender, first_derivative_sign_change, derivative_breakpoint),
                   by = c("GEOID", "Year", "Gender"))

growth_rate_male <- merged_data2 %>% filter(Gender == "Male")
growth_rate_female <- merged_data2 %>% filter(Gender == "Female")


# Step 5: Define Color Palette
num_colors <- length(unique(growth_rate_ci$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 6: Plot for Males
p3_male <- ggplot(data = growth_rate_male, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_male , derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Males",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 7: Plot for Females
p3_female <- ggplot(data = growth_rate_female, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Females",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 8: Save the Plots
ggsave(filename = "growth_rate_MDR_male.tiff", plot = p3_male, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_MDR_female.tiff", plot = p3_female, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")





merged_data$first_derivative_sign_change[merged_data$Year == 2004] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2004] <- 0
max_y<-80

ppred_male <- ggplot(data = filter(predictions, Gender == "Male"), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Gender == "Male"), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Gender == "Male"), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted MDR (%) - Male",
       x = "Year",
       y = "Predicted MDR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )

# Plot for Female
ppred_female <- ggplot(data = filter(predictions, Gender == "Female"), 
                       aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Gender == "Female"), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Gender == "Female"), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted MDR (%) - Female",
       x = "Year",
       y = "Predicted MDR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
# Display the plot
print(ppred_female)
ggsave(filename = "predictions_breakpoint_MDR_female.tiff", plot = ppred_female, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_MDR_male.tiff", plot = ppred_male, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


######
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#Subgroup analyses - Age/group-specific change points [Europe]
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#NEW, CARBAPENEM RESISTANCE GRAPH: ######
# Step 1: Prepare filtered dataset with Gender trimmed
country_resistance_carbap_eu <- country_resistance_carbap_euAgeg

# Step 2: Set up the model output
model_output_carbap_eu <- model_output_carbap_euAgeg

# Step 3: Prepare Year data
country_resistance_carbap_eu$Year <- as.numeric(as.character(country_resistance_carbap_eu$Year))
years_range <- seq(min(country_resistance_carbap_eu$Year), max(country_resistance_carbap_eu$Year), length.out = 100)

# Step 4: Prepare spatial data
states_spatial_filtered <- model_output_carbap_eu$states_spatial_filtered
states_data_df <- as.data.frame(states_spatial_filtered@data)
original_geo_levels <- unique(states_data_df$GEOID)

# Step 5: Prepare unique GEOID data
GEOID_nameC <- states_data_df %>%
  dplyr::select(GEOID, NAME) %>%
  distinct(GEOID, .keep_all = TRUE)

# Step 6: Expand prediction data including Agegroup
Ageg_levels <- c(0, 1, 2)

pred_data <- expand.grid(Year = seq(min(years_range), max(years_range), by = 1),
                         GEOID = original_geo_levels,
                         Agegroup = Ageg_levels)  # Include Gender in prediction

# Step 7: Predict from GAM using the created function
predictions <- gam_predictions(model_output_carbap_eu$fr, newdata = pred_data)
predictions <- merge(predictions, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions <- predictions %>% arrange(GEOID, Agegroup, Year)



# Step 8: Repeat for carbapenem resistance predictions
predictions_carb_EU <- gam_predictions(model_output_carbap_eu$fr, newdata = pred_data)
predictions_carb_EU <- merge(predictions_carb_EU, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions_carb_EU <- predictions_carb_EU %>% arrange(GEOID, Agegroup, Year)

# Step 9: Calculate growth derivatives
resultsGrowth <- derivatives_mh2(model_output_carbap_eu$fr, predictions_carb_EU)


library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Step 1: Calculate derivatives by Gender
derivatives_data <- derivatives_mh(model_output_carbap_eu$fr, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001, startpoint = 0)

# Step 2: Calculate growth rates
derivatives_data <- derivatives_data %>%
  mutate(growth_rate = first_derivative + (second_derivative / first_derivative))

# Step 3: Merge with GEOID data
merged_data <- merge(derivatives_data, unique_GEOID_data, by = "GEOID", all.x = TRUE)
merged_data$first_derivative_sign_change[merged_data$Year == 2004] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2004] <- 0
merged_data$first_derivative_sign_change[merged_data$Year == 2005] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2005] <- 0
# Step 4: Process predictions by Gender
predictions <- predictions %>%
  group_by(GEOID, Agegroup) %>%  # Grouping by Gender
  mutate(
    pred_lag = lag(pred, default = NA),
    pred_lagu = lag(pred_upper, default = NA),
    pred_lagl = lag(pred_lower, default = NA),
    growth_rate2 = if_else(is.na(pred_lag), NA_real_, 100 * (pred - pred_lag) / pred_lag),
    growth_rate2_up = if_else(is.na(pred_lagu), NA_real_, 100 * (pred_upper - pred_lagu) / pred_lagu),
    growth_rate2_lo = if_else(is.na(pred_lagl), NA_real_, 100 * (pred_lower - pred_lagl) / pred_lagl)
  ) %>%
  dplyr::select(-pred_lag, -pred_lagu, -pred_lagl)  # Clean up lag columns

# Step 5: Store results
Carb_predictions_grat <- predictions
Carb_changep_EU <- merged_data

# Step 6: Plotting by Agegroup
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Calculate y-axis limits
y_max <- max(merged_data$first_upper, na.rm = TRUE)
y_min <- min(merged_data$first_lower, na.rm = TRUE)
padding <- (y_max - y_min) * 0.05
y_max <- y_max + padding
y_min <- y_min - padding

# Plotting with facets for Gender
library(ggplot2)
library(dplyr)

# Step 1: Filter the merged data by Agegroup labels = c("≤18yo", "19≤ and ≤64", "≥65")
merged_data_age0 <- merged_data %>% filter(Agegroup == 0)
merged_data_age1 <- merged_data %>% filter(Agegroup == 1)
merged_data_age2 <- merged_data %>% filter(Agegroup == 2)

# Step 2: Plot for Age0
p_age0 <- ggplot(data = merged_data_age0, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age0, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age ≤18yo",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 3: Plot for Females
p_age1<- ggplot(data = merged_data_age1, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age 19≤ and ≤64",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

p_age2<- ggplot(data = merged_data_age2, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age ≥65",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")




# Step 4: Save the plots
#ggsave(filename = "first_derivat_eu_carb_male.tiff", plot = p_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "first_derivat_eu_carb_female.tiff", plot = p_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Step 1: Filter merged data by Gender
merged_data_age0 <- merged_data %>% filter(Agegroup == 0)
merged_data_age1 <- merged_data %>% filter(Agegroup == 1)
merged_data_age2 <- merged_data %>% filter(Agegroup == 2)

# Step 2: Define color palette
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 3: Plot for Males - Second Derivative
p2_age0_2 <- ggplot(data = merged_data_age0, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age0, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age ≤18yo",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Plot for Females - Second Derivative
p2_age1_2 <- ggplot(data = merged_data_age1, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age 19≤ and ≤64",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

p2_age2_2 <- ggplot(data = merged_data_age2, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age ≥65",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")


# Step 5: Save the plots
#ggsave(filename = "second_derivat_eu_carb_male.tiff", plot = p2_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "second_derivat_eu_carb_female.tiff", plot = p2_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Step 1: Calculate Growth Rate with agegroups
# Step 2: Merge with Predictions and GEOID Data (including Agegroups)
growth_rate_ci <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE) %>%
  arrange(GEOID, Year)

# Step 3: Calculate Doubling and Halving Times
growth_rate_ci <- growth_rate_ci %>%
  mutate(
    doubling_times = log(2) / growth_rate2,
    halving_times = log(0.5) / growth_rate2,
    NAME = NAME.x,  # Ensure proper naming
    growth_rate2 = if_else(Year == 2004, NA_real_, growth_rate2),
    growth_rate2_lo = if_else(Year == 2004, NA_real_, growth_rate2_lo),
    growth_rate2_up = if_else(Year == 2004, NA_real_, growth_rate2_up),
    doubling_times = if_else(Year == 2004, NA_real_, doubling_times),
    halving_times = if_else(Year == 2004, NA_real_, halving_times)
  )

# Step 4: Filter Data by Gender
merged_data2 <- growth_rate_ci %>%
  dplyr::left_join(merged_data %>% dplyr::select(GEOID, Year, Agegroup, first_derivative_sign_change, derivative_breakpoint),
                   by = c("GEOID", "Year", "Agegroup"))

growth_rate_age0 <- merged_data2 %>% filter(Agegroup == 0)
growth_rate_age1 <- merged_data2 %>% filter(Agegroup == 1)
growth_rate_age2 <- merged_data2 %>% filter(Agegroup == 2)

#≤18yo", "19≤ and ≤64", "≥65"
# Step 5: Define Color Palette
num_colors <- length(unique(growth_rate_ci$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 6: Plot for Males
p3_age0 <- ggplot(data = growth_rate_age0, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age0 , derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age ≤18yo",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 7: Plot for Females
p3_age1 <- ggplot(data = growth_rate_age1, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age 19≤ and ≤64",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

p3_age2 <- ggplot(data = growth_rate_age2, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age ≥65",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 8: Save the Plots
ggsave(filename = "growth_rate_carbEU_age0.tiff", plot = p3_age0, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_carbEU_age1.tiff", plot = p3_age1, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_carbEU_age2.tiff", plot = p3_age2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")



max_y<-40

ppred_age0 <- ggplot(data = filter(predictions, Agegroup == 0), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup==0), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup==0), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 10)) +
  labs(title = "Predicted Carbapenem-Resistance (%) - Age ≤18yo",
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
#≤18yo", "19≤ and ≤64", "≥65"

# Plot for Female
ppred_age1 <- ggplot(data = filter(predictions, Agegroup == 1), 
                       aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 10)) +
  labs(title = "Predicted Carbapenem-Resistance (%) - Age 19≤ and ≤64",
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )

ppred_age2 <- ggplot(data = filter(predictions, Agegroup == 2), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup == 2), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup == 2), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 10)) +
  labs(title = "Predicted Carbapenem-Resistance (%) - Age ≥65",
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
# Display the plot
print(ppred_female)
ggsave(filename = "predictions_breakpoint_carb_age0.tiff", plot = ppred_age0, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_carb_age1.tiff", plot = ppred_age1, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_carb_age2.tiff", plot = ppred_age2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


######
#------------------------------------------------------------------------------#
#NEW, CEPHALOSPORIN RESISTANCE GRAPH: ######
# Step 1: Prepare filtered dataset with Gender trimmed
country_resistance_cephalos_eu <- country_resistance_cephalos_euAgeg

# Step 2: Set up the model output
model_output_cephalos_eu <- model_output_cephalos_euAgeg

# Step 3: Prepare Year data
country_resistance_cephalos_eu$Year <- as.numeric(as.character(country_resistance_cephalos_eu$Year))
years_range <- seq(min(country_resistance_cephalos_eu$Year), max(country_resistance_cephalos_eu$Year), length.out = 100)

# Step 4: Prepare spatial data
states_spatial_filtered <- model_output_cephalos_eu$states_spatial_filtered
states_data_df <- as.data.frame(states_spatial_filtered@data)
original_geo_levels <- unique(states_data_df$GEOID)

# Step 5: Prepare unique GEOID data
GEOID_nameC <- states_data_df %>%
  dplyr::select(GEOID, NAME) %>%
  distinct(GEOID, .keep_all = TRUE)

# Step 6: Expand prediction data including Agegroup
Ageg_levels <- c(0, 1, 2)

pred_data <- expand.grid(Year = seq(min(years_range), max(years_range), by = 1),
                         GEOID = original_geo_levels,
                         Agegroup = Ageg_levels)  # Include Gender in prediction

# Step 7: Predict from GAM using the created function
predictions <- gam_predictions(model_output_cephalos_eu$fr, newdata = pred_data)
predictions <- merge(predictions, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions <- predictions %>% arrange(GEOID, Agegroup, Year)



# Step 8: Repeat for carbapenem resistance predictions
predictions_cephalos_EU <- gam_predictions(model_output_cephalos_eu$fr, newdata = pred_data)
predictions_cephalos_EU <- merge(predictions_cephalos_EU, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions_cephalos_EU <- predictions_cephalos_EU %>% arrange(GEOID, Agegroup, Year)

# Step 9: Calculate growth derivatives
resultsGrowth <- derivatives_mh2(model_output_cephalos_eu$fr, predictions_cephalos_EU)


library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Step 1: Calculate derivatives by Gender
derivatives_data <- derivatives_mh(model_output_cephalos_eu$fr, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001, startpoint = 0)

# Step 2: Calculate growth rates
derivatives_data <- derivatives_data %>%
  mutate(growth_rate = first_derivative + (second_derivative / first_derivative))

# Step 3: Merge with GEOID data
merged_data <- merge(derivatives_data, unique_GEOID_data, by = "GEOID", all.x = TRUE)
merged_data$first_derivative_sign_change[merged_data$Year == 2004] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2004] <- 0
merged_data$first_derivative_sign_change[merged_data$Year == 2005] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2005] <- 0
# Step 4: Process predictions by Gender
predictions <- predictions %>%
  group_by(GEOID, Agegroup) %>%  # Grouping by Gender
  mutate(
    pred_lag = lag(pred, default = NA),
    pred_lagu = lag(pred_upper, default = NA),
    pred_lagl = lag(pred_lower, default = NA),
    growth_rate2 = if_else(is.na(pred_lag), NA_real_, 100 * (pred - pred_lag) / pred_lag),
    growth_rate2_up = if_else(is.na(pred_lagu), NA_real_, 100 * (pred_upper - pred_lagu) / pred_lagu),
    growth_rate2_lo = if_else(is.na(pred_lagl), NA_real_, 100 * (pred_lower - pred_lagl) / pred_lagl)
  ) %>%
  dplyr::select(-pred_lag, -pred_lagu, -pred_lagl)  # Clean up lag columns

# Step 5: Store results
Carb_predictions_grat <- predictions
Carb_changep_EU <- merged_data

# Step 6: Plotting by Agegroup
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Calculate y-axis limits
y_max <- max(merged_data$first_upper, na.rm = TRUE)
y_min <- min(merged_data$first_lower, na.rm = TRUE)
padding <- (y_max - y_min) * 0.05
y_max <- y_max + padding
y_min <- y_min - padding

# Plotting with facets for Gender
library(ggplot2)
library(dplyr)

# Step 1: Filter the merged data by Agegroup labels = c("≤18yo", "19≤ and ≤64", "≥65")
merged_data_age0 <- merged_data %>% filter(Agegroup == 0)
merged_data_age1 <- merged_data %>% filter(Agegroup == 1)
merged_data_age2 <- merged_data %>% filter(Agegroup == 2)

# Step 2: Plot for age0
p_age0 <- ggplot(data = merged_data_age0, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age0, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age ≤18yo",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 3: Plot for age1
p_age1<- ggplot(data = merged_data_age1, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age 19≤ and ≤64",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

p_age2<- ggplot(data = merged_data_age2, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age ≥65",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")


# Step 4: Save the plots
#ggsave(filename = "first_derivat_eu_carb_male.tiff", plot = p_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "first_derivat_eu_carb_female.tiff", plot = p_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Step 1: Filter merged data by Gender
merged_data_age0 <- merged_data %>% filter(Agegroup == 0)
merged_data_age1 <- merged_data %>% filter(Agegroup == 1)
merged_data_age2 <- merged_data %>% filter(Agegroup == 2)

# Step 2: Define color palette
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 3: Plot for Males - Second Derivative
p2_age0_2 <- ggplot(data = merged_data_age0, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age0, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age ≤18yo",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Plot for Females - Second Derivative
p2_age1_2 <- ggplot(data = merged_data_age1, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age 19≤ and ≤64",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

p2_age2_2 <- ggplot(data = merged_data_age2, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age ≥65",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")


# Step 5: Save the plots
#ggsave(filename = "second_derivat_eu_carb_male.tiff", plot = p2_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "second_derivat_eu_carb_female.tiff", plot = p2_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Step 1: Calculate Growth Rate with agegroups
# Step 2: Merge with Predictions and GEOID Data (including Agegroups)
growth_rate_ci <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE) %>%
  arrange(GEOID, Year)

# Step 3: Calculate Doubling and Halving Times
growth_rate_ci <- growth_rate_ci %>%
  mutate(
    doubling_times = log(2) / growth_rate2,
    halving_times = log(0.5) / growth_rate2,
    NAME = NAME.x,  # Ensure proper naming
    growth_rate2 = if_else(Year == 2004, NA_real_, growth_rate2),
    growth_rate2_lo = if_else(Year == 2004, NA_real_, growth_rate2_lo),
    growth_rate2_up = if_else(Year == 2004, NA_real_, growth_rate2_up),
    doubling_times = if_else(Year == 2004, NA_real_, doubling_times),
    halving_times = if_else(Year == 2004, NA_real_, halving_times)
  )

# Step 4: Filter Data by Gender
merged_data2 <- growth_rate_ci %>%
  dplyr::left_join(merged_data %>% dplyr::select(GEOID, Year, Agegroup, first_derivative_sign_change, derivative_breakpoint),
                   by = c("GEOID", "Year", "Agegroup"))

growth_rate_age0 <- merged_data2 %>% filter(Agegroup == 0)
growth_rate_age1 <- merged_data2 %>% filter(Agegroup == 1)
growth_rate_age2 <- merged_data2 %>% filter(Agegroup == 2)

#≤18yo", "19≤ and ≤64", "≥65"
# Step 5: Define Color Palette
num_colors <- length(unique(growth_rate_ci$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 6: Plot for age0
p3_age0 <- ggplot(data = growth_rate_age0, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age0 , derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age ≤18yo",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 7: Plot for age1
p3_age1 <- ggplot(data = growth_rate_age1, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age 19≤ and ≤64",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

p3_age2 <- ggplot(data = growth_rate_age2, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age ≥65",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 8: Save the Plots
ggsave(filename = "growth_rate_cephalosEU_age0.tiff", plot = p3_age0, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_cephalosEU_age1.tiff", plot = p3_age1, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_cephalosEU_age2.tiff", plot = p3_age2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")



max_y<-80

ppred_age0 <- ggplot(data = filter(predictions, Agegroup == 0), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup==0), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup==0), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted 3GCR (%) - Age ≤18yo",
       x = "Year",
       y = "Predicted 3GCR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
#≤18yo", "19≤ and ≤64", "≥65"

# Plot for Female
ppred_age1 <- ggplot(data = filter(predictions, Agegroup == 1), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted 3GCR (%) - Age 19≤ and ≤64",
       x = "Year",
       y = "Predicted 3GCR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )

ppred_age2 <- ggplot(data = filter(predictions, Agegroup == 2), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup == 2), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup == 2), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted 3GCRe (%) - Age ≥65",
       x = "Year",
       y = "Predicted 3GCR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
# Display the plot
print(ppred_female)
ggsave(filename = "predictions_breakpoint_cephalos_age0.tiff", plot = ppred_age0, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_cephalos_age1.tiff", plot = ppred_age1, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_cephalos_age2.tiff", plot = ppred_age2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


######
#------------------------------------------------------------------------------#
#NEW, MDR RESISTANCE GRAPH: ######
# Step 1: Prepare filtered dataset with Gender trimmed
# Step 1: Prepare filtered dataset with Gender trimmed
country_resistance_mdr_eu <- country_resistance_mdr_euAgeg

# Step 2: Set up the model output
model_output_mdr_eu <- model_output_mdr_euAgeg

# Step 3: Prepare Year data
country_resistance_mdr_eu$Year <- as.numeric(as.character(country_resistance_mdr_eu$Year))
years_range <- seq(min(country_resistance_mdr_eu$Year), max(country_resistance_mdr_eu$Year), length.out = 100)

# Step 4: Prepare spatial data
states_spatial_filtered <- model_output_mdr_eu$states_spatial_filtered
states_data_df <- as.data.frame(states_spatial_filtered@data)
original_geo_levels <- unique(states_data_df$GEOID)

# Step 5: Prepare unique GEOID data
GEOID_nameC <- states_data_df %>%
  dplyr::select(GEOID, NAME) %>%
  distinct(GEOID, .keep_all = TRUE)

# Step 6: Expand prediction data including Agegroup
Ageg_levels <- c(0, 1, 2)

pred_data <- expand.grid(Year = seq(min(years_range), max(years_range), by = 1),
                         GEOID = original_geo_levels,
                         Agegroup = Ageg_levels)  # Include Gender in prediction

# Step 7: Predict from GAM using the created function
predictions <- gam_predictions(model_output_mdr_eu$sf, newdata = pred_data)
predictions <- merge(predictions, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions <- predictions %>% arrange(GEOID, Agegroup, Year)



# Step 8: Repeat for carbapenem resistance predictions
predictions_mdr_EU <- gam_predictions(model_output_mdr_eu$sf, newdata = pred_data)
predictions_mdr_EU <- merge(predictions_mdr_EU, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions_mdr_EU <- predictions_mdr_EU %>% arrange(GEOID, Agegroup, Year)

# Step 9: Calculate growth derivatives
resultsGrowth <- derivatives_mh2(model_output_mdr_eu$sf, predictions_mdr_EU)


library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Step 1: Calculate derivatives by Gender
derivatives_data <- derivatives_mh(model_output_mdr_eu$sf, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001, startpoint = 0)

# Step 2: Calculate growth rates
derivatives_data <- derivatives_data %>%
  mutate(growth_rate = first_derivative + (second_derivative / first_derivative))

# Step 3: Merge with GEOID data
merged_data <- merge(derivatives_data, unique_GEOID_data, by = "GEOID", all.x = TRUE)
merged_data$first_derivative_sign_change[merged_data$Year == 2004] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2004] <- 0
merged_data$first_derivative_sign_change[merged_data$Year == 2005] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2005] <- 0
# Step 4: Process predictions by Gender
predictions <- predictions %>%
  group_by(GEOID, Agegroup) %>%  # Grouping by Gender
  mutate(
    pred_lag = lag(pred, default = NA),
    pred_lagu = lag(pred_upper, default = NA),
    pred_lagl = lag(pred_lower, default = NA),
    growth_rate2 = if_else(is.na(pred_lag), NA_real_, 100 * (pred - pred_lag) / pred_lag),
    growth_rate2_up = if_else(is.na(pred_lagu), NA_real_, 100 * (pred_upper - pred_lagu) / pred_lagu),
    growth_rate2_lo = if_else(is.na(pred_lagl), NA_real_, 100 * (pred_lower - pred_lagl) / pred_lagl)
  ) %>%
  dplyr::select(-pred_lag, -pred_lagu, -pred_lagl)  # Clean up lag columns

# Step 5: Store results
Carb_predictions_grat <- predictions
Carb_changep_EU <- merged_data

# Step 6: Plotting by Agegroup
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Calculate y-axis limits
y_max <- max(merged_data$first_upper, na.rm = TRUE)
y_min <- min(merged_data$first_lower, na.rm = TRUE)
padding <- (y_max - y_min) * 0.05
y_max <- y_max + padding
y_min <- y_min - padding

# Plotting with facets for Gender
library(ggplot2)
library(dplyr)

# Step 1: Filter the merged data by Agegroup labels = c("≤18yo", "19≤ and ≤64", "≥65")
merged_data_age0 <- merged_data %>% filter(Agegroup == 0)
merged_data_age1 <- merged_data %>% filter(Agegroup == 1)
merged_data_age2 <- merged_data %>% filter(Agegroup == 2)

# Step 2: Plot for age0
p_age0 <- ggplot(data = merged_data_age0, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age0, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age ≤18yo",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 3: Plot for age1
p_age1<- ggplot(data = merged_data_age1, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age 19≤ and ≤64",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

p_age2<- ggplot(data = merged_data_age2, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age ≥65",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")


# Step 4: Save the plots
#ggsave(filename = "first_derivat_eu_carb_male.tiff", plot = p_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "first_derivat_eu_carb_female.tiff", plot = p_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Step 1: Filter merged data by Gender
merged_data_age0 <- merged_data %>% filter(Agegroup == 0)
merged_data_age1 <- merged_data %>% filter(Agegroup == 1)
merged_data_age2 <- merged_data %>% filter(Agegroup == 2)

# Step 2: Define color palette
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 3: Plot for Males - Second Derivative
p2_age0_2 <- ggplot(data = merged_data_age0, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age0, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age ≤18yo",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Plot for Females - Second Derivative
p2_age1_2 <- ggplot(data = merged_data_age1, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age 19≤ and ≤64",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

p2_age2_2 <- ggplot(data = merged_data_age2, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age ≥65",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")


# Step 5: Save the plots
#ggsave(filename = "second_derivat_eu_carb_male.tiff", plot = p2_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "second_derivat_eu_carb_female.tiff", plot = p2_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Step 1: Calculate Growth Rate with agegroups
# Step 2: Merge with Predictions and GEOID Data (including Agegroups)
growth_rate_ci <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE) %>%
  arrange(GEOID, Year)

# Step 3: Calculate Doubling and Halving Times
growth_rate_ci <- growth_rate_ci %>%
  mutate(
    doubling_times = log(2) / growth_rate2,
    halving_times = log(0.5) / growth_rate2,
    NAME = NAME.x,  # Ensure proper naming
    growth_rate2 = if_else(Year == 2004, NA_real_, growth_rate2),
    growth_rate2_lo = if_else(Year == 2004, NA_real_, growth_rate2_lo),
    growth_rate2_up = if_else(Year == 2004, NA_real_, growth_rate2_up),
    doubling_times = if_else(Year == 2004, NA_real_, doubling_times),
    halving_times = if_else(Year == 2004, NA_real_, halving_times)
  )

# Step 4: Filter Data by Gender
merged_data2 <- growth_rate_ci %>%
  dplyr::left_join(merged_data %>% dplyr::select(GEOID, Year, Agegroup, first_derivative_sign_change, derivative_breakpoint),
                   by = c("GEOID", "Year", "Agegroup"))

growth_rate_age0 <- merged_data2 %>% filter(Agegroup == 0)
growth_rate_age1 <- merged_data2 %>% filter(Agegroup == 1)
growth_rate_age2 <- merged_data2 %>% filter(Agegroup == 2)

#≤18yo", "19≤ and ≤64", "≥65"
# Step 5: Define Color Palette
num_colors <- length(unique(growth_rate_ci$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 6: Plot for age0
p3_age0 <- ggplot(data = growth_rate_age0, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age0 , derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age ≤18yo",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 7: Plot for age1
p3_age1 <- ggplot(data = growth_rate_age1, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age 19≤ and ≤64",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

p3_age2 <- ggplot(data = growth_rate_age2, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age ≥65",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 8: Save the Plots
ggsave(filename = "growth_rate_mdrEU_age0.tiff", plot = p3_age0, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_mdrEU_age1.tiff", plot = p3_age1, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_mdrEU_age2.tiff", plot = p3_age2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")



max_y<-80

ppred_age0 <- ggplot(data = filter(predictions, Agegroup == 0), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup==0), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup==0), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted MDR (%) - Age ≤18yo",
       x = "Year",
       y = "Predicted MDR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
#≤18yo", "19≤ and ≤64", "≥65"

# Plot for Female
ppred_age1 <- ggplot(data = filter(predictions, Agegroup == 1), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted MDR (%) - Age 19≤ and ≤64",
       x = "Year",
       y = "Predicted MDR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )

ppred_age2 <- ggplot(data = filter(predictions, Agegroup == 2), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup == 2), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup == 2), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted MDR (%) - Age ≥65",
       x = "Year",
       y = "Predicted MDR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
# Display the plot
print(ppred_female)
ggsave(filename = "predictions_breakpoint_mdr_age0.tiff", plot = ppred_age0, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_mdr_age1.tiff", plot = ppred_age1, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_mdr_age2.tiff", plot = ppred_age2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")



######
#------------------------------------------------------------------------------#




#                                                                              #


#. UNITED STATES, below.

#                                                                              #





#Subgroup analyses - Gender-specific change points [United States]

#------------------------------------------------------------------------------#
#NEW, CARBAPENEM RESISTANCE GRAPH: ######
# Step 1: Prepare filtered dataset with Gender trimmed
country_resistance_carbap_us <- state_resistance_carbapGen %>%
  mutate(Gender = str_trim(Gender)) %>%  # Trim whitespace
  filter(Gender == "Male")  # Filter for Male

# Step 2: Set up the model output
model_output_carbap_us <- model_output_carbap_usGen

# Step 3: Prepare Year data
country_resistance_carbap_us$Year <- as.numeric(as.character(country_resistance_carbap_us$Year))
years_range <- seq(min(country_resistance_carbap_us$Year), max(country_resistance_carbap_us$Year), length.out = 100)

# Step 4: Prepare spatial data
states_spatial_filtered <- model_output_carbap_us$states_spatial_filtered
states_data_df <- as.data.frame(states_spatial_filtered@data)
original_geo_levels <- unique(states_data_df$GEOID)

# Step 5: Prepare unique GEOID data
GEOID_nameC <- states_data_df %>%
  dplyr::select(GEOID, NAME) %>%
  distinct(GEOID, .keep_all = TRUE)

# Step 6: Expand prediction data including Gender
# Assuming Gender levels are "Male" and "Female"
gender_levels <- c("Male", "Female")

pred_data <- expand.grid(Year = seq(min(years_range), max(years_range), by = 1),
                         GEOID = original_geo_levels,
                         Gender = gender_levels)  # Include Gender in prediction

# Step 7: Predict from GAM using the created function
predictions <- gam_predictions(model_output_carbap_us$fr, newdata = pred_data)
predictions <- merge(predictions, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions <- predictions %>% arrange(GEOID, Gender, Year)



# Step 8: Repeat for carbapenem resistance predictions
predictions_carb_US <- gam_predictions(model_output_carbap_us$fr, newdata = pred_data)
predictions_carb_US <- merge(predictions_carb_US, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions_carb_US <- predictions_carb_US %>% arrange(GEOID, Gender, Year)

# Step 9: Calculate growth derivatives
resultsGrowth <- derivatives_mh2(model_output_carbap_us$fr, predictions_carb_US)


library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Step 1: Calculate derivatives by Gender
derivatives_data <- derivatives_mh(model_output_carbap_us$fr, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001, startpoint = 0)

# Step 2: Calculate growth rates
derivatives_data <- derivatives_data %>%
  mutate(growth_rate = first_derivative + (second_derivative / first_derivative))

# Step 3: Merge with GEOID data
merged_data <- merge(derivatives_data, unique_GEOID_data, by = "GEOID", all.x = TRUE)
merged_data$first_derivative_sign_change[merged_data$Year == 2004] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2004] <- 0
merged_data$first_derivative_sign_change[merged_data$Year == 2005] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2005] <- 0
# Step 4: Process predictions by Gender
predictions <- predictions %>%
  group_by(GEOID, Gender) %>%  # Grouping by Gender
  mutate(
    pred_lag = lag(pred, default = NA),
    pred_lagu = lag(pred_upper, default = NA),
    pred_lagl = lag(pred_lower, default = NA),
    growth_rate2 = if_else(is.na(pred_lag), NA_real_, 100 * (pred - pred_lag) / pred_lag),
    growth_rate2_up = if_else(is.na(pred_lagu), NA_real_, 100 * (pred_upper - pred_lagu) / pred_lagu),
    growth_rate2_lo = if_else(is.na(pred_lagl), NA_real_, 100 * (pred_lower - pred_lagl) / pred_lagl)
  ) %>%
  dplyr::select(-pred_lag, -pred_lagu, -pred_lagl)  # Clean up lag columns

# Step 5: Store results
Carb_predictions_grat <- predictions
Carb_changep_EU <- merged_data

# Step 6: Plotting by Gender
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Calculate y-axis limits
y_max <- max(merged_data$first_upper, na.rm = TRUE)
y_min <- min(merged_data$first_lower, na.rm = TRUE)
padding <- (y_max - y_min) * 0.05
y_max <- y_max + padding
y_min <- y_min - padding

# Plotting with facets for Gender
library(ggplot2)
library(dplyr)

# Step 1: Filter the merged data by Gender
merged_data_male <- merged_data %>% filter(Gender == "Male")
merged_data_female <- merged_data %>% filter(Gender == "Female")

# Step 2: Plot for Males
p_male <- ggplot(data = merged_data_male, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_male, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Males",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 3: Plot for Females
p_female <- ggplot(data = merged_data_female, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Females",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Save the plots
#ggsave(filename = "first_derivat_eu_carb_male.tiff", plot = p_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "first_derivat_eu_carb_female.tiff", plot = p_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Step 1: Filter merged data by Gender
merged_data_male <- merged_data %>% filter(Gender == "Male")
merged_data_female <- merged_data %>% filter(Gender == "Female")

# Step 2: Define color palette
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 3: Plot for Males - Second Derivative
p2_male2 <- ggplot(data = merged_data_male, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_male, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Males",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Plot for Females - Second Derivative
p2_female2 <- ggplot(data = merged_data_female, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Females",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 5: Save the plots
#ggsave(filename = "second_derivat_eu_carb_male.tiff", plot = p2_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "second_derivat_eu_carb_female.tiff", plot = p2_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")



###REVIEW GROWTH RATES, I presume logs are not calculated and bringing NA values due to 0-values.


library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Step 1: Calculate Growth Rate with Gender
# Step 2: Merge with Predictions and GEOID Data (including Gender)
growth_rate_ci <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE) %>%
  arrange(GEOID, Gender, Year)

# Step 3: Calculate Doubling and Halving Times
growth_rate_ci <- growth_rate_ci %>%
  mutate(
    doubling_times = log(2) / growth_rate2,
    halving_times = log(0.5) / growth_rate2,
    NAME = NAME.x,  # Ensure proper naming
    growth_rate2 = if_else(Year == 2004, NA_real_, growth_rate2),
    growth_rate2_lo = if_else(Year == 2004, NA_real_, growth_rate2_lo),
    growth_rate2_up = if_else(Year == 2004, NA_real_, growth_rate2_up),
    doubling_times = if_else(Year == 2004, NA_real_, doubling_times),
    halving_times = if_else(Year == 2004, NA_real_, halving_times)
  )

# Step 4: Filter Data by Gender
merged_data2 <- growth_rate_ci %>%
  dplyr::left_join(merged_data %>% dplyr::select(GEOID, Year, Gender, first_derivative_sign_change, derivative_breakpoint),
                   by = c("GEOID", "Year", "Gender"))

growth_rate_male <- merged_data2 %>% filter(Gender == "Male")
growth_rate_female <- merged_data2 %>% filter(Gender == "Female")

# Step 5: Define Color Palette
num_colors <- length(unique(growth_rate_ci$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)
growth_rate_male$NAME[growth_rate_male$NAME == "District of Columbia"] <- "D. Columbia"
growth_rate_female$NAME[growth_rate_female$NAME == "District of Columbia"] <- "D. Columbia"

# Step 6: Plot for Males
p3_male <- ggplot(data = growth_rate_male, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_male , derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Males",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 7: Plot for Females
p3_female <- ggplot(data = growth_rate_female, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Females",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 8: Save the Plots
ggsave(filename = "growth_rate_carbUS_male.tiff", plot = p3_male, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_carbUS_female.tiff", plot = p3_female, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")





max_y<-40
predictions$NAME[predictions$NAME == "District of Columbia"] <- "D. Columbia"
merged_data$NAME[merged_data$NAME == "District of Columbia"] <- "D. Columbia"



ppred_male <- ggplot(data = filter(predictions, Gender == "Male"), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Gender == "Male"), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Gender == "Male"), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 10)) +
  labs(title = "Predicted Carbapenem-Resistance (%) - Male",
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )

# Plot for Female
ppred_female <- ggplot(data = filter(predictions, Gender == "Female"), 
                       aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Gender == "Female"), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Gender == "Female"), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 10)) +
  labs(title = "Predicted Carbapenem-Resistance (%) - Female",
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
# Display the plot
print(ppred_female)
ggsave(filename = "predictions_breakpoint_carb_female_US.tiff", plot = ppred_female, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_carb_male_US.tiff", plot = ppred_male, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")



######
#------------------------------------------------------------------------------#
#NEW, CEPHALOSPORIN RESISTANCE GRAPH: ######
# Step 1: Prepare filtered dataset with Gender trimmed
country_resistance_cephalos_us <- state_resistance_cephalosGen %>%
  mutate(Gender = str_trim(Gender)) %>%  # Trim whitespace
  filter(Gender == "Male")  # Filter for Male

# Step 2: Set up the model output
model_output_cephalos_us <- model_output_cephalos_usGen

# Step 3: Prepare Year data
country_resistance_cephalos_us$Year <- as.numeric(as.character(country_resistance_cephalos_us$Year))
years_range <- seq(min(country_resistance_cephalos_us$Year), max(country_resistance_cephalos_us$Year), length.out = 100)

# Step 4: Prepare spatial data
states_spatial_filtered <- model_output_cephalos_us$states_spatial_filtered
states_data_df <- as.data.frame(states_spatial_filtered@data)
original_geo_levels <- unique(states_data_df$GEOID)

# Step 5: Prepare unique GEOID data
GEOID_nameC <- states_data_df %>%
  dplyr::select(GEOID, NAME) %>%
  distinct(GEOID, .keep_all = TRUE)

# Step 6: Expand prediction data including Gender
# Assuming Gender levels are "Male" and "Female"
gender_levels <- c("Male", "Female")

pred_data <- expand.grid(Year = seq(min(years_range), max(years_range), by = 1),
                         GEOID = original_geo_levels,
                         Gender = gender_levels)  # Include Gender in prediction

# Step 7: Predict from GAM using the created function
predictions <- gam_predictions(model_output_cephalos_us$fr, newdata = pred_data)
predictions <- merge(predictions, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions <- predictions %>% arrange(GEOID, Gender, Year)



# Step 8: Repeat for carbapenem resistance predictions
predictions_cephalos_US <- gam_predictions(model_output_cephalos_us$fr, newdata = pred_data)
predictions_cephalos_US <- merge(predictions_cephalos_US, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions_cephalos_US <- predictions_cephalos_US %>% arrange(GEOID, Gender, Year)

# Step 9: Calculate growth derivatives
resultsGrowth <- derivatives_mh2(model_output_cephalos_us$fr, predictions_cephalos_US)


library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Step 1: Calculate derivatives by Gender
derivatives_data <- derivatives_mh(model_output_cephalos_us$fr, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001, startpoint = 0)

# Step 2: Calculate growth rates
derivatives_data <- derivatives_data %>%
  mutate(growth_rate = first_derivative + (second_derivative / first_derivative))

# Step 3: Merge with GEOID data
merged_data <- merge(derivatives_data, unique_GEOID_data, by = "GEOID", all.x = TRUE)
merged_data$first_derivative_sign_change[merged_data$Year == 2004] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2004] <- 0
merged_data$first_derivative_sign_change[merged_data$Year == 2005] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2005] <- 0

# Step 4: Process predictions by Gender
predictions <- predictions %>%
  group_by(GEOID, Gender) %>%  # Grouping by Gender
  mutate(
    pred_lag = lag(pred, default = NA),
    pred_lagu = lag(pred_upper, default = NA),
    pred_lagl = lag(pred_lower, default = NA),
    growth_rate2 = if_else(is.na(pred_lag), NA_real_, 100 * (pred - pred_lag) / pred_lag),
    growth_rate2_up = if_else(is.na(pred_lagu), NA_real_, 100 * (pred_upper - pred_lagu) / pred_lagu),
    growth_rate2_lo = if_else(is.na(pred_lagl), NA_real_, 100 * (pred_lower - pred_lagl) / pred_lagl)
  ) %>%
  dplyr::select(-pred_lag, -pred_lagu, -pred_lagl)  # Clean up lag columns

# Step 5: Store results
Carb_predictions_grat <- predictions
Carb_changep_EU <- merged_data

# Step 6: Plotting by Gender
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Calculate y-axis limits
y_max <- max(merged_data$first_upper, na.rm = TRUE)
y_min <- min(merged_data$first_lower, na.rm = TRUE)
padding <- (y_max - y_min) * 0.05
y_max <- y_max + padding
y_min <- y_min - padding

# Plotting with facets for Gender
library(ggplot2)
library(dplyr)

# Step 1: Filter the merged data by Gender
merged_data_male <- merged_data %>% filter(Gender == "Male")
merged_data_female <- merged_data %>% filter(Gender == "Female")

# Step 2: Plot for Males
p_male <- ggplot(data = merged_data_male, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_male, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Males",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 3: Plot for Females
p_female <- ggplot(data = merged_data_female, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Females",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Save the plots
#ggsave(filename = "first_derivat_eu_carb_male.tiff", plot = p_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "first_derivat_eu_carb_female.tiff", plot = p_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Step 1: Filter merged data by Gender
merged_data_male <- merged_data %>% filter(Gender == "Male")
merged_data_female <- merged_data %>% filter(Gender == "Female")

# Step 2: Define color palette
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 3: Plot for Males - Second Derivative
p2_male2 <- ggplot(data = merged_data_male, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_male, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Males",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Plot for Females - Second Derivative
p2_female2 <- ggplot(data = merged_data_female, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Females",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 5: Save the plots
#ggsave(filename = "second_derivat_eu_carb_male.tiff", plot = p2_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "second_derivat_eu_carb_female.tiff", plot = p2_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")



###REVIEW GROWTH RATES, I presume logs are not calculated and bringing NA values due to 0-values.


library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Step 1: Calculate Growth Rate with Gender
# Step 2: Merge with Predictions and GEOID Data (including Gender)
growth_rate_ci <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE) %>%
  arrange(GEOID, Year)

# Step 3: Calculate Doubling and Halving Times
growth_rate_ci <- growth_rate_ci %>%
  mutate(
    doubling_times = log(2) / growth_rate2,
    halving_times = log(0.5) / growth_rate2,
    NAME = NAME.x,  # Ensure proper naming
    growth_rate2 = if_else(Year == 2004, NA_real_, growth_rate2),
    growth_rate2_lo = if_else(Year == 2004, NA_real_, growth_rate2_lo),
    growth_rate2_up = if_else(Year == 2004, NA_real_, growth_rate2_up),
    doubling_times = if_else(Year == 2004, NA_real_, doubling_times),
    halving_times = if_else(Year == 2004, NA_real_, halving_times)
  )

# Step 4: Filter Data by Gender
merged_data2 <- growth_rate_ci %>%
  dplyr::left_join(merged_data %>% dplyr::select(GEOID, Year, Gender, first_derivative_sign_change, derivative_breakpoint),
                   by = c("GEOID", "Year", "Gender"))

growth_rate_male <- merged_data2 %>% filter(Gender == "Male")
growth_rate_female <- merged_data2 %>% filter(Gender == "Female")


# Step 5: Define Color Palette
num_colors <- length(unique(growth_rate_ci$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

growth_rate_male$NAME[growth_rate_male$NAME == "District of Columbia"] <- "D. Columbia"
growth_rate_female$NAME[growth_rate_female$NAME == "District of Columbia"] <- "D. Columbia"


# Step 6: Plot for Males
p3_male <- ggplot(data = growth_rate_male, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_male , derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Males",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 7: Plot for Females
p3_female <- ggplot(data = growth_rate_female, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Females",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 8: Save the Plots
ggsave(filename = "growth_rate_3GCREU_male_US.tiff", plot = p3_male, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_3GCREU_female_US.tiff", plot = p3_female, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")





max_y<-70
predictions$NAME[predictions$NAME == "District of Columbia"] <- "D. Columbia"
merged_data$NAME[merged_data$NAME == "District of Columbia"] <- "D. Columbia"

ppred_male <- ggplot(data = filter(predictions, Gender == "Male"), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Gender == "Male"), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Gender == "Male"), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted 3GCR (%) - Male",
       x = "Year",
       y = "Predicted 3GCR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )

# Plot for Female
ppred_female <- ggplot(data = filter(predictions, Gender == "Female"), 
                       aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Gender == "Female"), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Gender == "Female"), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted 3GCR (%) - Female",
       x = "Year",
       y = "Predicted 3GCR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
# Display the plot
print(ppred_female)
ggsave(filename = "predictions_breakpoint_3GCR_female_US.tiff", plot = ppred_female, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_3GCR_male_US.tiff", plot = ppred_male, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


######
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#NEW, MDR RESISTANCE GRAPH: ######
# Step 1: Prepare filtered dataset with Gender trimmed
country_resistance_mdr_us <- state_resistance_mdrGen %>%
  mutate(Gender = str_trim(Gender)) %>%  # Trim whitespace
  filter(Gender == "Male")  # Filter for Male

# Step 2: Set up the model output
model_output_mdr_us <- model_output_mdr_usGen

# Step 3: Prepare Year data
country_resistance_mdr_us$Year <- as.numeric(as.character(country_resistance_mdr_us$Year))
years_range <- seq(min(country_resistance_mdr_us$Year), max(country_resistance_mdr_us$Year), length.out = 100)

# Step 4: Prepare spatial data
states_spatial_filtered <- model_output_mdr_us$states_spatial_filtered
states_data_df <- as.data.frame(states_spatial_filtered@data)
original_geo_levels <- unique(states_data_df$GEOID)

# Step 5: Prepare unique GEOID data
GEOID_nameC <- states_data_df %>%
  dplyr::select(GEOID, NAME) %>%
  distinct(GEOID, .keep_all = TRUE)

# Step 6: Expand prediction data including Gender
# Assuming Gender levels are "Male" and "Female"
gender_levels <- c("Male", "Female")

pred_data <- expand.grid(Year = seq(min(years_range), max(years_range), by = 1),
                         GEOID = original_geo_levels,
                         Gender = gender_levels)  # Include Gender in prediction

# Step 7: Predict from GAM using the created function
predictions <- gam_predictions(model_output_mdr_us$sf, newdata = pred_data)
predictions <- merge(predictions, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions <- predictions %>% arrange(GEOID, Gender, Year)



# Step 8: Repeat for carbapenem resistance predictions
predictions_mdr_US <- gam_predictions(model_output_mdr_us$sf, newdata = pred_data)
predictions_mdr_US <- merge(predictions_mdr_US, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions_mdr_US <- predictions_mdr_US %>% arrange(GEOID, Gender, Year)

# Step 9: Calculate growth derivatives
resultsGrowth <- derivatives_mh2(model_output_mdr_us$sf, predictions_mdr_US)


library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Step 1: Calculate derivatives by Gender
derivatives_data <- derivatives_mh(model_output_mdr_us$sf, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001, startpoint = 0)

# Step 2: Calculate growth rates
derivatives_data <- derivatives_data %>%
  mutate(growth_rate = first_derivative + (second_derivative / first_derivative))

# Step 3: Merge with GEOID data
merged_data <- merge(derivatives_data, unique_GEOID_data, by = "GEOID", all.x = TRUE)
merged_data$first_derivative_sign_change[merged_data$Year == 2004] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2004] <- 0
merged_data$first_derivative_sign_change[merged_data$Year == 2005] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2005] <- 0

# Step 4: Process predictions by Gender
predictions <- predictions %>%
  group_by(GEOID, Gender) %>%  # Grouping by Gender
  mutate(
    pred_lag = lag(pred, default = NA),
    pred_lagu = lag(pred_upper, default = NA),
    pred_lagl = lag(pred_lower, default = NA),
    growth_rate2 = if_else(is.na(pred_lag), NA_real_, 100 * (pred - pred_lag) / pred_lag),
    growth_rate2_up = if_else(is.na(pred_lagu), NA_real_, 100 * (pred_upper - pred_lagu) / pred_lagu),
    growth_rate2_lo = if_else(is.na(pred_lagl), NA_real_, 100 * (pred_lower - pred_lagl) / pred_lagl)
  ) %>%
  dplyr::select(-pred_lag, -pred_lagu, -pred_lagl)  # Clean up lag columns

# Step 5: Store results
Carb_predictions_grat <- predictions
Carb_changep_EU <- merged_data

# Step 6: Plotting by Gender
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Calculate y-axis limits
y_max <- max(merged_data$first_upper, na.rm = TRUE)
y_min <- min(merged_data$first_lower, na.rm = TRUE)
padding <- (y_max - y_min) * 0.05
y_max <- y_max + padding
y_min <- y_min - padding

# Plotting with facets for Gender
library(ggplot2)
library(dplyr)

# Step 1: Filter the merged data by Gender
merged_data_male <- merged_data %>% filter(Gender == "Male")
merged_data_female <- merged_data %>% filter(Gender == "Female")

# Step 2: Plot for Males
p_male <- ggplot(data = merged_data_male, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_male, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Males",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 3: Plot for Females
p_female <- ggplot(data = merged_data_female, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Females",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Save the plots
#ggsave(filename = "first_derivat_eu_carb_male.tiff", plot = p_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "first_derivat_eu_carb_female.tiff", plot = p_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(RColorBrewer)
merged_data$first_derivative_sign_change[merged_data$Year == 2004] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2004] <- 0


# Step 1: Filter merged data by Gender
merged_data_male <- merged_data %>% filter(Gender == "Male")
merged_data_female <- merged_data %>% filter(Gender == "Female")

# Step 2: Define color palette
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 3: Plot for Males - Second Derivative
p2_male2 <- ggplot(data = merged_data_male, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_male, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Males",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Plot for Females - Second Derivative
p2_female2 <- ggplot(data = merged_data_female, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Females",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 5: Save the plots
#ggsave(filename = "second_derivat_eu_carb_male.tiff", plot = p2_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "second_derivat_eu_carb_female.tiff", plot = p2_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")



###REVIEW GROWTH RATES, I presume logs are not calculated and bringing NA values due to 0-values.


library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Step 1: Calculate Growth Rate with Gender
# Step 2: Merge with Predictions and GEOID Data (including Gender)
growth_rate_ci <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE) %>%
  arrange(GEOID, Year)

# Step 3: Calculate Doubling and Halving Times
growth_rate_ci <- growth_rate_ci %>%
  mutate(
    doubling_times = log(2) / growth_rate2,
    halving_times = log(0.5) / growth_rate2,
    NAME = NAME.x,  # Ensure proper naming
    growth_rate2 = if_else(Year == 2004, NA_real_, growth_rate2),
    growth_rate2_lo = if_else(Year == 2004, NA_real_, growth_rate2_lo),
    growth_rate2_up = if_else(Year == 2004, NA_real_, growth_rate2_up),
    doubling_times = if_else(Year == 2004, NA_real_, doubling_times),
    halving_times = if_else(Year == 2004, NA_real_, halving_times)
  )

# Step 4: Filter Data by Gender
merged_data2 <- growth_rate_ci %>%
  dplyr::left_join(merged_data %>% dplyr::select(GEOID, Year, Gender, first_derivative_sign_change, derivative_breakpoint),
                   by = c("GEOID", "Year", "Gender"))

growth_rate_male <- merged_data2 %>% filter(Gender == "Male")
growth_rate_female <- merged_data2 %>% filter(Gender == "Female")


# Step 5: Define Color Palette
num_colors <- length(unique(growth_rate_ci$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)
growth_rate_male$NAME[growth_rate_male$NAME == "District of Columbia"] <- "D. Columbia"
growth_rate_female$NAME[growth_rate_female$NAME == "District of Columbia"] <- "D. Columbia"

# Step 6: Plot for Males
p3_male <- ggplot(data = growth_rate_male, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_male, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_male , derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Males",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 7: Plot for Females
p3_female <- ggplot(data = growth_rate_female, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_female, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_female, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Females",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 8: Save the Plots
ggsave(filename = "growth_rate_MDR_male_US.tiff", plot = p3_male, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_MDR_female_US.tiff", plot = p3_female, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")





merged_data$first_derivative_sign_change[merged_data$Year == 2004] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2004] <- 0
max_y<-80

predictions$NAME[predictions$NAME == "District of Columbia"] <- "D. Columbia"
merged_data$NAME[merged_data$NAME == "District of Columbia"] <- "D. Columbia"

ppred_male <- ggplot(data = filter(predictions, Gender == "Male"), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Gender == "Male"), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Gender == "Male"), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted MDR (%) - Male",
       x = "Year",
       y = "Predicted MDR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )

# Plot for Female
ppred_female <- ggplot(data = filter(predictions, Gender == "Female"), 
                       aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Gender == "Female"), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Gender == "Female"), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted MDR (%) - Female",
       x = "Year",
       y = "Predicted MDR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
# Display the plot
print(ppred_female)
ggsave(filename = "predictions_breakpoint_MDR_female_US.tiff", plot = ppred_female, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_MDR_male_US.tiff", plot = ppred_male, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


######
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#Subgroup analyses - Age/group-specific change points [United States]
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#NEW, CARBAPENEM RESISTANCE GRAPH: ######
# Step 1: Prepare filtered dataset with Gender trimmed
country_resistance_carbap_us <- state_resistance_carbapAgeg

# Step 2: Set up the model output
model_output_carbap_us <- model_output_carbap_usAgeg

# Step 3: Prepare Year data
country_resistance_carbap_us$Year <- as.numeric(as.character(country_resistance_carbap_us$Year))
years_range <- seq(min(country_resistance_carbap_us$Year), max(country_resistance_carbap_us$Year), length.out = 100)

# Step 4: Prepare spatial data
states_spatial_filtered <- model_output_carbap_us$states_spatial_filtered
states_data_df <- as.data.frame(states_spatial_filtered@data)
original_geo_levels <- unique(states_data_df$GEOID)

# Step 5: Prepare unique GEOID data
GEOID_nameC <- states_data_df %>%
  dplyr::select(GEOID, NAME) %>%
  distinct(GEOID, .keep_all = TRUE)

# Step 6: Expand prediction data including Agegroup
Ageg_levels <- c(0, 1, 2)

pred_data <- expand.grid(Year = seq(min(years_range), max(years_range), by = 1),
                         GEOID = original_geo_levels,
                         Agegroup = Ageg_levels)  # Include Gender in prediction

# Step 7: Predict from GAM using the created function
predictions <- gam_predictions(model_output_carbap_us$fr, newdata = pred_data)
predictions <- merge(predictions, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions <- predictions %>% arrange(GEOID, Agegroup, Year)



# Step 8: Repeat for carbapenem resistance predictions
predictions_carb_US <- gam_predictions(model_output_carbap_us$fr, newdata = pred_data)
predictions_carb_US <- merge(predictions_carb_US, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions_carb_US <- predictions_carb_US %>% arrange(GEOID, Agegroup, Year)

# Step 9: Calculate growth derivatives
resultsGrowth <- derivatives_mh2(model_output_carbap_us$fr, predictions_carb_US)


library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Step 1: Calculate derivatives by Gender
derivatives_data <- derivatives_mh(model_output_carbap_us$fr, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001, startpoint = 0)

# Step 2: Calculate growth rates
derivatives_data <- derivatives_data %>%
  mutate(growth_rate = first_derivative + (second_derivative / first_derivative))

# Step 3: Merge with GEOID data
merged_data <- merge(derivatives_data, unique_GEOID_data, by = "GEOID", all.x = TRUE)
merged_data$first_derivative_sign_change[merged_data$Year == 2004] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2004] <- 0
merged_data$first_derivative_sign_change[merged_data$Year == 2005] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2005] <- 0
# Step 4: Process predictions by Gender
predictions <- predictions %>%
  group_by(GEOID, Agegroup) %>%  # Grouping by Gender
  mutate(
    pred_lag = lag(pred, default = NA),
    pred_lagu = lag(pred_upper, default = NA),
    pred_lagl = lag(pred_lower, default = NA),
    growth_rate2 = if_else(is.na(pred_lag), NA_real_, 100 * (pred - pred_lag) / pred_lag),
    growth_rate2_up = if_else(is.na(pred_lagu), NA_real_, 100 * (pred_upper - pred_lagu) / pred_lagu),
    growth_rate2_lo = if_else(is.na(pred_lagl), NA_real_, 100 * (pred_lower - pred_lagl) / pred_lagl)
  ) %>%
  dplyr::select(-pred_lag, -pred_lagu, -pred_lagl)  # Clean up lag columns

# Step 5: Store results
Carb_predictions_grat <- predictions
Carb_changep_US <- merged_data

# Step 6: Plotting by Agegroup
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Calculate y-axis limits
y_max <- max(merged_data$first_upper, na.rm = TRUE)
y_min <- min(merged_data$first_lower, na.rm = TRUE)
padding <- (y_max - y_min) * 0.05
y_max <- y_max + padding
y_min <- y_min - padding

# Plotting with facets for Gender
library(ggplot2)
library(dplyr)

# Step 1: Filter the merged data by Agegroup labels = c("≤18yo", "19≤ and ≤64", "≥65")
merged_data_age0 <- merged_data %>% filter(Agegroup == 0)
merged_data_age1 <- merged_data %>% filter(Agegroup == 1)
merged_data_age2 <- merged_data %>% filter(Agegroup == 2)

# Step 2: Plot for Age0
p_age0 <- ggplot(data = merged_data_age0, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age0, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age ≤18yo",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 3: Plot for Females
p_age1<- ggplot(data = merged_data_age1, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age 19≤ and ≤64",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

p_age2<- ggplot(data = merged_data_age2, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age ≥65",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")




# Step 4: Save the plots
#ggsave(filename = "first_derivat_eu_carb_male.tiff", plot = p_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "first_derivat_eu_carb_female.tiff", plot = p_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Step 1: Filter merged data by Gender
merged_data_age0 <- merged_data %>% filter(Agegroup == 0)
merged_data_age1 <- merged_data %>% filter(Agegroup == 1)
merged_data_age2 <- merged_data %>% filter(Agegroup == 2)

# Step 2: Define color palette
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 3: Plot for Males - Second Derivative
p2_age0_2 <- ggplot(data = merged_data_age0, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age0, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age ≤18yo",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Plot for Females - Second Derivative
p2_age1_2 <- ggplot(data = merged_data_age1, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age 19≤ and ≤64",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

p2_age2_2 <- ggplot(data = merged_data_age2, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age ≥65",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")


# Step 5: Save the plots
#ggsave(filename = "second_derivat_eu_carb_male.tiff", plot = p2_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "second_derivat_eu_carb_female.tiff", plot = p2_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Step 1: Calculate Growth Rate with agegroups
# Step 2: Merge with Predictions and GEOID Data (including Agegroups)
growth_rate_ci <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE) %>%
  arrange(GEOID, Year)

# Step 3: Calculate Doubling and Halving Times
growth_rate_ci <- growth_rate_ci %>%
  mutate(
    doubling_times = log(2) / growth_rate2,
    halving_times = log(0.5) / growth_rate2,
    NAME = NAME.x,  # Ensure proper naming
    growth_rate2 = if_else(Year == 2004, NA_real_, growth_rate2),
    growth_rate2_lo = if_else(Year == 2004, NA_real_, growth_rate2_lo),
    growth_rate2_up = if_else(Year == 2004, NA_real_, growth_rate2_up),
    doubling_times = if_else(Year == 2004, NA_real_, doubling_times),
    halving_times = if_else(Year == 2004, NA_real_, halving_times)
  )

# Step 4: Filter Data by Gender
merged_data2 <- growth_rate_ci %>%
  dplyr::left_join(merged_data %>% dplyr::select(GEOID, Year, Agegroup, first_derivative_sign_change, derivative_breakpoint),
                   by = c("GEOID", "Year", "Agegroup"))

growth_rate_age0 <- merged_data2 %>% filter(Agegroup == 0)
growth_rate_age1 <- merged_data2 %>% filter(Agegroup == 1)
growth_rate_age2 <- merged_data2 %>% filter(Agegroup == 2)

#≤18yo", "19≤ and ≤64", "≥65"
# Step 5: Define Color Palette
num_colors <- length(unique(growth_rate_ci$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

growth_rate_age0$NAME[growth_rate_age0$NAME == "District of Columbia"] <- "D. Columbia"
growth_rate_age1$NAME[growth_rate_age1$NAME == "District of Columbia"] <- "D. Columbia"
growth_rate_age2$NAME[growth_rate_age2$NAME == "District of Columbia"] <- "D. Columbia"


# Step 6: Plot for Males
p3_age0 <- ggplot(data = growth_rate_age0, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age0 , derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age ≤18yo",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 7: Plot for Females
p3_age1 <- ggplot(data = growth_rate_age1, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age 19≤ and ≤64",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

p3_age2 <- ggplot(data = growth_rate_age2, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age ≥65",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 8: Save the Plots
ggsave(filename = "growth_rate_carbEU_age0_US.tiff", plot = p3_age0, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_carbEU_age1_US.tiff", plot = p3_age1, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_carbEU_age2_US.tiff", plot = p3_age2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")



max_y<-40
predictions$NAME[predictions$NAME == "District of Columbia"] <- "D. Columbia"
merged_data$NAME[merged_data$NAME == "District of Columbia"] <- "D. Columbia"

ppred_age0 <- ggplot(data = filter(predictions, Agegroup == 0), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup==0), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup==0), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 10)) +
  labs(title = "Predicted Carbapenem-Resistance (%) - Age ≤18yo",
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
#≤18yo", "19≤ and ≤64", "≥65"

# Plot for Female
ppred_age1 <- ggplot(data = filter(predictions, Agegroup == 1), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 10)) +
  labs(title = "Predicted Carbapenem-Resistance (%) - Age 19≤ and ≤64",
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )

ppred_age2 <- ggplot(data = filter(predictions, Agegroup == 2), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup == 2), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup == 2), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 10)) +
  labs(title = "Predicted Carbapenem-Resistance (%) - Age ≥65",
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
# Display the plot
print(ppred_female)
ggsave(filename = "predictions_breakpoint_carb_age0_US.tiff", plot = ppred_age0, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_carb_age1_US.tiff", plot = ppred_age1, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_carb_age2_US.tiff", plot = ppred_age2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


######
#------------------------------------------------------------------------------#
#NEW, CEPHALOSPORIN RESISTANCE GRAPH: ######
# Step 1: Prepare filtered dataset with Gender trimmed
country_resistance_cephalos_us <- state_resistance_cephalosAgeg

# Step 2: Set up the model output
model_output_cephalos_us <- model_output_cephalos_usAgeg

# Step 3: Prepare Year data
country_resistance_cephalos_us$Year <- as.numeric(as.character(country_resistance_cephalos_us$Year))
years_range <- seq(min(country_resistance_cephalos_us$Year), max(country_resistance_cephalos_us$Year), length.out = 100)

# Step 4: Prepare spatial data
states_spatial_filtered <- model_output_cephalos_us$states_spatial_filtered
states_data_df <- as.data.frame(states_spatial_filtered@data)
original_geo_levels <- unique(states_data_df$GEOID)

# Step 5: Prepare unique GEOID data
GEOID_nameC <- states_data_df %>%
  dplyr::select(GEOID, NAME) %>%
  distinct(GEOID, .keep_all = TRUE)

# Step 6: Expand prediction data including Agegroup
Ageg_levels <- c(0, 1, 2)

pred_data <- expand.grid(Year = seq(min(years_range), max(years_range), by = 1),
                         GEOID = original_geo_levels,
                         Agegroup = Ageg_levels)  # Include Gender in prediction

# Step 7: Predict from GAM using the created function
predictions <- gam_predictions(model_output_cephalos_us$fr, newdata = pred_data)
predictions <- merge(predictions, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions <- predictions %>% arrange(GEOID, Agegroup, Year)



# Step 8: Repeat for carbapenem resistance predictions
predictions_cephalos_US <- gam_predictions(model_output_cephalos_us$fr, newdata = pred_data)
predictions_cephalos_US <- merge(predictions_cephalos_US, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions_cephalos_US <- predictions_cephalos_US %>% arrange(GEOID, Agegroup, Year)

# Step 9: Calculate growth derivatives
resultsGrowth <- derivatives_mh2(model_output_cephalos_us$fr, predictions_cephalos_US)


library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Step 1: Calculate derivatives by Gender
derivatives_data <- derivatives_mh(model_output_cephalos_us$fr, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001, startpoint = 0)

# Step 2: Calculate growth rates
derivatives_data <- derivatives_data %>%
  mutate(growth_rate = first_derivative + (second_derivative / first_derivative))

# Step 3: Merge with GEOID data
merged_data <- merge(derivatives_data, unique_GEOID_data, by = "GEOID", all.x = TRUE)
merged_data$first_derivative_sign_change[merged_data$Year == 2004] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2004] <- 0
merged_data$first_derivative_sign_change[merged_data$Year == 2005] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2005] <- 0
# Step 4: Process predictions by Gender
predictions <- predictions %>%
  group_by(GEOID, Agegroup) %>%  # Grouping by Gender
  mutate(
    pred_lag = lag(pred, default = NA),
    pred_lagu = lag(pred_upper, default = NA),
    pred_lagl = lag(pred_lower, default = NA),
    growth_rate2 = if_else(is.na(pred_lag), NA_real_, 100 * (pred - pred_lag) / pred_lag),
    growth_rate2_up = if_else(is.na(pred_lagu), NA_real_, 100 * (pred_upper - pred_lagu) / pred_lagu),
    growth_rate2_lo = if_else(is.na(pred_lagl), NA_real_, 100 * (pred_lower - pred_lagl) / pred_lagl)
  ) %>%
  dplyr::select(-pred_lag, -pred_lagu, -pred_lagl)  # Clean up lag columns

# Step 5: Store results
Carb_predictions_grat <- predictions
Carb_changep_US <- merged_data

# Step 6: Plotting by Agegroup
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Calculate y-axis limits
y_max <- max(merged_data$first_upper, na.rm = TRUE)
y_min <- min(merged_data$first_lower, na.rm = TRUE)
padding <- (y_max - y_min) * 0.05
y_max <- y_max + padding
y_min <- y_min - padding

# Plotting with facets for Gender
library(ggplot2)
library(dplyr)

# Step 1: Filter the merged data by Agegroup labels = c("≤18yo", "19≤ and ≤64", "≥65")
merged_data_age0 <- merged_data %>% filter(Agegroup == 0)
merged_data_age1 <- merged_data %>% filter(Agegroup == 1)
merged_data_age2 <- merged_data %>% filter(Agegroup == 2)

# Step 2: Plot for age0
p_age0 <- ggplot(data = merged_data_age0, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age0, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age ≤18yo",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 3: Plot for age1
p_age1<- ggplot(data = merged_data_age1, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age 19≤ and ≤64",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

p_age2<- ggplot(data = merged_data_age2, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age ≥65",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")


# Step 4: Save the plots
#ggsave(filename = "first_derivat_eu_carb_male.tiff", plot = p_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "first_derivat_eu_carb_female.tiff", plot = p_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Step 1: Filter merged data by Gender
merged_data_age0 <- merged_data %>% filter(Agegroup == 0)
merged_data_age1 <- merged_data %>% filter(Agegroup == 1)
merged_data_age2 <- merged_data %>% filter(Agegroup == 2)

# Step 2: Define color palette
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 3: Plot for Males - Second Derivative
p2_age0_2 <- ggplot(data = merged_data_age0, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age0, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age ≤18yo",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Plot for Females - Second Derivative
p2_age1_2 <- ggplot(data = merged_data_age1, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age 19≤ and ≤64",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

p2_age2_2 <- ggplot(data = merged_data_age2, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age ≥65",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")


# Step 5: Save the plots
#ggsave(filename = "second_derivat_eu_carb_male.tiff", plot = p2_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "second_derivat_eu_carb_female.tiff", plot = p2_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Step 1: Calculate Growth Rate with agegroups
# Step 2: Merge with Predictions and GEOID Data (including Agegroups)
growth_rate_ci <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE) %>%
  arrange(GEOID, Year)

# Step 3: Calculate Doubling and Halving Times
growth_rate_ci <- growth_rate_ci %>%
  mutate(
    doubling_times = log(2) / growth_rate2,
    halving_times = log(0.5) / growth_rate2,
    NAME = NAME.x,  # Ensure proper naming
    growth_rate2 = if_else(Year == 2004, NA_real_, growth_rate2),
    growth_rate2_lo = if_else(Year == 2004, NA_real_, growth_rate2_lo),
    growth_rate2_up = if_else(Year == 2004, NA_real_, growth_rate2_up),
    doubling_times = if_else(Year == 2004, NA_real_, doubling_times),
    halving_times = if_else(Year == 2004, NA_real_, halving_times)
  )

# Step 4: Filter Data by Gender
merged_data2 <- growth_rate_ci %>%
  dplyr::left_join(merged_data %>% dplyr::select(GEOID, Year, Agegroup, first_derivative_sign_change, derivative_breakpoint),
                   by = c("GEOID", "Year", "Agegroup"))

growth_rate_age0 <- merged_data2 %>% filter(Agegroup == 0)
growth_rate_age1 <- merged_data2 %>% filter(Agegroup == 1)
growth_rate_age2 <- merged_data2 %>% filter(Agegroup == 2)

growth_rate_age0$NAME[growth_rate_age0$NAME == "District of Columbia"] <- "D. Columbia"
growth_rate_age1$NAME[growth_rate_age1$NAME == "District of Columbia"] <- "D. Columbia"
growth_rate_age2$NAME[growth_rate_age2$NAME == "District of Columbia"] <- "D. Columbia"


#≤18yo", "19≤ and ≤64", "≥65"
# Step 5: Define Color Palette
num_colors <- length(unique(growth_rate_ci$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 6: Plot for age0
p3_age0 <- ggplot(data = growth_rate_age0, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age0 , derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age ≤18yo",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 7: Plot for age1
p3_age1 <- ggplot(data = growth_rate_age1, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age 19≤ and ≤64",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

p3_age2 <- ggplot(data = growth_rate_age2, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age ≥65",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 8: Save the Plots
ggsave(filename = "growth_rate_cephalosEU_age0_US.tiff", plot = p3_age0, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_cephalosEU_age1_US.tiff", plot = p3_age1, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_cephalosEU_age2_US.tiff", plot = p3_age2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")



max_y<-80
predictions$NAME[predictions$NAME == "District of Columbia"] <- "D. Columbia"
merged_data$NAME[merged_data$NAME == "District of Columbia"] <- "D. Columbia"



ppred_age0 <- ggplot(data = filter(predictions, Agegroup == 0), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup==0), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup==0), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted 3GCR (%) - Age ≤18yo",
       x = "Year",
       y = "Predicted 3GCR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
#≤18yo", "19≤ and ≤64", "≥65"

# Plot for Female
ppred_age1 <- ggplot(data = filter(predictions, Agegroup == 1), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted 3GCR (%) - Age 19≤ and ≤64",
       x = "Year",
       y = "Predicted 3GCR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )

ppred_age2 <- ggplot(data = filter(predictions, Agegroup == 2), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup == 2), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup == 2), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted 3GCRe (%) - Age ≥65",
       x = "Year",
       y = "Predicted 3GCR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
# Display the plot
print(ppred_female)
ggsave(filename = "predictions_breakpoint_cephalos_age0_US.tiff", plot = ppred_age0, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_cephalos_age1_US.tiff", plot = ppred_age1, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_cephalos_age2_US.tiff", plot = ppred_age2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")


######
#------------------------------------------------------------------------------#
#NEW, MDR RESISTANCE GRAPH: ######
# Step 1: 
country_resistance_mdr_us <- state_resistance_mdrAgeg

# Step 2: Set up the model output
model_output_mdr_us <- model_output_mdr_usAgeg

# Step 3: Prepare Year data
country_resistance_mdr_us$Year <- as.numeric(as.character(country_resistance_mdr_us$Year))
years_range <- seq(min(country_resistance_mdr_us$Year), max(country_resistance_mdr_us$Year), length.out = 100)

# Step 4: Prepare spatial data
states_spatial_filtered <- model_output_mdr_us$states_spatial_filtered
states_data_df <- as.data.frame(states_spatial_filtered@data)
original_geo_levels <- unique(states_data_df$GEOID)

# Step 5: Prepare unique GEOID data
GEOID_nameC <- states_data_df %>%
  dplyr::select(GEOID, NAME) %>%
  distinct(GEOID, .keep_all = TRUE)

# Step 6: Expand prediction data including Agegroup
Ageg_levels <- c(0, 1, 2)

pred_data <- expand.grid(Year = seq(min(years_range), max(years_range), by = 1),
                         GEOID = original_geo_levels,
                         Agegroup = Ageg_levels)  # Include Gender in prediction

# Step 7: Predict from GAM using the created function
predictions <- gam_predictions(model_output_mdr_us$sf, newdata = pred_data)
predictions <- merge(predictions, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions <- predictions %>% arrange(GEOID, Agegroup, Year)



# Step 8: Repeat for carbapenem resistance predictions
predictions_mdr_US <- gam_predictions(model_output_mdr_us$sf, newdata = pred_data)
predictions_mdr_US <- merge(predictions_mdr_US, GEOID_nameC, by = "GEOID", all.x = TRUE)
predictions_mdr_US <- predictions_mdr_US %>% arrange(GEOID, Agegroup, Year)

# Step 9: Calculate growth derivatives
resultsGrowth <- derivatives_mh2(model_output_mdr_us$sf, predictions_mdr_EU)


library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Step 1: Calculate derivatives by Gender
derivatives_data <- derivatives_mh(model_output_mdr_us$sf, newdata = pred_data, type = "central", h1 = 0.001, h2 = 0.001, startpoint = 0)

# Step 2: Calculate growth rates
derivatives_data <- derivatives_data %>%
  mutate(growth_rate = first_derivative + (second_derivative / first_derivative))

# Step 3: Merge with GEOID data
merged_data <- merge(derivatives_data, unique_GEOID_data, by = "GEOID", all.x = TRUE)
merged_data$first_derivative_sign_change[merged_data$Year == 2004] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2004] <- 0
merged_data$first_derivative_sign_change[merged_data$Year == 2005] <- 0
merged_data$derivative_breakpoint[merged_data$Year == 2005] <- 0
# Step 4: Process predictions by Gender
predictions <- predictions %>%
  group_by(GEOID, Agegroup) %>%  # Grouping by Gender
  mutate(
    pred_lag = lag(pred, default = NA),
    pred_lagu = lag(pred_upper, default = NA),
    pred_lagl = lag(pred_lower, default = NA),
    growth_rate2 = if_else(is.na(pred_lag), NA_real_, 100 * (pred - pred_lag) / pred_lag),
    growth_rate2_up = if_else(is.na(pred_lagu), NA_real_, 100 * (pred_upper - pred_lagu) / pred_lagu),
    growth_rate2_lo = if_else(is.na(pred_lagl), NA_real_, 100 * (pred_lower - pred_lagl) / pred_lagl)
  ) %>%
  dplyr::select(-pred_lag, -pred_lagu, -pred_lagl)  # Clean up lag columns

# Step 5: Store results
Carb_predictions_grat <- predictions
Carb_changep_EU <- merged_data

# Step 6: Plotting by Agegroup
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Calculate y-axis limits
y_max <- max(merged_data$first_upper, na.rm = TRUE)
y_min <- min(merged_data$first_lower, na.rm = TRUE)
padding <- (y_max - y_min) * 0.05
y_max <- y_max + padding
y_min <- y_min - padding

# Plotting with facets for Gender
library(ggplot2)
library(dplyr)

# Step 1: Filter the merged data by Agegroup labels = c("≤18yo", "19≤ and ≤64", "≥65")
merged_data_age0 <- merged_data %>% filter(Agegroup == 0)
merged_data_age1 <- merged_data %>% filter(Agegroup == 1)
merged_data_age2 <- merged_data %>% filter(Agegroup == 2)

# Step 2: Plot for age0
p_age0 <- ggplot(data = merged_data_age0, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age0, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age ≤18yo",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 3: Plot for age1
p_age1<- ggplot(data = merged_data_age1, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age 19≤ and ≤64",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

p_age2<- ggplot(data = merged_data_age2, aes(x = Year, y = first_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "First Derivative - Age ≥65",
       x = "Year",
       y = "First Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")


# Step 4: Save the plots
#ggsave(filename = "first_derivat_eu_carb_male.tiff", plot = p_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "first_derivat_eu_carb_female.tiff", plot = p_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Step 1: Filter merged data by Gender
merged_data_age0 <- merged_data %>% filter(Agegroup == 0)
merged_data_age1 <- merged_data %>% filter(Agegroup == 1)
merged_data_age2 <- merged_data %>% filter(Agegroup == 2)

# Step 2: Define color palette
num_colors <- length(unique(merged_data$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

# Step 3: Plot for Males - Second Derivative
p2_age0_2 <- ggplot(data = merged_data_age0, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age0, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age ≤18yo",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 4: Plot for Females - Second Derivative
p2_age1_2 <- ggplot(data = merged_data_age1, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age 19≤ and ≤64",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")

p2_age2_2 <- ggplot(data = merged_data_age2, aes(x = Year, y = second_derivative, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = second_lower, ymax = second_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022) +
  labs(title = "Second Derivative - Age ≥65",
       x = "Year",
       y = "Second Derivative") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none"
  ) +
  facet_wrap(~NAME, scales = "free_y")


# Step 5: Save the plots
#ggsave(filename = "second_derivat_eu_carb_male.tiff", plot = p2_male, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")

#ggsave(filename = "second_derivat_eu_carb_female.tiff", plot = p2_female, device = "tiff", path = base_pathOut,
#       width = 11, height = 7, dpi = 500, units = "in")


library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Step 1: Calculate Growth Rate with agegroups
# Step 2: Merge with Predictions and GEOID Data (including Agegroups)
growth_rate_ci <- merge(predictions, unique_GEOID_data, by = "GEOID", all.x = TRUE) %>%
  arrange(GEOID, Year)

# Step 3: Calculate Doubling and Halving Times
growth_rate_ci <- growth_rate_ci %>%
  mutate(
    doubling_times = log(2) / growth_rate2,
    halving_times = log(0.5) / growth_rate2,
    NAME = NAME.x,  # Ensure proper naming
    growth_rate2 = if_else(Year == 2004, NA_real_, growth_rate2),
    growth_rate2_lo = if_else(Year == 2004, NA_real_, growth_rate2_lo),
    growth_rate2_up = if_else(Year == 2004, NA_real_, growth_rate2_up),
    doubling_times = if_else(Year == 2004, NA_real_, doubling_times),
    halving_times = if_else(Year == 2004, NA_real_, halving_times)
  )

# Step 4: Filter Data by Gender
merged_data2 <- growth_rate_ci %>%
  dplyr::left_join(merged_data %>% dplyr::select(GEOID, Year, Agegroup, first_derivative_sign_change, derivative_breakpoint),
                   by = c("GEOID", "Year", "Agegroup"))

growth_rate_age0 <- merged_data2 %>% filter(Agegroup == 0)
growth_rate_age1 <- merged_data2 %>% filter(Agegroup == 1)
growth_rate_age2 <- merged_data2 %>% filter(Agegroup == 2)

#≤18yo", "19≤ and ≤64", "≥65"
# Step 5: Define Color Palette
num_colors <- length(unique(growth_rate_ci$NAME))
palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)

growth_rate_age0$NAME[growth_rate_age0$NAME == "District of Columbia"] <- "D. Columbia"
growth_rate_age1$NAME[growth_rate_age1$NAME == "District of Columbia"] <- "D. Columbia"
growth_rate_age2$NAME[growth_rate_age2$NAME == "District of Columbia"] <- "D. Columbia"


# Step 6: Plot for age0
p3_age0 <- ggplot(data = growth_rate_age0, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age0, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age0 , derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age ≤18yo",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 7: Plot for age1
p3_age1 <- ggplot(data = growth_rate_age1, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age1, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age1, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age 19≤ and ≤64",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

p3_age2 <- ggplot(data = growth_rate_age2, aes(x = Year, y = growth_rate2, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_vline(data = filter(growth_rate_age2, first_derivative_sign_change == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(growth_rate_age2, derivative_breakpoint == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = 2004:2022, labels = as.character(2004:2022)) +
  labs(title = "Growth Rate - Age ≥65",
       x = "Year",
       y = "Growth Rate (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  facet_wrap(~NAME, scales = "free_y")

# Step 8: Save the Plots
ggsave(filename = "growth_rate_mdrEU_age0_US.tiff", plot = p3_age0, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_mdrEU_age1_US.tiff", plot = p3_age1, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")

ggsave(filename = "growth_rate_mdrEU_age2_US.tiff", plot = p3_age2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")



max_y<-80

predictions$NAME[predictions$NAME == "District of Columbia"] <- "D. Columbia"
merged_data$NAME[merged_data$NAME == "District of Columbia"] <- "D. Columbia"

ppred_age0 <- ggplot(data = filter(predictions, Agegroup == 0), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup==0), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup==0), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted MDR (%) - Age ≤18yo",
       x = "Year",
       y = "Predicted MDR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
#≤18yo", "19≤ and ≤64", "≥65"

# Plot for Female
ppred_age1 <- ggplot(data = filter(predictions, Agegroup == 1), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup == 1), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup == 1), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted MDR (%) - Age 19≤ and ≤64",
       x = "Year",
       y = "Predicted MDR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )

ppred_age2 <- ggplot(data = filter(predictions, Agegroup == 2), 
                     aes(x = Year, y = pred, group = NAME)) +
  geom_line(aes(color = NAME), size = 1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = NAME), alpha = 0.2) +
  geom_vline(data = filter(merged_data, first_derivative_sign_change == 1, Agegroup == 2), 
             aes(xintercept = Year), color = "#6baed6", linetype = "dashed", size = 0.5) +
  geom_vline(data = filter(merged_data, derivative_breakpoint == 1, Agegroup == 2), 
             aes(xintercept = Year), color = "#fed98e", linetype = "dashed", size = 0.5) +
  facet_wrap(~NAME, scales = "fixed") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = seq(min(predictions$Year), max(predictions$Year), by = 1)) +
  scale_y_continuous(limits = c(-5, max_y), breaks = seq(0, max_y, by = 20)) +
  labs(title = "Predicted MDR (%) - Age ≥65",
       x = "Year",
       y = "Predicted MDR (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )
# Display the plot
print(ppred_female)
ggsave(filename = "predictions_breakpoint_mdr_age0_US.tiff", plot = ppred_age0, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_mdr_age1_US.tiff", plot = ppred_age1, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")
ggsave(filename = "predictions_breakpoint_mdr_age2_US.tiff", plot = ppred_age2, device = "tiff", path = base_pathOut,
       width = 11, height = 7, dpi = 500, units = "in")



######
#------------------------------------------------------------------------------#




