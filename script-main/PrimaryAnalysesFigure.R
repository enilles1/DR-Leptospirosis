
# Load packages
```{r}
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggridges)
library(patchwork)
library(grid)
library(finalfit)
library(DescTools)
library(freqtables)
library(RColorBrewer)
library(finalfit)
library(regclass)
library(gridExtra)
```
#Loading dataset
```{r}
# Define  file path
file_path1 <- "[set path]df_lepto_deidentified_FINAL.csv" # Set path to file

# Read the csv file
df_lepto <- read_csv(file_path1)

```
# Demographics
```{r}
seroprev_table <- df_lepto %>%
  group_by(province) %>%
  dplyr::summarize(
    n = n(),
    sum_seropos = sum(seropos, na.rm = TRUE),
    mean_seroprev = mean(seropos, na.rm = TRUE)*100,
    se_seroprev = (sd(seropos, na.rm = TRUE) / sqrt(n))*100,
    lower_95_CI = mean_seroprev - (1.96 * se_seroprev),
    upper_95_CI = mean_seroprev + (1.96 * se_seroprev),
  ) %>%
  arrange(province)

seroprev_table

```
# Demographics
```{r}

# Counts and proportions for gender
gender_proportions <- df_lepto %>%
  group_by(gender) %>%
  dplyr::summarize(count = n()) %>%
  mutate(proportion = count / sum(count))

# Calculate counts and proportions for each province
province_proportions <- df_lepto %>%
  group_by(province) %>%
  dplyr::summarize(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  arrange(desc(proportion)) 

```
# Seroprevalence
```{r}
# Total seroprovalence
# Calculate overall seroprevalence for all serovars combined
overall_positive <- sum(c(8, 1, 4, 3, 38, 25, 1, 1, 20, 2, 38, 40, 5, 3, 4, 3, 41, 0))
overall_negative <- 2091 - overall_positive
overall_seroprevalence <- (overall_positive / 2091) * 100
overall_CI_lower <- overall_seroprevalence - 1.96 * sqrt((overall_seroprevalence/100 * (1 - overall_seroprevalence/100)) / 2091) * 100
overall_CI_upper <- overall_seroprevalence + 1.96 * sqrt((overall_seroprevalence/100 * (1 - overall_seroprevalence/100)) / 2091) * 100

```
#Table 1
```{r}

# List demographic factors
factors <- c("age_category", "gender", "province", "setting", "education", "occupation_grouped", "work_environment", 'animal_contact_rats',   'associated_serogroup')

summarize_seroprevalence <- function(data, factor) {
  # Exclude the specified serogroup
  data_filtered <- data 

  # Calculate seroprevalence
  data_filtered %>%
    group_by(!!sym(factor)) %>%
     dplyr::summarise(
      total = n(),
      positive = sum(seropos == 1, na.rm = TRUE),
      negative = total - positive,
      seroprevalence = (positive / total) * 100,
      CI_lower = seroprevalence - 1.96 * sqrt((seroprevalence / 100 * (1 - seroprevalence / 100)) / total) * 100,
      CI_upper = seroprevalence + 1.96 * sqrt((seroprevalence / 100 * (1 - seroprevalence / 100)) / total) * 100
    ) %>%
    ungroup()
}

# Apply function for factors
summary_tables <- lapply(factors, function(f) summarize_seroprevalence(df_lepto, f))
summary_tables

```
# Table 2 & Table S1
```{r}

df_lepto$num_serovar <- as.factor(df_lepto$num_serovar)

df_lepto$max_titer <- as.numeric(df_lepto$max_titer)

# Function to summarize seroprev
summarize_seroprevalence <- function(data, factor, total_population) {
  
  data %>%
    group_by(!!sym(factor)) %>%
    dplyr::summarise(
      positive = sum(seropos == 1),
      negative = total_population - positive, 
      total = positive + negative, 
      seroprevalence = (positive / total_population) * 100,
      CI_lower = ifelse(positive == 0, 0,
                        seroprevalence - 1.96 * sqrt((seroprevalence/100 * (1 - seroprevalence/100)) / total_population) * 100),
      CI_upper = ifelse(positive == 0, 0,
                        seroprevalence + 1.96 * sqrt((seroprevalence/100 * (1 - seroprevalence/100)) / total_population) * 100),
      mean_max_titer = mean(max_titer, na.rm = TRUE), 
      CI_lower_max_titer = mean_max_titer - 1.96 * sd(max_titer, na.rm = TRUE) / sqrt(length(na.omit(max_titer))), 
      CI_upper_max_titer = mean_max_titer + 1.96 * sd(max_titer, na.rm = TRUE) / sqrt(length(na.omit(max_titer))),
      Median_max_titer = quantile(max_titer, 0.5, na.rm = TRUE), 
      Q1_max_titer = quantile(max_titer, 0.25, na.rm = TRUE), 
      Q3_max_titer = quantile(max_titer, 0.75, na.rm = TRUE)  
    )
}

# Apply fx to  "associated_serogroup" factor
summary_table_serogroup <- summarize_seroprevalence(df_lepto, "associated_serogroup", 2091)
print(summary_table_serogroup)

# Apply fx to  "associated_serovar" factor
summary_table_serovar <- summarize_seroprevalence(df_lepto, "associated_serovar", 2091)
print(summary_table_serovar)


# Associated serogroup
serovar_table <- df_lepto %>%
  group_by(associated_serogroup) %>%
  dplyr::summarise(
    seropos = sum(seropos == 1),
    `Percent of all seropos` = (seropos / total_seropos) * 100
  )%>%
  arrange(desc(seropos)) 
print(serovar_table)

```
#Figure 1
## Fig 1A
```{r}

# Distribution of max titer
max_titer_proportions <- df_lepto %>%
  filter(max_titer != 'NA')%>%
  group_by(max_titer) %>%
  dplyr::summarize(count = n()) %>%
  mutate(proportion = (count / sum(count))*100) %>%
  arrange(desc(proportion))
         
# Define the color palette
lancet_colors <- c(
  "Australis" = "blue",
  "Canicola" = "#d95f02",
  "Djasiman" = "#7570b3",
  "Ictohaem" = "darkred",
  "Pyrogenes" = "#66a61e",
  "Mixed" = "#e6ab02",
  "Other" = "#a6761d"
)

# Filter for seropositive cases
df_lepto_seropos <- df_lepto %>%
  filter(seropos == '1')

# Define the threshold for low count serovars
threshold = 10

# Modify the ordered levels to include a ">=3200" category
ordered_levels_modified <- c(100, 200, 400, 800, 1600, ">=3200")

# Check if the max_titer column exists and is not NULL
if (!"max_titer" %in% names(df_lepto_seropos) || is.null(df_lepto_seropos$max_titer)) {
  stop("Column max_titer does not exist or is NULL in df_lepto_seropos")
}

# Identify serovars with count less than the threshold
low_count_serovars <- df_lepto_seropos %>%
  dplyr::count(associated_serogroup) %>%
  filter(n < threshold) %>%
  pull(associated_serogroup)

# Replace low count serovars with "Other"
df_lepto_seropos$grouped_serogroup <- ifelse(df_lepto_seropos$associated_serogroup %in% low_count_serovars, 
                                             "Other", 
                                             df_lepto_seropos$associated_serogroup)

# Group max_titer values
df_lepto_seropos$grouped_titer <- ifelse(df_lepto_seropos$max_titer >= 3200, ">=3200", as.character(df_lepto_seropos$max_titer))

# Filter for plot data
df_lepto_filtered <- df_lepto_seropos %>%
  filter(!is.na(grouped_titer) & !is.na(grouped_serogroup))

df_lepto_filtered <- df_lepto_filtered %>%
  mutate(grouped_serogroup = ifelse(grouped_serogroup == "Icterohaemorrhagiae", "Ictohaem", grouped_serogroup)) %>%
  mutate(grouped_serogroup = factor(grouped_serogroup, levels = c("Australis", "Canicola", "Djasiman", "Ictohaem", "Pyrogenes", "Mixed", "Other")))


# Generate plot
fig1a <- ggplot(df_lepto_filtered, aes(x = factor(grouped_titer, levels = ordered_levels_modified), 
                                       fill = factor(grouped_serogroup), group = factor(grouped_serogroup))) +
  geom_bar() +
  labs(
    title = 'A',
    x = "Highest titer",
    y = "Number of participants", fill = "Serogroup") +
  theme_classic() +
  scale_fill_manual(values = lancet_colors) +
  theme(legend.position = "right") +
  scale_y_continuous(limits = c(0, 130))


```
## Fig 1B
```{r}

# Create the ridge plot with the custom color palette
fig1b <- ggplot(df_lepto_filtered, aes(x = max_titer, y = grouped_serogroup, fill = grouped_serogroup)) +
  geom_density_ridges() +
  labs(
    title = 'B',
    x = "Highest titer",
    y = NULL,
    fill = "Serovar"
  ) +
  theme_classic() +
  scale_fill_manual(values = lancet_colors) +
  scale_x_continuous(trans = 'log2', 
                     breaks = c(100, 200, 400, 800, 1600, 3200), 
                     labels = c("100", "200", "400", "800", "1600", "3200")) +
  theme(
    legend.position = 'none'
  ) +
  coord_cartesian(xlim = c(120, 3200))

```
## Fig 1C 
```{r}

df_lepto_filtered <- df_lepto_filtered %>%
  mutate(
    age_category = factor(age_category, levels = c("5-19", "20-34", "35-49", "50-64", "65+", "Unknown"))
  )

# Calculate the total count per age category
age_totals <- df_lepto_filtered %>%
  group_by(age_category) %>%
  dplyr::summarise(total = n())

# Calculate the count and proportion for each grouped_serogroup within each age category
titer_summary <- df_lepto_filtered %>%
  group_by(age_category, grouped_serogroup) %>%
  dplyr::summarise(
    count = n(),
    proportion = (count / first(age_totals$total[age_totals$age_category == age_category])) * 100
  )

# Plot
fig1c <- ggplot(titer_summary, aes(fill = grouped_serogroup, y = proportion, x = age_category)) + 
  geom_bar(stat = "identity", position = "fill") +
  labs(title = 'C', y = NULL, x = "Age Category", fill = "Serovar") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = lancet_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

```
##Fig 1D
```{r}

# Count per setting and province
setting_province_totals <- df_lepto_filtered %>%
  group_by(setting, province) %>%
  dplyr::summarise(total = n())

# Count and proportion by associated_serogroup within each setting and province
titer_summary <- df_lepto_filtered %>%
  group_by(setting, province, grouped_serogroup) %>%
  dplyr::summarise(
    count = n(),
    proportion = (count / first(setting_province_totals$total[setting_province_totals$setting == setting & setting_province_totals$province == province])) * 100
  )

# Generate plot
fig1d <- ggplot(titer_summary, aes(x = setting, y = proportion, fill = grouped_serogroup)) + 
  geom_bar(stat = "identity", position = "fill") +
  labs(y = NULL, fill = "Serogroup", title = 'D') +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = lancet_colors) +
  facet_wrap(~ province) +  # Stratify by province
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

```
# Regression models
## Overall (Figure 2, Table S4)
```{r}

# Set reference values
df_lepto$occupation2 <- as.factor(df_lepto$occupation2)
df_lepto$occupation2 <- relevel(df_lepto$occupation2, ref = "Non-professional")
df_lepto$age_category <- as.factor(df_lepto$age_category)
df_lepto$age_category <- relevel(df_lepto$age_category, ref = "5-19")
df_lepto$setting <- as.factor(df_lepto$setting)
df_lepto$setting <- relevel(df_lepto$setting, ref = "Rural")
df_lepto$province <- as.factor(df_lepto$province)
df_lepto$province <- relevel(df_lepto$province, ref = "San Pedro de Macorís")
df_lepto$setting <- as.factor(df_lepto$setting)
df_lepto$setting <- relevel(df_lepto$setting, ref = "Urban")

df_lepto <- df_lepto%>%
  mutate(seropos_cat = as.factor(seropos))

# Explanatory variables
explanatory = c("age_category", "gender", "province", "setting", "occupation2",'animal_contact_rats')
dependent = 'seropos_cat'

df_lepto %>%
  finalfit(dependent, explanatory, metrics=TRUE)

```
## Ictohem (Figue 2, Table S5)
```{r}

df_lepto$Icterohaemorrhagiae_RGA_binary <- ifelse(df_lepto$associated_serogroup == 'Icterohaemorrhagiae', 1, 0)
df_lepto$Icterohaemorrhagiae_RGA_binary[is.na(df_lepto$Icterohaemorrhagiae_RGA_binary)] <- 0 # Replace NA values with 0 
df_lepto$Icterohaemorrhagiae_RGA_binary <- as.factor(df_lepto$Icterohaemorrhagiae_RGA_binary)

dependent = "Icterohaemorrhagiae_RGA_binary"
explanatory = c("age_category", "gender", "province", "setting", "occupation2", 'animal_contact_rats')
df_lepto %>%
  finalfit(dependent, explanatory, metrics=TRUE)

```
## Ictohaem (excluding non-ictohaem seropos)
```{r}

# DataFrame for Icterohaemorrhagiae seropositives
df_seropos_ictohaem <- df_lepto %>%
  filter(associated_serogroup == 'Icterohaemorrhagiae') %>%
  mutate(Icterohaemorrhagiae_RGA_binary = 1)

# DataFrame for seronegatives to all serogroups
df_seroneg <- df_lepto %>%
  filter(seropos == 0) %>%
  mutate(Icterohaemorrhagiae_RGA_binary = 0)

 # Combine the two DataFrames
df_lepto_clean <- bind_rows(df_seropos_ictohaem, df_seroneg)
df_lepto_clean$Icterohaemorrhagiae_RGA_binary <- as.factor(df_lepto_clean$Icterohaemorrhagiae_RGA_binary)

df_lepto_clean$occupation2 <- as.factor(df_lepto_clean$occupation2)
df_lepto_clean$occupation2 <- relevel(df_lepto_clean$occupation2, ref = "Non-professional")
df_lepto_clean$age_category <- as.factor(df_lepto_clean$age_category)
df_lepto_clean$age_category <- relevel(df_lepto_clean$age_category, ref = "5-19")
df_lepto_clean$setting <- as.factor(df_lepto_clean$setting)
df_lepto_clean$setting <- relevel(df_lepto_clean$setting, ref = "Rural")
df_lepto_clean$province <- as.factor(df_lepto_clean$province)
df_lepto_clean$province <- relevel(df_lepto_clean$province, ref = "San Pedro de Macorís")
df_lepto_clean$setting <- as.factor(df_lepto_clean$setting)
df_lepto_clean$setting <- relevel(df_lepto_clean$setting, ref = "Urban")
df_lepto_clean <- df_lepto_clean%>%
  mutate(seropos_cat = as.factor(seropos))

dependent <- "Icterohaemorrhagiae_RGA_binary"
explanatory <- c("age_category", "gender", "province", "setting", "occupation2", "animal_contact_rats")

df_lepto_clean %>%
  finalfit(dependent, explanatory, metrics=TRUE)

```
## Australis (Table S6)
```{r}

df_lepto$australis_M_binary <- ifelse(df_lepto$associated_serogroup == 'Australis', 1, 0)
df_lepto$australis_M_binary[is.na(df_lepto$australis_M_binary)] <- 0 # Replace NA values with 0 in Mankarso_M_binary
df_lepto$australis_M_binary <- as.factor(df_lepto$australis_M_binary)

dependent = "australis_M_binary"
explanatory = c("age_category", "gender", "province", "setting", "occupation2", 'animal_contact_rats')

df_lepto %>%
  finalfit(dependent, explanatory, metrics=TRUE)

```
## Australis (excluding non-australis seropos)
```{r}

# df for Australis
df_seropos_Australis<- df_lepto %>%
  filter(associated_serogroup == 'Australis') %>%
  mutate(Australis_binary = 1)

# df for seronegatives to all serogroups
df_seroneg <- df_lepto %>%
  filter(seropos == 0) %>%
  mutate(Australis_binary = 0)

# Combine df
df_lepto_clean <- bind_rows(df_seropos_Australis, df_seroneg)

df_lepto_clean$Australis_binary <- as.factor(df_lepto_clean$Australis_binary)

dependent <- "Australis_binary"
explanatory <- c("age_category", "gender", "province", "setting", "occupation2", "animal_contact_rats")

df_lepto_clean %>%
  finalfit(dependent, explanatory, metrics=TRUE)

```
## Canocola (Table S7)
```{r}

df_lepto$canicola_binary <- ifelse(df_lepto$associated_serogroup == 'Canicola', 1, 0)
df_lepto$canicola_binary[is.na(df_lepto$canicola_binary)] <- 0
df_lepto$canicola_binary <- as.factor(df_lepto$canicola_binary)

dependent = "canicola_binary"
explanatory = c("age_category", "gender", "province", "setting", "occupation2", 'animal_contact_rats')

df_lepto %>%
  finalfit(dependent, explanatory, metrics=TRUE)
```
## Canicola (excluding non-conacola seropos)
```{r}

df_seropos_Canicola <- df_lepto %>%
  filter(associated_serogroup == 'Canicola') %>%
  mutate(Canicola_binary = 1)

df_seroneg <- df_lepto %>%
  filter(seropos == 0) %>%
  mutate(Canicola_binary = 0)

df_lepto_clean <- bind_rows(df_seropos_Canicola, df_seroneg)
df_lepto_clean$Canicola_binary <- as.factor(df_lepto_clean$Canicola_binary)

dependent <- "Canicola_binary"
explanatory <- c("age_category", "gender", "province", "setting", "occupation2", "animal_contact_rats")

df_lepto_clean %>%
  finalfit(dependent, explanatory, metrics=TRUE)

```
#Figure 2
## Plot overall
```{r}

# Population characteristics
labels <- c("Age: 5-19 y (ref)",
            "        20-34 y", 
            "        35-49 y", 
            "        50-64 y",
            "        65+ y",
            "Gender: Female (ref)",
            "        Male",  
            "Region: SPM (ref)",
            "        Espaillat",
            "Setting: Urban (ref)", 
            "        Rural", 
            "Occup: Non-professional (ref)", 
            "        Farming",
            "        Professional",
            "Contact with rats: No (ref)",
            "         Yes")

# Extracted multivariable OR values
OR <- c(1, 
        3.87,  # from OR = 3.87 (2.05-7/98)
        5.57,  # from OR = 5.57 (2.97-11.41)
        5.07,  # from OR = 5.07 (2.67-10.48)
        6.77,  # from OR = 6.77 (3.58-13.96)
        1, 
        2.41,  # from OR = 2.41 (1.79-3.25)
        1, 
        1.84, # from OR = 1.84 (1.31-2.58)
        1, 
        1.19, # from OR = 1.19 (0.88-1.61)
        1, 
        1.63, # from OR = 1.63 (0.94-2.80)
        1.05, # from OR = 1.05 (0.47-2.08)
        1, 
        1.64) # from OR = 1.64 (1.06-2.52)


# Extracted confidence interval lower values, including '0' for reference categories
lower <- c(1, 2.05, 2.97, 2.67,3.58, 1,1.79, 1,1.31, 1,0.88, 1,0.94,0.47,1,1.06) 
            
upper <- c(1,7.98, 11.41,10.48,13.96,1,3.25,1,2.58, 1,1.61, 1,2.80, 2.08, 1,2.52) 

# Create a column for OR and 95% CI
OR_CI <- sprintf("%.2f (%.2f, %.2f)", OR, lower, upper)

# Create df
df <- data.frame(Variables = labels, OR_CI = OR_CI, OR = OR, Lower = lower, Upper = upper)

df$Variables <- factor(df$Variables, levels = rev(as.character(labels)))

plot <- ggplot(df, aes(x = OR, y = Variables)) +
  geom_vline(xintercept = 1, color = 'gray') +
  geom_point(size = 2.5, shape=22, color = "darkblue", fill="darkblue") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, size = 0.25, color = "black") +
  annotate("text", x = 20, y = df$Variables, label = df$OR_CI, hjust = "left", size = 2.5) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(), 
    axis.text.y = element_text(hjust = 0, size = 8), 
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.title.x = element_text(size = 8), axis.text.x = element_text(size = 7))+
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16), expand = expansion(mult = c(0.01, 0.1)))+
  coord_cartesian(xlim = c(0.25,30)) 

labels_grob <- textGrob(label = "Feature", x = unit(0.01, "npc"), y = unit(-1.8, "npc"), just = "left", gp = gpar(fontface = "bold", fontsize =9))
OR_CI_grob <- textGrob(label = "aOR (95% CI)", x = unit(0.97, "npc"), y = unit(0.1, "npc"), just = "right", gp = gpar(fontface = "bold", fontsize =9))
tit <- textGrob(label = "Overall", x = unit(0.5, "npc"), y = unit(0.1, "npc"), just = "center", gp = gpar(fontface = "bold", fontsize =12))

# Combine main plot and labels
grid.arrange(labels_grob, tit, OR_CI_grob, plot, heights = c(2, 2, 2, 50))

# Combine main plot and labels
combined_plot_overall <- arrangeGrob(labels_grob, tit, OR_CI_grob, fp, heights = c(2, 2, 2, 50))


```
## Plot Ictohem

```{r}

# Extracted OR values
OR <- c(1, 2.91, 3.77, 3.93, 5.83, 1, 2.69, 1, 2.82, 1, 1.62, 1, 0.74, 0.54, 1, 1.73)
lower <- c(1, 1.14,  1.50,  1.55,  2.37,  1, 1.71, 1, 1.65, 1, 1.02, 1, 0.31, 0.09, 1, 0.83)  
upper <- c(1, 8.96,  11.49,  12.02,  17.59,  1, 4.26,  1, 4.99, 1, 2.61, 1, 1.62, 1.85, 1, 3.51)  

# Create data frame 
df <- data.frame(Variables = labels, OR_CI = OR_CI, OR = OR, Lower = lower, Upper = upper)

df$Variables <- factor(df$Variables, levels = rev(as.character(labels)))

plot_icto <- ggplot(df, aes(x = OR, y = Variables)) +
  geom_vline(xintercept = 1, color = 'gray') +
  geom_point(size = 2.5, shape=22, color = "black", fill="red") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, size = 0.25, color = "black") +
  annotate("text", x = 20, y = df$Variables, label = df$OR_CI, hjust = "left", size = 2.5) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(), 
    axis.text.y = element_text(hjust = 0, size = 8), 
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.title.x = element_text(size = 8), axis.text.x = element_text(size = 7))+
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16), expand = expansion(mult = c(0.01, 0.1)))+
  coord_cartesian(xlim = c(0.25,30)) 

# Create labels for above the plot for the second plot
labels_grob_no_labels <- textGrob(label = NULL, x = unit(0.01, "npc"), y = unit(-1.8, "npc"), just = "left", gp = gpar(fontface = "bold", fontsize =10))
OR_CI_grob <- textGrob(label = "", x = unit(0.97, "npc"), y = unit(0.1, "npc"), just = "right", gp = gpar(fontface = "bold", fontsize =9))
tit_no_labels <- textGrob(label = "Icterohaemorrhagiae", x = unit(0.5, "npc"), y = unit(0.1, "npc"), just = "center", gp = gpar(fontface = "bold", fontsize =12))

# Combine main plot and labels
grid.arrange(labels_grob_no_labels, tit_no_labels, OR_CI_grob, plot_icto, heights = c(2, 2, 2, 50))

# Combine main plot and labels
combined_plot_icto <- arrangeGrob(labels_grob_no_labels, tit_no_labels, OR_CI_grob, forest_icto, heights = c(2, 2, 2, 50))


```
#Figure 3
## Seroprev for plots
```{r}

# Fxn to calculate seroprevalence
calculate_seroprevalence <- function(data) {
  age_levels <- c("5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39",
                  "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", 
                  "70-74", "75+")
  data %>%
    dplyr::mutate(
      age_category2 = factor(age_category2, levels = age_levels) 
    ) %>%
    dplyr::group_by(age_category2) %>%
    dplyr::summarise(
      seroprevalence = mean(seropos, na.rm = TRUE), 
      total = n(),                                
      seropos_count = sum(seropos, na.rm = TRUE), 
      .groups = "drop"                             
    ) %>%
    dplyr::mutate(
      std_error = sqrt(seroprevalence * (1 - seroprevalence) / total),
      ci_lower = seroprevalence - qnorm(0.975) * std_error,             
      ci_upper = seroprevalence + qnorm(0.975) * std_error            
    )
}

# overall data
overall_seroprevalence <- calculate_seroprevalence(df_lepto)

# San Pedro de Macorís
san_pedro_seroprevalence <- calculate_seroprevalence(filter(df_lepto, province == 'San Pedro de Macorís'))

# Espaillat
espaillat_seroprevalence <- calculate_seroprevalence(filter(df_lepto, province == 'Espaillat'))

# Female
female_seroprevalence <- calculate_seroprevalence(filter(df_lepto, gender == 'Female'))

# Male
male_seroprevalence <- calculate_seroprevalence(filter(df_lepto, gender == 'Male'))

# Urban
urban_seroprevalence <- calculate_seroprevalence(filter(df_lepto, setting == 'Urban'))

# Rural
rural_seroprevalence <- calculate_seroprevalence(filter(df_lepto, setting == 'Rural'))

```
## Plots
```{r}
# Overall
age_p1 <- ggplot(overall_seroprevalence, aes(x = age_category2, y = seroprevalence, weight = total)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0) +
    geom_smooth(method = "gam", formula = y ~ s(x),  # Using a cubic spline with k knots
              se = TRUE, color = "#1d3557", fill = "#1d3557", aes(group = 1)) +
  labs(x = "Age, years", y = "Proportion seropositive", title = "Overall") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold')) +
  coord_cartesian(ylim = c(0.015, 0.35)) 

# San Pedro de Macorís
age_p2 <- ggplot(san_pedro_seroprevalence, aes(x = age_category2, y = seroprevalence, weight = total)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0) +
    geom_smooth(method = "gam", formula = y ~ s(x),  
              se = TRUE, color = "#2a9d8f", fill = "#2a9d8f", aes(group = 1)) +
  labs(x = NULL, y = NULL, title = "Southeast") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5, face = 'bold')) +
  coord_cartesian(ylim = c(0.015, 0.4)) 

# Espaillat
age_p3 <- ggplot(espaillat_seroprevalence, aes(x = age_category2, y = seroprevalence, weight = total)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0) +
    geom_smooth(method = "gam", formula = y ~ s(x), 
              se = TRUE, color = "#2a9d8f", fill = "#2a9d8f", aes(group = 1)) +
  labs(x = "", y = NULL, title = "Northwest") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5, face = 'bold')) +
  coord_cartesian(ylim = c(0.015, 0.4)) 

# Female
gender_p1 <- ggplot(female_seroprevalence, aes(x = age_category2,y = seroprevalence, weight = total)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0) +
    geom_smooth(method = "gam", formula = y ~ s(x),  
              se = TRUE, color = "#457b9d", fill = "#457b9d", aes(group = 1)) +
  labs(x = "", y = NULL, title = "Female") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5, face = 'bold')) +
  coord_cartesian(ylim = c(0.015, 0.4)) 

# Male
gender_p2 <- ggplot(male_seroprevalence, aes(x = age_category2,y = seroprevalence, weight = total)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0) +
    geom_smooth(method = "gam", formula = y ~ s(x), 
              se = TRUE, color = "#457b9d", fill = "#457b9d", aes(group = 1)) +
  labs(x = "", y = NULL, title = "Male") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5, face = 'bold')) +
  coord_cartesian(ylim = c(0.015, 0.4)) 

# Urban
setting_p1 <- ggplot(urban_seroprevalence, aes(x = age_category2, y = seroprevalence, weight = total)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0) +
    geom_smooth(method = "gam", formula = y ~ s(x), 
              se = TRUE, color = "#F4A261", fill = "#F4A261",  aes(group = 1)) +
  labs(x = "Age, years", y = NULL, title = "Urban") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5, face = 'bold')) +
  coord_cartesian(ylim = c(0.015, 0.4)) 

# Rural
setting_p2 <- ggplot(rural_seroprevalence, aes(x = age_category2, y = seroprevalence, weight = total)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0) +
    geom_smooth(method = "gam", formula = y ~ s(x), 
              se = TRUE, color = "#F4A261", fill = "#F4A261",  aes(group = 1)) +
  labs(x = "Age, years", y = NULL, title = "Rural") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5, face = 'bold')) +
  coord_cartesian(ylim = c(0.015, 0.4)) 

```
