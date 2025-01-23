#BTC1877H: Data Science in Health II
#Kobayashi Maru 
#By Hasan Abdo

#loading the data set in 
patient_data <- read.csv("~/Documents/Data Science 2/Kobayashi Maru/study2_n1.csv", header = TRUE)

#taking out only the variables of interest 
library(tidyr)
library(dplyr)
patient_data1 <- patient_data %>%
  select(ptid, AGE, COHORT, DXSTAGE,
         T1CHEMO,  
         T1DIABET, T2DIAB, T3DIAB,
         T1EQTOT, T2EQTOT, T3EQTOT,
         T1KIDNEY, T2KIDNEY, T3KIDNEY,
         T1RP, t2rp, t3rp,
         T1RT, t2rt, t3rt,
         T1SPBONE, T2SPBONE, T3SPBONE,
         T1VIAGRA, T2VIAGRA, T3VIAGRA,
         t1ADT, t2ADT, t3ADT,
         t1arthritis, t2arthritis, t3arthritis,
         t1heart, T2HEART, T3HEART,
         t1sporg, T2SPORG, T3SPORG)

#Filtering out only the observations in Cohort A
patient_data1_cohA <- patient_data1 %>%
  filter(COHORT == "A")

#I want to run a mixed effects linear model for my analysis, since using a simple linear model wouldn't work due to the repeated measurements within same patients. 

#Need to change the data frame into a pivoted longer form, in order to get rid of different columns for different time-points, I'll have a time column with T1 = 1, T2= 2, ETC. 

#Matching "Tt" includes both cases of T since the variables are named differently for some
#The number after the "T" or "t" will be stored in a new "time" column, with the variable after it being split into its own column, with observations from the columns matching the "time" row. 
#This way, I have a time variable that I could use in my model down the line. 
patient_data1_arranged <- patient_data1_cohA %>%
  pivot_longer(cols = matches("^[Tt]\\d"),  
    names_to = c("time", ".value"),
    names_pattern = "[Tt]?(\\d)(.*)")

#this does result in some duplicate columns due to the nature of how the variables were initially coded
#using coalesce from the "dplyr" package 
patient_data1_combined <- patient_data1_arranged %>%
  mutate(DIABET = coalesce(DIABET, DIAB),      
    RP = coalesce(RP, rp),               
    RT = coalesce(RT, rt),               
    heart = coalesce(heart, HEART),      
    sporg = coalesce(sporg, SPORG)) %>%
  select(-DIAB, -rp, -rt, -SPORG, -HEART)
#Used select to remove duplicate columns 

#Now doing some data filtering and processing.. 
#Create the "Treatment" based on RP, RT, and ADT in long format

patient_data1_combined <- patient_data1_combined %>%
  mutate(Treatment = case_when(
    RP == 1 & RT == 0 & ADT == 0 ~ "RP Only",
    RP == 0 & RT == 1 & ADT == 0 ~ "RT Only",
    RP == 0 & RT == 0 & ADT == 1 ~ "ADT Only",
    RP == 1 & RT == 1 & ADT == 0 ~ "RP + RT",
    RP == 1 & RT == 0 & ADT == 1 ~ "RP + ADT",
    RP == 0 & RT == 1 & ADT == 1 ~ "RT + ADT",
    RP == 1 & RT == 1 & ADT == 1 ~ "RP + RT + ADT",
    RP == 0 & RT == 0 & ADT == 0 ~ "No Treatment"))

#Filtering:
#Filter out "RP + RT + ADT" from the data (little observations and it doesn't match T2)
#Filter out "ADT Only" from the data (this doesn't make sense according to the instructions)
#Remove patients who had treatment at baseline (time == "1", indicating T1)
#Need to group patients by ptid to ensure the data is even throughout (Meaning I am removign all observations from the patient that meets those conditions I want removed)

patient_data1_combined1 <- patient_data1_combined %>%
  group_by(ptid) %>%
  filter(
    all(Treatment != "RP + RT + ADT"),  
    all(Treatment != "ADT Only"),
    all(!(time == "1" & Treatment != "No Treatment")),
    all(!(time == "3" & Treatment == "RP + RT"))
  ) %>%
  ungroup()

### Doing some data cleaning/re-coding classes### 
library(Hmisc)
describe(patient_data1_combined1)
summary(patient_data1_combined1)

#re-coding variables 
patient_data1_combined1$time <- as.factor(patient_data1_combined1$time)
patient_data1_combined1$CHEMO <- as.factor(patient_data1_combined1$CHEMO)
patient_data1_combined1$DIABET <- as.factor(patient_data1_combined1$DIABET)
patient_data1_combined1$KIDNEY <- as.factor(patient_data1_combined1$KIDNEY)
patient_data1_combined1$RP <- as.factor(patient_data1_combined1$RP)
patient_data1_combined1$RT <- as.factor(patient_data1_combined1$RT)
patient_data1_combined1$SPBONE <- as.factor(patient_data1_combined1$SPBONE)
patient_data1_combined1$VIAGRA <- as.factor(patient_data1_combined1$VIAGRA)
patient_data1_combined1$ADT <- as.factor(patient_data1_combined1$ADT)
patient_data1_combined1$arthritis <- as.factor(patient_data1_combined1$arthritis)
patient_data1_combined1$heart <- as.factor(patient_data1_combined1$heart)
patient_data1_combined1$sporg <- as.factor(patient_data1_combined1$sporg)
patient_data1_combined1$Treatment <- as.factor(patient_data1_combined1$Treatment)


#remove negatives for EQTOT since the range is only from 0-1 (there's 1 observation that falls within this range)

patient_data1_combined$EQTOT[patient_data1_combined$EQTOT < 0] <- 0
patient_data1_combined1$EQTOT[patient_data1_combined1$EQTOT < 0] <- 0


## Creating mixed effects model ###
#loading package for linear mixed effects model
library(lme4)
#to get P value
library(lmerTest)

#with time interaction
mixedeff_model <- lmer(EQTOT ~ time * (Treatment) + AGE + DIABET + KIDNEY + (1 | ptid), data = patient_data1_combined)

summary(mixedeff_model)

#without time interaction
mixedeff_model1 <- lmer(EQTOT ~ Treatment + AGE + time + DIABET + SPBONE + (1 | ptid), data = patient_data1_combined1)

summary(mixedeff_model1)

#Model assumptions
#Model with time interaction: 
#Histogram of the residuals
hist(resid(mixedeff_model))
#Fitted vs residuals 
plot(fitted(mixedeff_model),resid(mixedeff_model))
#QQ-plot
qqnorm(resid(mixedeff_model))
qqline(resid(mixedeff_model), col=2)

#Model without time interaction: 
#Histogram of the residuals
hist(resid(mixedeff_model1))
#Fitted vs residuals 
plot(fitted(mixedeff_model1),resid(mixedeff_model1))
#QQ-plot
qqnorm(resid(mixedeff_model1))
qqline(resid(mixedeff_model1), col=2)

#According to the assumptions, the model(s) is not valid. 


### ---Friedman Test??? --- ### 
#(non parametric since normality isn't achieved)
#Friedman test wont work since all patients need to have observations for all treatment groups, repeated, which doesn't happen here.
friedman_results <- friedman.test(cbind("T1EQTOT", "T2EQTOT", "T3EQTOT") ~ cbind("treatment1","treatment2","treatment3") | ptid, data = patient_krus)
)

#Printing the results
print(friedman_results)

### Trying Kruskal Wallis to find differences at baseline, T2, and T3 ###

patient_krus <- patient_data1_cohA %>%
  select(ptid, 
         T1EQTOT, T2EQTOT, T3EQTOT,
         T1RP, t2rp, t3rp,
         T1RT, t2rt, t3rt,
         t1ADT, t2ADT, t3ADT)

patient_krus <- patient_krus %>%
  mutate(
    treatment1 = case_when(
      T1RP == 1 & T1RT == 0 & t1ADT == 0 ~ "RP Only",
      T1RP == 0 & T1RT == 1 & t1ADT == 0 ~ "RT Only",
      T1RP == 0 & T1RT == 0 & t1ADT == 1 ~ "ADT Only",
      T1RP == 1 & T1RT == 1 & t1ADT == 0 ~ "RP + RT",
      T1RP == 1 & T1RT == 0 & t1ADT == 1 ~ "RP + ADT",
      T1RP == 0 & T1RT == 1 & t1ADT == 1 ~ "RT + ADT",
      T1RP == 1 & T1RT == 1 & t1ADT == 1 ~ "RP + RT + ADT",
      T1RP == 0 & T1RT == 0 & t1ADT == 0 ~ "No Treatment"),
    
    treatment2 = case_when(
      t2rp == 1 & t2rt == 0 & t2ADT == 0 ~ "RP Only",
      t2rp == 0 & t2rt == 1 & t2ADT == 0 ~ "RT Only",
      t2rp == 0 & t2rt == 0 & t2ADT == 1 ~ "ADT Only",
      t2rp == 1 & t2rt == 1 & t2ADT == 0 ~ "RP + RT",
      t2rp == 1 & t2rt == 0 & t2ADT == 1 ~ "RP + ADT",
      t2rp == 0 & t2rt == 1 & t2ADT == 1 ~ "RT + ADT",
      t2rp == 1 & t2rt == 1 & t2ADT == 1 ~ "RP + RT + ADT",
      t2rp == 0 & t2rt == 0 & t2ADT == 0 ~ "No Treatment"),
    
    treatment3 = case_when(
      t3rp == 1 & t3rt == 0 & t3ADT == 0 ~ "RP Only",
      t3rp == 0 & t3rt == 1 & t3ADT == 0 ~ "RT Only",
      t3rp == 0 & t3rt == 0 & t3ADT == 1 ~ "ADT Only",
      t3rp == 1 & t3rt == 1 & t3ADT == 0 ~ "RP + RT",
      t3rp == 1 & t3rt == 0 & t3ADT == 1 ~ "RP + ADT",
      t3rp == 0 & t3rt == 1 & t3ADT == 1 ~ "RT + ADT",
      t3rp == 1 & t3rt == 1 & t3ADT == 1 ~ "RP + RT + ADT",
      t3rp == 0 & t3rt == 0 & t3ADT == 0 ~ "No Treatment"))

#Filter out observations where treatment1, treatment2, or treatment3 is "RP + RT + ADT" along with "ADT Only"
patient_krus_filtered <- patient_krus %>%
  filter(treatment1 != "RP + RT + ADT",
         treatment2 != "RP + RT + ADT",
         treatment3 != "RP + RT + ADT")

patient_krus_filtered <- patient_krus_filtered %>%
  filter(treatment1 != "ADT Only",
         treatment2 != "ADT Only",
         treatment3 != "ADT Only")

#Patients that received RP+RT at T3, these patients switched treatments between T2 and T3
patient_krus_filtered <- patient_krus_filtered %>%
  filter(treatment3 != "RP + RT")

#Also removed patients that had treatment at baseline, since it doesnt make sense for the "newly diagnosed" to be on medication
patient_krus_filtered <- patient_krus_filtered %>%
  filter(treatment1 == "No Treatment")
         

#Re-coding variables to correct classes to run Kruskal test:

patient_krus_filtered$T1RP <- as.factor(patient_krus_filtered$T1RP)
patient_krus_filtered$t2rp <- as.factor(patient_krus_filtered$t2rp)
patient_krus_filtered$t3rp <- as.factor(patient_krus_filtered$t3rp)
patient_krus_filtered$T1RT <- as.factor(patient_krus_filtered$T1RT)
patient_krus_filtered$t2rt <- as.factor(patient_krus_filtered$t2rt)
patient_krus_filtered$t3rt <- as.factor(patient_krus_filtered$t3rt)
patient_krus_filtered$t1ADT <- as.factor(patient_krus_filtered$t1ADT)
patient_krus_filtered$t2ADT <- as.factor(patient_krus_filtered$t2ADT)
patient_krus_filtered$t3ADT <- as.factor(patient_krus_filtered$t3ADT)

patient_krus_filtered$T2EQTOT[patient_krus_filtered$T2EQTOT < 0] <- 0



# Kruskal-Wallis test for T1 (comparing treatments at T1)
kruskal_T1 <- kruskal.test(T1EQTOT ~ treatment1, data = patient_krus_filtered)
print(kruskal_T1)
#-------------------------

# Kruskal-Wallis test for T2 (comparing treatments at T2)
kruskal_T2 <- kruskal.test(T2EQTOT ~ treatment2, data = patient_krus_filtered)
print(kruskal_T2)

# Kruskal-Wallis test for T3 (comparing treatments at T3)
kruskal_T3 <- kruskal.test(T3EQTOT ~ treatment3, data = patient_krus_filtered)
print(kruskal_T3)


#There is no significant difference between any of the treatment groups at any stage 
## Kruskal Plots ###
#Using ggstatsplot package to make the plots with p values
library(ggstatsplot)

#Baseline 
ggbetweenstats(
  data = patient_krus_filtered,
  x = treatment1,     
  y = T1EQTOT,        
  title = "Kruskal-Wallis Test for EQTOT at T1 (Baseline)",
  type = "nonparametric", 
  pairwise.comparisons = FALSE)

#T2
ggbetweenstats(
  data = patient_krus_filtered,
  x = treatment2,
  y = T2EQTOT,
  title = "Kruskal-Wallis Test for EQTOT at T2 (3 Months)",
  type = "nonparametric",
  pairwise.comparisons = FALSE)

#T3
ggbetweenstats(
  data = patient_krus_filtered,
  x = treatment3,
  y = T3EQTOT,
  title = "Kruskal-Wallis Test for EQTOT at T3 (1 Year)",
  type = "nonparametric",
  pairwise.comparisons = FALSE)






##### Plots ####
library(ggplot2)

#BASELINE plots/figures
## EQTOT histogram ##

ggplot(patient_data1_cohA, aes(x = T1EQTOT)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Baseline EQ5D Scores", x = "EQ5D Score", y = "Number of Patients") +
  theme_minimal()

## Treatment distribution ##
baseline_treat <- patient_data1_cohA %>%
  mutate(TreatmentGroup = case_when(
    T1RP == 1 & T1RT == 0 & t1ADT == 0 ~ "RP Only",
    T1RP == 0 & T1RT == 1 & t1ADT == 0 ~ "RT Only",
    T1RP == 0 & T1RT == 0 & t1ADT == 1 ~ "ADT Only",
    T1RP == 1 & T1RT == 1 & t1ADT == 0 ~ "RP + RT",
    T1RP == 1 & T1RT == 0 & t1ADT == 1 ~ "RP + ADT",
    T1RP == 0 & T1RT == 1 & t1ADT == 1 ~ "RT + ADT",
    T1RP == 1 & T1RT == 1 & t1ADT == 1 ~ "RP + RT + ADT",
    T1RP == 0 & T1RT == 0 & t1ADT == 0 ~ "No treatment"))

#Creating bar plot to show frequencies/counts
ggplot(baseline_treat, aes(x = TreatmentGroup)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = "Distribution of Treatment Combinations at Baseline", 
       x = "Treatment Group", 
       y = "Number of Patients") +
  theme_minimal() +
  ylim(0, 100) +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5)

#Creating TreatmentGroup for T2
T2_treat <- patient_data1_cohA %>%
  mutate(TreatmentGroup = case_when(
    t2rp == 1 & t2rt == 0 & t2ADT == 0 ~ "RP Only",
    t2rp == 0 & t2rt == 1 & t2ADT == 0 ~ "RT Only",
    t2rp == 0 & t2rt == 0 & t2ADT == 1 ~ "ADT Only",
    t2rp == 1 & t2rt == 1 & t2ADT == 0 ~ "RP + RT",
    t2rp == 1 & t2rt == 0 & t2ADT == 1 ~ "RP + ADT",
    t2rp == 0 & t2rt == 1 & t2ADT == 1 ~ "RT + ADT",
    t2rp == 1 & t2rt == 1 & t2ADT == 1 ~ "RP + RT + ADT",
    t2rp == 0 & t2rt == 0 & t2ADT == 0 ~ "No treatment"))

#Creating bar plot for T2
ggplot(T2_treat, aes(x = TreatmentGroup)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = "Distribution of Treatment Combinations at T2 (3 Months)", 
       x = "Treatment Group", 
       y = "Number of Patients") +
  theme_minimal() + 
  ylim(0, 100) +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5)

#Creating TreatmentGroup for T3
T3_treat <- patient_data1_cohA %>%
  mutate(TreatmentGroup = case_when(
    t3rp == 1 & t3rt == 0 & t3ADT == 0 ~ "RP Only",
    t3rp == 0 & t3rt == 1 & t3ADT == 0 ~ "RT Only",
    t3rp == 0 & t3rt == 0 & t3ADT == 1 ~ "ADT Only",
    t3rp == 1 & t3rt == 1 & t3ADT == 0 ~ "RP + RT",
    t3rp == 1 & t3rt == 0 & t3ADT == 1 ~ "RP + ADT",
    t3rp == 0 & t3rt == 1 & t3ADT == 1 ~ "RT + ADT",
    t3rp == 1 & t3rt == 1 & t3ADT == 1 ~ "RP + RT + ADT",
    t3rp == 0 & t3rt == 0 & t3ADT == 0 ~ "No treatment"))

# Creating bar plot for T3
ggplot(T3_treat, aes(x = TreatmentGroup)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = "Distribution of Treatment Combinations at T3 (1 Year)", 
       x = "Treatment Group", 
       y = "Number of Patients") +
  theme_minimal() +
  ylim(0, 100) +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5)




#OVERALL PLOTS
#Calculating mean values for EQTOT, but first grouping by time and treatment 
mean_eqtot <- patient_data1_combined %>%
  group_by(time, RP, RT, ADT) %>%
  summarise(mean_EQTOT = mean(EQTOT, na.rm = TRUE))

#Adding a label based on the combination of therapy taken 
mean_eqtot <- mean_eqtot %>%
  mutate(Treatment = case_when(
    RP == 1 & RT == 0 & ADT == 0 ~ "RP Only",
    RP == 0 & RT == 1 & ADT == 0 ~ "RT Only",
    RP == 0 & RT == 0 & ADT == 1 ~ "ADT Only",
    RP == 1 & RT == 1 & ADT == 0 ~ "RP + RT",
    RP == 1 & RT == 0 & ADT == 1 ~ "RP + ADT",
    RP == 0 & RT == 1 & ADT == 1 ~ "RT + ADT",
    RP == 1 & RT == 1 & ADT == 1 ~ "RP + RT + ADT",
    TRUE ~ "None"
  ))


#Plotting mean EQTOT scores over time, grouped by treatment combination
ggplot(mean_eqtot, aes(x = time, y = mean_EQTOT, group = Treatment, color = Treatment)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) + 
  labs(title = "Mean EQTOT Scores Over Time by Treatment Group", 
       x = "Time (T1: Baseline, T2: 3 Months, T3: 1 Year)", 
       y = "Mean EQTOT (Health Utility Score)", 
       color = "Treatment Group") +
  theme_classic() 


#Another plot, boxplot of treatments and how they influence EQTOT over time points
EQTOT_arranged <- patient_data1_combined %>%
  mutate(TreatmentGroup = case_when(
    RP == 1 & RT == 0 & ADT == 0 ~ "RP Only",
    RP == 0 & RT == 1 & ADT == 0 ~ "RT Only",
    RP == 0 & RT == 0 & ADT == 1 ~ "ADT Only",
    RP == 1 & RT == 1 & ADT == 0 ~ "RP + RT",
    RP == 1 & RT == 0 & ADT == 1 ~ "RP + ADT",
    RP == 0 & RT == 1 & ADT == 1 ~ "RT + ADT",
    RP == 1 & RT == 1 & ADT == 1 ~ "RP + RT + ADT",
    TRUE ~ "None"
  ))


count_data <- EQTOT_arranged %>%
  group_by(time, TreatmentGroup) %>%
  summarise(n = n())

#Boxplot for EQTOT, grouped 
ggplot(EQTOT_arranged, aes(x = time, y = EQTOT, fill = TreatmentGroup)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  geom_text(data = count_data, aes(x = time, y = max(EQTOT_arranged$EQTOT) + 0.05, label = paste0("n = ", n)), 
            position = position_dodge(width = 0.75), vjust = -0.5) +  # Add counts with "n = " label
  labs(title = "EQTOT Scores by Treatment Group and Time", 
       x = "Time (T1: Baseline, T2: 3 Months, T3: 1 Year)", 
       y = "EQTOT (Health Utility Score)", 
       fill = "Treatment Group") +
  theme_minimal()




#After processing 
EQTOT_arranged1 <- patient_data1_combined1 %>%
  mutate(TreatmentGroup = case_when(
    RP == 1 & RT == 0 & ADT == 0 ~ "RP Only",
    RP == 0 & RT == 1 & ADT == 0 ~ "RT Only",
    RP == 0 & RT == 0 & ADT == 1 ~ "ADT Only",
    RP == 1 & RT == 1 & ADT == 0 ~ "RP + RT",
    RP == 1 & RT == 0 & ADT == 1 ~ "RP + ADT",
    RP == 0 & RT == 1 & ADT == 1 ~ "RT + ADT",
    RP == 1 & RT == 1 & ADT == 1 ~ "RP + RT + ADT",
    TRUE ~ "No Treatment"
  ))

count_data1 <- EQTOT_arranged1 %>%
  group_by(time, TreatmentGroup) %>%
  summarise(n = n())


ggplot(EQTOT_arranged1, aes(x = time, y = EQTOT, fill = TreatmentGroup)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  geom_text(data = count_data1, aes(x = time, y = max(EQTOT_arranged1$EQTOT) + 0.05, label = paste0("n = ", n)), 
            position = position_dodge(width = 0.75), vjust = -0.5) +  # Add counts with "n = " label
  labs(title = "EQTOT Scores by Treatment Group and Time", 
       x = "Time (T1: Baseline, T2: 3 Months, T3: 1 Year)", 
       y = "EQTOT (Health Utility Score)", 
       fill = "Treatment Group") +
  theme_minimal()


#summary table:
library(gtsummary)

# Recode categorical variables inside the summary_data data frame
summary_data <- patient_data1_combined1 %>%
  select(AGE, Treatment, EQTOT, DIABET, KIDNEY, SPBONE, heart, arthritis, time) %>%
  mutate(
    # Recode Diabetes
    DIABET = factor(DIABET, levels = c(0, 1, 2, 3), labels = c("No", "Yes (Past)", "Yes (Now)", "Not Sure")),
    
    # Recode Kidney Problems
    KIDNEY = factor(KIDNEY, levels = c(0, 1, 2, 3), labels = c("No", "Yes (Past)", "Yes (Now)", "Not Sure")),
    
    # Recode Tumor Spread to Bones
    SPBONE = factor(SPBONE, levels = c(0, 1, 2), labels = c("No", "Yes", "Not Sure")),
    
    # Recode Heart Problems
    heart = factor(heart, levels = c(0, 1, 2, 3), labels = c("No", "Yes (Past)", "Yes (Now)", "Not Sure")),
    
    # Recode Arthritis
    arthritis = factor(arthritis, levels = c(0, 1, 2, 3), labels = c("No", "Yes (Past)", "Yes (Now)", "Not Sure")),
    
    # Treatment already recoded in previous steps
    Treatment = factor(Treatment)
  )

# Create the summary table using gtsummary
summary_table <- summary_data %>%
  tbl_summary(
    by = time,  # Grouping by the time variable (T1: Baseline, T2: 3 Months, T3: 1 Year)
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",  # Median and interquartile range for continuous variables
      all_categorical() ~ "{n} ({p}%)"               # Count and percentage for categorical variables (n (%))
    ),
    label = list(
      AGE ~ "Age",
      Treatment ~ "Treatment",
      EQTOT ~ "EQTOT Score",
      DIABET ~ "Diabetes",
      KIDNEY ~ "Kidney Issues",
      SPBONE ~ "Bone Metastasis",
      heart ~ "Heart Problems",
      arthritis ~ "Arthritis"
    ),
    missing = "no"  # Exclude missing values from the summary
  ) %>%
  add_stat_label() %>%  # Adds a "Statistic" column
  modify_header(
    label = "**Variable**",       # Keep only the variable name in this column
    stat_label = "**Statistic**",  # The column that describes the statistic (Median, Count, etc.)
    stat_1 = "**Baseline**",       # Rename for time == 1
    stat_2 = "**3 Months**",       # Rename for time == 2
    stat_3 = "**1 Year**"          # Rename for time == 3
  ) %>%
  modify_spanning_header(c(stat_1, stat_2, stat_3) ~ "Timepoints") %>%  # Optional: Add a spanning header
  modify_footnote(everything() ~ NA)

# Display the summary table
summary_table


# Recode categorical variables inside the summary_data data frame (non-filtered)
summary_data_non_filtered <- patient_data1_combined %>%
  select(AGE, Treatment, EQTOT, DIABET, KIDNEY, SPBONE, heart, arthritis, time) %>%
  mutate(
    # Recode Diabetes
    DIABET = factor(DIABET, levels = c(0, 1, 2, 3), labels = c("No", "Yes (Past)", "Yes (Now)", "Not Sure")),
    
    # Recode Kidney Problems
    KIDNEY = factor(KIDNEY, levels = c(0, 1, 2, 3), labels = c("No", "Yes (Past)", "Yes (Now)", "Not Sure")),
    
    # Recode Tumor Spread to Bones
    SPBONE = factor(SPBONE, levels = c(0, 1, 2), labels = c("No", "Yes", "Not Sure")),
    
    # Recode Heart Problems
    heart = factor(heart, levels = c(0, 1, 2, 3), labels = c("No", "Yes (Past)", "Yes (Now)", "Not Sure")),
    
    # Recode Arthritis
    arthritis = factor(arthritis, levels = c(0, 1, 2, 3), labels = c("No", "Yes (Past)", "Yes (Now)", "Not Sure")),
    
    # Treatment already recoded in previous steps
    Treatment = factor(Treatment)
  )

# Create the summary table using gtsummary for the non-filtered data
summary_table_non_filtered <- summary_data_non_filtered %>%
  tbl_summary(
    by = time,  # Grouping by the time variable (T1: Baseline, T2: 3 Months, T3: 1 Year)
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",  # Median and interquartile range for continuous variables
      all_categorical() ~ "{n} ({p}%)"               # Count and percentage for categorical variables (n (%))
    ),
    # Correctly matching the variables from the dataset to the labels
    label = list(
      AGE ~ "Age",
      Treatment ~ "Treatment",
      EQTOT ~ "EQTOT Score",
      DIABET ~ "Diabetes",
      KIDNEY ~ "Kidney Issues",
      SPBONE ~ "Bone Metastasis",
      heart ~ "Heart Problems",
      arthritis ~ "Arthritis"
    ),
    missing = "no"  # Exclude missing values from the summary
  ) %>%
  add_stat_label() %>%  # Adds a "Statistic" column for the type of statistic
  modify_header(
    label = "**Variable**",       # Only show variable names in the "Variable" column
    stat_label = "**Statistic**",  # The column that describes the statistic (Median, Count, etc.)
    stat_1 = "**Baseline**",       # Rename for time == 1 to "Baseline"
    stat_2 = "**3 Months**",       # Rename for time == 2 to "3 Months"
    stat_3 = "**1 Year**"          # Rename for time == 3 to "1 Year"
  ) %>%
  modify_spanning_header(c(stat_1, stat_2, stat_3) ~ "Timepoints") %>%  # Add the spanning header for Timepoints
  modify_footnote(everything() ~ NA) 


summary_table_non_filtered %>%
  as_gt() %>%  # Convert to gt table
  gt::tab_options(
    table.font.size = "small"  # Adjust font size to smaller
  )
# Display the summary table for non-filtered data
summary_table_non_filtered




# Create histogram for AGE (continuous variable) over time
ggplot(summary_data_non_filtered, aes(x = AGE)) +
  geom_histogram(binwidth = 5, color = "black", fill = "skyblue", alpha = 0.7) +
  labs(title = "Distribution of Age", x = "Age", y = "Number of Patients") +
  theme_minimal() +
  xlim(0,100)

# Create histogram for EQTOT (continuous variable) over time with modified time labels
ggplot(summary_data_non_filtered, aes(x = EQTOT, fill = time)) +
  geom_histogram(binwidth = 0.05, color = "black", alpha = 0.7) +
  labs(title = "Distribution of EQTOT Scores by Timepoint", x = "EQTOT Score", y = "Number of Patients", fill = "Timepoint") +  # Rename legend to "Timepoint"
  facet_wrap(~ time, ncol = 1) +  # Separate plots for each time point
  theme_minimal()

# Barplot for Treatment distribution over time
ggplot(summary_data_non_filtered, aes(x = Treatment, fill = time)) +
  geom_bar(color = "black", alpha = 0.7, position = "dodge") +
  stat_count(aes(label = paste0("n = ", ..count..)), geom = "text", 
             position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Distribution of Treatments by Timepoint", 
       x = "Treatment Type", y = "Number of Patients", fill = "Timepoint") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 14),   # Increase x-axis labels size
    axis.text.y = element_text(size = 14),   # Increase y-axis labels size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 12),   # Increase legend text size
    plot.title = element_text(size = 18, hjust = 0.5)  # Increase plot title size and center it
  )
###EQTOT PLOTS DISTRIBUTION ####
#Histogram for EQTOT at Time 2 (3 Months) with y-axis limit
ggplot(patient_data1_combined1 %>% filter(time == 2), aes(x = EQTOT)) +
  geom_histogram(bins = 15, fill = "skyblue", color = "black", alpha = 0.7) +  # Set number of bins to 30
  labs(title = "Distribution of EQTOT Scores at 3 Months (Time 2)", 
       x = "EQTOT Score", 
       y = "Number of Patients") +
  ylim(0, 60) +  # Keep y-axis limit at 60
  theme_minimal()

# Histogram for EQTOT at Time 3 (1 Year) without x-axis limit
ggplot(patient_data1_combined1 %>% filter(time == 3), aes(x = EQTOT)) +
  geom_histogram(bins = 15, fill = "lightcoral", color = "black", alpha = 0.7) +  # Set number of bins to 30
  labs(title = "Distribution of EQTOT Scores at 1 Year (Time 3)", 
       x = "EQTOT Score", 
       y = "Number of Patients") +
  ylim(0, 60) +  # Keep y-axis limit at 60
  theme_minimal()



# Barplot for DIABET (Diabetes) with n labels
ggplot(summary_data_non_filtered, aes(x = DIABET, fill = time)) +
  geom_bar(color = "black", alpha = 0.7, position = "dodge") +
  stat_count(aes(label = paste0("n = ", ..count..)), geom = "text", 
             position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Distribution of Diabetes Status by Timepoint", 
       x = "Diabetes Status", y = "Number of Patients", fill = "Timepoint") +
  theme_minimal()

# Barplot for KIDNEY problems with n labels
ggplot(summary_data_non_filtered, aes(x = KIDNEY, fill = time)) +
  geom_bar(color = "black", alpha = 0.7, position = "dodge") +
  stat_count(aes(label = paste0("n = ", ..count..)), geom = "text", 
             position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Distribution of Kidney Problems by Timepoint", 
       x = "Kidney Status", y = "Number of Patients", fill = "Timepoint") +
  theme_minimal()

# Barplot for SPBONE (Tumor spread to bones) with n labels
ggplot(summary_data_non_filtered, aes(x = SPBONE, fill = time)) +
  geom_bar(color = "black", alpha = 0.7, position = "dodge") +
  stat_count(aes(label = paste0("n = ", ..count..)), geom = "text", 
             position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Distribution of Tumor Spread to Bones by Timepoint", 
       x = "Tumor Spread to Bones", y = "Number of Patients", fill = "Timepoint") +
  theme_minimal()

# Barplot for heart problems with n labels
ggplot(summary_data_non_filtered, aes(x = heart, fill = time)) +
  geom_bar(color = "black", alpha = 0.7, position = "dodge") +
  stat_count(aes(label = paste0("n = ", ..count..)), geom = "text", 
             position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Distribution of Heart Problems by Timepoint", 
       x = "Heart Problems", y = "Number of Patients", fill = "Timepoint") +
  theme_minimal()

# Barplot for arthritis with n labels
ggplot(summary_data_non_filtered, aes(x = arthritis, fill = time)) +
  geom_bar(color = "black", alpha = 0.7, position = "dodge") +
  stat_count(aes(label = paste0("n = ", ..count..)), geom = "text", 
             position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Distribution of Arthritis by Timepoint", 
       x = "Arthritis Status", y = "Number of Patients", fill = "Timepoint") +
  theme_minimal()

