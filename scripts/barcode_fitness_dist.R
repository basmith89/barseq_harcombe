#Mereging barcode dataframes to plot fitness distributions

library(dplyr)
library(ggplot2)
library(gridExtra)

setwd("/home/brian/Documents/Work/barseq_proj/data/")

eo_t0_df <- read.csv("Eo_T0.codes", sep = "\t")
ES_t1_r1 <- read.csv("ES_R1.codes", sep = "\t")
ES_t1_r2 <- read.csv("ES_R2.codes", sep = "\t")
ES_t1_r3 <- read.csv("ES_R3.codes", sep = "\t")
ES_t1_r4 <- read.csv("ES_R4.codes", sep = "\t")
ES_t1_r5 <- read.csv("ES_R5.codes", sep = "\t")

#merged_df <- full_join(eo_t0_df, ES_t1_r1, ES_t1_r2, ES_t1_r3, ES_t1_r4, ES_t1_r5, by = c("barcode"))

#had to nest left_join to add all data sets together
merged_df <- left_join(eo_t0_df, ES_t1_r1, by="barcode") %>%
  left_join(., ES_t1_r2, by="barcode") %>%
  left_join(., ES_t1_r3, by="barcode") %>%
  left_join(., ES_t1_r4, by="barcode") %>%
  left_join(., ES_t1_r5, by="barcode")

#full join will append each df to each other -- not what I want here
#merged_df2 <- full_join(eo_t0_df, ES_t1_r1, by="barcode") %>%
#  full_join(., ES_t1_r2, by="barcode") %>%
#  full_join(., ES_t1_r3, by="barcode") %>%
#  full_join(., ES_t1_r4, by="barcode") %>%
#  full_join(., ES_t1_r5, by="barcode")

#replace all NA's with 0
merged_df <- merged_df %>% mutate_all(~replace(., is.na(.), 0))

#Store sum of each column
Eo_sum <- sum(merged_df$IT001)
ES_R1_sum <-sum(merged_df$IT007)
ES_R2_sum <-sum(merged_df$IT008)
ES_R3_sum <-sum(merged_df$IT009)
ES_R4_sum <-sum(merged_df$IT0010)
ES_R5_sum <-sum(merged_df$IT011)

#Create frequency columns
merged_df <- merged_df %>% mutate(Eo_T0_fq = IT001/Eo_sum)
merged_df <- merged_df %>% mutate(ES_R1_fq = IT007/ES_R1_sum)
merged_df <- merged_df %>% mutate(ES_R2_fq = IT008/ES_R2_sum)
merged_df <- merged_df %>% mutate(ES_R3_fq = IT009/ES_R3_sum)
merged_df <- merged_df %>% mutate(ES_R4_fq = IT0010/ES_R4_sum)
merged_df <- merged_df %>% mutate(ES_R5_fq = IT011/ES_R5_sum)

#create fitness columns
merged_df <- merged_df %>% mutate(R1ΔT0 = abs(ES_R1_fq - Eo_T0_fq))
merged_df <- merged_df %>% mutate(R2ΔT0 = abs(ES_R2_fq - Eo_T0_fq))
merged_df <- merged_df %>% mutate(R3ΔT0 = abs(ES_R3_fq - Eo_T0_fq))
merged_df <- merged_df %>% mutate(R4ΔT0 = abs(ES_R4_fq - Eo_T0_fq))
merged_df <- merged_df %>% mutate(R5ΔT0 = abs(ES_R5_fq - Eo_T0_fq))

merged_df <- merged_df %>% mutate(R1dbT0 = ES_R1_fq / Eo_T0_fq)
merged_df <- merged_df %>% mutate(R2dbT0 = ES_R2_fq / Eo_T0_fq)
merged_df <- merged_df %>% mutate(R3dbT0 = ES_R3_fq / Eo_T0_fq)
merged_df <- merged_df %>% mutate(R4dbT0 = ES_R4_fq / Eo_T0_fq)
merged_df <- merged_df %>% mutate(R5dbT0 = ES_R5_fq / Eo_T0_fq)

#log2
merged_df <- merged_df %>% mutate(R1dbT0 = log2(ES_R1_fq / Eo_T0_fq))
merged_df <- merged_df %>% mutate(R2dbT0 = log2(ES_R2_fq / Eo_T0_fq))
merged_df <- merged_df %>% mutate(R3dbT0 = log2(ES_R3_fq / Eo_T0_fq))
merged_df <- merged_df %>% mutate(R4dbT0 = log2(ES_R4_fq / Eo_T0_fq))
merged_df <- merged_df %>% mutate(R5dbT0 = log2(ES_R5_fq / Eo_T0_fq))

#trying to plot data
ggplot(data = merged_df, mapping = aes(x = R1ΔT0, y = ES_R1_fq)) + geom_point()

ggplot(data = merged_df, mapping = aes(x = R1dbT0, y = ES_R1_fq)) + geom_point()

#can subset data to get read of things with a value of 0 in one column
test_df <- subset(merged_df, R1dbT0 > 0)

#selecting specific columns
write_test <- select(test_df, ES_R1_fq, R1dbT0)

#plot multiple graphs into one fig. Below is example of storing plots to var for grid.arrange
#p11 <-ggplot(data = merged_df, mapping = aes(x = log(IT011))) + geom_histogram(binwidth = .1)
grid.arrange(p1, p7, p8, p9, p10, p11, ncol = 3)
