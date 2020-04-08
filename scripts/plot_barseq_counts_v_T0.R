#Plotting raw BarSeq counts vs T0

library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyr)


##inputs
setwd("/Users/briansmith/Documents/PostDoc_Harcombe/TnSeq/BarSeq/data/multicodes_out/")
min_filter = 4
max_filter = 40000

file_list <- list.files(".", pattern = "*.codes")
#
##df <- "df"
#finds the where T0 sample is in the list
T0_index <- grep(pattern = "*T0.codes", file_list)
#moves list index to top of the list
#T_0 <- file_list[c()]
file_names <- setdiff(file_list, file_list[T0_index])
file_names <- append(file_names, values = file_list[T0_index], after = match(file_list[1], file_names) - 1)
#T_0 <- append(list(file_list[T0_index]), file_list)


c <- 0
#loop for adding multiple files
for (f in file_names) {
  c <- c + 1
  fname <- f
  
  #msg to user
  print(paste0("Adding ", fname, " into a merged data frame"))
  
  #logic gate
  #first if loads intial data ##MAKE SURE THIS IS TIME 0##
  #second if adds samples to T0
  if (c == 1) {
    merged_df <- read.csv(f, sep = "\t")
    #adjust these filters as needed
    merged_df <- subset(merged_df, merged_df[ , 2] > min_filter)
    merged_df <- subset(merged_df, merged_df[ , 2] < max_filter)
  }
  else if (c >= 2)  {
    df2 <- read.csv(f, sep = "\t")
    #adjust these filters as needed
    df2 <- subset(df2, df2[ , 2] > min_filter)
    df2 <- subset(df2, df2[ , 2] < max_filter)
    merged_df <- left_join(merged_df, df2, by="barcode")
    
    #this one line will work if no threshold filtering is necessary
    #df <- left_join(df, read.csv(f, sep = "\t"), by="barcode")
  }
}

#fix Multicodes.pl bug that adds extra zero to index 10
if ("IT0010" %in% colnames(merged_df)) {
  merged_df <-rename(merged_df, IT010 = IT0010)
}

#replace all NA's with 0
merged_df <- merged_df  %>% mutate_all(~replace(., is.na(.), 0))

##store col names and remove the barcode character
#sample_ids <- tail(colnames(merged_df), -1)


#ptest <- ggplot(merged_df, aes(x = merged_df$IT001, y = merged_df$sample_ids[2])) +  geom_point() 

#plot all samples together against IT001
new_df <- merged_df %>%
  gather(index, counts, -barcode, -IT001)

new_df %>%
  ggplot(aes(x = IT001, y = counts)) +
  facet_wrap(~ index) +
  geom_point()



##############
#The long way#
##############

eo_t0_df <- read.csv("Eo_T0.codes", sep = "\t")
ES_t1_r1 <- read.csv("ES_R1.codes", sep = "\t")
ES_t1_r2 <- read.csv("ES_R2.codes", sep = "\t")
ES_t1_r3 <- read.csv("ES_R3.codes", sep = "\t")
ES_t1_r4 <- read.csv("ES_R4.codes", sep = "\t")
ES_t1_r5 <- read.csv("ES_R5.codes", sep = "\t")

#filter out low value counts
eo_t0_df <- subset(eo_t0_df, eo_t0_df[ , 2] > 4)
ES_t1_r1 <- subset(ES_t1_r1, ES_t1_r1[ , 2] > 4)
ES_t1_r2 <- subset(ES_t1_r2, ES_t1_r2[ , 2] > 4)
ES_t1_r3 <- subset(ES_t1_r3, ES_t1_r3[ , 2] > 4)
ES_t1_r4 <- subset(ES_t1_r4, ES_t1_r4[ , 2] > 4)
ES_t1_r5 <- subset(ES_t1_r5, ES_t1_r5[ , 2] > 4)

merged_df <- left_join(eo_t0_df, ES_t1_r1, by="barcode") %>%
  left_join(., ES_t1_r2, by="barcode") %>%
  left_join(., ES_t1_r3, by="barcode") %>%
  left_join(., ES_t1_r4, by="barcode") %>%
  left_join(., ES_t1_r5, by="barcode")

##replace all NA's with 0
merged_df <- merged_df %>% mutate_all(~replace(., is.na(.), 0))

#store plotes to variables
p1 <- ggplot(data = merged_df, mapping = aes(x = IT001, y = IT007)) + geom_point()
p2 <- ggplot(data = merged_df, mapping = aes(x = IT001, y = IT008)) + geom_point()
p3 <- ggplot(data = merged_df, mapping = aes(x = IT001, y = IT009)) + geom_point()
p4 <- ggplot(data = merged_df, mapping = aes(x = IT001, y = IT0010)) + geom_point()
p5 <- ggplot(data = merged_df, mapping = aes(x = IT001, y = IT011)) + geom_point()


#use varibles plot multiplot w/ gridarrange
grid.arrange(p1, p2, p3, p4, p5, ncol = 2)


##################
####Calc means####
##################

merged_df$E_mean <- rowMeans(subset(merged_df, select = c(IT002, IT003, IT004, IT005, IT006)), na.rm = TRUE)
merged_df$ES_mean <- rowMeans(subset(merged_df, select = c(IT007, IT008, IT009, IT010, IT011)), na.rm = TRUE)
merged_df$EM_mean <- rowMeans(subset(merged_df, select = c(IT012, IT013, IT014, IT015, IT016)), na.rm = TRUE)
merged_df$ESM_mean <- rowMeans(subset(merged_df, select = c(IT017, IT018, IT019, IT020, IT021)), na.rm = TRUE)


p1 <- ggplot(data = merged_df, mapping = aes(x = E_mean, y = ES_mean)) + geom_point()
p2 <- ggplot(data = merged_df, mapping = aes(x = E_mean, y = EM_mean)) + geom_point()
p3 <- ggplot(data = merged_df, mapping = aes(x = E_mean, y = ESM_mean)) + geom_point()
p4 <- ggplot(data = merged_df, mapping = aes(x = ES_mean, y = EM_mean)) + geom_point()
p5 <- ggplot(data = merged_df, mapping = aes(x = ES_mean, y = ESM_mean)) + geom_point()
p6 <- ggplot(data = merged_df, mapping = aes(x = EM_mean, y = ESM_mean)) + geom_point()

grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)


##################
####Calc ratios###
##################

#create ratios
merged_df$E_mean_div_IT001 <- round(as.numeric(merged_df$E_mean) / merged_df$IT001, digits =  4)
merged_df$ES_mean_div_IT001 <- round(as.numeric(merged_df$ES_mean) / merged_df$IT001, digits =  4)
merged_df$EM_mean_div_IT001 <- round(as.numeric(merged_df$EM_mean) / merged_df$IT001, digits =  4)
merged_df$ESM_mean_div_IT001 <- round(as.numeric(merged_df$ESM_mean) / merged_df$IT001, digits =  4)

#build ratio df with barcodes
ratio_df <- select(merged_df, barcode, E_mean_div_IT001, ES_mean_div_IT001, EM_mean_div_IT001, ESM_mean_div_IT001)


#filter values y/x > 2 and y/x < 1/2

ratio_df <- subset(ratio_df, ratio_df$E_mean_div_IT001 >= 2 | ratio_df$E_mean_div_IT001 <= 0.5)
ratio_df <- subset(ratio_df, ratio_df$ES_mean_div_IT001 >= 2 | ratio_df$ES_mean_div_IT001 <= 0.5)
ratio_df <- subset(ratio_df, ratio_df$EM_mean_div_IT001 >= 2 | ratio_df$EM_mean_div_IT001 <= 0.5)
ratio_df <- subset(ratio_df, ratio_df$ESM_mean_div_IT001 >= 2 | ratio_df$ESM_mean_div_IT001 <= 0.5)

#filter zeros out
ratio_df <- ratio_df[ratio_df$E_mean_div_IT001 != 0, ]
ratio_df <- ratio_df[ratio_df$ES_mean_div_IT001 != 0, ]
ratio_df <- ratio_df[ratio_df$EM_mean_div_IT001 != 0, ]
ratio_df <- ratio_df[ratio_df$ESM_mean_div_IT001 != 0, ]

##################
####Diagnostics###
##################

dep_means <- names(merged_df)[27:30]

setwd("../g/figures/")
for (d in dep_means) {
  print(d)
  #need get() to convert string list to an object
  linearMod <- lm(IT001 ~ get(d), data = merged_df)
  print(linearMod)
  summary(linearMod)
  #plot(linearMod, pch = 18, col = "red", sub.caption = paste("IT001 ~", d, sep = " "))
  #old <- par(oma = c(0,0,2,0), mfrow = c(3,2))
  png(filename = paste("IT001_vs_", d, ".png", sep = ""), width = 1920, height = 1080)
  par(oma = c(0,0,2,0), mfrow = c(3,2))
  plot(linearMod, pch = 18, col = "red", which = c(1,2,3,4,5,6), ask = FALSE, sub.caption = paste("Time0 ~", d, sep = " "), cex.id = 1.30, cex.oma.main = 1.75, cex.caption = 1.50)
  dev.off()
  #mtext("Test title", outer = TRUE, cex = 1.5)
  #par(old)
}



######The long way

linearMod <- lm (IT001 ~ IT003, data = merged_df)
print(linearMod)
summary(linearMod)
plot(linearMod, pch = 18, col = "red")
plot(linearMod, pch = 18, col = "red", which =c(4))