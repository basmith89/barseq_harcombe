#notes

library(Tnseq)
library(tidyr)

barseq_pool <- read.csv("~/Documents/PostDoc_Harcombe/TnSeq/BarSeq/data/mg1655/combined_mg1655_data.poolcount", sep = "\t")

bar_pool <- barseq_pool[,-c(1,2,3,4)]

bar_pool <- rename(bar_pool, IT010 = IT0010)

#reaarange columns in numerical order
bar_pool1 <- bar_pool[,c(1,26,20,21,22,23,24,14,15,16,17,18,2,3,4,5,6,8,9,10,11,12,25,19,7,13)]

#filter out rows that contain 0 counts for EVERY sample to reduce data size
bar_pool1 <- bar_pool1[rowSums((bar_pool1[,2:26])) > 0,]


#check for NAs in df
apply(bar_pool1, 2, function(x) any(is.na(x)))


#remove NA data
#need tidyr
bar_pool1 <- drop_na(bar_pool1)

#write this file out and use python script to create TnSeq_diff data format
write.table(bar_pool1, "~/Documents/PostDoc_Harcombe/TnSeq/BarSeq/data/mg1655/mg1655_pool_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE, dec = ".")

barseq_diff <- read.csv("~/Documents/PostDoc_Harcombe/TnSeq/BarSeq/data/mg1655/tnseq_diff/mg1655_pool_for_tnseqdiff.tsv", sep = "\t")


#seperate data
Eo_barseq <- bar_pool1[,c(2,3,4,5,6,7)]
ES_barseq <- bar_pool1[,c(2,8,9,10,11,12)]
EM_barseq <- bar_pool1[,c(2,15,14,15,16,17)]
ESM_barseq <- bar_pool1[,c(2,18,19,20,21,22)]

condition = c(rep("Input", 1), rep("Output", 5))
geneID = as.character(barseq_diff$GeneID)
location=barseq_diff$pos
pool = c(1,1,1,1,1,1)

Eo_diff <- TnseqDiff(Eo_barseq, geneID, location, pool, condition)
ES_diff <- TnseqDiff(ES_barseq, geneID, location, pool, condition)
EM_diff <- TnseqDiff(EM_barseq, geneID, location, pool, condition)
ESM_diff <- TnseqDiff(ESM_barseq, geneID, location, pool, condition)

#search for specific geneIDs
dplyr::filter(barseq_diff, grepl('BW25113_RS15410', geneID))

countData = barseq_diff[,c(4:26)]
condition = c(rep("Input", 1), rep("Output", 22))
geneID = as.character(barseq_diff$GeneID)
location=barseq_diff$pos
pool = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
test <- TnseqDiff(countData, geneID, location, pool, condition)

##################
###PCA Analysis###
##################

###notes
#test1 <- test[ , which(apply(test, 2, var) != 0)]
#test1[] <- lapply(test1, function(x) as.numeric(as.character(x)))
#class(test1$IT002)
#pca <- prcomp(t(test1), scale = TRUE)

#Sum all gene values otherwise prcomp fails on redundant gene names
barseq.agg <- aggregate(cbind(IT001, IT002, IT003, IT004, IT005, IT006, IT007, IT008, IT009, IT010, IT011, IT012, IT013, IT014, IT015, IT016, IT017, IT018, IT019, IT020, IT021, IT022,IT023, IT024,IT025) ~ GeneID, data = barseq_diff, FUN = sum)

#Gene names must be row names
row.names(barseq.agg) = barseq.agg$GeneID
barseq.agg <- barseq.agg[,c(-1)] #remove GeneID col
#barseq.agg <- barseq.agg %>% select(-GeneID)


barseq.agg$E_mean <- rowMeans(subset(barseq.agg, select = c(IT002, IT003, IT004, IT005, IT006)), na.rm = TRUE)
barseq.agg$ES_mean <- rowMeans(subset(barseq.agg, select = c(IT007, IT008, IT009, IT010, IT011)), na.rm = TRUE)
barseq.agg$EM_mean <- rowMeans(subset(barseq.agg, select = c(IT012, IT013, IT014, IT015, IT016)), na.rm = TRUE)
barseq.agg$ESM_mean <- rowMeans(subset(barseq.agg, select = c(IT017, IT018, IT019, IT020, IT021)), na.rm = TRUE)

#if necessary filter out any rows with all 0s other prcomp will fail
barseq.agg <- barseq.agg[rowSums((barseq.agg[,1:24])) > 0,]

barseq_means <- cbind(barseq.agg[,c(1, 26:29)])

#pca <- prcomp(t(subset(bar_pool1, select = -c(IT001))), scale = TRUE)

pca <- prcomp(t(barseq_means), scale = TRUE)

plot(pca$x[,1], pca$x[,2])

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab = "Principal Component", ylab = "Percent Variation")

pca.data <- data.frame(Sample = rownames(pca$x), X = pca$x[,1], Y = pca$x[,2])

ggplot(data = pca.data, aes(x=X, y = Y, label = Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep ="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_bw()


#######################
##Graph by Gene vs T0##
#######################

gene_fq <- barseq.agg %>%
  +     gather(index, counts, -GeneID, -IT001)

gene_fq %>%
  +     ggplot(aes(x = IT001, y = counts)) +
  +     facet_wrap(~ index) +
  +     geom_point()


#filter high T0 (x) point out
gene_fq_filtered <- subset(gene_fq, gene_fq[ ,2] < 206106)
gene_fq_filtered %>%
  +     ggplot(aes(x = IT001, y = counts)) +
  +     facet_wrap(~ index) +
  +     geom_point()
