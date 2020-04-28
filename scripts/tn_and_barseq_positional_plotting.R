#Plotting TnSeq map and positional infromation

library(ggplot2)

tnmap_df <- read.csv("~/Documents/PostDoc_Harcombe/TnSeq/BarSeq/data/g/Eo221/map_TnSeq/random_pool_files/KAPA_Eolib_NZ_CP009273_pool", sep = "\t")
barseq_pool <- read.csv("~/Documents/PostDoc_Harcombe/TnSeq/BarSeq/data/g/Eo221/KAPA_EoTn_NZ_CP009273/combined_BarSeq_data.poolcount", sep ="\t")

#plot tnseq map pos vs counts and color peak region with either 2749518 or 2727475.  
#The latter looks better for the start of the peak
tnmap_df %>% 
  mutate(Color = ifelse(pos >= 2727475, "black", "red")) %>% 
  ggplot(aes(x = pos, y = n, color = Color)) + 
  geom_point(size = 0.1, alpha = 0.1) + 
  scale_color_identity()

#plot a single barseq positional map
#could make for loop for this
ggplot(data = subset(barseq_pool, barseq_pool$IT022 < 100 & barseq_pool$IT022 > 5), mapping = aes(x = pos, y = IT022)) + 
  geom_point(size = 0.1, alpha = 0.1)

#widen barseq poolcount data
barseq_wide <- barseq_pool %>% gather(index, counts, -barcode, -rcbarcode, -scaffold, -strand, -pos, -IT001)

#filter out counts and plot multiple plots
barseq_wide <- subset(barseq_wide, counts < 40000)
barseq_wide <- subset(barseq_wide, counts > 4)
barseq_wide %>% 
  mutate(Color = ifelse(pos >= 2727475, "black", "red")) %>% 
  ggplot(aes(x = IT001, y = counts, color = Color)) + 
  facet_wrap(~ index) + geom_point() + 
  scale_color_identity()
