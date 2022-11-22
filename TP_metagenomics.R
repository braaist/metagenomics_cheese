library(ggplot2)
install.packages("vegan")
library(ggpubr)


#Download package
phylum_matrix <- read.csv("/Users/braaist/Documents/paris_metagenomics/mgRAST_mgp3362.txt", sep = "\t")
phylum_matrix<-as.matrix(phylum_matrix)
View(phylum_matrix)

#calculate OTU statistics
init_vec = c()
init_length = c()
for (col in colnames(data.frame(phylum_matrix))){
  init_vec = union(init_vec, rownames(data.frame(phylum_matrix[,col][phylum_matrix[,col] != 0])))
  init_length = c(init_length, length(init_vec))
}


#Plotting
out_df <- data.frame(samples=c(1:dim(phylum_matrix)[2]),
                 OTU=init_length)
ggplot(data=out_df, aes(x=samples, y=OTU, group=1)) +
  geom_line(color="red")+
  geom_point()


#Custom function for the Shannon Diversity
ShannonDiversity <- function(microbiome_vec){
  microbiome_vec <- microbiome_vec[microbiome_vec!=0]
  p <- microbiome_vec/sum(microbiome_vec)
  microbiome_vec[microbiome_vec!=0]
  SI <- -sum(p*log(p))
  return(SI)
}

ShannonDiversity(phylum_matrix[3,])
phylum_matrix = data.frame(phylum_matrix)

#Compare results of different functions
SD1 <- t(data.frame(lapply(phylum_matrix, diversity)))
SD2 <- t(data.frame(lapply(phylum_matrix, ShannonDiversity)))
SD1 == SD2


#Check the associations of factors with qualitative variables
metadata = read.csv("/Users/braaist/Documents/paris_metagenomics/MetadataCheese.csv", sep = "\t")
View(metadata)

metadata_quantitative <- metadata[,colnames(metadata) %in% c("Moisture", "pH", "NaCl")] 
metadata_rest <- metadata[,!colnames(metadata) %in% c("Moisture", "pH", "NaCl")] 

df_quantitative = data.frame(metadata_quantitative, SD[,1])
colnames(df_quantitative) = c("Moisture", "pH", "NaCl", "SD")

lmHeight = lm(SD~., data = df_quantitative) 
summary(lmHeight)

ggplot(df_quantitative, aes(x = SD, y = Moisture)) +
  geom_point() +
  stat_smooth(method = lm)


#Check the associations of factors with quantitative variables
df_rest = data.frame(metadata_rest, SD[,1])
colnames(df_rest) = c(colnames(df_rest)[-length(colnames(df_rest))], "SD")

my_comparisons <- list( c("natural", "washed"), c("natural", "bloomy"), c("washed", "bloomy") )
p <- ggboxplot(metadata_rest, x = "RindType", y = "SD",
               color = "RindType", palette = "jco",
               add = "jitter") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test")
p

my_comparisons <- list( c("Y", "N"))
p <- ggboxplot(metadata_rest, x = "Pasteurized", y = "SD",
               color = "Pasteurized", palette = "jco",
               add = "jitter") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test")
p


