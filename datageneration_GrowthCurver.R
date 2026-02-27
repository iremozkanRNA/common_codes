rm(list = ls())
lapply(c("growthcurver", "dplyr", "ggplot2"), library, character.only = T)

setwd("~/Desktop/Mutant/IOAnalysis/Vanco/")
d <- read.csv("Vanco_0.25ugml_GrowthCurver.csv", header = TRUE, stringsAsFactors = FALSE)
#For whatever reason, sometimes it doesn't read the time as numeric
d$time <- as.numeric(d$time)
#Change this based on data fit
d <- d[d$time >= 1 ,]
d <- d[d$time <= 15 ,]

gc_out <- SummarizeGrowthByPlate(d, plot_fit = T, plot_file = "GrowthCurverResults2/Vanco_0.25ugml_growth_summary.pdf")
hist(gc_out$sigma, main = "Histogram of sigma values", xlab = "sigma")

#plot(gc_out)
# Change file names as necessary
write.csv(gc_out, "GrowthCurverResults2/Vanco_0.25ugml_growth_summary.csv", quote=FALSE, row.names=FALSE)

gc_out %>% top_n(5, sigma) %>% arrange(desc(sigma))

pca_gc_out <- as_data_frame(gc_out) 

# Prepare the gc_out data for the PCA
rownames(pca_gc_out) <- pca_gc_out$sample
# Do the PCA
pca.res <- prcomp(pca_gc_out %>% select(k:sigma), center=TRUE, scale=TRUE)

# Plot the results
as_data_frame(list(PC1=pca.res$x[,1],
                   PC2=pca.res$x[,2],
                   samples = pca_gc_out$sample)) %>% 
  ggplot(aes(x=PC1,y=PC2, label=samples)) + 
  geom_text(size = 3)

# Create table (assign to object first)
pca_scores <- as_data_frame(list(PC1=pca.res$x[,1],
                                 PC2=pca.res$x[,2],
                                 samples = pca_gc_out$sample))

# Your exact same plotting code
pca_scores %>% 
  ggplot(aes(x=PC1,y=PC2, label=samples)) + 
  geom_text(size = 3)

write.csv(pca_scores, file = "GrowthCurverResults2/Vanco_0.25ugml_pcascores.csv", col.names = T)
