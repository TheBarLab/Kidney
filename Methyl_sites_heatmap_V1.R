# read excel table and show heatmap

library("readxl")
library("ggplot2")
library("gplots")
library("RColorBrewer")

setwd("/home/naorsagy/Desktop/R_Methyl")
my_data <- read_excel("multiGene2.xlsx")

data_dims=dim(my_data)
hm_dat=my_data[2:data_dims[1],10:data_dims[2]]
hm_dat_mat=as.matrix(hm_dat)
mat_num <- matrix(as.numeric(hm_dat_mat),    # Convert to numeric matrix
                  ncol = ncol(hm_dat_mat))

#jpeg(file="filename.jpg")
col_breaks <- seq(-3,3, by = 0.1)
colormaphm<- colorRampPalette(c("blue", "white", "yellow"))(n=299)
my_palette <- c(colorRampPalette(c("blue", "white"))(29), 
                colorRampPalette(c("white","white"))(2),
                colorRampPalette(c("white", "yellow"))(29))

                #colorRampPalette(c("gold", "red"))(length(col_breaks)-31))
# (optional) defines the color breaks manually for a "skewed" color transition


#heatmap.2(mat_num,Rowv=NA,Colv=NA,trace="none",col="bluered",scale="row",xlab="samples",ylab="genes")
heatmap.2(mat_num,Rowv=NA,Colv=NA,trace="none",col=my_palette, breaks=col_breaks,scale="row",xlab="samples",ylab="genes")
#dev.off()
#heatmap(hm_dat_mat,Rowv=NA)


# heatmap and pca for incidences of seqs in methylation sites

my_data <- read_excel("20220906_buckets.xlsx")

