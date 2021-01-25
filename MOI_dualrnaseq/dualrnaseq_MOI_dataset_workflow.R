############################################
#
#Dual RNA-Seq MOI dataset
#Regan Hayward - PhD
#
############################################

#Data description
#RNA-seq from Host (human) and bacteria (chlamydia trachomatis)
#Two time points: 1hr and 24hrs
#Three different MOI's: 0.1, 1 and 10
#Two different rRNA depletion methods: rRNA dep and rRNA dep + PolyA combined

#Pre-import steps
#Reads QC'ed and mapped with counts generated using feature counts




#---------------------------
#Install and load all required libraries
#---------------------------
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")

#Are packages already installed
list.of.packages <- c("easypackages","edgeR","scales","enrichR","rtracklayer",
                      "org.Hs.eg.db","PCAtools","genefilter","reshape2","Cairo","RColorBrewer",
                      "ggplot2","EDASeq","dendextend","EnsDb.Hsapiens.v86","dplyr","ggplot2",
                      "cowplot","pheatmap","beeswarm")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

#Load all required packages using easypackages
easypackages::libraries(list.of.packages)






#**************************************
#
# Load in the counts, filter and rename ----
#
#**************************************

#---
#load in
#---

#Chlamydial
setwd("Chlamydia")
counts_featureCounts_chlamyd <- readDGE(list.files(pattern = ".count"), columns = c(1,2))
counts_featureCounts_chlamyd$samples

#Human
setwd("Human")
counts_featureCounts_human <- readDGE(list.files(pattern = ".count"), columns = c(1,2))
counts_featureCounts_human$samples


#---
#Filter
#---

#Chlamydia
keep <- rowSums(cpm(counts_featureCounts_chlamyd)>1) >= 3 #filter
counts_featureCounts_chlamyd_v2 <- counts_featureCounts_chlamyd[keep,]
dim(counts_featureCounts_chlamyd);dim(counts_featureCounts_chlamyd_v2) #comparing final counts

#Human
keep <- rowSums(cpm(counts_featureCounts_human)>1) >= 3 #filter
counts_featureCounts_human_v2 <- counts_featureCounts_human[keep,]
dim(counts_featureCounts_human);dim(counts_featureCounts_human_v2) #comparing final counts



#---
#Rename the files names to be more readable
#---

a <- counts_featureCounts_chlamyd_v2$samples$files
counts_featureCounts_chlamyd_v2$samples$files <- substr(a,1,nchar(a)-6)

a <- counts_featureCounts_human_v2$samples$files
counts_featureCounts_human_v2$samples$files <- substr(a,1,nchar(a)-6)




#**************************************
#
# Library size plots ----
#
#**************************************

#---
#Prep colours
#---

plot_cols <- rep(brewer.pal(6,"Set3"),each=6)

#--
#Save plots
#--
width=10;height=5;dpi=300

#Split into 1 and 24 hours to get a better visual representation (as read #'s vary lots)

#Human 1hr
height=5;width=5
Cairo(file="library_size_human_1hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
ggplot(counts_featureCounts_human$samples[1:18,], aes(x = sort(files, decreasing = F), y=lib.size/10^6)) + #re-order values, small to large
  geom_bar(width=0.9, stat = "identity", fill=plot_cols[1:18] )+ #stacked and re-ordered
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7)) + #make labels vertical
  scale_y_continuous(name="Library size (millions)",expand = c(0, 0), limits = c(0, 70)) + #remove x-axis gap 
  theme(legend.title=element_blank()) +#remove title from legend
  scale_fill_manual(values=cbPalette) +#colour palette
  ggtitle("Human 1hr") +
  theme(plot.title = element_text(hjust = 0.5)) + #center title
  scale_x_discrete(labels=substr(counts_featureCounts_human_v2$samples$files[1:18],1,nchar(counts_featureCounts_human_v2$samples$files[1:18])-24))
dev.off()


#Human 24hr
Cairo(file="library_size_human_24hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
ggplot(counts_featureCounts_human$samples[19:36,], aes(x = sort(files, decreasing = F), y=lib.size/10^6)) + #re-order values, small to large
  geom_bar(width=0.9, stat = "identity", fill=plot_cols[1:18] )+ #stacked and re-ordered
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7)) + #make labels vertical
  scale_y_continuous(name="Library size (millions)",expand = c(0, 0), limits = c(0, 10)) + #remove x-axis gap 
  theme(legend.title=element_blank()) +#remove title from legend
  scale_fill_manual(values=cbPalette) +#colour palette
  ggtitle("Human 24hr") +
  theme(plot.title = element_text(hjust = 0.5)) + #center title
  scale_x_discrete(labels=substr(counts_featureCounts_human_v2$samples$files[19:36],1,nchar(counts_featureCounts_human_v2$samples$files[19:36])-24))
dev.off()


#CT 1hr
Cairo(file="library_size_CT_1hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
ggplot(counts_featureCounts_chlamyd$samples[1:18,], aes(x = sort(files, decreasing = F), y=lib.size)) + #re-order values, small to large
  geom_bar(width=0.9, stat = "identity", fill=plot_cols[1:18] )+ #stacked and re-ordered
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7)) + #make labels vertical
  scale_y_continuous(name="Library size",expand = c(0, 0), limits = c(0, 100000)) + #remove x-axis gap 
  theme(legend.title=element_blank()) +#remove title from legend
  scale_fill_manual(values=cbPalette) +#colour palette
  ggtitle("CT 1hr") +
  theme(plot.title = element_text(hjust = 0.5)) + #center title
  scale_x_discrete(labels=counts_featureCounts_chlamyd_v2$samples$files[1:18])
dev.off()

#CT 24hr
Cairo(file="library_size_CT_24hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
ggplot(counts_featureCounts_chlamyd$samples[19:36,], aes(x = sort(files, decreasing = F), y=lib.size/10^6)) + #re-order values, small to large
  geom_bar(width=0.9, stat = "identity", fill=plot_cols[1:18] )+ #stacked and re-ordered
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7)) + #make labels vertical
  scale_y_continuous(name="Library size (millions)",expand = c(0, 0), limits = c(0, 4)) + #remove x-axis gap 
  theme(legend.title=element_blank()) +#remove title from legend
  scale_fill_manual(values=cbPalette) +#colour palette
  ggtitle("CT 24hr") +
  theme(plot.title = element_text(hjust = 0.5)) + #center title
  scale_x_discrete(labels=counts_featureCounts_chlamyd_v2$samples$files[19:36])
dev.off()




#**************************************
#
# Plots showing the % of CT and Human reads ----
#
#**************************************

#--
#Prep colours and data structure
#--

#Prep colours
percent_plot_cols <- c("steelblue3","firebrick2")

#plotting % of human and CT and plasmid
total = counts_featureCounts_human_v2$samples$lib.size + counts_featureCounts_chlamyd_v2$samples$lib.size
human_chlamydia_df <- data.frame(Samples = counts_featureCounts_human_v2$samples$files,
                                 Human = counts_featureCounts_human_v2$samples$lib.size / total * 100,
                                 Chlamydia = counts_featureCounts_chlamyd_v2$samples$lib.size / total * 100)

#tidy up the names
human_chlamydia_df$Samples <- as.character(human_chlamydia_df$Samples)
human_chlamydia_df$Samples <- substr(human_chlamydia_df$Samples,1,nchar(human_chlamydia_df$Samples)-24)


#Melt data to be plot-ready 
human_chlamydia_df_melt <- melt(human_chlamydia_df, id='Samples')
human_chlamydia_df_melt$order <- 1:72
levels(human_chlamydia_df_melt$Samples) <- human_chlamydia_df_melt$order



#--
#Create plots
#--
width=10;height=5


Cairo(file="percent_mapped_human_CT.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
ggplot(data = human_chlamydia_df_melt, aes(x = Samples, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'stack', colour="white") +
  scale_fill_manual(values=percent_plot_cols)+ #colour palette
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7)) + #make labels vertical
  scale_y_continuous(name="% Mapped",expand = c(0, 0), limits = c(0, 101)) + #remove x-axis gap 
  ggtitle("Human and Chlamydial reads") +
  xlab("")+
  theme(plot.title = element_text(hjust = 0.5)) +#center title
  theme(legend.position="bottom") + #place bottom
  theme(legend.text = element_text(colour="black", size=7)) + #size and colour of text
  theme(legend.key.size = unit(0.5,"line")) + #size of boxes
  theme(legend.title=element_blank()) +   #remove the legend name
  scale_x_discrete(limits=human_chlamydia_df_melt$Samples[1:36])
dev.off()



#**************************************
#
# Working out the % of bacterial reads from the different depletion methods ----
#
#**************************************

#Note:
#There are more elegant ways of doing this, such as a loop assigning dynamic variables


#1hr
#Testing 1 hr MOI 0.1 - rRNA Dep
subset = 1:3
total_reads_moi_01_rRNA <- sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) #total
sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset]))
((sum(counts_featureCounts_human_v2$samples$lib.size[subset])) / total_reads_moi_01_rRNA) *100
((sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) / total_reads_moi_01_rRNA) *100

#Testing 1 hr MOI 0.1 - rRNA Dep + polyA dep
subset = 4:6
total_reads_moi_01_rRNA_PolyA <- sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) #total
sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset]))
((sum(counts_featureCounts_human_v2$samples$lib.size[subset])) / total_reads_moi_01_rRNA_PolyA) *100
((sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) / total_reads_moi_01_rRNA_PolyA) *100

#Testing 1 hr MOI 1 - rRNA Dep
subset = 7:9
total_reads_moi_1_rRNA <- sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) #total
sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset]))
((sum(counts_featureCounts_human_v2$samples$lib.size[subset])) / total_reads_moi_1_rRNA) *100
((sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) / total_reads_moi_1_rRNA) *100

#Testing 1 hr MOI 1 - rRNA Dep + polyA dep
subset = 10:12
total_reads_moi_1_rRNA_PolyA <- sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) #total
sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset]))
((sum(counts_featureCounts_human_v2$samples$lib.size[subset])) / total_reads_moi_1_rRNA_PolyA) *100
((sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) / total_reads_moi_1_rRNA_PolyA) *100

#Testing 1 hr MOI 10 - rRNA Dep
subset = 13:15
total_reads_moi_10_rRNA <- sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) #total
sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset]))
((sum(counts_featureCounts_human_v2$samples$lib.size[subset])) / total_reads_moi_10_rRNA) *100
((sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) / total_reads_moi_10_rRNA) *100

#Testing 1 hr MOI 10 - rRNA Dep + polyA dep
subset = 16:18
total_reads_moi_10_rRNA_PolyA <- sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) #total
sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset]))
((sum(counts_featureCounts_human_v2$samples$lib.size[subset])) / total_reads_moi_10_rRNA_PolyA) *100
((sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) / total_reads_moi_10_rRNA_PolyA) *100


#24hr
#Testing 24 hr MOI 0.1 - rRNA Dep
subset = 19:21
total_reads_moi_24_01_rRNA <- sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) #total
sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])/3)
((sum(counts_featureCounts_human_v2$samples$lib.size[subset])) / total_reads_moi_24_01_rRNA) *100
((sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) / total_reads_moi_24_01_rRNA) *100

#Testing 24 hr MOI 0.1 - rRNA Dep + polyA dep
subset = 22:24
total_reads_moi_24_01_rRNA_PolyA <- sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) #total
sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset]))
((sum(counts_featureCounts_human_v2$samples$lib.size[subset])) / total_reads_moi_24_01_rRNA_PolyA) *100
((sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) / total_reads_moi_24_01_rRNA_PolyA) *100

#Testing 24 hr MOI 1 - rRNA Dep
subset = 25:27
total_reads_moi_24_1_rRNA <- sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) #total
sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset]))
((sum(counts_featureCounts_human_v2$samples$lib.size[subset])) / total_reads_moi_24_1_rRNA) *100
((sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) / total_reads_moi_24_1_rRNA) *100

#Testing 24 hr MOI 1 - rRNA Dep + polyA dep
subset = 28:30
total_reads_moi_24_1_rRNA_PolyA <- sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) #total
sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset]))
((sum(counts_featureCounts_human_v2$samples$lib.size[subset])) / total_reads_moi_24_1_rRNA_PolyA) *100
((sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) / total_reads_moi_24_1_rRNA_PolyA) *100

#Testing 24 hr MOI 10 - rRNA Dep
subset = 31:33
total_reads_moi_24_10_rRNA <- sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) #total
sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset]))
((sum(counts_featureCounts_human_v2$samples$lib.size[subset])) / total_reads_moi_24_10_rRNA) *100
((sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) / total_reads_moi_24_10_rRNA) *100

#Testing 24 hr MOI 10 - rRNA Dep + polyA dep
subset = 34:36
total_reads_moi_24_10_rRNA_PolyA <- sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) #total
sum(sum(counts_featureCounts_human_v2$samples$lib.size[subset]) + sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset]))
((sum(counts_featureCounts_human_v2$samples$lib.size[subset])) / total_reads_moi_24_10_rRNA_PolyA) *100
((sum(counts_featureCounts_chlamyd_v2$samples$lib.size[subset])) / total_reads_moi_24_10_rRNA_PolyA) *100








#**************************************
#
# PCA plots (Before normalisation) ----
#
#**************************************

#--
#Set colours and shapes
#--

#first 18 are 1 and hollow
#second 18 are 24 and filled
plot.colors = c(rep(plot_cols[1], 6), rep(plot_cols[10],6),
                rep(plot_cols[13],6), rep(plot_cols[21], 6),
                rep(plot_cols[28], 6), rep(plot_cols[33], 6))

legend.colors = c(rep(plot_cols[1], 2), rep(plot_cols[10],2),
                  rep(plot_cols[13],2), rep(plot_cols[21], 2),
                  rep(plot_cols[28], 2), rep(plot_cols[33], 2))


#--
#Save plots
#--
width=8;height=8;dpi=300


#The split up plots - as thats what will be using downstream as the variation is too great between time points and MOI's

#Create function
pca_plots_before_norm <- function(filename, input_data, plot_cols, main, legend_names, legend_cols){
  Cairo(file=filename, type="png", units="in", width=8, height=8, dpi=300, bg = "white")
  par(xpd=T, mar=par()$mar+c(0,0,0,7)) # Expand right side of clipping rect to make room for the legend
  plotPCA(input_data, col=plot_cols, 
          pch = c(rep(1,3),rep(16,3),rep(1,3),rep(16,3),rep(1,3),rep(16,3)), 
          lwd=2 ,cex=2, labels=F, xlim=c(-0.5,0.5), main=main)
  legend(0.55,0.4, ncol = 1, bty = 'L', legend=legend_names, col=legend_cols, pch = c(1,16,1,16,1,16), cex = 0.9, horiz = F)
  dev.off()  
}

#Human - 1hr
pca_plots_before_norm(filename = "PCA_human_before_norm_1hr.png", input_data = counts_featureCounts_human_v2$counts[,c(1:18)], plot_cols = plot.colors,
                      main = "Human reads 1hr - before normalisation", 
                      legend_names = c("1hr 0.1 rRNAdep", "1hr 0.1 rRNApolyA", "1hr 1 rRNAdep", "1hr 1 rRNApolyA", "1hr 10 rRNAdep", "1hr 10 rRNApolyA"),
                      legend_cols = unique(plot.colors)[c(1,1,2,2,3,3)])

#Human - 24hr
pca_plots_before_norm(filename = "PCA_human_before_norm_24hr.png", input_data = counts_featureCounts_human_v2$counts[,c(19:36)], plot_cols = plot.colors,
                      main = "Human reads 24hr - before normalisation", 
                      legend_names = c("24hr 0.1 rRNAdep", "24hr 0.1 rRNApolyA", "24hr 1 rRNAdep", "24hr 1 rRNApolyA", "24hr 10 rRNAdep", "24hr 10 rRNApolyA"),
                      legend_cols = unique(plot.colors)[c(1,1,2,2,3,3)])




#CT - 1hr
pca_plots_before_norm(filename = "PCA_CT_before_norm_1hr.png", input_data = counts_featureCounts_chlamyd_v2$counts[,c(1:18)], plot_cols = plot.colors,
                      main = "CT reads 1hr - before normalisation", 
                      legend_names = c("1hr 0.1 rRNAdep", "1hr 0.1 rRNApolyA", "1hr 1 rRNAdep", "1hr 1 rRNApolyA", "1hr 10 rRNAdep", "1hr 10 rRNApolyA"),
                      legend_cols = unique(plot.colors)[c(1,1,2,2,3,3)])

#CT - 24hr
pca_plots_before_norm(filename = "PCA_CT_before_norm_24hr.png", input_data = counts_featureCounts_chlamyd_v2$counts[,c(19:36)], plot_cols = plot.colors,
                      main = "CT reads 24hr - before normalisation", 
                      legend_names = c("24hr 0.1 rRNAdep", "24hr 0.1 rRNApolyA", "24hr 1 rRNAdep", "24hr 1 rRNApolyA", "24hr 10 rRNAdep", "24hr 10 rRNApolyA"),
                      legend_cols = unique(plot.colors)[c(1,1,2,2,3,3)])






#**************************************
#
# RLE Plots (before norm) ----
#
#**************************************

#Plot params
width=10;height=5;dpi=300


#Tidy up the human names
human_rle_plots <- counts_featureCounts_human_v2$counts
colnames(human_rle_plots) <- substr(colnames(human_rle_plots),1,nchar(colnames(human_rle_plots))-24)

#Copy the CT names into the same var format (to be consistent)
chlamydial_rle_plots <- counts_featureCounts_chlamyd_v2$counts

#--
#RLE plots
#--

#Human - 1hr 
Cairo(file="RLE_plot_human_before_norm_1hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
par(mar=c(11,3,3,1))
plotRLE(human_rle_plots[,1:18], ylim=c(-4,4), las=2, cex=0.8, cex.axis = 0.7, outline=FALSE, col=plot_cols, main="RLE Plot 1hr - Human reads")
dev.off()

#Human - 24hr 
Cairo(file="RLE_plot_human_before_norm_24hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
par(mar=c(11,3,3,1))
plotRLE(human_rle_plots[,19:36], ylim=c(-4,4), las=2, cex=0.8, cex.axis = 0.7, outline=FALSE, col=plot_cols, main="RLE Plot 24hr - Human reads")
dev.off()

#CT - 1hr 
Cairo(file="RLE_plot_CT_before_norm_1hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
par(mar=c(11,3,3,1))
plotRLE(chlamydial_rle_plots[,1:18], ylim=c(-4,4), las=2, cex=0.8, cex.axis = 0.7, outline=FALSE, col=plot_cols, main="RLE Plot 1hr - CT reads")
dev.off()

#CT - 24hr 
Cairo(file="RLE_plot_CT_before_norm_24hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
par(mar=c(11,3,3,1))
plotRLE(chlamydial_rle_plots[,19:36], ylim=c(-4,4), las=2, cex=0.8, cex.axis = 0.7, outline=FALSE, col=plot_cols, main="RLE Plot 24hr - CT reads")
dev.off()





#**************************************
#
# Splitting into 1 and 24 hrs and removing outliers ----
#
#**************************************

#----Split out into 1 and 24 hours---
#Split out human 1 and 24 (as the variation is huge, especially for CT)
#Don't have to split human, but would be weird splitting one and not the other though

#outliers (in the first round)
#Human - 1HPI_MOI_10_rRNApolyAdep_rep2
#Human - 24HPI_MOI_10_rRNApolyAdep_rep3
#CT - 24HPI_MOI_0_1_rRNAployAdep_rep2

#Then in the second round, had to remove another one (second round)
#CT - 24HPI_MOI_0_1_rRNAployAdep_rep1


#So coming back and removing when splitting the counts into time points (step 1 + 2)


#Split out human 1 and 24hrs
edger_matrix_1hr_human <- counts_featureCounts_human_v2[,c(1:16,18)]
edger_matrix_24hr_human <- counts_featureCounts_human_v2[,c(19:35)]

#Split out CT 1 and 24
edger_matrix_1hr_CT <- counts_featureCounts_chlamyd_v2[,c(1:18)]
edger_matrix_24hr_CT <- counts_featureCounts_chlamyd_v2[,c(19:21,24:36)]


#**************************************
#
# Filtering again after splitting into 1 and 24 hrs and removing outliers ----
#
#**************************************

#using genefilter
#Filter based on the rows again (will help to make DE genes more accurate)

#function that implies our "expression measure is above 50 in at least 3 samples"
f1 <- kOverA(3, 50)
ffun <- filterfun(f1)

#Human
keep <- genefilter(edger_matrix_1hr_human, ffun)
edger_matrix_1hr_human_v2 <- edger_matrix_1hr_human[keep,]
dim(edger_matrix_1hr_human);dim(edger_matrix_1hr_human_v2)

keep <- genefilter(edger_matrix_24hr_human, ffun)
edger_matrix_24hr_human_v2 <- edger_matrix_24hr_human[keep,]
dim(edger_matrix_24hr_human);dim(edger_matrix_24hr_human_v2)

#function that implies our "expression measure is above 10 in at least 3 samples"
f1 <- kOverA(3, 10)
ffun <- filterfun(f1)

#CT
keep <- genefilter(edger_matrix_1hr_CT, ffun)
edger_matrix_1hr_CT_v2 <- edger_matrix_1hr_CT[keep,]
dim(edger_matrix_1hr_CT);dim(edger_matrix_1hr_CT_v2)

keep <- genefilter(edger_matrix_24hr_CT, ffun)
edger_matrix_24hr_CT_v2 <- edger_matrix_24hr_CT[keep,]
dim(edger_matrix_24hr_CT);dim(edger_matrix_24hr_CT_v2)


#**************************************
#
# Normalisation using edgeR ----
#
#**************************************


#--
#Create a df with metadata
#--

#Human
human_1_df <- data.frame("Fullname"=colnames(counts_featureCounts_human_v2[,c(1:18)]),
                         "RNA"=as.factor(rep(c(rep(c("rRNAdep","rRNA_polyA"),each=3)),3)),
                         "Treatment"=as.factor(rep(c("0.1","1","10"),each=6)))

human_24_df <- data.frame("Fullname"=colnames(counts_featureCounts_human_v2[,c(19:36)]),
                          "RNA"=as.factor(rep(c(rep(c("rRNAdep","rRNA_polyA"),each=3)),3)),
                          "Treatment"=as.factor(rep(c("0.1","1","10"),each=6)))

#CT
CT_1_df <- data.frame("Fullname"=colnames(counts_featureCounts_chlamyd_v2[,c(1:18)]),
                      "RNA"=as.factor(rep(c(rep(c("rRNAdep","rRNA_polyA"),each=3)),3)),
                      "Treatment"=as.factor(rep(c("0.1","1","10"),each=6)))

CT_24_df <- data.frame("Fullname"=colnames(counts_featureCounts_chlamyd_v2[,c(19:36)]),
                       "RNA"=as.factor(rep(c(rep(c("rRNAdep","rRNA_polyA"),each=3)),3)),
                       "Treatment"=as.factor(rep(c("0.1","1","10"),each=6)))


#Remove outlier rows
human_1_df <- human_1_df[-17,]
human_24_df <- human_24_df[-18,]

CT_1_df #not changed
CT_24_df <- CT_24_df[-5,]
CT_24_df <- CT_24_df[-4,]



#--
#Create the matrix
#--

#Human 1hr
the_design_human_1hr <- model.matrix(~human_1_df$RNA+human_1_df$Treatment)
rownames(the_design_human_1hr) <- human_1_df$Fullname
the_design_human_1hr
#Human 24hr
the_design_human_24hr <- model.matrix(~human_24_df$RNA+human_24_df$Treatment)
rownames(the_design_human_24hr) <- human_24_df$Fullname
the_design_human_24hr

#CT 1hr
the_design_CT_1hr <- model.matrix(~CT_1_df$RNA+CT_1_df$Treatment)
rownames(the_design_CT_1hr) <- CT_1_df$Fullname
the_design_CT_1hr
#CT 24hr
the_design_CT_24hr <- model.matrix(~CT_24_df$RNA+CT_24_df$Treatment)
rownames(the_design_CT_24hr) <- CT_24_df$Fullname
the_design_CT_24hr



#--
#Set up the counts into a var - doing this to recalculate the libray sizes after filtering
#--

#human
y1_human_1hr <- DGEList(edger_matrix_1hr_human_v2)
y1_human_24hr <- DGEList(edger_matrix_24hr_human_v2)
#CT
y1_CT_1hr <- DGEList(edger_matrix_1hr_CT_v2)
y1_CT_24hr <- DGEList(edger_matrix_24hr_CT_v2)


#--
#Normalisation
#--
#Note
#See version 30 for the different nomalisation methods and their settings, RLE, UQ etc

#Human
y2_human_1hr <- calcNormFactors(y1_human_1hr, method = "TMM") 
y2_human_24hr <- calcNormFactors(y1_human_24hr, method = "TMM")

#CT
y2_CT_1hr <- calcNormFactors(y1_CT_1hr, method = "TMM") #have tried all the methods including edger TMMwzp, also adjusting all the TMM params - no luck
y2_CT_24hr <- calcNormFactors(y1_CT_24hr, method = "TMM")



#--
#Estimate dispersions
#--

#Human
y3_human_1hr <-  estimateGLMCommonDisp(y2_human_1hr, the_design_human_1hr, verbose = TRUE)
y3_human_24hr <-  estimateGLMCommonDisp(y2_human_24hr, the_design_human_24hr, verbose = TRUE)


#CT
y3_CT_1hr <-  estimateGLMCommonDisp(y2_CT_1hr, the_design_CT_1hr, verbose = TRUE)
y3_CT_24hr <-  estimateGLMCommonDisp(y2_CT_24hr, the_design_CT_24hr, verbose = TRUE)



#--
#(BCV) Biological coefficient of variation numbers
#--
#around 0.1 is good (but normally around 0.15)

#Human
sqrt(y3_human_1hr$common.dispersion) 
sqrt(y3_human_24hr$common.dispersion)
#CT
sqrt(y3_CT_1hr$common.dispersion)
sqrt(y3_CT_24hr$common.dispersion)


#--
#BCV plots
#--

#Human
plotBCV(y3_human_1hr, main="edgeR - BCV plot - Human - 1hr")
plotBCV(y3_human_24hr, main="edgeR - BCV plot - Human - 24hr")
#CT
plotBCV(y3_CT_1hr, main="edgeR - BCV plot - CT - 1hr")
plotBCV(y3_CT_24hr, main="edgeR - BCV plot - CT - 24hr")


#--
#Fit the model
#--

#Human
fit_human_1hr <-  glmFit(y3_human_1hr, the_design_human_1hr)
fit_human_24hr <-  glmFit(y3_human_24hr, the_design_human_24hr)
#CT
fit_CT_1hr <-  glmFit(y3_CT_1hr, the_design_CT_1hr)
fit_CT_24hr <-  glmFit(y3_CT_24hr, the_design_CT_24hr)


#--
#Determine contrasts
#--

#Human
#1hr
lrt_human_1hr_treatment_1_vs_01 <- glmLRT(fit_human_1hr, coef = 3) 
lrt_human_1hr_treatment_10_vs_01 <- glmLRT(fit_human_1hr, coef = 4)
lrt_human_1hr_treatment_10_vs_1 <- glmLRT(fit_human_1hr, contrast=c(0,0,-1,1))
lrt_human_1hr_treatment_1_vs_10 <- glmLRT(fit_human_1hr, contrast=c(0,0,1,-1))
#24hr
lrt_human_24hr_treatment_1_vs_01 <- glmLRT(fit_human_24hr, coef = 3) 
lrt_human_24hr_treatment_10_vs_01 <- glmLRT(fit_human_24hr, coef = 4)
lrt_human_24hr_treatment_10_vs_1 <- glmLRT(fit_human_24hr, contrast=c(0,0,-1,1))
lrt_human_24hr_treatment_1_vs_10 <- glmLRT(fit_human_24hr, contrast=c(0,0,1,-1))

#CT
#1hr
lrt_CT_1hr_treatment_1_vs_01 <- glmLRT(fit_CT_1hr, coef = 3) 
lrt_CT_1hr_treatment_10_vs_01 <- glmLRT(fit_CT_1hr, coef = 4)
lrt_CT_1hr_treatment_10_vs_1 <- glmLRT(fit_CT_1hr, contrast=c(0,0,-1,1))
lrt_CT_1hr_treatment_1_vs_10 <- glmLRT(fit_CT_1hr, contrast=c(0,0,1,-1))
#24hr
lrt_CT_24hr_treatment_1_vs_01 <- glmLRT(fit_CT_24hr, coef = 3) 
lrt_CT_24hr_treatment_10_vs_01 <- glmLRT(fit_CT_24hr, coef = 4)
lrt_CT_24hr_treatment_10_vs_1 <- glmLRT(fit_CT_24hr, contrast=c(0,0,-1,1))
lrt_CT_24hr_treatment_1_vs_10 <- glmLRT(fit_CT_24hr, contrast=c(0,0,1,-1))




#--
#Count of top tags
#--

#Human
#1hr
dim(topTags(lrt_human_1hr_treatment_1_vs_01, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue"))
dim(topTags(lrt_human_1hr_treatment_1_vs_10, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue"))
#24hr
dim(topTags(lrt_human_24hr_treatment_1_vs_01, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue"))
dim(topTags(lrt_human_24hr_treatment_1_vs_10, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue"))

#CT
#1hr
dim(topTags(lrt_CT_1hr_treatment_1_vs_01, n=Inf, adjust.method = "fdr", sort.by = "PValue"))
dim(topTags(lrt_CT_1hr_treatment_1_vs_01, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue"))
dim(topTags(lrt_CT_1hr_treatment_1_vs_10, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue"))
#24hr
dim(topTags(lrt_CT_24hr_treatment_1_vs_01, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue"))
dim(topTags(lrt_CT_24hr_treatment_1_vs_10, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue"))




#--
#Top tags
#--

#Human
#1hr
DE_genes_human_1hr_1_vs_01 <- topTags(lrt_human_1hr_treatment_1_vs_01, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue")
DE_genes_human_1hr_10_vs_01 <- topTags(lrt_human_1hr_treatment_10_vs_01, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue")
DE_genes_human_1hr_10_vs_1 <- topTags(lrt_human_1hr_treatment_10_vs_1, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue")
DE_genes_human_1hr_1_vs_10 <- topTags(lrt_human_1hr_treatment_1_vs_10, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue")
#24hrs
DE_genes_human_24hr_1_vs_01 <- topTags(lrt_human_24hr_treatment_1_vs_01, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue")
DE_genes_human_24hr_10_vs_01 <- topTags(lrt_human_24hr_treatment_10_vs_01, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue")
DE_genes_human_24hr_10_vs_1 <- topTags(lrt_human_24hr_treatment_10_vs_1, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue")
DE_genes_human_24hr_1_vs_10 <- topTags(lrt_human_24hr_treatment_1_vs_10, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue")
#CT
#1hr
DE_genes_CT_1hr_1_vs_01 <- topTags(lrt_CT_1hr_treatment_1_vs_01, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue")
DE_genes_CT_1hr_10_vs_01 <- topTags(lrt_CT_1hr_treatment_10_vs_01, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue")
DE_genes_CT_1hr_10_vs_1 <- topTags(lrt_CT_1hr_treatment_10_vs_1, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue")
DE_genes_CT_1hr_1_vs_10 <- topTags(lrt_CT_1hr_treatment_1_vs_10, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue")
#24hr
DE_genes_CT_24hr_1_vs_01 <- topTags(lrt_CT_24hr_treatment_1_vs_01, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue")
DE_genes_CT_24hr_10_vs_01 <- topTags(lrt_CT_24hr_treatment_10_vs_01, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue")
DE_genes_CT_24hr_10_vs_1 <- topTags(lrt_CT_24hr_treatment_10_vs_1, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue")
DE_genes_CT_24hr_1_vs_10 <- topTags(lrt_CT_24hr_treatment_1_vs_10, p.value = 0.05, n=Inf, adjust.method = "fdr", sort.by = "PValue")



#--
#Save down the MD plots
#--
width=6;height=6

#Human
#1hr
Cairo(file="Human_MD_plot_TMM_1hr_moi_1_vs_01.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
plotMD(lrt_human_1hr_treatment_1_vs_01, main="Human - 1hr - MOI 1 vs 0.1 (TMM)") ; abline(h=c(-1,1), col="blue")
dev.off()

Cairo(file="Human_MD_plot_TMM_1hr_moi_1_vs_10.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
plotMD(lrt_human_1hr_treatment_1_vs_10, main="Human - 1hr - MOI 1 vs 10 (TMM)") ; abline(h=c(-1,1), col="blue")
dev.off()

#24hr
Cairo(file="Human_MD_plot_TMM_24hr_moi_1_vs_01.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
plotMD(lrt_human_24hr_treatment_1_vs_01, main="Human - 24hr - MOI 1 vs 0.1 (TMM)") ; abline(h=c(-1,1), col="blue")
dev.off()

Cairo(file="Human_MD_plot_TMM_24hr_moi_1_vs_10.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
plotMD(lrt_human_24hr_treatment_1_vs_10, main="Human - 24hr - MOI 1 vs 10 (TMM)") ; abline(h=c(-1,1), col="blue")
dev.off()


#CT
#1hr
Cairo(file="CT_MD_plot_TMM_1hr_moi_1_vs_01.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
plotMD(lrt_CT_1hr_treatment_1_vs_01, main="CT - 1hr - MOI 1 vs 0.1 (TMM)") ; abline(h=c(-1,1), col="blue")
dev.off()

Cairo(file="CT_MD_plot_TMM_1hr_moi_1_vs_10.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
plotMD(lrt_CT_1hr_treatment_1_vs_10, main="CT - 1hr - MOI 1 vs 10 (TMM)") ; abline(h=c(-1,1), col="blue")
dev.off()

#24hr
Cairo(file="CT_MD_plot_TMM_24hr_moi_1_vs_01.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
plotMD(lrt_CT_24hr_treatment_1_vs_01, main="CT - 24hr - MOI 1 vs 0.1 (TMM)") ; abline(h=c(-1,1), col="blue")
dev.off()

Cairo(file="CT_MD_plot_TMM_24hr_moi_1_vs_10.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
plotMD(lrt_CT_24hr_treatment_1_vs_10, main="CT - 24hr - MOI 1 vs 10 (TMM)") ; abline(h=c(-1,1), col="blue")
dev.off()




#**************************************
#
# Create beeswarm plots of the fold changes (just host) ----
#
#**************************************

   
#---
#Host
#---

#Make the dfs Human
#1hr
DE_genes_human_1hr_1_vs_01_bs_fc <- topTags(lrt_human_1hr_treatment_1_vs_01, n=Inf, adjust.method = "fdr", sort.by = "PValue")
DE_genes_human_1hr_10_vs_1_bs_fc <- topTags(lrt_human_1hr_treatment_10_vs_01, n=Inf, adjust.method = "fdr", sort.by = "PValue")
#24hr
DE_genes_human_24hr_1_vs_01_bs_fc <- topTags(lrt_human_24hr_treatment_1_vs_01, n=Inf, adjust.method = "fdr", sort.by = "PValue")
DE_genes_human_24hr_10_vs_1_bs_fc <- topTags(lrt_human_24hr_treatment_10_vs_01, n=Inf, adjust.method = "fdr", sort.by = "PValue")



width=5;height=3


#1hr - 1 vs 0.1 
DE_genes_human_1hr_1_vs_01_bs_fc_v2 <- DE_genes_human_1hr_1_vs_01_bs_fc$table[which(DE_genes_human_1hr_1_vs_01_bs_fc$table$FDR < 0.05),]
DE_genes_human_1hr_1_vs_01_bs_fc_v2$col <- "royalblue"
DE_genes_human_1hr_1_vs_01_bs_fc_v2$col[which(DE_genes_human_1hr_1_vs_01_bs_fc_v2$logFC > 0)] <- "firebrick1"
DE_genes_human_1hr_1_vs_01_bs_fc_v2$shape = 21
DE_genes_human_1hr_1_vs_01_bs_fc_v2$fill = grDevices::adjustcolor("royalblue", alpha.f = 0.5)
DE_genes_human_1hr_1_vs_01_bs_fc_v2$fill[which(DE_genes_human_1hr_1_vs_01_bs_fc_v2$logFC > 0)] <- grDevices::adjustcolor("firebrick1", alpha.f = 0.5)
DE_genes_human_1hr_1_vs_01_bs_fc_v2


Cairo(file="beeswarm_fc_human_1hrs_1_vs_01.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
beeswarm(DE_genes_human_1hr_1_vs_01_bs_fc_v2$logFC, horizontal = T, 
         pwcol = DE_genes_human_1hr_1_vs_01_bs_fc_v2$col, 
         pwpch = DE_genes_human_1hr_1_vs_01_bs_fc_v2$shape,
         pwbg = DE_genes_human_1hr_1_vs_01_bs_fc_v2$fill,
         spacing=0.45, cex=0.6, xlim = c(-3,5), method = 'swarm', cex.axis=1.2, main="Host - 1hr - 1 vs 0.1")
abline(v=c(-1,1), col="grey", lty=2)
dev.off()


#1hr - 10 vs 1
DE_genes_human_1hr_10_vs_1_bs_fc_v2 <- DE_genes_human_1hr_10_vs_1_bs_fc$table[which(DE_genes_human_1hr_10_vs_1_bs_fc$table$FDR < 0.05),]
DE_genes_human_1hr_10_vs_1_bs_fc_v2$col <- "royalblue"
DE_genes_human_1hr_10_vs_1_bs_fc_v2$col[which(DE_genes_human_1hr_10_vs_1_bs_fc_v2$logFC > 0)] <- "firebrick1"
DE_genes_human_1hr_10_vs_1_bs_fc_v2$shape = 21
DE_genes_human_1hr_10_vs_1_bs_fc_v2$fill = grDevices::adjustcolor("royalblue", alpha.f = 0.5)
DE_genes_human_1hr_10_vs_1_bs_fc_v2$fill[which(DE_genes_human_1hr_10_vs_1_bs_fc_v2$logFC > 0)] <- grDevices::adjustcolor("firebrick1", alpha.f = 0.5)
DE_genes_human_1hr_10_vs_1_bs_fc_v2
summary(DE_genes_human_1hr_10_vs_1_bs_fc_v2$logFC)

Cairo(file="beeswarm_fc_human_1hrs_10_vs_1.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
beeswarm(DE_genes_human_1hr_10_vs_1_bs_fc_v2$logFC, horizontal = T, 
         pwcol = DE_genes_human_1hr_10_vs_1_bs_fc_v2$col, 
         pwpch = DE_genes_human_1hr_10_vs_1_bs_fc_v2$shape, 
         pwbg = DE_genes_human_1hr_10_vs_1_bs_fc_v2$fill,
         spacing=0.2, cex=0.6, xlim = c(-8,10), method = 'swarm', cex.axis=1.2, main="Host - 1hr - 10 vs 1")
abline(v=c(-1,1), col="grey", lty=2)
dev.off()

#24hr - 1 vs 0.1 
DE_genes_human_24hr_1_vs_01_bs_fc_v2 <- DE_genes_human_24hr_1_vs_01_bs_fc$table[which(DE_genes_human_24hr_1_vs_01_bs_fc$table$FDR < 0.05),]
DE_genes_human_24hr_1_vs_01_bs_fc_v2$col <- "royalblue"
DE_genes_human_24hr_1_vs_01_bs_fc_v2$col[which(DE_genes_human_24hr_1_vs_01_bs_fc_v2$logFC > 0)] <- "firebrick1"
DE_genes_human_24hr_1_vs_01_bs_fc_v2$shape = 21
DE_genes_human_24hr_1_vs_01_bs_fc_v2$fill = grDevices::adjustcolor("royalblue", alpha.f = 0.5)
DE_genes_human_24hr_1_vs_01_bs_fc_v2$fill[which(DE_genes_human_24hr_1_vs_01_bs_fc_v2$logFC > 0)] <- grDevices::adjustcolor("firebrick1", alpha.f = 0.5)
DE_genes_human_24hr_1_vs_01_bs_fc_v2
summary(DE_genes_human_24hr_1_vs_01_bs_fc_v2$logFC)

Cairo(file="beeswarm_fc_human_24hrs_1_vs_01.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
beeswarm(DE_genes_human_24hr_1_vs_01_bs_fc_v2$logFC, horizontal = T, 
         pwcol = DE_genes_human_24hr_1_vs_01_bs_fc_v2$col, 
         pwpch = DE_genes_human_24hr_1_vs_01_bs_fc_v2$shape, 
         pwbg = DE_genes_human_24hr_1_vs_01_bs_fc_v2$fill, 
         spacing=0.7, cex=0.6, xlim = c(-2,3), method = 'swarm', cex.axis=1.2, main="Host - 24hr - 1 vs 0.1")
abline(v=c(-1,1), col="grey", lty=2)
dev.off()



#24hr - 10 vs 1
DE_genes_human_24hr_10_vs_1_bs_fc_v2 <- DE_genes_human_24hr_10_vs_1_bs_fc$table[which(DE_genes_human_24hr_10_vs_1_bs_fc$table$FDR < 0.05),]
DE_genes_human_24hr_10_vs_1_bs_fc_v2$col <- "royalblue"
DE_genes_human_24hr_10_vs_1_bs_fc_v2$col[which(DE_genes_human_24hr_10_vs_1_bs_fc_v2$logFC > 0)] <- "firebrick1"
DE_genes_human_24hr_10_vs_1_bs_fc_v2$shape = 21
DE_genes_human_24hr_10_vs_1_bs_fc_v2$fill = grDevices::adjustcolor("royalblue", alpha.f = 0.5)
DE_genes_human_24hr_10_vs_1_bs_fc_v2$fill[which(DE_genes_human_24hr_10_vs_1_bs_fc_v2$logFC > 0)] <- grDevices::adjustcolor("firebrick1", alpha.f = 0.5)
DE_genes_human_24hr_10_vs_1_bs_fc_v2
summary(DE_genes_human_24hr_10_vs_1_bs_fc_v2$logFC)

Cairo(file="beeswarm_fc_human_24hrs_10_vs_1.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
beeswarm(DE_genes_human_24hr_10_vs_1_bs_fc_v2$logFC, horizontal = T, 
         pwcol = DE_genes_human_24hr_10_vs_1_bs_fc_v2$col, 
         pwpch = DE_genes_human_24hr_10_vs_1_bs_fc_v2$shape,
         pwbg = DE_genes_human_24hr_10_vs_1_bs_fc_v2$fill,
         spacing=0.15, cex=0.6, xlim = c(-8,10), method = 'swarm', cex.axis=1.2, main="Host - 24hr - 10 vs 1")
abline(v=c(-1,1), col="grey", lty=2)
dev.off()




       




#**************************************
#
# Extracting out the normalised counts ----
#
#**************************************

#View the the normalised values
#human
cpm.DGEList(y2_human_1hr)[c(1:5),c(1:5)] 
cpm.DGEList(y2_human_24hr)[c(1:5),c(1:5)] 

#CT
cpm.DGEList(y2_CT_1hr)[c(1:5),c(1:5)] 
cpm.DGEList(y2_CT_24hr)[c(1:5),c(1:5)] 

#Extract out the normalised values
norm_counts_human_1hr_v1 <- cpm.DGEList(y2_human_1hr)
norm_counts_human_24hr_v1 <- cpm.DGEList(y2_human_24hr)

norm_counts_CT_1hr_v1 <- cpm.DGEList(y2_CT_1hr)
norm_counts_CT_24hr_v1 <- cpm.DGEList(y2_CT_24hr)




#**************************************
#
# PCA plots - after normalisation ----
#
#**************************************

#this could go in a function...

#Human - 1hr
Cairo(file="PCA_human_after_norm_TMM_1hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")

par(xpd=T, mar=par()$mar+c(0,0,0,7)) # Expand right side of clipping rect to make room for the legend
plotPCA(norm_counts_human_1hr_v1, col=unique(cbPalette)[c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6)], 
        pch = c(rep(16,3),rep(16,3),rep(16,3),rep(16,3),rep(16,3),rep(16,3)), 
        lwd=3 ,cex=3, labels=F, xlim=c(-0.5,0.5), main="Human reads 1hr - after norm. - (TMM)")
legend(0.55,0.4, ncol = 1, bty = 'L',
       legend=c("1hr 0.1 rRNAdep", "1hr 0.1 rRNApolyA", "1hr 1 rRNAdep", "1hr 1 rRNApolyA", "1hr 10 rRNAdep", "1hr 10 rRNApolyA"), 
       col=unique(cbPalette)[c(1:6)], pch = c(16,16,16,16,16,16), cex = 0.9, horiz = F)
dev.off()


#Human - 24hr
Cairo(file="PCA_human_after_norm_TMM_24hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
par(xpd=T, mar=par()$mar+c(0,0,0,7)) # Expand right side of clipping rect to make room for the legend
plotPCA(norm_counts_human_24hr_v1, col=unique(cbPalette)[c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6)], 
        pch = c(rep(16,3),rep(16,3),rep(16,3),rep(16,3),rep(16,3),rep(16,3)), 
        lwd=3 ,cex=3, labels=F, xlim=c(-0.5,0.5), main="Human reads 24hr - after norm. - (TMM)")
legend(0.55,0.25, ncol = 1, bty = 'L',
       legend=c("24hr 0.1 rRNAdep", "24hr 0.1 rRNApolyA", "24hr 1 rRNAdep", "24hr 1 rRNApolyA", "24hr 10 rRNAdep", "24hr 10 rRNApolyA"), 
       col=unique(cbPalette)[c(1:6)], pch = c(16,16,16,16,16,16), cex = 0.9, horiz = F)
dev.off()





#CT - 1hr
colnames(norm_counts_CT_1hr_v1)
Cairo(file="PCA_CT_after_norm_TMM_1hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
par(xpd=T, mar=par()$mar+c(0,0,0,7)) # Expand right side of clipping rect to make room for the legend
plotPCA(norm_counts_CT_1hr_v1, col=unique(cbPalette)[c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6)], 
        pch = c(rep(16,3),rep(16,3),rep(16,3),rep(16,3),rep(16,3),rep(16,3)), 
        lwd=3 ,cex=3, labels=F, xlim=c(-0.6,0.4), main="CT reads 1hr - after norm. - (TMM)")
legend(0.45,0.4, ncol = 1, bty = 'L',
       legend=c("1hr 0.1 rRNAdep", "1hr 0.1 rRNApolyA", "1hr 1 rRNAdep", "1hr 1 rRNApolyA", "1hr 10 rRNAdep", "1hr 10 rRNApolyA"), 
       col=unique(cbPalette)[c(1:6)], pch = c(16,16,16,16,16,16), cex = 0.9, horiz = F)
dev.off()

#CT - 24hr
colnames(norm_counts_CT_24hr_v1)
Cairo(file="PCA_CT_after_norm_TMM_24hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
par(xpd=T, mar=par()$mar+c(0,0,0,7)) # Expand right side of clipping rect to make room for the legend
plotPCA(norm_counts_CT_24hr_v1, col=unique(cbPalette)[c(1,1,1,2,3,3,3,4,4,4,5,5,5,6,6,6)], 
        pch = c(rep(16,3),rep(16,1),rep(16,3),rep(16,3),rep(16,3),rep(16,3)), 
        lwd=3 ,cex=3, labels=F, xlim=c(-0.4,0.8), main="CT reads 24hr - after norm. - (TMM)")
legend(0.86,0.43, ncol = 1, bty = 'L',
       legend=c("24hr 0.1 rRNAdep", "24hr 0.1 rRNApolyA", "24hr 1 rRNAdep", "24hr 1 rRNApolyA", "24hr 10 rRNAdep", "24hr 10 rRNApolyA"), 
       col=unique(cbPalette)[c(1:6)], pch = c(16,16,16,16,16,16), cex = 0.9, horiz = F)
dev.off()










#**************************************
#
# RLE plots - after normalisation ----
#
#**************************************

width=5;height=5;dpi=300

#---
#Human
#---

#Tidy up the names
human_rle_plots_norm_1hr <- norm_counts_human_1hr_v1
human_rle_plots_norm_24hr <- norm_counts_human_24hr_v1
colnames(human_rle_plots_norm_1hr) <- substr(colnames(human_rle_plots_norm_1hr),1,nchar(colnames(human_rle_plots_norm_1hr))-24)
colnames(human_rle_plots_norm_24hr) <- substr(colnames(human_rle_plots_norm_24hr),1,nchar(colnames(human_rle_plots_norm_24hr))-24)

#Human - 1hr 
Cairo(file="RLE_plot_human_after_norm_TMM_1hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
par(mar=c(11,3,3,1))
plotRLE(human_rle_plots_norm_1hr, ylim=c(-2,2), las=2, cex=0.8, cex.axis = 0.7, outline=FALSE, col=plot_cols[1:17], main="1hr - Human norm. (TMM)")
dev.off()

#Human - 24hr 
Cairo(file="RLE_plot_human_after_norm_TMM_24hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
par(mar=c(11,3,3,1))
plotRLE(human_rle_plots_norm_24hr, ylim=c(-2,2), las=2, cex=0.8, cex.axis = 0.7, outline=FALSE, col=plot_cols, main="24hr - Human norm. (TMM)")
dev.off()


#---
#CT
#---
#Keep the same variable naming scheme to be consistent for downstream
CT_rle_plots_norm_1hr <- norm_counts_CT_1hr_v1
CT_rle_plots_norm_24hr <- norm_counts_CT_24hr_v1
dim(norm_counts_CT_24hr_v1)
dim(CT_rle_plots_norm_24hr)

#CT - 1hr 
Cairo(file="RLE_plot_CT_after_norm_TMM_1hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
par(mar=c(11,3,3,1))
plotRLE(CT_rle_plots_norm_1hr, ylim=c(-10,8), las=2, cex=0.8, cex.axis = 0.7, outline=FALSE, col=plot_cols, main="TMM 1hr - CT norm. (TMM)")
dev.off()

#CT - 24hr 
colnames(CT_rle_plots_norm_24hr)
Cairo(file="RLE_plot_CT_after_norm_TMM_24hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
par(mar=c(11,3,3,1))
plotRLE(CT_rle_plots_norm_24hr, ylim=c(-1,1), las=2, cex=0.8, cex.axis = 0.7, outline=FALSE, col=
          unique(plot_cols)[c(rep(1,3),rep(2,1),rep(3,3),rep(4,3),rep(5,3),rep(6,3))], main="TMM 24hr - CT norm. (TMM)")
dev.off()






#**************************************
#
# Identifying outliers ----
#
#**************************************

#Note:
#Althogh the code is placed here, this is performed after normalisation, but DE genes are only determined once outliers have been removed 


#Three methods
#1) Looking at the eigen values from PC1 and PC2
#2) Looking at groupings from hierarchical clustering
#3) Looking at PCA plots using the normalised values and removing 10% of the lowest variance

#--
#1) Identifying if an outlier is > |3| sd's above the mean on PC1 and PC2
#--


#Calc and extract out the principal components
#Human
human_1hr_pcas <- prcomp(t(log2(human_rle_plots_norm_1hr +1)), scale. = F)
human_24hr_pcas <- prcomp(t(log2(human_rle_plots_norm_24hr + 1)), scale. = F)
#CT
CT_1hr_pcas <- prcomp(t(log2(CT_rle_plots_norm_1hr +1)), scale. = F)
CT_24hr_pcas <- prcomp(t(log2(CT_rle_plots_norm_24hr + 1)), scale. = F)


#Use the absolute values - get the mean, sd and calculate 3.sds
#Get mean
#human
human_1hr_mean <- mean(abs(as.numeric(human_1hr_pcas$x[,1])))
human_24hr_mean <- mean(abs(as.numeric(human_24hr_pcas$x[,1])))
#CT
CT_1hr_mean <- mean(abs(as.numeric(CT_1hr_pcas$x[,1])))
CT_24hr_mean <- mean(abs(as.numeric(CT_24hr_pcas$x[,1])))

#Get s.d
#Human
human_1hr_sd <- sd(abs(as.numeric(human_1hr_pcas$x[,1])))
human_24hr_sd <- sd(abs(as.numeric(human_24hr_pcas$x[,1])))
#CT
CT_1hr_sd <- sd(abs(as.numeric(CT_1hr_pcas$x[,1])))
CT_24hr_sd <- sd(abs(as.numeric(CT_24hr_pcas$x[,1])))

#get 2 and 3 s.ds
#Human
human_1hr_2_sd <- human_1hr_mean + human_1hr_sd + human_1hr_sd
human_1hr_3_sd <- human_1hr_mean + human_1hr_sd + human_1hr_sd + human_1hr_sd
human_24hr_2_sd <- human_24hr_mean + human_24hr_sd + human_24hr_sd
human_24hr_3_sd <- human_24hr_mean + human_24hr_sd + human_24hr_sd + human_24hr_sd
#CT
CT_1hr_2_sd <- CT_1hr_mean + CT_1hr_sd + CT_1hr_sd
CT_1hr_3_sd <- CT_1hr_mean + CT_1hr_sd + CT_1hr_sd + CT_1hr_sd
CT_24hr_2_sd <- CT_24hr_mean + CT_24hr_sd + CT_24hr_sd
CT_24hr_3_sd <- CT_24hr_mean + CT_24hr_sd + CT_24hr_sd + CT_24hr_sd


#--
#Are there any outliers on PC1
#--

#Human
row.names(human_1hr_pcas$x)[which(abs(as.numeric(human_1hr_pcas$x[,1]))>human_1hr_3_sd)]
row.names(human_1hr_pcas$x)[which(abs(as.numeric(human_1hr_pcas$x[,1]))>human_1hr_2_sd)]
row.names(human_24hr_pcas$x)[which(abs(as.numeric(human_24hr_pcas$x[,1]))>human_24hr_3_sd)]
row.names(human_24hr_pcas$x)[which(abs(as.numeric(human_24hr_pcas$x[,1]))>human_24hr_2_sd)]
#CT
row.names(CT_1hr_pcas$x)[which(abs(as.numeric(CT_1hr_pcas$x[,1]))>CT_1hr_3_sd)]
row.names(CT_1hr_pcas$x)[which(abs(as.numeric(CT_1hr_pcas$x[,1]))>CT_1hr_2_sd)]
row.names(CT_24hr_pcas$x)[which(abs(as.numeric(CT_24hr_pcas$x[,1]))>CT_24hr_3_sd)]
row.names(CT_24hr_pcas$x)[which(abs(as.numeric(CT_24hr_pcas$x[,1]))>CT_24hr_2_sd)]

#--
#Are there any outliers on PC2
#--

#Human
row.names(human_1hr_pcas$x)[which(abs(as.numeric(human_1hr_pcas$x[,2]))>human_1hr_3_sd)]
row.names(human_1hr_pcas$x)[which(abs(as.numeric(human_1hr_pcas$x[,2]))>human_1hr_2_sd)]
row.names(human_24hr_pcas$x)[which(abs(as.numeric(human_24hr_pcas$x[,2]))>human_24hr_3_sd)]
row.names(human_24hr_pcas$x)[which(abs(as.numeric(human_24hr_pcas$x[,2]))>human_24hr_2_sd)]
#CT
row.names(CT_1hr_pcas$x)[which(abs(as.numeric(CT_1hr_pcas$x[,2]))>CT_1hr_3_sd)]
row.names(CT_1hr_pcas$x)[which(abs(as.numeric(CT_1hr_pcas$x[,2]))>CT_1hr_2_sd)]
row.names(CT_24hr_pcas$x)[which(abs(as.numeric(CT_24hr_pcas$x[,2]))>CT_24hr_3_sd)]
row.names(CT_24hr_pcas$x)[which(abs(as.numeric(CT_24hr_pcas$x[,2]))>CT_24hr_2_sd)]



#--
#2) Now looking at hierarchical clustering
#--

#Need to understand the data shape
#Norm data will ne neg bionomial. Taking the log will give normal, from here: https://www.biostars.org/p/298155/
hist(norm_counts_human_24hr_v1)
hist(log2(norm_counts_human_24hr_v1))
hist(log2(norm_counts_human_24hr_v1+0.1))

#calc z-scores
zscore <- function(x){
  z<- (x - mean(x)) / sd(x)
  return(z)
}

#human
z.output_human_1hr <- zscore(log2(human_rle_plots_norm_1hr+0.1))
z.output_human_24hr <- zscore(log2(human_rle_plots_norm_24hr+0.1))
#CTn
z.output_CT_1hr <- zscore(log2(CT_rle_plots_norm_1hr+0.1))
z.output_CT_24hr <- zscore(log2(CT_rle_plots_norm_24hr+0.1))


width=6;height=6

#Human
Cairo(file="dendogram_human_1hr_samples.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
par(mfrow = c(1,1), mar = c(18,4,4,2))
dend <- as.dendrogram(hclust(dist(t(z.output_human_1hr),  method = "euclidean"), method = 'ward.D2' ))
dend %>% set("labels_colors", unique(plot_cols)[c(4,6,6,2,2,4,2,4,5,5,5,1,1,1,3,3,3)]) %>% plot(main="Human - 1 hr - removed 1 outlier")
dev.off()

Cairo(file="dendogram_human_24hr_samples.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
par(mfrow = c(1,1), mar = c(18,4,4,2))
dend <- as.dendrogram(hclust(dist(t(z.output_human_24hr),  method = "euclidean"), method = 'ward.D2' ))
dend %>% set("labels_colors", unique(plot_cols)[c(2,3,2,3,2,3,3,1,1,1,3,3,5,5,5,6,6)]) %>% plot(main="Human - 24 hr - removed 1 outlier")
dev.off()

#CT
Cairo(file="dendogram_CT_1hr_all_samples.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
par(mfrow = c(1,1), mar = c(18,4,4,2))
dend <- as.dendrogram(hclust(dist(t(z.output_CT_1hr),  method = "euclidean"), method = 'ward.D2' ))
dend %>% set("labels_colors", unique(plot_cols)[c(3,3,4,3,6,6,5,5,5,4,6,1,2,1,2,4,1,2)]) %>% plot(main="CT - 1 hr - all samples")
dev.off()

Cairo(file="dendogram_CT_24hr_samples.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
par(mfrow = c(1,1), mar = c(18,4,4,2))
dend <- as.dendrogram(hclust(dist(t(z.output_CT_24hr),  method = "euclidean"), method = 'ward.D2' ))
dend %>% set("labels_colors", unique(plot_cols)[c(2,6,6,6,4,4,4,1,1,1,5,5,5,3,3,3)]) %>% plot(main="CT - 24 hr - 2 outliers removed")
dev.off()




#--
#3) Looking at PCA plots removing the bottom 10% of variance
#--

#Set up plot settings
width=7;height=5


#Prepare the metadata
pca_metadata_1hr <-  data.frame("names" = counts_featureCounts_human_v2$samples$files[1:18],
                                "col" = unique(plot_cols)[c(1,1,1,2,2,2,5,5,5,6,6,6,3,3,3,4,4,4)])

pca_metadata_1hr$names <- substr(as.vector(pca_metadata_1hr$names),1,nchar(as.vector(pca_metadata_1hr$names))-24)


#Human
Cairo(file="pcaTools_biplot_human_1hr_all_samples.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
p <- pca(mat=log2(human_rle_plots_norm_1hr +1), removeVar = 0.1, scale = T) #by removing 0.1%, seem to get more variation in PCA 1 and 2
biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, 
       colkey = unique(plot_cols)[c(1,1,1,2,2,2,5,5,5,6,6,3,3,3,4,4,4)], pointSize = 5, title="Human 1hr\n(scaled + 10% var removed)\n1 outlier removed")
dev.off()

Cairo(file="pcaTools_biplot_human_24hr_all_samples.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
p <- pca(mat=log2(human_rle_plots_norm_24hr + 1), removeVar = 0.1, scale = T) #by removing 0.1%, seem to get more variation in PCA 1 and 2
biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, 
       colkey = unique(plot_cols)[c(1,1,1,2,2,2,5,5,5,6,6,3,3,3,4,4,4)], pointSize = 5, title="Human 24hr\n(scaled + 10% var removed)\n1 outlier removed")
dev.off()


#CT
Cairo(file="pcaTools_biplot_CT_1hr_all_samples.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
p <- pca(mat=log2(CT_rle_plots_norm_1hr + 1), removeVar = 0.1, scale = T)
biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, colkey = plot_cols, pointSize = 5, title="CT 1hr\n(scaled + 10% var removed)\n0 outliers removed")
dev.off()

Cairo(file="pcaTools_biplot_CT_24hr_all_samples.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
p <- pca(mat=log2(CT_rle_plots_norm_24hr + 1), removeVar = 0.1, scale = T) #by removing 0.1%, seem to get more variation in PCA 1 and 2
biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, 
       colkey = unique(plot_cols)[c(1,1,1,2,5,5,5,6,6,6,3,3,3,4,4,4)], pointSize = 5, title="CT 24hr\n(scaled + 10% var removed)\n1 outlier removed")
dev.off()




#Looking at genes and the amount of variation they drive
#all CT 1 hr
p <- pca(mat=log2(CT_rle_plots_norm_1hr + 1), removeVar = 0.1, scale = T)
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.01, absolute = T, labSize = 4.0, title = 'CT - 1hr', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 1% variables', shape = 21, shapeSizeRange = c(4,10), col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE)


#All CT 24 hours
p <- pca(mat=log2(CT_rle_plots_norm_24hr + 1), removeVar = 0.1, scale = T) #by removing 0.1%, seem to get more variation in PCA 1 and 2
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.005, absolute = T, labSize = 3.0, title = 'CT - 24hrs', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 0.5% variables', shape = 21, shapeSizeRange = c(2,8) ,col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE)


#all human 1hr
p <- pca(mat=log2(human_rle_plots_norm_1hr +1), removeVar = 0.1, scale = T)
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.0005, absolute = T, labSize = 3.0, title = 'CT - 24hrs', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 0.05% variables', shape = 21, shapeSizeRange = c(2,8) ,col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE)


#all human 24hrs
p <- pca(mat=log2(human_rle_plots_norm_24hr +1), removeVar = 0.1, scale = T)
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.0005, absolute = T, labSize = 3.0, title = 'CT - 24hrs', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 0.05% variables', shape = 21, shapeSizeRange = c(2,8) ,col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE)








#**************************************
#
# Annotate the human genes ----
#
#**************************************

#--
#Connect to the annotation db
#--
edb <- EnsDb.Hsapiens.v86 #connect
columns(edb) #The available columns

#Get the list of genes
ens_cols <- c("GENEID","GENENAME","GENEBIOTYPE")

#annotate
human_annotation_1hr_v1 <- select(edb,keys=as.vector(row.names(human_rle_plots_norm_1hr)),columns=ens_cols,keytype="GENEID")
human_annotation_24hr_v1 <- select(edb,keys=as.vector(row.names(human_rle_plots_norm_24hr)),columns=ens_cols,keytype="GENEID")

#Make a df and add in the gene annotation info
human_annotation_1hr_v2 <- as.data.frame(human_rle_plots_norm_1hr)
human_annotation_24hr_v2 <- as.data.frame(human_rle_plots_norm_24hr)

#Combine
human_annotation_1hr_v3 <- cbind(human_annotation_1hr_v2, human_annotation_1hr_v1)
human_annotation_24hr_v3 <- cbind(human_annotation_24hr_v2, human_annotation_24hr_v1)





#**************************************
#
# Determine the genetic variation between PC1 and PC2 (HUMAN)----
#
#**************************************


#---
#Want to deterine what genes are driving the variation between depletion methods
#---

#need to get the actual gene names as the row names, not the ENSG names
human_norm_1hr_annotated_df <- human_annotation_1hr_v3
human_norm_24hr_annotated_df <- human_annotation_24hr_v3

#change the row.names (allows for duplicate names, adds a .1 or .2 after them!)
rownames(human_norm_1hr_annotated_df) <- make.names(human_norm_1hr_annotated_df$GENENAME,unique = T)
rownames(human_norm_24hr_annotated_df) <- make.names(human_norm_24hr_annotated_df$GENENAME,unique = T)


#human 1hr
colnames(human_norm_1hr_annotated_df)
plot_loadings_human_1hr_moi_01 <- human_norm_1hr_annotated_df[,c(1:6)]
plot_loadings_human_1hr_moi_1 <- human_norm_1hr_annotated_df[,c(7:12)]
plot_loadings_human_1hr_moi_10 <- human_norm_1hr_annotated_df[,c(13:17)]
#human 24hr
colnames(human_norm_24hr_annotated_df)
plot_loadings_human_24hr_moi_01 <- human_norm_24hr_annotated_df[,c(1:6)]
plot_loadings_human_24hr_moi_1 <- human_norm_24hr_annotated_df[,c(7:12)]
plot_loadings_human_24hr_moi_10 <- human_norm_24hr_annotated_df[,c(13:17)]



#--Human - 1hr moi = 0.1
p <- pca(mat=log2(plot_loadings_human_1hr_moi_01 +1), removeVar = 0.1, scale = T) #convert into the pca format

#bi plot
Cairo(file="Human_1hr_moi_01_PCA.png", type="png", units="in", width=6, height=5, dpi=300, bg = "transparent")
biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, 
       colkey = unique(plot_cols)[c(1,1,1,2,2,2)], pointSize = 4, title="Human - 1hr - MOI 0.1\n(scaled + 10% var removed)")
dev.off()
#loadings plot
Cairo(file="Human_1hr_moi_01_loadings.png", type="png", units="in", width=5, height=5, dpi=300, bg = "transparent")
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.00005, absolute = T, labSize = 3.0, title = 'Human - 1hr - moi 0.1', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 0.005% variables', shape = 21, shapeSizeRange = c(2,8) ,col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE)
dev.off()


#--Human - 1hr moi = 1
p <- pca(mat=log2(plot_loadings_human_1hr_moi_1 +1), removeVar = 0.1, scale = T) #convert into the pca format
#bi plot
Cairo(file="Human_1hr_moi_1_PCA.png", type="png", units="in", width=6, height=5, dpi=300, bg = "transparent")
biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, 
       colkey = unique(plot_cols)[c(3,3,3,4,4,4)], pointSize = 5, title="Human - 1hr - MOI 1\n(scaled + 10% var removed)")
dev.off()
#loadings plot
Cairo(file="Human_1hr_moi_1_loadings.png", type="png", units="in", width=5, height=5, dpi=300, bg = "transparent")
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.00005, absolute = T, labSize = 3.0, title = 'Human - 1hr - moi 1', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 0.005% variables', shape = 21, shapeSizeRange = c(2,8) ,col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE)
dev.off()


#--Human - 1hr moi = 10
p <- pca(mat=log2(plot_loadings_human_1hr_moi_10 +1), removeVar = 0.1, scale = T) #convert into the pca format
#bi plot
Cairo(file="Human_1hr_moi_10_PCA.png", type="png", units="in", width=6, height=5, dpi=300, bg = "transparent")
biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, 
       colkey = unique(plot_cols)[c(5,5,5,6,6)], pointSize = 5, title="Human - 1hr - MOI 10\n(scaled + 10% var removed)")
dev.off()
#loadings plot
Cairo(file="Human_1hr_moi_10_loadings.png", type="png", units="in", width=5, height=5, dpi=300, bg = "transparent")
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.00005, absolute = T, labSize = 3.0, title = 'Human - 1hr - moi 10', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 0.005% variables', shape = 21, shapeSizeRange = c(2,8) ,col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE, ylim = c(0,0.020))
dev.off()


#Genes
human_1hr_moi_01_genes <- c("NT5DC3", "AC010642.1", "SRGN", "CDC42EP1", "ATAD2") #moi 0.1
human_1hr_moi_1_genes <- c("TRIP12", "TNKS2", "PA2G4", "TK1", "HSF2", "SFI1") #moi 1
human_1hr_moi_10_genes <- c("ARIH1", "TPM3", "ISCA1", "H3F3B", "TK1", "PRMT1", "VMA21", "ALKBH8", "GNB1") #moi 10

#Look at the biotypes by subsetting into the original df
human_annotation_1hr_v3[human_annotation_1hr_v3$GENENAME %in% human_1hr_moi_01_genes,18:20]
human_annotation_1hr_v3[human_annotation_1hr_v3$GENENAME %in% human_1hr_moi_1_genes,18:20]
human_annotation_1hr_v3[human_annotation_1hr_v3$GENENAME %in% human_1hr_moi_10_genes,18:20]




#--Human - 24hr moi = 0.1
p <- pca(mat=log2(plot_loadings_human_24hr_moi_01 +1), removeVar = 0.1, scale = T) #convert into the pca format
#bi plot
Cairo(file="Human_24hr_moi_01_PCA.png", type="png", units="in", width=6, height=5, dpi=300, bg = "transparent")
biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, 
       colkey = unique(plot_cols)[c(1,1,1,2,2,2)], pointSize = 4, title="Human - 24hr - MOI 0.1\n(scaled + 10% var removed)")
dev.off()
#loadings plot
Cairo(file="Human_24hr_moi_01_loadings.png", type="png", units="in", width=5, height=5, dpi=300, bg = "transparent")
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.00005, absolute = T, labSize = 3.0, title = 'Human - 24hr - moi 0.1', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 0.005% variables', shape = 21, shapeSizeRange = c(2,8) ,col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE)
dev.off()


#--Human - 24hr moi = 1
p <- pca(mat=log2(plot_loadings_human_24hr_moi_1 +1), removeVar = 0.1, scale = T) #convert into the pca format
#bi plot
Cairo(file="Human_24hr_moi_1_PCA.png", type="png", units="in", width=6, height=5, dpi=300, bg = "transparent")
biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, 
       colkey = unique(plot_cols)[c(3,3,3,4,4,4)], pointSize = 5, title="Human - 24hr - MOI 1\n(scaled + 10% var removed)")
dev.off()
#loadings plot
Cairo(file="Human_24hr_moi_1_loadings.png", type="png", units="in", width=5, height=5, dpi=300, bg = "transparent")
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.00005, absolute = T, labSize = 3.0, title = 'Human - 24hr - moi 1', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 0.005% variables', shape = 21, shapeSizeRange = c(2,8) ,col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE)
dev.off()


#--Human - 24hr moi = 10
p <- pca(mat=log2(plot_loadings_human_24hr_moi_10 +1), removeVar = 0.1, scale = T) #convert into the pca format
#bi plot
Cairo(file="Human_24hr_moi_10_PCA.png", type="png", units="in", width=6, height=5, dpi=300, bg = "transparent")
biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, 
       colkey = unique(plot_cols)[c(5,5,5,6,6)], pointSize = 5, title="Human - 24hr - MOI 10\n(scaled + 10% var removed)")
dev.off()
#loadings plot
Cairo(file="Human_24hr_moi_10_loadings.png", type="png", units="in", width=5, height=5, dpi=300, bg = "transparent")
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.00005, absolute = T, labSize = 3.0, title = 'Human - 24hr - moi 10', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 0.005% variables', shape = 21, shapeSizeRange = c(2,8) ,col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE, ylim = c(0,0.030))
dev.off()


#Genes
human_24hr_moi_01_genes <- c("FAT4", "RPS4X", "CBWD2", "STRN3") #moi 0.1
human_24hr_moi_1_genes <- c("FARP1", "TPM3", "RPL24", "RACK1", "DBF4B", "CLK1") #moi 1
human_24hr_moi_10_genes <- c("RTTN", "CDK5", "EMC7", "SPIN1", "PHF1") #moi 10

#Look at the biotypes by subsetting into the original df
human_annotation_24hr_v3[human_annotation_24hr_v3$GENENAME %in% human_24hr_moi_01_genes,18:20]
human_annotation_24hr_v3[human_annotation_24hr_v3$GENENAME %in% human_24hr_moi_1_genes,18:20]
human_annotation_24hr_v3[human_annotation_24hr_v3$GENENAME %in% human_24hr_moi_10_genes,18:20]







#**************************************
#
# Determine the genetic variation between PC1 and PC2 (CT)----
#
#**************************************


#---
#Want to deterine what genes are driving the variation between depletion methods
#---

#Getting the annotation from below 
CT_features3
CT_rle_plots_norm_1hr
CT_rle_plots_norm_24hr

#merge
CT_norm_1hr_annotated_df <- as.data.frame(CT_rle_plots_norm_1hr)
CT_norm_1hr_annotated_df$Gene <- row.names(CT_rle_plots_norm_1hr)
CT_norm_1hr_annotated_df_v2 <- merge(x=CT_norm_1hr_annotated_df, y=CT_features3, by.x="Gene", by.y="ID", all.x=T)
CT_norm_1hr_annotated <- CT_rle_plots_norm_1hr
row.names(CT_norm_1hr_annotated) <- CT_norm_1hr_annotated_df_v2$Name

CT_norm_24hr_annotated_df <- as.data.frame(CT_rle_plots_norm_24hr)
CT_norm_24hr_annotated_df$Gene <- row.names(CT_rle_plots_norm_24hr)
CT_norm_24hr_annotated_df_v2 <- merge(x=CT_norm_24hr_annotated_df, y=CT_features3, by.x="Gene", by.y="ID", all.x=T)
CT_norm_24hr_annotated <- CT_rle_plots_norm_24hr
row.names(CT_norm_24hr_annotated) <- CT_norm_24hr_annotated_df_v2$Name


#CT 1hr
colnames(CT_norm_1hr_annotated)
plot_loadings_CT_1hr_moi_01 <- CT_norm_1hr_annotated[,c(1:6)]
plot_loadings_CT_1hr_moi_1 <- CT_norm_1hr_annotated[,c(7:12)]
plot_loadings_CT_1hr_moi_10 <- CT_norm_1hr_annotated[,c(13:18)]
#CT 24hr
colnames(CT_norm_24hr_annotated)
plot_loadings_CT_24hr_moi_01 <- CT_norm_24hr_annotated[,c(1:4)]
plot_loadings_CT_24hr_moi_1 <- CT_norm_24hr_annotated[,c(5:10)]
plot_loadings_CT_24hr_moi_10 <- CT_norm_24hr_annotated[,c(11:16)]



#--CT - 1hr moi = 0.1
p <- pca(mat=log2(plot_loadings_CT_1hr_moi_01 +1), removeVar = 0.1, scale = T) #convert into the pca format
#bi plot
Cairo(file="CT_1hr_moi_01_PCA.png", type="png", units="in", width=6, height=5, dpi=300, bg = "transparent")
biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, 
       colkey = unique(plot_cols)[c(1,1,1,2,2,2)], pointSize = 4, title="CT - 1hr - MOI 0.1\n(scaled + 10% var removed)")
dev.off()
#loadings plot
Cairo(file="CT_1hr_moi_01_loadings.png", type="png", units="in", width=5, height=5, dpi=300, bg = "transparent")
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.00005, absolute = T, labSize = 3.0, title = 'CT - 1hr - moi 0.1', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 0.005% variables', shape = 21, shapeSizeRange = c(2,8) ,col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE)
dev.off()


#--CT - 1hr moi = 1
p <- pca(mat=log2(plot_loadings_CT_1hr_moi_1 +1), removeVar = 0.1, scale = T) #convert into the pca format
#bi plot
Cairo(file="CT_1hr_moi_1_PCA.png", type="png", units="in", width=6, height=5, dpi=300, bg = "transparent")
biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, 
       colkey = unique(plot_cols)[c(3,3,3,4,4,4)], pointSize = 5, title="CT - 1hr - MOI 1\n(scaled + 10% var removed)")
dev.off()
#loadings plot
Cairo(file="CT_1hr_moi_1_loadings.png", type="png", units="in", width=5, height=5, dpi=300, bg = "transparent")
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.00005, absolute = T, labSize = 3.0, title = 'CT - 1hr - moi 1', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 0.005% variables', shape = 21, shapeSizeRange = c(2,8) ,col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE)
dev.off()


#--CT - 1hr moi = 10
p <- pca(mat=log2(plot_loadings_CT_1hr_moi_10 +1), removeVar = 0.1, scale = T) #convert into the pca format
#bi plot
Cairo(file="CT_1hr_moi_10_PCA.png", type="png", units="in", width=6, height=5, dpi=300, bg = "transparent")
biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, 
       colkey = unique(plot_cols)[c(5,5,5,6,6,6)], pointSize = 5, title="CT - 1hr - MOI 10\n(scaled + 10% var removed)")
dev.off()
#loadings plot
Cairo(file="CT_1hr_moi_10_loadings.png", type="png", units="in", width=5, height=5, dpi=300, bg = "transparent")
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.00005, absolute = T, labSize = 3.0, title = 'CT - 1hr - moi 10', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 0.005% variables', shape = 21, shapeSizeRange = c(2,8) ,col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE )
dev.off()


#Genes
CT_1hr_moi_01_genes <- c("ct027v1_405", "dcd", "ct027v1_407", "ct027v1_946") #moi 0.1
CT_1hr_moi_1_genes <- c("ct027v1_415", "pyrH", "ct027v1_345", "proS") #moi 1
CT_1hr_moi_10_genes <- c("dnaN", "ct027v1_171", "ct027v1_771", "ct027v1_614") #moi 10

#Look at the biotypes by subsetting into the original df
CT_features3[CT_features3$Name %in% CT_1hr_moi_01_genes,]
CT_features3[CT_features3$Name %in% CT_1hr_moi_1_genes,]
CT_features3[CT_features3$Name %in% CT_1hr_moi_10_genes,]






#--CT - 24hr moi = 0.1
p <- pca(mat=log2(plot_loadings_CT_24hr_moi_01 +1), removeVar = 0.1, scale = T) #convert into the pca format
#bi plot
Cairo(file="CT_24hr_moi_01_PCA.png", type="png", units="in", width=6, height=5, dpi=300, bg = "transparent")
biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, 
       colkey = unique(plot_cols)[c(1,1,1,2,2,2)], pointSize = 4, title="CT - 24hr - MOI 0.1\n(scaled + 10% var removed)")
dev.off()
#loadings plot
Cairo(file="CT_24hr_moi_01_loadings.png", type="png", units="in", width=5, height=5, dpi=300, bg = "transparent")
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.00005, absolute = T, labSize = 3.0, title = 'CT - 24hr - moi 0.1', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 0.005% variables', shape = 21, shapeSizeRange = c(2,8) ,col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE)
dev.off()


#--CT - 24hr moi = 1
p <- pca(mat=log2(plot_loadings_CT_24hr_moi_1 +1), removeVar = 0.1, scale = T) #convert into the pca format
#bi plot
Cairo(file="CT_24hr_moi_1_PCA.png", type="png", units="in", width=6, height=5, dpi=300, bg = "transparent")
biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, 
       colkey = unique(plot_cols)[c(3,3,3,4,4,4)], pointSize = 5, title="CT - 24hr - MOI 1\n(scaled + 10% var removed)")
dev.off()
#loadings plot
Cairo(file="CT_24hr_moi_1_loadings.png", type="png", units="in", width=5, height=5, dpi=300, bg = "transparent")
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.00005, absolute = T, labSize = 3.0, title = 'CT - 24hr - moi 1', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 0.005% variables', shape = 21, shapeSizeRange = c(2,8) ,col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE)
dev.off()


#--CT - 24hr moi = 10
p <- pca(mat=log2(plot_loadings_CT_24hr_moi_10 +1), removeVar = 0.1, scale = T) #convert into the pca format
#bi plot
Cairo(file="CT_24hr_moi_10_PCA.png", type="png", units="in", width=6, height=5, dpi=300, bg = "transparent")
biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, 
       colkey = unique(plot_cols)[c(5,5,5,6,6,6)], pointSize = 5, title="CT - 24hr - MOI 10\n(scaled + 10% var removed)")
dev.off()
#loadings plot
Cairo(file="CT_24hr_moi_10_loadings.png", type="png", units="in", width=5, height=5, dpi=300, bg = "transparent")
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.00005, absolute = T, labSize = 3.0, title = 'CT - 24hr - moi 10', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 0.005% variables', shape = 21, shapeSizeRange = c(2,8) ,col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE )
dev.off()


#Genes
CT_24hr_moi_01_genes <- c("ct027v1_993", "rplM", "aroE", "ct027v1_365") #moi 0.1
CT_24hr_moi_1_genes <- c("ct027v1_820", "ct027v1_642", "ct027v1_649", "hemB") #moi 1
CT_24hr_moi_10_genes <- c("fumC", "ct027v1_556", "ct027v1_308", "rplU") #moi 10

#Look at the biotypes by subsetting into the original df
CT_features3[CT_features3$Name %in% CT_24hr_moi_01_genes,]
CT_features3[CT_features3$Name %in% CT_24hr_moi_1_genes,]
CT_features3[CT_features3$Name %in% CT_24hr_moi_10_genes,]
#



#Now looking at the top 5% of variance from 1hr - CT
#--CT - 1hr moi = 0.1
p <- pca(mat=log2(plot_loadings_CT_1hr_moi_01 +1), removeVar = 0.1, scale = T) #convert into the pca format
#bi plot

biplot(p, lab = FALSE, legendPosition = "right", legendLabSize = 6, legendIconSize = 3.0, 
       colkey = unique(plot_cols)[c(1,1,1,2,2,2)], pointSize = 4, title="CT - 1hr - MOI 0.1\n(scaled + 10% var removed)")

#loadings plot
plotloadings(p, components = getComponents(p,seq_len(2)), rangeRetain = 0.05, absolute = T, labSize = 3.0, title = 'CT - 1hr - moi 0.1', subtitle = 'Loadings of (PC1 and PC2)',
             caption = 'Top 5% variability', shape = 21, shapeSizeRange = c(2,8) ,col = c('limegreen', 'black', 'red3'), drawConnectors = TRUE)








#**************************************
#
# looking at the top 5% variance of genes, then making venn diagrams of the overlap ----
#
#**************************************


#1hr
plotload_CT_1hr_0.1_v1 <- pca(mat=log2(plot_loadings_CT_1hr_moi_01 +1), removeVar = 0.1, scale = T) #convert into the pca format
plotload_CT_1hr_0.1_v2 <- plotloadings(plotload_CT_1hr_0.1_v1, components = getComponents(plotload_CT_1hr_0.1_v1,seq_len(2)), rangeRetain = 0.05, absolute = T)

plotload_CT_1hr_1_v1 <- pca(mat=log2(plot_loadings_CT_1hr_moi_1 +1), removeVar = 0.1, scale = T) #convert into the pca format
plotload_CT_1hr_1_v2 <- plotloadings(plotload_CT_1hr_1_v1, components = getComponents(plotload_CT_1hr_1_v1,seq_len(2)), rangeRetain = 0.05, absolute = T)

plotload_CT_1hr_10_v1 <- pca(mat=log2(plot_loadings_CT_1hr_moi_10 +1), removeVar = 0.1, scale = T) #convert into the pca format
plotload_CT_1hr_10_v2 <- plotloadings(plotload_CT_1hr_10_v1, components = getComponents(plotload_CT_1hr_10_v1,seq_len(2)), rangeRetain = 0.05, absolute = T)

Reduce(intersect, list(plotload_CT_1hr_0.1_v2$data$var,plotload_CT_1hr_1_v2$data$var,plotload_CT_1hr_10_v2$data$var))
length(Reduce(intersect, list(plotload_CT_1hr_0.1_v2$data$var,plotload_CT_1hr_1_v2$data$var,plotload_CT_1hr_10_v2$data$var)))


area1 = length(plotload_CT_1hr_0.1_v2$data$var)
area2 = length(plotload_CT_1hr_1_v2$data$var)
area3 = length(plotload_CT_1hr_10_v2$data$var)
n12 = length(Reduce(intersect, list(plotload_CT_1hr_0.1_v2$data$var,plotload_CT_1hr_1_v2$data$var)))
n13 = length(Reduce(intersect, list(plotload_CT_1hr_0.1_v2$data$var,plotload_CT_1hr_10_v2$data$var)))
n23 = length(Reduce(intersect, list(plotload_CT_1hr_1_v2$data$var,plotload_CT_1hr_10_v2$data$var)))
n123 = length(Reduce(intersect, list(plotload_CT_1hr_0.1_v2$data$var,plotload_CT_1hr_1_v2$data$var,plotload_CT_1hr_10_v2$data$var)))


height=5;width=5
Cairo(file="venn.diagram_1hr_CT_overlaps.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
draw.triple.venn(area1, area2, area3, n12, n13, n23, n123, 
               category = c("0.1","1","10"), main="gg", alpha = 0.4, fill = unique(plot_cols)[c(2,4,6)], col = "white",
               cat.col = "black")
dev.off()



#24hr
plotload_CT_24hr_0.1_v1 <- pca(mat=log2(plot_loadings_CT_24hr_moi_01 +1), removeVar = 0.1, scale = T) #convert into the pca format
plotload_CT_24hr_0.1_v2 <- plotloadings(plotload_CT_24hr_0.1_v1, components = getComponents(plotload_CT_24hr_0.1_v1,seq_len(2)), rangeRetain = 0.05, absolute = T)

plotload_CT_24hr_1_v1 <- pca(mat=log2(plot_loadings_CT_24hr_moi_1 +1), removeVar = 0.1, scale = T) #convert into the pca format
plotload_CT_24hr_1_v2 <- plotloadings(plotload_CT_24hr_1_v1, components = getComponents(plotload_CT_24hr_1_v1,seq_len(2)), rangeRetain = 0.05, absolute = T)

plotload_CT_24hr_10_v1 <- pca(mat=log2(plot_loadings_CT_24hr_moi_10 +1), removeVar = 0.1, scale = T) #convert into the pca format
plotload_CT_24hr_10_v2 <- plotloadings(plotload_CT_24hr_10_v1, components = getComponents(plotload_CT_24hr_10_v1,seq_len(2)), rangeRetain = 0.05, absolute = T)

Reduce(intersect, list(plotload_CT_24hr_0.1_v2$data$var,plotload_CT_24hr_1_v2$data$var,plotload_CT_24hr_10_v2$data$var))
length(Reduce(intersect, list(plotload_CT_24hr_0.1_v2$data$var,plotload_CT_24hr_1_v2$data$var,plotload_CT_24hr_10_v2$data$var)))


area1 = length(plotload_CT_24hr_0.1_v2$data$var)
area2 = length(plotload_CT_24hr_1_v2$data$var)
area3 = length(plotload_CT_24hr_10_v2$data$var)
n12 = length(Reduce(intersect, list(plotload_CT_24hr_0.1_v2$data$var,plotload_CT_24hr_1_v2$data$var)))
n13 = length(Reduce(intersect, list(plotload_CT_24hr_0.1_v2$data$var,plotload_CT_24hr_10_v2$data$var)))
n23 = length(Reduce(intersect, list(plotload_CT_24hr_1_v2$data$var,plotload_CT_24hr_10_v2$data$var)))
n123 = length(Reduce(intersect, list(plotload_CT_24hr_0.1_v2$data$var,plotload_CT_24hr_1_v2$data$var,plotload_CT_24hr_10_v2$data$var)))



height=5;width=5

Cairo(file="venn.diagram_24hr_CT_overlaps.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
draw.triple.venn(area1, area2, area3, n12, n13, n23, n123, 
                 category = c("0.1","1","10"), main="gg", alpha = 0.4, fill = unique(plot_cols)[c(2,4,6)], col = "white",
                 cat.col = "black")
dev.off()










#**************************************
#
# Look at the % of protein coding genes (comparing rRNA methods) ----
#
#**************************************

#Need to annotate the normalised genes to get the biotypes (just human as no biotypes for CT genes unfortunately)


#--
#Get the average of the replicates and make new columns with their averages
#--

#Human
#1hr
colnames(human_annotation_1hr_v3)
human_annotation_1hr_v3$rRNAdep_01_1hr_Avg <- rowSums(human_annotation_1hr_v3[,1:3]) / 3
human_annotation_1hr_v3$rRNAPolyA_01_1hr_Avg <- rowSums(human_annotation_1hr_v3[,4:6]) /3
human_annotation_1hr_v3$rRNAdep_1_1hr_Avg <- rowSums(human_annotation_1hr_v3[,7:9]) / 3
human_annotation_1hr_v3$rRNAPolyA_1_1hr_Avg <- rowSums(human_annotation_1hr_v3[,10:12]) / 3
human_annotation_1hr_v3$rRNAdep_10_1hr_Avg <- rowSums(human_annotation_1hr_v3[,13:15]) / 3
human_annotation_1hr_v3$rRNAPolyA_10_1hr_Avg <- rowSums(human_annotation_1hr_v3[,16:17]) / 2

#24hr
colnames(human_annotation_24hr_v3)
human_annotation_24hr_v3$rRNAdep_01_24hr_Avg <- rowSums(human_annotation_24hr_v3[,1:3]) / 3
human_annotation_24hr_v3$rRNAPolyA_01_24hr_Avg <- rowSums(human_annotation_24hr_v3[,4:6]) / 3
human_annotation_24hr_v3$rRNAdep_1_24hr_Avg <- rowSums(human_annotation_24hr_v3[,7:9]) / 3
human_annotation_24hr_v3$rRNAPolyA_1_24hr_Avg <- rowSums(human_annotation_24hr_v3[,10:12]) / 3
human_annotation_24hr_v3$rRNAdep_10_24hr_Avg <- rowSums(human_annotation_24hr_v3[,13:15]) / 3
human_annotation_24hr_v3$rRNAPolyA_10_24hr_Avg <- rowSums(human_annotation_24hr_v3[,16:17]) / 2



#--
#Using dplyr to filter just protein coding genes
#--
just_protein_coding_1hr <- human_annotation_1hr_v3 %>% filter(GENEBIOTYPE =="protein_coding")
just_protein_coding_24hr <- human_annotation_24hr_v3 %>% filter(GENEBIOTYPE =="protein_coding")


#Save variables (%)
#1hr
human_pie_pc_1hr_rRNAdep_01 <- sum(just_protein_coding_1hr$rRNAdep_01_1hr_Avg) / sum(human_annotation_1hr_v3$rRNAdep_01_1hr_Avg) * 100
human_pie_pc_1hr_rRNAdep_1 <- sum(just_protein_coding_1hr$rRNAdep_1_1hr_Avg) / sum(human_annotation_1hr_v3$rRNAdep_1_1hr_Avg) * 100
human_pie_pc_1hr_rRNAdep_10 <- sum(just_protein_coding_1hr$rRNAdep_10_1hr_Avg) / sum(human_annotation_1hr_v3$rRNAdep_10_1hr_Avg) * 100
human_pie_pc_1hr_rRNApolyA_01 <- sum(just_protein_coding_1hr$rRNAPolyA_01_1hr_Avg) / sum(human_annotation_1hr_v3$rRNAPolyA_01_1hr_Avg) * 100
human_pie_pc_1hr_rRNApolyA_1 <- sum(just_protein_coding_1hr$rRNAPolyA_1_1hr_Avg) / sum(human_annotation_1hr_v3$rRNAPolyA_1_1hr_Avg) * 100
human_pie_pc_1hr_rRNApolyA_10 <- sum(just_protein_coding_1hr$rRNAPolyA_10_1hr_Avg) / sum(human_annotation_1hr_v3$rRNAPolyA_10_1hr_Avg) * 100

#24hr
human_pie_pc_24hr_rRNAdep_01 <- sum(just_protein_coding_24hr$rRNAdep_01_24hr_Avg) / sum(human_annotation_24hr_v3$rRNAdep_01_24hr_Avg) * 100
human_pie_pc_24hr_rRNAdep_1 <- sum(just_protein_coding_24hr$rRNAdep_1_24hr_Avg) / sum(human_annotation_24hr_v3$rRNAdep_1_24hr_Avg) * 100
human_pie_pc_24hr_rRNAdep_10 <- sum(just_protein_coding_24hr$rRNAdep_10_24hr_Avg) / sum(human_annotation_24hr_v3$rRNAdep_10_24hr_Avg) * 100
human_pie_pc_24hr_rRNApolyA_01 <- sum(just_protein_coding_24hr$rRNAPolyA_01_24hr_Avg) / sum(human_annotation_24hr_v3$rRNAPolyA_01_24hr_Avg) * 100
human_pie_pc_24hr_rRNApolyA_1 <- sum(just_protein_coding_24hr$rRNAPolyA_1_24hr_Avg) / sum(human_annotation_24hr_v3$rRNAPolyA_1_24hr_Avg) * 100
human_pie_pc_24hr_rRNApolyA_10 <- sum(just_protein_coding_24hr$rRNAPolyA_10_24hr_Avg) / sum(human_annotation_24hr_v3$rRNAPolyA_10_24hr_Avg) * 100




#--
#Create pie charts of the % of protein coding expression between depletion methods
#--

#create metadata
pie_labels <- c("Protein coding", "Non-Protein coding")
pie_cols <- c("deepskyblue","firebrick2")

#--
#Plot the legend
#--
Cairo(file="PC_legend.png", type="png", units="in", width=5, 5=height, dpi=300, bg = "transparent")
cairo_pdf(file="PC_legend.pdf", width=width, height=height, bg = "transparent")
pie(c(human_pie_pc_1hr_rRNAdep_01, 100-human_pie_pc_1hr_rRNAdep_01), labels = "", col = c("deepskyblue","firebrick2"), border = "gray18", radius = 0.1)
legend(0.25,0.30, pie_labels, cex = 0.8, fill = pie_cols, bty = "n")
dev.off()


#Make plot function
plot_pies_protein_coding <- function (filename, pie_value) {
  #Cairo(file=filename, type="png", units="in", width=5, height=5, dpi=300, bg = "transparent")
  cairo_pdf(file=filename, width=5, height=5, bg = "transparent")
  pie(c(pie_value, 100-pie_value), 
      labels = c(paste0(round(pie_value,digits = 1),"%"),paste0(round(100-pie_value, digits = 1),"%")), 
      col = c("deepskyblue","firebrick2"), border = "gray18")
  dev.off()  
}


#--
#Pie plots
#--
#1hr
plot_pies_protein_coding(filename ="PC_rRNAPolyA_01_1hr.pdf", pie_value = human_pie_pc_1hr_rRNApolyA_01)
plot_pies_protein_coding(filename ="PC_rRNAPolyA_1_1hr.pdf", pie_value = human_pie_pc_1hr_rRNApolyA_1)
plot_pies_protein_coding(filename ="PC_rRNAPolyA_10_1hr.pdf", pie_value = human_pie_pc_1hr_rRNApolyA_10)
plot_pies_protein_coding(filename ="PC_rRNAdep_01_1hr.pdf", pie_value = human_pie_pc_1hr_rRNAdep_01)
plot_pies_protein_coding(filename ="PC_rRNAdep_1_1hr.pdf", pie_value = human_pie_pc_1hr_rRNAdep_1)
plot_pies_protein_coding(filename ="PC_rRNAdep_10_1hr.pdf", pie_value = human_pie_pc_1hr_rRNAdep_10)

#24hr
plot_pies_protein_coding(filename ="PC_rRNAPolyA_01_24hr.pdf", pie_value = human_pie_pc_24hr_rRNApolyA_01)
plot_pies_protein_coding(filename ="PC_rRNAPolyA_1_24hr.pdf", pie_value = human_pie_pc_24hr_rRNApolyA_1)
plot_pies_protein_coding(filename ="PC_rRNAPolyA_10_24hr.pdf", pie_value = human_pie_pc_24hr_rRNApolyA_10)
plot_pies_protein_coding(filename ="PC_rRNAdep_01_24hr.pdf", pie_value = human_pie_pc_24hr_rRNAdep_01)
plot_pies_protein_coding(filename ="PC_rRNAdep_1_24hr.pdf", pie_value = human_pie_pc_24hr_rRNAdep_1)
plot_pies_protein_coding(filename ="PC_rRNAdep_10_24hr.pdf", pie_value = human_pie_pc_24hr_rRNAdep_10)



#**************************************
#
# Beeswarm plot of the bio-types of non-protein coding genes appearing ----
#
#**************************************

table(human_annotation_1hr_v3$GENEBIOTYPE)
table(human_annotation_24hr_v3$GENEBIOTYPE)

#Make a master list of the types of genes appearing
unique(c(names(table(human_annotation_1hr_v3$GENEBIOTYPE)),names(table(human_annotation_24hr_v3$GENEBIOTYPE)))) #There are 25

#Make a new df just for 1 hr
colnames(human_annotation_1hr_v3)
biotype_1hr <- human_annotation_1hr_v3[,c(21:26,20)]

#Make a new df just for 24 hr
colnames(human_annotation_24hr_v3)
biotype_24hr <- human_annotation_24hr_v3[,c(21:26,20)]

#sum up all the gene types - 1hr
biotype_1hr_v2 <- aggregate(cbind(biotype_1hr$rRNAdep_01_1hr_Avg,biotype_1hr$rRNAPolyA_01_1hr_Avg,biotype_1hr$rRNAdep_1_1hr_Avg,
                                  biotype_1hr$rRNAPolyA_1_1hr_Avg,biotype_1hr$rRNAdep_10_1hr_Avg,biotype_1hr$rRNAPolyA_10_1hr_Avg), 
                            by=list(Biotype=biotype_1hr$GENEBIOTYPE), FUN=sum)
colnames(biotype_1hr_v2) <- colnames(human_annotation_1hr_v3)[c(20:26)]

#sum up all the gene types - 24hr
biotype_24hr_v2 <- aggregate(cbind(biotype_24hr$rRNAdep_01_24hr_Avg,biotype_24hr$rRNAPolyA_01_24hr_Avg,biotype_24hr$rRNAdep_1_24hr_Avg,
                                   biotype_24hr$rRNAPolyA_1_24hr_Avg,biotype_24hr$rRNAdep_10_24hr_Avg,biotype_24hr$rRNAPolyA_10_24hr_Avg), 
                             by=list(Biotype=biotype_24hr$GENEBIOTYPE), FUN=sum)
colnames(biotype_24hr_v2) <- colnames(human_annotation_24hr_v3)[c(20:26)]


#Get the %s

#1hr
biotype_1hr_v3 <- mutate(biotype_1hr_v2, 
                         rRNAdep_01_1hr_Avg_pct = rRNAdep_01_1hr_Avg / sum(rRNAdep_01_1hr_Avg),
                         rRNAdep_1_1hr_Avg_pct = rRNAdep_1_1hr_Avg / sum(rRNAdep_1_1hr_Avg),
                         rRNAdep_10_1hr_Avg_pct = rRNAdep_10_1hr_Avg / sum(rRNAdep_10_1hr_Avg),
                         rRNAPolyA_01_1hr_Avg_pct = rRNAPolyA_01_1hr_Avg / sum(rRNAPolyA_01_1hr_Avg),
                         rRNAPolyA_1_1hr_Avg_pct = rRNAPolyA_1_1hr_Avg / sum(rRNAPolyA_1_1hr_Avg),
                         rRNAPolyA_10_1hr_Avg_pct = rRNAPolyA_10_1hr_Avg / sum(rRNAPolyA_10_1hr_Avg)
)

biotype_1hr_v4 <- biotype_1hr_v3[,c(1,8:13)]


#24hr
biotype_24hr_v3 <- mutate(biotype_24hr_v2, 
                          rRNAdep_01_24hr_Avg_pct = rRNAdep_01_24hr_Avg / sum(rRNAdep_01_24hr_Avg),
                          rRNAdep_1_24hr_Avg_pct = rRNAdep_1_24hr_Avg / sum(rRNAdep_1_24hr_Avg),
                          rRNAdep_10_24hr_Avg_pct = rRNAdep_10_24hr_Avg / sum(rRNAdep_10_24hr_Avg),
                          rRNAPolyA_01_24hr_Avg_pct = rRNAPolyA_01_24hr_Avg / sum(rRNAPolyA_01_24hr_Avg),
                          rRNAPolyA_1_24hr_Avg_pct = rRNAPolyA_1_24hr_Avg / sum(rRNAPolyA_1_24hr_Avg),
                          rRNAPolyA_10_24hr_Avg_pct = rRNAPolyA_10_24hr_Avg / sum(rRNAPolyA_10_24hr_Avg)
)

biotype_24hr_v4 <- biotype_24hr_v3[,c(1,8:13)]




#Plot

#prepare data - 1hr
biotype_1hr_v5 <- data.frame(t(biotype_1hr_v4[,2:7]))
colnames(biotype_1hr_v5) <- biotype_1hr_v4$GENEBIOTYPE
biotype_1hr_v5$other <- biotype_1hr_v5$`3prime_overlapping_ncRNA` + biotype_1hr_v5$non_coding + biotype_1hr_v5$polymorphic_pseudogene + biotype_1hr_v5$rRNA + biotype_1hr_v5$miRNA + biotype_1hr_v5$macro_lncRNA
biotype_1hr_v5 <- biotype_1hr_v5[,c(13,2,3,6:8,11:12,15:25,27,26)] #remove cols combined into other
biotype_1hr_v5$cols <- unique(plot_cols)[1:6]

colnames(biotype_1hr_v5)
#prepare data - 24hr
biotype_24hr_v5 <- data.frame(t(biotype_24hr_v4[,2:7]))
colnames(biotype_24hr_v5) <- biotype_24hr_v4$GENEBIOTYPE
biotype_24hr_v5$other <- biotype_24hr_v5$macro_lncRNA
biotype_24hr_v5 <- biotype_24hr_v5[,c(9,1:2,4:8,10:21)] #match that of 1hr
biotype_24hr_v5$cols <- unique(plot_cols)[1:6]

#Make sure the cols are now the same
colnames(biotype_24hr_v5) == colnames(biotype_1hr_v5)



#melt and plot
biotype_1hr_v7 <- melt(biotype_1hr_v5)
biotype_24hr_v7 <- melt(biotype_24hr_v5)

beeswarm(value ~ variable,data = biotype_1hr_v7, pch = c(16,16,16,16,16,16), las=2, pwcol=cols, method = 'swarm')
beeswarm(value ~ variable,data = biotype_24hr_v7, pch = 16, las=2, pwcol=cols, method = 'swarm')


#Now combine the data frames
biotype_1hr_v5
biotype_24hr_v5
#add in a shape variable
biotype_1hr_v6 <- biotype_1hr_v5
biotype_24hr_v6 <- biotype_24hr_v5
biotype_combined <- rbind(biotype_1hr_v6, biotype_24hr_v6)
#move cols around
biotype_combined_v2 <- biotype_combined[,c(1,5,12,4,3,2,13,8,7,9,6,10:11,14:21)]

biotype_combined_v3 <- melt(biotype_combined_v2)
#add in the shape
biotype_combined_v3$shape <- rep(rep(c(16,1),each=6),20)
biotype_combined_v3$shape <- rep(rep(c(16,1,15,0),each=3),20)

beeswarm(value ~ variable,data = biotype_combined_v3, pwpch = shape, las=2, pwcol=cols, method = 'swarm', 
         spacing=0.25, cex=0.8, side=0, main="Expression of gene biotypes", cex.axis=0.8, xlab="", ylab="percent (%) of expression")
legend("topright", legend = c("1hr MOI 0.1 rRNAdep","1hr MOI 0.1 rRNAdep + polyA",
                              "1hr MOI 1 rRNAdep","1hr MOI 1 rRNAdep + polyA",
                              "1hr MOI 10 rRNAdep","1hr MOI 10 rRNAdep + polyA",
                              "24hr MOI 0.1 rRNAdep","24hr MOI 0.1 rRNAdep + polyA",
                              "24hr MOI 1 rRNAdep","24hr MOI 1 rRNAdep + polyA",
                              "24hr MOI 10 rRNAdep","24hr MOI 10 rRNAdep + polyA"),
       col=unique(plot_cols)[c(1:6,1:6)], pch = c(16,16,16,1,1,1,15,15,15,0,0,0), cex = 0.5, horiz = F)



#Now remove the protein coding genes
biotype_combined_v4 <- biotype_combined_v3[c(13:240),]
beeswarm(value ~ variable,data = biotype_combined_v4, pwpch = shape, las=2, pwcol=cols, method = 'swarm', 
         spacing=0.25, cex=0.8, side=0, main="Expression of gene biotypes", cex.axis=0.8, xlab="", ylab="percent (%) of expression")
legend("topright", legend = c("1hr MOI 0.1 rRNAdep","1hr MOI 0.1 rRNAdep + polyA",
                              "1hr MOI 1 rRNAdep","1hr MOI 1 rRNAdep + polyA",
                              "1hr MOI 10 rRNAdep","1hr MOI 10 rRNAdep + polyA",
                              "24hr MOI 0.1 rRNAdep","24hr MOI 0.1 rRNAdep + polyA",
                              "24hr MOI 1 rRNAdep","24hr MOI 1 rRNAdep + polyA",
                              "24hr MOI 10 rRNAdep","24hr MOI 10 rRNAdep + polyA"),
       col=unique(plot_cols)[c(1:6,1:6)], pch = c(16,16,16,1,1,1,15,15,15,0,0,0), cex = 0.5, horiz = F)



Cairo(file="beeswarm_all_biotypes_expression_v2.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
par(xpd=T, mar=par()$mar+c(8,0,0,0)) # Expand bottom 
beeswarm(value ~ variable,data = biotype_combined_v3, pwpch = shape, las=2, pwcol=cols, method = 'swarm', 
         spacing=0.25, cex=0.8, side=0, main="Expression of gene biotypes", cex.axis=0.8, xlab="", ylab="percent (%) of expression")
legend("topright", legend = c("1hr MOI 0.1 rRNAdep","1hr MOI 0.1 rRNAdep + polyA",
                              "1hr MOI 1 rRNAdep","1hr MOI 1 rRNAdep + polyA",
                              "1hr MOI 10 rRNAdep","1hr MOI 10 rRNAdep + polyA",
                              "24hr MOI 0.1 rRNAdep","24hr MOI 0.1 rRNAdep + polyA",
                              "24hr MOI 1 rRNAdep","24hr MOI 1 rRNAdep + polyA",
                              "24hr MOI 10 rRNAdep","24hr MOI 10 rRNAdep + polyA"),
       col=unique(plot_cols)[c(1:6,1:6)], pch = c(16,1,16,1,16,1,15,0,15,0,15,0), cex = 0.5, horiz = F)
dev.off()


Cairo(file="beeswarm_all_biotypes_expression_v3.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
cairo_pdf(file="beeswarm_all_biotypes_expression_v3.pdf", width=width, height=height, bg = "transparent")
par(xpd=T, mar=par()$mar+c(8,0,0,0)) # Expand bottom 
beeswarm(value ~ variable,data = biotype_combined_v4, pwpch = shape, las=2, pwcol=cols, method = 'swarm', 
         spacing=0.25, cex=0.8, side=0, main="Expression of gene biotypes", cex.axis=0.8, xlab="", ylab="percent (%) of expression")
legend("topright", legend = c("1hr MOI 0.1 rRNAdep","1hr MOI 0.1 rRNAdep + polyA",
                              "1hr MOI 1 rRNAdep","1hr MOI 1 rRNAdep + polyA",
                              "1hr MOI 10 rRNAdep","1hr MOI 10 rRNAdep + polyA",
                              "24hr MOI 0.1 rRNAdep","24hr MOI 0.1 rRNAdep + polyA",
                              "24hr MOI 1 rRNAdep","24hr MOI 1 rRNAdep + polyA",
                              "24hr MOI 10 rRNAdep","24hr MOI 10 rRNAdep + polyA"),
       col=unique(plot_cols)[c(1:6,1:6)], pch = c(16,1,16,1,16,1,15,0,15,0,15,0), cex = 0.5, horiz = F)
dev.off()










#**************************************
#
# Where are the top expressed non-coding genes appearing/and what type of nc genes (plot 1/2) ----
#
#**************************************

#Set up the variables to use and the gene information
#1hr
subset = 19:20
top_express_nc_1hr_moi_01_rRNADep <- human_annotation_1hr_v3[,c(subset,21)]
top_express_nc_1hr_moi_01_rRNApolyA <- human_annotation_1hr_v3[,c(subset,22)]
top_express_nc_1hr_moi_1_rRNADep <- human_annotation_1hr_v3[,c(subset,23)]
top_express_nc_1hr_moi_1_rRNApolyA <- human_annotation_1hr_v3[,c(subset,24)]
top_express_nc_1hr_moi_10_rRNADep <- human_annotation_1hr_v3[,c(subset,25)]
top_express_nc_1hr_moi_10_rRNApolyA <- human_annotation_1hr_v3[,c(subset,26)]

#24hr
top_express_nc_24hr_moi_01_rRNADep <- human_annotation_24hr_v3[,c(subset,21)]
top_express_nc_24hr_moi_01_rRNApolyA <- human_annotation_24hr_v3[,c(subset,22)]
top_express_nc_24hr_moi_1_rRNADep <- human_annotation_24hr_v3[,c(subset,23)]
top_express_nc_24hr_moi_1_rRNApolyA <- human_annotation_24hr_v3[,c(subset,24)]
top_express_nc_24hr_moi_10_rRNADep <- human_annotation_24hr_v3[,c(subset,25)]
top_express_nc_24hr_moi_10_rRNApolyA <- human_annotation_24hr_v3[,c(subset,26)]



#--
#order and get the top 200
#--
ex_top <- 200

#1hr
top_express_nc_1hr_moi_01_rRNApolyA_v2 <- top_express_nc_1hr_moi_01_rRNApolyA[order(top_express_nc_1hr_moi_01_rRNApolyA$rRNAPolyA_01_1hr_Avg, decreasing = T),][1:ex_top,]
top_express_nc_1hr_moi_01_rRNADep_v2 <- top_express_nc_1hr_moi_01_rRNADep[order(top_express_nc_1hr_moi_01_rRNADep$rRNAdep_01_1hr_Avg, decreasing = T),][1:ex_top,]
top_express_nc_1hr_moi_1_rRNApolyA_v2 <- top_express_nc_1hr_moi_1_rRNApolyA[order(top_express_nc_1hr_moi_1_rRNApolyA$rRNAPolyA_1_1hr_Avg, decreasing = T),][1:ex_top,]
top_express_nc_1hr_moi_1_rRNADep_v2 <- top_express_nc_1hr_moi_1_rRNADep[order(top_express_nc_1hr_moi_1_rRNADep$rRNAdep_1_1hr_Avg, decreasing = T),][1:ex_top,]
top_express_nc_1hr_moi_10_rRNApolyA_v2 <- top_express_nc_1hr_moi_10_rRNApolyA[order(top_express_nc_1hr_moi_10_rRNApolyA$rRNAPolyA_10_1hr_Avg, decreasing = T),][1:ex_top,]
top_express_nc_1hr_moi_10_rRNADep_v2 <- top_express_nc_1hr_moi_10_rRNADep[order(top_express_nc_1hr_moi_10_rRNADep$rRNAdep_10_1hr_Avg, decreasing = T),][1:ex_top,]
#24hr
top_express_nc_24hr_moi_01_rRNApolyA_v2 <- top_express_nc_24hr_moi_01_rRNApolyA[order(top_express_nc_24hr_moi_01_rRNApolyA$rRNAPolyA_01_24hr_Avg, decreasing = T),][1:ex_top,]
top_express_nc_24hr_moi_01_rRNADep_v2 <- top_express_nc_24hr_moi_01_rRNADep[order(top_express_nc_24hr_moi_01_rRNADep$rRNAdep_01_24hr_Avg, decreasing = T),][1:ex_top,]
top_express_nc_24hr_moi_1_rRNApolyA_v2 <- top_express_nc_24hr_moi_1_rRNApolyA[order(top_express_nc_24hr_moi_1_rRNApolyA$rRNAPolyA_1_24hr_Avg, decreasing = T),][1:ex_top,]
top_express_nc_24hr_moi_1_rRNADep_v2 <- top_express_nc_24hr_moi_1_rRNADep[order(top_express_nc_24hr_moi_1_rRNADep$rRNAdep_1_24hr_Avg, decreasing = T),][1:ex_top,]
top_express_nc_24hr_moi_10_rRNApolyA_v2 <- top_express_nc_24hr_moi_10_rRNApolyA[order(top_express_nc_24hr_moi_10_rRNApolyA$rRNAPolyA_10_24hr_Avg, decreasing = T),][1:ex_top,]
top_express_nc_24hr_moi_10_rRNADep_v2 <- top_express_nc_24hr_moi_10_rRNADep[order(top_express_nc_24hr_moi_10_rRNADep$rRNAdep_10_24hr_Avg, decreasing = T),][1:ex_top,]



#--
#Prepare plots
#--
width=7;height=5


#Set up the background removal
remove_bkg <- theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
edit_legend_format <- theme(plot.title = element_text(size = 6),legend.title=element_text(size=0), legend.text=element_text(size=7), legend.background = element_blank(), legend.key.height = unit(0.4,"cm")) 
edit_legend_position <- theme(legend.justification = c(1, 1), legend.position = c(1, 1))


#There are a few different ways of plotting this.
#This idea shows 2 plots and the top 200 genes from each - highlighting the non-protein coding genes

#Function to plot
plot_pc_expression_plot1 <- function (the_data_a, the_data_b, y_val_a, y_val_b, x_scale_a, x_scale_b, manual_cols_a, manual_cols_b, filename) {
  plota <- ggplot(data=the_data_a, aes_string(x="GENENAME", y=y_val_a)) + geom_point(aes(colour=GENEBIOTYPE)) + 
    scale_x_discrete(limits=x_scale_a) + scale_color_manual(values=manual_cols_a) +
    remove_bkg + edit_legend_format + edit_legend_position
  
  plotb <- ggplot(data=the_data_b, aes_string(x="GENENAME", y=y_val_b)) + geom_point(aes(colour=GENEBIOTYPE)) + 
    scale_x_discrete(limits=x_scale_b) + scale_color_manual(values=manual_cols_b) + 
    remove_bkg + edit_legend_format + edit_legend_position
  
  Cairo(file=filename, type="png", units="in", width=7, height=5, dpi=300, bg = "transparent")
  print(plot_grid(plota, plotb, labels = c("A", "B"), nrow = 2))
  dev.off()
  
}


#--
#plots
#--

#1hr MOI 0.1
plot_pc_expression_plot1(the_data_a = top_express_nc_1hr_moi_01_rRNApolyA_v2, y_val_a = "rRNAPolyA_01_1hr_Avg", x_scale_a = top_express_nc_1hr_moi_01_rRNApolyA_v2$GENENAME, 
                         manual_cols_a = c("red","green","orange","grey90","blue","cyan"),
                         the_data_b = top_express_nc_1hr_moi_01_rRNADep_v2, y_val_b = "rRNAdep_01_1hr_Avg", x_scale_b = top_express_nc_1hr_moi_01_rRNADep_v2$GENENAME, 
                         manual_cols_b = c("red","green","orange","grey90","blue"),
                         filename = "expression_diff_1hr_moi_01.png")

#1hr MOI 1
plot_pc_expression_plot1(the_data_a = top_express_nc_1hr_moi_1_rRNApolyA_v2, y_val_a = "rRNAPolyA_1_1hr_Avg", x_scale_a = top_express_nc_1hr_moi_1_rRNApolyA_v2$GENENAME, 
                         manual_cols_a = c("red","green","orange","grey90","blue","cyan"),
                         the_data_b = top_express_nc_1hr_moi_1_rRNADep_v2, y_val_b = "rRNAdep_1_1hr_Avg", x_scale_b = top_express_nc_1hr_moi_1_rRNADep_v2$GENENAME, 
                         manual_cols_b = c("red","green","orange","grey90","blue"),
                         filename = "expression_diff_1hr_moi_1.png")

#1hr MOI 10
plot_pc_expression_plot1(the_data_a = top_express_nc_1hr_moi_10_rRNApolyA_v2, y_val_a = "rRNAPolyA_10_1hr_Avg", x_scale_a = top_express_nc_1hr_moi_10_rRNApolyA_v2$GENENAME, 
                         manual_cols_a = c("red","green","orange","grey90","blue","cyan"),
                         the_data_b = top_express_nc_1hr_moi_10_rRNADep_v2, y_val_b = "rRNAdep_10_1hr_Avg", x_scale_b = top_express_nc_1hr_moi_10_rRNADep_v2$GENENAME, 
                         manual_cols_b = c("red","green","grey90","blue"),
                         filename = "expression_diff_1hr_moi_10.png")


#24hr MOI 0.1
plot_pc_expression_plot1(the_data_a = top_express_nc_24hr_moi_01_rRNApolyA_v2, y_val_a = "rRNAPolyA_01_24hr_Avg", x_scale_a = top_express_nc_24hr_moi_01_rRNApolyA_v2$GENENAME, 
                         manual_cols_a = c("red","green","orange","grey90","blue","cyan"),
                         the_data_b = top_express_nc_24hr_moi_01_rRNADep_v2, y_val_b = "rRNAdep_01_24hr_Avg", x_scale_b = top_express_nc_24hr_moi_01_rRNADep_v2$GENENAME, 
                         manual_cols_b = c("red","green","orange","grey90","blue"),
                         filename = "expression_diff_24hr_moi_01.png")

#24hr MOI 1
plot_pc_expression_plot1(the_data_a = top_express_nc_24hr_moi_1_rRNApolyA_v2, y_val_a = "rRNAPolyA_1_24hr_Avg", x_scale_a = top_express_nc_24hr_moi_1_rRNApolyA_v2$GENENAME, 
                         manual_cols_a = c("red","green","orange","brown","grey90","blue","cyan"),
                         the_data_b = top_express_nc_24hr_moi_1_rRNADep_v2, y_val_b = "rRNAdep_1_24hr_Avg", x_scale_b = top_express_nc_24hr_moi_1_rRNADep_v2$GENENAME, 
                         manual_cols_b = c("red","green","orange","grey90","blue"),
                         filename = "expression_diff_24hr_moi_1.png")

#24hr MOI 10
plot_pc_expression_plot1(the_data_a = top_express_nc_24hr_moi_10_rRNApolyA_v2, y_val_a = "rRNAPolyA_10_24hr_Avg", x_scale_a = top_express_nc_24hr_moi_10_rRNApolyA_v2$GENENAME, 
                         manual_cols_a = c("red","green","orange","grey90","blue","cyan","brown"),
                         the_data_b = top_express_nc_24hr_moi_10_rRNADep_v2, y_val_b = "rRNAdep_10_24hr_Avg", x_scale_b = top_express_nc_24hr_moi_10_rRNADep_v2$GENENAME, 
                         manual_cols_b = c("red","green","magenta","grey90","blue"),
                         filename = "expression_diff_24hr_moi_10.png")







#**************************************
#
# Where are the top expressed non-coding genes appearing/and what type of nc genes (plot 2/2) ----
#
#**************************************

#Using rRNADep as the reference genes, get the expression values from rRNAPolyA

#1hr
vv_1hr_moi_01_v1 <- merge(x=top_express_nc_1hr_moi_01_rRNApolyA_v2, y=top_express_nc_1hr_moi_01_rRNADep, by.x="GENENAME", by.y="GENENAME", all.x=T)
vv_1hr_moi_1_v1 <- merge(x=top_express_nc_1hr_moi_1_rRNApolyA_v2, y=top_express_nc_1hr_moi_1_rRNADep, by.x="GENENAME", by.y="GENENAME", all.x=T)
vv_1hr_moi_10_v1 <- merge(x=top_express_nc_1hr_moi_10_rRNApolyA_v2, y=top_express_nc_1hr_moi_10_rRNADep, by.x="GENENAME", by.y="GENENAME", all.x=T)
#24hr
vv_24hr_moi_01_v1 <- merge(x=top_express_nc_24hr_moi_01_rRNApolyA_v2, y=top_express_nc_24hr_moi_01_rRNADep, by.x="GENENAME", by.y="GENENAME", all.x=T)
vv_24hr_moi_1_v1 <- merge(x=top_express_nc_24hr_moi_1_rRNApolyA_v2, y=top_express_nc_24hr_moi_1_rRNADep, by.x="GENENAME", by.y="GENENAME", all.x=T)
vv_24hr_moi_10_v1 <- merge(x=top_express_nc_24hr_moi_10_rRNApolyA_v2, y=top_express_nc_24hr_moi_10_rRNADep, by.x="GENENAME", by.y="GENENAME", all.x=T)


#--
#Just keep the cols of interest
#--

#1hr
cols_int <-  c(1,2,3,5)
vv_1hr_moi_01_v2 <- vv_1hr_moi_01_v1[,cols_int]
vv_1hr_moi_1_v2 <- vv_1hr_moi_1_v1[,cols_int]
vv_1hr_moi_10_v2 <- vv_1hr_moi_10_v1[,cols_int]
#24hr
vv_24hr_moi_01_v2 <- vv_24hr_moi_01_v1[,cols_int]
vv_24hr_moi_1_v2 <- vv_24hr_moi_1_v1[,cols_int]
vv_24hr_moi_10_v2 <- vv_24hr_moi_10_v1[,cols_int]

#--
#Melt into ggplot format
#--

#1hr
vv_1hr_moi_01_v3 <- melt(vv_1hr_moi_01_v2)
vv_1hr_moi_1_v3 <- melt(vv_1hr_moi_1_v2)
vv_1hr_moi_10_v3 <- melt(vv_1hr_moi_10_v2)
#24hr
vv_24hr_moi_01_v3 <- melt(vv_24hr_moi_01_v2)
vv_24hr_moi_1_v3 <- melt(vv_24hr_moi_1_v2)
vv_24hr_moi_10_v3 <- melt(vv_24hr_moi_10_v2)


#--
#Prepare the plots
#--

#Set up the background removal
remove_bkg <- theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
edit_legend_format <- theme(plot.title = element_text(size = 5),legend.title=element_text(size=1), legend.text=element_text(size=10), legend.background = element_blank(), legend.key.height = unit(0.4,"cm")) 
edit_legend_position <- theme(legend.justification = c(1, 1), legend.position = c(1, 1))

#Note
#These seem a little more difficult to look at than plot 1

#Create function
plot_pc_expression_plot2 <- function (filename, the_data, x_limit, col_manual, size_manual) {
  Cairo(file=filename, type="png", units="in", width=5, height=4, dpi=300, bg = "transparent")
  print(ggplot(data=the_data, aes(x=GENENAME, y=value))  + geom_line(aes(colour=variable, group=variable)) + 
          geom_point(aes(colour=GENEBIOTYPE.x, size=GENEBIOTYPE.x), show.legend = T) +
          scale_x_discrete(limits=x_limit) +  
          scale_color_manual(values=col_manual) + 
          scale_size_manual(values=size_manual) +
          remove_bkg + edit_legend_format + edit_legend_position)
  dev.off()  
}

#1hr
plot_pc_expression_plot2(filename = "expression_diff_of_biotypes_1hr_moi_01.png", the_data = vv_1hr_moi_01_v3, x_limit = top_express_nc_1hr_moi_01_rRNApolyA_v2$GENENAME,
                         col_manual = c("rRNAPolyA_01_1hr_Avg"="grey","rRNAdep_01_1hr_Avg"="black","protein_coding"="grey","snoRNA"="blue","snRNA"="cyan","Mt_rRNA"="red","lincRNA"="orange","misc_RNA"="pink"),
                         size_manual = c("rRNAPolyA_01_1hr_Avg"=0.5,"rRNAdep_01_1hr_Avg"=0.5,"protein_coding"=0.1,"snoRNA"=3,"snRNA"=3,"Mt_rRNA"=3,"lincRNA"=3,"misc_RNA"=3))

plot_pc_expression_plot2(filename = "expression_diff_of_biotypes_1hr_moi_1.png", the_data = vv_1hr_moi_1_v3, x_limit = top_express_nc_1hr_moi_1_rRNApolyA_v2$GENENAME,
                         col_manual = c("rRNAPolyA_1_1hr_Avg"="grey","rRNAdep_1_1hr_Avg"="black","protein_coding"="grey","snoRNA"="blue","snRNA"="cyan","Mt_rRNA"="red","lincRNA"="orange","misc_RNA"="pink"),
                         size_manual = c("rRNAPolyA_1_1hr_Avg"=0.5,"rRNAdep_1_1hr_Avg"=0.5,"protein_coding"=0.1,"snoRNA"=3,"snRNA"=3,"Mt_rRNA"=3,"lincRNA"=3,"misc_RNA"=3))

plot_pc_expression_plot2(filename = "expression_diff_of_biotypes_1hr_moi_10.png", the_data = vv_1hr_moi_10_v3, x_limit = top_express_nc_1hr_moi_10_rRNApolyA_v2$GENENAME,
                         col_manual = c("rRNAPolyA_10_1hr_Avg"="grey","rRNAdep_10_1hr_Avg"="black","protein_coding"="grey","snoRNA"="blue","snRNA"="cyan","Mt_rRNA"="red","lincRNA"="orange","misc_RNA"="pink"),
                         size_manual = c("rRNAPolyA_10_1hr_Avg"=0.5,"rRNAdep_10_1hr_Avg"=0.5,"protein_coding"=0.1,"snoRNA"=3,"snRNA"=3,"Mt_rRNA"=3,"lincRNA"=3,"misc_RNA"=3))

#24hr
plot_pc_expression_plot2(filename = "expression_diff_of_biotypes_24hr_moi_01.png", the_data = vv_24hr_moi_01_v3, x_limit = top_express_nc_24hr_moi_01_rRNApolyA_v2$GENENAME,
                         col_manual = c("rRNAPolyA_01_24hr_Avg"="grey","rRNAdep_01_24hr_Avg"="black","protein_coding"="grey","snoRNA"="blue","snRNA"="cyan","Mt_rRNA"="red","lincRNA"="orange","misc_RNA"="pink"),
                         size_manual = c("rRNAPolyA_01_24hr_Avg"=0.5,"rRNAdep_01_24hr_Avg"=0.5,"protein_coding"=0.1,"snoRNA"=3,"snRNA"=3,"Mt_rRNA"=3,"lincRNA"=3,"misc_RNA"=3))

plot_pc_expression_plot2(filename = "expression_diff_of_biotypes_24hr_moi_1.png", the_data = vv_24hr_moi_1_v3, x_limit = top_express_nc_24hr_moi_1_rRNApolyA_v2$GENENAME,
                         col_manual = c("rRNAPolyA_1_24hr_Avg"="grey","rRNAdep_1_24hr_Avg"="black","protein_coding"="grey","snoRNA"="blue","snRNA"="cyan","Mt_rRNA"="red","Mt_tRNA"="brown","lincRNA"="orange","misc_RNA"="pink"),
                         size_manual = c("rRNAPolyA_1_24hr_Avg"=0.5,"rRNAdep_1_24hr_Avg"=0.5,"protein_coding"=0.1,"snoRNA"=3,"snRNA"=3,"Mt_rRNA"=3,"Mt_tRNA"=3,"lincRNA"=3,"misc_RNA"=3))

plot_pc_expression_plot2(filename = "expression_diff_of_biotypes_24hr_moi_10.png", the_data = vv_24hr_moi_10_v3, x_limit = top_express_nc_24hr_moi_10_rRNApolyA_v2$GENENAME,
                         col_manual = c("rRNAPolyA_10_24hr_Avg"="grey","rRNAdep_10_24hr_Avg"="black","protein_coding"="grey","snoRNA"="blue","snRNA"="cyan","Mt_rRNA"="red","scaRNA"="magenta","lincRNA"="orange","misc_RNA"="pink"),
                         size_manual = c("rRNAPolyA_10_24hr_Avg"=0.5,"rRNAdep_10_24hr_Avg"=0.5,"protein_coding"=0.1,"snoRNA"=3,"snRNA"=3,"Mt_rRNA"=3,"scaRNA"=3,"lincRNA"=3,"misc_RNA"=3))










#**************************************
#
# What are the actual top expressed non-coding genes ----
#
#**************************************

#Can identify the genes names from each of the six graphs above

#1hr
vv_1hr_moi_01_v2[order(vv_1hr_moi_01_v2$rRNAPolyA_01_1hr_Avg, decreasing = T),][c(1:6),]
vv_1hr_moi_1_v2[order(vv_1hr_moi_1_v2$rRNAPolyA_1_1hr_Avg, decreasing = T),][c(1:6),]
vv_1hr_moi_10_v2[order(vv_1hr_moi_10_v2$rRNAPolyA_10_1hr_Avg, decreasing = T),][c(1:4,15),]

#24
vv_24hr_moi_01_v2[order(vv_24hr_moi_01_v2$rRNAPolyA_01_24hr_Avg, decreasing = T),][c(1:6),]
vv_24hr_moi_1_v2[order(vv_24hr_moi_1_v2$rRNAPolyA_1_24hr_Avg, decreasing = T),][c(1:5,7,8),]
vv_24hr_moi_10_v2[order(vv_24hr_moi_10_v2$rRNAPolyA_10_24hr_Avg, decreasing = T),][c(1:6),]


#--
#Save the full lists down
#--

write.table(vv_1hr_moi_01_v2[order(vv_1hr_moi_01_v2$rRNAPolyA_01_1hr_Avg, decreasing = T),], "1hr_top_nc_expression_moi_01.csv", sep = ",", row.names = F, col.names = T)
write.table(vv_1hr_moi_1_v2[order(vv_1hr_moi_1_v2$rRNAPolyA_1_1hr_Avg, decreasing = T),], "1hr_top_nc_expression_moi_1.csv", sep = ",", row.names = F, col.names = T)
write.table(vv_1hr_moi_10_v2[order(vv_1hr_moi_10_v2$rRNAPolyA_10_1hr_Avg, decreasing = T),], "1hr_top_nc_expression_moi_10.csv", sep = ",", row.names = F, col.names = T)

write.table(vv_24hr_moi_01_v2[order(vv_24hr_moi_01_v2$rRNAPolyA_01_24hr_Avg, decreasing = T),], "24hr_top_nc_expression_moi_01.csv", sep = ",", row.names = F, col.names = T)
write.table(vv_24hr_moi_1_v2[order(vv_24hr_moi_1_v2$rRNAPolyA_1_24hr_Avg, decreasing = T),], "24hr_top_nc_expression_moi_1.csv", sep = ",", row.names = F, col.names = T)
write.table(vv_24hr_moi_10_v2[order(vv_24hr_moi_10_v2$rRNAPolyA_10_24hr_Avg, decreasing = T),], "24hr_top_nc_expression_moi_10.csv", sep = ",", row.names = F, col.names = T)





#**************************************
#
#Find significant DE genes that increase  from 0.1, 1, 10 ----
#
#**************************************

#Are there any overlaps in general?
#Human (both DE tests)
length(intersect(row.names(DE_genes_human_1hr_1_vs_01),row.names(DE_genes_human_1hr_1_vs_10))) #118
intersect(row.names(DE_genes_human_1hr_1_vs_01),row.names(DE_genes_human_1hr_1_vs_10))


#What are the ranges of log fc's
boxplot(DE_genes_human_1hr_1_vs_01$table$logFC) 
boxplot(DE_genes_human_1hr_1_vs_10$table$logFC) 


#--
#How many genes appear with fold changes above/below 1, or even 2?
#--
#Human
#1hr
human_de_1hr_MOI_01_fc_1_up_v1 <- DE_genes_human_1hr_1_vs_01$table[DE_genes_human_1hr_1_vs_01$table$logFC >1,]
dim(human_de_1hr_MOI_01_fc_1_up_v1) #65
human_de_1hr_MOI_01_fc_1_dw_v1 <- DE_genes_human_1hr_1_vs_01$table[DE_genes_human_1hr_1_vs_01$table$logFC < -1,]
dim(human_de_1hr_MOI_01_fc_1_dw_v1) #14

human_de_1hr_MOI_01_fc_2_up_v1 <- DE_genes_human_1hr_1_vs_01$table[DE_genes_human_1hr_1_vs_01$table$logFC >2,]
dim(human_de_1hr_MOI_01_fc_2_up_v1) #9
human_de_1hr_MOI_01_fc_2_dw_v1 <- DE_genes_human_1hr_1_vs_01$table[DE_genes_human_1hr_1_vs_01$table$logFC < -2,]
dim(human_de_1hr_MOI_01_fc_2_dw_v1) #2

human_de_1hr_MOI_10_fc_1_up_v1 <- DE_genes_human_1hr_1_vs_10$table[DE_genes_human_1hr_1_vs_10$table$logFC >1,]
dim(human_de_1hr_MOI_10_fc_1_up_v1) #45
human_de_1hr_MOI_10_fc_1_dw_v1 <- DE_genes_human_1hr_1_vs_10$table[DE_genes_human_1hr_1_vs_10$table$logFC < -1,]
dim(human_de_1hr_MOI_10_fc_1_dw_v1) #250

human_de_1hr_MOI_10_fc_2_up_v1 <- DE_genes_human_1hr_1_vs_10$table[DE_genes_human_1hr_1_vs_10$table$logFC > 2,]
dim(human_de_1hr_MOI_10_fc_2_up_v1) #4
human_de_1hr_MOI_10_fc_2_dw_v1 <- DE_genes_human_1hr_1_vs_10$table[DE_genes_human_1hr_1_vs_10$table$logFC < -2,]
dim(human_de_1hr_MOI_10_fc_2_dw_v1) #39

#24hr
human_de_24hr_MOI_01_fc_1_up_v1 <- DE_genes_human_24hr_1_vs_01$table[DE_genes_human_24hr_1_vs_01$table$logFC >1,]
dim(human_de_24hr_MOI_01_fc_1_up_v1) #20
human_de_24hr_MOI_01_fc_1_dw_v1 <- DE_genes_human_24hr_1_vs_01$table[DE_genes_human_24hr_1_vs_01$table$logFC < -1,]
dim(human_de_24hr_MOI_01_fc_1_dw_v1) #2

human_de_24hr_MOI_01_fc_2_up_v1 <- DE_genes_human_24hr_1_vs_01$table[DE_genes_human_24hr_1_vs_01$table$logFC >2,]
dim(human_de_24hr_MOI_01_fc_2_up_v1) #1
human_de_24hr_MOI_01_fc_2_dw_v1 <- DE_genes_human_24hr_1_vs_01$table[DE_genes_human_24hr_1_vs_01$table$logFC < -2,]
dim(human_de_24hr_MOI_01_fc_2_dw_v1) #0

human_de_24hr_MOI_10_fc_1_up_v1 <- DE_genes_human_24hr_1_vs_10$table[DE_genes_human_24hr_1_vs_10$table$logFC >1,]
dim(human_de_24hr_MOI_10_fc_1_up_v1) #336
human_de_24hr_MOI_10_fc_1_dw_v1 <- DE_genes_human_24hr_1_vs_10$table[DE_genes_human_24hr_1_vs_10$table$logFC < -1,]
dim(human_de_24hr_MOI_10_fc_1_dw_v1) #400

human_de_24hr_MOI_10_fc_2_up_v1 <- DE_genes_human_24hr_1_vs_10$table[DE_genes_human_24hr_1_vs_10$table$logFC > 2,]
dim(human_de_24hr_MOI_10_fc_2_up_v1) #35
human_de_24hr_MOI_10_fc_2_dw_v1 <- DE_genes_human_24hr_1_vs_10$table[DE_genes_human_24hr_1_vs_10$table$logFC < -2,]
dim(human_de_24hr_MOI_10_fc_2_dw_v1) #62


#--
#How many genes appear with fold changes above/below 1, or even 2?
#--
#CT
length(intersect(row.names(DE_genes_CT_1hr_1_vs_01),row.names(DE_genes_CT_1hr_1_vs_10))) #4
intersect(row.names(DE_genes_CT_1hr_1_vs_01),row.names(DE_genes_CT_1hr_1_vs_10))

#What are the ranges of log fc's
boxplot(DE_genes_CT_1hr_1_vs_01$table$logFC) 
boxplot(DE_genes_CT_1hr_1_vs_10$table$logFC)


#1hr
CT_de_1hr_MOI_01_fc_1_up_v1 <- DE_genes_CT_1hr_1_vs_01$table[DE_genes_CT_1hr_1_vs_01$table$logFC >1,]
dim(CT_de_1hr_MOI_01_fc_1_up_v1) #33
CT_de_1hr_MOI_01_fc_1_dw_v1 <- DE_genes_CT_1hr_1_vs_01$table[DE_genes_CT_1hr_1_vs_01$table$logFC < -1,]
dim(CT_de_1hr_MOI_01_fc_1_dw_v1) #6

CT_de_1hr_MOI_01_fc_2_up_v1 <- DE_genes_CT_1hr_1_vs_01$table[DE_genes_CT_1hr_1_vs_01$table$logFC >2,]
dim(CT_de_1hr_MOI_01_fc_2_up_v1) #12
CT_de_1hr_MOI_01_fc_2_dw_v1 <- DE_genes_CT_1hr_1_vs_01$table[DE_genes_CT_1hr_1_vs_01$table$logFC < -2,]
dim(CT_de_1hr_MOI_01_fc_2_dw_v1) #5

CT_de_1hr_MOI_10_fc_1_up_v1 <- DE_genes_CT_1hr_1_vs_10$table[DE_genes_CT_1hr_1_vs_10$table$logFC >1,]
dim(CT_de_1hr_MOI_10_fc_1_up_v1) #1
CT_de_1hr_MOI_10_fc_1_dw_v1 <- DE_genes_CT_1hr_1_vs_10$table[DE_genes_CT_1hr_1_vs_10$table$logFC < -1,]
dim(CT_de_1hr_MOI_10_fc_1_dw_v1) #6

CT_de_1hr_MOI_10_fc_2_up_v1 <- DE_genes_CT_1hr_1_vs_10$table[DE_genes_CT_1hr_1_vs_10$table$logFC > 2,]
dim(CT_de_1hr_MOI_10_fc_2_up_v1) #1
CT_de_1hr_MOI_10_fc_2_dw_v1 <- DE_genes_CT_1hr_1_vs_10$table[DE_genes_CT_1hr_1_vs_10$table$logFC < -2,]
dim(CT_de_1hr_MOI_10_fc_2_dw_v1) #3

#24hr
CT_de_24hr_MOI_01_fc_1_up_v1 <- DE_genes_CT_24hr_1_vs_01$table[DE_genes_CT_24hr_1_vs_01$table$logFC >1,]
dim(CT_de_24hr_MOI_01_fc_1_up_v1) #1
CT_de_24hr_MOI_01_fc_1_dw_v1 <- DE_genes_CT_24hr_1_vs_01$table[DE_genes_CT_24hr_1_vs_01$table$logFC < -1,]
dim(CT_de_24hr_MOI_01_fc_1_dw_v1) #4

CT_de_24hr_MOI_01_fc_2_up_v1 <- DE_genes_CT_24hr_1_vs_01$table[DE_genes_CT_24hr_1_vs_01$table$logFC >2,]
dim(CT_de_24hr_MOI_01_fc_2_up_v1) #0
CT_de_24hr_MOI_01_fc_2_dw_v1 <- DE_genes_CT_24hr_1_vs_01$table[DE_genes_CT_24hr_1_vs_01$table$logFC < -2,]
dim(CT_de_24hr_MOI_01_fc_2_dw_v1) #1

CT_de_24hr_MOI_10_fc_1_up_v1 <- DE_genes_CT_24hr_1_vs_10$table[DE_genes_CT_24hr_1_vs_10$table$logFC >1,]
dim(CT_de_24hr_MOI_10_fc_1_up_v1) #6
CT_de_24hr_MOI_10_fc_1_dw_v1 <- DE_genes_CT_24hr_1_vs_10$table[DE_genes_CT_24hr_1_vs_10$table$logFC < -1,]
dim(CT_de_24hr_MOI_10_fc_1_dw_v1) #7

CT_de_24hr_MOI_10_fc_2_up_v1 <- DE_genes_CT_24hr_1_vs_10$table[DE_genes_CT_24hr_1_vs_10$table$logFC > 2,]
dim(CT_de_24hr_MOI_10_fc_2_up_v1) #0
CT_de_24hr_MOI_10_fc_2_dw_v1 <- DE_genes_CT_24hr_1_vs_10$table[DE_genes_CT_24hr_1_vs_10$table$logFC < -2,]
dim(CT_de_24hr_MOI_10_fc_2_dw_v1) #0





#--
#Any overlaps where 0.1 = down and 10 =up? (trending up-reg)
#--

#Human
#1hr
intersect(row.names(human_de_1hr_MOI_01_fc_1_up_v1),row.names(human_de_1hr_MOI_10_fc_1_dw_v1))
length(intersect(row.names(human_de_1hr_MOI_01_fc_1_up_v1),row.names(human_de_1hr_MOI_10_fc_1_dw_v1))) #46 = trend of up regulation fc = 1

intersect(row.names(human_de_1hr_MOI_01_fc_1_dw_v1),row.names(human_de_1hr_MOI_10_fc_1_up_v1))
length(intersect(row.names(human_de_1hr_MOI_01_fc_1_dw_v1),row.names(human_de_1hr_MOI_10_fc_1_up_v1))) #0 = trend of down regulation fc = 1

intersect(row.names(human_de_1hr_MOI_01_fc_2_up_v1),row.names(human_de_1hr_MOI_10_fc_2_dw_v1))
length(intersect(row.names(human_de_1hr_MOI_01_fc_2_up_v1),row.names(human_de_1hr_MOI_10_fc_2_dw_v1))) #7 = trend of up regulation fc = 2

intersect(row.names(human_de_1hr_MOI_01_fc_2_dw_v1),row.names(human_de_1hr_MOI_10_fc_2_up_v1))
length(intersect(row.names(human_de_1hr_MOI_01_fc_2_dw_v1),row.names(human_de_1hr_MOI_10_fc_2_up_v1))) #0 = trend of down regulation fc = 2

#24hr
intersect(row.names(human_de_24hr_MOI_01_fc_1_up_v1),row.names(human_de_24hr_MOI_10_fc_1_dw_v1))
length(intersect(row.names(human_de_24hr_MOI_01_fc_1_up_v1),row.names(human_de_24hr_MOI_10_fc_1_dw_v1))) #14 = trend of up regulation fc = 1

intersect(row.names(human_de_24hr_MOI_01_fc_1_dw_v1),row.names(human_de_24hr_MOI_10_fc_1_up_v1))
length(intersect(row.names(human_de_24hr_MOI_01_fc_1_dw_v1),row.names(human_de_24hr_MOI_10_fc_1_up_v1))) #1 = trend of down regulation fc = 1

intersect(row.names(human_de_24hr_MOI_01_fc_2_up_v1),row.names(human_de_24hr_MOI_10_fc_2_dw_v1))
length(intersect(row.names(human_de_24hr_MOI_01_fc_2_up_v1),row.names(human_de_24hr_MOI_10_fc_2_dw_v1))) #1 = trend of up regulation fc = 2

intersect(row.names(human_de_24hr_MOI_01_fc_2_dw_v1),row.names(human_de_24hr_MOI_10_fc_2_up_v1))
length(intersect(row.names(human_de_24hr_MOI_01_fc_2_dw_v1),row.names(human_de_24hr_MOI_10_fc_2_up_v1))) #0 = trend of down regulation fc = 2

#CT
#1hr
intersect(row.names(CT_de_1hr_MOI_01_fc_1_up_v1),row.names(CT_de_1hr_MOI_10_fc_1_dw_v1))
length(intersect(row.names(CT_de_1hr_MOI_01_fc_1_up_v1),row.names(CT_de_1hr_MOI_10_fc_1_dw_v1))) #0 = trend of up regulation fc = 1

intersect(row.names(CT_de_1hr_MOI_01_fc_1_dw_v1),row.names(CT_de_1hr_MOI_10_fc_1_up_v1))
length(intersect(row.names(CT_de_1hr_MOI_01_fc_1_dw_v1),row.names(CT_de_1hr_MOI_10_fc_1_up_v1))) #1 = trend of down regulation fc = 1

#24hr
intersect(row.names(CT_de_24hr_MOI_01_fc_1_up_v1),row.names(CT_de_24hr_MOI_10_fc_1_dw_v1))
length(intersect(row.names(CT_de_24hr_MOI_01_fc_1_up_v1),row.names(CT_de_24hr_MOI_10_fc_1_dw_v1))) #0 = trend of up regulation fc = 1

intersect(row.names(CT_de_24hr_MOI_01_fc_1_dw_v1),row.names(CT_de_24hr_MOI_10_fc_1_up_v1))
length(intersect(row.names(CT_de_24hr_MOI_01_fc_1_dw_v1),row.names(CT_de_24hr_MOI_10_fc_1_up_v1))) #1 = trend of down regulation fc = 1





#**************************************
#
#Trended MOI regulation plots ----
#
#**************************************

#--
#Trending up plots
#--

dim(human_rle_plots_norm_1hr)

#Human
human_1hr_trend_genes_v1 <- as.data.frame(human_rle_plots_norm_1hr)
human_24hr_trend_genes_v1 <- as.data.frame(human_rle_plots_norm_24hr)

human_1hr_trend_genes_v1$GENE_NAME <- row.names(human_1hr_trend_genes_v1)
human_24hr_trend_genes_v1$GENE_NAME <- row.names(human_24hr_trend_genes_v1)

#CT
CT_1hr_trend_genes_v1 <- as.data.frame(CT_rle_plots_norm_1hr)
CT_24hr_trend_genes_v1 <- as.data.frame(CT_rle_plots_norm_24hr)

CT_1hr_trend_genes_v1$GENE_NAME <- row.names(CT_1hr_trend_genes_v1)
CT_24hr_trend_genes_v1$GENE_NAME <- row.names(CT_24hr_trend_genes_v1)



#Extract out the genes that are trended up (fc 1)
#Human
human_1hr_trend_genes_up_v1 <- human_1hr_trend_genes_v1[human_1hr_trend_genes_v1$GENE_NAME %in% intersect(row.names(human_de_1hr_MOI_01_fc_1_up_v1),row.names(human_de_1hr_MOI_10_fc_1_dw_v1)),]
rownames(human_1hr_trend_genes_up_v1) <- c()

human_24hr_trend_genes_up_v1 <- human_24hr_trend_genes_v1[human_24hr_trend_genes_v1$GENE_NAME %in% intersect(row.names(human_de_24hr_MOI_01_fc_1_up_v1),row.names(human_de_24hr_MOI_10_fc_1_dw_v1)),]
rownames(human_24hr_trend_genes_up_v1) <- c()

#Extract out the genes that are trended down (fc 1)
#Human
human_24hr_trend_genes_dw_v1 <- human_24hr_trend_genes_v1[human_24hr_trend_genes_v1$GENE_NAME %in% intersect(row.names(human_de_24hr_MOI_01_fc_1_dw_v1),row.names(human_de_24hr_MOI_10_fc_1_up_v1)),]
rownames(human_24hr_trend_genes_dw_v1) <- c()

#CT
CT_1hr_trend_genes_dw_v1 <- CT_1hr_trend_genes_v1[CT_1hr_trend_genes_v1$GENE_NAME %in% intersect(row.names(CT_de_1hr_MOI_01_fc_1_dw_v1),row.names(CT_de_1hr_MOI_10_fc_1_up_v1)),]
rownames(CT_1hr_trend_genes_dw_v1) <- c()

CT_24hr_trend_genes_dw_v1 <- CT_24hr_trend_genes_v1[CT_24hr_trend_genes_v1$GENE_NAME %in% intersect(row.names(CT_de_24hr_MOI_01_fc_1_dw_v1),row.names(CT_de_24hr_MOI_10_fc_1_up_v1)),]
rownames(CT_24hr_trend_genes_dw_v1) <- c()

#what are the two chlamydial genes
CT_features3[CT_features3$ID %in% CT_1hr_trend_genes_dw_v1$GENE_NAME,]
CT_features3[CT_features3$ID %in% CT_24hr_trend_genes_dw_v1$GENE_NAME,]



#--
#Make the graph df of the top genes
#--
#Human
#1hr up
human_1hr_trend_genes_up_df <- data.frame("Gene1" = as.numeric(human_1hr_trend_genes_up_v1[1,1:17]),
                                          "Gene2" = as.numeric(human_1hr_trend_genes_up_v1[2,1:17]),
                                          "Gene3" = as.numeric(human_1hr_trend_genes_up_v1[3,1:17]),
                                          "Gene4" = as.numeric(human_1hr_trend_genes_up_v1[4,1:17]),
                                          "MOI" = c(rep("0.1",6),rep("1",6),rep("10",5)),
                                          "Depletion" = rep(c(rep("rRNA",3),rep("rRNA+PolyA",3)),3)[1:17])

#24 up
human_24hr_trend_genes_up_df <- data.frame("Gene1" = as.numeric(human_24hr_trend_genes_up_v1[1,1:17]),
                                           "Gene2" = as.numeric(human_24hr_trend_genes_up_v1[2,1:17]),
                                           "Gene3" = as.numeric(human_24hr_trend_genes_up_v1[3,1:17]),
                                           "Gene4" = as.numeric(human_24hr_trend_genes_up_v1[4,1:17]),
                                           "MOI" = c(rep("0.1",6),rep("1",6),rep("10",5)),
                                           "Depletion" = rep(c(rep("rRNA",3),rep("rRNA+PolyA",3)),3)[1:17])

#24 down
human_24hr_trend_genes_dw_df <- data.frame("Gene1" = as.numeric(human_24hr_trend_genes_dw_v1[1,1:17]),
                                           "Gene2" = as.numeric(human_24hr_trend_genes_dw_v1[2,1:17]),
                                           "Gene3" = as.numeric(human_24hr_trend_genes_dw_v1[3,1:17]),
                                           "Gene4" = as.numeric(human_24hr_trend_genes_dw_v1[4,1:17]),
                                           "MOI" = c(rep("0.1",6),rep("1",6),rep("10",5)),
                                           "Depletion" = rep(c(rep("rRNA",3),rep("rRNA+PolyA",3)),3)[1:17])

#CT
#1 down
CT_1hr_trend_genes_dw_df <- data.frame("Gene1" = as.numeric(CT_1hr_trend_genes_dw_v1[1,1:18]),
                                       "MOI" = c(rep("0.1",6),rep("1",6),rep("10",6)),
                                       "Depletion" = rep(c(rep("rRNA",3),rep("rRNA+PolyA",3)),3))

#24 down
CT_24hr_trend_genes_dw_df <- data.frame("Gene1" = as.numeric(CT_24hr_trend_genes_dw_v1[1,1:16]),
                                        "MOI" = c(rep("0.1",4),rep("1",6),rep("10",6)),
                                        "Depletion" = rep(c(rep("rRNA",3),rep("rRNA+PolyA",3)),3)[c(1:3,6:18)])



#---
#The jitter plots
#---


#Plot graphs function
trended_genes <- function (filename, the_data, title, no_genes) {
  #As there are lots of genes, make function to cycle through them
  for(i in 1:no_genes){
    #The plot part
    Cairo(file=paste0(filename,"_gene_",i,".png"), type="png", units="in", width=6, height=4, dpi=300, bg = "white")
    print(ggplot(the_data, aes(MOI, paste0("Gene",i), colour=Depletion)) + 
            geom_jitter(width = 0.1,alpha=0.3, size=1.5) + scale_x_discrete(limits=c("0.1","1","10")) + 
            ggtitle(paste0(title,"_gene",i)))
    dev.off()
  }
}


#Human
trended_genes(filename = "Trended_up_human_1hr", the_data = human_1hr_trend_genes_up_df, title = "Trended up - Human - 1hr", no_genes = 10)
trended_genes(filename = "Trended_up_human_24hr", the_data = human_24hr_trend_genes_up_df, title = "Trended up - Human - 24hr", no_genes = 10)
trended_genes(filename = "Trended_down_human_24hr", the_data = human_24hr_trend_genes_dw_df, title = "Trended down - Human - 24hr", no_genes = 1)
#CT
trended_genes(filename = "Trended_down_CT_1hr", the_data = CT_1hr_trend_genes_dw_df, title = "Trended down - CT - 1hr", no_genes = 1)
trended_genes(filename = "Trended_down_CT_24hr", the_data = CT_24hr_trend_genes_dw_df, title = "Trended down - CT - 24hr", no_genes = 1)







#**************************************
#
# Annotate all of the DE genes (Human) and save down ----
#
#**************************************

#The annotation db
edb <- EnsDb.Hsapiens.v86 #connect
columns(edb) #The available columns

#Get the list of gene types
ens_cols <- c("GENEID","GENENAME","GENEBIOTYPE")

#The lists of genes to use
#1hr
dim(DE_genes_human_1hr_1_vs_01) #304
dim(DE_genes_human_1hr_1_vs_10) #1423
#24hr
dim(DE_genes_human_24hr_1_vs_01) #363
dim(DE_genes_human_24hr_1_vs_10) #4606

#human
#1hr
human_1hr_moi_1_vs_01_pathways_v1 <- select(edb,keys=as.vector(row.names(DE_genes_human_1hr_1_vs_01)),columns=ens_cols,keytype="GENEID")
human_1hr_moi_1_vs_10_pathways_v1 <- select(edb,keys=as.vector(row.names(DE_genes_human_1hr_1_vs_10)),columns=ens_cols,keytype="GENEID")
#24hr
human_24hr_moi_1_vs_01_pathways_v1 <- select(edb,keys=as.vector(row.names(DE_genes_human_24hr_1_vs_01)),columns=ens_cols,keytype="GENEID")
human_24hr_moi_1_vs_10_pathways_v1 <- select(edb,keys=as.vector(row.names(DE_genes_human_24hr_1_vs_10)),columns=ens_cols,keytype="GENEID")


#Get the gene description
columns(org.Hs.eg.db)

#1hr
human_1hr_moi_1_vs_01_pathways_v2 <- ensembldb::select(org.Hs.eg.db, keys = as.character(human_1hr_moi_1_vs_01_pathways_v1$GENENAME), columns=c("ENTREZID","GENENAME","ENSEMBL"), keytype="SYMBOL", multiVals = "first")
human_1hr_moi_1_vs_10_pathways_v2 <- ensembldb::select(org.Hs.eg.db, keys = as.character(human_1hr_moi_1_vs_10_pathways_v1$GENENAME), columns=c("ENTREZID","GENENAME","ENSEMBL"), keytype="SYMBOL", multiVals = "first")
#24hr
human_24hr_moi_1_vs_01_pathways_v2 <- ensembldb::select(org.Hs.eg.db, keys = as.character(human_24hr_moi_1_vs_01_pathways_v1$GENENAME), columns=c("ENTREZID","GENENAME","ENSEMBL"), keytype="SYMBOL", multiVals = "first")
human_24hr_moi_1_vs_10_pathways_v2 <- ensembldb::select(org.Hs.eg.db, keys = as.character(human_24hr_moi_1_vs_10_pathways_v1$GENENAME), columns=c("ENTREZID","GENENAME","ENSEMBL"), keytype="SYMBOL", multiVals = "first")

#remove duplicates
#1hr
human_1hr_moi_1_vs_01_pathways_v2 <- human_1hr_moi_1_vs_01_pathways_v2[!duplicated(human_1hr_moi_1_vs_01_pathways_v2),]
human_1hr_moi_1_vs_10_pathways_v2 <- human_1hr_moi_1_vs_10_pathways_v2[!duplicated(human_1hr_moi_1_vs_10_pathways_v2),]
#24hr
human_24hr_moi_1_vs_01_pathways_v2 <- human_24hr_moi_1_vs_01_pathways_v2[!duplicated(human_24hr_moi_1_vs_01_pathways_v2),]
human_24hr_moi_1_vs_10_pathways_v2 <- human_24hr_moi_1_vs_10_pathways_v2[!duplicated(human_24hr_moi_1_vs_10_pathways_v2),]


#now combine them
#1hr
human_1hr_moi_1_vs_01_pathways_v3 <- merge(x=human_1hr_moi_1_vs_01_pathways_v1, y=human_1hr_moi_1_vs_01_pathways_v2, by.x="GENEID", by.y="ENSEMBL", all.x=T)
human_1hr_moi_1_vs_10_pathways_v3 <- merge(x=human_1hr_moi_1_vs_10_pathways_v1, y=human_1hr_moi_1_vs_10_pathways_v2, by.x="GENEID", by.y="ENSEMBL", all.x=T)
dim(human_1hr_moi_1_vs_01_pathways_v1);dim(human_1hr_moi_1_vs_01_pathways_v3);dim(DE_genes_human_1hr_1_vs_01)
dim(human_1hr_moi_1_vs_10_pathways_v1);dim(human_1hr_moi_1_vs_10_pathways_v3);dim(DE_genes_human_1hr_1_vs_10)
#24hr
human_24hr_moi_1_vs_01_pathways_v3 <- merge(x=human_24hr_moi_1_vs_01_pathways_v1, y=human_24hr_moi_1_vs_01_pathways_v2, by.x="GENEID", by.y="ENSEMBL", all.x=T)
human_24hr_moi_1_vs_10_pathways_v3 <- merge(x=human_24hr_moi_1_vs_10_pathways_v1, y=human_24hr_moi_1_vs_10_pathways_v2, by.x="GENEID", by.y="ENSEMBL", all.x=T)
dim(human_24hr_moi_1_vs_01_pathways_v1);dim(human_24hr_moi_1_vs_01_pathways_v3);dim(DE_genes_human_24hr_1_vs_01)
dim(human_24hr_moi_1_vs_10_pathways_v1);dim(human_24hr_moi_1_vs_10_pathways_v3);dim(DE_genes_human_24hr_1_vs_10)


#Just keep the cols of interest
#1hr
human_1hr_moi_1_vs_01_pathways_v4 <- human_1hr_moi_1_vs_01_pathways_v3[,c(1,2,5,6,3)]
human_1hr_moi_1_vs_10_pathways_v4 <- human_1hr_moi_1_vs_10_pathways_v3[,c(1,2,5,6,3)]
#24hr
human_24hr_moi_1_vs_01_pathways_v4 <- human_24hr_moi_1_vs_01_pathways_v3[,c(1,2,5,6,3)]
human_24hr_moi_1_vs_10_pathways_v4 <- human_24hr_moi_1_vs_10_pathways_v3[,c(1,2,5,6,3)]


#Now add in the fc info etc
#1hr
DE_genes_human_1hr_1_vs_01_v2 <- as.data.frame(DE_genes_human_1hr_1_vs_01)
DE_genes_human_1hr_1_vs_10_v2 <- as.data.frame(DE_genes_human_1hr_1_vs_10)
DE_genes_human_1hr_1_vs_01_v2$Gene_Name <- row.names(DE_genes_human_1hr_1_vs_01_v2)
DE_genes_human_1hr_1_vs_10_v2$Gene_Name <- row.names(DE_genes_human_1hr_1_vs_10_v2)
#24hr
DE_genes_human_24hr_1_vs_01_v2 <- as.data.frame(DE_genes_human_24hr_1_vs_01)
DE_genes_human_24hr_1_vs_10_v2 <- as.data.frame(DE_genes_human_24hr_1_vs_10)
DE_genes_human_24hr_1_vs_01_v2$Gene_Name <- row.names(DE_genes_human_24hr_1_vs_01_v2)
DE_genes_human_24hr_1_vs_10_v2$Gene_Name <- row.names(DE_genes_human_24hr_1_vs_10_v2)

#1hr
human_1hr_moi_1_vs_01_pathways_v5 <- merge(x=human_1hr_moi_1_vs_01_pathways_v4, y=DE_genes_human_1hr_1_vs_01_v2, by.x="GENEID", by.y="Gene_Name", all.x=T)
human_1hr_moi_1_vs_10_pathways_v5 <- merge(x=human_1hr_moi_1_vs_10_pathways_v4, y=DE_genes_human_1hr_1_vs_10_v2, by.x="GENEID", by.y="Gene_Name", all.x=T)
#24hr
human_24hr_moi_1_vs_01_pathways_v5 <- merge(x=human_24hr_moi_1_vs_01_pathways_v4, y=DE_genes_human_24hr_1_vs_01_v2, by.x="GENEID", by.y="Gene_Name", all.x=T)
human_24hr_moi_1_vs_10_pathways_v5 <- merge(x=human_24hr_moi_1_vs_10_pathways_v4, y=DE_genes_human_24hr_1_vs_10_v2, by.x="GENEID", by.y="Gene_Name", all.x=T)

#keep cols of interest
#1hr
human_1hr_moi_1_vs_01_pathways_v6 <- human_1hr_moi_1_vs_01_pathways_v5[,c(1:6,10)]
dim(human_1hr_moi_1_vs_01_pathways_v6)
human_1hr_moi_1_vs_10_pathways_v6 <- human_1hr_moi_1_vs_10_pathways_v5[,c(1:6,10)]
dim(human_1hr_moi_1_vs_10_pathways_v6)
#24hr
human_24hr_moi_1_vs_01_pathways_v6 <- human_24hr_moi_1_vs_01_pathways_v5[,c(1:6,10)]
dim(human_24hr_moi_1_vs_01_pathways_v6)
human_24hr_moi_1_vs_10_pathways_v6 <- human_24hr_moi_1_vs_10_pathways_v5[,c(1:6,10)]
dim(human_24hr_moi_1_vs_10_pathways_v6)


#Order by FDR
#1hr
human_1hr_moi_1_vs_01_pathways_v7 <- human_1hr_moi_1_vs_01_pathways_v6[order(human_1hr_moi_1_vs_01_pathways_v6$FDR),]
human_1hr_moi_1_vs_10_pathways_v7 <- human_1hr_moi_1_vs_10_pathways_v6[order(human_1hr_moi_1_vs_10_pathways_v6$FDR),]
#24hr
human_24hr_moi_1_vs_01_pathways_v7 <- human_24hr_moi_1_vs_01_pathways_v6[order(human_24hr_moi_1_vs_01_pathways_v6$FDR),]
human_24hr_moi_1_vs_10_pathways_v7 <- human_24hr_moi_1_vs_10_pathways_v6[order(human_24hr_moi_1_vs_10_pathways_v6$FDR),]


#--
#Save down the lists
#--

#human
#1hr
write.table(human_1hr_moi_1_vs_01_pathways_v7, "human_1hr_moi_1_vs_01_DE_genes.csv", row.names = F, col.names = T, sep = ",")
write.table(human_1hr_moi_1_vs_10_pathways_v7, "human_1hr_moi_1_vs_10_DE_genes.csv", row.names = F, col.names = T, sep = ",")
#24hr
write.table(human_24hr_moi_1_vs_01_pathways_v7, "human_24hr_moi_1_vs_01_DE_genes.csv", row.names = F, col.names = T, sep = ",")
write.table(human_24hr_moi_1_vs_10_pathways_v7, "human_24hr_moi_1_vs_10_DE_genes.csv", row.names = F, col.names = T, sep = ",")






#**************************************
#
# Annotate the CT DE genes ----
#
#**************************************


#--
#Load in the chlamydial GFF (using rtracklayer)
#--
features <- import("charm_genome.gff")

#just keep the three main types
table(features$type)
CT_features2 <- features[(features$type =="CDS" | features$type =="tRNA" | features$type =="rRNA") ,c(1,2,5,6,12:16)]
head(features2)

#Make a df of just the cols of interest (need for vlookup)
CT_features3 <- data.frame("ID" = CT_features2$ID,
                           "Name" = CT_features2$Name,
                           "Parent" = CT_features2$Parent,
                           "Product" = CT_features2$product,
                           "Bio_type" = CT_features2$type,
                           "Length" = CT_features2@ranges@width)

head(CT_features3)

#removing the last four characters .t01, .p01 etc from the ID
CT_features3$ID <- substr(as.character(CT_features3$ID),1,nchar(as.character(CT_features3$ID))-4)
CT_features3 <- CT_features3[,c(1,2,6:8)]
head(CT_features3)
dim(CT_features3)

#--
#Convert into DF
#--
#1hr
DE_genes_CT_1hr_1_vs_01_v2 <- as.data.frame(DE_genes_CT_1hr_1_vs_01)
DE_genes_CT_1hr_1_vs_10_v2 <- as.data.frame(DE_genes_CT_1hr_1_vs_10)
DE_genes_CT_1hr_1_vs_01_v2$Gene_name <- row.names(DE_genes_CT_1hr_1_vs_01_v2)
DE_genes_CT_1hr_1_vs_10_v2$Gene_name <- row.names(DE_genes_CT_1hr_1_vs_10_v2)
#24hr
DE_genes_CT_24hr_1_vs_01_v2 <- as.data.frame(DE_genes_CT_24hr_1_vs_01)
DE_genes_CT_24hr_1_vs_10_v2 <- as.data.frame(DE_genes_CT_24hr_1_vs_10)
DE_genes_CT_24hr_1_vs_01_v2$Gene_name <- row.names(DE_genes_CT_24hr_1_vs_01_v2)
DE_genes_CT_24hr_1_vs_10_v2$Gene_name <- row.names(DE_genes_CT_24hr_1_vs_10_v2)

#--
#Merge annotation with DE genes
#--
#1hr
DE_genes_CT_1hr_1_vs_01_v3 <- merge(x=DE_genes_CT_1hr_1_vs_01_v2, y=CT_features3, by.x="Gene_name", by.y="ID", all.x=T)
DE_genes_CT_1hr_1_vs_10_v3 <- merge(x=DE_genes_CT_1hr_1_vs_10_v2, y=CT_features3, by.x="Gene_name", by.y="ID", all.x=T)
#24hr
DE_genes_CT_24hr_1_vs_01_v3 <- merge(x=DE_genes_CT_24hr_1_vs_01_v2, y=CT_features3, by.x="Gene_name", by.y="ID", all.x=T)
DE_genes_CT_24hr_1_vs_10_v3 <- merge(x=DE_genes_CT_24hr_1_vs_10_v2, y=CT_features3, by.x="Gene_name", by.y="ID", all.x=T)

#Keep cols of interest
#1hr
DE_genes_CT_1hr_1_vs_01_v4 <- DE_genes_CT_1hr_1_vs_01_v3[,c(1,7,8,9,10,2,6)]
DE_genes_CT_1hr_1_vs_10_v4 <- DE_genes_CT_1hr_1_vs_10_v3[,c(1,7,8,9,10,2,6)]
#24hr
DE_genes_CT_24hr_1_vs_01_v4 <- DE_genes_CT_24hr_1_vs_01_v3[,c(1,7,8,9,10,2,6)]
DE_genes_CT_24hr_1_vs_10_v4 <- DE_genes_CT_24hr_1_vs_10_v3[,c(1,7,8,9,10,2,6)]

#order by FDR
#1hr
DE_genes_CT_1hr_1_vs_01_v5 <- DE_genes_CT_1hr_1_vs_01_v4[order(DE_genes_CT_1hr_1_vs_01_v4$FDR, decreasing = F),]
DE_genes_CT_1hr_1_vs_10_v5 <- DE_genes_CT_1hr_1_vs_10_v4[order(DE_genes_CT_1hr_1_vs_10_v4$FDR, decreasing = F),]
#24hr
DE_genes_CT_24hr_1_vs_01_v5 <- DE_genes_CT_24hr_1_vs_01_v4[order(DE_genes_CT_24hr_1_vs_01_v4$FDR, decreasing = F),]
DE_genes_CT_24hr_1_vs_10_v5 <- DE_genes_CT_24hr_1_vs_10_v4[order(DE_genes_CT_24hr_1_vs_10_v4$FDR, decreasing = F),]

#--
#Save down
#--
#1hr
write.table(DE_genes_CT_1hr_1_vs_01_v5, "DE_genes_CT_1hr_1_vs_01.csv", row.names = F, col.names = T, sep = ",")
write.table(DE_genes_CT_1hr_1_vs_10_v5, "DE_genes_CT_1hr_1_vs_10.csv", row.names = F, col.names = T, sep = ",")
#24hr
write.table(DE_genes_CT_24hr_1_vs_01_v5, "DE_genes_CT_24hr_1_vs_01.csv", row.names = F, col.names = T, sep = ",")
write.table(DE_genes_CT_24hr_1_vs_10_v5, "DE_genes_CT_24hr_1_vs_10.csv", row.names = F, col.names = T, sep = ",")





#**************************************
#
# Link up the annotated list just created to the trended genes, then save ----
#
#**************************************

#--
#Human
#--
#1hr
dim(human_1hr_trend_genes_up_v1)[1] == length(intersect(row.names(human_de_1hr_MOI_01_fc_1_up_v1),row.names(human_de_1hr_MOI_10_fc_1_dw_v1))) #Checking the numbers match above
final_trended_up_genes_human_1hr_annotated_v1 <-  human_1hr_moi_1_vs_01_pathways_v7[human_1hr_moi_1_vs_01_pathways_v7$GENEID %in% 
                                                                                      intersect(row.names(human_de_1hr_MOI_01_fc_1_up_v1),row.names(human_de_1hr_MOI_10_fc_1_dw_v1)),c(1:5)]

#24hr
dim(human_24hr_trend_genes_up_v1)[1] == length(intersect(row.names(human_de_24hr_MOI_01_fc_1_up_v1),row.names(human_de_24hr_MOI_10_fc_1_dw_v1))) #Checking the numbers match above
final_trended_up_genes_human_24hr_annotated_v1 <-  human_24hr_moi_1_vs_01_pathways_v7[human_24hr_moi_1_vs_01_pathways_v7$GENEID %in% 
                                                                                        intersect(row.names(human_de_24hr_MOI_01_fc_1_up_v1),row.names(human_de_24hr_MOI_10_fc_1_dw_v1)),c(1:5)]

dim(human_24hr_trend_genes_dw_v1)[1] == length(intersect(row.names(human_de_24hr_MOI_01_fc_1_dw_v1),row.names(human_de_24hr_MOI_10_fc_1_up_v1))) #Checking the numbers match above
final_trended_down_genes_human_24hr_annotated_v1 <-  human_24hr_moi_1_vs_01_pathways_v7[human_24hr_moi_1_vs_01_pathways_v7$GENEID %in% 
                                                                                          intersect(row.names(human_de_24hr_MOI_01_fc_1_dw_v1),row.names(human_de_24hr_MOI_10_fc_1_up_v1)),c(1:5)]
#--
#Save down
#--
write.table(final_trended_up_genes_human_1hr_annotated_v1, "final_trended_up_genes_human_1hr_annotated.csv", row.names = F, col.names = T, sep = ",")
write.table(final_trended_up_genes_human_24hr_annotated_v1, "final_trended_up_genes_human_24hr_annotated.csv", row.names = F, col.names = T, sep = ",")
write.table(final_trended_down_genes_human_24hr_annotated_v1, "final_trended_down_genes_human_24hr_annotated.csv", row.names = F, col.names = T, sep = ",")


#--
#CT
#--
#1hr
dim(CT_1hr_trend_genes_dw_v1)[1] == length(intersect(row.names(CT_de_1hr_MOI_01_fc_1_dw_v1),row.names(CT_de_1hr_MOI_10_fc_1_up_v1))) 
final_trended_down_genes_CT_1hr_annotated_v1 <- DE_genes_CT_1hr_1_vs_01_v5[DE_genes_CT_1hr_1_vs_01_v5$Gene_name %in% 
                                                                             intersect(row.names(CT_de_1hr_MOI_01_fc_1_dw_v1),row.names(CT_de_1hr_MOI_10_fc_1_up_v1)),c(1:5)]

#24hr
dim(CT_24hr_trend_genes_dw_v1)[1] == length(intersect(row.names(CT_de_24hr_MOI_01_fc_1_dw_v1),row.names(CT_de_24hr_MOI_10_fc_1_up_v1))) 
final_trended_down_genes_CT_24hr_annotated_v1 <- DE_genes_CT_24hr_1_vs_01_v5[DE_genes_CT_24hr_1_vs_01_v5$Gene_name %in% 
                                                                               intersect(row.names(CT_de_24hr_MOI_01_fc_1_dw_v1),row.names(CT_de_24hr_MOI_10_fc_1_up_v1)),]
#--
#Save down
#--
write.table(final_trended_down_genes_CT_1hr_annotated_v1, "final_trended_down_genes_CT_1hr_annotated.csv", row.names = F, col.names = T, sep = ",")
write.table(final_trended_down_genes_CT_24hr_annotated_v1, "final_trended_down_genes_CT_24hr_annotated.csv", row.names = F, col.names = T, sep = ",")




#**************************************
#
# Make heatmaps of the trended genes ----
#
#**************************************

#Create a dataframe of the normalised expression values to use for the heatmap
heatmap_expression_1hr_df_v1 <- as.data.frame(human_rle_plots_norm_1hr)
heatmap_expression_24hr_df_v1 <- as.data.frame(human_rle_plots_norm_24hr)
heatmap_expression_1hr_df_v1$gene <- row.names(heatmap_expression_1hr_df_v1)
heatmap_expression_24hr_df_v1$gene <- row.names(heatmap_expression_24hr_df_v1)

#Filter the values to just the subsets of interest
heatmap_expression_1hr_trended_up_df_v2 <- heatmap_expression_1hr_df_v1[heatmap_expression_1hr_df_v1$gene %in% final_trended_up_genes_human_1hr_annotated_v1$GENEID,]
heatmap_expression_24hr_trended_up_df_v2 <- heatmap_expression_24hr_df_v1[heatmap_expression_24hr_df_v1$gene %in% final_trended_up_genes_human_24hr_annotated_v1$GENEID,]
heatmap_expression_24hr_trended_down_df_v2 <- heatmap_expression_24hr_df_v1[heatmap_expression_24hr_df_v1$gene %in% final_trended_down_genes_human_24hr_annotated_v1$GENEID,]

#Combine up and down at 24 hours as there's only one at 24-down 
heatmap_expression_24hr_trended_both <- rbind(heatmap_expression_24hr_trended_up_df_v2, heatmap_expression_24hr_trended_down_df_v2)


#add in the readable gene names
heatmap_expression_1hr_trended_up_df_v3 <- merge(x=heatmap_expression_1hr_trended_up_df_v2, y=final_trended_up_genes_human_1hr_annotated_v1[,1:2], by.x="gene", by.y="GENEID", all.x=T)
row.names(heatmap_expression_1hr_trended_up_df_v3) <- heatmap_expression_1hr_trended_up_df_v3$GENENAME.x

heatmap_expression_24hr_trended_both_v2 <- merge(x=heatmap_expression_24hr_trended_both, y=final_trended_up_genes_human_24hr_annotated_v1[,1:2], by.x="gene", by.y="GENEID", all.x=T)
#Add in the one entry for the down-trend
heatmap_expression_24hr_trended_both_v2[14,19] <- "TXNIP"
row.names(heatmap_expression_24hr_trended_both_v2) <- heatmap_expression_24hr_trended_both_v2$GENENAME.x

#--
#Create the DF's to use for the plots
#--
#1 hr up
heatmap_annotation_trended_up_human <- data.frame(row.names = colnames(heatmap_expression_1hr_trended_up_df_v3)[2:18],
                                                  "MOI" = rep(c("0.1","1","10"),each=6)[1:17],
                                                  "Depletion" = rep(rep(c("rRNA","rRNA+polyA"), each=3),3)[1:17])


#Column annotation - 24 hr up
heatmap_annotation_trended_up_24_human <- data.frame(row.names = colnames(heatmap_expression_24hr_trended_up_df_v2)[1:17],
                                                     "MOI" = rep(c("0.1","1","10"),each=6)[c(1:17)],
                                                     "Depletion" = rep(rep(c("rRNA","rRNA+polyA"), each=3),3)[c(1:17)])


#--
#Create heatmap using pheatmap
#--

#Set plot params
width=6;height=6

#Human
#1hr - up
Cairo(file="Heatmap_Trended_up_human_1hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
pheatmap(heatmap_expression_1hr_trended_up_df_v3[,2:18], scale="row", cluster_cols = F, show_colnames = F,
         annotation_col = heatmap_annotation_trended_up_human, gaps_col = c(6,12), fontsize_row = 7, main = "Trended Genes - human 1hr - up")
dev.off()

#24hr - up
Cairo(file="Heatmap_Trended_up_human_24hr.png", type="png", units="in", width=width, height=height, dpi=dpi, bg = "white")
pheatmap(heatmap_expression_24hr_trended_both_v2[,2:18], scale="row", cluster_cols = F, cluster_rows = T, show_colnames = F,
         annotation_col = heatmap_annotation_trended_up_24_human, gaps_col = c(6,12), main = "Trended Genes - human 24hr - both ")
dev.off()





#**************************************
#
# Then of the top DE genes, enrich pathways ----
#
#**************************************

#--
#Show to terminal the numbers of up and down reg-genes
#--

#human
#1hr
dim(human_1hr_moi_1_vs_01_pathways_v7[human_1hr_moi_1_vs_01_pathways_v7$logFC > 0,]) 
dim(human_1hr_moi_1_vs_01_pathways_v7[human_1hr_moi_1_vs_01_pathways_v7$logFC < 0,])
dim(human_1hr_moi_1_vs_10_pathways_v7[human_1hr_moi_1_vs_10_pathways_v7$logFC > 0,])
dim(human_1hr_moi_1_vs_10_pathways_v7[human_1hr_moi_1_vs_10_pathways_v7$logFC < 0,])
#24hr
dim(human_24hr_moi_1_vs_01_pathways_v7[human_24hr_moi_1_vs_01_pathways_v7$logFC > 0,])
dim(human_24hr_moi_1_vs_01_pathways_v7[human_24hr_moi_1_vs_01_pathways_v7$logFC < 0,])
dim(human_24hr_moi_1_vs_10_pathways_v7[human_24hr_moi_1_vs_10_pathways_v7$logFC > 0,])
dim(human_24hr_moi_1_vs_10_pathways_v7[human_24hr_moi_1_vs_10_pathways_v7$logFC < 0,])


#--
#Separate out into up and down fold changes
#--
#1hr
human_1hr_moi_1_vs_01_enriched_pathways_up_v1 <- human_1hr_moi_1_vs_01_pathways_v7[human_1hr_moi_1_vs_01_pathways_v7$logFC > 0,]
human_1hr_moi_1_vs_01_enriched_pathways_down_v1 <- human_1hr_moi_1_vs_01_pathways_v7[human_1hr_moi_1_vs_01_pathways_v7$logFC < 0,]
human_1hr_moi_1_vs_10_enriched_pathways_up_v1 <- human_1hr_moi_1_vs_10_pathways_v7[human_1hr_moi_1_vs_10_pathways_v7$logFC > 0,]
human_1hr_moi_1_vs_10_enriched_pathways_down_v1 <- human_1hr_moi_1_vs_10_pathways_v7[human_1hr_moi_1_vs_10_pathways_v7$logFC < 0,]
#24hr
human_24hr_moi_1_vs_01_enriched_pathways_up_v1 <- human_24hr_moi_1_vs_01_pathways_v7[human_24hr_moi_1_vs_01_pathways_v7$logFC > 0,]
human_24hr_moi_1_vs_01_enriched_pathways_down_v1 <- human_24hr_moi_1_vs_01_pathways_v7[human_24hr_moi_1_vs_01_pathways_v7$logFC < 0,]
human_24hr_moi_1_vs_10_enriched_pathways_up_v1 <- human_24hr_moi_1_vs_10_pathways_v7[human_24hr_moi_1_vs_10_pathways_v7$logFC > 0,]
human_24hr_moi_1_vs_10_enriched_pathways_down_v1 <- human_24hr_moi_1_vs_10_pathways_v7[human_24hr_moi_1_vs_10_pathways_v7$logFC < 0,]


#--
#Use enrichr to annotate the genes
#--

#set and see the databases to choose from
dbs <- listEnrichrDbs()
head(dbs)

#Set the databases
dbs2 <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018","KEGG_2016")

#--
#Enrich
#--
#1hr
human_1hr_moi_1_vs_01_enriched_pathways_up_v2 <- enrichr(as.vector(human_1hr_moi_1_vs_01_enriched_pathways_up_v1$GENENAME.x),dbs2)
human_1hr_moi_1_vs_01_enriched_pathways_down_v2 <- enrichr(as.vector(human_1hr_moi_1_vs_01_enriched_pathways_down_v1$GENENAME.x),dbs2)
human_1hr_moi_1_vs_10_enriched_pathways_up_v2 <- enrichr(as.vector(human_1hr_moi_1_vs_10_enriched_pathways_up_v1$GENENAME.x),dbs2)
human_1hr_moi_1_vs_10_enriched_pathways_down_v2 <- enrichr(as.vector(human_1hr_moi_1_vs_10_enriched_pathways_down_v1$GENENAME.x),dbs2)
#24hr
human_24hr_moi_1_vs_01_enriched_pathways_up_v2 <- enrichr(as.vector(human_24hr_moi_1_vs_01_enriched_pathways_up_v1$GENENAME.x),dbs2)
human_24hr_moi_1_vs_01_enriched_pathways_down_v2 <- enrichr(as.vector(human_24hr_moi_1_vs_01_enriched_pathways_down_v1$GENENAME.x),dbs2)
human_24hr_moi_1_vs_10_enriched_pathways_up_v2 <- enrichr(as.vector(human_24hr_moi_1_vs_10_enriched_pathways_up_v1$GENENAME.x),dbs2)
human_24hr_moi_1_vs_10_enriched_pathways_down_v2 <- enrichr(as.vector(human_24hr_moi_1_vs_10_enriched_pathways_down_v1$GENENAME.x),dbs2)


#--
#set up the pathway plot function
#--
simple_plot <- function(filename, the_data, plot_colours, wrap_text, text_size, plot_title){
  Cairo(file=filename, type="png", units="in", width=6, height=8, dpi=300, bg = "white")
  print(ggplot(data=the_data, aes(y=Combined.Score ,x=reorder(Term, Combined.Score))) +
          geom_bar(stat="identity", fill=plot_colours) + coord_flip() +
          ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5)) + scale_x_discrete(name = "", labels = wrap_format(wrap_text)) +
          theme(axis.text.y = element_text(size=text_size)))
  dev.off()
}
wrap_text_size = 25
plot_text_size = 14


#Set colours
col_moi_01 <- "gray47"
col_moi_1 <- "tomato"
col_moi_10 <- "lightgreen"

#Select the top pathways
top_pathways <- 1:10
top_pathways <- 1:5

#1hr
human_1hr_moi_1_vs_01_enriched_pathways_up_v3 <- human_1hr_moi_1_vs_01_enriched_pathways_up_v2$KEGG_2016[top_pathways,]
human_1hr_moi_1_vs_01_enriched_pathways_down_v3 <- human_1hr_moi_1_vs_01_enriched_pathways_down_v2$KEGG_2016[top_pathways,]
human_1hr_moi_1_vs_10_enriched_pathways_up_v3 <- human_1hr_moi_1_vs_10_enriched_pathways_up_v2$KEGG_2016[top_pathways,]
human_1hr_moi_1_vs_10_enriched_pathways_down_v3 <- human_1hr_moi_1_vs_10_enriched_pathways_down_v2$KEGG_2016[top_pathways,]
#24hrs
human_24hr_moi_1_vs_01_enriched_pathways_up_v3 <- human_24hr_moi_1_vs_01_enriched_pathways_up_v2$KEGG_2016[top_pathways,]
human_24hr_moi_1_vs_01_enriched_pathways_down_v3 <- human_24hr_moi_1_vs_01_enriched_pathways_down_v2$KEGG_2016[top_pathways,]
human_24hr_moi_1_vs_10_enriched_pathways_up_v3 <- human_24hr_moi_1_vs_10_enriched_pathways_up_v2$KEGG_2016[top_pathways,]
human_24hr_moi_1_vs_10_enriched_pathways_down_v3 <- human_24hr_moi_1_vs_10_enriched_pathways_down_v2$KEGG_2016[top_pathways,]


#--
#Tidy up the names for KEGG pathways
#--
#1hr
human_1hr_moi_1_vs_01_enriched_pathways_up_v3$Term <- substr(human_1hr_moi_1_vs_01_enriched_pathways_up_v3$Term,1,nchar(human_1hr_moi_1_vs_01_enriched_pathways_up_v3$Term)-22)
human_1hr_moi_1_vs_01_enriched_pathways_down_v3$Term <- substr(human_1hr_moi_1_vs_01_enriched_pathways_down_v3$Term,1,nchar(human_1hr_moi_1_vs_01_enriched_pathways_down_v3$Term)-22)
human_1hr_moi_1_vs_10_enriched_pathways_up_v3$Term <- substr(human_1hr_moi_1_vs_10_enriched_pathways_up_v3$Term,1,nchar(human_1hr_moi_1_vs_10_enriched_pathways_up_v3$Term)-22)
human_1hr_moi_1_vs_10_enriched_pathways_down_v3$Term <- substr(human_1hr_moi_1_vs_10_enriched_pathways_down_v3$Term,1,nchar(human_1hr_moi_1_vs_10_enriched_pathways_down_v3$Term)-22)
#24hr
human_24hr_moi_1_vs_01_enriched_pathways_up_v3$Term <- substr(human_24hr_moi_1_vs_01_enriched_pathways_up_v3$Term,1,nchar(human_24hr_moi_1_vs_01_enriched_pathways_up_v3$Term)-22)
human_24hr_moi_1_vs_01_enriched_pathways_down_v3$Term <- substr(human_24hr_moi_1_vs_01_enriched_pathways_down_v3$Term,1,nchar(human_24hr_moi_1_vs_01_enriched_pathways_down_v3$Term)-22)
human_24hr_moi_1_vs_10_enriched_pathways_up_v3$Term <- substr(human_24hr_moi_1_vs_10_enriched_pathways_up_v3$Term,1,nchar(human_24hr_moi_1_vs_10_enriched_pathways_up_v3$Term)-22)
human_24hr_moi_1_vs_10_enriched_pathways_down_v3$Term <- substr(human_24hr_moi_1_vs_10_enriched_pathways_down_v3$Term,1,nchar(human_24hr_moi_1_vs_10_enriched_pathways_down_v3$Term)-22)



#--
#Plot the results
#--

col_up <- "firebrick1"
col_dw <- "royalblue"

#Human
#1hr
simple_plot(filename = "pathway_KEGG_human_1hr_moi_1_vs_01_up.png", the_data = human_1hr_moi_1_vs_01_enriched_pathways_up_v3, 
            plot_colours = col_up, wrap_text = wrap_text_size, text_size = plot_text_size, plot_title = "1hr - MOI 1 vs 0.1 - up")

simple_plot(filename = "pathway_KEGG_human_1hr_moi_1_vs_01_down.png", the_data = human_1hr_moi_1_vs_01_enriched_pathways_down_v3, 
            plot_colours = col_dw, wrap_text = wrap_text_size, text_size = plot_text_size, plot_title = "1hr - MOI 1 vs 0.1 - down")

simple_plot(filename = "pathway_KEGG_human_1hr_moi_1_vs_10_up.png", the_data = human_1hr_moi_1_vs_10_enriched_pathways_up_v3, 
            plot_colours = col_up, wrap_text = wrap_text_size, text_size = plot_text_size, plot_title = "1hr - MOI 1 vs 10 - up")

simple_plot(filename = "pathway_KEGG_human_1hr_moi_1_vs_10_down.png", the_data = human_1hr_moi_1_vs_10_enriched_pathways_down_v3, 
            plot_colours = col_dw, wrap_text = wrap_text_size, text_size = plot_text_size, plot_title = "1hr - MOI 1 vs 10 - down")

#24hr
simple_plot(filename = "pathway_KEGG_human_24hr_moi_1_vs_01_up.png", the_data = human_24hr_moi_1_vs_01_enriched_pathways_up_v3, 
            plot_colours = col_up, wrap_text = wrap_text_size, text_size = plot_text_size, plot_title = "24hr - MOI 1 vs 0.1 - up")

simple_plot(filename = "pathway_KEGG_human_24hr_moi_1_vs_01_down.png", the_data = human_24hr_moi_1_vs_01_enriched_pathways_down_v3, 
            plot_colours = col_dw, wrap_text = wrap_text_size, text_size = plot_text_size, plot_title = "24hr - MOI 1 vs 0.1 - down")

simple_plot(filename = "pathway_KEGG_human_24hr_moi_1_vs_10_up.png", the_data = human_24hr_moi_1_vs_10_enriched_pathways_up_v3, 
            plot_colours = col_up, wrap_text = wrap_text_size, text_size = plot_text_size, plot_title = "24hr - MOI 1 vs 10 - up")

simple_plot(filename = "pathway_KEGG_human_24hr_moi_1_vs_10_down.png", the_data = human_24hr_moi_1_vs_10_enriched_pathways_down_v3, 
            plot_colours = col_dw, wrap_text = wrap_text_size, text_size = plot_text_size, plot_title = "24hr - MOI 1 vs 10 - down")






#**************************************
#
# Create host pathway graphs that increase enrichment over MOI's ----
#
#**************************************
top_pathways_overlap <- 1:10

human_1hr_moi_1_vs_01_enriched_pathways_up_overlap_v1 <- human_1hr_moi_1_vs_01_enriched_pathways_up_v2$KEGG_2016[top_pathways_overlap,]
human_1hr_moi_1_vs_01_enriched_pathways_down_overlap_v1 <- human_1hr_moi_1_vs_01_enriched_pathways_down_v2$KEGG_2016[top_pathways_overlap,]

human_1hr_moi_1_vs_10_enriched_pathways_up_overlap_v1 <- human_1hr_moi_1_vs_10_enriched_pathways_up_v2$KEGG_2016[top_pathways_overlap,]
human_1hr_moi_1_vs_10_enriched_pathways_down_overlap_v1 <- human_1hr_moi_1_vs_10_enriched_pathways_down_v2$KEGG_2016[top_pathways_overlap,]



#24hrs
human_24hr_moi_1_vs_01_enriched_pathways_up_overlap_v1 <- human_24hr_moi_1_vs_01_enriched_pathways_up_v2$KEGG_2016[top_pathways_overlap,]
human_24hr_moi_1_vs_01_enriched_pathways_down_overlap_v1 <- human_24hr_moi_1_vs_01_enriched_pathways_down_v2$KEGG_2016[top_pathways_overlap,] 

human_24hr_moi_1_vs_10_enriched_pathways_up_overlap_v1 <- human_24hr_moi_1_vs_10_enriched_pathways_up_v2$KEGG_2016[top_pathways_overlap,]
human_24hr_moi_1_vs_10_enriched_pathways_down_overlap_v1 <- human_24hr_moi_1_vs_10_enriched_pathways_down_v2$KEGG_2016[top_pathways_overlap,]


#Look for overlaps that increase over time
#1hr
intersect(human_1hr_moi_1_vs_01_enriched_pathways_up_overlap_v1$Term,human_1hr_moi_1_vs_10_enriched_pathways_down_overlap_v1$Term)
#TNF signaling                                          115.62673 76.81002
#NF-kappa B signaling                                   62.58421  27.16542
#Cytokine-cytokine receptor interaction                 49.93421  30.48000
#NOD-like receptor signaling                            49.44465  33.07911
#AGE-RAGE signaling pathway in diabetic complications   49.44465  28.12171
#Rheumatoid arthritis                                   41.00843  24.51826
#Osteoclast differentiation                             36.25215  36.25215

#24hr
intersect(human_24hr_moi_1_vs_01_enriched_pathways_up_overlap_v1$Term,human_24hr_moi_1_vs_10_enriched_pathways_down_overlap_v1$Term)
#TNF signaling          47.40138 38.72732
#NF-kappa B signaling   36.75498 31.16056



#Look for overlaps that decrease over time
#1hr
intersect(human_1hr_moi_1_vs_01_enriched_pathways_down_overlap_v1$Term,human_1hr_moi_1_vs_10_enriched_pathways_up_overlap_v1$Term) #cant use as no p.value.adj is significant
#24hr
intersect(human_24hr_moi_1_vs_01_enriched_pathways_down_overlap_v1$Term,human_24hr_moi_1_vs_10_enriched_pathways_up_overlap_v1$Term) #Can use - will probably create a cutuff for the combined score though
#Carbon metabolism          8.301496 69.65291
#Citrate cycle (TCA cycle)  7.878103 36.32700


#Make the dataframes

#1hr increasing
increasing_1 <- data.frame("Pathway" = c("TNF signaling", "NF-kappa B signaling","Cytokine-cytokine receptor interaction","NOD-like receptor signaling",
                                         "AGE-RAGE signaling pathway in diabetic complications","Rheumatoid arthritis","Osteoclast differentiation"),
                           "Combined_score_1" = c(115.62673,62.58421,49.93421,49.44465,49.44465,41.00843,36.25215),
                           "Combined_score_2" = c(76.81002,27.16542,30.48,33.07911,28.12171,24.51826,36.25215))
increasing_1$Combined.Score = increasing_1$Combined_score_1 + increasing_1$Combined_score_2


#24hr increasing
increasing_24 <- data.frame("Pathway" = c("TNF signaling","NF-kappa B signaling"),
                            "Combined_score_1" = c(47.40138,36.75498),
                            "Combined_score_2" = c(38.72732,31.16056))
increasing_24$Combined.Score = increasing_24$Combined_score_1 + increasing_24$Combined_score_2


#24hr decreasing
decreasing_24 <- data.frame("Pathway" = c("Carbon metabolism","Citrate cycle (TCA cycle)"),
                            "Combined_score_1" = c(8.301496,7.878103),
                            "Combined_score_2" = c(69.65291,36.32700))
decreasing_24$Combined.Score = decreasing_24$Combined_score_1 + decreasing_24$Combined_score_2


grDevices::adjustcolor("firebrick1", alpha.f = 0.7)
grDevices::adjustcolor("royalblue", alpha.f = 0.7)

#1hr up
Cairo(file="enrichment_1hr_up.png", type="png", units="in", width=3, height=3, dpi=dpi, bg = "white")
ggplot(data=increasing_1, aes(y=Combined.Score ,x=reorder(Pathway, Combined.Score))) +
          geom_bar(stat="identity", fill="#FF3030B3") + coord_flip() +
          ggtitle("plot_title") + theme(plot.title = element_text(hjust = 0.5)) + scale_x_discrete(name = "", labels = wrap_format(25)) +
          theme(axis.text.y = element_text(size=6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey"))
dev.off()


#24hr up
Cairo(file="enrichment_24hr_up.png", type="png", units="in", width=3, height=2, dpi=dpi, bg = "white")
ggplot(data=increasing_24, aes(y=Combined.Score ,x=reorder(Pathway, Combined.Score))) +
  geom_bar(stat="identity", fill="#FF3030B3") + coord_flip() +
  ggtitle("plot_title") + theme(plot.title = element_text(hjust = 0.5)) + scale_x_discrete(name = "", labels = wrap_format(25)) +
  theme(axis.text.y = element_text(size=8)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey"))
dev.off()


#24hr down
Cairo(file="enrichment_24hr_down.png", type="png", units="in", width=3, height=2, dpi=dpi, bg = "white")
ggplot(data=decreasing_24, aes(y=Combined.Score ,x=reorder(Pathway, Combined.Score))) +
  geom_bar(stat="identity", fill="#4169E1B3") + coord_flip() +
  ggtitle("plot_title") + theme(plot.title = element_text(hjust = 0.5)) + scale_x_discrete(name = "", labels = wrap_format(25)) +
  theme(axis.text.y = element_text(size=8)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey"))
dev.off()