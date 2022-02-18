library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(data.table)
library(matrixStats)
library(gplots)
library(matrixTests)

#using filtered data from the first analysis
filename = "Noise_filtered_data.csv"
analysis.1.data <- read.csv(filename)

#converting probeids to rownames from a column
analysis.1.data.1 <- analysis.1.data[,-1]
rownames(analysis.1.data.1) <- analysis.1.data[,1]
analysis.1.data <- analysis.1.data.1
rm(analysis.1.data.1)

#transposing to have samples as rownames and probeids as column names
analysis.1.data_tp <- transpose(analysis.1.data)

# get row and col names in order
colnames(analysis.1.data_tp) <- rownames(analysis.1.data)
rownames(analysis.1.data_tp) <- colnames(analysis.1.data)

# making a distance matrix for clustering
dist_mat <- dist(analysis.1.data_tp, method = 'euclidean')

#Clustering performed
hclust_avg <- hclust(dist_mat)

#Cutting into two clusters
cut_avg <- cutree(hclust_avg, k = 2)
cut_vec <- as.vector(cut_avg)

#number of samples in each cluster
group_1 <- sum(cut_vec == 1)
group_2 <- sum(cut_vec == 2)

# HEATMAP

Y<-data.matrix(analysis.1.data_tp)
# rownames(Y)<analysis.1.data_tp[[1]]

#importing metadata
metadata<-read.csv("/projectnb/bf528/users/im_not_dead_yet/project_1/proj_metadata.csv")

#Applying conditon to cluster C3 as red and all else as blue
condition_colors <-
  transmute(
    metadata,
    color=if_else(cit.coloncancermolecularsubtype == "C3","red","blue")
  ) #coloring by subtypes


#plotting heatmap
png('heatmap.png', width=1920, height=1080, res=100)

heatmap.2(t(Y),ColSideColors=condition_colors$color, xlab="Patient Samples", ylab="Genes",
          main="Gene Expression Across Samples",trace="none", density.info = "none",
          key.xlab="Expression Level", scale="row", cexRow = 0.5, cexCol = 0.2, key = TRUE)
legend(x=0.9, y=1, xpd=TRUE, inset=c(-0.15,0),
       legend=c('C3', 'Other'), title='Molecular Cancer Subtype', fill=c('red', 'blue'))
dev.off()

#WELCH'S T-TEST

#Adding cluster assignment along the probeids
grouped_data <- analysis.1.data_tp
grouped_data$Group <- cut_vec

#dividing samples between two clusters
cluster_1 <- t(data.matrix(grouped_data[grouped_data$Group == 1,]))
cluster_2 <- t(data.matrix(grouped_data[grouped_data$Group == 2,]))

#Performing Welch's t-test between cluster 1 and cluster 2
welch.t.test <- row_t_welch(cluster_1, cluster_2)
p.values <- welch.t.test$pvalue
t_stat <- welch.t.test$statistic
p_adjusted = p.adjust(p.values,method = "fdr")
probeids <- rownames(cluster_1)

welch_test_df <- data.frame(probeids,t_stat,p.values, p_adjusted)
rownames(welch_test_df)<- welch_test_df$probeids
#Genes with maximum differential expression
max_diff_expr <- (head(rownames(welch_test_df[rev(order(abs(welch_test_df$t_stat))),])))

#The number of differentially genes
diff.expr.genes <- filter(welch_test_df, p_adjusted < 0.05)



#DELIVERABLES 

#report the number of samples in each cluster from Part 5.2
print(paste0('The number of samples in group 1 is ', group_1))
print(paste0('The number of samples in group 2 is ', group_2))

# a heatmap of the genes and samples with a color bar indicating which subtype each sample belongs to
#Refer to heatmap.png

# report the number of differentially expressed genes at p<0.05 between the two clusters
print(paste0('The number of differentially expressed genes at p<0.05 between the two clusters are ', nrow(diff.expr.genes)))

# a comma separated file containing the results of the Welch t-test for all genes irrespective of significance for each subtype comparison
write.csv(welch_test_df, "Clustered_and_subtype_disc_data.csv")

# report a list of the genes you feel best represent each cluster and explain how you came to your conclusion
print((head(rownames(welch_test_df[rev(order(abs(welch_test_df$t_stat))),]))))


