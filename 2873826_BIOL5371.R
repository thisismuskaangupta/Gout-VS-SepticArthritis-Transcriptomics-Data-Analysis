#This is the R script for my report.
#First we read our data into our environment.
annotations = read.table("C:/Users/2873826G/Documents/My Masters Degree UofG/Statistics/report/Annotations.csv",header=TRUE,row.names=1,sep = "\t")
de_gout_vs_hc = read.table("C:/Users/2873826G/Documents/My Masters Degree UofG/Statistics/report/DE_Gout_vs_HC.csv",header=TRUE,row.names=1,sep = "\t")
de_sa_vs_hc = read.table("C:/Users/2873826G/Documents/My Masters Degree UofG/Statistics/report/DE_SA_vs_HC.csv",header=TRUE,row.names=1,sep = "\t")
expression_table = read.table("C:/Users/2873826G/Documents/My Masters Degree UofG/Statistics/report/Expression_Table.csv",header=TRUE,row.names=1,sep = "\t")
sample_info = read.table("C:/Users/2873826G/Documents/My Masters Degree UofG/Statistics/report/Sample_Information.csv",header=TRUE,row.names=1,sep = "\t")
#Everything looks good so far.

#checking structure of sample_info so we can do background checks on the data.
str(sample_info)
#making the group and sex factors.
sample_info$GROUP = as.factor(sample_info$GROUP)
sample_info$SEX=as.factor(sample_info$SEX)
str(sample_info)
#now it looks good, and we can start our analysis.

#We check the DE tables for both comparisons - 
#questions - 
#are our groups well-matched clinically?
#checking the summary statistics.
summary(de_gout_vs_hc)
summary(de_sa_vs_hc)
summary(expression_table)
summary(sample_info)
#There are 4 females in total, 2 from HC and 2 from SA, none from gout.
#median age is 58.5, which is skewed.
#we can visualize the spread of our data using histograms.
ggp1 = ggplot(sample_info,aes(x=AGE)) + geom_histogram()
  +labs(title='Histogram of Age Across All Samples')
ggp1
#left skewed.
#we've checked the demographics; for neutrophils and monocytes, it does not make sense to check the summary statistics or plot a histogram because they are white blood cells and are affected by disease status, so the values are not independent like sex and age are. 

#separating the sample_info table into three components.
sample_info_hc = subset(sample_info,GROUP=='HEALTHY')
sample_info_sa = subset(sample_info,GROUP=='SA')
sample_info_gout = subset(sample_info,GROUP=='GOUT')

#checking the summary statistics of each group.
summary(sample_info_hc)
summary(sample_info_gout)
summary(sample_info_sa)
#write the 3 medians in the report, for all three tables.
#for all three, age is skewed towards older age, and for gout and SA, WBCs have higher medians compared to healthy, and roughly the same compared to each other.

#making histograms of the clinical measurements (by each clinical measurement).
ggp2 = ggplot(sample_info_hc,aes(x=AGE))+geom_histogram()+labs(title='Histogram of Age in Healthy Controls')
ggp2
ggp3 = ggplot(sample_info_sa,aes(x=AGE))+geom_histogram()+labs(title='Histogram of Age in SA Patients')
ggp3
ggp4 = ggplot(sample_info_gout,aes(x=AGE))+geom_histogram()+labs(title='Histogram of Age in Gout Patients')
ggp4

ggp5 = ggplot(sample_info_hc,aes(x=NEUTROPHILS))+geom_histogram()+labs(title='Histogram of Neutrophils in Healthy Controls')
ggp5
ggp6 = ggplot(sample_info_sa,aes(x=NEUTROPHILS))+geom_histogram()+labs(title='Histogram of Neutrophils in SA Patients')
ggp6
ggp7 = ggplot(sample_info_gout,aes(x=NEUTROPHILS))+geom_histogram()+labs(title='Histogram of Neutrophils in Gout Patients')
ggp7

ggp8 = ggplot(sample_info_hc,aes(x=MONOCYTES))+geom_histogram()+labs(title='Histogram of Monocytes in Healthy Controls')
ggp8
ggp9 = ggplot(sample_info_sa,aes(x=MONOCYTES))+geom_histogram()+labs(title='Histogram of Monocytes in SA Patients')
ggp9
ggp10 = ggplot(sample_info_gout,aes(x=MONOCYTES))+geom_histogram()+labs(title='Histogram of Monocytes in Gout Patients')
ggp10

#more background - on gene expression, we can use PCA plotting.
#first of all we can make a PCA plot of gene expression across our three groups.
PCA = prcomp(t(expression_table))
pca_coordinates = data.frame(PCA$x)
ggp11 = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = sample_info$GROUP)) +
  geom_point() + labs(title='PCA Plot of Gene Expression Across the Three Groups')
ggp11
#There is a no noticeable clustering.

#we can do the same for age.
PCA = prcomp(t(expression_table))
pca_coordinates = data.frame(PCA$x)
ggp12 = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = sample_info$AGE)) +
  geom_point() + labs(title='PCA Plot of Gene Expression Across Age (All Replicates)')
ggp12

#what variance is explained by each principal component?
library(factoextra)
ggp15 = fviz_eig(PCA)
ggp15
#refer to the lecture on PCA logic (lecture 8) to discuss this.

#annotating the DE files and expression table using the annotations file.
expression_table_annotated = merge(expression_table, annotations,by.x=0,by.y=0)
de_gout_vs_hc_annotated = merge(de_gout_vs_hc,annotations,by.x=0,by.y=0)
de_sa_vs_hc_annotated = merge(de_sa_vs_hc,annotations,by.x=0,by.y=0)

#making the ensembl IDs the row names for the annotated files.
row.names(de_gout_vs_hc_annotated) = de_gout_vs_hc_annotated[,'Row.names']
row.names(de_sa_vs_hc_annotated) = de_sa_vs_hc_annotated[,'Row.names']
row.names(expression_table_annotated) = expression_table_annotated[,'Row.names']

#renaming the first column in the annotated files to 'ensembl_ID'.
names(expression_table_annotated)[1]='ensembl_ID'
names(de_gout_vs_hc_annotated)[1]='ensembl_ID'
names(de_sa_vs_hc_annotated)[1]='ensembl_ID'

#getting significant genes for both gout vs hc DE and sa vs hc DE.
de_gout_vs_hc_annotated_sig = subset(de_gout_vs_hc_annotated, p.adj < 0.05)
de_sa_vs_hc_annotated_sig = subset(de_sa_vs_hc_annotated, p.adj < 0.05)

#how many significant genes are there?
nrow(de_gout_vs_hc_annotated_sig)
nrow(de_sa_vs_hc_annotated_sig)
#the answer - for gout vs hc, it is 5552.
#and for sa vs hc, it is 6956.

#from the tables we can read the top most significant gene for both comparisons. 
#I have noticed that the genes METTL7B and HP are both within the top 5 most significantly differentially expressed in both our DE comparisons. This means that both these genes' expression is affected in the cases of both SA and gout.

#making histograms of fold change values for both gout vs hc and sa vs hc DE.
#ggp13 = ggplot(de_gout_vs_hc_annotated_sig,aes(x=log2fold))+geom_histogram() + labs(title='Gout VS HC Fold Change Histogram')
#ggp14 = ggplot(de_sa_vs_hc_annotated_sig,aes(x=log2fold))+geom_histogram() + labs(title='SA VS HC Fold Change Histogram')
#for the sake of conciseness, I'm leaving this out of the report.
#in histograms of both datasets using the log2fold change value on x, we see bimodal histograms. this indicates that in both cases, upregulation as well as downregulation can be observed.
#we can depict these using volcano plots as well.
ggp16 = ggplot(de_gout_vs_hc_annotated_sig,aes(x=log2fold,y=-log10(p.adj)))+geom_point() + labs(title='Gout VS HC Volcano Plot')
ggp17 = ggplot(de_sa_vs_hc_annotated_sig,aes(x=log2fold,y=-log10(p.adj)))+geom_point() + labs(title='SA VS HC Volcano Plot')
ggp18 = ggplot(de_gout_vs_hc_annotated_sig,aes(x=log2fold,y=-log10(p.adj)))+geom_point()+geom_point(data=de_sa_vs_hc_annotated_sig,colour='red')+labs(title='Volcano Plot - SA VS HC in Red and Gout VS HC in Black')

#are the genes similar? 
#we can see using our volcano plot that most of the genes match and coincide with one another..
#we can test this in a numerical way by seeing which genes are common between the two comparisons (significant ones).
de_gout_vs_hc_sig_geneIDs = de_gout_vs_hc_annotated_sig$ensembl_ID
de_sa_vs_hc_sig_geneIDs = de_sa_vs_hc_annotated_sig$ensembl_ID
common_sig_geneIDs = intersect(de_gout_vs_hc_sig_geneIDs,de_sa_vs_hc_sig_geneIDs)
length(common_sig_geneIDs)
#By running this, we find that 3661 genes are common in both the significant differential analyses.

#we can create tables for gout vs hc and sa vs hc top 5000 genes.
de_gout_vs_hc_annotated_sig_top5000 = de_gout_vs_hc_annotated_sig[(order(de_gout_vs_hc_annotated_sig$p.adj,decreasing=FALSE)),]
de_gout_vs_hc_annotated_sig_top5000 = head(de_gout_vs_hc_annotated_sig_top5000,5000)

de_sa_vs_hc_annotated_sig_top5000 = de_sa_vs_hc_annotated_sig[(order(de_sa_vs_hc_annotated_sig$p.adj,decreasing=FALSE)),]
de_sa_vs_hc_annotated_sig_top5000 = head(de_sa_vs_hc_annotated_sig_top5000,5000)

#Making a heatmap comparing top 5000 genes.
#ggp6 = ggplot(de_gout_vs_hc_annotated_sig_top5000,aes(x=de_gout_vs_hc_annotated_sig_top5000$symbol,y=de_sa_vs_hc_annotated_sig_top5000$symbol,fill=log2fold))+geom_tile() + labs(x='Gout VS HC Symbols',y='SA VS HC Symbols',title = 'Top 5000 Most Significant Differentially Expressed Genes')

#Making a heatmap comparing top 10 genes.
#ggp7 = ggplot(de_gout_vs_hc_annotated_sig_top10,aes(x=de_gout_vs_hc_annotated_sig_top10$symbol,y=de_sa_vs_hc_annotated_sig_top10$symbol,fill=log2fold))+geom_tile() + labs(x='Gout VS HC Symbols',y='SA VS HC Symbols',title = 'Top 10 Most Significant Differentially Expressed Genes')

#I'm making gene data tables for METTL7B and HP.
METTL7B_gene_data = data.frame(t(expression_table['ENSG00000170439',]))
METTL7B_gene_data$GROUP = sample_info$GROUP
METTL7B_gene_data$SEX = sample_info$SEX
METTL7B_gene_data$AGE = sample_info$AGE
METTL7B_gene_data$NEUTROPHILS = sample_info$NEUTROPHILS
METTL7B_gene_data$MONOCYTES = sample_info$MONOCYTES
names(METTL7B_gene_data)[1] = 'GENE'

HP_gene_data = data.frame(t(expression_table['ENSG00000257017',]))
HP_gene_data$GROUP = sample_info$GROUP
HP_gene_data$SEX = sample_info$SEX
HP_gene_data$AGE = sample_info$AGE
HP_gene_data$NEUTROPHILS = sample_info$NEUTROPHILS
HP_gene_data$MONOCYTES = sample_info$MONOCYTES
names(HP_gene_data)[1]='GENE'

#Are these genes affected by any of the clinical measurements, such as Age, Sex, Neutrophils, Monocytes?
#since all of our data has a low sample size, and age and sex are skewed, and WBCs may or may not be skewed, we should use a non-parametric test like Spearman correlation.

cor(METTL7B_gene_data$AGE,METTL7B_gene_data$GENE,method='spearman')
#0.07138805 - weak positive correlation

#we don't run a test for sex because females are underrepresented.

cor(METTL7B_gene_data$NEUTROPHILS,METTL7B_gene_data$GENE,method='spearman')
#0.6892752 - moderate positive correlation

cor(METTL7B_gene_data$MONOCYTES,METTL7B_gene_data$GENE,method='spearman')
#-0.02917224 - weak negative correlation

#we can visualize these using scatterplots.
ggp19 = ggplot(METTL7B_gene_data,aes(x=AGE,y=GENE,colour=GROUP))+geom_point()+labs(title='Scatterplot of METTL7B Gene Expression VS Age')
#as age increases, gene expression increases, but in any case younger ages are not well-represented.
ggp20 = ggplot(METTL7B_gene_data,aes(x=NEUTROPHILS,y=GENE,colour=GROUP))+geom_point()+labs(title='Scatterplot of METTL7B Gene Expression VS Neutrophils')
#seemingly moderate positive correlation.
ggp21 = ggplot(METTL7B_gene_data,aes(x=MONOCYTES,y=GENE,colour=GROUP))+geom_point()+labs(title='Scatterplot of METTL7B Gene Expression VS Monocytes')
#weak negative correlation.

cor(HP_gene_data$AGE,HP_gene_data$GENE,method='spearman')
#0.127601 - weak positive correlation

#we don't run a test for sex because females are underrepresented.

cor(HP_gene_data$NEUTROPHILS,HP_gene_data$GENE,method='spearman')
#0.7172028 - moderate positive correlation

cor(HP_gene_data$MONOCYTES,HP_gene_data$GENE,method='spearman')
#-0.03784134 - weak negative correlation

#we can visualize these using scatterplots.
ggp22 = ggplot(HP_gene_data,aes(x=AGE,y=GENE,colour=GROUP))+geom_point()+labs(title='Scatterplot of HP Gene Expression VS Age')
#as age increases, gene expression increases, but in any case younger ages are not well-represented.
ggp23 = ggplot(HP_gene_data,aes(x=NEUTROPHILS,y=GENE,colour=GROUP))+geom_point()+labs(title='Scatterplot of HP Gene Expression VS Neutrophils')
#seemingly moderate positive correlation.
ggp24 = ggplot(HP_gene_data,aes(x=MONOCYTES,y=GENE,colour=GROUP))+geom_point()+labs(title='Scatterplot of HP Gene Expression VS Monocytes')
#weak negative correlation.

#Running gout vs SA analysis!
#first we make vectors of column names, due to a difference in case, we cannot directly take these from the sample_info file. I've chosen to just hard-code it.
expression_table_colnames_hc = c('Healthy_1','Healthy_2','Healthy_3','Healthy_4','Healthy_5','Healthy_6','Healthy_7','Healthy_8','Healthy_9','Healthy_10','Healthy_11','Healthy_12','Healthy_13','Healthy_14')
expression_table_colnames_gout = c('Gout_1','Gout_2','Gout_3','Gout_4','Gout_5','Gout_6','Gout_7','Gout_8','Gout_9','Gout_10','Gout_11','Gout_12','Gout_13','Gout_14')
expression_table_colnames_sa = row.names(sample_info_sa)

#making the expression tables for both gout and SA.
expression_table_gout = expression_table[,expression_table_colnames_gout]
expression_table_sa = expression_table[,expression_table_colnames_sa]

#creating an empty table to store our variables in.
de_gout_vs_sa = as.data.frame(matrix(0,ncol = 2,nrow = nrow(expression_table)))
names(de_gout_vs_sa) = c('Log2FC','p')
row.names(de_gout_vs_sa)=row.names(expression_table)

#making the loop and calculating fold change and p values.
for (row in 1:nrow(expression_table))
{
  gene_data_gout = as.numeric(expression_table_gout[row,])
  gene_data_sa = as.numeric(expression_table_sa[row,])
  log2fold = log2(mean(gene_data_gout))-log2(mean(gene_data_sa))
  p = (t.test(gene_data_sa,gene_data_gout))$p.value
  de_gout_vs_sa[row,1]=log2fold
  de_gout_vs_sa[row,2]=p
}

#adding p adjusted column using p.adjust from base R (stats package). specifically, the correction method of BH has been used.
de_gout_vs_sa$p.adjBH = p.adjust(as.numeric(de_gout_vs_sa$p),method = 'BH')
de_gout_vs_sa$p.adjHOLM = p.adjust(as.numeric(de_gout_vs_sa$p),method = 'holm')
de_gout_vs_sa$p.adjBONFERRONI = p.adjust(as.numeric(de_gout_vs_sa$p),method = 'bonferroni')
de_gout_vs_sa$p.adjBY = p.adjust(as.numeric(de_gout_vs_sa$p),method = 'BY')
summary(de_gout_vs_sa)
